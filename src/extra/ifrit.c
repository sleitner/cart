#include "defs.h"      


#include <mpi.h>
#include <stdio.h>

#include "../auxiliary.h"
#include "../tree.h"

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#endif

#include "ifrit.h"


float** extUniformGrid_Sample(int level, int nbin[3], double bb[6], int nvar, int *varid);


int extWriteIfritFile(int level, int *nbinIn, double *bbIn, int nvars, int *varid, const char *filename)
{
  int i, ntemp, nbin[3];
  double bb[6];
  float **vars;
  FILE *F;

  cart_assert(nvars > 0);
  for(i=0; i<nDim; i++)
    {
      cart_assert(nbinIn[i] > 0);
      cart_assert(bbIn[2*i+1] > bbIn[2*i]);

      nbin[i] = nbinIn[i];
      bb[2*i+0] = bbIn[2*i+0];
      bb[2*i+1] = bbIn[2*i+1];
    }

  for(i=nDim; i<3; i++)
    {
      nbin[i] = 1;
      bb[2*i+0] = 0.0;
      bb[2*i+1] = num_grid;
    }

  vars = extUniformGrid_Sample(level,nbin,bb,nvars,varid);
 
  if(local_proc_id == MASTER_NODE)
    {
      F = fopen(filename,"w");
      if(F == 0)
	{
	  cart_debug("Unable to open file %s for writing.",filename);
	  for(i=0; i<nvars; i++) cart_free(vars[i]);
	  cart_free(vars);
	  return 1;
	}
  
      ntemp = 12;
      fwrite(&ntemp,4,1,F);
      fwrite(nbin+0,4,1,F);
      fwrite(nbin+1,4,1,F);
      fwrite(nbin+2,4,1,F);
      fwrite(&ntemp,4,1,F);
      ntemp = 4*nbin[0]*nbin[1]*nbin[2];
      for(i=0; i<nvars; i++)
	{
	  fwrite(&ntemp,4,1,F); fwrite(vars[i],4,nbin[0]*nbin[1]*nbin[2],F); fwrite(&ntemp,4,1,F);
	}
      fclose(F);

      for(i=0; i<nvars; i++) cart_free(vars[i]);
      cart_free(vars);
    }

   return 0;
}


/*
// ************************************************
//
//              PRIVATE FUNCTIONS
// (should not be called from outside of this file)
//
// ************************************************
*/
long extUniformGrid_GetSize(int level, int nbin[3], double bb[6])
{
  int i, j, k, cell;
  long n;
  double pos[3];

  /* count number of buffer cells */
  n = 0;
  for(k=0; k<nbin[2]; k++)
    {
      pos[2] = bb[4] + (bb[5]-bb[4])*(k+0.5)/nbin[2];
      for(j=0; j<nbin[1]; j++)
	{
	  pos[1] = bb[2] + (bb[3]-bb[2])*(j+0.5)/nbin[1];
	  for(i=0; i<nbin[0]; i++)
	    {
	      pos[0] = bb[0] + (bb[1]-bb[0])*(i+0.5)/nbin[0];
	      cell = cell_find_position_above_level(level,pos);
	      if(cell!=-1 && root_cell_type(root_cell_sfc_index(cell_parent_root_cell(cell)))==1)
		{
		  n++;
		}
	    }
	}
    }
  
  return n;
}


void extUniformGrid_FillData(int level, int nbin[3], double bb[6], int nvars, int *varid, float **buf, long *loc)
{
  int i, j, k, var, cell;
  long offset, l;
  double pos[3];

  l = 0;
  for(k=0; k<nbin[2]; k++)
    {
      pos[2] = bb[4] + (bb[5]-bb[4])*(k+0.5)/nbin[2];
      for(j=0; j<nbin[1]; j++)
	{
	  pos[1] = bb[2] + (bb[3]-bb[2])*(j+0.5)/nbin[1];
	  offset = nbin[1]*(j+nbin[2]*k);
	  for(i=0; i<nbin[0]; i++)
	    {
	      pos[0] = bb[0] + (bb[1]-bb[0])*(i+0.5)/nbin[0];
	      cell = cell_find_position_above_level(level,pos);
	      if(cell!=-1 && root_cell_type(root_cell_sfc_index(cell_parent_root_cell(cell)))==1)
		{
		  for(var=0; var<nvars; var++)
		    {
		      if(varid[var]>=0 && varid[var]<num_vars)
			{
			  buf[var][l] = cell_var(cell,varid[var]);
			}
		      else if(varid[var]>=EXT_FRACTION && varid[var]<EXT_FRACTION+num_vars)
			{
			  buf[var][l] = cell_var(cell,varid[var]-EXT_FRACTION)/cell_gas_density(cell);
			}
		      else switch(varid[var])
			{
#ifdef RADIATIVE_TRANSFER
			case EXT_GAS_TEMPERATURE:
			  {
			    buf[var][l] = rtTemInK(cell);
			    break;
			  }
#endif  // RADIATIVE_TRANSFER
			case EXT_CELL_LEVEL:
			  {
			    buf[var][l] = cell_level(cell);
			    break;
			  }
			case EXT_LOCAL_PROC:
			  {
			    buf[var][l] = local_proc_id;
			    break;
			  }
			default:
			  {
			    buf[var][l] = 0.0;
			  }
			}
		    }
		  loc[l] = i + offset;
		  l++;
		}
	    }
	}
    }
}


float** extUniformGrid_Sample(int level, int nbin[3], double bb[6], int nvars, int *varid)
{
  int i, ip, done;
  long l, ncells;
  float **vars, **buf;
  long *loc;

  buf = cart_alloc(nvars*sizeof(float*));

  if(local_proc_id == MASTER_NODE)
    {
      vars = cart_alloc(nvars*sizeof(float*));
      for(i=0; i<nvars; i++)
	{
	  vars[i] = cart_alloc(nbin[0]*nbin[1]*nbin[2]*sizeof(float));
	}

      for(ip=0; ip<num_procs; ip++)
	{
	  /*
	  //  Measure buffer size
	  */
	  if(ip == 0)
	    {
	      ncells = extUniformGrid_GetSize(level,nbin,bb);
	    }
	  else
	    {
	      MPI_Recv(&ncells,1,MPI_LONG,ip,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }

	  /*
	  //  Allocate buffers
	  */
	  for(i=0; i<nvars; i++) buf[i] = cart_alloc(ncells*sizeof(float));
	  loc = cart_alloc(ncells*sizeof(long));

	  /*
	  //  Fill/transfer buffers
	  */
	  if(ip == 0)
	    {
	      extUniformGrid_FillData(level,nbin,bb,nvars,varid,buf,loc);
	    }
	  else
	    {
	      MPI_Recv(loc,ncells,MPI_LONG,ip,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      for(i=0; i<nvars; i++) MPI_Recv(buf[i],ncells,MPI_FLOAT,ip,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }

	  /*
	  //  Fill variable arrays
	  */
	  for(i=0; i<nvars; i++)
	    {
	      for(l=0; l<ncells; l++)
		{
		  vars[i][loc[l]] = buf[i][l];
		}
	    }
	  
	  /*
	  //  Free buffers
	  */
	  for(i=0; i<nvars; i++) cart_free(buf[i]);
	  cart_free(loc);
 	}

    }
  else
    {
      vars = 0;

      /*
      //  Measure & transfer buffer size
      */
      ncells = extUniformGrid_GetSize(level,nbin,bb);
      MPI_Send(&ncells,1,MPI_LONG,MASTER_NODE,0,MPI_COMM_WORLD);

      /*
      //  Allocate buffers
      */
      for(i=0; i<nvars; i++) buf[i] = cart_alloc(ncells*sizeof(float));
      loc = cart_alloc(ncells*sizeof(long));

      /*
      //  Fill & transfer buffers
      */
      extUniformGrid_FillData(level,nbin,bb,nvars,varid,buf,loc);

      MPI_Send(loc,ncells,MPI_LONG,MASTER_NODE,0,MPI_COMM_WORLD);
      for(i=0; i<nvars; i++) MPI_Send(buf[i],ncells,MPI_FLOAT,MASTER_NODE,0,MPI_COMM_WORLD);

      /*
      //  Free buffers
      */
      for(i=0; i<nvars; i++) cart_free(buf[i]);
      cart_free(loc);
    }

  cart_free(buf);

  return vars;
}

