#include "config.h"

#include <math.h>
#include <stdio.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "rt_solver.h"
#include "rt_utilities.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "halo_finder.h"
#include "ifrit.h"
#include "utils.h"


float** iUniformGrid_Sample(int level, int nbin[3], double bb[6], int nvar, int *varid);


#ifdef STARFORM
extern double sf_min_stellar_particle_mass;
#endif


int iOutputMesh(const char *filename, int level, int *nbinIn, double *bbIn, int nvars, int *varid)
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

  vars = iUniformGrid_Sample(level,nbin,bb,nvars,varid);
 
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


#ifdef PARTICLES

double gasdev()
{
  static int iset = 0;
  static double gset;
  double v1, v2, r2, fac;

  if(iset)
    {
      iset = 0;
      return gset;
    }
  else
    {
      do
	{
	  v1 = 2*cart_rand() - 1;
	  v2 = 2*cart_rand() - 1;
	  r2 = v1*v1 + v2*v2;
	}
      while(r2>1.0 || r2<1.0e-300);

      iset = 1;
      fac = sqrt(-2*log(r2)/r2);
      gset = v1*fac;
      return v2*fac;
    }
}


long iParticles_WriteArray(const char *filename, double *bb, int flags, float *arr, int idx, int nrec)
{
  const float dr = 0.33;
  int i, j, k, size, rank, ntemp, nsub;
  char s[999];
  double pos;
  float val;
  long ntot, nloc = 0L;
  double tot = 0.0;
  FILE *F;
#ifdef STARFORM
  float dm_star_min = sf_min_stellar_particle_mass * constants->Msun / units->mass;
#endif

  /*
  //  Write to a file in order of a proc rank
  */
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(local_proc_id == MASTER_NODE)
    {
      if(filename != NULL)
	{
	  F = fopen(filename,"a");
	  cart_assert(F != NULL);
	  if(idx > -1) ntemp = nrec*sizeof(double); else ntemp = nrec*sizeof(float);
	  fwrite(&ntemp,sizeof(int),1,F);
	  fclose(F);
	}
    }

  for(i=0; i<size; i++)
    {
      MPI_Barrier(MPI_COMM_WORLD);
      if(i == rank)
	{
	  if(filename != NULL)
	    {
	      F = fopen(filename,"a");
	      cart_assert(F != NULL);
	    }
	  for(j=0; j<num_particles; j++)
	    {
	      if(particle_id[j]!=NULL_PARTICLE &&
#ifdef STARFORM
		 particle_id_is_star(particle_id[j])==(flags&I_FLAG_STARS) &&
#endif
		 particle_x[j][0]>bb[0] && particle_x[j][0]<bb[1]
#if (nDim > 1)
		 && particle_x[j][1]>bb[2] && particle_x[j][1]<bb[3]
#if (nDim > 2)
		 && particle_x[j][2]>bb[4] && particle_x[j][2]<bb[5]
#endif
#endif
		 )
		{
		  nsub = 1;
#ifdef STARFORM
		  if(particle_id_is_star(particle_id[j])==(flags&I_FLAG_STARS) && (flags&I_FLAG_SPLIT_STARS))
		    {
		      nsub = (int)(0.5+particle_mass[j]/dm_star_min);
		      cart_assert(nsub > 0);
		    }
#endif
		  nloc += nsub;
		  if(filename != NULL)
		    {
		      if(idx > -1)
			{
			  fwrite(particle_x[j]+idx,sizeof(double),1,F);
			  for(k=1; k<nsub; k++)
			    {
			      pos = particle_x[j][idx] + dr*cell_size[particle_level[j]]*gasdev();
			      fwrite(&pos,sizeof(double),1,F);
			    }
			}
		      else
			{
			  val = arr[j];
			  tot += val;
			  if(flags & I_FLAG_ATTR_IS_MASS) val /= nsub;
			  for(k=0; k<nsub; k++)
			    {
			      fwrite(&val,sizeof(float),1,F);
			    }
			}
		    }
		}
	    }
	  if(filename != NULL)
	    {
	      fclose(F);
	    }
	}
    }

  if(flags & I_FLAG_ATTR_IS_MASS)
    {
      cart_debug("Total particle mass: %le",tot*units->mass/constants->Msun);
    }

  if(local_proc_id == MASTER_NODE)
    {
      if(filename != NULL)
	{
	  F = fopen(filename,"a");
	  cart_assert(F != NULL);
	  fwrite(&ntemp,sizeof(int),1,F);
	  fclose(F);
	}
    }

  MPI_Allreduce(&nloc,&ntot,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

  return ntot;

}


long iParticles_WriteBasicFile(const char *filename, double *bb, int flags)
{
  int i, j, k, size, rank, ntemp;
  FILE *F;
  float w;
  long ntot;

  ntot = iParticles_WriteArray(NULL,bb,flags,NULL,-1,0);

  cart_debug("Saving %ld particles...",ntot);

  if(ntot*4 != (int)(ntot*4))
    {
      cart_debug("Too many particles to output: %ld.",ntot);
      return 0L;
    }

  if(local_proc_id == MASTER_NODE)
    {
      F = fopen(filename,"w");
      if(F == 0)
	{
	  cart_debug("Unable to open file %s for writing.",filename);
	  return 0L;
	}
      ntemp = sizeof(int); fwrite(&ntemp,sizeof(int),1,F);
      ntemp = ntot;        fwrite(&ntemp,sizeof(int),1,F);
      ntemp = sizeof(int); fwrite(&ntemp,sizeof(int),1,F);

      ntemp = 6*sizeof(float); fwrite(&ntemp,sizeof(int),1,F);
      w = bb[0]; fwrite(&w,sizeof(float),1,F);
      w = bb[2]; fwrite(&w,sizeof(float),1,F);
      w = bb[4]; fwrite(&w,sizeof(float),1,F);
      w = bb[1]; fwrite(&w,sizeof(float),1,F);
      w = bb[3]; fwrite(&w,sizeof(float),1,F);
      w = bb[5]; fwrite(&w,sizeof(float),1,F);
      ntemp = 6*sizeof(float); fwrite(&ntemp,sizeof(int),1,F);
      fclose(F);
    }

  iParticles_WriteArray(filename,bb,flags,NULL,0,ntot);
#if (nDim > 1)
  iParticles_WriteArray(filename,bb,flags,NULL,1,ntot);
#if (nDim > 2)
  iParticles_WriteArray(filename,bb,flags,NULL,2,ntot);
#endif
#endif

  iParticles_WriteArray(filename,bb,flags|I_FLAG_ATTR_IS_MASS,particle_mass,-1,ntot);

  return ntot;
}


void iOutputParticles(const char *fileroot, double *bb)
{
  int j;
  float *arr;
  char str[999];
  long ntot;

  /*
  //  Dark matter particles
  */
  strcpy(str,fileroot);
  strcat(str,"-parts.bin");
  iParticles_WriteBasicFile(str,bb,0);

#ifdef STARFORM
  /*
  //  Stellar particles
  */
  strcpy(str,fileroot);
  strcat(str,"-stars.bin");
  ntot = iParticles_WriteBasicFile(str,bb,I_FLAG_STARS);

  iParticles_WriteArray(str,bb,I_FLAG_STARS,star_initial_mass,-1,ntot);

  arr = cart_alloc(float,num_particles);

  for(j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id_is_star(particle_id[j]))
    {
      arr[j] = tphys_from_tcode(particle_t[j]) - tphys_from_tcode(star_tbirth[j]);
    }

  iParticles_WriteArray(str,bb,I_FLAG_STARS,arr,-1,ntot);

  cart_free(arr);

#ifdef ENRICH
  iParticles_WriteArray(str,bb,I_FLAG_STARS,star_metallicity_II,-1,ntot);
#ifdef ENRICH_SNIa
  iParticles_WriteArray(str,bb,I_FLAG_STARS,star_metallicity_Ia,-1,ntot);
#endif /* ENRICH_SNIa */
#endif /* ENRICH */

#endif /* STARFORM */
}
#endif /* PARTICLES */


void iOutputHalo(const char *fileroot, int floor_level, float zoom, const halo *h, int nvars, int *varid)
{
  const int nbin1 = 256;
  int j, nbin[] = { nbin1, nbin1, nbin1 };
  double bb[6], pos[3], dbb;
  float dmax;
  char str[999];

  if(h == NULL)
    {
#ifdef HYDRO
      extFindMaxVar(HVAR_GAS_DENSITY,&dmax,pos,-1.0);
#else
      extFindMaxVar(VAR_DENSITY,&dmax,pos,-1.0);
#endif
    }
  else
    {
      for(j=0; j<nDim; j++) pos[j] = h->pos[j];
    }
  
  dbb = zoom*nbin1*pow(0.5,(double)floor_level);

  bb[0] = pos[0] - 0.5*dbb;
  bb[1] = pos[0] + 0.5*dbb;
  bb[2] = pos[1] - 0.5*dbb;
  bb[3] = pos[1] + 0.5*dbb;
  bb[4] = pos[2] - 0.5*dbb;
  bb[5] = pos[2] + 0.5*dbb;

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Rebinning for IFrIT at (%lf,%lf,%lf) with size %lf",pos[0],pos[1],pos[2],dbb);
    }

  strcpy(str,fileroot);
  strcat(str,"-mesh.bin");
  iOutputMesh(str,max_level,nbin,bb,nvars,varid);

#ifdef PARTICLES
  iOutputParticles(fileroot,bb);
#endif /* PARTICLES */
}


/*
// ************************************************
//
//              PRIVATE FUNCTIONS
// (should not be called from outside of this file)
//
// ************************************************
*/
long iUniformGrid_GetSize(int level, int nbin[3], double bb[6])
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


void iUniformGrid_FillData(int level, int nbin[3], double bb[6], int nvars, int *varid, float **buf, long *loc)
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
#ifdef HYDRO
		      else if(varid[var]>=I_FRACTION && varid[var]<I_FRACTION+num_vars)
			{
			  buf[var][l] = cell_var(cell,varid[var]-I_FRACTION)/cell_gas_density(cell);
			}
#endif /* HYDRO */
		      else switch(varid[var])
			{
#ifdef RADIATIVE_TRANSFER
			case I_GAS_TEMPERATURE:
			  {
			    buf[var][l] = units->temperature*cell_gas_temperature(cell);
			    break;
			  }
#endif /* RADIATIVE_TRANSFER */
			case I_CELL_LEVEL:
			  {
			    buf[var][l] = cell_level(cell);
			    break;
			  }
			case I_LOCAL_PROC:
			  {
			    buf[var][l] = local_proc_id;
			    break;
			  }
#ifdef HYDRO
			case I_GAS_NUMBER_DENSITY:
			  {
			    buf[var][l] = units->number_density*cell_gas_density(cell);
			    break;
			  }
#endif /* HYDRO */
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


float** iUniformGrid_Sample(int level, int nbin[3], double bb[6], int nvars, int *varid)
{
  int i, ip;
  long l, ncells;
  float **vars, **buf;
  long *loc;

  buf = cart_alloc(float*, nvars );

  if(local_proc_id == MASTER_NODE)
    {
      vars = cart_alloc(float*, nvars );
      for(i=0; i<nvars; i++)
	{
	  vars[i] = cart_alloc(float, nbin[0]*nbin[1]*nbin[2] );
	}

      for(ip=0; ip<num_procs; ip++)
	{
	  /*
	  //  Measure buffer size
	  */
	  if(ip == 0)
	    {
	      ncells = iUniformGrid_GetSize(level,nbin,bb);
	    }
	  else
	    {
	      MPI_Recv(&ncells,1,MPI_LONG,ip,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }

	  /*
	  //  Allocate buffers
	  */
	  for(i=0; i<nvars; i++) buf[i] = cart_alloc(float, ncells );
	  loc = cart_alloc(long, ncells );

	  /*
	  //  Fill/transfer buffers
	  */
	  if(ip == 0)
	    {
	      iUniformGrid_FillData(level,nbin,bb,nvars,varid,buf,loc);
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
      ncells = iUniformGrid_GetSize(level,nbin,bb);
      MPI_Send(&ncells,1,MPI_LONG,MASTER_NODE,0,MPI_COMM_WORLD);

      /*
      //  Allocate buffers
      */
      for(i=0; i<nvars; i++) buf[i] = cart_alloc(float, ncells );
      loc = cart_alloc(long, ncells );

      /*
      //  Fill & transfer buffers
      */
      iUniformGrid_FillData(level,nbin,bb,nvars,varid,buf,loc);

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

