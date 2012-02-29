#include "config.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "hydro.h"
#include "parallel.h"
#include "particle.h"
#include "rand.h"
#include "rt.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "halo_finder.h"
#include "ifrit.h"


float** ifritUniformGrid_Sample(int level, int nbin[3], double pcen[3], int nvar, const int *varid);


#ifdef STARFORM
extern double sf_min_stellar_particle_mass;
#endif


int ifritOutputMesh(const char *filename, int level, int nbinIn[], const double pcenIn[], int nvars, const int *varid)
{
  int i, ntemp, nbin[3];
  double pcen[3];
  float **vars;
  FILE *F;

  cart_assert(nvars > 0);
  for(i=0; i<nDim; i++)
    {
      cart_assert(nbinIn[i] > 0);

      nbin[i] = nbinIn[i];
      pcen[i] = pcenIn[i];
    }

  for(i=nDim; i<3; i++)
    {
      nbin[i] = 1;
      pcen[i] = 0.5;
    }


  if(level < 0)
    {
      level = min_level;
      ntemp = num_grid;
      while(ntemp < nbin[0])
	{
	  ntemp *= 2;
	  level++;
	}

      cart_debug("Level selected for IFrIT::OutputMesh: %d",level);
    }

  vars = ifritUniformGrid_Sample(level,nbin,pcen,nvars,varid);
 
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

double ifrit_gasdev()
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


long ifritParticles_WriteArray(const char *filename, double *bb, int flags, float *arr, int idx, int nrec)
{
  const float dr = 0.33;
  int i, j, k, size, rank, ntemp, nsub;
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
  MPI_Comm_size(mpi.comm.run,&size);
  MPI_Comm_rank(mpi.comm.run,&rank);

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
      MPI_Barrier(mpi.comm.run);
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
			      pos = particle_x[j][idx] + dr*cell_size[particle_level[j]]*ifrit_gasdev();
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

  MPI_Allreduce(&nloc,&ntot,1,MPI_LONG,MPI_SUM,mpi.comm.run);

  return ntot;

}


long ifritParticles_WriteBasicFile(const char *filename, double *bb, int flags)
{
  int ntemp;
  FILE *F;
  float w;
  long ntot;

  ntot = ifritParticles_WriteArray(NULL,bb,flags,NULL,-1,0);

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

  ifritParticles_WriteArray(filename,bb,flags,NULL,0,ntot);
#if (nDim > 1)
  ifritParticles_WriteArray(filename,bb,flags,NULL,1,ntot);
#if (nDim > 2)
  ifritParticles_WriteArray(filename,bb,flags,NULL,2,ntot);
#endif
#endif

  ifritParticles_WriteArray(filename,bb,flags|I_FLAG_ATTR_IS_MASS,particle_mass,-1,ntot);

  return ntot;
}


void ifritOutputParticles(const char *fileroot, double *bb)
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
  ifritParticles_WriteBasicFile(str,bb,0);

#ifdef STARFORM
  /*
  //  Stellar particles
  */
  strcpy(str,fileroot);
  strcat(str,"-stars.bin");
  ntot = ifritParticles_WriteBasicFile(str,bb,I_FLAG_STARS);

  ifritParticles_WriteArray(str,bb,I_FLAG_STARS,star_initial_mass,-1,ntot);

  arr = cart_alloc(float,num_particles);

  for(j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id_is_star(particle_id[j]))
    {
#ifdef COSMOLOGY
      arr[j] = tphys_from_tcode(particle_t[j]) - tphys_from_tcode(star_tbirth[j]);
#else
      arr[j] = particle_t[j] - star_tbirth[j];
#endif
    }

  ifritParticles_WriteArray(str,bb,I_FLAG_STARS,arr,-1,ntot);

  cart_free(arr);

#ifdef ENRICH
  ifritParticles_WriteArray(str,bb,I_FLAG_STARS,star_metallicity_II,-1,ntot);
#ifdef ENRICH_SNIa
  ifritParticles_WriteArray(str,bb,I_FLAG_STARS,star_metallicity_Ia,-1,ntot);
#endif /* ENRICH_SNIa */
#endif /* ENRICH */

#endif /* STARFORM */
}
#endif /* PARTICLES */


void ifritOutputBox(const char *fileroot, int pixel_level, int nbin[], const double pos[3], int nvars, const int *varid)
{
  double bb[6], dx;
  char str[999];

  dx = pow(0.5,(double)pixel_level);

  bb[0] = pos[0] - 0.5*dx*nbin[0];
  bb[1] = pos[0] + 0.5*dx*nbin[0];
  bb[2] = pos[1] - 0.5*dx*nbin[1];
  bb[3] = pos[1] + 0.5*dx*nbin[1];
  bb[4] = pos[2] - 0.5*dx*nbin[2];
  bb[5] = pos[2] + 0.5*dx*nbin[2];

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Rebinning for IFrIT at (%lf,%lf,%lf) with cell size %lf",pos[0],pos[1],pos[2],dx);
    }

  strcpy(str,fileroot);
  strcat(str,"-mesh.bin");
  ifritOutputMesh(str,pixel_level,nbin,pos,nvars,varid);

#ifdef PARTICLES
  ifritOutputParticles(fileroot,bb);
#endif /* PARTICLES */
}


void ifritOutputHalo(const char *fileroot, int pixel_level, int nbin[], const halo *h, int nvars, const int *varid)
{
  if(h == NULL)
    {
      cart_debug("No halo file is loaded. Skipping ifritOutputHalo.");
      return;
    }
  else
    {
      ifritOutputBox(fileroot,pixel_level,nbin,h->pos,nvars,varid);
    }
}


/*
// ************************************************
//
//              PRIVATE FUNCTIONS
// (should not be called from outside of this file)
//
// ************************************************
*/
long ifritUniformGrid_GetSize(int level, int nbin[3], double pcen[3])
{
  int i, j, k, cell;
  long n;
  double pos[3], dx;

  dx = pow(0.5,(double)level);

  /* count number of buffer cells */
  n = 0;
  for(k=0; k<nbin[2]; k++)
    {
      pos[2] = pcen[2] + dx*(k+0.5-0.5*nbin[2]);
      pos[2] -= num_grid*floor(pos[2]/num_grid);
      for(j=0; j<nbin[1]; j++)
	{
	  pos[1] = pcen[1] + dx*(j+0.5-0.5*nbin[1]);
	  pos[1] -= num_grid*floor(pos[1]/num_grid);
	  for(i=0; i<nbin[0]; i++)
	    {
	      pos[0] = pcen[0] + dx*(i+0.5-0.5*nbin[0]);
	      pos[0] -= num_grid*floor(pos[0]/num_grid);
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


void ifritUniformGrid_FillData(int level, int nbin[3], double pcen[3], int nvars, const int *varid, float **buf, long *loc)
{
  int i, j, k, var, cell;
  long offset, l;
  double pos[3], dx;

  dx = pow(0.5,(double)level);

  l = 0;
  for(k=0; k<nbin[2]; k++)
    {
      pos[2] = pcen[2] + dx*(k+0.5-0.5*nbin[2]);
      pos[2] -= num_grid*floor(pos[2]/num_grid);
      for(j=0; j<nbin[1]; j++)
	{
	  pos[1] = pcen[1] + dx*(j+0.5-0.5*nbin[1]);
	  pos[1] -= num_grid*floor(pos[1]/num_grid);
	  offset = nbin[1]*(j+nbin[2]*k);
	  for(i=0; i<nbin[0]; i++)
	    {
	      pos[0] = pcen[0] + dx*(i+0.5-0.5*nbin[0]);
	      pos[0] -= num_grid*floor(pos[0]/num_grid);
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
			case I_HI_FRACTION:
			  {
			    buf[var][l] = cell_HI_fraction(cell);
			    break;
			  }
			case I_H2_FRACTION:
			  {
			    buf[var][l] = cell_H2_fraction(cell);
			    break;
			  }
			case I_DMW:
			  {
			    buf[var][l] = rtDmw(cell);
			    break;
			  }
			case I_UMW:
			  {
			    buf[var][l] = rtUmw(cell);
			    break;
			  }
#endif
#ifdef HYDRO
			case I_GAS_NUMBER_DENSITY:
			  {
			    buf[var][l] = units->number_density*cell_gas_density(cell);
			    break;
			  }
#ifdef COSMOLOGY
			case I_GAS_OVERDENSITY:
			  {
			    buf[var][l] = cell_gas_density(cell)*cosmology->OmegaM/cosmology->OmegaB;
			    break;
			  }
#endif /* COSMOLOGY */
			case I_GAS_TEMPERATURE:
			  {
			    buf[var][l] = units->temperature*cell_gas_temperature(cell);
			    break;
			  }
			case I_GAS_TOVERMU:
			  {
			    buf[var][l] = units->temperature*(cell_gas_gamma(cell)-1)*cell_gas_internal_energy(cell)/cell_gas_density(cell);
			    break;
			  }
#endif /* HYDRO */
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
#ifdef ENRICH
			case I_GAS_METAL_DENSITY:
			  {
			    buf[var][l] = cell_gas_metal_density(cell);
			    break;
			  }
			case I_GAS_METALLICITY:
			  {
			    buf[var][l] = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
			    break;
			  }
#endif /* ENRICH */



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


float** ifritUniformGrid_Sample(int level, int nbin[3], double pcen[3], int nvars, const int *varid)
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
	      ncells = ifritUniformGrid_GetSize(level,nbin,pcen);
	    }
	  else
	    {
	      MPI_Recv(&ncells,1,MPI_LONG,ip,0,mpi.comm.run,MPI_STATUS_IGNORE);
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
	      ifritUniformGrid_FillData(level,nbin,pcen,nvars,varid,buf,loc);
	    }
	  else
	    {
	      MPI_Recv(loc,ncells,MPI_LONG,ip,0,mpi.comm.run,MPI_STATUS_IGNORE);
	      for(i=0; i<nvars; i++) MPI_Recv(buf[i],ncells,MPI_FLOAT,ip,0,mpi.comm.run,MPI_STATUS_IGNORE);
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
      ncells = ifritUniformGrid_GetSize(level,nbin,pcen);
      MPI_Send(&ncells,1,MPI_LONG,MASTER_NODE,0,mpi.comm.run);

      /*
      //  Allocate buffers
      */
      for(i=0; i<nvars; i++) buf[i] = cart_alloc(float, ncells );
      loc = cart_alloc(long, ncells );

      /*
      //  Fill & transfer buffers
      */
      ifritUniformGrid_FillData(level,nbin,pcen,nvars,varid,buf,loc);

      MPI_Send(loc,ncells,MPI_LONG,MASTER_NODE,0,mpi.comm.run);
      for(i=0; i<nvars; i++) MPI_Send(buf[i],ncells,MPI_FLOAT,MASTER_NODE,0,mpi.comm.run);

      /*
      //  Free buffers
      */
      for(i=0; i<nvars; i++) cart_free(buf[i]);
      cart_free(loc);
    }

  cart_free(buf);

  return vars;
}

const struct IFRIT_NAMESPACE ifrit = { ifritOutputMesh, ifritOutputHalo, ifritOutputBox };

