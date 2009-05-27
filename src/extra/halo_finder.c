#include "defs.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "parallel.h"
#include "timestep.h"
#include "units.h"

#include "halo_finder.h"


int hfReadHFINDHalos(const char* fname, int nmemMin, float massMin, float vmaxMin, float rvirMin, hfHalo **list, int *list_size)
{
  int i, j, ret, size, nread = 0;
  long np;
  FILE *f;
  char buffer[1024];
  float uPos, uVel;
  float a, OmM, OmL, OmB, h; 
  hfHalo q, *d, *tmp;

  ret = -1;

  while(ret==-1 && local_proc_id==MASTER_NODE)
    {

      f = fopen(fname,"r");
      if(f == NULL) break;

      /*
      //  Skip the job name
      */
      if(fgets(buffer,1024,f) == NULL) break;
      
      /*
      //  Check the scale factor
      */
      if(fscanf(f," A=%f A0=%*f Ampl=%*f Step=%*f\n",&a) != 1) break;
      if(fabs(a-auni[min_level]) > 1.0e-3)
	{
	  cart_debug("Scalar factor in HFIND file (%f) is different from the current value (%f)",a,auni[min_level]);
	}

      if(fgets(buffer,1024,f) == NULL) break;

      if(fscanf(f," Nrow=%*d Ngrid=%*d  Omega_0=%f OmLam_0=%f  Omegab_0=%f Hubble=%f\n",&OmM,&OmL,&OmB,&h) != 4) break;
      if(fabs(OmM/cosmology->OmegaM-1.0) > 1.0e-3)
	{
	  cart_debug("OmegaM in HFIND file (%f) is different from the current value (%f)",OmM,cosmology->OmegaM);
	}
      if(fabs(OmL/cosmology->OmegaL-1.0) > 1.0e-3)
	{
	  cart_debug("OmegaL in HFIND file (%f) is different from the current value (%f)",OmL,cosmology->OmegaL);
	}
      if(fabs(OmB/cosmology->OmegaB-1.0) > 1.0e-3)
	{
	  cart_debug("OmegaB in HFIND file (%f) is different from the current value (%f)",OmB,cosmology->OmegaB);
	}
      if(fabs(h/cosmology->h-1.0) > 1.0e-3)
	{
	  cart_debug("Hubble in HFIND file (%f) is different from the current value (%f)",h,cosmology->h);
	}

      /*
      //  Skip the rest
      */
      for(i=4; i<17; i++)
	{
	  if(fgets(buffer,1024,f) == NULL) break;
	}

      uPos = r0;
      uVel = v0/abox[min_level];

      size = 100;
      d = cart_alloc(hfHalo,size);
  
      while(fscanf(f,"%d %lf %lf %lf %f %f %f %f %*f %e %ld %f %*f %*f %*ld\n",&q.Id,q.Pos+0,q.Pos+1,q.Pos+2,q.Vel+0,q.Vel+1,q.Vel+2,&q.Rvir,&q.Mass,&np,&q.Vmax) == 11)
	{
	  if(np>nmemMin && q.Mass>massMin && q.Vmax>vmaxMin && q.Rvir>rvirMin)
	    {
	      if(nread == size)
		{
		  size *= 2;
		  tmp = cart_alloc(hfHalo,size);
		  memcpy(tmp,d,nread*sizeof(hfHalo));
		  cart_free(d);
		  d = tmp;
		}

	      for(j=0; j<3; j++)
		{
		  q.Pos[j] /= uPos;
		  q.Vel[j] /= uVel;
		}

	      d[nread++] = q;
	    }
	}

      fclose(f);
      ret = 0;
    }

  MPI_Bcast(&ret,1,MPI_INT,MASTER_NODE,MPI_COMM_WORLD);
  MPI_Bcast(&nread,1,MPI_INT,MASTER_NODE,MPI_COMM_WORLD);
  if(nread > 0)
    {
      if(local_proc_id != MASTER_NODE) d = cart_alloc(hfHalo,nread);
      MPI_Bcast(d,nread*sizeof(hfHalo),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);
    }
  else
    {
      if(local_proc_id == MASTER_NODE) cart_free(d);
      d = 0;
    }

  *list = d;
  *list_size = nread;

  return ret;
} 


int hfHaloLevel(const hfHalo *h)
{
  int llev, glev, cell;

  cart_assert(h);

  cell = cell_find_position(h->Pos);

  if(cell < 0)
    {
      llev = min_level - 1;
    }
  else
    {
      llev = cell_level(cell);
    }

  MPI_Allreduce(&llev,&glev,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  return glev;
}
