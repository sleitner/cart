#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(STAR_FORMATION) && defined(COSMOLOGY)

#include <math.h>
#include <stdio.h>

#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rt.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "fesc.h"

#include "cd.h"
#include "halo_tools.h"


extern double sf_min_stellar_particle_mass;
extern int current_step_level;


int fe_912A_Vars_[] = { RT_HVAR_HI_DENSITY, RT_HVAR_H2_DENSITY };
double fe_912A_Coll_(int num, double *vals)
{
  return exp(-6.3e-18*(vals[0]+2*vals[1])*units->number_density*units->length);
}


struct feInfo fe_912A = { 2, fe_912A_Vars_, fe_912A_Coll_ };


void feAddOneLocation(struct feInfo a, int nside, float w, int cell, double pos[3], double len, double s[], int l, int n)
{
  int npix = 12*nside*nside;
  int i, k;
  cdData outputs[npix];
  float tmp[a.num];
  double tnext;
  static double tprev = -1;

  s[npix] += w;

  for(i=0; i<a.num; i++) tmp[i] = cell_var(cell,a.var_ids[i]);
  for(i=0; i<a.num; i++) cell_var(cell,a.var_ids[i]) = 0.0;

  cdTraverseSky(a.num,a.var_ids,nside,pos,len,max_level,outputs);

  for(i=0; i<a.num; i++) cell_var(cell,a.var_ids[i]) = tmp[i];

  for(k=0; k<npix; k++)
    {
      s[k] += w*a.collector(a.num,outputs[k].val);
    }

  tnext = MPI_Wtime();
  if(tprev < 0.0) tprev = tnext;
  if(tnext > tprev+60)
    {
      tprev = tnext;
      cart_debug("%d out of %d (%4.1f%)",l,n,100.0*l/n);
    }

  if(l == n) tprev = -1;
}


void feSumOverParts(struct feInfo a, int nside, float scale, const halo *h, double sum[], int sample)
{
  int l, n;
  float w;
  int j, cell;

  cart_assert(nside > 0);
  cart_assert(h);

  rtSetupSource(min_level);

  for(n=j=0; j<num_particles; j++) if(particle_is_star(j) && compute_distance_periodic(particle_x[j],h->pos)<h->rhalo)
    {
      w = particle_mass[j]*rtSource(j);
      if(w > 0) n++;
    }

  if(n == 0) return;

  cart_debug("Ray-tracing %d particles...",n);
  current_step_level += 3;

  for(l=j=0; j<num_particles; j++) if(particle_is_star(j) && compute_distance_periodic(particle_x[j],h->pos)<h->rhalo)
    {
      w = particle_mass[j]*rtSource(j);
      if(w > 0)
	{
	  l++;
	  if(l%sample == 0)
	    {
	      cell = cell_find_position(particle_x[j]);
	      feAddOneLocation(a,nside,w,cell,particle_x[j],scale*h->rhalo,sum,l,n);
	    }
	}
    }

  current_step_level -= 3;
}


void feSumOverCells(struct feInfo a, int nside, float scale, const halo *h, double sum[], int sample)
{
  int l, n;
  float w;
  double pos[3];
  MESH_RUN_DECLARE(level,cell);

  cart_assert(nside > 0);
  cart_assert(h);

  n = 0;
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  cell_center_position(cell,pos);
  if(cell_is_leaf(cell) && compute_distance_periodic(pos,h->pos)<h->rhalo)
    {
      w = cell_rt_source(cell)*cell_volume[level];
      if(w > 0) n++;
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;

  if(n == 0) return;

  cart_debug("Ray-tracing %d cells...",n);
  current_step_level += 3;

  l = 0;
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  cell_center_position(cell,pos);
  if(cell_is_leaf(cell) && compute_distance_periodic(pos,h->pos)<h->rhalo)
    {
      w = cell_rt_source(cell)*cell_volume[level];
      if(w > 0)
	{
	  l++;
	  if(l%sample == 0) feAddOneLocation(a,nside,w,cell,pos,scale*h->rhalo,sum,l,n);
	}
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;

  current_step_level -= 3;
}


void extEscapeFraction(const char *fname, struct feInfo a, int nside, struct HALO_LIST *halos, int resolution_level, float outer_edge, int sample, int sum_over_parts)
{
  int npix = 12*nside*nside;
  int i, j, ih;
  double sum[npix+1], fesc[npix+1], sfr;
  FILE *f;
  double t;

  cart_assert(sample >= 1);

  f = fopen(fname,"w");
  cart_assert(f != NULL);

  fprintf(f,"# Columns:\n"); 
  fprintf(f,"# 1. Halo id\n"); 
  fprintf(f,"# 2. Halo mass (Msun)\n"); 
  fprintf(f,"# 3. SFR (Msun/yr)\n"); 
  fprintf(f,"# 4-6. <fesc>, fesc_{10%}, fesc_{90%}\n"); 
  fprintf(f,"#\n"); 

  for(ih=0; ih<halos->num_halos; ih++) if(halo_level(&(halos->list[ih]),mpi.comm.run)>=resolution_level)
    {
      t = MPI_Wtime();

      sfr = 0.0;
      for(j=0; j<num_particles; j++) if(particle_is_star(j) && units->time*(tl[min_level]-star_tbirth[j])<10*constants->Myr && compute_distance_periodic(particle_x[j],halos->list[ih].pos)<halos->list[ih].rhalo)
	{
	  sfr += particle_mass[j];
	}
      sfr = units->mass*sfr/constants->Msun/1.0e7;

      if(sfr < 1.0e-10) continue;

      for(i=0; i<=npix; i++) sum[i] = 0.0;

      if(sum_over_parts)
	{
	  feSumOverParts(a,nside,outer_edge,&(halos->list[ih]),sum,sample);
	}
      else
	{
	  feSumOverCells(a,nside,outer_edge,&(halos->list[ih]),sum,sample);
	}

      MPI_Allreduce(sum,fesc,npix+1,MPI_DOUBLE,MPI_SUM,mpi.comm.run);

      if(fesc[npix] > 0.0)
	{
	  for(i=0; i<npix; i++) fesc[i] /= fesc[npix];

	  for(i=0; i<npix-1; i++) for(j=i+1; j<npix; j++) if(fesc[i] > fesc[j])
	    {
	      fesc[npix] = fesc[i];
	      fesc[i] = fesc[j];
	      fesc[j] = fesc[npix];
	    }
	  
	  fesc[npix] = 0.0;
	  for(i=0; i<npix; i++) fesc[npix] += fesc[i];
	  fesc[npix] /= npix;

	  t = MPI_Wtime() - t;

	  cart_debug("Halo #%d, Fesc = %le %le %le (%7.1le s)",halos->list[ih].id,fesc[npix],fesc[nside*nside],fesc[11*nside*nside],t);

	  fprintf(f,"%d %le %le %le %le %le\n",halos->list[ih].id,units->mass*halos->list[ih].mvir/constants->Msun,sfr,fesc[npix],fesc[nside*nside],fesc[11*nside*nside]);
	  fflush(f);
	}
    }

  fclose(f);
}

#endif /* RADIATIVE_TRANSFER && STAR_FORMATION && COSMOLOGY */

