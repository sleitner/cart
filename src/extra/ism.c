#include "config.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rt_solver.h"
#include "rt_utilities.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "cd.h"
#include "halo_finder.h"
#include "ism.h"


/*
//  Set cell_var(c,var) with the halo id for each halo, or 0 if belongs to 
//  none; a cell belongs to a halo if it is inside its size_factor*Rtrunc, 
//  and satellites are accounted for properly.
*/
void extMapHaloCells(int var, int floor_level, const halo_list *halos, float size_factor)
{
  int j, ih, iold, *halo_levels;
  MESH_RUN_DECLARE(level,cell);
  double pos[3], dx, r2, r2old, r2Cut;

  cart_assert(halos != NULL);
  cart_assert(var>=0 && var<num_vars);

  halo_levels = cart_alloc(int,halos->num_halos);
  for(ih=0; ih<halos->num_halos; ih++) halo_levels[ih] = halo_level(&halos->list[ih],MPI_COMM_WORLD);

  /*
  //  Loop over levels first to avoid selecting cells multiple times
  */
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  cart_debug("Mapping level %d...",level);

  /*
  //  Zero map array
  */
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,var)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell)) cell_var(cell,var) = 0.0;
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

  for(ih=0; ih<halos->num_halos; ih++) if(halo_levels[ih] >= floor_level)
    {
      r2Cut = pow(size_factor*halos->list[ih].rhalo,2.0);

      /*
      //  Map halo indices(+1), not ids, first (to simplify inter-comparison)
      */
#pragma omp parallel for default(none), private(_Index,cell,j,dx,r2,iold,r2old,pos), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,var,r2Cut,halos,ih)
      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
      if(cell_is_leaf(cell))
	{
	  cell_position_double(cell,pos);
	  for(j=0, r2=0.0; j<nDim; j++)
	    {
	      dx = pos[j] - halos->list[ih].pos[j];
	      if(dx < -0.5*num_grid) dx += num_grid;
	      if(dx >  0.5*num_grid) dx -= num_grid;
	      r2 += dx*dx;
	    }
	  
	  if(r2 < r2Cut)
	    {
	      iold = (int)(0.5+cell_var(cell,var));
	      if(iold == 0)
		{
		  cell_var(cell,var) = ih + 1;
		}
	      else
		{
		  cart_assert(iold>=1 && iold<=halos->num_halos);
		  iold--;
		  for(j=0, r2old=0.0; j<nDim; j++)
		    {
		      dx = pos[j] - halos->list[iold].pos[j];
		      if(dx < -0.5*num_grid) dx += num_grid;
		      if(dx >  0.5*num_grid) dx -= num_grid;
		      r2old += dx*dx;
		    }
		  if(r2old > r2)
		    {
		      /*
		      //  This cells belongs to a satellite
		      */
		      cell_var(cell,var) = ih + 1;
		    }
		}
	    }
	}

      MESH_RUN_OVER_CELLS_OF_LEVEL_END;

    }

  MESH_RUN_OVER_LEVELS_END;

  cart_free(halo_levels);
}


#ifdef HYDRO
void extDumpLevels(const char *fname, DumpWorker worker, int level1, int level2, halo_list *halos)
{
  const int nout = 10;
  int i, j, ih, size, rank, select;
  MESH_RUN_DECLARE(level,cell);
  FILE *f;
  char str[999];
  float *buffer, *ptr;
  int nselproc, nseltot, ntot = 0;
  
#ifdef RADIATIVE_TRANSFER
  rtStepBegin();
  rtUpdateTables();
#endif

  /*
  //  Count number of cells
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct), reduction(+:ntot)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  ntot++;
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
  
  if(ntot == 0) return;

  cart_debug("Dumping %d cells on levels %d - %d",ntot,level1,level2);

  /*
  //  Create the buffer
  */
  buffer = cart_alloc(float, ntot*nout );

  if(halos != NULL)
    {
      /*
      //  Map cells
      */
      extMapHaloCells(VAR_ACCEL,min_level,halos,1.0);
    }

  /*
  //  Fill in the buffer for each halo
  */
  ih = 0;
  do
    {
      if(halos==NULL || halo_level(&halos->list[ih],MPI_COMM_WORLD)>=level1)
	{
	  nselproc = ntot = 0;
	  MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);

#pragma omp parallel for default(none), private(_Index,cell,ptr,select,i), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,buffer,ntot,halos,ih,worker,units), reduction(+:nselproc)
	  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
	  if(cell_is_leaf(cell))
	    {
	      if(halos != NULL)
		{
		  if(ih+1 == (int)(0.5+cell_var(cell,VAR_ACCEL))) select = 1; else select = 0;
		}
	      else select = 1;
	    }
	  else select = 0;

	  ptr = buffer + nout*(ntot+_Index);

	  if(select)
	    {
	      nselproc++;

	      ptr[0] = units->number_density*cell_gas_density(cell);
	      for(i=1; i<nout; i++) ptr[i] = 0.0;
	      worker(level,cell,nout-1,ptr+1);

	    }
	  else
	    {
	      ptr[0] = -1.0;
	    }
	  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

	  ntot += _Num_level_cells;

	  MESH_RUN_OVER_LEVELS_END;
  
	  MPI_Allreduce(&nselproc,&nseltot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	  if(nseltot > 0)
	    {
	      if(halos != NULL)
		{
		  sprintf(str,"%s.%04d",fname,halos->list[ih].id);
		}
	      else strcpy(str,fname);

	      /*
	      //  Write to a file in order of a proc rank
	      */
	      MPI_Comm_size(MPI_COMM_WORLD,&size);
	      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	      for(i=0; i<size; i++)
		{
		  MPI_Barrier(MPI_COMM_WORLD);
		  if(i == rank)
		    {
		      cart_debug("Writing file piece #%d",i);
		      f = fopen(str,(i==0?"w":"a"));
		      cart_assert(f != NULL);
		      for(j=0; j<ntot; j++)
			{
			  ptr = buffer + nout*j;
			  if(ptr[0] > 0.0)
			    {
#ifdef RADIATIVE_TRANSFER
			      fprintf(f,"%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",ptr[0],ptr[1],ptr[2],ptr[3],ptr[4],ptr[5],ptr[6],ptr[7],ptr[8],ptr[9]);
#else
			      fprintf(f,"%9.3e %9.3e %9.3e %9.3e\n",ptr[0],ptr[1],ptr[2],ptr[3]);
#endif
			    }
			}
		      fclose(f);
		    }
		}
	    }
	}
      ih++;
    }
  while(halos!=NULL && ih<halos->num_halos);

  cart_free(buffer);

}
#endif /* HYDRO */


void extDumpProfiles(const char *fname, DumpWorker worker, int floor_level, float rmin, float rmax, int ndex, halo_list *halos)
{
  const int nout = 10;
  int i, j, ih;
  MESH_RUN_DECLARE(level,cell);
  FILE *f;
  char str[999];
  float *buffer, *gbuffer, *ptr;
  double pos[3], dx, r2, uRad2;
  int ntot, ibin;
  float *rbin, lrmin;
  
#ifdef RADIATIVE_TRANSFER
  rtStepBegin();
  rtUpdateTables();
#endif

  cart_assert(ndex>0 && rmax>0.0 && rmin>0.0);

  if(halos == NULL)
    {
      cart_debug("No halo file is loaded. Skipping dumping profiles.");
      return;
    }

  lrmin = log10(rmin);
  ntot = (int)(0.5+(log10(rmax)-lrmin)*ndex);
  cart_assert(ntot > 0);

  cart_debug("Dumping profiles with %d bins...",ntot);

  rbin = cart_alloc(float,ntot);
  for(i=0; i<ntot; i++) rbin[i] = pow(10.0,lrmin+(i+0.5)/ndex); 

  /*
  //  Create the buffer
  */
  buffer = cart_alloc(float, ntot*nout );
  gbuffer = cart_alloc(float, ntot*nout );

  /*
  //  Map cells
  */
  extMapHaloCells(VAR_ACCEL,min_level,halos,1.0);

  /*
  //  Fill in the buffer for each halo
  */
  for(ih=0; ih<halos->num_halos; ih++) if(halo_level(&halos->list[ih],MPI_COMM_WORLD) >= floor_level)
    {

      cart_debug("Analysing halo #%d...",halos->list[ih].id);

      ptr = gbuffer;

      for(i=0; i<ntot; i++)
	{
	  for(j=0; j<nout; j++)
	    {
	      buffer[j+nout*i] = 0.0;
	    }
	}

      MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);

      uRad2 = pow(units->length/constants->kpc,2.0);

      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
      if(cell_is_leaf(cell) && ih+1==(int)(0.5+cell_var(cell,VAR_ACCEL)))
	{
	  cell_position_double(cell,pos);
	  for(j=0, r2=0.0; j<nDim; j++)
	    {
	      dx = pos[j] - halos->list[ih].pos[j];
	      if(dx < -0.5*num_grid) dx += num_grid;
	      if(dx >  0.5*num_grid) dx -= num_grid;
	      r2 += dx*dx;
	    }
	  
	  ibin = (int)(0.5+(0.5*log10(1.0e-35+uRad2*r2)-lrmin)*ndex);
	  if(ibin>=0 && ibin<ntot)
	    {
	      ptr[0] = units->mass/constants->Msun*cell_gas_density(cell)*cell_volume[level];
	      for(i=1; i<nout; i++) ptr[i] = 0.0;
	      worker(level,cell,nout-1,ptr+1);

	      buffer[ibin*nout] += ptr[0];
	      for(j=1; j<nout; j++)
		{
		  buffer[j+ibin*nout] += ptr[0]*ptr[j];
		}
	    }
	}
      MESH_RUN_OVER_CELLS_OF_LEVEL_END;
      MESH_RUN_OVER_LEVELS_END;
  
      MPI_Allreduce(buffer,gbuffer,ntot*nout,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

      for(i=0; i<ntot; i++)
	{
	  for(j=1; j<nout; j++) if(gbuffer[nout*i] > 0.0)
	    {
	      gbuffer[j+nout*i] /= gbuffer[nout*i];
	    }
	}

      if(local_proc_id == MASTER_NODE)
	{
	  sprintf(str,"%s.%04d",fname,halos->list[ih].id);
	  f = fopen(str,"w");
	  cart_assert(f != NULL);
	  for(i=0; i<ntot; i++)
	    {
	      ptr = gbuffer + nout*i;
#ifdef RADIATIVE_TRANSFER
	      fprintf(f,"%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",rbin[i],ptr[0],ptr[1],ptr[2],ptr[3],ptr[4],ptr[5],ptr[6],ptr[7],ptr[8],ptr[9]);
#else
	      fprintf(f,"%9.3e %9.3e %9.3e %9.3e %9.3e\n",rbin[i],ptr[0],ptr[1],ptr[2],ptr[3]);
#endif
	    }
	  fclose(f);
	}
    }

  cart_free(buffer);
  cart_free(gbuffer);

}


#if defined(PARTICLES) && defined(STARFORM)
void extStarFormationLaw(const char *fname, float spatial_scale, float time_scale, float stellar_age_limit, halo_list *halos)
{
  int i, j, k, ih, level, size, rank, select;
  int num_level_cells, *level_cells, *index, cell;
  float *sfr, *mst, uLen, uDen, uTime, uRate, dt, d2g;
  double dx, r2, r2Max;
  FILE *f;
  char str[999];
  int nselproc, nseltot;
  int level_halo, level_global;

  time_scale *= 1.0e6;  /* turn Myr into years */
  stellar_age_limit *= 1.0e6;

  uLen = units->length/constants->kpc; /* phys pc */
  level = nearest_int(-log(spatial_scale/uLen)/log(2.0));
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Using level: %d",level);
    }

  cart_assert(level>=min_level && level<=max_level);

  level_global = max_level_now_global(MPI_COMM_WORLD);

  /*
  // Units
  */
  uDen = units->density*pow(constants->kpc,3.0)/constants->Msun;    /* Msun/kpc^3 */
  uTime = units->time/constants->yr;         /* yr */
  uRate = uDen/time_scale;                     /* Msun/kpc^3/yr */

  /*
  //  Prepare forward and backward indicies
  */
  select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
  index = cart_alloc(int, num_cells );

#pragma omp parallel for default(none), private(i), shared(index)
  for(i=0; i<num_cells; i++)
    {
      index[i] = -1;
    }

#pragma omp parallel for default(none), private(i), shared(index,level_cells,num_level_cells)
  for(i=0; i<num_level_cells; i++)
    {
      index[level_cells[i]] = i;
    }

  /*
  //  Prepare the SFR buffer array
  */
  sfr = cart_alloc(float, num_level_cells );
  mst = cart_alloc(float, num_level_cells );

  ih = 0;
  do
    {
      if(halos==NULL || (level_halo=halo_level(&halos->list[ih],MPI_COMM_WORLD))>=level)
	{
	  if(halos != NULL)
	    {
	      r2Max = halos->list[ih].rhalo;
	      r2Max *= r2Max;
	    }
	  else
	    {
	      r2Max = (double)num_grid*num_grid;
	    }

#pragma omp parallel for default(none), private(i), shared(num_level_cells,sfr,mst)
	  for(i=0; i<num_level_cells; i++)
	    {
	      sfr[i] = mst[i] = 0.0;
	    }

	  /*
	  //  Measure time-averaged SFR
	  */
	  for(j=0; j<num_particles; j++) if(particle_level[j]>=level && particle_is_star(j))
	    {
	      if(halos != NULL)
		{
		  for(k=0, r2=0.0; k<nDim; k++)
		    {
		      dx = particle_x[j][k] - halos->list[ih].pos[k];
		      if(dx < -0.5*num_grid) dx += num_grid;
		      if(dx >  0.5*num_grid) dx -= num_grid;
		      r2 += dx*dx;
		    }
		  if(r2 < r2Max) select = 1; else select = 0;
		}
	      else select = 1;

	      if(select)
		{
		  cell = cell_find_position_level(level,particle_x[j]);
		  if(cell > -1)
		    {
		      cart_assert(index[cell]>=0 && index[cell]<num_level_cells);
		      dt = uTime*(tl[level]-star_tbirth[j]);
		      if(dt < time_scale)
			{
			  sfr[index[cell]] += particle_mass[j];
			}
		      if(dt < stellar_age_limit)
			{
			  mst[index[cell]] += particle_mass[j];
			}
		    }
		}
	    }

	  /*
	  //  Turn mass into SFR density
	  */
	  nselproc = 0;
#pragma omp parallel for default(none), private(i), shared(num_level_cells,sfr,mst,cell_volume,level), reduction(+:nselproc)
	  for(i=0; i<num_level_cells; i++)
	    {
	      sfr[i] /= cell_volume[level];
	      if(mst[i] > 1.0e-30) nselproc++;
	    }

	  MPI_Allreduce(&nselproc,&nseltot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	  if(nseltot > 0)
	    {
	      if(halos != NULL)
		{
		  sprintf(str,"%s.%04d",fname,halos->list[ih].id);
		}
	      else strcpy(str,fname);

	      /*
	      //  Write to a file in order of a proc rank
	      */
	      MPI_Comm_size(MPI_COMM_WORLD,&size);
	      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	      for(i=0; i<size; i++)
		{
		  MPI_Barrier(MPI_COMM_WORLD);
		  if(i == rank)
		    {
		      cart_debug("Writing file piece #%d",i);
		      f = fopen(str,(i==0?"w":"a"));
		      cart_assert(f != NULL);

		      /*
		      //  Manifest
		      */
		      if(i == 0)
			{
			  if(halos != NULL)
			    {
			      fprintf(f,"%d %9.3e %9.3e %9.3e %d\n",level,uLen*cell_size[level],spatial_scale,time_scale,level_halo);
			    }
			  else
			    {
			      fprintf(f,"%d %9.3e %9.3e %9.3e %d\n",level,uLen*cell_size[level],spatial_scale,time_scale,level_global);
			    }
			}

		      for(j=0; j<num_level_cells; j++) if(mst[j] > 1.0e-30)
			{
#ifdef RADIATIVE_TRANSFER
			  d2g = rtDustToGas(level_cells[j]);
			  fprintf(f,"%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",uRate*sfr[j],uDen*cell_gas_density(level_cells[j]),uDen*cell_HI_density(level_cells[j]),uDen*cell_H2_density(level_cells[j]),units->mass/constants->Msun*mst[j],uDen*d2g*cell_gas_density(level_cells[j]),uDen*d2g*cell_gas_density(level_cells[j])*(cell_HI_fraction(level_cells[j])+cell_H2_fraction(level_cells[j])));
#else
			  fprintf(f,"%9.3e %9.3e %9.3e %9.3e\n",uRate*sfr[j],uDen*cell_gas_density(level_cells[j]),units->mass/constants->Msun*mst[j],uDen*cell_gas_metal_density(level_cells[j])/constants->Zsun);
#endif
			}
		      fclose(f);
		    }
		}
	    }
	}
      ih++;
    }
  while(halos!=NULL && ih<halos->num_halos);

  cart_free(sfr);
  cart_free(mst);
  cart_free(index);
  cart_free(level_cells);

}
#endif


#if defined (HYDRO) && defined(STARFORM)
void extStarFormationLaw2(const char *fname, float spatial_scale, halo_list *halos)
{
  int i, j, k, ih, level, ll, size, rank, select;
  int num_level_cells, *level_cells, cell;
  int num_ll_cells, *ll_cells, *index;
  float *sfr, *sfr1, uLen, uDen, uTime, uRate, d2g;
  double dx, r2, r2Max;
  double pos[3];
  FILE *f;
  char str[999];
  int nselproc, nseltot;

  uLen = units->length/constants->kpc; /* phys pc */
  level = nearest_int(-log(spatial_scale/uLen)/log(2.0));
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Using level: %d",level);
    }

  cart_assert(level>=min_level && level<=max_level);

  /*
  // Units
  */
  uDen = units->density*pow(constants->kpc,3.0)/constants->Msun;    /* Msun/kpc^3 */
  uTime = units->time/constants->yr;      /* yr */
  uRate = uDen/uTime;                     /* Msun/kpc^3/yr */

  /*
  //  Prepare forward and backward indicies
  */
  select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
  index = cart_alloc(int, num_cells );

#pragma omp parallel for default(none), private(i), shared(index)
  for(i=0; i<num_cells; i++)
    {
      index[i] = -1;
    }

#pragma omp parallel for default(none), private(i), shared(index,level_cells,num_level_cells)
  for(i=0; i<num_level_cells; i++)
    {
      index[level_cells[i]] = i;
    }

  /*
  //  Prepare the SFR buffer array
  */
  sfr = cart_alloc(float, num_level_cells );

  ih = 0;
  do
    {
      if(halos==NULL || halo_level(&halos->list[ih],MPI_COMM_WORLD)>=level)
	{
	  if(halos != NULL)
	    {
	      r2Max = halos->list[ih].rhalo;
	      r2Max *= r2Max;
	    }
	  else
	    {
	      r2Max = (double)num_grid*num_grid;
	    }

#pragma omp parallel for default(none), private(i), shared(num_level_cells,sfr)
	  for(i=0; i<num_level_cells; i++)
	    {
	      sfr[i] = 0.0;
	    }

	  /*
	  //  Measure time-averaged SFR
	  */
	  for(ll=max_level_now(); ll>=level; ll--)
	    {
	      if(ll > level)
		{
		  select_level(ll,CELL_TYPE_LOCAL,&num_ll_cells,&ll_cells);
		}
	      else
		{
		  num_ll_cells = num_level_cells;
		  ll_cells = level_cells;
		}

	      sfr1 = cart_alloc(float,num_ll_cells);

	      star_formation_rate(ll,num_ll_cells,ll_cells,sfr1);

	      for(i=0; i<num_ll_cells; i++) if(cell_is_leaf(cell = ll_cells[i]))
		{
		  /*
		  // volume-weigh the SFR
		  */
		  sfr1[i] *= cell_volume[ll];

		  cell_position_double(cell,pos);

		  if(halos != NULL)
		    {
		      for(k=0, r2=0.0; k<nDim; k++)
			{
			  dx = pos[k] - halos->list[ih].pos[k];
			  if(dx < -0.5*num_grid) dx += num_grid;
			  if(dx >  0.5*num_grid) dx -= num_grid;
			  r2 += dx*dx;
			}
		      if(r2 < r2Max) select = 1; else select = 0;
		    }
		  else select = 1;

		  if(select)
		    {
		      while(cell_level(cell) > level)
			{
			  cell = cell_parent_cell(cell);
			  cart_assert(cell > -1);
			}
		      cart_assert(cell_level(cell)==level && index[cell]>=0 && index[cell]<num_level_cells);

		      sfr[index[cell]] += sfr1[i];
		    }
		}

	      cart_free(sfr1);
	      if(ll > level) cart_free(ll_cells);
	    }

	  /*
	  //  Finish volume-weighing of SFR density
	  */
	  nselproc = 0;
#pragma omp parallel for default(none), private(i), shared(num_level_cells,sfr,cell_volume,cell_vars,level_cells,level), reduction(+:nselproc)
	  for(i=0; i<num_level_cells; i++)
	    {
	      sfr[i] /= cell_volume[level];
	      if(sfr[i] > 1.0e-30) nselproc++;
	    }

	  MPI_Allreduce(&nselproc,&nseltot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	  if(nseltot > 0)
	    {
	      if(halos != NULL)
		{
		  sprintf(str,"%s.%04d",fname,halos->list[ih].id);
		}
	      else strcpy(str,fname);

	      /*
	      //  Write to a file in order of a proc rank
	      */
	      MPI_Comm_size(MPI_COMM_WORLD,&size);
	      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	      for(i=0; i<size; i++)
		{
		  MPI_Barrier(MPI_COMM_WORLD);
		  if(i == rank)
		    {
		      cart_debug("Writing file piece #%d",i);
		      f = fopen(str,(i==0?"w":"a"));
		      cart_assert(f != NULL);

		      /*
		      //  Manifest
		      */
		      if(i == 0)
			{
			  if(halos != NULL)
			    {
			      fprintf(f,"%d %9.3e %9.3e %9.3e %d\n",level,uLen*cell_size[level],spatial_scale,0.0,halo_level(&halos->list[ih],MPI_COMM_WORLD));
			    }
			  else
			    {
			      fprintf(f,"%d %9.3e %9.3e %9.3e %d\n",level,uLen*cell_size[level],spatial_scale,0.0,max_level_now_global(MPI_COMM_WORLD));
			    }
			}

		      for(j=0; j<num_level_cells; j++) if(sfr[j] > 1.0e-30)
			{
#ifdef RADIATIVE_TRANSFER
			  d2g = rtDustToGas(level_cells[j]);
			  fprintf(f,"%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",uRate*sfr[j],uDen*cell_gas_density(level_cells[j]),uDen*cell_HI_density(level_cells[j]),uDen*cell_H2_density(level_cells[j]),uDen*d2g*cell_gas_density(level_cells[j]),uDen*d2g*cell_gas_density(level_cells[j])*(cell_HI_fraction(level_cells[j])+cell_H2_fraction(level_cells[j])));
#else
			  fprintf(f,"%9.3e %9.3e %9.3e\n",uRate*sfr[j],uDen*cell_gas_density(level_cells[j]),uDen*cell_gas_metal_density(level_cells[j])/constants->Zsun);
#endif
			}
		      fclose(f);
		    }
		}
	    }
	}
      ih++;
    }
  while(halos!=NULL && ih<halos->num_halos);

  cart_free(sfr);
  cart_free(index);
  cart_free(level_cells);

}
#endif

#if defined (HYDRO) && defined(RADIATIVE_TRANSFER)
/*
//  Compute column densities and compare them with the Sobolev-like approximations
*/
void extCheckSobolevApproximations(const char *fname, int floor_level, int nside, double len)
{
  const float cs[] = { 1.0e-19, 1.0e-20, 1.0e-21, 1.0e-22 };
  const int ncs = sizeof(cs)/sizeof(float);

  int vars[2] = { RT_HVAR_OFFSET+0, RT_HVAR_OFFSET+5 };
  int npix = 12*nside*nside;
  int ipix, ics;
  double tau[ncs];
  double pos[3];
  FILE *f;
  MESH_RUN_DECLARE(level,cell);

  int i, nb[num_neighbors];
  float s1, s2, d;

  cdData *data;

  /*
  // Works only serially so far
  */
  cart_assert(num_procs == 1);

  /*
  // Open output file
  */
  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(fname,"w");
      if(f == NULL)
        {
          cart_error("Unable to open output file.");
        }
      fprintf(f,"# length... cdSob1... cdSob2... tau19.... tau20.... tau21.... tau22....\n");
      fprintf(f,"# chik      cm^{-2} ...\n");
      fprintf(f,"#\n");
    }

  len *= (constants->kpc/units->length);

  /*
  //  Prepare buffers for the LOS traversal
  */
  data = cart_alloc(cdData,npix);

  /*
  //  Loop over all cells of floor_level
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,floor_level,floor_level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

  cell_position_double(cell,pos);

  cdTraverseSky(2,vars,nside,pos,len,floor_level,data);

  cell_all_neighbors(cell,nb);

  /*
  //  Length factors
  */
  s1 = 0.0;
  for(i=0; i<nDim; i++)
    {
      d = cell_gas_density(nb[2*i+1]) - cell_gas_density(nb[2*i]);
      s1 += d*d;
    }
  s1 = cell_size[level]*max(0.0,cell_gas_density(cell)/(1.0e-30+sqrt(s1)));

  s2 = 0.0;
  for(i=0; i<num_neighbors; i++)
    {
      d = cell_gas_density(nb[i]) - cell_gas_density(cell);
      s2 += d*d;
    }
  s2 = 0.5*cell_size[level]*max(0.0,cell_gas_density(cell)/(1.0e-30+sqrt(s2/2)));

  for(ics=0; ics<ncs; ics++)
    {
      d = 0.0;
      for(ipix=0; ipix<npix; ipix++)
	{
	  d += exp(-cs[ics]*units->number_density*units->length*(data[ipix].val[0]+2*data[ipix].val[1]));
	}
      tau[ics] = -log(1.0e-35+d/npix);
    }

  if(local_proc_id == MASTER_NODE)
    {
      fprintf(f,"%9.3le %9.3le %9.3le %9.3le %9.3le %9.3le %9.3le\n",units->length/constants->kpc*data[0].len,units->number_density*units->length*(cell_HI_density(cell)+2*cell_H2_density(cell))*s1,units->number_density*units->length*(cell_HI_density(cell)+2*cell_H2_density(cell))*s2,tau[0],tau[1],tau[2],tau[3]);
    }

  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
 
  if(local_proc_id == MASTER_NODE)
    {
      fclose(f);
    }

  cart_free(data);
}
#endif  /* HYDRO && RADIATIVE_TRANSFER */




