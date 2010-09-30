#include "config.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

#ifdef RADIATIVE_TRANSFER
#include "F/frt_c.h"
#endif

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
void extDumpLevels(const char *fname, int nout, DumpWorker worker, int level1, int level2, halo_list *halos)
{
  int i, j, ih, node, size, rank, select;
  MESH_RUN_DECLARE(level,cell);
  FILE *f;
  char str[999];
  float *buffer, *ptr;
  int nselproc, nseltot, ntot = 0;
  int *selected;
  
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
  selected = cart_alloc(int, ntot );

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

#pragma omp parallel for default(none), private(_Index,cell,ptr,select,i), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,buffer,ntot,halos,ih,worker,units,selected,nout), reduction(+:nselproc)
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

	  selected[ntot+_Index] = select;

	  if(select)
	    {
	      nselproc++;
	      ptr = buffer + nout*(ntot+_Index);
	      for(i=0; i<nout; i++) ptr[i] = 0.0;
	      worker(level,cell,nout,ptr);

	    }
	  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

	  ntot += _Num_level_cells;

	  MESH_RUN_OVER_LEVELS_END;
  
	  MPI_Allreduce(&nselproc,&nseltot,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	  if(nseltot > 0)
	    {
	      if(halos != NULL)
		{
		  sprintf(str,"%s.%05d",fname,halos->list[ih].id);
		}
	      else strcpy(str,fname);

	      /*
	      //  Write to a file in order of a proc rank
	      */
	      MPI_Comm_size(MPI_COMM_WORLD,&size);
	      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	      for(node=0; node<size; node++)
		{
		  MPI_Barrier(MPI_COMM_WORLD);
		  if(node == rank)
		    {
		      cart_debug("Writing file piece #%d",node);
		      f = fopen(str,(node==0?"w":"a"));
		      cart_assert(f != NULL);
		      for(j=0; j<ntot; j++)
			{
			  ptr = buffer + nout*j;
			  if(selected[j])
			    {
			      fprintf(f,"%9.3e",ptr[0]);
			      for(i=1; i<nout; i++) fprintf(f," %9.3e",ptr[i]);
			      fprintf(f,"\n");
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

  cart_free(selected);
  cart_free(buffer);
}


void extDumpProfiles(const char *fname, int nout, DumpWorker worker, int floor_level, float rmin, float rmax, int ndex, halo_list *halos)
{
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
	      for(i=0; i<nout; i++) ptr[i] = 0.0;
	      worker(level,cell,nout,ptr);

	      ptr[0] = units->mass/constants->Msun*cell_gas_density(cell)*cell_volume[level];

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
	  sprintf(str,"%s.%05d",fname,halos->list[ih].id);
	  f = fopen(str,"w");
	  cart_assert(f != NULL);
	  for(i=0; i<ntot; i++)
	    {
	      ptr = gbuffer + nout*i;

              fprintf(f,"%9.3e ",pow(10.0,lrmin+(i+0.5)/ndex));
              for(j=0; j<nout; j++) fprintf(f," %9.3e",ptr[j]);
	      fprintf(f,"\n");
	    }
	  fclose(f);
	}
    }

  cart_free(buffer);
  cart_free(gbuffer);

}
#endif /* HYDRO */


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
		  sprintf(str,"%s.%05d",fname,halos->list[ih].id);
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
void extStarFormationLaw2(const char *fname, float spatial_scale, const halo_list *halos)
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
		  sprintf(str,"%s.%05d",fname,halos->list[ih].id);
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


#if defined (HYDRO) && defined(STARFORM) && defined(RADIATIVE_TRANSFER)

/*
//  Dump ISRF and a variable var[cell] with a hierarchy of levels
*/
void extRFvsSFR1(const char *froot, int top_level, float *var, const halo *h)
{
  const char *fext[] = { "sfr", "den", "dmw", "fh2" };
  const int nfiles = sizeof(fext)/sizeof(char*);
  MESH_RUN_DECLARE(level,cell);
  float rate[frtRATE_DIM], uLen;
  int i, j, l, parent, size, rank, save;
  double dx, pos[nDim], r = 0.0;
  FILE *f[nfiles];
  char str[999], fsuffix[99];

  cart_assert(top_level>=min_level && top_level<=max_level);

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /*
  // Units
  */
  uLen = units->length/constants->kpc; /* phys pc */

  /*
  //  Create a set of files for each requested level
  */
  for(level=max_level_now_global(MPI_COMM_WORLD); level>=top_level; level--)
    {
      cart_debug("Working on level %d...",level);

      /*
      //  Write to a file in order of a proc rank
      */
      if(h != NULL)
	{
	  if(halo_level(h,MPI_COMM_WORLD) < level) continue;
	  sprintf(fsuffix,"L=%02d.res.%05d",level,h->id);
	}
      else
	{
	  sprintf(fsuffix,"L=%02d.res",level);
	}

      select_level(level,CELL_TYPE_LOCAL,&_Num_level_cells,&_Level_cells);

      for(i=0; i<size; i++)
	{
	  MPI_Barrier(MPI_COMM_WORLD);
	  if(i == rank)
	    {
	      cart_debug("Writing file piece #%d",i);

	      for(j=0; j<nfiles; j++)
		{
		  sprintf(str,"%s-%s.%s",froot,fext[j],fsuffix);
		  f[j] = fopen(str,(i==0?"w":"a"));
		  cart_assert(f[j] != NULL);
		}

	      /*
	      //  Manifest
	      */
	      if(i == 0)
		{
		  for(j=0; j<nfiles; j++) fprintf(f[j],"%d %9.3e\n",level,uLen*cell_size[level]);
		}

	      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

	      if(cell_is_leaf(cell))
		{
		  save = 0;
		  for(parent=cell,l=level; l>=min_level; l--)
		    {
		      cart_assert(parent > -1);
		      if(var[parent] > 0.0) save = 1;
		      parent = cell_parent_cell(parent);
		    }

		  if(save && h!=NULL)
		    {
		      cell_position_double(cell,pos);
		      for(j=0, r=0.0; j<nDim; j++)
			{
			  dx = pos[j] - h->pos[j];
			  if(dx < -0.5*num_grid) dx += num_grid;
			  if(dx >  0.5*num_grid) dx -= num_grid;
			  r += dx*dx;
			}
		      r = sqrt(r);
		      if(r > h->rhalo) save = 0;
		    }

		  if(save)
		    {
		      rtGetPhotoRates(cell,rate);
		      for(j=0; j<nfiles; j++) fprintf(f[j],"%9.3e",rate[frtRATE_CiLW]*1.05e10);
		      for(parent=cell,l=level; l>=min_level; l--)
			{
			  if(nfiles > 0) fprintf(f[0]," %9.3e",var[parent]);
			  if(nfiles > 1) fprintf(f[1]," %9.3e",units->number_density*cell_gas_density(parent));
			  if(nfiles > 2) fprintf(f[2]," %9.3e",rtDustToGas(parent));
			  if(nfiles > 3) fprintf(f[3]," %9.3e",cell_H2_fraction(parent));
			  parent = cell_parent_cell(parent);
			}
		      if(h != NULL)
			{
			  for(j=0; j<nfiles; j++) fprintf(f[j]," %9.3le",uLen*r);
			}
		      for(j=0; j<nfiles; j++) fprintf(f[j],"\n");
		    }
		}
	      MESH_RUN_OVER_CELLS_OF_LEVEL_END;

	      for(j=0; j<nfiles; j++) fclose(f[j]);
	    }
	}

      cart_free(_Level_cells);
      _Level_cells = 0;
    }
}


void extRFvsSFR(const char *froot, int top_level, const halo_list *halos)
{
  MESH_RUN_DECLARE(level,cell);
  int j, ih;
  float *var, *sfr;
  float uLen, uDen, uTime, uRate;

  rtStepBegin();
  rtUpdateTables();

  /*
  // Units
  */
  uLen = units->length/constants->kpc; /* phys pc */
  uDen = units->density*pow(constants->kpc,3.0)/constants->Msun;    /* Msun/kpc^3 */
  uTime = units->time/constants->yr;      /* yr */
  uRate = uDen/uTime;                     /* Msun/kpc^3/yr */

  var = cart_alloc(float,num_cells);
  memset(var,0,sizeof(float)*num_cells);

  /*
  //  Fill in the data array
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,_MaxLevel,min_level);

  /*
  //  Prepare the SFR buffer array
  */
  sfr = cart_alloc(float,_Num_level_cells);
  star_formation_rate(level,_Num_level_cells,_Level_cells,sfr);

#pragma omp parallel for default(none), private(_Index,cell,j), shared(_Num_level_cells,_Level_cells,sfr,var,uRate,cell_child_oct)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

  if(cell_is_leaf(cell))
    {
      var[cell] = uRate*max(0.0,sfr[_Index]);
    }
  else
    {
      var[cell] = 0.0;
      for(j=0; j<num_children; j++)
	{
	  var[cell] += var[cell_child(cell,j)]/num_children;
	}
    }

  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  
  cart_free(sfr);

  MESH_RUN_OVER_LEVELS_END;

  if(halos == NULL)
    {
      extRFvsSFR1(froot,top_level,var,NULL);
    }
  else
    {
      for(ih=0; ih<halos->num_halos; ih++) if(halo_level(&halos->list[ih],MPI_COMM_WORLD) >= top_level)
	{
	  extRFvsSFR1(froot,top_level,var,&halos->list[ih]);
	}
    }

  cart_free(var);
}

#endif
