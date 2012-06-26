#include "config.h"

#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "hydro.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rt.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "halo_tools.h"
#include "ism.h"


/*
//  Return the value of the conserved cell variable var but only for the temperature
//  range Tmin < T < Tmax.
*/
#ifdef HYDRO
float extGetCellVarInTRange(int cell, int var, float Tmin, float Tmax)
{
  int j;
  float w;

  if(cell_is_refined(cell))
    {
      w = 0.0;
      for(j=0; j<num_children; j++) w += extGetCellVarInTRange(cell_child(cell,j),var,Tmin,Tmax);
      return w/num_children;
    }
  else
    {
      w = cell_gas_temperature(cell)*units->temperature;
      if(w>Tmin && w<Tmax) return cell_var(cell,var); else return 0.0; 
    }
}
#endif /* HYDRO */


void extDumpLevels(const char *fname, int nout, const DumpWorker *workers, int level1, int level2, struct HALO_LIST *halos)
{
  int i, j, ih, numh, node, size, rank, select;
  MESH_RUN_DECLARE(level,cell);
  FILE *f;
  char str[999];
  float *buffer, *ptr;
  int nselproc, nseltot, ntot = 0;
  int *selected;

  if(nout == 0)
    {
      cart_debug("There are no cell data to dump. Skipping DumpLevels...");
      return;
    }

  cart_assert(nout>0 && nout<100);
  cart_assert(workers != NULL);
  for(i=0; i<nout; i++)
    {
      cart_assert(workers[i].Value != NULL);
      cart_assert(workers[i].Header != NULL);
    }

  /*
  //  Count number of cells
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct), reduction(+:ntot)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  ntot++;
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
  
  if(ntot == 0)
    {
      cart_debug("There are no cells on levels %d to %d",level1,level2);
      return;
    }

  cart_debug("Dumping %d cells on levels %d - %d",ntot,level1,level2);

  /*
  //  Create the buffer
  */
  buffer = cart_alloc(float, ntot*nout );
  selected = cart_alloc(int, ntot );

#ifdef COSMOLOGY
  if(halos!=NULL && halos->map==NULL)
    {
      /*
      //  Map cells
      */
      map_halos(min_level,halos,1.0);
      numh = halos->num_halos;
    }
  else numh = 1;
#else
  halos = NULL;
  numh = 1;
#endif

  /*
  //  Fill in the buffer for each halo
  */
  for(ih=0; ih<numh; ih++)
#ifdef COSMOLOGY
    if(halos==NULL || halo_level(&halos->list[ih],mpi.comm.run)>=level1)
#endif
      {
	nselproc = ntot = 0;
	MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);

#pragma omp parallel for default(none), private(_Index,cell,ptr,select,i), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,buffer,ntot,halos,ih,workers,units,selected,nout,level2), reduction(+:nselproc)
	MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
	if(level==level2 || cell_is_leaf(cell))
	  {
#ifdef COSMOLOGY
	    if(halos != NULL)
	      {
		if(ih+1 == halos->map[cell]) select = 1; else select = 0;
	      }
	    else
#endif
	      {
		select = 1;
	      }
	  }
	else select = 0;

	selected[ntot+_Index] = select;

	if(select)
	  {
	    nselproc++;
	    ptr = buffer + nout*(ntot+_Index);
	    for(i=0; i<nout; i++)
	      {
#ifdef COSMOLOGY
		if(halos != NULL)
		  {
		    ptr[i] = workers[i].Value(level,cell,halos->list[ih].pos,halos->list[ih].vel);
		  }
		else
#endif
		  {
		    ptr[i] = workers[i].Value(level,cell,NULL,NULL);
		  }
	      }
	  }
	MESH_RUN_OVER_CELLS_OF_LEVEL_END;

	ntot += _Num_level_cells;

	MESH_RUN_OVER_LEVELS_END;
  
	MPI_Allreduce(&nselproc,&nseltot,1,MPI_INT,MPI_SUM,mpi.comm.run);
	if(nseltot > 0)
	  {
#ifdef COSMOLOGY
	    if(halos != NULL)
	      {
		sprintf(str,"%s.%05d",fname,halos->list[ih].id);
	      }
	    else
#endif
	      {
		strcpy(str,fname);
	      }

	    /*
	    //  Write to a file in order of a proc rank
	    */
	    MPI_Comm_size(mpi.comm.run,&size);
	    MPI_Comm_rank(mpi.comm.run,&rank);
	    for(node=0; node<size; node++)
	      {
		MPI_Barrier(mpi.comm.run);
		if(node == rank)
		  {
		    cart_debug("Writing file piece #%d",node);
		    f = fopen(str,(node==0?"w":"a"));
		    cart_assert(f != NULL);

		    if(node == 0)
		      {
			fprintf(f,"# Levels: from %d to %d\n",level1,level2);
			fprintf(f,"# Columns:\n");
			for(i=0; i<nout; i++) fprintf(f,"#   %2d. %s\n",i+1,workers[i].Header);
			fprintf(f,"#\n");
		      }

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

  cart_free(selected);
  cart_free(buffer);
}


void extDumpLevelsLowMemory(const char *fname, int nout, const DumpWorker *workers, int level1, int level2, struct HALO_LIST *halos)
{
  int i, ih, numh, node, size, rank, select;
  MESH_RUN_DECLARE(level,cell);
  FILE *f;
  char str[999];
  float *ptr;
  int nselproc, nseltot, ntot = 0;

  if(nout == 0)
    {
      cart_debug("There are no cell data to dump. Skipping DumpLevels...");
      return;
    }

  cart_assert(nout>0 && nout<100);
  cart_assert(workers != NULL);
  for(i=0; i<nout; i++)
    {
      cart_assert(workers[i].Value != NULL);
      cart_assert(workers[i].Header != NULL);
    }

  /*
  //  Count number of cells
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct), reduction(+:ntot)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  ntot++;
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
  
  if(ntot == 0)
    {
      cart_debug("There are no cells on levels %d to %d",level1,level2);
      return;
    }

  cart_debug("Dumping %d cells on levels %d - %d",ntot,level1,level2);

  /*
  //  Create the buffer
  */
  ptr = cart_alloc(float, nout );

#ifdef COSMOLOGY
  if(halos!=NULL && halos->map==NULL)
    {
      /*
      //  Map cells
      */
      map_halos(min_level,halos,1.0);
      numh = halos->num_halos;
    }
  else numh = 1;
#else
  halos = NULL;
  numh = 1;
#endif

  /*
  //  Fill in the buffer for each halo
  */
  for(ih=0; ih<numh; ih++)
#ifdef COSMOLOGY
    if(halos==NULL || halo_level(&halos->list[ih],mpi.comm.run)>=level1)
#endif
      {
	nselproc = ntot = 0;
	MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);

#pragma omp parallel for default(none), private(_Index,cell,select,i), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,halos,ih,level2), reduction(+:nselproc)
	MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
	if(level==level2 || cell_is_leaf(cell))
	  {
#ifdef COSMOLOGY
	    if(halos != NULL)
	      {
		if(ih+1 == halos->map[cell]) select = 1; else select = 0;
	      }
	    else
#endif
	      {
		select = 1;
	      }
	  }
	else select = 0;

	if(select) nselproc++;

	MESH_RUN_OVER_CELLS_OF_LEVEL_END;
	MESH_RUN_OVER_LEVELS_END;
  
	MPI_Allreduce(&nselproc,&nseltot,1,MPI_INT,MPI_SUM,mpi.comm.run);

	if(nseltot == 0) continue;

#ifdef COSMOLOGY
	if(halos != NULL)
	  {
	    sprintf(str,"%s.%05d",fname,halos->list[ih].id);
	  }
	else
#endif
	  {
	    strcpy(str,fname);
	  }

	/*
	//  Write to a file in order of a proc rank
	*/
	MPI_Comm_size(mpi.comm.run,&size);
	MPI_Comm_rank(mpi.comm.run,&rank);
	for(node=0; node<size; node++)
	  {
	    MPI_Barrier(mpi.comm.run);
	    if(node == rank)
	      {
		cart_debug("Writing file piece #%d",node);
		f = fopen(str,(node==0?"w":"a"));
		cart_assert(f != NULL);

		if(node == 0)
		  {
		    fprintf(f,"# Levels: from %d to %d\n",level1,level2);
		    fprintf(f,"# Columns:\n");
		    for(i=0; i<nout; i++) fprintf(f,"#   %2d. %s\n",i+1,workers[i].Header);
		    fprintf(f,"#\n");
		  }

		MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);
		MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
		if(level==level2 || cell_is_leaf(cell))
		  {
#ifdef COSMOLOGY
		    if(halos != NULL)
		      {
			if(ih+1 == halos->map[cell]) select = 1; else select = 0;
		      }
		    else
#endif
		      {
			select = 1;
		      }
		  }
		else select = 0;

		if(select)
		  {
		    for(i=0; i<nout; i++)
		      {
#ifdef COSMOLOGY
			if(halos != NULL)
			  {
			    ptr[i] = workers[i].Value(level,cell,halos->list[ih].pos,halos->list[ih].vel);
			  }
			else
#endif
			  {
			    ptr[i] = workers[i].Value(level,cell,NULL,NULL);
			  }
		      }

		    fprintf(f,"%9.3e",ptr[0]);
		    for(i=1; i<nout; i++) fprintf(f," %9.3e",ptr[i]);
		    fprintf(f,"\n");
		  }
		MESH_RUN_OVER_CELLS_OF_LEVEL_END;
		MESH_RUN_OVER_LEVELS_END;
  
		fclose(f);
	      }
	  }
      }

  cart_free(ptr);
}


#ifdef COSMOLOGY
void extDumpHaloProfiles(const char *fname, int nout, const DumpWorker *workers, float rmin, float rmax, int ndex, struct HALO_LIST *halos, int resolution_level, float outer_edge)
{
  const int num_weights = 4;
  int i, j, ih;
  MESH_RUN_DECLARE(level,cell);
  FILE *f;
  char str[999];
  float *buffer, *gbuffer, *weight, *gweight, *ptr, w[num_weights];
  double pos[3], r, uRad;
  int ntot, ibin;
  float *rbin, lrmin;
  
  if(nout == 0)
    {
      cart_debug("There are no cell data to dump. Skipping DumpHaloProfiles...");
      return;
    }

  if(halos == NULL)
    {
      cart_debug("No halo file is loaded. Skipping DumpHaloProfiles...");
      return;
    }

  cart_assert(nout>0 && nout<100);
  cart_assert(ndex>0 && rmax>0.0 && rmin>0.0);
  cart_assert(workers != NULL);
  for(i=0; i<nout; i++)
    {
      cart_assert(workers[i].Value != NULL);
      cart_assert(workers[i].Header != NULL);
    }

  for(j=0; j<nout; j++) 
    {
      if(workers[j].WeightId<0 || workers[j].WeightId>=num_weights) cart_error("Invalid weight id[%d] = %d (should be between 0 and %d)",j,workers[j].WeightId,num_weights-1);
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

  weight = cart_alloc(float, ntot*num_weights );
  gweight = cart_alloc(float, ntot*num_weights );

  uRad = units->length/constants->kpc;

  /*
  //  Map cells
  */
  if(halos->map == NULL)
    {
      map_halos(min_level,halos,outer_edge);
    }

  /*
  //  Fill in the buffer for each halo
  */
  for(ih=0; ih<halos->num_halos; ih++) if(halo_level(&halos->list[ih],mpi.comm.run) >= resolution_level)
    {

      cart_debug("Analysing halo #%d...",halos->list[ih].id);

      ptr = gbuffer;

      for(i=0; i<ntot; i++)
	{
	  for(j=0; j<num_weights; j++)
	    {
	      weight[j+i*num_weights] = 0.0;
	    }
	  for(j=0; j<nout; j++)
	    {
	      buffer[j+i*nout] = 0.0;
	    }
	}

      MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);

      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
      if(cell_is_leaf(cell) && ih+1==halos->map[cell])
	{
	  cell_center_position(cell,pos);
	  r = compute_distance_periodic(pos,halos->list[ih].pos);
	  
	  ibin = (int)((log10(1.0e-35+uRad*r)-lrmin)*ndex);
	  if(ibin < 0) ibin = 0;
	  if(ibin>=0 && ibin<ntot)
	    {
	      for(j=0; j<nout; j++)
		{
		  ptr[j] = workers[j].Value(level,cell,halos->list[ih].pos,halos->list[ih].vel);
		}

	      w[0] = cell_volume[level];
	      w[1] = (cell_total_mass(cell)+cell_volume[level]);
#ifdef HYDRO
	      w[2] = cell_gas_density(cell)*cell_volume[level];
#else
	      w[2] = 1;
#endif
#if defined(HYDRO) && defined(RADIATIVE_TRANSFER)
	      w[3] = cell_HI_density(cell)*cell_volume[level];
#else
	      w[3] = w[2];
#endif

	      for(j=0; j<num_weights; j++)
		{
		  weight[j+ibin*num_weights] += w[j];
		}
	      for(j=0; j<nout; j++)
		{
		  buffer[j+ibin*nout] += ptr[j]*w[workers[j].WeightId];
		}
	    }
	}
      MESH_RUN_OVER_CELLS_OF_LEVEL_END;
      MESH_RUN_OVER_LEVELS_END;
  
      MPI_Allreduce(buffer,gbuffer,ntot*nout,MPI_FLOAT,MPI_SUM,mpi.comm.run);
      MPI_Allreduce(weight,gweight,ntot*num_weights,MPI_FLOAT,MPI_SUM,mpi.comm.run);

      for(i=0; i<ntot; i++)
	{
	  for(j=0; j<nout; j++) if(gweight[workers[j].WeightId+i*num_weights] > 0.0)
	    {
	      gbuffer[j+i*nout] /= gweight[workers[j].WeightId+i*num_weights];
	    }
	}

      if(local_proc_id == MASTER_NODE)
	{
	  sprintf(str,"%s.%05d",fname,halos->list[ih].id);
	  f = fopen(str,"w");
	  cart_assert(f != NULL);

	  fprintf(f,"# Columns:\n");
	  for(j=0; j<nout; j++) fprintf(f,"#   %2d. %s\n",j+1,workers[j].Header);
	  fprintf(f,"#\n");

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

  cart_free(weight);
  cart_free(gweight);

  cart_free(buffer);
  cart_free(gbuffer);

}
#endif


void extDumpPointProfile(const char *fname, int nout, const DumpWorker *workers, float rmin, float rmax, int ndex, double center[3])
{
  const int num_weights = 4;
  int i, j;
  MESH_RUN_DECLARE(level,cell);
  FILE *f;
  char str[999];
  float *buffer, *gbuffer, *weight, *gweight, *ptr, w[num_weights];
  double pos[3], r, uRad;
  int ntot, ibin;
  float *rbin, lrmin;
  
  if(nout == 0)
    {
      cart_debug("There are no cell data to dump. Skipping DumpPointProfile...");
      return;
    }

  cart_assert(nout>0 && nout<100);
  cart_assert(ndex>0 && rmax>0.0 && rmin>0.0);
  cart_assert(workers != NULL);
  for(i=0; i<nout; i++)
    {
      cart_assert(workers[i].Value != NULL);
      cart_assert(workers[i].Header != NULL);
    }

  for(j=0; j<nout; j++) 
    {
      if(workers[j].WeightId<0 || workers[j].WeightId>=num_weights) cart_error("Invalid weight id[%d] = %d (should be between 0 and %d)",j,workers[j].WeightId,num_weights-1);
    }

  lrmin = log10(rmin);
  ntot = (int)(0.5+(log10(rmax)-lrmin)*ndex);
  cart_assert(ntot > 0);

  cart_debug("Dumping a profile with %d bins...",ntot);

  rbin = cart_alloc(float,ntot);
  for(i=0; i<ntot; i++) rbin[i] = pow(10.0,lrmin+(i+0.5)/ndex); 

  /*
  //  Create the buffer
  */
  buffer = cart_alloc(float, ntot*nout );
  gbuffer = cart_alloc(float, ntot*nout );

  weight = cart_alloc(float, ntot*num_weights );
  gweight = cart_alloc(float, ntot*num_weights );

  uRad = units->length/constants->kpc;

  /*
  //  Fill in the buffer
  */
  ptr = gbuffer;

  for(i=0; i<ntot; i++)
    {
      for(j=0; j<num_weights; j++)
	{
	  weight[j+i*num_weights] = 0.0;
	}
      for(j=0; j<nout; j++)
	{
	  buffer[j+i*nout] = 0.0;
	}
    }

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);

  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
  if(cell_is_leaf(cell))
    {
      cell_center_position(cell,pos);
      r = compute_distance_periodic(pos,center);
	  
      ibin = (int)((log10(1.0e-35+uRad*r)-lrmin)*ndex);
      if(ibin < 0) ibin = 0;
      if(ibin>=0 && ibin<ntot)
	{
	  for(j=0; j<nout; j++)
	    {
	      ptr[j] = workers[j].Value(level,cell,center,NULL);
	    }

	  w[0] = cell_volume[level];
#ifdef GRAVITY
	  w[1] = (cell_total_mass(cell)+cell_volume[level]);
#else
	  w[1] = 1;
#endif
#ifdef HYDRO
	  w[2] = cell_gas_density(cell)*cell_volume[level];
#else
	  w[2] = 1;
#endif
#if defined(HYDRO) && defined(RADIATIVE_TRANSFER)
	  w[3] = cell_HI_density(cell)*cell_volume[level];
#else
	  w[3] = w[2];
#endif

	  for(j=0; j<num_weights; j++)
	    {
	      weight[j+ibin*num_weights] += w[j];
	    }
	  for(j=0; j<nout; j++)
	    {
	      buffer[j+ibin*nout] += ptr[j]*w[workers[j].WeightId];
	    }
	}
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
  
  MPI_Allreduce(buffer,gbuffer,ntot*nout,MPI_FLOAT,MPI_SUM,mpi.comm.run);
  MPI_Allreduce(weight,gweight,ntot*num_weights,MPI_FLOAT,MPI_SUM,mpi.comm.run);

  for(i=0; i<ntot; i++)
    {
      for(j=0; j<nout; j++) if(gweight[workers[j].WeightId+i*num_weights] > 0.0)
	{
	  gbuffer[j+i*nout] /= gweight[workers[j].WeightId+i*num_weights];
	}
    }

  if(local_proc_id == MASTER_NODE)
    {
      sprintf(str,"%s",fname);
      f = fopen(str,"w");
      cart_assert(f != NULL);

      fprintf(f,"# Columns:\n");
      for(j=0; j<nout; j++) fprintf(f,"#   %2d. %s\n",j+1,workers[j].Header);
      fprintf(f,"#\n");

      for(i=0; i<ntot; i++)
	{
	  ptr = gbuffer + nout*i;

	  fprintf(f,"%9.3e ",pow(10.0,lrmin+(i+0.5)/ndex));
	  for(j=0; j<nout; j++) fprintf(f," %9.3e",ptr[j]);
	  fprintf(f,"\n");
	}
      fclose(f);
    }

  cart_free(weight);
  cart_free(gweight);

  cart_free(buffer);
  cart_free(gbuffer);

}


#if defined(PARTICLES) && defined(STAR_FORMATION)
void extStarFormationLaw(const char *fname, float spatial_scale, float time_scale, float stellar_age_limit, const struct HALO_LIST *halos)
{
  int i, j, ih, level, size, rank, lstarted, gstarted, select;
  int num_level_cells, *level_cells, *index, cell;
  float *sfr, *mst, *zst, uLen, uDen, uTime, uRate, dt;
  double pos[nDim];
  FILE *f;
  char str[999];
  int nselproc, nseltot;
  int level_halo, level_global;

  time_scale *= 1.0e6;  /* turn Myr into years */
  stellar_age_limit *= 1.0e6;

  uLen = units->length/constants->kpc; /* phys kpc */
  level = nearest_int(-log(spatial_scale/uLen)/log(2.0));
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Using level: %d",level);
    }

  cart_assert(level>=min_level && level<=max_level);

  level_global = max_level_now_global(mpi.comm.run);

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

#pragma omp parallel for default(none), private(i), shared(index,size_cell_array)
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
  zst = cart_alloc(float, num_level_cells );

#pragma omp parallel for default(none), private(i), shared(num_level_cells,sfr,mst,zst)
  for(i=0; i<num_level_cells; i++)
    {
      sfr[i] = mst[i] = zst[i] = 0.0;
    }

  /*
  //  Measure time-averaged SFR
  */
  for(j=0; j<num_particles; j++) if(particle_level[j]>=level && particle_is_star(j))
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
#ifdef ENRICHMENT
	      zst[index[cell]] += particle_mass[j]*star_metallicity_II[j];
#ifdef ENRICHMENT_SNIa
	      zst[index[cell]] += particle_mass[j]*star_metallicity_Ia[j];
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
	    }
	}
    }

  /*
  //  Turn mass into SFR density
  */
  nselproc = 0;
#pragma omp parallel for default(none), private(i), shared(num_level_cells,sfr,mst,zst,cell_volume,level), reduction(+:nselproc)
  for(i=0; i<num_level_cells; i++)
    {
      sfr[i] /= cell_volume[level];
      if(mst[i] > 1.0e-30)
	{
	  nselproc++;
	  zst[i] /= mst[i];
	}
    }

  MPI_Allreduce(&nselproc,&nseltot,1,MPI_INT,MPI_SUM,mpi.comm.run);

  if(nseltot > 0)
    {
      ih = -1;
      do
	{
	  if(ih<0 || (halos!=NULL && (level_halo=halo_level(&halos->list[ih],mpi.comm.run))>=level))
	    {
	      if(ih >= 0)
		{
		  cart_assert(halos != NULL);
		  sprintf(str,"%s.%05d",fname,halos->list[ih].id);
		}
	      else
		{
		  strcpy(str,fname);
		}

	      lstarted = 0;

	      /*
	      //  Write to a file in order of a proc rank
	      */
	      MPI_Comm_size(mpi.comm.run,&size);
	      MPI_Comm_rank(mpi.comm.run,&rank);
	      for(i=0; i<size; i++)
		{

		  MPI_Allreduce(&lstarted,&gstarted,1,MPI_INT,MPI_MAX,mpi.comm.run);

		  if(i == rank)
		    {
		      for(j=0; j<num_level_cells; j++) if(mst[j] > 1.0e-30)
			{
			  if(ih >= 0)
			    {
			      cell_center_position(level_cells[j],pos);
			      if(compute_distance_periodic(pos,halos->list[ih].pos) < halos->list[ih].rhalo) select = 1; else select = 0;
			    }
			  else select = 1;

			  if(select)
			    {
			      if(lstarted == 0)
				{
				  if(gstarted == 0)
				    {
				      cart_debug("Writing file %s",str);
				      f = fopen(str,"w");
				    }
				  else
				    {
				      f = fopen(str,"a");
				    }

				  cart_assert(f != NULL);

				  /*
				  //  Manifest
				  */
				  fprintf(f,"# Parameters:\n");
				  fprintf(f,"#   1. Level\n");
				  fprintf(f,"#   2. Actual scale (kpc)\n");
				  fprintf(f,"#   3. Requested scale (kpc)\n");
				  fprintf(f,"#   4. Time scale (Myr)\n");
				  fprintf(f,"#   5. Lowest level reached\n");

				  if(ih >= 0)
				    {
				      fprintf(f,"%2d %9.3e %9.3e %9.3e %d\n",level,uLen*cell_size[level],spatial_scale,time_scale,level_halo);
				    }
				  else
				    {
				      fprintf(f,"%2d %9.3e %9.3e %9.3e %d\n",level,uLen*cell_size[level],spatial_scale,time_scale,level_global);
				    }

				  fprintf(f,"# Columns:\n");
				  fprintf(f,"#   1. SFR density (Msun/yr/kpc^3)\n");
				  fprintf(f,"#   2. Total baryon density (Msun/kpc^3)\n");
				  fprintf(f,"#   3. Gas metallicity (solar units)\n");
#ifdef RADIATIVE_TRANSFER
				  fprintf(f,"#   4. Atomic hydrogen density (Msun/kpc^3)\n");
				  fprintf(f,"#   5. Ionized hydrogen density (Msun/kpc^3)\n");
				  fprintf(f,"#   6. Molecular hydrogen density (Msun/kpc^3)\n");
				  fprintf(f,"#   7. Dust-to-gas ratio (solar units)\n");
				  fprintf(f,"#   8. Instellar radiation field (solar units)\n");
				  fprintf(f,"#   9. Stellar mass (Msun)\n");
				  fprintf(f,"#  10. Stellar metallicity (solar units)\n");
#else
				  fprintf(f,"#   4. Stellar mass (Msun)\n");
				  fprintf(f,"#   5. Stellar metallicity (solar units)\n");
#endif
				  fprintf(f,"#\n");

				  lstarted = 1;
				}

			      fprintf(f,"%9.3e %9.3e %9.3e ",uRate*sfr[j],uDen*cell_gas_density(level_cells[j]),cell_gas_metal_density(level_cells[j])/constants->Zsun/cell_gas_density(level_cells[j]));
#ifdef RADIATIVE_TRANSFER
			      fprintf(f,"%9.3e %9.3e %9.3e %9.3e %9.3e ",uDen*cell_HI_density(level_cells[j]),uDen*cell_HII_density(level_cells[j]),uDen*2*cell_H2_density(level_cells[j]),rtDmw(level_cells[j]),rtUmw(level_cells[j]));
#endif
			      fprintf(f,"%9.3e %9.3e\n",units->mass/constants->Msun*mst[j],zst[j]/constants->Zsun);
			    }
			}

		      if(lstarted == 1) fclose(f);
		    }
		}
	    }
	  ih++;
	}
      while(halos!=NULL && ih<halos->num_halos);
    }

  cart_free(sfr);
  cart_free(mst);
  cart_free(zst);
  cart_free(index);
  cart_free(level_cells);

}
#endif /* PARTICLES && STAR_FORMATION */


#if defined (HYDRO) && defined(STAR_FORMATION)
void extStarFormationLaw2(const char *fname, float spatial_scale, const struct HALO_LIST *halos)
{
  int i, j, ih, level, ll, size, rank, lstarted, gstarted, select;
  int num_level_cells, *level_cells, cell;
  int num_ll_cells, *ll_cells, *index;
  float *sfr, *sfr1, uLen, uDen, uTime, uRate;
  double pos[nDim];
  FILE *f;
  char str[999];
  int nselproc, nseltot;
  int level_halo, level_global;

  uLen = units->length/constants->kpc; /* phys pc */
  level = nearest_int(-log(spatial_scale/uLen)/log(2.0));
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Using level: %d",level);
    }

  cart_assert(level>=min_level && level<=max_level);

  level_global = max_level_now_global(mpi.comm.run);

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

#pragma omp parallel for default(none), private(i), shared(index,size_cell_array)
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

	  while(cell_level(cell) > level)
	    {
	      cell = cell_parent_cell(cell);
	      cart_assert(cell > -1);
	    }
	  cart_assert(cell_level(cell)==level && index[cell]>=0 && index[cell]<num_level_cells);

	  sfr[index[cell]] += sfr1[i];
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

  MPI_Allreduce(&nselproc,&nseltot,1,MPI_INT,MPI_SUM,mpi.comm.run);
  if(nseltot > 0)
    {
      ih = -1;
      do
	{
	  if(ih<0 || (halos!=NULL && (level_halo=halo_level(&halos->list[ih],mpi.comm.run))>=level))
	    {
	      if(ih >= 0)
		{
		  cart_assert(halos != NULL);
		  sprintf(str,"%s.%05d",fname,halos->list[ih].id);
		}
	      else
		{
		  strcpy(str,fname);
		}

	      lstarted = 0;

	      /*
	      //  Write to a file in order of a proc rank
	      */
	      MPI_Comm_size(mpi.comm.run,&size);
	      MPI_Comm_rank(mpi.comm.run,&rank);
	      for(i=0; i<size; i++)
		{

		  MPI_Allreduce(&lstarted,&gstarted,1,MPI_INT,MPI_MAX,mpi.comm.run);

		  if(i == rank)
		    {
		      for(j=0; j<num_level_cells; j++) if(sfr[j] > 1.0e-30)
			{
			  if(ih >= 0)
			    {
			      cell_center_position(level_cells[j],pos);
			      if(compute_distance_periodic(pos,halos->list[ih].pos) < halos->list[ih].rhalo) select = 1; else select = 0;
			    }
			  else select = 1;

			  if(select)
			    {
			      if(lstarted == 0)
				{
				  if(gstarted == 0)
				    {
				      cart_debug("Writing file %s",str);
				      f = fopen(str,"w");
				    }
				  else
				    {
				      f = fopen(str,"a");
				    }

				  cart_assert(f != NULL);

				  /*
				  //  Manifest
				  */
				  fprintf(f,"# Parameters:\n");
				  fprintf(f,"#   1. Level\n");
				  fprintf(f,"#   2. Actual scale (kpc)\n");
				  fprintf(f,"#   3. Requested scale (kpc)\n");
				  fprintf(f,"#   4. Time scale (Myr)\n");
				  fprintf(f,"#   5. Lowest level reached\n");

				  if(ih >= 0)
				    {
				      fprintf(f,"%2d %9.3e %9.3e %9.3e %d\n",level,uLen*cell_size[level],spatial_scale,0.0,level_halo);
				    }
				  else
				    {
				      fprintf(f,"%2d %9.3e %9.3e %9.3e %d\n",level,uLen*cell_size[level],spatial_scale,0.0,level_global);
				    }

				  fprintf(f,"# Columns:\n");
				  fprintf(f,"#   1. SFR density (Msun/yr/kpc^3)\n");
				  fprintf(f,"#   2. Total baryon density (Msun/kpc^3)\n");
				  fprintf(f,"#   3. Gas metallicity (solar units)\n");
#ifdef RADIATIVE_TRANSFER
				  fprintf(f,"#   4. Atomic hydrogen density (Msun/kpc^3)\n");
				  fprintf(f,"#   5. Ionized hydrogen density (Msun/kpc^3)\n");
				  fprintf(f,"#   6. Molecular hydrogen density (Msun/kpc^3)\n");
				  fprintf(f,"#   7. Dust-to-gas ratio (solar units)\n");
				  fprintf(f,"#   8. Instellar radiation field (solar units)\n");
#endif
				  fprintf(f,"#\n");

				  lstarted = 1;
				}

			      fprintf(f,"%9.3e %9.3e %9.3e ",uRate*sfr[j],uDen*cell_gas_density(level_cells[j]),cell_gas_metal_density(level_cells[j])/constants->Zsun/cell_gas_density(level_cells[j]));
#ifdef RADIATIVE_TRANSFER
			      fprintf(f,"%9.3e %9.3e %9.3e %9.3e %9.3e ",uDen*cell_HI_density(level_cells[j]),uDen*cell_HII_density(level_cells[j]),uDen*2*cell_H2_density(level_cells[j]),rtDmw(level_cells[j]),rtUmw(level_cells[j]));
#endif
			      fprintf(f,"\n");
			    }
			}

		      if(lstarted == 1) fclose(f);
		    }
		}
	    }
	  ih++;
	}
      while(halos!=NULL && ih<halos->num_halos);
    }

  cart_free(sfr);
  cart_free(index);
  cart_free(level_cells);

}
#endif /* HYDRO && STAR_FORMATION */


#if defined(PARTICLES) && defined(STAR_FORMATION)
void extHaloStars(const char *fname, const halo *h, float rmax)
{
  int i, j, ntot, rank, size;
  double r, rad, uLen, uMass;
  char str[999];
  FILE *f;

  cart_assert(h != NULL);

  uLen = units->length/constants->kpc; /* phys kpc */
  uMass = units->mass/constants->Msun;

  rad = rmax*h->rhalo;

  /*
  //  Count our particles
  */
  ntot = 0;
  for(j=0; j<num_particles; j++) if(particle_is_star(j))
    {
      if(compute_distance_periodic(particle_x[j],(double *)h->pos) < rad) ntot++;
    }

  sprintf(str,"%s.%05d",fname,h->id);

  /*
  //  Write to a file in order of a proc rank
  */
  MPI_Comm_size(mpi.comm.run,&size);
  MPI_Comm_rank(mpi.comm.run,&rank);
  for(i=0; i<size; i++)
    {
      MPI_Barrier(mpi.comm.run);
      if(i==rank && ntot>0)
	{
	  cart_debug("Writing file piece #%d",i);
	  f = fopen(str,(i==0?"w":"a"));
	  cart_assert(f != NULL);

	  /*
	  //  Manifest
	  */
	  if(i == 0)
	    {
	      fprintf(f,"# Halo id, Mvir, Rvir, Rhalo, tNow(Gyr)\n");
	      fprintf(f,"%d %9.3e %9.3e %9.3e %7.4lf\n",h->id,h->mvir,h->rvir,h->rhalo,
#ifdef COSMOLOGY
		      tphys_from_tcode(tl[min_level])/1.0e9
#else
		      uTime*tl[min_level]
#endif /* COSMOLOGY */
		      );
	      fprintf(f,"#  R(kpc)   M(Msun)  Mi(Msun) tF(Gyr) ZII(Zsun) ZIa(Zsun)\n");
	      fprintf(f,"#\n");
	    }

	  for(j=0; j<num_particles; j++) if(particle_is_star(j))
	    {
	      r = compute_distance_periodic(particle_x[j],(double *)h->pos);
	      if(r < rad)
		{
		  fprintf(f,"%9.3e %9.3e %9.3e %7.4lf",uLen*r,uMass*particle_mass[j],uMass*star_initial_mass[j],
#ifdef COSMOLOGY
			  tphys_from_tcode(star_tbirth[j])/1.0e9
#else
			  uTime*star_tbirth[j]
#endif /* COSMOLOGY */
			  );
#ifdef ENRICHMENT
		  fprintf(f," %9.3e",star_metallicity_II[j]/constants->Zsun);
#ifdef ENRICHMENT_SNIa
		  fprintf(f," %9.3e",star_metallicity_Ia[j]/constants->Zsun);
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
		  fprintf(f,"\n");
		}
	    }
	  fclose(f);
	}
    }
}
#endif /* PARTICLES && STAR_FORMATION */

