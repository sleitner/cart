#include "defs.h"

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "auxiliary.h"
#include "particle.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"


#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#include "rt_utilities.h"
#include "F/frt_parameters.ch"
#ifdef RT_DEBUG
#include "rt_debug.h"
#endif
#endif



#ifdef RADIATIVE_TRANSFER
void rtGetSobolevFactors(int cell, int level, float *len, float *vel);
void extDumpChemicalState(int level1, int level2)
{
  const int nout = 8;
  int i, j, size, rank;
  MESH_RUN_DECLARE(level,cell);
  float uDen, uLen, uColumn;
  FILE *f;
  float *buffer, *ptr, soblen, sobvel, rate[IRATE_DIM];
  int ntot = 0;

  rtStepBegin();

  /*
  //  Count number of cells
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_child_oct), reduction(+:ntot)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell)) ntot++;
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
  
  cart_debug("Dumping %d cells",ntot);

  /*
  //  Create the buffer
  */
  buffer = (float *)cart_alloc(ntot*nout*sizeof(float));

  /*
  //  Fill in the buffer
  */
  ntot = 0;
  MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);

  uDen = 1.123e-5*Omega0*hubble*hubble/pow(aexp[level],3.0);
  uLen = 3.086e24*r0/hubble*aexp[level];

#pragma omp parallel for default(none), private(_Index,cell,soblen,sobvel,ptr,rate), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,uDen,uLen,buffer,ntot)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
      ptr = buffer + nout*(ntot+_Index);

#ifdef RT_CHEMISTRY
      rtGetSobolevFactors(cell,level,&soblen,&sobvel);
#else
      soblen = 0.0;
#endif

      rtGetPhotoRates(cell,rate);

      ptr[0] = uDen*cell_gas_density(cell);
      ptr[1] = rtTemInK(cell);
      ptr[2] = ptr[0]*uLen*soblen;
      ptr[3] = cell_gas_metallicity(cell)/(0.02*cell_gas_density(cell));
      ptr[4] = cell_HI_fraction(cell);
      ptr[5] = cell_HII_fraction(cell);
      ptr[6] = cell_H2_fraction(cell);
      ptr[7] = rate[12]*1.1e10;  /* UV field at 1000A in units of Draine field */
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

  ntot += _Num_level_cells;

  MESH_RUN_OVER_LEVELS_END;
  
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
	  f = fopen("chemistry-state.res",(i==0?"w":"a"));
	  cart_assert(f != NULL);
	  for(j=0; j<ntot; j++)
	    {
	      ptr = buffer + nout*j;
	      fprintf(f,"%9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e %9.3e\n",ptr[0],ptr[1],ptr[2],ptr[3],ptr[4],ptr[5],ptr[6],ptr[7]);
	    }
	  fclose(f);
	}
    }

  cart_free(buffer);

}
#endif


#if defined(PARTICLES) && defined(STARFORM)
void extDumpKennicuttLaw(float spatialScale, float timeScale)
{
  int i, j, level, size, rank;
  int num_level_cells, *level_cells, *index, cell;
  float *sfr, uLen, uDen, uTime, uRate, dt;
  FILE *f;

  timeScale *= 1.0e6;  /* turn Myr into years */

  uLen = r0*aexp[0]/hubble*1.0e3; /* phys pc */
  level = nearest_int(-log(spatialScale/uLen)/log(2.0));
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Using level: %d",level);
    }

  cart_assert(level>=min_level && level<=max_level);

  /*
  // Units
  */
  uLen = r0*aexp[level]/hubble*1.0e3;         /* phys kpc */
  uDen = rho0/pow(aexp[level],3.0)*1.0e-9;    /* Msun/kpc^3 */
  uTime = t0*aexp[level]*aexp[level];         /* yr */
  uRate = uDen/timeScale;                     /* Msun/kpc^3/yr */

  /*
  //  Prepare forward and backward indicies
  */
  select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
  index = cart_alloc(num_cells*sizeof(int));

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
  sfr = (float *)cart_alloc(num_level_cells*sizeof(float));

#pragma omp parallel for default(none), private(i), shared(num_level_cells,sfr)
  for(i=0; i<num_level_cells; i++)
    {
      sfr[i] = 0.0;
    }

  /*
  //  Measure time-averaged SFR
  */
  for(j=0; j<num_local_star_particles; j++)
    {
      dt = uTime*(tl[level]-star_tbirth[j]);
      if(dt < timeScale)
	{
	  cell = cell_find_position_level(level,particle_x[j]);
	  if(cell > -1)
	    {
	      cart_assert(index[cell]>=0 && index[cell]<num_level_cells);
	      sfr[index[cell]] += particle_mass[j];
	    }
	}
    }

  /*
  //  Turn mass into SFR density
  */
#pragma omp parallel for default(none), private(i), shared(num_level_cells,sfr,cell_volume,level)
  for(i=0; i<num_level_cells; i++)
    {
      sfr[i] /= cell_volume[level];
    }

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
	  f = fopen("kennicut.res",(i==0?"w":"a"));
	  cart_assert(f != NULL);

	  /*
	  //  Manifest
	  */
	  if(i == 0)
	    {
	      fprintf(f,"%d %9.3e %9.3e %9.3e\n",level,uLen*cell_size[level],spatialScale,timeScale);
	    }

	  for(j=0; j<num_level_cells; j++) if(sfr[j]*uRate > 1.0e-20)
	    {
#ifdef RADIATIVE_TRANSFER
	      fprintf(f,"%9.3e %9.3e %9.3e %9.3e\n",uRate*sfr[j],uDen*cell_gas_density(level_cells[j]),uDen*cell_HI_density(level_cells[j]),uDen*cell_H2_density(level_cells[j]));
#else
	      fprintf(f,"%9.3e %9.3e\n",uRate*sfr[j],uDen*cell_gas_density(level_cells[j]));
#endif
	    }
	  fclose(f);
	}
    }

  cart_free(sfr);
  cart_free(index);
  cart_free(level_cells);

}
#endif

