#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER)

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "auxiliary.h"
#include "iterators.h"
#include "logging.h"
#include "parallel.h"
#include "particle.h"
#include "rt_global.h"
#include "rt_transfer.h"
#include "rt_utilities.h"
#include "starformation.h"
#include "timestep.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "F/frt_c.h"


#ifdef RT_VAR_SOURCE
struct rtGlobalValue rtAvgSource;
#endif

#ifdef RT_VAR_OT_FIELD
struct rtGlobalValue rtAvgOTField;
#endif

struct rtGlobalValue rtAvgRF[rt_num_frequencies];


void rtTransferSplitUpdate(int level);
void rtTransferUpdateFields(int nvar, int var0, struct rtGlobalValue *avg, int level_begin, int level_end, MPI_Comm com); 
void rtSetGlobalAbsorption(struct rtGlobalValue *abcAvg, struct rtGlobalValue *abcMax);

#ifndef RT_VAR_SOURCE
void rtTransferUpdateUniformSource(struct rtGlobalValue *avg, MPI_Comm com);
#endif


#ifdef RT_SINGLE_SOURCE
float rtSingleSourceVal;
double rtSingleSourcePos[nDim];

#ifdef RT_VAR_SOURCE
void rtTransferAssignSingleSourceDensity(int level);
#endif
#endif


#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
#include "rt_otvet.h"
#endif


int rtIsThereWork()
{
#ifdef RT_TEST
  return 1;
#endif

#if defined(PARTICLES) && defined(STARFORM)

  return 1;

#else

  return 0;

#endif /* PARTICLES && STARFORM */
}


void rtInitRunTransfer()
{
#ifdef RT_SINGLE_SOURCE
  int i;
#endif
  frt_intg val = rt_num_frequencies;

  start_time(WORK_TIMER);

  frtCall(initruntransfer)(&val);

  end_time(WORK_TIMER);

#ifdef RT_SINGLE_SOURCE
  rtSingleSourceVal = 0.0;
  for(i=0; i<nDim; i++) rtSingleSourcePos[i] = 0.5*num_grid;
#endif

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

  rtInitRunTransferOtvet();

#endif
}


void rtStepBeginTransfer()
{
  int i;
  struct rtGlobalValue abcAvg[rt_num_frequencies], abcMax[rt_num_frequencies];
  frt_real frtAbcAvg[rt_num_frequencies];

  rtSetGlobalAbsorption(abcAvg,abcMax);

  start_time(WORK_TIMER);

  for(i=0; i<rt_num_frequencies; i++) frtAbcAvg[i] = abcAvg[i].Value;
  frtCall(stepbegintransfer)(frtAbcAvg);

  end_time(WORK_TIMER);

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

  rtStepBeginTransferOtvet(abcMax);

#endif

  /*
  //  Set all global values
  */
#ifdef RT_VAR_SOURCE
  rtGlobalValueInit(&rtAvgSource,0.0);
#endif

#ifdef RT_VAR_OT_FIELD
  rtGlobalValueInit(&rtAvgOTField,0.0);
#endif

  for(i=0; i<rt_num_frequencies; i++) rtGlobalValueInit(&rtAvgRF[i],0.0);

  /*
  //  No need to do that here, it is done in the first Superclass::UpdateTables() call
  */
  //  rtGlobalUpdateTransfer(min_level,MPI_COMM_WORLD);
}


void rtStepEndTransfer()
{
  int i;
  struct rtGlobalValue abcAvg[rt_num_frequencies], abcMax[rt_num_frequencies];
  frt_real frtAbcAvg[rt_num_frequencies];

  rtSetGlobalAbsorption(abcAvg,abcMax);

  start_time(WORK_TIMER);

  for(i=0; i<rt_num_frequencies; i++) frtAbcAvg[i] = abcAvg[i].Value;
  frtCall(stependtransfer)(frtAbcAvg);

  end_time(WORK_TIMER);
}


void rtLevelUpdateTransfer(int level)
{
  if(!rtIsThereWork()) return;

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

  rtLevelUpdateTransferOtvet(level);

#endif

  rtTransferSplitUpdate(level);
}


void rtGlobalUpdateTransfer(int top_level, MPI_Comm level_com)
{
  int iomp, j;
  MESH_RUN_DECLARE(level,cell);
#ifdef _OPENMP
  int nomp = omp_get_num_threads();
#else
  int nomp = 1;
#endif
  double sum[nomp][rt_num_frequencies+2];

  start_time(WORK_TIMER);

  /*
  //  Compute per-level averages
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,top_level,max_level);

  /*
  //  Because the reduction variable cannot be an array in C, doing
  //  reduction manually. Cannot re-arrange the loops because of the
  //  cache access pattern.
  */
  for(iomp=0; iomp<nomp; iomp++)
    {
      for(j=0; j<rt_num_frequencies+2; j++) sum[iomp][j] = 0.0;
    }

#pragma omp parallel for default(none), private(_Index,cell,iomp,j), shared(_Num_level_cells,_Level_cells,level,cell_vars,cell_child_oct,sum)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
#ifdef _OPENMP
      iomp = omp_get_thread_num();
#else
      iomp = 0;
#endif

      for(j=0; j<rt_num_frequencies; j++)
	{
	  sum[iomp][j] += cell_var(cell,rt_freq_offset+j);
	}

#ifdef RT_VAR_SOURCE
      sum[iomp][rt_num_frequencies+0] += cell_var(cell,RT_VAR_SOURCE);
#endif /* RT_VAR_SOURCE */

#ifdef RT_VAR_OT_FIELD
      sum[iomp][rt_num_frequencies+1] += cell_var(cell,RT_VAR_OT_FIELD);
#endif /* RT_VAR_OT_FIELD */
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

#ifdef _OPENMP
  for(iomp=1; iomp<nomp; iomp++)
    {
      for(j=0; j<rt_num_frequencies+2; j++) sum[0][j] += sum[iomp][j];
    }
#endif

  for(j=0; j<rt_num_frequencies; j++)
    {
      rtGlobalValueChange(&rtAvgRF[j],level,sum[0][j]);
    }

#ifdef RT_VAR_SOURCE
  rtGlobalValueChange(&rtAvgSource,level,sum[0][rt_num_frequencies+0]);
#endif /* RT_VAR_SOURCE */

#ifdef RT_VAR_OT_FIELD
  rtGlobalValueChange(&rtAvgSource,level,sum[0][rt_num_frequencies+1]);
#endif /* RT_VAR_OT_FIELD */
  
  MESH_RUN_OVER_LEVELS_END;

  end_time(WORK_TIMER);

  for(j=0; j<rt_num_frequencies; j++)
    {
      rtGlobalValueUpdate(&rtAvgRF[j],top_level,max_level,MPI_SUM,level_com);
    }

#ifdef RT_VAR_SOURCE
  rtGlobalValueUpdate(&rtAvgSource,top_level,max_level,MPI_SUM,level_com);
#endif /* RT_VAR_SOURCE */

#ifdef RT_VAR_OT_FIELD
  rtGlobalValueUpdate(&rtAvgSource,top_level,max_level,MPI_SUM,level_com);
#endif /* RT_VAR_OT_FIELD */
}


void rtAfterAssignDensityTransfer(int level, int num_level_cells, int *level_cells)
{
#ifdef RT_VAR_SOURCE
  int i;

#ifdef PARTICLES

  /*
  // If we have a source field that was set inside density(...), 
  // turn the mass per cell into density.
  */
  start_time( WORK_TIMER );

#pragma omp parallel for default(none), private(i), shared(level,num_level_cells,level_cells,cell_vars,cell_volume_inverse)
  for(i=0; i<num_level_cells; i++)
    {
      cell_rt_source(level_cells[i]) *= cell_volume_inverse[level];
    }

  end_time( WORK_TIMER );

#else  /* PARTICLES */

#ifdef RT_SINGLE_SOURCE
  /*
  // Set the source field from a single source
  */
  rtTransferAssignSingleSourceDensity(level);
#else  /* RT_SINGLE_SOURCE */

#error "Invalid set of switches: either PARTICLES or RT_SINGLE_SOURCE must be defined."

#endif /* RT_SINGLE_SOURCE */

#endif /* PARTICLES */

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
  /* Need only local cells */
  rtAfterAssignDensityTransferOtvet(level,num_cells_per_level[level],level_cells);
#endif

#endif /* RT_VAR_SOURCE */
}


void rtTransferSplitUpdate(int level)
{
  int i, j, k;
  int icell;
  int num_level_cells;
  int *level_cells;
  int children[num_children];
  double new_var;
  const double factor = ((double)(1.0/(1<<nDim)));

  start_time(WORK_TIMER);

  if(level < max_level)
    {
      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
#pragma omp parallel for default(none), private(i,icell,j,k,children,new_var), shared(num_level_cells,level_cells,cell_child_oct,cell_vars)
      for(i=0; i<num_level_cells; i++)
	{
	  icell = level_cells[i];
	  if(cell_is_refined(icell))
	    {
	      /*
	      // Average over children
	      */
	      cell_all_children(icell,children);
	      for(j=0; j<rt_num_frequencies; j++)
		{
		  new_var = 0.0;
		  for(k=0; k<num_children; k++)
		    {
		      new_var += cell_var(children[k],rt_freq_offset+j);
		    }
		  cell_var(icell,rt_freq_offset+j) = new_var*factor; 
		}
	    }
	}
      cart_free(level_cells);
    }

  end_time(WORK_TIMER);

}


/*
//  Computes the absorption coefficient at a single frequency bin 
//  ifreq on all cells from a supplied array and sends the result 
//  into abc[0] and (optionally) abc[1].
*/
void rtComputeAbsLevel(int level, int ncells, int *cells, int ifreq, float **abc)
{
  int i, cell;
  frt_real Zsol, dx;
  frt_real buffer[5], out[2];

  /* turn ifreq into a fortran index */
  ifreq++;

  dx = cell_size[level];

#pragma omp parallel for default(none), private(cell,i,Zsol,buffer,out), shared(ncells,cells,ifreq,abc,cell_vars,constants,dx)
  for(i=0; i<ncells; i++)
    {
      cell = cells[i];
#if defined(RT_UV) && defined(ENRICH)
      Zsol = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
      Zsol = 0.0;
#endif /* RT_UV && ENRICH */

      if(sizeof(frt_real) != sizeof(float))  /* Optimization */
	{
	  buffer[0] = cell_gas_density(cell);
	  buffer[1] = cell_HI_density(cell);
	  buffer[2] = cell_HeI_density(cell);
	  buffer[3] = cell_HeII_density(cell);
	  buffer[4] = cell_H2_density(cell);
#if (RT_CFI == 1)
	  frtCall(transfercomputecellabs)(&ifreq,&Zsol,buffer+0,buffer+1,buffer+2,buffer+3,buffer+4,&dx,out+0,out+1);
	  abc[0][i] = out[0];
	  abc[1][i] = out[1];
#else
	  frtCall(transfercomputecellabs)(&ifreq,&Zsol,buffer+0,buffer+1,buffer+2,buffer+3,buffer+4,&dx,out+0);
	  abc[0][i] = out[0];
#endif
	}
      else
	{
#if (RT_CFI == 1)
	  frtCall(transfercomputecellabs)(&ifreq,&Zsol,&(cell_gas_density(cell)),&(cell_HI_density(cell)),&(cell_HeI_density(cell)),&(cell_HeII_density(cell)),&(cell_H2_density(cell)),&dx,abc[0]+i,abc[1]+i);
#else
	  frtCall(transfercomputecellabs)(&ifreq,&Zsol,&(cell_gas_density(cell)),&(cell_HI_density(cell)),&(cell_HeI_density(cell)),&(cell_HeII_density(cell)),&(cell_H2_density(cell)),&dx,abc[0]+i);
#endif
	}

#ifdef RT_DEBUG
      if(abc[0][i]<0.0 || isnan(abc[0][i]))
	{
	  cart_debug("Oops: %d %d %g",i,cell,abc[0][i]);
	  cart_debug("Z:    %g",Zsol);
	  cart_debug("rho:  %g",cell_gas_density(cell));
	  cart_debug("Hi:   %g",cell_HI_density(cell));
	  cart_debug("HeI:  %g",cell_HeI_density(cell));
	  cart_debug("HeII: %g",cell_HeII_density(cell));
	  cart_debug("H2:   %g",cell_H2_density(cell));
	  cart_error("Negative absorption");
	}
#endif
    }
}


void rtSetGlobalAbsorption(struct rtGlobalValue *abcAvg, struct rtGlobalValue *abcMax)
{
  const int n = rt_num_frequencies/2;
  int i, ifield;
  struct rtGlobalValue tmp[2*n];
  double s, s0, s1;
  float amin, amax;
  MESH_RUN_DECLARE(level,cell);
  float *abc[2], *abc1, w;
#ifdef RT_DEBUG
  int j;
#endif

  start_time(WORK_TIMER);

  for(i=0; i<n; i++)
    {
      rtGlobalValueInit(&tmp[2*i+0],0.0);
      rtGlobalValueInit(&tmp[2*i+1],0.0);
      rtGlobalValueInit(&abcAvg[i],0.0);
      rtGlobalValueInit(&abcMax[i],0.0);
    }
  
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);

  abc[0] = cart_alloc(float,_Num_level_cells);
#if (RT_CFI == 1)
  abc[1] = cart_alloc(float,_Num_level_cells);
#else
  abc[1] = 0;
#endif

  for(i=0; i<n; i++)
    {
      /*
      //  Average by weighting with the far field only
      */
      ifield = rt_freq_offset + rt_num_frequencies/2 + i;

      rtComputeAbsLevel(level,_Num_level_cells,_Level_cells,i,abc);
#if (RT_CFI == 1)
      abc1 = abc[1];
#else
      abc1 = abc[0];
#endif

      rtuGetLinearArrayMaxMin(_Num_level_cells,abc1,&amax,&amin);
      rtGlobalValueChange(&abcMax[i],level,amax);

      s = s0 = s1 = 0.0;

#pragma omp parallel for default(none), private(_Index,cell,i,w), shared(_Num_level_cells,_Level_cells,level,n,cell_vars,cell_child_oct,ifield,abc1), reduction(+:s,s0,s1)
      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
      if(cell_is_leaf(cell))
	{
	  w = cell_var(cell,ifield);
	  s0 += w*cell_volume[level];
	  s1 += w*abc1[_Index]*cell_volume[level];
	  s += abc1[_Index]*cell_volume[level];

#ifdef RT_DEBUG
	  if(w<0.0 || isnan(w)) 
	    {
	      cart_debug("Oops: %d %d %d %d %g",i,ifield,_Index,cell,w);
	      for(j=0; j<num_vars; j++)
		{
		  cart_debug("Var: %d %g",j,cell_var(cell,j));
		}
	      cart_error("Negative radiation field");
	    }
#endif
	}
      MESH_RUN_OVER_CELLS_OF_LEVEL_END;

      rtGlobalValueChange(&tmp[2*i+0],level,s0);
      rtGlobalValueChange(&tmp[2*i+1],level,s1);
      rtGlobalValueChange(&abcAvg[i],level,s);

    }

  cart_free(abc[0]);
#if (RT_CFI == 1)
  cart_free(abc[1]);
#endif

  MESH_RUN_OVER_LEVELS_END;

  end_time(WORK_TIMER);

  rtGlobalValueUpdate(&tmp[2*i+0],min_level,max_level,MPI_SUM,MPI_COMM_WORLD);
  rtGlobalValueUpdate(&tmp[2*i+1],min_level,max_level,MPI_SUM,MPI_COMM_WORLD);
  rtGlobalValueUpdate(&abcAvg[i],min_level,max_level,MPI_SUM,MPI_COMM_WORLD);
  rtGlobalValueUpdate(&abcMax[i],min_level,max_level,MPI_MAX,MPI_COMM_WORLD);

  start_time(WORK_TIMER);

  /*
  //  Weighted average
  */
  for(i=0; i<n; i++)
    {
      if(tmp[2*i+0].Value > 1.0e-35)
	{
	  abcAvg[i].Value = tmp[2*i+1].Value/tmp[2*i+0].Value;
	}
    }

  end_time(WORK_TIMER);
}


#if defined(RT_SINGLE_SOURCE) && defined(RT_VAR_SOURCE)

void rtTransferAssignSingleSourceDensity(int level)
{
  int icell;
  double corner[nDim];
  double size2;
  float mass;
  double cornerx0, cornerx1, cornery0, cornery1, cornerz0, cornerz1;
  double x, y, z;
  double xs, ys, zs;
  double dx0, dx1, dy0, dy1, dz0, dz1;
  double d00, d01, d10, d11;

  icell = cell_find_position_level( level, rtSingleSourcePos );
  if(icell==-1 || cell_level(icell)<level) return;

  start_time( WORK_TIMER );

  size2 = 0.5*cell_size[level];

  x = rtSingleSourcePos[0];
  y = rtSingleSourcePos[1];
  z = rtSingleSourcePos[2];

  cornerx0 = x - size2;
  cornerx1 = x + size2;
  cornery0 = y - size2;
  cornery1 = y + size2;
  cornerz0 = z - size2;
  cornerz1 = z + size2;
		
  if ( cornerx0 < 0.0 ) cornerx0 += (double)num_grid;
  if ( cornerx1 >= (double)num_grid ) cornerx1 -= (double)num_grid;
  if ( cornery0 < 0.0 ) cornery0 += (double)num_grid;
  if ( cornery1 >= (double)num_grid ) cornery1 -= (double)num_grid;
  if ( cornerz0 < 0.0 ) cornerz0 += (double)num_grid;
  if ( cornerz1 >= (double)num_grid ) cornerz1 -= (double)num_grid;

  xs = x*cell_size_inverse[level] + 0.5;
  ys = y*cell_size_inverse[level] + 0.5;
  zs = z*cell_size_inverse[level] + 0.5;

  dx1 = xs - floor(xs);
  dy1 = ys - floor(ys);
  dz1 = zs - floor(zs);

  dx0 = 1.0 - dx1;
  dy0 = 1.0 - dy1;
  dz0 = 1.0 - dz1;

  dx0 *= rtSingleSourceVal*cell_volume_inverse[level];
  dx1 *= rtSingleSourceVal*cell_volume_inverse[level];

  d00 = dx0*dy0;
  d01 = dx0*dy1;
  d10 = dx1*dy0;
  d11 = dx1*dy1;

  /* child 0 */
  corner[0] = cornerx0;
  corner[1] = cornery0;
  corner[2] = cornerz0;
  
  icell = cell_find_position_level( level, corner );
  if ( icell != -1 ) {
    mass = d00*dz0;
    cell_rt_source(icell) += mass;
  }

  /* child 1 */
  corner[0] = cornerx1;

  icell = cell_find_position_level( level, corner );
  if ( icell != -1 ) {
    mass = d10*dz0;
    cell_rt_source(icell) += mass;
  }

  /* child 2 */
  corner[0] = cornerx0;
  corner[1] = cornery1;

  icell = cell_find_position_level( level, corner );
  if ( icell != -1 ) {
    mass = d01*dz0;
    cell_rt_source(icell) += mass;
  }

  /* child 3 */
  corner[0] = cornerx1;

  icell = cell_find_position_level( level, corner );
  if ( icell != -1 ) {
    mass = d11*dz0;
    cell_rt_source(icell) += mass;
  }

  /* child 4 */
  corner[0] = cornerx0;
  corner[1] = cornery0;
  corner[2] = cornerz1;
  
  icell = cell_find_position_level( level, corner );
  if ( icell != -1 ) {
    mass = d00*dz1;
    cell_rt_source(icell) += mass;
  }

  /* child 5 */
  corner[0] = cornerx1;

  icell = cell_find_position_level( level, corner );
  if ( icell != -1 ) {
    mass = d10*dz1;
    cell_rt_source(icell) += mass;
  }

  /* child 6 */
  corner[0] = cornerx0;
  corner[1] = cornery1;

  icell = cell_find_position_level( level, corner );
  if ( icell != -1 ) {
    mass = d01*dz1;
    cell_rt_source(icell) += mass;
  }

  /* child 7 */
  corner[0] = cornerx1;

  icell = cell_find_position_level( level, corner );
  if ( icell != -1 ) {
    mass = d11*dz1;
    cell_rt_source(icell) += mass;
  }

  end_time( WORK_TIMER );

}

#endif /* RT_SINGLE_SOURCE && RT_VAR_SOURCE */
#endif /* RADIATIVE_TRANSFER && RT_TRANSFER */

