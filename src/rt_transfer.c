#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER)

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "auxiliary.h"
#include "hydro.h"
#include "iterators.h"
#include "parallel.h"
#include "rt.h"
#include "rt_global.h"
#include "rt_transfer.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "F/frt_c.h"


struct rtGlobalValue rtAvgRF[rt_num_fields];
struct rtGlobalValue rtAvgACxRF[rt_num_fields];
struct rtGlobalValue rtMaxAC[rt_num_freqs];
struct rtGlobalValue rtAvgAC[rt_num_freqs];

float rtGlobalAC[rt_num_freqs];


#ifndef RT_VAR_SOURCE
void rtTransferUpdateUniformSource(struct rtGlobalValue *avg, MPI_Comm com);
#endif


#ifdef RT_SINGLE_SOURCE
int rtSingleSourceLevel;
float rtSingleSourceValue;
double rtSingleSourcePos[nDim];

#ifdef RT_VAR_SOURCE
void rtTransferAssignSingleSourceDensity(int level);
#endif
#endif


#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
#include "rt_otvet.h"
#endif


void rtInitRunTransfer()
{
  int cell, freq;
  frt_intg nfreq = rt_num_freqs;
  frt_intg ncomp = rt_num_fields_per_freq;
#ifdef RT_SINGLE_SOURCE
  int i;
#endif

  start_time(WORK_TIMER);

  frtCall(initruntransfer)(&nfreq,&ncomp);

  end_time(WORK_TIMER);

#ifdef RT_SINGLE_SOURCE
  rtSingleSourceLevel = min_level;
  rtSingleSourceValue = 0.0;
  for(i=0; i<nDim; i++) rtSingleSourcePos[i] = 0.5*num_grid;
#endif

  for(cell=0; cell<num_cells; cell++)
    {
      for(freq=0; freq<rt_num_freqs; freq++)
	{
	  cell_var(cell,rt_far_field_offset+freq) = 1.0;
	}
    }

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

  rtInitRunTransferOtvet();

#endif
}


void rtInitStepTransfer()
{
  int i;

  /*
  //  Set all global values
  */
  for(i=0; i<rt_num_fields; i++)
    {
      rtGlobalValueInit(&rtAvgRF[i],0.0);
      rtGlobalValueInit(&rtAvgACxRF[i],0.0);
    }
  for(i=0; i<rt_num_freqs; i++)
    {
      rtGlobalValueInit(&rtMaxAC[i],0.0);
      rtGlobalValueInit(&rtAvgAC[i],0.0);
    }

  rtGlobalUpdateTransfer(min_level,mpi.comm.run);

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

  rtInitStepTransferOtvet(rtMaxAC);

#endif
}


/*
//  Update all running averages
*/
void rtGlobalUpdateTransfer(int top_level, MPI_Comm level_com)
{
  int iomp, i, freq, field;
  int level, cell, *level_cells, num_level_cells, bottom_level = max_level_local();
  float amin, amax;
  float *abc[2];
#ifdef _OPENMP
  int nomp = omp_get_max_threads();
#else
  int nomp = 1;
#endif
  double s[nomp][rt_num_fields];
  double s1, sw[nomp][rt_num_fields_per_freq];
  frt_intg idx;
  frt_real fabc[rt_num_fields_per_freq];

  start_time(WORK_TIMER);

  /*
  //  Compute per-level averages
  */
  for(level=top_level; level<=bottom_level; level++)
    {
      select_level_with_condition(1,level,&num_level_cells,&level_cells);
      if(num_level_cells == 0) continue;

      /*
      //  Because the reduction variable cannot be an array in C, doing
      //  reduction manually. Cannot re-arrange the loops because of the
      //  cache access pattern.
      */
      for(i=0; i<nomp; i++)
	{
	  for(field=0; field<rt_num_fields; field++) s[i][field] = 0.0;
	}

#pragma omp parallel for default(none), private(i,field,cell,iomp), shared(num_level_cells,level_cells,level,cell_vars,cell_child_oct,nomp,s)
      for(i=0; i<num_level_cells; i++)
	{
	  cell = level_cells[i]; // No need to check for leaves, we selected only them!

#ifdef _OPENMP
	  iomp = omp_get_thread_num();
	  cart_assert(iomp>=0 && iomp<nomp);
#else
	  iomp = 0;
#endif

	  for(field=0; field<rt_num_fields; field++)
	    {
	      s[iomp][field] += cell_var(cell,rt_field_offset+field)*cell_volume[level]/num_root_cells;
	    }
	}

#ifdef _OPENMP
      for(i=1; i<nomp; i++)
	{
	  for(field=0; field<rt_num_fields; field++) s[0][field] += s[i][field];
	}
#endif

      for(field=0; field<rt_num_fields; field++)
	{
	  rtGlobalValueUpdate(&rtAvgRF[field],level,s[0][field]);
	}

      /*
      //  Now do absoprtion - since we need to recompute the abs. coefficient,
      //  loop over frequencies first
      */
      abc[0] = cart_alloc(float,num_level_cells);
#if (RT_CFI == 1)
      abc[1] = cart_alloc(float,num_level_cells);
#else
      abc[1] = abc[0];
#endif

      for(freq=0; freq<rt_num_freqs; freq++)
	{
	  /*
	  //  Average by weighting with the far field only
	  */
	  rtComputeAbsLevel(level,num_level_cells,level_cells,freq,abc);

	  linear_array_maxmin(num_level_cells,abc[1],&amax,&amin);
	  rtGlobalValueUpdate(&rtMaxAC[freq],level,amax);

	  /*
	  //  Because the reduction variable cannot be an array in C, doing
	  //  reduction manually. Cannot re-arrange the loops because of the
	  //  cache access pattern.
	  */
	  for(i=0; i<nomp; i++)
	    {
	      for(field=0; field<rt_num_fields_per_freq; field++) sw[i][field] = 0.0;
	    }

#pragma omp parallel for default(none), private(cell,i,iomp,field), shared(num_level_cells,level_cells,abc,level,cell_vars,freq,nomp,sw), reduction(+:s1)
	  for(i=0; i<num_level_cells; i++)
	    {
	      cell = level_cells[i]; // No need to check for leaves, we selected only them!

#ifdef _OPENMP
	      iomp = omp_get_thread_num();
	      cart_assert(iomp>=0 && iomp<nomp);
#else
	      iomp = 0;
#endif

	      for(field=0; field<rt_num_fields_per_freq; field++)
		{
		  sw[iomp][field] += cell_var(cell,rt_field_offset+rt_num_freqs*field+freq)*abc[1][i]*cell_volume[level]/num_root_cells;
		}
	      s1 += 1*abc[1][i]*cell_volume[level]/num_root_cells;
	    }

#ifdef _OPENMP
	  for(i=1; i<nomp; i++)
	    {
	      for(field=0; field<rt_num_fields_per_freq; field++)
		{
		  sw[0][field] += sw[i][field];
		}
	    }
#endif
	  rtGlobalValueUpdate(&rtAvgAC[freq],level,s1);
	  for(field=0; field<rt_num_fields_per_freq; field++) rtGlobalValueUpdate(&rtAvgACxRF[rt_num_freqs*field+freq],level,sw[0][field]);
	}

      cart_free(abc[0]);
#if (RT_CFI == 1)
      cart_free(abc[1]);
#endif

      cart_free(level_cells);
    }

  end_time(WORK_TIMER);

  for(field=0; field<rt_num_fields; field++)
    {
      rtGlobalValueCommunicate(&rtAvgRF[field],MPI_SUM,level_com);
      rtGlobalValueCommunicate(&rtAvgACxRF[field],MPI_SUM,level_com);
    }

  for(freq=0; freq<rt_num_freqs; freq++)
    {
      rtGlobalValueCommunicate(&rtMaxAC[freq],MPI_MAX,level_com);
      rtGlobalValueCommunicate(&rtAvgAC[freq],MPI_SUM,level_com);
    }

  start_time(WORK_TIMER);

  /*
  //  Weighted average
  */
  for(freq=0; freq<rt_num_freqs; freq++)
    {
      idx = freq + 1;

      for(field=0; field<rt_num_fields_per_freq; field++)
	{
	  if(rtAvgRF[rt_num_freqs*field+freq].Value > 1.0e-35)
	    {
	      fabc[field] = rtAvgACxRF[rt_num_freqs*field+freq].Value/rtAvgRF[rt_num_freqs*field+freq].Value;
	    }
	  else
	    {
	      fabc[field] = rtAvgAC[freq].Value;
	    }
	}
      
      rtGlobalAC[freq] = frtCall(transferglobalac)(&idx,fabc);
    }

  end_time(WORK_TIMER);


#ifdef RT_OUTPUT
  for(freq=0; freq<rt_num_freqs; freq++)
    {
      cart_debug("RT: Global abs[%d] = %10.3le, avg = %10.3le, max = %10.3le",freq,rtGlobalAC[freq],rtAvgAC[freq].Value,rtMaxAC[freq].Value);
    }
  for(field=0; field<rt_num_fields; field++)
    {
      cart_debug("RT: field=%d: <rf>=%10.3e, <abc>=%10.3e",field,rtAvgRF[field].Value,(rtAvgRF[field].Value>0.0)?rtAvgACxRF[field].Value/rtAvgRF[field].Value:0.0);
    }
#endif /* RT_OUTPUT */


  for(cell=0; cell<num_cells; cell++)
    {
      for(freq=0; freq<rt_num_freqs; freq++)
        {
//	  if(rtAvgRF[rt_num_freqs*rt_num_near_fields_per_freq+freq].Value > 0.0) cell_var(cell,rt_far_field_offset+freq) /= rtAvgRF[rt_num_freqs*rt_num_near_fields_per_freq+freq].Value;
        }
    }


#ifdef RT_SINGLE_SOURCE
  start_time(WORK_TIMER);
  cell = cell_find_position(rtSingleSourcePos);
  if(cell>-1 && cell_is_local(cell))
    {
      level = cell_level(cell);
    }
  else
    {
      level = -1;
    }
  end_time(WORK_TIMER);
 
  start_time(COMMUNICATION_TIMER);
  /*
  //  NG: I don't know why, but Bcast blocks here, hence using Allreduce
  */
  MPI_Allreduce(&level,&rtSingleSourceLevel,1,MPI_INT,MPI_MAX,level_com);
  end_time(COMMUNICATION_TIMER);
#endif /* RT_SINGLE_SOURCE */
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

#pragma omp parallel for default(none), private(i), shared(level,num_level_cells,level_cells,cell_vars)
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


/*
//  Computes the absorption coefficient at a single frequency bin 
//  freq on all cells from a supplied array and sends the result 
//  into abc[0] and (optionally) abc[1].
*/
void rtComputeAbsLevel(int level, int num_level_cells, int *level_cells, int freq, float **abc)
{
  int i, cell;
  frt_real dx, buffer[5], out[2];
#ifdef RT_ABSORPTION_CALLBACK_FULL
  frt_real var[FRT_DIM];
#else
  frt_real var[1];
#endif
  frt_intg idx;

  /* turn ifreq into a fortran index */
  idx = freq + 1;

#pragma omp parallel for default(none), private(cell,i,var,buffer,out,dx), shared(num_level_cells,level_cells,idx,abc,cell_vars,constants,level)
  for(i=0; i<num_level_cells; i++)
    {
      cell = level_cells[i];

      /*
      //  Not sure what to do here
      */
      //dx = 2*cell_size[level];
      dx = cell_sobolev_length2(cell,level,NULL);

#ifdef RT_ABSORPTION_CALLBACK_FULL
      rtPackCellData(level,cell,var,NULL);
#else  /* RT_ABSORPTION_CALLBACK_FULL */
#ifdef ENRICH
      var[0] = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
      var[0] = 0.0;
#endif
#endif /* RT_ABSORPTION_CALLBACK_FULL */

      if(sizeof(frt_real) == sizeof(float))  /* Optimization */
	{
#if (RT_CFI == 1)
	  frtCall(transfercomputecellabs)(&idx,&(cell_gas_density(cell)),&(cell_HI_density(cell)),&(cell_HeI_density(cell)),&(cell_HeII_density(cell)),&(cell_H2_density(cell)),&dx,out,var);
	  abc[0][i] = out[0];
	  abc[1][i] = out[1];
#else
	  frtCall(transfercomputecellabs)(&idx,&(cell_gas_density(cell)),&(cell_HI_density(cell)),&(cell_HeI_density(cell)),&(cell_HeII_density(cell)),&(cell_H2_density(cell)),&dx,abc[0]+i,var);
#endif
	}
      else
	{
	  buffer[0] = cell_gas_density(cell);
	  buffer[1] = cell_HI_density(cell);
	  buffer[2] = cell_HeI_density(cell);
	  buffer[3] = cell_HeII_density(cell);
	  buffer[4] = cell_H2_density(cell);
#if (RT_CFI == 1)
	  frtCall(transfercomputecellabs)(&idx,buffer+0,buffer+1,buffer+2,buffer+3,buffer+4,&dx,out,var);
	  abc[0][i] = out[0];
	  abc[1][i] = out[1];
#else
	  frtCall(transfercomputecellabs)(&idx,buffer+0,buffer+1,buffer+2,buffer+3,buffer+4,&dx,out,var);
	  abc[0][i] = out[0];
#endif
	}

#ifdef RT_DEBUG
      if(abc[0][i]<0.0 || isnan(abc[0][i]))
	{
	  cart_debug("Oops: %d %d %g",i,cell,abc[0][i]);
#ifdef RT_ABSORPTION_CALLBACK_FULL
	  cart_debug("Z:    %g",var[FRT_Metallicity]);
#else
	  cart_debug("Z:    %g",var[0]);
#endif
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

  dx0 *= rtSingleSourceValue*cell_volume_inverse[level];
  dx1 *= rtSingleSourceValue*cell_volume_inverse[level];

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

