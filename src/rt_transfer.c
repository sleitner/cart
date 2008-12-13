#include "defs.h"
#ifdef RADIATIVE_TRANSFER

#include <math.h>

#include "auxiliary.h"
#include "particle.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "rt_c2f.h"
#include "rt_utilities.h"


float rt_trad_on;


float rtSource(int ipart)
{
  int istar;
  float t1on, t2on;

#if defined(PARTICLES) && defined(STARFORM)
  if(!particle_is_star(ipart)) return 0.0;
#endif

#ifdef RT_TEST
  return 1.0;
#endif

#if defined(PARTICLES) && defined(STARFORM)
  /*
  //  The convention is different from HART
  */
  istar = ipart;
  /*
  // ******************************************************************
  //
  //     EXTREMELY IMPORTANT!!!!!
  //
  //  Specific form of this function depends on the order in which things
  //  happen in ART_Step. Right now this assumes that Move_Level is done
  //  before Assign_Density so that when we call this function, pt(is)
  //  is already pdt(is) after the real moment when this star starts emitting
  //  energy. So, the right quantity is the intergral from pt-pdt to pt
  //  (so that the total number of emitted photons does not depend on the time
  //  step).
  //
  // ******************************************************************
  */
  t2on = (float)(particle_t[istar]-star_tbirth[istar]); if(t2on < 0.0) t2on = 0.0;
  t1on = (float)(t2on-particle_dt[ipart]); if(t1on < 0.0) t1on = 0.0;

  if(t1on < 100*rt_trad_on)
    {
      return (exp(-t1on/rt_trad_on)-exp(-t2on/rt_trad_on))/particle_dt[ipart];
    }
  else
    {
      return 0.0;
    }
#else
  return 1.0;
#endif /* defined(PARTICLES) && defined(STARFORM) */
}


int rtIsThereWork()
{
  int num_local, num_global;

#ifdef RT_TEST
  return 1;
#endif

#if defined(PARTICLES) && defined(STARFORM)

  if(num_particle_species > 1)
    {
      num_local = num_local_star_particles;
    }
  else
    {
      /*
      //  This assumes that if we have just one particle species and
      //  RADIATIVE_TRANSFER is on, then this species is sources
      */
      num_local = particle_species_num[0];
    }

  MPI_Allreduce(&num_local,&num_global,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  return (num_global > 0);

#else

  return 0;

#endif /* defined(PARTICLES) && defined(STARFORM) */
}


/* ******************************************************* */

#ifdef RT_TRANSFER

#include "rt_transfer.h"


struct rtArrayAverageData rt_glob_Avg[2];

void rtTransferSplitUpdate(int level);
void rtTransferUpdateFields(int level1, int level2, int nvar, int var0, struct rtArrayAverageData *out); 
void rtTransferUpdateFieldAverages(int nvar, struct rtArrayAverageData *out); 

void f2c_wrapper(frtinitruntransfer)(f2c_intg *nfreq);
void f2c_wrapper(frtstepbegintransfer)();
void f2c_wrapper(frtupdatetablestransfer)(f2c_real *rfAvg);
#if (RT_CFI == 1)
void f2c_wrapper(frttransfercomputecellabs)(f2c_intg *L, f2c_real *Zsol, f2c_real *denB, f2c_real *denH1, f2c_real *denG1, f2c_real *denG2, f2c_real *denMH, f2c_real *abc, f2c_real *abc1);
#else
void f2c_wrapper(frttransfercomputecellabs)(f2c_intg *L, f2c_real *Zsol, f2c_real *denB, f2c_real *denH1, f2c_real *denG1, f2c_real *denG2, f2c_real *denMH, f2c_real *abc);
#endif


#ifndef RT_VAR_SOURCE
void rtTransferUpdateUniformSource(rtArrayAverageData *avg);
#endif


#ifdef RT_SINGLE_SOURCE
float rtSingleSourceVal;
double rtSingleSourcePos[nDim];

#ifdef RT_VAR_SOURCE
void rtTransferAssignSingleSourceDensity(int level);
#endif
#endif


#if RT_TRANSFER_METHOD == RT_METHOD_OTVET
#include "rt_otvet.h"
#endif


void rtInitRunTransfer()
{
  int i, level;
  f2c_intg val = rt_num_frequencies;

  f2c_wrapper(frtinitruntransfer)(&val);

  for(i=0; i<rt_num_glob; i++)
    {
      rt_glob_Avg[i].Value = 0.0;
      for(level=min_level; level<=max_level; level++)
	{
	  rt_glob_Avg[i].LevelSum[level] = 0.0;
	}
    }

#ifdef RT_SINGLE_SOURCE
  rtSingleSourceVal = 0.0;
  for(i=0; i<nDim; i++) rtSingleSourcePos[i] = 0.5*num_grid;
#endif

#if RT_TRANSFER_METHOD == RT_METHOD_OTVET

  rtInitRunTransferOtvet();

#endif
}


void rtStepBeginTransfer()
{
  int level;

  /* Time the source is on (20 Myr) */
  rt_trad_on = 2.0e7/(t0*aexp[0]*aexp[0]);

  f2c_wrapper(frtstepbegintransfer)();

#if RT_TRANSFER_METHOD == RT_METHOD_OTVET

  rtStepBeginTransferOtvet();

#endif

  rtTransferUpdateFields(min_level,max_level,rt_num_glob,RT_VAR_SOURCE,rt_glob_Avg);
}


void rtStepEndTransfer()
{
}


void rtLevelUpdateTransfer(int level, MPI_Comm local_comm)
{
  if(!rtIsThereWork()) return;

#if RT_TRANSFER_METHOD == RT_METHOD_OTVET

  rtLevelUpdateTransferOtvet(level,local_comm);

#endif

  rtTransferSplitUpdate(level);
}


void rtAfterAssignDensityTransfer(int level, int num_level_cells, int *level_cells)
{
  int i;
  double sum = 0.0;


#ifdef RT_VAR_SOURCE

#ifdef PARTICLES

  /*
  // If we have a source field that was set inside density(...), 
  // turn the mass per cell into density.
  */
#pragma omp parallel for default(none), private(i), shared(level,num_level_cells,level_cells,cell_vars,cell_volume_inverse)
  for(i=0; i<num_level_cells; i++)
    {
      cell_rt_source(level_cells[i]) *= cell_volume_inverse[level];
    }

#else  // PARTICLES

#ifdef RT_SINGLE_SOURCE
  /*
  // Set the source field from a single source
  */
  rtTransferAssignSingleSourceDensity(level);
#else  // RT_SINGLE_SOURCE

#error "Invalid set of switches: either PARTICLES or RT_SINGLE_SOURCE must be defined."

#endif  // RT_SINGLE_SOURCE

#endif  // PARTICLES

  //#pragma omp parallel for default(none), private(i), shared(level,num_level_cells,level_cells,cell_vars), reduction(+:sum)
  for(i=0; i<num_level_cells; i++)
    {
      sum += cell_rt_source(level_cells[i]);
    }
  rt_glob_Avg[0].LevelSum[level] = sum;
 
  rtTransferUpdateFieldAverages(1,rt_glob_Avg);

#else

  rtTransferUpdateUniformSource(&rt_glob_Avg[0]);  NOT IMPLEMENTED

#endif // RT_VAR_SOURCE
}


void rtUpdateTablesTransfer()
{
  int i;
  struct rtArrayAverageData tmp[rt_num_frequencies];
  f2c_real rfAvg[rt_num_frequencies];
  
  /*
  //  Compute global properties for all frequencies 
  */
  rtTransferUpdateFields(min_level,max_level,rt_num_frequencies,rt_freq_offset,tmp); 
  for(i=0; i<rt_num_frequencies; i++)
    {
      rfAvg[i] = tmp[i].Value;
    }

  f2c_wrapper(frtupdatetablestransfer)(rfAvg);
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
}


void rtTransferUpdateFields(int level1, int level2, int nvar, int var0, struct rtArrayAverageData *out)
{
  int i;
  double buffer[rt_num_vars];
  MESH_RUN_DECLARE(level,cell);

  /*
  //  Compute per-level averages
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2);

  for(i=0; i<nvar; i++) buffer[i] = 0.0;

  /*#pragma omp parallel for default(none), private(_Index,cell,i), shared(_Num_level_cells,_Level_cells,level,nvar,var0,cell_vars), reduction(+:buffer) */
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  for(i=0; i<nvar; i++) buffer[i] += cell_var(cell,var0+i);
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

  for(i=0; i<nvar; i++) out[i].LevelSum[level] = buffer[i];
  
  MESH_RUN_OVER_LEVELS_END;

  rtTransferUpdateFieldAverages(nvar,out);
}


void rtTransferUpdateFieldAverages(int nvar, struct rtArrayAverageData *out)
{
  int i, level, num_levels = 1 + max_level_local();
  double buffer[rt_num_vars];

  /*
  //  Compute local averages
  */
  for(i=0; i<nvar; i++)
    {
      buffer[i] = 0.0;
      for(level=min_level; level<num_levels; level++) buffer[i] += out[i].LevelSum[level]*cell_volume[level];
    }

  /*
  //  Compute global averages
  */
  rtuGlobalAverage(nvar,buffer);  
  for(i=0; i<nvar; i++) out[i].Value = buffer[i]/num_root_cells;
}


/*
//  Computes the absorption coefficient at a single frequency bin 
//  ifreq on all cells from a supplied array and sends the result 
//  into abc[0] and (optionally) abc[1].
*/
void rtComputeAbsLevel(int ncells, int *cells, int ifreq, float **abc)
{
  int i, cell;
  float rho, Zsol;
  f2c_real buffer[5];

  /* turn ifreq into a fortran index */
  ifreq++;

#pragma omp parallel for default(none), private(cell,i,rho,Zsol,buffer), shared(ncells,cells,ifreq,abc,cell_vars)
  for(i=0; i<ncells; i++)
    {
      cell = cells[i];
      rho = cell_gas_density(cell);

#ifdef RT_DUST
#ifdef ENRICH
      Zsol = cell_gas_metallicity_II(cell);
#ifdef ENRICH_SNIa
      Zsol += cell_gas_metallicity_Ia(cell);
#endif
      Zsol /= (0.02*rho);
#else
      Zsol = 0.0;
#endif
#else
      Zsol = 0.0;
#endif  // RT_DUST

      if(sizeof(f2c_real) != sizeof(float))  /* Optimization */
	{
	  buffer[0] = cell_gas_density(cell);
	  buffer[1] = cell_HI_density(cell);
	  buffer[2] = cell_HeI_density(cell);
	  buffer[3] = cell_HeII_density(cell);
	  buffer[4] = cell_H2_density(cell);
#if (RT_CFI == 1)
	  f2c_wrapper(frttransfercomputecellabs)(&ifreq,&Zsol,buffer+0,buffer+1,buffer+2,buffer+3,buffer+4,abc[0]+i,abc[1]+i);
#else
	  f2c_wrapper(frttransfercomputecellabs)(&ifreq,&Zsol,buffer+0,buffer+1,buffer+2,buffer+3,buffer+4,abc[0]+i);
#endif
	}
      else
	{
#if (RT_CFI == 1)
	  f2c_wrapper(frttransfercomputecellabs)(&ifreq,&Zsol,&(cell_gas_density(cell)),&(cell_HI_density(cell)),&(cell_HeI_density(cell)),&(cell_HeII_density(cell)),&(cell_H2_density(cell)),abc[0]+i,abc[1]+i);
#else
	  f2c_wrapper(frttransfercomputecellabs)(&ifreq,&Zsol,&(cell_gas_density(cell)),&(cell_HI_density(cell)),&(cell_HeI_density(cell)),&(cell_HeII_density(cell)),&(cell_H2_density(cell)),abc[0]+i);
#endif
	}
    }
}


void f2c_wrapper(rttransfergetglobalabs)(f2c_intg *nfields, f2c_intg *ifield, f2c_real *abcAvg)
{
  int i, n = *nfields;
  double buffer[3*rt_num_frequencies];
  MESH_RUN_DECLARE(level,cell);
  float *abc[2], *abc1;

  for(i=0; i<3*n; i++) buffer[i] = 0.0;
  
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);

  abc[0] = (float *)cart_alloc(_Num_level_cells*sizeof(float));
#if (RT_CFI == 1)
  abc[1] = (float *)cart_alloc(_Num_level_cells*sizeof(float));
#else
  abc[1] = 0;
#endif

  for(i=0; i<n; i++)
    {
      rtComputeAbsLevel(_Num_level_cells,_Level_cells,i,abc);
#if (RT_CFI == 1)
      abc1 = abc[1];
#else
      abc1 = abc[0];
#endif

  /*#pragma omp parallel for default(none), private(_Index,cell,i), shared(_Num_level_cells,_Level_cells,level,n,cell_vars), reduction(+:sum) */
      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

      buffer[3*i+0] += abc1[_Index]*cell_volume[level];
      buffer[3*i+1] += cell_var(cell,ifield[i])*cell_volume[level];
      buffer[3*i+2] += cell_var(cell,ifield[i])*abc1[_Index]*cell_volume[level];
  
      MESH_RUN_OVER_CELLS_OF_LEVEL_END;

    }

  MESH_RUN_OVER_LEVELS_END;

  rtuGlobalAverage(3*n,buffer);

  for(i=0; i<n; i++)
    {
      if(buffer[3*i+1] > 1.0e-35)
	{
	  abcAvg[i] = buffer[3*i+2]/buffer[3*i+1];
	}
      else
	{
	  abcAvg[i] = buffer[3*i+0]/num_root_cells;
	}
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
}

#endif  // defined(RT_SINGLE_SOURCE) && defined(RT_VAR_SOURCE)

#endif  // RT_TRANSFER
#endif  // RADIATIVE_TRANSFER
