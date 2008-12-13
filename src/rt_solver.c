#include "defs.h"      
#ifdef RADIATIVE_TRANSFER

#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */


#include "auxiliary.h"
#include "constants.h"
#include "hydro.h"
#include "timestep.h"
#include "timing.h"
#include "units.h"

#include "rt_c2f.h"
#include "rt_utilities.h"
#ifdef RT_TRANSFER
#include "rt_transfer.h"
#endif

#include "F/frt_parameters.ch"


float rt_XH;
float rt_XHe;

float rt_TemScale;


void rtCoolOff(int level, int cell, float *Ein, int thread);
void rtGetSobolevFactors(int cell, int level, float *len, float *vel);


/* 
//  Fortran interface 
*/
void f2c_wrapper(frtcooloff)(f2c_real *rPar, f2c_real *rRF0, f2c_real *rRF1, f2c_real *rTime, f2c_real *rVar, f2c_intg *iBuf);
void f2c_wrapper(frtinitrun)(f2c_real *Yp);
void f2c_wrapper(frtstepbeginart)(f2c_real *aexp, f2c_real *daexp, f2c_real *dt, f2c_real *Om0, f2c_real *h100, f2c_real *r0, f2c_real *t0);
void f2c_wrapper(frtstepend)(f2c_intg *id, f2c_real *vol, f2c_real *avg);

f2c_real f2c_wrapper(frtein)(f2c_real *tem, f2c_real *y);
f2c_real f2c_wrapper(frttem)(f2c_real *Ein, f2c_real *y);
f2c_real f2c_wrapper(frtgamma)(f2c_real *Ein, f2c_real *y);
f2c_real f2c_wrapper(frtqueryxh)(void);
f2c_real f2c_wrapper(frtqueryxhe)(void);
f2c_real f2c_wrapper(frtlogwithunits)(f2c_real *val, f2c_intg *pDen, f2c_intg *pLen, f2c_intg *pTime);

#ifdef RT_TEST
void f2c_wrapper(frttestmodifytimestep)(f2c_real *dt, f2c_intg *ns);
#endif

/*
d//  This size must be maintained to be always consistent with the parameter ioDim from file F/frt_base.inc
*/


/*
//  Applies cooling to all cells of a given level
*/
void rtApplyCooling(int level, int num_level_cells, int *level_cells)
{
  int i, j, icell, nchunk, thread;
  float Efact, Escale = rt_TemScale/(aexp[level]*aexp[level]);
  float vol = cell_volume[level];
  float tmp, Ein;

  /* 
  //  OpenMP chunk size 
  */
  nchunk = 128/(1<<level);
  if(nchunk < 1) nchunk = 1;

  /*
  //  Main loop
  */
#pragma omp parallel for default(none), private(i,icell,Ein,Efact,tmp,j,thread), shared(cell_vars,num_level_cells,cell_child_oct,level_cells,Escale,level,vol), schedule(dynamic,nchunk)
  for(i=0; i<num_level_cells; i++) if(cell_is_leaf((icell = level_cells[i])) && cell_gas_density(icell) > 0.0)  /* neg. density means a blow-up, let the code die gracefully in hydro_magic, not here */
    {
      /*
      //  Ein is the internal energy per baryon, = U/n_b = c_u*k_B*T/mu
      */
      Efact = Escale/cell_gas_density(icell);
      Ein = max(T_min,Efact*cell_gas_internal_energy(icell));

#ifdef _OPENMP
      thread = omp_get_thread_num();
#else
      thread = 0;
#endif /* _OPENMP */

      /*
      //  Wrapper over an actual (Fortran) worker
      */
      rtCoolOff(level,icell,&Ein,thread);

      cell_gas_internal_energy(icell) = max(T_min,Ein)/Efact;
      cell_gas_energy(icell) = cell_gas_kinetic_energy(icell) + cell_gas_internal_energy(icell);
    }
}


void rtInitRun()
{
#ifdef RT_TEST
  f2c_real Yp = 1.0e-10;
#else
  f2c_real Yp = Y_p;
#endif

  rtuInitRun();
  
  f2c_wrapper(frtinitrun)(&Yp);
  rt_XH = f2c_wrapper(frtqueryxh)();
  rt_XHe = f2c_wrapper(frtqueryxhe)();

#ifdef RT_TRANSFER
  rtInitRunTransfer();
#endif
}


/* This function can be called more than once per top level step */ 
void rtUpdateTables()
{
  start_time(RT_TABLES_TIMER);

#ifdef RT_TRANSFER
  rtUpdateTablesTransfer();
#endif

  /* Fill in the tables */
  f2c_wrapper(frtfillradiationtables)();

  end_time(RT_TABLES_TIMER);
}


void rtStepBegin()
{
  f2c_real rAStep = aexp[0];
  f2c_real rHStep = aexp[0] - aexp_old[0];
  f2c_real rTStep = dtl[0];
  f2c_real rOm0 = Omega0;
  f2c_real rH100 = hubble;
  f2c_real rR0 = r0;
  f2c_real rT0 = t0;

  rt_TemScale = T0/wmu; /* Need to set it here, since units may change after rtInitRun */

#ifdef COSMOLOGY
  rHStep = log(b2a(tl[0]+dtl[0])/b2a(tl[0]))/dtl[0];
#else
  rHStep = 0.0;
#endif

  f2c_wrapper(frtstepbeginart)(&rAStep,&rHStep,&rTStep,&rOm0,&rH100,&rR0,&rT0);

#ifdef RT_TRANSFER
  rtStepBeginTransfer();
#endif

  /* By default update tables once per step */
  rtUpdateTables();
}


void rtStepEnd()
{
  int i, n;
  f2c_intg id;
  f2c_real rVol = 1.0/num_root_cells;
  f2c_real *glob, rAvg[3];
  double *buffer;

  MESH_RUN_DECLARE(level,cell);
  double sumSrc, sumRho, sumRhoHI, sumRhoH2;

  /*
  //    2. Averages over mesh variables and their derivatives
  */
  sumSrc = sumRho = sumRhoHI = sumRhoH2 = 0.0;

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  /*#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_vars,cell_volume), reduction(+:sumSrc,sumRho,sumRhoHI,sumRhoH2) */
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
#ifdef RT_VAR_SOURCE
  sumSrc += cell_vars[cell][RT_VAR_SOURCE]*cell_volume[level];
#endif
#ifdef RT_CHEMISTRY
  sumRho += cell_gas_density(cell)*cell_volume[level];
  sumRhoHI += cell_HI_density(cell)*cell_volume[level];
  sumRhoH2 += cell_H2_density(cell)*cell_volume[level];
#endif

  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;

  buffer = (double *)cart_alloc(4*sizeof(double));
  
  buffer[0] = sumSrc;
  buffer[1] = sumRho;
  buffer[2] = sumRhoHI;
  buffer[3] = sumRhoH2;
  rtuGlobalAverage(4,buffer);

  rAvg[0] = buffer[0]/num_root_cells;
#ifdef RT_CHEMISTRY
  rAvg[1] = buffer[2]/buffer[1];
  rAvg[2] = buffer[3]/buffer[1];
#else
  rAvg[1] = rAvg[2] = 0.0;
#endif

  cart_free(buffer);

  id = local_proc_id;
  f2c_wrapper(frtstepend)(&id,&rVol,rAvg);

#ifdef RT_TRANSFER
  rtStepEndTransfer();
#endif  /* RT_TRANSFER */
}


void rtAfterAssignDensity1(int level)
{
  int num_level_cells;
  int *level_cells;

  start_time(RT_AFTER_DENSITY_TIMER);

#ifdef RT_TRANSFER
  /* assumes buffer gas density is up to date */
  select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
  rtAfterAssignDensityTransfer(level,num_level_cells,level_cells);
  cart_free(level_cells);
#endif

  end_time(RT_AFTER_DENSITY_TIMER);
}


void rtAfterAssignDensity2(int level, int num_level_cells, int *level_cells)
{
  start_time(RT_AFTER_DENSITY_TIMER);

#ifdef RT_TRANSFER
  rtAfterAssignDensityTransfer(level,num_level_cells,level_cells);
#endif

  end_time(RT_AFTER_DENSITY_TIMER);
}


void rtLevelUpdate(int level, MPI_Comm local_comm)
{
  start_time(RT_LEVEL_UPDATE_TIMER);

#ifdef RT_TRANSFER
  rtLevelUpdateTransfer(level,local_comm);
#endif

  end_time(RT_LEVEL_UPDATE_TIMER);
}


/*
//  Follow cooling evolution of one resolution element - convenient wrapper over a Fortran routine
*/
void rtCoolOff(int level, int cell, float *Ein, int thread)
{
  int i;
#ifdef RT_CHEMISTRY
  float soblen, sobvel;
#endif

  /* 
  //  Specify types for the Fortran interface
  */
  f2c_real rTime, rVar[IVAR_DIM], rPar[IPAR_DIM];
#ifdef RT_TRANSFER
  f2c_real rBuffer[rt_num_frequencies];
  f2c_real rRF0[2], *rRF1 = rBuffer;
#else
  f2c_real *rRF0 = 0, *rRF1 = 0;
#endif
  f2c_intg iBuffer = 1 + thread;

  for(i=0; i<IPAR_DIM; i++) rPar[i] = 0.0;

  /*
  //  Set parameters:
  //
  //    Density in code units
  */
  rPar[IPAR_RHOB] = cell_gas_density(cell);
  /*
  //    Metallicity in units of solar
  */
#ifdef METALCOOLING            
  rPar[IPAR_ZSOL] = cell_gas_metallicity(cell)/(0.02*cell_gas_density(cell));
#else
  rPar[IPAR_ZSOL] = 0.0;
#endif
  /*
  //    Sobolev length 
  */
#ifdef RT_CHEMISTRY
  rtGetSobolevFactors(cell,level,&soblen,&sobvel);
  rPar[IPAR_SOBL] = soblen;
  rPar[IPAR_NUMF] = sobvel*dtl[level]/cell_size[level];
#endif

#if defined(RT_DEBUG) && defined(RT_DEBUG_ONE_CELL_POSX) && defined(RT_DEBUG_ONE_CELL_POSY) && defined(RT_DEBUG_ONE_CELL_POSZ) && defined(RT_DEBUG_ONE_CELL_INFO)
  double pos[3];
  cell_position_double(cell,pos);
  pos[0] = (pos[0]-RT_DEBUG_ONE_CELL_POSX)*cell_size_inverse[level];
  pos[1] = (pos[1]-RT_DEBUG_ONE_CELL_POSY)*cell_size_inverse[level];
  pos[2] = (pos[2]-RT_DEBUG_ONE_CELL_POSZ)*cell_size_inverse[level];
  if(
     -0.5<pos[0] && pos[0]<=0.5 &&
     -0.5<pos[1] && pos[1]<=0.5 &&
     -0.5<pos[2] && pos[2]<=0.5)
    {
      rPar[IPAR_DEB] = RT_DEBUG_ONE_CELL_INFO;
    }
  else
    {
      rPar[IPAR_DEB] = 0.0;
    }
#endif

  /*
  //  The following may be a waste of time if the types are consistent, but compiler should take care of that 
  */
  rTime = dtl[level];

  /* 
  //  Pack elemental abundances 
  */
  rVar[IVAR_Ein] = *Ein;
  rVar[IVAR_XHI] = cell_HI_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XHII] = cell_HII_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XHeI] = cell_HeI_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XHeII] = cell_HeII_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XHeIII] = cell_HeIII_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XH2] = cell_H2_density(cell)/rPar[IPAR_RHOB];
#ifdef RT_8SPECIES
  rVar[IVAR_XH2p] = rVar[IVAR_XHm] = 0.0;
#endif

  /* 
  //  Pack radiation field 
  */
#ifdef RT_TRANSFER

  if(sizeof(f2c_real) != sizeof(float))  /* Optimization */
    {
      rRF1 = rBuffer;
      for(i=0; i<rt_num_frequencies; i++)
	{
	  rRF1[i] = cell_var(cell,rt_freq_offset+i);
	}
    }
  else
    {
      rRF1 = cell_vars[cell] + rt_freq_offset;
    }

#ifdef RT_VAR_OT_FIELD  
  rRF0[0] = cell_var(cell,RT_VAR_OT_FIELD);
  rRF0[1] = rt_ot_field_Avg.Value;
#else
  rRF0[0] = rRF0[1] = 0.0;
#endif

  /*
  //    Cell size for flux-conserving correction
  */
#ifdef RT_TRANSFER_FLUX_CONSERVING
  rPar[IPAR_CELL] = cell_size[level];
#endif

#endif  // RT_TRANSFER

  /*
  //  Call the Fortran worker 
  */
  f2c_wrapper(frtcooloff)(rPar,rRF0,rRF1,&rTime,rVar,&iBuffer);

  /*
  //  Unpack radiation field (if needed) 
  */
#if defined(RT_TRANSFER) && defined(RT_VARIABLE_RFIELD)
  if(sizeof(f2c_real) != sizeof(float))  /* Optimization */
    {
      for(i=0; i<rt_num_frequencies; i++)
	{
	  cell_var(cell,rt_freq_offset+i) = rRF1[i];
	}
    }
#endif
  
  /*
  //  Unpack elemental abundances 
  */
  *Ein = rVar[IVAR_Ein];
  cell_HI_density(cell) = rVar[IVAR_XHI]*rPar[IPAR_RHOB];
  cell_HII_density(cell) = rVar[IVAR_XHII]*rPar[IPAR_RHOB];
  cell_HeI_density(cell) = rVar[IVAR_XHeI]*rPar[IPAR_RHOB];
  cell_HeII_density(cell) = rVar[IVAR_XHeII]*rPar[IPAR_RHOB];
  cell_HeIII_density(cell) = rVar[IVAR_XHeIII]*rPar[IPAR_RHOB];
  cell_H2_density(cell) = rVar[IVAR_XH2]*rPar[IPAR_RHOB];

#ifndef RT_MONOATOMIC
#ifdef RT_HIGH_DENSITY
  cell_gas_gamma(cell) = f2c_wrapper(frtgamma)(Ein,rVar,rPar+1);
#else
  cell_gas_gamma(cell) = f2c_wrapper(frtgamma)(Ein,rVar);
#endif
#endif
}


/*
//  Sobolev approximation factors
*/
void rtGetSobolevFactors(int cell, int level, float *len, float *vel)
{
  int i, j, nb[num_neighbors];
  float s, d;
  float vc[nDim];
  
  cell_all_neighbors(cell,nb);

  /*
  //  Length factor
  */
  s = 0.0;
  for(i=0; i<nDim; i++)
    {
      d = cell_gas_density(nb[2*i+1]) - cell_gas_density(nb[2*i]);
      s += d*d;
    }

  /* 
  //  Factor of 1.0 is empirical, from detailed comparison of Sobolev
  //  approximation with a ray-tracer
  */
  *len = cell_size[level]*cell_gas_density(cell)/(1.0e-30+sqrt(s));
  if(*len > num_grid) *len = num_grid;

  /*
  //  Velocity factor
  */
  for(i=0; i<nDim; i++)
    {
      vc[i] = cell_var(cell,HVAR_MOMENTUM+i)/cell_gas_density(cell);
    }

  s = 0.0;
  for(j=0; j<num_neighbors; j++)
    {
      for(i=0; i<nDim; i++)
	{
	  d = cell_var(nb[j],HVAR_MOMENTUM+i)/cell_gas_density(nb[j]) - vc[i];
	  s += d*d;
	}
    }
  *vel = sqrt(s/num_neighbors);
}


/*
//  Used in testing
*/
void rtModifyTimeStep(double *dt)
{
#ifdef RT_TEST
  f2c_real rdt;
  f2c_intg ins;

  rdt = *dt;
  ins = step;

  f2c_wrapper(frttestmodifytimestep)(&rdt,&ins);

  *dt = rdt;
#endif
}


/*
//  Helper function to get gas temperature 7 other properties
*/
float rtTem(int cell)
{
  f2c_real buffer[7];

  if(sizeof(f2c_real) == sizeof(float))  /* Optimization */
    {
      return f2c_wrapper(frttem)(&(cell_gas_internal_energy(cell)),cell_vars[cell]+RT_HVAR_OFFSET-1);
    }
  else
    {
      buffer[0] = cell_gas_internal_energy(cell);
      buffer[1] = cell_HI_density(cell);
      buffer[2] = cell_HII_density(cell);
      buffer[3] = cell_HeI_density(cell);
      buffer[4] = cell_HeII_density(cell);
      buffer[5] = cell_HeIII_density(cell);
      buffer[6] = cell_H2_density(cell);
      return f2c_wrapper(frttem)(buffer+0,buffer);
    }
}


float rtTemInK(int cell)
{
  int level = cell_level(cell);
  return rt_TemScale*rtTem(cell)/(aexp[level]*aexp[level]);
}


#ifdef RT_TRANSFER
void f2c_wrapper(frttransferpackradiationfield)(f2c_real *par, f2c_real *y0, f2c_real *rawRF0, f2c_real *rawRF1, f2c_real *rf);
#endif
void f2c_wrapper(frtgetphotorates)(f2c_real *par, f2c_real *rf, f2c_intg *itab, f2c_real *y0, f2c_real *pRate);

void rtGetPhotoRates(int cell, float ionRates[3], float heatRates[3])
{
  int i;
  int level = cell_level(cell);
  float Escale = rt_TemScale/(aexp[level]*aexp[level]);

  /* 
  //  Specify types for the Fortran interface
  */
  f2c_real rVar[IVAR_DIM], rPar[IPAR_DIM], pRate[999];
#ifdef RT_TRANSFER
  f2c_real rBuffer[rt_num_frequencies];
  f2c_real rRF0[2], *rRF1 = rBuffer;
  f2c_real rf[1+2*rt_num_frequencies];
#else
  f2c_real *rRF0 = 0, *rRF1 = 0;
  f2c_real rf[1];
#endif
  f2c_intg iTab[2];

  for(i=0; i<IPAR_DIM; i++) rPar[i] = 0.0;

  /*
  //  Set parameters:
  //
  //    Density in code units
  */
  rPar[IPAR_RHOB] = cell_gas_density(cell);
  /*
  //    Metallicity in units of solar
  */
#ifdef METALCOOLING            
  rPar[IPAR_ZSOL] = cell_gas_metallicity(cell)/(0.02*cell_gas_density(cell));
#else
  rPar[IPAR_ZSOL] = 0.0;
#endif

  /* 
  //  Pack elemental abundances 
  */
  rVar[IVAR_Ein] = max(T_min,Escale*cell_gas_internal_energy(cell)/cell_gas_density(cell));
  rVar[IVAR_XHI] = cell_HI_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XHII] = cell_HII_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XHeI] = cell_HeI_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XHeII] = cell_HeII_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XHeIII] = cell_HeIII_density(cell)/rPar[IPAR_RHOB];
  rVar[IVAR_XH2] = cell_H2_density(cell)/rPar[IPAR_RHOB];
#ifdef RT_8SPECIES
  rVar[IVAR_XH2p] = rVar[IVAR_XHm] = 0.0;
#endif

  /* 
  //  Pack radiation field 
  */
#ifdef RT_TRANSFER

  if(sizeof(f2c_real) != sizeof(float))  /* Optimization */
    {
      rRF1 = rBuffer;
      for(i=0; i<rt_num_frequencies; i++)
	{
	  rRF1[i] = cell_var(cell,rt_freq_offset+i);
	}
    }
  else
    {
      rRF1 = cell_vars[cell] + rt_freq_offset;
    }

#ifdef RT_VAR_OT_FIELD  
  rRF0[0] = cell_var(cell,RT_VAR_OT_FIELD);
  rRF0[1] = rt_ot_field_Avg.Value;
#else
  rRF0[0] = rRF0[1] = 0.0;
#endif

#endif  // RT_TRANSFER

  /*
  //  Call the Fortran worker 
  */
#ifdef RT_TRANSFER
  f2c_wrapper(frttransferpackradiationfield)(rPar,rVar,rRF0,rRF1,rf);
#endif

  f2c_wrapper(frtgetphotorates)(rPar,rf,iTab,rVar,pRate);

  ionRates[0] = pRate[5];
  ionRates[1] = pRate[3];
  ionRates[2] = pRate[1];

  heatRates[0] = pRate[4];
  heatRates[1] = pRate[2];
  heatRates[2] = pRate[0];
}


#ifdef RT_DEBUG
void rtPrintValue(const char *name, float val, int pDen, int pLen, int pTime)
{
  f2c_real fval = val;
  f2c_intg fpDen = pDen;
  f2c_intg fpLen = pLen;
  f2c_intg fpTime = pTime;

  val = f2c_wrapper(frtlogwithunits)(&fval,&fpDen,&fpLen,&fpTime);
  cart_debug("Checking: lg(%s) = %g",name,val);
}
#endif

#endif  /* RADIATIVE_TRANSFER */
