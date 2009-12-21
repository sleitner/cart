#include "defs.h"      

#include <math.h>

#include "hydro.h"
#include "timestep.h"
#include "units.h"


float rt_TemScale;

float rtTem(int cell);


void rtSetTemUnits()
{
#ifdef RADIATIVE_TRANSFER
  rt_TemScale = T0/wmu;
#else
  rt_TemScale = T0;
#endif
}


float rtTemInK(int cell)
{
  int level = cell_level(cell);
  float uTem = rt_TemScale/(abox[level]*abox[level]);
#ifdef RADIATIVE_TRANSFER
  return uTem*rtTem(cell);
#else
  return uTem*cell_gas_internal_energy(cell)/cell_gas_density(cell);
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
//  RT-dependent code
*/
#ifdef RADIATIVE_TRANSFER

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */


#include "auxiliary.h"
#include "timing.h"

#include "rt_c2f.h"
#include "rt_utilities.h"
#ifdef RT_TRANSFER
#include "rt_transfer.h"
#endif

#include "F/frt_parameters.ch"

#ifdef RT_DEBUG
#include "rt_debug.h"
#else
int rt_debug = 0;  /* for OpenMP pragmas */
#endif


void rtPackCellData(int level, int cell, f2c_real rVar[], f2c_real rPar[], f2c_real *rRadField0, f2c_real **pRadField1);
void rtUnPackCellData(int level, int cell, f2c_real rVar[], f2c_real rPar[], f2c_real *rRadField1);


/* 
//  Fortran interface 
*/
void f2c_wrapper(frtcooloff)(f2c_real *rPar, f2c_real *rRadField0, f2c_real *rRadField1, f2c_real *rTime, f2c_real *rVar, f2c_intg *info);
void f2c_wrapper(frtinitrun)(f2c_real *Yp);
void f2c_wrapper(frtstepbeginart)(f2c_real *astep, f2c_real *dstep, f2c_real *dt, f2c_real *Om0, f2c_real *h100, f2c_real *r0, f2c_real *t0);
void f2c_wrapper(frtstepend)(f2c_real *vol, f2c_real *avg);
void f2c_wrapper(frtupdatetables)();

f2c_real f2c_wrapper(frtein)(f2c_real *tem, f2c_real *y);
f2c_real f2c_wrapper(frttem)(f2c_real *Ein, f2c_real *y);
f2c_real f2c_wrapper(frtgamma)(f2c_real *Ein, f2c_real *y);
f2c_real f2c_wrapper(frtlogwithunits)(f2c_real *val, f2c_intg *pDen, f2c_intg *pLen, f2c_intg *pTime);

#ifdef RT_TEST
void f2c_wrapper(frttestmodifytimestep)(f2c_real *dt, f2c_intg *ns);
#endif


/*
//  Applies cooling to all cells of a given level
*/
void rtApplyCooling(int level, int num_level_cells, int *level_cells)
{
  int i, cell, nchunk;
  float soblen, sobvel;

  /* 
  //  Specify types for the Fortran interface
  */
  f2c_real rTime, rVar[IVAR_DIM], rPar[IPAR_DIM];
#ifdef RT_TRANSFER
  f2c_real rBuffer[rt_num_frequencies];
  f2c_real rRadField0[2];
#else
  f2c_real *rBuffer = 0, *rRadField0 = 0;
#endif
  f2c_real *rRadField1 = rBuffer;
  f2c_intg info;

  /*
  //  The following may be a waste of time if the types are consistent, but compiler should take care of that 
  */
  rTime = dtl[level];

  /* 
  //  OpenMP chunk size 
  */
  nchunk = 128/(1<<level);
  if(nchunk < 1) nchunk = 1;

  /*
  //  Main loop
  */
#pragma omp parallel for default(none), private(i,cell,rVar,rPar,rBuffer,rRadField0,rRadField1,soblen,sobvel,info), shared(cell_vars,num_level_cells,cell_child_oct,level_cells,level,rTime,cell_size,rt_debug,nchunk), schedule(dynamic,nchunk)
  for(i=0; i<num_level_cells; i++) if(cell_is_leaf((cell = level_cells[i])) && cell_gas_density(cell) > 0.0)  /* neg. density means a blow-up, let the code die gracefully in hydro_magic, not here */
    {
      
      rtPackCellData(level,cell,rVar,rPar,rRadField0,&rRadField1);

      /*
      //    Cell size for flux-conserving correction
      */
#if defined(RT_TRANSFER) && defined(RT_TRANSFER_FLUX_CONSERVING)
      rPar[IPAR_CELL] = cell_size[level];
#endif  // defined(RT_TRANSFER) && defined(RT_TRANSFER_FLUX_CONSERVING)

      /*
      //    Sobolev length 
      */
#ifdef RT_CHEMISTRY
      rtGetSobolevFactors(cell,level,&soblen,&sobvel);
      rPar[IPAR_SOBL] = soblen;
      rPar[IPAR_NUMF] = sobvel*rTime/cell_size[level];
#endif

#ifdef RT_DEBUG
      if(rt_debug.Mode==1 && cell==cell_find_position(rt_debug.Pos))
	{
	  rPar[IPAR_DEB] = 1 + (rt_debug.Stop>1 ? 0.5 : 0.0);
	  cart_debug("In cell-level debug for cell %d#%d",cell,cell_level(cell));
	}
      else
	{
	  rPar[IPAR_DEB] = 0.0;
	}
#endif

      /*
      //  Call the Fortran worker 
      */
      f2c_wrapper(frtcooloff)(rPar,rRadField0,rRadField1,&rTime,rVar,&info);

      rtUnPackCellData(level,cell,rVar,rPar,rRadField1);
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

#ifdef RT_TRANSFER
  rtInitRunTransfer();
#endif
}


/* This function can be called more than once per top level step */ 
void rtUpdateTables()
{
  start_time(RT_TABLES_TIMER);

  /* Fill in the tables */
  f2c_wrapper(frtupdatetables)();

  end_time(RT_TABLES_TIMER);
}


void rtStepBegin()
{
  f2c_real rAStep = abox[0];
  f2c_real rHStep = abox[0] - abox_old[0];
  f2c_real rTStep = dtl[0];
  f2c_real rOmegaM = cosmology->OmegaM;
  f2c_real rHubble = cosmology->h;
  f2c_real rR0 = r0;
  f2c_real rT0 = t0;

  rtSetTemUnits();  /* Need to set it here, since units may change after rtInitRun */

#ifdef COSMOLOGY
  rHStep = log(abox_from_tcode(tl[0]+dtl[0])/abox_from_tcode(tl[0]))/dtl[0];
#else
  rHStep = 0.0;
#endif

  f2c_wrapper(frtstepbeginart)(&rAStep,&rHStep,&rTStep,&rOmegaM,&rHubble,&rR0,&rT0);

#ifdef RT_TRANSFER
  rtStepBeginTransfer();
#endif

#ifdef RT_DEBUG
  switch(rt_debug.Mode)
    {
    case 1:
      {
	int i, cell;
	cell = cell_find_position(rt_debug.Pos);
	cart_debug("In cell-level debug for cell %d/%d",cell,cell_level(cell));
	cart_debug("RT_HVAR_OFFSET: %d",RT_HVAR_OFFSET);
	cart_debug("RT_VAR_SOURCE: %d",RT_VAR_SOURCE);
	cart_debug("rt_grav_vars_offset: %d",rt_grav_vars_offset);
#ifdef RT_TRANSFER
	cart_debug("rt_num_vars: %d",rt_num_vars);
#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
	cart_debug("RT_VAR_OT_FIELD: %d",RT_VAR_OT_FIELD);
	cart_debug("rt_et_offset: %d",rt_et_offset);
	cart_debug("rt_freq_offset: %d",rt_freq_offset);
#endif
#endif
	for(i=0; i<num_vars; i++)
	  {
	    cart_debug("Var[%d] = %g",i,cell_var(cell,i));
	  }
	break;
      }
    }
#endif  /* RT_DEBUG */
}


void rtStepEnd()
{
  int i;
  f2c_real rVol = num_root_cells;
  f2c_real rAvg[3];
  double *buffer;

  MESH_RUN_DECLARE(level,cell);
  double sumSrc, sumRho, sumRhoHI, sumRhoH2;

  /*
  //    Averages over mesh variables and their derivatives
  */
  sumSrc = sumRho = sumRhoHI = sumRhoH2 = 0.0;

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_vars,cell_volume,cell_child_oct), reduction(+:sumSrc,sumRho,sumRhoHI,sumRhoH2)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
#ifdef RT_VAR_SOURCE
      sumSrc += cell_vars[cell][RT_VAR_SOURCE]*cell_volume[level];
#endif
#ifdef RT_CHEMISTRY
      sumRho += cell_gas_density(cell)*cell_volume[level];
      sumRhoHI += cell_HI_density(cell)*cell_volume[level];
      sumRhoH2 += cell_H2_density(cell)*cell_volume[level];
#endif
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;

  buffer = cart_alloc(double, 4 );
  
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

  f2c_wrapper(frtstepend)(&rVol,rAvg);

#ifdef RT_TRANSFER
  rtStepEndTransfer();
#endif  /* RT_TRANSFER */

#ifdef RT_DEBUG
  switch(rt_debug.Mode)
    {
    case 1:
      {
	cell = cell_find_position(rt_debug.Pos);
	cart_debug("In cell-level debug for cell %d/%d",cell,cell_level(cell));
	for(i=0; i<num_vars; i++)
	  {
	    cart_debug("Var[%d] = %g",i,cell_var(cell,i));
	  }
	break;
      }
    }
  if(rt_debug.Mode>0 && rt_debug.Stop)
    {
      cart_debug("Aborting on request...");
      MPI_Abort(MPI_COMM_WORLD,0);
      exit(0);
    }
#endif
}


void rtAfterAssignDensity1(int level)
{
  int num_level_cells;
  int *level_cells;

  start_time(RT_AFTER_DENSITY_TIMER);

#ifdef RT_TRANSFER
  /* assumes buffer gas density is up to date */
  select_level(level,CELL_TYPE_ANY,&num_level_cells,&level_cells);
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
//  Helper functions for Fortran workers
*/
void rtPackCellData(int level, int cell, f2c_real rVar[], f2c_real rPar[], f2c_real *rRadField0, f2c_real **pRadField1)
{
  int i;
  float uTem = rt_TemScale/(abox[level]*abox[level]);

  /*
  //  Default values of all parameters is 0
  */
  for(i=0; i<IPAR_DIM; i++) rPar[i] = 0.0;

  /*
  //  Set parameters:
  //
  //    Density in code units
  */
  rPar[IPAR_RHO] = cell_gas_density(cell);
  /*
  //    Cell volume in code units
  */
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
  rPar[IPAR_VOL] = cell_volume[level];
#endif
  /*
  //    Metallicity in units of solar
  */
#ifdef METALCOOLING            
  rPar[IPAR_ZSOL] = cell_gas_metallicity(cell)/(Zsolar*cell_gas_density(cell));
#else
  rPar[IPAR_ZSOL] = 0.0;
#endif

  /* 
  //  Pack elemental abundances 
  */
  rVar[IVAR_Ein] = uTem*cell_gas_internal_energy(cell)/rPar[IPAR_RHO];
  rVar[IVAR_XHI] = cell_HI_density(cell)/rPar[IPAR_RHO];
  rVar[IVAR_XHII] = cell_HII_density(cell)/rPar[IPAR_RHO];
  rVar[IVAR_XHeI] = cell_HeI_density(cell)/rPar[IPAR_RHO];
  rVar[IVAR_XHeII] = cell_HeII_density(cell)/rPar[IPAR_RHO];
  rVar[IVAR_XHeIII] = cell_HeIII_density(cell)/rPar[IPAR_RHO];
  rVar[IVAR_XH2] = cell_H2_density(cell)/rPar[IPAR_RHO];
#ifdef RT_8SPECIES
  rVar[IVAR_XH2p] = rVar[IVAR_XHm] = 0.0;
#endif

  /* 
  //  Pack radiation field 
  */
#ifdef RT_TRANSFER
  if(sizeof(f2c_real) != sizeof(float))  /* Optimization */
    {
      for(i=0; i<rt_num_frequencies; i++)
	{
	  (*pRadField1)[i] = cell_var(cell,rt_freq_offset+i);
	}
    }
  else
    {
      *pRadField1 = cell_vars[cell] + rt_freq_offset;
    }
#ifdef RT_VAR_OT_FIELD  
  rRadField0[0] = cell_var(cell,RT_VAR_OT_FIELD);
  rRadField0[1] = rt_ot_field_Avg.Value;
#else
  rRadField0[0] = rRadField0[1] = 0.0;
#endif
#endif  // RT_TRANSFER
}


void rtUnPackCellData(int level, int cell, f2c_real rVar[], f2c_real rPar[], f2c_real *rRadField1)
{
  float uTem = rt_TemScale/(abox[level]*abox[level]);

  /*
  //  Unpack radiation field (if needed) 
  */
#if defined(RT_TRANSFER) && defined(RT_VARIABLE_RFIELD)
  int i;
  if(sizeof(f2c_real) != sizeof(float))  /* Optimization */
    {
      for(i=0; i<rt_num_frequencies; i++)
	{
	  cell_var(cell,rt_freq_offset+i) = rRadField1[i];
	}
    }
#endif
  
  /*
  //  Unpack elemental abundances 
  */
  cell_gas_internal_energy(cell) = max(gas_temperature_floor,rVar[IVAR_Ein])*rPar[IPAR_RHO]/uTem;
  cell_gas_energy(cell) = cell_gas_kinetic_energy(cell) + cell_gas_internal_energy(cell);
  cell_HI_density(cell) = rVar[IVAR_XHI]*rPar[IPAR_RHO];
  cell_HII_density(cell) = rVar[IVAR_XHII]*rPar[IPAR_RHO];
  cell_HeI_density(cell) = rVar[IVAR_XHeI]*rPar[IPAR_RHO];
  cell_HeII_density(cell) = rVar[IVAR_XHeII]*rPar[IPAR_RHO];
  cell_HeIII_density(cell) = rVar[IVAR_XHeIII]*rPar[IPAR_RHO];
  cell_H2_density(cell) = rVar[IVAR_XH2]*rPar[IPAR_RHO];

#ifndef RT_MONOATOMIC
#ifdef RT_HIGH_DENSITY
  cell_gas_gamma(cell) = f2c_wrapper(frtgamma)(Ein,rVar,rPar+1);
#else
  cell_gas_gamma(cell) = f2c_wrapper(frtgamma)(Ein,rVar);
#endif
#endif
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
//  Helper function to get gas temperature
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


#ifdef RT_TRANSFER
void f2c_wrapper(frttransferpackradiationfield)(f2c_real *par, f2c_real *y0, f2c_real *rawRF0, f2c_real *rawRF1, f2c_real *rf);
#endif
void f2c_wrapper(frtgetphotorates)(f2c_real *par, f2c_real *rf, f2c_intg *itab, f2c_real *y0, f2c_real *pRate);
void f2c_wrapper(frtquerybackground)(f2c_intg *n, f2c_real *wlen, f2c_real *nxi);
f2c_real f2c_wrapper(frtgetradiationfield)(f2c_real *rf, f2c_intg *lr);


void rtGetPhotoRates(int cell, float rate[])
{
  int i;
  int level = cell_level(cell);

  /* 
  //  Specify types for the Fortran interface
  */
  f2c_real rVar[IVAR_DIM], rPar[IPAR_DIM], pRate[IRATE_DIM];
#ifdef RT_TRANSFER
  f2c_real rBuffer[rt_num_frequencies];
  f2c_real rRadField0[2];
#else
  f2c_real *rBuffer = 0, *rRadField0 = 0;
#endif
  f2c_real rf[1+2*rt_num_frequencies];
  f2c_real *rRadField1 = rBuffer;
  f2c_intg iTab[2];

  rtPackCellData(level,cell,rVar,rPar,rRadField0,&rRadField1);

  /*
  //  Call Fortran workers
  */
#ifdef RT_TRANSFER
  f2c_wrapper(frttransferpackradiationfield)(rPar,rVar,rRadField0,rRadField1,rf);
#endif

  if(sizeof(f2c_real) == sizeof(float))
    {
      f2c_wrapper(frtgetphotorates)(rPar,rf,iTab,rVar,(f2c_real *)rate);
    }
  else
    {
      f2c_wrapper(frtgetphotorates)(rPar,rf,iTab,rVar,pRate);
      for(i=0; i<IRATE_DIM; i++) rate[i] = pRate[i];
    }
}


void rtGetRadiationBackground(int *nPtr, float **wlenPtr, float **ngxiPtr)
{
  int i, n;
  float *wlen, *ngxi;
  f2c_intg iNum;
  f2c_real *pWlen, *pNgxi;

  iNum = 0;
  f2c_wrapper(frtquerybackground)(&iNum,NULL,NULL);
  *nPtr = n = iNum;

  cart_assert(n > 0);

  *wlenPtr = wlen = cart_alloc(float,n);
  *ngxiPtr = ngxi = cart_alloc(float,n);

  if(sizeof(f2c_real) == sizeof(float))
    {
      f2c_wrapper(frtquerybackground)(&iNum,(f2c_real *)wlen,(f2c_real *)ngxi);
    }
  else
    {
      pWlen = cart_alloc(f2c_real,n);
      pNgxi = cart_alloc(f2c_real,n);

      f2c_wrapper(frtquerybackground)(&iNum,pWlen,pNgxi);

      for(i=0; i<n; i++)
	{
	  wlen[i] = pWlen[i];
	  ngxi[i] = pNgxi[i]; 
	}

      cart_free(pWlen);
      cart_free(pNgxi);
    }
}


void rtGetRadiationField(int cell, int n, int lxi[], float ngxi[])
{
  int i;
  int level = cell_level(cell);

  /* 
  //  Specify types for the Fortran interface
  */
  f2c_real rVar[IVAR_DIM], rPar[IPAR_DIM], pRate[IRATE_DIM];
#ifdef RT_TRANSFER
  f2c_real rBuffer[rt_num_frequencies];
  f2c_real rRadField0[2];
#else
  f2c_real *rBuffer = 0, *rRadField0 = 0;
#endif
  f2c_real rf[1+2*rt_num_frequencies];
  f2c_real *rRadField1 = rBuffer;
  f2c_intg lr;

  rtPackCellData(level,cell,rVar,rPar,rRadField0,&rRadField1);

  /*
  //  Call Fortran workers
  */
#ifdef RT_TRANSFER
  f2c_wrapper(frttransferpackradiationfield)(rPar,rVar,rRadField0,rRadField1,rf);
#endif

  for(i=0; i<n; i++)
    {
      lr = lxi[i] + 1;  /* C-Fortran conversion!!! */
      ngxi[i] = f2c_wrapper(frtgetradiationfield)(rf,&lr);
    }
}

#endif  /* RADIATIVE_TRANSFER */
