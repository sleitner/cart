#include "config.h"

#include <math.h>

#include "auxiliary.h"
#include "tree.h"


/*
//  Sobolev approximation factors
*/
#ifdef HYDRO
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
      d = 0.5*(cell_gas_density(nb[2*i+1])-cell_gas_density(nb[2*i]));
      s += d*d;
    }

  /* 
  //  Factor of 0.5 is empirical, from detailed comparison of Sobolev
  //  approximation with a ray-tracer
  */
  *len = 0.5*cell_size[level]*cell_gas_density(cell)/(1.0e-30+sqrt(s));
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
#endif /* HYDRO */


/*
//  RT-dependent code
*/
#ifdef RADIATIVE_TRANSFER

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "iterators.h"
#include "particle.h"
#include "rt_global.h"
#include "rt_transfer.h"
#include "rt_utilities.h"
#include "starformation.h"
#include "timestep.h"
#include "timing.h"
#include "units.h"

#include "F/frt_c.h"

#ifdef RT_DEBUG
#include "rt_debug.h"
#else
int rt_debug = 0;  /* for OpenMP pragmas */
#endif


/*
//  Set the floor for the gast-to-dust ratio. May be needed to force the
//  switch from primordial, metal-free star formation to normal star
//  formation if the resolution is not high enough to resolve first stars
//  properly.
*/ 
float rt_dust_to_gas_floor = 0.001;

/*
//  Efficiencies for stellar and QSO components, relative to the
//  default value of EPSUV=5e-6 (~ Draine field)
*/
float rt_uv_emissivity_stars = 1.0;
float rt_uv_emissivity_quasars  = 0.0;

/*
//  Stellar population
*/
int rt_stellar_pop = 2;

/*
//  Clumping factor
*/
float rt_clumping_factor = 30.0;

/*
//  H2 coherence length (in pc)
*/
float rt_H2_coherence_length = 0.3;


void rtPackCellData(int level, int cell, frt_real rVar[], frt_real rPar[], frt_real *rRadField0, frt_real **pRadField1);
void rtUnPackCellData(int level, int cell, frt_real rVar[], frt_real rPar[], frt_real *rRadField1);


float rt_src_rate;


void rtInitSource(int level)
{
  /* Time the source is on (20 Myr) */
  const float ShiningTime = 2.0e7;

  rt_src_rate = units->time/(ShiningTime*constants->yr);
}


float rtSource(int ipart)
{
  int istar;
  float tPrev;

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
  tPrev = (float)(particle_t[istar]-star_tbirth[istar]-particle_dt[ipart]);
  if(tPrev < 0.0) tPrev = 0.0;

  if(rt_src_rate*tPrev < 100)
    {
      return exp(-rt_src_rate*tPrev)*(1.0-exp(-rt_src_rate*particle_dt[ipart]))/particle_dt[ipart];
    }
  else
    {
      return 0.0;
    }
#else
  return 1.0;
#endif /* PARTICLES && STARFORM */
}


void rtConfigInit()
{
  control_parameter_add2(control_parameter_float,&rt_dust_to_gas_floor,"rt:dust-to-gas-floor","rt_dust_to_gas_floor","the minimum value for the gas-to-dust ratio, in solar (i.e. MW) units");

  control_parameter_add2(control_parameter_float,&rt_uv_emissivity_stars,"rt:uv-emissivity-stars","rt_uv_emissivity_stars","the efficiency for the stellar UV component, relative to the default value (Draine field)");

  control_parameter_add2(control_parameter_float,&rt_uv_emissivity_quasars,"rt:uv-emissivity-quasars","rt_uv_emissivity_quasars","the efficiency for the quasar-like UV component, relative to the default value (Draine field)");

  control_parameter_add2(control_parameter_int,&rt_stellar_pop,"rt:stellar-pop","rt_stellar_pop","the dominant stellar populatio (2 or 3). This value is used to set the spectral shape for the stellar component.");

  control_parameter_add(control_parameter_float,&rt_clumping_factor,"@rt:clumping-factor","the clumping factor of the neutral gas");

  control_parameter_add(control_parameter_float,&rt_H2_coherence_length,"@rt:H2-coherence-length","the coherence length of molecular gas (in parsecs)");
}


void rtConfigVerify()
{
  cart_assert(!(rt_dust_to_gas_floor < 0.0));

  cart_assert(!(rt_uv_emissivity_stars < 0.0));

  cart_assert(!(rt_uv_emissivity_quasars < 0.0));

  cart_assert(rt_stellar_pop==2 || rt_stellar_pop==3);

  cart_assert(!(rt_clumping_factor < 1.0));

  cart_assert(rt_H2_coherence_length > 0.0);
}


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
  frt_real rTime, rVar[frtVAR_DIM], rPar[frtPAR_DIM];
#ifdef RT_TRANSFER
  frt_real rBuffer[rt_num_frequencies];
  frt_real rRadField0[2];
#else
  frt_real *rBuffer = 0, *rRadField0 = 0;
#endif
  frt_real *rRadField1 = rBuffer;
  frt_intg info;

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
      rPar[frtPAR_CELL] = cell_size[level];
#endif /* RT_TRANSFER && RT_TRANSFER_FLUX_CONSERVING */

      /*
      //    Sobolev length 
      */
#ifdef RT_CHEMISTRY
      rtGetSobolevFactors(cell,level,&soblen,&sobvel);
      rPar[frtPAR_SOBL] = soblen;
      rPar[frtPAR_NUMF] = sobvel*rTime/cell_size[level];
#endif

#ifdef RT_DEBUG
      if(rt_debug.Mode==1 && cell==cell_find_position(rt_debug.Pos))
	{
	  rPar[frtPAR_DEB] = 10 + (rt_debug.Stop>1 ? 1 : 0);
	  cart_debug("In cell-level debug for cell %d#%d",cell,cell_level(cell));
	}
      else
	{
	  rPar[frtPAR_DEB] = 0.0;
	}
#endif

      /*
      //  Call the Fortran worker 
      */
      frtCall(cooloff)(rPar,rRadField0,rRadField1,&rTime,rVar,&info);

      rtUnPackCellData(level,cell,rVar,rPar,rRadField1);
    }
}


void rtInitRun()
{
#ifdef RT_TEST
  frt_real Yp = 1.0e-10;
#else
  frt_real Yp = constants->Yp;
#endif
  frt_real Tmin = gas_temperature_floor;
  frt_real D2Gmin = rt_dust_to_gas_floor;
  frt_real ClumpH2 = rt_clumping_factor;
  frt_real CohLenH2 = rt_H2_coherence_length;
  frt_real fGal = rt_uv_emissivity_stars;
  frt_real fQSO = rt_uv_emissivity_quasars;
  frt_intg IPOP = rt_stellar_pop;
  frt_intg IREC = 0;

  rtuInitRun();

#ifdef RT_TEST
  frtCall(setrun).Yp = 1.0e-10;
#else
  frtCall(setrun).Yp = constants->Yp;
#endif
  frtCall(setrun).Tmin = gas_temperature_floor;
  frtCall(setrun).D2Gmin = rt_dust_to_gas_floor;
  frtCall(setrun).ClumpH2 = rt_clumping_factor;
  frtCall(setrun).CohLenH2 = rt_H2_coherence_length;
  frtCall(setrun).fGal = rt_uv_emissivity_stars;
  frtCall(setrun).fQSO = rt_uv_emissivity_quasars;
  frtCall(setrun).IPOP = rt_stellar_pop;
  frtCall(setrun).IREC = 0;
  
  //frtCall(initrun)();
  frtCall(initrun2)(&Yp,&Tmin,&D2Gmin,&ClumpH2,&CohLenH2,&fGal,&fQSO,&IPOP,&IREC);

#ifdef RT_TRANSFER
  rtInitRunTransfer();
#endif
}


/* This function can be called more than once per top level step */ 
void rtUpdateTables()
{
  start_time(RT_TABLES_TIMER);

  /* Fill in the tables */
  frtCall(updatetables)();

  end_time(RT_TABLES_TIMER);
}


void rtStepBegin()
{
  frt_real uDen = units->number_density;
  frt_real uLen = units->length;
  frt_real uTime = units->time;
  frt_real dtStep = dtl[min_level];

#ifdef COSMOLOGY
  frt_real aExp = abox[min_level];
  frt_real Hubble = log(abox_from_tcode(tl[0]+dtl[0])/abox_from_tcode(tl[0]))/dtl[0];
#else
  frt_real aExp = 1.0;
  frt_real Hubble = 0.0;
#endif

  frtCall(stepbegin)(&uDen,&uLen,&uTime,&dtStep,&aExp,&Hubble);

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
#endif /* RT_DEBUG */
}


void rtStepEnd()
{
#ifdef RT_DEBUG
  int i;
#endif
  frt_real rVol = num_root_cells;
  frt_real rAvg[3];
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

  frtCall(stepend)(&rVol,rAvg);

#ifdef RT_TRANSFER
  rtStepEndTransfer();
#endif /* RT_TRANSFER */

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
void rtPackCellData(int level, int cell, frt_real rVar[], frt_real rPar[], frt_real *rRadField0, frt_real **pRadField1)
{
  int i;

  /*
  //  Default values of all parameters is 0
  */
  for(i=0; i<frtPAR_DIM; i++) rPar[i] = 0.0;

  /*
  //  Set parameters:
  //
  //    Density in code units
  */
  rPar[frtPAR_RHO] = cell_gas_density(cell);
  /*
  //    Cell volume in code units
  */
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
  rPar[frtPAR_VOL] = cell_volume[level];
#endif
  /*
  //    Metallicity in units of solar
  */
#ifdef ENRICH
  rPar[frtPAR_ZSOL] = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  rPar[frtPAR_ZSOL] = 0.0;
#endif

  /* 
  //  Pack elemental abundances 
  */
  rVar[frtVAR_Ein] = units->temperature*cell_gas_internal_energy(cell)/rPar[frtPAR_RHO];
  rVar[frtVAR_XHI] = cell_HI_density(cell)/rPar[frtPAR_RHO];
  rVar[frtVAR_XHII] = cell_HII_density(cell)/rPar[frtPAR_RHO];
  rVar[frtVAR_XHeI] = cell_HeI_density(cell)/rPar[frtPAR_RHO];
  rVar[frtVAR_XHeII] = cell_HeII_density(cell)/rPar[frtPAR_RHO];
  rVar[frtVAR_XHeIII] = cell_HeIII_density(cell)/rPar[frtPAR_RHO];
  rVar[frtVAR_XH2] = cell_H2_density(cell)/rPar[frtPAR_RHO];
#ifdef RT_8SPECIES
  rVar[frtVAR_XH2p] = rVar[frtVAR_XHm] = 0.0;
#endif

  /* 
  //  Pack radiation field 
  */
#ifdef RT_TRANSFER
  if(sizeof(frt_real) != sizeof(float))  /* Optimization */
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
  rRadField0[1] = rtGlobals[RT_OT_FIELD_AVG].Value;
#else
  rRadField0[0] = rRadField0[1] = 0.0;
#endif
#endif /* RT_TRANSFER */
}


void rtUnPackCellData(int level, int cell, frt_real rVar[], frt_real rPar[], frt_real *rRadField1)
{
#ifdef RT_DEBUG
  int j, fail;
#endif

  /*
  //  Unpack radiation field (if needed) 
  */
#if defined(RT_TRANSFER) && defined(RT_VARIABLE_RFIELD)
  int i;
  if(sizeof(frt_real) != sizeof(float))  /* Optimization */
    {
      for(i=0; i<rt_num_frequencies; i++)
	{
	  cell_var(cell,rt_freq_offset+i) = rRadField1[i];
	}
    }
#endif

#ifdef RT_DEBUG
  for(fail=j=0; j<frtVAR_DIM; j++)
    {
      if(isnan(rVar[j])) fail = 1;
    }

  if(fail)
    {
      for(j=0; j<frtVAR_DIM; j++)
	{
	  cart_debug("Var[%d] = %g",j,rVar[j]);
	}
      cart_error("frtCoolOff returned NaN");
    }
#endif
  
  /*
  //  Unpack elemental abundances 
  */
  cell_gas_internal_energy(cell) = max(gas_temperature_floor,rVar[frtVAR_Ein])*rPar[frtPAR_RHO]/units->temperature;
  cell_gas_energy(cell) = cell_gas_kinetic_energy(cell) + cell_gas_internal_energy(cell);
  cell_HI_density(cell) = rVar[frtVAR_XHI]*rPar[frtPAR_RHO];
  cell_HII_density(cell) = rVar[frtVAR_XHII]*rPar[frtPAR_RHO];
  cell_HeI_density(cell) = rVar[frtVAR_XHeI]*rPar[frtPAR_RHO];
  cell_HeII_density(cell) = rVar[frtVAR_XHeII]*rPar[frtPAR_RHO];
  cell_HeIII_density(cell) = rVar[frtVAR_XHeIII]*rPar[frtPAR_RHO];
  cell_H2_density(cell) = rVar[frtVAR_XH2]*rPar[frtPAR_RHO];

#ifndef RT_MONOATOMIC
#ifdef RT_HIGH_DENSITY
  cell_gas_gamma(cell) = frtCall(gamma)(Ein,rVar,rPar+1);
#else
  cell_gas_gamma(cell) = frtCall(gamma)(Ein,rVar);
#endif
#endif
}


/*
//  Used in testing
*/
void rtModifyTimeStep(double *dt)
{
#ifdef RT_TEST
  frt_real rdt;
  frt_intg ins;

  rdt = *dt;
  ins = step;

  frtCall(testmodifytimestep)(&rdt,&ins);

  *dt = rdt;
#endif
}


/*
//  Helper function to get gas temperature
*/
float rtTem(int cell)
{
  frt_real buffer[7];

  if(sizeof(frt_real) == sizeof(float))  /* Optimization */
    {
      return frtCall(tem)(&(cell_gas_internal_energy(cell)),cell_vars[cell]+RT_HVAR_OFFSET-1);
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
      return frtCall(tem)(buffer+0,buffer);
    }
}


float rtDustToGas(int cell)
{
  int i;
  /* 
  //  Specify types for the Fortran interface
  */
  frt_real rPar[frtPAR_DIM];

  /*
  //  Default values of all parameters is 0
  */
  for(i=0; i<frtPAR_DIM; i++) rPar[i] = 0.0;

  /*
  //  Set parameters:
  //
  //    Density in code units
  */
  rPar[frtPAR_RHO] = cell_gas_density(cell);
  /*
  //    Cell volume in code units
  */
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
  rPar[frtPAR_VOL] = cell_volume[level];
#endif
  /*
  //    Metallicity in units of solar
  */
#ifdef ENRICH
  rPar[frtPAR_ZSOL] = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  rPar[frtPAR_ZSOL] = 0.0;
#endif

  /*
  //  Call Fortran worker
  */
  return frtCall(dusttogas)(rPar);
}


#ifdef RT_DEBUG
void rtPrintValue(const char *name, float val, int pDen, int pLen, int pTime)
{
  frt_real fval = val;
  frt_intg fpDen = pDen;
  frt_intg fpLen = pLen;
  frt_intg fpTime = pTime;

  val = frtCall(logwithunits)(&fval,&fpDen,&fpLen,&fpTime);
  cart_debug("Checking: lg(%s) = %g",name,val);
}
#endif


void rtGetPhotoRates(int cell, float rate[])
{
  int i, level;

  /* 
  //  Specify types for the Fortran interface
  */
  frt_real rVar[frtVAR_DIM], rPar[frtPAR_DIM], pRate[frtRATE_DIM];
#ifdef RT_TRANSFER
  frt_real rBuffer[rt_num_frequencies];
  frt_real rRadField0[2];
#else
  frt_real *rBuffer = 0, *rRadField0 = 0;
#endif
  frt_real rf[1+2*rt_num_frequencies];
  frt_real *rRadField1 = rBuffer;
  frt_intg iTab[2];

  if(cell < 0)
    {
      if(sizeof(frt_real) == sizeof(float))
	{
	  frtCall(getbackgroundphotorates)((frt_real *)rate);
	}
      else
	{
	  frtCall(getbackgroundphotorates)(pRate);
	  for(i=0; i<frtRATE_DIM; i++) rate[i] = pRate[i];
	}
    }
  else
    {
      level = cell_level(cell);
      rtPackCellData(level,cell,rVar,rPar,rRadField0,&rRadField1);

      /*
      //  Call Fortran workers
      */
#ifdef RT_TRANSFER
      frtCall(transferpackradiationfield)(rPar,rVar,rRadField0,rRadField1,rf);
#endif

      if(sizeof(frt_real) == sizeof(float))
	{
	  frtCall(getphotorates)(rPar,rf,iTab,rVar,(frt_real *)rate);
	}
      else
	{
	  frtCall(getphotorates)(rPar,rf,iTab,rVar,pRate);
	  for(i=0; i<frtRATE_DIM; i++) rate[i] = pRate[i];
	}
    }
}


void rtGetBinIds(int n, const float wlen[], int idxi[])
{
  int i;
  frt_real w;

  for(i=0; i<n; i++)
    {
      w = wlen[i];
      idxi[i] = frtCall(getbinid)(&w);
    }
}


void rtGetBinWavelengths(int n, const int idxi[], float wlen[])
{
  int i;
  frt_intg lr; 

  for(i=0; i<n; i++)
    {
      lr = idxi[i];
      wlen[i] = frtCall(getbinwavelength)(&lr);
    }
}


void rtGetRadiationField(int cell, int n, const int idxi[], float ngxi[])
{
  int i, level;

  /* 
  //  Specify types for the Fortran interface
  */
  frt_real rVar[frtVAR_DIM], rPar[frtPAR_DIM];
#ifdef RT_TRANSFER
  frt_real rBuffer[rt_num_frequencies];
  frt_real rRadField0[2];
#else
  frt_real *rBuffer = 0, *rRadField0 = 0;
#endif
  frt_real rf[1+2*rt_num_frequencies];
  frt_real *rRadField1 = rBuffer;
  frt_intg lr;

  if(cell < 0)
    {
      for(i=0; i<n; i++)
	{
	  lr = idxi[i];
	  ngxi[i] = frtCall(getbackgroundradiationfield)(&lr);
	}
    }
  else
    {
      level = cell_level(cell);
      frtCall(transferpackradiationfield)(rPar,rVar,rRadField0,rRadField1,rf);

      rtPackCellData(level,cell,rVar,rPar,rRadField0,&rRadField1);

      /*
      //  Call Fortran workers
      */
#ifdef RT_TRANSFER
      frtCall(transferpackradiationfield)(rPar,rVar,rRadField0,rRadField1,rf);
#endif

      for(i=0; i<n; i++)
	{
	  lr = idxi[i];
	  ngxi[i] = frtCall(getradiationfield)(&lr,rf);
	}
    }
}

#endif /* RADIATIVE_TRANSFER */

