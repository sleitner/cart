#include "config.h"
#ifdef RADIATIVE_TRANSFER

#include <math.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rt.h"
#include "rt_debug.h"
#include "rt_global.h"
#include "rt_transfer.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "frt/frt_c.h"


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
float rt_coherence_length = 0.3;

/*
//  Switch to limit the signal propagation speed to c
*/
int rt_limit_signal_speed_to_c = 1;


void rtPackCellData(int level, int cell, frt_real *var, frt_real **p_rawrf);
void rtUnPackCellData(int level, int cell, frt_real *var, frt_real *rawrf);


void rtSetupSource(int level)
{
#if defined(STAR_FORMATION) && !defined(RT_TEST)
  if(sf_feedback_particle->setup != NULL) sf_feedback_particle->setup(level);
#endif 
}


float rtConstantSource(int ipart)
{
  return 1.0;
}


float (*rtSource)(int ipart) = NULL;


void rtConfigInit()
{
#if defined(STAR_FORMATION) && !defined(RT_TEST)
  rtSource = sf_feedback_particle->rt_source;
#else
  rtSource = rtConstantSource;
#endif

  control_parameter_add2(control_parameter_float,&rt_dust_to_gas_floor,"rt:dust-to-gas-floor","rt_dust_to_gas_floor","the minimum value for the gas-to-dust ratio, in solar (i.e. MW) units.");

  control_parameter_add2(control_parameter_float,&rt_uv_emissivity_stars,"rt:uv-emissivity-stars","rt_uv_emissivity_stars","the efficiency for the stellar UV component, relative to the default value (Draine field).");

  control_parameter_add2(control_parameter_float,&rt_uv_emissivity_quasars,"rt:uv-emissivity-quasars","rt_uv_emissivity_quasars","the efficiency for the quasar-like UV component, relative to the default value (Draine field).");

  control_parameter_add2(control_parameter_int,&rt_stellar_pop,"rt:stellar-pop","rt_stellar_pop","the dominant stellar populatio (2 or 3). This value is used to set the spectral shape for the stellar component.");

  control_parameter_add2(control_parameter_bool,&rt_limit_signal_speed_to_c,"rt:limit-signal-speed-to-c","rt_limit_signal_speed_to_c","if set, limits the signal propagation speed to the spped of light.");

  control_parameter_add(control_parameter_float,&rt_clumping_factor,"@rt:clumping-factor","the clumping factor of the neutral gas.");

  control_parameter_add(control_parameter_float,&rt_coherence_length,"@rt:coherence-length","the coherence length of molecular gas (in parsecs).");

#ifdef RT_TRANSFER 
  rtConfigInitTransfer();
#endif
}


void rtConfigVerify()
{
  VERIFY(rt:dust-to-gas-floor, !(rt_dust_to_gas_floor < 0.0) );

  VERIFY(rt:uv-emissivity-stars, !(rt_uv_emissivity_stars < 0.0) );

  VERIFY(rt:uv-emissivity-quasars, !(rt_uv_emissivity_quasars < 0.0) );

  VERIFY(rt:stellar-pop, rt_stellar_pop==2 || rt_stellar_pop==3 );

  VERIFY(@rt:clumping-factor, !(rt_clumping_factor < 1.0) );

  VERIFY(@rt:coherence-length, rt_coherence_length > 0.0 );

#ifdef RT_TRANSFER 
  rtConfigVerifyTransfer();
#endif
}


void rtInitRun()
{
#ifdef RT_TEST
  frt_real Yp = 1.0e-10;
#else
  frt_real Yp = constants->Yp;
#endif
  frt_real Tmin = gas_temperature_floor;
  frt_real D2GminH2 = rt_dust_to_gas_floor*(constants->Dsun/constants->Zsun);
  frt_real CluFacH2 = rt_clumping_factor;
  frt_real CohLenH2 = rt_coherence_length;
  frt_real fGal = rt_uv_emissivity_stars;
  frt_real fQSO = rt_uv_emissivity_quasars;
  frt_intg IPOP = rt_stellar_pop;
  frt_intg IREC = 0;
  frt_intg IOUNIT = 81;

  start_time(WORK_TIMER);

  frtCall(initrun2)(&Yp,&Tmin,&D2GminH2,&CluFacH2,&CohLenH2,&fGal,&fQSO,&IPOP,&IREC,&IOUNIT);

  end_time(WORK_TIMER);

#ifdef RT_TRANSFER
  rtInitRunTransfer();
#endif
}


/* This function can be called more than once per top level step */ 
void rtUpdateTables()
{
#ifdef RT_TRANSFER
  int i;
  frt_real rfAvg[rt_num_fields];

  for(i=0; i<rt_num_fields; i++) rfAvg[i] = rtAvgRF[i].Value;

#else
  frt_real *rfAvg = NULL;
#endif

  start_time(RT_TABLES_TIMER);
  start_time(WORK_TIMER);

  /* Fill in the tables */
  frtCall(updatetables)(rfAvg);

  end_time(WORK_TIMER);
  end_time(RT_TABLES_TIMER);
}


void rtInitStep(double dt)
{
  frt_real uDen = units->number_density;
  frt_real uLen = units->length;
  frt_real uTime = units->time;
  frt_real dtCode = dt;

#ifdef COSMOLOGY
  frt_real aExp = abox[min_level];
  frt_real daExp = abox_from_tcode(tl[min_level]+dt) - abox[min_level];
  frt_real HExp = Hubble(aExp)/units->time;
#else
  frt_real aExp = 1.0;
  frt_real daExp = 0.0;
  frt_real HExp = 0.0;
#endif

  start_time(WORK_TIMER);

  frtCall(stepbegin)(&uDen,&uLen,&uTime,&aExp,&HExp,&daExp,&dtCode);

  end_time(WORK_TIMER);

#ifdef RT_TRANSFER
  rtInitStepTransfer();
#endif

  rtUpdateTables();

#ifdef RT_DEBUG
  switch(rt_debug.Mode)
    {
    case 1:
      {
	int i, cell;
	cell = cell_find_position(rt_debug.Pos);
	cart_debug("In cell-level debug for cell %d/%d",cell,cell_level(cell));
	cart_debug("RT_HVAR_OFFSET: %d",RT_HVAR_OFFSET);
#ifdef RT_VAR_SOURCE
	cart_debug("RT_VAR_SOURCE: %d",RT_VAR_SOURCE);
#endif
	cart_debug("rt_grav_vars_offset: %d",rt_grav_vars_offset);
#ifdef RT_TRANSFER
	cart_debug("rt_num_vars: %d",rt_num_vars);
#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
	cart_debug("RT_VAR_OT_FIELD: %d",RT_VAR_OT_FIELD);
	cart_debug("rt_et_offset: %d",rt_et_offset);
	cart_debug("rt_field_offset: %d",rt_field_offset);
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


void rtGlobalUpdate(int top_level, MPI_Comm level_com)
{
  start_time(RT_GLOBAL_UPDATE_TIMER);

#ifdef RT_TRANSFER
  rtGlobalUpdateTransfer(top_level,level_com);
#endif
  
  //rtUpdateTables(top_level,level_com);

  end_time(RT_GLOBAL_UPDATE_TIMER);
}


void rtAfterAssignDensity1(int level)
{
  int num_level_cells;
  int *level_cells;

  start_time(RT_AFTER_DENSITY_TIMER);

#ifdef RT_TRANSFER
  /* assumes buffer gas density is up to date */
  start_time( WORK_TIMER );
  select_level(level,CELL_TYPE_ANY,&num_level_cells,&level_cells);
  end_time( WORK_TIMER );

  rtAfterAssignDensityTransfer(level,num_level_cells,level_cells);

  start_time( WORK_TIMER );
  cart_free(level_cells);
  end_time( WORK_TIMER );
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


/*
//  Helper functions for Fortran workers
*/
void rtPackCellData(int level, int cell, frt_real *var, frt_real **p_rawrf)
{
  int i;
#ifdef RT_CHEMISTRY
  float sobvel;
#endif

  frtCall(initvar)(var);

  /* 
  //  Pack elemental abundances 
  */
  var[FRT_Ein] = units->temperature*cell_gas_internal_energy(cell)/cell_gas_density(cell);
  var[FRT_XHI] = cell_HI_density(cell)/cell_gas_density(cell);
  var[FRT_XHII] = cell_HII_density(cell)/cell_gas_density(cell);
  var[FRT_XHeI] = cell_HeI_density(cell)/cell_gas_density(cell);
  var[FRT_XHeII] = cell_HeII_density(cell)/cell_gas_density(cell);
  var[FRT_XHeIII] = cell_HeIII_density(cell)/cell_gas_density(cell);
  var[FRT_XH2] = cell_H2_density(cell)/cell_gas_density(cell);
#ifdef RT_8SPECIES
#error "The 8-species model is not currently implemented."
#endif

  /*
  //  Pack parameters (in the order they appear, to improve cache performance):
  //  Density in code units
  */
  var[FRT_Density] = cell_gas_density(cell);

  /*
  //  Metallicity in units of solar
  */
#ifdef ENRICHMENT
  var[FRT_Metallicity] = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  var[FRT_Metallicity] = 0.0;
#endif

  /*
  //  Dust-to-gas ratio
  */
#ifdef DUST_EVOLUTION
  var[FRT_DustToGas] = cell_dust_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  var[FRT_DustToGas] = constants->Dsun/constants->Zsun*var[FRT_Metallicity];
#endif

  /*
  //  Sobolev length 
  */
#ifdef RT_CHEMISTRY
  var[FRT_SobolevLength] = cell_sobolev_length2(cell,level,&sobvel);
  var[FRT_NumericalDiffusionFactor] = sobvel/cell_size[level];
#endif

  /* 
  //  OT radiation fields
  */
#ifdef RT_VAR_OT_FIELD  
  var[FRT_OTRadiationFieldLocal] = cell_var(cell,RT_VAR_OT_FIELD);
  var[FRT_OTRadiationFieldGlobal] = 1.0;
#endif

  /*
  //    Cell volume in code units
  */
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
  var[FRT_ResolutionElementVolume] = cell_volume[level];
#endif

  /*
  //    Cell size for flux-conserving correction
  */
#if defined(RT_TRANSFER) && defined(RT_TRANSFER_FLUX_CONSERVING)
  var[FRT_ResolutionElementSize] = cell_size[level];
#endif /* RT_TRANSFER && RT_TRANSFER_FLUX_CONSERVING */

  /* 
  //  Pack radiation field (if requested)
  */
#ifdef RT_TRANSFER
  if(p_rawrf != NULL)
    {
      if(sizeof(frt_real) != sizeof(float))  /* Optimization */
	{
	  for(i=0; i<rt_num_fields; i++)
	    {
	      (*p_rawrf)[i] = cell_var(cell,rt_field_offset+i);
	    }
	}
      else
	{
	  *p_rawrf = &(cell_var(cell,rt_field_offset));
	}
    }
#endif /* RT_TRANSFER */
}


void rtUnPackCellData(int level, int cell, frt_real *var, frt_real *rawrf)
{
#ifdef RT_DEBUG
  int j, fail;
#endif

  /*
  //  Unpack radiation field (if needed) 
  */
#if defined(RT_TRANSFER) && defined(RT_VARIABLE_RF)
  int i;
  if(sizeof(frt_real) != sizeof(float))  /* Optimization */
    {
      for(i=0; i<rt_num_fields; i++)
	{
	  cell_var(cell,rt_field_offset+i) = rawrf[i];
	}
    }
#endif

#ifdef RT_DEBUG
  for(fail=j=0; j<FRT_Debug; j++)
    {
      if(isnan(var[j])) fail = 1;
    }

  if(fail)
    {
      for(j=0; j<FRT_Debug; j++)
	{
	  cart_debug("Var[%d] = %g",j,var[j]);
	}
      double pos[nDim];
      cell_center_position(cell,pos);
      cart_error("frtCoolOff returned NaN, cell = %d, pos = (%lf,%lf,%lf)",cell,pos[0],pos[1],pos[2]);
    }
#endif
  
  /*
  //  Unpack elemental abundances 
  */
  cell_gas_internal_energy(cell) = MAX(gas_temperature_floor,var[FRT_Ein])*cell_gas_density(cell)/units->temperature;
  cell_gas_energy(cell) = cell_gas_kinetic_energy(cell) + cell_gas_internal_energy(cell);
  cell_HI_density(cell) = var[FRT_XHI]*cell_gas_density(cell);
  cell_HII_density(cell) = var[FRT_XHII]*cell_gas_density(cell);
  cell_HeI_density(cell) = var[FRT_XHeI]*cell_gas_density(cell);
  cell_HeII_density(cell) = var[FRT_XHeII]*cell_gas_density(cell);
  cell_HeIII_density(cell) = var[FRT_XHeIII]*cell_gas_density(cell);
  cell_H2_density(cell) = var[FRT_XH2]*cell_gas_density(cell);
#ifdef RT_8SPECIES
#error "The 8-species model is not currently implemented."
#endif

  /*
  //  Dust-to-gas ratio (only changes inside the FRT block if
  //  RT_DUST_EVOLUTION is set)
  */
#if defined(DUST_EVOLUTION) && defined(RT_DUST_EVOLUTION)
  cell_dust_density(cell) = var[FRT_DustToGas]*constants->Zsun*cell_gas_density(cell);
#endif

  cell_gas_gamma(cell) = var[FRT_Gamma];
}


/*
//  Helper function to get gas temperature
*/
float rtTem(int cell)
{
  frt_real var[FRT_DIM];

  rtPackCellData(cell_level(cell),cell,var,NULL);
  return frtCall(tem)(var)/units->temperature;
}


float rtDmw(int cell)
{
#ifdef DUST_EVOLUTION
  return cell_dust_density(cell)/(constants->Dsun*cell_gas_density(cell));
#else  /* DUST_EVOLUTION */
#ifdef ENRICHMENT
  return cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  return 0.0;
#endif
#endif /* DUST_EVOLUTION */
}


float rtDmwFL(int cell)
{
#ifdef RT_FIXED_ISM
  return rt_dust_to_gas_floor;
#else
  float d = rtDmw(cell);
  return MAX(rt_dust_to_gas_floor,d);
#endif
}


/*
//  UV field at 12.0eV in units of Draine field (1.0e6 phot/cm^2/s/ster/eV)
*/
void rtGetPhotoRates(int cell, float *rate);
float rtUmw(int cell)
{
  float rate[FRT_RATE_DIM];
  rtGetPhotoRates(cell,rate);
  return rate[FRT_RATE_DissociationLW]*1.05e10;
}


void rtGetPhotoRatesFS(int cell, float *rate);
float rtUmwFS(int cell)
{
  float rate[FRT_RATE_DIM];
  rtGetPhotoRatesFS(cell,rate);
  return rate[FRT_RATE_DissociationLW]*1.05e10;
}


void rtGetCoolingRate(int cell, float *cooling_rate, float *heating_rate)
{
  DEFINE_FRT_INTEFACE(var,rawrf);
  frt_real c, h;

  rtPackCellData(cell_level(cell),cell,var,&rawrf);

#ifdef RT_FIXED_ISM
  var[FRT_Metallicity] = rt_dust_to_gas_floor;
#endif

  /*
  //  Call the Fortran worker 
  */
  frtCall(getcoolingrate)(var,rawrf,&c,&h);

  *cooling_rate = c;
  *heating_rate = h;
}


void rtGetPhotoRatesWorker(int cell, float *rate, int free_space)
{
  int i;
  DEFINE_FRT_INTEFACE(var,rawrf);
  frt_real frate[FRT_RATE_DIM];

  if(cell < 0)
    {
      if(sizeof(frt_real) == sizeof(float))
	{
	  for(i=0; i<FRT_RATE_DIM; i++) rate[i] = 0.0;
	  frtCall(getbackgroundphotorates)((frt_real *)rate);
	}
      else
	{
	  for(i=0; i<FRT_RATE_DIM; i++) frate[i] = 0.0;
	  frtCall(getbackgroundphotorates)(frate);
	  for(i=0; i<FRT_RATE_DIM; i++) rate[i] = frate[i];
	}
    }
  else
    {
      rtPackCellData(cell_level(cell),cell,var,&rawrf);

      if(sizeof(frt_real) == sizeof(float))
	{
	  for(i=0; i<FRT_RATE_DIM; i++) rate[i] = 0.0;
	  if(free_space)
	    frtCall(getphotoratesfs)(var,rawrf,(frt_real *)rate);
	  else
	    frtCall(getphotorates)(var,rawrf,(frt_real *)rate);
	}
      else
	{
	  for(i=0; i<FRT_RATE_DIM; i++) frate[i] = 0.0;
	  if(free_space)
	    frtCall(getphotoratesfs)(var,rawrf,frate);
	  else
	    frtCall(getphotorates)(var,rawrf,frate);
	  for(i=0; i<FRT_RATE_DIM; i++) rate[i] = frate[i];
	}
    }
}


void rtGetPhotoRates(int cell, float *rate)
{
  rtGetPhotoRatesWorker(cell,rate,0);
}


void rtGetPhotoRatesFS(int cell, float *rate)
{
  rtGetPhotoRatesWorker(cell,rate,1);
}


void rtGetRadiationFieldWorker(int cell, int n, const float *wlen, float *ngxi, int free_space)
{
  int i;
  DEFINE_FRT_INTEFACE(var,rawrf);
  frt_real *fwlen, *fngxi;
  frt_intg nout = n;

  cart_assert(n > 0);

  if(cell < 0)
    {
      if(sizeof(frt_real) == sizeof(float))
	{
	  frtCall(getbackgroundradiationfield)(&nout,(frt_real *)wlen,(frt_real*)ngxi);
	}
      else
	{
	  fwlen = cart_alloc(frt_real,n);
	  fngxi = cart_alloc(frt_real,n);
	  for(i=0; i<n; i++) fwlen[i] = wlen[i];
	  frtCall(getbackgroundradiationfield)(&nout,fwlen,fngxi);
	  for(i=0; i<n; i++) ngxi[i] = fngxi[i];
	  cart_free(fwlen);
	  cart_free(fngxi);
	}
    }
  else
    {
      rtPackCellData(cell_level(cell),cell,var,&rawrf);

      if(sizeof(frt_real) == sizeof(float))
	{
	  if(free_space)
	    frtCall(getradiationfieldfs)(var,rawrf,&nout,(frt_real *)wlen,(frt_real*)ngxi);
	  else
	    frtCall(getradiationfield)(var,rawrf,&nout,(frt_real *)wlen,(frt_real*)ngxi);
	}
      else
	{
	  fwlen = cart_alloc(frt_real,n);
	  fngxi = cart_alloc(frt_real,n);
	  for(i=0; i<n; i++) fwlen[i] = wlen[i];
	  if(free_space)
	    frtCall(getradiationfieldfs)(var,rawrf,&nout,fwlen,fngxi);
	  else
	    frtCall(getradiationfield)(var,rawrf,&nout,fwlen,fngxi);
	  for(i=0; i<n; i++) ngxi[i] = fngxi[i];
	  cart_free(fwlen);
	  cart_free(fngxi);
	}
    }
}


void rtGetRadiationField(int cell, int n, const float *wlen, float *ngxi)
{
  rtGetRadiationFieldWorker(cell,n,wlen,ngxi,0);
}


void rtGetRadiationFieldFS(int cell, int n, const float *wlen, float *ngxi)
{
  rtGetRadiationFieldWorker(cell,n,wlen,ngxi,1);
}

#endif /* RADIATIVE_TRANSFER */

