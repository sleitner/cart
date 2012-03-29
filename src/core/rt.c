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


double rt_rp_amplt = 1.0;
double rt_rp_slope = 1.0;


void rtPackCellData(int level, int cell, frt_real *var, frt_real **p_rawrf);
void rtUnPackCellData(int level, int cell, frt_real *var, frt_real *rawrf);

void frtCall(getrpfactors)(frt_real *rf2Prad, frt_real *rho2abc);


float rt_src_rate;
struct rtRadiationPressureFactors
{
  float rf2Prad;
  float cd2tauUV;
  float cd2tauIR;
  float cd2sigma100;
}
rt_rp_factor;


float cell_radiation_pressure(int cell)
{
  float len = cell_sobolev_length(cell);
  float cd_gas = len*constants->XH*cell_gas_density(cell);
  float cd_dust = len*constants->XH*rtDmw(cell)*(cell_HI_density(cell)+2*cell_H2_density(cell));
  float tauUV = rt_rp_factor.cd2tauUV*cd_dust;

#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD==RT_METHOD_OTVET)
  float rfLoc;

#ifdef RT_UV
  float w;
  if(cell_var(cell,rt_field_offset+rt_num_freqs-1)>0.0 && cell_var(cell,rt_field_offset+rt_num_freqs-1)>1.0e-35*cell_var(cell,RT_VAR_OT_FIELD))
    {
      w = log(cell_var(cell,RT_VAR_OT_FIELD)/cell_var(cell,rt_field_offset+rt_num_freqs-1));
      if(tauUV > w) tauUV = w;
    }
  rfLoc = cell_var(cell,rt_field_offset+rt_num_freqs-1)*exp(tauUV);
#else /* RT_UV */
  rfLoc = cell_var(cell,RT_VAR_OT_FIELD);
#endif /* RT_UV */
  
  return (1-exp(-tauUV)+rt_rp_amplt*pow(rt_rp_factor.cd2sigma100*cd_gas,rt_rp_slope))*rt_rp_factor.rf2Prad*rfLoc;

#else /* RT_TRANSFER && RT_TRANSFER_METHOD==RT_METHOD_OTVET */
  cart_error("Radiation pressure without RT_TRANSFER and RT_TRANSFER_METHOD=RT_METHOD_OTVET is not implemented yet.");
  return 0.0;
#endif /* RT_TRANSFER && RT_TRANSFER_METHOD==RT_METHOD_OTVET */
}


void rtInitSource(int level)
{
  /* Time the source is on (20 Myr) */
#ifdef RT_OLDSTYLE_SOURCE_FUNCTION
  const float ShiningTime = 2.0e7;
#else
  const float ShiningTime = 3.0e6;
#endif /* RT_OLDSTYLE_SOURCE_FUNCTION */

  rt_src_rate = units->time/(ShiningTime*constants->yr);
}


float rtSource(int ipart)
{
  int istar;
#if defined(PARTICLES) && defined(STAR_FORMATION)
#ifdef RT_OLDSTYLE_SOURCE_FUNCTION
  float x1;
#else
  float x1, x2, dx;
#endif /* RT_OLDSTYLE_SOURCE_FUNCTION */
#endif /* PARTICLES && STAR_FORMATION */

#if defined(PARTICLES) && defined(STAR_FORMATION)
  if(!particle_is_star(ipart)) return 0.0;
#endif

#ifdef RT_TEST
  return 1.0;
#endif

#if defined(PARTICLES) && defined(STAR_FORMATION)
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
  x1 = rt_src_rate*(particle_t[istar]-star_tbirth[istar]-particle_dt[ipart]);
  if(x1 < 0.0) x1 = 0.0;

#ifdef RT_OLDSTYLE_SOURCE_FUNCTION
  if(x1 < 100)
    {
      return exp(-x1)*(1.0-exp(-rt_src_rate*particle_dt[ipart]))/particle_dt[ipart];
    }
#else
  if(x1 < 1.0e4)
    {
      /*
      //  This is a rough fit to Starburst99 evolving spectra
      */
      dx = rt_src_rate*particle_dt[ipart];
      if(dx > 1.0e-5)
	{
	  x2 = x1 + dx;
	  x1 *= (0.8+x1*x1);
	  x2 *= (0.8+x2*x2);
	  return (x2-x1)/(1+x1)/(1+x2)/particle_dt[ipart];
	}
      else
	{
	  x2 = x1*(0.8+x1*x1);
	  return (0.8+3*x1*x1)/(1+x2)/(1+x2)*rt_src_rate;
	}
    }
#endif /* RT_OLDSTYLE_SOURCE_FUNCTION */
  else
    {
      return 0.0;
    }
#else
  return 1.0;
#endif /* PARTICLES && STAR_FORMATION */
}


void rtConfigInit()
{
  control_parameter_add2(control_parameter_float,&rt_dust_to_gas_floor,"rt:dust-to-gas-floor","rt_dust_to_gas_floor","the minimum value for the gas-to-dust ratio, in solar (i.e. MW) units.");

  control_parameter_add2(control_parameter_float,&rt_uv_emissivity_stars,"rt:uv-emissivity-stars","rt_uv_emissivity_stars","the efficiency for the stellar UV component, relative to the default value (Draine field).");

  control_parameter_add2(control_parameter_float,&rt_uv_emissivity_quasars,"rt:uv-emissivity-quasars","rt_uv_emissivity_quasars","the efficiency for the quasar-like UV component, relative to the default value (Draine field).");

  control_parameter_add2(control_parameter_int,&rt_stellar_pop,"rt:stellar-pop","rt_stellar_pop","the dominant stellar populatio (2 or 3). This value is used to set the spectral shape for the stellar component.");

  control_parameter_add2(control_parameter_int,&rt_limit_signal_speed_to_c,"rt:limit-signal-speed-to-c","rt_limit_signal_speed_to_c","if set, limits the signal propagation speed to the spped of light.");

  control_parameter_add(control_parameter_float,&rt_clumping_factor,"@rt:clumping-factor","the clumping factor of the neutral gas.");

  control_parameter_add(control_parameter_float,&rt_coherence_length,"@rt:coherence-length","the coherence length of molecular gas (in parsecs).");

  control_parameter_add(control_parameter_double,&rt_rp_amplt,"@rt:rp-amplt","temporary control for testing.");
  control_parameter_add(control_parameter_double,&rt_rp_slope,"@rt:rp-slope","temporary control for testing.");
}


void rtConfigVerify()
{
  cart_assert(!(rt_dust_to_gas_floor < 0.0));

  cart_assert(!(rt_uv_emissivity_stars < 0.0));

  cart_assert(!(rt_uv_emissivity_quasars < 0.0));

  cart_assert(rt_stellar_pop==2 || rt_stellar_pop==3);

  cart_assert(rt_limit_signal_speed_to_c==0 || rt_limit_signal_speed_to_c==1);

  cart_assert(!(rt_clumping_factor < 1.0));

  cart_assert(rt_coherence_length > 0.0);
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
  frt_real CohLenH2 = rt_coherence_length;
  frt_real fGal = rt_uv_emissivity_stars;
  frt_real fQSO = rt_uv_emissivity_quasars;
  frt_intg IPOP = rt_stellar_pop;
  frt_intg IREC = 0;
  frt_intg IOUNIT = 81;

  start_time(WORK_TIMER);

#ifdef RT_TEST
  frtCall(setrun).Yp = 1.0e-10;
#else
  frtCall(setrun).Yp = constants->Yp;
#endif
  frtCall(setrun).Tmin = gas_temperature_floor;
  frtCall(setrun).D2Gmin = rt_dust_to_gas_floor;
  frtCall(setrun).ClumpH2 = rt_clumping_factor;
  frtCall(setrun).CohLenH2 = rt_coherence_length;
  frtCall(setrun).fGal = rt_uv_emissivity_stars;
  frtCall(setrun).fQSO = rt_uv_emissivity_quasars;
  frtCall(setrun).IPOP = rt_stellar_pop;
  frtCall(setrun).IREC = 0;
  frtCall(setrun).IOUNIT = 81;
  
  //frtCall(initrun)();
  frtCall(initrun2)(&Yp,&Tmin,&D2Gmin,&ClumpH2,&CohLenH2,&fGal,&fQSO,&IPOP,&IREC,&IOUNIT);

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
  frt_real rf2Prad, rho2abc;

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

  frtCall(stepbegin)(&uDen,&uLen,&uTime,&aExp,&HExp,&daExp);

  end_time(WORK_TIMER);

#ifdef RT_TRANSFER
  rtInitStepTransfer();
#endif

  rtUpdateTables();

  /*
  //  Don't forget that n_xi is in comoving units
  */
  frtCall(getrpfactors)(&rf2Prad,&rho2abc);
  /*
  //  Factor 0.57 is the ratio of Lbol to nu*Lnu at 1000 A, per Oscar Agertz e-mail
  */
  rt_rp_factor.rf2Prad = 0.57*(constants->k)/units->energy_density*rf2Prad;
  rt_rp_factor.cd2tauUV = rho2abc;

  /*
  //  Units for tauIR factor, assuming kIR = 5 cm^2/g, as in Hopkins et al 1101.4940
  */
  rt_rp_factor.cd2tauIR = 5*units->density*units->length;

  rt_rp_factor.cd2sigma100 = units->density*units->length/(100*constants->Msun/constants->pc/constants->pc);

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
  var[FRT_XH2p] = var[FRT_XHm] = 0.0;
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
  cell_gas_internal_energy(cell) = max(gas_temperature_floor,var[FRT_Ein])*cell_gas_density(cell)/units->temperature;
  cell_gas_energy(cell) = cell_gas_kinetic_energy(cell) + cell_gas_internal_energy(cell);
  cell_HI_density(cell) = var[FRT_XHI]*cell_gas_density(cell);
  cell_HII_density(cell) = var[FRT_XHII]*cell_gas_density(cell);
  cell_HeI_density(cell) = var[FRT_XHeI]*cell_gas_density(cell);
  cell_HeII_density(cell) = var[FRT_XHeII]*cell_gas_density(cell);
  cell_HeIII_density(cell) = var[FRT_XHeIII]*cell_gas_density(cell);
  cell_H2_density(cell) = var[FRT_XH2]*cell_gas_density(cell);
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
#ifdef RT_CUSTOM_DUST_TO_GAS
  frt_real var[FRT_DIM];

  rtPackCellData(cell_level(cell),cell,var,NULL);
  return frtCall(dusttogas)(var);
#else  /* RT_CUSTOM_DUST_TO_GAS */
#ifdef ENRICHMENT
  return cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  return 0.0;
#endif
#endif /* RT_CUSTOM_DUST_TO_GAS */
}


float rtDmw2(int cell)
{
#ifdef RT_FIXED_ISM
  return rt_dust_to_gas_floor;
#else
  float d = rtDmw(cell);
  return max(rt_dust_to_gas_floor,d);
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


void rtGetPhotoRates(int cell, float *rate)
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
	  frtCall(getphotorates)(var,rawrf,(frt_real *)rate);
	}
      else
	{
	  for(i=0; i<FRT_RATE_DIM; i++) frate[i] = 0.0;
	  frtCall(getphotorates)(var,rawrf,frate);
	  for(i=0; i<FRT_RATE_DIM; i++) rate[i] = frate[i];
	}
    }
}


void rtGetRadiationField(int cell, int n, const float *wlen, float *ngxi)
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
	  frtCall(getradiationfield)(var,rawrf,&nout,(frt_real *)wlen,(frt_real*)ngxi);
	}
      else
	{
	  fwlen = cart_alloc(frt_real,n);
	  fngxi = cart_alloc(frt_real,n);
	  for(i=0; i<n; i++) fwlen[i] = wlen[i];
	  frtCall(getradiationfield)(var,rawrf,&nout,fwlen,fngxi);
	  for(i=0; i<n; i++) ngxi[i] = fngxi[i];
	  cart_free(fwlen);
	  cart_free(fngxi);
	}
    }
}

#endif /* RADIATIVE_TRANSFER */

