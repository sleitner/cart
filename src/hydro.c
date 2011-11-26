#include "config.h"
#ifdef HYDRO 

#include <math.h>
#include <stdio.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "tree.h"
#include "units.h"


int pressure_floor_min_level = max_level; /* NG: that used to be MinL_Jeans define */
float pressure_floor_factor = 10.0;

int pressureless_fluid_eos = 0;           /* NG: that used to be PRESSURELESS_FLUID define */
int apply_lapidus_viscosity = 1;          /* NG: that used to be LAPIDUS define */
int smooth_density_gradients = 1;         /* NG: that used to be DENSGRADSMOOTH define */

float gas_density_floor = 1e-6;
float gas_temperature_floor = 3.0;        /* NG: that used to be T_min define */
/*
//  Radiation pressure fudge factor (Hopkins' parameter eta_p)
*/
float radiation_pressure_factor = 0.0;    /* by default, RP is off */


#ifdef BLASTWAVE_FEEDBACK
double blastwave_time_floor = 1.0e-30; 
double blastwave_time_cut = 1.0e-20;
#endif /* BLASTWAVE_FEEDBACK */


void config_init_hydro()
{
  control_parameter_add3(control_parameter_float,&gas_temperature_floor,"gas-temperature-floor","gas_temperature_floor","T_min","the minimum temperature for the gas (in K).");

  control_parameter_add2(control_parameter_float,&gas_density_floor,"gas-density-floor","gas_density_floor","the minimum densitye for the gas (in code units).");

  control_parameter_add2(control_parameter_int,&pressureless_fluid_eos,"pressureless-fluid-eos","pressureless_fluid_eos","use a pressureless fluid equation of state.");

  control_parameter_add2(control_parameter_int,&apply_lapidus_viscosity,"apply-lapidus-viscosity","apply_lapidus_viscosity","apply Lapidus viscosity in the hydro flux calculations (boolean value).");

  control_parameter_add2(control_parameter_int,&smooth_density_gradients,"smooth-density-gradients","smooth_density_gradients","smooth overly steep density gradients in the hydro flux calculations (boolean value).");

  control_parameter_add3(control_parameter_int,&pressure_floor_min_level,"pressure-floor-min-level","pressure_floor_min_level","MinL_Jeans","the level to apply the pressure floor. If this value is set to -1, the pressure floor correction is disabled.");

  control_parameter_add(control_parameter_float,&pressure_floor_factor,"@pressure-floor-factor","the factor to scale the pressure floor with. The default, thoroughly tested value is 10. If you change it, make sure you know what you are doing.");

  control_parameter_add2(control_parameter_float,&radiation_pressure_factor,"radiation-pressure-factor","radiation_pressure_factor","Hopkins' eta_p factor to scale the radiation pressure term; set it to zero to disable radiation pressure.");
}


void config_verify_hydro()
{
  cart_assert(gas_temperature_floor >= 0.0); 

  cart_assert(gas_density_floor > 0.0);

  cart_assert(pressureless_fluid_eos==0 || pressureless_fluid_eos==1);

  cart_assert(apply_lapidus_viscosity==0 || apply_lapidus_viscosity==1);

  cart_assert(smooth_density_gradients==0 || smooth_density_gradients==1);

  cart_assert(pressure_floor_min_level>=-1 && pressure_floor_min_level<=max_level);

  cart_assert(pressure_floor_factor > 0.0);

  cart_assert(!(radiation_pressure_factor < 0.0));

#ifdef BLASTWAVE_FEEDBACK

  cart_assert(blastwave_time_floor > 0.0);

  cart_assert(blastwave_time_cut > 0.0 && blastwave_time_cut > blastwave_time_floor);

#endif /* BLASTWAVE_FEEDBACK */
}


float cell_gas_kinetic_energy(int cell) {
	int j;
	double ke = 0.0;

	if(cell_gas_density(cell) > 0.0) {
		for(j=0; j<nDim; j++) ke += (double)cell_momentum(cell,j)*(double)cell_momentum(cell,j);
		return (float)(0.5*ke/cell_gas_density(cell));
	}
	else return 0.0;
}

#ifdef RADIATIVE_TRANSFER
float rtTem(int cell);
#endif /* RADIATIVE_TRANSFER */

float cell_gas_temperature(int cell) {
	if(cell_gas_density(cell) > 0.0) {
#ifdef RADIATIVE_TRANSFER
		return rtTem(cell);
#else
		return (cell_gas_gamma(cell)-1)*constants->wmu*cell_gas_internal_energy(cell)/cell_gas_density(cell);
#endif
	} else return 0.0;
}


/*
//  Sobolev approximation factors
*/
float cell_sobolev_length2(int cell, int level, float *vel)
{
  int i, j, nb[num_neighbors];
  float s, d, len;
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
  len = 0.5*cell_size[level]*cell_gas_density(cell)/(1.0e-30+sqrt(s));
  if(len > num_grid) len = num_grid;

  /*
  //  Velocity factor
  */
  if(vel!=NULL && cell_gas_density(cell)>0.0)
    {
      for(i=0; i<nDim; i++)
	{
	  vc[i] = cell_var(cell,HVAR_MOMENTUM+i)/cell_gas_density(cell);
	}

      s = 0.0;
      for(j=0; j<num_neighbors; j++) if(cell_gas_density(nb[j]) > 0.0)
	{
	  for(i=0; i<nDim; i++)
	    {
	      d = cell_var(nb[j],HVAR_MOMENTUM+i)/cell_gas_density(nb[j]) - vc[i];
	      s += d*d;
	    }
	}
      *vel = sqrt(s/num_neighbors);
    }

  return len;
}


#ifndef RADIATIVE_TRANSFER
float cell_radiation_pressure(int cell)
{
  cart_error("Radiation pressure without RADIATIVE_TRANSFER is not implemented yet.");
  return 0.0;
}
#endif /* RADIATIVE_TRANSFER */


float cell_gas_sound_speed( int icell ) {
	return sqrt(cell_gas_gamma(icell)*(cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell)/cell_gas_density(icell));
}

#endif /*HYDRO*/
