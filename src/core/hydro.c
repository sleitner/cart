#include "config.h"
#ifdef HYDRO 

#include <math.h>
#include <stdio.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "hydro.h"
#include "iterators.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

DEFINE_LEVEL_ARRAY(int,level_sweep_dir);


int pressure_floor_min_level = max_level; /* NG: that used to be MinL_Jeans define */
float pressure_floor_factor = 10.0;

int pressureless_fluid_eos = 0;           /* NG: that used to be PRESSURELESS_FLUID define */
int apply_lapidus_viscosity = 1;          /* NG: that used to be LAPIDUS define */
int smooth_density_gradients = 1;         /* NG: that used to be DENSGRADSMOOTH define */

float gas_density_floor = 1e-6;
float gas_temperature_floor = 3.0;        /* NG: that used to be T_min define */

#ifdef BLASTWAVE_FEEDBACK
double blastwave_time_floor = 1.0e-30; 
double blastwave_time_cut = 1.0e-20;
#endif /* BLASTWAVE_FEEDBACK */

#ifdef ISOTROPIC_TURBULENCE_ENERGY
double fix_isotropic_turbulence_dissipation_time = -1;
double fraction_SN_to_isotropic_turbulence  = 0.0;
#endif /* ISOTROPIC_TURBULENCE_ENERGY */

#if num_extra_energy_variables > 0
float extra_energy_gammas[num_extra_energy_variables];
#else
float extra_energy_gammas[1]; /* avoids extra flags around pragmas*/
#endif /* num_extra_energy_variables > 0 */

#if defined(ISOTROPIC_TURBULENCE_ENERGY) 
struct gamma_t
{
    float isotropic_turbulence;
};
extern struct gamma_t gamma_internal = {5.0/3.0};
#endif /* ISOTROPIC_TURBULENCE_ENERGY */


void config_init_hydro()
{
  int level;

  control_parameter_add3(control_parameter_float,&gas_temperature_floor,"gas-temperature-floor","gas_temperature_floor","T_min","the minimum temperature for the gas (in K).");

  control_parameter_add2(control_parameter_float,&gas_density_floor,"gas-density-floor","gas_density_floor","the minimum densitye for the gas (in code units).");

  control_parameter_add2(control_parameter_bool,&pressureless_fluid_eos,"pressureless-fluid-eos","pressureless_fluid_eos","use a pressureless fluid equation of state.");

  control_parameter_add2(control_parameter_bool,&apply_lapidus_viscosity,"apply-lapidus-viscosity","apply_lapidus_viscosity","apply Lapidus viscosity in the hydro flux calculations (boolean value).");

  control_parameter_add2(control_parameter_bool,&smooth_density_gradients,"smooth-density-gradients","smooth_density_gradients","smooth overly steep density gradients in the hydro flux calculations (boolean value).");

  control_parameter_add3(control_parameter_int,&pressure_floor_min_level,"pressure-floor-min-level","pressure_floor_min_level","MinL_Jeans","the level to apply the pressure floor. If this value is set to -1, the pressure floor correction is disabled.");

  control_parameter_add(control_parameter_float,&pressure_floor_factor,"@pressure-floor-factor","the factor to scale the pressure floor with. The default, thoroughly tested value is 10. If you change it, make sure you know what you are doing.");

#ifdef ISOTROPIC_TURBULENCE_ENERGY
  control_parameter_add2(control_parameter_float,&gamma_internal.isotropic_turbulence,"gamma:isotropic_turbulence","gamma.isotropic_turbulence","gamma for turbulent energy");
  
  control_parameter_add2(control_parameter_double,&fraction_SN_to_isotropic_turbulence,"isotropic_turbulence:fraction-SN-to-isotropic_turbulence","fraction_SN_to_isotropic_turbulence","fraction of SN energy that goes to isotropic_turbulence.");

  control_parameter_add3(control_parameter_time,&fix_isotropic_turbulence_dissipation_time,"isotropic_turbulence:fix-dissipation-time","fix_isotropic_turbulence_dissipation_time","isotropic_turbulence.fix_dissipation_time", "time scale over which isotropic_turbulence dissipates, if not set then defaults to cell crossing time.");
#endif /* ISOTROPIC_TURBULENCE_ENERGY */
  
  for(level=min_level; level<=max_level; level++)
    {
      level_sweep_dir[level] = 0;
    }

#ifdef ISOTROPIC_TURBULENCE_ENERGY
  isotropic_turbulence_gamma = gamma_internal.isotropic_turbulence;
#endif /* ISOTROPIC_TURBULENCE_ENERGY */

}


void config_verify_hydro()
{
  VERIFY(gas-temperature-floor, gas_temperature_floor >= 0.0 ); 

  VERIFY(gas-density-floor, gas_density_floor > 0.0 );

  VERIFY(pressure-floor-min-level, pressure_floor_min_level>=-1 && pressure_floor_min_level<=max_level );

  VERIFY(@pressure-floor-factor, pressure_floor_factor > 0.0 );

#ifdef BLASTWAVE_FEEDBACK

  cart_assert(blastwave_time_floor > 0.0);

  cart_assert(blastwave_time_cut > 0.0 && blastwave_time_cut > blastwave_time_floor);

#endif /* BLASTWAVE_FEEDBACK */

#ifdef ISOTROPIC_TURBULENCE_ENERGY
  VERIFY(fix-isotropic_turbulence-dissipation-time, fix_isotropic_turbulence_dissipation_time > 0.0 || fix_isotropic_turbulence_dissipation_time == -1.0 );

  VERIFY(isotropic_turbulence:fraction-SN-to-isotropic_turbulence, fraction_SN_to_isotropic_turbulence >= 0.0 && fraction_SN_to_isotropic_turbulence <= 1.0);
  
  VERIFY(gamma:isotropic_turbulence , gamma_internal.isotropic_turbulence  > 1.0);
#endif /* ISOTROPIC_TURBULENCE_ENERGY */

}


void hydro_magic_one_cell( int icell ) {
        int j;
	float average_density;
	int neighbors[num_neighbors];
	static int failed = 0;

	double kinetic_energy;
	double thermal_energy;

	/* do density floor stuff */
	if ( cell_gas_density(icell) < gas_density_floor ) {
		average_density = 0.0;
		cell_all_neighbors( icell, neighbors );

		for ( j = 0; j < num_neighbors; j++ ) {
			average_density += cell_gas_density( neighbors[j] );
		}

		if(!failed) {
			cart_debug("---------------------------------------------------------");
			cart_debug("HIT DENSITY FLOOR:");
			cart_debug("old density = %e g/cc", cell_gas_density(icell)*units->density/constants->gpercc );
			cart_debug("new density = %e g/cc", MAX( average_density/(float)num_neighbors, gas_density_floor ) );
			cart_debug("T  = %e K", cell_gas_temperature(icell)*units->temperature/constants->K );
			cart_debug("P  = %e ergs cm^-3", cell_gas_pressure(icell)*units->energy_density/constants->barye );
			cart_debug("v  = %e %e %e cm/s",
					cell_momentum(icell,0)/cell_gas_density(icell)*units->velocity/constants->cms,
					cell_momentum(icell,1)/cell_gas_density(icell)*units->velocity/constants->cms,
					cell_momentum(icell,2)/cell_gas_density(icell)*units->velocity/constants->cms );
			cart_debug("---------------------------------------------------------");
		}	

		failed = 1;

		cell_gas_density(icell) = MAX( average_density/(float)num_neighbors, gas_density_floor );
	}

	kinetic_energy = cell_gas_kinetic_energy(icell);
	thermal_energy = gas_temperature_floor/(units->temperature*constants->wmu*(constants->gamma-1)) * cell_gas_density(icell);

	cell_gas_internal_energy(icell) = MAX( cell_gas_internal_energy(icell), thermal_energy );
	cell_gas_energy(icell) = MAX( cell_gas_energy(icell), thermal_energy+kinetic_energy );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	cell_electron_internal_energy(icell) = MAX( cell_electron_internal_energy(icell), thermal_energy*constants->wmu/constants->wmu_e );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
	for ( j = 0; j < num_extra_energy_variables; j++ ) {
            cell_extra_energy_variables(icell,j) = MAX( 0, cell_extra_energy_variables(icell,j) );
        }

	for ( j = 0; j < num_chem_species; j++ ) {
		/* 
		   1e-15 may be too large a number for ionic species;
		   at least let's scale them with density and make 1e-20;
Gnedin: 1e-20 is not small enough for chemistry, making it 1e-30
		 */
		cell_advected_variable(icell,j) = MAX( 1e-30*cell_gas_density(icell), cell_advected_variable(icell,j) );
		/* Doug had it like that:
		   cell_advected_variable(icell,j) = MAX( 1e-15, cell_advected_variable(icell,j) );
		 */
	}
}

void hydro_magic( int level ) {
    int i;
    int num_level_cells;
    int *level_cells;

    start_time( WORK_TIMER );

    select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i), shared(num_level_cells,level_cells,cell_child_oct)
    for ( i = 0; i < num_level_cells; i++ ) {
		hydro_magic_one_cell(level_cells[i]);
	}
    cart_free( level_cells );

    end_time( WORK_TIMER );
} 

void hydro_split_update( int level ) {
	int i, j, k;
	int icell;
	int num_level_cells;
	int *level_cells;
	int children[num_children];
	double new_var;
	const double factor = ((double)(1.0/(1<<nDim)));

	if ( level < max_level ) {
		start_time( WORK_TIMER );

		select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_REFINED, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,j,k,children,new_var), shared(num_level_cells,level_cells,cell_child_oct,cell_vars)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			/* average over children */
			cell_all_children( icell, children );

			for ( j = 0; j < num_hydro_vars; j++ ) {
				new_var = 0.0;
				for ( k = 0; k < num_children; k++ ) {
					new_var += cell_hydro_variable(children[k], j);
				}

				cell_hydro_variable(icell, j) = new_var*factor; 
			}
		}
		cart_free( level_cells );

		end_time( WORK_TIMER );
	}
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

#ifdef ISOTROPIC_TURBULENCE_ENERGY
float cell_isotropic_turbulence_temperature(int cell){
	if(cell_gas_density(cell) > 0.0) {
	        return (isotropic_turbulence_gamma-1)*constants->wmu*cell_isotropic_turbulence_energy(cell)/cell_gas_density(cell);
	} else return 0.0;
}
#endif /* ISOTROPIC_TURBULENCE_ENERGY */

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


float cell_gas_sound_speed( int icell ) {
        float gPeff;
        int j;
        gPeff = cell_gas_gamma(icell)*(cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell);
        for ( j = 0; j < num_extra_energy_variables; j++ ) {
                gPeff += extra_energy_gamma(j)*cell_extra_energy_pressure(icell,j);
        }
        return sqrt(gPeff/cell_gas_density(icell));
}

#endif /*HYDRO*/
