#include "config.h"
#if defined(STAR_FORMATION) && defined(AGN)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "io.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "particle_buffer.h"
#include "starformation.h"
#include "starformation_feedback.h"
#include "times.h"
#include "tree.h"
#include "units.h"
#include "agn.h"

/* agn model parameters */
int agn_accretion_recipe = 1;              /* Set to 0 for pure eddington accretion, 1 for Bondi accretion (constant-alpha), 2 for Bondi accretion (constant-beta) */
int agn_feedback_recipe = 0;
int agn_feedback_storage = 1;              /* Set to 0 for no storage, 1 for storage via Booth and Schaye (2009) Eqn (7) */
int agn_dv_on_off = 1;                     /* Set dv to 0 in denominator of accretion rate or not (for Bondi-type accretion) */

double eddington_factor = 1.0;             /* normalization of Eddington accretion (1 for standard value) */
double radiative_efficiency = 0.1;         /* eta or eps_r */
double feedback_efficiency = 0.05;         /* eps_f */
double bondi_normalization = 100.;         /* alpha normalization of standard Bondi implementations */
double bondi_pivot_gas_density = 0.1;      /* gas density below which alpha=1 model applies, usually set to sf_min_gas_number_density */
double bondi_exponent = 2.0;               /* Beta for constant-Beta (density dependent) model of Bondi implementations */
double minimum_agn_feedback_temperature = 1e8; /* Tmin for critical energy in feedback */  
double sink_particle_delta = 4.0;          /* This determines maximum r_acc=sink_particle_delta*resolution_size */
double agn_merge_radius = 4.0;             /* maximum distance to merge two black hole particles */
double agn_merge_velocity = 1000.;         /* maximum velocity to merge in km/s */
double agn_seed_mass = 1e5;                /* initial seed black hole mass [Msun] */
double agn_halo_mass_threshold = 2e11;     /* minimum halo mass to seed into [Msun] */

void config_init_agn() {

	control_parameter_add(control_parameter_double,&eddington_factor,
			"agn:eddington_factor",
			"Normalization of Eddington accretion rate (typically a small number 1-few)");
	control_parameter_add2(control_parameter_double,&bondi_normalization,
			"agn:bondi_normalization","agn:alpha",
			"Normalization of Bondi-Hoyle accretion rate, used to compensate for not resolving r_bondi" );
	control_parameter_add3(control_parameter_double,&radiative_efficiency,
			"agn:radiative_efficiency","agn:eps_r","agn:eta",
			"Radiative efficiency, L = \\eta \\dot{M} c^2" );
	control_parameter_add2(control_parameter_double,&feedback_efficiency,
			"agn:feedback_efficiency","agn:eps_f",
			"Feedback efficiency, \\dot{E} = eps_f L = \\eps_f \\eps_r \\dot{M} c^2" );
	control_parameter_add(control_parameter_int, &agn_accretion_recipe,
			      "agn:accretion_recipe", "Accretion recipe for AGN: 0 = pure Eddington, 1 = Bondi Constant alpha model in single cell, 2= Bondi Constant beta model in single cell" ); 
	control_parameter_add(control_parameter_double,&bondi_pivot_gas_density,"agn:bondi_pivot_gas_density","Critical gas density above which the beta bondi model applies, \\alpha=(n_H/n_H^*)^\\beta");  
	control_parameter_add2(control_parameter_double,&bondi_exponent,"agn:bondi_exponent","agn:beta","exponential dependence of normalization of Bondi-Hoyle accretion rate on ISM gas density"); 
	control_parameter_add2(control_parameter_double,&minimum_agn_feedback_temperature,"agn:minimum_feedback_temperature","agn:Tmin","Minimum temperature that defines critical energy for an agn heating event to be triggered, E_crit=m_g*k*Temperature/(mu*Mh*(gamma-1))"); 
	control_parameter_add(control_parameter_int, &agn_feedback_storage, "agn:feedback_storage", "Feedback storage for AGN: 0 = no storage, 1 = storage with E_crit = m_g*k*Tmin/(gamma-1)/mu/mH"); 
	control_parameter_add2(control_parameter_double, &sink_particle_delta, "agn:sink_particle_delta", "agn:delta","Determines accretion radius for sink particle which is proportional to the resolution, r_{acc}=\\delta \\Delta x");
	control_parameter_add2(control_parameter_double, &agn_merge_radius, "agn:merge_radius", "agn:merge_radius","Determines maximum relative distance for two merging black holes with default, multiples of the resolution size"); 
	control_parameter_add2(control_parameter_double, &agn_merge_velocity, "agn:merge_velocity", "agn:merge_velocity","Determines maximum relative velocity for two merging black holes"); 
	control_parameter_add(control_parameter_int, &agn_dv_on_off, "agn:dv_on_off", "dv on or off for AGN: 0 = no dv in calculation, 1 = dv in calculation");
	control_parameter_add2(control_parameter_double,&agn_seed_mass,
            "agn:seed_mass", "agn_seed_mass", "The initial seed mass for black hole particles [Msun]");
    control_parameter_add2(control_parameter_double,&agn_halo_mass_threshold,
            "agn:seed_halo_min_mass", "agn_halo_mass_threshold", "The minimum halo mass to seed black holes into (strongly resolution dependent) [Msun]");
}

void config_verify_agn() {
	VERIFY(agn:eddington_factor, eddington_factor >= 1.0 );
	VERIFY(agn:bondi_normalization, bondi_normalization >= 1.0 );
	VERIFY(agn:radiative_efficiency, radiative_efficiency > 0.0 && radiative_efficiency <= 1.0 );
	VERIFY(agn:feedback_efficiency, feedback_efficiency >= 0.0 && feedback_efficiency <= 1.0 );
	VERIFY(agn:accretion_recipe, agn_accretion_recipe == 0 || agn_accretion_recipe == 1 || agn_accretion_recipe == 2 ); 
	VERIFY(agn:bondi_pivot_gas_density, bondi_pivot_gas_density > 0.0 ); 
	VERIFY(agn:bondi_exponent, bondi_exponent > 0.0 ); 
	VERIFY(agn:minimum_feedback_temperature, minimum_agn_feedback_temperature >= 1e6 ); 
	VERIFY(agn:feedback_storage, agn_feedback_storage == 0 || agn_feedback_storage == 1 ); 
	VERIFY(agn:sink_particle_delta, sink_particle_delta > 0.0); 
	VERIFY(agn:merge_radius, agn_merge_radius > 0.0 && agn_merge_radius < pow( 2.0, max_level - 1 ) );
	VERIFY(agn:merge_velocity, agn_merge_velocity > 0.0 );
	VERIFY(agn:dv_on_off, agn_dv_on_off == 0 || agn_dv_on_off == 1 );

	VERIFY(agn:seed_mass, agn_seed_mass > 0.0 );
	VERIFY(agn:seed_halo_min_mass, agn_halo_mass_threshold > 0.0 );
}

double sink_sample_weights[SINK_SAMPLES_3D];
double sink_sample_offset[SINK_SAMPLES_3D][nDim];
int num_sink_samples;

void init_agn() {
	int i, j, k;
	double x, y, z, r2;
	double weight_sum;

	/* compute sink formalism sampled points' weights and offsets */
	num_sink_samples = 0;
	for ( i = -SINK_RADIUS_SAMPLES; i <= SINK_RADIUS_SAMPLES; i++ ) {
		x = SINK_ACCRETION_RADIUS*(double)i;
		for ( j = -SINK_RADIUS_SAMPLES; j <= SINK_RADIUS_SAMPLES; j++ ) {
			y = SINK_ACCRETION_RADIUS*(double)j;
			for ( k = -SINK_RADIUS_SAMPLES; k <= SINK_RADIUS_SAMPLES; k++ ) {
				z = SINK_ACCRETION_RADIUS*(double)k;

				r2 = x*x+y*y+z*z;
				if ( r2 <= SINK_ACCRETION_RADIUS*SINK_ACCRETION_RADIUS ) {
					sink_sample_offset[num_sink_samples][0] = x;
					sink_sample_offset[num_sink_samples][1] = y;
					sink_sample_offset[num_sink_samples][2] = z;

					sink_sample_weights[num_sink_samples] = exp( -r2 );
					num_sink_samples++;
				}
			}
		}
	}

	weight_sum = 0.0;
	for ( i = 0; i < num_sink_samples; i++ ) {
		weight_sum += sink_sample_weights[i];
	}

	for ( i = 0; i < num_sink_samples; i++ ) {
		sink_sample_weights[i] /= weight_sum;
	}

}

double compute_r_K( int ipart, int icell ) {  
	double r_K;
	double dv_cell;
	double r_BH; /* Bondi-Hoyle radius */
	double Delta_x=cell_size[max_level];
	double Mbh = particle_mass[ipart];

	cart_assert( ipart >= 0 && ipart < num_star_particles );
	cart_assert( icell >= 0 && icell < num_cells );

	dv_cell =    (particle_v[ipart][0]-cell_momentum(icell,0)/cell_gas_density(icell)) *
		(particle_v[ipart][0]-cell_momentum(icell,0)/cell_gas_density(icell)) +
		(particle_v[ipart][1]-cell_momentum(icell,1)/cell_gas_density(icell)) *
		(particle_v[ipart][1]-cell_momentum(icell,1)/cell_gas_density(icell)) +
		(particle_v[ipart][2]-cell_momentum(icell,2)/cell_gas_density(icell)) *
		(particle_v[ipart][2]-cell_momentum(icell,2)/cell_gas_density(icell));   

	r_BH = constants->G*particle_mass[ipart]*units->mass/( cell_gas_sound_speed(icell)*cell_gas_sound_speed(icell) + dv_cell )/units->velocity/units->velocity/units->length;

	if ( r_BH < Delta_x/4 ) {
		r_K = Delta_x/4;
	} else if (r_BH <= sink_particle_delta * Delta_x/2) {
		r_K = r_BH;
	} else {
		r_K = sink_particle_delta * Delta_x/2;
	}

	return r_K;                                                                                                     
}

int sort_particles_by_mass( const void *a, const void *b ) {
	int index_a = *(int *)a;
	int index_b = *(int *)b;

	cart_assert( index_a >= 0 && index_a < num_particles );
	cart_assert( index_b >= 0 && index_b < num_particles );

	if ( particle_mass[index_a] < particle_mass[index_b] ) {
		return -1;
	} else if ( particle_mass[index_a] > particle_mass[index_b] ) {
		return 1;
	} else {
		/* decreasing order of id */
		return ( particle_id[index_b] - particle_id[index_a] );
	}
} 

#endif /* STAR_FORMATION && AGN */
