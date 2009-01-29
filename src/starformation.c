#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "tree.h"
#include "particle.h"
#include "starformation.h"
#include "timestep.h"
#include "units.h"
#include "constants.h"
#include "iterators.h"
#include "auxiliary.h"

#ifdef STARFORM

/*
//  NG: fall back to the old-style recipe if this is not set
*/
#ifndef SF_RECIPE
#define SF_RECIPE 1
#endif

int num_local_star_particles;
int last_star_id;
int num_new_stars;

double total_stellar_mass = 0.0;
double total_stellar_initial_mass = 0.0;

float star_tbirth[num_star_particles];
float star_initial_mass[num_star_particles];

#ifdef ENRICH
float star_metallicity_II[num_star_particles];
#ifdef ENRICH_SNIa
float star_metallicity_Ia[num_star_particles];
#endif /* ENRICH_SNIa */
#endif /* ENRICH */

float star_formation_volume_min[nDim];
float star_formation_volume_max[nDim];

int star_formation_frequency[max_level-min_level+1];

/* star formation parameters */
double alpha_SF	= 1.5;
double eps_SF = 1.5;
double dtmin_SF = 1e6;
double tau_SF = 1e7;
double dm_star_min = 0.0;
double rho_SF = 1e-1;
double T_SF = 2e4;
double a_IMF = 0.0;
double aM_stl = 0.1;
double aM_stu = 100.0;
double aM_SNII = 8.0;
double aM_SNIa1 = 3.0;
double aM_SNIa2 = 8.0;
double ejM_SNIa = 1.3;
double C_SNIa = 1.5e-2;
double t_SNIa = 0.2;
double E_51 = 2.0;
double t_fb = 1000.0;
double T0_ml = 5.0;
double c0_ml = 0.05;

double eps_SFH2 = 0.01;
double fH2_SFH2 = 0.1;
double den_SFH2_eff = 100.0;

double C_SFR;
double C_fb;
double C_fbIa;
double fmass_met;
double rho_SF_fact, rho_SFH2_eff, den_SFH2_fact;
double aMSN_ave;
double fmass_SN;
double RIaf;

double f_IMF( double amstar ) {
	/* Miller-Scalo (1979, ApJS 41, 513, eq. 30, Table 7) IMF */
	#define C_1     1.09
	#define C_2     -1.02
	return exp( -C_1 * ( (log10(amstar) - C_2)*(log10(amstar) - C_2) ) ) / amstar;
                                                                                                                                                            
	/* Chabrier, G. (2001, ApJ 554, 1274) */
	/* #define am0_Ch 716.4 */
	/* #define beta_Ch 0.25 */
	/* #define alpha_Ch (-3.3) */
	/* return exp( -pow( am0_Ch/amstar), beta_Ch) * pow(amstar,alpha_Ch); */
}

double f_SNIa_func( double xd ) {
	return max( exp( -xd*xd ) / sqrt( xd ), 1e-20 );
}
                                                                                                                                                            
double f_SNIa( double xd ) {
	return max( exp(-xd*xd) * sqrt(xd*xd*xd), 1e-20 );
}
                                                                                                                                                            
double f_IMF_plaw( double amstar ) {
	return pow( amstar, -a_IMF );
}
                                                                                                                                                            
double fm_IMF_plaw( double amstar ) {
	return pow( amstar, 1.0 - a_IMF );
}
                                                                                                                                                            
double fm_IMF( double amstar ) {
	return amstar * f_IMF(amstar);
}

double fmet_ej( double amstar ) {
        return min( 0.2, max( 0.01*amstar - 0.06, 1e-20 ) );
}

double fej_IMF( double amstar ) {
	return amstar * f_IMF(amstar) * fmet_ej(amstar);
}
                                                                                                                                                            
double fej_IMF_plaw( double amstar ) {
	return pow( amstar, 1.0 - a_IMF ) * fmet_ej( amstar );
}
                                                                                                                                                            
void init_star_formation() {
        int j;
	double Aprime;
	double aN_SNII, An_IMF;

#if (SF_RECIPE == 1)
	C_SFR = eps_SF * 2.5*pow(10.0, 6.0 - 16.0*alpha_SF) * t0 * pow(rho0,alpha_SF-1.0);
#else
	C_SFR = 49.53*2.5e-18*t0*sqrt(rho0);  /* This is from Kostas */
#endif
	/*
	// Code units
	*/  
	dm_star_min /= aM0; 

	if ( a_IMF > 0.0 ) {
		if ( a_IMF == 2.0 ) {
			Aprime = 1.0 / log(aM_stu/aM_stl);
		} else {
			Aprime = (2.0-a_IMF) / ( pow(aM_stu,2.0-a_IMF) - pow(aM_stl,2.0-a_IMF) );
		}

		An_IMF = Aprime;

		if ( a_IMF == 1.0 ) {
			Aprime *= log(aM_stu/aM_SNII);
		} else {
			Aprime /= (1.0-a_IMF)*( pow(aM_stu,2.0-a_IMF) - pow( aM_SNII,1.0-a_IMF) );
		}

#ifdef ENRICH
		fmass_met = integrate( fej_IMF_plaw, aM_SNII, aM_stu, 1e-6, 1e-9 );

		if ( a_IMF == 2.0 ) {
			fmass_met /= log(aM_stu/aM_stl);
		} else {
			fmass_met *= (2.0-a_IMF)/( pow(aM_stu,2.0-a_IMF) - pow(aM_stl,2.0-a_IMF) );
		}
#endif /* ENRICH */

		aN_SNII = integrate( f_IMF_plaw, aM_SNII, aM_stu, 1e-6, 1e-9 );
		aMSN_ave = integrate( fm_IMF_plaw, aM_SNII, aM_stu, 1e-6, 1e-9 ) / aN_SNII;
		aN_SNII *= An_IMF;

		fmass_SN = integrate( fm_IMF_plaw, aM_SNII, aM_stu, 1e-6, 1e-9 ) /
				integrate( fm_IMF_plaw, aM_stl, aM_stu, 1e-6, 1e-9 );

#ifdef FEEDBACK_SNIa
		RIaf = 1e-9 * t0 * aM0 * C_SNIa * integrate( f_IMF_plaw, aM_SNIa1, aM_SNIa2, 1e-6, 1e-9 ) /
			integrate( fm_IMF_plaw, aM_SNIa1, aM_SNIa2, 1e-6, 1e-9 ) /
			integrate( f_SNIa_func, 1e-2, 1e3, 1e-6, 1e-9 ) / t_SNIa;
#endif /* FEEDBACK_SNIa */

	} else {
		/* Miller-Scalo IMF */
		An_IMF = 1.0 / integrate( fm_IMF, aM_stl, aM_stu, 1e-6, 1e-9 );

		if ( An_IMF <= 0.0 ) {
			cart_error("An_IMF <= 0.0!");
		}

		Aprime = An_IMF * integrate( f_IMF, aM_SNII, aM_stu, 1e-6, 1e-9 );

		if ( Aprime <= 0.0 ) {
			cart_error("Aprime <= 0!");
		}

#ifdef ENRICH 
		fmass_met = integrate( fej_IMF, aM_SNII, aM_stu, 1e-6, 1e-9 ) /
			integrate( fm_IMF, aM_stl, aM_stu, 1e-6, 1e-9 );
#endif /* ENRICH */

		aN_SNII = integrate( f_IMF, aM_SNII, aM_stu, 1e-6, 1e-9 );
		aMSN_ave = integrate( fm_IMF, aM_SNII, aM_stu, 1e-6, 1e-9 ) / aN_SNII;
		aN_SNII *= An_IMF;

		fmass_SN = integrate( fm_IMF, aM_SNII, aM_stu, 1e-6, 1e-9 ) /
				integrate( fm_IMF, aM_stl, aM_stu, 1e-6, 1e-9 );

#ifdef FEEDBACK_SNIa
		RIaf = 1e-9 * t0 * aM0 * C_SNIa * integrate( f_IMF, aM_SNIa1, aM_SNIa2, 1e-6, 1e-9 ) /
			integrate( fm_IMF, aM_SNIa1, aM_SNIa2, 1e-6, 1e-9 ) /
			integrate( f_SNIa_func, 1e-3, 1e3, 1e-6, 1e-9 ) / t_SNIa;
#endif /* FEEDBACK_SNIa */
	}

	C_fb = 1e51 * E_51 * Aprime * aM0 / E0;
	C_fbIa = 1e51 * E_51 / E0;
	rho_SF_fact = rho_SF * 8.9e4 / (hubble*hubble) / Omega0 / ( 1.0 - Y_p );

	/*
	//  NG: This limit is very obscure, disable by default
	*/
	for(j=0; j<nDim; j++)
	  {
	    star_formation_volume_min[j] = 0.0;
	    star_formation_volume_max[j] = num_grid;
	  }

	if ( local_proc_id == MASTER_NODE ) {
		cart_debug("tau_SF = %e", tau_SF );
		cart_debug("T_SF = %e", T_SF );
		cart_debug("C_fb = %e", C_fb );
		cart_debug("C_sfr = %e", C_SFR );
		cart_debug("rho_SF_fact = %e", rho_SF_fact );
	}
}

#ifdef HYDRO
void star_formation( int level, int time_multiplier ) {
	int i, j;
	int icell;
	float pos[nDim];
	int num_level_cells;
	int *level_cells;
	int do_star_formation;
	double cell_fraction;
	double fmass;
	double T_SF_max;
	double Tcell;
	double dm_star;
	double rho_SF_min;
	double tau_SF_code;
	double tau_SF_eff;
	double dt_eff;
	double P_SF, fH2_cell, zSol_cell;

	cell_fraction = 0.667 * cell_volume[level];

	T_SF_max = T_SF / T0 * aexp[level]*aexp[level];
	rho_SF_min = max( rho_SF_fact * aexp[level]*aexp[level]*aexp[level], 200.0 * Omegab0 / Omega0 );
	tau_SF_code = tau_SF / ( t0 * aexp[level]*aexp[level] );
	dt_eff = dtl[level] * time_multiplier;

#if (SF_RECIPE == 1)
	fmass = pow(aexp[level], (5.0 - 3.0*alpha_SF)) * C_SFR * cell_volume[level] * tau_SF_code;
#else
	den_SFH2_fact = (1.123e-5*Omega0*hubble*hubble)/(aexp[level]*aexp[level]*aexp[level]);
	rho_SFH2_eff = den_SFH2_eff/den_SFH2_fact;
	fmass = eps_SFH2 * sqrt(aexp[level]) * C_SFR * cell_volume[level] * tau_SF_code;
#endif

	/* probability of forming a star is Poisson with <t> = tau_SF */
	P_SF = exp( -dt_eff / tau_SF_code );

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_leaf(icell) ) {
			/* check position */
			cell_position( icell, pos );

			do_star_formation = 1;
			for ( j = 0; j < nDim; j++ ) {
				if ( pos[j] < star_formation_volume_min[j] || pos[j] > star_formation_volume_max[j] ) {
					do_star_formation = 0;
				}
			}

			if ( do_star_formation ) {
				if ( cell_gas_density(icell) > rho_SF_min ) {
					Tcell = cell_gas_pressure(icell) / cell_gas_density(icell);
#ifdef RADIATIVE_TRANSFER
					fH2_cell = 2*cell_H2_density(icell)/(2*cell_H2_density(icell)+cell_HI_density(icell));
#else /* RADIATIVE_TRANSFER */
#ifdef ENRICH
					zSol_cell = cell_gas_metallicity(icell);
#else
					zSol_cell = 0.0;
#endif /* ENRICH */
					fH2_cell = (max(1.0e-3,zSol_cell)*den_SFH2_fact*cell_gas_density(icell) > 30.0) ? 1.0 : 0.0;
#endif /* RADIATIVE_TRANSFER */
					if ( Tcell < T_SF_max
#if (SF_RECIPE != 1)
					     && fH2_cell > fH2_SFH2
#endif
 ) {
						/* randomly generate particle on timescale tau_SF_eff */
						if ( cart_rand() > P_SF ) {
							/* how big of a star particle to create? */
#if (SF_RECIPE == 1)
						        dm_star = fmass * pow( cell_gas_density(icell), alpha_SF );
#else /* SF_RECIPE == 1 */
							dm_star = 0.0;  /* compiler will remove this if the setting is correct */
#if (SF_RECIPE == 2)
							dm_star = fmass*fH2_cell*cell_gas_density(icell)*sqrt(rho_SFH2_eff);
#endif
#if (SF_RECIPE == 3)
							dm_star = fmass*fH2_cell*cell_gas_density(icell)*sqrt(max(cell_gas_density(icell),rho_SFH2_eff));
#endif
#if (SF_RECIPE == 4)
							dm_star = fmass*fH2_cell*cell_gas_density(icell)*sqrt(cell_gas_density(icell));
#endif

#endif /* SF_RECIPE == 1 */

							dm_star = min( max(dm_star,dm_star_min), cell_fraction * cell_gas_density(icell) );

							/* create the new star */
							create_star_particle( icell, dm_star );
						}
					}
				}
			}
		}
	}

	cart_free( level_cells );
}

void create_star_particle( int icell, float mass ) {
	int i;
	int ipart;
	int id;
	int level;
	float pos[nDim];
	float new_density;
	float density_fraction;

	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( mass > 0.0 );

	id = last_star_id + local_proc_id + 1;
	last_star_id += num_procs;
	num_new_stars++;

	ipart = particle_alloc( id );
	cart_assert( ipart < num_star_particles );

	/* place particle at center of cell with cell momentum */
	cell_position(icell, pos );
	level = cell_level(icell);

	for ( i = 0; i < nDim; i++ ) {
		particle_x[ipart][i] = pos[i];
	}

	for ( i = 0; i < nDim; i++ ) {
		particle_v[ipart][i] = cell_momentum(icell,i) / cell_gas_density(icell);
	}

	particle_t[ipart] = tl[level];
	particle_dt[ipart] = dtl[level];

	star_tbirth[ipart] = tl[level];
	particle_mass[ipart] = mass;
	star_initial_mass[ipart] = mass;

#ifdef ENRICH
	star_metallicity_II[ipart] = cell_gas_metallicity_II(icell) / cell_gas_density(icell);
#ifdef ENRICH_SNIa
	star_metallicity_Ia[ipart] = cell_gas_metallicity_Ia(icell) / cell_gas_density(icell);
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
	
	/* insert particle into cell linked list */
	insert_particle( icell, ipart );

	/* adjust cell values */
	new_density = cell_gas_density(icell) - mass * cell_volume_inverse[level];
	density_fraction = new_density / cell_gas_density(icell);

	cell_gas_density(icell) = new_density;
	cell_gas_energy(icell) *= density_fraction;
	cell_gas_internal_energy(icell) *= density_fraction;
	cell_gas_pressure(icell) *= density_fraction;
	cell_momentum(icell,0) *= density_fraction;
	cell_momentum(icell,1) *= density_fraction;
	cell_momentum(icell,2) *= density_fraction;
		
#ifdef ENRICH
	cell_gas_metallicity_II(icell) = max( 1e-17,
		cell_gas_metallicity_II(icell) - star_metallicity_II[ipart]*mass*cell_volume_inverse[level] );
#ifdef ENRICH_SNIa
	cell_gas_metallicity_Ia(icell) = max( 1e-17,
		cell_gas_metallicity_Ia(icell) - star_metallicity_Ia[ipart]*mass*cell_volume_inverse[level] );
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
}

#endif /* HYDRO */

void remap_star_ids() {
	int i;
	int proc;
	int ipart;
	int block;
	int max_stars;	
	int new_id;
	int total_new_stars;
	int *block_ids;
	int proc_new_stars[MAX_PROCS];

	/* collect number of stars created */
	MPI_Allgather( &num_new_stars, 1, MPI_INT, proc_new_stars, 1, MPI_INT, MPI_COMM_WORLD );

	/* find how many "blocks" to expect */
	max_stars = 0;
	total_new_stars = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( proc_new_stars[proc] > max_stars ) {
			max_stars = proc_new_stars[proc];
		}

		total_new_stars += proc_new_stars[proc];
	}

	if ( total_new_stars > 0 ) {
		/* create lists of indices for each block */
		block_ids = cart_alloc( max_stars * sizeof(int) );

		block_ids[0] = 0;
		for ( block = 1; block < max_stars; block++ ) {
			block_ids[block] = block_ids[block-1];
			for ( proc = 0; proc < num_procs; proc++ ) {
				if ( proc_new_stars[proc] >= block ) {
					block_ids[block]++;
				}
			}
		}
	
		/* find all newly allocated stars and remap their id's (keeping order) */
		for ( ipart = 0; ipart < num_star_particles; ipart++ ) {
			if ( particle_level[ipart] != FREE_PARTICLE_LEVEL && 
					particle_id[ipart] >= particle_species_indices[num_particle_species] ) {
	
				block = ( particle_id[ipart] - particle_species_indices[num_particle_species] ) / num_procs;
					proc = ( particle_id[ipart] - particle_species_indices[num_particle_species] ) % num_procs;
				new_id = particle_species_indices[num_particle_species] + block_ids[block];
	
				for ( i = 0; i < proc; i++ ) {
					if ( proc_new_stars[i] > block ) {
							new_id++;
					}
				}
				
				cart_assert( new_id <= particle_id[ipart] && 
					new_id < particle_species_indices[num_particle_species]+total_new_stars );

				particle_id[ipart] = new_id;
			}
		}

		cart_free(block_ids);
	}

	particle_species_indices[num_particle_species] += total_new_stars;
	particle_species_num[num_particle_species-1] += total_new_stars;
	num_particles_total += total_new_stars;
	num_new_stars = 0;
}

#endif /* STARFORM */
