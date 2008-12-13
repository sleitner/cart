#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "io.h"
#include "hydro.h"
#include "timestep.h"
#include "tree.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "timing.h"
#include "viewdump.h"
#include "gravity.h"
#include "units.h"
#include "sfc.h"
#include "auxiliary.h"
#include "cooling.h"

#include "defs.h"
#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#endif


#ifdef HYDRO 

float pressure_floor_factor = 50.0;
float pressure_floor;
float gas_density_floor = 1e-6;

float backup_hvars[num_cells][num_hydro_vars-2];
float ref[num_cells];

#ifdef _OPENMP
double backup_fluxes[num_cells][num_hydro_vars-1];
#endif

#ifdef GRAVITY_IN_RIEMANN
void fluxh( double dtx, double dtx2, double v[num_hydro_vars-1][4], double g[2], double c[2], double f[num_hydro_vars-1] );
#else
void fluxh( double dtx, double dtx2, double v[num_hydro_vars-1][4], double c[2], double f[num_hydro_vars-1] );
#endif

#ifdef LAPIDUS
void lapidus( double dtx2, int L1, int R1, int sweep_direction, int mj3, int mj4, int mj5, double v[num_hydro_vars-1][4], double f[num_hydro_vars-1] );
#endif

int level_sweep_dir[max_level-min_level+1];
const int sweep_dir[2][nDim] = { { 1, 3, 5 }, { 5, 3, 1 } };

int sweep_direction;
int sweep_dimension;

int j3,j4,j5;
int mj3, mj4, mj5;
double dtx;
double dtx2;

void hydro_step( int level ) {
	int i, j;
	int dir;
	int L1, L2, R1, R2;
	double gravadd;
	int icell;
	int num_level_cells;
	int *level_cells;
	double dxi;
	double f[num_hydro_vars-1];

	const int momentum_permute[2*nDim][nDim] = {  { 0, 1, 2 }, { 0, 1, 2 },
						{ 1, 0, 2 }, { 1, 0, 2 },
						{ 2, 1, 0 }, { 2, 1, 0 } };


#ifdef PRESSURE_FLOOR
	if ( level >= MinL_Jeans ) {
		/* artificial pressure floor */
		pressure_floor = 0.47746 * pressure_floor_factor * aexp[level] * cell_size[max_level]*cell_size[max_level];
	} else {
		pressure_floor = 0.0;
	}
#else
	pressure_floor = 0.0;
#endif

        dtx = dtl[level] * cell_size_inverse[level];
	dxi = cell_size_inverse[level];
        dtx2 = 0.5*dtx;

	start_time( HYDRO_TIMER );

	/* zero ref's of cells on level with NO_SPLIT_NEIGHBORS */
	hydro_copy_vars( level, COPY_ZERO_REF, COPY_NO_SPLIT_NEIGHBORS );

	for ( dir = 0; dir < nDim; dir++ ) {
		/* H_Old_to_New( Level, 2 ) in Fortran ver */
		hydro_copy_vars( level, COPY, COPY_NO_SPLIT_NEIGHBORS );

		start_time( WORK_TIMER );

		sweep_direction = sweep_dir[level_sweep_dir[level]][dir];
		sweep_dimension = (sweep_direction-1)/2;

		j3 = momentum_permute[sweep_direction][0];
		j4 = momentum_permute[sweep_direction][1];
		j5 = momentum_permute[sweep_direction][2];

		mj3 = j3+2;
		mj4 = j4+2;
		mj5 = j5+2;

		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

#ifdef _OPENMP
#pragma omp parallel for default(none), private(i,icell,j), shared(num_level_cells,level_cells,backup_fluxes)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			for ( j = 0; j < num_hydro_vars-1; j++ ) {
				backup_fluxes[icell][j] = 0.0;
			}
		}
#endif /* openmp */

#pragma omp parallel for default(none), private(i,icell,R1,R2,L1,L2,f,j), shared(num_level_cells,level_cells,cell_child_oct,sweep_direction,backup_hvars,mj3,mj4,ref,dxi,level,mj5,backup_fluxes)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_leaf( icell ) ) {
				/* calculate neighbors */
				L1 = cell_neighbor( icell, reverse_direction[sweep_direction]);
				R1 = cell_neighbor( icell, sweep_direction );
				R2 = cell_neighbor( R1, sweep_direction );

				if ( cell_is_leaf( R1 ) ) {
					compute_hydro_fluxes( L1, icell, R1, R2, f );

					/* apply fluxes immediately for local, left cell */
					backup_hvars[icell][0] -= (float)f[0];
					backup_hvars[icell][1] -= (float)f[4];
					backup_hvars[icell][mj3] -= (float)f[1];
					backup_hvars[icell][mj4] -= (float)f[2];
					backup_hvars[icell][mj5] -= (float)f[3];
					backup_hvars[icell][5] -= (float)f[5];
					ref[icell] -= (float)f[6]*dxi;

#ifdef ADVECT_SPECIES
					for ( j = 0; j < num_chem_species; j++ ) {
						backup_hvars[icell][j+6] -= (float)f[j+7];
					}
#endif /* ADVECT_SPECIES */

#ifdef _OPENMP
					if ( cell_is_local(R1) ) {
						if ( cell_level(R1) < level ) {
							/* again apply fluxes immediately */
							backup_hvars[R1][0] += (float)f[0]*0.125;
							backup_hvars[R1][1] += (float)f[4]*0.125;
							backup_hvars[R1][mj3] += (float)f[1]*0.125;
							backup_hvars[R1][mj4] += (float)f[2]*0.125;
							backup_hvars[R1][mj5] += (float)f[3]*0.125;
							backup_hvars[R1][5] += (float)f[5]*0.125;
							ref[R1] += (float)f[6]*0.5*dxi*0.125;

#ifdef ADVECT_SPECIES
							for ( j = 0; j < num_chem_species; j++ ) {
								backup_hvars[R1][j+6] += (float)f[j+7]*0.125;	
							}
#endif /* ADVECT_SPECIES */

						} else {
							/* save fluxes for later */
							for ( j = 0; j < num_hydro_vars-1; j++ ) {
								backup_fluxes[R1][j] = f[j];
							}
						}
					}
#else
					if ( cell_is_local(R1) ) {
						if ( cell_level(R1) < level ) {
							backup_hvars[R1][0] += (float)f[0]*0.125;
							backup_hvars[R1][1] += (float)f[4]*0.125;
							backup_hvars[R1][mj3] += (float)f[1]*0.125;
							backup_hvars[R1][mj4] += (float)f[2]*0.125;
							backup_hvars[R1][mj5] += (float)f[3]*0.125;
							backup_hvars[R1][5] += (float)f[5]*0.125;
							ref[R1] += (float)f[6]*0.5*dxi*0.125;

#ifdef ADVECT_SPECIES
							for ( j = 0; j < num_chem_species; j++ ) {
								backup_hvars[R1][j+6] += (float)f[j+7]*0.125;
							}
#endif /* ADVECT_SPECIES */

						} else {
							backup_hvars[R1][0] += (float)f[0];
							backup_hvars[R1][1] += (float)f[4];
							backup_hvars[R1][mj3] += (float)f[1];
							backup_hvars[R1][mj4] += (float)f[2];
							backup_hvars[R1][mj5] += (float)f[3];
							backup_hvars[R1][5] += (float)f[5];
							ref[R1] += (float)f[6]*dxi;

#ifdef ADVECT_SPECIES
							for ( j = 0; j < num_chem_species; j++ ) {
								backup_hvars[R1][j+6] += (float)f[j+7];
							}
#endif /* ADVECT_SPECIES */

						}
					}
#endif /* OPENMP */
				}
                                                                                                                                                            
				if ( cell_level(L1) == level - 1 || ( cell_is_leaf(L1) && !cell_is_local(L1) ) ) {
					L2 = cell_neighbor( L1, reverse_direction[sweep_direction] );
					compute_hydro_fluxes( L2, L1, icell, R1, f );

					/* apply fluxes immediately for local, left cell (since level(L1) < level(cell)
					 * it won't ever be an R1, so it's okay to apply directly */
					backup_hvars[icell][0] += (float)f[0];
					backup_hvars[icell][1] += (float)f[4];
					backup_hvars[icell][mj3] += (float)f[1];
					backup_hvars[icell][mj4] += (float)f[2];
					backup_hvars[icell][mj5] += (float)f[3];
					backup_hvars[icell][5] += (float)f[5];
					ref[icell] += (float)f[6]*dxi;

#ifdef ADVECT_SPECIES
					for ( j = 0; j < num_chem_species; j++ ) {
						backup_hvars[icell][j+6] += (float)f[j+7];
					}
#endif /* ADVECT_SPECIES */


					/* apply directly since cell_level(L1) < level */
					if ( cell_is_local(L1) ) {
						backup_hvars[L1][0] -= (float)f[0]*0.125;
						backup_hvars[L1][1] -= (float)f[4]*0.125;
						backup_hvars[L1][mj3] -= (float)f[1]*0.125;
						backup_hvars[L1][mj4] -= (float)f[2]*0.125;
						backup_hvars[L1][mj5] -= (float)f[3]*0.125;
						backup_hvars[L1][5] -= (float)f[5]*0.125;
						ref[L1] -= (float)f[6]*0.5*dxi*0.125;

#ifdef ADVECT_SPECIES
						for ( j = 0; j < num_chem_species; j++ ) {
							backup_hvars[L1][j+6] -= (float)f[j+7]*0.125;
						}
#endif /* ADVECT_SPECIES */

					}
				}
			} 
		}

#ifdef _OPENMP
		/* apply fluxes for right cells */
#pragma omp parallel for default(none), private(i,icell,j), shared(num_level_cells,level_cells,backup_fluxes,mj3,mj4,mj5,ref,dxi,backup_hvars)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			backup_hvars[icell][0] += (float)backup_fluxes[icell][0];
			backup_hvars[icell][1] += (float)backup_fluxes[icell][4];
			backup_hvars[icell][mj3] += (float)backup_fluxes[icell][1];
			backup_hvars[icell][mj4] += (float)backup_fluxes[icell][2];
			backup_hvars[icell][mj5] += (float)backup_fluxes[icell][3];
			backup_hvars[icell][5] += (float)backup_fluxes[icell][5];
			ref[icell] += (float)backup_fluxes[icell][6]*dxi;

#ifdef ADVECT_SPECIES
			for ( j = 0; j < num_chem_species; j++ ) {
				backup_hvars[icell][j+6] += (float)backup_fluxes[icell][j+7];
			}
#endif /* ADVECT_SPECIES */

		}
#endif /* OPENMP */

		cart_free( level_cells );

		if ( level > min_level ) {
			select_level( level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,R1,R2,L1,L2,f,j), shared(num_level_cells,level_cells,cell_child_oct,sweep_direction,level,backup_hvars,mj3,mj4,mj5,dxi,ref)
			for ( i = 0; i < num_level_cells; i++ ) {
				icell = level_cells[i];

				if ( cell_is_leaf(icell) ) {
					R1 = cell_neighbor( icell, sweep_direction );
					L1 = cell_neighbor( icell, reverse_direction[sweep_direction]);

					if ( R1 != -1 && cell_level(R1) == level - 1 && cell_is_local(R1) ) {
						R2 = cell_neighbor( R1, sweep_direction );
						compute_hydro_fluxes( L1, icell, R1, R2, f );

						if ( cell_is_local(R1) ) {
							backup_hvars[R1][0] += (float)f[0]*0.125;
							backup_hvars[R1][1] += (float)f[4]*0.125;
							backup_hvars[R1][mj3] += (float)f[1]*0.125;
							backup_hvars[R1][mj4] += (float)f[2]*0.125;
							backup_hvars[R1][mj5] += (float)f[3]*0.125;
							backup_hvars[R1][5] += (float)f[5]*0.125;
							ref[R1] += (float)f[6]*0.5*dxi*0.125;

#ifdef ADVECT_SPECIES
							for ( j = 0; j < num_chem_species; j++ ) {
								backup_hvars[R1][j+6] += (float)f[j+7]*0.125;
							}
#endif /* ADVECT_SPECIES */

						}
					}

					if ( L1 != -1 && cell_level(L1) == level - 1 && cell_is_local(L1) ) {
						L2 = cell_neighbor( L1, reverse_direction[sweep_direction] );
						compute_hydro_fluxes( L2, L1, icell, R1, f );

						if ( cell_is_local(L1) ) {
							backup_hvars[L1][0] -= (float)f[0]*0.125;
							backup_hvars[L1][1] -= (float)f[4]*0.125;
							backup_hvars[L1][mj3] -= (float)f[1]*0.125;
							backup_hvars[L1][mj4] -= (float)f[2]*0.125;
							backup_hvars[L1][mj5] -= (float)f[3]*0.125;
							backup_hvars[L1][5] -= (float)f[5]*0.125;
							ref[L1] -= (float)f[6]*0.5*dxi*0.125;

#ifdef ADVECT_SPECIES
							for ( j = 0; j < num_chem_species; j++ ) {
								backup_hvars[L1][j+6] -= (float)f[j+7]*0.125;
							}
#endif /* ADVECT_SPECIES */
						}
					}
				}
			}
			cart_free( level_cells );
		}

#ifdef GRAVITY 
#ifndef GRAVITY_IN_RIEMANN
		/* now we need to apply a gravity correction */
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(icell,gravadd), shared(num_level_cells,level_cells,cell_child_oct,backup_hvars,cell_vars,sweep_dimension,mj3)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_leaf(icell) ) {
				gravadd = backup_hvars[icell][0] * cell_accel(icell,sweep_dimension);
				backup_hvars[icell][1] += cell_accel(icell,sweep_dimension) *
							( backup_hvars[icell][mj3] + 0.5 * gravadd );
				backup_hvars[icell][mj3] += gravadd;
			}
		}
		cart_free( level_cells );
#endif /* GRAVITY_IN_RIEMANN */
#endif /* GRAVITY */

		end_time( WORK_TIMER );

		if ( dir == nDim - 1 ) {
			hydro_copy_vars( level, RESTORE, COPY_ALL_LEAFS );
		} else {
			hydro_copy_vars( level, RESTORE, 
					COPY_NO_SPLIT_NEIGHBORS );

			hydro_magic( level );
			hydro_eos( level );

			start_time( HYDRO_UPDATE_TIMER );
			update_buffer_level( level, all_hydro_vars, num_hydro_vars );
			end_time( HYDRO_UPDATE_TIMER );
		}
	}

	hydro_advance_internalenergy( level );
	hydro_magic( level );
	hydro_eos( level );

	hydro_split_update( level );

	start_time( HYDRO_UPDATE_TIMER );
	update_buffer_level( level, all_hydro_vars, num_hydro_vars );
	end_time( HYDRO_UPDATE_TIMER );

	/* update sweep direction */
	level_sweep_dir[level] = (level_sweep_dir[level]+1)%2;

	end_time( HYDRO_TIMER );
}

void hydro_magic( int level ) {
	int i, j;
	int icell;
	int num_level_cells;
	int *level_cells;
	float average_density;
	int neighbors[num_neighbors];

	double Tminc;
	double kinetic_energy;
	double thermal_energy;

	Tminc = T_min * aexp[level]*aexp[level] / ( T0 * ( gamma - 1.0 ) );

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,j,kinetic_energy,thermal_energy,average_density,neighbors), shared(num_level_cells,level_cells,cell_child_oct,cell_vars,gas_density_floor,Tminc)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_leaf(icell) ) {
			/* do density floor stuff */
			if ( cell_gas_density(icell) < gas_density_floor ) {
				average_density = 0.0;
				cell_all_neighbors( icell, neighbors );
	
				for ( j = 0; j < num_neighbors; j++ ) {
					average_density += cell_gas_density( neighbors[j] );
				}

				cart_debug("cell %u hit density floor, old density %e, new density = %e",
					icell, cell_gas_density(icell), max( average_density/(float)num_neighbors,
					gas_density_floor ) );
#if defined(DEBUG) && (DEBUG-0 > 9)
				cart_error("aborting...");
#endif

				for ( i = 0; i < num_hydro_vars; i++ ) {
					cart_debug("hydro var %u: %e", i, cell_hydro_variable( icell, i ) );
				}	
	
				cell_gas_density(icell) = max( average_density/(float)num_neighbors, 
								gas_density_floor );
			}

			kinetic_energy = cell_gas_kinetic_energy(icell);
			thermal_energy = Tminc * cell_gas_density(icell);

			cell_gas_internal_energy(icell) = max( cell_gas_internal_energy(icell), thermal_energy );
			cell_gas_energy(icell) = max( cell_gas_energy(icell), thermal_energy+kinetic_energy );

#ifdef ADVECT_SPECIES
			for ( j = 0; j < num_chem_species; j++ ) {
			  /* 
			     1e-15 may be too large a number for ionic species;
			     at least let's scale them with density and make 1e-20;
                             Gnedin: 1e-20 is not small enough for chemistry, making it 1e-30
			  */
				cell_advected_variable(icell,j) = max( 1e-30*cell_gas_density(icell), cell_advected_variable(icell,j) );
				/* Doug had it like that:
				cell_advected_variable(icell,j) = max( 1e-15, cell_advected_variable(icell,j) );
				*/
			}
#endif /* ADVECT_SPECIES */
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
}

void hydro_eos( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;
	double kinetic_energy;

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,i,kinetic_energy), shared(num_level_cells,level_cells,cell_child_oct,cell_vars)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_leaf( icell) ) {
			kinetic_energy = cell_gas_kinetic_energy(icell);

#ifdef PRESSURELESS_FLUID
			cell_gas_pressure(icell) = 1e-20;
			cell_gas_internal_energy(icell) = cell_gas_pressure(icell) / (gamma-1.0);
			cell_gas_energy(icell) = cell_gas_internal_energy(icell) + kinetic_energy;
#else
			cell_gas_internal_energy(icell) = max( cell_gas_internal_energy(icell), 0.0 );
			cell_gas_energy(icell) = max( kinetic_energy, cell_gas_energy(icell) );
			cell_gas_pressure(icell) = max( (cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell), 0.0 );
#endif
		}
	}

	cart_free( level_cells );

	end_time( WORK_TIMER );
}


#ifdef COOLING

void hydro_cool_one_cell(int icell, double t_begin, double t_end, double a1, double afact, double Zdum, double rhogl, double rhog2, double fact_nH, double Tfact_cell, double e_small) {
	int continue_cooling;
	double t_curr, a_curr;
	double ai;
	double dE;
	double dt_e, ei1, T_gas;
	double Tfact_cell2;

#define dstep	(0.01)

	continue_cooling = 1;
	t_curr = t_begin;
	a_curr = a1;

	/* integrate cooling using smaller timestep */
	while ( ( t_curr < t_end ) && continue_cooling ) {
		ai = 1.0 / a_curr;
		Tfact_cell2 = Tfact_cell * ai * ai;
		T_gas = Tfact_cell2 * cell_gas_internal_energy(icell);

		/* compute new timestep */
		dE = cooling_rate( rhogl, T_gas, Zdum );
		dE *= -rhog2 * a_curr;
		dt_e = min( dstep * fabs( cell_gas_internal_energy(icell) / dE ), t_end - t_curr );

		ei1 = max( cell_gas_internal_energy(icell) + 0.5 * dE * dt_e, e_small );
		T_gas = Tfact_cell2 * ei1;

		dE = cooling_rate( rhogl, T_gas, Zdum );
		dE *= -rhog2 * a_curr * dt_e;

		/* adjust cell energies */
		cell_gas_internal_energy(icell) += dE;
		cell_gas_energy(icell) += dE;

		/* stop if we hit energy minimum */
		if ( cell_gas_internal_energy(icell) < e_small ) {
			continue_cooling = 0;
		}

		cell_gas_internal_energy(icell) = max( e_small, cell_gas_internal_energy(icell) );
		cell_gas_energy(icell) = max( e_small, cell_gas_energy(icell) );

		/* advance timestep */
		t_curr += dt_e;
		a_curr = a1 + afact * ( t_curr - t_begin );
	}
}

void hydro_apply_cooling(int level, int num_level_cells, int *level_cells) {
	int i;
	int icell;
	double t_begin, t_end;
	double ai, a1, a2, afact;
	double e_min, e_small;
	double Zdum;
	double fact_nH, Tfact, Tfact_cell;
	double rhogi, rhog2, rhogl;

	t_begin	= tl[level];
	t_end = tl[level] + dtl[level];

#ifdef COSMOLOGY
	ai = 1.0 / aexp[level];
	a1 = aexp[level];
	a2 = b2a( t_end );
	afact = ( a2 - a1 ) / dtl[level];
#else
	a1 = a2 = ai = 1.0;
	afact = 0.0;
#endif

	fact_nH = log10( 1.12e-5 * hubble * hubble * Omega0 * ( 1.0 - Y_p ) * ai*ai*ai );
	Tfact = T0 * ( gamma - 1.0 ) / 1e4;
	e_min = T_min * a1 * a1  / T0 / ( gamma - 1.0 );

#pragma omp parallel for default(none), private(icell,i,rhogi,rhog2,rhogl,Zdum,Tfact_cell,e_small), shared(num_level_cells,level_cells,t_begin,t_end,fact_nH,Tfact,e_min,cell_child_oct,cell_vars,a1,afact)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		if ( cell_is_leaf(icell) ) {
		
			cell_gas_gamma(icell) = gamma;

			rhogi = 1.0 / cell_gas_density(icell);
			rhog2 = cell_gas_density(icell)*cell_gas_density(icell);

#ifdef CLOUDY_COOLING
			/* take code density -> log10(n_H [cm^-3]) */
			rhogl = log10(cell_gas_density(icell)) + fact_nH;
#endif /* CLOUDY_COOLING */

#ifdef METALCOOLING
			Zdum = cell_gas_metallicity(icell);
			Zdum = max( Zdum, 1e-10 );
			Zdum = log10( Zdum * rhogi / Zsolar );
#else
			Zdum = 0.0;
#endif /* METALCOOLING*/

			Tfact_cell = Tfact * rhogi;
			e_small = e_min * cell_gas_density(icell);

			hydro_cool_one_cell(icell,t_begin,t_end,a1,afact,Zdum,rhogl,rhog2,fact_nH,Tfact_cell,e_small);
		}
	}
}

#endif /* COOLING */

void hydro_advance_internalenergy( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;
	double kinetic_energy;
	double energy;
	double gamma1, div, div_dt;

	start_time( WORK_TIMER );

	div_dt = dtl[level] / 3.0;

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(icell,i,kinetic_energy,energy,gamma1,div), shared(num_level_cells,level_cells,cell_child_oct,cell_vars,ref,div_dt)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_leaf(icell) ) {
			/* P dV term */
			gamma1 = cell_gas_gamma(icell) - 1.0;
			div = 1.0 + gamma1 * ref[icell] * div_dt;
			cell_gas_internal_energy(icell) = max( small, cell_gas_internal_energy(icell)*div*div*div);

			/* synchronize internal and total energy */
			kinetic_energy = cell_gas_kinetic_energy(icell);
			energy = cell_gas_energy(icell);

			/* we trust energy over internal energy since it's computed using
			 * the riemann solver rather than just advection equation, so if
			 * internal energy is sufficiently large then compute it from 
			 * e = E - rho * v**2 /2 */
			if ( ( energy - kinetic_energy) / energy > 1e-3 ) {
                                cell_gas_internal_energy(icell) = energy - kinetic_energy;
                        }
		}
	}

#ifdef COOLING
#ifdef RADIATIVE_TRANSFER
	start_time( RT_COOLING_TIMER );
	rtApplyCooling(level,num_level_cells,level_cells);
	end_time( RT_COOLING_TIMER );
#else
	start_time( COOLING_TIMER );
	hydro_apply_cooling(level,num_level_cells,level_cells);
	end_time( COOLING_TIMER );
#endif /* RADIATIVE_TRANSFER */
#endif /* COOLING */

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

		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,j,k,children,new_var), shared(num_level_cells,level_cells,cell_child_oct,cell_vars)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_refined(icell) ) {
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
		}
		cart_free( level_cells );

		end_time( WORK_TIMER );
	}
}

void compute_hydro_fluxes( int L2, int L1, int R1, int R2, double f[num_hydro_vars-1] ) {
        int j;
	double v[num_hydro_vars-1][4];
	double c[2];

#ifdef GRAVITY_IN_RIEMANN
	double g[2];
#endif

	cart_assert( cell_is_leaf(L1) && cell_is_leaf(R1) );

	/* L2 */
	v[0][0] = cell_gas_density(L2);
	v[1][0] = max( pressure_floor * v[0][0]*v[0][0], cell_gas_pressure(L2) );
	v[2][0] = cell_momentum(L2,j3)/cell_gas_density(L2);
	v[3][0] = cell_momentum(L2,j4)/cell_gas_density(L2);
	v[4][0] = cell_momentum(L2,j5)/cell_gas_density(L2);
	v[5][0] = cell_gas_gamma(L2);
	v[6][0] = gamma;

#ifdef ADVECT_SPECIES
	for ( j = 0; j < num_chem_species; j++ ) {
		v[j+7][0] = cell_advected_variable(L2,j);
	}
#endif /* ADVECT_SPECIES */
	
	/* L1 vars */
	v[0][1] = cell_gas_density(L1);
	v[1][1] = max( pressure_floor * v[0][1]*v[0][1], cell_gas_pressure(L1) );
	v[2][1] = cell_momentum(L1,j3)/cell_gas_density(L1);
	v[3][1] = cell_momentum(L1,j4)/cell_gas_density(L1);
	v[4][1] = cell_momentum(L1,j5)/cell_gas_density(L1);
	v[5][1] = cell_gas_gamma(L1);
	v[6][1] = gamma;

#ifdef ADVECT_SPECIES
	for ( j = 0; j < num_chem_species; j++ ) {
		v[j+7][1] = cell_advected_variable(L1,j);
	}
#endif /* ADVECT_SPECIES */
	
	/* R1 vars */
	v[0][2] = cell_gas_density(R1);
	v[1][2] = max( pressure_floor * v[0][2]*v[0][2], cell_gas_pressure(R1) );
	v[2][2] = cell_momentum(R1,j3)/cell_gas_density(R1);
	v[3][2] = cell_momentum(R1,j4)/cell_gas_density(R1);
	v[4][2] = cell_momentum(R1,j5)/cell_gas_density(R1);
	v[5][2] = cell_gas_gamma(R1);
	v[6][2] = gamma;

#ifdef ADVECT_SPECIES
	for ( j = 0; j < num_chem_species; j++ ) {
		v[j+7][2] = cell_advected_variable(R1,j);
	}
#endif /* ADVECT_SPECIES */
	
	/* R2 vars */
	v[0][3] = cell_gas_density(R2);
	v[1][3] = max( pressure_floor * v[0][3]*v[0][3], cell_gas_pressure(R2) );
	v[2][3] = cell_momentum(R2,j3)/cell_gas_density(R2);
	v[3][3] = cell_momentum(R2,j4)/cell_gas_density(R2);
	v[4][3] = cell_momentum(R2,j5)/cell_gas_density(R2);
	v[5][3] = cell_gas_gamma(R2);
	v[6][3] = gamma;

#ifdef ADVECT_SPECIES
	for ( j = 0; j < num_chem_species; j++ ) {
		v[j+7][3] = cell_advected_variable(R2,j);
	}
#endif /* ADVECT_SPECIES */


	if ( cell_level(R1) > cell_level(L1) ) {
		c[0] = 1.0/1.5;
		c[1] = 1.0/1.25;
	} else if ( cell_level(R1) < cell_level(L1) ) {
		c[0] = 1.0/1.25;
		c[1] = 1.0/1.5;
	} else {
		if ( cell_level( L2 ) == cell_level(L1) ) {
			c[0] = 1.0;
		} else {
			c[0] = 1.0/1.25;
		}

		if ( cell_level( R2 ) == cell_level(L1) ) {
			c[1] = 1.0;
		} else {
			c[1] = 1.0/1.25;
		}
	} 

#ifdef GRAVITY_IN_RIEMANN
	g[0] = 0.5*cell_accel( L1, j3 );
	g[1] = 0.5*cell_accel( R1, j3 );

	/* compute fluxes */
	fluxh( dtx, dtx2, v, g, c, f );
#else
	fluxh( dtx, dtx2, v, c, f );
#endif

#ifdef LAPIDUS
	lapidus( dtx2, L1, R1, sweep_direction, j3, j4, j5, v, f );
#endif

}
	
void hydro_copy_vars( int level, int direction, int copy_cells ) {
	int i, j;
	int neighbors[num_neighbors];
	int do_copy;
	int icell;
	int num_level_cells;
	int *level_cells;

#if nDim != 3
	#error	hydro_copy_vars only works for nDim = 3
#endif

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,j,neighbors,do_copy), shared(num_level_cells,level_cells,cell_child_oct,cell_vars,copy_cells,direction,backup_hvars,ref)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_leaf(icell) ) {
			do_copy = 1;

			if ( copy_cells != COPY_ALL_LEAFS ) {
				cell_all_neighbors( icell, neighbors );
	
				if ( copy_cells == COPY_SPLIT_NEIGHBORS ) {
					do_copy = 0;
					for ( j = 0; j < num_neighbors; j++ ) {
						if ( cell_is_refined(neighbors[j]) ) {
							do_copy = 1;
							break;
						}
					}
				} else if ( copy_cells == COPY_NO_SPLIT_NEIGHBORS ) {
					for ( j = 0; j < num_neighbors; j++ ) {
						if ( cell_is_refined(neighbors[j]) ) {
							do_copy = 0;
							break;
						}
					}
				}
			}

			if ( do_copy ) {
				if ( direction == COPY || direction == COPY_ZERO_REF ) {
					backup_hvars[icell][0] = cell_gas_density(icell);
					backup_hvars[icell][1] = cell_gas_energy(icell);
					backup_hvars[icell][2] = cell_momentum(icell,0);
					backup_hvars[icell][3] = cell_momentum(icell,1);
					backup_hvars[icell][4] = cell_momentum(icell,2);
					backup_hvars[icell][5] = cell_gas_internal_energy(icell);

#ifdef ADVECT_SPECIES
					for ( j = 0; j < num_chem_species; j++ ) {
						backup_hvars[icell][6+j] = cell_advected_variable(icell,j);
					}
#endif /* ADVECT_SPECIES */
		
					if ( direction == COPY_ZERO_REF ) {
						ref[icell] = 0.0;
					}
				} else {
					cell_gas_density(icell) = max( small, backup_hvars[icell][0] );
					cell_gas_energy(icell) = max( small, backup_hvars[icell][1] );
					cell_momentum(icell,0) = backup_hvars[icell][2];
					cell_momentum(icell,1) = backup_hvars[icell][3];
					cell_momentum(icell,2) = backup_hvars[icell][4];
					cell_gas_internal_energy(icell) = max( small, backup_hvars[icell][5] );
					
#ifdef ADVECT_SPECIES
					for ( j = 0; j < num_chem_species; j++ ) {
						cell_advected_variable(icell,j) = max( small, backup_hvars[icell][6+j] );
					}
#endif /* ADVECT_SPECIES */
				}
			}
		}
	}

	cart_free( level_cells );

	end_time( WORK_TIMER );
}

#endif /*HYDRO*/
