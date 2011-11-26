#include "config.h"
#ifdef GRAVITY 


#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "iterators.h"
#include "times.h"
#include "timing.h"
#include "tree.h"

#include "step.h"


#ifdef HYDRO

void copy_potential( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;

	start_time( GRAVITY_TIMER );
	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_ANY, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,cell_vars)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		cell_potential_hydro(icell) = cell_potential(icell);
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
	end_time( GRAVITY_TIMER );
}

void interpolate_potential( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;
	double dtdt2;

	start_time( GRAVITY_TIMER );
	start_time( WORK_TIMER );

	/*
	//  NG: dtl_old may not be set in the first time-step, but then
	//  cell_potential = cell_potential_hydro
	*/
	if(dtl_old[level] > 0.1*dtl[level])
	  {
	    dtdt2 = 0.5 * dtl[level]/dtl_old[level];
	  }
	else
	  {
	    dtdt2 = 0;
	  }

	select_level( level, CELL_TYPE_ANY, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,cell_vars,dtdt2)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		cell_potential_hydro(icell) = cell_potential(icell) +
			( cell_potential(icell) - cell_potential_hydro(icell) ) * dtdt2;
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
	end_time( GRAVITY_TIMER );
}

void compute_accelerations_hydro( int level ) {
	int i, j;
	double a2half;
	int neighbors[num_neighbors];
	int L1, R1;
	double phi_l, phi_r;
#ifdef GRAVITY_IN_RIEMANN
	const int accel_vars[nDim] = { VAR_ACCEL, VAR_ACCEL+1, VAR_ACCEL+2 };
#endif

	int icell;
	int num_level_cells;
	int *level_cells;

	start_time( GRAVITY_TIMER );
	start_time( HYDRO_ACCEL_TIMER );
	start_time( WORK_TIMER );

#ifdef COSMOLOGY
	a2half = abox_from_tcode( tl[level] + 0.5*dtl[level] );
	a2half = -0.5*dtl[level]*cell_size_inverse[level]*a2half*a2half;
#else
	a2half = -0.5*dtl[level]*cell_size_inverse[level];
#endif 

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#ifdef OPENMP_DECLARE_CONST
#pragma omp parallel for default(none), private(icell,j,neighbors,L1,R1,phi_l,phi_r), shared(num_level_cells,level_cells,level,cell_vars,a2half), shared(local)
#else  /* OPENMP_DECLARE_CONST */
#pragma omp parallel for default(none), private(icell,j,neighbors,L1,R1,phi_l,phi_r), shared(num_level_cells,level_cells,level,cell_vars,a2half)
#endif /* OPENMP_DECLARE_CONST */
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		cell_all_neighbors( icell, neighbors );
		for ( j = 0; j < nDim; j++ ) {
			L1 = neighbors[2*j];
			R1 = neighbors[2*j+1];

			if ( cell_level(L1) == level && cell_level(R1) == level ) {
				phi_l = cell_potential_hydro(L1);
				phi_r = cell_potential_hydro(R1);
			} else {
				if ( cell_level(L1) < level ) {
					phi_l = cell_interpolate( L1, local[cell_child_number(icell)][2*j], VAR_POTENTIAL );
					phi_r = cell_potential(R1);
				} else {
					phi_l = cell_potential(L1);
					phi_r = cell_interpolate( R1, local[cell_child_number(icell)][2*j], VAR_POTENTIAL );
				}
			}

			cell_accel( icell, j ) = (float)(a2half * ( phi_r - phi_l ) );
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );

#ifdef GRAVITY_IN_RIEMANN
	/* this only gets called if we pass gravity on to Riemann solver (and thus need accel in buffer cells) */
	start_time( HYDRO_ACCEL_UPDATE_TIMER );
	update_buffer_level( level, accel_vars, nDim );
	end_time( HYDRO_ACCEL_UPDATE_TIMER );
#endif

	end_time( HYDRO_ACCEL_TIMER );
	end_time( GRAVITY_TIMER );
}

#endif /* HYDRO */

#ifdef PARTICLES

void compute_accelerations_particles( int level ) {
	int i, j;
	double a2half;
	const int accel_vars[nDim] = { VAR_ACCEL, VAR_ACCEL+1, VAR_ACCEL+2 };
	int neighbors[num_neighbors];
	int L1, R1;
	double phi_l, phi_r;
	int icell;
	int num_level_cells;
	int *level_cells;

	start_time( GRAVITY_TIMER );
	start_time( PARTICLE_ACCEL_TIMER );
	start_time( WORK_TIMER );

#ifdef COSMOLOGY
	a2half = -0.5*abox[level]*abox[level]*cell_size_inverse[level];
#else
	a2half = -0.5 * cell_size_inverse[level];
#endif

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(icell,j,neighbors,L1,R1,phi_l,phi_r), shared(num_level_cells,level_cells,level,cell_vars,a2half)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		cell_all_neighbors( icell, neighbors );
		for ( j = 0; j < nDim; j++ ) {
			L1 = neighbors[2*j];
			R1 = neighbors[2*j+1];

			if ( cell_level(L1) < level ) {
				phi_l = 0.8*cell_potential(L1) + 0.2*cell_potential(icell);
			} else {
				phi_l = cell_potential(L1);
			}

			if ( cell_level(R1) < level ) {
				phi_r = 0.8*cell_potential(R1)+0.2*cell_potential(icell);
			} else {
				phi_r = cell_potential(R1);
			}

			cell_accel( icell, j ) = (float)(a2half * ( phi_r - phi_l ) );
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );

	/* update accelerations */
	start_time( PARTICLE_ACCEL_UPDATE_TIMER );
	update_buffer_level( level, accel_vars, nDim );
	end_time( PARTICLE_ACCEL_UPDATE_TIMER );

	end_time( PARTICLE_ACCEL_TIMER );
	end_time( GRAVITY_TIMER );
}

#endif /* PARTICLES */

#endif /* GRAVITY */
