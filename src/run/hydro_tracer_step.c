#include "config.h"
#if defined(HYDRO) && defined(HYDRO_TRACERS)

#include <math.h>

#include "auxiliary.h"
#include "hydro_tracer.h"
#include "iterators.h"
#include "timing.h"
#include "tree.h"

#include "step.h"


void move_hydro_tracers( int level ) {
	int i, j;
	int tracer;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	double vdt[nDim];
	int icell, icell_orig;
	int level1;
	int child;
	double pos[nDim];
	int found;
	int c[num_children];
	double diff1, diff2, diff3;
	double pt3, pd3;
	double t1,t2,t3,d1,d2,d3;
	double t2t1, t2d1, d2t1, d2d1;
	double t3t2t1, t3t2d1, t3d2t1, t3d2d1;
	double d3t2t1, d3t2d1, d3d2t1, d3d2d1;

	start_time( WORK_TIMER ); 

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		iter_cell = level_cells[i];

#ifdef HYDRO_TRACERS_NGP
		for ( j = 0; j < nDim; j++ ) {
			vdt[j] = cell_momentum(iter_cell,j)/cell_gas_density(iter_cell) * dtl[level];
		}
#endif /* HYDRO_TRACERS_NGP */

		tracer = cell_tracer_list[iter_cell];
		while ( tracer != NULL_TRACER ) {
			cart_assert( tracer >= 0 && tracer < num_tracers );

#ifndef HYDRO_TRACERS_NGP
			icell = iter_cell;
			level1 = level;

			do {
				found = 1;
				icell_orig = icell;
				cart_assert( icell != NULL_OCT );

				cell_center_position( icell, pos );

				/* find lower leftmost cell */
				child = 0;
				for ( j = 0; j < nDim; j++ ) {
					if ( tracer_x[tracer][j] >= pos[j] ) {
						child += (1<<j);
					}
				}

				cart_assert( child >= 0 && child < num_children );

				for ( j = 0; j < nDim; j++ ) {
					if ( neighbor_moves[child][j] == -1 ) {
						break;
					} else {
						icell = cell_neighbor(icell, neighbor_moves[child][j] );
						cart_assert( icell != NULL_OCT );

						if ( cell_level(icell) != level1 ) {
							icell = cell_parent_cell(icell_orig);
							cart_assert( icell != NULL_OCT );
							level1 = level1 - 1;
							found = 0;
							break;
						}
					}
				}

				if ( found ) {
					c[0] = icell;
					c[1] = cell_neighbor(icell,1);
					c[2] = cell_neighbor(icell,3);
					c[3] = cell_neighbor(c[1],3);
					c[4] = cell_neighbor(icell,5);
					c[5] = cell_neighbor(c[1],5);
					c[6] = cell_neighbor(c[2],5);
					c[7] = cell_neighbor(c[3],5);

					for ( j = 1; j < num_children; j++ ) {
						if ( cell_level(c[j]) != level1 ) {
							icell = cell_parent_cell(icell_orig);
							level1 = level1 - 1;
							cart_assert( icell != NULL_OCT );
							found = 0;
							break;
						}
					}
				}
			} while ( !found );

			cell_center_position( c[0], pos );

			/* now we have the level on which this particle will move */
			diff1 = pos[0] - tracer_x[tracer][0];
			if ( fabs(diff1) > (double)(num_grid/2) ) {
				if ( diff1 > 0.0 ) {
					diff1 -= (double)(num_grid);
				} else {
					diff1 += (double)(num_grid);
				}
			}
			d1 = fabs(diff1) * cell_size_inverse[level1];
			cart_assert( d1 >= 0.0 && d1 <= 1.0 );

			diff2 = pos[1] - tracer_x[tracer][1];
			if ( fabs(diff2) > (double)(num_grid/2) ) {
				if ( diff2 > 0.0 ) {
					diff2 -= (double)(num_grid);
				} else {
					diff2 += (double)(num_grid);
				}
			}
			d2 = fabs(diff2) * cell_size_inverse[level1];

			diff3 = pos[2] - tracer_x[tracer][2];
			if ( fabs(diff3) > (double)(num_grid/2) ) {
				if ( diff3 > 0.0 ) {
					diff3 -= (double)(num_grid);
				} else {
					diff3 += (double)(num_grid);
				}
			}
			d3 = fabs(diff3) * cell_size_inverse[level1];

			cart_assert( d1 >= 0.0 && d1 <= 1.0 );
			cart_assert( d2 >= 0.0 && d2 <= 1.0 );
			cart_assert( d3 >= 0.0 && d3 <= 1.0 );

			t1   = 1.0 - d1;
			t2   = 1.0 - d2;
			t3   = 1.0 - d3;

			cart_assert( t1 >= 0.0 && t1 <= 1.0 );
			cart_assert( t2 >= 0.0 && t2 <= 1.0 );
			cart_assert( t3 >= 0.0 && t3 <= 1.0 );

			t2t1 = t2 * t1;
			t2d1 = t2 * d1;
			d2t1 = d2 * t1;
			d2d1 = d2 * d1;

			pt3 = t3*dtl[level];
			pd3 = d3*dtl[level];

			t3t2t1 = pt3 * t2t1;
			t3t2d1 = pt3 * t2d1;
			t3d2t1 = pt3 * d2t1;
			t3d2d1 = pt3 * d2d1;
			d3t2t1 = pd3 * t2t1;
			d3t2d1 = pd3 * t2d1;
			d3d2t1 = pd3 * d2t1;
			d3d2d1 = pd3 * d2d1;

			for ( j = 0; j < nDim; j++ ) {
				vdt[j] =t3t2t1 * cell_momentum(c[0], j) / cell_gas_density(c[0]) +
					t3t2d1 * cell_momentum(c[1], j) / cell_gas_density(c[1]) +
					t3d2t1 * cell_momentum(c[2], j) / cell_gas_density(c[2]) +
					t3d2d1 * cell_momentum(c[3], j) / cell_gas_density(c[3]) +
					d3t2t1 * cell_momentum(c[4], j) / cell_gas_density(c[4]) +
					d3t2d1 * cell_momentum(c[5], j) / cell_gas_density(c[5]) +
					d3d2t1 * cell_momentum(c[6], j) / cell_gas_density(c[6]) +
					d3d2d1 * cell_momentum(c[7], j) / cell_gas_density(c[7]);
			}
#endif /* HYDRO_TRACERS_NGP */

			tracer_x[tracer][0] += vdt[0];
			tracer_x[tracer][1] += vdt[1];
			tracer_x[tracer][2] += vdt[2];

			/* enforce periodic boundaries */
			if ( tracer_x[tracer][0] < 0.0 ) {
				tracer_x[tracer][0] += (double)(num_grid);		
			}
				
			if ( tracer_x[tracer][0] >= (double)(num_grid) ) {
				tracer_x[tracer][0] -= (double)(num_grid);
			}

			if ( tracer_x[tracer][1] < 0.0 ) {
				tracer_x[tracer][1] += (double)(num_grid);
			}

			if ( tracer_x[tracer][1] >= (double)(num_grid) ) {
				tracer_x[tracer][1] -= (double)(num_grid);
			}

			if ( tracer_x[tracer][2] < 0.0 ) {
				tracer_x[tracer][2] += (double)(num_grid);
			}

			if ( tracer_x[tracer][2] >= (double)(num_grid) ) {
				tracer_x[tracer][2] -= (double)(num_grid);
			}

			tracer = tracer_list_next[tracer];
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
}

#endif /* HYDRO && HYDRO_TRACERS */
