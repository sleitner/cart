#include "config.h"

#include <math.h>

#include "auxiliary.h"
#include "min_heap.h"
#include "particle.h"
#include "particle_density.h"
#include "sfc.h"
#include "stack.h"
#include "timing.h"
#include "tree.h"

#ifdef PARTICLES

double smooth_kernel( double r, double h ) {
	double k;
	double rh = r/h;

	if ( rh < 1.0 ) {
		k = M_1_PI/(h*h*h)*(1-0.75*(2.-rh)*(rh*rh));
	} else if ( rh < 2.0 ) {
		k = 0.25*M_1_PI/(h*h*h)*(2.-rh)*(2.-rh)*(2.-rh);
	} else {
		k = 0;
	}

	return k;
}

min_heap *particle_heap;
stack *cell_list;

void particle_find_nearest_neighbors( int ipart, int *particle_flag, int num_nearest_neighbors, 
		double max_search_radius, int *num_neighbors_found, int *neighbors, double *distances ) {
	int i, j, k;
	int icell, child;
	int icell_leaf;
	int icell_root;
	int level;
	int iterpart;
	int ioct;
	double r, rmax, rprune;
	double delta, pos[nDim];
	int coords[nDim];
	double rmax_factor = 2.0;

	min_heap *particle_heap = min_heap_init( num_nearest_neighbors );
	stack *cell_list = stack_init();

	icell_leaf = cell_find_position( particle_x[ipart] );
	icell_root = cell_find_position_level( min_level, particle_x[ipart] );

	do {
		min_heap_flush(particle_heap);

		/* try pre-seeding local particles */
		iterpart = cell_particle_list[icell_leaf];
		while ( iterpart != NULL_PARTICLE ) {
			if ( particle_flag[iterpart] ) {
				r = compute_distance_periodic( particle_x[ipart], particle_x[iterpart] );
				min_heap_push( particle_heap, iterpart, r );
			}
			iterpart = particle_list_next[iterpart];
		}

		if ( min_heap_size( particle_heap ) < num_nearest_neighbors ) {
			/* estimate rmax */
			rmax = rmax_factor*pow( (double)num_nearest_neighbors /
					MAX( (double)min_heap_size( particle_heap ), 1.0 ) * 
					cell_volume[cell_level(icell_leaf)], 1./3.);
		} else {
			min_heap_peek( particle_heap, &j, &rmax );
		}

		rmax = MIN( rmax, max_search_radius );
		rprune = rmax*rmax;

		/* find root cells within rmax */
		for ( i = (int)floor(particle_x[ipart][0]-rmax); i <= (int)(particle_x[ipart][0]+rmax); i++ ) {
			coords[0] = ( i + num_grid ) % num_grid;
			for ( j = (int)floor(particle_x[ipart][1]-rmax); j <= (int)(particle_x[ipart][1]+rmax); j++ ) {
				coords[1] = ( j + num_grid ) % num_grid;
				for ( k = (int)floor(particle_x[ipart][2]-rmax); k <= (int)(particle_x[ipart][2]+rmax); k++ ) {
					coords[2] = ( k + num_grid ) % num_grid;
					icell = root_cell_location( sfc_index( coords ) );
					if ( icell != NULL_OCT && icell != icell_root ) {
						stack_push( cell_list, 	icell );
					}
				}
			}
		}
		stack_push( cell_list, icell_root );

		while ( stack_pop( cell_list, &icell ) ) {
			if ( icell == icell_leaf ) {
				continue;
			}

			iterpart = cell_particle_list[icell];
			while ( iterpart != NULL_PARTICLE ) {
				if ( particle_flag[iterpart] ) {
					r = compute_distance_periodic( particle_x[ipart], particle_x[iterpart] );
					if ( r <= rmax ) {
						min_heap_push( particle_heap, iterpart, r );
					}
				}
				iterpart = particle_list_next[iterpart];
			}

			if ( cell_is_refined( icell ) ) {
				if ( min_heap_size(particle_heap) == num_nearest_neighbors ) {
					min_heap_peek( particle_heap, &j, &r );
					rmax = r;
					rprune = rmax*rmax;
				}

                ioct = cell_child_oct[icell];
				level = oct_level[ioct];
				for ( j = 0; j < nDim; j++ ) {
					pos[j] = particle_x[ipart][j] - oct_pos[ioct][j];
				}

				for ( i = 0; i < num_children; i++ ) {
					child = oct_child(ioct,i);

					r = 0.0;
					for ( j = 0; j < nDim; j++ ) {
						delta = fabs( pos[j] - cell_size[level]*cell_delta[i][j] );

						if ( delta > (double)(num_grid/2) ) {
							delta = delta - (double)num_grid + 0.5*cell_size[level];
							r += delta*delta;
						} else if ( delta >= 0.5*cell_size[level] ) {
							delta -= 0.5*cell_size[level];
							r += delta*delta;
						}
					}

					if ( r < rprune ) {
						stack_push( cell_list, child );
					}
				}
			}
		}

		rmax_factor *= 2.0;
	} while ( min_heap_size(particle_heap) < num_nearest_neighbors && rmax < max_search_radius );

	*num_neighbors_found = min_heap_size(particle_heap);
	for ( i = 0; i < *num_neighbors_found; i++ ) {
		min_heap_pop( particle_heap, &neighbors[i], &distances[i] );
	}

	min_heap_destroy( particle_heap );
	stack_destroy( cell_list );
}

void find_neighbors_slow( int ipart, int *particle_flag, int *num_neighbors_found, 
		int *neighbors, double *distances ) {
	int i;
	double r;

	min_heap_flush( particle_heap );
	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_flag[i] ) {
			r = compute_distance_periodic( particle_x[ipart], particle_x[i] );
			min_heap_push( particle_heap, i, r );
		}
	}

	*num_neighbors_found = min_heap_size(particle_heap);
	for ( i = 0; i < *num_neighbors_found; i++ ) {
		min_heap_pop( particle_heap, &neighbors[i], &distances[i] );
	}
}

void compute_particle_densities( int num_nearest_neighbors, int *particle_flag, float *particle_density ) {
	int i;
	double h;
	int ipart;

	int num_neighbors_found;
	int *neighbors;
	double *distances;
	double W;

	start_time( SMOOTH_PARTICLE_DENSITY_TIMER );

	neighbors = cart_alloc(int, num_nearest_neighbors);
	distances = cart_alloc(double, num_nearest_neighbors);

#ifndef COMPILER_GCC
		/* Get compiler segfault under GCC */
#pragma omp parallel default(none) shared(particle_flag,particle_density,particle_mass,num_nearest_neighbors) private(ipart,num_neighbors_found,neighbors,distances,h,W,i)
#endif
	{
		neighbors = cart_alloc(int, num_nearest_neighbors);	
		distances = cart_alloc(double, num_nearest_neighbors);

		#pragma omp for schedule(dynamic)
		for ( ipart = 0; ipart < num_particles; ipart++ ) {
			if ( particle_flag[ipart] == 2 ) {
				/* add to particle list */
				particle_find_nearest_neighbors( ipart, particle_flag, num_nearest_neighbors, 
						2.0, &num_neighbors_found, neighbors, distances );

				if ( num_neighbors_found == num_nearest_neighbors ) {
					h = 0.5*distances[0];
				} else {
					h = 1.0;
				}

				/* mark additional particles needed to do densities */
				for ( i = 0; i < num_neighbors_found; i++ ) {
					W = smooth_kernel( distances[i], h ); 
					#pragma omp critical
					particle_density[ipart] += 0.5*particle_mass[neighbors[i]]*W;

					if ( particle_flag[neighbors[i]] == 2 ) {
						#pragma omp critical
						particle_density[neighbors[i]] += 0.5*particle_mass[ipart]*W;
					} else if ( particle_flag[neighbors[i]] == 1 ) {
						particle_flag[neighbors[i]] = 3;
					}
				}
			}
		}

		#pragma omp for schedule(dynamic)
		for ( ipart = 0; ipart < num_particles; ipart++ ) {
			if ( particle_flag[ipart] == 3 ) {
				particle_find_nearest_neighbors( ipart, particle_flag, num_nearest_neighbors, 
						2.0, &num_neighbors_found, neighbors, distances );

				if ( num_neighbors_found == num_nearest_neighbors ) {
					h = 0.5*distances[0];
				} else {
					h = 1.0;
				}

				for ( i = 0; i < num_neighbors_found; i++ ) {
					if ( particle_flag[neighbors[i]] == 2 ) {
						W = smooth_kernel( distances[i], h );
						#pragma omp critical
						particle_density[neighbors[i]] += 0.5*particle_mass[ipart]*W;
					}
				}
			}
		}

		cart_free( neighbors );
		cart_free( distances );
	}

	end_time( SMOOTH_PARTICLE_DENSITY_TIMER );

#ifdef DEBUG
	cart_debug("density time = %f", total_time( SMOOTH_PARTICLE_DENSITY_TIMER, -1 ) );
#endif
}

#ifdef PARTICLE_VDISP
/* this code is untested, won't compile, and shouldn't be used */
void compute_particle_vdisp() {
	int i, j;
	int ipart;
	double h, W, v2;
	int num_nearest_neighbors = 24;
	int num_neighbors_found;
	int neighbors[24];
	double distances[24];

	cell_list = stack_init();
	particle_heap = min_heap_init( num_nearest_neighbors );

	/* compute density */
	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_flag[ipart] ) {
			particle_find_nearest_neighbors( ipart, num_nearest_neighbors,
					2.0, &num_neighbors_found, neighbors, distances );

			if ( num_neighbors_found == num_nearest_neighbors ) {
				h = 0.5*distances[0];
			} else {
				h = 1.0;
			}

			for ( i = 0; i < num_neighbors_found; i++ ) {
				W = smooth_kernel( distances[i], h ); 
				particle_density[ipart] += 0.5*particle_mass[neighbors[i]]*W;
				particle_density[neighbors[i]] += 0.5*particle_mass[ipart]*W;
			}
		}
	}

	/* compute mean velocity at each particle position */
	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_flag[ipart] ) {
			particle_find_nearest_neighbors( ipart, num_nearest_neighbors,
					2.0, &num_neighbors_found, neighbors, distances );
			
			if ( num_neighbors_found == num_nearest_neighbors ) {
				h = 0.5*distances[0];
			} else {
				h = 1.0;
			}

			for ( i = 0; i < num_neighbors_found; i++ ) {
				W = 0.5*smooth_kernel( distances[i], h );

				for ( j = 0; j < nDim; j++ ) {
					particle_vmean[ipart][j] += particle_v[neighbors[i]][j]*
						particle_mass[neighbors[i]]*W / 
						particle_density[neighbors[i]];
				}

				for ( j = 0; j < nDim; j++ ) {
					particle_vmean[neighbors[i]][j] += particle_v[ipart][j]*
						particle_mass[ipart]*W /
						particle_density[ipart];
				}
			}
		}
	}

	/* compute velocity dispersion */
    for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_flag[ipart] ) {
			particle_find_nearest_neighbors( ipart, num_nearest_neighbors, 
					2.0, &num_neighbors_found, neighbors, distances );

			if ( num_neighbors_found == num_nearest_neighbors ) {
				h = 0.5*distances[0];
			} else {
				h = 1.0;
			}
			
			for ( i = 0; i < num_neighbors_found; i++ ) {
				W = 0.5*smooth_kernel( distances[i], h );

				v2 = 0.0;
				for ( j = 0; j < nDim; j++ ) {
					v2 += (particle_v[ipart][j] - particle_vmean[ipart][j] ) *
						(particle_v[ipart][j] - particle_vmean[ipart][j] );	
				}
				
				particle_vdisp[ipart] += v2*particle_mass[neighbors[i]]*W/
						particle_density[neighbors[i]];

				v2 = 0.0;
				for ( j = 0; j < nDim; j++ ) {
					v2 += (particle_v[neighbors[i]][j] - particle_vmean[neighbors[i]][j] ) *
						(particle_v[neighbors[i]][j] - particle_vmean[neighbors[i]][j] );  
				}
                
				particle_vdisp[neighbors[i]] += v2*particle_mass[ipart]*W/
						particle_density[ipart];	
			}
        }
    }
}
#endif /* PARTICLE_VDISP */

#endif /* PARTICLES */
