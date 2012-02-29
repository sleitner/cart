#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "particle.h"
#include "hydro_tracer.h"
#include "starformation.h"
#include "tree.h"


int compare_oct_parent_root_sfc( const void *a, const void *b ) {
	int ioct1 = *(int *)a;
	int ioct2 = *(int *)b;

	return ( oct_parent_root_sfc[ioct1] - oct_parent_root_sfc[ioct2] );
}

void cache_reorder_tree() {
	int i;
	int level;
	int ioct, ioct_new, ioct_old;
	int icell, icell_new;
	int parent, child;
	int ivar;
	int oct_level_indices[max_level+2];
	int *oct_index;
	int *oct_order;
	float *cell_scratch_vars;
#ifdef PARTICLES
	int *scratch_particle_list;
#endif
#ifdef HYDRO_TRACERS
	int *scratch_tracer_list;
#endif
	int *scratch_cell_child_oct;
	int num_local_octs;
	int oct_count;
	int first_oct;

	cart_debug("optimizing tree for cache efficiency");

	first_oct = cell_parent_oct( num_cells_per_level[min_level] + num_buffer_cells[min_level] ) + 1;

	num_local_octs = 0;
	oct_level_indices[min_level+1] = 0;
	for ( level = min_level+1; level <= max_level; level++ ) {
		num_local_octs += num_cells_per_level[level] / num_children;
		oct_level_indices[level+1] = num_local_octs;

		cart_assert( oct_level_indices[level+1] < num_octs );
	}

	/* since oct_index is last to be freed, allocate early */
	oct_index = cart_alloc(int, num_octs );
        for ( ioct = 0; ioct < num_octs; ioct++ ) {
                oct_index[ioct] = NULL_OCT;
        }

	oct_order = cart_alloc(int, num_local_octs );

	/* create mapping of new oct index to old oct index */
	oct_count = 0;
	for ( level = min_level+1; level <= max_level; level++ ) {
		ioct = local_oct_list[level];

		while ( ioct != NULL_OCT ) {
			oct_order[oct_count++] = ioct;
			ioct = oct_next[ioct];
		}

		cart_assert( oct_count == oct_level_indices[level+1] );
	}

	/* sort oct order */
	for ( level = min_level+1; level <= max_level; level++ ) {
		if ( num_cells_per_level[level] > 0 ) {
			qsort( &oct_order[ oct_level_indices[level] ], 
				oct_level_indices[level+1]-oct_level_indices[level], sizeof(int), 
				compare_oct_parent_root_sfc );
		}
	}

	for ( level = min_level; level <= max_level; level++ ) {
		oct_level_indices[level+1] += first_oct;
	}

	/* create oct_index mapping from old index to new packed oct index */
	for ( ioct_new = 0; ioct_new < num_local_octs; ioct_new++ ) {
		oct_index[ oct_order[ioct_new] ] = ioct_new + first_oct;
	}

	cart_free( oct_order );

	/* convert cell_child_oct pointers and copy over variables */
	for ( icell = 0; icell < num_cells_per_level[min_level]; icell++ ) {
		if ( cell_is_refined(icell) ) {
			ioct_new = oct_index[cell_child_oct[icell]];
			cell_child_oct[icell] = ioct_new;
		}
	}

	scratch_cell_child_oct = cart_alloc(int, num_children * num_local_octs );

	for ( level = min_level+1; level <= max_level; level++ ) {
		ioct = local_oct_list[level];

		while ( ioct != NULL_OCT ) {
			ioct_new = oct_index[ioct];

			for ( child = 0; child < num_children; child++ ) {
				icell = oct_child( ioct, child );
				icell_new = oct_child( ioct_new, child ) - num_children*first_oct;

				if ( cell_is_refined(icell) ) {
					scratch_cell_child_oct[icell_new] = oct_index[ cell_child_oct[icell] ];
				} else {
					scratch_cell_child_oct[icell_new] = UNREFINED_CELL;
				}
			}

			ioct = oct_next[ioct];
		}
	}

	for ( level = min_level+1; level <= max_level; level++ ) {
		for ( ioct = oct_level_indices[level]; ioct < oct_level_indices[level+1]; ioct++ ) {
			for ( child = 0; child < num_children; child++ ) {
				icell = oct_child( ioct, child );
				cell_child_oct[icell] = scratch_cell_child_oct[icell - num_children * first_oct];
			}
		}
	}

	for ( icell = num_children * ( first_oct + num_local_octs ); icell < num_cells; icell++ ) {
		cell_child_oct[icell] = UNREFINED_CELL;
	}

	cart_free( scratch_cell_child_oct );

	/* copy over variables one at a time (could do more, this ensures low memory usage) */
	cell_scratch_vars = cart_alloc(float, num_children * num_local_octs );

	for ( ivar = 0; ivar < num_vars; ivar++ ) {
		/* copy variable to scratch */
		for ( level = min_level+1; level <= max_level; level++ ) {
			ioct = local_oct_list[level];

			while ( ioct != NULL_OCT ) {
				ioct_new = oct_index[ioct];
	
				for ( child = 0; child < num_children; child++ ) {
					icell = oct_child( ioct, child );
					icell_new = oct_child( ioct_new, child );
					cell_scratch_vars[icell_new-num_children*first_oct] = cell_var( icell, ivar );
				}
		
				ioct = oct_next[ioct];
			}
		}

		/* now copy variable back */
		for ( level = min_level+1; level <= max_level; level++ ) {
			for ( ioct = oct_level_indices[level]; ioct < oct_level_indices[level+1]; ioct++ ) {
				for ( child = 0; child < num_children; child++ ) {
					icell = oct_child( ioct, child );
					cell_var( icell, ivar ) = cell_scratch_vars[icell-num_children*first_oct];
				}
			}
		}
	}

	cart_free( cell_scratch_vars );

#ifdef PARTICLES
	scratch_particle_list = cart_alloc(int, num_children * num_local_octs );

	for ( level = min_level+1; level <= max_level; level++ ) {
		ioct = local_oct_list[level];

		while ( ioct != NULL_OCT ) {
			ioct_new = oct_index[ioct];

			for ( child = 0; child < num_children; child++ ) {
				icell = oct_child( ioct, child );
				icell_new = oct_child( ioct_new, child );

				scratch_particle_list[icell_new-num_children*first_oct] = cell_particle_list[icell];
			}

			ioct = oct_next[ioct];
		}
	}

	/* now copy variable back */
	for ( level = min_level+1; level <= max_level; level++ ) {
		for ( ioct = oct_level_indices[level]; ioct < oct_level_indices[level+1]; ioct++ ) {
			for ( child = 0; child < num_children; child++ ) {
				icell = oct_child( ioct, child );
				cell_particle_list[icell] = scratch_particle_list[icell-num_children*first_oct];
			}
		}
	}

	cart_free( scratch_particle_list );

	for ( icell = num_children * ( first_oct + num_local_octs); icell < num_cells; icell++ ) {
		cell_particle_list[icell] = NULL_PARTICLE;
	}
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	scratch_tracer_list = cart_alloc(int, num_children * num_local_octs );

	for ( level = min_level+1; level <= max_level; level++ ) {
		ioct = local_oct_list[level];

		while ( ioct != NULL_OCT ) {
			ioct_new = oct_index[ioct];

			for ( child = 0; child < num_children; child++ ) {
				icell = oct_child( ioct, child );
				icell_new = oct_child( ioct_new, child );

				scratch_tracer_list[icell_new-num_children*first_oct] = cell_tracer_list[icell];
			}

			ioct = oct_next[ioct];
		}
	}

	/* now copy variable back */
	for ( level = min_level+1; level <= max_level; level++ ) {
		for ( ioct = oct_level_indices[level]; ioct < oct_level_indices[level+1]; ioct++ ) {
			for ( child = 0; child < num_children; child++ ) {
				icell = oct_child( ioct, child );
				cell_tracer_list[icell] = scratch_tracer_list[icell-num_children*first_oct];
			}
		}
	}

	cart_free( scratch_tracer_list );

	for ( icell = num_children * ( first_oct + num_local_octs ); icell < num_cells; icell++ ) {
		cell_tracer_list[icell] = NULL_TRACER;
	}
#endif /* HYDRO_TRACERS */

	/* pack oct values into oct_neighbors */
	for ( level = min_level+1; level <= max_level; level++ ) {
		ioct = local_oct_list[level];

		while ( ioct != NULL_OCT ) {
			ioct_new = oct_index[ioct];
                                                                                                                                                            
			/* use neighbors to store values we need to copy */
			oct_neighbors[ioct_new][0] = oct_parent_cell[ioct];
			oct_neighbors[ioct_new][1] = oct_parent_root_sfc[ioct];

			ioct = oct_next[ioct];
		}
	}

	/* copy over oct values
         * neighbors are recomputed in load balance */
	for ( level = min_level+1; level <= max_level; level++ ) {
		for ( ioct = oct_level_indices[level]; ioct < oct_level_indices[level+1]; ioct++ ) {
			/* adjust parent cell if necessary (level > root ) */
			parent = oct_neighbors[ioct][0];
			if ( cell_level(parent) > min_level ) {
				/* convert cell by converting parent cell's oct index */
				ioct_old = cell_parent_oct( parent );
				child = cell_child_number( parent );
				ioct_new = oct_index[ ioct_old ];
				parent = oct_child( ioct_new, child );
			}

			oct_parent_cell[ioct] = parent;
			oct_parent_root_sfc[ioct] = oct_neighbors[ioct][1];
			oct_level[ioct] = level;
			
			/* this requires levels fixed from min->max */
			cell_center_position( parent, oct_pos[ioct] );

			for ( i = 0; i < num_neighbors; i++ ) {
				oct_neighbors[ioct][i] = NULL_OCT;
			}

			/* initialize oct linked list */
			oct_next[ioct] = ioct+1;
			oct_prev[ioct] = ioct-1;
		}

		/* complete oct linked list */
		if ( oct_level_indices[level] == oct_level_indices[level+1] ) {
			local_oct_list[level] = NULL_OCT;
		} else {
			local_oct_list[level] = oct_level_indices[level];
			oct_prev[ oct_level_indices[level] ] = NULL_OCT;
			oct_next[ oct_level_indices[level+1] - 1 ] = NULL_OCT;
		}
	}

	cart_free( oct_index );

	/* reset rest of tree */
        for ( ioct = oct_level_indices[max_level+1]; ioct < num_octs; ioct++ ) {
		oct_parent_cell[ioct] = -1;
		oct_parent_root_sfc[ioct] = -1;
		oct_level[ioct] = FREE_OCT_LEVEL;
		oct_next[ioct] = NULL_OCT;
		oct_prev[ioct] = NULL_OCT;
        }

	next_free_oct = num_local_octs + first_oct;
	free_oct_list = NULL_OCT;
}

#ifdef PARTICLES
int *particle_parent_cell;

int sort_particles( const void *a, const void *b ) {
	int index_a = *(int *)a;
	int index_b = *(int *)b;

	cart_assert( index_a >= 0 && index_a < num_particles );
	cart_assert( index_b >= 0 && index_b < num_particles );

	if ( particle_parent_cell[index_a] < particle_parent_cell[index_b] ) {
		return -1;
	} else if ( particle_parent_cell[index_a] > particle_parent_cell[index_b] ) {
		return 1;
	} else {
		return ( particle_id[index_a] - particle_id[index_b] );
	}
}

void cache_reorder_particles() {
	int i, d;
	int index, star_index;
	int level;
	int icell;
	int ipart;
	int num_level_cells;
	int num_normal_particles;
	int *level_cells;
	int *backup_ints;
	float *backup_floats;
	double *backup_doubles;
	int *old_particle_index;
	int *particle_order;
	int *star_particle_order;

	cart_debug("beginning cache optimization for particles");

	old_particle_index = cart_alloc(int, num_particles );
	particle_parent_cell = cart_alloc(int, num_particles );

	for ( index = 0; index < num_particles; index++ ) {
		old_particle_index[index] = -1;
	}

#ifdef STARFORM
	num_normal_particles = num_local_particles-num_local_star_particles;
	star_particle_order = cart_alloc(int, num_local_star_particles );
	particle_order = cart_alloc(int, num_normal_particles );

	index = 0;
	star_index = 0;

	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells ); 

		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			ipart = cell_particle_list[icell];

			while ( ipart != NULL_PARTICLE ) {
				particle_parent_cell[ipart] = icell;

				if ( particle_is_star(ipart) ) {
					star_particle_order[star_index++] = ipart;
				} else {
					particle_order[index++] = ipart;
				}

				ipart = particle_list_next[ipart];
			}
		}

		cart_free( level_cells );
	}

	cart_assert( star_index == num_local_star_particles );
	cart_assert( index == num_normal_particles );

	/* sort star particles */
	qsort( star_particle_order, num_local_star_particles, sizeof(int), 
		sort_particles );

	/* assign star particles new positions at beginning of particle array */
	for ( index = 0; index < num_local_star_particles; index++ ) {
		old_particle_index[index] = star_particle_order[index];
	}

	cart_free( star_particle_order );

	/* sort normal particles */
	qsort( particle_order, num_normal_particles, sizeof(int), sort_particles );

	if ( num_normal_particles > ( num_particles - num_star_particles ) ) {
		/* need to use some star particle spots for normal particles */
		index = num_particles - num_normal_particles;
	} else {
		/* pack particles immediately after star spots */
		index = num_star_particles;
	}
	cart_assert( index >= num_local_star_particles );

	for ( i = 0; i < num_normal_particles; i++ ) {
		old_particle_index[index+i] = particle_order[i];
	}

	cart_free( particle_order );

	/* set up allocation variables */
	next_free_particle = index+num_normal_particles;
	free_particle_list = NULL_PARTICLE;

	next_free_star_particle = num_local_star_particles;
	free_star_particle_list = NULL_PARTICLE;
#else
	particle_order = cart_alloc(int, num_local_particles );

	index = 0;
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			ipart = cell_particle_list[icell];
			while ( ipart != NULL_PARTICLE ) {
				particle_parent_cell[ipart] = icell;
				particle_order[index++] = ipart;
				ipart = particle_list_next[ipart];
			}
		}

		cart_free( level_cells );
	}

	cart_assert( index == num_local_particles );

	/* sort particles */
	qsort( particle_order, num_local_particles, sizeof(int), sort_particles );

	for ( index = 0; index < num_local_particles; index++ ) {
		old_particle_index[index] = particle_order[index];
	}

	cart_free( particle_order );

	/* set up allocation variables */
	next_free_particle = num_local_particles;
	free_particle_list = NULL_PARTICLE;
#endif /* STARFORM */

	/* now put particles into proper order */
	backup_doubles = cart_alloc(double, num_particles );

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_doubles[ipart] = particle_t[ old_particle_index[ipart] ];
		} else {
			backup_doubles[ipart] = 0.0;
		}
	}

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		particle_t[ipart] = backup_doubles[ipart];
	}

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_doubles[ipart] = particle_dt[ old_particle_index[ipart] ];
		} else {
			backup_doubles[ipart] = 0.0;
		}
	}

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		particle_dt[ipart] = backup_doubles[ipart];
	}

	cart_free( backup_doubles );

	backup_floats = cart_alloc(float, num_particles );

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_floats[ipart] = particle_mass[ old_particle_index[ipart] ];
		} else {
			backup_floats[ipart] = 0.0;
		} 
	}

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		particle_mass[ipart] = backup_floats[ipart];
	}

#ifdef STARFORM
	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_floats[ipart] = star_tbirth[ old_particle_index[ipart] ];
		} else {
			backup_floats[ipart] = 0.0;
		}
	}

	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		star_tbirth[ipart] = backup_floats[ipart];
	}

	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_floats[ipart] = star_initial_mass[ old_particle_index[ipart] ];
		} else {
			backup_floats[ipart] = 0.0;
		}
	}

	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		star_initial_mass[ipart] = backup_floats[ipart];
	}

#ifdef ENRICH
	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_floats[ipart] = star_metallicity_II[ old_particle_index[ipart] ];
		} else {
			backup_floats[ipart] = 0.0;
		}
	}

	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		star_metallicity_II[ipart] = backup_floats[ipart];
	}
#endif /* ENRICH */

#ifdef ENRICH_SNIa
	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_floats[ipart] = star_metallicity_Ia[ old_particle_index[ipart] ];
		} else {
			backup_floats[ipart] = 0.0;
		}
	}

	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		star_metallicity_Ia[ipart] = backup_floats[ipart];
	}
#endif /* ENRICH_SNIa */

#endif /* STARFORM */

	cart_free( backup_floats );

	backup_ints = cart_alloc(int, num_particles );

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_ints[ipart] = particle_id[ old_particle_index[ipart] ];
		} else {
			backup_ints[ipart] = NULL_PARTICLE;
		}
	}

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		particle_id[ipart] = backup_ints[ipart];
	}

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_ints[ipart] = particle_level[ old_particle_index[ipart] ];
		} else {
			backup_ints[ipart] = FREE_PARTICLE_LEVEL;
		}
	}

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		particle_level[ipart] = backup_ints[ipart];
	}

#if defined(STARFORM) && defined(STAR_PARTICLE_TYPES)
	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		if ( old_particle_index[ipart] != -1 ) {
			backup_ints[ipart] = star_particle_type[ old_particle_index[ipart] ];
		} else {
			backup_ints[ipart] = STAR_TYPE_DELETED;
		}
	}

	for ( ipart = 0; ipart < num_local_star_particles; ipart++ ) {
		star_particle_type[ipart] = backup_ints[ipart];
	}
#endif /* STARFORM && STAR_PARTICLE_TYPES */

	cart_free( backup_ints );

	backup_doubles = cart_alloc(double, num_particles );

	for ( d = 0; d < nDim; d++ ) {
		for ( ipart = 0; ipart < num_particles; ipart++ ) {
			if ( old_particle_index[ipart] != -1 ) {
				backup_doubles[ipart] = particle_x[ old_particle_index[ipart] ][d];
			} else {
				backup_doubles[ipart] = 0.0;
			}
		}

		for ( ipart = 0; ipart < num_particles; ipart++ ) {
			particle_x[ipart][d] = backup_doubles[ipart];
		}
	}

	for ( d = 0; d < nDim; d++ ) {
		for ( ipart = 0; ipart < num_particles; ipart++ ) {
			if ( old_particle_index[ipart] != -1 ) {
				backup_doubles[ipart] = particle_v[ old_particle_index[ipart] ][d];
			} else {
				backup_doubles[ipart] = 0.0;
			}
		}

		for ( ipart = 0; ipart < num_particles; ipart++ ) {
			particle_v[ipart][d] = backup_doubles[ipart];
		}
	}

	cart_free( backup_doubles );

	for ( icell = 0; icell < num_cells; icell++ ) {
		cell_particle_list[icell] = NULL_PARTICLE;
	}

	for ( ipart = num_particles-1; ipart >= 0; ipart-- ) {
		particle_list_prev[ipart] = NULL_PARTICLE;

		if ( old_particle_index[ipart] != -1 ) {
			icell = particle_parent_cell[ old_particle_index[ipart] ];
			particle_list_next[ipart] = cell_particle_list[icell];

			if ( cell_particle_list[icell] != NULL_PARTICLE ) {
				cart_assert( cell_particle_list[icell] > ipart );
				particle_list_prev[ cell_particle_list[icell] ] = ipart;
			}

			cell_particle_list[icell] = ipart;
		} else {
			particle_list_next[ipart] = NULL_PARTICLE;
		}
	}

	cart_free( old_particle_index );
	cart_free( particle_parent_cell );
}
#endif /* PARTICLES */
