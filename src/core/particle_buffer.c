#include "config.h"

#include <math.h>                                                                                                                                                            
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "pack.h"
#include "parallel.h"
#include "sfc.h"
#include "skiplist.h"
#include "tree.h"
#include "particle.h"
#include "starformation.h"

#ifdef PARTICLES
int particle_buffer_enabled = 0;

void build_particle_buffer( int species, int subspecies ) {
	int i;
	int level, proc;
	int ipart, icell;
	int count;
	int num_level_cells, *level_cells;
	int num_parts_to_send[MAX_PROCS];
	int *particle_list_to_send[MAX_PROCS];

	cart_assert( !particle_buffer_enabled );

	/* determine number of particles to send */
	for ( proc = 0; proc < num_procs; proc++ ) {
		num_parts_to_send[proc] = 0;

		for ( level = min_level; level <= max_level; level++ ) {
			select_local_buffered_cells( level, proc, &num_level_cells, &level_cells );

			for ( i = 0; i < num_level_cells; i++ ) { 
				icell = level_cells[i];
				ipart = cell_particle_list[icell];
				while ( ipart != NULL_PARTICLE ) {
					if ( species == -1 || 
							( particle_species(particle_id[ipart]) == species
#ifdef STAR_PARTICLE_TYPES
							  && ( species != num_particle_species-1 || subspecies == -1 ||
								  star_particle_type[ipart] == subspecies ) 
#endif /* STAR_PARTICLE_TYPES */
							) ) {
						num_parts_to_send[proc]++;
					}
					ipart = particle_list_next[ipart];
				}
			}
			cart_free( level_cells );
		}

		if ( num_parts_to_send[proc] > 0 ) {
			particle_list_to_send[proc] = cart_alloc(int, num_parts_to_send[proc]);

			count = 0;
			for ( level = min_level; level <= max_level; level++ ) {
				select_local_buffered_cells( level, proc, &num_level_cells, &level_cells );

				for ( i = 0; i < num_level_cells; i++ ) { 
					icell = level_cells[i];
					ipart = cell_particle_list[icell];
					while ( ipart != NULL_PARTICLE ) {
						if ( species == -1 || 
								( particle_species(particle_id[ipart]) == species
#ifdef STAR_PARTICLE_TYPES
								  && ( species != num_particle_species-1 || subspecies == -1 ||
									  star_particle_type[ipart] == subspecies ) 
#endif /* STAR_PARTICLE_TYPES */
								) ) {
							particle_list_to_send[proc][count++] = ipart;
						}
						ipart = particle_list_next[ipart];
					}
				}
				cart_free( level_cells );
			}

			cart_assert( count == num_parts_to_send[proc] );
		}
    }

	/* communicate */
    trade_particle_lists( num_parts_to_send, particle_list_to_send, -1, SAVE_PARTICLE_LISTS );	

	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_parts_to_send[proc] > 0 ) {
			cart_free( particle_list_to_send[proc] );
		}
	}

	particle_buffer_enabled = 1;
}

void destroy_particle_buffer() {
	int icell, i, level, num_level_cells, *level_cells;

	cart_assert( particle_buffer_enabled );

	/* loop over buffer cells and remove any particles within */
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			if ( cell_particle_list[icell] != NULL_PARTICLE ) {
				particle_list_free( cell_particle_list[icell] );
			}
			cell_particle_list[icell] = NULL_PARTICLE;
		}
		cart_free( level_cells );
	}
	particle_buffer_enabled = 0;
}

#endif /* PARTICLES */
