#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "pack.h"
#include "parallel.h"
#include "sfc.h"
#include "timing.h"
#include "tree.h"


pack *pack_init( int cell_type ) {
	int proc;
	int level;
	pack *ret;

	ret = cart_alloc(pack, 1 );
	for ( proc = 0; proc < num_procs; proc++ ) {
		ret->tree_list[proc] = skiplist_init();
		ret->num_sending_cells_total[proc] = 0;
		ret->num_receiving_cells_total[proc] = 0;
		for ( level = min_level; level <= max_level; level++ ) {
			ret->num_sending_cells[proc][level-min_level] = 0;
			ret->num_receiving_cells[proc][level-min_level] = 0;
		}
	}

	ret->cell_type = cell_type;

	return ret;
}

void pack_destroy( pack *p ) {
	int proc;

	cart_assert( p != NULL );

	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( p->num_sending_cells[proc][min_level-min_level] > 0 ) {
			cart_free( p->root_cells[proc] );
			cart_free( p->cell_refined[proc] );
			cart_free( p->cell_vars[proc] );
		}
	}

	cart_free( p );	
}

int prune_proc;
int cell_count_with_pruning( int c ) {
	int i;
	int count;

        cart_assert( c >= 0 && c < num_cells );
	cart_assert( cell_is_local(c) );

	count = 1;

	if ( cell_is_refined(c) && !cell_can_prune( c, prune_proc ) ) {
		for ( i = 0; i < num_children; i++ ) {
			cart_assert( cell_is_local( cell_child(c,i) ) );
			cart_assert( cell_level( cell_child(c,i) ) > min_level );
			count += cell_count_with_pruning(cell_child(c,i));
		}
	}

	return count;
}

void pack_add_root_trees( pack *p, int *new_proc_sfc_index, int sfc1, int sfc2 ) {
	int j;
	int proc;
	int block_count;
	int a, b, c;
	int sfc = sfc1;

	cart_assert( sfc1 <= sfc2 );

	/* find processor for sfc1 */
	if ( sfc < new_proc_sfc_index[local_proc_id] ) {
		a = 0;
		b = local_proc_id-1;
	} else if ( sfc >= new_proc_sfc_index[local_proc_id+1] ) {
		a = local_proc_id + 1;
		b = num_procs-1;
	} else {
		cart_error("pack_add_root_trees called with new local sfc index!");
	}

	/* do binary search between procs a & b */
	while ( a != b ) {
		c = ( a + b + 1) / 2;

		if ( sfc < new_proc_sfc_index[c] ) {
			b = c-1;
		} else {
			a = c;
		}
	}

	proc = a;

	while ( sfc < sfc2 ) {
		block_count = MIN( sfc2, new_proc_sfc_index[proc+1] ) - sfc;
		
		for ( j = sfc; j < sfc+block_count; j++ ) {
			pack_add_root_tree(p, proc, j);
		}

		sfc += block_count;
		proc++;
	}

	cart_assert( sfc == sfc2 );
}


void pack_add_root_tree( pack *p, int proc, int sfc ) {
	cart_assert( p != NULL );
	cart_assert( proc >= 0 && proc < num_procs && proc != local_proc_id );
	cart_assert( sfc >= 0 && sfc < max_sfc_index );

	if ( skiplist_insert( p->tree_list[proc], sfc ) ) {
		p->num_sending_cells[proc][min_level-min_level]++;

		if ( p->cell_type == CELL_TYPE_BUFFER ) {
			/* prune cells if we're building a cell buffer */
			prune_proc = proc;
			cart_assert( root_cell_is_local(sfc) );
			cart_assert( cell_is_local( root_cell_location(sfc) ) );
			p->num_sending_cells_total[proc] += cell_count_with_pruning( root_cell_location(sfc) );
		} else {
			/* otherwise we're load balancing and send the entire tree */
			p->num_sending_cells_total[proc] += tree_traversal( root_cell_location( sfc ), cell_count_cells );
		}
	}
}

void pack_apply( pack *p ) {
	int i, j;
	int proc;
	int level;
	int icell, ioct;
	int child;
	int num_vars_packed;
	int level_count;
	int *root_cells;
	float *cell_packed_vars;
	int *cell_packed_child;
	int num_cells_packed;
	int old_offset;

	cart_assert( p != NULL );

	skiplist_destroy( p->tree_list[local_proc_id] );

	for ( proc = 0; proc < num_procs; proc++ ) {
		prune_proc = proc;

		if ( p->num_sending_cells[proc][min_level-min_level] > 0 ) {
			root_cells = cart_alloc(int, p->num_sending_cells[proc][min_level-min_level] );

			/* pack root sfc's into communication array */
			i = 0;
			skiplist_iterate( p->tree_list[proc] );
			while ( skiplist_next( p->tree_list[proc], &root_cells[i] ) ) i++;
			skiplist_destroy( p->tree_list[proc] );

			for ( i = 0; i < p->num_sending_cells[proc][min_level-min_level]; i++ ) {
				cart_assert( root_cells[i] >= 0 && root_cells[i] < max_sfc_index );
				cart_assert( root_cells[i] < proc_sfc_index[local_proc_id+1] );
			}

			p->root_cells[proc] = root_cells;

			/* allocate space for cell tree information */
			cell_packed_vars = cart_alloc(float, p->num_sending_cells_total[proc] * num_vars );
			cell_packed_child = cart_alloc(int, p->num_sending_cells_total[proc] );

			/* pack root cells */
			num_vars_packed = 0;
			for ( i = 0; i < p->num_sending_cells[proc][min_level-min_level]; i++ ) {
				/* pack tree into communication arrays */
				icell = root_cell_location(root_cells[i]);
				cart_assert( icell >= 0 && icell < num_cells_per_level[min_level] );

				/* pack root cell vars */
				for ( j = 0; j < num_vars; j++ ) {
					cell_packed_vars[num_vars_packed++] = cell_var(icell,j);
				}
			}

			/* pack root cell child ptr */
			if ( p->cell_type == CELL_TYPE_BUFFER ) {
				for ( i = 0; i < p->num_sending_cells[proc][min_level-min_level]; i++ ) {
					icell = root_cell_location(root_cells[i]);
					if ( cell_is_refined(icell) && !cell_can_prune(icell,proc) ) {
						cell_packed_child[i] = cell_child_oct[icell];
					} else {
						cell_packed_child[i] = NULL_OCT;
					}
				}
			} else {
				for ( i = 0; i < p->num_sending_cells[proc][min_level-min_level]; i++ ) {
					icell = root_cell_location(root_cells[i]);
					cell_packed_child[i] = cell_child_oct[icell];
				}
			}

			old_offset = 0;
			num_cells_packed = p->num_sending_cells[proc][min_level-min_level];

			/* pack each level in turn */
			for ( level = min_level+1; level <= max_level; level++ ) {
				level_count = 0;

				for ( i = old_offset; i < old_offset + p->num_sending_cells[proc][level-1-min_level]; i++ ) {
					ioct = cell_packed_child[i];

					if ( ioct != NULL_OCT ) {
						cart_assert( oct_level[ioct] != FREE_OCT_LEVEL );

						for ( child = 0; child < num_children; child++ ) {
							cart_assert( num_cells_packed < p->num_sending_cells_total[proc] );

							icell = oct_child( ioct, child );

							/* pack cell vars */
							for ( j = 0; j < num_vars; j++ ) { 
								cell_packed_vars[num_vars_packed++] = cell_var( icell, j );
							}

							/* pack child ptr */
							/*
							if ( cell_is_refined(icell) &&
									( p->cell_type != CELL_TYPE_BUFFER || 
									  !cell_can_prune( icell, proc) ) ) {
								cell_packed_child[num_cells_packed++] = 
										cell_child_oct[icell];
							} else {
								cell_packed_child[num_cells_packed++] = NULL_OCT;
							}
							*/
	
							level_count++;
						}
					}
				}

				/* pack child ptr's */
				if ( p->cell_type == CELL_TYPE_BUFFER ) {
					for ( i = old_offset; i < old_offset + p->num_sending_cells[proc][level-1-min_level]; i++ ) {
						ioct = cell_packed_child[i];

						if ( ioct != NULL_OCT ) {
							for ( child = 0; child < num_children; child++ ) {
								icell = oct_child( ioct, child );
								if ( cell_is_refined(icell) && !cell_can_prune(icell,proc) ) {
									cell_packed_child[num_cells_packed++] = 
										cell_child_oct[icell];
								} else {
									cell_packed_child[num_cells_packed++] = NULL_OCT;
								}
							}
						}
					}
				} else {
					for ( i = old_offset; i < old_offset + p->num_sending_cells[proc][level-1-min_level]; i++ ) {
						ioct = cell_packed_child[i];

						if ( ioct != NULL_OCT ) {
							for ( child = 0; child < num_children; child++ ) {
								icell = oct_child( ioct, child );
								cell_packed_child[num_cells_packed++] = cell_child_oct[icell];
							}
						}
					}
				}

				p->num_sending_cells[proc][level-min_level] = level_count;
				old_offset = num_cells_packed - level_count;

				cart_assert( old_offset >= 0 );
			}

			cart_assert( num_cells_packed == p->num_sending_cells_total[proc] );

			p->cell_refined[proc] = cell_packed_child;
			p->cell_vars[proc] = cell_packed_vars;
		}
	}

	/* do global communication for how many cells to expect */
	MPI_Alltoall( p->num_sending_cells, max_level-min_level+1, MPI_INT,
		p->num_receiving_cells, max_level-min_level+1, MPI_INT,
		mpi.comm.run );
	MPI_Alltoall( p->num_sending_cells_total, 1, MPI_INT,
		p->num_receiving_cells_total, 1, MPI_INT, mpi.comm.run );
}

void pack_communicate( pack *p ) {
	int i, j;
	int ret;
	int level;
	int proc;
	int icell, child;
	int next_level_count;
	int num_cells_unpacked;
	int num_cells_split;
	int num_vars_unpacked;
	int num_octs_unpacked;
	int cell_count;
	int num_recv_octs;
	int *root_cells[MAX_PROCS];
	int *cell_refined[MAX_PROCS];
	float *cell_recv_vars[MAX_PROCS];
	int *oct_indices[MAX_PROCS][2];
	int num_oct_indices[MAX_PROCS];
	int level_offset[MAX_PROCS];
	int *next_level_octs, *current_level_octs;
	int num_sends, num_receives;
	int num_send_requests, num_recv_requests;
	MPI_Request *sends;
	MPI_Request *receives;

#ifdef MPI_MAX_MESSAGE_SIZE
	int size;
#endif

	if ( p->cell_type == CELL_TYPE_BUFFER ) {
		for ( proc = 0; proc < num_procs; proc++ ) {
			if ( p->num_receiving_cells[proc][min_level-min_level] > 0 ) {
				num_recv_octs = ( p->num_receiving_cells_total[proc] -
					p->num_receiving_cells[proc][min_level-min_level] ) / num_children;

				oct_indices[proc][0] = cart_alloc(int, num_recv_octs );
				oct_indices[proc][1] = cart_alloc(int, num_recv_octs );
			}

			num_oct_indices[proc] = 0;
		}
	}

	for ( proc = 0; proc < num_procs; proc++ ) {
		level_offset[proc] = 0;
	}

	for ( level = min_level; level <= max_level; level++ ) {
		next_level_count = 0;
		num_receives = 0;
		num_sends = 0;

		/* determine number of pages if we have a limited send size */
#ifdef MPI_MAX_MESSAGE_SIZE
		num_send_requests = num_recv_requests = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			cell_count = p->num_receiving_cells[proc][level-min_level];
			if ( cell_count > 0 ) {
				if ( level == min_level ) {
					num_recv_requests += (sizeof(int)*cell_count-1)/MPI_MAX_MESSAGE_SIZE+1;
				}

				num_recv_requests += (sizeof(int)*cell_count-1)/MPI_MAX_MESSAGE_SIZE+1;
				num_recv_requests += (sizeof(float)*cell_count*num_vars-1)/MPI_MAX_MESSAGE_SIZE+1;
			}

			cell_count = p->num_sending_cells[proc][level-min_level];
			if ( cell_count > 0 ) {
				if ( level == min_level ) {
					num_send_requests += (sizeof(int)*cell_count-1)/MPI_MAX_MESSAGE_SIZE+1;
				}

				num_send_requests += (sizeof(int)*cell_count-1)/MPI_MAX_MESSAGE_SIZE+1;
				num_send_requests += (sizeof(float)*cell_count*num_vars-1)/MPI_MAX_MESSAGE_SIZE+1;
			}
		}
#else
		num_send_requests = num_recv_requests = 3*num_procs;
#endif

		sends = cart_alloc( MPI_Request, num_send_requests );
		receives = cart_alloc( MPI_Request, num_recv_requests );

		for ( proc = 0; proc < num_procs; proc++ ) {
			if ( level < max_level ) {
				next_level_count += p->num_receiving_cells[proc][level+1-min_level] / num_children;
			}

			cell_count = p->num_receiving_cells[proc][level-min_level];
			if ( cell_count > 0 ) {
				cell_refined[proc] = cart_alloc(int, cell_count );
				cell_recv_vars[proc] = cart_alloc(float, num_vars * cell_count );

				if ( level == min_level ) {
					root_cells[proc] = cart_alloc(int, cell_count );

#ifdef MPI_MAX_MESSAGE_SIZE
					i = 0;
					do {
						cart_assert( num_receives < num_recv_requests );
						size = MIN( MPI_MAX_MESSAGE_SIZE/sizeof(int), cell_count-i );
						MPI_Irecv( &root_cells[proc][i], size, MPI_INT, proc, i,
								mpi.comm.run, &receives[num_receives++] );
						i += size;
					} while ( i < cell_count );
#else 
					MPI_Irecv( root_cells[proc], cell_count, MPI_INT, proc, cell_count, 
						mpi.comm.run, &receives[num_receives++] );
#endif
				}

#ifdef MPI_MAX_MESSAGE_SIZE
				i = 0;
				do {
					cart_assert( num_receives < num_recv_requests );
					size = MIN( MPI_MAX_MESSAGE_SIZE/sizeof(int), cell_count-i );
					MPI_Irecv( &cell_refined[proc][i], size, MPI_INT, proc, i+cell_count,
							mpi.comm.run, &receives[num_receives++] );
					i += size;
				} while ( i < cell_count );

				i = 0;
				do {
					cart_assert( num_receives < num_recv_requests );
					size = MIN( MPI_MAX_MESSAGE_SIZE/sizeof(float), cell_count*num_vars-i );
					MPI_Irecv( &cell_recv_vars[proc][i], size, MPI_FLOAT, proc,
							i, mpi.comm.run, &receives[num_receives++] );
					i += size;
				} while ( i < cell_count*num_vars );
#else
				MPI_Irecv( cell_refined[proc], cell_count, MPI_INT, proc, level,
						mpi.comm.run, &receives[num_receives++] );
				MPI_Irecv( cell_recv_vars[proc], num_vars*cell_count, MPI_FLOAT, proc,
					level, mpi.comm.run, &receives[num_receives++] );
#endif
			}

			cell_count = p->num_sending_cells[proc][level-min_level];
			if ( cell_count > 0 ) {
				if ( level == min_level ) {
#ifdef MPI_MAX_MESSAGE_SIZE
					i = 0;
					do {
						cart_assert( num_sends < num_send_requests );
						size = MIN( MPI_MAX_MESSAGE_SIZE/sizeof(int), cell_count-i );
						MPI_Isend( &p->root_cells[proc][i], size, MPI_INT, proc, i,
								mpi.comm.run, &sends[num_sends++] );
						i += size;
					} while ( i < cell_count );
#else 
					MPI_Isend( p->root_cells[proc], cell_count, MPI_INT, proc, cell_count,
							mpi.comm.run, &sends[num_sends++] );
#endif
				}

#ifdef MPI_MAX_MESSAGE_SIZE
				i = 0;
				do {
					cart_assert( num_sends < num_send_requests );
					size = MIN( MPI_MAX_MESSAGE_SIZE/sizeof(int), cell_count-i );
					MPI_Isend( &p->cell_refined[proc][level_offset[proc]+i], size, MPI_INT, 
							proc, i+cell_count, mpi.comm.run, &sends[num_sends++] );
					i += size;
				} while ( i < cell_count );

				i = 0;
				do {
					cart_assert( num_sends < num_send_requests );
					size = MIN( MPI_MAX_MESSAGE_SIZE/sizeof(float), cell_count*num_vars-i );
					MPI_Isend( &p->cell_vars[proc][num_vars*level_offset[proc]+i], size, MPI_FLOAT,
							proc, i, mpi.comm.run, &sends[num_sends++] );
					i += size;
				} while ( i < cell_count*num_vars );

#else 
				MPI_Isend( &p->cell_refined[proc][level_offset[proc]], cell_count, MPI_INT, proc, 
						level, mpi.comm.run, &sends[num_sends++] );
				MPI_Isend( &p->cell_vars[proc][num_vars*level_offset[proc]], num_vars*cell_count, 
						MPI_FLOAT, proc, level, mpi.comm.run, &sends[num_sends++] );
#endif

				level_offset[proc] += cell_count;
			}
		}
	
		MPI_Waitall( num_receives, receives, MPI_STATUSES_IGNORE );

		if ( p->cell_type == CELL_TYPE_BUFFER ) {
			if ( level == min_level ) {
				for ( proc = 0; proc < num_procs; proc++ ) {
					for ( i = 0; i < p->num_receiving_cells[proc][min_level-min_level]; i++ ) {
						local_buffers[min_level][proc][i] = 
							cell_buffer_local_index( root_cells[proc][i] );
					}
				}
			} 

			if ( level < max_level ) {
				for ( proc = 0; proc < num_procs; proc++ ) {
					num_local_buffers[level+1][proc] = 0;
				}
			}
		}

		/* unpack level */
		num_cells_split = 0;
		next_level_octs = cart_alloc(int, next_level_count );
		num_octs_unpacked = 0;
		child = 0;

		for ( proc = 0; proc < num_procs; proc++ ) {
			cell_count = p->num_receiving_cells[proc][level-min_level];

			if ( cell_count > 0 ) {
				cart_assert( proc != local_proc_id );

				num_cells_unpacked = 0;
				num_vars_unpacked = 0;

				while ( num_cells_unpacked < cell_count ) {
					if ( level == min_level ) {
						/* root cells are pre-allocated */
						if ( p->cell_type == CELL_TYPE_BUFFER ) {
							icell = local_buffers[min_level][proc][num_cells_unpacked];
						} else {
							icell = root_cell_location( root_cells[proc][num_cells_unpacked] );
						}
					} else {
						if ( child == num_children ) {
							num_octs_unpacked++;
							child = 0;
						}

						icell = oct_child( current_level_octs[num_octs_unpacked], child );
						child++;
					}

					cart_assert( icell >= 0 && icell < num_cells );
					cart_assert( cell_is_leaf(icell) );

					for ( j = 0; j < num_vars; j++ ) {
						cell_var(icell,j) = cell_recv_vars[proc][num_vars_unpacked++];
					}

					if ( cell_refined[proc][num_cells_unpacked] != NULL_OCT ) {
						ret = split_cell( icell );
						if ( ret ) {
							if ( level == min_level ) {
								cart_debug("root_cells[%u][%u] = %d", proc, 
									num_cells_unpacked, 
									root_cells[proc][num_cells_unpacked] );
							}
							cart_error("Error in split in pack_communicate, ret = %d, icell = %u, level = %u, proc = %u, oct = %d", 
								ret, icell, level, proc, cell_child_oct[icell] );
						}

						cart_assert( cell_child_oct[icell] != NULL_OCT );
						cart_assert( oct_level[ cell_child_oct[icell] ] == level+1 );
						cart_assert( oct_parent_cell[cell_child_oct[icell]] == icell );
						cart_assert( num_cells_split < next_level_count );
						next_level_octs[num_cells_split++] = cell_child_oct[icell];

						if ( p->cell_type == CELL_TYPE_BUFFER ) {
							cart_assert( cell_refined[proc][num_cells_unpacked] >= 0 );
							cart_assert( cell_child_oct[icell] >= 0 );

							oct_indices[proc][0][num_oct_indices[proc]] = 
								cell_refined[proc][num_cells_unpacked];
							oct_indices[proc][1][num_oct_indices[proc]] = 
								cell_child_oct[icell];
							num_oct_indices[proc]++;

							local_buffers[level+1][proc][num_local_buffers[level+1][proc]++] = 
								cell_child_oct[icell];
						}
					}

					num_cells_unpacked++;
				}

				/* make sure we unpacked a full set of oct children */
				cart_assert( level == min_level || child == num_children );
	
				if ( level == min_level ) {
					cart_free( root_cells[proc] );
				}

				cart_free( cell_recv_vars[proc] );
				cart_free( cell_refined[proc] );
			}
		}

		cart_assert( num_cells_split == next_level_count );

		if ( level > min_level ) {
			cart_free( current_level_octs );
		}

		current_level_octs = next_level_octs;

		MPI_Waitall( num_sends, sends, MPI_STATUSES_IGNORE );

		cart_free( sends );
		cart_free( receives );
	}

	cart_free( current_level_octs );

	if ( p->cell_type == CELL_TYPE_BUFFER ) {
		for ( proc = 0; proc < num_procs; proc++ ) {
			if ( num_oct_indices[proc] > 0 ) {
				oct_hash_add_list( buffer_oct_hash[proc],  num_oct_indices[proc], 
					oct_indices[proc][0], oct_indices[proc][1] );
				oct_hash_add_list( buffer_oct_reverse_hash[proc], num_oct_indices[proc], 
					oct_indices[proc][1], oct_indices[proc][0] );
        		}

			if ( p->num_receiving_cells[proc][min_level-min_level] > 0 ) {
				cart_free( oct_indices[proc][0] );
				cart_free( oct_indices[proc][1] );
			}
		}
	}
}
