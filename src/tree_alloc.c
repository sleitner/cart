#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "hydro_tracer.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "sfc.h"
#include "tree.h"
#include "tree_linkedlist.h"

int next_free_oct = 0;
int free_oct_list = NULL_OCT;

/*******************************************************
 * cell_free
 ******************************************************/
void cell_free( int c ) 
/* purpose: frees a previously allocated cell and ensures
 * 	no octs have neighbors pointing to it.
 */
{
	int i;
	
	cart_assert( c >= 0 && c < num_cells );

	if ( cell_level(c) > min_level ) {
		if ( cell_is_local(c) ) {
			num_cells_per_level[ cell_level(c) ]--;
		} else {
			num_buffer_cells[ cell_level(c) ]--;
		}
	}

	if ( cell_is_refined(c) ) {
		for ( i = 0; i < num_children; i++ ) {
			cart_assert( cell_parent_cell( cell_child(c,i) ) == c );
			cart_assert( cell_is_local(c) == cell_is_local(cell_child(c,i) ) );
			cell_free( cell_child( c, i ) );
		}
		
		oct_free( cell_child_oct[c] );
		cell_child_oct[c] = UNREFINED_CELL;
	}

#ifdef PARTICLES
	cell_particle_list[c] = NULL_PARTICLE;
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	cell_tracer_list[c] = NULL_TRACER;
#endif /* HYDRO_TRACERS */
}

void cell_move( int icell_old, int icell_new ) {
	int i;
	int ioct;

	if ( cell_is_refined(icell_old) ) {
		ioct = cell_child_oct[icell_old];
		oct_parent_cell[ioct] = icell_new;
	}

        cell_child_oct[icell_new] = cell_child_oct[icell_old];
        cell_child_oct[icell_old] = UNREFINED_CELL;

	for ( i = 0; i < num_vars; i++ ) {
		cell_var(icell_new,i) = cell_var(icell_old,i);
	}

#ifdef PARTICLES
	cell_particle_list[icell_new] = cell_particle_list[icell_old];
	cell_particle_list[icell_old] = NULL_PARTICLE;
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	cell_tracer_list[icell_new] = cell_tracer_list[icell_old];
	cell_tracer_list[icell_old] = NULL_TRACER;
#endif /* HYDRO_TRACERS */
}

void oct_move( int oct_old, int oct_new ) {
	int i;

	cart_assert( oct_level[oct_old] != FREE_OCT_LEVEL );
	cart_assert( oct_level[oct_new] == FREE_OCT_LEVEL );
	cart_assert( oct_parent_cell[oct_old] >= 0 && oct_parent_cell[oct_old] < num_cells );
	cart_assert( cell_child_oct[ oct_parent_cell[oct_old] ] == oct_old );

	cell_child_oct[ oct_parent_cell[oct_old] ] = oct_new;
	oct_parent_cell[oct_new] = oct_parent_cell[oct_old];
	oct_level[oct_new] = oct_level[oct_old];
	oct_parent_root_sfc[oct_new] = oct_parent_root_sfc[oct_old];

	for ( i = 0; i < num_children; i++ ) {
		cell_move( oct_child( oct_old, i ), oct_child( oct_new, i ) );
	}

	for ( i = 0; i < num_neighbors; i++ ) {
		oct_neighbors[oct_new][i] = oct_neighbors[oct_old][i];
	}

	for ( i = 0; i < nDim; i++ ) {
		oct_pos[oct_new][i] = oct_pos[oct_old][i];
	}

	if ( root_cell_type(oct_parent_root_sfc[oct_new]) == CELL_TYPE_LOCAL ) {
		linked_list_insert( &local_oct_list[oct_level[oct_new]], oct_new );
	} else {
		linked_list_insert( &buffer_oct_list[oct_level[oct_new]], oct_new );
	}
}

/*******************************************************
 * oct_alloc
 *******************************************************/
int oct_alloc() 
/* purpose: takes a free oct from the octs array and
 * 	returns its index.
 *
 * returns: the index, or -1 if we have run out of free
 * 	octs.
 */
{
	int ioct;

	/* pull from the free_oct_list linked list 
	 * first, so we continually 
	 * refill the positions left by oct_free calls */
	if ( free_oct_list == NULL_OCT ) {
		if ( next_free_oct >= num_octs ) {
			cart_debug("oct_alloc: OUT OF OCTS!");
			return -1;
		}
	
		ioct = next_free_oct;
		next_free_oct++;
	} else {
		ioct = free_oct_list;
		free_oct_list = oct_next[ioct];

		if ( free_oct_list != NULL_OCT ) {
			oct_prev[free_oct_list] = NULL_OCT;
		}
	}

	cart_assert( ioct >= 0 && ioct < num_octs );
	cart_assert( oct_level[ioct] == FREE_OCT_LEVEL );
	cart_assert( ioct >= cell_parent_oct( num_cells_per_level[min_level] + num_buffer_cells[min_level] ) + 1 ); 
	return ioct;
}

/*******************************************************
 * oct_free
 *******************************************************/
void oct_free( int ioct ) 
/* purpose: returns a previously allocated oct to the
 * 	pool of free octs.  Sets its level to FREE_OCT_LEVEL
 * 	so it is skipped in the level iterators.  Also 
 * 	removes it from the proper linked list.
 */
{
	int i;
	int level;

	cart_assert( ioct >= 0 && ioct < num_octs );

	level = oct_level[ioct];
	cart_assert( level != FREE_OCT_LEVEL );

	if ( root_cell_type(oct_parent_root_sfc[ioct]) == CELL_TYPE_LOCAL ) {
		linked_list_remove( &local_oct_list[level], ioct );
	} else {
		linked_list_remove( &buffer_oct_list[level], ioct );
	}

	oct_parent_root_sfc[ioct] = NULL_OCT;
	oct_parent_cell[ioct] = NULL_OCT;
	oct_level[ioct] = FREE_OCT_LEVEL;

	for ( i = 0; i < num_neighbors; i++ ) {
		oct_neighbors[ioct][i] = NULL_OCT;
	}

	/* add new oct to free oct list */
	oct_next[ioct] = free_oct_list;
	oct_prev[ioct] = NULL_OCT;

	if ( free_oct_list != NULL_OCT ) {
		oct_prev[free_oct_list] = ioct;
	}
	free_oct_list = ioct;
}

 /*******************************************************
 * split_cell
 *******************************************************/
int split_cell( int icell ) 
/* purpose: splits the given cell into num_children and
 * 	handles tree variables (does NOT handle interpolation
 * 	of actual cell variables.
 *
 * returns: 	 0	if operation is successful
 * 		-1	if cell is already split
 * 		-2	if there are no free octs remaining
 * NO LONGER USED: we need to use this function on buffer
 * 	cells which may violate +/- 1 refinement criteria
 * 	since they are pruned, so the refinement check
 * 	is moved to split
 *  		-4	if the split would violate the 
 * 				+/- 1 refinement criteria
 */
{
	int i;
	int oct_ptr;
	int type;

	cart_assert ( icell >= 0 && icell < num_cells );
	cart_assert ( cell_level( icell ) < max_level );

	/* make sure cell isn't already split */
	if ( cell_is_refined(icell) ) {
		cart_debug("split_cell(%u) error: cell already refined!", icell );
		return -1;
	}

	/* allocate a new oct */
	oct_ptr = oct_alloc();

	/* make sure we haven't run out of octs */
	if ( oct_ptr == NULL_OCT ) {
		return -2;
	}

	oct_level[oct_ptr] = cell_level(icell)+1;
	cart_assert( oct_level[oct_ptr] > min_level && oct_level[oct_ptr] <= max_level );

	/* add oct to appropriate linked list */
	type = root_cell_type( cell_parent_root_sfc(icell) );
	if ( type == CELL_TYPE_LOCAL ) {
		linked_list_insert( &local_oct_list[oct_level[oct_ptr]], oct_ptr );
		num_cells_per_level[oct_level[oct_ptr]] += num_children;
	} else if ( type == CELL_TYPE_BUFFER ) {
		linked_list_insert( &buffer_oct_list[oct_level[oct_ptr]], oct_ptr );
		num_buffer_cells[oct_level[oct_ptr]] += num_children;
	} else {
		cart_error("Bad cell type in split_cell!");
	}

	/* set up oct neighbors */
	cell_all_neighbors( icell, oct_neighbors[oct_ptr] );

	/* set up oct position */
	cell_position( icell, oct_pos[oct_ptr] );

	/* set up new oct */
	oct_parent_cell[oct_ptr] = icell;
	oct_parent_root_sfc[oct_ptr] = cell_parent_root_sfc(icell);

	/* set all newly allocated children to leaves */
	for ( i = 0; i < num_children; i++ ) {
		cell_child_oct[ oct_child( oct_ptr, i ) ] = UNREFINED_CELL;
	}

	/* set cell child pointer */
	cell_child_oct[icell] = oct_ptr;

#ifdef PARTICLES
	for ( i = 0; i < num_children; i++ ) {
		cell_particle_list[ oct_child( oct_ptr, i ) ] = NULL_PARTICLE;
	}

	split_particle_list( icell );
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	for ( i = 0; i < num_children; i++ ) {
		cell_tracer_list[ oct_child( oct_ptr, i ) ] = NULL_TRACER;
	}

	split_tracer_list( icell );
#endif /* HYDRO_TRACERS */

	return 0;
}

int join_cell( int icell ) {
	int i;

	cart_assert( icell >= 0 && icell < num_cells );

	if ( cell_is_leaf(icell) ) {
		cart_debug("ERROR: called join_cell(%u) on a leaf", icell );
		return -1;
	}

	cart_assert( cell_child_oct[icell] != NULL_OCT );
	cart_assert( oct_level[cell_child_oct[icell]] != FREE_OCT_LEVEL );

#ifdef PARTICLES
	join_particle_list(icell);
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	join_tracer_list(icell);
#endif /* HYDRO_TRACERS */

	for ( i = 0; i < num_children; i++ ) {
		cell_free( cell_child( icell, i ) );
	}

	oct_free( cell_child_oct[icell] );
	cell_child_oct[icell] = UNREFINED_CELL; 

	return 0;
}
