#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"
#include "tree.h"
#include "particle.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "timing.h"
#include "auxiliary.h"

/* global tree variables */
float cell_vars[num_cells][num_vars];
int cell_child_oct[num_cells];

int oct_parent_cell[num_octs];
int oct_level[num_octs];
int oct_neighbors[num_octs][num_neighbors];
int oct_parent_root_sfc[num_octs];
int oct_next[num_octs];
int oct_prev[num_octs];
float oct_pos[num_octs][nDim];

int all_vars[num_vars];
#ifdef HYDRO
int all_hydro_vars[num_hydro_vars];
#endif

/* level linked list variables */
int num_cells_per_level[max_level-min_level+1];
int local_oct_list[max_level-min_level+1];

/* tables of stats for each level */
float cell_size[max_level-min_level+1];
float cell_size_inverse[max_level-min_level+1];
float cell_volume[max_level-min_level+1];
float cell_volume_inverse[max_level-min_level+1];

/*******************************************************
 * init_tree
 *******************************************************/
void init_tree() 
/* purpose: initializes the level tables (assumes proc_sfc_index set)
 *   needs to be called before a restart
 */
{
	int i;

	for ( i = 0; i < num_vars; i++ ) {
		all_vars[i] = i;
	}

#ifdef HYDRO
	for ( i = 0; i < num_hydro_vars; i++ ) {
		all_hydro_vars[i] = HVAR_GAS_DENSITY + i;
	}
#endif /* HYDRO */

	for ( i = min_level; i <= max_level; i++ ) {
		cell_size_inverse[i-min_level] = (float)(1<<i);
		cell_size[i-min_level] = 1.0 / cell_size_inverse[i-min_level];
		cell_volume_inverse[i-min_level] = pow(2.,3*i);
		cell_volume[i-min_level] = 1.0 / cell_volume_inverse[i-min_level];

		num_cells_per_level[i] = 0;
		local_oct_list[i] = NULL_OCT;
	}

#pragma omp parallel for default(none), private(i), shared(cell_child_oct)
	for ( i = 0; i < num_cells; i++ ) {
		cell_child_oct[i] = UNREFINED_CELL;
	}

	/* this now starts at first_oct due to root cells */
#pragma omp parallel for default(none), private(i), shared(oct_parent_cell,oct_parent_root_sfc,oct_level,oct_next,oct_prev)
	for ( i = 0; i < num_octs; i++ ) {
		oct_parent_cell[i] = -1;
		oct_parent_root_sfc[i] = -1;
		oct_level[i] = FREE_OCT_LEVEL;
		oct_next[i] = NULL_OCT;
		oct_prev[i] = NULL_OCT;
	}

	/* start with only root cells, more cells are created via refinement routines */
	num_cells_per_level[min_level] = proc_sfc_index[local_proc_id+1] - proc_sfc_index[local_proc_id];
	build_root_cell_buffer();

	next_free_oct = cell_parent_oct( num_cells_per_level[min_level] + num_buffer_cells[min_level] ) + 1;
	free_oct_list = NULL_OCT;
}

int max_level_now_global(MPI_Comm local_comm) {
	int level, max;

	level = max_level_local();

	start_time( COMMUNICATION_TIMER );
	start_time( MAX_LEVEL_TIMER );
	MPI_Allreduce( &level, &max, 1, MPI_INT, MPI_MAX, local_comm );
	end_time( MAX_LEVEL_TIMER );
	end_time( COMMUNICATION_TIMER );

	return max;
}

int max_level_now() {
	return max( max_level_local(), max_level_buffer() );
}

int max_level_local() {
	int i;
	int level;
	
	level = max_level;
	for ( i = min_level+1; i <= max_level; i++ ) {
		if ( num_cells_per_level[i] == 0 ) {
			level = i-1;
			break;
		}
	}

	return level;
}

int max_level_buffer() {
	int i;
	int level;

	level = max_level;
	for ( i = min_level+1; i <= max_level; i++ ) {
		if ( num_buffer_cells[i] == 0 ) {
			level = i-1;
			break;
		}
	}

	return level;
}

/*******************************************************
 * repair_neighbors
 *******************************************************/
void repair_neighbors() 
/* purpose: repairs all null neighbor pointers.  Works
 * 	from min_level -> max_level so the previous
 * 	level's neighbor pointers are always correct,
 * 	so they may be used to calculate sucessive levels
 */
{
	int i;
	int level;
	int ioct, icell;
	int num_level_cells;
	int *level_cells;

	/* repair neighbors level by level so each succeeding level can use
	 * previous for calculating neighbors */
	for ( level = min_level; level < max_level; level++ ) {
		select_level( level, CELL_TYPE_ANY, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,ioct), shared(num_level_cells,level_cells,cell_child_oct,oct_neighbors)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_refined(icell) ) {
				ioct = cell_child_oct[icell];
				cell_all_neighbors( icell, oct_neighbors[ioct] );
			}
		}
		cart_free( level_cells );
	}
}

/*******************************************************
 * root_cell_type
 *******************************************************/
int root_cell_type( int sfc ) 
/* purpose: identifies what type of root cell corresponds
 *  to the given index
 *
 *  returns: 1 if cell is local, 2 if cell is buffer, 
 *  	0 if cell is non-local
 */
{
	cart_assert( sfc >= 0 && sfc < max_sfc_index );
	
	if ( root_cell_is_local( sfc ) ) {
		return CELL_TYPE_LOCAL;
	} else if ( cell_buffer_exists(sfc) ) {
		return CELL_TYPE_BUFFER;
	} else {
		return CELL_TYPE_NONLOCAL;
	}	
}

int cell_is_local( int icell ) {
	int ioct;

	cart_assert( icell >= 0 && icell < num_cells );

	if ( cell_is_root_cell(icell) ) {
		return ( icell < num_cells_per_level[min_level] );
	} else {
		ioct = cell_parent_oct(icell);
		return root_cell_is_local( oct_parent_root_sfc[ioct] );
	}
}

/*******************************************************
 * root_cell_location
 *******************************************************/
int root_cell_location( int sfc )
/* purpose: determines where root cell is from its sfc index
 *
 * returns: index if cell is local or buffer, otherwise -1 
 * 	(passed by cell_buffer_local_index)
 */
{
	cart_assert( sfc >= 0 && sfc < max_sfc_index );
	
	if ( root_cell_is_local(sfc) ) { 
		return (sfc - proc_sfc_index[local_proc_id]);
	} else {
		return cell_buffer_local_index(sfc);
	}
}

/*******************************************************
 * cell_parent_cell
 *******************************************************/
int cell_parent_cell( int c ) 
/* equivalent to iPr in the fortran version of the code
 * 
 * purpose: finds parent for cell c
 * input: c -> integer offset of cell (0 <= c < num_cells)
 *        c <= num_root_cells -> no parent
 *        
 * returns: integer offset of parent cell or -1 if root cell or invalid
 */ 
{
	/* catch array bounds exception */
	cart_assert ( c >= 0 && c < num_cells );

	if ( cell_is_root_cell(c) ) {
		return -1;
	} else {
		return oct_parent_cell[cell_parent_oct(c)];
	}
}

/*******************************************************
 * cell_parent_root_cell
 *******************************************************/
int cell_parent_root_cell( int c )
/* purpose: determines the root cell which contains cell c
 *
 * Changed 11/8/02 to use oct field parent_root_cell
 * 
 * returns: c if is root cell, otherwise root cell parent
 */
{
	cart_assert ( c >= 0 && c < num_cells );

	if ( cell_is_root_cell(c) ) {
		return c;
	} else {
		return oct_parent_root_cell( cell_parent_oct(c) );
	}
}

/*******************************************************
 * oct_parent_root_cell
 *******************************************************/
int oct_parent_root_cell( int ioct ) 
/* purpose: determines the root cell which contains the
 * 	given oct.  Simply shortcut for running
 * 	cell_parent_root_cell on the oct's parent_cell
 *
 * 	Changed 11/8/02 to use oct field parent_root_cell
 */
{
	cart_assert ( ioct >= 0 && ioct < num_octs );

	return root_cell_location( oct_parent_root_sfc[ioct] );
}

int root_cell_sfc_index( int icell ) {
	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( cell_is_root_cell(icell) );
	cart_assert( root_buffer_enabled );

	if ( icell < num_cells_per_level[min_level] ) {
		return proc_sfc_index[local_proc_id] + icell;
	} else {
		return buffer_cell_sfc_index[ icell - num_cells_per_level[min_level] ];
	}
}

int cell_parent_root_sfc( int c ) {

	cart_assert( c >= 0 && c < num_cells ); 

	if ( cell_is_root_cell(c) ) {
		return root_cell_sfc_index(c);
	} else {
		return oct_parent_root_sfc[ cell_parent_oct(c) ];	
	}
}

/*******************************************************
 * cell_level
 *******************************************************/
int cell_level( int c ) 
/* equivalent to iLv in fortran version
 *
 * purpose: finds level for cell c
 *
 * returns: minlevel for root level cells
 * 		-1 for invalid cells
 * 		minlevel+1->maxlevel for normal cells
 */
{
	cart_assert ( c >= 0 && c < num_cells );

	if ( cell_is_root_cell(c) ) {
		return min_level;
	} else {
		return oct_level[cell_parent_oct(c)];
	}
}

/*******************************************************
 * cell_position
 *******************************************************/
void cell_position( int c, float position[nDim] ) {
/* equivalent to iPs in fortran version 
 *
 * purpose: finds coordinates of center of cell c
 */
	int i;
	int coords[nDim];
	int level, child;
	int parent;

	cart_assert ( c >= 0 && c < num_cells );

	if ( cell_is_root_cell(c) ) {
		/* convert sfc index of cell to 3-d coordinates */
		sfc_coords( cell_parent_root_sfc(c), coords);

		/* center of cell is actually cell_center_offset from
		 * its index (otherwise 0,0,0 cell straddles into negative
		 * coordinates) */
		for ( i = 0; i < nDim; i++ ) {
			position[i] = (float)coords[i] + cell_center_offset;
		}
	} else {
		parent = cell_parent_oct(c);
		level = oct_level[parent];
		child = cell_child_number(c); 

		/* use delta to compute offset from parent's center */
		for ( i = 0; i < nDim; i++ ) {
			position[i] = oct_pos[parent][i] + cell_size[level] * cell_delta[child][i];
		}
	}
}

void cell_position_double( int c, double position[nDim] ) {
/* equivalent to iPs2 in fortran version 
 *
 * purpose: finds coordinates of center of cell c
 */
	int i;
	int coords[nDim];
	int level, child;
	int parent;

	cart_assert ( c >= 0 && c < num_cells );

	if ( cell_is_root_cell(c) ) {
		/* convert sfc index of cell to 3-d coordinates */
		sfc_coords( cell_parent_root_sfc(c), coords);

		/* center of cell is actually cell_center_offset from
		 * its index (otherwise 0,0,0 cell straddles into negative
		 * coordinates) */
		for ( i = 0; i < nDim; i++ ) {
			position[i] = (double)coords[i] + cell_center_offset;
		}
	} else {
		parent = cell_parent_oct(c);
		level = oct_level[parent];
		child = cell_child_number(c); 

		/* use delta to compute offset from parent's center */
		for ( i = 0; i < nDim; i++ ) {
			position[i] = (double)oct_pos[parent][i] + (double)cell_size[level] * (double)cell_delta[child][i];
		}
	}
}

int cell_find_position( double position[nDim] ) {
	int i;
	int coords[nDim];
	int root_index;
	int c;
	int child;
	float pos[nDim];

	for ( i = 0; i < nDim; i++ ) {
		coords[i] = (int)(position[i]);

		if ( coords[i]<0 || coords[i]>=num_grid ) {
			return -1;
		}
	}
	root_index = sfc_index( coords );

	c = root_cell_location(root_index);
	if ( c != NULL_OCT ) {	
		while ( cell_is_refined(c) ) {
			cell_position(c, pos);

			/* determine which child cell contains the point */
			child = 0;
			for ( i = 0; i < nDim; i++ ) {
				if ( position[i] >= pos[i] ) {
					child += (1<<i);
				}
			}

			c = cell_child( c, child );
		}
	}

	return c;
}

int cell_find_position_level( int level, double position[nDim] ) {
	int i;
	int curlevel;
        int coords[nDim];
        int root_index;
        int c;
        int child;
        float pos[nDim];

	for ( i = 0; i < nDim; i++ ) {
		coords[i] = (int)(position[i]);

		if ( coords[i]<0 || coords[i]>=num_grid ) {
			return -1;
		}
	}
	root_index = sfc_index( coords );
                                                                                                                                                            
	c = root_cell_location(root_index);
	curlevel = min_level;
	if ( c != NULL_OCT ) {
		while ( curlevel != level && cell_is_refined(c) ) {
			cell_position(c, pos);

			/* determine which child cell contains the point */
			child = 0;
			for ( i = 0; i < nDim; i++ ) {
				if ( position[i] >= pos[i] ) {
					child += (1<<i);
				}
			}

			c = cell_child( c, child );
			curlevel++;
		}

		if ( cell_level(c) != level ) {
			c = -1;
		}
	}

	return c;
}

int cell_find_position_above_level( int level, double position[nDim] ) {
	int i;
	int curlevel;
	int coords[nDim];
	int root_index;
	int c;
	int child;
	float pos[nDim];

	for ( i = 0; i < nDim; i++ ) {
		coords[i] = (int)(position[i]);

		if ( coords[i]<0 || coords[i]>=num_grid ) {
			return -1;
		}
	}
	
	root_index = sfc_index( coords );
	c = root_cell_location(root_index);
	curlevel = min_level;

	if ( c != NULL_OCT ) {
		while ( curlevel != level && cell_is_refined(c) ) {
			cell_position(c, pos);

			/* determine which child cell contains the point */
			child = 0;
			for ( i = 0; i < nDim; i++ ) {
				if ( position[i] >= pos[i] ) {
					child += (1<<i);
				}
			}

			c = cell_child( c, child );
			curlevel++;
		}
	}

	return c;
}

int cell_contains_position( int cell, double position[nDim] ) {
	int i;
	float pos[nDim];
	float size;

	size = 0.5*cell_size[cell_level(cell)]+1e-6;
	cell_position( cell, pos );

	for ( i = 0; i < nDim; i++ ) {
		if ( position[i] > pos[i] + size || position[i] < pos[i] - size ) {
			return 0;
		}
	}

	return 1;
}

double compute_distance_periodic( double pos1[nDim], double pos2[nDim] ) {
	int d;
	double dx, r = 0.0;
	
	for ( d = 0; d < nDim; d++ ) {
		dx = fabs( pos1[d] - pos2[d] );

		if ( dx > (double)(num_grid/2) ) {
			dx -= num_grid;
		}

		r += dx*dx;
	}

	return sqrt(r);
}

float compute_distance_periodic_float( float pos1[nDim], float pos2[nDim] ) {
	int d;
	float dx, r = 0.0;

	for ( d = 0; d < nDim; d++ ) {
		dx = fabsf( pos1[d] - pos2[d] );

		if ( dx > (float)(num_grid/2) ) {
			dx -= num_grid;
		}
		r += dx*dx;
	}

	return sqrtf(r);
}


/*******************************************************
 * cell_child
 *******************************************************/
int cell_child( int c, int j ) 
/* purpose: returns the jth child of the given cell or
 * 	-1 if the cell is unrefined (handled by 
 * 	oct_child).
 */
{
	cart_assert ( c >= 0 && c < num_cells );
	cart_assert ( j >= 0 && j < num_children );

	if ( cell_is_refined(c) ) {
		return oct_child( cell_child_oct[c], j );
	} else {
		return UNREFINED_CELL;
	}
}

/*******************************************************
 * oct_all_children
 *******************************************************/
void oct_all_children( int oct, int children[num_children] ) 
/* purpose: returns all children of the given oct.  Simply
 * 	wraps around oct_child.
 */
{
	int i;

	cart_assert( oct >= 0 && oct < num_octs );

	for ( i = 0; i < num_children; i++ ) {
		children[i] = oct_child( oct, i );
	}	
}

/*******************************************************
 * cell_all_children
 *******************************************************/
void cell_all_children( int c, int children[num_children] ) 
/* purpose: returns all children of the given cell
 */
{
	cart_assert ( c >= 0 && c <= num_cells );
	cart_assert( cell_is_refined(c) );

	oct_all_children( cell_child_oct[c], children );
}

/*******************************************************
 * root_cell_neighbor
 *******************************************************/
int root_cell_neighbor( int sfc, int direction ) 
/* purpose: computes the sfc index of the neighbor of root cell sfc
 * 	does not check whether that cell is local, buffered or remote
 */
{
	int i;
	int coords[nDim];

	cart_assert( direction >= 0 && direction < num_neighbors );
	cart_assert( sfc >= 0 && sfc < max_sfc_index );

	/* calculate the coordinates of the given root cell */
	sfc_coords( sfc, coords );

	/* shift the coordinates in the given direction */
	for ( i = 0; i < nDim; i++ ) {
		coords[i] += ishift[direction][i];

		/* boundary check to ensure periodicity */
		if ( coords[i] >= num_grid ) {
			coords[i] = 0;
		} else if ( coords[i] < 0 ) {
			coords[i] = num_grid - 1;
		}
	}

	/* return the sfc index of the revised coordinates */
	return sfc_index( coords );	
}

void root_cell_all_dependent_neighbors( int sfc, int neighbors[num_dependent_neighbors] ) {
	int d;

	for ( d = 0; d < num_neighbors; d++ ) {
		neighbors[d] = root_cell_neighbor( sfc, d );
	}

	for ( d = 0; d < num_secondary_neighbors; d++ ) {
		neighbors[num_neighbors+d] = root_cell_neighbor( neighbors[ secondary_neighbors[d][0] ],
				secondary_neighbors[d][1] );
	}

	for ( d = 0; d < num_tertiary_neighbors; d++ ) {
		neighbors[num_neighbors+num_secondary_neighbors+d] = 
			root_cell_neighbor( neighbors[ num_neighbors + tertiary_neighbors[d][0] ],
				tertiary_neighbors[d][1] );		
	}
}

/*******************************************************
 * cell_neighbor
 *******************************************************/
int cell_neighbor( int c, int direction ) 
/* purpose: returns the neighbor of the given cell
 * 	in the given direction (either at the same
 * 	level or cell_level(c) - 1
 */
{
	int cellnum, parent, neighbor, neighbor_cell;

	cart_assert ( c >= 0 && c <= num_cells );
	cart_assert ( direction >= 0 && direction < num_neighbors );

	if ( cell_is_root_cell(c) ) {
		return root_cell_location( root_cell_neighbor( root_cell_sfc_index(c), direction) );
	} else {
		cellnum = cell_child_number(c);
		neighbor = local[cellnum][direction];

		/* if the cell is in the same oct, simply shift cell pointer,
		 * otherwise use oct's neighbor array */
		if ( in_local_oct[cellnum][direction] ) {
			return (c - cellnum) + neighbor;
		} else {
			parent = cell_parent_oct(c);
			cart_assert( parent >= 0 && parent < num_octs );
			
			neighbor_cell = oct_neighbors[parent][direction];
			cart_assert( neighbor_cell == NULL_OCT || ( neighbor_cell >= 0 && neighbor_cell < num_cells ) );
	
			if ( neighbor_cell != NULL_OCT && cell_is_refined(neighbor_cell) ) {
					return cell_child( neighbor_cell, neighbor );
			} else {
				return neighbor_cell;
			}
		}
	}
}

/*******************************************************
 * cell_all_neighbors
 *******************************************************/
void cell_all_neighbors( int c, int neighbors[num_neighbors] ) 
/* purpose: calculates all neighbors of the given cell c
 */
{
	int i;
	int cellnum, parent, neighbor, neighbor_cell;
	int direction, sfc;
	int coords[nDim], coords_neighbor[nDim];

	cart_assert ( c >= 0 && c < num_cells );

	if ( cell_is_root_cell(c) ) {
		sfc = root_cell_sfc_index(c);
		sfc_coords( sfc, coords );

		for ( direction = 0; direction < num_neighbors; direction++ ) {
			/* shift the coordinates in the given direction */
			for ( i = 0; i < nDim; i++ ) {
				coords_neighbor[i] = coords[i] + ishift[direction][i];

				/* boundary check to ensure periodicity */
				if ( coords_neighbor[i] >= num_grid ) {
					coords_neighbor[i] = 0;
				} else if ( coords_neighbor[i] < 0 ) {
					coords_neighbor[i] = num_grid - 1;
				}
			}
	
			neighbors[direction] = root_cell_location( sfc_index( coords_neighbor ) );
		}	
	} else {
		cellnum = cell_child_number(c);
		parent = cell_parent_oct(c);

		for ( direction = 0; direction < num_neighbors; direction++ ) {
			neighbor = local[cellnum][direction];
	
			if ( in_local_oct[cellnum][direction] ) {
				neighbors[direction] = (c - cellnum) + neighbor;
			} else {
				neighbor_cell = oct_neighbors[parent][direction];

				if ( neighbor_cell != NULL_OCT && cell_is_refined(neighbor_cell) ) {
					neighbors[direction] = cell_child( neighbor_cell, neighbor );
				} else {
					neighbors[direction] = neighbor_cell;
				}
			}
		}
	}
}


#ifdef HYDRO
float cell_gas_kinetic_energy(int icell)
{
  int j;
  double ke = 0.0;

  if(cell_gas_density(icell) > 0.0)
    {
      for(j=0; j<nDim; j++) ke += cell_momentum(icell,j)*cell_momentum(icell,j);
      return (float)(0.5*ke/cell_gas_density(icell));
    }
  else return 0.0;
}
#endif

