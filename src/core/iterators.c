#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "parallel.h"
#include "timing.h"
#include "tree.h"
#include "tree_linkedlist.h"


void select_level( int level, int cell_types, int *num_cells_selected, int **selection ) {
	int i;
	int num_selected;
	int *level_cells;
	int ioct, ichild;
	int refined_flag;

	start_time( SELECT_LEVEL_TIMER );

	cart_assert( level >= min_level && level <= max_level );

	/* hack to allow backwards compatibility for cell_types */
	if ( !((cell_types & CELL_TYPE_LEAF) || (cell_types & CELL_TYPE_REFINED)) ) {
		cell_types |= CELL_TYPE_LEAF | CELL_TYPE_REFINED;
	}
	
	num_selected = 0;

	if ( cell_types & CELL_TYPE_LOCAL ) {
		if ( cell_types & CELL_TYPE_LEAF ) {
			if ( level < max_level ) {
				num_selected += num_cells_per_level[level] - num_cells_per_level[level+1]/num_children;
			} else {
				num_selected += num_cells_per_level[level];
			}
		} 
		if ( (cell_types & CELL_TYPE_REFINED) && (level < max_level) ) {
			num_selected += num_cells_per_level[level+1]/num_children;
		}
	}

	if ( cell_types & CELL_TYPE_BUFFER ) {
		if ( cell_types & CELL_TYPE_LEAF ) {
			if ( level < max_level ) {
				num_selected += num_buffer_cells[level] - num_buffer_cells[level+1]/num_children;
			} else {
				num_selected += num_buffer_cells[level];
			}
		} 
        if ( (cell_types & CELL_TYPE_REFINED) && (level < max_level) ) {
            num_selected += num_buffer_cells[level+1]/num_children;
        }
    }

	level_cells = cart_alloc(int, num_selected );
	*num_cells_selected = num_selected;
	num_selected = 0;

	if ( (cell_types & CELL_TYPE_LEAF) && (cell_types & CELL_TYPE_REFINED) ) {
		if ( cell_types & CELL_TYPE_LOCAL ) {
			if ( level == min_level ) {
				for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
					level_cells[i] = i;
				}
				num_selected = num_cells_per_level[min_level];
			} else {
				ioct = local_oct_list[level];
				while ( ioct != NULL_OCT ) {
					for ( ichild = 0; ichild < num_children; ichild++ ) {
						level_cells[num_selected++] = oct_child(ioct,ichild);
					}
					ioct = oct_next[ioct];
				}
			}
		}

		if ( cell_types & CELL_TYPE_BUFFER ) {
			if ( level == min_level ) { 
				for ( i = num_cells_per_level[min_level]; 
						i < num_cells_per_level[min_level]+num_buffer_cells[min_level]; i++ ) {
					level_cells[num_selected++] = i;
				}
			} else {
				ioct = buffer_oct_list[level];
				while ( ioct != NULL_OCT ) {
					for ( ichild = 0; ichild < num_children; ichild++ ) {
						level_cells[num_selected++] = oct_child(ioct,ichild);
					}
					ioct = oct_next[ioct];
				}
			}
		}
	} else {
		refined_flag = ( cell_types & CELL_TYPE_REFINED ) ? 1 : 0;

		if ( cell_types & CELL_TYPE_LOCAL ) {
			if ( level == min_level ) {
				for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
					if ( cell_is_refined(i) == refined_flag ) {
						level_cells[num_selected++] = i;	
					}
				}
			} else {
				ioct = local_oct_list[level];
				while ( ioct != NULL_OCT ) {
					for ( ichild = 0; ichild < num_children; ichild++ ) {
						if ( cell_is_refined( oct_child(ioct,ichild) ) == refined_flag ) {
							level_cells[num_selected++] = oct_child(ioct,ichild);
						}
					}
					ioct = oct_next[ioct];
				}
			}
		}

		if ( cell_types & CELL_TYPE_BUFFER ) {
			if ( level == min_level ) { 
				for ( i = num_cells_per_level[min_level]; 
						i < num_cells_per_level[min_level]+num_buffer_cells[min_level]; i++ ) {
					if ( cell_is_refined(i) == refined_flag ) {
						level_cells[num_selected++] = i;
					}
				}
			} else {
				ioct = buffer_oct_list[level];
				while ( ioct != NULL_OCT ) {
					for ( ichild = 0; ichild < num_children; ichild++ ) {
						if ( cell_is_refined( oct_child(ioct,ichild) ) == refined_flag ) {
							level_cells[num_selected++] = oct_child(ioct,ichild);
						}
					}
					ioct = oct_next[ioct];
				}
			}
		}
	}

	cart_assert( num_selected == *num_cells_selected );
	*selection = level_cells;

	end_time( SELECT_LEVEL_TIMER );
}

void select_level_octs( int level, int oct_types, int *num_octs_selected, int **selection ) {
	int num_selected;
	int *level_octs;
	int ioct;

	start_time( SELECT_LEVEL_TIMER );

	cart_assert( level > min_level && level <= max_level );
	cart_assert( !((oct_types & CELL_TYPE_LEAF) || (oct_types & CELL_TYPE_REFINED) ) );

	num_selected = 0;
	if ( oct_types & CELL_TYPE_LOCAL ) {
		num_selected += num_cells_per_level[level]/num_children;
	}
	if ( oct_types & CELL_TYPE_BUFFER ) {
		num_selected += num_buffer_cells[level]/num_children;
	}

	level_octs = cart_alloc(int, num_selected );
	*num_octs_selected = num_selected;
	num_selected = 0;

	if ( oct_types & CELL_TYPE_LOCAL ) {
		ioct = local_oct_list[level];
		while ( ioct != NULL_OCT ) {
			level_octs[num_selected++] = ioct;
			ioct = oct_next[ioct];
		}
	}

	if ( oct_types & CELL_TYPE_BUFFER ) {
		ioct = buffer_oct_list[level];
		while ( ioct != NULL_OCT ) {
			level_octs[num_selected++] = ioct;
			ioct = oct_next[ioct];
		}
	}

	cart_assert( num_selected == *num_octs_selected );
	*selection = level_octs;

	end_time( SELECT_LEVEL_TIMER );
}

void select_local_buffered_cells( int level, int proc, int *num_cells_selected, int **selection ) {
	int i, j;
	int num_selected;
	int *level_cells;

	cart_assert( level >= min_level && level <= max_level );
	cart_assert( proc >= 0 && proc < num_procs );

	start_time( SELECT_LEVEL_TIMER );

	if ( level == min_level ) {
		*num_cells_selected = num_remote_buffers[min_level][proc];
	} else {
		*num_cells_selected = num_children*num_remote_buffers[level][proc];
	}

	level_cells = cart_alloc(int, *num_cells_selected);

	if ( level == min_level ) {
		for ( i = 0; i < num_remote_buffers[min_level][proc]; i++ ) {
			level_cells[i] = root_cell_location( remote_buffers[min_level][proc][i] );
		}
		num_selected = num_remote_buffers[min_level][proc];
	} else {
		num_selected = 0;
		for ( i = 0; i < num_remote_buffers[level][proc]; i++ ) {
			for ( j = 0; j < num_children; j++ ) {
				level_cells[num_selected++] = oct_child( remote_buffers[level][proc][i], j );
			}
		}
	}

	cart_assert( num_selected == *num_cells_selected );
	*selection = level_cells;

	end_time( SELECT_LEVEL_TIMER );
}

/*******************************************************
 * tree_traversal
 *******************************************************/
int tree_traversal( int cell, int workfunction( int, int ) ) 
/* purpose: recursively travels the tree, applying workfunction
 * 	to every cell.
 *
 * returns: the sum of all returns from the workfunction 
 * 	in this subtree.
 */
{
	int i;
	int count = 0;

	cart_assert( cell >= 0 && cell < num_cells );

	/* if this isn't a leaf cell, recursively traverse tree */	
	if ( cell_is_refined( cell ) ) {
		for ( i = 0; i < num_children; i++ ) {
			count += tree_traversal( cell_child( cell, i ), workfunction );
		}
	}

	count += workfunction( cell, cell_level(cell) );

	return count;
}

int tree_preorder_traversal( int cell, int workfunction(int, int) ) {
	int i;
	int count = 0;

	cart_assert( cell >= 0 && cell < num_cells );

	count = workfunction( cell, cell_level(cell) );

	if ( count > 0 ) {
		/* if this isn't a leaf cell, recursively traverse tree */      
		if ( cell_is_refined( cell ) ) {
			for ( i = 0; i < num_children; i++ ) {
				count += tree_traversal( cell_child( cell, i ), workfunction );
			}
		}
	}

	return count;
}

int tree_cell_count( int cell ) {
	int child;
	int count;

	cart_assert( cell >= 0 && cell < num_cells );

	count = 1;

	if ( cell_is_refined(cell) ) {
		for ( child = 0; child < num_children; child++ ) {
			count += tree_cell_count( cell_child( cell, child ) );
		}
	}

	return count;
}

int tree_oct_count( int icell ) {
	int i;
	int count = 0;
	if ( cell_is_refined(icell) ) {
		for ( i = 0; i < num_children; i++ ) {
			count += tree_oct_count( cell_child(icell,i));
		}
		count++;
	}

	return count;
}

void root_tree_level_oct_count( int icell, int level_oct_count[max_level-min_level] ) {
	int i;

	cart_assert( cell_level(icell) == min_level );

	for ( i = min_level; i < max_level; i++ ) {
		level_oct_count[i] = 0;
	}

	tree_level_oct_count( icell, level_oct_count );
}
    

void tree_level_oct_count( int icell, int level_oct_count[max_level-min_level] ) {
	int i;

	if ( cell_is_refined(icell) ) {
		level_oct_count[ cell_level(icell) ]++;

		for ( i = 0; i < num_children; i++ ) {
			tree_level_oct_count( cell_child( icell, i ), level_oct_count );
		}
	}
}

int tree_max_level( int icell ) {
	int i, level, max_tree_level;

	if ( cell_is_refined(icell) ) {
		i = 0;
		max_tree_level = tree_max_level( cell_child( icell, i ) );
		for ( ; i < num_children; i++ ) {
			level = tree_max_level( cell_child( icell, i ) );

			if ( level > max_tree_level ) {
				max_tree_level = level;
			}
		}
	} else {
		max_tree_level = cell_level(icell);
	}

	return max_tree_level;
}

/*******************************************************
 * tree_level_traversal
 *******************************************************/
int tree_level_traversal( int cell, int level, int workfunction( int, int ) ) 
/* purpose: recursively travels the tree, applying workfunction
 * 	to all cells on the given level.
 *
 * returns: the sum of all returns from the workfunction.
 */
{
	int i;
	int count = 0;

	cart_assert( level >= min_level && level <= max_level );
	cart_assert( cell_level(cell) <= level );

	if ( cell_level(cell) == level ) {
		count += workfunction(cell, level);
	} else {
		for ( i = 0; i < num_children; i++ ) {
                        count += tree_level_traversal( cell_child( cell, i ), level, workfunction );
                }
	}

	return count;
}
