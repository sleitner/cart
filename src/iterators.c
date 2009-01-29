#include <stdlib.h>
#include <stdio.h>

#include "defs.h"
#include "tree.h"
#include "iterators.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "tree_linkedlist.h"
#include "auxiliary.h"
#include "timing.h"

void select_level( int level, int cell_types, int *num_cells_selected, int **selection ) {
	int i;
	int num_selected;
	int *level_cells;
	int ioct, ichild;

	start_time( SELECT_LEVEL_TIMER );

	cart_assert( level >= min_level && level <= max_level );
	cart_assert( cell_types == CELL_TYPE_LOCAL || cell_types == CELL_TYPE_BUFFER || cell_types == CELL_TYPE_ANY );

	switch ( cell_types ) {
		case CELL_TYPE_LOCAL:
			*num_cells_selected = num_cells_per_level[level];
			break;
		case CELL_TYPE_BUFFER:
			*num_cells_selected = num_buffer_cells[level];
			break;
		case CELL_TYPE_ANY:
			*num_cells_selected = num_cells_per_level[level] + num_buffer_cells[level];
			break;
	}

	level_cells = cart_alloc( *num_cells_selected * sizeof(int) );
	num_selected = 0;

	if ( cell_types == CELL_TYPE_LOCAL || cell_types == CELL_TYPE_ANY ) {
		if ( level == min_level ) {
			cart_assert( num_cells_per_level[min_level] <= *num_cells_selected );
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

	if ( cell_types == CELL_TYPE_BUFFER || cell_types == CELL_TYPE_ANY ) {
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

	cart_assert( num_selected == *num_cells_selected );

	*selection = level_cells;

	end_time( SELECT_LEVEL_TIMER );
}

void select_level_with_condition( int select_leaves, int level, int *num_cells_selected, int **selection )
{
  int i;
  int num_selected, size;
  int *level_cells;
  int ioct, ichild;

  start_time( SELECT_LEVEL_TIMER );

  cart_assert( level >= min_level && level <= max_level );
  cart_assert( select_leaves == 0 || select_leaves == 1 );

  /*
  //  Measure array size
  */
  size = 0;
  if ( level == min_level )
    {
      for ( i = 0; i < num_cells_per_level[min_level]; i++ ) if(cell_is_leaf(i) == select_leaves) 
	{
	  size++;
	}
    }
  else
    {
      ioct = local_oct_list[level];
      while ( ioct != NULL_OCT )
	{
	  for ( ichild = 0; ichild < num_children; ichild++ )
	    {
	      i = oct_child(ioct,ichild);
	      if(cell_is_leaf(i) == select_leaves) 
		{
		  size++;
		}
	    }
	  ioct = oct_next[ioct];
	}
    }

  if(size == 0)
    {
      *selection = NULL;
      *num_cells_selected = 0;
    }
  else
    {
      level_cells = cart_alloc( size * sizeof(int) );
      num_selected = 0;

      if ( level == min_level )
	{
	  for ( i = 0; i < num_cells_per_level[min_level]; i++ ) if(cell_is_leaf(i) == select_leaves) 
	    {
	      level_cells[num_selected++] = i;
	    }
	}
      else
	{
	  ioct = local_oct_list[level];
	  while ( ioct != NULL_OCT )
	    {
	      for ( ichild = 0; ichild < num_children; ichild++ )
		{
		  i = oct_child(ioct,ichild);
		  if(cell_is_leaf(i) == select_leaves) 
		    {
		      level_cells[num_selected++] = i;
		    }
		}
	      ioct = oct_next[ioct];
	    }
	}

      *selection = level_cells;
      *num_cells_selected = num_selected;

    }

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
