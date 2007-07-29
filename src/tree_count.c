#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defs.h"
#include "tree.h"
#include "iterators.h"
#include "auxiliary.h"

/*******************************************************
 * cell_count_leaves
 ******************************************************/
int cell_count_leaves( int c, int level ) 
/* purpose: workfunction for use in counting all
 * 	leaves in a tree
 */
{
	cart_assert( c >= 0 && c < num_cells );
	
	if ( cell_is_leaf(c) ) {
		return 1;
	} else {
		return 0;
	}
}

/*******************************************************
 * cell_num_child_leaves
 ******************************************************/
int cell_num_child_leaves( int c, int level ) 
/* purpose: uses tree_traversal iterator to count the
 * 	number of leaves in a root cell tree
 */
{
	cart_assert( cell_is_root_cell(c) );
	return tree_traversal( c, cell_count_leaves );
}

/*******************************************************
 * cell_count_cells
 ******************************************************/
int cell_count_cells( int c, int level ) 
/* purpose: workfunction for use in counting all cells
 * 	in a tree.
 */
{
	cart_assert( c >= 0 && c < num_cells );
	
	return 1;
}

/*******************************************************
 * cell_num_child_cells
 ******************************************************/
int cell_num_child_cells( int c, int level ) 
/* purpose: uses tree_traversal iterator and 
 * 	cell_count_cells to count number of cells in
 * 	a root cell tree.
 */
{
	cart_assert( cell_is_root_cell(c) );
	
	return tree_traversal( c, cell_count_cells );
}

/*******************************************************
 * cell_count_octs
 ******************************************************/
int cell_count_octs( int c, int level ) 
/* purpose: workfunction for use in counting all octs
 * 	in a tree (essentially num_child_cells -
 * 	num_leaves).
 */
{
	cart_assert( c >= 0 && c < num_cells );
	
	if ( cell_is_refined(c) ) {
		return 1;
	} else {
		return 0;
	}
}

/*******************************************************
 * tree_num_octs
 ******************************************************/
int tree_num_octs( int c, int level ) 
/* purpose: counts the number of octs in a root cell
 * 	tree
 */
{
	cart_assert( c >= 0 && c < num_root_cells );

	return tree_traversal( c, cell_count_octs );
}

/*******************************************************
 * tree_num_cells
 ******************************************************/
int tree_num_cells( int c, int level )
/* purpose: counts the number of cells in a root cell
 *      tree
 */
{
        cart_assert( c >= 0 && c < num_root_cells );
        return tree_traversal( c, cell_count_cells );
}
