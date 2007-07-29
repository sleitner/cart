#ifndef __ITERATORS_H__
#define __ITERATORS_H__

#include "tree.h"
#include "cell_buffer.h"

void select_level( int level, int cell_types, int *num_cells_selected, int **selection );

int tree_cell_count( int cell );
int tree_traversal( int cell, int workfunction( int, int ) );
int tree_preorder_traversal( int cell, int workfunction(int, int) );
int tree_level_traversal( int cell, int level, int workfunction( int, int ) );

#endif
