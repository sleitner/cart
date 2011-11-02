#ifndef __ITERATORS_H__
#define __ITERATORS_H__


void select_level( int level, int cell_types, int *num_cells_selected, int **selection );
void select_level_octs( int level, int oct_types, int *num_octs_selected, int **selection );
void select_local_buffered_cells( int level, int proc, int *num_cells_selected, int **selection );

int tree_cell_count( int cell );
int tree_traversal( int cell, int workfunction( int, int ) );
int tree_preorder_traversal( int cell, int workfunction(int, int) );
int tree_level_traversal( int cell, int level, int workfunction( int, int ) );


#endif
