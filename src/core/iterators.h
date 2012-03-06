#ifndef __ITERATORS_H__
#define __ITERATORS_H__


void select_level( int level, int cell_types, int *num_cells_selected, int **selection );
void select_level_octs( int level, int oct_types, int *num_octs_selected, int **selection );
void select_local_buffered_cells( int level, int proc, int *num_cells_selected, int **selection );

int tree_cell_count( int cell );
int tree_traversal( int cell, int workfunction( int, int ) );
int tree_preorder_traversal( int cell, int workfunction(int, int) );
int tree_level_traversal( int cell, int level, int workfunction( int, int ) );

int tree_oct_count( int icell );
void root_tree_level_oct_count( int icell, int level_oct_count[max_level-min_level+1] );
void tree_level_oct_count( int icell, int level_oct_count[max_level-min_level+1] );
int tree_max_level( int icell );

/*
//  Useful macros
*/
#define MESH_RUN_DECLARE(level,cell) \
int cell, level, _Index, _MaxLevel = max_level_local(); \
int _Num_level_cells, *_Level_cells = 0


#define MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2) \
for(level=level1; (level1>level2)?level>=level2:level<=level2; (level1>level2)?level--:level++) \
{ \
  select_level(level,CELL_TYPE_LOCAL,&_Num_level_cells,&_Level_cells)

#define MESH_RUN_OVER_ALL_LEVELS_BEGIN(level) \
MESH_RUN_OVER_LEVELS_BEGIN(level,min_level,_MaxLevel) 

#define MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell) \
for(_Index=0; _Index<_Num_level_cells; _Index++) \
{ \
  cell = _Level_cells[_Index]

#define MESH_RUN_OVER_CELLS_OF_LEVEL_END \
}

#define MESH_RUN_OVER_LEVELS_END \
 cart_free(_Level_cells); \
 _Level_cells = 0; \
}

#endif
