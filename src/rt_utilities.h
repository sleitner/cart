#ifndef __RT_UTILITIES_H__
#define __RT_UTILITIES_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>

/*
//  Useful macros
*/
#define MESH_RUN_DECLARE(level,cell) \
int cell, level, _Index, _MaxLevel = max_level_local(); \
int _Num_level_cells, *_Level_cells = 0


#define MESH_RUN_OVER_LEVELS_BEGIN(level,level1,level2) \
for(level=level1; level<=level2; level++) \
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


#define rtuStencilSize     (2*nDim*nDim)
extern double rtuStencilDist2[rtuStencilSize];
extern double rtuStencilDelPos[rtuStencilSize][nDim];
extern double rtuStencilTensor[rtuStencilSize][nDim*(nDim+1)/2];


void rtuInitRun();

/*
// Helper functions
*/
void rtuGlobalAverage(int n, double *buffer);
void rtuGetStencil(int level, int cell, int nb[]);
void rtuGetLinearArrayMaxMin(int n, float *arr, float *max, float *min);

void rtuCopyArraysInt(int *dest, int *src, int size);
void rtuCopyArraysFloat(float *dest, float *src, int size);

#ifdef DEBUG
void rtuCheckGlobalValue(int val, char *name, MPI_Comm local_comm);
#endif


#endif /* __RT_UTILITIES_H__ */
