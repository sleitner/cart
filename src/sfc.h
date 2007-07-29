#ifndef __SFC_H__
#define __SFC_H__

#include "defs.h"
#include "tree.h"

#define	SLAB		0
#define MORTON		1
#define HILBERT		2

#ifndef SFC
#define SFC		HILBERT
#else
#if SFC == SLAB
#ifndef SLAB_DIM
#define SLAB_DIM	1
#endif
#endif
#endif


#define nBitsPerDim	num_root_grid_refinements
#define nBits		(nDim * nBitsPerDim)
#define max_sfc_index	(1<<nBits)

#define rollLeft(x,y,mask) ((x<<y) | (x>>(nDim-y))) & mask
#define rollRight(x,y,mask) ((x>>y) | (x<<(nDim-y))) & mask

int sfc_index( int coords[nDim] );
void sfc_coords( int index, int coords[nDim] );

#endif
