#ifndef __SFC_H__
#define __SFC_H__

#define nDim		3

#define	SLAB		0
#define MORTON		1
#define HILBERT		2

#define SLAB_DIM	0

#define num_refinement_levels         20 //needs to be higher than maxlevel-minlevel+1   

extern int num_grid;
extern int sfc_order;
extern int nBitsPerDim;
extern int nBits;
extern int max_sfc_index;

#define rollLeft(x,y,mask) ((x<<y) | (x>>(nDim-y))) & mask
#define rollRight(x,y,mask) ((x>>y) | (x<<(nDim-y))) & mask

void init_sfc();
int sfc_index( int coords[nDim] );
void sfc_coords( int index, int coords[nDim] );

#endif
