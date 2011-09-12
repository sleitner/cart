#ifndef __OUTPUTSLICE_H__
#define __OUTPUTSLICE_H__

#define OUT_CELL_DENSITY 1
#define OUT_CELL_INTERNAL_ENERGY (OUT_CELL_DENSITY+1)
#define OUT_CELL_MOMENTUM (OUT_CELL_INTERNAL_ENERGY+1)
#define OUT_CELL_SOUNDSPEED (OUT_CELL_MOMENTUM+3)
#define OUT_CELL_MACH (OUT_CELL_SOUNDSPEED+1)
#define OUT_CELL_LEVEL (OUT_CELL_MACH+1)

extern const int axis_direction[nDim][nDim-1];

void dump_plane(int iflag, int out_level, int slice_axis_z, double pos_zero[nDim], double slice_region_size, FILE *output);

#endif /* __OUTPUTSLICE_H__ */
