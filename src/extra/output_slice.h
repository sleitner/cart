#define OUT_CELL_DENSITY 1
#define OUT_CELL_INTERNAL_ENERGY (OUT_CELL_DENSITY+1)
#define OUT_CELL_MOMENTUM (OUT_CELL_INTERNAL_ENERGY+1)
#define OUT_CELL_PRESSURE (OUT_CELL_MOMENTUM+3)
#define OUT_CELL_LEVEL (OUT_CELL_PRESSURE+1)

void dump_plane(int iflag, int out_level, int slice_axis_z, double pos_zero[nDim], double slice_region_size, FILE *output);
    
