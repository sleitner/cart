#ifndef __VIEWDUMP_H__
#define __VIEWDUMP_H__

#define DUMP_HVARS		0
#define DUMP_GRAV_VARS		1
#define DUMP_REFINEMENT_VARS	2
#define DUMP_DIFFUSION_VARS	3
#define DUMP_LOCATION_VARS	4

void viewdump( const char *filename, int max_level_dumped, float slice, int d, int vars, int cell_type );
void viewdump_serial( const char *filename, int max_level_dumped, float slice, int d, int vars, int cell_type );

extern float drho[nDim];
extern int neighbors[num_neighbors];

#endif
