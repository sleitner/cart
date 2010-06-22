#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "defs.h"
#include "auxiliary.h"
#include "halos.h"

#define subhalo_exclude_radius	(1.5)
#define EXCLUDE_SUBHALOS

typedef struct {
} profile_struct;

void compute_halo_properties( char *output_directory, char *analysis_directory, halo_list *halos, halo_list *subhalos );

#endif
