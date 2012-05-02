#ifndef __HALO_FINDER_H__
#define __HALO_FINDER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#include "halos.h"
#include <mpi.h>

int halo_level( const halo *h, MPI_Comm local_comm );
void dump_region_around_halo(const char *filename, const halo *h, float size);

/*
//  Set halos->map with the halo id for each halo, or 0 if belongs to 
//  none; a cell belongs to a halo if it is inside its size_factor*Rtrunc, 
//  and satellites are accounted for properly.
*/
void map_halos(int resolution_level, halo_list *halos, float size_factor);

#endif
