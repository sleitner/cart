#ifndef __HALOS_H__
#define __HALOS_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef COSMOLOGY

#include <mpi.h>

typedef struct HALO {
	int id;
	double pos[nDim];
	double vel[nDim];
	double rhalo;
	double rvir;
	double mvir;
	double vmax;
	double rmax;
	int np;
	int *particles;
	int *binding_order;
	int flag;
} halo;

typedef struct HALO_LIST {
	int num_halos;
	int size;
	halo *list;
	int *map;  
} halo_list;

void load_halo_finder_epochs( char *filename, int *num_epochs, float **epoch );
halo_list *load_halo_finder_catalog( const char *filename, int nmem_min, float mvir_min, float vmax_min, float rvir_min, int max_num_halos);
void load_halo_particle_mapping( char *filename, halo_list *halos );
void destroy_halo_list( halo_list *halos );
halo *find_halo_by_id(halo_list *halos, int id);
halo_list *halo_list_alloc( int size );
halo *halo_list_add_halo( halo_list *halos );

#endif /* COSMOLOGY */
#endif
