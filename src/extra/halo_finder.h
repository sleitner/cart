#ifndef __HALO_FINDER_H__
#define __HALO_FINDER_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>


#define min_halo_particles      (1e3)
#define min_xray_removal_mass   (1e12)
#define xray_removal_radius     (1.0)

typedef struct HALO {
	int id;
	int proc;
	float vel[nDim];
	float rhalo;
	float rvir;
	float mvir;
	float vmax;
	int np;
	int *particles;
	int *binding_order;
	int section;
	double pos[nDim];
} halo;

typedef struct HALO_LIST {
	int num_halos;
	halo *list;
	int *map;  
} halo_list;

void load_halo_finder_epochs( char *filename, int *num_epochs, float **epoch );
halo_list *load_halo_finder_catalog( const char *filename, int nmem_min, float mvir_min, float vmax_min, float rvir_min);
void load_halo_particle_mapping( char *filename, halo_list *halos );
void destroy_halo_list( halo_list *halos );
void compute_halo_properties( char *analysis_directory, int halo_section, halo_list *halos );
void crude_stellar_mass_fractions( halo_list *halos );

int halo_level(const halo *h, MPI_Comm local_comm);
halo* find_halo_by_id(halo_list *halos, int id);

void dump_region_around_halo(const char *filename, const halo *h, float size);

/*
//  Set halos->map with the halo id for each halo, or 0 if belongs to 
//  none; a cell belongs to a halo if it is inside its size_factor*Rtrunc, 
//  and satellites are accounted for properly.
*/
void map_halos(int resolution_level, halo_list *halos, float size_factor);

#endif
