#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

typedef struct HALO {
	int id;
	int proc;
	int sfc_index;
	float pos[nDim];
	float vel[nDim];
	float rhalo;
	float rvir;
	float mvir;
	float vmax;
	int np;
	int *particles;
	int *binding_order;
	int section;
} halo;

typedef struct HALO_LIST {
	int num_halos;
	halo *list;
} halo_list;

void load_halo_finder_epochs( char *filename, int *num_epochs, float **epoch );
halo_list *load_halo_finder_catalog( char *filename );
void load_halo_particle_mapping( char *filename, halo_list *halos );
void destroy_halo_list( halo_list *halos );
void compute_halo_properties( char *analysis_directory, int halo_section, halo_list *halos );
void crude_stellar_mass_fractions( halo_list *halos );

#endif
