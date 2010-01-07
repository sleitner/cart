#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>


//#define rbinmax         (5e3/hubble)
//#define rbinmin         (10.0/hubble)
//#define max_bins        100

//#define rbinvirmin      0.05
//#define rbinvirmax      4.0
//#define num_vir_bins    50
//#define virial_radius_index     2       /* r500, 0 = rvir */

//#define points_per_cell 1

//#define Tcold           (1e5)
//#define new_star_age    (0.1)

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
} halo_list;

void load_halo_finder_epochs( char *filename, int *num_epochs, float **epoch );
halo_list *load_halo_finder_catalog( const char *filename, int nmem_min, float mvir_min, float vmax_min, float rvir_min);
void load_halo_particle_mapping( char *filename, halo_list *halos );
void destroy_halo_list( halo_list *halos );
void compute_halo_properties( char *analysis_directory, int halo_section, halo_list *halos );
void crude_stellar_mass_fractions( halo_list *halos );
int halo_level( const halo *h, MPI_Comm local_comm );

#endif
