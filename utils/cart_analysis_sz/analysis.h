#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "defs.h"
#include "auxiliary.h"
#include "skiplist.h"

#define PIXEL_NSIDE	5
#define PIXEL_NPIX	(6*PIXEL_NSIDE*PIXEL_NSIDE) /* really 12, but half are duplicates */
#define num_points_per_cell 100

//#define SCALED_PROJECTION_EXTENT

#define PROJECTED_PROFILES

/* in kpc/h comoving */
#define rbinmax         (4e3)
#define rbinmin         (5.0)
#define max_bins        50

extern double projection_extent;

/*
#define projection_extent (30./r0/aexpn)
*/

typedef struct {
	int id;
	double pos[nDim];
	double vel[nDim];
	double rvir;
	double mvir;

	double extent;

	double projected_sz_flux[PIXEL_NPIX];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double projected_electron_sz_flux[PIXEL_NPIX];
#endif

#ifdef PROJECTED_PROFILES
    double projected_sz_binned[PIXEL_NPIX][max_bins];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
    double projected_electron_sz_binned[PIXEL_NPIX][max_bins];
#endif
#endif

	double projection_volume[PIXEL_NPIX];
	

} halo_struct;

typedef struct {
	int num_halos;
	halo_struct *list;
} halo_list;

halo_list *load_halo_sample_file( char *filename );
void load_halo_finder_catalog( char *filename, halo_list * );
void load_baryon_catalog( char *filename, halo_list * );
void destroy_halo_list( halo_list *halos );

void compute_halo_properties( char *output_directory, char *analysis_directory, char *radius, halo_list *halos );

#endif
