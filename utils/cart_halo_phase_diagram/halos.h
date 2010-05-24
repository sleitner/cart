#ifndef __HALOS_H__
#define __HALOS_H__

#include "defs.h"
#include "auxiliary.h"

typedef struct {
	int id;
	double pos[nDim];
	double vel[nDim];
	double rvir;
	double mvir;
	double rhalo;
	double analysis_radius;
	void *data;
} halo_struct;

typedef struct {
	int num_halos;
	halo_struct *list;
} halo_list;

halo_list *load_halo_sample_file( char *filename );
void load_halo_finder_catalog( char *filename, halo_list * );
void load_baryon_catalog( char *filename, halo_list * );
void destroy_halo_list( halo_list *halos );

#endif
