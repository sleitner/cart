#ifndef __HALO_FINDER_H__
#define __HALO_FINDER_H__

#ifdef PARTICLES

#define MAX_HALO_BINS		150

typedef struct HALO {
	int hid;
    double x[nDim];
	double v[nDim];
    double rdelta;
    double Mdelta;
	double rmax;
    double vmax;
	int np;
	struct HALO *next;
} halo;

extern double min_halo_mass;
extern int halo_finder_frequency;
char halo_finder_output_directory[];

void config_init_halo_finder();
void config_verify_halo_finder();
void compute_halo_mass( halo *h );
void halo_recenter( halo *h );
halo *find_halos();
void write_halo_list( halo *list );
void write_halo_particle_list( halo *list );
void destroy_halo_list( halo *list );

#endif /* PARTICLES */

#endif /* __HALO_FINDER_H__ */

