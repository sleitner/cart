#ifndef __IO_H__
#define __IO_H__

#include "defs.h"
#include "analysis.h"
#include "halos.h"

#define MAX_OUTPUTS	256
#define MAX_PARTICLE_SPECIES	10

extern char output_directory[256];
extern char logfile_directory[256];
extern char jobname[256];
extern int num_output_files;

void reorder( char *buffer, int size );

#ifdef HYDRO
typedef struct {
	int level;
	double pos[nDim];
	float gas_density;
	float gas_energy;
	float gas_pressure;
	float gas_internal_energy;
	float momentum[nDim];
#ifdef ENRICH
	float metallicity_II;
#ifdef ENRICH_SNIa
	float metallicity_Ia;
#endif /* ENRICH_Ia */
#endif /* ENRICH */
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	float electron_internal_energy;
#endif
	int subhalo_flag;
} cell_struct;

void read_indexed_grid( char *, halo_list *, halo_list *, void callback( halo_struct *, cell_struct * ) );
#endif /* HYDRO */

#ifdef PARTICLES
typedef struct {
	float aexpn;
	float aexp0;
	float amplt;
	float astep;
	int   istep;
	float partw;
	float tintg;
	float ekin;
	float ekin1;
	float ekin2;
	float au0;
	float aeu0;
	int   Nrow;
	int   Ngrid;
	int   Nspecies;
	int   Nseed;
	float Om0;
	float Oml0;
	float hubble;
	float Wp5;
	float Ocurv;
	float Omb0;  
	float mass[10];
	unsigned int num[10];
	float fill[80];
} particle_header;

typedef struct {
	float aexpn;
	float aexp0;
	float amplt;
	float astep;
	int   istep;
	float partw;
	float tintg;
	float ekin;
	float ekin1;
	float ekin2;
	float au0;
	float aeu0;
	int   Nrow;
	int   Ngrid;
	int   Nspecies;
	int   Nseed;
	float Om0;
	float Oml0;
	float hubble;
	float Wp5;
	float Ocurv;
	float mass[10];
	unsigned int   num[10];
	float fill[80];
} nbody_particle_header;

typedef double particle_float;

typedef struct {
	double x[nDim];
	double v[nDim];
	float mass;
	int is_star;
#ifdef STARFORM
	float initial_mass;
	float star_tbirth;
	float star_metallicity_II;
	float star_metallicity_Ia;
#endif /* STARFORM */

	int subhalo_flag;
} particle_struct;

void read_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag );
void read_particles( char *header_filename, char *data_filename, char *stellar_filename,
		halo_list *, void callback( halo_struct *, particle_struct * ) );

void read_indexed_particles( char *header_filename, char *data_filename, char *stellar_filename,
		halo_list *, halo_list *, void callback( halo_struct *, particle_struct * ) );


#endif /* PARTICLES */

#endif
