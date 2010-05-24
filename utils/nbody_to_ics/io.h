#ifndef __IO_H__
#define __IO_H__

#define nDim		3
#define MAX_OUTPUTS	256
#define MAX_PARTICLE_SPECIES	10

extern char output_directory[256];
extern char logfile_directory[256];
extern char jobname[256];
extern int num_output_files;

extern double particle_species_mass[10];
extern int particle_species_indices[11];
extern int particle_species_num[10];
extern int num_particle_species;

void reorder( char *buffer, int size );

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
	unsigned int   num[10];
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
	particle_float x[nDim];
	particle_float v[nDim];
	float mass;
	int specie;
} particle_struct;

void read_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag );
void read_particles( char *header_filename, char *data_filename, void callback( particle_struct * ) );

#endif
