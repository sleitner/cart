#ifndef __PARTICLE_IO_H__
#define __PARTICLE_IO_H__

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
	int   num[10];
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
	int   num[10];
	float fill[80];
} nbody_particle_header;

void reorder( char *buffer, int size );
void read_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag ); 
void read_write_particle_header( char *header_filename, char *out_filename, particle_header *header, int *endian, int *nbody_flag );

#endif
