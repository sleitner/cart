#ifndef __IO_CART_H__
#define __IO_CART_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

void config_init_io_cart();
void config_verify_io_cart();
void reorder( char *buffer, int size );

#ifdef PARTICLES

typedef struct {
	float aunin;
	float auni0;
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
	float OmM0;
	float OmL0;
	float h100;
	float Wp5;
	float OmK0;
	float OmB0;  
	float mass[10];
	unsigned int   num[10];
	float magic1;
	float DelDC;
	float abox;   /* Scale factor in the box */
	float Hbox;   /* Hubble constant in the box */
	float magic2;
	float fill[75];
} particle_header;

typedef struct {
	float aunin;
	float auni0;
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
	float OmM0;
	float OmL0;
	float h100;
	float Wp5;
	float OmK0;
	float mass[10];
	unsigned int   num[10];
	float magic1;
	float DelDC;
	float abox;   /* Scale factor in the box */
	float Hbox;   /* Hubble constant in the box */
	float magic2;
	float fill[75];
} nbody_particle_header;

void write_cart_particles( char *header_filename, char *data_filename, char *timestep_filename, char *stellar_filename );
void read_cart_particles(  char *header_filename, char *data_filename, char *timestep_filename, char *stellar_filename, 
			int num_sfcs, int *sfc_list );
void read_cart_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag ); 
#endif /* PARTICLES */

void read_cart_restart(const char *);
void write_cart_restart( int grid_filename_flag, int particle_filename_flag, int tracer_filename_flag );

void restart_load_balance_cart( char *grid_filename, char *particle_header_filename, char *particle_data );
void write_cart_grid_binary( char *filename );
void read_cart_grid_binary( char *filename );

#ifdef HYDRO_TRACERS
void read_cart_hydro_tracers( char *filename );
void write_cart_hydro_tracers( char *filename );
#endif /* HYDRO_TRACERS */

/*
//  This is mostly for internal use. Setting a non-zero mode allows to read in 
//  a non-native file. The mode is automatically reset to 0 after each file 
//  read. 
//  Allowed modes for a grid file:
//     1. Read the file from a run with RT enabled.
//     2. Read the file from a run with RT+RT_UV enabled.
//  Allowed modes for particle files:
//     1. Read files with double-precision positions but single-precision times.
//     2. Read files with single-precision positions and times.
//  Other setting in the other run must be the same as in the current one (i.e.
//  no reading a cooing data file into an adiabatic run).
*/
void set_cart_grid_file_mode(int mode);
void set_cart_particle_file_mode(int mode);

extern int cart_particle_num_row;
extern int num_cart_output_files;
extern int num_cart_input_files;

#endif
