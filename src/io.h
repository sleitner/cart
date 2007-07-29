#ifndef __IO_H__
#define __IO_H__

#define MAX_OUTPUTS	256

extern char output_directory[256];
extern char logfile_directory[256];
extern char jobname[256];
extern char requeue_command[256];
extern float outputs[MAX_OUTPUTS];
extern int num_outputs;
extern int current_output;
extern int last_restart_step;

extern int num_output_files;

void reorder( char *buffer, int size );

void save_check();

#define WRITE_GENERIC		0
#define WRITE_SAVE		1

void write_restart( int gas_filename_flag, int particle_filename_flag, int tracer_filename_flag );
void read_restart( double aexpn );

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

void write_particles( char *header_filename, char *data_filename, char *timestep_filename, char *stellar_filename );
void read_particles(  char *header_filename, char *data_filename, char *timestep_filename, char *stellar_filename, int num_sfcs, int *sfc_list );
void read_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag ); 
#endif /* PARTICLES */

#ifdef HYDRO
void read_gas_ic( char *filename );
void write_grid_binary( char *filename );
void write_hart_gas_binary( char *filename );
void read_grid_binary( char *filename );
void read_indexed_grid( char *filename, int num_sfcs, int *sfc_list, int max_level_to_read );

#ifdef HYDRO_TRACERS
void read_hydro_tracers( char *filename );
void write_hydro_tracers( char *filename );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */

#endif
