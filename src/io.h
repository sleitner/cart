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
#define WRITE_SAVE			1
#define WRITE_BACKUP		2

void write_restart( int gas_filename_flag, int particle_filename_flag, int tracer_filename_flag );
void read_restart( double aexpn );
void restart_load_balance( char *grid_filename, char *particle_header_filename, char *particle_data );

#ifdef PARTICLES
#define PARTICLE_HEADER_MAGIC		(0.1234f)

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

void write_particles( char *header_filename, char *data_filename, char *timestep_filename, char *stellar_filename );
void read_particles(  char *header_filename, char *data_filename, char *timestep_filename, char *stellar_filename, int num_sfcs, int *sfc_list );
void read_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag ); 
#endif /* PARTICLES */

void write_grid_binary( char *filename );
void read_grid_binary( char *filename );

#ifdef HYDRO
void read_gas_ic( char *filename );
void read_indexed_grid( char *filename, int num_sfcs, int *sfc_list, int max_level_to_read );

#ifdef HYDRO_TRACERS
void read_hydro_tracers( char *filename );
void write_hydro_tracers( char *filename );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */

#endif
