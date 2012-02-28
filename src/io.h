#ifndef __IO_H__
#define __IO_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

extern char output_directory[];
extern char logfile_directory[];
extern char jobname[];
extern char requeue_command[];
extern int current_output;
extern int last_restart_step;

extern int output_frequency;

extern int num_outputs;
extern float *outputs;

void config_init_io();
void config_verify_io();

void write_restart( int grid_filename_flag, int particle_filename_flag, int tracer_filename_flag );
void read_restart( const char *label );
void save_check();

#define NO_WRITE		0
#define WRITE_GENERIC   1
#define WRITE_SAVE      2                                                                                                         
#define WRITE_BACKUP    3

#define WRITE_GRID		1
#define WRITE_PARTICLES	2
#endif
