#ifndef __IO_H__
#define __IO_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

extern const char* output_directory;
extern const char* logfile_directory;
extern const char* jobname;
extern const char* requeue_command;
extern int current_output;
extern int last_restart_step;

extern int current_restart_backup;
extern int num_restart_backups;

extern int output_frequency;

extern int num_outputs;
extern float *outputs;

void config_init_io();
void config_verify_io();

void read_restart( const char *label );
void set_jobname( const char *str );

#define NO_WRITE		0
#define WRITE_BACKUP	1
#define WRITE_SAVE      2                                                                                                         

#define WRITE_GRID		1
#define WRITE_PARTICLES	2
#define WRITE_TRACERS	4
#endif
