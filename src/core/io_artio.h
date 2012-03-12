#ifndef __IO_ARTIO_H__
#define __IO_ARTIO_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

void config_init_io_artio();
void config_verify_io_artio();

void artio_restart_load_balance( artio_file handle );
void read_artio_restart( const char *label );
void write_artio_restart( int grid_filename_flag, int particle_filename_flag, int tracer_filename_flag );

extern int num_artio_grid_files;
#ifdef PARTICLES
extern int num_artio_particle_files;
#endif /* PARTICLES */

#endif /* __IO_ARTIO_H__ */
