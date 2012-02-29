#ifndef __IO_CARTIO_H__
#define __IO_CARTIO_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

void config_init_io_cartio();
void config_verify_io_cartio();

void read_cartio_restart( const char *label );
void write_cartio_restart( int grid_filename_flag, int particle_filename_flag, int tracer_filename_flag );

extern int num_cartio_grid_files;
#ifdef PARTICLES
extern int num_cartio_particle_files;
#endif /* PARTICLES */

#endif /* __IO_CARTIO_H__ */
