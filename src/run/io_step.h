#ifndef __IO_STEP_H__
#define __IO_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

void write_restart( int grid_filename_flag, int particle_filename_flag, int tracer_filename_flag );
void save_check();

#endif
