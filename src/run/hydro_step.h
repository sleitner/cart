#ifndef __HYDRO_STEP_H__
#define __HYDRO_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef HYDRO

#define HYDRO_COPY_ALL      0
#define HYDRO_RESTORE_CLEAN 1
#define HYDRO_RESTORE_ALL   2

DECLARE_LEVEL_ARRAY(int,level_sweep_dir);

void hydro_step_init();
void hydro_step( int level );
void hydro_copy_vars( int level, int direction );

#ifdef ISOTROPIC_TURBULENCE_ENERGY
void hydro_turbulence_sources(int level);
#endif /* ISOTROPIC_TURBULENCE_ENERGY */
void hydro_eos(int level);

#endif /* HYDRO */
#endif /* __HYDRO_STEP_H__ */
