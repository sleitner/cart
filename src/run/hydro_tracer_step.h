#ifndef __HYDRO_TRACER_STEP_H__
#define __HYDRO_TRACER_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#if defined(HYDRO) && defined(HYDRO_TRACERS)

void move_hydro_tracers( int level );

#endif /* defined(HYDRO) && defined(HYDRO_TRACERS) */

#endif
