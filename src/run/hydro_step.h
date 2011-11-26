#ifndef __HYDRO_STEP_H__
#define __HYDRO_STEP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef HYDRO

#define COPY		0
#define RESTORE		1
#define COPY_ZERO_REF	2

#define	COPY_ALL_LEAFS		0
#define	COPY_SPLIT_NEIGHBORS	1
#define COPY_NO_SPLIT_NEIGHBORS	2

DECLARE_LEVEL_ARRAY(int,level_sweep_dir);

void hydro_step_init();
void hydro_step( int level );
void hydro_copy_vars( int level, int direction, int cell_type );

void hydro_eos(int level);
void hydro_magic(int level);

#endif /* HYDRO */
#endif /* __HYDRO_STEP_H__ */
