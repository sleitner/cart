#ifndef __RT_TREE_H__
#define __RT_TREE_H__


#include "rt_config.h"


/*
//  Some numbers that depend on the settings in rt_defs.h 
*/
#ifdef RT_TRANSFER

/*
//  Currently implemented methods for radiative transfer
*/
#define RT_METHOD_OTVET    1

#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

#define rt_num_frequencies      6

#define rt_num_et_vars          (nDim*(nDim+1)/2)
#define rt_num_vars             (2+rt_num_et_vars+rt_num_frequencies)
#define RT_VAR_OT_FIELD 	(rt_grav_vars_offset+1)
#define rt_et_offset            (rt_grav_vars_offset+2)
#define rt_freq_offset          (rt_grav_vars_offset+2+rt_num_et_vars)

#else
#error "There is no RT method beyond OTVET implemented yet."
#endif

#define rt_num_disk_vars        rt_num_frequencies
#define rt_disk_offset          rt_freq_offset

#else

#define rt_num_vars             0
#define rt_num_disk_vars        0
#define rt_disk_offset          rt_grav_vars_offset

#endif


/*
//  RT needs both COOLING and ADVECT_SPECIES to be set 
*/
#ifndef COOLING
#define COOLING
#endif

#ifndef ADVECT_SPECIES
#define ADVECT_SPECIES
#endif

#endif  /* __RT_TREE_H__ */
