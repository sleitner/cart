#ifndef __RT_TREE_H__
#define __RT_TREE_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

/*
//  Some numbers that depend on the settings in rt_defs.h 
*/
#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


#ifdef RT_TRANSFER

/*
//  Currently implemented framework for radiative transfer:
//  two-field ansatz
*/
#define rt_num_near_fields_per_freq      1

#ifdef RT_UV
#define rt_num_freqs      4
#else
#define rt_num_freqs      3
#endif

/*
//  Currently implemented methods for radiative transfer
*/
#if (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

#define rt_num_et_vars          (nDim*(nDim+1)/2)
#define rt_num_extra_vars       (2+rt_num_et_vars)

#define RT_VAR_OT_FIELD 	(rt_grav_vars_offset+1)
#define rt_et_offset            (rt_grav_vars_offset+2)

#endif /* RT_TRANSFER_METHOD == RT_METHOD_OTVET */

#if !defined(rt_num_extra_vars)
#error "A radiative transfer method must be specified.\n Please set RT_TRANSFER_METHOD to one of the supported methods (listed in rt_config.h)."
#endif

#define rt_field_offset         (rt_grav_vars_offset+rt_num_extra_vars)
#define rt_far_field_offset     (rt_field_offset+rt_num_near_fields_per_freq*rt_num_freqs)

#define rt_num_fields_per_freq  (1+rt_num_near_fields_per_freq)
#define rt_num_fields           (rt_num_fields_per_freq*rt_num_freqs)
#define rt_num_vars             (rt_num_extra_vars+rt_num_fields)
#define rt_num_disk_vars        rt_num_fields
#define rt_disk_offset          rt_field_offset

#else  /* RT_TRANSFER */

#define rt_num_vars             0
#define rt_num_disk_vars        0
#define rt_disk_offset          rt_grav_vars_offset

#endif /* RT_TRANSFER */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_TREE_H__ */
