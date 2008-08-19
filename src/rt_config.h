#ifndef __RT_CONFIG_H__
#define __RT_CONFIG_H__

/* 
// rt_defs.h file should only be included here 
*/
#define __RT_INCLUDING_DEFS
#include "rt_defs.h"
#undef __RT_INCLUDING_DEFS


#ifndef RT_TRANSFER

/*
//  There is no sense in filling in the tables in the OT mode
//  (in principle, this is independent of using RT, but
//  the mode in which this is useful is not implemented).
*/
#ifdef RT_TABLES
#undef RT_TABLES
#endif

#endif


#ifndef RT_CHEMISTRY

/*
//  There is no sense in using non-monoatomic thermodynamics
//  if chemistry is not included
*/
#ifndef RT_MONOATOMIC
#define RT_MONOATOMIC
#endif

/*
//  There is no sense in using high-density mode
//  if chemistry is not included
*/
#ifdef RT_HIGH_DENSITY
#undef RT_HIGH_DENSITY
#endif

/*
//  There is no sense in using LW bands
//  if chemistry is not included
*/
#ifdef RT_LWBANDS
#undef RT_LWBANDS
#endif

#endif


/*
//  Apply flux-conserving correction a-la Abel to photoionization rates.
//  There is no reason why it shouldn't be always on.
*/
#define RT_TRANSFER_FLUX_CONSERVING


/*
//  Use log interpolation for photo rates. 
*/
#define RT_INTERPOLLOG


/*
//  If the required photo rate is not found in the rate table, use direct
//  integration to compute it (which would be VERY SLOW). If unset, assume
//  that the table is wide enough (parameter acOmax from rt_tables.h is
//  large enough) for the source spectrum and set the rate to zero.
*/
#define RT_NARROWTABLE


/*
//  Allow for non-equilibrium abundances for H_2^+ and H^-. Tom Abel claimed
//  that those two are always in the equilibrium. This is introduced for 
//  testing purposes only, so it is hidden here.
*/
#undef RT_8SPECIES


/*
//  H2 cooling and formation/destruction rates
//  H2_RATE: 1 = Shapiro & Kang 
//           2 = Abel et all
//           3 = Galli & Palla
//  H2_COOL: 1 = Lepp & Shull 
//           2 = Galli & Palla
*/
#define RT_H2_RATE 3
#define RT_H2_COOL 2


/*
//  If we are running several specific tests, use a single source
*/
#if defined(RT_TEST) && (RT_TEST==1 || RT_TEST==5 || RT_TEST==6)
#define RT_SINGLE_SOURCE
#endif


#endif  /* __RT_CONFIG_H__ */
