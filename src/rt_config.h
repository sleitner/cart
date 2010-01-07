#ifndef __RT_CONFIG_H__
#define __RT_CONFIG_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


/*
//  Currently implemented methods for radiative transfer
*/
#define RT_METHOD_OTVET     1  /* Optically Thin Variable Eddington Tensor approximation */


/*
//  Currently implemented external backgrounds
*/
#define RT_BACKGROUND_SELFCONSISTENT         1
#define RT_BACKGROUND_HAARDT_MADAU           2


/* 
// rt_defs.h file should only be included here 
*/
#include "rt_defs.h"
#define RT_CONFIGURED


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
/* #define RT_8SPECIES */

/*
//  Apply flux-conserving correction a-la Abel to photoionization rates.
//  May slow down the cooling computation. So far, I found no effect, so it
//  is off by default
*/
/* #define RT_TRANSFER_FLUX_CONSERVING */


/*
//  H2 cooling and formation/destruction rates
//  H2_RATE: 0 = Glover & Abel
//           1 = Shapiro & Kang / Lepp & Shull
//           2 = Galli & Palla
*/
#define RT_H2_RATE 0

/*
//  By defult use the full chemical model -  it is the only one
//  that works in all regimes
*/
#define RT_CHEMISTRY_FULL_MODEL

/*
//  A helper switch to optimize calculations if photoionization and
//  photoheating rates are const during one cooling step. Should be set
//  only with RT_TRANSFER and RT_TRANSFER_FLUX_CONSERVING are both on.
//  (May be removed after the development is complete.)
*/
#if defined(RT_TRANSFER) && defined(RT_TRANSFER_FLUX_CONSERVING)
/* #define RT_VARIABLE_PRATES */
/* #define RT_VARIABLE_RFIELD */
#endif

/*
//  We need at least 1 global buffer
*/
#if (!defined(RT_PARALLEL_NUM_OPENMP_BUFFERS) || RT_PARALLEL_NUM_OPENMP_BUFFERS<1)
#define RT_PARALLEL_NUM_OPENMP_BUFFERS 1
#endif

/*
//  If we are running several specific tests, use a single source
*/
#if defined(RT_TEST) && (RT_TEST==1 || RT_TEST==5 || RT_TEST==6)
#define RT_SINGLE_SOURCE
#endif

/*
//  If we are running several specific tests, remove the background
*/
#if defined(RT_TEST) && (RT_TEST==1 || RT_TEST==5 || RT_TEST==6 || RT_TEST==11 || RT_TEST==15 || RT_TEST==16)
#ifdef RT_EXTERNAL_BACKGROUND
#undef RT_EXTERNAL_BACKGROUND
#endif
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


#endif /* __RT_CONFIG_H__ */
