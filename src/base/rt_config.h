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
#define RT_BACKGROUND_HAARDT_MADAU             1
#define RT_BACKGROUND_SELFCONSISTENT           2
#define RT_BACKGROUND_SELFCONSISTENT_AND_QSO   3


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
#ifdef RT_EXACT_EOS
#undef RT_EXACT_EOS
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

/*
//  There is no sense in using dust evolution
//  if chemistry is not included
*/
#ifdef RT_DUST_EVOLUTION
#undef RT_DUST_EVOLUTION
#endif

#endif


/*
//  Use log interpolation for photo rates. I cannot imagine it
//  should ever be off, hence it is here.
*/
#define RT_INTERPOLLOG


/*
//  Allow for non-equilibrium abundances for H_2^+ and H^-. Tom Abel claimed
//  that those two are always in the equilibrium. This is introduced for 
//  testing purposes only, so it is hidden here (in fact, it is not fully
//  implemented in the C part.
*/
/* #define RT_8SPECIES */


/*
//  H2 cooling and formation/destruction rates
//  H2_RATE: 0 = Glover & Abel
//           1 = Shapiro & Kang / Lepp & Shull
//           2 = Galli & Palla
*/
#ifndef RT_H2_RATE
#define RT_H2_RATE 0
#endif


/*
//  Dust absorption cross-section
//  0 = LMC-like dust
//  1 = SMC-like dust
//  (if underfined, defaults to LMC-like)
*/
#ifndef RT_DUST_CS
#define RT_DUST_CS 0
#endif


/*
//  Mode of ionization clumping factor:
//  0:  no clumping factors
//  1:  use clumping factors in external RF only
//  2:  use clumping factors in full RF 
//  (if underfined, defaults to no clumping)
*/
#ifndef RT_CFI
#define RT_CFI 0
#endif

/*
//  A helper switch to optimize calculations if photoionization and
//  photoheating rates are const during one cooling step. Should be set
//  only with RT_TRANSFER and RT_TRANSFER_FLUX_CONSERVING are both on.
//  (May be removed after the development is complete.)
*/
#if defined(RT_TRANSFER) && defined(RT_TRANSFER_FLUX_CONSERVING)
#define RT_VARIABLE_RF
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
//  RT needs COOLING to be set 
*/
#ifndef COOLING
#define COOLING
#endif

/*
//  There is no sense in using dust evolution
//  if the main code not using it
*/
#if !defined(DUST_EVOLUTION) && defined(RT_DUST_EVOLUTION)
#undef RT_DUST_EVOLUTION
#endif

#endif /* __RT_CONFIG_H__ */
