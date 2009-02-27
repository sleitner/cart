#ifdef __RT_INCLUDING_DEFS


#define RT_DEBUG


/*
//  Radiative transfer switches: only have effect when RADIATIVE_TRANSFER
//  switch in defs.h is switched on.
//  ---------------------------------------------------------------------
*/

/*
//  Produce diagnostic output 
*/
#define RT_OUTPUT


/*
//  Use tables for computing photoionization rates. 
*/
#define RT_TABLES


/*
//  Allow for spatially inhomogeneous radiation field.
*/
#define RT_TRANSFER


/*
//   Set the method for radiative transfer:
//   1.  OTVET approximation
*/
#define RT_TRANSFER_METHOD RT_METHOD_OTVET


/*
//  Allow for H2 checmistry & metal cooling. 
*/
#define RT_CHEMISTRY 


/*
//  Assume that gas behaves as monoatomic (i.e c_v is always 3/2).
//  This is valid if H2 fraction is always small (optimization switch).
*/
#define RT_MONOATOMIC


/*
//  Include additional processes if the gas density can be high:
//  H2 ortho-to-para computation, three-body formation of H2,
//  LTE cooling of H2, etc. The critical density is about 1e4 per cc,
//  although it depends somewhat on the metallicity and radiation field.
//  This is rarely needed for galaxy formation simulations.
*/
/* #define RT_HIGH_DENSITY */


/*
//  Allow for secondary electrons from X-rays. 
*/
/* #define RT_XRAYS */


/*
//  Allow for absorption by dust (disables the use of tables for comuting
//  photoionozation and photoheating rates, so the code would be VERY SLOW).
*/
/* #define RT_DUST */


/*
//  Allow for line transfer in Lyman-Werner bands 
*/
#define RT_LWBANDS


/*
//  Allow for heating by recoil in Lyman-alpha line using Tozzi et al formula.
//  Jordi Miralda-Escude claims that the formula is incorrect.
*/
/* #define RT_LYMAN_ALPHA_HEATING */


/*
//  Include PAH and cosmic ray ionizations and heating
//  (normally, not important)
*/
/* #define RT_PAH_CR */


/*
//  Keep the signal propagation speed equal c. 
*/
/* #define RT_SIGNALSPEED_TO_C */


/*
//  Set the floor for the gast-to-dust ratio. May be needed to force the
//  switch from primordial, metal-free star formation to normal star
//  formation if the resolution is not high enough to resolve first stars
//  properly.
*/
#define RT_DUST_TO_GAS_FLOOR 0.001


/*
//  Enable computing the evolving cosmic background
*/
/* #define RT_EXTERNAL_BACKGROUND RT_BACKGROUND_SELFCONSISTENT */


/*
//  Mode of ionization clumping factor:
//  0:  no clumping factors
//  1:  use clumping factors in external RF only
//  2:  use clumping factors in full RF 
*/
#define RT_CFI 0

/*
//  The number of OpenMP buffers to use for shared-memory-parallel 
//  accumulation of global quantities. It needs to be at least
//  the number of OpenMP threads you are using, but not too large to
//  avoid wasting storage (buffers are not that small).
*/
#define RT_PARALLEL_NUM_OPENMP_BUFFERS 8

/*
//  Define this if the program uses MPI
*/
#define RT_PARALLEL_USE_MPI

#else

#error "rt_def.h file should never be included directly"

#endif  /* __RT_INCLUDING_DEFS */