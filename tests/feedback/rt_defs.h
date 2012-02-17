/*#define RT_OTVET_SAVE_FLUX  */

/* //  Radiative transfer switches: only have effect when RADIATIVE_TRANSFER */
/* //  switch in defs.h is switched on. */
/* //  --------------------------------------------------------------------- */


/* //  Allow for spatially inhomogeneous radiation field. */
#define RT_TRANSFER



/* //   Set the method for radiative transfer: */
/* //   1.  OTVET approximation (RT_METHOD_OTVET) */
#define RT_TRANSFER_METHOD RT_METHOD_OTVET



/* //  Allow for H2 checmistry & metal cooling.  */
#define RT_CHEMISTRY



/* //  By default, if RT_CHEMISTRY is not set then the LW rate is not computed; */
/* //  if it is actually needed, use this switch. */
#define RT_COMPUTE_LW_RATE



/* //  Also include transfer of UV (non-ionizing) radiation. That requires */
/* //  RT_CHEMISTRY and adds two extra variables for each cell. */
#define RT_UV



/* //  Include additional processes if the gas density can be high: */
/* //  H2 ortho-to-para computation, three-body formation of H2, */
/* //  LTE cooling of H2, etc. The critical density is about 1e4 per cc, */
/* //  although it depends somewhat on the metallicity and radiation field. */
/* //  This is rarely needed for galaxy formation simulations. */
#define RT_HIGH_DENSITY



/* //  Allow for secondary electrons from X-rays.  */
#define RT_XRAYS



/* //  Allow for line transfer in Lyman-Werner bands  */
#define RT_LWBANDS



/* //  Allow for heating by recoil in Lyman-alpha line using Tozzi et al formula. */
/* //  Jordi Miralda-Escude claims that the formula is incorrect. */
#define RT_LYMAN_ALPHA_HEATING



/* //  Include PAH and cosmic ray ionizations and heating */
/* //  (normally, not important) */
#define RT_PAH_CR



/* //  Enable the external cosmic background, the valid values are */
/* //  1. Self-consistent from the sources inside the box (RT_BACKGROUND_SELFCONSISTENT) */
/* //  2. Haardt-Madau 2001 (RT_BACKGROUND_HAARDT_MADAU) */
#define RT_EXTERNAL_BACKGROUND RT_BACKGROUND_HAARDT_MADAU  
/* #define RT_EXTERNAL_BACKGROUND RT_BACKGROUND_SELFCONSISTENT */



/* //  The number of OpenMP buffers to use for shared-memory-parallel  */
/* //  accumulation of global quantities. It needs to be at least */
/* //  the number of OpenMP threads you are using, but not too large to */
/* //  avoid wasting storage (buffers are not that small). */
#define RT_PARALLEL_NUM_OPENMP_BUFFERS 8



/* //  Define this if the program uses MPI */
#define RT_PARALLEL_USE_MPI  

