#define RT_OUTPUT
#define RT_OTVET_CACHE_ET
#define RT_OTVET_CACHE_RF

#define RT_OTVET_SAVE_FLUX 
/* //  Allow for spatially inhomogeneous radiation field. */
#define RT_TRANSFER
#define RT_TRANSFER_METHOD RT_METHOD_OTVET

/* //  Allow for H2 checmistry & metal cooling.  */
#define RT_CHEMISTRY
/* //  By default, if RT_CHEMISTRY is not set then the LW rate is not computed; */
/* //  if it is actually needed, use this switch. */
#define RT_COMPUTE_LW_RATE
/* //  Also include transfer of UV (non-ionizing) radiation. That requires */
/* //  RT_CHEMISTRY and adds two extra variables for each cell. */
#define RT_UV
/* //  Allow for line transfer in Lyman-Werner bands  */
#define RT_LWBANDS


/* //  Enable the external cosmic background, the valid values are */
/* //  1. Self-consistent from the sources inside the box (RT_BACKGROUND_SELFCONSISTENT) */
/* //  2. Haardt-Madau 2001 (RT_BACKGROUND_HAARDT_MADAU) */
/* #define RT_EXTERNAL_BACKGROUND RT_BACKGROUND_HAARDT_MADAU   */
/* #define RT_EXTERNAL_BACKGROUND RT_BACKGROUND_SELFCONSISTENT */



/* //  The number of OpenMP buffers to use for shared-memory-parallel  */
/* //  accumulation of global quantities. It needs to be at least */
/* //  the number of OpenMP threads you are using, but not too large to */
/* //  avoid wasting storage (buffers are not that small). */
#define RT_PARALLEL_NUM_OPENMP_BUFFERS 8
/* //  Define this if the program uses MPI */
/* #define RT_PARALLEL_USE_MPI   */

