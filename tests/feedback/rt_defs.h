
#define RT_OLDSTYLE_SOURCE_FUNCTION
#define RT_OUTPUT


#define RT_TRANSFER
#define RT_TRANSFER_METHOD RT_METHOD_OTVET
#define RT_CHEMISTRY 
/*#define RT_LWBANDS*/ /*DE got rid for Z=1,U=100 run*/
#define RT_SIGNALSPEED_TO_C 
#define RT_EXTERNAL_BACKGROUND RT_BACKGROUND_HAARDT_MADAU
/*#define RT_EXTERNAL_BACKGROUND RT_BACKGROUND_SELFCONSISTENT/* /* This is the other extermal background allowed*/
/*#define RT_UV*/ /*Comment this out for the old runs (Milkyway and original 6MPC)*/

#define RT_PARALLEL_NUM_OPENMP_BUFFERS 8
#define RT_PARALLEL_USE_MPI

