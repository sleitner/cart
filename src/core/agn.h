#ifndef __AGN_H__
#define __AGN_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef AGN

#define SINK_ACCRETION_RADIUS   (3.0) /* in units of accretion kernel scale r_K */
#define SINK_RADIUS_SAMPLES     (3)   /* number of sampling points over accretion radius 
                                       * (this is what we'd change to 0 if we wanted to do single cell accretion)  
                                       * normally set to 3 */  
#define SINK_SAMPLES_1D         (2*SINK_RADIUS_SAMPLES+1)                         /* total number of sampling points in 1d */
#define SINK_SAMPLES_3D         (SINK_SAMPLES_1D*SINK_SAMPLES_1D*SINK_SAMPLES_1D) /* number of sample points per accretion radius */

void config_init_agn();
void config_verify_agn();
void init_agn();

#endif /* AGN */

#endif 
