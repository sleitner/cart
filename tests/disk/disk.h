#ifndef __STARTDISK_H__
#define __STARTDISK_H__

#ifndef STAR_FORMATION
#error need stars defined for feedback tests
#endif

/*
// fill in the mass of the missing particles with a uniform DM component
// #define UNIFORM_PARTICLES 
*/

#define nomore_starformation  0
/* #warning "no more star formation" */
/* #define nomore_starformation 1 */
#define n_ambient      (1e-5)
#define T_ambient      (1e6)
#define Tgasmodelat1cc  1.0e2
#define nbuild_mesh     2


int icell_wrt_central(double dispx,double dispy,double dispz);
void radial_average( int cell, int level );

#define omm0 0.30
#define oml0 0.7
#define omb0 0.05
#define hubble 0.7
#define deltadc 0.0

#define N_GAS   0 
#define N_HALO  0
#define N_STARS  0
#define NUM_VCIRC_POINTS 10000

#endif /*__STARTDISK_H__*/
