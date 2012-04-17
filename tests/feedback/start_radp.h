#ifndef __STARTRADP_H__
#define __STARTRADP_H__


#define omm0 1.0
#define oml0 0.0
#define omb0 1.0
#define hubble 1.0
#define deltadc 0.0
#define a0 0.9

#define NSTARS_1D   (1)
#define NSTARS (NSTARS_1D*NSTARS_1D*NSTARS_1D)

/* #define boxh (20.0e-3/a0*hubble) //20,000pc/2^9~2000/51=40pc  (4->8pc) */
/* #define box_traverse_time  (40e6)  */
#define box_traverse_time  (0) 

#define boxh (5.0e-3/a0*hubble/pow(2.,NSTARS_1D-1)) //20,000pc/2^9~2000/51=40pc  (4->8pc) 
/* #define box_traverse_time  (40e6/pow(2.,NSTARS_1D-1)) */



#define equil           (1e-7) 
/* #define equil           (1.0)  */
#define n_h2            (1e-5)
#define T_h2            (1.0e4)  
#define n_ambient       (n_h2/equil) 
#define T_ambient       (T_h2*equil)
//#define mstar_one_msun   1e5  //5e3; //10
#define mstar_one_msun   1e6  //5e3; //10
#define Radius_stargas   1000.  //5e3; //10



double advection_momentum(int icell,int idir);
extern double tot_momentum0;
//extern const int NSTARS;
extern float adv_velocity[3];
extern double tot_energy0;
extern int last_star_id;

int icell_central(double dispx,double dispy,double dispz);
    
void radial_average( int cell, int level );


#ifndef STARFORM
#error need stars defined for feedback tests
#endif

#endif
