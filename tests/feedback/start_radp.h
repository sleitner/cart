#ifndef __STARTRADP_H__
#define __STARTRADP_H__


#define slice_axis_z    2  
#define slice_hsize_pc  2000  //(->2000)

#define omm0 1.0
#define oml0 0.0
#define omb0 1.0
#define hubble 1.0
#define deltadc 0.0
#define a0 0.9

#define boxh (20.0e-3/a0*hubble) //20,000pc/2^9~2000/51=40pc  (4->8pc)


#define box_traverse_time  (40e6)
//#define box_traverse_time  (0)
#define equil 1.0
#define n_h2            (100.0)
#define n_ambient       (n_h2/equil) 
#define T_h2            (100.0)  
#define T_ambient       (T_h2*equil)
#define mstar_one_msun   1e5  //5e3; //10



double advection_momentum(int icell);
extern double tot_momentum0;
extern float adv_velocity[3];
extern double tot_energy0;
extern int last_star_id;
extern int central_cell;
int icell_central(double dispx,double dispy,double dispz);
    
void radial_average( int cell, int level );


#ifndef STARFORM
#error need stars defined for feedback tests
#endif

#endif
