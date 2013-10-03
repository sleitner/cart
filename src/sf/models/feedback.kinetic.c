#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "imf.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "tree.h"
#include "units.h"
#include "timing.h"

#include "rand.h"
#include "feedback.kinetic.h"
#include "plugin.h"

/*  
// This file distributes input momentum through kicks, internal energy, or the Riemann solver (Pextra)
*/

#ifdef STAR_FORMATION
extern double feedback_temperature_ceiling;
#endif /* STAR_FORMATION */
extern double feedback_turbulence_temperature_ceiling;
double kfb_boost_kicks = 1;

#define UNIT_VECTOR     -1
#define CUBE_VECTOR     -2
#define NEIGHBOR_VECTOR -3

#define KFB_METHOD_PRESSURIZE       0
#define KFB_METHOD_KICKS            1
#define KFB_METHOD_HYBRID           2
#define KFB_METHOD_HYBRIDGATHER     3
int kfb_internal_method = KFB_METHOD_PRESSURIZE;

#define KFB_SPREAD_CELL       0
#define KFB_SPREAD_OCT        1
#define KFB_SPREAD_CUBE       2
int kfb_internal_spread = KFB_SPREAD_CELL;

#define KFB_TURBULENCE_NONE   0
#define KFB_TURBULENCE_CANCEL 1
#define KFB_TURBULENCE_BOTH   2
int kfb_internal_turbulence = KFB_TURBULENCE_NONE;

extern double feedback_speed_time_ceiling;
extern double dvfact;

void kfb_kick_cube(double dp, int level, int icell);

void control_parameter_set_kfb_method(const char *value, void *ptr, int ind) {
	if ( strcmp(value,"pressurize") == 0 ) {
		kfb_internal_method = KFB_METHOD_PRESSURIZE;
	} else if ( strcmp(value, "hybrid") == 0 ) {
		kfb_internal_method = KFB_METHOD_HYBRID;
	} else if ( strcmp(value, "kicks") == 0 ) {
		kfb_internal_method = KFB_METHOD_KICKS;
	} else if ( strcmp(value, "hybrid_gather") == 0 ) {
		kfb_internal_method = KFB_METHOD_HYBRIDGATHER;
	} else {
		cart_error("String '%s' is not a valid method for kfb feedback.",value);
	}
}

void control_parameter_set_kfb_spread(const char *value, void *ptr, int ind) {
    if ( strcmp(value,"cell") == 0 ) {
        kfb_internal_spread = KFB_SPREAD_CELL;
    } else if ( strcmp(value, "oct") == 0 ) {
        kfb_internal_spread = KFB_SPREAD_OCT;
    } else if ( strcmp(value, "cube") == 0 ) {
        kfb_internal_spread = KFB_SPREAD_CUBE;
    } else {
        cart_error("String '%s' is not a valid spreading method for kfb feedback.",value);
    }
}

void control_parameter_set_kfb_turbulence(const char *value, void *ptr, int ind) {
	if ( strcmp(value,"none") == 0 ) {
		kfb_internal_turbulence = KFB_TURBULENCE_NONE;
	} else if ( strcmp(value, "cancel") == 0 ) {
		kfb_internal_turbulence = KFB_TURBULENCE_CANCEL;
	} else if ( strcmp(value, "both") == 0 ) {
		kfb_internal_turbulence = KFB_TURBULENCE_BOTH;
	} else {
		cart_error("String '%s' is not a valid spreading method for kfb feedback.",value);
	}
}

void control_parameter_kfb_opts_method(FILE *stream, const void *ptr){
}
void control_parameter_kfb_opts_spread(FILE *stream, const void *ptr){
}
void control_parameter_kfb_opts_turbulence(FILE *stream, const void *ptr){
}

void kfb_config_init(){
    ControlParameterOps control_parameter_kfb_method = { control_parameter_set_kfb_method, control_parameter_kfb_opts_method };
    ControlParameterOps control_parameter_kfb_spread = { control_parameter_set_kfb_spread, control_parameter_kfb_opts_spread };
    ControlParameterOps control_parameter_kfb_turbulence = { control_parameter_set_kfb_turbulence, control_parameter_kfb_opts_turbulence };

    control_parameter_add2(control_parameter_kfb_method,&kfb_internal_method,"kfb:method","kfb_method","the method for distributing kinetic feedback. Valid names are 'pressure', 'kicks', 'hybrid'.");
    control_parameter_add2(control_parameter_kfb_spread,&kfb_internal_spread,"kfb:spread","kfb_spread","the spread of distributed kinetic feedback. Valid names are 'cell', 'oct', 'cube' (27 cells).");
    control_parameter_add2(control_parameter_kfb_turbulence, &kfb_internal_turbulence,"kfb:turbulence","kfb_to_turbulence", "the sources of turbulent pressure from kinetic feedback. Valid names are 'none' 'cancel' (when kfb cancels->turbulence before dissipating) 'both' ");

    control_parameter_add2(control_parameter_double,&kfb_boost_kicks,"kfb:boost_kicks","kfb_boost_kicks","factor boosting momentum input.");
}

void kfb_config_verify()
{

	VERIFY( kfb:method,
#ifdef EXTRA_PRESSURE_SOURCE
	        kfb_internal_method == KFB_METHOD_PRESSURIZE || kfb_internal_method == KFB_METHOD_KICKS || kfb_internal_method == KFB_METHOD_HYBRID ||
	        kfb_internal_method == KFB_METHOD_HYBRIDGATHER
#else
	        kfb_internal_method == KFB_METHOD_KICKS
#endif /* EXTRA_PRESSURE_SOURCE */
	        );

	VERIFY( kfb:spread,
	        kfb_internal_spread == KFB_SPREAD_CELL || kfb_internal_spread == KFB_SPREAD_OCT || kfb_internal_spread == KFB_SPREAD_CUBE );

	VERIFY( kfb:turbulence,
#ifdef ISOTROPIC_TURBULENCE_ENERGY
	        kfb_internal_turbulence == KFB_TURBULENCE_NONE || kfb_internal_turbulence == KFB_TURBULENCE_CANCEL || kfb_internal_turbulence == KFB_TURBULENCE_BOTH
#else
	        kfb_internal_turbulence == KFB_TURBULENCE_NONE 
#endif
	        );

	VERIFY(kfb:boost_kicks, kfb_boost_kicks>=0);
}

/* ------- auxiliary routines ----------------------------------- */
#define num_side_children 4
/* for neighbor-oct j, which of its children are closest to the original cell? */
int neighbor_side_child[num_neighbors][num_side_children] = { 
    {1,3,5,7},
    {0,2,4,6},
    {0,1,4,5},
    {2,3,6,7},
    {0,1,2,3},
    {4,5,6,7}
};

/*******************************************************
 *                27 cell stencil                      *
 *******************************************************/
/*                       6  7  8  9 10 11 12 13 14 15 16 17  */
const int CubeDir1[] = { 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 2, 3 };
const int CubeDir2[] = { 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5 };
/* in direction +/-x from x=0 y=+/- and z=+/-  */  
/*                           18  19  20  21  22  23  24  25 */
const int CubeDir3a[]    = {  0,  1,  0,  1,  0,  1,  0,  1 };
const int CubeOrigin3a[] = { 12, 12, 13, 13, 16, 16, 17, 17 };
const int CubeDir3b[]    = {  2,  3,  2,  3,  2,  3,  2,  3 };
const int CubeOrigin3b[] = { 10, 10, 11, 11, 14, 14, 15, 15 };
/*
// This array returns 26 neighbors as if the whole mesh was uniform.
// In particular, one neighbor of higher level will appear more than once. 
// The first 6 neighbors are exactly as returned by cell_all_neighbors.
// The other 20 neighbors are packed as above:
*/
int CubeDelPos[][nDim] = {
    {-1, 0, 0},    { 1, 0, 0},
    { 0,-1, 0},    { 0, 1, 0},
    { 0, 0,-1},    { 0, 0, 1},
    //6-9
    {-1,-1, 0},    { 1,-1, 0},
    {-1, 1, 0},    { 1, 1, 0},
     //10-13
    {-1, 0,-1},    { 1, 0,-1},
    { 0,-1,-1},    { 0, 1,-1},
     //14-17
    {-1, 0, 1},    { 1, 0, 1},
    { 0,-1, 1},    { 0, 1, 1},
     //18-25
    {-1,-1,-1},    { 1,-1,-1},    {-1, 1,-1},    { 1, 1,-1},
    {-1,-1, 1},    { 1,-1, 1},    {-1, 1, 1},    { 1, 1, 1}

}; 
#define CubeStencilSize 26
#define num_corners 8
/* neighbors found are either 
(1) at the same level in the current oct
(2) at the same level (as a leaf or parent) in the neighboring oct
(3) at 1 higher level (neighbor_level < level) making up the neighboring oct
if at level or level+1 add 1 unit momentum to cell or all child leafs
if at level-1 add same momentum to that level 
*/

void GetCubeStencil(int level, int cell, int nb[CubeStencilSize])
{
  int j,k, levNb[num_neighbors];
  
  /*
  //  First level neighbors
  */
  cell_all_neighbors(cell,nb);
  for(j=0; j<num_neighbors; j++) levNb[j] = cell_level(nb[j]);
  
  /*
  //  Second level neighbors
  */
  for(j=0; j<CubeStencilSize-num_neighbors-num_corners; j++)
    {
      k = num_neighbors+j;
      if(levNb[CubeDir1[j]] == level)
	{
	  /*
	  // Our neighbor in the first direction is on the same level,
	  // it is sufficient to take its neighbor in the second direction.
	  */
	  nb[k] = cell_neighbor(nb[CubeDir1[j]],CubeDir2[j]);
	}
      else if(levNb[CubeDir2[j]] == level)
	{
	  /*
	  // Our neighbor in the second direction is on the same level,
	  // it is sufficient to take its neighbor in the first direction.
	  */
	  nb[k] = cell_neighbor(nb[CubeDir2[j]],CubeDir1[j]);
	}
      else /* tried to get to new cell from 2 directions, both higher level */
	{ 
	  /*
	  // Both our neighbors are on a higher level. The 2nd-neighbor must be 
	  // kitty-corner to the current oct.
	  */
	  nb[k] = cell_neighbor(nb[CubeDir1[j]],CubeDir2[j]);
          /* only include stencil w/in 1 level of cell */
	  if( cell_level(nb[k]) < level-1 ){
              nb[k] = -1;
          }
       }
  }
  /*
  // Only do the cube corners if there is a path there at the cell level
  */
  for(j=0; j<num_corners; j++){
      k = CubeStencilSize-num_corners+j; 
      if(        cell_level( nb[CubeOrigin3a[j]] ) == level ){
          nb[k] = cell_neighbor(nb[CubeOrigin3a[j]],CubeDir3a[j]);
          nb[k] = cell_level(nb[k]) != level ? -1 : nb[k];
      } else if( cell_level( nb[CubeOrigin3b[j]] ) == level ){
          nb[k] = cell_neighbor(nb[CubeOrigin3b[j]],CubeDir3b[j]);
          nb[k] = cell_level(nb[k]) != level ? -1 : nb[k];
      } else{
          nb[k] = -1;
      }
  }
}


void cart_rand_unit_vector_block(double uni[nDim], int idir[nDim]){
    double phi, r;
    /* 
    // returns a unit vector in direction of cube stencil member:
    // first 6 cell_all_neighbors get static direction kick, 
    // next 12 get 2D randomness (one direction fixed), 
    // last 8 corners get 3D randomness like an oct (randomness of same order as non-zero entries)
    */
    if(idir[0] == 0 && idir[1] == 0){
	    uni[2] = (double)idir[2];
    }else{
	    uni[2] = (double)idir[2]*cart_rand();
    }
    r = sqrt(fabs(1-uni[2]*uni[2]));
    if     (idir[0] == 1 && idir[1] == 1){  phi = 1*M_PI/2.+cart_rand()*M_PI/2.; } //+x.+y
    else if(idir[0] ==-1 && idir[1] == 1){  phi = 2*M_PI/2.+cart_rand()*M_PI/2.; } //-x,+y
    else if(idir[0] ==-1 && idir[1] ==-1){  phi = 3*M_PI/2.+cart_rand()*M_PI/2.; } //-x,-y
    else if(idir[0] == 1 && idir[1] ==-1){  phi = 0*M_PI/2.+cart_rand()*M_PI/2.; } //+x,-y

    else if(idir[0] == 1 && idir[1] == 0 ){  phi = 1*M_PI/2.                    ; }
    else if(idir[0] == 0 && idir[1] == 1 ){  phi = 2*M_PI/2.                    ; }
    else if(idir[0] ==-1 && idir[1] == 0 ){  phi = 3*M_PI/2.                    ; } 
    else if(idir[0] == 0 && idir[1] ==-1 ){  phi = 0*M_PI/2.                    ; } //-x,+y
    
    uni[0] = r * sin(phi);
    uni[1] = -r * cos(phi);
}

void cart_rand_unit_vector_oct(double uni[nDim], int ichild){
    double phi, r;
    /* returns a random unit in direction of oct child */
    /* children start at bottom left */
    /* 4 5   6 7 */
    /* 0 1   2 3   +x-> +z^ +y>> */
    if(ichild > 3){
		uni[2] = cart_rand();
    }else{
		uni[2] = -cart_rand();
    }
    r = sqrt(1-uni[2]*uni[2]);
    if(ichild == 3 || ichild == 7){  phi = 1*M_PI/2.+cart_rand()*M_PI/2.; } //+x,+y
    if(ichild == 2 || ichild == 6){  phi = 2*M_PI/2.+cart_rand()*M_PI/2.; } //-x,+y
    if(ichild == 0 || ichild == 4){  phi = 3*M_PI/2.+cart_rand()*M_PI/2.; } //-x,-y
    if(ichild == 1 || ichild == 5){  phi = 0*M_PI/2.+cart_rand()*M_PI/2.; } //+x.-y
    
    uni[0] = r * sin(phi);
    uni[1] = -r * cos(phi);
}
/* ------- end auxiliary routines ----------------------------------- */




/* 
// How to kick individual cells 
*/
void kinetic_to_internal(int icell, double p0, double p1){
    /* cancel momentum and convert to internal energy */
    double  p2_cancel, dU;
    if(p0*p1 >0){ return; } /* same direction */
    p2_cancel = MIN( p0*p0, p1*p1 );
    
	dU = p2_cancel / cell_gas_density(icell) ; /* \Delta_ie = 2*(dp)^2/2m (per volume) */

	if ( kfb_internal_turbulence == KFB_TURBULENCE_NONE ) {	
#ifdef STAR_FORMATION
        if(units->temperature*cell_gas_temperature(icell) < feedback_temperature_ceiling)
        {
            cell_gas_internal_energy(icell) += dU;
            cell_gas_energy(icell) += dU/2.;  /* half was already kinetic */
            cell_gas_pressure(icell) += dU*(cell_gas_gamma(icell)-1);
        }
#endif /* STAR_FORMATION */
	} else {
#ifdef ISOTROPIC_TURBULENCE_ENERGY
		if(units->temperature*cell_isotropic_turbulence_temperature(icell) < feedback_turbulence_temperature_ceiling)
		{
			cell_isotropic_turbulence_energy(icell) += dU;
			cell_gas_energy(icell) += dU/2.;  /* half was already kinetic */
			cell_gas_pressure(icell) += dU*(isotropic_turbulence_gamma-1);
		}
#else
		cart_error("ISOTROPIC_TURBULENCE_ENERGY must be defined");
#endif /* ISOTROPIC_TURBULENCE_ENERGY */
	}
}

void kfb_kick_cell(int icell, int ichild, int idir[nDim], double dp, int level){
	double uni[nDim], p1, dp_cell, dp_dir;
	double dp_prev;
	int i;
	if(ichild == UNIT_VECTOR ){
		cart_rand_unit_vector(uni);
	}else if(ichild >= 0 && ichild < num_children){
		cart_rand_unit_vector_oct(uni, ichild);
	}else if(ichild == CUBE_VECTOR){ 
		cart_rand_unit_vector_block(uni, idir);
	}else if(ichild == NEIGHBOR_VECTOR){ 
		uni[0]=(double)idir[0];
		uni[1]=(double)idir[1];
		uni[2]=(double)idir[2];
	}else{
		cart_error("bad child argument %d",ichild);
	}

	/* max momentum that can be added to cell corresponding to feedback_speed_time_ceiling */
	dp_cell = copysign(1.0,dp)*MIN( fabs(dp*cell_volume_inverse[level]), 
			dvfact*cell_gas_density(icell) ); 
	/////////////////////
	/* plot 'checkoctrandir' u 6:7:3:4 w vec */
	/* plot 'checkoctrandir' u 6:7:3:4 w vec, 'checkoctrandir' u 5:6:2:3 w vec, 'checkoctran\
	   dir' u 5:7:2:4 w vec */
/* 	double pos[nDim]; */
/* 	cell_center_position(icell, pos); */
/* 	printf("kicking %f %f %f  %f %f %f  %e\n", pos[0], pos[1], pos[2],uni[0], uni[1], uni[2],dp_cell); */
	/////////////////////
	for(i=0; i<nDim; i++){
		dp_dir = dp_cell * uni[i]; 
		p1 = cell_momentum(icell,i);
		kinetic_to_internal(icell, dp_dir, p1);

		dp_prev = cell_momentum(icell,i);
		cell_momentum(icell,i) += dp_dir;
		if( isnan(cell_momentum(icell,i)) ){
			cart_debug("dp_prev %e dp_dir %e  dv_prev %e dv_dir %e dvfact %e  nden %e",
					dp_prev,
					dp_dir,
					dp_prev/cell_gas_density(icell)*units->velocity/constants->kms,
					dp_dir/cell_gas_density(icell)*units->velocity/constants->kms,
					dvfact*units->velocity/constants->kms,
					cell_gas_density(icell)*units->number_density
					);
			cart_error("cell momentum is nan after kick in feedback.kinetic.c");
		}
	}
}
void kfb_kick_oct(double dp, int level, int icell){
	int ioct, ichild;
	int dum[]={0,0,0};
	if ( level == min_level ) {
		kfb_kick_cell(icell, UNIT_VECTOR, dum, dp, level);
	} else {
		ioct = cell_parent_oct(icell);
		for(ichild=0; ichild<num_children; ichild++){
			kfb_kick_cell(oct_child(ioct,ichild), ichild, dum, dp/num_children, level ); 
		}
	}
}

/* 
// Hybrid pressurize/kick cubes: gather cell weights and then kick all cells with a single velocity
// receiver cells at level-1 receive (m/8)*v 
*/
void kfb_gather_neighbor_stencil_weights(int level, int icell, double *weight){
	/* 
	// purpose: gather weighting such that kicks produce the same velocity across each A 
	// (note kicks into a larger cell are treated as kicks into an octant of the larger cell)
	*/
	int iPar, ichild,icell_child, iside;
	int nb[num_neighbors], num_local_cells=0, j;
	cell_all_neighbors(icell,nb);
	if(      kfb_internal_spread == KFB_METHOD_HYBRID || 
	         kfb_internal_spread == KFB_METHOD_HYBRIDGATHER ){ 
		*weight = cell_gas_density(icell)*6 ; /* P*Aside = \dot{p}/Acube*Aside*6 -- treat each face like a kick */
	}else if(kfb_internal_spread == KFB_METHOD_KICKS ){ 
		*weight = 0;
	}else{
		cart_error("gather not written for selected spread method");
	}
	for(j=0; j<num_neighbors; j++){
		if( cell_is_leaf(nb[j]) ){ /* receiver same or coarser */
			/* the kicker cell is finer or the same -- take 1 octant or 1 cell worth of weight respectively */
			*weight += cell_gas_density(nb[j]);
		}else{ /* receiver finer */
			/* the kicker cell is coarser -- take weights from the side of the receiver oct and they add to 1 side's worth */
			iPar = nb[j];
			for(ichild=0; ichild<num_side_children; ichild++){
				iside = neighbor_side_child[j][ichild];
				icell_child = cell_child(iPar,iside);
				/* normalizing by volume kicked relative to icell -- these get dv=constv/8 */ 
				*weight += cell_gas_density(icell_child)/num_side_children; 
			}
		}
	}
}
void kfb_gather_cube_stencil_weights(int level, int icell, double *weight){
	/* gather weights for kicks to all 26+1 cube-neighbors */
    int iPar, ichild,icell_child;
    int nb26[CubeStencilSize], num_local_cells=0, j;
    double dpi, dv;
    GetCubeStencil(level, icell, nb26);
	if(      kfb_internal_spread == KFB_METHOD_HYBRID || 
	         kfb_internal_spread == KFB_METHOD_HYBRIDGATHER ){ 
		*weight = cell_gas_density(icell)*6 ; /* P*Aside = \dot{p}/Acube*Aside*6 -- treat each face like a kick */
	}else if(kfb_internal_spread == KFB_METHOD_KICKS ){ 
		*weight = 0;
	}else{
		cart_error("gather not written for selected spread method");
	}
	for(j=0; j<CubeStencilSize; j++){
		if(nb26[j] != -1 &&  cell_is_local(nb26[j])){
			if( cell_is_leaf(nb26[j]) ){
				/* the kicker cell is finer or the same -- take 1 octant or 1 cell worth of weight respectively */
				*weight += cell_gas_density(nb26[j]); 
			}else{
				/* kick down only one level in cube */
				iPar = nb26[j];
				for(ichild=0; ichild<num_children; ichild++){
					icell_child = cell_child(iPar,ichild);
					if(cell_is_leaf(icell_child)){ 
						*weight += cell_gas_density(nb26[j])/num_children;
					}
				}
			}
		}
	}
}

void kfb_kick_cube_constv(double dp, int level, int icell, double weight){
	/* kick 26 cube-neighbors except for nonlocal cells and difficult to reach corners */
    int iPar, ichild,icell_child;
    int nb26[CubeStencilSize], j;
    double dpi, dv;
    GetCubeStencil(level, icell, nb26);
    
    dv = dp/(weight); /* constv*/
    for(j=0; j<CubeStencilSize; j++){
	    if( nb26[j] != -1 && cell_is_local(nb26[j]) ){
		    if( cell_is_leaf(nb26[j]) ){
			    dpi = dv*cell_gas_density( nb26[j] );
			    kfb_kick_cell(nb26[j], CUBE_VECTOR, CubeDelPos[j], dpi, level);
		    }else{
				/* kick down only one level in cube */
			    iPar = nb26[j];
			    for(ichild=0; ichild<num_children; ichild++){
				    icell_child = cell_child(iPar,ichild);
				    if(cell_is_leaf(icell_child)){
					    dpi = dv*cell_gas_density( icell_child );
					    kfb_kick_cell(icell_child, CUBE_VECTOR, CubeDelPos[j], dpi, level);
				    }
				}
			}
		}
	}
}


/* 
// Purely kick different stencils
*/
void kfb_kick(double dp, int level, int icell){
    int dum[]={0,0,0};

    switch ( kfb_internal_spread ) {
    case KFB_SPREAD_CELL:
			kfb_kick_cell(icell, UNIT_VECTOR, dum, dp, level);
			break;
    case KFB_SPREAD_OCT:
			kfb_kick_oct(dp, level, icell);
			break;
    case KFB_SPREAD_CUBE:
			kfb_kick_cube(dp, level, icell);
			break;
	default:
			cart_error("bad kfb spread option %s",kfb_internal_spread); 
    }
}

/* 
// Purely pressurizing different stencils
*/
#ifdef EXTRA_PRESSURE_SOURCE
void kfb_pressurize(double dPressure, int level, int icell) {
	int ioct, ichild;
	switch ( kfb_internal_spread ) {
	case KFB_SPREAD_CELL:
		cell_extra_pressure_source(icell) += dPressure;
		break;
	case KFB_SPREAD_OCT:
		if ( level == min_level ) {
			cell_extra_pressure_source(icell) += dPressure;
		} else {
			ioct = cell_parent_oct(icell);
			for ( ichild = 0; ichild < num_children; ichild++ ) {
				/* 4.0=24faces/6faces */
				cell_extra_pressure_source(oct_child(ioct,ichild)) += dPressure/4.0;
			}
		}
		break;
	case KFB_SPREAD_CUBE:
		if( kfb_internal_method == KFB_METHOD_HYBRID || 
		    kfb_internal_method == KFB_METHOD_HYBRIDGATHER ){
			cell_extra_pressure_source(icell) += dPressure;
		}else{
			cart_error(" Pressurizing across cube is not implemented");
		}
		break;
	default:
		cart_error("bad kfb spread option %s",kfb_internal_spread); 
	}
}
#endif /* EXTRA_PRESSURE_SOURCE */

/* 
// cube stencil stuff
*/
void kfb_kick_cube(double dp, int level, int icell){
	double weight;
    kfb_gather_cube_stencil_weights(level, icell, &weight ); 
    kfb_kick_cube_constv(dp, level, icell, weight ); 
}
/* Hybrid kick/pressurize different stencils */
#ifdef EXTRA_PRESSURE_SOURCE
void kfb_hybrid_pressurize_kick(double dp, double dPressure, int level, int icell){
	/* Hybrid approach: pressurize central cell +kick 26 surrounding cells */
	double dPi, weight;
    kfb_gather_cube_stencil_weights(level, icell, &weight ); 
    dPi = dPressure * 6*cell_gas_density(icell)/weight; /* 6 directions go into pressure (6 also in weight) */
    kfb_pressurize(dPi, level, icell); 
    kfb_kick_cube_constv(dp, level, icell, weight ); 
}
void kfb_hybrid_gather(double dPressure, int level, int icell){
	/* gather hybrid */
	double dPi, weight;
    kfb_gather_neighbor_stencil_weights( level, icell, &weight ); 
    dPi = dPressure * 6*cell_gas_density(icell)/weight; /* 6 directions go into pressure (6 also in weight) */
    kfb_pressurize( dPi, level, icell); 
    /* kicks happen through cell-based feedback function kfb_kick_from_pextra */
}
#endif /* EXTRA_PRESSURE_SOURCE */


/* ---------------- The external routine! -------------------- */
void distribute_momentum(double dp, int level, int icell, double dt){
	double dPressure;
	/* dp is in momentum code units -- NOT per volume (mass*velocity * dt/time )*/
	dp *= kfb_boost_kicks;
	PLUGIN_POINT(RecordDistributedMomentum)(dp, icell, level);

	switch ( kfb_internal_method ) {
	case KFB_METHOD_KICKS:
		kfb_kick(dp, level, icell);
		break;
#ifdef EXTRA_PRESSURE_SOURCE
	case KFB_METHOD_PRESSURIZE:
		dPressure = dp / dt / (6*cell_size[level]*cell_size[level]) ; /* Pr = \dot{p}/A */
		kfb_pressurize(dPressure, level, icell);
		break;
	case KFB_METHOD_HYBRID:
		dPressure = dp / dt / (6*cell_size[level]*cell_size[level]) ; /* Pr = \dot{p}/A */
		kfb_hybrid_pressurize_kick(dp, dPressure, level, icell);
		break;
	case KFB_METHOD_HYBRIDGATHER:
		dPressure = dp / dt / (6*cell_size[level]*cell_size[level]) ; /* Pr = \dot{p}/A */
		kfb_hybrid_gather(dPressure, level, icell);
		break;
#endif
	default:
		cart_error("Invalid kfb_internal_method in distribute_momentum");
	}
}



/* 
// pull kicks from extra pressure source rather than from stars individually
*/
#ifdef EXTRA_PRESSURE_SOURCE
double constv_from_extra_pressure(int kicker, int receiver_level, double dt){
	double Ar = cell_size[receiver_level]*cell_size[receiver_level];
	double dp1side = cell_extra_pressure_source(kicker) * Ar * dt ;
	return dp1side/cell_gas_density(kicker);
}
void kfb_kick_from_pextra(int level, int icell, double dt){
	int iPar, ichild,iside, icell_child, idir, isign;
	int nb[num_neighbors], levNb[num_neighbors], num_local_cells=0, j;
	int dirNb[nDim];
	double constv, dpi;
	cell_all_neighbors(icell,nb);
	for(j=0; j<num_neighbors; j++){
		if( cell_extra_pressure_source(nb[j])/cell_gas_pressure(nb[j]) > 1e-6 ){ /* could be change in gradient, but this should be fine */ 
			levNb[j] = cell_level(nb[j]);
			
			/* -(0:-x),-(1:+x)... ; -sign because kick comes FROM neighbor direction */
			dirNb[0]=0; dirNb[1]=0; dirNb[2]=0; 
			idir = j/2; isign = -( (j%2)*2-1 ); 
			dirNb[idir] = isign; 
			
			if( levNb[j] <= level && cell_is_leaf(nb[j]) ){ /* receiver same or fine */
				/* the kicker cell is coarser or the same -- give constv kick to receiver */
				constv = constv_from_extra_pressure(nb[j], level, dt);
				dpi = constv*cell_gas_density(icell);
				kfb_kick_cell(icell, NEIGHBOR_VECTOR, dirNb, dpi, level);
				//cell_gas_momentum[icell][idir] += isign*constv*cell_gas_density(icell);
			}else if( levNb[j] == level+1 ){  /* receiver coarser */
				/* 
				// the kicker cell is finer -- only receive an the octants worth
				// of momentum, consistent with the weight gather. Even if you fill the 
				// finer oct end up with half the velocity from the face because 
				// half of the oct (the other "side") is not touched.
				*/
				iPar = nb[j];
				for(ichild=0; ichild<num_side_children; ichild++){
					iside = neighbor_side_child[j][ichild];
					icell_child = cell_child(iPar,iside);
					constv = constv_from_extra_pressure(icell_child, level, dt);
					dpi = constv*cell_gas_density(icell);
					kfb_kick_cell(icell, NEIGHBOR_VECTOR, dirNb, dpi, level);
				}
			}else{
				cart_error("kicker cell is at the wrong level %d -- should be within 1 of %d", levNb[j], level );
			}
		}
	}
}
#endif /* EXTRA_PRESSURE_SOURCE */
/* ---------------- The external routine! -------------------- */
void gather_kicks(int level, int icell, double t_next, double dt ){
    /* 
    // This collects kicks assigned through a cell field to 6-neighbors 
    // and distributes kicks to current cell.
    */
	if(kfb_internal_method == KFB_METHOD_HYBRIDGATHER){
#ifdef EXTRA_PRESSURE_SOURCE
		kfb_kick_from_pextra(level, icell, dt);
#else
		cart_error("need EXTRA_PRESSURE_SOURCE for HYBRIDGATHER kicks currently");
#endif
	}
}

#endif /*STARFORM*/
