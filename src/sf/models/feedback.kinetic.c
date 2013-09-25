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
// This file tells how to distribute some input momentum to kicks and/or internal energy 
*/

#ifdef STAR_FORMATION
extern double feedback_temperature_ceiling;
#endif /* STAR_FORMATION */
extern double feedback_turbulence_temperature_ceiling;
double kfb_boost_kicks = 1;

#define KFB_METHOD_PRESSURIZE 0
#define KFB_METHOD_KICKS      1
#define KFB_METHOD_HYBRID     2
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

void control_parameter_set_kfb_method(const char *value, void *ptr, int ind) {
	if ( strcmp(value,"pressurize") == 0 ) {
		kfb_internal_method = KFB_METHOD_PRESSURIZE;
	} else if ( strcmp(value, "hybrid") == 0 ) {
		kfb_internal_method = KFB_METHOD_HYBRID;
	} else if ( strcmp(value, "kicks") == 0 ) {
		kfb_internal_method = KFB_METHOD_KICKS;
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
		kfb_internal_method == KFB_METHOD_PRESSURIZE || kfb_internal_method == KFB_METHOD_KICKS || kfb_internal_method == KFB_METHOD_HYBRID
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
/* first 6 get static direction kick, next 12 get 2D randomness (one direction fixed), last 8 corners get 3D randomness like an oct (randomness of same order as non-zero entries)*/
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

void GetCubeStencil(int level, int cell, int nb[CubeStencilSize])
{
  int j, levNb[num_neighbors];
  
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
      if(levNb[CubeDir1[j]] == level)
	{
	  /*
	  // Our neighbor in the first direction is on the same level,
	  // it is sufficient to take its neighbor in the second direction.
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[CubeDir1[j]],CubeDir2[j]);
	}
      else if(levNb[CubeDir2[j]] == level)
	{
	  /*
	  // Our neighbor in the second direction is on the same level,
	  // it is sufficient to take its neighbor in the first direction.
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[CubeDir2[j]],CubeDir1[j]);
	}
      else /* tried to get to new cell from 2 directions, both blocked, assume not low level (not necessarily true) */
	{ 
	  /*
	  // Both our neighbors are on a higher level. In that case the cell is kiddy-corner to 
	  // the local cell/oct...
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[CubeDir1[j]],CubeDir2[j]);
       }
  }
  /*
  //  Third level neighbors (skipping corners for now)
  */
   for(j=0; j<num_corners; j++){ 
       nb[CubeStencilSize-num_corners+j] = -1;
   }
/*   for(j=0; j<num_corners; j++){ */
/*       if(        cell_level( nb[CubeOrigin3a[j]] ) == level ){ */
/*           nb[CubeStencilSize-num_corners+j] = cell_neighbor(nb[CubeOrigin3a[j]],CubeDir3a[j]); */
/*       } else if( cell_level( nb[CubeOrigin3b[j]] ) == level ){ */
/*           nb[CubeStencilSize-num_corners+j] = cell_neighbor(nb[CubeOrigin3b[j]],CubeDir3b[j]); */
/*       } else{ */
/*           nb[CubeStencilSize-num_corners+j] = cell_neighbor(nb[CubeOrigin3a[j]],CubeDir3a[j]); */
/* // either way of getting to corner should be the same */
/* //          cart_assert( nb[CubeStencilSize-num_corners+j] == cell_neighbor(nb[CubeOrigin3b[j]],CubeDir3b[j]) ); */
/*       } */
/*   } */
      
}



void cart_rand_unit_vector_block(double uni[nDim], int idir[nDim]){
    double phi, r;
    /* returns a random unit in direction of oct child */
    /* children start at bottom left */
    /* 4 5   6 7 */
    /* 0 1   2 3   +x-> +z^ +y>> */
    if(idir[0] == 0 && idir[1] == 0){
	uni[2] = idir[2];
    }else{
	uni[2] = idir[2]*cart_rand();
    }
    r = sqrt(1-uni[2]*uni[2]);
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

#ifdef EXTRA_PRESSURE_SOURCE
void kfb_pressurize(double dPressure, int level, int icell) {
	int ioct, ichild;

	switch ( kfb_internal_spread ) {
		case KFB_SPREAD_CELL:
		case KFB_SPREAD_CUBE:
			cell_extra_pressure_source(icell) += dPressure;
			break;
		case KFB_SPREAD_OCT:
			if ( level == min_level ) {
				cell_extra_pressure_source(icell) += dPressure;
			} else {
				ioct = cell_parent_oct(cell);
				for ( ichild = 0; ichild < num_children; ichild++ ) {
					/* 4.0=24faces/6faces */
					cell_extra_pressure_source(oct_child(ioct,ichild)) += dPressure/4.0;
				}
			}
			break;
		default:
			cart_error("bad kfb spread option %s",kfb_internal_spread); 
	}
}

void kfb_hybrid_pressurize_kick(double dp, double dPressure, int level, int icell){
/*
 //  Hybrid approach: pressurize central cell +kick 26 surrounding cells
 */
    double dPi;
	double mall;

	cart_assert( kfb_internal_spread == KFB_SPREAD_CUBE );

    kfb_kick_cube_constv( dp, level, icell, &mall );
    dPi = dPressure * cell_gas_density(icell)/mall;
    kfb_pressurize( dPi, level, icell);
}
#endif /* EXTRA_PRESSURE_SOURCE */

void kinetic_to_internal(int icell, double p0, double p1 ){
	/* cancel momentum and convert to internal energy */
	double  p2_cancel, dU;
	p2_cancel = MIN( p0*p0, p1*p1 );
	if(p0*p1 >0){
		p2_cancel = 0;
	}

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

	if(ichild == -1){
		cart_rand_unit_vector(uni);
	}else if(ichild >= 0 && ichild < num_children){
		cart_rand_unit_vector_oct(uni, ichild);
	}else if(ichild == -2){
		cart_rand_unit_vector_block(uni, idir);
	}else{
		cart_error("bad child argument %d",ichild);
	}

	/* max momentum that can be added to cell corresponding to feedback_speed_time_ceiling */
	dp_cell = copysign(1.0,dp)*MIN( fabs(dp*cell_volume_inverse[level]), 
			dvfact*cell_gas_density(icell) ); 
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
		kfb_kick_cell(icell, -1, dum, dp, level);
	} else {
		ioct = cell_parent_oct(icell);
		for(ichild=0; ichild<num_children; ichild++){
			kfb_kick_cell(oct_child(ioct,ichild), ichild, dum, dp/num_children, level ); 
		}
	}
}

void kfb_kick_cube(double dp, int level, int icell){
	/* kick all cells with same momentum */
	int iPar, ichild,icell_child;
	int nb26[CubeStencilSize], j;
	double num_local_cells=0; 

	GetCubeStencil(level, icell, nb26);

	num_local_cells=1; /* This is not written to cross processor boundaries! */
	for(j=0; j<CubeStencilSize; j++){
		if( nb26[j] != -1 && cell_is_local(nb26[j]) ){
			if( cell_is_leaf(nb26[j]) ){
				num_local_cells++;
			}else{
				/* go down only one level in neighbors (only closest neighbors and neighbors at same level)*/
				iPar = nb26[j];
				for(ichild=0; ichild<num_children; ichild++){
					icell_child = cell_child(iPar,ichild);
					if(cell_is_leaf(icell_child)){ 
						num_local_cells += 1.0/num_children;
					}
				}
			}
		}
	}
	for(j=0; j<CubeStencilSize; j++){
		if( nb26[j] != -1 && cell_is_local(nb26[j]) ){
			if( cell_is_leaf(nb26[j]) ){
				kfb_kick_cell(nb26[j], -2, CubeDelPos[j], dp/num_local_cells, level);
			}else{
				/* go down only one level in neighbors (only closest neighbors and neighbors at same level)*/
				iPar = nb26[j];
				for(ichild=0; ichild<num_children; ichild++){
					icell_child = cell_child(iPar,ichild);
					if(cell_is_leaf(icell_child)){
						kfb_kick_cell(icell_child, -2, CubeDelPos[j], dp/num_local_cells, level);
					}
				}
			}
		}
	}
}

void kfb_kick_cube_constv(double dp, int level, int icell, double *mall_level){
	/* kick with a constant velocity corresponding to dp/mall */
	int iPar, ichild,icell_child;
	int nb26[CubeStencilSize], num_local_cells=0, j;
	double dpi, dv;
	GetCubeStencil(level, icell, nb26);

	*mall_level=cell_gas_density(icell) ; 
	for(j=0; j<CubeStencilSize; j++){
		if(nb26[j] != -1 &&  cell_is_local(nb26[j])){
			if( cell_is_leaf(nb26[j]) ){
				*mall_level += cell_gas_density(nb26[j]); 
			}else{
				/* go down only one level in neighbors (only closest neighbors and neighbors at same level)*/
				iPar = nb26[j];
				for(ichild=0; ichild<num_children; ichild++){
					icell_child = cell_child(iPar,ichild);
					if(cell_is_leaf(icell_child)){ 
						*mall_level += cell_gas_density(nb26[j])/num_children;
					}
				}
			}
		}
	}

	dv = dp/(*mall_level); /* constv*/
	for(j=0; j<CubeStencilSize; j++){
		if( nb26[j] != -1 && cell_is_local(nb26[j]) ){
			if( cell_is_leaf(nb26[j]) ){
				dpi=dv*cell_gas_density( nb26[j] );
				kfb_kick_cell(nb26[j], -2, CubeDelPos[j], dpi, level);
			}else{
				/* go down only one level in neighbors (only closest neighbors and neighbors at same level)*/
				iPar = nb26[j];
				for(ichild=0; ichild<num_children; ichild++){
					icell_child = cell_child(iPar,ichild);
					if(cell_is_leaf(icell_child)){
						dpi = dv*cell_gas_density( icell_child );
						kfb_kick_cell(icell_child, -2, CubeDelPos[j], dpi, level);
					}
				}
			}
		}
	}
}

void kfb_kick(double dp, int level, int icell){
    int dum[]={0,0,0};

    switch ( kfb_internal_spread ) {
        case KFB_SPREAD_CELL:
			kfb_kick_cell(icell, -1, dum, dp, level);
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

/* ---------------- The main routine! -------------------- */
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
		case KFB_METHOD_PRESSURIZE
			dPressure = dp / dt / (6*cell_size[level]*cell_size[level]) ; /* Pr = \dot{p}/A */
			kfb_pressurize(dPressure, level, icell);
			break;
		case KFB_METHOD_HYBRID:
			dPressure = dp / dt / (6*cell_size[level]*cell_size[level]) ; /* Pr = \dot{p}/A */
			kfb_hybrid_pressurize_kick(dp, dPressure, level, icell);
			break;
#endif
		default:
			cart_error("Invalid kfb_internal_method in distribute_momentum");
	}
}




#endif /*STARFORM*/

