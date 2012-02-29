#ifndef __REFINEMENT_H__
#define __REFINEMENT_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


/*
//  These are outside the ifdef clause since they are used in io.c
//  outside the clause too.
*/
extern int spatially_limited_refinement;
extern float refinement_volume_min[nDim];
extern float refinement_volume_max[nDim];


#ifdef REFINEMENT

#define OP_REFINE           1024
#define OP_DEREFINE         2048
#define OP_FORCE_DEREFINE   4096

void config_init_refinement();
void config_verify_refinement();

void modify( int level, int op );
void add_reaction( int icell );
void diffusion_step( int level, int icell );
void refine( int level );
void choose_cells_to_derefine( int cell ); 
void derefine( int level ); 

#endif /* REFINEMENT */

#endif
