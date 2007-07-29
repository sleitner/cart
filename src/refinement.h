#ifndef __REFINEMENT_H__
#define __REFINEMENT_H__

extern float split_tolerance;
extern float join_tolerance;
extern int num_diffusion_steps;
extern float reaction_increment;
extern float diffusion_coefficient;
extern float momentum_increment;

void modify( int level, int derefine );
void add_reaction( int cell );
void diffusion_step( int cell );
void refine( int level );
void choose_cells_to_derefine( int cell ); 
void derefine( int level ); 

#endif
