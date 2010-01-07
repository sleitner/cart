#include "config.h"


float refinement_volume_min[nDim] = { 0.0, 0.0, 0.0 };
float refinement_volume_max[nDim] = { num_grid, num_grid, num_grid };


#ifdef REFINEMENT

#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "iterators.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "timing.h"
#include "tree.h"


int cells_to_refine[num_octs];
int num_cells_to_refine;

float split_tolerance		= 0.8;
float join_tolerance		= 0.2;

int num_diffusion_steps		= 4;
float diffusion_coefficient	= 0.15;
float reaction_increment	= 0.1;
float momentum_increment	= 0.4;


void config_init_refinement()
{
  control_parameter_add3(control_parameter_float,&split_tolerance,"ref:split-tolerance","split_tolerance","wsplit","the dimensionless tolerance for adaptively splitting a given cell. Should be above <ref:join-tolerance> and below 1.");

  control_parameter_add3(control_parameter_float,&join_tolerance,"ref:join-tolerance","join_tolerance","wjoin","the dimensionless tolerance for adaptively un-splitting (joining) a given cell. Should be above 0 and below <ref:split-tolerance>.");

  control_parameter_add2(control_parameter_int,&num_diffusion_steps,"ref:diffusion-steps","num_diffusion_steps","number of diffusion steps for smoothing the refinement indicator field.");

  control_parameter_add2(control_parameter_float,&diffusion_coefficient,"ref:diffusion-coefficient","diffusion_coefficient","the diffusion coefficient for smoothing the refinement indicator field.");

  control_parameter_add2(control_parameter_float,&reaction_increment,"ref:reaction-increment","reaction_increment","--ask Doug--.");

  control_parameter_add2(control_parameter_float,&momentum_increment,"ref:momentum-increment","momentum_increment","--ask Doug--.");

  config_init_refinement_indicators();
}


void config_verify_refinement()
{
  cart_assert(split_tolerance>0.0 && split_tolerance<1.0);

  cart_assert(join_tolerance>0.0 && join_tolerance<1.0);

  cart_assert(num_diffusion_steps >= 0);

  cart_assert(diffusion_coefficient > 0.0);

  cart_assert(reaction_increment > 0.0);

  cart_assert(momentum_increment > 0.0);

  config_verify_refinement_indicators();
}


void modify( int level, int op ) {
	int i, j;
	float pos[nDim];
	int diff;
	int icell;
	int num_level_cells;
	int *level_cells;

	const int diffuse_vars[1] = { VAR_REFINEMENT_DIFFUSION };
        const int indicator_var[1] = { VAR_REFINEMENT_INDICATOR };

	cart_assert( level >= min_level && level <= max_level );

	if ( level == max_level ) {
		/* don't bother marking indicators, we won't
		 * be doing a refine or derefine call */
		return;
	}

	start_time( REFINEMENT_TIMER );

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,level)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		mark_refinement_indicators( icell, level );
	}

	end_time( WORK_TIMER );

	/* diffusion step */
	for ( diff = 0; diff < num_diffusion_steps; diff++ ) {
                start_time( DIFFUSION_STEP_TIMER );

		start_time( WORK_TIMER );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			add_reaction( icell );
		}
		end_time( WORK_TIMER );

		start_time( DIFFUSION_UPDATE_TIMER );
		update_buffer_level( level, diffuse_vars, 1 );
		end_time( DIFFUSION_UPDATE_TIMER );

		start_time( WORK_TIMER );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,level)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			diffusion_step( level, icell );
		}
		end_time( WORK_TIMER );

                end_time( DIFFUSION_STEP_TIMER );
        }

	/* check refinement mask */
#pragma omp parallel for default(none), private(i,icell,pos,j) shared(refinement_volume_min,refinement_volume_max,cell_vars,num_level_cells,level_cells)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		cell_position( icell, pos );

		for ( j = 0; j < nDim; j++ ) {
			if ( pos[j] < refinement_volume_min[j] || pos[j] > refinement_volume_max[j] ) {
				refinement_indicator( icell, 0 ) = 0.0;
			}
		}
	}

	cart_free( level_cells );

	start_time( MODIFY_UPDATE_TIMER );
        update_buffer_level( level, indicator_var, 1 );
	end_time( MODIFY_UPDATE_TIMER );

	/* do rest serially */
	refine( level );
	
	if ( op ) {
		derefine( level );
	}

	end_time( REFINEMENT_TIMER );
}

void add_reaction( int icell ) {
	/* cells which need splitting become sources, add reaction_increment to them */
	if ( refinement_indicator( icell, 0 ) >= split_tolerance ) {
		refinement_indicator( icell, 0 ) = min( 1.0, refinement_indicator( icell, 0 ) + reaction_increment );
	}

	/* we'll use refinement_indicator[1] during the diffusion step since
	 * we'll be changing refinement_indicator[0], and the diffusion step
	 * must be handled atomically */
	refinement_indicator( icell, 1 ) = refinement_indicator( icell, 0 );
}

#ifdef MOMENTUM_DIFFUSION
const float sign[num_neighbors] = {
	1.0, -1.0, 1.0, -1.0, 1.0, -1.0
};
                                                                                                                                                                                                     
const int dir[num_neighbors] = {
	0, 0, 1, 1, 2, 2
};
#endif

void diffusion_step( int level, int icell ) {
	int i;
	int neighbors[num_neighbors];

	float current_indicator = refinement_indicator( icell, 0 );
	float new_indicator = 0.0;

	cell_all_neighbors( icell, neighbors );

	for ( i = 0; i < num_neighbors; i++ ) {
		if ( cell_level( neighbors[i] ) == level ) {
			new_indicator += ( refinement_indicator( neighbors[i], 1 ) - current_indicator );

#if defined(HYDRO) && defined(MOMENTUM_DIFFUSION)
			if ( cell_is_refined( neighbors[i] ) && sign[i]*cell_momentum(neighbors[i],dir[i]) > 0.0  ) {
				new_indicator += momentum_increment*refinement_indicator( neighbors[i], 1 );
			}
#endif /* HYDRO && MOMENTUM_DIFFUSION */
		}
	}

	refinement_indicator( icell, 0 ) += diffusion_coefficient * new_indicator;
	refinement_indicator( icell, 0 ) = min( 1.0, refinement_indicator( icell, 0 ) );
}

void refine( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;
	
	start_time( WORK_TIMER );

	num_cells_to_refine = 0;
	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( refinement_indicator(icell,0) > split_tolerance &&
				cell_is_leaf(icell) ) {

			cells_to_refine[num_cells_to_refine++] = icell;

			if ( num_cells_to_refine >= num_octs ) {
				cart_error("Ran out of octs!  Increase num_octs and rerun");
			}
		}
	}
	cart_free( level_cells );

	for ( i = 0; i < num_cells_to_refine; i++ ) {
		if ( split( cells_to_refine[i] ) ) {
			/* cell was not split for some reason */
			cells_to_refine[i] = -1;
		}
	}

	end_time( WORK_TIMER );

	/* now worry about buffer cells */
	start_time( SPLIT_BUFFER_TIMER );
	split_buffer_cells( level, cells_to_refine, num_cells_to_refine );
	end_time( SPLIT_BUFFER_TIMER );
}

void choose_cells_wanting_derefinement( int cell ) {
	int i, j;
	int child;
	int neighbors[num_neighbors];

	cart_assert( cell >= 0 && cell < num_cells );
	
	if ( refinement_indicator(cell,0) < join_tolerance &&
			cell_is_refined(cell) ) {
		for ( i = 0; i < num_children; i++ ) {
			child = cell_child( cell, i );
			cart_assert( child >= 0 && child < num_cells );

			/* don't derefine if a child is split */
			if ( cell_is_refined( child ) ) {
				return;
			}

			/* don't derefine if any of the child's neighbors
			 * are refined */
			cell_all_neighbors( child, neighbors );
			for ( j = 0; j < num_neighbors; j++ ) {
				if ( cell_is_refined( neighbors[j] ) ) {
					return;
				}
			}
		}

		/* made it this far, mark for derefinement */
		refinement_indicator( cell, 0 ) = -1.0;
	}
}

void choose_cells_to_derefine( int cell ) {
	int i;
	int level;
	int neighbors[num_neighbors];

	cart_assert( cell >= 0 && cell < num_cells );

	if ( refinement_indicator(cell,0) < 0.0 ) {
		level = cell_level(cell);
		cell_all_neighbors( cell, neighbors );
	
		/* don't derefine if we're in the middle of a group
		 * of refined cells which aren't scheduled to be 
		 * derefined (check each dimension in turn) */	
		for ( i = 0; i < nDim; i++ ) {
			cart_assert( neighbors[2*i] >= 0 && neighbors[2*i] < num_cells );
			cart_assert( neighbors[2*i+1] >= 0 && neighbors[2*i+1] < num_cells );

			if ( 	cell_level( neighbors[2*i] ) == level &&
				cell_level( neighbors[2*i+1] ) == level &&
				cell_is_refined( neighbors[2*i] ) &&
				cell_is_refined( neighbors[2*i+1] ) &&

				refinement_indicator( neighbors[2*i], 0 ) > -0.5 &&
				refinement_indicator( neighbors[2*i+1],0) > -0.5 ) {

				return;
			}
		}
	
		/* add to list of cells to derefine */	
		cells_to_refine[num_cells_to_refine] = cell;
		num_cells_to_refine++;
		cart_assert( num_cells_to_refine < num_octs );
	}
}

void derefine( int level ) {
	int i;
	int *oct_list;
	int *parent_root_sfc;
	int num_octs_eliminated;
	int num_level_cells;
	int *level_cells;

	const int refine_vars[1] = { VAR_REFINEMENT_INDICATOR };

	start_time( WORK_TIMER );

	num_cells_to_refine = 0;
	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

#pragma omp parallel for default(none), private(i), shared(num_level_cells,level_cells)
	for ( i = 0; i < num_level_cells; i++ ) {
		choose_cells_wanting_derefinement( level_cells[i] );
	}
	
	/* update so we know which buffer cells are derefined */
	update_buffer_level( level, refine_vars, 1 );
	
	/* this loop needs to be serial for num_cells_to_refine++ */
	for ( i = 0; i < num_level_cells; i++ ) {
		choose_cells_to_derefine( level_cells[i] );
	}

	cart_free( level_cells );

	oct_list = cart_alloc(int, num_cells_to_refine );
	parent_root_sfc = cart_alloc(int, num_cells_to_refine );

	num_octs_eliminated = 0;
	
	for ( i = 0; i < num_cells_to_refine; i++ ) {
		oct_list[num_octs_eliminated] = cell_child_oct[ cells_to_refine[i] ];
		parent_root_sfc[num_octs_eliminated] = oct_parent_root_sfc[ oct_list[num_octs_eliminated] ];
		
		if ( join( cells_to_refine[i] ) == 0 ) {
			num_octs_eliminated++;	
		}
	}

	end_time( WORK_TIMER );

	/* now join buffer cells */
	start_time( JOIN_BUFFER_TIMER );
	join_buffer_cells( level, oct_list, parent_root_sfc, num_octs_eliminated );
	end_time( JOIN_BUFFER_TIMER );

	cart_free( oct_list );
	cart_free( parent_root_sfc );
}

#endif /* REFINEMENT */
