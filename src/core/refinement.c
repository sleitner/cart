#include "config.h"


float refinement_volume_min[nDim] = { 0.0, 0.0, 0.0 };
float refinement_volume_max[nDim] = { num_grid, num_grid, num_grid };


#ifdef REFINEMENT

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "hydro.h"
#include "iterators.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "timing.h"
#include "tree.h"
#include "units.h"


int cells_to_refine[num_octs];
int num_cells_to_refine;

int refinement_is_static        = 0;
int spatially_limited_refinement = 0;
int refinement_volume_level     = min_level;

float split_tolerance		= 0.8;
float join_tolerance		= 0.2;

int num_diffusion_steps		= 4;
float diffusion_coefficient	= 0.15;
float reaction_increment	= 0.1;

#ifdef MOMENTUM_DIFFUSION 
float momentum_increment	= 0.4;
#endif /* MOMENTUM_DIFFUSION */

#ifdef COSMOLOGY
double fixed_proper_resolution = 0.0;
#endif /* COSMOLOGY */

double G_code;

void config_init_refinement()
{
  control_parameter_add2(control_parameter_bool,&refinement_is_static,"ref:static","refinement_is_static","makes the refinement mesh static. All the refinement must be done manually in init_run() call.");

  control_parameter_add2(control_parameter_bool,&spatially_limited_refinement,"ref:spatial","spatially_limited_refinement","turns on/off the feature by which refinement_volume_min/max is used to limit refinement.");

  control_parameter_add2(control_parameter_int,&refinement_volume_level,"ref:volume-level","refinement_volume_level","sets the level below which refinement_volume_min/max restriction is tested.");

  control_parameter_add3(control_parameter_float,&split_tolerance,"ref:split-tolerance","split_tolerance","wsplit","the dimensionless tolerance for adaptively splitting a given cell. Should be above <ref:join-tolerance> and below 1.");

  control_parameter_add3(control_parameter_float,&join_tolerance,"ref:join-tolerance","join_tolerance","wjoin","the dimensionless tolerance for adaptively un-splitting (joining) a given cell. Should be above 0 and below <ref:split-tolerance>.");

  control_parameter_add2(control_parameter_int,&num_diffusion_steps,"ref:diffusion-steps","num_diffusion_steps","number of diffusion steps for smoothing the refinement indicator field.");

  control_parameter_add2(control_parameter_float,&diffusion_coefficient,"ref:diffusion-coefficient","diffusion_coefficient","the diffusion coefficient for smoothing the refinement indicator field.");

  control_parameter_add2(control_parameter_float,&reaction_increment,"ref:reaction-increment","reaction_increment","the constant source term in the refinement coefficient diffusion.");

#ifdef MOMENTUM_DIFFUSION
  control_parameter_add2(control_parameter_float,&momentum_increment,"ref:momentum-increment","momentum_increment","Adds an additional increment to the refinement diffusion in the direction of the cell velocity. This helps to pre-refine in front of shock fronts.");
#endif /* MOMENTUM_DIFFUSION */

#ifdef COSMOLOGY
  control_parameter_add2(control_parameter_double,&fixed_proper_resolution,"fixed-proper-resolution","fixed_proper_resolution","if set to a non-zero value, this parameter sets spatial resolution (in physical kpc, no h-ies) to be (approximately) fixed in proper units; refinement and de-refinement are adjusted to ensure that the cell size of finest resolved level is within a factor of 2 of this value, irrespectively of the values of refinement indicators.");
#endif /* COSMOLOGY */

  config_init_refinement_indicators();
}


void config_verify_refinement()
{
  VERIFY(ref:volume-level, refinement_volume_level>=min_level && refinement_volume_level<=max_level );

  VERIFY(ref:split-tolerance, split_tolerance>0.0 && split_tolerance<1.0 );

  VERIFY(ref:join-tolerance, join_tolerance>0.0 && join_tolerance<1.0 );

  VERIFY(ref:diffusion-steps, num_diffusion_steps >= 0 );

  VERIFY(ref:diffusion-coefficient, diffusion_coefficient > 0.0 );

  VERIFY(ref:reaction-increment, reaction_increment > 0.0 );

#ifdef MOMENTUM_DIFFUSION
  VERIFY(ref:momentum-increment, momentum_increment > 0.0 );
#endif 

#ifdef COSMOLOGY
  VERIFY(fixed-proper-resolution, !(fixed_proper_resolution < 0.0) );
#endif /* COSMOLOGY */

  config_verify_refinement_indicators();
}


void modify( int level, int op ) {
	int i, j;
	double pos[nDim];
	int diff;
	int icell;
	int num_level_cells;
	int *level_cells;
	float saved_join_tolerance;

	const int diffuse_vars[1] = { VAR_REFINEMENT_DIFFUSION };
	const int indicator_var[1] = { VAR_REFINEMENT_INDICATOR };

	cart_assert( level >= min_level && level <= max_level );

	if ( refinement_is_static || level==max_level ) {
		/* don't bother marking indicators, we won't
		 * be doing a refine or derefine call */
		return;
	}

	if ( op!=OP_REFINE && op!=OP_DEREFINE && op!=OP_FORCE_DEREFINE && op!=(OP_REFINE|OP_DEREFINE) ) {
		cart_error("Obsolete usage for modify(...); parameter <op> should be a combination of switches OP_REFINE and/or OP_DEREFINE");
	}

#ifdef COSMOLOGY
	if(op&OP_REFINE && fixed_proper_resolution>1.0e-100)
	  {
	    /*
	    //  First check that we can do that at all
	    */
	    if(units->length*cell_size[min_level]/constants->kpc < fixed_proper_resolution)
	      {
		cart_error("The size of the top grid is too small to maintain the fixed proper resolution of %le kpc",fixed_proper_resolution);
	      }
	    if(units->length*cell_size[max_level]/constants->kpc > fixed_proper_resolution)
	      {
		cart_error("The max_level is too small to maintain the fixed proper resolution of %le kpc",fixed_proper_resolution);
	      }
	    /*
	    //  If we are trying to refine below the resolution limit, 
	    //  unset the refinement flag
	    */
	    if(units->length*cell_size[level]/constants->kpc < fixed_proper_resolution)
	      {
		op = OP_FORCE_DEREFINE; /* if we, say, started with the over-refined restart... */
#ifdef HYDRO
		if(pressure_floor_min_level>-1) pressure_floor_min_level = level;
#endif /* HYDRO */
	      }
	  }
#endif /* COSMOLOGY */

	start_time( REFINEMENT_TIMER );

	start_time( WORK_TIMER );

	/*
	//  Gravitational constant in code units
	*/
	G_code = constants->G*units->mass*pow(units->time,2.0)/pow(units->length,3.0);

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

	if ( spatially_limited_refinement && level >= refinement_volume_level ) {
	  /* check refinement mask */
#pragma omp parallel for default(none), private(i,icell,pos,j) shared(refinement_volume_min,refinement_volume_max,cell_vars,num_level_cells,level_cells)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			cell_center_position( icell, pos );

			for ( j = 0; j < nDim; j++ )  if ( pos[j] < refinement_volume_min[j] || pos[j] > refinement_volume_max[j] ) {
				refinement_indicator( icell, 0 ) = 0.0;
			}
		}
	}

	cart_free( level_cells );

	start_time( MODIFY_UPDATE_TIMER );
	update_buffer_level( level, indicator_var, 1 );
	end_time( MODIFY_UPDATE_TIMER );

	/* do rest serially */
	if ( op & OP_REFINE ) refine( level );
	if ( op & OP_DEREFINE ) derefine( level );
	if ( op & OP_FORCE_DEREFINE )
	  {
	    saved_join_tolerance = join_tolerance;
	    join_tolerance = 1.1;
	    derefine( level );
	    join_tolerance = saved_join_tolerance;
	  }

	end_time( REFINEMENT_TIMER );
}

void add_reaction( int icell ) {
	/* cells which need splitting become sources, add reaction_increment to them */
	if ( refinement_indicator( icell, 0 ) >= split_tolerance ) {
		refinement_indicator( icell, 0 ) = MIN( 1.0, refinement_indicator( icell, 0 ) + reaction_increment );
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
	refinement_indicator( icell, 0 ) = MIN( 1.0, refinement_indicator( icell, 0 ) );
}

void refine( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;
	
	start_time( WORK_TIMER );

	num_cells_to_refine = 0;
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( refinement_indicator(icell,0) > split_tolerance ) {
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
	split_buffer_cells( level, cells_to_refine, num_cells_to_refine );
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

	end_time( WORK_TIMER );
	
	/* update so we know which buffer cells are derefined */
	start_time( DEREFINE_UPDATE_TIMER );
	update_buffer_level( level, refine_vars, 1 );
	end_time( DEREFINE_UPDATE_TIMER );

	start_time( WORK_TIMER );

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
	join_buffer_cells( level, oct_list, parent_root_sfc, num_octs_eliminated );

	cart_free( oct_list );
	cart_free( parent_root_sfc );
}

#endif /* REFINEMENT */
