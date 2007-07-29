#ifndef __TREE_H__
#define __TREE_H__

#include "defs.h"
#include "parallel.h"

#define		nDim			3
#define		num_grid		(1<<num_root_grid_refinements)		/* number of grid spaces in 1-d */
#define		num_root_cells		(1<<(nDim*num_root_grid_refinements))	/* total number of root cells */
#define		num_cells		(num_octs<<nDim)			/* number of cells in buffer */
#define		num_children		(1<<nDim)		
#define		num_neighbors		(2*nDim)

#if nDim == 3
#define		num_secondary_neighbors	18
#define		num_tertiary_neighbors	8
#define		num_dependent_neighbors (num_neighbors+num_secondary_neighbors+num_tertiary_neighbors)
#else
#error  "Unsupported number of directions for num_secondary_neighbors"
#endif

#define		min_level		0
#define		max_level		(num_refinement_levels)
#define		d_x			(1.0 / (float)( 1 << max_level ))		/* width of finest grid */

#define		cell_center_offset	0.5

#define		UNREFINED_CELL		-1
#define		NULL_OCT		-1
#define		FREE_OCT_LEVEL		-1

#define		CELL_TYPE_NONLOCAL	0
#define		CELL_TYPE_LOCAL		1
#define		CELL_TYPE_BUFFER	2
#define		CELL_TYPE_ANY		-1

#ifdef GRAVITY

#ifdef HYDRO
#ifdef PARTICLES
	#define num_grav_vars	(4+nDim)
#else
	#define num_grav_vars	(3+nDim)
#endif /* PARTICLES */
#else
#ifdef PARTICLES
	#define num_grav_vars	(3+nDim)
#else
	#define num_grav_vars	(2+nDim)
#endif /* PARTICLES */
#endif /* HYDRO */

	#define VAR_DENSITY		0

#ifdef PARTICLES
	#define VAR_FIRST_SPECIES_MASS	1
	#define VAR_POTENTIAL		2

	#define cell_first_species_mass(c)	(cell_vars[c][VAR_FIRST_SPECIES_MASS])
#else
	#define VAR_POTENTIAL		1
#endif /* PARTICLES */

#ifdef HYDRO
	#define VAR_POTENTIAL_HYDRO	(VAR_POTENTIAL+1)
	#define VAR_ACCEL		(VAR_POTENTIAL+2)
#else
	#define VAR_ACCEL		(VAR_POTENTIAL+1)
#endif /* HYDRO */

	#define cell_density(c) 		(cell_vars[c][VAR_DENSITY])
	#define cell_potential(c)		(cell_vars[c][VAR_POTENTIAL])

#ifdef HYDRO
	#define cell_potential_hydro(c)		(cell_vars[c][VAR_POTENTIAL_HYDRO])
#endif /* HYDRO */

	#define cell_accel(c,d)			(cell_vars[c][VAR_ACCEL+d])
#else

#ifdef PARTICLES
	#define num_grav_vars	2
	#define VAR_DENSITY             0
        #define VAR_FIRST_SPECIES_MASS  1
	#define cell_density(c)                 (cell_vars[c][VAR_DENSITY])
        #define cell_first_species_mass(c)      (cell_vars[c][VAR_FIRST_SPECIES_MASS])
#else 
	#define num_grav_vars	0
#endif /* PARTICLES */
#endif /* GRAVITY */

#ifdef HYDRO

#ifdef ADVECT_SPECIES
	#define HVAR_ADVECTED_VARIABLES		(num_grav_vars+5+nDim)

	#ifdef METALCOOLING
		#ifdef ENRICH /* turn on enrichment by stars */
			#define HVAR_METALLICITY_II		(num_grav_vars+5+nDim)
			#define cell_gas_metallicity_II(c)	(cell_vars[c][num_grav_vars+5+nDim])

			#ifdef ENRICH_SNIa
				#define	num_chem_species	(2)

				#define HVAR_METALLICITY_Ia		(num_grav_vars+5+nDim+1)
				#define cell_gas_metallicity_Ia(c)      (cell_vars[c][num_grav_vars+5+nDim+1])

				#define cell_gas_metallicity(c)		(cell_gas_metallicity_II(c)+cell_gas_metallicity_Ia(c))
			#else
				#define num_chem_species	(1)
				#define cell_gas_metallicity(c)		cell_gas_metallicity_II(c)
			#endif /* ENRICH_SNIa */
		#else
			#error "Metal cooling defined without enrichment!"
		#endif /* ENRICH */
	#else
		/* not sure how many to have here... */
		#define num_chem_species		(0)
	#endif /* METALCOOLING */

	#define cell_advected_variable(c,v)	(cell_vars[c][num_grav_vars+5+nDim+v])
#else
	#ifdef METALCOOLING
		#error "METALCOOLING specified without ADVECT_SPECIES set!"
	#endif /* METALCOOLING */

	#define num_chem_species	(0)
#endif /* ADVECT_SPECIES */

	#define	num_hydro_vars	(5+nDim+num_chem_species)
	#define HVAR_GAS_DENSITY	num_grav_vars
	#define HVAR_GAS_ENERGY		(num_grav_vars+1)
	#define HVAR_PRESSURE		(num_grav_vars+2)
	#define HVAR_GAMMA		(num_grav_vars+3)
	#define HVAR_INTERNAL_ENERGY	(num_grav_vars+4)
	#define HVAR_MOMENTUM		(num_grav_vars+5)

	#define cell_gas_density(c)		(cell_vars[c][HVAR_GAS_DENSITY])
	#define cell_gas_energy(c)		(cell_vars[c][HVAR_GAS_ENERGY])
	#define	cell_gas_pressure(c)		(cell_vars[c][HVAR_PRESSURE])
	#define cell_gas_gamma(c)		(cell_vars[c][HVAR_GAMMA])
	#define cell_gas_internal_energy(c)	(cell_vars[c][HVAR_INTERNAL_ENERGY])
	#define cell_momentum(c,d)		(cell_vars[c][HVAR_MOMENTUM+d])

	#define cell_hydro_variable(c,v)	(cell_vars[c][num_grav_vars+v])
#else
	#define num_hydro_vars	0
#endif /* HYDRO */

#ifdef GRAVITY
	/* reuse acceleration variables for Refinement */
	#define num_refinement_vars		0
	#define VAR_REFINEMENT_INDICATOR        (VAR_ACCEL)
	#define VAR_REFINEMENT_DIFFUSION	(VAR_ACCEL+1)
	#define refinement_indicator(c,x)	(cell_vars[c][VAR_ACCEL+x])
#else
	#define num_refinement_vars		2
	#define	VAR_REFINEMENT_INDICATOR	(num_grav_vars+num_hydro_vars)
	#define VAR_REFINEMENT_DIFFUSION	(num_grav_vars+num_hydro_vars+1)
	#define refinement_indicator(c,x)	(cell_vars[c][num_grav_vars+num_hydro_vars+x])
#endif /* GRAVITY */

#define num_vars	(num_grav_vars+num_hydro_vars+num_refinement_vars)

extern int all_vars[num_vars];
extern int all_hydro_vars[num_hydro_vars];

#define cell_var(c,v)	(cell_vars[c][v])

extern float cell_vars[num_cells][num_vars];
extern int cell_child_oct[num_cells];

extern int oct_parent_cell[num_octs];
extern int oct_level[num_octs];
extern int oct_neighbors[num_octs][num_neighbors];
extern int oct_parent_root_sfc[num_octs];
extern int oct_next[num_octs];
extern int oct_prev[num_octs];
extern float oct_pos[num_octs][nDim];

extern int next_free_oct;
extern int free_oct_list;

extern int num_cells_per_level[max_level-min_level+1];
extern int local_oct_list[max_level-min_level+1];
extern int oct_list_needs_ordering[max_level-min_level+1];

extern float cell_size[max_level-min_level+1];
extern float cell_size_inverse[max_level-min_level+1];
extern float cell_volume[max_level-min_level+1];
extern float cell_volume_inverse[max_level-min_level+1];

void init_tree();
int max_level_now_global();
int max_level_now();
int max_level_local();
int max_level_buffer();
double cell_interpolate( int cell, int cell_number, int variable );
void cell_interpolation_neighbors( int cell, int cell_number, int neighbor[nDim] );
double cell_interpolate_with_neighbors( int cell, int variable, int neighbor[nDim] );
double cell_interpolate_function_with_neighbors( int cell, double function(int), int neighbor[nDim] );
double cell_interpolate_function( int cell, int cell_number, double function(int) );

void repair_neighbors();
void check_map();
int cell_num_child_leaves( int c, int level );
int cell_num_child_cells( int c, int level );

int cell_count_octs( int c, int level );
int tree_num_octs( int c, int level );
int root_cell_type( int sfc );
int cell_is_local( int cell );
int root_cell_location( int index );
void cell_free( int c );
void oct_move( int oct_old, int oct_new );
void cell_move( int cell_old, int cell_new );
int cell_parent_cell( int c );
int cell_parent_root_cell( int c );
int oct_parent_root_cell( int oct );
int root_cell_sfc_index( int icell );
int cell_parent_root_sfc( int c);
int cell_level( int c );
void cell_position( int c, float position[nDim] );
int cell_find_position( double position[nDim] );
int cell_find_position_level( int level, double position[nDim] );
int cell_find_position_above_level( int level, double position[nDim] );
int cell_contains_position( int cell, double position[nDim] );
double compute_distance_periodic( double pos1[nDim], double pos2[nDim] );
float compute_distance_periodic_float( float pos1[nDim], float pos2[nDim] );
int cell_child( int c, int j );
void oct_all_children( int oct, int child_octs[num_children] );
void cell_all_children( int c, int child_octs[num_children] );
int root_cell_neighbor( int c, int direction );
void root_cell_all_dependent_neighbors( int sfc, int neighbors[num_dependent_neighbors] );
int cell_neighbor( int c, int direction );
void cell_all_neighbors( int c, int neighbors[num_neighbors] );
int split_cell( int cell );
int join_cell( int cell );

int oct_alloc();
void oct_free( int oct );

int cell_count_cells(int c, int level);
int tree_num_octs( int c, int level );
int tree_num_cells( int c, int level );

#define cell_is_leaf(c)			( cell_child_oct[c] == UNREFINED_CELL )
#define cell_is_refined(c)		( !cell_is_leaf(c) )
#define cell_is_root_cell(c)		( c < num_cells_per_level[min_level] + num_buffer_cells[min_level] )
#define root_cell_is_local(sfc)		( sfc >= proc_sfc_index[local_proc_id] && sfc < proc_sfc_index[local_proc_id+1] )
#define cell_parent_oct(c)		( c >> nDim )
#define oct_child( oct, j )		( oct * num_children + j )
#define cell_child_number(c)		( c % num_children)

#define min(x,y)        (((x) < (y)) ? (x): (y))
#define max(x,y)        (((x) > (y)) ? (x): (y))
#define sign(x,y)       ( (y>=0) ? fabs(x) : -fabs(x) )

/* public constant arrays (precomputed tables) */
extern const int external_direction[num_children][nDim];
extern const int secondary_neighbors[num_secondary_neighbors][2];
extern const int tertiary_neighbors[num_tertiary_neighbors][2];
extern const int secondary_external_neighbors[num_children][nDim];
extern const int pyramid_vertices[num_children][nDim];
extern const int local[num_children][num_neighbors];
extern const int in_local_oct[num_children][num_neighbors];
extern const int ishift[num_neighbors][nDim];
extern const float cell_delta[num_children][nDim];
extern const int reverse_direction[num_neighbors];
extern const int neighbor_moves[num_children][nDim];

#endif
