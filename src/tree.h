#ifndef __TREE_H__
#define __TREE_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>


#define		cell_center_offset	0.5

#define		UNREFINED_CELL		-1
#define		NULL_OCT			-1
#define		FREE_OCT_LEVEL		-1

#define		CELL_TYPE_NONLOCAL	0
#define		CELL_TYPE_LOCAL		1
#define		CELL_TYPE_BUFFER	2
#define		CELL_TYPE_ANY		-1


#ifdef RADIATIVE_TRANSFER
  #include "rt_tree.h"
#else
  #define rt_num_vars			0
#endif /* RADIATIVE_TRANSFER */


#ifdef GRAVITY

  #ifdef HYDRO
    #ifdef PARTICLES
      #define num_grav_vars		(4+nDim)
    #else
      #define num_grav_vars		(3+nDim)
    #endif /* PARTICLES */
  #else
    #ifdef PARTICLES
      #define num_grav_vars		(3+nDim)
    #else
      #define num_grav_vars		(2+nDim)
    #endif /* PARTICLES */
  #endif /* HYDRO */

  #define VAR_DENSITY			0

  #ifdef PARTICLES
    #define VAR_FIRST_SPECIES_MASS	1
    #define VAR_POTENTIAL		2
    #define cell_first_species_mass(c)	(cell_vars[c][VAR_FIRST_SPECIES_MASS])
  #else
    #define VAR_POTENTIAL		1
  #endif /* PARTICLES */

  #ifdef HYDRO
    #define VAR_POTENTIAL_HYDRO		(VAR_POTENTIAL+1)
    #define VAR_ACCEL			(VAR_POTENTIAL+2)
  #else
    #define VAR_ACCEL			(VAR_POTENTIAL+1)
  #endif /* HYDRO */

  #define cell_density(c) 		(cell_vars[c][VAR_DENSITY])
  #define cell_potential(c)		(cell_vars[c][VAR_POTENTIAL])

  #ifdef HYDRO
    #define cell_potential_hydro(c)	(cell_vars[c][VAR_POTENTIAL_HYDRO])
  #endif /* HYDRO */

  #define cell_accel(c,d)		(cell_vars[c][VAR_ACCEL+d])

#else /* GRAVITY */

  #ifdef PARTICLES
    #define num_grav_vars		2
    #define VAR_DENSITY			0
    #define VAR_FIRST_SPECIES_MASS	1
    #define cell_density(c)		(cell_vars[c][VAR_DENSITY])
    #define cell_first_species_mass(c)	(cell_vars[c][VAR_FIRST_SPECIES_MASS])
  #else
    #define num_grav_vars		0
  #endif /* PARTICLES */

#endif /* GRAVITY */


#ifdef RADIATIVE_TRANSFER
  #define rt_grav_vars_offset     (num_grav_vars)
  #if (rt_num_vars > 0)
    #define RT_VAR_SOURCE           (rt_grav_vars_offset)
    #define cell_rt_source(c)       (cell_vars[c][RT_VAR_SOURCE])
  #endif
#endif /* RADIATIVE_TRANSFER */


#ifdef HYDRO

  #ifdef ADVECT_SPECIES

    #ifdef ELECTRON_ION_NONEQUILIBRIUM
      #define HVAR_ADVECTED_VARIABLES		(num_grav_vars+rt_num_vars+6+nDim)
    #else
      #define HVAR_ADVECTED_VARIABLES		(num_grav_vars+rt_num_vars+5+nDim)
    #endif /* ELECTRON_ION_NONEQUILIBRIUM */

    #ifdef RADIATIVE_TRANSFER /* radiative transfer block */

      #define rt_num_chem_species		6
      #define RT_HVAR_OFFSET			(HVAR_ADVECTED_VARIABLES)

      #define cell_HI_density(c)		(cell_vars[c][RT_HVAR_OFFSET+0])
      #define cell_HII_density(c)		(cell_vars[c][RT_HVAR_OFFSET+1])
      #define cell_HeI_density(c)		(cell_vars[c][RT_HVAR_OFFSET+2])
      #define cell_HeII_density(c)		(cell_vars[c][RT_HVAR_OFFSET+3])
      #define cell_HeIII_density(c)		(cell_vars[c][RT_HVAR_OFFSET+4])
      #define cell_H2_density(c)		(cell_vars[c][RT_HVAR_OFFSET+5])
      #define cell_HI_fraction(c)		(cell_vars[c][RT_HVAR_OFFSET+0]/cell_gas_density(c)/constants->XH)
      #define cell_HII_fraction(c)		(cell_vars[c][RT_HVAR_OFFSET+1]/cell_gas_density(c)/constants->XH)
      #define cell_HeI_fraction(c)		(cell_vars[c][RT_HVAR_OFFSET+2]/cell_gas_density(c)/constants->XHe)
      #define cell_HeII_fraction(c)		(cell_vars[c][RT_HVAR_OFFSET+3]/cell_gas_density(c)/constants->XHe)
      #define cell_HeIII_fraction(c)		(cell_vars[c][RT_HVAR_OFFSET+4]/cell_gas_density(c)/constants->XHe)
      #define cell_H2_fraction(c)		(cell_vars[c][RT_HVAR_OFFSET+5]/cell_gas_density(c)/(0.5*constants->XH))

    #else

      #define rt_num_chem_species		0

    #endif /* RADIATIVE_TRANSFER */

    #ifdef ENRICH /* turn on enrichment by stars */

      #define HVAR_METAL_DENSITY_II		(HVAR_ADVECTED_VARIABLES+rt_num_chem_species)
      #define cell_gas_metal_density_II(c)	(cell_vars[c][HVAR_METAL_DENSITY_II])

      #ifdef ENRICH_SNIa
        #define	num_chem_species	        (2+rt_num_chem_species)
        #define HVAR_METAL_DENSITY_Ia		(HVAR_ADVECTED_VARIABLES+rt_num_chem_species+1)
        #define cell_gas_metal_density_Ia(c)    (cell_vars[c][HVAR_METAL_DENSITY_Ia])
        #define cell_gas_metal_density(c)	(cell_gas_metal_density_II(c)+cell_gas_metal_density_Ia(c))
      #else
        #define num_chem_species	        (1+rt_num_chem_species)
        #define cell_gas_metal_density(c)	cell_gas_metal_density_II(c)
      #endif /* ENRICH_SNIa */

    #else

      #define num_chem_species		 	(rt_num_chem_species)

    #endif /* ENRICH */

    #define cell_advected_variable(c,v)		(cell_vars[c][HVAR_ADVECTED_VARIABLES+v])

  #else  /* ADVECT_SPECIES */

    #ifdef ENRICH
      #error "ENRICH specified without ADVECT_SPECIES set!"
    #endif /* ENRICH */

    #ifdef RADIATIVE_TRANSFER /* radiative transfer block */
      #error "RADIATIVE_TRANSFER specified without ADVECT_SPECIES set!"
    #endif /* RADIATIVE_TRANSFER */

    #define num_chem_species			0

  #endif /* ADVECT_SPECIES */

  #define HVAR_GAS_DENSITY			(num_grav_vars+rt_num_vars)
  #define HVAR_GAS_ENERGY			(num_grav_vars+rt_num_vars+1)
  #define HVAR_PRESSURE				(num_grav_vars+rt_num_vars+2)
  #define HVAR_GAMMA				(num_grav_vars+rt_num_vars+3)
  #define HVAR_INTERNAL_ENERGY			(num_grav_vars+rt_num_vars+4)
  #define HVAR_MOMENTUM				(num_grav_vars+rt_num_vars+5)

  #define cell_gas_density(c)			(cell_vars[c][HVAR_GAS_DENSITY])
  #define cell_gas_energy(c)			(cell_vars[c][HVAR_GAS_ENERGY])
  #define cell_gas_pressure(c)			(cell_vars[c][HVAR_PRESSURE])
  #define cell_gas_gamma(c)			(cell_vars[c][HVAR_GAMMA])
  #define cell_gas_internal_energy(c)		(cell_vars[c][HVAR_INTERNAL_ENERGY])

  float cell_gas_kinetic_energy(int cell);
  float cell_gas_temperature(int cell);

  #ifdef ELECTRON_ION_NONEQUILIBRIUM
    #define HVAR_ELECTRON_INTERNAL_ENERGY	(num_grav_vars+rt_num_vars+5+nDim)
    #define cell_electron_internal_energy(c)	(cell_vars[c][HVAR_ELECTRON_INTERNAL_ENERGY])
    #define num_hydro_vars			(6+nDim+num_chem_species)
  #else
    #define num_hydro_vars			(5+nDim+num_chem_species)
  #endif /* ELECTRON_ION_NONEQUILIBRIUM */

  #define cell_momentum(c,d)			(cell_vars[c][HVAR_MOMENTUM+d])
  
  #define cell_hydro_variable(c,v)		(cell_vars[c][num_grav_vars+rt_num_vars+v])

#else

  #define num_hydro_vars			0

#endif /* HYDRO */

#ifdef GRAVITY
/* reuse acceleration variables for Refinement */
  #define num_refinement_vars			0
  #define VAR_REFINEMENT_INDICATOR		(VAR_ACCEL)
  #define VAR_REFINEMENT_DIFFUSION		(VAR_ACCEL+1)
  #define refinement_indicator(c,x)		(cell_vars[c][VAR_ACCEL+x])
#else
  #define num_refinement_vars			2
  #define	VAR_REFINEMENT_INDICATOR	(num_grav_vars+rt_num_vars+num_hydro_vars)
  #define VAR_REFINEMENT_DIFFUSION		(num_grav_vars+rt_num_vars+num_hydro_vars+1)
  #define refinement_indicator(c,x)		(cell_vars[c][num_grav_vars+rt_num_vars+num_hydro_vars+x])
#endif /* GRAVITY */

#define num_vars				(num_grav_vars+rt_num_vars+num_hydro_vars+num_refinement_vars)

extern int all_vars[num_vars];
#ifdef HYDRO
extern int all_hydro_vars[num_hydro_vars];
#endif /* HYDRO */

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

DECLARE_LEVEL_ARRAY(int,num_cells_per_level);
DECLARE_LEVEL_ARRAY(int,local_oct_list);

DECLARE_LEVEL_ARRAY(float,cell_size);
DECLARE_LEVEL_ARRAY(float,cell_size_inverse);
DECLARE_LEVEL_ARRAY(float,cell_volume);
DECLARE_LEVEL_ARRAY(float,cell_volume_inverse);

void init_tree();
int max_level_now_global(MPI_Comm local_comm);
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
void print_cell_values(int level);
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
void cell_position_double( int c, double position[nDim] );
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
void root_cell_uniform_stencil( int sfc, int neighbors[num_stencil] );
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
#define cell_is_root_cell(c)		((c)< num_cells_per_level[min_level] + num_buffer_cells[min_level] )
#define root_cell_is_local(sfc)		( sfc >= proc_sfc_index[local_proc_id] && sfc < proc_sfc_index[local_proc_id+1] )
#define cell_parent_oct(c)		((c)>> nDim )
#define oct_child( oct, j )		((oct)*num_children + (j) )
#define cell_child_number(c)		((c) % num_children)

/* public constant arrays (precomputed tables) */
extern const int external_direction[num_children][nDim];
extern const int uniform_stencil[num_stencil][nDim];
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
