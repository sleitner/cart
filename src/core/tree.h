#ifndef __TREE_H__
#define __TREE_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>


extern int size_oct_array;
extern int size_cell_array;


#define		UNREFINED_CELL		-1
#define		NULL_OCT		-1
#define		FREE_OCT_LEVEL		-1

#define		CELL_TYPE_NONLOCAL	0
#define		CELL_TYPE_LOCAL		1
#define		CELL_TYPE_BUFFER	2

#define 	CELL_TYPE_REFINED	4
#define		CELL_TYPE_LEAF		8

#define		CELL_TYPE_ANY           ( CELL_TYPE_LOCAL | CELL_TYPE_BUFFER | CELL_TYPE_REFINED | CELL_TYPE_LEAF )
#define		CELL_TYPE_ANY_LEAF      ( CELL_TYPE_LOCAL | CELL_TYPE_BUFFER | CELL_TYPE_LEAF )
#define		CELL_TYPE_ANY_REFINED   ( CELL_TYPE_LOCAL | CELL_TYPE_BUFFER | CELL_TYPE_REFINED )

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

  #define VAR_TOTAL_MASS		0

  #ifdef PARTICLES
    #define VAR_FIRST_SPECIES_MASS	1
    #define VAR_POTENTIAL		2
    #define cell_first_species_mass(c)	(cell_var(c,VAR_FIRST_SPECIES_MASS))
  #else
    #define VAR_POTENTIAL		1
  #endif /* PARTICLES */

  #ifdef HYDRO
    #define VAR_POTENTIAL_HYDRO		(VAR_POTENTIAL+1)
    #define VAR_ACCEL			(VAR_POTENTIAL+2)
  #else
    #define VAR_ACCEL			(VAR_POTENTIAL+1)
  #endif /* HYDRO */

  #define cell_total_mass(c) 		(cell_var(c,VAR_TOTAL_MASS))
  #define cell_potential(c)		(cell_var(c,VAR_POTENTIAL))

  #ifdef HYDRO
    #define cell_potential_hydro(c)	(cell_var(c,VAR_POTENTIAL_HYDRO))
  #endif /* HYDRO */

  #define cell_accel(c,d)		(cell_var(c,VAR_ACCEL+d))

#else /* GRAVITY */

  #ifdef PARTICLES
    #define num_grav_vars		2
    #define VAR_TOTAL_MASS		0
    #define VAR_FIRST_SPECIES_MASS	1
    #define cell_total_mass(c) 		(cell_var(c,VAR_TOTAL_MASS))
    #define cell_first_species_mass(c)	(cell_var(c,VAR_FIRST_SPECIES_MASS))
  #else
    #define num_grav_vars		0
  #endif /* PARTICLES */

#endif /* GRAVITY */


#ifdef RADIATIVE_TRANSFER
  #define rt_grav_vars_offset       (num_grav_vars)
  #if (rt_num_vars > 0)
    #define RT_VAR_SOURCE           (rt_grav_vars_offset)
    #define cell_rt_source(c)       (cell_var(c,RT_VAR_SOURCE))
  #endif
#endif /* RADIATIVE_TRANSFER */



#ifdef HYDRO

  #define HVAR_GAS_DENSITY			(num_grav_vars+rt_num_vars)
  #define HVAR_GAS_ENERGY			(num_grav_vars+rt_num_vars+1)
  #define HVAR_PRESSURE				(num_grav_vars+rt_num_vars+2)
  #define HVAR_GAMMA				(num_grav_vars+rt_num_vars+3)
  #define HVAR_INTERNAL_ENERGY			(num_grav_vars+rt_num_vars+4)
  #define HVAR_MOMENTUM				(num_grav_vars+rt_num_vars+5)

  #define cell_hydro_variable(c,v)		(cell_var(c,num_grav_vars+rt_num_vars+v))
  #define cell_gas_density(c)			(cell_var(c,HVAR_GAS_DENSITY))
  #define cell_gas_energy(c)			(cell_var(c,HVAR_GAS_ENERGY))
  #define cell_gas_pressure(c)			(cell_var(c,HVAR_PRESSURE))
  #define cell_gas_gamma(c)			(cell_var(c,HVAR_GAMMA))
  #define cell_gas_internal_energy(c)		(cell_var(c,HVAR_INTERNAL_ENERGY))
  #define cell_momentum(c,d)			(cell_var(c,HVAR_MOMENTUM+d))

  #define num_basic_hydro_vars      (5+nDim)

  #ifdef ELECTRON_ION_NONEQUILIBRIUM
    #define num_electronion_noneq_vars 1
    #define HVAR_ELECTRON_INTERNAL_ENERGY	(num_grav_vars+rt_num_vars+num_basic_hydro_vars)
    #define cell_electron_internal_energy(c)	(cell_var(c,HVAR_ELECTRON_INTERNAL_ENERGY))
  #else
    #define num_electronion_noneq_vars 0
  #endif /* ELECTRON_ION_NONEQUILIBRIUM */

  #define HVAR_EXTRA_ENERGY_VARIABLES           (num_grav_vars+rt_num_vars+num_basic_hydro_vars+num_electronion_noneq_vars)
  #define cell_extra_energy_variables(c,v)      (cell_var(c,HVAR_EXTRA_ENERGY_VARIABLES+v))
  #define cell_extra_energy_pressure(c,v)       ((extra_energy_gamma(v)-1.0)*cell_extra_energy_variables(c,v))
  #ifdef ISOTROPIC_TURBULENCE_ENERGY 
    #define num_turbulence_energy_vars           1
    #define HVAR_ISOTROPIC_TURBULENCE_ENERGY     (HVAR_EXTRA_ENERGY_VARIABLES+0)
    #define cell_isotropic_turbulence_energy(c)  (cell_var(c,HVAR_ISOTROPIC_TURBULENCE_ENERGY))
    #define isotropic_turbulence_gamma           (extra_energy_gamma(0))
  #else
    #define num_turbulence_energy_vars           0
  #endif /* ISOTROPIC_TURBULENCE_ENERGY */
  #define num_extra_energy_variables            (num_turbulence_energy_vars)

  #define num_extra_hydro_vars                  (num_electronion_noneq_vars+num_extra_energy_variables)

  #define HVAR_ADVECTED_VARIABLES		(num_grav_vars+rt_num_vars+num_basic_hydro_vars+num_extra_hydro_vars)
  #define cell_advected_variable(c,v)		(cell_var(c,HVAR_ADVECTED_VARIABLES+v))

  #ifdef RADIATIVE_TRANSFER /* radiative transfer block */
    #define rt_num_chem_species		6
    #define RT_HVAR_OFFSET			(HVAR_ADVECTED_VARIABLES)
    #define RT_HVAR_HI_DENSITY			(RT_HVAR_OFFSET+0)
    #define RT_HVAR_HII_DENSITY			(RT_HVAR_OFFSET+1)
    #define RT_HVAR_HeI_DENSITY			(RT_HVAR_OFFSET+2)
    #define RT_HVAR_HeII_DENSITY		(RT_HVAR_OFFSET+3)
    #define RT_HVAR_HeIII_DENSITY		(RT_HVAR_OFFSET+4)
    #define RT_HVAR_H2_DENSITY			(RT_HVAR_OFFSET+5)
    #define cell_HI_density(c)			(cell_var(c,RT_HVAR_OFFSET+0))
    #define cell_HII_density(c)			(cell_var(c,RT_HVAR_OFFSET+1))
    #define cell_HeI_density(c)			(cell_var(c,RT_HVAR_OFFSET+2))
    #define cell_HeII_density(c)		(cell_var(c,RT_HVAR_OFFSET+3))
    #define cell_HeIII_density(c)		(cell_var(c,RT_HVAR_OFFSET+4))
    #define cell_H2_density(c)			(cell_var(c,RT_HVAR_OFFSET+5))
    #define cell_HI_fraction(c)			(cell_var(c,RT_HVAR_OFFSET+0)/cell_gas_density(c)/constants->XH)
    #define cell_HII_fraction(c)		(cell_var(c,RT_HVAR_OFFSET+1)/cell_gas_density(c)/constants->XH)
    #define cell_HeI_fraction(c)		(cell_var(c,RT_HVAR_OFFSET+2)/cell_gas_density(c)/constants->XHe)
    #define cell_HeII_fraction(c)		(cell_var(c,RT_HVAR_OFFSET+3)/cell_gas_density(c)/constants->XHe)
    #define cell_HeIII_fraction(c)		(cell_var(c,RT_HVAR_OFFSET+4)/cell_gas_density(c)/constants->XHe)
    #define cell_H2_fraction(c)			(cell_var(c,RT_HVAR_OFFSET+5)/cell_gas_density(c)/(0.5*constants->XH))
  #else
    #define rt_num_chem_species			0
  #endif /* RADIATIVE_TRANSFER */

  #ifdef ENRICHMENT /* turn on enrichment by stars */
    #define HVAR_METAL_DENSITY_II		(HVAR_ADVECTED_VARIABLES+rt_num_chem_species)
    #define cell_gas_metal_density_II(c)	(cell_var(c,HVAR_METAL_DENSITY_II))

    #ifdef ENRICHMENT_SNIa
      #ifdef DUST_EVOLUTION
        #define num_enrichment_species          3
        #define HVAR_DUST_DENSITY               (HVAR_ADVECTED_VARIABLES+rt_num_chem_species+2)
      #else
        #define num_enrichment_species          2
      #endif
      #define HVAR_METAL_DENSITY_Ia		(HVAR_ADVECTED_VARIABLES+rt_num_chem_species+1)
      #define cell_gas_metal_density_Ia(c)	(cell_var(c,HVAR_METAL_DENSITY_Ia))
      #define cell_gas_metal_density(c)		(cell_gas_metal_density_II(c)+cell_gas_metal_density_Ia(c))
    #else
      #ifdef DUST_EVOLUTION
        #define num_enrichment_species          2
        #define HVAR_DUST_DENSITY               (HVAR_ADVECTED_VARIABLES+rt_num_chem_species+1)
      #else
        #define num_enrichment_species          1
      #endif
      #define cell_gas_metal_density(c)		cell_gas_metal_density_II(c)
    #endif /* ENRICHMENT_SNIa */

    #ifdef DUST_EVOLUTION
      #define cell_dust_density(c)              (cell_var(c,HVAR_DUST_DENSITY))
    #endif

  #else
    #define num_enrichment_species		0
  #endif /* ENRICHMENT */

  #ifdef BLASTWAVE_FEEDBACK
    #define HVAR_BLASTWAVE_TIME			(HVAR_ADVECTED_VARIABLES+rt_num_chem_species+num_enrichment_species)
    #define cell_blastwave_time(c)		(cell_var(c,HVAR_BLASTWAVE_TIME))
    #define num_feedback_species	        1
  #else
    #define num_feedback_species		0
  #endif /* BLASTWAVE_FEEDBACK*/

  #ifdef INERT_GAS_TRACER
    #define HVAR_INERT_GAS_TRACER		(HVAR_ADVECTED_VARIABLES+rt_num_chem_species+num_enrichment_species+num_feedback_species)
    #define cell_inert_gas_tracer(c)		(cell_var(c,HVAR_INERT_GAS_TRACER))
    #define num_inert_gas_tracers	        1
  #else
    #define num_inert_gas_tracers		0
  #endif /* INERT_GAS_TRACER */ 

  #define num_chem_species		 	(rt_num_chem_species+num_enrichment_species+num_feedback_species+num_inert_gas_tracers) 
  #define num_hydro_vars			(num_basic_hydro_vars+num_extra_hydro_vars+num_chem_species)

#else

  #define num_hydro_vars			0

#endif /* HYDRO */

#ifdef GRAVITY
/* reuse acceleration variables for Refinement */
  #define num_refinement_vars			0
  #define VAR_REFINEMENT_INDICATOR		(VAR_ACCEL)
  #define VAR_REFINEMENT_DIFFUSION		(VAR_ACCEL+1)
  #define refinement_indicator(c,x)		(cell_var(c,VAR_ACCEL+x))
#else
  #define num_refinement_vars			2
  #define VAR_REFINEMENT_INDICATOR	        (num_grav_vars+rt_num_vars+num_hydro_vars)
  #define VAR_REFINEMENT_DIFFUSION		(num_grav_vars+rt_num_vars+num_hydro_vars+1)
  #define refinement_indicator(c,x)		(cell_var(c,num_grav_vars+rt_num_vars+num_hydro_vars+x))
#endif /* GRAVITY */

#ifdef HYDRO
  #define VAR_EXTRA_SOURCE                      (num_grav_vars+rt_num_vars+num_hydro_vars+num_refinement_vars)
  #define cell_extra_source_variables(c,v)      (cell_var(c,VAR_EXTRA_SOURCE+v))
  #ifdef EXTRA_PRESSURE_SOURCE
    #define num_extra_pressure_source_vars      1
    #define VAR_EXTRA_PRESSURE_SOURCE           (VAR_EXTRA_SOURCE+0)
    #define cell_extra_pressure_source(c)       (cell_var(c,VAR_EXTRA_PRESSURE_SOURCE))
  #else                             
    #define num_extra_pressure_source_vars      0
  #endif /* EXTRA_PRESSURE_SOURCE */
  #define num_extra_source_vars                 (num_extra_pressure_source_vars)
#else 
  #define num_extra_source_vars                 0
#endif /* HYDRO */

#define num_vars				(num_grav_vars+rt_num_vars+num_hydro_vars+num_refinement_vars+num_extra_source_vars)

extern int all_vars[num_vars];
#ifdef HYDRO
extern int all_hydro_vars[num_hydro_vars];

#define extra_energy_gamma(v)	(extra_energy_gammas[v])
#if num_extra_energy_variables > 0 
extern float extra_energy_gammas[num_extra_energy_variables];
#else
extern float extra_energy_gammas[1]; /* avoids extra flags around pragmas*/
#endif /* num_extra_energy_variables > 0 */

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
extern double oct_pos[num_octs][nDim];

extern int next_free_oct;
extern int free_oct_list;

DECLARE_LEVEL_ARRAY(int,num_cells_per_level);
DECLARE_LEVEL_ARRAY(int,local_oct_list);

DECLARE_LEVEL_ARRAY(double,cell_size);
DECLARE_LEVEL_ARRAY(double,cell_size_inverse);
DECLARE_LEVEL_ARRAY(double,cell_volume);
DECLARE_LEVEL_ARRAY(double,cell_volume_inverse);

void init_tree();
int max_level_now_global(MPI_Comm level_com);
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
void cell_center_position( int c, double position[nDim] );
int cell_find_position( double position[nDim] );
int cell_find_position_sfc( int sfc, double position[nDim] );
int cell_find_position_level( int level, double position[nDim] );
int cell_find_position_level_sfc( int sfc, int level, double position[nDim] );
int cell_find_position_above_level( int level, double position[nDim] );
int cell_find_position_above_level_sfc( int sfc, int level, double position[nDim] );
int cell_contains_position( int cell, double position[nDim] );
double compute_distance_periodic( const double *pos1, const double *pos2 );
double compute_distance_periodic_1d( const double pos1, const double pos2 );
double compute_displacement_periodic_1d( const double pos1, const double pos2 );
void compute_displacement_periodic( const double *pos1, const double *pos2, double *dx12 );
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
extern const double cell_delta[num_children][nDim];
extern const int reverse_direction[num_neighbors];
extern const int neighbor_moves[num_children][nDim];

#endif
