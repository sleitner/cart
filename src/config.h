#ifndef __CONFIG_H__
#define __CONFIG_H__


/* 
// defs.h file should only be included here 
*/
#include "defs.h"
#define CONFIGURED


/*
//  The consistency of defs.h settings can be checked here...
*/
#if defined(COOLING) && !defined(HYDRO)
#error "COOLING cannot be defined without HYDRO"
#endif


#if defined(RADIATIVE_TRANSFER) && !defined(HYDRO)
#error "RADIATIVE_TRANSFER cannot be defined without HYDRO"
#endif


#if defined(ENRICH_SNIa) && !defined(ENRICH)
#error "ENRICH_SNIa cannot be defined without ENRICH"
#endif


#ifdef RADIATIVE_TRANSFER

#if defined(RT_UV) && !defined(RT_CHEMISTRY)
#error "RT_UV requires RT_CHEMISTRY"
#endif

#endif /* RADIATIVE_TRANSFER */


#ifdef ADVECT_SPECIES
#error "The ADVECT_SPECIES define is now obsolete; it is extraneous and is not needed at all."
#endif


#ifdef LAPIDUS
#error "The LAPIDUS define is now obsolete; use <apply-lapidus-viscosity> control parameter in the .cfg file."
#endif


#ifdef DENSGRADSMOOTH
#error "The DENSGRADSMOOTH define is now obsolete; use <smooth-density-gradients> control parameter in the .cfg file."
#endif


#ifdef PRESSURE_FLOOR
#error "The PRESSURE_FLOOR define is now obsolete; use <pressure-floor-min-level> control parameter in the .cfg file; set this to -1 to disable the pressure floor."
#endif


#ifdef PRESSURELESS_FLUID
#error "The PRESSURELESS_FLUID define is now obsolete; use <pressureless-fluid-eos> control parameter in the .cfg file."
#endif


#if defined(METALCOOLING) || defined(NO_METALCOOLING)
#error "Switches METALCOOLING and NO_METALCOOLING are now obsolete; metal cooling is on by default (as physically meaningful), to disable metal cooling do not use ENRICH define."
#endif


#if defined(FEEDBACK) || defined(FEEDBACK_SNIa)
#error "Switches FEEDBACK and FEEDBACK_SNIa are now obsolete; stellar feedback is on by default, set <snII:energy-per-explosion> and <snIa:energy-per-explosion> control parameters to zero in the .cfg file to disable stellar feedback of each kind."
#endif

#ifdef STELLARMASSLOSS
#error "The STELLARMASSLOSS define is now obsolete; stellar mass loss is on by default, set the <ml:loss-rate> control parameter to zero to disable the stellar mass loss."
#endif

#ifdef OLDSTYLE_PARTICLE_FILE_SINGLE_PRECISION
#error "The OLDSTYLE_PARTICLE_FILE_SINGLE_PRECISION define is now obsolete; -pfm=1: read files with double positions but single times;-pfm=2: single-precision both. "
#endif


/*
//  Check that minimum settings are actually set.
*/
#ifndef num_root_grid_refinements
#error "num_root_grid_refinements must be set."
#endif


/*
//  Maximum number of processors
*/
#ifndef MAX_PROCS
#define MAX_PROCS    512
#endif


/*
//  Default to 25 levels if the number of levels is not set.
*/
#ifndef num_refinement_levels 
#define num_refinement_levels 25
#endif


/*
//  Default to num_star_particles = num_particles if not set.
*/
#ifdef STARFORM
#if defined(num_particles) && !defined(num_star_particles)
#define num_star_particles num_particles
#endif
#endif /* STARFORM */


/*
//  Computational domain setup
*/
#define nDim		3
#define num_grid	(1<<num_root_grid_refinements)          /* number of grid spaces in 1-d */
#define num_root_cells	(1<<(nDim*num_root_grid_refinements))   /* total number of root cells */
#define num_children    (1<<nDim)		
#define num_neighbors   (2*nDim)

#ifdef num_octs

#define num_cells	(num_octs<<nDim)                        /* number of cells in buffer */
#define STATIC_MESH_DATA
#define STATIC_INIT      

#else

#define num_octs         size_oct_array
#define num_cells        size_cell_array
#ifdef STATIC_MESH_DATA
#undef STATIC_MESH_DATA
#endif

#define STATIC_INIT      = NULL

#endif


#ifdef STATIC_MESH_DATA

#define OCT_ARRAY(name)   name[num_octs]
#define CELL_ARRAY(name)  name[num_cells]
#define OCT_ARRAY2D(name,dim2)   name[num_octs][dim2]
#define CELL_VAR_ARRAY(name,nvars)  name[num_cells][nvars]

#else  /* STATIC_MESH_DATA */

#define OCT_ARRAY(name)   *name
#define CELL_ARRAY(name)  *name
#define OCT_ARRAY2D(name,dim2)   **name
#define CELL_VAR_ARRAY(name,nvars)  *name

#endif /* STATIC_MESH_DATA */


#if (nDim == 3)
#define	num_stencil		32
#define num_secondary_neighbors	12
#define num_tertiary_neighbors	8
#define num_dependent_neighbors (num_neighbors+num_secondary_neighbors+num_tertiary_neighbors)
#else
#error  "Unsupported number of dimensions."
#endif

#define min_level		0
#define max_level		(num_refinement_levels)


#define DEFINE_LEVEL_ARRAY(type,name) \
type name##_buffer[max_level-min_level+1]; \
type *const name = name##_buffer - min_level

#define DECLARE_LEVEL_ARRAY(type,name) \
extern type *const name


#ifdef RADIATIVE_TRANSFER
#include "rt_config.h"
#endif

/*
//  Some compiler-dependent switches
*/
#ifndef COMPILER_GCC
#define OPENMP_DECLARE_CONST
#endif


#endif
