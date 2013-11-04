#ifndef __CONFIG_H__
#define __CONFIG_H__


/* 
// defs.h file should only be included here 
*/
#include "defs.h"
#define CONFIGURED


/*
//  Keep both sets of names for some clauses
*/
#if defined(STARFORM) && !defined(STAR_FORMATION)
#define STAR_FORMATION
#endif
#if !defined(STARFORM) && defined(STAR_FORMATION)
#define STARFORM
#endif

#if defined(ENRICH) && !defined(ENRICHMENT)
#define ENRICHMENT
#endif
#if !defined(ENRICH) && defined(ENRICHMENT)
#define ENRICH
#endif

#if defined(ENRICH_SNIa) && !defined(ENRICHMENT_SNIa)
#define ENRICHMENT_SNIa
#endif
#if !defined(ENRICH_SNIa) && defined(ENRICHMENT_SNIa)
#define ENRICH_SNIa
#endif


/*
//  The consistency of defs.h settings can be checked here...
*/
#if defined(COOLING) && !defined(HYDRO)
#error "COOLING cannot be defined without HYDRO"
#endif


#if defined(STAR_FORMATION) && !defined(HYDRO)
#error "STAR_FORMATION cannot be defined without HYDRO"
#endif


#if defined(RADIATIVE_TRANSFER) && !defined(HYDRO)
#error "RADIATIVE_TRANSFER cannot be defined without HYDRO"
#endif


#if defined(ENRICHMENT_SNIa) && !defined(ENRICHMENT)
#error "ENRICHMENT_SNIa cannot be defined without ENRICHMENT"
#endif


#ifdef RADIATIVE_TRANSFER

#if defined(RT_UV) && !defined(RT_CHEMISTRY)
#error "RT_UV requires RT_CHEMISTRY"
#endif

#endif /* RADIATIVE_TRANSFER */


#if defined(AGN) && !defined(STAR_PARTICLE_TYPES)
#define STAR_PARTICLE_TYPES
#endif


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
#error "Switches METALCOOLING and NO_METALCOOLING are now obsolete; metal cooling is on by default (as physically meaningful), to disable metal cooling do not use ENRICHMENT define."
#endif

#if defined(CLOUDY_COOLING)
#error "The CLOUDY_COOLING define is now obsolete."
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
//  Note: changing MAX_PROCS mid-simulation will cause problems. Just don't.
*/
#ifndef MAX_PROCS
#define MAX_PROCS    512
#endif

#if MAX_PROCS >= 1000000 /* you never know */
#define ART_PROC_FORMAT   "%07u"
#elif MAX_PROCS >= 100000
#define ART_PROC_FORMAT   "%06u"
#elif MAX_PROCS >= 10000
#define ART_PROC_FORMAT   "%05u"
#elif MAX_PROCS >= 1000
#define ART_PROC_FORMAT   "%04u"
#elif MAX_PROCS >= 100
#define ART_PROC_FORMAT   "%03u"
#else
#define ART_PROC_FORMAT   "%u"
#endif


/*
//  Default to total of 2**30 spatial resolution 
*/
#ifndef num_refinement_levels 
#define num_refinement_levels (30-num_root_grid_refinements)
#endif

/*
//  Default to num_star_particles = num_particles if not set.
*/
#ifdef STAR_FORMATION
#if defined(num_particles) 
#if !defined(num_star_particles)
#define num_star_particles num_particles
#elif num_star_particles > num_particles
#error "num_star_particles must be less than or equal to num_particles"
#endif
#endif
#else
#if defined(num_star_particles)
#error "num_star_particles is defined, but not STAR_FORMATION!"
#endif
#endif /* STAR_FORMATION */


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
#else
#error "num_octs must be set."
#endif


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


#ifdef __cplusplus

#define DEFINE_LEVEL_ARRAY(type,name) \
type name##_buffer[max_level-min_level+1]; \
type * name = name##_buffer - min_level

#define DECLARE_LEVEL_ARRAY(type,name) \
extern type * name

#else  /* __cplusplus */

#define DEFINE_LEVEL_ARRAY(type,name) \
type name##_buffer[max_level-min_level+1]; \
type *const name = name##_buffer - min_level

#define DECLARE_LEVEL_ARRAY(type,name) \
extern type *const name

#endif /* __cplusplus */


#ifdef RADIATIVE_TRANSFER
#include "rt_config.h"
#endif

/*
//  Some compiler-dependent switches
*/
#ifndef COMPILER_GCC
#define OPENMP_DECLARE_CONST
#else
#define OPENMP_NO_COLLAPSE_CLAUSE
#endif

#endif
