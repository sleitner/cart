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

#ifndef STARFORM

#ifdef ENRICH
#error "ENRICH cannot be defined without STARFORM"
#endif

#ifdef FEEDBACK
#error "FEEDBACK cannot be defined without STARFORM"
#endif

#ifdef STELLARMASSLOSS
#error "STELLARMASSLOSS cannot be defined without STARFORM"
#endif

#endif /* STARFORM */

#if defined(FEEDBACK_SNIa) && !defined(FEEDBACK)
#error "FEEDBACK_SNIa cannot be defined without FEEDBACK"
#endif

#if defined(ENRICH_SNIa) && !defined(ENRICH)
#error "ENRICH_SNIa cannot be defined without ENRICH"
#endif


/*
//  Maximum number of processors
*/
#ifndef MAX_PROCS
#define MAX_PROCS    512
#endif


/*
//  Computational domain setup
*/
#define nDim			3
#define num_grid		(1<<num_root_grid_refinements)		/* number of grid spaces in 1-d */
#define num_root_cells		(1<<(nDim*num_root_grid_refinements))   /* total number of root cells */
#define num_cells		(num_octs<<nDim)			/* number of cells in buffer */
#define num_children		(1<<nDim)		
#define num_neighbors		(2*nDim)


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
// Some other, perhaps temporary, settings
*/
#ifdef COSMOLOGY
#define LEGACY_UNITS
#endif /* COSMOLOGY */


#endif
