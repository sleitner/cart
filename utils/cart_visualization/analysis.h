#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "defs.h"
#include "auxiliary.h"
#include "halos.h"

//#define SMOOTH_PARTICLES	
/* 16 Mpc slice for 80 box
#define NUM_PIXELS	1636
#define SLICE_MIN	600
#define SLICE_MAX	1000
*/

/* 12 Mpc slice for 80 Mpc box */
/*
#define NUM_PIXELS	1227
#define SLICE_MIN	527
#define SLICE_MAX	727
//#define SLICE_MIN	577
//#define SLICE_MAX	677
*/

/* 4 Mpc slice for the 80 Mpc box */
//#define NUM_PIXELS	410
//#define SLICE_MIN	99
//#define SLICE_MAX	106
//#define SLICE_MIN	77.5
//#define SLICE_MAX	127.5
//#define SLICE_MIN	51
//#define SLICE_MAX	358

/* 8 Mpc slice for 80 Mpc box, 12 Mpc for 120 */
#define NUM_PIXELS	820
#define SLICE_MIN	308
#define SLICE_MAX	513

/* 2 Mpc slice for 120 Mpc box 
#define NUM_PIXELS	546
#define SLICE_MIN	137
#define SLICE_MAX	410
*/

/* 12 Mpc slices for 120 Mpc box */
//#define NUM_PIXELS	3276
//#define SLICE_MIN	1501
//#define SLICE_MAX	1774
//#define SLICE_MIN	3003
//#define SLICE_MAX	3549
//#define SLICE_MIN	0
//#define SLICE_MAX	3276

/* 20 Mpc projection for 120 Mpc box */
/*
#define NUM_PIXELS	5460
#define SLICE_MIN	0
#define SLICE_MAX	5460
*/

/* 4 Mpc slice for 240 Mpc box *
#define NUM_PIXELS	137
#define SLICE_MIN	50
#define SLICE_MAX	87
*/

/* CL6csf visualization (16 Mpc) 
#define NUM_PIXELS	1638
#define SLICE_MIN	768
#define SLICE_MAX	870
#define PIXEL_LEVEL	6
#define PIXEL_SIZE	(cell_size[PIXEL_LEVEL])
*/

/* 4Mpc_1Mpc */
/*
#define NUM_PIXELS      1638
#define SLICE_MIN       614
#define SLICE_MAX       1024
#define PIXEL_LEVEL     8
#define PIXEL_SIZE      (cell_size[PIXEL_LEVEL])
*/

/* 1Mpc_1Mpc */
/*
#define NUM_PIXELS	819
#define SLICE_MIN       0
#define SLICE_MAX       819
#define PIXEL_LEVEL     9
#define PIXEL_SIZE      (cell_size[PIXEL_LEVEL])
#define SMOOTH_PARTICLES
*/

/* 12x12 Mpc */
/*
#define NUM_PIXELS      1225
#define SLICE_MIN       0
#define SLICE_MAX       1225
#define PIXEL_LEVEL     6
#define PIXEL_SIZE      (cell_size[PIXEL_LEVEL])
//#define SMOOTH_PARTICLES
*/

#define PIXEL_LEVEL     6
#define PIXEL_SIZE      (cell_size[PIXEL_LEVEL])

typedef struct {
	double sigma_dm[nDim][NUM_PIXELS][NUM_PIXELS];

#ifdef STARFORM
	double sigma_star[nDim][NUM_PIXELS][NUM_PIXELS];
#endif /* STARFORM */

	double sigma_gas[nDim][NUM_PIXELS][NUM_PIXELS];
	double sigma_Tm[nDim][NUM_PIXELS][NUM_PIXELS];
	double sigma_Tew[nDim][NUM_PIXELS][NUM_PIXELS];
	double sigma_ew[nDim][NUM_PIXELS][NUM_PIXELS];
	double sigma_S[nDim][NUM_PIXELS][NUM_PIXELS];
	
	double sigma_y[nDim][NUM_PIXELS][NUM_PIXELS];

	double sigma_refine[nDim][NUM_PIXELS][NUM_PIXELS];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double sigma_Te[nDim][NUM_PIXELS][NUM_PIXELS];
	double sigma_tei[nDim][NUM_PIXELS][NUM_PIXELS];

	double sigma_ye[nDim][NUM_PIXELS][NUM_PIXELS];
#endif

#ifdef SMOOTH_PARTICLES
	int leaf_level[NUM_PIXELS][NUM_PIXELS][NUM_PIXELS];
#endif /* SMOOTH_PARTICLES */
} visualization_struct;

void compute_halo_properties( char *output_directory, char *visualization_directory, halo_list *halos );

#endif
