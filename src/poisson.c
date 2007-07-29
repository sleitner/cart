#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tree.h"
#include "auxiliary.h"
#include "timestep.h"
#include "poisson.h"

#ifdef GRAVITY

#ifdef FFT_DOUBLE
#include <drfftw.h>
#include <dfftw.h>
#else
#include <srfftw.h>
#include <sfftw.h>
#endif /* FFT_DOUBLE */


void poisson( type_fft *density_grid, type_fft *potential_grid ) {
	int i, j, k;
	int index, index1, index2;
	fftwnd_plan forward, backward;
	double G, G_i, G_j, G_k;
	double green[num_grid];
	double lambda;
	fftw_complex *fft_density;
	double trphi;

	trphi = -6.0 / (4.0 * aexp[min_level] * (double)(num_grid*num_grid*num_grid) );

	fft_density = cart_alloc( num_grid*num_grid*(num_grid/2+1) * sizeof(fftw_complex) );

	forward = rfftw3d_create_plan(num_grid, num_grid, num_grid, 
			FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE );
	rfftwnd_one_real_to_complex( forward, (fftw_real *)density_grid, fft_density);
	rfftwnd_destroy_plan(forward);

	/* precompute G(k) */
	lambda = M_PI/(double)num_grid;
	for ( i = 0; i < num_grid; i++ ) {
		green[i] = sin(lambda*(double)i)*sin(lambda*(double)i);
	}

	for ( i = 0; i < num_grid; i++ ) {
		G_i = green[i];
		index1 = num_grid*i;

		for ( j = 0; j < num_grid; j++ ) {
			G_j = G_i + green[j];
			index2 = (num_grid/2+1)*(j+index1);

			for ( k = 0; k < num_grid/2 + 1; k++ ) {
				index = k+index2;

				G_k = G_j + green[k];

				if ( i == 0 && j == 0 && k == 0 ) {
					G = 0.0;
				} else {
					G = trphi / G_k;
				}

				fft_density[index].re *= G;
				fft_density[index].im *= G;
			}
		}
	}

	/* perform inverse transform to get potential */
	backward = rfftw3d_create_plan(num_grid, num_grid, num_grid, 
			FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE );
	rfftwnd_one_complex_to_real( backward, fft_density, (fftw_real *)potential_grid);
        rfftwnd_destroy_plan(backward);

	cart_debug("done with fft");

	cart_free( fft_density );
}

#endif
