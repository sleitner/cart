#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tree.h"
#include "auxiliary.h"
#include "timestep.h"
#include "poisson.h"

#ifdef GRAVITY

void poisson( int id, type_fft_complex *fft_density, type_fft_complex *dummy ) {
	int i, j, k;
	int index, index1, index2;
	double G, G_i, G_j, G_k;
	double green[num_grid];
	double lambda;
	double trphi;

	trphi = -6.0 / (4.0 * abox[min_level] * (double)(num_grid*num_grid*num_grid) );

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

				c_re(fft_density[index]) *= G;
				c_im(fft_density[index]) *= G;
			}
		}
	}
}

#endif
