#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <sfftw.h>
#include <srfftw.h>

#include "auxiliary.h"
#include "units.h"
#include "io.h"

int num_mesh;
double mesh_cell_size;
float *density;

void particle_callback( particle_struct *particle ) {
	double xs,ys,zs;
	double dx0,dx1,dy0,dy1,dz0,dz1;
	double d00, d01, d10, d11;
	int ix, ix1, iy, iy1, iz, iz1;

	xs = particle->x[0]/mesh_cell_size;
	ys = particle->x[1]/mesh_cell_size;
	zs = particle->x[2]/mesh_cell_size;

	ix = (int)(xs) % num_mesh;
	ix1 = (ix+1) % num_mesh;
	iy = (int)(ys) % num_mesh;
	iy1 = (iy+1) % num_mesh;
	iz = (int)(zs) % num_mesh;
	iz1 = (iz+1) % num_mesh;

	dx1 = xs + 0.5 - floor(xs+0.5);
	dy1 = ys + 0.5 - floor(ys+0.5);
	dz1 = zs + 0.5 - floor(zs+0.5);

	dx0 = 1.0 - dx1;
	dy0 = 1.0 - dy1;
	dz0 = 1.0 - dz1;

	dx0 *= particle->mass;
	dx1 *= particle->mass;

	d00 = dx0*dy0;
	d01 = dx0*dy1;
	d10 = dx1*dy0;
	d11 = dx1*dy1;

	density[iz+num_mesh*(iy+num_mesh*ix)] += d00*dz0;
	density[iz1+num_mesh*(iy+num_mesh*ix)] += d00*dz1;
	density[iz+num_mesh*(iy1+num_mesh*ix)] += d01*dz0;
	density[iz1+num_mesh*(iy1+num_mesh*ix)] += d01*dz1;
	density[iz+num_mesh*(iy+num_mesh*ix1)] += d10*dz0;
	density[iz1+num_mesh*(iy+num_mesh*ix1)] += d10*dz1;
	density[iz+num_mesh*(iy1+num_mesh*ix1)] += d11*dz0;
	density[iz1+num_mesh*(iy1+num_mesh*ix1)] += d11*dz1;
}

int main ( int argc, char *argv[]) {
	int i, j, k, m;
	double total, total2;
	int *num_modes;
	double *avg_k;
	double *power;
	double di, dj, d;
	double Pk, Pq, dk, wk;
	int bin, index;
	int num_power_foldings;

	int num_mesh_cells;
	int endian, nbody;
	particle_header header;

	FILE *output;

	fftwnd_plan forward;
	fftw_complex *density_fft;

	if ( argc < 6 ) {
		fprintf(stderr,"Usage: cart_power Lbox num_mesh:num_foldings PMcrd.DAT PMcrs.DAT output\n");
		exit(1);
	}

	Lbox = atof(argv[1]);

	num_power_foldings = 1;
	for ( i = 0; i < strlen( argv[2] ); i++ ) {
                if ( argv[2][i] == ':' ) {
                        argv[2][i] = '\0';
                        num_power_foldings = atoi( &argv[2][i+1] );
                        break;
                }
        }
	num_mesh = atoi(argv[2]);

	num_mesh_cells = num_mesh*num_mesh*num_mesh;

	cart_debug("num_mesh = %u", num_mesh );
	cart_debug("num_mesh_cells = %u", num_mesh_cells );

	density = cart_alloc( num_mesh_cells * sizeof(float) );

	num_modes = cart_alloc( num_mesh*sizeof(int) );
	avg_k = cart_alloc( num_mesh*sizeof(double) );
	power = cart_alloc( num_mesh*sizeof(double) );

	forward = rfftw3d_create_plan(num_mesh, num_mesh, num_mesh,
			FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE );

	/* load num_grid from particle header */
	read_particle_header( argv[3], &header, &endian, &nbody );
	num_grid = header.Ngrid;

	cart_debug("num_grid = %d", num_grid );
		
	output = fopen(argv[5], "w");
	if ( output == NULL ) {
		cart_error("Unable to open file %s\n", argv[5] );
	}

	/* write out header information */
	fprintf( output, "# k [mpc/h] P(k) mesh_folding num_modes\n" );

	/* compute density field for box folded m times */
	for ( m = 0; m < num_power_foldings; m++ ) {
		cart_debug("Folding %u", m );

		mesh_cell_size = (double)num_grid / (double)(num_mesh << m);

		for ( i = 0; i < num_mesh_cells; i++ ) {
			density[i] = 0.0;
		}

		read_particles( argv[3], argv[4], NULL, particle_callback );

		/* adjust total mass so <d> =  */
		total = 0.0;
		for ( i = 0; i < num_mesh_cells; i++ ) {
			total += density[i];
		}

		/* rescale total mass */
		total = (double)num_mesh_cells/total;
		for ( i = 0; i < num_mesh_cells; i++ ) {
			density[i] = density[i]*total - 1.0;
		}

		for ( i = 0; i < num_mesh_cells; i++ ) {
			total2 += density[i]*density[i];
		}

		printf("RMS = %e\n", sqrt( total2/(double)num_mesh_cells ) );

		printf("calculating fft...\n"); fflush(stdout);

		density_fft = cart_alloc(num_mesh*num_mesh*(num_mesh/2+1)*sizeof(fftw_complex));

		/* do fft */
		rfftwnd_one_real_to_complex( forward, (fftw_real *)density, density_fft );

		printf("binning power...\n"); fflush(stdout);

		/* now average over modes */
		for ( i = 0; i < num_mesh; i++ ) {
			num_modes[i] = 0;
			avg_k[i] = 0.0;
			power[i] = 0.0;
		}

		for ( i = 0; i < num_mesh; i++ ) {
			if ( i <= num_mesh/2 ) {
				di = i*i;
			} else {
				di = (i-num_mesh)*(i-num_mesh);
			}

			for ( j = 0; j < num_mesh; j++ ) {
				if ( j <= num_mesh/2 ) {
					dj = j*j;
				} else {
					dj = (j-num_mesh)*(j-num_mesh);
				}

				/* skip 0,0,0 mode */
				if ( i != 0 || j != 0 ) {
					d = sqrt( di + dj );

					bin = (int)(d-.5);
					index = (num_mesh/2+1) * ( j + num_mesh * i );
					Pq = (density_fft[index].re*density_fft[index].re +
							density_fft[index].im * density_fft[index].im);

					power[bin] += Pq;
					avg_k[bin] += d;
					num_modes[bin]++;
				}

				/* add power from critical mode */
				d = sqrt( di + dj + (double)(num_mesh*num_mesh/4) );
				bin = (int)(d-.5);
				index = num_mesh/2 + (num_mesh/2+1) * ( j + num_mesh * i );
				Pq = (density_fft[index].re*density_fft[index].re +
						density_fft[index].im * density_fft[index].im);

				power[bin] += Pq;
				avg_k[bin] += d;
				num_modes[bin]++;

				for ( k = 1; k < num_mesh/2; k++ ) {
					dk = k*k;
					d = sqrt( di + dj + dk );
					bin = (int)(d-.5);

					index = k + (num_mesh/2+1) * ( j + num_mesh * i );
					Pq = (density_fft[index].re*density_fft[index].re +
							density_fft[index].im * density_fft[index].im);

					power[bin] += 2.0*Pq;
					avg_k[bin] += 2.0*d;
					num_modes[bin] += 2;
				}
			}
		}

		/* now write out modes */
		for ( i = 0; i < num_mesh; i++ ) {
			//if ( num_modes[i] > 0 ) {
			if ( num_modes[i] > 0 && ( m == 0 || i > num_mesh/24 ) && i < num_mesh/12 ) {
				wk = avg_k[i]/(double)num_modes[i] * (float)(1<<m)*(2.*M_PI/Lbox);
				Pk = power[i]/(double)num_modes[i] * (Lbox*Lbox*Lbox) /
					((double)num_mesh_cells*(double)num_mesh_cells);
				fprintf(output, "%e %e %u %u\n", wk, Pk, m, num_modes[i] );
			}
		}
	}

	fclose(output);
	
        return 0;
}
