#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "auxiliary.h"
#include "analysis.h"
#include "units.h"
#include "constants.h"
#include "cooling.h"
#include "io.h"
#include "sfc.h"

double compute_distance_periodic_particle( double pos1[nDim], particle_float pos2[nDim] ) {
	int d;
	double dx, r = 0.0;

	for ( d = 0; d < nDim; d++ ) {
		dx = fabs( pos1[d] - pos2[d] );

		if ( dx > (double)(num_grid/2) ) {
			dx -= num_grid;
		}

		r += dx*dx;
	}

	return sqrt(r);
}

double compute_distance_periodic( double pos1[nDim], double pos2[nDim] ) {
        int d;
        double dx, r = 0.0;

        for ( d = 0; d < nDim; d++ ) {
                dx = fabs( pos1[d] - pos2[d] );

                if ( dx > (double)(num_grid/2) ) {
                        dx -= num_grid;
                }

                r += dx*dx;
        }

        return sqrt(r);
}

int position_to_pixel( double center, double pos ) {
	int dx;
	
	dx = ( (int)floor(pos/PIXEL_SIZE) - (int)floor(center/PIXEL_SIZE) );
	if ( abs(dx) >= num_grid/PIXEL_SIZE/2 ) {
		if ( dx > 0 ) {
			dx -= num_grid/PIXEL_SIZE;
		} else {
			dx += num_grid/PIXEL_SIZE;
		}
	}

	return dx+NUM_PIXELS/2;
}

double rfact, vfact;

#ifdef HYDRO
double fact_nH;
double Tfact, Sfact, nfact;
double dtfact, dEfact, szfact, Pfact;
double tnewstar;

void cell_callback( halo_struct *halo, cell_struct *cell ) {
	double r;
	int level;
	double rhogi;
	double Tcell, Tcell_kev;
	double Pcell, Scell, Ycell;
	double rhogl;
	double Zdum, Zldum;
	double dEcell, tcool;
	double frac;
	int pix;
	double mass;
	double ew;
	int i, j, k;
	int i1, i2, j1, j2, k1, k2;
	int ii, jj, kk;
	double Tecell, ne, loglambda, tei, Yecell;
	visualization_struct *vis;
	
	vis = (visualization_struct *)halo->data;
	level = cell->level;

	rhogi = 1.0 / cell->gas_density;

	/* grab cell properties */
	Tcell = Tfact * cell->gas_internal_energy * rhogi;
	Scell = Sfact * cell->gas_internal_energy*pow(rhogi,gamma);
	Ycell = szfact * cell->gas_internal_energy*cell_volume[level]/(cell_size[PIXEL_LEVEL]*cell_size[PIXEL_LEVEL]);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
        Tecell = Tfact * (wmu_e/wmu) * cell->electron_internal_energy * rhogi;
        ne = nfact*cell->gas_density;
        loglambda = max( 30.0, 37.8 + log(Tecell/1e7) - 0.5*log(ne/1e-5) );
        tei = 6.3e8*pow(Tecell/1e7,1.5)/(ne/1e-5)/(loglambda/40.0);
	Yecell = szfact * (wmu_e/wmu) * cell->electron_internal_energy*cell_volume[level]/(cell_size[PIXEL_LEVEL]*cell_size[PIXEL_LEVEL]);
#endif

	mass = cell->gas_density*cell_volume[level];
	ew = cell->gas_density*sqrt(Tcell)*mass;
	//mass = cell_volume[level];

	i = position_to_pixel( halo->pos[0], cell->pos[0]/*-0.5*cell_size[level]*/ );
	j = position_to_pixel( halo->pos[1], cell->pos[1]/*-0.5*cell_size[level]*/ );
	k = position_to_pixel( halo->pos[2], cell->pos[2]/*-0.5*cell_size[level]*/ );	

	if ( cell_size[level] > PIXEL_SIZE ) {
		pix = cell_size[level]/PIXEL_SIZE;

		i1 = max( 0, min( i, NUM_PIXELS ) );
		i2 = min( NUM_PIXELS, i+pix );

		j1 = max( 0, min( j, NUM_PIXELS ) );
		j2 = min( NUM_PIXELS, j+pix );

		k1 = max( 0, min( k, NUM_PIXELS ) );
		k2 = min( NUM_PIXELS, k+pix );

#ifdef SMOOTH_PARTICLES
		for ( i = i1; i < i2; i++ ) {
			for ( j = j1; j < j2; j++ ) {
				for ( k = k1; k < k2; k++ ) {
					cart_assert( vis->leaf_level[i][j][k] == -1 );
					vis->leaf_level[i][j][k] = level;
				}
			}
		}
#endif /* SMOOTH_PARTICLES */

		if ( i1 < SLICE_MAX && i2 > SLICE_MIN ) {
			frac = (double)( min(SLICE_MAX,i2) - max(SLICE_MIN,i1) )/(double)pix/(double)pix/(double)pix;
			cart_assert( frac >= 0.0 && frac <= 1.0 );

			for ( j = j1; j < j2; j++ ) {
				for ( k = k1; k < k2; k++ ) {
					cart_assert( j >= 0 && j < NUM_PIXELS );
					cart_assert( k >= 0 && k < NUM_PIXELS );

					vis->sigma_gas[0][j][k] += mass*frac;
					vis->sigma_Tm[0][j][k] += mass*Tcell*frac;
					vis->sigma_Tew[0][j][k] += ew*Tcell*frac;
					vis->sigma_ew[0][j][k] += ew*frac;
					vis->sigma_S[0][j][k] += mass*Scell*frac;
					vis->sigma_refine[0][j][k] += mass*(float)level*frac;

					vis->sigma_y[0][j][k] += Ycell*frac;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
					vis->sigma_Te[0][j][k] += mass*Tecell*frac;
					vis->sigma_tei[0][j][k] += mass*tei*frac;
					vis->sigma_ye[0][j][k] += Yecell*frac;
#endif
				}
			}
		}

		if ( j1 < SLICE_MAX && j2 > SLICE_MIN ) {
			frac = (double)( min(SLICE_MAX,j2) - max(SLICE_MIN,j1) )/(double)pix/(double)pix/(double)pix;
			cart_assert( frac >= 0.0 && frac < 1.0 );

			for ( i = i1; i < i2; i++ ) {
				for ( k = k1; k < k2; k++ ) {
					cart_assert( i >= 0 && i < NUM_PIXELS );
					cart_assert( k >= 0 && k < NUM_PIXELS );

					vis->sigma_gas[1][i][k] += mass*frac;
					vis->sigma_Tm[1][i][k] += mass*Tcell*frac;
					vis->sigma_Tew[1][i][k] += ew*Tcell*frac;
					vis->sigma_ew[1][i][k] += ew*frac;
					vis->sigma_S[1][i][k] += mass*Scell*frac;
					vis->sigma_refine[1][i][k] += mass*(float)level*frac;

					vis->sigma_y[1][i][k] += Ycell*frac;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
					vis->sigma_Te[1][i][k] += mass*Tecell*frac;
					vis->sigma_tei[1][i][k] += mass*tei*frac;
					vis->sigma_ye[1][i][k] += Yecell*frac;
#endif
				}
			}
		} 

		if ( k1 < SLICE_MAX && k2 > SLICE_MIN ) {
			frac = (double)( min(SLICE_MAX,k2) - max(SLICE_MIN,k1) )/(double)pix/(double)pix/(double)pix;
			cart_assert( frac >= 0.0 && frac < 1.0 );

			for ( i = i1; i < i2; i++ ) {
				for ( j = j1; j < j2; j++ ) {
					cart_assert( i >= 0 && i < NUM_PIXELS );
					cart_assert( j >= 0 && j < NUM_PIXELS );

					vis->sigma_gas[2][i][j] += mass*frac;
					vis->sigma_Tm[2][i][j] += mass*Tcell*frac;
					vis->sigma_Tew[2][i][j] += ew*Tcell*frac;
					vis->sigma_ew[2][i][j] += ew*frac;
					vis->sigma_S[2][i][j] += mass*Scell*frac;

					vis->sigma_refine[2][i][j] += mass*(float)level*frac;

					vis->sigma_y[2][i][j] += Ycell*frac;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
					vis->sigma_Te[2][i][j] += mass*frac*Tecell;
                                        vis->sigma_tei[2][i][j] += mass*tei*frac;
					vis->sigma_ye[2][i][j] += Yecell*frac;
#endif

                                }
                        }
                } 
	} else {
		if ( i >= 0 && i < NUM_PIXELS && j >= 0 && j < NUM_PIXELS && k >= 0 && k < NUM_PIXELS ) {
			if ( i >= SLICE_MIN && i < SLICE_MAX ) {
				vis->sigma_gas[0][j][k] += mass;
				vis->sigma_Tm[0][j][k] += mass*Tcell;
				vis->sigma_Tew[0][j][k] += ew*Tcell;
				vis->sigma_ew[0][j][k] += ew;
				vis->sigma_S[0][j][k] += mass*Scell;
				vis->sigma_refine[0][j][k] += (float)level;
				vis->sigma_y[0][j][k] += Ycell;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				vis->sigma_Te[0][j][k] += mass*Tecell;
				vis->sigma_tei[0][j][k] += mass*tei;
				vis->sigma_ye[0][j][k] += Yecell;
#endif

			}

			if ( j >= SLICE_MIN && j < SLICE_MAX ) {
				vis->sigma_gas[1][i][k] += mass;
				vis->sigma_Tm[1][i][k] += mass*Tcell;
				vis->sigma_Tew[1][i][k] += ew*Tcell;
				vis->sigma_ew[1][i][k] += ew;
				vis->sigma_S[1][i][k] += mass*Scell;
				vis->sigma_refine[1][i][k] += (float)level;
				vis->sigma_y[1][i][k] += Ycell;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				vis->sigma_Te[1][i][k] += mass*Tecell;
				vis->sigma_tei[1][i][k] += mass*tei;
				vis->sigma_ye[1][i][k] += Yecell;
#endif
			}

			if ( k >= SLICE_MIN && k < SLICE_MAX ) {
				vis->sigma_gas[2][i][j] += mass;
				vis->sigma_Tm[2][i][j] += mass*Tcell;
				vis->sigma_Tew[2][i][j] += ew*Tcell;
				vis->sigma_ew[2][i][j] += ew;
				vis->sigma_S[2][i][j] += mass*Scell;
				vis->sigma_refine[2][i][j] += (float)level;
				vis->sigma_y[2][i][j] += Ycell;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				vis->sigma_Te[2][i][j] += mass*Tecell;
				vis->sigma_tei[2][i][j] += mass*tei;
				vis->sigma_ye[2][i][j] += Yecell;
#endif
			}

#ifdef SMOOTH_PARTICLES
			vis->leaf_level[i][j][k] = level;
#endif /* SMOOTH_PARTICLES */
		}  
	}
}
#endif /* HYDRO */

void apply_particle( halo_struct *halo, particle_struct *particle ) {
	int i, j, k;
	visualization_struct *vis;

	vis = (visualization_struct *)halo->data;

	i = position_to_pixel( halo->pos[0], particle->x[0] );
	j = position_to_pixel( halo->pos[1], particle->x[1] );
	k = position_to_pixel( halo->pos[2], particle->x[2] );

	if ( particle->is_star ) {
#ifdef STARFORM
		if ( i >= SLICE_MIN && i < SLICE_MAX ) {
			if ( j >= 0 && j < NUM_PIXELS && k >= 0 && k < NUM_PIXELS ) {
				vis->sigma_star[0][j][k] += particle->mass;
			}
		}

		if ( j >= SLICE_MIN && j < SLICE_MAX ) {
			if ( i >= 0 && i < NUM_PIXELS && k >= 0 && k < NUM_PIXELS ) {
				vis->sigma_star[1][i][k] += particle->mass;
			}
		}

		if ( k >= SLICE_MIN  && k < SLICE_MAX ) {
			if ( i >= 0 && i < NUM_PIXELS && j >= 0 && j < NUM_PIXELS ) {
				vis->sigma_star[2][i][j] += particle->mass;
			}
		}
#endif /* STARFORM */
	} else {
		if ( i >= SLICE_MIN && i < SLICE_MAX ) {
			if ( j >= 0 && j < NUM_PIXELS && k >= 0 && k < NUM_PIXELS ) {
				vis->sigma_dm[0][j][k] += particle->mass;
			}
		}

		if ( j >= SLICE_MIN && j < SLICE_MAX ) {
			if ( i >= 0 && i < NUM_PIXELS && k >= 0 && k < NUM_PIXELS ) {
				vis->sigma_dm[1][i][k] += particle->mass;
			}
		}

		if ( k >= SLICE_MIN  && k < SLICE_MAX ) {
			if ( i >= 0 && i < NUM_PIXELS && j >= 0 && j < NUM_PIXELS ) {
				vis->sigma_dm[2][i][j] += particle->mass;
			}
		}
	}
}

void particle_callback( halo_struct *halo, particle_struct *particle ) {
	int i, j, k;
	int level;
	int num_parts_1d;
	particle_struct subpart;
	double step;

        visualization_struct *vis;
        vis = (visualization_struct *)halo->data;

#ifdef SMOOTH_PARTICLES	
	i = position_to_pixel( halo->pos[0], particle->x[0] );
	j = position_to_pixel( halo->pos[1], particle->x[1] );
	k = position_to_pixel( halo->pos[2], particle->x[2] );

	i = max( 0, min( i, NUM_PIXELS-1) );
	j = max( 0, min( j, NUM_PIXELS-1) );
	k = max( 0, min( k, NUM_PIXELS-1) );

	if ( particle->specie == 0 ) {
		level = vis->leaf_level[i][j][k];
		cart_assert( level >= 0 );
	} else {
		level = -1;
	}

	num_parts_1d = (1<<max(0,PIXEL_LEVEL-level));
	step = 1.0/(double)(2*num_parts_1d);

	subpart.is_star = particle->is_star; 
	subpart.mass = particle->mass/pow( (double)num_parts_1d, 3.);

	for ( i = 0; i < num_parts_1d; i++ ) {
		subpart.x[0] = particle->x[0]+((double)(2*i+1)*step-0.5)*cell_size[level];

		for ( j = 0; j < num_parts_1d; j++ ) {
			subpart.x[1] = particle->x[1]+((double)(2*j+1)*step-0.5)*cell_size[level];

			for ( k = 0; k < num_parts_1d; k++ ) {
				subpart.x[2] = particle->x[2]+((double)(2*k+1)*step-0.5)*cell_size[level];

				apply_particle( halo, &subpart );
			}
		}
	}
#else
	apply_particle( halo, particle );
#endif 

}

void write_image( char *filename, int d, double pos[nDim], double image[NUM_PIXELS][NUM_PIXELS] ) {
	FILE *output;
	int npix;
	int i, j;
	double dx;
	double physical[nDim];
	double maxval = -1e20;
	double minval = 1e20;

	output = fopen( filename, "w" );
	if ( output == NULL ) {
		fprintf( stderr, "Unable to open file %s for writing!\n", filename );
		exit(1);
	}

	npix = NUM_PIXELS;
	fwrite( &npix, sizeof(int), 1, output );	

	/* x,y,z for lower-left corner */
	for ( i = 0; i < nDim; i++ ) {
		physical[i] = (floor( pos[i] / PIXEL_SIZE ) - (double)(NUM_PIXELS/2))*PIXEL_SIZE*r0;
		if ( i == d ) {
			physical[i] += (double)SLICE_MIN*PIXEL_SIZE*r0;
		}

		if ( physical[i] < 0.0 ) {
			physical[i] += Lbox;
		}
	}

	dx = Lbox;
	fwrite( &dx, sizeof(double), 1, output );
	fwrite( physical, sizeof(double), nDim, output );
	
	/* image width, depth */
	dx = (double)NUM_PIXELS*PIXEL_SIZE*r0;
	fwrite( &dx, sizeof(double), 1, output );
	dx = (double)(SLICE_MAX-SLICE_MIN)*PIXEL_SIZE*r0;
	fwrite( &dx, sizeof(double), 1, output );
	
	for ( i = 0; i < NUM_PIXELS; i++ ) {
		for ( j = 0; j < NUM_PIXELS; j++ ) {
			if ( image[i][j] < minval ) {
				minval = image[i][j];
			} else if ( image[i][j] > maxval ) {
				maxval = image[i][j];
			}
		}
	}

	fwrite( image, sizeof(double), NUM_PIXELS*NUM_PIXELS, output );

	printf("%s: min = %e, max = %e\n", filename, minval, maxval );

	fclose(output);
}

void compute_halo_properties( char *output_directory, char *analysis_directory, halo_list *halos ) {
	int i, j, k;
	int d;
	FILE *output;
	char filename1[256];
	char filename2[256];
	char filename3[256];
	char dimension[nDim] = { 'x', 'y', 'z' };
	visualization_struct *vis;

	/* set up conversion constants */
#ifdef HYDRO
	Tfact = T0 * ( gamma - 1.0 ) / ( aexpn*aexpn );
	Pfact = P0/(aexpn*aexpn*aexpn*aexpn*aexpn);
	Sfact = S0;
	szfact = 3.4383e-15 * ( gamma - 1.0 ) * T0 * r0 * Omega0 * hubble / (aexpn*aexpn*aexpn*aexpn);
#endif

#ifdef ELECTRON_ION_NONEQUILIBRIUM
        nfact = 1.12e-5*hubble*hubble*Omega0 / wmu_e / (aexpn*aexpn*aexpn);
#endif

	rfact = 1000.0 * r0 * aexpn / hubble; /* proper kpc -> code units */
	vfact = v0 / aexpn;

#ifdef COOLING
	dtfact = 1e6 * t0 * aexpn;
	dEfact = t0 * aexpn/Pfact;

	fact_nH = log10( 1.12e-5 * hubble * hubble * Omega0 * ( 1.0 - Y_p ) / aexpn );
#endif

	for ( i = 0; i < halos->num_halos; i++ ) {
		vis = cart_alloc(sizeof(visualization_struct) );
		for ( j = 0; j < NUM_PIXELS; j++ ) {
			for ( k = 0; k < NUM_PIXELS; k++ ) {
				for ( d = 0; d < nDim; d++ ) {
					vis->sigma_dm[d][j][k] = 0.0;
					vis->sigma_gas[d][j][k] = 0.0;
					vis->sigma_Tm[d][j][k] = 0.0;
					vis->sigma_Tew[d][j][k] = 0.0;
					vis->sigma_ew[d][j][k] = 0.0;
					vis->sigma_S[d][j][k] = 0.0;
					vis->sigma_y[d][j][k] = 0.0;
					vis->sigma_refine[d][j][k] = 0.0;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
					vis->sigma_Te[d][j][k] = 0.0;
					vis->sigma_tei[d][j][k] = 0.0;
					vis->sigma_ye[d][j][k] = 0.0;
#endif

#ifdef STARFORM
					vis->sigma_star[d][j][k] = 0.0;
#endif /* STARFORM */
				}

#ifdef SMOOTH_PARTICLES
				for ( d = 0; d < NUM_PIXELS; d++ ) {
					vis->leaf_level[j][k][d] = -1;
				}
#endif
			}
		}
		halos->list[i].data = vis;
	}

#ifdef HYDRO
        sprintf(filename1, "%s/%s_a%6.4f.grid", output_directory, jobname, aexpn );
        read_indexed_grid( filename1, halos, cell_callback );
#endif /* HYDRO */

#ifdef PARTICLES
	/* open particle file */
	sprintf( filename1, "%s/PMcrda%6.4f.DAT", output_directory, aexpn );
	sprintf( filename2, "%s/PMcrs_indexed_a%6.4f.DAT", output_directory, aexpn );
	sprintf( filename3, "%s/stars_indexed_a%6.4f.dat", output_directory, aexpn );

	read_indexed_particles( filename1, filename2, filename3, halos, particle_callback );
#endif /* PARTICLES */

	for ( i = 0; i < halos->num_halos; i++ ) {
		vis = (visualization_struct *)halos->list[i].data;

		for ( d = 0; d < nDim; d++ ) {
#ifdef PARTICLES
			sprintf( filename1, "%s/sigma_dm_%c_%6.4f_%05u.dat", analysis_directory,
				dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_dm[d] );


#ifdef STARFORM
			sprintf( filename1, "%s/sigma_star_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_star[d] );
#endif /* STARFORM */
#endif /* PARTICLES */

#ifdef HYDRO
			sprintf( filename1, "%s/sigma_gas_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_gas[d] );

			sprintf( filename1, "%s/sigma_Tm_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_Tm[d] );

			sprintf( filename1, "%s/sigma_S_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_S[d] );

                        sprintf( filename1, "%s/sigma_ymap_%c_%6.4f_%05u.dat", analysis_directory,
                                        dimension[d], aexpn, halos->list[i].id );
                        write_image( filename1, d, halos->list[i].pos, vis->sigma_y[d] );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			sprintf( filename1, "%s/sigma_Te_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_Te[d] );

			sprintf( filename1, "%s/sigma_tei_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_tei[d] );

			sprintf( filename1, "%s/sigma_yemap_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
                        write_image( filename1, d, halos->list[i].pos, vis->sigma_ye[d] );
#endif

			sprintf( filename1, "%s/sigma_refine_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_refine[d] );

			for ( j = 0; j < NUM_PIXELS; j++ ) {
				for ( k = 0; k < NUM_PIXELS; k++ ) {
					if ( vis->sigma_ew[d][j][k] > 0.0 ) {
						vis->sigma_Tew[d][j][k] /= vis->sigma_ew[d][j][k];
					}
				}
			}

			sprintf( filename1, "%s/sigma_Tew_%c_%6.4f_%05u.dat", analysis_directory,
					dimension[d], aexpn, halos->list[i].id );
			write_image( filename1, d, halos->list[i].pos, vis->sigma_Tew[d] );
#endif /* HYDRO */
		}
	}	
}
