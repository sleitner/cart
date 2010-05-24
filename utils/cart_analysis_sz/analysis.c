#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "auxiliary.h"
#include "analysis.h"
#include "units.h"
#include "skiplist.h"
#include "constants.h"
#include "io.h"
#include "sfc.h"
#include "healpix.h"

#define dot(a,b) ((a[0]*b[0])+(a[1]*b[1])+(a[2]*b[2]))

double projection_extent = 1.0;
double pixel_axis[PIXEL_NPIX][3];

int compute_projected_radial_bin( double r, double rlmin, double drl ) {
	int bin;                                                                                                                                       
	double lgr = log10(r);

	if ( lgr < rlmin ) {
		bin = 0;
	} else {
		bin = (int)((lgr - rlmin)/drl) + 1;
	}

	return bin;
}

void cylindrical_coordinates( double p[nDim], double cylinder_axis[nDim], double *R, double *z ) {
	int i;
	double mag;
	double pproj[nDim];

	mag = dot(p,cylinder_axis);
	for ( i = 0; i < nDim; i++ ) {
		pproj[i] = p[i] - mag*cylinder_axis[i];
	}

	*z = fabs(mag);
	*R = sqrt(dot(pproj,pproj));
}

int point_cylinder_intersect( double p[nDim], double cylinder_axis[nDim], double cylinder_radius, 
		double cylinder_length) {

	int i;
	double pproj[nDim];
	double pperp[nDim];
	double mag;

	/* assumes cylinder_axis is a unit vector */
	mag = dot(p,cylinder_axis);
	if ( fabs(mag) < cylinder_length ) {
		for ( i = 0; i < nDim; i++ ) {
			pproj[i] = p[i] - mag*cylinder_axis[i];
		}

		if ( sqrt(dot(pproj,pproj)) < cylinder_radius ) {
			return 1;
		}
	}

	return 0;
}

int box_cylinder_intersect( double box_center[nDim], double box_size, double cylinder_axis[nDim], 
		double cylinder_radius, double cylinder_length ) {
	int i,j;
	double corner[nDim];

	/* check box center */
	if ( point_cylinder_intersect(box_center,cylinder_axis,cylinder_radius,cylinder_length) ) {
		return 1;
	}

	for ( i = 0; i < 8; i++ ) {
		for ( j = 0; j < nDim; j++ ) {
			corner[j] = box_center[j]+box_size*cell_delta[i][j];
		}

		if ( point_cylinder_intersect(corner,cylinder_axis,cylinder_radius,cylinder_length) ) {
			return 1;
		}		
	}

	return 0;
}

double rlmin, rlmax, drl;
int num_bins;
double szfact;

void cell_callback( halo_struct *halo, cell_struct *cell ) {
	int i, j, k;
	int level;
	int flag;
	int num_pixels;
	long pixels[PIXEL_NPIX];
	double box_center[nDim];
	double corner[nDim];
	double pos[nDim], r;
	double Ycell;
	int hits[PIXEL_NPIX];
#ifdef PROJECTED_PROFILES
	int profile_hits[PIXEL_NPIX][max_bins];
#endif
	double R, z;
	int bin;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	double Yecell;
#endif

	level = cell->level;
	Ycell = szfact * cell->gas_internal_energy*cell_volume[level];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	Yecell = szfact * (wmu_e/wmu) * cell->electron_internal_energy*cell_volume[level];
#endif

#ifndef PROJECTED_PROFILES
	/* if cell is completely within r_delta, just add it to all projections */
	flag = 0;
	for ( i = 0; i < 8; i++ ) {
		r = 0.0;
		for ( j = 0; j < nDim; j++ ) {
			corner[j] = cell->pos[j]+cell_size[level]*cell_delta[i][j] - halo->pos[j];
			if ( corner[j] > (double)num_grid/2. ) {
				corner[j] -= (double)num_grid;
			} else if ( corner[j] < -(double)num_grid/2. ) {
				corner[j] += (double)num_grid;
			}
			r += corner[j]*corner[j];
		}

		if ( sqrt(r) > halo->rvir ) {
			flag = 1;
			break;
		}
	}
#else
	flag = 1;
#endif

	if ( flag ) {
		/* loop over random points in the cell */
		for ( i = 0; i < nDim; i++ ) {
			box_center[i] = cell->pos[i] - halo->pos[i];
			if ( box_center[i] > (double)num_grid/2. ) {
				box_center[i] -= (double)num_grid;
			} else if ( box_center[i] < -(double)num_grid/2. ) {
				box_center[i] += (double)num_grid;
			}
		}

		num_pixels = 0;
		for ( i = 0; i < PIXEL_NPIX; i++ ) {
			hits[i] = 0;
#ifdef PROJECTED_PROFILES
			for ( j = 0; j < num_bins; j++ ) {
				profile_hits[i][j] = 0;
			}
#endif /* PROJECTED_PROFILES */

			if ( box_cylinder_intersect( box_center, cell_size[level], pixel_axis[i], 
#ifdef PROJECTED_PROFILES
					max( rbinmax/r0/1000.0, halo->rvir ),
#else
					halo->rvir + cell_size[level], 
#endif
					halo->extent + cell_size[level] ) ) {
				pixels[num_pixels++] = i;
			}
		}

		if ( num_pixels > 0 ) {
			for ( i = 0; i < num_points_per_cell; i++ ) {
				for ( j = 0; j < nDim; j++ ) {
					pos[j] = box_center[j] + cell_size[level]*(cart_rand()-0.5);
				}	

				for ( j = 0; j < num_pixels; j++ ) {
					cylindrical_coordinates( pos, pixel_axis[pixels[j]], &R, &z );

					if ( z < halo->extent ) {
						if ( R < halo->rvir ) {
							hits[pixels[j]]++;
						}

#ifdef PROJECTED_PROFILES
						bin = compute_projected_radial_bin( R, rlmin, drl );
						if ( bin < num_bins ) {
							profile_hits[pixels[j]][bin]++;
						}
#endif /* PROJECTED_PROFILES */
					}
				}
			}

			for ( j = 0; j < num_pixels; j++ ) {
				cart_assert( pixels[j] >= 0 && pixels[j] < PIXEL_NPIX );
				cart_assert( hits[pixels[j]] >= 0 && hits[pixels[j]] <= num_points_per_cell );

				halo->projected_sz_flux[pixels[j]] += Ycell*(double)hits[pixels[j]]/(double)num_points_per_cell;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
				halo->projected_electron_sz_flux[pixels[j]] += Yecell*(double)hits[pixels[j]]/(double)num_points_per_cell;
#endif
				halo->projection_volume[pixels[j]] += cell_volume[level]*(double)hits[pixels[j]]/(double)num_points_per_cell;

#ifdef PROJECTED_PROFILES
				for ( k = 0; k < num_bins; k++ ) {
					halo->projected_sz_binned[pixels[j]][k] += Ycell*(double)profile_hits[pixels[j]][k]/(double)num_points_per_cell;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
					halo->projected_electron_sz_binned[pixels[j]][k] += Yecell*(double)profile_hits[pixels[j]][k]/(double)num_points_per_cell;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
				}
#endif /* PROJECTED_PROFILES */
			}
		}
	} else {
		for ( j = 0; j < PIXEL_NPIX; j++ ) {
			halo->projected_sz_flux[j] += Ycell;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			halo->projected_electron_sz_flux[j] += Yecell;
#endif
			halo->projection_volume[j] += cell_volume[level];
		}
	}
}

void compute_halo_properties( char *output_directory, char *analysis_directory, 
		char *radius, halo_list *halos ) {
	long i;
	double axis[3];	
	int ihalo;
	int level;
	int bin;
	char filename[256];
	FILE *szlist, *eszlist, *healpix;
	FILE *szpro, *eszpro;
	FILE *volume_list;
	double total_sz_flux, total_electron_sz_flux;
	double rmid, rr;

    /* set up binning */
    rlmin = log10(rbinmin/r0/1000.0);
    rlmax = log10(rbinmax/r0/1000.0);
    drl = (rlmax - rlmin)/(float)(max_bins-1);
    num_bins = (rlmax - rlmin)/drl + 1;

	/* set up conversion constants */
	szfact = 3.4383e-15 * ( gamma - 1.0 ) * T0 * r0 * Omega0 * hubble / (aexpn*aexpn*aexpn*aexpn);
	cart_debug("szfact = %e", szfact );

	/* set up pixel axes */
	for ( i = 0; i < PIXEL_NPIX; i++ ) {
		hp_pix2vec_nest( PIXEL_NSIDE, i, pixel_axis[i] );
	}

	/* process all halos */
	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		/* clear out bins */
		for ( i = 0; i < PIXEL_NPIX; i++ ) {
			halos->list[ihalo].projected_sz_flux[i] = 0.0;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			halos->list[ihalo].projected_electron_sz_flux[i] = 0.0;
#endif
			halos->list[ihalo].projection_volume[i] = 0.0;

#ifdef PROJECTED_PROFILES
			for ( bin = 0; bin < num_bins; bin++ ) {
				halos->list[ihalo].projected_sz_binned[i][bin] = 0.0;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				halos->list[ihalo].projected_electron_sz_binned[i][bin] = 0.0;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
			}
#endif /* PROJECTED_PROFILES */
		}

#ifdef SCALED_PROJECTION_EXTENT	
		halos->list[ihalo].extent = projection_extent*halos->list[ihalo].rvir;
#else
		halos->list[ihalo].extent = projection_extent/r0;
#endif
	}

	sprintf(filename, "%s/%s_a%6.4f.grid", output_directory, jobname, aexpn );
	read_indexed_grid( filename, halos, cell_callback );

	cart_debug("done assigning gas to grid"); fflush(stdout);

	sprintf(filename, "%s/h_szlist_projected_%.2f_%s_a%6.4f.dat", 
			analysis_directory, projection_extent, radius, aexpn );
	szlist = fopen( filename, "w" );
	cart_assert( szlist != NULL );

	fprintf( szlist, "# PIXEL_NSIDE = %d\n", PIXEL_NSIDE );
	fprintf( szlist, "# PIXEL_NPIX = %d\n", PIXEL_NPIX );
#ifdef SCALED_PROJECTION_EXTENT
	fprintf( szlist, "# projected_extent = %e Rproj\n", projection_extent );
#else
	fprintf( szlist, "# projected_extent = %e [h^{-1} Mpc]\n", projection_extent );
#endif
	fprintf( szlist, "# Columns:\n" );
	fprintf( szlist, "# id Rproj Yproj values grouped by healpix projected angle\n" );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	sprintf( filename, "%s/h_electron_szlist_projected_%.2f_%s_a%6.4f.dat", 
			analysis_directory, projection_extent, radius, aexpn );
	eszlist = fopen( filename, "w" );
	cart_assert( eszlist != NULL );

	fprintf( eszlist, "# PIXEL_NSIDE = %d\n", PIXEL_NSIDE );
	fprintf( eszlist, "# PIXEL_NPIX = %d\n", PIXEL_NPIX );
#ifdef SCALED_PROJECTION_EXTENT
	fprintf( eszlist, "# projected_extent = %e Rproj\n", projection_extent );
#else
	fprintf( eszlist, "# projected_extent = %e [h^{-1} Mpc]\n", projection_extent );
#endif
	fprintf( eszlist, "# Columns:\n" );
	fprintf( eszlist, "# id Rproj Yproj values grouped by healpix projected angle\n" );
#endif

#ifdef PROJECTED_PROFILES
	sprintf(filename, "%s/h_szpro_projected_%.2f_%s_a%6.4f.dat", 
			analysis_directory, projection_extent, radius, aexpn );
	szpro = fopen( filename, "w" );
	cart_assert( szpro != NULL );

	fprintf( szpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
	fprintf( szpro, "# PIXEL_NSIDE = %d\n", PIXEL_NSIDE );
	fprintf( szpro, "# PIXEL_NPIX = %d\n", PIXEL_NPIX );
#ifdef SCALED_PROJECTION_EXTENT
	fprintf( szpro, "# projected_extent = %e Rproj\n", projection_extent );
#else
	fprintf( szpro, "# projected_extent = %e [h^{-1} Mpc]\n", projection_extent );
#endif
	fprintf( szpro, "# Profile Columns:\n");
    fprintf( szpro, "# rmid rr [/h kpc, projectd] Y(<rr) grouped by healpix projected angle\n");
    fprintf( szpro, "#########################################################################\n");

#ifdef ELECTRON_ION_NONEQUILIBRIUM
    sprintf(filename, "%s/h_electron_szpro_projected_%.2f_%s_a%6.4f.dat", 
            analysis_directory, projection_extent, radius, aexpn );
    eszpro = fopen( filename, "w" );
    cart_assert( eszpro != NULL );

    fprintf( eszpro, "# Binning: %u bins, %f to %f\n", num_bins, rbinmin, rbinmax );
    fprintf( eszpro, "# PIXEL_NSIDE = %d\n", PIXEL_NSIDE );
    fprintf( eszpro, "# PIXEL_NPIX = %d\n", PIXEL_NPIX );
#ifdef SCALED_PROJECTION_EXTENT
    fprintf( eszpro, "# projected_extent = %e Rproj\n", projection_extent );
#else
    fprintf( eszpro, "# projected_extent = %e [h^{-1} Mpc]\n", projection_extent );
#endif
    fprintf( eszpro, "# Profile Columns:\n");
    fprintf( eszpro, "# rmid rr [/h kpc, projectd] Y(<rr) grouped by healpix projected angle\n");
    fprintf( eszpro, "#########################################################################\n");
#endif

#endif /* PROJECTED_PROFILES */

	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		fprintf( szlist, "%u %e", halos->list[ihalo].id, halos->list[ihalo].rvir*r0 );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		fprintf( eszlist, "%u %e", halos->list[ihalo].id, halos->list[ihalo].rvir*r0 );
#endif

		/* process pixels */
		for ( bin = 0; bin < PIXEL_NPIX; bin++ ) {
			halos->list[ihalo].projected_sz_flux[bin] *= aexpn*aexpn*r0*r0/(hubble*hubble);
			fprintf( szlist, " %e", halos->list[ihalo].projected_sz_flux[bin] );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			halos->list[ihalo].projected_electron_sz_flux[bin] *= aexpn*aexpn*r0*r0/(hubble*hubble);
			fprintf( eszlist, " %e", halos->list[ihalo].projected_electron_sz_flux[bin] );
#endif
		}

		fprintf( szlist, "\n" );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		fprintf( eszlist, "\n" );
#endif

#ifdef PROJECTED_PROFILES
		/* make binned profiles cumulative */
		for ( i = 0; i < PIXEL_NPIX; i++ ) { 
			total_sz_flux = 0.0;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			total_electron_sz_flux = 0.0;
#endif
			for ( bin = 0; bin < num_bins; bin++ ) {
				halos->list[ihalo].projected_sz_binned[i][bin] += total_sz_flux;
				total_sz_flux = halos->list[ihalo].projected_sz_binned[i][bin];
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				halos->list[ihalo].projected_electron_sz_binned[i][bin] += total_electron_sz_flux;
				total_electron_sz_flux = halos->list[ihalo].projected_electron_sz_binned[i][bin];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
			}
		}

		fprintf( szpro, "# %u %e %e\n", halos->list[ihalo].id, 
				halos->list[ihalo].rvir, halos->list[ihalo].extent );
#ifdef ELECTRON_ION_NONEQUILIBRIUM
		fprintf( eszpro, "# %u %e %e\n", halos->list[ihalo].id, 
				halos->list[ihalo].rvir, halos->list[ihalo].extent );
#endif

		for ( bin = 0; bin < num_bins; bin++ ) {
			/*
			if ( bin == 0 ) {
				rmid = 0.5*pow( 10.0, rlmin )*1000.0*r0;
			} else {
			*/
				rmid = pow( 10.0, rlmin + ((float)bin-0.5)*drl )*1000.0*r0;
			//}

			rr = pow( 10.0, rlmin + (float)bin*drl )*1000.0*r0;

			fprintf( szpro, "%.3f %.3f", rmid, rr );
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			fprintf( eszpro, "%.3f %.3f", rmid, rr );
#endif
	
			for ( i = 0; i < PIXEL_NPIX; i++ ) {
				fprintf( szpro, " %e", halos->list[ihalo].projected_sz_binned[i][bin]*aexpn*aexpn*r0*r0/(hubble*hubble) );
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				fprintf( eszpro, " %e", halos->list[ihalo].projected_electron_sz_binned[i][bin]*aexpn*aexpn*r0*r0/(hubble*hubble) );
#endif
			}
			fprintf( szpro, "\n" );
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			fprintf( eszpro, "\n" );
#endif
		}	
#endif /* PROJECTED_PROFILES */
	}

	fclose(szlist);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	fclose(eszlist);
#endif

#ifdef PROJECTED_PROFILES
	fclose(szpro);
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	fclose(eszpro);
#endif
#endif

	sprintf( filename, "%s/h_healpix_projections.dat", analysis_directory );
	healpix = fopen( filename, "w" );
	cart_assert( healpix != NULL );

	fprintf( healpix, "# PIXEL_NSIDE = %d\n", PIXEL_NSIDE );
	fprintf( healpix, "# PIXEL_NPIX = %d\n", PIXEL_NPIX );

	for ( i = 0; i < PIXEL_NPIX; i++ ) {
		hp_pix2vec_nest( PIXEL_NSIDE, i, axis );
		fprintf( healpix, "%ld  %e %e %e\n", i, axis[0], axis[1], axis[2] );
	}
		 
	fclose(healpix);

	sprintf( filename, "%s/h_projection_volume_%s_a%6.4f.dat", analysis_directory,
			radius, aexpn  );
	volume_list = fopen( filename, "w" );
	cart_assert( volume_list != NULL );

	fprintf( volume_list, "# PIXEL_NSIDE = %d\n", PIXEL_NSIDE );
	fprintf( volume_list, "# PIXEL_NPIX = %d\n", PIXEL_NPIX );

	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		fprintf( volume_list, "%u %e", halos->list[ihalo].id, 
				M_PI*pow(halos->list[ihalo].rvir*r0,2.0)*2.0*halos->list[ihalo].extent*r0 );

		for ( bin = 0; bin < PIXEL_NPIX; bin++ ) {
			halos->list[ihalo].projection_volume[bin] *= r0*r0*r0;
			fprintf( volume_list, " %e", halos->list[ihalo].projection_volume[bin] );
		}
		fprintf( volume_list, "\n" );
	}

	fclose(volume_list);
}
