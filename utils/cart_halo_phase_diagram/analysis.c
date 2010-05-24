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
#include "halos.h"

double corner_offset[8][nDim] = {
                { -0.5, -0.5, -0.5 }, {  0.5, -0.5, -0.5 },
                { -0.5,  0.5, -0.5 }, {  0.5,  0.5, -0.5 },
                { -0.5, -0.5,  0.5 }, {  0.5, -0.5,  0.5 },
                { -0.5,  0.5,  0.5 }, {  0.5,  0.5,  0.5 }
};

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

int compute_radial_bin( double r, double rlmin, double drl ) {
	int bin;
	double lgr = log10(r);

	if ( lgr < rlmin ) {
		bin = 0;
	} else {
		bin = (int)((lgr - rlmin)/drl) + 1;
	}

	return bin;
}

double fact_nH;
double Tfact, vfact, rfact, Sfact;
double szfact, Pfact;
double rlmin, rlmax, drl;
double Tmin = 1e5;
int num_bins;

void cell_callback( halo_struct *halo, cell_struct *cell ) {
	int i, j;
	double rhogi, rhogl;
	double Tcell, Tecell;
	double Pcell, Pecell;
	double Scell;
	analysis_struct *data;
	int it, ig, ip, is, ite;
	int differ, bin, bin2, min_bin;
	int level;

	double r;
	double pos[nDim];
	double mass_fraction, volume_fraction;
	int bin_volume_fraction[max_bins];
        double min_bin_volume;

	level = cell->level;

	/* grab cell properties */
	rhogi = 1.0 / cell->gas_density;
	Tcell = Tfact * cell->gas_internal_energy*rhogi;
	Scell = Sfact * cell->gas_internal_energy*pow(rhogi,gamma);

	if ( Tcell < Tmin ) {
		Tmin = Tcell;
	}

	rhogl = log10(cell->gas_density) + fact_nH;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
        Tecell = Tfact * (wmu_e/wmu) * cell->electron_internal_energy * rhogi;
#endif

	it = (log10(Tcell) - phase_Tlmin)/phase_dlt;
	ig = (rhogl - phase_rhoglmin)/phase_drhol;
	is = (log10(Scell) - phase_Slmin)/phase_dls;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	ite = (log10(Tecell) - phase_Tlmin)/phase_dlt;
#endif

	data = (analysis_struct *)halo->data;

        r = compute_distance_periodic( halo->pos, cell->pos );

	if ( r < halo->rvir ) {
		if ( it >= 0 && ig >= 0 && it < phase_Tbins && ig < phase_rhobins ) {
			data->phase_table[it][ig] += cell_volume[level];
		}

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		if ( ite >= 0 && ig >= 0 && ite < phase_Tbins && ig < phase_rhobins ) {
			data->electron_phase_table[ite][ig] += cell_volume[level];
		}
#endif

		if ( is >= 0 && ig >= 0 && is < phase_Sbins && ig < phase_rhobins ) {
			data->entropy_table[is][ig] += cell_volume[level];
		}
	}

	r /= halo->rvir;
	bin = compute_radial_bin( r, rlmin, drl );
	min_bin = bin;

	differ = 0;
	for ( i = 0; i < 8; i++ ) {
		for ( j = 0; j < nDim; j++ ) {
			pos[j] = cell->pos[j] + corner_offset[i][j]*cell_size[level];
		}

		r = compute_distance_periodic( halo->pos, pos ) / halo->rvir;
		bin2 = compute_radial_bin( r, rlmin, drl );

                if ( bin2 < num_bins && bin != bin2 ) {
                        differ = 1;
                }

                min_bin = min( min_bin, bin2 );
        }

	if ( !differ && min_bin >= num_bins ) {
		return;
	}

        if ( !differ ) {
		cart_assert( bin >= 0 && bin < num_bins );

		if ( it >= 0 && it < phase_Tbins ) {
			data->radial_temperature[it][bin] += cell_volume[level];
		}
		if ( is >= 0 && is < phase_Sbins ) {
			data->radial_entropy[is][bin] += cell_volume[level];
		}
		if ( ig >= 0 && ig < phase_rhobins ) {
			data->radial_density[ig][bin] += cell_volume[level];
		}
#ifdef ELECTRON_ION_NONEQUILIBRIUM
		if ( ite >= 0 && ite < phase_Tbins ) {
			data->radial_electron_temperature[ite][bin] += cell_volume[level];
		}
#endif
	} else {
		for ( bin = 0; bin < num_bins; bin++ ) {
			bin_volume_fraction[bin] = 0;
		}

		/* throw random points */
		for ( i = 0; i < points_per_cell; i++ ) {
			for ( j = 0; j < nDim; j++ ) {
				pos[j] = cell->pos[j] + cell_size[level]*(cart_rand()-0.5);
			}

			r = compute_distance_periodic( halo->pos, pos ) / halo->rvir;
			bin = compute_radial_bin( r, rlmin, drl );

			if ( bin < num_bins ) {
				bin_volume_fraction[bin]++;
			}
		}

		for ( bin = 0; bin < num_bins; bin++ ){
			if ( bin_volume_fraction[bin] > 0 ) {
				volume_fraction = (double)bin_volume_fraction[bin] * cell_volume[level] /
					(double)points_per_cell;

				if ( it >= 0 && it < phase_Tbins ) {
					data->radial_temperature[it][bin] += volume_fraction;
				}
				if ( is >= 0 && is < phase_Sbins ) {
					data->radial_entropy[is][bin] += volume_fraction;
				}
				if ( ig >= 0 && ig < phase_rhobins ) {
					data->radial_density[ig][bin] += volume_fraction;
				}
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				if ( ite >= 0 && ite < phase_Tbins ) {
					data->radial_electron_temperature[ite][bin] += volume_fraction;
				}
#endif
			}
		}
	}
}

void write_table( char *filename, int num_bins1, float binlmin1, float binlmax1, 
			int num_bins2, float binlmin2, float binlmax2, float *data ) {
	int size;
	FILE *output;

	output = fopen( filename, "w" );
	size = sizeof(int);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &num_bins1, sizeof(int), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	size = 2*sizeof(float);
	fwrite( &size, sizeof(int), 1, output );
	fwrite( &binlmin1, sizeof(float), 1, output );
	fwrite( &binlmax1, sizeof(float), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	size = sizeof(int);
	fwrite( &size, sizeof(int), 1, output );
	fwrite( &num_bins2, sizeof(int), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	size = 2*sizeof(float);
	fwrite( &size, sizeof(int), 1, output );
	fwrite( &binlmin2, sizeof(float), 1, output );
	fwrite( &binlmax2, sizeof(float), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	size = num_bins1*num_bins2*sizeof(float);
	fwrite( &size, sizeof(int), 1, output );
	fwrite( data, sizeof(float), num_bins1*num_bins2, output );
	fwrite( &size, sizeof(int), 1, output );

	fclose(output);
}

void compute_halo_properties( char *output_directory, char *analysis_directory, halo_list *halos ) {
	int i, j;
	int ihalo;
	char filename1[256];
	char filename2[256];
	char filename3[256];
	FILE *output;
	float volume;
	int size;
	analysis_struct *data;

	Tfact = T0 * ( gamma - 1.0 ) / ( aexpn*aexpn );
        Pfact = P0/(aexpn*aexpn*aexpn*aexpn*aexpn);
        Sfact = S0; /* no trend with aexp if gamma = 5/3 */
	
	fact_nH = log10( 1.12e-5 * hubble * hubble * Omega0 * ( 1.0 - Y_p ) / aexpn );
        szfact = 3.4383e-15 * ( gamma - 1.0 ) * T0 * r0 * Omega0 * hubble / (aexpn*aexpn*aexpn*aexpn);

        rfact = 1000.0 * r0 * aexpn / hubble; /* proper kpc -> code units */
        vfact = v0 / aexpn;

        /* set up binning */
        rlmin = log10(rbinmin);
        rlmax = log10(rbinmax);
        drl = (rlmax - rlmin)/(float)(max_bins-1);
        num_bins = (rlmax - rlmin)/drl + 1;

	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		data = cart_alloc( sizeof(analysis_struct) );

		for ( i = 0; i < phase_Tbins; i++ ) {
			for ( j = 0; j < phase_rhobins; j++ ) {
				data->phase_table[i][j] = 0.0;
#ifdef ELECTRON_ION_NONEQULIBIRUM
				data->electron_phase_table[i][j] = 0.0;
#endif
			}
		}

		for ( i = 0; i < phase_Sbins; i++ ) {
			for ( j = 0; j < phase_rhobins; j++ ) {
				data->entropy_table[i][j] = 0.0;
			}
		}

		for ( i = 0; i < max_bins; i++ ) {
			for ( j = 0; j < phase_Tbins; j++ ) {
				data->radial_temperature[j][i] = 0.0;
			}
			for ( j = 0; j < phase_Sbins; j++ ) { 
				data->radial_entropy[j][i] = 0.0;
			}
			for ( j = 0; j < phase_rhobins; j++ ) {
				data->radial_density[j][i];
			}
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			for ( j = 0; j < phase_rhobins; j++ ) {
				data->radial_electron_temperature[j][i] = 0.0;
			}
#endif
		}

		halos->list[ihalo].data = data;
	}

	sprintf(filename1, "%s/%s_a%6.4f.grid", output_directory, jobname, aexpn );
	read_indexed_grid( filename1, halos, cell_callback );

	cart_debug("done reading grid");

	cart_debug("Tmin = %e", Tmin );

	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		data = (analysis_struct *)halos->list[ihalo].data;

		volume = (4./3.)*M_PI*pow(halos->list[ihalo].rvir,3.0);

		/* normalize to halo volumes */
		for ( i = 0; i < phase_Tbins; i++ ) {
			for ( j = 0; j < phase_rhobins; j++ ) {
				data->phase_table[i][j] /= volume;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				data->electron_phase_table[i][j] /= volume;
#endif
			}
		}

		for ( i = 0; i < phase_Sbins; i++ ) {
			for ( j = 0; j < phase_rhobins; j++ ) {
				data->entropy_table[i][j] /= volume;
			}
		}

		for ( i = 0; i < num_bins; i++ ) {
			if ( i == 0 ) {
				volume = 4.0*M_PI/3.0 * pow( rbinmin*halos->list[ihalo].rvir, 3.0);
			} else {
				volume = 4.0*M_PI/3.0 * pow( halos->list[ihalo].rvir, 3.0 ) * 
					( 	pow( 10., 3.0*(rlmin+(float)i*drl) ) - 
						pow( 10., 3.0*(rlmin+(float)(i-1)*drl) ) );
			}

			for ( j = 0; j < phase_Tbins; j++ ) {
				data->radial_temperature[j][i] /= volume;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				data->radial_electron_temperature[j][i] /= volume;
#endif
			}

			for ( j = 0; j < phase_Sbins; j++ ) {
				data->radial_entropy[j][i] /= volume;
			}

			for ( j = 0; j < phase_rhobins; j++ ) {
				data->radial_density[j][i] /= volume;
			}
		}
			

		sprintf( filename1, "%s/phase_a%6.4f_%04u.dat", analysis_directory,
				aexpn, halos->list[ihalo].id );
		write_table( filename1, phase_rhobins, phase_rhoglmin, phase_rhoglmax,
				phase_Tbins, phase_Tlmin, phase_Tlmax,
				(float *)data->phase_table );

		sprintf( filename1, "%s/entropy_a%6.4f_%04u.dat", analysis_directory,
				aexpn, halos->list[ihalo].id );
		write_table( filename1, phase_rhobins, phase_rhoglmin, phase_rhoglmax,
				phase_Sbins, phase_Slmin, phase_Slmax,
				(float *)data->entropy_table );

		sprintf( filename1, "%s/radial_temperature_a%6.4f_%04u.dat", analysis_directory,
                                aexpn, halos->list[ihalo].id );
                write_table( filename1, num_bins, rlmin, rlmax,
                                phase_Tbins, phase_Tlmin, phase_Tlmax, (float *)data->radial_temperature );

		sprintf( filename1, "%s/radial_entropy_%6.4f_%04u.dat", analysis_directory,
				aexpn, halos->list[ihalo].id );
		write_table( filename1, num_bins, rlmin, rlmax,
				phase_Sbins, phase_Slmin, phase_Slmax, (float *)data->radial_entropy);

		sprintf( filename1, "%s/radial_density_a%6.4f_%04u.dat", analysis_directory,
				aexpn, halos->list[ihalo].id );
		write_table( filename1, num_bins, rlmin, rlmax,
				phase_rhobins, phase_rhoglmin, phase_rhoglmax, (float *)data->radial_density );	

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		sprintf( filename1, "%s/electron_phase_a%6.4f_%04u.dat", analysis_directory,
				aexpn, halos->list[ihalo].id );
		write_table( filename1, phase_rhobins, phase_rhoglmin, phase_rhoglmax, 
				phase_Tbins, phase_Tlmin, phase_Tlmax,
				(float *)data->electron_phase_table );

		sprintf( filename1, "%s/radial_electron_temperature_a%6.4f_%04u.dat", analysis_directory,
				aexpn, halos->list[ihalo].id );
		write_table( filename1, num_bins, rlmin, rlmax,
				phase_Tbins, phase_Tlmin, phase_Tlmax, (float *)data->radial_electron_temperature );
#endif

	}
}
