#include "config.h"

#ifdef COSMOLOGY

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <srfftw.h>
#include <sfftw.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "tree.h"
#include "timestep.h"
#include "units.h"

#include "power.h"

#define num_power_foldings	4
#define power_mesh_refinements	9		/* 512 mesh */
#define power_mesh_size		(1<<power_mesh_refinements)
#define num_power_mesh		(power_mesh_size*power_mesh_size*power_mesh_size)


int bin_from_d(float d) {
	/* NG: old Doug's binning */
	// return (int)(d-.5);
	/* NG: logarithmically spaced binning */
	return (int)(10*log10(d)+0.5);
}

int mesh_index( int ix, int iy, int iz ) {
	cart_assert( ix >= 0 && ix < power_mesh_size );
	cart_assert( iy >= 0 && iy < power_mesh_size );
	cart_assert( iz >= 0 && iz < power_mesh_size );
	return iz + power_mesh_size * ( iy + power_mesh_size * ix );
}

void compute_power_spectrum( char *filename, int power_type ) {
	int i, j, k, m;
	FILE *output;
        fftwnd_plan forward;
        fftw_complex *density_fft;
	int num_level_cells;
	int *level_cells;
	float *local_mesh, *global_mesh;
	double pos[nDim];
	int icell, level;
	float mass;
	double total;
	int mesh_level;
	int mesh_count;
	int ix, iy, iz;
	int ix1, iy1, iz1;
	int ix2, iy2, iz2;
	float xs, ys, zs;
	float dx0, dx1, dy0, dy1, dz0, dz1;
	float d00, d01, d10, d11;
	int ipart;
	float mesh_offset;
	float mesh_cell_size, mesh_cell_volume;
	int num_modes[power_mesh_size];
	float power[power_mesh_size];
	float avg_k[power_mesh_size];
	float di, dj, dk;
	float d, Pq, Pk, wk;
	float mass_factor, fb;
	double stellar_mass, total_stellar_mass;
	int bin, index;

	fb = cosmology->OmegaB / cosmology->OmegaM;
	mass_factor = ((float)(num_power_mesh)/(float)(num_grid*num_grid*num_grid));

	if ( power_type == POWER_TYPE_GAS || power_type == POWER_TYPE_STARS || power_type == POWER_TYPE_BARYONS ) {
#ifdef STARFORM
		if ( power_type != POWER_TYPE_BARYONS ) {
			/* compute stellar mass */
			stellar_mass = 0.0;
			for ( ipart = 0; ipart < num_star_particles; ipart++ ) {
				if ( particle_level[ipart] != FREE_PARTICLE_LEVEL && particle_is_star(ipart) ) {
					stellar_mass += particle_mass[ipart];
				}
			}

			MPI_Allreduce( &stellar_mass, &total_stellar_mass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

			cart_debug("total_stellar_mass = %e", total_stellar_mass );

			if ( power_type == POWER_TYPE_GAS ) {
				fb -= total_stellar_mass / (double)(num_grid*num_grid*num_grid);
			} else {
				fb = total_stellar_mass / (double)(num_grid*num_grid*num_grid);
			}
		}
#endif

		cart_debug("fb = %e", fb );
		mass_factor /= fb;
	} else if ( power_type == POWER_TYPE_DARK ) {
		mass_factor /= ( 1.0 - fb );
	} else if ( power_type == POWER_TYPE_TOTAL ) {
#ifndef HYDRO
		/* account for not loading hydro data but fb != 0 */
	        /* NG: this is actually incorrect, removing it */
		//mass_factor /= ( 1.0 - fb );
#endif /* HYDRO */
	}

	local_mesh = cart_alloc(float, num_power_mesh );

	if ( local_proc_id == MASTER_NODE ) {
		global_mesh = cart_alloc(float, num_power_mesh );
		density_fft = cart_alloc(fftw_complex, power_mesh_size*power_mesh_size*(power_mesh_size/2+1) );

		forward = rfftw3d_create_plan(power_mesh_size, power_mesh_size, power_mesh_size,
				FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE );

		output = fopen( filename, "w" );
		if ( output == NULL ) {
			cart_error("Unable to open %s for writing.", output );
		}
		fprintf( output, "aexp: %le Dplus[aexp]: %le  Dplus[1]: %le\n", auni[min_level], dplus_from_auni(auni[min_level]), dplus_from_auni(1.0) );
	}

	cart_debug("mass_factor = %e", mass_factor );

	for ( m = 0; m < num_power_foldings; m++ ) {
		mesh_cell_size = (float)num_grid / (float)(power_mesh_size << m);
		mesh_level = max( min_level, power_mesh_refinements - num_root_grid_refinements + m );
		mesh_cell_volume = mesh_cell_size*mesh_cell_size*mesh_cell_size;
	
		cart_debug("fold %u, mesh_level = %u, mesh_cell_size = %e, mesh_cell_volume = %e", m, mesh_level, mesh_cell_size, mesh_cell_volume );
	
		/* initialize mesh */
		if ( local_proc_id == MASTER_NODE ) {
			for ( i = 0; i < num_power_mesh; i++ ) {
				local_mesh[i] = -1.0;
			}
		} else {
			for ( i = 0; i < num_power_mesh; i++ ) {
				local_mesh[i] = 0.0;
			}
		}

		/* assign density */
#ifdef PARTICLES
		if ( power_type == POWER_TYPE_TOTAL || power_type == POWER_TYPE_DARK || power_type == POWER_TYPE_BARYONS ||
				power_type == POWER_TYPE_STARS ) {

			cart_debug("assigning particle density to mesh %u", m);
			for ( ipart = 0; ipart < num_particles; ipart++ ) {
				if ( particle_level[ipart] != FREE_PARTICLE_LEVEL && 
						( power_type == POWER_TYPE_TOTAL ||
#ifdef STARFORM
						( power_type == POWER_TYPE_DARK && !particle_is_star(ipart) ) ||
						( ( power_type == POWER_TYPE_BARYONS || power_type == POWER_TYPE_STARS ) && particle_is_star(ipart) ) ) ) {
#else
						( power_type == POWER_TYPE_DARK ) ) ) {
#endif /* STARFORM */
					xs = particle_x[ipart][0]/mesh_cell_size - 0.5;
					ys = particle_x[ipart][1]/mesh_cell_size - 0.5;
					zs = particle_x[ipart][2]/mesh_cell_size - 0.5;

					if(xs < 0) xs += power_mesh_size;
					if(ys < 0) ys += power_mesh_size;
					if(zs < 0) zs += power_mesh_size;

					ix = (int)(xs) % power_mesh_size;
					ix1 = (ix+1) % power_mesh_size;
					iy = (int)(ys) % power_mesh_size;
					iy1 = (iy+1) % power_mesh_size;
					iz = (int)(zs) % power_mesh_size;
					iz1 = (iz+1) % power_mesh_size;

					dx1 = xs - floor(xs);
					dy1 = ys - floor(ys);
					dz1 = zs - floor(zs);

					cart_assert( dx1 >= 0.0 && dx1 <= 1.0 );
					cart_assert( dy1 >= 0.0 && dy1 <= 1.0 );
					cart_assert( dz1 >= 0.0 && dz1 <= 1.0 );

					dx0 = 1.0 - dx1;
					dy0 = 1.0 - dy1;
					dz0 = 1.0 - dz1;

					dx0 *= particle_mass[ipart]*mass_factor;
					dx1 *= particle_mass[ipart]*mass_factor;

					d00 = dx0*dy0;
					d01 = dx0*dy1;
					d10 = dx1*dy0;
					d11 = dx1*dy1;
				
					local_mesh[mesh_index(ix,iy,iz)] += d00*dz0;
					local_mesh[mesh_index(ix,iy,iz1)] += d00*dz1;
					local_mesh[mesh_index(ix,iy1,iz)] += d01*dz0;
					local_mesh[mesh_index(ix,iy1,iz1)] += d01*dz1;
					local_mesh[mesh_index(ix1,iy,iz)] += d10*dz0;
					local_mesh[mesh_index(ix1,iy,iz1)] += d10*dz1;
					local_mesh[mesh_index(ix1,iy1,iz)] += d11*dz0;
					local_mesh[mesh_index(ix1,iy1,iz1)] += d11*dz1;
				}
			}
		}
#endif /* PARTICLES */

#ifdef HYDRO
		if ( power_type == POWER_TYPE_TOTAL || power_type == POWER_TYPE_BARYONS || power_type == POWER_TYPE_GAS ) {
			cart_debug("assigning gas density to mesh %u", m );
			for ( level = min_level; level <= min( mesh_level, max_level ); level++ ) {
				select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

				if ( level < mesh_level ) {
					mesh_count = 1 << ( mesh_level - level );
					mesh_offset = 0.5*(cell_size[level] - mesh_cell_size );
			
					for ( i = 0; i < num_level_cells; i++ ) {
						icell = level_cells[i];
	
						if ( cell_is_leaf(icell) ) {
							cell_center_position(icell,pos);
	
							/* find lower left mesh position */
							for ( j = 0; j < nDim; j++ ) {
								pos[j] -= mesh_offset;
							}
	
							ix1 = (int)(pos[0]/mesh_cell_size) % power_mesh_size;
							iy1 = (int)(pos[1]/mesh_cell_size) % power_mesh_size;
							iz1 = (int)(pos[2]/mesh_cell_size) % power_mesh_size;
	
							mass = cell_gas_density(icell)*mesh_cell_volume*mass_factor;
	
							for ( ix = ix1; ix < ix1 + mesh_count; ix++ ) {
								ix2 = ix % power_mesh_size;
								for ( iy = iy1; iy < iy1 + mesh_count; iy++ ) {
									iy2 = iy % power_mesh_size;
									for ( iz = iz1; iz < iz1 + mesh_count; iz++ ) {
										iz2 = iz % power_mesh_size;
										local_mesh[mesh_index(ix2,iy2,iz2)] += mass;
									}
								}
							}
						}
					}
				} else {
					/* assign this density even if it's not a leaf */
					for ( i = 0; i < num_level_cells; i++ ) {
						icell = level_cells[i];
	
						/* assign cell gas mass to mesh cell */
						cell_center_position(icell,pos);
	
						ix = (int)(pos[0]/mesh_cell_size) % power_mesh_size;
						iy = (int)(pos[1]/mesh_cell_size) % power_mesh_size; 
						iz = (int)(pos[2]/mesh_cell_size) % power_mesh_size;
	
						local_mesh[mesh_index(ix,iy,iz)] += 
							cell_gas_density(icell)*mesh_cell_volume*mass_factor;
					}
				}

				cart_free( level_cells );
			}
		}
#endif /* HYDRO */
		
		cart_debug("sending densities to master_node");

		/* merge */
		MPI_Reduce( local_mesh, global_mesh, num_power_mesh, MPI_FLOAT, 
				MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );

		/* compute power spectrum */
		if ( local_proc_id == MASTER_NODE ) {
			total = 0.0;
			for ( i = 0; i < num_power_mesh; i++ ) {
				total += global_mesh[i];
			}

			rfftwnd_one_real_to_complex( forward, (fftw_real *)global_mesh, density_fft );

			/* now average over modes */
			for ( i = 0; i < power_mesh_size; i++ ) {
				num_modes[i] = 0;
				avg_k[i] = 0.0;
				power[i] = 0.0;
			}

			for ( i = 0; i < power_mesh_size; i++ ) {
				if ( i <= power_mesh_size/2 ) {
					di = i*i;
				} else {
					di = (i-power_mesh_size)*(i-power_mesh_size);
				}

				for ( j = 0; j < power_mesh_size; j++ ) {
					if ( j <= power_mesh_size/2 ) {
						dj = j*j;
					} else {
						dj = (j-power_mesh_size)*(j-power_mesh_size);
					}

					/* skip 0,0,0 mode */
					if ( i != 0 || j != 0 ) {
						d = sqrt( di + dj );
						bin = bin_from_d(d);
						index = (power_mesh_size/2+1) * ( j + power_mesh_size * i );

						Pq = (density_fft[index].re*density_fft[index].re +
							density_fft[index].im * density_fft[index].im);

						if ( bin == 0 ) {
							cart_debug("%u %u %u: Pq = %e, k = %e, density_fft[%u] = %e %e", i, j, k, Pq, d,
								index, density_fft[index].re, density_fft[index].im );
						}

						power[bin] += Pq;

						avg_k[bin] += d;
						num_modes[bin]++;
					}

					/* add power from critical mode */
					d = sqrt( di + dj + (float)(power_mesh_size*power_mesh_size/4) );
					bin = bin_from_d(d);
					index = power_mesh_size/2 + (power_mesh_size/2+1) * ( j + power_mesh_size * i );

					Pq = (density_fft[index].re*density_fft[index].re +
						density_fft[index].im * density_fft[index].im);
					power[bin] += Pq;

					avg_k[bin] += d;
					num_modes[bin]++;

					for ( k = 1; k < power_mesh_size/2; k++ ) {
						dk = k*k;
						d = sqrt( di + dj + dk );
						bin = bin_from_d(d);

						index = k + (power_mesh_size/2+1) * ( j + power_mesh_size * i );

						Pq = (density_fft[index].re*density_fft[index].re +
							density_fft[index].im * density_fft[index].im);
						power[bin] += 2.0*Pq;

						avg_k[bin] += 2.0*d;
						num_modes[bin] += 2;
					}
				}
			}

			/* now write out modes */
			for ( i = 0; i < power_mesh_size; i++ ) {
				if ( num_modes[i] > 0 ) {
					wk = avg_k[i]/(float)num_modes[i] * (2.*M_PI*(float)(1<<m)/box_size);
					Pk = power[i]/(float)num_modes[i] * (box_size*box_size*box_size) / 
						(float)(power_mesh_size*power_mesh_size*power_mesh_size)/
						(float)(power_mesh_size*power_mesh_size*power_mesh_size);

					if ( i == 0 ) {
						cart_debug("power[%u] = %e", i, power[i] );
						cart_debug("power_mesh_size^3 = %e", (float)(power_mesh_size*power_mesh_size*power_mesh_size) );
						cart_debug("box_size^3 = %e", (box_size*box_size*box_size) );
					}

					fprintf(output, "%u %e %e %u\n", m, wk, Pk, num_modes[i] );
				}
			}
		}
	}

	/* output computed power spectrum */
	if ( local_proc_id == MASTER_NODE ) {
		cart_free( global_mesh );
		cart_free( density_fft );

		rfftwnd_destroy_plan(forward);

		fclose(output);
	}

	cart_free( local_mesh );
}

#endif /* COSMOLOGY */
