#include "config.h"

#ifdef COSMOLOGY

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "tree.h"
#include "times.h"
#include "units.h"

#include "fft/fft3.h"

#include "power.h"

#define num_power_foldings	4
#define power_mesh_refinements	9		/* 512 mesh */
#define power_mesh_size		(1<<power_mesh_refinements)


int bin_from_d(double d) {
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
	int num_level_cells;
	int *level_cells;
	fft_t *local_mesh, *global_mesh;
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
	double power[power_mesh_size];
	double avg_k[power_mesh_size];
	double di, dj, dk;
	double d, Pq, Pk, wk;
	float mass_factor, fb;
	double stellar_mass, total_stellar_mass;
	int bin, index;
	int pads[] = { 1, 0, 0 };
	int dims[3], bbox[6];
	int num_power_mesh, jk[2];
	size_t offset;

	fb = cosmology->OmegaB / cosmology->OmegaM;

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

			MPI_Allreduce( &stellar_mass, &total_stellar_mass, 1, MPI_DOUBLE, MPI_SUM, mpi.comm.run );

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

	if ( local_proc_id == MASTER_NODE ) {
		dims[0] = dims[1] = dims[2] = power_mesh_size;
		num_power_mesh = (size_t)dims[0]*dims[1]*dims[2];

		fft3_init(MPI_COMM_SELF,dims,pads,bbox);
		global_mesh = fft3_allocate_data();
		cart_assert(global_mesh != NULL);
		memset(global_mesh,0,sizeof(fft_t)*num_power_mesh);

		output = fopen( filename, "w" );
		if ( output == NULL ) {
			cart_error("Unable to open %s for writing.", output );
		}
		fprintf( output, "aexp: %le Dplus[aexp]: %le  Dplus[1]: %le\n", auni[min_level], dplus_from_auni(auni[min_level]), dplus_from_auni(1.0) );
	}

	local_mesh = cart_alloc(fft_t, num_power_mesh );

    mass_factor = ((float)(num_power_mesh)/(float)(num_grid*num_grid*num_grid));
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
				MPI_SUM, MASTER_NODE, mpi.comm.run );

		/* compute power spectrum */
		if ( local_proc_id == MASTER_NODE ) {
			total = 0.0;
			for ( i = 0; i < num_power_mesh; i++ ) {
				total += global_mesh[i];
			}

			fft3_x2k(global_mesh,FFT3_FLAG_WACKY_K_SPACE);

			/* now average over modes */
			for ( i = 0; i < power_mesh_size; i++ ) {
				num_modes[i] = 0;
				avg_k[i] = 0.0;
				power[i] = 0.0;
			}

			for ( k = 0; k < dims[2]; k++ ) {
				for ( j = 0; j < dims[1]; j++ ) {

				        offset = fft3_jk_index(j,k,jk,FFT3_FLAG_WACKY_K_SPACE);
				        if(offset == -1) continue;

                                        if ( k <= power_mesh_size/2 ) {
                                                dk = k*k;
                                        } else {
                                                dk = (k-power_mesh_size)*(k-power_mesh_size);
                                        }

                                        if ( j <= power_mesh_size/2 ) {
                                                dj = j*j;
                                        } else {
                                                dj = (j-power_mesh_size)*(j-power_mesh_size);
                                        }

					for ( i = 0; i < power_mesh_size/2; i++ ) {

						if(i==0 && jk[0]==0 && jk[1]==0) continue;

						di = i*i;
						d = sqrt( di + dj + dk );
						bin = bin_from_d(d);

						Pq = global_mesh[2*i+0+offset*dims[0]]*global_mesh[2*i+0+offset*dims[0]] + global_mesh[2*i+1+offset*dims[0]]*global_mesh[2*i+1+offset*dims[0]];

						power[bin] += 2.0*Pq;
						avg_k[bin] += 2.0*d;
						num_modes[bin] += 2;
					}

					i = power_mesh_size/2; {

						di = i*i;
						d = sqrt( di + dj + dk );
						bin = bin_from_d(d);

						Pq = global_mesh[2*i+0+offset*dims[0]]*global_mesh[2*i+0+offset*dims[0]] + global_mesh[2*i+1+offset*dims[0]]*global_mesh[2*i+1+offset*dims[0]];

						power[bin] += 1.0*Pq;
						avg_k[bin] += 1.0*d;
						num_modes[bin] += 1;
					}
				}
			}

			/* now write out modes */
			for ( i = 0; i < power_mesh_size; i++ ) {
				if ( num_modes[i] > 0 ) {
					wk = avg_k[i]/(double)num_modes[i] * (2.*M_PI*(double)(1<<m)/box_size);
					Pk = power[i]/(double)num_modes[i] * (box_size*box_size*box_size) / 
						(double)(power_mesh_size*power_mesh_size*power_mesh_size)/
						(double)(power_mesh_size*power_mesh_size*power_mesh_size);

					fprintf(output, "%u %le %le %u\n", m, wk, Pk, num_modes[i] );
				}
			}
		}
	}

	/* output computed power spectrum */
	if ( local_proc_id == MASTER_NODE ) {
		cart_free( global_mesh );
		fft3_done();

		fclose(output);
	}

	cart_free( local_mesh );
}

#endif /* COSMOLOGY */
