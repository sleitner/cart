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

#include "../tools/fft/fft3.h"

#include "power.h"

int num_power_foldings = 1;
int power_mesh_size = num_grid;

int bin_from_d(double d) {
	/* NG: old Doug's binning */
	// return (int)(d-.5);
	/* NG: logarithmically spaced binning */
	return (int)(20*log10(d)+0.5);
}

size_t mesh_index( int ix, int iy, int iz ) {
  return iz + (power_mesh_size+2) * ( iy + (size_t)power_mesh_size * ix );
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
	double xs, ys, zs;
        float dx0, dx1, dy0, dy1, dz0, dz1;
	float d00, d01, d10, d11;
	int ipart;
	float mesh_offset;
	float mesh_cell_size, mesh_cell_volume;
	double di, dj, dk;
	double d, Pq, Pk, wk;
	int bin;
	int dims[nDim], bbox[2*nDim];
	int num_power_mesh, num_power_mesh2, jk[2];
	size_t offset;
	int num_bins;
	int *num_modes;
	double *power, *avg_k;
	int power_mesh_refinements = 0;

	i = power_mesh_size;
	while ( i > 1 ) {
		cart_assert( (i & 1) == 0 || (i >> 1) == 0 ); /* check for power of two */
		power_mesh_refinements++;
		i >>= 1;
	}

	dims[0] = dims[1] = dims[2] = power_mesh_size;
	num_power_mesh = (size_t)dims[0]*dims[1]*dims[2];
	num_power_mesh2 = (size_t)(dims[0]+2)*dims[1]*dims[2];

	num_bins = bin_from_d(1.0*power_mesh_size) + 1;
	cart_debug("Using %d bins",num_bins);

	num_modes = cart_alloc(int,num_bins);
	avg_k = cart_alloc(double,num_bins);
	power = cart_alloc(double,num_bins);

	if ( local_proc_id == MASTER_NODE ) {
		int pads[] = { 0, 0, 0 };
		fft3_init(MPI_COMM_SELF,dims,pads,bbox);
		global_mesh = fft3_allocate_data();
		cart_assert(global_mesh != NULL);
		memset(global_mesh,0,sizeof(fft_t)*num_power_mesh2);

		cart_debug("num_power_mesh = %d",num_power_mesh);
		cart_debug("dims = %d %d %d",dims[0],dims[1],dims[2]);
		cart_debug("bbox = %d %d %d %d %d %d",bbox[0],bbox[1],bbox[2],bbox[3],bbox[4],bbox[5]);

		output = fopen( filename, "w" );
		if ( output == NULL ) {
			cart_error("Unable to open %s for writing.", filename );
		}
		fprintf( output, "# aexp: %le Dplus[aexp]: %le  Dplus[1]: %le\n", auni[min_level], dplus_from_auni(auni[min_level]), dplus_from_auni(1.0) );
	}

	local_mesh = cart_alloc(fft_t, num_power_mesh2 );

	for ( m = 0; m < num_power_foldings; m++ ) {
		mesh_cell_size = (float)num_grid / (float)(power_mesh_size << m);
		mesh_level = MAX( min_level, power_mesh_refinements - num_root_grid_refinements + m );
		mesh_cell_volume = mesh_cell_size*mesh_cell_size*mesh_cell_size;
	
		cart_debug("fold %u, mesh_level = %u, mesh_cell_size = %e, mesh_cell_volume = %e", m, mesh_level, mesh_cell_size, mesh_cell_volume );
	
		/* initialize mesh */
		for ( i = 0; i < num_power_mesh2; i++ ) {
			local_mesh[i] = 0.0;
		}

		/* assign density */
#ifdef PARTICLES
		if ( power_type == POWER_TYPE_TOTAL || power_type == POWER_TYPE_DARK || power_type == POWER_TYPE_BARYONS ||
				power_type == POWER_TYPE_STARS ) {

			cart_debug("assigning particle density to mesh %u", m);
			for ( ipart = 0; ipart < num_particles; ipart++ ) {
				if ( particle_level[ipart] != FREE_PARTICLE_LEVEL && 
						( power_type == POWER_TYPE_TOTAL ||
#ifdef STAR_FORMATION
						( power_type == POWER_TYPE_DARK && !particle_is_star(ipart) ) ||
						( ( power_type == POWER_TYPE_BARYONS || power_type == POWER_TYPE_STARS ) && particle_is_star(ipart) ) ) 
#else
						( power_type == POWER_TYPE_DARK ) ) 
#endif /* STAR_FORMATION */
					) {

					xs = particle_x[ipart][0]/mesh_cell_size - 0.5;
					ys = particle_x[ipart][1]/mesh_cell_size - 0.5;
					zs = particle_x[ipart][2]/mesh_cell_size - 0.5;

					if(xs < 0) xs += (float)(power_mesh_size << m);
					if(ys < 0) ys += (float)(power_mesh_size << m);
					if(zs < 0) zs += (float)(power_mesh_size << m);

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

					dx0 *= particle_mass[ipart];
					dx1 *= particle_mass[ipart];

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
			for ( level = min_level; level <= MIN( mesh_level, max_level ); level++ ) {
				select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

				/*
				//  NG: for large mesh_level - level the original assignment was extremely
				//  inefficient, this is a more optimized version
				*/
#ifdef DEBUG
				cart_debug("... assigning density to level=%d at mesh_level=%d",level,mesh_level);
#endif

				if ( level < mesh_level ) {
				        float *wrap_buf;   // BTW: this is legal C89
					int wrap_size;

					mesh_count = 1 << ( mesh_level - level );
					mesh_offset = 0.5*(cell_size[level] - mesh_cell_size );
		    
					wrap_size = (power_mesh_size+mesh_count-1)/mesh_count;  // just in case
					wrap_buf = cart_alloc(fft_t,wrap_size*wrap_size*wrap_size);
					memset(wrap_buf,0,sizeof(float)*wrap_size*wrap_size*wrap_size);

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
	
							ix2 = ix1 / mesh_count;
							iy2 = iy1 / mesh_count;
							iz2 = iz1 / mesh_count;

							wrap_buf[iz2+wrap_size*(iy2+wrap_size*ix2)] += cell_gas_density(icell)*mesh_cell_volume;
						}
					}

					for ( ix = 0; ix < power_mesh_size; ix++ ) {
						ix2 = ix / mesh_count;
						for ( iy = 0; iy < power_mesh_size; iy++ ) {
							iy2 = iy / mesh_count;
							for ( iz = 0; iz < power_mesh_size; iz++ ) {
								iz2 = iz / mesh_count;
								local_mesh[mesh_index(ix,iy,iz)] += wrap_buf[iz2+wrap_size*(iy2+wrap_size*ix2)];
							}
						}
					}

					cart_free(wrap_buf);

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
							cell_gas_density(icell)*mesh_cell_volume;
					}
				}

				cart_free( level_cells );
			}
		}
#endif /* HYDRO */
		
		cart_debug("sending densities to master_node");

		/* merge */
		MPI_Reduce( local_mesh, global_mesh, num_power_mesh2, MPI_TYPE_FFT, 
				MPI_SUM, MASTER_NODE, mpi.comm.run );

		cart_debug("densities have been sent");

		/* compute power spectrum */
		if ( local_proc_id == MASTER_NODE ) {
			/* properly normalize mesh */
			total = 0.0;
			for ( ix = 0; ix < dims[2]; ix++ ) {
			  for ( iy = 0; iy < dims[1]; iy++ ) {
			    for ( iz = 0; iz < dims[0]; iz++ ) {
				total += global_mesh[mesh_index(ix,iy,iz)];
			    }
			  }
			}
			total /= (double)num_power_mesh; 

			if ( m == 0 ) {
			  fprintf( output, "# Omega (<rho>/rho_m) = %le\n", total*(double)num_power_mesh/pow(num_grid,3) );
			}

			/* convert to overdensity */
			for ( i = 0; i < num_power_mesh2; i++ ) {
				global_mesh[i] = global_mesh[i]/total - 1.0;
			}

			fft3_x2k(global_mesh,FFT3_FLAG_WACKY_K_SPACE);

			/* now average over modes */
			for ( i = 0; i < num_bins; i++ ) {
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

						cart_assert(bin>=0 && bin<num_bins);

						Pq = global_mesh[2*i+0+offset*dims[0]]*global_mesh[2*i+0+offset*dims[0]] + global_mesh[2*i+1+offset*dims[0]]*global_mesh[2*i+1+offset*dims[0]];

						power[bin] += 2.0*Pq;
						avg_k[bin] += 2.0*d;
						num_modes[bin] += 2;
					}

					i = power_mesh_size/2; {

						di = i*i;
						d = sqrt( di + dj + dk );
						bin = bin_from_d(d);

						cart_assert(bin>=0 && bin<num_bins);

						Pq = global_mesh[2*i+0+offset*dims[0]]*global_mesh[2*i+0+offset*dims[0]] + global_mesh[2*i+1+offset*dims[0]]*global_mesh[2*i+1+offset*dims[0]];

						power[bin] += 1.0*Pq;
						avg_k[bin] += 1.0*d;
						num_modes[bin] += 1;
					}
				}
			}

			/* now write out modes */
			for ( i = 0; i < num_bins; i++ ) {
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
	cart_free( num_modes );
	cart_free( power );
	cart_free( avg_k );
}

#endif /* COSMOLOGY */
