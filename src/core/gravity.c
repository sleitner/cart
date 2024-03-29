#include "config.h"
#ifdef GRAVITY 

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "gravity.h"
#include "iterators.h"
#include "parallel.h"
#include "sfc.h"
#include "root_grid_fft.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "../tools/fft/fft3.h"

extern int step;

int num_initial_smooth_iterations = 1000;
int num_initial_smooth_steps = 0; 
int num_smooth_iterations = 60;   // used to be called MAX_SOR_ITER
float spectral_radius = 0.95;     // used to be called rhoJ


void root_grid_fft_get_cell_ijk(int cell, int ijk[nDim]);
void root_grid_fft_gravity_worker(const root_grid_fft_t *config, int id, fft_t *fft_source, fft_t *fft_output, int flags);


void config_init_gravity()
{
  control_parameter_add2(control_parameter_int,&num_initial_smooth_iterations,"gravity:num-initial-iterations","num_initial_smooth_iterations","number of iterations in the gravity relaxation solver for the first num_initial_smooth_steps global timesteps");

  control_parameter_add2(control_parameter_int,&num_initial_smooth_steps,"gravity:num-initial-steps","num_initial_smooth_steps","number of global timesteps to use gravity:num-initial-iterations rather than gravity:num-iterations in the relaxation solver");

  control_parameter_add3(control_parameter_int,&num_smooth_iterations,"gravity:num-iterations","num_smooth_iterations","MAX_SOR_ITER","number of iterations in the gravity relation solver (smooth)");

  control_parameter_add3(control_parameter_float,&spectral_radius,"gravity:spectral-radius","spectral_radius","rhoJ","Jacobi spectra radius for successful overrelation iterations.");
}

void config_verify_gravity()
{
  VERIFY(gravity:num-initial-iterations, num_initial_smooth_iterations > 0);

  VERIFY(gravity:num-initial-steps, num_initial_smooth_steps >= 0);

  VERIFY(gravity:num-iterations, num_smooth_iterations > 0 ); 

  VERIFY(gravity:spectral-radius, spectral_radius>0.0 && spectral_radius<1.0 );
}


void solve_poisson( int level, int flag ) {
	const int potential_vars[1] = { VAR_POTENTIAL };

	cart_assert( level >= min_level && level <= max_level );

	start_time( GRAVITY_TIMER );

	if ( level == min_level ) {
		root_grid_fft_exec(VAR_TOTAL_MASS,1,potential_vars,root_grid_fft_gravity_worker);
	} else {
		relax( level, flag );
	}

	end_time( GRAVITY_TIMER );
}

void relax( int level, int flag ) {

	if ( level > min_level && flag == 0 ) {
		prolongate( level );
	}
		
	smooth( level );
}

void prolongate( int level ) {
	int i,j;
	const int prolongation_vars[1] = { VAR_POTENTIAL };
	int icell;
	int num_level_cells;
	int *level_cells;

	cart_assert( level >= min_level+1 );

	start_time( WORK_TIMER );

	select_level( level-1, CELL_TYPE_LOCAL | CELL_TYPE_REFINED, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,j), shared(num_level_cells,level_cells,cell_vars,cell_child_oct)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		for ( j = 0; j < num_children; j++ ) {
			cell_potential( cell_child(icell,j) ) = cell_interpolate( icell, j, VAR_POTENTIAL );
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );

	/* update buffers */
	start_time( PROLONGATE_UPDATE_TIMER );
	update_buffer_level( level, prolongation_vars, 1 );
	end_time( PROLONGATE_UPDATE_TIMER );
}

void smooth( int level ) {
	int niter, iter;
	int i, j, k;
	int num_local_blocks;
	int num_direct_blocks;
	int num_indirect_blocks;
	int num_blocks;
	double phi0, phi1, phi2, phi3, phi4, phi5, phi6, phi7;
	double wsor, wsor6, trfi2;
	double *phi_red, *phi_black;
	double *rho_red, *rho_black;
	int *ext_red, *ext_black;
	double phi_ext_neighbors;
	int coords[nDim], coords_neighbor[nDim];
	int neighbor, direction, direct;
	int oct_count;
	int num_level_cells;
	int *level_cells;
	int *ind;
	int *oct_list;
	int num_red_border_cells;
	int num_black_border_cells;
	int list_size, cur_oct, block, child, icell, ioct;
	int sfc, proc;
	int num_send_octs[MAX_PROCS];
	int num_recv_octs[MAX_PROCS];
	int *send_indices[MAX_PROCS];
	int *recv_indices[MAX_PROCS];
	int *send_recv_indices[MAX_PROCS];
	double *buffer_red[MAX_PROCS];
	double *buffer_black[MAX_PROCS];
	double *packed_red[MAX_PROCS];
	double *packed_black[MAX_PROCS];
	int num_sendrequests;
	int num_recvrequests;
	MPI_Request requests[2*MAX_PROCS];
	MPI_Request send_requests[2*MAX_PROCS];
	MPI_Request recv_requests[2*MAX_PROCS];

	const int block_index[num_children] = { 
		#if nDim == 3 
			0, 0, 1, 1, 2, 2, 3, 3 
		#else
			#error "Unknown nDim in block_index (smooth)"
		#endif
	};

	const int color[num_children] = {
		#if nDim == 3
			0, 1, 1, 0, 1, 0, 0, 1
		#else
			#error "Unknown nDim in color (smooth)"
		#endif
	};

	const int smooth_vars[1] = { VAR_POTENTIAL };
	
	cart_assert ( level > min_level );

	start_time( SMOOTH_TIMER );
	start_time( SMOOTH_SETUP_TIMER );
	start_time( WORK_TIMER );

	/* create list of cell blocks */
	list_size = num_buffer_cells[level-1] + num_cells_per_level[level-1];

	oct_list = cart_alloc(int, list_size );

	oct_count = 0;
	num_local_blocks = 0;
	num_direct_blocks = 0;
	num_indirect_blocks = 0;

	for ( proc = 0; proc < num_procs; proc++ ) {
		num_recv_octs[proc] = 0;
		
		if ( num_local_buffers[level][proc] > 0 ) {
			recv_indices[proc] = cart_alloc(int, num_local_buffers[level][proc] );
			send_recv_indices[proc] = cart_alloc(int, num_local_buffers[level][proc] );
		}

		num_send_octs[proc] = num_remote_buffers[level][proc];

		/* recieve send list (yeah, I know...) */
		if ( num_send_octs[proc] > 0 ) {
			send_indices[proc] = cart_alloc(int, num_send_octs[proc] );

			MPI_Irecv( send_indices[proc], num_send_octs[proc], MPI_INT,
				proc, 0, mpi.comm.run, &requests[proc] );
		} else {
			requests[proc] = MPI_REQUEST_NULL;
		}
	}

	ind = cart_alloc(int, num_octs );

#pragma omp parallel for default(none), private(i), shared(ind,size_oct_array)
	for ( i = 0; i < num_octs; i++ ) {
		ind[i] = -1;
	}

	select_level( level-1, CELL_TYPE_LOCAL | CELL_TYPE_REFINED, &num_level_cells, &level_cells );

	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		oct_list[oct_count] = cell_child_oct[icell];
		ind[cell_child_oct[icell]] = oct_count;
		oct_count++;
		cart_assert( oct_count <= list_size );
	}

	num_local_blocks = oct_count;

	cart_free( level_cells );

	/* loop over buffer cells */
	select_level( level-1, CELL_TYPE_BUFFER | CELL_TYPE_REFINED, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		ioct = cell_child_oct[icell];
		sfc = cell_parent_root_sfc(icell);
		proc = processor_owner( sfc );

		recv_indices[proc][num_recv_octs[proc]] = ioct;
		send_recv_indices[proc][num_recv_octs[proc]] = 
			oct_hash_lookup( buffer_oct_reverse_hash[proc], ioct );
		num_recv_octs[proc]++;

		direct = 0;
		if ( level == min_level+1 ) {
			sfc_coords( sfc, coords );

			for ( j = 0; j < num_neighbors; j++ ) {
				for ( k = 0; k < nDim; k++ ) {
					coords_neighbor[k] = coords[k] + ishift[j][k];

					/* boundary check to ensure periodicity */
					if ( coords_neighbor[k] >= num_grid ) {
						coords_neighbor[k] = 0;
					} else if ( coords_neighbor[k] < 0 ) {
						coords_neighbor[k] = num_grid - 1;
					}
				}

				if ( root_cell_is_local( sfc_index( coords_neighbor ) ) ) {
					direct = 1;
					break;
				}
			}
		} else {
			/* using the oct neighbors should be a safe 
			 * optimization since we only
			 * solve gravity when oct neighbors are ok */
			for ( k = 0; k < num_neighbors; k++ ) {
				if ( oct_neighbors[ioct][k] == NULL_OCT ) {
					break;
				} else if ( cell_is_local( oct_neighbors[ioct][k] ) ) {
					direct = 1;
					break;
				}
			}
		}

		if ( direct ) {
			oct_list[oct_count] = ioct;
			ind[ioct] = oct_count;
			oct_count++;
			num_direct_blocks++;
			cart_assert( oct_count < list_size );
		} else {
			/* pack at the end */
			oct_list[list_size - num_indirect_blocks - 1] = ioct;
			num_indirect_blocks++;
		}
	}

	cart_free( level_cells );

	/* send receive list */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_recv_octs[proc] > 0 ) {
			MPI_Isend( send_recv_indices[proc], num_recv_octs[proc], MPI_INT,
				proc, 0, mpi.comm.run, &requests[num_procs+proc] );
		} else {
			requests[num_procs+proc] = MPI_REQUEST_NULL;
		}
	}

	num_blocks = num_local_blocks + num_direct_blocks + num_indirect_blocks;

	if ( num_indirect_blocks + oct_count < list_size ) {
		for ( i = 0; i < num_indirect_blocks; i++ ) {
			cart_assert( oct_count+i < list_size - num_indirect_blocks + i );
			cart_assert( oct_count+i < list_size );
			cart_assert( list_size - num_indirect_blocks + i < list_size );

			oct_list[oct_count+i] = oct_list[list_size - num_indirect_blocks + i];
			ind[oct_list[oct_count+i]] = oct_count+i;
		}
	} else {
#pragma omp parallel for default(none), private(i), shared(num_indirect_blocks,ind,oct_list,oct_count)
		for ( i = 0; i < num_indirect_blocks; i++ ) {
			ind[oct_list[oct_count+i]] = oct_count+i;
		}
	}

	/* do some sanity checks */
	for ( proc = 0; proc < num_procs; proc++ ) {
		cart_assert( num_recv_octs[proc] == num_local_buffers[level][proc] );
	}

	/* compute neighbors for cell in blocks */
	num_red_border_cells = 0;
	num_black_border_cells = 0;
	
	ext_red = cart_alloc(int, 4*oct_count*nDim );
	ext_black = cart_alloc(int, 4*oct_count*nDim );
	
	/* compute neighbors for local blocks */	
	for ( i = 0; i < num_local_blocks; i++ ) {
		cur_oct = oct_list[i];

		for ( j = 0; j < num_children; j++ ) {
			icell = 4*i+block_index[j];
		
			for ( k = 0; k < nDim; k++ ) {
				direction = external_direction[j][k];
				neighbor = oct_neighbors[cur_oct][direction];
				cart_assert( neighbor != NULL_OCT );
				cart_assert( cell_level(neighbor) == level-1 );

				if ( cell_is_refined(neighbor) ) {
					block = ind[cell_child_oct[neighbor]];
					cart_assert( block >= 0 && block < num_blocks );
					
					if ( color[j] == 0 ) {
						ext_red[nDim*icell+k] = 4*block + block_index[ local[j][direction] ];
					} else {
						ext_black[nDim*icell+k] = 4*block + block_index[ local[j][direction] ];
					}
				} else {
					if ( color[j] == 0 ) {
						ext_red[nDim*icell+k] = 4*num_blocks + num_black_border_cells;
						num_black_border_cells++;
					} else {
						ext_black[nDim*icell+k] = 4*num_blocks + num_red_border_cells;
						num_red_border_cells++;
					}
				}
			}
		}
	}

	/* compute neighbors for red direct blocks */
	for ( i = num_local_blocks; i < num_local_blocks+num_direct_blocks; i++ ) {
		cur_oct = oct_list[i];

		for ( j = 0; j < num_children; j++ ) {
			if ( color[j] == 0 ) {
				icell = 4*i+block_index[j];

				for ( k = 0; k < nDim; k++ ) {
					direction = external_direction[j][k];
					neighbor = oct_neighbors[cur_oct][direction];
					cart_assert( neighbor != NULL_OCT );
					cart_assert( cell_level(neighbor) == level-1 );

					if ( cell_is_local(neighbor) ) {
						break;
					}
				}

				if ( k < nDim ) { 
					for ( k = 0; k < nDim; k++ ) {
						direction = external_direction[j][k];
						neighbor = oct_neighbors[cur_oct][direction];
						
						cart_assert( neighbor != NULL_OCT );
						cart_assert( cell_level(neighbor) == level-1 );

						if ( cell_is_refined(neighbor) ) {
							block = ind[cell_child_oct[neighbor]];
							cart_assert( block >= 0 && block < num_blocks );

							ext_red[nDim*icell+k] = 4*block + block_index[ local[j][direction] ];
						} else {
							ext_red[nDim*icell+k] = 4*num_blocks + num_black_border_cells;
							num_black_border_cells++;
						}
					}
				} else {
					for ( k = 0; k < nDim; k++ ) {
						ext_red[nDim*icell+k] = 0;
					}
				}
			}
		}
	}

	/* now allocate space for cell values (since we know the size exactly) */
	phi_red = cart_alloc(double, (4*num_blocks + num_red_border_cells) );
	phi_black = cart_alloc(double, (4*num_blocks + num_black_border_cells) );

	num_red_border_cells = 0;
	num_black_border_cells = 0;

	/* compute border cell values */
	for ( i = 0; i < num_local_blocks; i++ ) {
		cur_oct = oct_list[i];

		for ( j = 0; j < num_children; j++ ) {
			icell = 4*i+block_index[j];

			for ( k = 0; k < nDim; k++ ) {
				direction = external_direction[j][k];
				neighbor = oct_neighbors[cur_oct][direction];
				cart_assert( neighbor != NULL_OCT );
				cart_assert( cell_level(neighbor) == level-1 );

				if ( cell_is_leaf(neighbor) ) {
					if ( color[j] == 0 ) {
						phi_black[4*num_blocks+num_black_border_cells] =
							cell_interpolate( neighbor, local[j][direction], VAR_POTENTIAL );
						num_black_border_cells++;
					} else {
						phi_red[4*num_blocks+num_red_border_cells] =
							cell_interpolate( neighbor, local[j][direction], VAR_POTENTIAL );
						num_red_border_cells++;
					}
				}
			}
		}
	}

	/* compute values for neighbors of red direct blocks */
	for ( i = num_local_blocks; i < num_local_blocks+num_direct_blocks; i++ ) {
		cur_oct = oct_list[i];

		for ( j = 0; j < num_children; j++ ) {
			if ( color[j] == 0 ) {
				icell = 4*i+block_index[j];

				for ( k = 0; k < nDim; k++ ) {
					direction = external_direction[j][k];
					neighbor = oct_neighbors[cur_oct][direction];
					cart_assert( neighbor != NULL_OCT );
					cart_assert( cell_level(neighbor) == level-1 );

					if ( cell_is_local(neighbor) ) {
						break;
					}
				}

				if ( k < nDim ) {
					for ( k = 0; k < nDim; k++ ) {
						direction = external_direction[j][k];
						neighbor = oct_neighbors[cur_oct][direction];
						cart_assert( neighbor != NULL_OCT );
						cart_assert( cell_level(neighbor) == level-1 );

						if ( cell_is_leaf(neighbor) ) {
							phi_black[4*num_blocks+num_black_border_cells] =
								cell_interpolate( neighbor, local[j][direction], VAR_POTENTIAL );
							num_black_border_cells++;
						}
					}
				}
			}
		}
	}

	rho_red = cart_alloc(double, 4*num_blocks );
	rho_black = cart_alloc(double, 4*num_blocks );

	/* pack cell_blocks */
#pragma omp parallel for default(none), private(i,j,child), shared(num_blocks,oct_list,phi_red,phi_black,cell_vars,rho_red,rho_black)
	for ( i = 0; i < num_blocks; i++ ) {
		j = 4*i;
		child = oct_child(oct_list[i],0);

		/* child 0 (red) */
		phi_red[j]              = cell_potential(child);
		rho_red[j]              = cell_total_mass(child);

		/* child 1 (black) */
		phi_black[j]            = cell_potential(child+1);
		rho_black[j]		= cell_total_mass(child+1);

		/* child 2 (black) */
		phi_black[j+1]          = cell_potential(child+2);
		rho_black[j+1]          = cell_total_mass(child+2);

		/* child 3 (red) */
		phi_red[j+1]            = cell_potential(child+3);
		rho_red[j+1]            = cell_total_mass(child+3);

		/* child 4 (black) */
		phi_black[j+2]          = cell_potential(child+4);
		rho_black[j+2]          = cell_total_mass(child+4);

		/* child 5 (red) */
		phi_red[j+2]            = cell_potential(child+5);
		rho_red[j+2]            = cell_total_mass(child+5);

		/* child 6 (red) */
		phi_red[j+3]            = cell_potential(child+6);
		rho_red[j+3]            = cell_total_mass(child+6);

		/* child 7 (black) */
		phi_black[j+3]          = cell_potential(child+7);
		rho_black[j+3]          = cell_total_mass(child+7);
	}

	end_time( WORK_TIMER );

	/* finish preparing send/recv indices arrays */
	start_time( COMMUNICATION_TIMER );

	start_time( SMOOTH_COMMUNICATION_TIMER ); 
	MPI_Waitall( 2*num_procs, requests, MPI_STATUSES_IGNORE );
	end_time( SMOOTH_COMMUNICATION_TIMER );

	/* convert send_indices and recv_indices from oct index to index in packed arrays */
	for ( proc = 0; proc < num_procs; proc++ ) {
#pragma omp parallel 
		{

#pragma omp for nowait
		for ( i = 0; i < num_send_octs[proc]; i++ ) {
			send_indices[proc][i] = 4*ind[send_indices[proc][i]];
		}

#pragma omp for
		for ( i = 0; i < num_recv_octs[proc]; i++ ) {
			recv_indices[proc][i] = 4*ind[ recv_indices[proc][i] ];
		}

		} /* end parallel region */

		if ( num_local_buffers[level][proc] > 0 ) {
			cart_free( send_recv_indices[proc] );
		}
	}

	cart_free( ind );

	/* create communications objects for iteration steps */
	num_sendrequests = 0;
	num_recvrequests = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_send_octs[proc] > 0 ) {
			packed_red[proc] = cart_alloc(double, (num_children/2)*num_send_octs[proc] );
			packed_black[proc] = cart_alloc(double, (num_children/2)*num_send_octs[proc] );

			MPI_Send_init( packed_red[proc], (num_children/2)*num_send_octs[proc], MPI_DOUBLE, proc, 0,
					mpi.comm.run, &send_requests[num_sendrequests++] );
			MPI_Send_init( packed_black[proc], (num_children/2)*num_send_octs[proc], MPI_DOUBLE, proc, 1, 
					mpi.comm.run, &send_requests[num_sendrequests++]);
		}
                                                                                                                                                            
		if ( num_recv_octs[proc] > 0 ) {
			buffer_red[proc] = cart_alloc(double, (num_children/2)*num_recv_octs[proc] );
			buffer_black[proc] = cart_alloc(double, (num_children/2)*num_recv_octs[proc] );
			
			MPI_Recv_init( buffer_red[proc], (num_children/2)*num_recv_octs[proc], MPI_DOUBLE,
					proc, 0, mpi.comm.run, &recv_requests[num_recvrequests++] );
			MPI_Recv_init( buffer_black[proc], (num_children/2)*num_recv_octs[proc], MPI_DOUBLE,
					proc, 1, mpi.comm.run, &recv_requests[num_recvrequests++] );
		}
	}

	end_time( COMMUNICATION_TIMER );	
	end_time( SMOOTH_SETUP_TIMER );

	/* work loop */
	wsor = 2.0; /* done for wsor_(1/2) calc, really == 1.0 */
	wsor6 = 1.0 / 6.0;
	trfi2 = units->potential * cell_size_inverse[level];

	niter = ( step < num_initial_smooth_steps ) ? 
				num_initial_smooth_iterations : num_smooth_iterations;

	for ( iter = 0; iter < niter; iter++ ) {
		if ( iter > 0 ) {
			start_time( COMMUNICATION_TIMER );

			start_time( SMOOTH_COMMUNICATION_TIMER );
			MPI_Startall( num_recvrequests, recv_requests );
			end_time( SMOOTH_COMMUNICATION_TIMER );

			/* pack into arrays for communication */
			for ( proc = 0; proc < num_procs; proc++ ) {
#pragma omp parallel for default(none), private(i,j,k), shared(proc,num_send_octs,send_indices,packed_red,phi_red,packed_black,phi_black)
				for ( i = 0; i < num_send_octs[proc]; i++ ) {
					j = send_indices[proc][i];
					k = (num_children/2)*i;

					packed_red[proc][k]	= phi_red[j];
					packed_red[proc][k+1]	= phi_red[j+1];
					packed_red[proc][k+2]	= phi_red[j+2];
					packed_red[proc][k+3]	= phi_red[j+3];

					packed_black[proc][k]     = phi_black[j];
					packed_black[proc][k+1]   = phi_black[j+1];
					packed_black[proc][k+2]   = phi_black[j+2];
					packed_black[proc][k+3]   = phi_black[j+3];
				}
			}

			/* do communication here */
			start_time( SMOOTH_COMMUNICATION_TIMER );
			MPI_Startall( num_sendrequests, send_requests );
			MPI_Waitall( num_recvrequests, recv_requests, MPI_STATUSES_IGNORE );
			end_time( SMOOTH_COMMUNICATION_TIMER );

			/* now move from buffers to actual array */
			for ( proc = 0; proc < num_procs; proc++ ) {
#pragma omp parallel for default(none), private(i,j,k), shared(num_recv_octs,proc,recv_indices,phi_red,buffer_red,phi_black,buffer_black)
				for ( i = 0; i < num_recv_octs[proc]; i++ ) {
					j = recv_indices[proc][i];
					k = (num_children/2)*i;

					phi_red[j]	= buffer_red[proc][k];
					phi_red[j+1]	= buffer_red[proc][k+1];
					phi_red[j+2]	= buffer_red[proc][k+2];
					phi_red[j+3]	= buffer_red[proc][k+3];

					phi_black[j]      = buffer_black[proc][k];
					phi_black[j+1]    = buffer_black[proc][k+1];
					phi_black[j+2]    = buffer_black[proc][k+2];
					phi_black[j+3]    = buffer_black[proc][k+3];
				}
			}

			end_time( COMMUNICATION_TIMER );
		}

		start_time( WORK_TIMER );

		/* red (local red and buffered red) */
#pragma omp parallel for default(none), private(i,block,phi1,phi2,phi4,phi7,phi_ext_neighbors,icell), shared(num_local_blocks,num_direct_blocks,phi_black,ext_red,phi_red,wsor6,trfi2,rho_red)
		for ( i = 0; i < num_local_blocks + num_direct_blocks; i++ ) {
			block = i*4;
			
			phi1 = phi_black[block];
			phi2 = phi_black[block+1];
			phi4 = phi_black[block+2];
			phi7 = phi_black[block+3];

			/* child 0 */
			icell = block;
			phi_ext_neighbors = phi_black[ext_red[nDim*icell]]
					+   phi_black[ext_red[nDim*icell+1]] 
					+   phi_black[ext_red[nDim*icell+2]];
			phi_red[icell] = phi_red[icell] 
				+ wsor6 * ( phi_ext_neighbors + phi1 + phi2 + phi4 - 6.0*phi_red[icell] )
				- trfi2 * rho_red[icell];

			/* child 3 */
			icell = icell+1;  /* block + 1 */
			phi_ext_neighbors = phi_black[ext_red[nDim*icell]]
					+   phi_black[ext_red[nDim*icell+1]]
					+   phi_black[ext_red[nDim*icell+2]];
			phi_red[icell] = phi_red[icell]
				+ wsor6 * ( phi_ext_neighbors + phi1 + phi2 + phi7 - 6.0*phi_red[icell] )
				- trfi2 * rho_red[icell];

			/* child 5 */
			icell = icell+1;  /* block + 2 */
			phi_ext_neighbors = phi_black[ext_red[nDim*icell]]
					+   phi_black[ext_red[nDim*icell+1]]
					+   phi_black[ext_red[nDim*icell+2]];
			phi_red[icell] = phi_red[icell]
				+ wsor6 * ( phi_ext_neighbors + phi1 + phi4 + phi7 - 6.0*phi_red[icell] )
				- trfi2 * rho_red[icell];

			/* child 6 */
			icell = icell+1;  /* block + 3 */
			phi_ext_neighbors = phi_black[ext_red[nDim*icell]]
					+   phi_black[ext_red[nDim*icell+1]]
					+   phi_black[ext_red[nDim*icell+2]];
			phi_red[icell] = phi_red[icell]
				+ wsor6 * ( phi_ext_neighbors + phi2 + phi4 + phi7 - 6.0*phi_red[icell] )
				- trfi2 * rho_red[icell];
		}

		/* chebyshev acceleration */
		wsor = 1.0 / ( 1.0 - 0.25 * spectral_radius*spectral_radius*wsor );
		wsor6 = wsor / 6.0;
		trfi2 = wsor * units->potential * cell_size_inverse[level];

		/* black (just local black cells) */
#pragma omp parallel for default(none), private(i,block,phi0,phi3,phi5,phi6,phi_ext_neighbors,icell), shared(num_local_blocks,phi_red,ext_black,phi_black,wsor6,trfi2,rho_black)
		for ( i = 0; i < num_local_blocks; i++ ) {
			block = i*4;

			phi0 = phi_red[block];
			phi3 = phi_red[block+1];
			phi5 = phi_red[block+2];
			phi6 = phi_red[block+3];

			/* child 1 */
			icell = block;
			phi_ext_neighbors = phi_red[ext_black[nDim*icell]]
				+   phi_red[ext_black[nDim*icell+1]]
				+   phi_red[ext_black[nDim*icell+2]];
			phi_black[icell] = phi_black[icell]
				+ wsor6 * ( phi_ext_neighbors + phi0 + phi3 + phi5 - 6.0*phi_black[icell] )
				- trfi2 * rho_black[icell];

			/* child 2 */
			icell = icell+1;  /* block + 1 */
			phi_ext_neighbors = phi_red[ext_black[nDim*icell]]
				+   phi_red[ext_black[nDim*icell+1]]
				+   phi_red[ext_black[nDim*icell+2]];
			phi_black[icell] = phi_black[icell]
				+ wsor6 * ( phi_ext_neighbors + phi0 + phi3 + phi6 - 6.0*phi_black[icell] )
				- trfi2 * rho_black[icell];

			/* child 4 */
			icell = icell+1;  /* block + 2 */
			phi_ext_neighbors = phi_red[ext_black[nDim*icell]]
				+   phi_red[ext_black[nDim*icell+1]]
				+   phi_red[ext_black[nDim*icell+2]];
			phi_black[icell] = phi_black[icell]
				+ wsor6 * ( phi_ext_neighbors + phi0 + phi5 + phi6 - 6.0*phi_black[icell] )
				- trfi2 * rho_black[icell];

			/* child 7 */
			icell = icell+1;  /* block + 3 */
			phi_ext_neighbors = phi_red[ext_black[nDim*icell]]
				+   phi_red[ext_black[nDim*icell+1]]
				+   phi_red[ext_black[nDim*icell+2]];
			phi_black[icell] = phi_black[icell]
				+ wsor6 * ( phi_ext_neighbors + phi3 + phi5 + phi6 - 6.0*phi_black[icell] )
				- trfi2 * rho_black[icell];
		}

		/* chebyshev acceleration */
		wsor = 1.0 / ( 1.0 - 0.25 * spectral_radius*spectral_radius*wsor );
		wsor6 = wsor / 6.0;
		trfi2 = wsor * units->potential * cell_size_inverse[level];

		end_time( WORK_TIMER );

		if ( iter > 0 ) {
			/* wait for packed sends to complete */
			start_time( COMMUNICATION_TIMER );
			start_time( SMOOTH_COMMUNICATION_TIMER );
			MPI_Waitall( num_sendrequests, send_requests, MPI_STATUSES_IGNORE );
			end_time( SMOOTH_COMMUNICATION_TIMER );
			end_time( COMMUNICATION_TIMER );
		}
	}

	/* write back to cell array */
	start_time( WORK_TIMER );

#pragma omp parallel for default(none), private(i,j,child), shared(num_local_blocks,oct_list,cell_vars,phi_red,phi_black,size_oct_array,size_cell_array)
	for ( i = 0; i < num_local_blocks; i++ ) {
		j = 4*i;
		cart_assert( oct_list[i] >= 0 && oct_list[i] < num_octs );
		child = oct_child(oct_list[i],0);
		cart_assert( child >= 0 && child < num_cells );
	
		/* child 0 (red) */
		cell_potential(child) = phi_red[j];

		/* child 1 (black) */
		cell_potential(child+1) = phi_black[j];

		/* child 2 (black) */
		cell_potential(child+2) = phi_black[j+1];

		/* child 3 (red) */
		cell_potential(child+3) = phi_red[j+1];

		/* child 4 (black) */
		cell_potential(child+4) = phi_black[j+2];

		/* child 5 (red) */
		cell_potential(child+5) = phi_red[j+2];

		/* child 6 (red) */
		cell_potential(child+6) = phi_red[j+3];

		/* child 7 (black) */
		cell_potential(child+7) = phi_black[j+3];
	}

	/* free MPI_Requests */
	for ( i = 0; i < num_sendrequests; i++ ) {
		MPI_Request_free( &send_requests[i] );
	}

	for ( i = 0; i < num_recvrequests; i++ ) {
		MPI_Request_free( &recv_requests[i] );
	}

	/* cleanup allocated memory */
	cart_free( ext_red );
	cart_free( ext_black );

	cart_free( phi_red );
	cart_free( phi_black );
	cart_free( rho_red );
	cart_free( rho_black );
	cart_free( oct_list );

	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_send_octs[proc] > 0 ) {
			cart_free( send_indices[proc] );
			cart_free( packed_red[proc] );
			cart_free( packed_black[proc] );
		}

		if ( num_recv_octs[proc] > 0 ) {
			cart_free( recv_indices[proc] );
			cart_free( buffer_red[proc] );
			cart_free( buffer_black[proc] );
		}
	}

	end_time( WORK_TIMER );
	end_time( SMOOTH_TIMER );

	/* update buffer cells */
	start_time( SMOOTH_UPDATE_TIMER );
	update_buffer_level( level, smooth_vars, 1 ); 
	end_time( SMOOTH_UPDATE_TIMER );
}

void restrict_to_level( int level ) {
	int i,j;
	double sum;
	const double rfw = 1.0 / (double)num_children;
	const int restrict_vars[1] = { VAR_POTENTIAL };

	int icell;
	int num_level_cells;
	int *level_cells;

	start_time( GRAVITY_TIMER );
	start_time( WORK_TIMER );
                                                                                                                  
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_REFINED, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,sum,j), shared(num_level_cells,level_cells,cell_child_oct,cell_vars)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		sum = 0.0;
		for ( j = 0; j < num_children; j++ ) {
			sum += cell_potential(cell_child(icell,j));
		}

		cell_potential(icell) = sum*rfw;
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );

	/* update buffers */
	start_time( RESTRICT_UPDATE_TIMER );
	update_buffer_level( level, restrict_vars, 1 );
	end_time( RESTRICT_UPDATE_TIMER );

	end_time( GRAVITY_TIMER );
}


void root_grid_fft_gravity_worker(const root_grid_fft_t *config, int id, fft_t *fft_source, fft_t *fft_output, int flags)
{
  int i, j, k, jk[2];
  double G, G_jk;
  double green[num_grid];
  double lambda;
  double trphi;
  size_t offset;

  cart_assert(config->bbox[0]==0 && config->bbox[1]==num_grid);

  trphi = -1.5*units->potential/(num_grid*num_grid*num_grid);

  /* precompute G(k) */
  lambda = M_PI/num_grid;
#pragma omp parallel for default(none), private(i), shared(green,lambda)
  for(i=0; i<num_grid; i++)
    {
      green[i] = sin(lambda*i)*sin(lambda*i);
    }

#pragma omp parallel for default(none), private(i,j,k,G,G_jk,offset,jk), shared(config,green,fft_source,fft_output,flags,trphi)
  for(k=0; k<config->dims[2]; k++)
    {
      for(j=0; j<config->dims[1]; j++)
	{
	  offset = fft3_jk_index(j,k,jk,flags);
	  if(offset == (size_t)-1) continue;

	  G_jk = green[jk[0]] + green[jk[1]];
	  offset *= config->dims[0];

	  for(i=0; i<=num_grid/2; i++)
	    {
	      if(i==0 && jk[0]==0 && jk[1]==0)
		{
		  G = 0.0;
		}
	      else
		{
		  G = trphi/(G_jk+green[i]);
		}

	      fft_output[2*i+0+offset] = G*fft_source[2*i+0+offset];
	      fft_output[2*i+1+offset] = G*fft_source[2*i+1+offset];
	    }
	}
    }
}

#endif /* GRAVITY */

