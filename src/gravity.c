#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "defs.h"
#include "tree.h"
#include "timestep.h"
#include "timing.h"
#include "gravity.h"
#include "iterators.h"
#include "cell_buffer.h"
#include "parallel.h"
#include "sfc.h"
#include "poisson.h"
#include "units.h"
#include "auxiliary.h"

#ifdef GRAVITY 

void solve_poisson( int level, int flag ) {

	cart_assert( level >= min_level && level <= max_level );

	start_time( GRAVITY_TIMER );

	if ( level == min_level ) {
		potential();
	} else {
		relax( level, flag );
	}

	end_time( GRAVITY_TIMER );
}

#ifdef HYDRO

void copy_potential( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_ANY, &num_level_cells, &level_cells );
	#pragma omp parallel for private(icell)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		cell_potential_hydro(icell) = cell_potential(icell);
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
}

void interpolate_potential( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;
	double dtdt2;

	start_time( WORK_TIMER );

	dtdt2 = 0.5 * dtl[level]/dtl_old[level];

	select_level( level, CELL_TYPE_ANY, &num_level_cells, &level_cells );
	#pragma omp parallel for private(icell)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		cell_potential_hydro(icell) = cell_potential(icell) +
			( cell_potential(icell) - cell_potential_hydro(icell) ) * dtdt2;
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
}
#endif /* HYDRO */

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

	select_level( level-1, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	#pragma omp parallel for private(icell,j)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_refined(icell) ) {
			for ( j = 0; j < num_children; j++ ) {
				cell_potential( cell_child(icell,j) ) = cell_interpolate( icell, j, VAR_POTENTIAL );
			}
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );

	/* update buffers */
	start_time( PROLONGATE_UPDATE_TIMER );
	update_buffer_level( level, prolongation_vars, 1 );
	end_time( PROLONGATE_UPDATE_TIMER );
}

#define MAX_SOR_ITERS	60
#define rhoJ		0.95

void smooth( int level ) {
	int iter;
	int i, j, k, m, n;
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
	int neighbors[num_neighbors];
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
	int *block_size;
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

	oct_list = cart_alloc( list_size*sizeof(int) );

	oct_count = 0;
	num_local_blocks = 0;
	num_direct_blocks = 0;
	num_indirect_blocks = 0;

	for ( proc = 0; proc < num_procs; proc++ ) {
		num_recv_octs[proc] = 0;
		
		if ( num_local_buffers[level][proc] > 0 ) {
			recv_indices[proc] = cart_alloc( num_local_buffers[level][proc]*sizeof(int) );
			send_recv_indices[proc] = cart_alloc( num_local_buffers[level][proc]*sizeof(int) );
		}

		num_send_octs[proc] = num_remote_buffers[level][proc];

		/* recieve send list (yeah, I know...) */
		if ( num_send_octs[proc] > 0 ) {
			send_indices[proc] = cart_alloc( num_send_octs[proc]*sizeof(int) );

			MPI_Irecv( send_indices[proc], num_send_octs[proc], MPI_INT,
				proc, 0, MPI_COMM_WORLD, &requests[proc] );
		} else {
			requests[proc] = MPI_REQUEST_NULL;
		}
	}

	ind = cart_alloc( num_octs * sizeof(int) );

	#pragma omp parallel for
	for ( i = 0; i < num_octs; i++ ) {
		ind[i] = -1;
	}

	select_level( level-1, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_refined(icell) ) {
			oct_list[oct_count] = cell_child_oct[icell];
			ind[cell_child_oct[icell]] = oct_count;
			oct_count++;
			cart_assert( oct_count <= list_size );
		}
	}

	num_local_blocks = oct_count;

	cart_free( level_cells );

	/* loop over buffer cells */
	select_level( level-1, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		if ( cell_is_refined(icell) ) {
			ioct = cell_child_oct[icell];
			sfc = cell_parent_root_sfc(icell);
			proc = processor_owner( sfc );

			recv_indices[proc][num_recv_octs[proc]] = ioct;
			send_recv_indices[proc][num_recv_octs[proc]] = 
					index_hash_lookup( buffer_oct_reverse_hash[proc], ioct );
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
	}

	cart_free( level_cells );

        /* send receive list */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_recv_octs[proc] > 0 ) {
			MPI_Isend( send_recv_indices[proc], num_recv_octs[proc], MPI_INT,
				proc, 0, MPI_COMM_WORLD, &requests[num_procs+proc] );
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
		#pragma omp parallel for
		for ( i = 0; i < num_indirect_blocks; i++ ) {
			ind[oct_list[oct_count+i]] = oct_count+i;
		}
	}

	/* do some sanity checks */
	for ( proc = 0; proc < num_procs; proc++ ) {
		cart_assert( num_recv_octs[proc] == num_local_buffers[level][proc] );
	}

	for ( i = 0; i < num_blocks; i++ ) {
		cart_assert( ind[ oct_list[i] ] != NULL_OCT );
	}

	/* compute neighbors for cell in blocks */
	num_red_border_cells = 0;
	num_black_border_cells = 0;
	
	ext_red = cart_alloc( 4*oct_count*nDim * sizeof(int) );
	ext_black = cart_alloc( 4*oct_count*nDim * sizeof(int) );
	
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
        phi_red = cart_alloc( (4*num_blocks + num_red_border_cells) * sizeof(double) );
        phi_black = cart_alloc( (4*num_blocks + num_black_border_cells) * sizeof(double) );

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

	rho_red = cart_alloc( 4*num_blocks * sizeof(double) );
	rho_black = cart_alloc( 4*num_blocks * sizeof(double) );

	/* pack cell_blocks */
	#pragma omp parallel for private(j,child)
	for ( i = 0; i < num_blocks; i++ ) {
		j = 4*i;
		child = oct_child(oct_list[i],0);

		/* child 0 (red) */
		phi_red[j]              = cell_potential(child);
		rho_red[j]              = cell_density(child);

		/* child 1 (black) */
		phi_black[j]            = cell_potential(child+1);
		rho_black[j]		= cell_density(child+1);

		/* child 2 (black) */
		phi_black[j+1]          = cell_potential(child+2);
		rho_black[j+1]          = cell_density(child+2);

		/* child 3 (red) */
		phi_red[j+1]            = cell_potential(child+3);
		rho_red[j+1]            = cell_density(child+3);

		/* child 4 (black) */
		phi_black[j+2]          = cell_potential(child+4);
		rho_black[j+2]          = cell_density(child+4);

		/* child 5 (red) */
		phi_red[j+2]            = cell_potential(child+5);
		rho_red[j+2]            = cell_density(child+5);

		/* child 6 (red) */
		phi_red[j+3]            = cell_potential(child+6);
		rho_red[j+3]            = cell_density(child+6);

		/* child 7 (black) */
		phi_black[j+3]          = cell_potential(child+7);
		rho_black[j+3]          = cell_density(child+7);
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

		#pragma omp for, nowait
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
			packed_red[proc] = cart_alloc( (num_children/2)*num_send_octs[proc]*sizeof(double) );
			packed_black[proc] = cart_alloc( (num_children/2)*num_send_octs[proc]*sizeof(double) );

			MPI_Send_init( packed_red[proc], (num_children/2)*num_send_octs[proc], MPI_DOUBLE, proc, 0,
				MPI_COMM_WORLD, &send_requests[num_sendrequests++] );
			MPI_Send_init( packed_black[proc], (num_children/2)*num_send_octs[proc], MPI_DOUBLE, proc, 1, 
				MPI_COMM_WORLD, &send_requests[num_sendrequests++]);
                }
                                                                                                                                                            
                if ( num_recv_octs[proc] > 0 ) {
			buffer_red[proc] = cart_alloc( (num_children/2)*num_recv_octs[proc]*sizeof(double) );
			buffer_black[proc] = cart_alloc( (num_children/2)*num_recv_octs[proc]*sizeof(double) );
			
                        MPI_Recv_init( buffer_red[proc], (num_children/2)*num_recv_octs[proc], MPI_DOUBLE,
                                proc, 0, MPI_COMM_WORLD, &recv_requests[num_recvrequests++] );
                        MPI_Recv_init( buffer_black[proc], (num_children/2)*num_recv_octs[proc], MPI_DOUBLE,
                                proc, 1, MPI_COMM_WORLD, &recv_requests[num_recvrequests++] );
                }
        }

	end_time( COMMUNICATION_TIMER );	
	end_time( SMOOTH_SETUP_TIMER );

	/* work loop */
	wsor = 2.0; /* done for wsor_(1/2) calc, really == 1.0 */
	wsor6 = 1.0 / 6.0;
	trfi2 = cell_size_inverse[level] / aexp[level]; 

	for ( iter = 0; iter < MAX_SOR_ITERS; iter++ ) {
		if ( iter > 0 ) {
			start_time( COMMUNICATION_TIMER );

			start_time( SMOOTH_COMMUNICATION_TIMER );
			MPI_Startall( num_recvrequests, recv_requests );
			end_time( SMOOTH_COMMUNICATION_TIMER );

			/* pack into arrays for communication */
			for ( proc = 0; proc < num_procs; proc++ ) {
				#pragma omp parallel for private(j,k)
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
				#pragma omp parallel for private(j,k)
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
		#pragma omp parallel for private(block,phi1,phi2,phi4,phi7,phi_ext_neighbors,icell)
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
		wsor = 1.0 / ( 1.0 - 0.25 * rhoJ*rhoJ*wsor );
		wsor6 = wsor / 6.0;
		trfi2 = wsor * cell_size_inverse[level] / aexp[level]; 

		/* black (just local black cells) */
		#pragma omp parallel for private(block,phi0,phi3,phi5,phi6,phi_ext_neighbors,icell)
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
		wsor = 1.0 / ( 1.0 - 0.25 * rhoJ*rhoJ*wsor );
		wsor6 = wsor / 6.0;
		trfi2 = wsor * cell_size_inverse[level] / aexp[level]; 

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

	#pragma omp parallel for private(j,child)
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

	start_time( WORK_TIMER );
                                                                                                                  
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	#pragma omp parallel for private(icell,sum,j)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		if ( cell_is_refined(icell) ) {
			sum = 0.0;
			for ( j = 0; j < num_children; j++ ) {
				sum += cell_potential(cell_child(icell,j));
			}

			cell_potential(icell) = sum*rfw;
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );

	/* update buffers */
	start_time( RESTRICT_UPDATE_TIMER );
	update_buffer_level( level, restrict_vars, 1 );
	end_time( RESTRICT_UPDATE_TIMER );
}

/* compute potential on min_level (WILL BE REWRITTEN!) */
void potential() {
	const int potential_vars[1] = { VAR_POTENTIAL };
	
	int coords[3];
	int index;
	int p, i;
	type_fft *density, *potential;
	MPI_Status status;

	start_time( FFT_TIMER );

	/* gather all root densities to master node */
	if ( local_proc_id == MASTER_NODE ) {
		start_time( WORK_TIMER );

		potential = cart_alloc( num_root_cells * sizeof(type_fft) );
		density = cart_alloc( num_root_cells * sizeof(type_fft) );

		#pragma omp parallel for private(coords,index)
		for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
			sfc_coords( i, coords );
			index = num_grid*(num_grid*coords[0] + coords[1] ) + coords[2];
			density[index] = cell_density(i);
		}

		end_time( WORK_TIMER );

		start_time( COMMUNICATION_TIMER );

		for ( p = 1; p < num_procs; p++ ) {
			MPI_Recv( potential, proc_sfc_index[p+1]-proc_sfc_index[p], MPI_TYPE_FFT,
				p, 0, MPI_COMM_WORLD, &status );

			#pragma omp parallel for private(coords,index)
			for ( i = 0; i < proc_sfc_index[p+1]-proc_sfc_index[p]; i++ ) {
				sfc_coords(i+proc_sfc_index[p],coords);
				index = num_grid*(num_grid*coords[0] + coords[1] ) + coords[2];
				density[index] = potential[i];
			}
		}

		end_time( COMMUNICATION_TIMER );

		/* compute potential */
		start_time( WORK_TIMER );
	        poisson( density, potential );
		end_time( WORK_TIMER );

		start_time( COMMUNICATION_TIMER );
		for ( p = 1; p < num_procs; p++ ) {
			#pragma omp parallel for private(coords,index)
			for ( i = 0; i < proc_sfc_index[p+1]-proc_sfc_index[p]; i++ ) {
				sfc_coords(i+proc_sfc_index[p],coords);
				index = num_grid*(num_grid*coords[0] + coords[1] ) + coords[2];
				density[i] = potential[index];
			}

			MPI_Send( density, proc_sfc_index[p+1]-proc_sfc_index[p], MPI_TYPE_FFT, p, 1, MPI_COMM_WORLD );
		}
		end_time( COMMUNICATION_TIMER );

		start_time( WORK_TIMER );
		#pragma omp parallel for private(coords,index)
		for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
			sfc_coords( i, coords );
			index = num_grid*(num_grid*coords[0] + coords[1] ) + coords[2];

			cell_potential(i) = potential[index];
		}
		end_time( WORK_TIMER );

		cart_free(potential);
		cart_free(density);
	} else {
		start_time( COMMUNICATION_TIMER );
		potential = cart_alloc( num_cells_per_level[min_level] * sizeof(type_fft) );
		
		#pragma omp parallel for
		for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
			potential[i] = cell_density(i);
		}
		
		MPI_Send( potential, num_cells_per_level[min_level], MPI_TYPE_FFT, MASTER_NODE, 0, MPI_COMM_WORLD );
		MPI_Recv( potential, num_cells_per_level[min_level], MPI_TYPE_FFT, MASTER_NODE, 1, MPI_COMM_WORLD, &status );

		#pragma omp parallel for
		for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
			cell_potential(i) = potential[i];
		}

		cart_free(potential);
		end_time( COMMUNICATION_TIMER );
	}

	/* update cell buffer */
	start_time( FFT_UPDATE_TIMER );
	update_buffer_level( min_level, potential_vars, 1 );
	end_time( FFT_UPDATE_TIMER );

	end_time( FFT_TIMER );
}

#ifdef HYDRO 
void compute_accelerations_hydro( int level ) {
	int i, j;
	double a2half;
	const int accel_vars[nDim] = { VAR_ACCEL, VAR_ACCEL+1, VAR_ACCEL+2 };
	int neighbors[num_neighbors];
	int L1, R1;
	double phi_l, phi_r;

	int icell;
	int num_level_cells;
	int *level_cells;

	start_time( HYDRO_ACCEL_TIMER );

#ifdef COSMOLOGY
	a2half = b2a( tl[level] + 0.5*dtl[level] );
	a2half = a2half*a2half * -0.5*dtl[level]*cell_size_inverse[level];
#else
	a2half = 0.5 * dtl[level] * cell_size_inverse[level];
#endif 

	start_time( WORK_TIMER );

        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        #pragma omp parallel for private(icell,j,neighbors,L1,R1,phi_l,phi_r)
        for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		cell_all_neighbors( icell, neighbors );
		for ( j = 0; j < nDim; j++ ) {
			L1 = neighbors[2*j];
			R1 = neighbors[2*j+1];

			if ( cell_level(L1) == level && cell_level(R1) == level ) {
				phi_l = cell_potential_hydro(L1);
				phi_r = cell_potential_hydro(R1);
			} else {
				if ( cell_level(L1) < level ) {
					phi_l = cell_interpolate( L1, local[cell_child_number(icell)][2*j], VAR_POTENTIAL );
					phi_r = cell_potential(R1);
				} else {
					phi_l = cell_potential(L1);
					phi_r = cell_interpolate( R1, local[cell_child_number(icell)][2*j], VAR_POTENTIAL );
				}
			}

			cell_accel( icell, j ) = (float)(a2half * ( phi_r - phi_l ) );
                }
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );

#ifdef GRAVITY_IN_RIEMANN
	/* this only gets called if we pass gravity on to Riemann solver (and thus need accel in buffer cells) */
	start_time( HYDRO_ACCEL_UPDATE_TIMER );
	update_buffer_level( level, accel_vars, nDim );
	end_time( HYDRO_ACCEL_UPDATE_TIMER );
#endif

	end_time( HYDRO_ACCEL_TIMER );
}
#endif /* HYDRO */

#ifdef PARTICLES

void compute_accelerations_particles( int level ) {
	int i, j;
	double a2half;
	const int accel_vars[nDim] = { VAR_ACCEL, VAR_ACCEL+1, VAR_ACCEL+2 };
	int neighbors[num_neighbors];
	int L1, R1;
	double phi_l, phi_r;
	int icell;
	int num_level_cells;
	int *level_cells;

	start_time( PARTICLE_ACCEL_TIMER );
        
	start_time( WORK_TIMER );

#ifdef COSMOLOGY
	a2half = -0.5*aexp[level]*aexp[level]*cell_size_inverse[level];
#else
	a2half = -0.5 * cell_size_inverse[level];
#endif

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	#pragma omp parallel for private(icell,j,neighbors,L1,R1,phi_l,phi_r)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		cell_all_neighbors( icell, neighbors );
		for ( j = 0; j < nDim; j++ ) {
			L1 = neighbors[2*j];
			R1 = neighbors[2*j+1];

			if ( cell_level(L1) < level ) {
				phi_l = 0.8*cell_potential(L1) + 0.2*cell_potential(icell);
			} else {
				phi_l = cell_potential(L1);
			}

			if ( cell_level(R1) < level ) {
				phi_r = 0.8*cell_potential(R1)+0.2*cell_potential(icell);
			} else {
				phi_r = cell_potential(R1);
			}

			cell_accel( icell, j ) = (float)(a2half * ( phi_r - phi_l ) );
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );

	/* update accelerations */
	start_time( PARTICLE_ACCEL_UPDATE_TIMER );
	update_buffer_level( level, accel_vars, nDim );
	end_time( PARTICLE_ACCEL_UPDATE_TIMER );

	end_time( PARTICLE_ACCEL_TIMER );
}
#endif /* PARTICLES */

#endif /* GRAVITY */
