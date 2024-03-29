#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "hydro_tracer.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "sfc.h"
#include "starformation.h"
#include "tree.h"


void check_map() {
	int i, j, k;
	double pos[nDim];
	int neighbors[num_neighbors];
	int total_root_cells;
	int level;
	int icell, ioct;
	int num_level_cells;
	int *level_cells;
	int sfc;
	float max_var[num_vars];
	float min_var[num_vars];

#ifdef PARTICLES
	int ipart, ipart_next;
	int total_local_particles;
	particleid_t total_particles, tmp;
	particleid_t species_count[MAX_PARTICLE_SPECIES];
	particleid_t species_count_total[MAX_PARTICLE_SPECIES];
#endif /* PARTICLES */
#ifdef GRAVITY
	const int accel_vars[nDim] = { VAR_ACCEL, VAR_ACCEL+1, VAR_ACCEL+2 };
	const int color[num_children] = {
#if nDim == 3
		0, 1, 1, 0, 1, 0, 0, 1
#else
#error "Unknown nDim in color (smooth)"
#endif
	};
#endif

#ifdef HYDRO_TRACERS
	int itracer, itracer_next;
	int total_local_tracers;
	tracerid_t tmp_local_tracers;
	tracerid_t total_tracers;
#endif /* HYDRO_TRACERS */

	/* test root cells */
	cart_assert( num_cells_per_level[min_level] == proc_sfc_index[local_proc_id+1] - proc_sfc_index[local_proc_id] );

	cart_assert( proc_sfc_index[0] == 0 );
	cart_assert( proc_sfc_index[num_procs] == num_root_cells );

	MPI_Allreduce( &num_cells_per_level[min_level], &total_root_cells, 1, MPI_INT, MPI_SUM, mpi.comm.run );

	cart_assert( total_root_cells == num_root_cells );

	/* check level counts */
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			cart_assert( cell_is_local(icell) );
		}
		cart_free( level_cells );
		cart_assert( num_level_cells == num_cells_per_level[level] );

		select_level( level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			cart_assert( !cell_is_local(icell) );
		}
		cart_free( level_cells ); 
		cart_assert( num_level_cells == num_buffer_cells[level] );
	}

	/* check neighbors */
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#ifdef OPENMP_DECLARE_CONST
#pragma omp parallel for default(none), private(i,j,icell,sfc,neighbors,pos,ioct), shared(num_level_cells,level_cells,level,proc_sfc_index,local_proc_id,cell_child_oct,oct_neighbors,oct_level,oct_parent_cell,oct_parent_root_sfc,oct_pos,size_cell_array,size_oct_array), shared(reverse_direction)
#else  /* OPENMP_DECLARE_CONST */
#pragma omp parallel for default(none), private(i,j,icell,sfc,neighbors,pos,ioct), shared(num_level_cells,level_cells,level,proc_sfc_index,local_proc_id,cell_child_oct,oct_neighbors,oct_level,oct_parent_cell,oct_parent_root_sfc,oct_pos,size_cell_array,size_oct_array)
#endif /* OPENMP_DECLARE_CONST */
		for ( k = 0; k < num_level_cells; k++ ) {
			icell = level_cells[k];

			cart_assert( icell >= 0 && icell < num_cells );
			cart_assert( cell_level(icell) == level );
			cart_assert( level == min_level || cell_level( cell_parent_cell( icell ) ) == level - 1 );

			sfc = cell_parent_root_sfc(icell);
			cart_assert( sfc >= 0 && sfc < max_sfc_index );
			cart_assert( !cell_is_local(icell) || ( sfc >= proc_sfc_index[local_proc_id] && sfc < proc_sfc_index[local_proc_id+1] ) );

			cell_all_neighbors( icell, neighbors );
			cell_center_position( icell, pos );

			for ( i = 0; i < num_neighbors; i++ ) {
				if ( neighbors[i] == NULL_OCT && cell_is_local(icell) ) {
					cart_debug("level = %d", level );
					cart_debug("icell = %d", icell );
					cart_debug("cell_level = %d", cell_level(icell ) );
					cart_debug("sfc = %d", cell_parent_root_sfc(icell) );
					cart_debug("cell_is_local(icell) = %u", cell_is_local(icell) );
					cart_debug("neighbors[%u] = %d", i, neighbors[i] );
					cart_debug("oct = %d", cell_child_oct[icell] );
					if ( cell_is_refined(icell) ) {
						cart_debug("oct neighbor = %d", oct_neighbors[ cell_child_oct[icell] ][i] );
					}
				}
				cart_assert( neighbors[i] != NULL_OCT || !cell_is_local(icell) );

				if ( neighbors[i] != NULL_OCT ) {
					if ( cell_level(neighbors[i]) != level && cell_level(neighbors[i]) != level-1 ) {
						cart_debug("level = %u", level );
						cart_debug("icell = %u", icell );
						cart_debug("neighbor[%u] = %d, %u, %u", i, neighbors[i], cell_level(neighbors[i]), 
							cell_is_local(neighbors[i]) );

						cell_center_position(icell,pos);
						cart_debug("cell position = %le %le %le", pos[0], pos[1], pos[2] );
						cell_center_position(neighbors[i],pos);
						cart_debug("neighbors position = %le %le %le", pos[0], pos[1], pos[2] );

						cart_debug("parent oct = %u", cell_parent_oct(icell) );
						cart_debug("oct neighbor = %u", oct_neighbors[cell_parent_oct(icell)][i] );
					}

					cart_assert( cell_level( neighbors[i] ) == level || cell_level( neighbors[i] ) == level - 1 );

					if ( cell_level( neighbors[i] ) == level ) {
						cart_assert( cell_neighbor( neighbors[i], reverse_direction[i] ) == icell );
					}
				}
				cart_assert( neighbors[i] != NULL_OCT || !cell_is_local(icell) );
			}

			if ( cell_is_refined(icell) ) {
				ioct = cell_child_oct[icell];
				cart_assert( ioct >= 0 && ioct < num_octs );
				cart_assert( oct_level[ioct] == level+1 );
				cart_assert( oct_parent_cell[ioct] == icell );
				cart_assert( oct_parent_root_sfc[ioct] == sfc );

				for ( i = 0; i < num_neighbors; i++ ) {
					cart_assert( neighbors[i] == oct_neighbors[ioct][i] );
				}

				for ( i = 0; i < nDim; i++ ) {
					cart_assert( fabs(pos[i]-oct_pos[ioct][i])/fabs(pos[i]) < 1e-12 );
				}
			}
		}
		cart_free( level_cells );
	}

	for ( i = 0; i < num_vars; i++ ) {
		max_var[i] = -1e20;
		min_var[i] = 1e20;
	}

	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			for ( j = 0; j < num_vars; j++ ) {
				max_var[j] = MAX( max_var[j], cell_var(icell,j) );
				min_var[j] = MIN( min_var[j], cell_var(icell,j) );
			}
		}

		cart_free( level_cells );
	}

	for ( i = 0; i < num_vars; i++ ) {
		cart_debug("var %u: max = %e, min = %e", i, max_var[i], min_var[i] );
	}

#ifdef PARTICLES
	/* check consistency between particle species arrays */
	tmp = 0;
	cart_assert( particle_species_indices[0] == 0 );
	cart_assert( num_particle_species >= 0 && num_particle_species <= MAX_PARTICLE_SPECIES );
	for ( i = 0; i < num_particle_species; i++ ) {
		cart_assert( particle_species_num[i] >= 0 );
		if ( particle_species_indices[i+1] != tmp + particle_species_num[i] ) {
			for ( j = 0; j < num_particle_species; j++ ) {
				cart_debug("particle_species_num[%u] = %ld", j, particle_species_num[j]);
			}
			for ( j = 0; j < num_particle_species+1; j++ ) {	
				cart_debug("particle_species_indices[%u] = %ld", j, particle_species_indices[j] );
			}
			cart_debug("i = %d", i );
			cart_debug("tmp = %ld", tmp );
		}
		cart_assert( particle_species_indices[i+1] == tmp + particle_species_num[i] );
		tmp += particle_species_num[i];
	}
	cart_assert( tmp == num_particles_total );

	/* check particles */
	total_local_particles = 0;
	for ( i = 0; i < num_particle_species; i++ ) {
		species_count[i] = 0;
	}

	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL ) {
			if( particle_id[i] >= num_particles_total )
			  {
			    cart_error("Incorrect particle[%d] id=%ld, num_particles_total=%ld",i,particle_id[i],num_particles_total);
			  }
			total_local_particles++;
			species_count[ particle_species( particle_id[i] ) ]++;
		}
	}

	if ( total_local_particles != num_local_particles ) {
		cart_debug("total_local_particles = %u", total_local_particles );
		cart_debug("num_local_particles = %u", num_local_particles );
	}
	cart_assert( total_local_particles == num_local_particles );

	MPI_Allreduce( species_count, species_count_total, num_particle_species, MPI_PARTICLEID_T, MPI_SUM, mpi.comm.run );

	for ( i = 0; i < num_particle_species; i++ ) {
		cart_assert( species_count_total[i] == particle_species_num[i] );
	}

	total_local_particles = 0;
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			ipart = cell_particle_list[icell];
			cart_assert( ipart == NULL_PARTICLE || particle_list_prev[ipart] == NULL_PARTICLE );

			while ( ipart != NULL_PARTICLE ) {
				cart_assert( particle_level[ipart] != FREE_PARTICLE_LEVEL );
				cart_assert( particle_id[ipart] < num_particles_total );

				total_local_particles++;
				ipart_next = particle_list_next[ipart];
				cart_assert( ipart_next == NULL_PARTICLE || particle_list_prev[ipart_next] == ipart );
				ipart = ipart_next;
			}
		}
		cart_free( level_cells );
	}

	if ( total_local_particles != num_local_particles ) {
		cart_debug("num_local_particles = %u", num_local_particles );
		cart_debug("total_local_particles = %u", total_local_particles );
	}
	cart_assert( total_local_particles == num_local_particles );

	tmp = (particleid_t)total_local_particles;
	MPI_Allreduce( &tmp, &total_particles, 1, MPI_PARTICLEID_T, MPI_SUM, mpi.comm.run );

	cart_assert( total_particles == num_particles_total );

#ifdef STAR_FORMATION
	total_local_particles = 0;
	for ( i = 0; i < num_star_particles; i++ ) {
		if ( particle_is_star(i) ) {
			total_local_particles++;
		}
	}

	if ( total_local_particles != num_local_star_particles ) {
		cart_debug("total_local_particles= %u", total_local_particles );
		cart_debug("num_local_star_particles = %u", num_local_star_particles );
	}
	cart_assert( total_local_particles == num_local_star_particles );

	tmp = (particleid_t)total_local_particles;
	MPI_Allreduce( &tmp, &total_particles, 1, MPI_PARTICLEID_T, MPI_SUM, mpi.comm.run );

	cart_assert( total_particles == particle_species_num[num_particle_species-1] );

	for ( i = num_star_particles; i < num_particles; i++ ) {
		cart_assert( !particle_is_star(i) );
	}
#endif /* STAR_FORMATION */
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
	/* does the number of non-free tracers match num_local_tracers? */
	total_local_tracers = 0;
	for ( i = 0; i < num_tracers; i++ ) {
		if ( tracer_id[i] != NULL_TRACER ) {
			total_local_tracers++;
		}
	}

	cart_assert( total_local_tracers == num_local_tracers );

	/* does the number of tracers in cell lists match num_local_tracers? 
	 * also check the integrity of the cell tracer linked lists */
	total_local_tracers = 0;
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];
			itracer = cell_tracer_list[icell];
			cart_assert( itracer == NULL_TRACER || tracer_list_prev[itracer] == NULL_TRACER );

			while ( itracer != NULL_TRACER ) {
				cart_assert( tracer_id[itracer] != NULL_TRACER );
				cart_assert( tracer_id[itracer] < num_tracers_total );

				total_local_tracers++;
				itracer_next = tracer_list_next[itracer];
				if ( itracer_next != NULL_TRACER && tracer_list_prev[itracer_next] != itracer ) {
					cart_debug("itracer = %d", itracer );
					cart_debug("itracer_next = %d", itracer_next );
					cart_debug("tracer_list_prev[itracer_next] = %d", tracer_list_prev[itracer_next] );
					cart_debug("level = %d", level );
					cart_debug("icell = %d", icell );
					cart_debug("cell_tracer_list[icell] = %d", cell_tracer_list[icell] );
				}
				cart_assert( itracer_next == NULL_TRACER || tracer_list_prev[itracer_next] == itracer );
				itracer = itracer_next;
			}
		}
		cart_free( level_cells );
	}

	if ( total_local_tracers != num_local_tracers ) {
		cart_debug("total_local_tracers = %ld", total_local_tracers );
		cart_debug("num_local_tracers = %ld", num_local_tracers );
	}
	cart_assert( total_local_tracers == num_local_tracers );

	/* does the sum of num_local_tracers equal num_tracers_total? */
	tmp_local_tracers = num_local_tracers;
	MPI_Allreduce( &tmp_local_tracers, &total_tracers, 1, MPI_TRACERID_T, MPI_SUM, mpi.comm.run );

	cart_assert( total_tracers == num_tracers_total );
#endif /* HYDRO_TRACERS */

#ifdef GRAVITY
	for ( level = min_level+1; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );

#pragma omp parallel for private(icell,neighbors,i,j) 
		for ( k = 0; k < num_level_cells; k++ ) {
			icell = level_cells[k];
			cell_all_neighbors( icell, neighbors );

			for ( i = 0; i < num_neighbors; i++ ) {
				if ( neighbors[i] == NULL_OCT ) {
					break;
				} else if ( cell_is_local(neighbors[i]) ) {
					j = cell_child_number( icell );

					/* check that red direct neighbors have all neighbors */
					if ( color[j] ) {
						for ( j = i+1; j < num_neighbors; j++ ) {
							cart_assert( neighbors[i] != NULL_OCT );
						}
					}
					break;
				}
			}
		}
		cart_free( level_cells );
	}

	/* use cell_accel vars to test buffer cells */
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

#pragma omp parallel for default(none) private(icell,pos,j) shared(num_level_cells,level_cells,cell_vars)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			cell_center_position(icell, pos);

			for ( j = 0; j < nDim; j++ ) {
				cell_accel(icell,j) = pos[j];
			}
		}

		cart_free( level_cells );

		update_buffer_level( level, accel_vars, nDim );

		select_level( level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );

#ifdef OPENMP_DECLARE_CONST
#pragma omp parallel for default(none) private(icell,pos,j) shared(level,num_level_cells,level_cells,cell_vars,cell_size)
#else
#pragma omp parallel for default(none) private(icell,pos,j) shared(level,num_level_cells,level_cells,cell_vars)
#endif
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			cell_center_position( icell, pos );

			for ( j = 0; j < nDim; j++ ) {
				cart_assert( fabs( cell_accel(icell,j) - pos[j] ) < 0.5*cell_size[level] );
			}
		}

		cart_free( level_cells );
	}
#endif /* GRAVITY */

	cart_debug("done with check_map");
}

void print_cell_values(int level) {
	int i, j;
	int num_level_cells;
	int *level_cells;
	int icell;
	float max_var[num_vars];
	float min_var[num_vars];

    for ( i = 0; i < num_vars; i++ ) {
        max_var[i] = -1e20;
        min_var[i] = 1e20;
    }

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		for ( j = 0; j < num_vars; j++ ) {
			if ( isnan(cell_var(icell,j)) ) {
				cart_debug("NaN found in level %u, cell %u, var %u", level, icell, j );
			}

			max_var[j] = MAX( max_var[j], cell_var(icell,j) );
			min_var[j] = MIN( min_var[j], cell_var(icell,j) );
		}
	}

	cart_free( level_cells );

    for ( i = 0; i < num_vars; i++ ) {
        cart_debug("var %u: max = %e, min = %e", i, max_var[i], min_var[i] );
    }
}
