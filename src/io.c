#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "defs.h"
#include "particle.h"
#include "timestep.h"
#include "tree.h"
#include "io.h"
#include "units.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "sfc.h"
#include "iterators.h"
#include "particle.h"
#include "auxiliary.h"
#include "timing.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "starformation.h"
#include "load_balance.h"
#include "skiplist.h"
#include "index_hash.h"
#include "refinement_indicators.h"

#ifdef RADIATIVE_TRANSFER
#include "rt_io.h"
#endif

char output_directory[256];
char logfile_directory[256];
char jobname[256];
char requeue_command[256];
float outputs[MAX_OUTPUTS];
int num_outputs;
int current_output;
int last_restart_step;

int num_output_files = 1;

#ifdef OLDSTYLE_PARTICLE_FILE_SINGLE_PRECISION
typedef float particle_float;
#define MPI_PARTICLE_FLOAT	MPI_FLOAT
#else
typedef double particle_float;
#define MPI_PARTICLE_FLOAT	MPI_DOUBLE
#endif

void reorder( char *buffer, int size ) {
        int i;
        char tmp;
                                                                                
        for ( i = 0; i < (size/2); i++ ) {
                tmp = buffer[i];
                buffer[i] = buffer[size - i - 1];
                buffer[size - i - 1] = tmp;
        }
}

void write_restart( int gas_filename_flag, int particle_filename_flag, int tracer_filename_flag ) {
	FILE *restart;
	char filename[256];
	char filename1[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];
	char filename_gas[256];
	char filename_tracers[256];

	start_time( IO_TIMER );

#ifdef HYDRO
	cart_debug("Writing gas restart...");
	switch(gas_filename_flag)
	  {
	  case WRITE_SAVE:
	    {
	      sprintf( filename_gas, "%s/%s_a%06.4f.d", output_directory, jobname, aexp[min_level] );
	      break;
	    }
	  case WRITE_BACKUP:
	    {
	      sprintf( filename_gas, "%s/%s_2.d", output_directory, jobname );
	      break;
	    }
	  default:
	    {
	      sprintf( filename_gas, "%s/%s.d", output_directory, jobname );
	    }
	}

	start_time( GAS_IO_TIMER );
	write_grid_binary2( filename_gas );
	end_time( GAS_IO_TIMER );

#ifdef HYDRO_TRACERS
	cart_debug("Writing hydro tracer restart...");
	switch(tracer_filename_flag)
	  {
	  case WRITE_SAVE:
	    {
	      sprintf( filename_tracers, "%s/tracers_a%06.4f.dat", output_directory, aexp[min_level] );
	      break;
	    }
	  case WRITE_BACKUP:
	    {
	      sprintf( filename_tracers, "%s/tracers_2.dat", output_directory );
	      break;
	    }
	  default:
	    {
	      sprintf( filename_tracers, "%s/tracers.dat", output_directory );
	    }
	}

	start_time( PARTICLE_IO_TIMER );
	write_hydro_tracers( filename_tracers );
	end_time( PARTICLE_IO_TIMER );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */

#ifdef PARTICLES
	cart_debug("Writing particle restart...");
	switch(particle_filename_flag)
	  {
	  case WRITE_SAVE:
	    {
	      sprintf( filename1,"%s/PMcrda%06.4f.DAT", output_directory, aexp[min_level] );
	      sprintf( filename2, "%s/PMcrs0a%06.4f.DAT", output_directory, aexp[min_level] );
	      sprintf( filename3, "%s/pta%06.4f.dat", output_directory, aexp[min_level] );
	      sprintf( filename4, "%s/stars_a%06.4f.dat", output_directory, aexp[min_level] );
	      break;
	    }
	  case WRITE_BACKUP:
	    {
		sprintf( filename1, "%s/PMcrd_2.DAT", output_directory );
		sprintf( filename2, "%s/PMcrs_2.DAT", output_directory );
		sprintf( filename3, "%s/pt_2.dat", output_directory );
		sprintf( filename4, "%s/stars_2.dat", output_directory );
		break;
	    }
	  default:
	    {
		sprintf( filename1, "%s/PMcrd.DAT", output_directory );
		sprintf( filename2, "%s/PMcrs.DAT", output_directory );
		sprintf( filename3, "%s/pt.dat", output_directory );
		sprintf( filename4, "%s/stars.dat", output_directory );
	    }
	}
		
	start_time( PARTICLE_IO_TIMER );
#ifdef STARFORM
	write_particles( filename1, filename2, filename3, filename4 );
#else
	write_particles( filename1, filename2, filename3, NULL );
#endif /* STARFORM */
	end_time( PARTICLE_IO_TIMER );

#endif /* PARTICLES */

	if ( local_proc_id == MASTER_NODE ) {
		/* write out restart file */
		sprintf( filename, "%s/restart.dat", output_directory );
		restart = fopen( filename, "w" );

#ifdef HYDRO
		fprintf( restart, "%s\n", filename_gas );

#ifdef HYDRO_TRACERS
		fprintf( restart, "%s\n", filename_tracers );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */

#ifdef PARTICLES
		fprintf( restart, "%s\n", filename1 );
		fprintf( restart, "%s\n", filename2 );
		fprintf( restart, "%s\n", filename3 );
#ifdef STARFORM
		fprintf( restart, "%s\n", filename4 );
#endif /* STARFORM */
#endif /* PARTICLES */
	
		fclose(restart);
	}

	save_auxiliary();
	last_restart_step = step;

	end_time ( IO_TIMER );
}

void read_restart( double aexpn ) {
	int i,j;
	int icell, num_level_cells;
	int *level_cells;
	int level;
	FILE *restart;
	char filename[256];
	char filename_gas[256];
	char filename1[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];
	char filename_tracers[256];

	if ( buffer_enabled ) {
		destroy_cell_buffer();
	}

	if ( aexpn == 0.0 ) {
		/* read filenames from restart file */
		sprintf( filename, "%s/restart.dat", output_directory );
		restart = fopen( filename, "r" );

		if ( restart == NULL ) {
			cart_debug("Unable to locate restart.dat, trying default filenames!");

			/* try generic names */
			sprintf( filename_gas, "%s/%s.d", output_directory, jobname );
			sprintf( filename1,  "%s/PMcrd.DAT", output_directory );
			sprintf( filename2, "%s/PMcrs.DAT", output_directory );
			sprintf( filename3, "%s/pt.dat", output_directory );
			sprintf( filename4, "%s/stars.dat", output_directory );
			sprintf( filename_tracers, "%s/tracers.dat", output_directory );
		} else {
#ifdef HYDRO
			fscanf( restart, "%s\n", filename_gas );
#ifdef HYDRO_TRACERS
			fscanf( restart, "%s\n", filename_tracers );
#endif /* HYDRO_TRACERS */
#endif /* HYDRO */
			fscanf( restart, "%s\n", filename1 );
			fscanf( restart, "%s\n", filename2 );
			fscanf( restart, "%s\n", filename3 );
#ifdef STARFORM
			fscanf( restart, "%s\n", filename4 );
#endif /* STARFORM */
			fclose(restart);
		}
	} else {
		sprintf( filename_gas, "%s/%s_a%06.4f.d", output_directory, jobname, aexpn );
		sprintf( filename1,  "%s/PMcrda%06.4f.DAT", output_directory, aexpn );
		sprintf( filename2, "%s/PMcrs0a%06.4f.DAT", output_directory, aexpn );
		sprintf( filename3, "%s/pta%06.4f.dat", output_directory, aexpn );
		sprintf( filename4, "%s/stars_a%06.4f.dat", output_directory, aexpn );
		sprintf( filename_tracers, "%s/tracers_a%06.4f.dat", output_directory, aexpn );
	}

	/* do load balancing */
	restart_load_balance( filename_gas, filename1, filename2 );	

#ifdef HYDRO
	cart_debug("Reading gas restart...");
	read_grid_binary2( filename_gas );
	init_units();

#ifdef HYDRO_TRACERS
	read_hydro_tracers( filename_tracers );
#endif /* HYDRO_TRACERS */

#else
	build_root_cell_buffer();
#endif /* HYDRO */

#ifdef PARTICLES
	cart_debug("Reading particle restart...");
#ifdef STARFORM
	read_particles( filename1, filename2, filename3, filename4, 0, NULL );
#else
#ifdef HYDRO
	read_particles( filename1, filename2, filename3, NULL, 0, NULL );
#else
	cart_debug("reading in files %s and %s", filename1, filename2 );
	read_particles( filename1, filename2, NULL, NULL, 0, NULL );
#endif /* HYDRO */
#endif /* STARFORM */

#ifndef HYDRO
	init_units();
#endif /* HYDRO */
#endif /* PARTICLES */

	/* ensure current output points to correct place in output array */
	current_output = 0;
	while ( current_output < num_outputs &&
			aexp[min_level] >= outputs[current_output] ) {
		current_output++;
	}
}

void save_check() {
	char filename1[256];
	char filename2[256];
	char filename3[256];
	char filename4[256];

	int particle_save_flag = WRITE_GENERIC;
	int tracer_save_flag = WRITE_GENERIC;

	if ( particle_output_frequency != 0 && step % particle_output_frequency == 0 ) {
		particle_save_flag = WRITE_SAVE;
	}

	if ( tracer_output_frequency != 0 && step % tracer_output_frequency == 0 ) {
		tracer_save_flag = WRITE_SAVE;
	}

	if ( ( grid_output_frequency != 0 && step % grid_output_frequency == 0 ) ||
		( current_output < num_outputs && aexp[min_level] >= outputs[current_output] ) ) {
		/* we're saving information for permenant saving */
		write_restart( WRITE_SAVE, WRITE_SAVE, WRITE_SAVE );
		current_output++;
	} else if ( restart_frequency != 0 && step % restart_frequency == 0 ) {
		write_restart( WRITE_GENERIC, particle_save_flag, tracer_save_flag );
		/* A second restart file, in case the primary one is corrupted */
		if(step%(2*restart_frequency) == 0) write_restart(WRITE_BACKUP,WRITE_BACKUP,WRITE_BACKUP);
	} else {

#ifdef PARTICLES
		if ( particle_save_flag == WRITE_SAVE ) {
			start_time( IO_TIMER );
			start_time( PARTICLE_IO_TIMER );

			/* only write out particles, no restart */
			sprintf( filename1, "%s/PMcrda%06.4f.DAT", output_directory, aexp[min_level] );
			sprintf( filename2, "%s/PMcrs0a%06.4f.DAT", output_directory, aexp[min_level] );
			sprintf( filename3, "%s/pta%06.4f.dat", output_directory, aexp[min_level] );

#ifdef STARFORM
			sprintf( filename4, "%s/stars_a%06.4f.dat", output_directory, aexp[min_level] );
			write_particles( filename1, filename2, filename3, filename4 );
#else
			write_particles( filename1, filename2, filename3, NULL );
#endif /* STARFORM */

			end_time( PARTICLE_IO_TIMER );
			end_time( IO_TIMER );
		}
#endif /* PARTICLES */

#ifdef HYDRO_TRACERS
		if ( tracer_save_flag == WRITE_SAVE ) {
			start_time( IO_TIMER );
			start_time( PARTICLE_IO_TIMER );
			sprintf( filename1, "%s/tracers_a%06.4f.dat", output_directory, aexp[min_level] );
			write_hydro_tracers( filename1 );
			end_time( PARTICLE_IO_TIMER );
			end_time( IO_TIMER );
		}
#endif /* HYDRO_TRACERS */
	}
}

void restart_load_balance( char *grid_filename, char *particle_header_filename, char *particle_data ) {
	int i, j;
	int index;
	int coords[nDim];
	int page;
	int num_read;
	float *cell_work;	
	int *constrained_quantities;

#ifdef PARTICLES
	int num_parts_per_page;
	int num_parts_in_page;
	int num_pages;
	particle_float *x, *y, *z;
	particle_float *input_page;
	particle_header header;
#endif /* PARTICLES */

	FILE *input;
	int endian, nbody_flag;
	int grid_change_flag;
	int size, value;
	int *cellrefined;
	double rfact, vfact;
	double grid_shift;
	char filename[256];
	
	if ( num_procs == 1 ) {
		proc_sfc_index[0] = 0;
		proc_sfc_index[1] = num_root_cells;
		init_tree();
		return;
	}

	if ( local_proc_id == MASTER_NODE ) {
		/* do load balancing */
		constrained_quantities = cart_alloc( num_constraints*num_root_cells*sizeof(int) );
		cell_work = cart_alloc( num_root_cells * sizeof(float) );

		for ( i = 0; i < num_root_cells; i++ ) {
			cell_work[i] = 0.0;
		}

		for ( i = 0; i < num_constraints*num_root_cells; i++ ) {
			constrained_quantities[i] = 0;
		}

		/* load particle work information */
#ifdef PARTICLES
		if ( particle_header_filename != NULL ) {
			read_particle_header( particle_header_filename, &header, &endian, &nbody_flag );

			num_particles_total = header.num[ header.Nspecies - 1 ];
			num_parts_per_page = header.Nrow*header.Nrow;
			num_pages = (num_particles_total-1) / num_parts_per_page + 1;

			cart_debug("num_parts_per_page = %d", num_parts_per_page );
			cart_debug("num_particles_total = %d", num_particles_total );
			cart_debug("num_pages = %d", num_pages );
	
			if ( nbody_flag ) {
				grid_shift = 1.5;
			} else {
				grid_shift = 1.0;
			}

			if ( header.Ngrid != num_grid ) {
				rfact = (float)num_grid / (float)header.Ngrid;
				grid_change_flag = 1;
			} else {
				grid_change_flag = 0;
			}

			input_page = cart_alloc( nDim*num_parts_per_page*sizeof(particle_float) );

			x = input_page;
			y = &input_page[num_parts_per_page];
			z = &input_page[2*num_parts_per_page];

			input = fopen( particle_data, "r" );
			if ( input == NULL ) {
				cart_error( "Unable to load particle data file!\n");
			}

			for ( page = 0; page < num_pages; page++ ) {
				if ( page == num_pages - 1 ) {
					num_parts_in_page = num_particles_total - 
							num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				num_read = fread( input_page, sizeof(particle_float), 
							nDim*num_parts_per_page, input );

				if ( num_read != nDim*num_parts_per_page ) {
					cart_error("Error reading from particle file %s: insufficient data", particle_data );
				}

				if ( endian ) {
					for ( j = 0; j < num_parts_in_page; j++ ) {
						reorder( (char *)&x[j], sizeof(particle_float) );
						reorder( (char *)&y[j], sizeof(particle_float) );
						reorder( (char *)&z[j], sizeof(particle_float) );
					}
				}

				for ( j = 0; j < num_parts_in_page; j++ ) {
					x[j] -= grid_shift;
					y[j] -= grid_shift;
					z[j] -= grid_shift;

					if ( grid_change_flag ) {
						x[j] *= rfact;
						y[j] *= rfact;
						z[j] *= rfact;
					}

					/* enforce periodic boundary conditions */
					if ( x[j] < 0.0 ) {
						x[j] += (double)num_grid;
					} else if ( x[j] >= (double)num_grid ) {
						x[j] -= (double)num_grid;
					}

					if ( y[j] < 0.0 ) {
						y[j] += (double)num_grid;
					} else if ( y[j] >= (double)num_grid ) {
						y[j] -= (double)num_grid;
					}

					if ( z[j] < 0.0 ) {
						z[j] += (double)num_grid;
					} else if ( z[j] >= (double)num_grid ) {
						z[j] -= (double)num_grid;
					}

					coords[0] = (int)(x[j]);
					coords[1] = (int)(y[j]);
					coords[2] = (int)(z[j]);

					index = sfc_index( coords );
					cart_assert( index >= 0 && index < num_root_cells );

					constrained_quantities[num_constraints*index+1]++;
					cell_work[index] += cost_per_particle;
				}

				fseek( input, nDim*num_parts_per_page*sizeof(particle_float), SEEK_CUR );
			}

			cart_free( input_page );
			fclose(input);
		}
#endif /* PARTICLES */

		if ( grid_filename != NULL ) {
			index = 0;
			for ( i = 0; i < num_output_files; i++ ) {
				if ( num_output_files == 1 ) {
					sprintf( filename, "%s", grid_filename );
				} else {
					sprintf( filename, "%s.%03u", grid_filename, i );
				}

				input = fopen( filename, "r" );
				if ( input == NULL ) {
					cart_error( "Unable to open file %s for reading!", filename );
				}

				if ( i == 0 ) {
					/* skip over header information */
					fread( &size, sizeof(int), 1, input );
					endian = 0;
					if ( size != 256 ) {
						reorder( (char *)&size, sizeof(int) );
						if ( size != 256 ) {
							cart_error("Error: file %s is corrupted", filename );
						} else {
							endian = 1;
						}
					}

					fseek( input, 256*sizeof(char)+sizeof(int), SEEK_CUR );

					value = 0;
					while ( value != num_root_cells ) {
						fread( &size, sizeof(int), 1, input );

						if ( endian ) {
							reorder( (char *)&size, sizeof(int) );
						}

						if ( size == sizeof(int) ) {
							fread( &value, sizeof(int), 1, input );

							if ( endian ) {
								reorder( (char *)&value, sizeof(int) );
							}
						} else {
							fseek( input, size, SEEK_CUR );
						}

						fread( &size, sizeof(int), 1, input );
					}
				}

				/* read cellrefined */
				fread( &size, sizeof(int), 1, input );
			
				if ( endian ) {
					reorder( (char *)&size, sizeof(int) );
				}

				cellrefined = cart_alloc( size );
				size /= sizeof(int);
				fread( cellrefined, sizeof(int), size, input );

				fclose(input);

				if ( endian ) {
					for ( j = 0; j < size; j++ ) {
						reorder( (char *)&cellrefined, sizeof(int) );
					}
				}

				for ( j = 0; j < size; j++ ) {
					constrained_quantities[num_constraints*index] += cellrefined[j];
					cell_work[index] += cost_per_cell*cellrefined[j];
					index++;
				}
			}

			cart_free( cellrefined );
		} else {
			for ( index = 0; index < num_root_cells; index++ ) {
				constrained_quantities[num_constraints*index] = 1;
				cell_work[index] += cost_per_cell;
			}
		}	

		cart_debug("load balancing before i/o");
		load_balance_entire_volume( cell_work, constrained_quantities, proc_sfc_index );

		cart_free( cell_work );
		cart_free( constrained_quantities );
	}

        /* let all other processors know what their new workload is */
        MPI_Bcast( proc_sfc_index, num_procs+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
        init_tree();
}

#ifdef PARTICLES

int compare_particle_ids( const void *a, const void *b ) {
	return ( particle_id[*(int *)a] - particle_id[*(int *)b] );
}

void write_particles( char *header_filename, char *data_filename, char *timestep_filename, char *stellar_filename ) {
	int i, j;
	int num_parts;
	char desc[46];
	int processor_heap[MAX_PROCS];
	int *particle_order;
	int *page_ids[MAX_PROCS];
	particle_float *page[MAX_PROCS];
	particle_float *output_page;
	float *timestep_page[MAX_PROCS];
	int count[MAX_PROCS];
	int pos[MAX_PROCS];
	int pages[MAX_PROCS];
	int num_parts_in_page, num_parts_per_page;
	int num_parts_per_proc_page;
	int num_pages, index;
	int local_count;
	int current_id;
	int size;
	FILE *output;
	FILE *timestep_output;
	float *output_times;
	float *output_stars;
	particle_header header;
	MPI_Status status;
	int left, right, root, smallest;
	int leftproc, rightproc, rootproc, smallestproc;
	int proc;

#ifdef STARFORM
	FILE *stellar_output;
	int num_stars;
	int first_star;
	float *star_page[MAX_PROCS];
#endif /* STARFORM */

	num_parts_per_page = num_row*num_row;
	num_parts_per_proc_page = num_parts_per_page/num_procs;
	num_pages = (num_particles_total-1) / num_parts_per_page + 1;

	/* construct mapping of id -> particle index */
	particle_order = cart_alloc( num_local_particles * sizeof(int) );
        num_parts = 0;
        for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL ) {
			particle_order[num_parts++] = i;
		}
        }
	cart_assert( num_parts == num_local_particles );

	qsort( particle_order, num_parts, sizeof(int), compare_particle_ids );

	/* write file header */
	if ( local_proc_id == MASTER_NODE ) {
		output = fopen( header_filename, "w");
		if ( output == NULL ) {
			cart_error("Unable to open %s for writing", header_filename );
		}

		/* set up header and write it out */
		header.aexpn	= aexp[min_level];
		header.aexp0	= a_init;
		header.amplt	= 0.0;
		header.astep	= aexp[min_level] - aexp_old[min_level];
		header.istep	= step;
		header.partw	= 0.0;
		header.tintg	= tintg;
		header.ekin	= ekin;
		header.ekin1	= ekin1;
		header.ekin2	= 0.0;
		header.au0	= au0;
		header.aeu0	= aeu0;
		header.Nrow	= num_row;
		header.Ngrid	= num_grid;
		header.Nspecies	= num_particle_species;
		header.Nseed	= 0.0;
		header.Om0	= Omega0;
		header.Oml0	= OmegaL0;
		header.hubble	= hubble;
		header.Wp5	= 0.0;
		header.Ocurv	= 0.0;
		header.Omb0	= Omegab0;

		for ( i = 0; i < num_particle_species; i++ ) {
			header.mass[i] = particle_species_mass[i];
			header.num[i] = particle_species_indices[i+1];
		}

		/* write jobname to header desc (head of file) */
		for ( i = 0; i < 45; i++ ) {
			jobname[i] = " "; /* for fortran version */
		}
		snprintf( desc, 45, "%s", jobname );

		size = sizeof(particle_header)+45;
		fwrite( &size, sizeof(int), 1, output );
		fwrite( desc, sizeof(char), 45, output );
		fwrite( &header, sizeof(particle_header), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		fclose( output );

		/* now prepare to write data file */
		output = fopen( data_filename, "w" );
		if ( output == NULL ) {
			cart_error("Unable to open %s for writing.", data_filename );
		}

		if ( timestep_filename != NULL ) {
			timestep_output = fopen( timestep_filename, "w" );
			if ( timestep_output == NULL ) {
				cart_error("Unable to open %s for writing.", timestep_filename );
			}

			size = num_particles_total * sizeof(float);
			fwrite( &size, sizeof(int), 1, timestep_output );

			output_times = cart_alloc( num_parts_per_page*sizeof(float) );

			for ( i = 0; i < num_procs; i++ ) {
				timestep_page[i] = cart_alloc( num_parts_per_proc_page*sizeof(float) );
			}
		}

		/* allocate space to receive pages of particles from other procs */
		for ( i = 1; i < num_procs; i++ ) {
			page[i] = cart_alloc( 2*nDim*num_parts_per_proc_page*sizeof(particle_float) );
			page_ids[i] = cart_alloc( num_parts_per_proc_page * sizeof(int) );
		}

		/* allocate actual page which will be written */
		output_page = cart_alloc( 2*nDim*num_parts_per_page*sizeof(particle_float) );

		processor_heap[0] = MASTER_NODE;
		page_ids[0] = cart_alloc( sizeof(int) );
		pos[0] = 0;
		if ( num_local_particles > 0 ) {
			page_ids[0][0] = particle_id[particle_order[0]];
		} else {
			page_ids[0][0] = -1;
		}

		/* receive initial pages */
		for ( i = 1; i < num_procs; i++ ) {
			MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, i, 0, MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[i] );
			MPI_Recv( page[i], 2*nDim*num_parts_per_proc_page, MPI_PARTICLE_FLOAT, i, 0, MPI_COMM_WORLD, &status );
			
			if ( timestep_filename != NULL ) {
				MPI_Recv( timestep_page[i], num_parts_per_proc_page, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &status );
			}

			pos[i] = 0;
			processor_heap[i] = i;
			pages[i] = 1;

			/* fix for case where one processor has no particles */
			if ( count[i] == 0 ) {
				page_ids[i][0] = -1;
			}
		}

		/* build processor heap ordered by minimum particle id */
		for ( i = num_procs/2-1; i >= 0; i-- ) {
			root = i;
			left = 2*root+1;
			while ( left < num_procs ) {
				leftproc = processor_heap[left];
				rootproc = processor_heap[root];

				if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
						page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
					smallest = left;
				} else {
					smallest = root;
				}

				right = left+1;
				if ( right < num_procs && ( page_ids[processor_heap[smallest]][0] == -1 || 
						( page_ids[processor_heap[right]][0] != -1 && 
						  page_ids[processor_heap[right]][0] < 
						  page_ids[processor_heap[smallest]][0] ) ) ) {
					smallest = right;
				}

				if ( smallest != root ) {
					/* swap */
					processor_heap[root] = processor_heap[smallest];
					processor_heap[smallest] = rootproc;
					root = smallest;
					left = 2*root+1;
				} else {
					left = num_procs;
				}
			}
		}

		current_id = 0;
		local_count = 0;
		for ( i = 0; i < num_pages; i++ ) {
			/* construct page */
			num_parts = 0;

			if ( i == num_pages - 1 ) {
				num_parts_in_page = num_particles_total - num_parts_per_page*(num_pages-1);
			} else {
				num_parts_in_page = num_parts_per_page;
			}

			cart_assert( num_parts_in_page > 0 && num_parts_in_page <= num_parts_per_page );

			while ( num_parts < num_parts_in_page ) {
				/* pop minimum processor off the heap */
				proc = processor_heap[0];
				cart_assert( page_ids[proc][pos[proc]] == current_id );

				if ( proc == local_proc_id ) {
					/* add from our list */
					while ( num_parts < num_parts_in_page &&
							local_count < num_local_particles &&
							particle_id[particle_order[local_count]] == current_id ) {
						index = particle_order[local_count];

						output_page[num_parts] = particle_x[index][0] + 1.0;
						output_page[num_parts_per_page+num_parts] = particle_x[index][1] + 1.0;
						output_page[2*num_parts_per_page+num_parts] = particle_x[index][2] + 1.0;
						output_page[3*num_parts_per_page+num_parts] = particle_v[index][0];
						output_page[4*num_parts_per_page+num_parts] = particle_v[index][1];
						output_page[5*num_parts_per_page+num_parts] = particle_v[index][2];

						if ( timestep_filename != NULL ) {
							output_times[num_parts] = particle_dt[index];
						}

						local_count++;
						num_parts++;
						current_id++;
					}

					if ( local_count < num_local_particles ) {
						page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
					} else {
						page_ids[local_proc_id][0] = -1;
					}
				} else {
					while ( num_parts < num_parts_in_page &&
							pos[proc] < count[proc] &&
							page_ids[proc][pos[proc]] == current_id ) {

						output_page[num_parts] = 			page[proc][2*nDim*pos[proc]]   + 1.0;
						output_page[num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+1] + 1.0;
						output_page[2*num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+2] + 1.0;
						output_page[3*num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+3];
						output_page[4*num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+4];
						output_page[5*num_parts_per_page+num_parts] = 	page[proc][2*nDim*pos[proc]+5];

						if ( timestep_filename != NULL ) {
							output_times[num_parts] = timestep_page[proc][pos[proc]];
						}

						num_parts++;
						current_id++;
						pos[proc]++;

						/* if we've run out, refill page */
						if ( pos[proc] == count[proc] ) { 
							if ( count[proc] == num_parts_per_proc_page ) {
								MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT, 
									proc, pages[proc], MPI_COMM_WORLD, &status );
								MPI_Get_count( &status, MPI_INT, &count[proc] );
								MPI_Recv( page[proc], 2*nDim*num_parts_per_proc_page, 
									MPI_PARTICLE_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
								MPI_Recv( timestep_page[proc], num_parts_per_proc_page,
									MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );

								pos[proc] = 0;
								pages[proc]++;

								if ( count[proc] == 0 ) {
									page_ids[proc][0] = -1;
								}
							} else {
								pos[proc] = 0;
								count[proc] = 0;
								page_ids[proc][0] = -1;
							}
						}
					}
				}

				/* add back to heap (heapify) */
				root = 0;
				left = 1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][pos[rootproc]] == -1 ||
							( page_ids[leftproc][pos[leftproc]] != -1 && 
							  page_ids[leftproc][pos[leftproc]] < 
							  page_ids[rootproc][pos[rootproc]] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs ) {
						rightproc = processor_heap[right];
						smallestproc = processor_heap[smallest];

						if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
								( page_ids[rightproc][pos[rightproc]] != -1 &&
								  page_ids[rightproc][pos[rightproc]] < 
								  page_ids[smallestproc][pos[smallestproc]] ) ) {
							smallest = right;
						}
					}

					if ( smallest != root ) {
                                                /* swap */
                                                processor_heap[root] = processor_heap[smallest];
                                                processor_heap[smallest] = rootproc;
                                                root = smallest;
                                                left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			cart_assert( num_parts == num_parts_in_page );

			/* write page */
			fwrite( output_page, sizeof(particle_float), 2*nDim*num_parts_per_page, output );

			/* write particle timesteps */
			if ( timestep_filename != NULL ) {
				fwrite( output_times, sizeof(float), num_parts_per_page, timestep_output );
			}
		}

		fclose( output );
		cart_free( output_page );

		if ( timestep_filename != NULL ) {
			size = num_particles_total * sizeof(float);
                        fwrite( &size, sizeof(int), 1, timestep_output );

			fclose( timestep_output );
			cart_free( output_times );

			for ( i = 0; i < num_procs; i++ ) {
				cart_free( timestep_page[i] );
			}
		}

		cart_assert( local_count == num_local_particles );

		for ( i = 1; i < num_procs; i++ ) {
			cart_assert( pos[i] == 0 || pos[i] == count[i] );
			cart_free( page[i] );
		}

#ifdef STARFORM
		if ( stellar_filename != NULL ) {
			stellar_output = fopen(stellar_filename, "w");

			if ( stellar_output == NULL ) {
				cart_error("ERROR: unable to open %s for writing.", stellar_filename );
			}

			/* allocate buffers */
			for ( i = 1; i < num_procs; i++ ) {
				star_page[i] = cart_alloc( num_parts_per_proc_page * sizeof(float) );
			}

			output_stars = cart_alloc( num_parts_per_page * sizeof(float) );

			num_stars = particle_species_num[num_particle_species-1];

			if ( num_stars == 0 ) {
				num_pages = 0;
				first_star = 0;
			} else {
				num_pages = (num_stars-1) / num_parts_per_page + 1;

				first_star = 0;
				while ( first_star < num_local_particles && 
						!particle_is_star( particle_order[first_star] ) ) {
					first_star++;
				}
			}

			/* write header */
			size = 2*sizeof(double);
			fwrite( &size, sizeof(int), 1, stellar_output );
			fwrite( &tl[min_level], sizeof(double), 1, stellar_output );
			fwrite( &aexp[min_level], sizeof(double), 1, stellar_output );
			fwrite( &size, sizeof(int), 1, stellar_output );

			size = sizeof(int);
			fwrite( &size, sizeof(int), 1, stellar_output );
			fwrite( &num_stars, sizeof(int), 1, stellar_output );
			fwrite( &size, sizeof(int), 1, stellar_output );

			size = 2*sizeof(double);
			fwrite( &size, sizeof(int), 1, stellar_output );
			fwrite( &total_stellar_mass, sizeof(double), 1, stellar_output );
			fwrite( &total_stellar_initial_mass, sizeof(double), 1, stellar_output );
			fwrite( &size, sizeof(int), 1, stellar_output );

			/* particle_mass (pw) */
			size = num_stars * sizeof(float);
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
	                        page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				pages[i] = 1;
				processor_heap[i] = i;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || 
							( page_ids[leftproc][0] != -1 && 
							  page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = particle_mass[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									count[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < 
								  page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );

			/* star_initial_mass (pw0) */
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
	                        page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				processor_heap[i] = i;
				pages[i] = 1;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
								page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = star_initial_mass[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );

			/* star_tbirth (tbirth) */
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
	                        page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				pages[i] = 1;
				processor_heap[i] = i;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
								page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = star_tbirth[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );

#ifdef ENRICH
			/* star_metallicity_II (zstII) */
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
	                        page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				pages[i] = 1;
				processor_heap[i] = i;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
								page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = star_metallicity_II[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );
#endif /* ENRICH */
#ifdef ENRICH_SNIa
			/* star_metallicity_Ia (zstIa) */
			fwrite( &size, sizeof(int), 1, stellar_output );

			processor_heap[0] = MASTER_NODE;
			pos[0] = 0;
			if ( first_star == num_local_particles ) {
				page_ids[0][0] = -1;
			} else {
	                        page_ids[0][0] = particle_id[particle_order[first_star]];
			}

			/* receive initial pages */
			for ( i = 1; i < num_procs; i++ ) {
				MPI_Recv( page_ids[i], num_parts_per_proc_page, MPI_INT, 
						i, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_INT, &count[i] );
				MPI_Recv( star_page[i], num_parts_per_proc_page, MPI_FLOAT, 
						i, 0, MPI_COMM_WORLD, &status );
				pos[i] = 0;
				pages[i] = 1;
				processor_heap[i] = i;

				if ( count[i] == 0 ) {
					page_ids[i][0] = -1;
				}
			}

			/* build processor heap ordered by minimum particle id */
			for ( i = num_procs/2-1; i >= 0; i-- ) {
				root = i;
				left = 2*root+1;
				while ( left < num_procs ) {
					leftproc = processor_heap[left];
					rootproc = processor_heap[root];

					if ( page_ids[rootproc][0] == -1 || ( page_ids[leftproc][0] != -1 && 
								page_ids[leftproc][0] < page_ids[rootproc][0] ) ) {
						smallest = left;
					} else {
						smallest = root;
					}

					right = left+1;
					if ( right < num_procs && 
							( page_ids[processor_heap[smallest]][0] == -1 || 
							( page_ids[processor_heap[right]][0] != -1 && 
							  page_ids[processor_heap[right]][0] < 
							  page_ids[processor_heap[smallest]][0] ) ) ) {
						smallest = right;
					}

					if ( smallest != root ) {
						/* swap */
						processor_heap[root] = processor_heap[smallest];
						processor_heap[smallest] = rootproc;
						root = smallest;
						left = 2*root+1;
					} else {
						left = num_procs;
					}
				}
			}

			local_count = first_star;
			current_id = particle_species_indices[num_particle_species-1];
			for ( i = 0; i < num_pages; i++ ) {
				/* construct page */
				num_parts = 0;

				if ( i == num_pages - 1 ) {
					num_parts_in_page = num_stars - num_parts_per_page*(num_pages-1);
				} else {
					num_parts_in_page = num_parts_per_page;
				}

				while ( num_parts < num_parts_in_page ) {
					/* choose proc from min heap */
					proc = processor_heap[0];
					cart_assert( page_ids[proc][pos[proc]] == current_id );

					if ( proc == local_proc_id ) {
						/* add from our list */
						while ( num_parts < num_parts_in_page &&
								local_count < num_local_particles &&
								particle_id[particle_order[local_count]] == current_id ) {
							index = particle_order[local_count];
							output_stars[num_parts] = star_metallicity_Ia[index];
							local_count++;
							num_parts++;
							current_id++;
						}

						if ( local_count < num_local_particles ) {
							page_ids[local_proc_id][0] = particle_id[particle_order[local_count]];
						} else {
							page_ids[local_proc_id][0] = -1;
						}
					} else {
						while ( num_parts < num_parts_in_page &&
								pos[proc] < count[proc] &&
								page_ids[proc][pos[proc]] == current_id ) {

							output_stars[num_parts] = star_page[proc][pos[proc]];
							num_parts++;
							current_id++;
							pos[proc]++;

							/* if we've run out, refill page */
							if ( pos[proc] == count[proc] ) {
								if ( count[proc] == num_parts_per_proc_page ) {
									MPI_Recv( page_ids[proc], num_parts_per_proc_page, MPI_INT,
											proc, pages[proc], MPI_COMM_WORLD, &status );
									MPI_Get_count( &status, MPI_INT, &count[proc] );
									MPI_Recv( star_page[proc], num_parts_per_proc_page,
											MPI_FLOAT, proc, pages[proc], MPI_COMM_WORLD, &status );
									pos[proc] = 0;
									pages[proc]++;

									if ( count[proc] == 0 ) {
										page_ids[proc][0] = -1;
									}
								} else {
									pos[proc] = 0;
									page_ids[proc][0] = -1;
								}
							}
						}
					}

					/* re-heapify */
					root = 0;
					left = 1;
					while ( left < num_procs ) {
						leftproc = processor_heap[left];
						rootproc = processor_heap[root];

						if ( page_ids[rootproc][pos[rootproc]] == -1 ||
								( page_ids[leftproc][pos[leftproc]] != -1 && 
								  page_ids[leftproc][pos[leftproc]] < page_ids[rootproc][pos[rootproc]] ) ) {
							smallest = left;
						} else {
							smallest = root;
						}

						right = left+1;
						if ( right < num_procs ) {
							rightproc = processor_heap[right];
							smallestproc = processor_heap[smallest];

							if ( page_ids[smallestproc][pos[smallestproc]] == -1 ||
									( page_ids[rightproc][pos[rightproc]] != -1 &&
									  page_ids[rightproc][pos[rightproc]] < 
									  page_ids[smallestproc][pos[smallestproc]] ) ) {
								smallest = right;
							}
						}

						if ( smallest != root ) {
							/* swap */
							processor_heap[root] = processor_heap[smallest];
							processor_heap[smallest] = rootproc;
							root = smallest;
							left = 2*root+1;
						} else {
							left = num_procs;
						}
					}
				}

				/* write out page */
				cart_assert( num_parts == num_parts_in_page );
				fwrite( output_stars, sizeof(float), num_parts_in_page, stellar_output );
			}

			fwrite( &size, sizeof(int), 1, stellar_output );

#endif /* ENRICH_SNIa */

			cart_free( output_stars );

			for ( i = 1; i < num_procs; i++ ) {
				cart_free( star_page[i] );
			}

			fclose( stellar_output );
		}
#endif /* STARFORM */

		for ( i = 0; i < num_procs; i++ ) {
			cart_free( page_ids[i] );
		}
	} else {
		page[local_proc_id] = cart_alloc( 2*nDim*num_parts_per_proc_page*sizeof(particle_float) );
		page_ids[local_proc_id] = cart_alloc( num_parts_per_proc_page*sizeof(int) );
		pages[local_proc_id] = 0;

		if ( timestep_filename != NULL ) {
			timestep_page[local_proc_id] = cart_alloc( num_parts_per_page*sizeof(float) );
		}	

		num_parts = 0;
		do {
			/* pack page */
			local_count = 0;
			for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
				page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
				page[local_proc_id][local_count++] = particle_x[particle_order[num_parts]][0];
				page[local_proc_id][local_count++] = particle_x[particle_order[num_parts]][1];
				page[local_proc_id][local_count++] = particle_x[particle_order[num_parts]][2];
				page[local_proc_id][local_count++] = particle_v[particle_order[num_parts]][0];
				page[local_proc_id][local_count++] = particle_v[particle_order[num_parts]][1];
				page[local_proc_id][local_count++] = particle_v[particle_order[num_parts]][2];

				if ( timestep_filename != NULL ) {
					timestep_page[local_proc_id][i] = particle_dt[particle_order[num_parts]];
				}

				num_parts++;
			}

			/* send the page */
			MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
			MPI_Send( page[local_proc_id], local_count, MPI_PARTICLE_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );

			if ( timestep_filename != NULL ) {
				MPI_Send( timestep_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
			}

			pages[local_proc_id]++;
		} while ( i == num_parts_per_proc_page );

		cart_free( page[local_proc_id] );

		if ( timestep_filename != NULL ) {
			cart_free( timestep_page[local_proc_id] );
		}

#ifdef STARFORM
		if ( stellar_filename != NULL ) {
			first_star = 0;
			while ( first_star < num_local_particles && !particle_is_star( particle_order[first_star] ) ) {
				first_star++;
			}

			star_page[local_proc_id] = cart_alloc( num_parts_per_page * sizeof(float) );

			/* send stellar mass (pw) */
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = particle_mass[particle_order[num_parts]];
					num_parts++;
				}
	
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );

			/* send star_initial_mass (pw0) */
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = star_initial_mass[particle_order[num_parts]];
					num_parts++;
				}
	
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );
	
			/* send star_tbirth */
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = star_tbirth[particle_order[num_parts]];
					num_parts++;
				}
				
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );

#ifdef ENRICH
			/* send metallicities */
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = star_metallicity_II[particle_order[num_parts]];
					num_parts++;
				}
	
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );
	
#ifdef ENRICH_SNIa
			num_parts = first_star;
			pages[local_proc_id] = 0;
			do {
				for ( i = 0; i < num_parts_per_proc_page && num_parts < num_local_particles; i++ ) {
					page_ids[local_proc_id][i] = particle_id[particle_order[num_parts]];
					star_page[local_proc_id][i] = star_metallicity_Ia[particle_order[num_parts]];
					num_parts++;
				}
	
				/* send the page */
				MPI_Send( page_ids[local_proc_id], i, MPI_INT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				MPI_Send( star_page[local_proc_id], i, MPI_FLOAT, MASTER_NODE, pages[local_proc_id], MPI_COMM_WORLD );
				pages[local_proc_id]++;
			} while ( i == num_parts_per_proc_page );
#endif /* ENRICH_SNIa */
#endif /* ENRICH */

			cart_free( star_page[local_proc_id] );
		}
#endif /* STARFORM */

		cart_free( page_ids[local_proc_id] );
	}

	cart_free( particle_order );
}

void read_particle_header( char *header_filename, particle_header *header, int *endian, int *nbody_flag ) {
	int i;
	FILE *input;
	nbody_particle_header nbody_header;
	char desc[46];
	int size;

	*nbody_flag = 0;
	*endian = 0;

	/* read file header */
	input = fopen( header_filename, "r");
	if ( input == NULL ) {
		cart_error("Unable to open particle file %s", header_filename );
	}

	fread( &size, sizeof(int), 1, input );
	fread( desc, sizeof(char), 45, input );
	desc[45] = '\0';

	cart_debug( "Particle header file: %s", desc );

	if ( size != sizeof(particle_header)+45 ) {
		if ( size == sizeof(nbody_particle_header)+45 ) {
			*nbody_flag = 1;
		} else {
			reorder( (char *)&size, sizeof(int) );

			if ( size != sizeof(particle_header)+45 ) {
				if ( size == sizeof(nbody_particle_header)+45 ) {
					*endian = 1;
					*nbody_flag = 1;
				} else {
					cart_error("Size mismatch in reading particle file header %s (%u vs %u)",
							header_filename, size, sizeof(particle_header)+45 );
				}
			} else {
				*endian = 1;
			}
		}
	}

	if ( *nbody_flag ) {
		cart_debug("USING OLD NBODY FILE FORMAT (hope that's what you wanted...)");

		fread( &nbody_header, sizeof(nbody_particle_header), 1, input );

		header->aexpn = nbody_header.aexpn;
		header->aexp0 = nbody_header.aexp0;
		header->amplt = nbody_header.amplt;
		header->astep = nbody_header.astep;
		header->istep = nbody_header.istep;
		header->partw = nbody_header.partw;
		header->tintg = nbody_header.tintg;
		header->ekin = nbody_header.ekin;
		header->ekin1 = nbody_header.ekin1;
		header->ekin2 = nbody_header.ekin2;
		header->au0 = nbody_header.au0;
		header->aeu0 = nbody_header.aeu0;
		header->Nrow = nbody_header.Nrow;
		header->Ngrid = nbody_header.Ngrid;
		header->Nspecies = nbody_header.Nspecies;
		header->Nseed = nbody_header.Nseed;
		header->Om0 = nbody_header.Om0;
		header->Oml0 = nbody_header.Oml0;
		header->hubble = nbody_header.hubble;
		header->Wp5 = nbody_header.Wp5;
		header->Ocurv = nbody_header.Ocurv;

		for ( i = 0; i < 10; i++ ) {
			header->mass[i] = nbody_header.mass[i];
			header->num[i] = nbody_header.num[i];
		}

	} else {
		fread( header, sizeof(particle_header), 1, input );
	}

	fread( &size, sizeof(int), 1, input );
	fclose(input);

	if ( *endian ) {
		reorder( (char *)&header->aexpn, sizeof(float) );
		reorder( (char *)&header->aexp0, sizeof(float) );
		reorder( (char *)&header->amplt, sizeof(float) );
		reorder( (char *)&header->astep, sizeof(float) );
		reorder( (char *)&header->istep, sizeof(int) );
		reorder( (char *)&header->partw, sizeof(float) );
		reorder( (char *)&header->tintg, sizeof(float) );
		reorder( (char *)&header->ekin, sizeof(float) );
		reorder( (char *)&header->ekin1, sizeof(float) );
		reorder( (char *)&header->ekin2, sizeof(float) );
		reorder( (char *)&header->au0, sizeof(float) );
		reorder( (char *)&header->aeu0, sizeof(float) );
		reorder( (char *)&header->Nrow, sizeof(int) );
		reorder( (char *)&header->Ngrid, sizeof(int) );
		reorder( (char *)&header->Nspecies, sizeof(int) );
		reorder( (char *)&header->Nseed, sizeof(int) );
		reorder( (char *)&header->Om0, sizeof(float) );
		reorder( (char *)&header->Oml0, sizeof(float) );
		reorder( (char *)&header->hubble, sizeof(float) );
		reorder( (char *)&header->Wp5, sizeof(float) );
		reorder( (char *)&header->Ocurv, sizeof(float) );
		reorder( (char *)&header->Omb0, sizeof(float) );

		for ( i = 0; i < header->Nspecies; i++ ) {
			reorder( (char *)&header->mass[i], sizeof(float) );
			reorder( (char *)&header->num[i], sizeof(int) );
		}
	}
}

void read_particles( char *header_filename, char *data_filename, 
			char *timestep_filename, char *stellar_filename,
			int num_sfcs, int *sfc_list ) {
	int i, j, k;
	char desc[47];
	int proc;
	int num_parts;
	int *page_ids[MAX_PROCS];
	particle_float *page[MAX_PROCS];
	particle_float *input_page, *x, *y, *z, *vx, *vy, *vz;
	float *timestep_page[MAX_PROCS];
	float *pdt;
	double dt;
	int ipart;
	int sfc;
	int num_read;
	int count[MAX_PROCS];
	int coords[nDim];
	int num_parts_in_page, num_parts_per_page;
	int num_parts_per_proc_page;
	int num_pages, index;
	int current_page[MAX_PROCS];
	int current_id, current_type, local_count;
	int size, endian, dt_endian, stellar_endian;
	int *sfc_map;
	FILE *input;
	FILE *timestep_input;
	particle_header header;
	int nbody_flag;
	int grid_change_flag;
	float rfact, vfact;
	float grid_shift;
	MPI_Status status;

#ifdef STARFORM
	FILE *stellar_input;
        int num_stars;
        double st, sa;
        int first_star, first_star_index, num_stars_to_read;
	long seek_amount, var_first;
        float *pw, *pw0, *tbirth, *zstII, *zstIa;
	float *star_page[MAX_PROCS];
	int num_star_vars[MAX_PROCS];
#ifdef ENRICH
#ifdef ENRICH_SNIa
	#define num_star_variables	5
#else
	#define num_star_variables	4
#endif /* ENRICH_SNIa */
#else
	#define num_star_variables	3
#endif /* ENRICH */
#endif /* STARFORM */

	if ( local_proc_id == MASTER_NODE ) {
		read_particle_header( header_filename, &header, &endian, &nbody_flag );

		cart_debug("aexpn = %e", header.aexpn );
		cart_debug("aexp0 = %e", header.aexp0 );
		cart_debug("amplt = %e", header.amplt );
		cart_debug("astep = %e", header.astep );
		cart_debug("istep = %u", header.istep );
		cart_debug("Ngrid = %u", header.Ngrid );
		cart_debug("Nrow = %u", header.Nrow );
		cart_debug("Nspecies = %u", header.Nspecies );
		cart_debug("Omega0 = %e", header.Om0 );
		cart_debug("OmegaL0 = %e", header.Oml0 );
		cart_debug("hubble = %e", header.hubble );
		cart_debug("Omegab0 = %e", header.Omb0 );

		/* set units */
		Omega0          = header.Om0;
		OmegaL0         = header.Oml0;
		hubble          = header.hubble;

		aexp[min_level]	= header.aexpn;
		aexp_old[min_level] = aexp[min_level] - header.astep;
		tl[min_level]  = a2b( aexp[min_level] );

		if ( header.astep > 0.0 ) {
			dt = tl[min_level] - a2b( aexp[min_level] - header.astep );
		} else {
			dt = 0.0;
		}

		a_init		= header.aexp0;
		step		= header.istep;
	
		num_row		= header.Nrow;
	
		if ( nbody_flag ) {
			vfact = 2.0/sqrt(Omega0);
			grid_shift = 1.5;
		} else {
			Omegab0 = header.Omb0;
			vfact = 1.0;
			grid_shift = 1.0;
		}
	
		/* only root node needs to keep energy conservation variables */
		tintg		= header.tintg;
		ekin		= header.ekin;
		ekin1		= header.ekin1;
		ekin2		= header.ekin2;
		au0		= header.au0;
		aeu0		= header.aeu0;

		ap0 		= b2a( tl[min_level] - 0.5*dt );

#ifndef OLDSTYLE_PARTICLE_FILE_IGNORE_NGRID
		if ( header.Ngrid != num_grid ) {
			cart_debug( "Mismatch between particle file num_grid and compiled value!" );

			rfact = (float)num_grid / (float)header.Ngrid;
			vfact *= rfact;
			grid_change_flag = 1;

			for ( i = 0; i < header.Nspecies; i++ ) {
				header.mass[i] *= rfact*rfact*rfact;
			}
		} else
#endif
		{
			grid_change_flag = 0;
		}

		if ( header.Nspecies > MAX_PARTICLE_SPECIES ) {
			cart_error( "header.Nspecies > MAX_PARTICLE_SPECIES.  Increase and rerun.");
		}

		num_particle_species = header.Nspecies;
	
		particle_species_indices[0] = 0;
		for ( i = 0; i < num_particle_species; i++ ) {
			particle_species_mass[i] = header.mass[i];
			particle_species_indices[i+1] = header.num[i];
			particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];

			cart_debug("particle_species_mass[%u] = %e", i, particle_species_mass[i] );
			cart_debug("particle_species_num[%u] = %u", i, particle_species_num[i] );
		}

		for ( i = 0; i <= num_particle_species; i++ ) {
			cart_debug("particle_species_indices[%u] = %d", i, particle_species_indices[i] );
		}

		num_particles_total = particle_species_indices[num_particle_species];
		cart_debug("num_particles_total = %u", num_particles_total );

                num_parts_per_page = header.Nrow*header.Nrow;
                num_parts_per_proc_page = num_parts_per_page/num_procs;
		num_pages = (num_particles_total-1) / num_parts_per_page + 1;

		cart_assert( num_pages > 0 && num_parts_per_page*num_pages >= num_particles_total );

#ifdef STARFORM
		if ( stellar_filename == NULL ) {
			/* ICs don't include stars, initialize species to 0 */
			if ( num_particle_species+1 > MAX_PARTICLE_SPECIES ) {
				cart_error("header.Nspecies > MAX_PARTICLE_SPECIES.  Increase and rerun.");
			}
			particle_species_indices[num_particle_species+1] = particle_species_indices[num_particle_species];
			particle_species_num[num_particle_species] = 0;
			particle_species_mass[num_particle_species] = 0.0;
			num_particle_species++;

			total_stellar_mass = 0.0;
			total_stellar_initial_mass = 0.0;
		} else {
			stellar_input = fopen( stellar_filename, "r" );
			if ( stellar_input == NULL ) {
				cart_error("Unable to open file %s for reading.", stellar_filename );
			}

			/* read in header */
			fread( &size, sizeof(int), 1, stellar_input );
			if ( size != 2*sizeof(double) ) {
				reorder( (char *)&size, sizeof(int) );

				if ( size != 2*sizeof(double) ) {
					cart_error("Error reading from %s.", stellar_filename );
				}

				stellar_endian = 1;
			} else {
				stellar_endian = 0;
			}

			fread( &st, sizeof(double), 1, stellar_input );
			fread( &sa, sizeof(double), 1, stellar_input );
			fread( &size, sizeof(int), 1, stellar_input );

			fread( &size, sizeof(int), 1, stellar_input );
			fread( &num_stars, sizeof(int), 1, stellar_input );
			fread( &size, sizeof(int), 1, stellar_input );

			if ( stellar_endian ) {
				reorder( (char *)&num_stars, sizeof(int) );
			}

			if ( num_stars != particle_species_num[num_particle_species-1] ) {
				cart_error("num_stars in %s doesn't match last particle specie.", stellar_filename );
			}

			fread( &size, sizeof(int), 1, stellar_input );
			fread( &total_stellar_mass, sizeof(double), 1, stellar_input );
			fread( &total_stellar_initial_mass, sizeof(double), 1, stellar_input );
			fread( &size, sizeof(int), 1, stellar_input );

			if ( stellar_endian ) {
				reorder( (char *)&total_stellar_mass, sizeof(double) );
				reorder( (char *)&total_stellar_initial_mass, sizeof(double) );
			}

			pw = cart_alloc( num_parts_per_page * sizeof(float) );
			pw0 = cart_alloc( num_parts_per_page * sizeof(float) );
			tbirth = cart_alloc( num_parts_per_page * sizeof(float) );

#ifdef ENRICH
			zstII = cart_alloc( num_parts_per_page * sizeof(float) );
#ifdef ENRICH_SNIa
			zstIa = cart_alloc( num_parts_per_page * sizeof(float) );
#endif /* ENRICH_SNIa */
#endif /* ENRICH */

			/* allocate buffer space for particles on other processors */
			for ( i = 1; i < num_procs; i++ ) {
				star_page[i] = cart_alloc( num_star_variables*num_parts_per_proc_page*sizeof(float) );
				num_star_vars[i] = 0;
			}
		}
#endif /* STARFORM */

		MPI_Bcast( &aexp[min_level], 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &tl[min_level], 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &dt, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &a_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &Omega0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &OmegaL0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &Omegab0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &hubble, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &step, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &num_row, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &num_particle_species, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &num_particles_total, 1, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &num_parts_per_page, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &num_parts_per_proc_page, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &num_pages, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( particle_species_mass, num_particle_species, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( particle_species_num, num_particle_species, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( particle_species_indices, num_particle_species+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

		if ( num_sfcs > 0 ) {
			/* create hash for quick checking if we should read in particles */
			sfc_map = cart_alloc( num_root_cells * sizeof(int) );

			for ( i = 0; i < num_root_cells; i++ ) {
				sfc_map[i] = 0;
			}

			for ( sfc = 0; sfc < num_sfcs; sfc++ ) {
				sfc_map[ sfc_list[sfc] ] = 1;
			}
		}

		input_page = cart_alloc( 2*nDim*num_parts_per_page*sizeof(particle_float) );

		x = input_page;
		y = &input_page[num_parts_per_page];
		z = &input_page[2*num_parts_per_page];
		vx = &input_page[3*num_parts_per_page];
		vy = &input_page[4*num_parts_per_page];
		vz = &input_page[5*num_parts_per_page];

		/* allocate buffer space for particles on other processors */
		for ( i = 1; i < num_procs; i++ ) {
			page[i] = cart_alloc( 2*nDim*num_parts_per_proc_page*sizeof(particle_float) );
			page_ids[i] = cart_alloc( num_parts_per_proc_page * sizeof(int) );
			count[i] = 0;
			current_page[i] = 0;
		}

		/* start loading actual particle data */
		input = fopen( data_filename, "r" );
		if ( input == NULL ) {
			cart_error( "Unable to open particle file %s for reading!", data_filename );
		}

		if ( timestep_filename != NULL ) {
			for ( i = 1; i < num_procs; i++ ) {
				timestep_page[i] = cart_alloc( num_parts_per_proc_page*sizeof(float) );
			}
		
			timestep_input = fopen( timestep_filename, "r");
			if ( timestep_input == NULL ) {
				cart_error("Unable to open particle dt file %s for reading.", timestep_filename );
			}

			pdt = cart_alloc( num_parts_per_page * sizeof(float) );

			fread( &size, sizeof(int), 1, timestep_input );

			if ( size != num_particles_total * sizeof(int) ) {
				reorder( (char *)&size, sizeof(int) );

				if ( size != num_particles_total * sizeof(int) ) {
					cart_error("Error: particle dt file %s is corrupt.", timestep_filename );
				} else {
					dt_endian = 1;
				}
			} else {
				dt_endian = 0;
			}
		}

		current_id = 0;
		current_type = 0;

		for ( i = 0; i < num_pages; i++ ) {
			if ( i == num_pages - 1 ) {
				num_parts_in_page = num_particles_total - num_parts_per_page*(num_pages-1);
			} else {
				num_parts_in_page = num_parts_per_page;
			}

			num_read = fread( input_page, sizeof(particle_float), 2*nDim*num_parts_per_page, input );
			if ( num_read != 2*nDim*num_parts_per_page ) {
				cart_error("Error reading from particle file %s: insufficient data", data_filename );
			}

			if ( endian ) {
				for ( j = 0; j < num_parts_in_page; j++ ) {
					reorder( (char *)&x[j], sizeof(particle_float) );
					reorder( (char *)&y[j], sizeof(particle_float) );
					reorder( (char *)&z[j], sizeof(particle_float) );
					reorder( (char *)&vx[j], sizeof(particle_float) );
					reorder( (char *)&vy[j], sizeof(particle_float) );
					reorder( (char *)&vz[j], sizeof(particle_float) );
				}
			}

			if ( timestep_filename != NULL ) {
				num_read = fread( pdt, sizeof(float), num_parts_in_page, timestep_input );

				if ( num_read != num_parts_in_page ) {
					cart_error("Error: ran out of particle dt's in %s.", timestep_filename );
				}

				if ( dt_endian ) {
					for ( j = 0; j < num_parts_in_page; j++ ) {
						reorder( (char *)&pdt[j], sizeof(float) );
					}
				}
                        }

#ifdef STARFORM
			if ( stellar_filename != NULL ) {
				/* have we reached star particles yet? */
				if ( particle_id_is_star(current_id+num_parts_in_page-1) ) {
					/* load page values for stars */
					first_star = max( current_id - particle_species_indices[num_particle_species-1], 0 );
					cart_assert( first_star >= 0 && first_star < num_stars );

					first_star_index = first_star + particle_species_indices[num_particle_species-1] - current_id;
					cart_assert( first_star_index >= 0 && first_star_index < num_parts_per_page );

					num_stars_to_read = min( num_parts_in_page - first_star_index, 
						particle_species_indices[num_particle_species] - first_star );

					cart_assert( num_stars_to_read >= 0 && num_stars_to_read <= num_parts_in_page );

					var_first = 2*sizeof(int) + 2*sizeof(double) +		/* t, a */
						+ 2*sizeof(int) + sizeof(int) +			/* num_stars */
						+ 2*sizeof(int) + 2*sizeof(double) +		/* ws_old, ws_oldi */
						+ sizeof(int) + first_star*sizeof(float);	/* skip over first pw */

					seek_amount = (num_stars - num_stars_to_read)*sizeof(float) + 2*sizeof(int);

					/* weights */
					if ( fseek( stellar_input, var_first, SEEK_SET ) ) {
						cart_error("Error seeking to %u in file %s",
							num_stars*sizeof(float) + sizeof(int), stellar_filename );
					}

					num_read = fread( &pw[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
					if ( num_read != num_stars_to_read ) {
						cart_error("Error reading from stellar file %s", stellar_filename );
					}

					if ( stellar_endian ) {
						for ( j = 0; j < num_stars_to_read; j++ ) {
							reorder( (char *)&pw[first_star_index+j], sizeof(float) );
						}
					}
					
					/* original weights */
					if ( fseek( stellar_input, seek_amount, SEEK_CUR ) ) {
						cart_error("Error seeking %u bytes in file %s", seek_amount, stellar_filename );
					}

					num_read = fread( &pw0[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
					if ( num_read != num_stars_to_read ) {
						cart_error("Error reading from stellar file %s", stellar_filename );
					}

					if ( stellar_endian ) {
						for ( j = 0; j < num_stars_to_read; j++ ) {
							reorder( (char *)&pw0[first_star_index+j], sizeof(float) );
						}
					}

					/* birth times */
					if ( fseek( stellar_input, seek_amount, SEEK_CUR ) ) {
						cart_error("Error seeking %u bytes in file %s", seek_amount, stellar_filename );
					}

					num_read = fread( &tbirth[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
					if ( num_read != num_stars_to_read ) {
						cart_error("Error reading from stellar file %s", stellar_filename );
					}

					if ( stellar_endian ) {
						for ( j = 0; j < num_stars_to_read; j++ ) {
							reorder( (char *)&tbirth[first_star_index+j], sizeof(float) );
						}
					}

#ifdef ENRICH
					/* metallicity */
					if ( fseek( stellar_input, seek_amount, SEEK_CUR ) ) {
						cart_error("Error seeking %u bytes in file %s", seek_amount, stellar_filename );
					}

					num_read = fread( &zstII[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
					if ( num_read != num_stars_to_read ) {
						cart_error("Error reading from stellar file %s", stellar_filename );
					}
                                                                                                                                                            
					if ( stellar_endian ) {
						for ( j = 0; j < num_stars_to_read; j++ ) {
							reorder( (char *)&zstII[first_star_index+j], sizeof(float) );
						}
					}

#endif /* ENRICH */
#ifdef ENRICH_SNIa
					if ( fseek( stellar_input, seek_amount, SEEK_CUR ) ) {
						cart_error("Error seeking %u bytes in file %s", seek_amount, stellar_filename );
					}

					num_read = fread( &zstIa[first_star_index], sizeof(float), num_stars_to_read, stellar_input );
					if ( num_read != num_stars_to_read ) {
						cart_error("Error reading from stellar file %s", stellar_filename );
					}

					if ( stellar_endian ) {
						for ( j = 0; j < num_stars_to_read; j++ ) {
							reorder( (char *)&zstIa[first_star_index+j], sizeof(float) );
						}
					}
#endif /* ENRICH_SNIa */
				}
			}
#endif /* STARFORM */

			for ( j = 0; j < num_parts_in_page; j++ ) {
				/* convert to our coordinates 0->num_grid */
				x[j] -= grid_shift;
				y[j] -= grid_shift;
				z[j] -= grid_shift;

				if ( grid_change_flag ) {
					x[j] *= rfact;
					y[j] *= rfact;
					z[j] *= rfact;
				}

				if ( nbody_flag || grid_change_flag ) {
					vx[j] *= vfact;
					vy[j] *= vfact;
					vz[j] *= vfact;
				}
				
				/* enforce periodic boundary conditions */
				if ( x[j] < 0.0 ) {
					x[j] += (double)num_grid;
				} else if ( x[j] >= (double)num_grid ) {
					x[j] -= (double)num_grid;
				}

				if ( y[j] < 0.0 ) {
					y[j] += (double)num_grid;
				} else if ( y[j] >= (double)num_grid ) {
					y[j] -= (double)num_grid;
				}

				if ( z[j] < 0.0 ) {
					z[j] += (double)num_grid;
				} else if ( z[j] >= (double)num_grid ) {
					z[j] -= (double)num_grid;
				}

				coords[0] = (int)(x[j]);
				coords[1] = (int)(y[j]);
				coords[2] = (int)(z[j]);

				index = sfc_index( coords );
				cart_assert( index >= 0 && index < max_sfc_index );

				/* check if we're supposed to read in this particle */
				if ( num_sfcs > 0 && sfc_map[index] == 0 ) {
					proc = -1;
				} else if ( num_procs == 1 ) {
					proc = MASTER_NODE;
				} else {
					proc = processor_owner( index );
				}

				if ( current_id >= particle_species_indices[current_type+1] ) {
					current_type++;
				}

				if ( proc == MASTER_NODE ) {
					ipart = particle_alloc( current_id );
					cart_assert( ipart >= 0 && ipart < num_particles );

					particle_x[ipart][0] = x[j];
					particle_x[ipart][1] = y[j];
					particle_x[ipart][2] = z[j];
					particle_v[ipart][0] = vx[j];
					particle_v[ipart][1] = vy[j];
					particle_v[ipart][2] = vz[j];

					cart_assert( particle_specie( current_id ) == current_type );

					particle_t[ipart] = tl[min_level];
					if ( timestep_filename == NULL ) {
						particle_dt[ipart] = dt;
					} else {
						particle_dt[ipart] = pdt[j];
					}

#ifdef STARFORM
					if ( stellar_filename != NULL && particle_id_is_star( current_id ) ) {
						cart_assert( ipart >= 0 && ipart < num_star_particles );
						cart_assert( particle_is_star(ipart) );

						star_tbirth[ipart] = tbirth[j];
						star_initial_mass[ipart] = pw0[j];
						particle_mass[ipart] = pw[j];

						cart_assert( star_initial_mass[ipart] > 0.0 );

#ifdef ENRICH
						star_metallicity_II[ipart] = zstII[j];
#ifdef ENRICH_SNIa
						star_metallicity_Ia[ipart] = zstIa[j];
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
					} else {
						particle_mass[ipart] = particle_species_mass[particle_specie(current_id)];
					}
#else
					particle_mass[ipart] = particle_species_mass[particle_specie(current_id)];
#endif /* STARFORM */

				} else if ( proc >= 0 && proc < num_procs ) {
					/* add the particle to a processor page */
					page_ids[proc][count[proc]] = current_id;
					page[proc][2*nDim*count[proc]] = x[j];
					page[proc][2*nDim*count[proc]+1] = y[j];
					page[proc][2*nDim*count[proc]+2] = z[j];
					page[proc][2*nDim*count[proc]+3] = vx[j];
					page[proc][2*nDim*count[proc]+4] = vy[j];
					page[proc][2*nDim*count[proc]+5] = vz[j];

					if ( timestep_filename != NULL ) {
						timestep_page[proc][count[proc]] = pdt[j];
					}

#ifdef STARFORM
					if ( stellar_filename != NULL && particle_id_is_star( current_id ) ) {
						star_page[proc][num_star_vars[proc]++] = pw[j];
						star_page[proc][num_star_vars[proc]++] = pw0[j];
						star_page[proc][num_star_vars[proc]++] = tbirth[j];
					
#ifdef ENRICH
						star_page[proc][num_star_vars[proc]++] = zstII[j];
#endif /* ENRICH */
#ifdef ENRICH_SNIa
						star_page[proc][num_star_vars[proc]++] = zstIa[j];
#endif /* ENRICH_SNIa */
					}
#endif /* STARFORM */

					count[proc]++;

					if ( count[proc] == num_parts_per_proc_page ) {
						MPI_Send( page_ids[proc], num_parts_per_proc_page, MPI_INT, proc, 
							current_page[proc], MPI_COMM_WORLD );
						MPI_Send( page[proc], 2*nDim*num_parts_per_proc_page, 
							MPI_PARTICLE_FLOAT, proc, current_page[proc], MPI_COMM_WORLD );

						if ( timestep_filename != NULL ) {
							MPI_Send( timestep_page[proc], num_parts_per_proc_page,
								MPI_FLOAT, proc, current_page[proc], MPI_COMM_WORLD );
						}

#ifdef STARFORM
						if ( stellar_filename != NULL ) {
							MPI_Send( star_page[proc], num_star_vars[proc],
								MPI_FLOAT, proc, current_page[proc], MPI_COMM_WORLD );
							num_star_vars[proc] = 0;
						}
#endif /* STARFORM */
						count[proc] = 0;
						current_page[proc]++;
					}
				}

				current_id++;
			}
		}
	
		/* send final pages */
		for ( proc = 1; proc < num_procs; proc++ ) {
			MPI_Send( page_ids[proc], count[proc], MPI_INT, proc, current_page[proc], MPI_COMM_WORLD );
			MPI_Send( page[proc], 2*nDim*count[proc], MPI_PARTICLE_FLOAT, proc, current_page[proc], MPI_COMM_WORLD );

			if ( timestep_filename != NULL ) {
				MPI_Send( timestep_page[proc], count[proc], MPI_FLOAT, proc, current_page[proc], MPI_COMM_WORLD );
				cart_free( timestep_page[proc] );
			}

#ifdef STARFORM
			if ( stellar_filename != NULL ) {
				MPI_Send( star_page[proc], num_star_vars[proc], MPI_FLOAT, proc, current_page[proc], MPI_COMM_WORLD );
				cart_free( star_page[proc] );
			}
#endif /* STARFORM */

			cart_free( page_ids[proc] );
			cart_free( page[proc] );
		}

#ifdef STARFORM
		if ( stellar_filename != NULL ) {
#ifdef ENRICH_SNIa
			cart_free( zstIa );
#endif /* ENRICH_SNIa */
#ifdef ENRICH
			cart_free( zstII );
#endif /* ENRICH */

			cart_free( tbirth );
			cart_free( pw0 );
			cart_free( pw );
		}
#endif /* STARFORM */

		cart_free( input_page );

		fclose(input);

		if ( num_sfcs > 0 ) {
			cart_free( sfc_map );
		}
	} else {
                MPI_Bcast( &aexp[min_level], 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &tl[min_level], 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &dt, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &a_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &Omega0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &OmegaL0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &Omegab0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &hubble, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &step, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &num_row, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &num_particle_species, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &num_particles_total, 1, MPI_LONG, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &num_parts_per_page, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &num_parts_per_proc_page, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( &num_pages, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( particle_species_mass, num_particle_species, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( particle_species_num, num_particle_species, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
                MPI_Bcast( particle_species_indices, num_particle_species+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

		current_page[local_proc_id] = 0;

		page_ids[local_proc_id] = cart_alloc( num_parts_per_proc_page * sizeof(int) );
		page[local_proc_id] = cart_alloc( 2*nDim*num_parts_per_proc_page * sizeof(particle_float) );
		
		if ( timestep_filename != NULL ) {
			timestep_page[local_proc_id] = cart_alloc( num_parts_per_proc_page * sizeof(float) );
		}

#ifdef STARFORM
		if ( stellar_filename != NULL ) {
			star_page[local_proc_id] = cart_alloc( num_star_variables*num_parts_per_proc_page*sizeof(float) );
		}
#endif /* STARFORM */

		count[local_proc_id] = num_parts_per_proc_page;
		while ( count[local_proc_id] == num_parts_per_proc_page ) {
			MPI_Recv( page_ids[local_proc_id], num_parts_per_proc_page, MPI_INT, 
				MASTER_NODE, current_page[local_proc_id], MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[local_proc_id] );

			MPI_Recv( page[local_proc_id], 2*nDim*num_parts_per_proc_page, 
				MPI_PARTICLE_FLOAT, MASTER_NODE, current_page[local_proc_id], MPI_COMM_WORLD, &status );

			if ( timestep_filename != NULL ) {
				MPI_Recv( timestep_page[local_proc_id], num_parts_per_proc_page,
					MPI_FLOAT, MASTER_NODE, current_page[local_proc_id], MPI_COMM_WORLD, &status );
			}

#ifdef STARFORM
			if ( stellar_filename != NULL ) {
				MPI_Recv( star_page[local_proc_id], num_star_variables*num_parts_per_proc_page,
					MPI_FLOAT, MASTER_NODE, current_page[local_proc_id], MPI_COMM_WORLD, &status );
				num_star_vars[local_proc_id] = 0;
			}
#endif /* STARFORM */
			current_page[local_proc_id]++;

			for ( i = 0; i < count[local_proc_id]; i++ ) {
				ipart = particle_alloc( page_ids[local_proc_id][i] );
				cart_assert( ipart >= 0 && ipart < num_particles );

				particle_x[ipart][0] = page[local_proc_id][2*nDim*i];
				particle_x[ipart][1] = page[local_proc_id][2*nDim*i+1];
				particle_x[ipart][2] = page[local_proc_id][2*nDim*i+2];
				particle_v[ipart][0] = page[local_proc_id][2*nDim*i+3];
				particle_v[ipart][1] = page[local_proc_id][2*nDim*i+4];
				particle_v[ipart][2] = page[local_proc_id][2*nDim*i+5];

				particle_t[ipart] = tl[min_level];

				if ( timestep_filename == NULL ) {
					particle_dt[ipart] = dt;
				} else {
					particle_dt[ipart] = timestep_page[local_proc_id][i];
				}

#ifdef STARFORM
				if ( particle_id_is_star( particle_id[ipart] ) ) {
					cart_assert( ipart >= 0 && ipart < num_star_particles );
					cart_assert( particle_is_star(ipart) );

					particle_mass[ipart] = star_page[local_proc_id][num_star_vars[local_proc_id]++];
					star_initial_mass[ipart] = star_page[local_proc_id][num_star_vars[local_proc_id]++];
					star_tbirth[ipart] = star_page[local_proc_id][num_star_vars[local_proc_id]++];

#ifdef ENRICH
					star_metallicity_II[ipart] = star_page[local_proc_id][num_star_vars[local_proc_id]++];
#endif /* ENRICH */
#ifdef ENRICH_SNIa
					star_metallicity_Ia[ipart] = star_page[local_proc_id][num_star_vars[local_proc_id]++];
#endif /* ENRICH_SNIa */
				} else {
					particle_mass[ipart] = particle_species_mass[particle_specie(particle_id[ipart])];
				}
#else
				particle_mass[ipart] = particle_species_mass[particle_specie(particle_id[ipart])];
#endif /* STARFORM */
			}
		}

#ifdef STARFORM
		if ( stellar_filename != NULL ) {
			cart_free( star_page[local_proc_id] );
		}
#endif /* STARFORM */

		if ( timestep_filename != NULL ) {
			cart_free( timestep_page[local_proc_id] );
		}

		cart_free( page[local_proc_id] );
		cart_free( page_ids[local_proc_id] );
	}

	cart_debug("num_local_particles = %u", num_local_particles );
	cart_assert(num_local_particles >= 0 );
	build_particle_list();
}

#endif /* PARTICLES */

#ifdef HYDRO

#ifdef HYDRO_TRACERS
int compare_tracer_ids( const void *a, const void *b ) {
	return ( tracer_id[*(int *)a] - tracer_id[*(int *)b] );
}

void write_hydro_tracers( char *filename ) {
	int i, j, k;
	FILE *output;
	int size;
	int icell, index;
	float adum, ainit;
	float boxh, Om0, Oml0, Omb0, h;
	int *tracer_order;
	char label[32];
	int num_tracers_in_page;
	int num_tracers_per_page, num_tracers_per_proc_page, num_pages;
	int num_tracers_written, current_id, local_count;
	double *output_page;
	float *output_vars_page;
	int *page_ids[MAX_PROCS];
	float *vars_page[MAX_PROCS];
	double *page[MAX_PROCS];
	int num_pages_sent;
	int num_pages_received[MAX_PROCS];
        int count[MAX_PROCS];
        int pos[MAX_PROCS];
        MPI_Status status;

	cart_debug("writing hydro tracers");

	num_tracers_per_page = num_tracer_row*num_tracer_row;
	num_tracers_per_proc_page = num_tracers_per_page/num_procs;
	num_pages = (num_tracers_total-1) / num_tracers_per_page + 1;

	/* construct mapping of id -> tracer id */
	tracer_order = cart_alloc( num_local_tracers * sizeof(int) );
	j = 0;
	for ( i = 0; i < num_tracers; i++ ) {
		if ( tracer_id[i] != NULL_TRACER ) {
			tracer_order[j++] = i;
		}
	}
	cart_assert( j == num_local_tracers );

	qsort( tracer_order, num_local_tracers, sizeof(int), compare_tracer_ids );

	if ( local_proc_id == MASTER_NODE ) { 
		output = fopen( filename, "w" );
		if ( output == NULL ) {
			cart_error("Unable to open hydro tracer file %s.", filename );
		}

		/* write header to file */
		size = 256*sizeof(char);
		fwrite(&size, sizeof(int), 1, output );
		fwrite(jobname, sizeof(char), 256, output );
		fwrite(&size, sizeof(int), 1, output );
		
		/* istep, t, dt, adum, ainit */
		adum = aexp[min_level];
		ainit = a_init;
		size = sizeof(int) + 2*sizeof(double) + 2*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &step, sizeof(int), 1, output );
		fwrite( &tl[min_level], sizeof(double), 1, output );
		fwrite( &dtl[min_level], sizeof(double), 1, output );
		fwrite( &adum, sizeof(float), 1, output );
		fwrite( &ainit, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* boxh, Om0, Oml0, Omb0, hubble */
		boxh = Lbox;
		Om0 = Omega0;
		Oml0 = OmegaL0;
		Omb0 = Omegab0;
		h = hubble;

		size = 5*sizeof(float);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &boxh, sizeof(float), 1, output );
		fwrite( &Om0, sizeof(float), 1, output );
		fwrite( &Oml0, sizeof(float), 1, output );
		fwrite( &Omb0, sizeof(float), 1, output );
		fwrite( &h, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* num_grid */
		size = sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		size = num_grid;
		fwrite( &size, sizeof(int), 1, output );
		size = sizeof(int);
		fwrite( &size, sizeof(int), 1, output );

		/* number of tracers & page_size */
		size = 2*sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &num_tracers_total, sizeof(int), 1, output );
		fwrite( &num_tracer_row, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* hydro variables traced */
		size = sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &num_hydro_vars_traced, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		size = 32*num_hydro_vars_traced*sizeof(int);
		fwrite( &size, sizeof(int), 1, output );

		for ( i = 0; i < num_hydro_vars_traced; i++ ) {
			snprintf( label, 32, "%s", hydro_vars_traced_labels[i] );
			fwrite( label, sizeof(char), 32, output );
		}

		fwrite( &size, sizeof(int), 1, output );

		/* allocate space to receive pages of tracers from other procs */
		for ( i = 1; i < num_procs; i++ ) {
			page_ids[i] = cart_alloc( num_tracers_per_proc_page * sizeof(int) );
			page[i] = cart_alloc( nDim*num_tracers_per_proc_page*sizeof(double) );
			vars_page[i] = cart_alloc( num_hydro_vars_traced*num_tracers_per_proc_page*sizeof(float) );
		}

		/* allocate actual page which will be written */
		output_page = cart_alloc( nDim*num_tracers_per_page*sizeof(double) );

		/* receive initial pages */
		for ( i = 1; i < num_procs; i++ ) {
			num_pages_received[i] = 0;
			MPI_Recv( page_ids[i], num_tracers_per_proc_page, 
				MPI_INT, i, num_pages_received[i], MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[i] );
			cart_assert( count[i] >= 0 && count[i] <= num_tracers_per_proc_page );
			MPI_Recv( page[i], nDim*num_tracers_per_proc_page, MPI_DOUBLE, i, 
				num_pages_received[i], MPI_COMM_WORLD, &status );
			pos[i] = 0;
			num_pages_received[i]++;
		}

		current_id = 0;
		local_count = 0;
		for ( i = 0; i < num_pages; i++ ) {
			/* construct page */
			num_tracers_written = 0;
			if ( i == num_pages - 1 ) {
				num_tracers_in_page = num_tracers_total - num_tracers_per_page*(num_pages-1);
			} else {
				num_tracers_in_page = num_tracers_per_page;
			}

			cart_assert( num_tracers_in_page >= 0 && num_tracers_in_page <= num_tracers_per_page );

			while ( num_tracers_written < num_tracers_in_page ) {
				/* add from our list */
				while ( num_tracers_written < num_tracers_in_page &&
						local_count < num_local_tracers &&
						tracer_id[tracer_order[local_count]] == current_id ) {
					index = tracer_order[local_count];
					
					output_page[num_tracers_written]			= tracer_x[index][0] + 1.0;
					output_page[num_tracers_per_page+num_tracers_written] 	= tracer_x[index][1] + 1.0;
					output_page[2*num_tracers_per_page+num_tracers_written]	= tracer_x[index][2] + 1.0;

					local_count++;
					num_tracers_written++;
					current_id++;
				}

				/* add from each processor in turn */
				for ( j = 1; j < num_procs; j++ ) {
					while ( num_tracers_written < num_tracers_in_page &&
							pos[j] < count[j] &&
							page_ids[j][pos[j]] == current_id ) {

						output_page[num_tracers_written]			= page[j][nDim*pos[j]] + 1.0;
						output_page[num_tracers_per_page+num_tracers_written]	= page[j][nDim*pos[j]+1] + 1.0;
						output_page[2*num_tracers_per_page+num_tracers_written] = page[j][nDim*pos[j]+2] + 1.0;

						num_tracers_written++;
						current_id++;
						pos[j]++;

						/* if we've run out, refill page */
						if ( pos[j] == count[j] && count[j] == num_tracers_per_proc_page ) {
							MPI_Recv( page_ids[j], num_tracers_per_proc_page, MPI_INT,
								j, num_pages_received[j], MPI_COMM_WORLD, &status );
							MPI_Get_count( &status, MPI_INT, &count[j] );
							cart_assert( count[j] >= 0 && count[j] <= num_tracers_per_proc_page );
							MPI_Recv( page[j], nDim*num_tracers_per_proc_page,
								MPI_DOUBLE, j, num_pages_received[j], 
								MPI_COMM_WORLD, &status );
							pos[j] = 0;
							num_pages_received[j]++;
						}
					}
				}
			}

			cart_assert( num_tracers_written == num_tracers_in_page );

			/* write page */
			fwrite( output_page, sizeof(double), nDim*num_tracers_per_page, output );
		}

		cart_free( output_page );

                /* allocate actual page which will be written */
		output_vars_page = cart_alloc( num_hydro_vars_traced*num_tracers_per_page*sizeof(float) );

		/* receive initial pages */
		for ( i = 1; i < num_procs; i++ ) {
			MPI_Recv( page_ids[i], num_tracers_per_proc_page, MPI_INT, 
				i, num_pages_received[i], MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[i] );
			cart_assert( count[i] >= 0 && count[i] <= num_tracers_per_proc_page );
			MPI_Recv( page[i], num_hydro_vars_traced*num_tracers_per_proc_page, 
				MPI_FLOAT, i, num_pages_received[i], MPI_COMM_WORLD, &status );
			pos[i] = 0;
			num_pages_received[i]++;
		}

		current_id = 0;
		local_count = 0;
		for ( i = 0; i < num_pages; i++ ) {
			/* construct page */
			num_tracers_written = 0;
			if ( i == num_pages - 1 ) {
				num_tracers_in_page = num_tracers_total - num_tracers_per_page*(num_pages-1);
			} else {
				num_tracers_in_page = num_tracers_per_page;
			}

			cart_assert( num_tracers_in_page >= 0 && num_tracers_in_page <= num_tracers_per_page );

			while ( num_tracers_written < num_tracers_in_page ) {
				/* add from our list */
				while ( num_tracers_written < num_tracers_in_page &&
						local_count < num_local_tracers &&
						tracer_id[tracer_order[local_count]] == current_id ) {

					index = tracer_order[local_count];
					icell = cell_find_position( tracer_x[index] );

					for ( j = 0; j < num_hydro_vars_traced; j++ ) {
						cart_assert( j*num_tracers_per_page + num_tracers_written < num_hydro_vars_traced*num_tracers_per_page );
						output_vars_page[j*num_tracers_per_page + num_tracers_written] = 
								cell_vars[icell][ hydro_vars_traced[j] ];
					}

					local_count++;
					num_tracers_written++;
					current_id++;
				}

				/* add from each processor in turn */
				for ( j = 1; j < num_procs; j++ ) {
					while ( num_tracers_written < num_tracers_in_page &&
							pos[j] < count[j] &&
							page_ids[j][pos[j]] == current_id ) {

						for ( k = 0; k < num_hydro_vars_traced; k++ ) {
							output_vars_page[k*num_tracers_per_page + num_tracers_written ] =
								page[j][num_hydro_vars_traced*pos[j]+k];
						}
					
						num_tracers_written++;
						current_id++;
						pos[j]++;

						/* if we've run out, refill page */
						if ( pos[j] == count[j] && count[j] == num_tracers_per_proc_page ) {
							MPI_Recv( page_ids[j], num_tracers_per_proc_page, MPI_INT,
								j, num_pages_received[j], MPI_COMM_WORLD, &status );
							MPI_Get_count( &status, MPI_INT, &count[j] );
							cart_assert( count[j] >= 0 && count[j] <= num_tracers_per_proc_page );
							MPI_Recv( page[j], num_hydro_vars_traced*num_tracers_per_proc_page,
								MPI_FLOAT, j, num_pages_received[j], 
								MPI_COMM_WORLD, &status );
							pos[j] = 0;
							num_pages_received[j]++;
						}
					}
				}
			}

			cart_assert( num_tracers_written == num_tracers_in_page );

			/* write page */
			fwrite( output_vars_page, sizeof(float), num_hydro_vars_traced*num_tracers_per_page, output );
		}		

		cart_free( output_vars_page );

		for ( i = 1; i < num_procs; i++ ) {
			cart_free( vars_page[i] );
			cart_free( page[i] );
			cart_free( page_ids[i] );
		}
		fclose( output );
	} else {
		num_pages_sent = 0;
		page_ids[local_proc_id] = cart_alloc( num_tracers_per_proc_page*sizeof(int) );
		page[local_proc_id] = cart_alloc( nDim*num_tracers_per_proc_page*sizeof(double) );
		vars_page[local_proc_id] = cart_alloc( num_hydro_vars_traced*num_tracers_per_proc_page*sizeof(float) );

		num_tracers_written = 0;
		do {
			/* pack page */
			local_count = 0;
			for ( i = 0; i < num_tracers_per_proc_page && num_tracers_written < num_local_tracers; i++ ) {
				index = tracer_order[num_tracers_written];
				cart_assert( index >= 0 && index < num_tracers );
				cart_assert( tracer_id[index] >= 0 && tracer_id[index] < num_tracers_total );
				cart_assert( i >= 0 && i < num_tracers_per_proc_page );
				page_ids[local_proc_id][i] = tracer_id[index];
				page[local_proc_id][local_count++] = tracer_x[index][0];
				page[local_proc_id][local_count++] = tracer_x[index][1];
				page[local_proc_id][local_count++] = tracer_x[index][2];
				num_tracers_written++;
			}

			cart_assert( i <= num_tracers_per_proc_page );

			/* send the page */
			MPI_Send( page_ids[local_proc_id], i, MPI_INT, 
				MASTER_NODE, num_pages_sent, MPI_COMM_WORLD );
			MPI_Send( page[local_proc_id], local_count, MPI_DOUBLE, 
				MASTER_NODE, num_pages_sent, MPI_COMM_WORLD );
			num_pages_sent++;
		} while ( i == num_tracers_per_proc_page );

		cart_free( page[local_proc_id] );

		/* send hydro traced variables */
		num_tracers_written = 0;
		do {
			local_count = 0;
			for ( i = 0; i < num_tracers_per_proc_page && num_tracers_written < num_local_tracers; i++ ) {
				index = tracer_order[num_tracers_written];
				cart_assert( index >= 0 && index < num_tracers );
				cart_assert( tracer_id[index] >= 0 && tracer_id[index] < num_tracers_total );

				page_ids[local_proc_id][i] = tracer_id[index];
				icell = cell_find_position( tracer_x[index] );

				for ( j = 0; j < num_hydro_vars_traced; j++ ) {
					vars_page[local_proc_id][local_count++] = cell_vars[icell][ hydro_vars_traced[j] ];
				}

				num_tracers_written++;
			}

			cart_assert( i <= num_tracers_per_proc_page );

			MPI_Send( page_ids[local_proc_id], i, MPI_INT, 	
				MASTER_NODE, num_pages_sent, MPI_COMM_WORLD );
			MPI_Send( vars_page[local_proc_id], local_count, 
				MPI_FLOAT, MASTER_NODE, num_pages_sent, MPI_COMM_WORLD );
			num_pages_sent++;
		} while ( i == num_tracers_per_proc_page );

		cart_free( vars_page[local_proc_id] );
		cart_free( page_ids[local_proc_id] );
	}

	cart_free( tracer_order );
}

void read_hydro_tracers( char *filename ) {
	int i, j, k;
	FILE *input;
	int size, endian;
	char tracerjobname[256];
	double tmpt, tmpdt;
	int tmp_num_grid, tmp_num_hydro_vars_traced;
	int num_read;
	int proc;
	int current_id;
	int tracer;
	int index, icell;
	int coords[nDim];
	float ainit, adum;
	float boxh, Om0, Oml0, Omb0, h;
	int *tracer_order;
	char label[32];
	int num_tracers_in_page;
	int num_tracers_per_page, num_tracers_per_proc_page, num_pages;
	int num_tracers_written, local_count;
	double *input_page;
	double *x, *y, *z;
	int *page_ids[MAX_PROCS];
	double *page[MAX_PROCS];
	int count[MAX_PROCS];
	int current_page[MAX_PROCS];
	int pos[MAX_PROCS];
	MPI_Status status;

	if ( local_proc_id == MASTER_NODE ) {
		input = fopen( filename, "r" );
		if ( input == NULL ) {
			cart_error("Unable to open hydro tracer file %s.", filename );
		}

		endian = 0;

                /* write header to file */
		fread( &size, sizeof(int), 1, input );

		if ( size != 256*sizeof(char) ) {
			reorder( (char *)&size, sizeof(int) );

			if ( size == 256*sizeof(char) ) {
				endian = 1;
				cart_debug("switching endianness of tracer file");
			} else {
				cart_debug("hydro tracer file %s appears to be corrupt!", filename );
			}
		}

		fread( tracerjobname, sizeof(char), 256, input );
		fread(&size, sizeof(int), 1, input );

		cart_debug("tracer job name: %s", tracerjobname );

		fread( &size, sizeof(int), 1, input );
		fread( &step, sizeof(int), 1, input );
		fread( &tmpt, sizeof(double), 1, input );
		fread( &tmpdt, sizeof(double), 1, input );
		fread( &adum, sizeof(float), 1, input );
		fread( &ainit, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		/* boxh, Om0, Oml0, Omb0, hubble */
		fread( &size, sizeof(int), 1, input );
		fread( &boxh, sizeof(float), 1, input );
		fread( &Om0, sizeof(float), 1, input );
		fread( &Oml0, sizeof(float), 1, input );
		fread( &Omb0, sizeof(float), 1, input );
		fread( &h, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		/* num_grid */
		fread( &size, sizeof(int), 1, input );
		fread( &tmp_num_grid, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );	

		if ( endian ) {
			reorder( (char *)&tmp_num_grid, sizeof(int) );
		}

		if ( tmp_num_grid != num_grid ) {
			cart_error("num_grid in %s doesn't match compiled: %u vs %u\n", filename, tmp_num_grid, num_grid );
		}

		/* number of tracers & page_size */
		fread( &size, sizeof(int), 1, input );
		fread( &num_tracers_total, sizeof(int), 1, input );
		fread( &num_tracer_row, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&num_tracers_total, sizeof(int) );
			reorder( (char *)&num_tracer_row, sizeof(int) );
		}

		/* hydro variables traced */
		fread( &size, sizeof(int), 1, input );
		fread( &tmp_num_hydro_vars_traced, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		fread( &size, sizeof(int), 1, input );

		for ( i = 0; i < num_hydro_vars_traced; i++ ) {
			fread( label, sizeof(char), 32, input );
		}

		fread( &size, sizeof(int), 1, input );
	}

	MPI_Bcast( &num_tracer_row, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &num_tracers_total, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	num_tracers_per_page = num_tracer_row*num_tracer_row;
	num_tracers_per_proc_page = num_tracers_per_page/num_procs;
	num_pages = (num_tracers_total-1) / num_tracers_per_page + 1;

	if ( local_proc_id == MASTER_NODE ) {
		input_page = cart_alloc( nDim*num_tracers_per_page*sizeof(double) );

		x = input_page;
		y = &input_page[num_tracers_per_page];
		z = &input_page[2*num_tracers_per_page];

		/* allocate buffer space for tracers on other processors */
		for ( i = 1; i < num_procs; i++ ) {
			page_ids[i] = cart_alloc( num_tracers_per_proc_page*sizeof(int) );
			page[i] = cart_alloc( nDim*num_tracers_per_proc_page*sizeof(double) );
			count[i] = 0;
			current_page[i] = 0;
		}

		current_id = 0;

		for ( i = 0; i < num_pages; i++ ) {
			if ( i == num_pages - 1 ) {
				num_tracers_in_page = num_tracers_total - num_tracers_per_page*(num_pages-1);
			} else {
				num_tracers_in_page = num_tracers_per_page;
			}

			num_read = fread( input_page, sizeof(double), nDim*num_tracers_per_page, input );

			if ( num_read != nDim*num_tracers_per_page ) {
				cart_error("Error reading from tracer file %s: insufficient data", filename );
			}

			if ( endian ) {
				for ( j = 0; j < num_tracers_in_page; j++ ) {
					reorder( (char *)&x[j], sizeof(double) );
					reorder( (char *)&y[j], sizeof(double) );
					reorder( (char *)&z[j], sizeof(double) );
				}
			}
	
			for ( j = 0; j < num_tracers_in_page; j++ ) {
				/* enforce periodic boundary conditions */
				if ( x[j] < 0.0 ) {
					x[j] += (double)num_grid;
				} else if ( x[j] >= (double)num_grid ) {
					x[j] -= (double)num_grid;
				}

				if ( y[j] < 0.0 ) {
					y[j] += (double)num_grid;
				} else if ( y[j] >= (double)num_grid ) {
					y[j] -= (double)num_grid;
				}

				if ( z[j] < 0.0 ) {
					z[j] += (double)num_grid;
				} else if ( z[j] >= (double)num_grid ) {
					z[j] -= (double)num_grid;
				}

				coords[0] = (int)(x[j]);
				coords[1] = (int)(y[j]);
				coords[2] = (int)(z[j]);
	
				index = sfc_index( coords );
				cart_assert( index >= 0 && index < max_sfc_index );
				proc = processor_owner( index );
	
				if ( proc == MASTER_NODE ) {
					tracer = tracer_alloc( current_id );

					tracer_x[tracer][0] = x[j];
					tracer_x[tracer][1] = y[j];
					tracer_x[tracer][2] = z[j];
				} else {
					page_ids[proc][count[proc]] = current_id;
					page[proc][nDim*count[proc]] = x[j];
					page[proc][nDim*count[proc]+1] = y[j];
					page[proc][nDim*count[proc]+2] = z[j];
					count[proc]++;

					if ( count[proc] == num_tracers_per_proc_page ) {
						MPI_Send( page_ids[proc], num_tracers_per_proc_page, MPI_INT, proc,
							current_page[proc], MPI_COMM_WORLD );
						MPI_Send( page[proc], nDim*num_tracers_per_proc_page,
							MPI_DOUBLE, proc, current_page[proc], MPI_COMM_WORLD );
						count[proc] = 0;
						current_page[proc]++;
					}
				}

				current_id++;
			}
		}

		/* send final pages */
		for ( proc = 1; proc < num_procs; proc++ ) {
			MPI_Send( page_ids[proc], count[proc], MPI_INT, proc, current_page[proc], MPI_COMM_WORLD );
			MPI_Send( page[proc], nDim*count[proc], MPI_DOUBLE, proc, current_page[proc], MPI_COMM_WORLD );
			cart_free( page_ids[proc] );
			cart_free( page[proc] );
		}

		cart_free( input_page );

		fclose( input );
	} else {

		page_ids[local_proc_id] = cart_alloc( num_tracers_per_proc_page * sizeof(int) );
		page[local_proc_id] = cart_alloc( nDim*num_tracers_per_proc_page * sizeof(double) );

		count[local_proc_id] = num_tracers_per_proc_page;
		while ( count[local_proc_id] == num_tracers_per_proc_page ) {
			MPI_Recv( page_ids[local_proc_id], num_tracers_per_proc_page, MPI_INT,
				MASTER_NODE, current_page[local_proc_id], MPI_COMM_WORLD, &status );
			MPI_Get_count( &status, MPI_INT, &count[local_proc_id] );
			MPI_Recv( page[local_proc_id], nDim*num_tracers_per_proc_page,
				MPI_DOUBLE, MASTER_NODE, current_page[local_proc_id], MPI_COMM_WORLD, &status );

			current_page[local_proc_id]++;

			for ( i = 0; i < count[local_proc_id]; i++ ) {
				tracer = tracer_alloc( page_ids[local_proc_id][i] );

				tracer_x[tracer][0] = page[local_proc_id][nDim*i];
				tracer_x[tracer][1] = page[local_proc_id][nDim*i+1];
				tracer_x[tracer][2] = page[local_proc_id][nDim*i+2];
			}
		}

		cart_free( page[local_proc_id] );
		cart_free( page_ids[local_proc_id] );
        }

        cart_debug("num_local_tracers = %u", num_local_tracers );
	build_tracer_list();
}

#endif /* HYDRO_TRACERS */

void read_gas_ic( char *filename ) {
	int i;
	FILE *input;
	int size;
	float boxh, ainit, astep;
	int ncells;
	int endian;
	int coords[nDim];
	int index;
	int proc, icell;
	MPI_Status status;
	int page_count;
	float *input_page;
	int count[MAX_PROCS];
	float *page[MAX_PROCS];
	int *page_indices[MAX_PROCS];
	int var;
	const int num_gas_vars    = 6;
	const int var_index[] = { 	HVAR_GAS_DENSITY, HVAR_MOMENTUM, 
					HVAR_MOMENTUM+1, HVAR_MOMENTUM+2,
					HVAR_GAS_ENERGY, HVAR_INTERNAL_ENERGY };

	if ( local_proc_id == MASTER_NODE ) {
	        input = fopen(filename, "r");
		if ( input == NULL ) {
			cart_error("Unable to open %s for reading!", filename );
		}
	
		fread( &size, sizeof(int), 1, input );
		if ( size != sizeof(float) ) {
			reorder( (char *)&size, sizeof(int) );
			endian = 1;
			if ( size != sizeof(float) ) {
				cart_error("Bad file-format in read_cell_ic");
			} 
		} else {
			endian = 0;
		}
	
		fread( &boxh, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		fread( &size, sizeof(int), 1, input );
		fread( &ainit, sizeof(float), 1, input );
		fread( &astep, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		fread( &size, sizeof(int), 1, input );
		fread( &ncells, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );
	
		if ( endian ) {
			reorder( (char *)&boxh, sizeof(float) );
			reorder( (char *)&ainit, sizeof(float) );
			reorder( (char *)&astep, sizeof(float) );
			reorder( (char *)&ncells, sizeof(int) );
		}

		Lbox = boxh;
		a_init = ainit;

		MPI_Bcast( &Lbox, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &a_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	
		cart_debug("boxh = %f", boxh );
		cart_debug("ainit = %f", ainit );
		cart_debug("astep = %f", astep );
		cart_debug("ncells = %u", ncells );

		if ( ncells != num_root_cells ) {
			cart_error("ncells in %s does not match num_root_cells (%u vs %u)", 
					filename, ncells, num_root_cells );
		}

		input_page = cart_alloc( num_grid*num_grid*sizeof(float) );
                for ( proc = 1; proc < num_procs; proc++ ) {
			page[proc] = cart_alloc( num_grid*num_grid*sizeof(float) );
			page_indices[proc] = cart_alloc( num_grid*num_grid*sizeof(int) );
                }

		for ( var = 0; var < num_gas_vars; var++ ) {
			for ( proc = 1; proc < num_procs; proc++ ) {
				count[proc] = 0;
			}

			fread( &size, sizeof(int), 1, input );
			for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				size = fread( input_page, sizeof(float), num_grid*num_grid, input );
				if ( size != num_grid*num_grid ) {
					cart_error("Error reading from file %s", filename );
				}

				if ( endian ) {
					for ( i = 0; i < num_grid*num_grid; i++ ) {
						reorder( (char *)&input_page[i], sizeof(float) );
					}
				}
	
				page_count = 0;
				for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
					for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
						index = sfc_index( coords );
						proc = processor_owner( index );
	
						if ( proc == local_proc_id ) {
							icell = root_cell_location( index );
							cell_var( icell, var_index[var] ) = input_page[page_count];
						} else {
							/* add cell to send page */
							page[proc][count[proc]] = input_page[page_count];
							page_indices[proc][count[proc]] = index;
							count[proc]++;
	
							if ( count[proc] == num_grid*num_grid ) {
								MPI_Send( page[proc], num_grid*num_grid, MPI_FLOAT, 
									proc, 0, MPI_COMM_WORLD );
								MPI_Send( page_indices[proc], num_grid*num_grid, 
									MPI_INT, proc, 0, MPI_COMM_WORLD );
								count[proc] = 0;
							}
						}
	
						page_count++;
					}
				}
			}
			fread( &size, sizeof(int), 1, input );
			
			/* send last variables */
			for ( proc = 1; proc < num_procs; proc++ ) {
				MPI_Send( page[proc], count[proc], MPI_FLOAT, proc, 0, MPI_COMM_WORLD );
				MPI_Send( page_indices[proc], count[proc], MPI_INT, proc, 0, MPI_COMM_WORLD );
			}
		}
	
        	fclose(input);

		cart_free( input_page );
		for ( proc = 1; proc < num_procs; proc++ ) {
			cart_free( page[proc] );
			cart_free( page_indices[proc] );
		}
	} else {
		MPI_Bcast( &Lbox, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
		MPI_Bcast( &a_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );

		page[local_proc_id] = cart_alloc( num_grid*num_grid*sizeof(float) );
		page_indices[local_proc_id] = cart_alloc( num_grid*num_grid*sizeof(int) );

		for ( var = 0; var < num_gas_vars; var++ ) {
			page_count = num_grid*num_grid;
			while ( page_count == num_grid*num_grid ) {
				MPI_Recv( page[local_proc_id], num_grid*num_grid, MPI_FLOAT, MASTER_NODE, 0, MPI_COMM_WORLD, &status );
				MPI_Get_count( &status, MPI_FLOAT, &page_count );
				MPI_Recv( page_indices[local_proc_id], num_grid*num_grid, MPI_INT, MASTER_NODE, 0, MPI_COMM_WORLD, &status );

				for ( i = 0; i < page_count; i++ ) {
					icell = root_cell_location( page_indices[local_proc_id][i] );
					cell_var( icell, var_index[var] ) = page[local_proc_id][i];
				}
			}
		}

		cart_free( page[local_proc_id] );
		cart_free( page_indices[local_proc_id] );
	}

	/* set gas gamma on root level */
	for ( i = 0; i < num_cells_per_level[min_level]; i++ ) {
		cell_gas_gamma(i) = gamma;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		cell_electron_internal_energy(i) = cell_gas_internal_energy(i)*wmu/wmu_e;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
	}
}

void write_grid_binary( char *filename ) {
	int i, j, k, m;
	int size;
	FILE *output;
	float adum, ainit;
	float boxh, Om0, Oml0, Omb0, h;
	int minlevel, maxlevel;
	int nextras;
	int *cellrefined;
	float *cellhvars;
	float *cellvars;
	int *order;
	int *current_level;
	int proc;
	int level;
	int icell, ioct;
	int page_size;
	int ncell0, sfc_order;
	int current_level_count;
	int next_level_count;
	int page;
	int page_count;
	long total_cells[max_level-min_level+1];
	MPI_Status status;
	char parallel_filename[256];
	int proc_num_cells[MAX_PROCS*(max_level-min_level+1)];
	int file_index, file_parent, file_num_procs;

	/* ensure consistency of num_output_files */
	num_output_files = min( num_output_files, num_procs );
	num_output_files = max( num_output_files, 1 );

	page_size = num_grid*num_grid;

	/* determine parallel output options */
	file_index = local_proc_id * num_output_files / num_procs;

	file_parent = local_proc_id;
	while ( file_parent > 0 && 
			(file_parent-1)*num_output_files / num_procs == file_index ) {
		file_parent--;
	}	

	file_num_procs = 1;
	while ( file_parent + file_num_procs < num_procs &&
			(file_parent+file_num_procs)*num_output_files / num_procs == file_index ) {
		file_num_procs++;
	}	

	cellrefined = cart_alloc( page_size * sizeof(int) );
	cellhvars = cart_alloc( num_hydro_vars * page_size * sizeof(float) );
	cellvars = cart_alloc( 2 * page_size * sizeof(float) );

	/* get number of cells on each level */
	MPI_Allgather( num_cells_per_level, max_level-min_level+1, MPI_INT,
		proc_num_cells, max_level-min_level+1, MPI_INT, MPI_COMM_WORLD );

	minlevel = min_level;
	maxlevel = max_level_now_global();

	/* open file handle if parent of parallel file */
	if ( local_proc_id == file_parent ) {
		if ( num_output_files == 1 ) {
			output = fopen(filename,"w");
		} else {
			sprintf( parallel_filename, "%s.%03u", filename, file_index );
			output = fopen(parallel_filename, "w");
		}

		if ( output == NULL ) {
			cart_error( "Unable to open file %s for writing!", filename );
		}

		for ( level = min_level; level <= max_level; level++ ) {
			total_cells[level] = 0;
			for ( proc = local_proc_id; proc < local_proc_id+file_num_procs; proc++ ) {
				total_cells[level] += proc_num_cells[(max_level-min_level+1)*proc+level];
			}
		}
	}

	/* only write one copy of the header information */
	if ( local_proc_id == MASTER_NODE ) {
		size = 256*sizeof(char);
		fwrite(&size, sizeof(int), 1, output );
		fwrite(&jobname, sizeof(char), 256, output );
		fwrite(&size, sizeof(int), 1, output );

		/* istep, t, dt, adum, ainit */
		adum = aexp[min_level];
		ainit = a_init;
		size = sizeof(int) + 2*sizeof(double) + 2*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &step, sizeof(int), 1, output );
		fwrite( &tl[min_level], sizeof(double), 1, output );
		fwrite( &dtl[min_level], sizeof(double), 1, output );
		fwrite( &adum, sizeof(float), 1, output );
		fwrite( &ainit, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* boxh, Om0, Oml0, Omb0, hubble */
		boxh = Lbox;
		Om0 = Omega0;
		Oml0 = OmegaL0;
		Omb0 = Omegab0;
		h = hubble;
		size = 5*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &boxh, sizeof(float), 1, output );
		fwrite( &Om0, sizeof(float), 1, output );
		fwrite( &Oml0, sizeof(float), 1, output );
		fwrite( &Omb0, sizeof(float), 1, output );
		fwrite( &h, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* nextra (no evidence extras are used...) extra lextra */
		size = sizeof(int);
		nextras = 0;

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &nextras, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* extra */
		size = nextras * sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* lextra */
		size = nextras * 256 * sizeof(char);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );
	
		/* Minlevel, MaxLevelNow */
		size = 2 * sizeof(int);
		fwrite(&size, sizeof(int), 1, output );
		fwrite(&minlevel, sizeof(int), 1, output );
		fwrite(&maxlevel, sizeof(int), 1, output );
		fwrite(&size, sizeof(int), 1, output );

		size = (maxlevel-minlevel+1) * sizeof(double);

		/* tl */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &tl, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* dtl */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &dtl, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* tl_old */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &tl_old, sizeof(double), maxlevel-minlevel+1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* dtl_old */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &dtl_old, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* iSO */
		size = (maxlevel-minlevel+1) * sizeof(int);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &level_sweep_dir, sizeof(int), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* sfc ordering used */
		sfc_order = SFC;
		size = sizeof(int);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &sfc_order, sizeof(int), 1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* refinement volume */
		size = 2*nDim*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( refinement_volume_min, sizeof(float), nDim, output );
		fwrite( refinement_volume_max, sizeof(float), nDim, output );
		fwrite( &size, sizeof(int), 1, output );

#ifdef STARFORM
		/* refinement volume */
		size = 2*nDim*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( star_formation_volume_min, sizeof(float), nDim, output );
		fwrite( star_formation_volume_max, sizeof(float), nDim, output );
		fwrite( &size, sizeof(int), 1, output );
#endif /* STARFORM */

		/* ncell0 */
		ncell0 = num_grid*num_grid*num_grid;
		size = sizeof(int);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &ncell0, sizeof(int), 1, output);
		fwrite( &size, sizeof(int), 1, output );
	}

	/* now start writing pages of root level cell children */
	if ( local_proc_id == file_parent ) {
		size = total_cells[min_level] * sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
	}

	/* holds list of next level octs to write */
	order = cart_alloc( (num_cells_per_level[min_level+1] / num_children) * sizeof(int) );
	next_level_count = 0;
	current_level_count = 0;
	page = 0;

	while ( current_level_count < num_cells_per_level[min_level] ) {
		page_count = min( page_size, num_cells_per_level[min_level] - current_level_count );

		for ( i = 0; i < page_count; i++ ) {
			if ( cell_is_refined(current_level_count) ) {
				cellrefined[i] = tree_cell_count(current_level_count);
				order[next_level_count++] = cell_child_oct[current_level_count];
			} else {
				cellrefined[i] = 1;
			}
			current_level_count++;
		}

		if ( local_proc_id == file_parent ) {
			fwrite( cellrefined, sizeof(int), page_count, output );
		} else {
			MPI_Send( cellrefined, page_count, MPI_INT, file_parent, page, MPI_COMM_WORLD );
			page++;
		}
	}

	if ( local_proc_id == file_parent ) {
		for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
			i = 0;
			page = 0;
			while ( i < proc_num_cells[(max_level-min_level+1)*proc+min_level] ) {
				page_count = min( page_size,  proc_num_cells[(max_level-min_level+1)*proc+min_level] - i );
				MPI_Recv( cellrefined, page_count, MPI_INT, proc, page, MPI_COMM_WORLD, &status );
				fwrite( cellrefined, sizeof(int), page_count, output );
				i += page_count;
				page++;
			}
		}

		fwrite( &size, sizeof(int), 1, output );
	}

#ifdef HYDRO
	/* now write pages of root level hydro variables */
	if ( local_proc_id == file_parent ) {
		size = num_hydro_vars * total_cells[min_level] * sizeof(float);
		fwrite( &size, sizeof(int), 1, output );
	}

	current_level_count = 0;
	page = 0;
	while ( current_level_count < num_cells_per_level[min_level] ) {
		page_count = min( page_size, num_cells_per_level[min_level] - current_level_count );

		i = 0;
		while ( i < num_hydro_vars*page_count ) {
			icell = current_level_count;
			cellhvars[i++] = cell_gas_density(icell);
			cellhvars[i++] = cell_gas_energy(icell);
			cellhvars[i++] = cell_momentum(icell,0);
			cellhvars[i++] = cell_momentum(icell,1);
			cellhvars[i++] = cell_momentum(icell,2);
			cellhvars[i++] = cell_gas_pressure(icell);
			cellhvars[i++] = cell_gas_gamma(icell);
			cellhvars[i++] = cell_gas_internal_energy(icell);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			cellhvars[i++] = cell_electron_internal_energy(icell);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef ADVECT_SPECIES
			for ( j = 0; j < num_chem_species; j++ ) {
				cellhvars[i++] = cell_advected_variable(icell,j);
			}
#endif /* ADVECT_SPECIES */

			current_level_count++;
		}

		if ( local_proc_id == file_parent ) {
			fwrite( cellhvars, sizeof(float), num_hydro_vars*page_count, output );
		} else {
			MPI_Send( cellhvars, num_hydro_vars*page_count, MPI_FLOAT,
				file_parent, page, MPI_COMM_WORLD );
			page++;
		}
	}

	if ( local_proc_id == file_parent ) {
		for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
			i = 0;
			page = 0;
			while ( i <  proc_num_cells[(max_level-min_level+1)*proc+min_level] ) {
				page_count = min( page_size,  proc_num_cells[(max_level-min_level+1)*proc+min_level] - i );

				MPI_Recv( cellhvars, num_hydro_vars*page_count, MPI_FLOAT, 
					proc, page, MPI_COMM_WORLD, &status );
				fwrite( cellhvars, sizeof(float), num_hydro_vars*page_count, output );
				i += page_count;
				page++;
			}
		}

		fwrite( &size, sizeof(int), 1, output );
	}
#endif /* HYDRO */

#ifdef GRAVITY
	/* finally, write potential variables */
	if ( local_proc_id == file_parent ) {
		/* finally write potential variables */
		size = 2 * total_cells[min_level] * sizeof(float);
		fwrite( &size, sizeof(int), 1, output );
	}

	current_level_count = 0;
	page = 0;
	while ( current_level_count < num_cells_per_level[min_level] ) {
		page_count = min( page_size, num_cells_per_level[min_level] - current_level_count );
	
		i = 0;
		while ( i < 2*page_count ) {
			icell = current_level_count;
			cellvars[i++] = cell_potential(icell);
			cellvars[i++] = cell_potential_hydro(icell);
			current_level_count++;
		}

		if ( local_proc_id == file_parent ) {
			fwrite( cellvars, sizeof(float), 2*page_count, output );
		} else {
			MPI_Send( cellvars, 2*page_count, MPI_FLOAT,
				file_parent, page, MPI_COMM_WORLD );
			page++;
		}
	}

	if ( local_proc_id == file_parent ) {
		for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
			i = 0;
			page = 0;
			while ( i <  proc_num_cells[(max_level-min_level+1)*proc+min_level] ) {
				page_count = min( page_size,  proc_num_cells[(max_level-min_level+1)*proc+min_level] - i );
				MPI_Recv( cellvars, 2*page_count, MPI_FLOAT,
					proc, page, MPI_COMM_WORLD, &status );
				fwrite( cellvars, sizeof(float), 2*page_count, output );
				i += page_count;
				page++;
			}
		}

		fwrite ( &size, sizeof(int), 1, output );
	}
#endif /* GRAVITY */

	/* then write each level's cells in turn */
	for ( level = min_level+1; level <= maxlevel; level++ ) {
		if ( local_proc_id == file_parent ) {
			/* write size */
			size = sizeof(long);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &total_cells[level], sizeof(long), 1, output );
			fwrite( &size, sizeof(int), 1, output );

			size = total_cells[level] * sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
		}

		current_level_count = next_level_count;
		cart_assert( num_cells_per_level[level] == current_level_count*num_children );
		next_level_count = 0;
		current_level = order;

		if ( level < maxlevel ) {
			order = cart_alloc( ( num_cells_per_level[level+1]/num_children) * sizeof(int) );
		} else {
			order = cart_alloc( 0 );
		}

		i = 0;
		page = 0;
		while ( i < current_level_count ) {
			page_count = min( page_size, num_cells_per_level[level] - i*num_children );

			j = 0;
			while ( j < page_count) {
				ioct = current_level[i];

				for ( k = 0; k < num_children; k++ ) {
					icell = oct_child( ioct, k );

					if ( cell_is_refined(icell) ) {
						cellrefined[j++] = 1;
						order[next_level_count++] = cell_child_oct[icell];
					} else {
						cellrefined[j++] = 0;
					}
				}

				i++;
			}

			if ( local_proc_id == file_parent ) {
				fwrite( cellrefined, sizeof(int), page_count, output );
			} else {
				MPI_Send( cellrefined, page_count, MPI_INT, file_parent, page, MPI_COMM_WORLD );
				page++;
			}
		}

		if ( local_proc_id == file_parent ) {
			for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
				i = 0;
				page = 0;
				while ( i < proc_num_cells[(max_level-min_level+1)*proc+level] ) {
					page_count = min( page_size, 
							proc_num_cells[(max_level-min_level+1)*proc+level] - i );
					MPI_Recv( cellrefined, page_count, MPI_INT, proc, page, MPI_COMM_WORLD, &status );
					fwrite( cellrefined, sizeof(int), page_count, output );
					i += page_count;
					page++;
				}
			}

			fwrite( &size, sizeof(int), 1, output );
		}

#ifdef HYDRO
		/* now write pages of root level hydro variables */
		if ( local_proc_id == file_parent ) {
			size = num_hydro_vars * total_cells[level] * sizeof(float);
			fwrite( &size, sizeof(int), 1, output );
		}

		i = 0;
		page = 0;
		while ( i < current_level_count ) {
			page_count = min( page_size, num_cells_per_level[level] - i*num_children );

			j = 0;
			while ( j < num_hydro_vars*page_count ) { 
				ioct = current_level[i];

				for ( k = 0; k < num_children; k++ ) {
					icell = oct_child( ioct, k );

					cellhvars[j++] = cell_gas_density(icell);
					cellhvars[j++] = cell_gas_energy(icell);
					cellhvars[j++] = cell_momentum(icell,0);
					cellhvars[j++] = cell_momentum(icell,1);
					cellhvars[j++] = cell_momentum(icell,2);
					cellhvars[j++] = cell_gas_pressure(icell);
					cellhvars[j++] = cell_gas_gamma(icell);
					cellhvars[j++] = cell_gas_internal_energy(icell);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
					cellhvars[j++] = cell_electron_internal_energy(icell);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef ADVECT_SPECIES
					for ( m = 0; m < num_chem_species; m++ ) {
						cellhvars[j++] = cell_advected_variable(icell,m);
					}
#endif /* ADVECT_SPECIES */

				}
				i++;
			}

			if ( local_proc_id == file_parent ) {
				fwrite( cellhvars, sizeof(float), num_hydro_vars*page_count, output );
			} else {
				MPI_Send( cellhvars, num_hydro_vars*page_count, 
						MPI_FLOAT, file_parent, page, MPI_COMM_WORLD );
				page++;
			}
		}

		if ( local_proc_id == file_parent ) {
			for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
				i = 0;
				page = 0;
				while ( i < proc_num_cells[(max_level-min_level+1)*proc+level] ) {
					page_count = min( page_size, 
						proc_num_cells[(max_level-min_level+1)*proc+level] - i );
					MPI_Recv( cellhvars, num_hydro_vars*page_count, MPI_FLOAT,
						proc, page, MPI_COMM_WORLD, &status );
					fwrite( cellhvars, sizeof(float), num_hydro_vars*page_count, output );
					i += page_count;
					page++;
				}
			}

			fwrite( &size, sizeof(int), 1, output );
		}
#endif /* HYDRO */

#ifdef GRAVITY
		/* finally write potential variables */
		if ( local_proc_id == file_parent ) {
			size = 2 * total_cells[level] * sizeof(float);
			fwrite( &size, sizeof(int), 1, output );
		}

		i = 0;
		page = 0;
		while ( i < current_level_count ) {
			page_count = min( page_size, num_cells_per_level[level] - i*num_children );

			j = 0;
			while ( j < 2*page_count ) {
				ioct = current_level[i];

				for ( k = 0; k < num_children; k++ ) {
					icell = oct_child( ioct, k );

					cellvars[j++] = cell_potential(icell);
					cellvars[j++] = cell_potential_hydro(icell);
				}

				i++;
			}

			if ( local_proc_id == file_parent ) {
				fwrite( cellvars, sizeof(float), 2*page_count, output );
			} else {
				MPI_Send( cellvars, 2*page_count, MPI_FLOAT, 
						file_parent, page, MPI_COMM_WORLD );
				page++;
			}
		}

		if ( local_proc_id == file_parent ) {
			for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
				i = 0;
				page = 0;
				while ( i < proc_num_cells[(max_level-min_level+1)*proc+level] ) {
					page_count = min( page_size, 
						proc_num_cells[(max_level-min_level+1)*proc+level] - i );
					MPI_Recv( cellvars, 2*page_count, MPI_FLOAT,
						proc, page, MPI_COMM_WORLD, &status );
					fwrite( cellvars, sizeof(float), 2*page_count, output );
					i += page_count;
					page++;
				}
			}

			fwrite ( &size, sizeof(int), 1, output );
		}
#endif /* GRAVITY */

		cart_free( current_level );
	}

	if ( local_proc_id == file_parent ) {
		fclose( output );
	}

	cart_free( order );
	cart_free( cellrefined );
	cart_free( cellhvars );
	cart_free( cellvars );
}

void read_grid_binary( char *filename ) {
	int i, j, k, m;
	int size;
	int num_read, flag;
	FILE *input;
	char job[256];
	int minlevel, maxlevel;
	double t, dt;
	float adum, ainit;
        float boxh, Om0, Oml0, Omb0, h;
	int nextras;
	float extra[10];
	char lextra[10][256];
	int *cellrefined[MAX_PROCS], *cellrefinedbuffer;
	float *cellhvars[MAX_PROCS], *cellhvars_buffer;
	float *cellvars[MAX_PROCS], *cellvars_buffer;
	int *order;
	int *current_level;
	int endian;
	int proc;
	int ret;
	int level;
	int icell, ioct;
	int cell_counts;
	int page_size;
	int ncell0, sfc_order;
	long current_level_count, current_read_count;
	long next_level_count;
	int file_index, file_parent;
	int local_file_root_cells;
	int page_count;
	long count, start;
	long total_cellrefined;
	long total_cells[max_level-min_level+1];
	int file_parent_proc[MAX_PROCS];
	int proc_num_cells[MAX_PROCS];
	int proc_cur_cells[MAX_PROCS];
	int proc_page_count[MAX_PROCS];
	long first_page_count[MAX_PROCS];
	long proc_next_level_octs[MAX_PROCS];
	long proc_cell_index[MAX_PROCS];
	long proc_first_index[MAX_PROCS+1];
	long next_proc_index[MAX_PROCS+1];
	int file_sfc_index[MAX_PROCS+1];
	char parallel_filename[256];
	int num_requests, num_send_requests;
	int continue_reading, ready_to_read;
	MPI_Request requests[2*MAX_PROCS];
	MPI_Request send_requests[MAX_PROCS];
	MPI_Status status;

	cart_assert( num_output_files >= 1 && num_output_files <= num_procs );

	page_size = num_grid*num_grid;

	/* set up global file information */
	proc = 0;	
	for ( i = 0; i < num_output_files; i++ ) {
		while ( proc*num_output_files / num_procs != i ) {
			proc++;
		}
		file_parent_proc[i] = proc;
	}

	/* determine parallel output options */
	file_index = local_proc_id * num_output_files / num_procs;
	file_parent = file_parent_proc[file_index];

	/* open file handle if parent of parallel file */
	if ( local_proc_id == file_parent ) {
		if ( num_output_files == 1 ) {
			input = fopen(filename,"r");
		} else {
			sprintf( parallel_filename, "%s.%03u", filename, file_index );
			input = fopen(parallel_filename, "r");
		}

		if ( input == NULL ) {
			cart_error( "Unable to open file %s for reading!", filename );
		}
	}

	/* the header exists only in the first file */
	if ( local_proc_id == MASTER_NODE ) {
                fread(&size, sizeof(int), 1, input );
		endian = 0;
		if ( size != 256 ) {
			reorder( (char *)&size, sizeof(int) );
			if ( size != 256 ) {
				cart_error("Error: file %s is corrupted", filename );
			} else {
				endian = 1;
				cart_debug("Reordering bytes (file endianness is opposite program)");
			}
		}

                fread(&job, sizeof(char), 256, input );
                fread(&size, sizeof(int), 1, input );

		/* istep, t, dt, adum, ainit */
		fread( &size, sizeof(int), 1, input );
		fread( &step, sizeof(int), 1, input );
		fread( &t, sizeof(double), 1, input );
		fread( &dt, sizeof(double), 1, input );
		fread( &adum, sizeof(float), 1, input );
		fread( &ainit, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&step, sizeof(int) );
			reorder( (char *)&ainit, sizeof(float) );
		}

		a_init = ainit;

#ifndef COSMOLOGY
		if ( endian ) {
			reorder( (char *)&adum, sizeof(float) );
		}
		for(i=min_level; i<=max_level; i++) aexp[i] = adum;
#endif

                /* boxh, Om0, Oml0, Omb0, hubble */
		fread( &size, sizeof(int), 1, input );
		fread( &boxh, sizeof(float), 1, input );
		fread( &Om0, sizeof(float), 1, input );
		fread( &Oml0, sizeof(float), 1, input );
		fread( &Omb0, sizeof(float), 1, input );
		fread( &h, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&boxh, sizeof(float) );
			reorder( (char *)&Om0, sizeof(float) );
			reorder( (char *)&Oml0, sizeof(float) );
			reorder( (char *)&Omb0, sizeof(float) );
			reorder( (char *)&h, sizeof(float) );
		}

		Lbox = boxh;
		Omega0 = Om0;
		OmegaL0 = Oml0;
		Omegab0 = Omb0;
		hubble = h;

		/* nextra (no evidence extras are used...) extra lextra */
		fread( &size, sizeof(int), 1, input );
		fread( &nextras, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&nextras, sizeof(int) );
		}

		/* extra */
		fread( &size, sizeof(int), 1, input );
		fread( extra, sizeof(float), nextras, input );
		fread( &size, sizeof(int), 1, input );

		/* lextra */
		fread( &size, sizeof(int), 1, input );
		fread( lextra, 256*sizeof(char), nextras, input );
		fread( &size, sizeof(int), 1, input );

		/* Minlevel, MaxLevelNow */
		fread( &size, sizeof(int), 1, input );
		fread( &minlevel, sizeof(int), 1, input );
		fread( &maxlevel, sizeof(int), 1, input);

		if ( endian ) {
			reorder( (char *)&minlevel, sizeof(int) );
			reorder( (char *)&maxlevel, sizeof(int) );
		}

		if ( maxlevel > max_level ) {
			cart_error("File %s has more levels than compiled program (%u)", filename, maxlevel );
		}

		cart_assert( minlevel == min_level );

		fread( &size, sizeof(int), 1, input );

		/* tl */
		fread( &size, sizeof(int), 1, input );
		fread( &tl, sizeof(double), maxlevel-minlevel+1, input );
		fread( &size, sizeof(int), 1, input);

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&tl[i], sizeof(double) );
			}
		}

		/* dtl */
		fread( &size, sizeof(int), 1, input );
		fread( &dtl, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&dtl[i], sizeof(double) );
			}
		}

		/* tl_old */
		fread( &size, sizeof(int), 1, input );
		fread( &tl_old, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&tl_old[i], sizeof(double) );
			}
		}

		/* dtl_old */
		fread( &size, sizeof(int), 1, input );
		fread( &dtl_old, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );	

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&dtl_old[i], sizeof(double) );
			}
		}

		/* iSO */
		fread( &size, sizeof(int), 1, input );
		fread( &level_sweep_dir, sizeof(int), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&level_sweep_dir[i], sizeof(int) );
			}
		}

		/* sfc ordering used */
		fread( &size, sizeof(int), 1, input );
		fread( &sfc_order, sizeof(int), 1, input);

		if ( endian ) {
			reorder( (char *)&sfc_order, sizeof(int) );
		}

		if ( sfc_order != SFC ) {
			cart_error("File has different sfc indexing than program");
		}
		fread( &size, sizeof(int), 1, input );

		/* refinement volume */
		fread( &size, sizeof(int), 1, input );
		fread( refinement_volume_min, sizeof(float), nDim, input );
		fread( refinement_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&refinement_volume_min[i], sizeof(float) );
				reorder( (char *)&refinement_volume_max[i], sizeof(float) );
			}
		}

#ifdef STARFORM
		/* star formation volume */
		fread( &size, sizeof(int), 1, input );
		fread( star_formation_volume_min, sizeof(float), nDim, input );
		fread( star_formation_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&star_formation_volume_min[i], sizeof(float) );
				reorder( (char *)&star_formation_volume_max[i], sizeof(float) );
			}
		}
#endif /* STARFORM */

		/* ncell0 */
		fread( &size, sizeof(int), 1, input );
		fread( &ncell0, sizeof(int), 1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&ncell0, sizeof(int) );
		}

		if ( ncell0 != (num_grid*num_grid*num_grid) ) {
			cart_error("File has different num_grid than compiled program");
		}
	}

	/* send header information to all other processors */
	MPI_Bcast( &endian, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &minlevel, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &maxlevel, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &step, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &a_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Lbox, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Omega0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &OmegaL0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Omegab0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &hubble, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( tl, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( dtl, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( dtl_old, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( level_sweep_dir, max_level-min_level+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	MPI_Bcast( refinement_volume_min, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( refinement_volume_max, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );

#ifdef STARFORM
	MPI_Bcast( star_formation_volume_min, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( star_formation_volume_max, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
#endif /* STARFORM */

	if ( local_proc_id == file_parent ) {
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&size, sizeof(int) );
		}
		
		local_file_root_cells = size / sizeof(int);
		cellrefinedbuffer = cart_alloc( local_file_root_cells * sizeof(int) );

		num_read = fread( cellrefinedbuffer, sizeof(int), local_file_root_cells, input );

		if ( num_read != local_file_root_cells ) {
			cart_error("I/O error in read_grid_binary: num_read = %d, local_file_root_cells = %d",
				num_read, local_file_root_cells );
		}

		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < local_file_root_cells; i++ ) {
				reorder( (char *)&cellrefinedbuffer[i], sizeof(int) );
			}
		}

		if ( local_proc_id == MASTER_NODE ) {
			cell_counts = local_file_root_cells;

			file_sfc_index[0] = 0;
			for ( i = 1; i < num_output_files; i++ ) {
				/* need to know how many root cells are in each file */
				/* receive cellrefined */
				file_sfc_index[i] = cell_counts;
				MPI_Recv( &file_sfc_index[i+1], 1, MPI_INT, file_parent_proc[i],
					0, MPI_COMM_WORLD, &status );	
				cell_counts += file_sfc_index[i+1];
			}
			file_sfc_index[num_output_files] = cell_counts;

			cart_assert( cell_counts == num_root_cells );
		} else {
			/* send cellrefined array to MASTER_NODE */
			MPI_Send( &local_file_root_cells, 1, MPI_INT, MASTER_NODE, 0, MPI_COMM_WORLD );
		}
	}

	/* send block information */
	MPI_Bcast( file_sfc_index, num_output_files+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	/* determine how many cells to expect from each processor */
	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_num_cells[proc] = 0;
	}

	for ( i = 0; i < num_output_files; i++ ) {
		if ( proc_sfc_index[local_proc_id] < file_sfc_index[i+1] &&
				proc_sfc_index[local_proc_id+1] >= file_sfc_index[i] ) {
			proc_num_cells[ file_parent_proc[i] ] = 
					min( proc_sfc_index[local_proc_id+1], file_sfc_index[i+1] ) - 
					max( proc_sfc_index[local_proc_id], file_sfc_index[i] );
		}
	}

	order = cart_alloc( num_cells_per_level[min_level] * sizeof(int) );
	cellrefined[local_proc_id] = cart_alloc( num_cells_per_level[min_level] * sizeof(int) );

	num_requests = 0;

	/* set up non-blocking receives */
	count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( proc_num_cells[proc] > 0 ) { 
			if ( proc == local_proc_id ) {
				/* copy data directly */
				start = max( 0, proc_sfc_index[local_proc_id] - file_sfc_index[file_index] );
				cart_assert( start >= 0 && start <= local_file_root_cells - proc_num_cells[local_proc_id] );
				for ( i = 0; i < proc_num_cells[local_proc_id]; i++ ) {
					cellrefined[local_proc_id][count+i] = cellrefinedbuffer[start+i];
				}
			} else {
				MPI_Irecv( &cellrefined[local_proc_id][count], proc_num_cells[proc], MPI_INT, proc,
					proc_num_cells[proc], MPI_COMM_WORLD, &requests[num_requests++] );
			}

			count += proc_num_cells[proc];
		}
	}

	cart_assert( count == num_cells_per_level[min_level] );

	/* need to avoid deadlocks! senders may also be receivers! */
	if ( local_proc_id == file_parent ) {
		start = 0;
		next_proc_index[proc] = 0;

		/* send all necessary information */
		for ( proc = 0; proc < num_procs; proc++ ) {
			next_proc_index[proc+1] = 0;

			if ( proc == local_proc_id ) { 
				for ( i = 0; i < proc_num_cells[local_proc_id]; i++ ) {
					if ( cellrefinedbuffer[start+i] > 1 ) {
						next_proc_index[proc+1]++;
					}
				}
				start += proc_num_cells[local_proc_id];
			} else if ( start < local_file_root_cells && proc != local_proc_id && 
					file_sfc_index[file_index]+start < proc_sfc_index[proc+1] &&
					file_sfc_index[file_index]+start >= proc_sfc_index[proc] ) {
				count = min( proc_sfc_index[proc+1], file_sfc_index[file_index+1] ) -
						max( file_sfc_index[file_index] + start, proc_sfc_index[proc] ); 
				cart_assert( count > 0 && start + count <= local_file_root_cells );
				for ( i = 0; i < count; i++ ) {
					if ( cellrefinedbuffer[start+i] > 1 ) {
						next_proc_index[proc+1]++;
					}
				}
				MPI_Isend( &cellrefinedbuffer[start], count, MPI_INT, proc, 
					count, MPI_COMM_WORLD, &requests[num_requests++] );
				start += count;
			}
		}

		cart_assert( start == local_file_root_cells );
	}

	/* wait for sends/receives to complete */
	MPI_Waitall( num_requests, requests, MPI_STATUSES_IGNORE );

	if ( local_proc_id == file_parent ) {
		cart_free( cellrefinedbuffer );
	}

	current_level_count = 0;
	next_level_count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_cell_index[proc] = proc_sfc_index[local_proc_id] + current_level_count;
		proc_next_level_octs[proc] = 0;
	
		for ( i = 0; i < proc_num_cells[proc]; i++ ) {
			if ( cellrefined[local_proc_id][current_level_count] > 1 ) {
				ret = split_cell(current_level_count);
				if ( ret ) {
					cart_error("Unable to finish splitting root cells, ran out of octs?");
				}
				cart_assert( next_level_count < num_cells_per_level[min_level] );
				order[next_level_count++] = cell_child_oct[current_level_count];
				proc_next_level_octs[proc]++;
			}
			current_level_count++;
		}
	}

	cart_assert( next_level_count <= num_cells_per_level[min_level] );

	cart_free( cellrefined[local_proc_id] );

#ifdef HYDRO
	if ( local_proc_id == file_parent ) {
		fread( &size, sizeof(int), 1, input );

		cellhvars_buffer = cart_alloc( 
				num_hydro_vars*min( local_file_root_cells, page_size ) * sizeof(float) );
	}

	num_requests = 0;

	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_cur_cells[proc] = 0;

		if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
			/* set up receive */
			proc_page_count[proc] = min( proc_num_cells[proc], page_size - ( proc_cell_index[proc] - 
				file_sfc_index[(int)((proc * num_output_files)/ num_procs)]) % page_size);
			cellhvars[proc] = cart_alloc( num_hydro_vars*min( page_size, proc_num_cells[proc] )*sizeof(float) );
			MPI_Irecv( cellhvars[proc], num_hydro_vars*proc_page_count[proc], MPI_FLOAT, proc, 
				proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
			num_requests++;
		} else {
			proc_page_count[proc] = 0;
			requests[proc] = MPI_REQUEST_NULL;
		}
	}


	current_level_count = 0;
	current_read_count = 0;
	continue_reading = 1;
	ready_to_read = 1;

	while ( continue_reading ) {
		flag = 0;

		if ( local_proc_id == file_parent && current_read_count < local_file_root_cells && ready_to_read ) {
			/* read in a page */
			page_count = min( page_size, local_file_root_cells - current_read_count );
			num_read = fread( cellhvars_buffer, sizeof(float), num_hydro_vars*page_count, input );

			if ( num_read != num_hydro_vars*page_count ) {
				cart_error("I/O Error in read_grid_binary: num_read = %u, num_hydro_vars*page_count = %u",
					num_read, num_hydro_vars*page_count );
			}

			if ( endian ) {
				for ( i = 0; i < num_hydro_vars*page_count; i++ ) {
					reorder( (char *)&cellhvars_buffer[i], sizeof(float) );
				}
			}

			/* send info to other processors */
			start = 0;
			num_send_requests = 0;
			for ( proc = 0; proc < num_procs; proc++ ) {
				if ( start < page_count && 
					file_sfc_index[file_index]+current_read_count+start < proc_sfc_index[proc+1] &&
					file_sfc_index[file_index]+current_read_count+start >= proc_sfc_index[proc] ) {

					count = min( proc_sfc_index[proc+1], 
							file_sfc_index[file_index]+current_read_count+page_count ) -
						max( proc_sfc_index[proc], 
							file_sfc_index[file_index] + current_read_count + start );
					cart_assert( count > 0 && count <= page_count );

					if ( proc == local_proc_id ) {
						/* copy into local buffer */
						proc_page_count[local_proc_id] = count;

						cellhvars[local_proc_id] = 
							cart_alloc( num_hydro_vars*
							proc_page_count[local_proc_id]*sizeof(float));

						for ( i = 0; i < num_hydro_vars*proc_page_count[local_proc_id]; i++ ) {
							cart_assert( num_hydro_vars*start + i < num_hydro_vars*page_count );
							cellhvars[local_proc_id][i] = 
									cellhvars_buffer[num_hydro_vars*start+i];
						}
					} else {
						MPI_Isend( &cellhvars_buffer[num_hydro_vars*start], 
							num_hydro_vars*count, MPI_FLOAT, proc,
							count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
					}

					start += count;
				}
			}

			cart_assert( start == page_count );
			current_read_count += page_count;
			ready_to_read = 0;

			/* unpack local cells first */
			if ( proc_page_count[local_proc_id] > 0 ) {
				flag = 1;
				proc = local_proc_id;
			} else if ( num_requests > 0 ) {
				/* see if we've received anything */
				MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
			}
		} else if ( num_requests > 0 ) {
			/* wait for a receive to complete */
			MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
			num_requests--;
			flag = 1;
		}

		if ( flag == 1 && proc != MPI_UNDEFINED ) {
			cart_assert( proc >= 0 && proc < num_procs );
			
			/* unpack received page */
			i = 0;
			while ( i < num_hydro_vars*proc_page_count[proc] ) {
				icell = root_cell_location( proc_cell_index[proc] + proc_cur_cells[proc]);
				cart_assert( icell >= 0 && icell < num_cells_per_level[min_level] );

				cell_gas_density(icell) = cellhvars[proc][i++];
				cell_gas_energy(icell) = cellhvars[proc][i++];
				cell_momentum(icell,0) = cellhvars[proc][i++];
				cell_momentum(icell,1) = cellhvars[proc][i++];
				cell_momentum(icell,2) = cellhvars[proc][i++];
				cell_gas_pressure(icell) = cellhvars[proc][i++];
				cell_gas_gamma(icell) = cellhvars[proc][i++];
				cell_gas_internal_energy(icell) = cellhvars[proc][i++];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
				cell_electron_internal_energy(icell) = cellhvars[proc][i++];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef ADVECT_SPECIES
				for ( m = 0; m < num_chem_species; m++ ) {
					cell_advected_variable(icell,m) = cellhvars[proc][i++];
				}
#endif /* ADVECT_SPECIES */

				current_level_count++;
				proc_cur_cells[proc]++;
			}

			/* if necessary, start a new receive for that proc */
			if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id ) {
				proc_page_count[proc] = min( page_size, proc_num_cells[proc] - proc_cur_cells[proc] );
				cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
				MPI_Irecv( cellhvars[proc], num_hydro_vars*proc_page_count[proc], MPI_FLOAT, 
					proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
				num_requests++;
			} else {
				proc_page_count[proc] = 0;
				cart_free( cellhvars[proc] );
				requests[proc] = MPI_REQUEST_NULL;
			}
		}

		if ( local_proc_id == file_parent ) {
			if ( num_requests > 0 ) {
				/* see if sends have completed, if they have we can continue to read */
				MPI_Testall( num_send_requests, send_requests, &ready_to_read, MPI_STATUSES_IGNORE );
			} else {
				/* not waiting on any receives, so we can't do anything until
				 * we read in additional data */
				MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
				num_send_requests = 0;
				ready_to_read = 1;
			}

			if ( current_read_count >= local_file_root_cells && 
					current_level_count >= num_cells_per_level[min_level] ) {
				continue_reading = 0;
			}
		} else {
			if ( current_level_count >= num_cells_per_level[min_level] ) {
				continue_reading = 0;
			}
		}
	}

	if ( local_proc_id == file_parent ) {
		if ( !ready_to_read ) {
			MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
		}

		cart_free( cellhvars_buffer );
		fread( &size, sizeof(int), 1, input );
	}
#endif /* HYDRO */

#ifdef GRAVITY
	if ( local_proc_id == file_parent ) {
		fread( &size, sizeof(int), 1, input );

		cellvars_buffer = cart_alloc( 2*min( local_file_root_cells, page_size ) * sizeof(float) );
	}

	num_requests = 0;

	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_cur_cells[proc] = 0;

		if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
			/* set up receive */
			proc_page_count[proc] = min( proc_num_cells[proc], page_size - ( proc_cell_index[proc] - 
				file_sfc_index[(int)((proc * num_output_files)/ num_procs)]) % page_size );
			cellvars[proc] = cart_alloc( 2*min( page_size, proc_num_cells[proc] )*sizeof(float) );
			MPI_Irecv( cellvars[proc], 2*proc_page_count[proc], MPI_FLOAT, proc, 
				proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
			num_requests++;
		} else {
			proc_page_count[proc] = 0;
			requests[proc] = MPI_REQUEST_NULL;
		}
	}

	current_level_count = 0;
	current_read_count = 0;
	continue_reading = 1;
	ready_to_read = 1;

	while ( continue_reading ) {
		flag = 0;

		if ( local_proc_id == file_parent && current_read_count < local_file_root_cells && ready_to_read ) {
			/* read in a page */
			page_count = min( page_size, local_file_root_cells - current_read_count );
			num_read = fread( cellvars_buffer, sizeof(float), 2*page_count, input );

			if ( num_read != 2*page_count ) {
				cart_error("I/O Error in read_grid_binary: num_read = %u, 2*page_count = %u",
					num_read, 2*page_count );
			}

			if ( endian ) {
				for ( i = 0; i < 2*page_count; i++ ) {
					reorder( (char *)&cellvars_buffer[i], sizeof(float) );
				}
			}

			/* send info to other processors */
			start = 0;
			num_send_requests = 0;
			for ( proc = 0; proc < num_procs; proc++ ) {
				if ( start < page_count && 
					file_sfc_index[file_index]+current_read_count+start < proc_sfc_index[proc+1] &&
					file_sfc_index[file_index]+current_read_count+start >= proc_sfc_index[proc] ) {

					count = min( proc_sfc_index[proc+1], 
							file_sfc_index[file_index]+current_read_count+page_count ) -
						max( proc_sfc_index[proc], 
							file_sfc_index[file_index] + current_read_count + start );
					cart_assert( count > 0 && count <= page_count );

					if ( proc == local_proc_id ) {
						/* copy into local buffer */
						proc_page_count[local_proc_id] = count;

						cellvars[local_proc_id] = 
							cart_alloc( 2*proc_page_count[local_proc_id]*sizeof(float));

						for ( i = 0; i < 2*proc_page_count[local_proc_id]; i++ ) {
							cart_assert( 2*start + i < 2*page_count );
							cellvars[local_proc_id][i] = cellvars_buffer[2*start+i];
						}
					} else {
						MPI_Isend( &cellvars_buffer[2*start], 2*count, MPI_FLOAT, proc,
							count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
					}

					start += count;
				}
			}

			cart_assert( start == page_count );
			current_read_count += page_count;
			ready_to_read = 0;

			/* unpack local cells first */
			if ( proc_page_count[local_proc_id] > 0 ) {
				flag = 1;
				proc = local_proc_id;
			} else if ( num_requests > 0 ) {
				/* see if we've received anything */
				MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
			}
		} else if ( num_requests > 0 ) {
			/* wait for a receive to complete */
			MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
			num_requests--;
			flag = 1;
		}

		if ( flag == 1 && proc != MPI_UNDEFINED ) {
			cart_assert( proc >= 0 && proc < num_procs );

			/* unpack received page */
			i = 0;

			while ( i < 2*proc_page_count[proc] ) {
				icell = root_cell_location( proc_cell_index[proc] + proc_cur_cells[proc]);
				cart_assert( icell >= 0 && icell < num_cells_per_level[min_level] );

				cell_potential(icell) = cellvars[proc][i++];
				cell_potential_hydro(icell) = cellvars[proc][i++];

				current_level_count++;
				proc_cur_cells[proc]++;
			}

			/* if necessary, start a new receive for that proc */
			if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id ) {
				proc_page_count[proc] = min( page_size, proc_num_cells[proc] - proc_cur_cells[proc] );
				cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
				MPI_Irecv( cellvars[proc], 2*proc_page_count[proc], MPI_FLOAT, 
					proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
				num_requests++;
			} else {
				proc_page_count[proc] = 0;
				cart_free( cellvars[proc] );
				requests[proc] = MPI_REQUEST_NULL;
			}

			flag = 0;
		}

		if ( local_proc_id == file_parent ) {
			if ( num_requests > 0 ) {
				/* see if sends have completed */
				MPI_Testall( num_send_requests, send_requests, &ready_to_read, MPI_STATUSES_IGNORE );
			} else {
				MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
				num_send_requests = 0;
				ready_to_read = 1;
			}

			if ( current_read_count >= local_file_root_cells && 
					current_level_count >= num_cells_per_level[min_level] ) {
				continue_reading = 0;
			}
		} else {
			if ( current_level_count >= num_cells_per_level[min_level] ) {
				continue_reading = 0;
			}
		}
	}

	if ( local_proc_id == file_parent ) {
		if ( !ready_to_read ) {
			MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
		}

		cart_free( cellvars_buffer );
		fread( &size, sizeof(int), 1, input );
	}
#endif /* GRAVITY */

	/* now read levels */
	for ( level = min_level+1; level <= maxlevel; level++ ) {
		num_requests = 0;

		count = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			proc_num_cells[proc] = proc_next_level_octs[proc]*num_children;

			if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
				MPI_Irecv( &first_page_count[proc], 1, MPI_LONG, proc,
						0, MPI_COMM_WORLD, &requests[num_requests++] );
			} else {
				first_page_count[proc] = 0;
			}

			proc_cell_index[proc] = count;
			count += proc_next_level_octs[proc];
			proc_next_level_octs[proc] = 0;
		}

		if ( local_proc_id == file_parent ) {
			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
			}

			if ( size == sizeof(int) ) {
				fread( &size, sizeof(int), 1, input );

				if ( endian ) {
					reorder( (char *)&size, sizeof(int) );
				}

				total_cells[level] = (long)size;
			} else {
				fread( &total_cells[level], sizeof(long), 1, input );
				if ( endian ) {
					reorder( (char *)&total_cells[level], sizeof(long) );
				}
			}
			fread( &size, sizeof(int), 1, input );
			fread( &size, sizeof(int), 1, input );

			cart_debug("total_cells[%u] = %d", level, total_cells[level] );

			proc_first_index[0] = 0;
			for ( proc = 1; proc <= num_procs; proc++ ) {
				proc_first_index[proc] = (long)num_children*next_proc_index[proc] + 
								proc_first_index[proc-1];
				next_proc_index[proc] = 0;
			}

			cart_assert( proc_first_index[num_procs] == total_cells[level] );

			for ( proc = 0; proc < num_procs; proc++ ) {
				if ( proc_first_index[proc+1]-proc_first_index[proc] > 0 && proc != local_proc_id ) {
					MPI_Isend( &proc_first_index[proc], 1, MPI_LONG, proc,
						0, MPI_COMM_WORLD, &requests[num_requests++] );
				}
			}

			cellrefinedbuffer = cart_alloc( min( total_cells[level], page_size ) * sizeof(int) );
		}

		MPI_Waitall( num_requests, requests, MPI_STATUS_IGNORE );

		num_requests = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			proc_cur_cells[proc] = 0;

			if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
				/* set up receive */
				proc_page_count[proc] = min( page_size - first_page_count[proc] % page_size, 
								proc_num_cells[proc] );
				cellrefined[proc] = cart_alloc( min( page_size, proc_num_cells[proc] )*sizeof(float) );
				MPI_Irecv( cellrefined[proc], proc_page_count[proc], MPI_INT, proc,
					proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
				num_requests++;
			} else {
				proc_page_count[proc] = 0;
				requests[proc] = MPI_REQUEST_NULL;
			}
		}

		current_level = order;
		order = cart_alloc( num_children*next_level_count * sizeof(int) );

		current_level_count = 0;
		current_read_count = 0;
		next_level_count = 0;
		continue_reading = 1;
		ready_to_read = 1;
	
		while ( continue_reading ) {
			flag = 0;

			if ( local_proc_id == file_parent && current_read_count < total_cells[level] && ready_to_read ) {
				/* read in a page */
				page_count = min( page_size, total_cells[level] - current_read_count );
				num_read = fread( cellrefinedbuffer, sizeof(int), page_count, input );

				if ( num_read != page_count ) {
					cart_error("I/O Error in read_grid_binary: num_read = %u, page_count = %u",
						num_read, page_count );
				}

				if ( endian ) {
					for ( i = 0; i < page_count; i++ ) {
						reorder( (char *)&cellrefinedbuffer[i], sizeof(int) );
					}
				}

				/* send info to other processors */
				start = 0;
				num_send_requests = 0;

				for ( proc = 0; proc < num_procs; proc++ ) {
					if ( start < page_count && current_read_count + start >= proc_first_index[proc] &&
							current_read_count + start < proc_first_index[proc+1] ) {
						
						count = min( proc_first_index[proc+1], current_read_count+page_count ) -
							max( proc_first_index[proc], current_read_count + start );

						cart_assert( count > 0 && count <= page_count );

						for ( i = 0; i < count; i++ ) {
							cart_assert( start + i < page_count );
							if ( cellrefinedbuffer[start+i] ) {
								next_proc_index[proc+1]++;
							}
						}

						if ( proc == local_proc_id ) {
							/* copy into local buffer */
							proc_page_count[local_proc_id] = count;
							cellrefined[local_proc_id] =
								cart_alloc( proc_page_count[local_proc_id]*sizeof(int));
							for ( i = 0; i < proc_page_count[local_proc_id]; i++ ) {
								cart_assert( start + i < page_count );
								cellrefined[local_proc_id][i] = cellrefinedbuffer[start+i];
							}
						} else {
							MPI_Isend( &cellrefinedbuffer[start], count, MPI_INT, proc,
								count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
						}

						start += count;
					}
				}

				cart_assert( start == page_count );

				current_read_count += page_count;
				ready_to_read = 0;

				/* split local cells first */
				if ( proc_page_count[local_proc_id] > 0 ) {
					flag = 1;
					proc = local_proc_id;
				} else if ( num_requests > 0 ) {
					/* see if we've received anything */
					MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
				}

			} else if ( num_requests > 0 ) {
				/* wait for a receive to complete */
				MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
				num_requests--;
				flag = 1;
			}

			if ( flag == 1 && proc != MPI_UNDEFINED ) {
				cart_assert( proc >= 0 && proc < num_procs );

				/* split cells in received page */
				i = 0;
				while ( i < proc_page_count[proc] ) {
					cart_assert( proc_cell_index[proc] + proc_cur_cells[proc]/num_children 
							< num_cells_per_level[level]/num_children );
					ioct = current_level[ proc_cell_index[proc] + proc_cur_cells[proc]/num_children ]; 

					cart_assert( ioct >= 0 && ioct < num_octs );
					cart_assert( oct_level[ioct] == level );

					for ( j = 0; j < num_children; j++ ) {
						icell = oct_child( ioct, j );
						cart_assert( icell >= 0 && icell < num_cells );
						cart_assert( cell_level(icell) == level );

						if ( cellrefined[proc][i] == 1 ) {
							ret = split_cell(icell);

							if ( ret ) {
								cart_error("Unable to finish splitting cells, ran out of octs?");
							}

							order[ proc_cell_index[proc]*num_children +
								proc_next_level_octs[proc] ] = cell_child_oct[icell];
							next_level_count++;
							proc_next_level_octs[proc]++;

						}

						i++;
						proc_cur_cells[proc]++;
						current_level_count++;
					}
				}

				/* if necessary, start a new receive from that proc */
				if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id ) {
					proc_page_count[proc] = min( page_size, 
						proc_num_cells[proc] - proc_cur_cells[proc] );
					cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
					MPI_Irecv( cellrefined[proc], proc_page_count[proc], MPI_INT,
						proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
					num_requests++;
				} else {
					proc_page_count[proc] = 0;
					cart_free( cellrefined[proc] );
					requests[proc] = MPI_REQUEST_NULL;
				}

				flag = 0;
			}

			/* check if we're done reading */
			if ( local_proc_id == file_parent ) {
				if ( num_requests > 0 ) {
					/* see if sends have completed */
					MPI_Testall( num_send_requests, send_requests, 
							&ready_to_read, MPI_STATUSES_IGNORE );
				} else {
					MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
					num_send_requests = 0;
					ready_to_read = 1;
				}

				if ( current_read_count >= total_cells[level] &&
						current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			} else {
				if ( current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			}
		}

		if ( local_proc_id == file_parent ) {
			if ( !ready_to_read ) {
				MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
			}

			cart_free( cellrefinedbuffer );

			fread( &size, sizeof(int), 1, input );
		}

		/* repack the order array */
		count = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			for ( i = 0; i < proc_next_level_octs[proc]; i++ ) {
				cart_assert( proc_cell_index[proc]*num_children + i  >= count );
				cart_assert( order[ proc_cell_index[proc]*num_children + i ] >= 0 &&
						order[ proc_cell_index[proc]*num_children + i ] < num_octs );
				cart_assert( oct_level[ order[ proc_cell_index[proc]*num_children + i ] ] == level+1 );

				order[count++] = order[ proc_cell_index[proc]*num_children + i ];
			}
		}

		cart_assert( count == next_level_count );

#ifdef HYDRO
		if ( local_proc_id == file_parent ) {
			fread( &size, sizeof(int), 1, input );

			cellhvars_buffer = cart_alloc( 
					num_hydro_vars*min( total_cells[level], page_size ) * sizeof(float) );
		}
	
		num_requests = 0;

		for ( proc = 0; proc < num_procs; proc++ ) {
			proc_cur_cells[proc] = 0;

			if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
				/* set up receive */
				proc_page_count[proc] = min( page_size - first_page_count[proc] % page_size, proc_num_cells[proc] );
				cellhvars[proc] = cart_alloc( num_hydro_vars*min( page_size, proc_num_cells[proc] )*sizeof(float) );
				MPI_Irecv( cellhvars[proc], num_hydro_vars*proc_page_count[proc], MPI_FLOAT, proc, 
					proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
				num_requests++;
			} else {
				proc_page_count[proc] = 0;
				requests[proc] = MPI_REQUEST_NULL;
			}
		}
	
		current_level_count = 0;
		current_read_count = 0;
		continue_reading = 1;
		ready_to_read = 1;

		while ( continue_reading ) {
			flag = 0;

			if ( local_proc_id == file_parent && current_read_count < total_cells[level] && ready_to_read ) {
				/* read in a page */
				page_count = min( page_size, total_cells[level] - current_read_count );
				num_read = fread( cellhvars_buffer, sizeof(float), num_hydro_vars*page_count, input );
	
				if ( num_read != num_hydro_vars*page_count ) {
					cart_error("I/O Error in read_grid_binary: num_read = %u, num_hydro_vars*page_count = %u",
						num_read, num_hydro_vars*page_count );
				}
	
				if ( endian ) {
					for ( i = 0; i < num_hydro_vars*page_count; i++ ) {
						reorder( (char *)&cellhvars_buffer[i], sizeof(float) );
					}
				}
	
				/* send info to other processors */
				start = 0;
				num_send_requests = 0;
				for ( proc = 0; proc < num_procs; proc++ ) {
					if ( start < page_count && current_read_count + start >= proc_first_index[proc] &&
							current_read_count + start < proc_first_index[proc+1] ) {
				
						count = min( proc_first_index[proc+1], current_read_count+page_count ) -
							max( proc_first_index[proc], current_read_count + start );
						cart_assert( count > 0 && count <= page_count );

						if ( proc == local_proc_id ) {
							/* copy into local buffer */
							proc_page_count[local_proc_id] = count;

							cellhvars[local_proc_id] = 
								cart_alloc( num_hydro_vars*proc_page_count[local_proc_id]*sizeof(float));
	
							for ( i = 0; i < num_hydro_vars*proc_page_count[local_proc_id]; i++ ) {
								cart_assert( num_hydro_vars*start + i < num_hydro_vars*page_count );
								cellhvars[local_proc_id][i] = cellhvars_buffer[num_hydro_vars*start+i];
							}
						} else {
							MPI_Isend( &cellhvars_buffer[num_hydro_vars*start], num_hydro_vars*count, MPI_FLOAT, proc,
								count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
						}
	
						start += count;
					}
				}
	
				cart_assert( start == page_count );
				current_read_count += page_count;
				ready_to_read = 0;
	
				/* unpack local cells first */
				if ( proc_page_count[local_proc_id] > 0 ) {
					flag = 1;
					proc = local_proc_id;
				} else if ( num_requests > 0 ) {
					/* see if we've received anything */
					MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
				}
			} else if ( num_requests > 0 ) {
				/* wait for a receive to complete */
				MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
				num_requests--;
				flag = 1;
			}
	
			if ( flag == 1 && proc != MPI_UNDEFINED ) {
				cart_assert( proc >= 0 && proc < num_procs );

				/* unpack received page */
				i = 0;
	
				while ( i < num_hydro_vars*proc_page_count[proc] ) {
					ioct = current_level[ proc_cell_index[proc] + proc_cur_cells[proc]/num_children];
					cart_assert( ioct >= 0 && ioct < num_octs );

					for ( j = 0; j < num_children; j++ ) {	
						icell = oct_child( ioct, j );
						cart_assert( icell >= 0 && icell < num_cells );
						cart_assert( cell_level(icell) == level );

						cell_gas_density(icell) = cellhvars[proc][i++];
						cell_gas_energy(icell) = cellhvars[proc][i++];
						cell_momentum(icell,0) = cellhvars[proc][i++];
						cell_momentum(icell,1) = cellhvars[proc][i++];
						cell_momentum(icell,2) = cellhvars[proc][i++];
						cell_gas_pressure(icell) = cellhvars[proc][i++];
						cell_gas_gamma(icell) = cellhvars[proc][i++];
						cell_gas_internal_energy(icell) = cellhvars[proc][i++];

#ifdef ELECTRON_ION_NONEQUILIBRIUM
						cell_electron_internal_energy(icell) = cellhvars[proc][i++];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
	
#ifdef 	ADVECT_SPECIES
						for ( m = 0; m < num_chem_species; m++ ) {
							cell_advected_variable(icell,m) = cellhvars[proc][i++];
						}
#endif /* ADVECT_SPECIES */
	
						current_level_count++;
						proc_cur_cells[proc]++;
					}
				}

				/* if necessary, start a new receive for that proc */
				if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id ) {
					proc_page_count[proc] = min( page_size, proc_num_cells[proc] - proc_cur_cells[proc] );
					cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
					MPI_Irecv( cellhvars[proc], num_hydro_vars*proc_page_count[proc], MPI_FLOAT, 
						proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
					num_requests++;
				} else {
					proc_page_count[proc] = 0;
					cart_free( cellhvars[proc] );
					requests[proc] = MPI_REQUEST_NULL;
				}
	
				flag = 0;
			}

			if ( local_proc_id == file_parent ) {
				if ( num_requests > 0 ) {
					/* see if sends have completed */
					MPI_Testall( num_send_requests, send_requests, &ready_to_read, MPI_STATUSES_IGNORE );
				} else {
					MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
					num_send_requests = 0;
					ready_to_read = 1;
				}
	
				if ( current_read_count >= total_cells[level] && 
						current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			} else {
				if ( current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			}
		}

		if ( local_proc_id == file_parent ) {
			if ( !ready_to_read ) {
				MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
			}
	
			cart_free( cellhvars_buffer );
			fread( &size, sizeof(int), 1, input );
		}
#endif /* HYDRO */

#ifdef GRAVITY
		if ( local_proc_id == file_parent ) {
			fread( &size, sizeof(int), 1, input );

			cellvars_buffer = cart_alloc( 
					2*min( total_cells[level], page_size ) * sizeof(float) );
		}
	
		num_requests = 0;

		for ( proc = 0; proc < num_procs; proc++ ) {
			proc_cur_cells[proc] = 0;

			if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
				/* set up receive */
				proc_page_count[proc] = min( page_size - first_page_count[proc] % page_size, proc_num_cells[proc] );
				cellvars[proc] = cart_alloc( 2*min( page_size, proc_num_cells[proc] )*sizeof(float) );
				MPI_Irecv( cellvars[proc], 2*proc_page_count[proc], MPI_FLOAT, proc, 
					proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
				num_requests++;
			} else {
				proc_page_count[proc] = 0;
				requests[proc] = MPI_REQUEST_NULL;
			}
		}
	
		current_level_count = 0;
		current_read_count = 0;
		continue_reading = 1;
		ready_to_read = 1;

		while ( continue_reading ) {
			flag = 0;

			if ( local_proc_id == file_parent && current_read_count < total_cells[level] && ready_to_read ) {
				/* read in a page */
				page_count = min( page_size, total_cells[level] - current_read_count );
				num_read = fread( cellvars_buffer, sizeof(float), 2*page_count, input );
	
				if ( num_read != 2*page_count ) {
					cart_error("I/O Error in read_grid_binary: num_read = %u, 2*page_count = %u",
						num_read, 2*page_count );
				}
	
				if ( endian ) {
					for ( i = 0; i < 2*page_count; i++ ) {
						reorder( (char *)&cellvars_buffer[i], sizeof(float) );
					}
				}
	
				/* send info to other processors */
				start = 0;
				num_send_requests = 0;
				for ( proc = 0; proc < num_procs; proc++ ) {
					if ( start < page_count && current_read_count + start >= proc_first_index[proc] &&
							current_read_count + start < proc_first_index[proc+1] ) {
				
						count = min( proc_first_index[proc+1], current_read_count+page_count ) -
							max( proc_first_index[proc], current_read_count + start );
						cart_assert( count > 0 && count <= page_count );

						if ( proc == local_proc_id ) {
							/* copy into local buffer */
							proc_page_count[local_proc_id] = count;

							cellvars[local_proc_id] = 
								cart_alloc( 2*proc_page_count[local_proc_id]*sizeof(float));
	
							for ( i = 0; i < 2*proc_page_count[local_proc_id]; i++ ) {
								cart_assert( 2*start + i < 2*page_count );
								cellvars[local_proc_id][i] = cellvars_buffer[2*start+i];
							}
						} else {
							MPI_Isend( &cellvars_buffer[2*start], 2*count, MPI_FLOAT, proc,
								count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
						}
	
						start += count;
					}
				}
	
				cart_assert( start == page_count );
				current_read_count += page_count;
				ready_to_read = 0;
	
				/* unpack local cells first */
				if ( proc_page_count[local_proc_id] > 0 ) {
					flag = 1;
					proc = local_proc_id;
				} else if ( num_requests > 0 ) {
					/* see if we've received anything */
					MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
				}
			} else if ( num_requests > 0 ) {
				/* wait for a receive to complete */
				MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
				num_requests--;
				flag = 1;
			}
	
			if ( flag == 1 && proc != MPI_UNDEFINED ) {
				cart_assert( proc >= 0 && proc < num_procs );

				/* unpack received page */
				i = 0;
	
				while ( i < 2*proc_page_count[proc] ) {
					ioct = current_level[ proc_cell_index[proc] + proc_cur_cells[proc]/num_children];
					cart_assert( ioct >= 0 && ioct < num_octs );

					for ( j = 0; j < num_children; j++ ) {	
						icell = oct_child( ioct, j );
						cart_assert( icell >= 0 && icell < num_cells );
						cart_assert( cell_level(icell) == level );

						cell_potential(icell) = cellvars[proc][i++];
						cell_potential_hydro(icell) = cellvars[proc][i++];
	
						current_level_count++;
						proc_cur_cells[proc]++;
					}
				}

				/* if necessary, start a new receive for that proc */
				if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id ) {
					proc_page_count[proc] = min( page_size, proc_num_cells[proc] - proc_cur_cells[proc] );
					cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
					MPI_Irecv( cellvars[proc], 2*proc_page_count[proc], MPI_FLOAT, 
						proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
					num_requests++;
				} else {
					proc_page_count[proc] = 0;
					cart_free( cellvars[proc] );
					requests[proc] = MPI_REQUEST_NULL;
				}
	
				flag = 0;
			}

			if ( local_proc_id == file_parent ) {
				if ( num_requests > 0 ) {
					/* see if sends have completed */
					MPI_Testall( num_send_requests, send_requests, &ready_to_read, MPI_STATUSES_IGNORE );
				} else {
					MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
					num_send_requests = 0;
					ready_to_read = 1;
				}
	
				if ( current_read_count >= total_cells[level] && 
						current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			} else {
				if ( current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			}
		}

		if ( local_proc_id == file_parent ) {
			if ( !ready_to_read ) {
				MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
			}
	
			cart_free( cellvars_buffer );
			fread( &size, sizeof(int), 1, input );
		}
#endif /* GRAVITY */

		cart_free( current_level );		
	}

	if ( local_proc_id == file_parent ) {
		fclose( input );
	}

	cart_free( order );
}

void read_indexed_grid( char *filename, int num_sfcs, int *sfc_list, int max_level_to_read ) {
	int i, j, k;
	int size;
	int num_read, flag;
	FILE *input;
	char job[256];
	int minlevel, maxlevel;
	double t, dt;
	float adum, ainit;
        float boxh, Om0, Oml0, Omb0, h;
	int nextras;
	float extra[10];
	char lextra[10][256];
	int *order;
	int *current_level;
	int endian;
	int proc, sfc;
	int ret;
	int level;
	int icell;
	int last_sfc;
	int ncell0, sfc_order;
	long current_level_count, current_read_count;
	long next_level_count;
	int local_file_root_cells;
	long *root_cell_index;
	int *cell_refined;
	float *vars;
	int *cells_per_level;
	int num_cells_unpacked, num_refined_unpacked;
	char varname[32];
	int *varindex;
	long celloffset;
	double total_work;
	int *cell_count;
	float *cell_work;
	long total_cells, total_refined;
	int cells_to_read, cell_refined_to_read;
	int num_file_vars;
	int neighbors[num_dependent_neighbors];
	int *buffer_indices;
	skiplist *buffer_list;
	int index, hash;

	max_level_to_read = min( max_level, max_level_to_read );

	/* open file handle if parent of parallel file */
	if ( local_proc_id == MASTER_NODE ) {
		input = fopen(filename,"r");
		if ( input == NULL ) {
			cart_error( "Unable to open file %s for reading!", filename );
		}

                fread(&size, sizeof(int), 1, input );
		endian = 0;
		if ( size != 256 ) {
			reorder( (char *)&size, sizeof(int) );
			if ( size != 256 ) {
				cart_error("Error: file %s is corrupted", filename );
			} else {
				endian = 1;
				cart_debug("Reordering bytes (file endianness is opposite program)");
			}
		}

                fread(&job, sizeof(char), 256, input );
                fread(&size, sizeof(int), 1, input );

		/* istep, t, dt, adum, ainit */
		fread( &size, sizeof(int), 1, input );
		fread( &step, sizeof(int), 1, input );
		fread( &t, sizeof(double), 1, input );
		fread( &dt, sizeof(double), 1, input );
		fread( &adum, sizeof(float), 1, input );
		fread( &ainit, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&step, sizeof(int) );
			reorder( (char *)&ainit, sizeof(float) );
		}

		a_init = ainit;

                /* boxh, Om0, Oml0, Omb0, hubble */
		fread( &size, sizeof(int), 1, input );
		fread( &boxh, sizeof(float), 1, input );
		fread( &Om0, sizeof(float), 1, input );
		fread( &Oml0, sizeof(float), 1, input );
		fread( &Omb0, sizeof(float), 1, input );
		fread( &h, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&boxh, sizeof(float) );
			reorder( (char *)&Om0, sizeof(float) );
			reorder( (char *)&Oml0, sizeof(float) );
			reorder( (char *)&Omb0, sizeof(float) );
			reorder( (char *)&h, sizeof(float) );
		}

		Lbox = boxh;
		Omega0 = Om0;
		OmegaL0 = Oml0;
		Omegab0 = Omb0;
		hubble = h;

		/* nextra (no evidence extras are used...) extra lextra */
		fread( &size, sizeof(int), 1, input );
		fread( &nextras, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&nextras, sizeof(int) );
		}

		/* extra */
		fread( &size, sizeof(int), 1, input );
		fread( extra, sizeof(float), nextras, input );
		fread( &size, sizeof(int), 1, input );

		/* lextra */
		fread( &size, sizeof(int), 1, input );
		fread( lextra, 256*sizeof(char), nextras, input );
		fread( &size, sizeof(int), 1, input );

		/* Minlevel, MaxLevelNow */
		fread( &size, sizeof(int), 1, input );
		fread( &minlevel, sizeof(int), 1, input );
		fread( &maxlevel, sizeof(int), 1, input);

		if ( endian ) {
			reorder( (char *)&minlevel, sizeof(int) );
			reorder( (char *)&maxlevel, sizeof(int) );
		}

		if ( maxlevel > max_level ) {
			cart_error("File %s has more levels than compiled program (%u)", filename, maxlevel );
		}

		cart_assert( minlevel == min_level );

		fread( &size, sizeof(int), 1, input );

		/* tl */
		fread( &size, sizeof(int), 1, input );
		fread( &tl, sizeof(double), maxlevel-minlevel+1, input );
		fread( &size, sizeof(int), 1, input);

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&tl[i], sizeof(double) );
			}
		}

		/* dtl */
		fread( &size, sizeof(int), 1, input );
		fread( &dtl, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&dtl[i], sizeof(double) );
			}
		}

		/* tl_old */
		fread( &size, sizeof(int), 1, input );
		fread( &tl_old, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&tl_old[i], sizeof(double) );
			}
		}

		/* dtl_old */
		fread( &size, sizeof(int), 1, input );
		fread( &dtl_old, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );	

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&dtl_old[i], sizeof(double) );
			}
		}

		/* iSO */
		fread( &size, sizeof(int), 1, input );
		fread( &level_sweep_dir, sizeof(int), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&level_sweep_dir[i], sizeof(int) );
			}
		}

		/* sfc ordering used */
		fread( &size, sizeof(int), 1, input );
		fread( &sfc_order, sizeof(int), 1, input);

		if ( endian ) {
			reorder( (char *)&sfc_order, sizeof(int) );
		}

		if ( sfc_order != SFC ) {
			cart_error("File has different sfc indexing than program");
		}
		fread( &size, sizeof(int), 1, input );

		/* refinement volume */
		fread( &size, sizeof(int), 1, input );
		fread( refinement_volume_min, sizeof(float), nDim, input );
		fread( refinement_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&refinement_volume_min[i], sizeof(float) );
				reorder( (char *)&refinement_volume_max[i], sizeof(float) );
			}
		}

#ifdef STARFORM
		/* star formation volume */
		fread( &size, sizeof(int), 1, input );
		fread( star_formation_volume_min, sizeof(float), nDim, input );
		fread( star_formation_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&star_formation_volume_min[i], sizeof(float) );
				reorder( (char *)&star_formation_volume_max[i], sizeof(float) );
			}
		}
#endif /* STARFORM */

		/* read in variable strings */
		size = sizeof(int);
		fread( &size, sizeof(int), 1, input );
		fread( &num_file_vars, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&num_file_vars, sizeof(int) );
		}

		cart_debug("num_file_vars = %d", num_file_vars );

		varindex = cart_alloc( num_file_vars*sizeof(int) );

		size = 32*sizeof(char);
		for ( i = 0; i < num_file_vars; i++ ) {
			fread( &size, sizeof(int), 1, input );
			fread( varname, sizeof(char), 32, input );
			fread( &size, sizeof(int), 1, input );

			cart_debug("variable[%u] = %s", i, varname );

			if ( strncmp( varname, "hydro_gas_density",32 ) == 0 ) {
				varindex[i] = HVAR_GAS_DENSITY;
			} else if ( strncmp( varname, "hydro_gas_energy", 32 ) == 0 ) {
				varindex[i] = HVAR_GAS_ENERGY;
			} else if ( strncmp( varname, "hydro_gas_pressure", 32 ) == 0 ) {
				varindex[i] = HVAR_PRESSURE;
			} else if ( strncmp( varname, "hydro_gas_gamma", 32 ) == 0 ) {
				varindex[i] = HVAR_GAMMA;
			} else if ( strncmp( varname, "hydro_gas_internal_energy", 32 ) == 0 ) {
				varindex[i] = HVAR_INTERNAL_ENERGY;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			} else if ( strncmp( varname, "hydro_electron_internal_energy", 32 ) == 0 ) {
				varindex[i] = HVAR_ELECTRON_INTERNAL_ENERGY;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
			} else if ( strncmp( varname, "hydro_momentum_x", 32 ) == 0 ) {
				varindex[i] = HVAR_MOMENTUM;
			} else if ( strncmp( varname, "hydro_momentum_y", 32 ) == 0 ) {
				varindex[i] = HVAR_MOMENTUM+1;
			} else if ( strncmp( varname, "hydro_momentum_z", 32 ) == 0 ) {
				varindex[i] = HVAR_MOMENTUM+2;
#ifdef ENRICH
			} else if ( strncmp( varname, "hydro_metallicity_II", 32 ) == 0 ) {
				varindex[i] = HVAR_METALLICITY_II;
			} else if ( strncmp( varname, "hydro_metallicity_Ia", 32 ) == 0 ) {
				varindex[i] = HVAR_METALLICITY_Ia;
#endif /* ENRICH */
			} else {
				varindex[i] = -1;
			}
		}

		/* ncell0 */
		fread( &size, sizeof(int), 1, input );
		fread( &ncell0, sizeof(int), 1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&ncell0, sizeof(int) );
		}

		if ( ncell0 != (num_grid*num_grid*num_grid) ) {
			cart_error("File has different num_grid than compiled program");
		}

		cart_debug("ncell0 = %d", ncell0 );
	}

	/* send header information to all other processors */
	MPI_Bcast( &endian, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &minlevel, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &maxlevel, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &step, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &a_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Lbox, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Omega0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &OmegaL0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Omegab0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &hubble, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( tl, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( dtl, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( dtl_old, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( level_sweep_dir, max_level-min_level+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	MPI_Bcast( refinement_volume_min, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( refinement_volume_max, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );

#ifdef STARFORM
	MPI_Bcast( star_formation_volume_min, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( star_formation_volume_max, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
#endif /* STARFORM */

	MPI_Bcast( &num_file_vars, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	
	if ( local_proc_id != MASTER_NODE ) {
		varindex = cart_alloc( num_file_vars*sizeof(int) );
	}

	MPI_Bcast( varindex, num_file_vars, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	cells_per_level = cart_alloc( (maxlevel-minlevel+1)*sizeof(int) );

	if ( local_proc_id == MASTER_NODE ) {
		/* pick out indices we need */
		fread( &size, sizeof(int), 1, input );
		celloffset = ftell(input);

		root_cell_index = cart_alloc( num_sfcs * sizeof(long) );
		
		last_sfc = -1;
		for ( i = 0; i < num_sfcs; i++ ) {
			if ( sfc_list[i] != (last_sfc+1) ) {
				fseek( input, celloffset+(long)sfc_list[i]*sizeof(long), SEEK_SET );
			}

			fread( &root_cell_index[i], sizeof(long), 1, input );
			last_sfc = sfc_list[i];
		}
	
		if ( endian ) {
			for ( i = 0; i < num_sfcs; i++ ) {
				reorder( (char *)&root_cell_index[i], sizeof(long) );
			}
		}

	}

	cart_debug("now reading in selected sfc trees...");

	last_sfc = -2;
	for ( i = 0; i < num_sfcs; i++ ) {
		if ( local_proc_id == MASTER_NODE ) {
			if ( sfc_list[i] != last_sfc+1 ) {
				fseek( input, root_cell_index[i], SEEK_SET );
			}

			/* read in root tree */
			fread( &size, sizeof(int), 1, input );
			cart_assert( size == (maxlevel-minlevel+1)*sizeof(int));
			fread( cells_per_level, sizeof(int), maxlevel-minlevel+1, input );
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				for ( level = minlevel; level <= maxlevel; level++ ) {
					reorder( (char *)&cells_per_level, sizeof(int) );
				}
			}

			cells_to_read = 0;
			for ( level = 0; level <= min( maxlevel, max_level_to_read ); level++ ) {
				cells_to_read += cells_per_level[level];
			}

			total_cells = 0;
			for ( level = 0; level <= maxlevel; level++ ) {
				total_cells += cells_per_level[level];
			}

			total_refined = total_cells - cells_per_level[maxlevel];
			cell_refined_to_read = cells_to_read-cells_per_level[min( maxlevel, max_level_to_read)];

			cell_refined = cart_alloc( cell_refined_to_read*sizeof(int) );
			vars = cart_alloc( cells_to_read*num_file_vars*sizeof(float) );

			fread( &size, sizeof(int), 1, input );
			fread( cell_refined, sizeof(int), cell_refined_to_read, input );

			if ( total_cells > cells_to_read ) {
				fseek( input, (long)(total_refined-cell_refined_to_read)*sizeof(int), SEEK_CUR );
			}
			fread( &size, sizeof(int), 1, input );

			fread( &size, sizeof(int), 1, input );
			fread( vars, sizeof(float), num_file_vars*cells_to_read, input );
			if ( total_cells > cells_to_read ) {
				fseek( input, (long)(total_cells-cells_to_read)*
						(long)num_file_vars*sizeof(float), SEEK_CUR );
			}
			fread( &size, sizeof(int), 1, input );
			
			if ( num_procs > 1 && !root_cell_is_local( sfc_list[i] ) ) {
				/* send to owner */
				proc = processor_owner( sfc_list[i] );
				MPI_Send( &cell_refined_to_read, 1, MPI_INT, proc, sfc_list[i], MPI_COMM_WORLD );
				MPI_Send( &cells_to_read, 1, MPI_INT, proc, sfc_list[i], MPI_COMM_WORLD );
				MPI_Send( cell_refined, cell_refined_to_read, MPI_INT, proc, sfc_list[i], MPI_COMM_WORLD );
				MPI_Send( vars, cells_to_read*num_file_vars, MPI_FLOAT, proc, sfc_list[i], MPI_COMM_WORLD );

				cart_free( cell_refined );
				cart_free( vars );
			}

			last_sfc = sfc_list[i];
		} else {
			if ( root_cell_is_local( sfc_list[i] ) ) {
				/* receive from root */
				MPI_Recv( &cell_refined_to_read, 1, MPI_INT, MASTER_NODE, sfc_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				MPI_Recv( &cells_to_read, 1, MPI_INT, MASTER_NODE, sfc_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE );

				cell_refined = cart_alloc( cell_refined_to_read*sizeof(int) );
				vars = cart_alloc( cells_to_read*num_file_vars*sizeof(float) );

				MPI_Recv( cell_refined, cell_refined_to_read, MPI_INT, MASTER_NODE, sfc_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				MPI_Recv( vars, cells_to_read*num_file_vars, MPI_FLOAT, MASTER_NODE, sfc_list[i], MPI_COMM_WORLD, MPI_STATUS_IGNORE );
			}
		}

		if ( num_procs == 1 || root_cell_is_local( sfc_list[i] ) ) {
			/* unpack cells */
			level = minlevel;
			current_level_count = 1;
			num_refined_unpacked = 0;
			num_cells_unpacked = 0;

			order = cart_alloc( cells_to_read * sizeof(int) );
			order[0] = root_cell_location( sfc_list[i] );
			cart_assert( order[0] >= 0 && order[0] < num_cells );

			while ( current_level_count > 0 ) {
				next_level_count = 0;

				if ( level < min( maxlevel, max_level_to_read ) ) {
					/* refine this level */
					for ( j = 0; j < current_level_count; j++ ) {
						cart_assert( num_refined_unpacked < cell_refined_to_read );
						if ( cell_refined[num_refined_unpacked++] ) {
							ret = split_cell( order[num_cells_unpacked+j] );

							if ( ret ) {
								cart_error("Unable to split cell!");
							}

							for ( k = 0; k < num_children; k++ ) {
								order[num_cells_unpacked+current_level_count+next_level_count++] =
									cell_child( order[num_cells_unpacked+j], k );
							}
						}
					}
				}
				
				/* unpack cell variables */
				for ( k = 0; k < num_file_vars; k++ ) {
					for ( j = 0; j < current_level_count; j++ ) {
						cell_var( order[num_cells_unpacked+j], varindex[k] ) = 
							vars[k*cells_to_read+(num_cells_unpacked+j)];
					}
				}

				num_cells_unpacked += current_level_count;
				current_level_count = next_level_count;
				level++;
			}

			cart_assert( num_cells_unpacked == cells_to_read );
			
			cart_free( order );
			cart_free( cell_refined );
			cart_free( vars );
		}
	}

	cart_free( cells_per_level );
	cart_free( varindex );

	if ( local_proc_id == MASTER_NODE ) {
		cart_free( root_cell_index );
		fclose(input);
	}

	buffer_enabled = 1;
}

void write_hart_gas_binary( char *filename ) {
	int i, j, k, m;
	int size;
	FILE *output;
	float adum, ainit;
	float boxh, Om0, Oml0, Omb0, h;
	int minlevel, maxlevel;
	int nextras;
	int *cellrefined;
	float *cellhvars;
	float *cellvars;
	int level;
	int coords[nDim];
	int icell, ioct;
	int ncell0;
	int current_level_count;
	int page_size;
	int iOctFree, nOct;
	int iNOLL, iHOLL;
	int pos[nDim];
	int inext, iprev;
	int icellnum;
	int neighbors[num_neighbors];
	int parent_cell;
	int iOctCh;

	/* this is a single-processor only function */
	cart_assert( num_procs == 1 );

	/* should replace this with a memory parameters and use it in trade_particle_lists as well */	
	page_size = num_grid*num_grid;

	/* allocate pages for writing */
        cellrefined = cart_alloc( page_size * sizeof(int) );
        cellhvars = cart_alloc( num_hydro_vars * page_size * sizeof(float) );
        cellvars = cart_alloc( 2 * page_size * sizeof(float) );

	minlevel = min_level;
	maxlevel = max_level_now();

	/* open file and write header */
	output = fopen(filename,"w");
	if ( output == NULL ) {
		cart_error( "Unable to open file %s for writing!", filename );
	}

	size = 256*sizeof(char);
	fwrite(&size, sizeof(int), 1, output );
	fwrite(&jobname, sizeof(char), 256, output );
	fwrite(&size, sizeof(int), 1, output );

	/* istep, t, dt, adum, ainit */
	adum = aexp[min_level];
	ainit = a_init;
	size = sizeof(int) + 2*sizeof(double) + 2*sizeof(float);

	fwrite(&size, sizeof(int), 1, output );
	fwrite( &step, sizeof(int), 1, output );
	fwrite( &tl[min_level], sizeof(double), 1, output );
	fwrite( &dtl[min_level], sizeof(double), 1, output );
	fwrite( &adum, sizeof(float), 1, output );
	fwrite( &ainit, sizeof(float), 1, output );
	fwrite(&size, sizeof(int), 1, output );

	/* boxh, Om0, Oml0, Omb0, hubble */
	boxh = Lbox;
	Om0 = Omega0;
	Oml0 = OmegaL0;
	Omb0 = Omegab0;
	h = hubble;
	size = 5*sizeof(float);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &boxh, sizeof(float), 1, output );
	fwrite( &Om0, sizeof(float), 1, output );
	fwrite( &Oml0, sizeof(float), 1, output );
	fwrite( &Omb0, sizeof(float), 1, output );
	fwrite( &h, sizeof(float), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	/* nextra (no evidence extras are used...) extra lextra */
	size = sizeof(int);
	nextras = 0;

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &nextras, sizeof(int), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	/* extra */
	size = nextras * sizeof(float);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	/* lextra */
	size = nextras * 256 * sizeof(char);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &size, sizeof(int), 1, output );
	
	/* Minlevel, MaxLevelNow */
	size = 2 * sizeof(int);
	fwrite(&size, sizeof(int), 1, output );
	fwrite(&minlevel, sizeof(int), 1, output );
	fwrite(&maxlevel, sizeof(int), 1, output );
	fwrite(&size, sizeof(int), 1, output );

	/* tl */
	size = (maxlevel-minlevel+1) * sizeof(double);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &tl, sizeof(double), maxlevel-minlevel+1, output);
	fwrite( &size, sizeof(int), 1, output );

	/* dtl */
	size = (maxlevel-minlevel+1) * sizeof(double);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &dtl, sizeof(double), maxlevel-minlevel+1, output);
	fwrite( &size, sizeof(int), 1, output );

        /* tl_old */
        size = (maxlevel-minlevel+1) * sizeof(double);

        fwrite( &size, sizeof(int), 1, output );
        fwrite( &tl_old, sizeof(double), maxlevel-minlevel+1, output);
        fwrite( &size, sizeof(int), 1, output );

	/* dtl_old */
	size = (maxlevel-minlevel+1) * sizeof(double);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &dtl_old, sizeof(double), maxlevel-minlevel+1, output);
	fwrite( &size, sizeof(int), 1, output );

	/* iSO */
	size = (maxlevel-minlevel+1) * sizeof(int);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &level_sweep_dir, sizeof(int), maxlevel-minlevel+1, output);
	fwrite( &size, sizeof(int), 1, output );

	/* ncell0 */
	ncell0 = num_grid*num_grid*num_grid;
	size = sizeof(int);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &ncell0, sizeof(int), 1, output);
	fwrite( &size, sizeof(int), 1, output );

	/* now we start writing pages of root level cell children */
	size = ncell0 * sizeof(int);
	fwrite( &size, sizeof(int), 1, output );

	/* holds list of next level octs to write */
	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
		i = 0;
		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				icell = sfc_index( coords );
				cellrefined[i++] = cell_child_oct[icell] + 1;
			}
		}

		fwrite( cellrefined, sizeof(int), page_size, output );
	}

	fwrite( &size, sizeof(int), 1, output );

	/* now write pages of root level hydro variables */
	size = num_hydro_vars * ncell0 * sizeof(float);
	fwrite( &size, sizeof(int), 1, output );

	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
		i = 0;
		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				icell = sfc_index( coords );

				cellhvars[i++] = cell_gas_density(icell);
				cellhvars[i++] = cell_gas_energy(icell);
				cellhvars[i++] = cell_momentum(icell,0);
				cellhvars[i++] = cell_momentum(icell,1);
				cellhvars[i++] = cell_momentum(icell,2);
				cellhvars[i++] = cell_gas_pressure(icell);
				cellhvars[i++] = cell_gas_gamma(icell);
				cellhvars[i++] = cell_gas_internal_energy(icell);

#ifdef ADVECT_SPECIES
				for ( j = 0; j < num_chem_species; j++ ) {
					cellhvars[i++] = cell_advected_variable(icell,j);
				}
#endif /* ADVECT_SPECIES */

			}
		}

		fwrite( cellhvars, sizeof(float), num_hydro_vars*page_size, output );
	}

	fwrite( &size, sizeof(int), 1, output );

#ifdef GRAVITY
	/* finally write potential variables */
	size = 2 * ncell0 * sizeof(float);
	fwrite( &size, sizeof(int), 1, output );

	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
		i = 0;
		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				icell = sfc_index( coords );

				cellvars[i++] = cell_potential(icell);
				cellvars[i++] = cell_potential_hydro(icell);
			}
		}

		fwrite( cellvars, sizeof(float), 2*page_size, output );
	}

	fwrite ( &size, sizeof(int), 1, output );
#endif /* GRAVITY */

	size = 2*sizeof(int);

	fwrite( &size, sizeof(int), 1, output );

	nOct = 0;
	iOctFree = 0;
	for ( ioct = 0; ioct < num_octs; ioct++ ) {
		if ( oct_level[ioct] != FREE_OCT_LEVEL ) {
			nOct++;
			iOctFree = ioct+2;
		}
	}
	fwrite( &iOctFree, sizeof(int), 1, output );
	fwrite( &nOct, sizeof(int), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	/* then write each level's cells in turn */
	for ( level = min_level+1; level <= maxlevel; level++ ) {
		/* write size */
		size = 3*sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &level, sizeof(int), 1, output );

		iHOLL = local_oct_list[level] + 1;
		iNOLL = num_cells_per_level[level] / num_children;

		fwrite( &iNOLL, sizeof(int), 1, output );
		fwrite( &iHOLL, sizeof(int), 1, output );

		fwrite( &size, sizeof(int), 1, output );

		ioct = local_oct_list[level];
		for ( i = 0; i < iNOLL; i++ ) {
			size = 	3*sizeof(int) + 		/* pos */
				num_neighbors*sizeof(int) + 	/* neighbors */
				sizeof(int) +			/* parent */
				sizeof(int) +			/* level */
				sizeof(int) +			/* next */
				sizeof(int);			/* prev */

			for ( j = 0; j < nDim; j++ ) {
				pos[j] = (int)((oct_pos[ioct][j]+1.0)*(1<<(max_level+1)));
			}

			for ( j = 0; j < num_neighbors; j++ ) {
				/* if neighbor is root cell need to convert to index */
				if ( cell_level(oct_neighbors[ioct][j]) == min_level ) {
					sfc_coords( cell_parent_root_sfc( oct_neighbors[ioct][j] ), coords );
					neighbors[j] = coords[2]+num_grid*(coords[1]+num_grid*coords[0])+1;
				} else {
					neighbors[j] = oct_neighbors[ioct][j] + 1;
				}
			}

			if ( oct_level[ioct] == min_level+1 ) {
				sfc_coords( oct_parent_root_sfc[ioct], coords );
				parent_cell = coords[2]+num_grid*(coords[1]+num_grid*coords[0])+1;
			} else {
				parent_cell = oct_parent_cell[ioct] + 1;
			}

			fwrite( &size, sizeof(int), 1, output );
			fwrite( pos, sizeof(int), 3, output );
			fwrite( neighbors, sizeof(int), num_neighbors, output );
			fwrite( &parent_cell, sizeof(int), 1, output );
			fwrite( &oct_level[ioct], sizeof(int), 1, output );

			inext = oct_next[ioct] + 1;
			iprev = oct_prev[ioct] + 1;

			fwrite( &inext, sizeof(int), 1, output );
			fwrite( &iprev, sizeof(int), 1, output );
			fwrite( &size, sizeof(int), 1, output );

			ioct = oct_next[ioct];
		}

		ioct = local_oct_list[level];
		for ( i = 0; i < iNOLL; i++ ) {
			for ( j = 0; j < num_children; j++ ) {
				icell = oct_child( ioct, j );

				size = 	sizeof(int) +			/* idc */
					sizeof(int) + 			/* iOctCh */
					num_hydro_vars*sizeof(float) +	/* hvar */
					2*sizeof(float);		/* var */

				icellnum = icell + 1;
				iOctCh = cell_child_oct[icell] + 1;

				cellhvars[0] = cell_gas_density(icell);
				cellhvars[1] = cell_gas_energy(icell);
				cellhvars[2] = cell_momentum(icell,0);
				cellhvars[3] = cell_momentum(icell,1);
				cellhvars[4] = cell_momentum(icell,2);
				cellhvars[5] = cell_gas_pressure(icell);
				cellhvars[6] = cell_gas_gamma(icell);
				cellhvars[7] = cell_gas_internal_energy(icell);

#ifdef ADVECT_SPECIES
				for ( m = 0; m < num_chem_species; m++ ) {
					cellhvars[8+m] = cell_advected_variable(icell,m);
				}
#endif /* ADVECT_SPECIES */

#ifdef GRAVITY
				cellvars[0] = cell_potential(icell);
				cellvars[1] = cell_potential_hydro(icell);
#endif /* GRAVITY */

				fwrite( &size, sizeof(int), 1, output );
				fwrite( &icellnum, sizeof(int), 1, output );
				fwrite( &iOctCh, sizeof(int), 1, output );
				fwrite( cellhvars, sizeof(float), num_hydro_vars, output );

#ifdef GRAVITY
				fwrite( cellvars, sizeof(float), 2, output );
#endif /* GRAVITY */

				fwrite( &size, sizeof(int), 1, output );
			}

			ioct = oct_next[ioct];
		}
	}

	fclose( output );

	cart_free( cellrefined );
	cart_free( cellhvars );
	cart_free( cellvars );
}



/*
// NG: Slightly re-written form of grid binary I/O that (a) eliminates 
// code duplication and (b) allows complete flexibility in what is written
// to disk (the latter is needed for adding RT block I/O).
*/


/* two helpers */
void write_grid_binary_top_level_vars(int num_out_vars, int *out_var, FILE *output, int file_parent, int file_num_procs, long *total_cells, int page_size, int *proc_num_cells);
void write_grid_binary_lower_level_vars(int num_out_vars, int *out_var, FILE *output, int file_parent, int file_num_procs, long *total_cells, int page_size, int *proc_num_cells, int level, int current_level_count, int *current_level);


void write_grid_binary2( char *filename ) {
        int i, j, k;
        int size;
	FILE *output;
	float adum, ainit;
	float boxh, Om0, Oml0, Omb0, h;
	int minlevel, maxlevel;
	int nextras;
	int *cellrefined;
	int *order;
	int *current_level;
	int proc;
	int level;
	int icell, ioct;
	int page_size;
	int ncell0, sfc_order;
	int current_level_count;
	int next_level_count;
	int page;
	int page_count;
	long total_cells[max_level-min_level+1];
	MPI_Status status;
	char parallel_filename[256];
	int proc_num_cells[MAX_PROCS*(max_level-min_level+1)];
	int file_index, file_parent, file_num_procs;

	int hydro_vars[num_hydro_vars];
	int num_other_vars = 0;
	int *other_vars = 0;

	/*
	// Maintain the same order as in a previous version
	*/
	hydro_vars[0] = HVAR_GAS_DENSITY;
	hydro_vars[1] = HVAR_GAS_ENERGY;
	hydro_vars[2] = HVAR_MOMENTUM + 0;
	hydro_vars[3] = HVAR_MOMENTUM + 1;
	hydro_vars[4] = HVAR_MOMENTUM + 2;
	hydro_vars[5] = HVAR_PRESSURE;
	hydro_vars[6] = HVAR_GAMMA;
	hydro_vars[7] = HVAR_INTERNAL_ENERGY;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	hydro_vars[8] = HVAR_ELECTRON_INTERNAL_ENERGY;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

#ifdef ADVECT_SPECIES
	for(j=0; j<num_chem_species; j++)
	  {
	    hydro_vars[num_hydro_vars-num_chem_species+j] = HVAR_ADVECTED_VARIABLES+j;
	  }
#endif /* ADVECT_SPECIES */

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
#ifdef GRAVITY
	num_other_vars += 2;
#endif
#ifdef RADIATIVE_TRANSFER
	num_other_vars += rt_num_disk_vars;
#endif

	other_vars = cart_alloc(num_other_vars*sizeof(int));

#ifdef GRAVITY
	other_vars[0] = VAR_POTENTIAL;
	other_vars[1] = VAR_POTENTIAL_HYDRO;
	k = 2;
#else
	k = 0;
#endif
#ifdef RADIATIVE_TRANSFER
	for(j=0; j<rt_num_disk_vars; j++) other_vars[k+j] = rt_disk_offset + j; 
#endif
#endif /* defined(GRAVITY) || defined(RADIATIVE_TRANSFER) */


	/* ensure consistency of num_output_files */
	num_output_files = min( num_output_files, num_procs );
	num_output_files = max( num_output_files, 1 );

	page_size = num_grid*num_grid;

	/* determine parallel output options */
	file_index = local_proc_id * num_output_files / num_procs;

	file_parent = local_proc_id;
	while ( file_parent > 0 && 
			(file_parent-1)*num_output_files / num_procs == file_index ) {
		file_parent--;
	}	

	file_num_procs = 1;
	while ( file_parent + file_num_procs < num_procs &&
			(file_parent+file_num_procs)*num_output_files / num_procs == file_index ) {
		file_num_procs++;
	}	

	cellrefined = cart_alloc( page_size * sizeof(int) );

	/* get number of cells on each level */
	MPI_Allgather( num_cells_per_level, max_level-min_level+1, MPI_INT,
		proc_num_cells, max_level-min_level+1, MPI_INT, MPI_COMM_WORLD );

	minlevel = min_level;
	maxlevel = max_level_now_global();

	/* open file handle if parent of parallel file */
	if ( local_proc_id == file_parent ) {
		if ( num_output_files == 1 ) {
			output = fopen(filename,"w");
		} else {
			sprintf( parallel_filename, "%s.%03u", filename, file_index );
			output = fopen(parallel_filename, "w");
		}

		if ( output == NULL ) {
			cart_error( "Unable to open file %s for writing!", filename );
		}

		for ( level = min_level; level <= max_level; level++ ) {
			total_cells[level] = 0;
			for ( proc = local_proc_id; proc < local_proc_id+file_num_procs; proc++ ) {
				total_cells[level] += proc_num_cells[(max_level-min_level+1)*proc+level];
			}
		}
	}

	/* only write one copy of the header information */
	if ( local_proc_id == MASTER_NODE ) {
		size = 256*sizeof(char);
		fwrite(&size, sizeof(int), 1, output );
		fwrite(&jobname, sizeof(char), 256, output );
		fwrite(&size, sizeof(int), 1, output );

		/* istep, t, dt, adum, ainit */
		adum = aexp[min_level];
		ainit = a_init;
		size = sizeof(int) + 2*sizeof(double) + 2*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &step, sizeof(int), 1, output );
		fwrite( &tl[min_level], sizeof(double), 1, output );
		fwrite( &dtl[min_level], sizeof(double), 1, output );
		fwrite( &adum, sizeof(float), 1, output );
		fwrite( &ainit, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* boxh, Om0, Oml0, Omb0, hubble */
		boxh = Lbox;
		Om0 = Omega0;
		Oml0 = OmegaL0;
		Omb0 = Omegab0;
		h = hubble;
		size = 5*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &boxh, sizeof(float), 1, output );
		fwrite( &Om0, sizeof(float), 1, output );
		fwrite( &Oml0, sizeof(float), 1, output );
		fwrite( &Omb0, sizeof(float), 1, output );
		fwrite( &h, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* nextra (no evidence extras are used...) extra lextra */
		size = sizeof(int);
		nextras = 0;

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &nextras, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* extra */
		size = nextras * sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* lextra */
		size = nextras * 256 * sizeof(char);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );
	
		/* Minlevel, MaxLevelNow */
		size = 2 * sizeof(int);
		fwrite(&size, sizeof(int), 1, output );
		fwrite(&minlevel, sizeof(int), 1, output );
		fwrite(&maxlevel, sizeof(int), 1, output );
		fwrite(&size, sizeof(int), 1, output );

		size = (maxlevel-minlevel+1) * sizeof(double);

		/* tl */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &tl, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* dtl */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &dtl, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* tl_old */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &tl_old, sizeof(double), maxlevel-minlevel+1, output );
		fwrite( &size, sizeof(int), 1, output );

		/* dtl_old */
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &dtl_old, sizeof(double), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* iSO */
		size = (maxlevel-minlevel+1) * sizeof(int);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &level_sweep_dir, sizeof(int), maxlevel-minlevel+1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* sfc ordering used */
		sfc_order = SFC;
		size = sizeof(int);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &sfc_order, sizeof(int), 1, output);
		fwrite( &size, sizeof(int), 1, output );

		/* refinement volume */
		size = 2*nDim*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( refinement_volume_min, sizeof(float), nDim, output );
		fwrite( refinement_volume_max, sizeof(float), nDim, output );
		fwrite( &size, sizeof(int), 1, output );

#ifdef STARFORM
		/* refinement volume */
		size = 2*nDim*sizeof(float);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( star_formation_volume_min, sizeof(float), nDim, output );
		fwrite( star_formation_volume_max, sizeof(float), nDim, output );
		fwrite( &size, sizeof(int), 1, output );
#endif /* STARFORM */

		/* ncell0 */
		ncell0 = num_grid*num_grid*num_grid;
		size = sizeof(int);

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &ncell0, sizeof(int), 1, output);
		fwrite( &size, sizeof(int), 1, output );
	}

	/* now start writing pages of root level cell children */
	if ( local_proc_id == file_parent ) {
		size = total_cells[min_level] * sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
	}

	/* holds list of next level octs to write */
	order = cart_alloc( (num_cells_per_level[min_level+1] / num_children) * sizeof(int) );
	next_level_count = 0;
	current_level_count = 0;
	page = 0;

	while ( current_level_count < num_cells_per_level[min_level] ) {
		page_count = min( page_size, num_cells_per_level[min_level] - current_level_count );

		for ( i = 0; i < page_count; i++ ) {
			if ( cell_is_refined(current_level_count) ) {
				cellrefined[i] = tree_cell_count(current_level_count);
				order[next_level_count++] = cell_child_oct[current_level_count];
			} else {
				cellrefined[i] = 1;
			}
			current_level_count++;
		}

		if ( local_proc_id == file_parent ) {
			fwrite( cellrefined, sizeof(int), page_count, output );
		} else {
			MPI_Send( cellrefined, page_count, MPI_INT, file_parent, page, MPI_COMM_WORLD );
			page++;
		}
	}

	if ( local_proc_id == file_parent ) {
		for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
			i = 0;
			page = 0;
			while ( i < proc_num_cells[(max_level-min_level+1)*proc+min_level] ) {
				page_count = min( page_size,  proc_num_cells[(max_level-min_level+1)*proc+min_level] - i );
				MPI_Recv( cellrefined, page_count, MPI_INT, proc, page, MPI_COMM_WORLD, &status );
				fwrite( cellrefined, sizeof(int), page_count, output );
				i += page_count;
				page++;
			}
		}

		fwrite( &size, sizeof(int), 1, output );
	}

	/* now write pages of root level hydro variables */
	write_grid_binary_top_level_vars(num_hydro_vars,hydro_vars,output,file_parent,file_num_procs,total_cells,page_size,proc_num_cells);

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
	write_grid_binary_top_level_vars(num_other_vars,other_vars,output,file_parent,file_num_procs,total_cells,page_size,proc_num_cells);
#endif /* defined(GRAVITY) || defined(RADIATIVE_TRANSFER) */

	/* then write each level's cells in turn */
	for ( level = min_level+1; level <= maxlevel; level++ ) {
		if ( local_proc_id == file_parent ) {
			/* write size */
			size = sizeof(long);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &total_cells[level], sizeof(long), 1, output );
			fwrite( &size, sizeof(int), 1, output );

			size = total_cells[level] * sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
		}

		current_level_count = next_level_count;
		cart_assert( num_cells_per_level[level] == current_level_count*num_children );
		next_level_count = 0;
		current_level = order;

		if ( level < maxlevel ) {
			order = cart_alloc( ( num_cells_per_level[level+1]/num_children) * sizeof(int) );
		} else {
			order = cart_alloc( 0 );
		}

		i = 0;
		page = 0;
		while ( i < current_level_count ) {
			page_count = min( page_size, num_cells_per_level[level] - i*num_children );

			j = 0;
			while ( j < page_count) {
				ioct = current_level[i];

				for ( k = 0; k < num_children; k++ ) {
					icell = oct_child( ioct, k );

					if ( cell_is_refined(icell) ) {
						cellrefined[j++] = 1;
						order[next_level_count++] = cell_child_oct[icell];
					} else {
						cellrefined[j++] = 0;
					}
				}

				i++;
			}

			if ( local_proc_id == file_parent ) {
				fwrite( cellrefined, sizeof(int), page_count, output );
			} else {
				MPI_Send( cellrefined, page_count, MPI_INT, file_parent, page, MPI_COMM_WORLD );
				page++;
			}
		}

		if ( local_proc_id == file_parent ) {
			for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ ) {
				i = 0;
				page = 0;
				while ( i < proc_num_cells[(max_level-min_level+1)*proc+level] ) {
					page_count = min( page_size, 
							proc_num_cells[(max_level-min_level+1)*proc+level] - i );
					MPI_Recv( cellrefined, page_count, MPI_INT, proc, page, MPI_COMM_WORLD, &status );
					fwrite( cellrefined, sizeof(int), page_count, output );
					i += page_count;
					page++;
				}
			}

			fwrite( &size, sizeof(int), 1, output );
		}

		/* now write pages of lower level hydro variables */
		write_grid_binary_lower_level_vars(num_hydro_vars,hydro_vars,output,file_parent,file_num_procs,total_cells,page_size,proc_num_cells,level,current_level_count,current_level);

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
		write_grid_binary_lower_level_vars(num_other_vars,other_vars,output,file_parent,file_num_procs,total_cells,page_size,proc_num_cells,level,current_level_count,current_level);
#endif /* defined(GRAVITY) || defined(RADIATIVE_TRANSFER) */

		cart_free( current_level );
	}

	if ( local_proc_id == file_parent ) {
		fclose( output );
	}

	if(other_vars != 0) cart_free(other_vars);

	cart_free( order );
	cart_free( cellrefined );

#ifdef RADIATIVE_TRANSFER
	/* Save RF data */
	rtWriteRadiationFieldData(filename,1);
#endif
}


void write_grid_binary_top_level_vars(int num_out_vars, int *out_var, FILE *output, int file_parent, int file_num_procs, long *total_cells, int page_size, int *proc_num_cells)
{
  int i, j;
  int size;
  float *cellvars;
  int current_level_count;
  int page, page_count;
  int icell;
  int proc;
  MPI_Status status;

  if(num_out_vars < 1) return;

  cellvars = cart_alloc( num_out_vars * page_size * sizeof(float) );

  /* now write pages of root level hydro variables */
  if ( local_proc_id == file_parent )
    {
      size = num_out_vars * total_cells[min_level] * sizeof(float);
      fwrite( &size, sizeof(int), 1, output );
    }

  current_level_count = 0;
  page = 0;
  while ( current_level_count < num_cells_per_level[min_level] )
    {
      page_count = min( page_size, num_cells_per_level[min_level] - current_level_count );
    
      i = 0;
      while ( i < num_out_vars*page_count )
	{
	  icell = current_level_count;

	  for(j=0; j<num_out_vars; j++)
	    {
	      cellvars[i++] = cell_var(icell,out_var[j]);
	    }

	  current_level_count++;
	}

      if ( local_proc_id == file_parent )
	{
	  fwrite( cellvars, sizeof(float), num_out_vars*page_count, output );
	} 
      else
	{
	  MPI_Send( cellvars, num_out_vars*page_count, MPI_FLOAT, file_parent, page, MPI_COMM_WORLD );
	  page++;
	}
    }

  if ( local_proc_id == file_parent )
    {
      for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ )
	{
	  i = 0;
	  page = 0;
	  while ( i <  proc_num_cells[(max_level-min_level+1)*proc+min_level] )
	    {
	      page_count = min( page_size,  proc_num_cells[(max_level-min_level+1)*proc+min_level] - i );
	
	      MPI_Recv( cellvars, num_out_vars*page_count, MPI_FLOAT, proc, page, MPI_COMM_WORLD, &status );
	      fwrite( cellvars, sizeof(float), num_out_vars*page_count, output );
	      i += page_count;
	      page++;
	    }
	}
      
      fwrite( &size, sizeof(int), 1, output );
    }

  cart_free( cellvars );
}


void write_grid_binary_lower_level_vars(int num_out_vars, int *out_var, FILE *output, int file_parent, int file_num_procs, long *total_cells, int page_size, int *proc_num_cells, int level, int current_level_count, int *current_level)
{
  int i, j, k, m;
  int size;
  float *cellvars;
  int page, page_count;
  int icell, ioct;
  int proc;
  MPI_Status status;

  if(num_out_vars < 1) return;

  cellvars = cart_alloc( num_out_vars * page_size * sizeof(float) );

  /* now write pages of root level hydro variables */
  if ( local_proc_id == file_parent )
    {
      size = num_out_vars * total_cells[level] * sizeof(float);
      fwrite( &size, sizeof(int), 1, output );
    }

  i = 0;
  page = 0;
  while ( i < current_level_count )
    {
      page_count = min( page_size, num_cells_per_level[level] - i*num_children );
      
      j = 0;
      while ( j < num_out_vars*page_count ) { 
	ioct = current_level[i];

	for ( k = 0; k < num_children; k++ ) {
	  icell = oct_child( ioct, k );

	  for ( m = 0; m < num_out_vars; m++ )
	    {
	      cellvars[j++] = cell_var(icell,out_var[m]);
	    }
	  
	}
	i++;
      }

      if ( local_proc_id == file_parent )
	{
	  fwrite( cellvars, sizeof(float), num_out_vars*page_count, output );
	}
      else
	{
	  MPI_Send( cellvars, num_out_vars*page_count, MPI_FLOAT, file_parent, page, MPI_COMM_WORLD );
	  page++;
	}
    }

  if ( local_proc_id == file_parent )
    {
      for ( proc = local_proc_id+1; proc < local_proc_id+file_num_procs; proc++ )
	{
	  i = 0;
	  page = 0;
	  while ( i < proc_num_cells[(max_level-min_level+1)*proc+level] )
	    {
	      page_count = min( page_size, proc_num_cells[(max_level-min_level+1)*proc+level] - i );
	      MPI_Recv( cellvars, num_out_vars*page_count, MPI_FLOAT, proc, page, MPI_COMM_WORLD, &status );
	      fwrite( cellvars, sizeof(float), num_out_vars*page_count, output );
	      i += page_count;
	      page++;
	    }
	}
      
      fwrite( &size, sizeof(int), 1, output );
    }

  cart_free( cellvars );
}


/* two helpers */
void read_grid_binary_top_level_vars(int num_out_vars, int *out_var, FILE *input, int endian, int file_parent, int file_index, int local_file_root_cells, int page_size, int *proc_num_cells, long *proc_cell_index, int *file_sfc_index);
void read_grid_binary_lower_level_vars(int num_out_vars, int *out_var, FILE *input, int endian, int file_parent, int file_index, long *total_cells, int page_size, int *proc_num_cells, int level, long *first_page_count, long *proc_first_index, long *proc_cell_index, int *current_level);


void read_grid_binary2( char *filename ) {
        int i, j, k;
	int size;
	int num_read, flag;
	FILE *input;
	char job[256];
	int minlevel, maxlevel;
	double t, dt;
	float adum, ainit;
        float boxh, Om0, Oml0, Omb0, h;
	int nextras;
	float extra[10];
	char lextra[10][256];
	int *cellrefined[MAX_PROCS], *cellrefinedbuffer;
	int *order;
	int *current_level;
	int endian;
	int proc;
	int ret;
	int level;
	int icell, ioct;
	int cell_counts;
	int page_size;
	int ncell0, sfc_order;
	long current_level_count, current_read_count;
	long next_level_count;
	int file_index, file_parent;
	int local_file_root_cells;
	int page_count;
	long count, start;
	long total_cellrefined;
	long total_cells[max_level-min_level+1];
	int file_parent_proc[MAX_PROCS];
	int proc_num_cells[MAX_PROCS];
	int proc_cur_cells[MAX_PROCS];
	int proc_page_count[MAX_PROCS];
	long first_page_count[MAX_PROCS];
	long proc_next_level_octs[MAX_PROCS];
	long proc_cell_index[MAX_PROCS];
	long proc_first_index[MAX_PROCS+1];
	long next_proc_index[MAX_PROCS+1];
	int file_sfc_index[MAX_PROCS+1];
	char parallel_filename[256];
	int num_requests, num_send_requests;
	int continue_reading, ready_to_read;
	MPI_Request requests[2*MAX_PROCS];
	MPI_Request send_requests[MAX_PROCS];
	MPI_Status status;

        int hydro_vars[num_hydro_vars];
        int num_other_vars = 0;
        int *other_vars = 0;


        /*
        // Maintain the same order as in a previous version
        */
        hydro_vars[0] = HVAR_GAS_DENSITY;
        hydro_vars[1] = HVAR_GAS_ENERGY;
        hydro_vars[2] = HVAR_MOMENTUM + 0;
        hydro_vars[3] = HVAR_MOMENTUM + 1;
        hydro_vars[4] = HVAR_MOMENTUM + 2;
        hydro_vars[5] = HVAR_PRESSURE;
        hydro_vars[6] = HVAR_GAMMA;
        hydro_vars[7] = HVAR_INTERNAL_ENERGY;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	hydro_vars[8] = HVAR_ELECTRON_INTERNAL_ENERGY;
#endif
#ifdef ADVECT_SPECIES
	for(j=0; j<num_chem_species; j++) {
		hydro_vars[8+j] = HVAR_ADVECTED_VARIABLES+j;
	}
#endif /* ADVECT_SPECIES */

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
#ifdef GRAVITY
        num_other_vars += 2;
#endif
#ifdef RADIATIVE_TRANSFER
        num_other_vars += rt_num_disk_vars;
#endif

        other_vars = cart_alloc(num_other_vars*sizeof(int));

#ifdef GRAVITY
        other_vars[0] = VAR_POTENTIAL;
        other_vars[1] = VAR_POTENTIAL_HYDRO;
        k = 2;
#else
        k = 0;
#endif
#ifdef RADIATIVE_TRANSFER
        for(j=0; j<rt_num_disk_vars; j++) other_vars[k+j] = rt_disk_offset + j; 
#endif
#endif /* defined(GRAVITY) || defined(RADIATIVE_TRANSFER) */


	cart_assert( num_output_files >= 1 && num_output_files <= num_procs );

	page_size = num_grid*num_grid;

	/* set up global file information */
	proc = 0;	
	for ( i = 0; i < num_output_files; i++ ) {
		while ( proc*num_output_files / num_procs != i ) {
			proc++;
		}
		file_parent_proc[i] = proc;
	}

	/* determine parallel output options */
	file_index = local_proc_id * num_output_files / num_procs;
	file_parent = file_parent_proc[file_index];

	/* open file handle if parent of parallel file */
	if ( local_proc_id == file_parent ) {
		if ( num_output_files == 1 ) {
			input = fopen(filename,"r");
		} else {
			sprintf( parallel_filename, "%s.%03u", filename, file_index );
			input = fopen(parallel_filename, "r");
		}

		if ( input == NULL ) {
			cart_error( "read_grid_binary2: unable to open file %s for reading!", filename );
		}
	}

	/* the header exists only in the first file */
	if ( local_proc_id == MASTER_NODE ) {
                fread(&size, sizeof(int), 1, input );
		endian = 0;
		if ( size != 256 ) {
			reorder( (char *)&size, sizeof(int) );
			if ( size != 256 ) {
				cart_error("Error: file %s is corrupted", filename );
			} else {
				endian = 1;
				cart_debug("Reordering bytes (file endianness is opposite program)");
			}
		}

                fread(&job, sizeof(char), 256, input );
                fread(&size, sizeof(int), 1, input );

		/* istep, t, dt, adum, ainit */
		fread( &size, sizeof(int), 1, input );
		fread( &step, sizeof(int), 1, input );
		fread( &t, sizeof(double), 1, input );
		fread( &dt, sizeof(double), 1, input );
		fread( &adum, sizeof(float), 1, input );
		fread( &ainit, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&step, sizeof(int) );
			reorder( (char *)&ainit, sizeof(float) );
		}

		a_init = ainit;

#ifndef COSMOLOGY
		if ( endian ) {
			reorder( (char *)&adum, sizeof(float) );
		}
		for(i=min_level; i<=max_level; i++) aexp[i] = adum;
#endif

                /* boxh, Om0, Oml0, Omb0, hubble */
		fread( &size, sizeof(int), 1, input );
		fread( &boxh, sizeof(float), 1, input );
		fread( &Om0, sizeof(float), 1, input );
		fread( &Oml0, sizeof(float), 1, input );
		fread( &Omb0, sizeof(float), 1, input );
		fread( &h, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&boxh, sizeof(float) );
			reorder( (char *)&Om0, sizeof(float) );
			reorder( (char *)&Oml0, sizeof(float) );
			reorder( (char *)&Omb0, sizeof(float) );
			reorder( (char *)&h, sizeof(float) );
		}

		Lbox = boxh;
		Omega0 = Om0;
		OmegaL0 = Oml0;
		Omegab0 = Omb0;
		hubble = h;

		/* nextra (no evidence extras are used...) extra lextra */
		fread( &size, sizeof(int), 1, input );
		fread( &nextras, sizeof(int), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&nextras, sizeof(int) );
		}

		/* extra */
		fread( &size, sizeof(int), 1, input );
		fread( extra, sizeof(float), nextras, input );
		fread( &size, sizeof(int), 1, input );

		/* lextra */
		fread( &size, sizeof(int), 1, input );
		fread( lextra, 256*sizeof(char), nextras, input );
		fread( &size, sizeof(int), 1, input );

		/* Minlevel, MaxLevelNow */
		fread( &size, sizeof(int), 1, input );
		fread( &minlevel, sizeof(int), 1, input );
		fread( &maxlevel, sizeof(int), 1, input);

		if ( endian ) {
			reorder( (char *)&minlevel, sizeof(int) );
			reorder( (char *)&maxlevel, sizeof(int) );
		}

		if ( maxlevel > max_level ) {
			cart_error("File %s has more levels than compiled program (%u)", filename, maxlevel );
		}

		cart_assert( minlevel == min_level );

		fread( &size, sizeof(int), 1, input );

		/* tl */
		fread( &size, sizeof(int), 1, input );
		fread( &tl, sizeof(double), maxlevel-minlevel+1, input );
		fread( &size, sizeof(int), 1, input);

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&tl[i], sizeof(double) );
			}
		}

		/* dtl */
		fread( &size, sizeof(int), 1, input );
		fread( &dtl, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&dtl[i], sizeof(double) );
			}
		}

		/* tl_old */
		fread( &size, sizeof(int), 1, input );
		fread( &tl_old, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&tl_old[i], sizeof(double) );
			}
		}

		/* dtl_old */
		fread( &size, sizeof(int), 1, input );
		fread( &dtl_old, sizeof(double), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );	

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&dtl_old[i], sizeof(double) );
			}
		}

		/* iSO */
		fread( &size, sizeof(int), 1, input );
		fread( &level_sweep_dir, sizeof(int), maxlevel-minlevel+1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = minlevel; i <= maxlevel; i++ ) {
				reorder( (char *)&level_sweep_dir[i], sizeof(int) );
			}
		}

		/* sfc ordering used */
		fread( &size, sizeof(int), 1, input );
		fread( &sfc_order, sizeof(int), 1, input);

		if ( endian ) {
			reorder( (char *)&sfc_order, sizeof(int) );
		}

		if ( sfc_order != SFC ) {
			cart_error("File has different sfc indexing than program");
		}
		fread( &size, sizeof(int), 1, input );

		/* refinement volume */
		fread( &size, sizeof(int), 1, input );
		fread( refinement_volume_min, sizeof(float), nDim, input );
		fread( refinement_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&refinement_volume_min[i], sizeof(float) );
				reorder( (char *)&refinement_volume_max[i], sizeof(float) );
			}
		}

#ifdef STARFORM
		/* star formation volume */
		fread( &size, sizeof(int), 1, input );
		fread( star_formation_volume_min, sizeof(float), nDim, input );
		fread( star_formation_volume_max, sizeof(float), nDim, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < nDim; i++ ) {
				reorder( (char *)&star_formation_volume_min[i], sizeof(float) );
				reorder( (char *)&star_formation_volume_max[i], sizeof(float) );
			}
		}
#endif /* STARFORM */

		/* ncell0 */
		fread( &size, sizeof(int), 1, input );
		fread( &ncell0, sizeof(int), 1, input);
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&ncell0, sizeof(int) );
		}

		if ( ncell0 != (num_grid*num_grid*num_grid) ) {
			cart_error("File has different num_grid than compiled program");
		}
	}

	/* send header information to all other processors */
	MPI_Bcast( &endian, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &minlevel, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &maxlevel, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &step, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &a_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Lbox, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Omega0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &OmegaL0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Omegab0, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &hubble, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( tl, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( dtl, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( dtl_old, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( level_sweep_dir, max_level-min_level+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	MPI_Bcast( refinement_volume_min, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( refinement_volume_max, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );

#ifdef STARFORM
	MPI_Bcast( star_formation_volume_min, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( star_formation_volume_max, nDim, MPI_FLOAT, MASTER_NODE, MPI_COMM_WORLD );
#endif /* STARFORM */

#ifndef COSMOLOGY
	MPI_Bcast( aexp, maxlevel-minlevel+1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
#endif

	if ( local_proc_id == file_parent ) {
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&size, sizeof(int) );
		}
		
		local_file_root_cells = size / sizeof(int);
		cellrefinedbuffer = cart_alloc( local_file_root_cells * sizeof(int) );

		num_read = fread( cellrefinedbuffer, sizeof(int), local_file_root_cells, input );

		if ( num_read != local_file_root_cells ) {
			cart_error("I/O error in read_grid_binary: num_read = %d, local_file_root_cells = %d",
				num_read, local_file_root_cells );
		}

		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < local_file_root_cells; i++ ) {
				reorder( (char *)&cellrefinedbuffer[i], sizeof(int) );
			}
		}

		if ( local_proc_id == MASTER_NODE ) {
			cell_counts = local_file_root_cells;

			file_sfc_index[0] = 0;
			for ( i = 1; i < num_output_files; i++ ) {
				/* need to know how many root cells are in each file */
				/* receive cellrefined */
				file_sfc_index[i] = cell_counts;
				MPI_Recv( &file_sfc_index[i+1], 1, MPI_INT, file_parent_proc[i],
					0, MPI_COMM_WORLD, &status );	
				cell_counts += file_sfc_index[i+1];
			}
			file_sfc_index[num_output_files] = cell_counts;

			cart_assert( cell_counts == num_root_cells );
		} else {
			/* send cellrefined array to MASTER_NODE */
			MPI_Send( &local_file_root_cells, 1, MPI_INT, MASTER_NODE, 0, MPI_COMM_WORLD );
		}
	}

	/* send block information */
	MPI_Bcast( file_sfc_index, num_output_files+1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	/* determine how many cells to expect from each processor */
	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_num_cells[proc] = 0;
	}

	for ( i = 0; i < num_output_files; i++ ) {
		if ( proc_sfc_index[local_proc_id] < file_sfc_index[i+1] &&
				proc_sfc_index[local_proc_id+1] >= file_sfc_index[i] ) {
			proc_num_cells[ file_parent_proc[i] ] = 
					min( proc_sfc_index[local_proc_id+1], file_sfc_index[i+1] ) - 
					max( proc_sfc_index[local_proc_id], file_sfc_index[i] );
		}
	}

	order = cart_alloc( num_cells_per_level[min_level] * sizeof(int) );
	cellrefined[local_proc_id] = cart_alloc( num_cells_per_level[min_level] * sizeof(int) );

	num_requests = 0;

	/* set up non-blocking receives */
	count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( proc_num_cells[proc] > 0 ) { 
			if ( proc == local_proc_id ) {
				/* copy data directly */
				start = max( 0, proc_sfc_index[local_proc_id] - file_sfc_index[file_index] );
				cart_assert( start >= 0 && start <= local_file_root_cells - proc_num_cells[local_proc_id] );
				for ( i = 0; i < proc_num_cells[local_proc_id]; i++ ) {
					cellrefined[local_proc_id][count+i] = cellrefinedbuffer[start+i];
				}
			} else {
				MPI_Irecv( &cellrefined[local_proc_id][count], proc_num_cells[proc], MPI_INT, proc,
					proc_num_cells[proc], MPI_COMM_WORLD, &requests[num_requests++] );
			}

			count += proc_num_cells[proc];
		}
	}

	cart_assert( count == num_cells_per_level[min_level] );

	/* need to avoid deadlocks! senders may also be receivers! */
	if ( local_proc_id == file_parent ) {
		start = 0;
		next_proc_index[proc] = 0;

		/* send all necessary information */
		for ( proc = 0; proc < num_procs; proc++ ) {
			next_proc_index[proc+1] = 0;

			if ( proc == local_proc_id ) { 
				for ( i = 0; i < proc_num_cells[local_proc_id]; i++ ) {
					if ( cellrefinedbuffer[start+i] > 1 ) {
						next_proc_index[proc+1]++;
					}
				}
				start += proc_num_cells[local_proc_id];
			} else if ( start < local_file_root_cells && proc != local_proc_id && 
					file_sfc_index[file_index]+start < proc_sfc_index[proc+1] &&
					file_sfc_index[file_index]+start >= proc_sfc_index[proc] ) {
				count = min( proc_sfc_index[proc+1], file_sfc_index[file_index+1] ) -
						max( file_sfc_index[file_index] + start, proc_sfc_index[proc] ); 
				cart_assert( count > 0 && start + count <= local_file_root_cells );
				for ( i = 0; i < count; i++ ) {
					if ( cellrefinedbuffer[start+i] > 1 ) {
						next_proc_index[proc+1]++;
					}
				}
				MPI_Isend( &cellrefinedbuffer[start], count, MPI_INT, proc, 
					count, MPI_COMM_WORLD, &requests[num_requests++] );
				start += count;
			}
		}

		cart_assert( start == local_file_root_cells );
	}

	/* wait for sends/receives to complete */
	MPI_Waitall( num_requests, requests, MPI_STATUSES_IGNORE );

	if ( local_proc_id == file_parent ) {
		cart_free( cellrefinedbuffer );
	}

	current_level_count = 0;
	next_level_count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		proc_cell_index[proc] = proc_sfc_index[local_proc_id] + current_level_count;
		proc_next_level_octs[proc] = 0;
	
		for ( i = 0; i < proc_num_cells[proc]; i++ ) {
			if ( cellrefined[local_proc_id][current_level_count] > 1 ) {
				ret = split_cell(current_level_count);
				if ( ret ) {
					cart_error("Unable to finish splitting root cells, ran out of octs?");
				}
				cart_assert( next_level_count < num_cells_per_level[min_level] );
				order[next_level_count++] = cell_child_oct[current_level_count];
				proc_next_level_octs[proc]++;
			}
			current_level_count++;
		}
	}

	cart_assert( next_level_count <= num_cells_per_level[min_level] );

	cart_free( cellrefined[local_proc_id] );

	read_grid_binary_top_level_vars(num_hydro_vars,hydro_vars,input,endian,file_parent,file_index,local_file_root_cells,page_size,proc_num_cells,proc_cell_index,file_sfc_index);
#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
	read_grid_binary_top_level_vars(num_other_vars,other_vars,input,endian,file_parent,file_index,local_file_root_cells,page_size,proc_num_cells,proc_cell_index,file_sfc_index);
#endif /* defined(GRAVITY) || defined(RADIATIVE_TRANSFER) */

	/* now read levels */
	for ( level = min_level+1; level <= maxlevel; level++ ) {
		num_requests = 0;

		count = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			proc_num_cells[proc] = proc_next_level_octs[proc]*num_children;

			if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
				MPI_Irecv( &first_page_count[proc], 1, MPI_LONG, proc,
						0, MPI_COMM_WORLD, &requests[num_requests++] );
			} else {
				first_page_count[proc] = 0;
			}

			proc_cell_index[proc] = count;
			count += proc_next_level_octs[proc];
			proc_next_level_octs[proc] = 0;
		}

		if ( local_proc_id == file_parent ) {
			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
			}

			if ( size == sizeof(int) ) {
				fread( &size, sizeof(int), 1, input );

				if ( endian ) {
					reorder( (char *)&size, sizeof(int) );
				}

				total_cells[level] = (long)size;
			} else {
				fread( &total_cells[level], sizeof(long), 1, input );
				if ( endian ) {
					reorder( (char *)&total_cells[level], sizeof(long) );
				}
			}
			fread( &size, sizeof(int), 1, input );
			fread( &size, sizeof(int), 1, input );

			cart_debug("total_cells[%u] = %d", level, total_cells[level] );

			proc_first_index[0] = 0;
			for ( proc = 1; proc <= num_procs; proc++ ) {
				proc_first_index[proc] = (long)num_children*next_proc_index[proc] + 
								proc_first_index[proc-1];
				next_proc_index[proc] = 0;
			}

			cart_assert( proc_first_index[num_procs] == total_cells[level] );

			for ( proc = 0; proc < num_procs; proc++ ) {
				if ( proc_first_index[proc+1]-proc_first_index[proc] > 0 && proc != local_proc_id ) {
					MPI_Isend( &proc_first_index[proc], 1, MPI_LONG, proc,
						0, MPI_COMM_WORLD, &requests[num_requests++] );
				}
			}

			cellrefinedbuffer = cart_alloc( min( total_cells[level], page_size ) * sizeof(int) );
		}

		MPI_Waitall( num_requests, requests, MPI_STATUS_IGNORE );

		num_requests = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			proc_cur_cells[proc] = 0;

			if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) {
				/* set up receive */
				proc_page_count[proc] = min( page_size - first_page_count[proc] % page_size, 
								proc_num_cells[proc] );
				cellrefined[proc] = cart_alloc( min( page_size, proc_num_cells[proc] )*sizeof(float) );
				MPI_Irecv( cellrefined[proc], proc_page_count[proc], MPI_INT, proc,
					proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
				num_requests++;
			} else {
				proc_page_count[proc] = 0;
				requests[proc] = MPI_REQUEST_NULL;
			}
		}

		current_level = order;
		order = cart_alloc( num_children*next_level_count * sizeof(int) );

		current_level_count = 0;
		current_read_count = 0;
		next_level_count = 0;
		continue_reading = 1;
		ready_to_read = 1;
	
		while ( continue_reading ) {
			flag = 0;

			if ( local_proc_id == file_parent && current_read_count < total_cells[level] && ready_to_read ) {
				/* read in a page */
				page_count = min( page_size, total_cells[level] - current_read_count );
				num_read = fread( cellrefinedbuffer, sizeof(int), page_count, input );

				if ( num_read != page_count ) {
					cart_error("I/O Error in read_grid_binary: num_read = %u, page_count = %u",
						num_read, page_count );
				}

				if ( endian ) {
					for ( i = 0; i < page_count; i++ ) {
						reorder( (char *)&cellrefinedbuffer[i], sizeof(int) );
					}
				}

				/* send info to other processors */
				start = 0;
				num_send_requests = 0;

				for ( proc = 0; proc < num_procs; proc++ ) {
					if ( start < page_count && current_read_count + start >= proc_first_index[proc] &&
							current_read_count + start < proc_first_index[proc+1] ) {
						
						count = min( proc_first_index[proc+1], current_read_count+page_count ) -
							max( proc_first_index[proc], current_read_count + start );

						cart_assert( count > 0 && count <= page_count );

						for ( i = 0; i < count; i++ ) {
							cart_assert( start + i < page_count );
							if ( cellrefinedbuffer[start+i] ) {
								next_proc_index[proc+1]++;
							}
						}

						if ( proc == local_proc_id ) {
							/* copy into local buffer */
							proc_page_count[local_proc_id] = count;
							cellrefined[local_proc_id] =
								cart_alloc( proc_page_count[local_proc_id]*sizeof(int));
							for ( i = 0; i < proc_page_count[local_proc_id]; i++ ) {
								cart_assert( start + i < page_count );
								cellrefined[local_proc_id][i] = cellrefinedbuffer[start+i];
							}
						} else {
							MPI_Isend( &cellrefinedbuffer[start], count, MPI_INT, proc,
								count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
						}

						start += count;
					}
				}

				cart_assert( start == page_count );

				current_read_count += page_count;
				ready_to_read = 0;

				/* split local cells first */
				if ( proc_page_count[local_proc_id] > 0 ) {
					flag = 1;
					proc = local_proc_id;
				} else if ( num_requests > 0 ) {
					/* see if we've received anything */
					MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
				}

			} else if ( num_requests > 0 ) {
				/* wait for a receive to complete */
				MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
				num_requests--;
				flag = 1;
			}

			if ( flag == 1 && proc != MPI_UNDEFINED ) {
				cart_assert( proc >= 0 && proc < num_procs );

				/* split cells in received page */
				i = 0;
				while ( i < proc_page_count[proc] ) {
					cart_assert( proc_cell_index[proc] + proc_cur_cells[proc]/num_children 
							< num_cells_per_level[level]/num_children );
					ioct = current_level[ proc_cell_index[proc] + proc_cur_cells[proc]/num_children ]; 

					cart_assert( ioct >= 0 && ioct < num_octs );
					cart_assert( oct_level[ioct] == level );

					for ( j = 0; j < num_children; j++ ) {
						icell = oct_child( ioct, j );
						cart_assert( icell >= 0 && icell < num_cells );
						cart_assert( cell_level(icell) == level );

						if ( cellrefined[proc][i] == 1 ) {
							ret = split_cell(icell);

							if ( ret ) {
								cart_error("Unable to finish splitting cells, ran out of octs?");
							}

							order[ proc_cell_index[proc]*num_children +
								proc_next_level_octs[proc] ] = cell_child_oct[icell];
							next_level_count++;
							proc_next_level_octs[proc]++;

						}

						i++;
						proc_cur_cells[proc]++;
						current_level_count++;
					}
				}

				/* if necessary, start a new receive from that proc */
				if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id ) {
					proc_page_count[proc] = min( page_size, 
						proc_num_cells[proc] - proc_cur_cells[proc] );
					cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
					MPI_Irecv( cellrefined[proc], proc_page_count[proc], MPI_INT,
						proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
					num_requests++;
				} else {
					proc_page_count[proc] = 0;
					cart_free( cellrefined[proc] );
					requests[proc] = MPI_REQUEST_NULL;
				}

				flag = 0;
			}

			/* check if we're done reading */
			if ( local_proc_id == file_parent ) {
				if ( num_requests > 0 ) {
					/* see if sends have completed */
					MPI_Testall( num_send_requests, send_requests, 
							&ready_to_read, MPI_STATUSES_IGNORE );
				} else {
					MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
					num_send_requests = 0;
					ready_to_read = 1;
				}

				if ( current_read_count >= total_cells[level] &&
						current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			} else {
				if ( current_level_count >= num_cells_per_level[level] ) {
					continue_reading = 0;
				}
			}
		}

		if ( local_proc_id == file_parent ) {
			if ( !ready_to_read ) {
				MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
			}

			cart_free( cellrefinedbuffer );

			fread( &size, sizeof(int), 1, input );
		}

		/* repack the order array */
		count = 0;
		for ( proc = 0; proc < num_procs; proc++ ) {
			for ( i = 0; i < proc_next_level_octs[proc]; i++ ) {
				cart_assert( proc_cell_index[proc]*num_children + i  >= count );
				cart_assert( order[ proc_cell_index[proc]*num_children + i ] >= 0 &&
						order[ proc_cell_index[proc]*num_children + i ] < num_octs );
				cart_assert( oct_level[ order[ proc_cell_index[proc]*num_children + i ] ] == level+1 );

				order[count++] = order[ proc_cell_index[proc]*num_children + i ];
			}
		}

		cart_assert( count == next_level_count );

		read_grid_binary_lower_level_vars(num_hydro_vars,hydro_vars,input,endian,file_parent,file_index,total_cells,page_size,proc_num_cells,level,first_page_count,proc_first_index,proc_cell_index,current_level);

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)
		read_grid_binary_lower_level_vars(num_other_vars,other_vars,input,endian,file_parent,file_index,total_cells,page_size,proc_num_cells,level,first_page_count,proc_first_index,proc_cell_index,current_level);
#endif /* defined(GRAVITY) || defined(RADIATIVE_TRANSFER) */

		cart_free( current_level );		
	}

	if ( local_proc_id == file_parent ) {
		fclose( input );
	}

        if(other_vars != 0) cart_free(other_vars);

	cart_free( order );

#ifdef RADIATIVE_TRANSFER
	/* Load RF data */
	rtReadRadiationFieldData(filename,1);
#endif
}


void read_grid_binary_top_level_vars(int num_out_vars, int *out_var, FILE *input, int endian, int file_parent, int file_index, int local_file_root_cells, int page_size, int *proc_num_cells, long *proc_cell_index, int *file_sfc_index)
{
  int i, m;
  int size;
  float *cellvars[MAX_PROCS], *cellvars_buffer;
  int num_requests, num_send_requests;
  int proc;
  int proc_cur_cells[MAX_PROCS];
  int proc_page_count[MAX_PROCS];
  long current_level_count, current_read_count;
  int continue_reading, ready_to_read;
  int num_read, flag;
  int page_count;
  long start, count;
  int icell;
  MPI_Request requests[2*MAX_PROCS];
  MPI_Request send_requests[MAX_PROCS];


  if(num_out_vars < 1) return;


  if ( local_proc_id == file_parent )
    {
      fread( &size, sizeof(int), 1, input );
      cellvars_buffer = cart_alloc( num_out_vars*min( local_file_root_cells, page_size ) * sizeof(float) );
    }

  num_requests = 0;

  for ( proc = 0; proc < num_procs; proc++ )
    {
      proc_cur_cells[proc] = 0;

      if ( proc_num_cells[proc] > 0 && proc != local_proc_id ) 
	{
	  /* set up receive */
	  proc_page_count[proc] = min( proc_num_cells[proc], page_size - ( proc_cell_index[proc] - file_sfc_index[(int)((proc * num_output_files)/ num_procs)]) % page_size);
	  cellvars[proc] = cart_alloc( num_out_vars*min( page_size, proc_num_cells[proc] )*sizeof(float) );
	  MPI_Irecv( cellvars[proc], num_out_vars*proc_page_count[proc], MPI_FLOAT, proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
	  num_requests++;
	} 
      else 
	{
	  proc_page_count[proc] = 0;
	  requests[proc] = MPI_REQUEST_NULL;
	}
    }


  current_level_count = 0;
  current_read_count = 0;
  continue_reading = 1;
  ready_to_read = 1;

  while ( continue_reading )
    {
      flag = 0;

      if ( local_proc_id == file_parent && current_read_count < local_file_root_cells && ready_to_read )
	{
	  /* read in a page */
	  page_count = min( page_size, local_file_root_cells - current_read_count );
	  num_read = fread( cellvars_buffer, sizeof(float), num_out_vars*page_count, input );

	  if ( num_read != num_out_vars*page_count )
	    {
	      cart_error("I/O Error in read_grid_binary: num_read = %u, num_out_vars*page_count = %u", num_read, num_out_vars*page_count );
	    }

	  if ( endian )
	    {
	      for ( i = 0; i < num_out_vars*page_count; i++ )
		{
		  reorder( (char *)&cellvars_buffer[i], sizeof(float) );
		}
	    }

	  /* send info to other processors */
	  start = 0;
	  num_send_requests = 0;
	  for ( proc = 0; proc < num_procs; proc++ )
	    {
	      if ( start < page_count && 
		   file_sfc_index[file_index]+current_read_count+start < proc_sfc_index[proc+1] &&
		   file_sfc_index[file_index]+current_read_count+start >= proc_sfc_index[proc] )
		{
		  
		  count = min( proc_sfc_index[proc+1], file_sfc_index[file_index]+current_read_count+page_count ) - max( proc_sfc_index[proc],  file_sfc_index[file_index] + current_read_count + start );
		  cart_assert( count > 0 && count <= page_count );

		  if ( proc == local_proc_id )
		    {
		      /* copy into local buffer */
		      proc_page_count[local_proc_id] = count;

		      cellvars[local_proc_id] = cart_alloc( num_out_vars*proc_page_count[local_proc_id]*sizeof(float));

		      for ( i = 0; i < num_out_vars*proc_page_count[local_proc_id]; i++ )
			{
			  cart_assert( num_out_vars*start + i < num_out_vars*page_count );
			  cellvars[local_proc_id][i] = cellvars_buffer[num_out_vars*start+i];
			}
		    } 
		  else
		    {
		      MPI_Isend( &cellvars_buffer[num_out_vars*start], num_out_vars*count, MPI_FLOAT, proc, count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
		    }
		  
		  start += count;
		}
	    }

	  cart_assert( start == page_count );
	  current_read_count += page_count;
	  ready_to_read = 0;

	  /* unpack local cells first */
	  if ( proc_page_count[local_proc_id] > 0 )
	    {
	      flag = 1;
	      proc = local_proc_id;
	    } 
	  else if ( num_requests > 0 ) 
	    {
	      /* see if we've received anything */
	      MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
	    }
	} 
      else if ( num_requests > 0 )
	{
	  /* wait for a receive to complete */
	  MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
	  num_requests--;
	  flag = 1;
	}

      if ( flag == 1 && proc != MPI_UNDEFINED )
	{
	  cart_assert( proc >= 0 && proc < num_procs );
			
	  /* unpack received page */
	  i = 0;
	  while ( i < num_out_vars*proc_page_count[proc] ) 
	    {
	      icell = root_cell_location( proc_cell_index[proc] + proc_cur_cells[proc]);
	      cart_assert( icell >= 0 && icell < num_cells_per_level[min_level] );
	      
	      for ( m = 0; m < num_out_vars; m++ )
		{
		  cell_var(icell,out_var[m]) = cellvars[proc][i++];
		}

	      current_level_count++;
	      proc_cur_cells[proc]++;
	    }

	  /* if necessary, start a new receive for that proc */
	  if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id )
	    {
	      proc_page_count[proc] = min( page_size, proc_num_cells[proc] - proc_cur_cells[proc] );
	      cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
	      MPI_Irecv( cellvars[proc], num_out_vars*proc_page_count[proc], MPI_FLOAT, proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
	      num_requests++;
	    } 
	  else
	    {
	      proc_page_count[proc] = 0;
	      cart_free( cellvars[proc] );
	      requests[proc] = MPI_REQUEST_NULL;
	    }
	}

      if ( local_proc_id == file_parent )
	{
	  if ( num_requests > 0 )
	    {
	      /* see if sends have completed, if they have we can continue to read */
	      MPI_Testall( num_send_requests, send_requests, &ready_to_read, MPI_STATUSES_IGNORE );
	    } 
	  else
	    {
	      /* not waiting on any receives, so we can't do anything until
	       * we read in additional data */
	      MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	      num_send_requests = 0;
	      ready_to_read = 1;
	    }

	  if ( current_read_count >= local_file_root_cells && current_level_count >= num_cells_per_level[min_level] )
	    {
	      continue_reading = 0;
	    }
	} 
      else
	{
	  if ( current_level_count >= num_cells_per_level[min_level] )
	    {
	      continue_reading = 0;
	    }
	}
    }

  if ( local_proc_id == file_parent )
    {
      if ( !ready_to_read ) 
	{
	  MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	}

      cart_free( cellvars_buffer );
      fread( &size, sizeof(int), 1, input );

    }
}


void read_grid_binary_lower_level_vars(int num_out_vars, int *out_var, FILE *input, int endian, int file_parent, int file_index, long *total_cells, int page_size, int *proc_num_cells, int level, long *first_page_count, long *proc_first_index, long *proc_cell_index, int *current_level)
{
  int i, j, m;
  int size;
  float *cellvars[MAX_PROCS], *cellvars_buffer;
  int num_requests, num_send_requests;
  int proc;
  int proc_cur_cells[MAX_PROCS];
  int proc_page_count[MAX_PROCS];
  long current_level_count, current_read_count;
  int continue_reading, ready_to_read;
  int num_read, flag;
  int page_count;
  long start, count;
  int icell, ioct;
  MPI_Request requests[2*MAX_PROCS];
  MPI_Request send_requests[MAX_PROCS];


  if(num_out_vars < 1) return;


  if ( local_proc_id == file_parent )
    {
      fread( &size, sizeof(int), 1, input );
      cellvars_buffer = cart_alloc( num_out_vars*min( total_cells[level], page_size ) * sizeof(float) );
    }
	
  num_requests = 0;

  for ( proc = 0; proc < num_procs; proc++ )
    {
      proc_cur_cells[proc] = 0;

      if ( proc_num_cells[proc] > 0 && proc != local_proc_id )
	{
	  /* set up receive */
	  proc_page_count[proc] = min( page_size - first_page_count[proc] % page_size, proc_num_cells[proc] );
	  cellvars[proc] = cart_alloc( num_out_vars*min( page_size, proc_num_cells[proc] )*sizeof(float) );
	  MPI_Irecv( cellvars[proc], num_out_vars*proc_page_count[proc], MPI_FLOAT, proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
	  num_requests++;
	}
      else
	{
	  proc_page_count[proc] = 0;
	  requests[proc] = MPI_REQUEST_NULL;
	}
    }
	
  current_level_count = 0;
  current_read_count = 0;
  continue_reading = 1;
  ready_to_read = 1;

  while ( continue_reading )
    {
      flag = 0;

      if ( local_proc_id == file_parent && current_read_count < total_cells[level] && ready_to_read )
	{
	  /* read in a page */
	  page_count = min( page_size, total_cells[level] - current_read_count );
	  num_read = fread( cellvars_buffer, sizeof(float), num_out_vars*page_count, input );
	
	  if ( num_read != num_out_vars*page_count )
	    {
	      cart_error("I/O Error in read_grid_binary: num_read = %u, num_out_vars*page_count = %u", num_read, num_out_vars*page_count );
	    }
	
	  if ( endian )
	    {
	      for ( i = 0; i < num_out_vars*page_count; i++ )
		{
		  reorder( (char *)&cellvars_buffer[i], sizeof(float) );
		}
	    }
	
	  /* send info to other processors */
	  start = 0;
	  num_send_requests = 0;
	  for ( proc = 0; proc < num_procs; proc++ )
	    {
	      if ( start < page_count && current_read_count + start >= proc_first_index[proc] && current_read_count + start < proc_first_index[proc+1] )
		{
				
		  count = min( proc_first_index[proc+1], current_read_count+page_count ) - max( proc_first_index[proc], current_read_count + start );
		  cart_assert( count > 0 && count <= page_count );

		  if ( proc == local_proc_id )
		    {
		      /* copy into local buffer */
		      proc_page_count[local_proc_id] = count;

		      cellvars[local_proc_id] = cart_alloc( num_out_vars*proc_page_count[local_proc_id]*sizeof(float));
	
		      for ( i = 0; i < num_out_vars*proc_page_count[local_proc_id]; i++ )
			{
			  cart_assert( num_out_vars*start + i < num_out_vars*page_count );
			  cellvars[local_proc_id][i] = cellvars_buffer[num_out_vars*start+i];
			}
		    }
		  else
		    {
		      MPI_Isend( &cellvars_buffer[num_out_vars*start], num_out_vars*count, MPI_FLOAT, proc, count, MPI_COMM_WORLD, &send_requests[num_send_requests++] );
		    }
	
		  start += count;
		}
	    }
	
	  cart_assert( start == page_count );
	  current_read_count += page_count;
	  ready_to_read = 0;
	
	  /* unpack local cells first */
	  if ( proc_page_count[local_proc_id] > 0 )
	    {
	      flag = 1;
	      proc = local_proc_id;
	    }
	  else if ( num_requests > 0 )
	    {
	      /* see if we've received anything */
	      MPI_Testany( num_procs, requests, &proc, &flag, MPI_STATUS_IGNORE );
	    }
	}
      else if ( num_requests > 0 )
	{
	  /* wait for a receive to complete */
	  MPI_Waitany( num_procs, requests, &proc, MPI_STATUS_IGNORE );
	  num_requests--;
	  flag = 1;
	}
	
      if ( flag == 1 && proc != MPI_UNDEFINED )
	{
	  cart_assert( proc >= 0 && proc < num_procs );

	  /* unpack received page */
	  i = 0;
	
	  while ( i < num_out_vars*proc_page_count[proc] )
	    {
	      ioct = current_level[ proc_cell_index[proc] + proc_cur_cells[proc]/num_children];
	      cart_assert( ioct >= 0 && ioct < num_octs );

	      for ( j = 0; j < num_children; j++ )
		{	
		  icell = oct_child( ioct, j );
		  cart_assert( icell >= 0 && icell < num_cells );
		  cart_assert( cell_level(icell) == level );

		  for ( m = 0; m < num_out_vars; m++ )
		    {
		      cell_var(icell,out_var[m]) = cellvars[proc][i++];
		    }

		  
		  current_level_count++;
		  proc_cur_cells[proc]++;
		}
	    }

	  /* if necessary, start a new receive for that proc */
	  if ( proc_cur_cells[proc] < proc_num_cells[proc] && proc != local_proc_id )
	    {
	      proc_page_count[proc] = min( page_size, proc_num_cells[proc] - proc_cur_cells[proc] );
	      cart_assert( proc_page_count[proc] > 0 && proc_page_count[proc] <= page_size );
	      MPI_Irecv( cellvars[proc], num_out_vars*proc_page_count[proc], MPI_FLOAT, proc, proc_page_count[proc], MPI_COMM_WORLD, &requests[proc] );
	      num_requests++;
	    }
	  else
	    {
	      proc_page_count[proc] = 0;
	      cart_free( cellvars[proc] );
	      requests[proc] = MPI_REQUEST_NULL;
	    }
	
	  flag = 0;
	}

      if ( local_proc_id == file_parent )
	{
	  if ( num_requests > 0 )
	    {
	      /* see if sends have completed */
	      MPI_Testall( num_send_requests, send_requests, &ready_to_read, MPI_STATUSES_IGNORE );
	    }
	  else
	    {
	      MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	      num_send_requests = 0;
	      ready_to_read = 1;
	    }
	  
	  if ( current_read_count >= total_cells[level] && current_level_count >= num_cells_per_level[level] )
	    {
	      continue_reading = 0;
	    }
	}
      else
	{
	  if ( current_level_count >= num_cells_per_level[level] )
	    {
	      continue_reading = 0;
	    }
	}
    }

  if ( local_proc_id == file_parent )
    {
      if ( !ready_to_read )
	{
	  MPI_Waitall( num_send_requests, send_requests, MPI_STATUSES_IGNORE );
	}
	
      cart_free( cellvars_buffer );
      fread( &size, sizeof(int), 1, input );
    }
}


#endif /* HYDRO */

