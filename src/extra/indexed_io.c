#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "defs.h"
#include "timestep.h"
#include "tree.h"
#include "io.h"
#include "units.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "sfc.h"
#include "iterators.h"
#include "auxiliary.h"
#include "hydro.h"
#include "skiplist.h"
#include "index_hash.h"

#ifdef HYDRO

void read_indexed_grid( char *filename, int num_sfcs, int *sfc_list, int max_level_to_read ) {
	int i, j, k;
	int size;
	FILE *input;
	char job[256];
	int minlevel, maxlevel;
	double t, dt;
	float adum, ainit;
	float boxh, OmM0, OmL0, OmB0, h100;
	int nextras;
	float extra[10];
	char lextra[10][256];
	int *order;
	int endian;
	int proc;
	int ret;
	int level;
	int last_sfc;
	int ncell0, sfc_order;
	long current_level_count;
	long next_level_count;
	long *root_cell_index;
	int *cell_refined;
	float *vars;
	int *cells_per_level;
	int num_cells_unpacked, num_refined_unpacked;
	char varname[32];
	int *varindex;
	long celloffset;
	long total_cells, total_refined;
	int cells_to_read, cell_refined_to_read;
	int num_file_vars;

	struct CosmologyParameters temp_cosmo;

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

		auni_init = ainit;

                /* boxh, Om0, Oml0, Omb0, hubble */
		fread( &size, sizeof(int), 1, input );
		fread( &boxh, sizeof(float), 1, input );
		fread( &OmM0, sizeof(float), 1, input );
		fread( &OmL0, sizeof(float), 1, input );
		fread( &OmB0, sizeof(float), 1, input );
		fread( &h100, sizeof(float), 1, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			reorder( (char *)&boxh, sizeof(float) );
			reorder( (char *)&OmM0, sizeof(float) );
			reorder( (char *)&OmL0, sizeof(float) );
			reorder( (char *)&OmB0, sizeof(float) );
			reorder( (char *)&h100, sizeof(float) );
		}

		Lbox = boxh;
		cosmology_set(OmegaM,OmM0);
		cosmology_set(OmegaL,OmL0);
		cosmology_set(OmegaB,OmB0);
		cosmology_set(h,h100);
		temp_cosmo = *cosmology;

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

		varindex = cart_alloc(int, num_file_vars );

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
	MPI_Bcast( &auni_init, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Bcast( &Lbox, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );

	MPI_Bcast( (char *)&temp_cosmo, sizeof(struct CosmologyParameters), MPI_BYTE, MASTER_NODE, MPI_COMM_WORLD );
	cosmology_copy(&temp_cosmo);

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
		varindex = cart_alloc(int, num_file_vars );
	}

	MPI_Bcast( varindex, num_file_vars, MPI_INT, MASTER_NODE, MPI_COMM_WORLD );

	cells_per_level = cart_alloc(int, (maxlevel-minlevel+1) );

	if ( local_proc_id == MASTER_NODE ) {
		/* pick out indices we need */
		fread( &size, sizeof(int), 1, input );
		celloffset = ftell(input);

		root_cell_index = cart_alloc(long, num_sfcs );
		
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
					reorder( (char *)&cells_per_level[level], sizeof(int) );
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

			cell_refined = cart_alloc(int, cell_refined_to_read );
			vars = cart_alloc(float, cells_to_read*num_file_vars );

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

				cell_refined = cart_alloc(int, cell_refined_to_read );
				vars = cart_alloc(float, cells_to_read*num_file_vars );

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

			order = cart_alloc(int, cells_to_read );
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

#endif  /* HYDRO */
