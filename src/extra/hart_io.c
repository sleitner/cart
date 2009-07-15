#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "defs.h"
#include "tree.h"
#include "io.h"
#include "units.h"
#include "timestep.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "sfc.h"
#include "iterators.h"
#include "auxiliary.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "starformation.h"
#include "load_balance.h"
#include "skiplist.h"
#include "index_hash.h"

#ifdef HYDRO

typedef struct {
	int size_start;
	int cell;
	int refined;
	float vars[num_hydro_vars+2];
	int size_end;
} cell_file_struct;

void read_hart_grid_binary( char *filename ) {
	int i, j, k, m;
	int size;
	FILE *input;
	char job[256];
	int minlevel, maxlevel;
	double t, dt;
	float adum, ainit;
        float boxh, Om0, Oml0, Omb0, h;
	int nextras;
	float extra[10];
	char lextra[10][256];
	int refined;
	int *oct_list;
	int  *cellrefined_buffer;
	float *vars_buffer;
	int endian;
	int proc;
	int ret;
	int level;
	int index;
	int next, prev;
	int iNOLL, iHOLL;
	int icell, ioct;
	int cell_counts;
	int page_size;
	int coords[nDim];
	int ncell0, sfc_order;
	cell_file_struct *child_cells;
	index_hash *hart_oct_hash;
	int num_next_level_octs;
	int *next_level_cart_octs, *next_level_hart_octs;
	int iOctFree, nOct;

	/* this function is intended to be run serially */
	cart_assert( num_procs == 1 );

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
		reorder( (char *)&adum, sizeof(float) );
	}

	auni[min_level] = adum;
	auni_init = ainit;

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
	cosmology_set(OmegaM, Om0 );
	cosmology_set(OmegaL, Oml0 );
	cosmology_set(OmegaB, Omb0 );
	cosmology_set( h, h );

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
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&minlevel, sizeof(int) );
		reorder( (char *)&maxlevel, sizeof(int) );
	}

	if ( maxlevel > max_level ) {
		cart_error("File %s has more levels than compiled program (%u)", filename, maxlevel );
	}

	cart_assert( minlevel == min_level );

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

	proc_sfc_index[0] = 0;
	proc_sfc_index[1] = ncell0;

	init_tree();
	init_units();

	/* create hash for hart oct index -> our oct index */
	hart_oct_hash = index_hash_create( ncell0 );

	cellrefined_buffer = cart_alloc( int, num_grid*num_grid );

	fread( &size, sizeof(int), 1, input );
	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
		fread( cellrefined_buffer, sizeof(int), num_grid*num_grid, input );

		if ( endian ) {
			for ( i = 0; i < num_grid*num_grid; i++ ) {
				reorder( (char *)&cellrefined_buffer[i], sizeof(int) );
			}
		}
		i = 0;
	
		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				index = sfc_index( coords );
				icell = root_cell_location( index );
	
				if ( cellrefined_buffer[i] ) {
					cart_assert( cell_is_leaf( icell ) );
					ret = split_cell( icell );
					if ( ret ) {
						cart_error("Error splitting cell %d\n", icell );
					}
					index_hash_add( hart_oct_hash, cellrefined_buffer[i], cell_child_oct[icell] );
				}

				i++;
			}
		}
	}
	fread( &size, sizeof(int), 1, input );

	cart_free( cellrefined_buffer );

#ifdef HYDRO
	vars_buffer = cart_alloc( float, num_grid*num_grid*num_hydro_vars );

	fread( &size, sizeof(int), 1, input );
	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
		fread( vars_buffer, sizeof(float), num_hydro_vars*num_grid*num_grid, input );

		if ( endian ) {
			for ( i = 0; i < num_hydro_vars*num_grid*num_grid; i++ ) {
				reorder( (char *)&vars_buffer[i], sizeof(float) );
			}
		}
		i = 0;

		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				index = sfc_index( coords );
				icell = root_cell_location( index );

				cell_gas_density(icell) = vars_buffer[i++];
				cell_gas_energy(icell) = vars_buffer[i++];
				cell_momentum(icell,0) = vars_buffer[i++];
				cell_momentum(icell,1) = vars_buffer[i++];
				cell_momentum(icell,2) = vars_buffer[i++];
				cell_gas_pressure(icell) = vars_buffer[i++];
				cell_gas_gamma(icell) = vars_buffer[i++];
				cell_gas_internal_energy(icell) = vars_buffer[i++];

#ifdef ADVECT_SPECIES
				for ( m = 0; m < num_chem_species; m++ ) {
					cell_advected_variable(icell,m) = vars_buffer[i++];
				}
#endif /* ADVECT_SPECIES */
			}
		}
	}
	fread( &size, sizeof(int), 1, input );

	cart_free( vars_buffer );
#endif /* HYDRO */

#ifdef GRAVITY 
	vars_buffer = cart_alloc( float, 2*num_grid*num_grid );

	fread( &size, sizeof(int), 1, input );
	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
		fread( vars_buffer, sizeof(float), 2*num_grid*num_grid, input );

		if ( endian ) {
			for ( i = 0; i < 2*num_grid*num_grid; i++ ) {
				reorder( (char *)&vars_buffer[i], sizeof(float) );
			}
		}
		i = 0;

		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				index = sfc_index( coords );
				icell = root_cell_location( index );
				cell_potential(icell) = vars_buffer[i++];
				cell_potential_hydro(icell) = vars_buffer[i++];
                        }
                }
        }
	fread( &size, sizeof(int), 1, input );

	cart_free( vars_buffer );
#endif /* GRAVITY */

	/* iOctfree and nOct */
	fread( &size, sizeof(int), 1, input );
	fread( &iOctFree, sizeof(int), 1, input );
	fread( &nOct, sizeof(int), 1, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&iOctFree, sizeof(int) );
		reorder( (char *)&nOct, sizeof(int) );
	}

	cart_debug("iOctFree = %d", iOctFree );
	cart_debug("nOct = %d", nOct );

	child_cells = cart_alloc( cell_file_struct, num_children );

	for ( level = 1; level <= maxlevel; level++ ) {
		fread( &size, sizeof(int), 1, input );
		if ( endian ) {
			reorder( (char *)&size, sizeof(int) );
		}

		cart_assert( size == 3*sizeof(int) );

		fread( &size, sizeof(int), 1, input );
		fread( &iNOLL, sizeof(int), 1, input );
		fread( &iHOLL, sizeof(int), 1, input );
		
		if ( endian ) {
			reorder( (char *)&iNOLL, sizeof(int) );
			reorder( (char *)&iHOLL, sizeof(int) );
			reorder( (char *)&size, sizeof(int) );
		}

		cart_debug("level %d, %d octs", level, iNOLL );

		fread( &size, sizeof(int), 1, input );

		i = 0;
		oct_list = cart_alloc( int, iNOLL );

		while ( iHOLL != 0 ) {
			oct_list[i] = index_hash_lookup( hart_oct_hash, iHOLL );
			if ( oct_list[i] == NULL_OCT ) {
				cart_debug("oct_list[%u] = %d", i, oct_list[i] );
				cart_debug("iHOLL = %d", iHOLL );
			}
			cart_assert( oct_list[i] != NULL_OCT );
			cart_assert( oct_level[oct_list[i]] == level );
			
			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
			}

			cart_assert( size == 13*sizeof(int) );

			/* skip oct information that's not needed */
			fseek( input, 11*sizeof(int), SEEK_CUR );
			fread( &next, sizeof(int), 1, input );
			fread( &prev, sizeof(int), 1, input );
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				reorder( (char *)&next, sizeof(int) );
			}

			/* goto next oct */
			iHOLL = next;
			i++;
		}

		if ( i != iNOLL ) {
			cart_debug("i = %d, iNOLL = %d", i, iNOLL );
		}
		cart_assert( i == iNOLL );

		index_hash_free( hart_oct_hash );
		next_level_hart_octs = cart_alloc( int, num_children*iNOLL );
		next_level_cart_octs = cart_alloc( int, num_children*iNOLL );
		num_next_level_octs = 0;

		for ( i = 0; i < iNOLL; i++ ) {
			fread( child_cells, sizeof(cell_file_struct), num_children, input );

			for ( j = 0; j < num_children; j++ ) {
				icell = oct_child( oct_list[i], j );
				cart_assert( cell_level(icell) == level );
				
				if ( endian ) {
					reorder( (char *)&child_cells[j].refined, sizeof(int) );
				}

			//	if ( child_cells[j].refined > 0 && child_cells[j].refined < nOct ) {
				if ( child_cells[j].refined > 0 ) {
					ret = split_cell( icell );
					if ( ret ) {
						cart_error("Error splitting cell %d", icell );
					}
					next_level_hart_octs[num_next_level_octs] = child_cells[j].refined;
					next_level_cart_octs[num_next_level_octs] = cell_child_oct[icell];
					num_next_level_octs++;
				}

				if ( endian ) {
					for ( k = 0; k < num_hydro_vars+2; k++ ) {
						reorder( (char *)&child_cells[j].vars[k], sizeof(float) );
					}
				}

				k = 0;
                                cell_gas_density(icell) = child_cells[j].vars[k++];
                                cell_gas_energy(icell) = child_cells[j].vars[k++];
                                cell_momentum(icell,0) = child_cells[j].vars[k++];
                                cell_momentum(icell,1) = child_cells[j].vars[k++];
                                cell_momentum(icell,2) = child_cells[j].vars[k++];
                                cell_gas_pressure(icell) = child_cells[j].vars[k++];
                                cell_gas_gamma(icell) = child_cells[j].vars[k++];
                                cell_gas_internal_energy(icell) = child_cells[j].vars[k++];

#ifdef ADVECT_SPECIES
                                for ( m = 0; m < num_chem_species; m++ ) {
                                        cell_advected_variable(icell,m) = child_cells[j].vars[k++];
                                }
#endif /* ADVECT_SPECIES */

				cell_potential(icell) = child_cells[j].vars[k++];
				cell_potential_hydro(icell) = child_cells[j].vars[k++];
			}
		}

		cart_debug("creating list for level %u, %d", level, num_next_level_octs );
/*
		for ( i = 1; i < num_next_level_octs; i++ ) {
			for ( j = 0; j < i; j++ ) {
				if ( next_level_hart_octs[i] == next_level_hart_octs[j] ) {
					cart_debug("%u %u = %d", i, j, next_level_hart_octs[i] );
				}
			}
		}
*/
		hart_oct_hash = index_hash_create( 2*num_next_level_octs );
		index_hash_add_list( hart_oct_hash, num_next_level_octs, 
					next_level_hart_octs, next_level_cart_octs );

		cart_free( next_level_hart_octs );
		cart_free( next_level_cart_octs );
		cart_free( oct_list );
	}

	cart_free( child_cells );
	index_hash_free( hart_oct_hash );

	cart_debug("building cell buffer");
	build_cell_buffer();
	repair_neighbors();
	
	fclose( input );
}

void write_hart_grid_binary( char *filename ) {
	int i, j, m;
	int size;
	FILE *output;
	float adum, ainit;
	float boxh, OmM0, OmL0, OmB0, h100;
	int minlevel, maxlevel;
	int nextras;
	int *cellrefined;
	float *cellhvars;
	float *cellvars;
	int level;
	int coords[nDim];
	int icell, ioct;
	int ncell0;
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
        cellrefined = cart_alloc(int, page_size );
        cellhvars = cart_alloc(float, num_hydro_vars * page_size );
        cellvars = cart_alloc(float, 2 * page_size );

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
	adum = auni[min_level];
	ainit = auni_init;
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
	OmM0 = cosmology->OmegaM;
	OmL0 = cosmology->OmegaL;
	OmB0 = cosmology->OmegaB;
	h100 = cosmology->h;
	size = 5*sizeof(float);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( &boxh, sizeof(float), 1, output );
	fwrite( &OmM0, sizeof(float), 1, output );
	fwrite( &OmL0, sizeof(float), 1, output );
	fwrite( &OmB0, sizeof(float), 1, output );
	fwrite( &h100, sizeof(float), 1, output );
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
        fwrite( &tl_old, sizeof(double), maxlevel-minlevel+1, output);
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

#endif  /* HYDRO */