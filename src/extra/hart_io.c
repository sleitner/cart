#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "index_hash.h"
#include "io.h"
#include "iterators.h"
#include "load_balance.h"
#include "parallel.h"
#include "sfc.h"
#include "skiplist.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"
#include "refinement.h"



#ifdef HYDRO

#ifndef GRAVITY
#error  warning: without GRAVITY define, io will be wrong 
#endif

/*#define HART_NVAR_L6N*/
/*#define DEBUG*/

#define CONVERT_FOR_IFRIT

#define HART_rt_num_chem_species  8
#define HART_num_enrichment_species  2

#ifdef RADIATIVE_TRANSFER
#ifndef HART_NVAR_L6N 
const int HART_nvarMax = (2 + (rt_num_et_vars+1+rt_num_disk_vars));
#else
const int HART_nvarMax = (2 + (rt_num_et_vars+rt_num_disk_vars));
#endif /* HART_NVAR_L6N */
#else
const int HART_nvarMax = 2 ; //potential vars
#endif /* RADIATIVE_TRANSFER */

const int HART_num_hydro_vars = (num_hydro_vars + HART_num_enrichment_species - num_enrichment_species + HART_rt_num_chem_species - rt_num_chem_species);


typedef struct {
	int cell;
	int refined;
	float *vars;
} cell_file_struct;

void read_hart_grid_binary( char *filename ) {
	int i, j, k, m, idx;
	int size,size2;
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
	if ( !control_parameter_is_set("jobname") ) {
	  cart_debug("setting jobname to header value");
	  strcpy( jobname, job );
	}
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

	box_size = boxh;
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
	fread( tl, sizeof(double), maxlevel-minlevel+1, input );
	fread( &size, sizeof(int), 1, input);

	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&tl[i], sizeof(double) );
		}
	}

	/* dtl */
	fread( &size, sizeof(int), 1, input );
	fread( dtl, sizeof(double), maxlevel-minlevel+1, input);
	fread( &size, sizeof(int), 1, input );
	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&dtl[i], sizeof(double) );
		}
	}

	/* tl_old */
	fread( &size, sizeof(int), 1, input );
	fread( tl_old, sizeof(double), maxlevel-minlevel+1, input);
	fread( &size, sizeof(int), 1, input );
	
	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&tl_old[i], sizeof(double) );
		}
	}

	/* dtl_old */
	fread( &size, sizeof(int), 1, input );
	fread( dtl_old, sizeof(double), maxlevel-minlevel+1, input);
	fread( &size, sizeof(int), 1, input );	

	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&dtl_old[i], sizeof(double) );
		}
	}

	/* iSO */
	fread( &size, sizeof(int), 1, input );
	fread( level_sweep_dir, sizeof(int), maxlevel-minlevel+1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		for ( i = minlevel; i <= maxlevel; i++ ) {
			reorder( (char *)&level_sweep_dir[i], sizeof(int) );
		}
	}
	
	for(i=0;i<nDim;i++){
#ifdef STARFORM
	  star_formation_volume_min[i] = 0;
	  star_formation_volume_max[i] = num_grid;
#endif
	  refinement_volume_min[i] = 0;
	  refinement_volume_max[i] = num_grid;
	  /* warning: hart format does not retain star_formation_volume or refinement_volume arrays */
	  /*  this will lead to a 48byte difference between reverted files */
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

	/* create hash for hart oct index -> our oct index */
	hart_oct_hash = index_hash_create( ncell0 );

	cellrefined_buffer = cart_alloc( int, num_grid*num_grid );

	fread( &size, sizeof(int), 1, input );

	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
	  cart_debug("0-level: %d/%d hash",coords[0],num_grid);
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
	fread( &size2, sizeof(int), 1, input ); 

	cart_free( cellrefined_buffer );

#ifdef HYDRO
	fread( &size, sizeof(int), 1, input );
	if ( endian ) {
	  reorder( (char *)&size, sizeof(int) );
	}
#ifdef DEBUG
	cart_debug("Record size: %d",size);
#endif

	int num_hydro_vars_file = size/num_grid/num_grid/num_grid/sizeof(float);

	cart_debug("num_chem_species_code: %d",num_chem_species);
	cart_debug("num_hydro_vars_file: %d",num_hydro_vars_file);
	cart_assert(num_hydro_vars_file == HART_num_hydro_vars)

	vars_buffer = cart_alloc( float, num_grid*num_grid*num_hydro_vars_file );

	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
	  cart_debug("0-level: %d/%d",coords[0],num_grid);
		fread( vars_buffer, sizeof(float), num_hydro_vars_file*num_grid*num_grid, input );

		if ( endian ) {
			for ( i = 0; i < num_hydro_vars_file*num_grid*num_grid; i++ ) {
				reorder( (char *)&vars_buffer[i], sizeof(float) );
			}
		}
		idx = 0;

		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				index = sfc_index( coords );
				icell = root_cell_location( index );
				
				i = idx*num_hydro_vars_file;

				cell_gas_density(icell) = vars_buffer[i++];
				cell_gas_energy(icell) = vars_buffer[i++];
				cell_momentum(icell,0) = vars_buffer[i++];
				cell_momentum(icell,1) = vars_buffer[i++];
				cell_momentum(icell,2) = vars_buffer[i++];
				cell_gas_pressure(icell) = vars_buffer[i++];
				cell_gas_gamma(icell) = vars_buffer[i++];
				cell_gas_internal_energy(icell) = vars_buffer[i++];
#if defined(ELECTRON_ION_NONEQUILIBRIUM) && defined(CONVERT_FOR_IFRIT)
				cell_electron_internal_energy(icell) = vars_buffer[i++]; /* not used in ART */
#endif 
#ifdef ENRICH
				cell_gas_metal_density_II(icell) = vars_buffer[i++];
#ifdef ENRICH_SNIa              
				cell_gas_metal_density_Ia(icell) = vars_buffer[i++];
#else 
				i++; //skip -- ART always writes enrich
#endif /*  ENRICH_SNIa  */             
#else
				i++; //skip -- ART always writes enrich
#endif /*  ENRICH */             
#ifdef RADIATIVE_TRANSFER
				cell_HI_density(icell) = vars_buffer[i++];
				cell_HII_density(icell) = vars_buffer[i++];
				cell_HeI_density(icell) = vars_buffer[i++];
				cell_HeII_density(icell) = vars_buffer[i++];
				cell_HeIII_density(icell) = vars_buffer[i++];
				cell_H2_density(icell) = vars_buffer[i++];
				for(m=0; m<HART_rt_num_chem_species - rt_num_chem_species; m++){ 
				  i++; //skip -- ART always writes H2+ and H- density
				}
#endif /* RADIATIVE TRANSFER */
#if defined(BLASTWAVE_FEEDBACK) && defined(CONVERT_FOR_IFRIT)
				  cell_blastwave_time(icell) = vars_buffer[i++]; /* not used in ART */
#endif 
				idx++;

			}
		}
	}
	fread( &size, sizeof(int), 1, input );
#ifdef DEBUG
	if ( endian ) {
	  reorder( (char *)&size, sizeof(int) );
	}
	cart_debug("... matching size: %d",size);
#endif

	cart_free( vars_buffer );
#endif /* HYDRO */

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER) 

	fread( &size, sizeof(int), 1, input );
	if ( endian ) {
	  reorder( (char *)&size, sizeof(int) );
	}
#ifdef DEBUG
	cart_debug("Record size: %d",size);
#endif
	int nvar_file = size/num_grid/num_grid/num_grid/sizeof(float); 

	if(HART_nvarMax != nvar_file){
#ifdef RADIATIVE_TRANSFER
	  cart_debug("HART_nvarMax = (2 + (rt_num_et_vars=%d + 1 + rt_num_fields_per_freq=%d * rt_num_freqs=%d ))", rt_num_et_vars,rt_num_fields_per_freq,rt_num_freqs);
	  cart_debug("def RT_UV ? rt_num_freqs=4 : 3 ", rt_num_et_vars,rt_num_fields_per_freq,rt_num_freqs);
#else
	  cart_debug("HART_nvarMax = 2");
#endif /* RADIATIVE_TRANSFER */
	  cart_debug("you may need to alter the rt defs, define HART_NVAR_L6N, or play with HART_nvarMax manually");
	  cart_error("nvar_file=%d != HART_nvarMax=%d",nvar_file,HART_nvarMax);
	}

	vars_buffer = cart_alloc( float, nvar_file*num_grid*num_grid ); 

	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
		fread( vars_buffer, sizeof(float), nvar_file*num_grid*num_grid, input );

		if ( endian ) {
			for ( i = 0; i < nvar_file*num_grid*num_grid; i++ ) {
				reorder( (char *)&vars_buffer[i], sizeof(float) );
			}
		}
		idx = 0;

		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				index = sfc_index( coords );
				icell = root_cell_location( index );

				i = idx*nvar_file;
#ifdef GRAVITY
				cell_potential(icell) = vars_buffer[i++];
				cell_potential_hydro(icell) = vars_buffer[i++];
#endif
#ifdef RADIATIVE_TRANSFER
				for(j=0; j<rt_num_disk_vars; j++) {
				   cell_vars[icell][j+rt_disk_offset] = vars_buffer[i++];
				}
				for(j=0; j<rt_num_et_vars; j++) {
				  cell_vars[icell][j+rt_et_offset] = vars_buffer[i++];
				}
#ifndef HART_NVAR_L6N 
				cell_rt_source(icell) = vars_buffer[i++];
#endif
#endif
				idx++;
                        }
                }
        }
	fread( &size, sizeof(int), 1, input );
#ifdef DEBUG
	if ( endian ) {
	  reorder( (char *)&size, sizeof(int) );
	}
	cart_debug("... matching size: %d",size);
#endif

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

	child_cells = cart_alloc( cell_file_struct, num_children );
	for(i=0; i<num_children; i++) child_cells[i].vars = cart_alloc(float, num_hydro_vars_file+nvar_file);

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

		fread( &size, sizeof(int), 1, input );

		i = 0;
		oct_list = cart_alloc( int, iNOLL );

		while ( iHOLL != 0 ) {
			oct_list[i] = index_hash_lookup( hart_oct_hash, iHOLL );
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

		cart_assert( i == iNOLL );

		index_hash_free( hart_oct_hash );
		next_level_hart_octs = cart_alloc( int, num_children*iNOLL );
		next_level_cart_octs = cart_alloc( int, num_children*iNOLL );
		num_next_level_octs = 0;

		for ( i = 0; i < iNOLL; i++ ) {
		  for(j=0; j<num_children; j++) {
			fread( &size, sizeof(int), 1, input );
			fread( &child_cells[j].cell, sizeof(int), 1, input );
			fread( &child_cells[j].refined, sizeof(int), 1, input );
			fread( child_cells[j].vars, sizeof(float), num_hydro_vars_file+nvar_file, input );
			fread( &size, sizeof(int), 1, input );
		  }

			for ( j = 0; j < num_children; j++ ) {
				icell = oct_child( oct_list[i], j );
				cart_assert( cell_level(icell) == level );
				
				if ( endian ) {
					reorder( (char *)&child_cells[j].refined, sizeof(int) );
				}

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
					for ( k = 0; k < num_hydro_vars_file+nvar_file; k++ ) {
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
#if defined(ELECTRON_ION_NONEQUILIBRIUM) && defined(CONVERT_FOR_IFRIT)
				cell_electron_internal_energy(icell) = child_cells[j].vars[k++]; /* not used in ART */
#endif 
#ifdef ENRICH
				cell_gas_metal_density_II(icell) = child_cells[j].vars[k++];
#ifdef ENRICH_SNIa              
				cell_gas_metal_density_Ia(icell) = child_cells[j].vars[k++];
#else 
				k++; //skip -- ART always writes enrich
#endif /*  ENRICH_SNIa  */             
#else
				k++; //skip -- ART always writes enrich
#endif /*  ENRICH */             
#ifdef RADIATIVE_TRANSFER
				cell_HI_density(icell) = child_cells[j].vars[k++];
				cell_HII_density(icell) = child_cells[j].vars[k++];
				cell_HeI_density(icell) = child_cells[j].vars[k++];
				cell_HeII_density(icell) = child_cells[j].vars[k++];
				cell_HeIII_density(icell) = child_cells[j].vars[k++];
				cell_H2_density(icell) = child_cells[j].vars[k++];
				for(m=0; m<HART_rt_num_chem_species - rt_num_chem_species; m++){ 
				  k++; //skip -- ART always writes H2+ and H- density
				}
#endif /* RADIATIVE TRANSFER */
#if defined(BLASTWAVE_FEEDBACK) && defined(CONVERT_FOR_IFRIT)
				cell_blastwave_time(icell) = child_cells[j].vars[k++]; /* not used in ART */
#endif 
				cart_assert(k == num_hydro_vars_file);
#ifdef GRAVITY
				cell_potential(icell) = child_cells[j].vars[k++];
				cell_potential_hydro(icell) = child_cells[j].vars[k++];
#endif /* GRAVITY */

#ifdef RADIATIVE_TRANSFER
				for(m=0; m<rt_num_disk_vars; m++) {
				  cell_vars[icell][m+rt_disk_offset] = child_cells[j].vars[k++];
				}
				for(m=0; m<rt_num_et_vars; m++) {
				  cell_vars[icell][m+rt_et_offset] = child_cells[j].vars[k++];
				}
#ifndef HART_NVAR_L6N 
				cell_rt_source(icell) = child_cells[j].vars[k++];
#endif
				cart_assert(k == num_hydro_vars_file+nvar_file);
#endif
			}
		}

		hart_oct_hash = index_hash_create( 2*num_next_level_octs );
		index_hash_add_list( hart_oct_hash, num_next_level_octs, 
					next_level_hart_octs, next_level_cart_octs );

		cart_free( next_level_hart_octs );
		cart_free( next_level_cart_octs );
		cart_free( oct_list );
	}

	fclose(input);

	for(i=0; i<num_children; i++) cart_free(child_cells[i].vars);
	cart_free( child_cells );
	index_hash_free( hart_oct_hash );

	build_cell_buffer();
	repair_neighbors();
}

void write_hart_grid_binary( char *filename ) {
        int i, j, k, m, idx;
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
	int first_oct;

	/* this is a single-processor only function */
	cart_assert( num_procs == 1 );

	/* should replace this with a memory parameters and use it in trade_particle_lists as well */	
	page_size = num_grid*num_grid;

	/* allocate pages for writing */
        cellrefined = cart_alloc( int, page_size );
	
	/* assign cellhvars to the right HART size */
#ifdef RADIATIVE_TRANSFER
	cart_debug("Note ART wrote rt_num_chem_species=8 instead of CART's 6." );
	cart_debug("Adding two zeroed float fields for ART ivarHp ivarHm");
#endif
	cart_debug("HART_num_hydro_vars = %d ; num_hydro_vars = %d",HART_num_hydro_vars,num_hydro_vars);
	cellhvars = cart_alloc( float, HART_num_hydro_vars*page_size );

	/* assign nvars to the right HART size */
	cart_debug("HART_nvarMax = %d ; #potential vars= %d",HART_nvarMax,2);
	cellvars = cart_alloc( float, HART_nvarMax*page_size );

	minlevel = min_level;
	maxlevel = max_level_now();

	first_oct = cell_parent_oct( num_cells_per_level[min_level] + num_buffer_cells[min_level] );

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
	boxh = box_size;
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
	fwrite( tl, sizeof(double), maxlevel-minlevel+1, output);
	fwrite( &size, sizeof(int), 1, output );

	/* dtl */
	fwrite( &size, sizeof(int), 1, output );
	fwrite( dtl, sizeof(double), maxlevel-minlevel+1, output);
	fwrite( &size, sizeof(int), 1, output );

	/* tl_old */
	fwrite( &size, sizeof(int), 1, output );
	fwrite( tl_old, sizeof(double), maxlevel-minlevel+1, output);
	fwrite( &size, sizeof(int), 1, output );

	/* dtl_old */
	fwrite( &size, sizeof(int), 1, output );
	fwrite( dtl_old, sizeof(double), maxlevel-minlevel+1, output);
	fwrite( &size, sizeof(int), 1, output );

	/* iSO */
	size = (maxlevel-minlevel+1) * sizeof(int);

	fwrite( &size, sizeof(int), 1, output );
	fwrite( level_sweep_dir, sizeof(int), maxlevel-minlevel+1, output);
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
				if ( cell_is_refined(icell) ) {
					cellrefined[i++] = cell_child_oct[icell] - first_oct + 1;
				} else {
					cellrefined[i++] = 0;
				}
			}
		}

		fwrite( cellrefined, sizeof(int), page_size, output );
	}

	fwrite( &size, sizeof(int), 1, output );

	/* now write pages of root level hydro variables */
	size = HART_num_hydro_vars * ncell0 * sizeof(float);
	fwrite( &size, sizeof(int), 1, output );

	i = size/sizeof(float)/num_grid;
	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
	        cart_debug("0-level: %d/%d",coords[0],num_grid);
	        cart_assert(i == size/sizeof(float)/num_grid);
		idx = 0 ;
		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				icell = sfc_index( coords );

				i = idx*HART_num_hydro_vars;

				cellhvars[i++] = cell_gas_density(icell);
				cellhvars[i++] = cell_gas_energy(icell);
				cellhvars[i++] = cell_momentum(icell,0);
				cellhvars[i++] = cell_momentum(icell,1);
				cellhvars[i++] = cell_momentum(icell,2);
				cellhvars[i++] = cell_gas_pressure(icell);
				cellhvars[i++] = cell_gas_gamma(icell);
				cellhvars[i++] = cell_gas_internal_energy(icell);
#if defined(ELECTRON_ION_NONEQUILIBRIUM) && defined(CONVERT_FOR_IFRIT)
				cellhvars[i++] = cell_electron_internal_energy(icell); /* not used in ART */
#endif 
#ifdef ENRICH
				cellhvars[i++] = cell_gas_metal_density_II(icell);
#ifdef ENRICH_SNIa              
				cellhvars[i++] = cell_gas_metal_density_Ia(icell);
#else 
				cellhvars[i++] = 0;
#endif /*  ENRICH_SNIa  */             
#else
				cellhvars[i++] = 0;
#endif /*  ENRICH */             
#ifdef RADIATIVE_TRANSFER
				cellhvars[i++] = cell_HI_density(icell);
				cellhvars[i++] = cell_HII_density(icell);
				cellhvars[i++] = cell_HeI_density(icell);
				cellhvars[i++] = cell_HeII_density(icell);
				cellhvars[i++] = cell_HeIII_density(icell);
				cellhvars[i++] = cell_H2_density(icell);
				for(j=0; j<HART_rt_num_chem_species - rt_num_chem_species; j++){
				  cellhvars[i++] = 0; // H2+ density H- density
				}
#endif /* RADIATIVE TRANSFER */
#if defined(BLASTWAVE_FEEDBACK) && defined(CONVERT_FOR_IFRIT)
				cellhvars[i++] = cell_blastwave_time(icell); /* not used in ART */
#endif 
				idx++;
			}
		}
		
		fwrite( cellhvars, sizeof(float), HART_num_hydro_vars*page_size, output );
	}
	fwrite( &size, sizeof(int), 1, output );

	/* write grav and RT variables */
	cart_debug("write grav and RT variables");
	size = HART_nvarMax * ncell0 * sizeof(float);
	fwrite( &size, sizeof(int), 1, output );

	i = size/sizeof(float)/num_grid;
	for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
	        cart_assert(i == size/sizeof(float)/num_grid);
		i = 0;
		for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
			for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
				icell = sfc_index( coords );
#ifdef GRAVITY
				cellvars[i++] = cell_potential(icell);
				cellvars[i++] = cell_potential_hydro(icell);
#endif /* GRAVITY */
#ifdef RADIATIVE_TRANSFER
				for(m=0; m<rt_num_disk_vars; m++) {
				  cellvars[i++] = cell_vars[icell][m+rt_disk_offset];
/* 				  parameter ( ifrOTf   = 3 + 1 ) //Far RF at 0 */
/* 				  parameter ( ifrH1f   = 3 + 2 ) //Far RF at HI */
/* 				  parameter ( ifrG1f   = 3 + 3 ) //Far RF at HeI */
/* 				  parameter ( ifrG2f   = 3 + 4 ) //Far RF at HeII */
/* 				  parameter ( ifrOTn   = 3 + 5 ) //Near RF at 0 */
/* 				  parameter ( ifrH1n   = 3 + 6 ) // */
/* 				  parameter ( ifrG1n   = 3 + 7 ) // */
/* 				  parameter ( ifrG2n   = 3 + 8 ) // */
				}
				for(m=0; m<rt_num_et_vars; m++) {
				  cellvars[i++] = cell_vars[icell][m+rt_et_offset];
/* 				  parameter ( irtET1   = 11 + 1) //Eddington Tensor 11 */
/* 				  parameter ( irtET2   = 11 + 2) //Eddington Tensor 12 */
/* 				  parameter ( irtET3   = 11 + 3) //Eddington Tensor 13 */
/* 				  parameter ( irtET4   = 11 + 4) //Eddington Tensor 22 */
/* 				  parameter ( irtET5   = 11 + 5) //Eddington Tensor 23 */
/* 				  parameter ( irtET6   = 11 + 6) //Eddington Tensor 33 */
				}
#ifndef HART_NVAR_L6N 
				cellvars[i++] = cell_rt_source(icell);
/*				parameter ( irtSor   = 11 + 7) // */
#endif
#endif /* RADIATIVE_TRANSFER */
			}
		}

		fwrite( cellvars, sizeof(float), HART_nvarMax*page_size, output );
	}
	fwrite ( &size, sizeof(int), 1, output );
	cart_debug("done writing RT+grav");

	size = 2 * sizeof(int);
	fwrite( &size, sizeof(int), 1, output );

	nOct = 0;
	iOctFree = 0;
	for ( ioct = 0; ioct < num_octs; ioct++ ) {
		if ( oct_level[ioct] != FREE_OCT_LEVEL ) {
			nOct++;
			iOctFree = ioct-first_oct+2;
		}
	}
	fwrite( &iOctFree, sizeof(int), 1, output );
	fwrite( &nOct, sizeof(int), 1, output );
	fwrite( &size, sizeof(int), 1, output );

	/* then write each level's cells in turn */
	for ( level = minlevel+1; level <= maxlevel; level++ ) {
		/* write size */
		size = 3*sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &level, sizeof(int), 1, output );

		if ( local_oct_list[level] == NULL_OCT ) {
			iHOLL = 0;
		} else {
			iHOLL = local_oct_list[level] - first_oct + 1;
		}
		iNOLL = num_cells_per_level[level] / num_children;

		fwrite( &iNOLL, sizeof(int), 1, output );
		fwrite( &iHOLL, sizeof(int), 1, output );

		fwrite( &size, sizeof(int), 1, output );

		ioct = local_oct_list[level];
		for ( i = 0; i < iNOLL; i++ ) {
			size = 	nDim*sizeof(int) + 		/* pos */
				num_neighbors*sizeof(int) + 	/* neighbors */
				sizeof(int) +			/* parent */
				sizeof(int) +			/* level */
				sizeof(int) +			/* next */
				sizeof(int);			/* prev */

			for ( j = 0; j < nDim; j++ ) {
				pos[j] = (int)((oct_pos[ioct][j]+1.0)*(double)(1<<(maxlevel+1)));
			}

			for ( j = 0; j < num_neighbors; j++ ) {
				/* if neighbor is root cell need to convert to index */
				if ( cell_level(oct_neighbors[ioct][j]) == min_level ) {
					sfc_coords( cell_parent_root_sfc( oct_neighbors[ioct][j] ), coords );
					neighbors[j] = coords[2]+num_grid*(coords[1]+num_grid*coords[0])+1;
				} else {
					neighbors[j] = ( cell_parent_oct( oct_neighbors[ioct][j] ) - first_oct ) * num_children +
						cell_child_number( oct_neighbors[ioct][j] ) + ncell0 + 1;
					//neighbors[j] = oct_neighbors[ioct][j] + 1;
				}
			}

			if ( oct_level[ioct] == min_level+1 ) {
				sfc_coords( oct_parent_root_sfc[ioct], coords );
				parent_cell = coords[2]+num_grid*(coords[1]+num_grid*coords[0])+1;
			} else {
				parent_cell = ( cell_parent_oct( oct_parent_cell[ioct] ) - first_oct ) * num_children + 
					cell_child_number( oct_parent_cell[ioct] ) + ncell0 + 1;
				//parent_cell = oct_parent_cell[ioct] + 1;
			}

			fwrite( &size, sizeof(int), 1, output );
			fwrite( pos, sizeof(int), nDim, output );
			fwrite( neighbors, sizeof(int), num_neighbors, output );
			fwrite( &parent_cell, sizeof(int), 1, output );
			fwrite( &oct_level[ioct], sizeof(int), 1, output );

			if ( oct_next[ioct] == NULL_OCT ) {
				inext = 0;
			} else {
				inext = oct_next[ioct] - first_oct + 1;
			}

			if ( oct_prev[ioct] == NULL_OCT ) {
				iprev = 0;
			} else {
				iprev = oct_prev[ioct] - first_oct + 1;
			}

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
					HART_num_hydro_vars*sizeof(float) +	/* hvar */
					HART_nvarMax*sizeof(float);		/* var */

				icellnum = ( ioct - first_oct ) * num_children + cell_child_number( icell ) + ncell0 + 1;
				//icellnum = icell + 1 + ncell0;

				if ( cell_is_refined(icell) ) {
					iOctCh = cell_child_oct[icell] - first_oct + 1;
				} else {
					iOctCh = 0;
				}

				cellhvars[0] = cell_gas_density(icell);
				cellhvars[1] = cell_gas_energy(icell);
				cellhvars[2] = cell_momentum(icell,0);
				cellhvars[3] = cell_momentum(icell,1);
				cellhvars[4] = cell_momentum(icell,2);
				cellhvars[5] = cell_gas_pressure(icell);
				cellhvars[6] = cell_gas_gamma(icell);
				cellhvars[7] = cell_gas_internal_energy(icell);
				k=8;
#if defined(ELECTRON_ION_NONEQUILIBRIUM) && defined(CONVERT_FOR_IFRIT)
				cellhvars[k++] = cell_electron_internal_energy(icell); /* not used in ART */
#endif 
#ifdef ENRICH
				cellhvars[k++] = cell_gas_metal_density_II(icell);
#ifdef ENRICH_SNIa              
				cellhvars[k++] = cell_gas_metal_density_Ia(icell);
#else 
				cellhvars[k++] = 0;
#endif /*  ENRICH_SNIa  */             
#else
				cellhvars[k++] = 0;
#endif /*  ENRICH */             
#ifdef RADIATIVE_TRANSFER
				cellhvars[k++] = cell_HI_density(icell);
				cellhvars[k++] = cell_HII_density(icell);
				cellhvars[k++] = cell_HeI_density(icell);
				cellhvars[k++] = cell_HeII_density(icell);
				cellhvars[k++] = cell_HeIII_density(icell);
				cellhvars[k++] = cell_H2_density(icell);
				for(m=0; m<HART_rt_num_chem_species - rt_num_chem_species; m++){
				  cellhvars[k++] = 0; // H2+ density H- density
				}
#endif /* RADIATIVE TRANSFER */
#if defined(BLASTWAVE_FEEDBACK) && defined(CONVERT_FOR_IFRIT)
				cellhvars[k++] = cell_blastwave_time(icell); /* not used in ART */
#endif 

				k = 0; 
#ifdef GRAVITY
				cellvars[k++] = cell_potential(icell);
				cellvars[k++] = cell_potential_hydro(icell);
#endif /* GRAVITY */
#ifdef RADIATIVE_TRANSFER
				for(m=0; m<rt_num_disk_vars; m++) {
				  cellvars[k++] = cell_vars[icell][m+rt_disk_offset];
				}
				for(m=0; m<rt_num_et_vars; m++) {
				  cellvars[k++] = cell_vars[icell][m+rt_et_offset];
				}
#ifndef HART_NVAR_L6N 
				cellvars[k++] = cell_rt_source(icell);
#endif
#endif

				fwrite( &size, sizeof(int), 1, output );
				fwrite( &icellnum, sizeof(int), 1, output );
				fwrite( &iOctCh, sizeof(int), 1, output );
 				fwrite( cellhvars, sizeof(float), HART_num_hydro_vars, output ); 
 				fwrite( cellvars, sizeof(float), HART_nvarMax, output ); 

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
