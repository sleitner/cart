#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "defs.h"
#include "io.h"
#include "units.h"
#include "sfc.h"
#include "auxiliary.h"
#include "skiplist.h"
#include "healpix.h"
#include "analysis.h"

char output_directory[256];
char logfile_directory[256];
char jobname[256];
int num_output_files = 1;

const double cell_delta[8][3] = {
	{ -0.5, -0.5, -0.5 }, { 0.5, -0.5, -0.5 },
	{ -0.5, 0.5, -0.5 }, { 0.5, 0.5, -0.5 },
	{ -0.5, -0.5, 0.5 }, { 0.5, -0.5, 0.5 },
	{ -0.5, 0.5, 0.5 }, { 0.5, 0.5, 0.5 }
};

void reorder( char *buffer, int size ) {
	int i;
	char tmp;

	for ( i = 0; i < (size/2); i++ ) {
		tmp = buffer[i];
		buffer[i] = buffer[size - i - 1];
		buffer[size - i - 1] = tmp;
	}
}

int compare_indices( const void *a, const void *b ) {
        return *(int *)a - *(int *)b;
}

#ifdef PARTICLES

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

#endif /* PARTICLES */

void read_indexed_grid( char *filename, halo_list *halos, void callback( halo_struct *, cell_struct * ) ) {
	int i, j, k, m, n;
	FILE *input;
	int size;
	int endian;
	int step;
	int sfc;
	int ihalo;
	double aexpn;
	double t, dt;
	float adum, ainit;
	float Lbox, Om0, Oml0, Omb0, h;
	char job[256];
	int nextras;
	float extra[10];
	char lextra[10][256];
	int maxlevel, minlevel;
	double tl[10], dtl[10], tl_old[10], dtl_old[10];
	int level_sweep_dir[10];
	int icell, ncell0;
	int count;
	float refinement_volume_min[3], refinement_volume_max[3];
	float star_formation_volume_min[3], star_formation_volume_max[3];
	long lsize;
	int level;
	long level_count, next_level_count;
	int num_in_page;
	int local_file_root_cells;
	long total_cells;
	int index;
	int num_cells;
	int num_file_vars;
	int current_root_index;
	int *refined;
	float *vars;

	cell_struct cell;
	int density_index;
	int energy_index;
	int pressure_index;
	int internal_energy_index;
	int momentum_x_index;
	int momentum_y_index;
	int momentum_z_index;
	int metallicity_II_index;
	int metallicity_Ia_index;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	int electron_internal_energy_index;
#endif 	

	long *root_index;
	int *root_level_count;
	int *root_next_level_count;
	int cells_per_level[10];
	int *cell_refined;
	int cell_refined_to_read;
	int num_read;
	char varname[32];
	int dim;
	float *pos_level;
	float *pos_next_level;
	int coords[3];
	double ppos[3];
	double axis[3];
	long pixel;
	double r;
	int *sfc_list;

	skiplist *sfc_skiplist;
	
	input = fopen( filename, "r" );
	if ( input == NULL ) {
		cart_error("Unable to open %s for reading!", filename );
	}

	fread(&size, sizeof(int), 1, input );
	endian = 0;
	if ( size != 256*sizeof(char) ) {
		reorder( (char *)&size, sizeof(int) );
		if ( size != 256*sizeof(char) ) {
			cart_error("Error: file %s is corrupted", filename );
		} else {
			endian = 1;
			cart_debug("Reordering bytes (file endianness is opposite program)");
		}
	}

	fread(job, sizeof(char), 256, input );
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
		reorder( (char *)&t, sizeof(double) );
		reorder( (char *)&dt, sizeof(double) );
		reorder( (char *)&adum, sizeof(float) );
		reorder( (char *)&ainit, sizeof(float) );
	}

	/* boxh, Om0, Oml0, Omb0, hubble */
	fread( &size, sizeof(int), 1, input );
	fread( &Lbox, sizeof(float), 1, input );
	fread( &Om0, sizeof(float), 1, input );
	fread( &Oml0, sizeof(float), 1, input );
	fread( &Omb0, sizeof(float), 1, input );
	fread( &h, sizeof(float), 1, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&Lbox, sizeof(float) );
		reorder( (char *)&Om0, sizeof(float) );
		reorder( (char *)&Oml0, sizeof(float) );
		reorder( (char *)&Omb0, sizeof(float) );
		reorder( (char *)&h, sizeof(float) );
	}

	size = 5*sizeof(float);

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
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&sfc_order, sizeof(int) );
	}

	/* refinement volume */
	fread( &size, sizeof(int), 1, input );
	fread( refinement_volume_min, sizeof(float), 3, input );
	fread( refinement_volume_max, sizeof(float), 3, input );
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		for ( i = 0; i < 3; i++ ) {
			reorder( (char *)&refinement_volume_min[i], sizeof(float) );
			reorder( (char *)&refinement_volume_max[i], sizeof(float) );
		}
	}

	/* star formation volume (when included) */
	fread( &size, sizeof(int), 1, input );

	if ( size == 6*sizeof(float) ) {
		fread( star_formation_volume_min, sizeof(float), 3, input );
		fread( star_formation_volume_max, sizeof(float), 3, input );
		fread( &size, sizeof(int), 1, input );

		if ( endian ) {
			for ( i = 0; i < 3; i++ ) {
				reorder( (char *)&star_formation_volume_min[i], sizeof(float) );
				reorder( (char *)&star_formation_volume_max[i], sizeof(float) );
			}
		}

		/* ncell0 */
		fread( &size, sizeof(int), 1, input );
	}

	/* read in variable strings */
	fread( &num_file_vars, sizeof(int), 1, input );
	fread( &size, sizeof(int), 1, input );

	cart_debug("num_file_vars = %d", num_file_vars );

	for ( i = 0; i < num_file_vars; i++ ) {
		fread( &size, sizeof(int), 1, input );
		fread( varname, sizeof(char), 32, input );
		fread( &size, sizeof(int), 1, input );

		cart_debug("file_var[%d] = %s", i, varname );

		if ( strcmp( varname, "hydro_gas_density" ) == 0 ) {
			density_index = i;
		} else if ( strcmp( varname, "hydro_gas_energy" ) == 0 ) {
			energy_index = i;
		} else if ( strcmp( varname, "hydro_momentum_x" ) == 0 ) {
			momentum_x_index = i;
		} else if ( strcmp( varname, "hydro_momentum_y" ) == 0 ) {
			momentum_y_index = i;
		} else if ( strcmp( varname, "hydro_momentum_z" ) == 0 ) {
			momentum_z_index = i;
		} else if ( strcmp( varname, "hydro_gas_pressure" ) == 0 ) {
			pressure_index = i;
		} else if ( strcmp( varname, "hydro_gas_internal_energy" ) == 0 ) {
			internal_energy_index = i;
		} else if ( strcmp( varname, "hydro_metallicity_II" ) == 0 ) {
			metallicity_II_index = i;
		} else if ( strcmp( varname, "hydro_metallicity_Ia" ) == 0 ) {
			metallicity_Ia_index = i;
#ifdef ELECTRON_ION_NONEQUILIBRIUM
		} else if ( strcmp( varname, "hydro_electron_internal_energy" ) == 0 ) {
			electron_internal_energy_index = i;
#endif
		}
	}

	fread( &size, sizeof(int), 1, input );
	fread( &ncell0, sizeof(int), 1, input);
	fread( &size, sizeof(int), 1, input );

	if ( endian ) {
		reorder( (char *)&ncell0, sizeof(int) );
	}

	init_sfc();

	root_index = cart_alloc( ncell0 * sizeof(long) );

	fread( &size, sizeof(int), 1, input );
	fread( root_index, sizeof(long), ncell0, input );
	fread( &size, sizeof(int), 1, input );

	for ( ihalo = 0; ihalo < halos->num_halos; ihalo++ ) {
		sfc_skiplist = skiplist_init();

		/* loop over root level cellsi */
		for ( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
			ppos[0] = ( (double)coords[0]+0.5 - halos->list[ihalo].pos[0] );
			if ( ppos[0] > (double)num_grid/2. ) {
				ppos[0] -= (double)num_grid;
			} else if ( ppos[0] < -(double)num_grid/2. ) {
				ppos[0] += (double)num_grid;
			}

			for ( coords[1] = 0; coords[1] < num_grid; coords[1]++ ) {
				ppos[1] = ( (double)coords[1]+0.5 - halos->list[ihalo].pos[1] );
				if ( ppos[1] > (double)num_grid/2. ) {
					ppos[1] -= (double)num_grid;
				} else if ( ppos[1] < -(double)num_grid/2. ) {
					ppos[1] += (double)num_grid;
				}
			
				for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
					ppos[2] = ( (double)coords[2]+0.5 - halos->list[ihalo].pos[2] );
					if ( ppos[2] > (double)num_grid/2. ) {
						ppos[2] -= (double)num_grid;
					} else if ( ppos[2] < -(double)num_grid/2. ) {
						ppos[2] += (double)num_grid;
					}

					/* -1 takes care of center->corner differences and possible overlap */
					r = sqrt(ppos[0]*ppos[0]+ppos[1]*ppos[1]+ppos[2]*ppos[2]) - 1.0;
#ifdef PROJECTED_PROFILES
					if ( r < max( halos->list[ihalo].extent, rbinmax/r0/1000.0 ) ) {
#else
					if ( r < halos->list[ihalo].extent ) {
#endif
						index = sfc_index(coords);
						skiplist_insert(sfc_skiplist, index);
					}
				}
			}
		}

		/* extract list of unique sfc indices */
		sfc_list = cart_alloc( skiplist_size(sfc_skiplist)*sizeof(int) );

		count = 0;
		skiplist_iterate( sfc_skiplist );
		while ( skiplist_next( sfc_skiplist, &index ) ) {
			sfc_list[count++] = index;
		}

		skiplist_destroy( sfc_skiplist );

		/* probably not necessary.... */
		qsort( sfc_list, count, sizeof(int), compare_indices );

		cart_debug("selected %u root trees...", count );

		for ( sfc = 0; sfc < count; sfc++ ) {
			if ( sfc % 1024 == 0 ) {
				cart_debug("%7d/%7d", sfc, count );
			}
			fseek( input, root_index[sfc_list[sfc]], SEEK_SET );

			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
				if ( size != (maxlevel-minlevel+1)*sizeof(int)) {
					cart_error("Error, file is corrupt!");
				}
			}

			fread( cells_per_level, sizeof(int), maxlevel-minlevel+1, input );
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				for ( level = minlevel; level <= maxlevel; level++ ) {
					reorder( (char *)&cells_per_level[level], sizeof(int) );
				}
			}

			total_cells = 0;
			for ( level = 0; level <= maxlevel; level++ ) {
				total_cells += cells_per_level[level];
			}

			cell_refined_to_read = total_cells-cells_per_level[maxlevel];

			cell_refined = cart_alloc( cell_refined_to_read*sizeof(int) );
			vars = cart_alloc( num_file_vars*total_cells*sizeof(float) );

			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );

				if ( size != cell_refined_to_read*sizeof(int)) {
					cart_error("Error: file is corrupt reading cell_refined!");
				}
			}
			fread( cell_refined, sizeof(int), cell_refined_to_read, input );
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				for ( i = 0; i < cell_refined_to_read; i++ ) {
					reorder( (char *)&cell_refined[i], sizeof(int) );
				}
			}

			fread( &size, sizeof(int), 1, input );
			fread( vars, sizeof(float), num_file_vars*total_cells, input );
			fread( &size, sizeof(int), 1, input );

			pos_level = cart_alloc( nDim*sizeof(float) );

			sfc_coords( sfc_list[sfc], coords );

			for ( dim = 0; dim < nDim; dim++ ) {
				pos_level[dim] = (double)coords[dim]+0.5;
			}

			icell = 0;

			for ( level = minlevel; level <= maxlevel; level++ ) {
				if ( level < maxlevel ) {
					pos_next_level = cart_alloc( nDim*cells_per_level[level+1]*sizeof(float) );
					next_level_count = 0;
				} else {
					pos_next_level = NULL;
				}

				for ( j = 0; j < cells_per_level[level]; j++ ) {
					if ( level == maxlevel || !cell_refined[icell+j] ) {
						cell.level = level;
						for ( k = 0; k < nDim; k++ ) {
							cell.pos[k] = pos_level[nDim*j+k];
						}

						cell.gas_density = vars[density_index*total_cells + icell + j];
						cell.gas_energy = vars[energy_index*total_cells + icell + j];
						cell.gas_pressure = vars[pressure_index*total_cells + icell + j];
						cell.gas_internal_energy = vars[internal_energy_index*total_cells + icell + j];
						cell.momentum[0] = vars[momentum_x_index*total_cells + icell + j];
						cell.momentum[1] = vars[momentum_y_index*total_cells + icell + j];
						cell.momentum[2] = vars[momentum_z_index*total_cells + icell + j];

#ifdef ENRICH
						cell.metallicity_II = vars[metallicity_II_index*total_cells + icell + j];
#ifdef ENRICH_SNIa
						cell.metallicity_Ia = vars[metallicity_Ia_index*total_cells + icell + j];
#endif
#endif

#ifdef ELECTRON_ION_NONEQUILIBRIUM
						cell.electron_internal_energy = vars[electron_internal_energy_index*total_cells + icell + j];
#endif

						/* callback here */
						callback( &halos->list[ihalo], &cell );

					} else if ( cells_per_level[level+1] > 0 ) {
						for ( k = 0; k < 8; k++ ) {
							for ( dim = 0; dim < nDim; dim++ ) {
								pos_next_level[nDim*next_level_count+nDim*k+dim] =
									pos_level[nDim*j+dim] + 0.5*cell_size[level]*cell_delta[k][dim];
							}
						}

						next_level_count += 8;
					}
				}

				icell += cells_per_level[level];

				cart_free(pos_level);
				pos_level = pos_next_level;
			}

			cart_free( cell_refined );
			cart_free( vars );
		}

		cart_free( sfc_list ); 
	}

	cart_free( root_index );
	fclose(input);
}

