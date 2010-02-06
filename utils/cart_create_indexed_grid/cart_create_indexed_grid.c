#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "sfc.h"

/* #define ELECTRON_ION_NONEQUILIBRIUM  */
/* #define ENRICH */

#define num_children	(1<<nDim)
#define num_grav_vars	2


#define min(x,y)        (((x) < (y)) ? (x): (y))
#define max(x,y)        (((x) > (y)) ? (x): (y))

const float cell_delta[num_children][nDim] = {
        #if nDim == 1
                { -0.5 }, { 0.5 }
        #elif nDim == 2
                { -0.5, -0.5 }, { 0.5, -0.5 }, { -0.5, 0.5 }, { 0.5, 0.5 }
        #elif nDim == 3
                { -0.5, -0.5, -0.5 }, {  0.5, -0.5, -0.5 }, { -0.5,  0.5, -0.5 },
                {  0.5,  0.5, -0.5 }, { -0.5, -0.5,  0.5 }, {  0.5, -0.5,  0.5 },
                { -0.5,  0.5,  0.5 }, {  0.5,  0.5,  0.5 }
        #else
                #error "No valid cell_delta for that number of dimensions!"
        #endif
};

const char *variable_name[] = {
	"hydro_gas_density",
	"hydro_gas_energy",
	"hydro_momentum_x",
	"hydro_momentum_y",
	"hydro_momentum_z",
	"hydro_gas_pressure",
	"hydro_gas_gamma",
	"hydro_gas_internal_energy",
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	"hydro_electron_internal_energy",
#endif
#ifdef ENRICH
	"hydro_metallicity_II",
	"hydro_metallicity_Ia"
#endif
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

int main ( int argc, char *argv[] ) {
	int i, j;
	int file_index;
	int num_output_grid_files = 1;
	FILE *input, *output;
	FILE *positions, *newpositions;
	char job[256];
	char filename[256];
	int size;
	long lsize;
	int endian = 0;
	int step;
	double t, dt;
	float adum, ainit;
	float boxh, Om0, Oml0, Omb0, h;
	int nextras;
	float extra[10];
	char lextra[10][256];
	int maxlevel, minlevel;
	double tl[num_refinement_levels], dtl[num_refinement_levels], tl_old[num_refinement_levels], dtl_old[num_refinement_levels];
	int level_sweep_dir[num_refinement_levels];
	int ncell0;
	int count;
	float refinement_volume_min[3], refinement_volume_max[3];
	float star_formation_volume_min[3], star_formation_volume_max[3];
	int level;
	long level_count, next_level_count;
	int num_in_page;
	int local_file_root_cells;
	long total_cells[20];
	int index;
	int num_cells;
	int num_hydro_vars;
	long root_file_index_ptr;
	int current_root_index;
	int *refined;
	float *vars;
	float *cell_vars;
	long *root_file_index;
	int *root_index;
	int *root_level_count;
	int *root_next_level_count;
	int *cells_per_level;
	long int total_local_cells;
	int num_read;
	char varname[32];
	float var_max[10], var_min[10];
	float global_var_max[10], global_var_min[10];
	int sfc;
	int coords[3], min_coords[3], max_coords[3];

	if ( argc != 3 ) {
		fprintf(stderr,"Usage: cart_create_indexed_grid input[:num_input_files] output_grid_file\n");
		exit(1);
	}

	for ( index = 0; index < strlen( argv[1] ); index++ ) {
		if ( argv[1][index] == ':' ) {
			argv[1][index] = '\0';
			num_output_grid_files = atoi( &argv[1][index+1] );
			break;
		}
	}

	output = fopen( argv[2], "w" );
	if ( output == NULL ) {
		fprintf(stderr,"Unable to open %s for writing!\n", argv[2] );
		exit(1);
	}
	
	for ( i = 0; i < 10; i++ ) {
		global_var_max[i] = -1e20;
		global_var_min[i] = 1e20;
	}

	for ( file_index = 0; file_index < num_output_grid_files; file_index++ ) {
		if ( num_output_grid_files == 1 ) {
			sprintf(filename,"%s", argv[1] );
		} else {
			sprintf(filename,"%s.%03u", argv[1], file_index );
		}

		printf("file %u\n", file_index ); fflush(stdout);

		/* open file */
		input = fopen( filename,"r");
		if ( input == NULL && file_index != 11 ) {
			printf( "Unable to open file %s for input!\n", filename);
			exit(1);
		}

		if ( file_index == 0 ) {
			fread(&size, sizeof(int), 1, input );
			if ( size != 256 ) {
				reorder( (char *)&size, sizeof(int) );
				if ( size != 256 ) {
					printf("Error: file %s is corrupted\n", argv[1] );
				} else {
					endian = 1;
					printf("Reordering bytes (file endianness is opposite program)\n");
				}
			}

			fread(&job, sizeof(char), 256, input );
			fread(&size, sizeof(int), 1, input );

			printf("job: %s\n", job );

			size = 256*sizeof(char);
			fwrite(&size, sizeof(int), 1, output );
			fwrite(&job, sizeof(char), 256, output );
			fwrite(&size, sizeof(int), 1, output );

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

			printf("step = %u\n", step );
			printf("t = %e\n", t );
			printf("dt = %e\n", dt );
			printf("adum = %e\n", adum );
			printf("ainit = %e\n", ainit );

			size = sizeof(int) + 2*sizeof(double) + 2*sizeof(float);

			fwrite( &size, sizeof(int), 1, output );
			fwrite( &step, sizeof(int), 1, output );
			fwrite( &t, sizeof(double), 1, output );
			fwrite( &dt, sizeof(double), 1, output );
			fwrite( &adum, sizeof(float), 1, output );
			fwrite( &ainit, sizeof(float), 1, output );
			fwrite( &size, sizeof(int), 1, output );

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

			printf("boxh = %f\n", boxh );
			printf("Om0 = %f\n", Om0 );
			printf("Oml0 = %f\n", Oml0 );
			printf("Omb0 = %e\n", Omb0 );
			printf("h = %f\n", h );

			size = 5*sizeof(float);

			fwrite( &size, sizeof(int), 1, output );
			fwrite( &boxh, sizeof(float), 1, output );
			fwrite( &Om0, sizeof(float), 1, output );
			fwrite( &Oml0, sizeof(float), 1, output );
			fwrite( &Omb0, sizeof(float), 1, output );
			fwrite( &h, sizeof(float), 1, output );
			fwrite( &size, sizeof(int), 1, output );

			/* nextra (no evidence extras are used...) extra lextra */
			fread( &size, sizeof(int), 1, input );
			fread( &nextras, sizeof(int), 1, input );
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				reorder( (char *)&nextras, sizeof(int) );
			}

			size = sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &nextras, sizeof(int), 1, output );
			fwrite( &size, sizeof(int), 1, output );

			/* extra */
			fread( &size, sizeof(int), 1, input );
			fread( extra, sizeof(float), nextras, input );
			fread( &size, sizeof(int), 1, input );

			size = nextras*sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( extra, sizeof(int), nextras, output );
			fwrite( &size, sizeof(int), 1, output );

			/* lextra */
			fread( &size, sizeof(int), 1, input );
			fread( lextra, 256*sizeof(char), nextras, input );
			fread( &size, sizeof(int), 1, input );

			size = 256*sizeof(char)*nextras;
			fwrite( &size, sizeof(int), 1, output );
			fwrite( lextra, 256*sizeof(char), nextras, output );
			fwrite( &size, sizeof(int), 1, output );

			/* Minlevel, MaxLevelNow */
			fread( &size, sizeof(int), 1, input );
			fread( &minlevel, sizeof(int), 1, input );
			fread( &maxlevel, sizeof(int), 1, input);
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				reorder( (char *)&minlevel, sizeof(int) );
				reorder( (char *)&maxlevel, sizeof(int) );
			}

			printf("minlevel = %u\n", minlevel );
			printf("maxlevel = %u\n", maxlevel );

			size = 2 * sizeof(int);
			fwrite(&size, sizeof(int), 1, output );
			fwrite(&minlevel, sizeof(int), 1, output );
			fwrite(&maxlevel, sizeof(int), 1, output );
			fwrite(&size, sizeof(int), 1, output );

			for ( level = minlevel; level <= maxlevel; level++ ) {
				total_cells[level] = 0;
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

			size = (maxlevel-minlevel+1) * sizeof(double);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &tl, sizeof(double), maxlevel-minlevel+1, output);
			fwrite( &size, sizeof(int), 1, output );

			/* dtl */
			fread( &size, sizeof(int), 1, input );
			fread( &dtl, sizeof(double), maxlevel-minlevel+1, input);
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				for ( i = minlevel; i <= maxlevel; i++ ) {
					reorder( (char *)&dtl[i], sizeof(double) );
				}
			}

			size = (maxlevel-minlevel+1) * sizeof(double);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &dtl, sizeof(double), maxlevel-minlevel+1, output);
			fwrite( &size, sizeof(int), 1, output );

			/* tl_old */
			fread( &size, sizeof(int), 1, input );
			fread( &tl_old, sizeof(double), maxlevel-minlevel+1, input);
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				for ( i = minlevel; i <= maxlevel; i++ ) {
					reorder( (char *)&tl_old[i], sizeof(double) );
				}
			}

			size = (maxlevel-minlevel+1) * sizeof(double);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &tl_old, sizeof(double), maxlevel-minlevel+1, output );
			fwrite( &size, sizeof(int), 1, output );

			/* dtl_old */
			fread( &size, sizeof(int), 1, input );
			fread( &dtl_old, sizeof(double), maxlevel-minlevel+1, input);
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				for ( i = minlevel; i <= maxlevel; i++ ) {
					reorder( (char *)&dtl_old[i], sizeof(double) );
				}
			}

			size = (maxlevel-minlevel+1) * sizeof(double);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &dtl_old, sizeof(double), maxlevel-minlevel+1, output);
			fwrite( &size, sizeof(int), 1, output );

			printf("tl dtl tl_old dtl_old\n");
			for ( i = minlevel; i <= maxlevel; i++ ) {
				printf("%f %f %f %f\n", tl[i], dtl[i], tl_old[i], dtl_old[i] );
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

			size = (maxlevel-minlevel+1) * sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &level_sweep_dir, sizeof(int), maxlevel-minlevel+1, output);
			fwrite( &size, sizeof(int), 1, output );

			/* sfc ordering used */
			fread( &size, sizeof(int), 1, input );
			fread( &sfc_order, sizeof(int), 1, input);
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				reorder( (char *)&sfc_order, sizeof(int) );
			}

			printf("sfc_order = %u\n", sfc_order );

			size = sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &sfc_order, sizeof(int), 1, output);
			fwrite( &size, sizeof(int), 1, output );

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

			for ( i = 0; i < 3; i++ ) {
				printf("refinement_volume[%u] = %f to %f\n", i, refinement_volume_min[i], 
					refinement_volume_max[i] );
			}

			size = 2*nDim*sizeof(float);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( refinement_volume_min, sizeof(float), nDim, output );
			fwrite( refinement_volume_max, sizeof(float), nDim, output );
			fwrite( &size, sizeof(int), 1, output );

			/* star formation volume (when included) */
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
			}

			if ( size == 6*sizeof(float) ) {
				printf("reading star formation information...\n");

				fread( star_formation_volume_min, sizeof(float), 3, input );
				fread( star_formation_volume_max, sizeof(float), 3, input );
				fread( &size, sizeof(int), 1, input );

				if ( endian ) {
					for ( i = 0; i < 3; i++ ) {
						reorder( (char *)&star_formation_volume_min, sizeof(float) );
						reorder( (char *)&star_formation_volume_max, sizeof(float) );
					}
				}

				size = 2*nDim*sizeof(float);
				fwrite( &size, sizeof(int), 1, output );
				fwrite( star_formation_volume_min, sizeof(float), nDim, output );
				fwrite( star_formation_volume_max, sizeof(float), nDim, output );
				fwrite( &size, sizeof(int), 1, output );

#ifdef ENRICH
				num_hydro_vars = 5+nDim+2;
#else
				num_hydro_vars = 5+nDim;
#endif

				/* ncell0 */
				fread( &size, sizeof(int), 1, input );
			} else {
				num_hydro_vars = 5+nDim;
			}

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			num_hydro_vars++;
#endif

			printf("num_hydro_vars = %u\n", num_hydro_vars );

			fread( &ncell0, sizeof(int), 1, input);
			fread( &size, sizeof(int), 1, input );

			if ( endian ) {
				reorder( (char *)&ncell0, sizeof(int) );
			}

			printf("num_root_cells = %d\n", ncell0 );

			/* write out variable strings */
			size = sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &num_hydro_vars, sizeof(int), 1, output );
			fwrite( &size, sizeof(int), 1, output );

			size = 32*sizeof(char);
			for ( i = 0; i < num_hydro_vars; i++ ) {
				snprintf( varname, 32, "%s", variable_name[i] );
				fwrite( &size, sizeof(int), 1, output );
				fwrite( varname, sizeof(char), 32, output );
				fwrite( &size, sizeof(int), 1, output );
			}

			size = sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &ncell0, sizeof(int), 1, output );
			fwrite( &size, sizeof(int), 1, output );

			/* figure out num_grid */
			num_grid = 1;
			size = ncell0;
			while ( size >>= 3 ) {
				num_grid <<= 1;
			}

			printf("num_grid = %d\n", num_grid );

			init_sfc();

			current_root_index = 0;
			root_file_index = malloc( ncell0*sizeof(long) );
			if ( root_file_index == NULL ) {
				fprintf(stderr, "Ran out of memory!\n");
				exit(1);
			}

			root_file_index_ptr = ftell(output);

			size = ncell0*sizeof(long);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( root_file_index, sizeof(long), ncell0, output );
			fwrite( &size, sizeof(int), 1, output );

			root_file_index[0] = ftell(output);

			min_coords[0] = min_coords[1] = min_coords[2] = num_grid+1;
			max_coords[0] = max_coords[1] = max_coords[2] = -1;
		}

		for ( i = 0; i < num_hydro_vars; i++ ) {
			var_min[i] = 1e20;
			var_max[i] = -1e20;
		}

		/* figure out how many total cells */
		fread( &size, sizeof(int), 1, input );
		if ( endian ) {
			reorder( (char *)&size, sizeof(int) );
		}
		local_file_root_cells = size / sizeof(int);
		root_index = malloc( local_file_root_cells * sizeof(int) );
		if ( root_index == NULL ) {
			fprintf(stderr,"Ran out of memory!\n");
			exit(1);
		}
		printf("local_file_root_cells = %d\n", local_file_root_cells );
		fread( root_index, sizeof(int), local_file_root_cells, input );
		fread( &size, sizeof(int), 1, input );
		
		if ( endian ) {
			for ( i = 0; i < local_file_root_cells; i++ ) {
				reorder( (char *)&root_index[i], sizeof(int) );
			}
		}

		root_level_count = malloc( local_file_root_cells*sizeof(int) );
		root_next_level_count = malloc( local_file_root_cells*sizeof(int) );
		cells_per_level = malloc( local_file_root_cells*(maxlevel+1)*sizeof(int) );

		if ( root_level_count == NULL || root_next_level_count == NULL || cells_per_level == NULL ) {
			fprintf(stderr, "Ran out of memory!\n");
			exit(1);
		}

		next_level_count = 0;
		total_local_cells = 0;
		for ( i = 0; i < local_file_root_cells; i++ ) {
			if ( root_index[i] > 1 ) {
				root_next_level_count[i] = 8;
				next_level_count += 8;
			} else {
				root_next_level_count[i] = 0;
			}

			total_local_cells += root_index[i];
			root_index[i] = total_local_cells - root_index[i];
			root_level_count[i] = 1;

			for ( level = minlevel; level <= maxlevel; level++ ) {
				cells_per_level[(maxlevel-minlevel+1)*i+level] = 0;
			}
		}

		printf("total_local_cells = %ld\n", total_local_cells );

		/* read in cells in pages */
		cell_vars = malloc( total_local_cells*num_hydro_vars*sizeof(float) ); 
		refined = malloc( total_local_cells*sizeof(int) );
		if ( cell_vars == NULL || refined == NULL ) {
			fprintf(stderr,"Ran out of memory!\n");
			exit(1);
		}

		for ( i = 0; i < local_file_root_cells; i++ ) {
			if ( root_next_level_count[i] > 0 ) {
				refined[root_index[i]] = 1;
			} else {
				refined[root_index[i]] = 0;
			}
		}

		level_count = local_file_root_cells;

		for ( level = minlevel; level <= maxlevel && level_count > 0; level++ ) {
			printf("level = %u, level_count = %ld, next_level_count = %ld\n", 
				level, level_count, next_level_count ); fflush(stdout);

			total_cells[level] += level_count;

			/* read in hydro vars */
			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
			}

			count = 0;
			for ( index = 0; index < local_file_root_cells; index++ ) {
				cells_per_level[(maxlevel-minlevel+1)*index + level] = root_level_count[index];

				count += root_level_count[index]*num_hydro_vars;

				/* read in each root cell's vars */
				fread( &cell_vars[ num_hydro_vars*root_index[index] ], 
					sizeof(float), num_hydro_vars*root_level_count[index], input );

				if ( endian ) {
					for ( i = 0; i < num_hydro_vars*root_level_count[index]; i++ ) {
						reorder( (char *)&cell_vars[ num_hydro_vars*root_index[index]+i], 
								sizeof(float) );
					}
				}
			}

			if ( count*sizeof(float) != size ) {
				printf("size = %d, count = %d, %ld\n", size, count, count*sizeof(float) );
			}

			if ( size > count*sizeof(float) ) {
				fseek( input, size-count*sizeof(float), SEEK_CUR );
			}

			fread( &size, sizeof(int), 1, input );

			/* skip gravity variables */
			fread( &size, sizeof(int), 1, input );
			if ( endian ) {
				reorder( (char *)&size, sizeof(int) );
			}
			if ( size != num_grav_vars*level_count*sizeof(float) ) {
				printf("size = %d, num_grav_vars = %d, level_count = %ld, total = %ld\n",
					size, num_grav_vars, level_count, 
					num_grav_vars*level_count*sizeof(float) );
			}		

			fseek( input, (long)(size), SEEK_CUR );
			fread( &size, sizeof(int), 1, input );

			if ( level < maxlevel && next_level_count > 0 ) {
				/* read next level count */
				fread( &size, sizeof(int), 1, input );
				if ( endian ) {
					reorder( (char *)&size, sizeof(int) );
				}
				fread( &lsize, sizeof(long), 1, input );
				if ( endian ) {
					reorder( (char *)&lsize, sizeof(long) );
				}
				printf("%ld vs %ld\n", lsize, next_level_count );
				fread( &size, sizeof(int), 1, input );

				level_count = next_level_count;
				next_level_count = 0;
				fread( &size, sizeof(int), 1, input );
				if ( endian ) {
					reorder( (char *)&size, sizeof(int) );
				}

				count = 0;
				for ( index = 0; index < local_file_root_cells; index++ ) {
					root_index[index] += root_level_count[index];

					count += root_next_level_count[index];

					if ( root_next_level_count[index] > 0 ) {
						num_read = fread( &refined[root_index[index]], sizeof(int), 
							root_next_level_count[index], input );

						if ( num_read != root_next_level_count[index] ) {
							printf("Error reading from file!\n" );
							exit(1);
						}

						if ( endian ) {
							for ( i = 0; i < root_next_level_count[index]; i++ ) {
								reorder( (char *)&refined[root_index[index]+i], sizeof(int) );
							}
						}
					}

					root_level_count[index] = root_next_level_count[index];
					root_next_level_count[index] = 0;

					for ( i = 0; i < root_level_count[index]; i++ ) {
						if ( refined[root_index[index]+i] > 0 ) {
							root_next_level_count[index] += 8;
							next_level_count += 8;
						}
					}
				}

				if ( size != count*sizeof(int) ) {
					printf("after refined: size = %d, count = %d, %ld\n", size, 
						count, count*sizeof(int) );
					exit(1);
				}

				fread( &size, sizeof(int), 1, input );
			} else {
				level_count = next_level_count;
			}
		}

		fclose(input);

		root_index[0] = 0;

		printf("writing data to grid file...\n"); fflush(stdout);

		for ( index = 0; index < local_file_root_cells; index++ ) {
			
			num_cells = 0;
			for ( level = minlevel; level <= maxlevel; level++ ) {
				num_cells += cells_per_level[(maxlevel-minlevel+1)*index+level];
			}

			if ( index < local_file_root_cells - 1 ) {
				root_index[index+1] = root_index[index]+num_cells;
			}

			size = (maxlevel-minlevel+1)*sizeof(int);
		
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &cells_per_level[(maxlevel-minlevel+1)*index], 
					sizeof(int), maxlevel-minlevel+1, output );
			fwrite( &size, sizeof(int), 1, output );
			
			size = (num_cells - cells_per_level[(maxlevel-minlevel+1)*index+maxlevel])*sizeof(int);
			fwrite( &size, sizeof(int), 1, output );
			fwrite( &refined[root_index[index]],
					sizeof(int),
					num_cells-cells_per_level[(maxlevel-minlevel+1)*index+maxlevel],
					output );
			fwrite( &size, sizeof(int), 1, output );

			size = num_cells*num_hydro_vars*sizeof(float);
			fwrite( &size, sizeof(int), 1, output );

			/* reorder from cell-ordered to variable ordered */
			vars = malloc( num_cells*sizeof(float) );
			if ( vars == NULL ) {
				fprintf(stderr,"Error allocating memory for vars!\n");
				exit(1);
			}

			for ( i = 0; i < num_hydro_vars; i++ ) {
				for ( j = 0; j < num_cells; j++ ) {
					vars[j] = cell_vars[ (root_index[index]+j)*num_hydro_vars + i ];
					var_max[i] = ( vars[j] > var_max[i] ) ? vars[j] : var_max[i];
					var_min[i] = ( vars[j] < var_min[i] ) ? vars[j] : var_min[i];
				}

				fwrite( vars, sizeof(float), num_cells, output );
			}

			free( vars );

			fwrite( &size, sizeof(int), 1, output );
			
			current_root_index++;
			if ( current_root_index < ncell0 ) {
				root_file_index[current_root_index] = ftell(output);
			}
		}

		for ( i = 0; i < num_hydro_vars; i++ ) {
			global_var_max[i] = ( var_max[i] > global_var_max[i] ) ? var_max[i] :
						global_var_max[i];
			global_var_min[i] = ( var_min[i] < global_var_min[i] ) ? var_min[i] :
						global_var_min[i];

			printf("var %u: %s %e %e\n", i, variable_name[i], var_min[i], var_max[i] );
		}

		printf("Done writing, freeing allocated memory...\n"); fflush(stdout);

		free( cell_vars ); 
		free( refined ); 
		free( root_level_count ); 
		free( root_next_level_count );
		free( cells_per_level ); 
		free( root_index ); 
	}

	fseek( output, root_file_index_ptr, SEEK_SET );

	size = ncell0*sizeof(long);
	fwrite( &size, sizeof(int), 1, output );
	fwrite( root_file_index, sizeof(long), ncell0, output );
	fwrite( &size, sizeof(int), 1, output );

	fclose(output);

	for ( i = 0; i < num_hydro_vars; i++ ) {
		printf("var[%u] = %s = %e %e\n", i, variable_name[i], global_var_min[i], global_var_max[i] );
	}

	return 0;
}
