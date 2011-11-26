#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "refinement_indicators.h"
#include "tree.h"
#include "viewdump.h"


FILE *output;
float plane;
int dimension;
int lev;
int corner_count;

#if nDim != 3
#error Viewdump designed for 3 dimensions only!
#endif

const int other_directions[nDim][nDim-1] = {
	{ 1, 2 }, { 0, 2 }, { 0, 1 }
};

int is_intersected( int cell ) {
	double position[nDim];

	if ( cell_is_leaf(cell) || cell_level(cell) == lev ) {
		cell_center_position( cell, position );

		return ( fabs( position[dimension] - plane ) <= (cell_size[cell_level(cell)] / 2.0) ); 
	} else {
		return 0;
	}
}

void print_cell_corners( int cell, int level ) {
	double position[nDim];
	float x,y;
	float cellsize;

	if ( is_intersected(cell) ) {
		cell_center_position(cell, position );
		cellsize = cell_size[level] / 2.0;

		x = position[ other_directions[dimension][0] ];
		y = position[ other_directions[dimension][1] ];

		fprintf(output, "%e %e %e %e\n", x - cellsize + 1.0, y - cellsize + 1.0, 0.0, 0.0 );
		fprintf(output, "%e %e %e %e\n", x - cellsize + 1.0, y + cellsize + 1.0, 0.0, 0.0 );
		fprintf(output, "%e %e %e %e\n", x + cellsize + 1.0, y + cellsize + 1.0, 0.0, 0.0 );
		fprintf(output, "%e %e %e %e\n", x + cellsize + 1.0, y - cellsize + 1.0, 0.0, 0.0 );
	}
}

void print_cell_location_vars( int cell, int level ) {
	        double position[nDim];

        if ( is_intersected(cell) ) {
                cell_center_position( cell, position);

#ifdef HYDRO
				fprintf(output, "%u %u %u %u 1 %u %u %u %u %e %e %e\n", corner_count, corner_count+1,
                        corner_count+2, corner_count+3, level, local_proc_id, cell_is_local(cell),
						processor_owner(cell_parent_root_sfc(cell)),
						position[0], position[1], position[2] );
#endif

		corner_count += 4;
		cart_assert(corner_count >= 0 );
	}
}

void print_cell_hvars( int cell, int level ) {
	double position[nDim];

	if ( is_intersected(cell) ) {
		cell_center_position( cell, position);

#ifdef HYDRO
		fprintf(output, "%u %u %u %u 1 %u %u %d %d %d %e %e %e %e %e %e %e %e %e %e\n", corner_count, corner_count+1,
			corner_count+2, corner_count+3, level, local_proc_id, 
			cell, cell_parent_cell(cell), cell_parent_root_sfc(cell),
			position[0], position[1], position[2], cell_gas_density(cell),
		        cell_gas_internal_energy(cell),
		        cell_gas_pressure(cell), cell_gas_gamma(cell), 
			cell_momentum(cell,0), cell_momentum(cell,1), cell_momentum(cell,2) );
#endif

		corner_count += 4;
		cart_assert(corner_count >= 0 );
	}
}

void print_cell_grav_vars( int cell, int level ) {
	double position[nDim];
	
	if ( is_intersected(cell) ) {
		cell_center_position( cell, position );

#ifdef GRAVITY
		fprintf( output, "%u %u %u %u 1 %u %u %u %e %e %e %e %e %e %e %e\n", corner_count, corner_count+1,
			corner_count+2, corner_count+3, level, local_proc_id, cell,
			position[0], position[1], position[2], 
			cell_total_mass(cell), cell_potential(cell), cell_accel(cell,0),
			cell_accel(cell,1), cell_accel(cell,2) );
#endif

		corner_count += 4;
	}
}

void print_cell_refinement_vars( int cell, int level ) {
	int i;
	int neighbors[num_neighbors];
	float drho[nDim];

	if ( is_intersected(cell) ) {
		cell_all_neighbors( cell, neighbors );

#ifdef HYDRO
		for ( i = 0; i < nDim; i++ ) {
			drho[i] = fabs( cell_gas_density( neighbors[2*i] )
					- cell_gas_density( neighbors[2*i+1] ) ) /
				min ( cell_gas_density( neighbors[2*i] ),
					cell_gas_density( neighbors[2*i+1] ) );
		}

/*
		fprintf(output, "%u %u %u %u 1 %u %u %u %f %f %f %f %f %f\n", 
			corner_count, corner_count+1, corner_count+2, corner_count+3,
			level, local_proc_id, cell,
			mass_indicator( cell, level ),
			shock_indicator( cell, level ),
			contact_discontinuity_indicator( cell, level ),
			density_gradient_indicator( cell, level ),
			pressure_gradient_indicator( cell, level ),
			entropy_gradient_indicator( cell, level ) );
*/
#endif

		corner_count += 4;
	}
}

void print_diffusion_vars( int cell, int level ) {
	if ( is_intersected(cell) ) {
		fprintf(output, "%u %u %u %u 1 %u %u %u %f %f\n", corner_count,
			corner_count+1, corner_count+2, corner_count+3,
			level, local_proc_id, cell,
			refinement_indicator(cell,0),
			refinement_indicator(cell,1) );
		corner_count += 4;		
	}
}

void viewdump( const char *filename, int max_level_dumped, float slice, int d, int vars, int cell_type ) {
	int i, j, p;
	int num_cells_intersected;
	int total_intersected;
	int num_level_cells;
	int *level_cells;
	int cell_count[MAX_PROCS];

	lev = max_level_dumped;
	dimension = d;
	plane = slice;
	num_cells_intersected = 0;

	if ( local_proc_id == MASTER_NODE ) {
		cart_debug("viewdump %s", filename );
	}
	
	for ( i = min_level; i <= max_level_dumped; i++ ) {
		select_level( i, cell_type, &num_level_cells, &level_cells );
		for ( j = 0; j < num_level_cells; j++ ) {
			if ( is_intersected( level_cells[j] ) ) {
				num_cells_intersected++;
			}
		}
		cart_free( level_cells );
	}

	MPI_Reduce( &num_cells_intersected, &total_intersected, 1, MPI_INT, MPI_SUM, MASTER_NODE, mpi.comm.run );
	MPI_Allgather( &num_cells_intersected, 1, MPI_INT, cell_count, 1, MPI_INT, mpi.comm.run );

	for ( p = 0; p < num_procs; p++ ) {
		if ( p == local_proc_id ) {
			cart_debug("dumping %u cells", num_cells_intersected );
		} 

		MPI_Barrier(mpi.comm.run );
	}

	if ( local_proc_id == MASTER_NODE ) {
	        output = fopen(filename, "w");
	
		fprintf(output, "testname 0.0 0.0 0\n");
		fprintf(output, "zones      1 1 1\n");

		if ( vars == DUMP_HVARS ) {
			fprintf(output, "15 level processor cellnum parent parent_sfc x y z density energy pressure gamma px py pz\n");
		} else if ( vars == DUMP_GRAV_VARS ) {
			fprintf(output, "11 level processor cellnum x y z density potential ax ay az\n" );
		} else if ( vars == DUMP_REFINEMENT_VARS ) {
			fprintf(output, "9 level processor cellnum mass_ind shock_ind contact_int density_ind pressure_ind entropy_ind\n");
		} else if ( vars == DUMP_DIFFUSION_VARS ) {
			fprintf(output, "5 level processor cellnum refinement1 refinement2\n" );
		} else if ( vars == DUMP_LOCATION_VARS ) {
			fprintf(output, "7 level print_processor cell_is_local processor_owner x y z\n");	
		}

		fprintf(output, "0\n");
		fprintf(output, "%u\n", 4*total_intersected );

		fclose(output);
	}

	MPI_Barrier(mpi.comm.run);

	for ( p = 0; p < num_procs; p++ ) {
		if ( p == local_proc_id ) {
			output = fopen(filename,"a");

			for ( i = min_level; i <= max_level_dumped; i++ ) {
				select_level( i, cell_type, &num_level_cells, &level_cells );
				for ( j = 0; j < num_level_cells; j++ ) {
					print_cell_corners( level_cells[j], i );
				}
				cart_free( level_cells );
			}

			fclose(output);
		}

		MPI_Barrier(mpi.comm.run);
	}

	if ( local_proc_id == MASTER_NODE ) {
		output = fopen(filename, "a");
		fprintf(output, "%u\n", total_intersected );
		fclose(output);
        }

	corner_count = 1;
	for ( p = 0; p < num_procs; p++ ) {
		if ( p == local_proc_id ) {
			output = fopen(filename,"a");

			for ( i = min_level; i <= max_level_dumped; i++ ) {
				select_level( i, cell_type, &num_level_cells, &level_cells );
				if ( vars == DUMP_HVARS ) {
					for ( j = 0; j < num_level_cells; j++ ) {
						print_cell_hvars( level_cells[j], i );
					}
				} else if ( vars == DUMP_GRAV_VARS ) {
					for ( j = 0; j < num_level_cells; j++ ) {
						print_cell_grav_vars( level_cells[j], i );
					}
				} else if ( vars == DUMP_REFINEMENT_VARS ) {
					for ( j = 0; j < num_level_cells; j++ ) {
						print_cell_refinement_vars( level_cells[j], i );
					}
				} else if ( vars == DUMP_DIFFUSION_VARS ) {
					for ( j = 0; j < num_level_cells; j++ ) {
						print_diffusion_vars( level_cells[j], i );
					}
				} else if ( vars == DUMP_LOCATION_VARS ) {
					for ( j = 0; j < num_level_cells; j++ ) {
						print_cell_location_vars( level_cells[j], i );
					}
				}

                                cart_free( level_cells );
			}

			fclose(output);
		} else {
			corner_count += cell_count[p] * 4;
		}

                MPI_Barrier(mpi.comm.run);
        }

	if ( local_proc_id == MASTER_NODE ) {
                output = fopen(filename, "a");
                fprintf(output, "interfaces 0\n");
                fclose(output);
        }
}

void viewdump_serial( const char *filename, int max_level_dumped, float slice, int d, int vars, int cell_type ) {
	int i, j;
	int num_cells_intersected;
	int total_intersected;
	int num_level_cells;
	int *level_cells;

	lev = max_level_dumped;
	dimension = d;
	plane = slice;
	num_cells_intersected = 0;

        for ( i = min_level; i <= max_level_dumped; i++ ) {
                select_level( i, cell_type, &num_level_cells, &level_cells );
                for ( j = 0; j < num_level_cells; j++ ) {
                        if ( is_intersected( level_cells[j] ) ) {
                                num_cells_intersected++;
                        }
                }
                cart_free( level_cells );
        }

	total_intersected = num_cells_intersected;

	output = fopen(filename, "w");

	fprintf(output, "testname 0.0 0.0 0\n");
	fprintf(output, "zones      1 1 1\n");

	if ( vars == DUMP_HVARS ) {
		fprintf(output, "15 level processor cellnum parent parent_sfc x y z density energy pressure gamma px py pz\n");
	} else if ( vars == DUMP_GRAV_VARS ) {
		fprintf(output, "11 level processor cellnum x y z density potential ax ay az\n" );
	} else if ( vars == DUMP_REFINEMENT_VARS ) {
		fprintf(output, "9 level processor cellnum mass_ind shock_ind contact_int density_ind pressure_ind entropy_ind\n");
	} else if ( vars == DUMP_DIFFUSION_VARS ) {
		fprintf(output, "5 level processor cellnum refinement1 refinement2\n" );
	} else if ( vars == DUMP_LOCATION_VARS ) {
		fprintf(output, "7 level print_processor cell_is_local processor_owner x y z\n");
	}

	fprintf(output, "0\n");
	fprintf(output, "%u\n", 4*total_intersected );

	for ( i = min_level; i <= max_level_dumped; i++ ) {
		select_level( i, cell_type, &num_level_cells, &level_cells );
		for ( j = 0; j < num_level_cells; j++ ) {
			print_cell_corners( level_cells[j], i );
		}
		cart_free( level_cells );
	}

	fprintf(output, "%u\n", total_intersected );

        corner_count = 1;
	for ( i = min_level; i <= max_level_dumped; i++ ) {
		select_level( i, cell_type, &num_level_cells, &level_cells );

		if ( vars == DUMP_HVARS ) {
			for ( j = 0; j < num_level_cells; j++ ) {
				print_cell_hvars( level_cells[j], i );
			}
		} else if ( vars == DUMP_GRAV_VARS ) {
			for ( j = 0; j < num_level_cells; j++ ) {
				print_cell_grav_vars( level_cells[j], i );
			}
		} else if ( vars == DUMP_REFINEMENT_VARS ) {
			for ( j = 0; j < num_level_cells; j++ ) {
				print_cell_refinement_vars( level_cells[j], i );
			}
		} else if ( vars == DUMP_DIFFUSION_VARS ) {
			for ( j = 0; j < num_level_cells; j++ ) {
				print_diffusion_vars( level_cells[j], i );
			}
		} else if ( vars == DUMP_LOCATION_VARS ) {
			for ( j = 0; j < num_level_cells; j++ ) {
				print_cell_location_vars( level_cells[j], i );
			}
		}
		
		cart_free( level_cells );
	}

	fprintf(output, "interfaces 0\n");
	fclose(output);
}

