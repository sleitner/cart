#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include <sys/types.h>
#include <unistd.h>

#include "tree.h"
#include "cosmology.h"
#include "particle.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "timestep.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "gravity.h"
#include "density.h"
#include "io.h"
#include "auxiliary.h"

#ifdef VIEWDUMP
#include "extras/viewdump.h"
#endif

#include "extra/ifrit.h"

#define refine_radius	(4.5)
#define r0		(1.0/num_grid)
#define sedov_radius	(cell_size[max_level]/num_grid)
#define delta_r		(1.0)
#define rho0		(1.0)
#define p0		(1.0)
#define E0		(1.0)
#define E		(1e7)

void refine_level( int cell, int level ) {
	float pos[nDim];
	float r;

	cart_assert( cell >= 0 && cell < num_cells );
	cart_assert( cell_level(cell) == level );
	
	cell_position(cell, pos);

	pos[0] -= 0.5*num_grid;
	pos[1] -= 0.5*num_grid;
	pos[2] -= 0.5*num_grid;

	r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);

	if ( r < refine_radius*cell_size[level] ) {
		refinement_indicator(cell,0) = 1.0;	
	} else {
		refinement_indicator(cell,0) = 0.0;
	}
}
	
void sedov_initial_conditions( int icell ) {
	float r;
	float pos[nDim];
	
	cell_gas_density(icell) = 1.0/rho0;
	cell_momentum(icell,0) = 0.0;
	cell_momentum(icell,1) = 0.0;
	cell_momentum(icell,2) = 0.0;
	cell_gas_gamma(icell) = constants->gamma;

	/* now add some energy  */
	cell_position(icell, pos);

	pos[0] -= 0.5*num_grid;
	pos[1] -= 0.5*num_grid;
	pos[2] -= 0.5*num_grid;

	r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])*r0;

	if ( r <= 4.0*r0*cell_size[min_level] ) {
		cell_gas_pressure(icell) = exp( -r*r / (2.0*sedov_radius*sedov_radius ) );
	} else {
		cell_gas_pressure(icell) = 0.0;
	}

	cell_gas_internal_energy(icell) = cell_gas_pressure(icell) / (cell_gas_gamma(icell)-1.0);
	cell_gas_energy(icell) = cell_gas_internal_energy(icell);
}

void set_sedov_initial_conditions() {
	int i;
	int icell;
	int level;
	int num_level_cells;
	int *level_cells;
	double scale;
	double scale_energy;

	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			sedov_initial_conditions( level_cells[i] );
		}
		cart_free( level_cells );
	}

	scale_energy = 0.0;
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_leaf(icell) ) {
				scale_energy += cell_gas_internal_energy(icell)*cell_volume[level];
			}
		}
		cart_free( level_cells );
	}
	scale = scale_energy;
	
	MPI_Allreduce( &scale, &scale_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	scale_energy = (E/E0) / scale_energy;
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			cell_gas_internal_energy(icell) = cell_gas_internal_energy(icell)*scale_energy + 
								p0 / ( (cell_gas_gamma(icell)-1.0) );
			cell_gas_pressure(icell) = cell_gas_internal_energy(icell) * ( (cell_gas_gamma(icell)-1.0) );
			cell_gas_energy(icell) = cell_gas_internal_energy(icell);
		}
	
		cart_free( level_cells );
	}

	for ( level = max_level - 1; level >= min_level; level-- ) {
		hydro_split_update(level);
	}
}

/* radial binning for analysis */
#define num_bins        (num_grid*(1<<(max_level-min_level)))
#define bin_width       ((float)num_grid/((float)num_bins))

float radii[num_bins];
float vel[num_bins];
float pressure[num_bins];
float rho[num_bins];
float avgs[num_bins];

void radial_average( int cell, int level ) {
        float pos[nDim];
        float r1, r, v;
        float amt, cur_r;
        int bin;
                                                                                                                                                            
        if ( cell_is_leaf( cell ) ) {
                cell_position(cell, pos);
                r = sqrt( (pos[0]-((float)num_grid/2.0))*(pos[0]-((float)num_grid/2.0))+
                        (pos[1]-((float)num_grid/2.0))*(pos[1]-((float)num_grid/2.0))+
                        (pos[2]-((float)num_grid/2.0))*(pos[2]-((float)num_grid/2.0)));
                                                                                                                                                            
                /* this needs to be converted to v_r */
                v = sqrt( cell_momentum(cell,0)*cell_momentum(cell,0) +
                        cell_momentum(cell,1)*cell_momentum(cell,1) +
                        cell_momentum(cell,2)*cell_momentum(cell,2) ) / cell_gas_density(cell);
                                                                                                                                                            
                /* 1) determine which bin left hand edge of cell is in
                 * 2) go outwards in bins until right hand edge of bin > right hand edge of cell
                 * 3) apply proper amounts to each cell */
                                                                                                                                                            
                bin = (int)((r-0.5*cell_size[level])/bin_width);
                cart_assert( bin >= 0 && bin < num_bins );

                cur_r = r-0.5*cell_size[level];
                amt = 1.0;

                while ( (float)(bin+1)*bin_width < (r+0.5*cell_size[level]) ) {
                        r1 = ( (float)(bin+1)*bin_width - cur_r )/cell_size[level];

                        vel[bin] += r1*v;
                        rho[bin] += r1*cell_gas_density(cell);
                        pressure[bin] += r1*cell_gas_pressure(cell);
                        avgs[bin] += r1;
                                                                                                                                                            
                        amt -= r1;
                        cart_assert( amt >= 0.0 );
                                                                                                                                                            
                        cur_r = (float)(bin+1)*bin_width;
                        bin++;
                }

                if ( amt > 0.0 && bin < num_bins - 1 ) {
                        /* apply amt to last bin */
                        r1 = amt;
                        vel[bin] += r1*v;
                        rho[bin] += r1*cell_gas_density(cell);
                        pressure[bin] += r1*cell_gas_pressure(cell);
                        avgs[bin] += r1;
                }
        }
}

void run_output() {
	int i, j;
	char filename[128];
	FILE *SEDOV;
	float reduced_rho[num_bins];
	float reduced_pressure[num_bins];
	float reduced_vel[num_bins];
	float reduced_avgs[num_bins];
	int level;
	int num_level_cells;
	int *level_cells;
	float pos[nDim];
	int cell;
	int min_index, max_index;

	int icell, sfc, size;
	float value;
	float *slice;
	FILE *output;
	int coords[nDim];

	const int nvars = 4;
	const int nbin1 = 128;
	int varid[] = { HVAR_PRESSURE, HVAR_GAS_DENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC };
	int nbin[] = { nbin1, nbin1, nbin1 };
	double bb[6];


	bb[0] = bb[2] = bb[4] = num_grid*(0.5-0.125);
	bb[1] = bb[3] = bb[5] = num_grid*(0.5+0.125);

	sprintf(filename,"%s/out.%05d.bin",output_directory,step);
	ifrit.OutputMesh(filename,max_level,nbin,bb,nvars,varid);

	/* now dump the radial profiles */
	for ( i = 0; i < num_bins; i++ ) {
		radii[i] = ((float)i + 0.5) * bin_width;
		vel[i] = 0.0;
		pressure[i] = 0.0;
		rho[i] = 0.0;
		avgs[i] = 0.0;
	}

	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			radial_average( level_cells[i], level );
		}
		cart_free( level_cells );
	}

	MPI_Reduce( rho, reduced_rho, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( pressure, reduced_pressure, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( vel, reduced_vel, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( avgs, reduced_avgs, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );

	if ( local_proc_id == MASTER_NODE ) {
		sprintf(filename, "%s/%s_[%5.4f].dat", output_directory, jobname, tl[min_level]-t_init );

		SEDOV = fopen(filename,"w");
		for ( i = 0; i < num_bins/2; i++ ) {
			if ( reduced_avgs[i] > 0.0 ) {
				reduced_vel[i] /= reduced_avgs[i];
				reduced_pressure[i] /= reduced_avgs[i];
				reduced_rho[i] /= reduced_avgs[i];

				fprintf(SEDOV, "%e %e %e %e\n", radii[i]*r0, reduced_rho[i], reduced_pressure[i], reduced_vel[i] );
			}
		}
		fclose(SEDOV);
	}

#ifdef VIEWDUMP
	sprintf(filename, "%s/%s_%04u.v", output_directory, jobname, step );
        viewdump( filename, max_level, ((float)(num_grid)/2.0)+0.5*cell_size[max_level], 2, DUMP_HVARS, CELL_TYPE_LOCAL );
#endif

#ifdef HYDRO_TRACERS
	sprintf( filename, "%s/tracers_%04u.dat", output_directory, step );
	write_hydro_tracers( filename );
#endif /* HYDRO_TRACERS */

	/* output a 2-d slice through the center of the box */
	if ( local_proc_id == MASTER_NODE ) {
		sprintf( filename, "%s/%s_slice_%04u.dat", output_directory, jobname, step );
		output = fopen( filename, "w" );

		size = num_grid;
		fwrite( &size, sizeof(int), 1, output );
		size = num_grid;
		fwrite( &size, sizeof(int), 1, output );

		slice = cart_alloc(float, num_grid*num_grid );

		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_gas_density(icell);
				} else {
					MPI_Recv( &value, 1, MPI_FLOAT, processor_owner(sfc), sfc,
							MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}

				slice[ coords[2]*num_grid + coords[0] ] = value;
			}
		}

		fwrite( slice, sizeof(float), num_grid*num_grid, output );

		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_gas_internal_energy(icell)/cell_gas_density(icell);
				} else {
					MPI_Recv( &value, 1, MPI_FLOAT, processor_owner(sfc), sfc,
							MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}

				slice[ coords[2]*num_grid + coords[0] ] = value;
			}
		}

		fwrite( slice, sizeof(float), num_grid*num_grid, output );

		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_momentum(icell,0)/cell_gas_density(icell);
				} else {
					MPI_Recv( &value, 1, MPI_FLOAT, processor_owner(sfc), sfc,
							MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}

				slice[ coords[2]*num_grid + coords[0] ] = value;
			}
		}

		fwrite( slice, sizeof(float), num_grid*num_grid, output );

		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_momentum(icell,2)/cell_gas_density(icell);
				} else {
					MPI_Recv( &value, 1, MPI_FLOAT, processor_owner(sfc), sfc,
							MPI_COMM_WORLD, MPI_STATUS_IGNORE );
				}

				slice[ coords[2]*num_grid + coords[0] ] = value;
			}
		}

		fwrite( slice, sizeof(float), num_grid*num_grid, output );

		cart_free( slice );
		fclose(output);
	} else {
		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_gas_density(icell);
					MPI_Send( &value, 1, MPI_FLOAT, MASTER_NODE, sfc, MPI_COMM_WORLD );
				}
			}
		}

		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_gas_internal_energy(icell)/cell_gas_density(icell);
					MPI_Send( &value, 1, MPI_FLOAT, MASTER_NODE, sfc, MPI_COMM_WORLD );
				}
			}
		}

		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_momentum(icell,0)/cell_gas_density(icell);
					MPI_Send( &value, 1, MPI_FLOAT, MASTER_NODE, sfc, MPI_COMM_WORLD );
				}
			}
		}

		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_momentum(icell,1)/cell_gas_density(icell);
					MPI_Send( &value, 1, MPI_FLOAT, MASTER_NODE, sfc, MPI_COMM_WORLD );
				}
			}
		}
	}
}

void init_run() {
	int i;
	int level;
	int ioct;
	int num_level_cells;
	int *level_cells;
	char filename[128];
	int min_index, max_index;

	for ( i = 0; i < nDim; i++ ) {
		refinement_volume_min[i] = 0.0;
		refinement_volume_max[i] = (double)num_grid;
	}

	cart_debug("in init");

 	units_set(1.0,1.0,1.0);
 	units_reset();

	/* build buffer */
	build_cell_buffer();
	cart_debug("built cell buffer");
	repair_neighbors();

	check_map();

	cart_debug("repaired neighbors");

	/* do initial refinements */
	for ( level = min_level; level < max_level; level++ ) {
		cart_debug("refining level %u", level );

		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		cart_debug("num_level_cells = %u", num_level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			refine_level( level_cells[i], level );
		}
		cart_free( level_cells );
		cart_debug("about to refine level %u", level );
		refine(level);
	}

	cart_debug("setting initial conditions");
	set_sedov_initial_conditions();

#ifdef HYDRO_TRACERS
	cart_debug("setting hydro tracers");
	set_hydro_tracers( min_level+1 );
#endif /* HYDRO_TRACERS */

	cart_debug("set initial conditions");

	for ( level = min_level; level <= max_level; level++ ) {
		cart_debug("updating level %u", level );
		update_buffer_level( level, all_hydro_vars, num_hydro_vars );
	}

	cart_debug("done updating initial conditions");

	/* set time variables */
	tl[min_level] = t_init;

	cart_debug("done with initialization");

	check_map();
}
