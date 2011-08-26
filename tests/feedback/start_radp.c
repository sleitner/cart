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

#define omm0 1.0
#define oml0 0.0
#define omb0 1.0
#define hubble 1.0
#define deltadc 0.0
#define a0 0.9
#define boxh (0.001/a0*hubble) //1pkpc
#define refine_radius	(0.1*constants->kpc/units->length)
#define radp_radius	(1.0*constants->pc/units->length)
#define delta_r		(1.0*constants->pc/units->length)
#define rho0		(1.0*constants->gpercc/units->density)
#define p0		(1.0*constants->kms*constants->Msun/(units->velocity*units->mass))
#define E0		(1.0*constants->ergs/units->energy)
#define E		(1e4*constants->ergs/units->energy)

void units_set_art(double OmegaM, double h, double Lbox);

void refine_level( int cell, int level ) {
	double pos[nDim];
	float r;

	cart_assert( cell >= 0 && cell < num_cells );
	cart_assert( cell_level(cell) == level );
	
	cell_center_position(cell, pos);

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
	
void radp_initial_conditions( int icell ) {
	float r;
	double pos[nDim];
	
	cell_gas_density(icell) = rho0;
	cell_momentum(icell,0) = 0.0;
	cell_momentum(icell,1) = 0.0;
	cell_momentum(icell,2) = 0.0;
	cell_gas_gamma(icell) = constants->gamma;

	/* now add some energy  */
/* 	cell_center_position(icell, pos); */

/* 	pos[0] -= 0.5*num_grid; */
/* 	pos[1] -= 0.5*num_grid; */
/* 	pos[2] -= 0.5*num_grid; */

/* 	r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])*r0; */

/* 	if ( r <=  ) { */
/* 		cell_gas_pressure(icell) = exp( -r*r / (2.0*radp_radius*radp_radius ) ); */
/* //	if ( r <= cell_size[max_level]*r0 ) { */
/* //		cell_gas_pressure(icell) = 1.0; */
/* //              cart_debug("picked a cell!!!!!!!!!!!!!!!!!"); */
 	pos[0] = 0.5*num_grid; 
 	pos[1] = 0.5*num_grid; 
 	pos[2] = 0.5*num_grid; 
        
        if(cell_find_position(pos) == icell)
                cell_gas_pressure(icell) = p0;
	} else {
		cell_gas_pressure(icell) = p0*1e-6;
	}

	cell_gas_internal_energy(icell) = cell_gas_pressure(icell) / (cell_gas_gamma(icell)-1.0);
	cell_gas_energy(icell) = cell_gas_internal_energy(icell);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
        cell_electron_internal_energy(icell) = cell_gas_internal_energy(icell)*constants->wmu/constants->wmu_e;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
        
#ifdef ENRICH
        cell_gas_metal_density_II(icell) = 1e-30;
#ifdef ENRICH_SNIa
        cell_gas_metal_density_Ia(icell) = 1e-30;
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
}

void set_radp_initial_conditions() {
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
			radp_initial_conditions( level_cells[i] );
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

			cell_gas_internal_energy(icell) = cell_gas_internal_energy(icell)*scale_energy +  //this scales energy to init pressure 
								p0 / (cell_gas_gamma(icell)-1.0);
			cell_gas_pressure(icell) = cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0); //this adds energy back to pressure to preserve eos?
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

float intPt;
float radii[num_bins];
float vel[num_bins];
float pressure[num_bins];
float rho[num_bins];
float avgs[num_bins];

void radial_average( int cell, int level ) {
        double pos[nDim];
        float r1, r, v;
        float amt, cur_r;
        int bin;

        if ( cell_is_leaf( cell ) ) {
                cell_center_position(cell, pos);
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
        int const endian_test=-99;
	char filename[128];
	FILE *RADP;
	float reduced_rho[num_bins];
	float reduced_pressure[num_bins];
	float reduced_vel[num_bins];
	float reduced_avgs[num_bins];
	int level;
	int num_level_cells;
	int *level_cells;

	int icell, sfc, size;
	float value;
	float *slice;
	FILE *output;
	int coords[nDim];

        int max_child_level;

	const int nvars = 4;
	const int nbin1 = 128;
	int varid[] = { HVAR_PRESSURE, HVAR_GAS_DENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC };
	int nbin[] = { nbin1, nbin1, nbin1 };
	double bb[6];
        float tot_momentum=0,dPt=0 ;

        cart_debug("run_output:");
        
	bb[0] = bb[2] = bb[4] = num_grid*(0.5-0.125);
	bb[1] = bb[3] = bb[5] = num_grid*(0.5+0.125);
        
	sprintf(filename,"%s/out.%05d.bin",output_directory,step);
	ifrit.OutputMesh(filename,max_level,nbin,bb,nvars,varid);
        
        
	/* now dump the total gas momentum */
        tot_momentum=0, tot_momentum=0;
	for ( level = min_level; level <= max_level; level++ ) {
            select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
            for ( i = 0; i < num_level_cells; i++ ) {
                j=level_cells[i];
                tot_momentum += sqrt( cell_momentum(j,0)*cell_momentum(j,0) +
                                      cell_momentum(j,1)*cell_momentum(j,1) +
                                      cell_momentum(j,2)*cell_momentum(j,2) )*cell_volume[level];
                dPt += cell_gas_pressure(j)* 6*cell_size[level]*cell_size[level];
            }
            cart_free( level_cells );
	}
//        int MPI_Reduce(void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) 
        MPI_Reduce( MPI_IN_PLACE, &tot_momentum, 1, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
        MPI_Reduce( MPI_IN_PLACE, &dPt, 1, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
            
	if ( local_proc_id == MASTER_NODE ) {
		sprintf(filename, "%s/mom_%s.dat", output_directory, jobname );
		RADP = fopen(filename,"a+");
                fprintf(RADP, "%e %e %e\n", tl[min_level]-t_init , tot_momentum,dPt*(tl[min_level]-t_init) );
		fclose(RADP);
	}


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
                sprintf(filename, "%s/%s_stp%04u_tl%5.4f.txt", output_directory, jobname, step, tl[min_level]-t_init );

		RADP = fopen(filename,"w");
		for ( i = 0; i < num_bins/2; i++ ) {
			if ( reduced_avgs[i] > 0.0 ) {
				reduced_vel[i] /= reduced_avgs[i];
				reduced_pressure[i] /= reduced_avgs[i];
				reduced_rho[i] /= reduced_avgs[i];

				fprintf(RADP, "%e %e %e %e\n", radii[i]*units_length/constants->pc, reduced_rho[i], reduced_pressure[i], reduced_vel[i] );
			}
		}
		fclose(RADP);
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

		fwrite( &endian_test, sizeof(int), 1, output );
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

		coords[1] = num_grid/2;
		for ( coords[2] = 0; coords[2] < num_grid; coords[2]++ ) {
			for( coords[0] = 0; coords[0] < num_grid; coords[0]++ ) {
				sfc = sfc_index( coords );

				if ( root_cell_is_local(sfc) ) {
					icell = root_cell_location(sfc);
					value = cell_gas_pressure(icell);
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
                                        level=max_level;
                                        max_child_level=min_level;
                                        while(level>min_level && max_child_level == min_level){
                                                select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
                                                for (j=0; j<num_level_cells; j++){
                                                        if(cell_parent_root_cell(level_cells[j]) == icell){
                                                                max_child_level = level;
//   if(level>min_level){cart_debug("%d %d %d", coords[0],coords[2],max_child_level);}
                                                        }
                                                }
                                                cart_free( level_cells );
                                                level--;
                                        }
                                        value = max_child_level;
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
                cart_error("need to insert level output before you can be here");
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
					value = cell_gas_pressure(icell);
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
	int num_level_cells;
	int *level_cells;
        
	for ( i = 0; i < nDim; i++ ) {
		refinement_volume_min[i] = 0.0;
		refinement_volume_max[i] = (double)num_grid;
	}
        
        cosmology_set(OmegaM,omm0);
        cosmology_set(OmegaL,oml0);
        cosmology_set(OmegaB,omb0);
        cosmology_set(h,hubble);
        cosmology_set(DeltaDC,deltadc);
        
        box_size = boxh;
        auni_init = a0;
        auni[min_level] = a0;
        abox[min_level] = a0;
        t_init = tcode_from_auni(auni_init);
        
        units_set_art(omm0,hubble,box_size);
 	//units_set(1.0,1.0,1.0); //mass time length
 	units_reset();
	units_update( min_level );
        
////////////////////////////////////////////////////////////////////////////////////
	/* build buffer */
	build_cell_buffer();
	cart_debug("built cell buffer");
	repair_neighbors();
	cart_debug("repaired neighbors");
        
	check_map();
        
	cart_debug("setting initial conditions on root level");
	set_radp_initial_conditions(); 
        
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
        
	cart_debug("re-setting initial conditions on all levels");
	set_radp_initial_conditions(); 

        
#ifdef HYDRO_TRACERS
	cart_debug("setting hydro tracers");
	set_hydro_tracers( min_level+1 );
#endif /* HYDRO_TRACERS */
        
	cart_debug("set initial conditions");
        
        hydro_magic( min_level );
        hydro_eos( min_level );
////////////////////////////////////////////////////////////////////////////////////
        
        /* set time variables */
        tl[min_level] = t_init;
        for(level=min_level+1; level<=max_level; level++)
        {
                tl[level] = tl[min_level];
                auni[level] = auni[min_level];
                abox[level] = abox[min_level];
        }
        cosmology_set_fixed();
        
        cart_debug("tl[min_level] = %f", tl[min_level] );
        cart_debug("au[min_level] = %f", auni[min_level] );
        cart_debug("ab[min_level] = %f", abox[min_level] );
        cart_debug("DC mode = %f", cosmology->DeltaDC );
        

#ifdef PARTICLES
        for(i=0; i<num_particles; i++) if(particle_level[i] != FREE_PARTICLE_LEVEL)
        {
//  We set the step to 0 so that the first leapfrog step is correct
                particle_t[i] = tl[min_level];
                particle_dt[i] = 0.0;
        }
#endif
        
        
        
        for(level = min_level; level <= max_level; level++){
                cart_debug("updating level %u", level );
                update_buffer_level(level, all_hydro_vars, num_hydro_vars);
        }
        
        cart_debug("done updating initial conditions");
        
#ifdef PARTICLES
        build_mesh(); //defines refinement volume
        
        for ( i = 0; i < nDim; i++ ) {
                refinement_volume_min[i] = 0;
                refinement_volume_max[i] = (double)num_grid;
        }
        for ( i = 0; i < nDim; i++ ) {
                star_formation_volume_min[i] = refinement_volume_min[i];
                star_formation_volume_max[i] = refinement_volume_max[i];
        }
#endif /* PARTICLES */
        
        if ( !buffer_enabled ) {
                cart_debug("building cell buffer");
                build_cell_buffer();
                repair_neighbors();
        }
        
        cart_debug("done with initialization");
        
	check_map();
}





