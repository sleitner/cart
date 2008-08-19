#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defs.h"
#include "io.h"
#include "timestep.h"
#include "iterators.h"
#include "timing.h"
#include "units.h"
#include "tree.h"
#include "hydro.h"
#include "particle.h"
#include "auxiliary.h"
#include "cell_buffer.h"
#include "starformation.h"
#include "logging.h"

FILE *steptimes;
FILE *timing;
FILE *energy;
FILE *workload;
FILE *dependency;

#ifdef STARFORM
FILE *star_log;
#endif /* STARFORM */

void init_logging( int restart ) {
	int i;
	char mode[2];
	char filename[256];

	if ( logfile_directory[0] == '\0' ) {
		strcpy( logfile_directory, output_directory );
	}

	if ( restart ) {
		mode[0] = 'a';
	} else {
		mode[0] = 'w';
	}
	mode[1] = '\0';
	
	if ( local_proc_id == MASTER_NODE ) {
		/* open log files */
		sprintf(filename,"%s/times.log", logfile_directory );
		steptimes = fopen(filename,mode);
                                                                                                                                                            
		if ( steptimes == NULL ) {
			cart_error("Unable to open %s for writing!", filename );
		}

		if ( !restart || restart == 2 ) {
			fprintf(steptimes,"# step tl dtl aexp daexp\n");
		}

                sprintf(filename, "%s/energy.log", logfile_directory );
                energy = fopen(filename,mode);

                if ( energy == NULL ) {
                        cart_error("Unable to open %s for writing!", filename );
                }

		if ( !restart || restart == 2 ) {
			fprintf(energy, "# step tl aexp gas_thermal gas_kinetic gas_potential total_gas_energy baryon_mass particle_kinetic particle_potential total_particle_energy error\n");
		}

#ifdef STARFORM
		sprintf( filename, "%s/sf.log", logfile_directory );
		star_log = fopen(filename, mode);

		if ( star_log == NULL ) {
			cart_error("Unable to open %s for writing!", filename );
		}

		if ( !restart ) {
			fprintf(star_log, "# step t dt aexp t [Gyrs] dt [yrs] N* M* dM* Mi* dMi* [Msun] SFR [Msun/yr/Mpc^3]\n");  
		}
#endif /* STARFORM */
	}

	sprintf(filename, "%s/timing.%03u.log", logfile_directory, local_proc_id );
	timing = fopen(filename,mode);

	if ( timing == NULL ) {
		cart_error("Unable to open %s for writing!", filename );
	}

	if ( !restart || restart == 2 ) {
		/* write header for timing file (total time treated specially) */
		fprintf( timing, "# %u levels %u timers\n", num_refinement_levels+1, NUM_TIMERS-1 );
		fprintf( timing, "# step t aexp total_time" );
		for ( i = 1; i < NUM_TIMERS; i++ ) {
			fprintf( timing, " %s", timer_name[i] );
		}
		fprintf( timing, "\n" );
		fflush(timing);
	}

	sprintf(filename, "%s/workload.%03u.dat", logfile_directory, local_proc_id );
	workload = fopen(filename, mode);

	if ( workload == NULL ) {
		cart_error("Unable to open %s for writing!", filename );
	}

	if ( !restart || restart == 2 ) {
		/* write header for workload file */
		fprintf( workload, "# %u levels\n", num_refinement_levels+1 );
		fprintf( workload, "# step t aexp num_local_particles max_level" );
		for ( i = min_level; i <= max_level; i++ ) {
			fprintf( workload, " num_local_cells[%u] num_buffer_cells[%u]", i, i );
		}
		fprintf( workload, "\n" );
		fflush(workload);
	}

	sprintf(filename, "%s/dependency.%03u.dat", logfile_directory, local_proc_id );
	dependency = fopen( filename, "a" );

	if ( dependency == NULL ) {
		cart_error( "Unable to open %s for writing!", filename );
	}


#ifdef DEBUG
	log_in_debug(-1,0,__FILE__,__LINE__);
#endif
}

void finalize_logging() {
	if ( local_proc_id == MASTER_NODE ) {
		/* close log files */
		fclose(steptimes);
		fclose(energy);
	}

	fclose(timing);
	fclose(workload);
}

void log_diagnostics() {
	int i,j;
	int level;
	int icell;
	int num_level_cells;
	int *level_cells;
	double kinetic_energy;
	double gas_kinetic, gas_thermal, gas_potential, gas_mass;
	double total_gas_kinetic, total_gas_thermal, total_gas_potential, total_gas_mass;
	double particle_kinetic, particle_potential;
	double total_particle_kinetic, total_particle_potential;
	double error;
	double da;

#ifdef STARFORM
	double stellar_mass, stellar_initial_mass;
	double old_stellar_mass, old_stellar_initial_mass;
	double d_stellar_mass, d_stellar_initial_mass;
	double resolved_volume[max_level-min_level+1];
	double local_resolved_volume[max_level-min_level+1];
	double total_resolved_volume;
	double dtyears, current_age;
	double sfr;
#endif /* STARFORM */

#ifdef PARTICLES
#ifdef COSMOLOGY
	ap1 = b2a( tl[min_level] - 0.5*dtl[min_level] );
	da = aexp[min_level] - aexp_old[min_level];
#else
	ap1 = 1.0;
	ap0 = 0.0;
	da = 0.0;
#endif
#endif
	
	/* log profiling information */
	fprintf( timing, "%u %e %e %e", step, tl[min_level], aexp[min_level], current_time( TOTAL_TIME, min_level ) );
	for ( level = min_level; level <= max_level; level++ ) {
		for ( i = 1; i < NUM_TIMERS; i++ ) {
			fprintf( timing, " %e", total_time(i, level) );
		}
	}
	fprintf( timing, "\n" );
	fflush(timing);

	/* log workload information */
#ifdef PARTICLES
	fprintf( workload, "%u %e %e %u %u", step, tl[min_level], aexp[min_level], num_local_particles, max_level_now() );
#else
	fprintf( workload, "%u %e %e 0 %u", step, tl[min_level], aexp[min_level], max_level_now() );
#endif /* PARTICLES */

	for ( i = min_level; i <= max_level; i++ ) {
		fprintf(workload, " %u %u", num_cells_per_level[i], num_buffer_cells[i] );
	}
	fprintf(workload, "\n");
	fflush(workload);

	/* log dependency information */
	fprintf( dependency, "%u %e %e", step, tl[min_level], aexp[min_level] );
	for ( level = min_level; level <= max_level; level++ ) {
		for ( i = 0; i < num_procs; i++ ) {
			if ( level == min_level ) {
				fprintf( dependency, " %u %u", num_remote_buffers[level][i],
					num_local_buffers[level][i] );
			} else {
				fprintf( dependency, " %u %u", num_children*num_remote_buffers[level][i],
					num_children*num_local_buffers[level][i] );
			}
		}
	}
	fprintf(dependency, "\n");
	fflush(dependency);

	/* compute energies */
	gas_kinetic = gas_thermal = gas_potential = gas_mass = 0.0;
	total_gas_kinetic = total_gas_thermal = total_gas_potential = total_gas_mass = 0.0;
#ifdef GRAVITY
#ifdef HYDRO
	for ( level = min_level; level <= max_level; level++ ) {
		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_leaf( icell ) ) {
				gas_thermal += cell_volume[level]*cell_gas_pressure(icell)/(gamma-1.0);

				kinetic_energy = 0.0;
				for ( j = 0; j < nDim; j++ ) {
					kinetic_energy += cell_momentum(icell,j)*cell_momentum(icell,j);
				}
				kinetic_energy *= 0.5*cell_volume[level]/cell_gas_density(icell);

				gas_kinetic += kinetic_energy;
				gas_potential += cell_gas_density(icell)*cell_volume[level]*cell_potential(icell);
				gas_mass += cell_gas_density(icell)*cell_volume[level];
			}
		}
		cart_free( level_cells );
	}

	/* add stellar mass to gas mass */
#ifdef STARFORM
	stellar_mass = 0.0;
	stellar_initial_mass = 0.0;
	for ( i = 0; i < num_star_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL && particle_is_star(i) ) {
			gas_mass += particle_mass[i];
			stellar_mass += particle_mass[i];
			stellar_initial_mass += star_initial_mass[i];
		}
	}
#endif /* STARFORM */

	MPI_Reduce( &gas_thermal, &total_gas_thermal, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( &gas_kinetic, &total_gas_kinetic, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( &gas_potential, &total_gas_potential, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( &gas_mass, &total_gas_mass, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
#endif /* HYDRO */
#endif /* GRAVITY */

#ifdef STARFORM
	dtyears = dtl[min_level] * t0 * aexp[min_level]*aexp[min_level];
	current_age = age(aexp[min_level]);
	
	old_stellar_mass = total_stellar_mass;
	old_stellar_initial_mass = total_stellar_initial_mass;

	cart_debug("stellar_mass = %e", stellar_mass );
	cart_debug("stellar_initial_mass = %e", stellar_initial_mass );

	MPI_Reduce( &stellar_mass, &total_stellar_mass, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( &stellar_initial_mass, &total_stellar_initial_mass, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );

	if ( local_proc_id == MASTER_NODE ) {
		cart_debug("total_stellar_mass = %e", total_stellar_mass );
	}

	d_stellar_mass = total_stellar_mass - old_stellar_mass;
	d_stellar_initial_mass = total_stellar_initial_mass - old_stellar_initial_mass;

	/* compute resolved volume */
	local_resolved_volume[min_level] = 0.0;
	for ( level = min_level+1; level <= max_level; level++ ) {
		local_resolved_volume[level] = 0.0;

		select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_leaf( icell ) ) {
				local_resolved_volume[level] += cell_volume[level];
			}
		}

		cart_free(level_cells);
	}

	MPI_Reduce( local_resolved_volume, resolved_volume, max_level-min_level+1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );

	if ( local_proc_id == MASTER_NODE ) {
		/* sum resolved volume over all levels except min_level (works for MMR simulations) */
		total_resolved_volume = 0.0;
		for ( level = min_level; level <= max_level; level++ ) {
			cart_debug("resolved_volume[%u] = %e", level, resolved_volume[level] );
			total_resolved_volume += resolved_volume[level];
		}
		total_resolved_volume *= (Lbox*Lbox*Lbox) / num_root_cells;

		cart_debug("total_resolved_volume = %e", total_resolved_volume);

		if ( total_resolved_volume > 0.0 ) {
			sfr = d_stellar_initial_mass * aM0 * hubble * hubble * hubble / dtyears / total_resolved_volume;
		} else {
			sfr = 0.0;
		}

		fprintf( star_log, "%u %e %e %e %e %e %u %e %e %e %e %e\n", step, tl[min_level],
			dtl[min_level], aexp[min_level], current_age, dtyears, 
			particle_species_num[num_particle_species-1],
			total_stellar_mass*aM0, d_stellar_mass*aM0,
			total_stellar_initial_mass*aM0, d_stellar_initial_mass*aM0,
			sfr );

		fflush(star_log);
	}
#endif /* STARFORM */


	particle_kinetic = particle_potential = 0.0;
	total_particle_kinetic = total_particle_potential = 0.0;
#ifdef PARTICLES
	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL ) {
			kinetic_energy = 0.0;
			for ( j = 0; j < nDim; j++ ) {
				kinetic_energy += particle_v[i][j]*particle_v[i][j];
			}
			kinetic_energy *= particle_mass[i];

			particle_kinetic += kinetic_energy;
#ifdef GRAVITY
			particle_potential += particle_pot[i];
#endif /* GRAVITY */
		}
	}

	MPI_Reduce( &particle_kinetic, &total_particle_kinetic, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( &particle_potential, &total_particle_potential, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );

	total_particle_kinetic *= 0.5/(ap1*ap1);
	total_particle_potential *= 0.5;
#endif

	if ( local_proc_id == MASTER_NODE ) {
#ifdef PARTICLES
		ekin = total_particle_kinetic;

#ifdef COSMOLOGY
		if ( step == 1 ) {
			au0 = aexp_old[min_level]*total_particle_potential;
			aeu0 = au0 + aexp_old[min_level]*total_particle_kinetic;
			tintg = 0.0;
			error = 0.0;
		} else {
			tintg += 2.0*(aexp_old[min_level] - ap0)*ekin;
			kinetic_energy = (ekin*(ap1-ap0-0.5*da) + ekin1*0.5*da)/(ap1-ap0);
			error = ( aexp_old[min_level]*(kinetic_energy+total_particle_potential) - aeu0 + tintg ) /
					( aexp_old[min_level]*total_particle_potential - au0 );
		}
#else
		if ( step == 1 ) {
			au0 = total_particle_potential;
			aeu0 = au0 + total_particle_kinetic;
			tintg = 0.0;
			error = 0.0;
		} else {
 		        kinetic_energy = 0.5*(ekin+ekin1);
			if(aeu0 != 0.0) error = (kinetic_energy+total_particle_potential)/aeu0 - 1.0; else error = 0.0;
		}
#endif  /* COSMOLOGY */

#endif /* PARTICLES */

		fprintf(energy, "%u %e %e %e %e %e %e %e %e %e %e %e\n",
			step, tl[min_level], aexp[min_level], 
			total_gas_thermal, total_gas_kinetic, total_gas_potential, 
			total_gas_thermal+total_gas_kinetic+total_gas_potential,
			total_gas_mass,
			total_particle_kinetic, total_particle_potential,
			total_particle_kinetic+total_particle_potential,
			error );
		fflush(energy);

		fprintf(steptimes, "%u %e %e %e %e\n", step, tl[min_level], dtl[min_level], aexp[min_level],
			aexp[min_level]-aexp_old[min_level] );
		fflush(steptimes);
        }

#ifdef PARTICLES
	ekin1 = total_particle_kinetic;
	ap0 = ap1;
#endif

}


#ifdef DEBUG

extern const char *timer_name[];
double offset = 0.0;

void log_in_debug(int timerid, int start, const char *file, int line)
{
#if (DEBUG-0 > 1)
  char filename[256];
  FILE *f = 0;
  
  sprintf(filename,"%s/debug.%03u.log",logfile_directory,local_proc_id);
  if(timerid < 0)
    {
      offset = MPI_Wtime();
      f = fopen(filename,"w");
      if(f != 0) fclose(f);
      return;
    }
  
  f = fopen(filename,"a");
  if(f != 0)
    {
      switch(start)
	{
	case 0:
	  {
	    fprintf(f,"%s @ %d: %s done at %f sec.\n",file,line,timer_name[timerid],MPI_Wtime()-offset);
	    break;
	  }
	case 1:
	  {
	    fprintf(f,"%s @ %d: %s started at %f sec.\n",file,line,timer_name[timerid],MPI_Wtime()-offset);
	    break;
	  }
	default:
	  {
	    fprintf(f,"%s @ %d: marker #%d set at %f sec.\n",file,line,timerid,MPI_Wtime()-offset);
	  }
	}
      fclose(f);
    }
#endif
}
#endif
