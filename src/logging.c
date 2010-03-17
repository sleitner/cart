#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "hydro.h"
#include "io.h"
#include "iterators.h"
#include "logging.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "control_parameter.h"
#include "timestep.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#ifdef MPE_LOG
#include <mpe.h>
#endif
#ifdef LOG_STAR_CREATION
#include "starformation_recipes.h"
#endif

FILE *steptimes;
FILE *timing;
FILE *energy;
FILE *workload;
FILE *dependency;

#ifdef STARFORM
FILE *star_log;
#endif /* STARFORM */


#ifdef DEBUG_MEMORY_USE
unsigned long dmuReportAllocatedMemory();
#endif /* DEBUG_MEMORY_USE */

#ifdef STARFORM
#ifdef LOG_STAR_CREATION

void output_star_creation( int icell, double mass, FILE *f ){
  int level_KS, icell_KS;
  double side_KS;
  int level;
  int id;
  int i;

  id = last_star_id + local_proc_id + 1;

  level_KS = ( log( 1.0*auni[level]*(constants->kpc/units->length)/cell_size[min_level] ) / log(2.0) ); //want KS to be ~1pKpc
  level = cell_level(icell);
  icell_KS = icell;
  if( level > level_KS ){
    for( i = level; i < level_KS; i++ ){
      icell_KS = cell_parent_cell(icell_KS);
    }
  }
  side_KS = cell_size[ cell_level(icell_KS) ]*(units->length/constants->kpc)*auni[level];
  
  //starid, auni,tcosmo, level, over_density,cell_gas density, consumption timescale, mass, icell_KS,KSlength(Kpc), KS gas density
  fprintf(f,"%d %d  %f %e  %d %e %e  %e  %e  %d %f %e\n",
	  local_proc_id, id, 
	  auni_from_tcode(tl[level]) , tphys_from_tcode(tl[level]) ,
	  level , cell_density(icell) , cell_gas_density(icell)*units->density/(constants->Msun/pow(constants->pc,3)) ,
	  ( cell_gas_density(icell)*cell_volume[level]/sf_rate(icell) )*(units->time/constants->yr),
	  mass*units->mass/constants->Msun , 
	  icell_KS , side_KS , cell_gas_density(icell_KS)*units->density/(constants->Msun/pow(constants->pc,3)) 
	  );
}


void log_star_creation( int icell, double mass, int fileop_temp ){
  static char filename[256];
  static FILE *temp_file;

  switch(fileop_temp)
    {
    case FILE_OPEN:
      {
	sprintf(filename,"%s/star_temp.%03u.log",output_directory,local_proc_id);
	temp_file = fopen(filename,"w");
	cart_debug("opening %s",filename);
	if ( temp_file == NULL ){cart_error("Unable to open %s!", filename );}
	break;
      }
    case FILE_RECORD:
      {
	cart_assert( icell >= 0 && icell < num_cells );
	cart_assert( mass > 0.0 );
	if ( temp_file == NULL ) {cart_error("Unable to use %s for writing!", filename );}
	output_star_creation( icell, mass, temp_file );
	break;
      }
    case FILE_CLOSE:
      {
	cart_debug("closing %s",filename);
	if ( temp_file == NULL ){cart_error("Unable to close %s!", filename );}
	fclose(temp_file);
	break;
      }
    default:
      cart_error("bad option");
    } 
}

void check_restart_star_creation(){
  /* Don't want ot restart from jobname.d and have a corrupt/empty/inaccurate restart_sc file.
   * Match the number of stars since the last output ae in the code to the number of stars listed in screst.
   */
  double last_ae;
  float fdummy;
  int local_count_stars,count_stars,count_recent_stars, nsclog, iout,i;
  int max_len=256;
  char str_buf[max_len+1];
  char filename_screst[256];
  FILE *fp;
  cart_debug("checking the restart_star_creation file for consistency with stars in memory");
  sprintf(filename_screst,"%s/restart_star_creation.log",output_directory);

  local_count_stars = 0;  count_stars = 0;
  for ( i = 0; i < num_star_particles; i++ ) {
    if ( particle_level[i] != FREE_PARTICLE_LEVEL && particle_is_star(i) ) {
      local_count_stars++ ;
    }
  }
  MPI_Allreduce( &local_count_stars, &count_stars, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  cart_debug("At check_restart_sc %d stars exist",count_stars); 

  /*read in the last labeled restart*/ 
  if ( local_proc_id == MASTER_NODE ) {
    iout=1;while ( iout < num_outputs && auni[min_level] >= outputs[iout] ) {iout++;}
    iout--;

    if( count_stars == 0 ){ // no stars: wipe restart
      last_ae = auni[min_level];
      wipe_restart_star_creation( last_ae );
    } 
    else{ // starting from labeled restart: wipe restart (last stars were recorded) 
      if(auni_from_abox(abox_old[min_level]) < outputs[iout] && auni[min_level] > outputs[iout]){  
	last_ae = auni[min_level];
	wipe_restart_star_creation( last_ae );
      } 
      else {
	if((fp=fopen(filename_screst,"r")) == NULL) {
	  cart_error("restart_sc must contain last_ae if (not just initialized after a labeled output or no stars exist)");
	} else {
	  if(fgets(str_buf, max_len + 1, fp) != NULL){/* read last_ae if it exists */
	    i = sscanf(str_buf, "%f", &fdummy);
	    last_ae=fdummy;
	    cart_debug("found last_ae=%f restart_sc file",last_ae);
	  }else {
	    cart_error("restart_sc must contain last_ae if it has not just been initialized after a labeled output or no stars exist");
	  }
	  fclose(fp);
	}
      }
    }
    /* check last_ae is inside relevant interval */
    if( (last_ae < outputs[iout] || last_ae > auni[min_level]) && iout > 0 ){
      cart_error("bad last_ae: %f < %f || %f > %f",last_ae,outputs[iout],last_ae,auni[min_level]);
    }
    cart_debug("star creation outputs check: iout=%d a=%f last_ae=%f outi=%f outi+1=%f",i, auni[min_level],last_ae,outputs[iout],outputs[iout+1]);
  }
  MPI_Bcast(&last_ae, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD );

  if(last_ae != auni[min_level] && count_stars > 0){
    cart_debug("counting recent star creations after last_ae=%f",last_ae);
    /*count stars since the last restart*/ 
    local_count_stars = 0; count_recent_stars = 0;
    for ( i = 0; i < num_star_particles; i++ ) {
      if ( particle_level[i] != FREE_PARTICLE_LEVEL && particle_is_star(i) ) {
	if ( auni_from_tcode(star_tbirth[i]) > last_ae ) {
	  local_count_stars++ ;
	}
      }
    }
    MPI_Reduce( &local_count_stars, &count_recent_stars, 1, MPI_INT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
    
    /* make sure your restart_sc file matches the number of stars since the last backup to a labeled file */
    if ( local_proc_id == MASTER_NODE ) {
      if( count_recent_stars > 0 ){
	nsclog = count_lines( filename_screst );
	nsclog = nsclog - 1;// one line is a header
	if( nsclog != count_recent_stars ){
	  cart_error( "restart_star_creation.log has a bad star count: (log:%d) != (dst:%d)",nsclog,count_recent_stars);
	}else{
	  cart_debug("total number of recent star creations since last output=%d matches restart line count=%d",count_recent_stars, nsclog);
	}
      }else{
	cart_debug("no stars newer at a=%f than last_ae=%f",auni[min_level],last_ae);
      }
    }
  }
  cart_debug("done checking restart_star_creation.log");


}

void combine_star_creation_log(){
  char filename_screst[256];
  char filename_local[256];
  int proc;
  FILE *fp;

  sprintf(filename_screst,"%s/restart_star_creation.log",output_directory);
  for ( proc = 0; proc < num_procs; proc++ ) {
    sprintf(filename_local,"%s/star_temp.%03d.log",output_directory,proc);
    fp = fopen(filename_local,"r");
    if( fp != NULL ) {fclose(fp); 
      cart_debug("combine_star_creation_log %s to %s",filename_local,filename_screst);
      append_file( filename_local, filename_screst);
    }
  }
}

void wipe_restart_star_creation( double last_ae ){
  char filename_screst[256];
  int proc;
  FILE *fp;
  sprintf(filename_screst,"%s/restart_star_creation.log",output_directory);
  if((fp=fopen(filename_screst,"w")) == NULL) {cart_error("Cannot open file %s",filename_screst);}
  fprintf(fp,"%f\n",last_ae);
  fclose(fp);
}

void finalize_star_creation_log( char *filename_sclog ){
  char filename_screst[256];
  FILE *fp;

  cart_debug("finalize_star_creation");
  sprintf(filename_screst,"%s/restart_star_creation.log",output_directory);

  fp = fopen(filename_screst,"r"); //make sure restart_sc exists
  if( fp != NULL ) {
    fclose(fp);
    copy_file( filename_screst, filename_sclog);
    wipe_restart_star_creation( auni[min_level] );
  } 
}


int count_lines(char *filename){
  int i=0;
  char ch='\0';
  FILE *fp;

  fp=fopen(filename,"r");
  if ( fp == NULL ) {cart_error("Cannot open count file: %s",filename);}
  while(ch!=EOF) {
    ch=fgetc(fp);
    if(ch=='\n')  i++;
  }
  fclose(fp);
  return i;
}

void append_file(char *file_path_from, char *file_path_to){
  int max_line_len=1000; 
  FILE *f_from;
  FILE *f_to;
  char buf[max_line_len+1];  

  /* open the source and the target files. */
  f_from = fopen(file_path_from, "r");
  if ( f_from == NULL ) {cart_error("Cannot open source file: %s",file_path_from);}
  f_to = fopen(file_path_to, "a");
  if ( f_to == NULL ) {cart_error("Cannot open target file: %s",file_path_to);}

  /* copy source file to target file, line by line. */
  while (fgets(buf, max_line_len+1, f_from)) {
    if (fputs(buf, f_to) == EOF) {  cart_error("Error writing to target file: %s",file_path_to);}
  }
  /* fgets failed _not_ due to encountering EOF */
  if (!feof(f_from)) { cart_error("Error reading from source file: %s",file_path_from);}

  /* close source and target file streams. */
  if (fclose(f_from) == EOF) {cart_error("Error when closing source file: %s",file_path_from);}
  if (fclose(f_to) == EOF) {cart_error("Error when closing target file: %s",file_path_to);}
}

void copy_file(char *file_path_from, char *file_path_to){
  int max_line_len=1000; 
  FILE *f_from;
  FILE *f_to;
  char buf[max_line_len+1];  

  /* open the source and the target files. */
  f_from = fopen(file_path_from, "r");
  if (!f_from) {cart_error("Cannot open source file: %s",file_path_from);}
  f_to = fopen(file_path_to, "w");
  if (!f_to) {cart_error("Cannot open target file: %s",file_path_to);}

  /* copy source file to target file, line by line. */
  while (fgets(buf, max_line_len+1, f_from)) {
    if (fputs(buf, f_to) == EOF) {  cart_error("Error writing to target file: %s",file_path_to);}
  }
  /* fgets failed _not_ due to encountering EOF */
  if (!feof(f_from)) { cart_error("Error reading from source file: %s",file_path_from);}

  /* close source and target file streams. */
  if (fclose(f_from) == EOF) {cart_error("Error when closing source file: %s",file_path_from);}
  if (fclose(f_to) == EOF) {cart_error("Error when closing target file: %s",file_path_to);}
}


#endif //LOG_STAR_CREATION
#endif //STARFORM

void init_logging( int restart ) {
	int i;
	char mode[2];
	char filename[256];

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
#ifdef COSMOLOGY
			fprintf(steptimes,"# step tl dtl auni abox dabox\n");
#else
			fprintf(steptimes,"# step tl dtl\n");
#endif /* COSMOLOGY */
		}

		sprintf(filename, "%s/energy.log", logfile_directory );
		energy = fopen(filename,mode);

		if ( energy == NULL ) {
			cart_error("Unable to open %s for writing!", filename );
		}

		if ( !restart || restart == 2 ) {
#ifdef COSMOLOGY
			fprintf(energy, "# step tl a gas_thermal gas_kinetic gas_potential total_gas_energy baryon_mass particle_kinetic particle_potential total_particle_energy error\n");
#else
			fprintf(energy, "# step tl gas_thermal gas_kinetic gas_potential total_gas_energy baryon_mass particle_kinetic particle_potential total_particle_energy error\n");
#endif /* COSMOLOGY */
		}

#ifdef STARFORM
		sprintf( filename, "%s/sf.log", logfile_directory );
		star_log = fopen(filename, mode);

		if ( star_log == NULL ) {
			cart_error("Unable to open %s for writing!", filename );
		}

		if ( !restart || restart == 2 ) {
#ifdef COSMOLOGY
			fprintf(star_log, "# step t dt a t [Gyrs] dt [yrs] N* M* dM* Mi* dMi* [Msun] SFR [Msun/yr/Mpc^3]\n");  
#else
			fprintf(star_log, "# step t dt t [Gyrs] dt [yrs] N* M* dM* Mi* dMi* [Msun] SFR [Msun/yr/Mpc^3]\n");  
#endif /* COSMOLOGY */
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
		fprintf( timing, "# %u levels %u timers\n", num_refinement_levels+2, NUM_TIMERS-1 );
#ifdef COSMOLOGY
		fprintf( timing, "# step t a total_time" );
#else
		fprintf( timing, "# step t total_time" );
#endif /* COSMOLOGY */
		for ( i = 1; i < NUM_TIMERS; i++ ) {
			fprintf( timing, " %s", timer_name[i][0] );
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
#ifdef COSMOLOGY
		fprintf( workload, "# step t a num_local_particles max_level" );
#else
		fprintf( workload, "# step t num_local_particles max_level" );
#endif /* COSMOLOGY */
		for ( i = min_level; i <= max_level; i++ ) {
			fprintf( workload, " num_local_cells[%u] num_buffer_cells[%u]", i, i );
		}
#ifdef DEBUG_MEMORY_USE
		fprintf( workload, " allocated_memory" );
#endif /* DEBUG_MEMORY_USE */
		fprintf( workload, "\n" );
		fflush(workload);
	}

	sprintf(filename, "%s/dependency.%03u.dat", logfile_directory, local_proc_id );
	dependency = fopen( filename, "a" );

	if ( dependency == NULL ) {
		cart_error( "Unable to open %s for writing!", filename );
	}


#ifdef DEBUG
	debug_breakpoint(-1,0,__FILE__,__LINE__);
#endif
}

void finalize_logging() {
#ifdef MPE_LOG
	char filename[256];
	sprintf( filename, "%s/mpe_log.dat", logfile_directory );
	MPE_Finish_log(filename);
#endif

#ifdef LOG_STAR_CREATION
	log_star_creation(-1,-1.0,FILE_CLOSE);
#endif

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
	ap1 = abox_from_tcode( tl[min_level] - 0.5*dtl[min_level] );
	da = abox[min_level] - abox_old[min_level];
#else
	ap1 = 1.0;
	ap0 = 0.0;
	da = 0.0;
#endif
#endif
	
	/* log profiling information */
#ifdef COSMOLOGY
	fprintf( timing, "%u %e %e %e", step, tl[min_level], auni[min_level], current_time( TOTAL_TIME, min_level-1 ) );
#else
	fprintf( timing, "%u %e %e", step, tl[min_level], current_time( TOTAL_TIME, min_level-1 ) );
#endif /* COSMOLOGY */
	for ( level = min_level-1; level <= max_level; level++ ) {
		for ( i = 1; i < NUM_TIMERS; i++ ) {
			fprintf( timing, " %e", total_time(i, level) );
		}
	}
	fprintf( timing, "\n" );
	fflush(timing);

	/* log workload information */
#ifdef PARTICLES
#ifdef COSMOLOGY
	fprintf( workload, "%u %e %e %u %u", step, tl[min_level], auni[min_level], num_local_particles, max_level_now() );
#else
	fprintf( workload, "%u %e %u %u", step, tl[min_level], num_local_particles, max_level_now() );
#endif /* COSMOLOGY */
#else
#ifdef COSMOLOGY
	fprintf( workload, "%u %e %e 0 %u", step, tl[min_level], auni[min_level], max_level_now() );
#else
	fprintf( workload, "%u %e 0 %u", step, tl[min_level], max_level_now() );
#endif /* COSMOLOGY */
#endif /* PARTICLES */

	for ( i = min_level; i <= max_level; i++ ) {
		fprintf(workload, " %u %u", num_cells_per_level[i], num_buffer_cells[i] );
	}
#ifdef DEBUG_MEMORY_USE
	fprintf( workload, " %lu",dmuReportAllocatedMemory());
#endif /* DEBUG_MEMORY_USE */
	fprintf(workload, "\n");
	fflush(workload);

	/* log dependency information */
#ifdef COSMOLOGY
	fprintf( dependency, "%u %e %e", step, tl[min_level], auni[min_level] );
#else
	fprintf( dependency, "%u %e", step, tl[min_level] );
#endif /* COSMOLOGY */
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
				gas_thermal += cell_volume[level]*cell_gas_pressure(icell)/(constants->gamma-1);

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
	dtyears = dtl[min_level] * units->time / constants->yr;
#ifdef COSMOLOGY
	current_age = tphys_from_abox(abox[min_level]);
#else
	current_age = tl[min_level];
#endif
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
#ifdef COSMOLOGY
		total_resolved_volume *= pow(units->length/constants->Mpc/abox[min_level],3.0);
#else
		total_resolved_volume *= pow(units->length/constants->Mpc,3.0);
#endif /* COSMOLOGY */


		cart_debug("total_resolved_volume = %e", total_resolved_volume);

		if ( total_resolved_volume > 0.0 ) {
			sfr = d_stellar_initial_mass * units->mass / constants->Msun / dtyears / total_resolved_volume;
		} else {
			sfr = 0.0;
		}

#ifdef COSMOLOGY
		fprintf( star_log, "%u %e %e %e %e %e %u %e %e %e %e %e\n", step, tl[min_level],
			dtl[min_level], auni[min_level], current_age, dtyears, 
#else
		fprintf( star_log, "%u %e %e %e %e %u %e %e %e %e %e\n", step, tl[min_level],
			dtl[min_level], current_age, dtyears, 
#endif /* COSMOLOGY */
			particle_species_num[num_particle_species-1],
			total_stellar_mass*units->mass/constants->Msun, 
			d_stellar_mass*units->mass/constants->Msun,
			total_stellar_initial_mass*units->mass/constants->Msun,
			d_stellar_initial_mass*units->mass/constants->Msun, sfr );

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
			au0 = abox_old[min_level]*total_particle_potential;
			aeu0 = au0 + abox_old[min_level]*total_particle_kinetic;
			tintg = 0.0;
			error = 0.0;
		} else {
			tintg += 2.0*(abox_old[min_level] - ap0)*ekin;
			kinetic_energy = (ekin*(ap1-ap0-0.5*da) + ekin1*0.5*da)/(ap1-ap0);
			error = ( abox_old[min_level]*(kinetic_energy+total_particle_potential) - aeu0 + tintg ) /
					( abox_old[min_level]*total_particle_potential - au0 );
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
#endif /* COSMOLOGY */
#else
		error = 0.0;
#endif /* PARTICLES */

#ifdef COSMOLOGY
		fprintf(energy, "%u %e %e %e %e %e %e %e %e %e %e %e\n",
			step, tl[min_level], auni[min_level], 
#else
		fprintf(energy, "%u %e %e %e %e %e %e %e %e %e %e\n",
			step, tl[min_level], 
#endif /* COSMOLOGY */
			total_gas_thermal, total_gas_kinetic, total_gas_potential, 
			total_gas_thermal+total_gas_kinetic+total_gas_potential,
			total_gas_mass,
			total_particle_kinetic, total_particle_potential,
			total_particle_kinetic+total_particle_potential,
			error );
		fflush(energy);

#ifdef COSMOLOGY
		fprintf(steptimes, "%u %e %e %e %e %e\n", step, tl[min_level], dtl[min_level], auni[min_level], abox[min_level], abox[min_level]-abox_old[min_level] );
#else
		fprintf(steptimes, "%u %e %e\n", step, tl[min_level], dtl[min_level] );
#endif /* COSMOLOGY */
		fflush(steptimes);
	}

#ifdef PARTICLES
	ekin1 = total_particle_kinetic;
	ap0 = ap1;
#endif

}


#ifdef DEBUG

double offset = 0.0;
unsigned long record = 0;

void debug_breakpoint(int timerid, int start, const char *file, int line)
{
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
	    fprintf(f,"%10lu: %s @ %d: %s done at %f sec.\n",record++,file,line,timer_name[timerid][0],MPI_Wtime()-offset);
	    break;
	  }
	case 1:
	  {
	    fprintf(f,"%10lu: %s @ %d: %s started at %f sec.\n",record++,file,line,timer_name[timerid][0],MPI_Wtime()-offset);
	    break;
	  }
	default:
	  {
	    fprintf(f,"%10lu: %s @ %d: marker #%d set at %f sec.\n",record++,file,line,timerid,MPI_Wtime()-offset);
	  }
	}
      fclose(f);
    }
}
#endif
