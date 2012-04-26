#include "config.h"

#include <math.h>
#include <stdio.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "io.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "starformation.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "logging.h"
#include "step.h"

#ifdef MPE_LOG
#include <mpe.h>
#endif
#ifdef LOG_STAR_CREATION
#include "starformation_recipes.h"
#endif

FILE *steptimes;
FILE *levelsteptimes;
FILE *timing;
FILE *energy;
FILE *workload;
FILE *dependency;

#ifdef STAR_FORMATION
FILE *star_log;
#endif /* STAR_FORMATION */

#ifdef STAR_FORMATION
#ifdef LOG_STAR_CREATION

void output_star_creation( int icell, double mass, FILE *f ){
  int level_UP, icell_UP;
  double side_UP;
  double side;
  int level;
  int id;
  int i;
  double zSol_cell,zSol_cell_UP;
  int high_levels,children_at_level, k,ic,icell0,icid,icicell,ilev,cells_below_ilev, index_to_child;
  int icipar;
  double sfrbelow_UP;
  double tem_max , rho_min ;
  double sf_max_gas_temperature = 2.0e4;
  tem_max = sf_max_gas_temperature/(constants->wmu*units->temperature); 
  rho_min = 0.5/(constants->XH*units->number_density);

  //true, but useless id = last_star_id + local_proc_id + 1;

  level = cell_level(icell);
  level_UP = ( log( cell_size[min_level] / (1.0*constants->kpc/units->length) ) / log(2.0) ); //want UP to be ~1pKpc
    //level_UP = level-1;
  if(level_UP<=min_level){level_UP=min_level;}
  icell_UP = icell;
  if( level > level_UP ){
    for( i = level; i > level_UP; i-- ){
      icell_UP = cell_parent_cell(icell_UP);
    }
  } else {
    level_UP = level;
  }
  side    = cell_size[ level    ]*(units->length/constants->kpc);
  side_UP = cell_size[ level_UP ]*(units->length/constants->kpc);
 


  high_levels = max_level_now() - level_UP  ; 
  children_at_level = pow( 2, 3*high_levels );
  sfrbelow_UP = 0; 
  ic = 0;
  icell0 = icell_UP ;
  while ( ic < children_at_level ){ //for each of the children that could potentially exist at max_level
    icid = ic;
    icicell = icell0;
    ilev = 0; 
    while ( ilev < high_levels  && icicell !=-1 ){ //find where the child lives in each higher level, and progress downward as long as refined
      cells_below_ilev = pow(  2 , 3*( high_levels-(ilev+1) ) );  //below means ilev+1
      index_to_child = (int)icid / cells_below_ilev ;
      icid -= cells_below_ilev * index_to_child;
      icipar=icicell;
      icicell = cell_child ( icicell, index_to_child );
      ilev++;
    }
    ilev--; // ilev should not have been incremented if you are not going to a child
    if ( icicell != -1 ) {  //if you reach maxlevel add it
      //cart_debug("snl %d %d %d ",ic,icicell, cell_level(icicell));
        sfrbelow_UP += sf_recipe->rate(icicell)*cell_volume[cell_level(icicell)]; 
      ic++;
    }else { //if you dont reach maxlevel then the leaf is not at maxlevel, skip ic to the next cell at this level
      //cart_debug("snl %d %d %d ",ic,icipar, cell_level(icipar) );
      if( cell_gas_density(icipar)>rho_min && cell_gas_pressure(icipar)/cell_gas_density(icipar)<tem_max){
          sfrbelow_UP += sf_recipe->rate(icipar)*cell_volume[cell_level(icipar)]; 
      }//*pow( 2,3*(high_levels-ilev));
      ic += pow( 2,3*(high_levels-ilev) ) ; 
    }
  }
  sfrbelow_UP=sfrbelow_UP/cell_volume[level_UP];
  zSol_cell_UP = cell_gas_metal_density(icell_UP)/(constants->Zsun*cell_gas_density(icell_UP));
  zSol_cell_UP = MAX(1.0e-3,zSol_cell_UP);
  zSol_cell = cell_gas_metal_density(icell)/(constants->Zsun*cell_gas_density(icell));
  zSol_cell = MAX(1.0e-3,zSol_cell);
  
  //abox, auni, tcosmo, level, side(pkpc), over_density, rho_H[msun/pc3], Z, temp, SFRD Msun/pkpc^3, mass, icell_UP,levelUP, sideUP(pkpc)...temp
  fprintf(f,"%.10lf %.6lf %e  %d %f  %e %e %e %e  %e | %e %d   %d %f  %e %e  %e %e %e  %e\n",
	  abox[level] , auni[level], tphys_from_tcode(tl[level]),
	  level, side,
	  cell_total_mass(icell), constants->XH*cell_gas_density(icell)*units->density/(constants->Msun/pow(constants->pc,3)) ,
	  zSol_cell,
	  cell_gas_pressure(icell)/cell_gas_density(icell)*(constants->wmu*units->temperature), //snl
	  sf_recipe->rate(icell)/(units->time/constants->yr)* units->density/(constants->Msun/pow(constants->kpc,3)) ,

	  mass*units->mass/constants->Msun , icell_UP, 
	  level_UP, side_UP ,
	  cell_total_mass(icell_UP),  constants->XH*cell_gas_density(icell_UP)*units->density/(constants->Msun/pow(constants->pc,3)) , 
	  sfrbelow_UP/(units->time/constants->yr)* units->density/(constants->Msun/pow(constants->kpc,3)), 
	  zSol_cell_UP,
	  sf_recipe->rate(icell_UP)/(units->time/constants->yr)* units->density/(constants->Msun/pow(constants->kpc,3))
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
  /* Don't want to restart from jobname.d and have a corrupt/empty/inaccurate restart_sc file.
   * Match the number of stars since the last output ae in the code to the number of stars listed in screst.
   */
  double last_ae;
  double fstar_ae;
  int local_count_stars,count_recent_stars; 
  int low_local_count_stars,low_count_recent_stars; 
  int high_local_count_stars,high_count_recent_stars; 
  int count_stars;
  int nsclog, iout,i;
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
  MPI_Allreduce( &local_count_stars, &count_stars, 1, MPI_INT, MPI_SUM, mpi.comm.run );
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
	    i = sscanf(str_buf, "%lf", &last_ae);
	    cart_debug("found last_ae=%f restart_sc file",last_ae);
	  }else {
	    cart_error("restart_sc must contain last_ae if it has not just been initialized after a labeled output or no stars exist");
	  }
	  fclose(fp);
	}
      }
    }
    /* check last_ae is inside relevant interval */
    if( (last_ae < outputs[iout] || last_ae > auni[min_level]+1e-6) && iout > 0 ){
      cart_error("bad last_ae: %f < %f || %f > %f",last_ae,outputs[iout],last_ae,auni[min_level]);
    }
    cart_debug("star creation outputs check: iout=%d a=%f last_ae=%f outi=%f outi+1=%f",i, auni[min_level],last_ae,outputs[iout],outputs[iout+1]);
  }
  MPI_Bcast(&last_ae, 1, MPI_DOUBLE, MASTER_NODE, mpi.comm.run );

  if(last_ae != auni[min_level] && count_stars > 0){
    cart_debug("counting recent star creations after last_ae=%f (low=%f, high=%f)", last_ae, last_ae+1e-6 , last_ae-1e-6 );
    /*count stars since the last restart*/ 
    local_count_stars = 0; count_recent_stars = 0;
    low_local_count_stars = 0; low_count_recent_stars = 0;
    high_local_count_stars = 0; high_count_recent_stars = 0;
    for ( i = 0; i < num_star_particles; i++ ) {
      if ( particle_level[i] != FREE_PARTICLE_LEVEL && particle_is_star(i) ) {
	/*precision of star_tbirth is low so stars scatter in and out of last_ae and count is bad*/
	fstar_ae = auni_from_tcode((double)star_tbirth[i]);
	if (  fstar_ae > last_ae ) {
	  local_count_stars++ ;
	  //cart_debug("snl %.10lf %d",auni_from_tcode((double)star_tbirth[i]),particle_id[i]);
	}
	if ( fstar_ae > last_ae + 1e-6 ) {
	  low_local_count_stars++ ;
	}
	if ( fstar_ae > last_ae - 1e-6 ) {
	  high_local_count_stars++ ;
	}
      }
    }
    MPI_Reduce( &local_count_stars, &count_recent_stars, 1, MPI_INT, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( &low_local_count_stars, &low_count_recent_stars, 1, MPI_INT, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( &high_local_count_stars, &high_count_recent_stars, 1, MPI_INT, MPI_SUM, MASTER_NODE, mpi.comm.run );
    
    /* make sure your restart_sc file matches the number of stars since the last backup to a labeled file */
    if ( local_proc_id == MASTER_NODE ) {
      if( count_recent_stars > 0 ){
	nsclog = count_lines( filename_screst );
	nsclog = nsclog - 1;// one line is a header
	if( nsclog != count_recent_stars ){
	  cart_debug( "star count in restart_star_creation.log and star_tbirth don't match: (log:%d) != (dst:%d)",nsclog,count_recent_stars);
	  if( nsclog < low_count_recent_stars || nsclog > high_count_recent_stars ){
	    cart_error( "even accounting for low precision star_tbirth there is an error in the file: (log:%d)  < (dst_low:%d) || > (dst_high:%d)",nsclog,low_count_recent_stars,high_count_recent_stars );
	  }else{
	    cart_debug( "accounting for low precision star_tbirth: (dst_low:%d) <= (log:%d) <= (dst_high:%d)",low_count_recent_stars,nsclog,high_count_recent_stars );
	  }
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
  fprintf(fp,"%.16lf\n",last_ae);
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
#endif //STAR_FORMATION

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
			fprintf(steptimes,"# step tl [Gyr] dtl [Myr] auni abox dabox\n");
#else
			fprintf(steptimes,"# step tl dt\n");
#endif /* COSMOLOGY */
		}

		sprintf(filename,"%s/level_times.log", logfile_directory );
		levelsteptimes = fopen(filename,mode);

		if ( levelsteptimes == NULL ) {
			cart_error("Unable to open %s for writing!", filename );
		}

		if ( !restart || restart == 2 ) {
#ifdef COSMOLOGY
			fprintf(levelsteptimes,"# step dtl [Myr] dtl[levels] [code units]\n");
#else
			fprintf(levelsteptimes,"# step dt dtl[levels] [code units]\n");
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

#ifdef STAR_FORMATION
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
#endif /* STAR_FORMATION */
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
		fclose(levelsteptimes);
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
	double dtyears, current_age, current_dt;

#ifdef STAR_FORMATION
	double stellar_mass, stellar_initial_mass;
	double old_stellar_mass, old_stellar_initial_mass;
	double d_stellar_mass, d_stellar_initial_mass;
	double resolved_volume[max_level-min_level+1];
	double local_resolved_volume[max_level-min_level+1];
	double total_resolved_volume;
	double sfr;
#endif /* STAR_FORMATION */

	dtyears = dtl[min_level]*units->time/constants->yr;
#ifdef COSMOLOGY
	current_age = tphys_from_abox(abox[min_level]);
	current_dt = dtyears;
#else
	current_age = tl[min_level]*units->time;
	current_dt = dtl[min_level]*units->time;
#endif

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
#ifdef GRAVITY
				gas_potential += cell_gas_density(icell)*cell_volume[level]*cell_potential(icell);
#endif
				gas_mass += cell_gas_density(icell)*cell_volume[level];
			}
		}
		cart_free( level_cells );
	}

	/* add stellar mass to gas mass */
#ifdef STAR_FORMATION
	stellar_mass = 0.0;
	stellar_initial_mass = 0.0;
	for ( i = 0; i < num_star_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL && particle_is_star(i) ) {
			gas_mass += particle_mass[i];
			stellar_mass += particle_mass[i];
			stellar_initial_mass += star_initial_mass[i];
		}
	}
#endif /* STAR_FORMATION */

	MPI_Reduce( &gas_thermal, &total_gas_thermal, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
	MPI_Reduce( &gas_kinetic, &total_gas_kinetic, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
#ifdef GRAVITY
	MPI_Reduce( &gas_potential, &total_gas_potential, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
#endif /* GRAVITY */
	MPI_Reduce( &gas_mass, &total_gas_mass, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
#endif /* HYDRO */

#ifdef STAR_FORMATION
	old_stellar_mass = total_stellar_mass;
	old_stellar_initial_mass = total_stellar_initial_mass;

	MPI_Reduce( &stellar_mass, &total_stellar_mass, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
	MPI_Reduce( &stellar_initial_mass, &total_stellar_initial_mass, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );

	d_stellar_mass = MAX( 0.0, total_stellar_mass - old_stellar_mass );
	d_stellar_initial_mass = MAX( 0.0, total_stellar_initial_mass - old_stellar_initial_mass );

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

	MPI_Reduce( local_resolved_volume, resolved_volume, max_level-min_level+1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );

	if ( local_proc_id == MASTER_NODE ) {
		/* sum resolved volume over all levels except min_level (works for MMR simulations) */
		total_resolved_volume = 0.0;
		for ( level = min_level; level <= max_level; level++ ) {
			total_resolved_volume += resolved_volume[level];
		}
#ifdef COSMOLOGY
		total_resolved_volume *= pow(units->length/constants->Mpc/abox[min_level],3.0);
#else
		total_resolved_volume *= pow(units->length/constants->Mpc,3.0);
#endif /* COSMOLOGY */

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
#endif /* STAR_FORMATION */


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

	MPI_Reduce( &particle_kinetic, &total_particle_kinetic, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
	MPI_Reduce( &particle_potential, &total_particle_potential, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );

	total_particle_kinetic *= 0.5;
	total_particle_potential *= 0.5;
#endif

	if ( local_proc_id == MASTER_NODE ) {
#ifdef PARTICLES
		ekin = total_particle_kinetic/(ap1*ap1);

#ifdef COSMOLOGY
		if ( step == 1 ) {
			au0 = abox_old[min_level]*total_particle_potential;
			aeu0 = au0 + abox_old[min_level]*ekin;
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
			aeu0 = au0 + ekin;
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
		fprintf(steptimes, "%u %e %e %e %e %e\n", 
				step, current_age, current_dt, auni[min_level], abox[min_level], 
				abox[min_level]-abox_old[min_level] );
#else
		fprintf(steptimes, "%u %e %e\n", step, current_age, current_dt );
#endif /* COSMOLOGY */
		fflush(steptimes);

		fprintf(levelsteptimes, "%u %e", step, current_dt );
		for ( i = min_level; i <= max_level; i++ ) {
			fprintf(levelsteptimes, " %e", dtl[i] );
		}
		fprintf(levelsteptimes, "\n" );
		fflush(levelsteptimes);
	}

#ifdef PARTICLES
	ekin1 = ekin;
	ap0 = ap1;
#endif

}

