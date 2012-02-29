#include "config.h"


#include "auxiliary.h"
#include "parallel.h"
#include "timing.h"

#ifdef MPE_LOG
#include <mpe.h>
#endif


timer timers[num_refinement_levels+2][NUM_TIMERS];
int current_timer_level;

#ifdef DEBUG_TIMING
double timing_synch_precision = 0.01;  /* fraction precision */
double timing_synch_gap_width = 1;     /* absolute difference in seconds */

int report_gaps = 0;                   /* need to block gap reporting between levels - those gaps always exist and are ok */

int balanced_timers_ids[] = { WORK_TIMER, COMMUNICATION_TIMER, LOWER_LEVEL_TIMER };
const int num_balanced_timers = sizeof(balanced_timers_ids)/sizeof(int);

#endif

#ifdef MPE_LOG
int event[2*NUM_TIMERS];
#endif 

const char *timer_name[][2] = {
{ "total", "blue" }, 
{ "init", "blue" },
{ "restart", "blue" }, 
{ "output", "orange" },
{ "io", "orange" }, 
{ "particle_write_io", "orange" },
{ "particle_read_io", "orange" },
{ "gas_write_io", "orange" },
{ "gas_read_io", "orange" },
{ "load_balance", "blue" }, 
{ "load_balance_communication", "red" },
{ "choose_timestep", "blue" },
{ "choose_timestep_communication", "red" },
{ "hydro", "blue" }, 
{ "hydro_update", "red" }, 
{ "cooling", "blue" },
{ "accel_parts", "blue" },
{ "move_parts", "blue" }, 
{ "stellar_feedback_update", "red" },
{ "update_parts", "blue" }, 
{ "update_parts_communication", "red" },
{ "trade_particles", "red" },
{ "trade_particle_timer", "red" },
{ "gravity", "blue" }, 
{ "fft", "blue" },
{ "fft_communication", "red" },
{ "fft_update", "red" },
{ "prolongate_update", "red" },
{ "restrict_update", "red" },
{ "smooth", "blue" }, 
{ "smooth_setup", "blue" }, 
{ "smooth_update", "red" },
{ "smooth_communication", "red" }, 
{ "hydro_accel", "blue" },
{ "hydro_accel_update", "red" },
{ "particle_accel", "blue" },
{ "particle_accel_update", "red" },
{ "density", "blue" },
{ "merge_density", "red" },
{ "merge_density_update", "red" },
{ "merge_density_communication", "red" },
{ "refinement", "blue" },
{ "diffusion_step", "blue" },
{ "diffusion_update", "red" },
{ "modify_update", "red" },
{ "derefine_update", "red" },
{ "split_buffer", "red" },
{ "split_buffer_communication", "red" },
{ "join_buffer", "red" }, 
{ "join_buffer_communication", "red" },
{ "update", "red" }, 
{ "update_send", "red" }, 
{ "update_recv", "red" },
{ "build_cell_buffer", "green" }, 
{ "select_level", "green" }, 
{ "max_level", "red" },
{ "communication", "red" },
{ "lower_level", "white" },  /* this is for internal accouting only */
#ifdef USER_PLUGIN
{ "plugin", "gray" },
#endif
{ "work", "green" },
#ifdef RADIATIVE_TRANSFER
{ "RT_tables", "green" },
{ "RT_cooling", "green" },
{ "RT_level_update", "red" },
{ "RT_global_update", "blue" },
{ "RT_after_density", "blue" },
{ "RT_tree_emulator_update", "red" },
{ "RT_single_source_update", "red" },
{ "RT_solve_equation_update", "red" },
#endif
{ "level_total", "blue" }
};

void init_timers() {
	int i, j;

	for ( i = min_level-1; i <= max_level; i++ ) {
		for ( j = 0; j < NUM_TIMERS; j++ ) {
			timers[i+1][j].num_calls = 0;
			timers[i+1][j].current_time = -1.0;
			timers[i+1][j].total_time = 0.0;
			timers[i+1][j].last_time = 0.0;
#ifdef DEBUG_TIMING 
			timers[i+1][j].last_wtime = 0.0;
#endif
		}
	}

	current_timer_level = min_level-1;

#ifdef MPE_LOG
	MPE_Init_log();
	for ( i = 0; i < NUM_TIMERS; i++ ) {
        MPE_Log_get_state_eventIDs( &event[2*i], &event[2*i+1]);
		if ( local_proc_id == MASTER_NODE ) {
			MPE_Describe_state( event[2*i], event[2*i+1], timer_name[i][0], timer_name[i][1] );
		}
	}
	MPI_Barrier(mpi.comm.run);
	MPE_Start_log();
#endif /* MPE_LOG */
}

void start_time_at_location( int timerid, const char *file, int line ) {
	double wtime = MPI_Wtime();
#ifdef DEBUG_TIMING
	double d;
	timer *prev;
	int i, check_balance;
#endif

	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	if ( timers[current_timer_level+1][timerid].current_time != -1.0 ) {
		cart_error("Timer already started @ %s:%u: level = %d, timerid = %u", 
				file, line, current_timer_level, timerid );
	}

	timers[current_timer_level+1][timerid].current_time = wtime;

#ifdef DEBUG_TIMING

	if(timerid==WORK_TIMER && timers[current_timer_level+1][COMMUNICATION_TIMER].current_time>0.0)
	  {
	    cart_error("TIMING: WORK and COMMUNICATION overlap at level %d in %s:%u",current_timer_level,file,line);
	  }

	if(timerid==COMMUNICATION_TIMER && timers[current_timer_level+1][WORK_TIMER].current_time>0.0)
	  {
	    cart_error("TIMING: COMMUNICATION and WORK overlap at level %d in %s:%u",current_timer_level,file,line);
	  }

	/*
	//  Check for gaps inside a step
	*/
	if(report_gaps && current_timer_level>-1 && timers[current_timer_level+1][LEVEL_TIMER].current_time>0.0)
	  {
	    check_balance = 0;
	    for(i=0; i<num_balanced_timers; i++) if(timerid == balanced_timers_ids[i])
	      {
		check_balance = 1;
		break;
	      }

	    if(check_balance)
	      {
		prev = &timers[current_timer_level+1][balanced_timers_ids[0]];
		for(i=1; i<num_balanced_timers; i++) if(prev->last_wtime < timers[current_timer_level+1][balanced_timers_ids[i]].last_wtime)
		  {
		    prev = &timers[current_timer_level+1][balanced_timers_ids[i]];
		  }

		d = wtime - prev->last_wtime;
		if(prev->last_time>0.0 && d>timing_synch_gap_width && d>timing_synch_precision*prev->last_time)
		  {
		    cart_debug("TIMING: gap from (%s:%u) to (%s:%u) [%lg s,%lg %]",prev->last_file,prev->last_line,file,line,d,100*d/prev->last_time);
		  }
	      }
	  }

#endif /* DEBUG_TIMING */

#ifdef MPE_LOG
	MPE_Log_event( event[2*timerid], current_timer_level, timer_name[timerid][0] );
#endif

#ifdef DEBUG
	debug_breakpoint(timerid,1,file,line);
#endif
}

double end_time_at_location( int timerid, const char *file, int line ) {
	double elapsed;
	double wtime = MPI_Wtime();
#ifdef DEBUG_TIMING
	int i;
#endif

	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( timers[current_timer_level+1][timerid].current_time > 0.0 );

	elapsed = wtime - timers[current_timer_level+1][timerid].current_time;

#ifdef MPE_LOG
	MPE_Log_event( event[2*timerid+1], current_timer_level, timer_name[timerid][0] );
#endif /* MPE_LOG */

	timers[current_timer_level+1][timerid].last_time = elapsed;	
	timers[current_timer_level+1][timerid].total_time += elapsed;
	timers[current_timer_level+1][timerid].current_time = -1.0;
	timers[current_timer_level+1][timerid].num_calls++;
#ifdef DEBUG_TIMING
	timers[current_timer_level+1][timerid].last_wtime = wtime;	
	timers[current_timer_level+1][timerid].last_file = file;
	timers[current_timer_level+1][timerid].last_line = line;

	for(i=0; i<num_balanced_timers; i++) if(timerid == balanced_timers_ids[i])
	  {
	    report_gaps = 1;
	    break;
	  }
#endif /* DEBUG_TIMING */

#ifdef DEBUG
	debug_breakpoint(timerid,0,file,line);
#endif

	return elapsed;
}

double current_time( int timerid, int level ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( level >= min_level-1 && level <= max_level );
	cart_assert( timers[level+1][timerid].current_time > 0.0 );

	return MPI_Wtime() - timers[level+1][timerid].current_time;
}

double average_time( int timerid, int level ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( level >= min_level-1 && level <= max_level );

	return timers[level+1][timerid].total_time / (double) timers[level+1][timerid].num_calls;
}

double total_time( int timerid, int level ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( level >= min_level-1 && level <= max_level );
	return timers[level+1][timerid].total_time;
}

double last_time( int timerid, int level ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( level >= min_level-1 && level <= max_level );
	return timers[level+1][timerid].last_time;
}

void start_timing_level( int level ) {
	current_timer_level = level;
#ifdef DEBUG_TIMING
	report_gaps = 0;
#endif /* DEBUG_TIMING */
}

void end_timing_level( int level ) {
	cart_assert( level >= min_level && level <= max_level );
	current_timer_level = level-1;
}


#ifdef DEBUG

#include <stdio.h>
#include "io.h"

double offset = 0.0;
unsigned long record = 0;
extern int current_step_level;

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
	    fprintf(f,"%10lu/L=%2d: %s @ %d: %s done at %f sec.\n",record++,current_step_level,file,line,timer_name[timerid][0],MPI_Wtime()-offset);
	    break;
	  }
	case 1:
	  {
	    fprintf(f,"%10lu/L=%2d: %s @ %d: %s started at %f sec.\n",record++,current_step_level,file,line,timer_name[timerid][0],MPI_Wtime()-offset);
	    break;
	  }
	default:
	  {
	    fprintf(f,"%10lu/L=%2d: %s @ %d: marker #%d set at %f sec.\n",record++,current_step_level,file,line,timerid,MPI_Wtime()-offset);
	  }
	}
      fclose(f);
    }
}
#endif
