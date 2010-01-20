#include "config.h"

#include <mpi.h>

#ifdef MPE_LOG
#include <mpe.h>
#endif

#include "parallel.h"
#include "auxiliary.h"
#include "logging.h"
#include "timing.h"

timer timers[num_refinement_levels+2][NUM_TIMERS];
int current_timer_level;

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
{ "work", "green" },
#ifdef RADIATIVE_TRANSFER
{ "RT_tables", "green" },
{ "RT_cooling", "green" },
{ "RT_level_update", "red" },
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
	MPI_Barrier(MPI_COMM_WORLD);
	MPE_Start_log();
#endif /* MPE_LOG */
}

void start_time_at_location( int timerid, const char *file, int line ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	if ( timers[current_timer_level+1][timerid].current_time != -1.0 ) {
		cart_error("Timer already started @ %s:%u: level = %d, timerid = %u", 
				file, line, current_timer_level, timerid );
	}

	timers[current_timer_level+1][timerid].current_time = MPI_Wtime();

#ifdef MPE_LOG
	MPE_Log_event( event[2*timerid], current_timer_level, timer_name[timerid][0] );
#endif

#ifdef DEBUG
	log_in_debug(timerid,1,file,line);
#endif
}

double end_time_at_location( int timerid, const char *file, int line ) {
	double elapsed;
	
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( timers[current_timer_level+1][timerid].current_time > 0.0 );

	elapsed = MPI_Wtime() - timers[current_timer_level+1][timerid].current_time;

#ifdef MPE_LOG
	MPE_Log_event( event[2*timerid+1], current_timer_level, timer_name[timerid][0] );
#endif /* MPE_LOG */

	timers[current_timer_level+1][timerid].last_time = elapsed;	
	timers[current_timer_level+1][timerid].total_time += elapsed;
	timers[current_timer_level+1][timerid].current_time = -1.0;
	timers[current_timer_level+1][timerid].num_calls++;

#ifdef DEBUG
	log_in_debug(timerid,0,file,line);
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
}

void end_timing_level( int level ) {
	cart_assert( level >= min_level && level <= max_level );
	current_timer_level = level-1;
}

