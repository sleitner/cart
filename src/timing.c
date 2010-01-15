#include "config.h"

#include <mpi.h>

#ifdef MPE_LOG
#include <mpe.h>
#endif

#include "auxiliary.h"
#include "logging.h"
#include "timing.h"

timer timers[num_refinement_levels+1][NUM_TIMERS];
int current_timer_level;

#ifdef MPE_LOG
int event[2*NUM_TIMERS];
#endif 

const char *timer_name[][2] = {
{	"total", "blue" }, 
{	"init", "blue" },
{	"restart", "blue" }, 
{	"io", "orange" }, 
{	"load_balance", "blue" }, 
{	"output", "orange" }, 
{	"gravity", "blue" }, 
{	"particle_accel", "blue" }, 
{	"hydro", "blue" }, 
{	"hydro_update", "red" }, 
{	"hydro_accel", "blue" }, 
{	"density", "blue" }, 
{	"move_parts", "blue" }, 
{	"update_parts", "blue" }, 
{	"trade_particles", "red" }, 
{	"refinement", "red" }, 
{	"smooth", "blue" }, 
{	"smooth_setup", "blue" }, 
{	"smooth_communication", "red" }, 
{	"fft", "blue" }, 
{	"update", "red" }, 
{	"update_send", "red" }, 
{	"update_recv", "red" },
{	"build_cell_buffer", "green" }, 
{	"diffusion_step", "blue" }, 
{	"split_buffer", "red" },
{	"join_buffer", "red" }, 
{	"select_level", "green" }, 
{	"select_level_optimize_level", "green" },
{	"max_level", "red" },
{	"hydro_accel_update", "red" },
{	"diffusion_update", "red" },
{	"modify_update", "red" },
{	"prolongate_update", "red" },
{	"smooth_update", "red" },
{	"restrict_update", "red" },
{	"fft_update", "red" },
{	"particle_accel_update", "red" },
{	"merge_densities_update", "red" },
{	"hydro_particle_update", "red" },
{	"work", "green" },
{	"communication", "red" },
{	"particle_write_io", "orange" },
{	"particle_read_io", "orange" },
{	"gas_write_io", "orange" },
{	"gas_read_io", "orange" },
{	"cooling", "blue" },
{	"merge_density", "red" },
#ifdef RADIATIVE_TRANSFER
{	"RT_tables", "green" },
{	"RT_cooling", "green" },
{	"RT_level_update", "blue" },
{	"RT_after_density", "blue" },
#endif
{	"level_total" "blue" },
};

void init_timers() {
	int i, j;

	for ( i = min_level; i <= max_level; i++ ) {
		for ( j = 0; j < NUM_TIMERS; j++ ) {
			timers[i][j].num_calls = 0;
			timers[i][j].current_time = -1.0;
			timers[i][j].total_time = 0.0;
			timers[i][j].last_time = 0.0;
		}
	}

	current_timer_level = min_level;

#ifdef MPE_LOG
	MPE_Init_log();
	for ( i = 0; i < NUM_TIMERS; i++ ) {
        MPE_Log_get_state_eventIDs( &event[2*i], &event[2*i+1]);
		MPE_Describe_state( event[2*i], event[2*i+1], timer_name[i][0], timer_name[i][1] );
	}
	MPE_Start_log();
#endif /* MPE_LOG */
}

void start_time_at_location( int timerid, const char *file, int line ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( timers[current_timer_level][timerid].current_time == -1.0 );

	timers[current_timer_level][timerid].current_time = MPI_Wtime();

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
	cart_assert( timers[current_timer_level][timerid].current_time > 0.0 );

	elapsed = MPI_Wtime() - timers[current_timer_level][timerid].current_time;

#ifdef MPE_LOG
	MPE_Log_event( event[2*timerid+1], current_timer_level, timer_name[timerid][0] );
#endif /* MPE_LOG */

	timers[current_timer_level][timerid].last_time = elapsed;	
	timers[current_timer_level][timerid].total_time += elapsed;
	timers[current_timer_level][timerid].current_time = -1.0;
	timers[current_timer_level][timerid].num_calls++;

#ifdef DEBUG
	log_in_debug(timerid,0,file,line);
#endif

	return elapsed;
}

double current_time( int timerid, int level ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( level >= min_level && level <= max_level );
	cart_assert( timers[level][timerid].current_time > 0.0 );

	return MPI_Wtime() - timers[level][timerid].current_time;
}

double average_time( int timerid, int level ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( level >= min_level && level <= max_level );

	return timers[level][timerid].total_time / (double) timers[level][timerid].num_calls;
}

double total_time( int timerid, int level ) {
	return timers[level][timerid].total_time;
}

double last_time( int timerid, int level ) {
	return timers[level][timerid].last_time;
}

void start_timing_level( int level ) {
	current_timer_level = level;
}

void end_timing_level( int level ) {
	current_timer_level = max( level-1, min_level );
}

