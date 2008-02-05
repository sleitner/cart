#include <mpi.h>

#include "defs.h"
#include "tree.h"
#include "timing.h"
#include "auxiliary.h"

double timelimit = 0.0;
timer timers[num_refinement_levels+1][NUM_TIMERS];
int current_timer_level;

const char *timer_name[] = {
        "total",
        "init",
        "restart",
        "io",
        "load_balance",
        "output",
        "gravity",
	"particle_accel",
        "hydro",
	"hydro_update",
	"hydro_accel",
        "density",
	"merge_particle_densities",
        "move_parts",
        "update_parts",
	"trade_particles",
        "refinement",
        "smooth",
        "smooth_setup",
        "smooth_communication",
        "fft",
        "update",
	"update_send",
	"update_recv",
        "build_cell_buffer",
        "diffusion_step",
        "split_buffer",
        "join_buffer",
        "select_level",
        "select_level_optimize_level",
	"max_level",
	"hydro_accel_update",
	"diffusion_update",
	"modify_update",
	"prolongate_update",
	"smooth_update",
	"restrict_update",
	"fft_update",
	"particle_accel_update",
	"merge_densities_update",
	"hydro_particle_update",
	"work",
	"communication",
	"particle_io",
	"gas_io",
        "level_total"
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
}

void start_time( int timerid ) {
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( timers[current_timer_level][timerid].current_time == -1.0 );

	timers[current_timer_level][timerid].current_time = MPI_Wtime();
}

double end_time( int timerid ) {
	double elapsed;
	
	cart_assert( timerid >= 0 && timerid < NUM_TIMERS );
	cart_assert( timers[current_timer_level][timerid].current_time > 0.0 );

	elapsed = MPI_Wtime() - timers[current_timer_level][timerid].current_time;

	timers[current_timer_level][timerid].last_time = elapsed;	
	timers[current_timer_level][timerid].total_time += elapsed;
	timers[current_timer_level][timerid].current_time = -1.0;
	timers[current_timer_level][timerid].num_calls++;

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

