#ifndef __TIMING_H__
#define __TIMING_H__

typedef struct TIMER {
	int num_calls;
	double current_time;
	double last_time;
	double total_time;
} timer;

extern double timelimit;

#define		TOTAL_TIME			0
#define		INIT_TIMER			1
#define		RESTART_TIMER			2
#define		IO_TIMER			3
#define		LOAD_BALANCE_TIMER		4
#define		OUTPUT_TIMER			5
#define		GRAVITY_TIMER			6
#define		PARTICLE_ACCEL_TIMER		7
#define		HYDRO_TIMER			8
#define		HYDRO_UPDATE_TIMER		9
#define		HYDRO_ACCEL_TIMER		10
#define		DENSITY_TIMER			11
#define		MOVE_PARTS_TIMER		12
#define		UPDATE_PARTS_TIMER		13
#define		TRADE_PARTICLE_TIMER		14
#define		REFINEMENT_TIMER		15
#define		SMOOTH_TIMER			16
#define		SMOOTH_SETUP_TIMER		17
#define		SMOOTH_COMMUNICATION_TIMER	18
#define		FFT_TIMER			19
#define 	UPDATE_TIMER			20
#define		UPDATE_SEND_TIMER		21
#define		UPDATE_RECV_TIMER		22
#define		BUILD_CELL_BUFFER_TIMER		23
#define		DIFFUSION_STEP_TIMER		24
#define		SPLIT_BUFFER_TIMER		25
#define		JOIN_BUFFER_TIMER		26
#define		SELECT_LEVEL_TIMER		27
#define		SELECT_LEVEL_OPTIMIZE_LEVEL	28
#define		MAX_LEVEL_TIMER			29
#define		HYDRO_ACCEL_UPDATE_TIMER	30
#define		DIFFUSION_UPDATE_TIMER		31
#define		MODIFY_UPDATE_TIMER		32
#define		PROLONGATE_UPDATE_TIMER		33
#define		SMOOTH_UPDATE_TIMER		34
#define		RESTRICT_UPDATE_TIMER 		35
#define		FFT_UPDATE_TIMER		36
#define		PARTICLE_ACCEL_UPDATE_TIMER	37
#define		MERGE_DENSITIES_UPDATE_TIMER	38
#define		HYDRO_PARTICLE_UPDATE_TIMER	39
#define		WORK_TIMER                      40
#define		COMMUNICATION_TIMER             41
#define		PARTICLE_WRITE_IO_TIMER		42
#define		PARTICLE_READ_IO_TIMER		43
#define		GAS_WRITE_IO_TIMER			44
#define		GAS_READ_IO_TIMER			45
#define		COOLING_TIMER			46
#define		MERGE_DENSITY_TIMER		47

#ifdef RADIATIVE_TRANSFER

#define		RT_TABLES_TIMER			48
#define		RT_COOLING_TIMER		49
#define		RT_LEVEL_UPDATE_TIMER	        50
#define		RT_AFTER_DENSITY_TIMER		51
#define		LEVEL_TIMER			52

#else

#define		LEVEL_TIMER			48

#endif

#define		NUM_TIMERS			(LEVEL_TIMER+1)


#define start_time(timerid) start_time_at_location(timerid,__FILE__,__LINE__)
#define   end_time(timerid)   end_time_at_location(timerid,__FILE__,__LINE__)


extern const char *timer_name[];

void init_timers();
void start_timing_level( int level );
void end_timing_level( int level );
void start_time_at_location( int timerid, const char *file, int line );
double end_time_at_location( int timerid, const char *file, int line );

double current_time( int timerid, int level );
double average_time( int timerid, int level );
double total_time( int timerid, int level );
double last_time( int timerid, int level );

#endif
