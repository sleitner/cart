#ifndef __TIMING_H__
#define __TIMING_H__

typedef struct TIMER {
	int num_calls;
	double current_time;
	double last_time;
	double total_time;
} timer;

extern double timelimit;

#define		NUM_TIMERS			(46)

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
#define		MERGE_DENSITY_TIMER		12
#define		MOVE_PARTS_TIMER		13
#define		UPDATE_PARTS_TIMER		14
#define		TRADE_PARTICLE_TIMER		15
#define		REFINEMENT_TIMER		16
#define		SMOOTH_TIMER			17
#define		SMOOTH_SETUP_TIMER		18
#define		SMOOTH_COMMUNICATION_TIMER	19
#define		FFT_TIMER			20
#define 	UPDATE_TIMER			21
#define		UPDATE_SEND_TIMER		22
#define		UPDATE_RECV_TIMER		23
#define		BUILD_CELL_BUFFER_TIMER		24
#define		DIFFUSION_STEP_TIMER		25
#define		SPLIT_BUFFER_TIMER		26
#define		JOIN_BUFFER_TIMER		27
#define		SELECT_LEVEL_TIMER		28
#define		SELECT_LEVEL_OPTIMIZE_LEVEL	29
#define		MAX_LEVEL_TIMER			30
#define		HYDRO_ACCEL_UPDATE_TIMER	31
#define		DIFFUSION_UPDATE_TIMER		32
#define		MODIFY_UPDATE_TIMER		33
#define		PROLONGATE_UPDATE_TIMER		34
#define		SMOOTH_UPDATE_TIMER		35
#define		RESTRICT_UPDATE_TIMER 		36
#define		FFT_UPDATE_TIMER		37
#define		PARTICLE_ACCEL_UPDATE_TIMER	38
#define		MERGE_DENSITIES_UPDATE_TIMER	39
#define		HYDRO_PARTICLE_UPDATE_TIMER	40
#define         WORK_TIMER                      41
#define         COMMUNICATION_TIMER             42
#define		PARTICLE_IO_TIMER		43
#define		GAS_IO_TIMER			44
#define		LEVEL_TIMER			45

extern const char *timer_name[];

void init_timers();
void start_timing_level( int level );
void end_timing_level( int level );
void start_time( int timerid );
double end_time( int timerid );

double current_time( int timerid, int level );
double average_time( int timerid, int level );
double total_time( int timerid, int level );
double last_time( int timerid, int level );

#endif
