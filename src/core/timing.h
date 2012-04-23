#ifndef __TIMING_H__
#define __TIMING_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


typedef struct TIMER {
	int num_calls;
	double current_time;
	double last_time;
	double total_time;
#ifdef DEBUG_TIMING
	double last_wtime;
	const char* last_file;
	int last_line;
#endif
} timer;

#define	    TOTAL_TIME                              0
#define	    INIT_TIMER                              1
#define	    RESTART_TIMER                           2
#define     OUTPUT_TIMER                            3
#define	    IO_TIMER                                4
#define     PARTICLE_WRITE_IO_TIMER                 5
#define     PARTICLE_READ_IO_TIMER                  6
#define     GAS_WRITE_IO_TIMER                      7
#define     GAS_READ_IO_TIMER                       8
#define	    LOAD_BALANCE_TIMER                      9
#define	    LOAD_BALANCE_COMMUNICATION_TIMER        10
#define     CHOOSE_TIMESTEP_TIMER                   11
#define     CHOOSE_TIMESTEP_COMMUNICATION_TIMER     12
#define	    HYDRO_TIMER                             13
#define	    HYDRO_UPDATE_TIMER                      14
#define     COOLING_TIMER                           15
#define     ACCEL_PARTS_TIMER                       16
#define	    MOVE_PARTS_TIMER                        17
#define     STELLAR_FEEDBACK_UPDATE_TIMER           18
#define	    UPDATE_PARTS_TIMER                      19
#define	    UPDATE_PARTS_COMMUNICATION_TIMER        20
#define	    TRADE_PARTICLE_TIMER                    21
#define	    TRADE_PARTICLE_COMMUNICATION_TIMER      22
#define     GRAVITY_TIMER                           23
#define	    FFT_TIMER                               24
#define	    FFT_COMMUNICATION_TIMER                 25
#define     FFT_UPDATE_TIMER                        26
#define     PROLONGATE_UPDATE_TIMER                 27
#define	    RESTRICT_UPDATE_TIMER                   28
#define     SMOOTH_TIMER                            29
#define     SMOOTH_SETUP_TIMER                      30
#define	    SMOOTH_UPDATE_TIMER	                    31
#define     SMOOTH_COMMUNICATION_TIMER              32
#define     HYDRO_ACCEL_TIMER                       33
#define     HYDRO_ACCEL_UPDATE_TIMER                34
#define     PARTICLE_ACCEL_TIMER                    35
#define	    PARTICLE_ACCEL_UPDATE_TIMER	            36
#define     DENSITY_TIMER                           37
#define     MERGE_DENSITY_TIMER                     38
#define     MERGE_DENSITIES_UPDATE_TIMER            39
#define     MERGE_DENSITIES_COMMUNICATION_TIMER     40
#define     REFINEMENT_TIMER                        41
#define     DIFFUSION_STEP_TIMER                    42
#define     DIFFUSION_UPDATE_TIMER                  43
#define     MODIFY_UPDATE_TIMER                     44
#define     DEREFINE_UPDATE_TIMER                   45
#define     SPLIT_BUFFER_TIMER                      46
#define     JOIN_BUFFER_TIMER                       47
#define     SPLIT_BUFFER_COMMUNICATION_TIMER        48
#define     JOIN_BUFFER_COMMUNICATION_TIMER         49
#define     UPDATE_TIMER                            50
#define	    UPDATE_SEND_TIMER                       51
#define	    UPDATE_RECV_TIMER                       52
#define	    BUILD_CELL_BUFFER_TIMER                 53
#define	    SELECT_LEVEL_TIMER                      54
#define	    MAX_LEVEL_TIMER                         55
#define     SMOOTH_PARTICLE_DENSITY_TIMER           56
#define     HALO_FINDER_TIMER                       57
#define     HALO_FINDER_MASS_TIMER                  58
#define     HALO_FINDER_RECENTER_TIMER              59
#define		HALO_FINDER_WRITE_PARTICLES_TIMER		60
#define	    COMMUNICATION_TIMER                     61
#define	    LOWER_LEVEL_TIMER                       62
#define	    WORK_TIMER                              63

#ifdef RADIATIVE_TRANSFER

#define	    RT_TABLES_TIMER                         (WORK_TIMER+1)
#define	    RT_COOLING_TIMER                        (WORK_TIMER+2)
#define	    RT_LEVEL_UPDATE_TIMER                   (WORK_TIMER+3)
#define	    RT_GLOBAL_UPDATE_TIMER                  (WORK_TIMER+4)
#define	    RT_AFTER_DENSITY_TIMER                  (WORK_TIMER+5)
#define	    RT_TREE_EMULATOR_UPDATE_TIMER           (WORK_TIMER+6)
#define	    RT_SINGLE_SOURCE_UPDATE_TIMER           (WORK_TIMER+7)
#define	    RT_SOLVE_EQUATION_UPDATE_TIMER          (WORK_TIMER+8)
#define	    LEVEL_TIMER                             (WORK_TIMER+9)

#else

#define	    LEVEL_TIMER                             (WORK_TIMER+1)

#endif

#define		NUM_TIMERS                              (LEVEL_TIMER+1)


#define start_time(timerid) start_time_at_location(timerid,__FILE__,__LINE__)
#define end_time(timerid)   end_time_at_location(timerid,__FILE__,__LINE__)

#ifdef MPE_LOG
extern int event[2*NUM_TIMERS];
#endif

extern const char *timer_name[][2];

void init_timers();
void start_timing_level( int level );
void end_timing_level( int level );
void start_time_at_location( int timerid, const char *file, int line );
double end_time_at_location( int timerid, const char *file, int line );

double current_time( int timerid, int level );
double average_time( int timerid, int level );
double total_time( int timerid, int level );
double last_time( int timerid, int level );

#ifdef DEBUG
void debug_breakpoint(int timerid, int start, const char *file, int line);
#define SET_MARKER(id)  debug_breakpoint(id,-1,__FILE__,__LINE__)
#endif

#endif
