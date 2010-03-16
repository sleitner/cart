#ifndef __LOGGING_H__
#define __LOGGING_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


void init_logging( int restart );
void finalize_logging();
void log_diagnostics();
#ifdef LOG_STAR_CREATION 
void log_star_creation( int icell, double mass );
void combine_star_creation_log();
void finalize_star_creation_log( char *filename_sclog );
void append_file(char *file_path_from, char *file_path_to);
void copy_file(char *file_path_from, char *file_path_to);
void check_restart_star_creation();
void wipe_restart_star_creation( double aexpn );
void wipe_temp();
int count_lines(char *filename);
#endif 

#ifdef DEBUG
void debug_breakpoint(int timerid, int start, const char *file, int line);
#define SET_MARKER(id)  debug_breakpoint(id,-1,__FILE__,__LINE__)
#endif

#endif
