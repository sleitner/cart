#ifndef __LOGGING_H__
#define __LOGGING_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


void init_logging( int restart );
void finalize_logging();
void log_diagnostics();

#ifdef DEBUG
void log_in_debug(int timerid, int start, const char *file, int line);
#define SET_MARKER(id)  log_in_debug(id,-1,__FILE__,__LINE__)
#endif

#endif
