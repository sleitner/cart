#ifndef __LOGGING_H__
#define __LOGGING_H__

void init_logging( int restart );
void finalize_logging();
void log_diagnostics();

#ifdef DEBUG
void log_in_debug(int timerid, int start, const char *file, int line);
#endif

#endif
