#ifndef __SYSTEM_H__
#define __SYSTEM_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


/*
//  This file contains wrappers for system portability only
*/
void system_get_host_name( char *, int );
int system_get_pid();
char* system_get_time_stamp(int utc);

#endif

