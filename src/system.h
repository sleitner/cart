#ifndef __SYSTEM_H__
#define __SYSTEM_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


/*
//  This file contains wrappers for system portability only
*/
char* system_get_host_name();
char* system_get_time_stamp(int utc);

#endif

