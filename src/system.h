#ifndef __SYSTEM_H__
#define __SYSTEM_H__


/*
//  This file contains wrappers for system portability only
*/
void system_get_host_name( char *, int );
int system_get_pid();
char* system_get_time_stamp(int utc);
double system_get_available_memory();

int system_mkdir(const char *name);
int system_chdir(const char *name);

int file_exists(const char *name);

#endif

