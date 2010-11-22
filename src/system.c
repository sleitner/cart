#include "config.h"


/*
//  Only ported to Unix at present
*/
#include <errno.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

void system_get_host_name( char *name, int len )
{
  if(gethostname(name,len) != 0)
    {
      strcpy(name,"unknown");
    }
  else
    {
      name[len-1] = 0; /* just in case, gethostname does not guarantee null-terminated string. */
    }
}

int system_get_pid() {
	return getpid();
}	

char* system_get_time_stamp(int utc)
{
  struct timeval tv;

  if(gettimeofday(&tv,NULL) != 0)
    {
      return "unknown";
    }
  else
    {
      if(utc) return asctime(gmtime(&tv.tv_sec)); else return asctime(localtime(&tv.tv_sec));
    }
}

int system_mkdir(const char *name)
{
  switch(mkdir(name,S_IRWXU | S_IRGRP | S_IXGRP))
    {
    case 0:
      {
	return 0;
      }
    default:
      {
	return (errno == EEXIST) ? 0 : errno;
      }
    }
}


int system_chdir(const char *name)
{
  switch(chdir(name))
    {
    case 0:
      {
	return 0;
      }
    default:
      {
	return errno;
      }
    }
}
