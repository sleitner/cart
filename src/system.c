#include "config.h"


/*
//  Only ported to Unix at present
*/
#include <time.h>
#include <unistd.h>
#include <sys/time.h>


char* system_get_host_name()
{
  const int len = 1024;
  char name[len];

  if(gethostname(name,len) != 0)
    {
      strcpy(name,"unknown");
    }
  else
    {
      name[len-1] = 0; /* just in case, gethostname does not guarantee null-terminated string. */
    }

  return name;
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
