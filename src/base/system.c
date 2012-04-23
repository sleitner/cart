
/*
//  Only ported to Unix at present
*/
#ifdef _WIN32
#include <windows.h>
#else
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
#endif

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

/*
//  Gives the available memory for this process
*/
double system_get_available_memory()
{
#ifdef _WIN32
  MEMORYSTATUS r;
  ::GlobalMemoryStatus(&r);
  return (double)r.dwTotalPhys;
#else
  static struct rlimit r;
  int d = getrlimit(RLIMIT_DATA,&r);
  double v1 = (double)r.rlim_cur;
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
  unsigned long np = sysconf(_SC_PHYS_PAGES);
  unsigned long ps = sysconf(_SC_PAGESIZE);
  double v2 = (double)np*(double)ps;
  if(v1 > v2) v1 = v2;
#endif
  if(d == 0) return v1; else return (double)ULONG_MAX;
#endif
}

const char *system_getcwd()
{
  const int size = 1024;
  static char buf[1024];

  return getcwd(buf,size);
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

/*
//  Is there a better implementation?
*/
int file_exists(const char *name)
{
  FILE *f;

  f = fopen(name,"r");
  if(f != NULL)
    {
      fclose(f);
      return 1;
    }
  else
    {
      return 0;
    }
}
