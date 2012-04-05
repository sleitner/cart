
#include <stdio.h>

#include <mpi.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "system.h"

#define MAX_PROCS (1024*256)


int main(int argc, char *argv[])
{
  int proc, num_procs, local_proc_id;
  char task_hostnames[MAX_PROCS][256]; 
  int local_pid, pid[MAX_PROCS];
  char hostname[256];

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD,&local_proc_id);

#ifdef _OPENMP
  printf("task=%d, num_omp_threads=%u\n",local_proc_id,omp_get_max_threads());
#endif

  system_get_host_name(hostname,256);
  local_pid = system_get_pid(); 

  MPI_Gather(hostname,256,MPI_CHAR,task_hostnames,256,MPI_CHAR,0,MPI_COMM_WORLD);
  MPI_Gather(&local_pid,1,MPI_INT,pid,1,MPI_INT,0,MPI_COMM_WORLD);

  if(local_proc_id != 0) return;

  for(proc=0; proc<num_procs; proc++)
    {
      printf("mpi task %3u hostname %s:%d\n",proc,task_hostnames[proc],pid[proc]);
    }
	
  MPI_Finalize();

  return 0;
}
