#ifndef __PARALLEL_CONFIG_H__
#define __PARALLEL_CONFIG_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>


#define MPI_TASK_TYPE_UNDEFINED     0
#define MPI_TASK_TYPE_RUN           1
#define MPI_TASK_TYPE_FFT           2


struct ART_COMM_TYPE
{
  MPI_Comm world;
  MPI_Comm run;
  MPI_Comm fft;
};


struct ART_TASK_TYPE
{
  int rank, size;
};


struct ART_MPI_TYPE
{
  struct ART_COMM_TYPE comm;
  struct ART_TASK_TYPE world;
  struct ART_TASK_TYPE run;
  struct ART_TASK_TYPE fft;
  int task_type;
};


#endif
