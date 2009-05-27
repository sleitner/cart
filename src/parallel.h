#ifndef __PARALLEL_H__
#define __PARALLEL_H__

#include <mpi.h>

#define		MAX_PROCS		512

#define		MASTER_NODE		0

extern int num_procs;
extern int local_proc_id;
extern int proc_sfc_index[MAX_PROCS+1];

void init_parallel_grid();
int processor_owner( int sfc );

#define MPI_CUSTOM_NONE      0x0U
#define MPI_CUSTOM_SYNC      0x1U

extern unsigned int mpi_customization_mode;

#endif
