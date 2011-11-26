#ifndef __PARALLEL_H__
#define __PARALLEL_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include "parallel_config.h"


/*
//  This is a trick to ensure communicator compliance
*/
#ifdef MPI_COMM_WORLD
#undef MPI_COMM_WORLD
#endif

#ifdef MPI_COMM_SELF
#undef MPI_COMM_SELF
#endif

#define MPI_COMM_WORLD
#define MPI_COMM_SELF


#define	MASTER_NODE			0


extern struct ART_MPI_TYPE mpi;


extern int num_procs;
extern int local_proc_id;
extern int tasks_per_node;
extern int proc_sfc_index[MAX_PROCS+1];

void configure_runtime_setup(); // Should be called by ALL tasks

void config_init_parallel();
void config_verify_parallel();

void init_parallel_grid();
int processor_owner( int sfc );

#define MPI_CUSTOM_NONE			0x0U

extern unsigned int mpi_custom_flags;

void print_comm_contents(MPI_Comm com, const char *name);

#endif
