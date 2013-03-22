#include "config.h"

#include <mpi.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "auxiliary.h"
#include "control_parameter.h"


extern int num_procs;
extern int local_proc_id;

extern const char* executable_name;


int drive();
void constants_init();


int main ( int argc, char *argv[] ) {
	int i, ret;
	const char **buffer;

	MPI_Init( &argc, &argv );

	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&local_proc_id);

	/*
	//  Command-line options need to be set in the very beginning, as
	//  configure_runtime_setup() uses some of them (-mpi, -omp).
	*/
	num_options = argc - 1;
	/*
	//  We copy, rather than assign, argv pointers so that
	//  we could manipulated them later without messing up the
	//  original argv[] array.
	*/
	buffer = options = cart_alloc(const char*,argc);
	for(i=0; i<num_options; i++) options[i] = argv[i+1];

	executable_name = argv[0];

	/*
	//  Initialize various physical constants as they can be used anywhere
	*/
	constants_init();

	/* 
	//  Run the actual driver
	*/
	ret = drive();

	/* 
	//  Clean-up and exit
	*/
	cart_free(buffer);

	MPI_Finalize();

	return ret;
}
