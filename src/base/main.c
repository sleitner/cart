#include "config.h"

#include <mpi.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "auxiliary.h"
#include "control_parameter.h"


int num_procs;
int local_proc_id;

const char* logfile_directory = NULL;
const char* executable_name;


int drive();


int main ( int argc, char *argv[] ) {
	int i, ret;

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
	options = cart_alloc(const char*,argc);
	for(i=0; i<num_options; i++) options[i] = argv[i+1];

	executable_name = argv[0];

	/* Run the code */
	ret = drive();


	cart_free(options);

	MPI_Finalize();

	return ret;
}