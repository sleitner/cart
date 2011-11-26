#include "config.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "parallel_config.h"

/*
//  This is a trick to ensure communicator comliance
*/
struct ART_MPI_TYPE mpi = { { MPI_COMM_WORLD, MPI_COMM_WORLD, MPI_COMM_NULL }, MPI_TASK_TYPE_UNDEFINED };


#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "parallel.h"
#include "root_grid_fft.h"
#include "sfc.h"
#include "tree.h"


int num_procs;
int local_proc_id;
int tasks_per_node;
int proc_sfc_index[MAX_PROCS+1];


unsigned int mpi_custom_flags = MPI_CUSTOM_NONE;


void configure_runtime_setup()
{
  const char *str;

  MPI_Group world_grp, run_grp, fft_grp;
  int i, irun, nrun, ifft, nfft, sel[MAX_PROCS];
  char **str1, **str2;
  char *buf1, *buf2;

  MPI_Comm_size(mpi.comm.world,&mpi.world.size);
  MPI_Comm_rank(mpi.comm.world,&mpi.world.rank);

  /*
  //  for cart_error to work
  */
  num_procs = mpi.world.size;
  local_proc_id = mpi.world.rank;

  MPI_Comm_group(mpi.comm.world,&world_grp);

  /*
  //  ************************************************
  //
  //  Main configuration of the runtime setup
  //
  //  ------------------------------------------------
  */
  str = extract_option1("mpi-setup","mpi",NULL);
  if(str == NULL)
    {
      /* 
      //  Default behaviour: all tasks are run tasks, at least one K-slice per fft task 
      */
      irun = 0;
      nrun = mpi.world.size;

      ifft = 0;
      nfft = min(num_grid,mpi.world.size);
    }
  else
    {
      if(sscanf(str,"run:%d-%d,fft:%d-%d",&irun,&nrun,&ifft,&nfft) != 4)
	{
	  cart_error("A valid format for the --mpi-setup option argument is run:N1-N2,fft:N3-N4, where N1-N2 is the range of run tasks ids, and N3-N4 is the range of fft tasks ids.");
	}

      if(irun<0 || irun>nrun)
	{
	  cart_error("Invalid range %d - %d",irun,nrun);
	}
      if(ifft<0 || ifft>nfft)
	{
	  cart_error("Invalid range %d - %d",ifft,nfft);
	}

      if(nrun >= mpi.world.size)
	{
	  cart_error("The range of run tasks overflows the available number of tasks %d",mpi.world.size);
	}
      if(nfft >= mpi.world.size)
	{
	  cart_error("The range of fft tasks overflows the available number of tasks %d",mpi.world.size);
	}

      nrun = nrun - irun + 1;
      nfft = nfft - ifft + 1;
    }
 
  cart_assert(irun>=0 && irun+nrun<=mpi.world.size);
  cart_assert(ifft>=0 && ifft+nfft<=mpi.world.size);

  for(i=0; i<nrun; i++) sel[i] = irun + i;
  MPI_Group_incl(world_grp,nrun,sel,&run_grp);

  for(i=0; i<nfft; i++) sel[i] = ifft + i;
  MPI_Group_incl(world_grp,nfft,sel,&fft_grp);

  /*
  //  ************************************************
  //
  //  Create our communicators, etc (no customization here)
  */
  MPI_Comm_create(mpi.comm.world,run_grp,&mpi.comm.run);
  MPI_Comm_create(mpi.comm.world,fft_grp,&mpi.comm.fft);

  mpi.task_type = 0;

  /*
  //  Sizes and ranks can only be safely querued from a group, 
  //  not a communicator!!!
  */
  MPI_Group_rank(run_grp,&i);
  if(i != MPI_UNDEFINED)
    {
      mpi.task_type += MPI_TASK_TYPE_RUN;
      mpi.run.size = nrun;
      mpi.run.rank = i;
    }
  else
    {
      mpi.run.size = 0;
      mpi.run.rank = i;
    }

  MPI_Group_rank(fft_grp,&i);
  if(i != MPI_UNDEFINED)
    {
      mpi.task_type += MPI_TASK_TYPE_FFT;
      mpi.fft.size = nrun;
      mpi.fft.rank = i;
    }
  else
    {
      mpi.fft.size = 0;
      mpi.fft.rank = i;
    }

  MPI_Group_free(&world_grp);

  str = extract_option1("num-omp-threads","omp",NULL);
  if(str != NULL)
    {
#ifdef _OPENMP
      if(sscanf(str,"%d",&i)!=1 || i<1 || i>256)
	{
	  cart_error("--num-omp-threads=<num> option requires a positive integer <num> as an argument");
	}
      omp_set_num_threads(i);
      cart_debug("num openmp threads = %u", omp_get_max_threads() );
#else
      cart_debug("OpenMP support is not compiled in; ignoring --num-omp-threads option.");
#endif
    }

  root_grid_fft_init(run_grp,fft_grp);

  MPI_Group_free(&run_grp);
  MPI_Group_free(&fft_grp);

  /*
  //  Measure tasks per node
  */
  buf1 = cart_alloc(char,mpi.world.size*MPI_MAX_PROCESSOR_NAME);
  buf2 = cart_alloc(char,mpi.world.size*MPI_MAX_PROCESSOR_NAME);
  str1 = cart_alloc(char*,mpi.world.size);
  str2 = cart_alloc(char*,mpi.world.size);
  for(i=0; i<mpi.world.size; i++)
    {
      str1[i] = buf1 + i*MPI_MAX_PROCESSOR_NAME;
      str2[i] = buf2 + i*MPI_MAX_PROCESSOR_NAME;
    }

  MPI_Get_processor_name(str1[0],&i);

  for(i=1; i<mpi.world.size; i++)
    {
      strcpy(str1[i],str1[0]);
    }

  MPI_Alltoall(buf1,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,buf2,MPI_MAX_PROCESSOR_NAME,MPI_CHAR,mpi.comm.world);

  tasks_per_node = 0;
  for(i=0; i<mpi.world.size; i++)
    {
      if(strcmp(str2[i],str2[mpi.world.rank]) == 0) tasks_per_node++;
    }

  cart_debug("Tasks per node: %d",tasks_per_node);

  cart_assert(tasks_per_node > 0);

  cart_free(buf1);
  cart_free(buf2);
  cart_free(str1);
  cart_free(str2);
}


void config_init_parallel()
{
  control_parameter_add2(control_parameter_int,&mpi_custom_flags,"@mpi:custom-flags","mpi_custom_flags","flags that can be set to customize MPI performance. This parameter is experimental and may be removed in the future.");
}


void config_verify_parallel()
{
}


/*******************************************************  
 * init_parallel_grid
 ******************************************************/
void init_parallel_grid() 
/* purpose: initializes proc_sfc_index array 
 * requires: local_proc_id and num_procs are already set
 */
{
	int i;
	int index;

	/* to start with, divide all root octs equally amongst all nodes */
	for ( i = 0, index = 0; i < num_procs; i++, index += num_root_cells / num_procs ) {
		proc_sfc_index[i] = index;
	}
	proc_sfc_index[num_procs] = num_root_cells;
}

/*******************************************************
 * processor_owner
 ******************************************************/
int processor_owner( int sfc ) 
/* purpose: determines which processor currently 
 * 	owns the given sfc index 
 */
{
	int a, b, c;

	cart_assert( sfc >= 0 && sfc < max_sfc_index );

	/* determine if sfc is local */
	if ( sfc < proc_sfc_index[local_proc_id] ) {
		a = 0;
		b = local_proc_id-1;
	} else if ( sfc >= proc_sfc_index[local_proc_id+1] ) {
		a = local_proc_id + 1;
		b = num_procs-1;
	} else {
		return local_proc_id;
	}

	/* do binary search between procs a & b */
	while ( a != b ) {
		c = ( a + b + 1) / 2;
		
		if ( sfc < proc_sfc_index[c] ) {
			b = c-1;
		} else {
			a = c;
		}
	}

	cart_assert( a >= 0 && a < num_procs );
	cart_assert( sfc >= proc_sfc_index[a] && sfc < proc_sfc_index[a+1] );

	return a;
}


/*
//  This function must be called by all members of communicator com.
*/
void print_comm_contents(MPI_Comm com, const char *name)
{
  MPI_Group world, local;
  int i, n, *ranks_local, *ranks_world;

  MPI_Comm_group(mpi.comm.world,&world);
  MPI_Comm_group(com,&local);

  MPI_Group_size(local,&n);
  MPI_Group_rank(local,&i);

  if(i == 0)
    {
      ranks_local = cart_alloc(int,n);
      ranks_world = cart_alloc(int,n);
  
      for(i=0; i<n; i++) ranks_local[i] = i;

      MPI_Group_translate_ranks(local,n,ranks_local,world,ranks_world);

      cart_debug("Communicator %s (%p), size = %d:",name,com,n);
      for(i=0; i<n; i++) cart_debug("id = %d -> world id = %d",i,ranks_world[i]);

      cart_free(ranks_local);
      cart_free(ranks_world);
    }

  MPI_Group_free(&local);
  MPI_Group_free(&world);
}
