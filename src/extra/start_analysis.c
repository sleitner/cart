#include "defs.h"

#include <mpi.h>

#include "parallel.h"


void analysis();


void run_output()
{
  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Entering analysis mode...");
    }

  analysis();

  MPI_Finalize();
  exit(0);
}


void init_run()
{
  if(local_proc_id == MASTER_NODE)
    {
      cart_error("This is the analysis utility. It must be run in a restart mode.");
    }
}

