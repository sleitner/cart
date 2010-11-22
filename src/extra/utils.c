#include "config.h"

#include <math.h>
#include <stdio.h>

#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "rt_utilities.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"



void extFindMaxVar(int var, float *val, double *pos, double dist)
{
  MESH_RUN_DECLARE(level,cell);
  float vMax = -1.0e35;
  int select, cellMax = 0;
  double p[3];
  struct
  {
    float val;
    int rank;
  } in, out;

  cart_assert(var>=0 && var<num_vars);

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
      if(pos!=NULL && dist>0.0)
	{
	  cell_center_position(cell,p);
	  select = (compute_distance_periodic(p,pos) < dist);
	}
      else select = 1;

      if(select && vMax<cell_var(cell,var))
	{
	  vMax = cell_var(cell,var);
	  cellMax = cell;
	}
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
  
  in.val = vMax;
  MPI_Comm_rank(MPI_COMM_WORLD,&in.rank);
  MPI_Allreduce(&in,&out,1,MPI_FLOAT_INT,MPI_MAXLOC,MPI_COMM_WORLD);

  *val = out.val;
  if(pos != NULL)
    {
      if(local_proc_id == out.rank)
	{
	  cell_center_position(cellMax,pos);
	}
      MPI_Bcast(pos,3,MPI_DOUBLE,out.rank,MPI_COMM_WORLD);
    }
}
