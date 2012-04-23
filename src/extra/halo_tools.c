#include "config.h"

#ifdef COSMOLOGY

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cooling.h"
#include "cosmology.h"
#include "iterators.h"
#include "hydro.h"
#include "io_cart.h"
#include "parallel.h"
#include "particle.h"
#include "sfc.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "halos.h"
#include "halo_tools.h"

int halo_level( const halo *h, MPI_Comm local_comm )
{
  int llev, glev, cell;

  cart_assert(h);

  cell = cell_find_position((double *)h->pos);

  if(cell < 0)
    {
      llev = min_level - 1;
    }
  else
    {
      llev = cell_level(cell);
    }

  MPI_Allreduce(&llev,&glev,1,MPI_INT,MPI_MAX,local_comm);

  return glev;
}


void dump_region_around_halo(const char *filename, const halo *h, float size)
{
  int i, j, n, nbuf;
  int *ids;
  FILE *f;

  cart_assert(h != NULL);

  for(n=j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id[j]<particle_species_indices[1])
    {
      if(compute_distance_periodic(particle_x[j],(double *)h->pos) < h->rvir*size)
        {
          n++;
        }
    }

  nbuf = n + 1;  /* in case n is zero */
  ids = cart_alloc(int,nbuf);
  
  for(n=j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id[j]<particle_species_indices[1])
    {
      /*
      //  If particle is within the size*Rvir, save its id
      */
      if(compute_distance_periodic(particle_x[j],(double *)h->pos) < h->rvir*size)
        {
          ids[n++] = particle_id[j];
        }
    }

  /*
  //  Write from the master node
  */
  if(local_proc_id == MASTER_NODE)
    {

      cart_debug("Dumping region for halo %d",h->id);

      f = fopen(filename,"w");
      cart_assert(f != NULL);

      for(j=0; j<n; j++)
	{
	  fprintf(f,"%d\n",ids[j]);
	}

      for(i=1; i<num_procs; i++)
        {
          MPI_Recv(&n,1,MPI_INT,i,0,mpi.comm.run,MPI_STATUS_IGNORE);
          if(n > nbuf)
            {
              cart_free(ids);
	      nbuf = n;
              ids = cart_alloc(int,nbuf);
            }
          MPI_Recv(ids,n,MPI_INT,i,0,mpi.comm.run,MPI_STATUS_IGNORE);
	  for(j=0; j<n; j++)
	    {
	      fprintf(f,"%d\n",ids[j]);
	    }
        }

      fclose(f);
    }
  else
    {
      MPI_Send(&n,1,MPI_INT,MASTER_NODE,0,mpi.comm.run);
      MPI_Send(ids,n,MPI_INT,MASTER_NODE,0,mpi.comm.run);
   }

  cart_free(ids);
}


/*
//  Set halos->map with the halo id for each halo, or 0 if belongs to 
//  none; a cell belongs to a halo if it is inside its size_factor*Rtrunc, 
//  and satellites are accounted for properly.
//
//  MG: needs to be rewritten, very inefficient.
*/
void map_halos(int resolution_level, halo_list *halos, float size_factor)
{
  int j, ih, iold, *halo_levels, *map;
  MESH_RUN_DECLARE(level,cell);
  double pos[3], dx, r2, r2old, r2Cut;

  cart_assert(halos != NULL);
  if(halos->map != NULL) return; /* Already mapped */

  map = cart_alloc(int,num_cells);
  /*
  //  Zero map array
  */
  memset(map,0,num_cells*sizeof(int));

  halo_levels = cart_alloc(int,halos->num_halos);
  for(ih=0; ih<halos->num_halos; ih++) halo_levels[ih] = halo_level(&halos->list[ih],mpi.comm.run);

  /*
  //  Loop over levels first to avoid selecting cells multiple times
  */
  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  cart_debug("Mapping level %d...",level);

  for(ih=0; ih<halos->num_halos; ih++) if(halo_levels[ih] >= resolution_level)
    {
      r2Cut = pow(size_factor*halos->list[ih].rhalo,2.0);

      /*
      //  Map halo indices(+1), not ids, first (to simplify inter-comparison)
      */
#pragma omp parallel for default(none), private(_Index,cell,j,dx,r2,iold,r2old,pos), shared(_Num_level_cells,_Level_cells,level,cell_child_oct,cell_vars,r2Cut,halos,ih,map)
      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  
      if(cell_is_leaf(cell))
	{
	  cell_center_position(cell,pos);
	  for(j=0, r2=0.0; j<nDim; j++)
	    {
	      dx = pos[j] - halos->list[ih].pos[j];
	      if(dx < -0.5*num_grid) dx += num_grid;
	      if(dx >  0.5*num_grid) dx -= num_grid;
	      r2 += dx*dx;
	    }
	  
	  if(r2 < r2Cut)
	    {
	      iold = map[cell];
	      if(iold == 0)
		{
		   map[cell] = ih + 1;
		}
	      else
		{
		  cart_assert(iold>=1 && iold<=halos->num_halos);
		  iold--;
		  for(j=0, r2old=0.0; j<nDim; j++)
		    {
		      dx = pos[j] - halos->list[iold].pos[j];
		      if(dx < -0.5*num_grid) dx += num_grid;
		      if(dx >  0.5*num_grid) dx -= num_grid;
		      r2old += dx*dx;
		    }
		  if(r2old > r2)
		    {
		      /*
		      //  This cells belongs to a satellite
		      */
		      map[cell] = ih + 1;
		    }
		}
	    }
	}

      MESH_RUN_OVER_CELLS_OF_LEVEL_END;

    }

  MESH_RUN_OVER_LEVELS_END;

  cart_free(halo_levels);

  halos->map = map;
}

#endif /* COSMOLOGY */
