#include "config.h"

#include <stdio.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "parallel.h"
#include "particle.h"
#include "sfc.h"
#include "tree.h"


/*
//  Data structure
*/
typedef struct MESH_STRUCTURE
{
  int NumCells;
  int NumParticles;
  int NumBorderCells[num_neighbors];
}
MeshStructure;


/*
//  Recursive worker
*/
void extCountRootCellWorker(int top_level, int top_cell, int cell, MeshStructure mesh[])
{
  int j, l, level, top_neighbor;
  int children[num_children], neighbors[num_neighbors];
#ifdef PARTICLES
  int part, num_parts;
#endif /* PARTICLES */

  cart_assert(cell>=0 && cell<num_cells);

  level = cell_level(cell);
  cell_all_neighbors(cell,neighbors);

  mesh[level].NumCells++;

#ifdef PARTICLES
  /* count number of particles in cell */
  num_parts = 0;
  part = cell_particle_list[cell];
  while(part != NULL_PARTICLE)
    {
    num_parts++;
    part = particle_list_next[part];              
  }
  mesh[level].NumParticles += num_parts;
#endif /* PARTICLES */

  for(j=0; j<num_neighbors; j++)
    {
      for(l=cell_level(neighbors[j]), top_neighbor=neighbors[j]; l>top_level; l--) 
	{
	  top_neighbor = cell_parent_cell(top_neighbor);
	}	  

      if(top_cell != top_neighbor)
	{
	  mesh[level].NumBorderCells[j]++;
	}
    }
      
  if(cell_is_refined(cell))
    {
      cell_all_children(cell,children);
      for(j=0; j<num_children; j++) extCountRootCellWorker(top_level,top_cell,children[j],mesh);
    }
}

/*
//  Driver for recursive worker
*/
void extCountRootCell(int top_cell, MeshStructure mesh[])
{
  int j, level;
  int top_level = cell_level(top_cell);

  for(level=top_level; level<=max_level; level++)
    {
      mesh[level].NumCells = 0;
      mesh[level].NumParticles = 0;
      for(j=0; j<num_neighbors; j++) mesh[level].NumBorderCells[j] = 0;
    }

  extCountRootCellWorker(top_level,top_cell,top_cell,mesh);
}


void extDumpMeshStructure(const char *filename, int level_shift)
{
  FILE *f;
  int j, i, level, ijk[3], coords[3], shifts[3];
  int sfc, root_cell, floor_level, cell_floor_level, subsample;
  DEFINE_LEVEL_ARRAY(MeshStructure,mesh);

  cart_assert(level_shift>=0 && level_shift<=1);
  subsample = 1<<level_shift;

  floor_level = max_level_now_global(mpi.comm.run);

  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(filename,"w");
      if(f == 0)
        {
          cart_debug("Unable to open file %s for writing.",filename);
          return;
        }

      fprintf(f,"# Root dimensions, floor level:\n");
      fprintf(f,"%d %d %d %d\n",subsample*num_grid,subsample*num_grid,subsample*num_grid,floor_level-level_shift);
      fprintf(f,"# Format of each cell entry:\n");
      fprintf(f,"#   Cell <I> <J> <K> <sfc index> <floor level>\n");
      fprintf(f,"#   Then, for each level between %d and <floor level>:\n",min_level);
      fprintf(f,"#   <Number of cells at this level> <Number of particle at this level> <Number of boundary cells at this level for 6 directions (-X,+X,-Y,+Y,-Z,+Z)>\n");
    }

  for(ijk[2]=0; nDim>2 && ijk[2]<subsample*num_grid; ijk[2]++)
    {
      cart_debug("Doing Z-slice #%d...",ijk[2]);
      coords[2] = ijk[2]/subsample;
      shifts[2] = ijk[2]%subsample;
      for(ijk[1]=0; nDim>1 && ijk[1]<subsample*num_grid; ijk[1]++)
	{
	  coords[1] = ijk[1]/subsample;
	  shifts[1] = ijk[1]%subsample;
	  for(ijk[0]=0; ijk[0]<subsample*num_grid; ijk[0]++)
	    {
	      coords[0] = ijk[0]/subsample;
	      shifts[0] = ijk[0]%subsample;

	      sfc = sfc_index(coords);

	      if(root_cell_type(sfc) == CELL_TYPE_LOCAL)
		{
		  root_cell = root_cell_location(sfc);

		  switch(level_shift)
		    {
		    case 0:
		      {
			extCountRootCell(root_cell,mesh);
			break;
		      }
		    case 1:
		      {
			if(cell_is_refined(root_cell))
			  {
			    root_cell = cell_child(root_cell,shifts[0]+2*shifts[1]+4*shifts[2]);
			    extCountRootCell(root_cell,mesh);
			  }
			else
			  {
			    extCountRootCell(root_cell,mesh);
			    mesh[min_level+1].NumCells = 1;
			    mesh[min_level+1].NumParticles = mesh[min_level].NumParticles;
			    for(i=0; i<num_neighbors; i++)
			      {
				mesh[min_level+1].NumBorderCells[i] = 1;
			      }
			  }

			/*
			// Shift by one level
			*/
			for(level=min_level; level<floor_level; level++)
			  {
			    mesh[level] = mesh[level+1];
			  }

			break;
		      }
		    default:
		      {
			cart_error("Should not be here!");
		      }
		    }

		  /*
		  //  Send to master node
		  */
		  if(local_proc_id != MASTER_NODE)
		    {
		      MPI_Send(mesh+min_level,sizeof(MeshStructure)*(floor_level-min_level+1-level_shift),MPI_BYTE,MASTER_NODE,0,mpi.comm.run);
		    }
		}
	      else
		{
		  if(local_proc_id == MASTER_NODE)
		    {
		      MPI_Recv(mesh+min_level,sizeof(MeshStructure)*(floor_level-min_level+1-level_shift),MPI_BYTE,processor_owner(sfc),0,mpi.comm.run,MPI_STATUS_IGNORE);
		    }
		}
	      
	      if(local_proc_id == MASTER_NODE)
		{
		  for(cell_floor_level=level=min_level; level<=floor_level-level_shift; level++)
		    {
		      if(mesh[level].NumCells > 0) cell_floor_level = level;
		    }

		  fprintf(f,"%d %d %d %d %d\n",ijk[0],ijk[1],ijk[2],sfc,cell_floor_level);
		  for(level=min_level; level<=cell_floor_level; level++)
		    {
		      fprintf(f,"%d %d",mesh[level].NumCells,mesh[level].NumParticles);
		      for(j=0; j<num_neighbors; j++) fprintf(f," %d",mesh[level].NumBorderCells[j]); 
		      fprintf(f,"\n");
		    }
		}		  
	    }
	}
    }

  if(local_proc_id == MASTER_NODE)
    {
      fclose(f);
    }
}
