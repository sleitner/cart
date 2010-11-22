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
void extCountRootCellWorker(int root_parent, int cell, MeshStructure mesh[])
{
  int j, level;
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
      if(root_parent != cell_parent_root_cell(neighbors[j]))
	{
	  mesh[level].NumBorderCells[j]++;
	}
    }
      
  if(cell_is_refined(cell))
    {
      cell_all_children(cell,children);
      for(j=0; j<num_children; j++) extCountRootCellWorker(root_parent,children[j],mesh);
    }
}

/*
//  Driver for rRecursive worker
*/
void extCountRootCell(int root_cell, MeshStructure mesh[])
{
  int j, level;

  cart_assert(cell_is_root_cell(root_cell)); 

  for(level=min_level; level<=max_level; level++)
    {
      mesh[level].NumCells = 0;
      mesh[level].NumParticles = 0;
      for(j=0; j<num_neighbors; j++) mesh[level].NumBorderCells[j] = 0;
    }

  extCountRootCellWorker(root_cell,root_cell,mesh);
}


void extDumpMeshStructure(const char *filename, int mode)
{
  FILE *f;
  int j, level, coords[3];
  int sfc, root_cell, floor_level, cell_floor_level;
  DEFINE_LEVEL_ARRAY(MeshStructure,mesh);
  DEFINE_LEVEL_ARRAY(int,ncells);
  int nparts;

  floor_level = max_level_now_global(MPI_COMM_WORLD);

  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(filename,"w");
      if(f == 0)
        {
          cart_debug("Unable to open file %s for writing.",filename);
          return;
        }

      fprintf(f,"# Root dimensions, floor level:\n");
      fprintf(f,"%d %d %d %d\n",num_grid,num_grid,num_grid,floor_level);
      fprintf(f,"# Format of each cell entry:\n");
      switch(mode)
	{
	case 1:
	  {
	    fprintf(f,"#   Cell <I> <J> <K> <sfc index> <floor level>\n");
	    fprintf(f,"#   Then, for each level between %d and <floor level>:\n",min_level);
	    break;
	  }
	default:
	  {
	    fprintf(f,"#   Cell <I> <J> <K> <sfc index>\n");
	    fprintf(f,"#   Then, for each level between %d and %d:\n",min_level,floor_level);
	  }
	}
      fprintf(f,"#   <Number of cells at this level> <Number of particle at this level> <Number of boundary cells at this level for 6 directions (-X,+X,-Y,+Y,-Z,+Z)>\n");
    }

  /*
  //  This for internal checking
  */
  nparts = 0;
  for(level=min_level; level<=floor_level; level++) ncells[level] = 0;

  for(coords[2]=0; nDim>2 && coords[2]<num_grid; coords[2]++)
    {
      cart_debug("Doing Z-slice #%d...",coords[2]);
      for(coords[1]=0; nDim>1 && coords[1]<num_grid; coords[1]++)
	{
	  for(coords[0]=0; coords[0]<num_grid; coords[0]++)
	    {
	      sfc = sfc_index(coords);

	      if(root_cell_type(sfc) == CELL_TYPE_LOCAL)
		{
		  root_cell = root_cell_location(sfc);
		  extCountRootCell(root_cell,mesh);

		  for(level=min_level; level<=floor_level; level++)
		    {
		      nparts += mesh[level].NumParticles;
		      ncells[level] += mesh[level].NumCells;
		    }

		  /*
		  //  Send to master node
		  */
		  if(local_proc_id != MASTER_NODE)
		    {
		      MPI_Send(mesh+min_level,sizeof(MeshStructure)*(floor_level-min_level+1),MPI_BYTE,MASTER_NODE,0,MPI_COMM_WORLD);
		    }
		}
	      else
		{
		  if(local_proc_id == MASTER_NODE)
		    {
		      MPI_Recv(mesh+min_level,sizeof(MeshStructure)*(floor_level-min_level+1),MPI_BYTE,processor_owner(sfc),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		    }
		}
	      
	      if(local_proc_id == MASTER_NODE)
		{
		  for(cell_floor_level=level=min_level; level<=floor_level; level++)
		    {
		      if(mesh[level].NumCells > 0) cell_floor_level = level;
		    }

		  switch(mode)
		    {
		    case 1:
		      {
			fprintf(f,"%d %d %d %d %d\n",coords[0],coords[1],coords[2],sfc,cell_floor_level);
			for(level=min_level; level<=cell_floor_level; level++)
			  {
			    fprintf(f,"%d %d",mesh[level].NumCells,mesh[level].NumParticles);
			    for(j=0; j<num_neighbors; j++) fprintf(f," %d",mesh[level].NumBorderCells[j]); 
			    fprintf(f,"\n");
			  }
			break;
		      }
		    default:
		      {
			fprintf(f,"%d %d %d %d\n",coords[0],coords[1],coords[2],sfc);
			for(level=min_level; level<=floor_level; level++)
			  {
			    fprintf(f,"%d %d",mesh[level].NumCells,mesh[level].NumParticles);
			    for(j=0; j<num_neighbors; j++) fprintf(f," %d",mesh[level].NumBorderCells[j]); 
			    fprintf(f,"\n");
			  }
		      }
		    }
		}		  
	    }
	}
    }

  if(local_proc_id == MASTER_NODE)
    {
      fclose(f);
    }

  /*
  //  Internal check
  */
#ifdef PARTICLES 
  cart_assert(nparts == num_local_particles);
#endif /* PARTICLES */

  for(level=min_level; level<=floor_level; level++)
    {
      cart_assert(ncells[level] == num_cells_per_level[level]);
    }
}
