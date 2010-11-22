#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "parallel.h"
#include "particle.h"
#include "tree.h"

#include "chull.h"
#include "halo_finder.h"

#include "gic_tools.h"


void gicMakeMaskAddPoints(int mode, int n, const int *pos, char *mask)
{
  int j;
  const int *idx;

  switch(mode)
    {
    case GIC_MASK_MODE_PLAIN:
      {
	for(j=0; j<n; j++)
	  {
	    idx = pos + nDim*j;
	    mask[idx[0]+num_grid*(idx[1]+num_grid*idx[2])] = 1;
	  }
	break;
      }
    case GIC_MASK_MODE_CHULL:
      {
	chAddPoints(n,pos);
	break;
      }
    default:
      {
	cart_error("gicMakeMask mode %d is not defined.",mode);
      }
    }
}


void gicMakeMask(const char *filename, int num_halos, const struct HALO **halos, float size, int mode)
{
  int ih;
  int i, j, k, n, nbuf;
  int *pos, min[nDim], max[nDim], shift[nDim], box[nDim], p[3];
  char *mask = 0;
  const halo *h;
  FILE *f;

  cart_assert(halos != NULL);
  cart_assert(nDim == 3);

  if(particle_species_indices[1] != num_root_cells)
    {
      cart_debug("The first particle species is not an exact cube, num_particles[1]=%d, num_root_cells=%d",particle_species_indices[1],num_root_cells);
      cart_debug("gicMakeMask call is skipped.");
      return;
    }

  if(local_proc_id == MASTER_NODE)
    {
      mask = cart_alloc(char,num_root_cells);
      memset(mask,0,sizeof(char)*num_root_cells);
    }

  /*
  //  Do a union of masks for each halo
  */
  for(ih=0; ih<num_halos; ih++)
    {
      h = halos[ih];
      cart_assert(h != NULL);
  
      for(n=j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id[j]<particle_species_indices[1])
	{
	  if(compute_distance_periodic(particle_x[j],(double*)h->pos) < h->rvir*size)
	    {
	      n++;
	    }
	}

      nbuf = n + 1;  /* in case n is zero */
      pos = cart_alloc(int,nDim*nbuf);
  
      /*
      //  Shift the halo to the box center - chull does not do peridoc BC
      */
      if(mode == GIC_MASK_MODE_CHULL)
	{
	  for(k=0; k<nDim; k++)
	    {
	      box[k] = num_grid;
	      shift[k] = box[k] + num_grid/2 - h->pos[k];
	    }
	}

      for(n=j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id[j]<particle_species_indices[1])
	{
	  /*
	  //  If particle is NOW within the size*Rvir, get it INITIAL position
	  */
	  if(compute_distance_periodic(particle_x[j],(double*)h->pos) < h->rvir*size)
	    {
	      pos[nDim*n+0] = particle_id[j] % num_grid;
#if (nDim > 1)
	      pos[nDim*n+1] = (particle_id[j]/num_grid) % num_grid;
#if (nDim > 2)
	      pos[nDim*n+2] = (particle_id[j]/(num_grid*num_grid)) % num_grid;
#endif
#endif
	      if(mode == GIC_MASK_MODE_CHULL)
		{
		  for(k=0; k<nDim; k++)
		    {
		      pos[nDim*n+k] = (pos[nDim*n+k]+shift[k]) % box[k];
		    }
		}

	      n++;
	    }
	}

      /*
      //  Since a union of several convex hulls is not a convex hull, 
      //  we need to make a hull on the master node.
      */
      if(local_proc_id == MASTER_NODE)
	{
	  gicMakeMaskAddPoints(mode,n,pos,mask);
	  for(i=1; i<num_procs; i++)
	    {
	      MPI_Recv(&n,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      if(n > nbuf)
		{
		  cart_free(pos);
		  nbuf = n;
		  pos = cart_alloc(int,nDim*nbuf);
		}
	      MPI_Recv(pos,nDim*n,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      gicMakeMaskAddPoints(mode,n,pos,mask);
	    }
	  cart_free(pos);

	  if(mode == GIC_MASK_MODE_CHULL)
	    {
	      /*
	      // Check that the whole region fits inside the box
	      */
	      chGetLimits(min,max);
	      for(k=0; k<nDim; k++)
		{
		  cart_assert(min[k]>=0 && min[k]<num_grid);
		  cart_assert(max[k]>=0 && max[k]<num_grid);
		  if(min[k]==0 || max[k]==num_grid-1)
		    {
		      cart_debug("Unable to account for the periodic BC.");
		      cart_debug("gicMakeMask call is skipped.");
		      cart_free(mask);
		      return;
		    }
		}

	      chConstruct();
	  
	      for(k=0; k<num_grid; k++)
		{
		  for(j=0; j<num_grid; j++)
		    {
		      for(i=0; i<num_grid; i++)
			{
			  p[0] = (i+shift[0]) % box[0];
			  p[1] = (j+shift[1]) % box[1];
			  p[2] = (k+shift[2]) % box[2];
			  if(chIsPointInside(p)) mask[i+num_grid*(j+num_grid*k)] = 1;
			}
		    }
		}

	      chReset();
	    }
	}
      else
	{
	  MPI_Send(&n,1,MPI_INT,MASTER_NODE,0,MPI_COMM_WORLD);
	  MPI_Send(pos,nDim*n,MPI_DOUBLE,MASTER_NODE,0,MPI_COMM_WORLD);
	  cart_free(pos);
	}

      cart_debug("Done halo #%d, id=%d",ih,h->id);
    }

  /*
  //  Write plain ascii file. Should consider adding writing a GIC .rfm file directly.
  */
  if(local_proc_id == MASTER_NODE)
    {
      f = fopen(filename,"w");
      cart_assert(f != NULL);

      for(k=0; k<num_grid; k++)
	{
	  for(j=0; j<num_grid; j++)
	    {
	      for(i=0; i<num_grid; i++)
		{
		  if(mask[i+num_grid*(j+num_grid*k)])
		    {
		      fprintf(f,"%d %d %d\n",i+1,j+1,k+1);
		    }
		}
	    }
	}
  
      fclose(f);

      cart_free(mask);
    }
}

