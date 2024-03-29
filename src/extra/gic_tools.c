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
#include "halos.h"

#include "gic_tools.h"


struct gicCHullConfig gictol = { 0.01, 0.05, 3, 2 };


#ifdef PARTICLES

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


void gicMakeMask(const char *filename, const halo_list *halos, float size, int mode, int level, int width)
{
  int ih;
  int i, j, k, i1, j1, k1, ioff, joff, koff, done, n, nbuf, lev, width2;
  int *pos, min[nDim], max[nDim], shift[nDim], box[nDim], p[3];
  char *mask = NULL;
  char *halo_mask = NULL;
  halo *h;
  FILE *f;
  char str[256];

  cart_assert(halos != NULL);
  cart_assert(nDim == 3);
  cart_assert(level > 0);
  cart_assert(width >= 0);

  if(particle_species_indices[1] != num_root_cells)
    {
      cart_debug("The first particle species is not an exact cube, num_particles[1]=%d, num_root_cells=%d",particle_species_indices[1],num_root_cells);
      cart_debug("gicMakeMask call is skipped.");
      return;
    }

  if(local_proc_id == MASTER_NODE)
    {
      halo_mask = cart_alloc(char,num_root_cells);
      mask = cart_alloc(char,num_root_cells);
      memset(mask,0,sizeof(char)*num_root_cells);
    }

  width2 = width*width;

  /*
  //  Do a union of masks for each halo
  */
  for(ih=0; ih<halos->num_halos; ih++)
    {
      h = &halos->list[ih];
      cart_assert(h != NULL);
      if ( !h->flag ) continue;
  
      for(n=j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id[j]<particle_species_indices[1])
	{
	  if(compute_distance_periodic(particle_x[j],h->pos) < h->rvir*size)
	    {
	      n++;
	    }
	}

      cart_debug("Found %d particles.",n);

      nbuf = n + 1;  /* in case n is zero */
      pos = cart_alloc(int,nDim*nbuf);
  
      /*
      //  Shift the halo to the box center - chull does not do periodic BC
      */
      for(k=0; k<nDim; k++)
	{
	  box[k] = num_grid;
	  shift[k] = num_grid/2 - h->pos[k];
	  shift[k] = (box[k]+shift[k]) % box[k];
	}

      for(n=j=0; j<num_particles; j++) if(particle_id[j]!=NULL_PARTICLE && particle_id[j]<particle_species_indices[1])
	{
	  /*
	  //  If particle is NOW within the size*Rvir, get it INITIAL position
	  */
	  if(compute_distance_periodic(particle_x[j],h->pos) < h->rvir*size)
	    {
	      pos[nDim*n+0] = particle_id[j] % num_grid;
#if (nDim > 1)
	      pos[nDim*n+1] = (particle_id[j]/num_grid) % num_grid;
#if (nDim > 2)
	      pos[nDim*n+2] = (particle_id[j]/(num_grid*num_grid)) % num_grid;
#endif
#endif
	      for(k=0; k<nDim; k++)
		{
		  pos[nDim*n+k] = (pos[nDim*n+k]+shift[k]) % box[k];
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
	  memset(halo_mask,0,sizeof(char)*num_root_cells);
	  gicMakeMaskAddPoints(mode,n,pos,halo_mask);
	  for(i=1; i<num_procs; i++)
	    {
	      MPI_Recv(&n,1,MPI_INT,i,0,mpi.comm.run,MPI_STATUS_IGNORE);
	      if(n > nbuf)
		{
		  cart_free(pos);
		  nbuf = n;
		  pos = cart_alloc(int,nDim*nbuf);
		}
	      MPI_Recv(pos,nDim*n,MPI_INT,i,0,mpi.comm.run,MPI_STATUS_IGNORE);
	      gicMakeMaskAddPoints(mode,n,pos,halo_mask);
	    }
	  cart_free(pos);

	  switch(mode)
	    {
	    case GIC_MASK_MODE_PLAIN:
	      {
		for(k=0; k<num_grid; k++)
		  {
		    for(j=0; j<num_grid; j++)
		      {
			for(i=0; i<num_grid; i++)
			  {
			    halo_mask[i+num_grid*(j+num_grid*k)] *= level;
			  }
		      }
		  }
		break;
	      }
	    case GIC_MASK_MODE_CHULL:
	      {
		/*
		// Check that the whole region fits inside the box
		*/
		chGetLimits(min,max);
		for(k=0; k<nDim; k++)
		  {
		    cart_debug("Region[%1d]: %d %d",k,min[k],max[k]);
		    cart_assert(min[k]>=0 && min[k]<num_grid);
		    cart_assert(max[k]>=0 && max[k]<num_grid);
		    if(min[k]==0 || max[k]==num_grid-1)
		      {
			cart_debug("Unable to account for the periodic BC.");
			cart_debug("gicMakeMask call is skipped.");
			cart_free(halo_mask);
			cart_free(mask);
			return;
		      }
		  }

		chMakeHull(gictol.Particles,gictol.Volume,gictol.NumSteps,gictol.Loud);

		for(k=0; k<num_grid; k++)
		  {
		    for(j=0; j<num_grid; j++)
		      {
			for(i=0; i<num_grid; i++)
			  {
			    p[0] = i;
			    p[1] = j;
			    p[2] = k;
			    if(chIsPointInside(p)) halo_mask[i+num_grid*(j+num_grid*k)] = level;
			  }
		      }
		  }

		chReset();
	      }
	    }

	  /*
	  // Merge halo mask into complete mask
	  */
	  for(k=0; k<num_grid; k++)
	  {
              k1 = (k+shift[2])%box[2];
	      for(j=0; j<num_grid; j++)
	      {
		  j1 = (j+shift[1])%box[1];
		  for(i=0; i<num_grid; i++)
		  {
		      i1 = (i+shift[0])%box[0];
		      mask[i+num_grid*(j+num_grid*k)] |= halo_mask[i1+num_grid*(j1+num_grid*k1)];
		  }
	      }
	  }
	}
      else
	{
	  MPI_Send(&n,1,MPI_INT,MASTER_NODE,0,mpi.comm.run);
	  MPI_Send(pos,nDim*n,MPI_INT,MASTER_NODE,0,mpi.comm.run);
	  cart_free(pos);
	}

      cart_debug("Done halo #%d, id=%d",ih,h->id);
    }


  if ( local_proc_id == MASTER_NODE ) {
	  /*
	  //  Diffuse to create lower level buffers
	  */
	  for(lev=level-1; lev>0; lev--)
	  {
		  for(k=0; k<num_grid; k++)
		  {
			  for(j=0; j<num_grid; j++)
			  {
				  for(i=0; i<num_grid; i++)
				  {
					  if(mask[i+num_grid*(j+num_grid*k)] == 0)
					  {
						  done = 0;
						  for(koff=-width; done==0 && koff<=width; koff++) for(joff=-width; done==0 && joff<=width; joff++) for(ioff=-width; done==0 && ioff<=width; ioff++)
						  {
							  if(ioff*ioff+joff*joff+koff*koff <= width2)
							  {
								  i1 = (i+ioff+num_grid) % num_grid;
								  j1 = (j+joff+num_grid) % num_grid;
								  k1 = (k+koff+num_grid) % num_grid;
								  if(mask[i1+num_grid*(j1+num_grid*k1)] > lev)
								  {
									  mask[i+num_grid*(j+num_grid*k)] = lev;
									  done = 1;
								  }
							  }
						  }
					  }
				  }
			  }
		  }
	  }

	  /*
	  //  Write plain ascii files. 
	  */
	  for(lev=1; lev<=level; lev++)
	  {
		  sprintf(str,"%s.%-u",filename,lev);
		  f = fopen(str,"w");
		  cart_assert(f != NULL);

		  /* fprintf(f,"%d %d %d\n",shift[0],shift[1],shift[2]); */
		  /* DHR: masks are now always unshifted */
		  fprintf(f, "0 0 0\n");

		  for(k=0; k<num_grid; k++)
		  {
			  for(j=0; j<num_grid; j++)
			  {
				  for(i=0; i<num_grid; i++)
				  {
					  if(mask[i+num_grid*(j+num_grid*k)] == lev)
					  {
						  fprintf(f,"%d %d %d\n",i+1,j+1,k+1);
					  }
				  }
			  }
		  }

		  fclose(f);
	  }

	  cart_free(halo_mask);
	  cart_free(mask);
  }
}

#endif /* PARTICLES */
