#include "config.h"

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "healpix.h"
#include "los.h"
#include "parallel.h"
#include "tree.h"


/*
//  Traverse a segment of a LOS located on a single processor.
*/
void losTraverseSegment(int id, double pos0[3], double theta, double phi, double len, int floor_level, losWorkerCallback worker, losSegment *segment)
{
  int j, cell, level, ret = 0;
  double e[3], posBox[nDim], posCell[nDim];
  double x, h, r = 0.0;
  int start = 1;

  cart_assert(worker != 0);
  cart_assert(segment != NULL);

  segment->Id = id;
  segment->Range[0] = segment->Range[1] = 0.0;

  e[0] = sin(theta)*cos(phi);
  e[1] = sin(theta)*sin(phi);
  e[2] = cos(theta);

  while(r <= len)
    {
      for(j=0; j<nDim; j++)
	{
	  posBox[j] = pos0[j] + r*e[j];
	  posBox[j] = posBox[j] - num_grid*floor(posBox[j]/num_grid);
	}

      cell = cell_find_position_above_level(floor_level,posBox);
      if(cell >= 0)
	{
	  cell_center_position(cell,posCell);
	  level = cell_level(cell);
	}
      else
	{
	  for(j=0; j<nDim; j++) posCell[j] = 0.5 + floor(posBox[j]);
	  level = min_level;
	}

      h = 1.75; /* just a bit above sqrt(3) */
      for(j=0; j<nDim; j++)
	{
	  x = (posBox[j]-posCell[j])/cell_size[level];

	  if(e[j] > 0.0) 
	    x = (0.5-x)/e[j];
	  else if(e[j] < 0.0) 
	    x = (-0.5-x)/e[j];
	  else
	    x = 10.0;

	  if(h > x) h = x;
	}
      h += 0.01;
      h *= cell_size[level];

      if(cell>=0 && cell_is_local(cell))
	{
	  if(start)
	    {
	      segment->Range[0] = r;
	      start = 0;
	    }
	  ret = worker(id,cell,r,r+h,segment->Buffer);
	  r += h;
	  segment->Range[1] = r;
	  if(ret != 0) break;
	}
      else r += h;
    }

  segment->WorkerReturnCode = ret;
}


/*
//  Collect all LOS segments from different processors on the master node
//  and broadcast them back.
*/
void losCollectSegments(const losBuffer *result, const losSegment *segment, losCollectorCallback collector)
{
  int i;
  losSegment *line;

  if(local_proc_id == MASTER_NODE)
    {
      cart_assert(collector != 0);

      line = cart_alloc(losSegment,num_procs);

      line[0] = *segment;
      for(i=1; i<num_procs; i++)
	{
	  MPI_Recv(line+i,sizeof(losSegment),MPI_BYTE,i,0,mpi.comm.run,MPI_STATUS_IGNORE);

	  if(line[i].Buffer.Size > 0)
	    {
	      line[i].Buffer.Data = cart_alloc(char,line[i].Buffer.Size);
	      MPI_Recv(line[i].Buffer.Data,line[i].Buffer.Size,MPI_BYTE,i,0,mpi.comm.run,MPI_STATUS_IGNORE);
	    }
	}
      
      collector(result,num_procs,line);

      for(i=1; i<num_procs; i++) if(line[i].Buffer.Size > 0)
	{
	  cart_free(line[i].Buffer.Data);
	}
      cart_free(line);

      if(result->Size > 0)
	{
	  MPI_Bcast(result->Data,result->Size,MPI_BYTE,MASTER_NODE,mpi.comm.run);
	}
    }
  else
    {
      MPI_Send((void *)segment,sizeof(losSegment),MPI_BYTE,MASTER_NODE,0,mpi.comm.run);
      if(segment->Buffer.Size > 0)
	{
	  MPI_Send(segment->Buffer.Data,segment->Buffer.Size,MPI_BYTE,MASTER_NODE,0,mpi.comm.run);
	}

      if(result->Size > 0)
	{
	  MPI_Bcast(result->Data,result->Size,MPI_BYTE,MASTER_NODE,mpi.comm.run);
	}
    }
}


/*
//  Traverse LOS over the whole sky, sampled by HealPIX
*/
void losTraverseSky(int nside, double pos0[3], double len, int floor_level, losBuffer *lines, losWorkerCallback worker, losCollectorCallback collector)
{
  int npix = 12*nside*nside;
  int ipix;
  double *theta, *phi;
  losSegment *segments;

  cart_assert(lines != NULL);
  cart_assert(worker != NULL);
  cart_assert(collector != NULL);

  theta = cart_alloc(double,npix);
  phi = cart_alloc(double,npix);

  /*
  //  HealPIX is not thread-safe, so we fill in angle arrays serially.
  */
  for(ipix=0; ipix<npix; ipix++)
    {
      hp_pix2ang_nest(nside,ipix,theta+ipix,phi+ipix);
    }

  /*
  //  Create worker segments. We only need to assign buffers, other data are set
  //  inside losTraverseSegment
  */
  segments = cart_alloc(losSegment,npix);
  for(ipix=0; ipix<npix; ipix++)
    {
      segments[ipix].Buffer.Size = lines[ipix].Size;
      if(lines[ipix].Size > 0)
	{
	  segments[ipix].Buffer.Data = cart_alloc(char,lines[ipix].Size);
	  memcpy(segments[ipix].Buffer.Data,lines[ipix].Data,lines[ipix].Size);
	}
    }

  /*
  //  Loop over directions
  */
#pragma omp parallel for default(none), private(ipix), shared(npix,theta,phi,segments,pos0,len,floor_level,worker), schedule(dynamic,1)
  for(ipix=0; ipix<npix; ipix++)
    {
      losTraverseSegment(ipix,pos0,theta[ipix],phi[ipix],len,floor_level,worker,segments+ipix);
    }

  /*
  //  Collect segments serially
  */
  for(ipix=0; ipix<npix; ipix++)
    {
      losCollectSegments(lines+ipix,segments+ipix,collector);
    }

  /*
  // Clean-up
  */
  for(ipix=0; ipix<npix; ipix++)
    {
      if(lines[ipix].Size > 0) cart_free(segments[ipix].Buffer.Data);
    }
  cart_free(segments);
  cart_free(theta);
  cart_free(phi);
}


