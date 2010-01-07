#include "config.h"

#include <mpi.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "parallel.h"
#include "rt_solver.h"
#include "rt_utilities.h"
#include "tree.h"


/*
// **************************************************************
//
// General routines - used in RT block, but not requiring RT 
// include files (can be taken outside ifdef RADIATIVE_TRANSFER).
//
// **************************************************************
*/


/*                          6  7  8  9 10 11 12 13 14 15 16 17  */
const int StencilDir1[] = { 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 2, 3 };
const int StencilDir2[] = { 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5 };

double rtuStencilDist2[rtuStencilSize];
double rtuStencilDelPos[rtuStencilSize][nDim];
double rtuStencilTensor[rtuStencilSize][nDim*(nDim+1)/2];


void rtuInitRun()
{
  int i, j, k, l;
  double r2;

  /*
  //  Fill in Stencil positions
  */
  for(l=0; l<rtuStencilSize; l++)
    {
      for(j=0; j<nDim; j++) rtuStencilDelPos[l][j] = 0.0;
    }
     
  for(l=0; l<num_neighbors; l++) rtuStencilDelPos[l][l/2] = 2*(l%2) - 1;
  for(l=0; l<rtuStencilSize-num_neighbors; l++)
    { 
      rtuStencilDelPos[num_neighbors+l][StencilDir1[l]/2] = 2*(StencilDir1[l]%2) - 1;
      rtuStencilDelPos[num_neighbors+l][StencilDir2[l]/2] = 2*(StencilDir2[l]%2) - 1;
    }

  for(l=0; l<rtuStencilSize; l++)
    {
      r2 = 0.0;
      for(j=0; j<nDim; j++) r2 += rtuStencilDelPos[l][j]*rtuStencilDelPos[l][j];
      rtuStencilDist2[l] = r2;

      for(k=j=0; j<nDim; j++)
	{
	  for(i=0; i<=j; i++)
	    {
	      rtuStencilTensor[l][k++] = rtuStencilDelPos[l][i]*rtuStencilDelPos[l][j]/r2;
	    }
	}
    }
}


/*
//  Compute a global average of a buffer
*/
void rtuGlobalAverage(int n, double *lBuffer)
{
  double gBuffer[n];
  MPI_Allreduce(lBuffer,gBuffer,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  memcpy(lBuffer,gBuffer,n*sizeof(double));
}


/*
// This routine returns 18 neighbors as if the whole mesh was uniform.
// In particular, one neighbor of higher level will appear more than
// once. The first 6 neighbors are exactly as returned by cell_all_neighbors.
// The other 12 neighbors are packed as follows:
//  6:  --0
//  7:  +-0
//  8:  -+0
//  9:  ++0
// 10:  -0-
// 11:  +0-
// 12:  -0+
// 13:  +0+
// 14:  0--
// 15:  0+-
// 16:  0-+
// 17:  0++
*/
void rtuGetStencil(int level, int cell, int nb[])
{
  int j, levNb[num_neighbors];
  
  /*
  //  First level neighbors
  */
  cell_all_neighbors(cell,nb);
  for(j=0; j<num_neighbors; j++) levNb[j] = cell_level(nb[j]);
  
  /*
  //  Second level neighbors
  */
  for(j=0; j<rtuStencilSize-num_neighbors; j++)
    {
      if(levNb[StencilDir1[j]] == level)
	{
	  /*
	  // Our neighbor in the first direction is on the same level,
	  // it is sufficient to take its neighbor in the second direction.
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[StencilDir1[j]],StencilDir2[j]);
	}
      else if(levNb[StencilDir2[j]] == level)
	{
	  /*
	  // Our neighbor in the second direction is on the same level,
	  // it is sufficient to take its neighbor in the first direction.
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[StencilDir2[j]],StencilDir1[j]);
	}
      else
	{
	  /*
	  // Both our neighbors are on a higher level. In that case the corner cell
	  // cannot be our immediate neighbor, it must be a common neighbor of
	  // two higher level cells.
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[StencilDir1[j]],StencilDir2[j]);
#ifdef DEBUG
	  cart_assert(nb[num_neighbors+j] == cell_neighbor(nb[StencilDir2[j]],StencilDir1[j]));
#endif      
	}
    }


#ifdef DEBUG
  double p0[nDim], p[nDim];
  int k;
 
  cell_position_double(cell,p0);

  for(j=0; j<rtuStencilSize-num_neighbors; j++)
    {
      for(k=0; k<nDim; k++)
	{
	  p[k] = p0[k] + cell_size[level]*rtuStencilDelPos[num_neighbors+j][k];
	  if(p[k] < 0.0) p[k] += num_grid;
	  if(p[k] > num_grid) p[k] -= num_grid;
	}

      cart_assert(cell_contains_position(nb[num_neighbors+j],p));
    }
#endif

}


/*
//  Compute the maxium and minimum of a (cached) 1-component array.
*/
void rtuGetLinearArrayMaxMin(int n, float *arr, float *max, float *min)
{
  int j, i, ibeg, iend;
  float *vmax, *vmin;
#ifdef _OPENMP
  int num_pieces = omp_get_num_threads();
#else
  int num_pieces = 1;
#endif 
  int len_piece = (n+num_pieces-1)/num_pieces;
  
  vmax = cart_alloc(float, num_pieces );
  vmin = cart_alloc(float, num_pieces );

#pragma omp parallel for default(none), private(j,i,ibeg,iend), shared(arr,vmin,vmax,n,len_piece,num_pieces)
  for(j=0; j<num_pieces; j++)
    {
      ibeg = j*len_piece;
      iend = ibeg + len_piece;
      if(iend > n) iend = n;
  
      vmin[j] = vmax[j] = arr[ibeg];

      for(i=ibeg+1; i<iend; i++)
	{
	  if(arr[i] > vmax[j]) vmax[j] = arr[i];
	  if(arr[i] < vmin[j]) vmin[j] = arr[i];
	}
    }


  *min = vmin[0];
  *max = vmax[0];
  for(j=1; j<num_pieces; j++)
    {
      if(*max < vmax[j]) *max = vmax[j];
      if(*min > vmin[j]) *min = vmin[j];
    }

  cart_free(vmax);
  cart_free(vmin);
}


void rtuCopyArraysInt(int *dest, int *src, int size)
{
  int i;

  /*
  // Hard-code memcpy for now, but in general we need to check whether
  // doing an OpenMP-parallel loop is faster.
  */
  if(1)
    {
      memcpy(dest,src,sizeof(int)*size);
    }
  else
    {
#pragma omp parallel for default(none), private(i), shared(dest,src,size)
      for(i=0; i<size; i++)
	{
	  dest[i] = src[i];
	}
    }
}


void rtuCopyArraysFloat(float *dest, float *src, int size)
{
  int i;

  /*
  // Hard-code memcpy for now, but in general we need to check whether
  // doing an OpenMP-parallel loop is faster.
  */
  if(1)
    {
      memcpy(dest,src,sizeof(float)*size);
    }
  else
    {
#pragma omp parallel for default(none), private(i), shared(dest,src,size)
      for(i=0; i<size; i++)
	{
	  dest[i] = src[i];
	}
    }
}


#ifdef DEBUG
void rtuCheckGlobalValue(int val, char *name, MPI_Comm local_comm)
{
  int gMin, gMax;

  MPI_Allreduce(&val,&gMin,1,MPI_INT,MPI_MIN,local_comm);
  MPI_Allreduce(&val,&gMax,1,MPI_INT,MPI_MAX,local_comm);

  if(val!=gMin || val!=gMax)
    {
      cart_error("Incorrect global value %s: %d (min: %d, max:%d)",name,val,gMin,gMax);
    }
}
#endif

