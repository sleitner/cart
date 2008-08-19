#include "defs.h"
#ifdef RADIATIVE_TRANSFER


#include "auxiliary.h"
#include "constants.h"
#include "timestep.h"
#include "units.h"

#include "rt_tree.h"
#include "rt_utilities.h"


#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)


#include "rt_transfer.h"

#include <stdio.h>
#include <string.h>


float *BufferFactor, BufferFactor1[max_level-min_level+1];
const int NumVars = 1 + rt_num_et_vars;
const int Vars[] = { RT_VAR_OT_FIELD,
		     rt_et_offset + 0,
		     rt_et_offset + 1,
		     rt_et_offset + 2,
		     rt_et_offset + 3,
		     rt_et_offset + 4,
		     rt_et_offset + 5 };


void rtOtvetComputeGreenFunctions();
void rtOtvetEddingtonTensor(int level);
void rtOtvetLaplacian(int mode, int num_level_cells, int level, float *abc, int *neib, int *info, float abcNum, float *ddL, float *d0L, int maskBit, int *indL2G, int ivarL, int ivarG);

/*
// Tensor stencils
*/
void rtOtvetElementUnitaryTensorCross(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *dd, float *d0);
void rtOtvetElementGenericTensorSplit(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *dd, float *d0);
void rtOtvetElementGenericTensorShell(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *dd, float *d0);
void rtOtvetElementGenericTensorMixed(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *dd, float *d0);


void rtInitRunTransferOtvet()
{
  int i, j;

  BufferFactor = BufferFactor1 - min_level;
  for(i=min_level; i<=max_level; i++) BufferFactor[i] = 1.5;

  if(local_proc_id == MASTER_NODE)
    {
      rtOtvetComputeGreenFunctions();
    }
}


void rtStepBeginTransferOtvet()
{
  static int start = 1;
  int level;

  /*
  // We need to recompute tensors the first time we get here, since
  // they are not saved into a dump file, nor they are required to be
  // set at a start-up.
  */
  if(start)
    {
      start = 0;

      cart_debug("recreating tensors...");

      for(level=min_level; level<=max_level_local(); level++)
	{
	  rtOtvetEddingtonTensor(level);
	}
    }
}


/* 
// Helper macros
*/
#define varL(iL)    cell_var(indL2G[iL],ivarL)
#define varG(iL)    cell_var(indL2G[iL],ivarG)

#define srcL(iL)    cell_rt_source(indL2G[iL])
#define srcG(iL)    rt_source_Avg.Value

#define otfL(iL)    cell_var(indL2G[iL],RT_VAR_OT_FIELD)
#define otfG(iL)    rt_ot_field_Avg.Value


void rtLevelUpdateTransferOtvet(int level, MPI_Comm local_comm)
{
  /*
  // Numerical parameters
  */
  float epsAlp = 0.1;
  float safety = 0.8;
  float tauMin = 1.0e-3;
  int numIterMax = 100;
  int numIterMin = 1;

  /*
  // Extra arrays
  */
  int *indL2G, *indG2L, *neib, *info, *tmp;
  float *abc[2], *ddL, *ddG, *d0L, *d0G;

  /*
  // Work variables
  */
  int work, *level_cells;
  float xiUnit, xiFin, dxi, dx2;
  int num_level_cells, num_all_cells, num_total_cells, *nb;
  int iL, iG, j, offset;
  int var_indicies[2];
  int nit, nit0, ifreq, ivarL, ivarG, it;
  float abcMin, abcMax, abcNum;
  float abcMin1, abcMax1;
  float abcL, abcG, facL, facG;
  float abcLmid[num_neighbors], abcGmid[num_neighbors];

#ifdef RT_DEBUG
  double pos[3];
#endif


  xiUnit = (mpc/clight/year)*r0/hubble/(t0*aexp[level]);
  xiFin = 0.5*dtl[level]/xiUnit;

  /*
  // If there are no local cells, we do a dry run only
  */
  work = (num_cells_per_level[level] > 0);
  if(work)
    {
      /*
      //  Insure that our tensors are up-to-date; out parent tensors must be up-to-date.
      */
      rtOtvetEddingtonTensor(level);

      /* 
      //  Allocate memory for index arrays
      */
      indG2L = cart_alloc(num_cells*sizeof(int));  /* NG: I do not know how to avoid using this large array */
      select_level(level,CELL_TYPE_LOCAL,&num_all_cells,&level_cells);

      /*
      //  We got all cells, but we only need leaves!!!
      //  (If we solve levels separately, we can never insure that 
      //  light fronts on both levels propagate at the same speed,
      //  even if I-fronts do.)
      */
      info = (int *)cart_alloc(num_all_cells*sizeof(int));
      num_level_cells = 0;
      for(iL=0; iL<num_all_cells; iL++) if(cell_is_leaf(level_cells[iL]))
	{
	  info[num_level_cells++] = level_cells[iL];
	}
      cart_free(level_cells);

      /*
      // Are we wasting too much memory?
      */
      if(num_level_cells < num_all_cells/10)
	{
	  level_cells = (int *)cart_alloc(num_level_cells*sizeof(int));
	  rtuCopyArraysInt(level_cells,info,num_level_cells);
	  cart_free(info);
	  info = level_cells;
	}

      /*
      // Allocate the rest of arrays
      */
      neib = (int *)cart_alloc(rtuStencilSize*num_level_cells*sizeof(int));

      num_all_cells = 1 + (int)(BufferFactor[level]*num_level_cells);
      indL2G = (int *)cart_alloc(num_all_cells*sizeof(int));
      rtuCopyArraysInt(indL2G,info,num_level_cells);

      /*
      // Initialize the global-to-local index array
      */
#pragma omp parallel for default(none), private(iG), shared(indG2L)
      for(iG=0; iG<num_cells; iG++)
	{
	  indG2L[iG] = -1;
	}

      /*
      // Compute global-to-local mapping indices
      */
#pragma omp parallel for default(none), private(iL), shared(indL2G,indG2L,num_level_cells)
      for(iL=0; iL<num_level_cells; iL++)
	{
	  indG2L[indL2G[iL]] = iL;
	}

      /*
      //  Compute neighbors as global indices
      */
#pragma omp parallel for default(none), private(iL,offset,j), shared(neib,info,num_level_cells,indL2G,level)
      for(iL=0; iL<num_level_cells; iL++ )
	{
	  offset = rtuStencilSize*iL;
	  
	  rtuGetStencil(level,indL2G[iL],neib+offset);
	  
	  info[iL] = 0;
	  for(j=0; j<2*nDim; j++)
	    {
	      if(cell_level(neib[offset+j]) < level) info[iL] |= (1<<j);
	    }
	}

      /*
      // Distribute members as local indicies
      // THIS IS A SERIAL LOOP!!!!
      */
      num_total_cells = num_level_cells;
      for(iL=0; iL<num_level_cells; iL++ )
	{
	  nb = neib + rtuStencilSize*iL;
	  
	  for(j=0; j<rtuStencilSize; j++)
	    {
	      it = indG2L[nb[j]];
	      if(it >= 0)
		{
		  nb[j] = it;
		}
	      else
		{
		  
		  cart_assert(num_total_cells <= num_all_cells);
		  
		  /*
		  // Check that we are not overfilling the allocated memory - a-la IFRIT's self-extendable arrays
		  */
		  if(num_total_cells == num_all_cells)
		    {
		      nit = num_all_cells;
		      BufferFactor[level] *= 1.5;
		      num_all_cells = 1 + (int)(BufferFactor[level]*num_level_cells);
		      
		      cart_debug("Extending local-to-global index buffer (%d -> %d) to margin %f",nit,num_all_cells,BufferFactor[level]);
		      
		      tmp = (int *)cart_alloc(num_all_cells*sizeof(int));
		      rtuCopyArraysInt(tmp,indL2G,nit);
		      cart_free(indL2G);
		      indL2G = tmp;
		    }
		  
		  indG2L[nb[j]] = num_total_cells;
		  indL2G[num_total_cells] = nb[j];
		  nb[j] = num_total_cells;
		  num_total_cells++;
		}
	    }
	}
      
      ddL = (float *)cart_alloc(num_level_cells*sizeof(float));
      ddG = (float *)cart_alloc(num_level_cells*sizeof(float));
      d0L = (float *)cart_alloc(num_level_cells*sizeof(float));
      d0G = (float *)cart_alloc(num_level_cells*sizeof(float));

      /*
      // Are we using too much memory?
      */
      if(num_total_cells<0.5*num_all_cells && BufferFactor[level]>1.5/0.8)
	{
	  BufferFactor[level] *= 0.8;
	}
      
#ifdef RTO_CACHE_VARS
  /*
C
C  Cache the Eddington Tensor, the OT solution, and the source
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(iL),
C$OMP+SHARED(var,varET,varOT,varSR,indLG,nTotCells)
      do iL=1,nTotCells
         varOT(iL) = var(ifrOTn,indLG(iL))
         varET(1,iL) = var(irtET1,indLG(iL))
         varET(2,iL) = var(irtET2,indLG(iL))
         varET(3,iL) = var(irtET3,indLG(iL))
         varET(4,iL) = var(irtET4,indLG(iL))
         varET(5,iL) = var(irtET5,indLG(iL))
         varET(6,iL) = var(irtET6,indLG(iL))
         varSR(iL) = var(irtSor,indLG(iL))
      enddo
  */
#endif

    }
  else
    {

      select_level(level,CELL_TYPE_ANY,&num_total_cells,&indL2G);
      cart_assert(num_total_cells > 0);

    }

  /*
  // Allocate memory for absorption coefficient and other arrays
  */
  abc[0] = (float *)cart_alloc(num_total_cells*sizeof(float));
#if (RT_CFI == 1)
  abc[1] = (float *)cart_alloc(num_total_cells*sizeof(float));
#else
  abc[1] = 0;
#endif

  /*
  // Iterate over all non-zero frequencies, two at a time (local and global fields)
  */
  for(ifreq=0; ifreq<rt_num_frequencies/2; ifreq++)
    {

      ivarL = rt_freq_offset + ifreq;
      ivarG = ivarL + rt_num_frequencies/2;

      var_indicies[0] = ivarL;
      var_indicies[1] = ivarG;

#ifdef RTO_CACHE_VARS
  /*
C
C  Cache the radiation fields
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(iL,j),
C$OMP+SHARED(var,varLN,varLF,ivarN,ivarF,indLG,nTotCells)
      do iL=1,nTotCells
         varLN(iL) = var(ivarN,indLG(iL))
         varLF(iL) = var(ivarF,indLG(iL))
      enddo
  */
#endif

      /*
      // Compute the absorption coefficient at this level
      */
      rtComputeAbsLevel(num_total_cells,indL2G,ifreq,abc);

      /*
      // Numerical "diffusion" (actually, inverse one)
      */
      abcNum = epsAlp/cell_size[0];

      /*
      // Compute the time-step
      */
      rtuGetLinearArrayMaxMin(num_total_cells,abc[0],&abcMax,&abcMin);
#if (RT_CFI == 1)
      rtuGetLinearArrayMaxMin(num_total_cells,abc[1],&abcMax1,&abcMin1);
      if(abcMin > abcMin1) abcMin = abcMin1;
      if(abcMax < abcMax1) abcMax = abcMax1;
#endif

      /*
      //  Make sure all procs have the same abcMin/Max
      */
      MPI_Allreduce(&abcMin,&abcMin1,1,MPI_FLOAT,MPI_MIN,local_comm);
      MPI_Allreduce(&abcMax,&abcMax1,1,MPI_FLOAT,MPI_MAX,local_comm);
      abcMin = abcMin1;
      abcMax = abcMax1;

      /*
      // Is the box optically thin?
      */
      if(abcMax*num_grid < tauMin)
	{
	  if(work)
	    {
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,rt_glob_Avg)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  varL(iL) = otfL(iL);
		  varG(iL) = otfG(iL);
		}
	    }

	  it = 0;
	  
	}
      else
	{

	  abcMin = abcMin + abcNum;

	  dxi = safety*0.5*cell_size[level]*cell_size[level]*abcMin/3.0;

	  nit = 1 + (int)(xiFin/dxi);
	  dxi = xiFin/nit;
	  dx2 = 1.0/(cell_size[level]*cell_size[level]);

#ifdef RT_SIGNALSPEED_TO_C
	  /*
	  // Number of iterations needed for the signal propagation speed to be equal c.
	  */
	  nit0 = 1 + (int)(2*xiFin/cell_size[level]);
	  if(nit > nit0) nit = nit0;
#endif

	  nit0 = nit;
	  if(nit < numIterMin) nit = numIterMin;
	  if(nit > numIterMax) nit = numIterMax;

#ifdef RT_OUTPUT
	  MPI_Comm_rank(local_comm,&it);
	  if(it == 0)
	    { 
	      cart_debug("OTVET step: level=%2d, freq=%d, nit=%d, nit0=%d",level,ifreq,nit,nit0);
	    }
#endif

#ifdef DEBUG
	  rtuCheckGlobalValue(ifreq,"Otvet:ifreq");
	  rtuCheckGlobalValue(nit0,"Otvet:nit0");
	  rtuCheckGlobalValue(nit,"Otvet:nit");
#endif

	  for(it=0; it<nit; it++)
	    {

	      if(work)
		{

		  /*
		  // Compute Laplacian term
		  */
#pragma omp parallel for default(none), private(iL,nb,abcLmid,abcGmid,j), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,abc,neib,level,abcNum,info,d0L,ddL,d0G,ddG,dx2)
		  for(iL=0; iL<num_level_cells; iL++)
		    {

		      nb = neib + rtuStencilSize*iL;

		      /*
		      // Flux-conserving diffusion coefficients a-la Oran & Boris
		      */
		      if(level == min_level)
			{
			  for(j=0; j<num_neighbors; j++)
			    {
#if (RT_CFI == 1)
			      abcLmid[j] = 0.5*(abc[0][iL]+abc[0][nb[j]]) + abcNum;
			      abcGmid[j] = 0.5*(abc[1][iL]+abc[1][nb[j]]) + abcNum;
#else
			      abcLmid[j] = abcGmid[j] = 0.5*(abc[0][iL]+abc[0][nb[j]]) + abcNum;
#endif
			    }
			}
		      else
			{
			  for(j=0; j<num_neighbors; j++)
			    {
			      if(info[iL] & (1<<j))
				{
#if (RT_CFI == 1)
				  abcLmid[j] = 0.333*abc[0][iL] + 0.667*abc[0][nb[j]] + abcNum;
				  abcGmid[j] = 0.333*abc[1][iL] + 0.667*abc[1][nb[j]] + abcNum;
#else
				  abcLmid[j] = abcGmid[j] = 0.333*abc[0][iL] + 0.667*abc[0][nb[j]] + abcNum;
#endif
				}
			      else
				{
#if (RT_CFI == 1)
				  abcLmid[j] = 0.5*(abc[0][iL]+abc[0][nb[j]]) + abcNum;
				  abcGmid[j] = 0.5*(abc[1][iL]+abc[1][nb[j]]) + abcNum;
#else
				  abcLmid[j] = abcGmid[j] = 0.5*(abc[0][iL]+abc[0][nb[j]]) + abcNum;
#endif
				}
			    }
			}

		      /*
		      // Laplacian and its diagonal element
		      */
		      rtOtvetElementGenericTensorMixed(ivarL,iL,indL2G,nb,abcLmid,ddL+iL,d0L+iL);
		      rtOtvetElementUnitaryTensorCross(ivarG,iL,indL2G,nb,abcGmid,ddG+iL,d0G+iL);

		      ddL[iL] *= dx2;
		      d0L[iL] *= dx2;
		      ddG[iL] *= dx2;
		      d0G[iL] *= dx2;
		    }

		  /*
		  // Update the radiation fields
		  */
#pragma omp parallel for default(none), private(iL,abcL,abcG,facL,facG), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,abc,dxi,d0L,d0G,ddL,ddG,rt_glob_Avg)
		  for(iL=0; iL<num_level_cells; iL++)
		    {
		      abcL = abc[0][iL];
#if (RT_CFI == 1)
		      abcG = abc[1][iL];
#else
		      abcG = abc[0][iL];
#endif

		      facL = dxi/(1.0+dxi*(abcL+d0L[iL]));
		      facG = dxi/(1.0+dxi*(abcG+d0G[iL]));

		      varL(iL) += facL*(ddL[iL]-abcL*varL(iL)+srcL(iL)) /* *(1.0-varL(iL)/otfL(iL)) */;
		      varG(iL) += facG*(ddG[iL]-abcG*varG(iL)+srcG(iL)) /* *(1.0-varG(iL)/otfG(iL)) */;

		      /*
		      // Limit
		      */
		      if(varL(iL) < 0.0) varL(iL) = 0.0;
		      if(varL(iL) > otfL(iL)) varL(iL) = otfL(iL);
		      if(varG(iL) < 0.0) varG(iL) = 0.0;
		      if(varG(iL) > otfG(iL)) varG(iL) = otfG(iL);

		    }

		}
#ifdef DEBUG
	      rtuCheckGlobalValue(it,"Otvet:it");
#endif
	      update_buffer_level(level,var_indicies,2);
	    }

#ifdef RTO_CACHE_VARS
  /*
C
C  Un-cache the radiation fields
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(iL,j),
C$OMP+SHARED(nLevCells,var,varLN,varLF,ivarN,ivarF,indLG,nTotCells)
         do iL=1,nLevCells
            var(ivarN,indLG(iL)) = varLN(iL)
            var(ivarF,indLG(iL)) = varLF(iL)
         enddo
  */
#endif
	}
    }

  cart_free(indL2G);
  cart_free(abc[0]);
#if (RT_CFI == 1)
  cart_free(abc[1]);
#endif

  if(work)
    {
      cart_free(indG2L);
      cart_free(neib);
      cart_free(info);
      cart_free(ddL);
      cart_free(ddG);
      cart_free(d0L);
      cart_free(d0G);
    }

#ifdef RT_DEBUG
  /* for(j=0; j<=4*num_grid/2; j++)
    {
      pos[0] = 0.25*(j+0.5);
      pos[1] = pos[2] = 0.5*num_grid - 0.125;
      iL = cell_find_position(pos);
      cart_debug("%d %15g %15g %15g",j,cell_var(iL,rt_freq_offset+0),cell_var(iL,rt_freq_offset+1),cell_var(iL,rt_freq_offset+2));
    }
  /*   exit(0); */
#endif

}


/*
// Compute the Laplacian and its diagonal element in several ways
*/
#define var(iL)    cell_var(indL2G[iL],ivar)
#define varET(ind,iL)    cell_var(indL2G[iL],rt_et_offset+ind)
#define RPT(ind,iL)      (var(iL)*varET(ind,iL))


void rtOtvetElementUnitaryTensorCross(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *dd, float *d0)
{
  int j;
  float s, vmax;

#ifndef RT_DEBUG_BLOCK_MASKING
  /*
  // Compute update mask
  */
  vmax = var(iL);
  for(j=0; j<rtuStencilSize; j++)
    {
      if(vmax < var(nb[j])) vmax = var(nb[j]);
    }
  if(vmax > 1.0e-35)
#endif
    {
      s = 0.0;
      for(j=0; j<num_neighbors; j++) s += (var(nb[j])-var(iL))/abcLoc[j];
      *dd = s/nDim;
    }
#ifndef RT_DEBUG_BLOCK_MASKING
  else *dd = 0.0;
#endif
  
  /* Do always */
    {
      s = 0.0;
      for(j=0; j<num_neighbors; j++) s += 1.0/abcLoc[j];
      *d0 = s/nDim;
    }
}


void rtOtvetElementGenericTensorSplit(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *dd, float *d0)
{
  int j;
  float vmax;

#ifndef RT_DEBUG_BLOCK_MASKING
  /*
  // Compute update mask
  */
  vmax = var(iL);
  for(j=0; j<rtuStencilSize; j++)
    {
      if(vmax < var(nb[j])) vmax = var(nb[j]);
    }
  if(vmax > 1.0e-35)
#endif
    {
      /*
      // Diagonal elements
      */
      *dd = 
	(RPT(0,nb[0])-RPT(0,iL))/abcLoc[0] +
	(RPT(0,nb[1])-RPT(0,iL))/abcLoc[1]
#if (nDim > 1)
	+
	(RPT(2,nb[2])-RPT(2,iL))/abcLoc[2] +
	(RPT(2,nb[3])-RPT(2,iL))/abcLoc[3]
#if (nDim > 2)
#endif /* nDim > 2 */
	+
	(RPT(5,nb[4])-RPT(5,iL))/abcLoc[4] +
	(RPT(5,nb[5])-RPT(5,iL))/abcLoc[5]
#endif /* nDim > 1 */
	;

	/*
	// Off-diagonal elements
	*/
#if (nDim > 1)
	*dd += 1.0*0.25*(
		     (RPT(1,nb[ 9])+RPT(1,nb[ 1])-RPT(1,nb[ 8])-RPT(1,nb[ 0]))/abcLoc[3] -
		     (RPT(1,nb[ 7])+RPT(1,nb[ 1])-RPT(1,nb[ 6])-RPT(1,nb[ 0]))/abcLoc[2] +
		     (RPT(1,nb[ 9])+RPT(1,nb[ 3])-RPT(1,nb[ 7])-RPT(1,nb[ 2]))/abcLoc[1] -
		     (RPT(1,nb[ 8])+RPT(1,nb[ 3])-RPT(1,nb[ 6])-RPT(1,nb[ 2]))/abcLoc[0]
#if (nDim > 2)
		     +
		     (RPT(3,nb[15])+RPT(3,nb[ 1])-RPT(3,nb[14])-RPT(3,nb[ 0]))/abcLoc[5] -
		     (RPT(3,nb[11])+RPT(3,nb[ 1])-RPT(3,nb[10])-RPT(3,nb[ 0]))/abcLoc[4] +
		     (RPT(3,nb[15])+RPT(3,nb[ 5])-RPT(3,nb[11])-RPT(3,nb[ 4]))/abcLoc[1] -
		     (RPT(3,nb[14])+RPT(3,nb[ 5])-RPT(3,nb[10])-RPT(3,nb[ 4]))/abcLoc[0] +
		     (RPT(4,nb[17])+RPT(4,nb[ 3])-RPT(4,nb[16])-RPT(4,nb[ 2]))/abcLoc[5] -
		     (RPT(4,nb[13])+RPT(4,nb[ 3])-RPT(4,nb[12])-RPT(4,nb[ 2]))/abcLoc[4] +
		     (RPT(4,nb[17])+RPT(4,nb[ 5])-RPT(4,nb[13])-RPT(4,nb[ 4]))/abcLoc[3] -
		     (RPT(4,nb[16])+RPT(4,nb[ 5])-RPT(4,nb[12])-RPT(4,nb[ 4]))/abcLoc[2]
#endif /* nDim > 2 */
		     );
#endif /* nDim > 1 */
    }
#ifndef RT_DEBUG_BLOCK_MASKING
  else *dd = 0.0;
#endif
  
  /* Do always */
    {
      *d0 =  varET(0,iL)*(1.0/abcLoc[0]+1.0/abcLoc[1])
#if (nDim > 1)
	+ varET(2,iL)*(1.0/abcLoc[2]+1.0/abcLoc[3])
#if (nDim > 2)
	+ varET(5,iL)*(1.0/abcLoc[4]+1.0/abcLoc[5])
#endif /* nDim > 2 */
#endif /* nDim > 1 */
	;
    }
}



#if (nDim == 3)

#define SXp(p,i) (p[1][i]+p[3][i]+p[5][i]+p[7][i])
#define SXm(p,i) (p[0][i]+p[2][i]+p[4][i]+p[6][i])

#define SYp(p,i) (p[2][i]+p[3][i]+p[6][i]+p[7][i])
#define SYm(p,i) (p[0][i]+p[1][i]+p[4][i]+p[5][i])

#define SZp(p,i) (p[4][i]+p[5][i]+p[6][i]+p[7][i])
#define SZm(p,i) (p[0][i]+p[1][i]+p[2][i]+p[3][i])

#define DX(p,i) (SXp(p,i)-SXm(p,i))
#define DY(p,i) (SYp(p,i)-SYm(p,i))
#define DZ(p,i) (SZp(p,i)-SZm(p,i))

void rtOtvetElementGenericTensorShell(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *dd, float *d0)
{
  const int index[8][6] = {
    { 0, 2, 6, 4, 10, 12 },
    { 1, 2, 7, 4, 11, 12 },
    { 0, 3, 8, 4, 10, 13 },
    { 1, 3, 9, 4, 11, 13 },
    { 0, 2, 6, 5, 14, 16 },
    { 1, 2, 7, 5, 15, 16 },
    { 0, 3, 8, 5, 14, 17 },
    { 1, 3, 9, 5, 15, 17 } };
    
  const int dir[8][3] = { /* 0x1, 0x2, 1x2 */
    {  1,  1,  1 },    /* { -1, -1, -1 } */
    { -1, -1,  1 },    /* {  1, -1, -1 } */
    { -1,  1, -1 },    /* { -1,  1, -1 } */
    {  1, -1, -1 },    /* {  1,  1, -1 } */
    {  1, -1, -1 },    /* { -1, -1,  1 } */
    { -1,  1, -1 },    /* {  1, -1,  1 } */
    { -1, -1,  1 },    /* { -1,  1,  1 } */
    {  1,  1,  1 } };  /* {  1,  1,  1 } */

  int i, j, k;
  float p[8][rt_num_et_vars];
  float vmax, s[8][nDim];
 
#ifndef RT_DEBUG_BLOCK_MASKING
  /*
  // Compute update mask
  */
  vmax = var(iL);
  for(j=0; j<rtuStencilSize; j++)
    {
      if(vmax < var(nb[j])) vmax = var(nb[j]);
    }
  if(vmax > 1.0e-35)
#endif
    {
      for(j=0; j<rt_num_et_vars; j++) p[0][j] = RPT(j,iL);

      for(i=0; i<8; i++)
	{
	  for(j=0; j<rt_num_et_vars; j++)
	    {
	      for(k=0; k<6; k++) p[k+1][j] = RPT(j,nb[index[i][k]]);
	      p[7][j] = p[0][j] + p[3][j] + p[5][j] + p[6][j] - p[1][j] - p[2][j] - p[4][j];
	    }
	  
	  s[i][0] = 0.25*(DX(p,0)+dir[i][0]*DY(p,1)+dir[i][1]*DZ(p,3));
	  s[i][1] = 0.25*(dir[i][0]*DX(p,1)+DY(p,2)+dir[i][2]*DZ(p,4));
	  s[i][2] = 0.25*(dir[i][1]*DX(p,3)+dir[i][2]*DY(p,4)+DZ(p,5));
	}

      *dd = 0.25*(
		   SXp(s,0)/abcLoc[1] + SXm(s,0)/abcLoc[0] +
		   SYp(s,1)/abcLoc[3] + SYm(s,1)/abcLoc[2] +
		   SZp(s,2)/abcLoc[5] + SZm(s,2)/abcLoc[4]);
    }
#ifndef RT_DEBUG_BLOCK_MASKING
  else *dd = 0.0;
#endif

  /* Do always */
    {
      *d0 = 0.25*(
	varET(0,iL)*(1.0/abcLoc[1]+1.0/abcLoc[0]) +
	varET(2,iL)*(1.0/abcLoc[3]+1.0/abcLoc[2]) +
	varET(5,iL)*(1.0/abcLoc[5]+1.0/abcLoc[4]));
    }
}


void rtOtvetElementGenericTensorMixed(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *dd, float *d0)
{
  const float eps = 3.0;

  const int index[8][6] = {
    { 0, 2, 6, 4, 10, 12 },
    { 1, 2, 7, 4, 11, 12 },
    { 0, 3, 8, 4, 10, 13 },
    { 1, 3, 9, 4, 11, 13 },
    { 0, 2, 6, 5, 14, 16 },
    { 1, 2, 7, 5, 15, 16 },
    { 0, 3, 8, 5, 14, 17 },
    { 1, 3, 9, 5, 15, 17 } };
    
  const int dir[8][3] = { /* 0x1, 0x2, 1x2 */
    {  1,  1,  1 },    /* { -1, -1, -1 } */
    { -1, -1,  1 },    /* {  1, -1, -1 } */
    { -1,  1, -1 },    /* { -1,  1, -1 } */
    {  1, -1, -1 },    /* {  1,  1, -1 } */
    {  1, -1, -1 },    /* { -1, -1,  1 } */
    { -1,  1, -1 },    /* {  1, -1,  1 } */
    { -1, -1,  1 },    /* { -1,  1,  1 } */
    {  1,  1,  1 } };  /* {  1,  1,  1 } */

  const float diag[6] = { 1.0, 0.0, 1.0, 0.0, 0.0, 1.0 };

  int i, j, k;
  float p[8][rt_num_et_vars];
  float vmax, q, s[8][nDim];

#ifndef RT_DEBUG_BLOCK_MASKING
  /*
  // Compute update mask
  */
  vmax = var(iL);
  for(j=0; j<rtuStencilSize; j++)
    {
      if(vmax < var(nb[j])) vmax = var(nb[j]);
    }
  if(vmax > 1.0e-35)
#endif
    {
      q = 0.0;
      for(j=0; j<num_neighbors; j++) q += (var(nb[j])-var(iL))/abcLoc[j];
      *dd = eps*q/nDim;

      for(j=0; j<rt_num_et_vars; j++) p[0][j] = RPT(j,iL) - eps*diag[j]/3.0;

      for(i=0; i<8; i++)
	{
	  for(j=0; j<rt_num_et_vars; j++)
	    {
	      for(k=0; k<6; k++) p[k+1][j] = RPT(j,nb[index[i][k]]) - eps*diag[j]/3.0;
	      p[7][j] = p[0][j] + p[3][j] + p[5][j] + p[6][j] - p[1][j] - p[2][j] - p[4][j];
	    }
	  
	  s[i][0] = 0.25*(DX(p,0)+dir[i][0]*DY(p,1)+dir[i][1]*DZ(p,3));
	  s[i][1] = 0.25*(dir[i][0]*DX(p,1)+DY(p,2)+dir[i][2]*DZ(p,4));
	  s[i][2] = 0.25*(dir[i][1]*DX(p,3)+dir[i][2]*DY(p,4)+DZ(p,5));
	}

      *dd += 0.25*(
		   SXp(s,0)/abcLoc[1] + SXm(s,0)/abcLoc[0] +
		   SYp(s,1)/abcLoc[3] + SYm(s,1)/abcLoc[2] +
		   SZp(s,2)/abcLoc[5] + SZm(s,2)/abcLoc[4]);
    }
#ifndef RT_DEBUG_BLOCK_MASKING
  else *dd = 0.0;
#endif

  /* Do always */
    {
      q = 0.0;
      for(j=0; j<num_neighbors; j++) q += 1.0/abcLoc[j];
      *d0 = eps*q/nDim;

      *d0 += 0.25*(
	(varET(0,iL)-eps/3.0)*(1.0/abcLoc[1]+1.0/abcLoc[0]) +
	(varET(2,iL)-eps/3.0)*(1.0/abcLoc[3]+1.0/abcLoc[2]) +
	(varET(5,iL)-eps/3.0)*(1.0/abcLoc[5]+1.0/abcLoc[4]));
    }
}

#endif /* nDim == 3 */



#endif  // defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
#endif  // RADIATIVE_TRANSFER
