#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "parallel.h"
#include "rt_global.h"
#include "rt_transfer.h"
#include "rt_utilities.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"



int rtOtvetMaxNumIter = 10;


int NumVars = 1 + rt_num_et_vars;
int Vars[] = { RT_VAR_OT_FIELD,
		       rt_et_offset + 0,
		       rt_et_offset + 1,
		       rt_et_offset + 2,
		       rt_et_offset + 3,
		       rt_et_offset + 4,
		       rt_et_offset + 5 };

DEFINE_LEVEL_ARRAY(float,BufferFactor);


void rtOtvetComputeGreenFunctions();
void rtOtvetEddingtonTensor(int level, int num_level_cells, int *level_cells);

void rtOtvetMidPointAbsorptionCoefficients(int level, int iL, int *info, int *nb, float *abc, float *abcMid);


typedef float (*LaplacianDiag)(int iL, int *indL2G, float *abcLoc);
typedef float (*LaplacianFull)(int ivar, int iL, int *indL2G, int *nb, float *abcLoc);

void rtOtvetSolveFieldEquation(int ivar, int level, int num_level_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd, int nit, int work, LaplacianDiag dWorker, LaplacianFull fWorker);

/*
// Tensor stencils
*/
float rtOtvetLaplacian_UnitaryTensor_CrossStencil_Diag(int iL, int *indL2G, float *abcLoc);
float rtOtvetLaplacian_GenericTensor_SplitStencil_Diag(int iL, int *indL2G, float *abcLoc);

float rtOtvetLaplacian_UnitaryTensor_CrossStencil_Full(int ivar, int iL, int *indL2G, int *nb, float *abcLoc);
float rtOtvetLaplacian_GenericTensor_SplitStencil_Full(int ivar, int iL, int *indL2G, int *nb, float *abcLoc);


void rtInitRunTransferOtvet()
{
  int i;

  for(i=min_level; i<=max_level; i++)
    {
      BufferFactor[i] = 2.0;
    }

  if(local_proc_id == MASTER_NODE)
    {
      rtOtvetComputeGreenFunctions();
    }
}


void rtStepBeginTransferOtvet()
{
}


void rtAfterAssignDensityTransferOtvet(int level, int num_level_cells, int *level_cells)
{
  rtOtvetEddingtonTensor(level,num_level_cells,level_cells);
}


/* 
// Helper macros
*/
#define varL(iL)    cell_var(indL2G[iL],ivarL)
#define varG(iL)    cell_var(indL2G[iL],ivarG)

#define srcL(iL)    cell_rt_source(indL2G[iL])
#define srcG(iL)    rtGlobals[RT_SOURCE_AVG].Value

#define otfL(iL)    cell_var(indL2G[iL],RT_VAR_OT_FIELD)
#define otfG(iL)    rtGlobals[RT_OT_FIELD_AVG].Value


void rtLevelUpdateTransferOtvet(int level, MPI_Comm local_comm)
{
  const float tauMin = 1.0e-2;

  /*
  // Extra arrays
  */
  int *indL2G, *indG2L, *neib, *info, *tmp;
  float *abc[2], *dd, *jac, *rhs;

  /*
  // Work variables
  */
  int work;
  int num_level_cells, num_all_cells, num_total_cells, *nb;
  int iL, iG, j, offset, index;
  int nit, ifreq, ivarL, ivarG;
  float abcMin, abcMax;
  float abcMin1, abcMax1;
#ifdef RT_SIGNALSPEED_TO_C
  int nit0;
  float xiUnit;
#endif


  /*
  // If there are no local cells, we do a dry run only
  */
  work = (num_cells_per_level[level] > 0);
  if(work)
    {
      /* 
      //  Allocate memory for index arrays
      */
      indG2L = cart_alloc(int, num_cells );  /* NG: I do not know how to avoid using this large array */

      /*
      //  Find all leaves (if we solve levels separately, we can never insure that 
      //  light fronts on both levels propagate at the same speed, even if I-fronts do.)
      */
      select_level_with_condition(1,level,&num_level_cells,&info);

      /*
      // Allocate the rest of arrays
      */
      neib = cart_alloc(int, rtuStencilSize*num_level_cells );
      num_all_cells = 100 + (int)(BufferFactor[level]*num_level_cells);
      indL2G = cart_alloc(int, num_all_cells );
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
	      index = indG2L[nb[j]];
	      if(index >= 0)
		{
		  nb[j] = index;
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
		      num_all_cells = 100 + (int)(BufferFactor[level]*num_level_cells);
		      
		      cart_debug("Extending local-to-global index buffer (%d -> %d) to margin %f",nit,num_all_cells,BufferFactor[level]);
		      
		      tmp = cart_alloc(int, num_all_cells );
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

      /*
      //  It is not needed any more
      */
      cart_free(indG2L);

      rhs = cart_alloc(float, num_level_cells );
      jac = cart_alloc(float, num_level_cells );
      dd  = cart_alloc(float, num_level_cells );

      /*
      // Are we using too much memory?
      */
      if(num_total_cells<0.5*num_all_cells && BufferFactor[level]>2.0/0.8)
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
  abc[0] = cart_alloc(float, num_total_cells );
#if (RT_CFI == 1)
  abc[1] = cart_alloc(float, num_total_cells );
#else
  abc[1] = abc[0];
#endif

  /*
  // Iterate over all non-zero frequencies, two at a time (local and global fields)
  */
  for(ifreq=0; ifreq<rt_num_frequencies/2; ifreq++)
    {

      ivarL = rt_freq_offset + ifreq;
      ivarG = ivarL + rt_num_frequencies/2;

      /*
      // Compute the absorption coefficient at this level
      */
      rtComputeAbsLevel(num_total_cells,indL2G,ifreq,abc);

      /*
      // Compute the range of the abc array(s)
      */
      if(num_total_cells > 0)
	{
	  rtuGetLinearArrayMaxMin(num_total_cells,abc[0],&abcMax,&abcMin);
#if (RT_CFI == 1)
	  rtuGetLinearArrayMaxMin(num_total_cells,abc[1],&abcMax1,&abcMin1);
	  if(abcMin > abcMin1) abcMin = abcMin1;
	  if(abcMax < abcMax1) abcMax = abcMax1;
#endif
	}
      else
	{
	  abcMin = 1.0e35;
	  abcMax = 0.0;
	}

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
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,rtGlobals)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  varL(iL) = otfL(iL);
		  varG(iL) = otfG(iL);
		}
	    }
	}
      else
	{
	  /*
	  //  Number of iterations
	  */
	  nit = rtOtvetMaxNumIter;

#ifdef RT_SIGNALSPEED_TO_C
	  /*
	  // Number of iterations needed for the signal propagation speed to less or equal c.
	  */
	  xiUnit = units->length/(constants->c*units->time);
	  nit0 = 1 + (int)(dtl[level]/xiUnit/cell_size[level]);
	  if(nit > nit0)
	    {
	      cart_debug("OTVET: reducing iteration count from %d to %d",nit,nit0);
	      nit = nit0;
	    }
#endif

#ifdef DEBUG
	  rtuCheckGlobalValue(ifreq,"Otvet:ifreq",local_comm);
#ifdef RT_SIGNALSPEED_TO_C
	  rtuCheckGlobalValue(nit0,"Otvet:nit0",local_comm);
#endif
	  rtuCheckGlobalValue(nit,"Otvet:nit",local_comm);
#endif
	  /*
	  // Update local field
	  */
	  if(work)
	    {
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,level,rtGlobals,rhs)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  rhs[iL] = srcL(iL);
		}
	    }

	  rtOtvetSolveFieldEquation(ivarL,level,num_level_cells,indL2G,neib,info,abc[0],rhs,jac,dd,nit,work,rtOtvetLaplacian_GenericTensor_SplitStencil_Diag,rtOtvetLaplacian_GenericTensor_SplitStencil_Full);

#ifdef RT_EXTERNAL_BACKGROUND
	  /*
	  // Update global field
	  */
	  if(work)
	    {
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,level,rtGlobals,rhs)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  rhs[iL] = srcG(iL);
		}
	    }

	  rtOtvetSolveFieldEquation(ivarG,level,num_level_cells,indL2G,neib,info,abc[1],rhs,jac,dd,nit,work,rtOtvetLaplacian_UnitaryTensor_CrossStencil_Diag,rtOtvetLaplacian_UnitaryTensor_CrossStencil_Full);

#endif
	  /*
	  // Limit
	  */
	  if(work)
	    {
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,rtGlobals)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  if(varL(iL) < 0.0) varL(iL) = 0.0;
		  if(varL(iL) > otfL(iL)) varL(iL) = otfL(iL);
#ifdef RT_EXTERNAL_BACKGROUND
		  if(varG(iL) < 0.0) varG(iL) = 0.0;
		  if(varG(iL) > otfG(iL)) varG(iL) = otfG(iL);
#else
		  varG(iL) = 0.0;
#endif
		}
	    }
	}
    }

  cart_free(indL2G);
  cart_free(abc[0]);
#if (RT_CFI == 1)
  cart_free(abc[1]);
#endif

  if(work)
    {
      cart_free(neib);
      cart_free(info);
      cart_free(rhs);
      cart_free(jac);
      cart_free(dd);
    }
}


/* 
// Helper macros
*/
#define var(iL)    cell_var(indL2G[iL],ivar)
#define varET(ind,iL)    cell_var(indL2G[iL],rt_et_offset+ind)
#define RPT(ind,iL)      (varET(ind,iL)*var(iL))


void rtOtvetSolveFieldEquation(int ivar, int level, int num_level_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd, int nit, int work, LaplacianDiag dWorker, LaplacianFull fWorker)
{
  /*
  // Numerical parameters
  // alpha should be: 
  //   (1) sufficiently smaller than 1 for stability
  //   (2) not too small for efficiency
  */
  const float alpha = 0.5;
  const float beta  = 1.0;
  const float gamma = 1.0;
  const float epsNum = 1.0e-3;

  int j, it, iL, *nb;
  float abcMid[num_neighbors];
  float dx = cell_size[level];

#ifdef RTO_CACHE_VARS
  /*
C
C  Cache the radiation fields
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(iL,j),
C$OMP+SHARED(var,varL,ivar,indLG,nTotCells)
      do iL=1,nTotCells
         varL(iL) = var(ivar,indLG(iL))
      enddo
  */
#endif

  if(work)
    {

#pragma omp parallel for default(none), private(iL,nb,abcMid,j), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,jac,dx,dWorker,rhs,dd)
      for(iL=0; iL<num_level_cells; iL++)
	{
      
	  nb = neib + rtuStencilSize*iL;

	  rtOtvetMidPointAbsorptionCoefficients(level,iL,info,nb,abc,abcMid);
	  /*
	  //  Make everything dimensionsless
	  */
	  for(j=0; j<num_neighbors; j++) abcMid[j] = abcMid[j]*dx + epsNum;

	  /*
	  // Inverse Jacobian
	  */
	  jac[iL] = gamma/(1.0+gamma*(abc[iL]*dx+dWorker(iL,indL2G,abcMid)));

	  /*
	  // RHS: absorption must act on the initial field only to insure
	  // photon number conservation
	  */
	  rhs[iL] = dx*rhs[iL] - (1-beta)*0.5*(1.0-exp(-2*abc[iL]*dx))*var(iL);
	  //dd[iL] = rhs[iL];
	}
    }

  for(it=0; it<nit; it++)
    {
      if(work)
	{

	  /*
	  // Compute Laplacian term
	  */
#pragma omp parallel for default(none), private(iL,nb,abcMid,j), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,rhs,jac,dd,dx,fWorker)
	  for(iL=0; iL<num_level_cells; iL++)
	    {

	      nb = neib + rtuStencilSize*iL;

	      rtOtvetMidPointAbsorptionCoefficients(level,iL,info,nb,abc,abcMid);
	      /*
	      //  Make everything dimensionsless
	      */
	      for(j=0; j<num_neighbors; j++) abcMid[j] = abcMid[j]*dx + epsNum;

	      /*
	      // Main operator element
	      */
	      dd[iL] = fWorker(ivar,iL,indL2G,nb,abcMid) - beta*0.5*(1.0-exp(-2*abc[iL]*dx))*var(iL) + rhs[iL];
	    }

	  /*
	  // Update the radiation field
	  */
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivar,dd,jac)
	  for(iL=0; iL<num_level_cells; iL++)
	    {
	      var(iL) += alpha*jac[iL]*dd[iL];
	      if(var(iL) < 0.0) var(iL) = 0;
	    }
	}
      update_buffer_level(level,&ivar,1);
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



/*
// Flux-conserving diffusion coefficients a-la Oran & Boris
*/
void rtOtvetMidPointAbsorptionCoefficients(int level, int iL, int *info, int *nb, float *abc, float *abcMid)
{
  int j;

  if(level == min_level)
    {
      for(j=0; j<num_neighbors; j++)
	{
	  abcMid[j] = 0.5*(abc[iL]+abc[nb[j]]);
	}
    }
  else
    {
      for(j=0; j<num_neighbors; j++)
	{
	  if(info[iL] & (1<<j))
	    {
	      abcMid[j] = (1.0/3.0)*abc[iL] + (2.0/3.0)*abc[nb[j]];
	    }
	  else
	    {
	      abcMid[j] = 0.5*(abc[iL]+abc[nb[j]]);
	    }
	}
    }
}


/*
// Compute the Laplacian and its diagonal element in several ways
*/
float rtOtvetLaplacian_UnitaryTensor_CrossStencil_Full(int ivar, int iL, int *indL2G, int *nb, float *abcLoc)
{
  int j;
  float s, vmax;

#ifndef RT_DEBUG_BLOCK_MASKING
  /*
  // Compute update mask
  */
  vmax = fabs(var(nb[0]));
  for(j=1; j<num_neighbors; j++)
    {
      vmax = (vmax > fabs(var(nb[j]))) ? vmax : fabs(var(nb[j]));
    }
  if(vmax > 1.0e-35)
#endif
    {
      s = 0.0;
      for(j=0; j<num_neighbors; j++) s += (var(nb[j])-var(iL))/abcLoc[j];
      return s/nDim;
    }
#ifndef RT_DEBUG_BLOCK_MASKING
  else return 0.0;
#endif
}


float rtOtvetLaplacian_UnitaryTensor_CrossStencil_Diag(int iL, int *indL2G, float *abcLoc)
{
  int j;
  float s;

  s = 0.0;
  for(j=0; j<num_neighbors; j++) s += 1.0/abcLoc[j];
  return s/nDim;
}


float rtOtvetLaplacian_GenericTensor_SplitStencil_Full(int ivar, int iL, int *indL2G, int *nb, float *abcLoc)
{
  int j;
  float q, vmax;
  float flux[num_neighbors];

#ifndef RT_DEBUG_BLOCK_MASKING
  /*
  // Compute update mask
  */
  vmax = fabs(var(nb[0]));
  for(j=1; j<rtuStencilSize; j++)
    {
      vmax = (vmax > fabs(var(nb[j]))) ? vmax : fabs(var(nb[j]));
    }
  if(vmax > 1.0e-35)
#endif
    {
      /*
      // Co-axial elements
      */
      flux[0] = RPT(0,nb[0]) - RPT(0,iL);
      flux[1] = RPT(0,nb[1]) - RPT(0,iL);
#if (nDim > 1)
      flux[2] = RPT(2,nb[2]) - RPT(2,iL);
      flux[3] = RPT(2,nb[3]) - RPT(2,iL);
#if (nDim > 2)
      flux[4] = RPT(5,nb[4]) - RPT(5,iL);
      flux[5] = RPT(5,nb[5]) - RPT(5,iL);
#endif /* nDim > 2 */
#endif /* nDim > 1 */

      /*
      // Off-axial elements
      */
#if (nDim > 1)
      flux[0] -= 0.25*(RPT(1,nb[ 8])+RPT(1,nb[ 3])-RPT(1,nb[ 6])-RPT(1,nb[ 2]));
      flux[1] += 0.25*(RPT(1,nb[ 9])+RPT(1,nb[ 3])-RPT(1,nb[ 7])-RPT(1,nb[ 2]));
      flux[2] -= 0.25*(RPT(1,nb[ 7])+RPT(1,nb[ 1])-RPT(1,nb[ 6])-RPT(1,nb[ 0]));
      flux[3] += 0.25*(RPT(1,nb[ 9])+RPT(1,nb[ 1])-RPT(1,nb[ 8])-RPT(1,nb[ 0]));
#if (nDim > 2)
      flux[0] -= 0.25*(RPT(3,nb[14])+RPT(3,nb[ 5])-RPT(3,nb[10])-RPT(3,nb[ 4]));
      flux[1] += 0.25*(RPT(3,nb[15])+RPT(3,nb[ 5])-RPT(3,nb[11])-RPT(3,nb[ 4]));
      flux[2] -= 0.25*(RPT(4,nb[16])+RPT(4,nb[ 5])-RPT(4,nb[12])-RPT(4,nb[ 4]));
      flux[3] += 0.25*(RPT(4,nb[17])+RPT(4,nb[ 5])-RPT(4,nb[13])-RPT(4,nb[ 4]));
      flux[4] -= 0.25*(RPT(3,nb[11])+RPT(3,nb[ 1])-RPT(3,nb[10])-RPT(3,nb[ 0])+RPT(4,nb[13])+RPT(4,nb[ 3])-RPT(4,nb[12])-RPT(4,nb[ 2]));
      flux[5] += 0.25*(RPT(3,nb[15])+RPT(3,nb[ 1])-RPT(3,nb[14])-RPT(3,nb[ 0])+RPT(4,nb[17])+RPT(4,nb[ 3])-RPT(4,nb[16])-RPT(4,nb[ 2]));
#endif /* nDim > 2 */
#endif /* nDim > 1 */

      /*
      // Flux limiter OFF (breaks down spherical symmetry!!!)
      */
      /*
      for(j=0; j<num_neighbors; j++)
        {
          q = 0.5*abcLoc[j]*(var(iL)+var(nb[j]));
          if(flux[j] > q)
            {
              flux[j] = q;
            }
          else if(flux[j] < -q)
            {
              flux[j] = -q;
            }
        }
      */

      q = 0.0;
      for(j=0; j<num_neighbors; j++)
        {
	  q += flux[j]/abcLoc[j];
        }
      return q;

    }
#ifndef RT_DEBUG_BLOCK_MASKING
  else return 0.0;
#endif
}


float rtOtvetLaplacian_GenericTensor_SplitStencil_Diag(int iL, int *indL2G, float *abcLoc)
{
  return varET(0,iL)*(1.0/abcLoc[0]+1.0/abcLoc[1])
#if (nDim > 1)
    + varET(2,iL)*(1.0/abcLoc[2]+1.0/abcLoc[3])
    //+ 1.0*fabs(varET(1,iL))*(1.0/abcLoc[0]+1.0/abcLoc[1]+1.0/abcLoc[2]+1.0/abcLoc[3])
#if (nDim > 2)
    + varET(5,iL)*(1.0/abcLoc[4]+1.0/abcLoc[5])
    //+ 1.0*fabs(varET(3,iL))*(1.0/abcLoc[0]+1.0/abcLoc[1]+1.0/abcLoc[4]+1.0/abcLoc[5])
    //+ 1.0*fabs(varET(4,iL))*(1.0/abcLoc[2]+1.0/abcLoc[3]+1.0/abcLoc[4]+1.0/abcLoc[5])
#endif /* nDim > 2 */
#endif /* nDim > 1 */
    ;
}

#endif /* RADIATIVE_TRANSFER && RT_TRANSFER && (RT_TRANSFER_METHOD == RT_METHOD_OTVET) */
