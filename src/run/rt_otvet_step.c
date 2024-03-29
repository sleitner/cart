#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

#include <math.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "parallel.h"
#include "rt.h"
#include "rt_global.h"
#include "rt_otvet.h"
#include "rt_transfer.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "frt/frt_c.h"

#include "step.h"


extern int rtOtvetMaxNumIter;


#ifdef RT_OTVET_SAVE_FLUX
#ifdef RT_UV
int rt_flux_field = rt_num_freqs - 1;  /* Local UV by default */
#else
int rt_flux_field = 0;                 /* Local HI by default */
#endif /* RT_UV */
float rt_flux[num_cells][num_neighbors];
#endif /* RT_OTVET_SAVE_FLUX */


#ifdef RT_EXTERNAL_BACKGROUND

extern struct rtGlobalValue rtAvgRF[];
extern struct rtGlobalValue rtAvgAC[];
extern struct rtGlobalValue rtAvgACxRF[];

const float rtConvFac = 1.0;
const float rtFMaxFac = 3.0;
const float rtFMinFac = 0.0e-6;

float rtBarF[rt_num_freqs];
float rtBarK[rt_num_freqs];

#endif /* RT_EXTERNAL_BACKGROUND */


DEFINE_LEVEL_ARRAY(float,BufferFactor);


float *cache_var;
float *cache_varET;


extern int rt_limit_signal_speed_to_c;
extern int rtNumOtvetETVars;
extern int rtOtvetETVars[];
extern int rtOtvetOTBox[];
extern float rtGlobalAC[];


void rtOtvetMidPointAbsorptionCoefficients(int level, int iL, int *info, int *nb, float *abc, float *abcMid);

typedef struct 
{
  float (*Diag)(int iL, int *indL2G, float *abcLoc);
  float (*Full)(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *flux);
}
rt_laplacian_t;


void rtOtvetSolveFieldEquation(int ivar, int level, int num_level_cells, int num_total_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd2, float *dfx, int nit, int work, rt_laplacian_t lap);

extern rt_laplacian_t rt_unitary;
extern rt_laplacian_t rt_generic;


void rtStepBeginTransferOtvet(struct rtGlobalValue *maxAC)
{
  const float tauMin = 1.0e-2;
  int level, freq;

  for(level=min_level; level<=max_level; level++)
    {
      if(BufferFactor[level] < 1.0) BufferFactor[level] = 1.0;
    }

  for(freq=0; freq<rt_num_freqs; freq++)
    {
      rtOtvetOTBox[freq] = (maxAC[freq].Value*num_grid < tauMin);

#ifdef RT_EXTERNAL_BACKGROUND
      /*
      //  This is so that these values are the same for all nodes.
      */
      rtBarF[freq] = rtAvgRF[rt_num_freqs+freq].Value;
      if(rtBarF[freq] > 0.0)
	{
	  rtBarK[freq] = rtAvgACxRF[rt_num_freqs+freq].Value/rtBarF[freq];
	}
      else
	{
	  rtBarK[freq] = rtAvgAC[rt_num_freqs+freq].Value;
	}

#ifdef RT_OUTPUT
      cart_debug("RT: OTVET Far Field %d, <f> = %10.3le, <k> = %10.3le",freq,rtBarF[freq],rtBarK[freq]);
#endif

#endif /* RT_EXTERNAL_BACKGROUND */
    }
}


/* 
// Helper macros
*/
#define varL(iL)    cell_var(indL2G[iL],ivarL)
#define varG(iL)    cell_var(indL2G[iL],ivarG)

#define otfL(iL)    cell_var(indL2G[iL],RT_VAR_OT_FIELD)
#define otfG(iL)    1.0


void rtLevelUpdateTransferOtvet(int level)
{
  /*
  // Extra arrays
  */
  int *indL2G, *neib, *info, *tmp;
  float *abc[2], *dd2, *jac, *rhs, *dfx;

  /*
  // Work variables
  */
  int work, num_hashed;
  int num_level_cells, num_all_cells, num_total_cells, *nb;
  int iL, iG, j, offset, index;
  int nit, freq, ivarL, ivarG, nit0;
  float xiUnit;
  int nvars, vars[rt_num_fields];

  start_time(WORK_TIMER);

  /*
  // If there are no local cells, we do a dry run only
  */
  work = (num_cells_per_level[level] > 0);
  if(work)
    {
      /* 
      //  Allocate memory for index arrays
      */
#ifdef RT_OTVET_NO_GLOBAL_ARRAY
      int *indG2L = cart_alloc(int,num_cells);
#else
      static int indG2L[num_cells];
#endif

      /*
      //  Find all leaves (if we solve levels separately, we can never insure that 
      //  light fronts on both levels propagate at the same speed, even if I-fronts do.)
      */
      select_level(level,CELL_TYPE_LOCAL | CELL_TYPE_LEAF,&num_level_cells,&info);
      /*
      // Allocate the rest of arrays
      */
      neib = cart_alloc(int, (size_t)rtStencilSize*num_level_cells );
      num_all_cells = 100 + (int)(BufferFactor[level]*num_level_cells);
      indL2G = cart_alloc(int, num_all_cells );
      linear_array_copy_int(indL2G,info,num_level_cells);

      /*
      // Find the minimum size of the hash table
      */
      if(num_level_cells > 0)
        {
          linear_array_max_int(num_level_cells,indL2G,&num_hashed);
          num_hashed++;
        }
      else num_hashed = 0;

      /*
      // Initialize the global-to-local index array
      */
#pragma omp parallel for default(none), private(iG), shared(indG2L,size_cell_array,num_hashed)
      for(iG=0; iG<num_hashed; iG++)
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
	  offset = rtStencilSize*iL;
	  
	  rtGetStencil(level,indL2G[iL],neib+offset);
	  
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
	  nb = neib + rtStencilSize*iL;
	  
	  for(j=0; j<rtStencilSize; j++)
	    {
	      if(nb[j] < num_hashed)
		{
		  index = indG2L[nb[j]];
		}
	      else
		{
		  /*
		  // Initialize the next segment of the global-to-local index array
		  */
#pragma omp parallel for default(none), private(iG), shared(indG2L,size_cell_array,num_hashed,nb,j)
		  for(iG=num_hashed; iG<=nb[j]; iG++)
		    {
		      indG2L[iG] = -1;
		    }
		  num_hashed = nb[j] + 1;
		  index = -1;

		}

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
		      linear_array_copy_int(tmp,indL2G,nit);
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

#ifdef RT_OTVET_NO_GLOBAL_ARRAY
      /*
      //  It is not needed any more
      */
      cart_free(indG2L);
#endif

      rhs = cart_alloc(float, num_level_cells );
      jac = cart_alloc(float, num_level_cells );
      dd2 = cart_alloc(float, num_level_cells );
      dfx = cart_alloc(float, num_level_cells );

      /*
      // Are we using too much memory?
      */
      if(num_total_cells<0.5*num_all_cells && BufferFactor[level]>2.0/0.8)
	{
	  BufferFactor[level] *= 0.8;
	}
      
#ifdef RT_OTVET_CACHE_RF
      cache_var = cart_alloc(float, num_total_cells );
#endif

#ifdef RT_OTVET_CACHE_ET
      cache_varET = cart_alloc(float, 6*num_total_cells );

#pragma omp parallel for default(none), private(iL,j), shared(num_total_cells,cell_vars,indL2G,cache_varET)
      for(iL=0; iL<num_total_cells; iL++)
	{
	  for(j=0; j<6; j++)
	    {
	      cache_varET[j+6*iL] = cell_var(indL2G[iL],rt_et_offset+j);
	    }
	}
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

  end_time(WORK_TIMER);

  /*
  // Iterate over all non-zero frequencies, two at a time (local and global fields)
  */
  for(freq=0; freq<rt_num_freqs; freq++)
    {

      start_time(WORK_TIMER);

      ivarL = rt_field_offset + freq;
      ivarG = ivarL + rt_num_freqs;

      /*
      // Compute the absorption coefficient at this level
      */
      rtComputeAbsLevel(level,num_total_cells,indL2G,freq,abc);

      end_time(WORK_TIMER);

      /*
      // Is the box optically thin?
      */
      if(rtOtvetOTBox[freq])
	{
	  
	  start_time(WORK_TIMER);

	  if(work)
	    {
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  varL(iL) = otfL(iL);
		  varG(iL) = otfG(iL);
		}
	    }

	  end_time(WORK_TIMER);

	}
      else
	{

	  start_time(WORK_TIMER);

	  /*
	  //  Number of iterations
	  */
	  nit = rtOtvetMaxNumIter;

	  /*
	  // Number of iterations needed for the signal propagation speed to less or equal c.
	  */
	  if(rt_limit_signal_speed_to_c)
	    {
	      xiUnit = units->length/(constants->c*units->time);
	      nit0 = 1 + (int)(dtl[level]/xiUnit/cell_size[level]);
	      if(nit > nit0)
		{
		  nit = nit0;
		}
	    }

	  /*
	  // Update local field
	  */
	  if(work)
	    {
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,level,rhs)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  rhs[iL] = cell_rt_source(indL2G[iL]);
		}

#ifdef RT_OTVET_CACHE_RF
#pragma omp parallel for default(none), private(iL), shared(num_total_cells,cell_vars,indL2G,ivarL,cache_var)
	      for(iL=0; iL<num_total_cells; iL++)
		{
		  cache_var[iL] = varL(iL);
		}
#endif /* RT_OTVET_CACHE_RF */
	    }
	  end_time(WORK_TIMER);

	  rtOtvetSolveFieldEquation(ivarL,level,num_level_cells,num_total_cells,indL2G,neib,info,abc[0],rhs,jac,dd2,dfx,nit,work,rt_generic);

	  if(work)
	    {
	      start_time(WORK_TIMER);

#ifdef RT_OTVET_CACHE_RF
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,cache_var)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  if(cache_var[iL] < 0.0) cache_var[iL] = 0.0;
		  if(cache_var[iL] > otfL(iL)) cache_var[iL] = otfL(iL);
		  varL(iL) = cache_var[iL];
		}
#else  /* RT_OTVET_CACHE_RF */
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  if(varL(iL) < 0.0) varL(iL) = 0.0;
		  if(varL(iL) > otfL(iL)) varL(iL) = otfL(iL);
		}
#endif /* RT_OTVET_CACHE_RF */

	      end_time(WORK_TIMER);
	    }

#ifdef RT_EXTERNAL_BACKGROUND
	  /*
	  // Update global field
	  */
	  if(work)
	    {
	      frt_intg idx;
	      float QTerm, CFac;

	      start_time(WORK_TIMER);

	      if(rtBarF[freq] > 0.0)
		{
		  CFac = rtBarK[freq]/pow(rtBarF[freq],rtConvFac);
		}
	      else
		{
		  CFac = rtBarK[freq];
		}

	      idx = freq + 1;
              QTerm = frtCall(transferglobalqterm)(&idx);

#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarL,ivarG,level,rhs,abc,QTerm,CFac,freq)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  rhs[iL] = varG(iL)*CFac + QTerm;
		  abc[1][iL] += QTerm;
		}

#ifdef RT_OTVET_CACHE_RF
#pragma omp parallel for default(none), private(iL,j), shared(num_total_cells,cell_vars,indL2G,ivarG,cache_var)
	      for(iL=0; iL<num_total_cells; iL++)
		{
		  cache_var[iL] = varG(iL);
		}
#endif /* RT_OTVET_CACHE_RF */
	      end_time(WORK_TIMER);
	    }

	  rtOtvetSolveFieldEquation(ivarG,level,num_level_cells,num_total_cells,indL2G,neib,info,abc[1],rhs,jac,dd2,dfx,nit,work,rt_unitary);

	  if(work)
	    {
	      start_time(WORK_TIMER);

#ifdef RT_OTVET_CACHE_RF
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarG,cache_var)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  if(cache_var[iL] < 0.0) cache_var[iL] = 0.0;
		  if(cache_var[iL] > rtFMaxFac*otfG(iL)) cache_var[iL] = rtFMaxFac*otfG(iL);
		  varG(iL) = cache_var[iL];
		}
#else  /* RT_OTVET_CACHE_RF */
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarG)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  if(varG(iL) < rtFMinFac*otfG(iL)) varG(iL) = rtFMinFac*otfG(iL);
		  if(varG(iL) > rtFMaxFac*otfG(iL)) varG(iL) = rtFMaxFac*otfG(iL);
		}
#endif /* RT_OTVET_CACHE_RF */
	      end_time(WORK_TIMER);
	    }

#else  /* RT_EXTERNAL_BACKGROUND */
	  /*
	  // Set global field to zero
	  */
	  if(work)
	    {
	      start_time(WORK_TIMER);

#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivarG)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  varG(iL) = 0.0;
		}

	      end_time(WORK_TIMER);
	    }
#endif /* RT_EXTERNAL_BACKGROUND */
	}
    }

  /*
  //  Last buffer update - can combine for all fields
  */
  nvars = rt_num_freqs;
#ifdef RT_EXTERNAL_BACKGROUND
  nvars *= 2;
#endif /* RT_EXTERNAL_BACKGROUND */
  for(j=0; j<nvars; j++)
    {
      vars[j] = rt_field_offset + j;
    }

  start_time(RT_SOLVE_EQUATION_UPDATE_TIMER);
  update_buffer_level(level,vars,nvars);
  end_time(RT_SOLVE_EQUATION_UPDATE_TIMER);

  start_time(WORK_TIMER);

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
      cart_free(dd2);
      cart_free(dfx);
#ifdef RT_OTVET_CACHE_RF
      cart_free(cache_var);
#endif
#ifdef RT_OTVET_CACHE_ET
      cart_free(cache_varET);
#endif
    }

  end_time(WORK_TIMER);

}


/* 
// Helper macros
*/
#ifdef RT_OTVET_CACHE_RF
#define var(iL)          cache_var[iL]
#else
#define var(iL)          cell_var(indL2G[iL],ivar)
#endif

#ifdef RT_OTVET_CACHE_ET
#define varET(ind,iL)    cache_varET[ind+6*iL]
#else
#define varET(ind,iL)    cell_var(indL2G[iL],rt_et_offset+ind)
#endif

#define RPT(ind,iL)      (varET(ind,iL)*var(iL))


void rtOtvetIterate1(int it, int nit, int ivar, int level, int num_level_cells, int num_total_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd2, float *dfx, rt_laplacian_t lap);
void rtOtvetIterate2(int it, int nit, int ivar, int level, int num_level_cells, int num_total_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd2, float *dfx, rt_laplacian_t lap);
void rtOtvetIterate3(int it, int nit, int ivar, int level, int num_level_cells, int num_total_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd2, float *dfx, rt_laplacian_t lap);


void rtOtvetSolveFieldEquation(int ivar, int level, int num_level_cells, int num_total_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd2, float *dfx, int nit, int work, rt_laplacian_t lap)
{
  int it;
  
  for(it=1; it<=nit; it++)
    {
      if(work)
	{
	  start_time(WORK_TIMER);

	  rtOtvetIterate1(it,nit,ivar,level,num_level_cells,num_total_cells,indL2G,neib,info,abc,rhs,jac,dd2,dfx,lap);

	  end_time(WORK_TIMER);
	}

      /*
      //  The last update done simultaneously for all fields
      */
      if(it < nit)
	{
#ifdef RT_OTVET_CACHE_RF
	  int iL;

	  /*
	  //  Update level cells from the cache for the buffer update
	  */
	  if(work)
	    {
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivar,cache_var)
	      for(iL=0; iL<num_level_cells; iL++)
		{
		  cell_var(indL2G[iL],ivar) = cache_var[iL];
		}
	    }
#endif /* RT_OTVET_CACHE_RF */

	  start_time(RT_SOLVE_EQUATION_UPDATE_TIMER);
	  update_buffer_level(level,&ivar,1);
	  end_time(RT_SOLVE_EQUATION_UPDATE_TIMER);

#ifdef RT_OTVET_CACHE_RF
	  /*
	  //  Update cache with buffer cells, no need to update level cells
	  //  as they are already updated.
	  */
	  if(work)
	    {
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,num_total_cells,cell_vars,indL2G,ivar,cache_var)
	      for(iL=num_level_cells; iL<num_total_cells; iL++)
		{
		  cache_var[iL] = cell_var(indL2G[iL],ivar);
		}
	    }
#endif /* RT_OTVET_CACHE_RF */
	}
    }
}


void rtOtvetIterate1(int it, int nit, int ivar, int level, int num_level_cells, int num_total_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd2, float *dfx, rt_laplacian_t lap)
{
  /*
  // Numerical parameters
  // alpha should be: 
  //   (1) sufficiently smaller than 1 for stability
  //   (2) not too small for efficiency
  */
  const float CFL = 0.8;
  const float gamma = 1.0;
  const float epsNum = 1.0e-3;
  /*
  //  NG: I am not sure why this parameter is needed, but it is found empirically to work well. 
  */
  const float eta = 0.1;

  int j, iL, *nb;
  float abcMid[num_neighbors], flux[num_neighbors];
  float dx = cell_size[level];

  /*
  // Compute Laplacian term
  */
#ifdef RT_OTVET_SAVE_FLUX
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,flux), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,rhs,jac,dd2,dx,lap,cache_var,it,rt_flux_field,rt_flux,nit)
#else
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,flux), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,rhs,jac,dd2,dx,lap,cache_var,it)
#endif /* RT_OTVET_SAVE_FLUX */
  for(iL=0; iL<num_level_cells; iL++)
    {

      nb = neib + rtStencilSize*iL;

      rtOtvetMidPointAbsorptionCoefficients(level,iL,info,nb,abc,abcMid);
      /*
      //  Make everything dimensionsless
      */
      for(j=0; j<num_neighbors; j++) abcMid[j] = abcMid[j]*dx + epsNum;

      if(it == 1)
	{
	  /*
	  // Jacobian
	  */
	  jac[iL] = gamma/(1.0+gamma*(eta*abc[iL]*dx+lap.Diag(iL,indL2G,abcMid)));

	  /*
	  // RHS: absorption must act on the initial field only to insure
	  // photon number conservation
	  */
	  rhs[iL] = dx*rhs[iL] - (1-eta)*abc[iL]*dx*var(iL);
	}

      /*
      // Main operator element
      */
      dd2[iL] = lap.Full(ivar,iL,indL2G,nb,abcMid,flux) - eta*abc[iL]*dx*var(iL) + rhs[iL];

#ifdef RT_OTVET_SAVE_FLUX
      /* save flux at the last iteration */
      if(it==nit && rt_field_offset+rt_flux_field==ivar)
	{
	  for(j=0; j<num_neighbors; j++) rt_flux[indL2G[iL]][j] = -flux[j];
	}
#endif /* RT_OTVET_SAVE_FLUX */
    }

  /*
  // Update the radiation field
  */
#pragma omp parallel for default(none), private(iL), shared(num_level_cells,cell_vars,indL2G,ivar,dd2,jac,cache_var)
  for(iL=0; iL<num_level_cells; iL++)
    {
      var(iL) += CFL*jac[iL]*dd2[iL];
      if(var(iL) < 0.0) var(iL) = 0;
    }
}


void rtOtvetIterate2(int it, int nit, int ivar, int level, int num_level_cells, int num_total_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd2, float *dfx, rt_laplacian_t lap)
{
  /*
  // Numerical parameters:
  // gamma <= 2/(2alpha-1)
  // linear stability analysis is ignorant of the value of beta
  */
  const float CFL = 0.8;
  const float alpha = 1;
  const float beta  = alpha;
  const float gamma = 1;
  const float epsNum = 1.0e-3;
  const float eta = 0.15;

  int j, iL, *nb;
  float abcMid[num_neighbors], flux[num_neighbors];
  float dx = cell_size[level];
  float abcCen, s1, s2, b1, b2, det;

  /*
  // Compute Laplacian term
  */
#ifdef RT_OTVET_SAVE_FLUX
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,flux), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,jac,rhs,dd2,dfx,dx,lap,cache_var,it,rt_flux_field,rt_flux,nit)
#else
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,flux), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,jac,rhs,dd2,dfx,dx,lap,cache_var,it)
#endif /* RT_OTVET_SAVE_FLUX */
  for(iL=0; iL<num_level_cells; iL++)
    {

      nb = neib + rtStencilSize*iL;

      rtOtvetMidPointAbsorptionCoefficients(level,iL,info,nb,abc,abcMid);

      /*
      //  Make everything dimensionsless
      */
      for(j=0; j<num_neighbors; j++)
	{
	  abcMid[j] = abcMid[j]*dx + epsNum;
	}

      /*
      // Main operator element
      */
      dd2[iL] = lap.Full(ivar,iL,indL2G,nb,abcMid,flux);

      if(it == 1)
	{
	  /*
	  // Jacobian
	  */
	  jac[iL] = lap.Diag(iL,indL2G,abcMid);

	  /*
	  // Initial value for the flux divergence
	  */
	  dfx[iL] = dd2[iL];

	  /*
	  // RHS: absorption must act on the initial field only to insure
	  // photon number conservation
	  */
	  rhs[iL] = dx*rhs[iL] - (1-eta)*abc[iL]*dx*var(iL);
	}

#ifdef RT_OTVET_SAVE_FLUX
      /* save flux at the last iteration */
      if(it==1 && rt_field_offset+rt_flux_field==ivar)
	{
	  for(j=0; j<num_neighbors; j++) rt_flux[indL2G[iL]][j] = -flux[j];
	}
#endif /* RT_OTVET_SAVE_FLUX */
    }

  /*
  // Update the radiation field
  */
#ifdef OPENMP_DECLARE_CONST
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,abcCen,det,s1,s2,b1,b2), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,dx,dd2,jac,dfx,rhs,cache_var,eta,epsNum)
#else
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,abcCen,det,s1,s2,b1,b2), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,dx,dd2,jac,dfx,rhs,cache_var)
#endif /* OPENMP_DECLARE_CONST */
  for(iL=0; iL<num_level_cells; iL++)
    {
      nb = neib + rtStencilSize*iL;

      rtOtvetMidPointAbsorptionCoefficients(level,iL,info,nb,abc,abcMid);

      /*
      //  Make everything dimensionsless
      */
      abcCen = 0.0;
      for(j=0; j<num_neighbors; j++)
	{
	  abcMid[j] = abcMid[j]*dx + epsNum;
	  abcCen += abcMid[j];
	}
      abcCen /= num_neighbors;
      abcCen = abc[iL]*dx + epsNum;

      s1 = 1 + gamma*eta*abc[iL]*dx;
      s2 = 1 + gamma*beta*abcCen;
      b1 = dfx[iL] - eta*abc[iL]*dx*var(iL) + rhs[iL];
      b2 = -abcCen*dfx[iL] + abcCen*dd2[iL];

      det = s1*s2 + alpha*gamma*gamma*abcCen*jac[iL];

      var(iL) += CFL*gamma*(s2*b1+alpha*gamma*b2)/det;
      dfx[iL] += CFL*gamma*(s1*b2-gamma*abcCen*jac[iL]*b1)/det;

      if(var(iL) < 0.0) var(iL) = 0;
    }
}


void rtOtvetIterate3(int it, int nit, int ivar, int level, int num_level_cells, int num_total_cells, int *indL2G, int *neib, int *info, float *abc, float *rhs, float *jac, float *dd2, float *dfx, rt_laplacian_t lap)
{
  /*
  // Numerical parameters:
  // gamma <= 2/(2alpha-1)
  // linear stability analysis is ignorant of the value of beta
  */
  const float CFL = 0.8;
  const float alpha = 1;
  const float beta  = alpha;
  const float gamma = 1;
  const float epsNum = 1.0e-3;
  const float eta = 0.15;
  float var1, var2;

  int j, iL, *nb;
  float abcMid[num_neighbors], flux[num_neighbors];
  float dx = cell_size[level];
  float abcCen, s1, s2, b1, b2, det;

  /*
  // Compute Laplacian term
  */
#ifdef RT_OTVET_SAVE_FLUX
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,flux), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,jac,rhs,dd2,dfx,dx,lap,cache_var,it,rt_flux_field,rt_flux,nit)
#else
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,flux), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,jac,rhs,dd2,dfx,dx,lap,cache_var,it)
#endif /* RT_OTVET_SAVE_FLUX */
  for(iL=0; iL<num_level_cells; iL++)
    {

      nb = neib + rtStencilSize*iL;

      rtOtvetMidPointAbsorptionCoefficients(level,iL,info,nb,abc,abcMid);

      /*
      //  Make everything dimensionsless
      */
      for(j=0; j<num_neighbors; j++)
	{
	  abcMid[j] = abcMid[j]*dx + epsNum;
	}

      /*
      // Main operator element
      */
      dd2[iL] = lap.Full(ivar,iL,indL2G,nb,abcMid,flux);

      if(it == 1)
	{
	  /*
	  // Jacobian
	  */
	  jac[iL] = lap.Diag(iL,indL2G,abcMid);

	  /*
	  // Initial value for the flux divergence
	  */
	  dfx[iL] = dd2[iL];

	  /*
	  // RHS: absorption must act on the initial field only to insure
	  // photon number conservation
	  */
	  rhs[iL] = dx*rhs[iL] - (1-eta)*abc[iL]*dx*var(iL);
	}

#ifdef RT_OTVET_SAVE_FLUX
      /* save flux at the last iteration */
      if(it==1 && rt_field_offset+rt_flux_field==ivar)
	{
	  for(j=0; j<num_neighbors; j++) rt_flux[indL2G[iL]][j] = -flux[j];
	}
#endif /* RT_OTVET_SAVE_FLUX */
    }

  /*
  // Update the radiation field
  */
#pragma omp parallel for default(none), private(iL,nb,abcMid,j,abcCen,det,s1,s2,b1,b2,var1,var2), shared(num_level_cells,cell_vars,indL2G,ivar,abc,neib,level,info,dx,dd2,jac,dfx,rhs,cache_var,it)
  for(iL=0; iL<num_level_cells; iL++)
    {
      nb = neib + rtStencilSize*iL;

      rtOtvetMidPointAbsorptionCoefficients(level,iL,info,nb,abc,abcMid);

      /*
      //  Make everything dimensionsless
      */
      abcCen = 0.0;
      for(j=0; j<num_neighbors; j++)
	{
	  abcMid[j] = abcMid[j]*dx + epsNum;
	  abcCen += abcMid[j];
	}
      abcCen /= num_neighbors;
      abcCen = abc[iL]*dx + epsNum;

      var1 = var(iL) + CFL*gamma/(1.0+gamma*(eta*abc[iL]*dx+jac[iL]))*(dd2[iL]-eta*abc[iL]*dx*var(iL)+rhs[iL]);

      s1 = 1 + gamma*eta*abc[iL]*dx;
      s2 = 1 + gamma*beta*abcCen;
      b1 = dfx[iL] - eta*abc[iL]*dx*var(iL) + rhs[iL];
      b2 = -abcCen*dfx[iL] + abcCen*dd2[iL];

      det = s1*s2 + alpha*gamma*gamma*abcCen*jac[iL];

      var2 = var(iL) + CFL*gamma*(s2*b1+alpha*gamma*b2)/det;
      dfx[iL] += CFL*gamma*(s1*b2-gamma*abcCen*jac[iL]*b1)/det;

      var(iL) = 0.95*var1 + 0.05*var2;

      if(var(iL) < 0.0) var(iL) = 0;
    }
}


void rtOtvetSplitUpdate(int level, int num_level_cells, int *level_cells)
{
#ifdef RT_OTVET_SAVE_FLUX
  int i, j, k;
  int cell;
  int children[num_children];
  double new_var;
  const double factor = ((double)(1.0/(1<<nDim)));

  if(level < max_level)
    {
#pragma omp parallel for default(none), private(i,cell,j,k,children,new_var), shared(num_level_cells,level_cells,cell_child_oct,rt_flux)
      for(i=0; i<num_level_cells; i++)
	{
	  cell = level_cells[i];
	  if(cell_is_refined(cell))
	    {
	      /*
	      // Average over children
	      */
	      cell_all_children(cell,children);
	      for(j=0; j<num_neighbors; j++)
		{
		  new_var = 0.0;
		  for(k=0; k<num_children; k++)
		    {
		      new_var += rt_flux[children[k]][j];
		    }
		  rt_flux[cell][j] = new_var*factor; 
		}
	    }
	}
    }
#endif /* RT_OTVET_SAVE_FLUX */
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
float rtOtvet_UnitaryTensorFull(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *flux)
{
  int j;
  float q, vmax;

#ifndef RT_DEBUG_BLOCK_MASKING
  /*
  // Compute update mask
  */
  vmax = fabs(var(iL));
  for(j=0; j<num_neighbors; j++)
    {
      vmax = (vmax > fabs(var(nb[j]))) ? vmax : fabs(var(nb[j]));
    }
  if(vmax > 1.0e-35)
#endif
    {
      for(j=0; j<num_neighbors; j++) flux[j] = (var(nb[j])-var(iL))/abcLoc[j]/nDim;

      q = 0.0;
      for(j=0; j<num_neighbors; j++) q += flux[j];
      return q;
    }
#ifndef RT_DEBUG_BLOCK_MASKING
  else
    {
#ifdef RT_OTVET_SAVE_FLUX
      for(j=0; j<num_neighbors; j++) flux[j] = 0.0;
#endif /* RT_OTVET_SAVE_FLUX */

      return 0.0;
    }
#endif
}


float rtOtvet_UnitaryTensorDiag(int iL, int *indL2G, float *abcLoc)
{
  int j;
  float s;

  s = 0.0;
  for(j=0; j<num_neighbors; j++) s += 1.0/abcLoc[j];
  return s/nDim;
}


float rtOtvet_GenericTensorFull(int ivar, int iL, int *indL2G, int *nb, float *abcLoc, float *flux)
{
  int j;
  float q, vmax;

#ifndef RT_DEBUG_BLOCK_MASKING
  /*
  // Compute update mask
  */
  vmax = fabs(var(iL));
  for(j=0; j<rtStencilSize; j++)
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

      for(j=0; j<num_neighbors; j++) flux[j] /= abcLoc[j];

      q = 0.0;
      for(j=0; j<num_neighbors; j++) q += flux[j];
      return q;

    }
#ifndef RT_DEBUG_BLOCK_MASKING
  else
    {
#ifdef RT_OTVET_SAVE_FLUX
      for(j=0; j<num_neighbors; j++) flux[j] = 0.0;
#endif /* RT_OTVET_SAVE_FLUX */

      return 0.0;
    }
#endif
}


float rtOtvet_GenericTensorDiag(int iL, int *indL2G, float *abcLoc)
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

rt_laplacian_t rt_unitary = { rtOtvet_UnitaryTensorDiag, rtOtvet_UnitaryTensorFull };
rt_laplacian_t rt_generic = { rtOtvet_GenericTensorDiag, rtOtvet_GenericTensorFull };


#ifdef RT_OTVET_SAVE_FLUX
#include "frt/frt_c.h"

#ifdef __cplusplus
extern "C" {
#endif
void frtCall(getrfunits)(frt_intg *freq, frt_real *uNear, frt_real *uFar);
void frtCall(getrfunits2)(frt_real *var, frt_real *rawrf, frt_intg *freq, frt_real *uNear, frt_real *uFar);
#ifdef __cplusplus
}
#endif


void rtGetRadiationFlux(int cell, float flux[num_neighbors])
{
  int j;
  double fac;
  frt_intg freq = rt_flux_field + 1;
  frt_real uNear, uFar;
  int level = cell_level(cell);
  
#ifdef RT_UV_OLDSTYLE_3x1
  DEFINE_FRT_INTEFACE(var,rawrf);
  rtPackCellData(level,cell,var,&rawrf);
  frtCall(getrfunits2)(var,rawrf,&freq,&uNear,&uFar);
#else
  frtCall(getrfunits)(&freq,&uNear,&uFar);
#endif

  if(rt_flux_field < rt_far_freq_offset)
    {
      fac = constants->c*6.626e-27*uNear;
    }
  else
    {
      fac = constants->c*6.626e-27*uFar;
    }

  for(j=0; j<num_neighbors; j++)
    {
      flux[j] = fac*rt_flux[cell][j];
    }
}
#endif /* RT_OTVET_SAVE_FLUX */


#endif /* RADIATIVE_TRANSFER && RT_TRANSFER && (RT_TRANSFER_METHOD == RT_METHOD_OTVET) */
