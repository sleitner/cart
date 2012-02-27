#include "config.h"
#ifdef RADIATIVE_TRANSFER

#ifdef _OPENMP
#include <omp.h>
#endif /* _OPENMP */

#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "rt.h"
#include "rt_debug.h"
#include "rt_global.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "F/frt_c.h"

#include "rt_transfer_step.h"
#include "step.h"


#ifndef RT_DEBUG
int rt_debug = 0;
#endif


void rtPackCellData(int level, int cell, frt_real *var, frt_real **p_rawrf);
void rtUnPackCellData(int level, int cell, frt_real *var, frt_real *rawrf);


#ifdef BLASTWAVE_FEEDBACK
extern double blastwave_time_floor;
extern double blastwave_time_cut;
extern double blastwave_time_code_max;
#endif

/*
//  Applies cooling to all cells of a given level
*/
void rtApplyCooling(int level, int num_level_cells, int *level_cells)
{
  int i, cell, nchunk;
#ifdef BLASTWAVE_FEEDBACK
  double blastwave_time;
#endif /* BLASTWAVE_FEEDBACK */
  /* 
  //  Specify types for the Fortran interface
  */
  DEFINE_FRT_INTEFACE(var,rawrf);
  frt_real time;
  frt_intg info;

  /*
  //  The following may be a waste of time if the types are consistent, but compiler should take care of that 
  */
  time = dtl[level];

  /* 
  //  OpenMP chunk size 
  */
  nchunk = 128/(1<<level);
  if(nchunk < 1) nchunk = 1;

  /*
  //  Main loop
  */
#ifdef BLASTWAVE_FEEDBACK
#pragma omp parallel for default(none), private(i,cell,var,rawrf,rawrfBuffer,info), shared(cell_vars,num_level_cells,cell_child_oct,level_cells,level,time,rt_debug,nchunk), schedule(dynamic,nchunk), private(blastwave_time), shared(blastwave_time_cut,blastwave_time_floor,units,constants)
#else      
#pragma omp parallel for default(none), private(i,cell,var,rawrf,rawrfBuffer,info), shared(cell_vars,num_level_cells,cell_child_oct,level_cells,level,time,rt_debug,nchunk), schedule(dynamic,nchunk)
#endif /* BLASTWAVE_FEEDBACK */
  for(i=0; i<num_level_cells; i++) if(cell_is_leaf((cell = level_cells[i])) && cell_gas_density(cell) > 0.0)  /* neg. density means a blow-up, let the code die gracefully in hydro_magic, not here */
    {
      rtPackCellData(level,cell,var,&rawrf);

#ifdef RT_DEBUG
      if(rt_debug.Mode!=0 && cell==cell_find_position(rt_debug.Pos))
	{
	  var[FRT_Debug] = rt_debug.Mode;
	  cart_debug("In cell-level debug for cell %d, level %d",cell,cell_level(cell));
	}
      else
	{
	  var[FRT_Debug] = 0.0;
	}
#endif

#ifdef BLASTWAVE_FEEDBACK
      blastwave_time = cell_blastwave_time(cell)/cell_gas_density(cell);
      if(blastwave_time > blastwave_time_cut)
	{
	  var[FRT_CoolingSuppressionFactor] = 0.0;

	  blastwave_time -= dtl[level]*units->time/constants->yr;
	  if(blastwave_time < blastwave_time_cut) blastwave_time = blastwave_time_floor;
	  cell_blastwave_time(cell) = cell_gas_density(cell)*blastwave_time;
	}
#endif /* BLASTWAVE_FEEDBACK */
      
      /*
      //  Call the Fortran worker 
      */
      frtCall(cooloff)(var,rawrf,&time,&info);

      rtUnPackCellData(level,cell,var,rawrf);
    }
}


void rtStepBegin()
{
  rtInitStep(dtl[min_level]);

#ifdef RT_TRANSFER
  rtStepBeginTransfer();
#endif
}


void rtStepEnd()
{
  frt_real dt = dtl[min_level];
  frt_real vol = num_root_cells;
  frt_real par[3];
  int i;
  struct rtGlobalValue tmp[4];
  MESH_RUN_DECLARE(level,cell);
  double sumSrc, sumRho, sumRhoHI, sumRhoH2;

  /*
  //  Compute global parameters
  */
  start_time(WORK_TIMER);

  for(i=0; i<4; i++) rtGlobalValueInit(&tmp[i],0.0);

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);

  /*
  //    Averages over mesh variables and their derivatives
  */
  sumSrc = sumRho = sumRhoHI = sumRhoH2 = 0.0;

#pragma omp parallel for default(none), private(_Index,cell), shared(_Num_level_cells,_Level_cells,level,cell_vars,cell_child_oct), reduction(+:sumSrc,sumRho,sumRhoHI,sumRhoH2)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell))
    {
#ifdef RT_VAR_SOURCE
      sumSrc += cell_var(cell,RT_VAR_SOURCE)*cell_volume[level]/num_root_cells;
#endif
#ifdef RT_CHEMISTRY
      sumRho += cell_gas_density(cell)*cell_volume[level]/num_root_cells;
      sumRhoHI += cell_HI_density(cell)*cell_volume[level]/num_root_cells;
      sumRhoH2 += cell_H2_density(cell)*cell_volume[level]/num_root_cells;
#endif
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

  rtGlobalValueUpdate(&tmp[0],level,sumSrc);
  rtGlobalValueUpdate(&tmp[1],level,sumRho);
  rtGlobalValueUpdate(&tmp[2],level,sumRhoHI);
  rtGlobalValueUpdate(&tmp[3],level,sumRhoH2);

  MESH_RUN_OVER_LEVELS_END;

  end_time(WORK_TIMER);

  for(i=0; i<4; i++) rtGlobalValueCommunicate(&tmp[i],MPI_SUM,mpi.comm.run);

  start_time(WORK_TIMER);

  par[0] = tmp[0].Value;
#ifdef RT_CHEMISTRY
  par[1] = tmp[2].Value/tmp[1].Value;
  par[2] = tmp[3].Value/tmp[1].Value;
#else
  par[1] = par[2] = 0.0;
#endif

  end_time(WORK_TIMER);

  /*
  //  End step in the reverse order of its beginning
  */
#ifdef RT_TRANSFER
  rtStepEndTransfer();
#endif /* RT_TRANSFER */

  start_time(WORK_TIMER);
  frtCall(stepend)(&dt,&vol,par);
  end_time(WORK_TIMER);

#ifdef RT_DEBUG
  switch(rt_debug.Mode)
    {
    case  1:
    case -1:
      {
	int i, cell;
	cell = cell_find_position(rt_debug.Pos);
	cart_debug("In cell-level debug for cell %d/%d",cell,cell_level(cell));
	for(i=0; i<num_vars; i++)
	  {
	    cart_debug("Var[%d] = %g",i,cell_var(cell,i));
	  }
	break;
      }
    }
  if(rt_debug.Mode < 0)
    {
      cart_error("Aborting on request...");
    }
#endif
}


void rtLevelUpdate(int level)
{
  start_time(RT_LEVEL_UPDATE_TIMER);

#ifdef RT_TRANSFER
  rtLevelUpdateTransfer(level);
#endif

  end_time(RT_LEVEL_UPDATE_TIMER);
}

#ifdef RT_TEST
extern int step;
#endif

/*
//  Used in testing
*/
void rtModifyTimeStep(double *dt)
{
#ifdef RT_TEST
  frt_real rdt;
  frt_intg ins;

  rdt = *dt;
  ins = step;

  frtCall(testmodifytimestep)(&rdt,&ins);

  *dt = rdt;
#endif
}

#endif /* RADIATIVE_TRANSFER */

