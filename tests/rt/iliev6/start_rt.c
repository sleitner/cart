
#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "hydro.h"
#include "iterators.h"
#include "logging.h"
#include "parallel.h"
#include "refinement.h"
#include "rt_utilities.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "extra/healpix.h"
#include "extra/ifrit.h"


const float N50 = 1.0;
const float T_i = 1.0e2;
const int BottomLevel = 2;
const float ExtraScale = num_grid/16.0;


extern float rtSingleSourceVal;
extern double rtSingleSourcePos[nDim];

double tStart;
extern int rtOtvetMaxNumIter;


void units_set_art(double OmegaM, double h, double Lbox);


void refine_level( int cell, int level )
{
  int j;
  float pos[nDim];
  float r;

  cart_assert( cell >= 0 && cell < num_cells );
  cart_assert( cell_level(cell) == level );
	
  cell_position(cell, pos);

  for(j=0; j<nDim; j++)
    {
      pos[j] -= 0.5*num_grid;
    }

  r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);

  if ( r < num_grid/ExtraScale )
    {
      refinement_indicator(cell,0) = 1.0;	
    }
  else
    {
      refinement_indicator(cell,0) = 0.0;
    }
}
	

void rt_initial_conditions( int cell )
{
  int j;
  float pos[nDim];
  float r;
  float rho, dxscal = 0.0915/0.8*num_grid/ExtraScale;

  cell_position(cell, pos);

  for(j=0; j<nDim; j++)
    {
      pos[j] -= 0.5*num_grid;
    }

  r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
                  
  rho = dxscal*dxscal/(1.0e-30+r*r);
  if(rho > 1.0) rho = 1.0;

  cell_gas_density(cell) = rho;
  cell_momentum(cell,0) = 0.0;
  cell_momentum(cell,1) = 0.0;
  cell_momentum(cell,2) = 0.0;
  cell_gas_gamma(cell) = (5.0/3.0);

  cell_gas_internal_energy(cell) = T_i/units->temperature/(constants->gamma-1);

  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)*(constants->gamma-1);
  cell_gas_energy(cell) = cell_gas_internal_energy(cell);

  cell_HI_density(cell) = 1.0*rho;
  cell_HII_density(cell) = 0.0;
  cell_HeI_density(cell) = 1.0-10*rho;
  cell_HeII_density(cell) = 0.0;
  cell_HeIII_density(cell) = 0.0;
  cell_H2_density(cell) = 0.0;
}


void set_rt_initial_conditions(void)
{
  MESH_RUN_DECLARE(level,cell);

  _MaxLevel = BottomLevel;

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  rt_initial_conditions(cell);
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;

  for(level=_MaxLevel; level>=min_level; level--)
    {
      hydro_split_update(level);
    }
}


void run_output()
{
  const int nvars = 8;
  const int nbin1 = 32 * (1 << BottomLevel);
  int varid[] = { I_FRACTION+RT_HVAR_OFFSET+0, HVAR_GAS_DENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC, rt_freq_offset+0, rt_freq_offset+1, rt_freq_offset+2 };
  int nbin[] = { nbin1, nbin1, nbin1 };
  double bb[6];
  int done;
  double tPhys;
  char filename[99];

  if(step>0 && step%2==0) return;

  bb[0] = bb[2] = bb[4] = num_grid*(0.5-1/ExtraScale);
  bb[1] = bb[3] = bb[5] = num_grid*(0.5+1/ExtraScale);
 
  sprintf(filename,"OUT/out.%05d.bin",step);
  ifrit.OutputMesh(filename,max_level,nbin,bb,nvars,varid);

  done = 0;
  if(local_proc_id == MASTER_NODE)
    {

      tPhys = units->time*(tl[0]-tStart)/constants->Myr;
      printf("Output: %d,  Time: %lg\n",step,tPhys);
      if(tPhys > 24.9) done = 1;
    }

  MPI_Bcast(&done,1,MPI_INT,MASTER_NODE,MPI_COMM_WORLD);

  if(done)
    {
      finalize_logging();
      MPI_Finalize();
      exit(0);
    }
}


void init_run()
{
   int i, level;
   int num_level_cells;
   int *level_cells;
   float astart, hubble;
   const float n0 = 3.2;

   /* set units */
   astart = 1;
   hubble = 1;
   units_set_art(n0*pow(astart,3)/(1.123e-5*hubble*hubble),hubble,ExtraScale*0.8e-3/(astart*hubble));

   units_reset();
   units_update(min_level);

   cart_debug("in init");

   /* build buffer */
   build_cell_buffer();
   cart_debug("built cell buffer");
   repair_neighbors();

   check_map();

   cart_debug("repaired neighbors");

   set_rt_initial_conditions();

   /* do initial refinements */
   for ( level = min_level; level < BottomLevel; level++ )
     {
       cart_debug("refining level %u", level );

       select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
       cart_debug("num_level_cells = %u", num_level_cells );
       for ( i = 0; i < num_level_cells; i++ )
	 {
	   refine_level( level_cells[i], level );
	 }
       cart_free( level_cells );
       cart_debug("about to refine level %u", level );
       refine(level);
     }

   cart_debug("setting initial conditions");
   set_rt_initial_conditions();

#ifdef HYDRO_TRACERS
   cart_debug("setting hydro tracers");
   set_hydro_tracers( min_level+1 );
#endif /* HYDRO_TRACERS */

   cart_debug("set initial conditions");

   for ( level = min_level; level <= max_level; level++ )
     {
       cart_debug("updating level %u", level );
       update_buffer_level( level, all_hydro_vars, num_hydro_vars );
     }

   cart_debug("done updating initial conditions");

   /* set time variables */
   tStart = tl[min_level] = 0.0;

   dtl[min_level] = 0.5*constants->Myr/units->time;
   choose_timestep( &dtl[min_level] );

   for ( level = min_level+1; level <= max_level; level++ )
     {
       dtl[level] = 0.5*dtl[level-1];
       tl[level] = tl[min_level];
     }

   /* source */
   rtSingleSourceVal = N50*(units->time/constants->yr)*pow(constants->Mpc/units->length,3)/9.35e15/n0;
   rtSingleSourcePos[0] = rtSingleSourcePos[1] = rtSingleSourcePos[2] = 0.5*num_grid;
   
   //rtOtvetMaxNumIter = 30;

   cart_debug("done with initialization");
   
   check_map();
   
   run_output();
}
