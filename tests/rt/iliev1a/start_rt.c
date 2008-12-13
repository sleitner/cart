#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include <sys/types.h>
#include <unistd.h>

#include "defs.h"
#include "tree.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "timestep.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "viewdump.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "gravity.h"
#include "density.h"
#include "io.h"
#include "auxiliary.h"
#include "particle.h"
#include "starformation.h"

#include "rt_solver.h"
#include "rt_utilities.h"


#define N50             0.05
#define T_i             1.0e4
#define BottomLevel     3


double tStart;

void pix2ang_nest( long nside, long ipix, double *theta, double *phi);


#define refine_radius   (0.25*num_grid)

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

  if ( r < refine_radius )
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
  cell_gas_density(cell) = 1.0;
  cell_momentum(cell,0) = 0.0;
  cell_momentum(cell,1) = 0.0;
  cell_momentum(cell,2) = 0.0;
  cell_gas_gamma(cell) = (5.0/3.0);

  cell_gas_internal_energy(cell) = T_i*wmu/T0*aexp[0]*aexp[0]/(gamma-1)*(rt_XH+rt_XHe);

  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)*(gamma-1);
  cell_gas_energy(cell) = cell_gas_internal_energy(cell);

  cell_HI_density(cell) = (1.0-1.2e-3)*rt_XH;
  cell_HII_density(cell) = 1.2e-3*rt_XH;
  cell_HeI_density(cell) = rt_XHe;
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


float fHI(float r, double e[])
{
  int j, cell;
  float f, ff;
  double pos[3];

  for(j=0; j<3; j++)
    {
      pos[j] = 0.5*num_grid + r*e[j];
      if(pos[j] > num_grid) pos[j] -= num_grid;
      if(pos[j] < 0.0) pos[j] += num_grid;
    }
  cell = cell_find_position(pos);
  if(cell > -1)
    {
      f = cell_HI_fraction(cell);
    }
  else
    {
      f = -1.0;
    }
  MPI_Allreduce(&f,&ff,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);

  return ff;
}


void FindIFront(float val, float *riAvg, float *riMin, float *riMax)
{
  const long nside = 16;
  const long npix = 12*nside*nside;
  long ipix;
  int cell, j;
  double theta, phi, ravg, rmin, rmax;
  int navg;
  double ra, rb, rc, pos[3], e[3];
  float f, fa, fb, fc, f0 = val*rt_XH;

  ravg = 0.0;
  navg = 0;
  rmin = 1.0e20;
  rmax = 0.0;

  for(ipix=0; ipix<npix; ipix++)
    {
      pix2ang_nest(nside,ipix,&theta,&phi);
      e[0] = sin(theta)*cos(phi);
      e[1] = sin(theta)*sin(phi);
      e[2] = cos(theta);
        
      ra = 0.0;
      fa = fHI(ra,e);

      rb = 0.25*num_grid;
      fb = fHI(rb,e);

      if(fa<f0 && fb>f0)
	{
	  while(fabs(rb-ra)/(1.0e-10+ra+rb) > 0.001)
	    {
	      rc = 0.5*(ra+rb);
	      fc = fHI(rc,e);
	      if(fc < f0)
		{
		  ra = rc;
		}
	      else if(fc > f0)
		{
		  rb = rc;
		}
	      else
		{
		  ra = rb = rc;
		}
	    }            

	  rc = 0.5*(ra+rb);

	  rmin = min(rmin,rc);
	  rmax = max(rmax,rc);
	  ravg += rc;
	  navg++;
	}
    }

  if(navg > 0)
    {
      *riAvg = 1.0e3/hubble*r0*aexp[0]*ravg/navg;
      *riMin = 1.0e3/hubble*r0*aexp[0]*rmin;
      *riMax = 1.0e3/hubble*r0*aexp[0]*rmax;
    }
  else
    {
      *riAvg = *riMin = *riMax = 0.0;
    }

}


void run_output()
{
  const int nvars = 15;
  const int nbin1 = 32 * (1 << BottomLevel);
  int varid[] = { RTU_FRACTION+RT_HVAR_OFFSET+0, HVAR_GAS_DENSITY, RTU_GAS_TEMPERATURE, RTU_CELL_LEVEL, RTU_LOCAL_PROC, RT_VAR_OT_FIELD, rt_freq_offset+0, rt_freq_offset+1, rt_freq_offset+2, rt_et_offset+0, rt_et_offset+1, rt_et_offset+2, rt_et_offset+3, rt_et_offset+4, rt_et_offset+5 };
  int nbin[] = { nbin1, nbin1, nbin1 };
  double bb[6];
  int done;
  float riAvg[3], riMin, riMax;
  double tPhys;
  char filename[99];
  FILE *f;

  tPhys = 1.0e-6*pow(aexp[0],2)*t0*(tl[0]-tStart);

  bb[0] = bb[2] = bb[4] = num_grid*(0.5-0.25);
  bb[1] = bb[3] = bb[5] = num_grid*(0.5+0.25);

  sprintf(filename,"OUT/out.%05d.bin",step);
  rtuWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);

  FindIFront(0.01,riAvg+2,&riMin,&riMax);
  FindIFront(0.1,riAvg+1,&riMin,&riMax);
  FindIFront(0.5,riAvg+0,&riMin,&riMax);

  done = 0;
  if(local_proc_id == MASTER_NODE)
    {
      if(step == 0)
	{
	  f = fopen("OUT/rf.res","w");
	}
      else
	{
	  f = fopen("OUT/rf.res","a");
	}
      cart_assert(f != 0);
      fprintf(f,"%6.2f %6.3f %6.3f %6.3f %6.3f %6.3f\n",(float)tPhys,riAvg[0],riMin,riMax,riAvg[1],riAvg[2]);
      fclose(f);

      printf("Output: %d,  Time: %lg\n",step,tPhys);
      if(tPhys > 999.0) done = 1;
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
   int i, j, species, id, level, cell;
   int num_level_cells;
   int *level_cells;
   float astart;
   double pos[3];

   /* set units */
   astart = 1;
   hubble = 1.0;
   Lbox = 4*6.6e-3/(astart*hubble);
   Omega0 = 1.0e-3*pow(astart,3)/(1.123e-5*hubble*hubble);
   Omegab0 = Omega0;
   OmegaL0 = 0.0;
   aexp[min_level] = astart;

   init_units();

   /* create array with all hydro variable indices */
   for ( i = 0; i < num_hydro_vars; i++ )
     {
       all_hydro_vars[i] = HVAR_GAS_DENSITY + i;
     }
   
   for ( i = 0; i < nDim; i++ )
     {
       refinement_volume_min[i] = 0.0;
       refinement_volume_max[i] = num_grid;
     }

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

   dtl[min_level] = 1.0e7/(t0*astart*astart);
   choose_timestep( &dtl[min_level] );

   for ( level = min_level+1; level <= max_level; level++ )
     {
       dtl[level] = 0.5*dtl[level-1];
       tl[level] = tl[min_level];
       aexp[level] = aexp[min_level];		
     }

   /* source */
   num_row = num_grid;
   num_particle_species = 2;
   particle_species_mass[0] = 1.0;
   particle_species_mass[1] = N50*t0*pow(astart,2)/(1.05e11*Omega0/hubble*pow(r0
,3));

   num_particles_total = 10000;
   particle_species_num[0] = num_particles_total - num_star_particles;
   particle_species_num[1] = num_star_particles;
   particle_species_indices[0] = 0;
   particle_species_indices[1] = particle_species_num[0];
   particle_species_indices[2] = num_particles_total;
  
   num_local_particles = 0;
   num_local_star_particles = 0;
  
   for(i=0; i<num_particles; i++)
     {
       particle_level[num_local_particles] = -1;
     }

   for(i=0; i<num_particles_total; i++)
     {
       if(i < num_star_particles)
         {
           species = 1;
           id = particle_species_num[0] + i;
           for(j=0; j<3; j++) pos[j] = 0.5*num_grid;
         }
       else
         {
           species = 0;
           id = i - num_star_particles;
           pos[0] = 0.5*num_grid + 0*0.49*num_grid*sin(i*1.0);
           pos[1] = 0.5*num_grid + 0*0.49*num_grid*sin(i*3.0);
           pos[2] = 0.5*num_grid + 0*0.49*num_grid*sin(i*7.0);
         }
      
       cell = cell_find_position(pos);
      
       /* purpose: identifies what type of root cell corresponds
        *  to the given index
        *
        *  returns: 1 if cell is local, 2 if cell is buffer, 
        *      0 if cell is non-local
        */
       if(cell!=-1 && root_cell_type(cell_parent_root_sfc(cell))==1)
         {
          
           for(j=0; j<nDim; j++)
             {
               particle_x[num_local_particles][j] = pos[j];
               particle_v[num_local_particles][j] = 0.0;
             }
          
           particle_id[num_local_particles] = id;
           particle_mass[num_local_particles] = particle_species_mass[species];
          
           particle_t[num_local_particles] = 0.0;
           particle_dt[num_local_particles] = dtl[cell_level(cell)];

           particle_level[num_local_particles] = cell_level(cell);

           if(species == 1)
             {
               star_tbirth[num_local_star_particles] = 0.0;
               star_initial_mass[num_local_star_particles] = particle_mass[num_local_particles];
               num_local_star_particles++;
             }
          
           num_local_particles++;
         }
     }
  
   build_particle_list();
   
   cart_debug("done with initialization");
   
   check_map();
   
   run_output();
}
