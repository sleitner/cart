
#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "hydro.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "../extra/healpix.h"
#include "../extra/ifrit.h"
#include "../run/step.h"
#include "../run/logging.h"

#include "../et/oldstyle_units.h"


const float N50 = 0.05;
const float T_i = 1.0e4;
const int BottomLevel = 2;


extern float rtSingleSourceVal;
extern double rtSingleSourcePos[nDim];

double tStart;
const float refine_radius = 0.25*num_grid;
extern int rtOtvetMaxNumIter;

void units_set_art(double OmegaM, double h, double Lbox);


void refine_level( int cell, int level )
{
  int j;
  double pos[nDim];
  float r;

  cart_assert( cell >= 0 && cell < num_cells );
  cart_assert( cell_level(cell) == level );
	
  cell_center_position(cell, pos);

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

  cell_gas_internal_energy(cell) = T_i/units->temperature/(constants->gamma-1);

  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)*(constants->gamma-1);
  cell_gas_energy(cell) = cell_gas_internal_energy(cell);

  cell_HI_density(cell) = (1.0-1.2e-3);
  cell_HII_density(cell) = 1.2e-3;
  cell_HeI_density(cell) = 1.0e-10;
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
  MPI_Allreduce(&f,&ff,1,MPI_FLOAT,MPI_MAX,mpi.comm.run);

  return ff;
}


void FindIFront(float val, float *riAvg, float *riMin, float *riMax)
{
  const long nside = 16;
  const long npix = 12*nside*nside;
  long ipix;
  double theta, phi, ravg, rmin, rmax;
  int navg;
  double ra, rb, rc, e[3];
  float fa, fb, fc, f0 = val;

  ravg = 0.0;
  navg = 0;
  rmin = 1.0e20;
  rmax = 0.0;

  for(ipix=0; ipix<npix; ipix++)
    {
      hp_pix2ang_nest(nside,ipix,&theta,&phi);
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

	  rmin = MIN(rmin,rc);
	  rmax = MAX(rmax,rc);
	  ravg += rc;
	  navg++;
	}
    }

  if(navg > 0)
    {
      *riAvg = units->length/constants->kpc*ravg/navg;
      *riMin = units->length/constants->kpc*rmin;
      *riMax = units->length/constants->kpc*rmax;
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
  int varid[] = { I_FRACTION+RT_HVAR_OFFSET+0, HVAR_GAS_DENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC, RT_VAR_OT_FIELD };
  int nbin[] = { nbin1, nbin1, nbin1 };
  double bb[6];
  int done;
  float riAvg[3], riMin, riMax;
  double tPhys;
  char filename[99];
  FILE *f;

  tPhys = units->time*(tl[0]-tStart)/constants->Myr;

  bb[0] = bb[2] = bb[4] = num_grid*(0.5-0.25);
  bb[1] = bb[3] = bb[5] = num_grid*(0.5+0.25);

  sprintf(filename,"OUT/out.%05d.bin",step);
  ifrit.OutputMesh(filename,max_level,nbin,bb,nvars,varid);

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

  MPI_Bcast(&done,1,MPI_INT,MASTER_NODE,mpi.comm.run);

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
   float astart, hubble;
   const float n0 = 1.0e-3;
   double pos[3];

   /* set units */
   astart = 1;
   hubble = 1;
   oldstyle_units_set(n0*pow(astart,3)/(1.123e-5*hubble*hubble),hubble,4*6.6e-3/(astart*hubble));

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

   max_dt = 10*constants->Myr/units->time;

   for ( level = min_level+1; level <= max_level; level++ )
     {
       tl[level] = tl[min_level];
     }

   /* source */
   num_particle_species = 2;
   particle_species_mass[0] = 1.0;
   particle_species_mass[1] = N50*(units->time/constants->yr)*pow(constants->Mpc/units->length,3)/9.35e15/n0;

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
           particle_dt[num_local_particles] = 0.0;

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
   
   //rtOtvetMaxNumIter = 30;

   cart_debug("done with initialization");
   
   check_map();
   
   run_output();
}
