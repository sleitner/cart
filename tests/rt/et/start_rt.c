#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include <sys/types.h>
#include <unistd.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "density.h"
#include "hydro.h"
#include "iterators.h"
#include "logging.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "rt_utilities.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "extra/healpix.h"
#include "extra/ifrit.h"

const int BottomLevel = 6;

extern float rtSingleSourceVal;
extern double rtSingleSourcePos[nDim];

void units_set_art(double OmegaM, double h, double Lbox);


void refine_level( int cell, int level )
{
  int j;
  double pos[nDim];
  float r = 0.25*pow(0.5,level)*num_grid;
  
  cart_assert( cell >= 0 && cell < num_cells );
  cart_assert( cell_level(cell) == level );
  
  cell_center_position(cell, pos);

  for(j=0; j<nDim; j++)
    {
      pos[j] -= rtSingleSourcePos[j];
      if(pos[j] >  num_grid/2) pos[j] -= num_grid;
      if(pos[j] < -num_grid/2) pos[j] += num_grid;
    }

  if(fabs(pos[0])<r && fabs(pos[1])<r && fabs(pos[2])<r)
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

  cell_gas_internal_energy(cell) = 1.0e4/units->temperature/(constants->gamma-1)*(constants->XH+constants->XHe);

  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)*(constants->gamma-1);
  cell_gas_energy(cell) = cell_gas_internal_energy(cell);

  cell_HI_density(cell) = constants->XH;
  cell_HII_density(cell) = 0.0;
  cell_HeI_density(cell) = constants->XHe;
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


void SaveProfile(char *filename)
{
  const long nside = 16;
  const long npix = 12*nside*nside;

  const float lrmin = log10(0.3*pow(0.5,BottomLevel));
  const float lrmax = log10(1.0*num_grid);
  const float dlr = 0.05;
  const int n = 1 + (int)((lrmax-lrmin)/dlr);

  int i, cell;
  long ipix;
  float lr;
  double theta, phi, s1, s2, w;
  int ns;
  double r, pos[3];
  FILE *f;

  f = fopen(filename,"w");
  cart_assert(f != 0);

  for(i=0; i<n; i++)
    {
      s1 = s2 = 0.0;
      ns = 0;

      lr = lrmin + dlr*i;
      r = pow(10.0,(double)lr);

      for(ipix=0; ipix<npix; ipix++)
	{
	  hp_pix2ang_nest(nside,ipix,&theta,&phi);
	  pos[0] = rtSingleSourcePos[0] + r*sin(theta)*cos(phi);
	  pos[1] = rtSingleSourcePos[1] + r*sin(theta)*sin(phi);
	  pos[2] = rtSingleSourcePos[2] + r*cos(theta);

	  cell = cell_find_position(pos);
	  if(cell > -1)
	    {
	      w = (double)cell_var(cell,RT_VAR_OT_FIELD);
	      s1 += w;
	      s2 += w*w;
	      ns++;
	    }
	}

      if(ns > 0)
	{
	  s1 /= ns;
	  s2 = s2/ns - s1*s1;
	  if(s2 > 0.0) s2 = sqrt(s2); else s2 = 0.0;
	}

      pos[0] = rtSingleSourcePos[0] + r;
      pos[1] = rtSingleSourcePos[1];
      pos[2] = rtSingleSourcePos[2];
      
      cell = cell_find_position(pos);
      if(cell > -1)
	{
	  w = cell_var(cell,rt_et_offset+0);
	}
      else w = 0.0;

      fprintf(f,"%g %g %g %g\n",lr,(float)s1,(float)s2,(float)w);
    }

  fclose(f);
}


void rtOtvetTreeEmulatorEddingtonTensor(int level);
void rtOtvetSingleSourceEddingtonTensor(int level, float srcVal, double *srcPos);


void run_output()
{
  const int nvars = 8;
  int i, level, varid[nvars], nbin[3];
  double bb[6];

  for(i=0; i<6; i++) varid[i] = rt_et_offset + i;
  varid[6] = RT_VAR_OT_FIELD;
  varid[7] = I_CELL_LEVEL;

  nbin[0] = num_grid;
  for(i=0; i<=BottomLevel && nbin[0]<256; i++) nbin[0] *= 2;
  nbin[1] = nbin[2] = nbin[0];
  
  bb[0] = bb[2] = bb[4] = 0.5*num_grid - pow(0.5,1.0+BottomLevel)*nbin[0];
  bb[1] = bb[3] = bb[5] = 0.5*num_grid + pow(0.5,1.0+BottomLevel)*nbin[0];
  
  for(level=min_level; level<=BottomLevel; level++)
    {
      cart_debug("assigning density on level %u", level );
      assign_density( level );
    }
  /*
  //  Create outputs: single source mode
  */
  for(level=min_level; level<=BottomLevel; level++)
    {
      cart_debug("computing Single Source ET on level %u", level );
      rtOtvetSingleSourceEddingtonTensor(level,rtSingleSourceVal,rtSingleSourcePos);
    }

//  rtuWriteIfritFile(max_level,nbin,bb,nvars,varid,"OUT/out_ss.bin");
  SaveProfile("OUT/prof_ss.res");

  /*
  //  Create outputs: tree emulator mode
  */
  for(level=min_level; level<=BottomLevel; level++)
    {
      cart_debug("computing TreeEmulator ET on level %u", level );
      rtOtvetTreeEmulatorEddingtonTensor(level);
    }

  ifrit.OutputMesh("OUT/out_te.bin",max_level,nbin,bb,nvars,varid);
  //SaveProfile("OUT/prof_te.res");

  finalize_logging();
  MPI_Finalize();
  exit(0);
}


void init_run()
{
  int i, j, species, id, level, cell;
  int num_level_cells;
  int *level_cells;
  float astart, hubble;
  double pos[3];

  /* set units */
  astart = 1;
  hubble = 1;
  units_set_art(1.0e-3*pow(astart,3)/(1.123e-5*hubble*hubble),hubble,4*6.6e-3/(astart*hubble));

  units_reset();
  units_update(min_level);

  /* source */
  rtSingleSourceVal = 1.0;
  rtSingleSourcePos[0] = rtSingleSourcePos[1] = rtSingleSourcePos[2] = 0.5*num_grid - 0.0;

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
  
  cart_debug("set initial conditions");
  
  for ( level = min_level; level <= max_level; level++ )
    {
      cart_debug("updating level %u", level );
      update_buffer_level( level, all_hydro_vars, num_hydro_vars );
    }

  cart_debug("done updating initial conditions");

  /* set time variables */
  tl[min_level] = 0.0;
  
  dtl[min_level] = 10*constants->Myr/units->time;

  for ( level = min_level+1; level <= max_level; level++ )
    {
      tl[level] = tl[min_level];
    }

  /* particles */
  num_row = num_grid;
  num_particle_species = 2;
  particle_species_mass[0] = 1.0;
  particle_species_mass[1] = 1.0;

  num_particles_total = num_particles;
  particle_species_num[0] = num_particles - num_star_particles;
  particle_species_num[1] = num_star_particles;
  particle_species_indices[0] = 0;
  particle_species_indices[1] = particle_species_num[0];
  particle_species_indices[2] = num_particles_total;
  
  num_local_particles = 0;
  num_local_star_particles = 0;
  
  for(i=0; i<num_particles_total; i++)
    {
      if(i < num_star_particles)
	{
	  species = 1;
	  id = particle_species_num[0] + i;
	  for(j=0; j<3; j++) pos[j] = rtSingleSourcePos[j];
	}
      else
	{
	  species = 0;
	  id = i - num_star_particles;
	  pos[0] = 0.5*num_grid + 0.2*num_grid*cos(1.0*i);
	  pos[1] = 0.5*num_grid + 0.2*num_grid*sin(1.0*i);
	  pos[2] = 0.5*num_grid;
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
