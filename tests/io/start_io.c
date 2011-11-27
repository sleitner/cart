#include "config.h"

#include <stdio.h>
#include <math.h>

#include "auxiliary.h"
#include "hydro.h"
#include "iterators.h"
#include "particle.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"
#include "cosmology.h"


#define BottomLevel     1



#define refine_radius   (0.4*num_grid)

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

  if ( r < refine_radius*(0.5+0.5/(level+1.0)) )
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
  double pos[nDim];
  float r;
  float rho, dxscal = 0.1*num_grid;

  cell_center_position(cell, pos);

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

  cell_gas_internal_energy(cell) = (1.0e3/units->temperature)*rho;

  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)*(cell_gas_gamma(cell)-1);
  cell_gas_energy(cell) = cell_gas_internal_energy(cell);
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
}


void init_run()
{
   int i, j, species, id, level, cell;
   int num_level_cells;
   int *level_cells;
   float astart;
   double pos[3];
   float hubble;

   /* set units */
   astart = 0.1;
   hubble = 1.0;
   cosmology_set(h,hubble);
   cosmology_set(OmegaL,0.0);
   cosmology_set(OmegaM,3.2*pow(astart,3)/(1.123e-5*hubble*hubble));
   cosmology_set(OmegaB,cosmology->OmegaM);
   box_size = 1.6e-3/(astart*hubble);
   units_set_art(cosmology->OmegaM,cosmology->h,box_size);
   abox[min_level] = astart;

   units_reset();
   
   for ( i = 0; i < nDim; i++ )
     {
       refinement_volume_min[i] = 0.0;
       refinement_volume_max[i] = num_grid;
     }

   cart_debug("in init");

   set_rt_initial_conditions();

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

   max_dt = 1.0e7/(units->time*constants->yr);

   for ( level = min_level+1; level <= max_level; level++ )
     {
       tl[level] = tl[min_level];
       abox[level] = abox[min_level];		
     }

   /* particles */
   num_row = num_grid;
   num_particle_species = 2;
   particle_species_mass[0] = 1.0;
   particle_species_mass[1] = 0.1;

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
	   pos[0] = 0.5*num_grid + 0.3*num_grid*cos(1.0*i);
	   pos[1] = 0.5*num_grid + 0.3*num_grid*sin(1.0*i);
	   pos[2] = 0.5*num_grid;
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

   cart_debug("done with initialization");
   
   check_map();
   
   run_output();
}
