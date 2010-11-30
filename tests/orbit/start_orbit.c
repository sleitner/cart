#include "config.h"


#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include "auxiliary.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "rt_utilities.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#include "extra/ifrit.h"


#define depth  (max_level-min_level-3)

const float phi0 = M_PI/180*0;
const float vel0 = 1.3;
const float fac0 = 2.0;

const int force_free = 0;


const double scale = 1 << depth;


void refine_level(int cell, int level)
{
  int j;
  double pos[nDim];
  float r, rmax;

  cart_assert( cell >= 0 && cell < num_cells );
  cart_assert( cell_level(cell) == level );
	
  cell_center_position(cell, pos);

  for(j=0; j<nDim; j++)
    {
      pos[j] -= 0.5*num_grid;
    }

  r = scale*max(fabs(pos[0]),max(fabs(pos[1]),fabs(pos[2])));

  if(level < min_level+depth)
    {
      rmax = 8.0 + 2*scale*cell_size[level];
    }
  else
    {
      level -= (min_level+depth);
      rmax = 8;
      rmax *= pow(0.5,(double)(1+level));
    }

  if(r < rmax)
    {
      refinement_indicator(cell,0) = 1.0;	
    }
  else
    {
      refinement_indicator(cell,0) = 0.0;
    }
}


void save_pos(int level)
{
  static int init = 1;
  int j;
  FILE *f;

  f = fopen("orbit.res",(init==1)?"w":"a");
  cart_assert(f != NULL);
  init = 0;

  for(j=0; j<num_particles; j++)
    {
      if(particle_id[j]==1 && particle_level[j]==level)
	{
	  fprintf(f,"%9.3le %6.4lf %6.4lf %6.4lf %le %le %le %e %d\n",pow(scale,1.5)*tl[level],scale*(particle_x[j][0]-0.5*num_grid),scale*(particle_x[j][1]-0.5*num_grid),scale*(particle_x[j][2]-0.5*num_grid),particle_v[j][0],particle_v[j][1],particle_v[j][2],particle_pot[j],particle_level[j]);
	}
    }

  fclose(f);
}


void run_output()
{
  const int prolongation_vars[1] = { VAR_POTENTIAL };
  const int nvars = 3;
  const int nbin1 = 256;
  int varid[] = { I_CELL_LEVEL, I_LOCAL_PROC, VAR_DENSITY };
  int nbin[] = { nbin1, nbin1, nbin1 };
  double bb[6];
  MESH_RUN_DECLARE(level,cell);
  FILE *f;
  double r, pos[3], p0;
  int j;

  if(step != 0) return;


  cart_debug("Boundary is at %lg AU",0.5*num_grid*scale);


  bb[0] = bb[2] = bb[4] = 0.5*num_grid - 8.0/scale; // 16 AU
  bb[1] = bb[3] = bb[5] = 0.5*num_grid + 8.0/scale;
  //ifrit.OutputMesh("mesh1.bin",max_level,nbin,bb,nvars,varid);

  bb[0] = bb[2] = bb[4] = 0.0;
  bb[1] = bb[3] = bb[5] = num_grid;
  //ifrit.OutputMesh("mesh2.bin",max_level,nbin,bb,nvars,varid);

  f = fopen("pot.res","w");
  cart_assert(f);

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

  //if(cell_is_leaf(cell))
    {
      cell_center_position(cell,pos);
      for(j=0; j<nDim; j++)
	{
	  pos[j] -= 0.5*num_grid;
	}

      r = scale*sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);

      p0 = scale*4*pow(M_PI,3.0)/(3*r);

      fprintf(f,"%d %le %le %le\n",level,r,cell_potential(cell),p0);

      //cell_potential(cell) = -p0;
    }

  MESH_RUN_OVER_CELLS_OF_LEVEL_END;

  update_buffer_level(level,prolongation_vars,1);

  MESH_RUN_OVER_LEVELS_END;

  fclose(f);
}


void init_run()
{
  int i, j, cell, level;
  int num_level_cells, *level_cells;
  double r, pos[3];

  units_set(constants->Msun/num_root_cells,pow(scale,1.5)*constants->yr,scale*constants->pc/206264.806);
  units_reset();

  cart_debug("in init");

  /* build buffer */
  build_cell_buffer();
  cart_debug("built cell buffer");

  repair_neighbors();
  check_map();
  cart_debug("repaired neighbors");

  /* do initial refinements */
  for(level=min_level; level<max_level; level++)
    {
      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells );
      for(i=0; i<num_level_cells; i++)
	{
	  refine_level(level_cells[i],level);
	}
      cart_free( level_cells );
      refine(level);
    }

  num_level_cells = 0;
  for(level=min_level; level<=max_level; level++)
    {
      cart_debug("num_level_cells[%u] = %u",level,num_cells_per_level[level]);
      num_level_cells += num_cells_per_level[level];
    }
  cart_debug("Total cells: %u",num_level_cells);

  cart_debug("set initial conditions");

  /* set time variables */
  tl[min_level] = 0.0;
   
  num_row = num_grid;
  num_particle_species = 2;
  particle_species_mass[0] = force_free ? 0 : num_root_cells;
  particle_species_mass[1] = force_free ? num_root_cells : 0;

  num_particles_total = 2;
  particle_species_num[0] = 1;
  particle_species_num[1] = 1;
  particle_species_indices[0] = 0;
  particle_species_indices[1] = 1;
  particle_species_indices[2] = 2;
  
  num_local_particles = 0;

  pos[0] = 0.5*num_grid;
  pos[1] = 0.5*num_grid;
  pos[2] = 0.5*num_grid;
      
  for(i=0; i<num_particles; i++)
    {
      particle_level[i] = -1;
    }

  for(i=0; i<num_particles; i++)
    {
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

	  particle_id[num_local_particles] = i;
	  particle_mass[num_local_particles] = particle_species_mass[i];

	  particle_t[num_local_particles] = 0.0;
	  particle_dt[num_local_particles] = 0.0;

	  particle_level[num_local_particles] = cell_level(cell);

	  if(i == 1)
	    {
	      particle_v[num_local_particles][0] = 6.2830*vel0*sin(phi0)*sqrt(scale*fac0);
	      particle_v[num_local_particles][1] = 6.2830*vel0*cos(phi0)*sqrt(scale*fac0);
	      particle_v[num_local_particles][2] = 0;
	      for(j=0; j<nDim; j++)
		{
		  pos[j] -= 0.5*num_grid;
		}
	      r = sqrt(1.0e-4+pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);

	      particle_pot[num_local_particles] = -4*pow(M_PI,3.0)/(3*r);
	    }

	  num_local_particles++;
	  pos[0] += 1.05/(scale*fac0)*cos(phi0);
	  pos[1] -= 1.05/(scale*fac0)*sin(phi0);

	}
    }
  
  build_particle_list();
   
  cart_debug("done with initialization");
   
  check_map();

  units_update(min_level);

  cart_debug("units->potential = %le",units->potential);

  for(level=min_level; level<=max_level; level++) save_pos(level);
}
