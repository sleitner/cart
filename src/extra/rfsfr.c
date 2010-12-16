#include "config.h"


#include <math.h>
#include <stdio.h>
#include <string.h>


#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "rt_solver.h"
#include "rt_utilities.h"
#include "starformation.h"
#include "tree.h"
#include "units.h"

#include "extra/halo_finder.h"


#if defined (HYDRO) && defined(STARFORM) && defined(RADIATIVE_TRANSFER)
/*
//  Dump ISRF and a variable var[cell] with a hierarchy of levels
*/
void extRFvsSFR1(const char *froot, int top_level, float *var, const halo *h)
{
  const char *fext[] = { "sfr", "den", "dmw", "fh2" };
  const int nfiles = sizeof(fext)/sizeof(char*);
  MESH_RUN_DECLARE(level,cell);
  float uLen;
  int i, j, l, parent, size, rank, save;
  double dx, pos[nDim], r = 0.0;
  FILE *f[nfiles];
  char str[999], fsuffix[99];

  cart_assert(top_level>=min_level && top_level<=max_level);

  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  /*
  // Units
  */
  uLen = units->length/constants->kpc; /* phys pc */

  /*
  //  Create a set of files for each requested level
  */
  for(level=max_level_now_global(MPI_COMM_WORLD); level>=top_level; level--)
    {
      cart_debug("Working on level %d...",level);

      /*
      //  Write to a file in order of a proc rank
      */
      if(h != NULL)
	{
	  if(halo_level(h,MPI_COMM_WORLD) < level) continue;
	  sprintf(fsuffix,"L=%02d.res.%05d",level,h->id);
	}
      else
	{
	  sprintf(fsuffix,"L=%02d.res",level);
	}

      select_level(level,CELL_TYPE_LOCAL,&_Num_level_cells,&_Level_cells);

      for(i=0; i<size; i++)
	{
	  MPI_Barrier(MPI_COMM_WORLD);
	  if(i == rank)
	    {
	      cart_debug("Writing file piece #%d",i);

	      for(j=0; j<nfiles; j++)
		{
		  sprintf(str,"%s-%s.%s",froot,fext[j],fsuffix);
		  f[j] = fopen(str,(i==0?"w":"a"));
		  cart_assert(f[j] != NULL);
		}

	      /*
	      //  Manifest
	      */
	      if(i == 0)
		{
		  for(j=0; j<nfiles; j++) fprintf(f[j],"%d %9.3e\n",level,uLen*cell_size[level]);
		}

	      MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

	      if(cell_is_leaf(cell))
		{
		  save = 0;
		  for(parent=cell,l=level; l>=min_level; l--)
		    {
		      cart_assert(parent > -1);
		      if(var[parent] > 0.0) save = 1;
		      parent = cell_parent_cell(parent);
		    }

		  if(save && h!=NULL)
		    {
		      cell_center_position(cell,pos);
		      for(j=0, r=0.0; j<nDim; j++)
			{
			  dx = pos[j] - h->pos[j];
			  if(dx < -0.5*num_grid) dx += num_grid;
			  if(dx >  0.5*num_grid) dx -= num_grid;
			  r += dx*dx;
			}
		      r = sqrt(r);
		      if(r > h->rhalo) save = 0;
		    }

		  if(save)
		    {
		      for(j=0; j<nfiles; j++) fprintf(f[j],"%9.3e",rtUmw(cell));
		      for(parent=cell,l=level; l>=min_level; l--)
			{
			  if(nfiles > 0) fprintf(f[0]," %9.3e",var[parent]);
			  if(nfiles > 1) fprintf(f[1]," %9.3e",units->number_density*cell_gas_density(parent));
			  if(nfiles > 2) fprintf(f[2]," %9.3e",rtDmw(parent));
			  if(nfiles > 3) fprintf(f[3]," %9.3e",cell_H2_fraction(parent));
			  parent = cell_parent_cell(parent);
			}
		      if(h != NULL)
			{
			  for(j=0; j<nfiles; j++) fprintf(f[j]," %9.3le",uLen*r);
			}
		      for(j=0; j<nfiles; j++) fprintf(f[j],"\n");
		    }
		}
	      MESH_RUN_OVER_CELLS_OF_LEVEL_END;

	      for(j=0; j<nfiles; j++) fclose(f[j]);
	    }
	}

      cart_free(_Level_cells);
      _Level_cells = 0;
    }
}


void extRFvsSFR(const char *froot, int top_level, const halo_list *halos)
{
  MESH_RUN_DECLARE(level,cell);
  int j, ih;
  float *var, *sfr;
  float uDen, uTime, uRate;

  /*
  // Units
  */
  uDen = units->density*pow(constants->kpc,3.0)/constants->Msun;    /* Msun/kpc^3 */
  uTime = units->time/constants->yr;      /* yr */
  uRate = uDen/uTime;                     /* Msun/kpc^3/yr */

  var = cart_alloc(float,num_cells);
  memset(var,0,sizeof(float)*num_cells);

  /*
  //  Fill in the data array
  */
  MESH_RUN_OVER_LEVELS_BEGIN(level,_MaxLevel,min_level);

  /*
  //  Prepare the SFR buffer array
  */
  sfr = cart_alloc(float,_Num_level_cells);
  star_formation_rate(level,_Num_level_cells,_Level_cells,sfr);

#pragma omp parallel for default(none), private(_Index,cell,j), shared(_Num_level_cells,_Level_cells,sfr,var,uRate,cell_child_oct)
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);

  if(cell_is_leaf(cell))
    {
      var[cell] = uRate*max(0.0,sfr[_Index]);
    }
  else
    {
      var[cell] = 0.0;
      for(j=0; j<num_children; j++)
	{
	  var[cell] += var[cell_child(cell,j)]/num_children;
	}
    }

  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  
  cart_free(sfr);

  MESH_RUN_OVER_LEVELS_END;

  if(halos == NULL)
    {
      extRFvsSFR1(froot,top_level,var,NULL);
    }
  else
    {
      for(ih=0; ih<halos->num_halos; ih++) if(halo_level(&halos->list[ih],MPI_COMM_WORLD) >= top_level)
	{
	  extRFvsSFR1(froot,top_level,var,&halos->list[ih]);
	}
    }

  cart_free(var);
}

#endif
