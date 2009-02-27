#include "defs.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "auxiliary.h"
#include "density.h"
#include "hydro.h"
#include "io.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "sfc.h"
#include "starformation.h"
#include "timestep.h"
#include "tree.h"
#include "units.h"

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#ifdef RT_DEBUG
#include "rt_debug.h"
#endif
#endif

#include "rt_utilities.h"
#include "extra/ifrit.h"

#include "gic_reader.h"


extern int num_options;
extern char **options;


/*
//  Helper functions
*/
void FindMaxVar(int var, float *val, double *pos)
{
  MESH_RUN_DECLARE(level,cell);
  float vMax = -1.0e35;
  int cellMax = 0;

  cart_assert(var>=0 && var<num_vars);

  MESH_RUN_OVER_ALL_LEVELS_BEGIN(level);
  MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(cell);
  if(cell_is_leaf(cell) && vMax<cell_var(cell,var))
    {
      vMax = cell_var(cell,var);
      cellMax = cell;
    }
  MESH_RUN_OVER_CELLS_OF_LEVEL_END;
  MESH_RUN_OVER_LEVELS_END;
  
  *val = vMax;
  if(pos != NULL) cell_position_double(cellMax,pos);
}


void gicStartFile(char *filename, struct gicFile *input, struct gicManifest *manifest, struct gicFileHeader *fileHeader)
{
  int ret;

  input->File = fopen(filename,"r");
  if(input->File == NULL)
    {
      cart_error("Unable to open file %s for reading",filename);
    }

  ret = gicReadManifest(input,manifest);
  if(ret != 0)
    {
      cart_error("Unable to read GIC manifest of file %s, error code: %d",filename,ret);
    }

  ret = gicReadFileHeader(input,fileHeader);
  if(ret != 0)
    {
      cart_error("Unable to read GIC file header of file %s  error code: %d",filename,ret);
    }

  /*
  //  Only low-res files are currently supported for the single-res case
  */
  if(num_grid!=fileHeader->dims[0] || num_grid!=fileHeader->dims[1] || num_grid!=fileHeader->dims[2])
    {
      cart_error("Input grid at the top level (%d,%d,%d) of file %s is not a cube of size %d\n",fileHeader->dims[0],fileHeader->dims[1],fileHeader->dims[2],filename,num_grid);
    }

  /*
  //  Set units
  */
  Omega0 = manifest->OmegaX + manifest->OmegaB;
  Omegab0 = manifest->OmegaB;
  OmegaL0 = manifest->OmegaL;
  hubble = manifest->h100;
  Lbox = manifest->dx*num_grid;

  aexp[min_level] = fileHeader->aBegin;
  tl[min_level] = a2b(aexp[min_level]);
}


long gicSkipToLevel(int Level, int lMax, struct gicFile *input, struct gicFileHeader *fileHeader, struct gicLevelHeader *levelHeader, int narr, int indexed)
{
  int i, j, l, L, ret;
  int page, num_pages;
  long offset = 0;

  /*
  //  This is not very efficient, but easier to program and to read files backward.
  //  First, count the length to skip.
  */
  offset += (2*sizeof(GIC_RECORD)+GIC_MANIFEST_SIZE);
  offset += (2*sizeof(GIC_RECORD)+GIC_FILEHEADER_SIZE);
  if(fileHeader->Lmax > 0)  /* Skip mask (if needed) */
    {
      num_pages = (num_root_cells+fileHeader->Nrec-1)/fileHeader->Nrec;
      offset += num_pages*(2*sizeof(GIC_RECORD)+fileHeader->Nrec*sizeof(GIC_INTEGER));
    }

  /*
  //  Actually skip this distance along all input streams.
  */
  for(l=0; l<lMax; l++)
    {
      ret = fseek(input[l].File,offset,SEEK_SET);
      if(ret != 0)
	{
	  cart_error("Error in reposition stream %d, fseek error %d",l,ret);
	}
    }

  /*
  //  Skip lower levels (work is done only for a multi-res file)
  */
  offset = 0;
  for(L=0; L<Level; L++)
    {
      /*
      //  Read level header for stream 0 - we need to know the amount of data to skip
      */
      ret = gicReadLevelHeader(input,levelHeader);
      if(ret != 0)
	{
	  cart_error("Unable to read GIC level header for input stream 0, level %d, error code: %d",L,ret);
	}
      
      /*
      //  Skip other streams
      */
      for(l=1; l<lMax; l++)
	{
	  ret = fseek(input[l].File,2*sizeof(GIC_RECORD)+GIC_LEVELHEADER_SIZE,SEEK_CUR);
	  if(ret != 0)
	    {
	      cart_error("Error in reposition stream %d, fseek error %d",l,ret);
	    }
	}

      offset += levelHeader->Nlev;
      num_pages = (levelHeader->Nlev+fileHeader->Nrec-1)/fileHeader->Nrec;

      /*
      //  Skip narr arrays
      */
      for(i=0; i<narr; i++)
	{
	  for(page=0; page<num_pages; page++)
	    {
	      for(l=0; l<lMax; l++)
		{
		  if(gicSkipFortranRecordReal(input+l) != 0)
		    {
		      cart_error("Error in reposition stream %d, fseek error %d",l,ret);
		    }
		  offset += fileHeader->Nrec*sizeof(GIC_REAL);
		  if(indexed)
		    {
		      for(j=0; j<3; j++) if(gicSkipFortranRecordInteger(input+l) != 0)
			{
			  cart_error("Error in reposition stream %d, fseek error %d",l,ret);
			}
		      offset += 3*fileHeader->Nrec*sizeof(GIC_INTEGER);
		    }
		}
	    }
	}
    }

  return offset;
}


/* 
//  This function is a modified version of restart_load_balance from io.c
*/
void gicBalanceLoad(char *rootname, char *type)
{
  int i, j;
  int L, nref;
  int index;
  int coords[3];
  int page, num_pages;
  float *cell_work;	
  int *constrained_quantities;

  struct gicFile input;
  struct gicManifest manifest;
  struct gicFileHeader fileHeader;
  char filename[257];   /* GIC filenames are limited to 256 bytes */
  GIC_INTEGER *mask;

	
  if(num_procs == 1)
    {
      proc_sfc_index[0] = 0;
      proc_sfc_index[1] = num_root_cells;
      init_tree();
      return;
    }

  if(strlen(rootname) > 250)
    {
      cart_error("GIC filenames are limited to 256 bytes; rootnames to 250 bytes");
    }

  if(local_proc_id == MASTER_NODE)
    {
      /* do load balancing */
      constrained_quantities = (int *)cart_alloc(num_constraints*num_root_cells*sizeof(int));
      cell_work = (float *)cart_alloc(num_root_cells*sizeof(float));

      for(i=0; i<num_root_cells; i++)
	{
	  cell_work[i] = 0.0;
	}

      for(i=0; i<num_constraints*num_root_cells; i++)
	{
	  constrained_quantities[i] = 0;
	}

      /*
      //  load mask information and compute work
      */
      strcpy(filename,rootname);
      strcat(filename,"_");
      strcat(filename,type);
      strcat(filename,".vel");

      gicStartFile(filename,&input,&manifest,&fileHeader);

      /*
      //  Only low-res files are currently supported for the single-res case
      */
      L = 0;
#ifdef PARTICLES
      num_particles_total = fileHeader.Ntot;
#endif

      if(fileHeader.Lmax == 0)
	{
	  /*
	  //  Single-resolution file - all work is uniform
	  */
	  nref = 1 << (nDim*L);
	  for(index=0; index<num_root_cells; index++)
	    {
	      constrained_quantities[num_constraints*index+0] = nref;
	      cell_work[index] += cost_per_cell*nref;
#ifdef PARTICLES
	      constrained_quantities[num_constraints*index+1] = nref;
	      cell_work[index] += cost_per_particle*nref;
#endif
	    }
	}
      else
	{

	  mask = (int *)cart_alloc(fileHeader.Nrec*sizeof(GIC_INTEGER));
	  num_pages = (num_root_cells+fileHeader.Nrec-1)/fileHeader.Nrec;

	  coords[0] = coords[1] = coords[2] = 0;
	  for(page=0; page<num_pages; page++)
	    {
	      if(gicReadFortranRecordInteger(&input,mask) != 0)
		{
		  cart_error("File is corrupted, error in reading mask array on page %d",page);
		}
	      
	      for(j=0; coords[2]<num_grid && j<fileHeader.Nrec; j++)
		{
		  
		  index = sfc_index(coords);
		  cart_assert(index>=0 && index<num_root_cells);
		  
		  nref = 1 << mask[j];

		  constrained_quantities[num_constraints*index+0] = nref;
		  cell_work[index] += cost_per_cell*nref;
#ifdef PARTICLES
		  constrained_quantities[num_constraints*index+1] = nref;
		  cell_work[index] += cost_per_particle*nref;
#endif

		  coords[0]++;
		  if(coords[0] == num_grid)
		    {
		      coords[0] = 0;
		      coords[1]++;
		    }
		  if(coords[1] == num_grid)
		    {
		      coords[1] = 0;
		      coords[2]++;
		    }
		}
	    }
	}

      cart_free(mask);
      fclose(input.File);

      cart_debug("load balancing before i/o");
      load_balance_entire_volume(cell_work,constrained_quantities,proc_sfc_index);

      cart_free(cell_work);
      cart_free(constrained_quantities);
    
    }

  /* let all other processors know what their new workload is */
  MPI_Bcast(proc_sfc_index,num_procs+1,MPI_INT,MASTER_NODE,MPI_COMM_WORLD);
  init_tree();
}


/*
//  ******************
//  *
//  *    PARTICLES
//  *
//  ******************
*/
#ifdef PARTICLES

void gicReadParticleLevel(int species, int Level, int lMax, struct gicFile *input, struct gicFileHeader *fileHeader, struct gicLevelHeader *levelHeader, long *idOffset, GIC_REAL *buffer[])
{
  int j, n, l, ret;
  int page, num_pages;
  float xFac, vFac;
  GIC_REAL *x, *y, *z, *vx, *vy, *vz;
  long noff, ntot, id, offset;
  int ipart, index;
  int coords[nDim];

  x = buffer[0];
  y = buffer[1];
  z = buffer[2];
  vx = buffer[3];
  vy = buffer[4];
  vz = buffer[5];

  gicSkipToLevel(Level,lMax,input,fileHeader,levelHeader,3,0);

  id = *idOffset;

  /*
  //  Read level data (if needed)
  */
  if(fileHeader->Lmax > 0)
    {
      for(l=0; l<lMax; l++)
	{
	  ret = gicReadLevelHeader(input+l,levelHeader+l);
	  if(ret != 0)
	    {
	      cart_error("Unable to read GIC level header for input stream %d, level %d, error code: %d",l,Level,ret);
	    }
	  if(Level!=levelHeader[l].L || fileHeader->Lmax!=levelHeader[l].Lmax)
	    {
	      cart_error("Input stream %d is corrupted, error code 010",l);
	    }
	  if(l!=0 && memcmp(levelHeader,levelHeader+l,sizeof(struct gicLevelHeader))!=0)
	    {
	      cart_error("Incompatible input streams 0 and %d at level %d",l,Level);
	    }
	}

      /*
      //  Set species arrays
      */
      particle_species_mass[species] = levelHeader->Mlev/aM0;
      ntot = levelHeader->Nlev;
    }
  else
    {
      /*
      //  Set species arrays
      */
#ifdef HYDRO
      particle_species_mass[species] = 1.0 - Omegab0/Omega0;
#else
      particle_species_mass[species] = 1.0;
#endif
      ntot = fileHeader->Ntot;
    }

  particle_species_num[species] = ntot;
  particle_species_indices[species+1] = ntot + particle_species_indices[species];
  
  cart_debug("particle_species_mass[%d] = %e",species,particle_species_mass[species]);
  cart_debug("particle_species_num[%d] = %u",species,particle_species_num[species]);

  num_pages = (ntot+fileHeader->Nrec-1)/fileHeader->Nrec;

  xFac = 1.0/r0;
  vFac = fileHeader->aBegin/v0;

  /*
  // Skip to the appropriate record
  */
  for(l=0; l<lMax; l++)
    {
      offset = (l%3)*num_pages*(2*sizeof(GIC_RECORD)+fileHeader->Nrec*sizeof(GIC_REAL));
      ret = fseek(input[l].File,offset,SEEK_CUR);
      if(ret != 0)
	{
	  cart_error("Error in reposition stream %d, fseek error %d",l,ret);
	}
    }

  noff = 0;
  for(page=0; page<num_pages; page++)
    {
      n = fileHeader->Nrec;
      if(noff+n > ntot) n = ntot - noff;
      noff += fileHeader->Nrec;

      for(l=0; l<lMax; l++)
	{
	  if(gicReadFortranRecordReal(input+l,buffer[l]) != 0)
	    {
	      cart_error("Input stream %d is corrupted, error code 020",l);
	    }
	}

      for(j=0; j<n; j++)
	{
	  /*
	  //  Units
	  */
	  x[j] *= xFac;
	  y[j] *= xFac;
	  z[j] *= xFac;

	  vx[j] *= vFac;
	  vy[j] *= vFac;
	  vz[j] *= vFac;

	  /* enforce periodic boundary conditions */
	  if(x[j] < 0.0)
	    {
	      x[j] += (double)num_grid;
	    } 
	  else if(x[j] >= (double)num_grid)
	    {
	      x[j] -= (double)num_grid;
	    }

	  if(y[j] < 0.0)
	    {
	      y[j] += (double)num_grid;
	    }
	  else if(y[j] >= (double)num_grid)
	    {
	      y[j] -= (double)num_grid;
	    }

	  if(z[j] < 0.0)
	    {
	      z[j] += (double)num_grid;
	    }
	  else if(z[j] >= (double)num_grid)
	    {
	      z[j] -= (double)num_grid;
	    }

	  coords[0] = (int)(x[j]);
	  coords[1] = (int)(y[j]);
	  coords[2] = (int)(z[j]);

	  index = sfc_index( coords );
	  cart_assert( index >= 0 && index < num_root_cells );
          
	  /* check if we're supposed to read in this particle */
	  if(local_proc_id == processor_owner(index))
	    {
	      ipart = particle_alloc(id);
	      cart_assert(ipart>=0 && ipart<num_particles );

	      particle_x[ipart][0] = x[j];
	      particle_x[ipart][1] = y[j];
	      particle_x[ipart][2] = z[j];
	      particle_v[ipart][0] = vx[j];
	      particle_v[ipart][1] = vy[j];
	      particle_v[ipart][2] = vz[j];

	      particle_id[ipart] = id;
	      particle_mass[ipart] = particle_species_mass[species];
	      particle_level[ipart] = min_level;
	    }

	  id++;

	}
    }

  *idOffset = id;
}


void gicReadParticleData(char *rootname, char *type)
{
  const int lMax = 6;
  int l, L;
  long idOffset = 0;

  struct gicFile input[lMax];
  struct gicManifest manifest[lMax];
  struct gicFileHeader fileHeader[lMax];
  struct gicLevelHeader levelHeader[lMax];
  char filename[257];   /* GIC filenames are limited to 256 bytes */

  GIC_REAL *buffer[lMax];


  if(strlen(rootname) > 250)
    {
      cart_error("GIC filenames are limited to 256 bytes; rootnames to 250 bytes");
    }

  /*
  //  Loop over all open streams
  */
  for(l=0; l<lMax; l++)
    {
      strcpy(filename,rootname);
      strcat(filename,"_");
      strcat(filename,type);
      if(l < 3) 
	{
	  strcat(filename,".pos");
	}
      else
	{
	  strcat(filename,".vel");
	}

      gicStartFile(filename,input+l,manifest+l,fileHeader+l);

      if(l!=0 && (memcmp(manifest,manifest+l,sizeof(struct gicManifest))!=0 || memcmp(fileHeader,fileHeader+l,sizeof(struct gicFileHeader))!=0))
	{
	  cart_error("Incompatible input streams 0 and %d (%s)",l,filename);
	}
    }

  init_units();

  num_row = fileHeader->dims[0];
  num_particles_total = fileHeader->Ntot;
  num_particle_species = fileHeader->Lmax + 1;
  particle_species_indices[0] = 0;
  particle_species_indices[num_particle_species] = num_particles_total;

#ifdef STARFORM
  if(num_particle_species+1 > MAX_PARTICLE_SPECIES)
    {
      cart_error("header.Nspecies > MAX_PARTICLE_SPECIES.  Increase and rerun.");
    }
  particle_species_num[num_particle_species] = 0;
  particle_species_mass[num_particle_species] = 0.0;
  particle_species_indices[num_particle_species+1] = particle_species_indices[num_particle_species];

  num_particle_species++;

  total_stellar_mass = 0.0;
  total_stellar_initial_mass = 0.0;
#endif

  cart_debug("num_particle_species = %d",num_particle_species);
  cart_debug("num_particles_total = %d",num_particles_total);

  /*
  //  Create buffers
  */
  for(l=0; l<lMax; l++)
    {
      buffer[l] = (GIC_REAL *)cart_alloc(fileHeader->Nrec*sizeof(GIC_REAL));
    }

  /*
  //  Loop over all levels in REVERSE ORDER (!!!) since the highest level must be
  //  the first species for refining over it.
  */
  idOffset = 0L;
  for(L=fileHeader->Lmax; L>=0; L--)
    {
      gicReadParticleLevel(fileHeader->Lmax-L,L,lMax,input,fileHeader,levelHeader,&idOffset,buffer);
    }

  /*
  //  Delete buffers
  */
  for(l=0; l<lMax; l++)
    {
      cart_free(buffer[l]);
    }

  build_particle_list();
}

#endif  /* PARTICLES */


/*
//  *************************
//  *
//  *    GAS, alias HYDRO
//  *
//  *************************
*/
#ifdef HYDRO

float fracB, temIn, fracHII;


void gicAssignCellData(int cell, float w, float delta, float vx, float vy, float vz)
{
  int l;
  float d;
  int children[num_children];

  cart_assert(cell > -1);

  d = w*fracB*(1.0+delta);

  cell_gas_density(cell) += d;
  cell_momentum(cell,0) += d*vx;
  cell_momentum(cell,1) += d*vy;
  cell_momentum(cell,2) += d*vz;
  cell_gas_gamma(cell) = gamma;

  cell_gas_internal_energy(cell) += d*temIn/(gamma-1)*(1.0-Y_p+0.25*Y_p);

  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)*(gamma-1);
  cell_gas_energy(cell) = cell_gas_internal_energy(cell) + cell_gas_kinetic_energy(cell);

#ifdef RADIATIVE_TRANSFER
  cell_HI_density(cell) += d*rtXH*(1.0-fracHII);
  cell_HII_density(cell) += d*rtXH*fracHII;
  cell_HeI_density(cell) += d*rtXHe;
  cell_HeII_density(cell) += d*0.0;
  cell_HeIII_density(cell) += d*0.0;
  cell_H2_density(cell) += d*rtXH*2.0e-6;
#endif

  if(!cell_is_leaf(cell))
    {
      cart_assert(fabs(w-1.0)<1.0e-4);

      cell_all_children(cell,children);
      for(l=0; l<num_children; l++)
	{
	  gicAssignCellData(children[l],w,delta,vx,vy,vz);
	}
    }
}


void gicReadGasData(char *rootname, char *type)
{
  const int lMax = 4;
  int i, j, l, n, L, Lmax, ret;
  int save_ri0[num_refinement_levels+1];
  int children[num_children];
  float q, w, vFac;
  int page, num_pages;
  long noff, ntot, offset;
  double pos[3], fRef;
  int coords[3];
  int level, cell, num_level_cells;
  int *level_cells;

  struct gicFile input[lMax];
  struct gicManifest manifest[lMax];
  struct gicFileHeader fileHeader[lMax];
  struct gicLevelHeader levelHeader[lMax];
  char filename[257];   /* GIC filenames are limited to 256 bytes */

  GIC_REAL *buffer[lMax];
  GIC_REAL *delta, *vx, *vy, *vz;
  GIC_INTEGER *idx, *jdx, *kdx;
  int mask[num_root_cells];


  if(strlen(rootname) > 250)
    {
      cart_error("GIC filenames are limited to 256 bytes; rootnames to 250 bytes");
    }

  /*
  //  Loop over all open streams
  */
  for(l=0; l<lMax; l++)
    {
      strcpy(filename,rootname);
      strcat(filename,"_");
      strcat(filename,type);
      if(l > 0) 
	{
	  strcat(filename,".vel");
	}
      else
	{
	  strcat(filename,".den");
	}

      gicStartFile(filename,input+l,manifest+l,fileHeader+l);

      if(l!=0 && (memcmp(manifest,manifest+l,sizeof(struct gicManifest))!=0 || memcmp(fileHeader,fileHeader+l,sizeof(struct gicFileHeader))!=0))
	{
	  cart_error("Incompatible input streams 0 and %d (%s)",l,filename);
	}
    }

  if(min_level+fileHeader->Lmax < max_level)
    {
      Lmax = fileHeader->Lmax;
    }
  else
    {
      Lmax = max_level - min_level;
    }

  init_units();

  /*
  //  Thermal state of the primordial gas
  */
  fracB = Omegab0/Omega0;
  fracHII = 1.2e-5*sqrt(Omega0)/(Omegab0*hubble);
  q = fileHeader->aBegin*137.0*pow(Omegab0*hubble*hubble/0.022,0.4);
  temIn = T_CMB0/aexp[min_level]*q/pow(pow(q,1.73)+1,1.0/1.73)*wmu/T0;

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("f_HII: %e, T_in: %e",fracHII,temIn);
    }

  idx = (GIC_INTEGER *)cart_alloc(fileHeader->Nrec*sizeof(GIC_INTEGER));

  if(fileHeader->Lmax > 0)
    {
      /*
      //  Read the mask
      */
      num_pages = (num_root_cells+fileHeader->Nrec-1)/fileHeader->Nrec;
      noff = 0;
      for(page=0; page<num_pages; page++)
	{
	  n = fileHeader->Nrec;
	  if(noff+n > ntot) n = ntot - noff;
	  if(gicReadFortranRecordInteger(input,idx) != 0)
	    {
	      cart_error("File is corrupted, error in reading mask array on page %d",page);
	    }
	  for(j=0; j<n; j++)
	    {
	      mask[noff+j] = idx[j];
	    }
	  for(l=1; l<lMax; l++)
	    {
	      if(gicSkipFortranRecordInteger(input+l) != 0)
		{
		  cart_error("File is corrupted, error in skipping mask array in stream %d on page %d",l,page);
		}
	    }
	  noff += fileHeader->Nrec;
	}

      /*
      //  Refine by the mask
      */
      for(j=0; j<=num_refinement_levels; j++)
	{
	  save_ri0[j] = use_refinement_indicator[0][j];
	  use_refinement_indicator[0][j] = 1;
	}

      for(level=min_level; level<min_level+Lmax && level<max_level; level++)
	{
	  cart_debug("refining level %u",level);
	  select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
	  cart_debug("num_level_cells = %u", num_level_cells );
	  for(i=0; i<num_level_cells; i++)
	    {
	      cell_position_double(level_cells[i],pos);
	      cell = (int)pos[0] + num_grid*((int)pos[1]+num_grid*(int)pos[2]);
	      if(min_level+mask[cell] > level)
		{
		  refinement_indicator(level_cells[i],0) = 1.0;       
		}
	      else
		{
		  refinement_indicator(level_cells[i],0) = 0.0;
		}
	      /*
	      //  Zero cell variables - just in case
	      */
	      for(j=0; j<num_hydro_vars; j++)
		{
		  cell_var(level_cells[i],all_hydro_vars[j]) = 0.0;
		}
	      cell_gas_density(level_cells[i]) = 1.0e-30;  /* workaround a bug in cell energy refinement */
	    }
	  cart_free( level_cells );

	  cart_debug("about to refine level %u", level );
	  refine(level);
	}

      for(j=0; j<=num_refinement_levels; j++)
	{
	  use_refinement_indicator[0][j] = save_ri0[j];
	}
    }

  /*
  //  Create buffers
  */
  for(l=0; l<lMax; l++)
    {
      buffer[l] = (GIC_REAL *)cart_alloc(fileHeader->Nrec*sizeof(GIC_REAL));
    }
  jdx = (GIC_INTEGER *)cart_alloc(fileHeader->Nrec*sizeof(GIC_INTEGER));
  kdx = (GIC_INTEGER *)cart_alloc(fileHeader->Nrec*sizeof(GIC_INTEGER));

  delta = buffer[0];
  vx = buffer[1];
  vy = buffer[2];
  vz = buffer[3];

  vFac = fileHeader->aBegin/v0;

#ifdef DEBUG
  double dPos[3];
  cell_position_double(37212,dPos);
  cart_debug("That cell: %lf %lf %lf",dPos[0],dPos[1],dPos[2]);
#endif
  /*
  //  Now read the data
  */
  for(L=0; L<=Lmax; L++)
    {
      if(fileHeader->Lmax > 0)
	{
	  for(l=0; l<lMax; l++)
	    {
	      ret = gicReadLevelHeader(input+l,levelHeader+l);
	      if(ret != 0)
		{
		  cart_error("Unable to read GIC level header for input stream %d, level %d, error code: %d",l,L,ret);
		}
	      if(L!=levelHeader[l].L || fileHeader->Lmax!=levelHeader[l].Lmax)
		{
		  cart_error("Input stream %d is corrupted, error code 110",l);
		}
	      if(l == 0)
		{
		  if(levelHeader->ind == 0)
		    {
		      cart_error("Multi-resolution file must be indexed on level %d",L);
		    }
		  levelHeader->ind = 0; /* for comparison below */
		}
	      else if(memcmp(levelHeader,levelHeader+l,sizeof(struct gicLevelHeader))!=0 && levelHeader)
		{
		  cart_error("Incompatible input streams 0 and %d at level %d",l,L);
		}
	    }
	  levelHeader->ind = 1;
	  ntot = levelHeader->Nlev;
	}
      else
	{
	  ntot = num_root_cells;
	}

      num_pages = (ntot+fileHeader->Nrec-1)/fileHeader->Nrec;

      /*
      // Skip velocity streams to the appropriate record
      */
      for(l=2; l<lMax; l++)
	{
	  offset = ((l-1)%3)*num_pages*(2*sizeof(GIC_RECORD)+fileHeader->Nrec*sizeof(GIC_REAL));
	  ret = fseek(input[l].File,offset,SEEK_CUR);
	  if(ret != 0)
	    {
	      cart_error("Error in reposition stream %d, fseek error %d",l,ret);
	    }
	}

      fRef = pow(0.5,L);

      noff = 0;
      coords[0] = coords[1] = coords[2] = 0;
      for(page=0; page<num_pages; page++)
	{
	  n = fileHeader->Nrec;
	  if(noff+n > ntot) n = ntot - noff;
	  noff += fileHeader->Nrec;

	  for(l=0; l<lMax; l++)
	    {
	      if(gicReadFortranRecordReal(input+l,buffer[l]) != 0)
		{
		  cart_error("Input stream %d is corrupted, error code 120",l);
		}
	    }

	  /*
	  //  Read indecies if needed
	  */
	  if(levelHeader->ind != 0)
	    {
	      if(gicReadFortranRecordInteger(input,idx) != 0)
		{
		  cart_error("Input stream 0 is corrupted, error code 130");
		}
	      if(gicReadFortranRecordInteger(input,jdx) != 0)
		{
		  cart_error("Input stream 0 is corrupted, error code 140");
		}
	      if(gicReadFortranRecordInteger(input,kdx) != 0)
		{
		  cart_error("Input stream 0 is corrupted, error code 150");
		}
	    }
	  else
	    {
	      for(j=0; j<n; j++)
		{
		  idx[j] = coords[0];
		  jdx[j] = coords[1];
		  kdx[j] = coords[2];
                  coords[0]++;
                  if(coords[0] == num_grid)
                    {
                      coords[0] = 0;
                      coords[1]++;
                    }
                  if(coords[1] == num_grid)
                    {
                      coords[1] = 0;
                      coords[2]++;
                    }
		}
	    }

	  /*
	  //  Assign data now
	  */
	  for(j=0; j<n; j++)
	    {
	      /*
	      //  Units
	      */
	      vx[j] *= vFac;
	      vy[j] *= vFac;
	      vz[j] *= vFac;
	      
	      pos[0] = fRef*(idx[j]-0.5);
	      pos[1] = fRef*(jdx[j]-0.5);
	      pos[2] = fRef*(kdx[j]-0.5);

#ifdef DEBUG
	      if(fabs(pos[0]-dPos[0])<1.0e-6 && fabs(pos[1]-dPos[1])<1.0e-6 && fabs(pos[2]-dPos[2])<1.0e-6)
		{
		  cart_debug("Found that cell!");
		}
#endif

	      cell = cell_find_position_above_level(min_level+L,pos);
	      if(cell > -1)
		{
		  w = pow(0.125,min_level+L-cell_level(cell));
		  gicAssignCellData(cell,w,delta[j],vx[j],vy[j],vz[j]);
		}
	    }
	}
      /*
      // Skip velocity streams to the end of this level
      */
      for(l=1; l<lMax-1; l++)
	{
	  offset = ((3-l)%3)*num_pages*(2*sizeof(GIC_RECORD)+fileHeader->Nrec*sizeof(GIC_REAL));
	  ret = fseek(input[l].File,offset,SEEK_CUR);
	  if(ret != 0)
	    {
	      cart_error("Error in reposition stream %d, fseek error %d",l,ret);
	    }
	}
    }

  /*
  //  Delete buffers
  */
  for(l=0; l<lMax; l++)
    {
      cart_free(buffer[l]);
    }
  cart_free(idx);
  cart_free(jdx);
  cart_free(kdx);

  /*
  //  Fill in the holes left at higher levels
  */
  for(level=min_level+Lmax-1; level>=min_level; level--)
    {
      select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
      for(i=0; i<num_level_cells; i++) if(cell_gas_gamma(cell = level_cells[i]) < 1.0e-10)
	{
	  cell_all_children(cell,children);
	  for(j=0; j<num_hydro_vars; j++)
	    {
	      q = 0.0;
	      for(l=0; l<num_children; l++)
		{
		  q += cell_var(children[l],j);
		}
	      cell_var(cell,j) = q/num_children;
	    }
	}
      cart_free( level_cells );

      update_buffer_level(level,all_hydro_vars,num_hydro_vars);
    }
}

#endif  /* HYDRO */


void run_output()
{
#ifdef HYDRO

  const int nbin1 = 256;
#ifdef RADIATIVE_TRANSFER
  int varid[] = { EXT_FRACTION+RT_HVAR_OFFSET+0, HVAR_GAS_DENSITY, EXT_GAS_TEMPERATURE, EXT_FRACTION+RT_HVAR_OFFSET+5, EXT_CELL_LEVEL, EXT_LOCAL_PROC };
#else
  int varid[] = { HVAR_PRESSURE, HVAR_GAS_DENSITY, EXT_CELL_LEVEL };
#endif
  int nbin[] = { nbin1, nbin1, nbin1 };
  int nvars = sizeof(varid)/sizeof(int);
  double bb[6], pos[3], dbb;
  float dmax;
  char filename[99];

#ifdef RADIATIVE_TRANSFER
  rtSetTemUnits();
#endif

  bb[0] = 0.0;
  bb[1] = num_grid;
  bb[2] = 0.0;
  bb[3] = num_grid;
  bb[4] = 0.0;
  bb[5] = num_grid;

  sprintf(filename,"%s/out-box.%04d.bin",output_directory,(int)(aexp[0]*1.0e4));
  extWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);

  FindMaxVar(HVAR_GAS_DENSITY,&dmax,pos);

  dbb = max(1.0,nbin1*pow(0.5,(double)max_level));

  bb[0] = pos[0] - 0.5*dbb;
  bb[1] = pos[0] + 0.5*dbb;
  bb[2] = pos[1] - 0.5*dbb;
  bb[3] = pos[1] + 0.5*dbb;
  bb[4] = pos[2] - 0.5*dbb;
  bb[5] = pos[2] + 0.5*dbb;

  sprintf(filename,"%s/out-zoom.%04d.bin",output_directory,(int)(aexp[0]*1.0e4));
  extWriteIfritFile(max_level,nbin,bb,nvars,varid,filename);

#endif  /* HYDRO */
}


void init_run()
{
  int i, j;
  char *rootname, type[2];

  /*
  //  Where do we get the root name? Use options for now
  */
  if(local_proc_id == MASTER_NODE)
    {
      if(num_options < 1)
	{
          cart_error("An option -root=<name> is required, where <name> is the root name for a set of GIC input files.\n");
        }

      if(strncmp(options[0],"-root",5)!=0 || options[0][5]!='=' || strlen(options[0])<10)
	{
	  cart_error("-root option format is: -root=<name>");
	}
    }

  rootname = options[0] + 6;
  type[1] = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  /*
  //  Begin with load balancing 
  */
#ifdef HYDRO
  type[0] = 'B';
#else
  type[0] = 'M';
#endif
  gicBalanceLoad(rootname,type);


#ifdef PARTICLES

#ifdef HYDRO
  type[0] = 'D';
#else
  type[0] = 'M';
#endif

  gicReadParticleData(rootname,type);
  cart_debug("read in particles");
#endif

#ifdef HYDRO

  type[0] = 'B';
  gicReadGasData(rootname,type);
  cart_debug("read in gas");

  hydro_magic(min_level);
  hydro_eos(min_level);
#endif /* HYDRO */

  cart_debug("tl[min_level] = %f", tl[min_level] );
  cart_debug("aexp[min_level] = %f", aexp[min_level] );

  for(j=min_level+1; j<=max_level; j++)
    {
      tl[j] = tl[min_level];
      aexp[j] = aexp[min_level];
    }

  dtl[min_level] = 0.0;
  choose_timestep( &dtl[min_level] );

#ifdef PARTICLES
  for(j=0; j<num_particles; j++) if(particle_level[j] != FREE_PARTICLE_LEVEL)
    {
      particle_t[j] = tl[min_level];
      particle_dt[j] = dtl[min_level];
    }

  build_mesh();
#endif /* PARTICLES */

#ifdef STARFORM
  for ( i = 0; i < nDim; i++ )
    {
      star_formation_volume_min[i] = refinement_volume_min[i];
      star_formation_volume_max[i] = refinement_volume_max[i];
    }
#endif

  /*
  //  Debugging parameters
  */
#ifdef RADIATIVE_TRANSFER
#ifdef RT_DEBUG
  rt_debug.Mode = 0;
  rt_debug.Stop = 1;
  rt_debug.Pos[0] = 29.619681;
  rt_debug.Pos[1] = 14.352527;
  rt_debug.Pos[2] = 32.880561;
#endif
#endif
}
