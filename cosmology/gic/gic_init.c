#include "config.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "hydro.h"
#include "iterators.h"
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


#include "gic_reader.h"


extern int num_options;
extern char **options;


#if (!defined(PARTICLES) || !defined(COSMOLOGY))
#error "Both COSMOLOGY and PARTICLES must be set."
#endif


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
void gicBalanceLoad(const char *rootname, char *type)
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
      constrained_quantities = cart_alloc(int,num_constraints*num_root_cells);
      cell_work = cart_alloc(float,num_root_cells);

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
      num_particles_total = fileHeader.Ntot;

      if(fileHeader.Lmax == 0)
	{
	  /*
	  //  Single-resolution file - all work is uniform
	  */
	  nref = 1 << (nDim*L);
	  for(index=0; index<num_root_cells; index++)
	    {
	      constrained_quantities[num_constraints*index+0] = nref;
	      constrained_quantities[num_constraints*index+1] = nref;
	      cell_work[index] += (cost_per_cell+cost_per_particle)*nref;
	    }
	}
      else
	{
	  mask = cart_alloc(GIC_INTEGER,fileHeader.Nrec);
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
		  
		  nref = 1 << (nDim*mask[j]);

		  constrained_quantities[num_constraints*index+0] = nref;
		  constrained_quantities[num_constraints*index+1] = nref;
		  cell_work[index] += (cost_per_cell+cost_per_particle)*nref;

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
	  cart_free(mask);
	}

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

  /*
  //  All file operations done on the master node only
  */
  if(local_proc_id == MASTER_NODE)
    {
      gicSkipToLevel(Level,lMax,input,fileHeader,levelHeader,3,0);
    }

  id = *idOffset;

  /*
  //  Read level data (if needed)
  */
  if(fileHeader->Lmax > 0)
    {
      if(local_proc_id == MASTER_NODE)
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
	}

      MPI_Bcast(levelHeader,sizeof(struct gicLevelHeader),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);

      /*
      //  Set species arrays
      */
      particle_species_mass[species] = levelHeader->Mlev*constants->Msun/units->mass;
      ntot = levelHeader->Nlev;
    }
  else
    {
      /*
      //  Set species arrays
      */
#ifdef HYDRO
      particle_species_mass[species] = 1.0 - cosmology->OmegaB/cosmology->OmegaM;
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

  xFac = 1.0/units->length_in_chimps;
  /*
  //  GIC velocities are for auni, not abox
  */
  vFac = constants->kms*auni[min_level]/(units->velocity*abox[min_level]);

  /*
  // Skip to the appropriate record
  */
  if(local_proc_id == MASTER_NODE)
    {
      for(l=0; l<lMax; l++)
	{
	  offset = (l%3)*num_pages*(2*sizeof(GIC_RECORD)+fileHeader->Nrec*sizeof(GIC_REAL));
	  ret = fseek(input[l].File,offset,SEEK_CUR);
	  if(ret != 0)
	    {
	      cart_error("Error in reposition stream %d, fseek error %d",l,ret);
	    }
	}
    }

  noff = 0;
  for(page=0; page<num_pages; page++)
    {
      n = fileHeader->Nrec;
      if(noff+n > ntot) n = ntot - noff;
      noff += fileHeader->Nrec;

      if(local_proc_id == MASTER_NODE)
	{
	  for(l=0; l<lMax; l++)
	    {
	      if(gicReadFortranRecordReal(input+l,buffer[l]) != 0)
		{
		  cart_error("Input stream %d is corrupted, error code 020",l);
		}
	    }
	}

      for(l=0; l<lMax; l++)
	{
	  MPI_Bcast(buffer[l],fileHeader->Nrec*sizeof(GIC_REAL),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);
	}

      /*
      //  We need a barrier here to avoid overfilling MPI buffers 
      //  with too many asynchronized broadcasts
      */
      if(page%100 == 99) MPI_Barrier(MPI_COMM_WORLD);

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
	      particle_t[ipart] = particle_dt[ipart] = 0.0; /* used in load_balance, set to something to avoid junk in the buffers */
	    }

	  id++;

	}
    }

  *idOffset = id;
}


void gicReadParticleData(const char *rootname, char *type, int dc_off)
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
  if(local_proc_id == MASTER_NODE)
    {
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

      cart_debug("Number of data elements to load: %ld",fileHeader->Ntot);
    }

  MPI_Bcast(manifest,sizeof(struct gicManifest),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);
  MPI_Bcast(fileHeader,sizeof(struct gicFileHeader),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);

  /*
  //  Set units
  */
  cosmology_set(OmegaM,manifest->OmegaX+manifest->OmegaB);
  cosmology_set(OmegaB,manifest->OmegaB);
  cosmology_set(OmegaL,manifest->OmegaL);
  cosmology_set(h,manifest->h100);

  if(dc_off)
    cosmology_set(DeltaDC,0.0);
  else
    cosmology_set(DeltaDC,fileHeader->DeltaDC);

  box_size = manifest->dx*num_grid;

  auni[min_level] = fileHeader->aBegin;
  tl[min_level] = tcode_from_auni(auni[min_level]);
  abox[min_level] = abox_from_auni(auni[min_level]);

  units_reset();
  units_update(min_level);

  num_row = 4*fileHeader->dims[0];  /* Empirical number to reduce the # of communications */
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
      buffer[l] = cart_alloc(GIC_REAL,fileHeader->Nrec);
    }

  /*
  //  Loop over all levels in REVERSE ORDER (!!!) since the highest level must be
  //  the first species for refining over it.
  */
  idOffset = 0L;
  for(L=fileHeader->Lmax; L>=0; L--)
    {
      gicReadParticleLevel(fileHeader->Lmax-L,L,lMax,input,fileHeader,levelHeader,&idOffset,buffer);
      if(local_proc_id == MASTER_NODE)
	{
	  cart_debug("Finished reading level %d",L);
	}
    }

  /*
  //  Delete buffers
  */
  for(l=0; l<lMax; l++)
    {
      cart_free(buffer[l]);
    }

  if(local_proc_id == MASTER_NODE)
    {
      for(l=0; l<lMax; l++)
	{
	  fclose(input[l].File);
	}
    }

  build_particle_list();
  build_mesh();
}


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

  cell_gas_gamma(cell) = constants->gamma;
  cell_gas_internal_energy(cell) =  cell_gas_density(cell)*temIn/(constants->gamma-1)*(1.0-constants->Yp+0.25*constants->Yp);
  cell_gas_pressure(cell) = cell_gas_internal_energy(cell)*(constants->gamma-1);
  cell_gas_energy(cell) = cell_gas_internal_energy(cell) + cell_gas_kinetic_energy(cell);

#ifdef RADIATIVE_TRANSFER
  cell_HI_density(cell) = cell_gas_density(cell)*constants->XH*(1.0-fracHII);
  cell_HII_density(cell) = cell_gas_density(cell)*constants->XH*fracHII;
  cell_HeI_density(cell) = cell_gas_density(cell)*constants->XHe;
  cell_HeII_density(cell) = cell_gas_density(cell)*0.0;
  cell_HeIII_density(cell) = cell_gas_density(cell)*0.0;
  cell_H2_density(cell) = cell_gas_density(cell)*constants->XH*2.0e-6;
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


void gicReadGasData(const char *rootname, char *type)
{
  const int lMax = 4;
  int i, j, l, n, L, Lmax, ret;
  int children[num_children];
  float q, w, vFac, dMax, dMin, delMax, delMin;
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


  if(strlen(rootname) > 250)
    {
      cart_error("GIC filenames are limited to 256 bytes; rootnames to 250 bytes");
    }

  /*
  //  Loop over all open streams
  */
  if(local_proc_id == MASTER_NODE)
    {
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
    }

  MPI_Bcast(manifest,sizeof(struct gicManifest),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);
  MPI_Bcast(fileHeader,sizeof(struct gicFileHeader),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);

  if(min_level+fileHeader->Lmax < max_level)
    {
      Lmax = fileHeader->Lmax;
    }
  else
    {
      Lmax = max_level - min_level;
    }

  /*
  //  Thermal state of the primordial gas
  */
  fracB = cosmology->OmegaB/cosmology->OmegaM;
  fracHII = 1.2e-5*sqrt(cosmology->Omh2)/cosmology->Obh2;
  q = auni[min_level]*137.0*pow(cosmology->Obh2/0.022,0.4);
  temIn = 2.728/auni[min_level]*q/pow(pow(q,1.73)+1,1.0/1.73);

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Initial temperature: %f",temIn);
    }

  temIn /= units->temperature;

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("f_HII: %e, T_in: %e",fracHII,temIn);
    }

  /*
  //  Skip the mask
  */
  if(local_proc_id == MASTER_NODE)
    {
      if(fileHeader->Lmax > 0)
	{
	  num_pages = (num_root_cells+fileHeader->Nrec-1)/fileHeader->Nrec;
	  for(page=0; page<num_pages; page++)
	    {
	      for(l=0; l<lMax; l++)
		{
		  if(gicSkipFortranRecordInteger(input+l) != 0)
		    {
		      cart_error("File is corrupted, error in skipping mask array in stream %d on page %d",l,page);
		    }
		}
	    }
	}
    }

  /*
  //  Create buffers
  */
  for(l=0; l<lMax; l++)
    {
      buffer[l] = cart_alloc(GIC_REAL,fileHeader->Nrec);
    }
  idx = cart_alloc(GIC_INTEGER,fileHeader->Nrec);
  jdx = cart_alloc(GIC_INTEGER,fileHeader->Nrec);
  kdx = cart_alloc(GIC_INTEGER,fileHeader->Nrec);

  delta = buffer[0];
  vx = buffer[1];
  vy = buffer[2];
  vz = buffer[3];

  /*
  //  GIC velocities are for auni, not abox
  */
  vFac = constants->kms*auni[min_level]/(units->velocity*abox[min_level]);

  dMin =  1.0e35;
  dMax = -1.0e35;

  /*
  //  Now read the data
  */
  for(L=0; L<=Lmax; L++)
    {
      if(fileHeader->Lmax > 0)
	{
	  if(local_proc_id == MASTER_NODE)
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
	    }

	  MPI_Bcast(levelHeader,sizeof(struct gicLevelHeader),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);

	  ntot = levelHeader->Nlev;
	}
      else
	{
	  ntot = num_root_cells;
	  levelHeader->ind = 0;
	}

      num_pages = (ntot+fileHeader->Nrec-1)/fileHeader->Nrec;

      /*
      // Skip velocity streams to the appropriate record
      */
      if(local_proc_id == MASTER_NODE)
	{
	  for(l=2; l<lMax; l++)
	    {
	      offset = ((l-1)%3)*num_pages*(2*sizeof(GIC_RECORD)+fileHeader->Nrec*sizeof(GIC_REAL));
	      ret = fseek(input[l].File,offset,SEEK_CUR);
	      if(ret != 0)
		{
		  cart_error("Error in reposition stream %d, fseek error %d",l,ret);
		}
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

	  if(local_proc_id == MASTER_NODE)
	    {
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
		      idx[j] = 1 + coords[0];
		      jdx[j] = 1 + coords[1];
		      kdx[j] = 1 + coords[2];
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

	  for(l=0; l<lMax; l++)
	    {
	      MPI_Bcast(buffer[l],fileHeader->Nrec*sizeof(GIC_REAL),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);
	    }
	  MPI_Bcast(idx,fileHeader->Nrec*sizeof(GIC_INTEGER),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);
	  MPI_Bcast(jdx,fileHeader->Nrec*sizeof(GIC_INTEGER),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);
	  MPI_Bcast(kdx,fileHeader->Nrec*sizeof(GIC_INTEGER),MPI_BYTE,MASTER_NODE,MPI_COMM_WORLD);

	  /*
	  //  We need a barrier here to avoid overfilling MPI buffers 
	  //  with too many asynchronized broadcasts
	  */
	  if(page%100 == 99) MPI_Barrier(MPI_COMM_WORLD);

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

	      cell = cell_find_position_above_level(min_level+L,pos);
	      if(cell > -1)
		{
		  if(dMin > delta[j]) dMin = delta[j];
		  if(dMax < delta[j]) dMax = delta[j];
		  w = pow(0.125,min_level+L-cell_level(cell));
		  gicAssignCellData(cell,w,delta[j],vx[j],vy[j],vz[j]);
		}
	    }
	}
      /*
      // Skip velocity streams to the end of this level
      */
      if(local_proc_id == MASTER_NODE)
	{
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

  MPI_Allreduce(&dMax,&delMax,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(&dMin,&delMin,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);

  if(local_proc_id == MASTER_NODE)
    {
      cart_debug("Overdensity range: %f %f",delMin,delMax);
      cart_debug("Gas density range: %f %f",fracB*(1+delMin),fracB*(1+delMax));
    }

  if(local_proc_id == MASTER_NODE)
    {
      for(l=0; l<lMax; l++)
        {
          fclose(input[l].File);
        }
    }

  /*
  //  Fill in the holes left at higher levels
  */
  for(level=min_level+Lmax-1; level>=min_level; level--)
    {
      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
      for(i=0; i<num_level_cells; i++) if(cell_gas_gamma(cell = level_cells[i]) < 1.0e-10)
	{
	  cell_all_children(cell,children);
	  for(j=0; j<num_hydro_vars; j++)
	    {
	      q = 0.0;
	      for(l=0; l<num_children; l++)
		{
		  q += cell_var(children[l],all_hydro_vars[j]);
		}
	      cell_var(cell,all_hydro_vars[j]) = q/num_children;
	    }
	}
      cart_free(level_cells);
    }

  /*
  //  Fill in the holes left at lower levels
  */
  for(level=min_level+1; level<=min_level+Lmax; level++)
    {
      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
      for(i=0; i<num_level_cells; i++) if(cell_gas_gamma(cell = level_cells[i]) < 1.0e-10)
	{
	  for(j=0; j<num_hydro_vars; j++)
	    {
	      cell_var(cell,all_hydro_vars[j]) = cell_var(cell_parent_cell(cell),all_hydro_vars[j]);
	    }
	}
      cart_free(level_cells);
    }

  /*
  //  Update the buffer everywhere
  */
  for(level=min_level; level<=max_level; level++)
    {
      update_buffer_level(level,all_hydro_vars,num_hydro_vars);
    }

#ifdef DEBUG
  for(level=min_level; level<=max_level; level++)
    {
      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
      for(i=0; i<num_level_cells; i++)
	{
	  cell = level_cells[i];

	  if(cell_gas_gamma(cell) < 1.0e-10)
	    {
	      cart_debug("Error in cell %d#%d",cell,cell_level(cell));
	      cell_center_position(cell,pos);
	      cart_debug("Position: %lf %lf %lf",pos[0],pos[1],pos[2]);
	      for(i=0; i<num_vars; i++)
		{
		  cart_debug("Var[%d] = %g",i,cell_var(cell,i));
		}
	      cart_error("Not all vars are filled!!!");
	    }

	  if(cell_is_refined(cell))
	    {
	      cell_all_children(cell,children);

	      for(j=0; j<num_hydro_vars; j++)
		{
		  q = 0.0;
		  for(l=0; l<num_children; l++)
		    {
		      q += cell_var(children[l],all_hydro_vars[j]);
		    }
		  q /= num_children;
		  if(fabs(cell_var(cell,all_hydro_vars[j])-q)/(1.0e-35+fabs(cell_var(cell,all_hydro_vars[j]))+fabs(q)) > 1.0e-4)
		    {
		      cart_debug("Error in cell %d#%d",cell,cell_level(cell));
		      cell_center_position(cell,pos);
		      cart_debug("Position: %lf %lf %lf",pos[0],pos[1],pos[2]);
		      for(l=0; l<num_children; l++)
			{
			  cart_debug("Child %d: %d",l,children[l]);
			}
		      cart_debug("Var[%d] = %g",all_hydro_vars[j],cell_var(cell,all_hydro_vars[j]));
		      cart_debug("Sum over children: = %g",q);
		      cart_error("Conservation is broken!!!");
		    }
		}
	    }
	}
      cart_free(level_cells);
    }
#endif
}

#endif  /* HYDRO */


void gic_init()
{
  int i, level;
  const char *tmp, *rootname;
  char type[2];
  int dc_off = 0;

  /*
  //  Make sure we have a blank slate
  */
  if(cosmology_is_set())
    {
      cart_error("Cosmology is set before GIC files are read. Cosmological parameters should NOT be set in the config file when using the GIC reader.");
    }
  
  /*
  //  Where do we get the root name? Use options for now
  */
  if(local_proc_id == MASTER_NODE)
    {
      if(num_options < 1)
	{
          cart_error("An option --root=<name> is required, where <name> is the root name for a set of GIC input files.\n");
        }
    }

  for(i=0; i<num_options; i++)
    {
      /*
      //  Root name for data files
      */
      tmp = check_option1(options[i],"-root",NULL);
      if(tmp != NULL)
	{
	  rootname = tmp;
	  continue;
	}

      /*
      //  Switch off the DC mode
      */
      tmp = check_option0(options[i],"-no-dc");
      if(tmp != NULL)
	{
	  dc_off = 1;
	  continue;
	}

      cart_error("Unrecognized option: %s",options[i]);
    }

  type[1] = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  /*
  //  Begin with load balancing 
  */
#ifdef HYDRO
  type[0] = 'D';
#else
  type[0] = 'M';
#endif

  gicBalanceLoad(rootname,type);
  gicReadParticleData(rootname,type,dc_off);
  cart_debug("read in particles");

#ifdef HYDRO

  type[0] = 'B';
  gicReadGasData(rootname,type);
  cart_debug("read in gas");

  hydro_magic(min_level);
  hydro_eos(min_level);

#endif /* HYDRO */

  cart_debug("tl[min_level] = %f", tl[min_level] );
  cart_debug("au[min_level] = %f", auni[min_level] );
  cart_debug("ab[min_level] = %f", abox[min_level] );
  cart_debug("DC mode = %f", cosmology->DeltaDC );

  cosmology_set_fixed();

  for(level=min_level+1; level<=max_level; level++)
    {
      tl[level] = tl[min_level];
      auni[level] = auni[min_level];
      abox[level] = abox[min_level];
    }

  for(i=0; i<num_particles; i++) if(particle_level[i] != FREE_PARTICLE_LEVEL)
    {
      particle_t[i] = tl[min_level];
      /*
      //  We set the step to 0 so that the first leapfrog step is correct
      */
      particle_dt[i] = 0.0;
    }

#ifdef STARFORM
  for(i=0; i<nDim; i++)
    {
      star_formation_volume_min[i] = refinement_volume_min[i];
      star_formation_volume_max[i] = refinement_volume_max[i];
    }
#endif
}

