#include "config.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "hydro.h"
#include "io.h"
#include "io_cart.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "sfc.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "../run/hydro_step.h"

#include "../gic/gic_reader.h"


#if (!defined(PARTICLES) || !defined(COSMOLOGY) || !defined(HYDRO))
#error "All of COSMOLOGY, PARTICLES, and HYDRO must be set."
#endif


void cosmics_init()
{
  static const int page_size = 262144;
  FILE *f[6];
  int i, j, n, index, ipart, coords[3];
  int level, levelMax;
  int l, wrong_order[6], page, num_pages;
  int cell, num_level_cells, *level_cells;
  int slice, num_slices, num_data_per_slice, num_data_done;
  long id;
  int ng1, ng2;
  double x0[3], x[3], fRef;
  float q, xFac, vFac;
  float *buffer[6], *vxc, *vyc, *vzc, *vxb, *vyb, *vzb;
  float fracB, temIn, fracHII;
  int children[num_children];

  GIC_RECORD s1, s2;
  const char *tmp, *dir;
  char filename[999];
  struct cosmics_header
  {
    int n[3];
    float dx, abeg, OmegaM, OmegaL, H0;
  } header[6];

  
  /*
  //  Where do we get the root name? Use options for now
  */
  tmp = extract_option1("dir","dir",NULL);
  if(tmp != NULL)
    {
      dir = tmp;
    }
  else
    {
      cart_error("An option --dir=<name> is required, where <name> is the directory name for a set of COSMICS input files.");
    }
  
  /*
  //  No more options are allowed.
  */
  if(num_options > 0)
    {
      cart_error("Unrecognized option: %s",options[0]);
    }

  MPI_Barrier(mpi.comm.run);

  if(local_proc_id == MASTER_NODE)
    {
      for(l=0; l<6; l++)
	{
	  strcpy(filename,dir);
	  switch(l)
	    {
	    case 0:
	      {
		strcat(filename,"/ic_vxc.dat");
		break;
	      }
	    case 1:
	      {
		strcat(filename,"/ic_vyc.dat");
		break;
	      }
	    case 2:
	      {
		strcat(filename,"/ic_vzc.dat");
		break;
	      }
	    case 3:
	      {
		strcat(filename,"/ic_vxb.dat");
		break;
	      }
	    case 4:
	      {
		strcat(filename,"/ic_vyb.dat");
		break;
	      }
	    case 5:
	      {
		strcat(filename,"/ic_vzb.dat");
		break;
	      }
	    }

	  f[l] = fopen(filename,"r");
	  cart_assert(f[l] != NULL);

	  if(gicReadRecordHelper(f[l],sizeof(struct cosmics_header),header+l,wrong_order+l) != 0)
	    {
	      cart_error("Error in reading the header for stream %d, file %s",l,filename);
	    }
	  if(l!=0 && memcmp(header,header+l,sizeof(struct cosmics_header))!=0)
	    {
	      cart_error("Incompatible input streams 0 and %d",l);
	    }
	}

      if(wrong_order[0])
	{
	  reorder((char*)&header->n[0],sizeof(int));
	  reorder((char*)&header->n[1],sizeof(int));
	  reorder((char*)&header->n[2],sizeof(int));
	  reorder((char*)&header->dx,sizeof(float));
	  reorder((char*)&header->abeg,sizeof(float));
	  reorder((char*)&header->OmegaM,sizeof(float));
	  reorder((char*)&header->OmegaL,sizeof(float));
	  reorder((char*)&header->H0,sizeof(float));
	}


      if(header->n[0]!=header->n[1] || header->n[1]!=header->n[2])
	{
	  cart_error("Only a cubic input mesh is supported.");
	}
    }

  MPI_Bcast(header,sizeof(struct cosmics_header),MPI_BYTE,MASTER_NODE,mpi.comm.run);

  levelMax = 0;
  while(header->n[0] > num_grid)
    {
      header->n[0] = header->n[0] >> 1;
      levelMax++;
    }

  if(num_grid != header->n[0])
    {
      cart_error("The input grid size (=%d) is not a power-of-two multiple of num_grid (=%d).",header->n[1],num_grid);
    }

  cart_assert(header->n[1] == num_grid << levelMax);
  levelMax += min_level;

  /*
  //  Set units
  */
  cosmology_set(OmegaM,header->OmegaM);
  cosmology_set(OmegaB,0.04);
  cosmology_set(OmegaL,header->OmegaL);
  cosmology_set(h,header->H0/100.0);
  cosmology_set(DeltaDC,0.0);
  box_size = header->dx*cosmology->h*num_grid;

  auni[min_level] = header->abeg;
  tl[min_level] = tcode_from_auni(auni[min_level]);
  abox[min_level] = abox_from_auni(auni[min_level]);

  units_init();
  units_update(min_level);

  /*
  //  Particle parameters
  */
  num_particles_total = (particleid_t)header->n[1]*(particleid_t)header->n[1]*(particleid_t)header->n[1];
  num_particle_species = 1;
  particle_species_num[0] = num_particles_total;
  particle_species_mass[0] = 1.0 - cosmology->OmegaB/cosmology->OmegaM;
  particle_species_indices[0] = 0;
  particle_species_indices[1] = num_particles_total;

#ifdef STARFORM
  if(MAX_PARTICLE_SPECIES < 2)
    {
      cart_error("MAX_PARTICLE_SPECIES should be at least 2. Increase and rerun.");
    }

  num_particle_species = 2;

  particle_species_num[1] = 0;
  particle_species_mass[1] = 0.0;
  particle_species_indices[2] = particle_species_indices[1];

  total_stellar_mass = 0.0;
  total_stellar_initial_mass = 0.0;
#endif

  cart_debug("num_particle_species = %d",num_particle_species);
  cart_debug("num_particles_total = %d",num_particles_total);


  /*
  //  Balance load - split uniformly
  */
  for(i=0; i<=num_procs; i++)
    {
      proc_sfc_index[i] = ((unsigned long)num_root_cells*(unsigned long)i)/num_procs;
    }
  init_tree();

  for(i=0; i<nDim; i++)
    {
      refinement_volume_min[i] = 0.0;
      refinement_volume_max[i] = num_grid;
    }

  /*
  // Refine grid uniformly to levelMax
  */
  for(level=min_level; level<levelMax; level++)
    {
      cart_debug("refining level %d",level);

      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
      cart_debug("num_level_cells = %d",num_level_cells);
      
      for(i=0; i<num_level_cells; i++)
	{
	  refinement_indicator(level_cells[i],0) = 1.0;       
	}
       cart_free( level_cells );

       refine(level);
    }

  /*
  //  Read in the data
  */
  for(l=0; l<6; l++)
    {
      buffer[l] = cart_alloc(float,page_size);
    }

  vxc = buffer[0];
  vyc = buffer[1];
  vzc = buffer[2];
  vxb = buffer[3];
  vyb = buffer[4];
  vzb = buffer[5];

  /*
  //  Unit conversion factors
  */
  vFac = constants->kms/units->velocity;
  xFac = abox[min_level]*abox[min_level]*constants->Mpc/(100*cosmology->h*units->length)*dPlus(abox[min_level])/qPlus(abox[min_level]);

  if(header->n[1] > 256)
    {
      num_slices = header->n[1];
      num_data_per_slice = header->n[1]*header->n[1];
    }
  else
    {
      num_slices = 1;
      num_data_per_slice = header->n[1]*header->n[1]*header->n[1];
    }

  num_pages = (num_data_per_slice+page_size-1)/page_size;

  id = 0L;

  fRef = pow(0.5,levelMax);
  ng1 = num_grid << levelMax;
  ng2 = ng1*ng1;

  for(slice=0; slice<num_slices; slice++)
    {

      num_data_done = 0;

      if(local_proc_id == MASTER_NODE)
	{
	  for(l=0; l<6; l++)
	    {
	      if(fread(&s1,sizeof(GIC_RECORD),1,f[l]) != 1)
		{
		  cart_error("Error in reading header for file %d, record %d",l,slice);
		}
	      if(wrong_order[l]) reorder((char *)&s1,sizeof(s1));
	      if(s1 != sizeof(float)*num_data_per_slice)
		{
		  cart_error("Header for file %d, record %d is corrupted: %d, should be %d",l,slice,s1,num_data_per_slice);
		}
	    }
	}

      for(page=0; page<num_pages; page++)
	{

	  n = page_size;
	  if(num_data_done+n > num_data_per_slice)
	    {
	      n = num_data_per_slice - num_data_done;
	      cart_assert(page == (num_pages-1));
	    }
	  num_data_done += n;

	  if(local_proc_id == MASTER_NODE)
	    {
	      for(l=0; l<6; l++)
		{
		  if(fread(buffer[l],sizeof(float),n,f[l]) != n)
		    {
		      cart_error("Error in reading data for file %d, record %d, page %d",l,slice,page);
		    }
		  if(wrong_order[l])
		    {
		      for(j=0; j<n; j++) reorder((char *)(buffer[l]+j),sizeof(float));
		    }
		}
	    }

	  for(l=0; l<6; l++)
	    {
	      MPI_Bcast(buffer[l],n,MPI_FLOAT,MASTER_NODE,mpi.comm.run);
	    }

	  /*
	  //  We need a barrier here to avoid overfilling MPI buffers 
	  //  with too many asynchronized broadcasts
	  */
	  if(page%100 == 99) MPI_Barrier(mpi.comm.run);

	  for(j=0; j<n; j++)
	    {

	      /*
	      //  Particle position
	      */
	      x0[0] = fRef*(0.5+(id % ng1));
	      x0[1] = fRef*(0.5+(id/ng1 % ng1));
	      x0[2] = fRef*(0.5+(id/ng2 % ng1));
	      
	      x[0] = xFac*vxc[j] + x0[0];
	      x[1] = xFac*vyc[j] + x0[1];
	      x[2] = xFac*vzc[j] + x0[2];

	      /* enforce periodic boundary conditions */
	      for(i=0; i<3; i++)
		{
		  if(x[i] < 0.0)
		    {
		      x[i] += (double)num_grid;
		    } 
		  else if(x[i] >= (double)num_grid)
		    {
		      x[i] -= (double)num_grid;
		    }
		  coords[i] = (int)(x[i]);
		}

	      index = sfc_index( coords );
	      cart_assert( index >= 0 && index < num_root_cells );
          
	      /* check if we're supposed to read in this particle */
	      if(local_proc_id == processor_owner(index))
		{
		  ipart = particle_alloc(id);
		  cart_assert(ipart>=0 && ipart<num_particles );

		  particle_x[ipart][0] = x[0];
		  particle_x[ipart][1] = x[1];
		  particle_x[ipart][2] = x[2];
		  particle_v[ipart][0] = vFac*vxc[j];
		  particle_v[ipart][1] = vFac*vyc[j];
		  particle_v[ipart][2] = vFac*vzc[j];

		  particle_id[ipart] = id;
		  particle_mass[ipart] = particle_species_mass[0];
		  particle_level[ipart] = min_level + levelMax;
		}

	      for(i=0; i<3; i++)
		{
		  coords[i] = (int)(x0[i]);
		}

	      index = sfc_index( coords );
	      cart_assert( index >= 0 && index < num_root_cells );
	      if(local_proc_id == processor_owner(index))
		{
		  cell = cell_find_position(x0);
#ifdef DEBUG
		  if(cell == -1)
		    {
		      cart_debug("%lf %lf %lf",x0[0],x0[1],x0[2]);
		      cart_debug("%ld %d %g",id,ng1,fRef);
		    }
#endif
		  cart_assert(cell != -1);

		  cell_var(cell,HVAR_MOMENTUM+0) = vFac*vxb[j];
		  cell_var(cell,HVAR_MOMENTUM+1) = vFac*vyb[j];
		  cell_var(cell,HVAR_MOMENTUM+2) = vFac*vzb[j];
		}

	      id++;
	    }
	}

      if(local_proc_id == MASTER_NODE)
	{
	  for(l=0; l<6; l++)
	    {
	      if(fread(&s2,sizeof(GIC_RECORD),1,f[l]) != 1)
		{
		  cart_error("Error in reading footer for file %d, record %d",l,slice);
		}
	      if(wrong_order[l]) reorder((char *)&s2,sizeof(s2));
	      if(s2 != sizeof(float)*num_data_per_slice)
		{
		  cart_error("Footer for file %d, record %d is corrupted: %d, should be %d",l,slice,s2,num_data_per_slice);
		}
	    }
	}
    }

  if(local_proc_id == MASTER_NODE)
    {
      for(l=0; l<6; l++) fclose(f[l]);
    }

  
  for(l=0; l<6; l++) cart_free(buffer[l]);
  
  build_particle_list();

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
  //  Finish filling in the lowest level
  */
  select_level(min_level+levelMax,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
  for(i=0; i<num_level_cells; i++)
    {
      cell = level_cells[i];

      cell_gas_density(cell) = fracB;
      cell_momentum(cell,0) *= fracB;
      cell_momentum(cell,1) *= fracB;
      cell_momentum(cell,2) *= fracB;

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
#ifdef EXTRA_PRESSURE_SOURCE
      cell_extra_pressure_source(cell) = 0;
#endif /* EXTRA_PRESSURE_SOURCE */
#ifdef ISOTROPIC_TURBULENCE_ENERGY
      cell_isotropic_turbulence_energy(cell) = 0;
#endif /* ISOTROPIC_TURBULENCE_ENERGY */
    }
  cart_free(level_cells);


  /*
  //  Finish filling in the grid
  */
  for(level=min_level+levelMax-1; level>=min_level; level--)
    {
      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);
      for(i=0; i<num_level_cells; i++)
        {
	  cell = level_cells[i];

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

  build_cell_buffer();
  repair_neighbors();

  /*
  //  Update the buffer everywhere
  */
  for(level=min_level; level<=max_level; level++)
    {
      update_buffer_level(level,all_hydro_vars,num_hydro_vars);
    }

  hydro_magic(min_level);
  hydro_eos(min_level);

  cart_debug("tl[min_level] = %f", tl[min_level] );
  cart_debug("au[min_level] = %f", auni[min_level] );
  cart_debug("ab[min_level] = %f", abox[min_level] );

  for(level=min_level+1; level<=max_level; level++)
    {
      tl[level] = tl[min_level];
      auni[level] = auni[min_level];
      abox[level] = abox[min_level];
    }

  for(i=0; i<num_particles; i++) if(particle_level[i] != FREE_PARTICLE_LEVEL)
    {
      particle_t[i] = tl[min_level];
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








