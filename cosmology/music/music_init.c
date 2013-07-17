#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>

#include <mpi.h>

#include "config.h"
#include "io.h"
#include "io_artio.h"
#include "io_cart.h"
#include "tree.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "timing.h"
#include "times.h"
#include "auxiliary.h"
#include "cosmology.h"
#include "particle.h"

#include "hydro.h"
#include "iterators.h"
#include "load_balance.h"
#include "units.h"
#include "hydro.h"
#include "starformation.h"

#include "music_init.h"
#include "io_music.h"


#ifdef HYDRO
void error_check_music_io();
void zero_hydro( int icell ){
    cell_gas_density(icell) = 0;
    cell_momentum(icell,0) = 0;
    cell_momentum(icell,1) = 0;
    cell_momentum(icell,2) = 0;
    cell_gas_gamma(icell) = 0;
    cell_gas_internal_energy(icell) =  0;
    cell_gas_pressure(icell) =  0;
    cell_gas_energy(icell) = 0;
#ifdef ENRICHMENT
    cell_gas_metal_density_II(icell) = 0;
#ifdef ENRICHMENT_SNIa
    cell_gas_metal_density_Ia(icell) = 0;
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef RADIATIVE_TRANSFER
    cell_HI_density(icell) = 0;
    cell_HII_density(icell) = 0;
    cell_HeI_density(icell) = 0;
    cell_HeII_density(icell) = 0;
    cell_HeIII_density(icell) = 0;
    cell_H2_density(icell) = 0;
#endif
#ifdef EXTRA_PRESSURE_SOURCE
    cell_extra_pressure_source(icell) = 0;
#endif /* EXTRA_PRESSURE_SOURCE */
#ifdef ISOTROPIC_TURBULENCE_ENERGY
    cell_isotropic_turbulence_energy(icell) = 0;
#endif /* ISOTROPIC_TURBULENCE_ENERGY */
}

void set_zero_hydro() {
    int i, level, icell;
    int num_level_cells;
    int *level_cells;

    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, (CELL_TYPE_BUFFER || CELL_TYPE_LOCAL) | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            icell = level_cells[i] ;
            zero_hydro( icell );
	}
        cart_free( level_cells );
    }


    for ( level = max_level - 1; level >= min_level; level-- ) {
        hydro_split_update(level); /* update all non-leafs with their child's value */
        hydro_eos(level);
    }
}
#endif /* HYDRO */


void music_init() {
  char filename_header[257];  
  char filename_particles[257];
  char filename_hydro[257];  
  int i, level;
  const char *rootname, *tmp;
  char type[2];
  /*  Make sure we have a blank slate */
  if(cosmology_is_set()){
      cart_error("Cosmology is set before MUSIC files are read. Cosmological parameters should NOT be set in the config file when using the MUSIC reader.");
  }
  /*  Where do we get the root name? Use options for now */
  tmp = extract_option1("root-file-name","root",NULL);
  if(tmp != NULL){
      rootname = tmp;
  } else {
      cart_error("An option --root-file-name=<name> is required, where <name> is the root name for a set of MUSIC input files.");
  }
  /*  No more options are allowed. */
  if(num_options > 0){ cart_error("Unrecognized option: %s",options[0]);}
  MPI_Barrier(mpi.comm.run);


#ifdef HYDRO
  type[0] = 'H';
#else
  type[0] = 'D';
#endif
  type[1] = '\0';

#ifdef HYDRO
  /*
  //  If HYDRO is on, refinement now requires that the mesh contained
  //  the valid data. Hence, we fill the root cells with some dummy
  //  data.
  */
  for(i=0; i<num_root_cells; i++)
    {
      cell_gas_density(i) = 1e-10;
      cell_momentum(i,0) = 0;
      cell_momentum(i,1) = 0;
      cell_momentum(i,2) = 0;

      cell_gas_gamma(i) = constants->gamma;
      cell_gas_internal_energy(i) =  1;
      cell_gas_pressure(i) = cell_gas_internal_energy(i)*(constants->gamma-1);
      cell_gas_energy(i) = cell_gas_internal_energy(i);

#ifdef RADIATIVE_TRANSFER
      cell_HI_density(i) = cell_gas_density(i)*constants->XH;
      cell_HII_density(i) = 0;
      cell_HeI_density(i) = cell_gas_density(i)*constants->XHe;
      cell_HeII_density(i) = 0;
      cell_HeIII_density(i) = 0;
      cell_H2_density(i) = 0;
#endif
#ifdef EXTRA_PRESSURE_SOURCE
      cell_extra_pressure_source(i) = 0;
#endif /* EXTRA_PRESSURE_SOURCE */
#ifdef ISOTROPIC_TURBULENCE_ENERGY
      cell_isotropic_turbulence_energy(i) = 0;
#endif /* ISOTROPIC_TURBULENCE_ENERGY */
    }
#endif /* HYDRO */


#ifdef STARFORM
  for(i=0; i<nDim; i++)
    {
	star_formation_volume_min[i] = 0;
	star_formation_volume_max[i] = num_grid;
    }
#endif /* STARFORM */

      strcpy(filename_header,rootname);
      strcat(filename_header,"_");
      strcat(filename_header,type);
      strcat(filename_header,".mdh");

#ifdef PARTICLES
      strcpy(filename_particles,rootname);
      strcat(filename_particles,"_");
      strcat(filename_particles,type);
      strcat(filename_particles,".mdxv");
#endif

#ifdef HYDRO
      strcpy(filename_hydro,rootname);
      strcat(filename_hydro,"_");
      strcat(filename_hydro,type);
      strcat(filename_hydro,".md");
#endif


#ifdef PARTICLES
	cart_debug("load balance based on particles only");
	restart_load_balance_cart( NULL, filename_header, filename_particles );
#else
	cart_error("need particles for MUSIC ICs");
#endif /* PARTICLES */


#ifdef PARTICLES
	read_cart_particles( filename_header, filename_particles, NULL, NULL, 0, NULL );
#endif /* PARTICLES */

	read_cart_header_to_units(filename_header);

	units_update(min_level); /* aexpn only set on min_level */
	cosmology_set_fixed();
 	build_mesh();  

#ifdef HYDRO
	set_zero_hydro();
	read_music_gas_particles( filename_header, filename_hydro, 0, NULL ); 
#endif /* HYDRO */
	
	/* build_mesh only load balanced particles, now add gas */
  	load_balance(); 

#ifdef HYDRO
	for(level=max_level; level>=min_level; level--){
	    hydro_split_update(level);	    /* update non-leafs */
	}

	cart_debug("music io: error check");
	error_check_music_io();

	for(level=min_level; level<max_level; level++){
	    hydro_magic( level );
	    hydro_eos( level );
	    update_buffer_level(level, all_hydro_vars, num_hydro_vars);
	}	    
#endif /* HYDRO */

	check_map();

#ifdef HYDRO_TRACERS
	cart_debug("setting hydro tracers");
	set_hydro_tracers( min_level+1 );
#endif /* HYDRO_TRACERS */

	write_artio_restart(WRITE_SAVE,WRITE_SAVE,WRITE_SAVE);

	cart_debug("done reading data...");
	
}

#ifdef HYDRO
void error_check_music_io(){
    /* 
    // check that ART refinement doesn't over-resolve MUSIC 
    // refinement map (no empty cells after assign density) 
    // and make sure correct initial gas mass 
    */
    int num_level_cells, icell;
    int *level_cells;
    double sum=0, checksum=0;
    int level, i;
    double pos[nDim]; //snl
    for ( level = max_level; level >= min_level; level-- ) {
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for reduction(+ : sum) default(none), private(i,icell,pos), shared(num_level_cells,level_cells,cell_vars,level,constants)
	for ( i = 0; i < num_level_cells; i++ ) {
	    icell = level_cells[i] ;
	    sum+=cell_gas_density(icell)*cell_volume[level];

            if(cell_gas_gamma(icell)<=1){ //snl
                cell_gas_gamma(icell) = constants->gamma;
                hydro_magic_one_cell(icell);
#ifdef RADIATIVE_TRANSFER
                cell_HI_density(i) = cell_gas_density(i)*constants->XH;
                cell_HII_density(i) = 0;
                cell_HeI_density(i) = cell_gas_density(i)*constants->XHe;
                cell_HeII_density(i) = 0;
                cell_HeIII_density(i) = 0;
                cell_H2_density(i) = 0;
#endif
                cell_center_position(icell, pos);
                cart_error("location of unitialized cells %f %f %f %d %e",pos[0], pos[1], pos[2], level, cell_gas_density(icell));
            }
	}
	cart_free( level_cells );
    }

    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_DOUBLE, MPI_SUM, mpi.comm.run ); 

    if(local_proc_id==MASTER_NODE){
	checksum = sum * cosmology->OmegaM/cosmology->OmegaB / pow(num_grid,3.0) ;
	checksum = fabs(checksum)-1;
	if( checksum > 1.0e-6 ){
	    cart_error("baryon mass is wrong by %e",checksum);
	}
    }

    cart_debug("done checking hydro in music");
}
#endif /* HYDRO */


        
