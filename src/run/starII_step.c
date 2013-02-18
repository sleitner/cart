#include "config.h"
#ifdef STAR_FORMATION
#ifdef STAR_PARTICLE_TYPES

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "iterators.h" 
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "io.h"
#include "particle.h"
#include "tree.h"
#include "starformation.h"

#include "units.h"
#include "onestarfits.h"

void starII_delete_snII(int level){
    int i;
    int ipart, ipart_prev;
    int iter_cell;
    int num_level_cells;
    int *level_cells;
    double ini_mass_sol, Zsol;
    float star_age;

    select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
    
#pragma omp parallel for default(none), private(iter_cell,ipart,star_age,ini_mass_sol,Zsol,ipart_prev), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,star_tbirth,particle_id,star_particle_type,particle_species_indices,num_particle_species,particle_list_next, particle_list_prev, star_initial_mass, star_metallicity_II, units, constants), schedule(dynamic)
    for ( i = 0; i < num_level_cells; i++ ) {
	iter_cell = level_cells[i];
	
	ipart = cell_particle_list[iter_cell];
	    while ( ipart != NULL_PARTICLE ) {
		star_age = particle_t[ipart] - star_tbirth[ipart] ;
		ini_mass_sol = star_initial_mass[ipart]*units->mass/constants->Msun;
		Zsol = star_metallicity_II[ipart]/constants->Zsun;
		ipart_prev = particle_list_prev[ipart] ;

		if ( particle_is_star(ipart) &&
		     star_particle_type[ipart] == STAR_TYPE_STARII &&
		     star_age > OneStar_stellar_lifetime(ini_mass_sol, Zsol)
		    ){
		    /* delete star*/      
		    delete_particle(iter_cell,ipart);
		    particle_free(ipart); 
		    /* deleted current particle so go back one */
		    ipart = ipart_prev; 
		}

		/* go to next particle in list */
		if( ipart != NULL_PARTICLE ){  ipart = particle_list_next[ipart]; }
	    }
    }
    cart_free(level_cells);
    
}


#endif /* STAR_PARTICLE_TYPES */
#endif /* STARFORM */
