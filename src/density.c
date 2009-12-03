#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "timestep.h"
#include "density.h"
#include "iterators.h"
#include "tree.h"
#include "cell_buffer.h"
#include "particle.h"
#include "constants.h"
#include "auxiliary.h"
#include "timing.h"

#ifndef DENSITY_CHUNK_SIZE
#define DENSITY_CHUNK_SIZE	16384
#endif /* DENSITY_CHUNK_SIZE */

#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#endif


#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER)

void initialize_density( int level ) {
	int i;
	int icell;
        int num_level_cells;
        int *level_cells;

	start_time( WORK_TIMER );

#ifdef RT_VAR_SOURCE
	rtInitSource(level);
#endif

#ifdef PARTICLES
	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,cell_volume,cell_vars,level)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		cell_density(icell) = -cell_volume[level]; 
		cell_first_species_mass(icell) = 0.0;
#ifdef RT_VAR_SOURCE
		cell_rt_source(icell) = 0.0;
#endif	
	}

	cart_free( level_cells );

	select_level( level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,cell_vars)
        for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		cell_density(icell) = 0.0; 
		cell_first_species_mass(icell) = 0.0;
#ifdef RT_VAR_SOURCE
		cell_rt_source(icell) = 0.0;
#endif	
	}
	
	cart_free( level_cells );
#else
	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,cell_vars,level,cell_volume)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
#ifdef GRAVITY
		cell_density(icell) = -cell_volume[level];
#endif
#ifdef RT_VAR_SOURCE
		cell_rt_source(icell) = 0.0;
#endif	
	}

	cart_free( level_cells );
#endif /* PARTICLES */

	end_time( WORK_TIMER );
}

void assign_density( int level ) {

	start_time( DENSITY_TIMER );

	initialize_density(level);

#ifdef PARTICLES

#ifdef MAX_LEVEL_DARK_DENSITY
        if ( level <= MAX_LEVEL_DARK_DENSITY ) {
          assign_particle_density( level );
        } else {
          assign_particle_density_smoothed( level );
        }
#else
	assign_particle_density( level );
#endif

#endif /* PARTICLES */

#ifdef HYDRO
	assign_hydro_density( level );
#else /* HYDRO */
#ifdef RADIATIVE_TRANSFER
	rtAfterAssignDensity1(level);
#endif
#endif /* HYDRO */

	end_time( DENSITY_TIMER );
}


#ifdef HYDRO
void assign_hydro_density( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;

	start_time( WORK_TIMER );

	/* assumes buffer gas density is up to date */
	select_level( level, CELL_TYPE_ANY, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,cell_vars,level,cell_volume)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
#ifdef GRAVITY
		cell_density(icell) += cell_gas_density(icell)*cell_volume[level];
#endif
	}

#ifdef RADIATIVE_TRANSFER
	rtAfterAssignDensity2(level,num_level_cells,level_cells);
#endif

	cart_free( level_cells );

	end_time( WORK_TIMER );
}
#endif /* HYDRO */

#ifdef PARTICLES
void find_single_particle_density( int level, double size2, double size_inverse, int ipart, 
		int *cell_list, float *mass_assigned ) {
        double corner[nDim];
        double cornerx0, cornerx1, cornery0, cornery1, cornerz0, cornerz1;
        double x, y, z;
        double xs, ys, zs;
        double dx0, dx1, dy0, dy1, dz0, dz1;
        double d00, d01, d10, d11;

	x = particle_x[ipart][0];
	y = particle_x[ipart][1];
	z = particle_x[ipart][2];

	cornerx0 = x - size2;
	cornerx1 = x + size2;
	cornery0 = y - size2;
	cornery1 = y + size2;
	cornerz0 = z - size2;
	cornerz1 = z + size2;

	if ( cornerx0 < 0.0 ) cornerx0 += (double)num_grid;
	if ( cornerx1 >= (double)num_grid ) cornerx1 -= (double)num_grid;
	if ( cornery0 < 0.0 ) cornery0 += (double)num_grid;
	if ( cornery1 >= (double)num_grid ) cornery1 -= (double)num_grid;
	if ( cornerz0 < 0.0 ) cornerz0 += (double)num_grid;
	if ( cornerz1 >= (double)num_grid ) cornerz1 -= (double)num_grid;

	xs = x*size_inverse + 0.5;
	ys = y*size_inverse + 0.5;
	zs = z*size_inverse + 0.5;

	dx1 = xs - floor(xs);
	dy1 = ys - floor(ys);
	dz1 = zs - floor(zs);

	dx0 = 1.0 - dx1;
	dy0 = 1.0 - dy1;
	dz0 = 1.0 - dz1;

	dx0 *= particle_mass[ipart];
	dx1 *= particle_mass[ipart];

	d00 = dx0*dy0;
	d01 = dx0*dy1;
	d10 = dx1*dy0;
	d11 = dx1*dy1;

	/* child 0 */
	corner[0] = cornerx0;
	corner[1] = cornery0;
	corner[2] = cornerz0;

	cell_list[0] = cell_find_position_level( level, corner );
	mass_assigned[0] = d00*dz0;

	/* child 1 */
	corner[0] = cornerx1;

	cell_list[1] = cell_find_position_level( level, corner );
	mass_assigned[1] = d10*dz0;

	/* child 2 */
	corner[0] = cornerx0;
	corner[1] = cornery1;

	cell_list[2] = cell_find_position_level( level, corner );
	mass_assigned[2] = d01*dz0;

	/* child 3 */
	corner[0] = cornerx1;

	cell_list[3] = cell_find_position_level( level, corner );
	mass_assigned[3] = d11*dz0;

	/* child 4 */
	corner[0] = cornerx0;
	corner[1] = cornery0;
	corner[2] = cornerz1;

	cell_list[4] = cell_find_position_level( level, corner );
	mass_assigned[4] = d00*dz1;

	/* child 5 */
	corner[0] = cornerx1;

	cell_list[5] = cell_find_position_level( level, corner );
	mass_assigned[5] = d10*dz1;

	/* child 6 */
	corner[0] = cornerx0;
	corner[1] = cornery1;

	cell_list[6] = cell_find_position_level( level, corner );
	mass_assigned[6] = d01*dz1;

	/* child 7 */
	corner[0] = cornerx1;

	cell_list[7] = cell_find_position_level( level, corner );
	mass_assigned[7] = d11*dz1;
}


void assign_particle_density( int level ) {
        int i, j, k;
	int count;
        int ipart;
        int icell;
        int is_first;
	double size2, size_inverse;
#ifdef RT_VAR_SOURCE
        float sor;
#endif
	int particle_list[DENSITY_CHUNK_SIZE];
	int cell_list[num_children*DENSITY_CHUNK_SIZE];
	float mass_assigned[num_children*DENSITY_CHUNK_SIZE];

        start_time( WORK_TIMER );

        size2 = 0.5*cell_size[level];
        size_inverse = cell_size_inverse[level];
	count = 0;

        for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] >= level ) {
			particle_list[count++] = i;
		}

		if ( i == num_particles-1 || count == DENSITY_CHUNK_SIZE ) {
#pragma omp parallel for default(none), private(j,ipart), shared(count,level,size2,size_inverse,particle_list,cell_list,mass_assigned)
			for ( j = 0; j < count; j++ ) {
				find_single_particle_density( level, size2, size_inverse, 
					particle_list[j], &cell_list[num_children*j],
					&mass_assigned[num_children*j] );
			}

			/* must execute serially! */
			for ( j =  0; j < count; j++ ) {
				ipart = particle_list[j];

#ifdef STARFORM
				is_first = ( particle_id[ipart] < particle_species_indices[1] || particle_is_star(ipart) );
#else
				is_first = ( particle_id[ipart] < particle_species_indices[1] );
#endif /* STARFORM */

#ifdef RT_VAR_SOURCE
				sor = rtSource(ipart);
#endif

				if ( is_first ) {
					for ( k = 0; k < num_children; k++ ) {
						icell = cell_list[num_children*j+k];

						if ( icell != -1 ) {
							cell_density(icell) += mass_assigned[num_children*j+k];
							cell_first_species_mass(icell) += mass_assigned[num_children*j+k];
#ifdef RT_VAR_SOURCE
							cell_rt_source(icell) += sor*mass_assigned[num_children*j+k];
#endif
						}
					}
				} else {
					for ( k = 0; k < num_children; k++ ) {
						icell = cell_list[num_children*j+k];

						if ( icell != -1 ) {
							cell_density(icell) += mass_assigned[num_children*j+k];
#ifdef RT_VAR_SOURCE
							cell_rt_source(icell) += sor*mass_assigned[num_children*j+k];
#endif
						}
					}
				}
			}

			count = 0;
		}
	}

	end_time(WORK_TIMER);

	/* now collect buffer densities (after this function we
	 *  have correct densities on all cells/buffer cells) */
	merge_buffer_cell_densities( level );
}

#ifdef MAX_LEVEL_DARK_DENSITY
/* snl With a few more ifdefs all this could be done in the original routine */
void assign_particle_density_smoothed( int level ) { 
        int i, j, k;
	int count;
        int ipart;
        int icell, icell0;
        int is_first;
        int is_star;
	double size2, size_inverse;
	double size2_star, size_star_inverse;
#ifdef RT_VAR_SOURCE
        float sor;
#endif
	int particle_list[DENSITY_CHUNK_SIZE];
	int cell_list[num_children*DENSITY_CHUNK_SIZE];
	float mass_assigned[num_children*DENSITY_CHUNK_SIZE];

	int index_to_child; 
	int high_levels;
	int children_at_level;
	int cells_below_ilev;
	int ic,icid;
	int ilev;
	
        start_time( WORK_TIMER );

	high_levels = level - MAX_LEVEL_DARK_DENSITY  ;
	cart_assert ( high_levels > 0 );
	children_at_level = pow( 2, 3*high_levels );

	/* smoothed interpolation is done on MAX_LEVEL_DARK_DENSITY so sizes reflect this*/
        size2 = 0.5*cell_size[MAX_LEVEL_DARK_DENSITY];
        size_inverse = cell_size_inverse[MAX_LEVEL_DARK_DENSITY];

        size2_star = 0.5*cell_size[level];
        size_star_inverse = cell_size_inverse[level];

	count = 0;
        for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] >= level ) {
		  particle_list[count++] = i; 
		}

		if ( i == num_particles-1 || count == DENSITY_CHUNK_SIZE ) {
#pragma omp parallel for default(none), private(j,ipart,is_star), shared(count,level,size2,size_inverse,size2_star,size_star_inverse,particle_list,cell_list,mass_assigned,particle_species_indices,num_particle_species,particle_id)
			for ( j = 0; j < count; j++ ) {
				ipart = particle_list[j];
#ifdef STARFORM				
				is_star = ( particle_is_star(ipart) );
#else 
				is_star = 0;
#endif
				if ( is_star ) {
				  /* stars put density on their own leaves */
				  find_single_particle_density( level, size2_star, size_star_inverse, 
								ipart, &cell_list[num_children*j],
								&mass_assigned[num_children*j] );
				} else {
				  /* assign density of particles in list at MAX_LEVEL_DARK_DENSITY*/
				  find_single_particle_density( MAX_LEVEL_DARK_DENSITY, size2, size_inverse, 
								ipart, &cell_list[num_children*j],
								&mass_assigned[num_children*j] );
				}
			}

			/* must execute serially! */
			for ( j =  0; j < count; j++ ) {
				ipart = particle_list[j];

				is_first = particle_is_first(ipart);
#ifdef STARFORM				
				is_star =  particle_is_star(ipart);
#else 
				is_star = 0;
#endif

#ifdef RT_VAR_SOURCE
				sor = rtSource(ipart);
#endif
	
				if ( is_star ) {
      				        for ( k = 0; k < num_children; k++ ) {
					        icell = cell_list[num_children*j+k];

						if ( icell != -1 ) {
						      cell_density(icell) += mass_assigned[num_children*j+k] ;
						      cell_first_species_mass(icell) += mass_assigned[num_children*j+k];
#ifdef STARFORM
						      cell_stellar_particle_mass(icell) += mass_assigned[num_children*j+k];
#endif
#ifdef RT_VAR_SOURCE
						      cell_rt_source(icell) += sor*mass_assigned[num_children*j+k];
#endif
						}
					}
				} else {
				  //cart_assert( is_first || particle_mass[ipart]==0 ); //highspecies or stars only in these high levels, otherwise contaminated.
					for ( k = 0; k < num_children; k++ ) {
						icell0 = cell_list[num_children*j+k];
						
						if ( icell0 != -1 ) {
						  mass_assigned[num_children*j+k] /= pow( 2., 3*(high_levels));
						  
						  /* now find all youngest relatives of these cells to level "level" and give the fraction of mass_assigned to
						   * cell_density[icell] to cells that are refined at or below the particle level =level */
						  ic = 0;
						  while ( ic < children_at_level ) {
						    icid = ic;
						    icell = icell0;
						    ilev = 0;

						    while ( ilev < high_levels && icell !=-1 ){
						      cells_below_ilev = pow(  2 , 3*( high_levels-(ilev+1) )  );
                                                      index_to_child = (int)icid / cells_below_ilev ;
                                                      icid -= cells_below_ilev * index_to_child;
						      icell = cell_child ( icell, index_to_child );
						      ilev++;
						    }

						    if ( icell != -1 ) {
						      cell_density(icell) += mass_assigned[num_children*j+k] ;
						      if(is_first){
							cell_first_species_mass(icell) += mass_assigned[num_children*j+k];
						      }
#ifdef RT_VAR_SOURCE
						      cell_rt_source(icell) += sor*mass_assigned[num_children*j+k];
#endif
						      ic++; 
						    }else {
						      ic += pow( 2,3*(high_levels-ilev) ) ;
						    }
						    
						  }
						  
						}
					}
					
				}
			}

			count = 0;
		}
	}

	end_time(WORK_TIMER);

	/* now collect buffer densities (after this function we
	 *  have correct densities on all cells/buffer cells) */
	merge_buffer_cell_densities( level );
}
#endif /* MAX_LEVEL_DARK_DENSITY */

#endif /* PARTICLES */

#endif /* GRAVITY */
