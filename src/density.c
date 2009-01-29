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
	assign_particle_density( level );
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

	end_time( WORK_TIMER );

#ifdef RADIATIVE_TRANSFER
	rtAfterAssignDensity2(level,num_level_cells,level_cells);
#endif

	start_time( WORK_TIMER );
	cart_free( level_cells );
	end_time( WORK_TIMER );
}
#endif /* HYDRO */

#ifdef PARTICLES

void assign_particle_density( int level ) {
	int i, j;
	int ipart;
	int icell;
	double corner[nDim];
	double size2, size_inverse;
	float mass;
	double cornerx0, cornerx1, cornery0, cornery1, cornerz0, cornerz1;
	double x, y, z;
	double xs, ys, zs;
	double dx0, dx1, dy0, dy1, dz0, dz1;
	double d00, d01, d10, d11;
	int is_first;
#ifdef RT_VAR_SOURCE
	float sor;
#endif

	start_time( WORK_TIMER );

	size2 = 0.5*cell_size[level];
	size_inverse = cell_size_inverse[level];

	for ( ipart = 0; ipart < num_particles; ipart++ ) {
		if ( particle_level[ipart] >= level ) {
#ifdef STARFORM
			is_first = ( particle_id[ipart] < particle_species_indices[1] || particle_is_star(ipart) );
#else
			is_first = ( particle_id[ipart] < particle_species_indices[1] );
#endif /* STARFORM */

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

#ifdef RT_VAR_SOURCE
			sor = rtSource(ipart);
#endif
			/* child 0 */
			corner[0] = cornerx0;
			corner[1] = cornery0;
			corner[2] = cornerz0;

			icell = cell_find_position_level( level, corner );
			if ( icell != -1 ) {
				mass = d00*dz0;
				cell_density(icell) += mass;
#ifdef RT_VAR_SOURCE
				cell_rt_source(icell) += mass*sor;
#endif	
				if ( is_first ) {
					cell_first_species_mass(icell) += mass;
				}
			}

			/* child 1 */
			corner[0] = cornerx1;

			icell = cell_find_position_level( level, corner );
			if ( icell != -1 ) {
				mass = d10*dz0;
				cell_density(icell) += mass;
#ifdef RT_VAR_SOURCE
				cell_rt_source(icell) += mass*sor;
#endif	
				if ( is_first ) {
					cell_first_species_mass(icell) += mass;
				}
			}

			/* child 2 */
			corner[0] = cornerx0;
			corner[1] = cornery1;

			icell = cell_find_position_level( level, corner );
			if ( icell != -1 ) {
				mass = d01*dz0;
				cell_density(icell) += mass;
#ifdef RT_VAR_SOURCE
				cell_rt_source(icell) += mass*sor;
#endif	
				if ( is_first ) {
					cell_first_species_mass(icell) += mass;
				}
			}

			/* child 3 */
			corner[0] = cornerx1;

			icell = cell_find_position_level( level, corner );
			if ( icell != -1 ) {
				mass = d11*dz0;
				cell_density(icell) += mass;
#ifdef RT_VAR_SOURCE
				cell_rt_source(icell) += mass*sor;
#endif	
				if ( is_first ) {
					cell_first_species_mass(icell) += mass;
				}
			}

			/* child 4 */
			corner[0] = cornerx0;
			corner[1] = cornery0;
			corner[2] = cornerz1;

			icell = cell_find_position_level( level, corner );
			if ( icell != -1 ) {
				mass = d00*dz1;
				cell_density(icell) += mass;
#ifdef RT_VAR_SOURCE
				cell_rt_source(icell) += mass*sor;
#endif	
				if ( is_first ) {
					cell_first_species_mass(icell) += mass;
				}
			}

			/* child 5 */
			corner[0] = cornerx1;

			icell = cell_find_position_level( level, corner );
			if ( icell != -1 ) {
				mass = d10*dz1;
				cell_density(icell) += mass;
#ifdef RT_VAR_SOURCE
				cell_rt_source(icell) += mass*sor;
#endif	
                                if ( is_first ) {
                                        cell_first_species_mass(icell) += mass;
                                }
			}

                        /* child 6 */
			corner[0] = cornerx0;
			corner[1] = cornery1;

			icell = cell_find_position_level( level, corner );
			if ( icell != -1 ) {
				mass = d01*dz1;
				cell_density(icell) += mass;
#ifdef RT_VAR_SOURCE
				cell_rt_source(icell) += mass*sor;
#endif	
				if ( is_first ) {
					cell_first_species_mass(icell) += mass;
				}
			}

                        /* child 7 */
			corner[0] = cornerx1;

			icell = cell_find_position_level( level, corner );
			if ( icell != -1 ) {
				mass = d11*dz1;
				cell_density(icell) += mass;
#ifdef RT_VAR_SOURCE
				cell_rt_source(icell) += mass*sor;
#endif	
				if ( is_first ) {
					cell_first_species_mass(icell) += mass;
				}
			}
		}
	}

	end_time( WORK_TIMER );

	/* now collect buffer densities (after this function we
	 *  have correct densities on all cells/buffer cells) */
	merge_buffer_cell_densities( level );
}

#endif /* PARTICLES */

#endif /* GRAVITY */
