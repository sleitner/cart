#include "config.h"
#ifdef PARTICLES

#include <math.h>

#include "auxiliary.h"
#include "iterators.h"
#include "particle.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"

#include "step.h"


void move_particles( int level ) {
	int i, m;
	int ipart;
	int iter_cell;
	int num_level_cells;
	int *level_cells;
	int icell, icell_orig;
	int level1;
	int child;
	double pos[nDim];
	int found;
	double t_next;
	int c[num_children];
	double diff1, diff2, diff3;
	double pt3, pd3;
	double t1,t2,t3,d1,d2,d3;
	double t2t1, t2d1, d2t1, d2d1;
	double x, y, z, vx, vy, vz, ax, ay, az;
	double t3t2t1, t3t2d1, t3d2t1, t3d2d1;
	double d3t2t1, d3t2d1, d3d2t1, d3d2d1;
	double pconst;
	double delta_t;

	cart_assert( level >= min_level && level <= max_level );
        
	start_time( MOVE_PARTS_TIMER );
	start_time( WORK_TIMER ); 

#ifdef STARFORM
	setup_star_formation_feedback(level);
#endif /* STARFORM */

	t_next = tl[level] + dtl[level];

	select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );

#ifdef OPENMP_DECLARE_CONST

#ifdef GRAVITY
#pragma omp parallel for default(none), private(iter_cell,ipart,x,y,z,icell,level1,found,icell_orig,pos,child,i,c,diff1,diff2,diff3,d1,d2,d3,t1,t2,t3,pconst,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1,d3t2t1,d3t2d1,d3d2t1,d3d2d1,ax,ay,az,vx,vy,vz,delta_t,t2t1,t2d1,d2t1,d2d1), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,particle_x,t_next,particle_dt,particle_v,particle_id,particle_species_indices,num_particle_species,particle_list_next,particle_pot,particle_mass,cell_vars), schedule(dynamic), shared(neighbor_moves)
#else
#pragma omp parallel for default(none), private(iter_cell,ipart,x,y,z,icell,level1,found,icell_orig,pos,child,i,c,diff1,diff2,diff3,d1,d2,d3,t1,t2,t3,pconst,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1,d3t2t1,d3t2d1,d3d2t1,d3d2d1,ax,ay,az,vx,vy,vz,delta_t), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,particle_x,t_next,particle_dt,particle_v,particle_id,particle_species_indices,num_particle_species,particle_list_next), shared(neighbor_moves)
#endif

#else  /* OPENMP_DECLARE_CONST */

#ifdef GRAVITY
#pragma omp parallel for default(none), private(iter_cell,ipart,x,y,z,icell,level1,found,icell_orig,pos,child,i,c,diff1,diff2,diff3,d1,d2,d3,t1,t2,t3,pconst,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1,d3t2t1,d3t2d1,d3d2t1,d3d2d1,ax,ay,az,vx,vy,vz,delta_t,t2t1,t2d1,d2t1,d2d1), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,particle_x,t_next,particle_dt,particle_v,particle_id,particle_species_indices,num_particle_species,particle_list_next,particle_pot,particle_mass,cell_vars), schedule(dynamic)
#else
#pragma omp parallel for default(none), private(iter_cell,ipart,x,y,z,icell,level1,found,icell_orig,pos,child,i,c,diff1,diff2,diff3,d1,d2,d3,t1,t2,t3,pconst,pt3,pd3,t3t2t1,t3t2d1,t3d2t1,t3d2d1,d3t2t1,d3t2d1,d3d2t1,d3d2d1,ax,ay,az,vx,vy,vz,delta_t), shared(num_level_cells,level_cells,cell_particle_list,particle_level,level,particle_t,particle_x,t_next,particle_dt,particle_v,particle_id,particle_species_indices,num_particle_species,particle_list_next)
#endif

#endif /* OPENMP_DECLARE_CONST */
	for ( m = 0; m < num_level_cells; m++ ) {
		iter_cell = level_cells[m];

		ipart = cell_particle_list[iter_cell];
		while ( ipart != NULL_PARTICLE ) {
			cart_assert( ipart >= 0 && ipart < num_particles );
			cart_assert( particle_level[ipart] == level );
		
			if ( particle_t[ipart] < t_next - 0.5*dtl[max_level] ) {
				x = particle_x[ipart][0];
				y = particle_x[ipart][1];
				z = particle_x[ipart][2];

				icell = iter_cell;
				level1 = level;
				do {
					found = 1;
					icell_orig = icell;
					cart_assert( icell != NULL_OCT );

					cell_center_position( icell, pos );

					/* find lower leftmost cell */
					child = 0;
					for ( i = 0; i < nDim; i++ ) {
						if ( particle_x[ipart][i] >= pos[i] ) {
							child += (1<<i);
						}
					}
					cart_assert( child >= 0 && child < num_children );

					for ( i = 0; i < nDim; i++ ) {
						if ( neighbor_moves[child][i] == -1 ) {
							break;
						} else {
							icell = cell_neighbor(icell, neighbor_moves[child][i] );
							cart_assert( icell != NULL_OCT );

							if ( cell_level(icell) != level1 ) {
								icell = cell_parent_cell(icell_orig);
								cart_assert( icell != NULL_OCT );
								level1 = level1 - 1;
								found = 0;
								break;
							}
						}
					}
				
					if ( found ) {
						c[0] = icell;
						c[1] = cell_neighbor(icell,1);
						c[2] = cell_neighbor(icell,3);
						c[3] = cell_neighbor(c[1],3);
						c[4] = cell_neighbor(icell,5);
						c[5] = cell_neighbor(c[1],5);
						c[6] = cell_neighbor(c[2],5);
						c[7] = cell_neighbor(c[3],5);
	
						for ( i = 1; i < num_children; i++ ) {
							if ( cell_level(c[i]) != level1 ) {
								icell = cell_parent_cell(icell_orig);
								level1 = level1 - 1;
								cart_assert( icell != NULL_OCT );
								found = 0;
								break;
							}
						}
					}
				} while ( !found );
	
				for ( i = 0; i < num_children; i++ ) {
					cart_assert( cell_level(c[i]) == level1 );
				}

				cart_assert( c[0] != NULL_OCT );
				cell_center_position( c[0], pos );

#ifdef GRAVITY
				/* now we have the level on which this particle will move */
				diff1 = pos[0] - x;
				if ( fabs(diff1) > (double)(num_grid/2) ) {
					if ( diff1 > 0.0 ) {
						diff1 -= (double)(num_grid);
					} else {
						diff1 += (double)(num_grid);
					}
				}
				d1 = fabs(diff1) * cell_size_inverse[level1];
				cart_assert( d1 >= 0.0 && d1 <= 1.0 );

				diff2 = pos[1] - y;
				if ( fabs(diff2) > (double)(num_grid/2) ) {
					if ( diff2 > 0.0 ) {
						diff2 -= (double)(num_grid);
					} else {
						diff2 += (double)(num_grid);
					}
				}
				d2 = fabs(diff2) * cell_size_inverse[level1];

				diff3 = pos[2] - z;
				if ( fabs(diff3) > (double)(num_grid/2) ) {
					if ( diff3 > 0.0 ) {
						diff3 -= (double)(num_grid);
					} else {
						diff3 += (double)(num_grid);
					}
				}
				d3 = fabs(diff3) * cell_size_inverse[level1];

				cart_assert( d1 >= 0.0 && d1 <= 1.0 );
				cart_assert( d2 >= 0.0 && d2 <= 1.0 );
				cart_assert( d3 >= 0.0 && d3 <= 1.0 );

				t1   = 1.0 - d1;
				t2   = 1.0 - d2;
				t3   = 1.0 - d3;

				cart_assert( t1 >= 0.0 && t1 <= 1.0 );
				cart_assert( t2 >= 0.0 && t2 <= 1.0 );
				cart_assert( t3 >= 0.0 && t3 <= 1.0 );

				t2t1 = t2 * t1;
				t2d1 = t2 * d1;
				d2t1 = d2 * t1;
				d2d1 = d2 * d1;

				/* a2cs term in HART is computed in accel */
				pconst = tl[level] - particle_t[ipart] + 0.5*( dtl[level] + particle_dt[ipart] );

				pt3 = pconst * t3;
				pd3 = pconst * d3;

				t3t2t1 = pt3 * t2t1;
				t3t2d1 = pt3 * t2d1;
				t3d2t1 = pt3 * d2t1;
				t3d2d1 = pt3 * d2d1;
				d3t2t1 = pd3 * t2t1;
				d3t2d1 = pd3 * t2d1;
				d3d2t1 = pd3 * d2t1;
				d3d2d1 = pd3 * d2d1;

				particle_pot[ipart] = particle_mass[ipart] * 
						(	t3t2t1 * cell_potential(c[0]) +
							t3t2d1 * cell_potential(c[1]) +
							t3d2t1 * cell_potential(c[2]) +
							t3d2d1 * cell_potential(c[3]) +
							d3t2t1 * cell_potential(c[4]) +
							d3t2d1 * cell_potential(c[5]) +
							d3d2t1 * cell_potential(c[6]) +
							d3d2d1 * cell_potential(c[7]) ) / pconst;

				ax = 	t3t2t1 * cell_accel(c[0], 0) +
					t3t2d1 * cell_accel(c[1], 0) +
					t3d2t1 * cell_accel(c[2], 0) +
					t3d2d1 * cell_accel(c[3], 0) + 
					d3t2t1 * cell_accel(c[4], 0) + 
					d3t2d1 * cell_accel(c[5], 0) + 
					d3d2t1 * cell_accel(c[6], 0) + 
					d3d2d1 * cell_accel(c[7], 0);
	
				ay =    t3t2t1 * cell_accel(c[0], 1) +
					t3t2d1 * cell_accel(c[1], 1) +
					t3d2t1 * cell_accel(c[2], 1) +	
					t3d2d1 * cell_accel(c[3], 1) +
					d3t2t1 * cell_accel(c[4], 1) +
					d3t2d1 * cell_accel(c[5], 1) +
					d3d2t1 * cell_accel(c[6], 1) +
					d3d2d1 * cell_accel(c[7], 1);
	
				az =    t3t2t1 * cell_accel(c[0], 2) +
					t3t2d1 * cell_accel(c[1], 2) +
					t3d2t1 * cell_accel(c[2], 2) +
					t3d2d1 * cell_accel(c[3], 2) +
					d3t2t1 * cell_accel(c[4], 2) +
					d3t2d1 * cell_accel(c[5], 2) +
					d3d2t1 * cell_accel(c[6], 2) +
					d3d2d1 * cell_accel(c[7], 2);

				vx = particle_v[ipart][0] + ax;
				vy = particle_v[ipart][1] + ay;
				vz = particle_v[ipart][2] + az;
#else
				vx = particle_v[ipart][0];
				vy = particle_v[ipart][1];
				vz = particle_v[ipart][2];
#endif /* GRAVITY */
	
				delta_t = t_next - particle_t[ipart];

#ifdef STARFORM
				/* do feedback, enrichment, etc */
				if ( particle_is_star(ipart) ) {
				    stellar_feedback(level,iter_cell,ipart,delta_t,t_next,vx,vy,vz);
				}
#endif /* STARFORM */

				x += vx * delta_t;
				y += vy * delta_t;
				z += vz * delta_t;
	
				if ( x < 0.0 ) {
					x += (double)(num_grid);		
				} else if ( x >= (double)(num_grid) ) {
					x -= (double)(num_grid);
				}
	
				if ( y < 0.0 ) {
					y += (double)(num_grid);
				} else if ( y >= (double)(num_grid) ) {
					y -= (double)(num_grid);
				}
	
				if ( z < 0.0 ) {
					z += (double)(num_grid);
				} else if ( z >= (double)(num_grid) ) {
					z -= (double)(num_grid);
				}
	
				particle_x[ipart][0] = x;
				particle_x[ipart][1] = y;
				particle_x[ipart][2] = z;
	
				particle_v[ipart][0] = vx;
				particle_v[ipart][1] = vy;
				particle_v[ipart][2] = vz;
	
				particle_t[ipart] = t_next;
				particle_dt[ipart] = dtl[level];
			}
		
			ipart = particle_list_next[ipart];
		}
	}
	cart_free( level_cells );

	end_time( WORK_TIMER );
	end_time( MOVE_PARTS_TIMER );
}

#endif /* PARTICLES */
