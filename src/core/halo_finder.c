#include "config.h"
#ifdef PARTICLES

#include <dirent.h>
#include <stdio.h>
#include <math.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "halos.h"
#include "halo_finder.h"
#include "io.h"
#include "parallel.h"
#include "particle.h"
#include "particle_buffer.h"
#include "particle_density.h"
#include "stack.h"
#include "timing.h"
#include "times.h"
#include "tree.h"
#include "units.h"

/* define binning */
int num_bins = 80;
double rmin_physical = 1e-3;
double rmax_physical = 2.0;
double delta_vir = 180.;
double delta_halo_center = 1000;
double min_halo_center_level = max_level-2;
double min_halo_mass = 0.0; 

int halo_finder_frequency = 0;              /* repeated halo finding disabled by default */
int halo_particle_list_flag = 1;			/* write hp_*.dat files mapping particles to halos and computing binding energy */
int halo_center_definition = 0;             /* 0 = cm convergence, 1 density peak */
int halo_num_nearest_neighbors = 24;
double cm_radius_initial_reduction = 0.333;
double cm_radius_freduce = 0.05;
double cm_convergence_ftol = 1e-4;
double cm_convergence_abs = 5.0;            /* in units of cell_size[max_level] */

double rlmin, rlmax, drl;
double rl[MAX_HALO_BINS];
double rr[MAX_HALO_BINS];
double bin_volume[MAX_HALO_BINS];
double bin_volume_cumulative[MAX_HALO_BINS];

extern int step;

char halo_finder_output_directory_d[CONTROL_PARAMETER_STRING_LENGTH] = ".";
const char* halo_finder_output_directory = halo_finder_output_directory_d;

void config_init_halo_finder() {
	control_parameter_add2(control_parameter_string,halo_finder_output_directory_d,"directory:halo-finder","halo_finder_output_directory","directory for output halo catalog files." );
	control_parameter_add(control_parameter_int,&halo_finder_frequency,"frequency:halo-finder","Sets number of root level steps between halo finding steps");

	control_parameter_add(control_parameter_int,&halo_center_definition,"halo:center_definition","Halo centers are defined as 0=iterative cm convergence or 1=location of highest density particle");
	control_parameter_add(control_parameter_int,&halo_num_nearest_neighbors,"halo:num_nearest_neighbors","Number of nearest neighbors used to define smoothed density for initial peak finding");
	control_parameter_add(control_parameter_double,&cm_radius_initial_reduction,"halo:cm_radius_initial_reduction","Initial fraction of rvir to use for iterative center of mass calculation");
	control_parameter_add(control_parameter_double,&cm_radius_freduce,"halo:cm_radius_freduce","Fraction of center of mass search radius to decrease during each iteration.");
	control_parameter_add(control_parameter_double,&cm_convergence_ftol,"halo:cm_convergence_ftol","Fractional movement in halo center of mass necessary for convergence.");
	control_parameter_add(control_parameter_double,&cm_convergence_abs,"halo:cm_convergence_abs","Absolute movement in halo center of mass (units of cell_size[max_level]) necessary for convergence.");
    control_parameter_add(control_parameter_double,&rmin_physical,"halo:rhalo_bin_rmin", "Minimum radius to search for overdensity radius rvir [comoving Mpc/h]");
	control_parameter_add(control_parameter_double,&rmax_physical,"halo:rhalo_bin_rmax", "Maximum radius to serach for overdensity radius rvir [comoving Mpc/h]");
	control_parameter_add2(control_parameter_double,&delta_vir,"halo:overdensity","halo:delta_vir", "Overdensity with respect to critical which defines halo rvir and Mvir");
	control_parameter_add2(control_parameter_double,&delta_halo_center,"halo:center_overdensity","halo:deltamin", "Overdensity necessary for a particle to be considered a potential halo center");
	control_parameter_add(control_parameter_int,&min_halo_center_level,"halo:center_min_level","Minimum level to compute smoothed densities and therefore be eligible to be a halo center");
	control_parameter_add(control_parameter_double,&min_halo_mass,"halo:min_halo_mass","Minimum mass within chosen overdensity to consider a halo");
}

void config_verify_halo_finder() {
	DIR *d;

	if((d = opendir(halo_finder_output_directory)) == NULL) {
		cart_error("Directory %s does not exist.",halo_finder_output_directory);
	} else closedir(d);

	VERIFY( frequency:halo-finder, halo_finder_frequency >= 0 );
	VERIFY( halo:center_definition, halo_center_definition == 0 || halo_center_definition == 1 );
	VERIFY( halo:num_nearest_neighbors, halo_num_nearest_neighbors > 0 );
	VERIFY( halo:cm_radius_initial_reduction, cm_radius_initial_reduction > 0.0 && cm_radius_initial_reduction <= 1.0 );
	VERIFY( halo:cm_radius_freduce, cm_radius_freduce > 0.0 && cm_radius_freduce < 1.0 );
	VERIFY( halo:cm_convergence_ftol, cm_convergence_ftol > 0.0 && cm_convergence_ftol < 1.0 );
	VERIFY( halo:cm_convergence_abs, cm_convergence_abs >= 0.0 );
	VERIFY( halo:rhalo_bin_rmin, rmin_physical > 0 && rmin_physical < rmax_physical );
	VERIFY( halo:rhalo_bin_rmax, rmax_physical > 0.0 );
	VERIFY( halo:overdensity, delta_vir > 0.0 );
	VERIFY( halo:center_overdensity, delta_halo_center >= delta_vir );
	VERIFY( halo:center_min_level, min_halo_center_level >= min_level && min_halo_center_level <= max_level );
	VERIFY( halo:min_halo_mass, min_halo_mass >= 0.0 );
}

double loglin_interpolate( double *binned_var, int bin, double rlout, double rri, double rll ) {
	double ah, bh;

	ah = ( binned_var[bin] - binned_var[bin-1] ) * rri;
	bh = binned_var[bin] - ah * rll;
	return ah*rlout + bh;
}

double log_interpolate( double *binned_var, int bin, double rlout, double rri, double rll ) {
	double aM1, aM2, ah, bh;

	aM1 = ( binned_var[bin-1] > 0.0 ) ? log10( binned_var[bin-1] ) : -15.0;
	aM2 = ( binned_var[bin] > 0.0 ) ? log10( binned_var[bin] ) : -15.0;
	ah = ( aM2 - aM1 ) * rri;
	bh = aM1 - ah * rll;
	return pow( 10.0, ah*rlout + bh );
}

int compute_radial_bin( double r ) {
	double lgr = log10(r);
	return ( lgr < rlmin ) ? 0 : (int)((lgr - rlmin)/drl) + 1;
}

void compute_halo_mass( halo *h ) {
	int i, j, k, m;
	int icell;
	int ioct;
	int ipart;
	int bin;
	double r;
	int coords[nDim];
	stack *cell_list;
	double dbi1, dbi2;
	double dlbi1, dlbi2;
	double rrl, rll, rri;
	double bin_mass[MAX_HALO_BINS];
	double thread_local_mass[MAX_HALO_BINS];
	double local_mass[MAX_HALO_BINS];
	double bin_total_mass[MAX_HALO_BINS];
	double thread_local_cv[nDim][MAX_HALO_BINS];
	double local_cv[nDim][MAX_HALO_BINS];
	double bin_cv[nDim][MAX_HALO_BINS];
	double bin_total_cv[nDim][MAX_HALO_BINS];
	double total_cv[nDim];
	double total_mass;
	double Gcode;
	double vcirc, vmax, rmax;
	double dx, lgr;
	int level, child, parent;

	start_time( HALO_FINDER_MASS_TIMER );

	for ( bin = 0; bin < num_bins; bin++ ) {
		local_mass[bin] = 0.0;
		for ( i = 0; i < nDim; i++ ) {
			local_cv[i][bin] = 0.0;
		}
	}

#pragma omp parallel default(none) shared(h,rr,cell_child_oct,cell_particle_list,particle_list_next,cell_vars,particle_mass,particle_x,particle_v,num_bins,local_mass,local_cv,oct_pos,oct_level,cell_delta,rlmin,drl) private(bin,i,j,k,m,coords,cell_list,icell,ipart,r,thread_local_mass,thread_local_cv,lgr,dx,parent,child,level,ioct)
	{
		cell_list = stack_init();

		for ( bin = 0; bin < num_bins; bin++ ) {
			thread_local_mass[bin] = 0.0;
			for ( i = 0; i < nDim; i++ ) {
				thread_local_cv[i][bin] = 0.0;
			}
		}

		/* select root cells, assumes rmax ~ cell_size[min_level] so no pruning */
#pragma omp for collapse(3) schedule(dynamic) nowait
		for ( i = (int)floor(h->pos[0]-rr[num_bins-1]); i <= (int)(h->pos[0]+rr[num_bins-1]); i++ ) {
			for ( j = (int)floor(h->pos[1]-rr[num_bins-1]); j <= (int)(h->pos[1]+rr[num_bins-1]); j++ ) {
				for ( k = (int)floor(h->pos[2]-rr[num_bins-1]); k <= (int)(h->pos[2]+rr[num_bins-1]); k++ ) {
					coords[0] = ( i + num_grid ) % num_grid;
					coords[1] = ( j + num_grid ) % num_grid;
					coords[2] = ( k + num_grid ) % num_grid;
					icell = root_cell_location( sfc_index( coords ) );
					if ( icell != NULL_OCT && cell_is_local(icell) ) {
						stack_push( cell_list, icell );

						while ( stack_pop( cell_list, &icell ) ) {
							/* inlined cell_position_double+compute_distance_periodic+compute_radial_bin */
							if ( cell_is_root_cell(icell) ) {
								level = min_level;
								r = 0.0;
								for ( m = 0; m < nDim; m++ ) {
									dx = (double)coords[m] + 0.5 - h->pos[m];
									if ( dx < -num_grid/2 ) {
										dx += num_grid;
									} else if ( dx > num_grid/2 ) {
										dx -= num_grid;
									}
									r += dx*dx;
								}
							} else {
								parent = cell_parent_oct(icell);
								level = oct_level[parent];
								child = cell_child_number(icell);

								r = 0.0;
								for ( m = 0; m < nDim; m++ ) {
									dx = (double)oct_pos[parent][m] + 
										(double)cell_size[level]*(double)cell_delta[child][m] - 
										h->pos[m];
									if ( dx < -num_grid/2 ) {
										dx += num_grid;
									} else if ( dx > num_grid/2 ) {
										dx -= num_grid;
									}
									r += dx*dx;
								}
							}
							lgr = 0.5*log10(r);
							bin = ( lgr < rlmin ) ? 0 : (int)((lgr - rlmin)/drl) + 1;

							if ( cell_is_refined(icell) ) {
								ioct = cell_child_oct[icell];
								for ( m = 0; m < num_children; m++ ) {
									stack_push( cell_list, oct_child( ioct, m ) );
								}
#ifdef HYDRO
							} else {
								if ( bin < num_bins ) { 
									thread_local_mass[bin] += cell_gas_density(icell)*cell_volume[level];
									for ( m = 0; m < nDim; m++ ) {
										thread_local_cv[m][bin] += cell_momentum(icell,m)*cell_volume[level];
									}
								}
#endif /* HYDRO */
							}

							if ( bin < num_bins ) {
								ipart = cell_particle_list[icell];
								while ( ipart != NULL_PARTICLE ) {
									thread_local_mass[bin] += particle_mass[ipart];
									for ( m = 0; m < nDim; m++ ) {
										thread_local_cv[m][bin] += particle_mass[ipart]*particle_v[ipart][m];
									}
									ipart = particle_list_next[ipart];
								}
							}
						}
					}
				}
			}
		}

#pragma omp critical
		{
			for ( bin = 0; bin < num_bins; bin++ ) {
				local_mass[bin] += thread_local_mass[bin];
				for ( i = 0; i < nDim; i++ ) {
					local_cv[i][bin] += thread_local_cv[i][bin];
				}
			}
		}

		stack_destroy(cell_list); 
	}

	MPI_Allreduce( local_mass, bin_mass, num_bins, MPI_DOUBLE, MPI_SUM, mpi.comm.run );
	for ( i = 0; i < nDim; i++ ) {
		MPI_Allreduce( local_cv[i], bin_cv[i], num_bins, MPI_DOUBLE, MPI_SUM, mpi.comm.run );
	}

	/* find rvir, mvir, rmax, vmax */
	Gcode = constants->G/units->length/units->length/units->length*units->mass*units->time*units->time;
	vmax = 0.0;
	total_mass = 0.0;
	for ( i = 0; i < nDim; i++ ) {
		total_cv[i] = 0.0;
	}

	for ( bin = 0; bin < num_bins; bin++ ) {
		total_mass += bin_mass[bin];
		bin_total_mass[bin] = total_mass;

		for ( i = 0; i < nDim; i++ ) {
			total_cv[i] += bin_cv[i][bin];

			if ( bin_total_mass[bin] > 0.0 ) {
				bin_total_cv[i][bin] = total_cv[i]/bin_total_mass[bin];
			} else {
				bin_total_cv[i][bin] = 0.0;
			}
		}

		if ( bin == 0 ) {
			dbi1 = 0.0;
		} else {
			dbi1 = ( bin_total_mass[bin-1] / bin_volume_cumulative[bin-1] );
		}

		dbi2 = bin_total_mass[bin] / bin_volume_cumulative[bin];

		vcirc = sqrt( Gcode*bin_total_mass[bin] / rr[bin] );
		if ( vcirc >= vmax ) {
			rmax = rr[bin];
			vmax = vcirc;
		}

		/* compute virial radius */
		if ( bin > 0 && ( dbi1 >= delta_vir && dbi2 < delta_vir ) ) {
			rrl = log10(rr[bin]);
			rll = log10(rl[bin]);
			rri = 1.0/(rrl-rll);
			dlbi1 = log10(dbi1);
			dlbi2 = log10(dbi2);
			h->rvir = pow( 10.0, (log10(delta_vir)*(rrl-rll) + 
							rll*dlbi2 - rrl*dlbi1)/(dlbi2-dlbi1));
			h->mvir = log_interpolate( bin_total_mass, bin, log10(h->rvir), rri, rll );

			for ( i = 0; i < nDim; i++ ) {
				h->vel[i] = loglin_interpolate( bin_total_cv[i], bin, log10(h->rvir), rri, rll );
			}

			break;
		}
	}

	if ( bin == num_bins ) {
		h->rvir = rr[num_bins-1];
		h->mvir = bin_total_mass[num_bins-1];

		for ( i = 0; i < nDim; i++ ) {
			h->vel[i] = bin_total_cv[i][num_bins-1];
		}
	}

	h->rmax = rmax;
	h->vmax = vmax;

	end_time( HALO_FINDER_MASS_TIMER );
}

void halo_recenter( halo *h ) {
	int i, j, k, m;
	int icell, ipart;
	double r, rcm, dr;
	int niter;
	double dx[nDim];
	double cm[nDim], cm_total[nDim], local_cm[nDim];
	double cm_mass, cm_mass_total, local_cm_mass;
	float mass;
	int coords[nDim];
	stack *cell_list;
	int ioct;
	int parent, child, level;

	start_time( HALO_FINDER_RECENTER_TIMER );

	rcm = h->rvir * cm_radius_initial_reduction;
	niter = 0;

	do {
		cm_mass = 0.0;
		for ( i = 0; i < nDim; i++ ) {
			cm[i] = 0.0;
		}

#pragma omp parallel default(none) shared(h,cm,cm_mass,rcm,cell_child_oct,cell_particle_list,particle_list_next,cell_vars,particle_mass,particle_x,num_bins,cell_volume,oct_pos,oct_level,cell_delta) private(cell_list,i,j,k,m,local_cm,local_cm_mass,coords,r,icell,ipart,dx,mass,level,ioct,parent,child)
		{
			cell_list = stack_init();
			local_cm_mass = 0.0;
			for ( i = 0; i < nDim; i++ ) {
				local_cm[i] = 0.0;
			}

			/* select root cells, assumes rmax ~ cell_size[min_level] so no pruning */
#pragma omp for collapse(3) schedule(dynamic) nowait
			for ( i = (int)floor(h->pos[0]-rcm); i <= (int)(h->pos[0]+rcm); i++ ) {
				for ( j = (int)floor(h->pos[1]-rcm); j <= (int)(h->pos[1]+rcm); j++ ) {
					for ( k = (int)floor(h->pos[2]-rcm); k <= (int)(h->pos[2]+rcm); k++ ) {
						coords[0] = ( i + num_grid ) % num_grid;
						coords[1] = ( j + num_grid ) % num_grid;
						coords[2] = ( k + num_grid ) % num_grid;
						icell = root_cell_location( sfc_index( coords ) );
						if ( icell != NULL_OCT && cell_is_local(icell) ) {
							stack_push( cell_list, icell );

							while ( stack_pop( cell_list, &icell ) ) {
								if ( cell_is_root_cell(icell) ) {
									level = min_level;
									r = 0.0;
									for ( m = 0; m < nDim; m++ ) {
										dx[m] = (double)coords[m] + 0.5 - h->pos[m];
										if ( dx[m] < -num_grid/2 ) {
											dx[m] += num_grid;
										} else if ( dx[m] > num_grid/2 ) {
											dx[m] -= num_grid;
										}
										r += dx[m]*dx[m];
									}
								} else {
									parent = cell_parent_oct(icell);
									level = oct_level[parent];
									child = cell_child_number(icell);

									r = 0.0;
									for ( m = 0; m < nDim; m++ ) {
										dx[m] = (double)oct_pos[parent][m] + 
											(double)cell_size[level]*(double)cell_delta[child][m] - 
											h->pos[m];
										if ( dx[m] < -num_grid/2 ) {
											dx[m] += num_grid;
										} else if ( dx[m] > num_grid/2 ) {
											dx[m] -= num_grid;
										}
										r += dx[m]*dx[m];
									}
								}

								if ( cell_is_refined(icell) ) { 
									ioct = cell_child_oct[icell];
									for ( m = 0; m < num_children; m++ ) {
										stack_push( cell_list, oct_child( ioct, m ) );
									}
#ifdef HYDRO
								} else {
									if ( r < rcm*rcm ) {
										mass = cell_gas_density(icell)*cell_volume[level];
										local_cm_mass += mass;
										for ( m = 0; m < nDim; m++ ) {
											local_cm[m] += mass*dx[m];
										}
									}
#endif /* HYDRO */
								}

								ipart = cell_particle_list[icell];
								while ( ipart != NULL_PARTICLE ) {
									local_cm_mass += particle_mass[ipart];
									r = 0.0;
									for ( m = 0; m < nDim; m++ ) {
										dx[m] = particle_x[ipart][m] - h->pos[m];
										if ( dx[m] < -num_grid/2 ) {
											dx[m] += num_grid;
										} else if ( dx[m] > num_grid/2 ) {
											dx[m] -= num_grid;
										}
										r += dx[m]*dx[m];
									}

									if ( r < rcm*rcm ) {
										local_cm_mass += particle_mass[ipart];
										for ( m = 0; m < nDim; m++ ) {
											local_cm[m] += particle_mass[ipart]*dx[m];
										}
									}

									ipart = particle_list_next[ipart];
								}
							}
						}
					}
				}
			}

#pragma omp critical
			{
				cm_mass += local_cm_mass;
				for ( i = 0; i < nDim; i++ ) {
					cm[i] += local_cm[i];
				}
			}

			stack_destroy( cell_list );
		}

		MPI_Allreduce( cm, cm_total, nDim, MPI_DOUBLE, MPI_SUM, mpi.comm.run );
		MPI_Allreduce( &cm_mass, &cm_mass_total, 1, MPI_DOUBLE, MPI_SUM, mpi.comm.run );

		for ( i = 0; i < nDim; i++ ) {
			cm_total[i] /= cm_mass_total;
		}

		dr = 0.0;
		for ( i = 0; i < nDim; i++ ) {
			dr += cm_total[i]*cm_total[i];
            h->pos[i] += cm_total[i];
			if ( h->pos[i] < 0.0 ) {
				h->pos[i] += (double)num_grid;
			} else if ( h->pos[i] >= (double)num_grid ) {
				h->pos[i] -= (double)num_grid;
			}
        }

		dr = sqrt(dr)/rcm;
		rcm *= 1.0-cm_radius_freduce;
		niter++;

	} while ( dr > cm_convergence_ftol && 
				rcm > cm_convergence_abs*cell_size[max_level] && niter < 100 );

	end_time( HALO_FINDER_RECENTER_TIMER );
}

float *particle_density;
int sort_particles_by_density_desc( const void *a, const void *b ) {
    int index_a = *(int *)a;
    int index_b = *(int *)b;

	if ( particle_density[index_a] < particle_density[index_b] ) {
		return 1;
	} else {
		return -1;
	}
}

struct {
	float density;
	int proc;
} local_particle, global_particle;

halo_list *find_halos() {
	int i, j;
	int ih;
	int num_centers_local;
	int *order;
	int *particle_flag;
	double r;
	halo *h, *h2;
	halo_list *halos;
	int hid;
	double min_halo_mass_code;

#ifdef COSMOLOGY
	min_halo_mass_code = min_halo_mass*constants->Msun*cosmology->h/units->mass;
#else
	min_halo_mass_code = min_halo_mass*constants->Msun/units->mass;
#endif

	start_time( HALO_FINDER_TIMER );

	cart_debug("Finding halos...");

	/* set up binning */
	rlmin = log10( rmin_physical/units->length_in_chimps );
	rlmax = log10( rmax_physical/units->length_in_chimps );
	drl = (rlmax-rlmin)/(double)(num_bins-1);

	rl[0] = 0.0;
	rr[0] = rmin_physical/units->length_in_chimps;
	bin_volume[0] = 4.0*M_PI/3.0 * rr[0]*rr[0]*rr[0];
	bin_volume_cumulative[0] = bin_volume[0];

    for ( i = 1; i < num_bins; i++ ) {
        rl[i] = rr[i-1];
        rr[i] = pow( 10.0, rlmin + (double)i*drl );
        bin_volume[i] = 4.0*M_PI/3.0*( rr[i]*rr[i]*rr[i] - rl[i]*rl[i]*rl[i] );
        bin_volume_cumulative[i] = 4.0*M_PI/3.0*rr[i]*rr[i]*rr[i];
    }

	/* allocate particle_flag and particle_density */
	particle_flag = cart_alloc(int, num_particles);

	/* mark species 0 particles for density calculation */
	num_centers_local = 0;
	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL && 
				particle_id[i] < particle_species_indices[1] ) {
			if ( particle_level[i] >= min_halo_center_level ) {
				particle_flag[i] = 2;
				num_centers_local++;
			} else {
				particle_flag[i] = 1;
			}
		} else {
			particle_flag[i] = 0;
		}
	}

	cart_debug("computing density for %d particles", num_centers_local );

	/* create particle buffer of highest-res dark matter */
	build_particle_buffer( 0, -1 );

	/* mark buffered particles */
	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL && 
				particle_id[i] < particle_species_indices[1] && 
				particle_flag[i] == 0 ) {
			particle_flag[i] = 1;
		}
	}

	particle_density = cart_alloc(float, num_particles);
	compute_particle_densities( halo_num_nearest_neighbors, particle_flag, particle_density );
	destroy_particle_buffer();

	/* apply density threshold and density max constraint for halo centers */
	num_centers_local = 0;
	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_flag[i] == 2 && particle_density[i] > delta_halo_center ) {
			num_centers_local++;
		}
	}

	cart_debug("identified %u potential centers meeting density criterion", num_centers_local );

	order = cart_alloc(int, num_centers_local);

	num_centers_local = 0;
	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_flag[i] == 2 && particle_density[i] > delta_halo_center ) {
			order[num_centers_local++] = i;
		}
	}

	cart_free( particle_flag );

	qsort( order, num_centers_local, sizeof(int), sort_particles_by_density_desc );

	halos = halo_list_alloc(100);

	i = 0;
	hid = 1;
	while (1) {
		/* find next highest density local particle */
		while ( i < num_centers_local && order[i] < 0 ) i++;

		if ( i < num_centers_local ) {
			local_particle.proc = local_proc_id;
			local_particle.density = particle_density[order[i]];
		} else {
			local_particle.proc = -1;
			local_particle.density = 0.0;
		}

		/* gather max density particle */
		MPI_Allreduce( &local_particle, &global_particle, 1, 
				MPI_FLOAT_INT, MPI_MAXLOC, mpi.comm.run );

		if ( global_particle.proc == -1 ) {
			break;
		}

		h = halo_list_add_halo(halos);
		h->id = hid;

		if ( global_particle.proc == local_proc_id ) {
			for ( j = 0; j < nDim; j++ ) {
				h->pos[j] = particle_x[order[i]][j];
			}
			i++;
		}

		/* broadcast position to other processors */
		MPI_Bcast( h->pos, nDim, MPI_DOUBLE, global_particle.proc, mpi.comm.run );

		/* grow sphere */
		compute_halo_mass(h);

		if ( halo_center_definition == 0 ) {
			halo_recenter(h);
			compute_halo_mass(h);
		}

		/* eliminate lower density centers */
#pragma omp parallel for default(none) shared(order,num_centers_local,i,particle_x,h) private(j,r)
		for ( j = i; j < num_centers_local; j++ ) {
			if ( order[j] >= 0 ) {
				r = compute_distance_periodic( particle_x[order[j]], h->pos );
				if ( r < h->rvir ) {
					order[j] = -1;
				}
			}
		}

		if ( h->mvir < min_halo_mass_code ) {
			cart_free(h);
			continue;
		}

		if ( halo_center_definition == 0 ) {
			/* eliminate halos that have wandered into virial radius of previously defined halo */
			for ( ih = 0; ih < halos->num_halos-1; ih++ ) {
				h2 = &halos->list[ih];

				r = compute_distance_periodic( h2->pos, h->pos );
				if ( r < h2->rvir ) {
					halos->num_halos--;
					h = NULL;
					break;
				}
			}

			if ( h == NULL ) {
				continue;
			}
		}

#ifdef DEBUG
		if ( local_proc_id == MASTER_NODE ) {
			cart_debug("halo %d %e %e %e %e %e", 
					hid,
					h->pos[0]*units->length_in_chimps, 
					h->pos[1]*units->length_in_chimps,
					h->pos[2]*units->length_in_chimps,
					h->rvir*units->length_in_chimps, 
					h->mvir*units->mass/constants->Msun );
		}
#endif /* DEBUG */

		/* eliminate lower density centers */
#pragma omp parallel for default(none) shared(order,num_centers_local,i,particle_x,h) private(j,r)
		for ( j = i; j < num_centers_local; j++ ) {
			if ( order[j] >= 0 ) {
				r = compute_distance_periodic( particle_x[order[j]], h->pos );
				if ( r < h->rvir ) {
					order[j] = -1;
				}
			}
		}

		hid++;
	}

	cart_free( order );

	end_time( HALO_FINDER_TIMER );

#ifdef DEBUG
	cart_debug( "mass timer = %f", total_time( HALO_FINDER_MASS_TIMER, -1 ) );
	cart_debug( "cm timer = %f", total_time( HALO_FINDER_RECENTER_TIMER, -1 ) );
	cart_debug( "halo finder = %f", total_time( HALO_FINDER_TIMER, -1 ) );
#endif /* DEBUG */

	return halos;
}

void write_halo_list( halo_list *halos ) {
	int i;
	int ih;
	char filename[512];
	FILE *output;
	halo *h;

	if ( halo_particle_list_flag ) {
		/* also computes Np */
		write_halo_particle_list(halos);
	}

	if ( local_proc_id == MASTER_NODE ) {
#ifdef COSMOLOGY
		sprintf( filename, "%s/halo_catalog_a%6.4f.dat", halo_finder_output_directory, abox[min_level] );
#else
		sprintf( filename, "%s/halo_catalog_%05u.dat", halo_finder_output_directory, step );
#endif /* COSMOLOGY */

		output = fopen( filename, "w" );
		if ( output == NULL ) {
			cart_error("Unable to open halo catalog output file %s", filename );
		}

		fprintf( output, "# %s\n", jobname );
#ifdef COSMOLOGY
		fprintf( output, "# step = %u, auni = %8.6f, abox = %8.6f\n", step, auni[min_level], abox[min_level] );
		fprintf( output, "# Cosmology: OmM = %.2f, OmL = %.2f, OmB = %.4f, h = %.2f, DeltaDC = %.3f\n",
				cosmology->OmegaM, cosmology->OmegaL, cosmology->OmegaB, cosmology->h, cosmology->DeltaDC );
#endif /* COSMOLOGY */
		fprintf( output, "# Lbox = %.2f [Mpc/h comoving]\n", box_size );
		fprintf( output, "# num_neighbors = %u, Deltamin = %.2f, Deltavir = %.2f (mean)\n", 
			halo_num_nearest_neighbors, delta_halo_center, delta_vir);
		fprintf( output, "# Binning: rmin = %.3f, rmax = %.3f [kpc/h comoving], num_bins = %u\n", rmin_physical, rmax_physical, num_bins );

		if ( halo_center_definition == 0 ) {
			fprintf( output, "# Halo centering: center of mass, freduce = %.3f, conv. ftol = %.2e, abs = %.2f [kpc/h comoving]\n", 
					cm_radius_freduce, cm_convergence_ftol, cm_convergence_abs );	
		} else {
			fprintf( output, "# Halo centering: density peak\n" );
		}

#ifdef COSMOLOGY
		fprintf( output, "# Columns: halo_id x y z [Mpc/h comoving] vx vy vz [peculiar km/s] Rvir [kpc/h comoving]\n#     Mvir [Msun/h] vmax [km/s] rmax [kpc/h comoving]\n" );
#else
		fprintf( output, "# Columns: halo_id x y z [Mpc] vx vy vz [peculiar km/s] Rvir [kpc]\n#        Mvir [Msun] Np [specie 0] vmax [km/s] rmax [kpc]\n" );
#endif
		fprintf( output, "#################################################################################################################\n");

		for ( ih = 0; ih < halos->num_halos; ih++ ) {
			h = &halos->list[ih];
#ifdef COSMOLOGY
			fprintf( output, "%5u %10.5f %10.5f %10.5f %8.2f %8.2f %8.2f %9.4f %.5e %7u %7.2f %9.4f\n", 
					h->id, 
					h->pos[0]*units->length_in_chimps,
					h->pos[1]*units->length_in_chimps,
					h->pos[2]*units->length_in_chimps,
					h->vel[0]*units->velocity/constants->kms,
					h->vel[1]*units->velocity/constants->kms,
					h->vel[2]*units->velocity/constants->kms,
					h->rvir*units->length_in_chimps*1000.,
					h->mvir*cosmology->h*units->mass/constants->Msun,
					h->np,
					h->vmax*units->velocity/constants->kms,
					h->rmax*units->length_in_chimps*1000. );
#else
			fprintf("%5u %10.5f %10.5f %10.5f %8.2f %8.2f %8.2f %9.4f %.5e %7u %7.2f %9.4f\n", h->id, 
					h->pos[0]*units->length / constants->Mpc,
					h->pos[1]*units->length / constants->Mpc,
					h->pos[2]*units->length / constants->Mpc,
					h->vel[0]*units->velocity/constants->kms,
					h->vel[1]*units->velocity/constants->kms,
					h->vel[2]*units->velocity/constants->kms,
					h->rvir*units->length / constants->kpc,
					h->mvir*units->mass/constants->Msun,
					h->np,
					h->vmax*units->velocity/constants->kms,
					h->rmax*units->length / constants->kpc );
#endif
		}

		fclose(output);
	}
}

void write_halo_particle_list( halo_list *halos ) {
    int i, j, k, m;
	int ih;
    char filename[512];
    FILE *output;
    halo *h;
	stack *cell_list;
	int coords[nDim];
	double rvir2;
	float aexp;
	int size;
	int icell, ioct, ipart;
	int bin, level, parent, child;
	int nh, np;
	int *particle_counts;
	int thread_particle_count;
	int local_particle_count;
	int total_particle_count;
	int proc;
	int *ids;
	float *bind;
	double dx, r, v, lgr;
	double phi, v2kms2, phi2kms2;
	int plocal, pindex;
	double thread_radial_potential[MAX_HALO_BINS];
	double thread_bin_volume[MAX_HALO_BINS];
	double local_radial_potential[MAX_HALO_BINS];
	double local_bin_volume[MAX_HALO_BINS];
	double radial_potential[MAX_HALO_BINS];
	double actual_bin_volume[MAX_HALO_BINS];

	start_time( HALO_FINDER_WRITE_PARTICLES_TIMER );

	if ( local_proc_id == MASTER_NODE ) {
#ifdef COSMOLOGY
		sprintf( filename, "%s/halo_particles_a%6.4f.dat", halo_finder_output_directory, abox[min_level] );
#else
		sprintf( filename, "%s/halo_particles_%05u.dat", halo_finder_output_directory, step );
#endif /* COSMOLOGY */

		output = fopen( filename, "w" );

		size = sizeof(float);
		aexp = auni[min_level];

		fwrite( &size, sizeof(int), 1, output );
		fwrite( &aexp, sizeof(float), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		nh = halos->num_halos;
		np = particle_species_num[0];
		
		size = 2*sizeof(int);
		fwrite( &size, sizeof(int), 1, output );
		fwrite( &nh, sizeof(float), 1, output );
		fwrite( &np, sizeof(int), 1, output );
		fwrite( &size, sizeof(int), 1, output );

		particle_counts = cart_alloc(int, num_procs);
	}

	for ( ih = 0; ih < halos->num_halos; ih++ ) {
		h = &halos->list[ih];
		rvir2 = h->rvir*h->rvir;

		for ( bin = 0; bin < num_bins; bin++ ) {
			local_bin_volume[bin] = 0.0;
			local_radial_potential[bin] = 0.0;
		}
		local_particle_count = 0;

		/* construct radial average potential */
#pragma omp parallel default(none) shared(h,cell_child_oct,cell_vars,num_bins,cell_volume,oct_pos,oct_level,cell_delta,rvir2,particle_id,particle_x,cell_particle_list,particle_list_next,rlmin,rr,drl,local_radial_potential,local_bin_volume,local_particle_count) private(cell_list,i,j,k,m,coords,r,icell,dx,level,ioct,parent,child,lgr,bin,thread_radial_potential,thread_bin_volume,thread_particle_count,ipart)
        {
            cell_list = stack_init();

			for ( bin = 0; bin < num_bins; bin++ ) {
				thread_radial_potential[bin] = 0.0;
				thread_bin_volume[bin] = 0.0;
			}
			thread_particle_count = 0;

			/* construct radial average potential */
#pragma omp for collapse(3) schedule(dynamic) nowait
			for ( i = (int)floor(h->pos[0]-rr[num_bins-1]); i <= (int)(h->pos[0]+rr[num_bins-1]); i++ ) {
				for ( j = (int)floor(h->pos[1]-rr[num_bins-1]); j <= (int)(h->pos[1]+rr[num_bins-1]); j++ ) {
					for ( k = (int)floor(h->pos[2]-rr[num_bins-1]); k <= (int)(h->pos[2]+rr[num_bins-1]); k++ ) {
						coords[0] = ( i + num_grid ) % num_grid;
						coords[1] = ( j + num_grid ) % num_grid;
						coords[2] = ( k + num_grid ) % num_grid;
						icell = root_cell_location( sfc_index( coords ) );
						if ( icell != NULL_OCT && cell_is_local(icell) ) {
							stack_push( cell_list, icell );

							while ( stack_pop( cell_list, &icell ) ) {
								if ( cell_is_root_cell(icell) ) {
									level = min_level;
									r = 0.0;
									for ( m = 0; m < nDim; m++ ) {
										dx = (double)coords[m] + 0.5 - h->pos[m];
										if ( dx < -num_grid/2 ) {
											dx += num_grid;
										} else if ( dx > num_grid/2 ) {
											dx -= num_grid;
										}
										r += dx*dx;
									}
								} else {
									parent = cell_parent_oct(icell);
									level = oct_level[parent];
									child = cell_child_number(icell);

									r = 0.0;
									for ( m = 0; m < nDim; m++ ) {
										dx = (double)oct_pos[parent][m] + 
												(double)cell_size[level]*(double)cell_delta[child][m] - h->pos[m];
										if ( dx < -num_grid/2 ) {
											dx += num_grid;
										} else if ( dx > num_grid/2 ) {
											dx -= num_grid;
										}
										r += dx*dx;
									}
								}

								lgr = 0.5*log10(r);
								bin = ( lgr < rlmin ) ? 0 : (int)((lgr - rlmin)/drl) + 1;

								if ( cell_is_refined(icell) ) { 
									ioct = cell_child_oct[icell];
									for ( m = 0; m < num_children; m++ ) {
										stack_push( cell_list, oct_child( ioct, m ) );
									}
								} else {
									if ( bin < num_bins ) {
										thread_radial_potential[bin] += cell_potential(icell)*cell_volume[level];
										thread_bin_volume[bin] += cell_volume[level];
									}
								}

								ipart = cell_particle_list[icell];
								while ( ipart != NULL_PARTICLE ) {
									if ( particle_species( particle_id[ipart] ) == 0 ) {
										r = 0.0;
										for ( m = 0; m < nDim; m++ ) {
											dx = particle_x[ipart][m] - h->pos[m];
											if ( dx < -num_grid/2 ) {
												dx += num_grid;
											} else if ( dx > num_grid/2 ) {
												dx -= num_grid;
											}
											r += dx*dx;
										} 

										if ( r < rvir2 ) {
											thread_particle_count++;
										}
									}
									ipart = particle_list_next[ipart];
								}
							}
						}
					}
				}
			}

#pragma omp critical
			{
				for ( bin = 0; bin < num_bins; bin++ ) {
					local_radial_potential[bin] += thread_radial_potential[bin];
					local_bin_volume[bin] += thread_bin_volume[bin];
				}
				local_particle_count += thread_particle_count;
			}

            stack_destroy( cell_list );
		}

		MPI_Allreduce( local_radial_potential, radial_potential, num_bins, MPI_DOUBLE, MPI_SUM, mpi.comm.run );
		MPI_Allreduce( local_bin_volume, actual_bin_volume, num_bins, MPI_DOUBLE, MPI_SUM, mpi.comm.run );

		for ( bin = 0; bin < num_bins; bin++ ) {
			if ( actual_bin_volume[bin] > 0 ) {
				radial_potential[bin] /= actual_bin_volume[bin];
			} else {
				radial_potential[bin] = 0.0;
			}
		}

		/* ensure non-zero values, assumes possible 0 volumes more likely at low radius */
		for ( bin = num_bins-1; bin >= 0; bin-- ) {
			if ( actual_bin_volume[bin] == 0.0 ) {
				if ( bin == num_bins-1 ) {
					i = num_bins-2;
					while ( i >= 0 && radial_potential[i] == 0.0 ) i--;
					if ( i == 0 ) {
						cart_error("Potential appears to be 0, did you forget to solve the poisson equation?");
					}
					radial_potential[bin] = radial_potential[i];
				} else {
					radial_potential[bin] = radial_potential[bin+1];
				}
			}
		}

		/* measure binding energy for local particles */
		ids = cart_alloc(int, local_particle_count);
		bind = cart_alloc(float, local_particle_count);
		pindex = 0;
		v2kms2 = pow( units->velocity/constants->kms, 2.0 );
		phi2kms2 = pow( units->velocity*abox[min_level]/constants->kms, 2.0 );

#pragma omp parallel default(none) shared(h,cell_child_oct,cell_vars,num_bins,cell_volume,oct_pos,oct_level,cell_delta,rvir2,particle_x,particle_id,particle_v,particle_list_next,cell_particle_list,rlmin,drl,pindex,ids,bind,radial_potential,rl,phi2kms2,v2kms2) private(cell_list,i,j,k,m,coords,r,lgr,v,icell,dx,level,ioct,parent,child,plocal,ipart,phi,bin)
		{
			cell_list = stack_init();

			/* construct radial average potential */
#pragma omp for collapse(3) schedule(dynamic) nowait
			for ( i = (int)floor(h->pos[0]-h->rvir); i <= (int)(h->pos[0]+h->rvir); i++ ) {
				for ( j = (int)floor(h->pos[1]-h->rvir); j <= (int)(h->pos[1]+h->rvir); j++ ) {
					for ( k = (int)floor(h->pos[2]-h->rvir); k <= (int)(h->pos[2]+h->rvir); k++ ) {
						coords[0] = ( i + num_grid ) % num_grid;
						coords[1] = ( j + num_grid ) % num_grid;
						coords[2] = ( k + num_grid ) % num_grid;
						icell = root_cell_location( sfc_index( coords ) );
						if ( icell != NULL_OCT && cell_is_local(icell) ) {
							stack_push( cell_list, icell );

							while ( stack_pop( cell_list, &icell ) ) {
								if ( cell_is_refined(icell) ) { 
									ioct = cell_child_oct[icell];
									for ( m = 0; m < num_children; m++ ) {
										stack_push( cell_list, oct_child( ioct, m ) );
									}
								}

								ipart = cell_particle_list[icell];
								while ( ipart != NULL_PARTICLE ) {
									if ( particle_species( particle_id[ipart] ) == 0 ) {
										r = 0.0;
										for ( m = 0; m < nDim; m++ ) {
											dx = particle_x[ipart][m] - h->pos[m];
											if ( dx < -num_grid/2 ) {
												dx += num_grid;
											} else if ( dx > num_grid/2 ) {
												dx -= num_grid;
											}
											r += dx*dx;
										}

										lgr = 0.5*log10(r);
										bin = ( lgr < rlmin ) ? 0 : (int)((lgr - rlmin)/drl) + 1;

										if ( r < rvir2 && bin < num_bins ) {
											/* bad: introduces false cache sharing between threads, oh well */
											#pragma omp critical 
											{
												plocal = pindex;
												pindex++;	
											}
											ids[plocal] = particle_id[ipart];

											v = 0.0;
											for ( m = 0; m < nDim; m++ ) {
												v += (particle_v[ipart][m]-h->vel[m])*(particle_v[ipart][m]-h->vel[m]);
											}

											if ( bin == 0 ) {
												phi = radial_potential[0];
											} else {
												phi = loglin_interpolate( radial_potential, bin, lgr, 1.0/drl, log10(rl[bin]) );
											}

											/* phi = cell_potential(icell); */
											bind[plocal] = phi*phi2kms2 + 0.5*v*v2kms2;
										}
									}
									ipart = particle_list_next[ipart];
								}
							}
						}
					}
				}
			}

			stack_destroy( cell_list );
		}

		cart_assert( pindex == local_particle_count );

		/* gather binding energies */
		MPI_Gather( &local_particle_count, 1, MPI_INT, particle_counts, 1, MPI_INT, MASTER_NODE, mpi.comm.run );

		if ( local_proc_id == MASTER_NODE ) {
			total_particle_count = 0;
			for ( proc = 0; proc < num_procs; proc++ ) {
				total_particle_count += particle_counts[proc];
			}

			h->np = total_particle_count;
			
			size = total_particle_count*sizeof(int) + 
					total_particle_count*sizeof(float) +
					2*sizeof(int);

			fwrite( &size, sizeof(int), 1, output );
			fwrite( &h->id, sizeof(int), 1, output );
			fwrite( &total_particle_count, sizeof(int), 1, output );
			fwrite( ids, sizeof(int), local_particle_count, output );
			cart_free( ids );

			for ( proc = 1; proc < num_procs; proc++ ) {
				ids = cart_alloc(int, particle_counts[proc] );
				MPI_Recv( ids, particle_counts[proc], MPI_INT, proc, 0, mpi.comm.run, MPI_STATUS_IGNORE );
				fwrite( ids, sizeof(int), particle_counts[proc], output );
				cart_free( ids );
			}

			fwrite( bind, sizeof(float), local_particle_count, output );
			cart_free( bind );

			for ( proc = 1; proc < num_procs; proc++ ) {
				bind = cart_alloc(float, particle_counts[proc]);
				MPI_Recv( bind, particle_counts[proc], MPI_FLOAT, proc, 0, mpi.comm.run, MPI_STATUS_IGNORE );
				fwrite( bind, sizeof(float), particle_counts[proc], output );
				cart_free( bind );
			}

			fwrite( &size, sizeof(int), 1, output );
		} else {
			/* send binding energies to root processor */
			MPI_Send( ids, local_particle_count, MPI_INT, MASTER_NODE, 0, mpi.comm.run );
			MPI_Send( bind, local_particle_count, MPI_FLOAT, MASTER_NODE, 0, mpi.comm.run );

			cart_free( ids );
			cart_free( bind );
		}
	}

	if ( local_proc_id == MASTER_NODE ) {
		cart_free( particle_counts );
		fclose(output);
	}

	end_time( HALO_FINDER_WRITE_PARTICLES_TIMER );
}

#endif /* PARTICLES */
