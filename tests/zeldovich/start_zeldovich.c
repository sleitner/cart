#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include <sys/types.h>
#include <unistd.h>

#include "defs.h"
#include "auxiliary.h"
#include "tree.h"
#include "particle.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "timestep.h"
#include "refinement.h"
#include "refinement_operations.h"
#include "viewdump.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "gravity.h"
#include "density.h"
#include "io.h"

#define a_cross		(1.0/(1.0+10.0))

#define lambda		((double)num_grid)

double p0;
double rhogas0;
double dgrowth;
double ddgrowthdt;
double ampl;
double ak;

void compute_omegas( double a, double *Omega, double *OmegaL ) {
	double f;
                                                                                                                                                            
	f = a + Omega0 * (1.0-a) + OmegaL0 * (a*a*a - a);
                                                                                                                                                            
	*Omega = Omega0 / f;
	*OmegaL = OmegaL0 / f;
}


double g_CPT( double Omega, double OmegaL ) {
	return 0.25* Omega / (
		pow( Omega, 4.0/7.0 ) -
		OmegaL +
		(1.0-0.5*Omega)*(1.0+OmegaL/70.0) );
}

double growth( double a ) 
/* Calculates growth factor using formula
 *  from Caroll, Press & Turner 1992, ARA&A 30, 499
 */
{
	double Omega, OmegaL;

	if ( Omega0 == 1.0 && OmegaL0 == 0.0 ) {
		return a;
	} else {
		compute_omegas( a, &Omega, &OmegaL );

		return a * g_CPT( Omega, OmegaL ) / 
			g_CPT( Omega0, OmegaL0 );
	}
}

double dgrowthdt( double a ) {
	return 2.0 * sqrt( a*a*a );
}

double x0;

double qsolve( double q ) {
	double x;

	x = q + dgrowth*ampl*sin(ak*q);

	return x-x0;
}

#ifdef HYDRO
void initial_conditions( int cell, int level ) {
	double q;
	float pos[nDim];

	double q1, q2;

	cell_position( cell, pos );

	x0 = pos[0];
	q = root_finder( qsolve, 0.0, num_grid, 1e-9, 1e-9 );

	x0 = pos[0] - 0.5*cell_size[level];
	cart_assert( x0 >= 0.0 && x0 < num_grid );
	q1 = root_finder( qsolve, 0.0, num_grid, 1e-9, 1e-9 );
	
	x0 = pos[0] + 0.5*cell_size[level];
	if ( x0 == num_grid ) {
		q2 = num_grid;
	} else {
		q2 = root_finder( qsolve, 0.0, num_grid, 1e-9, 1e-9 );
	}

	/* cell_gas_density(cell) = rhogas0*(q2-q1)/cell_size[level]; */
	cell_gas_density(cell) = rhogas0 / ( 1.0 + ak*dgrowth*ampl*cos(ak*q) );
	cart_assert( cell_gas_density(cell) > 0.0 );

	cell_momentum(cell,0) = cell_gas_density(cell)*ddgrowthdt * ampl * sin( ak * q );
        cell_momentum(cell,1) = 0.0;
        cell_momentum(cell,2) = 0.0;
        cell_gas_gamma(cell) = (5.0/3.0);

	cell_gas_pressure(cell) = p0;
	cell_gas_internal_energy(cell) = p0/(cell_gas_gamma(cell)-1.0);
	cell_gas_energy(cell) = cell_gas_internal_energy(cell) + 
		0.5*(cell_momentum(cell,0)*cell_momentum(cell,0))/cell_gas_density(cell);
}
#endif /* HYDRO */

double particle_q_init( int id ) {
	int i, j, k;
	double qfact;

	k = id % num_row;
	j = ((id-k)/num_row) % num_row;
	i = (((id-k)/num_row) - j ) / num_row;

	qfact = (double)num_grid/(double)num_row;

	return qfact*((double)i + 0.5);
}

FILE *output;
FILE *particles;
double c1,c2,c3;
double phi0;

void output_cell( int cell, int level ) {
	float pos[nDim];
	double q;
	double kq;
	double rho, v, g, phi;
	double x, x_a;
	double aexp_a, dgrowth_a;

#ifdef PARTICLES
	int ipart;
	int npart;
	double q_i;
	double x_z;
	double v_x;
	double rho_x;
#endif

	if ( cell_is_leaf(cell) ) {
		cell_position(cell,pos);

		x0 = pos[0];
		q = root_finder( qsolve, 0.0, num_grid, 1e-9, 1e-9 );
		kq = ak*q;

		if ( tl_old[min_level] == 0.0 ) {
			aexp_a = b2a( tl[level] - 0.5*dtl[level] );
		} else {
			aexp_a = b2a( tl[level] + 0.5*dtl[level] );
		}
		dgrowth_a = growth(aexp_a);

		x_a = q + dgrowth_a*ampl*sin(kq);
                                                                                                                                                            
		rho = Omega0 / ( 1.0 + ak*dgrowth*ampl*cos(kq) );
		v = ddgrowthdt*ampl*sin(kq);
		g = -6.0*aexp_a*( rhogas0*q - x_a );

		phi = c1*(c2*(0.5*kq*kq + c3*(kq*sin(kq)+cos(kq)-1.0)) - 0.5*(x0*x0)) + phi0;

#ifdef HYDRO
/*
		fprintf(output, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %u\n",
			q, x0, 
			(cell_density(cell)/cell_volume[level]+1.0), rho, fabs(rho-(cell_density(cell)/cell_volume[level])+1.0)/rho,
			cell_momentum(cell,0)/cell_gas_density(cell), v, fabs(v-cell_momentum(cell,0)/cell_gas_density(cell))/v,
			cell_accel(cell, 0)/dtl[level], g, fabs(g-cell_accel(cell, 0)/dtl[level])/g,
			cell_potential(cell), phi, fabs(phi-cell_potential(cell))/cell_potential(cell), 
			cell_gas_energy(cell), cell_gas_pressure(cell),
			cell_gas_internal_energy(cell), cell_level(cell) );
*/
		fprintf(output, "%e %e %u %e %e %e %e %e %e %e %e %u %e\n", q, x0, cell_level(cell),
			cell_density(cell), cell_potential(cell), cell_potential_hydro(cell),
			cell_gas_density(cell), cell_momentum(cell,0)/cell_gas_density(cell),
			cell_gas_pressure(cell), cell_momentum(cell,0), ref[cell], cell, phi );
#else 
/*
		fprintf(output, "%e %e %e %e %e %e %e %e %e %e %e\n",
			q, x0,
			(cell_density(cell)+cell_volume[level]), rho, fabs(rho-(cell_density(cell)+cell_volume[level]))/rho,
			cell_accel(cell, 0), g, fabs(g-cell_accel(cell, 0))/g,
			cell_potential(cell), phi, fabs(phi-cell_potential(cell))/cell_potential(cell) );
*/
/*		fprintf(output, "%e %e %e %e\n", q, x0, cell_density(cell), cell_potential(cell) ); */
#endif /* HYDRO */

#ifdef PARTICLES
/*
		ipart = cell_particle_list[cell];

		dgrowth = growth(aexp[min_level]);
		ddgrowthdt = dgrowthdt( b2a( tl[min_level] - 0.5*dtl[min_level] ) );

		while ( ipart != NULL_PARTICLE ) {
			x0 = particle_x[ipart][0];
			q = root_finder( qsolve, 0.0, num_grid, 1e-9, 1e-9 );
			kq = ak*q;

			q_i = particle_q_init( particle_id[ipart] );

			x_z = q_i + dgrowth*ampl*sin(ak*q_i);
			v = ddgrowthdt*ampl*sin(ak*q_i);

			fprintf(particles, "%u %e %e %e %e %e %e %e %e %e\n", ipart,
				particle_x[ipart][0], particle_x[ipart][1], particle_x[ipart][2],
				particle_v[ipart][0], particle_v[ipart][1], particle_v[ipart][2],
				particle_a[ipart][0], particle_a[ipart][1], particle_a[ipart][2] );

			ipart = particle_list_next[ipart];
		}
*/
#endif /* PARTICLES */
	}
}

void run_output() {
	int i, j;
	int ipart;
	char filename[128];
	double pos[nDim];
	double phi;
	double x_rms, x_norm;
	double v_rms, v_norm;
	double x_rms_total, v_rms_total;
	double x_norm_total, v_norm_total;
	double x_z, v_z;
	double q, q_i;
	int num_level_cells;
	int *level_cells;

	sprintf(filename, "dumps/zeldovich_[%5.4f]_%03d.dat", aexp[min_level], local_proc_id );
	output = fopen(filename, "w");

	sprintf(filename, "dumps/zeldovich_particles_[%5.4f]_%03d.dat", aexp[min_level], local_proc_id );
	particles = fopen(filename, "w");

	dgrowth = growth(aexp[min_level]);
        ddgrowthdt = dgrowthdt(aexp[min_level]);

	c1 = 6.0/aexp[min_level];
	c2 = 1.0 / (ak*ak);
	c3 = dgrowth*ampl*ak;

	pos[1] = pos[2] = (double)num_grid / 2.0 - 0.5*cell_size[max_level];
	pos[0] = 0.0;
	i = cell_find_position(pos);
	if ( i != -1 ) {
		phi = cell_potential(i);
	} else { 
		phi = 0.0;
	}

	MPI_Allreduce( &phi, &phi0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

	/* uncommment to print all cells */
/*
	for ( i = min_level; i <= max_level; i++ ) {
		select_level( i, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( j = 0; j < num_level_cells; j++ ) {
			output_cell( level_cells[j], i );
		}
		cart_free(level_cells);
	}
*/

	for ( pos[0] = 0.0; pos[0] < num_grid; pos[0] += 0.5*cell_size[max_level] ) {
		i = cell_find_position(pos);

		if ( i != -1 ) {
			output_cell( i, cell_level(i) );
		}
	}

	fclose(output);

	/* compute rms fluctuations */
	x_rms = x_norm = 0.0;
	v_rms = v_norm = 0.0;

	ddgrowthdt = dgrowthdt( b2a( tl[min_level] - 0.5*dtl[min_level] ) );

	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL ) {
			ipart = i;
			x0 = particle_x[i][0];
			q = root_finder( qsolve, 0.0, num_grid, 1e-9, 1e-9 );

			fprintf(particles, "%u %e %e %e %e %e %e %e\n", particle_id[i],
				particle_x[i][0], particle_x[i][1], particle_x[i][2],
				particle_v[i][0], particle_v[i][1], particle_v[i][2],
				particle_pot[i] );

			q_i = particle_q_init( particle_id[i] );

			x_z = q_i + dgrowth*ampl*sin(ak*q_i);
			v_z = ddgrowthdt*ampl*sin(ak*q_i);

			x_rms += (x0-x_z)*(x0-x_z);
			x_norm += (x_z-q_i)*(x_z-q_i);

			v_rms += (particle_v[i][0]-v_z)*(particle_v[i][0]-v_z);
			v_norm += v_z*v_z;
		}
	}

	fclose(particles);

	MPI_Reduce( &x_rms, &x_rms_total, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( &v_rms, &v_rms_total, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( &x_norm, &x_norm_total, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
	MPI_Reduce( &v_norm, &v_norm_total, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );

	if ( local_proc_id == MASTER_NODE ) {	

		x_rms_total /= x_norm_total;
		v_rms_total /= v_norm_total;

		x_rms_total = sqrt(x_rms_total);
		v_rms_total = sqrt(v_rms_total);

		particles = fopen("dumps/particle_rms.dat", "a");
		fprintf(particles, "%e %e %e %e %e\n", aexp[min_level], x_rms_total, v_rms_total, x_norm_total, v_norm_total );
		fclose(particles);
	}
}

void init_run() {
	int i, j, k;
	int index;
	double t, a;
	double a_th;
	double TinitK, Tinit;

	int ipart;
	int icell;
	double qi, qj, qk;
	double xcons, vcons;
	double dx, dvx;
	float pw;
	double a_vel;
	double qfact;

	int num_level_cells;
	int *level_cells;

        build_cell_buffer();
        repair_neighbors();

	a = a_init;
	t = a2b( a_init );

	for ( i = min_level; i <= max_level; i++ ) {
                tl[i] = t;
                aexp[i] = a;
	}

	rhogas0 = Omegab0/Omega0;
	cart_debug("rhogas0 = %e", rhogas0 );

	/* calculate temperature */
	a_th = 1.0 / ( 1e3 * pow(Omegab0*hubble*hubble, 0.4) );

	if ( a < a_th ) {
		TinitK = 2.726 / a;
	} else {
		TinitK = 2.726 / a_th * (a_th/a)*(a_th/a);
	}

	Tinit = TinitK * a*a / (0.31*T0/0.6);

#ifdef PRESSURELESS_FLUID
	p0 = 1e-20;
#else
	//p0 = Tinit * Omega0 * (Omegab0/Omega0);
	p0 = 3.2e-9;
#endif

	ak = 2.0*M_PI / lambda;
	dgrowth = growth(a);
	ddgrowthdt = dgrowthdt(a);
	ampl = 1.0 / ( growth(a_cross) * ak );

#ifdef HYDRO
	for ( i = min_level; i <= max_level; i++ ) {
		cart_debug("generating initial conditions on level %u", i );
	
		select_level( i, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( j = 0; j < num_level_cells; j++ ) {
			initial_conditions( level_cells[j], i );
		}
		cart_free( level_cells );	
	}

	for ( i = min_level; i <= max_level; i++ ) {
		update_buffer_level( i, all_hydro_vars, num_hydro_vars );
	}
#endif /* HYDRO */

        dtl[min_level] = 0.0;
        choose_timestep( &dtl[min_level] );

#ifdef PARTICLES
	num_row = num_grid;
	qfact = (double)num_grid / (double)num_row;
	pw = (1.0-rhogas0)*qfact*qfact*qfact;

	cart_debug("particle weight = %e", pw );

	xcons = dgrowth*ampl;
	a_vel = b2a( t - 0.5*dtl[min_level] );
	vcons = ampl * dgrowthdt( a_vel );

	ipart = 0;
	for ( i = 0; i < num_row; i++ ) {
		qi = qfact*((double)i + 0.5);
		dx = xcons * sin( ak * qi );
		dvx = vcons * sin( ak * qi );
		for ( j = 0; j < num_row; j++ ) {
			qj = qfact*((double)j + 0.5);
			for ( k = 0; k < num_row; k++ ) {
				qk = qfact*((double)k + 0.5);

				particle_x[ipart][0] = qi + dx;
				particle_x[ipart][1] = qj;
				particle_x[ipart][2] = qk;

				if ( particle_x[ipart][0] >= (double)num_grid ) {
					particle_x[ipart][0] -= num_grid;
				}

				if ( particle_x[ipart][1] >= (double)num_grid ) {
					particle_x[ipart][1] -= num_grid;
				}

				if ( particle_x[ipart][2] >= (double)num_grid ) {
                                      	particle_x[ipart][2] -= num_grid;
                                }

				icell = cell_find_position( particle_x[ipart] );

				if ( icell != -1 && cell_is_local(icell) ) {
					particle_v[ipart][0] = dvx;
					particle_v[ipart][1] = 0.0;
					particle_v[ipart][2] = 0.0;

					particle_id[ipart] = num_row*num_row*i + num_row*j + k;

					cart_assert( qi == particle_q_init( particle_id[ipart] ) );

					particle_t[ipart] = t;
					particle_dt[ipart] = dtl[min_level];

					ipart++;
				}
			}
		}
	}

	cart_debug("created %u particles", ipart );

	num_local_particles = ipart;
	num_particles_total = num_row*num_row*num_row;
	num_particle_species = 1;
	particle_species_mass[0] = pw;
	particle_species_num[0] = num_particles_total;
	particle_species_indices[0] = 0;
	particle_species_indices[1] = num_particles_total;

	build_particle_list();

/*
	assign_density( min_level, min_level );
	modify( min_level, 0 );
*/

	if ( local_proc_id == MASTER_NODE ) {
		particles = fopen("dumps/particle_rms.dat", "w");
		fclose(particles);
	}

	

#endif

	for ( i = min_level+1; i < max_level; i++ ) {
		dtl[i] = 0.5*dtl[i-1];
	}

	cart_debug("done with initialization");

	check_map();
}
