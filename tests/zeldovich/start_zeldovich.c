#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cosmology.h"
#include "density.h"
#include "gravity.h"
#include "hydro.h"
#include "io.h"
#include "iterators.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "refinement_operations.h"
#include "sfc.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "run/step.h"

#include "extra/viewdump.h"

/* setup as in Springel 2009 */
/* No shocking until a_cross */
/* #define a_cross		(0.50)  */
#define a_cross		(0.025)
#define OmM0            1.0
#define OmB0            1.0
#define OmL0            0.0
#define h0              0.5
#define dDC             0.0
#define Lbox0           64.0
#define TinitK          (100.0)

#define lambda		((double)num_grid)

double Tinit;
double rhogas0;
double dgrowth;
double ddgrowthdt;
double ampl;
double ak;

double dgrowth0;
void compute_omegas( double a, double *OmegaM, double *OmegaL ) {
	double f;
                                                                                                                                                            
	f = a + cosmology->OmegaM * (1.0-a) + cosmology->OmegaL * (a*a*a - a);
                                                                                                                                                            
	*OmegaM = cosmology->OmegaM / f;
	*OmegaL = cosmology->OmegaL / f;
}


double g_CPT( double OmegaM, double OmegaL ) {
	return 0.25* OmegaM / (
		pow( OmegaM, 4.0/7.0 ) -
		OmegaL +
		(1.0-0.5*OmegaM)*(1.0+OmegaL/70.0) );
}

double growth( double a ) 
/* Calculates growth factor using formula
 *  from Caroll, Press & Turner 1992, ARA&A 30, 499
 */
{
	double OmegaM, OmegaL;

	if ( cosmology->OmegaM == 1.0 && cosmology->OmegaL == 0.0 ) {
		return a;
	} else {
		compute_omegas( a, &OmegaM, &OmegaL );

		return a * g_CPT( OmegaM, OmegaL ) / 
			g_CPT( cosmology->OmegaM, cosmology->OmegaL );
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
#ifdef GRAVITY
void initial_gravity( int cell, int level ) {
	double q;
        double grav_accel, abox_a, x_a, dgrowth_a,kq ;
	double q1, q2;
	double g1, g2, g3, phi;
	double pos[nDim];
        
	cell_center_position( cell, pos );

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
        
        
        kq = ak*q;
        abox_a = abox_from_tcode( tl[level] - 0.5*dtl[level] );
        dgrowth_a = growth(abox_a);
        x_a = q + dgrowth_a*ampl*sin(kq);
        grav_accel = -6.0*abox_a*( rhogas0*q - x_a )*dtl[level];
        
	g1 = 6.0/abox[min_level]; 
	g2 = 1.0 / (ak*ak);
	g3 = dgrowth*ampl*ak;
        phi = g1*(g2*(0.5*kq*kq + g3*(kq*sin(kq)+cos(kq)-1.0)) - 0.5*(x0*x0)); //+ phi0; 
        
 	cell_accel(cell,0) = grav_accel;
 
 	cell_potential(cell) = phi;
 	cell_potential_hydro(cell) = phi;
 	cell_accel(cell,1) = 0;
 	cell_accel(cell,2) = 0;
}
#endif /* GRAVITY */

void initial_conditions( int cell, int level ) {
	double q;
	double q1, q2;
	double pos[nDim];
        
	cell_center_position( cell, pos );

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
        
        
/*         cell_gas_density(cell) = rhogas0*(q2-q1)/cell_size[level];  */
 	cell_gas_density(cell) = rhogas0 / ( 1.0 + ak*dgrowth*ampl*cos(ak*q) );
 	cell_total_mass(cell) = rhogas0 / ( 1.0 + ak*dgrowth*ampl*cos(ak*q) )*cell_volume[level];
	cart_assert( cell_gas_density(cell) > 0.0 );

 	cell_momentum(cell,0) = cell_gas_density(cell)*ddgrowthdt * ampl * sin( ak * q ); 
        cell_momentum(cell,1) = 0.0;
        cell_momentum(cell,2) = 0.0;
        cell_gas_gamma(cell) = (constants->gamma);

 	cell_gas_pressure(cell) = Tinit*cell_gas_density(cell)/constants->wmu;
 	cell_gas_internal_energy(cell) = cell_gas_pressure(cell)/(cell_gas_gamma(cell)-1.0);  
        
	cell_gas_energy(cell) = cell_gas_internal_energy(cell) + cell_gas_kinetic_energy(cell);
}
#endif /* HYDRO */

#ifdef PARTICLES
double particle_q_init( int id ) {
	int i, j, k;
	double qfact;

	k = id % num_row;
	j = ((id-k)/num_row) % num_row;
	i = (((id-k)/num_row) - j ) / num_row;

	qfact = (double)num_grid/(double)num_row;

	return qfact*((double)i + 0.5);
}
#endif

FILE *output;
FILE *particles;
double c1,c2,c3;
double phi0;

void output_cell( int cell, int level ) {
	double q;
	double kq;
	double rho0, rho, v, grav_accel, phi, TK, cell_tempPK, cell_tempK;
	double x, x_a;
	double abox_a, dgrowth_a;
        double pos[nDim];
        

#ifdef PARTICLES
	int ipart;
	int npart;
	double q_i;
	double x_z;
	double v_x;
	double rho_x;
#endif

	if ( cell_is_leaf(cell) ) {

                cell_center_position( cell, pos );

		x0 = pos[0];
		q = root_finder( qsolve, 0.0, num_grid, 1e-9, 1e-9 );
		kq = ak*q;

		if ( tl_old[min_level] == 0.0 ) {
			abox_a = abox_from_tcode( tl[level] - 0.5*dtl[level] );
		} else {
			abox_a = abox_from_tcode( tl[level] + 0.5*dtl[level] );
		}
		dgrowth_a = growth(abox_a);

		x_a = q + dgrowth_a*ampl*sin(kq);
                                                                                                                                                            
		rho = cosmology->OmegaM / ( 1.0 + ak*dgrowth*ampl*cos(kq) );
		v = ddgrowthdt*ampl*sin(kq);
/* 		grav_accel = -6.0*abox_a*( rhogas0*q - x_a ); */
		grav_accel = -6.0*abox[level]*( rhogas0*q - x0 );
		phi = c1*(c2*(0.5*kq*kq + c3*(kq*sin(kq)+cos(kq)-1.0)) - 0.5*(x0*x0)) + phi0; 
#ifdef HYDRO
                rho0 = rhogas0 / ( 1.0 + ak*dgrowth0*ampl*cos(kq) );
		TK = TinitK * 
                    pow(
                        pow(auni_init/abox[level],3.0)
                        *rho/rho0
                        ,constants->gamma-1.0);

                //~T=~P/~rho
                cell_tempK = cell_gas_temperature(cell)*units->temperature;
                cell_tempPK = cell_gas_pressure(cell)/cell_gas_density(cell)*units->temperature*constants->wmu;
                    
		fprintf(output, "%e %e   %e %e %e   %e %e %e   %e %e %e   %e %e %e   %e %e %e   %e %e %e %e   %e %e %d %d\n",
			q, x0, 
//3
			(cell_gas_density(cell)/cell_volume[level]+1.0), rho, fabs(rho-(cell_gas_density(cell)/cell_volume[level]+1.0))/rho,
//6
			cell_momentum(cell,0)/cell_gas_density(cell), v, fabs((v-cell_momentum(cell,0)/cell_gas_density(cell))/v),
//9
#ifdef GRAVITY
                        cell_accel(cell, 0)/dtl[level], grav_accel, fabs(grav_accel-cell_accel(cell, 0)/dtl[level])/grav_accel,
#else
                        -99.00,-99.00,-99.00,
#endif                        
//12                        
                        cell_tempK, TK, fabs(TK-cell_tempK)/TK,
//15                        
#ifdef GRAVITY
			cell_potential(cell), phi, fabs((phi-cell_potential(cell))/cell_potential(cell)),
#else
                        -99.00,-99.00,-99.00,
#endif                        
//18                        
                        cell_gas_internal_energy(cell),
                        cell_gas_density(cell)*TK/units->temperature/(cell_gas_gamma(cell)-1)/constants->wmu,
                        cell_gas_energy(cell)-cell_gas_kinetic_energy(cell),
                        cell_gas_pressure(cell)/(cell_gas_gamma(cell)-1),
//22                        
			cell_gas_energy(cell), cell_gas_kinetic_energy(cell),
                        cell_level(cell), cell );

/* 		fprintf(output, "%e %e %u %e %e %e %e %e %e %e %e %u %e %e %e\n", q, x0, cell_level(cell), */
/* 			cell_gas_density(cell), cell_potential(cell), cell_potential_hydro(cell), */
/* 			cell_gas_density(cell), cell_momentum(cell,0)/cell_gas_density(cell), */
/* 			cell_gas_pressure(cell), cell_momentum(cell,0), ref[cell], cell, phi, rho, v ); */
#else

		fprintf(output, "%e %e %e %e %e %e %e %e %e %e %e\n",
			q, x0,
			(cell_total_mass(cell)+cell_volume[level]), rho, fabs(rho-(cell_total_mass(cell)+cell_volume[level]))/rho,
			cell_accel(cell, 0), grav_accel, fabs(grav_accel-cell_accel(cell, 0))/grav_accel,
			cell_potential(cell), phi, fabs(phi-cell_potential(cell))/cell_potential(cell) );
#endif /* HYDRO */
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

	sprintf(filename, "dumps/zeldovich_[%5.4f]_%03d.dat", abox[min_level], local_proc_id );
	output = fopen(filename, "w");

#ifdef PARTICLES
	sprintf(filename, "dumps/zeldovich_particles_[%5.4f]_%03d.dat", abox[min_level], local_proc_id );
	particles = fopen(filename, "w");
#endif
	dgrowth0 = growth(auni_init);
	dgrowth = growth(abox[min_level]);
        ddgrowthdt = dgrowthdt(abox[min_level]);

	c1 = 6.0/abox[min_level];
	c2 = 1.0 / (ak*ak);
	c3 = dgrowth*ampl*ak;
#ifdef GRAVITY
	pos[1] = pos[2] = (double)num_grid / 2.0 - 0.5*cell_size[max_level];
	pos[0] = 0.0;
	i = cell_find_position(pos);
	if ( i != -1 ) {
		phi = cell_potential(i);
	} else { 
		phi = 0.0;
	}
#endif

	MPI_Allreduce( &phi, &phi0, 1, MPI_DOUBLE, MPI_MAX, mpi.comm.run );

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
            i=cell_find_position(pos);
		if ( i != -1 ) {
			output_cell( i, cell_level(i) );
		}
	}

	fclose(output);

	/* compute rms fluctuations */
	x_rms = x_norm = 0.0;
	v_rms = v_norm = 0.0;

	ddgrowthdt = dgrowthdt( abox_from_tcode( tl[min_level] - 0.5*dtl[min_level] ) );

#ifdef PARTICLES
	for ( i = 0; i < num_particles; i++ ) {
		if ( particle_level[i] != FREE_PARTICLE_LEVEL ) {
			ipart = i;
			x0 = particle_x[i][0];
			q = root_finder( qsolve, 0.0, num_grid, 1e-9, 1e-9 );

			q_i = particle_q_init( particle_id[i] );

			x_z = q_i + dgrowth*ampl*sin(ak*q_i);
			v_z = ddgrowthdt*ampl*sin(ak*q_i);

			x_rms += (x0-x_z)*(x0-x_z);
			x_norm += (x_z-q_i)*(x_z-q_i);

			v_rms += (particle_v[i][0]-v_z)*(particle_v[i][0]-v_z);
			v_norm += v_z*v_z;
                        
			fprintf(particles, "%u %e %e %e %e %e \n",
                                particle_id[i],
				particle_x[i][0],
				particle_v[i][0], 
				particle_pot[i],
                                x_z,v_z ); 
                        
                        
		}
	}
        
	fclose(particles);
	MPI_Reduce( &x_rms, &x_rms_total, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
	MPI_Reduce( &v_rms, &v_rms_total, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
	MPI_Reduce( &x_norm, &x_norm_total, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
	MPI_Reduce( &v_norm, &v_norm_total, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );

	if ( local_proc_id == MASTER_NODE ) {	

		x_rms_total /= x_norm_total;
		v_rms_total /= v_norm_total;

		x_rms_total = sqrt(x_rms_total);
		v_rms_total = sqrt(v_rms_total);

		particles = fopen("dumps/particle_rms.dat", "a");
		fprintf(particles, "%e %e %e %e %e\n", abox[min_level], x_rms_total, v_rms_total, x_norm_total, v_norm_total );
		fclose(particles);
	}
#endif /* PARTICLES */
}

void units_set_art(double OmegaM, double h, double Lbox);

void init_run() {
	int i, j, k;
	int index;
	double a_th;

	int ipart;
	int icell;
	double qi, qj, qk;
	double xcons, vcons;
	double dx, dvx;
	double pw;
	double a_vel;
	double qfact;

	int num_level_cells;
	int *level_cells;

        cosmology_set(OmegaM,OmM0);
        cosmology_set(OmegaB,OmB0);
        cosmology_set(OmegaL,OmL0);
        cosmology_set(h,h0);
        cosmology_set(DeltaDC,dDC);
        box_size = Lbox0;
        
	units_set_art(cosmology->OmegaM,cosmology->h,box_size);
        units_reset();
        build_cell_buffer();
        repair_neighbors();
        
	auni[min_level] = auni_init;
	tl[min_level] = tcode_from_auni( auni_init );
	for ( i = min_level; i <= max_level; i++ ) { tl[i] = tl[min_level]; }
	abox[min_level] = auni_init;

        for(i=min_level+1; i<=max_level; i++)
        {
            tl[i] = tl[min_level];
            auni[i] = auni[min_level];
            abox[i] = abox[min_level];
        }
        
        units_reset();
        units_update(min_level);
        cart_debug("tl[min_level] = %f", tl[min_level] );
        cart_debug("au[min_level] = %f", auni[min_level] );
        cart_debug("ab[min_level] = %f", abox[min_level] );
        cart_debug("DC mode = %f", cosmology->DeltaDC );
        cosmology_set_fixed();


	rhogas0 = cosmology->OmegaB/cosmology->OmegaM;
	cart_debug("rhogas0 = %e", rhogas0 );

	Tinit = TinitK/units->temperature;
        

	ak = 2.0*M_PI / lambda;
	dgrowth = growth(abox[min_level]);
	ddgrowthdt = dgrowthdt(abox[min_level]);
	ampl = 1.0 / ( growth(a_cross) * ak );
	cart_debug("Tinit,TinitK = %e %e", Tinit,TinitK );

#ifdef HYDRO
	for ( i = min_level; i <= max_level; i++ ) {
		cart_debug("generating initial conditions on level %u", i );
	
		select_level( i, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( j = 0; j < num_level_cells; j++ ) {
//                    cart_debug("%d %d",level_cells[j],num_cells);
			initial_conditions( level_cells[j], i );
		}
		cart_free( level_cells );
	}

	for ( i = min_level; i <= max_level; i++ ) {
		update_buffer_level( i, all_hydro_vars, num_hydro_vars );
	}
#endif /* HYDRO */

        cart_debug("choose timestep and set velocity on the half step");
        dtl[min_level] = 0.0;
        set_timestepping_scheme();
        dtl[min_level]=.125;
        cart_debug("=======================%e",dtl[min_level]);
        
        dtl_old[min_level] = dtl[min_level];
        tl_old[min_level] = tl[min_level]-dtl[min_level];
        abox_old[min_level] = abox_from_tcode(tl_old[min_level]);
        dtl_old[min_level] = dtl[min_level];
        
	for ( i = min_level+1; i <= max_level; i++ ) {
            tl_old[i] = tl[i]-dtl[i];
            abox_old[i] = abox_from_tcode(tl_old[i]);
            dtl_old[i] = dtl[i];
        }
                
#ifdef GRAVITY
#ifdef HYDRO
	for ( i = min_level; i <= max_level; i++ ) {
		cart_debug("generating gravity on level %u", i );
	
//                cart_assert(dtl[i]==0);
		select_level( i, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
		for ( j = 0; j < num_level_cells; j++ ) {
			initial_gravity( level_cells[j], i );
		}
		cart_free( level_cells );
	}

	for ( i = min_level; i <= max_level; i++ ) {
		update_buffer_level( i, all_hydro_vars, num_hydro_vars );
	}
#endif /* GRAVITY */
#endif /* HYDRO */
        
        
#ifdef PARTICLES
	num_row = num_grid;
	qfact = (double)num_grid / (double)num_row;
	pw = (1.0-rhogas0)*qfact*qfact*qfact;

	cart_debug("particle weight = %e", pw );

	xcons = dgrowth*ampl;
	a_vel = abox_from_tcode( tl[min_level] - 0.5*dtl[min_level] );
	vcons = ampl * dgrowthdt( a_vel);

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
                                        particle_mass[ipart] = pw;

					cart_assert( qi == particle_q_init( particle_id[ipart] ) );

					particle_t[ipart] = tl[min_level];
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


/* 	assign_density( min_level, min_level ); */ //for refinement
/* 	modify( min_level, 0 ); */
 	assign_density( min_level, min_level );  //for refinement
 	modify( min_level, 0 ); 


	if ( local_proc_id == MASTER_NODE ) {
		particles = fopen("dumps/particle_rms.dat", "w");
		fclose(particles);
	}

	

#endif

        cart_debug("snl10 %e ",cell_accel(0,0));

	check_map();
	cart_debug("done with initialization");
}
