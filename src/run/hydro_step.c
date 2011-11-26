#include "config.h"
#ifdef HYDRO 

#include <math.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "cooling.h"
#include "cosmology.h"
#include "hydro.h"
#include "iterators.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "hydro_step.h"
#include "rt_step.h"
#include "step.h"


#ifndef HYDRO_CHUNK_SIZE
#define HYDRO_CHUNK_SIZE        65536
#endif /* HYDRO_CHUNK_SIZE */

extern int pressure_floor_min_level; /* NG: that used to be MinL_Jeans define */
extern float pressure_floor_factor;
float pressure_floor;

extern int pressureless_fluid_eos;           /* NG: that used to be PRESSURELESS_FLUID define */
extern int apply_lapidus_viscosity;          /* NG: that used to be LAPIDUS define */
extern int smooth_density_gradients;         /* NG: that used to be DENSGRADSMOOTH define */

extern float gas_density_floor;
extern float gas_temperature_floor;          /* NG: that used to be T_min define */
/*
//  Radiation pressure fudge factor (Hopkins' parameter eta_p)
*/
extern float radiation_pressure_factor;


#ifdef BLASTWAVE_FEEDBACK
extern double blastwave_time_floor; 
extern double blastwave_time_cut;
#endif /* BLASTWAVE_FEEDBACK */

#ifdef STATIC_MESH_DATA
#define backup_hvar(c,v)	(backup_hvars[c][v])
#else  /* STATIC_MESH_DATA */
#define backup_hvar(c,v)	(backup_hvars[(c)*(num_hydro_vars-2)+(v)])
#endif /* STATIC_MESH_DATA */

float CELL_VAR_ARRAY(backup_hvars,num_hydro_vars-2) STATIC_INIT;
float CELL_ARRAY(ref) STATIC_INIT;
int CELL_ARRAY(backup_dirty) STATIC_INIT;

#ifdef GRAVITY_IN_RIEMANN
void fluxh( double dtx, double dtx2, double v[num_hydro_vars-1][4], double g[2], double c[2], double f[num_hydro_vars-1] );
#else
void fluxh( double dtx, double dtx2, double v[num_hydro_vars-1][4], double c[2], double f[num_hydro_vars-1] );
#endif

void lapidus( double dtx2, int L1, int R1, int sweep_direction, int mj3, int mj4, int mj5, double v[num_hydro_vars-1][4], double f[num_hydro_vars-1] );

float cell_radiation_pressure(int cell);

DEFINE_LEVEL_ARRAY(int,level_sweep_dir);
const int sweep_dir[2][nDim] = { { 1, 3, 5 }, { 5, 3, 1 } };

int sweep_direction;
int sweep_dimension;

int j3,j4,j5;
int mj3, mj4, mj5;
double dtx;
double dtx2;
double dxi;
double dxi2;

const int momentum_permute[2*nDim][nDim] = {  { 0, 1, 2 }, { 0, 1, 2 },
	{ 1, 0, 2 }, { 1, 0, 2 },
	{ 2, 1, 0 }, { 2, 1, 0 } };


void apply_hydro_fluxes( int icell, double factor, double dxi_factor, double f[ /* num_hydro_vars-1 */ ] );
void hydro_sweep_1d( int level );
#ifdef GRAVITY
void hydro_apply_gravity( int level );
#endif /* GRAVITY */
void compute_hydro_fluxes( int cell_list[4], double f[ /* num_hydro_vars-1 */ ] );
void hydro_eos(int level);
void hydro_magic(int level);
void hydro_advance_internalenergy(int level);
void hydro_split_update( int level );


void hydro_step_init() {
#ifndef STATIC_MESH_DATA
	if(backup_hvars == NULL)
	  {
	    backup_hvars = cart_alloc(float,(size_t)num_cells*(num_hydro_vars-2));
	    ref = cart_alloc(float,num_cells);
	  }
#endif  /* STATIC_MESH_DATA */
}


void hydro_step( int level ) {
	int dir;

	if ( pressure_floor_min_level >= 0 && level >= pressure_floor_min_level ) {
		/* artificial pressure floor */
		pressure_floor = pressure_floor_factor * constants->G*pow(units->density*units->length,2.0)/units->energy_density * cell_size[level] * cell_size[level];
	} else {
		pressure_floor = 0.0;
	}

	dtx = dtl[level] * cell_size_inverse[level];
	dxi = cell_size_inverse[level];
	dxi2 = 0.5*cell_size_inverse[level];
	dtx2 = 0.5*dtx;

	start_time( HYDRO_TIMER );

	for ( dir = 0; dir < nDim; dir++ ) {
		start_time( WORK_TIMER );

		sweep_direction = sweep_dir[level_sweep_dir[level]][dir];
		sweep_dimension = (sweep_direction-1)/2;

		j3 = momentum_permute[sweep_direction][0];
		j4 = momentum_permute[sweep_direction][1];
		j5 = momentum_permute[sweep_direction][2];

		mj3 = j3+2;
		mj4 = j4+2;
		mj5 = j5+2;

		/* compute fluxes across cell interfaces long sweep_dimension */
		hydro_sweep_1d( level );

#if defined(GRAVITY) && (!defined(GRAVITY_IN_RIEMANN))
		hydro_apply_gravity( level );
#endif /* GRAVITY && !GRAVITY_IN_RIEMANN */

		end_time( WORK_TIMER );

		if ( dir == nDim - 1 ) {
			hydro_copy_vars( level, HYDRO_RESTORE_ALL );
		} else {

			hydro_copy_vars( level, HYDRO_RESTORE_CLEAN );
					
			hydro_magic( level );
			hydro_eos( level );
			
			start_time( HYDRO_UPDATE_TIMER );
			update_buffer_level( level, all_hydro_vars, num_hydro_vars );
			end_time( HYDRO_UPDATE_TIMER );
		}
	}

	hydro_advance_internalenergy( level );
	hydro_magic( level );
	hydro_eos( level );

	hydro_split_update( level );

	start_time( HYDRO_UPDATE_TIMER );
	update_buffer_level( level, all_hydro_vars, num_hydro_vars );
	end_time( HYDRO_UPDATE_TIMER );
	
	/* update sweep direction */
	level_sweep_dir[level] = (level_sweep_dir[level]+1)%2;

	end_time( HYDRO_TIMER );
}

void hydro_sweep_1d( int level ) {
	int i, j;
	int icell;
	int L2, L1, R1, R2;
	int count;
	int num_level_cells;
	int *level_cells;

	int cell_list[HYDRO_CHUNK_SIZE][4];
	double f[HYDRO_CHUNK_SIZE][num_hydro_vars-1];

	count = 0;

	if ( level == min_level ) {
		select_level( min_level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
	} else {
		select_level( level, CELL_TYPE_ANY_LEAF, &num_level_cells, &level_cells );
	}

	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		/* calculate neighbors */
		L1 = cell_neighbor( icell, reverse_direction[sweep_direction]);
		R1 = cell_neighbor( icell, sweep_direction );

		if ( ( cell_is_local(icell) && cell_is_leaf(R1) ) ||
				(!cell_is_local(icell) && R1 != -1 && cell_is_local(R1) && cell_is_leaf(R1) ) ) {

			R2 = cell_neighbor( R1, sweep_direction );

			cell_list[count][0] = L1;
			cell_list[count][1] = icell;
			cell_list[count][2] = R1;
			cell_list[count][3] = R2;
			count++;
		}

		if ( ( level == min_level && !cell_is_local(L1) && cell_is_leaf(L1) ) || 
				( L1 != -1 && cell_level(L1) == level - 1 && 
				( cell_is_local(icell) || cell_is_local(L1) ) ) ) {

			L2 = cell_neighbor( L1, reverse_direction[sweep_direction] );

			cell_list[count][0] = L2;
			cell_list[count][1] = L1;
			cell_list[count][2] = icell;
			cell_list[count][3] = R1;
			count++;
		}

		/* must be HYDRO_CHUNK_SIZE-1 since each cell can add 2 interfaces */
		if ( i == num_level_cells-1 || count >= HYDRO_CHUNK_SIZE-1 ) {
#pragma omp parallel for default(none), private(j), shared(cell_list,f,count,dtx,dtx2,sweep_direction,j3,j4,j5)
			for ( j = 0; j < count; j++ ) {
				compute_hydro_fluxes( cell_list[j], f[j] );
			}

 			/* apply fluxes to left cells */
#pragma omp parallel for default(none), private(j,icell), shared(cell_list,count,level,f,dxi,dxi2)
			for ( j = 0; j < count; j++ ) {
				icell = cell_list[j][1];

				if ( cell_is_local(icell) && cell_level(icell) == level ) {
					apply_hydro_fluxes( icell, -1.0, dxi, f[j] );
				}
			}

#pragma omp parallel for default(none), private(j,icell), shared(count,cell_list,level,f,dxi,dxi2)
			for ( j = 0; j < count; j++ ) {
				icell = cell_list[j][2];

				if ( cell_is_local(icell) && cell_level(icell) == level ) {
					apply_hydro_fluxes( icell, 1.0, dxi, f[j] );
				}
			}

			/* Apply fluxes to higher level cells on both left and right interfaces 
			 * MUST execute serially since lower level cells will appear multiple times! */
			for ( j = 0; j < count; j++ ) {
				icell = cell_list[j][1];

				if ( cell_is_local(icell) && cell_level(icell) < level ) {
					apply_hydro_fluxes( icell, -0.125, dxi2, f[j] );
					backup_dirty[icell] = 1;
				}

				icell = cell_list[j][2];
				if ( cell_is_local(icell) && cell_level(icell) < level ) {
					apply_hydro_fluxes( icell, 0.125, dxi2, f[j] );
					backup_dirty[icell] = 1;
				}
			}

			count = 0;
		}
	}

	cart_free( level_cells );
}

#if defined(GRAVITY) && (!defined(GRAVITY_IN_RIEMANN))
void hydro_apply_gravity( int level ) {
	int i, icell;
	int num_level_cells;
	int *level_cells;
	double gravadd;

	/* now we need to apply a gravity correction */
	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(icell,gravadd), shared(num_level_cells,level_cells,cell_child_oct,backup_hvars,cell_vars,sweep_dimension,mj3)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		gravadd = backup_hvar(icell,0) * cell_accel(icell,sweep_dimension);
		backup_hvar(icell,1) += cell_accel(icell,sweep_dimension) *
		  ( backup_hvar(icell,mj3) + 0.5 * gravadd );
		backup_hvar(icell,mj3) += gravadd;
	}

	cart_free( level_cells );
}
#endif /* GRAVITY && !GRAVITY_IN_RIEMANN */

void hydro_magic_one_cell( int icell ) {
	int j;
	float average_density;
	int neighbors[num_neighbors];
	static int failed = 0;

	double kinetic_energy;
	double thermal_energy;

	/* do density floor stuff */
	if ( cell_gas_density(icell) < gas_density_floor ) {
		average_density = 0.0;
		cell_all_neighbors( icell, neighbors );

		for ( j = 0; j < num_neighbors; j++ ) {
			average_density += cell_gas_density( neighbors[j] );
		}

		if(!failed) {
			cart_debug("---------------------------------------------------------");
			cart_debug("HIT DENSITY FLOOR:");
			cart_debug("old density = %e g/cc", cell_gas_density(icell)*units->density/constants->gpercc );
			cart_debug("new density = %e g/cc", max( average_density/(float)num_neighbors, gas_density_floor ) );
			cart_debug("T  = %e K", cell_gas_temperature(icell)*units->temperature/constants->K );
			cart_debug("P  = %e ergs cm^-3", cell_gas_pressure(icell)*units->energy_density/constants->barye );
			cart_debug("v  = %e %e %e cm/s",
					cell_momentum(icell,0)/cell_gas_density(icell)*units->velocity/constants->cms,
					cell_momentum(icell,1)/cell_gas_density(icell)*units->velocity/constants->cms,
					cell_momentum(icell,2)/cell_gas_density(icell)*units->velocity/constants->cms );
			cart_debug("---------------------------------------------------------");
		}	

		failed = 1;

		cell_gas_density(icell) = max( average_density/(float)num_neighbors, gas_density_floor );
	}

	kinetic_energy = cell_gas_kinetic_energy(icell);
	thermal_energy = units->Emin * cell_gas_density(icell);

	cell_gas_internal_energy(icell) = max( cell_gas_internal_energy(icell), thermal_energy );
	cell_gas_energy(icell) = max( cell_gas_energy(icell), thermal_energy+kinetic_energy );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	cell_electron_internal_energy(icell) = max( cell_electron_internal_energy(icell), thermal_energy*constants->wmu/constants->wmu_e );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = 0; j < num_chem_species; j++ ) {
		/* 
		   1e-15 may be too large a number for ionic species;
		   at least let's scale them with density and make 1e-20;
Gnedin: 1e-20 is not small enough for chemistry, making it 1e-30
		 */
		cell_advected_variable(icell,j) = max( 1e-30*cell_gas_density(icell), cell_advected_variable(icell,j) );
		/* Doug had it like that:
		   cell_advected_variable(icell,j) = max( 1e-15, cell_advected_variable(icell,j) );
		 */
	}
}

void hydro_magic( int level ) {
    int i;
    int icell;
    int num_level_cells;
    int *level_cells;

    start_time( WORK_TIMER );

    select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,cell_child_oct)
    for ( i = 0; i < num_level_cells; i++ ) {
		hydro_magic_one_cell(level_cells[i]);
	}
    cart_free( level_cells );

    end_time( WORK_TIMER );
} 

void hydro_eos( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;
	double kinetic_energy;

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,kinetic_energy), shared(num_level_cells,level_cells,cell_child_oct,cell_vars,constants,pressureless_fluid_eos)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		kinetic_energy = cell_gas_kinetic_energy(icell);

		if(pressureless_fluid_eos)
		  {
		    cell_gas_pressure(icell) = 1e-20;
		    cell_gas_internal_energy(icell) = cell_gas_pressure(icell) / (constants->gamma-1.0);
		    cell_gas_energy(icell) = cell_gas_internal_energy(icell) + kinetic_energy;
		  }
		else
		  {
		    cell_gas_internal_energy(icell) = max( cell_gas_internal_energy(icell), 0.0 );
		    cell_gas_energy(icell) = max( kinetic_energy, cell_gas_energy(icell) );
		    cell_gas_pressure(icell) = max( (cell_gas_gamma(icell)-1.0)*cell_gas_internal_energy(icell), 0.0 );

		    /*
		    // NG: this is done this way to ensure that 
		    // cell_radiation_pressure(...) is not called 
		    // at all when RP is off, since it may not be
		    // implemented yet.
		    */
		    if(radiation_pressure_factor > 0.0)
		      {
			cell_gas_pressure(icell) += radiation_pressure_factor*cell_radiation_pressure(icell);
		      }

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		    cell_electron_internal_energy(icell) = min( cell_electron_internal_energy(icell), cell_gas_internal_energy(icell)*constants->wmu/constants->wmu_e );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
		  }
	}

	cart_free( level_cells );

	end_time( WORK_TIMER );
}


#if defined(COOLING) && !defined(RADIATIVE_TRANSFER)

#ifndef OLDSTYLE_COOLING_EXPLICIT_SOLVER
void qss_getcooling ( double t, double *y, void *params, double *w, double *a) {
	double nHlog = ((double *)params)[0];
	double Tfac_cell = ((double *)params)[1];
	double Zlog = ((double *)params)[2];
	double rhog2 = ((double *)params)[3];
	double unit_cl = ((double *)params)[6];
	double t0 = ((double *)params)[7];
	double Hdum = ((double *)params)[8];
	double f_curr = 1 + Hdum*(t-t0);
	double etmp = rhog2*f_curr*unit_cl;

	cooling_t coolrate = cooling_rate(nHlog,y[0]*Tfac_cell/(f_curr*f_curr),Zlog);
	
	a[0] = etmp*coolrate.Cooling/ y[0];
	w[0] = etmp*coolrate.Heating;
}

void adjust_internalenergy( double t, double *y, void *params ) {
  /* RL: put temperature/internal energy floor in here??? */
	double Emin_cell = ((double *)params)[5];
	if (y[0] < Emin_cell) y[0] = Emin_cell;
}

void hydro_apply_cooling(int level, int num_level_cells, int *level_cells) {
	int i;
	int icell;
	double t_begin, t_end;
	double Zlog, Hdum;
	double Tfac, Tfac_cell, Emin_cell, blastwave_time;
	double rhog2, nHlog;
	double e_curr;
	double unit_cl = units->time*pow(constants->XH*units->number_density,2.0)/units->energy_density;
	double err[1] = { 1e-2 };
	double params[9];
	qss_system *sys;

	t_begin	= tl[level];
	t_end = tl[level] + dtl[level];

#ifdef COSMOLOGY
	Hdum = (abox_from_tcode(t_end)/abox[level] - 1) / dtl[level];
#else
	Hdum = 0.0;
#endif

	/* Note: removed 10^4 K normalization since it wasn't used - DHR */
	Tfac = units->temperature*constants->wmu*( constants->gamma-1 )/constants->K;

#ifdef BLASTWAVE_FEEDBACK
#pragma omp parallel default(none), shared(num_level_cells,level_cells,level,t_begin,Tfac,units,t_end,cell_child_oct,err,constants,cell_vars,Hdum,unit_cl,blastwave_time_cut,blastwave_time_floor), private(i,icell,rhog2,nHlog,Zlog,Tfac_cell,e_curr,Emin_cell,params,sys,blastwave_time)
#else
#pragma omp parallel default(none), shared(num_level_cells,level_cells,t_begin,Tfac,units,t_end,cell_child_oct,err,constants,cell_vars,Hdum,unit_cl), private(i,icell,rhog2,nHlog,Zlog,Tfac_cell,e_curr,Emin_cell,params,sys)
#endif /* BLASTWAVE_FEEDBACK*/
	{
	  sys = qss_alloc( 1, &qss_getcooling, &adjust_internalenergy );

#pragma omp for
	  for ( i = 0; i < num_level_cells; i++ ) {
		  icell = level_cells[i];
		  if ( cell_is_leaf(icell) ) {
#ifdef BLASTWAVE_FEEDBACK
		    blastwave_time = cell_blastwave_time(icell) / cell_gas_density(icell);
		    if(blastwave_time <= blastwave_time_cut){ 
#endif /* BLASTWAVE_FEEDBACK */

		    cell_gas_gamma(icell) = constants->gamma;
		    rhog2 = cell_gas_density(icell)*cell_gas_density(icell);
		    /* take code density -> log10(n_H [cm^-3]) */
		    nHlog = log10(constants->XH*units->number_density*cell_gas_density(icell)/constants->cc);
#ifdef ENRICH
		    Zlog = log10(max(1.0e-10,cell_gas_metal_density(icell)/(constants->Zsun*cell_gas_density(icell))));
#else
		    Zlog = -10.0;
#endif /* ENRICH */
			  
		    Tfac_cell = Tfac/cell_gas_density(icell);
		    Emin_cell = units->Emin*cell_gas_density(icell);
		    
		    e_curr = cell_gas_internal_energy(icell);

		    params[0] = nHlog; //
		    params[1] = Tfac_cell; // to get the cooling rate...
		    params[2] = Zlog; //
		    params[3] = rhog2;
		    params[4] = e_curr;
		    params[5] = Emin_cell;
		    params[6] = unit_cl;
		    params[7] = t_begin;
		    params[8] = Hdum;
	
		    qss_solve( sys, t_begin, t_end, &e_curr, err, &params );

		    cell_gas_internal_energy(icell) = max(Emin_cell,e_curr);
		    cell_gas_energy(icell) = cell_gas_kinetic_energy(icell) + cell_gas_internal_energy(icell);
#ifdef BLASTWAVE_FEEDBACK
		  }else { 
		    blastwave_time -= dtl[level]*units->time/constants->yr; 
		    if(blastwave_time < blastwave_time_cut ){
		      blastwave_time = blastwave_time_floor;
		    }
		    cell_blastwave_time(icell) = cell_gas_density(icell) * blastwave_time;
		  }
#endif /* BLASTWAVE_FEEDBACK */
		  }
		}

		qss_free(sys);
	} /* END omp parallel block */
}

#else /* OLDSTYLE_COOLING_EXPLICIT_SOLVER */

void hydro_cool_one_cell(int icell, double t_begin, double t_end, double Hdum, double Zlog, double nHlog, double rhog2, double Tfact_cell, double Emin_cell, double unit_cl) {
	int continue_cooling;
	double t_curr, f_curr;
	double dE;
	double dt_e, ei1, T_gas;

#define dstep	(0.01)

	continue_cooling = 1;
	t_curr = t_begin;

	/* integrate cooling using smaller timestep */
	while ( ( t_curr < t_end ) && continue_cooling ) {
		f_curr = 1 + Hdum*(t_curr-t_begin);
		T_gas = Tfact_cell * cell_gas_internal_energy(icell) / ( f_curr * f_curr );

		/* compute new timestep */
		dE = unit_cl*cooling_rate( nHlog, T_gas, Zlog );
		dE *= -rhog2 * f_curr;
		dt_e = min( dstep * fabs( cell_gas_internal_energy(icell) / dE ), t_end - t_curr );

		ei1 = max( cell_gas_internal_energy(icell) + 0.5 * dE * dt_e, Emin_cell );
		T_gas = Tfact_cell * ei1 / ( f_curr * f_curr );

		dE = unit_cl*cooling_rate( nHlog, T_gas, Zlog );
		dE *= -rhog2 * f_curr * dt_e;
		/* adjust cell energies */
		cell_gas_internal_energy(icell) += dE;
		cell_gas_energy(icell) += dE;

		/* stop if we hit energy minimum */
		if ( cell_gas_internal_energy(icell) < Emin_cell ) {
			continue_cooling = 0;
		}

		cell_gas_internal_energy(icell) = max( Emin_cell, cell_gas_internal_energy(icell) );
		cell_gas_energy(icell) = max( Emin_cell, cell_gas_energy(icell) );

		/* advance timestep */
		t_curr += dt_e;
	}
}

void hydro_apply_cooling(int level, int num_level_cells, int *level_cells) {
	int i;
	int icell;
	double t_begin, t_end;
	double Zlog, Hdum;
	double Tfac, Tfac_cell, Emin_cell;
	double rhog2, nHlog;
	double unit_cl = units->time*pow(constants->XH*units->number_density,2.0)/units->energy_density;

	t_begin	= tl[level];
	t_end = tl[level] + dtl[level];

#ifdef COSMOLOGY
	Hdum = (abox_from_tcode(t_end)/abox[level] - 1) / dtl[level];
#else
	Hdum = 0.0;
#endif

	/* Note: removed 10^4 K term since it isn't used - DHR */
	Tfac = units->temperature*constants->wmu*(constants->gamma-1)/constants->K;

#ifdef BLASTWAVE_FEEDBACK
#pragma omp parallel for default(none), private(icell,i,rhog2,nHlog,Zlog,Tfac_cell,Emin_cell,blastwave_time), shared(num_level_cells,level_cells,t_begin,t_end,Tfac,units,constants,cell_child_oct,cell_vars,Hdum,unit_cl,blastwave_time_cut,blastwave_time_floor)
#else
#pragma omp parallel for default(none), private(icell,i,rhog2,nHlog,Zlog,Tfac_cell,Emin_cell), shared(num_level_cells,level_cells,t_begin,t_end,Tfac,units,constants,cell_child_oct,cell_vars,Hdum,unit_cl)
#endif /* BLASTWAVE_FEEDBACK*/ 
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];
		if ( cell_is_leaf(icell) ) {
#ifdef BLASTWAVE_FEEDBACK
		  blastwave_time = cell_blastwave_time(icell) / cell_gas_density(icell);
		  if(blastwave_time <= blastwave_time_cut){
#endif /* BLASTWAVE_FEEDBACK */
		
			cell_gas_gamma(icell) = constants->gamma;
			rhog2 = cell_gas_density(icell)*cell_gas_density(icell);
			/* take code density -> log10(n_H [cm^-3]) */
			nHlog = log10(constants->XH*units->number_density*cell_gas_density(icell)/constants->cc);

#ifdef ENRICH
			Zlog = log10(max(1.0e-10,cell_gas_metal_density(icell)/(constants->Zsun*cell_gas_density(icell))));
#else
			Zlog = -10.0;
#endif /* ENRICH */

			Tfac_cell = Tfac/cell_gas_density(icell);
			Emin_cell = units->Emin*cell_gas_density(icell);

			hydro_cool_one_cell(icell,t_begin,t_end,Hdum,Zlog,nHlog,rhog2,Tfac_cell,Emin_cell,unit_cl);
			
#ifdef BLASTWAVE_FEEDBACK
		  }else { 
		    blastwave_time -= dtl[level]*units->time/constants->yr; 
		    if(blastwave_time < blastwave_time_cut ){
		      blastwave_time = blastwave_time_floor;
		    }
		    cell_blastwave_time(icell) = cell_gas_density(icell) * blastwave_time;
		  }
#endif /* BLASTWAVE_FEEDBACK */
		}
	}
}
#endif /* OLDSTYLE_COOLING_EXPLICIT_SOLVER */

#endif /* COOLING && !RADIATIVE_TRANSFER */

#ifdef ELECTRON_ION_NONEQUILIBRIUM
void heating_rates ( double t, double *y, void *params, double *w, double *a) {
	double dEfact = ((double *)params)[0];
	double e_equil = ((double *)params)[1];
	double e_init = ((double *)params)[2];
	double t0 = ((double *)params)[3];
	double Hdum = ((double *)params)[4];
	double f_curr = 1 + Hdum*(t-t0);
	double e_curr = max( e_init, y[0] );

	a[0] = dEfact*f_curr*f_curr*pow(e_curr,-1.5);
	w[0] = a[0]*e_equil;
}

void adjust_temperatures( double t, double *y, void *params ) {
	double e_equil = ((double *)params)[1];
	double e_init = ((double *)params)[2];

	if ( y[0] < e_init ) {
		y[0] = e_init;
	} else if ( y[0] > e_equil ) {
		y[0] = e_equil;
	}
}

void hydro_apply_electron_heating(int level, int num_level_cells, int *level_cells) {
	int i;
	int icell;
	double t_begin, t_end, Hdum;
	double e_equil, e_curr;
	double nfact, Tefact, dEfact;
	double n_5, Te7, dEfact_cell;
	double logcoulomb;

	double err[1] = { 1e-2 };
	double params[5];
	qss_system *sys;

	t_begin	= tl[level];
	t_end = tl[level] + dtl[level];

#ifdef COSMOLOGY
	Hdum = (abox_from_tcode(t_end)/abox[level] - 1) / dtl[level];
#else
	Hdum = 0.0;
#endif

	nfact = 1.0e5*units->number_density*(1.0/constants->wmu - 1.0/constants->wmu_e)/constants->cc;
	Tefact = units->temperature*constants->wmu_e*(constants->gamma-1)/1.0e7/constants->K;
	dEfact = pow(Tefact,-1.5)*units->time/(constants->yr*6.3e8)/40.0/(1.0-constants->wmu/constants->wmu_e)/constants->erg; 

#pragma omp parallel default(none), shared(nfact,Tefact,dEfact,Hdum,cell_vars,t_begin,t_end,num_level_cells,level_cells,cell_child_oct,err,constants), private(i,icell,e_equil,e_curr,Te7,n_5,logcoulomb,dEfact_cell,params,sys)
	{
		sys = qss_alloc( 1, &heating_rates, &adjust_temperatures );

#pragma omp for
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( cell_is_leaf(icell) ) {
				e_equil = cell_gas_internal_energy(icell)*constants->wmu/constants->wmu_e;
				e_curr = cell_electron_internal_energy(icell);
				Te7 = Tefact * cell_electron_internal_energy(icell) / cell_gas_density(icell); /* a^2 Te/10^7 K */
		
				n_5 = nfact * cell_gas_density(icell); 
				logcoulomb = max( 30.0, 37.8 + log(Te7) - 0.5*log(n_5) );
				dEfact_cell = dEfact*n_5*logcoulomb*pow( cell_gas_density(icell), 1.5);
		
				params[0] = dEfact_cell;
				params[1] = e_equil;
				params[2] = e_curr;
				params[3] = t_begin;
				params[4] = Hdum;
	
				qss_solve( sys, t_begin, t_end, &e_curr, err, &params );
	
				cell_electron_internal_energy(icell) = min( e_curr, e_equil );
			}
		}

		qss_free(sys);
	} /* END omp parallel block */
}
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

void hydro_advance_internalenergy( int level ) {
	int i;
	int icell;
	int num_level_cells;
	int *level_cells;
	double kinetic_energy;
	double energy;
	double gamma1, div, div_dt;

	start_time( WORK_TIMER );

	div_dt = dtl[level] / 3.0;

	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(icell,i,kinetic_energy,energy,gamma1,div), shared(num_level_cells,level_cells,cell_child_oct,cell_vars,ref,div_dt)
	for ( i = 0; i < num_level_cells; i++ ) {
		icell = level_cells[i];

		/* P dV term */
		gamma1 = cell_gas_gamma(icell) - 1.0;
		div = 1.0 + gamma1 * ref[icell] * div_dt;
		cell_gas_internal_energy(icell) = max( 1.0e-30, cell_gas_internal_energy(icell)*div*div*div);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
		cell_electron_internal_energy(icell) = max( 1.0e-30, cell_electron_internal_energy(icell)*div*div*div );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

		/* synchronize internal and total energy */
		kinetic_energy = cell_gas_kinetic_energy(icell);
		energy = cell_gas_energy(icell);

		/* we trust energy over internal energy since it's computed using
		 * the riemann solver rather than just advection equation, so if
		 * internal energy is sufficiently large then compute it from 
		 * e = E - rho * v**2 /2 */
		if ( ( energy - kinetic_energy) / energy > 1e-3 ) {
			cell_gas_internal_energy(icell) = energy - kinetic_energy;
		}
	}

#ifdef COOLING
#ifdef RADIATIVE_TRANSFER
	start_time( RT_COOLING_TIMER );
	rtApplyCooling(level,num_level_cells,level_cells);
	end_time( RT_COOLING_TIMER );
#else
	start_time( COOLING_TIMER );
	hydro_apply_cooling(level,num_level_cells,level_cells);
	end_time( COOLING_TIMER );
#endif /* RADIATIVE_TRANSFER */
#else
#ifdef ELECTRON_ION_NONEQUILIBRIUM
	start_time(COOLING_TIMER);
	hydro_apply_electron_heating(level,num_level_cells,level_cells); 
	end_time(COOLING_TIMER);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#endif /* COOLING */

	cart_free( level_cells );

	end_time( WORK_TIMER );
}

void hydro_split_update( int level ) {
	int i, j, k;
	int icell;
	int num_level_cells;
	int *level_cells;
	int children[num_children];
	double new_var;
	const double factor = ((double)(1.0/(1<<nDim)));

	if ( level < max_level ) {
		start_time( WORK_TIMER );

		select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_REFINED, &num_level_cells, &level_cells );
#pragma omp parallel for default(none), private(i,icell,j,k,children,new_var), shared(num_level_cells,level_cells,cell_child_oct,cell_vars)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			/* average over children */
			cell_all_children( icell, children );

			for ( j = 0; j < num_hydro_vars; j++ ) {
				new_var = 0.0;
				for ( k = 0; k < num_children; k++ ) {
					new_var += cell_hydro_variable(children[k], j);
				}

				cell_hydro_variable(icell, j) = new_var*factor;	
			}
		}
		cart_free( level_cells );

		end_time( WORK_TIMER );
	}
}

void apply_hydro_fluxes( int icell, double factor, double dxi_factor, double f[num_hydro_vars-1] ) {
	int j;

        backup_hvar(icell,0) += factor*f[0];
        backup_hvar(icell,1) += factor*f[4];
        backup_hvar(icell,mj3) += factor*f[1];
        backup_hvar(icell,mj4) += factor*f[2];
        backup_hvar(icell,mj5) += factor*f[3];
        backup_hvar(icell,5) += factor*f[5];
	ref[icell] += factor*f[6]*dxi_factor;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
        backup_hvar(icell,6) += factor*f[7];
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = num_hydro_vars-num_chem_species-2; j < num_hydro_vars-2; j++ ) {
                backup_hvar(icell,j) += factor*f[j+1];
	}
}

void compute_hydro_fluxes( int cell_list[4], double f[num_hydro_vars-1] ) {
	int j;
	double v[num_hydro_vars-1][4];
	double c[2];

#ifdef GRAVITY_IN_RIEMANN
	double g[2];
#endif

	int L2 = cell_list[0];
	int L1 = cell_list[1];
	int R1 = cell_list[2];
	int R2 = cell_list[3];

	cart_assert( cell_is_leaf(L1) && cell_is_leaf(R1) );

	/* L2 */
	v[0][0] = cell_gas_density(L2);
	v[1][0] = max( pressure_floor * v[0][0]*v[0][0], cell_gas_pressure(L2) );
	v[2][0] = cell_momentum(L2,j3)/cell_gas_density(L2);
	v[3][0] = cell_momentum(L2,j4)/cell_gas_density(L2);
	v[4][0] = cell_momentum(L2,j5)/cell_gas_density(L2);
	v[5][0] = cell_gas_gamma(L2);
	v[6][0] = constants->gamma;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	v[7][0] = cell_electron_internal_energy(L2);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = 0; j < num_chem_species; j++ ) {
		v[num_hydro_vars-num_chem_species-1+j][0] = cell_advected_variable(L2,j)/cell_gas_density(L2);
	}
	
	/* L1 vars */
	v[0][1] = cell_gas_density(L1);
	v[1][1] = max( pressure_floor * v[0][1]*v[0][1], cell_gas_pressure(L1) );
	v[2][1] = cell_momentum(L1,j3)/cell_gas_density(L1);
	v[3][1] = cell_momentum(L1,j4)/cell_gas_density(L1);
	v[4][1] = cell_momentum(L1,j5)/cell_gas_density(L1);
	v[5][1] = cell_gas_gamma(L1);
	v[6][1] = constants->gamma;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	v[7][1] = cell_electron_internal_energy(L1);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = 0; j < num_chem_species; j++ ) {
		v[num_hydro_vars-num_chem_species-1+j][1] = cell_advected_variable(L1,j)/cell_gas_density(L1);
	}
	
	/* R1 vars */
	v[0][2] = cell_gas_density(R1);
	v[1][2] = max( pressure_floor * v[0][2]*v[0][2], cell_gas_pressure(R1) );
	v[2][2] = cell_momentum(R1,j3)/cell_gas_density(R1);
	v[3][2] = cell_momentum(R1,j4)/cell_gas_density(R1);
	v[4][2] = cell_momentum(R1,j5)/cell_gas_density(R1);
	v[5][2] = cell_gas_gamma(R1);
	v[6][2] = constants->gamma;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	v[7][2] = cell_electron_internal_energy(R1);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = 0; j < num_chem_species; j++ ) {
		v[num_hydro_vars-num_chem_species-1+j][2] = cell_advected_variable(R1,j)/cell_gas_density(R1);
	}
	
	/* R2 vars */
	v[0][3] = cell_gas_density(R2);
	v[1][3] = max( pressure_floor * v[0][3]*v[0][3], cell_gas_pressure(R2) );
	v[2][3] = cell_momentum(R2,j3)/cell_gas_density(R2);
	v[3][3] = cell_momentum(R2,j4)/cell_gas_density(R2);
	v[4][3] = cell_momentum(R2,j5)/cell_gas_density(R2);
	v[5][3] = cell_gas_gamma(R2);
	v[6][3] = constants->gamma;

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	v[7][3] = cell_electron_internal_energy(R2);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = 0; j < num_chem_species; j++ ) {
		v[num_hydro_vars-num_chem_species-1+j][3] = cell_advected_variable(R2,j)/cell_gas_density(R2);
	}

	if ( cell_level(R1) > cell_level(L1) ) {
		c[0] = 1.0/1.5;
		c[1] = 1.0/1.25;
	} else if ( cell_level(R1) < cell_level(L1) ) {
		c[0] = 1.0/1.25;
		c[1] = 1.0/1.5;
	} else {
		if ( cell_level( L2 ) == cell_level(L1) ) {
			c[0] = 1.0;
		} else {
			c[0] = 1.0/1.25;
		}

		if ( cell_level( R2 ) == cell_level(L1) ) {
			c[1] = 1.0;
		} else {
			c[1] = 1.0/1.25;
		}
	} 

#ifdef GRAVITY_IN_RIEMANN
	g[0] = 0.5*cell_accel( L1, j3 );
	g[1] = 0.5*cell_accel( R1, j3 );

	/* compute fluxes */
	fluxh( dtx, dtx2, v, g, c, f );
#else
	fluxh( dtx, dtx2, v, c, f );
#endif

	if(apply_lapidus_viscosity) lapidus( dtx2, L1, R1, sweep_direction, j3, j4, j5, v, f );
}
	
void hydro_copy_vars( int level, int direction ) {
	int i, j;
	int neighbors[num_neighbors];
	int icell;
	int num_level_cells;
	int *level_cells;

#if nDim != 3
	#error	hydro_copy_vars only works for nDim = 3
#endif

	start_time( WORK_TIMER );

	select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );

	if ( direction == HYDRO_COPY_ALL ) {
#pragma omp parallel for default(none), private(i,icell,j,neighbors), shared(num_level_cells,level_cells,cell_child_oct,cell_vars,direction,backup_hvars,backup_dirty,ref)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			backup_hvar(icell,0) = cell_gas_density(icell);
			backup_hvar(icell,1) = cell_gas_energy(icell);
			backup_hvar(icell,2) = cell_momentum(icell,0);
			backup_hvar(icell,3) = cell_momentum(icell,1);
			backup_hvar(icell,4) = cell_momentum(icell,2);
			backup_hvar(icell,5) = cell_gas_internal_energy(icell);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
			backup_hvar(icell,6) = cell_electron_internal_energy(icell);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

			for ( j = 0; j < num_chem_species; j++ ) {
			  backup_hvar(icell,num_hydro_vars-num_chem_species-2+j) = cell_advected_variable(icell,j);
			}

			ref[icell] = 0.0;
			backup_dirty[icell] = 0;
		}
	} else if ( direction == HYDRO_RESTORE_ALL ) {
#pragma omp parallel for default(none), private(i,icell,j), shared(num_level_cells,level_cells,cell_child_oct,cell_vars,direction,backup_hvars,ref)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			cell_gas_density(icell) = max( 1.0e-30, backup_hvar(icell,0) );
			cell_gas_energy(icell) = max( 1.0e-30, backup_hvar(icell,1) );
			cell_momentum(icell,0) = backup_hvar(icell,2);
			cell_momentum(icell,1) = backup_hvar(icell,3);
			cell_momentum(icell,2) = backup_hvar(icell,4);
			cell_gas_internal_energy(icell) = max( 1.0e-30, backup_hvar(icell,5) );
                                        
#ifdef ELECTRON_ION_NONEQUILIBRIUM
			cell_electron_internal_energy(icell) = backup_hvar(icell,6);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

			for ( j = 0; j < num_chem_species; j++ ) {
				cell_advected_variable(icell,j) = max( 1.0e-30, 
						backup_hvar(icell,num_hydro_vars-num_chem_species-2+j) );
			}
		}
	} else if ( direction == HYDRO_RESTORE_CLEAN ) {
#pragma omp parallel for default(none), private(i,icell,j), shared(num_level_cells,level_cells,cell_child_oct,cell_vars,direction,backup_hvars,ref,backup_dirty)
		for ( i = 0; i < num_level_cells; i++ ) {
			icell = level_cells[i];

			if ( !backup_dirty[icell] ) {                                                                                                           
				cell_gas_density(icell) = max( 1.0e-30, backup_hvar(icell,0) );
				cell_gas_energy(icell) = max( 1.0e-30, backup_hvar(icell,1) );
				cell_momentum(icell,0) = backup_hvar(icell,2);
				cell_momentum(icell,1) = backup_hvar(icell,3);
				cell_momentum(icell,2) = backup_hvar(icell,4);
				cell_gas_internal_energy(icell) = max( 1.0e-30, backup_hvar(icell,5) );
                                        
#ifdef ELECTRON_ION_NONEQUILIBRIUM
				cell_electron_internal_energy(icell) = backup_hvar(icell,6);
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

				for ( j = 0; j < num_chem_species; j++ ) {
					cell_advected_variable(icell,j) = max( 1.0e-30, 
							backup_hvar(icell,num_hydro_vars-num_chem_species-2+j) );
				}
			}
		}
	}		

	cart_free( level_cells );

	end_time( WORK_TIMER );
}

#endif /*HYDRO*/
