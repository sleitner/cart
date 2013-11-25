#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#include <sys/types.h>
#include <unistd.h>

#include "tree.h"
#include "cosmology.h"
#include "particle.h"
#include "sfc.h"
#include "parallel.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "load_balance.h"
#include "run/step.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "iterators.h"
#include "times.h"
#include "timing.h"
#include "units.h"

#include "hydro.h"
#include "hydro_tracer.h"
#include "gravity.h"
#include "density.h"
#include "io.h"
#include "auxiliary.h"
#include "starformation.h"
#include "plugin.h"

#include "rt_debug.h"

#include "disk.h"
void merge_buffer_cell_gas_density_momentum( int level ) ;
int create_star_particle( int icell, float mass, double dt, int type );
int place_star_particle( int icell, float mass, double Zsol, double pos[nDim], double vel[nDim], double pdt, double age, int type );

const int nbuild_mesh = 2;
double pos_central[nDim]={num_grid/2.,num_grid/2.,num_grid/2.};
extern int last_star_id;
double uniform_particles_mass;
int icell_wrt_central(double dispx,double dispy,double dispz){
    double pos[nDim];
    pos[0] = pos_central[0]+dispx;
    pos[1] = pos_central[1]+dispy;
    pos[2] = pos_central[2]+dispz;
    return cell_find_position(pos);
}
void units_set_art(double OmegaM, double h, double Lbox);

float convert_temp_to_ie(float temp,int cell){
    return temp/units->temperature*cell_gas_density(cell)/
        ((cell_gas_gamma(cell)-1)*constants->wmu);
}
float convert_n_to_density(float n){
    return n/constants->cc/units->number_density;
}

/* UNITS */
/* ----------- */
/* Velocity: km/s */
/* Mass: 10^9 Msun */
/* Length: kpc */

/* Disk properties (in above units) */
/* ---------------------------- */
/* Disk scale length: 3.57927 */
/* Disk scale height: 0.357927 */
/* M_DISK: 44.8431 */
/* M_GAS: 8.96862 */
/* M200: 1121.08 */
/* R200: 214.286 */
/* Total mass will be: 1358.5 */
void find_gas_particle_density( int level, double size2, double size_inverse, double pos[nDim], double mass,
				int cell_list[num_children], float mass_assigned[num_children] ) ;
#define N_GAS   1e5 
#define N_HALO  1e5
#define N_DISK  1e5
#define N_BULGE 1e4
const float Tgasmodelat1cc = 1.0e2;
    //double gasdx = 4*constants->kpc/sqrt(N_DISK);
void assign_gas_particle_density( int level, double pos[nDim], double vel[nDim], double mass ) {
    int  k;
    int icell;
    double size2, size_inverse;
    int cell_list[num_children];
    float mass_assigned[num_children];
    double newrho;

    size2 = 0.5*cell_size[level];
    size_inverse = cell_size_inverse[level];
    find_gas_particle_density( level, size2, size_inverse, pos, mass, 
			       cell_list, mass_assigned );
    for ( k = 0; k < num_children; k++ ) {
	icell = cell_list[k];
	if ( icell > -1 ) {
	    newrho = mass_assigned[k] *cell_volume_inverse[level];
	    cell_gas_density(icell) += newrho;
	    cell_momentum(icell,0) += newrho*vel[0];
	    cell_momentum(icell,1) += newrho*vel[1];
	    cell_momentum(icell,2) += newrho*vel[2];
	    cell_gas_gamma(icell) = constants->gamma;
/*   	    cell_gas_internal_energy(icell) = convert_temp_to_ie(Tgasmodelat1cc, icell);   */
	    /* isobaric disk */
 	    cell_gas_internal_energy(icell) = convert_temp_to_ie(Tgasmodelat1cc, icell) * (1/units->number_density) /cell_gas_density(icell); 
/* 	    cell_gas_internal_energy(icell) = u_gas; */
/* 	    cart_error("ugas units?"); */
	    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
	    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
	    cart_assert( cell_gas_internal_energy(icell) > 0);
#ifdef ENRICHMENT 
	    cell_gas_metal_density_II(icell) = constants->Zsun*cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
	    cell_gas_metal_density_Ia(icell) = constants->Zsun*cell_gas_density(icell);
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef RADIATIVE_TRANSFER
	    cell_HI_density(icell) = cell_gas_density(icell)*constants->XH;
	    cell_HII_density(icell) = 0;
	    cell_HeI_density(icell) = cell_gas_density(icell)*constants->XHe;
	    cell_HeII_density(icell) = 0;
	    cell_HeIII_density(icell) = 0;
	    cell_H2_density(icell) = 0;
#endif
	}
    }
}

void find_gas_particle_density( int level, double size2, double size_inverse, double pos[nDim], double mass,
				int cell_list[num_children], float mass_assigned[num_children] ) {
	double corner[nDim];
	double cornerx0, cornerx1, cornery0, cornery1, cornerz0, cornerz1;
	double x, y, z;
	double xs, ys, zs;
	double dx0, dx1, dy0, dy1, dz0, dz1;
	double d00, d01, d10, d11;

	x = pos[0];
	y = pos[1];
	z = pos[2];

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

	dx0 *= mass;
	dx1 *= mass;

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


void assign_gas_model(FILE *fd){
    int i, icell, level;
    float xp,yp,zp, vxp,vyp,vzp, mp, u_gas;
    double pos[nDim], vel[nDim];
    double Zsol = constants->Zsun;
    double mass; 

    for(i=1;i<=N_GAS;i++)
	{
	    fscanf(fd," %g",&xp);
	    fscanf(fd," %g",&yp);
	    fscanf(fd," %g",&zp);
	    fscanf(fd," %g",&vxp);
	    fscanf(fd," %g",&vyp);
	    fscanf(fd," %g",&vzp);
	    fscanf(fd," %g",&mp);
	    fscanf(fd," %g",&u_gas);
	    pos[0] = xp*constants->kpc/units->length + pos_central[0]; 
	    pos[1] = yp*constants->kpc/units->length + pos_central[1];
	    pos[2] = zp*constants->kpc/units->length + pos_central[2];
	    cart_assert(pos[0]<num_grid && pos[0]>0 );
	    cart_assert(pos[1]<num_grid && pos[1]>0 );
	    cart_assert(pos[2]<num_grid && pos[2]>0 );
	    icell = cell_find_position(pos);
	    if(icell>-1 && cell_is_local(icell)){
		vel[0] = vxp * constants->kms / units->velocity ;
		vel[1] = vyp * constants->kms / units->velocity ;
		vel[2] = vzp * constants->kms / units->velocity ;
		mass    = mp*1e9*constants->Msun / units->mass ;
		icell = cell_find_position(pos);
		level = cell_level(icell);
		assign_gas_particle_density( level, pos, vel, mass );
	    }
	}      
}

void init_dm_particles(float mp){
    int i;

    /* setup first specie*/
    num_particle_species = 1; /* not including stars */
    i=0;
    particle_species_indices[i] = 0;
    particle_species_mass[i] = mp ;
    particle_species_indices[i+1] = N_HALO; /* first index of specie */
    particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];

#ifdef UNIFORM_PARTICLES
    num_particle_species++; 
    i++;
    particle_species_mass[i] = uniform_particles_mass ;
    particle_species_indices[i+1] = particle_species_indices[i] + num_grid*num_grid*num_grid ; /* first index of specie */
    particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];
#endif
    
    for ( i = 0; i <= num_particle_species; i++ ) {
	cart_debug("particle_species_mass[%u] = %e", i, particle_species_mass[i] );
	cart_debug("particle_species_indices[%u] = %d", i, particle_species_indices[i] );
	cart_debug("particle_species_num[%u] = %u", i, particle_species_num[i] );
    }

    num_particles_total = particle_species_indices[num_particle_species];
    cart_debug("num_particles_total (nonstar) = %u", num_particles_total );
#ifdef STAR_FORMATION    
    if ( num_particle_species+1 > MAX_PARTICLE_SPECIES ) {
	cart_error("header.Nspecies > MAX_PARTICLE_SPECIES.  Increase and rerun.");
    }

    /* NOW add a specie for stars*/
    num_particle_species++;
    particle_species_mass[num_particle_species-1] = 0.0;
    particle_species_indices[num_particle_species] = particle_species_indices[num_particle_species-1];
    particle_species_num[num_particle_species-1] = 0;
#endif
}


void assign_darkmatter_model(FILE *fd){
    int icell, level, ipart;
	particleid_t i;
    float xp,yp,zp, vxp,vyp,vzp, mp, mp1, u_gas;
    double pos[nDim], vel[nDim];
    float Zsol = constants->Zsun;
    double pdt;
    int current_type;
    
    
    current_type = 0;
    for(i=0;i<N_HALO;i++)
	{
	    fscanf(fd," %e",&xp);
	    fscanf(fd," %e",&yp);
	    fscanf(fd," %e",&zp);
	    fscanf(fd," %e",&vxp);
	    fscanf(fd," %e",&vyp);
	    fscanf(fd," %e",&vzp);
	    fscanf(fd," %e",&mp);
	    
	    pos[0] = xp*constants->kpc/units->length + pos_central[0]; 
	    pos[1] = yp*constants->kpc/units->length + pos_central[1];
	    pos[2] = zp*constants->kpc/units->length + pos_central[2];
	    cart_assert(pos[0]<num_grid && pos[0]>0 );
	    cart_assert(pos[1]<num_grid && pos[1]>0 );
	    cart_assert(pos[2]<num_grid && pos[2]>0 );
	    if(i==0){
		mp1 = mp;
		cart_debug("1st specie's mass %e",mp);
		mp = mp * 1e9*constants->Msun / units->mass ;
		init_dm_particles( mp);
/* 	    }else{ */
/* 		    cart_assert(mp == mp1); */
	    }
	    icell = cell_find_position(pos); //returns -1 if not local or on buffer
	    if(icell>=0 && cell_is_local(icell)){
		vel[0] = vxp * constants->kms / units->velocity ;
		vel[1] = vyp * constants->kms / units->velocity ;
		vel[2] = vzp * constants->kms / units->velocity ;
		
		
		level = cell_level(icell);
		pdt = 0; /* dtl[level];*/

		
		/* Do this for dark matter -- taken from io_cart2.def     */
		ipart = particle_alloc( i );
		cart_assert( ipart > -1 && ipart < num_particles );
		particle_x[ipart][0] = pos[0];
		particle_x[ipart][1] = pos[1];
		particle_x[ipart][2] = pos[2];
		particle_v[ipart][0] = vel[0];
		particle_v[ipart][1] = vel[1];
		particle_v[ipart][2] = vel[2];
		
		particle_t[ipart] = tl[min_level];
		particle_dt[ipart] = pdt;
		particle_mass[ipart] = particle_species_mass[particle_species(i)];
		/* particle_id[ipart] in particle_alloc; id is local
		// particle_level are and insert_particle respectively
		*/

		/* insert particle into cell linked list */
		insert_particle( icell, ipart );
		cart_assert(particle_x[ipart][0]<num_grid && particle_x[ipart][0]>0 );
		cart_assert(particle_x[ipart][1]<num_grid && particle_x[ipart][1]>0 );
		cart_assert(particle_x[ipart][2]<num_grid && particle_x[ipart][2]>0 );
	    }
	}       

    if( particle_species( i-1 ) != current_type ) 
	{ 
	    cart_error("Assertion failed: particle_species(%d)=%d, current_type=%d",i-1,particle_species(i-1),current_type); 
	} 
	
}
    
#ifdef UNIFORM_PARTICLES
void assign_darkmatter_uniform(){
    int i, icell, level, ipart;
    double pos[nDim], vel[nDim];
    float Zsol = constants->Zsun;
    double pdt;
    particleid_t current_id;
    int current_type;
    
    current_id = N_HALO;
    current_type = 1;
	/* dhr - this is 64-bit unsafe for large num_grid */
    for(i=0;i<num_grid*num_grid*num_grid;i++)
	{
	    cell_center_position(i, pos);
	    icell = cell_find_position(pos);
	    vel[0] = 0;
	    vel[1] = 0;
	    vel[2] = 0;

	    level = cell_level(icell);
	    ipart = particle_alloc( current_id );
	    cart_assert( ipart > -1 && ipart < num_particles );
	    particle_x[ipart][0] = pos[0];
	    particle_x[ipart][1] = pos[1];
	    particle_x[ipart][2] = pos[2];
	    particle_v[ipart][0] = vel[0];
	    particle_v[ipart][1] = vel[1];
	    particle_v[ipart][2] = vel[2];

	    if( particle_species( current_id ) != current_type )
		{
		    cart_error("Assertion failed: particle_species(%d)=%d, current_type=%d",current_id,particle_species(current_id),current_type);
		}
	    particle_t[ipart] = tl[min_level];
	    particle_dt[ipart] = 0.0; /* dtl[level]; */
	    particle_mass[ipart] = particle_species_mass[particle_species(current_id)];

	    /* insert particle into cell linked list */
	    insert_particle( icell, ipart );
	    if(particle_x[ipart][0]>=num_grid || particle_x[ipart][0]<=0 ||
	       particle_x[ipart][1]>=num_grid || particle_x[ipart][1]<=0 ||
	       particle_x[ipart][2]>=num_grid || particle_x[ipart][2]<=0 ){
		cart_debug("%e %d %e ", particle_x[ipart][0],num_grid, particle_x[ipart][1]);
	    }

	    cart_assert(particle_x[ipart][0]<num_grid && particle_x[ipart][0]>0 );
	    cart_assert(particle_x[ipart][1]<num_grid && particle_x[ipart][1]>0 );
	    cart_assert(particle_x[ipart][2]<num_grid && particle_x[ipart][2]>0 );
	    
	    current_id++;
	}       
}
#endif /* UNIFORM_PARTICLES */

void assign_star_model(FILE *fd){
    int i, level, icell, ipart ;
    float xp,yp,zp, vxp,vyp,vzp, mp, u_gas;
    double pos[nDim], vel[nDim];
    double Zsol = constants->Zsun;
    double pdt, age;

    int add_particles=0;
    int Alladded_particles=0;

    cart_debug("snl starting assign star");
    last_star_id = particle_species_indices[num_particle_species-1];
    for(i=0;i<N_DISK+N_BULGE;i++)
	{
	    fscanf(fd,"%e",&xp);
	    fscanf(fd,"%e",&yp);
	    fscanf(fd,"%e",&zp);
	    fscanf(fd,"%e",&vxp);
	    fscanf(fd,"%e",&vyp);
	    fscanf(fd,"%e",&vzp);
	    fscanf(fd,"%e",&mp);
	    pos[0] = xp*constants->kpc/units->length + pos_central[0]; 
	    pos[1] = yp*constants->kpc/units->length + pos_central[1];
	    pos[2] = zp*constants->kpc/units->length + pos_central[2];
	    icell = cell_find_position(pos); /* -1 for nonlocal but still buffered cells */
	    if(icell>=0 && cell_is_local(icell)){
		vel[0] = vxp * constants->kms / units->velocity ;
		vel[1] = vyp * constants->kms / units->velocity ;
		vel[2] = vzp * constants->kms / units->velocity ;
		mp = mp * 1e9*constants->Msun / units->mass ;
		level = cell_level(icell);
		pdt = 0.0; /* TIMESTEP IS IMPORTANT -- it cannot be the min_level timestep!! dtl[level]; */
		age = 1e9*constants->yr/units->time; /* initial stars are old */
		ipart = place_star_particle(icell, mp, Zsol, pos, vel, pdt, age, STAR_TYPE_NORMAL );
		cart_assert(particle_x[ipart][0]<num_grid && particle_x[ipart][0]>0 );
		cart_assert(particle_x[ipart][1]<num_grid && particle_x[ipart][1]>0 );
		cart_assert(particle_x[ipart][2]<num_grid && particle_x[ipart][2]>0 );
		add_particles++;
	    }
	}

    cart_debug("stars %d on node ",add_particles);
    MPI_Allreduce(&add_particles, &Alladded_particles, 1, MPI_INT, MPI_SUM, mpi.comm.run );
    cart_debug("all stars %d ",Alladded_particles);
    cart_assert(particle_is_star(ipart));
    cart_debug("snl finished assign star\n");
    
}

#include <string.h>
void read_model_particles(int iter_number){
    FILE *fd;
    char filename[256];
    char a[1000];
    int nlines, level;

    sprintf(filename, "%s/1e5.dat",output_directory);
    fd = fopen( filename, "r" );
    if ( fd == NULL ) {cart_error("Unable to open %s", filename );}

    if(iter_number == 0){
	nlines=0;
	while(!feof(fd))
	    {
		fgets(a,1000,fd);
		nlines++;
	    }
	rewind(fd);
	nlines--;
	cart_debug("reading %d lines from disk file",nlines);
	cart_assert(nlines -N_GAS -N_HALO -N_DISK -N_BULGE == 0);
    }
	
    /* ORDERING IS IMPORTANT: Assignment MUST be done in this order for the file format + iteration*/
    cart_debug("assigning gas, ignoring gas internal energy");
    assign_gas_model(fd);
  
    if(iter_number == 0){
	cart_debug("assigning darkmatter");
	assign_darkmatter_model(fd);

#ifdef UNIFORM_PARTICLES
	cart_assert(uniform_particles_mass>0);
	assign_darkmatter_uniform();
#endif

	cart_debug("assigning stars");
	assign_star_model(fd);
	remap_star_ids();
////////////////////////////////////////////////////// from tree_debug.c
/* 	int count=0; */
/* 	int species_count_total[100]; */
/*         int species_count[100]; */
/* 	int i; */

/*         for ( i = 0; i < num_particle_species; i++ ) { */
/* 	    species_count[i] = 0; */
/* 	    species_count_total[i] = 0; */
/*         } */
/*         for ( i = 0; i < num_particles; i++ ) { */
/* 	    if ( particle_level[i] != FREE_PARTICLE_LEVEL ) { */
/* 		if( particle_id[i] >= num_particles_total ) */
/* 		    { */
/* 			cart_debug("Incorrect particle[%d] id=%d, num_particles_total=%d",i,particle_id[i],num_particles_total); */
/* 		    } */
/* 		count++; */
/* 		species_count[ particle_species( particle_id[i] ) ]++; //snl   */
/* 	    } */
/*         } */

/* 	MPI_Allreduce( species_count, species_count_total, num_particle_species, MPI_INT, MPI_SUM, mpi.comm.run ); */
/* 	MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, mpi.comm.run ); */

/* 	cart_debug(" count %d particles total %d ",count,num_particles_total);  */
/*         for ( i = 0; i < num_particle_species; i++ ) { //snl */
/* 	    cart_debug( "spec %d counted %d supposedly there %d",i, species_count_total[i], particle_species_num[i] ); */
/* 	    cart_assert( species_count_total[i] == particle_species_num[i] ); */
/*         } */
//	exit(1);
//////////////////////////////////////////////////////


    }
	
    fclose(fd);
}

#define NUM_VCIRC_POINTS 10000
double vcirc_dat[NUM_VCIRC_POINTS+1];
double rad_vcirc_dat[NUM_VCIRC_POINTS+1];
void  read_vcirc_file(){
    FILE *fd;
    char filename[256];
    char a[1000];
    int nlines,i;
    

    sprintf(filename, "%s/vcirc.dat",output_directory);
    fd = fopen( filename, "r" );
    if ( fd == NULL ) {cart_error("Unable to open %s", filename );}

    nlines=0;
    while(!feof(fd))
	{
	    fgets(a,1000,fd);
	    nlines++;
	}
    rewind(fd);
    nlines--;
    cart_assert(nlines-NUM_VCIRC_POINTS  == 0);
    
    float rad_kpc, vcirc_dat_kms;
    rad_vcirc_dat[0] = 0;
    vcirc_dat[0] = 0;
    for(i=1;i<=NUM_VCIRC_POINTS;i++){
	fscanf(fd," %g",&rad_kpc);
	fscanf(fd," %g",&vcirc_dat_kms);
	rad_vcirc_dat[i] = rad_kpc*constants->kpc/units->length;
	vcirc_dat[i] = vcirc_dat_kms*constants->kms/units->velocity;
    }
    fclose(fd);
}
double vcirc_from_file( double rad ){
    double dr_bin, dr;
    int bin, bin2;
    dr_bin = rad_vcirc_dat[1]-rad_vcirc_dat[0]; /* constant dr bins */
    bin = (int) (rad / dr_bin);
    /*the upper limit shouldn't happen because of rmax*/
    cart_assert( bin >= 0 && bin < NUM_VCIRC_POINTS ); 

    dr = (rad-rad_vcirc_dat[bin])/dr_bin ;
    bin2 = bin+1 < NUM_VCIRC_POINTS ? bin+1 : bin ;
    return vcirc_dat[bin]*(1-dr) + dr*vcirc_dat[bin2];
}

void analytic_gas_model(int icell){

    double pos[nDim];
    double rpos[nDim];
    const double r0_kpc = 3.43218;
    double r0=r0_kpc*constants->kpc/units->length;
    const double z0_kpc = 0.343218;
    double z0=z0_kpc*constants->kpc/units->length;
    const double rmax = 20*constants->kpc/units->length;
    const double zmax = 3*constants->kpc/units->length;

    const double Mgas_msun = 8.59322e9;
    double Mgas = Mgas_msun * constants->Msun/units->mass;
    
    const double T_analytic = 1e4;

    double rho0 = Mgas/4/M_PI/(r0*r0)/(z0);

    double rad, height, rcyl;
    double vcirc, vx, vy, ex, ey, hypo;

    cell_center_position(icell, pos);

    height = compute_distance_periodic_1d(pos[2], pos_central[2]);
    
    rpos[0]=pos[0];
    rpos[1]=pos[1];
    rpos[2]=pos_central[2];
    rcyl = compute_distance_periodic(rpos, pos_central);

    if(height<zmax && rcyl<rmax){
	rad = compute_distance_periodic(pos, pos_central);
    
	cell_gas_density(icell) = rho0*exp(-rcyl/r0)*exp(-fabs(height)/z0);
	
	vcirc = vcirc_from_file(rad);
	hypo = compute_distance_periodic(rpos,pos_central);
	ex= compute_displacement_periodic_1d(pos_central[0],pos[0])/hypo;
	ey= compute_displacement_periodic_1d(pos_central[1],pos[1])/hypo;
	
	vx= -ey*vcirc;
	vy= ex*vcirc;
	    
	cell_momentum(icell,0) = cell_gas_density(icell) * vx;
	cell_momentum(icell,1) = cell_gas_density(icell) * vy;
	cell_momentum(icell,2) = 0;

	cell_gas_gamma(icell) = constants->gamma;
	cell_gas_internal_energy(icell) =  convert_temp_to_ie(T_analytic, icell);
	cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
	cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
	cart_assert( cell_gas_internal_energy(icell) > 0);
#ifdef ENRICHMENT 
	cell_gas_metal_density_II(icell) = constants->Zsun*cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
	cell_gas_metal_density_Ia(icell) = constants->Zsun*cell_gas_density(icell);
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef RADIATIVE_TRANSFER
	/* start all singly ionized */
	cell_HI_density(icell) = 0;
	cell_HII_density(icell) = cell_gas_density(icell)*constants->XH;
	cell_HeI_density(icell) = 0;
	cell_HeII_density(icell) = cell_gas_density(icell)*constants->XHe;
	cell_HeIII_density(icell) = 0;
	cell_H2_density(icell) = 0;
#endif
    }else{ /* ambient reset */
	cell_gas_density(icell) = convert_n_to_density(1e-6);
	cell_momentum(icell,0) = 0;
	cell_momentum(icell,1) = 0;
	cell_momentum(icell,2) = 0;

	cell_gas_gamma(icell) = constants->gamma;
	cell_gas_internal_energy(icell) =  convert_temp_to_ie(T_analytic, icell);
	cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
	cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
	cart_assert( cell_gas_internal_energy(icell) > 0);
#ifdef ENRICHMENT 
	cell_gas_metal_density_II(icell) = 1e-20*constants->Zsun*cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
	cell_gas_metal_density_Ia(icell) = 1e-20*constants->Zsun*cell_gas_density(icell);
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef RADIATIVE_TRANSFER
	/* start all singly ionized I guess */
	cell_HI_density(icell) = 0;
	cell_HII_density(icell) = cell_gas_density(icell)*constants->XH;
	cell_HeI_density(icell) = 0;
	cell_HeII_density(icell) = cell_gas_density(icell)*constants->XHe;
	cell_HeIII_density(icell) = 0;
	cell_H2_density(icell) = 0;
#endif
	
    }
    
}
void set_analytic_gas_model(){
    int i, level, icell;
    int num_level_cells;
    int *level_cells;
    
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
	    icell = level_cells[i] ;
	    analytic_gas_model( icell ); 
	    cart_assert(cell_gas_density(icell) > 0);
        }
        cart_free( level_cells );
    }

    
    for ( level = max_level - 1; level >= min_level; level-- ) {
        hydro_split_update(level); /* update all non-leafs with their child's value */
        hydro_eos(level);
    }
}



void ambient_conditions( int icell ){
    cell_gas_density(icell) = convert_n_to_density(n_ambient);
    cell_momentum(icell,0) = 0;
    cell_momentum(icell,1) = 0;
    cell_momentum(icell,2) = 0;
    cell_gas_gamma(icell) = constants->gamma;
    cell_gas_internal_energy(icell) =  convert_temp_to_ie(T_ambient, icell);
    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
    cart_assert( cell_gas_internal_energy(icell) > 0);
#ifdef ENRICHMENT 
    cell_gas_metal_density_II(icell) = 1e-20*constants->Zsun*cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
    cell_gas_metal_density_Ia(icell) = 1e-20*constants->Zsun*cell_gas_density(icell);
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef RADIATIVE_TRANSFER
    cell_HI_density(icell) = 0;
    cell_HII_density(icell) = cell_gas_density(icell)*constants->XH;
    cell_HeI_density(icell) = 0;
    cell_HeII_density(icell) = cell_gas_density(icell)*constants->XHe;
    cell_HeIII_density(icell) = 0;
    cell_H2_density(icell) = 0;
#endif
}

void set_ambient_initial_conditions() {
    int i, level, icell;
    int num_level_cells;
    int *level_cells;
    
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
	    icell = level_cells[i] ;
	    ambient_conditions( icell ); 
	    cart_assert(cell_gas_density(icell) > 0);
        }
        cart_free( level_cells );
    }

    
    for ( level = max_level - 1; level >= min_level; level-- ) {
        hydro_split_update(level); /* update all non-leafs with their child's value */
        hydro_eos(level);
    }
}


	
void zero_buffer( int icell ){
    cell_gas_density(icell) = 0;
    cell_momentum(icell,0) = 0;
    cell_momentum(icell,1) = 0;
    cell_momentum(icell,2) = 0;
    cell_gas_gamma(icell) = constants->gamma;
    cell_gas_internal_energy(icell) =  convert_temp_to_ie(T_ambient, icell);
    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
//    cart_assert( cell_gas_internal_energy(icell) > 0);
#ifdef ENRICHMENT 
    cell_gas_metal_density_II(icell) = 1e-20*constants->Zsun*cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
    cell_gas_metal_density_Ia(icell) = 1e-20*constants->Zsun*cell_gas_density(icell);
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
#ifdef RADIATIVE_TRANSFER
    cell_HI_density(icell) = 0;
    cell_HII_density(icell) = cell_gas_density(icell)*constants->XH;
    cell_HeI_density(icell) = 0;
    cell_HeII_density(icell) = cell_gas_density(icell)*constants->XHe;
    cell_HeIII_density(icell) = 0;
    cell_H2_density(icell) = 0;
#endif
}


void set_zero_buffer(){
    int i, level, icell;
    int num_level_cells;
    int *level_cells;
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_BUFFER, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
	    icell = level_cells[i] ;
	    zero_buffer( icell ); 
        }
        cart_free( level_cells );
    }
}


void init_run() {
    int i,j,k;
    int level;
    
    int num_level_cells;
    int *level_cells;
    int icell,idir;
    float frame_momentum[nDim];
    double pos[nDim];
    double rstargas;
    double a0, boxh;
    char *substring;
	
    for ( i = 0; i < nDim; i++ ) {
        refinement_volume_min[i] = 0.0;
        refinement_volume_max[i] = (double)num_grid;
    }

#ifdef COSMOLOGY    
    substring = strstr(jobname, "hiz");
    if(substring != NULL){
	a0=0.25;
	cart_debug("hiz %s",jobname);
	boxh = 4.0;
    }else{
	a0=0.9000;
	boxh = 2.0;
	substring = strstr(jobname, "loz");
	cart_assert(substring != NULL);

	cart_debug("loz %s",jobname);
    }

    cosmology_set(OmegaM,omm0);
    cosmology_set(OmegaL,oml0);
    cosmology_set(OmegaB,omb0);
    cosmology_set(h,hubble);
    cosmology_set(DeltaDC,deltadc);
        
    box_size = boxh;
    auni[min_level] = a0;
    abox[min_level] = a0;
    abox_old[min_level] = a0;
    tl[min_level]=tcode_from_auni(a0);
    tl_old[min_level]=tl[min_level];
        
    units_set_art(omm0,hubble,box_size);
#else
#warning disk is untested out of COSMOLOGY
    units_set(1.0,1.0,1.0); /* mass time length */
#endif
    units_init();
    units_update( min_level );

    /* 
    // set time variables 
    */
#ifdef COSMOLOGY
    for(level=min_level+1; level<=max_level; level++)
    {
        tl[level]       = tl[min_level];
        tl_old[level]   = tl[min_level];
        auni[level]     = auni[min_level];
        abox[level]     = abox[min_level];
        abox_old[level] =  abox[min_level] ; 
    }
        
    cart_debug("tl[min_level] = %f", tl[min_level] );
    cart_debug("au[min_level] = %f", auni[min_level] );
    cart_debug("ab[min_level] = %f", abox[min_level] );
    cart_debug("abo[min_level] = %f", abox_old[min_level] );
    cart_debug("DC mode = %f", cosmology->DeltaDC );
    cosmology_set_fixed();

    set_timestepping_scheme(); //set tl as above and dtl =max_dt...
#endif /* COSMOLOGY */    


    /* describe expected refinement */
    double split_on_8gas,split_on_1dm,code_tot_mass,dm1_mass,gas_mass,star_mass,uniform_particles;
    const double model_gas_particle_mass =  8.96862e-05; /* units 1e9*Msun*/
    const double model_star_particle_mass =  4.484309e-04;
    const double model_dm_particle_mass =  1.309238e-02;

    split_on_8gas = 8 * 8.96862e-05*1e9*constants->Msun/units->mass;
    split_on_1dm = 1.309238e-02*1e9*constants->Msun/units->mass;

    code_tot_mass =  pow(num_grid,3.0);
    gas_mass  = N_GAS   * model_gas_particle_mass  * 1e9*constants->Msun/units->mass;
    dm1_mass  = N_HALO   * model_dm_particle_mass * 1e9*constants->Msun/units->mass;
    star_mass = (N_BULGE+N_DISK) * model_star_particle_mass * 1e9*constants->Msun/units->mass;

#ifdef UNIFORM_PARTICLES
    uniform_particles_mass = (code_tot_mass - dm1_mass - gas_mass - star_mass)/pow(num_grid,3.0);
#else
    uniform_particles_mass = 0;
#endif
    
    
    cart_debug( "gas refinement=8x%e=%e ;  DM refinement=%e ; particles mass to place in every root cell =%e", 
		split_on_8gas/8.0, split_on_8gas,
		split_on_1dm,
		uniform_particles_mass
	);
    
   
    cart_debug("read vcirc");
    read_vcirc_file();
    cart_debug("done reading vcirc");

    /* 
    //set hydro 
    */
    init_particles();
    cart_debug("set particles NULL initially");
    build_particle_list(); /* double check particles */
    
    get_refinement_region();
    build_cell_buffer();
    repair_neighbors();

    cart_debug("refine then read iteratively until refinement conditions matches input");
    cart_assert(nbuild_mesh>1); /* once to set refinement structure, once to populate cells (at least)*/
    
    for(i = 0; i < nbuild_mesh; i++){  
	cart_debug("read particles %d times ==========================",i);
	set_ambient_initial_conditions(); /* resets all cell values */

	set_zero_buffer(); /* zero the buffer */
	read_model_particles(i); 
	for(level=min_level; level<max_level; level++){
	    merge_buffer_cell_gas_density_momentum(level);
	}
	build_refinement_region(-1);
    }
    load_balance();

    /* now reset cell values with analytic model */
    cart_debug("reset hydro to analytic model");
    set_analytic_gas_model();
    
    //   cart_debug( "ambient mass[min_level]=%e ",cell_volume[min_level] * n_ambient*units->mass/constants->Msun );


    hydro_magic( min_level );
    hydro_eos( min_level );
    
    for(level = min_level; level <= max_level; level++){
        cart_debug("updating level %u", level );
        update_buffer_level(level, all_hydro_vars, num_hydro_vars);
    }

    if(nomore_starformation == 1){
	/* make it impossible to produce more stars */
	for ( i = 0; i < nDim; i++ ) {
	    star_formation_volume_min[i] = 0;
	    star_formation_volume_max[i] = 0;
	}
	cart_debug("no star formation after ICs !");
    }else{
	for ( i = 0; i < nDim; i++ ) {
            star_formation_volume_min[i] = 0;
            star_formation_volume_max[i] = num_grid;
        }
    }
	
    /* hydro stuff */
    if ( !buffer_enabled ) {
        cart_debug("building cell buffer");
        build_cell_buffer();
        repair_neighbors();
    }

    hydro_eos( min_level );
    hydro_magic( min_level );
    check_map();
    cart_debug("done with initialization");

    
#ifdef RADIATIVE_TRANSFER
#ifdef RT_DEBUG
  rt_debug.Mode = -1;
  rt_debug.Pos[0] = 0.5;
  rt_debug.Pos[1] = 0.5;
  rt_debug.Pos[2] = 0.5;
#endif
#endif

#ifdef HYDRO_TRACERS
    cart_debug("setting hydro tracers");
    set_hydro_tracers( min_level+1 );
#endif /* HYDRO_TRACERS */

}

int place_star_particle( int icell, float mass, double Zsol, double pos[nDim], double vel[nDim], double pdt, double age, int type) {
	int i;
	int ipart;
	particleid_t id;
	int level;
	float new_density;
	float density_fraction, thermal_pressure;
	float true_mass;

	cart_assert( icell > -1 && icell < num_cells );
	cart_assert( mass > 0.0 );
	cart_assert(age > 0 );

	id = last_star_id + local_proc_id + 1;
	last_star_id += num_procs;
	num_new_stars++;

	ipart = particle_alloc( id );
/* 	cart_debug("%d %d ", ipart, num_star_particles); */
	cart_assert( ipart < num_star_particles );

#ifdef STAR_PARTICLE_TYPES
	star_particle_type[ipart] = type;
#endif

	level = cell_level(icell);
	for ( i = 0; i < nDim; i++ ) {
		particle_x[ipart][i] = pos[i];
	}
	for ( i = 0; i < nDim; i++ ) {
		particle_v[ipart][i] = vel[i];
	}
	
	particle_t[ipart] = tl[level]; 

	particle_dt[ipart] = pdt;

	star_tbirth[ipart] = tl[level] - age; 
	particle_mass[ipart] = mass;
	star_initial_mass[ipart] = mass;
	/* don't count these stars as new star formation in logging */
	total_stellar_initial_mass += mass;
	total_stellar_mass += mass;

#ifdef ENRICHMENT
	star_metallicity_II[ipart] = Zsol;
#ifdef ENRICHMENT_SNIa
	star_metallicity_Ia[ipart] = Zsol;
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
	
	/* insert particle into cell linked list */
	insert_particle( icell, ipart );

	return ipart;
}

#ifdef HYDRO
void merge_buffer_cell_gas_density_momentum( int level );
#endif
#ifdef HYDRO
void merge_buffer_cell_gas_density_momentum( int level ) {
	int i;
	int index, child;
	int icell, proc;

	const int vars_per_cell = 4; /* total vars to send */
	const int updatedbuffer_vars_size = 4; /* vars where buffer is updated after the merge */
	const int updatedbuffer_vars[4] = { 
	    HVAR_GAS_DENSITY, HVAR_MOMENTUM+0,HVAR_MOMENTUM+1, HVAR_MOMENTUM+2 };

	MPI_Request sends[MAX_PROCS];
	MPI_Request receives[MAX_PROCS];
	MPI_Status status;
	MPI_Status statuses[MAX_PROCS];

	float *send_buffer;
	float *recv_buffer;
	int recv_offset[MAX_PROCS];
	int num_send_vars[MAX_PROCS];
	int num_recv_vars[MAX_PROCS];
	int total_send_vars, total_recv_vars;
	int send_offset;
	int send_count, recv_count;

	start_time( COMMUNICATION_TIMER );

	cart_assert( buffer_enabled );

	total_send_vars = total_recv_vars = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		cart_assert( num_remote_buffers[level][proc] >= 0 );
		cart_assert( num_local_buffers[level][proc] >= 0 );

		if ( level == min_level ) {
			num_recv_vars[proc] = vars_per_cell*num_remote_buffers[min_level][proc];
			num_send_vars[proc] = vars_per_cell*num_local_buffers[min_level][proc];
		} else {
			num_recv_vars[proc] = vars_per_cell*num_children*num_remote_buffers[level][proc];
			num_send_vars[proc] = vars_per_cell*num_children*num_local_buffers[level][proc];
		}

		cart_assert( num_send_vars[proc] >= 0 && num_recv_vars[proc] >= 0 );

		total_send_vars += num_send_vars[proc];
		total_recv_vars += num_recv_vars[proc];
	}

	/* set up receives */
	recv_buffer = cart_alloc(float, total_recv_vars );

	recv_count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_recv_vars[proc] > 0 ) {
			recv_offset[proc] = recv_count;
			MPI_Irecv( &recv_buffer[recv_count], num_recv_vars[proc], MPI_FLOAT,
				proc, 0, mpi.comm.run, &receives[proc] );
			recv_count += num_recv_vars[proc];
		} else {
			receives[proc] = MPI_REQUEST_NULL;
		}
	}	

	/* pack cell ids and densities */
	send_buffer = cart_alloc(float, total_send_vars );

	send_offset = send_count = 0;
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( num_send_vars[proc] > 0 ) {
			if ( level == min_level ) {
				cart_assert( num_send_vars[proc] == vars_per_cell*num_local_buffers[level][proc] );
				for ( i = 0; i < num_local_buffers[min_level][proc]; i++ ) {
					icell = local_buffers[min_level][proc][i];

					send_buffer[send_count++] = cell_gas_density(icell);
					send_buffer[send_count++] = cell_momentum(icell,0);
					send_buffer[send_count++] = cell_momentum(icell,1);
					send_buffer[send_count++] = cell_momentum(icell,2);
				}
			} else {
				cart_assert( num_send_vars[proc] == vars_per_cell*num_children*num_local_buffers[level][proc] );

				for ( i = 0; i < num_local_buffers[level][proc]; i++ ) {
					index = local_buffers[level][proc][i];
					cart_assert( index >= 0 && index < num_octs );
					cart_assert( oct_level[index] == level );

					for ( child = 0; child < num_children; child++ ) {
						icell = oct_child( index, child );

						cart_assert( icell >= 0 && icell < num_cells );
						cart_assert( cell_level(icell) == level );

						send_buffer[send_count++] = cell_gas_density(icell);
						send_buffer[send_count++] = cell_momentum(icell,0);
						send_buffer[send_count++] = cell_momentum(icell,1);
						send_buffer[send_count++] = cell_momentum(icell,2);

					}
				}
			}

#ifdef DEBUG
			if ( send_offset + num_send_vars[proc] != send_count ) {
				for ( proc = 0; proc < num_procs; proc++ ) {
					cart_debug("proc = %d, num_send = %d, num_local = %d", 
						proc, num_send_vars[proc], num_local_buffers[level][proc] );
				}
				cart_debug("level = %d", level );
				cart_debug("total_send_vars = %d", total_send_vars );
				cart_debug("i = %d", i );
				cart_debug("send_offset = %d", send_offset );
				cart_debug("num_send_vars[%u] = %d", proc, num_send_vars[proc] );
				cart_debug("num_local_buffers[%d][%d] = %d", level, proc );
				cart_debug("send_count = %d", send_count );
			}
#endif

			cart_assert( send_count <= total_send_vars );
			cart_assert( send_offset + num_send_vars[proc] == send_count );

			MPI_Isend( &send_buffer[send_offset], num_send_vars[proc], MPI_FLOAT,
				proc, 0, mpi.comm.run, &sends[proc] );

			send_offset = send_count;
		} else {
			sends[proc] = MPI_REQUEST_NULL;
		}
	}

	/* process receives as they come in */
	do {
		MPI_Waitany( num_procs, receives, &proc, &status );

		if ( proc != MPI_UNDEFINED ) {
			recv_count = recv_offset[proc];

			if ( level == min_level ) {
				for ( i = 0; i < num_remote_buffers[min_level][proc]; i++ ) {
					icell = root_cell_location( remote_buffers[min_level][proc][i] );

					cell_gas_density(icell) += recv_buffer[recv_count++];
					cell_momentum(icell,0)  += recv_buffer[recv_count++];
					cell_momentum(icell,1)  += recv_buffer[recv_count++];
					cell_momentum(icell,2)  += recv_buffer[recv_count++];
				}
			} else {
				for ( i = 0; i < num_remote_buffers[level][proc]; i++ ) {
					index = remote_buffers[level][proc][i];

					for ( child = 0; child < num_children; child++ ) {
						icell = oct_child( index, child );

						cell_gas_density(icell) += recv_buffer[recv_count++];
						cell_momentum(icell,0)  += recv_buffer[recv_count++];
						cell_momentum(icell,1)  += recv_buffer[recv_count++];
						cell_momentum(icell,2)  += recv_buffer[recv_count++];
					}
				}
			}

			cart_assert( recv_offset[proc] + num_recv_vars[proc] == recv_count );
		}
	} while ( proc != MPI_UNDEFINED );

	cart_free( recv_buffer );

	end_time(COMMUNICATION_TIMER );

	/* now update density variables */
	update_buffer_level( level, updatedbuffer_vars, updatedbuffer_vars_size );

	/* wait for sends */
	start_time( COMMUNICATION_TIMER );
	MPI_Waitall( num_procs, sends, statuses );
	end_time( COMMUNICATION_TIMER );

	cart_free( send_buffer );

}
#endif /* HYDRO */

