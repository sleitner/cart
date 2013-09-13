#include "config.h"

#include <string.h>
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
#include "run/starformation_step.h"
#include "run/hydro_step.h"
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

extern int last_star_id;
double pos_central[nDim]={num_grid/2.,num_grid/2.,num_grid/2.};
char fname_allparticles[256];
char fname_vcirc[256];
int use_analytic_gas_distribution = 0;
static double ZDISK = 0.1; //can be changed by filename

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

void find_gas_particle_density( int level, double size2, double size_inverse, 
                                double pos[nDim], double mass,
				int cell_list[num_children], 
                                float mass_assigned[num_children] ) ;
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
	    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
	    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
	    cart_assert( cell_gas_internal_energy(icell) > 0);
#ifdef ENRICHMENT 
	    cell_gas_metal_density_II(icell) = ZDISK * constants->Zsun*cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
	    cell_gas_metal_density_Ia(icell) = ZDISK * constants->Zsun*cell_gas_density(icell);
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

void find_gas_particle_density( int level, double size2, double size_inverse, 
                                double pos[nDim], double mass,
				int cell_list[num_children], 
                                float mass_assigned[num_children] ) {
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


void assign_gas_model(FILE *fd, int ngas){
    int i, icell, level;
    float xp,yp,zp, vxp,vyp,vzp, mp, u_gas;
    double pos[nDim], vel[nDim];
    double mass; 

    for(i=0;i<ngas;i++)
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

void init_dm_particles(float mp, int nhalo){
    int i;

    /* setup first specie*/
    num_particle_species = 1; /* not including stars */
    i=0;
    particle_species_indices[i] = 0;
    particle_species_mass[i] = mp ;
    particle_species_indices[i+1] = nhalo; /* first index of specie */
    particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];
/*
    num_particle_species++; 
    i++;
    particle_species_mass[i] = next_particles_mass ;
    particle_species_indices[i+1] = particle_species_indices[i] + next_count ; 
    particle_species_num[i] = particle_species_indices[i+1] - particle_species_indices[i];
*/
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


void assign_darkmatter_model(FILE *fd, int nhalo){
    static float mp1;
    int i, icell, level, ipart;
    float xp,yp,zp, vxp,vyp,vzp, mp;
    double pos[nDim], vel[nDim];
    double pdt;
    int current_type;
    current_type = 0;
    for(i=0;i<nhalo;i++)
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
		cart_debug("1st specie's mass %f e9Msun",mp);
		mp = mp * 1e9*constants->Msun / units->mass ;
		init_dm_particles( mp, nhalo);
 	    }else{
                cart_assert(mp == mp1); 
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
		/* particle_id[ipart] in particle_alloc; id is local */
                /* particle_level in insert_particle */
		particle_t[ipart] = tl[min_level];
		particle_dt[ipart] = pdt;
		particle_mass[ipart] = particle_species_mass[particle_species(i)];
		/* insert particle into cell linked list */
		insert_particle( icell, ipart );
		cart_assert(particle_x[ipart][0]<num_grid && particle_x[ipart][0]>0 );
		cart_assert(particle_x[ipart][1]<num_grid && particle_x[ipart][1]>0 );
		cart_assert(particle_x[ipart][2]<num_grid && particle_x[ipart][2]>0 );
	    }
	}       
    if(particle_species(i-1) != current_type) { 
	    cart_error("particle_species(%d)=%d != current_type=%d",
                       i-1,particle_species(i-1),current_type); 
    } 
	
}
    

void assign_star_model(FILE *fd, int nstars){
    int i, level, icell, ipart ;
    float xp,yp,zp, vxp,vyp,vzp, mp;
    double pos[nDim], vel[nDim];
    double Zsol = constants->Zsun;
    double pdt, age;

    int add_particles=0;
    int Alladded_particles=0;

    cart_debug("starting assign star");
    /* place_star adds one before assigning id*/
    last_star_id = particle_species_indices[num_particle_species-1]-1;
    for(i=0;i<nstars;i++)
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
		cart_assert(particle_is_star(ipart));
		cart_assert(particle_x[ipart][0]<num_grid && particle_x[ipart][0]>0 );
		cart_assert(particle_x[ipart][1]<num_grid && particle_x[ipart][1]>0 );
		cart_assert(particle_x[ipart][2]<num_grid && particle_x[ipart][2]>0 );
		add_particles++;
	    }
	}

    cart_debug("stars %d on node ",add_particles);
    MPI_Allreduce(&add_particles, &Alladded_particles, 1, MPI_INT, MPI_SUM, mpi.comm.run );
    cart_debug("all stars %d ",Alladded_particles);
    if(nstars!=0 && !particle_is_star(ipart)){
        cart_error("ipart %d should be a star",ipart );
    }
    cart_debug("finished assign star\n");
    
}

void count_particle_types(int *ngas, int *nhalo, int *nstars){
//
#define NF_GAS 8
#define NF_PARTICLES 7
#define NF_PARTICLES_M 6
#define NF_GAS_M 6
    FILE *fd;
    char buffer[1000];
    char *pch;
    char mfirst[20], mass_str[20];
    int nlines, nf;
    char filename[256];
    int GF_SPEC_COUNT=0;

    sprintf(filename,"%s",fname_allparticles);
    fd = fopen( filename, "r" );
    if ( fd == NULL ) {cart_error("Unable to open %s", filename );}
 
    nlines = 0;
    *ngas = 0;
    *nhalo = 0;
    *nstars = 0;

	nlines = 0;
	while(!feof(fd)){
            nlines++;
            fgets(buffer,1000,fd);
            nf = 0; 
            pch = strtok(buffer," ");
            while(pch != NULL){
                if(GF_SPEC_COUNT == 0 && nf == NF_GAS_M ){
                    sprintf(mass_str,"%s",pch);
                }else if(GF_SPEC_COUNT != 0 && nf == NF_PARTICLES_M ){
                    sprintf(mass_str,"%s",pch);
                }
                pch = strtok(NULL," ");
                nf++;
            }

            /* check which specie we are on */
            if(  nf == NF_GAS ){ /* GAS has an extra (energy) field */
                if(strcmp(mfirst, mass_str)!=0){ /* new mass --> next specie*/
                    GF_SPEC_COUNT=0;
                    cart_debug("spec%d @line %d: mlast %s-> mnew %s ", 
                               GF_SPEC_COUNT,nlines, mfirst, mass_str);
                    sprintf(mfirst,"%s",mass_str);
                }
            }else if(nf == NF_PARTICLES){ /* others have same nf, but different masses */ 
                if(strcmp(mfirst, mass_str)!=0){ /* new mass --> next specie*/
                    GF_SPEC_COUNT++ ;
                    cart_debug("spec%d @line %d: mlast %s-> mnew %s ", 
                               GF_SPEC_COUNT,nlines, mfirst, mass_str);
                    sprintf(mfirst,"%s",mass_str);
                }
            }else if(nf==1){
                cart_debug("done reading");
                break;
            }else{
                cart_error("number of fields=%d is unexpected: %d or %d", NF_GAS, NF_PARTICLES);
            }

            /* increment the count on that spec */
            if(      GF_SPEC_COUNT == 0){
                (*ngas)++;
            }else if(GF_SPEC_COUNT == 1){
                (*nhalo)++;
            }else if(GF_SPEC_COUNT == 2 || GF_SPEC_COUNT == 3){
                (*nstars)++;
            }else{
                cart_error("too many species! %d+1", GF_SPEC_COUNT);
            }
        }
	rewind(fd);
	nlines--;
	cart_debug("reading %d lines from disk file",nlines);
        printf("nlines %d ngas %d nhalo %d nstars %d \n", nlines, *ngas, *nhalo, *nstars );
	cart_assert(nlines -*ngas -*nhalo -*nstars == 0);
	if(!(*ngas   % 100 == 0 ||
             *nhalo  % 100 == 0 || 
             *nstars % 100 == 0)){
            cart_error("why assign less than 100 particles of any species? %d %d %d",*ngas, *nhalo, *nstars);
        }
    fclose(fd);
}
	

void read_model_particles(int iter_number){
    FILE *fd;
    char filename[256];
    static int ngas, nhalo, nstars;

    if(iter_number == 0){
        if( N_GAS==0 && N_STARS==0 && N_HALO==0){
            /* if you don't knoow particle counts try to detect them */
            count_particle_types(&ngas, &nhalo, &nstars);
        }else{
            ngas = N_GAS; 
            nstars = N_STARS; 
            nhalo = N_HALO;
        }
    }
    
    sprintf(filename,"%s",fname_allparticles);
    fd = fopen( filename, "r" );
    if ( fd == NULL ) {cart_error("Unable to open %s", filename );}
    /* 
    // Assignment MUST be this order for the file format + iteration 
    */
    cart_debug("assigning gas, ignoring gas internal energy");
    assign_gas_model(fd, ngas);
    if(iter_number == 0){
	cart_debug("assigning darkmatter");
	assign_darkmatter_model(fd, nhalo);
	/*assign_second_darkmatter_component(nhalo);*/
	cart_debug("assigning stars");
	assign_star_model(fd, nstars);
	remap_star_ids();
    }
    fclose(fd);
}

double vcirc_dat[NUM_VCIRC_POINTS+1];
double rad_vcirc_dat[NUM_VCIRC_POINTS+1];
void  read_vcirc_file(){
    FILE *fd;
    char filename[256];
    char buffer[1000];
    int nlines,i;
    float rad_kpc, vcirc_dat_kms;
    sprintf(filename,"%s", fname_vcirc);
    fd = fopen( filename, "r" );
    if ( fd == NULL ) {cart_error("Unable to open %s", filename );}
    nlines=0;
    while(!feof(fd)){
	    fgets(buffer,1000,fd);
	    nlines++;
    }
    rewind(fd);
    nlines--;
    cart_assert(nlines-NUM_VCIRC_POINTS  == 0);
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
    const double T_analytic = 1e4;
    const double r0_kpc = 3.43218;
    const double z0_kpc = 0.343218;
    const double Mgas_msun = 8.59322e9;

    double Mgas = Mgas_msun * constants->Msun/units->mass;
    double rho0, rmax, zmax ;
    double rad, height, rcyl;
    double vcirc, vx, vy, ex, ey, hypo, r0, z0;

    double pos[nDim];
    double rpos[nDim];

    r0 = r0_kpc*constants->kpc/units->length;
    z0 = z0_kpc*constants->kpc/units->length;
    rho0 = Mgas/(4*M_PI* r0*r0 *z0 ); /* comes from readme */
    rmax = 10*r0; /* truncates disk model*/
    zmax = 10*z0;

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
	cell_gas_metal_density_II(icell) = ZDISK * constants->Zsun*cell_gas_density(icell);
#ifdef ENRICHMENT_SNIa
	cell_gas_metal_density_Ia(icell) = ZDISK * constants->Zsun*cell_gas_density(icell);
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
	/* start all singly ionized */
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
    //   cart_assert( cell_gas_internal_energy(icell) > 0);
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

void cfg_options_fromname(double *a0, double *boxh){
    char prefix[256];
    ZDISK = 0.1;
    /* high/lowz */
    if(      strstr(jobname, "hiz")!= NULL){
        sprintf(prefix, "%s/hiz",output_directory);
	*a0=0.25;
	*boxh = 4.0;
    }else if(strstr(jobname, "loz")!= NULL){
        sprintf(prefix, "%s/loz",output_directory);
        *a0=0.9000;
	*boxh = 2.0;
    }else{
        cart_error("bad jobname %s -- need loz/hiz substring", jobname);
    }

    /* galaxy mass */
    if(      strstr(jobname, "h12")!= NULL){
        sprintf(fname_allparticles,"%s",strcat(prefix, "mh12.dat\0"));
    }else if( strstr(jobname, "h11")!= NULL){
        sprintf(fname_allparticles,"%s",strcat(prefix, "mh11.dat\0"));
        *boxh /= pow(10,.333);
    }else if( strstr(jobname, "h10")!= NULL){
        sprintf(fname_allparticles,"%s",strcat(prefix, "mh10.dat\0"));
        *boxh /= pow(10,.666);
    }else if( strstr(jobname, "h09")!= NULL){
        sprintf(fname_allparticles,"%s",strcat(prefix, "mh09.dat\0"));
        *boxh /= 10.0;
    }else if( strstr(jobname, "h08")!= NULL){
        sprintf(fname_allparticles,"%s",strcat(prefix, "mh08.dat\0"));
        *boxh /= pow(10,10.333);
    }
    cart_debug("jobname %s ; IC name %s",jobname, fname_allparticles);

    /* vicrc defined by file? */
    fname_vcirc[0] = 0;
    if(strstr(jobname,"h12")!=NULL && strstr(jobname,"loz")!=NULL){
        sprintf(fname_vcirc,"%s/vcirc_lozmh12.dat",output_directory);
        use_analytic_gas_distribution = 1; 
	ZDISK = 1.0;
    }
}
void init_run() {
    int i;
    int level;
    double a0, boxh;
    
    cfg_options_fromname(&a0, &boxh);
#ifdef COSMOLOGY    
    
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

    for ( i = 0; i < nDim; i++ ) {
        refinement_volume_min[i] = 0.0;
        refinement_volume_max[i] = (double)num_grid;
    }
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
   
    cart_debug("read vcirc");
    if(strlen(fname_vcirc) != 0 ){
        read_vcirc_file();
    }
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
    if(strlen(fname_vcirc) != 0 && use_analytic_gas_distribution == 1 ){
        cart_debug("reset hydro to analytic model");
        set_analytic_gas_model();
    }

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

    
#ifdef HYDRO_TRACERS
    cart_debug("setting hydro tracers");
    set_hydro_tracers( min_level+1 );
#endif /* HYDRO_TRACERS */

}

int place_star_particle( int icell, float mass, double Zsol, double pos[nDim], double vel[nDim], double pdt, double age, int type) {
	int i;
	int ipart;
	int id;
	int level;

	cart_assert( icell > -1 && icell < num_cells );
	cart_assert( mass > 0.0 );
	cart_assert(age > 0 );

	id = last_star_id + local_proc_id + 1;
	last_star_id += num_procs;
	num_new_stars++;

	ipart = particle_alloc( id );
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
	star_metallicity_II[ipart] = ZDISK * Zsol;
#ifdef ENRICHMENT_SNIa
	star_metallicity_Ia[ipart] = ZDISK * Zsol;
#endif /* ENRICHMENT_SNIa */
#endif /* ENRICHMENT */
	
	/* insert particle into cell linked list */
	insert_particle( icell, ipart );

	return ipart;
}

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

