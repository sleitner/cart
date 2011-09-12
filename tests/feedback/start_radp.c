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
#include "timestep.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "rt_utilities.h"
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

#include "extra/output_slice.h"

#ifdef VIEWDUMP
#include "extra/viewdump.h"
#endif

#include "extra/ifrit.h"

#define OUTLEVEL (max_level)
#define refine_radius	(4.0)

#ifdef COSMOLOGY
#define omm0 1.0
#define oml0 0.0
#define omb0 1.0
#define hubble 1.0
#define deltadc 0.0
#define a0 0.9
//#define boxh (1.0e-6/a0*hubble) //1ppc
#define boxh (1.0e-2/a0*hubble) //10pkpc //for level4
//#define boxh (2.5e-3/a0*hubble) //2.5pkpc //for level2

#define blast_radius    (cell_size[max_level]) 
#define rho0            (1e7)

#ifdef STARFORM
#define E_ambient       (1.0e10)
#else /* STARFORM */
#define E_ambient       (1.0)
#define E_det           (1.0e7*E_ambient)
#define P_ambient	(1.0) 
#endif /* STARFORM */

#else /* !COSMOLOGY */

#define blast_radius    (cell_size[max_level]) 
#define rho0            (1.0e7) 
#define P_ambient	(1.0) 
#define refine_radius	(4.0)

#endif  /* COSMOLOGY */

#ifdef ONE_CELL_IN_THE_CENTER
const int central_cell=num_grid/2;
int icell_central(){
    double pos[nDim];
    pos[0] = central_cell;
    pos[1] = central_cell;
    pos[2] = central_cell;
    return cell_find_position(pos);
}
#ifdef STARFORM 
const double mstar_one_msun = 10;
extern int last_star_id=-1;
#endif
#endif



///////////////////////////////////////////////////////////////
double rad_sedov_pc(double Edum,double rhodum,double timedum){
    return (0.868*pow(Edum/rhodum,0.2)*pow(timedum,0.4))*constants->cm/constants->pc;
}


void units_set_art(double OmegaM, double h, double Lbox);

void refine_level( int cell, int level ) {
    double pos[nDim];
    float r;

    cart_assert( cell >= 0 && cell < num_cells );
    cart_assert( cell_level(cell) == level );
	
    cell_center_position(cell, pos);

    pos[0] -= 0.5*num_grid;
    pos[1] -= 0.5*num_grid;
    pos[2] -= 0.5*num_grid;

    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);

    if ( r < refine_radius*cell_size[level] ) { 
        refinement_indicator(cell,0) = 1.0;	 
    } else { 
        refinement_indicator(cell,0) = 0.0; 
    } 
}


void radp_initial_conditions_one_cell( int icell ) {
    cell_gas_density(icell) = rho0;
    cell_momentum(icell,0) = 0.0;
    cell_momentum(icell,1) = 0.0;
    cell_momentum(icell,2) = 0.0;
    cell_gas_gamma(icell) = constants->gamma;
    
#ifdef STARFORM
    if ( icell == icell_central()){
        cell_gas_density(icell) += mstar_one_msun*constants->Msun/units->mass * cell_volume_inverse[cell_level(icell)];
    } 
    cell_gas_internal_energy(icell) =  E_ambient;
    
#else
    
    if ( icell == icell_central()){
        cell_gas_internal_energy(icell) = 1.0;
    } else { 
        cell_gas_internal_energy(icell) = E_ambient ;
    }
#endif
    
    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
#ifdef ELECTRON_ION_NONEQUILIBRIUM
    cell_electron_internal_energy(icell) = cell_gas_internal_energy(icell)*constants->wmu/constants->wmu_e;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
        
#ifdef ENRICH
    cell_gas_metal_density_II(icell) = 1e-30;
#ifdef ENRICH_SNIa
    cell_gas_metal_density_Ia(icell) = 1e-30;
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
}


void radp_initial_conditions( int icell ) {
    double pos[nDim];
    float r;
    
    cell_gas_density(icell) = rho0;
    cell_momentum(icell,0) = 0.0;
    cell_momentum(icell,1) = 0.0;
    cell_momentum(icell,2) = 0.0;
    cell_gas_gamma(icell) = constants->gamma;
    
    /* now add some energy  */
    cell_center_position(icell, pos); 
    pos[0] -= 0.5*num_grid;
    pos[1] -= 0.5*num_grid;
    pos[2] -= 0.5*num_grid; 
    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    if ( r <= refine_radius*cell_size[min_level] ) {
        cell_gas_internal_energy(icell) = exp( -r*r / (2.0*blast_radius*blast_radius ) );
    } else { 
        cell_gas_internal_energy(icell) =  E_ambient;
    }
    
    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);

#ifdef ELECTRON_ION_NONEQUILIBRIUM
    cell_electron_internal_energy(icell) = cell_gas_internal_energy(icell)*constants->wmu/constants->wmu_e;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
        
#ifdef ENRICH
    cell_gas_metal_density_II(icell) = 1e-30;
#ifdef ENRICH_SNIa
    cell_gas_metal_density_Ia(icell) = 1e-30;
#endif /* ENRICH_SNIa */
#endif /* ENRICH */
}


void set_radp_initial_conditions() {
    int i;
    int level;
    int num_level_cells;
    int *level_cells;
    
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
#ifdef ONE_CELL_IN_THE_CENTER
            radp_initial_conditions_one_cell( level_cells[i] );
#else 
            radp_initial_conditions( level_cells[i] );
#endif /* ONE_CELL_IN_THE_CENTER */
        }
        cart_free( level_cells );
    }
    
#ifndef STARFORM
    double scale_energy;
    
///////////// get E_det normalization
    scale_energy = 0.0;
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            icell = level_cells[i];

            if ( cell_is_leaf(icell) ) {
                scale_energy += cell_gas_internal_energy(icell)*cell_volume[level];
            }
        }
        cart_free( level_cells );
    }
    scale = scale_energy;
    MPI_Allreduce( &scale, &scale_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    
///////////// normalize
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            icell = level_cells[i];
            cell_gas_internal_energy(icell) =
                E_det*cell_gas_internal_energy(icell)/scale_energy
                + P_ambient / (cell_gas_gamma(icell)-1.0);
               
            
            //add a uniform pressure component
            cell_gas_pressure(icell) = cell_gas_internal_energy(icell)
                * (cell_gas_gamma(icell)-1.0);
            cell_gas_energy(icell) = cell_gas_internal_energy(icell)
                + cell_gas_kinetic_energy(icell);
        }
	
        cart_free( level_cells );
    }
#endif /* STARFORM */
    
    for ( level = max_level - 1; level >= min_level; level-- ) {
        hydro_split_update(level);
    }
}



void init_run() {
    int i;
    int level;
    int num_level_cells;
    int *level_cells;
        
    for ( i = 0; i < nDim; i++ ) {
        refinement_volume_min[i] = 0.0;
        refinement_volume_max[i] = (double)num_grid;
    }

#ifdef COSMOLOGY    
    cosmology_set(OmegaM,omm0);
    cosmology_set(OmegaL,oml0);
    cosmology_set(OmegaB,omb0);
    cosmology_set(h,hubble);
    cosmology_set(DeltaDC,deltadc);
        
    box_size = boxh;
    auni_init = a0;
    auni[min_level] = a0;
    abox[min_level] = a0;
    t_init = tcode_from_auni(auni_init);
        
    units_set_art(omm0,hubble,box_size);
#else
    units_set(1.0,1.0,1.0); //mass time length
#endif
    
    units_reset();
    units_update( min_level );
        
////////////////////////////////////////////////////
    /* build buffer */
    build_cell_buffer();
    cart_debug("built cell buffer");
    repair_neighbors();
    cart_debug("repaired neighbors");
        
        
    /* setup initial refinement for ICs to be placed on */
    for ( level = min_level; level < max_level; level++ ) {
        cart_debug("refining level %u", level );

        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        cart_debug("num_level_cells = %u", num_level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            refine_level( level_cells[i], level );
        }
        cart_free( level_cells );
        cart_debug("about to refine level %u", level );
        refine(level);
    }

    build_cell_buffer();
    repair_neighbors();
        
    cart_debug("setting initial conditions on root level");
    set_radp_initial_conditions();
    
    cart_debug( "mass[min_level]=%e ",cell_volume[min_level] * rho0*units->mass/constants->Msun );
        
#ifdef HYDRO_TRACERS
    cart_debug("setting hydro tracers");
    set_hydro_tracers( min_level+1 );
#endif /* HYDRO_TRACERS */
        
    cart_debug("set initial conditions");
        
    hydro_magic( min_level );
    hydro_eos( min_level );
////////////////////////////////////////////////////////
        
    /* set time variables */
    tl[min_level] = t_init;
    cart_debug("t_init=%e", t_init);
    
#ifdef COSMOLOGY
    for(level=min_level+1; level<=max_level; level++)
    {
        tl[level] = tl[min_level];
        auni[level] = auni[min_level];
        abox[level] = abox[min_level];
    }
    cosmology_set_fixed();
        
    cart_debug("tl[min_level] = %f", tl[min_level] );
    cart_debug("au[min_level] = %f", auni[min_level] );
    cart_debug("ab[min_level] = %f", abox[min_level] );
    cart_debug("DC mode = %f", cosmology->DeltaDC );
#endif /* COSMOLOGY */    

    
#ifdef PARTICLES
    for(i=0; i<num_particles; i++) if(particle_level[i] != FREE_PARTICLE_LEVEL)
    {
//  We set the step to 0 so that the first leapfrog step is correct
        particle_t[i] = tl[min_level];
        particle_dt[i] = 0.0;
    }
#endif /* PARTICLES */
        
    for(level = min_level; level <= max_level; level++){
        cart_debug("updating level %u", level );
        update_buffer_level(level, all_hydro_vars, num_hydro_vars);
    }
        
    cart_debug("done updating initial conditions");
        
#ifdef PARTICLES
//need zoom+DM for buildmesh:    build_mesh(); //defines refinement volume 
    for ( i = 0; i < nDim; i++ ) {
        refinement_volume_min[i] = 0;
        refinement_volume_max[i] = (double)num_grid;
    }

    
#ifdef STARFORM
#ifdef ONE_CELL_IN_THE_CENTER
    int icell;
    float dm_star;
    double old_density, density_fraction;
    
    icell = icell_central();
    level = cell_level(icell);
    dm_star = mstar_one_msun*constants->Msun/units->mass;
    
    /* for ( i = 0; i < nDim; i++ ) {  */
    /*     star_formation_volume_min[i] = pos[i];   */
    /*     star_formation_volume_max[i] = pos[i]+cell_size(level);  */
    /* }  */
#ifdef LOG_STAR_CREATION
    log_star_creation( icell, dm_star, FILE_RECORD);
#endif
    num_particle_species=1;
    num_particles_total=0;
    cart_debug("snl1 %e ", cell_gas_density(icell));
    if(dm_star!=0){
        create_star_particle(icell, dm_star);
    }
    cart_debug("snl1 %e ", cell_gas_density(icell));

    remap_star_ids();
    
    /* reset ICs where star formed -- may have removed way too much mass*/
    if(dm_star < 0.6*cell_gas_density(icell)*cell_volume[cell_level(icell)]){
        old_density = cell_gas_density(icell) +  dm_star * cell_volume_inverse[level];
        density_fraction = old_density / cell_gas_density(icell);
        cart_debug("snl1density_fraction=%e %e %e",density_fraction, old_density,cell_gas_energy(icell+1)  );
        cell_gas_energy(icell) *= density_fraction;
        cell_gas_internal_energy(icell) *= density_fraction;
        cell_gas_pressure(icell) *= density_fraction;
        cell_momentum(icell,0) *= density_fraction;
        cell_momentum(icell,1) *= density_fraction;
        cell_momentum(icell,2) *= density_fraction;
        for ( i = 0; i < num_chem_species; i++ ) {
            cell_advected_variable(icell,i) *= density_fraction;
        }
    }else{
        cell_gas_density(icell) = rho0;
        cell_momentum(icell,0) = 0.0;
        cell_momentum(icell,1) = 0.0;
        cell_momentum(icell,2) = 0.0;
        cell_gas_internal_energy(icell) =  E_ambient;
        cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
        cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
        cart_debug("snl1=%e %e ",cell_gas_energy(icell),cell_gas_pressure(icell) );
    }
    
    cell_gas_density(icell) = rho0;
    
    /* make it impossible to produce more stars */
    for ( i = 0; i < nDim; i++ ) {
        star_formation_volume_min[i] = 0;
        star_formation_volume_max[i] = 0;
    }
    

#ifdef BLASTWAVE_FEEDBACK
    init_blastwave(icell);
#endif /* BLASTWAVE_FEEDBACK */

#else /* ONE_CELL_IN_THE_CENTER */
    for ( i = 0; i < nDim; i++ ) {
        star_formation_volume_min[i] = refinement_volume_min[i];
        star_formation_volume_max[i] = refinement_volume_max[i];
    }
    
#endif /* ONE_CELL_IN_THE_CENTER */
#endif /* STARFORM */
    
#endif /* PARTICLES */
        
    if ( !buffer_enabled ) {
        cart_debug("building cell buffer");
        build_cell_buffer();
        repair_neighbors();
    }
    
        
    cart_debug("done with initialization");
        
    check_map();
    /* for ( level = min_level; level < max_level; level++ ) { */
    /*     select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells ); */
    /*     for ( i = 0; i < num_level_cells; i++ ) { */
    /*         cart_debug("blah %e",cell_gas_density( level_cells[i])); */
    /*     } */
    /*     cart_free( level_cells ); */
    /* } */
    /* exit(-1); */


}




/////////////////////////////////////////////////////////////////////////////////
/////////////////////////   OUTPUT    ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
int const endian_test=-99;
int const slice_axis_x=0, slice_axis_y=1, slice_axis_z=2;



/* radial binning for analysis */
#define num_bins        (num_grid*(1<<(max_level-min_level)))
#define bin_width       ((float)num_grid/((float)num_bins))

float intPt;
float radii[num_bins];
float vel[num_bins];
float pressure[num_bins];
float rho[num_bins];
float avgs[num_bins];
float sedov_analytic[num_bins];

void radial_average( int cell, int level ) {
    double pos[nDim];
    float r1, r, v;
    float amt, cur_r;
    int bin;

    if ( cell_is_leaf( cell ) ) {
        cell_center_position(cell, pos);
        r = sqrt( (pos[0]-((float)num_grid/2.0))*(pos[0]-((float)num_grid/2.0))+
                  (pos[1]-((float)num_grid/2.0))*(pos[1]-((float)num_grid/2.0))+
                  (pos[2]-((float)num_grid/2.0))*(pos[2]-((float)num_grid/2.0)));

        /* this needs to be converted to v_r */
        v = sqrt( cell_momentum(cell,0)*cell_momentum(cell,0) +
                  cell_momentum(cell,1)*cell_momentum(cell,1) +
                  cell_momentum(cell,2)*cell_momentum(cell,2) ) / cell_gas_density(cell);

        /* 1) determine which bin left hand edge of cell is in
         * 2) go outwards in bins until right hand edge of bin > right hand edge of cell
         * 3) apply proper amounts to each cell */

        bin = (int)((r-0.5*cell_size[level])/bin_width);

        cart_assert( bin >= 0 && bin < num_bins );

        cur_r = r-0.5*cell_size[level];
        amt = 1.0;

        while ( (float)(bin+1)*bin_width < (r+0.5*cell_size[level]) ) {
            r1 = ( (float)(bin+1)*bin_width - cur_r )/cell_size[level];

            vel[bin] += r1*v;
            rho[bin] += r1*cell_gas_density(cell);
            pressure[bin] += r1*cell_gas_pressure(cell);
            avgs[bin] += r1;
                                                                                
            amt -= r1;
            cart_assert( amt >= 0.0 );
                              
            cur_r = (float)(bin+1)*bin_width;
            bin++;
        }

        if ( amt > 0.0 && bin < num_bins - 1 ) {
            /* apply amt to last bin */
            r1 = amt;
            vel[bin] += r1*v;
            rho[bin] += r1*cell_gas_density(cell);
            pressure[bin] += r1*cell_gas_pressure(cell);
            avgs[bin] += r1;
        }
    }
}





#ifdef USER_PLUGIN
void outslice();
plugin_t outslicePlugin = {NULL};
const plugin_t* add_plugin(int id){
    if(id==0){
        outslicePlugin.RunBegin = outslice;
        outslicePlugin.GlobalStepEnd = outslice;
        return &outslicePlugin;
    }else{
        return NULL;
    }
}

void outslice() {
    int i, j;
    char filename[128];
    FILE *RADP;
    float reduced_rho[num_bins];
    float reduced_pressure[num_bins];
    float reduced_vel[num_bins];
    float reduced_avgs[num_bins];
    int level;
    int num_level_cells;
    int *level_cells;

    FILE *output;
    int icell;
    const int nvars = 4;
    const int nbin1 = 128;
    int varid[] = { HVAR_PRESSURE, HVAR_GAS_DENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC };
    int nbin[] = { nbin1, nbin1, nbin1 };
    double bb[6];
    float tot_momentum=0,dPt=0 ;

    cart_debug("output feedback test:");
        
    bb[0] = bb[2] = bb[4] = num_grid*(0.5-0.125);
    bb[1] = bb[3] = bb[5] = num_grid*(0.5+0.125);
        
    sprintf(filename,"%s/out.%05d.bin",output_directory,step);
    ifrit.OutputMesh(filename,max_level,nbin,bb,nvars,varid);
        
        
    /* now dump the total gas momentum */
    tot_momentum=0, tot_momentum=0;
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            j=level_cells[i];
            tot_momentum += sqrt( cell_momentum(j,0)*cell_momentum(j,0) +
                                  cell_momentum(j,1)*cell_momentum(j,1) +
                                  cell_momentum(j,2)*cell_momentum(j,2) )*cell_volume[level];
            dPt += cell_gas_pressure(j)* 6*cell_size[level]*cell_size[level];
        }
        cart_free( level_cells );
    }
    MPI_Reduce( MPI_IN_PLACE, &tot_momentum, 1, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
    MPI_Reduce( MPI_IN_PLACE, &dPt, 1, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
            
    if ( local_proc_id == MASTER_NODE ) {
        sprintf(filename, "%s/mom_%s.dat", output_directory, jobname );
        RADP = fopen(filename,"a+");
        fprintf(RADP, "%e %e %e\n", tl[min_level]-t_init , tot_momentum,dPt*(tl[min_level]-t_init) );
        fclose(RADP);
    }


    /* now dump the radial profiles */
    for ( i = 0; i < num_bins; i++ ) {
        radii[i] = ((float)i + 0.5) * bin_width;
        vel[i] = 0.0;
        pressure[i] = 0.0;
        rho[i] = 0.0;
        avgs[i] = 0.0;
        sedov_analytic[i] = 0.0;
    }

#ifndef STARFORM
    double rad_an, curr_time, init_E, init_rho;
    int bin_an;
    cart_debug("code E: %e rho: %e t_code: %e; dm0=%e;", E_det, rho0, (tl[min_level]-t_init), rho0*cell_volume[min_level] );
    curr_time = (tphys_from_tcode(tl[min_level]) - tphys_from_tcode(t_init))
        *constants->yr/constants->s;
    init_E = E_det*units->energy/constants->erg;
    init_rho = rho0*units->density/constants->gpercc;
    cart_debug("physical: t=%es ; E=%eerg ; rho=%e#/cc",
               curr_time, init_E, rho0*units->number_density/constants->cc);
    rad_an = rad_sedov_pc(init_E, init_rho, curr_time);
    rad_an = rad_an*constants->pc/units->length;
    bin_an = (int)(rad_an/bin_width);
    cart_debug("time=%eyr sedov radius=%epc ",
               curr_time*constants->s/constants->yr, rad_an*units->length/constants->pc);
    cart_debug("bin_an:%d, %f %f\n",bin_an, rad_an*units->length/constants->pc,bin_width);
    cart_assert( bin_an >= 0 && bin_an < num_bins );
    sedov_analytic[bin_an]=1.0;
#endif
    
    cart_debug("\n\n");
    icell=icell_central();
    cart_debug("--=physical conditions(cent cell=%d)=--",icell);
    cart_debug("rho=%e[#/cc] %e[g/cc] %e[code]", rho0*units->number_density/constants->cc, rho0*units->density/constants->gpercc, rho0);
    cart_debug("T=%e[K] %e[code]", cell_gas_temperature(icell)*units->temperature/constants->K,cell_gas_temperature(icell) );
    cart_debug("P=%e[ergs/cc] %e[code]", cell_gas_pressure(icell)*units->energy_density/constants->barye, cell_gas_pressure(icell));
    cart_debug("crho[%d]=%e[#/cc] %e[g/cc] %e[code]", icell,cell_gas_density(icell)*units->number_density/constants->cc, cell_gas_density(icell)*units->density/constants->gpercc, cell_gas_density(icell));
    cart_debug("mstar=%e[Msun] %e[code] ",mstar_one_msun,mstar_one_msun/units->mass*constants->Msun);
    cart_debug("cell_max[pc]=%e cell_min[pc]%e ",cell_size[min_level]*units->length/constants->pc,cell_size[cell_level(icell)]*units->length/constants->pc);
    
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            radial_average( level_cells[i], level );
        }
        cart_free( level_cells );
    }

    MPI_Reduce( rho, reduced_rho, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
    MPI_Reduce( pressure, reduced_pressure, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
    MPI_Reduce( vel, reduced_vel, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );
    MPI_Reduce( avgs, reduced_avgs, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, MPI_COMM_WORLD );

    if ( local_proc_id == MASTER_NODE ) {
        sprintf(filename, "%s/%s_stp%04u.txt", output_directory, jobname, step );
        RADP = fopen(filename,"w");
        fprintf(RADP, "#time=%e \n",tl[min_level]-t_init);
        for ( i = 0; i < num_bins/2; i++ ) {
            if ( reduced_avgs[i] > 0.0 ) {
                reduced_vel[i] /= reduced_avgs[i];
                reduced_pressure[i] /= reduced_avgs[i];
                reduced_rho[i] /= reduced_avgs[i];
                
                fprintf(RADP, "%e %e %e %e %e\n", radii[i]*units->length/constants->pc, reduced_rho[i], reduced_pressure[i], reduced_vel[i], sedov_analytic[i] );
            }
        }
        fclose(RADP);
    }


        
#ifdef VIEWDUMP
    sprintf(filename, "%s/%s_%04u.v", output_directory, jobname, step );
    viewdump( filename, max_level, ((float)(num_grid)/2.0)+0.5*cell_size[max_level], 2, DUMP_HVARS, CELL_TYPE_LOCAL );
#endif

#ifdef HYDRO_TRACERS
    sprintf( filename, "%s/tracers_%04u.dat", output_directory, step );
    write_hydro_tracers( filename );
#endif /* HYDRO_TRACERS */

    /* output a 2-d slice through the center of the box */
    if ( local_proc_id == MASTER_NODE ) {
        double slice_region_hsize=500; //100pc
        slice_region_hsize *= constants->pc/units->length;
        cart_debug("original ouput size=%e [pc]",slice_region_hsize*units->length/constants->pc);
        int nsgrid_half=(int)slice_region_hsize/cell_size[OUTLEVEL];
        slice_region_hsize=nsgrid_half*cell_size[OUTLEVEL]+cell_size[OUTLEVEL]/2.;
        cart_debug("adjusted ouput size=%e [code] %d[cells]",2*slice_region_hsize,2*nsgrid_half+1);
        cart_debug("adjusted ouput size=%e [pc]",2*slice_region_hsize*units->length/constants->pc);
        
#ifdef ONE_CELL_IN_THE_CENTER
        double pos[3];
        icell=icell_central();
        cell_center_position(icell,pos);
#else
        for(i=0;i<nDim; i++){
            pos[i]=num_grid/2;
        }
        cell_center_position(cell_find_position(pos),pos);
#endif
        /* for(i=0;i<nDim; i++){ */
        /*     //now everything in a given cell will be assigned to that cell by pos/dx */
        /*     // |00x00|1111| */
        /*     pos[i]-=cell_size[OUTLEVEL]/2.0; */
        /* } */
        cart_debug("snl1:%d%d%d",axis_direction[slice_axis_z][0],axis_direction[slice_axis_z][1], slice_axis_z);
        
        sprintf( filename, "%s/%s_slice_%04u.dat", output_directory, jobname, step );
        output = fopen( filename, "w" );
        cart_debug("dumping plane");
        dump_plane(OUT_CELL_DENSITY, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_INTERNAL_ENERGY, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_MOMENTUM+axis_direction[slice_axis_z][0], OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_MOMENTUM+axis_direction[slice_axis_z][1], OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_MOMENTUM+slice_axis_z, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_SOUNDSPEED, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_MACH, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_LEVEL, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        fclose(output);
    }else{
        cart_error("output not MPI_parallel!!");
    }

}


#endif /*  USER_PLUGIN */

void run_output(){
}
