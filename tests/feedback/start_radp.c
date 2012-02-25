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

#include "start_radp.h"
void create_star_particle( int icell, float mass, int type );
    

const int central_cell=num_grid/2;
double tot_energy0=0;
double tot_momentum0=0;
/* const int NSTARS=1; */
/* const int NSTARS=8; */
/* const int NSTARS=27; */
/* const int NSTARS=64; */

float adv_velocity[3] ;

void ic_star();
void ic_star_spread(int icell, float dmstar);
void ic_refine_levels();
double advection_momentum(int icell,int idir);


int icell_central(double dispx,double dispy,double dispz){
    double pos[nDim];
    pos[0] = central_cell+dispx;
    pos[1] = central_cell+dispy;
    pos[2] = central_cell+dispz;
    return cell_find_position(pos);
}


void units_set_art(double OmegaM, double h, double Lbox);

float convert_temp_to_ie(float temp,int cell){
    return temp*constants->K/units->temperature*cell_gas_density(cell)/
        ((cell_gas_gamma(cell)-1)*constants->wmu);
}
float convert_n_to_density(float n){
    return n/constants->cc/units->number_density;
}


double advection_momentum(int icell,int idir)
{
  double btt; 
  btt = box_traverse_time;
    if(btt  > 0){
        if(idir==0){
            return cell_gas_density(icell) *
                boxh/btt*constants->Mpc/constants->yr/units->velocity;
        }else if(idir==1){
            return 0;
/*             return cell_gas_density(icell)/2 * */
/*                 boxh/btt*constants->Mpc/constants->yr/units->velocity; */
        }else{
            return 0;
        }
    }else{
        return 0;
    }
}
void star_cell_conditions( int icell ) {
    cell_gas_density(icell) = convert_n_to_density(n_h2);
    cart_debug("snl1 1e8 ish: %e", cell_gas_density(icell));
    cell_momentum(icell,0) = advection_momentum(icell,0);
    cell_momentum(icell,1) = advection_momentum(icell,1);
    cell_momentum(icell,2) = advection_momentum(icell,2);
    cell_gas_gamma(icell) = constants->gamma;
    
    cell_gas_internal_energy(icell) =  convert_temp_to_ie(T_h2,icell);
    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);

#ifdef ENRICH 
    cell_gas_metal_density_II(icell) = constants->Zsun*cell_gas_density(icell);
#ifdef ENRICH_SNIa
    cell_gas_metal_density_Ia(icell) = constants->Zsun*cell_gas_density(icell);
#endif
#endif
}

void ambient_conditions( int icell ) {
    cell_gas_density(icell) = convert_n_to_density(n_ambient);
    cell_momentum(icell,0) = advection_momentum(icell,0);
    cell_momentum(icell,1) = advection_momentum(icell,1);
    cell_momentum(icell,2) = advection_momentum(icell,2);
    cell_gas_gamma(icell) = constants->gamma;
    
    cell_gas_internal_energy(icell) =  convert_temp_to_ie(T_ambient,icell);
    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);

#ifdef ENRICH 
    cell_gas_metal_density_II(icell) = constants->Zsun*cell_gas_density(icell);
#ifdef ENRICH_SNIa
    cell_gas_metal_density_Ia(icell) = constants->Zsun*cell_gas_density(icell);
#endif
#endif
}



void set_radp_initial_conditions() {
    int i;
    int level;
    int num_level_cells;
    int *level_cells;
    
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            ambient_conditions( level_cells[i] );
        }
        cart_free( level_cells );
    }
    
    for ( level = max_level - 1; level >= min_level; level-- ) {
        hydro_split_update(level); //update all non-leafs with their child's value
        hydro_eos(level);
    }
}

void    ic_star_spread(int icell, float dm_star){
    double density_fraction;
    double new_density;
    int i,level;
    level = cell_level(icell);
    int neighbors[num_neighbors];
    
#ifdef LOG_STAR_CREATION
    log_star_creation( icell, dm_star, FILE_RECORD);
#endif

    /* star particle can form and takes properties from gas */
    star_cell_conditions( icell );
    // particle_v[ipart][i] = cell_momentum(icell,i) / cell_gas_density(icell);
    new_density = cell_gas_density(icell) + dm_star * cell_volume_inverse[cell_level(icell)];
    density_fraction = new_density / cell_gas_density(icell);
    cell_gas_density(icell) = new_density;
    cell_momentum(icell,0) *= density_fraction;
    cell_momentum(icell,1) *= density_fraction;
    cell_momentum(icell,2) *= density_fraction;
#ifdef ENRICH
    cell_gas_metal_density_II(icell) *= density_fraction;
#ifdef ENRICH_SNIa
    cell_gas_metal_density_Ia(icell) *= density_fraction;
#endif
#endif
    
    if(dm_star!=0){
        cart_debug("star mass %e",dm_star);
        create_star_particle(icell, dm_star, STAR_TYPE_NORMAL);
    }
    // cart_debug("snl1:npspec npt %d %d", num_particle_species,num_particles_total);
    remap_star_ids();

    /* set cells neighboring star to have the star-cell properties*/
    cell_all_neighbors( icell, neighbors );
    star_cell_conditions( icell );
    for(i=0;i<num_neighbors;i++){
        star_cell_conditions( neighbors[i] );
    }
        
    
    
    hydro_split_update(cell_level(icell)); //update all non-leafs with their child's value
    
#ifdef BLASTWAVE_FEEDBACK
    init_blastwave(icell);
#endif /* BLASTWAVE_FEEDBACK */
}

void ic_refine_levels(){
  int i, j, level;
    int num_level_cells;
    int *level_cells;

  /*
  //  If HYDRO is on, refinement now requires that the mesh contained
  //  the valid data. Hence, we fill the root cells with some dummy
  //  data.
  */
  for(i=0; i<num_root_cells; i++)
    {
      cell_gas_density(i) = 1;
      cell_momentum(i,0) = 0;
      cell_momentum(i,1) = 0;
      cell_momentum(i,2) = 0;

      cell_gas_gamma(i) = constants->gamma;
      cell_gas_internal_energy(i) =  1;
      cell_gas_pressure(i) = cell_gas_internal_energy(i)*(constants->gamma-1);
      cell_gas_energy(i) = cell_gas_internal_energy(i);

#ifdef RADIATIVE_TRANSFER
      cell_HI_density(i) = cell_gas_density(i)*constants->XH;
      cell_HII_density(i) = 0;
      cell_HeI_density(i) = cell_gas_density(i)*constants->XHe;
      cell_HeII_density(i) = 0;
      cell_HeIII_density(i) = 0;
      cell_H2_density(i) = 0;
#endif
    }

    for(level=min_level; level<max_level; level++){
        modify( level, OP_REFINE );
        if((!refinement_indicator[SPATIAL_INDICATOR].use[level])){
            cart_error("need spatial refinement");
        }
    }
    build_cell_buffer();
    repair_neighbors();
    load_balance();
}


////////////////////////////////////////////////////
                   /* init_run */
//////////////////////////////////////////////////// 
void init_run() {
    int i,j,k;
    int level;
    
    int num_level_cells;
    int *level_cells;
    int icell,idir;
    float frame_momentum[nDim];

    
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
    t_init = tcode_from_auni(auni_init);
    auni[min_level] = a0;
    abox[min_level] = a0;
    abox_old[min_level] = a0;
    tl_old[min_level]   = t_init;
        
    units_set_art(omm0,hubble,box_size);
#else
    units_set(1.0,1.0,1.0); //mass time length
#endif
    
    units_reset();
    units_update( min_level );
        
////////////////////////////////////////////////////
    /* set hydro */
    cart_debug("built cell buffer");
    cart_debug("repaired neighbors");
    build_cell_buffer();
    repair_neighbors();
    ic_refine_levels(); 
    
        
    cart_debug("setting initial conditions on root level");
    set_radp_initial_conditions();
    
    cart_debug( "ambient mass[min_level]=%e ",cell_volume[min_level] * n_ambient*units->mass/constants->Msun );
        
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
        tl[level]       = tl[min_level];
        tl_old[level]   = tl[min_level];
        auni[level]     = auni[min_level];
        abox[level]     = abox[min_level];
        abox_old[level] =  abox[min_level] ; 
    }
    auni_init = abox[min_level];
    cosmology_set_fixed();
        
    cart_debug("tl[min_level] = %f", tl[min_level] );
    cart_debug("au[min_level] = %f", auni[min_level] );
    cart_debug("ab[min_level] = %f", abox[min_level] );
    cart_debug("abo[min_level] = %f", abox_old[min_level] );
    cart_debug("DC mode = %f", cosmology->DeltaDC );
#endif /* COSMOLOGY */    


    num_row=512;
    for(i=0; i<num_particles; i++) if(particle_level[i] != FREE_PARTICLE_LEVEL)
    {
        particle_t[i] = tl[min_level];
        particle_dt[i] = 0.0;
    }
        
    for(level = min_level; level <= max_level; level++){
        cart_debug("updating level %u", level );
        update_buffer_level(level, all_hydro_vars, num_hydro_vars);
    }
        
//need zoom_DM for buildmesh:    build_mesh(); //defines refinement volume
    for ( i = 0; i < nDim; i++ ) {
        refinement_volume_min[i] = 0;
        refinement_volume_max[i] = (double)num_grid;
    }

    num_particle_species=1;
    num_particles_total=0;
    last_star_id=-1;
    float dm_star; double dx;
    int ic[NSTARS];
    ic[0] = icell_central(0,0,0);
    level = cell_level(ic[0]);
    dx=cell_size[level];
    if(NSTARS==64){
        for(k=0;k<4;k++){
            for(j=0;j<4;j++){
                for(i=0;i<4;i++){
                    ic[k*4*4+ j*4+ i] = icell_central(i*dx, j*dx, k*dx);
                }
            }
        }
    }
    if(NSTARS==27){
        for(k=0;k<3;k++){
            for(j=0;j<3;j++){
                for(i=0;i<3;i++){
                    ic[k*3*3+ j*3+ i] = icell_central(i*dx, j*dx, k*dx);
                }
            }
        }
    }
    if(NSTARS==8){
        ic[1] = icell_central(dx,  0, 0);
        ic[2] = icell_central(0,  dx, 0);
        ic[3] = icell_central(dx, dx, 0);
        ic[4] = icell_central(0,   0, dx);
        ic[5] = icell_central(dx,  0, dx);
        ic[6] = icell_central(0,  dx, dx);
        ic[7] = icell_central(dx, dx, dx);
    }
    if(NSTARS==7){
        ic[1] = icell_central( dx,0,0);
        ic[2] = icell_central(-dx,0,0);
        ic[3] = icell_central(0, dx,0);
        ic[4] = icell_central(0,-dx,0);
        ic[5] = icell_central(0, 0 , dx);
        ic[6] = icell_central(0, 0 ,-dx);
    }

/*     cart_debug(""); */
/*     cart_debug(""); */
/*     cart_debug("Checking direction scheme"); */
/*     cart_debug(""); */
/*     int ioct, ichild[num_children]; */
/*     int iplus_direction ; //-x=0,+x=1;+y=3;+z=5 */
/*     int mj3, ii; */
/*     double pos[nDim]; */
/*     icell = icell_central(dx,-dx,dx); */

/*     for(mj3=0;mj3<nDim;mj3++){ */
/*         iplus_direction = mj3*2+1; //-x=0,+x=1;+y=3;+z=5 */
/*         cart_debug("%d-direction",mj3); */
/*         ioct = cell_parent_oct(icell); */
/*         oct_all_children(ioct,ichild); */
/*         for ( j=0 ; j<num_children; j++ ){ */
/*             i = in_local_oct[j][iplus_direction]*(-2) + 1; */
/*             ii = ichild[j]; */
/*             if(ii==icell){cart_debug("DOHNE!");} */
/*             cell_center_position(ii,pos); */
/*             cart_debug("sign=%d,position=%f",i,pos[mj3]); */
/*         } */
/*     } */
/*     exit(1); */
    
    for(i=0;i<NSTARS;i++){
        cart_debug("ic: %d",ic[i]);
    }
    dm_star = mstar_one_msun*constants->Msun/units->mass/(NSTARS_1D*NSTARS_1D*1.0); //snl resolution scaling /NSTARS;
    for(i=0; i<NSTARS; i++){
        ic_star_spread(ic[i],dm_star);
    }
    
    /* make it impossible to produce more stars */
    for ( i = 0; i < nDim; i++ ) {
        star_formation_volume_min[i] = 0;
        star_formation_volume_max[i] = 0;
    }
        
    if ( !buffer_enabled ) {
        cart_debug("building cell buffer");
        build_cell_buffer();
        repair_neighbors();
    }

    hydro_magic( min_level );
    hydro_eos( min_level );
    
    cart_debug("done with initialization");
        
    check_map();

    
    tot_momentum0=0;
    icell = icell_central(0,0,0);
    for(idir=0;idir<nDim;idir++){
        adv_velocity[idir] = cell_momentum(icell,idir)/cell_gas_density(icell);
        cart_assert(adv_velocity[idir] == cell_momentum(0,idir)/cell_gas_density(0));
        cart_debug("advection dir %d velocity %e ",idir,adv_velocity[idir]*units->velocity/constants->kms);
    }
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL && CELL_TYPE_LEAF, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            icell=level_cells[i];
            tot_energy0 += cell_gas_energy(icell)*cell_volume[level];
            
            if(cell_is_leaf(icell)){
                for(idir=0;idir<nDim;idir++){
                    frame_momentum[idir]=
                        (cell_momentum(icell,idir) - adv_velocity[idir]*cell_gas_density(icell));
                    cart_assert(frame_momentum[idir]==0);
                }
                
                tot_momentum0 += sqrt( cell_momentum(icell,0)*cell_momentum(icell,0) +
                                       cell_momentum(icell,1)*cell_momentum(icell,1) +
                                       cell_momentum(icell,2)*cell_momentum(icell,2) )
                    *cell_volume[level];
/*                 tot_momentum0 += sqrt( */
/*                     frame_momentum[0]*frame_momentum[0] + */
/*                     frame_momentum[1]*frame_momentum[1] + */
/*                     frame_momentum[2]*frame_momentum[2] ) */
/*                     *cell_volume[level]; */
            }
        }
        cart_free( level_cells );
    }
    cart_debug("momentum %e",tot_momentum0);
    hydro_magic( min_level );
    hydro_eos( min_level );


#ifdef RADIATIVE_TRANSFER
#ifdef RT_DEBUG
  rt_debug.Mode = -1;
  rt_debug.Pos[0] = 0.5*num_grid - pow(0.5,1.0+max_level);
  rt_debug.Pos[1] = 0.5*num_grid - pow(0.5,1.0+max_level);
  rt_debug.Pos[2] = 0.5*num_grid - pow(0.5,1.0+max_level);
#endif
#endif


}


