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
void create_star_particle( int icell, float mass, double dt, int type );
    

double pos_central[nDim]={num_grid/2.,num_grid/2.,num_grid/2.};
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
    cell_momentum(icell,0) = advection_momentum(icell,0);
    cell_momentum(icell,1) = advection_momentum(icell,1);
    cell_momentum(icell,2) = advection_momentum(icell,2);
    cell_gas_gamma(icell) = constants->gamma;
    
    cell_gas_internal_energy(icell) =  convert_temp_to_ie(T_h2,icell);
    cell_gas_pressure(icell) =  cell_gas_internal_energy(icell) * (cell_gas_gamma(icell)-1.0);
    cell_gas_energy(icell) = cell_gas_internal_energy(icell)+cell_gas_kinetic_energy(icell);
    cart_assert( cell_gas_internal_energy(icell) > 0);

#ifdef ENRICH 
    cell_gas_metal_density_II(icell) = constants->Zsun*cell_gas_density(icell);
#ifdef ENRICH_SNIa
    cell_gas_metal_density_Ia(icell) = constants->Zsun*cell_gas_density(icell);
#endif
#endif
#ifdef RADIATIVE_TRANSFER
    if(T_ambient < T_h2){
	cell_HI_density(icell) = 0;
	cell_HII_density(icell) = cell_gas_density(icell)*constants->XH;
	cell_HeI_density(icell) = cell_gas_density(icell)*constants->XHe;
	cell_HeII_density(icell) = 0;
	cell_HeIII_density(icell) = 0;
	cell_H2_density(icell) = 0;
    }else{
	cell_HI_density(icell) = cell_gas_density(icell)*constants->XH;
	cell_HII_density(icell) = 0;
	cell_HeI_density(icell) = cell_gas_density(icell)*constants->XHe;
	cell_HeII_density(icell) = 0;
	cell_HeIII_density(icell) = 0;
	cell_H2_density(icell) = 0;
    }

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
    cart_assert( cell_gas_internal_energy(icell) > 0);

#ifdef ENRICH 
    cell_gas_metal_density_II(icell) = constants->Zsun*cell_gas_density(icell);
#ifdef ENRICH_SNIa
    cell_gas_metal_density_Ia(icell) = constants->Zsun*cell_gas_density(icell);
#endif
#endif
#ifdef RADIATIVE_TRANSFER
    if(T_ambient < T_h2){
	cell_HI_density(icell) = cell_gas_density(icell)*constants->XH;
	cell_HII_density(icell) = 0;
	cell_HeI_density(icell) = cell_gas_density(icell)*constants->XHe;
	cell_HeII_density(icell) = 0;
	cell_HeIII_density(icell) = 0;
	cell_H2_density(icell) = 0;
    }else{
	cell_HI_density(icell) = 0;
	cell_HII_density(icell) = cell_gas_density(icell)*constants->XH;
	cell_HeI_density(icell) = cell_gas_density(icell)*constants->XHe;
	cell_HeII_density(icell) = 0;
	cell_HeIII_density(icell) = 0;
	cell_H2_density(icell) = 0;
    }
      

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
        hydro_split_update(level); /* update all non-leafs with their child's value */
        hydro_eos(level);
    }
}

void ic_star_spread(int icell, float dm_star){
    double density_fraction;
    double new_density;
    int i,level, j, k;
    level = cell_level(icell);
    int neighbors[num_neighbors];
    int nneighbors[num_neighbors];
    int nnneighbors[num_neighbors];
    
#ifdef LOG_STAR_CREATION
    log_star_creation( icell, dm_star, FILE_RECORD);
#endif

    /* star particle can form and takes properties from gas */
    star_cell_conditions( icell );
    /* particle_v[ipart][i] = cell_momentum(icell,i) / cell_gas_density(icell); */
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
        create_star_particle(icell, dm_star, (double)dtl[level], STAR_TYPE_NORMAL);
    }
    remap_star_ids();

    /* set cells neighboring star to have the star-cell properties*/
    cell_all_neighbors( icell, neighbors );
    star_cell_conditions( icell );
    for(i=0;i<num_neighbors;i++){
        star_cell_conditions( neighbors[i] );
        
/* ////////////// neighbors of neighbors (or just choose a region) */
/*         cell_all_neighbors( neighbors[i], nneighbors ); */
/*         for(j=0;j<num_neighbors;j++){ */
/*             star_cell_conditions( nneighbors[j] ); */
            
/*             cell_all_neighbors( nneighbors[i], nnneighbors ); */
/*             for(k=0;k<num_neighbors;k++){ */
/*                 star_cell_conditions( nnneighbors[j] ); */
                
/*             } */
/*         } */
                
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
  cart_debug("refining");
    for(level=min_level; level<max_level; level++){
	cart_debug("%d level",level);
        modify( level, OP_REFINE );
	cart_debug("%d level",level);
        if((!refinement_indicator[SPATIAL_INDICATOR].use[level])){
            cart_error("need spatial refinement");
        }
    }
    cart_debug("done refining");
    build_cell_buffer();
    repair_neighbors();
    load_balance();
}


/* 
//init_run 
*/
void init_run() {
    int i,j,k;
    int level;
    
    int num_level_cells;
    int *level_cells;
    int icell,idir;
    float frame_momentum[nDim];
    double pos[nDim];
    double rstargas;

    
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
    auni[min_level] = a0;
    abox[min_level] = a0;
    abox_old[min_level] = a0;
    tl[min_level]=tcode_from_auni(a0);
    tl_old[min_level]=tl[min_level];
        
    units_set_art(omm0,hubble,box_size);
#else
    units_set(1.0,1.0,1.0); /* mass time length */
#endif
    units_init();
    units_update( min_level );
        
    /* 
    //set hydro 
    */
    cart_debug("start building cell buffers ");
    build_cell_buffer();
    cart_debug("built cell buffer");
    repair_neighbors();
    cart_debug("repaired neighbors");
    ic_refine_levels(); 
    cart_debug("refined levels ");
        
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
    cosmology_set_fixed();
        
    cart_debug("tl[min_level] = %f", tl[min_level] );
    cart_debug("au[min_level] = %f", auni[min_level] );
    cart_debug("ab[min_level] = %f", abox[min_level] );
    cart_debug("abo[min_level] = %f", abox_old[min_level] );
    cart_debug("DC mode = %f", cosmology->DeltaDC );
#endif /* COSMOLOGY */    


    for(i=0; i<num_particles; i++) if(particle_level[i] != FREE_PARTICLE_LEVEL)
    {
        particle_t[i] = tl[min_level];
        particle_dt[i] = 0.0;
    }
        
    for(level = min_level; level <= max_level; level++){
        cart_debug("updating level %u", level );
        update_buffer_level(level, all_hydro_vars, num_hydro_vars);
    }
    /* build_mesh(); need zoomed dark matter for buildmesh -- defines refinement volume */

    /* position stars */
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

    for(i=0;i<NSTARS;i++){
        cart_debug("ic: %d",ic[i]);
    }
    
    /* resolution scaling /NSTARS: */
    dm_star = mstar_one_msun*constants->Msun/units->mass/(NSTARS_1D*NSTARS_1D*1.0); 
    for(i=0; i<NSTARS; i++){
        ic_star_spread(ic[i],dm_star);
    }
    
    /* make it impossible to produce more stars */
    for ( i = 0; i < nDim; i++ ) {
        star_formation_volume_min[i] = 0;
        star_formation_volume_max[i] = 0;
    }

    rstargas = Radius_stargas*constants->pc/units->length;
    /* establish stellar conditions around the central region*/
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            icell=level_cells[i];
            cell_center_position( icell, pos );
            
            if(compute_distance_periodic(pos,pos_central) < rstargas ){
                star_cell_conditions( icell );
            }
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


    /* record initial momentum */
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
        cart_free( level_cells );
    }
    cart_debug("momentum %e",tot_momentum0);
    

#ifdef RADIATIVE_TRANSFER
#ifdef RT_DEBUG
  rt_debug.Mode = -1;
  rt_debug.Pos[0] = 0.5;
  rt_debug.Pos[1] = 0.5;
  rt_debug.Pos[2] = 0.5;
#endif
#endif


}
