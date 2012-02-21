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

#include "extra/output_slice.h"

#include "start_radp.h"


#include "run/gravity_step.h"
#include "run/hydro_step.h"
#include "run/hydro_tracer_step.h"
#include "run/particle_step.h"
#include "run/rt_step.h"
#include "run/starformation_step.h"
#include "run/starformation_feedback_step.h"
#include "run/step.h"


void Out4IFrIT();


#ifdef VIEWDUMP
#include "extra/viewdump.h"
#endif
#include "extra/ifrit.h"

#define OUTLEVEL (max_level)
#ifndef STARFORM
#error need stars defined for feedback tests
#endif

extern double totmomentum0;
/* extern const int NSTARS; */


void record_momentum_input(int level, int cell, double dU);
void outslice();
void outslicebegin();
void outsliceCFL();
plugin_t outslicePlugin = {NULL};
const plugin_t* add_plugin(int id){
    if(id==0){
        outslicePlugin.AfterCFLRestart = outsliceCFL;
        outslicePlugin.RunBegin = outslicebegin;
        outslicePlugin.GlobalStepEnd = outslice;
//snl1
//        outslicePlugin.StarformationFeedbackEnd = record_momentum_input;
        return &outslicePlugin;
    }else{
        return NULL;
    }
}

/* double dmom_from_press =0; */
/* /\* double dE_from_press =0; *\/ */
/* static double told=-99; */
/* void record_momentum_input(int level, int icell, double dU){ */
/*     double dmom; */
/*     if(told != tl[level]){ */
/* /\*     double dE;  *\/ */
/* /\*     dE =  (dU - cell_gas_pressure(icell)/(cell_gas_gamma(icell)-1))* *\/ */
/* /\*         cell_volume[level]; //energy input *\/ */
/* #ifdef STAR_PRESSURE_IN_INTERNAL_ENERGY */
/*         dmom =   cell_gas_pressure(icell)  *6 *(cell_size[level]*cell_size[level])*dtl[level]; */
/* #else */
/* #ifdef  STAR_PRESSURE_FROM_PARTICLES */
/*         dmom = dU*(cell_gas_gamma(icell)-1)*6 *(cell_size[level]*cell_size[level])*dtl[level]; */
/* #else */
/* #error "Star pressure must come from somewhere" */
/* #endif */
/* #endif */
/*         if(NSTARS ==8){ */
/* //multiply pressure by new block size area; pressure */
/*             dmom = (2*2)*dmom;  */
/*         } */
/*         if(NSTARS ==27){ */
/*             dmom = 3*3*dmom;  */
/*         } */
/*         if(NSTARS ==64){ */
/*             dmom = 4*4*dmom;  */
/*         } */
/*         MPI_Reduce( MPI_IN_PLACE, &dmom, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run ); */
/* /\*     MPI_Reduce( MPI_IN_PLACE, &dE  , 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run ); *\/ */
/*         if ( local_proc_id == MASTER_NODE ) { */
/*             dmom_from_press +=dmom; */
/* /\*         dE_from_press +=dE; *\/ */
/*         } */
/*         told=tl[level]; */
/*     } */
        
    
/* } */

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////   OUTPUT    ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



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
    double momentum;
    float r1, r, v;
    float amt, cur_r;
    int bin;

    if ( cell_is_leaf( cell ) ) {
        cell_center_position(cell, pos);
        r = sqrt( (pos[0]-((float)num_grid/2.0))*(pos[0]-((float)num_grid/2.0))+
                  (pos[1]-((float)num_grid/2.0))*(pos[1]-((float)num_grid/2.0))+
                  (pos[2]-((float)num_grid/2.0))*(pos[2]-((float)num_grid/2.0)));

        /* this needs to be converted to v_r */
        momentum = sqrt( cell_momentum(cell,0)*cell_momentum(cell,0) +
                               cell_momentum(cell,1)*cell_momentum(cell,1) +
                               cell_momentum(cell,2)*cell_momentum(cell,2) ) ;
        v = momentum / cell_gas_density(cell);
        
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
 void outsliceCFL(){ 
     ic_refine_levels();
/*      dmom_from_press=0; */
/* /\*      dE_from_press = 0; *\/ */
} 
void outslicebegin(){
    write_restart( WRITE_SAVE, WRITE_SAVE, WRITE_SAVE, NULL );
    outslice(); 
}
double mom_from_press=0;
/* double E_from_press=0; */
void outslice() {
    int i;
    char filename[128];
    FILE *RADP;
    float reduced_rho[num_bins];
    float reduced_pressure[num_bins];
    float reduced_vel[num_bins];
    float reduced_avgs[num_bins];
    int level;
    int num_level_cells;
    int *level_cells;
    double pos[nDim], pabs[nDim], center_pos[nDim];

    FILE *output;
    int icell, j;
    const int nvars = 4;
    const int nbin1 = 128;
    int varid[] = { HVAR_PRESSURE, HVAR_GAS_DENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC };
    int nbin[] = { nbin1, nbin1, nbin1 };
    double bb[6];
    float frame_momentum[nDim];
    double some_momentum;
    double some_momentum_dir[nDim];
    double tot_momentum;
    double tot_momentum_dir[nDim];
    double tot_energy=0; 
    int icell_level;
    int icenter = icell_central(0,0,0);
    int icenter_level=cell_level(icenter);
    cart_debug("output feedback test:");

    
    cell_center_position(icenter,center_pos);
    center_pos[0]=particle_x[0][0];
    center_pos[1]=particle_x[0][1];
    center_pos[2]=particle_x[0][2];
        
    bb[0] = bb[2] = bb[4] = num_grid*(0.5-0.125);
    bb[1] = bb[3] = bb[5] = num_grid*(0.5+0.125);
    
    sprintf( filename, "%s/pabs_%04u.dat", output_directory, step );
    output = fopen( filename, "w" );
         
    /* now dump the total gas momentum */
    some_momentum=0;
    some_momentum_dir[0]=0;
    some_momentum_dir[1]=0;
    some_momentum_dir[2]=0;
    tot_momentum_dir[0]=0;
    tot_momentum_dir[1]=0;
    tot_momentum_dir[2]=0;
    tot_energy=0;
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LEAF && CELL_TYPE_LOCAL , &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            icell=level_cells[i];
            if(cell_is_leaf(icell)){
                for(j=0;j<nDim;j++){
                    frame_momentum[j]=
                        (cell_momentum(icell,j) - adv_velocity[j]*cell_gas_density(icell));
                    pabs[j]=fabs(frame_momentum[j]);
                }
                if(pabs[0] > 0 || pabs[1] > 0 || pabs[2] > 0){
/*                 if(pabs[0] > 1e5 || pabs[1] > 1e5 || pabs[2] > 1e5){ */
                    cell_center_position(icell,pos);
                    if( compute_distance_periodic(pos,center_pos)< 5){
                        fprintf(output,"%f %f %f %e %e %e %e\n",
                                pos[0],pos[1],pos[2],
                                frame_momentum[0],frame_momentum[1],frame_momentum[2],
                                pabs[0]+pabs[1]+pabs[2]);
                        some_momentum += 
                            (pabs[0]+pabs[1]+pabs[2])*cell_volume[level];
                        for(j=0;j<nDim;j++){
                            some_momentum_dir[j] += pabs[j]*cell_volume[level];
                        }
                    }
                }
/*                 some_momentum += sqrt( */
/*                     frame_momentum[0]*frame_momentum[0]+ */
/*                     frame_momentum[1]*frame_momentum[1]+ */
/*                     frame_momentum[2]*frame_momentum[2] */
/*                     )*cell_volume[level]; */
                tot_energy += cell_gas_energy(icell)
                    *cell_volume[level];
            }
        }
        cart_free( level_cells );
    }
    fclose(output);
    
/*     MPI_Reduce( &stellar_mass, &total_stellar_mass, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );        */
    MPI_Reduce( &some_momentum, &tot_momentum, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( &some_momentum_dir[0], &tot_momentum_dir[0], 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( &some_momentum_dir[1], &tot_momentum_dir[1], 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( &some_momentum_dir[2], &tot_momentum_dir[2], 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
    cart_debug("\n\ntot_momentum(t),(0)=%e %e",tot_momentum,tot_momentum0);
    cart_debug("momentum(center)=%e momentum(5)=%e advection_velocity=%e \n\n",cell_momentum(icell_central(0,0,0),0),cell_momentum(5,0),adv_velocity[0]);
    

    if ( local_proc_id == MASTER_NODE ) {
/*         mom_from_press += dmom_from_press; */
/*         dmom_from_press=0; */
/* /\*         E_from_press += dmom_from_press; *\/ */
/* /\*         dE_from_press=0; *\/ */
        sprintf(filename, "totmompress.dat" );
        RADP = fopen(filename,"a+");
        fprintf(RADP, "%e %e %e  %e %e %e %e  \n",
                (tphys_from_tcode(tl[min_level])-tphys_from_auni(auni_init))*1e-6,
                (tot_momentum)*units->velocity*units->mass/constants->kms/constants->Msun ,
                mom_from_press       *units->velocity*units->mass/constants->kms/constants->Msun,
                (tot_momentum_dir[0])*units->velocity*units->mass/constants->kms/constants->Msun ,
                (tot_momentum_dir[1])*units->velocity*units->mass/constants->kms/constants->Msun ,
                (tot_momentum_dir[2])*units->velocity*units->mass/constants->kms/constants->Msun ,
                mom_from_press/3.0   *units->velocity*units->mass/constants->kms/constants->Msun 
/*                  (tot_energy-tot_energy0)*units->energy/constants->erg  */
/*                 E_from_press*units->energy/constants->erg */
                );
                
        fclose(RADP);
    }

    cart_debug("--=physical conditions(cent cell=%d)=--",icell);
    cart_debug("IC:rho=%e[#/cc] %e[g/cc] %e[code]", n_h2, cell_gas_density(icell)*units->density/constants->gpercc, cell_gas_density(icell));
    cart_debug("rho=%e[#/cc] %e[g/cc] %e[code]", cell_gas_density(icell)*units->number_density/constants->cc, cell_gas_density(icell)*units->density/constants->gpercc, cell_gas_density(icell));
    cart_debug("T=%e[K] %e[code]", cell_gas_temperature(icell)*units->temperature,cell_gas_temperature(icell) );
    cart_debug("P=%e[ergs/cc] %e[code]", cell_gas_pressure(icell)*units->energy_density/constants->barye, cell_gas_pressure(icell));
    cart_debug("mstar=%e[Msun] %e[code] ",mstar_one_msun,mstar_one_msun/units->mass*constants->Msun);
    cart_debug("cell_max[pc]=%e cell_min[pc]%e ",cell_size[min_level]*units->length/constants->pc,cell_size[cell_level(icell)]*units->length/constants->pc);
    
    /* now dump the radial profiles */
    for ( i = 0; i < num_bins; i++ ) {
        radii[i] = ((float)i + 0.5) * bin_width;
        vel[i] = 0.0;
        pressure[i] = 0.0;
        rho[i] = 0.0;
        avgs[i] = 0.0;
        sedov_analytic[i] = 0.0;
    }
    for ( level = min_level; level <= max_level; level++ ) {
        select_level( level, CELL_TYPE_LOCAL, &num_level_cells, &level_cells );
        for ( i = 0; i < num_level_cells; i++ ) {
            radial_average( level_cells[i], level );
        }
        cart_free( level_cells );
    }

    MPI_Reduce( rho, reduced_rho, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( pressure, reduced_pressure, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( vel, reduced_vel, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( avgs, reduced_avgs, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, mpi.comm.run );

    if ( local_proc_id == MASTER_NODE ) {
        sprintf(filename, "%s/%s_stp%04u.txt", output_directory, jobname, step );
        RADP = fopen(filename,"w");
        fprintf(RADP, "#time=%e \n",tl[min_level]-t_init);
        for ( i = 0; i < num_bins/2; i++ ) {
            if ( reduced_avgs[i] > 0.0 ) {
                reduced_vel[i] /= reduced_avgs[i];
                reduced_pressure[i] /= reduced_avgs[i];
                reduced_rho[i] /= reduced_avgs[i];
                
                fprintf(RADP, "%e %e %e %e \n",
                        radii[i]*units->length/constants->pc,
                        reduced_rho[i], reduced_pressure[i],
                        reduced_vel[i] );
//                , sedov_analytic[i]
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
        double slice_region_hsize=slice_hsize_pc; 
        slice_region_hsize *= constants->pc/units->length;
        cart_debug("original ouput size=%e [pc]",slice_region_hsize*units->length/constants->pc);
        int nsgrid_half=(int)slice_region_hsize/cell_size[OUTLEVEL];
        slice_region_hsize=nsgrid_half*cell_size[OUTLEVEL]+cell_size[OUTLEVEL]/2.;
        cart_debug("adjusted ouput size=%e [code] %d[cells]",2*slice_region_hsize,2*nsgrid_half+1);
        cart_debug("adjusted ouput size=%e [pc]",2*slice_region_hsize*units->length/constants->pc);
        
        double pos[3];
        icell=icell_central(0,0,0);
        cell_center_position(icell,pos);
        /* for(i=0;i<nDim; i++){ */
        /*     //now everything in a given cell will be assigned to that cell by pos/dx */
        /*     // |00x00|1111| */
        /*     pos[i]-=cell_size[OUTLEVEL]/2.0; */
        /* } */
        
        sprintf( filename, "%s/%s_slice_%04u.dat", output_directory, jobname, step );
        output = fopen( filename, "w" );
        cart_debug("dumping plane");
        cart_debug("in output_slice; axes (xyz)=%d%d%d nsgrid=%d out level=%d" ,
                   axis_direction[slice_axis_z][0],
                   axis_direction[slice_axis_z][1],slice_axis_z, 2*nsgrid_half, OUTLEVEL);
        
        dump_plane(OUT_CELL_DENSITY, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_INTERNAL_ENERGY, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_MOMENTUM+axis_direction[slice_axis_z][0], OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_MOMENTUM+axis_direction[slice_axis_z][1], OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_MOMENTUM+slice_axis_z, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
/*         dump_plane(OUT_CELL_SOUNDSPEED, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); */
/*         dump_plane(OUT_CELL_MACH, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); */
        dump_plane(OUT_CELL_LEVEL, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
        dump_plane(OUT_CELL_PRESSURE, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);
/*         dump_plane(OUT_CELL_TAUUV, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); */
/*         dump_plane(OUT_CELL_RADIATION_PRESSURE, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); */
/*         dump_plane(OUT_CELL_URAD, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); */
        cart_debug("\n\n");
        fclose(output);
    }else{
        cart_error("output not MPI_parallel!!");
    }

    //Out4IFrIT();

}



#endif /*  USER_PLUGIN */

void run_output(){
}



#include "extra/ifrit.h"


void Out4IFrIT()
{
  const int nbin1 = 256;
  int ids[] = { I_HI_FRACTION, I_GAS_OVERDENSITY, I_GAS_TEMPERATURE, I_CELL_LEVEL, I_LOCAL_PROC, rt_field_offset, rt_field_offset+rt_num_freqs };
  int nbin[] = { nbin1, nbin1, nbin1 };
  double pos[] = { 0.5*num_grid, 0.5*num_grid, 0.5*num_grid };
  char fname[999];

  sprintf(fname,"%s/ifrit-box-a=%05d.bin",output_directory,step);

  ifrit.OutputMesh(fname,-1,nbin,pos,sizeof(ids)/sizeof(int),ids);
}

