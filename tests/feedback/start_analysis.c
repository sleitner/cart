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


#include "start_radp.h"


#include "../../src/run/gravity_step.h"
#include "../../src/run/hydro_step.h"
#include "../../src/run/hydro_tracer_step.h"
#include "../../src/run/particle_step.h"
#include "../../src/run/rt_step.h"
#include "../../src/run/starformation_step.h"
#include "../../src/run/starformation_feedback_step.h"
#include "../../src/run/step.h"

#include "../../src/extra/output_slice.h"

#define sqr(a) ((a)*(a))
void Out4IFrIT();

#ifdef RT_OTVET_SAVE_FLUX
extern int rt_flux_frequency;
extern float rt_flux[num_cells][num_neighbors]; // //-x=0,+x=1;+y=3;+z=5 
#endif

#ifdef VIEWDUMP
#include "../../src/extra/viewdump.h"
#endif
#include "../../src/extra/ifrit.h"

#define OUTLEVEL (max_level)
#ifndef STARFORM
#error need stars defined for feedback tests
#endif

extern double totmomentum0;
/* extern const int NSTARS; */

float cell_Urad(int cell);
float cell_tauUV(int cell);


/*
//  Plugin
*/
void record_momentum_input(int level, int cell, double dU);
void outslice();
void outslicebegin();
void outsliceCFL();

#include "../../src/extra/output_slice.h"

plugin_t outslicePlugin = {NULL};

const plugin_t* add_plugin(int id){
    if(id==0){
/*         outslicePlugin.AfterCFLRestart = outsliceCFL; */
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
float wgt_avg[num_bins];
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
            wgt_avg[bin] += r1;
                                                                                
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
            wgt_avg[bin] += r1;
        }
    }
}





 void outsliceCFL(){ 
     ic_refine_levels();
/*      dmom_from_press=0; */
/* /\*      dE_from_press = 0; *\/ */
} 
void outslicebegin(){
    write_restart( WRITE_SAVE, WRITE_SAVE, WRITE_SAVE, NULL );
//    outslice();  
}
double mom_from_press=0;
/* double E_from_press=0; */
void outslice() {

    int i;
    char filename[128];
    FILE *RADP;
    FILE *fp1;
    FILE *fp2;
    FILE *fp3, *fp4;
    float reduced_rho[num_bins];
    float reduced_pressure[num_bins];
    float reduced_vel[num_bins];
    float reduced_wgt_avg[num_bins];

    float *fluxr;
    float *fluxl;
    
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
//snl1               (tphys_from_tcode(tl[min_level])-tphys_from_auni(auni_init))*1e-6,
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
        wgt_avg[i] = 0.0;
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
    MPI_Reduce( wgt_avg, reduced_wgt_avg, num_bins, MPI_FLOAT, MPI_SUM, MASTER_NODE, mpi.comm.run );

    if ( local_proc_id == MASTER_NODE ) {
        sprintf(filename, "%s/%s_stp%04u.txt", output_directory, jobname, step );
        RADP = fopen(filename,"w");
//snl       fprintf(RADP, "#time=%e \n",tl[min_level]-t_init);
        fprintf(RADP, "#time=%e \n",tl[min_level]);
        for ( i = 0; i < num_bins/2; i++ ) {
            if ( reduced_wgt_avg[i] > 0.0 ) {
                reduced_vel[i] /= reduced_wgt_avg[i];
                reduced_pressure[i] /= reduced_wgt_avg[i];
                reduced_rho[i] /= reduced_wgt_avg[i];
                
                fprintf(RADP, "%e %e %e %e \n",
                        radii[i]*units->length/constants->pc,
                        reduced_rho[i],
                        reduced_pressure[i],
                        reduced_vel[i]
                    );
            }
        }
        fclose(RADP);


#ifdef RT_OTVET_SAVE_FLUX
#define num_wlen  5
        //6eV, 13.8eV, 24eV,120eV, 1.2KeV
        float wlen[num_wlen]={2000.0,1000.0,500.0,100.0,10.0};
        float ngxi[num_wlen];
        double Z_cell, dust_atten_factor, dust_atten_tauUV;
        double tau_tauMW, tau_tauUV;
            
        double pos_central[nDim]={num_grid/2.,num_grid/2.,num_grid/2.};
        int idir=0;
        int nline;
        double fact_hi_level=pow(2.0,max_level-min_level);
    
        sprintf(filename, "%s/%s_flux%04u.txt", output_directory, jobname, step );
        fp1 = fopen(filename,"w");
//snl        fprintf(fp1, "#time=%e \n",tl[min_level]-t_init );
        fprintf(fp1, "#time=%e \n",tl[min_level]);
        
        sprintf(filename, "%s/%s_inten%04u.txt", output_directory, jobname, step );
        fp2 = fopen(filename,"w");
        fprintf(fp2, "#time=%e \n",tl[min_level]);
        
        sprintf(filename, "%s/%s_attn%04u.txt", output_directory, jobname, step );
        fp3 = fopen(filename,"w");
        fprintf(fp3, "#time=%e \n",tl[min_level]);
        
        sprintf(filename, "%s/%s_kicks%04u.txt", output_directory, jobname, step );
        fp4 = fopen(filename,"w");
        fprintf(fp4, "#time=%e \n",tl[min_level]);
        

        cart_debug("ul=%e ut=%e ue=%e c=%e ul/ut=%e dx=%e",
                   units->length, units->time, units->energy_density,
                   constants->c, units->length/units->time, cell_size[max_level]
            );
        
#ifdef SNLUPDATE             //need cell_tauUV...
        nline=num_grid/cell_size[max_level];
        pos[1]=num_grid/2;
        pos[2]=num_grid/2;
/*         fluxr = cart_alloc(float, nline); */
/*         fluxl = cart_alloc(float, nline); */
/*             fluxl[i] = rt_flux[icell][2*idir]; */
/*             fluxr[i] = rt_flux[icell][2*idir+1]; */
/*         cart_erase(fluxl, fluxr); */
        double tauMW_sum=0;
        double tauUV_sum=0;
        double fiMW;
        double fiUV, fir2;
        double rabs;
        for(i=0;i<nline;i++){
            pos[idir] = (i-nline/2)/fact_hi_level + pos_central[idir];
            
            if(pos[idir]>=num_grid){pos[idir]-=num_grid;}
            if(pos[idir]<0        ){pos[idir]+=num_grid;}
            
            icell=cell_find_position_above_level(max_level,pos);
            level=cell_level(icell);
            cart_assert(icell!=-1);
            
            
//  UV field at 12.0eV in units of Draine field (1.0e6 phot/cm^2/s/ster/eV)
            //rtUmw(icell); // needs chemistry or LW_RATE
            //rtDmw(icell);      //cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell))
        
            //cell_gas_metal_density(cell)/cell_gas_density(cell)/constants->Zsun;
            Z_cell = rtDmw(icell); // dust or metals in units of zsun
            tau_tauMW = (cell_HI_density(icell)+2.0*cell_H2_density(icell))
                                 *units->number_density * 2.0e-21 * Z_cell
                                 *units->length*cell_sobolev_length(icell);
            tau_tauUV =   cell_tauUV(icell) ;
            cart_assert(tau_tauUV>0);
            
            dust_atten_factor = exp( -tau_tauMW );
            dust_atten_tauUV = exp( -tau_tauUV );
            
            
            rtGetRadiationField( icell, num_wlen, wlen,ngxi);

            if(i>nline/2){
                tauMW_sum +=tau_tauMW;
                tauUV_sum +=tau_tauUV;

                rabs=fabs((pos[idir]-pos_central[idir]) * units->length/constants->pc);
                
                fiMW=400/(rabs*rabs)*exp(-tauMW_sum);
                fiUV=400/(rabs*rabs)*exp(-tauUV_sum);
                fir2=400/rabs/rabs;
            }else{
                tauMW_sum=0;
                tauUV_sum=0;
                fiMW=0;
                fiUV=0;
                fir2=0;
            }
            fprintf(fp1,"%e   %e %e %e  %e %e  %e %e %e   %e %e %e\n",
                    (pos[idir]-pos_central[idir]) * units->length/constants->pc,
/*                     rt_flux[icell][2*idir] */
/*                      *units->energy_density*units->length/units->time,  */
/*                     rt_flux[icell][2*idir+1] */
/*                      *units->energy_density*units->length/units->time,  //ergs/cm^3*cm/s */
                    cell_radiation_flux(icell,2*idir)  
                      *units->energy_density*units->length/units->time,  
                    cell_radiation_flux(icell,2*idir+1)
                      *units->energy_density*units->length/units->time,  
                    cell_Urad(icell)
                     *units->energy_density*constants->c,  //ergs/cm^3*c[cm/s]
                    
                    tau_tauMW,tau_tauUV,
                    
                    fiMW,fiUV, fir2,
                    
                    cell_var(icell,RT_VAR_OT_FIELD)
		    *units->energy_density*constants->c,
                    cell_var(icell,RT_VAR_OT_FIELD+1)
		    *units->energy_density*constants->c,
                    cell_var(icell,RT_VAR_OT_FIELD+2)
		    *units->energy_density*constants->c
                );
                    
            fprintf(fp2,"%e   %e %e  %e %e  %e %e  %e %e  %e %e  \n",
                    (pos[idir]-pos_central[idir]) * units->length/constants->pc,
                    wlen[0],ngxi[0],
                    wlen[1],ngxi[1],
                    wlen[2],ngxi[2],
                    wlen[3],ngxi[3],
                    wlen[4],ngxi[4]
                );
                    
            fprintf(fp3,"%e   %e %e %e \n",
                    (pos[idir]-pos_central[idir]) * units->length/constants->pc,
                    dust_atten_factor,dust_atten_tauUV, 1-cell_HI_fraction(icell)
                );
            
            fprintf(fp4,"%e   %e \n",
                    (pos[idir]-pos_central[idir]) * units->length/constants->pc,
                    cell_radiation_pdot(icell,2,level)/cell_gas_density(icell)*units->velocity/constants->kms
                );
        }
        fclose(fp1);
        fclose(fp2);
        fclose(fp3);
        fclose(fp4);
#endif
    }
    cart_debug("length %e erg %e cc %e ",units->length,constants->erg, constants->cc);
#endif



        
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
#ifdef RADIATIVE_TRANSFER
#ifdef SNLUPDATE
        dump_plane(OUT_CELL_TAUUV, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); 
/*         dump_plane(OUT_CELL_RADIATION_PRESSURE, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);  */
        dump_plane(OUT_CELL_URAD, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); 
#endif
        dump_plane(OUT_CELL_FHI, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); 
        dump_plane(OUT_CELL_FH2, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output); 
#ifdef RT_OTVET_SAVE_FLUX
         dump_plane(OUT_CELL_FLUX0, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, output);  
#endif        
#endif
        cart_debug("\n\n");
        fclose(output);
    }else{
        cart_error("output not MPI_parallel!!");
    }

    //Out4IFrIT();
}



void run_output(){
    cart_debug("");   
}



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

