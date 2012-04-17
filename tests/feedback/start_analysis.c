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
#include "frt/frt_c.h"
#include "rt.h"

#include "start_radp.h"

#include "run/gravity_step.h"
#include "run/hydro_step.h"
#include "run/hydro_tracer_step.h"
#include "run/particle_step.h"
#include "run/rt_step.h"
#include "run/starformation_step.h"
#include "run/starformation_feedback_step.h"
#include "run/step.h"

#include "extra/output_slice.h"

#define sqr(a) ((a)*(a))
void Out4IFrIT();

#ifdef RT_OTVET_SAVE_FLUX
extern int rt_flux_frequency;
extern float rt_flux[num_cells][num_neighbors]; /* -x=0,+x=1;+y=3;+z=5 */
#endif

#ifdef VIEWDUMP
#include "extra/viewdump.h"
#endif
#include "extra/ifrit.h"

#define OUTLEVEL (max_level)
#ifndef STARFORM
#error need stars defined for feedback tests
#endif

extern double p_tot0;
double p_input=0;

/*
//  Plugin
*/
void analysis();
void check_momentum();
void check_radiation();
void writebegin();
void record_momentum_input(int level, int cell, float pressure);
void radial_profiles();
void outslice();
void startCFL();

#include "extra/output_slice.h"

plugin_t outslicePlugin = {NULL};

const plugin_t* add_plugin(int id){
    if(id==0){
	/* If hitting CFL use writebegin to choose shorter timestep */
        outslicePlugin.RunBegin = writebegin;
        outslicePlugin.GlobalStepEnd = analysis;
/* 	outslicePlugin.AfterCFLRestart = startCFL;  */ /* not currently implimented */
/*      outslicePlugin.RadiationFeedbackEnd = record_momentum_input; *//* not currently implimented */
        return &outslicePlugin;
    }else{
        return NULL;
    }
}

double dp_input =0;
void startCFL(){ 
/*   ic_refine_levels(); */
     dp_input=0; 
} 
void writebegin(){
    write_restart( WRITE_SAVE, WRITE_SAVE, WRITE_SAVE, NULL ); 
    /*  analysis();  */
}

void analysis() {
    int icenter = icell_central(0,0,0);
    cart_debug("--=physical conditions(cent cell=%d)=--",icenter);
    cart_debug("IC:rho=%e[#/cc] %e[g/cc] %e[code]", n_h2, cell_gas_density(icenter)*units->density/constants->gpercc, cell_gas_density(icenter));
    cart_debug("rho=%e[#/cc] %e[g/cc] %e[code]", cell_gas_density(icenter)*units->number_density/constants->cc, cell_gas_density(icenter)*units->density/constants->gpercc, cell_gas_density(icenter));
    cart_debug("T=%e[K] %e[code]", cell_gas_temperature(icenter)*units->temperature,cell_gas_temperature(icenter) );
    cart_debug("P=%e[ergs/cc] %e[code]", cell_gas_pressure(icenter)*units->energy_density/constants->barye, cell_gas_pressure(icenter));
    cart_debug("mstar=%e[Msun] %e[code] ",mstar_one_msun,mstar_one_msun/units->mass*constants->Msun);
    cart_debug("cell_max[pc]=%e cell_min[pc]%e ",cell_size[min_level]*units->length/constants->pc,cell_size[cell_level(icenter)]*units->length/constants->pc);
    
#ifdef RT_OTVET_SAVE_FLUX
    cart_debug("check_radiation");
    check_radiation();
#endif
    cart_debug("check_momentum");
    check_momentum();
    cart_debug("radial_profiles");
    radial_profiles();
    if ( local_proc_id == MASTER_NODE ) {
    cart_debug("outslice");
	outslice();
    }
}

/* check on radiation-related quantities*/
#ifdef RT_OTVET_SAVE_FLUX
void rtGetRadiationFlux(int cell, float flux[num_neighbors]);
#define num_wlen  2
#define mabs(c) (c < 0 ? -c : c)
void check_radiation(){
    FILE *fp1;
    char filename[128];

    int level;
    int num_level_cells;
    int *level_cells;
    double pos[nDim], pabs[nDim], center_pos[nDim];

    float wlen[num_wlen]={911.0,1000.0}; /* 13.6eV = HI, 1000A = rtUV*/
    float ngxi[num_wlen];
    double Z_cell, dust_atten_factor, dust_atten_tauUV;
    double tau_tauMW, tauMW_sum, fiMW, rabs;

    int i,icell,nline, idir=0;
    double pos_central[nDim]={num_grid/2.,num_grid/2.,num_grid/2.};
    double fact_hi_level=pow(2.0,max_level-min_level);
    
    
    float fluxneighbors[num_neighbors];
    float fluxat1000;
    float ng911,ng1000, unitgetrf;


 
    sprintf(filename, "%s/%s_flux%04u.txt", output_directory, jobname, step );
    fp1 = fopen(filename,"w");
    
    tauMW_sum=0;
    nline=num_grid/cell_size[max_level];
    pos[1]=particle_x[0][1];
    pos[2]=particle_x[0][2];
    for(i=0;i<nline;i++){
	pos[idir] = (i-nline/2)/fact_hi_level + pos_central[idir];
        
	if(pos[idir]>=num_grid){pos[idir]-=num_grid;}
	if(pos[idir]<0        ){pos[idir]+=num_grid;}
        
	icell=cell_find_position_above_level(max_level,pos);
	level=cell_level(icell);
	cart_assert(icell!=-1);
        
	Z_cell = rtDmw(icell); 
	tau_tauMW = (cell_HI_density(icell)+2.0*cell_H2_density(icell))
	    *units->number_density * 2.0e-21 * Z_cell
	    *units->length*cell_sobolev_length(icell);
	
	if(i>nline/2){
	    tauMW_sum +=tau_tauMW;
	    rabs=fabs((pos[idir]-pos_central[idir]) * units->length/constants->pc);
	    fiMW=400/(rabs*rabs)*exp(-tauMW_sum);
	}else{
	    tauMW_sum=0;
	    fiMW=0;
	}

	rtGetRadiationFlux(icell,fluxneighbors);
	fluxat1000= 0.5*(mabs(fluxneighbors[2*idir])+mabs(fluxneighbors[2*idir+1]));
	rtGetRadiationField( icell, num_wlen, wlen,ngxi); // photon number density per Hz at wlen
        ng911=ngxi[0];
        ng1000=ngxi[1];
	
	/*convert ngxi to flux in [ergs/cm^2/s/Hz]*/
	unitgetrf  = 6.626e-27*constants->c;
        ng911          *=  unitgetrf;
        ng1000         *=  unitgetrf;

	fprintf(fp1,"%e %e %e %e \n",
		(pos[idir]-pos_central[idir]) * units->length/constants->pc,
		fluxat1000, ng1000, ng911
                );
    }
    fclose(fp1);
}
#endif




static double told=-99;
void record_momentum_input(int level, int icell, float pressure ){
    double dmom;
    if(told != tl[level]){
	dmom =   pressure *6 *(cell_size[level]*cell_size[level])*dtl[level]; 
/*         dmom = dU*(cell_gas_gamma(icell)-1)*6 *(cell_size[level]*cell_size[level])*dtl[level]; */
        if(NSTARS ==8){ /* multiply pressure by new block size area; pressure */
            dmom = (2*2)*dmom;
        }
        if(NSTARS ==27){
            dmom = 3*3*dmom;
        }
        if(NSTARS ==64){
            dmom = 4*4*dmom;
        }
        MPI_Reduce( MPI_IN_PLACE, &dmom, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
        if ( local_proc_id == MASTER_NODE ) {
            dp_input +=dmom;
        }
        told=tl[level];
    }
        
    
}


/* radial binning for analysis */
#define  num_bins   (num_grid*(1<<(max_level-min_level)))
#define  bin_width  ((float)num_grid/((float)num_bins))
float radii[num_bins];
float vel[num_bins];
float pressure[num_bins];
float rho[num_bins];
float wgt_avg[num_bins];
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





void check_momentum(){
    int i;
    char filename[128];
    FILE *fp1;
    int level;
    int num_level_cells;
    int *level_cells;
    double pos[nDim], pabs[nDim], center_pos[nDim];

    FILE *output;
    int icell, j;
    const int nvars = 4;
    float frame_momentum[nDim];
    double some_momentum;
    double some_momentum_dir[nDim];
    double p_tot;
    double p_tot_dir[nDim];
    int icell_level;
    int icenter = icell_central(0,0,0);
    int icenter_level=cell_level(icenter);
    
    cell_center_position(icenter,center_pos);
    center_pos[0]=particle_x[0][0];
    center_pos[1]=particle_x[0][1];
    center_pos[2]=particle_x[0][2];
        
    sprintf( filename, "%s/pabs_%04u.dat", output_directory, step );
    fp1 = fopen( filename, "w" );
         
    /* now dump the total gas momentum */
    some_momentum=0;
    some_momentum_dir[0]=0;
    some_momentum_dir[1]=0;
    some_momentum_dir[2]=0;
    p_tot_dir[0]=0;
    p_tot_dir[1]=0;
    p_tot_dir[2]=0;
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
                    cell_center_position(icell,pos);
                    if( compute_distance_periodic(pos,center_pos)< 5){
                        fprintf(fp1,"%f %f %f %e %e %e %e\n",
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
            }
        }
        cart_free( level_cells );
    }
    fclose(fp1);
    
/*     MPI_Reduce( &stellar_mass, &total_stellar_mass, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );        */
    MPI_Reduce( &some_momentum, &p_tot, 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( &some_momentum_dir[0], &p_tot_dir[0], 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( &some_momentum_dir[1], &p_tot_dir[1], 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
    MPI_Reduce( &some_momentum_dir[2], &p_tot_dir[2], 1, MPI_DOUBLE, MPI_SUM, MASTER_NODE, mpi.comm.run );
    cart_debug("\n\np_tot(t),(0)=%e %e",p_tot,p_tot0);
    cart_debug("momentum(center)=%e momentum(5)=%e advection_velocity=%e \n\n",cell_momentum(icell_central(0,0,0),0),cell_momentum(5,0),adv_velocity[0]);
    

    if ( local_proc_id == MASTER_NODE ) {
         p_input += dp_input; 
         dp_input=0; 
	 sprintf(filename, "momentum_output.dat" );
	 fp1 = fopen(filename,"a+");
	 fprintf(fp1, "%e %e %e  %e %e %e %e  \n",
		 (tphys_from_tcode(tl[min_level]))*1e-6,
		 (p_tot)*units->velocity*units->mass/constants->kms/constants->Msun ,
		 p_input       *units->velocity*units->mass/constants->kms/constants->Msun,
		 (p_tot_dir[0])*units->velocity*units->mass/constants->kms/constants->Msun ,
		 (p_tot_dir[1])*units->velocity*units->mass/constants->kms/constants->Msun ,
		 (p_tot_dir[2])*units->velocity*units->mass/constants->kms/constants->Msun ,
		 p_input/3.0   *units->velocity*units->mass/constants->kms/constants->Msun 
	     );
	 
	 fclose(fp1);
    }
}
void radial_profiles()
{
    float reduced_rho[num_bins];
    float reduced_pressure[num_bins];
    float reduced_vel[num_bins];
    float reduced_wgt_avg[num_bins];
    int i;
    char filename[128];
    FILE *fp1;
    int level;
    int num_level_cells;
    int *level_cells;

    /* now dump the radial profiles */
    for ( i = 0; i < num_bins; i++ ) {
        radii[i] = ((float)i + 0.5) * bin_width;
        vel[i] = 0.0;
        pressure[i] = 0.0;
        rho[i] = 0.0;
        wgt_avg[i] = 0.0;
        
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
        fp1 = fopen(filename,"w");
        fprintf(fp1, "#time=%e \n",tl[min_level]);
        for ( i = 0; i < num_bins/2; i++ ) {
            if ( reduced_wgt_avg[i] > 0.0 ) {
                reduced_vel[i] /= reduced_wgt_avg[i];
                reduced_pressure[i] /= reduced_wgt_avg[i];
                reduced_rho[i] /= reduced_wgt_avg[i];
                
                fprintf(fp1, "%e %e %e %e \n",
                        radii[i]*units->length/constants->pc,
                        reduced_rho[i],
                        reduced_pressure[i],
                        reduced_vel[i]
                    );
            }
        }
        fclose(fp1);
    }
}



#define slice_axis_z    2  
#define slice_hsize_pc  2000  //(->2000)

void outslice(){
    char filename[128];
    FILE *fp1;
    
    double slice_region_hsize;
    int icell,nsgrid_half=(int)slice_region_hsize/cell_size[OUTLEVEL];
    double pos[3];
    
#ifdef VIEWDUMP
    sprintf(filename, "%s/%s_%04u.v", output_directory, jobname, step );
    viewdump( filename, max_level, ((float)(num_grid)/2.0)+0.5*cell_size[max_level], 2, DUMP_HVARS, CELL_TYPE_LOCAL );
#endif

#ifdef HYDRO_TRACERS
    sprintf( filename, "%s/tracers_%04u.dat", output_directory, step );
    write_hydro_tracers( filename );
#endif /* HYDRO_TRACERS */

    /* output a 2-d slice through the center of the box */
    slice_region_hsize = slice_hsize_pc*constants->pc/units->length;
    slice_region_hsize = min(num_grid/2,slice_region_hsize);
    cart_debug("original ouput size=%e [pc]",slice_region_hsize*units->length/constants->pc);
    slice_region_hsize=nsgrid_half*cell_size[OUTLEVEL]+cell_size[OUTLEVEL]/2.;
    cart_debug("adjusted ouput size=%e [code] %d[cells]",2*slice_region_hsize,2*nsgrid_half+1);
    cart_debug("adjusted ouput size=%e [pc]",2*slice_region_hsize*units->length/constants->pc);
    
    icell=icell_central(0,0,0);
    cell_center_position(icell,pos);
    cart_debug("%e %e %e",pos[0],pos[1],pos[2]);
    
    sprintf( filename, "%s/%s_slice_%04u.dat", output_directory, jobname, step );
    fp1 = fopen( filename, "w" );
    cart_debug("dumping plane");
    cart_debug("in output_slice; axes (xyz)=%d%d%d nsgrid=%d out level=%d" ,
	       axis_direction[slice_axis_z][0],
	       axis_direction[slice_axis_z][1],slice_axis_z, 2*nsgrid_half, OUTLEVEL);
    
    dump_plane(OUT_CELL_DENSITY, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);
    dump_plane(OUT_CELL_INTERNAL_ENERGY, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);
    dump_plane(OUT_CELL_MOMENTUM+axis_direction[slice_axis_z][0], OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);
    dump_plane(OUT_CELL_MOMENTUM+axis_direction[slice_axis_z][1], OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);
    dump_plane(OUT_CELL_MOMENTUM+slice_axis_z, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);
/*         dump_plane(OUT_CELL_SOUNDSPEED, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1); */
/*         dump_plane(OUT_CELL_MACH, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1); */
    dump_plane(OUT_CELL_LEVEL, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);
    dump_plane(OUT_CELL_PRESSURE, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);
#ifdef RADIATIVE_TRANSFER
/*         dump_plane(OUT_CELL_TAUUV, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);  */
/*         dump_plane(OUT_CELL_RADIATION_PRESSURE, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);  */
/*         dump_plane(OUT_CELL_URAD, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);  */
    dump_plane(OUT_CELL_FHI, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1); 
    dump_plane(OUT_CELL_FH2, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1); 
#ifdef RT_OTVET_SAVE_FLUX
/*          dump_plane(OUT_CELL_FLUX0, OUTLEVEL, slice_axis_z, pos, slice_region_hsize, fp1);   */
#endif        
#endif
    cart_debug("\n\n");
    fclose(fp1);
	
    /* Out4IFrIT(); */
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

