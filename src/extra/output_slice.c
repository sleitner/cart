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
#include "times.h"
#include "refinement.h"
#include "refinement_indicators.h"
#include "refinement_operations.h"
#include "timing.h"
#include "units.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "gravity.h"
#include "density.h"
#include "io.h"
#include "auxiliary.h"
#include "starformation.h"

#ifdef RT_OTVET_SAVE_FLUX
extern int rt_flux_frequency;
extern float rt_flux[num_cells][num_neighbors]; // //-x=0,+x=1;+y=3;+z=5 
#endif

#include "output_slice.h"
const int string_size=256;

const int axis_direction[nDim][nDim-1] = {
#if nDim == 2
    {1},{0}
#elif nDim == 3
    {1,2},{2,0},{0,1}
#else
    #error  "Unsupported number of dimensions for other_axis"
#endif
};

void name_cell_property(int iflag, char varname[]){
    switch(iflag){
    case OUT_CELL_DENSITY:
        sprintf(varname,"density_numbercc\n");
        break;
    case OUT_CELL_INTERNAL_ENERGY:
        sprintf(varname,"temperature_kelvin\n");
        break;
    case OUT_CELL_MOMENTUM+0:
        sprintf(varname,"vx_kms\n");
        break;
    case OUT_CELL_MOMENTUM+1: 
        sprintf(varname,"vy_kms\n");
        break;
    case OUT_CELL_MOMENTUM+2: 
        sprintf(varname,"vz_kms\n");
        break;
    case OUT_CELL_SOUNDSPEED: 
        sprintf(varname,"cs_kms\n");
        break;
    case OUT_CELL_MACH: 
        sprintf(varname,"mach_number\n");
        break;
    case  OUT_CELL_PRESSURE:
        sprintf(varname,"pressure_ergcc\n");
        break;
        
#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD==RT_METHOD_OTVET)
    case  OUT_CELL_TAUUV: 
        sprintf(varname,"tauUV_number\n");
        break;
    case  OUT_CELL_RADIATION_PRESSURE:
        sprintf(varname,"radPoverP_number\n");
        break;
    case  OUT_CELL_URAD:
        sprintf(varname,"Urad_ergcc\n");
        break;
    case  OUT_CELL_FHI:
        sprintf(varname,"fh1_number\n");
        break;
    case  OUT_CELL_FH2:
        sprintf(varname,"fh2_number\n");
        break;
#ifdef RT_OTVET_SAVE_FLUX
    case  OUT_CELL_FLUX0:
        sprintf(varname,"flux0_ergcm2\n");
        break;
#endif        
#endif
    case OUT_CELL_LEVEL: 
        sprintf(varname,"level_number\n");
        break;
    default:
        cart_error("bad flag");
        break;
    }
}

float value_cell_property(int icell,int iflag){
    double px2, py2, pz2, vmag,cs;
    switch(iflag){
    case OUT_CELL_DENSITY: 
        return (float)(cell_gas_density(icell)
                       *units->number_density*constants->cc);
        break;
    case OUT_CELL_INTERNAL_ENERGY:
        return (float)(cell_gas_temperature(icell)
                       *units->temperature/constants->K);
        break;
    case OUT_CELL_MOMENTUM+0:
        return (float)(cell_momentum(icell,0)/cell_gas_density(icell)
                       *units->velocity/constants->kms);
        break;
    case OUT_CELL_MOMENTUM+1:
        return (float)(cell_momentum(icell,1)/cell_gas_density(icell)
                       *units->velocity/constants->kms);
        break;
    case OUT_CELL_MOMENTUM+2:
        return (float)(cell_momentum(icell,2)/cell_gas_density(icell)
                       *units->velocity/constants->kms);
        break;
    case OUT_CELL_SOUNDSPEED:
        return (float)( sqrt(
            cell_gas_pressure(icell)/cell_gas_density(icell)*constants->gamma )
                        *units->velocity/constants->kms); 
        break;
    case  OUT_CELL_PRESSURE:
        return (float)
            ( cell_gas_pressure(icell)*units->energy_density*constants->cc/constants->erg );
        break;
#if defined(RADIATIVE_TRANSFER) && (RT_TRANSFER_METHOD==RT_METHOD_OTVET)
#ifdef SNLUPDATE
    case  OUT_CELL_TAUUV:  
        return (float) cell_tauUV(icell);
        break;
    case  OUT_CELL_RADIATION_PRESSURE:
        return (float)
            ( cell_radiation_pressure(icell)
              /(cell_gas_pressure(icell)-cell_radiation_pressure(icell)) );
        break;
    case  OUT_CELL_URAD:
        return (float)
            ( cell_Urad(icell)*units->energy_density*constants->cc/constants->erg );
        break;
#endif
    case  OUT_CELL_FHI:
        return (float)
            cell_HI_fraction(icell);
        break;
    case  OUT_CELL_FH2:
        return (float)
            cell_H2_fraction(icell);
        break;
#ifdef RT_OTVET_SAVE_FLUX
    case  OUT_CELL_FLUX0:
        return (float)
            rt_flux[icell][0]*units->energy_density*units->length/constants->erg;
        break;
#endif        
#endif
    case OUT_CELL_MACH: 
        px2 = cell_momentum(icell,0)*cell_momentum(icell,0);
        py2 = cell_momentum(icell,1)*cell_momentum(icell,1);
        pz2 = cell_momentum(icell,2)*cell_momentum(icell,2);
        vmag=sqrt(px2+py2+pz2)/cell_gas_density(icell);
        cs=sqrt( cell_gas_pressure(icell)/cell_gas_density(icell)*constants->gamma );
        return (float) vmag/cs;
            /* return sqrt(px2+py2+pz2)/cell_gas_density(icell) */
            /* / sqrt( cell_gas_pressure(icell)/cell_gas_density(icell)*constants->gamma ); */
        break;
    case OUT_CELL_LEVEL: return (float)cell_level(icell);
        break;
    default:
        cart_error("bad flag");
        return -1;
        break;
    }
}


void dump_plane(int iflag, int out_level, int slice_axis_z, double pos_zero[nDim], double slice_region_hsize, FILE *output){
    int i, idummy ;
    float fdummy;
    int islice, icell;
    int slice_indx, slice_indy ;
    double pos[nDim];
    float *slice;
    double fact_hi_level=pow(2.0,out_level-min_level);
    int axis[2];
    const int endian_test=-99;
    const float unassigned=-98;
    const int nsgrid=slice_region_hsize*2/cell_size[out_level]; //odd
    char varname[string_size];
    axis[0] = axis_direction[slice_axis_z][0];
    axis[1] = axis_direction[slice_axis_z][1] ;
//    FILE *fp;
//    fp=fopen("checkrp.txt","w");

    
    fwrite( &endian_test, sizeof(int), 1, output );
    
    fwrite( &string_size, sizeof(int), 1, output );
    name_cell_property(iflag, varname);
    fwrite(varname, string_size*sizeof(char), 1, output );
    
    idummy  = nsgrid;
    fwrite( &idummy, sizeof(int), 1, output );
    idummy  = nsgrid;
    fwrite( &idummy, sizeof(int), 1, output );
    fdummy  = 2*slice_region_hsize*units->length/constants->kpc;
    fwrite( &fdummy, sizeof(float), 1, output );
#ifdef COSMOLOGY
    fdummy  = auni[min_level];
    fwrite( &fdummy, sizeof(float), 1, output );
    fdummy  = tphys_from_tcode(tl[min_level])*1e-6;
    fwrite( &fdummy, sizeof(float), 1, output );
#else
    cart_error("pick a different time variable to output");
#endif
    fdummy  = 0*units->time/constants->yr*1e-6;
    /* fdummy  = dtl[min_level]*units->time/constants->yr*1e-6; # Not available in analysis mode */
    fwrite( &fdummy, sizeof(float), 1, output );
        
    slice = cart_alloc(float, nsgrid*nsgrid);
        
    for(i=0; i<nsgrid*nsgrid; i++){
        slice[i] = unassigned;
    }

    pos[slice_axis_z] = pos_zero[slice_axis_z];
    
    for(slice_indy=0;slice_indy<nsgrid;slice_indy++){
        pos[axis[1]] = (slice_indy-nsgrid/2)/fact_hi_level + pos_zero[axis[1]]; 
            
        if(pos[axis[1]]>=num_grid){pos[axis[1]]-=num_grid;}
        if(pos[axis[1]]<0        ){pos[axis[1]]+=num_grid;}
            
        for(slice_indx=0;slice_indx<nsgrid;slice_indx++){
            pos[axis[0]] = (slice_indx-nsgrid/2)/fact_hi_level + pos_zero[axis[0]];
                
            if(pos[axis[0]]>=num_grid){pos[axis[0]]-=num_grid;}
            if(pos[axis[0]]<0        ){pos[axis[0]]+=num_grid;}
            
            islice = slice_indy*nsgrid + slice_indx;
            icell = cell_find_position_above_level(out_level+1,pos);
            cart_assert(icell!=-1);
            slice[islice] = value_cell_property(icell,iflag);
/*             if(slice[islice]==unassigned){ */
/*                 fprintf(fp,"%d %f %f %e\n",iflag, pos[axis[0]],pos[axis[1]],slice[i]); */
/*             } */
        }
    }
//    fclose(fp);
    
    int    debug_unassigned=0;
    for(i=0; i<nsgrid*nsgrid; i++){
        if(slice[i]==unassigned){
            debug_unassigned=1;
        }
        if(debug_unassigned==1){
            cart_debug("%d %e %e",iflag,slice[i],unassigned);
            exit(1);
        }
    }

    fwrite( slice, sizeof(float), nsgrid*nsgrid, output );
        
    cart_free(slice);
}




