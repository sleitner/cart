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

#include "output_slice.h"


const int axis_direction[nDim][nDim-1] = {
#if nDim == 2
    {1},{0}
#elif nDim == 3
    {1,2},{2,0},{0,1}
#else
    #error  "Unsupported number of dimensions for other_axis"
#endif
};

float value_cell_property(int icell,int iflag){
    switch(iflag){
    case OUT_CELL_DENSITY: return (float)cell_gas_density(icell)*units->number_density/constants->cc;
        break;
    case OUT_CELL_INTERNAL_ENERGY: return (float)cell_gas_temperature(icell)*units->temperature/constants->K;
        break;
    case OUT_CELL_MOMENTUM+0: return (float)cell_momentum(icell,0)/cell_gas_density(icell)*units->velocity/constants->kms;
        break;
    case OUT_CELL_MOMENTUM+1: return (float)cell_momentum(icell,1)/cell_gas_density(icell)*units->velocity/constants->kms;
        break;
    case OUT_CELL_MOMENTUM+2: return (float)cell_momentum(icell,2)/cell_gas_density(icell)*units->velocity/constants->kms;
        break;
    case OUT_CELL_PRESSURE: return (float)cell_gas_pressure(icell)*units->energy_density/constants->erg/pow(constants->pc,3);
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
    const int unassigned=-99;
    const int nsgrid=slice_region_hsize*2/cell_size[out_level]; //odd
    axis[0] = axis_direction[slice_axis_z][0];
    axis[1] = axis_direction[slice_axis_z][1] ;
    cart_debug("permuting axes, %d %d %d", axis[0],axis[1],slice_axis_z);

    cart_debug("in output_slice; axes (xyz)=%d%d%d",axis[0],axis[1],axis[2]);

    
    fwrite( &endian_test, sizeof(int), 1, output );
    idummy  = nsgrid;
    fwrite( &idummy, sizeof(int), 1, output );
    idummy  = nsgrid;
    fwrite( &idummy, sizeof(int), 1, output );
    fdummy  = 2*slice_region_hsize*units->length/constants->kpc;
    fwrite( &fdummy, sizeof(float), 1, output );
    fdummy  = tphys_from_tcode(tl[min_level])*1e-6;
    fwrite( &fdummy, sizeof(float), 1, output );
//    fdummy  = (tphys_from_tcode(tl[min_level])-tphys_from_tcode(t_init))*1e-6;
    fdummy  = dtl[min_level]*units->time/constants->yr*1e-6;
    fwrite( &fdummy, sizeof(float), 1, output );
    cart_debug("nsgrid=%d,out_level=%d, time=%e",nsgrid,out_level, fdummy);
        
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
        }
    }

    for(i=0; i<nsgrid*nsgrid; i++){
        cart_assert(slice[i]!=unassigned);
    }

    fwrite( slice, sizeof(float), nsgrid*nsgrid, output );
        
    cart_free(slice);
}




