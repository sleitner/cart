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
    {1,2},{0,2},{0,1}
#else
    #error  "Unsupported number of dimensions for other_axis"
#endif
};

int is_cell_in_slice_region(double pos_zero[nDim], int slice_axis_x, int slice_axis_y, double slice_region_hsize, double pos[nDim],int level ) {
        return (compute_distance_periodic_1d(
                    pos[slice_axis_x] , pos_zero[slice_axis_x])
                < slice_region_hsize &&
                compute_distance_periodic_1d(
                    pos[slice_axis_y] , pos_zero[slice_axis_y])
                < slice_region_hsize ) ;
}

int is_cell_in_slice_plane(double pos_zero[nDim], int slice_axis_z, double pos[nDim], int level) {
    double disp;
    disp=pos[slice_axis_z]-pos_zero[slice_axis_z];
        return (disp >=0 && disp <= cell_size[level]/2.0) ;
}

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


void dump_plane(int iflag, int out_level, int slice_axis_z, double pos_zero[nDim], double slice_region_hsize, FILE *output){//pos[slice_axis_z]=galaxy center, out_level
    int i, idummy ;
    float fdummy;
    //also need to adjust to take cells level with your selected cell in cell_size need cell_size at eachlevel so block_sign[level]
    int islice, icell;
    int slice_indx, slice_indy ;
    double pos[nDim];
    float *slice;
    float fact_hi_level=pow(2.0,out_level-min_level);
    int axis[2] = { axis_direction[slice_axis_z][0],  axis_direction[slice_axis_z][1] };
//lower cell_delta block, output everything within 80kpc? for each projection axis. Around each cell with some critical value of density that is the most extreme in its rootgrid? just the first for now... could do that, eliminate all cells from the list then step to the next...
    const int endian_test=-99;
    const int unassigned=-99;
    const int nsgrid=slice_region_hsize*2/cell_size[out_level]; //should be odd
    
    
    fwrite( &endian_test, sizeof(int), 1, output );
//    nsgrid = num_grid*fact_hi_level;
    idummy  = nsgrid;
    fwrite( &idummy, sizeof(int), 1, output );
    idummy  = nsgrid;
    fwrite( &idummy, sizeof(int), 1, output );
    fdummy  = 2*slice_region_hsize*units->length/constants->kpc;
    fwrite( &fdummy, sizeof(float), 1, output );
    fdummy  = (tphys_from_tcode(tl[min_level])-tphys_from_tcode(t_init))*1e-6;
    fwrite( &fdummy, sizeof(float), 1, output );
    cart_debug("nsgrid=%d,out_level=%d, time=%e",nsgrid,out_level, fdummy);
        
    slice = cart_alloc(float, nsgrid*nsgrid);
        
    for(i=0; i<nsgrid*nsgrid; i++){
        slice[i] = unassigned;
    }

    pos[slice_axis_z] = pos_zero[slice_axis_z];
    
    for(slice_indy=0;slice_indy<nsgrid;slice_indy++){
        pos[axis[1]] = (slice_indy-nsgrid/2)/fact_hi_level +
            pos_zero[axis[1]]; 
        if(pos[axis[1]]>=num_grid){pos[axis[1]]-=num_grid;}
        if(pos[axis[1]]<0        ){pos[axis[1]]+=num_grid;}
            
        for(slice_indx=0;slice_indx<nsgrid;slice_indx++){
            pos[axis[0]] = (slice_indx-nsgrid/2)/fact_hi_level +
                pos_zero[axis[axis[0]]];
            if(pos[axis[0]]>=num_grid){pos[axis[0]]-=num_grid;}
            if(pos[axis[0]]<0        ){pos[axis[0]]+=num_grid;}
            
            islice = slice_indy*nsgrid + slice_indx;
            icell = cell_find_position_above_level(out_level,pos);
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




