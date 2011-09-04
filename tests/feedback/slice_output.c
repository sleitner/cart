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

#include "slice_output.h"


const int axis_direction[nDim][nDim-1] = {
#if nDim == 2
    {1},{0}
#elif nDim == 3
    {1,2},{0,2},{0,1}
#else
    #error  "Unsupported number of dimensions for other_axis"
#endif
};


int is_cell_in_slice_plane(double pos_zero[nDim], int slice_axis_z, double pos[nDim], int level) {
    return (fabs(pos[slice_axis_z] + 1e-20 - pos_zero[slice_axis_z]) <= cell_size[level]/2.0) ;
}

float value_cell_property(int icell,int iflag){
    switch(iflag){
    case OUT_CELL_DENSITY: return (float)cell_gas_density(icell);
        break;
    case OUT_CELL_INTERNAL_ENERGY: return (float)cell_gas_internal_energy(icell);
        break;
    case OUT_CELL_MOMENTUM+0: return (float)cell_momentum(icell,0);
        break;
    case OUT_CELL_MOMENTUM+1: return (float)cell_momentum(icell,1);
        break;
    case OUT_CELL_MOMENTUM+2: return (float)cell_momentum(icell,2);
        break;
    case OUT_CELL_PRESSURE: return (float)cell_gas_pressure(icell);
        break;
    case OUT_CELL_LEVEL: return (float)cell_level(icell);
        break;
    default:
        cart_error("bad flag");
        return -1;
        break;
    }
}
void dump_plane(int iflag, int out_level, int slice_axis_z, double pos_zero[nDim], FILE *output){//pos[slice_axis_z]=galaxy center, out_level
    int i, size, nsgrid;
    //also need to adjust to take cells level with your selected cell in cell_size need cell_size at eachlevel so block_sign[level]
    int slice_indx, slice_indy, fill_level;
    int ix,iy;
    double pos[nDim];
    float block_delta;
    float *slice;
    float fact_hi_level=pow(2.0,out_level-min_level);
    int slice_axis_x=axis_direction[slice_axis_z][0];
    int slice_axis_y=axis_direction[slice_axis_z][1];
//lower cell_delta block, output everything within 80kpc? for each projection axis. Around each cell with some critical value of density that is the most extreme in its rootgrid? just the first for now... could do that, eliminate all cells from the list then step to the next...
    const int endian_test=-99;
        
        
    fwrite( &endian_test, sizeof(int), 1, output );
    nsgrid = num_grid*fact_hi_level;
    size  = nsgrid;
    fwrite( &size, sizeof(int), 1, output );
    size  = nsgrid;
    fwrite( &size, sizeof(int), 1, output );
        
    slice = cart_alloc(float, nsgrid*nsgrid);
    for(i=0; i<nsgrid*nsgrid; i++){ slice[i] = 0;}
        
    MESH_RUN_DECLARE(level, icell);
    MESH_RUN_OVER_LEVELS_BEGIN(level,min_level,out_level);
    MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(icell);
    if(cell_is_leaf(icell) || level == out_level){
        cell_center_position(icell,pos);
        if(is_cell_in_slice_plane(pos_zero,slice_axis_z,pos,level)) { 
            slice_indx = (int)((pos[slice_axis_x]-cell_size[level]/2.0)*fact_hi_level);
            slice_indy = (int)((pos[slice_axis_y]-cell_size[level]/2.0)*fact_hi_level) ;
            cart_assert(slice_indx ==
                        (pos[slice_axis_x]-cell_size[level]/2.0)*fact_hi_level);
            cart_assert(slice_indy ==
                        (pos[slice_axis_y]-cell_size[level]/2.0)*fact_hi_level);
            //fill holes (if you want)
            fill_level=out_level-level;
            for(iy=0;iy<pow(2,fill_level);iy++){
                for(ix=0;ix<pow(2,fill_level);ix++){
                    cart_assert(slice[(slice_indy+iy)*nsgrid + slice_indx+ix] == 0) ;
                    slice[(slice_indy+iy)*nsgrid + slice_indx+ix] = value_cell_property(icell,iflag);
                }
            }
        }
    }
    MESH_RUN_OVER_CELLS_OF_LEVEL_END;
    MESH_RUN_OVER_LEVELS_END;
    
    fwrite( slice, sizeof(float), nsgrid*nsgrid, output );
//    cart_debug("slice in cell prop %d output",iflag);
        
    cart_free(slice);
}




/* void dump_slice(int iflag, FILE *output, int out_level, int slice_axis_z){//pos[slice_axis_z]=galaxy center, out_level */
/*     int i, size, nsgrid; */
/*     double pos[nDim]; */
/*     int const block_sign=-1; */
/* //lower cell_delta block, output everything within 80kpc? for each projection axis. Around each cell with some critical value of density that is the most extreme in its rootgrid? just the first for now... could do that, eliminate all cells from the list then step to the next... */
/*     //also need to adjust to take cells level with your selected cell in cell_size need cell_size at eachlevel so block_sign[level] */
    
/*     int block_level, block_cell, slice_indx, slice_indy, fill_level; */
/*     int ix,iy; */
/*     float block_delta; */
/*     float *slice; */
/*     float fact_hi_level=pow(2.0,out_level-min_level); */
        
        
/*     fwrite( &endian_test, sizeof(int), 1, output ); */
/*     nsgrid = num_grid*fact_hi_level; */
/*     size  = nsgrid; */
/*     fwrite( &size, sizeof(int), 1, output ); */
/*     size  = nsgrid; */
/*     fwrite( &size, sizeof(int), 1, output ); */
        
/*     slice = cart_alloc(float, nsgrid*nsgrid); */
/*     for(i=0; i<nsgrid*nsgrid; i++){ slice[i] = 0;} */
        
/*     MESH_RUN_DECLARE(level, icell); */
/*     MESH_RUN_OVER_LEVELS_BEGIN(level,min_level,out_level); */
/*     MESH_RUN_OVER_CELLS_OF_LEVEL_BEGIN(icell); */
/*     if(cell_is_leaf(icell) || level == out_level){ */
/*         cell_center_position(icell,pos); */
/*         if((int)pos[slice_axis_z] == num_grid/2){  */
/*             block_level = level; */
/*             block_cell = icell; */
/*             block_delta = cell_delta[cell_child_number(block_cell)][slice_axis_z]; */
/*             // check that cell_delta is negative for current and all parents-> */
/*             //(take the cells closest to the axis) */
/*             while(block_delta*block_sign>0 && block_level>min_level){ */
/*                 block_level--; */
/*                 block_cell = cell_parent_cell(block_cell)  ; */
/*                 block_delta = cell_delta[cell_child_number(block_cell)][slice_axis_z]; */
/*             } */
/*             if(block_level==min_level){ */
/*                 //made it to root and were in right cell the whole way */
/*                 slice_indx = (int)((pos[slice_axis_x]-cell_size[level]/2.0)*fact_hi_level); */
/*                 slice_indy = (int)((pos[slice_axis_y]-cell_size[level]/2.0)*fact_hi_level) ; */
/*                 cart_assert(slice_indx == */
/*                             (pos[slice_axis_x]-cell_size[level]/2.0)*fact_hi_level); */
/*                 cart_assert(slice_indy == */
/*                             (pos[slice_axis_y]-cell_size[level]/2.0)*fact_hi_level); */
/*                 //fill holes (if you want) */
/*                 fill_level=out_level-level; */
/*                 for(iy=0;iy<pow(2,fill_level);iy++){ */
/*                     for(ix=0;ix<pow(2,fill_level);ix++){ */
/*                         cart_assert(slice[(slice_indy+iy)*nsgrid + slice_indx+ix] == 0) ; */
/*                         slice[(slice_indy+iy)*nsgrid + slice_indx+ix] = cell_vars[icell][iflag]; */
/*                     } */
/*                 } */
                
/*             } */
/*         } */
/*     } */
/*     MESH_RUN_OVER_CELLS_OF_LEVEL_END; */
/*     MESH_RUN_OVER_LEVELS_END; */
        
/*     fwrite( slice, sizeof(float), nsgrid*nsgrid, output ); */
/* //    cart_debug("slice in cell prop %d output",iflag); */
        
/*     cart_free(slice); */
/* } */




