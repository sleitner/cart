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

#include "disk.h"

#ifndef STARFORM
#error need stars defined for feedback tests
#endif

void writebegin();

plugin_t outPlugin = {NULL};

void everystep();

const plugin_t* add_plugin(int id){
    if(id==0){
	/* If hitting CFL use writebegin to choose shorter timestep */
        outPlugin.RunBegin = writebegin;
        return &outPlugin;
    }else{
        return NULL;
    }
}


void writebegin(){
    int icell;
    write_restart( WRITE_SAVE, WRITE_SAVE, WRITE_SAVE );
}

void run_output(){
    cart_debug("");   
}

