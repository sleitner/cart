#include "config.h"

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "hydro.h"
#include "hydro_tracer.h"
#include "io.h"
#include "io_cart.h"
#include "iterators.h"
#include "load_balance.h"
#include "parallel.h"
#include "particle.h"
#include "refinement.h"
#include "rt_io.h"
#include "sfc.h"
#include "starformation.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "io_music.h"

/* mimicking io_cart */
#ifdef HYDRO
#define FUNCTION                      read_music_gas_particles_double
#define PARTICLE_FLOAT                double
#define MPI_PARTICLE_FLOAT            MPI_DOUBLE
#include "io_music.def"
void read_music_gas_particles( char *header_filename, char *data_filename, int num_sfcs, int *sfc_list ) {
    read_music_gas_particles_double(header_filename,data_filename,num_sfcs,sfc_list);
}
#endif /* HYDRO */
