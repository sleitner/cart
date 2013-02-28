#include "config.h"

#include "../base/auxiliary.h"
#include "../base/cosmology.h"
#include "../base/rand.h"
#include "../base/times.h" 
#include "../base/units.h"

#include "../core/cell_buffer.h"
#include "../core/parallel.h"
#include "../core/particle.h"
#include "../core/plugin.h"
#include "../core/timing.h"
#include "../core/tree.h"


extern const char* executable_name;

void config_init();
const plugin_t* add_plugin(int id)
{
  return NULL;
}


void init();
void read_file();
void write_file();


int drive() {

	init();

	init_timers();

	//config_init();

	init_rand();
	init_parallel_grid();
	init_cell_buffer();
	init_tree();
#ifdef PARTICLES
	init_particles();
#endif
#ifdef HYDRO_TRACERS
    init_hydro_tracers();
#endif /* HYDRO_TRACERS */

	read_file();
	
	units_init();
#ifdef COSMOLOGY
	abox[min_level] = abox_from_tcode(tl[min_level]);
	auni[min_level] = auni_from_tcode(tl[min_level]);
#else
	abox[min_level] = auni[min_level];
#endif
	units_update(min_level);

	cart_debug("tl = %e", tl[min_level] );
	cart_debug("abox = %e", abox[min_level] );
	cart_debug("done reading data...");

	write_file();

	return 0;
}
