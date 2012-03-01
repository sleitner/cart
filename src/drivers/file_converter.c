#include "config.h"

#include "../base/auxiliary.h"
#include "../base/cosmology.h"
#include "../base/times.h" 
#include "../base/units.h"

#include "../core/cell_buffer.h"
#include "../core/parallel.h"
#include "../core/timing.h"
#include "../core/tree.h"


extern const char* executable_name;

void config_init();


void read_file(const char* fname);
void write_file(const char* fname);


int drive() {

	if ( num_options != 2 ) {
	  cart_error("Usage: %s <input> <output>",executable_name);
	}

	init_timers();

	config_init();

	init_parallel_grid();
	init_tree();
	init_cell_buffer();

	read_file(options[0]);
	
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

	write_file(options[1]);

	return 0;
}
