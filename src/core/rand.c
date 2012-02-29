#include "config.h"

#include <math.h>
#include <stdio.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>

#include "auxiliary.h"


unsigned long int rng_seed = 0L;
gsl_rng *cart_random_generator;


void init_rand() {
	FILE *state;
	char filename[256];

	cart_random_generator = gsl_rng_alloc (gsl_rng_mt19937);

	/* attempt to reload state information (just use natural seeding otherwise) */
	sprintf( filename, "%s/rng_state_%03u.dat", logfile_directory, local_proc_id );
	state = fopen( filename, "r" );
	if ( state != NULL ) {
		gsl_rng_fread( state, cart_random_generator );
		fclose(state);
	} else {
		gsl_rng_set( cart_random_generator, rng_seed );
	}
}

void save_rand() {
	FILE *state;
	char filename[256];

#ifdef STARFORM
	sprintf( filename, "%s/rng_state_%03u.dat", logfile_directory, local_proc_id );
	state = fopen( filename, "w" );
	if ( state == NULL ) {
		cart_error("Unable to open %s for saving rng state", filename );
	}

	gsl_rng_fwrite ( state, cart_random_generator );
	fclose(state);
#endif
}

double cart_rand() {
	double ret;

#ifdef UNIQUE_RAND
	int proc;

	/* ensure unique random number sequence generated on each proc */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( proc == local_proc_id ) {
			ret = gsl_rng_uniform( cart_random_generator );
		} else {
			gsl_rng_uniform( cart_random_generator );
		}
	}
#else
	ret = gsl_rng_uniform( cart_random_generator );
#endif /* UNIQUE_RAND */

	return ret;
}

double cart_rand_lognormal(double sigma) {
	double ret;
	double zeta = -0.5*sigma*sigma;

#ifdef UNIQUE_RAND
	int proc;

	/* ensure unique random number sequence generated on each proc */
	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( proc == local_proc_id ) {
		  ret = gsl_ran_lognormal( cart_random_generator, zeta, sigma );
		} else {
			gsl_ran_lognormal( cart_random_generator, zeta, sigma );
		}
	}
#else
	ret = gsl_ran_lognormal( cart_random_generator, zeta, sigma );
#endif /* UNIQUE_RAND */

	return ret;
}
