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

gsl_rng *cart_random_generator;

void init_rand() {
	cart_random_generator = gsl_rng_alloc (gsl_rng_mt19937);
}

void cart_rand_set_seed( unsigned long int seed ) {
	cart_assert( cart_random_generator != NULL );
	gsl_rng_set( cart_random_generator, seed );
}

void cart_rand_load_state( const char *state_filename, int fail_on_error ) {
	FILE *state = fopen( state_filename, "r" );
	if ( state != NULL ) {
		gsl_rng_fread( state, cart_random_generator );
		fclose(state);
	} else if ( fail_on_error ) {
		cart_error("Unable to load random number generator state file %s!", state_filename );
	} /* otherwise silently skip over the error */
}

void cart_rand_save_state( const char *state_filename ) {
	FILE *state;

	state = fopen( state_filename, "w" );
	if ( state == NULL ) {
		cart_error("Unable to open %s for saving randomnumber generator state!", state_filename );
	}
	gsl_rng_fwrite ( state, cart_random_generator );
	fclose(state);
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

unsigned int cart_rand_poisson( double mu ) {
	int ret;

#ifdef UNIQUE_RAND
	int proc;

	for ( proc = 0; proc < num_procs; proc++ ) {
		if ( proc == local_proc_id ) {
			ret = gsl_ran_poisson( cart_random_generator, mu );
		} else {
			gsl_ran_poisson( cart_random_generator, mu );
		}
	}
#else
    ret = gsl_ran_poisson( cart_random_generator, mu );
#endif /* UNIQUE_RAND */

	return ret;
}

