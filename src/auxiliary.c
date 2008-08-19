#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <math.h>

#include <mpi.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

#include "defs.h"
#include "tree.h"
#include "auxiliary.h"
#include "parallel.h"
#include "io.h"

unsigned long int rng_seed = 0L;
gsl_rng *cart_random_generator;

extern int current_step_level;


void init_auxiliary() {
	FILE *state;
	char filename[256];

	cart_random_generator = gsl_rng_alloc (gsl_rng_mt19937);

	/* attempt to reload state information (just use natural seeding otherwise) */
	sprintf( filename, "%s/rng_state_%03u.dat", output_directory, local_proc_id );
	state = fopen( filename, "r" );
	if ( state != NULL ) {
		gsl_rng_fread( state, cart_random_generator );
		fclose(state);
	} else {
		gsl_rng_set( cart_random_generator, rng_seed );
	}
}

void save_auxiliary() {
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

double gsl_function_wrapper( double x, void *params ) {
        double (*f)(double) = (double (*)(double))params;
        return  f(x);
}

double integrate( double (*f)(double), double a, double b, double epsrel, double epsabs ) {
        double result, error;
        gsl_function F;
        gsl_integration_workspace *w;

	if ( a == b ) {
		return 0.0;
	}

	w = gsl_integration_workspace_alloc(1000);

        F.function = gsl_function_wrapper;
        F.params = f;

	gsl_integration_qag(&F, a, b, epsrel, epsabs, 1000, 6,
                w, &result, &error);
                                                                                                                                                            
        gsl_integration_workspace_free(w);
                                                                                                                                                            
        return result;
}

#define MAX_ITER 1000
                                                                                                                                                            
double root_finder( double (*f)(double), double a, double b, double epsrel, double epsabs ) {
        int status,i;
        double root;
        const gsl_root_fsolver_type *T;
        gsl_root_fsolver *s;
        gsl_function F;
                                                                                                                                                            
        F.function = gsl_function_wrapper;
        F.params = f;
                                                                                                                                                            
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc (T);
        gsl_root_fsolver_set (s, &F, a, b);
                                                                                                                                                            
        for ( i = 0; i < MAX_ITER; i++ ) {
                status = gsl_root_fsolver_iterate (s);
                status = gsl_root_test_interval (gsl_root_fsolver_x_lower(s),
                                                 gsl_root_fsolver_x_upper(s),
                                                 epsrel, epsabs);

                if (status == GSL_SUCCESS) {
                        root = gsl_root_fsolver_root(s);
                        gsl_root_fsolver_free(s);
                        return root;
                }
        }

	cart_error("Did not reach root after %u iterations!", MAX_ITER);
	return 0.0;
}

#undef MAX_ITER

int compare_ints( const void *a, const void *b ) {
        if ( *(int *)a == -1 ) {
                return 1;
        } else if ( *(int *)b == -1 ) {
                return -1;
        } else {
                return ( *(int *)a - *(int *)b );
        }
}

int compare_floats( const void *a, const void *b ) {
	if ( *(float *)a > *(float *)b ) {
		return 1;
	} else if ( *(float *)a < *(float *)b ) {
		return -1;
	} else {
		return 0;
	}
}

int nearest_int( double x ) {
	int ix;
	double frac;

	ix = (int)x;
	frac = x - (double)ix;

	if ( frac < 0.5 ) {
		return ix;
	} else {
		return ix+1;
	}
}

double cart_rand() {
	int proc;
	double ret;

#ifdef UNIQUE_RAND
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
	
void *cart_alloc( size_t size ) {
	void *ptr;

	if ( size > 0 ) {
		ptr = malloc( size );

		if ( ptr == NULL ) {
			cart_error( "Failure allocating %d bytes", size );
		}

		return ptr;
	} else {
		return NULL;
	}
}

void cart_free( void *ptr ) {
	if ( ptr != NULL ) {
		free( ptr );
	}
}

void cart_error( const char *fmt, ... ) {
	char message[256];

	va_list args;

	va_start( args, fmt );
	vsnprintf( message, 256, fmt, args );
	fprintf(stderr, "%u: %s\n", local_proc_id, message );
	fflush(stderr);
	va_end( args );
	MPI_Abort( MPI_COMM_WORLD, 1 );
}

#ifndef NDEBUG
void cart_debug( const char *fmt, ... ) {
        int i;
	char message[256], prompt[256];
	va_list args;

	/* prompt */
	if(current_step_level > -1)
	  {
	    strcpy(prompt,"> ");
	    for(i=1; i<=current_step_level; i++)
	      {
		if(i%5 == 0) strcat(prompt,": "); else strcat(prompt,". ");
	      }
	  }
	else
	  {
	    prompt[0] = 0;
	  }

	va_start( args, fmt );
	vsnprintf( message, 256, fmt, args );
	fprintf( stdout, "%u: %s%s\n", local_proc_id, prompt, message );
	va_end(args);
}
#endif
