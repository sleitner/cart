#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
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

qss_system *qss_alloc( size_t num_eqn,  
		void (* rates)(double, double *, void *, double *, double *),
		void (* adjust)(double, double *, void *) ) {

	qss_system *sys = (qss_system *)cart_alloc( sizeof(qss_system) );
	sys->buf = cart_alloc( 7*num_eqn*sizeof(double) );
	sys->y0 = sys->buf;
	sys->y1 = &sys->buf[num_eqn];
	sys->rs = &sys->buf[2*num_eqn];
	sys->a0 = &sys->buf[3*num_eqn];
	sys->w0 = &sys->buf[4*num_eqn];
	sys->a1 = &sys->buf[5*num_eqn];
	sys->w1 = &sys->buf[6*num_eqn];
	sys->num_eqn = num_eqn;
	sys->rates = rates;
	sys->adjust = adjust;
	return sys;
}

void qss_free( qss_system *sys ) {
	cart_free( sys->buf );
	cart_free( sys );
}

double qss_alpha( double tau ) {
	if ( tau < 0.01 ) {
		return 0.5 + tau/12 - tau*tau*tau/720;
	} else {
		return (1-(1-exp(-tau))/tau)/(1-exp(-tau));
	}
}

void qs1_step( qss_system *sys, double t, double dt, double yf[], void *params ) {
	int i;
	double tau, alp, pbar;

	for ( i = 0; i < sys->num_eqn; i++ ) {
		tau = dt*sys->a0[i];
		alp = qss_alpha(tau);
		sys->y1[i] = sys->y0[i] +
			dt*(sys->w0[i]-sys->a0[i]*sys->y0[i])/(1+alp*tau);
	}

	/* compute rates for corrector step */
	sys->rates( t, sys->y1, params, sys->w1, sys->a1 );

	for ( i = 0; i < sys->num_eqn; i++ ) {
		pbar = 0.5*(sys->a0[i]+sys->a1[i]);
		tau = pbar*dt;
		alp = qss_alpha(tau);
		yf[i] = sys->y0[i] +
			dt*(alp*sys->w1[i]+(1-alp)*sys->w0[i]- pbar*sys->y0[i]) /
                                (1+alp*tau);
                sys->rs[i] = yf[i] - sys->y1[i];
        }
}


void qsn_step( qss_system *sys, double t, double dt, double yf[], void *params ) {
        int i;
	int nIter, nCor;
	double rCor;
        double tau, alp, pbar;

	for ( i = 0; i < sys->num_eqn; i++ ) {
		tau = dt*sys->a0[i];
		alp = qss_alpha(tau);
		sys->y1[i] = sys->y0[i] +
			dt*(sys->w0[i]-sys->a0[i]*sys->y0[i])/(1+alp*tau);
		sys->rs[i] = sys->y1[i];
	}

	nIter = 0;
	nCor = 2;
	rCor = 1.0;

	do {
		/* compute rates for corrector step */
		sys->rates( t, sys->y1, params, sys->w1, sys->a1 );

		for ( i = 0; i < sys->num_eqn; i++ ) {
			pbar = 0.5*(sys->a0[i]+sys->a1[i]);
			tau = pbar*dt;
			alp = qss_alpha(tau);
			yf[i] = sys->y0[i] +
				dt*(alp*sys->w1[i]+(1-alp)*sys->w0[i]-pbar*sys->y0[i]) /
				(1+alp*tau);
		}

		/* Nick's black magic, designed to advance further when mode is oscillating */
		if ( (yf[0]-sys->y0[0])*(sys->rs[0]-sys->y0[0]) < 0.0 ) {
			nCor = 4;
			rCor = 0.2;
		}

		nIter++;

		if ( nIter < nCor ) {
			for ( i = 0; i < sys->num_eqn; i++ ) {
				sys->y1[i] = (1-rCor)*sys->y1[i] + rCor*yf[i];
			}
		}
	} while ( nIter < nCor );

	for ( i = 0; i < sys->num_eqn; i++ ) {
		sys->rs[i] = yf[i] - sys->rs[i];
        }
}

void qss_solve( qss_system *sys, double t_begin, double t_end, double y[], const double err[], void *params ) {
	int i;
	int nstep = 0;
	double errmax;
	double erri;

	double t = t_begin;
	double dt = t_end-t_begin;

	do {
		for ( i = 0; i < sys->num_eqn; i++ ) {
			sys->y0[i] = y[i];
		}

		/* compute initial rates */
		sys->rates( t, sys->y0, params, sys->w0, sys->a0 );

		while ( 1 ) {
			/* take step */
			if ( nstep < 300 ) {
				qs1_step( sys, t, dt, y, params );
			} else {
				qsn_step( sys, t, dt, y, params );
			}

			/* check for errors */
			errmax = fabs(sys->rs[0])/max(y[0],1e-30)/err[0];
			for ( i = 1; i < sys->num_eqn; i++ ) {
				erri = fabs(sys->rs[i])/max(y[i],1e-30)/err[i];
				if ( erri > errmax ) {
					errmax = erri;
				}
			}

			if ( errmax > 1.0 ) {
				dt = min( 0.9*dt/sqrt(2.+errmax), t_end-t );
			} else {
				t += dt;
				nstep++;
				dt = min( 0.9*dt/pow(0.01+errmax,0.3), t_end-t );
				break;
			}
		}

		/* adjust variables for max/min */
		sys->adjust( t, y, params );
	} while ( t < t_end ); 
}


/*
//  Changed by Gnedin to catch memory leaks
*/

#ifdef DEBUG_MEMORY_USE
void dmuRegister(void *ptr, unsigned long size, const char *file, int line);
void dmuUnRegister(void *ptr);
void dmuPrintRegistryContents();
#endif


void* cart_alloc_at_location(size_t size, const char *file, int line)
{
  void *ptr;

  if(size > 0)
    {
      ptr = malloc( size );

      if(ptr == NULL)
	{
#ifdef DEBUG_MEMORY_USE
	  dmuPrintRegistryContents();
#endif
	  cart_error( "Failure allocating %d bytes", size );
	}

#ifdef DEBUG_MEMORY_USE
      dmuRegister(ptr,size,file,line);
#endif

      return ptr;
    }
  else
    {
      return NULL;
    }
}


void cart_free(void *ptr)
{
  if(ptr != NULL)
    {
      free(ptr);
#ifdef DEBUG_MEMORY_USE
      dmuUnRegister(ptr);
#endif
    }
}

#ifdef DEBUG_MEMORY_USE

int dmuRegistrySize = 0;
int dmuNumItems = 0;

struct dmuItem
{
  void *Ptr;
  unsigned long Size;
  const char *File;
  int Line;
};

struct dmuItem *dmuRegistry = NULL;


void dmuRegister(void *ptr, unsigned long size, const char *file, int line)
{
  struct dmuItem *tmp;

  if(dmuNumItems == dmuRegistrySize)
    {
      /*
      //  Extend registry
      */
      dmuRegistrySize += 1000;
      tmp = malloc(dmuRegistrySize*sizeof(struct dmuItem));
      if(tmp == NULL)
	{
	  cart_error( "DMU: failure extendion allocation registry to %d items, %d bytes",dmuRegistrySize,dmuRegistrySize*sizeof(struct dmuItem));
	}
      if(dmuNumItems > 0)
	{
	  memcpy(tmp,dmuRegistry,dmuNumItems*sizeof(struct dmuItem));
	  free(dmuRegistry);
	}
      dmuRegistry = tmp;
    }
      
  dmuRegistry[dmuNumItems].Ptr = ptr;
  dmuRegistry[dmuNumItems].Size = size;
  dmuRegistry[dmuNumItems].File = file;
  dmuRegistry[dmuNumItems].Line = line;
  dmuNumItems++;
}


void dmuUnRegister(void *ptr)
{
  int i, j;
  
  for(i=0; i<dmuNumItems; i++)
    {
      if(dmuRegistry[i].Ptr == ptr) break;
    }

  if(i == dmuNumItems)
    {
      cart_error("DMU: freeing unregistered pointer %x",ptr);
    }

  for(j=i; j<dmuNumItems-1; j++)
    {
      dmuRegistry[j] = dmuRegistry[j+1];
    }
  dmuNumItems--;
}


void dmuPrintRegistryContents()
{
  int i;
  unsigned long tot = 0UL;

  for(i=0; i<dmuNumItems; i++)
    {
      cart_debug("DMU: ptr=%p, size=%lu, file=%s, line=%d",dmuRegistry[i].Ptr,dmuRegistry[i].Size,dmuRegistry[i].File,dmuRegistry[i].Line);
      tot += dmuRegistry[i].Size;
    }

  cart_debug("DMU: total allocated size=%lu",tot);
}


unsigned long dmuReportAllocatedMemory()
{
  int i;
  unsigned long tot = 0UL;

  for(i=0; i<dmuNumItems; i++)
    {
      tot += dmuRegistry[i].Size;
    }

  return tot;
}

#endif  /* DEBUG_MEMORY_USE */
