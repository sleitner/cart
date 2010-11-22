#include "config.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_rng.h>

#include "auxiliary.h"
#include "parallel.h"
#include "io.h"
#include "timing.h"


unsigned long int rng_seed = 0L;
gsl_rng *cart_random_generator;

extern int current_step_level;


int num_options = 0;
char **options = NULL;


void init_auxiliary() {
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
	F.params = (void *)f;  /* NG: BAD!!! Unsafe type cast */

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
	F.params = (void *)f;

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
	

void cart_error( const char *fmt, ... ) {
	char message[256];
	char filename[256];
	FILE *f;

	va_list args;

	va_start( args, fmt );
	vsnprintf( message, 256, fmt, args );
	fprintf(stderr, "%u: %s\n", local_proc_id, message );
	fflush(stderr);
	va_end( args );

	sprintf(filename,"%s/stdout.%03u.log",logfile_directory,local_proc_id);
	f = fopen(filename,"a");
	if (f != 0) {
		fprintf(f,"ERROR: %s\n",message);
		fclose(f);
	}

	MPI_Abort( MPI_COMM_WORLD, 1 );
}

#ifndef NDEBUG
void cart_debug( const char *fmt, ... ) {
	int i;
	char message[256], prompt[256];
	va_list args;
	static int first_time = 1;
	char filename[256];
	FILE *f;

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

	sprintf(filename,"%s/stdout.%03u.log",logfile_directory,local_proc_id);
	if(first_time)
	  {
	    first_time = 0;
	    f = fopen(filename,"w");
	  }
	else f = fopen(filename,"a");
	if(f != 0)
	  {
	    fprintf(f,"%s%s\n",prompt,message);
	    fclose(f);
	  }
}
#endif

qss_system *qss_alloc( size_t num_eqn,  
		void (* rates)(double, double *, void *, double *, double *),
		void (* adjust)(double, double *, void *) ) {

	qss_system *sys = cart_alloc(qss_system, 1 );
	sys->buf = cart_alloc(double, 7*num_eqn );
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
	} while ( t < t_end && nstep < 100000 ); 
}


const char* check_option0(const char* option, const char* name)
{
  static const char *empty = "";

  int len = strlen(name);

  if(option[0]!='-' || strcmp(option+1,name)!=0) return NULL;

  return empty;
}


const char* check_option1(const char* option, const char* name, const char *default_value)
{
  int len = strlen(name);

  if(option[0]!='-' || strncmp(option+1,name,len)!=0) return NULL;

  if(strlen(option+1) == len)
    {
      if(default_value == NULL)
	{
	  cart_error("Option %s has no default value; a value must be set as -%s=<value>",name,name);
	  return NULL;
	}
      return default_value;
    }

  if(option[len+1]!='=' || strlen(option)<len+3)
    {
      if(default_value == NULL)
	{
	  cart_error("Valid format for option %s is: -%s=<value>",name,name);
	}
      else
	{
	  cart_error("Valid format for option %s is: -%s[=<value>]",name,name);
	}
      return NULL;
    }
  else
    {
      return option + len + 2;
    }
}


/*
//  Changed by Gnedin to catch memory leaks
*/

#ifdef DEBUG_MEMORY_USE
void dmuRegister(void *ptr, unsigned long size, const char *file, int line);
void dmuUnRegister(void *ptr, const char *file, int line);
void dmuPrintRegistryContents();
#endif


void* cart_alloc_worker(size_t size, const char *file, int line)
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
	  cart_error( "Failure allocating %d bytes in file %s, line %d", size, file, line );
	}

#ifdef DEBUG_MEMORY_USE
      dmuRegister(ptr,size,file,line);

#ifdef DEBUG_MEMORY_USE_VERBOSE
      /* if allocating significant chunk (i.e. not skiplist nodes, report it */
      if( size > 65536)
	{
	  cart_debug("allocated %d bytes, %ld total",size,dmuReportAllocatedMemory() );
	}
#endif
#endif

      return ptr;
    }
  else
    {
      return NULL;
    }
}


void cart_free_worker(void *ptr, const char *file, int line)
{
  if(ptr != NULL)
    {
      free(ptr);
#ifdef DEBUG_MEMORY_USE
      dmuUnRegister(ptr,file,line);
#endif
    }
}

#ifdef DEBUG_MEMORY_USE


struct dmuItem
{
  void *Ptr;
  unsigned long Size;
  const char *File;
  int Line;
};


struct dmuRegistry
{
  int Size;
  int NumItems;
  struct dmuItem *Data;
};

int dmuSmallLimit = 100;
int dmuMaxSmallChunk = 0;
unsigned long dmuNumChunks[2] = { 0UL, 0UL };

struct dmuRegistry dmuSmall = { 0, 0, NULL};
struct dmuRegistry dmuLarge = { 0, 0, NULL};


void dmuCheckRegistry(struct dmuRegistry *registry)
{
  struct dmuItem *tmp;

  cart_assert(registry);

  /*
  //  Create registry initially
  */
  if(registry->Size == 0)
    {
      registry->Size = 1000;
      registry->Data = (struct dmuItem *)malloc(registry->Size*sizeof(struct dmuItem));
      cart_assert(registry->Data != NULL);
    }
      
  /*
  //  Extend registry if needed
  */
  if(registry->NumItems == registry->Size)
    {
      registry->Size *= 2;
      tmp = (struct dmuItem *)malloc(registry->Size*sizeof(struct dmuItem));
      if(tmp == NULL)
	{
	  dmuPrintRegistryContents();
	  cart_error( "DMU: failure extendion allocation registry to %d items, %d bytes",registry->Size,registry->Size*sizeof(struct dmuItem));
	}
      memcpy(tmp,registry->Data,registry->NumItems*sizeof(struct dmuItem));
      free(registry->Data);
      registry->Data = tmp;
    }
}

      
void dmuRegister(void *ptr, unsigned long size, const char *file, int line)
{
  int i;

  if(size > dmuSmallLimit)
    {
      dmuCheckRegistry(&dmuLarge);

      dmuLarge.Data[dmuLarge.NumItems].Ptr = ptr;
      dmuLarge.Data[dmuLarge.NumItems].Size = size;
      dmuLarge.Data[dmuLarge.NumItems].File = file;
      dmuLarge.Data[dmuLarge.NumItems].Line = line;
      dmuLarge.NumItems++;
    }
  else
    {
      dmuCheckRegistry(&dmuSmall);

      if(size > dmuMaxSmallChunk) dmuMaxSmallChunk = size;

      for(i=0; i<dmuSmall.NumItems; i++)
	{
	  if(strcmp(dmuSmall.Data[i].File,file)==0 && dmuSmall.Data[i].Line==line)
	    {
	      break;
	    }
	}
      if(i == dmuSmall.NumItems)
	{
	  dmuSmall.Data[i].Ptr = NULL;
	  dmuSmall.Data[i].Size = 0UL;
	  dmuSmall.Data[i].File = file;
	  dmuSmall.Data[i].Line = line;
	  dmuSmall.NumItems++;
	}

      dmuSmall.Data[0].Size++;
    }

  dmuNumChunks[0]++;
  if(dmuNumChunks[0] > dmuNumChunks[1]) dmuNumChunks[1] = dmuNumChunks[0];
}


void dmuUnRegister(void *ptr, const char *file, int line)
{
  int i, j;
  
  /*
  //  Check large buffer registry first
  */
  for(i=0; i<dmuLarge.NumItems; i++)
    {
      if(dmuLarge.Data[i].Ptr == ptr) break;
    }

  if(i < dmuLarge.NumItems)
    {
      dmuLarge.NumItems--;
      for(j=i; j<dmuLarge.NumItems; j++)
	{
	  dmuLarge.Data[j] = dmuLarge.Data[j+1];
	}
    }
  else
    {
      /*
      //  Assume it was a small-size allocation
      */
      if(dmuSmall.Data[0].Size == 0)
	{
	  cart_error("DMU: freeing unregistered pointer %p, file %s, line %d",ptr,file,line);
	}
      else
	{
	  dmuSmall.Data[0].Size--;
	}
    }

  dmuNumChunks[0]--;
}


void dmuPrintRegistryContents()
{
  int i;
  unsigned long tot = dmuLarge.Size*sizeof(struct dmuItem) + dmuSmall.Size*sizeof(struct dmuItem);

  cart_debug("DMU: Large allocated chunks:");
  for(i=0; i<dmuLarge.NumItems; i++)
    {
      cart_debug("DMU: ptr=%p, size=%lu, file=%s, line=%d",dmuLarge.Data[i].Ptr,dmuLarge.Data[i].Size,dmuLarge.Data[i].File,dmuLarge.Data[i].Line);
      tot += dmuLarge.Data[i].Size;
    }

  if(dmuSmall.Data[0].Size > 0)
    {
      cart_debug("DMU: There are also %lu small allocated chunks (of size at most %d), created in files:",dmuSmall.Data[0].Size,dmuMaxSmallChunk);
      for(i=0; i<dmuSmall.NumItems; i++)
	{
	  cart_debug("DMU: file=%s, line=%d",dmuSmall.Data[i].File,dmuSmall.Data[i].Line);
	}
      cart_debug("DMU: total allocated memory does not exceed %lu",tot+dmuMaxSmallChunk*dmuSmall.Data[0].Size);
    }
  else
    {
      cart_debug("DMU: total allocated memory is %lu",tot);
    }
  cart_debug("DMU: current number of chunks is %lu",dmuNumChunks[0]);
  cart_debug("DMU: maximum number of chunks is %lu",dmuNumChunks[1]);
}


unsigned long dmuReportAllocatedMemory()
{
  int i;
  unsigned long tot = dmuLarge.Size*sizeof(struct dmuItem) + dmuSmall.Size*sizeof(struct dmuItem);

  for(i=0; i<dmuLarge.NumItems; i++)
    {
      tot += dmuLarge.Data[i].Size;
    }

  if(dmuSmall.NumItems > 0)
    {
      tot += dmuMaxSmallChunk*dmuSmall.Data[0].Size;
    }

  return tot;
}

#endif /* DEBUG_MEMORY_USE */
