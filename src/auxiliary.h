#ifndef __AUXILIARY_H__
#define __AUXILIARY_H__

#include <stdlib.h>
#include <stdarg.h>

void init_auxiliary();
void save_auxiliary();

extern unsigned long int rng_seed;

double integrate( double (*f)(double), double a, double b, double epsrel, double epsabs );
double root_finder( double (*f)(double), double a, double b, double epsrel, double epsabs );
int compare_ints( const void *a, const void *b );
int compare_floats( const void *a, const void *b );
int nearest_int( double x );

double cart_rand();

void *cart_alloc( size_t size );
void cart_free( void *ptr );

void cart_error( const char *fmt, ... );

#ifndef NDEBUG
void cart_debug( const char *fmt, ... );
#else
#define cart_debug( ... )   
#endif

#ifdef __STDC__
#ifndef NDEBUG
#define cart_assert( x ) if (!(x)) { cart_error( "Assertion (%s) failed: %s line %u", #x, __FILE__, __LINE__ ); }
#else
#define cart_assert( x ) 
#endif
#else
#define cart_assert( x ) assert(x)
#endif

typedef struct {
        int num_eqn;
        double *y0;
	double *y1;
        double *rs; 
        double *a0; 
        double *a1;
        double *w0;
        double *w1;
        double *buf;
        void (* rates)(double t, double *y, void *params, double *w, double *a);
        void (* adjust)(double t, double *y, void *params);
} qss_system;

qss_system *qss_alloc( size_t num_eqn,  
		void (* rates)(double, double *, void *, double *, double *),
		void (* adjust)(double, double *, void *) ); 
void qss_free( qss_system *sys );
void qs1_step( qss_system *sys, double t, double dt, double yf[], void *params );
void qsn_step( qss_system *sys, double t, double dt, double yf[], void *params );
void qss_solve( qss_system *sys, double t_begin, double delta_t, double y[], const double err[], void *params );

#endif
