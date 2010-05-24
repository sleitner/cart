#ifndef __AUXILIARY_H__
#define __AUXILIARY_H__

#include <stdlib.h>
#include <stdarg.h>

#include "sfc.h"

#define min(x,y)        (((x) < (y)) ? (x): (y))
#define max(x,y)        (((x) > (y)) ? (x): (y))
#define sign(x,y)       ( (y>=0) ? fabs(x) : -fabs(x) )

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


#endif
