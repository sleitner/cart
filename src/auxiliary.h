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

extern int num_options;
extern char **options;

const char* check_option0(const char* option, const char* name);
const char* check_option1(const char* option, const char* name, const char *default_value);

/*
//  Useful macros
*/
#ifndef min
#define min(x,y)        (((x) < (y)) ? (x): (y))
#endif
#ifndef max
#define max(x,y)        (((x) > (y)) ? (x): (y))
#endif
#define sign(x,y)       ( (y>=0) ? fabs(x) : -fabs(x) )

/*
//  Macros for memory leak locating
*/
#define cart_alloc_bytes(size)     cart_alloc_worker(size,__FILE__,__LINE__)
#define cart_alloc(type,size)      (type *)cart_alloc_worker((size)*sizeof(type),__FILE__,__LINE__)   /* A properly cast version for strict type checking */
#define cart_free(ptr)             { cart_free_worker(ptr,__FILE__,__LINE__); ptr = NULL; }

void *cart_alloc_worker(size_t size, const char *file, int line);
void cart_free_worker(void *ptr, const char *file, int line);

#ifdef DEBUG_MEMORY_USE
void dmuPrintRegistryContents();
unsigned long dmuReportAllocatedMemory();
#endif

#endif
