#ifndef __AUXILIARY_H__
#define __AUXILIARY_H__

#include <stdlib.h>


extern int num_procs;
extern int local_proc_id;


double integrate( double (*f)(double), double a, double b, double epsrel, double epsabs );
double root_finder( double (*f)(double), double a, double b, double epsrel, double epsabs );
int compare_ints( const void *a, const void *b );
int compare_floats( const void *a, const void *b );
int nearest_int( double x );


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

extern int num_options;
extern const char **options;

int is_option_present(const char* full_name, const char* short_name, int with_argument);
const char* extract_option0(const char* full_name, const char* short_name);
const char* extract_option1(const char* full_name, const char* short_name, const char *default_value);
void die_on_unknown_options();

/* 
// Endian helper function
*/
void reorder( char *buffer, int size );

/*
// Helper functions
*/
void linear_array_maxmin(int n, float *arr, float *max, float *min);
void linear_array_copy_int(int *dest, int *src, int size);
void linear_array_copy_float(float *dest, float *src, int size);


/*
//  Useful macros
*/
#ifndef MIN
#define MIN(x,y)        (((x) < (y)) ? (x): (y))
#endif
#ifndef MAX
#define MAX(x,y)        (((x) > (y)) ? (x): (y))
#endif
#define SIGN(x,y)       ( (y>=0) ? fabs(x) : -fabs(x) )

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
size_t dmuReportAllocatedMemory();
#endif

#endif
