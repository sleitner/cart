#ifndef __MEMORY_H__
#define __MEMORY_H__


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
size_t dmuReportAllocatedMemory();
#endif

#endif
