#ifndef __MIN_HEAP_H__
#define __MIN_HEAP_H__

typedef struct MIN_HEAP {
	int size;
	int cursize;
	int *keys;
	double *values;
} min_heap;

min_heap *min_heap_init( size_t min_heap_size );
void min_heap_destroy( min_heap * );
void min_heap_flush( min_heap * );
int min_heap_size( min_heap * );
void min_heap_print( min_heap *h );
void min_heap_push( min_heap *, int key, double value );
int min_heap_peek( min_heap *, int *key, double *value );
int min_heap_pop( min_heap *, int *key, double *value );
void min_heap_reform( min_heap *h );

#endif /* __min_heap_H__ */
