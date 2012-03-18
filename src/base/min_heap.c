
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "min_heap.h"

min_heap *min_heap_init( size_t heap_size ) {
	min_heap *newheap;

	newheap = cart_alloc(min_heap,1);
	newheap->keys = cart_alloc( int, heap_size );
	newheap->values = cart_alloc( double, heap_size );

	newheap->size = heap_size;
	newheap->cursize = 0;

	return newheap;
}

void min_heap_destroy( min_heap *h ) {
	cart_free( h->keys );
	cart_free( h->values );
	cart_free( h );
}

void min_heap_flush( min_heap *h ) {
	cart_assert( h != NULL );
	h->cursize = 0;
}

int min_heap_size( min_heap *h ) {
	return h->cursize;
}

void min_heap_push( min_heap *h, int key, double value ) {
	int i;

	if ( h->cursize == h->size ) {
		if ( h->values[0] > value ) {
			/* replace extreme element and reheapify */
			h->keys[0] = key;	
			h->values[0] = value;
            min_heap_reform(h);	
		}
		return;
	}

	i = h->cursize;
	h->cursize++;

	while ( i > 0 && h->values[(i-1)>>1] < value ) {
		h->keys[i] = h->keys[(i-1)>>1];
		h->values[i] = h->values[(i-1)>>1];
		i = (i-1)>>1;
	}

	h->keys[i] = key;
	h->values[i] = value;
}

int min_heap_peek( min_heap *h, int *key, double *value ) {
	if ( h->cursize == 0 ) {
		return 0;
	} else {
		*key = h->keys[0];
		*value = h->values[0];
	}

	return 1;
}

void min_heap_print( min_heap *h ) {
	int i;

	for ( i = 0; i < h->cursize; i++ ) {
		cart_debug("%u %d %e", i, h->keys[i], h->values[i] );
	}
}

int min_heap_pop( min_heap *h, int *key, double *value ) {
	if ( h->cursize == 0 ) {
		return 0;
	} else {
		*key = h->keys[0];
		*value = h->values[0];

		h->cursize--;

		if ( h->cursize > 0 ) {
			h->keys[0] = h->keys[h->cursize];
			h->values[0] = h->values[h->cursize];
			min_heap_reform(h);
		}
	}

	return 1;
}

void min_heap_reform( min_heap *h ) {
	int i, left, right, extremum;
	int temp_key;
	double temp_value;

	i = 0; 
	while (1) {
		left = 2*i+1;
		right = left+1;

		if ( left < h->cursize && h->values[left] > h->values[i] ) {
			extremum = left;
		} else {
			extremum = i;
		}

		if ( right < h->cursize && h->values[right] > h->values[extremum] ) {
			extremum = right;
		}

		if ( extremum != i ) {
			/* extremum -> tmp */
			temp_key = h->keys[extremum];
			temp_value = h->values[extremum];

			/* i -> extremum */
			h->keys[extremum] = h->keys[i];
			h->values[extremum] = h->values[i];

			/* tmp -> i */
			h->keys[i] = temp_key;
			h->values[i] = temp_value;
			i = extremum;
		} else {
			break;
		}
	}
}
