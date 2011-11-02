#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "oct_hash.h"
#include "parallel.h"
#include "tree.h"


/*******************************************************
 * next_largest_prime
 *******************************************************/
int next_largest_prime( int count )
/* purpose: calculates the next largest prime by brute 
 * 	force, used since primes make good hash sizes
 * 	(keeps the number of hash collisions low)
 * returns: the next prime larger than count, or 3 if count <= 2
 */
{
	int i;
	int is_prime;
	int next_prime;

	if ( count > 2 ) {
		next_prime = ( count % 2 == 0 ) ? count + 1: count + 2;
	} else {
		next_prime = 3;
	}
	
	do {
		is_prime = 1;

		for ( i = 3; i <= sqrt(next_prime); i += 2 ) {
			if ( next_prime % i == 0 ) {
				is_prime = 0;
				next_prime += 2;
				break;
			}
		}
	} while ( !is_prime );

	return next_prime;
}

/*******************************************************
 * oct_hash_create 
 *******************************************************/
oct_hash *oct_hash_create( int size ) 
/* purpose: allocates space and initializes a new 
 * 	index hash which will contain at most size
 * 	entries
 *
 * returns: a pointer to the newly allocated oct hash
 */
{
	int i;
	oct_hash *hash;

	/* allocate space for the new hash */
	hash = cart_alloc(oct_hash, 1 );
	
	/* select a prime size larger than 2*size (leaving room for refinement) */
	hash->hash_size = next_largest_prime( 2*size );
	hash->num_entries = 0;
	
	/* allocate space for the actual hash array */
	hash->hash_array = cart_alloc(oct_hash_entry, hash->hash_size );

	/* initialize the hash array to empty entries */
	for ( i = 0; i < hash->hash_size; i++ ) {
		hash->hash_array[i].remote_index = NULL_OCT;
		hash->hash_array[i].local_index = NULL_OCT;
	}

	return hash;
}

/*******************************************************
 * oct_hash_free
 *******************************************************/
void oct_hash_free( oct_hash *hash ) 
/* purpose: deallocates the given oct hash
 */
{
	cart_assert( hash != NULL );
	
	if ( hash->hash_size > 0 ) {
		cart_assert( hash->hash_array != NULL );
		cart_free( hash->hash_array );
	}

	cart_free( hash );
}

void oct_hash_resize( oct_hash *hash, int new_size ) {
	int i;
	int count;
	int old_entries;
	oct_hash_entry *old_hash;
	int *old_remote, *old_local;

	cart_assert( hash->hash_size < new_size );

	old_remote = cart_alloc(int, hash->num_entries );
	old_local = cart_alloc(int, hash->num_entries );
	old_hash = hash->hash_array;

	count = 0;
	for ( i = 0; i < hash->hash_size; i++ ) {
		if ( old_hash[i].remote_index != NULL_OCT &&
				old_hash[i].remote_index != DELETED_ENTRY ) {
			old_remote[count] = old_hash[i].remote_index;
			old_local[count] = old_hash[i].local_index;
			count++;
		}
	}

	cart_assert( count == hash->num_entries );
	cart_free( hash->hash_array );

	hash->hash_size = next_largest_prime( new_size );
	hash->hash_array = cart_alloc(oct_hash_entry, hash->hash_size );

	/* initialize the hash array to empty entries */
	for ( i = 0; i < hash->hash_size; i++ ) {
		hash->hash_array[i].remote_index = NULL_OCT;
		hash->hash_array[i].local_index = NULL_OCT;
	}

	old_entries = hash->num_entries;
	hash->num_entries = 0;

	/* re-add all the old entries */
	oct_hash_add_list( hash, old_entries, old_remote, old_local );

	count = 0;
	for ( i = 0; i < hash->hash_size; i++ ) {
		if ( hash->hash_array[i].remote_index != DELETED_ENTRY &&
			hash->hash_array[i].remote_index != NULL_OCT ) {
			count++;
		}
	}
	
	cart_assert( count == hash->num_entries );
	cart_assert( hash->num_entries == old_entries );

	/* free arrays */
	cart_free( old_remote );
	cart_free( old_local );
}

/*******************************************************
 * oct_hash_add
 ********************************************************/
void oct_hash_add( oct_hash *hash, int remote_index, int local_index )	
/* purpose: adds the given remote_index->local_index mapping to
 * 	the given oct hash
 *
 * Uses Brent's method for reducing probing times
 * Brent, Richard.  Reducing the Retrieval Time of Scatter Storage Techniques.
 * CACM, Vol. 16 #2.  1973.
 * -or-
 * Knuth.  The Art of Computer Programming Vol 3 pg 532-533
 */
{
	int c, h, i, j, p, q, r, s, t;
	int placed;

	cart_assert( hash != NULL );
	cart_assert( remote_index >= 0 );
	cart_assert( local_index >= 0 );

	/* we may later want to lower the threshold for increasing
	 * hash size to keep the number of probes down */
	if ( hash->num_entries >= hash->hash_size ) {
		oct_hash_resize( hash, 2*hash->hash_size );
	}

	cart_assert( hash->num_entries < hash->hash_size );

	r = remote_index % hash->hash_size;
	q = (remote_index % (hash->hash_size-2)) + 1;
	h = r;

	/* probe for the proper hash position */
	for ( s = 0; s < hash->hash_size; s++ ) {
		cart_assert( h >= 0 && h < hash->hash_size );
		cart_assert( hash->hash_array[h].remote_index != remote_index );

		if ( hash->hash_array[h].remote_index == NULL_OCT ||
				hash->hash_array[h].remote_index == DELETED_ENTRY ) {
			break;
		}
		h = ( h + q ) % hash->hash_size;
	}

	cart_assert( s < hash->hash_size );

	/* look for a better place to put this entry to reduce number of probes */
	placed = 0;
	for ( i = 0; i < s - 1 && !placed; i++ ) {
		p = ( r + i*q) % hash->hash_size;
		c = (hash->hash_array[p].remote_index % ( hash->hash_size - 2 ) ) + 1;
		
		for ( j = 1; j < s - i - 1; j++ ) {
			t = ( p + j * c) % hash->hash_size;
			
			/* if this is a better place, switch entry at p to t and place in p */
			if ( hash->hash_array[t].remote_index == NULL_OCT || 
					hash->hash_array[t].remote_index == DELETED_ENTRY ) {
				hash->hash_array[t].remote_index = hash->hash_array[p].remote_index;
				hash->hash_array[t].local_index = hash->hash_array[p].local_index;
				hash->hash_array[p].remote_index = remote_index;
				hash->hash_array[p].local_index = local_index;

				placed = 1;
				break;
			}
		}
	}

	/* we couldn't find a better place to put this entry, so place it at h */
	if ( !placed ) {
		hash->hash_array[h].remote_index = remote_index;
		hash->hash_array[h].local_index = local_index;
	}

	hash->num_entries++;

	placed = 0;
	for ( i = 0; i < hash->hash_size; i++ ) {
		if ( hash->hash_array[i].remote_index != NULL_OCT &&
			hash->hash_array[i].remote_index != DELETED_ENTRY ) {

			placed++;
		}
	}
	cart_assert( placed == hash->num_entries );	
}

/*******************************************************
 * oct_hash_add_list
 ********************************************************/
void oct_hash_add_list( oct_hash *hash, int num_indices, int *remote_index, int *local_index )
/* purpose: adds the given list of remote_index->local_index mappings to
 *      the given oct hash
 *      resizes the hash 1 time (major savings when we're adding large numbers of indices)
 */
{
	int c, h, i, j, p, q, r, s, t;
	int index;
	int placed;

	cart_assert( hash != NULL );

	/* we may later want to lower the threshold for increasing
	 * hash size to keep the number of probes down */
	if ( hash->num_entries+num_indices >= hash->hash_size ) {
		oct_hash_resize( hash, 2*hash->hash_size+num_indices );
	}

	cart_assert( hash->num_entries < hash->hash_size );

	for ( index = 0; index < num_indices; index++ ) {
		cart_assert( remote_index[index] >= 0 );
		cart_assert( local_index[index] >= 0 );

		r = remote_index[index] % hash->hash_size;
		q = (remote_index[index] % (hash->hash_size-2)) + 1;
		h = r;

		/* probe for the proper hash position */
		for ( s = 0; s < hash->hash_size; s++ ) {
			cart_assert( h >= 0 && h < hash->hash_size );

			/* ensure no duplicate entries */
			cart_assert( hash->hash_array[h].remote_index != remote_index[index] );

			if ( hash->hash_array[h].remote_index == NULL_OCT ||
					hash->hash_array[h].remote_index == DELETED_ENTRY ) {
				break;
			}
			h = ( h + q ) % hash->hash_size;
		}

		cart_assert( s < hash->hash_size );

		/* look for a better place to put this entry to reduce number of probes */
		placed = 0;
		for ( i = 0; i < s - 1 && !placed; i++ ) {
			p = ( r + i*q) % hash->hash_size;
			c = (hash->hash_array[p].remote_index % ( hash->hash_size - 2 ) ) + 1;

			for ( j = 1; j < s - i - 1; j++ ) {
				t = ( (long)p + (long)j*(long)c) % hash->hash_size;
				cart_assert( t >= 0 && t < hash->hash_size );

				/* if this is a better place, switch entry at p to t and place in p */
				if ( hash->hash_array[t].remote_index == NULL_OCT ||
						hash->hash_array[t].remote_index == DELETED_ENTRY ) {

					hash->hash_array[t].remote_index = hash->hash_array[p].remote_index;
					hash->hash_array[t].local_index = hash->hash_array[p].local_index;
					hash->hash_array[p].remote_index = remote_index[index];
					hash->hash_array[p].local_index = local_index[index];
					hash->num_entries++;

					placed = 1;
					break;
				}
			}
		}

		/* we couldn't find a better place to put this entry, so place it at h */
		if ( !placed ) {
			hash->hash_array[h].remote_index = remote_index[index];
			hash->hash_array[h].local_index = local_index[index];
			hash->num_entries++;
		}
	}
}

/*******************************************************
 * oct_hash_delete
 *******************************************************/
void oct_hash_delete( oct_hash *hash, int local_index ) 
/* purpose: removes the given local_index from the oct 
 * 	hash.
 */
{
	int i,j,count, deleted;
	
	cart_assert( hash != NULL );

	count = 0;
	deleted = 0;
	for ( j = 0; j < hash->hash_size; j++ ) {
		if ( hash->hash_array[j].remote_index != DELETED_ENTRY &&
				hash->hash_array[j].remote_index != NULL_OCT ) {
			count++;
		}

		if ( hash->hash_array[j].remote_index == DELETED_ENTRY ) {
			deleted++;
		}
	}

	if ( count != hash->num_entries ) {
		cart_debug("count = %u", count );
		cart_debug("deleted = %d", deleted );
		cart_debug("hash->num_entries = %d", hash->num_entries );
		cart_debug("hash->hash_size = %d", hash->hash_size );
	}

	cart_assert( count == hash->num_entries );
	
	for ( i = 0; i < hash->hash_size; i++ ) {
		if ( hash->hash_array[i].local_index == local_index ) {
			hash->hash_array[i].remote_index = DELETED_ENTRY;
			hash->hash_array[i].local_index = NULL_OCT;
			hash->num_entries--;

			count = 0;
			for ( j = 0; j < hash->hash_size; j++ ) {
				if ( hash->hash_array[j].remote_index != DELETED_ENTRY &&
					hash->hash_array[j].remote_index != NULL_OCT ) {
					count++;
				}
			}

			cart_assert( count == hash->num_entries );

			return;
		}
	}

	cart_error("oct_hash_delete unable to find %u", local_index );
}

/*******************************************************
 * oct_hash_lookup
 *******************************************************/
int oct_hash_lookup ( oct_hash *hash, int remote_index ) 
/* purpose: finds the local_index corresponding to the
 * 	given remote index in the given oct hash.
 *
 * returns: the corresponding local index or -1 if the
 * 	remote index was not found.
 */
{
	int h, q, i;

	if ( remote_index != NULL_OCT ) {
		h = remote_index % hash->hash_size;
		q = (remote_index % (hash->hash_size - 2)) + 1;

		for ( i = 0; i < hash->hash_size; i++ ) {
			if ( hash->hash_array[h].remote_index == NULL_OCT ) {
				return NULL_OCT;
			} else if ( hash->hash_array[h].remote_index == remote_index ) {
				return hash->hash_array[h].local_index;
			} else {
				/* note we're skipping over DELETED_ENTRY spaces */
				h = (h+q) % hash->hash_size;
			}
		}
	}

	return NULL_OCT;
}
