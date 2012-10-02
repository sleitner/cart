#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>

#include "auxiliary.h"
#include "index_hash.h"
#include "rand.h"

const index_hash_entry null_entry = { -1, -1 };

/*******************************************************
 * next_largest_prime
 *******************************************************/
int64_t next_largest_prime64( int64_t count )
/* purpose: calculates the next largest prime by brute 
 * 	force, used since primes make good hash sizes
 * 	(keeps the number of hash collisions low)
 * returns: the next prime larger than count, or 3 if count <= 2
 */
{
	int64_t i;
	int is_prime;
	int64_t next_prime;

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
 * index_hash_create 
 *******************************************************/
index_hash *index_hash_create( int size, int64_t max_key, int *remote_index, int *local_index ) 
/* purpose: allocates space and initializes a new 
 * 	index hash which will contain at most size
 * 	entries
 *
 * returns: a pointer to the newly allocated index hash
 */
{
	int i, j, k, h;
	int oldj;
	index_hash *hash;
	int iter = 0;
	size_t total;
	int rekey;
	index_hash_entry tmp_entry, new_entry;
	int64_t a0, b0, min_a0, min_b0;
	size_t min_total;

	/* allocate space for the new hash */
	hash = cart_alloc( index_hash, 1 );

	/* create universal hash function for first layer */
	hash->hash_size = size;
	hash->p = next_largest_prime64( 3*max_key );

	hash->a = cart_alloc( int64_t, hash->hash_size );
	hash->b = cart_alloc( int64_t, hash->hash_size );
	hash->s = cart_alloc( int64_t, hash->hash_size );
	hash->n = cart_alloc( int, hash->hash_size );

	min_total = SIZE_MAX;

	do { 
		/* choose a first-level hash function */
		a0 = (int64_t)(hash->p*cart_rand()) + 1;
		b0 = (int64_t)(hash->p*cart_rand());
		
		for ( i = 0; i < size; i++ ) {
			hash->n[i] = 0;
		}

		for ( i = 0; i < size; i++ ) {
			h = (((int64_t)remote_index[i]*a0 + b0) % hash->p) % hash->hash_size;
			hash->n[h]++;
		}

		total = 0;
		for ( i = 0; i < size; i++ ) {
			hash->n[i] *= hash->n[i];
			total += hash->n[i];
		}

		if ( total < min_total ) {
			min_total = total;
			min_a0 = a0;
			min_b0 = b0;
		}

		iter++;
	} while ( min_total >= 2*size && iter < 20 );

	hash->a0 = min_a0;
	hash->b0 = min_b0;

	total = 0;
	for ( i = 0; i < size; i++ ) {
		hash->s[i] = total;
		total += hash->n[i];
	}

	/* allocate space for the actual hash array */
	hash->hash_array = cart_alloc( index_hash_entry, total );

	for ( i = 0; i < total; i++ ) {
		hash->hash_array[i] = null_entry;
	}

	for ( i = 0; i < size; i++ ) {
		hash->n[i] = 0;
	}

	for ( i = 0; i < size; i++ ) {
		h = (((int64_t)remote_index[i]*hash->a0 + hash->b0) % hash->p) % hash->hash_size;
		hash->hash_array[ hash->s[h] + hash->n[h] ].remote_index = remote_index[i];
		hash->hash_array[ hash->s[h] + hash->n[h] ].local_index = local_index[i];
		hash->n[h]++;
	}

	for ( i = 0; i < size; i++ ) {
		hash->n[i] *= hash->n[i];
	}	
	
	/* hash each secondary array */
	#pragma omp parallel for default(none) shared(hash,null_entry,size) private(i,j,k,oldj,rekey,tmp_entry,new_entry) schedule(dynamic)
	for ( i = 0; i < size; i++ ) {
		do {
			rekey = 0;
	
			/* choose secondary hash function */
			hash->a[i] = (int64_t)(hash->p*cart_rand()) + 1;
			hash->b[i] = (int64_t)(hash->p*cart_rand());

			for ( j = 0; j < hash->n[i]; j++ ) {
				if ( hash->hash_array[ hash->s[i] + j ].remote_index != -1 ) {
					/* find new hash location */
					new_entry = hash->hash_array[ hash->s[i] + j ];
					hash->hash_array[ hash->s[i] + j ] = null_entry;
					oldj = j;

					while (1) {
						k = ((hash->a[i]*(int64_t)new_entry.remote_index + hash->b[i]) % hash->p) % hash->n[i];
						if ( hash->hash_array[ hash->s[i] + k ].remote_index == -1 ) {
							hash->hash_array[ hash->s[i] + k ] = new_entry;
							break;
						} else if ( k > oldj ) {
							tmp_entry = hash->hash_array[ hash->s[i] + k ];
							hash->hash_array[ hash->s[i] + k ] = new_entry;
							new_entry = tmp_entry;
							oldj = k;
						} else {
							/* unresolvable collision, need to rekey, put new_entry in first open spot */
							for ( k = 0; k < hash->n[i]; k++ ) {
								if ( hash->hash_array[ hash->s[i] + k ].remote_index == -1 ) {
									hash->hash_array[ hash->s[i] + k ] = new_entry;
									break;
								}
							}
							rekey = 1;
							break;	
						}
					}
				}

				if ( rekey ) {
					break;
				}
			}
		} while ( j < hash->n[i] );
	}

	return hash;
}

/*******************************************************
 * index_hash_lookup
 *******************************************************/
int index_hash_lookup ( index_hash *hash, int remote_index ) 
/* purpose: finds the local_index corresponding to the
 * 	given remote index in the given index hash.
 *
 * returns: the corresponding local index or -1 if the
 * 	remote index was not found.
 */
{
	int h, k;

	if ( remote_index != INDEX_HASH_NULL_ENTRY ) {
		h = (((int64_t)remote_index*hash->a0 + hash->b0) % hash->p) % hash->hash_size;      

		if ( hash->n[h] > 0 ) {
			k = ((hash->a[h]*(int64_t)remote_index+hash->b[h]) % hash->p) % hash->n[h];
			if ( hash->hash_array[ hash->s[h] + k ].remote_index == remote_index ) {
				return hash->hash_array[ hash->s[h] + k ].local_index;
			}
		}
	}

	return INDEX_HASH_NULL_ENTRY;
}

void index_hash_free( index_hash *hash ) {
	cart_assert( hash != NULL );

	if ( hash->hash_size > 0 ) {
		cart_free( hash->n );
		cart_free( hash->s );
		cart_free( hash->a );
		cart_free( hash->b );
		cart_free( hash->hash_array );
	}

	cart_free( hash );
}
