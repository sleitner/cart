#ifndef __INDEX_HASH_H__
#define __INDEX_HASH_H__

#include <sys/types.h>
                                                                                                                                              
typedef struct INDEX_HASH_ENTRY {
	int local_index;
	int remote_index;
} index_hash_entry;

typedef struct INDEX_HASH {
	int hash_size;
	int64_t a0, b0;
	int64_t p;
	int *n, *s;
	int64_t *a, *b;
    index_hash_entry *hash_array;
} index_hash;

index_hash *index_hash_create( int size, int64_t max_key, int *remote_index, int *local_index );
void index_hash_free( index_hash *hash );
int index_hash_lookup ( index_hash *hash, int remote_index );

#endif /* __INDEX_HASH_H__ */
