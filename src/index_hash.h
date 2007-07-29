#ifndef __INDEX_HASH_H__
#define __INDEX_HASH_H__

typedef struct INDEX_HASH_ENTRY {
	int local_index;
	int remote_index;
} index_hash_entry;

typedef struct INDEX_HASH {
	int hash_size;
	int num_entries;
	index_hash_entry *hash_array;
} index_hash;

#define DELETED_ENTRY		-2

index_hash *index_hash_create( int size );
void index_hash_free( index_hash *hash );
void index_hash_add( index_hash *hash, int remote_index, int local_index );
void index_hash_add_list( index_hash *hash, int num_indices, int *remote_index, int *local_index );
void index_hash_delete( index_hash *hash, int local_index );
int index_hash_lookup ( index_hash *hash, int remote_index );

#endif
