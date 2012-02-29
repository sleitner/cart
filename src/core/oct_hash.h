#ifndef __OCT_HASH_H__
#define __OCT_HASH_H__

typedef struct OCT_HASH_ENTRY {
	int local_index;
	int remote_index;
} oct_hash_entry;

typedef struct OCT_HASH {
	int hash_size;
	int num_entries;
	oct_hash_entry *hash_array;
} oct_hash;

#define DELETED_ENTRY		-2

oct_hash *oct_hash_create( int size );
void oct_hash_free( oct_hash *hash );
void oct_hash_add( oct_hash *hash, int remote_index, int local_index );
void oct_hash_add_list( oct_hash *hash, int num_indices, int *remote_index, int *local_index );
void oct_hash_delete( oct_hash *hash, int local_index );
int oct_hash_lookup ( oct_hash *hash, int remote_index );

#endif
