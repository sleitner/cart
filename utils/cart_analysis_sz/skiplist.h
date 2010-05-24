#ifndef __SKIPLIST_H__
#define __SKIPLIST_H__

#define MAX_SKIPLIST_LEVEL	16

typedef struct SKIPLIST_NODE {
	int value;
	struct SKIPLIST_NODE *next[1];
} skiplist_node;

typedef struct SKIPLIST {
	int num_nodes;
	int level;
	int random_buffer;
	int num_random;
	skiplist_node *head;
	skiplist_node *ptr;
} skiplist;

skiplist *skiplist_init();
int skiplist_size( skiplist *list );
void skiplist_destroy( skiplist *list );
int skiplist_insert( skiplist *list, int value );
int skiplist_contains( skiplist *list, int value );
void skiplist_iterate( skiplist *list );
int skiplist_next( skiplist *list, int *value );

#endif
