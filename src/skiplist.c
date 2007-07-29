#include <stdlib.h>
#include <stdio.h>

#include "auxiliary.h"
#include "skiplist.h"

skiplist *skiplist_init() {
	int i;
	skiplist *list;

	list = (skiplist *)cart_alloc( sizeof( skiplist ) );

	list->num_nodes = 0;
	list->level = 0;
	list->head = (skiplist_node *)cart_alloc( sizeof(skiplist_node) + MAX_SKIPLIST_LEVEL*sizeof(skiplist_node *) );
	list->head->value = -1;
	for ( i = 0; i <= MAX_SKIPLIST_LEVEL; i++ ) {
		list->head->next[i] = NULL;
	}
	list->ptr = NULL;

	list->random_buffer = rand();
	list->num_random = sizeof(int)*4;

	return list;
}

void skiplist_destroy( skiplist *list ) {
	skiplist_node *n, *p;

	cart_assert( list != NULL );

	n = list->head;
	while ( n != NULL ) {
		p = n;
		n = n->next[0];
		cart_free(p);
	}

	cart_free( list );
}

int skiplist_insert( skiplist *list, int value ) {
	skiplist_node *n, *p;
	skiplist_node *update[MAX_SKIPLIST_LEVEL+1];
	int level;
	int i, r;

	cart_assert( list != NULL );
	cart_assert( value >= 0 );

	n = list->head;
	level = list->level;

	do {
		/* find position at this level */
		p = n->next[level];
		while ( p != NULL && p->value <= value ) {
			n = p;
			p = p->next[level];
		}
		update[level] = n;
		level--;
	} while ( level >= 0 );
	
	/* silently drop duplicate values */
	if ( n->value == value ) {
		return 0;
	}

	/* generate a random level for the new node */
	level = 0;
	do {
		r = list->random_buffer & 3;
		list->random_buffer >>= 2;
		list->num_random--;
		if ( list->num_random == 0 ) {
			list->random_buffer = rand();
			list->num_random = sizeof(int)*4;
		}
		
		if ( !r ) {
			level++;
		}
	} while ( !r && level < MAX_SKIPLIST_LEVEL );

	cart_assert( level < MAX_SKIPLIST_LEVEL );

	if ( level > list->level ) {
		for ( i = list->level + 1; i <= level; i++ ) {
			update[i] = list->head;
		}

		list->level = level;
	}

	/* allocate the new node and insert into the list */
	n = (skiplist_node *)cart_alloc( sizeof(skiplist_node) + level*sizeof(skiplist_node *) );
	n->value = value;
	for ( i = 0; i <= level; i++ ) {
		p = update[i];
		n->next[i] = p->next[i];
		p->next[i] = n;
	}

	list->num_nodes++;
	return 1;
}

int skiplist_size( skiplist *list ) {
	return list->num_nodes;
}

void skiplist_iterate( skiplist *list ) {
	cart_assert( list != NULL );

	/* skip over first blank node */
	list->ptr = list->head->next[0];
}

int skiplist_next( skiplist *list, int *value ) {
	int val;

	cart_assert( list != NULL );

	if ( list->ptr == NULL ) {
		return 0;
	} else {
		*value = list->ptr->value;
		list->ptr = list->ptr->next[0];
		return 1;
	}
}
