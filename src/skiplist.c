#include <stdlib.h>
#include <stdio.h>

#include "auxiliary.h"
#include "skiplist.h"

skiplist *skiplist_init() {
	int i;
	skiplist *list;
	skiplist_page *page;

	list = cart_alloc(skiplist, 1 );

	list->num_nodes = 0;
	list->level = 0;
	page = cart_alloc(skiplist_page, 1 );
	list->page_head = page;
	
	list->head = (skiplist_node *)page->buffer;
	page->current_ptr = (void *)((skiplist_node **)((skiplist_node *)page->buffer+1)+(MAX_SKIPLIST_LEVEL));
	page->bytes_remaining = SKIPLIST_PAGE_SIZE*sizeof(int) - sizeof(skiplist_node) - MAX_SKIPLIST_LEVEL*sizeof(skiplist *);
	page->next = NULL;
	
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
	skiplist_page *n, *p;

	cart_assert( list != NULL );

	n = list->page_head;
	while ( n != NULL ) {
		p = n;
		n = n->next;
		cart_free(p);
	}

	cart_free( list );
}

int skiplist_insert( skiplist *list, int value ) {
	skiplist_node *n, *p;
	skiplist_node *update[MAX_SKIPLIST_LEVEL+1];
	skiplist_page *page;
	int level;
	int i, r;
	size_t size;

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

	cart_assert( level <= MAX_SKIPLIST_LEVEL );

	if ( level > list->level ) {
		for ( i = list->level + 1; i <= level; i++ ) {
			update[i] = list->head;
		}

		list->level = level;
	}

	/* allocate the new node and insert into the list */
	size = sizeof(skiplist_node)+level*sizeof(skiplist_node *);
	if ( list->page_head->bytes_remaining < size ) {
		page = cart_alloc(skiplist_page, 1 );
		page->current_ptr = page->buffer;
		page->bytes_remaining = SKIPLIST_PAGE_SIZE*sizeof(int);
		page->next = list->page_head;
		list->page_head = page;
	} else {
		page = list->page_head;
	}

	n = (skiplist_node *)page->current_ptr;
	page->current_ptr = (void *)((skiplist_node **)((skiplist_node *)page->current_ptr+1)+level);
	page->bytes_remaining -= size;
	cart_assert( page->bytes_remaining >= 0 );

	n->value = value;
	for ( i = 0; i <= level; i++ ) {
		p = update[i];
		n->next[i] = p->next[i];
		p->next[i] = n;
	}

	list->num_nodes++;
	return 1;
}

int skiplist_contains( skiplist *list, int value ) {
	skiplist_node *n, *p;
	int level;

	n = list->head;
	level = list->level;

	do {
		/* find position at this level */
		p = n->next[level];
		while ( p != NULL && p->value <= value ) {
			n = p;
			p = p->next[level];
		}
		level--;
	} while ( level >= 0 );

	/* silently drop duplicate values */
	if ( n->value == value ) {
		return 1;
        } else {
                return 0;
        }
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

	cart_assert( list != NULL );

	if ( list->ptr == NULL ) {
		return 0;
	} else {
		*value = list->ptr->value;
		list->ptr = list->ptr->next[0];
		return 1;
	}
}
