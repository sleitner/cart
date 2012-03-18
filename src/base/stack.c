
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "stack.h"

stack *stack_init() {
	stack *newstack;
	stack_page *page;

	newstack = cart_alloc(stack,1);
	page = cart_alloc(stack_page,1);

	page->index = 0;
    page->next = NULL;

	newstack->head = page;
	newstack->size = 0;

	return newstack;
}

void stack_destroy( stack *s ) {
	stack_page *tmp;
	stack_page *page = s->head;

	while ( page != NULL ) {
		tmp = page->next;
		cart_free( page );
		page = tmp;
	}

	cart_free( s );
}

void stack_print( stack *s ) {
	int i;
	stack_page *p;

	for ( i = s->head->index-1; i >= 0; i-- ) {
		cart_debug("%d", s->head->data[i] );
	}

	p = s->head->next;	
	while ( p != NULL ) {
		for ( i = p->index-1; i >= 0; i-- )  {
			cart_debug("%d", p->data[i] );
		}
		p = p->next;
	}
}

void stack_push( stack *s, int value ) {
	stack_page *newpage;
	stack_page *page = s->head;

	if ( page->index == STACK_PAGE_SIZE ) {
		newpage = cart_alloc(stack_page,1);

		newpage->index = 1;
		newpage->data[0] = value;
		newpage->next = s->head;
		s->head = newpage;
	} else {
		page->data[page->index++] = value;
	}

	s->size++;
}

int stack_pop( stack *s, int *value ) {
	stack_page *page;

	page = s->head;
	if ( page->index == 0 ) {
		return 0;
	} else {
		*value = page->data[--page->index];
		s->size--;

		if ( page->index == 0 && page->next != NULL ) {
			s->head = page->next;
			cart_free(page);
		}

		return 1;
	}
}
