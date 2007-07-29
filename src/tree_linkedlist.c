#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "tree.h"
#include "tree_linkedlist.h"
#include "auxiliary.h"

void linked_list_insert( int *head, int ioct ) {
	cart_assert( head != NULL );
	cart_assert( ioct >= 0 && ioct < num_octs );

	/* push onto the front O(1) insert, then reordering on iteration */
	oct_next[ioct] = *head;
	oct_prev[ioct] = NULL_OCT;

	if ( *head != NULL_OCT ) {
		oct_prev[*head] = ioct;
	}

	*head = ioct;
}

void linked_list_remove( int *head, int ioct ) {
	cart_assert( head != NULL );
	cart_assert( ioct >= 0 && ioct < num_octs );

	if ( *head == ioct ) {
		*head = oct_next[ioct];
	} else {
		oct_next[oct_prev[ioct]] = oct_next[ioct];
	}
	
	if ( oct_next[ioct] != NULL_OCT ) {	
		oct_prev[oct_next[ioct]] = oct_prev[ioct];
	}

	oct_next[ioct] = NULL_OCT;
	oct_prev[ioct] = NULL_OCT;
}

void linked_list_optimize( int *head, int num_entries ) {
	int count;
	int *order;
	int cur;

	if ( num_entries > 0 ) {
		order = cart_alloc( num_entries * sizeof(int) );
		count = 0;
		cur = *head;

		while ( cur != NULL_OCT ) {
			order[count++] = cur;
			cur = oct_next[cur];
		}

		cart_assert( count == num_entries );

		/* sort the array */
		qsort( order, count, sizeof(int), compare_ints );

		*head = order[0];
		oct_prev[order[0]] = NULL_OCT;
		for ( cur = 1; cur < count; cur++ ) {
			oct_next[order[cur-1]] = order[cur];
			oct_prev[order[cur]] = order[cur-1];
		}
		oct_next[order[count-1]] = NULL_OCT;

		cart_free(order);
	}
}

