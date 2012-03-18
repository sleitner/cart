#ifndef __STACK_H__
#define __STACK_H__

#define STACK_PAGE_SIZE  1024

typedef struct STACK_PAGE {
	struct STACK_PAGE* next;
	int index;
	int data[STACK_PAGE_SIZE];
} stack_page;

typedef struct STACK {
    int size;
	stack_page *head;
} stack;

stack *stack_init();
int stack_size( stack * );
void stack_print( stack * );
void stack_destroy( stack * );
void stack_push( stack *, int );
int stack_pop( stack *, int * );

#endif /* __STACK_H__ */
