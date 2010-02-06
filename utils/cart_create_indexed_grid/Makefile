CC = icc
CFLAGS = -O2
LIBS = 
INCLUDES = 

cart_create_indexed_grid: cart_create_indexed_grid.c sfc.h sfc.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) $(LIBS) sfc.c cart_create_indexed_grid.c -o cart_create_indexed_grid

clean:
	rm -rf *.il core cart_create_indexed_grid
