CC = icc
CFLAGS = -O2
LIBS = 
INCLUDES = 

cart_create_indexed_particles: cart_create_indexed_particles.c sfc.h sfc.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) $(LIBS) sfc.c cart_create_indexed_particles.c -o cart_create_indexed_particles

clean:
	rm -rf *.il core cart_create_indexed_particles
