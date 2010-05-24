CC = gcc
CFLAGS = -O2 -m64 
LIBS = -lm -lgsl -lgslcblas
INCLUDES = -I${HOME}/include -L${HOME}/lib

all: defs.h cart_halo_phase_diagram.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) *.c \
		-o cart_halo_phase_diagram \
		$(LIBS)

clean: 
	@rm -f cart_halo_phase_diagram 
