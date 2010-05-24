CC = icc
CFLAGS = -O2
LIBS = -lm -lgsl -lgslcblas -lsrfftw -lsfftw
INCLUDES = -I${HOME}/include -L${HOME}/lib

cart_power: defs.h cart_power.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) *.c \
		-o cart_power \
		$(LIBS)

clean:
	@rm -f cart_power 
