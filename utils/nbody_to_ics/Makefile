CC = icc
CFLAGS = -O2 
LIBS = -lm -lgsl -lgslcblas 
#INCLUDES = -L/home/dhr9/gsl-1.12/lib -I/home/dhr9/gsl-1.12/include

nbody_to_ics: nbody_to_ics.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) *.c \
		-o nbody_to_ics \
		$(LIBS)

clean:
	@rm -f nbody_to_ics
