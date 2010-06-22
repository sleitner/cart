CC = gcc
CFLAGS = -O2
LIBS = -lm -lgsl -lgslcblas
INCLUDES = -L/sw/lib -I/sw/include

halo_analysis: defs.h halo_analysis.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) *.c \
		-o halo_analysis \
		$(LIBS)

clean:
	rm -f halo_analysis
