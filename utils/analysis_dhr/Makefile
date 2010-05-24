CC = icc
CFLAGS = -O2
LIBS = -lm -lgsl -lgslcblas
INCLUDES = 

cart_galaxy_finder: defs.h cart_analysis_serial.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) *.c \
		-o cart_analysis_serial \
		$(LIBS)

clean:
	rm -f cart_analysis_serial
