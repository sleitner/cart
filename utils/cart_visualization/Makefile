CC = icc
CFLAGS = -O2 
LIBS = -lm -lgsl -lgslcblas
INCLUDES = -I${HOME}/include -L${HOME}/lib

cart_analysis_generic: defs.h cart_analysis_generic.c analysis.c analysis.h
	$(CC) $(CFLAGS) -I. $(INCLUDES) *.c \
		-o cart_analysis_generic \
		$(LIBS)

clean:
	rm -f cart_analysis_generic
