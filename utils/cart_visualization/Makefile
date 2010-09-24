CC = icc
CFLAGS = -O2 
LIBS = -lm -lgsl -lgslcblas
INCLUDES = -I${HOME}/include -L${HOME}/lib

cart_visualization: defs.h cart_visualization.c analysis.c analysis.h
	$(CC) $(CFLAGS) -I. $(INCLUDES) *.c \
		-o cart_visualization \
		$(LIBS)

clean:
	rm -f cart_visualization
