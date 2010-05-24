CC = icc
CFLAGS = -O2
LIBS = -lm -lgsl -lgslcblas
INCLUDES =

cart_analysis_sz: defs.h cart_analysis_sz.c analysis.h analysis.c io.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) *.c \
		-o cart_analysis_sz \
		$(LIBS)

clean:
	@rm -f cart_analysis_sz
