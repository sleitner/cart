TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

gic_reader.o: ../gic/gic_reader.c
	$(CC) $(CFLAGS) -I. -I../gic -I$(SRC_DIR) $(INCLUDES) -c ../gic/gic_reader.c

cosmics_init.o: $(DEFS) cosmics_init.c
	$(CC) $(CFLAGS) -I. -I../gic -I$(SRC_DIR) $(INCLUDES) -c cosmics_init.c

art: .obj $(DEFS) $(OBJS) start_cosmology.c gic_reader.o cosmics_init.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) cosmics_init.o gic_reader.o start_cosmology.c -o art $(LIBS) $(EXTRALIBS)
	which mpirun

analysis: .obj $(DEFS) $(OBJS) $(EXT_OBJS) analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) $(EXT_OBJS) analysis.c -o analysis $(LIBS) $(EXTRALIBS)
	which mpirun

clean:
	rm -rf .obj *.il core start_cosmology.o gic_reader.o cosmics_init.o art analysis
