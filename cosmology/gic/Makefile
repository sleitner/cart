TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

gic_reader.o: gic_reader.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) -c gic_reader.c

gic_init.o: $(DEFS) gic_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) -c gic_init.c

art: .obj $(DEFS) $(OBJS) start_cosmology.c gic_init.o gic_reader.o .obj/extra/ifrit.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) .obj/extra/ifrit.o gic_reader.o gic_init.o start_cosmology.c -o art $(LIBS) $(EXTRALIBS)
	which mpirun

analysis: .obj $(DEFS) $(OBJS) $(EXT_OBJS) analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) $(EXT_OBJS) analysis.c -o analysis $(LIBS) $(EXTRALIBS)
	which mpirun

clean:
	rm -rf .obj *.il core start_cosmology.o gic_init.o gic_reader.o art analysis
