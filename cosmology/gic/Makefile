TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
PLATFORM = fulla

include $(TOP_DIR)/Make.config

compile: art

run:

art: .obj  $(OBJS) defs.h rt_defs.h start_cosmology_gic.c gic_reader.o .obj/extra/ifrit.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) .obj/extra/ifrit.o gic_reader.o start_cosmology_gic.c -o art $(LIBS) $(EXTRALIBS)

analysis: .obj $(OBJS) $(EXT_OBJS) defs.h rt_defs.h analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) $(EXT_OBJS) analysis.c -o analysis $(LIBS) $(EXTRALIBS)

clean:
	rm -rf $(OBJS) .obj .obj/F *.il core start_cosmology_gic.o art analysis
