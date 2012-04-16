TOP_DIR = $(ART)/CART/trunk
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

.obj/oldstyle_units.o: ../et/oldstyle_units.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) -o .obj/oldstyle_units.o -c ../et/oldstyle_units.c

art: .obj/.phony $(SIM_OBJS) defs.h rt_defs.h start_rt.c .obj/extra/healpix.o .obj/extra/ifrit.o .obj/extra/utils.o .obj/extra/ism.o .obj/extra/halo_finder.o .obj/oldstyle_units.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_rt.c -o art $(LIBS) $(EXTRALIBS) .obj/extra/healpix.o .obj/extra/ifrit.o .obj/extra/utils.o .obj/extra/ism.o .obj/extra/halo_finder.o .obj/oldstyle_units.o

CLEAN := start_rt.o art
