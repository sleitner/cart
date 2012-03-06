TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

.obj/oldstyle_units.o: ../et/oldstyle_units.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) -o .obj/oldstyle_units.o -c ../et/oldstyle_units.c

art: .obj/.phony $(SIM_OBJS) $(DEFS) start_rt.c .obj/healpix.o .obj/ifrit.o .obj/utils.o .obj/oldstyle_units.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_rt.c -o art $(LIBS) $(EXTRALIBS) .obj/healpix.o .obj/ifrit.o .obj/utils.o .obj/oldstyle_units.o

clean:
	rm -rf .obj start_rt.o art core
