TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

.obj/oldstyle_units.o: ../et/oldstyle_units.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) -o .obj/oldstyle_units.o -c ../et/oldstyle_units.c

art: .obj/.phony $(SIM_OBJS) defs.h start_rt.c .obj/oldstyle_units.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_rt.c .obj/oldstyle_units.o -o art $(LIBS) $(EXTRALIBS)

CLEAN := start_rc.o art
