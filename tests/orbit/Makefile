TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile: orbit

include $(TOP_DIR)/Make.config

orbit: .obj/.phony $(SIM_OBJS) defs.h start_orbit.c .obj/extra/ifrit.o .obj/extra/utils.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_orbit.c -o orbit $(LIBS) $(EXTRALIBS) .obj/extra/ifrit.o .obj/extra/utils.o

CLEAN := start_orbit.o orbit
