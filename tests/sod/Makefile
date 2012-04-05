TOP_DIR = ../..

compile:	sod

include $(TOP_DIR)/Make.config

sod: .obj/.phony $(SIM_OBJS) defs.h start_sod.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_sod.c -o sod $(LIBS) $(EXTRALIBS)

CLEAN := start_sod.o sod
