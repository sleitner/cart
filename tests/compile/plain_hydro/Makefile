TOP_DIR = ../../..
WITH_EXTRA = 1

compile: art

include $(TOP_DIR)/Make.config

art: .obj/.phony $(DEFS) $(SIM_OBJS) start.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start.c -o art $(LIBS) $(EXTRALIBS)

ALL += art
CLEAN += art start.o
