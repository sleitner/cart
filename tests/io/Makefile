TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile: art

run:

include $(TOP_DIR)/Make.config

art: .obj/.phony $(SIM_OBJS) defs.h start_io.c 
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_io.c -o art $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj *.il core start_io.o art
