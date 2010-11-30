TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile: art

run:

include $(TOP_DIR)/Make.config

art: .obj $(OBJS) defs.h start_io.c 
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_io.c -o art $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj *.il core start_io.o art list_defines.h
