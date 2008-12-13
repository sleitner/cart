TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
PLATFORM = fulla

include $(TOP_DIR)/Make.config

compile: art

run:

art: .obj  $(OBJS) defs.h rt_defs.h start.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start.c -o art $(LIBS) $(EXTRALIBS)

clean:
	rm -rf $(OBJS) .obj .obj/F *.il core start.o art
