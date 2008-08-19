TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src
PLATFORM = idunn

include $(TOP_DIR)/Make.config

compile: art

run:

art: .obj .obj/F $(OBJS) defs.h start_rt.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_rt.c -o art $(LIBS) $(EXTRALIBS)

clean:
	rm -rf $(OBJS) .obj .obj/F *.il core start_rc.o art
