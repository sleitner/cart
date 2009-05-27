TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src
PLATFORM = fulla

compile: art

include $(TOP_DIR)/Make.config

run:

art: .obj  $(OBJS) defs.h start_rt.c .obj/extra/ifrit.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) .obj/extra/ifrit.o start_rt.c -o art $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj *.il core start_rc.o art
