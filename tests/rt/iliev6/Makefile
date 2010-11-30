TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

art: .obj  $(OBJS) defs.h rt_defs.h start_rt.c .obj/extra/ifrit.o .obj/extra/utils.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) .obj/extra/ifrit.o start_rt.c -o art $(LIBS) $(EXTRALIBS) .obj/extra/utils.o

clean:
	rm -rf .obj *.il core start_rc.o art list_defines.h
