TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

art: .obj $(OBJS) defs.h start_orbit.c .obj/extra/ifrit.o .obj/extra/utils.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_orbit.c -o art $(LIBS) $(EXTRALIBS) .obj/extra/ifrit.o .obj/extra/utils.o

clean:
	rm -rf .obj core start_orbit.o art list_defines.h
