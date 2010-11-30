TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

art: .obj $(OBJS) $(DEFS) start_rt.c .obj/extra/healpix.o .obj/extra/ifrit.o .obj/extra/utils.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_rt.c -o art $(LIBS) $(EXTRALIBS) .obj/extra/healpix.o .obj/extra/ifrit.o .obj/extra/utils.o

clean:
	rm -rf .obj start_rt.o art core *.il list_defines.h
