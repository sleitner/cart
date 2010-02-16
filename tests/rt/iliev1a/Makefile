TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

art: .obj $(OBJS) defs.h rt_defs.h start_rt.c .obj/extra/healpix.o .obj/extra/ifrit.o .obj/extra/utils.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_rt.c -o art $(LIBS) $(EXTRALIBS) .obj/extra/healpix.o .obj/extra/ifrit.o .obj/extra/utils.o

clean:
	rm -rf .obj *.o art core *.il
