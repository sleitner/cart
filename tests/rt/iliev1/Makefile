TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

art: .obj .obj/rt .obj/rt/F $(OBJS) defs.h rt_defs.h start_rt.c ../common/pix2ang_nest.o ../common/mk_pix2xy.o .obj/extra/ifrit.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_rt.c -o art $(LIBS) $(EXTRALIBS) ../common/pix2ang_nest.o ../common/mk_pix2xy.o .obj/extra/ifrit.o

clean:
	rm -rf .obj *.o art core *.il
