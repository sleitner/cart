TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src
PLATFORM = idunn

include $(TOP_DIR)/Make.config

compile: art

run:

art: .obj .obj/F $(OBJS) defs.h start_rt.c pix2ang_nest.o mk_pix2xy.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_rt.c -o art $(LIBS) $(EXTRALIBS) pix2ang_nest.o mk_pix2xy.o

clean:
	rm -rf $(OBJS) .obj .obj/F start_rt.o art core *.il
