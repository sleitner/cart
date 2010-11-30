TOP_DIR = ../../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

art: .obj $(OBJS) defs.h start_rt.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_rt.c -o art $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj core start_rc.o art list_defines.h
