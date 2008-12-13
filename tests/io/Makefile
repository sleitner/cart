TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
NO_RT = 
#PLATFORM = idunn.debug

include $(TOP_DIR)/Make.config

compile: art

run:

art: .obj $(OBJS) defs.h start_io.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_io.c -o art $(LIBS) $(EXTRALIBS)

clean:
	rm -rf $(OBJS) .obj *.il core start_io.o art
