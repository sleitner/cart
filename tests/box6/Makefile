TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

art: .obj defs.h rt_defs.h $(OBJS) start.c 
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start.c -o art $(LIBS) $(EXTRALIBS) 
	which mpirun

analysis: .obj $(OBJS) $(EXT_OBJS) defs.h rt_defs.h analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) $(EXT_OBJS) analysis.c -o analysis $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj *.il core start.o art analysis
