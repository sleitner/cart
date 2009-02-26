TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
NO_RT = 1

include $(TOP_DIR)/Make.config

compile:	sod

run:
	
sod: .obj $(OBJS) defs.h start_sod.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_sod.c -o sod $(LIBS) $(EXTRALIBS)

clean:
	rm -f $(OBJS) *.il core start_sod.o sod
