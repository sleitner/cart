TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
NO_RT = 1

compile:	sod

run:

include $(TOP_DIR)/Make.config

sod: .obj $(OBJS) defs.h start_sod.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_sod.c -o sod $(LIBS) $(EXTRALIBS)

clean:
	rm -f $(OBJS) *.il core start_sod.o sod
