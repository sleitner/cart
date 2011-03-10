TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile:	sod

run:

include $(TOP_DIR)/Make.config

sod: .obj $(OBJS) defs.h start_sod.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_sod.c -o sod $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj *.il core start_sod.o sod
