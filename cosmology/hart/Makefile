TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

include $(TOP_DIR)/Make.config

compile:	cart

run:
	
cart: .obj $(OBJS) .obj/extra/viewdump.o defs.h start_cosmology.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_cosmology.c -o cart $(LIBS) $(EXTRALIBS)

clean:
	rm -f $(OBJS) .obj *.il core start_cosmology.o cart
