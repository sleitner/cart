TOP_DIR = ../../
SRC_DIR = $(TOP_DIR)/src

include $(TOP_DIR)/Make.config

compile:	cart

run:
	
cart: .obj $(OBJS) defs.h start_cosmology.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_cosmology.c -o cart $(LIBS) $(EXTRALIBS)

clean:
	rm -rf $(OBJS) .obj *.il core start_cosmology.o cart logs/* output/*
