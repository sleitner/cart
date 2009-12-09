TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
NO_RT = 1

compile:	art

include $(TOP_DIR)/Make.config

hart_init.o: $(DEFS) hart_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) -c hart_init.c

hart_reader.o: $(DEFS) hart_reader.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) -c hart_reader.c

art: .obj $(DEFS) $(OBJS) hart_init.o hart_reader.o start_cosmology.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) hart_init.o hart_reader.o start_cosmology.c -o cart $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj *.il core* start_cosmology.o hart_init.o hart_reader.o art
