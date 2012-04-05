TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile: art

include $(TOP_DIR)/Make.config

.obj/hart_init.o: $(DEFS) $(TOP_DIR)/cosmology/hart/hart_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/hart_init.o -c $(TOP_DIR)/cosmology/hart/hart_init.c

.obj/hart_reader.o: $(DEFS) $(TOP_DIR)/cosmology/hart/hart_reader.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/hart_reader.o -c $(TOP_DIR)/cosmology/hart/hart_reader.c

art: .obj/.phony $(DEFS) $(SIM_OBJS) start.c .obj/hart_init.o .obj/hart_reader.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start.c .obj/hart_init.o .obj/hart_reader.o -o art $(LIBS) $(EXTRALIBS) 

CLEAN := start.o art
