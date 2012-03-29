TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile:	cart

include $(TOP_DIR)/Make.config

.obj/hart_init.o: $(DEFS) $(TOP_DIR)/cosmology/hart/hart_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/hart_init.o -c $(TOP_DIR)/cosmology/hart/hart_init.c

.obj/hart_reader.o: $(DEFS) $(TOP_DIR)/cosmology/hart/hart_reader.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/hart_reader.o -c $(TOP_DIR)/cosmology/hart/hart_reader.c

cart: .obj/.phony $(DEFS) $(SIM_OBJS) .obj/hart_init.o .obj/hart_reader.o start_cosmology.c 
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) .obj/hart_init.o .obj/hart_reader.o start_cosmology.c -o cart $(LIBS) $(EXTRALIBS)

DAT:
	@mkdir DAT

logs:
	@mkdir logs

ICs:
	@svn checkout https://www.drudd.com/svn/cart/data/cluster6/ICs ICs

prep: DAT logs ICs 

CLEAN := start_cosmology.o cart
