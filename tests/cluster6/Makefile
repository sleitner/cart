TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile:	cart

include $(TOP_DIR)/Make.config

hart_reader.o: $(DEFS) hart_reader.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) -c hart_reader.c

cart: .obj $(OBJS) hart_reader.o defs.h start_cosmology.c 
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) hart_reader.o start_cosmology.c -o cart $(LIBS) $(EXTRALIBS)

DAT:
	@mkdir DAT

logs:
	@mkdir logs

ICs:
	@svn checkout https://www.drudd.com/svn/cart/data/cluster6/ICs ICs

prep: DAT logs ICs 
	
clean:
	rm -rf $(OBJS) .obj *.il core start_cosmology.o cart
