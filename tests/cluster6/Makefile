TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
NO_RT = 1

compile:	cart

include $(TOP_DIR)/Make.config

cart: .obj $(OBJS) defs.h start_cosmology.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_cosmology.c -o cart $(LIBS) $(EXTRALIBS)

DAT:
	@mkdir DAT

logs:
	@mkdir logs

ICs:
	@svn checkout https://www.drudd.com/svn/cart/data/cluster6 ICs

prep: DAT logs ICs 
	
clean:
	rm -rf $(OBJS) .obj *.il core start_cosmology.o cart 
