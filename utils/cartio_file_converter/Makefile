TOP_DIR = ../../
SRC_DIR = $(TOP_DIR)/src

all: cart_to_cartio 

include $(TOP_DIR)/Make.config

cart_to_cartio: .obj $(OBJS2) defs.h rt_defs.h cart_to_cartio.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS2) \
                cart_to_cartio.c -o cart_to_cartio \
                $(LIBS) $(EXTRALIBS)

#cartio_to_cart: .obj $(OBJS2) defs.h rt_defs.h cartio_to_cart.c
#	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS2) \
#			cartio_to_cart.c -o cartio_to_cart \
#			$(LIBS) $(EXTRALIBS)

clean:
	@rm -rf .obj *.o *.il core cart_to_cartio
