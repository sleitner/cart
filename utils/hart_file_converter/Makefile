TOP_DIR = ../../
SRC_DIR = $(TOP_DIR)/src
NO_RT = 1

all: hart_to_cart cart_to_hart

include $(TOP_DIR)/Make.config

NEWSRCS = $(filter-out main.c,$(SRCS)) 
NEWOBJS = $(patsubst %.c,.obj/./%.o,$(NEWSRCS)) .obj/extra/hart_io.o

hart_to_cart: .obj $(NEWOBJS) defs.h hart_to_cart.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(NEWOBJS) \
		hart_to_cart.c -o hart_to_cart \
		$(LIBS) $(EXTRALIBS)

cart_to_hart: .obj $(NEWOBJS) defs.h cart_to_hart.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(NEWOBJS) \
                cart_to_hart.c -o cart_to_hart \
                $(LIBS) $(EXTRALIBS)

clean:
	@rm -rf .obj *.o *.il core hart_to_cart cart_to_hart
