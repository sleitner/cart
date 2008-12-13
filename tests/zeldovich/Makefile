TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
NO_RT = 

include $(TOP_DIR)/Make.config

compile:	zeldovich

run:
	
zeldovich: .obj $(OBJS) defs.h start_zeldovich.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(LIBS) $(OBJS) start_zeldovich.c -o zeldovich

clean:
	rm -f $(OBJS) *.il core start_seldovich.o zeldovich
