TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile:	zeldovich

include $(TOP_DIR)/Make.config

zeldovich: .obj $(OBJS) defs.h start_zeldovich.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_zeldovich.c -o zeldovich $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj *.il core start_seldovich.o zeldovich
