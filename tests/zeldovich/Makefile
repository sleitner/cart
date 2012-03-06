TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile:	zeldovich

include $(TOP_DIR)/Make.config

zeldovich: .obj/.phony $(SIM_OBJS) defs.h start_zeldovich.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_zeldovich.c -o zeldovich $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj core start_seldovich.o zeldovich
