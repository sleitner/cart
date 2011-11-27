TOP_DIR = ../../
SRC_DIR = $(TOP_DIR)/src

all: check_placement

include $(TOP_DIR)/Make.config

check_placement: .obj .obj/./system.o check_placement.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) .obj/system.o check_placement.c -o check_placement

clean:
	@rm -rf .obj *.o *.il core check_placement
