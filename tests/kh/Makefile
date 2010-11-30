TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

compile:	kh	

run: DAT logs
	mpirun -np 1 ./kh kh.cfg

DAT:
	@mkdir DAT

logs:
	@mkdir logs

include $(TOP_DIR)/Make.config

kh: .obj $(OBJS) defs.h start_kh.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_kh.c -o kh $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj *.il core start_kh.o kh list_defines.h
