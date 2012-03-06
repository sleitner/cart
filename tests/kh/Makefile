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

kh: .obj/.phony $(SIM_OBJS) defs.h start_kh.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_kh.c -o kh $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj core start_kh.o kh
