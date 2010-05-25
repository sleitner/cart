TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
NO_RT = 1

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
	rm -f $(OBJS) *.il core start_kh.o kh
