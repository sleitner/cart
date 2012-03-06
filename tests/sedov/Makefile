TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

.phony: compile run clean

compile:	sedov

run: dumps sedov
	mpirun -np 1 ./sedov sedov.cfg

dumps:
	@mkdir dumps

include $(TOP_DIR)/Make.config

sedov: .obj/.phony $(SIM_OBJS) defs.h start_sedov.c .obj/viewdump.o .obj/ifrit.o .obj/utils.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) start_sedov.c -o sedov .obj/viewdump.o .obj/ifrit.o .obj/utils.o $(LIBS) $(EXTRALIBS)

clean:
	@rm -rf .obj sedov core start_sedov.o
