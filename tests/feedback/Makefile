TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

.phony: compile run clean

compile:	radp

run: dumps radp
	mpirun -np 1 ./radp radp.cfg

dumps:
	@mkdir dumps

include $(TOP_DIR)/Make.config

radp: .obj $(OBJS) defs.h start_radp.c .obj/extra/viewdump.o .obj/extra/ifrit.o .obj/extra/utils.o slice_output.c slice_output.h
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_radp.c -o radp .obj/extra/viewdump.o .obj/extra/ifrit.o .obj/extra/utils.o slice_output.c $(LIBS) $(EXTRALIBS)

clean:
	@rm -rf .obj radp core
