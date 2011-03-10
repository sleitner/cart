TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src

.phony: compile run clean

compile:	sedov

run: dumps sedov
	mpirun -np 1 ./sedov sedov.cfg

dumps:
	@mkdir dumps

include $(TOP_DIR)/Make.config

sedov: .obj $(OBJS) defs.h start_sedov.c .obj/extra/viewdump.o .obj/extra/ifrit.o .obj/extra/utils.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_sedov.c -o sedov .obj/extra/viewdump.o .obj/extra/ifrit.o .obj/extra/utils.o $(LIBS) $(EXTRALIBS)

clean:
	@rm -rf .obj sedov core
