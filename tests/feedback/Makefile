TOP_DIR = ../../
SRC_DIR = $(TOP_DIR)/src

.phony: compile run clean

compile:	radp clearfile

run: dumps radp
	mpirun -np 1 ./radp radp.cfg

dumps:
	@mkdir dumps

include $(TOP_DIR)/Make.config

#all: .obj/start_analysis.o radp

VIEWOBJ = .obj/extra/viewdump.o .obj/extra/ifrit.o .obj/extra/utils.o .obj/extra/output_slice.o 

start_analysis.o: start_radp.h start_analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(VIEWOBJ) -c start_analysis.c 

radp: .obj $(OBJS) defs.h start_radp.h start_radp.c start_analysis.c $(VIEWOBJ) start_analysis.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) $(VIEWOBJ) start_analysis.o start_radp.c -o radp  $(LIBS) $(EXTRALIBS)

clean:
	@rm -rf .obj radp core *.o *.log

clearfile: 
	@touch totmompress.dat
	@cp -f totmompress.dat totmompress.dat.old
	@rm -f totmompress.dat
