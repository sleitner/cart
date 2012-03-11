TOP_DIR = ../../
SRC_DIR = $(TOP_DIR)/src

.phony: compile run clean

compile:	radp clearfile

run: dumps radp
	mpirun -np 1 ./radp radp.cfg

dumps:
	@mkdir dumps

PLUGIN = /no_default

include $(TOP_DIR)/Make.config

#all: .obj/start_analysis.o radp

VIEWOBJ = .obj/extra/viewdump.o .obj/extra/ifrit.o .obj/extra/utils.o .obj/extra/output_slice.o 

start_analysis.o: start_radp.h start_analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(VIEWOBJ) -c start_analysis.c 

radp: .obj/.phony $(SIM_OBJS) defs.h start_radp.h start_radp.c start_analysis.c $(VIEWOBJ) start_analysis.o
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) $(VIEWOBJ) start_analysis.o start_radp.c -o radp  $(LIBS) $(EXTRALIBS)

CLEAN += radp *.o *.log

clearfile: 
	@touch totmompress.dat
	@rm -f totmompress.dat.old
	@cp totmompress.dat totmompress.dat.old
	@rm -f totmompress.dat
