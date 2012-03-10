TOP_DIR = ../
SRC_DIR = $(TOP_DIR)/src

compile: demo-base demo-analysis demo-file-converter demo-run-gic demo-run-hart demo-run-cosmics

include $(TOP_DIR)/Make.config

demo-base: .obj/.phony $(DEFS) $(BAS_OBJS) demo-base.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) $(BAS_OBJS) demo-base.c -o demo-base $(LIBS) $(EXTRALIBS)

demo-analysis: .obj/.phony $(DEFS) $(ANL_OBJS) demo-analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) $(ANL_OBJS) demo-analysis.c -o demo-analysis $(LIBS) $(EXTRALIBS)

demo-file-converter: .obj/.phony $(DEFS) $(FIC_OBJS) .obj/extra/hart_io.o demo-file-converter.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(FIC_OBJS) .obj/extra/hart_io.o demo-file-converter.c -o demo-file-converter $(LIBS) $(EXTRALIBS)

demo-run-hart: .obj/.phony $(DEFS) $(SIM_OBJS) .obj/hart_reader.o .obj/hart_init.o  demo-run-gic.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(TOP_DIR)/cosmology/hart $(INCLUDES) $(SIM_OBJS) .obj/hart_reader.o .obj/hart_init.o demo-run-hart.c -o demo-run-hart $(LIBS) $(EXTRALIBS)

.obj/hart_init.o: $(DEFS) $(TOP_DIR)/cosmology/hart/hart_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/hart_init.o -c $(TOP_DIR)/cosmology/hart/hart_init.c

.obj/hart_reader.o: $(DEFS) $(TOP_DIR)/cosmology/hart/hart_reader.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/hart_reader.o -c $(TOP_DIR)/cosmology/hart/hart_reader.c

demo-run-cosmics: .obj/.phony $(DEFS) $(SIM_OBJS) .obj/gic_reader.o .obj/cosmics_init.o  demo-run-gic.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(TOP_DIR)/cosmology/cosmics $(INCLUDES) $(SIM_OBJS) .obj/gic_reader.o .obj/cosmics_init.o demo-run-cosmics.c -o demo-run-cosmics $(LIBS) $(EXTRALIBS)

.obj/cosmics_init.o: $(DEFS) $(TOP_DIR)/cosmology/cosmics/cosmics_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/cosmics_init.o -c $(TOP_DIR)/cosmology/cosmics/cosmics_init.c

EXEC = demo-run-gic
include $(TOP_DIR)/cosmology/gic/Make.include


CLEAN += *.o demo-base demo-analysis demo-file-converter demo-run-gic demo-run-hart demo-run-cosmics
