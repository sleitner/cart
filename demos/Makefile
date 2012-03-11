TOP_DIR = ../
SRC_DIR = $(TOP_DIR)/src

compile: demo-base demo-analysis demo-file-converter demo-run

include $(TOP_DIR)/Make.config

demo-base: .obj/.phony $(DEFS) $(BAS_OBJS) demo-base.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) $(BAS_OBJS) demo-base.c -o demo-base $(LIBS) $(EXTRALIBS)

demo-analysis: .obj/.phony $(DEFS) $(ANL_OBJS) demo-analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) $(ANL_OBJS) demo-analysis.c -o demo-analysis $(LIBS) $(EXTRALIBS)

demo-file-converter: .obj/.phony $(DEFS) $(FIC_OBJS) .obj/extra/hart_io.o demo-file-converter.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(FIC_OBJS) .obj/extra/hart_io.o demo-file-converter.c -o demo-file-converter $(LIBS) $(EXTRALIBS)

demo-run: .obj/.phony $(DEFS) $(SIM_OBJS) demo-run.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(SIM_OBJS) demo-run.c -o demo-run $(LIBS) $(EXTRALIBS)

CLEAN += *.o demo-base demo-analysis demo-file-converter demo-run
