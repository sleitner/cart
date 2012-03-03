TOP_DIR = ../../
SRC_DIR = $(TOP_DIR)/src

compile: test-base test-analysis test-file-converter test-run-gic test-run-hart test-run-cosmics

include $(TOP_DIR)/Make.config

test-base: .obj/.phony $(DEFS) $(BAS_OBJS) test-base.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) $(BAS_OBJS) test-base.c -o test-base $(LIBS) $(EXTRALIBS)

test-analysis: .obj/.phony $(DEFS) $(ANL_OBJS) test-analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) $(ANL_OBJS) test-analysis.c -o test-analysis $(LIBS) $(EXTRALIBS)

test-file-converter: .obj/.phony $(DEFS) $(FIC_OBJS) .obj/hart_io.o test-file-converter.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) $(FIC_OBJS) .obj/hart_io.o test-file-converter.c -o test-file-converter $(LIBS) $(EXTRALIBS)

test-run-gic: .obj/.phony $(DEFS) $(SIM_OBJS) .obj/gic_reader.o .obj/gic_init.o  test-run-gic.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(TOP_DIR)/cosmology/gic $(INCLUDES) $(SIM_OBJS) .obj/gic_reader.o .obj/gic_init.o test-run-gic.c -o test-run-gic $(LIBS) $(EXTRALIBS)

.obj/gic_reader.o: $(TOP_DIR)/cosmology/gic/gic_reader.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) -o .obj/gic_reader.o -c $(TOP_DIR)/cosmology/gic/gic_reader.c

.obj/gic_init.o: $(DEFS) $(TOP_DIR)/cosmology/gic/gic_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/gic_init.o -c $(TOP_DIR)/cosmology/gic/gic_init.c

test-run-hart: .obj/.phony $(DEFS) $(SIM_OBJS) .obj/hart_reader.o .obj/hart_init.o  test-run-gic.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(TOP_DIR)/cosmology/hart $(INCLUDES) $(SIM_OBJS) .obj/hart_reader.o .obj/hart_init.o test-run-hart.c -o test-run-hart $(LIBS) $(EXTRALIBS)

.obj/hart_init.o: $(DEFS) $(TOP_DIR)/cosmology/hart/hart_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/hart_init.o -c $(TOP_DIR)/cosmology/hart/hart_init.c

.obj/hart_reader.o: $(DEFS) $(TOP_DIR)/cosmology/hart/hart_reader.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/hart_reader.o -c $(TOP_DIR)/cosmology/hart/hart_reader.c

test-run-cosmics: .obj/.phony $(DEFS) $(SIM_OBJS) .obj/gic_reader.o .obj/cosmics_init.o  test-run-gic.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core -I$(TOP_DIR)/cosmology/cosmics $(INCLUDES) $(SIM_OBJS) .obj/gic_reader.o .obj/cosmics_init.o test-run-cosmics.c -o test-run-cosmics $(LIBS) $(EXTRALIBS)

.obj/cosmics_init.o: $(DEFS) $(TOP_DIR)/cosmology/cosmics/cosmics_init.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base -I$(SRC_DIR)/core $(INCLUDES) -o .obj/cosmics_init.o -c $(TOP_DIR)/cosmology/cosmics/cosmics_init.c

clean:
	rm -rf .obj core *.o test-base test-analysis test-file-converter test-run-gic test-run-hart
