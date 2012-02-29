TOP_DIR = ../../
SRC_DIR = $(TOP_DIR)/src

compile: test-base test-analysis

include $(TOP_DIR)/Make.config

test-base: .obj/base/.phony $(DEFS) $(BAS_OBJS) test-base.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) $(BAS_OBJS) test-base.c -o test-base $(LIBS) $(EXTRALIBS)

test-analysis: .obj/base/.phony $(DEFS) $(ANL_OBJS) test-analysis.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR)/base $(INCLUDES) $(ANL_OBJS) test-analysis.c -o test-analysis $(LIBS) $(EXTRALIBS)

clean:
	rm -rf .obj core *.o test-base test-analysis
