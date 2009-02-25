TOP_DIR = ../..
SRC_DIR = $(TOP_DIR)/src
NO_RT = 1

include $(TOP_DIR)/Make.config

.phony: compile run clean

compile:	sedov

run:

sedov: .obj $(OBJS) defs.h start_sedov.c
	$(CC) $(CFLAGS) -I. -I$(SRC_DIR) $(INCLUDES) $(OBJS) start_sedov.c -o sedov $(LIBS) $(EXTRALIBS)

clean:
	@rm -f $(OBJS) sedov core
