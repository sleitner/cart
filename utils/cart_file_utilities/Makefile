TOP_DIR = ../../
ifdef PLATFORM                                                                                                                                                                
include $(TOP_DIR)/platforms/$(PLATFORM)
else
@echo "A variable PLATFORM must be set. Available platforms:"
@ls -m $(TOP_DIR)/platforms
@cat $(TOP_DIR)/platforms >& /dev/null
endif

all: print_particle_info print_gas_info

print_particle_info: print_particle_info.c particle_io.c particle_io.h
	$(CC) $(CFLAGS) -I. $(INCLUDES) $(LIBS) $(EXTRALIBS) print_particle_info.c particle_io.c -o print_particle_info

print_gas_info: print_gas_info.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) $(LIBS) $(EXTRALIBS) print_gas_info.c -o print_gas_info

convert_particle_endianness: convert_particle_endianness.c
	$(CC) $(CFLAGS) -I. $(INCLUDES) $(LIBS) $(EXTRALIBS) convert_particle_endianness.c -o convert_particle_endianness

clean:
	rm -rf *.il core print_particle_info print_tracer_info print_gas_info convert_particle_endianness
