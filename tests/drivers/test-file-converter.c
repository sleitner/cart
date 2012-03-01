#include "config.h"


#include "io_art.h"
#include "../extra/hart_io.h"


void read_file(const char* fname)
{
  read_hart_grid_binary((char *)fname);
}


void write_file(const char* fname)
{
  write_art_grid_binary((char *)fname);
}
