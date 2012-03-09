#include "config.h"


#include "auxiliary.h"
#include "io_cart.h"
#include "../extra/hart_io.h"


extern const char* executable_name;


void init()
{
  if(num_options != 2)
    {
      cart_error("Usage: %s <input> <output>",executable_name);
    }
}


void read_file(const char* fname)
{
  read_hart_grid_binary((char *)fname);
}


void write_file(const char* fname)
{
  write_cart_grid_binary((char *)fname);
}
