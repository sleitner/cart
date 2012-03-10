#include "config.h"


#include "auxiliary.h"


extern const char* executable_name;

const char *input;
const char *output;

void init()
{
  const char *str;

  cart_debug("First, we check some command-line options.");

  /*
  // For example, like that:
  */
  input = extract_option1("input","i",NULL);
  if(input == NULL)
    {
      cart_error("Command-line option -i/--input=<name> must be specified.");
    }

  output = extract_option1("output","o",NULL);
  if(output == NULL)
    {
      cart_error("Command-line option -o/--output=<name> must not specified.");
    }
}


void read_file()
{
  cart_debug("Second, we read the input data file %s.",input);
}


void write_file(const char* fname)
{
  cart_debug("Last, we write the output data file %s.",output);
}
