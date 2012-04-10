#include "defs.h"
#include <stdio.h>
#define STR_VALUE(arg)      #arg
#define to_string(name)     STR_VALUE(name)
int main()
{
#ifdef SF_RECIPE
  printf("%s\n",to_string(SF_RECIPE));
#endif
  return 0;
}

