#include "defs.h"
#include <stdio.h>
#define STR_VALUE(arg)      #arg
#define to_string(name)     STR_VALUE(name)
int main()
{
#ifdef SF_FEEDBACK
  printf("%s\n",to_string(SF_FEEDBACK));
#endif
  return 0;
}

