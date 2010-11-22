#include "config.h"


#include <stdlib.h>

#include "plugin.h"


#ifdef USER_PLUGIN

#define MAX_COUNT          100
const plugin_t* plugins_internal[MAX_COUNT+1];

struct PLUGIN_LIST_TYPE plugins = { 0, plugins_internal };


void config_plugins()
{
  int i, j;
  const plugin_t *ptr;
  
  for(i=0; i<=MAX_COUNT; i++) plugins_internal[i] = NULL;

  i = 0;
  do
    {
      ptr = add_plugin(i);
      if(ptr == NULL) break;

      for(j=0; j<i; j++) 
	{
	  if(plugins_internal[j] == ptr) break;
	}
      if(j < i) break;

      plugins_internal[i++] = ptr;
    }
  while(i < MAX_COUNT);
}

#endif
