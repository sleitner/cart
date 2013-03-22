#include "config.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "auxiliary.h"
#include "plugin.h"


#define MAX_COUNT          100


plugin_t plugins_internal[MAX_COUNT+1]; /* the last one is all NULLs */
struct PluginList plugins; // = { { NULL }, NULL, (const plugin_t*)&plugins_internal };


typedef void (*member_t)();


void config_init_plugins()
{
  const int num_members = sizeof(plugin_t)/sizeof(member_t);
  int j, k, m, num_plugins;
  const plugin_t *ptr;
  member_t *w1, *w2;
  int max_member[num_members];
  
  for(m=0; m<num_members; m++) max_member[m] = 0;

  for(k=0; k<MAX_COUNT; k++)
    {
      ptr = add_plugin(k);
      if(ptr == NULL) break;

      /*
      //  Store members indepedently, since plugins are expected to
      //  be sparse; that will economize on looping over plugins - 
      //  we loop over non-NULL members instead.
      */
      w1 = (member_t*)ptr;
      for(m=0; m<num_members; m++) if(w1[m] != NULL)
	{
	  for(j=0; j<max_member[m]; j++) 
	    {
	      w2 = (member_t*)&plugins_internal[j];
	      if(w2[m] == w1[m]) /* already exists */
		{
		  break;
		}
	    }
	  if(j == max_member[m]) /* a new member */
	    {
	      w2 = (member_t*)&plugins_internal[max_member[m]++];
	      w2[m] = w1[m];
	    }
	}
    }

  num_plugins = 0;
  for(m=0; m<num_members; m++) if(max_member[m] > num_plugins) num_plugins = max_member[m];

  /*
  //  Assign the list members. We force the assignment, since
  //  we declared .active as const to avoid user messing with
  //  it.
  */
  *((plugin_t*)&plugins.active) = plugins_internal[0];
  *((const plugin_t**)&plugins.head) = (const plugin_t*)&plugins_internal;

#ifdef DEBUG
  cart_debug("Plugin map:");

  for(j=0; j<num_plugins; j++)
    {
      cart_debug("i=%d ptr=%p",j,plugins_internal[j]);
      w1 = (member_t*)&plugins_internal[j];
      for(m=0; m<num_members; m++)
	{
	  cart_debug("   m=%d ptr[m]=%p",m,w1[m]);
	}

    }
#endif

  PLUGIN_POINT(ConfigInit)();
}

void config_verify_plugins()
{
  PLUGIN_POINT(ConfigVerify)();
}

#ifdef __cplusplus
PluginList::PluginList() : active(plugin_t()), next(0), head(0){}
#endif
