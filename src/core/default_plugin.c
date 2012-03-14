#include "config.h"


#include <stdio.h>

#include "plugin.h"
#include "timing.h"


void run_output();


void default_plugin_worker()
{
  start_time( OUTPUT_TIMER );
  run_output();
  end_time( OUTPUT_TIMER );
}


plugin_t default_plugin = { NULL };



const plugin_t* add_plugin(int id)
{
  if(id == 0)
    {
      default_plugin.RunBegin = default_plugin.GlobalStepEnd = default_plugin_worker;
      return &default_plugin;
    }
  else return NULL;
}

