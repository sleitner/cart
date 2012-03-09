#include "config.h"


#include "plugin.h"
#include "timing.h"


void default_plugin_worker()
{
  start_time( OUTPUT_TIMER );
  run_output();
  end_time( OUTPUT_TIMER );
}


void set_plugin(struct Plugin *plugin)
{
  plugin->RunBegin = default_plugin_worker;
  plugin->GlobalStepEnd = default_plugin_worker;
}

