#include "config.h"


#include <stdio.h>

#include "plugin.h"


struct Plugin plugin_internal = { NULL };
const struct Plugin *plugin = &plugin_internal;

void config_plugin()
{
  set_plugin(&plugin_internal);
}
