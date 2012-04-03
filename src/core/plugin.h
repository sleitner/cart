#ifndef __PLUGIN_H__
#define __PLUGIN_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>

/*
//  This file contains prototypes for user plugin functions.
//  These functions are combined in a single structure to
//  simplify extensibility: if a new plugin point is created,
//  it will be set to NULL by default. To insure that old codes
//  remains functional, ALWAYS ADD NEW POINTS TO THE END of this
//  structure. All non-NULL functions pointed to by the structure
//  members need to be defined and linked in at compile time.
//
//  After you add a plugin point here, you need to insert it into
//  the appropriate place in the code using PLUGIN_POINT macro, c.f.
//    PLUGIN_POINT(LevelStepBegin)(level,level_com);
*/
typedef struct Plugin
{
  void (*ConfigInit)();
  void (*ConfigVerify)();
  void (*RunBegin)();
  void (*RunEnd)();
  void (*GlobalStepBegin)();
  void (*GlobalStepEnd)();
  void (*LevelStepBegin)(int level, MPI_Comm level_com);
  void (*LevelStepEnd)(int level, MPI_Comm level_com);
  void (*LevelStepFail)(int level, MPI_Comm level_com);
  /* void (**AfterCFLRestart)(); -- not used? */
  /* void (**StarformationFeedbackEnd)(int level, int cell); -- not used? */
}
plugin_t;


/*
//  This function must be implemented by your plugin code. When being 
//  called subsequently with different id arguments, it should return 
//  const pointers to all user plugins or NULL when no more plugins are
//  needed. It is a responsibility of user's code to keep the actual 
//  plugin_t structures and populate them properly. However, only 
//  functions that are used need to be set in the structure; any 
//  function pointer inside plugin_t can be set to NULL, and then 
//  that plugin point will not be called. 
*/
const plugin_t* add_plugin(int id);


/*
//  The rest of this file is internal implementation
*/
void config_init_plugins();
void config_verify_plugins();

struct PluginList
{
  const plugin_t active;
  const plugin_t* next;
  const plugin_t* const head;
#ifdef __cplusplus
  PluginList();  // Not implemented
#endif
};

extern struct PluginList plugins;

#ifdef DISABLE_PLUGINS
#define PLUGIN_POINT(call)
#else  /* DISABLE_PLUGINS */
#define PLUGIN_POINT(call) \
  if(plugins.active.call != NULL) for(plugins.next=plugins.head; plugins.next->call!=NULL; plugins.next++) plugins.next->call
#endif /* DISABLE_PLUGINS */

#endif

