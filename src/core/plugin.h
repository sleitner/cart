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
//  it will be set to NULL by default. That way old codes
//  do NOT need to be modified. All non-NULL functions need to be 
//  defined and linked in at compile time.
//  The order in which plugin points appear in this structure is
//  irrelevant, so it is recommended that you add your points at the end.
*/
struct Plugin
{
  void (*RunBegin)();
  void (*RunEnd)();
  void (*GlobalStepBegin)();
  void (*GlobalStepEnd)();
  void (*LevelStepBegin)(int level, MPI_Comm level_com);
  void (*LevelStepEnd)(int level, MPI_Comm level_com);
  void (*LevelStepFail)(int level, MPI_Comm level_com);
  /* void (*AfterCFLRestart)(); -- not used? */
  /* void (*StarformationFeedbackEnd)(int level, int cell); -- not used? */
}
;

extern const struct Plugin *plugin;

/*
//  This function must be implemented by your plugin code. The pointer
//  plugin comes non-NULL, but all its members are set to NULL. Reset
//  those of them that you actually need to your own defined funtions.
*/
void set_plugin(struct Plugin *plugin);


#define PLUGIN_POINT(call) \
  if(plugin->call != NULL) plugin->call

#endif

