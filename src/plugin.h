#ifndef __PLUGIN_H__
#define __PLUGIN_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef USER_PLUGIN

#include <mpi.h>

/*
//  This file contains prototypes for user plugin functions.
//  These functions are combined in a single structure to
//  simplify extensibility: if a new plugin point is created,
//  it needs to be set to NULL by default. That way old codes
//  do NOT need to be modified. All non-NULL functions need to be 
//  defined and linked in at compile time.
*/
typedef struct PLUGIN_TYPE
{
  void (*RunBegin)();
  void (*RunEnd)();
  void (*GlobalStepBegin)();
  void (*GlobalStepEnd)();
  void (*AfterCFLRestart)();
  void (*LevelStepBegin)(int level, MPI_Comm level_com);
  void (*LevelStepEnd)(int level, MPI_Comm level_com);
  void (*LevelStepFail)(int level, MPI_Comm level_com);
    void (*StarformationFeedbackEnd)(int level, int cell, double dU);
}
plugin_t;

struct PLUGIN_LIST_TYPE
{
  int index;
  const plugin_t **list;
};

extern struct PLUGIN_LIST_TYPE plugins;

#define PLUGIN_POINT(call) \
  for(plugins.index=0; plugins.list[plugins.index]!=NULL; plugins.index++) if(plugins.list[plugins.index]->call != NULL) plugins.list[plugins.index]->call

const plugin_t* add_plugin(int id);

#endif /* USER_PLUGIN */

#endif

