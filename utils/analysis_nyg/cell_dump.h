#ifndef __CELL_DUMP_H__
#define __CELL_DUMP_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include "extra/ism.h"


/*
//  Class used for representing what we dump for a given cell
*/
struct ngCellDump
{
  int Size;
  DumpWorker Worker;
  const int *Weight;
  const char **Header;
};

#endif  /* __CELL_DUMP_H__ */
