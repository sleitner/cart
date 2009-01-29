#ifndef __RT_DEBUG_H__
#define __RT_DEBUG_H__


#include <stdio.h>


/*
//  Debugging extensions
*/
struct rtDebugData
{
  int Mode;
  int Stop;
  double Pos[3];
};

extern struct rtDebugData rt_debug;

#endif  /* __RT_DEBUG_H__ */

