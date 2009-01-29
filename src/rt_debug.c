#include "defs.h"      
#ifdef RADIATIVE_TRANSFER
#include "rt_config.h"

#ifdef RT_DEBUG

#include "rt_debug.h"

/*
//  Debug modes:
//  1. Cell-level debug
//  2. Dump atomic properties at max_level
//
*/

struct rtDebugData rt_debug = { 0, 0, { 0.0, 0.0, 0.0 } };


#endif  /* RT_DEBUG */
#endif  /* RADIATIVE_TRANSFER */
