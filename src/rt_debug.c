#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(RT_DEBUG)

#include "rt_debug.h"

/*
//  Debug modes:
//  1. Cell-level debug
//  2. Dump atomic properties at max_level
//
*/

struct rtDebugData rt_debug = { 0, 0, { 0.0, 0.0, 0.0 } };


#endif /* RADIATIVE_TRANSFER && RT_DEBUG */
