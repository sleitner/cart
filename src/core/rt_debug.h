#ifndef __RT_DEBUG_H__
#define __RT_DEBUG_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef RADIATIVE_TRANSFER

#ifndef RT_CONFIGURED
#error "Missing rt_config.h include."
#endif


#ifdef RT_DEBUG
/*
//  Debugging extensions
*/
struct rtDebugData
{
  int Mode;
  double Pos[3];
};

extern struct rtDebugData rt_debug;

#endif /* RT_DEBUG */
#endif /* RADIATIVE_TRANSFER */

#endif /* __RT_DEBUG_H__ */

