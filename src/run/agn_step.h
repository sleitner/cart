#ifndef __AGN_STEP_H__
#define __AGN_STEP_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef AGN 

#include "halos.h"

void agn_feedback(int level);
void agn_find_mergers();
void agn_seed( halo_list *list );

#endif /* AGN */

#endif /* __AGN_STEP_H__ */
