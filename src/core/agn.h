#ifndef __AGN_H__
#define __AGN_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef AGN

void config_init_agn();
void config_verify_agn();
void init_agn();
void agn_feedback(int level );
void agn_find_mergers();

#endif /* AGN */

#endif 
