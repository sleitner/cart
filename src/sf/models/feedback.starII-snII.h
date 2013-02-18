#ifndef __FEEDBACK_STARII_SNII_H__
#define __FEEDBACK_STARII_SNII_H__

#ifdef STARFORM
void starII_explosion_config_init();
void starII_explosion_config_verify();
void starII_explosion_feedback_init();
void starII_explosion_setup(int level);

void starII_explosion_mass(int level, int icell, int ipart);
void starII_explosion_thermal(int level, int cell, int ipart);
void starII_explosion_kicks(int level, int cell, int ipart);
#endif /* STARFORM */

#endif /* __FEEDBACK_STARII_SNII_H__ */
