#ifndef __SF_STARFORM_H__
#define __SF_STARFORM_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

#ifdef STAR_FORMATION
extern double tdelay_popM_feedback;
extern int poissonRF12_starformation_indicator;
extern int continuous_starformation_indicator;

void star_form_config_init();
void star_form_config_verify();
void star_form_setup(int level);
void star_form_particles(int level, int icell, double dtl, double dt, float sfr);
#endif /* STAR_FORMATION */

#endif
