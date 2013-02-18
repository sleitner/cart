#ifndef __FEEDBACK_SNII_H__
#define __FEEDBACK_SNII_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

void snII_config_init();
void snII_config_verify();
void snII_init();
void snII_setup(int level);

#if defined(HYDRO) && defined(PARTICLES)
void snII_thermal_feedback(int level, int cell, int ipart, double t_next );
void snII_kinetic_feedback(int level, int cell, int ipart, double t_next );
#endif /* HYDRO && PARTICLES */

struct SNII_t
{
  double energy_per_explosion;     /* used to be called E_51 */
  double kinetic_energy_per_explosion;   
  double time_duration;            /* used to be called t_fb */
  double time_delay;               
  double yield_factor;             /* fraction yield relative to the one coded in */
  double min_mass;                 /* used to be called aM_SNII */
  double max_mass;                 
};
struct SNII_PROP_t
{
  double energy;
  double kinetic_energy;
  double metals;
  double teject;
  double tdelay;
  double fmass;
};

extern struct SNII_t snII;
extern struct SNII_PROP_t snII_phys, snII_code;
#endif /* STAR_FORMATION */
#endif
