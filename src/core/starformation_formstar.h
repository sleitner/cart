#ifndef __STAR_FORMATION_FORMSTAR_H__
#define __STAR_FORMATION_FORMSTAR_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#if defined(HYDRO) && defined(STAR_FORMATION)

/*
//  ATTENTION DEVELOPERS:
//  ONLY add new members at the end of the structure!!!
//  ONLY add new members if they are inserted in a new place in the code!!!
*/ 
struct FormStar
{
  const char *name;
  void (*form_star_particles)(int level, int icell, double dtl, double dt, float sfr); 
  void (*config_init)();           /* can be NULL */
  void (*config_verify)();         /* can be NULL */
  void (*init)();                  /* can be NULL */
  void (*setup)(int level);        /* can be NULL */
};

extern const struct FormStar *sf_formstar;

void config_init_formstar();
void config_verify_formstar();
void init_formstar();
void setup_formstar(int level);

#endif /* HYDRO && STAR_FORMATION */

#endif
