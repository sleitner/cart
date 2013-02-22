#ifndef __STAR_FORMATION_RECIPES_H__
#define __STAR_FORMATION_RECIPES_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#if defined(HYDRO) && defined(STAR_FORMATION)

/*
//  ATTENTION DEVELOPERS:
//  ONLY add new members at the end of the structure!!!
//  ONLY add new members if they are inserted in a new place in the code!!!
*/ 
struct StarFormationRecipe
{
  const char *name;
  double (*rate)(int cell);
  void (*config_init)();           /* can be NULL */
  void (*config_verify)();         /* can be NULL */
  void (*setup)(int level);        /* can be NULL */
  void (*form_star_particles)(int level, int icell, double dt, float sfr); /* NULL */
};

extern const struct StarFormationRecipe *sf_recipe;

void config_init_star_formation_recipe();
void config_verify_star_formation_recipe();

#endif /* HYDRO && STAR_FORMATION */

#endif
