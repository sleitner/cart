#ifndef __STAR_FORMATION_RECIPES_H__
#define __STAR_FORMATION_RECIPES_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#if defined(HYDRO) && defined(STAR_FORMATION)

struct StarFormationFeedback;

struct StarFormationRecipe
{
  const char *name;
  double (*rate)(int cell);
  void (*config_init)();           /* can be NULL */
  void (*config_verify)();         /* can be NULL */
  void (*setup_feedback)();        /* can be NULL */
  void (*setup_level)(int level);  /* can be NULL */
};

extern const struct StarFormationRecipe *sf_recipe;

void config_init_star_formation_recipe();
void config_verify_star_formation_recipe();

#endif /* HYDRO && STAR_FORMATION */

#endif
