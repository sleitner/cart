#ifndef __STARFORMATION_RECIPES_H__
#define __STARFORMATION_RECIPES_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#if defined(HYDRO) && defined(STARFORM)

struct StarFormationRecipe
{
  const char *name;
  void (*config_init)();
  void (*config_verify)();
  void (*level_setup)(int level);
  double (*rate)(int cell);
};

extern const struct StarFormationRecipe *sf_recipe;

void config_init_star_formation_recipe();
void config_verify_star_formation_recipe();

#endif /* HYDRO && STARFORM */

#endif
