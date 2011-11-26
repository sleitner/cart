#ifndef __STARFORMATION_RECIPES_H__
#define __STARFORMATION_RECIPES_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#if defined(HYDRO) && defined(STARFORM)

struct StarFormationRecipe
{
  void (*setup)(int level);
  double (*rate)(int cell);
  int id;
  const char *name;
};

extern const struct StarFormationRecipe *sf_recipe;

void config_init_star_formation_recipes();
void config_verify_star_formation_recipes();

#endif /* HYDRO && STARFORM */

#endif
