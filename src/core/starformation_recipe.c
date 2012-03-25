#include "config.h"
#if defined(HYDRO) && defined(STARFORM)


#include "starformation_recipe.h"



extern struct StarFormationRecipe sf_recipe_internal;
const struct StarFormationRecipe *sf_recipe = &sf_recipe_internal;


/*
//  Configuration
*/
void config_init_star_formation_recipe()
{
  sf_recipe_internal.config_init();
}


void config_verify_star_formation_recipe()
{
  sf_recipe_internal.config_verify();
}

#endif /* HYDRO && STARFORM */
