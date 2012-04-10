#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)


#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "starformation_recipe.h"



extern struct StarFormationRecipe sf_recipe_internal;
const struct StarFormationRecipe *sf_recipe = &sf_recipe_internal;


/*
//  Configuration
*/
void control_parameter_set_recipe(const char *value, void *ptr, int ind)
{
  cart_error("Star formation recipe is specified at compile-time and cannot be set in the configuration file.");
}


void control_parameter_list_recipe(FILE *stream, const void *ptr)
{
  fprintf(stream,"<%s>",sf_recipe->name);
}


void config_init_star_formation_recipe()
{
  ControlParameterOps r = { control_parameter_set_recipe, control_parameter_list_recipe };

  control_parameter_add(r,(void *)sf_recipe_internal.name,"sf:recipe","a recipe for star formation. See /src/sf for available recipes.");

  if(sf_recipe_internal.config_init != NULL) sf_recipe_internal.config_init();
}


#define STR_VALUE(arg)      #arg
#define to_string(name)     STR_VALUE(name)

void config_verify_star_formation_recipe()
{
  char recipe_internal_name[99];
#ifdef SF_RECIPE
  const char *recipe_external_name = to_string(SF_RECIPE);
#else
  const char *recipe_external_name = "";
#endif

  cart_assert(sf_recipe_internal.name != NULL);
  cart_assert(sf_recipe_internal.rate != NULL);

  sprintf(recipe_internal_name,"<%s>",sf_recipe_internal.name);
  cart_assert(strcmp(recipe_internal_name,recipe_external_name) == 0);

  VERIFY(sf:recipe, 1 );

  if(sf_recipe_internal.config_verify != NULL) sf_recipe_internal.config_verify();
}

#endif /* HYDRO && STAR_FORMATION */
