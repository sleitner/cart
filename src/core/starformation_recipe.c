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
void control_parameter_list_recipe(FILE *stream, const void *ptr)
{
  fprintf(stream,"<%s>",sf_recipe->name);
}


void config_init_star_formation_recipe()
{
  static char *ptr;
  ControlParameterOps r = { NULL, control_parameter_list_recipe };

  ptr = cart_alloc(char,strlen(sf_recipe_internal.name)+1);
  strcpy(ptr,sf_recipe_internal.name);
  control_parameter_add(r,ptr,"sf:recipe","a recipe for star formation. This parameter is for listing only, and must be set with SF_RECIPE define in defs.h. See /src/sf for available recipes.");

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
  if(strcmp("<custom>",recipe_external_name)!=0 && strcmp(recipe_internal_name,recipe_external_name)!=0)
    {
      cart_error("Misconfiguration: the internal SF recipe name (%s) does not match the name set in defs.h (%s)",recipe_internal_name,recipe_external_name);
    }

  VERIFY(sf:recipe, 1 );

  if(sf_recipe_internal.config_verify != NULL) sf_recipe_internal.config_verify();
}
#undef STR_VALUE
#undef to_string

#endif /* HYDRO && STAR_FORMATION */
