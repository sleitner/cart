#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)


#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "starformation_recipe.h"



extern struct StarFormationRecipe sf_recipe_internal;
const struct StarFormationRecipe *sf_recipe = &sf_recipe_internal;

int poissonRF12_starformation_indicator = 1;
int continuous_starformation_indicator = 0;

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

  control_parameter_add2(control_parameter_int,&poissonRF12_starformation_indicator,"sfRF12:indicator","poissonRF12_starformation_indicator","create star particles using a poisson sampling around a mass that is typically formed over sfRF12_timescale for the given star formation rate");

  control_parameter_add2(control_parameter_int,&continuous_starformation_indicator,"sfcontinuous:indicator","continuous_starformation_indicator","create star particles and then grow them using the star formation rate.");

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

  VERIFY(sfcontinuous:indicator, continuous_starformation_indicator==0 || continuous_starformation_indicator==1 );
  VERIFY(sfRF12:indicator, poissonRF12_starformation_indicator==0 || poissonRF12_starformation_indicator==1 );
  if(   !(poissonRF12_starformation_indicator || poissonRF12_starformation_indicator) 
      || (poissonRF12_starformation_indicator && poissonRF12_starformation_indicator) ) 
      {
	  cart_error("One and only one star creation function should be used");
      }

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

#endif /* HYDRO && STAR_FORMATION */
