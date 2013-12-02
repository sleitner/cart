#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)


#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "starformation_formstar.h"

extern struct FormStar sf_formstar_internal;
const struct FormStar *sf_formstar = &sf_formstar_internal;

/*
//  Configuration
*/
void control_parameter_list_formstar(FILE *stream, const void *ptr)
{
  fprintf(stream,"<%s>",sf_formstar->name);
}


void config_init_formstar()
{
  static char *ptr;
  ControlParameterOps r = { NULL, control_parameter_list_formstar };

  ptr = cart_alloc(char,strlen(sf_formstar_internal.name)+1);
  strcpy(ptr,sf_formstar_internal.name);
  control_parameter_add(r,ptr,"sf:formstar","a method for forming star particles. This parameter is for listing only, and must be set with SF_FORMSTAR define in defs.h. See /src/sf for available methods.");

  if(sf_formstar_internal.config_init != NULL) sf_formstar_internal.config_init();
}


#define STR_VALUE(arg)      #arg
#define to_string(name)     STR_VALUE(name)
void config_verify_formstar()
{
  char formstar_internal_name[99];
#ifdef SF_FORMSTAR
  const char *formstar_external_name = to_string(SF_FORMSTAR);
#else
  const char *formstar_external_name = "";
#endif

  cart_assert(sf_formstar_internal.name != NULL);
  cart_assert(sf_formstar_internal.form_star_particles != NULL);

  sprintf(formstar_internal_name,"<%s>",sf_formstar_internal.name);
  if(strcmp("<custom>",formstar_external_name)!=0 && strcmp(formstar_internal_name,formstar_external_name)!=0)
    {
      cart_error("Misconfiguration: the internal SF formstar name (%s) does not match the name set in defs.h (%s)",formstar_internal_name,formstar_external_name);
    }

  VERIFY(sf:formstar, 1 );

  if(sf_formstar_internal.config_verify != NULL) sf_formstar_internal.config_verify();
}
#undef STR_VALUE
#undef to_string

void init_formstar()
{
  if(sf_formstar_internal.init != NULL) sf_formstar_internal.init ();
}

#endif /* HYDRO && STAR_FORMATION */
