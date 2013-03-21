#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include "starformation.h"

#include "form_star.hart.h"

void star_form_config_init()
{
  hart_config_init();
}


void star_form_config_verify()
{
  hart_config_verify();
}


void star_form_setup(int level)
{
  hart_setup(level);
}

void star_form_particles(int level, int icell, double dtl, double dt, float sfr){
  hart_star_formation( level, icell, dtl, dt, sfr );
}

struct FormStar sf_formstar_internal = 
  {
    "hart\0",
    star_form_particles,
    star_form_config_init,
    star_form_config_verify,
    star_form_setup
  };

#endif /* HYDRO && STAR_FORMATION */
