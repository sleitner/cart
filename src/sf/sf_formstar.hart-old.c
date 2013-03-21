#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include "starformation.h"

#include "form_star.hart-old.h"

void star_form_config_init()
{
  hart_old_config_init();
}


void star_form_config_verify()
{
  hart_old_config_verify();
}


void star_form_setup(int level)
{
  hart_old_setup(level);
}

void star_form_particles(int level, int icell, double dtl, double dt, float sfr){
  hart_old_star_formation( level, icell, dtl, dt, sfr );
}

struct FormStar sf_formstar_internal = 
  {
    "hart-old\0",
    star_form_particles,
    star_form_config_init,
    star_form_config_verify,
    star_form_setup
  };

#endif /* HYDRO && STAR_FORMATION */
