#include "config.h"
#if defined(HYDRO) && defined(STAR_FORMATION)

#include "form_star.poissonRF12.h"
#include "form_star.continuous.h"
#include "form_star.all.h"

void star_form_config_init()
{
  poissonRF12_config_init();
  continuous_config_init();
}


void star_form_config_verify()
{
  poissonRF12_config_verify();
  continuous_config_verify();
}


void star_form_setup(int level)
{
  poissonRF12_setup(level);
  continuous_setup(level);
}

void star_form_particles(int level, int icell, double dtl, double dt, float sfr){
    poissonRF12_star_formation( level, icell, dtl, dt, sfr );
    continuous_star_formation( level, icell, dtl, dt, sfr );
}

#endif /* HYDRO && STAR_FORMATION */
