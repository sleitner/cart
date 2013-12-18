#include "config.h"
#if defined(STAR_FORMATION) && defined(HYDRO)

#include <math.h>

#include "auxiliary.h"
#include "iterators.h"
#include "parallel.h"
#include "particle.h"
#include "rand.h"
#include "starformation.h"
#include "starformation_recipe.h"
#include "starformation_formstar.h"
#include "starformation_feedback.h"
#include "times.h"
#include "timing.h"
#include "tree.h"
#include "units.h"

#include "starformation_step.h"
#include "step.h"

#ifdef LOG_STAR_CREATION
#include "logging.h"
#endif

/* star formation parameters */
extern int sf_min_level;

void star_formation(int level, int time_multiplier) {
  int i, j;
  int icell;
  int num_level_cells;
  int *level_cells;
  double dt_eff;
  float *sfr;

  if ( level < sf_min_level ) return;

  start_time( WORK_TIMER );

  select_level( level, CELL_TYPE_LOCAL | CELL_TYPE_LEAF,  &num_level_cells, &level_cells );

  sfr = cart_alloc(float,num_level_cells);
  star_formation_rate(level,num_level_cells,level_cells,sfr);
  dt_eff = dtl[level] * time_multiplier;

  if(sf_formstar->setup != NULL) sf_formstar->setup(level);

/* cannot be done shared-memory parallel */
/* #pragma omp parallel for default(none), private(i,icell), shared(num_level_cells,level_cells,level,sf_formstar,dt_eff, sfr, dtl), schedule(dynamic) */
  for(i=0;i<num_level_cells; i++)
    {
        if ( sfr[i] <=0 ) continue;
        icell = level_cells[i];
        sf_formstar->form_star_particles(level,icell,dtl[level],dt_eff,sfr[i]);
    }

  cart_free(level_cells);
  cart_free(sfr);
  end_time( WORK_TIMER );
}

#endif /* STAR_FORMATION && HYDRO */
