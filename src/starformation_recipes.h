#ifndef __STARFORMATION_RECIPES_H__
#define __STARFORMATION_RECIPES_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STARFORM

typedef double (*StarFormationRate)(int cell);

extern StarFormationRate sf_rate;

void config_init_star_formation_recipes();
void config_verify_star_formation_recipes();

void setup_star_formation_recipes(int level);

#endif /* STARFORM */

#endif
