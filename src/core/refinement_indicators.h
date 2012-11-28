#ifndef __REFINEMENT_INDICATORS_H__
#define __REFINEMENT_INDICATORS_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef REFINEMENT

#define num_refinement_indicators	11

#define DARK_MASS_INDICATOR             0
#define GAS_MASS_INDICATOR              1
#define SHOCK_INDICATOR                 2
#define CONTACT_DISCONTINUITY_INDICATOR 3
#define DENSITY_GRADIENT_INDICATOR      4
#define PRESSURE_GRADIENT_INDICATOR     5
#define ENTROPY_GRADIENT_INDICATOR      6
#define SPATIAL_INDICATOR               7
#define JEANS_LENGTH_INDICATOR          8
#define DARK_1STSPEC_INDICATOR          9
#define GAS_1STSPEC_INDICATOR           10


typedef struct REFINEMENT_INDICATOR_TYPE
{
  int *use;
  int use_buffer[max_level-min_level+1];
  float weight;
  float *threshold;
  float threshold_buffer[max_level-min_level+1];
}
refinement_t;
extern refinement_t refinement_indicator[num_refinement_indicators];


void config_init_refinement_indicators();
void config_verify_refinement_indicators();

void mark_refinement_indicators( int cell, int level );

float dark_mass_indicator( int cell, int level );
float gas_mass_indicator( int cell, int level );
float dark_1stspec_indicator( int cell, int level );
float gas_1stspec_indicator( int cell, int level );
float spatial_indicator( int cell, int level );
float shock_indicator( int cell, int level, int neighbors[] );
float contact_discontinuity_indicator( int cell, int level, int neighbors[], float drho[] );
float density_gradient_indicator( int cell, int level, int neighbors[], float drho[] );
float pressure_gradient_indicator( int cell, int level, int neighbors[] );
float entropy_gradient_indicator( int cell, int level, int neighbors[] );
float jeans_length_indicator( int cell, int level );

#endif /* REFINEMENT */

#endif
