#ifndef __REFINEMENT_INDICATORS_H__
#define __REFINEMENT_INDICATORS_H__

#define num_refinement_indicators	7

#define DARK_MASS_INDICATOR             0
#define GAS_MASS_INDICATOR              1
#define SHOCK_INDICATOR                 2
#define CONTACT_DISCONTINUITY_INDICATOR 3
#define DENSITY_GRADIENT_INDICATOR      4
#define PRESSURE_GRADIENT_INDICATOR     5
#define ENTROPY_GRADIENT_INDICATOR      6

void mark_refinement_indicators( int cell, int level );

float dark_mass_indicator( int cell, int level );
float gas_mass_indicator( int cell, int level );
float shock_indicator( int cell, int level );
float contact_discontinuity_indicator( int cell, int level );
float density_gradient_indicator( int cell, int level );
float pressure_gradient_indicator( int cell, int level );
float entropy_gradient_indicator( int cell, int level );

extern float refinement_indicator_threshold[num_refinement_indicators][num_refinement_levels];
extern float refinement_indicator_weight[num_refinement_indicators];
extern int use_refinement_indicator[num_refinement_indicators][num_refinement_levels];
extern float refinement_volume_min[nDim];
extern float refinement_volume_max[nDim];

#endif
