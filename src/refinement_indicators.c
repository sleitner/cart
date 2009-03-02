#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "defs.h"
#include "tree.h"
#include "refinement_indicators.h"

/*
//  NG: +1 here is so that num_refinement_levels can be set to zero 
//      for a uniform mesh mode
*/
float refinement_indicator_threshold[num_refinement_indicators][num_refinement_levels+1];
float refinement_indicator_weight[num_refinement_indicators+1];
int use_refinement_indicator[num_refinement_indicators][num_refinement_levels+1];

float refinement_volume_min[nDim];
float refinement_volume_max[nDim];

int neighbors[num_neighbors];
float drho[nDim];

void mark_refinement_indicators( int cell, int level ) {
	int i;
        float indicator = 0.0;

	if ( use_refinement_indicator[DARK_MASS_INDICATOR][level] ) {
	        indicator = max( dark_mass_indicator(cell, level), indicator );
	}

#ifdef HYDRO
	if ( use_refinement_indicator[GAS_MASS_INDICATOR][level] ) {
		indicator = max( gas_mass_indicator( cell, level ), indicator );
	}

	cell_all_neighbors( cell, neighbors );
	for ( i = 0; i < nDim; i++ ) {
		drho[i] = fabs( cell_gas_density( neighbors[2*i] ) 
			- cell_gas_density( neighbors[2*i+1] ) ) / 
			min ( cell_gas_density( neighbors[2*i] ), 
				cell_gas_density( neighbors[2*i+1] ) );
	}

	if ( use_refinement_indicator[SHOCK_INDICATOR][level] ) {
		indicator = max( shock_indicator(cell, level), indicator );
	}

	if ( use_refinement_indicator[CONTACT_DISCONTINUITY_INDICATOR][level] ) {
		indicator = max( contact_discontinuity_indicator(cell,level), indicator );
	}

	if ( use_refinement_indicator[DENSITY_GRADIENT_INDICATOR][level] ) {
		indicator = max( density_gradient_indicator(cell,level), indicator );
	}

	if ( use_refinement_indicator[PRESSURE_GRADIENT_INDICATOR][level] ) {
		indicator = max( pressure_gradient_indicator(cell,level), indicator );
	}

	if ( use_refinement_indicator[ENTROPY_GRADIENT_INDICATOR][level] ) {
		indicator = max( entropy_gradient_indicator(cell,level), indicator );
	}
#endif /* HYDRO */

        refinement_indicator(cell, 0) = indicator;
}

float dark_mass_indicator( int cell, int level ) {
	float ave_mass;

#ifdef PARTICLES
	ave_mass = (cell_first_species_mass(cell)) /
			refinement_indicator_threshold[DARK_MASS_INDICATOR][level];
#else
	ave_mass = 0.0;
#endif /* PARTICLES */

	return min( ave_mass, refinement_indicator_weight[DARK_MASS_INDICATOR] );
}

#ifdef HYDRO

float gas_mass_indicator( int cell, int level ) {
	float ave_mass;

	ave_mass = ( cell_volume[level] * cell_gas_density(cell) ) / refinement_indicator_threshold[GAS_MASS_INDICATOR][level];
	return min( ave_mass, refinement_indicator_weight[GAS_MASS_INDICATOR] );
}

float shock_indicator( int cell, int level ) {
	int i;
	float dp;
	float dv;
	float indicator = 0.0;

	for ( i = 0; i < nDim; i++ ) {
		if ( cell_level( neighbors[2*i] ) == level 
				&& cell_level( neighbors[2*i+1] ) == level
				&& cell_level( cell_neighbor( neighbors[2*i], i ) ) == level
				&& cell_level( cell_neighbor( neighbors[2*i+1], i+1 ) ) == level ) {
			dp = fabs( cell_gas_pressure( neighbors[2*i] ) 
				- cell_gas_pressure( neighbors[2*i+1] ) ) / 
				min( cell_gas_pressure( neighbors[2*i] ),
					cell_gas_pressure( neighbors[2*i+1] ) );

			dv = fabs( cell_momentum(neighbors[2*i],i) ) / cell_gas_density( neighbors[2*i] )
				- cell_momentum(neighbors[2*i+1],i) / cell_gas_density( neighbors[2*i+1] );

			if ( dp > 3.0 && dv > 0.0 && 
					cell_gas_density(cell) > refinement_indicator_threshold[SHOCK_INDICATOR][level] ) {
				indicator = max( refinement_indicator_weight[SHOCK_INDICATOR], indicator );
			}
		}
	}

	return indicator;
}

float contact_discontinuity_indicator( int cell, int level ) {
        int i;
        float dp;
        float indicator = 0.0;
                                                                                             
        for ( i = 0; i < nDim; i++ ) {
		if ( cell_level( neighbors[2*i] ) == level
				&& cell_level( neighbors[2*i+1] ) == level
				&& cell_level( cell_neighbor( neighbors[2*i], i ) ) == level
				&& cell_level( cell_neighbor( neighbors[2*i+1], i+1 ) ) == level ) {
	                dp = fabs( cell_gas_pressure( neighbors[2*i] )
        	                - cell_gas_pressure( neighbors[2*i+1] ) ) /
                	        min( cell_gas_pressure( neighbors[2*i] ),
                        	        cell_gas_pressure( neighbors[2*i+1] ) );
                                                                                             
	                if ( dp < 10.0 && drho[i] >= 0.6 ) {
        	                indicator = max( refinement_indicator_weight[CONTACT_DISCONTINUITY_INDICATOR], indicator );
	                }
		}
        }
                                                                                             
        return indicator;
}

float density_gradient_indicator( int cell, int level ) {
	int i;
	float rho_ind1, rho_ind2;
	float drho_indicator;
	float indicator = 0.0;

	for ( i = 0; i < nDim; i++ ) {
		rho_ind1 = fabs( cell_gas_density(cell) - cell_gas_density(neighbors[2*i]) ) /
			min ( cell_gas_density(cell), cell_gas_density(neighbors[2*i] ) );

		rho_ind2 = fabs(cell_gas_density(cell) - cell_gas_density(neighbors[2*i+1])) /
			min ( cell_gas_density(cell), cell_gas_density(neighbors[2*i+1] ) );

		drho_indicator = max( drho[i], max( rho_ind1, rho_ind2 ) );
		if ( drho_indicator > 1.0 && 
			cell_gas_density(cell) > refinement_indicator_threshold[DENSITY_GRADIENT_INDICATOR][level] ) {
			indicator = max( indicator, refinement_indicator_weight[DENSITY_GRADIENT_INDICATOR]*max(rho_ind1, rho_ind2) );
		}
	}

	return indicator;
}

float pressure_gradient_indicator( int cell, int level ) {
	int i;
        float press_ind1, press_ind2;
        float indicator = 0.0;
                                                                                             
        for ( i = 0; i < nDim; i++ ) {
		if ( cell_level( neighbors[2*i] ) == level 
				&& cell_level( neighbors[2*i+1] ) == level
				&& cell_level( cell_neighbor( neighbors[2*i], 2*i ) ) == level 
				&& cell_level( cell_neighbor( neighbors[2*i+1], 2*i+1 ) ) == level ) {
			press_ind1 = fabs( cell_gas_pressure(cell) - cell_gas_pressure(neighbors[2*i]) )
				/ ( cell_gas_pressure(cell) + cell_gas_pressure(neighbors[2*i]) );
                                                                                             
	                press_ind2 = fabs(cell_gas_pressure(cell) - cell_gas_pressure(neighbors[2*i+1])) 
				/ ( cell_gas_pressure(cell) + cell_gas_pressure(neighbors[2*i+1] ) );
                                                                                             
	                if ( cell_gas_density(cell) > refinement_indicator_threshold[PRESSURE_GRADIENT_INDICATOR][level] ) {
        	                indicator = max( indicator, 
					refinement_indicator_weight[PRESSURE_GRADIENT_INDICATOR]*max(press_ind1, press_ind2) );
                	}
		}
        }
                                                                                             
        return indicator;

}

float entropy_gradient_indicator( int cell, int level ) {
	return 0.0;
}

#endif /* HYDRO */
