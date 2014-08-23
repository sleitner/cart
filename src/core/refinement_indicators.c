#include "config.h"
#ifdef REFINEMENT

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "refinement_indicators.h"
#include "tree.h"
#include "plugin.h"

#ifdef PARTICLES
#include "particle.h"
extern float particle_species_mass[MAX_PARTICLE_SPECIES];
#endif /* PARTICLES */

/*
//  List of units for the 
*/


refinement_t refinement_indicator[num_refinement_indicators];


void control_parameter_set_refinement_indicator(const char *value, void *ptr, int ind)
{
  char *str, *tok;
  int id, from_level, to_level, level, n;
  float w;
  char c, unit[10];

  str = cart_alloc(char,strlen(value)+1);
  strcpy(str,value); /* strtok destroys the input string */

  /*
  //  indicator id
  */
  tok = strtok(str," ");
  if(sscanf(tok,"%d",&id)!=1 && sscanf(tok,"id=%d",&id)!=1)
    {
      cart_error("Unable to read refinement ID parameter from token '%s' of string '%s'",tok,value);
    }
  
  if(id<-1 || id>= num_refinement_indicators)
    {
      cart_error("Invalid indicator id %d in '%s'",id,value);
    }

  if(id == -1)
    {
      /*
      //  That means the indicator is not set - just skip
      */
      return;
    }
 
  /*
  // weight
  */
  tok = strtok(NULL," ");
  if(sscanf(tok,"%g",&w)!=1 && sscanf(tok,"weight=%g",&w)!=1)
    {
      cart_error("Unable to read refinement WEIGHT parameter from token '%s' of string '%s'",tok,value);
    }
  if(w < 0.0)
    {
      cart_error("Refinement weight for indicator %d must be positive",id,w);
    }
  refinement_indicator[id].weight = w;

  /*
  // levels
  */
  tok = strtok(NULL," ");
  if(sscanf(tok,"%d",&from_level)!=1 && sscanf(tok,"from-level=%d",&from_level)!=1)
    {
      cart_error("Unable to read refinement FROM-LEVEL parameter from token '%s' of string '%s'",tok,value);
    }
  if(from_level < min_level)
    {
      cart_error("<from-level> for indicator %d cannot be less than min_level (%d)",id,min_level);
    }
  if(from_level >= max_level && max_level > min_level )
    {
      cart_error("<from-level> for indicator %d must be less than max_level (%d)",id,max_level);
    }

  tok = strtok(NULL," ");
  if(sscanf(tok,"%d",&to_level)!=1 && sscanf(tok,"to-level=%d",&to_level)!=1)
    {
      cart_error("Unable to read refinement TO-LEVEL parameter from token '%s' of string '%s'",tok,value);
    }
  if(to_level < from_level)
    {
      cart_error("<to-level> for indicator %d cannot be less than <from-level> (%d)",id,from_level);
    }
  if(to_level > max_level) to_level = max_level;

  /*
  // threshold values: at least one value is required
  */
  for(level=from_level; level<=to_level; )
    {
      tok = strtok(NULL," ");
      if(tok==NULL && level>from_level)
	{
	  /*
	  // Use the last value for all missing tokens: the last value must be present
	  */
	  w = refinement_indicator[id].threshold[level-1];
	  for(; level<=to_level; level++)
	    {
	      refinement_indicator[id].use[level] = 1;
	      refinement_indicator[id].threshold[level] = w;
	    }
	}
      else
	{
	  n = 1;
	  unit[0] = 0;
	  if(sscanf(tok,"%g%c",&w,&c)!=1 && sscanf(tok,"%d*%g%c",&n,&w,&c)!=2 && sscanf(tok,"%d*%g[%s]%c",&n,&w,unit,&c)!=3)
	    {
	      cart_error("Unable to read refinement THRESHOLD from token '%s' of string '%s'",tok,value);
	    }
	  if(n < 1)
	    {
	      cart_error("Numerical multiplier %d in the threshold value for indicator %d must be positive",n,id);
	    }
	  /*
	  //  Parse unit identifier
	  */
	  if(unit[0] != 0)
	    {
	      if(strcmp(unit,"M0") == 0)
		{
		  /*
		  //  Multiply the value by the mean specie mass
		  */
		  switch(id)
		    {
		    default:
		      {
			cart_error("Unit <M0> is not valid for refinement indicator %d",id);
		      }
		    }
		}
	      else
		{
		  cart_error("String '%s' is not a valid refinement indicator unit",unit);
		}
	    }

	  for(; level<=to_level && n>0; level++, n--)
	    {
	      refinement_indicator[id].use[level] = 1;
	      refinement_indicator[id].threshold[level] = w;
	    }
	}
    }

  cart_free(str);
}


void control_parameter_print_name(FILE *f, const char *name);

void control_parameter_list_refinement_indicator(FILE *stream, const void *ptr)
{
  int i, from_level, to_level, level;
  int newline = 0;

  for(i=0; i<num_refinement_indicators; i++)
    {
      from_level = min_level;
      while(from_level < max_level)
	{
	  if(refinement_indicator[i].use[from_level])
	    {
	      to_level = from_level + 1;
	      while(to_level<max_level && refinement_indicator[i].use[to_level]) to_level++;
	      to_level--;
	      if(newline)
		{
		  fprintf(stream,"\n");
		  control_parameter_print_name(stream,"ref:indicator");
		}
	      fprintf(stream,"id=%-d weight=%-g from-level=%-d to-level=%-d ",i,refinement_indicator[i].weight,from_level,to_level);
	      for(level=from_level; level<=to_level; level++)
		{
		  fprintf(stream,"%g ",refinement_indicator[i].threshold[level]);
		}
	      from_level = to_level;
	      newline = 1;
	    }
	  from_level++;
	}
    }

  if(newline == 0)
    {
      /*
      //  Refinement indicators are not set
      */
      fprintf(stream,"-1 (NOT SET)");
    }
}


void config_init_refinement_indicators()
{
  int i, level;
  ControlParameterOps control_parameter_refinement_indicator = { control_parameter_set_refinement_indicator, control_parameter_list_refinement_indicator };

  /*
  //  Init everything first
  */
  for(i=0; i<num_refinement_indicators; i++)
    {
      refinement_indicator[i].weight = 0.0;
      refinement_indicator[i].use = refinement_indicator[i].use_buffer - min_level;
      refinement_indicator[i].threshold = refinement_indicator[i].threshold_buffer - min_level;

      for(level=min_level; level<=max_level; level++)
	{
	  refinement_indicator[i].use[level] = 0;
	  refinement_indicator[i].threshold[level] = 0.0;
	}
    }

  control_parameter_add2(control_parameter_refinement_indicator,refinement_indicator,"ref:indicator","refinement","set one of the refinement indicators. This command may appear multiple times in the config file (to set several refinement indicators or to set a complex refinement pattern to some of them). The full format of this command is\n  refinement id=<id> weight=<weight> from-level=<from-level>\n  to-level=<to-level> [multiplier*]value [ [multiplier*]value ... ].\nHere <id> is the index of the indicator (check the include file 'refinement_indicators.h' for the list of available indicators), <weight> is the relative weight of this indicator (must be between 0 and 1), <from-level> and <to-level> is the range of levels to which the current setting applies, and the rest of the line lists the values of refinement thresholds. There should be exactly (<to-level>-<from-level>+1) values listed on the line, with 2 exceptions: (i) a single value may be prefixed with an integer multiplier to short-hand several identical values ('3*0.5' is equivalent to '0.5 0.5 0.5'), and (ii) if there are not enough values listed, the last value is used to set the missing thresholds. For example, the following command:\n  refinement id=0 weight=0.5 from-level=0 to-level=9 2*0.5 0.2\nsets 10 threshold values for levels from 0 to 9: levels 0 and 1 are set the value of 0.5, and the value of 0.2 is used to set all the values on levels 2 to 9. Be advised that setting the refinement indicator thresholds is a delicate task - thresholds for some indicators are dimensional and understanding of the approriate units is required to set the thresholds correctly.");

}


void config_verify_refinement_indicators()
{
  /*
  //  All checking is already done in control_parameter_set_refinement_indicator
  */
  VERIFY(ref:indicator, 1 );
}


void mark_refinement_indicators( int cell, int level ) {
	int i;
	float indicator = 0.0;
	int neighbors[num_neighbors];
	float drho[nDim];

	if ( refinement_indicator[DARK_MASS_INDICATOR].use[level] ) {
		indicator = MAX( dark_mass_indicator(cell, level), indicator );
	}

	if ( refinement_indicator[DARK_1STSPEC_INDICATOR].use[level] ) {
		indicator = MAX( dark_1stspec_indicator(cell, level), indicator );
	}

	if ( refinement_indicator[DARK_PARTICLE_COUNT_INDICATOR].use[level] ) {
		indicator = MAX( dark_particle_count_indicator( cell, level ), indicator );
	}

	if ( refinement_indicator[PLUGIN_INDICATOR].use[level] ) {
		indicator = MAX( plugin_indicator( cell, level ), indicator );
	}

#ifdef HYDRO
	if ( refinement_indicator[GAS_MASS_INDICATOR].use[level] ) {
		indicator = MAX( gas_mass_indicator( cell, level ), indicator );
	}

	if ( refinement_indicator[GAS_1STSPEC_INDICATOR].use[level] ) {
		indicator = MAX( gas_1stspec_indicator( cell, level ), indicator );
	}

	cell_all_neighbors( cell, neighbors );
	for ( i = 0; i < nDim; i++ ) {
		drho[i] = fabs( cell_gas_density( neighbors[2*i] ) 
			- cell_gas_density( neighbors[2*i+1] ) ) / 
			MIN(1.0e-35+cell_gas_density( neighbors[2*i] ), 
				1.0e-35+cell_gas_density( neighbors[2*i+1] ) );
	}

	if ( refinement_indicator[SHOCK_INDICATOR].use[level] ) {
		indicator = MAX( shock_indicator(cell, level,neighbors), indicator );
	}

	if ( refinement_indicator[CONTACT_DISCONTINUITY_INDICATOR].use[level] ) {
		indicator = MAX( contact_discontinuity_indicator(cell,level,neighbors,drho), indicator );
	}

	if ( refinement_indicator[DENSITY_GRADIENT_INDICATOR].use[level] ) {
		indicator = MAX( density_gradient_indicator(cell,level,neighbors,drho), indicator );
	}

	if ( refinement_indicator[PRESSURE_GRADIENT_INDICATOR].use[level] ) {
		indicator = MAX( pressure_gradient_indicator(cell,level,neighbors), indicator );
	}

	if ( refinement_indicator[ENTROPY_GRADIENT_INDICATOR].use[level] ) {
		indicator = MAX( entropy_gradient_indicator(cell,level,neighbors), indicator );
	}
        
	if ( refinement_indicator[SPATIAL_INDICATOR].use[level] ) {
		indicator = MAX( spatial_indicator(cell,level), indicator );
	}

	if ( refinement_indicator[JEANS_LENGTH_INDICATOR].use[level] ) {
		indicator = MAX( jeans_length_indicator(cell,level), indicator );
	}
#endif /* HYDRO */

	refinement_indicator(cell, 0) = indicator;
}

float dark_mass_indicator( int cell, int level ) {
	float ave_mass;

#ifdef PARTICLES
	ave_mass = (cell_first_species_mass(cell)) /
			refinement_indicator[DARK_MASS_INDICATOR].threshold[level];
#else
	ave_mass = 0.0;
#endif /* PARTICLES */

	return MIN( ave_mass, refinement_indicator[DARK_MASS_INDICATOR].weight );
}

float dark_1stspec_indicator( int cell, int level ) {
	float ave_mass;

#ifdef PARTICLES
	ave_mass = (cell_first_species_mass(cell)) /
	    (particle_species_mass[0] 
	     * refinement_indicator[DARK_1STSPEC_INDICATOR].threshold[level] );
#else
	ave_mass = 0.0;
	cart_error("need PARTICLES defined for gas_1stspec_indicator");
#endif /* PARTICLES */
	return MIN( ave_mass, refinement_indicator[DARK_1STSPEC_INDICATOR].weight );
}

#ifdef PARTICLES
int recursive_dark_particle_count( int cell, int num_highres_dark_particles ) {
	int i, ipart, count = 0;

	ipart = cell_particle_list[cell];
	while ( ipart != NULL_PARTICLE ) {
		if ( particle_id[ipart] < num_highres_dark_particles ) {
            count++;
        }
        ipart = particle_list_next[ipart];
    }

	if ( cell_is_refined(cell) ) {
		for ( i = 0; i < num_children; i++ ) {
			count += recursive_dark_particle_count( cell_child(cell,i), num_highres_dark_particles );
		}
	}

	return count;
}
#endif /* PARTICLES */

float dark_particle_count_indicator( int cell, int level ) {
#ifdef PARTICLES
	int num_highres_dark_particles;
	int count;

#ifdef STARFORM
	if ( num_particle_species > 2 ) {
		num_highres_dark_particles = particle_species_indices[num_particle_species-2];
	} else {
		num_highres_dark_particles = particle_species_indices[num_particle_species-1];
	}
#else
	if ( num_particle_species > 1 ) {
		num_highres_dark_particles = particle_species_indices[num_particle_species-1];
	} else {
		num_highres_dark_particles = particle_species_indices[num_particle_species];
	}	
#endif

	count = recursive_dark_particle_count(cell, num_highres_dark_particles);
/* Hard threshold works poorly with our diffusion scheme, better to use a ratio
	if ( count >= refinement_indicator[DARK_PARTICLE_COUNT_INDICATOR].threshold[level] ) {
		return 1.0;
	} 
*/
	return (float)count/refinement_indicator[DARK_PARTICLE_COUNT_INDICATOR].threshold[level]; 
#else 
	return 0.0;
#endif /* PARTICLES */
}

float plugin_indicator( int cell, int level ) {
        float indicator;
        PLUGIN_POINT(RefinementIndicator,(cell, level, &indicator));
        return indicator;
}

#ifdef HYDRO

float spatial_indicator( int cell, int level ) {
        /*
        // hack: refine around the central_cell -- initialized here to the center of the grid.
        // threshold[level] sets the code unit-distances to different levels of refinement 
        // the selected level gets fully refined (assuming weight==1) if r<distance(level) 
        */
	float in_region = 0;
	const double central_cell[nDim] = {num_grid/2,num_grid/2,num_grid/2};
	double pos[3];

	cell_center_position( cell, pos );
	if(refinement_indicator[SPATIAL_INDICATOR].threshold[level] > compute_distance_periodic(pos,central_cell)){
		in_region = refinement_indicator[SPATIAL_INDICATOR].weight;
	}
	return in_region ;
}

float gas_mass_indicator( int cell, int level ) {
	float ave_mass;

	ave_mass = ( cell_volume[level] * cell_gas_density(cell) ) / refinement_indicator[GAS_MASS_INDICATOR].threshold[level];
	return MIN( ave_mass, refinement_indicator[GAS_MASS_INDICATOR].weight );
}

float gas_1stspec_indicator( int cell, int level ) {
	float ave_mass;

#ifdef PARTICLES 
	ave_mass = ( cell_volume[level] * cell_gas_density(cell) ) / 
	    ( particle_species_mass[0] 
#ifdef COSMOLOGY
	      * cosmology->OmegaB/cosmology->OmegaM
#endif
	      * refinement_indicator[GAS_1STSPEC_INDICATOR].threshold[level] );
#else
	ave_mass = 0.0;
	cart_error("need particle defined for gas_1stspec_indicator");
#endif /* PARTICLES */
	return MIN( ave_mass, refinement_indicator[GAS_1STSPEC_INDICATOR].weight );
}

float shock_indicator( int cell, int level, int neighbors[] ) {
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
				MIN( cell_gas_pressure( neighbors[2*i] ),
					cell_gas_pressure( neighbors[2*i+1] ) );

			dv = fabs( cell_momentum(neighbors[2*i],i) ) / cell_gas_density( neighbors[2*i] )
				- cell_momentum(neighbors[2*i+1],i) / cell_gas_density( neighbors[2*i+1] );

			if ( dp > 3.0 && dv > 0.0 && 
					cell_gas_density(cell) > refinement_indicator[SHOCK_INDICATOR].threshold[level] ) {
				indicator = MAX( refinement_indicator[SHOCK_INDICATOR].weight, indicator );
			}
		}
	}

	return indicator;
}

float contact_discontinuity_indicator( int cell, int level, int neighbors[], float drho[] ) {
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
                	        MIN( cell_gas_pressure( neighbors[2*i] ),
                        	        cell_gas_pressure( neighbors[2*i+1] ) );
                                                                                             
	                if ( dp < 10.0 && drho[i] >= 0.6 ) {
        	                indicator = MAX( refinement_indicator[CONTACT_DISCONTINUITY_INDICATOR].weight, indicator );
	                }
		}
        }
                                                                                             
        return indicator;
}

float density_gradient_indicator( int cell, int level, int neighbors[], float drho[] ) {
	int i;
	float rho_ind1, rho_ind2;
	float drho_indicator;
	float indicator = 0.0;

	for ( i = 0; i < nDim; i++ ) {
		rho_ind1 = fabs( cell_gas_density(cell) - cell_gas_density(neighbors[2*i]) ) /
			MIN( cell_gas_density(cell), cell_gas_density(neighbors[2*i] ) );

		rho_ind2 = fabs(cell_gas_density(cell) - cell_gas_density(neighbors[2*i+1])) /
			MIN( cell_gas_density(cell), cell_gas_density(neighbors[2*i+1] ) );

		drho_indicator = MAX( drho[i], MAX( rho_ind1, rho_ind2 ) );
		if ( drho_indicator > 1.0 && 
			cell_gas_density(cell) > refinement_indicator[DENSITY_GRADIENT_INDICATOR].threshold[level] ) {
			indicator = MAX( indicator, refinement_indicator[DENSITY_GRADIENT_INDICATOR].weight*MAX(rho_ind1, rho_ind2) );
		}
	}

	return indicator;
}

float pressure_gradient_indicator( int cell, int level, int neighbors[] ) {
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

			if ( cell_gas_density(cell) > refinement_indicator[PRESSURE_GRADIENT_INDICATOR].threshold[level] ) {
				indicator = MAX( indicator, 
						refinement_indicator[PRESSURE_GRADIENT_INDICATOR].weight*MAX(press_ind1, press_ind2) );
			}
		}
	}

	return indicator;
}


extern double G_code;

float entropy_gradient_indicator( int cell, int level, int neighbors[] ) {
	return 0.0;
}

float jeans_length_indicator(int cell, int level)
{
  int i, j, nb[num_neighbors];
  float s, d, q, dv2, cs2, ljeans;
  float vc[nDim];

  if(cell_gas_density(cell) < 1.0e-30) return 0.0;

  /*
  //  Sound speed squared
  */
  cs2 = cell_gas_gamma(cell)*MAX((cell_gas_gamma(cell)-1.0)*cell_gas_internal_energy(cell),0.0)/cell_gas_density(cell);

  /*
  //  Local velocity dispersion squared
  */
  cell_all_neighbors(cell,nb);
  for(i=0; i<nDim; i++)
    {
      vc[i] = cell_var(cell,HVAR_MOMENTUM+i)/cell_gas_density(cell);
    }

  s = 0.0;
  for(j=0; j<num_neighbors; j++) if(cell_gas_density(nb[j]) > 0.0)
    {
      for(i=0; i<nDim; i++)
	{
	  d = cell_var(nb[j],HVAR_MOMENTUM+i)/cell_gas_density(nb[j]) - vc[i];
	  s += d*d;
	}
    }
  dv2 = s/num_neighbors;

  ljeans = sqrt(M_PI*(cs2+dv2)/(G_code*cell_gas_density(cell)));

  q = cell_size[level]*refinement_indicator[JEANS_LENGTH_INDICATOR].threshold[level]/ljeans;

  return MIN( q, refinement_indicator[JEANS_LENGTH_INDICATOR].weight );
}

#endif /* HYDRO */
#endif /* REFINEMENT */
