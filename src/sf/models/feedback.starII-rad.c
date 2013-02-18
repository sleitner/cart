#include "config.h"
#ifdef STARFORM
#ifdef STAR_PARTICLE_TYPES

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "particle.h"
#include "starformation.h"
#include "units.h"

#include "onestarfits.h"
#include "feedback.rad.h"
#include "feedback.starII-rad.h"

/* provides the luminosity used in density calculation*/

extern double tdelay_popM_feedback;
/* RT source function for individual stars */
float rad_luminosity_popM_starII0(int ipart){
    float Lion ; 
    double Zsol;
    if(!particle_is_star(ipart)) return 0.0;
    
    if(star_particle_type[ipart] == STAR_TYPE_STARII){
	Lion = (float) OneStar_Lion_RT(ipart);  
	return Lion;
    }else{
	if( particle_t[ipart] - star_tbirth[ipart] > tdelay_popM_feedback ){  
	    Lion = rad_luminosity_popM(ipart);
	    return Lion;
	}else{
	/* should return "\Int Lion{0,min_starII_mass}" at early times, but min_starII_mass should be small enough anyway*/
	    return 0; 
	}
    }
}




#endif /* STAR_PARTICLE_TYPES */
#endif /* STARFORM */
