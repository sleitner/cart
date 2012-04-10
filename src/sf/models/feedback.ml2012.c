#include "config.h"
#ifdef STARFORM

#include <math.h>
#include <string.h>

#include "auxiliary.h"
#include "control_parameter.h"
#include "cosmology.h"
#include "imf.h"
#include "particle.h"
#include "starformation.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "feedback.ml.h"


/*
//  Stellar Mass Loss, new form by Sam Leitner
*/

extern struct
{
  double loss_rate;       /* used to be called c0_ml */
  double time_interval;   /* used to be called T0_ml */
}
  ml;


void ml2012_config_init()
{
  ml.loss_rate = -1;
  ml.time_interval = -1;

  ml_config_init();
}


void ml2012_config_verify()
{
  /*
  //  if mass loss was not set in parameters then align it with the appropriate IMF.
  */
  /*    Leitner & Kravtsov 2011:  */
  if(ml.loss_rate == -1.0 ){
	  if(      strcmp("Salpeter",imf->name) == 0){
		  ml.loss_rate = 0.032 ;
	  }else if(strcmp("Miller-Scalo",imf->name) == 0){ //1979 
		  ml.loss_rate = 0.05 ; //left as old default -- 0.058 in paper 
	  }else if(strcmp("Chabrier",imf->name) == 0){ //Chabrier 2001~2003
		  ml.loss_rate = 0.046 ;
	  }else if(strcmp("Kroupa",imf->name) == 0){ //2001  
		  ml.loss_rate = 0.046 ;
	  }else{
		  cart_debug("IMF '%s' does not have associated IMF parameters.",imf->name);
		  cart_error("ART is terminating.");
	  }
  }

  if(ml.time_interval == -1.0 ){
	  if(      strcmp("Salpeter",imf->name) == 0){
		  ml.time_interval = 5.13e5 ;
	  }else if(strcmp("Miller-Scalo",imf->name) == 0){ //1979 
		  ml.time_interval = 5.0e6 ; //left as old default -- 6.04e6 in paper
	  }else if(strcmp("Chabrier",imf->name) == 0){ //Chabrier 2001~2003
		  ml.time_interval = 2.76e5 ;
	  }else if(strcmp("Kroupa",imf->name) == 0){ //2001  
		  ml.time_interval = 2.76e5 ;
	  }else{
		  cart_debug("IMF '%s' does not have associated IMF parameters.",imf->name);
		  cart_error("ART is terminating.");
	  }
  }

  ml_config_verify();
}


void ml2012_init()
{
  ml_init();
}


void ml2012_setup(int level)
{
  ml_setup(level);
}


#if defined(HYDRO) && defined(PARTICLES)
/*
//  Used with a #define for efficiency
//  void ml_hydrodynamic_feedback(int level, int cell, int ipart, double t_next )
*/
#endif /* HYDRO && PARTICLES */

#endif /* STARFORM */
