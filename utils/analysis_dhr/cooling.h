#ifndef __COOLING_H__
#define __COOLING_H__

#ifdef COOLING

void init_cooling();
void set_cooling_redshift( double a );
double cooling_rate( double rhogl, double T_g, double Z_met );
double ionization_fraction( double rhogl, double T_g, double Z_met );

#endif /* COOLING */

#endif
