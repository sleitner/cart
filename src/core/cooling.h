#ifndef __COOLING_H__
#define __COOLING_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef COOLING

void init_cooling();
void set_cooling_redshift( double a );

#ifdef OLDSTYLE_COOLING_EXPLICIT_SOLVER
typedef double cooling_t;
#else
typedef struct COOLING_TYPE
{
  double Cooling;
  double Heating;
}
cooling_t;
#endif

cooling_t cooling_rate( double nHlog, double T_g, double Zlog );
double cooling_fion( double rhogl, double T_g, double Zlog );

#endif /* COOLING */

#endif
