#ifndef __COSMOLOGY_H__
#define __COSMOLOGY_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef COSMOLOGY

/*
//  INTERNAL DECLARATIONS
*/

struct CosmologyParameters
{
  double OmegaM;
  double OmegaD;
  double OmegaB;
  double OmegaL;
  double OmegaK;
  double OmegaR;
  double h;
  double DeltaDC;
  int flat;
  double Omh2;
  double Obh2;
};
extern const struct CosmologyParameters *cosmology;


#define COSMOLOGY_DECLARE_PRIMARY_PARAMETER(name) \
void cosmology_set_##name(double value)

#define cosmology_set(name,value)	\
cosmology_set_##name(value)

COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaM);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaB);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaL);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(h);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(DeltaDC);


/*
//  USER INTEFACE:
//
//  To access parameters, use cosmology-><parameter>, for example
//    cosmology->OmegaM for the matter density parameter.
//
//  To set the parameter, call cosmology_set(<parameter>,<value>), for example
//    cosmology_set(OmegaM,0.3);
*/

/*
//  Copy constructor
*/
void cosmology_copy(const struct CosmologyParameters *c);


/*
//  Check that all required cosmological parameters have been set.
//  The minimum set is OmegaM, OmegaB, and h. By default, zero OmegaL,
//  OmegaK, and the DC mode are assumed.
*/
int cosmology_is_set();


/*
//  Freeze the cosmology and forbid any further changes to it.
//  In codes that include user-customizable segments (like plugins),
//  this function van be used for insuring that a user does not
//  change the cosmology in mid-run.
*/
void cosmology_set_fixed();


/*
//  Manual initialization. This does not need to be called, 
//  the initialization is done automatically on the first call
//  to a relevant function.
*/
void cosmology_init();


/*
//  Set the range of global scale factors for thread-safe
//  calls to direct functions until the argument leaves the range.
*/
void cosmology_set_thread_safe_range(double amin, double amax);

/*
//  Direct functions take the global cosmological scale factor as the argument.
//  These functionsare are thread-safe if called with the argument in the
//  range set by a prior call to cosmology_set_thread_safe_range(...).
//  Calling them with the argument outside that range is ok, but breaks
//  thread-safety assurance.
*/
double  aBox(double a);
double tCode(double a);
double tPhys(double a);
double dPlus(double a);
double qPlus(double a); /* Q+ = a^2 dD+/(H0 dt) */

/*
//  Inverse conversions (converting other variables to the global 
//  scale factor aUni). These functions are NOT generally thread-safe.
*/
double inv_aBox(double abox);
double inv_tCode(double tcode);
double inv_tPhys(double tphys);
double inv_dPlus(double dplus);

/*
//  Conversion macros
*/
#define  abox_from_auni(a)   aBox(a)
#define tcode_from_auni(a)  tCode(a)
#define tphys_from_auni(a)  tPhys(a)
#define dplus_from_auni(a)  dPlus(a)

#define auni_from_abox(v)   inv_aBox(v)
#define auni_from_tcode(v)  inv_tCode(v)
#define auni_from_tphys(v)  inv_tPhys(v)
#define auni_from_dplus(v)  inv_dPlus(v)

#define abox_from_tcode(tcode)   aBox(inv_tCode(tcode))
#define tcode_from_abox(abox)    tCode(inv_aBox(abox))

#define tphys_from_abox(abox)    tPhys(inv_aBox(abox))
#define tphys_from_tcode(tcode)  tPhys(inv_tCode(tcode))
#define dplus_from_tcode(tcode)  dPlus(inv_tCode(tcode))

#endif /* COSMOLOGY */

#endif /* __COSMOLOGY_H__ */

