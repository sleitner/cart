#ifndef __UNITS_H__
#define __UNITS_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct Constants
{
  double cm;
  double g;
  double s;
  double K;
  double yr;
  double Myr;
  double Gyr;
  double pc;
  double kpc;
  double Mpc;
  double kms;
  double mp;
  double k;
  double G;
  double c;
  double eV;
  double amu;
  double mH;
  double mHe;
  double Msun;
  double Zsun;
  double Yp;
  double wmu;
  double wmu_e;
  double XH;
  double XHe;
  double gamma;
  double sigmaT;
};
extern const struct Constants *constants;


struct PrimaryUnits
{
  int set;
  double mass;
  double time;
  double length;
};
extern const struct PrimaryUnits *primary_units;


struct Units
{
  double mass;
  double time;
  double length;
  double energy;
  double density;
  double velocity;
  double potential;
  double temperature;
  double energy_density;
  double number_density;
  double Emin;
#ifdef COSMOLOGY
  double length_in_chimps;
#endif /* COSMOLOGY */
};
extern const struct Units *units;


#ifdef CHECK_LEGACY_UNITS
struct LegacyUnits {
  double H0;
  double r0;
  double t0;
  double v0;
  double rho0;
  double den0;
  double P0;
  double T0;
  double E0;
  double M0;
  double S0;
  double AL_SD;
};
extern const struct LegacyUnits *legacy_units;
#endif /* CHECK_LEGACY_UNITS */


extern double box_size;


#ifndef COSMOLOGY
void units_set(double mass, double time, double length);
#endif /* COSMOLOGY */

void config_init_units();
void config_verify_units();

void units_reset();
void units_update(int level);

#endif
