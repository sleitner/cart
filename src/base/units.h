#ifndef __UNITS_H__
#define __UNITS_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct Constants
{
  double K;
  double cms;
  double erg;
  double barye;
  double dyne;
  double gpercc;
  double cc;
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
  double Dsun;
  double Lsun;
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
#ifdef COSMOLOGY
  double length_in_chimps;
#endif /* COSMOLOGY */
};
extern const struct Units *units;


extern double box_size;


#ifndef COSMOLOGY
void units_set(double mass, double time, double length);
#endif /* COSMOLOGY */


void units_init();
void units_update(int level);

#endif
