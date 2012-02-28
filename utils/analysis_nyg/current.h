#ifndef __NG_CURRENT_H__
#define __NG_CURRENT_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include "profiles.h"


namespace ng
{
  class CellData;

  namespace ifrit
  {
    class VarSet;
  };

  namespace h2
  {
    const CellData& cd();
    const ifrit::VarSet& vs();
  };

  namespace btfr
  {
    const CellData& cd();
  };

  namespace rei
  {
    const CellData& cd();
    const ifrit::VarSet& vs();
  };


  //
  //  Compute gas and other mass fractions for halos.
  //
  class MassFractions : public Algorithm
  {
  public:

    //MassFractions();

    virtual void Exec(const char *path);
  };

  //
  //  Write the Kennicutt-Schmidth relation on scale <LengthScale> kpc,
  //  averaged over <TimeScale> Myr, and limit stellar age in output M*
  //  to <StellarAgeLimit> Myr. If <TimeScale> is zero, the instanteneous KSR
  //  is returned.
  //
  class KSRelation : public Algorithm
  {
  public:

    KSRelation();
    KSRelation(float length, float time, float limit = 1.0e9);

    void SetTimeScale(float time);
    void SetLengthScale(float length);
    void SetStellarAgeLimit(float limit);

    virtual void Exec(const char *path);

    Shorthand<float> dt;
    Shorthand<float> dl;
    Shorthand<float> lim;

  protected:

    struct Pars
    {
      float TimeScale;
      float LengthScale;
      float StellarAgeLimit;
    }
    p;
  };
};


#endif // __NG_CURRENT_H__
