#include "config.h"
#include "framework.h"


extern "C"
{
#include "../extra/ifrit.h"
#include "../extra/igm.h"
#include "../extra/ism.h"
#include "../drivers/start_analysis.h"
}


#include "current.h"

#include "celldata.h"
#include "ifrit.h"


namespace ng
{
  namespace h2
  {
    const CellData& cd()
    {
      static CellData d;

      if(d.Size() == 0)
	{
	  d += cd::baryon_number_density;
	  d += cd::temperature;
	  d += cd::baryon_column_density;
	  d += cd::gas_metallicity;
	  d += cd::HI_fraction;
	  d += cd::HII_fraction;
	  d += cd::H2_fraction;
	  d += cd::interstellar_radiation_field;
	  d += cd::dust_to_gas_ratio;
	  d += cd::local_gas_overdensity;
	}

      return d;
    }

    const ifrit::VarSet& vs()
    {
      static ifrit::VarSet d;

      if(d.Size() == 0)
	{
#ifdef RADIATIVE_TRANSFER
	  d += I_HI_FRACTION;
	  d += I_GAS_NUMBER_DENSITY;
	  d += I_GAS_TEMPERATURE;
	  d += I_FRACTION+RT_HVAR_OFFSET+5;
	  d += RT_HVAR_OFFSET+5;
	  d += I_CELL_LEVEL;
	  d += I_LOCAL_PROC;
#endif /* RADIATIVE_TRANSFER */
	}

      return d;
    }
  };

  namespace btfr
  {
    const CellData& cd()
    {
      static CellData d;

      if(d.Size() == 0)
	{
	  d += cd::HI_fraction;
	  d += cd::HII_fraction;
	  d += cd::H2_fraction;
	  d += cd::total_density;
	  d += cd::gas_density;
	  d += cd::stellar_density;
	  d += cd::rotational_velocity;
	}

      return d;
    }
  };

  namespace rei
  {
    const CellData& cd()
    {
      static CellData d;

      if(d.Size() == 0)
	{
	  d += cd::baryon_number_density;
	  d += cd::temperature;
	  d += cd::baryon_column_density;
	  d += cd::gas_metallicity;
	  d += cd::HI_fraction;
	  d += cd::HII_fraction;
	  d += cd::H2_fraction;
	  d += cd::interstellar_radiation_field;
	  d += cd::dust_to_gas_ratio;
	  d += cd::total_density;
	  d += cd::stellar_density;
	}

      return d;
    }

    const ifrit::VarSet& vs()
    {
      static ifrit::VarSet d;

      if(d.Size() == 0)
	{
#ifdef RADIATIVE_TRANSFER
	  d += I_HI_FRACTION;
	  d += I_GAS_NUMBER_DENSITY;
	  d += I_GAS_TEMPERATURE;
	  d += I_FRACTION+RT_HVAR_OFFSET+0;
	  d += RT_HVAR_OFFSET+0;
	  d += I_CELL_LEVEL;
	  d += I_LOCAL_PROC;
#endif /* RADIATIVE_TRANSFER */
	}

      return d;
    }
  };
};

using namespace ng;


//
//  Class MassFractions
//
//MassFractions::MassFractions(){}


void MassFractions::Exec(const char *path)
{
#if defined(COSMOLOGY) && defined(HYDRO)
  init->TotalDensity();

  extMassFractions(this->File(path),this->Halos());
#else
  cart_debug("COSMOLOGY and HYDRO are not set. Skipping writing mass fractions...");
#endif /* COSMOLOGY && HYDRO */
}


//
//  Class KSR
//
KSRelation::KSRelation() : dt(p.TimeScale), dl(p.LengthScale), lim(p.StellarAgeLimit)
{
  p.TimeScale = 20.0;
  p.LengthScale = 500.0;
  p.StellarAgeLimit = 1.0e9;
}


KSRelation::KSRelation(float length, float time, float limit) : dt(p.TimeScale), dl(p.LengthScale), lim(p.StellarAgeLimit)
{
  p.TimeScale = time;
  p.LengthScale = length;
  p.StellarAgeLimit = limit;
}


void KSRelation::SetTimeScale(float time)
{
  cart_assert(time>=0.0 && time<1.0e6);

  p.TimeScale = time;
}


void KSRelation::SetLengthScale(float length)
{
  cart_assert(length>0.0 && length<1.0e6);

  p.LengthScale = length;
}


void KSRelation::SetStellarAgeLimit(float limit)
{
  cart_assert(limit > 0.0);

  p.StellarAgeLimit = limit;
}


void KSRelation::Exec(const char *path)
{
#if defined(HYDRO) && defined(STARFORM)
  init->RadiativeTransfer();

  if(p.TimeScale > 0.0)
    {
#if defined(PARTICLES)
      extStarFormationLaw(this->File(path),p.LengthScale,p.TimeScale,p.StellarAgeLimit,this->Halos());
#else
      cart_debug("PARTICLES are not set. Skipping writing Kennicutt-Schmidt relation for non-zero time-scale...");
#endif
    }
  else
    {
      extStarFormationLaw2(this->File(path),p.LengthScale,this->Halos());
    }
#else
  cart_debug("HYDRO and STARFORM are not set. Skipping writing Kennicutt-Schmidt relation...");
#endif
}

