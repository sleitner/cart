#include "config.h"
#include "framework.h"


extern "C"
{
#include "hydro.h"
#include "rt.h"
#include "../extra/ism.h"
#include "../extra/utils.h"
#include "../drivers/start_analysis.h"
}


#include "celldata.h"


using namespace ng;


//
//  CellData class
//
const CellDataWorker& CellData::Worker(int i) const
{
  return this->Data()[i];
}


bool CellData::IsValid(const CellDataWorker& w) const
{
  return (w.Value != NULL);  // C-style NULL
}


#define DECL(var,header,weight) \
  float var##_fun(int level, int cell, double *ref_pos, double *ref_vel); \
  CellDataWorker var##_impl = { var##_fun, header, weight }; \
  const CellDataWorker& var(var##_impl); \
  float var##_fun(int level, int cell, double *ref_pos, double *ref_vel)


namespace ng
{
  bool operator==(const CellDataWorker& w1, const CellDataWorker& w2)
  {
    return w1.Value==w2.Value && w1.Header==w2.Header && w1.WeightId==w2.WeightId;
  }

  namespace cd
  {
    CellDataWorker null = { NULL, NULL, 0 };  // Need C-style NULL, not 0

#ifdef GRAVITY

    DECL( total_density, "total density (g/cm^3)", 0 )
    {
      return units->density*cell_total_density(cell);
    }

#else

    const CellDataWorker& total_density(null);

#endif // GRAVITY


#if defined(GRAVITY) && defined(STARFORM)

    DECL( stellar_density, "stellar density (g/cm^3)", 0 )
    {
      return units->density*cell_stellar_density(cell);
    }

#else

    const CellDataWorker& stellar_density(null);

#endif // GRAVITY && STARFORM


#ifdef HYDRO

    DECL( gas_density, "gas density (g/cm^3)", 0 )
    {
      return units->density*cell_gas_density(cell);
    }

    DECL( baryon_number_density, "baryon number density (cm^{-3})", 0 )
    {
      return units->number_density*cell_gas_density(cell);
    }

    DECL( temperature, "temperature (K)", 2 )
    {
      return units->temperature*cell_gas_temperature(cell);
    }

    DECL( baryon_column_density, "baryon column density (cm^{-2})", 0 )
    {
      return cell_gas_density(cell)*cell_sobolev_length2(cell,level,NULL)*units->number_density*units->length;
    }

    DECL( thermal_pressure, "thermal pressure (erg/cm^{-3})", 0 )
    {
      return units->energy_density*(cell_gas_gamma(cell)-1)*cell_gas_internal_energy(cell);
    }

    DECL( rotational_velocity, "rotational velocity (km/s)", 3 )
    {
      int j, cc, vars[4];
      double vr, r2, v2;
      double dr[nDim];
      float dv[nDim];
      float H0 = 100.0*constants->kms/constants->Mpc*units->time;
      float vals[4];

      cart_assert(ref_pos != NULL);

      for(j=0; j<nDim; j++) vars[j] = HVAR_MOMENTUM + j;
      vars[nDim] = HVAR_GAS_DENSITY;

      cc = cell_find_position(ref_pos);
      cell_interpolate_at_position(cc,ref_pos,nDim+1,vars,vals);
      for(j=0; j<nDim; j++)
	{
	  dv[j] = vals[j]/vals[nDim];
      }

      cell_center_position(cell,dr);
      compute_displacement_periodic(ref_pos,dr,dr);
      for(j=0; j<nDim; j++)
	{
	  dv[j] = cell_momentum(cell,j)/cell_gas_density(cell) - dv[j] + H0*dr[j];
	}

      for(vr=r2=v2=0.0, j=0; j<nDim; j++)
	{
	  r2 += dr[j]*dr[j];
	  vr += dr[j]*dv[j];
	  v2 += dv[j]*dv[j];
	}

      return sqrt(v2-vr*vr/(1.0e-10+r2))*units->velocity/constants->kms;
    }

    DECL( local_gas_overdensity, "gas density in this cell over the average density in neighboring cells", 2 )
    {
      float s;
      int i, nb[num_neighbors];

      cell_all_neighbors(cell,nb);
      //
      //  Mean of all neighbors
      //
      s = 0.0;
      for(i=0; i<num_neighbors; i++)
	{
	  s += cell_gas_density(nb[i]);
	}
      s /= num_neighbors;
      return cell_gas_density(cell)/s;
    }

#else

    const CellDataWorker& gas_density(null);
    const CellDataWorker& baryon_number_density(null);
    const CellDataWorker& temperature(null);
    const CellDataWorker& baryon_column_density(null);
    const CellDataWorker& thermal_pressure(null);
    const CellDataWorker& radiation_pressure(null);
    const CellDataWorker& rotational_velocity(null);
    const CellDataWorker& local_gas_overdensity(null);

#endif // HYDRO

#if defined(HYDRO) && defined(ENRICH)

    DECL( gas_metallicity, "gas metallicity (solar units)", 2 )
    {
      return cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
    }

#else

    const CellDataWorker& gas_metallicity(null);

#endif // HYDRO && ENRICH


#if defined(HYDRO) && defined(RADIATIVE_TRANSFER)

    DECL( HI_fraction, "HI fraction", 2 )
    {
      return cell_HI_fraction(cell);
    }

    DECL( HII_fraction, "HII fraction", 2 )
    {
      return cell_HII_fraction(cell);
    }

    DECL( HeI_fraction, "HeI fraction", 2 )
    {
      return cell_HeI_fraction(cell);
    }

    DECL( HeII_fraction, "HeII fraction", 2 )
    {
      return cell_HeII_fraction(cell);
    }

    DECL( HeIII_fraction, "HeIII fraction", 2 )
    {
      return cell_HeIII_fraction(cell);
    }

    DECL( H2_fraction, "H2 fraction", 2 )
    {
      return cell_H2_fraction(cell);
    }

    DECL( HI_number_density, "HI number density", 0 )
    {
      return units->number_density*cell_HI_density(cell);
    }

    DECL( HII_number_density, "HII number density", 0 )
    {
      return units->number_density*cell_HII_density(cell);
    }

    DECL( HeI_number_density, "HeI number density", 0 )
    {
      return units->number_density*cell_HeI_density(cell);
    }

    DECL( HeII_number_density, "HeII number density", 0 )
    {
      return units->number_density*cell_HeII_density(cell);
    }

    DECL( HeIII_number_density, "HeIII number density", 0 )
    {
      return units->number_density*cell_HeIII_density(cell);
    }

    DECL( H2_number_density, "H2 number density", 0 )
    {
      return units->number_density*cell_H2_density(cell);
    }

    DECL( dust_to_gas_ratio, "dust-to-gas ratio in MW units (D_{MW})", 2 )
    {
      return rtDmw(cell);
    }

    DECL( interstellar_radiation_field, "radiation field in MW units (U_{MW})", 0 )
    {
      return rtUmw(cell);
    }

    DECL( cooling_rate, "cooling rate per baryon in ergs/s", 2 )
    {
      float cfun, hfun;
      rtGetCoolingRate(cell,&cfun,&hfun);
      return cfun;
    }

    DECL( heating_rate, "heating rate per baryon in ergs/s", 2 )
    {
      float cfun, hfun;
      rtGetCoolingRate(cell,&cfun,&hfun);
      return hfun;
    }

#else

    const CellDataWorker& HI_fraction(null);
    const CellDataWorker& HII_fraction(null);
    const CellDataWorker& HeI_fraction(null);
    const CellDataWorker& HeII_fraction(null);
    const CellDataWorker& HeIII_fraction(null);
    const CellDataWorker& H2_fraction(null);
    const CellDataWorker& HI_number_density(null);
    const CellDataWorker& HII_number_density(null);
    const CellDataWorker& HeI_number_density(null);
    const CellDataWorker& HeII_number_density(null);
    const CellDataWorker& HeIII_number_density(null);
    const CellDataWorker& H2_number_density(null);
    const CellDataWorker& dust_to_gas_ratio(null);
    const CellDataWorker& interstellar_radiation_field(null);
    const CellDataWorker& cooling_rate(null);
    const CellDataWorker& heating_rate(null);

#endif // HYDRO && RADIATIVE_TRANSFER


  };
};
