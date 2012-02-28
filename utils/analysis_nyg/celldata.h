#ifndef __NG_CELLDATA_H__
#define __NG_CELLDATA_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


extern "C"
{
  struct DUMP_WORKER;
}


#include "list.h"


namespace ng
{
  typedef struct DUMP_WORKER CellDataWorker;

  class CellData : public List<CellDataWorker>
  {
  public:

    const CellDataWorker& Worker(int i) const;

  protected:

    virtual bool IsValid(const CellDataWorker& w) const;
  };

  bool operator==(const CellDataWorker&,const CellDataWorker&);

  namespace cd
  {
    extern const CellDataWorker& total_density;
    extern const CellDataWorker& stellar_density;
    extern const CellDataWorker& gas_density;
    extern const CellDataWorker& baryon_number_density;
    extern const CellDataWorker& temperature;
    extern const CellDataWorker& baryon_column_density;
    extern const CellDataWorker& thermal_pressure;
    extern const CellDataWorker& radiation_pressure;
    extern const CellDataWorker& rotational_velocity;
    extern const CellDataWorker& local_gas_overdensity;
    extern const CellDataWorker& gas_metallicity;
    extern const CellDataWorker& HI_fraction;
    extern const CellDataWorker& HII_fraction;
    extern const CellDataWorker& HeI_fraction;
    extern const CellDataWorker& HeII_fraction;
    extern const CellDataWorker& HeIII_fraction;
    extern const CellDataWorker& H2_fraction;
    extern const CellDataWorker& HI_number_density;
    extern const CellDataWorker& HII_number_density;
    extern const CellDataWorker& HeI_number_density;
    extern const CellDataWorker& HeII_number_density;
    extern const CellDataWorker& HeIII_number_density;
    extern const CellDataWorker& H2_number_density;
    extern const CellDataWorker& dust_to_gas_ratio;
    extern const CellDataWorker& interstellar_radiation_field;
    extern const CellDataWorker& cooling_rate;
    extern const CellDataWorker& heating_rate;
  };
};


#define NEW_CELLDATA_WORKER(var,header,weight) \
  float var##_fun(int level, int cell, double *ref_pos, float *ref_vel); \
  CellDataWorker var = { var##_fun, header, weight }; \
  float var##_fun(int level, int cell, double *ref_pos, float *ref_vel)

#endif  // __NG_CELLDATA_H__
