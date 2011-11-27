#include "config.h"


#include <math.h>
#include <stdio.h>
#include <string.h>


#include "auxiliary.h"
#include "cosmology.h"
#include "hydro.h"
#include "rt.h"
#include "tree.h"
#include "units.h"

#include "extra/start_analysis.h"

#include "cell_dump.h"


#ifdef RADIATIVE_TRANSFER
const char *ngCDH_p_[] = {
  "baryon number density (cm^{-3})",
  "temperature (K)",
  "baryon column density (cm^{-2})",
  "gas metallicity (solar units)",
  "HI fraction",
  "HII fraction",
  "H_2 fraction",
  "radiation field in MW units (U_{MW})",
  "dust-to-gas ratio in MW units (D_{MW}",
  "density in this cell over the average density in neighboring cells",
  "cooling rate (erg/s)",
  "heating rate (erg/s)"
};
#else
const char *ngCDH_p_[] = {
  "baryon number density (cm^{-3})",
  "temperature (K)",
  "baryon column density (cm^{-2})",
  "gas metallicity (solar units)",
  "total density (g/cm^3)",
};
#endif


#ifdef RADIATIVE_TRANSFER
const int ngCDW_p_[] = {
  0,
  2,
  2,
  2,
  2,
  2,
  2,
  0,
  2,
  0,
  2,
  2
};
#else
const int ngCDW_p_[] = {
  0,
  2,
  2,
  2,
  0
};
#endif


void ngCDF_p_(int level, int cell, int num, float *ptr, double *halo_pos, float *halo_vel)
{
  float s, cr, hr;
  int i, nb[num_neighbors];

#ifdef RADIATIVE_TRANSFER
  cart_assert(num >= 12);
#else
  cart_assert(num >= 5);
#endif

#ifdef HYDRO
  ptr[0] = units->number_density*cell_gas_density(cell);
  ptr[1] = units->temperature*cell_gas_temperature(cell);

  ptr[2] = cell_gas_density(cell)*cell_sobolev_length2(cell,level,NULL)*units->number_density*units->length;
#else
  ptr[0] = ptr[1] = ptr[2] = 0.0;
#endif /* HYDRO */

#ifdef ENRICH
  ptr[3] = cell_gas_metal_density(cell)/(constants->Zsun*cell_gas_density(cell));
#else
  ptr[3] = 0.0;
#endif /* ENRICH */

#ifdef RADIATIVE_TRANSFER
  ptr[4] = cell_HI_fraction(cell);
  ptr[5] = cell_HII_fraction(cell);
  ptr[6] = cell_H2_fraction(cell);
  ptr[7] = rtUmw(cell);
  ptr[8] = rtDmw(cell);

  cell_all_neighbors(cell,nb);
  /*
  //  Mean of all neighbors
  */
  s = 0.0;
  for(i=0; i<num_neighbors; i++)
    {
      s += cell_gas_density(nb[i]);
    }
  s /= num_neighbors;
  ptr[9] = cell_gas_density(cell)/s;

  rtGetCoolingRate(cell,&cr,&hr);
  ptr[10] = cr;
  ptr[11] = hr;

#else  /* RADIATIVE_TRANSFER */

  ptr[4] = units->density*cell_total_density(cell)/(constants->g/pow(constants->cm,3.0));

#endif /* RADIATIVE_TRANSFER */
}


struct ngCellDump ngCD_p_ ;//= { sizeof(ngCDW_p_)/sizeof(int), ngCDF_p_, ngCDW_p_, ngCDH_p_ };
