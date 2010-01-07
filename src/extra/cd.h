#ifndef __EXT_CD_H__
#define __EXT_CD_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


/*
//  Buffer for keeping information used during the LOS traversal
*/
#define CD_MAX_NVARS 10
typedef struct cdDataType
{
  double val[CD_MAX_NVARS];
  double len;
}
cdData;


/*
//  Compute column densities of <nvars> grid variables (with indicies 
//  supplied in <vars[]> array, sampled over the sky using the HealPIX 
//  binning with <nside> bins per section (total number of rays is 
//  12*nside^2) with common origin <pos>, of maximum length <len>, not deeper 
//  than <floor_level>. The final result is returned in <output[]>, whose
//  dimension must be 12*nside^2 and whose val[] arrays must be allocated to 
//  hold <nvars> values.
*/
void cdTraverseSky(int nvars, int vars[], int nside, double pos[3], double len, int floor_level, cdData output[]);

#endif  /* __EXT_CD_H__ */
