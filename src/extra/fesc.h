#ifndef __EXT_FESC_H__
#define __EXT_FESC_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO_LIST;


/*
//  Information for what is actually escaping - column densities for
//  <num> grid variables specified in <vars[]> are computed and
//  combined together into a single escape fraction in the collector function.
*/
struct feInfo
{
  int num;
  int *var_ids;
  double (*collector)(int num, double *vals);
};

extern struct feInfo fe_912A;


/*
//  Compute escape fractions specified by feInfo<a> using the HealPIX 
//  binning with <nside> bins per section (total number of rays is 
//  12*nside^2) for all halos resolved at least to <resolution_level>;
//  escape fractions are integrated to <outer_edge> fraction of the virial
//  radius; sources can be downsamples by <sample> to spped-up the 
//  calculation; sources are stellar particles or the grid-assigned radiation
//  source (usually faster) dependending on <sum_over_parts>.
*/
void extEscapeFraction(const char *fname, struct feInfo a, int nside, struct HALO_LIST *halos, int resolution_level, float outer_edge, int sample, int sum_over_parts);

#endif  /* __EXT_FESC_H__ */
