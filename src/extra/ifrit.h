#ifndef __EXT_IFRIT_H__
#define __EXT_IFRIT_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO;


/*
// Additional options for cell variables; used in binning on a uniform mesh.
*/
#define I_FRACTION            1000
#define I_GAS_TEMPERATURE     2000
#define I_CELL_LEVEL          2001
#define I_LOCAL_PROC          2002
#define I_GAS_NUMBER_DENSITY  2003
#define I_GAS_METAL_DENSITY   2004
#define I_GAS_METALLICITY     2005
#define I_GAS_TOVERMU         2006
#define I_GAS_OVERDENSITY     2007
#define I_HI_FRACTION         2008
#define I_H2_FRACTION         2009
#define I_DMW                 2010
#define I_UMW                 2011

#define I_FLAG_STARS          1
#define I_FLAG_SPLIT_STARS    2
#define I_FLAG_ATTR_IS_MASS   4


struct IFRIT_NAMESPACE
{
  int  (*OutputMesh)(const char *filename, int floor_level, int nbin[], const double pcen[], int nvars, const int *varid);
  void (*OutputHalo)(const char *filename, int floor_level, int nbin[], const struct HALO *h, int nvars, const int *varid);
  void (*OutputBox)(const char *fileroot, int floor_level, int nbin[], const double pos[], int nvars, const int *varid);
};

extern const struct IFRIT_NAMESPACE ifrit;

#endif  /* __EXT_IFRIT_H__ */
