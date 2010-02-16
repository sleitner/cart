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

#define I_FLAG_STARS          1
#define I_FLAG_SPLIT_STARS    2
#define I_FLAG_ATTR_IS_MASS   4


struct IFRIT_NAMESPACE
{
  int (*OutputMesh)(const char *filename, int level, int *nbinIn, double *bbIn, int nvars, int *varid);
  void (*OutputHalo)(const char *fname, int floor_level, float zoom, const struct HALO *h, int nvars, int *varid);
};

extern const struct IFRIT_NAMESPACE ifrit;

#endif  /* __EXT_IFRIT_H__ */
