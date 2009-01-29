#ifndef __EXT_IFRIT_H__
#define __EXT_IFRIT_H__


/*
// Additional options for cell variables; used in binning on a uniform mesh.
*/
#define EXT_FRACTION            1000
#define EXT_GAS_TEMPERATURE     2000
#define EXT_CELL_LEVEL          2001
#define EXT_LOCAL_PROC          2002

int extWriteIfritFile(int level, int *nbinIn, double *bbIn, int nvars, int *varid, const char *filename);


#endif  /* __EXT_IFRIT_H__ */
