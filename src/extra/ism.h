#ifndef __EXT_ISM_H__
#define __EXT_ISM_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO_LIST;

typedef void (*DumpWorker)(int level, int cell, int num, float *ptr);


/*
//  Set cell_var(c,var) with the halo index+1 for each halo, or 0 if belongs to none;
//  a cell belongs to a halo if it is inside its size_factor*Rtrunc, and satellites are
//  accounted for properly.
*/
void extMapHaloCells(int var, int floor_level, const struct HALO_LIST *halos, float size_factor);

void extDumpLevels(const char *fname, DumpWorker worker, int level1, int level2, const struct HALO_LIST *halos);

void extDumpProfiles(const char *fname, DumpWorker worker, int floor_level, float rmin, float rmax, int ndex, const struct HALO_LIST *halos);

/*
//  SF law from stellar particles
*/
void extStarFormationLaw(const char *fname, float spatial_scale, float time_scale, float stellar_age_limit, const struct HALO_LIST *halos);

/*
//  SF law for instantaneous star formation from gas only
*/
void extStarFormationLaw2(const char *fname, float spatial_scale, const struct HALO_LIST *halos);

/*
//  Compute column densities and compare them with the Sobolev-like approximations
*/
void extCheckSobolevApproximations(const char *fname, int floor_level, int nside, double len);

#endif  /* __EXT_ISM_H__ */
