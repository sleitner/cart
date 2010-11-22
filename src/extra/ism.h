#ifndef __EXT_ISM_H__
#define __EXT_ISM_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO_LIST;

typedef void (*DumpWorker)(int level, int cell, int num, float *ptr);


/*
//  Dump ASCII files with cell data and spherically averaged profiles.
//  If halos is not NULL, these two calls overwrite VAR_ACCEL variable.
//  <weight_id> array determined as the profiles are weighted: 
//    0: volume average
//    1: total mass average
//    2: baryonic mass average
*/
void extDumpLevels(const char *fname, int nout, DumpWorker worker, const char **header, int level1, int level2, struct HALO_LIST *halos);

void extDumpProfiles(const char *fname, int nout, DumpWorker worker, const char **header, const int *weight_id, int resolution_level, float rmin, float rmax, int ndex, struct HALO_LIST *halos);

/*
//  SF law from stellar particles
*/
#if defined(PARTICLES) && defined(STARFORM)
void extStarFormationLaw(const char *fname, float spatial_scale, float time_scale, float stellar_age_limit, const struct HALO_LIST *halos);
#endif

/*
//  SF law for instantaneous star formation from gas only
*/
#if defined (HYDRO) && defined(STARFORM)
void extStarFormationLaw2(const char *fname, float spatial_scale, const struct HALO_LIST *halos);
#endif

/*
//  Dump stellar particles for a given halo within a given fraction of Rvir
*/
void extHaloStars(const char *fname, const struct HALO *h, float rmax);

#endif  /* __EXT_ISM_H__ */
