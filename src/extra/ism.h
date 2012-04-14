#ifndef __EXT_ISM_H__
#define __EXT_ISM_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO;
struct HALO_LIST;


/*
//  Dump ASCII files with cell data and spherically averaged profiles.
//  <WeightId> values determine how the profiles are weighted: 
//    0: volume average
//    1: total mass average
//    2: baryonic mass average
//    3: HI mass average
*/
typedef struct DUMP_WORKER
{
  float (*Value)(int level, int cell, double *ref_pos, float *ref_vel);
  const char *Header;
  int WeightId;
}
DumpWorker;

void extDumpLevels(const char *fname, int nout, const DumpWorker *workers, int level1, int level2, struct HALO_LIST *halos);

void extDumpLevelsLowMemory(const char *fname, int nout, const DumpWorker *workers, int level1, int level2, struct HALO_LIST *halos);

#ifdef COSMOLOGY
void extDumpHaloProfiles(const char *fname, int nout, const DumpWorker *workers, float rmin, float rmax, int ndex, struct HALO_LIST *halos, int resolution_level, float outer_edge);
#endif

void extDumpPointProfile(const char *fname, int nout, const DumpWorker *workers, float rmin, float rmax, int ndex, double center[3]);

/*
//  SF law from stellar particles
*/
#if defined(PARTICLES) && defined(STAR_FORMATION)
void extStarFormationLaw(const char *fname, float spatial_scale, float time_scale, float stellar_age_limit, const struct HALO_LIST *halos);
#endif

/*
//  SF law for instantaneous star formation from gas only
*/
#if defined (HYDRO) && defined(STAR_FORMATION)
void extStarFormationLaw2(const char *fname, float spatial_scale, const struct HALO_LIST *halos);
#endif

/*
//  Dump stellar particles for a given halo within a given fraction of Rvir
*/
void extHaloStars(const char *fname, const struct HALO *h, float rmax);

#endif  /* __EXT_ISM_H__ */
