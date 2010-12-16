#ifndef __ENGINE_NG_ANALYSIS_H__
#define __ENGINE_NG_ANALYSIS_H__


typedef struct NG_VAR_LIST_TYPE
{
  int Num;
  int *Ids;
}
ng_var_list_t;


const char *ngOutputFile(const char *filepath);


/* Load hlist file for other command to use. */
void ngLoadHalos(const char *dirname);

/* Dump chemical state of the gas at levels <level> and below. */
void ngDumpLevels(const char *filename, int level);

/* Dump profiles of gas quantities from <rmin> to <rmax> (in kpc) 
// for halos resolved to at least <resolution_level>. */
void ngDumpProfiles(const char *filename, int resolution_level, float rmin, float rmax);

/* Dump gas fractions for halos. */
void ngGasFractions(const char *filename);

/* Produce IFrIR uniform scalars file of <nbin>^3 size, centered at <pos>, resolved to <floor_level>.
// If <list> is NULL, use the default variables to save in the file. */
void ngRegion2IFrIT(const char *filename, int nbin, int floor_level, const double pos[], const ng_var_list_t *list);

/* As above, but now center at a given halo with id=<id>. */
void ngHalo2IFrIT(const char *filename, int floor_level, int id, const ng_var_list_t *list);

/* Compute proximity zones for a specific halo with id=<id>, or for all halos, if <id>=0.
// Only consider halos resolved to at least <resolution_level>. */
void ngProximityZone(const char *filename, int id, int resolution_level);

/* Output radiation field as a function of SFR for all levels at or below <top_level>. */
void ngRFvsSFR(const char *fileroot, int top_level);

/* Dump star formation law (i.e. Kennicutt-Schmidth relation) on scale <spatial_scale> kpc, averaged 
// over <time_scale> Myr, and limit stellar age in output M* to <stellar_age_limit> Myr. */
void ngStarFormationLaw(const char *filename, float spatial_scale, float time_scale, float stellar_age_limit);

/* Same as above, but compute instantaneous SFR. */
void ngStarFormationLaw2(const char *filename, float spatial_scale);

/* Dump ids of all particles within <rmax>*Rvir of a given halo with id <id>. */
void ngHaloParticles(const char *filename, int id, float rmax);

/* Dump properties of stellar particles within <rmax>*Rvir of a given halo with id <id>. */
void ngHaloStars(const char *filename, int id, float rmax);

#endif  /* __ENGINE_NG_ANALYSIS_H__ */
