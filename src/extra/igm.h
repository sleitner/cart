#ifndef __EXT_IGM_H__
#define __EXT_IGM_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO_LIST;


void extDumpRadiationBackground();
void extExtractRadiationField(int nbins, const float wbins[], float *mean_rf);
void extFindProximityZones(const char *fname, int floor_level, int nside, int halo_id, const struct HALO_LIST *halos);

void extDumpGasFractions(const char *fname, const struct HALO_LIST *halos);

#endif  /* __EXT_IGM_H__ */
