#ifndef __EXT_IGM_H__
#define __EXT_IGM_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO_LIST;

#ifdef COSMOLOGY

#ifdef RADIATIVE_TRANSFER
/*
//  Save cosmic radiation background into a file. 
//  Radiation field must be initialized before this is called.
*/
void extDumpRadiationBackground(const char *filename);

/*
//  Extract radiation field at given frequencies and put it into
//  mesh variables from rt_disk_offset to rt_disk_offset+nbins-1.
//  Radiation field must be initialized before this is called.
*/
void extExtractRadiationField(int nout, const float *wlen, float *mean_rf);

/*
//  Extract photo rates field given by idx (as set in F/frt_index.h) and 
//  put them into mesh variables from rt_disk_offset to rt_disk_offset+nbins-1.
//  Radiation field must be initialized before this is called.
*/
void extExtractPhotoRates(int nout, const int *idx, float *mean_rate);

/*
//  Measure proximity zones around halos for LW band radiation and ionizing
//  radiation for HI, HeI, and HeII.
//  Radiation field must be initialized before this is called.
*/
#if defined(RT_TRANSFER) && defined(RT_EXTERNAL_BACKGROUND)
void extProximityZones(const char *fname, int resolution_level, int nside, int halo_id, const struct HALO_LIST *halos);
#endif  /* RT_TRANSFER && RT_EXTERNAL_BACKGROUND */
#endif  /* RADIATIVE_TRANSFER */

/*
//  Produce gas fractions of various gas species.
//  Density must be assigned before this is called.
//  If halos is not NULL, this call overwrites VAR_ACCEL variable.
*/
#ifdef HYDRO
void extMassFractions(const char *fname, struct HALO_LIST *halos);
#endif /* HYDRO */
#endif /* COSMOLOGY */

#endif  /* __EXT_IGM_H__ */
