#ifndef __EXT_SPS_H__
#define __EXT_SPS_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO_LIST;


/*
//  Compute luminosities in a given band from an FSPS table
*/
void extGalLums(const char *fname, const struct HALO_LIST *halos, float rmax, float tSFR, const char *spsfile);

#endif  /* __EXT_SPS_H__ */
