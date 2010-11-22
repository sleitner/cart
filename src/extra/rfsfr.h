#ifndef __EXT_RFSFR_H__
#define __EXT_RFSFR_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO_LIST;

/*
//  Dump ISRF and some variables with a hierarchy of levels
*/
void extRFvsSFR(const char *froot, int top_level, const struct HALO_LIST *halos);

#endif  /* __EXT_RFSFR_H__ */
