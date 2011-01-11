#ifndef __EXT_GIC_TOOLS_H__
#define __EXT_GIC_TOOLS_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct HALO;


#define GIC_MASK_MODE_PLAIN  0
#define GIC_MASK_MODE_CHULL  1


void gicMakeMask(const char *filename, int num_halos, const struct HALO **halos, float size, int mode);


#endif  /* __EXT_GIC_TOOLS_H__ */