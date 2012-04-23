#ifndef __EXT_GIC_TOOLS_H__
#define __EXT_GIC_TOOLS_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#define GIC_MASK_MODE_PLAIN  0
#define GIC_MASK_MODE_CHULL  1


struct gicCHullConfig
{
  float Particles;
  float Volume;
  int   NumSteps;
  int   Loud;
};

extern struct gicCHullConfig gictol;


#ifdef PARTICLES
void gicMakeMask(const char *filename, const halo_list *halos, float size, int mode, int level, int width);
#endif

#endif  /* __EXT_GIC_TOOLS_H__ */
