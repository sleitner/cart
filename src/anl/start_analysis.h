#ifndef __START_ANALYSIS_H__
#define __START_ANALYSIS_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


struct init_t
{
  void (*All)();
  void (*TotalDensity)();
  void (*StellarDensity)();
  void (*RadiativeTransfer)();
};
extern const struct init_t *init;


struct snapshot_t
{
  int Current;
  int Number;
};
extern const struct snapshot_t *snapshot;


#ifdef GRAVITY
#define VAR_TOTAL_DENSITY            VAR_TOTAL_MASS
#define cell_total_density(cell)     cell_total_mass(cell)
#ifdef STARFORM
#define VAR_STELLAR_DENSITY          VAR_FIRST_SPECIES_MASS
#define cell_stellar_density(cell)   cell_first_species_mass(cell)
#endif /* STARFORM */
#endif /* GRAVITY */

#endif  /* __START_ANALYSIS_H__ */
