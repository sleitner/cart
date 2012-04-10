#ifndef __IMF_H__
#define __IMF_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#ifdef STAR_FORMATION

struct InitialMassFunction
{
  const char *name;      /* name of IMF function */
  double min_mass;       /* used to be called aM_stl */
  double max_mass;       /* used to be called aM_stu */
  double min_SNII_mass;  /* used to be called aM_SNII */
  double min_SNIa_mass;  /* used to be called aM_SNIa1 */
  double max_SNIa_mass;  /* used to be called aM_SNIa2 */
  double (*f)( double mstar );
  double (*fm)( double mstar );
  double (*fmz)( double mstar );
};


extern const struct InitialMassFunction *imf;

typedef double(*fimf)(double);
struct IMF_t
{
  char* name;
  fimf  f;
};

void config_init_imf();
void config_verify_imf();

#endif /* STAR_FORMATION */

#endif
