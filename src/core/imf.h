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
