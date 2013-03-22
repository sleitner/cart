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

double tlf( double logM, double logZ ); /* log_10(lifetime/yr) of a star of mass log_10(M/Msun) */
                                             /*  and log_10(Z) (absolute ox. abundance) */
double mlf( double logt, double logZ ); /* log_10(mass/Msun) of a star with lifetime log_10(t/yr) */
                                             /*  and log_10(Z) (absolute ox. abundance) */

#endif /* STAR_FORMATION */

#endif
