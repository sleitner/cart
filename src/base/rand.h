#ifndef __RAND_H__
#define __RAND_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif

void init_rand();
void cart_rand_set_seed( unsigned long int seed );
void cart_rand_load_state( const char *state_filename, int fail_on_error );
void cart_rand_save_state( const char *state_filename );

double cart_rand();
unsigned long int cart_rand_ulong(unsigned long int nmax);
double cart_rand_lognormal(double sigma);
unsigned int cart_rand_poisson(double mu);
void cart_rand_unit_vector(double uni[nDim]);


#endif
