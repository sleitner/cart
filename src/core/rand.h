#ifndef __RAND_H__
#define __RAND_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


void init_rand();
void save_rand();


double cart_rand();
double cart_rand_lognormal(double sigma);


#endif
