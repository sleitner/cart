#ifndef __EXT_UTILS_H__
#define __EXT_UTILS_H__


#ifndef CONFIGURED
#error "Missing config.h include."
#endif


/*
//  Miscelaneous helper functions
*/
void find_max_var(int var, float *val, double *pos, double dist);

void cell_interpolate_at_position(int cell, double pos[], int nvars, int vars[], float vals[]);

#endif  /* __EXT_UTILS_H__ */
