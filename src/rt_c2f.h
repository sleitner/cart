#ifndef __RT_C2F_H__
#define __RT_C2F_H__

/*  Types for C - Fortran interface */

#define c2f_cell	integer*4
#define c2f_int		integer*4
#define c2f_float	real*4
#define c2f_double	real*8

#define f2c_intg	int
#define f2c_real	float

#define f2c_wrapper(fun) fun##_


#endif  /* __RT_C2F_H__ */
