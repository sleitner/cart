#include <stdlib.h>
#include <stdio.h>

/**************************/
/* R. Feldmann 2011       */
/**************************/

typedef float real;

#define IND2(X,Y,NY) ( (Y) + (NY)*(X) )

#define IND3(X,Y,Z,NY,NZ) ( (Z) + (NZ)*( (Y) + (NY)*(X) ) )

#define IND4(X,Y,Z,T,NY,NZ,NT) (  (T) + (NT)*( (Z) + (NZ)*( (Y) + (NY)*(X) ) ) )

int inter_bisect(real x, real *xx, int n, int mm);

real inter_linear_1d(real x, 
                     real *xx, int n, real *yy);

real inter_linear_2d(real x, real y,
                     real *xx, int nx, real *yy, int ny, real *mat);

real inter_linear_3d(real x, real y, real z,
                     real *xx, int nx, real *yy, int ny, real *zz, int nz, 
		     real *mat);

real inter_linear_4d(real x, real y, real z, real t,
                     real *xx, int nx, real *yy, int ny, real *zz, int nz, real *tt, int nt, 
		     real *mat);

real inter_linear_nd(real *x,   int n,
                     real **xx, int *nx, 
		     real *mat);	     
		     
