#include "interpolation.h"
#include "auxiliary.h"

/**************************/
/* R. Feldmann 2011       */
/**************************/

int inter_bisect(real x, real *xx, int n, int mm) {
  int ju,jm,jl;
  
  if (n<2 || mm<2 || mm>n) {
    printf("error: locate");
    exit(1);
  }
  
  int ascnd = (xx[n-1] >= xx[0]);
  jl=0; ju=n-1;
  while (ju-jl>1) {
    jm  = (ju+jl)>>1;
    if (x>=xx[jm] ==ascnd)
      jl=jm;
    else
      ju=jm;
  }
  return MAX( 0, MIN(n-mm, jl-((mm-2)>>1)) );
}

real inter_linear_1d(real x, 
                     real *xx, int n, real *yy) {
  int i = inter_bisect(x,xx,n,2);
  double a;
  
  a = (x-xx[i])/(xx[i+1]-xx[i]);
  
  /* no extrapolation */
  if (a<0) a=0;
  if (a>1) a=1;
  
  double na = a*yy[i+1];
  double ia = (1.-a)*yy[i];
  
  if (a<=0) na=0.;
  if (a>=1) ia=0.;
  
  return ia + na;
}

real inter_linear_2d(real x, real y,
                     real *xx, int nx, real *yy, int ny, real *mat) {
  int i,j;
  double a,b;
  i = inter_bisect(x,xx,nx,2);
  j = inter_bisect(y,yy,ny,2);
  
  a = (x-xx[i])/(xx[i+1]-xx[i]);
  b = (y-yy[j])/(yy[j+1]-yy[j]);
  
  if (a<0) a=0.;
  if (a>1) a=1.;
  if (b<0) b=0.;
  if (b>1) b=1.;
  
  double iaib = (1.-a)*(1.-b)*mat[IND2(i,j,ny)];
  double aib  = a*(1.-b)*mat[IND2(i+1,j,ny)];
  double iab  = (1.-a)*b*mat[IND2(i,j+1,ny)];
  double  ab  = a*b*mat[IND2(i+1,j+1,ny)];
 
  if (a<=0) aib=ab=0.;
  if (a>=1) iaib=iab=0.;
  if (b<=0) iab=ab=0.;
  if (b>=1) iaib=aib=0.;
  
  return iaib + aib + iab + ab;  
}

real inter_linear_3d(real x, real y, real z,
                     real *xx, int nx, real *yy, int ny, real *zz, int nz, 
		     real *mat) {
  int i,j,k;
  double a,b,c;
  i = inter_bisect(x,xx,nx,2);
  j = inter_bisect(y,yy,ny,2);
  k = inter_bisect(z,zz,nz,2);
  
  a = (x-xx[i])/(xx[i+1]-xx[i]);
  b = (y-yy[j])/(yy[j+1]-yy[j]);
  c = (z-zz[k])/(zz[k+1]-zz[k]);
  
  if (a<0) a=0.;
  if (a>1) a=1.;
  if (b<0) b=0.;
  if (b>1) b=1.;
  if (c<0) c=0.;
  if (c>1) c=1.;
  
  double iaibic = (1.-a)*(1.-b)*(1.-c)*mat[IND3(i,j,k,ny,nz)];
  double aibic = a*(1.-b)*(1.-c)*mat[IND3(i+1,j,k,ny,nz)];
  double iabic = (1.-a)*b*(1.-c)*mat[IND3(i,j+1,k,ny,nz)];
  double iaibc = (1.-a)*(1.-b)*c*mat[IND3(i,j,k+1,ny,nz)];
  double aibc  =  a*(1.-b)*c*mat[IND3(i+1,j,k+1,ny,nz)];
  double iabc  = (1.-a)*b*c*mat[IND3(i,j+1,k+1,ny,nz)];
  double abic  = a*b*(1.-c)*mat[IND3(i+1,j+1,k,ny,nz)];
  double abc   = a*b*c*mat[IND3(i+1,j+1,k+1,ny,nz)];
  
  if (a<=0) aibic=aibc=abic=abc=0.;
  if (a>=1) iaibic=iabic=iaibc=iabc=0.;
  if (b<=0) iabic=iabc=abic=abc=0.;
  if (b>=1) iaibic=aibic=iaibc=aibc=0.;
  if (c<=0) iaibc=aibc=iabc=abc=0.;
  if (c>=1) iaibic=aibic=iabic=abic=0.;

  return  iaibic + aibic + iabic + iaibc + aibc + iabc + abic + abc;
}

real inter_linear_4d(real x, real y, real z, real t,
                     real *xx, int nx, real *yy, int ny, real *zz, int nz, real *tt, int nt, 
		     real *mat) {
  int i,j,k,l;
  double a,b,c,d;
  i = inter_bisect(x,xx,nx,2);
  j = inter_bisect(y,yy,ny,2);
  k = inter_bisect(z,zz,nz,2);
  l = inter_bisect(t,tt,nt,2);
  
  a = (x-xx[i])/(xx[i+1]-xx[i]);
  b = (y-yy[j])/(yy[j+1]-yy[j]);
  c = (z-zz[k])/(zz[k+1]-zz[k]);
  d = (t-tt[l])/(tt[l+1]-tt[l]);
  
  if (a<0) a=0.;
  if (a>1) a=1.;
  if (b<0) b=0.;
  if (b>1) b=1.;
  if (c<0) c=0.;
  if (c>1) c=1.;
  if (d<0) d=0.;
  if (d>1) d=1.;  
  
  double iaibicid = (1.-a)*(1.-b)*(1.-c)*(1.-d)*mat[IND4(i,j,k,l,ny,nz,nt)];  
  
  double aibicid =    a   *(1.-b)*(1.-c)*(1.-d)*mat[IND4(i+1,j,k,l,ny,nz,nt)];
  double iabicid = (1.-a)*   b   *(1.-c)*(1.-d)*mat[IND4(i,j+1,k,l,ny,nz,nt)];  
  double iaibcid = (1.-a)*(1.-b)*   c   *(1.-d)*mat[IND4(i,j,k+1,l,ny,nz,nt)]; 
  double iaibicd = (1.-a)*(1.-b)*(1.-c)  * d   *mat[IND4(i,j,k,l+1,ny,nz,nt)];
   
  double abicid =      a  *   b  *(1.-c)*(1.-d)*mat[IND4(i+1,j+1,k,l,ny,nz,nt)];  
  double aibcid =      a  *(1.-b)*   c  *(1.-d)*mat[IND4(i+1,j,k+1,l,ny,nz,nt)];  
  double aibicd =      a  *(1.-b)*(1.-c)*   d  *mat[IND4(i+1,j,k,l+1,ny,nz,nt)];  
  
  double iabcid =   (1.-a)*  b   *   c  *(1.-d)*mat[IND4(i,j+1,k+1,l,ny,nz,nt)];  
  double iabicd =   (1.-a)*  b   *(1.-c)*  d   *mat[IND4(i,j+1,k,l+1,ny,nz,nt)];  
  double iaibcd =   (1.-a)*(1.-b)*   c  *  d   *mat[IND4(i,j,k+1,l+1,ny,nz,nt)];  
  
  double abcid =      a   *   b  *   c  *(1.-d)*mat[IND4(i+1,j+1,k+1,l,ny,nz,nt)];  
  double abicd =      a   *   b  *(1.-c)*   d  *mat[IND4(i+1,j+1,k,l+1,ny,nz,nt)];  
  double aibcd =      a   *(1.-b)*   c  *   d  *mat[IND4(i+1,j,k+1,l+1,ny,nz,nt)];  
  double iabcd =    (1.-a)*   b  *   c  *   d  *mat[IND4(i,j+1,k+1,l+1,ny,nz,nt)];  
  
  double abcd =       a   *   b  *   c  *   d  *mat[IND4(i+1,j+1,k+1,l+1,ny,nz,nt)];  
  
  if (a<=0) aibicid=abicid=aibcid=aibicd=abcid=abicd=aibcd=abcd=0;
  if (a>=1) iaibicid=iabicid=iaibcid=iaibicd=iabcid=iabicd=iaibcd=iabcd=0;
  if (b<=0) iabicid=abicid=iabcid=iabicd=abcid=abicd=iabcd=abcd=0;
  if (b>=1) iaibicid=aibicid=iaibcid=iaibicd=aibcid=aibicd=iaibcd=aibcd=0;
  if (c<=0) iaibcid=aibcid=iabcid=iaibcd=abcid=aibcd=iabcd=abcd=0;
  if (c>=1) iaibicid=aibicid=iabicid=iaibicd=abicid=aibicd=iabicd=abicd=0;
  if (d<=0) iaibicd=aibicd=iabicd=iaibcd=abicd=aibcd=iabcd=abcd=0;
  if (d>=1) iaibicid=aibicid=iabicid=iaibcid=abicid=aibcid=iabcid=abcid=0;

  return  iaibicid + aibicid + iabicid + iaibcid + iaibicd + abicid + aibcid + aibicd 
          + iabcid + iabicd + iaibcd + abcid + abicd + aibcd + iabcd + abcd;
}

/* linear interpolation on a regular, cartesian grid in arbitrary dimension */
/* Note: this code should only be used for moderately high n, since it invokes 2^(n-1) 1-dim interpolations */
real inter_linear_nd(real *x,   int n,
                     real **xx, int *nx, 
		     real *mat) {

  /* x is a pointer to the n variables (x^1,x^2,x^3,..,x^n) */
  /* xx is a pointer to the n 1-dimensional fields that specify the grid locations */
  /* nx is a pointer to the number of grid points in each dimension */
  /* mat is a pointer to the function values on the grid points */
  /* (same ordering as for any of the inter_linear_?d functions) */

  /* bisect in dimension n */
  int i = inter_bisect(x[0],xx[0],nx[0],2);
  double a = (x[0]-xx[0][i])/(xx[0][i+1]-xx[0][i]);
  int j;
  int prod_nx;
  double I0,I1;
  
  /* no extrapolation */
  if (a<0) a=0;
  if (a>1) a=1;

  /* interpolate in n-1 dimension */
  if (n>1) {
    prod_nx=1;
    for (j=1;j<n;j++) {
      prod_nx*=nx[j];
    }
  
    I0 = inter_linear_nd(x+1, n-1, xx+1, nx+1, mat + prod_nx*i );
    I1 = inter_linear_nd(x+1, n-1, xx+1, nx+1, mat + prod_nx*(i+1) );
  } else {
    I0 = inter_linear_1d(x[0], xx[0], nx[0], mat);
    I1 = inter_linear_1d(x[0], xx[0], nx[0], mat);
  }

  return a*I1 + (1.-a)*I0;
}
