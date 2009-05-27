/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2008 Krzysztof M. Gorski, Eric Hivon, 
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */

/*
//  05202009 Modified by NG: names changed to avoid name clash
*/

/* Standard Includes */
#include <math.h>

void hp_mk_pix2xy(int *pix2x, int *pix2y) {

  /* =======================================================================
   * subroutine mk_pix2xy
   * =======================================================================
   * constructs the array giving x and y in the face from pixel number
   * for the nested (quad-cube like) ordering of pixels
   *
   * the bits corresponding to x and y are interleaved in the pixel number
   * one breaks up the pixel number by even and odd bits
   * =======================================================================
   */

  int i, kpix, jpix, IX, IY, IP, ID;
  for (i = 0; i < 1023; i++) pix2x[i]=0;
  
  for( kpix=0;kpix<1024;kpix++ ) {
    jpix = kpix;
    IX = 0;
    IY = 0;
    IP = 1 ;//              ! bit position (in x and y)
    while( jpix!=0 ){// ! go through all the bits
      ID = (int)fmod(jpix,2);//  ! bit value (in kpix), goes in ix
      jpix = jpix/2;
      IX = ID*IP+IX;
      
      ID = (int)fmod(jpix,2);//  ! bit value (in kpix), goes in iy
      jpix = jpix/2;
      IY = ID*IP+IY;
      
      IP = 2*IP;//         ! next bit (in x and y)
    }
    
    pix2x[kpix] = IX;//     ! in 0,31
    pix2y[kpix] = IY;//     ! in 0,31
  }
  
  /* Later */
  return;
}


void hp_mk_xy2pix(int *x2pix, int *y2pix) {
  /* =======================================================================
   * subroutine mk_xy2pix
   * =======================================================================
   * sets the array giving the number of the pixel lying in (x,y)
   * x and y are in {1,128}
   * the pixel number is in {0,128**2-1}
   *
   * if  i-1 = sum_p=0  b_p * 2^p
   * then ix = sum_p=0  b_p * 4^p
   * iy = 2*ix
   * ix + iy in {0, 128**2 -1}
   * =======================================================================
   */
  int i, K,IP,I,J,ID;
  
  for(i = 0; i < 127; i++) x2pix[i] = 0;
  for( I=1;I<=128;I++ ) {
    J  = I-1;//            !pixel numbers
    K  = 0;//
    IP = 1;//
    truc : if( J==0 ) {
      x2pix[I-1] = K;
      y2pix[I-1] = 2*K;
    }
    else {
      ID = (int)fmod(J,2);
      J  = J/2;
      K  = IP*ID+K;
      IP = IP*4;
      goto truc;
    }
  }     
  
}


/* Standard Includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void hp_pix2ang_nest( long nside, long ipix, double *theta, double *phi) {

  /*
    c=======================================================================
    subroutine pix2ang_nest(nside, ipix, theta, phi)
    c=======================================================================
    c     gives theta and phi corresponding to pixel ipix (NESTED) 
    c     for a parameter nside
    c=======================================================================
  */
    
      int npix, npface, face_num;
      int  ipf, ip_low, ip_trunc, ip_med, ip_hi;
      int     ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
      double z, fn, fact1, fact2;
      double piover2=0.5*M_PI;
      int ns_max=8192;
      
      static int pix2x[1024], pix2y[1024];
      //      common /pix2xy/ pix2x, pix2y
      
      int jrll[12], jpll[12];// ! coordinate of the lowest corner of each face
      //      data jrll/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/ ! in unit of nside
      //      data jpll/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/ ! in unit of nside/2
      jrll[0]=2;
      jrll[1]=2;
      jrll[2]=2;
      jrll[3]=2;
      jrll[4]=3;
      jrll[5]=3;
      jrll[6]=3;
      jrll[7]=3;
      jrll[8]=4;
      jrll[9]=4;
      jrll[10]=4;
      jrll[11]=4;
      jpll[0]=1;
      jpll[1]=3;
      jpll[2]=5;
      jpll[3]=7;
      jpll[4]=0;
      jpll[5]=2;
      jpll[6]=4;
      jpll[7]=6;
      jpll[8]=1;
      jpll[9]=3;
      jpll[10]=5;
      jpll[11]=7;
      
      
      if( nside<1 || nside>ns_max ) {
	fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
	exit(0);
      }
      npix = 12 * nside*nside;
      if( ipix<0 || ipix>npix-1 ) {
	fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
	exit(0);
      }

      /* initiates the array for the pixel number -> (x,y) mapping */
      if( pix2x[1023]<=0 ) hp_mk_pix2xy(pix2x,pix2y);

      fn = 1.*nside;
      fact1 = 1./(3.*fn*fn);
      fact2 = 2./(3.*fn);
      nl4   = 4*nside;

      //c     finds the face, and the number in the face
      npface = nside*nside;

      face_num = ipix/npface;//  ! face number in {0,11}
      ipf = (int)fmod(ipix,npface);//  ! pixel number in the face {0,npface-1}

      //c     finds the x,y on the face (starting from the lowest corner)
      //c     from the pixel number
      ip_low = (int)fmod(ipf,1024);//       ! content of the last 10 bits
      ip_trunc =   ipf/1024 ;//       ! truncation of the last 10 bits
      ip_med = (int)fmod(ip_trunc,1024);//  ! content of the next 10 bits
      ip_hi  =     ip_trunc/1024   ;//! content of the high weight 10 bits

      ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
      iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

      //c     transforms this in (horizontal, vertical) coordinates
      jrt = ix + iy;//  ! 'vertical' in {0,2*(nside-1)}
      jpt = ix - iy;//  ! 'horizontal' in {-nside+1,nside-1}

      //c     computes the z coordinate on the sphere
      //      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
      jr =  jrll[face_num]*nside - jrt - 1;
      //      cout << "face_num=" << face_num << endl;
      //      cout << "jr = " << jr << endl;
      //      cout << "jrll(face_num)=" << jrll[face_num] << endl;
      //      cout << "----------------------------------------------------" << endl;
      nr = nside;//                  ! equatorial region (the most frequent)
      z  = (2*nside-jr)*fact2;
      kshift = (int)fmod(jr - nside, 2);
      if( jr<nside ) { //then     ! north pole region
         nr = jr;
         z = 1. - nr*nr*fact1;
         kshift = 0;
      }
      else {
	if( jr>3*nside ) {// then ! south pole region
         nr = nl4 - jr;
         z = - 1. + nr*nr*fact1;
         kshift = 0;
	}
      }
      *theta = acos(z);
      
      //c     computes the phi coordinate on the sphere, in [0,2Pi]
      //      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
      jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;
      if( jp>nl4 ) jp = jp - nl4;
      if( jp<1 )   jp = jp + nl4;

      *phi = (jp - (kshift+1)*0.5) * (piover2 / nr);

}


void hp_ang2pix_nest( const long nside, double theta, double phi, long *ipix) {

  /* =======================================================================
   * subroutine ang2pix_nest(nside, theta, phi, ipix)
   * =======================================================================
   * gives the pixel number ipix (NESTED) corresponding to angles theta and phi
   *
   * the computation is made to the highest resolution available (nside=8192)
   * and then degraded to that required (by integer division)
   * this doesn't cost more, and it makes sure that the treatement of round-off 
   * will be consistent for every resolution
   * =======================================================================
   */
  
  double z, za, z0, tt, tp, tmp;
  int    face_num,jp,jm;
  long   ifp, ifm;
  int    ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  double piover2 = 0.5*M_PI, pi = M_PI, twopi = 2.0*M_PI;
  int    ns_max = 8192;
  static int x2pix[128], y2pix[128];
  static char setup_done = 0;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  if( theta<0. || theta>pi ) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
  if( !setup_done ) {
    hp_mk_xy2pix(x2pix,y2pix);
    setup_done = 1;
  }
  
  z  = cos(theta);
  za = fabs(z);
  z0 = 2./3.;
  if( phi>=twopi ) phi = phi - twopi;
  if( phi<0. )    phi = phi + twopi;
  tt = phi / piover2; /* in [0,4[ */
  
  if( za<=z0 ) { /* equatorial region */
    
    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floor(ns_max*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floor(ns_max*(0.5 + tt + z*0.75)); /* descending edge line index */
    
    /* finds the face */
    ifp = jp / ns_max; /* in {0,4} */
    ifm = jm / ns_max;
    
    if( ifp==ifm ) face_num = (int)fmod(ifp,4) + 4; /* faces 4 to 7 */
    else if( ifp<ifm ) face_num = (int)fmod(ifp,4); /* (half-)faces 0 to 3 */
    else face_num = (int)fmod(ifm,4) + 8;           /* (half-)faces 8 to 11 */
    
    ix = (int)fmod(jm, ns_max);
    iy = ns_max - (int)fmod(jp, ns_max) - 1;
  }
  else { /* polar region, za > 2/3 */
    
    ntt = (int)floor(tt);
    if( ntt>=4 ) ntt = 3;
    tp = tt - ntt;
    tmp = sqrt( 3.*(1. - za) ); /* in ]0,1] */
    
    /* (the index of edge lines increase when distance from the closest pole
     * goes up)
     */
    /* line going toward the pole as phi increases */
    jp = (int)floor( ns_max * tp          * tmp ); 

    /* that one goes away of the closest pole */
    jm = (int)floor( ns_max * (1. - tp) * tmp );
    jp = (int)(jp < ns_max-1 ? jp : ns_max-1);
    jm = (int)(jm < ns_max-1 ? jm : ns_max-1);
    
    /* finds the face and pixel's (x,y) */
    if( z>=0 ) {
      face_num = ntt; /* in {0,3} */
      ix = ns_max - jm - 1;
      iy = ns_max - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
    }
  }
  
  ix_low = (int)fmod(ix,128);
  ix_hi  =     ix/128;
  iy_low = (int)fmod(iy,128);
  iy_hi  =     iy/128;

  ipf = (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)+ (x2pix[ix_low]+y2pix[iy_low]);
  ipf = (long)(ipf / pow(ns_max/nside,2));     /* in {0, nside**2 - 1} */
  *ipix =(long)( ipf + face_num*pow(nside,2)); /* in {0, 12*nside**2 - 1} */
}
