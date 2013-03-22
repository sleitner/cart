#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "auxiliary.h"
#include "cosmology.h"
#include "hydro.h"
#include "rand.h"
#include "system.h"
#include "times.h"
#include "tree.h"
#include "units.h"

#include "ism.h"
#include "los.h"
#include "synthetic_spectra.h"


#if defined(HYDRO) && defined(RADIATIVE_TRANSFER) && defined(COSMOLOGY)


#ifndef SS_MAX_VARS
#define SS_MAX_VARS  2
#endif


typedef struct
{
  double r;
  double vlos;
  float vals[SS_MAX_VARS];
}
ssPoint;


typedef struct
{
  int NumVars, NumPts, Size;
  double e[3];
  DumpWorker Vars[SS_MAX_VARS];
  ssPoint* Data;
}
ssLine;


int ssWorker(int id, int cell, double r1, double r2, losBuffer data)
{
  int j;
  ssLine *d = (ssLine*)data.Data;

  ssPoint p;
  for(j=0; j<d->NumVars; j++)
    {
      p.vals[j] = d->Vars[j].Value(cell_level(cell),cell,0,0);
    }
  p.r = 0.5*(r2+r1);
  p.vlos = 0.0;
  for(j=0; j<nDim; j++) p.vlos += cell_momentum(cell,j)*d->e[j];
  p.vlos /= cell_gas_density(cell);

  if(d->Size == d->NumPts)
    {
      ssPoint *tmp;

      d->Size = (int)(d->Size*1.5+100);
      tmp = cart_alloc(ssPoint,d->Size);
      if(d->NumPts > 0)
	{
	  memcpy(tmp,d->Data,d->NumPts*sizeof(ssPoint));
	  cart_free(d->Data);
	}
      d->Data = tmp;
    }
  
  d->Data[d->NumPts++] = p;

  return 0;
}


int ssPointCmp(const void *p1, const void *p2)
{
  const ssPoint *v1 = (const ssPoint*)p1;
  const ssPoint *v2 = (const ssPoint*)p2;
  if(v1->r > v2->r) return 1; else if(v1->r < v2->r) return -1; else return 0;
}


void ssCollector1(const losBuffer *result, int num_segments, const losSegment *segments)
{
  int i, npts = 0;

  for(i=0; i<num_segments; i++)
    {
      const ssLine *d = (const ssLine*)segments[i].Buffer.Data;
      cart_assert(d != 0);
      npts += d->NumPts;
    }

  *((int *)result->Data) = npts;
}


void ssCollector2(const losBuffer *result, int num_segments, const losSegment *segments)
{
  int i, j, n, npts = 0;
  ssPoint *arr = (ssPoint*)result->Data;

  for(i=0; i<num_segments; i++)
    {
      const ssPoint *d = (const ssPoint *)segments[i].Buffer.Data;
      n = segments[i].Buffer.Size/sizeof(ssPoint);
      cart_assert(n==0 || d!=0);

      for(j=0; j<n; j++)
	{
	  arr[npts++] = d[j];
	}
    }

  qsort(arr,npts,sizeof(ssPoint),ssPointCmp);
}


float ss_tem(int level, int cell, double *ref_pos, double *ref_vel)
{
  return units->temperature*cell_gas_temperature(cell);
}


float ss_nHI(int level, int cell, double *ref_pos, double *ref_vel)
{
  return units->number_density*cell_HI_density(cell);
}


ssSpectrum ssMakeSpectrum_HI1216(double pos[3], double theta, double phi, double len, double vel_pixel, int floor_level)
{
  static const double BWID = 5.0;
  static const double wb0 = 0.12856*1.5e2; /* 1.5e2 = sqrt(2.0e4) */

  ssLine line;
  int i, n, nn;
  double vb, vd, fact0 = 4.478e-18*2.998e5/cosmology->h*3.081e24*auni[min_level];
  double xtl, xtr, vtl, vtr, dx, dv, facts;
  int l, ls, ll, lu;
  double aH = auni[min_level]*Hubble(auni[min_level])/cosmology->h;  // this H per CHIMP
  double *vt, wstp, wmin, wmax;
  double *w, *f;
  losSegment seg;
  losBuffer res;
  ssSpectrum spec;
  ssPoint *arr;

  /*
  //  Spectrum specialization
  */
  line.NumVars = 2;
  line.Vars[1].Value = ss_tem;
  line.Vars[0].Value = ss_nHI;

  line.e[0] = sin(theta)*cos(phi);
  line.e[1] = sin(theta)*sin(phi);
  line.e[2] = cos(theta);

  line.NumPts = line.Size = 0;
  line.Data = NULL;

  seg.Buffer.Size = sizeof(ssLine);
  seg.Buffer.Data = &line;

  losTraverseSegment(0,pos,theta,phi,len,floor_level,ssWorker,&seg);
 
  res.Size = sizeof(int);
  res.Data = &n;
  losCollectSegments(&res,&seg,ssCollector1);

  if(n <= 0)
    {
      cart_debug("ss:n = %d",n);
      cart_debug("ss:pos = %lg %lg %lg",pos[0],pos[1],pos[2]);
      cart_debug("ss:ang = %lg %lg",theta,phi);
      cart_debug("ss:len = %lg",len);
      spec.Length = 0;
      spec.Vel = NULL;
      spec.Flux = NULL;
      return spec;
    }
  arr = cart_alloc(ssPoint,n);

  seg.Buffer.Size = line.NumPts*sizeof(ssPoint);
  seg.Buffer.Data = line.Data;
  res.Size = sizeof(ssPoint)*n;
  res.Data = arr;
  losCollectSegments(&res,&seg,ssCollector2);
  if(line.Data != NULL) cart_free(line.Data);

  vt = cart_alloc(double,n);
  for(i=0; i<n; i++)
    {
      vt[i] = aH*arr[i].r + arr[i].vlos*units->velocity/constants->kms;
    }

  wstp = vel_pixel;
  wmin = vt[0];
  wmax = vt[0];
  for(i=1; i<n; i++)
    {
      wmin = MIN(wmin,vt[i]);
      wmax = MAX(wmax,vt[i]);
    }

  wmin += BWID*wb0;
  wmax -= BWID*wb0;
  nn = 1 + (int)(0.5+(wmax-wmin)/wstp);
  if(nn < 1) nn = 1;

  w = cart_alloc(double,nn);
  f = cart_alloc(double,nn);
  for(i=0; i<nn; i++)
    {
      w[i] = wmin + wstp*i;
      f[i] = 0.0;
    }

  for(i=1; i<n; i++)
    {
      vb = 0.12856*sqrt(0.5*(arr[i-1].vals[1]+arr[i].vals[1]));
      vd = 3.0e5*2.02e-8;

      xtl = arr[i-1].r;
      xtr = arr[i].r;

      vtl = MIN(vt[i-1],vt[i]);
      vtr = MAX(vt[i-1],vt[i]);

      dx = xtr - xtl;
      dv = vtr - vtl;

      ls = 1 + (int)(BWID*vb/wstp);
      ll = 1 + (int)(0.5+(vtl-wmin)/wstp) - ls;
      lu = 1 + (int)(0.5+(vtr-wmin)/wstp) + ls;
      if(ll <  0) ll = 0;
      if(ll >= nn) ll = nn - 1;
      if(lu <  0) lu = 0;
      if(lu >= nn) lu = nn - 1;

      facts = fact0*0.5*(arr[i-1].vals[0]+arr[i].vals[0]);

      if(dv > 0.001*vb)
	{
	  /*
	  // Not a caustic
	  */
	  facts = facts*(dx/dv);
	  for(l=ll; l<=lu; l++)
	    {
	      f[l] += facts*voigtInt(vtl-w[l],vtr-w[l],vb,vd);
	    }
	}
      else
	{
	  /*
	  //  Velocity caustics
	  */
	  facts = facts*dx/(vb*1.77245385);
	  for(l=ll; l<=lu; l++)
	    {
	      f[l] += facts*voigt(w[l]-vtl,vb,vd);
            }
	}
    }

  for(l=0; l<nn; l++)
    {
      f[l] = exp(-f[l]);
    }

  cart_free(arr);
  cart_free(vt);

  spec.Length = nn;
  spec.Vel = w;
  spec.Flux = f;

  return spec;
}

#endif /* HYDRO && RADIATIVE_TRANSFER && COSMOLOGY */


double voigt(double w, double wb, double wl)
{
  // 0.832554611 = sqrt(ln(2))

  double wd = 0.832554611*wb;
  double wv = 0.5346*wl + sqrt(0.2166*wl*wl+wd*wd);
  double d = (wl-wd)/(wl+wd);

  double cl = 0.68188 + d*( 0.61293 + d*(-0.18384 - 0.11568*d));
  double cg = 0.32460 + d*(-0.61825 + d*( 0.17681 + 0.12109*d));

  double x = (w/wv);
  x *= x;

  return (cl/3.1415927/(x+1)+cg*0.832554611/1.77245385*exp(-0.693147181*x))/wv;
}


double voigtInt(double w1, double w2, double wb, double wl)
{
  // 0.832554611 = sqrt(ln(2))

  double wd = 0.832554611*wb;
  double wv = 0.5346*wl + sqrt(0.2166*wl*wl+wd*wd);
  double d = (wl-wd)/(wl+wd);

  double cl = 0.68188 + d*(0.61293 + d*(-0.18384 - 0.11568*d));
  double cg = 0.32460 + d*(-0.61825 + d*(0.17681 + 0.12109*d));

  double x1 = w1/wv;
  double x2 = w2/wv;

  return cl*3.1415927*(atan(x2)-atan(x1)) + cg*0.5*(erf(0.832554611*x2)-erf(0.832554611*x1));
}

