/*=========================================================================

Copyright (c) 2011 Nick Gnedin
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Nick Gnedin nor the names of any contributors may be used
   to endorse or promote products derived from this software without specific
   prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif


#include "ffti.h"
#include "ffts.h"



/*
//  Types for C - Fortran interface
*/
#define f77_call(fun)   fun##_

#ifdef FFT_DOUBLE
#define rffti1 dffti1
#define rfftf1 dfftf1
#define rfftb1 dfftb1
#define cffti1 zffti1
#define cfftf1 zfftf1
#define cfftb1 zfftb1
#endif

void f77_call(rffti1)(int *n, fft_t *work, int *facs);
void f77_call(rfftf1)(int *n, fft_t *data, fft_t *temp, fft_t *work, int *facs);
void f77_call(rfftb1)(int *n, fft_t *data, fft_t *temp, fft_t *work, int *facs);

void f77_call(cffti1)(int *n, fft_t *work, int *facs);
void f77_call(cfftf1)(int *n, fft_t *data, fft_t *temp, fft_t *work, int *facs);
void f77_call(cfftb1)(int *n, fft_t *data, fft_t *temp, fft_t *work, int *facs);


struct fftiInternal
{
  int dims[3];
  int *facs[3];
  fft_t *work[3];
  fft_t **tomp;
  int nomp;
}
ffti_internal_data;
#define d ffti_internal_data


void ffti_init(int dims[3])
{
  int i, len;

  for(i=0; i<3; i++) d.dims[i] = dims[i];

  d.facs[0] = (int *)malloc(15*sizeof(int));
  d.facs[1] = (int *)malloc(15*sizeof(int));
  d.facs[2] = (int *)malloc(15*sizeof(int));

  if(d.facs[0]==NULL || d.facs[1]==NULL || d.facs[2]==NULL)
    {
      FAIL("ffti: unable to allocate temporary storage.");
    }

  d.work[0] = NEW(1*d.dims[0]);
  d.work[1] = NEW(2*d.dims[1]);
  d.work[2] = NEW(2*d.dims[2]);

  if(d.work[0]==NULL || d.work[1]==NULL || d.work[2]==NULL)
    {
      FAIL("ffti: unable to allocate temporary storage.");
    }

  f77_call(rffti1)(d.dims+0,d.work[0],d.facs[0]);
  f77_call(cffti1)(d.dims+1,d.work[1],d.facs[1]);
  f77_call(cffti1)(d.dims+2,d.work[2],d.facs[2]);

#ifdef _OPENMP
  d.nomp = omp_get_max_threads();
#else
  d.nomp = 1;
#endif

  d.tomp = (fft_t**)malloc(sizeof(fft_t*)*d.nomp);
  if(d.tomp == NULL)
    {
      FAIL("ffti: unable to allocate temporary storage.");
    }

  len = d.dims[0];
  if(len < 2*d.dims[1]) len = 2*d.dims[1];
  if(len < 2*d.dims[2]) len = 2*d.dims[2];

  for(i=0; i<d.nomp; i++)
    {
      d.tomp[i] = NEW(len);
      if(d.tomp[i] == NULL)
	{
	  FAIL("ffti: unable to allocate temporary storage.");
	}
    }
}


void ffti_done()
{
  int i;

  for(i=0; i<3; i++)
    {
      DEL(d.work[i]);
      free(d.facs[i]);
    }

  for(i=0; i<d.nomp; i++) DEL(d.tomp[i]);
  free(d.tomp);
}


void ffti_r2c_x(int dir, fft_t *array)
{
  int i;
#ifdef _OPENMP
  int iomp = omp_get_thread_num();
#else
  int iomp = 0;
#endif
  
  if(iomp<0 || iomp>=d.nomp)
    {
      FAIL("ffti: invalid OpenMP thread.");
    }

  if(dir)
    {
      f77_call(rfftf1)(d.dims+0,array,d.tomp[iomp],d.work[0],d.facs[0]);

      /*
      //  Unwrap the k-space representation
      */
      array[d.dims[0]+1] = 0.0;
      for(i=d.dims[0]; i>1; i--)
	{
	  array[i] = array[i-1];
	}
      array[1] = 0.0;
    }
  else
    {
      /*
      //  Wrap the k-space representation
      */
      for(i=1; i<=d.dims[0]; i++)
	{
	  array[i] = array[i+1];
	}

      f77_call(rfftb1)(d.dims+0,array,d.tomp[iomp],d.work[0],d.facs[0]);
    }
}


void ffti_internal_c2c(int dir, int axis, fft_t *array);


void ffti_c2c_y(int dir, fft_t *array)
{
  ffti_internal_c2c(dir,1,array);
}


void ffti_c2c_z(int dir, fft_t *array)
{
  ffti_internal_c2c(dir,2,array);
}


void ffti_internal_c2c(int dir, int axis, fft_t *array)
{
#ifdef _OPENMP
  int iomp = omp_get_thread_num();
#else
  int iomp = 0;
#endif

  if(iomp<0 || iomp>=d.nomp)
    {
      FAIL("ffti: invalid OpenMP thread.");
    }

  if(dir)
    {
      f77_call(cfftf1)(d.dims+axis,array,d.tomp[iomp],d.work[axis],d.facs[axis]);
    }
  else
    {
      f77_call(cfftb1)(d.dims+axis,array,d.tomp[iomp],d.work[axis],d.facs[axis]);
    }
}


