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



/*
//  Types for C - Fortran interface
*/
#define f77_call(fun)   fun##_

#ifdef FFT_DOUBLE
#define cffti1 zffti1
#define cfftf1 zfftf1
#define cfftb1 zfftb1
#endif

void f77_call(cffti1)(int *n, fft_t *work, int *facs);
void f77_call(cfftf1)(int *n, fft_t *data, fft_t *temp, fft_t *work, int *facs);
void f77_call(cfftb1)(int *n, fft_t *data, fft_t *temp, fft_t *work, int *facs);


int ffti_get_work_size(int n)
{
  return 4*n + 15;
}


void ffti_init(int n, fft_t *work)
{
  f77_call(cffti1)(&n,work+2*n,(int *)(work+4*n));
}


void ffti_fc2c(int n, fft_t *array, fft_t *work)
{
  f77_call(cfftf1)(&n,array,work,work+2*n,(int *)(work+4*n));
}


void ffti_bc2c(int n, fft_t *array, fft_t *work)
{
  f77_call(cfftb1)(&n,array,work,work+2*n,(int *)(work+4*n));
}

