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
#ifndef __FFTI_H__
#define __FFTI_H__


#ifdef FFT_DOUBLE
typedef double		fft_t;
#define MPI_TYPE_FFT	MPI_DOUBLE
#else
typedef float		fft_t;
#define MPI_TYPE_FFT    MPI_FLOAT
#endif


/*
//  If the implementation needs some internal memory, it should
//  return the number of fft_t values it needs; if the memory
//  is not needed, the implementation MUST return 0.
*/
int ffti_get_work_size(int n);

/*
//
//  ffti_init will be called multiple times (3*num_omp_threads). 
*/
void ffti_init(int n, fft_t *work);

/*
//  1D complex-to-complex FFTs (forward and backward)
*/
void ffti_fc2c(int n, fft_t *array, fft_t *work);
void ffti_bc2c(int n, fft_t *array, fft_t *work);

#endif
