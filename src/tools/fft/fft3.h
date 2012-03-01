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
#ifndef __FFT3_H__
#define __FFT3_H__


#include <mpi.h>
#include "ffti.h"


#define FFT3_FLAG_FILL_IN_ARRAY  0x1
#define FFT3_FLAG_DOUBLE_MEMORY  0x2
#define FFT3_FLAG_WACKY_K_SPACE  0x4


typedef struct fft3Times
{
  double total;
  double work, comm, serv;
  void (*reset)();
}
fft3_times_t;
extern const fft3_times_t *fft3_times;


void fft3_init(MPI_Comm fft_com, int dims[3], int pads[3], int bbox[6]);
void fft3_done();

fft_t* fft3_allocate_data();

void fft3_x2k(fft_t *data, int flags);
void fft3_k2x(fft_t *data, int flags);

size_t fft3_jk_index(int array_j, int array_k, int *jk, int flags); 

#endif
