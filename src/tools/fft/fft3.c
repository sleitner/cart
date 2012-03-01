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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fft3.h"
#include "ffts.h"


/*
//  Useful macros
*/
#ifndef MIN
#define MIN(x,y)        (((x) < (y)) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x,y)        (((x) > (y)) ? (x) : (y))
#endif


void fft3_reset_times_internal();
fft3_times_t fft3_times_internal = { 0.0, 0.0, 0.0, 0.0, fft3_reset_times_internal };
const fft3_times_t* fft3_times = &fft3_times_internal;


struct fft3Internal
{
  MPI_Comm com;
  int rank, size;
  int fft_dims[3];
  int arr_dims[3];
  int jrange[2];
  int krange[2];
  int jstride;
  int kstride;
  fft_t **tomp;
  int nomp;
  fft_t **work[3];
  int init;
}
fft3_internal_data;
#define d fft3_internal_data


void fft3_fail(const char *fmt, ... )
{
  char str[256];
  va_list args;

  va_start(args,fmt);
  vsnprintf(str,256,fmt,args);

  FAIL(str);

  va_end(args);
}


void fft3_init(MPI_Comm fft_com, int dims[3], int pads[3], int bbox[6])
{
  int i, j, n;

  d.com = fft_com;

  MPI_Comm_size(d.com,&d.size);
  MPI_Comm_rank(d.com,&d.rank);

  if(d.size<0 && d.rank<0 && d.rank>=d.size)
    {
      fft3_fail("fft3: Invalid MPI Communicator.");
    }

  if((dims[0] % 2) == 1)
    {
      fft3_fail("fft3: First dimension MUST be even");
    }

  for(j=0; j<3; j++)
    {
      if(dims[j] <= 0)
	{
	  fft3_fail("fft3: Invalid input dimension[%d]=%d.",j,dims[j]);
	}
      d.fft_dims[j] = dims[j];
    }

  if(dims[2] < d.size)
    {
      fft3_fail("fft3: Input K dimension is shorter than the number of MPI tasks.");
    }

  d.arr_dims[0] = d.fft_dims[0] + MAX(2,pads[0]);

  d.jstride = (d.fft_dims[1]+pads[1]+d.size-1)/d.size;
  d.kstride = (d.fft_dims[2]+pads[2]+d.size-1)/d.size;

  d.arr_dims[1] = d.jstride*d.size;
  d.arr_dims[2] = d.kstride;

  d.jrange[0] = d.jstride*d.rank;
  d.jrange[1] = d.jstride*(d.rank+1);

  d.krange[0] = d.kstride*d.rank;
  d.krange[1] = d.kstride*(d.rank+1);

  for(j=0; j<3; j++) dims[j] = d.arr_dims[j];

  bbox[0] = 0;
  bbox[1] = d.fft_dims[0];
  bbox[2] = 0;
  bbox[3] = d.fft_dims[1];
  bbox[4] = d.krange[0];
  bbox[5] = MIN(d.krange[1],d.fft_dims[2]);

#ifdef _OPENMP
  d.nomp = omp_get_max_threads();
#else
  d.nomp = 1;
#endif

  d.tomp = (fft_t**)malloc(sizeof(fft_t*)*d.nomp);
  if(d.tomp == NULL)
    {
      fft3_fail("fft3: unable to allocate temporary storage.");
    }
  for(i=0; i<d.nomp; i++)
    {
      d.tomp[i] = NEW(2*MAX(2+d.fft_dims[0],MAX(d.fft_dims[1],d.fft_dims[2])));
      if(d.tomp[i] == NULL)
	{
	  fft3_fail("fft3: unable to allocate temporary storage.");
	}
    }

  for(j=0; j<3; j++)
    {
      d.work[j] = (fft_t**)malloc(sizeof(fft_t*)*d.nomp);
      if(d.work[j] == NULL)
	{
	  FAIL("fft3: unable to allocate temporary storage.");
	}

      n = ffti_get_work_size(d.fft_dims[j]);
      for(i=0; i<d.nomp; i++)
	{
	  if(n > 0)
	    {
	      d.work[j][i] = NEW(n);
	      if(d.work[j][i] == NULL)
		{
		  fft3_fail("fft3: unable to allocate temporary storage (%d values).",n);
		}
	    }
	  else d.work[j][i] = NULL;
	  ffti_init(d.fft_dims[j],d.work[j][i]);
	}
    }

  d.init = 1;
}


void fft3_done()
{
  int i, j;

  if(d.init != 1)
    {
      fft3_fail("fft3: has not been initialized.");
    }

  for(i=0; i<d.nomp; i++) DEL(d.tomp[i]);
  free(d.tomp);

  for(j=0; j<3; j++)
    {
      for(i=0; i<d.nomp; i++)
	{
	  if(d.work[j][i] != NULL) DEL(d.work[j][i]);
	}
      free(d.work[j]);
    }

  d.init = 0;
}


fft_t* fft3_allocate_data()
{
  if(d.init != 1)
    {
      fft3_fail("fft3: has not been initialized.");
    }

  return NEW((size_t)d.arr_dims[0]*d.arr_dims[1]*d.kstride);
}


size_t fft3_jk_index(int jarr, int karr, int *jk, int flags)
{
  int jfft, kfft;

  if(jarr<0 || jarr>=d.arr_dims[1] || karr<0 || karr>=d.arr_dims[2])
    {
      return (size_t)-1;
    }

  if(d.size>1 && (flags & FFT3_FLAG_WACKY_K_SPACE))
    {
      if(flags & FFT3_FLAG_DOUBLE_MEMORY)
	{
	  jfft = d.jrange[0] + (jarr % d.jstride);
	  kfft = jarr/d.jstride + d.size*karr;
	}
      else
	{
	  jfft = d.jrange[0] + (jarr % d.jstride);
	  kfft = karr + d.kstride*(jarr/d.jstride);
	}
    }
  else
    {
      jfft = jarr;
      kfft = d.krange[0] + karr;
    }

  if(jfft>=d.fft_dims[1] || kfft>=d.fft_dims[2])
    {
      return (size_t)-1;
    }
    
  if(jk != NULL)
    {
      jk[0] = jfft;
      jk[1] = kfft;
    }
  
  return jarr + (size_t)d.arr_dims[1]*karr;
}


void fft3_reset_times_internal()
{
  fft3_times_internal.total = 0.0;
  fft3_times_internal.work = 0.0;
  fft3_times_internal.comm = 0.0;
  fft3_times_internal.serv = 0.0;
}


void fft3_internal_worker(fft_t *data, char dir, int flags);


void fft3_x2k(fft_t *data, int flags)
{
  fft3_internal_worker(data,1,flags);
}


void fft3_k2x(fft_t *data, int flags)
{
  fft3_internal_worker(data,0,flags);
}


/*
//  ----------  INTERNAL IMPLEMENTATION  ------------
*/
void fft3_internal_step_x(fft_t *data, int mode, int pad)
{
  int i, j, k, iomp;
  size_t offset;
  fft_t *ptr;
  double t = MPI_Wtime();

  /*
  //  Transform X direction
  */
  //#pragma omp parallel for default(none), private(i,j,k,ptr,iomp,offset), shared(data,d,mode,pad)
  for(k=d.krange[0]; k<MIN(d.krange[1],d.fft_dims[2]); k++)
    {
#ifdef _OPENMP
      iomp = omp_get_thread_num();
      if(iomp<0 || iomp>=d.nomp)
	{
	  fft3_fail("fft3: invalid OpenMP thread %d out of %d.",iomp,d.nomp);
	}
#else
      iomp = 0;
#endif
      ptr = d.tomp[iomp];

      offset = (size_t)d.arr_dims[1]*(k-d.krange[0]);
      for(j=0; j<d.fft_dims[1]; j++)
	{
	  if(mode & 0x1)
	    {
	      /*
	      //  First, cache X direction
	      */
	      for(i=0; i<d.fft_dims[0]; i++)
		{
		  ptr[2*i+0] = data[i+d.arr_dims[0]*(j+offset)];
		  ptr[2*i+1] = 0.0;
		}
	      
	      ffti_fc2c(d.fft_dims[0],ptr,d.work[0][iomp]);

	      /*
	      //  De-cache it back
	      */
	      for(i=0; i<=d.fft_dims[0]/2; i++)
		{
		  data[2*i+0+d.arr_dims[0]*(j+offset)] = ptr[2*i+0];
		  data[2*i+1+d.arr_dims[0]*(j+offset)] = ptr[2*i+1];
		}
	    }
	  else
	    {
	      /*
	      //  First, cache X direction
	      */
	      for(i=0; i<=d.fft_dims[0]/2; i++)
		{
		  ptr[2*i+0] = data[2*i+0+d.arr_dims[0]*(j+offset)];
		  ptr[2*i+1] = data[2*i+1+d.arr_dims[0]*(j+offset)];
		}
	      for(i=d.fft_dims[0]/2+1; i<d.fft_dims[0]; i++)
		{
		  ptr[2*i+0] =  ptr[2*(d.fft_dims[0]-i)+0];
		  ptr[2*i+1] = -ptr[2*(d.fft_dims[0]-i)+1];
		}
	      ptr[1] = ptr[d.fft_dims[0]+1] = 0.0;

	      ffti_bc2c(d.fft_dims[0],ptr,d.work[0][iomp]);

	      /*
	      //  De-cache it back
	      */
	      for(i=0; i<d.fft_dims[0]; i++)
		{
		  data[i+d.arr_dims[0]*(j+offset)] = ptr[2*i+0];
		}
	    }
	  /*
	  //  Pad data with zeros to avoid garbage in MPI buffers
	  */
	  if(pad)
	    {
	      for(i=d.fft_dims[0]+2; i<d.arr_dims[0]; i++)
		{
		  data[i+d.arr_dims[0]*(j+offset)] = 0.0;
		}
	    }
	}
      if(pad)
	{
	  for(j=d.fft_dims[1]; j<d.arr_dims[1]; j++)
	    {
	      memset(data+d.arr_dims[0]*(j+offset),0,d.arr_dims[0]*sizeof(fft_t));
	    }
	}
    }

  if(pad)
    {
      for(k=d.fft_dims[2]; k<d.krange[1]; k++)
	{
	  offset = (size_t)d.arr_dims[1]*(k-d.krange[0]);

#pragma omp parallel for default(none), private(i,j), shared(k,data,d,offset)
	  for(j=0; j<d.arr_dims[1]; j++)
	    {
	      memset(data+d.arr_dims[0]*(j+offset),0,d.arr_dims[0]*sizeof(fft_t));
	    }
	}
    }

  fft3_times_internal.work += (MPI_Wtime()-t);
}


void fft3_internal_step_y(fft_t *data, int mode)
{
  int i, j, k, iomp;
  size_t offset;
  fft_t *ptr;
  double t = MPI_Wtime();

  /*
  //  Transform Y direction
  */
#pragma omp parallel for default(none), private(i,j,k,ptr,iomp,offset), shared(data,d,mode)
  for(k=d.krange[0]; k<MIN(d.krange[1],d.fft_dims[2]); k++)
    {
#ifdef _OPENMP
      iomp = omp_get_thread_num();
      if(iomp<0 || iomp>=d.nomp)
	{
	  fft3_fail("fft3: invalid OpenMP thread %d out of %d.",iomp,d.nomp);
	}
#else
      iomp = 0;
#endif
      ptr = d.tomp[iomp];

      offset = (size_t)d.arr_dims[1]*(k-d.krange[0]);
      for(i=0; i<=d.fft_dims[0]/2; i++)
	{
	  /*
	  //  First, cache Y direction
	  */
	  for(j=0; j<d.fft_dims[1]; j++)
	    {
	      ptr[2*j+0] = data[2*i+0+d.arr_dims[0]*(j+offset)];
	      ptr[2*j+1] = data[2*i+1+d.arr_dims[0]*(j+offset)];
	    }
	     
	  if(mode&0x1)
	    {
	      ffti_fc2c(d.fft_dims[1],ptr,d.work[1][iomp]);
	    }
	  else
	    {
	      ffti_bc2c(d.fft_dims[1],ptr,d.work[1][iomp]);
	    }

	  /*
	  //  De-cache it back
	  */
	  for(j=0; j<d.fft_dims[1]; j++)
	    {
	      data[2*i+0+d.arr_dims[0]*(j+offset)] = ptr[2*j+0];
	      data[2*i+1+d.arr_dims[0]*(j+offset)] = ptr[2*j+1];
	    }
	}
    }

  fft3_times_internal.work += (MPI_Wtime()-t);
}


void fft3_internal_step_z_no_communication(fft_t *data, int mode)
{
  int i, j, k, iomp;
  size_t offset;
  fft_t *ptr;
  double tb, te;

  /*
  //  Transform Z direction (which is now weirdly packed)
  */
  tb = MPI_Wtime();
#pragma omp parallel for default(none), private(i,j,k,ptr,iomp,offset), shared(data,d,mode)
  for(j=0; j<d.fft_dims[1]; j++)
    {
#ifdef _OPENMP
      iomp = omp_get_thread_num();
      if(iomp<0 || iomp>=d.nomp)
	{
	  fft3_fail("fft3: invalid OpenMP thread %d out of %d.",iomp,d.nomp);
	}
#else
      iomp = 0;
#endif
      ptr = d.tomp[iomp];

      for(i=0; i<=d.fft_dims[0]/2; i++)
	{
	  /*
	  //  First, cache Z direction from the first array
	  */
	  for(k=0; k<d.fft_dims[2]; k++)
	    {
	      offset = d.arr_dims[0]*(j+(size_t)d.arr_dims[1]*k);
	      ptr[2*k+0] = data[2*i+0+offset];
	      ptr[2*k+1] = data[2*i+1+offset];
	    }
	     
	  if(mode&0x1)
	    {
	      ffti_fc2c(d.fft_dims[2],ptr,d.work[2][iomp]);
	    }
	  else
	    {
	      ffti_bc2c(d.fft_dims[2],ptr,d.work[2][iomp]);
	    }

	  /*
	  //  De-cache it back into the second array
	  */
	  for(k=0; k<d.fft_dims[2]; k++)
	    {
	      offset = d.arr_dims[0]*(j+(size_t)d.arr_dims[1]*k);
	      data[2*i+0+offset] = ptr[2*k+0];
	      data[2*i+1+offset] = ptr[2*k+1];
	    }
	}
    }
  te = MPI_Wtime();
  fft3_times_internal.work += (te-tb);
}


void fft3_internal_step_z_single_communication(fft_t *data, int mode)
{
  int i, j, k, l, iomp, err;
  size_t offset, block_size;
  fft_t *ptr;
  fft_t *tmp;
  double tb, te;

  /*
  //  Create storage for out-of-place data transfer (forward transform)
  */
  block_size = (size_t)d.arr_dims[0]*d.jstride*d.kstride;

  tb = MPI_Wtime();
  tmp = NEW(block_size*d.size);
  if(tmp == NULL)
    {
      fft3_fail("fft3: unable to allocate temporary storage of size %ld.",block_size*d.size);
    }

  te = MPI_Wtime();
  fft3_times_internal.serv += (te-tb);

  /*
  //  For forward transform, we transpose first
  */
  if((mode&0x1) || (mode&0x2))
    {
      tb = te;

      /*
      //  Transpose a local piece
      */
      for(l=0; l<d.size; l++)
	{
	  for(j=0; j<d.jstride; j++)
	    {
	      for(k=0; k<d.kstride; k++)
		{
		  memcpy(tmp+d.arr_dims[0]*(j+(size_t)d.jstride*k)+block_size*l,data+d.arr_dims[0]*(j+l*d.jstride+(size_t)d.arr_dims[1]*k),d.arr_dims[0]*sizeof(fft_t));
		}
	    }
	}

      /*
      //  Transpose pieces in a global array
      */
      err = MPI_Alltoall(tmp,block_size,MPI_TYPE_FFT,data,block_size,MPI_TYPE_FFT,d.com);

      if(err != MPI_SUCCESS)
	{
	  fft3_fail("fft3: MPI_Alltoall return an error code %d",err);
	}

      te = MPI_Wtime();
      fft3_times_internal.comm += (te-tb);
    }

  /*
  //  Transform Z direction
  */
  tb = te;
#pragma omp parallel for default(none), private(i,j,k,ptr,iomp,offset), shared(data,d,mode)
  for(j=d.jrange[0]; j<MIN(d.jrange[1],d.fft_dims[1]); j++)
    {
#ifdef _OPENMP
      iomp = omp_get_thread_num();
      if(iomp<0 || iomp>=d.nomp)
	{
	  fft3_fail("fft3: invalid OpenMP thread %d out of %d.",iomp,d.nomp);
	}
#else
      iomp = 0;
#endif
      ptr = d.tomp[iomp];

      for(i=0; i<=d.fft_dims[0]/2; i++)
	{
	  /*
	  //  First, cache Z direction from the first array
	  */
	  for(k=0; k<d.fft_dims[2]; k++)
	    {
	      offset = d.arr_dims[0]*(j-d.jrange[0]+(size_t)d.jstride*k);
	      ptr[2*k+0] = data[2*i+0+offset];
	      ptr[2*k+1] = data[2*i+1+offset];
	    }
	     
	  if(mode&0x1)
	    {
	      ffti_fc2c(d.fft_dims[2],ptr,d.work[2][iomp]);
	    }
	  else
	    {
	      ffti_bc2c(d.fft_dims[2],ptr,d.work[2][iomp]);
	    }

	  /*
	  //  De-cache it back into the second array
	  */
	  for(k=0; k<d.fft_dims[2]; k++)
	    {
	      offset = d.arr_dims[0]*(j-d.jrange[0]+(size_t)d.jstride*k);
	      data[2*i+0+offset] = ptr[2*k+0];
	      data[2*i+1+offset] = ptr[2*k+1];
	    }
	}
    }
  te = MPI_Wtime();
  fft3_times_internal.work += (te-tb);

  /*
  //  For backward transform, we transpose last
  */
  if(!(mode&0x1) || (mode&0x2))
    {
      tb = te;

      /*
      //  Transpose pieces in a global array
      */
      err = MPI_Alltoall(data,block_size,MPI_TYPE_FFT,tmp,block_size,MPI_TYPE_FFT,d.com);

      if(err != MPI_SUCCESS)
	{
	  fft3_fail("fft3: MPI_Alltoall return an error code %d",err);
	}

      /*
      //  Transpose a local piece
      */
      for(l=0; l<d.size; l++)
	{
	  for(j=0; j<d.jstride; j++)
	    {
	      for(k=0; k<d.kstride; k++)
		{
		  memcpy(data+d.arr_dims[0]*(j+l*d.jstride+(size_t)d.arr_dims[1]*k),tmp+d.arr_dims[0]*(j+(size_t)d.jstride*k)+block_size*l,d.arr_dims[0]*sizeof(fft_t));
		}
	    }
	}

      te = MPI_Wtime();
      fft3_times_internal.comm += (te-tb);
    }

  tb = te;
  DEL(tmp);
  te = MPI_Wtime();
  fft3_times_internal.serv += (te-tb);
}


void fft3_internal_step_z_multiple_communication(fft_t *data, int mode)
{
  int i, j, k, l, iomp, err;
  size_t offset, block_size;
  fft_t *ptr;
  fft_t *tmp1, *tmp2;
  fft_t *arr1, *arr2;
  double tb, te;

  /*
  //  Create storage for out-of-place data transfer
  */
  block_size = (size_t)d.arr_dims[0]*d.kstride;

  tb = MPI_Wtime();
  tmp1 = NEW(block_size*d.size);
  if(tmp1 == NULL)
    {
      fft3_fail("fft3: unable to allocate temporary storage of size %ld.",block_size*d.size);
    }
  tmp2 = NEW(block_size*d.size);
  if(tmp2 == NULL)
    {
      DEL(tmp1);
      fft3_fail("fft3: unable to allocate temporary storage of size %ld.",block_size*d.size);
    }

  /*
  //  Set-up array pointers
  */
  if(mode&0x2)
    {
      arr1 = tmp2;
      arr2 = tmp2;
    }
  else if(mode&0x1)
    {
      arr1 = tmp2;
      arr2 = tmp1;
    }
  else
    {
      arr1 = tmp1;
      arr2 = tmp2;
    }
  te = MPI_Wtime();
  fft3_times_internal.serv += (te-tb);

  /*
  //  Transform Z direction
  */
  for(j=0; j<d.jstride; j++)
    {
      tb = te;
      /*
      //  Buffer
      */
      for(l=0; l<d.size; l++)
	{
	  for(k=0; k<d.kstride; k++)
	    {
	      memcpy(tmp1+(size_t)d.arr_dims[0]*k+block_size*l,data+d.arr_dims[0]*(j+l*d.jstride+(size_t)d.arr_dims[1]*k),d.arr_dims[0]*sizeof(fft_t));
	    }
	}
      te = MPI_Wtime();
      fft3_times_internal.serv += (te-tb);

      /*
      //  For forward transform, we transpose first
      */
      if((mode&0x1) || (mode&0x2))
	{
	  tb = te;
	  err = MPI_Alltoall(tmp1,block_size,MPI_TYPE_FFT,tmp2,block_size,MPI_TYPE_FFT,d.com);
	  if(err != MPI_SUCCESS)
	    {
	      fft3_fail("fft3: MPI_Alltoall return an error code %d",err);
	    }
	  te = MPI_Wtime();
	  fft3_times_internal.comm += (te-tb);
	}

      tb = te;
#pragma omp parallel for default(none), private(i,k,ptr,iomp,offset), shared(j,data,d,arr1,arr2,mode)
      for(i=0; i<=d.fft_dims[0]/2; i++)
	{
#ifdef _OPENMP
	  iomp = omp_get_thread_num();
	  if(iomp<0 || iomp>=d.nomp)
	    {
	      fft3_fail("fft3: invalid OpenMP thread %d out of %d.",iomp,d.nomp);
	    }
#else
	  iomp = 0;
#endif
	  ptr = d.tomp[iomp];

	  /*
	  //  First, cache Z direction
	  */
	  for(k=0; k<d.fft_dims[2]; k++)
	    {
	      offset = (size_t)d.arr_dims[0]*k;
	      ptr[2*k+0] = arr1[2*i+0+offset];
	      ptr[2*k+1] = arr1[2*i+1+offset];
	    }
	     
	  if(mode&0x1)
	    {
	      ffti_fc2c(d.fft_dims[2],ptr,d.work[2][iomp]);
	    }
	  else
	    {
	      ffti_bc2c(d.fft_dims[2],ptr,d.work[2][iomp]);
	    }

	  /*
	  //  De-cache it back into the data array in proper order
	  */
	  for(k=0; k<d.fft_dims[2]; k++)
	    {
	      offset = (size_t)d.arr_dims[0]*k;
	      arr2[2*i+0+offset] = ptr[2*k+0];
	      arr2[2*i+1+offset] = ptr[2*k+1];
	    }
	}
      te = MPI_Wtime();
      fft3_times_internal.work += (te-tb);

      /*
      //  For backward transform, we transpose last
      */
      if(!(mode&0x1) || (mode&0x2))
	{
	  tb = te;
	  err = MPI_Alltoall(tmp2,block_size,MPI_TYPE_FFT,tmp1,block_size,MPI_TYPE_FFT,d.com);
	  if(err != MPI_SUCCESS)
	    {
	      fft3_fail("fft3: MPI_Alltoall return an error code %d.",err);
	    }
	  te = MPI_Wtime();
	  fft3_times_internal.comm += (te-tb);
	}

      /*
      //  Un-buffer
      */
      tb = te;
      for(l=0; l<d.size; l++)
	{
	  for(k=0; k<d.kstride; k++)
	    {
	      memcpy(data+(size_t)d.arr_dims[0]*(j+l*d.jstride+(size_t)d.arr_dims[1]*k),tmp1+(size_t)d.arr_dims[0]*k+block_size*l,d.arr_dims[0]*sizeof(fft_t));
	    }
	}
      te = MPI_Wtime();
      fft3_times_internal.serv += (te-tb);
    }

  tb = te;
  DEL(tmp1);
  DEL(tmp2);
  te = MPI_Wtime();
  fft3_times_internal.serv += (te-tb);
}


void fft3_internal_worker(fft_t *data, char dir, int flags)
{
  int mode = 0;

  if(d.init != 1)
    {
      fft3_fail("fft3: has not been initialized.");
    }

  if(data == NULL)
    {
      fft3_fail("NULL data at the entry point."); 
    } 

  if((dir&0x1) != 0) mode += 0x1;
  if((flags&FFT3_FLAG_WACKY_K_SPACE) == 0) mode += 0x2;

  if(dir)
    {
      fft3_internal_step_x(data,mode,((flags&FFT3_FLAG_FILL_IN_ARRAY)==0) ? 0 : 1);
      fft3_internal_step_y(data,mode);
      if(d.size == 1)
	{
	  fft3_internal_step_z_no_communication(data,mode);
	}
      else if((flags&FFT3_FLAG_DOUBLE_MEMORY) == 0)
	{
	  fft3_internal_step_z_multiple_communication(data,mode);
	}
      else
	{
	  fft3_internal_step_z_single_communication(data,mode);
	}
    }
  else
    {
      if(d.size == 1)
	{
	  fft3_internal_step_z_no_communication(data,mode);
	}
      else if((flags&FFT3_FLAG_DOUBLE_MEMORY) == 0)
	{
	  fft3_internal_step_z_multiple_communication(data,mode);
	}
      else
	{
	  fft3_internal_step_z_single_communication(data,mode);
	}
      fft3_internal_step_y(data,mode);
      fft3_internal_step_x(data,mode,0);
    }

  fft3_times_internal.total = fft3_times_internal.work + fft3_times_internal.comm + fft3_times_internal.serv;
}


