#ifndef __TOP_LEVEL_FFT_H__
#define __TOP_LEVEL_FFT_H__

#ifdef FFT_DOUBLE
typedef double		type_fft;
#define MPI_TYPE_FFT    MPI_DOUBLE
#include <drfftw.h>
#include <dfftw.h>
#else
typedef float		type_fft;
#define MPI_TYPE_FFT    MPI_FLOAT
#include <srfftw.h>
#include <sfftw.h>
#endif

typedef void (*top_level_fft_op)(int id, fftw_complex *fft_source, fftw_complex *fft_output);

void top_level_fft(int in_var, int num_out_vars, const int *out_vars, top_level_fft_op worker);

#endif
