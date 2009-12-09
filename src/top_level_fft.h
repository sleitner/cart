#ifndef __TOP_LEVEL_FFT_H__
#define __TOP_LEVEL_FFT_H__

#include "defs.h"

#ifdef FFTW3
	#include "fftw3.h"

	#define c_re(c) ((c)[0])
	#define c_im(c) ((c)[1])
#else /* default to fftw 2.1.5 */
	#ifdef FFT_OPENMP
		#ifdef FFT_DOUBLE
			#include <drfftw_threads.h>
			#include <dfftw_threads.h>
		#else
			#include <srfftw_threads.h>
			#include <sfftw_threads.h>
		#endif
	#else
		#ifdef FFT_DOUBLE
			#include <drfftw.h>
			#include <dfftw.h>
		#else
			#include <srfftw.h>
			#include <sfftw.h>
		#endif
	#endif
#endif

#ifdef FFT_DOUBLE
typedef fftw_complex	type_fft_complex;
typedef double			type_fft;
#define MPI_TYPE_FFT	MPI_DOUBLE
#else
#ifdef FFTW3
typedef fftwf_complex	type_fft_complex;
#else
typedef fftw_complex	type_fft_complex;
#endif
typedef float			type_fft;
#define MPI_TYPE_FFT    MPI_FLOAT
#endif

void init_fft();
typedef void (*top_level_fft_op)(int id, type_fft_complex *fft_source, type_fft_complex *fft_output);
void top_level_fft(int in_var, int num_out_vars, const int *out_vars, top_level_fft_op worker);

#endif
