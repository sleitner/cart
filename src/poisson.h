#ifndef __POISSON_H__
#define __POISSON_H__

#ifdef FFT_DOUBLE
        typedef double		type_fft;
#define MPI_TYPE_FFT            MPI_DOUBLE
#else
        typedef float		type_fft;
#define MPI_TYPE_FFT            MPI_FLOAT
#endif

void poisson( type_fft *density_grid, type_fft *potential_grid );

#endif
