#ifndef __ROOT_GRID_FFT_H__
#define __ROOT_GRID_FFT_H__

#ifndef CONFIGURED
#error "Missing config.h include."
#endif


#include <mpi.h>
#include "../fft/ffti.h"


typedef struct RootGridFFTData
{
  int dims[nDim];
  int bbox[2*nDim];
}
root_grid_fft_t;


typedef struct RootGridFFTTune
{
  int report_times;
  int use_recv_map;
  int cache_cell_ijk;
  int num_fft_buffers;
  int num_recv_slices;
}
root_grid_fft_tune_t;
extern root_grid_fft_tune_t* root_grid_fft_tune;


typedef struct RootGridFFTTimes
{
  double prep, send, recv, work;
}
root_grid_fft_times_t;
extern const root_grid_fft_times_t* root_grid_fft_times;


typedef void (*root_grid_fft_op)(const root_grid_fft_t *config, int id, fft_t *fft_source, fft_t *fft_output, int flags);


void root_grid_fft_init(MPI_Group run_grp, MPI_Group fft_grp);
void root_grid_fft_done();
void root_grid_fft_exec(int in_var, int num_out_vars, const int *out_vars, root_grid_fft_op worker);

#endif
