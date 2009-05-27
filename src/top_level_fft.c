#include "defs.h"

#include <math.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "sfc.h"
#include "timing.h"
#include "top_level_fft.h"

#if defined(GRAVITY) || defined(RADIATIVE_TRANSFER) 

/* compute FFT on min_level (WILL BE REWRITTEN!) */
void top_level_fft(int in_var, int num_out_vars, const int *out_vars, top_level_fft_op worker)
{
  int coords[3];
  int index;
  int p, i, var;
  type_fft *data, *buffer;
  fftw_complex *fft_source, *fft_output;
  fftwnd_plan forward, backward;
  MPI_Status status;

  cart_assert(in_var>=0 && in_var<num_vars);
  cart_assert(num_out_vars > 0);
  for(var=0; var<num_out_vars; var++)
    {
      cart_assert(out_vars[var]>=0 && out_vars[var]<num_vars);
    }

  start_time( FFT_TIMER );

  buffer = cart_alloc(type_fft, num_root_cells );

  /* gather all root sources to master node */
  if ( local_proc_id == MASTER_NODE )
    {
      start_time( WORK_TIMER );

      data = cart_alloc(type_fft, num_root_cells );

#pragma omp parallel for default(none), private(i,coords,index), shared(num_cells_per_level,data,cell_vars,in_var)
      for ( i = 0; i < num_cells_per_level[min_level]; i++ )
	{
	  sfc_coords( i, coords );
	  index = num_grid*(num_grid*coords[0] + coords[1] ) + coords[2];
	  data[index] = cell_var(i,in_var);
	}

      end_time( WORK_TIMER );

      start_time( COMMUNICATION_TIMER );

      for ( p = 1; p < num_procs; p++ )
	{
	  MPI_Recv( buffer, proc_sfc_index[p+1]-proc_sfc_index[p], MPI_TYPE_FFT, p, 0, MPI_COMM_WORLD, &status );

#pragma omp parallel for default(none), private(i,coords,index), shared(proc_sfc_index,p,data,buffer)
	  for ( i = 0; i < proc_sfc_index[p+1]-proc_sfc_index[p]; i++ )
	    {
	      sfc_coords(i+proc_sfc_index[p],coords);
	      index = num_grid*(num_grid*coords[0] + coords[1] ) + coords[2];
	      data[index] = buffer[i];
	    }
	}

      end_time( COMMUNICATION_TIMER );

      /*
      //  In principle, memory can be saved by overlaying this array with data[] and doing in-place FFT.
      */
      fft_source = cart_alloc(fftw_complex, num_grid*num_grid*(num_grid/2+1) );

      forward = rfftw3d_create_plan( num_grid, num_grid, num_grid, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE );
      rfftwnd_one_real_to_complex( forward, data, fft_source );
      rfftwnd_destroy_plan(forward);

      if(num_out_vars > 1)
	{
	  fft_output = cart_alloc(fftw_complex, num_grid*num_grid*(num_grid/2+1) );
	}
      else
	{
	  fft_output = fft_source;
	}

      backward = rfftw3d_create_plan( num_grid, num_grid, num_grid, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_OUT_OF_PLACE );

      /* compute FFT */
      for(var=0; var<num_out_vars; var++)
	{
				     
	  start_time( WORK_TIMER );
	  (*worker)(var,fft_source,fft_output);
	  end_time( WORK_TIMER );

	  /* perform inverse transform to get result */
	  rfftwnd_one_complex_to_real( backward, fft_output, data);

	  start_time( COMMUNICATION_TIMER );
	  for ( p = 1; p < num_procs; p++ )
	    {
#pragma omp parallel for default(none), private(i,coords,index), shared(proc_sfc_index,p,data,buffer)
	      for ( i = 0; i < proc_sfc_index[p+1]-proc_sfc_index[p]; i++ )
		{
		  sfc_coords(i+proc_sfc_index[p],coords);
		  index = num_grid*(num_grid*coords[0] + coords[1] ) + coords[2];
		  buffer[i] = data[index];
		}

	      MPI_Send( buffer, proc_sfc_index[p+1]-proc_sfc_index[p], MPI_TYPE_FFT, p, 1, MPI_COMM_WORLD );
	    }
	  end_time( COMMUNICATION_TIMER );

	  start_time( WORK_TIMER );
#pragma omp parallel for default(none), private(i,coords,index), shared(num_cells_per_level,cell_vars,data,var,out_vars)
	  for ( i = 0; i < num_cells_per_level[min_level]; i++ )
	    {
	      sfc_coords( i, coords );
	      index = num_grid*(num_grid*coords[0] + coords[1] ) + coords[2];
	      
	      cell_var(i,out_vars[var]) = data[index];
	    }
	  end_time( WORK_TIMER );

	}

      rfftwnd_destroy_plan(backward);

      if(num_out_vars > 1)
	{
	  cart_free(fft_output);
	}
      cart_free(fft_source);
      cart_free(data);
    }
  else
    {
      start_time( COMMUNICATION_TIMER );
      
#pragma omp parallel for default(none), private(i), shared(num_cells_per_level,cell_vars,buffer,in_var)
      for ( i = 0; i < num_cells_per_level[min_level]; i++ )
	{
	  buffer[i] = cell_var(i,in_var);
	}
		
      MPI_Send( buffer, num_cells_per_level[min_level], MPI_TYPE_FFT, MASTER_NODE, 0, MPI_COMM_WORLD );

      for(var=0; var<num_out_vars; var++)
	{
	  MPI_Recv( buffer, num_cells_per_level[min_level], MPI_TYPE_FFT, MASTER_NODE, 1, MPI_COMM_WORLD, &status );

#pragma omp parallel for default(none), private(i), shared(num_cells_per_level,cell_vars,buffer,var,out_vars)
	  for ( i = 0; i < num_cells_per_level[min_level]; i++ )
	    {
	      cell_var(i,out_vars[var]) = buffer[i];
	    }
	}
      end_time( COMMUNICATION_TIMER );
    }
     
  cart_free(buffer);

  /* update cell buffer */
  start_time( FFT_UPDATE_TIMER );
  update_buffer_level( min_level, out_vars, num_out_vars );
  end_time( FFT_UPDATE_TIMER );

  end_time( FFT_TIMER );
}

#endif /* defined(GRAVITY) || defined(RADIATIVE_TRANSFER) */
