#include "config.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "auxiliary.h"
#include "cell_buffer.h"
#include "iterators.h"
#include "parallel.h"
#include "root_grid_fft.h"
#include "sfc.h"
#include "timing.h"
#include "tree.h"

#include "../tools/fft/fft3.h"


#define __CHECK


#ifndef ROOT_GRID_FFT_IJK_TYPE
#define ROOT_GRID_FFT_IJK_TYPE int
#endif

typedef ROOT_GRID_FFT_IJK_TYPE ijk_t;


root_grid_fft_t root_grid_fft_internal_data;
#define d root_grid_fft_internal_data

root_grid_fft_tune_t root_grid_fft_internal_tune = { 0, 1, 1, 1, 4 };
root_grid_fft_tune_t* root_grid_fft_tune = &root_grid_fft_internal_tune;


root_grid_fft_times_t root_grid_fft_times_internal = { 0.0, 0.0, 0.0, 0.0 };
const root_grid_fft_times_t* root_grid_fft_times = &root_grid_fft_times_internal;


typedef struct RootGridFFTCache
{
  int *tmp;
  MPI_Request *reqs;
  int num_tasks, *tasks;
  int *sizes;
}
cache_t;


typedef struct RootGridFFTValue
{
  float value;
  ijk_t index;
}
value_t;


typedef struct RootGridFFTMap
{
  ijk_t **lines;
}
map_t;


typedef struct RootGridFFTCellIJK
{
  int ijk[nDim];
}
cell_ijk_t;


struct RootGridFFTConfig
{
  MPI_Comm com;
  int rank, size;
  int is_run, is_fft, is_run_head, is_fft_head;
  int run_size, fft_size;
  int *state;
  map_t *map;
  cell_ijk_t *cell;
}
root_grid_fft_internal_config;
#define c root_grid_fft_internal_config


void root_grid_fft_init(MPI_Group run_grp, MPI_Group fft_grp)
{
  int pads[] = { 1, 0, 0 };
  int i, run_rank, fft_rank;
  MPI_Group all_grp;
  int *state;
  size_t l;

  /*
  //  Create a communicator for FFT data transfer.
  //  FFT group goes first, so that we know that rank=0 task is an FFT one.
  */
  MPI_Group_union(fft_grp,run_grp,&all_grp);
  MPI_Comm_create(mpi.comm.world,all_grp,&c.com);
  MPI_Group_free(&all_grp);

  /*
  //  Sizes and ranks can only be safely querued from a group, 
  //  not a communicator!!!
  */
  MPI_Group_rank(run_grp,&run_rank);
  c.is_run = (run_rank != MPI_UNDEFINED);

  MPI_Group_rank(fft_grp,&fft_rank);
  c.is_fft = (fft_rank != MPI_UNDEFINED);

  if(!c.is_run && !c.is_fft) return;

  MPI_Group_size(run_grp,&c.run_size);
  MPI_Group_size(fft_grp,&c.fft_size);

  MPI_Comm_size(c.com,&c.size);
  MPI_Comm_rank(c.com,&c.rank);
  cart_assert(c.rank>=0 && c.rank<c.size);

  if(fft_rank == 0)
    {
      cart_assert(c.rank == 0);
    }

  c.is_run_head = (run_rank == 0);
  c.is_fft_head = (fft_rank == 0);
  
  /*
  //  Create a state
  */
  state = cart_alloc(int,c.size);
  c.state = cart_alloc(int,c.size);
  for(i=0; i<c.size; i++)
    {
      state[i] = (c.is_fft?1:0) + (c.is_run?2:0);
    }
  
  MPI_Alltoall(state,1,MPI_INT,c.state,1,MPI_INT,c.com);

  cart_free(state);

  /*
  //  Make sure all fft tasks are upfront.
  */
  for(i=0; i<c.fft_size; i++) cart_assert(c.state[i] & 1);

  /*
  //  Init FFT
  */
  if(c.is_fft)
    {
      d.dims[0] = d.dims[1] = d.dims[2] = num_grid;
      fft3_init(mpi.comm.fft,d.dims,pads,d.bbox);
      cart_assert(c.rank*d.dims[2] == d.bbox[4]);
    }

  /*
  //  All tasks must know the data dimensions; rank=0 is an FFT task.
  */
  MPI_Bcast(d.dims,3,MPI_INT,0,c.com);

  /*
  //  Check that typedef ijk_t is large enough to address this space
  */
  l = (size_t)d.dims[0]*d.dims[1]*d.dims[2];
  i = 0;
  while(l > 0)
    {
      l = (l >> 8);
      i++;
    }
  if(sizeof(ijk_t) < i)
    {
      cart_error("root_grid_fft: type ijk_t needs to be upgraded to a wider int type, of at least %d bytes. Please define ROOT_GRID_FFT_IJK_TYPE with a wide enough integer type.",i);
    }
}


void root_grid_fft_done()
{
  MPI_Comm_free(&c.com);
  cart_free(c.state);
}


/*
//  Domain mapping interface.
*/
void root_grid_fft_get_cell_ijk(int cell, int ijk[nDim]);
int root_grid_fft_get_task(int ijk[nDim]);
int root_grid_fft_get_task_idx(int ijk[nDim], const cache_t *cache);
size_t root_grid_fft_get_page_index(int ijk[nDim], int task, int num_slices, int page);
void root_grid_fft_print_map(const char *title, cache_t *cache);


void root_grid_fft_internal_set_cache(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells);
void root_grid_fft_internal_send_data(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells, int in_var, fft_t *data);
void root_grid_fft_internal_recv_data(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells, int out_var, fft_t *data);


/* 
//  Compute FFT on the root grid
*/
void root_grid_fft_exec(int in_var, int num_out_vars, const int *out_vars, root_grid_fft_op worker)
{
  int i, var;
  int num_level_cells, *level_cells;
  cache_t run2fft, fft2run;
  fft_t *data, *dout;
  double tb, ttl[4], ttg[4];
  int fft_flags = FFT3_FLAG_WACKY_K_SPACE;


  if(!c.is_run && !c.is_fft) return;

  cart_assert(in_var>=0 && in_var<num_vars);
  cart_assert(num_out_vars > 0);
  for(var=0; var<num_out_vars; var++)
    {
      cart_assert(out_vars[var]>=0 && out_vars[var]<num_vars);
    }

  start_time( FFT_TIMER );

  if(c.is_run)
    {
      select_level(min_level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);

      if(root_grid_fft_tune->cache_cell_ijk)
	{
	  c.cell = cart_alloc(cell_ijk_t,num_level_cells);
	}
      else c.cell = NULL;
    }
  else
    {
      num_level_cells = 0;
      level_cells = NULL;
      c.cell = NULL;
    }

  if(c.is_fft)
    {
      data = fft3_allocate_data();
      if(num_out_vars > 1)
	{
	  /*
	  //  Out of place transform
	  */
	  dout = fft3_allocate_data();
	}
      else dout = data;
      if(root_grid_fft_tune->use_recv_map)
	{
	  c.map = cart_alloc(map_t,1);
	}
      else c.map = NULL;
    }
  else
    {
      data = NULL;
      dout = NULL;
      c.map = NULL;
    }

  /*
  //  Pass our local data from in_var to the respective sections of fft nodes.
  */
  start_time( COMMUNICATION_TIMER );
  start_time( FFT_COMMUNICATION_TIMER );

  tb = MPI_Wtime();
  root_grid_fft_internal_set_cache(&run2fft,&fft2run,num_level_cells,level_cells);
  root_grid_fft_times_internal.prep = MPI_Wtime() - tb;

  tb = MPI_Wtime();
  root_grid_fft_internal_send_data(&run2fft,&fft2run,num_level_cells,level_cells,in_var,data);
  root_grid_fft_times_internal.send = MPI_Wtime() - tb;

  end_time( FFT_COMMUNICATION_TIMER );
  end_time( COMMUNICATION_TIMER );

  /*
  //  Do the actual transform: 
  //    1. transform input forward
  //    2. then for each output field:
  //       3. call the worker
  //       4. transform result backward
  //       5. send the result back to run tasks
  //
  //  TO DO: is there any work run tasks can do while fft tasks work?
  */
  if(c.is_fft)
    {
      start_time( WORK_TIMER );

      tb = MPI_Wtime();
      fft3_x2k(data,fft_flags);
      root_grid_fft_times_internal.work = MPI_Wtime() - tb;

      end_time( WORK_TIMER );
    }

  root_grid_fft_times_internal.recv = 0.0;
  for(var=0; var<num_out_vars; var++)
    {
      if(c.is_fft)
	{
	  start_time( WORK_TIMER );

	  tb = MPI_Wtime();
	  worker(&d,var,data,dout,fft_flags);
	  fft3_k2x(dout,fft_flags);
	  root_grid_fft_times_internal.work += (MPI_Wtime()-tb);

	  end_time( WORK_TIMER );
	}

      start_time( COMMUNICATION_TIMER );
      start_time( FFT_COMMUNICATION_TIMER );

      tb = MPI_Wtime();
      root_grid_fft_internal_recv_data(&run2fft,&fft2run,num_level_cells,level_cells,out_vars[var],dout);
      root_grid_fft_times_internal.recv += (MPI_Wtime()-tb);

      end_time( FFT_COMMUNICATION_TIMER );
      end_time( COMMUNICATION_TIMER );
    }

  cart_free(level_cells);
  cart_free(c.cell);
  cart_free(data);
  if(num_out_vars > 1) cart_free(dout);

  if(c.map != NULL)
    {
      for(i=0; i<fft2run.num_tasks; i++)
	{
	  cart_free(c.map->lines[i]);
	}
      cart_free(c.map->lines);
      cart_free(c.map);
    }
 
  cart_free(run2fft.tasks);
  cart_free(run2fft.sizes);
  cart_free(fft2run.tasks);
  cart_free(fft2run.sizes);

  /*
  //  Update cell buffer
  */
  start_time( FFT_UPDATE_TIMER );
  update_buffer_level(min_level,out_vars,num_out_vars);
  end_time( FFT_UPDATE_TIMER );

  if(c.is_fft)
    {
      ttl[0] = root_grid_fft_times->prep;
      ttl[1] = root_grid_fft_times->send;
      ttl[2] = root_grid_fft_times->recv;
      ttl[3] = root_grid_fft_times->work;

      switch(root_grid_fft_tune->report_times)
	{
	case 1:
	  {
	    MPI_Allreduce(ttl,ttg,4,MPI_DOUBLE,MPI_SUM,mpi.comm.fft);
	    
	    for(i=0; i<4; i++) ttl[i] = ttg[i]/c.fft_size;
	    i = c.is_fft_head;
	    break;
	  }
	case 2:
	  {
	    i = 1;
	    break;
	  }
	default:
	  {
	    i = 0;
	  }
	}

      if(i)
	{
	  cart_debug("FFT times: prep=%le, send=%le, recv=%le, work=%le",ttl[0],ttl[1],ttl[2],ttl[3]);
	}
    }

  end_time( FFT_TIMER );
}


void root_grid_fft_internal_set_cache(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells)
{
#ifdef _OPENMP
  int num_omp = omp_get_max_threads();
#else
  int num_omp = 1;
#endif 
  int iomp, *tmp;
  int i, j, k, ijk[nDim];

  /*
  //  Establish a maping of who is talking to whom.
  //  The domain decomposition of fft tasks is already
  //  abstracted.
  */
  if(c.is_fft)
    {
      /*
      //  We receive the size list and cache it just like run tasks do.
      //  convenient to cache recipient tasks, and we save memory too.
      */
      fft2run->tmp = cart_alloc(int,c.run_size);
      fft2run->reqs = cart_alloc(MPI_Request,c.run_size);
      for(i=j=0; i<c.size; i++) if(c.state[i] & 2)
	{
	  MPI_Irecv(fft2run->tmp+j,1,MPI_INT,i,10,c.com,fft2run->reqs+j);
	  j++;
	}
      cart_assert(j == c.run_size);
    }

  if(c.is_run)
    {
      run2fft->tmp = cart_alloc(int,c.fft_size);
      tmp = cart_alloc(int,c.fft_size*num_omp);
      memset(tmp,0,c.fft_size*num_omp*sizeof(int));

#pragma omp parallel for default(none), private(i,k,iomp,ijk), shared(num_level_cells,level_cells,tmp,run2fft,num_omp,root_grid_fft_internal_config)
      for(i=0; i<num_level_cells; i++)
	{
#ifdef _OPENMP
	  iomp = omp_get_thread_num();
	  cart_assert(iomp>=0 && iomp<num_omp);
#else
	  iomp = 0;
#endif
	  root_grid_fft_get_cell_ijk(level_cells[i],ijk);
	  k = root_grid_fft_get_task(ijk);
	  cart_assert(k>=0 && k<c.fft_size);

	  tmp[iomp+k*num_omp]++;

	  if(c.cell != NULL)
	    {
	      for(k=0; k<nDim; k++) c.cell[i].ijk[k] = ijk[k];
	    }
	}

      for(k=0; k<c.fft_size; k++)
	{
	  run2fft->tmp[k] = tmp[0+k*num_omp];
	  for(j=1; j<num_omp; j++)
	    {
	      run2fft->tmp[k] += tmp[j+k*num_omp];
	    }
	}

      cart_free(tmp);

      run2fft->reqs = cart_alloc(MPI_Request,c.fft_size);
      for(i=j=0; i<c.size; i++) if(c.state[i] & 1)
	{
	  MPI_Isend(run2fft->tmp+j,1,MPI_INT,i,10,c.com,run2fft->reqs+j);
	  j++;
	}
      cart_assert(j == c.fft_size);

      /*
      //  Now run2fft contains the sizes of overlapping regions for all fft 
      //  tasks for each run task and vice versa for fft2run. Now it is
      //  convenient to cache recipient tasks, and we save memory too.
      //
      //  We are not yet in a hurry to complete sends on run tasks, we can 
      //  overlay some work with them.
      */
      for(run2fft->num_tasks=j=0; j<c.fft_size; j++) if(run2fft->tmp[j] > 0) run2fft->num_tasks++;
      cart_assert(run2fft->num_tasks > 0);
      
      run2fft->tasks = cart_alloc(int,run2fft->num_tasks);
      run2fft->sizes = cart_alloc(int,run2fft->num_tasks);
      for(i=j=k=0; i<c.size; i++) if(c.state[i] & 1)
	{
	  if(run2fft->tmp[j] > 0)
	    {
	      run2fft->tasks[k] = i;
	      run2fft->sizes[k] = run2fft->tmp[j];
	      k++;
	    }
	  j++;
	}
      cart_assert(j == c.fft_size);
      cart_assert(k == run2fft->num_tasks);

     }
  else
    {
      run2fft->tasks = NULL;
      run2fft->sizes = NULL;
    }

  if(c.is_fft)
    {
      /*
      //  We must complete our receives to ensure we have the correct
      //  size list array.
      */
      MPI_Waitall(c.run_size,fft2run->reqs,MPI_STATUSES_IGNORE);
      cart_free(fft2run->reqs);

      /*
      //  Cache it now.
      */
      for(fft2run->num_tasks=j=0; j<c.run_size; j++) if(fft2run->tmp[j] > 0) fft2run->num_tasks++;
      /*
      //  fft2run->num_tasks can be zero, so we don't assert it
      //  (example: dims[2]=64, fft_size=12, task=11). 
      */

      fft2run->tasks = cart_alloc(int,fft2run->num_tasks);
      fft2run->sizes = cart_alloc(int,fft2run->num_tasks);
      for(i=j=k=0; i<c.size; i++) if(c.state[i] & 2)
	{
	  if(fft2run->tmp[j] > 0)
	    {
	      fft2run->tasks[k] = i;
	      fft2run->sizes[k] = fft2run->tmp[j];
	      k++;
	    }
	  j++;
	}
      cart_assert(j == c.run_size);
      cart_assert(k == fft2run->num_tasks);

      /*
      //  tmps are no longer needed - we cached all communication relations.
      */
      cart_free(fft2run->tmp);

#ifdef TEST
      root_grid_fft_print_map("fft->run:",fft2run);
#endif
    }
  else
    {
      fft2run->tasks = NULL;
      fft2run->sizes = NULL;
    }

  if(c.is_run)
    {
      /*
      //  Complete the first wave of sends now - this is the latest we can do so.
      //  If we wanted to save even more memory, we would do that before we
      //  created data buffers.
      */
      MPI_Waitall(c.fft_size,run2fft->reqs,MPI_STATUSES_IGNORE);
      cart_free(run2fft->reqs);

#ifdef TEST
      root_grid_fft_print_map("run->fft:",run2fft);
#endif
    }
}


void root_grid_fft_internal_send_data(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells, int in_var, fft_t *data)
{
  int i, j, k, idx, ijk[nDim];
  value_t **sbuffers, **rbuffers, *ptr;
  int *counts, size;
  MPI_Status status;
  int l, lidx, kb, num_fft_buffers = MAX(1,root_grid_fft_tune->num_fft_buffers);
  int sender[num_fft_buffers];

  /*
  //  Since we send data from a random distribution, a recipient 
  //  fft task has no way of knowing from where the value came, hence
  //  we need to send indicies as well as values. We send the array index 
  //  of type ijk_t rather than 3 int indicies to reduce the amount of data
  //  sent; in most cases ijk_t can be just an int. We also may consider
  //  saving the source run task for each value, so that we know where
  //  to send it back. Alternatively, we will save memory by not storing
  //  the source task, but then we will need to send all the data from each
  //  fft task to all connected run tasks. Optional map is used for that.
  //
  //  We can either transfer to one fft node at a time, requiring only one
  //  additional buffer array. Alternatively, we can allocate buffers for
  //  all fft tasks. That will require only one pass over root cells, but 
  //  will usenum_ more memory (4 fields for all root cells). It seems to me
  //  that the first approach will cause major non-synchronization of data
  //  transfers, so at present the second approach is adopted. 
  //
  //  TO DO: It would be nice in the future to implement sending by pages, 
  //  to minimize memory use. That is not easy, though, as sends and recvs
  //  need to be overlayed to ensure page integrity.
  */
  if(c.is_run)
    {
      sbuffers = cart_alloc(value_t*,run2fft->num_tasks);
      for(k=0; k<run2fft->num_tasks; k++)
	{
	  sbuffers[k] = cart_alloc(value_t,run2fft->sizes[k]);
	}

      counts = cart_alloc(int,run2fft->num_tasks);
      memset(counts,0,run2fft->num_tasks*sizeof(int));

      /*
      //  Serial loop. It is not easy to make it OpenMP parallel.
      */
      for(i=0; i<num_level_cells; i++)
	{
	  if(c.cell == NULL)
	    {
	      root_grid_fft_get_cell_ijk(level_cells[i],ijk);
	      k = root_grid_fft_get_task_idx(ijk,run2fft);
	    }
	  else
	    {
	      k = root_grid_fft_get_task_idx(c.cell[i].ijk,run2fft);
	    }
	  cart_assert(k>=0 && k<run2fft->num_tasks);

	  sbuffers[k][counts[k]].value = cell_var(level_cells[i],in_var);
	  /*
	  //  Pack the index
	  */
	  if(c.cell == NULL)
	    {
	      sbuffers[k][counts[k]].index = root_grid_fft_get_page_index(ijk,run2fft->tasks[k],1,0);
	    }
	  else
	    {
	      sbuffers[k][counts[k]].index = root_grid_fft_get_page_index(c.cell[i].ijk,run2fft->tasks[k],1,0);
	    }

	  counts[k]++;
	}

      for(k=0; k<run2fft->num_tasks; k++)
	{
	  cart_assert(counts[k] == run2fft->sizes[k]);
	}

      cart_free(counts);

      /*
      //  tmps are no longer needed - we cached all communication relations.
      */
      cart_free(run2fft->tmp);

      run2fft->reqs = cart_alloc(MPI_Request,run2fft->num_tasks);
      for(k=0; k<run2fft->num_tasks; k++)
	{
	  MPI_Isend(sbuffers[k],run2fft->sizes[k]*sizeof(value_t),MPI_BYTE,run2fft->tasks[k],30,c.com,run2fft->reqs+k);
	}
    }

  if(c.is_fft)
    {
      /*
      //  Use num_fft_buffers buffers of the largest size.
      //  We start with 1, because fft2run->num_tasks can be 0.
      */
      size = 1;
      for(k=0; k<fft2run->num_tasks; k++)
	{
	  if(size < fft2run->sizes[k]) size = fft2run->sizes[k];
	}

      rbuffers = cart_alloc(value_t*,num_fft_buffers);
      for(k=0; k<num_fft_buffers; k++)
	{
	  rbuffers[k] = cart_alloc(value_t,size);
	}

      /*
      //  Allocate map if present.
      */
      if(c.map != NULL)
	{
	  c.map->lines = cart_alloc(ijk_t*,fft2run->num_tasks);
	  for(k=0; k<fft2run->num_tasks; k++)
	    {
	      c.map->lines[k] = cart_alloc(ijk_t,fft2run->sizes[k]);
	    }
	}

      fft2run->reqs = cart_alloc(MPI_Request,num_fft_buffers);
      for(k=l=kb=0; k<fft2run->num_tasks; k++)
	{
	  while(l<num_fft_buffers && kb<fft2run->num_tasks)
	    {
	      kb++;

	      MPI_Probe(MPI_ANY_SOURCE,30,c.com,&status);

	      for(idx=0; idx<fft2run->num_tasks; idx++)
		{
		  if(fft2run->tasks[idx] == status.MPI_SOURCE) break;
		}
	      cart_assert(idx < fft2run->num_tasks);

	      MPI_Irecv(rbuffers[l],fft2run->sizes[idx]*sizeof(value_t),MPI_BYTE,fft2run->tasks[idx],30,c.com,fft2run->reqs+l);
	      sender[l] = idx;
	      l++;
	    }

	  MPI_Waitany(l,fft2run->reqs,&lidx,MPI_STATUS_IGNORE);
	  idx = sender[lidx];
	  cart_assert(idx>=0 && idx<fft2run->num_tasks);

	  /*
	  //  Some data arrived, put it in and re-use rbuffers
	  */
#pragma omp parallel for default (none), private(i,j,ptr), shared(fft2run,idx,data,rbuffers,root_grid_fft_internal_data,lidx,root_grid_fft_internal_config)
	  for(i=0; i<fft2run->sizes[idx]; i++)
	    {
	      ptr = rbuffers[lidx] + i;

#ifdef __CHECK
	      ijk_t index = ptr->index;
	      for(j=0; j<nDim; j++)
		{
		  if((index%d.dims[j])<0 || (index%d.dims[j])>=d.bbox[2*j+1]-d.bbox[2*j+0])
		    {
		      cart_error("Received value is outside BBox, d=%d, val=%d, bbox=%d, from=%d",j,index%d.dims[j],d.bbox[2*j+1]-d.bbox[2*j+0],fft2run->tasks[idx]);
		    }
		  index = index/d.dims[j];
		}
#endif
	      data[ptr->index] = ptr->value;
	      if(c.map != NULL)
		{
		  c.map->lines[idx][i] = ptr->index;
		}
	    }

	  /*
	  // Free buffer #lidx
	  */
	  l--;
	  ptr = rbuffers[lidx];
	  rbuffers[lidx] = rbuffers[l];
	  rbuffers[l] = ptr;
	  sender[lidx] = sender[l];
	  fft2run->reqs[lidx] = fft2run->reqs[l];
	}

      cart_free(fft2run->reqs);

      for(k=0; k<root_grid_fft_tune->num_fft_buffers; k++)
	{
	  cart_free(rbuffers[k]);
	}
      cart_free(rbuffers);
    }

  if(c.is_run)
    {
      MPI_Waitall(run2fft->num_tasks,run2fft->reqs,MPI_STATUSES_IGNORE);
      cart_free(run2fft->reqs);

      for(k=0; k<run2fft->num_tasks; k++)
	{
	  cart_free(sbuffers[k]);
	}
      cart_free(sbuffers);
    }
}


void root_grid_fft_internal_recv_data_no_map(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells, int out_var, fft_t *data);
void root_grid_fft_internal_recv_data_mapped(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells, int out_var, fft_t *data);


void root_grid_fft_internal_recv_data(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells, int out_var, fft_t *data)
{
  /*
  //  We need to send the data back to run tasks. If there is no map set,
  //  we have just one choice - for each fft task to send the whole bbox
  //  to every run task the fft task is connected to, and let the run tasks
  //  to select the data they need. If there is a map, then we can send each 
  //  run task only the data it needs.
  */
  if(root_grid_fft_tune->use_recv_map)
    {
      root_grid_fft_internal_recv_data_mapped(run2fft,fft2run,num_level_cells,level_cells,out_var,data);
    }
  else
    {
      root_grid_fft_internal_recv_data_no_map(run2fft,fft2run,num_level_cells,level_cells,out_var,data);
    }
}


void root_grid_fft_internal_recv_data_no_map(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells, int out_var, fft_t *data)
{
  int num_slices = (root_grid_fft_tune->num_recv_slices > 0) ? root_grid_fft_tune->num_recv_slices : d.dims[nDim-1];
  int i, k, l, task_idx, page_idx, ijk[nDim];
  fft_t *buffer;
  int num_pages;
  size_t slice_size, page_size, index;
  MPI_Status status;

  slice_size = 1;
  for(k=0; k<nDim-1; k++) slice_size *= d.dims[k];

  num_pages = (d.dims[nDim-1]+num_slices-1)/num_slices;

  if(c.is_fft)
    {
      /*
      //  Send our whole data array to all run tasks. No need to 
      //  create any buffers. To avoid forcing a run task to 
      //  duplicate the whole data array, we send a few k-slices at 
      //  a time;
      */
      fft2run->reqs = cart_alloc(MPI_Request,fft2run->num_tasks*num_pages);

      for(l=0; l<num_pages; l++)
	{
	  if(l == num_pages-1)
	    {
	      page_size = slice_size*(d.dims[nDim-1]-num_slices*l);
	    }
	  else
	    {
	      page_size = slice_size*num_slices;
	    }
	  for(k=0; k<fft2run->num_tasks; k++)
	    {
	      MPI_Isend(data+l*slice_size*num_slices,page_size*sizeof(fft_t),MPI_BYTE,fft2run->tasks[k],40+l,c.com,fft2run->reqs+k+fft2run->num_tasks*l);
	    }
	}
    }

  if(c.is_run)
    {
      buffer = cart_alloc(fft_t,slice_size*num_slices);

      for(l=0; l<num_pages; l++)
	{
	  for(k=0; k<run2fft->num_tasks; k++)
	    {
	      MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,c.com,&status);
	  
	      for(task_idx=0; task_idx<run2fft->num_tasks; task_idx++)
		{
		  if(run2fft->tasks[task_idx] == status.MPI_SOURCE) break;
		}
	      cart_assert(task_idx < run2fft->num_tasks);

 	      page_idx = status.MPI_TAG - 40;
	      if(page_idx == num_pages-1)
		{
		  page_size = slice_size*(d.dims[nDim-1]-num_slices*page_idx);
		}
	      else
		{
		  page_size = slice_size*num_slices;
		}

	      MPI_Recv(buffer,page_size*sizeof(fft_t),MPI_BYTE,run2fft->tasks[task_idx],status.MPI_TAG,c.com,MPI_STATUSES_IGNORE);

	      /*
	      //  We need to scan our whole root grid segment to find a match
	      */
#pragma omp parallel for default(none), private(i,ijk,index), shared(num_level_cells,level_cells,buffer,out_var,run2fft,task_idx,page_idx,page_size,num_slices,cell_vars,root_grid_fft_internal_config)
	      for(i=0; i<num_level_cells; i++)
		{
		  if(c.cell == NULL)
		    {
		      root_grid_fft_get_cell_ijk(level_cells[i],ijk);
		      index = root_grid_fft_get_page_index(ijk,run2fft->tasks[task_idx],num_slices,page_idx);
		    }
		  else
		    {
		      index = root_grid_fft_get_page_index(c.cell[i].ijk,run2fft->tasks[task_idx],num_slices,page_idx);
		    }

		  if(index < page_size) // size_t is unsigned
		    {
		      cell_var(level_cells[i],out_var) = buffer[index];
		    }
		}
	    }
	}

      cart_free(buffer);
    }

  if(c.is_fft)
    {
      MPI_Waitall(fft2run->num_tasks*num_pages,fft2run->reqs,MPI_STATUSES_IGNORE);
      cart_free(fft2run->reqs);
    }
}


void root_grid_fft_internal_recv_data_mapped(cache_t *run2fft, cache_t *fft2run, int num_level_cells, const int *level_cells, int out_var, fft_t *data)
{
  int i, k, ijk[nDim];
  fft_t **sbuffers, **rbuffers;
  int *counts, size;
  int l, num_fft_buffers = MAX(1,root_grid_fft_tune->num_fft_buffers);

  /*
  //  This is just send_data in the reverse order. We just need to be
  //  careful to set recvs before sends.
  */
  if(c.is_run)
    {
      sbuffers = cart_alloc(fft_t*,run2fft->num_tasks);
      for(k=0; k<run2fft->num_tasks; k++)
	{
	  sbuffers[k] = cart_alloc(fft_t,run2fft->sizes[k]);
	}

      run2fft->reqs = cart_alloc(MPI_Request,run2fft->num_tasks);
      for(k=0; k<run2fft->num_tasks; k++)
	{
	  MPI_Irecv(sbuffers[k],run2fft->sizes[k]*sizeof(fft_t),MPI_BYTE,run2fft->tasks[k],50,c.com,run2fft->reqs+k);
	}
    }

  if(c.is_fft)
    {
      /*
      //  Use num_fft_buffers buffers of the largest size.
      //  We start with 1, because fft2run->num_tasks can be 0.
      */
      size = 1;
      for(k=0; k<fft2run->num_tasks; k++)
	{
	  if(size < fft2run->sizes[k]) size = fft2run->sizes[k];
	}

      rbuffers = cart_alloc(fft_t*,num_fft_buffers);
      for(k=0; k<num_fft_buffers; k++)
	{
	  rbuffers[k] = cart_alloc(fft_t,size);
	}

      fft2run->reqs = cart_alloc(MPI_Request,num_fft_buffers);
      for(k=0; k<fft2run->num_tasks; k+=num_fft_buffers)
	{
	  for(l=0; l<num_fft_buffers && k+l<fft2run->num_tasks; l++)
	    {
	      /*
	      //  Pack the data
	      */
#pragma omp parallel for default (none), private(i), shared(fft2run,k,l,data,rbuffers,root_grid_fft_internal_data,root_grid_fft_internal_config)
	      for(i=0; i<fft2run->sizes[k+l]; i++)
		{
		  rbuffers[l][i] = data[c.map->lines[k+l][i]];
		}

	      MPI_Isend(rbuffers[l],fft2run->sizes[k+l]*sizeof(fft_t),MPI_BYTE,fft2run->tasks[k+l],50,c.com,fft2run->reqs+l);
	    }

	  MPI_Waitall(l,fft2run->reqs,MPI_STATUS_IGNORE);
	}
      cart_free(fft2run->reqs);

      for(k=0; k<num_fft_buffers; k++)
	{
	  cart_free(rbuffers[k]);
	}
      cart_free(rbuffers);
    }

  if(c.is_run)
    {
      MPI_Waitall(run2fft->num_tasks,run2fft->reqs,MPI_STATUSES_IGNORE);
      cart_free(run2fft->reqs);

      counts = cart_alloc(int,run2fft->num_tasks);
      memset(counts,0,run2fft->num_tasks*sizeof(int));

      /*
      //  Serial loop. It is not easy to make it OpenMP parallel.
      */
      for(i=0; i<num_level_cells; i++)
	{
	  if(c.cell == NULL)
	    {
	      root_grid_fft_get_cell_ijk(level_cells[i],ijk);
	      k = root_grid_fft_get_task_idx(ijk,run2fft);
	    }
	  else
	    {
	      k = root_grid_fft_get_task_idx(c.cell[i].ijk,run2fft);
	    }
	  cart_assert(k>=0 && k<run2fft->num_tasks);

	  cell_var(level_cells[i],out_var) = sbuffers[k][counts[k]];
	  counts[k]++;
	}

      for(k=0; k<run2fft->num_tasks; k++)
	{
	  cart_assert(counts[k] == run2fft->sizes[k]);
	}

      cart_free(counts);
      for(k=0; k<run2fft->num_tasks; k++)
	{
	  cart_free(sbuffers[k]);
	}
      cart_free(sbuffers);
    }
}


/*
//  Domain mapping for fft tasks. In principle, it can be arbitrary, 
//  if bboxes from all fft tasks are stored at each run task.
//  Current implementation assumes slab decomposition in Z
//  direction as used by the fft3 library.
*/
void root_grid_fft_get_cell_ijk(int cell, int ijk[nDim])
{
  sfc_coords( root_cell_sfc_index(cell), ijk );
}


int root_grid_fft_get_task(int ijk[nDim])
{
  return ijk[nDim-1]/d.dims[nDim-1];
}


int root_grid_fft_get_task_idx(int ijk[nDim], const cache_t *cache)
{
  int i, task_id = root_grid_fft_get_task(ijk);

  for(i=0; i<cache->num_tasks; i++)
    {
      if(cache->tasks[i] == task_id) break;
    }
  cart_assert(i < cache->num_tasks);

  return i;
}


size_t root_grid_fft_get_page_index(int ijk[nDim], int task, int num_slices, int page)
{
  return ijk[0]
#if (nDim > 1)
    +d.dims[0]*(ijk[1]
#if (nDim > 2)
		+(size_t)d.dims[1]*(ijk[2]-task*d.dims[2]-page*num_slices)
#else
		-task*d.dims[1]-page*num_slices
#endif
		)
#else
    - task*d.dims[0] - page*num_slices;
#endif
    ;
}


void root_grid_fft_print_map(const char *title, cache_t *cache)
{
  int i;
  char str[999], str1[24];

  strncpy(str,title,990);

  for(i=0; i<cache->num_tasks; i++)
    {
      if(strlen(str) > 999-24) break;
      snprintf(str1,23," %d[%d]",cache->tasks[i],cache->sizes[i]);
      strcat(str,str1);
    }
  cart_debug(str);
}

