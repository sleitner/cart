#include "config.h"
#if defined(RADIATIVE_TRANSFER) && defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)

#include <math.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "control_parameter.h"
#include "iterators.h"
#include "parallel.h"
#include "root_grid_fft.h"
#include "rt_global.h"
#include "rt_otvet.h"
#include "rt_transfer.h"
#include "timing.h"
#include "tree.h"

#include "../tools/fft/fft3.h"


/*
//  These can be changed to avoid excessive cache missing
*/
#define NUM_LEAD      (num_grid+2)
#define NUM_TRAIL     num_grid

#define GREEN_SIZE    ((size_t)num_grid*NUM_TRAIL*NUM_LEAD/2)

fft_t *rtGreenET[] = { 0, 0, 0, 0, 0, 0 };

int rtNumOtvetETVars = 1 + rt_num_et_vars;
int rtOtvetETVars[] = { RT_VAR_OT_FIELD,
                        rt_et_offset + 0,
                        rt_et_offset + 1,
                        rt_et_offset + 2,
                        rt_et_offset + 3,
                        rt_et_offset + 4,
                        rt_et_offset + 5 };

int rtOtvetOTBox[rt_num_freqs];


/*                          6  7  8  9 10 11 12 13 14 15 16 17  */
const int StencilDir1[] = { 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 2, 3 };
const int StencilDir2[] = { 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5 };

double rtStencilDist2[rtStencilSize];
double rtStencilDelPos[rtStencilSize][nDim];
double rtStencilTensor[rtStencilSize][nDim*(nDim+1)/2];


int rtOtvetMaxNumIter = 10;


void rtOtvetEddingtonTensor(int level, int num_level_cells, int *level_cells);
void rtOtvetTreeEmulatorEddingtonTensor(int level, int num_level_cells, int *level_cells);
void rtOtvetSingleSourceEddingtonTensor(int level);
void rtOtvetTopLevelFFTWorker(const root_grid_fft_t *config, int id, fft_t *fft_source, fft_t *fft_output, int flags);

void root_grid_fft_get_cell_ijk(int cell, int ijk[nDim]);


void rtConfigInitTransferOtvet()
{
  control_parameter_add2(control_parameter_int,&rtOtvetMaxNumIter,"rt:otvet:num-iterations","rtOtvetMaxNumIter","the number of iterations in the OTVET solver.");
}


void rtConfigVerifyTransferOtvet()
{
  VERIFY(rt:otvet:num-iterations, (rtOtvetMaxNumIter > 0) );
}


void rtInitRunTransferOtvet()
{
  int i, j, k, l;
  double r2;

  /*
  //  Fill in Stencil positions
  */
  for(l=0; l<rtStencilSize; l++)
    {
      for(j=0; j<nDim; j++) rtStencilDelPos[l][j] = 0.0;
    }
     
  for(l=0; l<num_neighbors; l++) rtStencilDelPos[l][l/2] = 2*(l%2) - 1;
  for(l=0; l<rtStencilSize-num_neighbors; l++)
    { 
      rtStencilDelPos[num_neighbors+l][StencilDir1[l]/2] = 2*(StencilDir1[l]%2) - 1;
      rtStencilDelPos[num_neighbors+l][StencilDir2[l]/2] = 2*(StencilDir2[l]%2) - 1;
    }

  for(l=0; l<rtStencilSize; l++)
    {
      r2 = 0.0;
      for(j=0; j<nDim; j++) r2 += rtStencilDelPos[l][j]*rtStencilDelPos[l][j];
      rtStencilDist2[l] = r2;

      for(k=j=0; j<nDim; j++)
        {
          for(i=0; i<=j; i++)
            {
	      rtStencilTensor[l][k++] = rtStencilDelPos[l][i]*rtStencilDelPos[l][j]/r2;
            }
        }
    }
}


void rtInitStepTransferOtvet(struct rtGlobalValue *maxAC)
{
  const float tauMin = 1.0e-2;
  int freq;

  for(freq=0; freq<rt_num_freqs; freq++) rtOtvetOTBox[freq] = (maxAC[freq].Value*num_grid < tauMin);
}


void rtAfterAssignDensityTransferOtvet(int level, int num_level_cells, int *level_cells)
{
  rtOtvetEddingtonTensor(level,num_level_cells,level_cells);
}


/*
// This routine returns 18 neighbors as if the whole mesh was uniform.
// In particular, one neighbor of higher level will appear more than
// once. The first 6 neighbors are exactly as returned by cell_all_neighbors.
// The other 12 neighbors are packed as follows:
//  6:  --0
//  7:  +-0
//  8:  -+0
//  9:  ++0
// 10:  -0-
// 11:  +0-
// 12:  -0+
// 13:  +0+
// 14:  0--
// 15:  0+-
// 16:  0-+
// 17:  0++
*/
void rtGetStencil(int level, int cell, int nb[])
{
  int j, levNb[num_neighbors];
  
  /*
  //  First level neighbors
  */
  cell_all_neighbors(cell,nb);
  for(j=0; j<num_neighbors; j++) levNb[j] = cell_level(nb[j]);
  
  /*
  //  Second level neighbors
  */
  for(j=0; j<rtStencilSize-num_neighbors; j++)
    {
      if(levNb[StencilDir1[j]] == level)
	{
	  /*
	  // Our neighbor in the first direction is on the same level,
	  // it is sufficient to take its neighbor in the second direction.
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[StencilDir1[j]],StencilDir2[j]);
	}
      else if(levNb[StencilDir2[j]] == level)
	{
	  /*
	  // Our neighbor in the second direction is on the same level,
	  // it is sufficient to take its neighbor in the first direction.
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[StencilDir2[j]],StencilDir1[j]);
	}
      else
	{
	  /*
	  // Both our neighbors are on a higher level. In that case the corner cell
	  // cannot be our immediate neighbor, it must be a common neighbor of
	  // two higher level cells.
	  */
	  nb[num_neighbors+j] = cell_neighbor(nb[StencilDir1[j]],StencilDir2[j]);
#ifdef DEBUG
	  cart_assert(nb[num_neighbors+j] == cell_neighbor(nb[StencilDir2[j]],StencilDir1[j]));
#endif      
	}
    }


#ifdef DEBUG
  double p0[nDim], p[nDim];
  int k;
 
  cell_center_position(cell,p0);

  for(j=0; j<rtStencilSize-num_neighbors; j++)
    {
      for(k=0; k<nDim; k++)
	{
	  p[k] = p0[k] + cell_size[level]*rtStencilDelPos[num_neighbors+j][k];
	  if(p[k] < 0.0) p[k] += num_grid;
	  if(p[k] > num_grid) p[k] -= num_grid;
	}

      cart_assert(cell_contains_position(nb[num_neighbors+j],p));
    }
#endif

}


void rtOtvetEddingtonTensor(int level, int num_level_cells, int *level_cells)
{
#ifdef RT_SINGLE_SOURCE

  rtOtvetSingleSourceEddingtonTensor(level);

#else

  rtOtvetTreeEmulatorEddingtonTensor(level,num_level_cells,level_cells);

#endif
}


void rtOtvetTreeEmulatorEddingtonTensor(int level, int num_level_cells, int *level_cells)
{
  const int NumSmooth = 2;
  const double S1 = 1.0;
  const double S2 = (S1-1)/nDim;

  int i, j, k, l, cell, parent;
  int nb3[nDim], nb18[rtStencilSize];
  int num_parent_cells, *parent_cells;
  float norm, h2, eps2, *tmp;
  float ot, et[rt_num_et_vars], sor;
  double r2, q;

  if(level == min_level)
    {
      root_grid_fft_exec(RT_VAR_SOURCE,rtNumOtvetETVars-1,rtOtvetETVars+1,rtOtvetTopLevelFFTWorker);

      start_time(WORK_TIMER);

      /*
      //  Normalize
      */
#pragma omp parallel for default(none), private(i,j,cell), shared(num_level_cells,level_cells,cell_vars)
      for(i=0; i<num_level_cells; i++)
	{
	  cell = level_cells[i];

	  cell_var(cell,RT_VAR_OT_FIELD) = cell_var(cell,rt_et_offset+0)
#if (nDim > 1)
	    + cell_var(cell,rt_et_offset+2)
#if (nDim > 2)
	    + cell_var(cell,rt_et_offset+5)
#endif /* nDim > 2 */
#endif /* nDim > 1 */
	    ;

	  if(cell_var(cell,RT_VAR_OT_FIELD) > 0.0)
	    {
	      for(j=0; j<rt_num_et_vars; j++)
		{
		  cell_var(cell,rt_et_offset+j) /= cell_var(cell,RT_VAR_OT_FIELD);
		}
	    }
	  else
	    {
	      cell_var(cell,RT_VAR_OT_FIELD) = 0.0;
	      cell_var(cell,rt_et_offset+0) = 1.0/nDim;
#if (nDim > 1)
	      cell_var(cell,rt_et_offset+1) = 0.0;
	      cell_var(cell,rt_et_offset+2) = 1.0/nDim;
#if (nDim > 2)
	      cell_var(cell,rt_et_offset+3) = 0.0;
	      cell_var(cell,rt_et_offset+4) = 0.0;
	      cell_var(cell,rt_et_offset+5) = 1.0/nDim;
#endif /* nDim > 2 */
#endif /* nDim > 1 */
	    }
	}

      end_time(WORK_TIMER);

    }
  else
    {

      start_time(WORK_TIMER);

      /*
      // We start with interpolating from parents
      */
      select_level(level-1,CELL_TYPE_LOCAL | CELL_TYPE_REFINED,&num_parent_cells,&parent_cells);

      h2 = cell_size[level]*cell_size[level];
      norm = cell_volume[level]/(4*M_PI*h2);

      if(level == 1)
	{
	  eps2 = 0.1;
	}
      else
	{
	  eps2 = 0.05;
	}

#pragma omp parallel for default(none), private(i,j,k,parent,cell,ot,et,nb3,nb18,sor,l,r2,q), shared(level,num_parent_cells,parent_cells,cell_vars,cell_child_oct,norm,h2,rtStencilDist2,rtStencilTensor,eps2)
      for(i=0; i<num_parent_cells; i++)
	{
	  parent = parent_cells[i];

	  /*
	  //  Loop over all children
	  */
	  for(j=0; j<num_children; j++)
	    {
	      /*
	      //  Interpolate from parents and turn ET into OT radiation pressure tensor
	      */
	      cell_interpolation_neighbors(parent,j,nb3);

	      /*
	      //  NG: this is the best interpolation, I did check two other forms
	      */
	      ot = cell_interpolate_with_neighbors(parent,RT_VAR_OT_FIELD,nb3);
	      for(k=0; k<rt_num_et_vars; k++)
		{
		  et[k] = ot*cell_interpolate_with_neighbors(parent,rt_et_offset+k,nb3);
		}

	      cell = cell_child(parent,j);
	      
	      /*
	      //  Add local contributions from 18 neighbors
	      */
	      rtGetStencil(level,cell,nb18);
	      for(l=0; l<rtStencilSize; l++) if((sor = cell_rt_source(nb18[l])) > 0.0)
		{
		  r2 = rtStencilDist2[l];

		  sor *= norm/(r2+eps2);

		  ot += sor;

		  et[0] += sor*(S1*rtStencilTensor[l][0]-S2);
#if (nDim > 1)
		  et[1] += sor*rtStencilTensor[l][1];
		  et[2] += sor*(S1*rtStencilTensor[l][2]-S2);
#if (nDim > 2)
		  et[3] += sor*rtStencilTensor[l][3];
		  et[4] += sor*rtStencilTensor[l][4];
		  et[5] += sor*(S1*rtStencilTensor[l][5]-S2);
#endif /* nDim > 2 */
#endif /* nDim > 1 */
		}

	      cell_var(cell,RT_VAR_OT_FIELD) = ot;
	      if(ot > 0.0)
		{
		  q = et[0]
#if (nDim > 1)
		    + et[2]
#if (nDim > 2)
		    + et[5]
#endif /* nDim > 2 */
#endif /* nDim > 1 */
		    ;
		  for(k=0; k<rt_num_et_vars; k++)
		    {
		      cell_var(cell,rt_et_offset+k) = et[k]/q;
		    }
		}
	      else
		{
		  cell_var(cell,rt_et_offset+0) = 1.0/nDim;
#if (nDim > 1)
		  cell_var(cell,rt_et_offset+1) = 0.0;
		  cell_var(cell,rt_et_offset+2) = 1.0/nDim;
#if (nDim > 2)
		  cell_var(cell,rt_et_offset+3) = 0.0;
		  cell_var(cell,rt_et_offset+4) = 0.0;
		  cell_var(cell,rt_et_offset+5) = 1.0/nDim;
#endif /* nDim > 2 */
#endif /* nDim > 1 */
		}
	    }
	}
      
      cart_free(parent_cells);

      end_time(WORK_TIMER);

    }
 
  start_time(RT_TREE_EMULATOR_UPDATE_TIMER);
  update_buffer_level(level,rtOtvetETVars,rtNumOtvetETVars);
  end_time(RT_TREE_EMULATOR_UPDATE_TIMER);

  /*
  // Smooth a few times
  */
  start_time(WORK_TIMER);
  tmp = cart_alloc(float, num_level_cells*rt_num_et_vars);
  end_time(WORK_TIMER);

  for(l=0; l<NumSmooth; l++)
    {

      start_time(WORK_TIMER);

#pragma omp parallel for default(none), private(i,j,k,cell,et,nb18), shared(level,num_level_cells,level_cells,cell_vars,tmp)
      for(i=0; i<num_level_cells; i++)
	{
	  cell = level_cells[i];

	  rtGetStencil(level,cell,nb18);
	  //cell_all_neighbors(cell,nb18);

	  for(k=0; k<rt_num_et_vars; k++) et[k] = 2*cell_var(cell,rt_et_offset+k);

	  for(j=0; j<rtStencilSize; j++)
	    //for(j=0; j<num_neighbors; j++)
	    {
	      for(k=0; k<rt_num_et_vars; k++)
		{
		  et[k] += cell_var(nb18[j],rt_et_offset+k);
		}
	    }

	  for(k=0; k<rt_num_et_vars; k++) tmp[i*rt_num_et_vars+k] = et[k]/(2+rtStencilSize);
	}

#pragma omp parallel for default(none), private(i,k,cell), shared(level,num_level_cells,level_cells,cell_vars,tmp)
      for(i=0; i<num_level_cells; i++)
	{
	  cell = level_cells[i];
	  
	  for(k=0; k<rt_num_et_vars; k++)
	    {
	      cell_var(cell,rt_et_offset+k) = tmp[i*rt_num_et_vars+k];
	    }
	}

      end_time(WORK_TIMER);

      start_time(RT_TREE_EMULATOR_UPDATE_TIMER);
      update_buffer_level(level,rtOtvetETVars+1,rtNumOtvetETVars-1);
      end_time(RT_TREE_EMULATOR_UPDATE_TIMER);
    }

  start_time(WORK_TIMER);
  cart_free(tmp);
  end_time(WORK_TIMER);

}


void rtOtvetTopLevelFFTWorker(const root_grid_fft_t *config, int id, fft_t *fft_source, fft_t *fft_output, int flags)
{
  static int never_run = 1;
  int i, j, k, m;
  size_t size, offset;
  double dx, dy, dz, r[6], r2, q;
  double norm;
  fft_t *gf;

  cart_assert(config->bbox[0]==0 && config->bbox[1]==num_grid);

  if(never_run)
    {
      never_run = 0;

      size = (size_t)config->dims[0]*config->dims[1]*config->dims[2];

      for(m=0; m<6; m++)
	{
	  rtGreenET[m] = cart_alloc(fft_t,size);
	}

      norm = 1.0/(4*M_PI*num_grid*num_grid*num_grid);

#pragma omp parallel for default(none), private(i,j,k,m,offset,dx,dy,dz,r,r2,q), shared(rtGreenET,norm,config)
      for(k=config->bbox[4]; k<config->bbox[5]; k++)
	{
	  dz = k;
	  if(dz > num_grid/2) dz -= num_grid;
	  for(j=config->bbox[2]; j<config->bbox[3]; j++)
	    {
	      dy = j;
	      if(dy > num_grid/2) dy -= num_grid;

	      offset = config->dims[0]*(j-config->bbox[2]+(size_t)config->dims[1]*(k-config->bbox[4]));

	      for(i=config->bbox[0]; i<config->bbox[1]; i++)
		{
		  dx = i;
		  if(dx > num_grid/2) dx -= num_grid;

		  if(k==0 && j==0 && i==0)
		    {
		      r[0] = r[2] = r[5] = 1.0/nDim;
		      r[1] = r[3] = r[4] = 0.0;
		    }
		  else
		    {
		      r[0] = dx*dx;
		      r[1] = dx*dy;
		      r[2] = dy*dy;
		      r[3] = dx*dz;
		      r[4] = dy*dz;
		      r[5] = dz*dz;
		    }

		  r2 = r[0]
#if (nDim > 1)
		    + r[2]
#if (nDim > 2)
		    + r[5]
#endif /* nDim > 2 */
#endif /* nDim > 1 */
		    ;

		  /*
		  //  Correction to the Green function
		  */
		  if(r2 < 8.0)
		    {
		      if(r2 > 1.0)
			{
			  q = 0.7 + 0.3*(r2/8.0);
			}
		      else q = 1.25;
		    }
		  else
		    {
		      q = 1.0;
		    }
		  
		  for(m=0; m<6; m++)
		    {
		      rtGreenET[m][i+offset] = norm*r[m]/(r2*r2)*q;
		    }
		}
	    }
	}

      for(m=0; m<6; m++)
	{
	  fft3_x2k(rtGreenET[m],flags);
	}

    }

  gf = rtGreenET[id];

#pragma omp parallel for default(none), private(i,j,k,offset), shared(fft_source,fft_output,gf,config,flags)
  for(k=0; k<config->dims[2]; k++)
    {
      for(j=0; j<config->dims[1]; j++)
	{
	  offset = fft3_jk_index(j,k,NULL,flags);
	  if(offset == (size_t)-1) continue;
	  
	  offset *= config->dims[0];

	  for(i=0; i<=num_grid/2; i++)
	    {
	      fft_output[2*i+0+offset] = 
		fft_source[2*i+0+offset]*gf[2*i+0+offset] - 
		fft_source[2*i+1+offset]*gf[2*i+1+offset];
	      fft_output[2*i+1+offset] = 
		fft_source[2*i+1+offset]*gf[2*i+0+offset] +
		fft_source[2*i+0+offset]*gf[2*i+1+offset];
	    }
	}
    }
}


#ifdef RT_SINGLE_SOURCE

void rtOtvetSingleSourceEddingtonTensor(int level)
{
  int i, j, l, index, cell;
  int num_level_cells, *level_cells;
  double eps1, eps2, dr2, pos[nDim];

  start_time(WORK_TIMER);

  eps1 = 0.01*cell_size[rtSingleSourceLevel]*cell_size[rtSingleSourceLevel];
  eps2 = 9*cell_size[rtSingleSourceLevel]*cell_size[rtSingleSourceLevel];

  select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);

#pragma omp parallel for default(none), private(index,cell,pos,dr2,i,j,l), shared(level,num_level_cells,level_cells,rtSingleSourceValue,rtSingleSourcePos,eps1,eps2,cell_vars)
  for(index=0; index<num_level_cells; index++)
    {
      cell = level_cells[index];
      cell_center_position(cell,pos);

      dr2 = eps1;
      for(i=0; i<nDim; i++)
	{
	  pos[i] -= rtSingleSourcePos[i];
	  if(pos[i] >  0.5*num_grid) pos[i] -= num_grid;
	  if(pos[i] < -0.5*num_grid) pos[i] += num_grid;
	  dr2 += pos[i]*pos[i];
	}

      cell_var(cell,RT_VAR_OT_FIELD) = rtSingleSourceValue/(4*M_PI*dr2);
      
      dr2 += nDim*eps2;
      for(l=j=0; j<nDim; j++)
	{
	  for(i=0; i<j; i++)
	    {
	      cell_var(cell,rt_et_offset+l++) = pos[i]*pos[j]/dr2;
	    }
	  cell_var(cell,rt_et_offset+l++) = (eps2+pos[j]*pos[j])/dr2;
	}
    }

  cart_free(level_cells);

  end_time(WORK_TIMER);

  start_time(RT_SINGLE_SOURCE_UPDATE_TIMER);
  update_buffer_level(level,rtOtvetETVars,rtNumOtvetETVars);
  end_time(RT_SINGLE_SOURCE_UPDATE_TIMER);
}

#endif /* RT_SINGLE_SOURCE */

#endif /* RADIATIVE_TRANSFER && RT_TRANSFER && (RT_TRANSFER_METHOD == RT_METHOD_OTVET) */
