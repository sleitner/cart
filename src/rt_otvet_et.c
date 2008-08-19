#include "defs.h"
#ifdef RADIATIVE_TRANSFER

#include <math.h>

#include "auxiliary.h"
#include "cell_buffer.h"
#include "top_level_fft.h"
#include "tree.h"

#include "rt_tree.h"
#include "rt_utilities.h"


#if defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)


#include "rt_transfer.h"


/*
//  These can be changed to avoid excessive cache missing
*/
#define NUM_LEAD      (num_grid+2)
#define NUM_TRAIL     num_grid

#define GREEN_SIZE    ((long)num_grid*(long)NUM_TRAIL*(long)NUM_LEAD/2)

fftw_complex *rtGreenET[] = { 0, 0, 0, 0, 0, 0 };

extern const int NumVars;
extern const int Vars[];


void rtOtvetTopLevelEddingtonTensor(int id, fftw_complex *fft_source, fftw_complex *fft_output);
void rtOtvetTreeEmulatorEddingtonTensor(int level);
void rtOtvetSingleSourceEddingtonTensor(int level, float srcVal, double *srcPos, double srcEps);

#ifdef RT_VAR_SOURCE
float rtSource(int ipart);
#endif


void rtOtvetEddingtonTensor(int level)
{
#ifdef RT_SINGLE_SOURCE

  rtOtvetSingleSourceEddingtonTensor(level,rtSingleSourceVal,rtSingleSourcePos,rtSingleSourceEps);

#else

  rtOtvetTreeEmulatorEddingtonTensor(level);

#endif

  update_buffer_level(level,Vars,NumVars);
  rtTransferUpdateFields(level,level,1,RT_VAR_OT_FIELD,rt_glob_Avg+1);
}


void rtOtvetTreeEmulatorEddingtonTensor(int level)
{
  const double S1 = 1.0;
  const double S2 = (S1-1)/nDim;

  int i, j, k, l, cell, parent;
  int nb3[nDim], nb6[num_neighbors], nb18[rtuStencilSize];
  int num_level_cells, *level_cells;
  float norm, h2, eps2;
  float ot, et[rt_num_et_vars], sor;
  double r2, q;

  if(level == min_level)
    {
      top_level_fft(RT_VAR_SOURCE,NumVars-1,Vars+1,rtOtvetTopLevelEddingtonTensor);

      /*
      //  Normalize
      */
      select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);

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
    }
  else
    {
      /*
      // We start with interpolating from parents
      */
      select_level(level-1,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);

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

#pragma omp parallel for default(none), private(i,j,k,parent,cell,ot,et,nb3,nb18,sor,l,r2,q), shared(level,num_level_cells,level_cells,cell_vars,cell_child_oct,norm,h2,rtuStencilDist2,rtuStencilTensor,eps2)
      for(i=0; i<num_level_cells; i++) if(!cell_is_leaf(level_cells[i]))
	{
	  parent = level_cells[i];

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
	      rtuGetStencil(level,cell,nb18);
	      for(l=0; l<-rtuStencilSize; l++) if((sor = cell_rt_source(nb18[l])) > 0.0)
		{
		  r2 = rtuStencilDist2[l];

		  sor *= norm/(r2+eps2);

		  ot += sor;

		  et[0] += sor*(S1*rtuStencilTensor[l][0]-S2);
#if (nDim > 1)
		  et[1] += sor*rtuStencilTensor[l][1];
		  et[2] += sor*(S1*rtuStencilTensor[l][2]-S2);
#if (nDim > 2)
		  et[3] += sor*rtuStencilTensor[l][3];
		  et[4] += sor*rtuStencilTensor[l][4];
		  et[5] += sor*(S1*rtuStencilTensor[l][5]-S2);
#endif /* nDim > 2 */
#endif /* nDim > 1 */
		}

	      cell_var(cell,RT_VAR_OT_FIELD) = ot;
	      q = et[0] + et[2] + et[5];
	      for(k=0; k<rt_num_et_vars; k++)
		{
		  //cell_var(cell,rt_et_offset+k) = min(1.0,et[k]/q);
		  cell_var(cell,rt_et_offset+k) = et[k]/q;
		}
	    }
	}
    }

  cart_free(level_cells);
}


void rtOtvetComputeGreenFunctions()
{
  int i, j, k, m;
  long offset;
  float dx, dy, dz, r[6], r2, q;
  float norm;
  fftwnd_plan forward, backward;
  fftw_real *gf[6];

  for(m=0; m<6; m++)
    {
      rtGreenET[m] = (fftw_complex *)cart_alloc(GREEN_SIZE*sizeof(fftw_complex));
      gf[m] = (fftw_real *)rtGreenET[m];
    }

  norm = 1.0/(4*M_PI*num_grid*num_grid*num_grid);

#pragma omp parallel for default(none), private(i,j,k,m,offset,dx,dy,dz,r,r2,q), shared(gf,norm)
  for(i=0; i<num_grid; i++)
    {
      dx = i;
      if(dx > num_grid/2) dx -= num_grid;
      for(j=0; j<num_grid; j++)
	{
	  dy = j;
	  if(dy > num_grid/2) dy -= num_grid;

	  offset = NUM_LEAD*(j+NUM_TRAIL*i);
	  for(k=0; k<num_grid; k++)
	    {
	      dz = k;
	      if(dz > num_grid/2) dz -= num_grid;

	      if(offset==0 && k==0)
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

	      for(m=0; m<6; m++) gf[m][k+offset] = norm*r[m]/(r2*r2)*q;
	    }
	}
    }

  forward = rfftw3d_create_plan(num_grid,num_grid,num_grid,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE | FFTW_IN_PLACE);
  cart_assert(forward != 0);

  for(m=0; m<6; m++)
    {
      rfftwnd_one_real_to_complex(forward,gf[m],rtGreenET[m]);
    }
  
  rfftwnd_destroy_plan(forward);
}


void rtOtvetTopLevelEddingtonTensor(int id, fftw_complex *fft_source, fftw_complex *fft_output)
{
  int i, j, k;
  long offset;
  fftw_complex *gf = rtGreenET[id];

#pragma omp parallel for default(none), private(i,j,k,offset), shared(fft_source,fft_output,gf)
  for(i=0; i<num_grid; i++)
    {
      for(j=0; j<num_grid; j++)
	{
	  offset = NUM_LEAD*(j+NUM_TRAIL*i)/2;
	      
	  for(k=0; k<num_grid/2+1; k++)
	    {
	      fft_output[k+offset].re = 
		fft_source[k+offset].re*gf[k+offset].re - 
		fft_source[k+offset].im*gf[k+offset].im;
	      fft_output[k+offset].im = 
		fft_source[k+offset].im*gf[k+offset].re +
		fft_source[k+offset].re*gf[k+offset].im;
	    }
	}
    }
}


void rtOtvetSingleSourceEddingtonTensor(int level, float srcVal, double *srcPos, double srcEps)
{
  int i, j, l;
  int num_level_cells, *level_cells;
  int index, cell;
  double eps2, dr2, pos[nDim];

  eps2 = srcEps*srcEps;

  select_level(level,CELL_TYPE_LOCAL,&num_level_cells,&level_cells);

#pragma omp parallel for default(none), private(index,cell,pos,dr2,i,j,l), shared(level,num_level_cells,level_cells,srcVal,srcPos,eps2,cell_vars)
  for(index=0; index<num_level_cells; index++)
    {
      cell = level_cells[index];
      cell_position_double(cell,pos);

      dr2 = eps2;
      for(i=0; i<nDim; i++)
	{
	  pos[i] -= srcPos[i];
	  if(pos[i] >  0.5*num_grid) pos[i] -= num_grid;
	  if(pos[i] < -0.5*num_grid) pos[i] += num_grid;
	  dr2 += pos[i]*pos[i];
	}

      cell_var(cell,RT_VAR_OT_FIELD) = srcVal/(4*M_PI*dr2);

      for(l=j=0; j<nDim; j++)
	{
	  for(i=0; i<j; i++)
	    {
	      cell_var(cell,rt_et_offset+l++) = pos[i]*pos[j]/dr2;
	    }
	  cell_var(cell,rt_et_offset+l++) = (eps2/nDim+pos[j]*pos[j])/dr2;
	}
    }

  cart_free(level_cells);
}

#endif  // defined(RT_TRANSFER) && (RT_TRANSFER_METHOD == RT_METHOD_OTVET)
#endif  // RADIATIVE_TRANSFER
