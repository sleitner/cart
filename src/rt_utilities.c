#include "defs.h"      


#include <mpi.h>
#include <stdio.h>
#include <string.h>

#include "auxiliary.h"
#include "parallel.h"
#include "tree.h"

#include "rt_utilities.h"
#ifdef RADIATIVE_TRANSFER
#include "rt_solver.h"
#endif

/*
// **************************************************************
//
// General routines - used in RT block, but not requiring RT 
// include files (can be taken outside ifdef RADIATIVE_TRANSFER).
//
// **************************************************************
*/


/*                          6  7  8  9 10 11 12 13 14 15 16 17  */
const int StencilDir1[] = { 0, 1, 0, 1, 0, 1, 2, 3, 0, 1, 2, 3 };
const int StencilDir2[] = { 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5 };

double rtuStencilDist2[rtuStencilSize];
double rtuStencilDelPos[rtuStencilSize][nDim];
double rtuStencilTensor[rtuStencilSize][nDim*(nDim+1)/2];


void rtuInitRun()
{
  int i, j, k, l;
  double r2;

  /*
  //  Fill in Stencil positions
  */
  for(l=0; l<rtuStencilSize; l++)
    {
      for(j=0; j<nDim; j++) rtuStencilDelPos[l][j] = 0.0;
    }
     
  for(l=0; l<num_neighbors; l++) rtuStencilDelPos[l][l/2] = 2*(l%2) - 1;
  for(l=0; l<rtuStencilSize-num_neighbors; l++)
    { 
      rtuStencilDelPos[num_neighbors+l][StencilDir1[l]/2] = 2*(StencilDir1[l]%2) - 1;
      rtuStencilDelPos[num_neighbors+l][StencilDir2[l]/2] = 2*(StencilDir2[l]%2) - 1;
    }

  for(l=0; l<rtuStencilSize; l++)
    {
      r2 = 0.0;
      for(j=0; j<nDim; j++) r2 += rtuStencilDelPos[l][j]*rtuStencilDelPos[l][j];
      rtuStencilDist2[l] = r2;

      for(k=j=0; j<nDim; j++)
	{
	  for(i=0; i<=j; i++)
	    {
	      rtuStencilTensor[l][k++] = rtuStencilDelPos[l][i]*rtuStencilDelPos[l][j]/r2;
	    }
	}
    }
}


/*
//  Compute a global average of a buffer
*/
void rtuGlobalAverage(int n, double *lBuffer)
{
  double gBuffer[n];

  /*
    MPI_Reduce(lBuffer,gBuffer,n,MPI_DOUBLE,MPI_SUM,MASTER_NODE,MPI_COMM_WORLD);
    MPI_Bcast(gBuffer,n,MPI_DOUBLE,MASTER_NODE,MPI_COMM_WORLD);
  */
  MPI_Allreduce(lBuffer,gBuffer,n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  memcpy(lBuffer,gBuffer,n*sizeof(double));
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
void rtuGetStencil(int level, int cell, int nb[])
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
  for(j=0; j<rtuStencilSize-num_neighbors; j++)
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
 
  cell_position_double(cell,p0);

  for(j=0; j<rtuStencilSize-num_neighbors; j++)
    {
      for(k=0; k<nDim; k++)
	{
	  p[k] = p0[k] + cell_size[level]*rtuStencilDelPos[num_neighbors+j][k];
	  if(p[k] < 0.0) p[k] += num_grid;
	  if(p[k] > num_grid) p[k] -= num_grid;
	}

      cart_assert(cell_contains_position(nb[num_neighbors+j],p));
    }
#endif

}


/*
//  Compute the maxium and minimum of a (cached) 1-component array.
*/
void rtuGetLinearArrayMaxMin(int n, float *arr, float *max, float *min)
{
  int j, i, ibeg, iend;
  float *vmax, *vmin;
#ifdef _OPENMP
  int num_pieces = omp_get_num_threads();
#else
  int num_pieces = 1;
#endif 
  int len_piece = (n+num_pieces-1)/num_pieces;
  
  vmax = cart_alloc(num_pieces*sizeof(float));
  vmin = cart_alloc(num_pieces*sizeof(float));

#pragma omp parallel for default(none), private(j,i,ibeg,iend), shared(arr,vmin,vmax,n,len_piece,num_pieces)
  for(j=0; j<num_pieces; j++)
    {
      ibeg = j*len_piece;
      iend = ibeg + len_piece;
      if(iend > n) iend = n;
  
      vmin[j] = vmax[j] = arr[ibeg];

      for(i=ibeg+1; i<iend; i++)
	{
	  if(arr[i] > vmax[j]) vmax[j] = arr[i];
	  if(arr[i] < vmin[j]) vmin[j] = arr[i];
	}
    }


  *min = vmin[0];
  *max = vmax[0];
  for(j=1; j<num_pieces; j++)
    {
      if(*max < vmax[j]) *max = vmax[j];
      if(*min > vmin[j]) *min = vmin[j];
    }

  cart_free(vmax);
  cart_free(vmin);
}


void rtuCopyArraysInt(int *dest, int *src, int size)
{
  int i;

  /*
  // Hard-code memcpy for now, but in general we need to check whether
  // doing an OpenMP-parallel loop is faster.
  */
  if(1)
    {
      memcpy(dest,src,sizeof(int)*size);
    }
  else
    {
#pragma omp parallel for default(none), private(i), shared(dest,src,size)
      for(i=0; i<size; i++)
	{
	  dest[i] = src[i];
	}
    }
}


void rtuCopyArraysFloat(float *dest, float *src, int size)
{
  int i;

  /*
  // Hard-code memcpy for now, but in general we need to check whether
  // doing an OpenMP-parallel loop is faster.
  */
  if(1)
    {
      memcpy(dest,src,sizeof(float)*size);
    }
  else
    {
#pragma omp parallel for default(none), private(i), shared(dest,src,size)
      for(i=0; i<size; i++)
	{
	  dest[i] = src[i];
	}
    }
}


#ifdef DEBUG
void rtuCheckGlobalValue(int val, char *name)
{
  int gMin, gMax;

  MPI_Allreduce(&val,&gMin,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&val,&gMax,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

  if(val!=gMin || val!=gMax)
    {
      cart_error("Incorrect global value %s: %d (min: %d, max:%d)",name,val,gMin,gMax);
    }
}
#endif


float** rtupUniformGrid_Sample(int level, int nbin[3], double bb[6], int nvar, int *varid);

int rtuWriteIfritFile(int level, int *nbinIn, double *bbIn, int nvars, int *varid, const char *filename)
{
  int i, ntemp, nbin[3];
  double bb[6];
  float **vars;
  FILE *F;

  cart_assert(nvars > 0);
  for(i=0; i<nDim; i++)
    {
      cart_assert(nbinIn[i] > 0);
      cart_assert(bbIn[2*i+1] > bbIn[2*i]);

      nbin[i] = nbinIn[i];
      bb[2*i+0] = bbIn[2*i+0];
      bb[2*i+1] = bbIn[2*i+1];
    }

  for(i=nDim; i<3; i++)
    {
      nbin[i] = 1;
      bb[2*i+0] = 0.0;
      bb[2*i+1] = num_grid;
    }

  vars = rtupUniformGrid_Sample(level,nbin,bb,nvars,varid);
 
  if(local_proc_id == MASTER_NODE)
    {
      F = fopen(filename,"w");
      if(F == 0)
	{
	  cart_debug("Unable to open file %s for writing.",filename);
	  for(i=0; i<nvars; i++) cart_free(vars[i]);
	  cart_free(vars);
	  return 1;
	}
  
      ntemp = 12;
      fwrite(&ntemp,4,1,F);
      fwrite(nbin+0,4,1,F);
      fwrite(nbin+1,4,1,F);
      fwrite(nbin+2,4,1,F);
      fwrite(&ntemp,4,1,F);
      ntemp = 4*nbin[0]*nbin[1]*nbin[2];
      for(i=0; i<nvars; i++)
	{
	  fwrite(&ntemp,4,1,F); fwrite(vars[i],4,nbin[0]*nbin[1]*nbin[2],F); fwrite(&ntemp,4,1,F);
	}
      fclose(F);

      for(i=0; i<nvars; i++) cart_free(vars[i]);
      cart_free(vars);
    }

   return 0;
}


/*
// ************************************************
//
//              PRIVATE FUNCTIONS
// (should not be called from outside of this file)
//
// ************************************************
*/
long rtupUniformGrid_GetSize(int level, int nbin[3], double bb[6])
{
  int i, j, k, cell;
  long n;
  double pos[3];

  /* count number of buffer cells */
  n = 0;
  for(k=0; k<nbin[2]; k++)
    {
      pos[2] = bb[4] + (bb[5]-bb[4])*(k+0.5)/nbin[2];
      for(j=0; j<nbin[1]; j++)
	{
	  pos[1] = bb[2] + (bb[3]-bb[2])*(j+0.5)/nbin[1];
	  for(i=0; i<nbin[0]; i++)
	    {
	      pos[0] = bb[0] + (bb[1]-bb[0])*(i+0.5)/nbin[0];
	      cell = cell_find_position_above_level(level,pos);
	      if(cell!=-1 && root_cell_type(root_cell_sfc_index(cell_parent_root_cell(cell)))==1)
		{
		  n++;
		}
	    }
	}
    }
  
  return n;
}


void rtupUniformGrid_FillData(int level, int nbin[3], double bb[6], int nvars, int *varid, float **buf, long *loc)
{
  int i, j, k, var, cell;
  long offset, l;
  double pos[3];

  l = 0;
  for(k=0; k<nbin[2]; k++)
    {
      pos[2] = bb[4] + (bb[5]-bb[4])*(k+0.5)/nbin[2];
      for(j=0; j<nbin[1]; j++)
	{
	  pos[1] = bb[2] + (bb[3]-bb[2])*(j+0.5)/nbin[1];
	  offset = nbin[1]*(j+nbin[2]*k);
	  for(i=0; i<nbin[0]; i++)
	    {
	      pos[0] = bb[0] + (bb[1]-bb[0])*(i+0.5)/nbin[0];
	      cell = cell_find_position_above_level(level,pos);
	      if(cell!=-1 && root_cell_type(root_cell_sfc_index(cell_parent_root_cell(cell)))==1)
		{
		  for(var=0; var<nvars; var++)
		    {
		      if(varid[var]>=0 && varid[var]<num_vars)
			{
			  buf[var][l] = cell_var(cell,varid[var]);
			}
		      else if(varid[var]>=RTU_FRACTION && varid[var]<RTU_FRACTION+num_vars)
			{
			  buf[var][l] = cell_var(cell,varid[var]-RTU_FRACTION)/cell_gas_density(cell);
			}
		      else switch(varid[var])
			{
#ifdef RADIATIVE_TRANSFER
			case RTU_GAS_TEMPERATURE:
			  {
			    buf[var][l] = rtTemInK(cell);
			    break;
			  }
#endif  // RADIATIVE_TRANSFER
			case RTU_CELL_LEVEL:
			  {
			    buf[var][l] = cell_level(cell);
			    break;
			  }
			case RTU_LOCAL_PROC:
			  {
			    buf[var][l] = local_proc_id;
			    break;
			  }
			default:
			  {
			    buf[var][l] = 0.0;
			  }
			}
		    }
		  loc[l] = i + offset;
		  l++;
		}
	    }
	}
    }
}


float** rtupUniformGrid_Sample(int level, int nbin[3], double bb[6], int nvars, int *varid)
{
  int i, ip, done;
  long l, ncells;
  float **vars, **buf;
  long *loc;

  buf = cart_alloc(nvars*sizeof(float*));

  if(local_proc_id == MASTER_NODE)
    {
      vars = cart_alloc(nvars*sizeof(float*));
      for(i=0; i<nvars; i++)
	{
	  vars[i] = cart_alloc(nbin[0]*nbin[1]*nbin[2]*sizeof(float));
	}

      for(ip=0; ip<num_procs; ip++)
	{
	  /*
	  //  Measure buffer size
	  */
	  if(ip == 0)
	    {
	      ncells = rtupUniformGrid_GetSize(level,nbin,bb);
	    }
	  else
	    {
	      MPI_Recv(&ncells,1,MPI_LONG,ip,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }

	  /*
	  //  Allocate buffers
	  */
	  for(i=0; i<nvars; i++) buf[i] = cart_alloc(ncells*sizeof(float));
	  loc = cart_alloc(ncells*sizeof(long));

	  /*
	  //  Fill/transfer buffers
	  */
	  if(ip == 0)
	    {
	      rtupUniformGrid_FillData(level,nbin,bb,nvars,varid,buf,loc);
	    }
	  else
	    {
	      MPI_Recv(loc,ncells,MPI_LONG,ip,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      for(i=0; i<nvars; i++) MPI_Recv(buf[i],ncells,MPI_FLOAT,ip,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	    }

	  /*
	  //  Fill variable arrays
	  */
	  for(i=0; i<nvars; i++)
	    {
	      for(l=0; l<ncells; l++)
		{
		  vars[i][loc[l]] = buf[i][l];
		}
	    }
	  
	  /*
	  //  Free buffers
	  */
	  for(i=0; i<nvars; i++) cart_free(buf[i]);
	  cart_free(loc);
 	}

    }
  else
    {
      vars = 0;

      /*
      //  Measure & transfer buffer size
      */
      ncells = rtupUniformGrid_GetSize(level,nbin,bb);
      MPI_Send(&ncells,1,MPI_LONG,MASTER_NODE,0,MPI_COMM_WORLD);

      /*
      //  Allocate buffers
      */
      for(i=0; i<nvars; i++) buf[i] = cart_alloc(ncells*sizeof(float));
      loc = cart_alloc(ncells*sizeof(long));

      /*
      //  Fill & transfer buffers
      */
      rtupUniformGrid_FillData(level,nbin,bb,nvars,varid,buf,loc);

      MPI_Send(loc,ncells,MPI_LONG,MASTER_NODE,0,MPI_COMM_WORLD);
      for(i=0; i<nvars; i++) MPI_Send(buf[i],ncells,MPI_FLOAT,MASTER_NODE,0,MPI_COMM_WORLD);

      /*
      //  Free buffers
      */
      for(i=0; i<nvars; i++) cart_free(buf[i]);
      cart_free(loc);
    }

  cart_free(buf);

  return vars;
}


/*
C
C **************************************************************
C
C   General routines - used in RT block, but not requiring RT 
C   include files. 
C
C **************************************************************
C
C  Compute the average of an array component
C  Attempt to fold loops over the top level and lower levels into one
C  to avoid code duplication
C
      subroutine ArrAvg(mode,ncom,arr,iarr,vAvg)
C
C  mode = 0: volume weighted
C  mode = 1: mass weighted
C  mode = 2: volume weighted average of arr(iarr,*)/hvar(1,*)
C  mode = 3: mass weighted average of arr(iarr,*)/hvar(1,*)
C
      include 'a_tree.h'
      include 'a_control.h'
      dimension arr(ncom,*)
      double precision sum1k(NCPUs)
      double precision sum2k(NCPUs)

      if(iarr.le.0 .or. iarr.gt.ncom) then
         write(0,*) 'Error in ArrAvg: incorrect iarr, ', iarr
         stop
      endif
      
      if(mode.lt.0 .or. mode.gt.3) then
         write(0,*) 'Error in ArrAvg: incorrect mode, ', mode
         stop
      endif
      
      call Get_MaxLevelNow()

      sum1 = 0.0
      sum2 = 0.0

      do Level = MinLevel, MaxLevelNow

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(sum1k,iOctCh,arr,iarr,mode,iSelect,nLevel,maxChild,icd,
C$OMP+ Level,hvar,sum2k)
         DO LL=1,NCPUs
            sum1k(LL) = 0.0d0
            sum2k(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild 
                  i1 = icell + iChild
                  if ( iOctCh(i1) .eq. nil ) then ! only leaves 
                     
                     if(mode .eq. 0) then
                        sum1k(LL) = sum1k(LL) + arr(iarr,i1)
                     endif

                     if(mode .eq. 1) then
                        sum1k(LL) = sum1k(LL) + arr(iarr,i1)*hvar(1,i1)
                        sum2k(LL) = sum2k(LL) + hvar(1,i1)
                     endif
                        
                     if(mode .eq. 2) then
                        sum1k(LL) = sum1k(LL) + arr(iarr,i1)/hvar(1,i1)
                     endif

                     if(mode .eq. 3) then
                        sum1k(LL) = sum1k(LL) + arr(iarr,i1)
                        sum2k(LL) = sum2k(LL) + hvar(1,i1)
                     endif
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         do LL=1,NCPUs
            sum1 = sum1 + sum1k(LL)*CellSize(Level)**3
            sum2 = sum2 + sum2k(LL)*CellSize(Level)**3
         enddo

      enddo

      if(mode.eq.0 .or. mode.eq.2) sum2 = ncell0

      vAvg = sum1/sum2

      return
      end
C
C  Compute the masked average of an array component
C  Attempt to fold loops over the top level and lower levels into one
C  to avoid code duplication
C
      subroutine MaskedArrAvg(mode,ncom,arr,iarr,vAvg,maskBits)
C
C  mode = 0: volume weighted
C  mode = 1: mass weighted
C  mode = 2: volume weighted average of arr(iarr,*)/hvar(1,*)
C  mode = 3: mass weighted average of arr(iarr,*)/hvar(1,*)
C
      include 'a_tree.h'
      include 'a_mask.h'
      include 'a_control.h'
      dimension arr(ncom,*)
      double precision sum1k(NCPUs)
      double precision sum2k(NCPUs)

      if(.not. lUseMask) then
         call ArrAvg(mode,ncom,arr,iarr,vAvg)
         return
      endif

      if(iarr.le.0 .or. iarr.gt.ncom) then
         write(0,*) 'Error in ArrAvg: incorrect iarr, ', iarr
         stop
      endif
      
      if(mode.lt.0 .or. mode.gt.3) then
         write(0,*) 'Error in ArrAvg: incorrect mode, ', mode
         stop
      endif
      
      call Get_MaxLevelNow()

      sum1 = 0.0
      sum2 = 0.0

      do Level = MinLevel, MaxLevelNow

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(sum1k,iOctCh,arr,iarr,mode,iSelect,nLevel,maxChild,icd,
C$OMP+ Level,hvar,sum2k,iMask,maskBits)
         DO LL=1,NCPUs
            sum1k(LL) = 0.0d0
            sum2k(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild 
                  i1 = icell + iChild
                  if(iOctCh(i1).eq.nil .and. 
     .                 and(iMask(i1),maskBits).ne.0) then ! only leaves 
                     
                     if(mode .eq. 0) then
                        sum1k(LL) = sum1k(LL) + arr(iarr,i1)
                     endif

                     if(mode .eq. 1) then
                        sum1k(LL) = sum1k(LL) + arr(iarr,i1)*hvar(1,i1)
                        sum2k(LL) = sum2k(LL) + hvar(1,i1)
                     endif
                        
                     if(mode .eq. 2) then
                        sum1k(LL) = sum1k(LL) + arr(iarr,i1)/hvar(1,i1)
                     endif

                     if(mode .eq. 3) then
                        sum1k(LL) = sum1k(LL) + arr(iarr,i1)
                        sum2k(LL) = sum2k(LL) + hvar(1,i1)
                     endif
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         do LL=1,NCPUs
            sum1 = sum1 + sum1k(LL)*CellSize(Level)**3
            sum2 = sum2 + sum2k(LL)*CellSize(Level)**3
         enddo

      enddo

      if(mode.eq.0 .or. mode.eq.2) call MaskedVolume(sum2,maskBits)

      if(sum2 .gt. 0.5) then
         vAvg = sum1/sum2
      else
         vAvg = 0.0
      endif

      return
      end
C
C  Compute the volume of the masked part of the mesh
C
      subroutine MaskedVolume(vol,maskBits)
      include 'a_tree.h'
      include 'a_mask.h'
      include 'a_control.h'
      double precision sum(NCPUs)

      if(.not. lUseMask) then
         vol = ncell0
         return
      endif

      call Get_MaxLevelNow()

      icd = (ncell0+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic1,ic2),
C$OMP+SHARED(sum,iOctCh,iMask,maskBits,icd)
      do LL=1,NCPUs
         sum(LL) = 0.0
         ic1 = icd*(LL-1) + 1
         ic2 = icd*LL
         if(ic2 .gt. ncell0) ic2 = ncell0
         do i1=ic1,ic2
            if(and(iMask(i1),maskBits) .ne. 0) then
               sum(LL) = sum(LL) + 1.0
            endif
         enddo
      enddo

      vol = 0.0
      do LL=1,NCPUs
         vol = vol + sum(LL)
      enddo

      return
      end
C
C  Compute the maximum and minimum of an array component
C  Attempt to fold loops over the top level and lower levels into one
C  to avoid code duplication
C
      subroutine ArrMaxMin(ncom,arr,iarr,vMax,vMin)
C
      include 'a_tree.h'
      include 'a_control.h'
      dimension arr(ncom,*)
      dimension vmaxk(NCPUs), vmink(NCPUs)

      if(iarr.le.0 .or. iarr.gt.ncom) then
         write(0,*) 'Error in ArrMaxMin: incorrect iarr, ', iarr
         stop
      endif
      
      call Get_MaxLevelNow()

      vMax = arr(iarr,1)
      vMin = arr(iarr,1)
      do LL=1,NCPUs
         vmaxk(LL) = vMax
         vmink(LL) = vMin
      enddo

      do Level = MinLevel, MaxLevelNow

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(iOctCh,arr,iarr,iSelect,nLevel,maxChild,icd,
C$OMP+ Level,vmaxk,vmink)
         DO LL=1,NCPUs
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild 
                  i1 = icell + iChild
                  if ( iOctCh(i1) .eq. nil ) then ! only leaves 
                     
                     vmaxk(LL) = max(vmaxk(LL),arr(iarr,i1))
                     vmink(LL) = min(vmink(LL),arr(iarr,i1))
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         do LL=1,NCPUs
            vMax = max(vMax,vmaxk(LL))
            vMin = min(vMin,vmink(LL))
         enddo

      enddo

      return
      end
C
C  Compute the maxium and minimum of an array component, but only
C  for cells belonging to a given level (all cells, not only leaves)
C  Attempt to fold loops over the top level and lower levels into one
C  to avoid code duplication
C
      subroutine LevArrMaxMin(Level,ncom,arr,iarr,vMax,vMin)
C
      include 'a_tree.h'
      include 'a_control.h'
      dimension vmaxk(NCPUs), vmink(NCPUs)
      dimension arr(ncom,*)

      if(iarr.le.0 .or. iarr.gt.ncom) then
         write(0,*) 'Error in LevArrMaxMin: incorrect iarr, ', iarr
         stop
      endif
      
      call Get_MaxLevelNow()

      if(Level.lt.MinLevel .or. Level.gt.MaxLevelNow) then
         write(0,*) 'Error in LevArrMaxMin: incorrect Level, ', Level
         stop
      endif
      
      vMax = -1.0e35
      vMin =  1.0e35
      do LL=1,NCPUs
         vmaxk(LL) = vMax
         vmink(LL) = vMin
      enddo

      if(Level .eq. MinLevel) then
         nLevel = ncell0
         maxChild = 1
      else
         nLevel = iNOLL(Level)
         call Select_Cells(Level,nLevel) 
         maxChild = 8
      endif

      icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(iOctCh,arr,iarr,iSelect,nLevel,maxChild,icd,
C$OMP+ Level,vmaxk,vmink)
      DO LL=1,NCPUs
         ic1 = icd*(LL-1) + 1
         ic2 = icd*LL
         if(ic2 .gt. nLevel) ic2 = nLevel
         do ic=ic1,ic2

            if(Level .eq. MinLevel) then
               icell = ic - 1
            else
               icell = iSelect(ic) - 1
            endif

            do iChild = 1,maxChild 
               i1 = icell + iChild
                     
               vmaxk(LL) = max(vmaxk(LL),arr(iarr,i1))
               vmink(LL) = min(vmink(LL),arr(iarr,i1))
                        
            enddo
         enddo
      ENDDO

      do LL=1,NCPUs
         vMax = max(vMax,vmaxk(LL))
         vMin = min(vMin,vmink(LL))
      enddo

      return
      end
C
C  Compute the average of a product of two array components
C  Attempt to fold loops over the top level and lower levels into one
C  to avoid code duplication
C
      subroutine Arr2Avg(ncom1,arr1,iarr1,ncom2,arr2,iarr2,vAvg)
C
      include 'a_tree.h'
      include 'a_control.h'
      dimension arr1(ncom1,*), arr2(ncom2,*)
      double precision sumk(NCPUs)

      if(iarr1.le.0 .or. iarr1.gt.ncom1) then
         write(0,*) 'Error in Arr2Avg: incorrect iarr1, ', iarr1
         stop
      endif
      
      if(iarr2.le.0 .or. iarr2.gt.ncom2) then
         write(0,*) 'Error in Arr2Avg: incorrect iarr2, ', iarr2
         stop
      endif
      
      call Get_MaxLevelNow()
      
      sum = 0.0

      do Level = MinLevel, MaxLevelNow

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(sumk,iOctCh,arr1,arr2,iSelect,nLevel,maxChild,icd,
C$OMP+ Level,iarr1,iarr2)
         DO LL=1,NCPUs
            sumk(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild 
                  i1 = icell + iChild
                  if ( iOctCh(i1) .eq. nil ) then ! only leaves 
                     
                     sumk(LL) = sumk(LL) + 
     .                    arr1(iarr1,i1)*arr2(iarr2,i1)
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         do LL=1,NCPUs
            sum = sum + sumk(LL)*CellSize(Level)**3
         enddo

      enddo

      vAvg = sum/ncell0

      return
      end
C
C  Compute the masked average of a product of two array components
C  Attempt to fold loops over the top level and lower levels into one
C  to avoid code duplication
C
      subroutine MaskedArr2Avg(ncom1,arr1,iarr1,ncom2,arr2,iarr2,
     .     vAvg,maskBits)
C
      include 'a_tree.h'
      include 'a_mask.h'
      include 'a_control.h'
      dimension arr1(ncom1,*), arr2(ncom2,*)
      double precision sumk(NCPUs)

      if(.not. lUseMask) then
         call Arr2Avg(ncom1,arr1,iarr1,ncom2,arr2,iarr2,vAvg)
         return
      endif

      if(iarr1.le.0 .or. iarr1.gt.ncom1) then
         write(0,*) 'Error in Arr2Avg: incorrect iarr1, ', iarr1
         stop
      endif
      
      if(iarr2.le.0 .or. iarr2.gt.ncom2) then
         write(0,*) 'Error in Arr2Avg: incorrect iarr2, ', iarr2
         stop
      endif
      
      call Get_MaxLevelNow()
      
      sum = 0.0

      do Level = MinLevel, MaxLevelNow

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(sumk,iOctCh,arr1,arr2,iSelect,nLevel,maxChild,icd,
C$OMP+ Level,iarr1,iarr2,iMask,maskBits)
         DO LL=1,NCPUs
            sumk(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild 
                  i1 = icell + iChild
                  if(iOctCh(i1).eq.nil .and.
     .                 and(iMask(i1),maskBits).ne.0) then ! only leaves 
                     
                      sumk(LL) = sumk(LL) + 
     .                    arr1(iarr1,i1)*arr2(iarr2,i1)
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         do LL=1,NCPUs
            sum = sum + sumk(LL)*CellSize(Level)**3
         enddo

      enddo

      call MaskedVolume(vol,maskBits)

      if(vol .gt. 0.5) then
         vAvg = sum/vol
      else
         vAvg = 0.0
      endif

      return
      end
C
C  Turn the density into mass in a cell and back - needed for
C  solving Poisson equation - but only for a range of levels
C
      subroutine Reformat_ArrMassToDensity(mode,ncom,arr,iarr,
     .     MinModify,MaxModify)
C
C  mode = 0:  format density to mass
C  mode = 1:  format mass to density
C
      include 'a_tree.h'
      include 'a_control.h'
      dimension arr(ncom,*)

      if(iarr.le.0 .or. iarr.gt.ncom) then
         write(0,*) 'Error in Reformat_ArrMassToDensity: ',
     .        'incorrect iarr, ', iarr
         stop
      endif
      
      if(mode.lt.0 .or. mode.gt.1) then
         write(0,*) 'Error in Reformat_ArrMassToDensity: ',
     .        'incorrect mode, ', mode
         stop
      endif
      
      call Get_MaxLevelNow()

      LevelStart = max(Minlevel,MinModify)
      LevelEnd = min(MaxModify,MaxLevelNow)

      do Level = LevelStart, LevelEnd

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         if(mode .eq. 0) then
            formatFac = CellSize(Level)**3
         else
            formatFac = 1.0/CellSize(Level)**3
         endif

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild),
C$OMP+SHARED(sumk,iOctCh,arr,iarr,iSelect,nLevel,maxChild,
C$OMP+ Level,formatFac)
         do ic = 1,nLevel

            if(Level .eq. MinLevel) then
               icell = ic - 1
            else
               icell = iSelect(ic) - 1
            endif

            do iChild = 1,maxChild 
               i1 = icell + iChild
                  
               arr(iarr,i1) = formatFac*arr(iarr,i1)
               
            enddo
         enddo

      enddo

      return
      end
C
C  Project the value from the parent to the child with interpolation
C
      function Project1(iCell,iChild,iarr,ncom,arr)  
C
C  a-la Pyramide
C
      include 'a_tree.h'
      include 'a_control.h'
C
      integer iCell, iChild, iarr
      dimension arr(ncom,*)
C
      dimension iPyr(nchild,3)   ! interpolation pyramid vertices 
c
      data iPyr / 1, 2, 1, 2, 1, 2, 1, 2,   
     &            3, 3, 4, 4, 3, 3, 4, 4,
     &            5, 5, 5, 5, 6, 6, 6, 6  /

      save iPyr
c
        iNb1   = iNb(iCell,iPyr(iChild,1))
        iNb2   = iNb(iCell,iPyr(iChild,2))
        iNb3   = iNb(iCell,iPyr(iChild,3))

        Project1 =  0.1*arr(iarr,iCell) + 
     &                      0.3*arr(iarr,iNb1) +
     &                      0.3*arr(iarr,iNb2) +         
     &                      0.3*arr(iarr,iNb3) 
c
      return
      end
C
C  Project the value from the parent to the child with interpolation
C
      function Project2(iCell,iChild,iarr,ncom,arr)  
C
      include 'a_tree.h'
      include 'a_control.h'
C     
      integer iCell, iChild, iarr
      dimension arr(ncom,*)
C
      dimension pCh(3,nchild)
      dimension nb(6)
C
      data pCh/ 
     .     -0.25, -0.25, -0.25,
     .      0.25, -0.25, -0.25,
     .     -0.25,  0.25, -0.25,
     .      0.25,  0.25, -0.25,
     .     -0.25, -0.25,  0.25,
     .      0.25, -0.25,  0.25,
     .     -0.25,  0.25,  0.25,
     .      0.25,  0.25,  0.25 /

      call iNbAll(iCell,nb)

      a = arr(iarr,iCell)

      b1 = 0.5*(arr(iarr,nb(2))-arr(iarr,nb(1)))
      b2 = 0.5*(arr(iarr,nb(4))-arr(iarr,nb(3)))
      b3 = 0.5*(arr(iarr,nb(6))-arr(iarr,nb(5)))

      c1 = arr(iarr,nb(2)) + arr(iarr,nb(1)) - 2*a
      c2 = arr(iarr,nb(4)) + arr(iarr,nb(3)) - 2*a
      c3 = arr(iarr,nb(6)) + arr(iarr,nb(5)) - 2*a
      
      Project2 = a + 
     .     b1*pCh(1,iChild) + b2*pCh(2,iChild) + b3*pCh(3,iChild) + 
     .     0.5*(c1*pCh(1,iChild)**2+c2*pCh(2,iChild)**2+ 
     .     c3*pCh(3,iChild)**2)

      return
      end
C
C
C
      subroutine UpdateArrAvg(Level1,Level2,ncom,arr,iarr,
     .     arrAvg,arrAvgLevel)
C
C  Update the average of an array without recomputing the average of the whole
C  array.
C
      include 'a_tree.h'
      include 'a_control.h'
      dimension arr(ncom,*), arrAvgLevel(MinLevel:MaxLevel)
      double precision sumk(NCPUs)

      call Get_MaxLevelNow()

      do Level = max(MinLevel,Level1),min(MaxLevelNow,Level2)

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(sumk,iOctCh,iSelect,nLevel,maxChild,Level,arr,iarr,icd)
         DO LL=1,NCPUs
            sumk(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild
                  i1 = icell + iChild
                  if ( iOctCh(i1) .eq. nil ) then ! only leaves 
                     
                     sumk(LL) = sumk(LL) + arr(iarr,i1)
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         sum = 0.0
         do LL=1,NCPUs
            sum = sum + sumk(LL)*CellSize(Level)**3
         enddo
         arrAvgLevel(Level) = sum

      enddo

      arrAvg = 0.0
      do Level=MinLevel,MaxLevelNow
         arrAvg = arrAvg + arrAvgLevel(Level)
      enddo
      arrAvg = arrAvg/ncell0

      return
      end
C
C
C
      subroutine UpdateMaskedArrAvg(Level1,Level2,ncom,arr,iarr,
     .     arrAvg,arrAvgLevel,maskBits)
C
C  Update the average source without recomputing the average of the whole
C  array.
C
      include 'a_tree.h'
      include 'a_mask.h'
      include 'a_control.h'
      dimension arr(ncom,*), arrAvgLevel(MinLevel:MaxLevel)
      double precision sumk(NCPUs)

      if(.not. lUseMask) then
         call UpdateArrAvg(Level1,Level2,ncom,arr,iarr,arrAvg,
     .        arrAvgLevel)
         return
      endif

      call Get_MaxLevelNow()

      do Level = max(MinLevel,Level1),min(MaxLevelNow,Level2)

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(sumk,iOctCh,iSelect,nLevel,maxChild,Level,arr,iarr,
C$OMP+ iMask,maskBits,icd)
         DO LL=1,NCPUs
            sumk(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild
                  i1 = icell + iChild
                  if(iOctCh(i1).eq.nil .and.
     .                 and(iMask(i1),maskBits).ne.0) then ! only leaves

                      sumk(LL) = sumk(LL) + arr(iarr,i1)

                  endif         ! end if iOctCh
               enddo
            enddo
         ENDDO

         sum = 0.0
         do LL=1,NCPUs
            sum = sum + sumk(LL)*CellSize(Level)**3
         enddo
         arrAvgLevel(Level) = sum

      enddo

      call MaskedVolume(sum2,iMaskRTBits)

      arrAvg = 0.0
      do Level=MinLevel,MaxLevelNow
         arrAvg = arrAvg + arrAvgLevel(Level)
      enddo

      if(sum2 .gt. 0.5) then
         arrAvg = arrAvg/sum2
      else
         arrAvg = 0.0
      endif

      return
      end
C
C
C
      subroutine UpdateArr2Avg(Level1,Level2,ncom1,arr1,iarr1,
     .     ncom2,arr2,iarr2,arrAvg,arrAvgLevel)
C
C  Update the average of an array without recomputing the average of the whole
C  array.
C
      include 'a_tree.h'
      include 'a_control.h'
      dimension arr1(ncom1,*), arr2(ncom2,*), 
     .     arrAvgLevel(MinLevel:MaxLevel)
      double precision sumk(NCPUs)

      call Get_MaxLevelNow()

      do Level = max(MinLevel,Level1),min(MaxLevelNow,Level2)

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(sumk,iOctCh,iSelect,nLevel,maxChild,Level,arr1,iarr1,
C$OMP+ arr2,iarr2,icd)
         DO LL=1,NCPUs
            sumk(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild
                  i1 = icell + iChild
                  if ( iOctCh(i1) .eq. nil ) then ! only leaves 
                     
                     sumk(LL) = sumk(LL) + arr1(iarr1,i1)*arr2(iarr2,i1)
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         sum = 0.0
         do LL=1,NCPUs
            sum = sum + sumk(LL)*CellSize(Level)**3
         enddo
         arrAvgLevel(Level) = sum

      enddo

      arrAvg = 0.0
      do Level=MinLevel,MaxLevelNow
         arrAvg = arrAvg + arrAvgLevel(Level)
      enddo
      arrAvg = arrAvg/ncell0

      return
      end
C
C
C
      subroutine UpdateMaskedArr2Avg(Level1,Level2,ncom1,arr1,iarr1,
     .     ncom2,arr2,iarr2,arrAvg,arrAvgLevel,maskBits)
C
C  Update the average source without recomputing the average of the whole
C  array.
C
      include 'a_tree.h'
      include 'a_mask.h'
      include 'a_control.h'
      dimension arr1(ncom1,*), arr2(ncom2,*), 
     .     arrAvgLevel(MinLevel:MaxLevel)
      double precision sumk(NCPUs)

      if(.not. lUseMask) then
         call UpdateArr2Avg(Level1,Level2,ncom1,arr1,iarr1,
     .        ncom2,arr2,iarr2,arrAvg,arrAvgLevel)
         return
      endif

      call Get_MaxLevelNow()

      do Level = max(MinLevel,Level1),min(MaxLevelNow,Level2)

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,ic1,ic2),
C$OMP+SHARED(sumk,iOctCh,iSelect,nLevel,maxChild,Level,arr1,iarr1,
C$OMP+ iMask,maskBits,icd,arr2,iarr2)
         DO LL=1,NCPUs
            sumk(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild
                  i1 = icell + iChild
                  if(iOctCh(i1).eq.nil .and.
     .                 and(iMask(i1),maskBits).ne.0) then ! only leaves

                     sumk(LL) = sumk(LL) + arr1(iarr1,i1)*arr2(iarr2,i1)

                  endif         ! end if iOctCh
               enddo
            enddo
         ENDDO

         sum = 0.0
         do LL=1,NCPUs
            sum = sum + sumk(LL)*CellSize(Level)**3
         enddo
         arrAvgLevel(Level) = sum

      enddo

      call MaskedVolume(sum2,iMaskRTBits)

      arrAvg = 0.0
      do Level=MinLevel,MaxLevelNow
         arrAvg = arrAvg + arrAvgLevel(Level)
      enddo

      if(sum2 .gt. 0.5) then
         arrAvg = arrAvg/sum2
      else
         arrAvg = 0.0
      endif

      return
      end
#ifdef RADIATIVE_TRANSFER
C
C **************************************************************
C
C   General routines - used in RT block and requiring RT include files
C
C **************************************************************
C
C  SplitUpdate for rad. field
C
      subroutine rtSplitUpdate ( Level ) 
      include 'a_tree.h'
      include 'a_control.h'
      integer Level 
      integer i1, i2, iv, jch, icell, idcell, iCh1, iChild, nLevel
C
      real    factor
      parameter ( factor = 0.5**ndim )
C
      IF ( Level .eq. MinLevel ) THEN
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iCh1,iv,jch,iChild)
         do i1 = 1 , ncell0
            iCh1 = iOctCh(i1)
            if ( iCh1 .gt. nil ) then ! only split cells
               iCh1 = iCh1*nchild + mbshift
               do iv = 4 , nvarConservedMax ! 1st child separate to avoid zeroing 
                  var(iv,i1) = var(iv,iCh1)
               enddo            
               do jch = 1 , nchild-1
                  iChild = iCh1 + jch
                  do iv = 4 , nvarConservedMax
                     var(iv,i1) = var(iv,i1) + var(iv,iChild)
                  enddo
               enddo
               do iv = 4 , nvarConservedMax
                  var(iv,i1) = factor * var(iv,i1)
               enddo
            endif
         enddo
      ELSE
         nLevel = iNOLL(Level)
         call Select_Cells ( Level , nLevel ) 
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(icell,i2,idcell,iCh1,iv,jch,iChild)
         do i1 = 1 , nLevel
            icell = iSelect(i1)
            do i2 = 0 , 7
               idcell = icell + i2
               iCh1 = iOctCh(idcell) 
               if ( iCh1 .gt. nil ) then ! only split cells
                  iCh1 = iCh1*nchild + mbshift
                  do iv = 4 , nvarConservedMax ! 1st child separate to avoid zeroing 
                     var(iv,idcell) = var(iv,iCh1)
                  enddo
                  do jch = 1 , nchild-1
                     iChild = iCh1 + jch
                     do iv = 4 , nvarConservedMax
                        var(iv,idcell) = var(iv,idcell) + var(iv,iChild)
                     enddo
                  enddo
                  do iv = 4 , nvarConservedMax
                     var(iv,idcell) = factor * var(iv,idcell)
                  enddo              
               endif
            enddo
         enddo
      ENDIF

      return
      end
C
C  Compute the average of the temperature.
C  Attempt to fold loops over the top level and lower levels into one
C  to avoid code duplication
C
      subroutine TemAvg(mode,tAvg)
C
C  mode = 0: volume weighted
C  mode = 1: mass weighted
C
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_rt_rad.h'
      double precision sum1k(NCPUs)
      double precision sum2k(NCPUs)

      if(mode.lt.0 .or. mode.gt.1) then
         write(0,*) 'Error in TemAvg: incorrect mode, ', mode
         stop
      endif
      
      call Get_MaxLevelNow()

      sum1 = 0.0
      sum2 = 0.0
      Tfactor = (gamma-1)*T_0/aa**2/wmu

      do Level = MinLevel, MaxLevelNow

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,Tgas,ic1,ic2),
C$OMP+SHARED(sum1k,iOctCh,arr,iarr,mode,iSelect,nLevel,maxChild,icd,
C$OMP+ Level,hvar,sum2k,Tfactor)
         DO LL=1,NCPUs
            sum1k(LL) = 0.0d0
            sum2k(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild 
                  i1 = icell + iChild
                  if ( iOctCh(i1) .eq. nil ) then ! only leaves 
                     
                     Tgas = Tfactor*hvar(8,i1)/((XH+XG)*hvar(1,i1)+
     .                    hvar(ivarH2,i1)+hvar(ivarG2,i1)+
     .                    2*hvar(ivarG3,i1))

                     if(mode .eq. 0) then
                        sum1k(LL) = sum1k(LL) + Tgas
                     endif

                     if(mode .eq. 1) then
                        sum1k(LL) = sum1k(LL) + Tgas*hvar(1,i1)
                        sum2k(LL) = sum2k(LL) + hvar(1,i1)
                     endif
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         do LL=1,NCPUs
            sum1 = sum1 + sum1k(LL)*CellSize(Level)**3
            sum2 = sum2 + sum2k(LL)*CellSize(Level)**3
         enddo

      enddo

      if(mode .eq. 0) sum2 = ncell0

      tAvg = sum1/sum2

      return
      end
C
C  Compute the masked average of the temperature.
C  Attempt to fold loops over the top level and lower levels into one
C  to avoid code duplication
C
      subroutine MaskedTemAvg(mode,tAvg,maskBits)
C
C  mode = 0: volume weighted
C  mode = 1: mass weighted
C
      include 'a_tree.h'
      include 'a_mask.h'
      include 'a_control.h'
      include 'a_rt_rad.h'
      double precision sum1k(NCPUs)
      double precision sum2k(NCPUs)

      if(.not. lUseMask) then
         call TemAvg(mode,tAvg)
         return
      endif

      if(mode.lt.0 .or. mode.gt.1) then
         write(0,*) 'Error in TemAvg: incorrect mode, ', mode
         stop
      endif
      
      call Get_MaxLevelNow()

      sum1 = 0.0
      sum2 = 0.0
      Tfactor = (gamma-1)*T_0/aa**2/wmu

      do Level = MinLevel, MaxLevelNow

         if(Level .eq. MinLevel) then
            nLevel = ncell0
            maxChild = 1
         else
            nLevel = iNOLL(Level)
            call Select_Cells(Level,nLevel) 
            maxChild = 8
         endif

         icd = (nLevel+NCPUs-1)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(LL,i1,ic,icell,iChild,Tgas,ic1,ic2),
C$OMP+SHARED(sum1k,iOctCh,arr,iarr,mode,iSelect,nLevel,maxChild,icd,
C$OMP+ Level,hvar,sum2k,Tfactor,iMask,maskBits)
         DO LL=1,NCPUs
            sum1k(LL) = 0.0d0
            sum2k(LL) = 0.0d0
            ic1 = icd*(LL-1) + 1
            ic2 = icd*LL
            if(ic2 .gt. nLevel) ic2 = nLevel
            do ic=ic1,ic2

               if(Level .eq. MinLevel) then
                  icell = ic - 1
               else
                  icell = iSelect(ic) - 1
               endif

               do iChild = 1,maxChild 
                  i1 = icell + iChild
                  if(iOctCh(i1).eq.nil .and. 
     .                 and(iMask(i1),maskBits).ne.0) then ! only leaves 
                     
                     Tgas = Tfactor*hvar(8,i1)/((XH+XG)*hvar(1,i1)+
     .                    hvar(ivarH2,i1)+hvar(ivarG2,i1)+
     .                    2*hvar(ivarG3,i1))

                     if(mode .eq. 0) then
                        sum1k(LL) = sum1k(LL) + Tgas
                     endif

                     if(mode .eq. 1) then
                        sum1k(LL) = sum1k(LL) + Tgas*hvar(1,i1)
                        sum2k(LL) = sum2k(LL) + hvar(1,i1)
                     endif
                        
                  endif         ! end if iOctCh            
               enddo
            enddo
         ENDDO

         do LL=1,NCPUs
            sum1 = sum1 + sum1k(LL)*CellSize(Level)**3
            sum2 = sum2 + sum2k(LL)*CellSize(Level)**3
         enddo

      enddo

      if(mode .eq. 0) then

         call MaskedVolume(sum2,maskBits)

      endif

      if(sum2 .gt. 0.5) then
         tAvg = sum1/sum2
      else
         tAvg = 0.0
      endif

      return
      end
C
C  Compute average stellar density
C
      subroutine ComputeStellarDensity(starAvg)
      include 'a_tree.h'
      include 'a_control.h'
      double precision sump(NCPUs)

      isd = (nsp(nspecies,2)-nsp(nspecies,1)+NCPUs)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(l,is,is1,is2),
C$OMP+SHARED(nsp,nspecies,pw,sump,isd)
      do l=1,NCPUs
         sump(l) = 0.0d0
         is1 = isd*(l-1) + nsp(nspecies,1)
         is2 = is1 + isd - 1
         if(is2 .gt. nsp(nspecies,2)) is2 = nsp(nspecies,2)
         do is=is1,is2
            sump(l) = sump(l) + pw(is)
         enddo
      enddo

      starAvg = 0.0
      do l=1,NCPUs
         starAvg = starAvg + sump(l)
      enddo
      starAvg = starAvg/ncell0

      return
      end
C
C  Compute masked average stellar density
C
      subroutine ComputeMaskedStellarDensity(starAvg,maskBits)
      include 'a_tree.h'
      include 'a_mask.h'
      include 'a_control.h'
      double precision sump(NCPUs)

      if(.not. lUseMask) then
         call ComputeStellarDensity(starAvg)
         return
      endif

      isd = (nsp(nspecies,2)-nsp(nspecies,1)+NCPUs)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(l,is,is1,is2,ic),
C$OMP+SHARED(nsp,nspecies,pw,sump,isd,iMask,maskBits,x,y,z)
      do l=1,NCPUs
         sump(l) = 0.0d0
         is1 = isd*(l-1) + nsp(nspecies,1)
         is2 = is1 + isd - 1
         if(is2 .gt. nsp(nspecies,2)) is2 = nsp(nspecies,2)
         do is=is1,is2
            ic = iFindCell(MinLevel,x(is),y(is),z(is))
            if(and(iMask(ic),maskBits) .ne. 0) then
               sump(l) = sump(l) + pw(is)
            endif
         enddo
      enddo

      starAvg = 0.0
      do l=1,NCPUs
         starAvg = starAvg + sump(l)
      enddo

      call MaskedVolume(starVol,maskBits)

      if(starVol .gt. 0.5) then
         starAvg = starAvg/starVol
      else
         starAvg = 0.0
      endif

      return
      end
C
C  Compute average source density
C
      subroutine ComputeSourceDensity(srcAvg)
      include 'a_tree.h'
      include 'a_control.h'
      include 'a_rt_rad.h'
      double precision sump(NCPUs)

      isd = (nsp(nspecies,2)-nsp(nspecies,1)+NCPUs)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(l,is,is1,is2),
C$OMP+SHARED(nsp,nspecies,pw,sump,isd)
      do l=1,NCPUs
         sump(l) = 0.0d0
#ifdef RT_SOURCES
         is1 = isd*(l-1) + nsp(nspecies,1)
         is2 = is1 + isd - 1
         if(is2 .gt. nsp(nspecies,2)) is2 = nsp(nspecies,2)
         do is=is1,is2
            sump(l) = sump(l) + pw(is)*rtSor(is)
         enddo
#endif
      enddo

      srcAvg = 0.0
      do l=1,NCPUs
         srcAvg = srcAvg + sump(l)
      enddo
      srcAvg = srcAvg/ncell0

      return
      end
C
C  Compute masked average source density
C
      subroutine ComputeMaskedSourceDensity(srcAvg,maskBits)
      include 'a_tree.h'
      include 'a_mask.h'
      include 'a_control.h'
      include 'a_rt_rad.h'
      double precision sump(NCPUs)

      if(.not. lUseMask) then
         call ComputeSourceDensity(srcAvg)
         return
      endif

      isd = (nsp(nspecies,2)-nsp(nspecies,1)+NCPUs)/NCPUs

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(l,is,is1,is2,ic),
C$OMP+SHARED(nsp,nspecies,pw,sump,isd,iMask,maskBits,x,y,z)
      do l=1,NCPUs
         sump(l) = 0.0d0
#ifdef RT_SOURCES
         is1 = isd*(l-1) + nsp(nspecies,1)
         is2 = is1 + isd - 1
         if(is2 .gt. nsp(nspecies,2)) is2 = nsp(nspecies,2)
         do is=is1,is2
            ic = iFindCell(MinLevel,x(is),y(is),z(is))
            if(and(iMask(ic),maskBits) .ne. 0) then
               sump(l) = sump(l) + pw(is)*rtSor(is)
            endif
         enddo
#endif
      enddo

      srcAvg = 0.0
      do l=1,NCPUs
         srcAvg = srcAvg + sump(l)
      enddo

      call MaskedVolume(srcVol,maskBits)

      if(srcVol .gt. 0.5) then
         srcAvg = srcAvg/srcVol
      else
         srcAvg = 0.0
      endif

      return
      end
C
C
C
      subroutine DoNothing(a)
C
C  This subroutine does what it claims to do.
C
      return
      end
*/
