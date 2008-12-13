#include "defs.h"      
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
C
C  Add one global contributions to a specific buffer
C
      subroutine frtIncrementGlobalFields(iBuf,avg,itab,info)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_global.inc'
      dimension avg(iavgDim), itab(2)
C
      real*8 dBremBuf(NRAD,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      real*8 dDataBuf(iavgMax1,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      real*8 dInfoBuf(infoDim,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      common/RT_GLBBUF/ dBremBuf, dDataBuf, dInfoBuf

      if(iBuf .lt. 1) then
         write(0,*) 'RT_Global: fatal error, ',
     .        'invalid buffer ID; buffer ID must be > 0.'
         stop
      endif

      if(iBuf .gt. RT_PARALLEL_NUM_OPENMP_BUFFERS) then
         write(0,*) 'RT_Global: fatal error, ',
     .        'not enough global buffers. Set the variable ',
     .        'RT_PARALLEL_NUM_OPENMP_BUFFERS in rt_defs.h ',
     .        'to at least ', iBuf
         stop
      endif

#ifndef RT_NO_BACKGROUND
C
C  Bremstrahlung
C
      do lr=1,NRAD
         teml = avg(iavgTemp)
         if(teml .gt. 1.0e-2*Txi(lr)) then
            dBremBuf(lr,iBuf) = dBremBuf(lr,iBuf) + 
     .           avg(iavgBrem)*exp(-Txi(lr)/teml)/teml
         endif
      enddo
C
C Recombination photons and other quantities
C
      do j=1,iavgMax1
         dDataBuf(j,iBuf) = dDataBuf(j,iBuf) + avg(j)
      enddo
#endif

#ifdef RT_OUTPUT
C
C  Diagnostics
C
      dInfoBuf(infoTabH,iBuf) = dInfoBuf(infoTabH,iBuf) + itab(1)
      dInfoBuf(infoTabM,iBuf) = dInfoBuf(infoTabM,iBuf) + itab(2)
      if(info .lt. 0) then
         dInfoBuf(infoFail,iBuf) = dInfoBuf(infoFail,iBuf) + 1.0d0
         dInfoBuf(infoCall,iBuf) = dInfoBuf(infoCall,iBuf) + 1.0d0
      else if(info .gt. 0) then
         dInfoBuf(infoNAvg,iBuf) = dInfoBuf(infoNAvg,iBuf) + info
         dInfoBuf(infoNMax,iBuf) = max(dInfoBuf(infoNMax,iBuf),
     .        dble(info))
         dInfoBuf(infoCall,iBuf) = dInfoBuf(infoCall,iBuf) + 1.0d0
      endif
#endif

      return
      end
C
C Initialize global field buffers
C
      subroutine frtInitGlobalFields
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_global.inc'
C
      real*8 dBremBuf(NRAD,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      real*8 dDataBuf(iavgMax1,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      real*8 dInfoBuf(infoDim,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      common/RT_GLBBUF/ dBremBuf, dDataBuf, dInfoBuf
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(i,j),
C$OMP+SHARED(dBremBuf,dDataBuf,dInfoBuf)
      do i=1,RT_PARALLEL_NUM_OPENMP_BUFFERS
         do j=1,NRAD
            dBremBuf(j,i) = 0.0
         enddo
         do j=1,iavgMax1
            dDataBuf(j,i) = 0.0
         enddo
         do j=1,infoDim
            dInfoBuf(j,i) = 0.0
         enddo
      enddo

      return
      end
C
C  Sum-up global contributions
C
      subroutine frtFinalizeGlobalFields(vol)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_global.inc'
C
      real*8 dBremBuf(NRAD,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      real*8 dDataBuf(iavgMax1,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      real*8 dInfoBuf(infoDim,RT_PARALLEL_NUM_OPENMP_BUFFERS)
      common/RT_GLBBUF/ dBremBuf, dDataBuf, dInfoBuf
C
      real*8 s
      real*8 dBremBufLoc(NRAD)
      real*8 dDataBufLoc(iavgMax1)
      real*8 dInfoBufLoc(infoDim)
      real*8 dBremBufGlob(NRAD)
      real*8 dDataBufGlob(iavgMax1)
      real*8 dInfoBufGlob(infoDim)
C
      weight = vol/tStep
C
C  Bremstrahlung
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr,ib,s),
C$OMP+SHARED(dBremBufLoc,dBremBuf,weight)
      do lr=1,NRAD
         s = 0.0d0
         do ib=1,RT_PARALLEL_NUM_OPENMP_BUFFERS
            s = s + dBremBuf(lr,ib)
         enddo
         dBremBufLoc(lr) = weight*s
      enddo
C
C Recombination photons and other quantities
C
      do j=1,iavgMax1
         s = 0.0d0
         do ib=1,RT_PARALLEL_NUM_OPENMP_BUFFERS
            s = s + dDataBuf(j,ib)
         enddo
         dDataBufLoc(j) = weight*s
      enddo
C
C  Diagnostics
C
      do j=1,infoMax1
         s = 0.0d0
         do ib=1,RT_PARALLEL_NUM_OPENMP_BUFFERS
            s = s + dInfoBuf(j,ib)
         enddo
         dInfoBufLoc(j) = s
      enddo

      s = dInfoBuf(infoNMax,1)
      do ib=2,RT_PARALLEL_NUM_OPENMP_BUFFERS
         s = max(s,dInfoBuf(infoNMax,ib))
      enddo
      dInfoBufLoc(infoNMax) = s
C
C  Final summation (if needed)
C
      call frtGlobalAverage(0,NRAD,dBremBufLoc,dBremBufGlob)
      call frtGlobalAverage(0,iavgMax1,dDataBufLoc,dDataBufGlob)
      call frtGlobalAverage(0,infoMax1,dInfoBufLoc,dInfoBufGlob)
      call frtGlobalAverage(1,1,dInfoBufLoc(infoNMax),
     .     dInfoBufGlob(infoNMax))

      do lr=1,NRAD
         gBrem(lr) = dBremBufGlob(lr)
      enddo

      do j=1,iavgMax1
         gData(j) = dDataBufGlob(j)
      enddo

      do j=1,infoDim
         gInfo(j) = dInfoBufGlob(j)
      enddo

      return
      end

#endif

      subroutine frtDummyGlobal(i)
      i = 0
      end
