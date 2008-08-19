#include "defs.h"      
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
C
C  Sum-up global contributions
C
      subroutine frtAccumulateGlobalFields(ntot,rGlob)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_global.inc'
      dimension rGlob(ioDim,*)
C
      real*8 dTabH, dTabM, dNAvg, dFail, dCall
      real*8 dCXH1, dCXG1, dCXG2, dDust, dBrem
      real*8 dRAH2, dRAG2, dRAG3, dRBH2, dRBG2, dRBG3
C
C Recombination and bremstrahlung photons that need to be accumulated
C during the CoolOff calls
C
      dCXH1 = 0.0d0
      dCXG1 = 0.0d0
      dCXG2 = 0.0d0
      dRAH2 = 0.0d0
      dRAG2 = 0.0d0
      dRAG3 = 0.0d0
      dRBH2 = 0.0d0
      dRBG2 = 0.0d0
      dRBG3 = 0.0d0
      dDust = 0.0d0

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(ic),
C$OMP+SHARED(ntot,rGlob),
C$OMP+REDUCTION(+: dCXH1,dCXG1,dCXG2,dDust,dRAH2,dRAG2,dRAG3,dRBH2,
C$OMP+ dRBG2,dRBG3)
      do ic=1,ntot
         dCXH1 = dCXH1 + rGlob(ioCXH1,ic)
         dCXG1 = dCXG1 + rGlob(ioCXG1,ic)
         dCXG2 = dCXG2 + rGlob(ioCXG2,ic)
         dRAH2 = dRAH2 + rGlob(ioRAH2,ic)
         dRAG2 = dRAG2 + rGlob(ioRAG2,ic)
         dRAG3 = dRAG3 + rGlob(ioRAG3,ic)
         dRBH2 = dRBH2 + rGlob(ioRBH2,ic)
         dRBG2 = dRBG2 + rGlob(ioRBG2,ic)
         dRBG3 = dRBG3 + rGlob(ioRBG3,ic)
         dDust = dDust + rGlob(ioDust,ic)
      enddo

      gData(ioCXH1) = gData(ioCXH1) + dCXH1
      gData(ioCXG1) = gData(ioCXG1) + dCXG1
      gData(ioCXG2) = gData(ioCXG2) + dCXG2
      gData(ioRAH2) = gData(ioRAH2) + dRAH2
      gData(ioRAG2) = gData(ioRAG2) + dRAG2
      gData(ioRAG3) = gData(ioRAG3) + dRAG3
      gData(ioRBH2) = gData(ioRBH2) + dRBH2
      gData(ioRBG2) = gData(ioRBG2) + dRBG2
      gData(ioRBG3) = gData(ioRBG3) + dRBG3
      gData(ioDust) = gData(ioDust) + dDust
C
C  Not every OMP implementation supports array reduction, so
C  we do this loop separately. 
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(ic,lr,teml,dBrem),
C$OMP+SHARED(ntot,rGlob,Txi,gBrem)
      do lr=1,NRAD

         dBrem = 0.0d0

         do ic=1,ntot
            teml = rGlob(ioTemp,ic)
            if(teml .gt. 1.0e-2*Txi(lr)) then
               dBrem = dBrem + rGlob(ioBrem,ic)*exp(-Txi(lr)/teml)/teml
            endif
         enddo

         gBrem(lr) = gBrem(lr) + dBrem

      enddo
C
C  Diagnostics
C
      dTabH = 0.0d0
      dTabM = 0.0d0
      dNAvg = 0.0d0
      dFail = 0.0d0
      dCall = 0.0d0

      rNMax = 0.0

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(ic),
C$OMP+SHARED(ntot,rGlob),
C$OMP+REDUCTION(+: dTabH,dTabM,dNAvg,dFail,dCall), 
C$OMP+REDUCTION(max: rNMax)
      do ic=1,ntot
         dTabH = dTabH + rGlob(ioTabH,ic)
         dTabM = dTabM + rGlob(ioTabM,ic)
         if(rGlob(ioInfo,ic) .lt. 0.0) then
            dFail = dFail + 1.0d0
            dCall = dCall + 1.0d0
         else if(rGlob(ioInfo,ic) .gt. 0.0) then
            dNAvg = dNAvg + rGlob(ioInfo,ic)
            rNMax = max(rNMax,rGlob(ioInfo,ic))
            dCall = dCall + 1.0d0
         endif
      enddo

      gInfo(infoTabH) = gInfo(infoTabH) + dTabH
      gInfo(infoTabM) = gInfo(infoTabM) + dTabM
      gInfo(infoNAvg) = gInfo(infoNAvg) + dNAvg
      gInfo(infoNMax) = max(gInfo(infoNMax),rNMax)
      gInfo(infoFail) = gInfo(infoFail) + dFail
      gInfo(infoCall) = gInfo(infoCall) + dCall

      return
      end
C
C
C
      subroutine frtResetGlobalFields
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_global.inc'

      do j=1,iglobBDim
         gBrem(j) = 0.0
      enddo
      do j=1,iglobDDim
         gData(j) = 0.0
      enddo
      do j=1,iglobIDim
         gInfo(j) = 0.0
      enddo

      return
      end
C
C
C
      function frtQueryGlobalFieldSize()
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_global.inc'
      integer frtQueryGlobalFieldSize

      frtQueryGlobalFieldSize = iglobBDim + iglobDDim + iglobIDim

      return
      end
C
C
C
      subroutine frtQueryGlobalFieldData(data)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_global.inc'
      dimension data(*)

      do j=1,iglobBDim
         data(j) = gBrem(j)
      enddo
      do j=1,iglobDDim
         data(iglobBDim+j) = gData(j)
      enddo
      do j=1,iglobIDim
         data(iglobBDim+iglobDDim+j) = gInfo(j)
      enddo

      return
      end
C
C
C
      subroutine frtUpdateGlobalFieldData(data)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_global.inc'
      dimension data(*)

      do j=1,iglobBDim
         gBrem(j) = data(j)
      enddo
      do j=1,iglobDDim
         gData(j) = data(iglobBDim+j)
      enddo
      do j=1,iglobIDim
         gInfo(j) = data(iglobBDim+iglobDDim+j)
      enddo

      return
      end

#endif

      subroutine frtDummyGlobal(i)
      i = 0
      end
