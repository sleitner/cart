#include "../config.h"
#ifdef RADIATIVE_TRANSFER
#ifdef RT_TRANSFER
C
C
C
      subroutine frtInitRunTransfer(nfreq)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
C
C  Set sampling points here (except the last point)
C
#ifdef RT_UV
      if(NSETS .ne. 4) then
         write(0,*) 'RT_Transfer: fatal error, ',
     .        'only NSETS=4 is supported.'
         call frtAbort
      endif
#else
      if(NSETS .ne. 3) then
         write(0,*) 'RT_Transfer: fatal error, ',
     .        'only NSETS=3 is supported.'
         call frtAbort
      endif
#endif

      if(nfreq .ne. 2*NSETS) then
         write(0,*) 'RT_Transfer: fatal error, ',
     .        'only nfreq=2*NSETS is supported.'
         call frtAbort
      endif

      LRSETS(1) = lrH1tr
      LRSETS(2) = lrG1tr
      LRSETS(3) = lrG2tr
#ifdef RT_UV
      LRSETS(4) = lrUVtr
#endif

      return
      end
C
C
C
      subroutine frtSetGlobalRadiationFieldsTransfer(rfAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension w(NSETS), rfAvg(2*NSETS)
      include 'frt_inline.inc'
C
      do L=1,NSETS

         lr = LRSETS(L)

         aqbar = facLX*guvEff(lr)
         ajbar = angAvg(lr)

         q = ajbar - aqbar*rfAvg(L)
         if(q .lt. 0.0) q = 0.0
         p = ajbar*rfAvg(NSETS+L)
         if(p .gt. 1.0e-35) then
            w(L) = q/p
         else
            w(L) = 0.0
         endif

      enddo

      tauH1 = QIuni(w(1))
      do lr=1,NRAD
         tauLoc(lr) = tauH1*csH1n(lr)
      enddo

      tauG1 = QIuni(w(2)) - tauH1*csH1n(lrG1tr)
      do lr=lrG1tr,NRAD
         tauLoc(lr) = tauLoc(lr) + tauG1*csG1n(lr)
      enddo

      tauG2 = QIuni(w(3)) - tauG1*csG1n(lrG2tr) - 
     .     tauH1*csH1n(lrG2tr)
      do lr=lrG2tr,NRAD
         tauLoc(lr) = tauLoc(lr) + tauG2*csG2n(lr)
      enddo

#ifdef RT_UV
C
C  Ignoring average absorption by dust
C
#endif

      return
      end
C
C
C
      subroutine frtStepBeginTransfer(abcAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension abcAvg(*)
C
      common /RT_Transfer_Save1/ abcBeg(NSETS)
C
      do i=1,NSETS
         abcBeg(i) = abcAvg(i)
      enddo

#ifdef RT_DEBUG
      if(ioNode) then
         abe1 = abcAvg(1)
         abe2 = abcAvg(2) - abe1*csH1n(lrG1tr)
         abe3 = abcAvg(3) - abe1*csH1n(lrG2tr) - abe2*csG1n(lrG2tr)
         write(6,*) 'RT: <abcH1> = ', abcAvg(1), abe1/csH1tr/uColumn
         write(6,*) 'RT: <abcG1> = ', abcAvg(2), abe2/csG1tr/uColumn
         write(6,*) 'RT: <abcG2> = ', abcAvg(3), abe3/csG2tr/uColumn
      endif
#endif

      do i=1,NSETS
         if(abcAvg(i) .lt. 0.0) then
            write(0,*) 'RT_Transfer: fatal error, ',
     .           'negative global absoprtion: ', i, abcAvg(i)
            call frtAbort
         endif
      enddo

      return
      end
C
C
C
      subroutine frtStepEndTransfer(abcAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension abcAvg(*)
C
      common /RT_Transfer_Save2/ abcEnd(NSETS)
C
      do i=1,NSETS
         abcEnd(i) = abcAvg(i)
      enddo

      do i=1,NSETS
         if(abcAvg(i) .lt. 0.0) then
            write(0,*) 'RT_Transfer: fatal error, ',
     .           'negative global absoprtion: ', i, abcAvg(i)
            call frtAbort
         endif
      enddo

      return
      end
C
C  This is hard-coded for NSETS=3 or 4 case
C
      subroutine frtTransferComputeAverageAbs(axi)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension axi(*)
C
      common /RT_Transfer_Save1/ abcBeg(NSETS)
      common /RT_Transfer_Save2/ abcEnd(NSETS)
C
      do i=1,NSETS
         abcEnd(i) = 0.5*(abcEnd(i)+abcBeg(i))
      enddo

      axiH1 = max(0.0,abcEnd(1)/facLX)
      do lr=lrH1tr,lrmax
         axi(lr) = axi(lr) + axiH1*csH1n(lr)
      enddo

      axiG1 = max(0.0,abcEnd(2)/facLX-axiH1*csH1n(lrG1tr))
      do lr=lrG1tr,lrmax
         axi(lr) = axi(lr) + axiG1*csG1n(lr)
      enddo

      axiG2 = max(0.0,abcEnd(3)/facLX-axiG1*csG1n(lrG2tr)
     .     -axiH1*csH1n(lrG2tr))
      do lr=lrG2tr,lrmax
         axi(lr) = axi(lr) + axiG2*csG2n(lr)
      enddo

#ifdef RT_UV
C
C  Ignoring average absorption by dust
C
#endif

      return
      end
C
C  Do some preprocessing
C
      subroutine frtTransferPackRadiationField(par,y,rawRF0,rawRF1,rf)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension par(*), y(*), rawRF0(*), rawRF1(*), rf(*)
#ifdef RT_TRANSFER_FLUX_CONSERVING
      parameter(wFlux = 0.0, wCell = 0.0)
      dimension tau(NSETS)
#endif

      rf(irfNorm) = rawRF0(1)
C
C  1. Local field
C
      call frtTransferPackHelper(0,par,y,rawRF0,rawRF1,rf)
C
C  2. Global field
C
#ifdef RT_EXTERNAL_BACKGROUND
      call frtTransferPackHelper(1,par,y,rawRF0,rawRF1,rf)
#else
      do j=1,NSETS
         rf(irfOgf+j) = -1.0
      enddo
#endif

#ifdef RT_TRANSFER_FLUX_CONSERVING
      tau(1) = csH1tr*par(iparDelN)
      tau(2) = csG1tr*par(iparDelN)
      tau(3) = csG2tr*par(iparDelN)
#ifdef RT_UV
      tau(4) = frtDustToGas(par(iparZsol))*csUVtr*par(iparDelN)
#endif

      do j=1,NSETS
         q = 1.0/max(ySmall(1+j),y(1+j))
         if(rf(irfOlf+j) .lt. 0.0) then
            rf(irfOlf+j) = -1.0
            rf(irfOlv+j) = 0.0
         else
            rf(irfOlf+j) = max(0.0,rf(irfOlf+j)-wFlux*tau(j))
            rf(irfOlv+j) = wFlux*tau(j)*q
         endif
         if(rf(irfOgf+j) .lt. 0.0) then
            rf(irfOgf+j) = -1.0
            rf(irfOgv+j) = 0.0
         else
            rf(irfOgf+j) = max(0.0,rf(irfOgf+j)-wFlux*tau(j))
            rf(irfOgv+j) = wFlux*tau(j)*q
         endif
         rf(irfOwf+j) = (1-wCell)*tau(j)
         rf(irfOwv+j) = wCell*tau(j)*q
      enddo
#endif

      return
      end
C
C Do some postprocessing
C
#ifdef RT_VARIABLE_RFIELD
      subroutine frtTransferUnPackRadiationField(y,rawRF0,rawRF1,rf)
      dimension y(*), rawRF0(*), rawRF1(*), rf(*)
C
C  1. Local field
C
      call frtTransferUnPackHelper(0,y,rawRF0,rawRF1,rf)
C
C  2. Global field
C
      call frtTransferUnPackHelper(1,y,rawRF0,rawRF1,rf)

      return
      end
#endif /* RT_TRANSFER_FLUX_CONSERVING */
C
C  *******************
C
C    Helpers
C
C  *******************
C
      subroutine frtTransferPackHelper(k,par,y,w0,w,rf)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension par(*), y(*), w0(*), w(*), rf(*)
      include 'frt_inline.inc'
C
C  Works only for the NSETS=3 or 4 case
C
      if(w0(k+1) .gt. 1.0e-35) then

         if(w(k*NSETS+1) .gt. 1.0e-35) then
            t1 = min(1.0,w(k*NSETS+1)/w0(k+1))
            if(k .eq. 0) then
               ta1 = QIloc(t1)
            else
               ta1 = QIuni(t1)
            endif
         else
            ta1 = -1.0
         endif

         if(w(k*NSETS+2) .gt. 1.0e-35) then
            t2 = min(1.0,w(k*NSETS+2)/w0(k+1))
            if(k .eq. 0) then
               ta2 = QIloc(t2)
            else
               ta2 = QIuni(t2)
            endif
            if(ta1 .gt. 0.0) then
               ta2 = ta2 - ta1*csH1n(lrG1tr)
            endif
         else
            ta2 = -1.0
         endif

         if(w(k*NSETS+3) .gt. 1.0e-35) then
            t3 = min(1.0,w(k*NSETS+3)/w0(k+1))
            if(k .eq. 0) then
               ta3 = QIloc(t3)
            else
               ta3 = QIuni(t3)
            endif
            if(ta1 .gt. 0.0) then
               ta3 = max(1.0e-10,ta3-ta1*csH1n(lrG2tr))
            endif
            if(ta2 .gt. 0.0) then
               ta3 = max(1.0e-10,ta3-ta2*csG1n(lrG2tr))
            endif
         else
            ta3 = -1.0
         endif

#ifdef RT_UV
         if(w(k*NSETS+4) .gt. 1.0e-35) then
            t4 = min(1.0,w(k*NSETS+4)/w0(k+1))
            if(k .eq. 0) then
               ta4 = QIloc(t4)
            else
               ta4 = QIuni(t4)
            endif
         else
            ta4 = -1.0
         endif
#endif

         rf(k*NSETS+irfH1lf) = ta1
         rf(k*NSETS+irfG1lf) = ta2
         rf(k*NSETS+irfG2lf) = ta3
#ifdef RT_UV
         rf(k*NSETS+irfUVlf) = ta4
#endif

      else

         rf(k*NSETS+irfH1lf) = -1.0
         rf(k*NSETS+irfG1lf) = -1.0
         rf(k*NSETS+irfG2lf) = -1.0
#ifdef RT_UV
         rf(k*NSETS+irfUVlf) = -1.0
#endif

      endif

      return
      end
C
C
C
#ifdef RT_VARIABLE_RFIELD
      subroutine frtTransferUnPackHelper(k,y,w0,w,rf)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension y(*), w0(*), w(*), rf(*)
      include 'frt_inline.inc'
C
C  Works only for the NSETS=3 or 4 case
C
#ifdef RT_TRANSFER_FLUX_CONSERVING
      ta1 = y(ieqXH1)*rf(k*NSETS+irfH1lv) + rf(k*NSETS+irfH1lf)
      ta2 = y(ieqXG1)*rf(k*NSETS+irfG1lv) + rf(k*NSETS+irfG1lf)
      ta3 = y(ieqXG2)*rf(k*NSETS+irfG2lv) + rf(k*NSETS+irfG2lf)
#ifdef RT_UV
      ta4 = (y(ieqXH1)+2*y(ieqXMH))*rf(k*NSETS+irfUVlv) + 
     .     rf(k*NSETS+irfUVlf)
#endif
#else
      ta1 = rf(k*NSETS+irfH1lf)
      ta2 = rf(k*NSETS+irfG1lf)
      ta3 = rf(k*NSETS+irfG2lf)
#ifdef RT_UV
      ta4 = rf(k*NSETS+irfUVlf)
#endif
#endif
            
      if(ta2 .gt. 0.0) then
         if(ta3 .ge. 0.0) ta3 = csG1n(lrG2tr)*ta2 + ta3
      endif

      if(ta1 .gt. 0.0) then
         if(ta2 .ge. 0.0) ta2 = csH1n(lrG1tr)*ta1 + ta2
         if(ta3 .ge. 0.0) ta3 = csH1n(lrG2tr)*ta1 + ta3
      endif

      if(ta1 .lt. 0.0) then
         w(k*NSETS+1) = 0.0
      else if(k .eq. 0) then
         w(k*NSETS+1) = w0(k+1)*QFloc(ta1)
      else
         w(k*NSETS+1) = w0(k+1)*QFuni(ta1)
      endif

      if(ta2 .lt. 0.0) then
         w(k*NSETS+2) = 0.0
      else if(k .eq. 0) then
         w(k*NSETS+2) = w0(k+1)*QFloc(ta2)
      else
         w(k*NSETS+2) = w0(k+1)*QFuni(ta2)
      endif

      if(ta3 .lt. 0.0) then
         w(k*NSETS+3) = 0.0
      else if(k .eq. 0) then
         w(k*NSETS+3) = w0(k+1)*QFloc(ta3)
      else
         w(k*NSETS+3) = w0(k+1)*QFuni(ta3)
      endif

#ifdef RT_UV
      if(ta4 .lt. 0.0) then
         w(k*NSETS+4) = 0.0
      else if(k .eq. 0) then
         w(k*NSETS+4) = w0(k+1)*QFloc(ta4)
      else
         w(k*NSETS+4) = w0(k+1)*QFuni(ta4)
      endif
#endif

      return
      end
#endif /* RT_TRANSFER_FLUX_CONSERVING */
C
C  Absoprtion coefficient for 1 cell at frequency bin lr.
C
#if (RT_CFI == 1)
      subroutine frtTransferComputeCellAbs(L,Zsol,denB,denH1,denG1,
     .     denG2,denMH,dx,abc,abc1)
#else
      subroutine frtTransferComputeCellAbs(L,Zsol,denB,denH1,denG1,
     .     denG2,denMH,dx,abc)
#endif
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
C
      lr = LRSETS(L)
C
C  Adopt an approximation that H2 absorbs ionizing radiation as two HI
C
      denHn = denH1 + 2*denMH
C
#ifdef RT_UV
C
C  This is differnt from HI,HeI,HeII cases - we only include dust
C  absoprtion in 1 frequency in the RT solver, and add it later
C  approximately, because the RT solver seem to be unable to handle
C  it directly (loss of precision?).
C
      if(lr .eq. lrUVtr) then
         denHd = denH1 + 2*denMH/ClumpH2
         abc = uColumn*frtDustToGas(Zsol)*csUV(lr)*denHd
         abc = (1.0-exp(-abc*2*dx))/(2*dx)
      else
         abc = 0.0
      endif
#else
      abc = 0.0
#endif

#if (RT_CFI == 0)
      abc = abc + uColumn*(csH1(lr)*denHn+csG1(lr)*denG1+
     .     csG2(lr)*denG2)
#endif

#if (RT_CFI == 1)
      abc1 = abc + uColumn*(csH1(lr)*denHn+csG1(lr)*denG1+
     .     csG2(lr)*denG2)
#endif

#if (RT_CFI==1 || RT_CFI==2)
      abc = abc + uColumn*(
     .     csH1(lr)*denH1*CIonH1(denB,denHn/denB,denG1/denB,denG2/denB)
     .     csG1(lr)*denG1*CIonG1(denB,denHn/denB,denG1/denB,denG2/denB)
     .     csG2(lr)*denG2*CIonG2(denB,denHn/denB,denG1/denB,denG2/denB)
#endif

      return
      end
#endif /* RT_TRANSFER */
#endif /* RADIATIVE_TRANSFER */

      subroutine frtDummyTransfer(i)
      i = 0
      end
