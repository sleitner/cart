#include "../config.h"
#ifdef RADIATIVE_TRANSFER
#ifdef RT_TRANSFER
C
C
C
      subroutine frtInitRunTransfer(nfreq,ncomp)
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

      if(nfreq .ne. NSETS) then
         write(0,*) 'RT_Transfer: fatal error, ',
     .        'only nfreq=NSETS is supported.'
         call frtAbort
      endif

      if(ncomp .ne. 2) then
         write(0,*) 'RT_Transfer: fatal error, ',
     .        'only ncomp=2 is supported.'
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
      subroutine frtSetGlobalRadiationFieldsTransfer(rfAvg,otAvg)
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
         ajbar = angLoc(lr)

         if(ajbar .gt. 0.0) then
            q = 1.0 - aqbar*rfAvg(L)/ajbar
         else
            q = 0.0
         endif

         w(L) = max(0.0,q)

      enddo

      tauH1 = QIuni(w(1))
      do lr=1,NRAD
         tauLoc(lr) = tauH1*csH1n(lr)
      enddo

      tauG1 = max(0.0,QIuni(w(2))-tauH1*csH1n(lrG1tr))
      do lr=lrG1tr,NRAD
         tauLoc(lr) = tauLoc(lr) + tauG1*csG1n(lr)
      enddo

      tauG2 = max(0.0,QIuni(w(3))-tauG1*csG1n(lrG2tr)- 
     .     tauH1*csH1n(lrG2tr))
      do lr=lrG2tr,NRAD
         tauLoc(lr) = tauLoc(lr) + tauG2*csG2n(lr)
      enddo

#ifdef RT_OUTPUT
      if(ioNode) then
         open(unit=IOUNIT, file='rt-run.log', access='append')
         write(IOUNIT,*) '  Photon balance:'
         write(IOUNIT,*) '    nu_H1: ', w(1), tauH1, tauLoc(lrH1tr)
         write(IOUNIT,*) '    nu_G1: ', w(2), tauG1, tauLoc(lrG1tr)
         write(IOUNIT,*) '    nu_G2: ', w(3), tauG2, tauLoc(lrG2tr)
         close(IOUNIT)
      endif
#endif /* RT_OUTPUT */

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
      subroutine frtTransferPackRadiationField(var,rawrf,rf)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension var(*), rawrf(*), rf(*)
#if defined(RT_TRANSFER_FLUX_CONSERVING) || defined(RT_VARIABLE_RF)
      parameter(wFlux = 0.0, wCell = 0.0)
      dimension tau(NSETS)
#endif

      rf(irfNorm) = var(ivarRF0l)
C
C  1. Local field
C
      call frtTransferPackHelper(0,var,rawrf,rf)
C
C  2. Global field
C
#ifdef RT_EXTERNAL_BACKGROUND
      call frtTransferPackHelper(1,var,rawrf,rf)
#else
      do j=1,NSETS
         rf(irfOgf+j) = -1.0
      enddo
#endif

#if defined(RT_TRANSFER_FLUX_CONSERVING) || defined(RT_VARIABLE_RF)
      tau(1) = csH1tr*var(ivarDelN)
      tau(2) = csG1tr*var(ivarDelN)
      tau(3) = csG2tr*var(ivarDelN)
#ifdef RT_UV
#ifdef RT_CUSTOM_DUST_TO_GAS
      tau(4) = frtDustToGas(var)*csUVtr*var(ivarDelN)
#else
      tau(4) = var(ivarZsol)*csUVtr*var(ivarDelN)
#endif
#endif

      do j=1,NSETS
         q = 1.0/max(varLow(1+j),var(1+j))
         if(rf(irfOlf+j) .lt. 0.0) then
            rf(irfOlf+j) = -1.0
            rf(irfOlv+j) = 0.0
         else
            rf(irfOlf+j) = max(0.0,rf(irfOlf+j)-wFlux*tau(j))
            rf(irfOlv+j) = wFlux*tau(j)*q
         endif
#ifdef RT_EXTERNAL_BACKGROUND
         if(rf(irfOgf+j) .lt. 0.0) then
            rf(irfOgf+j) = -1.0
            rf(irfOgv+j) = 0.0
         else
            rf(irfOgf+j) = max(0.0,rf(irfOgf+j)-wFlux*tau(j))
            rf(irfOgv+j) = wFlux*tau(j)*q
         endif
#endif
         rf(irfOwf+j) = (1-wCell)*tau(j)
         rf(irfOwv+j) = wCell*tau(j)*q
      enddo
#endif

      return
      end
C
C Do some postprocessing
C
#ifdef RT_VARIABLE_RF
      subroutine frtTransferUnPackRadiationField(var,rawrf,rf)
      dimension var(*), rawrf(*), rf(*)
C
C  1. Local field
C
      call frtTransferUnPackHelper(0,var,rawrf,rf)
C
C  2. Global field
C
      call frtTransferUnPackHelper(1,var,rawrf,rf)

      return
      end
#endif /* RT_VARIABLE_RF */
C
C  *******************
C
C    Helpers
C
C  *******************
C
      subroutine frtTransferPackHelper(k,var,w,rf)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension var(*), w(*), rf(*)
      include 'frt_inline.inc'
C
C  Works only for the NSETS=3 or 4 case
C
      w0 = var(ivarRF0l+k)
      if(w0 .gt. 1.0e-35) then

         if(w(k*NSETS+1) .gt. 1.0e-35) then
            t1 = min(1.0,w(k*NSETS+1)/w0)
            if(k .eq. 0) then
               ta1 = QIloc(t1)
            else
               ta1 = QIuni(t1)
            endif
         else
            ta1 = -1.0
         endif

         if(w(k*NSETS+2) .gt. 1.0e-35) then
            t2 = min(1.0,w(k*NSETS+2)/w0)
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
            t3 = min(1.0,w(k*NSETS+3)/w0)
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
            t4 = min(1.0,w(k*NSETS+4)/w0)
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
#ifdef RT_VARIABLE_RF
      subroutine frtTransferUnPackHelper(k,var,w,rf)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension var(*), w(*), rf(*)
      include 'frt_inline.inc'
C
C  Works only for the NSETS=3 or 4 case
C
      ta1 = var(ivarXH1)*rf(k*NSETS+irfH1lv) + rf(k*NSETS+irfH1lf)
      ta2 = var(ivarXG1)*rf(k*NSETS+irfG1lv) + rf(k*NSETS+irfG1lf)
      ta3 = var(ivarXG2)*rf(k*NSETS+irfG2lv) + rf(k*NSETS+irfG2lf)
#ifdef RT_UV
      ta4 = (var(ivarXH1)+2*var(ivarXMH))*rf(k*NSETS+irfUVlv) + 
     .     rf(k*NSETS+irfUVlf)
#endif
            
      if(ta2 .gt. 0.0) then
         if(ta3 .ge. 0.0) ta3 = csG1n(lrG2tr)*ta2 + ta3
      endif

      if(ta1 .gt. 0.0) then
         if(ta2 .ge. 0.0) ta2 = csH1n(lrG1tr)*ta1 + ta2
         if(ta3 .ge. 0.0) ta3 = csH1n(lrG2tr)*ta1 + ta3
      endif

      w0 = var(ivarRF0l+k)

      if(ta1 .lt. 0.0) then
         w(k*NSETS+1) = 0.0
      else if(k .eq. 0) then
         w(k*NSETS+1) = w0*QFloc(ta1)
      else
         w(k*NSETS+1) = w0*QFuni(ta1)
      endif

      if(ta2 .lt. 0.0) then
         w(k*NSETS+2) = 0.0
      else if(k .eq. 0) then
         w(k*NSETS+2) = w0*QFloc(ta2)
      else
         w(k*NSETS+2) = w0*QFuni(ta2)
      endif

      if(ta3 .lt. 0.0) then
         w(k*NSETS+3) = 0.0
      else if(k .eq. 0) then
         w(k*NSETS+3) = w0*QFloc(ta3)
      else
         w(k*NSETS+3) = w0*QFuni(ta3)
      endif

#ifdef RT_UV
      if(ta4 .lt. 0.0) then
         w(k*NSETS+4) = 0.0
      else if(k .eq. 0) then
         w(k*NSETS+4) = w0*QFloc(ta4)
      else
         w(k*NSETS+4) = w0*QFuni(ta4)
      endif
#endif

      return
      end
#endif /* RT_VARIABLE_RF */
C
C  Absorption coefficient for 1 cell at frequency bin lr.
C
#ifdef RT_ABSORPTION_CALLBACK_FULL
      subroutine frtTransferComputeCellAbs(L,denB,denH1,denG1,
     .     denG2,denMH,dx,abc,var)
#else
      subroutine frtTransferComputeCellAbs(L,denB,denH1,denG1,
     .     denG2,denMH,dx,abc,Zsol)
#endif
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension abc(*)
#ifdef RT_ABSORPTION_CALLBACK_FULL
      dimension var(*)
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
#ifdef RT_CUSTOM_DUST_TO_GAS
         abc(1) = uColumn*frtDustToGas(var)*csUV(lr)*denHd
#else
         abc(1) = uColumn*Zsol*csUV(lr)*denHd
#endif
         abc(1) = (1.0-exp(-abc(1)*2*dx))/(2*dx)
      else
         abc(1) = 0.0
      endif
#else
      abc(1) = 0.0
#endif

#if (RT_CFI == 0)
      abc(1) = abc(1) + uColumn*(csH1(lr)*denHn+csG1(lr)*denG1+
     .     csG2(lr)*denG2)
#endif

#if (RT_CFI == 1)
      abc(2) = abc(1) + uColumn*(csH1(lr)*denHn+csG1(lr)*denG1+
     .     csG2(lr)*denG2)
#endif

#if (RT_CFI==1 || RT_CFI==2)
      abc(1) = abc(1) + uColumn*(csH1(lr)*denH1*CIonH1(var)+
     .     csG1(lr)*denG1*CIonG1(var)+csG2(lr)*denG2*CIonG2(var))
#endif

      return
      end
#endif /* RT_TRANSFER */
#endif /* RADIATIVE_TRANSFER */

      subroutine frtDummyTransfer(i)
      i = 0
      end
