#include "frt_config.h"
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
         write(0,*) 'RT::InitRunTransfer: fatal error, ',
     .        'only NSETS=4 is supported.'
         call frtAbort
      endif
#else
      if(NSETS .ne. 3) then
         write(0,*) 'RT::InitRunTransfer: fatal error, ',
     .        'only NSETS=3 is supported.'
         call frtAbort
      endif
#endif

      if(nfreq .ne. NSETS) then
         write(0,*) 'RT::InitRunTransfer: fatal error, ',
     .        'only nfreq=NSETS is supported.'
         call frtAbort
      endif

      if(ncomp .ne. 2) then
         write(0,*) 'RT::InitRunTransfer: fatal error, ',
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
      subroutine frtSetGlobalRadiationFieldsTransfer(rfAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension rfAvg(2*NSETS)
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

         WSETS(L) = max(0.0,q)

      enddo

      tauH1 = QIuni(WSETS(1))
      do lr=1,NRAD
         tauLoc(lr) = tauH1*csH1n(lr)
      enddo

      tauG1 = max(0.0,QIuni(WSETS(2))-tauH1*csH1n(lrG1tr))
      do lr=lrG1tr,NRAD
         tauLoc(lr) = tauLoc(lr) + tauG1*csG1n(lr)
      enddo

      tauG2 = max(0.0,QIuni(WSETS(3))-tauG1*csG1n(lrG2tr)- 
     .     tauH1*csH1n(lrG2tr))
      do lr=lrG2tr,NRAD
         tauLoc(lr) = tauLoc(lr) + tauG2*csG2n(lr)
      enddo

#ifdef RT_OUTPUT
      if(ioNode) then
         open(unit=IOUNIT, file='rt-run.log', access='append')
         write(IOUNIT,*) '  Photon balance:'
         do L=1,NSETS
            lr = LRSETS(L)
            aqbar = facLX*guvEff(lr)
            ajbar = angLoc(lr)
            write(IOUNIT,901) L, WSETS(L), rfAvg(L), 
     .           aqbar*rfAvg(L), ajbar, tauLoc(lr)
         enddo
         close(IOUNIT)
      endif
 901  format('  L=',I1,1X,F5.3,4(1X,E11.3))
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
      subroutine frtStepEndTransfer(abc)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension abc(*)
C
      common /RT_Transfer_Save/ abcAvg(NSETS)
C
      do i=1,NSETS
         abcAvg(i) = abc(i)
      enddo

      do i=1,NSETS
         if(abcAvg(i) .lt. 0.0) then
            write(0,*) 'RT::StepEndTransfer: fatal error, ',
     .           'negative global absoprtion: ', i, abcAvg(i)
            call frtAbort
         endif
      enddo

      return
      end
C
C  This is hard-coded for NSETS=3 or 4 case
C
      function frtTransferGlobalAC(L,abc)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
      dimension abc(*)

      frtTransferGlobalAC = abc(1)*(1.0-WSETS(L)) + WSETS(L)*abc(2)
*      frtTransferGlobalAC = abc(2)

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
      common /RT_Transfer_Save/ abcAvg(NSETS)

      axiH1 = max(0.0,abcAvg(1)/facLX)
      do lr=lrH1tr,lrmax
         axi(lr) = axi(lr) + axiH1*csH1n(lr)
      enddo

      axiG1 = max(0.0,abcAvg(2)/facLX-axiH1*csH1n(lrG1tr))
      do lr=lrG1tr,lrmax
         axi(lr) = axi(lr) + axiG1*csG1n(lr)
      enddo

      axiG2 = max(0.0,abcAvg(3)/facLX-axiG1*csG1n(lrG2tr)
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
      dimension tau(NSETS), x(NSETS)
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
      do L=1,NSETS
         rf(irfOgf+L) = 10*VLARGE
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

      x(1) = max(varLow(ivarXH1),var(ivarXH1))
      x(2) = max(varLow(ivarXG1),var(ivarXG1))
      x(3) = max(varLow(ivarXG2),var(ivarXG2))
#ifdef RT_UV
      x(4) = max(varLow(ivarXH1),var(ivarXH1)+2*var(ivarXMH))
#endif

      do L=1,NSETS
         if(rf(irfOlf+L) .gt. VLARGE) then
            rf(irfOlf+L) = 10*VLARGE
            rf(irfOlv+L) = 0.0
         else if(rf(irfOlf+j) .gt. 0.0) then
            rf(irfOlf+L) = rf(irfOlf+L) - wFlux*tau(L)
            rf(irfOlv+L) = wFlux*tau(L)/x(L)
         endif
#ifdef RT_EXTERNAL_BACKGROUND
         if(rf(irfOgf+L) .gt. VLARGE) then
            rf(irfOgf+L) = 10*VLARGE
            rf(irfOgv+L) = 0.0
         else if(rf(irfOgf+L) .gt. 0.0)
            rf(irfOgf+L) = rf(irfOgf+L) - wFlux*tau(L)
            rf(irfOgv+L) = wFlux*tau(L)/x(L)
         endif
#endif
         rf(irfOwf+L) = (1-wCell)*tau(L)
         rf(irfOwv+L) = wCell*tau(L)/x(L)
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
      dimension tau(NSETS)
      include 'frt_inline.inc'
C
C  Works only for the NSETS=3 or 4 case
C
      w0 = var(ivarRF0l+k)

      if(w0 .lt. VSMALL) then
         do L=1,NSETS
            rf(k*NSETS+irfOlf+L) = 10*VLARGE
         enddo
         return
      endif

      do L=1,NSETS
         if(w(k*NSETS+L) .gt. VSMALL) then
            q = w(k*NSETS+L)/w0
            if(k .eq. 0) then
               tau(L) = QIloc(q)
            else
               tau(L) = QIuni(q)
            endif
         else
            tau(L) = 10*VLARGE
         endif
      enddo
C
C  Now invert matrix.
C  We don't need to check whether we subtract from 10*VLARGE since
C  no matter how much we subtract, the value will still be above VLARGE
C
#if defined(RT_UV) && !defined(RT_UV_OLDSTYLE_3x1)
      if(tau(4) .lt. VLARGE) then
         do L=1,3
            tau(L) = tau(L) - tau(4)*csUVn(LRSETS(L))
         enddo
      endif
#endif /* RT_UV && !RT_UV_OLDSTYLE_3x1 */

      if(tau(1) .lt. VLARGE) then
         tau(2) = tau(2) - tau(1)*csH1n(lrG1tr)
         tau(3) = tau(3) - tau(1)*csH1n(lrG2tr)
      endif

      if(tau(2) .lt. VLARGE) then
         tau(3) = tau(3) - tau(2)*csG1n(lrG2tr)
      endif

      do L=1,NSETS
         rf(k*NSETS+irfOlf+L) = tau(L)
      enddo

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
      dimension tau(NSETS), x(NSETS)
      include 'frt_inline.inc'
C
C  Works only for the NSETS=3 or 4 case
C
      x(1) = max(varLow(ivarXH1),var(ivarXH1))
      x(2) = max(varLow(ivarXG1),var(ivarXG1))
      x(3) = max(varLow(ivarXG2),var(ivarXG2))
#ifdef RT_UV
      x(4) = max(varLow(ivarXH1),var(ivarXH1)+2*var(ivarXMH))
#endif

      do L=1,NSETS
         tau(L) = rf(k*NSETS+irfOlf+L) + x(L)*rf(k*NSETS+irfOlv+L)
      enddo
C
C  Now apply matrix.
C
      if(tau(2) .lt. VLARGE) then
         tau(3) = tau(3) + csG1n(lrG2tr)*tau(2)
      endif

      if(tau(1) .lt. VLARGE) then
         tau(2) = tau(2) + csH1n(lrG1tr)*tau(1)
         tau(3) = tau(3) + csH1n(lrG2tr)*tau(1)
      endif

#if defined(RT_UV) && !defined(RT_UV_OLDSTYLE_3x1)
      if(tau(4) .lt. VLARGE) then
         do L=1,3
            tau(L) = tau(L) + tau(4)*csUVn(LRSETS(L))
         enddo
      endif
#endif /* RT_UV && !RT_UV_OLDSTYLE_3x1 */

      w0 = var(ivarRF0l+k)

      do L=1,NSETS
         if(tau(L) .lt. VLARGE) then
            if(k .eq. 0) then
               w(k*NSETS+L) = w0*QFloc(tau(L))
            else
               w(k*NSETS+L) = w0*QFuni(tau(L))
            endif
         else
            w(k*NSETS+L) = 0.0
         endif
      enddo

      return
      end
#endif /* RT_VARIABLE_RF */
C
C  Scaling factors (= units) for the radiation field L
C  (only works with new UV)
C     
      subroutine frtGetRFUnits(L,uLoc,uGlob)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
 
      if(L.ge.1 .and. L.le.NSETS) then

         lr = LRSETS(L)

         uLoc = qExp3*facLX*guvEff(lr)

#ifdef RT_EXTERNAL_BACKGROUND
         uGlob = qExp3*angEff(lr)
#else
         uGlob = 0.0
#endif

      else

         uLoc = 0.0
         uGlob = 0.0

      endif

      return
      end
C
C  This works with old and new UV model
C
      subroutine frtGetRFUnits2(var,rawrf,L,uLoc,uGlob)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension var(*), rawrf(*)
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
#ifdef RT_TRANSFER_FLUX_CONSERVING
      dimension rf(irfDim)
#else
      dimension rf(irfMax)
#endif
      include 'frt_inline.inc'

      if(L.ge.1 .and. L.le.NSETS) then

#ifdef RT_TRANSFER_FLUX_CONSERVING
         var(ivarDelN) = uColumn*var(ivarRho)*var(ivarDelX)
#endif
         call frtTransferPackRadiationField(var,rawrf,rf)

         call frtGetRFHelper(LRSETS(L),rf,uLoc,uGlob)
         
         if(rawrf(L+0*NSETS) .gt. 1.0e-30) then
            uLoc = uLoc*(rf(irfNorm)/rawrf(L+0*NSETS))
         else
            uLoc = 0.0
         endif
#ifdef RT_EXTERNAL_BACKGROUND
         if(rawrf(L+1*NSETS) .gt. 1.0e-30) then
            uGlob = uGlob*(1.0/rawrf(L+1*NSETS))
         else
            uGlob = 0.0
         endif
#endif /* RT_EXTERNAL_BACKGROUND */

      else

         uLoc = 0.0
         uGlob = 0.0

      endif

      return
      end
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
#ifdef RT_UV_OLDSTYLE_3x1
C
C  This is differnt from HI,HeI,HeII cases - we only include dust
C  absoprtion in 1 frequency in the RT solver, and add it later
C  approximately, because the RT solver seem to be unable to handle
C  it directly (loss of precision?).
C
      if(lr .eq. lrUVtr) then
#endif /* RT_UV_OLDSTYLE_3x1 */
         denHd = denH1 + 2*denMH/ClumpH2
#ifdef RT_CUSTOM_DUST_TO_GAS
         abc(1) = uColumn*frtDustToGas(var)*csUV(lr)*denHd
#else
         abc(1) = uColumn*Zsol*csUV(lr)*denHd
#endif
         abc(1) = (1.0-exp(-abc(1)*dx))/dx
#ifdef RT_UV_OLDSTYLE_3x1
      else
         abc(1) = 0.0
      endif
#endif /* RT_UV_OLDSTYLE_3x1 */
#else  /* RT_UV */
      abc(1) = 0.0
#endif /* RT_UV */

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
