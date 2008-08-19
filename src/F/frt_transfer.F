#include "defs.h"
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
#ifdef RT_TRANSFER
C
C
C
      subroutine frtInitRunTransfer(nfreq)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_transfer.inc'
C
C  Set sampling points here (except the last point)
C
      if(NSETS .ne. 3) then
         write(0,*) 'frtInitRunTransfer: fatal error, ',
     .        'only NSETS=3 is supported.'
         stop
      endif

      if(nfreq .ne. 2*NSETS) then
         write(0,*) 'frtInitRunTransfer: fatal error, ',
     .        'only nfreq=2*NSETS is supported.'
         stop
      endif

      LRSETS(1) = lrH1tr
      LRSETS(2) = lrG1tr
      LRSETS(3) = lrG2tr

      return
      end
C
C
C
      subroutine frtUpdateTablesTransfer(rfAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_transfer.inc'
      dimension rfAvg(*)
C
      dimension w(NSETS)

      do L=1,NSETS

         lr = LRSETS(L)

         aqbar = facLX*guveff(lr)
         ajbar = ang(lr)

         q = ajbar - aqbar*rfAvg(L)
         if(q .lt. 0.0) q = 0.0
         p = ajbar*rfAvg(NSETS+L)
         if(p .gt. 1.0e-35) then
            w(L) = q/p
         else
            w(L) = 0.0
         endif

      enddo

      tauH1 = QIhat(w(1))
      do lr=1,NRAD
         tauhat_xi(lr) = tauH1*csH1n(lr)
      enddo

      tauG1 = QIhat(w(2)) - tauH1*csH1n(lrG1tr)
      do lr=lrG1tr,NRAD
         tauhat_xi(lr) = tauhat_xi(lr) + tauG1*csG1n(lr)
      enddo

      tauG2 = QIhat(w(3)) - tauG1*csG1n(lrG2tr) - 
     .     tauH1*csH1n(lrG2tr)
      do lr=lrG2tr,NRAD
         tauhat_xi(lr) = tauhat_xi(lr) + tauG2*csG2n(lr)
      enddo

      do lr=1,NRAD
         qfdhat_xi(lr) = QFhat(tauhat_xi(lr),1)
      enddo

      return
      end
C
C
C
      subroutine frtStepBeginTransfer
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_transfer.inc'
C
      common /RT_Transfer_Save1/ abcBeg(NSETS)
C
      call frtTransferGetGlobalAbs(abcBeg)

      return
      end
C
C  This is hard-coded for NSETS=3 case
C
      subroutine frtTransferComputeAverageAbs(axi)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_transfer.inc'
      dimension axi(*)
C
      common /RT_Transfer_Save1/ abcBeg(NSETS)
      dimension abcEnd(NSETS)
C
      call frtTransferGetGlobalAbs(abcEnd)

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

      return
      end
C
C  Using a C-callback
C
      subroutine frtTransferGetGlobalAbs(abc)
      include 'frt_base.inc'
      include 'frt_transfer.inc'
      dimension abc(*)
      dimension ifield(NSETS)
C
C  Average global fields
C
      do i=1,NSETS
         ifield(i) = NSETS + i
      enddo 

      call rtTransferGetGlobalAbs(NSETS,ifield,abc)

      return
      end
C
C  Do some preprocessing
C
      subroutine frtTransferPackRadiationField(par,y,rawRF0,rawRF1,rf)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_transfer.inc'
      dimension par(*), y(*), rawRF0(*), rawRF1(*), rf(*)

      rf(irfNorm) = rawRF0(1)
C
C  1. Local field
C
      call frtTransferPackHelper(0,par,y,rawRF0,rawRF1,rf)
C
C  2. Global field
C
      call frtTransferPackHelper(1,par,y,rawRF0,rawRF1,rf)

#ifdef RT_TRANSFER_FLUX_CONSERVING
      rf(irfH1hw) = 0.5*csH1tr*par(iparDelN)
      rf(irfG1hw) = 0.5*csG1tr*par(iparDelN)
      rf(irfG2hw) = 0.5*csG2tr*par(iparDelN)
#endif

      return
      end
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
      include 'frt_transfer.inc'
      dimension par(*), y(*), w0(*), w(*), rf(*)
      dimension tauFront(3), tauWidth(3)
C
C  Works only for the NSETS=3 case
C
      if(w0(k+1) .gt. 1.0e-35) then

         if(w(k*NSETS+1) .gt. 1.0e-35) then
            t1 = min(1.0,w(k*NSETS+1)/w0(k+1))
            if(k .eq. 0) then
               ta1 = QIsor(t1)
            else
               ta1 = QIhat(t1)
            endif
         else
            ta1 = -1.0
         endif

         if(w(k*NSETS+2) .gt. 1.0e-35) then
            t2 = min(1.0,w(k*NSETS+2)/w0(k+1))
            if(k .eq. 0) then
               ta2 = QIsor(t2)
            else
               ta2 = QIhat(t2)
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
               ta3 = QIsor(t3)
            else
               ta3 = QIhat(t3)
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
         
         rf(k*NSETS+irfH1lf) = ta1
         rf(k*NSETS+irfG1lf) = ta2
         rf(k*NSETS+irfG2lf) = ta3

      else

         rf(k*NSETS+irfH1lf) = -1.0
         rf(k*NSETS+irfG1lf) = -1.0
         rf(k*NSETS+irfG2lf) = -1.0

      endif

      return
      end
#endif
#endif

      subroutine frtDummyTransfer(i)
      i = 0
      end
