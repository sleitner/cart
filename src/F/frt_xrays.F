#include "../config.h"
#ifdef RADIATIVE_TRANSFER
#ifdef RT_XRAYS
C
C
C
      subroutine frtInitTablesSE
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_xrays.inc'

      do lr=1,NRAD

         wpTabSE(istabPiH1H1w1,lr) = csH1(lr)*se1H1i(lr)*
     .        max(0.0,Txi(lr)-TH1-TH1)/TH1
         wpTabSE(istabPiH1G1w1,lr) = csG1(lr)*se1H1i(lr)*
     .        max(0.0,Txi(lr)-TH1-TG1)/TH1
         wpTabSE(istabPiH1G2w1,lr) = csG2(lr)*se1H1i(lr)*
     .        max(0.0,Txi(lr)-TH1-TG2)/TH1
         wpTabSE(istabPiG1H1w1,lr) = csH1(lr)*se1G1i(lr)*
     .        max(0.0,Txi(lr)-TG1-TH1)/TG1
         wpTabSE(istabPiG1G1w1,lr) = csG1(lr)*se1G1i(lr)*
     .        max(0.0,Txi(lr)-TG1-TG1)/TG1
         wpTabSE(istabPiG1G2w1,lr) = csG2(lr)*se1G1i(lr)*
     .        max(0.0,Txi(lr)-TG1-TG2)/TG1
         
         wpTabSE(istabPiH1H1w2,lr) = csH1(lr)*se2H1i(lr)*
     .        max(0.0,Txi(lr)-TH1-TH1)/TH1
         wpTabSE(istabPiH1G1w2,lr) = csG1(lr)*se2H1i(lr)*
     .        max(0.0,Txi(lr)-TH1-TG1)/TH1
         wpTabSE(istabPiH1G2w2,lr) = csG2(lr)*se2H1i(lr)*
     .        max(0.0,Txi(lr)-TH1-TG2)/TH1
         wpTabSE(istabPiG1H1w2,lr) = csH1(lr)*se2G1i(lr)*
     .        max(0.0,Txi(lr)-TG1-TH1)/TG1
         wpTabSE(istabPiG1G1w2,lr) = csG1(lr)*se2G1i(lr)*
     .        max(0.0,Txi(lr)-TG1-TG1)/TG1
         wpTabSE(istabPiG1G2w2,lr) = csG2(lr)*se2G1i(lr)*
     .        max(0.0,Txi(lr)-TG1-TG2)/TG1

         wpTabSE(istabPhH1w1,lr) = csH1(lr)*se1H1h(lr)*(Txi(lr)-TH1)
         wpTabSE(istabPhG1w1,lr) = csG1(lr)*se1G1h(lr)*(Txi(lr)-TG1)
         wpTabSE(istabPhH1w2,lr) = csH1(lr)*se2H1h(lr)*(Txi(lr)-TH1)
         wpTabSE(istabPhG1w2,lr) = csG1(lr)*se2G1h(lr)*(Txi(lr)-TG1)

      enddo

      return
      end
C
C
C
      subroutine frtFillRadiationTablesSE
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_xrays.inc'
C
      dimension pH1(istabDim), qH1(istabDim)
      dimension pG1(istabDim), qG1(istabDim)
      dimension pG2(istabDim), qG2(istabDim)
C
C  Main tables
C
#ifndef RT_NO_TABLE

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lo1,lo2,lo3,tauH1l,tauG1l,tauG2l,itind,itind2,j,
C$OMP+ pH1,pG1,pG2,qH1,qG1,qG2)
C$OMP+SHARED(tauH1m,acTAU,csH1tr,csG1tr,csG2tr,tauG1m,tauG2m,
C$OMP+ itloc2,pTabSE,qTabSE,wpTabSE,facPhot,facLX)
      do lo1=1,NOPT1

         do j=1,istabDim
#ifdef RT_EXTERNAL_BACKGROUND
            pH1(j) = 0.0
#endif
#ifdef RT_TRANSFER
            qH1(j) = 0.0
#endif
         enddo

         if(lo1 .lt. NOPT1) then

            tauH1l = tauH1m*acTAU(lo1)
#ifdef RT_EXTERNAL_BACKGROUND
            call frtIntegrateRates_H1(0,istabDim,wpTabSE,
     .           istabH1,istabDim,pH1,0.0,tauH1l)
#endif
#ifdef RT_TRANSFER
            call frtIntegrateRates_H1(1,istabDim,wpTabSE,
     .           istabH1,istabDim,qH1,0.0,tauH1l)
#endif

         else

            tauH1l = -1.0

         endif

         do lo2=1,NOPT1

            do j=1,istabDim
#ifdef RT_EXTERNAL_BACKGROUND
               pG1(j) = pH1(j)
#endif
#ifdef RT_TRANSFER
               qG1(j) = qH1(j)
#endif
            enddo

            if(lo2 .lt. NOPT1) then

               tauG1l = tauG1m*acTAU(lo2)
#ifdef RT_EXTERNAL_BACKGROUND
               call frtIntegrateRates_G1(0,istabDim,wpTabSE,
     .              istabG1,istabDim,pG1,0.0,tauH1l,tauG1l)
#endif
#ifdef RT_TRANSFER
               call frtIntegrateRates_G1(1,istabDim,wpTabSE,
     .              istabG1,istabDim,qG1,0.0,tauH1l,tauG1l)
#endif

            else

               tauG1l = -1.0

            endif

            itind2 = itloc2(lo2,lo1)

            do lo3=1,NOPT1
C
C  No heating and ionizing HeII
C
               do j=1,istabDim
#ifdef RT_EXTERNAL_BACKGROUND
                  pG2(j) = pG1(j)
#endif
#ifdef RT_TRANSFER
                  qG2(j) = qG1(j)
#endif
               enddo

               if(lo3 .lt. NOPT1) then

                  tauG2l = tauG2m*acTAU(lo3)
#ifdef RT_EXTERNAL_BACKGROUND
                  call frtIntegrateRates_G2(0,istabDim,wpTabSE,
     .                 istabG2,istabDim,pG2,0.0,tauH1l,tauG1l,tauG2l)
#endif
#ifdef RT_TRANSFER
                  call frtIntegrateRates_G2(1,istabDim,wpTabSE,
     .                 istabG2,istabDim,qG2,0.0,tauH1l,tauG1l,tauG2l)
#endif

               endif

               itind = lo3 + itind2
               do j=1,istabDim
#ifdef RT_INTERPOLLOG
#ifdef RT_EXTERNAL_BACKGROUND
                  if(pG2(j) .gt. 0.0) then
                     pTabSE(j,itind) = log(1.0e-35+facPhot*pG2(j))
                  else
                     pTabSE(j,itind) = -1000.0
                  endif                  
#endif
#ifdef RT_TRANSFER
                  if(qG2(j) .gt. 0.0) then
                     qTabSE(j,itind) = log(1.0e-35+facLX*facPhot*qG2(j))
                  else
                     qTabSE(j,itind) = -1000.0
                  endif
#endif
#else  // RT_INTERPOLLOG
#ifdef RT_EXTERNAL_BACKGROUND
                  pTabSE(j,itind) = max(0.0,facPhot*pG2(j))
#endif
#ifdef RT_TRANSFER
                  qTabSE(j,itind) = max(0.0,facLX*facPhot*qG2(j))
#endif
#endif // RT_INTERPOLLOG
               enddo

            enddo
         enddo
      enddo

#endif /* RT_NO_TABLE */

      return
      end
C
C
C     
      subroutine frtPhotoRatesSE(var,rf,itab,pRate,qRate)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_TRANSFER
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
#endif
      include 'frt_xrays.inc'
      dimension var(*), rf(*)
      dimension itab(*)
      dimension pRate(*), qRate(*)
      dimension pRate1(istabDim), qRate1(istabDim)
C
      include 'frt_xrays_fun.inc'
C
C  Global field
C
#ifdef RT_EXTERNAL_BACKGROUND
#ifdef RT_TRANSFER
      tauH1l = rf(irfH1gf)
      tauG1l = rf(irfG1gf)
      tauG2l = rf(irfG2gf)
#ifdef RT_UV
      tauUVl = rf(irfUVgf)
#else
      tauUVl = 0.0
#endif
#else
      tauH1l = 0.0
      tauG1l = 0.0
      tauG2l = 0.0
      tauUVl = 0.0
#endif

#ifndef RT_NO_TABLE

      call frtSampleTable3(tauH1l,tauG1l,tauG2l,
     .     istabDim,pTabSE,istabDim,pRate1,iOk)

      if(iOk .eq. 1) then
         itab(1) = itab(1) + 1

#ifdef RT_UV
         do j=1,istabDim
            pRate1(j) = pRate1(j)*exp(-tauUVl*uvFac(j))
         enddo
#endif
      else

#endif /* RT_NO_TABLE */

         do j=1,istabDim
            pRate1(j) = 0.0
         enddo

#if defined(RT_NARROW_TABLE) || defined(RT_NO_TABLE)

         itab(2) = itab(2) + 1

         call frtIntegrateRates_H1(0,istabDim,wpTabSE,
     .        istabH1,istabDim,pRate1,tauUVl,tauH1l)

         call frtIntegrateRates_G1(0,istabDim,wpTabSE,
     .        istabG1,istabDim,pRate1,tauUVl,tauH1l,tauG1l)

         call frtIntegrateRates_G2(0,istabDim,wpTabSE,
     .        istabG2,istabDim,pRate1,tauUVl,tauH1l,tauG1l,tauG2l)

         do j=1,istabDim
            pRate1(j) = facPhot*max(0.0,pRate1(j))
         enddo

#endif /* RT_NARROW_TABLE || RT_NO_TABLE */

#ifndef RT_NO_TABLE
      endif
#endif
#endif // RT_EXTERNAL_BACKGROUND
C
C  Local field
C
#ifdef RT_TRANSFER
      tauH1l = rf(irfH1lf)
      tauG1l = rf(irfG1lf)
      tauG2l = rf(irfG2lf)
#ifdef RT_UV
      tauUVl = rf(irfUVlf)
#else
      tauUVl = 0.0
#endif

#ifndef RT_NO_TABLE

      call frtSampleTable3(tauH1l,tauG1l,tauG2l,
     .     istabDim,qTabSE,istabDim,qRate1,iOk)

      if(iOk .eq. 1) then
         itab(1) = itab(1) + 1

#ifdef RT_UV
         do j=1,istabDim
            qRate1(j) = qRate1(j)*exp(-tauUVl*uvFac(j))
         enddo
#endif
      else

#endif /* RT_NO_TABLE */

         do j=1,istabDim
            qRate1(j) = 0.0
         enddo

#if defined(RT_NARROW_TABLE) || defined(RT_NO_TABLE)

         itab(2) = itab(2) + 1

         call frtIntegrateRates_H1(1,istabDim,wpTabSE,
     .        istabH1,istabDim,qRate1,tauUVl,tauH1l)

         call frtIntegrateRates_G1(1,istabDim,wpTabSE,
     .        istabG1,istabDim,qRate1,tauUVl,tauH1l,tauG1l)

         call frtIntegrateRates_G2(1,istabDim,wpTabSE,
     .        istabG2,istabDim,qRate1,tauUVl,tauH1l,tauG1l,tauG2l)

         do j=1,istabDim
            qRate1(j) = facLX*facPhot*max(0.0,qRate1(j))
         enddo

#endif /* RT_NARROWTABLE || RT_NO_TABLE */

#ifndef RT_NO_TABLE
      endif
#endif /* RT_NO_TABLE */

#endif /* RT_TRANSFER */

      x = max(0.0,min(1.0,1.0-(var(ivarXH1)+var(ivarXG2)+
     .     2*var(ivarXG1))/XFEMAX))
      fh1 = se1fun(x,SE1H1A,SE1H1B,SE1H1C)
      fh2 = se2fun(x,SE2H1A,SE2H1B,SE2H1C)
      fg1 = se1fun(x,SE1G1A,SE1G1B,SE1G1C)
      fg2 = se2fun(x,SE2G1A,SE2G1B,SE2G1C)

      do j=1,istabDim2
#ifdef RT_EXTERNAL_BACKGROUND
         pRate(j) = fh1*pRate1(2*j-1) + fh2*pRate1(2*j)
#endif
#ifdef RT_TRANSFER
         qRate(j) = fh1*qRate1(2*j-1) + fh2*qRate1(2*j)
#endif
      enddo

      return
      end
#endif
#endif
      subroutine frtDummyXrays(i)
      i = 0
      end
