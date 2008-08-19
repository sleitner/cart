#include "defs.h"
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
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
      dimension pH1(istabDim), sgpH1(istabDim)
      dimension pG1(istabDim), sgpG1(istabDim)
      dimension pG2(istabDim), sgpG2(istabDim)
C
C  Integration normalization
C
      upTabSE(istabPiH1H1w1) = facPi
      upTabSE(istabPiH1G1w1) = facPi
      upTabSE(istabPiH1G2w1) = facPi
      upTabSE(istabPiG1H1w1) = facPi
      upTabSE(istabPiG1G1w1) = facPi
      upTabSE(istabPiG1G2w1) = facPi

      upTabSE(istabPiH1H1w2) = facPi
      upTabSE(istabPiH1G1w2) = facPi
      upTabSE(istabPiH1G2w2) = facPi
      upTabSE(istabPiG1H1w2) = facPi
      upTabSE(istabPiG1G1w2) = facPi
      upTabSE(istabPiG1G2w2) = facPi

      upTabSE(istabPhH1w1) = facPh
      upTabSE(istabPhG1w1) = facPh
      upTabSE(istabPhH1w2) = facPh
      upTabSE(istabPhG1w2) = facPh
C
C  Main tables
C
#if defined(RT_TABLES) && !defined(RT_DUST)

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lo1,lo2,lo3,tauH1l,tauG1l,tauG2l,itind,itind2,j,
C$OMP+ pH1,pG1,pG2,sgpH1,sgpG1,sgpG2)
C$OMP+SHARED(tauH1m,acTAU,csH1tr,csG1tr,csG2tr,tauG1m,tauG2m,
C$OMP+ itloc2,pTabSE,sgpTabSE,pOT,sgpOT,wpTabSE,upTabSE,facLX)
      do lo1=1,NOPT1

         do j=1,istabDim
            pH1(j) = 0.0
#ifdef RT_TRANSFER
            sgpH1(j) = 0.0
#endif
         enddo

         if(lo1 .lt. NOPT1) then

            tauH1l = tauH1m*acTAU(lo1)
            call frtIntegrateRates_H1(istabDim,wpTabSE,
     .           istabH1,istabDim,pH1,
     .           0.0,tauH1l)
#ifdef RT_TRANSFER
            call frtLocalIntegrateRates_H1(istabDim,wpTabSE,
     .           istabH1,istabDim,sgpH1,
     .           0.0,tauH1l)
#endif

         else

            tauH1l = -1.0

         endif

         do lo2=1,NOPT1

            do j=1,istabDim
               pG1(j) = pH1(j)
#ifdef RT_TRANSFER
               sgpG1(j) = sgpH1(j)
#endif
            enddo

            if(lo2 .lt. NOPT1) then

               tauG1l = tauG1m*acTAU(lo2)
               call frtIntegrateRates_G1(istabDim,wpTabSE,
     .              istabG1,istabDim,pG1,
     .              0.0,tauH1l,tauG1l)
#ifdef RT_TRANSFER
               call frtLocalIntegrateRates_G1(istabDim,wpTabSE,
     .              istabG1,istabDim,sgpG1,
     .              0.0,tauH1l,tauG1l)
#endif

            else

               tauG1l = -1.0

            endif
jjj
            itind2 = itloc2(lo2,lo1)

            do lo3=1,NOPT1
C
C  No heating and ionizing HeII
C
               do j=1,istabDim
                  pG2(j) = pG1(j)
#ifdef RT_TRANSFER
                  sgpG2(j) = sgpG1(j)
#endif
               enddo

               if(lo3 .lt. NOPT1) then

                  tauG2l = tauG2m*acTAU(lo3)
                  call frtIntegrateRates_G2(istabDim,wpTabSE,
     .                 istabG2,istabDim,pG2,
     .                 0.0,tauH1l,tauG1l,tauG2l)
#ifdef RT_TRANSFER
                  call frtLocalIntegrateRates_G2(istabDim,wpTabSE,
     .                 istabG2,istabDim,sgpG2,
     .                 0.0,tauH1l,tauG1l,tauG2l)
#endif

               endif

               itind = lo3 + itind2
               do j=1,istabDim
#ifdef RT_INTERPOLLOG
                  if(pG2(j) .gt. 0.0) then
                     pTabSE(j,itind) = log(1.0e-35+upTabSE(j)*pG2(j))
                  else
                     pTabSE(j,itind) = -1000.0
                  endif                     
#ifdef RT_TRANSFER
                  if(sgpG2(j) .gt. 0.0) then
                     sgpTabSE(j,itind) = log(1.0e-35+
     .                    facLX*upTabSE(j)*sgpG2(j))
                  else
                     sgpTabSE(j,itind) = -1000.0
                  endif
#endif
#else  // RT_INTERPOLLOG
                  pTabSE(j,itind) = max(0.0,upTabSE(j)*pG2(j))
#ifdef RT_TRANSFER
                  sgpTabSE(j,itind) = 
     .                 max(0.0,facLX*upTabSE(j)*sgpG2(j))
#endif
#endif // RT_INTERPOLLOG
               enddo

            enddo
         enddo
      enddo

#endif  // defined(RT_TABLES) && !defined(RT_DUST)

      return
      end
C
C
C     
      subroutine frtGetPhotoRatesSE(par,rf,itab,y,pRate,sgpRate)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_TRANSFER
      include 'frt_transfer.inc'
#endif
      include 'frt_xrays.inc'
      dimension par(*), rf(*)
      dimension itab(*)
      dimension y(*), pRate(*), sgpRate(*)
      dimension pRate1(istabDim), sgpRate1(istabDim)
C
      include 'frt_xrays_fun.inc'
C
C  Global field
C
#ifdef RT_TRANSFER
      tauH1l = y(ieqXH1)*rf(irfH1gv) + rf(irfH1gf)
      tauG1l = y(ieqXG1)*rf(irfG1gv) + rf(irfG1gf)
      tauG2l = y(ieqXG2)*rf(irfG2gv) + rf(irfG2gf)
#else
      tauH1l = 0.0
      tauG1l = 0.0
      tauG2l = 0.0
#endif

#if defined(RT_TABLES) && !defined(RT_DUST)

      call frtSampleTable(tauH1l,tauG1l,tauG2l,
     .     istabDim,pTabSE,istabDim,pRate1,iOk)

      if(iOk .eq. 1) then
         itab(1) = itab(1) + 1
      else

#endif  // defined(RT_TABLES) && !defined(RT_DUST)

         do j=1,istabDim
            pRate1(j) = 0.0
         enddo

#if defined(RT_NARROWTABLE) || !defined(RT_TABLES) || defined(RT_DUST)

         itab(2) = itab(2) + 1

         call frtIntegrateRates_H1(istabDim,wpTabSE,
     .        istabH1,istabDim,pRate1,
     .        par(iparDGas),tauH1l)

         call frtIntegrateRates_G1(istabDim,wpTabSE,
     .        istabG1,istabDim,pRate1,
     .        par(iparDGas),tauH1l,tauG1l)

         call frtIntegrateRates_G2(istabDim,wpTabSE,
     .        istabG2,istabDim,pRate1,
     .        par(iparDGas),tauH1l,tauG1l,tauG2l)

         do j=1,istabDim
            pRate1(j) = upTabSE(j)*max(0.0,pRate1(j))
         enddo

#endif  // defined(RT_NARROWTABLE) || !defined(RT_TABLES) || defined(RT_DUST)

#if defined(RT_TABLES) && !defined(RT_DUST)
      endif
#endif
C
C  Local field
C
#ifdef RT_TRANSFER

      tauH1l = y(ieqXH1)*rf(irfH1lv) + rf(irfH1lf)
      tauG1l = y(ieqXG1)*rf(irfG1lv) + rf(irfG1lf)
      tauG2l = y(ieqXG2)*rf(irfG2lv) + rf(irfG2lf)

#if defined(RT_TABLES) && !defined(RT_DUST)

      call frtSampleTable(tauH1l,tauG1l,tauG2l,
     .     istabDim,sgpTabSE,istabDim,sgpRate1,iOk)

      if(iOk .eq. 1) then
         itab(1) = itab(1) + 1
      else

#endif  // defined(RT_TABLES) && !defined(RT_DUST)

         do j=1,istabDim
            sgpRate1(j) = 0.0
         enddo

#if defined(RT_NARROWTABLE) || !defined(RT_TABLES) || defined(RT_DUST)

         itab(2) = itab(2) + 1

         call frtLocalIntegrateRates_H1(istabDim,wpTabSE,
     .        istabH1,istabDim,sgpRate1,
     .        par(iparDGas),tauH1l)

         call frtLocalIntegrateRates_G1(istabDim,wpTabSE,
     .        istabG1,istabDim,sgpRate1,
     .        par(iparDGas),tauH1l,tauG1l)

         call frtLocalIntegrateRates_G2(istabDim,wpTabSE,
     .        istabG2,istabDim,sgpRate1,
     .        par(iparDGas),tauH1l,tauG1l,tauG2l)

         do j=1,istabDim
            sgpRate1(j) = facLX*upTabSE(j)*max(0.0,sgpRate1(j))
         enddo

#endif  // defined(RT_NARROWTABLE) || !defined(RT_TABLES) || defined(RT_DUST)

#if defined(RT_TABLES) && !defined(RT_DUST)
      endif
#endif  // defined(RT_TABLES) && !defined(RT_DUST)

#endif  //  RT_TRANSFER

      x = max(0.0,min(1.0,1.0-(y(ieqXH1)+y(ieqXG2)+2*y(ieqXG1))/XFEMAX))
      fh1 = se1fun(x,SE1H1A,SE1H1B,SE1H1C)
      fh2 = se2fun(x,SE2H1A,SE2H1B,SE2H1C)
      fg1 = se1fun(x,SE1G1A,SE1G1B,SE1G1C)
      fg2 = se2fun(x,SE2G1A,SE2G1B,SE2G1C)

      do j=1,istabDim2
         pRate(j) = fh1*pRate1(2*j-1) + fh2*pRate1(2*j)
#ifdef RT_TRANSFER
         sgpRate(j) = fh1*sgpRate1(2*j-1) + fh2*sgpRate1(2*j)
#endif
      enddo

      return
      end
#endif
#endif
      subroutine frtDummyXrays(i)
      i = 0
      end
