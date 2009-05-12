#include "defs.h"      
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
C
C
C
      subroutine frtInitTables
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
#ifdef RT_XRAYS
      include 'frt_xrays.inc'
#endif
C
      common/RT_Z_LYSS/ oscstr(48), oscstr_other
      common/RT_Z_SD93/ sd93lt(90), sd93c0(90), sd93c1(90)
C
C  The fractional cross-section at the last point
C
      parameter (CSNMIN = 1.0e-7)
C
#ifdef RT_XRAYS
      parameter (EM5 = 1.0e-5)
      include 'frt_xrays_fun.inc'
#endif
C
C  Dust cross section fitting formula
C
      dustfit(w,wi,ai,bi,pw1,pw2) = ai/((w/wi)**pw1+(wi/w)**pw2+bi)
C
C  SD93 re-interpolation
C
      dsd93 = sd93lt(2) - sd93lt(1)
      nsd93 = 90

      temZLT = 10.0**sd93lt(1)
C
C  Wolfire et al (ApJ 587 278, 2003) phiPAH parameter and
C  the effective cross-section (from W&D 2003, ApJ 528, 296)
C
      d_pPAH = 0.5
      d_sPAH = 3.4e-21/CSSTAR
C
C  Radiation tables
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,fnu,x,cs0,fn0,y0,y1,yw,ya,p,y,w)
C$OMP+SHARED(Txi,xig,angUni,angLoc,angAvg,tauLoc,csH1,csG1,csG2,csDust,
C$OMP+ csMH27,csMH28,csMH29,csMH30,csMH31,csMH32,csMHLW)
      do lr = 1,NRAD

         xig(lr) = XISTEP*(real(lr)-0.5) + log(TXI0/TH1)
         Txi(lr) = TH1*exp(xig(lr))
         angUni(lr) = 0.0
         angLoc(lr) = 0.0
         angAvg(lr) = 0.0
         tauLoc(lr) = 0.0
         fnu = Txi(lr)/T1eV

         if(Txi(lr) .lt. TH1) then
            csH1(lr) = 0.0
         else 
            cs0 = 5.475e-14
            fn0 = 4.298e-01
            y0  = 0.0
            y1  = 0.0
            yw  = 0.0
            ya  = 3.288e+01
            p   = 2.963
            x = fnu/fn0 - y0
            y = sqrt(x**2+y1**2)
            csH1(lr) = cs0*((x-1.0)**2+yw**2)*y**(0.5*p-5.5)/
     .           (1.0+sqrt(y/ya))**p/CSSTAR
         endif
          
         if(Txi(lr) .lt. TG1) then
            csG1(lr) = 0.0
         else
            cs0 = 9.492e-16
            fn0 = 1.361e+01
            y0  = 4.434e-01 
            y1  = 2.136
            yw  = 2.039
            ya  = 1.469
            p   = 3.188
            x = fnu/fn0 - y0
            y = sqrt(x**2+y1**2)
            csG1(lr) = cs0*((x-1.0)**2+yw**2)*y**(0.5*p-5.5)/
     .           (1.0+sqrt(y/ya))**p/CSSTAR
         endif

         if(Txi(lr) .lt. TG2) then
            csG2(lr) = 0.0
         else
            cs0 = 1.369e-14
            fn0 = 1.720
            y0  = 0.0 
            y1  = 0.0
            yw  = 0.0
            ya  = 3.288e+01
            p   = 2.963
            x = fnu/fn0 - y0
            y = sqrt(x**2+y1**2)
            csG2(lr) = cs0*((x-1.0)**2+yw**2)*y**(0.5*p-5.5)/
     .           (1.0+sqrt(y/ya))**p/CSSTAR
         endif

         w = 1.24/fnu
C
C  SMC dust - from Weingartner & Draine 2001, ApJ 548, 296
C  Using fitting functions from Pei 1992, ApJ 395, 130
C  My fudging: 5.50->15.50 in term #2, plus my term #7
C  Assuming ZSMC = -0.7dex from Welty et al 1997, ApJ 489, 672
C
         csDust(lr) = 5.6e-22/CSSTAR*(
     .        dustfit(w, 0.042, 185.0  , 90.0 , 2.0, 2.0) +
     .        dustfit(w, 0.08 ,  27.0  , 15.50, 4.0, 4.0) +  !! 15.5, not 5.5
     .        dustfit(w, 0.22 ,   0.005, -1.95, 2.0, 2.0) +
     .        dustfit(w, 9.7  ,   0.010, -1.95, 2.0, 2.0) +
     .        dustfit(w,18.0  ,   0.012, -1.80, 2.0, 2.0) +
     .        dustfit(w,25.0  ,   0.030,  0.0 , 2.0, 2.0) )
         if(w .gt. 1.0e-3) then
            csDust(lr) = 5.6e-22/CSSTAR*(
     .           dustfit(w, 0.067 ,  10.0 ,  1.90, 4.0, 15.0) )
         endif
C
C  LMC dust
C  My fudging: 175->90 in term #1, 5.5->21 in term #2, plus my term #7
C  The same term #7 for both SMC and LMC
C  Assuming ZLMC=-0.5dex from Welty et al 1999, ApJ, 512, 636 
C
*         csDust(lr) = 1.1e-21/CSSTAR*(
*     .        dustfit(w, 0.046, 90.0  , 90.0 , 2.0, 2.0) +
*     .        dustfit(w, 0.08 ,  19.0  , 21.0, 4.5, 4.5) +  !! 21, not 5.5
*     .        dustfit(w, 0.22 ,   0.023, -1.95, 2.0, 2.0) +
*     .        dustfit(w, 9.7  ,   0.005, -1.95, 2.0, 2.0) +
*     .        dustfit(w,18.0  ,   0.006, -1.80, 2.0, 2.0) +
*     .        dustfit(w,25.0  ,   0.020,  0.0 , 2.0, 2.0) )
*         if(w .gt. 1.0e-3) then
*            csDust(lr) = 1.1e-21/CSSTAR*(
*     .           dustfit(w, 0.067 ,  10.0 ,  1.90, 4.0, 15.0) )
*         endif

#ifdef RT_CHEMISTRY

         if(fnu .lt. 0.755) then
            csMH27(lr) = 0.0
         else
            csMH27(lr) = 2.16e-16*(fnu-0.755)**1.5/fnu**3/CSSTAR
         endif

         if(fnu .lt. 2.65) then
            csMH28(lr) = 0.0
         else if(fnu .lt. 11.27) then
            csMH28(lr) = 10.0**(-40.97+6.03*fnu-0.504*fnu**2+1.387e-2*
     .           fnu**3)/CSSTAR
         else if(fnu .lt. 21.0) then
            csMH28(lr) = 10.0**(-30.26+2.79*fnu-0.184*fnu**2+3.535e-3*
     .           fnu**3)/CSSTAR
         else
            csMH28(lr) = 8.39e-21*(21.0/fnu)**3/CSSTAR
         endif

         if(fnu .lt. 15.42) then
            csMH29(lr) = 0.0
         else if(fnu .lt. 16.5) then
            csMH29(lr) = (6.2e-18*fnu-9.4e-17)/CSSTAR
         else if(fnu .lt. 17.7) then
            csMH29(lr) = (1.4e-18*fnu-1.48e-17)/CSSTAR
         else 
            csMH29(lr) = 2.5e-14/fnu**2.71/CSSTAR
         endif

         if(fnu .lt. 30.0) then
            csMH30(lr) = 0.0
         else if(fnu .lt. 70.0) then
            csMH30(lr) = 10.0**(-16.926-4.528e-2*fnu+2.238e-4*fnu**2+
     .           4.245e-7*fnu**3)/CSSTAR
         else
            csMH30(lr) = 4.0e-14/fnu**3/CSSTAR
         endif

         if(fnu .lt. 14.675) then
            csMH31(lr) = 0.0
         else 
            if(fnu .lt. 16.820) then
               csMH31(lr) = 10.0**(-18+15.129-1.051*fnu)
            else
               csMH31(lr) = 10.0**(-18-31.41+1.804e-2*fnu**3-
     .              4.234e-5*fnu**5)
            endif
            csMH31(lr) = csMH31(lr) + 10.0**(-18+13.531-0.918*fnu)
         endif

         if(fnu .lt. 14.159) then
            csMH32(lr) = 0.0
         else 
            if(fnu .lt. 15.302) then
               csMH32(lr) = 10.0**(-18+12.022-1.051*fnu)
            else
               csMH32(lr) = 10.0**(-18+16.046-1.0824*fnu)
            endif
            csMH32(lr) = csMH32(lr) + 10.0**(-18+12.874-0.851*fnu)
         endif

         if(fnu.lt.11.3 .or. fnu.gt.13.6) then
            csMHLW(lr) = 0.0
         else
C
C  R = 10^8 s^{-1} I_\nu  (that is what Tom Abel told me)
C
            csMHLW(lr) = 3.6e-18/CSSTAR
         endif

#endif

      enddo

#ifdef RT_XRAYS

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,s),
C$OMP+SHARED(Txi,
C$OMP+ se1H1i,se2H1i,se1H1h,se2H1h,se1G1i,se2G1i,se1G1h,se2G1h)
      do lr = 1,NRAD

         s = (28.0*T1eV/max(EM5,Txi(lr)-TH1))**0.4
         if(s .gt. 1.0) then
            se1H1i(lr) = 0.0
            se2H1i(lr) = 0.0
         else
            se1H1i(lr) = 1.0
            se2H1i(lr) = s
         endif
               
         s = (11.0*T1eV/max(EM5,Txi(lr)-TH1))**0.7
         if(s .gt. 1.0) then
            se1H1h(lr) = 0.0
            se2H1h(lr) = 0.0
         else
            se1H1h(lr) = 1.0
            se2H1h(lr) = s
         endif
               
         s = (28.0*T1eV/max(EM5,Txi(lr)-TG1))**0.4
         if(s .gt. 1.0) then
            se1G1i(lr) = 0.0
            se2G1i(lr) = 0.0
         else
            se1G1i(lr) = 1.0
            se2G1i(lr) = s
         endif
               
         s = (11.0*T1eV/max(EM5,Txi(lr)-TG1))**0.7
         if(s .gt. 1.0) then
            se1G1h(lr) = 0.0
            se2G1h(lr) = 0.0
         else
            se1G1h(lr) = 1.0
            se2G1h(lr) = s
         endif
               
      enddo

#endif
C
      csH1AB = 1.18e-16*sqrt(1.0)*0.416/10.2/CSSTAR
      csG1AB = 1.18e-16*sqrt(4.0)*0.285/21.2/CSSTAR
      csG2AB = 1.18e-16*sqrt(4.0)*0.416/40.8/CSSTAR
C
      lrH1tr = 0
      lrG1tr = 0
      lrG2tr = 0
      do lr=1,NRAD
         if(lrH1tr.eq.0 .and. csH1(lr).gt.1.0e-35) lrH1tr = lr
         if(lrG1tr.eq.0 .and. csG1(lr).gt.1.0e-35) lrG1tr = lr
         if(lrG2tr.eq.0 .and. csG2(lr).gt.1.0e-35) lrG2tr = lr
      enddo
      if(lrH1tr .eq. 0) lrH1tr = 1
      if(lrG1tr .eq. 0) lrG1tr = 1
      if(lrG2tr .eq. 0) lrG2tr = 1
      csH1tr = csH1(lrH1tr)
      csG1tr = csG1(lrG1tr)
      csG2tr = csG2(lrG2tr)

      lrH1rB = 0
      lrG1rB = 0
      lrG2rB = 0
      do lr=1,NRAD
         if(lrH1rB.eq.0 .and. Txi(lr).gt.0.75*TH1) lrH1rB = lr
         if(lrG1rB.eq.0 .and. Txi(lr).gt.0.86*TG1) lrG1rB = lr
         if(lrG2rB.eq.0 .and. Txi(lr).gt.0.75*TG2) lrG2rB = lr
      enddo
      if(lrH1rB .eq. 0) lrH1rB = 1
      if(lrG1rB .eq. 0) lrG1rB = 1
      if(lrG2rB .eq. 0) lrG2rB = 1
C
#ifdef RT_DUST
      lrDUtr = 0
      do lr=1,NRAD
         if(lrDUtr.eq.0 .and. csDust(lr).gt.1.0e-35) lrDUtr = lr
      enddo
      if(lrDUtr .eq. 0) lrDUtr = NRAD
#else
      lrDUtr = NRAD
#endif
C
#ifdef RT_XRAYS
      lrSEtr = 0
      do lr=1,NRAD
         if(lrSEtr.eq.0 .and. Txi(lr).gt.2*TH1) lrSEtr = lr
      enddo
      if(lrSEtr .eq. 0) lrSEtr = NRAD
#else
      lrSEtr = NRAD
#endif
C
      lrH1co = 0
      do lr=NRAD,lrH1tr,-1
         if(lrH1co.eq.0 .and. csH1n(lr).gt.CSNMIN) lrH1co = lr
      enddo
      if(lrH1co .eq. 0) lrH1co = NRAD

      lrG1co = 0
      do lr=NRAD,lrG1tr,-1
         if(lrG1co.eq.0 .and. csG1n(lr).gt.CSNMIN) lrG1co = lr
      enddo
      if(lrG1co .eq. 0) lrG1co = NRAD

      lrG2co = 0
      do lr=NRAD,lrG2tr,-1
         if(lrG2co.eq.0 .and. csG2n(lr).gt.CSNMIN) lrG2co = lr
      enddo
      if(lrG2co .eq. 0) lrG2co = NRAD

      lrmax = max(lrH1co,lrG1co,lrG2co)
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(csH1,csG1,csG2,csH1n,csG1n,csG2n,csH1tr,csG1tr,csG2tr,
C$OMP+ csDust,csDustn)
      do lr=1,NRAD
         csH1n(lr) = csH1(lr)/csH1tr
         csG1n(lr) = csG1(lr)/csG1tr
         csG2n(lr) = csG2(lr)/csG2tr
         csDustn(lr) = csDust(lr)/csH1tr
      enddo

#ifdef RT_CHEMISTRY

      lr27tr = 0
      lr28tr = 0
      lr29tr = 0
      lr30tr = 0
      lr31tr = 0
      lr32tr = 0
      lrLWtr = 0
      do lr=1,NRAD
         if(lr27tr.eq.0 .and. csMH27(lr).gt.1.0e-35) lr27tr = lr
         if(lr28tr.eq.0 .and. csMH28(lr).gt.1.0e-35) lr28tr = lr
         if(lr29tr.eq.0 .and. csMH29(lr).gt.1.0e-35) lr29tr = lr
         if(lr30tr.eq.0 .and. csMH30(lr).gt.1.0e-35) lr30tr = lr
         if(lr31tr.eq.0 .and. csMH31(lr).gt.1.0e-35) lr31tr = lr
         if(lr32tr.eq.0 .and. csMH32(lr).gt.1.0e-35) lr32tr = lr
         if(lrLWtr.eq.0 .and. csMHLW(lr).gt.1.0e-35) lrLWtr = lr
      enddo
      if(lr27tr .eq. 0) lr27tr = 1
      if(lr28tr .eq. 0) lr28tr = 1
      if(lr29tr .eq. 0) lr29tr = 1
      if(lr30tr .eq. 0) lr30tr = 1
      if(lr31tr .eq. 0) lr31tr = 1
      if(lr32tr .eq. 0) lr32tr = 1
      if(lrLWtr .eq. 0) lrLWtr = 1

#ifdef RT_LWBANDS
C
C  Lyman-Werner band with its own frequency dependence
C
      XILWST = log(fnuLW2/fnuLW1)/real(NRLW-1)
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(xigLW,angLW,XILWST)
      do lr=1,NRLW
         xigLW(lr) = XILWST*(lr-1)
         angLW(lr) = 0.0
      enddo
#endif
#endif
C
C  Temperature table
C
      acLTmin = log(acTmin)
      acLTmax = log(acTmax)
      acLTstp = (acLTmax-acLTmin)/(real(NTEM)-1)

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lc,tem,altem,q,p,alt10,srt,alt,alH1,alG1,alG2,s,
C$OMP+ tev,x,itl,itu,rC1a,rC2e,rC2a,rN1e,
C$OMP+ rO1a,rMg2e,rSi2e,rSi2a,rS2e,rFe2e,rFe2a,rO1e)
C$OMP+SHARED(acLTmin,acLTstp,acLT,acTEM,acSQT,RateAtom,RateChem,
C$OMP+ dsd93,nsd93,sd93lt,sd93c0,sd93c1,RateHigh)
      do lc = 1,NTEM

         altem = acLTmin + acLTstp*real(lc-1)
         acLT(lc) = altem

         tem = exp(altem)
         srt = sqrt(tem)

         acTEM(lc) = tem
         acSQT(lc) = srt

         alt10 = log10(tem)
C
C  Atomic rates
C
         alH1 = 2.0*TH1/tem
         alG1 = 2.0*TG1/tem
         alG2 = 2.0*TG2/tem
         q = 1.0/tem**1.5
         p = 1.0/(1.0+sqrt(tem/1.0e5))
         alt = alt10 - 1.0

         if(alH1 .gt. 2.0e2) then
            s = 0.0
         else 
            s = 21.11*q*exp(-0.5*alH1)*
     .           alH1**(-1.089)/(1.0+(alH1/0.354)**0.874)**1.101
         endif
         RateAtom(irateCIcH1,lc) = (1.381e-16/CFSTAR)*TH1*s
         RateAtom(irateCIiH1,lc) = s

         if(alG1 .gt. 2.0e2) then
            s = 0.0
         else 
            s = 32.38*q*exp(-0.5*alG1)*
     .           alG1**(-1.146)/(1.0+(alG1/0.416)**0.987)**1.056
         endif
         RateAtom(irateCIcG1,lc) = (1.381e-16/CFSTAR)*TG1*s
         RateAtom(irateCIiG1,lc) = s

         if(alG2 .gt. 2.0e2) then
            s = 0.0
         else 
            s = 19.95*q*exp(-0.5*alG2)*
     .           alG2**(-1.089)/(1.0+(alG2/0.553)**0.735)**1.275
         endif
         RateAtom(irateCIcG2,lc) = (1.381e-16/CFSTAR)*TG2*s
         RateAtom(irateCIiG2,lc) = s
C
C  Recombination case A
C
         RateAtom(irateRAiH2,lc) = 1.269e-13*alH1**1.503/
     .        (1.0+(alH1/0.522)**0.470)**1.923
         RateAtom(irateRAcH2,lc) = (1.778e-29/CFSTAR)*tem*alH1**1.965/
     .        (1.0+(alH1/0.541)**0.502)**2.697
C
C  HeII rates are from Hummer & Storey, MNRAS, 297, 1073 with
C  correction for altem>4.4 by Massimo Ricotti
C
         if(alt .lt. 3.4) then
            RateAtom(irateRAiG2,lc) = 1.0e-11/srt*(10.85-2.34*alt+
     .           0.019*alt**3)
            RateAtom(irateRAcG2,lc) = 1.0e-11/srt*( 9.92-2.37*alt+
     .           0.030*alt**3)*(1.381e-16/CFSTAR)*tem
         else
            RateAtom(irateRAiG2,lc) = 2.31e-10/tem**0.681
            RateAtom(irateRAcG2,lc) = 1.96e-10/tem**0.669*
     .           (1.381e-16/CFSTAR)*tem
         endif

         RateAtom(irateRAiG3,lc) = 2.0*1.269e-13*alG2**1.503/
     .        (1.0+(alG2/0.522)**0.470)**1.923
         RateAtom(irateRAcG3,lc) = 8.0*(1.778e-29/CFSTAR)*tem*
     .        alG2**1.965/(1.0+(alG2/0.541)**0.502)**2.697
C
C  Recombination case B
C
         RateAtom(irateRBiH2,lc) = 2.753e-14*alH1**1.500/
     .        (1.0+(alH1/2.740)**0.407)**2.242
         RateAtom(irateRBcH2,lc) = (3.435e-30/CFSTAR)*tem*alH1**1.970/
     .        (1.0+(alH1/2.250)**0.376)**3.720

         if(alt .lt. 3.4) then
            RateAtom(irateRBiG2,lc) = 1.0e-11/srt*(9.28-2.31*alt+
     .           0.016*alt**3)
            RateAtom(irateRBcG2,lc) = 1.0e-11/srt*(8.35-2.37*alt+
     .           0.028*alt**3)*(1.381e-16/CFSTAR)*tem
         else
            RateAtom(irateRBiG2,lc) = 5.13e-10/tem**0.817
            RateAtom(irateRBcG2,lc) = 5.93e-10/tem**0.870*
     .           (1.381e-16/CFSTAR)*tem
         endif

         RateAtom(irateRBiG3,lc) = 2.0*2.753e-14*alG2**1.500/
     .        (1.0+(alG2/2.740)**0.407)**2.242
         RateAtom(irateRBcG3,lc) = 8.0*(3.435e-30/CFSTAR)*tem*
     .        alG2**1.970/
     .        (1.0+(alG2/2.250)**0.376)**3.720

         if(tem .lt. 4.7e+3) then
            s = 0.0
         else 
            s = exp(-4.7e+5/tem)*
     .           (1.0+0.3*exp(-9.4e+4/tem))/(tem*srt)
         endif
         RateAtom(irateDRcG2,lc) = (1.24e-13/CFSTAR)*s
         RateAtom(irateDRiG2,lc) = 1.90e-03*s

         if(tem .lt. 1183.0) then
            s = 0.0
         else 
            s = p*exp(-118348.0/tem)
         endif
         RateAtom(irateLCcH1,lc) = (7.50e-19/CFSTAR)*s
         
         if(tem .lt. 4736.0) then
            s = 0.0
         else 
            s = p*tem**(-0.397)*exp(-473638.0/tem)
         endif
         RateAtom(irateLCcG2,lc) = (5.54e-17/CFSTAR)*s
C
C  Chemical rates
C
#ifdef RT_CHEMISTRY

         tev = tem/T1eV
         alt = alt10

#if (RT_H2_RATE == 1)
C
C  Molecular Hydrogen (MH) rates (from Shapiro & Kang, 1987)
C
         if(tem .gt. 1.5e4) then
            RateChem(irateMHi07,lc) = 10.0**(-14.10+0.1175*alt-
     .           9.813e-3*alt**2)
         else
            RateChem(irateMHi07,lc) = 1.0e-18*tem
         endif

         if(tem .gt. 1.0e4) then
            RateChem(irateMHi08,lc) = 10.0**(-8.78+0.113*alt-
     .           3.475e-2*alt**2)
         else
            RateChem(irateMHi08,lc) = 1.3e-9
         endif

         if(tem .gt. 6.7e3) then
            RateChem(irateMHi09,lc) = 5.81e-16*(tem/56200.0)**
     .           (3.1619-0.6657*alt)
         else
            RateChem(irateMHi09,lc) = 1.85e-23*tem**1.8
         endif

         RateChem(irateMHi10,lc) = 6.0e-10
            
         RateChem(irateMHi11,lc) = 1.12e-10*exp(-7.035e4/tem)
            
         RateChem(irateMHi12,lc) = 2.4e-9*exp(-21200.0/tem)
            
         if(tem .gt. 2500.0) then
            RateChem(irateMHi13,lc) = 9.69e-13*exp(-11.323/
     .           (altem-7.28))
         else
            RateChem(irateMHi13,lc) = 0.0
         endif
C
C  This is Shapiro & Kang 98 rate
C
*         RateChem(irateMHi14,lc) = 4.38e-10*exp(-102000.0/tem)*
*     .          (tem/(1.0+tem/1.0e5))**0.35
C
C  This is Donahue & Shull 91 rate (a smaller one)
C
         RateChem(irateMHi14,lc) = 5.6e-11*exp(-102000.0/tem)*
     .        (tem/(1.0+tem/1.0e5))**0.5

         RateChem(irateMHi15,lc) = 1.18e-10*exp(-6.95e4/tem)

         RateChem(irateMHi16,lc) = 4.0e-12*exp(-8750/tem)*
     .        (2.0*tem/(1.0+tem/16000.0))

         RateChem(irateMHi17,lc) = 5.3e-20*exp(-8750/tem)*
     .        (2.0*tem/(1.0+tem/16000.0))**2.17
            
         RateChem(irateMHi18,lc) = 4.0e-6/srt

         if(tem .gt. 1.0e4) then
            RateChem(irateMHi19,lc) = 4.0e-4/tem**1.4*
     .           exp(-15100.0/tem)
         else
            RateChem(irateMHi19,lc) = 1.0e-8/tem**0.4
         endif

         RateChem(irateMHi20,lc) = 2.25e-6/tem**0.4
         
         RateChem(irateMHi21,lc) = 4.0e-6/srt

#define RT_H2_RATE_SET
#endif
#if (RT_H2_RATE == 2)
C
C  Molecular Hydrogen (MH) rates (from Abel et al.)
C
         if(tem .gt. 6.e3) then
            RateChem(irateMHi07,lc) = 3.802e-17*tem**(0.1998*
     .           log10(tem))*exp(4.0415e-5*(log10
     .           (tem))**6-5.447e-3*(log10(tem))**4)
         else
            RateChem(irateMHi07,lc) = 1.429e-18*tem**0.7620*
     .           tem**(0.1523*log10(tem))*
     .           tem**(-3.274e-2*(log10(tem))**2)
         endif

         if(tev .gt. 0.1) then
            RateChem(irateMHi08,lc) = exp(-20.06913897+0.22898*
     .           log(tev)+3.5998377e-2*(log(tev))**2-
     .           4.55512e-3*(log(tev))**3-3.10511544e-4*(log
     .           (tev))**4+1.0732940e-5*(log(tev))**5-
     .           8.36671960e-6*(log(tev))**6+2.23830623e-7*
     .           (log(tev))**7)
         else
            RateChem(irateMHi08,lc) = 1.428e-9
         endif
            
         if(tev .gt. 0.577) then
            RateChem(irateMHi09,lc) = 5.81e-16*(0.20651*
     .           tev)**(-0.2891*log10(0.20651*tev))
         else
            RateChem(irateMHi09,lc) = 3.833e-16*tev**1.8
         endif
            
         RateChem(irateMHi10,lc) = 6.4e-10

         RateChem(irateMHi11,lc) = 1.067e-10*tev**2.012*
     .        exp(-4.463/tev*(1.+0.2472*tev)**3.512)

         RateChem(irateMHi12,lc) = exp(-24.24914687+3.40082444*
     .        log(tev)-3.89800396*(log(tev))**2+2.04558782*
     .        (log(tev))**3-0.541618285*(log(tev))**4+
     .        8.41077503e-2*(log(tev))**5-7.87902615e-3*
     .        (log(tev))**6+4.13839842e-4*(log(tev))**7-
     .        9.36345888e-6*(log(tev))**8)

         RateChem(irateMHi13,lc) = 0.0

         RateChem(irateMHi14,lc) = 5.6e-11*tem**0.5*
     .        exp(-102124/tem)

         RateChem(irateMHi15,lc) = 0.0

         RateChem(irateMHi16,lc) = exp(-18.01849334+2.3608522*
     .        log(tev)-0.28274430*(log(tev))**2+1.62331664e-2*
     .        (log(tev))**3-3.36501203e-2*(log(tev))**4+
     .        1.17832978e-2*(log(tev))**5-1.65619470e-3*
     .        (log(tev))**6+1.06827520e-4*(log(tev))**7-
     .        2.63128581e-6*(log(tev))**8)
         
         if(tev .gt. 0.1) then
            RateChem(irateMHi17,lc) = exp(-20.37260896+
     .           1.13944933*log(tev)-0.14210135*(log(tev))**2+
     .           8.4644554e-3*(log(tev))**3-1.4327641e-3*
     .           (log(tev))**4+2.0122503e-4*(log(tev))**5+
     .           8.6639632e-5*(log(tev))**6-2.5850097e-5*
     .           (log(tev))**7+2.4555012e-6*(log(tev))**8-
     .           8.0683825e-8*(log(tev))**9)
         else
            RateChem(irateMHi17,lc) =2.5634e-9*tev**1.78186
         endif

         RateChem(irateMHi18,lc) = 7.e-8*(tem/1.d2)**(-0.5)

         if(tev .gt. 1.719) then
            RateChem(irateMHi19,lc) = 8.4258e-10*tev**(-1.4)*
     .           exp(-1.301/tev)
         else
            RateChem(irateMHi19,lc) = 2.291e-10*tev**(-0.4)
         endif

         if(tem .ge. 617) then
            RateChem(irateMHi20,lc) = 1.32e-6*tem**(-0.76)
         else
            RateChem(irateMHi20,lc) = 1.e-8
         endif

         RateChem(irateMHi21,lc) = 5.e-7*(1.e2/tem)**0.5

#define RT_H2_RATE_SET
#endif
#if (RT_H2_RATE == 3)
C
C  Molecular Hydrogen (MH) rates (from Galli & Palla)
C
         RateChem(irateMHi07,lc) = 1.4e-18*tem**0.928*
     .        exp(-tem/16200)
        
         if(tem .gt. 300.) then
            RateChem(irateMHi08,lc) = 4.0e-9*tem**(-0.17)
         else
            RateChem(irateMHi08,lc) = 1.5e-9
         endif

         RateChem(irateMHi09,lc) = 10.**(-19.38 - 1.523*alt + 
     .        1.118*alt**2 - 0.1269*alt**3)

         RateChem(irateMHi10,lc) = 6.4e-10

         RateChem(irateMHi11,lc) = 0.0

         if(tem .gt. 10000.0) then
            RateChem(irateMHi12,lc) = 1.5e-10*exp(-14000/tem)
         else
            RateChem(irateMHi12,lc) = 3.0e-10*exp(-21050/tem)
         endif

         RateChem(irateMHi13,lc) = 2.7e-8/tem**1.27*
     .        exp(-43000/tem)      

         RateChem(irateMHi14,lc) = 4.4e-10*tem**0.35*
     .        exp(-102000/tem)

         RateChem(irateMHi15,lc) = 0.0

         RateChem(irateMHi16,lc) = 0.0
            
         RateChem(irateMHi17,lc) = 0.0

         RateChem(irateMHi18,lc) = 5.7e-6/tem**0.5 + 6.3e-8 - 
     .        9.2e-11*tem**0.5 + 4.4e-13*tem

         if(tem .gt. 8000.0) then
            RateChem(irateMHi19,lc) = 9.6e-7/tem**0.9
         else
            RateChem(irateMHi19,lc) = 6.9e-9/tem**0.35
         endif

         RateChem(irateMHi20,lc) = 2.0e-7/tem**0.5
         
         RateChem(irateMHi21,lc) = 0.0

#define RT_H2_RATE_SET
#endif

#ifndef RT_H2_RATE_SET
#error "Incorrect value for RT_H2_RATE: must be 1, 2, or 3"
#endif

#if (RT_H2_COOL == 1)
C
C     Lepp & Shull cooling (something wrong)
C
         if(tem .gt. 1635.0) then
            RateChem(irateMHcLV,lc) = (8.18e-25/CFSTAR)*srt*
     .           exp(-1000.0/tem)
         else
            RateChem(irateMHcLV,lc) = (1.15e-25/CFSTAR)*
     .           exp(tem/125.0-(tem/577.0)**2)
         endif

         if(tem .gt. 4031.0) then
            RateChem(irateMHcLR,lc) = (1.66e-22/CFSTAR)*
     .           exp(-9243.0/tem)
         else
            x = log10(tem/1e4)
            RateChem(irateMHcLR,lc) = 1.2*10.0**(-22.9-0.553*x-
     .           1.148*x**2)/CFSTAR
         endif

#define RT_H2_COOL_SET
#endif
#if (RT_H2_COOL == 2)
C
C  Galli & Palla cooling (low density limit n<10^3 cm^-3)
C
         RateChem(irateMHcLV,lc) = 0.
         
         RateChem(irateMHcLR,lc) = 10.**(-103.0+97.59*alt-
     .        48.05*alt**2+10.80*alt**3-0.9032*alt**4)/CFSTAR
         
#define RT_H2_COOL_SET
#endif

#ifndef RT_H2_COOL_SET
#error "Incorrect value for RT_H2_COOL: must be 1 or 2"
#endif
C
C  Other H2 rates
C
         RateChem(irateMHcUE,lc) = 4.476*(1.602e-12/CFSTAR)*
     .        (RateChem(irateMHi13,lc)+RateChem(irateMHi14,lc))
      
         RateChem(irateMHcUH,lc) = 4.476*(1.602e-12/CFSTAR)*
     .        RateChem(irateMHi11,lc)

         RateChem(irateMHcUM,lc) = 4.476*(1.602e-12/CFSTAR)*
     .        RateChem(irateMHi15,lc)

         RateChem(irateMHc16,lc) = 0.75*(1.602e-12/CFSTAR)*
     .        RateChem(irateMHi16,lc)

         RateChem(irateMHc17,lc) = 0.75*(1.602e-12/CFSTAR)*
     .        RateChem(irateMHi17,lc)
C
C  heavy elements (HE) rates
C
C  1. SD93 re-interpolation
C
         itl = 1 + int((alt-sd93lt(1))/dsd93)
         if(itl .lt. 1) then
            RateChem(irateZHcTA,lc) = 0.0
            RateChem(irateZHcTB,lc) = 0.0
            RateChem(irateZHcTC,lc) = 0.0
         else
            if(itl .ge. nsd93) itl = nsd93 - 1
            itu = itl + 1
            RateChem(irateZHcTA,lc) = (sd93c0(itl) + 
     .           (sd93c0(itu)-sd93c0(itl))/dsd93*(alt-sd93lt(itl)))/
     .           CFSTAR
            RateChem(irateZHcTB,lc) = (sd93c1(itl) + 
     .           (sd93c1(itu)-sd93c1(itl))/dsd93*(alt-sd93lt(itl)))/
     .           CFSTAR
            RateChem(irateZHcTC,lc) = log(1.0+tem/1.0e8)
         endif
C
C  2. low T regime: build the CF from the data from Penston, ApJ, 162, 771
C  and Dalgarno & McCray, ARAA, 10, 375
C
C  1. Collisions with electrons
C
         rC2e = 1.22e-17/tem**0.5*exp(-61900/tem) + 
     .        6.67e-20/tem**0.5*exp(-92/tem)         !! D&MC: 7.9
         rN1e = 1.09e-21*tem**0.5*exp(-27700/tem)    !! D&MC: 0.82
         rO1e = 1.40e-22*tem**0.5*exp(-22500/tem) +  !! D&MC: 0.94
     .        1.17e-24*tem**0.333*exp(-326/tem)      !! D&MC: 0.76
         rMg2e = 6.08e-17/tem**0.5*exp(-51300/tem)
         rSi2e = 2.36e-17/tem**0.5*exp(-79600/tem) + 
     .        1.67e-18/tem**0.5*exp(-413/tem)        !! D&MC: 1.9
         rS2e = 8.44e-18/tem**0.5*exp(-21400/tem)
         rFe2e = 1.50e-18/tem**0.5*exp(-961/tem) +   !! D&MC: 1.3
     .        1.05e-18/tem**0.5*exp(-554/tem)        !! D&MC: 1.1
C
C  2. Collisions with H atoms
C
         rC1a = 3.4e-24*tem**0.333*exp(-62/tem) + 
     .        0.81e-24*tem**0.333*exp(-23/tem)
         rO1a = 2.53e-24*tem**0.333*exp(-326/tem) + 
     .        6.43e-24*tem**0.333*exp(-228/tem)
         rFe2a = 1.41e-22*exp(-961/tem) +            !! D&MC: 1.54 
     .        1.38e-22*exp(-554/tem)                 !! D&MC: 1.1
C
C  The next two are from Barinovs et al, ApJ, 620, 537
C
         rC2a = 1.0e-24*exp(-92/tem)*max(0.0,16+0.344*sqrt(tem)-
     .        27.7/tem)
         rSi2a = 1.0e-24*exp(-413/tem)*(43.5+1.78*sqrt(tem)+
     .        0.005*tem)
C
C  The suppresion factor is somewhat arbitrary. I am not sure that
C  there rates are valid even at T=1.0e4. The supression factor will
C  make sure that SD93 rates are used above 1.0e4, not these ones.
C
         RateChem(irateZLcTA,lc) = (3.55e-5*rSi2a + 3.31e-4*rC2a + 
     .        6.76e-4*rO1a + 3.16e-5*rFe2a + 3.31e-4*rC1a)/
     .        (1.0+(tem/2.0e4)**4)/CFSTAR
         RateChem(irateZLcTx,lc) = (3.55e-5*rSi2e + 3.31e-4*rC2e + 
     .        3.80e-5*rMg2e + 8.32e-5*rN1e + 6.76e-4*rO1e + 
     .        2.14e-5*rS2e +  3.16e-5*rFe2e)/
     .        (1.0+(tem/2.0e4)**4)/CFSTAR

C
C  Dust: from Draine, ApJ, 245, 880
C
         RateChem(irateMHcDU,lc) = (2.0e-22/CFSTAR)*(tem/1.0e6)**1.5/
     .        (1+(tem/1.36e7)**1.5)
C
C  H2 formation on dust: from Cazauz & Spaans (ApJ, 611, 40)
C
*         RateChem(irateMHiDU,lc) = 3.5e-18*sqrt(tem/100.0)/
*     .        (1+0.4*sqrt((tem+d_temGr)/100)+0.2*tem/100+
*     .        0.08*(tem/100)**2)
C
C  H2 formation on dust: empirical rate from Wolfire et al (arXiv:0803.0138)
C
         RateChem(irateMHiDU,lc) = 3.5e-17
C
C  High density regime
C
#ifdef RT_HIGH_DENSITY

         alt = alt10 - 4
C
C  Shapiro & Kang 87
C
         RateHigh(irateHDi11,lc) = 1.20e-9*exp(-5.24e4/tem)*
     .        10.0**(4.0-0.416*alt-0.327*alt**2)

         RateHigh(irateHDi15,lc) = 1.30e-9*exp(-5.33e4/tem)*
     .        10.0**(4.845-1.3*alt+1.62*alt**2)

         RateHigh(irateHDi22,lc) = 5.5e-29/tem

         RateHigh(irateHDc08,lc) = 1.0e6/srt*(1.4+
     .        1.6*exp(-(400/tem)**2))/CFSTAR
C
C  Hollenbach & McKee, 1979, ApJS, 41, 555
C
         RateHigh(irateHDcLV,lc) = (6.7e-19*exp(-5.86e3/tem) +
     .        1.6e-18*exp(-11.7e3/tem))/CFSTAR

         RateHigh(irateHDcLR,lc) = (9.5e-22*(tem/1.0e3)**3.76/
     .        (1+0.12*(tem/1.0e3)**2.1)*exp(-(0.13e3/tem)**3) +
     .        3.0e-24*exp(-0.51e3/tem))/CFSTAR

#endif  // RT_HIGH_DENSITY
#endif  // RT_CHEMISTRY

      enddo
C
C  Radiation field tables
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lo1,lo2)
C$OMP+SHARED(itloc2)
      do lo1=1,NOPT1
         do lo2=1,NOPT1
            itloc2(lo2,lo1) = NOPT1*(lo2-1+NOPT1*(lo1-1))
         enddo
      enddo
C
C  Optical depth table (using nonuniform mapping)
C
      acLOmin = log(acOmin)
      acLOmax = log(acOmax)
      acLOstp = (acLOmax-acLOmin)/real(NOPT-2)
      do lo = 2,NOPT
         acLO(lo) = acLOmin+acLOstp*real(lo-2)
         acOPT(lo) = exp(acLO(lo))
      enddo
      acLO(1) = -1.0e+2
      acOPT(1) = 0.0

      acLO100 = log(1.0e2)

      acTOstp = (acLOmax-acLOmin)/real(NTOL-1)
      acT2O(1) = acLOmin
      acO2T(1) = acLOmin
      acT2O(NTOL) = acLOmax
      acO2T(NTOL) = acLOmax
      toalp = 6.0
      tobet = 0.7
      toenn = 0.5

      do lo=2,NTOL-1

         altaul = acLOmin + acTOstp*real(lo-1)

         q = abs(altaul)
         if(q .gt. acLO100) then
            toy = 0.0
         else
            toy = (1.0-q/acLO100)*
     .           (toalp*q/(1.0+tobet*q))**toenn
         endif

         if(altaul .gt. 0.0) then
            aloptl = altaul + toy
         else
            aloptl = altaul - toy
         endif
         
         if(aloptl .lt. 1.00001*acT2O(lo-1)) then
            write(0,*) 'RT_Tables: fatal error, ',
     .           'opacity table is not monotonic.'
            call frtAbort
         endif

         acT2O(lo) = aloptl

      enddo

      lol = 1
      do lo=2,NTOL-1
         aloptl = acLOmin + acTOstp*real(lo-1)
 10      lor = lol + 1
         if(acT2O(lor) .lt. aloptl) then
            lol = lor
            goto 10
         endif
         altaul = acLOmin + acTOstp*real(lol-1) +
     .        acTOstp/(acT2O(lor)-acT2O(lol))*(aloptl-acT2O(lol))
         acO2T(lo) = altaul
      enddo

      acTAU(1) = acOPT(1)
      acTAU(NOPT) = acOPT(NOPT)
      do lo=2,NOPT-1
         alopt = acLO(lo)
         loptl = int((alopt-acLOmin)/acTOstp*0.99999)+1
         loptu = loptl + 1
         altau = acO2T(loptl) + (acO2T(loptu)-acO2T(loptl))*
     .        (alopt-acLOmin-acTOstp*real(loptl-1))/acTOstp
         acTAU(lo) = exp(altau)
      enddo
C
C  Oscillator strengths for all Ly transition above 48
C  (for Ly-series absorption)
C
      sum1 = 0.0
      do l=1,48
         sum1 = sum1 + oscstr(l)
      enddo
      oscstr_other = 1.0 - sum1
C
C  Integration weights
C
      do lr=1,NRAD
         wpTab(iptabPiH1,lr) = csH1(lr)
         wpTab(iptabPhH1,lr) = csH1(lr)*(Txi(lr)-TH1)
         wpTab(iptabPiG1,lr) = csG1(lr)
         wpTab(iptabPhG1,lr) = csG1(lr)*(Txi(lr)-TG1)
         wpTab(iptabPiG2,lr) = csG2(lr)
         wpTab(iptabPhG2,lr) = csG2(lr)*(Txi(lr)-TG2)
#ifdef RT_CHEMISTRY
         wpTab(iptabCi27,lr) = csMH27(lr)
         wpTab(iptabCi28,lr) = csMH28(lr)
         wpTab(iptabCi29,lr) = csMH29(lr)
         wpTab(iptabCi30,lr) = csMH30(lr)
         wpTab(iptabCi31,lr) = csMH31(lr)
         wpTab(iptabCi32,lr) = csMH32(lr)
         wpTab(iptabCiLW,lr) = csMHLW(lr)
#endif  // RT_CHEMISTRY
      enddo
 
#ifdef RT_XRAYS
      call frtInitTablesSE
#endif

      return
      end
C
C
C
      subroutine frtFillRadiationTables
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
C
      dimension pOT(iptabDim), sgpOT(iptabDim)
      dimension pH1(iptabDim), sgpH1(iptabDim)
      dimension pG1(iptabDim), sgpG1(iptabDim)
      dimension pG2(iptabDim), sgpG2(iptabDim)
C
      do lr=1,NRAD
         angEff(lr) = angUni(lr) + angLoc(lr)*QFuni(tauLoc(lr))
      enddo
C
      iptabMax = 6

      sumc1 = 0.0
      sumc2 = 0.0
      sumc3 = 0.0
      sumi1 = 0.0
      sumi2 = 0.0
      sumi3 = 0.0

      do lr=1,NRAD
         sumc1 = sumc1 + angAvg(lr)*csH1(lr)*(Txi(lr)-TH1)
         sumc2 = sumc2 + angAvg(lr)*csG1(lr)*(Txi(lr)-TG1)
         sumc3 = sumc3 + angAvg(lr)*csG2(lr)*(Txi(lr)-TG2)
         sumi1 = sumi1 + angAvg(lr)*csH1(lr)
         sumi2 = sumi2 + angAvg(lr)*csG1(lr)
         sumi3 = sumi3 + angAvg(lr)*csG2(lr)
      enddo

      pRate0(iptabPhH1) = max(0.0,facPhot*sumc1)
      pRate0(iptabPhG1) = max(0.0,facPhot*sumc2)
      pRate0(iptabPhG2) = max(0.0,facPhot*sumc3)
      pRate0(iptabPiH1) = max(0.0,facPhot*sumi1)
      pRate0(iptabPiG1) = max(0.0,facPhot*sumi2)
      pRate0(iptabPiG2) = max(0.0,facPhot*sumi3)

#ifdef RT_CHEMISTRY

      iptabMax = 12

      sumi27 = 0.0
      sumi28 = 0.0
      sumi29 = 0.0
      sumi30 = 0.0
      sumi31 = 0.0
      sumi32 = 0.0

      do lr=1,NRAD
         sumi27 = sumi27 + angAvg(lr)*csMH27(lr)
         sumi28 = sumi28 + angAvg(lr)*csMH28(lr)
         sumi29 = sumi29 + angAvg(lr)*csMH29(lr)
         sumi30 = sumi30 + angAvg(lr)*csMH30(lr)
         sumi31 = sumi31 + angAvg(lr)*csMH31(lr)
         sumi32 = sumi32 + angAvg(lr)*csMH32(lr)
      enddo

      pRate0(iptabCi27) = max(0.0,facPhot*sumi27)
      pRate0(iptabCi28) = max(0.0,facPhot*sumi28)
      pRate0(iptabCi29) = max(0.0,facPhot*sumi29)
      pRate0(iptabCi30) = max(0.0,facPhot*sumi30)
      pRate0(iptabCi31) = max(0.0,facPhot*sumi31)
      pRate0(iptabCi32) = max(0.0,facPhot*sumi32)
C
C  Do Lyman-Werner bands separately
C
      sumiLW = 0.0
#ifdef RT_LWBANDS
      do lr=1,NRLW
         sumiLW = sumiLW + angLW(lr)*csLWMH(lr)
      enddo
      sumiLW = sumiLW*XILWST/XISTEP
#else
      do lr=1,NRAD
         sumiLW = sumiLW + angAvg(lr)*csMHLW(lr)
      enddo
#endif  // RT_LWBANDS
      pRate0(iptabCiLW) = max(0.0,facPhot*sumiLW)

#endif  // RT_CHEMISTRY

#ifdef RT_TRANSFER

      sumc1 = 0.0
      sumc2 = 0.0
      sumc3 = 0.0
      sumi1 = 0.0
      sumi2 = 0.0
      sumi3 = 0.0

      do lr=1,NRAD
         sumc1 = sumc1 + guvEff(lr)*csH1(lr)*(Txi(lr)-TH1)
         sumc2 = sumc2 + guvEff(lr)*csG1(lr)*(Txi(lr)-TG1)
         sumc3 = sumc3 + guvEff(lr)*csG2(lr)*(Txi(lr)-TG2)
         sumi1 = sumi1 + guvEff(lr)*csH1(lr)
         sumi2 = sumi2 + guvEff(lr)*csG1(lr)
         sumi3 = sumi3 + guvEff(lr)*csG2(lr)
      enddo

      sgpRate0(iptabPhH1) = max(0.0,facPhot*sumc1*facLX)
      sgpRate0(iptabPhG1) = max(0.0,facPhot*sumc2*facLX)
      sgpRate0(iptabPhG2) = max(0.0,facPhot*sumc3*facLX)
      sgpRate0(iptabPiH1) = max(0.0,facPhot*sumi1*facLX)
      sgpRate0(iptabPiG1) = max(0.0,facPhot*sumi2*facLX)
      sgpRate0(iptabPiG2) = max(0.0,facPhot*sumi3*facLX)

#ifdef RT_CHEMISTRY

      sumi27 = 0.0
      sumi28 = 0.0
      sumi29 = 0.0
      sumi30 = 0.0
      sumi31 = 0.0
      sumi32 = 0.0

      do lr=1,NRAD
         sumi27 = sumi27 + guvEff(lr)*csMH27(lr)
         sumi28 = sumi28 + guvEff(lr)*csMH28(lr)
         sumi29 = sumi29 + guvEff(lr)*csMH29(lr)
         sumi30 = sumi30 + guvEff(lr)*csMH30(lr)
         sumi31 = sumi31 + guvEff(lr)*csMH31(lr)
         sumi32 = sumi32 + guvEff(lr)*csMH32(lr)
      enddo

      sgpRate0(iptabCi27) = max(0.0,facPhot*sumi27*facLX)
      sgpRate0(iptabCi28) = max(0.0,facPhot*sumi28*facLX)
      sgpRate0(iptabCi29) = max(0.0,facPhot*sumi29*facLX)
      sgpRate0(iptabCi30) = max(0.0,facPhot*sumi30*facLX)
      sgpRate0(iptabCi31) = max(0.0,facPhot*sumi31*facLX)
      sgpRate0(iptabCi32) = max(0.0,facPhot*sumi32*facLX)
C
C  Do Lyman-Werner bands separately
C
      sumiLW = 0.0
#ifdef RT_LWBANDS
      do lr=1,NRLW
         sumiLW = sumiLW + guvLWe(lr)*csLWMH(lr)
      enddo
      sumiLW = sumiLW*XILWST/XISTEP
#else
      do lr=1,NRAD
         sumiLW = sumiLW + guvEff(lr)*csMHLW(lr)
      enddo
#endif  // RT_LWBANDS
      sgpRate0(iptabCiLW) = max(0.0,facPhot*sumiLW*facLX)

#endif  // RT_CHEMISTRY

#endif  // RT_TRANSFER
C
C  Lyman-series absorption, averaged over a fraction of unit redshift
C  (currently the XISTEP: over one bin)
C
      do lr=1,NRAD
         csLyS(lr) = 0.0
      enddo
C
C  Switch it off. Not clear if it works at all. 
C  Not a continuous medium limit.
C
*      do l=2,NLYSS
*
*         fnuL = fnuH1*(1.0-1.0/real(l+1)**2)
*         cs2 = 1.098e-16*oscstr(l)/CSSTAR
*
*         do lr=2,lrH1tr-1
*            fnu1 = 0.5*(Txi(lr)+Txi(lr-1))/T1eV
*            fnu2 = 0.5*(Txi(lr)+Txi(lr+1))/T1eV
*            if(fnuL.gt.fnu1 .and. fnuL.le.fnu2) then
*               csLyS(lr) = csLyS(lr) + cs2/(fnu2-fnu1)
*            endif
*         enddo
*
*      enddo
*
*      fnu1 = 0.5*(Txi(lrH1tr-1)+Txi(lrH1tr-2))/T1eV
*      fnu2 = 0.5*(Txi(lrH1tr-1)+Txi(lrH1tr+0))/T1eV
*      cs2 = 1.098e-16*oscstr_other/CSSTAR
*      csLyS(lrH1tr-1) = csLyS(lrH1tr-1) + cs2/(fnu2-fnu1)
C
      do lr=1,lrH1tr
         csLySn(lr) = csLyS(lr)/csH1tr
      enddo

      lrLStr = 0
      do lr=1,NRAD
         if(lrLStr.eq.0 .and. csLyS(lr).gt.1.0e-35) lrLStr = lr
      enddo
      if(lrLStr .eq. 0) lrLStr = NRAD

      lrmin = min(lrH1tr,lrDUtr,lrLStr)
C
C  Main tables
C
#if defined(RT_TABLES) && !defined(RT_DUST)

      tauH1m = 1.0
      tauG1m = 1.0
      tauG2m = 1.0

      altH1m = log(tauH1m)
      altG1m = log(tauG1m)
      altG2m = log(tauG2m)

      tH1min = tauH1m*acOmin
      tG1min = tauG1m*acOmin
      tG2min = tauG2m*acOmin

      tH1max = tauH1m*acOmax
      tG1max = tauG1m*acOmax
      tG2max = tauG2m*acOmax
C
C  OT part
C
      do j=1,iptabMax
#ifdef RT_EXTERNAL_BACKGROUND
         pOT(j) = 0.0
#endif
#ifdef RT_TRANSFER
         sgpOT(j) = 0.0
#endif
      enddo

#ifdef RT_CHEMISTRY
#ifdef RT_EXTERNAL_BACKGROUND
      call frtIntegrateRates_OT(iptabDim,wpTab,
     .     iptabOT,iptabMax,pOT)
#endif
#ifdef RT_TRANSFER
      call frtLocalIntegrateRates_OT(iptabDim,wpTab,
     .     iptabOT,iptabMax,sgpOT)
#endif
#endif  // RT_CHEMISTRY

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lo1,lo2,lo3,tauH1l,tauG1l,tauG2l,itind,itind2,j,
C$OMP+ pH1,pG1,pG2,sgpH1,sgpG1,sgpG2)
C$OMP+SHARED(tauH1m,acTAU,csH1tr,csG1tr,csG2tr,tauG1m,tauG2m,
C$OMP+ itloc2,pTab,sgpTab,pOT,sgpOT,
C$OMP+ wpTab,facPhot,facLX,iptabMax)
      do lo1=1,NOPT1

         do j=1,iptabMax
#ifdef RT_EXTERNAL_BACKGROUND
            pH1(j) = pOT(j)
#endif
#ifdef RT_TRANSFER
            sgpH1(j) = sgpOT(j)
#endif
         enddo

         if(lo1 .lt. NOPT1) then

            tauH1l = tauH1m*acTAU(lo1)
#ifdef RT_EXTERNAL_BACKGROUND
            call frtIntegrateRates_H1(iptabDim,wpTab,
     .           iptabH1,iptabMax,pH1,
     .           0.0,tauH1l)
#endif
#ifdef RT_TRANSFER
            call frtLocalIntegrateRates_H1(iptabDim,wpTab,
     .           iptabH1,iptabMax,sgpH1,
     .           0.0,tauH1l)
#endif

         else

            tauH1l = -1.0

         endif

         do lo2=1,NOPT1

            do j=1,iptabMax
#ifdef RT_EXTERNAL_BACKGROUND
               pG1(j) = pH1(j)
#endif
#ifdef RT_TRANSFER
               sgpG1(j) = sgpH1(j)
#endif
            enddo

            if(lo2 .lt. NOPT1) then

               tauG1l = tauG1m*acTAU(lo2)
#ifdef RT_EXTERNAL_BACKGROUND
               call frtIntegrateRates_G1(iptabDim,wpTab,
     .              iptabG1,iptabMax,pG1,
     .              0.0,tauH1l,tauG1l)
#endif
#ifdef RT_TRANSFER
               call frtLocalIntegrateRates_G1(iptabDim,wpTab,
     .              iptabG1,iptabMax,sgpG1,
     .              0.0,tauH1l,tauG1l)
#endif

            else

               tauG1l = -1.0

            endif

            itind2 = itloc2(lo2,lo1)

            do lo3=1,NOPT1

               do j=1,iptabMax
#ifdef RT_EXTERNAL_BACKGROUND
                  pG2(j) = pG1(j)
#endif
#ifdef RT_TRANSFER
                  sgpG2(j) = sgpG1(j)
#endif
               enddo

               if(lo3 .lt. NOPT1) then

                  tauG2l = tauG2m*acTAU(lo3)
#ifdef RT_EXTERNAL_BACKGROUND
                  call frtIntegrateRates_G2(iptabDim,wpTab,
     .                 iptabG2,iptabMax,pG2,
     .                 0.0,tauH1l,tauG1l,tauG2l)
#endif
#ifdef RT_TRANSFER
                  call frtLocalIntegrateRates_G2(iptabDim,wpTab,
     .                 iptabG2,iptabMax,sgpG2,
     .                 0.0,tauH1l,tauG1l,tauG2l)
#endif

               endif

               itind = lo3 + itind2
               do j=1,iptabMax
#ifdef RT_INTERPOLLOG
#ifdef RT_EXTERNAL_BACKGROUND
                  if(pG2(j) .gt. 0.0) then
                     pTab(j,itind) = log(1.0e-35+facPhot*pG2(j))
                  else
                     pTab(j,itind) = -1000.0
                  endif
#endif
#ifdef RT_TRANSFER
                  if(sgpG2(j) .gt. 0.0) then
                     sgpTab(j,itind) = log(1.0e-35+
     .                    facLX*facPhot*sgpG2(j))
                  else
                     sgpTab(j,itind) = -1000.0
                  endif
#endif
#else  // RT_INTERPOLLOG
#ifdef RT_EXTERNAL_BACKGROUND
                  pTab(j,itind) = max(0.0,facPhot*pG2(j))
#endif
#ifdef RT_TRANSFER
                  sgpTab(j,itind) = max(0.0,facLX*facPhot*sgpG2(j))
#endif
#endif // RT_INTERPOLLOG
               enddo

            enddo
         enddo
      enddo

#endif  // defined(RT_TABLES) && !defined(RT_DUST)

#ifdef RT_XRAYS
      call frtFillRadiationTablesSE
#endif

      return
      end
C
C
C     
      subroutine frtPhotoRatesTable(par,rf,itab,y,pRate,jmax1)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_TRANSFER
      include 'frt_transfer.inc'
#endif
      dimension par(*), rf(*)
      dimension itab(*)
      dimension y(*), pRate(*)
C
      dimension sgpRate(iptabDim)
#ifdef RT_XRAYS
      include 'frt_xrays.inc'
      dimension pRateSE(istabDim), sgpRateSE(istabDim)
#endif
C
C  This is just to placate ftnchek
C
      jmax = jmax1
C
C  Global field
C
#ifdef RT_EXTERNAL_BACKGROUND

#if defined(RT_TRANSFER)
      tauH1l = rf(irfH1gf)
      tauG1l = rf(irfG1gf)
      tauG2l = rf(irfG2gf)
#else
      tauH1l = 0.0
      tauG1l = 0.0
      tauG2l = 0.0
#endif

#if defined(RT_TABLES) && !defined(RT_DUST)

      call frtSampleTable(tauH1l,tauG1l,tauG2l,
     .     iptabDim,pTab,jmax,pRate,iOk)

      if(iOk .eq. 1) then
         itab(1) = itab(1) + 1
      else

#endif  // defined(RT_TABLES) && !defined(RT_DUST)

         do j=1,jmax
            pRate(j) = 0.0
         enddo

#if defined(RT_NARROWTABLE) || !defined(RT_TABLES) || defined(RT_DUST)

         itab(2) = itab(2) + 1

#ifdef RT_CHEMISTRY
         call frtIntegrateRates_OT(iptabDim,wpTab,
     .        iptabOT,jmax,pRate)
#endif  // RT_CHEMISTRY

         call frtIntegrateRates_H1(iptabDim,wpTab,
     .        iptabH1,jmax,pRate,
     .        par(iparDGas),tauH1l)

         call frtIntegrateRates_G1(iptabDim,wpTab,
     .        iptabG1,jmax,pRate,
     .        par(iparDGas),tauH1l,tauG1l)

         call frtIntegrateRates_G2(iptabDim,wpTab,
     .        iptabG2,jmax,pRate,
     .        par(iparDGas),tauH1l,tauG1l,tauG2l)

         do j=1,jmax
            pRate(j) = facPhot*max(0.0,pRate(j))
         enddo

#endif  // defined(RT_NARROWTABLE) || !defined(RT_TABLES) || defined(RT_DUST)

#if defined(RT_TABLES) && !defined(RT_DUST)
      endif
#endif  // defined(RT_TABLES) && !defined(RT_DUST)

#ifdef RT_CHEMISTRY
      pRate(iptabCiLW) = pRate0(iptabCiLW)
#endif  // RT_CHEMISTRY

#else   // RT_EXTERNAL_BACKGROUND

      do j=1,iptabDim
         pRate(j) = 0.0
      enddo

#endif  // RT_EXTERNAL_BACKGROUND
C
C  Local field
C
#ifdef RT_TRANSFER
      tauH1l = rf(irfH1lf)
      tauG1l = rf(irfG1lf)
      tauG2l = rf(irfG2lf)

#if defined(RT_TABLES) && !defined(RT_DUST)

      call frtSampleTable(tauH1l,tauG1l,tauG2l,
     .     iptabDim,sgpTab,jmax,sgpRate,iOk)

      if(iOk .eq. 1) then
         itab(1) = itab(1) + 1
      else

#endif  // defined(RT_TABLES) && !defined(RT_DUST)

         do j=1,jmax
            sgpRate(j) = 0.0
         enddo

#if defined(RT_NARROWTABLE) || !defined(RT_TABLES) || defined(RT_DUST)

         itab(2) = itab(2) + 1

#ifdef RT_CHEMISTRY
         call frtLocalIntegrateRates_OT(iptabDim,wpTab,
     .        iptabOT,jmax,sgpRate)
#endif  // RT_CHEMISTRY

         call frtLocalIntegrateRates_H1(iptabDim,wpTab,
     .        iptabH1,jmax,sgpRate,
     .        par(iparDGas),tauH1l)

         call frtLocalIntegrateRates_G1(iptabDim,wpTab,
     .        iptabG1,jmax,sgpRate,
     .        par(iparDGas),tauH1l,tauG1l)

         call frtLocalIntegrateRates_G2(iptabDim,wpTab,
     .        iptabG2,jmax,sgpRate,
     .        par(iparDGas),tauH1l,tauG1l,tauG2l)

         do j=1,jmax
            sgpRate(j) = facLX*facPhot*max(0.0,sgpRate(j))
         enddo

#endif  // defined(RT_NARROWTABLE) || !defined(RT_TABLES) || defined(RT_DUST)

#if defined(RT_TABLES) && !defined(RT_DUST)
      endif
#endif  // defined(RT_TABLES) && !defined(RT_DUST)

#ifdef RT_CHEMISTRY
      sgpRate(iptabCiLW) = sgpRate0(iptabCiLW)
#endif  // RT_CHEMISTRY

#endif  // RT_TRANSFER


#ifdef RT_XRAYS

      call frtGetPhotoRatesSE(par,rf,itab,y,pRateSE,sgpRateSE)

#ifdef RT_EXTERNAL_BACKGROUND
      pRate(iptabPiH1) = pRate(iptabPiH1) + 
     .     y(ieqXH1)*pRateSE(istabPiH1H1w1) +
     .     y(ieqXG1)*pRateSE(istabPiH1G1w1) +
     .     y(ieqXG2)*pRateSE(istabPiH1G2w1)
      pRate(iptabPiG1) = pRate(iptabPiG1) +
     .     y(ieqXH1)*pRateSE(istabPiG1H1w1) +
     .     y(ieqXG2)*pRateSE(istabPiG1G1w1) +
     .     y(ieqXG2)*pRateSE(istabPiG1G2w1)
      pRate(iptabPhH1) = max(0.0,
     .     pRate(iptabPhH1)-pRateSE(istabPhH1w1))
      pRate(iptabPhG1) = max(0.0,
     .     pRate(iptabPhG1)-pRateSE(istabPhG1w1))
#endif

#ifdef RT_TRANSFER
      sgpRate(iptabPiH1) = sgpRate(iptabPiH1) + 
     .     y(ieqXH1)*sgpRateSE(istabPiH1H1w1) +
     .     y(ieqXG1)*sgpRateSE(istabPiH1G1w1) +
     .     y(ieqXG2)*sgpRateSE(istabPiH1G2w1)
      sgpRate(iptabPiG1) = pRate(iptabPiG1) +
     .     y(ieqXH1)*sgpRateSE(istabPiG1H1w1) +
     .     y(ieqXG2)*sgpRateSE(istabPiG1G1w1) +
     .     y(ieqXG2)*sgpRateSE(istabPiG1G2w1)
      sgpRate(iptabPhH1) = max(0.0,sgpRate(iptabPhH1) -
     .     sgpRateSE(istabPhH1w1))
      sgpRate(iptabPhG1) = max(0.0,sgpRate(iptabPhG1) -
     .     sgpRateSE(istabPhG1w1))
#endif  // RT_TRANSFER

#endif  // RT_XRAYS

      par(iparCFH1) = 1.0
      par(iparCFG1) = 1.0
      par(iparCFG2) = 1.0

#if (RT_CFI == 1)
#ifdef RT_EXTERNAL_BACKGROUND
      ch = CHeat(par(iparRhoB),y)
      pRate(iptabPhH1) = pRate(iptabPhH1)*ch*par(iparBias)
      pRate(iptabPhG1) = pRate(iptabPhG1)*ch*par(iparBias)
      pRate(iptabPhG2) = pRate(iptabPhG2)*ch*par(iparBias)
      par(iparCFH1) = CIonH1(par(iparRhoB),y)
      par(iparCFG1) = CIonG1(par(iparRhoB),y)
      par(iparCFG2) = CIonG2(par(iparRhoB),y)
      pRate(iptabPiH1) = pRate(iptabPiH1)*par(iparCFH1)*par(iparBias)
      pRate(iptabPiG1) = pRate(iptabPiG1)*par(iparCFG1)*par(iparBias)
      pRate(iptabPiG2) = pRate(iptabPiG2)*par(iparCFG2)*par(iparBias)
#endif
#endif

#ifdef RT_TRANSFER
      do j=1,jmax
         pRate(j) = pRate(j) + rf(irfNorm)*sgpRate(j)
      enddo
#ifdef RT_CHEMISTRY
      pRate(iptabCiLW) = pRate(iptabCiLW) + rf(irfNorm)*
     .     sgpRate(iptabCiLW)
#endif
#endif  // RT_TRANSFER

#if (RT_CFI == 2)
      ch = CHeat(par(iparRhoB),y)
      pRate(iptabPhH1) = pRate(iptabPhH1)*ch*par(iparBias)
      pRate(iptabPhG1) = pRate(iptabPhG1)*ch*par(iparBias)
      pRate(iptabPhG2) = pRate(iptabPhG2)*ch*par(iparBias)
      par(iparCFH1) = CIonH1(par(iparRhoB),y)
      par(iparCFG1) = CIonG1(par(iparRhoB),y)
      par(iparCFG2) = CIonG2(par(iparRhoB),y)
      pRate(iptabPiH1) = pRate(iptabPiH1)*par(iparCFH1)*par(iparBias)
      pRate(iptabPiG1) = pRate(iptabPiG1)*par(iparCFG1)*par(iparBias)
      pRate(iptabPiG2) = pRate(iptabPig2)*par(iparCFG2)*par(iparBias)
#endif

      return
      end
C
C  ***************************************
C
C  Tools for manipulating tables
C
C  ***************************************
C
      subroutine frtIntegrateRates_OT(nmax,wi,j1,jn,s)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)

      do lr=1,lrmin-1
         fac = angEff(lr)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      do lr=lrmax+1,NRAD
         fac = angEff(lr)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return 
      end
C
C
C
      subroutine frtIntegrateRates_H1(nmax,wi,j1,jn,s,
     .     zdust,tauH1l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)

      if(tauH1l .lt. 0.0) return

      do lr=lrmin,lrG1tr-1
#ifdef RT_DUST
         csH1l = zdust*csDustn(lr) + csLySn(lr) + csH1n(lr)
#else
         csH1l = csLySn(lr) + csH1n(lr)
#endif
         taul = csH1l*tauH1l

         fac = angEff(lr)*QFuni(taul)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C
      subroutine frtIntegrateRates_G1(nmax,wi,j1,jn,s,
     .     zdust,tauH1l,tauG1l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)

      if(tauG1l .lt. 0.0) return

      do lr=lrG1tr,lrG2tr-1
         if(tauH1l .gt. 0.0) then
#ifdef RT_DUST
            csH1l = zdust*csDustn(lr) +  csH1n(lr)
#else
            csH1l = csH1n(lr)
#endif
            taul = csH1l*tauH1l
         else
            taul = 0.0
         endif
         taul = taul + csG1n(lr)*tauG1l

         fac = angEff(lr)*QFuni(taul)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C
      subroutine frtIntegrateRates_G2(nmax,wi,j1,jn,s,
     .     zdust,tauH1l,tauG1l,tauG2l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)

      if(tauG2l .lt. 0.0) return

      do lr=lrG2tr,lrmax
         if(tauH1l .gt. 0.0) then
#ifdef RT_DUST
            csH1l = zdust*csDustn(lr) +  csH1n(lr)
#else
            csH1l = csH1n(lr)
#endif
            taul = csH1l*tauH1l
         else
            taul = 0.0
         endif
         if(tauG1l .gt. 0.0) then
            taul = taul + csG1n(lr)*tauG1l
         endif
         taul = taul + csG2n(lr)*tauG2l

         fac = angEff(lr)*QFuni(taul)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C
      subroutine frtLocalIntegrateRates_OT(nmax,wi,j1,jn,s)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)

      do lr=1,lrmin-1
         fac = guvEff(lr)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      do lr=lrmax+1,NRAD
         fac = guvEff(lr)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return 
      end
C
C
C
      subroutine frtLocalIntegrateRates_H1(nmax,wi,j1,jn,s,
     .     zdust,tauH1l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)

      if(tauH1l .lt. 0.0) return

      do lr=lrmin,lrG1tr-1
#ifdef RT_DUST
         csH1l = zdust*csDustn(lr) + csLySn(lr) + csH1n(lr)
#else
         csH1l = csLySn(lr) + csH1n(lr)
#endif
         taul = csH1l*tauH1l

         fac = guvEff(lr)*QFloc(taul)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C
      subroutine frtLocalIntegrateRates_G1(nmax,wi,j1,jn,s,
     .     zdust,tauH1l,tauG1l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)

      if(tauG1l .lt. 0.0) return

      do lr=lrG1tr,lrG2tr-1
         if(tauH1l .gt. 0.0) then
#ifdef RT_DUST
            csH1l = zdust*csDustn(lr) +  csH1n(lr)
#else
            csH1l = csH1n(lr)
#endif
            taul = csH1l*tauH1l
         else
            taul = 0.0
         endif
         taul = taul + csG1n(lr)*tauG1l

         fac = guvEff(lr)*QFloc(taul)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C
      subroutine frtLocalIntegrateRates_G2(nmax,wi,j1,jn,s,
     .     zdust,tauH1l,tauG1l,tauG2l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)

      if(tauG2l .lt. 0.0) return

      do lr=lrG2tr,lrmax
         if(tauH1l .gt. 0.0) then
#ifdef RT_DUST
            csH1l = zdust*csDustn(lr) +  csH1n(lr)
#else
            csH1l = csH1n(lr)
#endif
            taul = csH1l*tauH1l
         else
            taul = 0.0
         endif
         if(tauG1l .gt. 0.0) then
            taul = taul + csG1n(lr)*tauG1l
         endif
         taul = taul + csG2n(lr)*tauG2l

         fac = guvEff(lr)*QFloc(taul)
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C     
      subroutine frtSampleTable(tauH1l,tauG1l,tauG2l,
     .     nmax,table,nout,out,iOk)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension table(nmax,*), out(*)

      if(tauH1l .lt. 0.0) then
         lH1tab = 2
      else if(tauH1l .lt. tH1max) then
         lH1tab = 1
      else
         lH1tab = 0
      endif

      if(tauG1l .lt. 0.0) then
         lG1tab = 2
      else if(tauG1l .lt. tG1max) then
         lG1tab = 1
      else
         lG1tab = 0
      endif

      if(tauG2l .lt. 0.0) then
         lG2tab = 2
      else if(tauG2l .lt. tG2max) then
         lG2tab = 1
      else
         lG2tab = 0
      endif

      if(lH1tab.gt.0 .and. lG1tab.gt.0 .and. lG2tab.gt.0) then

         if(lH1tab .eq. 2) then
            j1l = NOPT
            j1u = NOPT1
            p1l = 0.0
            p1u = 1.0
         else if(tauH1l .gt. tH1min) then
            altau = log(tauH1l)-altH1m
            ltaul = int((altau-acLOmin)/acTOstp*0.99999)+1
            ltauu = ltaul + 1
            alopt = acT2O(ltaul) + (acT2O(ltauu)-acT2O(ltaul))*
     .           (altau-acLOmin-acTOstp*real(ltaul-1))/acTOstp
            j1l = int((alopt-acLOmin)/acLOstp*0.99999)+2
            if(j1l .lt. 2   ) j1l = 2
            if(j1l .ge. NOPT) j1l = NOPT-1
            j1u = j1l + 1
            p1l = max(0.0,min(1.0,(acLO(j1u)-alopt)/acLOstp))
            p1u = 1.0-p1l
         else
            j1l = 1
            j1u = 2
            p1u = max(0.0,tauH1l/tH1min)
            p1l = 1.0-p1u
         endif

         if(lG1tab .eq. 2) then
            j2l = NOPT
            j2u = NOPT1
            p2l = 0.0
            p2u = 1.0
         else if(tauG1l .gt. tG1min) then
            altau = log(tauG1l)-altG1m
            ltaul = int((altau-acLOmin)/acTOstp*0.99999)+1
            ltauu = ltaul + 1
            alopt = acT2O(ltaul) + (acT2O(ltauu)-acT2O(ltaul))*
     .           (altau-acLOmin-acTOstp*real(ltaul-1))/acTOstp
            j2l = int((alopt-acLOmin)/acLOstp*0.99999)+2
            if(j2l .lt. 2   ) j2l = 2
            if(j2l .ge. NOPT) j2l = NOPT-1
            j2u = j2l + 1
            p2l = max(0.0,min(1.0,(acLO(j2u)-alopt)/acLOstp))
            p2u = 1.0-p2l
         else
            j2l = 1
            j2u = 2
            p2u = max(0.0,tauG1l/tG1min)
            p2l = 1.0-p2u
         endif

         if(lG2tab .eq. 2) then
            j3l = NOPT
            j3u = NOPT1
            p3l = 0.0
            p3u = 1.0
         else if(tauG2l .gt. tG2min) then
            altau = log(tauG2l)-altG2m
            ltaul = int((altau-acLOmin)/acTOstp*0.99999)+1
            ltauu = ltaul + 1
            alopt = acT2O(ltaul) + (acT2O(ltauu)-acT2O(ltaul))*
     .           (altau-acLOmin-acTOstp*real(ltaul-1))/acTOstp
            j3l = int((alopt-acLOmin)/acLOstp*0.99999)+2
            if(j3l .lt. 2   ) j3l = 2
            if(j3l .ge. NOPT) j3l = NOPT-1
            j3u = j3l + 1
            p3l = max(0.0,min(1.0,(acLO(j3u)-alopt)/acLOstp))
            p3u = 1.0-p3l
         else
            j3l = 1
            j3u = 2
            p3u = max(0.0,tauG2l/tG2min)
            p3l = 1.0-p3u
         endif

         itlll = j3l + itloc2(j2l,j1l)
         itull = j3u + itloc2(j2l,j1l)
         itlul = j3l + itloc2(j2u,j1l)
         ituul = j3u + itloc2(j2u,j1l)
         itllu = j3l + itloc2(j2l,j1u)
         itulu = j3u + itloc2(j2l,j1u)
         itluu = j3l + itloc2(j2u,j1u)
         ituuu = j3u + itloc2(j2u,j1u)
         
         wlll = p3l*p2l*p1l
         wull = p3u*p2l*p1l
         wlul = p3l*p2u*p1l
         wuul = p3u*p2u*p1l
         wllu = p3l*p2l*p1u
         wulu = p3u*p2l*p1u
         wluu = p3l*p2u*p1u
         wuuu = p3u*p2u*p1u

         do j=1,nout
            out(j) = 
     .           wlll*table(j,itlll) + wull*table(j,itull) +
     .           wlul*table(j,itlul) + wuul*table(j,ituul) +
     .           wllu*table(j,itllu) + wulu*table(j,itulu) +
     .           wluu*table(j,itluu) + wuuu*table(j,ituuu)
#ifdef RT_INTERPOLLOG
            out(j) = exp(out(j))
#endif
         enddo

         iOk = 1

      else
         
         iOk = 0
         
      endif

      return
      end
C
C
C
      block data frt_LYSS_DATA_SET
      common/RT_Z_LYSS/ oscstr(48), oscstr_other
      data (oscstr(i), i=1,48)
     .     /4.161e-1, 7.910e-2, 2.899e-2, 1.393e-2, 7.799e-3, 4.814e-3,
     .      3.183e-3, 2.216e-3, 1.605e-3, 1.201e-3, 9.214e-4, 7.227e-4,
     .      5.774e-4, 4.686e-4, 3.856e-4, 3.211e-4, 2.702e-4, 2.296e-4,
     .      1.967e-4, 1.698e-4, 1.476e-4, 1.291e-4, 1.136e-4, 1.005e-4,
     .      8.928e-5, 7.970e-5, 7.144e-5, 6.429e-5, 5.806e-5, 5.261e-5,
     .      4.782e-5, 4.360e-5, 3.986e-5, 3.653e-5, 3.357e-5, 3.092e-5,
     .      2.854e-5, 2.639e-5, 2.446e-5, 2.271e-5, 2.113e-5, 1.969e-5,
     .      1.837e-5, 1.717e-5, 1.608e-5, 1.507e-5, 1.415e-5, 1.330e-5/
      end
C
C
C
      block data frt_SD93COOLING
      common/RT_Z_SD93/ sd93lt(90), sd93c0(90), sd93c1(90)
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=1,10)
     .  /   4,   1.327e-26,   4.666e-27,
     .   4.05,   3.329e-25,   4.478e-26,
     .    4.1,   2.297e-24,   3.933e-25,
     .   4.15,   1.044e-23,   2.523e-24,
     .    4.2,   2.647e-23,   2.799e-24,
     .   4.25,   3.983e-23,   9.227e-24,
     .    4.3,   4.793e-23,   2.003e-23,
     .   4.35,   6.077e-23,   2.104e-23,
     .    4.4,     8.7e-23,   3.303e-23,
     .   4.45,   1.231e-22,   4.527e-23/
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=11,20)
     .  / 4.5,   1.799e-22,   6.018e-23,
     .   4.55,   2.431e-22,   8.119e-23,
     .    4.6,   3.285e-22,   1.131e-22,
     .   4.65,   4.342e-22,    1.39e-22,
     .    4.7,   5.531e-22,    1.89e-22,
     .   4.75,   6.975e-22,   2.174e-22,
     .    4.8,   8.523e-22,   2.549e-22,
     .   4.85,   1.014e-21,   2.832e-22,
     .    4.9,   1.161e-21,   3.184e-22,
     .   4.95,   1.281e-21,   3.669e-22/
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=21,30)
     .  /   5,   1.293e-21,   4.349e-22,
     .   5.05,   1.215e-21,   5.088e-22,
     .    5.1,   1.196e-21,   6.002e-22,
     .   5.15,   1.261e-21,   6.939e-22,
     .    5.2,   1.359e-21,   7.601e-22,
     .   5.25,   1.394e-21,    8.39e-22,
     .    5.3,    1.43e-21,   9.207e-22,
     .   5.35,   1.433e-21,   1.039e-21,
     .    5.4,   1.338e-21,   9.176e-22,
     .   5.45,   9.675e-22,   6.999e-22/
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=31,40)
     .  / 5.5,   5.799e-22,   3.919e-22,
     .   5.55,   3.634e-22,   2.418e-22,
     .    5.6,    2.81e-22,   1.853e-22,
     .   5.65,   2.561e-22,   1.687e-22,
     .    5.7,   2.504e-22,   1.652e-22,
     .   5.75,   2.507e-22,   1.559e-22,
     .    5.8,   2.231e-22,   1.294e-22,
     .   5.85,   1.761e-22,   8.714e-23,
     .    5.9,   1.491e-22,   5.655e-23,
     .   5.95,   1.423e-22,    4.41e-23/
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=41,50)
     .  /   6,   1.389e-22,   3.384e-23,
     .   6.05,   1.389e-22,   2.527e-23,
     .    6.1,   1.423e-22,    2.19e-23,
     .   6.15,   1.457e-22,   1.845e-23,
     .    6.2,   1.457e-22,   2.269e-23,
     .   6.25,   1.355e-22,   2.008e-23,
     .    6.3,   1.063e-22,   1.961e-23,
     .   6.35,   7.729e-23,   1.891e-23,
     .    6.4,   5.868e-23,   1.603e-23,
     .   6.45,   4.754e-23,   1.425e-23/
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=51,60)
     .  / 6.5,   4.155e-23,   1.094e-23,
     .   6.55,   3.821e-23,   1.055e-23,
     .    6.6,   3.508e-23,   8.592e-24,
     .   6.65,   3.199e-23,   7.796e-24,
     .    6.7,   2.923e-23,   7.269e-24,
     .   6.75,   2.726e-23,   5.531e-24,
     .    6.8,   2.535e-23,   5.056e-24,
     .   6.85,   2.517e-23,   5.945e-24,
     .    6.9,   2.636e-23,   4.885e-24,
     .   6.95,   2.677e-23,    3.23e-24/
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=61,70)
     .  /   7,   2.657e-23,   4.131e-24,
     .   7.05,   2.534e-23,   4.031e-24,
     .    7.1,   2.259e-23,   2.744e-24,
     .   7.15,   1.928e-23,   3.316e-24,
     .    7.2,   1.683e-23,   2.992e-24,
     .   7.25,   1.508e-23,   3.139e-24,
     .    7.3,   1.337e-23,    3.18e-24,
     .   7.35,   1.224e-23,   3.805e-24,
     .    7.4,   1.166e-23,   3.832e-24,
     .   7.45,   1.074e-23,   3.149e-24/
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=71,80)
     .  / 7.5,   1.066e-23,   2.511e-24,
     .   7.55,   9.983e-24,    2.41e-24,
     .    7.6,   9.855e-24,   3.122e-24,
     .   7.65,   1.032e-23,   3.269e-24,
     .    7.7,   1.016e-23,   4.065e-24,
     .   7.75,   1.064e-23,   2.608e-24,
     .    7.8,   1.001e-23,   2.458e-24,
     .   7.85,   1.049e-23,   2.574e-24,
     .    7.9,   1.098e-23,   2.695e-24,
     .   7.95,   1.073e-23,   3.593e-24/
      data (sd93lt(i),sd93c0(i),sd93c1(i),i=81,90)
     .  /   8,   1.069e-23,   2.582e-24,
     .   8.05,   1.119e-23,   2.704e-24,
     .    8.1,   1.172e-23,   5.107e-24,
     .   8.15,   1.227e-23,   5.348e-24,
     .    8.2,   1.217e-23,   4.148e-24,
     .   8.25,   1.275e-23,   4.343e-24,
     .    8.3,   1.335e-23,   4.548e-24,
     .   8.35,   1.319e-23,   3.055e-24,
     .    8.4,   1.381e-23,     3.2e-24,
     .   8.45,   1.446e-23,   6.566e-24/
      end
#endif
      subroutine frtDummyTables(i)
      i = 0
      end
