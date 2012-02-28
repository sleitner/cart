#include "../config.h"
#ifdef RADIATIVE_TRANSFER
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
C  Effective dust cross-section for the MW in the LW band
C  (for H2 formation)
C
      d_sDUST = 2.0e-21/CSSTAR
C
C  Radiation tables
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,fnu,x,cs0,fn0,y0,y1,yw,ya,p,y,w)
C$OMP+SHARED(Txi,xig,angUni,angLoc,angAvg,tauLoc,csH1,csG1,csG2,csDust,
C$OMP+ csMH27,csMH28,csMH29,csMH30,csMH31,csMH32,csMHLW,csUV,
C$OMP+ csO8,csS16,csFe26)
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

         if(fnu .lt. 871.4) then
            csO8(lr) = 0.0
         else
            cs0 = 8.554e-16
            fn0 = 27.54
            y0  = 0.0 
            y1  = 0.0
            yw  = 0.0
            ya  = 3.288e+01
            p   = 2.963
            x = fnu/fn0 - y0
            y = sqrt(x**2+y1**2)
            csO8(lr) = cs0*((x-1.0)**2+yw**2)*y**(0.5*p-5.5)/
     .           (1.0+sqrt(y/ya))**p/CSSTAR
         endif

         if(fnu .lt. 3494) then
            csS16(lr) = 0.0
         else
            cs0 = 2.139e-16
            fn0 = 110.4
            y0  = 0.0 
            y1  = 0.0
            yw  = 0.0
            ya  = 3.288e+01
            p   = 2.963
            x = fnu/fn0 - y0
            y = sqrt(x**2+y1**2)
            csS16(lr) = cs0*((x-1.0)**2+yw**2)*y**(0.5*p-5.5)/
     .           (1.0+sqrt(y/ya))**p/CSSTAR
         endif

         if(fnu .lt. 9278) then
            csFe26(lr) = 0.0
         else
            cs0 = 8.099e-17
            fn0 = 293.2
            y0  = 0.0 
            y1  = 0.0
            yw  = 0.0
            ya  = 3.288e+01
            p   = 2.963
            x = fnu/fn0 - y0
            y = sqrt(x**2+y1**2)
            csFe26(lr) = cs0*((x-1.0)**2+yw**2)*y**(0.5*p-5.5)/
     .           (1.0+sqrt(y/ya))**p/CSSTAR
         endif

         w = 1.24/fnu

#if (RT_DUST_CS == 1)
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
            csDust(lr) = csDust(lr) + 5.6e-22/CSSTAR*(
     .           dustfit(w, 0.067 ,  10.0 ,  1.90, 4.0, 15.0) )
         endif
#else
C
C  LMC dust
C  My fudging: 175->90 in term #1, 5.5->21 in term #2, plus my term #7
C  The same term #7 for both SMC and LMC
C  Assuming ZLMC=-0.5dex from Welty et al 1999, ApJ, 512, 636 
C
         csDust(lr) = 1.1e-21/CSSTAR*(
     .        dustfit(w, 0.046, 90.0  , 90.0 , 2.0, 2.0) +
     .        dustfit(w, 0.08 ,  19.0  , 21.0, 4.5, 4.5) +  !! 21, not 5.5
     .        dustfit(w, 0.22 ,   0.023, -1.95, 2.0, 2.0) +
     .        dustfit(w, 9.7  ,   0.005, -1.95, 2.0, 2.0) +
     .        dustfit(w,18.0  ,   0.006, -1.80, 2.0, 2.0) +
     .        dustfit(w,25.0  ,   0.020,  0.0 , 2.0, 2.0) )
         if(w .gt. 1.0e-3) then
            csDust(lr) = csDust(lr) + 1.1e-21/CSSTAR*(
     .           dustfit(w, 0.067 ,  10.0 ,  1.90, 4.0, 15.0) )
         endif
#endif

         csUV(lr) = csDust(lr)

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

#endif
C
C  We need the LW rate even if CHEMISTRY = off
C
         if(fnu.lt.11.3 .or. fnu.gt.13.6) then
            csMHLW(lr) = 0.0
         else
C
C  R = 10^8 s^{-1} I_\nu  (that is what Tom Abel told me)
C
            csMHLW(lr) = 3.6e-18/CSSTAR
         endif

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
         if(lrH1tr.eq.0 .and. csH1(lr).gt.VSMALL) lrH1tr = lr
         if(lrG1tr.eq.0 .and. csG1(lr).gt.VSMALL) lrG1tr = lr
         if(lrG2tr.eq.0 .and. csG2(lr).gt.VSMALL) lrG2tr = lr
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
      lrDUtr = 0
      do lr=1,NRAD
         if(lrDUtr.eq.0 .and. csDust(lr).gt.VSMALL) lrDUtr = lr
      enddo
      if(lrDUtr .eq. 0) lrDUtr = NRAD
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
C$OMP+SHARED(csH1,csG1,csG2,csH1n,csG1n,csG2n,csH1tr,csG1tr,csG2tr)
      do lr=1,NRAD
         csH1n(lr) = csH1(lr)/csH1tr
         csG1n(lr) = csG1(lr)/csG1tr
         csG2n(lr) = csG2(lr)/csG2tr
      enddo
C
      lrUVtr = 0
      do lr=1,NRAD
         if(lrUVtr.eq.0 .and. Txi(lr).gt.12.0/13.6*TH1) lrUVtr = lr
      enddo
      if(lrUVtr .eq. 0) lrUVtr = NRAD

      csUVtr = csUV(lrUVtr)

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(csUV,csUVn,csUVtr)
      do lr=1,NRAD
         if(csUVtr .gt. 1.0e-30) then
            csUVn(lr) = csUV(lr)/csUVtr
         else
            csUVn(lr) = 0.0
         endif
      enddo
C
      lrLWtr = 0
      do lr=1,NRAD
         if(lrLWtr.eq.0 .and. csMHLW(lr).gt.VSMALL) lrLWtr = lr
      enddo
      if(lrLWtr .eq. 0) lrLWtr = 1

#ifdef RT_CHEMISTRY

      lr27tr = 0
      lr28tr = 0
      lr29tr = 0
      lr30tr = 0
      lr31tr = 0
      lr32tr = 0
      do lr=1,NRAD
         if(lr27tr.eq.0 .and. csMH27(lr).gt.VSMALL) lr27tr = lr
         if(lr28tr.eq.0 .and. csMH28(lr).gt.VSMALL) lr28tr = lr
         if(lr29tr.eq.0 .and. csMH29(lr).gt.VSMALL) lr29tr = lr
         if(lr30tr.eq.0 .and. csMH30(lr).gt.VSMALL) lr30tr = lr
         if(lr31tr.eq.0 .and. csMH31(lr).gt.VSMALL) lr31tr = lr
         if(lr32tr.eq.0 .and. csMH32(lr).gt.VSMALL) lr32tr = lr
      enddo
      if(lr27tr .eq. 0) lr27tr = 1
      if(lr28tr .eq. 0) lr28tr = 1
      if(lr29tr .eq. 0) lr29tr = 1
      if(lr30tr .eq. 0) lr30tr = 1
      if(lr31tr .eq. 0) lr31tr = 1
      if(lr32tr .eq. 0) lr32tr = 1

#ifdef RT_LWBANDS
C
C  Lyman-Werner band with its own frequency dependence
C
      XILWST = log(fnuLW2/fnuLW1)/real(NRLW-1)
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(xigLW,angLW,XILWST)
      do lr=1,NRLW
         xigLW(lr) = XILWST*(lr-NRLW)
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
C$OMP+PRIVATE(lc,tem,altem,q,p,alt10,srt,alH1,alG1,alG2,s,
C$OMP+ tev,x,itl,itu,so,sp,rC1a,rC2e,rC2a,rN1e,
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
         tev = tem/T1eV
C
C  Atomic rates
C
         alH1 = 2.0*TH1/tem
         alG1 = 2.0*TG1/tem
         alG2 = 2.0*TG2/tem
         q = 1.0/tem**1.5
         p = 1.0/(1.0+sqrt(tem/1.0e5))

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
         x = alt10 - 1.0
         if(x .lt. 3.4) then
            RateAtom(irateRAiG2,lc) = 1.0e-11/srt*(10.85-2.34*x+
     .           0.019*x**3)
            RateAtom(irateRAcG2,lc) = 1.0e-11/srt*( 9.92-2.37*x+
     .           0.030*x**3)*(1.381e-16/CFSTAR)*tem
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

         x = alt10 - 1.0
         if(x .lt. 3.4) then
            RateAtom(irateRBiG2,lc) = 1.0e-11/srt*(9.28-2.31*x+
     .           0.016*x**3)
            RateAtom(irateRBcG2,lc) = 1.0e-11/srt*(8.35-2.37*x+
     .           0.028*x**3)*(1.381e-16/CFSTAR)*tem
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
C
C   1. H + e -> Hm + g
C
#if (RT_H2_RATE == 1)
         if(tem .gt. 1.5e4) then
            s = 10.0**(-14.10+0.1175*alt10-9.813e-3*alt10**2)
         else
            s = 1.0e-18*tem
         endif
#elif (RT_H2_RATE == 2)
         s = 1.4e-18*tem**0.928*exp(-tem/16200)
#else
         if(tem .lt. 6.0e3) then
            s = 10.0**(-17.845+7.62e-1*alt10+1.52e-1*alt10**2-
     .           3.27e-2*alt10**3)
         else
            s = 10.0**(-16.420+2.00e-1*alt10-5.45e-3*alt10**4+
     .           4.04e-5*alt10**6)
         endif
#endif
         RateChem(irateMHi01,lc) = s
C
C   2. Hm + H -> H2 + e
C
#if (RT_H2_RATE == 1)
         if(tem .gt. 1.0e4) then
            s = 10.0**(-8.78+1.13e-1*alt10-3.48e-2*alt10**2)
         else
            s = 1.3e-9
         endif
#elif (RT_H2_RATE == 2)
         if(tem .gt. 300.) then
            s = 4.0e-9*tem**(-0.17)
         else
            s = 1.5e-9
         endif
#else
         s = 1.3e-9
#endif
         RateChem(irateMHi02,lc) = s
C
C   3. H + Hp -> H2p + g
C
#if (RT_H2_RATE == 1)
         if(tem .gt. 6.7e3) then
            s = 5.81e-16*(tem/5.62e4)**(3.16-0.67*alt10)
         else
            s = 1.85e-23*tem**1.8
         endif
#elif (RT_H2_RATE == 2)
         s = 10.0**(-19.38-1.52*alt10+1.12*alt10**2-0.127*alt10**3)
#else
         s = 10.0**(-19.38-1.52*alt10+1.12*alt10**2-0.127*alt10**3)
#endif
         RateChem(irateMHi03,lc) = s
C
C   4. H + H2p -> H2 + Hp
C
#if (RT_H2_RATE == 1)
         s = 6.0e-10
#elif (RT_H2_RATE == 2)
         s = 6.4e-10
#else
         s = 6.4e-10
#endif
         RateChem(irateMHi04,lc) = s
C
C   5. H + H2p -> H2 + Hp
C
#if (RT_H2_RATE == 1)
         s = 4.0e-6/srt
#elif (RT_H2_RATE == 2)
         s = 5.7e-6/tem**0.5 + 6.3e-8 - 9.2e-11*tem**0.5 + 4.4e-13*tem
#else
         s = 2.4e-6/srt*(1+tem/2.0e4)
#endif
         RateChem(irateMHi05,lc) = s
C
C   6. H + H2p -> H2 + Hp
C
#if (RT_H2_RATE == 1)
         s = 2.25e-6/tem**0.4
#elif (RT_H2_RATE == 2)
         s = 2.0e-7/tem**0.5
#else
         s = 1.32e-6/max(617.0,tem)**0.76
#endif
         RateChem(irateMHi06,lc) = s
C
C   7. H2 + Hp -> H2p + H
C
#if (RT_H2_RATE == 1)
         s = 2.4e-9*exp(-2.12e4/tem)
#elif (RT_H2_RATE == 2)
         if(tem .gt. 1.0e4) then
            s = 1.5e-10*exp(-1.40e4/tem)
         else
            s = 3.0e-10*exp(-2.11e4/tem)
         endif
#else
C
C This rate from G&A08 is a bad fitting form and causes some weird behavior
C  for very low T. Use G&P instead.
C
*         s = (-3.32e-7+3.37e-7*alt10-1.45e-7*alt10**2+3.42e-8*alt10**3-
*     .        4.78e-9*alt10**4+3.97e-10*alt10**5-1.82e-11*alt10**6+
*     .        3.53e-13*alt10**7)*exp(-2.12e4/tem)
         if(tem .gt. 1.0e4) then
            s = 1.5e-10*exp(-1.40e4/tem)
         else
            s = 3.0e-10*exp(-2.11e4/tem)
         endif
#endif
         RateChem(irateMHi07,lc) = s
C
C   8. H2 + e -> H + H + e
C
#if (RT_H2_RATE == 1)
C
C  This is Shapiro & Kang 98 rate
C
*         s = 4.38e-10*exp(-1.02e5/tem)*(tem/(1.0+tem/1.0e5))**0.35
C
C  This is Donahue & Shull 91 rate (a smaller one)
C
         s = 5.6e-11*exp(-1.02e5/tem)*(tem/(1.0+tem/1.0e5))**0.5
#elif (RT_H2_RATE == 2)
         s = 4.4e-10*tem**0.35*exp(-1.02e5/tem)
#else
         s = 4.49e-9*tem**0.11*exp(-1.02e5/tem)
#endif
         RateChem(irateMHi08,lc) = s
C
C   9. H2 + H -> H + H + H
C
#if (RT_H2_RATE == 1)
         s = 1.12e-10*exp(-7.035e4/tem)
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = 6.67e-12*tem**0.5*exp(-(1+6.36e4/tem))
#endif
         RateChem(irateMHi09,lc) = s
C
C  10. H2 + H2 -> H2 + H + H
C
#if (RT_H2_RATE == 1)
         s = 1.18e-10*exp(-6.95e4/tem)
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = 6.0e-30*tem**4.19/(1+6.76e-6*tem)**5.69*exp(-5.47e4/tem)
#endif
         RateChem(irateMHi10,lc) = s
C
C  11. H2 + He -> H + H + He
C
#if (RT_H2_RATE == 1)
         s = 0.0
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = 10.0**(-27.029+3.80*alt10-2.95e4/tem)
#endif
         RateChem(irateMHi11,lc) = s
C
C  14. Hm + e -> H + e + e
C
#if (RT_H2_RATE == 1)
         s = 4.0e-12*exp(-8.75e3/tem)*(2*tem/(1+tem/1.6e4))
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = exp(-18.018+2.36*altem-2.83e-1*altem**2+1.62e-2*altem**3-
     .        3.37e-2*altem**4+1.18e-2*altem**5-1.66e-3*altem**6+
     .        1.07e-4*altem**7-2.63e-6*altem**8)
#endif
         RateChem(irateMHi14,lc) = s
C
C  15. Hm + H -> H + H + e
C
#if (RT_H2_RATE == 1)
         s = 5.3e-20*exp(-8.75e3/tem)*(2*tem/(1+tem/1.6e4))**2.17
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         if(tev .lt. 0.1) then
            s = 2.56e-9*tem**1.78
         else
            s = exp(-20.372+1.14*altem-1.42e-1*altem**2+
     .           8.46e-3*altem**3-1.43e-3*altem**4+2.01e-4*altem**5+
     .           8.66e-5*altem**6-2.59e-5*altem**7+2.46e-6*altem**8-
     .           8.07e-8*altem**9)
         endif
#endif
         RateChem(irateMHi15,lc) = s
C
C  16. Hp + Hm -> H2p + e
C
#if (RT_H2_RATE == 1)
         if(tem .gt. 1.0e4) then
            s = 4.0e-4/tem**1.4*exp(-1.51e4/tem)
         else
            s = 1.0e-8/tem**0.4
         endif
#elif (RT_H2_RATE == 2)
         if(tem .gt. 8.0e3) then
            s = 9.6e-7/tem**0.9
         else
            s = 6.9e-9/tem**0.35
         endif
#else
         if(tem .gt. 8.0e3) then
            s = 9.6e-7/tem**0.9
         else
            s = 6.9e-9/tem**0.35
         endif
#endif
         RateChem(irateMHi16,lc) = s
C
C  21. Hm + H2p -> H2 + H
C
#if (RT_H2_RATE == 1)
         s = 4.0e-6/srt
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = 4.0e-6/srt
#endif
         RateChem(irateMHi21,lc) = s
C
C  22. Hm + H2p -> H + H + H
C
#if (RT_H2_RATE == 1)
         s = 0.0
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = 1.4e-7*(300/tem)**0.5
#endif
         RateChem(irateMHi22,lc) = s
C
C  23. H2 + e -> Hm + H
C
#if (RT_H2_RATE == 1)
         if(tem .gt. 2.5e3) then
            s = 9.69e-13*exp(-11.32/(altem-7.28))
         else
            s = 0.0
         endif
#elif (RT_H2_RATE == 2)
         s = 2.7e-8/tem**1.27*exp(-4.3e4/tem)      
#else
         s = 2.7e-8/tem**1.27*exp(-4.3e4/tem)      
#endif
         RateChem(irateMHi23,lc) = s
C
C  *****************
C
C  Reaction 24 and 25 destroy H2 at NH>1e23 way too efficiently.
C  Something must be wrong with them (probably rates go to zero
C  at low T).
C
C  *****************
C
C  24. H2 + Hep -> He + H + Hp
C
#if (RT_H2_RATE == 1)
         s = 0.0
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
C
C  G&A have a plus sign under the exponent, but that must be a typo, it causes 
C  truly weird results.
C
         s = 3.7e-14*exp(-35/tem)
#endif
         RateChem(irateMHi24,lc) = s
C
C  25. H2 + Hep -> H2p + He
C
#if (RT_H2_RATE == 1)
         s = 0.0
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = 7.2e-15
#endif
         RateChem(irateMHi25,lc) = s
C
C  26. Hep + H -> He + Hp
C
#if (RT_H2_RATE == 1)
         s = 0.0
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
C
C  I have no idea where this formula is coming from. G&A refer to 
C  Zygelman et al. 1989, PhRA, 40, 2340, but results there are not 
C  consistent with that fit. Their table is consistent with a flat rate
C  of 8e-15 until T=1e3. With no data above T=1e3, I chop the rate off
C  just in case.
C
*         s = 1.2e-15*(tem/300)**0.25
         s = 8.0e-15*exp(-tem/3.0e3)
#endif
         RateChem(irateMHi26,lc) = s
C
C  27. He + Hp -> Hep + H
C
#if (RT_H2_RATE == 1)
         s = 0.0
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         if(tem .lt. 1.0e4) then
            s = 1.26e-9/tem**0.75*exp(-1.28e5/tem)
         else
C
C This grows too fast, clump it at T=1e7
C
*            s = 4.0d-37*tem**4.74
            s = 6.05e-4*min(1.0,tem/1.0e7)**4.74
         endif
#endif
         RateChem(irateMHi27,lc) = s
C
C  28. Hep + Hm -> He + H
C
#if (RT_H2_RATE == 1)
         s = 0.0
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = 2.32e-7*(300/tem)**0.52*exp(-tem/2.24e4)
#endif
         RateChem(irateMHi28,lc) = s
C
C  29. He + Hm -> He + H + e
C
#if (RT_H2_RATE == 1)
         s = 0.0
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         s = 4.1e-17*tem**2*exp(-1.99e4/tem)
#endif
         RateChem(irateMHi29,lc) = s
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
C  Dust: from Draine, ApJ, 245, 880
C
         RateChem(irateMHcDU,lc) = (2.0e-22/CFSTAR)*(tem/1.0e6)**1.5/
     .        (1+(tem/1.36e7)**1.5)
C
C  H2 cooling
C
C  1. collisions with H
C
#if (RT_H2_RATE == 1)
         if(tem .gt. 1635.0) then
            so = (8.18e-25/CFSTAR)*srt*exp(-1.0e3/tem)
         else
            so = (1.15e-25/CFSTAR)*exp(tem/125-(tem/577)**2)
         endif

         if(tem .gt. 4031.0) then
            so = so + (1.66e-22/CFSTAR)*exp(-9.24e3/tem)
         else
            x = log10(tem/1.0e4)
            so = so + 1.2*10.0**(-22.9-0.553*x-1.148*x**2)/CFSTAR
         endif
         sp = so
#elif (RT_H2_RATE == 2)
         so = 10.0**(-103.0+97.59*alt10-48.05*alt10**2+10.80*alt10**3-
     .        0.9032*alt10**4-log10(CFSTAR))
         sp = so
#else
C
C  Factor q is my correction to reproduce G&P behaviour at T>1e4 - the
C  G&A rates vlow up unphysically at >1e4
C
         q = exp(-(tem/17000)**2)
         x = log10(min(10.0,tem/1.0e3))
         if(tem .lt. 100) then
            so = 5.09e-27/CFSTAR*sqrt(tem/1.0e3)*exp(-853/tem)
         else if(tem .lt. 1.0e3) then
            so = 10.0**(-24.331+4.440*x-4.046*x**2-1.139*x**3+
     .           9.809*x**4+8.627*x**5-log10(CFSTAR))
         else
            so = 10.0**(-24.329+4.611*x-3.950*x**2+12.363*x**3-
     .           32.403*x**4+48.853*x**5-38.542*x**6+12.067*x**7-
     .           log10(CFSTAR))*q
         endif
         if(tem .lt. 100) then
            sp = 8.16e-26/CFSTAR*sqrt(tem/1.0e3)*exp(-510/tem)
         else if(tem .lt. 1.0e3) then
            sp = 10.0**(-24.216+3.323*x-11.642*x**2-35.553*x**3-
     .           35.106*x**4-10.922*x**5-log10(CFSTAR))
         else
            sp = 10.0**(-24.216+4.205*x-1.316*x**2-1.655*x**3+
     .           4.178*x**4-0.569*x**5-3.382*x**6+1.090*x**7-
     .           log10(CFSTAR))*q
         endif
#endif
         RateChem(irateMHcH1o,lc) = so
         RateChem(irateMHcH1p,lc) = sp
C
C  2. collisions with ortho-H2
C
#if (RT_H2_RATE == 1)
         so = 0.0
         sp = 0.0
#elif (RT_H2_RATE == 2)
         so = 0.0
         sp = 0.0
#else
         so = 10.0**(-24.020+2.269*x-1.020*x**2+0.835*x**3-0.408*x**4+
     .        0.096*x**5-log10(CFSTAR))*q
         sp = 10.0**(-24.126+2.325*x-1.008*x**2+0.548*x**3-0.337*x**4+
     .        0.207*x**5-log10(CFSTAR))*q
#endif
         RateChem(irateMHcMOo,lc) = so
         RateChem(irateMHcMOp,lc) = sp
C
C  3. collisions with para-H2
C
#if (RT_H2_RATE == 1)
         so = 0.0
         sp = 0.0
#elif (RT_H2_RATE == 2)
         so = 0.0
         sp = 0.0
#else
         so = 10.0**(-23.749+1.767*x-0.586*x**2+0.311*x**3-0.175*x**4+
     .        0.185*x**5-log10(CFSTAR))*q
         sp = 10.0**(-23.890+1.855*x-0.556*x**2+0.284*x**3-0.205*x**4+
     .        0.131*x**5-log10(CFSTAR))*q
#endif
         RateChem(irateMHcMPo,lc) = so
         RateChem(irateMHcMPp,lc) = sp
C
C  4. collisions with He
C
#if (RT_H2_RATE == 1)
         so = 0.0
         sp = 0.0
#elif (RT_H2_RATE == 2)
         so = 0.0
         sp = 0.0
#else
         so = 10.0**(-23.775+2.407*x-1.234*x**2+0.740*x**3-0.259*x**4+
     .        0.121*x**5-log10(CFSTAR))*q
         sp = 10.0**(-23.489+1.821*x-0.591*x**2+0.423*x**3-0.302*x**4+
     .        0.129*x**5-log10(CFSTAR))*q
#endif
         RateChem(irateMHcG1o,lc) = so
         RateChem(irateMHcG1p,lc) = sp
C
C  5. collisions with Hp
C
#if (RT_H2_RATE == 1)
         so = 0.0
         sp = 0.0
#elif (RT_H2_RATE == 2)
         so = 0.0
         sp = 0.0
#else
         so = 10.0**(-21.707+1.390*x-0.350*x**2+0.075*x**3-0.232*x**4+
     .        0.069*x**5-log10(CFSTAR))*q
         sp = 10.0**(-21.757+1.400*x-0.372*x**2+0.062*x**3-0.372*x**4+
     .        0.233*x**5-log10(CFSTAR))*q
#endif
         RateChem(irateMHcH2o,lc) = so
         RateChem(irateMHcH2p,lc) = sp
C
C  6. collisions with e
C
#if (RT_H2_RATE == 1)
         so = 0.0
         sp = 0.0
#elif (RT_H2_RATE == 2)
         so = 0.0
         sp = 0.0
#else
         so = 10.0**(-21.703+0.761*x+0.506*x**2+0.050*x**3-0.104*x**4-
     .        0.036*x**5-log10(CFSTAR))*exp(-845/tem)*q
         if(tem .lt. 1.0e3) then
            sp = 10.0**(-22.818+0.957*x+0.793*x**2+0.568*x**3-
     .           0.279*x**4+0.056*x**5-log10(CFSTAR))*exp(-510/tem)
         else
            sp = q*10.0**(-22.818+0.669*x-7.119*x**2-11.177*x**3+
     .           7.047*x**4-1.647*x**5-log10(CFSTAR))*exp(-510/tem)
         endif
#endif
         RateChem(irateMHcFEo,lc) = so
         RateChem(irateMHcFEp,lc) = sp
C
C  Other H2 rates
C
         RateChem(irateMHc08,lc) = 4.476*(1.602e-12/CFSTAR)*
     .        (RateChem(irateMHi08,lc)+RateChem(irateMHi23,lc))
      
         RateChem(irateMHc09,lc) = 4.476*(1.602e-12/CFSTAR)*
     .        RateChem(irateMHi09,lc)

         RateChem(irateMHc10,lc) = 4.476*(1.602e-12/CFSTAR)*
     .        RateChem(irateMHi10,lc)

         RateChem(irateMHc14,lc) = 0.75*(1.602e-12/CFSTAR)*
     .        RateChem(irateMHi14,lc)

         RateChem(irateMHc15,lc) = 0.75*(1.602e-12/CFSTAR)*
     .        RateChem(irateMHi15,lc)
C
C  High density regime
C
#ifdef RT_HIGH_DENSITY

         x = alt10 - 4
C
C  Shapiro & Kang 87
C
         if(tem .gt. 1.0e3) then
            RateHigh(irateHDi09,lc) = 1.20e-9*exp(-5.24e4/tem)*
     .           10.0**(4.0-0.416*x-0.327*x**2)
            RateHigh(irateHDi10,lc) = 1.30e-9*exp(-5.33e4/tem)*
     .           10.0**(4.845-1.3*x+1.62*x**2)
         else
            RateHigh(irateHDi09,lc) = 0
            RateHigh(irateHDi10,lc) = 0
         endif
C
C  30. H + H + H -> H2 + H
C
#if (RT_H2_RATE == 1)
         s = 5.5e-29/tem
#elif (RT_H2_RATE == 2)
         s = 0.0
#else
         if(tem .lt. 300) then
            s = 1.14e-31/tem**0.38
         else
            s = 3.9e-30/tem
         endif
#endif
         RateHigh(irateHDi30,lc) = s
         RateHigh(irateHDi31,lc) = s/8.0
         RateHigh(irateHDi32,lc) = 6.9e-32/tem**0.4

         RateHigh(irateHDc02,lc) = 1.0e6/srt*(1.4+
     .        1.6*exp(-(400/tem)**2))/CFSTAR
C
C  Critical density (Hollenbach & McKee, 1979, ApJS, 41, 555)
C
         RateHigh(irateHDcNC,lc) = 0.7e6/srt*exp(1.2e3/(tem+1.2e3))

#endif /* RT_HIGH_DENSITY */
#endif /* RT_CHEMISTRY */

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
            write(0,*) 'RT::InitTables: fatal error, ',
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
         wpTab(iptabCiLW,lr) = csMHLW(lr)
#ifdef RT_CHEMISTRY
         wpTab(iptabCi27,lr) = csMH27(lr)
         wpTab(iptabCi28,lr) = csMH28(lr)
         wpTab(iptabCi29,lr) = csMH29(lr)
         wpTab(iptabCi30,lr) = csMH30(lr)
         wpTab(iptabCi31,lr) = csMH31(lr)
         wpTab(iptabCi32,lr) = csMH32(lr)
#endif /* RT_CHEMISTRY */
      enddo

#ifdef RT_UV
      uvFac(iptabPiH1) = csUVn(lrH1tr)
      uvFac(iptabPhH1) = csUVn(lrH1tr)
      uvFac(iptabPiG1) = csUVn(lrG1tr)
      uvFac(iptabPhG1) = csUVn(lrG1tr)
      uvFac(iptabPiG2) = csUVn(lrG2tr)
      uvFac(iptabPhG2) = csUVn(lrG2tr)
      uvFac(iptabCiLW) = 1.0
#ifdef RT_CHEMISTRY
      uvFac(iptabCi27) = 1.0
      uvFac(iptabCi28) = 1.0
      uvFac(iptabCi29) = csUVn(lrH1tr)
      uvFac(iptabCi30) = csUVn(lrG1tr)
      uvFac(iptabCi31) = csUVn(lrH1tr)
      uvFac(iptabCi32) = csUVn(lrH1tr)
#endif /* RT_CHEMISTRY */
#endif /* RT_UV */

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
      dimension pUV(iptabDim), qUV(iptabDim)
      dimension pH1(iptabDim), qH1(iptabDim)
      dimension pG1(iptabDim), qG1(iptabDim)
      dimension pG2(iptabDim), qG2(iptabDim)
C
      include 'frt_inline.inc'
C
      iptabMax = 7

      sumc1 = 0.0
      sumc2 = 0.0
      sumc3 = 0.0
      sumi1 = 0.0
      sumi2 = 0.0
      sumi3 = 0.0
      sumi4 = 0.0

      do lr=1,NRAD
         sumc1 = sumc1 + angAvg(lr)*csH1(lr)*(Txi(lr)-TH1)
         sumc2 = sumc2 + angAvg(lr)*csG1(lr)*(Txi(lr)-TG1)
         sumc3 = sumc3 + angAvg(lr)*csG2(lr)*(Txi(lr)-TG2)
         sumi1 = sumi1 + angAvg(lr)*csH1(lr)
         sumi2 = sumi2 + angAvg(lr)*csG1(lr)
         sumi3 = sumi3 + angAvg(lr)*csG2(lr)
         sumi4 = sumi4 + angAvg(lr)*csO8(lr)
      enddo

      pRate0(iptabPhH1) = max(0.0,facPhot*sumc1)
      pRate0(iptabPhG1) = max(0.0,facPhot*sumc2)
      pRate0(iptabPhG2) = max(0.0,facPhot*sumc3)
      pRate0(iptabPiH1) = max(0.0,facPhot*sumi1)
      pRate0(iptabPiG1) = max(0.0,facPhot*sumi2)
      pRate0(iptabPiG2) = max(0.0,facPhot*sumi3)
      pRate0(iptabPiO8) = max(0.0,facPhot*sumi4)

      sumiLW = 0.0
      do lr=1,NRAD
         sumiLW = sumiLW + angAvg(lr)*csMHLW(lr)
      enddo
      pRate0(iptabCiLW) = max(0.0,facPhot*sumiLW)

#ifdef RT_CHEMISTRY

      iptabMax = 13

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
C  Redo Lyman-Werner bands if needed
C
#ifdef RT_LWBANDS
      sumiLW = 0.0
      do lr=1,NRLW
         sumiLW = sumiLW + angLW(lr)*csLWMH(lr)
      enddo
      sumiLW = sumiLW*XILWST/XISTEP
      pRate0(iptabCiLW) = max(0.0,facPhot*sumiLW)
#endif /* RT_LWBANDS */

#endif /* RT_CHEMISTRY */

#ifdef RT_TRANSFER

      sumc1 = 0.0
      sumc2 = 0.0
      sumc3 = 0.0
      sumi1 = 0.0
      sumi2 = 0.0
      sumi3 = 0.0
      sumi4 = 0.0

      do lr=1,NRAD
         sumc1 = sumc1 + guvEff(lr)*csH1(lr)*(Txi(lr)-TH1)
         sumc2 = sumc2 + guvEff(lr)*csG1(lr)*(Txi(lr)-TG1)
         sumc3 = sumc3 + guvEff(lr)*csG2(lr)*(Txi(lr)-TG2)
         sumi1 = sumi1 + guvEff(lr)*csH1(lr)
         sumi2 = sumi2 + guvEff(lr)*csG1(lr)
         sumi3 = sumi3 + guvEff(lr)*csG2(lr)
         sumi4 = sumi4 + guvEff(lr)*csO8(lr)
      enddo

      qRate0(iptabPhH1) = max(0.0,facPhot*sumc1*facLX)
      qRate0(iptabPhG1) = max(0.0,facPhot*sumc2*facLX)
      qRate0(iptabPhG2) = max(0.0,facPhot*sumc3*facLX)
      qRate0(iptabPiH1) = max(0.0,facPhot*sumi1*facLX)
      qRate0(iptabPiG1) = max(0.0,facPhot*sumi2*facLX)
      qRate0(iptabPiG2) = max(0.0,facPhot*sumi3*facLX)
      qRate0(iptabPiO8) = max(0.0,facPhot*sumi4*facLX)

      sumiLW = 0.0
      do lr=1,NRAD
         sumiLW = sumiLW + guvEff(lr)*csMHLW(lr)
      enddo
      qRate0(iptabCiLW) = max(0.0,facPhot*sumiLW*facLX)

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

      qRate0(iptabCi27) = max(0.0,facPhot*sumi27*facLX)
      qRate0(iptabCi28) = max(0.0,facPhot*sumi28*facLX)
      qRate0(iptabCi29) = max(0.0,facPhot*sumi29*facLX)
      qRate0(iptabCi30) = max(0.0,facPhot*sumi30*facLX)
      qRate0(iptabCi31) = max(0.0,facPhot*sumi31*facLX)
      qRate0(iptabCi32) = max(0.0,facPhot*sumi32*facLX)
C
C  Redo Lyman-Werner bands if needed
C
#ifdef RT_LWBANDS
      sumiLW = 0.0
      do lr=1,NRLW
         sumiLW = sumiLW + guvLWe(lr)*csLWMH(lr)
      enddo
      sumiLW = sumiLW*XILWST/XISTEP
      qRate0(iptabCiLW) = max(0.0,facPhot*sumiLW*facLX)
#endif /* RT_LWBANDS */

#endif /* RT_CHEMISTRY */

#endif /* RT_TRANSFER */

C
C  Main tables: we work in the approximation that the UV absorption
C  can be added posteriori to a 3D table.
C
#ifndef RT_NO_TABLE

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
C  ***********************************************
C
      do j=1,iptabMax
#ifdef RT_EXTERNAL_BACKGROUND
         pUV(j) = 0.0
#endif
#ifdef RT_TRANSFER
         qUV(j) = 0.0
#endif
      enddo

#ifdef RT_CHEMISTRY
#ifdef RT_EXTERNAL_BACKGROUND
      call frtIntegrateRates_UV(0,iptabDim,wpTab,
     .     iptabOT,iptabMax,pUV,0.0)
#endif
#ifdef RT_TRANSFER
      call frtIntegrateRates_UV(1,iptabDim,wpTab,
     .     iptabOT,iptabMax,qUV,0.0)
#endif
#endif /* RT_CHEMISTRY */

C
C  3D table for Ionizing radiation
C  ***********************************************
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lo1,lo2,lo3,tauH1l,tauG1l,tauG2l,itind,itind2,j,
C$OMP+ pH1,pG1,pG2,qH1,qG1,qG2)
C$OMP+SHARED(tauH1m,acTAU,csH1tr,csG1tr,csG2tr,tauG1m,tauG2m,
C$OMP+ itloc2,pTab3,qTab3,pUV,qUV,wpTab,facPhot,facLX,iptabMax)
      do lo1=1,NOPT1

         do j=1,iptabMax
#ifdef RT_EXTERNAL_BACKGROUND
            pH1(j) = pUV(j)
#endif
#ifdef RT_TRANSFER
            qH1(j) = qUV(j)
#endif
         enddo

         if(lo1 .lt. NOPT1) then

            tauH1l = tauH1m*acTAU(lo1)
#ifdef RT_EXTERNAL_BACKGROUND
            call frtIntegrateRates_H1(0,iptabDim,wpTab,
     .           iptabH1,iptabMax,pH1,0.0,tauH1l)
#endif
#ifdef RT_TRANSFER
            call frtIntegrateRates_H1(1,iptabDim,wpTab,
     .           iptabH1,iptabMax,qH1,0.0,tauH1l)
#endif

         else

            tauH1l = 2*VLARGE

         endif

         do lo2=1,NOPT1

            do j=1,iptabMax
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
               call frtIntegrateRates_G1(0,iptabDim,wpTab,
     .              iptabG1,iptabMax,pG1,0.0,tauH1l,tauG1l)
#endif
#ifdef RT_TRANSFER
               call frtIntegrateRates_G1(1,iptabDim,wpTab,
     .              iptabG1,iptabMax,qG1,0.0,tauH1l,tauG1l)
#endif

            else

               tauG1l = 2*VLARGE

            endif

            itind2 = itloc2(lo2,lo1)

            do lo3=1,NOPT1

               do j=1,iptabMax
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
                  call frtIntegrateRates_G2(0,iptabDim,wpTab,
     .                 iptabG2,iptabMax,pG2,0.0,tauH1l,tauG1l,tauG2l)
#endif
#ifdef RT_TRANSFER
                  call frtIntegrateRates_G2(1,iptabDim,wpTab,
     .                 iptabG2,iptabMax,qG2,0.0,tauH1l,tauG1l,tauG2l)
#endif

               endif

               itind = lo3 + itind2
               do j=1,iptabMax
#ifdef RT_INTERPOLLOG
#ifdef RT_EXTERNAL_BACKGROUND
                  if(pG2(j) .gt. 0.0) then
                     pTab3(j,itind) = log(VSMALL+facPhot*pG2(j))
                  else
                     pTab3(j,itind) = -1000.0
                  endif
#endif
#ifdef RT_TRANSFER
                  if(qG2(j) .gt. 0.0) then
                     qTab3(j,itind) = log(VSMALL+facLX*facPhot*qG2(j))
                  else
                     qTab3(j,itind) = -1000.0
                  endif
#endif
#else  /* RT_INTERPOLLOG */
#ifdef RT_EXTERNAL_BACKGROUND
                  pTab3(j,itind) = max(0.0,facPhot*pG2(j))
#endif
#ifdef RT_TRANSFER
                  qTab3(j,itind) = max(0.0,facLX*facPhot*sgpG2(j))
#endif
#endif /* RT_INTERPOLLOG */
               enddo

            enddo
         enddo
      enddo

#endif /* RT_NO_TABLE */

#ifdef RT_XRAYS
      call frtFillRadiationTablesSE
#endif

      return
      end
C
C
C     
      subroutine frtPhotoRatesTable(var,rf,itab,pRate,jmax1)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_TRANSFER
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
#endif
      dimension var(*), rf(*)
      dimension itab(*)
      dimension pRate(*)
C
      dimension qRate(iptabDim)
#ifdef RT_XRAYS
      include 'frt_xrays.inc'
      dimension pRateSE(istabDim), qRateSE(istabDim)
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
     .     iptabDim,pTab3,jmax,pRate,iOk)

      if(iOk .eq. 1) then
         itab(1) = itab(1) + 1

#ifdef RT_UV
C
C  This is temporary fix - disable dust opacity for external background
C
*         do j=1,jmax
*            pRate(j) = pRate(j)*exp(-tauUVl*uvFac(j))
*         enddo
#endif
      else

#endif /* RT_NO_TABLE */

         do j=1,jmax
            pRate(j) = 0.0
         enddo

#if defined(RT_NARROW_TABLE) || defined(RT_NO_TABLE)

         itab(2) = itab(2) + 1

#ifdef RT_CHEMISTRY
         call frtIntegrateRates_UV(0,iptabDim,wpTab,
     .        iptabOT,jmax,pRate,tauUVl)
#endif /* RT_CHEMISTRY */

         call frtIntegrateRates_H1(0,iptabDim,wpTab,
     .        iptabH1,jmax,pRate,tauUVl,tauH1l)

         call frtIntegrateRates_G1(0,iptabDim,wpTab,
     .        iptabG1,jmax,pRate,tauUVl,tauH1l,tauG1l)

         call frtIntegrateRates_G2(0,iptabDim,wpTab,
     .        iptabG2,jmax,pRate,tauUVl,tauH1l,tauG1l,tauG2l)

         do j=1,jmax
            pRate(j) = facPhot*max(0.0,pRate(j))
         enddo

#endif /* RT_NARROWTABLE || RT_NO_TABLE */

#ifndef RT_NO_TABLE
      endif
#endif /* RT_NO_TABLE */

#ifdef RT_UV
      pRate(iptabCiLW) = pRate0(iptabCiLW)*exp(-tauUVl)
#else
      pRate(iptabCiLW) = pRate0(iptabCiLW)
#endif

#else  /* RT_EXTERNAL_BACKGROUND */

      do j=1,iptabDim
         pRate(j) = 0.0
      enddo

#endif /* RT_EXTERNAL_BACKGROUND */
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
     .     iptabDim,qTab3,jmax,qRate,iOk)

      if(iOk .eq. 1) then
         itab(1) = itab(1) + 1

#ifdef RT_UV
         do j=1,jmax
            qRate(j) = qRate(j)*exp(-tauUVl*uvFac(j))
         enddo
#endif
      else

#endif /* RT_NO_TABLE */

         do j=1,jmax
            qRate(j) = 0.0
         enddo

#if defined(RT_NARROW_TABLE) || defined(RT_NO_TABLE)

         itab(2) = itab(2) + 1

#ifdef RT_CHEMISTRY
         call frtIntegrateRates_UV(1,iptabDim,wpTab,
     .        iptabOT,jmax,qRate,tauUVl)
#endif /* RT_CHEMISTRY */

         call frtIntegrateRates_H1(1,iptabDim,wpTab,
     .        iptabH1,jmax,qRate,tauUVl,tauH1l)

         call frtIntegrateRates_G1(1,iptabDim,wpTab,
     .        iptabG1,jmax,qRate,tauUVl,tauH1l,tauG1l)

         call frtIntegrateRates_G2(1,iptabDim,wpTab,
     .        iptabG2,jmax,qRate,tauUVl,tauH1l,tauG1l,tauG2l)

         do j=1,jmax
            qRate(j) = facLX*facPhot*max(0.0,qRate(j))
         enddo

#endif /* RT_NARROWTABLE || RT_NO_TABLE */

#ifndef RT_NO_TABLE
      endif
#endif /* RT_NO_TABLE */

#ifdef RT_UV
      qRate(iptabCiLW) = qRate0(iptabCiLW)*exp(-tauUVl)
#else
      qRate(iptabCiLW) = qRate0(iptabCiLW)
#endif

#endif /* RT_TRANSFER */


#ifdef RT_XRAYS

      call frtPhotoRatesSE(var,rf,itab,pRateSE,qRateSE)

#ifdef RT_EXTERNAL_BACKGROUND
      pRate(iptabPiH1) = pRate(iptabPiH1) + 
     .     var(ivarXH1)*pRateSE(istabPiH1H1) +
     .     var(ivarXG1)*pRateSE(istabPiH1G1) +
     .     var(ivarXG2)*pRateSE(istabPiH1G2)
      pRate(iptabPiG1) = pRate(iptabPiG1) +
     .     var(ivarXH1)*pRateSE(istabPiG1H1) +
     .     var(ivarXG2)*pRateSE(istabPiG1G1) +
     .     var(ivarXG2)*pRateSE(istabPiG1G2)
      pRate(iptabPhH1) = max(0.0,
     .     pRate(iptabPhH1)-pRateSE(istabPhH1))
      pRate(iptabPhG1) = max(0.0,
     .     pRate(iptabPhG1)-pRateSE(istabPhG1))
#endif

#ifdef RT_TRANSFER
      qRate(iptabPiH1) = qRate(iptabPiH1) + 
     .     var(ivarXH1)*qRateSE(istabPiH1H1) +
     .     var(ivarXG1)*qRateSE(istabPiH1G1) +
     .     var(ivarXG2)*qRateSE(istabPiH1G2)
      qRate(iptabPiG1) = pRate(iptabPiG1) +
     .     var(ivarXH1)*qRateSE(istabPiG1H1) +
     .     var(ivarXG2)*qRateSE(istabPiG1G1) +
     .     var(ivarXG2)*qRateSE(istabPiG1G2)
      qRate(iptabPhH1) = max(0.0,qRate(iptabPhH1) -
     .     qRateSE(istabPhH1))
      qRate(iptabPhG1) = max(0.0,qRate(iptabPhG1) -
     .     qRateSE(istabPhG1))
#endif /* RT_TRANSFER */

#endif /* RT_XRAYS */

      var(ivarCFH1) = 1.0
      var(ivarCFG1) = 1.0
      var(ivarCFG2) = 1.0

#if (RT_CFI == 1)
#ifdef RT_EXTERNAL_BACKGROUND
      ch = CHeat(var)
      pRate(iptabPhH1) = pRate(iptabPhH1)*ch*var(ivarBias)
      pRate(iptabPhG1) = pRate(iptabPhG1)*ch*var(ivarBias)
      pRate(iptabPhG2) = pRate(iptabPhG2)*ch*var(ivarBias)
      var(ivarCFH1) = CIonH1(var)
      var(ivarCFG1) = CIonG1(var)
      var(ivarCFG2) = CIonG2(var)
      pRate(iptabPiH1) = pRate(iptabPiH1)*var(ivarCFH1)*var(ivarBias)
      pRate(iptabPiG1) = pRate(iptabPiG1)*var(ivarCFG1)*var(ivarBias)
      pRate(iptabPiG2) = pRate(iptabPiG2)*var(ivarCFG2)*var(ivarBias)
#endif
#endif

#ifdef RT_TRANSFER
      do j=1,jmax
         pRate(j) = pRate(j) + rf(irfNorm)*qRate(j)
      enddo
      pRate(iptabCiLW) = pRate(iptabCiLW) + rf(irfNorm)*
     .     qRate(iptabCiLW)
#endif /* RT_TRANSFER */

#if (RT_CFI == 2)
      ch = CHeat(var)
      pRate(iptabPhH1) = pRate(iptabPhH1)*ch*var(ivarBias)
      pRate(iptabPhG1) = pRate(iptabPhG1)*ch*var(ivarBias)
      pRate(iptabPhG2) = pRate(iptabPhG2)*ch*var(ivarBias)
      var(ivarCFH1) = CIonH1(var)
      var(ivarCFG1) = CIonG1(var)
      var(ivarCFG2) = CIonG2(var)
      pRate(iptabPiH1) = pRate(iptabPiH1)*var(ivarCFH1)*var(ivarBias)
      pRate(iptabPiG1) = pRate(iptabPiG1)*var(ivarCFG1)*var(ivarBias)
      pRate(iptabPiG2) = pRate(iptabPig2)*var(ivarCFG2)*var(ivarBias)
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
      subroutine frtIntegrateRates_UV(id,nmax,wi,j1,jn,s,tauUVl)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)
      include 'frt_inline.inc'

      do lr=lrmax+1,NRAD
         if(id .eq. 0) then
            fac = angEff(lr)
         else
            fac = guvEff(lr)
         endif
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      if(tauUVl .gt. VLARGE) return

      do lr=1,lrH1tr-1
         if(id .eq. 0) then
#ifdef RT_UV
            taul = csUVn(lr)*tauUVl
            fac = angEff(lr)*QFuni(taul)
#else
            fac = angEff(lr)
#endif
         else
#ifdef RT_UV
            taul = csUVn(lr)*tauUVl
            fac = guvEff(lr)*QFloc(taul)
#else
            fac = guvEff(lr)
#endif
         endif
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return 
      end
C
C
C
      subroutine frtIntegrateRates_H1(id,nmax,wi,j1,jn,s,tauUVl,
     .     tauH1l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)
      include 'frt_inline.inc'

      if(tauH1l .gt. VLARGE) return

      do lr=lrH1tr,lrG1tr-1
         taul = csH1n(lr)*tauH1l
#ifdef RT_UV
         if(tauUVl .lt. VLARGE) then
            taul = taul + csUVn(lr)*tauUVl
         endif
#endif
         if(id .eq. 0) then
            fac = angEff(lr)*QFuni(taul)
         else
            fac = guvEff(lr)*QFloc(taul)
         endif
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C
      subroutine frtIntegrateRates_G1(id,nmax,wi,j1,jn,s,tauUVl,
     .     tauH1l,tauG1l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)
      include 'frt_inline.inc'

      if(tauG1l .gt. VLARGE) return

      do lr=lrG1tr,lrG2tr-1
         taul = csG1n(lr)*tauG1l
         if(tauH1l .lt. VLARGE) then
            taul = taul + csH1n(lr)*tauH1l
         endif
#ifdef RT_UV
         if(tauUVl .lt. VLARGE) then
            taul = taul + csUVn(lr)*tauUVl
         endif
#endif
         if(id .eq. 0) then
            fac = angEff(lr)*QFuni(taul)
         else
            fac = guvEff(lr)*QFloc(taul)
         endif
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C
      subroutine frtIntegrateRates_G2(id,nmax,wi,j1,jn,s,tauUVl,
     .     tauH1l,tauG1l,tauG2l)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension s(*), wi(nmax,*)
      include 'frt_inline.inc'

      if(tauG2l .gt. VLARGE) return

      do lr=lrG2tr,lrmax
         taul = csG2n(lr)*tauG2l
         if(tauH1l .lt. VLARGE) then
            taul = taul + csH1n(lr)*tauH1l
         endif
         if(tauG1l .lt. VLARGE) then
            taul = taul + csG1n(lr)*tauG1l
         endif
#ifdef RT_UV
         if(tauUVl .lt. VLARGE) then
            taul = taul + csUVn(lr)*tauUVl
         endif
#endif
         if(id .eq. 0) then
            fac = angEff(lr)*QFuni(taul)
         else
            fac = guvEff(lr)*QFloc(taul)
         endif
         do j=j1,jn
            s(j) = s(j) + fac*wi(j,lr)
         enddo
      enddo

      return
      end
C
C
C     
#ifndef RT_NO_TABLE
      subroutine frtSampleTable3(tauH1l,tauG1l,tauG2l,
     .     nmax,table,nout,out,iOk)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension table(nmax,*), out(*)

      if(tauH1l .gt. VLARGE) then
         lH1tab = 2
      else if(tauH1l .lt. tH1max) then
         lH1tab = 1
      else
         lH1tab = 0
      endif

      if(tauG1l .gt. VLARGE) then
         lG1tab = 2
      else if(tauG1l .lt. tG1max) then
         lG1tab = 1
      else
         lG1tab = 0
      endif

      if(tauG2l .gt. VLARGE) then
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
#endif /* RT_NO_TABLE */
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
