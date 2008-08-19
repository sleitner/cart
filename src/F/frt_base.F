#include "defs.h"      
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
C
C
C
      subroutine frtInitRun(Yp)
      include 'frt_base.inc'
      include 'frt_tables.inc'
C
C  Default configuration
C
      call frtInitRunGeneric(Yp,10.0,10.0)

      return
      end
C
C
C
      subroutine frtInitRunGeneric(Yp,TemFloor,ClumpH2)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
#ifdef RT_XRAYS
      include 'frt_xrays.inc'
#endif
C
#ifdef RT_DUST
      write(0,*) 'frtInitRun: fatal error, ',
     .   'RT_DUST is on - this is not fully debugged yet.'
      stop
#endif
C
      if(Yp.lt.1.0e-10 .or. Yp.gt.0.99) then
         write(0,*) 'frtInitRun: fatal error, ',
     .        'Yp must be between 1e-10 and 0.99.'
         stop
      endif
      XH = 1 - Yp
      XG = 0.25*Yp
      XFEMAX = XH + 2*XG
C
C  Someparameters:
C  1. H2 formation model: clumping factor in molecular clouds
C  2. Temperature floor
C
      d_clump = ClumpH2
      Tmin = TemFloor
C     
C  Table initialization
C
      call frtInitTables
      call frtInitSourceSpectrum

#ifdef RT_CHEMISTRY
#ifdef RT_LWBANDS
C
C  Call Massimo's routine
C
      call frtLW_InitRun()
      call frtLW_ReadData(guvLW0)
#endif
#endif
C
C  Init radiation field factors
C
      do lr=1,NRAD
         tauhat_xi(lr) = 0.0
         qfdhat_xi(lr) = 1.0
      enddo

      return
      end
C
C
C
      subroutine frtStepBeginART(as,hs,ts,Om0,h100,r0,t0)
      include 'frt_base.inc'
      include 'frt_tables.inc'
C
#ifdef RT_TEST
      call frtTestStepBegin(Om0,h100,r0,t0)
#endif
C
C  physical n_b/rho
C
      ud = 1.123e-5*Om0*h100**2/as**3
C
C  Time unit (in s)
C
      ut = t0*3.15e7*as**2
C
C  a*r0*1Mpc/h
C
      ul = 3.086e24*r0/h100*as

      call frtStepBegin(ud,ul,ut,as,hs,ts,0)

      return
      end
C
C
C
      subroutine frtStepBegin(ud,ul,ut,as,hs,ts,mab)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      dimension guvloc(NRADL)
C
C  Parameters of the step
C
      uDen = ud
      uLen = ul
      uTime = ut
      aStep = as
      hStep = hs
      tStep = ts
      modeAB = mab

#ifdef RT_TEST
C
C  RT-TEST +++++++++++++++++
C
      modeAB = 2
C
C  RT-TEST -----------------
C
#endif
C
C  Some constants
C
C ----------------------------------------------------
C
C  Additional units (and conversions)
C
      uRate = uDen*uTime
      uColumn = CSSTAR*uDen*uLen
C
C  CMB temperature (in K)
C
      TemCMB = 2.728/aStep
C
C  Destruction of H2 by CMB photons
C
      facH2CMB = 1.1e-1*TemCMB**2.13*exp(-8823.0/TemCMB)
C
C   dEin/dt = facCool * CoolFun * rho 
C   Only need to switch from proper time unit and from ergs to K
C
      facCool = uRate*(CFSTAR/1.381e-16)
C
C  P_i = facPi*Sum_j sigma_j ang_j/aStep**3 
C  facPi = c * CSSTAR *dxi
C
      facPi = 2.880e-13*XISTEP/aStep**3
C
C  P_h = facPh*t0 * Sum_j sigma_j T_j ang_j
C
      facPh = facPi*uTime
C
C  Compton cooling: dEin/dt (K/code units) = facCC*Xfe*(T-TCMB)
C
      facCC = 7.372e-22*uTime*TemCMB**4
C
C  Cell crossing time in units of c [uLen/(c*uTime)]
C
      facLX = (uLen/2.998e10)/uTime
C
C  Lyman-alpha pumping - Tozzi et al formula. Jordi says it is wrong.
C  dEa/dt = 1.735e-11 K cm^3/c ang(Lya)
C  Here t0 is in years.
C
#ifdef RT_LYMAN_ALPHA_HEATING
      facLya = 1.735e-11*uTime*ang(lrH1rB)
#else
      facLya = 0.0
#endif
C
C  Recombination case A/B:
C    tau = XH1*fOSC/hnu_eV*tau0
C    tau0 = c sigmaAB CSTAR nb / H / sqrt(Pi)
C
      if(modeAB.eq.0 .and. hStep.gt.1.0e-20) then
         facRecB = 2.998e10*CSSTAR*uDen/(uTime*hStep)
      else
         facRecB = 0.0
      endif
C
C Get the source spectrum
C
      call frtComputeSourceSpectrum(aStep,guvloc)
C
C  Normalize
C
      uComovingDen = uDen*aStep**3
      do lr=1,NRAD
         guveff(lr) = uComovingDen*guvloc(lr)
      enddo

#ifdef RT_CHEMISTRY
#ifdef RT_LWBANDS

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(guvLW0,guvLWe,guveff,lr99tr)
      do lr=1,NRLW
         guvLWe(lr) = guvLW0(lr)*guveff(lr99tr)
      enddo
C     
C  Call Massimo's function
C         
      call  frtLW_SigmaH2(hStep,tStep,csLWMH)
C
C  Renormalize
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr)
C$OMP+SHARED(csLWMH,csLWHI)
      do lr=1,NRLW
         csLWMH(lr) = csLWMH(lr)/CSSTAR
         csLWHI(lr) = csLWHI(lr)/CSSTAR
      enddo

#endif
#endif
C
C  Lower limit and pivot points for cooling ODE solvers
C
      yLimit(ieqEin) = 1.0e-2*TemCMB
      yLimit(ieqXMH) = 1.0e-20*XH
      yLimit(ieqXH1) = 1.0e-15*XH
      yLimit(ieqXH2) = 1.0e-15*XH
      yLimit(ieqXG1) = 1.0e-15*XG
      yLimit(ieqXG2) = 1.0e-15*XG
      yLimit(ieqXG3) = 1.0e-15*XG
      yLimit(ieqXHp) = 1.0e-25*XH
      yLimit(ieqXHm) = 1.0e-25*XH

      ySmall(ieqEin) = TemCMB
      ySmall(ieqXMH) = 1.0e-8*XH
      ySmall(ieqXH1) = 1.0e-6*XH
      ySmall(ieqXH2) = 1.0e-6*XH
      ySmall(ieqXG1) = 1.0e-6*XG
      ySmall(ieqXG2) = 1.0e-6*XG
      ySmall(ieqXG3) = 1.0e-6*XG
      ySmall(ieqXHp) = 1.0e-8*XH
      ySmall(ieqXHm) = 1.0e-8*XH
C
C  Reset global field counters
C
      call frtResetGlobalFields

      return
      end
C
C
C
      subroutine frtStepEnd(id,vol,avg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      include 'frt_global.inc'
      dimension avg(*)
C
      dimension axi(NRAD), sxi(NRAD), xib(NRAD), anb(NRAD)
#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)
      dimension axiLW(NRLW), sxiLW(NRLW), xibLW(NRLW), anbLW(NRLW)
#endif
      character*(99) str1, str2, str3
C
      common/RT_RT_UV/ duvin, axiH1tr, angH1tr
C
C  Unpack global quantities
C
      duvin = avg(1)
#ifdef RT_CHEMISTRY
      XH1_va = avg(2)
      XMH_va = avg(3)
#endif

      s = vol/tStep
      do j=1,NRAD
         gBrem(j) = gBrem(j)*s
      enddo
      do j=1,ioMax1
         gData(j) = gData(j)*s
      enddo
C
C  Diagnostics
C
#ifdef RT_OUTPUT

      if(id .ge. 0) then

         ftabHit = gInfo(infoTabH)
         ftabMiss = gInfo(infoTabM)

#ifdef RT_XRAYS
         ftabFill = 4.0*NOPT**3
#else
         ftabFill = 2.0*NOPT**3
#endif

         if(id .gt. 0) then
            write(str3,*) id
            ibeg = 1
            iok = 0
            do i=1,99
               if(iok.eq.0 .and. str3(i:i).eq.' ') then
                  ibeg = ibeg + 1
               else
                  iok = 1
               endif
            enddo
            iend = ibeg + index(str3(ibeg:99),' ') - 2
            str1 = 'RT_Tables['//str3(ibeg:iend)//'] '
            str2 = 'RT_Cooling['//str3(ibeg:iend)//'] '
         else
            str1 = 'RT_Tables '
            str2 = 'RT_Cooling '
         endif
         iend1 = index(str1,' ') - 1
         iend2 = index(str2,' ') - 1
         
         if(ftabHit .gt. 0.5) then
            write(6,*) str1(1:iend1)//': miss/hit ', 
     .           ftabMiss/ftabHit, ', fill/hit ', ftabFill/ftabHit
         endif

         if(gInfo(infoCall) .gt. gInfo(infoFail)) then
            write(6,*) str2(1:iend2)//': avg # of steps ', 
     .           gInfo(infoNAvg)/(gInfo(infoCall)-gInfo(infoFail))
            write(6,*) str2(1:iend2)//': max # of steps ', 
     .           gInfo(infoNMax)
         endif

         if(gInfo(infoFail) .gt. 0.0d0) then
            write(6,*) str2(1:iend2)//': fail/call ', 
     .           gInfo(infoFail)/gInfo(infoCall)
         endif

      endif
#endif

#ifndef RT_TEST
C
C  Coefficient for Bremstrahlung: 
C     dn_xi/dt = 1/V e^{-h\nu/kT} 1/kT dE_tot/dt
C     gBrem(xi) = Sum dV rho dT/dt 1/T e^{-T_xi/T}
C     CRAD1 = n_b/V
C
      CRAD1 = uDen
C
C  Coefficient for absorption:
C    c*CSSTAR*n_b/uTime = uColumn/facLX
C
      CRAD2 = uColumn/facLX
C
C  Advance the average radiation field:
C  Source function
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(CRAD1,gBrem,duvin,guveff,sxi)
      do lr=1,NRAD
         sxi(lr) = CRAD1*gBrem(lr) + duvin*guveff(lr)
      enddo
C
C  Absorption
C
      do lr=1,NRAD
         axi(lr) = 0.0
      enddo
C
C Optically thin absortion from Lyman Series and dust
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(axi,CRAD2,csDust,lrmin,lrmax,csLyS,gData)
      do lr=lrmin,lrmax
         axi(lr) = axi(lr) + CRAD2*csLyS(lr)*gData(ioCXH1)
#ifdef RT_DUST
         axi(lr) = axi(lr) + CRAD2*csDust(lr)*gData(ioDust)
#endif
      enddo

#ifdef RT_TRANSFER
      
      call frtTransferComputeAverageAbs(axi)

#else

C
C Optically thin absortion
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(axi,CRAD2,lrmin,lrmax,csH1,csG1,csG2,gData)
      do lr=lrmin,lrmax
               
         axi(lr) = axi(lr) + CRAD2*(csH1(lr)*gData(ioCXH1)+
     .        csG1(lr)*gData(ioCXG1)+csG2(lr)*gData(ioCXG2))
      enddo

#endif

#ifdef RT_CHEMISTRY
C
C  Absorption of Lyman-Wienner band photons (only included in the volume 
C  averaged equation to fix Olber's paradox).
C
      do lr=lr99tr,lrH1tr-1
         axi(lr) = axi(lr) + csMH99(lr)*XMH_va*CRAD2
      enddo
#endif

      axiH1tr = axi(lrH1tr)
      angH1tr = ang(lrH1tr)
C
C  Hydrogen and Helium recombination photons:
C
C  Recombination case A
C
      sxi(lrH1tr) = sxi(lrH1tr) + CRAD1*gData(ioRAH2)/XISTEP
      sxi(lrG1tr) = sxi(lrG1tr) + CRAD1*gData(ioRAG2)/XISTEP
      sxi(lrG2tr) = sxi(lrG2tr) + CRAD1*gData(ioRAG3)/XISTEP
C
C  Recombination case B
C
      sxi(lrH1rB) = sxi(lrH1rB) + CRAD1*gData(ioRBH2)/XISTEP
      sxi(lrG1rB) = sxi(lrG1rB) + CRAD1*gData(ioRBG2)/XISTEP
      sxi(lrG2rB) = sxi(lrG2rB) + CRAD1*gData(ioRBG3)/XISTEP
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,t12,p12,q12),
C$OMP+SHARED(axi,tStep,anb,ang,sxi)
      do lr=1,NRAD

         t12 = axi(lr)*tStep
         if(t12 .gt. 20.0) then
            p12 = 0.0
            q12 = 1.0/axi(lr)
         else if(t12 .lt. 0.001) then
            p12 = 1.0-t12
            q12 = tStep
         else
            p12 = exp(-t12)
            q12 = (1.0-p12)/axi(lr)
         endif

         anb(lr) = ang(lr)*p12 + sxi(lr)*q12 

      enddo

      shift = hStep*tStep
      lshift = int(shift/XISTEP)
      do lr=1,NRAD
         xib(lr) = xig(lr) + shift
      enddo

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,lbl,lbr,anbl,xibl,anbr,xibr),
C$OMP+SHARED(lshift,ang,anb,xib,xig)
      do lr=1,NRAD
         lbl = lr + lshift
         lbr = lbl+1
         if(lbl .gt. NRAD) then
            ang(lr) = 0.0
         else 
            anbl = anb(lbl)
            xibl = xig(lbl)
            if(lbr .gt. NRAD) then
               anbr = 0.0
            else
               anbr = anb(lbr)
            endif
            xibr = xibl + XISTEP
            ang(lr) = ((xib(lr)-xibl)*anbr+(xibr-xib(lr))*anbl)/XISTEP
         endif
      enddo
C
C  This is just in case...
C
      do lr=1,NRAD
         ang(lr) = max(0.0,ang(lr))
      enddo
C
C  Do Lyman-Werner bands separately
C
#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)

      fMH = XMH_va*CRAD2/CSSTAR
      fH1 = XH1_va*CRAD2/CSSTAR
      call frtLW_ComputeLineShifts(hStep,fMH)
      call frtLW_ComputeAbs(hStep,fH1,tStep,0.0,axiLW)
      call frtLW_ComputeSource(hStep,fH1,tStep,guvLWe,sxiLW)

      do lr=1,NRLW
         anbLW(lr) = angLW(lr)*axiLW(lr) + duvin*sxiLW(lr)*tStep
      enddo
C
      shift = 0  !!  already included 
      lshift = int(shift/XILWST)
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr)
C$OMP+SHARED(xibLW,xigLW,shift)
      do lr=1,NRLW
         xibLW(lr) = xigLW(lr) + shift
      enddo

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,lbl,lbr,anbl,xibl,anbr,xibr),
C$OMP+SHARED(lshift,angLW,anbLW,xibLW,xigLW,angH1tr,XILWST)
      do lr=1,NRLW
         lbl = lr + lshift
         lbr = lbl+1
         if(lbl .gt. NRLW) then
            angLW(lr) = 0.0
         else 
            anbl = anbLW(lbl)
            xibl = XILWST*(real(lbl)-1.0)
            if(lbr .gt. NRLW) then
               anbr = angH1tr
            else
               anbr = anbLW(lbr)
            endif
            xibr = xibl+XILWST
            angLW(lr) = ((xibLW(lr)-xibl)*anbr+(xibr-xibLW(lr))*anbl)/
     .           XILWST
         endif
      enddo
C
C  This is just in case...
C
      do lr=1,NRLW
         angLW(lr) = max(0.0,angLW(lr))
      enddo

#endif /*  defined(RT_CHEMISTRY) && defined(RT_LWBANDS) */

#endif /* RT_TEST */

      return
      end
C
C
C
      function frtQueryXH()
      include 'frt_base.inc'

      frtQueryXH = XH

      return
      end
C
C
C
      function frtQueryXHe()
      include 'frt_base.inc'

      frtQueryXHe = XG

      return
      end
#endif
      subroutine frtDummyBase(i)
      i = 0
      end
