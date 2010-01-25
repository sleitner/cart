#include "../config.h"      
#ifdef RADIATIVE_TRANSFER
C
C
C
      subroutine frtInitRun2(Ypl,Tminl,D2Gminl,ClumpH2l,CohLenH2l,
     .     fGall,fQSOl,IPOPl,IRECl)
      include 'frt_f77.inc'
C
      Yp = Ypl
      Tmin = Tminl
      D2Gin = D2Gminl
      ClumpH2 = ClumpH2l
      CohLenH2 = CohLenH2l
      fGal = fGall
      fQSO = fQSOl
      IPOP = IPOPl
      IREC = IRECl
C
      call frtInitRun
C
      return
      end
C
C
C
      subroutine frtInitRun
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
#ifdef RT_XRAYS
      include 'frt_xrays.inc'
#endif
C
      data Yp/0.24/, Tmin/10.0/, D2Gmin/1.0e-3/
      data ClumpH2/30.0/, CohLenH2/0.3/
      data fGal /1.0/, fQSO/1.0/
      data IPOP /2/, IREC/0/
C
#ifdef RT_DUST
      write(0,*) 'RT_Base: fatal error, ',
     .   'RT_DUST is on - this is not fully debugged yet.'
      call frtAbort
#endif
C
      if(Yp.lt.1.0e-10 .or. Yp.gt.0.99) then
         write(0,*) 'RT_Base: fatal error, ',
     .        'Yp must be between 1e-10 and 0.99.'
         call frtAbort
      endif
C
#ifdef RT_OUTPUT
      call frtQueryId(id)
      if(id .eq. 0) then
         open(unit=81, file='rt-run.log')
         write(81,*) 'Run:'
         write(81,*) 'Yp:       ', Yp
         write(81,*) 'Tmin:     ', Tmin
         write(81,*) 'D2Gmin:   ', D2Gmin
         write(81,*) 'ClumpH2:  ', ClumpH2
         write(81,*) 'CohLenH2: ', CohLenH2
         write(81,*) 'fGal:     ', fGal
         write(81,*) 'fQSO:     ', fQSO
         write(81,*) 'IPOP:     ', IPOP
         write(81,*) 'IREC:     ', IREC
         close(81)
      endif
#endif /* RT_OUTPUT */
C
      XH = 1 - Yp
      XG = 0.25*Yp
      XFEMAX = XH + 2*XG
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
C  Some extra initialization to avoid division by zero
C
      aExp = 1

#ifdef RT_TEST
C
C  RT-TEST +++++++++++++++++
C
      IREC = 2
C
C  RT-TEST -----------------
C
#endif

      return
      end
C
C
C
      subroutine frtStepBegin(uDenl,uLenl,uTimel,dtStepl,aExpl,Hubblel,
     .     rfAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      dimension guv(NRAD), rfAvg(*)
C
      include 'frt_inline.inc'
C
      uDen = uDenl
      uLen = uLenl
      uTime = uTimel
      dtStep = dtStepl
      aExp = aExpl
      Hubble = Hubblel
C
#ifdef RT_OUTPUT
      call frtQueryId(id)
      if(id .eq. 0) then
         open(unit=81, file='rt-run.log', access='append')
         write(81,*) 'Step:'
         write(81,*) '   uDen:   ', uDen
         write(81,*) '   uLen:   ', uLen
         write(81,*) '   uTime:  ', uTime
         write(81,*) '   aExp:   ', aExp
         write(81,*) '   Hubble: ', Hubble
         write(81,*) '   dtStep: ', dtStep
         close(81)
      endif
#endif /* RT_OUTPUT */

#ifdef RT_TEST
      call frtTestStepBegin
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

#ifdef RT_DEBUG
      call frtQueryId(id)
      if(id .eq. 0) then
         write(6,*) 'RT_Base: uDen = ', uDen
         write(6,*) 'RT_Base: uLen = ', uLen
         write(6,*) 'RT_Base: uTime = ', uTime
         write(6,*) 'RT_Base: aExp = ', aExp
         write(6,*) 'RT_Base: Hubble = ', Hubble
         write(6,*) 'RT_Base: dtStep = ', dtStep
         write(6,*) 'RT_Base: uRate = ', uRate
         write(6,*) 'RT_Base: uColumn = ', uColumn
      endif
#endif
C
C  CMB temperature (in K)
C
      TemCMB = 2.728/aExp
C
C  Destruction of H- by CMB photons
C  NG 20091001: I do not recall where this is coming from, and it is
C  completely negligible for z<100, so ignore it.
C
C      facH2CMB = 1.1e-1*TemCMB**2.13*exp(-8823.0/TemCMB)
      facH2CMB = 0.0
C
C   dEin/dt = facCool * CoolFun * rho 
C   Only need to switch from proper time unit and from ergs to K
C
      facCool = uRate*(CFSTAR/1.381e-16)
C
C  P_i = facPhot * Sum_j sigma_j ang_j/aExp**3 
C  facPhot = c * CSSTAR * XISTEP / aExp**3 
C
      facPhot = 2.998e10*CSSTAR*XISTEP/aExp**3
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
C  dEa/dt = 1.735e-11 K cm^3/c angAvg(Lya)
C  Here t0 is in years.
C
#ifdef RT_LYMAN_ALPHA_HEATING
      facLya = 1.735e-11*angAvg(lrH1rB)
#else
      facLya = 0.0
#endif
C
C  Recombination case A/B:
C    tau = XH1*fOSC/hnu_eV*tau0
C    tau0 = c sigmaAB CSTAR nb / H / sqrt(Pi)
C
      if(IREC.eq.0 .and. Hubble.gt.1.0e-20) then
         facRecB = 2.998e10*CSSTAR*uDen/(uTime*Hubble)
      else
         facRecB = 0.0
      endif
C
C Set the source spectrum
C
      call frtSetSourceSpectrum(aExp,guv)
C
C  Normalize
C
      fac = uDen*aExp**3
      do lr=1,NRAD
         guvEff(lr) = fac*guv(lr)
      enddo
C
#ifdef RT_CHEMISTRY
#ifdef RT_LWBANDS

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(guvLW0,guvLWe,guvEff,lrLWtr)
      do lr=1,NRLW
         guvLWe(lr) = guvLW0(lr)*guvEff(lrLWtr)
      enddo
C     
C  Call Massimo's function
C         
      call  frtLW_SigmaH2(Hubble,dtStep,csLWMH)
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
      call frtInitGlobalFields

#ifdef RT_EXTERNAL_BACKGROUND
#if (RT_EXTERNAL_BACKGROUND != RT_BACKGROUND_SELFCONSISTENT)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(angLoc)
      do lr=1,NRAD
         angLoc(lr) = 0.0
      enddo

      call frtSetExternalBackground(aExp,angUni,angLW)

#endif /* RT_EXTERNAL_BACKGROUND != RT_BACKGROUND_SELFCONSISTENT */
#endif /* RT_EXTERNAL_BACKGROUND */

      call frtSetGlobalRadiationFields(rfAvg)

      return
      end
C
C
C
      subroutine frtSetGlobalRadiationFields(rfAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension rfAvg(*)
C
      include 'frt_inline.inc'
C
#ifdef RT_TRANSFER
      call frtSetGlobalRadiationFieldsTransfer(rfAvg)
#endif

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(angUni,angLoc,angAvg,angEff,tauLoc)
      do lr=1,NRAD
         angEff(lr) = angUni(lr) + angLoc(lr)*QFuni(tauLoc(lr))
         angAvg(lr) = angUni(lr) + angLoc(lr)
      enddo

#ifdef RT_OUTPUT
      call frtQueryId(id)
      if(id .eq. 0) then
         write(6,*) 'RT: properties of global radiation fields:'
         write(6,*) '  nu_H1: ', angUni(lrH1tr), angAvg(lrH1tr), 
     .        angEff(lrH1tr), facLX*guvEff(lrH1tr)
         write(6,*) '  nu_G1: ', angUni(lrG1tr), angAvg(lrG1tr), 
     .        angEff(lrG1tr), facLX*guvEff(lrG1tr)
         write(6,*) '  nu_G2: ', angUni(lrG2tr), angAvg(lrG2tr), 
     .        angEff(lrG2tr), facLX*guvEff(lrG2tr)
      endif
#endif

      return
      end
C
C
C
      subroutine frtUpdateTables(rfAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension rfAvg(*)

      call frtSetGlobalRadiationFields(rfAvg)
      call frtFillRadiationTables

      return
      end
C
C
C
      subroutine frtStepEnd(vol,parAvg,abcAvg)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      include 'frt_global.inc'
      dimension parAvg(*), abcAvg(*)
C
      dimension axi(NRAD), sxi(NRAD), xib(NRAD)
#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)
      dimension axiLW(NRLW), sxiLW(NRLW), xibLW(NRLW), anbLW(NRLW)
#endif
C
      common/RT_RT_UV/ duvin, axiH1tr, angH1tr
C
C  Diagnostics
C
#ifdef RT_OUTPUT
      call frtQueryId(id)
      if(id .eq. 0) then

         ftabHit = gInfo(infoTabH)
         ftabMiss = gInfo(infoTabM)

#ifdef RT_XRAYS
         ftabFill = 4.0*NOPT**3
#else
         ftabFill = 2.0*NOPT**3
#endif
         
         if(ftabHit .gt. 0.5) then
            write(6,*) 'RT_Tables: miss/hit ', 
     .           ftabMiss/ftabHit, ', fill/hit ', ftabFill/ftabHit
         endif

         if(gInfo(infoCall) .gt. gInfo(infoFail)) then
            write(6,*) 'RT_Cooling: avg # of steps ', 
     .           gInfo(infoNAvg)/(gInfo(infoCall)-gInfo(infoFail))
            write(6,*) 'RT_Cooling: max # of steps ', 
     .           gInfo(infoNMax)
         endif

         if(gInfo(infoFail) .gt. 0.0d0) then
            write(6,*) 'RT_Coolin: fail/call ', 
     .           gInfo(infoFail)/gInfo(infoCall)
         endif
      endif
#endif

#ifdef RT_EXTERNAL_BACKGROUND
#if (RT_EXTERNAL_BACKGROUND == RT_BACKGROUND_SELFCONSISTENT)
C
C  Unpack global quantities
C
      duvin = parAvg(1)
#ifdef RT_CHEMISTRY
      XH1_va = parAvg(2)
      XMH_va = parAvg(3)
#endif

      call frtFinalizeGlobalFields(vol)
C
C  Coefficient for Bremstrahlung: 
C     d\tilde{n}_xi/dt = 1/V e^{-h\nu/kT} 1/kT dE_tot/dt
C     gBrem(xi) = Sum dV rho dT/dt 1/T e^{-T_xi/T} / V
C     CRAD1 = \tilde{n}_b = comoving n_b
C
      CRAD1 = uDen*aExp**3
C
C  Coefficient for absorption:
C    c*CSSTAR*n_b/uTime = uColumn/facLX
C
      CRAD2 = uColumn/facLX
C
C  Advance the average radiation field:
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
         axi(lr) = axi(lr) + CRAD2*csLyS(lr)*gData(iavgCXH1)
#ifdef RT_DUST
         axi(lr) = axi(lr) + CRAD2*csDust(lr)*gData(iavgDust)
#endif
      enddo
 
#ifdef RT_TRANSFER
     
      call frtTransferComputeAverageAbs(axi,abcAvg)

#else

C
C Optically thin absortion
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(axi,CRAD2,lrmin,lrmax,csH1,csG1,csG2,gData)
      do lr=lrmin,lrmax
               
         axi(lr) = axi(lr) + CRAD2*(csH1(lr)*gData(iavgCXH1)+
     .        csG1(lr)*gData(iavgCXG1)+csG2(lr)*gData(iavgCXG2))
      enddo

#endif

#ifdef RT_CHEMISTRY
C
C  Absorption of Lyman-Wienner band photons (only included in the volume 
C  averaged equation to fix Olber's paradox).
C
      do lr=lrLWtr,lrH1tr-1
         axi(lr) = axi(lr) + csMHLW(lr)*XMH_va*CRAD2
      enddo
#endif
C
C  Frequency shift
C
      shift = Hubble*dtStep
      lshift = int(shift/XISTEP)
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(xib,xig,shift)
      do lr=1,NRAD
         xib(lr) = xig(lr) + shift
      enddo
C
C  Step 1: background due to uniform components
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(CRAD1,gBrem,sxi)
      do lr=1,NRAD
         sxi(lr) = CRAD1*gBrem(lr)
      enddo
C
C  Hydrogen and Helium recombination photons:
C  Recombination case A
C
      sxi(lrH1tr) = sxi(lrH1tr) + CRAD1*gData(iavgRAH2)/XISTEP
      sxi(lrG1tr) = sxi(lrG1tr) + CRAD1*gData(iavgRAG2)/XISTEP
      sxi(lrG2tr) = sxi(lrG2tr) + CRAD1*gData(iavgRAG3)/XISTEP
C
C  Recombination case B
C
      sxi(lrH1rB) = sxi(lrH1rB) + CRAD1*gData(iavgRBH2)/XISTEP
      sxi(lrG1rB) = sxi(lrG1rB) + CRAD1*gData(iavgRBG2)/XISTEP
      sxi(lrG2rB) = sxi(lrG2rB) + CRAD1*gData(iavgRBG3)/XISTEP
C
      call frtAdvanceBackground(lshift,xib,axi,sxi,angUni)
C
C  Step 2: background due to local sources
C
#ifdef RT_TRANSFER
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(duvin,guvEff,sxi)
      do lr=1,NRAD
         sxi(lr) = duvin*guvEff(lr)
      enddo
C
      call frtAdvanceBackground(lshift,xib,axi,sxi,angLoc)
#endif
C
C  Do Lyman-Werner bands separately
C
#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)

      fMH = XMH_va*CRAD2/CSSTAR
      fH1 = XH1_va*CRAD2/CSSTAR
      call frtLW_ComputeLineShifts(Hubble,fMH)
      call frtLW_ComputeAbs(Hubble,fH1,dtStep,0.0,axiLW)
      call frtLW_ComputeSource(Hubble,fH1,dtStep,guvLWe,sxiLW)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(anbLW,angLW,axiLW,duvin,sxiLW,dtStep)
      do lr=1,NRLW
         anbLW(lr) = angLW(lr)*axiLW(lr) + duvin*sxiLW(lr)*dtStep
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
         lbr = lbl + 1
         if(lbl .gt. NRLW) then
            angLW(lr) = 0.0
         else 
            anbl = anbLW(lbl)
            xibl = xigLW(lbl)
            if(lbr .gt. NRLW) then
               anbr = angH1tr
            else
               anbr = anbLW(lbr)
            endif
            xibr = xibl + XILWST
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

      axiH1tr = axi(lrH1tr)

#endif /* RT_EXTERNAL_BACKGROUND == RT_BACKGROUND_SELFCONSISTENT */
C
C  Total
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(angUni,angLoc,angAvg)
      do lr=1,NRAD
         angAvg(lr) = angUni(lr) + angLoc(lr)
      enddo

      angH1tr = angAvg(lrH1tr)

#ifdef RT_DEBUG
      call frtQueryId(id)
      if(id .eq. 0) then
         fac = 1.59e-17/aExp**3
         write(6,*) 'RT: <abs_HI> = ', axiH1tr, axiH1tr/csH1tr/CRAD2
         write(6,*) 'RT: <PIR_HI> = ', 
     .        angUni(lrH1tr)*csH1tr*facPhot, 
     .        angLoc(lrH1tr)*csH1tr*facPhot
         write(6,*) 'RT: <Jnu_HI> = ', 
     .        angUni(lrH1tr)*fac, angLoc(lrH1tr)*fac
         write(6,*) 'RT: <Jnu_LW> = ', 
     .        angUni(lrH1tr-1)*fac, angLoc(lrH1tr-1)*fac
      endif
#endif
#endif /* RT_EXTERNAL_BACKGROUND */

      return
      end
C
C
C
      subroutine frtAdvanceBackground(lshift,xib,axi,sxi,ang)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension xib(NRAD), axi(NRAD), sxi(NRAD), ang(NRAD)
C
      dimension anb(NRAD)
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,t12,p12,q12),
C$OMP+SHARED(axi,dtStep,anb,ang,sxi)
      do lr=1,NRAD

         t12 = axi(lr)*dtStep
         if(t12 .gt. 20.0) then
            p12 = 0.0
            q12 = 1.0/axi(lr)
         else if(t12 .lt. 0.001) then
            p12 = 1.0-t12
            q12 = dtStep
         else
            p12 = exp(-t12)
            q12 = (1.0-p12)/axi(lr)
         endif

         anb(lr) = ang(lr)*p12 + sxi(lr)*q12 

      enddo

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,lbl,lbr,anbl,xibl,anbr,xibr),
C$OMP+SHARED(lshift,ang,anb,xib,xig)
      do lr=1,NRAD
         lbl = lr + lshift
         lbr = lbl + 1
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
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(ang)
      do lr=1,NRAD
         ang(lr) = max(0.0,ang(lr))
      enddo

      return
      end
C
C
C
      function frtLogWithUnits(val,iDen,iLen,iTime)
      include 'frt_base.inc'
      include 'frt_tables.inc'

      if(uDen.gt.0.0 .and. uLen.gt.0.0 .and. uTime.gt.0.0 .and.
     .     val.gt.0.0) then
         frtLogWithUnits = log10(val) + iDen*log10(uDen) + 
     .        iLen*log10(uLen) + iTime*log10(uTime)
      else
         frtLogWithUnits = -1000
      endif

      return
      end
C
C
C     
      subroutine frtGetBackgroundPhotoRates(pRate)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension pRate(*)

      do i=1,iptabDim
         pRate(i) = pRate0(i)
      enddo

      return
      end
C
C
C
      function frtGetBinId(wlen)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      integer frtGetBinId

      if(wlen .lt. 1.0e-35) then
         frtGetBinId = 0
         return
      endif

      xi = log(911.75/wlen)

      if(xi.lt.xig(1) .or. xi.gt.xig(NRAD)) then
         lr = 0
      else
         lr = 1 + nint((xi-xig(1))/XISTEP)
         if(lr .gt. NRAD) lr = NRAD
      endif

      frtGetBinId = lr

      return
      end
C
C
C
      function frtGetBinWavelength(lr)
      include 'frt_base.inc'
      include 'frt_tables.inc'

      if(lr.ge.1 .and. lr.le.NRAD) then

         wlen = 911.75*exp(-xig(lr))

      else

         wlen = 0.0

      endif

      frtGetBinWavelength = wlen

      return
      end
C
C
C
      function frtGetBackgroundRadiationField(lr)
      include 'frt_base.inc'
      include 'frt_tables.inc'

      if(lr.ge.1 .and. lr.le.NRAD) then

         anxi = angAvg(lr)

      else

         anxi = 0.0

      endif

      frtGetBackgroundRadiationField = anxi

      return
      end
C
C
C     
      function frtGetRadiationField(lr,rf)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_TRANSFER
      include 'frt_transfer.inc'
#endif
      dimension rf(*)
      include 'frt_inline.inc'

      if(lr.ge.1 .and. lr.le.NRAD) then

#ifdef RT_EXTERNAL_BACKGROUND
#ifdef RT_TRANSFER
         tauH1l = rf(irfH1gf)
         tauG1l = rf(irfG1gf)
         tauG2l = rf(irfG2gf)
#else
         tauH1l = 0.0
         tauG1l = 0.0
         tauG2l = 0.0
#endif
         if(tauH1l .gt. 0.0) then
#ifdef RT_DUST
            csH1l = zdust*csDustn(lr) + csLySn(lr) + csH1n(lr)
#else
            csH1l = csLySn(lr) + csH1n(lr)
#endif
            taul = csH1l*tauH1l
         else
            taul = 0.0
         endif
         if(tauG1l .gt. 0.0) then
            taul = taul + csG1n(lr)*tauG1l
         endif
         if(tauG2l .gt. 0.0) then
            taul = taul + csG2n(lr)*tauG2l
         endif

         if(
     .        (tauH1l.lt.0 .and. lr.ge.lrmin .and. lr.lt.lrG1tr) .or.
     .        (tauG1l.lt.0 .and. lr.ge.lrG1tr .and. lr.lt.lrG2tr) .or.
     .        (tauG2l.lt.0 .and. lr.ge.lrG2tr)) then
            anxi = 0.0
         else
            anxi = angEff(lr)*QFuni(taul)
         endif
#else
         anxi = 0.0
#endif

#ifdef RT_TRANSFER
         tauH1l = rf(irfH1lf)
         tauG1l = rf(irfG1lf)
         tauG2l = rf(irfG2lf)

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
         if(tauG2l .gt. 0.0) then
            taul = taul + csG2n(lr)*tauG2l
         endif

         if(
     .        (tauH1l.lt.0 .and. lr.ge.lrmin .and. lr.lt.lrG1tr) .or.
     .        (tauG1l.lt.0 .and. lr.ge.lrG1tr .and. lr.lt.lrG2tr) .or.
     .        (tauG2l.lt.0 .and. lr.ge.lrG2tr)) then
         else
            anxi = anxi + facLX*rf(irfNorm)*guvEff(lr)*QFloc(taul)
         endif
#endif

      else

         anxi = 0.0

      endif

      frtGetRadiationField = anxi

      return
      end
#endif
      subroutine frtDummyBase(i)
      i = 0
      end
