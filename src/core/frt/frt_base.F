#include "frt_config.h"      
#ifdef RADIATIVE_TRANSFER
C
C
C
      subroutine frtInitRun2(y,tm,d2g,ch2,cl,fg,fq,ip,ir,io)
      include 'frt_f77.inc'
C
      Yp = y
      Tmin = tm
      D2Gin = d2g
      ClumpH2 = ch2
      CohLenH2 = cl
      fGal = fg
      fQSO = fq
      IPOP = ip
      IREC = ir
      IOUNIT = io
C
      DataPath = './drt'
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
      data IPOP /2/, IREC/0/, IOUNIT/81/
C
      if(Yp.lt.1.0e-10 .or. Yp.gt.0.99) then
         write(0,*) 'RT::InitRun: fatal error, ',
     .        'Yp must be between 1e-10 and 0.99.'
         call frtAbort
      endif
C
      call frtInterfaceCheck
C
#ifdef RT_OUTPUT
      call frtSetIONode(ioNode)
      if(ioNode) then
         open(unit=IOUNIT, file='rt-run.log')
         write(IOUNIT,*) 'Run:'
         write(IOUNIT,*) 'Yp:       ', Yp
         write(IOUNIT,*) 'Tmin:     ', Tmin
         write(IOUNIT,*) 'D2Gmin:   ', D2Gmin
         write(IOUNIT,*) 'ClumpH2:  ', ClumpH2
         write(IOUNIT,*) 'CohLenH2: ', CohLenH2
         write(IOUNIT,*) 'fGal:     ', fGal
         write(IOUNIT,*) 'fQSO:     ', fQSO
         write(IOUNIT,*) 'IPOP:     ', IPOP
         write(IOUNIT,*) 'IREC:     ', IREC
         write(IOUNIT,*) 'IOUNIT:   ', IOUNIT
         write(IOUNIT,*) 'DataPath: ', DataPath
         close(IOUNIT)
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

      ib = index(DataPath,' ') - 1
      call frtInitCF(IOUNIT,DataPath(1:ib)//'/drt_cf.LHO13.dat')

#ifdef RT_CHEMISTRY
#ifdef RT_LWBANDS
C
C  Call Massimo's routine
C
      call frtLW_InitRun()
      call frtLW_ReadData(DataPath(1:ib),guvLW0)
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
      subroutine frtStepBegin(ud,ul,ut,a,h,da)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      dimension guv(NRAD)
C
      include 'frt_inline.inc'
C
      uDen = ud
      uLen = ul
      uTime = ut
      aExp = a
      Hubble = h
      daExp = da
C
C  Additional units (and conversions)
C
      qExp3 = 1.0/aExp**3
      uRate = uDen*uTime
      uColumn = CSSTAR*uDen*uLen
C
#ifdef RT_OUTPUT
      if(ioNode) then
         open(unit=IOUNIT, file='rt-run.log', access='append')
         write(IOUNIT,*) 'Step:'
         write(IOUNIT,*) '  aExp:    ', aExp
         write(IOUNIT,*) '  uDen:    ', uDen
         write(IOUNIT,*) '  uLen:    ', uLen
         write(IOUNIT,*) '  uTime:   ', uTime
         write(IOUNIT,*) '  uRate:   ', uRate
         write(IOUNIT,*) '  uColumn: ', uColumn
         write(IOUNIT,*) '  Hubble:  ', Hubble
         write(IOUNIT,*) '  daExp:   ', daExp
         close(IOUNIT)
      endif
#endif /* RT_OUTPUT */

      if(aExp .lt. 1.0e-10) then
         write(0,*) 'RT::StepBegin: fatal error, ',
     .        'aExp has to be positive.'
         call frtAbort
      endif
C
C  Some constants
C
C ----------------------------------------------------
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
      call  frtLW_SigmaH2(daExp/aExp,csLWMH)
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
      varMin(ivarEin) = 1.0e-2*TemCMB
      varMin(ivarXMH) = 1.0e-20*XH
      varMin(ivarXH1) = 1.0e-15*XH
      varMin(ivarXH2) = 1.0e-15*XH
      varMin(ivarXG1) = 1.0e-15*XG
      varMin(ivarXG2) = 1.0e-15*XG
      varMin(ivarXG3) = 1.0e-15*XG
      varMin(ivarXHp) = 1.0e-25*XH
      varMin(ivarXHm) = 1.0e-25*XH

      varLow(ivarEin) = TemCMB
      varLow(ivarXMH) = 1.0e-8*XH
      varLow(ivarXH1) = 1.0e-6*XH
      varLow(ivarXH2) = 1.0e-6*XH
      varLow(ivarXG1) = 1.0e-6*XG
      varLow(ivarXG2) = 1.0e-6*XG
      varLow(ivarXG3) = 1.0e-6*XG
      varLow(ivarXHp) = 1.0e-8*XH
      varLow(ivarXHm) = 1.0e-8*XH
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

#ifdef RT_OUTPUT
      if(ioNode) then
         fac = 1.59e-17/aExp**3
         write(6,*) 'RT: <PIR_HI> = ', 
     .        angUni(lrH1tr)*csH1tr*facPhot, 
     .        angLoc(lrH1tr)*csH1tr*facPhot
         write(6,*) 'RT: <Jnu_HI> = ', 
     .        angUni(lrH1tr)*fac, angLoc(lrH1tr)*fac
         write(6,*) 'RT: <Jnu_LW> = ', 
     .        angUni(lrH1tr-1)*fac, angLoc(lrH1tr-1)*fac
      endif
#endif

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
      if(ioNode) then
         open(unit=IOUNIT, file='rt-run.log', access='append')
         write(IOUNIT,*) '  Properties of global radiation fields:'
         write(IOUNIT,*) '    nu_H1: ', angUni(lrH1tr), angAvg(lrH1tr),
     .        angEff(lrH1tr), facLX*guvEff(lrH1tr)
         write(IOUNIT,*) '    nu_G1: ', angUni(lrG1tr), angAvg(lrG1tr),
     .        angEff(lrG1tr), facLX*guvEff(lrG1tr)
         write(IOUNIT,*) '    nu_G2: ', angUni(lrG2tr), angAvg(lrG2tr),
     .        angEff(lrG2tr), facLX*guvEff(lrG2tr)
         close(IOUNIT)
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
      subroutine frtStepEnd(dt,vol,par)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
#ifdef RT_ADD_EXTERNAL_QSO_BACKGROUND
      include 'frt_spectrum.inc'
#endif
      include 'frt_global.inc'
      dimension par(*)
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
      if(ioNode) then

         ftabHit = gInfo(infoTabH)
         ftabMiss = gInfo(infoTabM)

#ifdef RT_XRAYS
         ftabFill = 4.0*NOPT**3
#else
         ftabFill = 2.0*NOPT**3
#endif
         
         open(unit=IOUNIT, file='rt-run.log', access='append')

         if(ftabHit .gt. 0.5) then
            write(IOUNIT,*) '  Tables: miss/hit ', 
     .           ftabMiss/ftabHit, ', fill/hit ', ftabFill/ftabHit
         endif

         if(gInfo(infoCall) .gt. gInfo(infoFail)) then
            write(IOUNIT,*) '  Cooling: avg # of steps ', 
     .           gInfo(infoNAvg)/(gInfo(infoCall)-gInfo(infoFail))
            write(IOUNIT,*) '  Cooling: max # of steps ', 
     .           gInfo(infoNMax)
         endif

         if(gInfo(infoFail) .gt. 0.0d0) then
            write(IOUNIT,*) '  Cooling: fail/call ', 
     .           gInfo(infoFail)/gInfo(infoCall)
         endif

         close(IOUNIT)

      endif
#endif

      axiH1tr = 0.0
C
C  Coefficient for absorption:
C    c*CSSTAR*n_b/uTime = uColumn/facLX
C
      CRAD2 = uColumn/facLX

#ifdef RT_EXTERNAL_BACKGROUND
#if (RT_EXTERNAL_BACKGROUND == RT_BACKGROUND_SELFCONSISTENT)
C
C  Unpack global quantities
C
#ifdef RT_CHEMISTRY
      XH1_va = par(2)
      XMH_va = par(3)
#endif
C
C  This must be so for photon number conservation
C
      duvin = par(1)

      call frtFinalizeGlobalFields(vol,dt)
C
C  Coefficient for Bremstrahlung: 
C     d\tilde{n}_xi/dt = 1/V e^{-h\nu/kT} 1/kT dE_tot/dt
C     gBrem(xi) = Sum dV rho dT/dt 1/T e^{-T_xi/T} / V
C     CRAD1 = \tilde{n}_b = comoving n_b
C
      CRAD1 = uDen*aExp**3
C
C  Advance the average radiation field:
C  Absorption
C
      do lr=1,NRAD
         axi(lr) = 0.0
      enddo
C
C Optically thin absortion from dust
C
#ifdef RT_DUST_ABSORBS_BACKGROUND
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(axi,CRAD2,csDust,lrmax,gData)
      do lr=1,lrmax
         axi(lr) = axi(lr) + CRAD2*csDust(lr)*gData(iavgDust)
      enddo
#endif

#ifdef RT_TRANSFER
     
      call frtTransferComputeAverageAbs(axi)

#else

C
C Optically thin absortion
C
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(axi,CRAD2,lrH1tr,lrmax,csH1,csG1,csG2,gData)
      do lr=lrH1tr,lrmax
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
      axiH1tr = axi(lrH1tr)
C
C  Frequency shift
C
      shift = Hubble*dt
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
C  Recombination case AaLumDen
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
C  Background from external (outside of the box) QSOs
C
#ifdef RT_ADD_EXTERNAL_QSO_BACKGROUND

      z = 1/(aExp+0.5*daExp) - 1
C
C  My fit to Hopkins et all 2007 (ApJ, 654, 731) comoving 
C  bolometric luminosity density from QSOs (in Lsun/Mpc^3)
C
      aLumDen = 10.0**7.7/((exp(1.2*(z-4)))**0.33+
     .     (3*exp(-2.7*z))**0.33)**3
      aKion = 3.9
C
C  Change units:
C    j = \int_0^\infty E_\nu (dn_\xi/dt) d\xi
C    j = A \int_0^\infty E_\nu g_\xi d\xi
C      = A * Kion * m_pc^2
C  where Kion = 3.9 is the bolometric correction from the ionizing to
C  the bolometric luminosity. Hence
C    (dn_\xi/dt) = (j/(Kion*m_pc^2)) g_\xi
C  or
C    (dn_\xi/dt)_{code} = uTime/uDen*(j/(Kion*m_pc^2)) g_\xi
C    q = uTime/uDen*(j/(Kion*m_pc^2)) = (uTime/1e18s)*(1cm^{-3}/uDen)*
C        1e18s*cm^3*(1Mpc/3.086e24cm)^3*(j)*3.9e33 erg/s/Mpc^3/
C        (Kion*1.674e-24*9e20 erg)
C
C  Transform to code units, being careful with not overfilling 
C  a single precision number.
C
      q = (uTime/1.0e18)*(1.0/uDen)*aLumDen*8.81e-20/aKion

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(q,guv_Q,sxi)
      do lr=1,NRAD
         sxi(lr) = sxi(lr) + q*guv_Q(lr)
      enddo
#endif /* RT_ADD_EXTERNAL_QSO_BACKGROUND */

      call frtAdvanceBackground(dt,lshift,xib,axi,sxi,angUni)
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
      call frtAdvanceBackground(dt,lshift,xib,axi,sxi,angLoc)
#endif
C
C  Do Lyman-Werner bands separately
C
#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)

      fMH = XMH_va*CRAD2/CSSTAR
      fH1 = XH1_va*CRAD2/CSSTAR
      call frtLW_ComputeLineShifts(Hubble,fMH)
      call frtLW_ComputeAbs(Hubble,fH1,dt,0.0,axiLW)
      call frtLW_ComputeSource(Hubble,fH1,dt,guvLWe,sxiLW)

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(lr),
C$OMP+SHARED(anbLW,angLW,axiLW,duvin,sxiLW,dt)
      do lr=1,NRLW
         anbLW(lr) = angLW(lr)*axiLW(lr) + duvin*sxiLW(lr)*dt
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

#ifdef RT_OUTPUT
      if(ioNode) then
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
      subroutine frtAdvanceBackground(dt,lshift,xib,axi,sxi,ang)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension xib(NRAD), axi(NRAD), sxi(NRAD), ang(NRAD)
C
      dimension anb(NRAD)
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,t12,p12,q12),
C$OMP+SHARED(axi,dt,anb,ang,sxi)
      do lr=1,NRAD

         t12 = axi(lr)*dt
         if(t12 .gt. 20.0) then
            p12 = 0.0
            q12 = 1.0/axi(lr)
         else if(t12 .lt. 0.001) then
            p12 = 1.0-t12
            q12 = dt
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

#if defined(RT_OUTPUT) && defined(RT_DEBUG)
      if(ioNode) then
         lr = lrH1tr
         write(6,*) 'RT::Advance Background:', lr, ang(lr), axi(lr),
     .        axi(lr)*dt, sxi(lr)
      endif
#endif

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
C  Return cooling & heating rates of one resolution element
C
      subroutine frtGetCoolingRate(var,rawrf,rateCool,rateHeat)
      include 'frt_base.inc'
      dimension var(*), rawrf(*)

      call frtCoolOff(var,rawrf,-1.0,info)

      rateCool = 1.381e-16*var(ivarCool)/uTime
      rateHeat = 1.381e-16*var(ivarHeat)/uTime

      return
      end
C
C
C     
      subroutine frtGetPhotoRates(var,rawrf,pRate)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension var(*), rawrf(*)
      dimension pRate(*)
#ifdef RT_TRANSFER
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
#else
      dimension rf(1)
#endif
      dimension itab(2)

#ifdef RT_TRANSFER
#ifdef RT_TRANSFER_FLUX_CONSERVING
      var(ivarDelN) = uColumn*var(ivarRho)*var(ivarDelX)
#endif
      call frtTransferPackRadiationField(var,rawrf,rf)
#endif

      call frtPhotoRates(var,rf,itab,pRate)

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
      subroutine frtGetRadiationField(var,rawrf,nout,wlen,anxi)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension var(*), rawrf(*)
      dimension wlen(nout), anxi(nout)
#ifdef RT_TRANSFER
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
#else
      dimension rf(1)
#endif
      integer frtBinId
      external frtBinId
      include 'frt_inline.inc'

#ifdef RT_TRANSFER
#ifdef RT_TRANSFER_FLUX_CONSERVING
      var(ivarDelN) = uColumn*var(ivarRho)*var(ivarDelX)
#endif
      call frtTransferPackRadiationField(var,rawrf,rf)
#endif

      do iout=1,nout

         lr = frtBinId(wlen(iout))

         if(lr.ge.1 .and. lr.le.NRAD) then

            call frtGetRFHelper(lr,rf,uLoc,uGlob)
            anxi(iout) = uLoc*rf(irfNorm) + uGlob*1.0

         else

            anxi(iout) = 0.0

         endif
 
      enddo

      return
      end
C
C
C     
      subroutine frtGetRFHelper(lr,rf,uLoc,uGlob)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension rf(*)
#ifdef RT_TRANSFER
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
#endif
      include 'frt_inline.inc'

#ifdef RT_EXTERNAL_BACKGROUND
#ifdef RT_TRANSFER
      tauH1l = rf(irfH1gf)
      tauG1l = rf(irfG1gf)
      tauG2l = rf(irfG2gf)
#ifdef RT_UV
      tauUVl = rf(irfUVgf)
#endif
#else
      tauH1l = 0.0
      tauG1l = 0.0
      tauG2l = 0.0
#ifdef RT_UV
      tauUVl = 0.0
#endif
#endif

#ifdef RT_UV
      taul = csUVn(lr)*tauUVl
#else
      taul = 0.0
#endif
      if(tauH1l .lt. VLARGE) then
         taul = taul + csH1n(lr)*tauH1l
      endif
      if(tauG1l .lt. VLARGE) then
         taul = taul + csG1n(lr)*tauG1l
      endif
      if(tauG2l .lt. VLARGE) then
         taul = taul + csG2n(lr)*tauG2l
      endif

      if(
#ifdef RT_UV
     .     (tauUVl.gt.VLARGE .and. lr.ge.1 .and. lr.lt.lrH1tr) .or.
#endif
     .     (tauH1l.gt.VLARGE .and. lr.ge.lrH1tr .and. lr.lt.lrG1tr)
     .     .or.
     .     (tauG1l.gt.VLARGE .and. lr.ge.lrG1tr .and. lr.lt.lrG2tr)
     .     .or. (tauG2l.gt.VLARGE .and. lr.ge.lrG2tr)) then
         uGlob = 0.0
      else
C
C  Turn to physical units
C
         uGlob = angEff(lr)*QFuni(taul)*qExp3
      endif
#else
      uGlob = 0.0
#endif

#ifdef RT_TRANSFER
      tauH1l = rf(irfH1lf)
      tauG1l = rf(irfG1lf)
      tauG2l = rf(irfG2lf)
#ifdef RT_UV
      tauUVl = rf(irfUVlf)
#endif

#ifdef RT_UV
      taul = csUVn(lr)*tauUVl
#else
      taul = 0.0
#endif
      if(tauH1l .lt. VLARGE) then
         taul = taul + csH1n(lr)*tauH1l
      endif
      if(tauG1l .lt. VLARGE) then
         taul = taul + csG1n(lr)*tauG1l
      endif
      if(tauG2l .lt. VLARGE) then
         taul = taul + csG2n(lr)*tauG2l
      endif

      if(
#ifdef RT_UV
     .     (tauUVl.gt.VLARGE .and. lr.ge.1 .and. lr.lt.lrH1tr) .or.
#endif
     .     (tauH1l.gt.VLARGE .and. lr.ge.lrH1tr .and. lr.lt.lrG1tr)
     .     .or.
     .     (tauG1l.gt.VLARGE .and. lr.ge.lrG1tr .and. lr.lt.lrG2tr)
     .     .or. (tauG2l.gt.VLARGE .and. lr.ge.lrG2tr)) then
         uLoc = 0.0
      else
C
C  Turn to physical units
C
         uLoc = facLX*guvEff(lr)*QFloc(taul)*qExp3
      endif
#endif

      return
      end
C
C
C
      subroutine frtGetBackgroundRadiationField(nout,wlen,anxi)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension wlen(nout), anxi(nout)
      integer frtBinId
      external frtBinId

      do iout=1,nout

         lr = frtBinId(wlen(iout))

         if(lr.ge.1 .and. lr.le.NRAD) then
            anxi(iout) = angAvg(lr)*qExp3
         else
            anxi(iout) = 0.0
         endif

      enddo

      return
      end
C
C
C
      function frtBinId(wlen)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      integer frtBinId

      if(wlen .lt. VSMALL) then
         frtBinId = 0
         return
      endif

      xi = log(911.75/wlen)

      if(xi.lt.xig(1) .or. xi.gt.xig(NRAD)) then
         lr = 0
      else
         lr = 1 + nint((xi-xig(1))/XISTEP)
         if(lr .gt. NRAD) lr = NRAD
      endif

      frtBinId = lr

      return
      end
C
C  This is just for adding radiation pressure;
C  rf2Prad is returned in K, 
C  rho2mfp is in code units
C     
      subroutine frtGetRPFactors(rf2Prad,rho2abc)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      integer frtBinId
      external frtBinId

#ifdef RT_TRANSFER
C
C  n_xi is in comoving units, hence aExp
C
      rf2Prad = Txi(lrUVtr)*facLX*guvEff(lrUVtr)*qExp3
      rho2abc = uColumn*csUV(lrUVtr)
#else
      rf2Prad = 0.0
      rho2abc = 0.0
#endif

      return
      end
C
C
C
#include "frt_f77.h"
#define CHECK(v1,v2) if(v1 .ne. v2) ierr = v1
      subroutine frtInterfaceCheck
      include 'frt_base.inc'

      ierr = 0

      CHECK(ivarEin,FRT_Ein)
      CHECK(ivarXH1,FRT_XHI)
      CHECK(ivarXH2,FRT_XHII)
      CHECK(ivarXG1,FRT_XHeI)
      CHECK(ivarXG2,FRT_XHeII)
      CHECK(ivarXG3,FRT_XHeIII)
      CHECK(ivarXMH,FRT_XH2)
      CHECK(ivarXHp,FRT_XH2p)
      CHECK(ivarXHm,FRT_XHm)
      CHECK(ivarRho,FRT_Density)
      CHECK(ivarZsol,FRT_Metallicity)
      CHECK(ivarSobL,FRT_SobolevLength)
      CHECK(ivarNumF,FRT_NumericalDiffusionFactor)
      CHECK(ivarRF0l,FRT_OTRadiationFieldLocal)
      CHECK(ivarRF0g,FRT_OTRadiationFieldGlobal)
      CHECK(ivarDTex,FRT_ExternalHeatingRate)
      CHECK(ivarBias,FRT_SourceBias)
      CHECK(ivarCSup,FRT_CoolingSuppressionFactor)
      CHECK(ivarDelX,FRT_ResolutionElementSize)
      CHECK(ivarDelV,FRT_ResolutionElementVolume)
      CHECK(ivarLTE,FRT_LTEFlag)
      CHECK(ivarDeb,FRT_Debug)
      CHECK(ivarGam,FRT_Gamma)
      CHECK(ivarDGas,FRT_DustToGasRatio)
      CHECK(ivarDim,FRT_DIM)

      if(ierr .ne. 0) then
         write(0,*) 'RT: Internal Bug (interface error) ', ierr
         stop
      endif

      return
      end
#endif
      subroutine frtDummyBase(i)
      i = 0
      end
