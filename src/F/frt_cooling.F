#include "../config.h"
#ifdef RADIATIVE_TRANSFER
C
C  Return cooling & heating rates of one resolution element
C
      subroutine frtCoolingRate(par,rawRF0,rawRF1,y0,rateCool,rateHeat)
      include 'frt_base.inc'
      dimension par(*), rawRF0(*), rawRF1(*), y0(*)

      call frtCoolOff(par,rawRF0,rawRF1,0.0,y0,info)

      rateCool = 1.381e-16*par(iparCool)/uTime
      rateHeat = 1.381e-16*par(iparHeat)/uTime

      return
      end
C
C  Follow cooling evolution of one resolution element
C
      subroutine frtCoolOff(par,rawRF0,rawRF1,tFin,y0,info)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_TRANSFER
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
#endif
      dimension par(*), rawRF0(*), rawRF1(*), y0(*)
C
C  This function solves the equation
C   
C  dEin/dt = fEinSh*delta(t) + fEinAd + fEinHeat - fEinCool from t=0 to t=tFin
C
      parameter(errlim = 1.0e-2)
#ifdef RT_8SPECIES
      parameter(nEqs = ieqDim)
#else
#ifdef RT_CHEMISTRY
      parameter(nEqs = 7)
#else
      parameter(nEqs = 6)
#endif
#endif
C
      dimension yf(ieqDim), rs(ieqDim)
      dimension w0(ieqDim), a0(ieqDim)
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
      dimension avg(iavgDim), avg0(iavgMax1), avg1(iavgMax1)
#else
      dimension avg(1), avg0(1), avg1(1)
#endif
C
      dimension itab(2), pRate(iptabDim)
#ifdef RT_TRANSFER
#ifdef RT_TRANSFER_FLUX_CONSERVING
      dimension rf(irfDim)
#else
      dimension rf(irfMax)
#endif
#else
      dimension rf(1)
#endif
C
#if defined(RT_TEST) && (RT_TEST==1 || RT_TEST==11)
      y0save = y0(ieqEin)
#endif

      y0(ieqEin) = y0(ieqEin) + par(iparDTsh)*tFin

#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
      do j=1,iavgMax1
         avg(j) = 0.0
      enddo
#endif
C
C  Sub-cell bias factor & cooling suppressing factor:
C  if not set, set them both to 1
C
      if(par(iparBias) .lt. 2.0e-38) par(iparBias) = 1.0
      if(par(iparDFac) .lt. 2.0e-38) par(iparDFac) = 1.0
C
C  Metallicity-dependent cooling: parameters are saved in par(*) array.
C  par(iparSDp[1-3]) - parameters for the Sutherland & Dopita 93, ApJS,
C                  88, 253 cooling function as parametrized by Gnedin 
C  par(iparDGas) - dust-to-gas ratio. We assume it is proportional to 
C                  metallicity in solar units (more precisely, it is 
C                  10*SMC value is SMC metallicity is 0.1 solar).
C  Assignment order matters here, since par(iparDGas) is aliased.
C
      par(iparSDp1) = par(iparZsol)
      par(iparSDp2) = 2.7*exp(-(par(iparZsol)/0.316)**3)
      par(iparSDp3) = 0.04*par(iparZsol)**2
      par(iparDGas) = max(D2Gmin,frtDustToGas(par(iparZsol)))

      par(iparRecB) = facRecB*par(iparRhoB)

      itab(1) = 0
      itab(2) = 0

#ifdef RT_CHEMISTRY
      par(iparSobL) = uColumn*par(iparRhoB)*par(iparSobL)
C
C  Correct for numerical diffusion, if needed
C
      if(par(iparNumF) .gt. 1.0e-10) then
         XH1ll = max(yLimit(ieqXH1),min(XH,y0(ieqXH1)))
         XMHll = max(yLimit(ieqXMH),min(0.5*XH,y0(ieqXMH)))
         cdDUl = par(iparDGas)*(XH1ll+2*XMHll)*par(iparSobL)
         sFact1 = exp(-0.3*d_sDUST*cdDUl)
         fac = (1.0-sFact1)*par(iparNumF)
         y0(ieqEin) = y0(ieqEin)/(1+fac*y0(ieqEin)/30.0)
         y0(ieqXH1) = y0(ieqXH1)/(1+fac/max(1.0e-6,y0(ieqXH1)/XH))
         y0(ieqXH2) = y0(ieqXH2)/(1+fac/max(1.0e-6,y0(ieqXH2)/XH))
      endif
#endif

      call frtAdjust(y0,nEqs)

#ifdef RT_DEBUG
      ideb1 = nint(par(iparDeb))/10
      ideb2 = mod(nint(par(iparDeb)),10)
      if(ideb1 .eq. 1) then
         open(unit=81, file='rt-debug.log', access='append')
         write(81,*) 'In cell-level debug: ', uDen*par(iparRhoB)
         write(6,*) 'In cell-level debug: ', uDen*par(iparRhoB)
      endif
#endif

#ifdef RT_TRANSFER
#ifdef RT_TRANSFER_FLUX_CONSERVING
      par(iparDelN) = uColumn*par(iparRhoB)*par(iparCell)
#endif
      call frtTransferPackRadiationField(par,y0,rawRF0,rawRF1,rf)
#endif

#ifndef RT_VARIABLE_PRATES
      call frtGetPhotoRates(par,rf,itab,y0,pRate)
#endif
C
C  Cooling rate mode
C
      if(tFin .lt. 1.0e-30) then
         call frtGetRates(par,rf,itab,y0,w0,a0,avg0,pRate)
         par(iparCool) = a0(ieqEin)*y0(ieqEin)
         par(iparHeat) = w0(ieqEin)
         return
      endif
C
C  Evolution mode
C
      dt = tFin
      t = 0.0
      iLast = 1
      nStep = 0

      info = 9999

#ifdef RT_DEBUG
      if(ideb1 .eq. 1) then
#ifdef RT_VARIABLE_PRATES
         call frtGetPhotoRates(par,rf,itab,y0,pRate)
#endif
         write(81,*) 'Photo rates: ', 
     .        pRate(iptabPiH1), pRate(iptabPiG1), pRate(iptabPiG2),
     .        pRate(iptabPhH1), pRate(iptabPhG1), pRate(iptabPhG2)
#ifdef RT_CHEMISTRY
         write(81,*) 'Chemistry rates:',
     .        (pRate(j),j=iptabCi27,iptabCiLW)
#endif
         call frtDebugPrint1(81,nStep,t,dt,tFin,y0,1) 
      endif
#endif

 10   nStep = nStep + 1

      call frtGetRates(par,rf,itab,y0,w0,a0,avg0,pRate)
C
C  ODE Solver
C
      nFail = 0

 20   continue

      if(nStep .lt. 300) then
         call frtCoolOffStepQS1(nEqs,dt,par,rf,itab,
     .        y0,w0,a0,yf,rs,d1,avg1,pRate)
      else
         call frtCoolOffStepQSN(nEqs,dt,par,rf,itab,
     .        y0,w0,a0,yf,rs,d1,avg1,pRate)
      endif

      ieqmax = 1
      errmax = abs(rs(1))/max(yf(1),ySmall(1))
      do j=2,nEqs
         err = abs(rs(j))/max(yf(j),ySmall(j))
         if(err .gt. errmax) then
            errmax = err
            ieqmax = j
         endif
      enddo
      errmax = errmax/errlim

      if(errmax .gt. 1.0) then

         iLast = 0
         nFail = nFail + 1
         if(nFail .gt. 20) then
            info = -1
         else
            dt = 0.9*dt/(2+errmax)**0.5
            goto 20
         endif

      else

#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
         do j=1,iavgMax1
            avg(j) = avg(j) + dt*(d1*avg0(j)+(1-d1)*avg1(j))
         enddo
#endif

         t = t + dt

         dt = 0.9*dt/(0.01+errmax)**0.3
 
         if(iLast.eq.0 .and. t+dt.gt.tFin) then
            dt = tFin - t
            iLast = 1
         endif
C
C  In general, schemes are not conserving
C
         call frtAdjust(yf,nEqs)
         do j=1,nEqs
            y0(j) = yf(j)
         enddo

#ifdef RT_DEBUG
         if(ideb1 .eq. 1) then
            call frtDebugPrint1(81,nStep,t,dt,tFin,yf,ieqmax) 
         endif
#endif

      endif

      if(tFin-t.gt.0.5*dt .and. nStep.lt.info) goto 10

      if(info .ge. 0) then
         info = nStep
      endif

#ifdef RT_DEBUG
      if(ideb1 .eq. 1) then
         close(81)
         if(ideb2 .gt. 0) call frtAbort
      endif
#endif

#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
      avg(iavgBrem) = avg(iavgBrem)*par(iparRhoB)

      do j=1,iavgMax1
         avg(j) = avg(j)*par(iparRhoB)
      enddo

      avg(iavgNorm) = tFin

#if defined(RT_EXACT_EOS) && defined(RT_HIGH_DENSITY)
      avg(iavgTemp) = frtTem(y0(ieqEin),y0,0)
#else
      avg(iavgTemp) = frtTem(y0(ieqEin),y0)
#endif

      call frtIncrementGlobalFields(avg,par(iparVolB),itab,info)
#endif

#if defined(RT_TRANSFER) && defined(RT_VARIABLE_RFIELD)
      call frtTransferUnPackRadiationField(y0,rawRF0,rawRF1,rf)
#endif

#if defined(RT_TEST) && (RT_TEST==1 || RT_TEST==11)
      y0(ieqEin) = y0save
#endif

      return
      end
C
C ---------------------------------------------
C
C  Internal routines of the cooling block
C
C ---------------------------------------------
C
      subroutine frtAdjust(y,nEqs)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension y(*)

      do i=1,nEqs
         y(i) = max(y(i),yLimit(i))
      enddo

#ifdef RT_CHEMISTRY
      fact = XH/(y(ieqXH1)+y(ieqXH2)+2*y(ieqXMH))
      y(ieqXMH) = fact*y(ieqXMH)
#else
      fact = XH/(y(ieqXH1)+y(ieqXH2))
#endif
      y(ieqXH1) = fact*y(ieqXH1)
      y(ieqXH2) = fact*y(ieqXH2)

      fact = XG/(y(ieqXG1)+y(ieqXG2)+y(ieqXG3))
      y(ieqXG1) = fact*y(ieqXG1)
      y(ieqXG2) = fact*y(ieqXG2)
      y(ieqXG3) = fact*y(ieqXG3)

      if(Tmin.gt.0.0 .and. y(ieqEin).lt.10*Tmin) then
#if defined(RT_EXACT_EOS) && defined(RT_HIGH_DENSITY)
         y(ieqEin) = max(y(ieqEin),frtEin(Tmin,y,0))
#else
         y(ieqEin) = max(y(ieqEin),frtEin(Tmin,y))
#endif
      endif

      return
      end
C
C
C
      subroutine frtGetRates(par,rf,itab,
     .     y,w,a,avg,pRate)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      dimension par(*), rf(*), itab(*)
      dimension y(*), w(*), a(*)
      dimension avg(*), pRate(*)
C
      dimension aRate(irateAtomN)
#ifdef RT_CHEMISTRY
      dimension cRate(irateChemN)
#ifdef RT_HIGH_DENSITY
      dimension hRate(irateHighN)
#endif
#endif
      double precision BH1, BG1, BG2, AH2, AG2, AG3, dfac
C      
      den = uDen*par(iparRhoB)

      Einll = max(yLimit(ieqEin),y(ieqEin))
      XH1ll = max(yLimit(ieqXH1),min(XH,y(ieqXH1)))
      XH2ll = max(yLimit(ieqXH2),min(XH,y(ieqXH2)))
      XG1ll = max(yLimit(ieqXG1),min(XG,y(ieqXG1)))
      XG2ll = max(yLimit(ieqXG2),min(XG,y(ieqXG2)))
      XG3ll = max(yLimit(ieqXG3),min(XG,y(ieqXG3)))
#ifdef RT_CHEMISTRY
      XMHll = max(yLimit(ieqXMH),min(0.5*XH,y(ieqXMH)))
#ifdef RT_8SPECIES
      XHpll = max(yLimit(ieqXHp),min(0.5*XH,y(ieqXHp)))
      XHmll = max(yLimit(ieqXHm),min(0.5*XH,y(ieqXHm)))
#endif
#endif

      Xfell = XH2ll + XG2ll + 2*XG3ll  !! ignore contributions of H2+ and H-

      de = Xfell*den
#if defined(RT_EXACT_EOS) && defined(RT_HIGH_DENSITY)
      iLTE = 0
      tem = max(1.0e-3,Tmin,frtTem(Einll,y,iLTE))
#else
      tem = max(1.0e-3,Tmin,frtTem(Einll,y))
#endif
      altem = log(max(acTmin,tem))
      il = int((altem-acLTmin)/acLTstp*0.99999) + 1
      if(il .lt.    1) il = 1
      if(il .ge. NTEM) il = NTEM-1
      iu = il + 1
      ql = max(0.0,min(1.0,(acLT(iu)-altem)/acLTstp))
      qu = 1.0 - ql
C
      do j=1,irateAtomN
         aRate(j) = ql*RateAtom(j,il) + qu*RateAtom(j,iu)
      enddo

#ifdef RT_VARIABLE_PRATES
      call frtGetPhotoRates(par,rf,itab,y,pRate)
#endif
C
#ifdef RT_CHEMISTRY
C
C  H2 dust-shielding model
C
      Xdust = par(iparDGas)*(XH1ll+2*XMHll)*ClumpH2
*      Xdust = par(iparDGas)*(XH1ll+2*XMHll*ClumpH2)
      cdDUl = par(iparDGas)*(XH1ll+2*XMHll)*par(iparSobL)
      sFact1 = exp(-d_sDUST*cdDUl)
#ifdef RT_UV
      sFactI = exp(-d_sDUST*cdDUl)
#else
      sFactI = exp(-d_sDUST*cdDUl)
#endif
C
C  H2 self-shielding model
C
      cdMHl = XMHll*CSSTAR*3.086e18*den*CohLenH2
      sFact2 = max(1.0,cdMHl/(CSSTAR*1e14))**(-0.75)

*      x = cdMHl/(CSSTAR*5e14)
*      sFact2 = 0.965/(1+x)**2 + 
*     .     0.035/(1+x)**0.5*exp(-8.5e-4*(1+x)**0.5)

#else  /* RT_CHEMISTRY */
      sFact1 = 1.0
      sFact2 = 1.0
#endif /* RT_CHEMISTRY */

C
#ifdef RT_XLF_BUG_FIX1
C
C *****************************
C
C  This is to safeguard against that weird XLF bug
C
      wa = 0.0
      wb = 0.0
C
C *****************************
C
#endif
C
C  Recombination rates
C
      if(IREC .eq. 1) then
         RiH2l = aRate(irateRAiH2)
         RcH2l = aRate(irateRAcH2)
         RiG2l = aRate(irateRAiG2)
         RcG2l = aRate(irateRAcG2)
         RiG3l = aRate(irateRAiG3)
         RcG3l = aRate(irateRAcG3)
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
         avg(iavgRAH2) = uRate*XH2ll*de*aRate(irateRAiH2)
         avg(iavgRAG2) = uRate*XG2ll*de*aRate(irateRAiG2)
         avg(iavgRAG3) = uRate*XG3ll*de*aRate(irateRAiG3)
         avg(iavgRBH2) = 0.0
         avg(iavgRBG2) = 0.0
         avg(iavgRBG3) = 0.0
#endif
      else if(IREC .eq. 2) then
         RiH2l = aRate(irateRBiH2)
         RcH2l = aRate(irateRBcH2)
         RiG2l = aRate(irateRBiG2)
         RcG2l = aRate(irateRBcG2)
         RiG3l = aRate(irateRBiG3)
         RcG3l = aRate(irateRBcG3)
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
         avg(iavgRAH2) = 0.0
         avg(iavgRAG2) = 0.0
         avg(iavgRAG3) = 0.0
         avg(iavgRBH2) = uRate*XH2ll*de*aRate(irateRBiH2)
         avg(iavgRBG2) = uRate*XG2ll*de*aRate(irateRBiG2)
         avg(iavgRBG3) = uRate*XG3ll*de*aRate(irateRBiG3)
#endif
      else
         call frtGetRecCaseH2(XH1ll,par(iparRecB),wa,wb)
         RiH2l = wa*aRate(irateRAiH2) + wb*aRate(irateRBiH2)
         RcH2l = wa*aRate(irateRAcH2) + wb*aRate(irateRBcH2)
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
         avg(iavgRAH2) = uRate*XH2ll*de*wa*aRate(irateRAiH2)
         avg(iavgRBH2) = uRate*XH2ll*de*wb*aRate(irateRBiH2)
#endif

         call frtGetRecCaseG2(XG1ll,par(iparRecB),wa,wb)
         RiG2l = wa*aRate(irateRAiG2) + wb*aRate(irateRBiG2)
         RcG2l = wa*aRate(irateRAcG2) + wb*aRate(irateRBcG2)
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
         avg(iavgRAG2) = uRate*XG2ll*de*wa*aRate(irateRAiG2)
         avg(iavgRBG2) = uRate*XG2ll*de*wb*aRate(irateRBiG2)
#endif

         call frtGetRecCaseG3(XG2ll,par(iparRecB),wa,wb)
         RiG3l = wa*aRate(irateRAiG3) + wb*aRate(irateRBiG3)
         RcG3l = wa*aRate(irateRAcG3) + wb*aRate(irateRBcG3)
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
         avg(iavgRAG3) = uRate*XG3ll*de*wa*aRate(irateRAiG3)
         avg(iavgRBG3) = uRate*XG3ll*de*wb*aRate(irateRBiG3)
#endif
      endif
C
C  H ionization rates
C
      BH1 = aRate(irateCIiH1)*de + sFactI*pRate(iptabPiH1)
      AH2 = RiH2l*de
C
C  He ionization rates
C
      BG1 = aRate(irateCIiG1)*de + sFactI*pRate(iptabPiG1)
      BG2 = aRate(irateCIiG2)*de + sFactI*pRate(iptabPiG2)
      AG2 = (RiG2l+aRate(irateDRiG2))*de
      AG3 = RiG3l*de
C
C  Hydrogen rates can be modified by H2 reactions below
C
      a(ieqXH1) = uTime*BH1
      w(ieqXH1) = uTime*XH2ll*AH2

      a(ieqXH2) = uTime*AH2
      w(ieqXH2) = uTime*XH1ll*BH1
C
C  Helium rates are final, we do not include He moleculae
C
      a(ieqXG1) = uTime*BG1
      w(ieqXG1) = uTime*XG2ll*AG2

      a(ieqXG2) = uTime*(AG2+BG2)
      w(ieqXG2) = uTime*(XG1ll*BG1+XG3ll*AG3)

      a(ieqXG3) = uTime*AG3
      w(ieqXG3) = uTime*XG2ll*BG2
C
C  Cooling rates
C
      CF_PPl = Xfell*(
     .     aRate(irateCIcH1)*XH1ll + aRate(irateCIcG1)*XG1ll + 
     .     aRate(irateCIcG2)*XG2ll +
     .     RcH2l*XH2ll + RcG2l*XG2ll + RcG3l*XG3ll +
     .     aRate(irateDRcG2)*XG2ll +
     .     aRate(irateLCcH1)*XH1ll + aRate(irateLCcG2)*XG2ll)

      CF_FFl = facFF*sqrt(tem)*Xfell*(XH2ll+XG2ll+4*XG3ll)
C
C  Return parameters
C
#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
      avg(iavgCXH1) = par(iparCFH1)*XH1ll
      avg(iavgCXG1) = par(iparCFG1)*XG1ll
      avg(iavgCXG2) = par(iparCFG2)*XG2ll
      avg(iavgBrem) = facCool*CF_FFl
#endif

      aEinll = facCool*par(iparRhoB)*(CF_PPl+CF_FFl) !! Atomic cooling and Bremstrahlung
      wEinll = uTime*sFact1*(XH1ll*(facLya+pRate(iptabPhH1))+  !! PhotoHeating
     .     XG1ll*pRate(iptabPhG1)+XG2ll*pRate(iptabPhG2))

      aEinl2 = facCC*Xfell*tem  !! Compton cooling
      wEinl2 = facCC*Xfell*TemCMB  !! Compton heating

#ifdef RT_CHEMISTRY
C
C ***********************************************
C
C   Molecular hydrogen
C
C ***********************************************
C
      do j=irateChemC+1,irateChemN
         cRate(j) = ql*RateChem(j,il) + qu*RateChem(j,iu)
      enddo
#ifdef RT_HIGH_DENSITY
      do j=irateHighC+1,irateHighN
         hRate(j) = ql*RateHigh(j,il) + qu*RateHigh(j,iu)
      enddo
#endif

      do j=1,irateChemC
         cRate(j) = den*(ql*RateChem(j,il)+qu*RateChem(j,iu))
      enddo
#ifdef RT_HIGH_DENSITY
      do j=1,irateHighC
         hRate(j) = den*(ql*RateHigh(j,il)+qu*RateHigh(j,iu))
      enddo
#endif
C
C  H-minus rates
C
      aXHml = (cRate(irateMHi02)+cRate(irateMHi15))*XH1ll +
     .     (cRate(irateMHi05)+cRate(irateMHi16))*XH2ll + 
     .     cRate(irateMHi14)*Xfell + cRate(irateMHi28)*XG2ll +
     .     cRate(irateMHi29)*XG1ll + facH2CMB +
     .     sFact1*pRate(iptabCi27)

      wXHml = cRate(irateMHi01)*XH1ll*Xfell +
     .     cRate(irateMHi23)*XMHll*Xfell
#ifdef RT_8SPECIES
C     
C  H-minus non-equilibrium rates
C
      a(ieqXHm) = uTime*aXHml
      w(ieqXHm) = uTime*wXHml
#else
C
C  H-minus equilibrium
C
      XHmll = min(0.5*XH1ll,max(yLimit(ieqXHm),
     .     wXHml/(1.0e-35+aXHml)))
#endif
C
C  H2-plus rates
C
      pRate29 = sFact1*sFact2*pRate(iptabCi29)

      aXHpl = cRate(irateMHi04)*XH1ll + cRate(irateMHi06)*Xfell + 
     .     sFact1*(pRate(iptabCi28)+pRate(iptabCi30))

      wXHpl = cRate(irateMHi03)*XH1ll*XH2ll +
     .     cRate(irateMHi07)*XMHll*XH2ll + 
     .     cRate(irateMHi16)*XH2ll*XHmll + 
     .     cRate(irateMHi25)*XMHll*XG2ll +
     .     pRate29*XMHll
C
#ifdef RT_8SPECIES
C
C  H2-plus non-equilibrium rates
C
      a(ieqXHp) = uTime*aXHpl
      w(ieqXHp) = uTime*wXHpl
#else
C
C  H2-plus equilibrium
C
      XHpll = min(0.5*XH1ll,max(yLimit(ieqXHp),
     .     wXHpl/(1.0e-35+aXHpl)))
#endif
C
C  H2 ortho-to-para ratio
C
#ifdef RT_EXACT_EOS
#ifdef RT_HIGH_DENSITY
      rOP = frtH2_op(tem,iLTE)
#else
      rOP = frtH2_op(tem)
#endif /* RT_HIGH_DENSITY */
#else
      rOP = 3.0
#endif /* RT_EXACT_EOS */
C
C  MH gas-phase reactions and correction for HI. Because the abundance
C  of H2+ and H- are always small, we ignore them as depositories of 
C  hydrogen atoms.
C
      pRateH2 = sFact1*sFact2*(pRate(iptabCiLW)+
     .     (pRate(iptabCi31)+rOP*pRate(iptabCi32))/(1+rOP))

#ifdef RT_TEST
C
C  RT-TEST +++++++++++++++++
C
      a(ieqXMH) = 0.0
      w(ieqXMH) = 0.0
C
C  RT-TEST -----------------
C
#else  /* RT_TEST */
#ifdef RT_HIGH_DENSITY
      cRate3B = den*XH1ll*(hRate(irateHDi30)*XH1ll +  !! 3-body reaction
     .     hRate(irateHDi31)*XMHll + hRate(irateHDi32)*XG1ll)
      cRate09 = cRate(irateMHi09)*XH1ll/
     .     (1+cRate(irateMHi09)*XH1ll*den/hRate(irateHDi09)) +
     .     cRate(irateMHi10)*XMHll/
     .     (1+cRate(irateMHi10)*XMHll*den/hRate(irateHDi10))
#else
      cRate3B = 0.0
      cRate09 = cRate(irateMHi09)*XH1ll + cRate(irateMHi10)*XMHll
#endif
C
C  Combining this all together
C
#ifdef RT_CHEMISTRY_MINIMAL_MODEL
C
C  Minimal mode:
C  production rate of HI = 2*destruction rate of H2;
C  production rate of H2 = 0.5*destruction rate of HI;
C
      aXMHl = pRateH2
      aXH1l = 2*cRate(irateMHiDU)*Xdust

      wXH1l = 2*aXMHl*XMHll

      wXMHl = 0.5*aXH1l*XH1ll

      a(ieqXH1) = a(ieqXH1) + uTime*aXH1l
      w(ieqXH1) = w(ieqXH1) + uTime*wXH1l

      a(ieqXMH) = uTime*aXMHl
      w(ieqXMH) = uTime*wXMHl

#else
C
C  Full (6- or 8-) species model
C
      aXH1l = cRate(irateMHi01)*Xfell + cRate(irateMHi02)*XHmll +
     .     cRate(irateMHi03)*XH2ll + cRate(irateMHi04)*XHpll +
     .     cRate(irateMHi26)*XG2ll + 2*cRate3B +
     .     2*cRate(irateMHiDU)*Xdust

      aXH2l = cRate(irateMHi03)*XH1ll + cRate(irateMHi07)*XMHll + 
     .     (cRate(irateMHi05)+cRate(irateMHi16))*XHmll +
     .     cRate(irateMHi27)*XG1ll

      aXG1l = cRate(irateMHi27)*XH2ll + cRate(irateMHi29)*XHmll

      aXG2l = (cRate(irateMHi24)+cRate(irateMHi25))*XMHll +
     .     cRate(irateMHi26)*XH1ll + cRate(irateMHi28)*XHmll

      aXMHl = cRate(irateMHi07)*XH2ll + cRate(irateMHi11)*XG1ll + 
     .     (cRate(irateMHi08)+cRate(irateMHi23))*Xfell + 
     .     cRate(irateMHi24)*XG2ll + cRate(irateMHi25)*XG2ll +
     .     cRate09 + pRate29 + pRateH2

      wXH1l = sFact1*(pRate(iptabCi27)*XHmll+pRate(iptabCi28)*XHpll) +
     .     2*pRateH2*XMHll + 2*cRate(irateMHi05)*XH2ll*XHMll +
     .     2*cRate(irateMHi06)*Xfell*XHpll + 
     .     cRate(irateMHi07)*XH2ll*XMHll +
     .     2*cRate(irateMHi08)*Xfell*XMHll + 2*cRate09*XMHll +
     .     2*cRate(irateMHi11)*XG1ll*XMHll +
     .     cRate(irateMHi14)*Xfell*XHmll +
     .     cRate(irateMHi15)*XH1ll*XHmll +
     .     cRate(irateMHi23)*Xfell*XMHll +
     .     cRate(irateMHi24)*XG2ll*XMHll +
     .     cRate(irateMHi27)*XG1ll*XH2ll +
     .     cRate(irateMHi28)*XG2ll*XHmll +
     .     cRate(irateMHi29)*XG1ll*XHmll

      wXH2l = cRate(irateMHi04)*XH1ll*XHpll + 
     .     cRate(irateMHi24)*XMHll*XG2ll +
     .     cRate(irateMHi26)*XH1ll*XG2ll

      wXG1l = aXG2l*XG2ll

      wXG2l = aXG1l*XG1ll

      wXMHl = XH1ll*(cRate(irateMHi02)*XHmll+cRate(irateMHi04)*XHpll+ 
     .     cRate3B+cRate(irateMHiDU)*Xdust)

      a(ieqXH1) = a(ieqXH1) + uTime*aXH1l
      w(ieqXH1) = w(ieqXH1) + uTime*wXH1l

      a(ieqXH2) = a(ieqXH2) + uTime*aXH2l
      w(ieqXH2) = w(ieqXH2) + uTime*wXH2l

      a(ieqXG1) = a(ieqXG1) + uTime*aXG1l
      w(ieqXG1) = w(ieqXG1) + uTime*wXG1l

      a(ieqXG2) = a(ieqXG2) + uTime*aXG2l
      w(ieqXG2) = w(ieqXG2) + uTime*wXG2l

      a(ieqXMH) = uTime*aXMHl
      w(ieqXMH) = uTime*wXMHl
#endif /* RT_CHEMISTRY_MINIMAL_MODEL */
#endif /* RT_TEST */
C
C  Heating due to gas reactions
C
      wEinll = wEinll + uTime*T1eV*0.4*XMHll*pRateH2
C
C  Cooling function - exchange reactions
C
      CF_MHl = XMHll*(cRate(irateMHc08)*Xfell+cRate(irateMHc09)*XH1ll+
     .     cRate(irateMHc10)*XMHll) + XHmll*(cRate(irateMHc14)*Xfell+
     .     cRate(irateMHc15)*XH1ll)
C
C  From Shapiro & Kang
C
#ifdef RT_HIGH_DENSITY
      CF_MHl = CF_MHl + 
     .     (1.602e-12/CFSTAR)*XH1ll*(3.53*XHmll*cRate(irateMHi02)+
     .     1.38*XHpll*cRate(irateMHi04))*den/(den+hRate(irateHDc02))
#endif
C
C  Cooling function - level excitations
C
      qo = cRate(irateMHcH1o)*XH1ll + cRate(irateMHcH2o)*XH2ll +
     .     cRate(irateMHcG1o)*XG1ll + cRate(irateMHcFEo)*Xfell +
     .     XMHll*(rOP*cRate(irateMHcMOo)+cRate(irateMHcMPo))/(1+rOP)

      qp = cRate(irateMHcH1p)*XH1ll + cRate(irateMHcH2p)*XH2ll +
     .     cRate(irateMHcG1p)*XG1ll + cRate(irateMHcFEp)*Xfell +
     .     XMHll*(rOP*cRate(irateMHcMOp)+cRate(irateMHcMPp))/(1+rOP)

      CF_MHl = CF_MHl + XMHll*(rOP*qo+qp)/(1+rOP)
#ifdef RT_HIGH_DENSITY
     .     /(1+den/hRate(irateHDcNC))
#endif
C
C ***********************************************
C
C  Heavy element (metals+dust)
C
C ***********************************************
C
C
C  1. High temperature part (SD93 fit)
C
      if(tem .gt. temZLT) then
         cfzHT = par(iparSDp1)*(cRate(irateZHcTA)+par(iparSDp2)*
     .        cRate(irateZHcTB))/(1.0+par(iparSDp3)*cRate(irateZHcTC))*
     .        Xfell*(XH2ll+XG2ll+XG3ll)
      else
         cfzHT = 0.0
      endif
C
C  2. Low temperature part (Preston, Dalgarno&McCray)
C
      cfzLT = par(iparSDp1)*XH1ll*(XH1ll*cRate(irateZLcTA)+
     .     Xfell*cRate(irateZLcTX))

      CF_HEl = cfzLT + cfzHT
C
C  3. dust
C
      CF_HEl = CF_HEl + cRate(irateMHcDU)*Xdust*Xfell

#ifdef RT_PAH_CR
C
C  4. PAHs (chi is in comoving units since de is!!!)
C
      chi = sFact1*1.05e10*pRate(iptabCiLW)
      psi = chi*sqrt(tem)/(max(de,1.0e-6*den)*d_pPAH)
      eps = 4.9e-2/(1+4.0e-3*psi**0.73) + 
     .     3.7e-2*(1.0e-4*tem)**0.7/(1+2.0e-4*psi)
      
      pRatePAHh = par(iparSDp1)*285*eps*exp(-d_sPAH*cdDUl)*
     .     pRate(iptabCiLW)
C
C  5. Cosmic rays heating and ionizations
C
      pRateCRi = 1.0e-17*par(iparDGas)*den
      pRateCRh = pRateCRi*232080  !! 232080 = 20eV in K

      a(ieqXH1) = a(ieqXH1) + uTime*pRateCRi
      w(ieqXH2) = w(ieqXH2) + uTime*XH1ll*pRateCRi

      wEinll = wEinll + uTime*((XH1ll+2*XMHll)*pRatePAHh+
     .     XH1ll*pRateCRh)
#endif
C
C  Final aseembly
C
      aEinll = aEinll + facCool*par(iparRhoB)*(CF_MHl+CF_HEl)

#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
      avg(iavgDust) = par(iparDGas)*(XH1ll+2*XMHll)
#endif

#else  /* RT_CHEMISTRY */

      a(ieqXMH) = 0.0
      w(ieqXMH) = 0.0

#if defined(RT_EXTERNAL_BACKGROUND) && (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_SELFCONSISTENT)
      avg(iavgDust) = 0.0
#endif

#endif /* RT_CHEMISTRY */
C
C  Turn cooling into heating below CMB (induced emission)
C
#if defined(RT_TEST) && (RT_TEST==1 || RT_TEST==11)

      a(ieqEin) = 0.0
      w(ieqEin) = 0.0

#else

      fac = 1.0/(tem+TemCMB)*(1.0-Tmin/tem)

      a(ieqEin) = (fac*aEinll*tem+aEinl2*par(iparDFac))/Einll
      w(ieqEin) = max(0.0,wEinll+wEinl2+fac*aEinll*TemCMB+par(iparDTad))

#endif

      return
      end
C
C
C
      function frtDustToGas(Zsol)
      include 'frt_base.inc'

#ifdef RT_FIXED_ISM
      frtDustToGas = D2Gmin
#else
      frtDustToGas = Zsol
#endif /* RT_FIXED_ISM */

      return
      end
C
C
C
      subroutine frtGetRecCaseH2(XH1ll,cdRecB,wa,wb)
      include 'frt_base.inc'
      include 'frt_tables.inc'

      taul = XH1ll*cdRecB*csH1AB
      wa = exp(-taul)
      wb = 1.0 - wa

      return
      end
C
C
C
      subroutine frtGetRecCaseG2(XG1ll,cdRecB,wa,wb)
      include 'frt_base.inc'
      include 'frt_tables.inc'

      taul = XG1ll*cdRecB*csG1AB
      wa = exp(-taul)
      wb = 1.0 - wa

      return
      end
C
C
C
      subroutine frtGetRecCaseG3(XG2ll,cdRecB,wa,wb)
      include 'frt_base.inc'
      include 'frt_tables.inc'

      taul = XG2ll*cdRecB*csG2AB
      wa = exp(-taul)
      wb = 1.0 - wa

      return
      end
C
C
C     
      subroutine frtGetPhotoRates(par,rf,itab,y,pRate)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_TRANSFER
#ifdef RT_UV
      include 'frt_transfer4.inc'
#else
      include 'frt_transfer3.inc'
#endif
#endif
#ifdef RT_FIXED_ISM
      include 'frt_spectrum.inc'
#endif /* RT_FIXED_ISM */
      dimension par(*), rf(*)
      dimension itab(*)
      dimension y(*), pRate(*)
#if defined(RT_TRANSFER) && defined(RT_TRANSFER_FLUX_CONSERVING)
      parameter (wFront = 0.5)
      parameter (TAUMIN = 1.0e-2)
      dimension rfFrnt(irfMax), rfBack(irfMax), drf(NSETS)
      dimension pRateBack(iptabDim)
#endif

#ifdef RT_1ZONE
C
C  Single zone plug-in
C
      call frtPhotoRates1Zone(par,rf,itab,y,pRate)
#else

#if defined(RT_TRANSFER) && defined(RT_TRANSFER_FLUX_CONSERVING)
C
C  Flux-conserving correction
C
      drf(1) = rf(irfH1wf) + y(ieqXH1)*rf(irfH1wv)
      drf(2) = rf(irfG1wf) + y(ieqXG1)*rf(irfG1wv)
      drf(3) = rf(irfG2wf) + y(ieqXG2)*rf(irfG2wv)

      if(max(drf(1),drf(2),drf(3)) .gt. TAUMIN) then

         rfFrnt(irfNorm) = rf(irfNorm)
         rfBack(irfNorm) = rf(irfNorm)
         do j=1,NSETS
            if(rf(irfOlf+j) .lt. 0.0) then
               rfFrnt(irfOlf+j) = -1.0
               rfBack(irfOlf+j) = -1.0
            else
               rfFrnt(irfOlf+j) = max(0.0,rf(irfOlf+j)+
     .              y(1+j)*rf(irfOlv+j)-wFront*drf(j))
               rfBack(irfOlf+j) = rfFrnt(irfOlf+j) + drf(j)
            endif
            if(rf(irfOgf+j) .lt. 0.0) then
               rfFrnt(irfOgf+j) = -1.0
               rfBack(irfOgf+j) = -1.0
            else
               rfFrnt(irfOgf+j) = max(0.0,rf(irfOgf+j)+
     .              y(1+j)*rf(irfOgv+j)-wFront*drf(j))
               rfBack(irfOgf+j) = rfFrnt(irfOgf+j) + drf(j)
            endif
         enddo

         call frtPhotoRatesTable(par,rfFrnt,itab,y,pRate,iptabMax)
         call frtPhotoRatesTable(par,rfBack,itab,y,pRateBack,iptabIon)
      
         if(drf(1) .gt. TAUMIN) then
            pRate(iptabPhH1) = max(0.0,pRate(iptabPhH1)-
     .           pRateBack(iptabPhH1))/drf(1)
            pRate(iptabPiH1) = max(0.0,pRate(iptabPiH1)-
     .           pRateBack(iptabPiH1))/drf(1)
         endif
         if(drf(2) .gt. TAUMIN) then
            pRate(iptabPhG1) = max(0.0,pRate(iptabPhG1)-
     .           pRateBack(iptabPhG1))/drf(2)
            pRate(iptabPiG1) = max(0.0,pRate(iptabPiG1)-
     .           pRateBack(iptabPiG1))/drf(2)
         endif
         if(drf(3) .gt. TAUMIN) then
            pRate(iptabPhG2) = max(0.0,pRate(iptabPhG2)-
     .           pRateBack(iptabPhG2))/drf(3)
            pRate(iptabPiG2) = max(0.0,pRate(iptabPiG2)-
     .           pRateBack(iptabPiG2))/drf(3)
         endif

      else
#endif /* RT_TRANSFER && RT_TRANSFER_FLUX_CONSERVING */
C
C  Generic rates
C
         call frtPhotoRatesTable(par,rf,itab,y,pRate,iptabMax)

#if defined(RT_TRANSFER) && defined(RT_TRANSFER_FLUX_CONSERVING)
      endif
#endif /* RT_TRANSFER && RT_TRANSFER_FLUX_CONSERVING */

#endif

#ifdef RT_FIXED_ISM
      umw = EPSG_UV/5.0e-6
      urf = pRate(iptabCiLW)*1.05e10
      if(urf .gt. 0) then
         q = umw/urf
         do j=1,iptabDim
            pRate(j) = pRate(j)*q
         enddo
      endif
#endif /* RT_FIXED_ISM */

      return
      end
C
C  -------------------------
C
C  Thermodynamic functions
C
C  -------------------------
C
#ifdef RT_EXACT_EOS
C
C  Exacts (non-monoatomic) EOS
C
#ifdef RT_HIGH_DENSITY
      function frtEin(tem,y,iLTE)
#else
      function frtEin(tem,y)
#endif
      include 'frt_base.inc'
      dimension y(*)

      Xtot = y(ieqXH1) + y(ieqXMH) + 2*y(ieqXH2) +
     .     y(ieqXG1) + 2*y(ieqXG2) + 3*y(ieqXG3)

#ifdef RT_HIGH_DENSITY
      frtEin = tem*(1.5*Xtot+y(ieqXMH)*frtH2_Ce(tem,iLTE))
#else
      frtEin = tem*(1.5*Xtot+y(ieqXMH)*frtH2_Ce(tem))
#endif

      return
      end
C
C
C
#ifdef RT_HIGH_DENSITY
      function frtTem(Ein,y,iLTE)
#else
      function frtTem(Ein,y)
#endif
      include 'frt_base.inc'
      dimension y(*)
C
C Use bi-section to insure convergence
C      
      Xtot = y(ieqXH1) + y(ieqXMH) + 2*y(ieqXH2) +
     .     y(ieqXG1) + 2*y(ieqXG2) + 3*y(ieqXG3)
      fH2 = y(ieqXMH)/Xtot

      if(Ein.lt.15.0 .or. fH2.lt.1.0e-3) then
         frtTem = Ein/(1.5*Xtot)
         return
      endif

      tem1 = Ein/(1.5*Xtot)    !!  atomic gas
      tem3 = tem1*(3.0/7.0)    !!  diatomic gas
C
C  This assumes that Ein is monotonic with T, i.e. 
C  f(tem1) >= Ein and f(tem3) <= Ein
C
      iter = 0

 10   continue

      tem2 = 0.5*(tem1+tem3)
C
#ifdef RT_HIGH_DENSITY
      Ein2 = Xtot*tem2*(1.5+fH2*frtH2_Ce(tem2,iLTE))
#else
      Ein2 = Xtot*tem2*(1.5+fH2*frtH2_Ce(tem2))
#endif

      if(Ein2 .ge. Ein) then
         tem1 = tem2
      endif
      
      if(Ein2 .le. Ein) then
         tem3 = tem2
      endif
      
      iter = iter + 1

      if(iter.lt.10 .and. 
     .     abs(tem1-tem2)/(1.0e-9+tem1+tem2).gt.1.0e-3) goto 10

      frtTem = 0.5*(tem1+tem3)

      return
      end
C
C
C
#ifdef RT_HIGH_DENSITY
      function frtGamma(Ein,y,iLTE)
#else
      function frtGamma(Ein,y)
#endif
      include 'frt_base.inc'
      dimension y(*)

      Xtot = y(ieqXH1) + y(ieqXMH) + 2*y(ieqXH2) +
     .     y(ieqXG1) + 2*y(ieqXG2) + 3*y(ieqXG3)
      fH2 = y(ieqXMH)/Xtot

      if(fH2 .lt. 1.0e-3) then
         cv = 1.5
      else
#ifdef RT_HIGH_DENSITY
         cv = 1.5 + fH2*frtH2_Cv(frtTem(Ein,y,iLTE),iLTE)
#else
         cv = 1.5 + fH2*frtH2_Cv(frtTem(Ein,y))
#endif
      endif

      frtGamma = (cv+1)/cv

      return
      end
C
C
C
#ifdef RT_HIGH_DENSITY
      function frtH2_Cv(tem,iLTE)
#else
      function frtH2_Cv(tem)
#endif

      if(tem .gt. 10.0) then
         x = tem/85.4
#ifdef RT_HIGH_DENSITY
         o2p = frtH2_op(tem,iLTE)
         if(o2p .lt. 0.0) then
#endif
C
C LTE
C
            crot = exp(-(0.7/x)**2)*(1+6.4/x**5)*(1.0-exp(-0.28*x**2))
#ifdef RT_HIGH_DENSITY
         else
C     
C     Fixed value
C
            qp = exp(-(2.1/x)**2)*(1+29/x**3.1)
            qo = exp(-(2.0/x)**4)
            crot = (qo*o2p+qp)/(o2p+1)
         endif
#endif
      else
         crot = 0
      endif

      if(tem .gt. 500.0) then
         x = 6100/tem
         cvib = x**2*exp(x)/(exp(x)-1)**2
      else
         cvib = 0
      endif

      frtH2_Cv = crot + cvib

      return
      end
C
C
C
#ifdef RT_HIGH_DENSITY
      function frtH2_Ce(tem,iLTE)
#else
      function frtH2_Ce(tem)
#endif

      if(tem .gt. 10.0) then
         x = tem/85.4
#ifdef RT_HIGH_DENSITY
         o2p = frtH2_op(tem,iLTE)
         if(o2p .lt. 0.0) then
#endif
C
C LTE
C
            crot = exp(-(1.2/x)**1.5)*(1+3.2/x**2.65)
#ifdef RT_HIGH_DENSITY
         else
C
C  Fixed value
C
            qp = exp(-(1.5/x)**2.5)
            qo = exp(-(3.5/x)**1.2) !!+ 2.0/x
            crot = (qo*o2p+qp)/(o2p+1)
         endif
#endif
      else
         crot = 0
      endif

      if(tem .gt. 500.0) then
         x = 6100/tem
         cvib = x/(exp(x)-1)
      else
         cvib = 0
      endif

      frtH2_Ce = crot + cvib

      return
      end
C
C
C
#ifdef RT_HIGH_DENSITY
      function frtH2_op(tem,iLTE)
#else
      function frtH2_op(tem)
#endif

#ifdef RT_HIGH_DENSITY
C
C  Is H2 formation on dust slower than ortho-para proton 
C  exchange reaction?
C
      if(iLTE .eq. 1) then
C
C  Assume LTE
C
#endif
         if(tem .gt. 10.0) then
            x = tem/85.4
            frtH2_op = 9*exp(-2/x)/(1+2*exp(-3/x-(1.75/x)**2.5))
         else
            frtH2_op = 0.0
         endif
#ifdef RT_HIGH_DENSITY
C
C  Assume a ratio of 3
C
      else
         frtH2_op = 3.0
      endif
#endif

      return
      end
#else  /* RT_EXACT_EOS */
C
C  EOS for monoatomic gas
C
      function frtEin(tem,y)
      include 'frt_base.inc'
      dimension y(*)
      
#ifdef RT_CHEMISTRY
      Xtot = y(ieqXH1) + y(ieqXMH) + 2*y(ieqXH2) +
     .     y(ieqXG1) + 2*y(ieqXG2) + 3*y(ieqXG3)
#else
      Xtot = y(ieqXH1) + 2*y(ieqXH2) +
     .     y(ieqXG1) + 2*y(ieqXG2) + 3*y(ieqXG3)
#endif

      frtEin = tem*1.5*Xtot

      return
      end
C
C
C
      function frtTem(Ein,y)
      include 'frt_base.inc'
      dimension y(*)

#if defined(RT_TEST) && (RT_TEST==1 || RT_TEST==11)

      Xtot = XH

#else

#ifdef RT_CHEMISTRY
      Xtot = y(ieqXH1) + y(ieqXMH) + 2*y(ieqXH2) +
     .     y(ieqXG1) + 2*y(ieqXG2) + 3*y(ieqXG3)
#else
      Xtot = y(ieqXH1) + 2*y(ieqXH2) +
     .     y(ieqXG1) + 2*y(ieqXG2) + 3*y(ieqXG3)
#endif

#endif

      frtTem = Ein/(1.5*Xtot)

      return
      end
C
C
C
      function frtGamma(Ein,y)
      include 'frt_base.inc'
      dimension y(*)

      frtGamma = 5.0/3.0

      return
      end
#endif /* RT_EXACT_EOS */
#ifdef RT_DEBUG
C
C  Debug helpers
C
      subroutine frtDebugPrint1(io,nStep,t,dt,tFin,y,ieqmax)
      include 'frt_base.inc'
      dimension y(*)

      write(io,9) nStep, t, dt, tFin, frtTem(y(ieqEin),y),
     .     y(ieqXH1)/XH, y(ieqXH2)/XH, y(ieqXMH)/(0.5*XH), 
     .     y(ieqmax), ieqmax

 9    format(i4,1p,8(1x,e11.3),1x,i1)

      return
      end
#endif

#endif
      subroutine frtDummyCooling(i)
      i = 0
      end
