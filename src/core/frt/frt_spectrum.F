#include "frt_config.h"      
#ifdef RADIATIVE_TRANSFER
C
C
C
      subroutine frtInitSourceSpectrum
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_spectrum.inc'
C
C  Parameters of QSO spectral shape
C
      parameter(aeQSO = -1.3)
C
C  Parameters of stellar spectral shape
C  **MASSI**
C  IPOP = 3 Jason pop3 IPOP = 2 nick pop2
C
      parameter(G_Abrk = 5.5)
      parameter(G_Bbrk = 2.5)
      parameter(G_Xray = 2.0e-3)
      parameter(G_pAbrk = 2.5)
      parameter(G_pBbrk = -4.0)
      parameter(G_pXray = 2.0e1)
C
#ifndef RT_TEST
      EPSG_UV = fGal
      EPSQ_UV = fQSO
C
C  Start initialization
C
      if(IPOP.ne.2 .and. IPOP.ne.3) then
         write(0,*) 'RT::InitSourceSpectrum: ',
     .        'fatal error, invalid IPOP value: ', IPOP
         call frtAbort
      endif
C
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,x,fnu,guvl_G,guvl_Q)
C$OMP+SHARED(Txi,IPOP,guv_G,guv_Q)
      do lr = 1,NRAD

         x = Txi(lr)/TH1
         fnu = Txi(lr)/T1eV

	 IF(IPOP .eq. 2) then 
            if(x .lt. 1.0) then
               guvl_G = G_Abrk
            else if(x .lt. 2.5) then
               guvl_G = 1.0/x**1.8
            else if(x .lt. 4.0) then
               guvl_G = 1.0/(G_Bbrk*x**1.8)
            else if(x .lt. 100.0) then
               guvl_G = G_xray*x**3/(exp(x/1.4)-1.0)
            else
               guvl_G = 0.0
            endif
            guvl_G = 0.09517*guvl_G
	 ELSE IF(IPOP .eq. 3) then 
            if(x .lt. 1.0) then
               guvl_G = x**1.4*G_pAbrk
            else if(x .lt. 4.) then
               guvl_G = exp(G_pBbrk/x)/exp(G_pBbrk)/x**2.
            else
               guvl_G = G_pXray*exp(-x/4.)/x**2.
            endif
            guvl_G = 0.013283*guvl_G
         ELSE
            guvl_G = 0.0
         ENDIF
C
C  Average QSO spectral shape used by Hopkins et al 2007, ApJ, 654, 731
C  (i.e. my fit to it).
C
         guvl_Q = 9.5e-5/(1+(fnu/300))**0.8*exp(-fnu/5e5) +
     .        0.1*(10/fnu)**2.2/(1+(9/fnu)**5)**0.4 +
     .        5.8*exp(-fnu/0.2)

         guv_G(lr) = guvl_G/T1eV
         guv_Q(lr) = guvl_Q/T1eV

      enddo
C
C  Normalize spectral shapes:
C  int_1^\infty E_\nu g_\xi d\xi = m_p c^2 = 938.3MeV
C
      sum_Q = 0.0
      sum_G = 0.0
      do lr = lrH1tr,lrmax
         fnu = Txi(lr)/T1eV
         sum_Q = sum_Q + guv_Q(lr)*fnu
         sum_G = sum_G + guv_G(lr)*fnu
      enddo
      sum_Q = sum_Q*XISTEP
      sum_G = sum_G*XISTEP
C
C  Now we have integrals \int E_\nu g_\xi d\xi in eV. Normalize them
C
      fac_Q = 9.383e8/sum_Q
      fac_G = 9.383e8/sum_G
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(guv_G,guv_Q,fac_G,fac_Q)
      do lr = 1,NRAD
         guv_Q(lr) = fac_Q*guv_Q(lr)
         guv_G(lr) = fac_G*guv_G(lr)
      enddo

#else  /* RT_TEST */

      call frtTestInitSourceSpectrum

#endif /* RT_TEST */

#ifdef RT_OUTPUT
C
C  Compute the average photon energies for diagnostics
C
      sum0_Q = 0.0
      sum0_G = 0.0
      sum1_Q = 0.0
      sum1_G = 0.0
      do lr = lrH1tr,lrmax
         fnu = Txi(lr)/T1eV
         sum0_Q = sum0_Q + guv_Q(lr)
         sum0_G = sum0_G + guv_G(lr)
         sum1_Q = sum1_Q + guv_Q(lr)*fnu
         sum1_G = sum1_G + guv_G(lr)*fnu
      enddo
      EAvg_Q = sum1_Q/(VSMALL+sum0_Q)
      EAvg_G = sum1_G/(VSMALL+sum0_G)

      if(ioNode) then
         open(unit=IOUNIT, file='rt-run.log', access='append')
         write(IOUNIT,*) 'Spectrum: <E>_G, eps_G = ', EAvg_G, EPSG_UV
         write(IOUNIT,*) 'Spectrum: <E>_Q, eps_Q = ', EAvg_Q, EPSQ_UV
         close(IOUNIT)
      endif
#endif /* RT_OUTPUT */

      return
      end
C
C
C
      subroutine frtSetSourceSpectrum(a,guv)
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_spectrum.inc'
      dimension guv(*)
C
C  Units of the source function:
C     dn_\xi dt = guveff drho_*/dt
C     guveff = eps_uv * n_b * guv, where guv is normalized:
C     int_o^\infty guv_\xi E dxi = m_p c^2
C     Then eps_uv is the fraction of rest mass emitted above 1 Ryd
C
      do lr=1,NRAD
         guv(lr) = EPSG_UV*guv_G(lr) + EPSQ_UV*guv_Q(lr)
      enddo

      return
      end
#endif
      subroutine frtDummySpectrum(i)
      i = 0
      end
