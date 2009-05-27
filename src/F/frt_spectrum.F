#include "defs.h"      
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
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
C  IFPOP = 3 Jason pop3 IFPOP = 2 nick pop2
C
      parameter(G_Abrk = 5.5)
      parameter(G_Bbrk = 2.5)
      parameter(G_Xray = 2.0e-3)
      parameter(G_pAbrk = 2.5)
      parameter(G_pBbrk = -4.0)
      parameter(G_pXray = 2.0e1)
C
      dimension guv_G_esc(NRAD), guv_Q_esc(NRAD)
C
#ifndef RT_TEST
C
C  Set switches - must be specified somewhere else???
C  
C  The default value: 5.0e-6 is obtained as follows:
C    Draine field (1.0e6 phot/cm^2/s/ster/eV) at 12.0 eV is
C      n_\xi (=\nu n_\nu) = 0.005/cm^3.
C    In the code, for LW band (in the main rad. grid),
C      n_\xi = epsUV*SFR3/r1^2/E30*1.0e3 /cm^3,
C    where
C      SFR3 = SFR/(3 Msun/yr)
C      r1   = <r>/(1 kpc) (averaged as 1/<r>^2 = \int dV 1/r^2)
C      E30  = <E_\nu>/(30 eV)
C    Assuming SFR3 = r1 = E30 = 1, we get epsUV = 5.0e-6
C
      EPSG_UV = 5.0e-6
      EPSQ_UV = 0.0
C
C  Spectral shape parameters
C
      IPOP = 2

      FESCSET = 1.0
      XG1HESC = 0.1
      XG2HESC = 0.01
C
C  Start initialization
C
      if(IPOP.ne.2 .and. IPOP.ne.3) then
         write(0,*) 'RT_Spectrum: ',
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

         if(fnu .lt. 1.0e+3) then
            guvl_Q = 0.0123*fnu**3/(exp(min(50.0,real(fnu)/2.59))-1.0)+
     .           fnu**aeQSO
         else
            if(fnu .lt. 1.0e+5) then
               guvl_Q = 0.5*fnu**aeQSO*(1.0+(fnu/1.0e+3)**0.7)
            else
               guvl_Q = 0.5*fnu**aeQSO*(1.0+(fnu/1.0e+3)**0.7)/
     .              (1.0+(fnu/1.0e+5))*2.0
            endif
         endif
         guvl_Q = 0.2675*guvl_Q

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
C
C  **MASSI**
C  f_esc and epsilon_uv
C
      IF(FESCSET.gt.-1.0 .and. FESCSET.lt.1.0) then

         if(abs(FESCSET) .lt. 1.0e-5) then
            write(0,*) 'RT_Spectrum: fatal error, FESCSET = ', FESCSET,
     .           '; such a small an escape fraction makes no sense.'
            call frtAbort
         endif

         iter = 0

         sum_Q = 0.0
         sum_G = 0.0
COMPILER_DIRECTIVE_SERIAL_LOOP
         do lr = lrH1tr,lrmax
            sum_Q = sum_Q + guv_Q(lr)
            sum_G = sum_G + guv_G(lr)
         enddo
         sum0_Q = sum_Q*XISTEP
         sum0_G = sum_G*XISTEP

         xg0 = XG1HESC
         xg1 = XG2HESC
         fesc = abs(FESCSET)

         aN_a = 0.0
         aN_b = 1.0e5/csH1tr
         
 7777    iter = iter + 1

         aN_c = 0.5*(aN_a+aN_b)
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr,tau),
C$OMP+SHARED(aN_c,csH1,xg0,csG1,xg1,csG2,guv_Q_esc,guv_G_esc,
C$OMP+ guv_G,guv_Q,lrH1tr,lrmax)
         do lr = lrH1tr,lrmax
            tau = aN_c*(csH1(lr)+xg0*csG1(lr)+xg1*csG2(lr))
            guv_Q_esc(lr) = guv_Q(lr)*exp(-tau)  
            guv_G_esc(lr) = guv_G(lr)*exp(-tau)  
         enddo

         sum_Q = 0.0
         sum_G = 0.0
         do lr = lrH1tr,lrmax
            sum_Q = sum_Q + guv_Q_esc(lr)
            sum_G = sum_G + guv_G_esc(lr)
         enddo
         sum_Q = sum_Q*XISTEP
         sum_G = sum_G*XISTEP

         fesc_Q = sum_Q/sum0_Q
         fesc_G = sum_G/sum0_G

         fe_c = fesc_G*1.0 + fesc_Q*0.0  

         if(fesc .le. fe_c) aN_a = aN_c
         if(fesc .ge. fe_c) aN_b = aN_c

         if(abs(fe_c/fesc-1.0) .gt. 1.0e-3) then

            if(iter .gt. 99) then
               write(0,*) 'RT_Spectrum: fatal error, FESC does not ',
     .              'converge ', fesc, iter
               call frtAbort
            endif
            goto 7777

         else

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(guv_Q_esc,guv_G_esc,
C$OMP+ guv_G,guv_Q,lrH1tr)
            do lr = lrH1tr,NRAD
               guv_Q(lr) = guv_Q_esc(lr)
               guv_G(lr) = guv_G_esc(lr)
            enddo

         endif

      ENDIF

#else

      call frtTestInitSourceSpectrum

#endif  /* RT_TEST */

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
