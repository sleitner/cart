#include "defs.h"      
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
#ifdef RT_TEST
C
C
C
      subroutine frtTestInitSourceSpectrum
      include 'frt_base.inc'
      include 'frt_tables.inc'
      include 'frt_spectrum.inc'
C
C  Check that Yp is set to 0
C
      Yp = 4*XG
      if(Yp .gt. 1.0e-10) then
         write(0,*) 'RT_Testing: fatal error, ',
     .        'supported tests only work with Yp<1e-10.'
         call frtAbort
      endif

#if (RT_TEST==1 || RT_TEST==11)
C
C  Monochromatic
C
      do lr = 1,NRAD
         guv_Q(lr) = 0.0
         guv_G(lr) = 0.0
      enddo
      guv_G(lrH1tr) = 1.0/XISTEP
#else
C
C  Black body
C
      do lr = 1,NRAD
         guv_Q(lr) = 0.0
         guv_G(lr) = Txi(lr)**3/(exp(Txi(lr)/1.0e5)-1.0)
      enddo
#endif

      sum_G = 0.0
      do lr = lrH1tr,lrmax
         sum_G = sum_G + guv_G(lr)
      enddo
      sum_G = sum_G*XISTEP
C
C  Now we have integrals \int g_\xi d\xi. Normalize them
C
      fac_G = 1.0/sum_G
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(lr),
C$OMP+SHARED(guv_G,guv_Q,fac_G,fac_Q)
      do lr = 1,NRAD
         guv_G(lr) = fac_G*guv_G(lr)
      enddo

#if (RT_TEST == -1)
C
C  Test #0
C
      do lr = 1,NRAD
         angUni(lr) = 1.0e12/3e10*guv_G(lr)
      enddo

      EPSG_UV = 0.0
      EPSQ_UV = 0.0
#else
C
C  Other tests
C
      EPSG_UV = 1.0
      EPSQ_UV = 0.0
#endif

      return
      end
C
C
C
      subroutine frtTestStepBegin(Om0,h100,r0,t0)
      include 'frt_base.inc'
      include 'frt_tables.inc'
C
C  Test #0
C
#if (RT_TEST == -1)
      if(t0*tStep .gt. 0.7e5*(10**0.1-1.0)) then
         do lr = 1,NRAD
            angUni(lr) = 0.0
         enddo
      endif

      write(6,*) angUni(lrH1tr)

      if(t0*tStep .gt. 4.0e7) then
         call frtAbort
      endif
#endif

      return
      end
C
C
C
      subroutine frtTestModifyTimeStep(dt,ns)

#if (RT_TEST == -1)
C
C  Test #0
C
      if(ns .eq. 0) then
         dt = dt
      else if(ns .eq. 1) then
         dt = dt*(10.0**0.1-1)
      else
         dt = dt*10.0**0.1
      endif
#else
C
C  Other tests
C
      save nsRamp
      save dtMarch
      data nsRamp/10/, dtMarch/-1.0/

      if(ns .lt. nsRamp) then
         if(ns .eq. 0) then
            r = 10.0**(0.2-0.2*nsRamp)
            dtMarch = dt
         else
            r = 10.0**(0.2*(ns-nsRamp))*(10.0**0.2-1.0)
            if(dtMarch .lt. 0.0) then !!  check for a restart
               dtMarch = dt/r*10**0.2
            endif
         endif
         dt = dtMarch*r
      else
         if(dtMarch .lt. 0.0) then !!  check for a restart
            dtMarch = dt
         endif
         dt = dtMarch
      endif
#endif
  
      return
      end
#ifdef RT_EXTERNAL_BACKGROUND
C
C
C
      subroutine frtTestOutputBackground(z)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)
      include 'frt_tablesLW.inc'
#endif /*  defined(RT_CHEMISTRY) && defined(RT_LWBANDS) */
      character*99 str

      aa = 1.0/(1.0+z)
      call frtSetExternalBackground(aa,angAvg,angLW)

      uJnu = 1.5808e-17/aa**3
      
      write(str,'(f4.2)') z
      is = index(str,' ') - 1
      open(unit=81, file='test-bkg-'//str(1:is)//'.res')
      do lr=1,NRAD
         write(81,*) 911.75*exp(-xig(lr)), angAvg(lr)*uJnu
      enddo
      close(81)

#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)
      open(unit=81, file='test-bkgLW-'//str(1:is)//'.res')
      do lr=1,NRLW
         write(81,*) 911.75*exp(-xigLW(lr)), angLW(lr)*uJnu
      enddo
      close(81)
#endif /*  defined(RT_CHEMISTRY) && defined(RT_LWBANDS) */

      return
      end
#endif  /* RT_EXTERNAL_BACKGROUND */
#endif  /* RT_TEST */
#endif

      subroutine frtDummyTesting(i)
      i = 0
      end
