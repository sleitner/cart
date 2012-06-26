#include "frt_config.h"      
#ifdef RADIATIVE_TRANSFER
#ifdef RT_EXTERNAL_BACKGROUND
C
C  Add external background
C
      subroutine frtSetExternalBackground(aa,ang,ang2)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      dimension ang(NRAD)
#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)
      dimension ang2(NRLW)
#endif /*  defined(RT_CHEMISTRY) && defined(RT_LWBANDS) */

#if (RT_EXTERNAL_BACKGROUND == RT_BACKGROUND_SELFCONSISTENT)
#elif (RT_EXTERNAL_BACKGROUND == RT_BACKGROUND_HAARDT_MADAU)
      call frtBackgroundHM(aa,ang,ang2)
#else
#ifndef RT_TEST
#error "Invalid value for RT_EXTERNAL_BACKGROUND"
#endif
#endif

      return
      end
#if (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_HAARDT_MADAU)
C
C
C
      subroutine frtBackgroundHM(aa,ang,ang2)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      dimension ang(NRAD)
#if defined(RT_CHEMISTRY) && defined(RT_LWBANDS)
      dimension ang2(NRLW)
#endif /*  defined(RT_CHEMISTRY) && defined(RT_LWBANDS) */
      include 'frt_files.inc'
C
      parameter(NHMZ = 50)
      parameter(NHMW = 432)
      parameter(NHMP = 2)
      common /RT_BKG_HM_0/ ip0
      common /RT_BKG_HM_1/ z(NHMZ), d(NRAD,NHMZ)
      dimension q(NHMP), w(NHMW), xiw(NHMW), dw(NHMW)
      dimension iw1(NRAD), iw2(NRAD)
C
      logical inited
      save inited
      data inited/.false./, ip0/2/

      if(.not. inited) then
C
C  Load the data on the first call
C
         inited = .true.

         open(unit=IOUNIT, file=fileHM, status='old')
         read(IOUNIT,*) np1, nw1, nz1
         if(np1.ne.NHMP .or. nw1.ne.NHMW .or. nz1.ne.NHMZ) then
            write(0,*) 'RT::BackgroundHM: fatal error, ',
     .           'corrupted bkg_HM.dat file:'
            write(0,*) '> NP= in file: ', np1, ' in code: ', NHMP
            write(0,*) '> NW= in file: ', nw1, ' in code: ', NHMW
            write(0,*) '> NZ= in file: ', nz1, ' in code: ', NHMZ
            close(IOUNIT)
            call frtAbort
         endif

         do k=1,NHMZ
            do j=1,NHMW

               read(IOUNIT,*) z1, w1, (q(i),i=1,NHMP)

               if(j .eq. 1) then
                  z(k) = z1
                  uJnu = 1.5808e-17*(1.0+z1)**3
               else
                  if(abs((1.0+z1)/(1.0+z(k))-1.0) .gt. 1.0e-4) then
                     write(0,*) 'RT::BackgroundHM: fatal error, ',
     .                    'corrupted bkg_HM.dat file:'
                     write(0,*) '> Wrong z: ', k, j, z1, z(k)
                     close(IOUNIT)
                     call frtAbort
                  endif
               endif

               if(k .eq. 1) then
                  w(NHMW+1-j) = w1
                  xiw(NHMW+1-j) = log(912.0/w1)
               else
                  if(abs(w1/w(NHMW+1-j)-1.0) .gt. 1.0e-4) then
                     write(0,*) 'RT::BackgroundHM: fatal error, ',
     .                    'corrupted bkg_HM.dat file:'
                     write(0,*) '> Wrong w: ', k, j, w1, w(NHMW+1-j)
                     close(IOUNIT)
                     call frtAbort
                  endif
               endif

               dw(NHMW+1-j) = q(ip0)/uJnu

            enddo
C
C  Interpolate onto FRT mesh
C
            if(k .eq. 1) then
               ju = 1
               do lr=1,NRAD
                  do while(ju.lt.NHMW .and. xiw(ju).lt.xig(lr))
                     ju = ju + 1
                  enddo
                  if(ju .eq. 1) then
                     iw1(lr) = 1
                     iw2(lr) = 2
                  else
                     iw1(lr) = ju - 1
                     iw2(lr) = ju
                  endif
               enddo
            endif

            do lr=1,NRAD
               d(lr,k) = dw(iw1(lr)) + (dw(iw2(lr))-dw(iw1(lr)))/
     .              (xiw(iw2(lr))-xiw(iw1(lr)))*(xig(lr)-xiw(iw1(lr)))
            enddo

         enddo

         close(IOUNIT)

      endif
C
C  Find a proper z
C
      z1 = 1.0/aa - 1.0
      if(z1 .lt. z(1)) then
         z1 = z(1)
         kl = 1
         ku = 2
      else if(z1 .gt. z(NHMZ)) then
         z1 = z(NHMZ)
         kl = NHMZ - 1
         ku = NHMZ
      else
         ku = 0
         do k=1,NHMZ
            if(ku.eq.0 .and. z1.lt.z(k)) ku = k
         enddo
         if(ku .eq. 1) then
            kl = 1
            ku = 2
         else
            kl = ku - 1
         endif
      endif

      qu = (z1-z(kl))/(z(ku)-z(kl))
      ql = (z(ku)-z1)/(z(ku)-z(kl))
      do lr=1,NRAD
         ang(lr) = max(0.0,ql*d(lr,kl)+qu*d(lr,ku))
      enddo

#ifdef RT_LWBANDS
      if(lrLWtr .ge. lrH1tr) then
         do lr=1,NRLW
            angLW(lr) = ang(lrLWtr)
         enddo
      else
         dxi = xig(lrH1tr-1) - xig(lrLWtr-1)
         r = (ang(lrH1tr-1)-ang(lrLWtr-1))/dxi
         do lr=1,NRLW
            ang2(lr) = ang(lrLWtr-1) + r*
     .           min(dxi,xigLW(lr)-xig(lrLWtr-1))
         enddo
      endif
#endif /* RT_LWBANDS */

      return
      end
#endif /* (RT_EXTERNAL_BACKGROUND==RT_BACKGROUND_HAARDT_MADAU) */

#endif /* RT_EXTERNAL_BACKGROUND */
#endif /* RADIATIVE_TRANSFER */
      subroutine frtDummyBackground(i)
      i = 0
      end
