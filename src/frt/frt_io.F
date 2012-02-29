#include "../config.h"      
#ifdef RADIATIVE_TRANSFER
C
C  IO routines
C
      subroutine frtWriteRadiationFields(iunit)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
C
C  Write into existing Fortran binary file
C
      write(iunit) (angUni(j),j=1,NRAD)
      write(iunit) (angLoc(j),j=1,NRAD)
#ifdef RT_LWBANDS
      write(iunit) (angLW(j),j=1,NRLW)
#endif

      return
      end
C
C
C
      subroutine frtReadRadiationFields(iunit)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
C
C  Read from existing Fortran binary file
C
      read(iunit,err=10) (angUni(j),j=1,NRAD)
      read(iunit,err=10) (angLoc(j),j=1,NRAD)
#ifdef RT_LWBANDS
      read(iunit,err=20) (angLW(j),j=1,NRLW)
#endif

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(j),
C$OMP+SHARED(angUni,angLoc,angAvg)
      do j=1,NRAD
         angAvg(j) = angUni(j) + angLoc(j)
      enddo

      return

 10   write(6,*) 'Error reading radiation field from a file.'
      return

#ifdef RT_LWBANDS
 20   write(6,*) 'Error reading LW bands from a file.'
      return
#endif

      end
C
C
C
      subroutine frtPackRadiationFields(n,data)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      dimension data(*)

      n0 = 2*NRAD
#ifdef RT_LWBANDS
      n0 = n0 + NRLW
#endif
      if(n .ne. n0) then
         n = n0
      else
         do j=1,NRAD
            data(j) = angUni(j)
         enddo
         do j=1,NRAD
            data(NRAD+j) = angLoc(j)
         enddo
#ifdef RT_LWBANDS
         do j=1,NRLW
            data(2*NRAD+j) = angLW(j)
         enddo
#endif
      endif

      return
      end
C
C
C
      subroutine frtUnpackRadiationFields(n,data)
      include 'frt_base.inc'
      include 'frt_tables.inc'
#ifdef RT_LWBANDS
      include 'frt_tablesLW.inc'
#endif
      dimension data(*)

      n0 = 2*NRAD
#ifdef RT_LWBANDS
      n0 = n0 + NRLW
#endif
      if(n .ne. n0) then
         n = n0
      else
         do j=1,NRAD
            angUni(j) = data(j)
         enddo
         do j=1,NRAD
            angLoc(j) = data(NRAD+j)
         enddo
#ifdef RT_LWBANDS
         do j=1,NRLW
            angLW(j) = data(2*NRAD+j)
         enddo
#endif

C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(j),
C$OMP+SHARED(angUni,angLoc,angAvg)
         do j=1,NRAD
            angAvg(j) = angUni(j) + angLoc(j)
         enddo

      endif

      return
      end

#endif
      subroutine frtDummyIO(i)
      i = 0
      end
