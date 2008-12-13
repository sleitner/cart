#include "defs.h"      
#ifdef RADIATIVE_TRANSFER
#include "../rt_config.h"
C
C  Sum/max local buffers and send the result to the global buffer
C
      subroutine frtGlobalAverage(mode,nBuf,bufLoc,bufGlob)
      real*8 bufLoc(*), bufGlob(*)
#ifdef RT_PARALLEL_USE_MPI
      include 'mpif.h'
#endif /* RT_PARALLEL_USE_MPI */

#ifdef RT_PARALLEL_USE_MPI

      if(mode .eq. 1) then
         call MPI_ALLREDUCE(bufLoc,bufGlob,nBuf,
     .        MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      else
         call MPI_ALLREDUCE(bufLoc,bufGlob,nBuf,
     .        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif

      if(ierr .ne. MPI_SUCCESS) then
         write(0,*) 'RT_Mpi: fatal error, ',
     .        'MPI communication error #',ierr
         stop
      endif

#else

      do j=1,nBuf
         bufGlob(j) = bufLoc(j)
      enddo

#endif /* RT_PARALLEL_USE_MPI */

      return
      end

#endif

      subroutine frtDummyMpi(i)
      i = 0
      end
