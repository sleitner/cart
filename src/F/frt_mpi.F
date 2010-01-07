#include "../config.h"      
#ifdef RADIATIVE_TRANSFER
C
C  Sum/max local buffers and send the result to the global buffer
C
      subroutine frtGlobalAverage(iOp,nBuf,bufLoc,bufGlob)
      real*8 bufLoc(*), bufGlob(*)
#ifdef RT_PARALLEL_USE_MPI
      include 'mpif.h'
#endif /* RT_PARALLEL_USE_MPI */

#ifdef RT_PARALLEL_USE_MPI

      if(iOp .eq. 1) then
         call MPI_ALLREDUCE(bufLoc,bufGlob,nBuf,
     .        MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      else
         call MPI_ALLREDUCE(bufLoc,bufGlob,nBuf,
     .        MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      endif

      if(ierr .ne. MPI_SUCCESS) then
         write(0,*) 'RT_MPI: fatal error, ',
     .        'MPI communication error #',ierr
         call frtAbort
      endif

#else

      do j=1,nBuf
         bufGlob(j) = bufLoc(j)
      enddo

#endif /* RT_PARALLEL_USE_MPI */

      return
      end
C
C
C
      subroutine frtQueryId(id)
#ifdef RT_PARALLEL_USE_MPI
      include 'mpif.h'
#endif /* RT_PARALLEL_USE_MPI */

#ifdef RT_PARALLEL_USE_MPI

      call MPI_COMM_RANK(MPI_COMM_WORLD,id,ierr)
      if(ierr .ne. MPI_SUCCESS) then
         write(0,*) 'RT_MPI: fatal error, ',
     .        'MPI communication error #',ierr
         call frtAbort
      endif

#else

      id = 0

#endif /* RT_PARALLEL_USE_MPI */

      return
      end
C
C
C
      subroutine frtAbort()
#ifdef RT_PARALLEL_USE_MPI
      include 'mpif.h'
#endif /* RT_PARALLEL_USE_MPI */

C$OMP CRITICAL

      write(0,*) 'Abort.'
#ifdef RT_PARALLEL_USE_MPI
      call MPI_ABORT(MPI_COMM_WORLD,81,ierr)
#endif /* RT_PARALLEL_USE_MPI */
      stop

C$OMP END CRITICAL

      return
      end
#endif

      subroutine frtDummyMpi(i)
      i = 0
      end
