#include "../config.h"
#ifdef RADIATIVE_TRANSFER

#undef WA2_METHOD2

C
C  ODE solvers for derCoolOff routine
C
      subroutine frtCoolOffStepRK2(nEqs,dt,var,rf,itab,
     .     w0,a0,rs,d1,out1,pRate,chf)
      dimension var(*), rf(*), itab(*)
      dimension w0(*), a0(*), rs(*)
      dimension out1(*), pRate(*), chf(2,*)
C
      dimension y0(9), y1(9), w1(9), a1(9)
C
      d1 = 0.5
      
      do j=1,nEqs
         y0(j) = var(j)
         y1(j) = y0(j) + dt*(w0(j)-a0(j)*y0(j))
      enddo

      do j=1,nEqs
         var(j) = y1(j)
      enddo
      call frtRates(var,rf,itab,w1,a1,out1,pRate,chf)

      do j=1,nEqs
         var(j) = y0(j) + 0.5*dt*(w0(j)-a0(j)*y0(j)+w1(j)-a1(j)*y1(j))
      enddo

      do j=1,nEqs
         rs(j) = 0.5*dt*(w0(j)-a0(j)*y0(j)-w1(j)+a1(j)*y1(j))
      enddo
 
      return
      end
C
C
C
      subroutine frtCoolOffStepWA1(nEqs,dt,var,rf,itab,
     .     w0,a0,rs,d1,out1,pRate,chf)
      dimension var(*), rf(*), itab(*)
      dimension w0(*), a0(*), rs(*)
      dimension out1(*), pRate(*), chf(2,*)
C
      dimension y0(9), y1(9), y2(9), w1(9), a1(9)
C
      parameter(alp = 1.0/3.0)
      parameter(bet = 2.0/3.0)
C
      d1 = 0.25

      do j=1,nEqs
         y0(j) = var(j)
         tau = bet*dt*a0(j)
         if(tau .lt. 0.01) then
            fun2 = 1 - tau/2 + tau**2/6
         else
            fun2 = (1-exp(-tau))/tau
         endif
         y1(j) = y0(j)*exp(-tau) + bet*dt*w0(j)*fun2
         tau = dt*a0(j)
         if(tau .lt. 0.01) then
            fun2 = 1 - tau/2 + tau**2/6
         else
            fun2 = (1-exp(-tau))/tau
         endif
         y2(j) = y0(j)*exp(-tau) + dt*w0(j)*fun2
      enddo

      do j=1,nEqs
         var(j) = y1(j)
      enddo
      call frtRates(var,rf,itab,w1,a1,out1,pRate,chf)

      do j=1,nEqs
         tau = dt*(a1(j)*y1(j)+alp*a0(j)*y0(j))/(y1(j)+alp*y0(j))
         if(tau .lt. 0.01) then
            fun2 = 1 - tau/2 + tau**2/6
         else
            fun2 = (1-exp(-tau))/tau
         endif
         var(j) = y0(j)*exp(-tau) + dt*(0.25*w0(j)+0.75*w1(j))*fun2
      enddo

      do j=1,nEqs
         dy = y1(j) - y0(j)
         if(dy .ne. 0.0) then
            rs(j) = 4*(var(j)-y2(j))**2/(9*dy) + 
     .           dt*a0(j)*(var(j)-y2(j))/6
         else
            rs(j) = 0.0
         endif
      enddo

      return
      end
C
C
C
      subroutine frtCoolOffStepWA2(nEqs,dt,var,rf,itab,
     .     w0,a0,rs,d1,out1,pRate,chf)
      dimension var(*), rf(*), itab(*)
      dimension w0(*), a0(*), rs(*)
      dimension out1(*), pRate(*), chf(2,*)
C
      dimension y0(9), y1(9), y2(9), w1(9), a1(9)
C
      parameter(alp = 0.25)
      parameter(bet = 2.0/3.0)
C
      d1 = 0.25

      do j=1,nEqs
         y0(j) = var(j)
         tau = bet*dt*a0(j)
         if(tau .lt. 0.01) then
            fun2 = 1 - tau/2 + tau**2/6
         else
            fun2 = (1-exp(-tau))/tau
         endif
         y1(j) = y0(j)*exp(-tau) + bet*dt*w0(j)*fun2
         tau = dt*a0(j)
         if(tau .lt. 0.01) then
            fun2 = 1 - tau/2 + tau**2/6
         else
            fun2 = (1-exp(-tau))/tau
         endif
         y2(j) = y0(j)*exp(-tau) + dt*w0(j)*fun2
      enddo

      do j=1,nEqs
         var(j) = y1(j)
      enddo
      call frtRates(var,rf,itab,w1,a1,out1,pRate,chf)

      do j=1,nEqs
         tau = dt*(alp*a0(j)+(1-alp)*a1(j))
         if(tau .lt. 0.01) then
            fun2 = 1 - tau/2 + tau**2/6
         else
            fun2 = (1-exp(-tau))/tau
         endif
         var(j) = y0(j)*exp(-tau) + dt*(0.25*w0(j)+0.75*w1(j))*fun2
      enddo

      do j=1,nEqs
         dy = y1(j) - y0(j)
         if(dy .ne. 0.0) then
            rs(j) = 4*(var(j)-y2(j))**2/(9*dy) + 
     .           0.125*dt**2*(a0(j)*w1(j)-a1(j)*w0(j))
         else
            rs(j) = 0.0
         endif
      enddo

      return
      end
C
C
C
      subroutine frtCoolOffStepQS1(nEqs,dt,var,rf,itab,
     .     w0,a0,rs,d1,out1,pRate,chf)
      dimension var(*), rf(*), itab(*)
      dimension w0(*), a0(*), rs(*)
      dimension out1(*), pRate(*), chf(2,*)
C
      dimension y0(9), y1(9), w1(9), a1(9), ys(9)
C
      d1 = 0.5
C
      do j=1,nEqs
         y0(j) = var(j)
         tau = dt*a0(j)
         if(tau .lt. 0.01) then
            alp = 0.5 + tau/12 - tau**3/720
         else
            alp = (1-(1-exp(-tau))/tau)/(1-exp(-tau))
         endif
         y1(j) = y0(j) + dt*(w0(j)-a0(j)*y0(j))/(1+alp*tau)
         ys(j) = y1(j)
      enddo

      do j=1,nEqs
         var(j) = y1(j)
      enddo
      call frtRates(var,rf,itab,w1,a1,out1,pRate,chf)

      do j=1,nEqs
         tau = 0.5*dt*(a0(j)+a1(j))
         if(tau .lt. 0.01) then
            alp = 0.5 + tau/12 - tau**3/720
         else
            alp = (1-(1-exp(-tau))/tau)/(1-exp(-tau))
         endif
         var(j) = y0(j) + dt*(alp*w1(j)+(1-alp)*w0(j)-
     .        0.5*(a0(j)+a1(j))*y0(j))/(1+alp*tau)
      enddo

      do j=1,nEqs
         rs(j) = var(j) - ys(j)
      enddo

      return
      end
C
C
C
      subroutine frtCoolOffStepQSN(nEqs,dt,var,rf,itab,
     .     w0,a0,rs,d1,out1,pRate,chf)
      dimension var(*), rf(*), itab(*)
      dimension w0(*), a0(*), rs(*)
      dimension out1(*), pRate(*), chf(2,*)
C
      dimension y0(9), y1(9), w1(9), a1(9), ys(9)
C
      d1 = 0.5
C
      nCor = 2
      rCor = 1.0

      do j=1,nEqs
         y0(j) = var(j)
         tau = dt*a0(j)
         if(tau .lt. 0.01) then
            alp = 0.5 + tau/12 - tau**3/720
         else
            alp = (1-(1-exp(-tau))/tau)/(1-exp(-tau))
         endif
         y1(j) = y0(j) + dt*(w0(j)-a0(j)*y0(j))/(1+alp*tau)
         ys(j) = y1(j)
      enddo

      iCor = 0
 20   iCor = iCor + 1

      do j=1,nEqs
         var(j) = y1(j)
      enddo
      call frtRates(var,rf,itab,w1,a1,out1,pRate,chf)

      do j=1,nEqs
         tau = 0.5*dt*(a0(j)+a1(j))
         if(tau .lt. 0.01) then
            alp = 0.5 + tau/12 - tau**3/720
         else
            alp = (1-(1-exp(-tau))/tau)/(1-exp(-tau))
         endif
         var(j) = y0(j) + dt*(alp*w1(j)+(1-alp)*w0(j)-
     .        0.5*(a0(j)+a1(j))*y0(j))/(1+alp*tau)
      enddo
C
C  Black magic: oscillating vs evolving mode
C
      if((var(1)-y0(1))*(ys(1)-y0(1)) .lt. 0.0) then
         nCor = 4
         rCor = 0.2
      endif

      if(iCor .lt. nCor) then
         do j=1,nEqs
            y1(j) = (1-rCor)*y1(j) + rCor*var(j)
         enddo
         goto 20
      endif

      do j=1,nEqs
         rs(j) = var(j) - ys(j)
      enddo

      return
      end
#endif
      subroutine frtDummyCoolingSolvers(i)
      i = 0
      end
