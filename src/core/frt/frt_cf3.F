#include "frt_config.h"      
#ifdef RADIATIVE_TRANSFER
      subroutine frtInitCF(io,filename)
      include 'frt_cf3.inc'
      character*(*) filename

      open(unit=io, file=filename, status='old',
     .     form='unformatted')

      read(io) lt, ld, lp1, lp2, lp3, qstp(1), qstp(2), qstp(3),
     .     qmin(1), qmin(2), qmin(3), qmax(1), qmax(2), qmax(3)

      if(lt.ne.NT .or. ld.ne.ND .or. lp1.ne.NP1 .or. 
     .     lp2.ne.NP2 .or. lp3.ne.NP3) then
         write(0,*) 'RT::InitCF: fatal error, corrupted table:'
         write(0,*) '> NT= in file: ', lt, ' in code: ', NT
         write(0,*) '> ND= in file: ', ld, ' in code: ', ND
         write(0,*) '> NP1= in file: ', lp1, ' in code: ', NP1
         write(0,*) '> NP2= in file: ', lp2, ' in code: ', NP2
         write(0,*) '> NP3= in file: ', lp3, ' in code: ', NP3
         close(io)
         call frtAbort
      endif

      np(1) = NP1
      np(2) = NP2
      np(3) = NP3

      do i=1,3
         if(abs(qmax(i)-qmin(i)-(np(i)-1)*qstp(i)) .gt. 
     .        0.1*qstp(i)) then
            write(0,*) 'RT::InitCF: fatal error, corrupted table:'
            write(0,*) '> np=', np(i), ' qmin=', qmin(i), 'qmax=', 
     .           qmax(i), ' qstp=', qstp(i)
            close(io)
            call frtAbort
         endif
      enddo

      if(NX .gt. 1) then
         read(io) lx, xmin, xmax
         if(lx.ne.NX .or. NN.ne.1) then
            write(0,*) 'RT::InitCF: fatal error, corrupted table:'
            write(0,*) '> NX= in file: ', lx, ' in code: ', NX
            write(0,*) '> NX= in file: ',  1, ' in code: ', NN
            close(io)
            call frtAbort
         endif
         xstp = (xmax-xmin)/(NX-1)
      endif

      read(io) (alt(i),i=1,NT)
      altmin = alt(1)
      altstp = alt(2) - alt(1)

      read(io) (((indx(i,j,k), i=1,NP1), j=1,NP2), k=1,NP3)

      do ix=1,NX
         do id=1,ND
            read(io) (((data(j,ic,it,id,ix),j=1,NN),ic=1,NC),it=1,NT)
         enddo
      enddo

      if(NX .gt. 1) then
         do ix=1,NX
C$OMP PARALLEL DO DEFAULT(NONE),
C$OMP+PRIVATE(id,it,ic,j)
C$OMP+SHARED(data,ix)
            do id=1,ND
               do it=1,NT
                  do ic=1,NC
                     do j=1,NN
                        data(j,ic,it,id,ix) = log(1.0e-37+abs(
     .                       data(j,ic,it,id,ix)))
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif
 
      close(io)

      return
      end
C
C  Decode the interpolated function
C
      subroutine frtCFDecode1(p,it,id,wl,ws,jl,ju,ql,qu)
      include 'frt_cf3.inc'
      dimension p(NC), id(8), wl(3), ws(3)
      dimension vl(8), vu(8)
C
      do ic=1,NC
C
C  Already turned to log at init
C
         do j=1,8
            vl(j) = data(1,ic,it,id(j),jl)
            vu(j) = data(1,ic,it,id(j),ju)
         enddo

         p(ic) = exp(
     .        ql*(wl(3)*(wl(2)*(wl(1)*vl(1)+ws(1)*vl(2))+
     .                   ws(2)*(wl(1)*vl(3)+ws(1)*vl(4))) +
     .            ws(3)*(wl(2)*(wl(1)*vl(5)+ws(1)*vl(6))+
     .                   ws(2)*(wl(1)*vl(7)+ws(1)*vl(8))))+
     .        qu*(wl(3)*(wl(2)*(wl(1)*vu(1)+ws(1)*vu(2))+
     .                   ws(2)*(wl(1)*vu(3)+ws(1)*vu(4))) +
     .            ws(3)*(wl(2)*(wl(1)*vu(5)+ws(1)*vu(6))+
     .                   ws(2)*(wl(1)*vu(7)+ws(1)*vu(8)))))

      enddo

      return
      end
C
C  Decode the fitted function
C
      subroutine frtCFDecode2(p,it,id,wl,ws,dl,ds)
      include 'frt_cf3.inc'
      dimension p(NC), id(8), wl(3), ws(3)
      dimension v(8)
C
      do ic=1,NC

         do j=1,8
            a = data(1,ic,it,id(j),1)
            x1 = data(2,ic,it,id(j),1)*ds**data(3,ic,it,id(j),1)
            x2 = data(4,ic,it,id(j),1)*ds**data(5,ic,it,id(j),1)

            v(j) = (a**2+a*x1+x1**2)/(1+x2+x2**2)
            if(data(6,ic,it,id(j),1) .gt. 0.0) then
               v(j) = log(1.0e-37+abs(v(j)+data(6,ic,it,id(j),1)*dl))
            else
               v(j) = log(1.0e-37+abs(v(j)+data(6,ic,it,id(j),1)*ds))
            endif
         enddo

         p(ic) = exp(
     .        wl(3)*(wl(2)*(wl(1)*v(1)+ws(1)*v(2))+
     .               ws(2)*(wl(1)*v(3)+ws(1)*v(4))) +
     .        ws(3)*(wl(2)*(wl(1)*v(5)+ws(1)*v(6))+
     .               ws(2)*(wl(1)*v(7)+ws(1)*v(8))))

      enddo

      return
      end
C
C  Fill in the array chf(2,nmax) with cooling and heating functions
C
      subroutine frtSetCF(Z,den,pRate,chf)
      include 'frt_tables.inc'
      include 'frt_cf3.inc'
      dimension pRate(*), chf(2,*)
C
      dimension q(3), il(3), is(3), wl(3), ws(3)
      dimension id(8), p(NC)
C
C  Convert from nb to nH from Cloudy models 
C
      dl = max(1.0e-10,den*(1.0-0.02*Z)/1.4)

      if(pRate(iptabCiLW) .gt. 0.0) then 
         q(1) = log10(1.0e-37+abs(pRate(iptabCiLW)/dl))
         q(2) = log10(1.0e-37+abs(pRate(iptabPiH1)/pRate(iptabCiLW)))
         q(3) = log10(1.0e-37+abs(pRate(iptabPiO8)/pRate(iptabCiLW)))
C
C  qmin, qstp, etc are boundries of cells, not their centers
C
         do j=1,3
            w = 0.5 + (q(j)-qmin(j))/qstp(j)
            il(j) = int(w) + 1
            if(w .gt. il(j)-0.5) then
               is(j) = il(j) + 1
            else
               is(j) = il(j) - 1
            endif
            ws(j) = abs(il(j)-0.5-w)
            wl(j) = 1 - ws(j)
            if(il(j) .lt. 1) il(j) = 1
            if(il(j) .gt. np(j)) il(j) = np(j)
            if(is(j) .lt. 1) is(j) = 1
            if(is(j) .gt. np(j)) is(j) = np(j)
         enddo
      else
         do j=1,3
            il(j) = 1
            is(j) = 2
            wl(j) = 1
            ws(j) = 0
         enddo
      endif
C
C  Do not forget C-to-F77 index conversion
C
      id(1) = 1 + indx(il(1),il(2),il(3))
      id(2) = 1 + indx(is(1),il(2),il(3))
      id(3) = 1 + indx(il(1),is(2),il(3))
      id(4) = 1 + indx(is(1),is(2),il(3))
      id(5) = 1 + indx(il(1),il(2),is(3))
      id(6) = 1 + indx(is(1),il(2),is(3))
      id(7) = 1 + indx(il(1),is(2),is(3))
      id(8) = 1 + indx(is(1),is(2),is(3))

      if(NX .gt. 1) then
C
C  Density interpolation is still CIC
C
         w = (log10(dl)-xmin)/xstp
         jl = int(w) + 1
         if(jl .lt.  1) jl = 1
         if(jl .ge. NX) jl = NX-1
         ju = jl + 1
         ql = max(0.0,min(1.0,jl-w))
         qu = 1.0 - ql
      else
         ds = min(1.0e6,dl)
      endif

      do it=1,NT

         if(NX .gt. 1) then
            call frtCFDecode1(p,it,id,wl,ws,jl,ju,ql,qu)
         else
            call frtCFDecode2(p,it,id,wl,ws,dl,ds)
         endif

         a0 = p(1)
         a1 = p(2)
         a2 = p(3)
         p(2) = 2*a1 - 0.5*a2 - 1.5*a0
         p(3) = 0.5*(a0+a2) - a1

         a0 = p(4)
         a1 = p(5)
         a2 = p(6)
         p(5) = 2*a1 - 0.5*a2 - 1.5*a0
         p(6) = 0.5*(a0+a2) - a1

         if(MetalPartOnly) then
            chf(1,it) = (Z*p(3)+p(2))*Z
            chf(2,it) = (Z*p(6)+p(5))*Z
         else
            chf(1,it) = (Z*p(3)+p(2))*Z + p(1)
            chf(2,it) = (Z*p(6)+p(5))*Z + p(4)
         endif

      enddo

      return
      end
C
C  Get the cooling and heating functions for given T from the
C  pre-computed array chf(2,nmax)
C
      subroutine frtGetCF(tem,chf,cfun,hfun)
      include 'frt_cf3.inc'
      dimension chf(2,*)

      call frtGetCFln(log(max(1.0,tem)),chf,cfun,hfun)

      return
      end
C
C  Get the cooling and heating functions for given ln(T) from the
C  pre-computed array chf(2,nmax)
C
      subroutine frtGetCFln(altem,chf,cfun,hfun)
      include 'frt_cf3.inc'
      dimension chf(2,*)

      x = altem/2.3025851  !! log10(T)

      il = int((x-altmin)/altstp*0.99999) + 1
      if(il .lt.  1) il = 1
      if(il .ge. NT) il = NT-1
      iu = il + 1
      ql = max(0.0,min(1.0,(alt(iu)-x)/altstp))
      qu = 1.0 - ql
      
      cfun = ql*chf(1,il) + qu*chf(1,iu)
      hfun = ql*chf(2,il) + qu*chf(2,iu)

      return
      end
C
C  Get the cooling and heating functions for given ln(T) from the
C  pre-computed array chf(2,nmax)
C
      subroutine frtGetCFlg(x,chf,cfun,hfun)
      include 'frt_cf3.inc'
      dimension chf(2,*)

      il = int((x-altmin)/altstp*0.99999) + 1
      if(il .lt.  1) il = 1
      if(il .ge. NT) il = NT-1
      iu = il + 1
      ql = max(0.0,min(1.0,(alt(iu)-x)/altstp))
      qu = 1.0 - ql
      
      cfun = ql*chf(1,il) + qu*chf(1,iu)
      hfun = ql*chf(2,il) + qu*chf(2,iu)

      return
      end
#endif
