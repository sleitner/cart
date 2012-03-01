#include "frt_config.h"
#ifdef RADIATIVE_TRANSFER
#ifdef RT_LWBANDS
C  
C  Block data....define common blocks
C  
      block data frt_LWBandsData
      include 'frt_lwbands.inc'
      integer vLy, vW, vH
      parameter(vLy=18)
      dimension R0osc_str(0:vLy),R0diss(0:vLy),R0lam(0:vLy)
      dimension R1osc_str(0:vLy),R1diss(0:vLy),R1lam(0:vLy)
      dimension P1osc_str(0:vLy),P1diss(0:vLy),P1lam(0:vLy)
      common /RT_LW_data/R0osc_str,R0diss,R0lam,R1osc_str,R1diss,R1lam,
     + P1osc_str,P1diss,P1lam
      parameter(vW=5)
      dimension WR0osc_str(0:vW),WR0lam(0:vW),WR0diss(0:vW)
      dimension WQ1osc_str(0:vW),WQ1lam(0:vW),WQ1diss(0:vW)
      dimension WR1osc_str(0:vW),WR1lam(0:vW),WR1diss(0:vW)
      common /RT_LW_dataW/WR0osc_str,WR0diss,WR0lam,WR1osc_str,
     +     WR1diss,WR1lam,WQ1osc_str,WQ1diss,WQ1lam
      parameter(vH=14)
      dimension LyHosc_str(0:vH),LyHlam(0:vH)
      common /RT_LW_dataH/LyHlam,LyHosc_str

      data R0lam/1108.1303,1092.2,1077.142,1062.8825,1049.3674,1036.546,
     +1024.3759,1012.8174,1001.8284,991.3838,981.4426,971.9887,962.9799,
     + 954.4168, 946.1735, 938.4711, 931.0648, 923.9848, 917.2499/
      data R1lam/1108.6359,1092.7373,1077.7021,1063.4604,1049.9596,
     + 1037.1502,1024.99,1013.4412,1002.4565,992.0209,982.078,972.6347,
     + 963.6102,955.0691,946.9807,939.1275,931.7344,924.6434,917.9195/
      data P1lam/1110.066,1094.057,1078.929,1064.606,1051.033,1038.158,
     + 1025.9376,1014.332,1003.3015,992.815,982.839,973.347,964.3125,
     + 955.712, 947.5179, 939.7106, 932.2694, 925.1749, 918.4116/
c ... Lyb 1025.7223 .723 .7219 P1(6)
c ... Lyc 972.5368 .5371 .5367 R1(11)
c ... Ly11 918.1294 R1(18)
c ... Ly12 917.1806 R0(18)
      
      data R0osc_str/.17e-2,.58e-2,.12e-1,.18e-1,.23e-1,.27e-1,.29e-1,
     + .30e-1,.27e-1,.26e-1,.21e-1,.20e-1,.13e-1,.14e-1,.14e-2,.94e-2,
     + .11e-1,.61e-2,.61e-2/
      data R1osc_str/.11e-2,.38e-2,.77e-2,.12e-1,.16e-1,.18e-1,.20e-1,
     + .20e-1,.18e-1,.18e-1,.14e-1,.14e-1,.70e-2,.98e-2,.19e-1,.65e-2,
     + .74e-2,.41e-2,.44e-2/ 
      data P1osc_str/.58e-3,.20e-2,.39e-2,.60e-2,.76e-2,.87e-2,.91e-2,
     + .90e-2,.85e-2,.77e-2,.68e-2,.59e-2,.60e-2,.42e-2,.35e-2,.29e-2,
     + .24e-2,.19e-2,.16e-2/

      data P1diss/1.580e-6,3.250e-7,1.183e-5,1.308e-4,1.325e-3,1.779e-2,
     + 3.123e-2,2.074e-1,3.016e-1,4.133e-1,4.103e-1,4.063e-1,5.233e-1,
     + 5.111e-1,5.397e-1,5.603e-1,5.559e-1,6.062e-1,6.111e-1/
      data R0diss/1.486e-6,3.478e-7,1.198e-5,1.506e-4,1.282e-3,1.861e-2,
     + 3.331e-2,2.054e-1,3.121e-1,4.125e-1,4.100e-1,4.071e-1,5.199e-1,
     + 5.097e-1,3.481e-1,5.596e-1,5.507e-1,6.064e-1,6.090e-1/
      data R1diss/1.470e-6,3.630e-7,1.225e-5,1.951e-4,1.222e-3,2.009e-2,
     + 3.895e-2,2.004e-1,3.331e-1,4.122e-1,4.096e-1,4.112e-1,5.084e-1,
     + 5.072e-1,4.497e-1,5.580e-1,5.461e-1,6.061e-1,6.057e-1/
      
      data WR0lam/1008.5511,985.6294,964.9772,946.4237,929.5350,
     +     914.3972/
      data WQ1lam/1009.7703,986.7939,966.0895,947.4208,930.5792,
     +     915.4011/
      data WR1lam/1008.4972,985.6398,965.0577,946.3868,929.6911,
     +     914.6071/

      data WR0osc_str/.44e-1,.69e-1,.69e-1,.62e-1,.34e-1,.24e-1/
      data WQ1osc_str/.24e-1,.36e-1,.35e-1,.27e-1,.19e-1,.13e-1/
      data WR1osc_str/.21e-1,.33e-1,.35e-1,.13e-1,.16e-1,.12e-1/
      data WR0diss/1.353e-4,3.891e-4,1.599e-3,1.396e-1,2.873e-3,
     +     1.334e-3/
      data WQ1diss/1.218e-9,2.197e-10,3.205e-10,1.838e-8,2.277e-6,
     +     8.583e-6/
      data WR1diss/3.680e-4,1.278e-4,8.201e-3,5.735e-2,5.284e-3,
     +     3.713e-3/

      data LyHlam/1025.7223,1025.7230,1025.7219,972.5368,972.5371,
     +     972.5367, 949.7431,937.8035,930.7483,926.2257,923.1504,
     +     920.9631,919.3514,918.1294,917.1806/
      data LyHosc_str/7.912e-2,2.637e-2,5.274e-2,2.900e-2,9.668e-3,
     +     1.934e-2,1.934e-2,7.804e-3,4.816e-3,3.183e-3,2.216e-3,
     +     1.605e-3,1.200e-3,9.213e-4,7.226e-4/
      end
C
C
C
      subroutine frtLW_ReadData(path,JLWH2)
      include 'frt_lwbands.inc'
      integer n_lines_file
      parameter (n_lines_file=399)
      parameter (pi=3.14159265e0)
      dimension JLWH2(0:NRH2), tempJ(0:n_lines_file)
      integer jl, ii, i, IOUNIT
      parameter(IOUNIT = 99)
      character*(*) path

      XISH2=(lnu_b-lnu_a)/NRH2
      clight=3.e10   ![cm s-1]
      Atocm=1.e-8
      
      open(unit=IOUNIT, file=path//'/drt_lwbands.dat',
     .    status='old')
      do i=1,n_lines_file
         read(IOUNIT,*) lambd,ttpop3 ![erg s-1 A-1] 
         nu=clight/lambd/Atocm
         if (i.eq.1) lamin=lambd
         if (i.eq.n_lines_file) lamax=lambd
         tempJ(i-1)=ttpop3*lambd/nu/4./pi   !*1.e30 [erg s-1 Hz-1 sr-1]
      end do
      close(IOUNIT)
      tempJdla=(lamax-lamin)/(n_lines_file-1)

      
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(jl,lnu,nu,la,ii,q1,delta_la),
C$OMP+SHARED(XISH2,clight,Atocm,lamin,tempJdla,tempJ,JLWH2)
      do jl=0,NRH2
         
         lnu=lnu_a+jl*XISH2

         nu=exp(lnu)
         la=clight/nu/Atocm
         ii=int((la-lamin)/tempJdla+0.5)

         delta_la=la-(lamin+ii*tempJdla)
         q1=tempJ(ii)+(tempJ(ii+1)-tempJ(ii))/tempJdla*delta_la
         JLWH2(jl)=q1

      end do

      fnor = JLWH2(NRH2/2)

C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(jl),
C$OMP+SHARED(fnor,JLWH2)
      do jl=0,NRH2
         JLWH2(jl)=JLWH2(jl)/fnor
      end do

      return
      end
C
C
C
      subroutine frtLW_InitRun()
      include 'frt_lwbands.inc'
      integer nlines, vLy, vW, vH
      integer il,of,i_line, i_Hline,tr
      parameter (pi=3.14159265e0)
      parameter (nlines=24)
      parameter (vLy=18)
      dimension R0osc_str(0:vLy),R0diss(0:vLy),R0lam(0:vLy)
      dimension R1osc_str(0:vLy),R1diss(0:vLy),R1lam(0:vLy)
      dimension P1osc_str(0:vLy),P1diss(0:vLy),P1lam(0:vLy)
      common /RT_LW_data/R0osc_str,R0diss,R0lam,R1osc_str,R1diss,R1lam,
     +  P1osc_str,P1diss,P1lam
      parameter(vW=5)
      dimension WR0osc_str(0:vW),WR0lam(0:vW),WR0diss(0:vW)
      dimension WQ1osc_str(0:vW),WQ1lam(0:vW),WQ1diss(0:vW)
      dimension WR1osc_str(0:vW),WR1lam(0:vW),WR1diss(0:vW)
      common /RT_LW_dataW/WR0osc_str,WR0diss,WR0lam,WR1osc_str,
     +     WR1diss,WR1lam,WQ1osc_str,WQ1diss,WQ1lam
      parameter(vH=14)
      dimension LyHosc_str(0:vH),LyHlam(0:vH)
      common /RT_LW_dataH/LyHlam,LyHosc_str

      dimension nu_l(0:nlines,3), fs(0:nlines,3),fdis(0:nlines,3),
     +     nu_lH(0:vH)
c ... lines
      dimension i_line(0:nlines,3), alph_i(0:nlines,3),
     +     i_Hline(0:vH), fsH(0:vH), aH(0:vH)
      common /RT_LW_RAD_TR/ alph_i, fs, fdis, fsH, aH, i_line, i_Hline
c ... lines

      XISH2=(lnu_b-lnu_a)/NRH2
c     lnu=lnu_a+i*XISH
      clight=3.e10   ![cm s-1]
      Atocm=1.e-8
      para=0.25
      orto=0.75
C ....nlines=18+6+15
c ....=54+18+15=72+15=87

      do il=0,18
c ... remember to multiply by the para-orto ratio

         nu_l(il,1)=clight/P1lam(il)/Atocm
         fs(il,1)=orto*0.02654*P1osc_str(il)/nu_l(il,1)  !cm^2
         fdis(il,1)=fs(il,1)*P1diss(il)                  !cm^2
         nu_l(il,2)=clight/R1lam(il)/Atocm                 
         fs(il,2)=orto*0.02654*R1osc_str(il)/nu_l(il,2)  !cm^2
         fdis(il,2)=fs(il,2)*R1diss(il)                  !cm^2
         nu_l(il,3)=clight/R0lam(il)/Atocm    
         fs(il,3)=para*0.02654*R0osc_str(il)/nu_l(il,3)  !cm^2
         fdis(il,3)=fs(il,3)*R0diss(il)                  !cm^2
 
      end do

      of=-19
      do il=19,nlines
c ... remember to multiply by the para-orto ratio

         nu_l(il,1)=clight/WQ1lam(il+of)/Atocm
         fs(il,1)=orto*0.02654*WQ1osc_str(il+of)/nu_l(il,1)  !cm^2
         fdis(il,1)=fs(il,1)*WQ1diss(il+of)                  !cm^2
         nu_l(il,2)=clight/WR0lam(il+of)/Atocm
         fs(il,2)=para*0.02654*WR0osc_str(il+of)/nu_l(il,2)  !cm^2
         fdis(il,2)=fs(il,2)*WR0diss(il+of)                  !cm^2
         nu_l(il,3)=clight/WR1lam(il+of)/Atocm
         fs(il,3)=orto*0.02654*WR1osc_str(il+of)/nu_l(il,3)  !cm^2
         fdis(il,3)=fs(il,3)*WR1diss(il+of)                  !cm^2

      end do

c ... A21=6.669e15*f_osc/lamb(A)^2
c ... A_value=0.0449*A21/nu_l/4.3e-7
      do il=0,14

      nu_lH(il)=clight/LyHlam(il)/Atocm
      fsH(il)=0.02654*LyHosc_str(il)/nu_lH(il)               !cm^2 
      A21=6.669e15*LyHosc_str(il)/LyHlam(il)**2
      aH(il)=A21/nu_lH(il)/(2.*pi)
      i_Hline(il)=int((log(nu_lH(il))-lnu_a)/XISH2 + 0.5)    
      end do

      do il=0,nlines
         do tr=1,3
            i_line(il,tr)=int((log(nu_l(il,tr))-lnu_a)/XISH2 + 0.5)
         end do
      end do

      return
      end
C
C
C
      subroutine frtLW_SigmaH2(dLoga, SIG)
      include 'frt_lwbands.inc'
      integer i_line, i_Hline,i1,i2, nlines, vH
      integer i,il,tr,delta_i,delta_i1
      dimension SIG(0:NRH2)
c ... lines    
      parameter (nlines=24, vH=14)
      dimension i_line(0:nlines,3),alph_i(0:nlines,3),fdis(0:nlines,3),
     +     i_Hline(0:vH), fsH(0:vH), aH(0:vH), fs(0:nlines,3)
      common /RT_LW_RAD_TR/ alph_i, fs, fdis, fsH, aH, i_line, i_Hline
c ... lines

      XISH2=(lnu_b-lnu_a)/NRH2
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(i),
C$OMP+SHARED(SIG)
      do i=0,NRH2
         SIG(i)=0.
      end do

      delta_i = int(dLoga/XISH2 + 0.5)
      delta_i1 = 0

      do il=0,nlines
         do tr=1,3
            
            i1=max(i_line(il,tr)-delta_i,0)
            i2=max(i_line(il,tr)-delta_i1,0)
            i1=min(i1,NRH2)
            i2=min(i2,NRH2)

            do i=i1,i2
               
               SIG(i)=SIG(i)+fdis(il,tr)/max(XISH2,dLoga)

            end do
         end do
      end do

      return
      end
C
C
C
      subroutine frtLW_ComputeLineShifts(HUBBLE, CCC)
      include 'frt_lwbands.inc'
      integer il, tr, i_line, i_Hline, nlines, vH
      parameter (nlines=24, vH=14)
c ... lines
      dimension i_line(0:nlines,3),alph_i(0:nlines,3),fdis(0:nlines,3),
     +     i_Hline(0:vH), fsH(0:vH), aH(0:vH), fs(0:nlines,3)
      common /RT_LW_RAD_TR/ alph_i, fs, fdis, fsH, aH, i_line, i_Hline
c ... lines

       do il=0,nlines
         do tr=1,3
            alph_i(il,tr)=exp(-CCC/HUBBLE*(1.-0.14)*fs(il,tr))
         end do
      end do

      return
      end
C
C
C
      subroutine frtLW_ComputeAbs(HUBBLE, CCCH, DELTA_t, STAR_t, ALPH)
      include 'frt_lwbands.inc'
      integer i_line, i_Hline,i1,i2, nlines, vH
      integer i,il,tr,delta_i,delta_i1,wing
      parameter (pi=3.14159265e0)
      dimension ALPH(0:NRH2)
c ... lines    
      parameter (nlines=24, vH=14)
      dimension i_line(0:nlines,3),alph_i(0:nlines,3),fdis(0:nlines,3),
     +     i_Hline(0:vH), fsH(0:vH), aH(0:vH), fs(0:nlines,3)
      common /RT_LW_RAD_TR/ alph_i, fs, fdis, fsH, aH, i_line, i_Hline
c ... lines

      XISH2=(lnu_b-lnu_a)/NRH2
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(i),
C$OMP+SHARED(alph)
      do i=0,NRH2
         alph(i)=1.
      end do

      delta_i = int(HUBBLE*DELTA_t/XISH2 + 0.5)
      delta_i1= int(HUBBLE*STAR_t /XISH2 + 0.5)
      wing= int( 4.e-4 / XISH2 + 0.5)

c      CCC=?*CRAD2/CSSTAR

      do il=0,nlines
         do tr=1,3
            
            i1=max(i_line(il,tr)-delta_i,0)
            i2=max(i_line(il,tr)-delta_i1,0)
            i1=min(i1,NRH2)
            i2=min(i2,NRH2)
	    
            do i=i1,i2
              
               alph(i)=alph(i)*alph_i(il,tr)

            end do
         end do
      end do

      do il=0,14
         ff=XISH2/aH(il)

         i1=max(i_Hline(il)-delta_i  - wing,0)
         i2=max(i_Hline(il)-delta_i1 + wing,0)
         i1=min(i1,NRH2)
         i2=min(i2,NRH2)

         do i=i1,i2
            y1=(i-i_Hline(il)-delta_i1)*ff
            y2=(i-i_Hline(il)+ delta_i)*ff
            p=CCCH/HUBBLE*fsH(il)*(atan(y2)-atan(y1))/pi

            if (p.gt.10.) then
               pex=0.
            else
               pex=exp(-p)
            end if
            alph(i)=alph(i)*pex

         end do
      end do

      return
      end
C
C
C      
      subroutine frtLW_ComputeSource(HUBBLE, CCCH, DELTA_t, JLWH2, BETA)
      include 'frt_lwbands.inc'
      integer i
      dimension ALPH(0:NRH2), BETA(0:NRH2), JLWH2(0:NRH2)

      XISH2=(lnu_b-lnu_a)/NRH2
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(i),
C$OMP+SHARED(BETA)
      do i=0,NRH2
      BETA(i)=0.
      end do
      
      STEP = XISH2/HUBBLE

      do STAR_t = 0., DELTA_t, STEP

         call frtLW_ComputeAbs(HUBBLE, CCCH, DELTA_t, STAR_t, ALPH)
      
C$OMP PARALLEL DO DEFAULT(NONE), 
C$OMP+PRIVATE(i),
C$OMP+SHARED(BETA,STEP,ALPH,DELTA_t,JLWH2)
         do i=0,NRH2
            BETA(i)=BETA(i)+STEP*ALPH(i)/DELTA_t*JLWH2(i)
         end do
      end do

      return
      end
#endif
#endif
      subroutine frtDummyLWBands(i)
      i = 0
      end
