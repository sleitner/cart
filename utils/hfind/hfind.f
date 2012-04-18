c     ------------------------------------------------------------
c
c     HALO FINDER (fall 2002)
c 
c     Finds virial radius correctly. Removes unbound particles.
c
c     Escape velocity at a given halo centric radius is computed via
c     analytical expression for a NFW halo
c
c     Andrey Kravtsov 
c     
c     Updated: April 2007 
c
c     * reads in smoothed binary density file
c     * reads in both hydro and nbody particle formats
c     * writes new hp format which includes binding energy
c     * code cleanups
c     * added to cart utilities
c
c     ------------------------------------------------------------

      include 'hfind.h'
      character*40 FileName
      character*5 aname
      character*256 fname1, fname2, fname3, datpath, denpath, respath
      integer ldpath, ldenpath, lrpath
      common / abin / alpha
      character*256 ctmp
      character*10 chaexp
      character*250 froot

      alpha = 1.5
C
      if (IARGC().lt.4) then
        write (*,*)
     &'usage: exe aexp[e.g. 0.999] Dmin Dvir r_find[/h comoving kpc]'
        stop
      else
        call GETARG(1,ctmp)
        read (ctmp,*) chaexp
        call GETARG(2,ctmp)
        read (ctmp,*) DeltaMin
        call GETARG(3,ctmp)
        read (ctmp,*) Deltavir
        call GETARG(4,ctmp)
        read (ctmp,*) rfind

        msort = 0
        froot = '.'

        do iarg=5,iargc()

           call getarg(iarg,ctmp)

           if(ctmp .eq. '-vmax') then
              msort = 1
           else if(ctmp .eq. '-mvir') then
              msort = 2
           else if(ctmp(1:6) .eq. '-root=') then
              ib = index(ctmp,' ') - 1
              froot = ctmp(7:ib)
           else 
              write(0,*) 'Invalid option: ', ctmp
              stop
           endif

        enddo

      end if

c
c.... read control data
c
      open ( 15 , file = 'hfind.inp', status = 'old' )
      read  ( unit = 15 , fmt = '(10x,    a256)')  datpath
      read  ( unit = 15 , fmt = '(10x,    a256)')  denpath
      read  ( unit = 15 , fmt = '(10x,    a256)')  respath
      read  ( unit = 15 , fmt = '(70x,   e10.2)')  boxh
      read  ( unit = 15 , fmt = '(70x,     i10)')  ibtype
      read  ( unit = 15 , fmt = '(70x,   e10.2)')  rmin
      read  ( unit = 15 , fmt = '(70x,   e10.2)')  rmax
      read  ( unit = 15 , fmt = '(70x,     i10)')  nbins
      read  ( unit = 15 , fmt = '(70x,     i10)')  niter
      read  ( unit = 15 , fmt = '(70x,   e10.2)')  rvel
      read  ( unit = 15 , fmt = '(70x,     i10)')  nmin      
      close ( 15 ) 
c
c.... write it out for control
c
      write  ( unit = * , fmt ='(/'' data path:'',/A)')  datpath
      write  ( unit = * , fmt ='(/'' density path:'',/A)') denpath
      write  ( unit = * , fmt ='(/'' results path:'',/A)') respath
      write  ( unit = * , fmt = '( 1x,''Lbox:   '', 1pe10.2)') boxh
      write  ( unit = * , fmt = '( 1x,''ibtype:   '', i2,4x,
     &               ''niter:  '', i2)')  ibtype, niter
      write  ( unit = * , fmt = '( 1x,''rmin:   '', 1pe10.2,4x,
     &               ''rmax:   '', 1pe10.2,4x, 
     &               ''nbins:  '', i4)')  rmin, rmax, nbins
      write  ( unit = * , fmt = '( 1x,''rvel:   '', 1pe10.2,4x,
     &               ''nmin:  '', i4)') rvel, nmin

      ldpath = index( datpath , ' ' ) - 1 
      ldenpath = index( denpath , ' ' ) - 1 
      ldaexp = index( chaexp, ' ' ) - 1

      if(froot .eq. '.') then
         fname1 = datpath(1:ldpath)//
     .        '/PMcrda'//chaexp(1:ldaexp)//'.DAT '
         fname2 = datpath(1:ldpath)//
     .        '/PMcrs0a'//chaexp(1:ldaexp)//'.DAT '
      else
         ib = index(froot,' ') - 1
         fname1 = datpath(1:ldpath)//'/'//froot(1:ib)//'_a'//
     .        chaexp(1:ldaexp)//'.dph '
         fname2 = datpath(1:ldpath)//'/'//froot(1:ib)//'_a'//
     .        chaexp(1:ldaexp)//'.dxv '
      endif

      fname3 = denpath(1:ldenpath)//'/rho_smooth_a'
     &                  //chaexp(1:ldaexp)//'.dat '

      call InitArrays    ()
      call Read_Particles_Binary ( fname1, fname2, fname3) 

      call Conversions   ( boxh )

      lrpath = index( respath , ' ' ) - 1 
      fname1 = respath(1:lrpath)//'/hlist_'//chaexp(1:ldaexp)//'.dat '
      fname2 = respath(1:lrpath)//'/hpro_'//chaexp(1:ldaexp)//'.dat '
      fname3 = respath(1:lrpath)//'/hp_'//chaexp(1:ldaexp)//'.dat '

      rhalo     = rfind * rkpc2g
      rmin      = rmin  * rkpc2g 
      rmax      = rmax  * rkpc2g 

      if ( ibtype .eq. 2 ) then 
        nbins = int(2.*(rmax/rmin-1.0)**(2./3.)-1.0) 
      endif
      if ( nbins .gt. nbmax ) then 
        write(*,*) '* error: nbins =',nbins,' > nbmax in header: ',nbmax
        write(*,*) '* increase nbmax or decrease nbins and rerun...'
        stop
      else
        write(*,*) 'will construct profiles with nbins=',nbins
      endif

      write(*,*) 'Constructing Linked List...'
      call LL_Construct  ()
      write(*,*) 'Done.'
      write(*,*) 'Sorting particles...'
      call SortParticles ()
      write(*,*) 'Done.'
      write(*,*) 'Finding Halos...'
c     Warning: using a small value of rhalo or a large value of 
c     nmin parameter of FindHaloes can create an unintended density 
c     threshold for halo centers - DHR 11/6/2010
      call FindHaloes ( rhalo , nmin )
c      call RemoveSmall   ( nmin )
      write(*,*) 'Done.'
c      write(*,*) 'Iterating halos...'
c      call IterateHaloCM ( rhalo , rminshift )      
c      write(*,*) 'Done.'
c      write(*,*) 'Eating...'
c      call Cannibalism   ()
c      write(*,*) 'Done.'
      write(*,*) 'Computing halo properties...'
      call Compute_Halo_Properties ( ibtype, rmin, rmax, rvel, 
     &                               nbins, niter )
      call Remove_Velocity_Duplicates ()

      if(msort .eq. 1) then
         do ih=1,nhalo
            vhmax(ih) = -vhmax(ih)
         enddo
         call indexx(nhalo,vhmax,idx)
         do ih=1,nhalo
            vhmax(ih) = -vhmax(ih)
         enddo
      else if(msort .eq. 2) then
         do ih=1,nhalo
            amh(ih) = -amh(ih)
         enddo
         call indexx(nhalo,amh,idx)
         do ih=1,nhalo
            amh(ih) = -amh(ih)
         enddo
      else
         do ih=1,nhalo
            idx(ih) = ih
         enddo
      endif

      write(*,*) 'Done.'
      write(*,*) 'Writing halo catalogs...'
      call Write_Halo_Catalogs ( fname1, fname2 ,
     &                           ibtype , rmin , rmax , nbins , nmin,
     &                           rvel, rhalo )
      call Write_Halo_Particles ( fname3, 
     &                            ibtype, rmin, rmax, rvel, 
     &                            nbins, nmin )
      write(*,*) 'Done.'

      STOP
      END

c     -----------------------------------------------------------------
      subroutine FindRBin ( rd, ibtype, rmind, rmaxd, drd, ibin ) 
c     -----------------------------------------------------------------
c
c     given radius rd, and bin type ibtype find the corresponding
c     bin number ibin
c
c
      integer ibtype, ibin
      real rd, rmind, rmaxd, drd
      common / abin / alpha
c
      if ( rd .le. rmind ) then 
        ibin = 0 
      else
        if ( ibtype .eq. 0 ) then 
          ibin = int((log10(rd)-rmind)/drd)+1
        elseif ( ibtype .eq. 1 ) then ! linear 
          ibin = int((rd-rmind)/drd) + 1
        elseif ( ibtype .eq. 2 ) then ! custom
          ibin = int(2.*(rd/rmind)**(1./alpha)-1)
        else
          write(*,*) '* error in FindRBin unknown bin type =',ibtype
          return
        endif
      endif

      return
      end
c
c     -----------------------------------------------------------------
      subroutine FindIBin (ibin, ibtype, rmind, rmaxd, 
     &                     drd, rd1, rm, rd2, volr, voli  ) 
c     -----------------------------------------------------------------
c     given a bin number and binning info, return left, middle, and
c     right edges of the bin: rd1, rm, rd2, as well as differential volume
c     of the shell corresponding to the radial bin and integral volume
c     voli = V(<rr)
c
      include 'hfind.h'
      parameter ( pi43 = 4.0 * pi / 3.0 ) 
      common / abin / alpha

      IF ( ibin .eq. 0 ) then 
        if ( ibtype .eq. 0 ) then 
          rd2 = 10.0**(rmind)
        else
          rd2 = rmind
        endif
        rd1 = 0
        rm  = rd2*0.5
        voli = pi43 * rd2**3 * rhoaver
        volr = voli 
      ELSE 
        if ( ibtype .eq. 0 ) then ! log10 bins
          rd1 = 10.0**(rmind + float(ibin-1)*drd)
          rm  = rd1 * 10.0**(drd*0.5)
          rd2 = 10.0**(rmind + float(ibin)*drd)
        elseif ( ibtype .eq. 1 ) then 
          rd1 = rmind + float(ibin-1)*drd
          rm  = rd1 + 0.5*drd 
          rd2 = rmind + float(ibin)*drd
        elseif ( ibtype .eq. 2 ) then ! custom 
          rd1 = rmind * (float(ibin-1)/2. + 1)**alpha 
          rm  = rmind * ((float(ibin)-0.5)/2. + 1)**alpha
          rd2 = rmind * (float(ibin)/2. + 1)**alpha
        else
          write(*,*) '* error in FindIBin unknown bin type =',ibtype
          return
        endif
        voli = pi43 * rd2**3 * rhoaver
        volr = pi43 * (rd2**3 - rd1**3) * rhoaver
      ENDIF
c
      return
      end

c     ------------------------------------------------------------------
      subroutine Compute_Halo_Properties( ibtype, rpmin, rpmax, rvel, 
     &                                    npbin, niter )
c     ------------------------------------------------------------------
c
c     ibtype = 0 - log10 bins, 1 - even bins
c     rpmin,max - min/max. radius for profile construction in grid units
c     rvel = radius within which to measure mean halo velocity 
c     niter = number of iterations for unbound particles 
c

      include 'hfind.h'
      parameter ( pi43 = 4.0 * pi / 3.0 ) 
      real toohot 
      common / TOOHOT1 / toohot
      parameter ( dstopslope = -0.5 ) ! stop if density profile becomes flatter than this slope
      integer ibtype, npbin, niter
      dimension rr(1000) , vr(1000) , er(1000), dir(1000)
      real di(nil:nbmax), dd(nil:nbmax), ddi(nil:nbmax)
      real*8 hmassi 
      real hmold
      real EPSILON

      EPSILON = 1e-4

      toohot = 2.0  ! should really be in hfind.inp

      if     ( ibtype .eq. 0 ) then 
        rmax   = rpmax
        rmin   = rpmin
        rlmax  = log10(rmax)
        rlmin  = log10(rmin)
        drl    = (rlmax - rlmin)/float(npbin)
        nbins  = int((rlmax - rlmin)/drl) + 1
        rmind  = rlmin
        rmaxd  = rlmax
        drd = drl 
      elseif ( ibtype .eq. 1 ) then 
        rmax   = rpmax
        rmin   = rpmin
        drbin  = (rmax - rmin)/float(npbin) ! bin width in grid units
        nbins = npbin
        rmind  = rmin
        rmaxd  = rmax
        drd = drbin
      elseif ( ibtype .eq. 2 ) then 
        rmax  = rpmax
        rmin  = rpmin 
        rmaxd = rpmax
        rmind = rpmin
        nbins = npbin
        drd = 0.0
      else
        write(*,*) 'Compute_Halo_Properties: unknown bin type:',ibtype
        write(*,*) 'exiting...'
        return
      endif
      
      if ( ibtype .eq. 0 ) then 
      write(*,*) 'ibtype =',ibtype
      write(*,*) 'rmin =',rmind,' rmax=',rmaxd,' dr =',drd
      do ic2 = 0 , nbins
        call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                  rd1, rm, rd2, volr, voli ) 
        write(*,88) ic2, 0, 
     &       rd1*rg2kpc, rm*rg2kpc, rd2*rg2kpc, volr, voli
      enddo
      endif

      if ( ibtype .eq. 1 ) then   
      write(*,*) 'ibtype =',ibtype
      write(*,*) 'rmin =',rmind,' rmax=',rmaxd,' dr =',drd
      do ic2 = 0 , nbins
        call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                  rd1, rm, rd2, volr, voli ) 
        write(*,88) ic2, rd1*rg2kpc, rm*rg2kpc, rd2*rg2kpc, volr, voli
      enddo
      endif

      if ( ibtype .eq. 2 ) then 
      write(*,*) 'ibtype =',ibtype
      write(*,*) 'rmin =',rmind,' rmax=',rmaxd,' dr =',drd
      do ic2 = 0 , nbins
        call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                  rd1, rm, rd2, volr, voli ) 
        call FindRBin ( rm , ibtype, rmind, rmaxd, drd,
     &                  ibin ) 

        write(*,88) ic2, ibin, 
     &     rd1*rg2kpc, rm*rg2kpc, rd2*rg2kpc, volr, voli
      enddo
      endif
 88   format(2(i3,1x),3(f8.3,1x),1x,2(e11.6,1x))     

      if ( nbins .gt. nbmax ) then 
        write(*,*) 'error : Compute_Halo_Properties : nbins > nbmax'
        write(*,*) 'nbins =',nbins,' nbmax =',nbmax
        stop
      endif

      rhalo   = rmax
      rhalo2  = rmax * rmax
      rmax2   = rmax * rmax
      dvirlog = log10(Deltavir)

      do ic1 = 1 , npmax 
        iSp(ic1) = nil 
      enddo

      nr = 5 ! # of points to use for linear regression
      nr2 = 1
      if ( nr .lt. 2 ) then 
        write(*,*) 'number of bins for linear regression too small!!!'
        write(*,*) 'nr = nbins/4 =',nr
        write(*,*) 'stopping...'
        stop
      endif

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP+PRIVATE(iter,nhpold,ic2,ic3,ddi,di,dd,ptot,ivir,xhalo,yhalo,zhalo)
C$OMP+PRIVATE(vxhalo,vyhalo,vzhalo,vxhnew,vyhnew,vzhnew,vmassnew)
C$OMP+PRIVATE(wdummy,ndummy,rhalo,rdvel,xd,yd,zd,hmold)
C$OMP+PRIVATE(imin,jmin,kmin,imax,jmax,kmax)
C$OMP+PRIVATE(i,ic,j,jc,k,kc,idm)
C$OMP+PRIVATE(vxd,vyd,vzd,diff_x,diff_y,diff_z,corr_x,corr_y,corr_z)
C$OMP+PRIVATE(rd,dv,vesc,resc,ibin,irmaxflag,dex)
C$OMP+PRIVATE(rnp,vmn,rd1,rd2,rm,volr,voli,rld1,rld2,dil1,dil2,rvir)
C$OMP+PRIVATE(rhdum,ntot,ph,rvmax,vmax,imaxflag,irsflag,nrg,er,rr,dir)
C$OMP+PRIVATE(vr,a,ae,b,be,csq,q,hmassi,rcirc,rcircd,pdum,ptotd)
C$OMP+PRIVATE(ph1,ph2,irflag,pl1,pl2,ah,bh), schedule(dynamic,10)
      do ic1 = 1 , nhalo
        if ( mod(ic1-1,100) .eq. 0 ) then 
          write(*,*) 'processing halo',ic1,'...',nhalo-ic1,' left'
        endif

c....   iterate to remove unbound particles 
c
        do iter = 0 , niter
          nhpold = nhp(ic1)
          hmold = amh(ic1)
          do ic2 = nil , nbins
            na   (ic2,ic1) = nil 
            pn   (ic2,ic1) = zero
            pnt  (ic2,ic1) = zero
            pvx  (ic2,ic1) = zero
            pvy  (ic2,ic1) = zero
            pvz  (ic2,ic1) = zero
            vcirc(ic2,ic1) = zero
            vrms (ic2,ic1) = zero
            di(ic2) = zero
            dd(ic2) = zero 
          enddo
          ptot = 0.
          ivir = 0 
          xhalo = xh(ic1) 
          yhalo = yh(ic1) 
          zhalo = zh(ic1) 

          
          vxhalo = vxh(ic1)
          vyhalo = vyh(ic1)
          vzhalo = vzh(ic1)
          vxhnew = 0.
          vyhnew = 0.
          vzhnew = 0.
          vmassnew = 0.0 
          wdummy = 0.0
          ndummy = nil                           ! # of halo particles
c
c....     define the radius for vhalo estimate
c
          if ( iter .eq. 0 ) then 
            rhalo = rmax
            rdvel = rh(ic1)/5.0
          else
            rhalo = min(2.*rh(ic1),rhvir(ic1))
            rdvel = rhmax(ic1) * rvel 
          endif
             
c....     convert from original grid to chaining mesh          
          xd    = xhalo - rhalo
          yd    = yhalo - rhalo
          zd    = zhalo - rhalo
          imin  = int(float(nll) * (xd - 1) / ngrid + 1)
          jmin  = int(float(nll) * (yd - 1) / ngrid + 1)
          kmin  = int(float(nll) * (zd - 1) / ngrid + 1)
          xd    = xhalo + rhalo
          yd    = yhalo + rhalo
          zd    = zhalo + rhalo
          imax  = int(float(nll) * (xd - 1) / ngrid + 1)
          jmax  = int(float(nll) * (yd - 1) / ngrid + 1)
          kmax  = int(float(nll) * (zd - 1) / ngrid + 1)
c....     sweep over neighbors checking periodic boundary conditions
          do i = imin , imax
            ic = i 
            if ( ic .lt. 1   ) ic = ic + nll
            if ( ic .gt. nll ) ic = ic - nll
            do j = jmin , jmax 
              jc = j 
              if ( jc .lt. 1   ) jc = jc + nll
              if ( jc .gt. nll ) jc = jc - nll
              do k = kmin , kmax               
                kc = k 
                if ( kc .lt. 1   ) kc = kc + nll
                if ( kc .gt. nll ) kc = kc - nll 
                idm = iCL(ic,jc,kc)                 ! read LL head
                do while ( idm .ne. nil )
                  xd = xp(idm)
                  yd = yp(idm)
                  zd = zp(idm)
                  vxd = vxp(idm)
                  vyd = vyp(idm)
                  vzd = vzp(idm)
c....             enforce periodic boundary conditions
                  diff_x = xd - xhalo                  
                  diff_y = yd - yhalo 
                  diff_z = zd - zhalo 
                  corr_x = zero 
                  corr_y = zero 
                  corr_z = zero 

                  if ( abs(diff_x) .gt. ngrid/2 ) then 
                    if ( diff_x .gt. zero ) then 
                      corr_x = -ngrid 
                    else
                      corr_x =  ngrid
                    endif
                  endif
                  if ( abs(diff_y) .gt. ngrid/2 ) then 
                    if ( diff_y .gt. zero ) then 
                      corr_y = -ngrid 
                    else
                      corr_y =  ngrid
                    endif
                  endif
                  if ( abs(diff_z) .gt. ngrid/2 ) then 
                    if ( diff_z .gt. zero ) then 
                      corr_z = -ngrid 
                    else
                      corr_z =  ngrid
                    endif
                  endif

c....             compute distance particle - halo center squared
c                 corr_x , corr_y , corr_z - take care of periodicity
                  rd  = (diff_x + corr_x)**2 + 
     &                  (diff_y + corr_y)**2 + 
     &                  (diff_z + corr_z)**2   
c
c....             count particle if it falls inside halo radius and v < vesc 
c                   
                  if ( rd .lt. rmax2 ) then 
 
                    rd = max(sqrt(rd),1.e-7)

                    if ( iter .eq. 0 ) then 
                      dv  = 0.       
                      vesc = 1.e10
                    else
                      dv = sqrt( (vxd - vxhalo)**2 +
     &                           (vyd - vyhalo)**2 + 
     &                           (vzd - vzhalo)**2   ) * vg2kms         
                      resc = rd * rg2pkpc
                      vesc = ve ( resc , rhmax(ic1) , vhmax(ic1) )
                      vesc = vesc * toohot * 1.15**(niter-iter)
                    endif

                    if ( dv .lt. vesc ) then             
                     call FindRBin(rd, ibtype, rmind, rmaxd, drd, ibin ) 
                     ndummy = ndummy + 1
                     wdummy = wdummy + pw(idm) 
                     if ( rd .lt. min(rdvel,rh(ic1)) ) then 
                       vxhnew = vxhnew + pw(idm)*vxp(idm)
                       vyhnew = vyhnew + pw(idm)*vyp(idm)
                       vzhnew = vzhnew + pw(idm)*vzp(idm)
                       vmassnew = vmassnew + pw(idm)
                     endif
                     na (ibin,ic1) =  na(ibin,ic1) + 1
                     pn (ibin,ic1) =  pn(ibin,ic1) + pw(idm) 
                     pvx(ibin,ic1) = pvx(ibin,ic1) + vxp(idm)*pw(idm) 
                     pvy(ibin,ic1) = pvy(ibin,ic1) + vyp(idm)*pw(idm)
                     pvz(ibin,ic1) = pvz(ibin,ic1) + vzp(idm)*pw(idm)
                     vrms(ibin,ic1) = vrms(ibin,ic1) + 
     &                                  (vxp(idm)**2 + 
     &                                   vyp(idm)**2 + 
     &                                   vzp(idm)**2)*pw(idm)
                    endif
                  endif
                 idm = iLL(idm) ! next particle from linked list
                enddo  ! end do while 
              enddo  ! end k 
            enddo  ! end j 
          enddo  ! end i
c
c....     compute halo velocity 
c          
          if ( vmassnew .gt. 0.0 ) then 
            vxh(ic1) = vxhnew / vmassnew 
            vyh(ic1) = vyhnew / vmassnew 
            vzh(ic1) = vzhnew / vmassnew 
          else
            vxh(ic1) = 0.
            vyh(ic1) = 0.
            vzh(ic1) = 0.
          endif
c
c....     compute velocity components, and rms of radial shells 
c
          ptot = 0.
          irmaxflag = 0 
          if ( iter .eq. 0 ) then 
            rh(ic1) = rmax          
            rhvir(ic1) = rmax
          endif

          do ic2 = 0 , nbins 
            if ( na(ic2,ic1) .gt. 0 ) then 
              rnp = 1.0/pn(ic2,ic1)
              ptot = ptot + pn(ic2,ic1)
              pvx(ic2,ic1) = pvx(ic2,ic1) * rnp
              pvy(ic2,ic1) = pvy(ic2,ic1) * rnp
              pvz(ic2,ic1) = pvz(ic2,ic1) * rnp         
              vmn = pvx(ic2,ic1)**2 + pvy(ic2,ic1)**2 + pvz(ic2,ic1)**2
              vrms(ic2,ic1) = sqrt(abs(vrms(ic2,ic1)*rnp - vmn))
            endif
            pnt(ic2,ic1) = ptot
            call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rd2, volr, voli ) 
            di(ic2) = ptot / voli            
            dd(ic2) = pn(ic2,ic1) / volr
c
c....       estimate the virial radius 
c                      

            if ( ic2 .gt. 0 .and. ivir .eq. 0 .and. 
     &         di(ic2) .lt. Deltavir .and.
     &         di(ic2-1) .ge. Deltavir .and. iter .eq. 0 ) then
              ivir = ic2
              rld1 = log10(rd1)
              rld2 = log10(rd2)
              if ( di(ic2-1) .gt. 0.0 ) then 
                dil1 = log10(di(ic2-1))
              else
                dil1 = -3.0 
              endif
              if ( di(ic2) .gt. 0.0 ) then 
                dil2 = log10(di(ic2))
              else
                dil2 = -3.0 
              endif
              rvir = 10.0**(
     &                      (dvirlog * (rld2 - rld1) + rld1*dil2 - 
     &                       rld2*dil1) / (dil2 - dil1)
     &                     ) ! interpolate linearly in log r - log delta
              rhvir(ic1) = rvir 
            endif
            if ( iter .gt. 0 ) then 
              dex = di(ic2)/max((ddi(ic2)-di(ic2)),1.e-10)
            endif
            if ( iter .eq. 0 ) then 
              ddi(ic2) = di(ic2)
              rh(ic1) = rhvir(ic1)
            elseif (       (ddi(ic2) .gt. 0.0)
     &               .and. dex .le. 2.0 
     &               .and. irmaxflag .eq. 0 ) then 
              irmaxflag = 1 
              rh(ic1) = rd2
            endif
          enddo
c
c....     compute halo's extent
c
          rh(ic1) = min(rh(ic1),rhvir(ic1))
          if ( iter .eq. 0 ) then 
            rhdum = rhvir(ic1)
          else
            rhdum = rh(ic1)
          endif
c
c....     compute components of halo's velocity 
c
          vxhalo = 0.
          vyhalo = 0.
          vzhalo = 0.
          ntot = 0 
          ph = 0.
          rd = zero
          do ic2 = 0 , nbins
            call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rd, volr, voli ) 
c            if ( rhdum .gt. rd ) then 
            if ( rhdum .le. rd ) then
              ph = ph + pn(ic2,ic1)
            endif
          enddo
          if ( ph .gt. 0 ) then 
            nhp(ic1) = int(ph/pw(1))
            amh(ic1) = ph
          else
            nhp(ic1) = 0
            amh(ic1) = 0.
          endif
c
c....     construct circular velocity profile 
c 
          rd   = zero
          vmax = -1000.
          imaxflag = 0 
          irmaxflag = 0 
          irsflag = 0.0
          if ( iter .eq. 0 ) then 
            rsh(ic1) = rh(ic1)
          else
            rsh(ic1) = rhmax(ic1)/2.2
          endif

          do ic2 = 0 , nbins 
            call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rd, volr, voli ) 
            ptot = pnt(ic2,ic1)
            ptotd = ptot - pn(ic2,ic1)
            
            if ( irmaxflag .eq. 0 .and. rd .le. rh(ic1) ) then 
              if ( ic2 .lt. nbins-nr ) then 
                nrg = 0
                do ic3 = 0 , nr-1
                  nrg = nrg + 1
                  call FindIBin ( ic2+ic3 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rr(nrg), volr, voli ) 
                  er(nrg) = 1.0
                  rr(nrg) = log10(rr(nrg) * rg2pkpc)
                  if ( dd(ic2+ic3) .gt. 0. ) then 
                    dir(nrg) = log10(dd(ic2+ic3))
                  else
                    dir(nrg) = -3.
                  endif
                enddo
                call fit ( rr, dir, nrg, er, 0, a, b, ae, be, csq, q )
                if ( b-be .ge. dstopslope ) then 
                  rh(ic1) = rd
                  irmaxflag = 1
                endif
                if ( irsflag .eq. 0 ) then 
                  if ( b-be .le. -2.0 ) then 
                    rsh(ic1) = 10.0**(rr(nr2))
                    irsflag = 1
                  endif
                endif  
              endif 
            endif
c
c....  now work on vmax
c
            hmassi     = ptot * pmmsun
            rcirc      = rd * rg2pkpc * akpc
            vcirc(ic2,ic1) = 1.d-5 * dsqrt(G * hmassi * solarm/rcirc)
            if ( imaxflag .eq. 0 ) then 
              if ( ic2 .lt. nbins-nr ) then 
                nrg = 0
                do ic3 = 0 , nr-1
                  nrg = nrg + 1
                  call FindIBin ( ic2+ic3 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rr(nrg), volr, voli ) 
                  er(nrg) = 1.0
                  ptotd   = ptotd + pn(ic2+ic3,ic1)
                  hmassi  = ptotd * pmmsun
                  rr(nrg) = rr(nrg) * rg2pkpc ! r -> (physical) kpc
                  rcirc   = rr(nrg) * akpc ! r -> cm            
                  vr(nrg) = 1.d-5 * dsqrt(G * hmassi * solarm/rcirc) 
                  if ( vr(nrg) .gt. vmax ) then
                    vmax = vr(nrg)
                    rvmax = rr(nrg) 
                  endif
                enddo
                call fit ( rr, vr, nr, er, 0, a, b, ae, be, csq, q )
                if ( rd .le. rh(ic1) ) then
                  rhmax(ic1) = rvmax
                  vhmax(ic1) = vmax
                  if ( b-be .lt. 0 ) then 
                    imaxflag = 1
                  endif
                else                  
                  rcircd = min(2.16*rsh(ic1)/rg2pkpc,rh(ic1))
                  rcirc = rcircd*rg2pkpc
                  ic3 = 0 
                  pdum = 0.0 
                  do while ( pdum .eq. 0.0 .and. ic3 .le. nbins ) 
                    call FindIBin ( ic3 , ibtype, rmind, rmaxd, drd,
     &                      rd1, rm, rd2, volr, voli ) 
                    if ( rcircd .gt. rd1 .and. rcircd .le. rd2 ) then
                      pdum = pnt(ic3,ic1)
                      rcirc = rd2 * rg2pkpc
                    endif                        
                    ic3 = ic3 + 1
                  enddo        
                  rhmax(ic1) = rcirc
                  rcirc = rcirc * akpc
                  hmassi = pdum*pmmsun
                  vhmax(ic1) = 1.d-5 * dsqrt(G * hmassi * solarm/rcirc)
                  imaxflag = 1                   
                endif
              endif 
            endif
          enddo ! ic2
c
c....     compute halo radius, mass and the # of particles it contains
c
          rhdum = min(rh(ic1),rhvir(ic1))         
          if ( iter .gt. 0 ) then 
            rd = 0.
            ph1 = 0.      
            ph2 = 0.
            irflag = 0 
            do ic2 = 0 , nbins
              ph2 = ph2 + pn(ic2,ic1)
              call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                        rd1, rm, rd2, volr, voli ) 
              if ( rd1 .le. rhdum .and. rd2 .ge. rhdum ) then 
                if ( ic2 .eq. 0 ) then 
                  ph = pn(ic2,ic1)
                else
                  rld1 = log10(rd1)
                  rld2 = log10(rd2)
                  if ( ph1 .gt. 0.0 ) then 
                    pl1 = log10(ph1)
                  else
                    pl1 = -15.0 
                  endif
                  if ( ph2 .gt. 0.0 ) then 
                    pl2 = log10(ph2)
                  else
                    pl2 = -15.0
                  endif
                  ah = (pl2 - pl1)/(rld2 - rld1) 
                  bh = pl1 - ah * rld1
                  ph = 10.0**(ah*log10(rhdum) + bh) ! interpolate linearly in log r - log M
                endif
              endif            
              ph1 = ph2 
            enddo
            nhp(ic1) = int(ph/pw(1))
            amh(ic1) = ph
          endif
c
          if ( (nhp(ic1) .eq. 0) .or. 
     &         (abs((amh(ic1)/hmold - 1.0)) .le. EPSILON) ) goto 20
        enddo ! iter

 20     continue

      enddo ! end ic1

      return
      end

c     ------------------------------------------------------------------
      subroutine Write_Halo_Catalogs ( fname1, fname2, 
     &                               ibtype, rpmin, rpmax, npbin, nmin,
     &                               rvel, rhalo )
c     ------------------------------------------------------------------
      include 'hfind.h'
      parameter ( pi43 = pi * 4.0 / 3.0 ) 
      real toohot 
      common / TOOHOT1 / toohot
      real dvdup
      common / VDUP / dvdup 
      integer ibtype, npbin
      character*256 fname1, fname2

      nfn1 = index ( fname1 , ' ' ) - 1
      nfn2 = index ( fname2 , ' ' ) - 1 

      if     ( ibtype .eq. 0 ) then 
        rmax   = rpmax
        rmin   = rpmin
        rmaxd  = log10(rmax)
        rmind  = log10(rmin)
        drd    = (rlmax - rlmin)/float(npbin)
        nbins  = int((rlmax - rlmin)/drd) + 1
      elseif ( ibtype .eq. 1 ) then 
        rmaxd   = rpmax
        rmind   = rpmin
        drd  = drb 
        nbins = npbin
      elseif ( ibtype .eq. 2 ) then 
        rmaxd = rpmax
        rmind = rpmin
        nbins = npbin
      else
        write(*,*) 'Compute_Halo_Properties: unknown bin type:',ibtype
        write(*,*) 'exiting...'
        return
      endif
      if ( nbins .gt. nbmax ) then 
        write(*,*) 'error : Compute_Halo_Properties : nbins > nbmax'
        write(*,*) 'nbins =',nbins,' nbmax =',nbmax
        stop
      endif

      open ( 33 , file = fname1(1:nfn1) )
      open ( 34 , file = fname2(1:nfn2) )
      write (33,105) HEADER,
     &               AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &               EKIN,
     &               NROWC,NGRID,Om0,Oml0,Omb0,hubble
      write (34,105) HEADER,
     &               AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &               EKIN,
     &               NROWC,NGRID,Om0,Oml0,Omb0,hubble

 105  format (' ',A45,/
     &        ' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &        ' I =',I4,' WEIGHT=',F8.3,' Ekin=',E12.3,/
     &        ' Nrow=',I4,' Ngrid=',I4,
     &        '  Omega_0=',F7.3,' OmLam_0=',F7.4,
     &        '  Omegab_0=',F7.3,' Hubble=',f7.3)
      pm1 = pw(1) * pmmsun
      write(33,106) box0,rg2Mpc, vg2kms
      write(33,107) pmmsun,pm1
      write(33,120) DeltaMin
      write(33,121) Deltavir
      write(33,122) rhalo*rg2kpc
      write(33,108) rpmin*rg2kpc, rpmax*rg2kpc, npbin+1, ibtype
      write(33,109) rvel, toohot,dvdup
      write(33,110) nmin
      write(33,111) 
      write(33,112) 
      write(33,113) 
      write(33,114) 
      write(33,115) 
      write(34,106) box0,rg2Mpc, vg2kms
      write(34,107) pmmsun,pm1
      write(34,120) DeltaMin
      write(34,121) Deltavir
      write(34,122) rhalo*rg2kpc
      write(34,108) rpmin*rg2kpc, rpmax*rg2kpc, npbin+1, ibtype
      write(34,109) rvel, toohot, dvdup
      write(34,110) nmin
      write(34,116) 
      write(34,112) 
      write(34,113) 
      write(34,114) 
      write(34,115) 
      write(34,117) 
      write(34,118) 
      write(34,119) 

 120  format(' Deltamin=',e10.3,' [min. overd. for density maxima]')
 121  format(' Deltavir=',e10.3,' [virial overdensity]')
 122  format(' rhalo=',e9.2,' [/h kpc comoving=search radius]')
 106  format(' Lbox =',f7.2,' [/h Mpc comoving];',
     &       ' rg2Mpc=',e9.4,' [grid units->comoving /h Mpc];',
     &       ' vg2kms=',e9.4,' [grid units-> km/s peculiar velocity]')
 107  format(' M0 =',e9.4,' [/h Msun mass unit];',
     &       ' M(1)=',e9.4,' [/h Msun = mass of the 1st specie]')
 108  format(' rmin=',f8.3,'; rmax=',f8.3,' [/h kpc comoving];',
     &       ' nbins=',i4,' ibtype=',i2)
 109  format(' rvel=',f5.3,' [r/rhmax]; toohot=',f6.2,' dvdup=',f6.2)
 110  format(' min. number of particles =',i6)
 111  format(' column description:')
 112  format(' index x y z [/h Mpc comoving] vx vy vz [km/s peculiar]')
 113  format(' min(rvir,rmax) min(rt,rvir,rmax) [/h kpc comoving]')
 114  format(' M(<rh) [/h Msun] Np(<rh)')
 115  format(' vmax [km/s] rmax r_s [/h kpc proper] Central particle')
 116  format(' column description for halo header:')
 117  format(' column description for halo profiles:')
 118  format(' rm rr [/h kpc comoving] d(<r) d(r) Vc(r) Vrms(r) [km/s]')
 119  format(' m(<r) [/h Msun] Np(r)=# of particles in the shell')
c
      ih = 0 
      do ic0 = 1 , nhalo
        ic1 = idx(ic0)
        if ( amh(ic1) .gt. nmin*pw(1) ) then
        ih = ih + 1
        rkpc   = rh(ic1) * rg2kpc
        amassh = pmmsun * amh(ic1)
        xmpc  = (xh(ic1)-1.) * rg2Mpc
        ympc  = (yh(ic1)-1.) * rg2Mpc
        zmpc  = (zh(ic1)-1.) * rg2Mpc
        vxkms = vxh(ic1) * vg2kms
        vykms = vyh(ic1) * vg2kms
        vzkms = vzh(ic1) * vg2kms
        rvir = rhvir(ic1) * rg2kpc

        write(33,94)  ih,xmpc ,  ympc ,  zmpc ,
     &                vxkms   ,  vykms   , vzkms , 
     &                rvir, rkpc, amassh, nhp(ic1), 
     &                vhmax(ic1), rhmax(ic1), rsh(ic1),hc(ic1)-1
        write(34,94)  ih,xmpc ,  ympc ,  zmpc ,
     &                vxkms   ,  vykms   , vzkms , 
     &                rvir, rkpc, amassh, nhp(ic1), 
     &                vhmax(ic1), rhmax(ic1), rsh(ic1),hc(ic1)-1
       
        rd = zero
        do ic2 = 0 , nbins
          call FindIBin ( ic2 , ibtype, rmind, rmaxd, drd,
     &                    rd1, rdm, rd2, volr, voli ) 
          rdkpc1 = rdm * rg2kpc
          rdkpc2 = rd2 * rg2kpc
          odi = pnt(ic2,ic1)/voli
          od = pn(ic2,ic1)/volr
          write(34,95) rdkpc1, rdkpc2, odi, od, vcirc(ic2,ic1),
     &       vrms(ic2,ic1)*vg2kms, pnt(ic2,ic1)*pmmsun, int(na(ic2,ic1))
        enddo
        endif
      enddo

 94   format(i6,3(1x,f9.5),3(1x,f8.2),2(1x,f7.2),
     &        2x,e10.3,1x,i8,1x,f8.3,2(1x,f8.3),1x,i10)
 95   format(2(f8.3,1x),2x,2(e13.6,1x),1x,
     &       2(1x,f9.3),2x,e13.6,2x,i10)

      close ( 33 ) 
      close ( 34 ) 

      return
      end

c     ------------------------------------------------------------------
      subroutine Write_Halo_Particles( fname, 
     &                                 ibtype, rpmin, rpmax, rvel, 
     &                                 npbin, nmin )
c     ------------------------------------------------------------------
c
c     ibtype = 0 - log10 bins, 1 - even bins
c     rpmin,max - min/max. radius for profile construction in grid units
c     rvel = radius within which to measure mean halo velocity 
c     niter = number of iterations for unbound particles 
c
      include 'hfind.h'
      character*256 fname
      parameter ( pi43 = 4.0 * pi / 3.0 ) 
      real toohot 
      common / TOOHOT1 / toohot
      integer ibtype, npbin, niter

      nhd = 0 
      do ic1 = 1 , nhalo 
        if ( amh(ic1) .gt. nmin*pw(1) ) then
          nhd = nhd + 1
        endif
      enddo

      nfn1 = index ( fname , ' ' ) - 1
      open (21 , file = fname(1:nfn1), form='unformatted' )
      write(21) aexpn
      write(21) nhd, np

      do ic1 = 1 , npmax 
        iSp(ic1) = nil 
      enddo

      ih = 0 
      do ic0 = 1 , nhalo 
        ic1 = idx(ic0)
        if ( amh(ic1) .gt. nmin*pw(1) ) then
          ih = ih + 1
          rhalo  = rh(ic1)
          rhmaxd = rhmax(ic1)/rg2pkpc
          inp = 0
          xhalo = xh(ic1) 
          yhalo = yh(ic1) 
          zhalo = zh(ic1) 

          vxhalo = vxh(ic1)
          vyhalo = vyh(ic1)
          vzhalo = vzh(ic1)
c....     convert from original grid to chaining mesh          
          xd    = xhalo - rhalo
          yd    = yhalo - rhalo
          zd    = zhalo - rhalo
          imin  = int(float(nll) * (xd - 1) / ngrid + 1)
          jmin  = int(float(nll) * (yd - 1) / ngrid + 1)
          kmin  = int(float(nll) * (zd - 1) / ngrid + 1)
          xd    = xhalo + rhalo
          yd    = yhalo + rhalo
          zd    = zhalo + rhalo
          imax  = int(float(nll) * (xd - 1) / ngrid + 1)
          jmax  = int(float(nll) * (yd - 1) / ngrid + 1)
          kmax  = int(float(nll) * (zd - 1) / ngrid + 1)
c....     sweep over neighbors checking periodic boundary conditions
          do i = imin , imax
            ic = i 
            if ( ic .lt. 1   ) ic = ic + nll
            if ( ic .gt. nll ) ic = ic - nll
            do j = jmin , jmax 
              jc = j 
              if ( jc .lt. 1   ) jc = jc + nll
              if ( jc .gt. nll ) jc = jc - nll
              do k = kmin , kmax               
                kc = k 
                if ( kc .lt. 1   ) kc = kc + nll
                if ( kc .gt. nll ) kc = kc - nll 
                idm = iCL(ic,jc,kc)                 ! read LL head
                do while ( idm .ne. nil )
                  if ( iSp(idm) .eq. nil ) then 
                   xd = xp(idm)
                   yd = yp(idm)
                   zd = zp(idm)
                   vxd = vxp(idm)
                   vyd = vyp(idm)
                   vzd = vzp(idm)
c....              enforce periodic boundary conditions
                   diff_x = xd - xhalo                  
                   diff_y = yd - yhalo 
                   diff_z = zd - zhalo 
                   corr_x = zero 
                   corr_y = zero 
                   corr_z = zero 
 
                   if ( abs(diff_x) .gt. ngrid/2 ) then 
                     if ( diff_x .gt. zero ) then 
                       corr_x = -ngrid 
                     else
                       corr_x =  ngrid
                     endif
                   endif
                   if ( abs(diff_y) .gt. ngrid/2 ) then 
                     if ( diff_y .gt. zero ) then 
                       corr_y = -ngrid 
                     else
                       corr_y =  ngrid
                     endif
                   endif
                   if ( abs(diff_z) .gt. ngrid/2 ) then 
                     if ( diff_z .gt. zero ) then 
                       corr_z = -ngrid 
                     else
                       corr_z =  ngrid
                     endif
                   endif

c....              compute distance particle - halo center squared
c                  corr_x , corr_y , corr_z - take care of periodicity
                   rd  = sqrt( (diff_x + corr_x)**2 + 
     &                         (diff_y + corr_y)**2 + 
     &                         (diff_z + corr_z)**2   )   
                   rd = max(rd,1.e-7)

                   if ( rd .lt. rhalo ) then
                     resc = rd * rg2pkpc
                     dv = (vxd - vxhalo)**2 +
     &                    (vyd - vyhalo)**2 + 
     &                    (vzd - vzhalo)**2  
                     vesc = ve ( resc , rhmax(ic1) , vhmax(ic1) )
                     be = dv - vesc**2
                     dv = sqrt(dv) * vg2kms
                     vesc = vesc * toohot 
c
c....                count particle if it falls inside halo radius and v < vesc 
c                   
                     if ( dv .lt. vesc ) then
                       bind(inp) = be
                       ind(inp) = idm-1  
                       inp = inp + 1  
                     endif
                   endif
                  endif 
                  idm = iLL(idm) ! next particle from linked list
                enddo  ! end do while 
              enddo  ! end k 
            enddo  ! end j 
          enddo  ! end i

          write(21) ih, inp, (ind(i), i=1,inp), (bind(i), i=1,inp)
        endif
      enddo ! end ic1

      close ( 21 ) 

      return
      end
c
c     ------------------------
      subroutine InitArrays ()
c     ------------------------
c     purpose: to initialize working arrays
c     -------------------------------------
      include 'hfind.h'

      do ic1 = 1 , nll
        do ic2 = 1 , nll
          do ic3 = 1 , nll 
            iCL(ic3,ic2,ic1) = nil 
          enddo
        enddo
      enddo

      do ic1 = 1 , npmax 
        iLL(ic1) = nil 
        iSp(ic1) = nil 
      enddo

      do ic1 = 1 , nh 
        iSh(ic1) = nil 
        nhp(ic1) = nil 
        rhvir(ic1) = 0.0
        rh(ic1) = 0.0
        rsh(ic1) = 0.0 
        amh(ic1) = 0.0
      enddo

      return
      end

c     --------------------------------
      subroutine Conversions ( boxh0 ) 
c     --------------------------------
c     purpose: sets up various grid <-> physical 
c              conversion factors
c     input  : boxh0 - box size in h^{-1} Mpc at z=0
c     -------------------------------------------------------
      include 'hfind.h'
      real boxh0   
c
c.... particle mass in h^{-1} solar masses (CDM,LCDM)
c     divide by hubble to get physical mass
c
      pmmsun = 2.760368e11 * Om0 * 
     &         (boxh0)**3 / (ngridc**3) ! must be in /h Msun !
      fb = Omb0 / Om0 
      pm1 = pw(1) * pmmsun
      box0 = boxh0
c
c.... grid scale <-> physical scale conversions
c     using boxh will give units in h^{-1} 
c
      rMpc2g  = float(ngrid)/boxh0   ! conversion /h Mpc -> grid units
      rkpc2g  = rMpc2g / 1000.0     ! conversion /h kpc -> grid units
      rg2Mpc  = boxh0 / ngrid        ! conversion grid units -> /h Mpc
      rg2kpc  = 1000. * boxh0 / ngrid ! conversion grid units -> /h kpc
      rg2pMpc = boxh0 * aexpn / ngrid ! conversion grid units -> proper /h Mpc
      rg2pkpc = 1000. * boxh0 * aexpn / ngrid ! conversion grid units -> proper /h kpc
c
c.... grid velocity <-> km/s = (x0*Ho)/aexpn
c
      if ( hydro_flag .eq. 1 ) then
        vg2kms = 50.0*boxh0 / (aexpn*ngridc) * sqrt(Om0)
      else
        vg2kms = 100.0*boxh0 / (aexpn*ngridc)
      endif

      write(*,*) 'aexp =',aexpn
      
      write(*,*) 'Lbox =',boxh0,' [/h Mpc]'
      write(*,*) 'M0 =',pmmsun,' [/h Msun]'
      write(*,*) 'M(ispec=1) =',pm1,' [/h Msun]'
      write(*,*) 'pw(1) =',pw(1)
      write(*,*) 'rkpc2g =',rkpc2g
      write(*,*) 'vg2kms =',vg2kms
      return
      end

c     --------------------------
      subroutine LL_Construct ()
c     --------------------------
c     constructs particle linked list 
c     see description of algorithms in  
c       T.H.Cormen , C.E.Leiserson , R.L.Rivest
c       "Introduction to algorithms" pp. 204-208
c     ------------------------------------------
      include 'hfind.h'
      
c.... construct linked list
      do ic1 = 1 , np
        ip            = int(float(nll)*(xp(ic1)-1.0)/float(ngrid)+1)
        jp            = int(float(nll)*(yp(ic1)-1.0)/float(ngrid)+1)
        kp            = int(float(nll)*(zp(ic1)-1.0)/float(ngrid)+1)
        if ( ip .gt. nll ) ip = 1
        if ( jp .gt. nll ) jp = 1
        if ( kp .gt. nll ) kp = 1
        if ( ip .gt. nll .or. jp .gt. nll .or. kp .gt. nll .or.
     &       ip .le. 0 .or. jp .le. 0 .or. kp .le. 0 ) then
           write(*,*) 'ic1 = ', ic1
           write(*,*) 'x = ', xp(ic1), yp(ic1), zp(ic1)
           write(*,*) 'ip = ', ip, jp, kp
           write(*,*) 'nll = ', nll
           write(*,*) 'ngrid = ', ngrid
           stop
        endif
        iLL(ic1)      = iCL(ip,jp,kp)
        iCL(ip,jp,kp) = ic1        
      enddo
      return
      end


c     ---------------------------
      subroutine SortParticles ()
c     ---------------------------
c     sort particles according to their local density
c     ind(np) has index of particle with largest local density
c
c     uses: indexx subroutine from NR
c     -----------------------------------------------------
      include 'hfind.h'

      call indexx ( np , dnb , ind )

      return
      end

c     -------------------------------
      subroutine FindHaloes ( rhalo , npmin )
c     -------------------------------
c     purpose  : finds initial halo centers 
c
c     algorithm: particle with the largest local density (in dnb)
c                determines the approximate location of the first 
c                halo center. Particles located inside a sphere
c                of radius rhalo centered at the halo center are assigned 
c                to the same halo and are excluded from  the list of 
c                particles used to identify halos (iSp is set to 1). 
c                The radius rhalo is an adjustable parameter; its values
c                should be of order of spatial resolution of a simulation. 
c                The procedure repeats for the particle with the largest 
c                local density in the list of remaining particles. The peaks 
c                are identified until there are no particles in the list.
c     ----------------------------------------------------------------------
      include 'hfind.h'

      nhalo = nil       
c
c.... the particles are sorted so that ind(np) has index of particle 
c     with largest local density - this determines the inverse order 
c     in the loop below

      rhalo2  = rhalo**2

      do ic1 = np , 1 , -1
        if (mod(ic1,500000).eq.nil) then 
          write(*,*) ic1,' processed delta=', dnb(ind(ic1))
        endif
c....   determine center of the next halo
        if ( iSp(ind(ic1)) .eq. nil .and. 
     &       dnb(ind(ic1)) .gt. Deltamin ) then 

          iCenter      = ind(ic1)     
          xhalo        = xp(iCenter)
          yhalo        = yp(iCenter)
          zhalo        = zp(iCenter)
c....     convert from original grid to chaining mesh          
          ihalo        = int(float(nll) * (xhalo - 1) / ngrid + 1)
          jhalo        = int(float(nll) * (yhalo - 1) / ngrid + 1)
          khalo        = int(float(nll) * (zhalo - 1) / ngrid + 1)
          xd           = xhalo - rhalo
          yd           = yhalo - rhalo
          zd           = zhalo - rhalo
          imin         = int(float(nll) * (xd - 1) / ngrid + 1)
          jmin         = int(float(nll) * (yd - 1) / ngrid + 1)
          kmin         = int(float(nll) * (zd - 1) / ngrid + 1)
          xd           = xhalo + rhalo
          yd           = yhalo + rhalo
          zd           = zhalo + rhalo
          imax         = int(float(nll) * (xd - 1) / ngrid + 1)
          jmax         = int(float(nll) * (yd - 1) / ngrid + 1)
          kmax         = int(float(nll) * (zd - 1) / ngrid + 1)
c....     sweep over neighbors checking periodic boundary conditions
          xcm = 0.
          ycm = 0. 
          zcm = 0. 
          wdummy = nil                           ! # of halo particles
          n10    = nil 
          dmax = -1.e8
          ndummy = nil 
          do i = imin , imax
            ic = i 
            if ( ic .lt. 1   ) ic = ic + nll
            if ( ic .gt. nll ) ic = ic - nll
            do j = jmin , jmax 
              jc = j 
              if ( jc .lt. 1   ) jc = jc + nll
              if ( jc .gt. nll ) jc = jc - nll
              do k = kmin , kmax               
                kc = k 
                if ( kc .lt. 1   ) kc = kc + nll
                if ( kc .gt. nll ) kc = kc - nll 
                idummy = iCL(ic,jc,kc)                 ! read LL head
                do while ( idummy .ne. nil )
                  if ( iSp(idummy) .eq. nil ) then
                    xd = xp(idummy)
                    yd = yp(idummy)
                    zd = zp(idummy)
                    diff_x = xd - xhalo 
                    diff_y = yd - yhalo 
                    diff_z = zd - zhalo 
                    corr_x = zero 
                    corr_y = zero 
                    corr_z = zero 
                    if ( abs(diff_x) .gt. ngrid/2 ) then 
                      if ( diff_x .gt. zero ) then 
                        corr_x = -ngrid 
                      else
                        corr_x = ngrid
                      endif
                    endif
                    if ( abs(diff_y) .gt. ngrid/2 ) then 
                      if ( diff_y .gt. zero ) then 
                        corr_y = -ngrid 
                      else
                        corr_y = ngrid
                      endif
                    endif
                    if ( abs(diff_z) .gt. ngrid/2 ) then 
                      if ( diff_z .gt. zero ) then 
                        corr_z = -ngrid 
                      else
                        corr_z = ngrid
                      endif
                    endif

                    rd = (xd - xhalo + corr_x)**2 + 
     &                   (yd - yhalo + corr_y)**2 + 
     &                   (zd - zhalo + corr_z)**2 

                    if ( rd .le. rhalo2 ) then
                      ndummy = ndummy + 1
                      wdummy     = wdummy + pw(idummy)
                      xcm = xcm + xp(idummy)*pw(idummy)
                      ycm = ycm + yp(idummy)*pw(idummy)
                      zcm = zcm + zp(idummy)*pw(idummy)
                      dmax = max(dmax,dnb(idummy))
                      iSp(idummy) = iCenter                     
                    endif
                  endif
                  idummy = iLL(idummy)
                enddo  ! end do while 
              enddo  ! end k 
            enddo  ! end j 
          enddo  ! and i
          vol   = 4.0 * pi * rhalo**3 / 3.0
          overd = wdummy / vol / rhoaver

          if ( (wdummy.gt.0.) .and. 
     &          dnb(iCenter) .ge. dmax .and. 
     &          ndummy .gt. npmin ) then
            nhalo      = nhalo + 1
            if ( nhalo .gt. nh ) then 
              write(*,*) '* error: nh is too small: nh =',nhalo
              write(*,*) '* while more halos are found. stopping...'
              stop
            endif

            hc(nhalo)  = iCenter
            xh(nhalo)  = xhalo
            yh(nhalo)  = yhalo
            zh(nhalo)  = zhalo
            rh(nhalo)  = rhalo 
            nhp(nhalo) = int(wdummy/pw(1))
            amh(nhalo) = wdummy

            xcm = xcm / wdummy
            ycm = ycm / wdummy
            zcm = zcm / wdummy
          else 
c           clear particles marked for this center since it lacked
c           sufficient particles within rhalo
            do i = imin , imax
              ic = i 
              if ( ic .lt. 1   ) ic = ic + nll
              if ( ic .gt. nll ) ic = ic - nll
              do j = jmin , jmax 
                jc = j 
                if ( jc .lt. 1   ) jc = jc + nll
                if ( jc .gt. nll ) jc = jc - nll
                do k = kmin , kmax               
                  kc = k 
                  if ( kc .lt. 1   ) kc = kc + nll
                  if ( kc .gt. nll ) kc = kc - nll 
                  idummy = iCL(ic,jc,kc)                 ! read LL head
                  do while ( idummy .ne. nil )
                    if ( iSp(idummy) .eq. iCenter ) then
                      iSp(idummy) = nil
                    endif
                    idummy = iLL(idummy)
                  enddo  ! end do while 
                enddo  ! end k 
              enddo  ! end j 
            enddo  ! and i
          endif
        endif
      enddo

      write(*,*) ' nhalo =', nhalo 

      return
      end

c     ------------------------------------------
      subroutine IterateHaloCM ( r , rminshift )
c     ------------------------------------------
c     purpose: finds centers of mass for haloes iteratively by 
c              finding the center of mass of all particles inside r and 
c              displacing the center of the sphere to the center of mass. 
c              The procedure is iterated until convergence (displacement 
c              is less than rminshift or if number of particles inside r 
c              starts to decrease).
c     -----------------------------------------------------
      include 'hfind.h'
      
      real r , rminshift 
      real xt , yt , zt , xnew , ynew , znew , riter
      integer nsave(nh)

      do ic1 = 1 , nhalo 
        iSh(ic1) = nil 
         rh(ic1) = r
        nhp(ic1) = 0 
      enddo

      riter = r
      niter = nil 
      nstop = nil 
 7    iter  = nil 
      niter = niter + 1
      do ic1 = 1 , nhalo 
        xt = xh(ic1) 
        yt = yh(ic1) 
        zt = zh(ic1) 
        nt = nhp(ic1)
        it = iSh(ic1)
        
        if ( it .eq. nil ) then 
          call CM ( xt , yt , zt , xnew , ynew , znew , 
     &              vxnew , vynew , vznew , nnew , riter )
          rshift = max (abs(xnew-xt) , abs(ynew-yt) , abs(znew-zt))
          if ((rshift.lt.rminshift) .or. (nnew.le.nt)) then 
            iSh(ic1) = 1
          else
            iter = iter + 1
            if ( xnew .lt. 1.0 ) xnew = xnew + ngrid 
            if ( ynew .lt. 1.0 ) ynew = ynew + ngrid 
            if ( znew .lt. 1.0 ) znew = znew + ngrid 
            if ( xnew .gt. xn  ) xnew = xnew - ngrid 
            if ( ynew .gt. xn  ) ynew = ynew - ngrid 
            if ( znew .gt. xn  ) znew = znew - ngrid 
            xh(ic1)  = xnew 
            yh(ic1)  = ynew 
            zh(ic1)  = znew 
            if ( nnew .gt. 0 ) then 
              vxh(ic1) = vxnew 
              vyh(ic1) = vynew 
              vzh(ic1) = vznew 
              nhp(ic1) = nnew       
            endif
          endif
        endif
      enddo
      write(*,*) 'iterating halo mass centers ', iter, 'left'
      if ((iter .gt. 0).and.(niter.lt.100)) goto 7 
      
      return
      end

c     -------------------------------------------------------------
      subroutine CM ( xold , yold , zold , 
     &                xnew , ynew , znew , vxnew , vynew , vznew , 
     &                nnew , r )
c     -------------------------------------------------------------
      include 'hfind.h'
      real xold , yold , zold , xnew , ynew , znew , r
      
      r2 = r**2
      
      xnew   = 0.0
      ynew   = 0.0
      znew   = 0.0 
      vxnew  = 0.0 
      vynew  = 0.0 
      vznew  = 0.0 

      wdummy = 0.0

c.... convert from original grid to chaining mesh          
      xd    = xold - r
      yd    = yold - r
      zd    = zold - r
      imin  = int(float(nll) * (xd - 1) / ngrid + 1)
      jmin  = int(float(nll) * (yd - 1) / ngrid + 1)
      kmin  = int(float(nll) * (zd - 1) / ngrid + 1)
      xd    = xold + r
      yd    = yold + r
      zd    = zold + r
      imax  = int(float(nll) * (xd - 1) / ngrid + 1)
      jmax  = int(float(nll) * (yd - 1) / ngrid + 1)
      kmax  = int(float(nll) * (zd - 1) / ngrid + 1)
c.... sweep over neighbors checking periodic boundary conditions
      ndummy = nil                           ! # of halo particles
      do i = imin , imax
        ic = i 
        if ( ic .lt. 1   ) ic = ic + nll
        if ( ic .gt. nll ) ic = ic - nll
        do j = jmin , jmax 
          jc = j 
          if ( jc .lt. 1   ) jc = jc + nll
          if ( jc .gt. nll ) jc = jc - nll
          do k = kmin , kmax               
            kc = k 
            if ( kc .lt. 1   ) kc = kc + nll
            if ( kc .gt. nll ) kc = kc - nll 
            idummy = iCL(ic,jc,kc)                 ! read LL head
            do while ( idummy .ne. nil )
              xd = xp(idummy)
              yd = yp(idummy)
              zd = zp(idummy)
              diff_x = xd - xold
              diff_y = yd - yold
              diff_z = zd - zold
              corr_x = zero 
              corr_y = zero 
              corr_z = zero 
              if ( abs(diff_x) .gt. ngrid/2 ) then 
                if ( diff_x .gt. zero ) then 
                  corr_x = -ngrid 
                else
                  corr_x = ngrid
                endif
              endif
              if ( abs(diff_y) .gt. ngrid/2 ) then 
                if ( diff_y .gt. zero ) then 
                  corr_y = -ngrid 
                else
                  corr_y = ngrid
                endif
              endif
              if ( abs(diff_z) .gt. ngrid/2 ) then 
                if ( diff_z .gt. zero ) then 
                  corr_z = -ngrid 
                else
                  corr_z = ngrid
                endif
              endif

              rd = (xd - xold + corr_x)**2 + 
     &             (yd - yold + corr_y)**2 + 
     &             (zd - zold + corr_z)**2 

              if ( rd .le. r2 ) then
                 wdummy = wdummy + pw(idummy)
                 xnew   = xnew  +  xp(idummy)*pw(idummy)
                 ynew   = ynew  +  yp(idummy)*pw(idummy)
                 znew   = znew  +  zp(idummy)*pw(idummy)
                 vxnew  = vxnew + vxp(idummy)*pw(idummy)
                 vynew  = vynew + vyp(idummy)*pw(idummy)
                 vznew  = vznew + vzp(idummy)*pw(idummy)
              endif
              idummy = iLL(idummy)
            enddo  ! end do while 
          enddo  ! end k 
        enddo  ! end j 
      enddo  ! and i            
      
      if (wdummy .eq. 0.0 ) then 
        xnew  =  xold
        ynew  =  yold
        znew  =  zold
        nnew  = 0
      else
       xnew  =  xnew / wdummy 
       ynew  =  ynew / wdummy 
       znew  =  znew / wdummy 
       vxnew = vxnew / wdummy
       vynew = vynew / wdummy
       vznew = vznew / wdummy
       nnew  = int(wdummy/pw(1))
      endif
      

      return
      end

c     -------------------------------
      subroutine RemoveSmall ( nmin )
c     -------------------------------
c     purpose: removes small (with nhp < nmin) halos
c     input  : nmin - keep only halos with nhp > nmin 
c     ---------------------------------------------
      include 'hfind.h'
      integer nmin 

      nnew = nil 
      do ic1 = 1 , nhalo
        if ( nhp(ic1) .gt. nmin ) then 
          nnew      = nnew + 1
          call HOld2New ( ic1 , nnew ) 
        endif        
      enddo

      write(*,*) nhalo-nnew, ' small haloes were removed'

      nhalo = nnew

      return
      end

c     ----------------------------------------
      subroutine Remove_Velocity_Duplicates ()
c     ----------------------------------------
c     neighboring haloes eat each other 
c     ---------------------------------
      include 'hfind.h'
      common / VDUP / dvdup 
      dvdup = 0.1 
c
      do ic1 = 1 , nhalo 
        iSh(ic1) = nil 
      enddo 

      do ic1 = 1 , nhalo-1
        if ( iSh(ic1) .eq. nil ) then
          vxhd1 = vxh(ic1)
          vyhd1 = vyh(ic1)
          vzhd1 = vzh(ic1)
          vhdumi = 1.0 / sqrt(vxhd1**2 + vyhd1**2 + vzhd1**2)
          do ic2 = ic1+1 , nhalo 
            if ( iSh(ic2) .eq. nil ) then 
              vxhd2 = vxh(ic2)
              vyhd2 = vyh(ic2)
              vzhd2 = vzh(ic2)
              diff_x = xh(ic2) - xh(ic1) 
              diff_y = yh(ic2) - yh(ic1) 
              diff_z = zh(ic2) - zh(ic1)
              corr_x = zero 
              corr_y = zero 
              corr_z = zero 
              if ( abs(diff_x) .gt. ngrid/2 ) then 
                if ( diff_x .gt. zero ) then 
                  corr_x = -ngrid 
                else
                  corr_x = ngrid
                endif
              endif
              if ( abs(diff_y) .gt. ngrid/2 ) then 
                if ( diff_y .gt. zero ) then 
                  corr_y = -ngrid 
                else
                  corr_y = ngrid
                endif
              endif
              if ( abs(diff_z) .gt. ngrid/2 ) then 
                if ( diff_z .gt. zero ) then 
                  corr_z = -ngrid 
                else
                  corr_z = ngrid
                endif
              endif
              rd = (xh(ic2) - xh(ic1) + corr_x)**2 +
     &             (yh(ic2) - yh(ic1) + corr_y)**2 +
     &             (zh(ic2) - zh(ic1) + corr_z)**2 
              rhmaxd = min(rhmax(ic1),rhmax(ic2))
              
              dv = sqrt((vxhd1-vxhd2)**2 + 
     &                  (vyhd1-vyhd2)**2 + 
     &                  (vzhd1-vzhd2)**2 ) * vhdumi 
              if ( sqrt(rd)*rg2pkpc .lt. rhmaxd .and. 
     &             dv .lt. dvdup ) then
                if ( vhmax(ic2) .le. vhmax(ic1) ) then
                  iSh(ic2) = 1 
                else
                  iSh(ic1) = 1 
                  go to 5 
                endif
              endif
            endif
          enddo
        endif
 5      continue 
      enddo

      nnew = nil 
      do ic1 = 1 , nhalo
        if ( iSh(ic1) .eq. nil ) then 
          nnew = nnew + 1
          call HOld2New ( ic1 , nnew ) 
        endif        
      enddo

      nhalo = nnew 

      write(*,*) 'removing vel. duplicates: ',nhalo, ' haloes survived'

      return
      end

c     -------------------------
      subroutine Cannibalism ()
c     -------------------------
c     neighboring haloes eat each other 
c     ---------------------------------
      include 'hfind.h'

      do ic1 = 1 , nhalo 
        iSh(ic1) = nil 
      enddo 

      do ic1 = 1 , nhalo-1
        if ( iSh(ic1) .eq. nil ) then
          do ic2 = ic1+1 , nhalo 
            if ( iSh(ic2) .eq. nil ) then 
              diff_x = xh(ic2) - xh(ic1) 
              diff_y = yh(ic2) - yh(ic1) 
              diff_z = zh(ic2) - zh(ic1)
              corr_x = zero 
              corr_y = zero 
              corr_z = zero 
              if ( abs(diff_x) .gt. ngrid/2 ) then 
                if ( diff_x .gt. zero ) then 
                  corr_x = -ngrid 
                else
                  corr_x = ngrid
                endif
              endif
              if ( abs(diff_y) .gt. ngrid/2 ) then 
                if ( diff_y .gt. zero ) then 
                  corr_y = -ngrid 
                else
                  corr_y = ngrid
                endif
              endif
              if ( abs(diff_z) .gt. ngrid/2 ) then 
                if ( diff_z .gt. zero ) then 
                  corr_z = -ngrid 
                else
                  corr_z = ngrid
                endif
              endif
              rd = (xh(ic2) - xh(ic1) + corr_x)**2 +
     &             (yh(ic2) - yh(ic1) + corr_y)**2 +
     &             (zh(ic2) - zh(ic1) + corr_z)**2 
              rmax = max(rh(ic1),rh(ic2))
              if ( sqrt(rd) .lt. rmax ) then
                if ( nhp(ic2) .le. nhp(ic1) ) then
                  iSh(ic2) = 1 
                else
                  iSh(ic1) = 1 
                  go to 5 
                endif
              endif
            endif
          enddo
        endif
 5      continue 
      enddo

      nnew = nil 
      do ic1 = 1 , nhalo
        if ( iSh(ic1) .eq. nil ) then 
          nnew = nnew + 1
          call HOld2New ( ic1 , nnew ) 
        endif        
      enddo

      nhalo = nnew 

      write(*,*) 'cannibalism: ',nhalo, ' haloes survived'

      return
      end

c     ----------------------------
      subroutine Cannibalism200 ()
c     ----------------------------
c     neighboring haloes eat each other after growing 
c     -----------------------------------------------
      include 'hfind.h'

      do ic1 = 1 , nhalo 
        iSh(ic1) = nil 
      enddo 

      do ic1 = 1 , nhalo
        if ( iSh(ic1) .eq. nil ) then
          do ic2 = 1 , nhalo 
            if ( (iSh(ic2) .eq. nil)
     &                 .and. 
     &           ( ic1 .ne. ic2   ) ) then 
              diff_x = xh(ic2) - xh(ic1) 
              diff_y = yh(ic2) - yh(ic1) 
              diff_z = zh(ic2) - zh(ic1)
              corr_x = zero 
              corr_y = zero 
              corr_z = zero 
              if ( abs(diff_x) .gt. ngrid/2 ) then 
                if ( diff_x .gt. zero ) then 
                  corr_x = -ngrid 
                else
                  corr_x = ngrid
                endif
              endif
              if ( abs(diff_y) .gt. ngrid/2 ) then 
                if ( diff_y .gt. zero ) then 
                  corr_y = -ngrid 
                else
                  corr_y = ngrid
                endif
              endif
              if ( abs(diff_z) .gt. ngrid/2 ) then 
                if ( diff_z .gt. zero ) then 
                  corr_z = -ngrid 
                else
                  corr_z = ngrid
                endif
              endif
              rd = (xh(ic2) - xh(ic1) + corr_x)**2 +
     &             (yh(ic2) - yh(ic1) + corr_y)**2 +
     &             (zh(ic2) - zh(ic1) + corr_z)**2 
              if ( rd .lt. rhvir(ic1)**2 ) then
                if ( amh(ic2) .le. amh(ic1) ) then
                  iSh(ic2) = 1 
                else
                  iSh(ic1) = 1 
                  go to 5 
                endif
              endif
            endif
          enddo
        endif
 5      continue 
      enddo

      nnew = nil 
      do ic1 = 1 , nhalo
        if ( iSh(ic1) .eq. nil ) then 
          nnew = nnew + 1
          call HOld2New ( ic1 , nnew ) 
        endif        
      enddo

      nhalo = nnew 

      write(*,*) 'cannibalism: ',nhalo, ' haloes survived'

      close (40) 

      return
      end

c     -----------------------------------
      subroutine HOld2New ( nold , nnew )
c     -----------------------------------
c     copy halo variables from nold to nnew location
c     -----------------------------------------------
      include 'hfind.h'
c
          hc(nnew)  = hc(nold)
          xh(nnew)  = xh(nold) 
          yh(nnew)  = yh(nold) 
          zh(nnew)  = zh(nold) 
          vxh(nnew) = vxh(nold) 
          vyh(nnew) = vyh(nold) 
          vzh(nnew) = vzh(nold) 
          rh(nnew)  = rh(nold)
          rsh(nnew)  = rsh(nold)
          rhvir(nnew) = rhvir(nold)
          amh(nnew) = amh(nold)
          rst(nnew)  = rst(nold)
          nhp(nnew) = nhp(nold) 
          vhmax(nnew) = vhmax(nold) 
          rhmax(nnew) = rhmax(nold) 
          iHP(nnew) = iHP(nold) 
          iLH(nnew) = iLH(nold) 
       
      do i = 0 , nbmax
        pn(i,nnew) = pn(i,nold)
        pnt(i,nnew) = pnt(i,nold)
        na(i,nnew) = na(i,nold)
        pvx(i,nnew) = pvx(i,nold)
        pvy(i,nnew) = pvy(i,nold)
        pvz(i,nnew) = pvz(i,nold)
        vcirc(i,nnew) = vcirc(i,nold)
        vrms(i,nnew) = vrms(i,nold)
      enddo
c
      return
      end

c     -------------------------------
      function ve ( r , rmax , vmax )
c     -------------------------------
c
c     computes escape velocity for a NFW profile
c     c is approximated roughly
c
c.... compute escape velocity
c
      x = r/rmax
      ve = 2.15 * vmax 
     &          * sqrt(
     &                  log(1.0+2.0*x) / x
     &                )       
      return
      end

c     --------------------
      function con ( vmass )
c     --------------------
c     
c     concentration for mass vmass
c
      dimension cm(9)
c
c     approximation for c(M):
c
c     M < 10**   11.0  11.5  12.0  12.5  13.0  13.5  14.0  14.5  15.0
      data cm /  20.0, 17.8, 16.0, 14.3, 12.3, 11.0,  9.6,  8.3,  7.2 /
c
c.... approximate c
c
      if ( vmass .le. 0.5e11 ) then
        con  = 20.0
      else
        im = int(log10(vmass)/0.5) - 20
        if ( im .gt. 9 ) then
          con = 7.0
        else
          con = cm(im)
        endif
      endif

      return
      end

c     --------------------------------------
      subroutine Write_Short_Halo_Catalog ()
c     --------------------------------------
      include 'hfind.h'

      open ( 33 , file = 'haloshort.dat' )
      write( 33 , * ) Header
      do ic1 = 1 , nhalo 
        rkpc   = rh(ic1) * Cg2Mpc * 1000.0 * hubble
        amassh = pmmsun * nhp(ic1) 

        write(33,94)  xh(ic1) ,  yh(ic1) ,  zh(ic1) ,
     &                vxkms   ,  vykms   , vzkms , 
     &                rh(ic1)    ,  amassh  , nhp(ic1)
      enddo

 94   format ( 3(1x,f6.2), 2x, 3(1x,f7.2) , 3x , f10.4,
     &           1x, g13.5,1x,i6 )

      close ( 33 ) 

      return
      end

c     ------------------------------------------------------------
      subroutine Read_Particles_Binary ( fname1 , fname2, fname3 )
c     ------------------------------------------------------------
c
c     purpose: opens files with control information and particle data
c              reads in particle coordinates and momenta
c
      include 'hfind.h'
      integer ret
      character*256 fname1, fname2, fname3
      real maxdnb, mindnb

      hydro_flag = 1

      nfn1 = index ( fname1 , ' ' ) - 1 
      nfn2 = index ( fname2 , ' ' ) - 1
      nfn3 = index ( fname3 , ' ' ) - 1

      write(*,'(A)') fname1(1:nfn1)
      write(*,'(A)') fname2(1:nfn2)
      write(*,'(A)') fname3(1:nfn3)

      open ( 23 ,file =fname1(1:nfn1), form = 'unformatted')

c.... read control information and check whether it has proper structure
      read      (23,IOSTAT=ret) HEADER, 
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,
     &                  hubble,Wp5,Ocurv,Omb0,extras

      write(*,*) 'Omb0 = ', Omb0
      write(*,*) 'aexpn = ', AEXPN

      if ( ret .ne. 0 ) then
        write(*,*) 'Not hydro, trying N-body format...',ret
        rewind(23)
        hydro_flag = 0
        read      (23,IOSTAT=ret) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     &                  NROWC,NGRIDC,nspecies,Nseed,Om0,Oml0,
     &                  hubble,Wp5,Ocurv,extras
      endif

      if ( ret .ne. 0 ) then
       write(*,*) 'Error reading from header file!'
       stop
      endif

      if(nrow .ne. NROWC) then
         write(0,*) 'Code NROW value: ', nrow
         write(0,*) 'File NROW value: ', NROWC
         write(0,*) 'I don''t like that, chao.'
         stop
      endif

      if((( zero1 .eq. 0.0 ) .or. ( zero1 .eq. 0.1234 )) .and.
     &     (( zero2 .eq. 0.0 ) .or. ( zero2 .eq. 0.1234 )) .and.
     &     ( DelDC .ne. 0.0 )) then
         AEXPN = abox
         write(6,*) 'Detected DC mode, setting AEXPN to: ', abox
      endif

      write (*,100) HEADER,
     &                  AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     &                  EKIN,EKIN1,EKIN2,
     &                  NROWC,NGRID,NRECL,Om0,Oml0,Omb0,hubble
 100  format (1X,'Header=>',A45,/
     &           1X,' A=',F8.3,' A0=',F8.3,' Ampl=',F8.3,' Step=',F8.3,/
     &           1X,' I =',I4,' WEIGHT=',F8.3,' Ekin=',3E12.3,/
     &           1X,' Nrow=',I4,' Ngrid=',I4,' Nrecl=',I8,/
     &           1x,' Omega_0=',F7.3,' OmLam_0=',F7.4,
     &              ' Omegab_0=',F7.3,' Hubble=',f7.3)
      if ( NGRIDC .ne. NGRID ) then 
        write(*,*) 'NGRIDC .ne. NGRID :',NGRIDC, NGRID
        write(*,*) 'hope this is ok...'
      endif
      if(nspecies .eq. 0 ) nspecies = 1
      If( npmax .lt. lspecies(1) ) then
      write (*,*) ' Wrong number of particles: '
      write (*,*) ' must be at least=',lspecies(nspecies),' (lspecies)'
      write (*,*) ' but is set to ',npmax,' in a_setup.h...'
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec)
         enddo
         STOP
      Endif 

      if ( hydro_flag .eq. 1 ) then
c       rescale masses for purposes of computing properties
        do ispec = 1, nspecies
           wspecies(ispec) = wspecies(ispec) / (1.0-Omb0/Om0)
           write(*,*) 'Rescaled particle mass wpecies(', ispec, ') = ',
     &            wspecies(ispec)
        enddo
      endif

      
      nbyte  = nrecl * floatsize
      nacces = nbyte / nbyteword
 
      open ( 22 , file = fname2(1:nfn2), access = 'direct',
     &          status = 'unknown', recl = nacces )

      rewind 23

c     Only the first specie is used for halo finding
      N_particles = lspecies(1)
      np = lspecies(1)

      Npages      = (N_particles -1)/npage + 1
      N_in_last   = N_particles - npage*(Npages-1)
      write (*,*) ' Pages=',Npages,' Species=',Nspecies
      write (*,*) 'N_particles =',N_particles, ' N_in_last=',N_in_last
         do ispec = 1 , nspecies
           write(*,*) ispec, lspecies(ispec)
         enddo

      do ip = 1 , lspecies(1)
        pw(ip) = wspecies(1)
      enddo 

c     Read in densities (computed by smooth)
      open (21,file = fname3(1:nfn3),
     &         form='unformatted' )
      read(21) npd
      if ( npd .ne. np ) then 
        write(*,*) '# of particles in density and particle files'
        write(*,*) ' do not match:'
        write(*,*) 'npd =',npd
        write(*,*) 'np =',np
        write(*,*) 'aexpnd =',aexpnd
        stop
      endif
      
      read(21) (dnb(i),i=1,np)
      close (21)

c     rescale dm densities
      maxdnb = -1.e10
      mindnb = 1.e10
      if ( hydro_flag .eq. 1 ) then
        write(*,*) 'Scaling densities by 1-fb'
        do i=1,np
          dnb(i) = dnb(i) / (1.0 - Omb0/Om0)
          maxdnb = max( dnb(i), maxdnb )
          mindnb = min( dnb(i), mindnb )
        enddo
      endif

      write(*,*) 'Density range: ', mindnb, maxdnb

      ip = 0 
      
      do irow = 1 , Npages         ! loop over particle pages
        In_page = npage
        if ( irow .eq. Npages ) In_page = N_in_last
        iL = npage * (irow-1)
        CALL GetRow(irow,22) ! read in a page of particles
        do in = 1 , In_page          ! Loop over particles
          ipp = in + iL                     ! current particle number
          if ( ipp .le. np ) then 
            ip = ip + 1
          xp(ip) = xpar(in)
          yp(ip) = ypar(in)
          zp(ip) = zpar(in)
          vxp(ip) = vxx(in)
          vyp(ip) = vyy(in)
          vzp(ip) = vzz(in)
          pw(ip) = pw(ipp)
          dnb(ip) = dnb(ipp)
          endif
        enddo
      enddo
      np = ip 
      write(*,*) 'will use np=',np,' particles...'

      close (22)
      close (23)

      return
      end

C*********************************************************************
C                       Read current PAGE from disk
C
C                       NRECL - length of ROW block in words
 
      SUBROUTINE GETROW(IROW,Ifile)
      INCLUDE 'hfind.h'
         READ  (Ifile,REC=IROW) RECDAT
      RETURN
      END


c     -----------------------------
      real function rinput ( text )
c     -----------------------------
      character text*(*)
       
      write (*,'(A,$)') text
      read (*,*) x
      rinput = x

      return
      end

c     -------------------------------------
      character*40 function tinput ( text )
c     -------------------------------------
      character text*(*)
      character*40 textline
       
      write (*,'(A,$)') text
      read (*, 90 ) textline
      tinput = textline

 90   format (A)
      return
      end

c     -----------------------------
      SUBROUTINE indexx(n,arr,indx)
c     -----------------------------
c     (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
c     ----------------------------------------------------------
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

      SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
      INTEGER mwt,ndata
      REAL a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
CU    USES gammq
      INTEGER i
      REAL sigdat,ss,st2,sx,sxoss,sy,t,wt,gammq
      sx=0.
      sy=0.
      st2=0.
      b=0.
      if(mwt.ne.0) then
        ss=0.
        do 11 i=1,ndata
          wt=1./(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do 13 i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((1.+sx*sx/(ss*st2))/ss)
      sigb=sqrt(1./st2)
      chi2=0.
      if(mwt.eq.0) then
        do 15 i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
15      continue
        q=1.
        sigdat=sqrt(chi2/(ndata-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do 16 i=1,ndata
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
16      continue
        q=gammq(0.5*(ndata-2),0.5*chi2)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
      FUNCTION gammq(a,x)
      REAL a,gammq,x
CU    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.
      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .-35?421.1-9.

