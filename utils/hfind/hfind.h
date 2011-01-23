c     ----------------------------------------------------------------
c
c     hfind.h - definitions for hfind.f (Andrey Kravtsov, december 1996)
c
c             - parameters suitable for CLS_ICs (Doug Rudd, April 2007)
c     ----------------------------------------------------------------

      parameter ( pi     = 3.1415926 )           !
      parameter ( G      = 6.6732e-8 )           ! grav const [cm^3 g^{-1} s^{-2}]
      parameter ( Solarm = 1.989e33  )           ! solar mass in grams
      parameter ( pc     = 3.0856e18 )           ! parsec in cm
      parameter ( akpc   = 3.0856e21 )           ! kiloparsec in cm


      parameter ( nil    = 0         )           ! integer zero 
      parameter ( zero   = 0.0       )           ! real zero
      parameter ( nrow   = 128       )            ! # of particles  in 1D
      parameter ( npmax  = 18000000  )            !
      parameter ( ncold  = nrow**3   )           !
      parameter ( npage  = nrow**2   )           !
      parameter ( nrecl  = npage * 6 )           !
      parameter ( ngrid  = 128       )           ! # of grid cells in 1D
      parameter ( nll    = 128       )           ! # of chain mesh cells in 1D
      parameter ( ncell  = ngrid**3  )           ! total # of grid cells
      parameter ( ncell0 = ngrid**3  )           ! # of 0-level cells in the original sim.
      parameter ( nh     = 400000    )            ! max # of haloes
      parameter ( floatsize = 8      )
      parameter ( nbyteword = 1      )  
      parameter ( rhoaver= 1.0       )           ! average density of particles
      parameter ( rinit  = 100.0     )           ! initial halo radius in kpc

      parameter ( xn = (1.0*ngrid)+1.-1.E-7 )    ! boundary parameters
      parameter ( yn = (1.0*ngrid))       

c.... grid arrays 
      common / GRID1 / iCL(nll,nll,nll)

c.... particle arrays

      common / PART0 / np 
      common / PART1  /  xp(npmax)               ! particle x-coordinates
      common / PART2  /  yp(npmax)               ! particle y-coordinates
      common / PART3  /  zp(npmax)               ! particle z-coordinates
      common / PART4  / vxp(npmax)               ! particle x-velocities
      common / PART5  / vyp(npmax)               ! particle y-velocities
      common / PART6  / vzp(npmax)               ! particle z-velocities
      common / PART12  / pw(npmax)               ! particle mass
      common / PART8  / dnb(npmax)               ! distance to the Nth clos.nb
      common / PART9  / ind(npmax)               ! index array for sorting     
      common / PART10 / iSp(npmax)               ! 0 - free, 1 - belongs to a halo
      common / PART11 / iLL(npmax)               ! particle-mesh linked list 
      common / PART12 /bind(npmax)               ! binding energy of particle
c
c.... halo arrays 
c
      integer hc(nh)
      common / HALO00 /  hc		         ! particle at center of halo
      common / HALO01 /  xh(nh)                  ! halo x-coordinates
      common / HALO02 /  yh(nh)                  ! halo y-coordinates
      common / HALO03 /  zh(nh)                  ! halo z-coordinates
      common / HALO04 / vxh(nh)                  ! halo x-velocities
      common / HALO05 / vyh(nh)                  ! halo y-velocities
      common / HALO06 / vzh(nh)                  ! halo z-velocities
      common / HALO07 /  rh(nh)                  ! halo radius up to DeltaMin
      common / HALO07 / rsh(nh)                  ! scale radius r_s
      common / HALO08 / rst(nh)                  ! stop radius (profile flattens)
      common / HALO09 / nhp(nh)                  ! # of parts inside rh
      common / HALO09 / amh(nh)                  ! halo mass 
      common / HALO10 / iHP(nh)                  ! halo-particle LL headers
      common / HALO11 / iLH(nh)                  ! particle-halo linked list
      common / HALO12 / iSh(nh)                  ! 1 - radius reached 0 - not
      common / HALO13 / aMc(nh)                  ! core mass
      common / HALO14 / DeltaMin , Deltavir      ! halo overdensity
      common / HALO15 / nhalo                    ! current # of haloes 
      common / HALO16 / vhmax(nh)                ! 
      common / HALO17 / rhmax(nh)                ! 
      common / HALO18 / rhvir(nh)                ! 

      parameter ( nbmax    = 150 )               ! size of the profile arrays
      real pn   (nil:nbmax,nh)                   ! mass  profile
      real pnt  (nil:nbmax,nh)                   ! tot mass within r 
      real na   (nil:nbmax,nh)                   ! auxiliary arrays:
      real pvx  (nil:nbmax,nh)                   ! x-velocity profile 
      real pvy  (nil:nbmax,nh)                   ! y-velocity profile 
      real pvz  (nil:nbmax,nh)                   ! z-velocity profile 
      real vcirc(nil:nbmax,nh)                   ! circular velocity profile
      real vrms (nil:nbmax,nh)                   ! 3D velocity dispersion profile
c
      common / PROFILES / pn, pnt, na, pvx, pvy, pvz, vcirc, vrms
c
c.... variables
c
      common / VARS / dr
c
c.... conversion factors
c
      common / FACTORS / box0 , box , pmmsun , 
     &                   rMpc2g , rkpc2g , rg2Mpc , rg2kpc , vg2kms,
     &                   rg2pMpc, rg2pkpc
c
c.... control variables, arrays for reading in binary format
c
      COMMON /CONTROL/ AEXPN,AEXP0,AMPLT,ASTEP,ISTEP,PARTW,
     +                 TINTG,EKIN,EKIN1,EKIN2,AU0,AEU0,
     +                 NROWC,NGRIDC,Nspecies,Nseed,Om0,Oml0,hubble,Wp5,
     +                 Ocurv,Omb0,extras(100)

      real               ns      
      COMMON / TRUNCOM / Om,Omb,Omc,Omnu,Par(6),ns,qqscaleb,
     &                   QSCALE, SCALEL

c     Change to real*4 for single-precision particle files
      real*8        XPAR(npage), YPAR(npage), ZPAR(npage),
     &                  VXX(npage), VYY(npage), VZZ(npage)
      real*8		RECDAT(NRECL)
	  common / ROW / xpar, ypar, zpar, vxx, vyy, vzz

      DIMENSION      wspecies(10),lspecies(10)
      EQUIVALENCE    (RECDAT(1),XPAR(1)),
     +               (wspecies(1),extras(1)),
     +               (lspecies(1),extras(11)),
     +               (zero1,extras(21)),
     +               (DelDC,extras(22)),
     +               (abox,extras(23)),
     +               (Hbox,extras(24)),
     +               (zero2,extras(25))

      CHARACTER*45      HEADER
      COMMON / HEADDR/  HEADER

      integer		hydro_flag
      common / HFLAG /  hydro_flag
