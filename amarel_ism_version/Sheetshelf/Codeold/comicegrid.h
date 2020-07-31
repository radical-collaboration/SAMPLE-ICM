c comicegrid.h

      parameter (pi     = 3.14159265358979)
      parameter (radius = 6.371220e6) ! old, std 
c     parameter (radius = 6.378137e6) ! WGS84 semi-major axis
c     parameter (radius = 6.356752e6) ! WGS84 semi-minor axis

c===================
#if defined (LONLAT)
c===================

c----------------------
#  if defined (EISLINE)
c----------------------

#    if defined (SLMODEL)
c     parameter (dlatd = 0.05,        dlond = dlatd)           ! ~5   km
      parameter (dlatd = 0.1,         dlond = dlatd)           ! ~10  km
c     parameter (dlatd = 0.01,        dlond = dlatd)           ! ~1   km
c     parameter (dlatd = 0.001,       dlond = dlatd)           ! ~0.1 km
      parameter (alon1 = 0.,          alon2 = alon1 + dlond)   !    nx=1
      parameter (alat1 =-90.,         alat2 =-60.)
c     test of numerics, with ny=1, around equator(!):
c     parameter (alon1 = 0.,          alon2 = 30.)
c     parameter (alat1 = 0.,          alat2 = alat1 + dlatd)
#    endif

c--------------------
#  elif defined (NHA)
c--------------------

c++++++++++++++++++++++++++
#    if defined (GREENLAND)
c++++++++++++++++++++++++++

c............................
#      if ! defined (NESTING)
c............................

#        if defined (RESOL80)
      parameter (dlatd =   0.80,      dlond =  1.60)
#        elif defined (RESOL40)
      parameter (dlatd =   0.40,      dlond =  0.80)
#        elif defined (RESOL20)
      parameter (dlatd =   0.20,      dlond =  0.40)
#        elif defined (RESOL10)
      parameter (dlatd =   0.10,      dlond =  0.20)
#        elif defined (RESOL5)
      parameter (dlatd =   0.05,      dlond =  0.10)
#        else
      parameter (dlatd =   0.20,      dlond =  0.40)
#        endif

      parameter (alon1 =  -76.,       alon2 = -10.)
      parameter (alat1 =   58.,       alat2 =  85.)

c............................
#      elif defined (NESTING)
c............................

#        if defined (RESOL1)
      parameter (dlatd =   0.01,      dlond =  0.02)
#        elif defined (RESOL2)
      parameter (dlatd =   0.02,      dlond =  0.04)
#        else 
      parameter (dlatd =   0.02,      dlond =  0.04)
#        endif

c     Jakobshavn: 69.17N, 49.83W
c old:
c     parameter (alon1 = -49.83 - 2.50,  alon2 = -49.83 + 2.50)  ! Jakob
c     parameter (alat1 =  69.17 - 1.25,  alat2 =  69.17 + 1.25)  ! Jakob
c cur:
      parameter (alon1 = -49.83 - 2.00,  alon2 = -49.83 + 2.00)  ! Jakob
      parameter (alat1 =  69.17 - 0.75,  alat2 =  69.17 + 0.75)  ! Jakob
c too restrictive, but used in debplotice3dnet_mel (nlonin, etc): 
c     parameter (alon1 = -50.10 - 1.30,  alon2 = -50.10 + 1.30)  ! Jakob
c     parameter (alat1 =  69.17 - 0.25,  alat2 =  69.17 + 0.25)  ! Jakob

c     Helheim: 66.35N, 50.22W (bombs):
c     parameter (alon1 = -38.20 - 2.50,  alon2 = -38.20 + 2.50)  ! Helhe
c     parameter (alat1 =  66.35 - 1.25,  alat2 =  66.35 + 1.25)  ! Helhe
cc    parameter (alon1 = -38.20 - 5.00,  alon2 = -38.20 + 5.00)  ! Helhe
cc    parameter (alat1 =  66.35 - 2.50,  alat2 =  66.35 + 2.50)  ! Helhe

c     Kangerlussuaq (Russell Gl.): 67.10N, 50.22W (bombs):
c     parameter (alon1 = -50.22 - 5.00,  alon2 = -50.22 + 5.00)  ! Russ
c     parameter (alat1 =  67.10 - 2.50,  alat2 =  67.10 + 2.50)  ! Russ

c...........
#      endif
c...........

c+++++++++++++++++++++++++++
#    elif defined (LOVECLIP)
c+++++++++++++++++++++++++++

      parameter (dlatd =   0.5,       dlond =  1.0)       ! loveclip std
      parameter (alon1 = -180.,       alon2 = 180.)       ! loveclip std
      parameter (alat1 =   20.,       alat2 =  85.)       ! loveclip std

c++++++++++++++++++++++++++
#    elif defined (NHAGLAC)
c++++++++++++++++++++++++++

#      if defined (NHAGLOB)
      parameter (dlatd =   1.0,       dlond =  2.0)
      parameter (alon1 = -180.,       alon2 =  180.)
      parameter (alat1 =  -90.,       alat2 =  90.)
#      else
      parameter (dlatd =   1.0,       dlond =  2.0)
c     parameter (alon1 = -180.,       alon2 = -10.)
c     parameter (alat1 =   30.,       alat2 =  85.)
      parameter (alon1 = -180.,       alon2 =  180.)
      parameter (alat1 =   20.,       alat2 =  90.)
#      endif

c++++++++++++++++++++++++++
#    elif defined (GLACBAY)
c++++++++++++++++++++++++++

#      if defined (RESOL05)
      parameter (dlatd =    .005,     dlond =     .01)
      parameter (alon1 = -137.85,     alon2 = -135.50)  ! 235
      parameter (alat1 =   58.25,     alat2 =   59.50)  ! 250
#      elif defined (RESOL1)
      parameter (dlatd =     .01,     dlond =     .02)
      parameter (alon1 = -137.85,     alon2 = -135.49)  ! 118
      parameter (alat1 =   58.25,     alat2 =   59.50)  ! 125
#      endif

c++++++++
#    else
c++++++++

      parameter (dlatd =   2.0,       dlond =  2.0)
c     parameter (dlatd =   1.0,       dlond =  2.0)
c     parameter (dlatd =   0.5,       dlond =  1.0)
      parameter (alon1 = -180.,       alon2 = -10.)
      parameter (alat1 =   10.,       alat2 =  86.)

c+++++++++
#    endif
c+++++++++

c---------------------------------------
#  elif defined (CARB) || defined (ORDO)
c---------------------------------------

      parameter (dlatd = 1.0,         dlond = 2.0)
      parameter (alon1 = -180.,       alon2 = 180.)
      parameter (alat1 = -90.,        alat2 = 90.) 

c------
#  else
c------

      parameter (dlatd = 0.5,         dlond = 0.5)
      parameter (alon1 = -0.5*dlond,  alon2 = 0.5*dlond) 
      parameter (alat1 = 30.,         alat2 = 70.) 

c-------
#  endif
c-------

      parameter (dlat = dlatd*pi/180., 
     *           dlon = dlond*pi/180.,
     *           nx   = (alon2-alon1)/dlond + .001,  
     *           ny   = (alat2-alat1)/dlatd + .001,
     *           dd0  = dlat*radius, dx0=dd0, dy0=dd0)
      parameter (nxp=nx+1, nyp=ny+1)

      common /lonlat/
     *  alon(nx),     clon(nx),
     *  alat(ny),     clat(ny),
     *  alatu(0:nyp), clatu(0:nyp),
     *  alatv(0:nyp), clatv(0:nyp)

c=====================
#elif defined (STEREO)
c=====================

c----------------------
#  if defined (EISLINE)
c----------------------

#    if defined (LINEB)
      parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
      parameter (nx =80,    ny=1)
#    elif defined (LINEC)
c     parameter (dx0=100.e3,dy0=dx0,  dd0=dx0)
      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0=10.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0= 5.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0= 2.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0= 1.e3, dy0=dx0,  dd0=dx0)

#      if defined (FLATANT)
      parameter (nx =2100.e3/dx0 + .001,   ny=1)
#      else
      parameter (nx =1600.e3/dx0 + .001,   ny=1)
#      endif

#    elif defined (LINED)
      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)
      parameter (nx =40,    ny=1)
#    elif defined (LINEE)
c     parameter (dx0=0.1e3, dy0=dx0,  dd0=dx0)
      parameter (dx0=0.5e3, dy0=dx0,  dd0=dx0)
      parameter (nx =50.e3/dx0 + .001,  ny=1)
#    elif defined (LINEF)
c     parameter (dx0=0.1e3, dy0=dx0, dd0=dx0)
c     parameter (dx0=0.2e3, dy0=dx0, dd0=dx0)
c     parameter (dx0=0.5e3, dy0=dx0, dd0=dx0)
c     parameter (dx0=1.e3,  dy0=dx0, dd0=dx0)
c     parameter (dx0=2.e3,  dy0=dx0, dd0=dx0)
      parameter (dx0=5.e3,  dy0=dx0, dd0=dx0)
c     parameter (dx0=10.e3, dy0=dx0, dd0=dx0)
c     parameter (dx0=40.e3, dy0=dx0, dd0=dx0)
c     parameter (nx =1000.e3/dx0 + .001,  ny=1)
c     parameter (nx =1000.e3/dx0 + .001,  ny=5)
c     parameter (ny =1000.e3/dx0 + .001,  nx=1)
c     parameter (ny =1000.e3/dx0 + .001,  nx=5)
      parameter (nx =1500.e3/dx0 + .001,  ny=1)   ! schoof
#    elif defined (LINEG)
      parameter (dx0=10.e3, dy0=dx0,  dd0=dx0)
c     parameter (nx =10.e3/dx0   + .001,  ny=1)     ! 1pt
      parameter (nx =1000.e3/dx0 + .001,  ny=1)
#    elif defined (LINEH)
      parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
      parameter (nx =4000.e3/dx0 + .001,  ny=1)
#    elif defined (LINEH2)
      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0= 1.e3, dy0=dx0,  dd0=dx0)
      parameter (nx =2000.e3/dx0 + .001,  ny=1)
#    elif defined (LINEH3)
      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
      parameter (nx =2000.e3/dx0 + .001,  ny=1)
#    elif defined (LINEH4)
      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)
c     parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
      parameter (nx =4000.e3/dx0 + .001,  ny=nx)

#    elif defined (LINEM)

#      if defined (RESOL1)
      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)
#      elif defined (RESOL2)
      parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
#      elif defined (RESOL3)
      parameter (dx0=10.e3, dy0=dx0,  dd0=dx0)
#      elif defined (RESOL4)
      parameter (dx0= 5.e3, dy0=dx0,  dd0=dx0)
#      elif defined (RESOL5)
      parameter (dx0= 1.e3, dy0=dx0,  dd0=dx0)
#      elif defined (RESOL6)
      parameter (dx0= 0.1e3,dy0=dx0,  dd0=dx0)
#      endif
      parameter (nx = 1000.e3/dx0 + .001,  ny=1)

#    elif defined (SLMODEL)
      parameter (dx0=10.e3, dy0=dx0,  dd0=dx0)
      parameter (nx = 3300.e3/dx0 + .001,  ny=1)

#    endif

c------------------------
#  elif defined (EISANTA)
c------------------------

c...................................................
#    if ! defined (NESTING) && ! defined (TRANSECTA)
c...................................................

#      if defined (RESOL160)
      parameter (dx0=160.e3, dy0=dx0,  dd0=dx0)
      parameter (nx=35,     ny=35)

#      elif defined (RESOL80)
      parameter (dx0=80.e3, dy0=dx0,  dd0=dx0)
      parameter (nx=70,     ny=70)

#      elif defined (RESOL40)
      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)
      parameter (nx=140,    ny=140)                        ! was 141,141

#      elif defined (RESOL20)
      parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
      parameter (nx=280,    ny=280)                        ! was 281,281

#      elif defined (RESOL16)
      parameter (dx0=16.e3,  dy0=dx0,  dd0=dx0)
      parameter (nx=381,   ny=381)

#      elif defined (RESOL12)
      parameter (dx0=12.e3, dy0=dx0,  dd0=dx0)
      parameter (nx=466,    ny=466) ! nominal

#      elif defined (RESOL10)
      parameter (dx0=10.e3, dy0=dx0,  dd0=dx0)
      parameter (nx=560,    ny=560) ! too big with melange (500 ok)

#      elif defined (RESOL5)
      parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
      parameter (nx=1120,   ny=1120)

#      else 
      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)             ! default=40
      parameter (nx=140,    ny=140)
c     parameter (dx0=20.e3, dy0=dx0,  dd0=dx0)
c     parameter (nx=280,    ny=280)
c     parameter (dx0=10.e3, dy0=dx0,  dd0=dx0)
c     parameter (nx=560,    ny=560)

#      endif

c...................................................
#    elif ! defined (NESTING) && defined (TRANSECTA)
c...................................................
 
c~~~~~~~~~~~~~~~~~~~~~~~~~
#      if defined (ZOOMGL)
c~~~~~~~~~~~~~~~~~~~~~~~~~

      parameter (dx0=5.e3, dy0=dx0, dd0=dx0)
      parameter (nx0=303,  ny=1)
      parameter (dzoom = 0.1e3)
      parameter (npzoom = dx0/dzoom + .001) 
      parameter (iwidzoom = 1)
      parameter (nx = nx0 + (2*iwidzoom+1)*(npzoom-1))

      parameter (nzoom = nx0*npzoom)  
      common /czoomgl/ 
     *  xhzoom(nzoom),    yhzoom(nzoom),
     *  alondzoom(nzoom), alatdzoom(nzoom),
     *  distzoom(nzoom),  distlen(nx), totlen

c~~~~~~~~~~~
#       else
c~~~~~~~~~~~

      parameter (dx0=5.e3, dy0=dx0,  dd0=dx0)
      parameter (nx=303,   ny=1)
      common /ctransecta/ distlen(nx), totlen 

c~~~~~~~~~~~
#      endif
c~~~~~~~~~~~

c..........................
#    elif defined (NESTING)
c..........................

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      if defined (RESOLASE10)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=151,    ny=151)
      parameter (dx0=10.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -75.1, alondnest0 = -104.0) ! ase  10km
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (RESOLASE5)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=201,    ny=201)
      parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -75.1, alondnest0 = -104.0) ! ase  5km
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (RESOLASE2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     parameter (nx=221,    ny=221)
c     parameter (dx0=2.e3,  dy0=dx0,  dd0=dx0)
c     parameter (alatdnest0 = -75.1, alondnest0 = -104.0) ! ase  2km
c     parameter (alatdnest0 = -75.1, alondnest0 = -100.0) ! pig  2km
c     parameter (alatdnest0 = -75.6, alondnest0 = -106.7) ! thw  2km
c     parameter (alatdnest0 = -74.5, alondnest0 = -111.0) ! bear 2km

      parameter (nx=221,    ny=221)
      parameter (dx0=2.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -75.6, alondnest0 = -106.7) ! thw  2km
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (RESOLASE1)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     parameter (nx=221,    ny=221)
c     parameter (nx=331,    ny=331)  ! A7
      parameter (nx=401,    ny=221)  ! A8a,b, A8a
c     parameter (nx=801,    ny=221)  ! A8c
c     parameter (nx=801,    ny=221)  ! A9a

      parameter (dx0=1.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -74.5, alondnest0 = -111.0) ! bear 1km

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (RESOLPEN5)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     use with -DCLEMLIM:
      parameter (nx=300,    ny=300)
      parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -70.0, alondnest0 = -65.0)  !peninsula 5km

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (RESOLPEN2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     use with -DCLEMLIM:
      parameter (nx=350,    ny=350)
      parameter (dx0=2.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -70.0, alondnest0 = -65.0)  !peninsula 2km

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (PIRRIT5)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=200,    ny=200)
      parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -81.13, alondnest0 = -85.39) !Pirrit Hills

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (WILKES20)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=100,    ny=75)
      parameter (dx0=20.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -73.0, alondnest0 =  151.) ! Wilkes 20km 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (WILKES10)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=200,    ny=150)
      parameter (dx0=10.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -73.0, alondnest0 =  151.) ! Wilkes 10km 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (WILKES5)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=400,    ny=300)
      parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -73.0, alondnest0 =  151.) ! Wilkes 5km 
c     parameter (nx=200,    ny=200)                      ! new_dec12
c     parameter (dx0=10.e3,  dy0=dx0,  dd0=dx0)          ! new_dec12
c     parameter (alatdnest0 = -73.0, alondnest0 =  150.) ! Wilkes 10km

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (AURORA5)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=330,    ny=330)
      parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -70.5, alondnest0 =  114.) ! Aurora 5km 

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (RECOVERY5)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=300,    ny=30.)
      parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -84.0, alondnest0 =  -16.) ! Recovery 5km 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (WAIS20)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=171,    ny=171)
      parameter (dx0=20.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -82.,  alondnest0 = -85.)  ! all WAIS 20km

c~~~~~~~~~~~~~~~~~~~~~~~~~~~
#      elif defined (WAIS10)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~
      parameter (nx=341,    ny=341)
      parameter (dx0=10.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -82.,  alondnest0 = -85.)  ! all WAIS 10km

c~~~~~~~~~~
#      else
c~~~~~~~~~~
      parameter (nx=341,    ny=341)
      parameter (dx0=10.e3,  dy0=dx0,  dd0=dx0)
      parameter (alatdnest0 = -82.,  alondnest0 = -85.)  ! all WAIS 10km
c pite:
c     parameter (nx=171,    ny=171)
c     parameter (dx0=20.e3,  dy0=dx0,  dd0=dx0)
c     parameter (alatdnest0 = -82.,  alondnest0 = -85.)  ! all WAIS 20km
c+++
c     parameter (nx=201,    ny=201)
c     parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
c     parameter (alatdnest0 = -82.,  alondnest0 =-155.)  ! Siple 5km
c+++
c     parameter (nx=201,    ny=201)
c     parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)
c     parameter (alatdnest0 = -78.,  alondnest0 = 166.)  ! Ross Is. 5km

c     parameter (nx=201,    ny=201)
c     parameter (dx0=10.e3,  dy0=dx0,  dd0=dx0)
c     parameter (alatdnest0 = -80., alondnest0 = -160.)  ! Ross Is. 10km
c+++
c     parameter (nx=221,    ny=221)                      ! pitk 5km
c     parameter (dx0=5.e3,  dy0=dx0,  dd0=dx0)           ! pitk 5km
c     parameter (alatdnest0 = -75.0, alondnest0 = -105.) ! pitk 5km

c     parameter (nx=111,    ny=111)                      ! pig  2km
c     parameter (dx0=2.e3,  dy0=dx0,  dd0=dx0)           ! pig  2km
c     parameter (alatdnest0 = -75.1, alondnest0 = -100.0)! pig  2km

cc    parameter (nx=111,    ny=111)                      ! pitk 10 km
c     parameter (nx=166,    ny=111)                      ! pitk 10 km
c     parameter (dx0=10.e3,  dy0=dx0,  dd0=dx0)          ! pitk 10 km
cc    parameter (alatdnest0 = -75.0, alondnest0 = -105.) ! pitk 10 km
c     parameter (alatdnest0 = -72.5, alondnest0 = -105.) ! pitk 10 km
c~~~~~~~~~~~
#      endif
c~~~~~~~~~~~

c.........
#    endif
c.........

c-----------------------
#  elif defined (TEST2D)
c-----------------------

c     for mismipplus:
      parameter (xdomain = 640.e3, ydomain = 80.e3)
#    if defined (RESOL1) 
      parameter (dx0=1.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL2) 
      parameter (dx0=2.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL5) 
      parameter (dx0=5.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL10) 
      parameter (dx0=10.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL20) 
      parameter (dx0=20.e3, dy0=dx0, dd0=dx0)
#    else
      parameter (dx0=1.e3, dy0=dx0, dd0=dx0)   ! 1 km default
#    endif
      parameter (nx=xdomain/dx0 + .001, ny=ydomain/dy0 + .001)

c--------------------------
#  elif defined (TESTMEL2D)
c--------------------------

#    if defined (RESOL1) 
      parameter (dx0=1.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL2) 
      parameter (dx0=2.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL2P5) 
      parameter (dx0=2.5e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL5) 
      parameter (dx0=5.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL10) 
      parameter (dx0=10.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL20) 
        parameter (dx0=20.e3, dy0=dx0, dd0=dx0)
#    elif defined (RESOL40) 
      parameter (dx0=40.e3, dy0=dx0, dd0=dx0)
#    else
      parameter (dx0=20.e3, dy0=dx0, dd0=dx0)
#    endif

c............................
#    if defined (TESTMELPROF)
c............................
c     parameter (xdomain = 100.e3, ydomain = dy0)       ! 1-D
#      if defined (YDOMAINHALF) 
      parameter (xdomain = 300.e3, ydomain = 50.e3)     ! 2-D
#      else
      parameter (xdomain = 300.e3, ydomain = 100.e3)    ! 2-D
#      endif
c.......
#    else
c.......
c     parameter (xdomain = 700.e3, ydomain = dy0)       ! 1-D
      parameter (xdomain = 700.e3, ydomain = 200.e3)    ! 2-D
c.........
#    endif
c.........

      parameter (nx=xdomain/dx0 + .001, ny=ydomain/dy0 + .001)

c------
#  else
c------

      parameter (dx0=40.e3, dy0=dx0,  dd0=dx0)
      parameter (nx=40,     ny=1)
c     parameter (nx=10,     ny=10)
c     parameter (nx=40,     ny=40)
c     parameter (nx=150,    ny=150)

c-------
#  endif
c-------
      parameter (nxp=nx+1,  nyp=ny+1)

      common /cxyoffc/  xoffa, yoffa

c=====
#endif
c=====

c((((((((((((((((((((((((((((((((((((((((
#if defined (NESTING) && defined (LONLAT)
c((((((((((((((((((((((((((((((((((((((((
c     greenland, orig resol20:
c     parameter (nxa=nint((-10.+76.)/.40), nya=nint((85.-58.)/.20))
c     greenland, orig resol10:
      parameter (nxa=nint((-10.+76.)/.20), nya=nint((85.-58.)/.10))
c     greenland, orig resol5:
c     parameter (nxa=nint((-10.+76.)/.10), nya=nint((85.-58.)/.05))

c((((((((((((((((((((((((((((((((((((((((((
#elif defined (NESTING) && defined (STEREO)
c((((((((((((((((((((((((((((((((((((((((((
#  if defined (PIRRIT5)
      parameter (dxa=40.e3,  dya=dxa,  dda=dxa)
      parameter (nxa=140,    nya=140)
#  elif defined  (WILKES20) || defined (WILKES10) || defined (WILKES5) || defined (AURORA5) || defined (RECOVERY5) || defined (WAIS20)  || defined (WAIS10)
      parameter (dxa=10.e3,  dya=dxa,  dda=dxa)
      parameter (nxa=560,    nya=560)
#  else
      parameter (dxa=20.e3,  dya=dxa,  dda=dxa)
      parameter (nxa=280,    nya=280)
c     parameter (dxa=10.e3,  dya=dxa,  dda=dxa)     ! for cliff
c     parameter (nxa=560,    nya=560)               ! for cliff
#  endif
c(((((
#endif
c(((((

      common /hgrid/
     *  darea(nx,ny),     totarea,
     *  xh (nx,ny),       yh (nx,ny),
     *  dx (nx,ny),       dy (nx,ny), 
     *  alond(nx,ny),     alatd(nx,ny),
     *  dxu(0:nxp,0:nyp), dyu(0:nxp,0:nyp),
     *  dxv(0:nxp,0:nyp), dyv(0:nxp,0:nyp),
     *  dxc(0:nxp,0:nyp), dyc(0:nxp,0:nyp)

      parameter (nlev=10, nlevp=nlev+1)
c     parameter (nlev=20, nlevp=nlev+1)
c     parameter (nlev=30, nlevp=nlev+1)

#if defined (SEDIMENT) 
      parameter (nsed = 10)
#else
      parameter (nsed = 1)
#endif
      parameter (nsedp = nsed+1)
       
      parameter (nbed = 1)
c     parameter (nbed = 6)
      parameter (nbedp = nbed+1)

      common /vgrid/
     *  dzeta(nlev),    zeta(0:nlevp),
     *  dzetah(0:nlev), zetah(0:nlev),

     *  zsed(0:nsedp), zsedm(0:nsed), dzsed(nsed),
     *  zbed(0:nbedp), zbedm(0:nbed), dzbed(nbed),

     *  bedthick

#if defined (CALVDAMAGE)
c     tracer(*,*,*,2) = "r" (Bassis and Ma,EPSL,2015):
      parameter (ntraca=1, ntrace=2)
#else
      parameter (ntraca=1, ntrace=1)
#endif

      parameter (interm=5, ioterm=6, iunamel=10, iuokend=11,
     *           iuout1dzoom=18, iuout2d=19,iuout1d=20, iuplot1d=21, 
     *           iutab=22, iusedbud=23, iusedtrack=24,
     *           iutabwais=25, iupts=26,
     *           iutabbud=28, iutabmel=280, iutabmelforce=281,  
     *           iutabtro=29, iutabtrohb=290,
     *           iures=30, iunest=31, iucrhmel=32,
     *           iueis=40, iueistab=41,
     *           iubed=45, iubedout=46, 
     *           iunh=47, iubas=48, iutsin=49, iuaga=50,
     *           iuto=60, iupr=61, iuta=62, iuru=63,
     *           iuclim=64,
     *           iuslg=71, iusli=72, iuslb=73, iuslout=74, iuslit=75,
     *           iuslout2=76, iusle=77, 
     *           iuhis=92, iulag=93,
     *           iumismipplus_a=98, iumismipplus_b=99)

c     for gcm (ascii) or rcm (netcdf) meteorol forcing (iceclimgcm.F)

#if defined (GCMMATRIX) || defined (CO2INTER) || defined (RCMMATRIX) || defined (RACMOCLIM) || defined (RACMOANTSEAS) || defined (LOVECLIP) || defined (CSMCLIM) || defined (CLEMCLIM) || defined (CSMORDO)

#  if defined (NHAOBS)
      parameter (nlong=96, nlatg=48)                           ! GCM T31

#  elif defined (NHAGLAC)
      parameter (nlong=96, nlatg=48)                           ! GCM T31

#  elif defined (RCMMATRIX)
c     parameter (nlong= 91, nlatg= 91, dd0rcm=80.e3)           ! Ant 80
      parameter (nlong=181, nlatg=181, dd0rcm=40.e3)           ! Ant 40

#  elif defined (RACMOCLIM) || defined (RACMOANTSEAS) || defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM)

#    if defined (GREENLAND)
      parameter (nlong=306, nlatg=312, 
     *           nlonrac=nlong, nlatrac=nlatg)
#    elif defined (EISANTA)
      parameter (nlong=262, nlatg=240, 
     *           nlonrac=nlong, nlatrac=nlatg)
#    endif
      common /gcmgridrac/ 
     *  alondrac(nlonrac,nlatrac), alatdrac(nlonrac,nlatrac)

#  elif defined (LOVECLIP)
      parameter (nlong=64, nlatg=32)                           ! T21

#  elif defined (CSMCLIM)

#    if defined (CSMHIRES)
      parameter (nlong=288, nlatg=192) 
#    else
      parameter (nlong=96,  nlatg=48)
#    endif

#    if defined (NTCSM)
      parameter (ntcsm=NTCSM)  ! J. Tsai
#    else
      parameter  ntcsm=5000)   ! J. Tsai all years
c     parameter (ntcsm=1)      ! J. Tsai
#    endif

#  elif defined (CLEMCLIM)

      parameter (nlong=192, nlatg=96) 

#  elif defined (CSMORDO)

      parameter (nlong=96, nlatg=48) 

#  else
      parameter (nlong=180, nlatg=90)                          ! GCM 2X2
#  endif

      common /gcmgrid/ along(nlong), alatg(nlatg)
#  if defined (LOVECLIP)
     *             ,wgauss(nlatg), darlc(nlong,nlatg)              ! T21
#  endif

#endif

c     longitudinal shift to Genesis gcm grid (deg eastward), passed 
c     to setinterp:
#if defined (EISANTA)
c     parameter (rotate_to_gcm = -15.) ! for 34 Ma  Genesis gcm 
      parameter (rotate_to_gcm = 0.)   ! for modern Genesis gcm
#endif

#if defined (LOVECLIP)
      parameter (nmat=1,  nmon=12)

#elif defined (NHA) 
c---
#  if defined (NHAOBS)
      parameter (nmat=1,  nmon=12)      ! modern observed climatol.
#  elif defined (NHAGLAC)
      parameter (nmat=7,  nmon=12)      ! gen2: [No,11k,Fu][co,wa], Obs
#  elif defined (CSMCLIM)
      parameter (nmat=1,  nmon=12)
#  else
      parameter (nmat=3,  nmon=12)      ! 3 orbs, modern 
#  endif
c---

#elif defined (ASYNCH)
      parameter (nmat=1,  nmon=12)

#elif defined (RCMMATRIX)

#  if defined (RCMMATCO2)
#    if defined (RAMPPLIO) || defined (RAMP2X) || defined (RAMP4X) || defined (RAMP8X) || defined (RAMP8XW)
      parameter (nmat=2,  nmon=12)      ! 1x, 2x/4x/8x/8xw 
#    else
      parameter (nmat=4,  nmon=12)      ! 1x, 2x, 4x, 8x
#    endif
#  else
      parameter (nmat=2,  nmon=12)      ! RCMFILE, RCMFILE or 1x, 1x
#  endif

#elif defined (RACMOCLIM) || defined (CSMCLIM) || defined (CLEMCLIM) || defined (CSMORDO)
      parameter (nmat=1,  nmon=12)

#else
      parameter (nmat=1,  nmon=12)
c     parameter (nmat=18, nmon=12)      ! 9x2 matrix, mis31 gcm
c     parameter (nmat=27, nmon=12)      ! 3x3x3 matrix
c     parameter (nmat=36, nmon=12)      ! 3x3x4 matrix
#endif

c     intra-annual timestep for surfmodel[1-4] (secs): 

#if defined (DTIMED)
      parameter (dtimed = DTIMED*3600.)             ! DTIMED is in hours
#else
      parameter (dtimed = 5.*86400.)
c     parameter (dtimed = 1.*86400.)
#endif
      parameter (nspy = nint((86400.*365)/dtimed)) ! # of timesteps/yr

c     minimum limits to handle small velocities in icedyn 

      parameter (gmin   = (0.01/dd0)**2)
      parameter (ubmin  =  0.30)            ! 0.01

c        For land-ocean outlines in printmap (maskcur), display for
c        balance calcs (maskbal)

      common /cmaskdisp/ maskcur(nx,ny), maskbal(nx,ny)

c        For 1-D flowline strip output in iceshow1d and iceplot1d

      common /c1ddisp/ i1s(nx+ny), j1s(nx+ny), d1s(nx+ny), n1s 

c        For sediment tracking in sedtracking_eul (and on restart file).
c        ntr*dtr must be >= the thickest sediment encountered in the run

#if defined (SEDTRACK)
      parameter (ntr = 2000, dtr = 1.)
#else
      parameter (ntr =    1, dtr = 1.)
#endif
      common /csedtrack/ sedtrack(nx,ny,0:ntr), timebot, itrtop(nx,ny)

      parameter ( npoimax = max(1,(nx*ny)/5), ! max # of pts in one lake
     *            nlakemax = 50        )      ! max # of lakes

c        For reading restart files, interp to different resolutions 
    
#if defined (EISLINE) 
      parameter (ninmx = 15002*3)
#elif defined (TEST2D) 
c     parameter (ninmx = 50000)
      parameter (ninmx = 250000)
#else
      parameter (ninmx = 702*702)
#endif
      common /cgridinterp/ nxin, nyin, dxin, dyin, 
     *                     xhin(ninmx), yhin(ninmx), xoffin, yoffin,
     *                     alondin(ninmx), alatdin(ninmx), 
     *                     versresin 

c       For nestdrive filename list (namelist, for nesting only)

      character*80 cnestpath, cnestlist(1000)
      common /comnest/ cnestpath, cnestlist, nnestlist

c       Number of ocean sectors for ocean melt (OCMARTIN)

#if defined (OCMARTIN) 
      parameter (nsecmax=4)
c#elif defined (OCOH)
c     parameter (nsecmax=7)
#else
      parameter (nsecmax=1)
#endif
      character*16 csector(nsecmax)
      common /comsec/ csector, nsec

c^^^^^^^^^^^^^^^^^^
#if defined (RWRAP)
c^^^^^^^^^^^^^^^^^^
c        Input variables passed from R-wrapper (see icectl.F),
c        copied from argument array to these common variables for 
c        use throughout code

      common /comrwrap/ 
     *  rwrap_ocfacmult,
     *  rwrap_ocfacmultase,
     *  rwrap_calvnick,
     *  rwrap_crhshelf,
     *  rwrap_tauasth,
     *  rwrap_calvliq,
     *  rwrap_cliffvmax,
     *  rwrap_facemeltrate,
     *  rwrap_enhancesheet,
     *  rwrap_enhanceshelf,
     *  rwrap_crhfac
c^^^^^
#endif
c^^^^^

