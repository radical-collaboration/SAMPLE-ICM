c comiceprojection.h

c     Used in subrs ps_xy, ps_latlon, ps_scale (iceprojection.F),
c     to set properties of ellipsoid and polar stereographic projection,
c     selected by label cproj.

c     For settings below and equations in iceprojection.F, see: 
c     "OGP Publication 373-7-2 – Geomatics Guidance Note number 7, 
c     part 2 – June 2013. To facilitate improvement, this document is 
c     subject to revision. The current version is available at 
c     www.epsg.org"

      parameter (pi = 3.14159265358979)

c======

c----------------------------------------
#if defined (EISANTA) && defined (STEREO)
c----------------------------------------

c     WGS84 ellipsoid, southern hemisphere:
      ecc = .081819191
      arad = 6378137.0
      stdparallel = -71.0
      alons = 90. * pi/180.      ! longitude line defined by y=0, x->pos

c     sphere:
c     ecc = 0.
c     arad = 6371220.
c     stdparallel = -71.0
c     alons = 90. * pi/180.      ! longitude line defined by y=0, x->pos

c-----------------------------------------
#elif defined (NHA) && defined (GREENLAND)
c-----------------------------------------

c     WGS84 ellipsoid, projection for Bamber et al. Greenland dataset,
c     used only to interpolate that data to model lon,lat grid
      ecc = .081819191
      arad = 6378137.0
      stdparallel = 71.0
c     x=0, y->neg in data is at -39E, so 51E is y=0, x->pos:
      alons = 51. * pi/180.      ! longitude line defined by y=0, x->pos

c------------------------
#elif ! defined (EISLINE)
c------------------------

      write (6,'(3(/a))') '*** Error:',
     *  '    Unexpected call to ps_[xy,latlon,scale]', 
     *  '    without -DEISANTA and -DSTEREO, or -DNHA and -DGREENLAND'  
      stop

c------
# endif
c------

# if ! defined (EISLINE) 
      alats = abs (stdparallel) * pi/180.

      ts = tan (0.25*pi - 0.5*alats)
     *   / ( ((1.-ecc*sin(alats)) / (1.+ecc*sin(alats))) ** (0.5*ecc) )

      ams = cos(alats) / sqrt (1. - (ecc*sin(alats))**2)
#endif

