c     program ps_test
c
c     Basic driver to test projection subrs ps_[xy,latlon,scale]
c
c     parameter (pi = 3.14159265358979)
c
c     xoffa = 0.
c     yoffa = 0.
c
c     dd0 = 10.e3
c
c-----------------------
c     do ilat=-85,-65,10
c     do ilon=-180,180,45
c------------------------
c     zlatd = ilat
c     zlond = ilon
c     zlat = zlatd*pi/180.
c     zlon = zlond*pi/180.
c
c     call ps_xy (zx, zy, zlat, zlon, xoffa, yoffa)
c
c     call ps_latlon (zlatchk, zlonchk, zx, zy, xoffa, yoffa)
c
c     call ps_scale (zscalex, zscaley, zx, zy, xoffa, yoffa, dd0)
c
c     write (6,'(a,2f10.2,2f16.10)') 'lat,lon, zscale[x,y]=', 
c    *                        zlatd,zlond, zscalex,zscaley
c----------
c     enddo
c     enddo
c----------
c
c     end

c-----------------------------------------------------------------------

      subroutine ps_xy (zx, zy, zlat, zlon, xoffa, yoffa)

c     Returns zx,zy (m) polar-stereo coordinates given
c     zlat,zlon (rad). xoffa,yoffa = offset from pole.

c     Code segment in comiceprojection.h sets ecc,arad,stdparallel,
c     alons,stdparallel,ts,ams, for hard-coded ellipsoid and projection.

#include <comiceprojection.h>

      zlata = abs(zlat)
      t  = tan (0.25*pi - 0.5*zlata)
     *   / ( ((1.-ecc*sin(zlata)) / (1.+ecc*sin(zlata))) ** (0.5*ecc) )

      zdlon = sign(1.,stdparallel) * (zlon-alons)
      zx = arad * ams * (t/ts) * cos (zdlon) - xoffa
      zy = arad * ams * (t/ts) * sin (zdlon) - yoffa

      return
      end

c-----------------------------------------------------------------------

      subroutine ps_latlon (zlat, zlon, zx, zy, xoffa, yoffa)

c     Returns zlat,zlon (rad) given polar-stereo coordinates 
c     zx,zy (m). xoffa,yoffa = offset from pole.

c     Code segment in comiceprojection.h sets ecc,arad,stdparallel,
c     alons,stdparallel,ts,ams, for hard-coded ellipsoid and projection.

#include <comiceprojection.h>

      tp = sqrt ( (zx+xoffa)**2 + (zy+yoffa)**2 ) * ts / (arad*ams)

      xsi = sign(1.,stdparallel) * (0.5*pi - 2.*atan(tp))

      zlat = xsi + (         ecc**2 / 2.
     *               +   5.* ecc**4 / 24.
     *               +       ecc**6 / 12.
     *               +  13.* ecc**8 / 360.   ) * sin(2.*xsi)
  
     *           + (     7.* ecc**4 / 48.
     *               +  29.* ecc**6 / 240.
     *               + 811.* ecc**8 / 11520. ) * sin(4.*xsi)

     *           + (     7.* ecc**6 / 120.
     *               +  81.* ecc**8 / 1120.  ) * sin(6.*xsi)
 
     *           + (  4279.* ecc**8 / 161280.) * sin(8.*xsi)
      
      zlon = alons + sign(1.,stdparallel) * atan2 (zy+yoffa, zx+xoffa) 
      if (zlon.gt.pi)  zlon = zlon - 2.*pi
      if (zlon.lt.-pi) zlon = zlon + 2.*pi

      return
      end

c-----------------------------------------------------------------------

      subroutine ps_scale (zscalex,zscaley, zx,zy, xoffa,yoffa, zdelin)

c     Returns scale factors zscalex, zscaley (eareth distance / p.s. 
c     plane distance) for small displacements in x,y directions,
c     given polar-stereo coordinates zx,zy (m). xoffa,yoffa = offset 
c     from pole. zdel = length of displacement (m).  

c     Code segment in comiceprojection.h sets ecc,arad,stdparallel,
c     alons,stdparallel,ts,ams, for given ellipsoid and projection.

#include <comiceprojection.h>

c---
c     my old: exact only for ecc=0 (sphere)
c     zlambda = 1. + cos(-abs(stdparallel)*pi/180. + 0.5*pi)
c     zr = sqrt ((zx+xoffa)**2 + (zy + yoffa)**2) / arad
c     zscalex = 2. / ( zlambda * (1. + (zr/zlambda)**2) )
c     zscaley = zscalex
c---
c     call ps_latlon (zlat, zlon, zx, zy, xoffa, yoffa)
c     zlata = abs(zlat)
c     zscalex = (1. + sin(zlata)) / (1. + sin(abs(stdparallel)*pi/180.))
c     zscaley = zscalex
c---

c     zdel = zdelin
      zdel = 1.e3

c        Calculate scale factor crudely by finding lats,lons of
c        slightly displaced points (for both x and y directions),
c        then calculating great-circle distance between them.

      call ps_latlon (zlat, zlon, zx, zy, xoffa, yoffa)

      call ps_latlon (zlat1, zlon1, zx-.5*zdel, zy, xoffa, yoffa)

      call ps_latlon (zlat2, zlon2, zx+.5*zdel, zy, xoffa, yoffa)

      call ps_latlon (zlat3, zlon3, zx, zy-.5*zdel, xoffa, yoffa)

      call ps_latlon (zlat4, zlon4, zx, zy+.5*zdel, xoffa, yoffa)

c        Convert geodetic to geocentric latitudes

      zlatgeod = zlat
      zlat  = atan (tan(zlat) *(1.-ecc**2))
      zlat1 = atan (tan(zlat1)*(1.-ecc**2))
      zlat2 = atan (tan(zlat2)*(1.-ecc**2))
      zlat3 = atan (tan(zlat3)*(1.-ecc**2))
      zlat4 = atan (tan(zlat4)*(1.-ecc**2))

c        Compute Cartesian coordinates and arcs between vectors

      xcart1 = cos(zlat1)*cos(zlon1)
      ycart1 = cos(zlat1)*sin(zlon1)
      zcart1 = sin(zlat1)
      xcart2 = cos(zlat2)*cos(zlon2)
      ycart2 = cos(zlat2)*sin(zlon2)
      zcart2 = sin(zlat2)
      xcart3 = cos(zlat3)*cos(zlon3)
      ycart3 = cos(zlat3)*sin(zlon3)
      zcart3 = sin(zlat3)
      xcart4 = cos(zlat4)*cos(zlon4)
      ycart4 = cos(zlat4)*sin(zlon4)
      zcart4 = sin(zlat4)
      zarcx= acos (min (1., xcart1*xcart2+ycart1*ycart2+zcart1*zcart2) )
      zarcy= acos (min (1., xcart3*xcart4+ycart3*ycart4+zcart3*zcart4) )

c     zarcx= sqrt ((xcart1-xcart2)**2 + (ycart1-ycart2)**2 
c    *           + (zcart1-zcart2)**2)
c     zarcy= sqrt ((xcart3-xcart4)**2 + (ycart3-ycart4)**2 
c    *           + (zcart3-zcart4)**2)

c       zr is radius from geocentric center 

      brad = arad * sqrt (1.-ecc**2)
      zr = arad*brad / sqrt ((brad*cos(zlat))**2 + (arad*sin(zlat))**2)

c        Allow for angle (beta) between tangent planes of ellipsoid vs. 
c        sphere. Small effect, O(ecc**2), but needed to get good results
c        (nearly isotropic, not quite).

c~~~~~~~~
c Either: 
c~~~~~~~~
c     see "angle factor" notes 2/28/14:
c     ztanbeta = -(arad**2 - brad**2) * cos(zlat) * sin(zlat)
c    *         / ((brad*cos(zlat))**2 + (arad*sin(zlat))**2)
c     zbeta = abs (atan(ztanbeta))
c~~~~
c Or (gives virtually the same as above):
c~~~~
      zbeta = abs (zlatgeod-zlat)

c        Factor only applies to "radial" direction from pole, and is 
c        zero for "circumference" direction. 
c        For general x-direction: cos ( sin-1 (cos(theta) sin(beta) ) )
c        For general y-direction: cos ( sin-1 (sin(theta) sin(beta) ) )
c        where theta = atan(y/x). Results are nearly isotropic
c        for all theta (differences between zscalex and zscaley
c        are ~1.e-10), and exactly isotropic for theta = +/- 45 deg.

      zsinbeta = sin (zbeta)
      ztheta = atan2 (zy+yoffa, zx+xoffa)
      zangx = asin (abs(cos(ztheta)*zsinbeta))
      zangy = asin (abs(sin(ztheta)*zsinbeta))

      zdearthx = zarcx * zr / cos(zangx)
      zdearthy = zarcy * zr / cos(zangy)

      zscalex = zdearthx / zdel
      zscaley = zdearthy / zdel

c     zscale = sqrt (zscalex*zscaley)

      return
      end
