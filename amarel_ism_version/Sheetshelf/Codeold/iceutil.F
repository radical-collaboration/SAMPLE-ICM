      subroutine itor (iarr, arr, n)
      dimension iarr(n), arr(n)
      do i=1,n
        arr(i) = iarr(i)
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine rtoi (arr, iarr, n)
      dimension arr(n), iarr(n)
      do i=1,n
        iarr(i) = nint(arr(i))
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      FUNCTION LENCHR(CH)
c     Returns length of character variable ch, less trailing blanks
      CHARACTER*(*) CH
      DO 10 I = LEN(CH), 1, -1
        LENCHR = I
        IF (CH(I:I).NE.' ' .AND. CH(I:I).NE.CHAR(0)) RETURN
   10 CONTINUE
      RETURN
      END
c
      SUBROUTINE ZERO (ARR, NAR)
c     Zeros nar words starting at arr(1)
      DIMENSION ARR(NAR)
      DO 10 J=1,NAR
        ARR(J) = 0.
   10 CONTINUE
      RETURN
      END
c
      SUBROUTINE IZERO (IARR, NAR)
c     Zeros nar words starting at iarr(1)
      DIMENSION IARR(NAR)
      DO 10 J=1,NAR
        IARR(J) = 0
   10 CONTINUE
      RETURN
      END
c
c-----------------------------------------------------------------------
c
      subroutine reseti (ia, n, ival)
 
c        Sets integer array ia(n) to ival
 
      dimension ia(*)
 
      do 10 i=1,n
        ia(i) = ival
   10 continue
 
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine resetr (a, n, val)

c        Sets real array a(n) to val

      dimension a(*)

      do 10 i=1,n
        a(i) = val
   10 continue

      return
      end

c-----------------------------------------------------------------------

      subroutine scopy_i (n, iarr, inca, ibrr, incb)
c     Like Cray scopy, except for integers
      dimension iarr(*), ibrr(*)
      ib = 1
      do ia=1,n,inca
        ibrr(ib) = iarr(ia)
        ib = ib + incb
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
c#if defined (SUN)
cc     to enable nint() with large values (eg, showice) with -i4:
cc     also declare in comdriveice.h:
c      integer*8 function nint(x)
c      nint = x  + sign(0.5, x)
c      return
c      end
c#endif
c
c-----------------------------------------------------------------------
c
#if defined (SGI) || defined (SUN)

      integer function flush(iu)
c     Don't do this for std input or output:
      if (iu.eq.0 .or. iu.eq.5 .or. iu.eq.6) then
c     if (iu.eq.5) then
         flush = 0
         return
      endif
      close (iu, err = 1000, iostat = ios1)
#if defined (SUN) && defined (F90)
 1000 open  (iu, status='old', position='append', err=2000, iostat=ios2)
#else
 1000 open  (iu, status='old', access  ='append', err=2000, iostat=ios2)
#endif
 2000 if (ios1.eq.0 .and. ios2.eq.0) then
        flush = 0
      else
        flush = -1
        write (6,*) 
     *    '*** Warning (flush): iu=',iu,' ios1=',ios1,' ios2=',ios2
      endif
      return
      end

#elif defined (AIX)

      integer function flush(lun)
      ! - call XLF general service and utility procedure:
      integer*4 lun4
      lun4 = lun
      call flush_(lun4)
      flush = 0
      return
      end

#endif
c
c-----------------------------------------------------------------------
c
#if ! defined (CRAY)
      function cvmgt (x,y,l)
c     Duplicates Cray-Vector-Merge-GT function.
c     Only safe when first two arguments are reals (especially
c     if integers are 32 bits and reals are 64 bits).
      logical l
      if (l) then
        cvmgt = x
      else
        cvmgt = y
      endif
      return
      end

      subroutine scopy (n, a, inca, b, incb)
c     Duplicates Cray scopy 
      dimension a(*), b(*)
      ib = 1
      do ia=1,n,inca
        b(ib) = a(ia)
        ib = ib + incb
      enddo
      return
      end

      integer function ishell (a)
c     Duplicates Cray ishell
      character*(*) a
#if defined (AIX)
      integer*4 result
      ! - result is the return status of a waitpid(2) call and hence
      ! - if status is non-zero it will be different from the return
      ! - status if the same command were exectued from command line
      ! - as waitpid(2) encodes additional information -- see man pg.
      call system (a,result)
      ishell = result
#elif defined (SUN)
      integer*4 system
      ishell = system (a)
#else
      integer system
      ishell = system (a)
#endif
      return
      end
#endif
c
c-----------------------------------------------------------------------
c
#if defined (SUN)
      integer function hand1()
      write (6,*) '*** STOP (exception: fpe) in hand1'
c     write (6,*) '*** STOP (exception: division) in hand1'
      hand1 = 0.
      stop
      end

c     integer function hand2()
c     write (6,*) '*** STOP (exception: overflow) in hand2'
c     hand2 = 0.
c     stop
c     end

c     integer function hand3()
c     write (6,*) '*** STOP (exception: invalid) in hand3'
c     hand3 = 0.
c     stop
c     end
#endif

c-----------------------------------------------------------------------

      subroutine insolhalf (iyear, alorb, fluxh, ecc,obliq, prec)

c        Computes summer half-year insolation fluxh (W/m2)
c        at latitude alorb (radians), for calendar year iyear.
c        iyear will be interpreted by zenorb* as "BP", ie, relative 
c        to 1950 AD. ecc,obliq,prec are passed.

      parameter (nflux=365)  ! changed (from 180), 11/3/2011
c     parameter (nflux=180)
c     parameter (nflux=24)   ! old fur.f
      dimension flux(nflux)

      parameter (pi=3.14159265358979)

      parameter (solcon = 1367.)          ! W/m2

      dimension
     *  cosq(1), fraq(1), cosq24(1), fraq24(1), alatq(1)

      alatq(1) = alorb

c        Calculate vern as in old zenorb

      zs = sin (0.5*prec) / sqrt ((1.+ecc)/(1.-ecc))
      zc = cos (0.5*prec) 
      ze = 2. * atan2 (zs,zc)
      vern = ze - ecc * sin(ze)

c        Loop through year in equal time increments, storing
c        24-hour mean fluxes at latitude alorb in flux().
c        Then bubble sort flux and take the mean of the largest
c        half of its members, to get the half-year insolation
c        (fluxh) at alorb.

      isecdy = 86400/2
      dt  = 86400.

      do m=1,nflux
        isecyr = nint (86400.*365.*(m-.5)/nflux)
        call zencal (iyear, isecyr, isecdy, dt,
     *               ecc, obliq, prec, vern, dist, eccf,
     *               cosq, fraq, cosq24, fraq24, alatq, 1)
        flux(m) = cosq(1)*fraq(1)*eccf*solcon
      enddo

      call bubblesort (flux, nflux)

      fluxh = 0.
c     nwid = nint (182.5/(365./nflux))   ! half-year
      nwid = nint (90./(365./nflux))     ! 90-day (changed 11/3/2011)
      do m = nflux-nwid+1, nflux
        fluxh = fluxh + flux(m)
      enddo
      fluxh = fluxh / nwid

      return
      end

c----------------------------------------------------------------------

      subroutine bubblesort (a, n)

c        Bubble sorts a(n) into ascending order

      dimension a(n)

      do ibub = n-1,1,-1
        do j=ibub,n-1
          if (a(j).gt.a(j+1)) then
            z = a(j+1)
            a(j+1) = a(j)
            a(j) = z
          else
            go to 10
          endif
        enddo
   10   continue
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine bubblemean (iarr, a, w, zmean, ilon, ilat, nlon,nlat,n)

c        (1) Bubble sorts a(n) into ascending order,sorting weights
c            w(n) to match.
c        (2) Subtracts weighted mean from a, and normalize w.
c        (3) Packs a and normalized w values into iarr(ilon,ilat):
c            iarr = 10000*nint(a+10000) + min(9999,nint(10000*w))
c            (So a+10000 must be >= 0 and w must be between 0 and 1).

      dimension iarr(nlon,nlat,n), a(n), w(n)

c        Bubble sort a(n) into ascending order, with w(n) following.

      do ibub = n-1,1,-1
        do j=ibub,n-1
          if (a(j).gt.a(j+1)) then
            z = a(j+1)
            a(j+1) = a(j)
            a(j) = z

            z = w(j+1)
            w(j+1) = w(j)
            w(j) = z
          else
            go to 10
          endif
        enddo
   10   continue
      enddo

c        Subtract out the mean of a (weighted by w), and normalize w

      zmean = 0.
      zarea = 0.
      do i=1,n
        zmean = zmean + a(i)*w(i)
        zarea = zarea +      w(i)
      enddo
      zmean = zmean/zarea
      do i=1,n
        a(i) = a(i) - zmean
        w(i) = w(i) / zarea
      enddo

c        Pack values and weights into iarr

      do i=1,n
        iarr(ilon,ilat,i) =
     *    10000*nint(a(i)+10000.) + min (9999, nint(10000.*w(i)))
      enddo

      return
      end

c----------------------------------------------------------------------

      SUBROUTINE ZENCAL (IYEAR, ISECYR, ISECDY, DT,
     *                   ECC, OBLIQ, PREC,VERN,DIST, ECCF,
     *                   COSQ, FRAQ, COSQ24, FRAQ24, ALAT, NLAT)
c
c        Calculates COSQ(J) and FRAQ(J), cos (zenith angle), 
c        and daylight fraction in interval DT, versus latitude.
c        Also calculates 24-hr means in COSQ24 and FRAQ24.
c        Also calculates DIST, current earth-sun distance relative to
c        the semi-major axis, and eccentricity factor ECCF = 1/DIST**2.
c
c     IYEAR = year number (not used) (supplied)
c     ISECYR= secs into current calendar year from 00:00 Jan 1st (supp)
c     ISECDY= secs into current calendar day  from 00:00 (supplied)
c     DT    = time interval (must be .le. 1 day) (seconds) (supplied)
c     ECC   = eccentricity (supplied)
c     OBLIQ = obliquity (radians) (supplied)
c     PREC  = precession (prograde angle between perihelion and n.h.
c             vernal equinox) (0 to 2*pi) (supplied)
c     VERN  = 2*pi *  time from perihelion to n.h. vernal equinox 
c             / perihelion-to-perihelion year (supplied)
c     DIST  = current earth-sun distance / semi-major axis (returned)
c     ECCF  = eccentricity factor (1/DIST**2) (returned)
c     COSQ  = cos (zenith angle) averaged over DT vs lat (returned)
c     FRAQ  = fraction of DT with daylight vs lat (returned)
c     COSQ24= 24-hour mean of cos(zenith angle) vs lat (returned)
c     FRAQ24= 24-hr fraction with daylight vs lat (returned)
c     ALAT  = latitudes (radians, cannot be pi/2 or -pi/2) (supplied)
c     NLAT  = latitudinal grid size (supplied)
c
c     Local variables:
c     SECPD = day length (constant, approx midnight-to-midnight) (secs)
c     AYEAR = perihelion-to-perihelion year (anomalistic) (days)
c     TYEAR = equinox-to-equinox year (tropical) (days)
c 
      DIMENSION COSQ(NLAT), FRAQ(NLAT), COSQ24(NLAT), FRAQ24(NLAT),
     *          ALAT(NLAT)
      DIMENSION NDAYPM(12)
      PARAMETER (PI = 3.14159265358979)
c
      DATA NDAYPM /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
c
      DATA SECPD, AYEAR, TYEAR / .864E5, 365.25964, 365.24219 /
c
      DATA EPSORB / 1.E-7 /    ! or 1.e-5 for real*4
c
      TPI = 2.*PI
      PIH = PI/2.
      RPD = PI/180.
c
c        Calculate time from present to n.h.vernal equinox, by computing
c        time T from Jan 1st to present, then subtracting time TVE 
c        from Jan 1st to vernal equinox. TVE is incremented each year
c        from 1989 by the difference between the tropical year TYEAR 
c        (equinox-to-equinox) and 365 days, and decremented by 1 day
c        following each leap year. We define the calendar by setting 
c        N.H. vernal equinox at 12:00 GMT March 21, 79.5 days after 
c        00:00 GMT Jan 1 (as recommended for 6K BP by PMIP)
c
      VERNCAL = 79.5
c
      T = ISECYR
      TVE = VERNCAL
c     However, model ignores leap years, so comment out the
c     adjustment to TVE from 1989.
c     TVE = TVE + (IYEAR-1989)*(TYEAR-365.) - (IYEAR-1989)/4
      T = T - SECPD*TVE
c
c        Change to "time" from perihelion, expressed as
c        2*pi * fraction of perihelion-to-perihelion year
c
      T = VERN + TPI*T/(AYEAR*SECPD)
      IF (T.LT.0.) T = T + TPI
c
c        Solve for THETA, angle between perihelion and current position,
c        using Newton's method for ecc anomaly E (Danby, Fundamentals
c        of Celestial Mechanics, 1962, Eqs. (6.3.12) and (6.3.19).)
c
      E = T
      DO 20 ILOOP=1,100
        DE = - (E-ECC*SIN(E)-T) / (1.-ECC*COS(E))
        IF (ABS(DE).LE.EPSORB) GOTO 22
        E = E + DE
   20 CONTINUE
   22 CONTINUE
c
      ZS = SIN (0.5*E) * SQRT((1.+ECC)/(1.-ECC))
      ZC = COS (0.5*E)
      THETA = 2. * ATAN2 (ZS,ZC)
c
c        Calculate DIST (earth-sun distance / semi-major axis) and ECCF
c
      DIST = (1.-ECC*ECC) /  (1. + ECC*COS(THETA))
      ECCF = 1./(DIST*DIST)
c
c        Calculate B (angle between Earth spin axis and Earth-Sun line)
c        and CB = COS(B), SB = SIN(B)
c
      CB = SIN(THETA-PREC) * COS(PIH-OBLIQ)
c     Need ACOS to return 0 to pi for SB to have the right sign
      B = ACOS(CB)
      SB = SIN(B)
c
c        Loop over latitudes (nb: ALAT cannot be pi/2 or -pi/2)
c
      DO 100 J=1,NLAT

        CL = COS (PIH-ALAT(J))
        SL = SIN (PIH-ALAT(J))
c
c          Set RISE and SET, hour angles (between 0 and 2*pi, from 
c          local midnight) of sunrise and sunset
c
        IF (CL*CB .GT. SL*SB)  THEN
c         Polar day
          X = PI
        ELSE IF (CL*CB .LT. -SL*SB)  THEN
c         Polar night
          X = 0.
        ELSE
c         sunrise and sunset occur
c         (need ACOS to return 0 to pi) 
          X = ACOS(-CL*CB/(SL*SB))
        ENDIF
        RISE = PI - X
        SET  = PI + X
c
c          Set 24-hr total cos(zen) TCOSQ24 and daylight time TIMLQ24
c
        TCOSQ24 = 2.* (CL*CB*X + SL*SB*SIN(X))
        TIMLQ24 = SET - RISE
c
c          For an arbitrary longitude, set hour angles (0 to 2*pi,
c          from local  midnight) of T1 and T2, the start and end of
c          interval DT.
c
        ALON = 0. ! (greenwich, arbitrarily)
        T1 = ALON + ISECDY*TPI/SECPD
        T2 = T1 + DT*TPI/SECPD
        T1 = MOD (T1+TPI, TPI)
        T2 = MOD (T2+TPI, TPI)
c
c          Flag if local midnight occurs during DT,
c          and if so switch T1 and T2. Then will do calcs below for
c          24 hrs - DT, and convert back at end. First test for daily
c          mean case (DT = 24 hrs) and handle similarly.
c
        IFLAG = 1
        X = ABS ( MOD (T2-T1, TPI) )
        IF (X.LT.EPSORB .OR. ABS(X-TPI).LT.EPSORB) THEN
c         daily mean (DT = 24hrs)
          IFLAG = -1
          T2 = T1
        ELSE IF (T1. GT. T2) THEN
c         midnight occurs during DT (and not daily mean)
          IFLAG = -1
          X = T1
          T1 = T2
          T2 = X
        ENDIF
c
c          Constrain integ limits A1, A2 to be between RISE and SET,
c          and integrate to get total cos(zen) TCOSQ and daylight time
c          TIMLQ. Subtract PI from A1,A2 to get in range -PI to PI.
c
        A1 = MAX (RISE, MIN (SET, T1)) - PI
        A2 = MAX (RISE, MIN (SET, T2)) - PI
        TCOSQ = CL*CB*(A2-A1) + SL*SB*(SIN(A2)-SIN(A1))
        TIMLQ = A2 - A1
c
c          If local midnight occurs within DT, subtract from 24-hr tots
c
        IF (IFLAG.EQ.-1) THEN
          TCOSQ = TCOSQ24 - TCOSQ
          TIMLQ = TIMLQ24 - TIMLQ
        ENDIF
c
c          Set avg daylight cos(zen) COSQ and daylight fraction FRAQ
c          for time interval DT
c
        IF (TIMLQ.GT.0.) THEN
          COSQ(J) = TCOSQ / TIMLQ
        ELSE
          COSQ(J) = 0.
        ENDIF
        FRAQ(J) = TIMLQ / (TPI*DT/SECPD)
c
c          Set 24-hr mean cos(zen) and 24-hr daylight fraction
c
        COSQ24(J) = TCOSQ24 / TPI
        FRAQ24(J) = TIMLQ24 / TPI
c
  100 CONTINUE
c
      RETURN
      END
   
c-----------------------------------------------------------------------
   
      subroutine checkrcode (rcode, cfilin, cvar, ilab)
      character*(*) cfilin, cvar
      integer*4 rcode
#include <netcdf.inc>
      character cfile*240, msg*80
   
      if (rcode.ne.0) then
        cfile = cfilin
        msg = nf_strerror(rcode)
        write (6,'(/a/2x,a/2a, a,i6, a,i6/2a/)')
     *    'Error (netcdf file):',
     *    cfile(1:lenchr(cfile)),
     *    '  var: ',cvar, '   label:',ilab, '   rcode:',rcode,
     *    '  explanation: ',msg(1:lenchr(msg))
        stop    
      endif   
   
      return  
      end  
   
