c@@@@@@@@@@@@@@@@@@@@@
#if defined (LOVECLIP)
c@@@@@@@@@@@@@@@@@@@@@

c-----------------------------------------------------------------------

      subroutine readloveclim (arr, cvarin, nmonin)

c     Reads nc file written by Loveclim before each ice model segment 
c     of coupled runs, for one field, returned in arr on LC grid.
c     LC grid is global spectral T21 (nlong=64, nlatg=32). 
c     Netcdf field names (cvarin) are: 
c       surfelev  = surface elevation (m)
c       tair2m    = 2-m surface air temperature (C)
c       precip    = precipitation (m/y liquid equiv)
c       evap      = evaporation (m/y liquid equiv)
c       solarflux = incident downward solar flux (W/m2), sfc or TOA(?)
c       tocean    = ocean temperature at ~400 m depth (C)
c       icefrac   = fraction of LC cell with ice
c       oceanfrac = fraction of LC cell with ocean
c     nmonin = number of months to be read (1 or 12)
c     Null data value is vnull_lc (set in comiceloveclim.h)

c     Can be called for anomaly method, using entry "readlovectl"(below)

#include <comicephys.h>
#include <comicegrid.h>
#include <comiceloveclim.h>
#include <netcdf.inc>

      dimension arr(nlong,nlatg,nmonin)
      character*(*) cvarin

      character cvar*24
      integer*4 ncid, rcode, varid, start(3), count(3)

      real*4 arrin(nlong,nlatg,12), var4(nlong+nlatg)

c     from .../Nhicemodel/fort.6_genesis_t21:
      dimension alatt21(32), wgausst21(32)
      data alatt21 /
     * -85.7605858535, -80.2687778865, -74.7445392644, -69.2129751469,
     * -63.6786346204, -58.1429531903, -52.6065252572, -47.0696413643,
     * -41.5324606331, -35.9950778795, -30.4575535112, -24.9199282618,
     * -19.3822310601, -13.8444835299,  -8.3067027338,  -2.7689029668,
     *   2.7689029668,   8.3067027338,  13.8444835299,  19.3822310601,
     *  24.9199282618,  30.4575535112,  35.9950778795,  41.5324606331,
     *  47.0696413643,  52.6065252572,  58.1429531903,  63.6786346204,
     *  69.2129751469,  74.7445392644,  80.2687778865,  85.7605858535 /
      data wgausst21 /
     *   0.0070186100,   0.0162743947,   0.0253920653,   0.0342738629,
     *   0.0428358980,   0.0509980593,   0.0586840935,   0.0658222228,
     *   0.0723457941,   0.0781938958,   0.0833119242,   0.0876520930,
     *   0.0911738787,   0.0938443991,   0.0956387201,   0.0965400885,
     *   0.0965400885,   0.0956387201,   0.0938443991,   0.0911738787,
     *   0.0876520930,   0.0833119242,   0.0781938958,   0.0723457941,
     *   0.0658222228,   0.0586840935,   0.0509980593,   0.0428358980,
     *   0.0342738629,   0.0253920653,   0.0162743947,   0.0070186100 /
c     nb: sum of wgausst21(1:32) = 2

      logical firstlc
      save firstlc
      data firstlc /.true./

      character cfile*240, clab*24

c        Set input Loveclim file. If "readlovectl" called (for control
c        climate in anomaly method, -DLOVEANOM, entry statement below,
c        use climfilectl. If not, use climfilein. 
c        climfilein and climfilectl are in comicegrid.h, from namelist. 

c---
      cfile = climfilein
      clab = 'current'
      goto 10

      entry readlovectl (arr, cvarin, nmonin)
      cfile = climfilectl
      clab = 'modern control'
   10 continue
c---

      cvar = cvarin

c        If fake filename given, use simple hard-coded settings, return

c???????????????????????????????????????????????????
      if (cfile.eq.'FAKE' .or. cfile.eq.'fake') then
c???????????????????????????????????????????????????

        write (ioterm,'(2a)') 'Setting fake Loveclim input for ', cvar
  
        if (firstlc) then
          do i=1,nlong
            along(i) = 2.*pi*(i-1)/nlong  ! 1st grid center at Greenwich
          enddo
          do j=1,nlatg
            alatg(j) = alatt21(j)*pi/180.
            wgauss(j) = wgausst21(j)
          enddo
          zta = 4.*pi*(radius**2) 
          do j=1,nlatg
            do i=1,nlong
              darlc(i,j) = zta * 0.5*wgauss(j) / nlong
            enddo
          enddo
          firstlc = .false.
        endif
  
        do m=1,nmonin
          do j=1,nlatg
            do i=1,nlong
c---
              zelev = -100.
              if (j.gt.nlatg/2 .and. i.gt.nlong/2) zelev = 100.
              if (j.le.nlatg/2 .and. i.le.nlong/2) zelev = 100.
c---
c             zelev = -100.
c             zy = (j -.3*nlatg)/(.7*nlatg)
c             zx = (i -.3*nlong)/(.7*nlong)
c             if (zx.gt.0. .and. zy.gt.0.) 
c    *          zelev = 800.*abs(zx-zy) + 200.*sqrt(zx**2 + zy**2)
c---
              if (cvar.eq.'surfelev') then
                arr(i,j,m) = zelev
              else if (cvar.eq.'tair2m') then
                arr(i,j,m) = -5.
              else if (cvar.eq.'precip') then
                arr(i,j,m) = 0.20
              else if (cvar.eq.'evap') then
                arr(i,j,m) = 0.01
              else if (cvar.eq.'solarflux') then
                arr(i,j,m) = 50.
              else if (cvar.eq.'tocean') then
                arr(i,j,m) = vnull_lc
                if (zelev.lt.0.) arr(i,j,m) = -0.5
              else if (cvar.eq.'icefrac') then
                arr(i,j,m) = 0.1
              else if (cvar.eq.'oceanfrac') then
                arr(i,j,m) = 0.
                if (zelev.lt.0.) arr(i,j,m) = 1.
              else 
                write (ioterm,'(/2a)') 
     *            '*** unknown fake Loveclim variable ', cvar
                stop
            endif
          enddo
        enddo
      enddo

      return

c??????????
      endif
c??????????

c        Open existing netcdf data file and read field

      write (6,'(4a)') 
     *  'reading ', clab(1:lenchr(clab)), ' Loveclim input file for ', 
     *  cvar
      rcode = nf_open(cfile, nf_nowrite, ncid)
      call checkrcode (rcode, cfile, 'open', 6000)

c        First call, set Loveclim (T21) lons, lats, gaussian weights,
c        area of cells (in common in comicegrid.h). 

c======================
      if (firstlc) then
c======================
        start(1) = 1
        count(1) = nlong
        start(2) = 1
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lon', varid)
        call checkrcode (rcode, cfile, 'lon', 7750)
        rcode = nf_get_vara_real (ncid, varid, start, count, var4)
        call checkrcode (rcode, cfile, 'lon', 7751)
        do i=1,nlong
          along(i) = var4(i)*pi/180.
        enddo

        start(1) = 1
        count(1) = nlatg
        start(2) = 1
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lat', varid)
        call checkrcode (rcode, cfile, 'lat', 7760)
        rcode = nf_get_vara_real (ncid, varid, start, count, var4)
        call checkrcode (rcode, cfile, 'lat', 7761)
        do j=1,nlatg
          alatg(j) = var4(j)*pi/180.
          wgauss(j) = wgausst21(j)
        enddo

c           Set area of each loveclim cell (like darea for ice model).
c           Sum of darlc's over all i,j = total global area (zta).

        zta = 4.*pi*(radius**2) 
        do j=1,nlatg
          do i=1,nlong
            darlc(i,j) = zta * 0.5*wgauss(j) / nlong
          enddo
        enddo

        firstlc = .false.
c==========
      endif
c==========

c        Read data field

      start(1) = 1
      count(1) = nlong
      start(2) = 1
      count(2) = nlatg
      start(3) = 1
      count(3) = nmonin

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfile, cvar, 7770)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfile, cvar, 7771)

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 8000)

      do im=1,nmonin
        do j=1,nlatg
          do i=1,nlong
            arr(i,j,im) = arrin(i,j,im)
          enddo
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine globcorrect (prge, evge, pmon, emon, h)

c     Does domain-wide correction to mean precip, evap on ice surfaces, 
c     to (approx) conserve water between Loveclim and ice model. 
c     Called from climate_gcm.

c     If first call, compute ratios of average precip and
c     evap over Loveclim's ice area, vs. bilin interp fields
c     over ice model's ice area (ratiopr, ratioev, saved). 
c     For subsequent calls, do simple domain-wide correction by 
c     multiplying interp fields uniformly by this ratio. (This ignores 
c     issues due to changes in ice area during the ice model run).

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  prge(nmon,nlong,nlatg),  evge(nmon,nlong,nlatg),
     *  pmon(nmon,nx,ny),        emon(nmon,nx,ny),
     *  h(nx,ny)

      dimension 
     *  ficelc(nlong,nlatg), icedomainlc(nlong,nlatg),
     *  inearlc(nx,ny), jnearlc(nx,ny)

      dimension ndaypm(12)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm

      save ratiopr, ratioev

      logical firstnormlc
      save firstnormlc
      data firstnormlc /.true./

c==========================
      if (firstnormlc) then
c==========================

c          Read LC ice fraction on its global grid (ficelc)

        call readloveclim (ficelc, 'icefrac', 1)  

c          Set LC sub-area containing ice model domain (icedomainlc),
c          using LC indices of boxes that contain ice grid cell centers 

        call setintnear (along, alatg, nlong, nlatg, 
     *                   alond, alatd, nx, ny,
     *                   inearlc, jnearlc)

        call izero (icedomainlc, nlong*nlatg) 
        do j=1,ny
          do i=1,nx
            icedomainlc (inearlc(i,j), jnearlc(i,j)) = 1 
          enddo
        enddo

c          Calculate total precip, evap on ice fractions in LC grid 
         
        totprlc = 0.
        totevlc = 0.
        totarlc = 0.
        do j=1,nlatg
          do i=1,nlong
            if (icedomainlc(i,j).eq.1) then
              do m=1,12
                totprlc = totprlc  + prge(m,i,j)*ficelc(i,j)
     *                               * darlc(i,j) * ndaypm(m)/365.
                totevlc = totevlc  + evge(m,i,j)*ficelc(i,j)
     *                               * darlc(i,j) * ndaypm(m)/365.
              enddo
              totarlc = totarlc + ficelc(i,j)*darlc(i,j)
            endif
          enddo
        enddo

c          Calculate total precip, evap on ice points (h>0) in ice grid
         
        totprim = 0.
        totevim = 0.
        totarim = 0.
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              do m=1,12
                totprim = totprim  + pmon(m,i,j)
     *                               * darea(i,j) * ndaypm(m)/365.
                totevim = totevim  + emon(m,i,j)
     *                               * darea(i,j) * ndaypm(m)/365.
              enddo
              totarim = totarim + darea(i,j)
            endif
          enddo
        enddo

        if (totprim.ne.0. and. totprlc.ne.0.) then 
          ratiopr = totprlc / totprim
        else
          ratiopr = 1.
        endif

        if (totevim.ne.0. and. totevlc.ne.0.) then 
          ratioev = totevlc / totevim
        else
          ratioev = 1.
        endif

        write ( ioterm, '(/ a / a,2i10, 2(/ a,3f10.4))' )
     *    'Loveclim-to-Ice model precip and evap corrections:',
     *    '  ice area, km^2 (LC,IM)         :',
     *         nint(totarlc/1.e6), nint(totarim/1.e6),
     *    '  mean precip, m/y (LC,IM), ratio:',
     *         totprlc/max(totarlc,.001),  totprim/max(totarim,.001),
     *         ratiopr,
     *    '  mean evap,   m/y (LC,IM), ratio:',
     *         totevlc/max(totarlc,.001),  totevim/max(totarim,.001),
     *         ratioev

        firstnormlc = .false.

c=========
      else
c=========

c        Apply saved ratios uniformly to interpolated precip, evap.
c        Apply over whole ice model domain, even though only calculated
c        over ice area (and at beginning of run, so ignores issues with
c        changing ice area during the ice model run)

        do j=1,ny
          do i=1,nx
            do m=1,12
              pmon(m,i,j) = pmon(m,i,j) * ratiopr
              emon(m,i,j) = emon(m,i,j) * ratioev
            enddo
          enddo
        enddo

c==========
      endif
c==========

      return
      end

c----------------------------------------------------------------------

      subroutine setintnear (along, alatg, nlong, nlatg,
     *                       alond, alatd, nx, ny,
     *                       indi, indj)

c        Like setinterp, except sets gcm grid box containing each 
c        ice grid  cell center. Returns gcm indices (indi,indj).
c
c        along,alatg supplied in radians, alond,alatd in degrees.
c        Works for any range (cut) of along - just must increase 
c        eastward, and be uniformly spaced.

      dimension
     *  along(nlong),  alatg(nlatg),
     *  alond(nx,ny),  alatd(nx,ny),
     *  indi(nx,ny),   indj(nx,ny)

      parameter (pi=3.14159265358979, tpi=(2.*pi))

c        Define range of gcm longitudes (starting from along(1), 
c        increasing through 2*pi, regardless of other along's)

      dlong = tpi/nlong
      zlong1 = along(1) - 0.5*dlong
      zlong2 = zlong1 + tpi

c        Bilinear interp...calculate indices into gcm grid

      do j=1,ny
        do i=1,nx
          zlat = alatd(i,j) * pi/180. 
          zlon = alond(i,j) * pi/180.

c         shift ice grid longitude into gcm range [zlong1,zlong2]
          do ishif = 1,10
            if (zlon.lt.zlong1) then 
              zlon = zlon + 2.*pi
            else if (zlon.gt.zlong2) then 
              zlon = zlon - 2.*pi
            else
              goto 50
            endif
          enddo
          write (6,*) 'should not get here: (setintnear, ishif loop)'
          stop
   50     continue

c         compute gcm longitude box index assuming uniform spacing:
          val = (zlon-zlong1)/dlong
          ia = int(val) + 1
          ia = mod (ia -1 + nlong, nlong) + 1
          indi(i,j) = ia

c         compute gcm latitude box index:  
          val = zlat
          do jo = 1,nlatg-1
            if (val.lt.0.5*(alatg(jo)+alatg(jo+1))) then
              indj(i,j) = jo
              goto 34
            endif
          enddo
          indj(i,j) = nlatg
   34     continue

        enddo
      enddo

c no longer: LC now changed to S->N (also see setinterp, iceclimgcm.F):
c#if defined (LOVECLIP)
c      do j=1,ny
c        do i=1,nx 
c          indj(i,j) = nlatg - indj(i,j) + 1       ! switch N->S to S->N
c        enddo
c      enddo
c#endif

      return
      end

c----------------------------------------------------------------------

      subroutine writeloveclim (h, hs, hb, fedge, maskwater, 
     *                          runliqav, runfrozav, timefluxprev,
     *                          sealev, timeice, nyearstart)

c        Write new Netcdf file containing ice-model output fields, 
c        tiling to Loveclim grid, to be read by Loveclim at start of 
c        its next asynchronous segment. Fields written are:
c          falc  = fraction of LC cell occupied by ice-model domain
c          zslc  = mean sfc elevation of ice-model (all sfc types) (m)
c          zilc  = mean sfc elevation of ice-model ice (m)
c          hilc  = mean thickness of ice-model ice (m)
c          filc  = fraction of LC cell occupied by ice-model ice
c          rullc = liquid runoff from i.s., av over LC cell (m.w.e./yr)
c          ruflc = ice discharge from i.s., av over LC cell (m.w.e./yr)

c        For discharge from ice sheet ("runoff" runliqav, runfrozav):
c          (1) route to ice-model ocean or domain edge, on ice grid.
c          (2) tile to LC grid.
c          (3) route to LC ocean, on LC grid.
c        runliqav is from oceanmelt, facemelt, baseperc, basefrml.
c        runfrozav is from calving and cliff failure. Both have
c        been accumulated in icetabbud, corrected for fedge and
c        zero-h overshoots, since timefluxprev (which for LOVECLIP,
c        is the beginning of the run (chunk)).  

#include <comicephys.h>
#include <comicegrid.h>
#include <comiceloveclim.h>

c     passed:
      dimension
     *  h(nx,ny),          hs(nx,ny),         hb(nx,ny), 
     *  fedge(nx,ny),      maskwater(nx,ny),
     *  runliqav(nx,ny),   runfrozav(nx,ny)

c     local, ice grid: 
      dimension 
     *  runliqavo(nx,ny),    runfrozavo(nx,ny),
     *  hypot(nx,ny),        indroute(nx,ny),
     *  inearlc(nx,ny),      jnearlc(nx,ny),
     *  work(nx,ny)

c     local, LC grid: 
      dimension
     *  falc(nlong,nlatg),
     *  zslc(nlong,nlatg),    zslc_old(nlong,nlatg),
     *  zilc(nlong,nlatg),    folc_old(nlong,nlatg),
     *  hilc(nlong,nlatg),    filc(nlong,nlatg),     
     *  rullc(nlong,nlatg),   ruflc(nlong,nlatg),
     *  rullco(nlong,nlatg),  ruflco(nlong,nlatg),
     *  arealc(nlong,nlatg),  areaicelc(nlong,nlatg),
     *  hypotlc(nlong,nlatg), indroutelc(nlong,nlatg),
     *  worklc(nlong,nlatg)

c     for printmaplc output only:
      dimension masklc(nlong,nlatg)

c++++++++++++++++++++++++
c     for Netcdf writing:
c++++++++++++++++++++++++
#include <netcdf.inc>
      integer*4
     * lenattr, nx4, ny4, nlev4, ival4, jval4,
     * varid, rcode, ncerr, ncid,
     * xid1, yid1, xid0, yid0, levid, timid, minmaxid,
     * minmax(2), ndim, idim(4), start(4), count(4),
     * n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls.
     * ncglobal4, ncunlim4,
     * nf_float4, nf_short4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_short4 = nf_short,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)
      save varid, rcode, ncid,
     *     xid1, yid1, xid0, yid0, levid, timid, minmaxid,
     *     minmax, idim, start, count,
     *     n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13
      data minmax /1,2/
      data  n1, n2, n3, n4, n5, n6, n7, n8, n9,n10,
     *     n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,
     *     n21,n22,n23,n24,n25,n26,n27,n28,n29,n30
     *     / 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     *      11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     *      21, 22, 23, 24, 25, 26, 27, 28, 29, 30/
      character cruntitle*80, cvar*24, cvarl*80, cunits*24
      dimension ztmp(20000)
c++++++++++++++++++++++++

c        Set hydraulic potential for ice-model routing
c        = bed elevation + (rhoice/rholiq)*ice thickness
c        (see comments in routerun)

      do j=1,ny
        do i=1,nx
          hypot(i,j) = hb(i,j) + (rhoice/rholiq)*h(i,j)
        enddo
      enddo

c        Find ocn+lakes and route ice-model runoff to ice-model ocean  
c        (moves runliqav to runliqavo, runfrozav to runfrozavo) 

      call routerun (hypot, indroute,
     *               runliqav, runfrozav, runliqavo, runfrozavo,  
     *               darea, sealev, 0., .false., nx, ny) 

c        Calculate fraction of each LC cell occupied by IM domain
c        (This is almost exact, unlike approximate arealc:darlc below).

      call calcfracim (falc)

c        Calculate LC box indices containing IM grid centers

      call setintnear (along, alatg, nlong, nlatg, 
     *                 alond, alatd, nx, ny,
     *                 inearlc, jnearlc)

c        Simple tiling to LC grid, ignoring overlaps of ice boxes 
c        across LC borders. For intrinsic quantities like elevation,
c        ice thickness, etc, weight by appropriate area of each IM cell,
c        and accumulate and normalize by that area (to not distort the
c        intrinsic value, regardless of whether each IM cell is entirely
c        within the LC cell or not).

      call zero (zslc,  nlong*nlatg)
      call zero (zilc,  nlong*nlatg)
      call zero (hilc,  nlong*nlatg)
      call zero (filc,  nlong*nlatg)
      call zero (rullc, nlong*nlatg)
      call zero (ruflc, nlong*nlatg)
      call zero (arealc,nlong*nlatg)
      call zero (areaicelc,nlong*nlatg)

c     ztott = max (.001, timeice-nyearstart)
      ztott = max (.001, timeice-timefluxprev)
      do j=1,ny 
        do i=1,nx 
          ic = inearlc(i,j) 
          jc = jnearlc(i,j) 
          zda = darea(i,j)
          if (h(i,j).gt.0.) then
            zilc(ic,jc) = zilc(ic,jc) + fedge(i,j)*hs(i,j)*zda
            hilc(ic,jc) = hilc(ic,jc) + fedge(i,j)*h(i,j)*zda
            areaicelc(ic,jc) = areaicelc(ic,jc) + fedge(i,j)*zda
          endif
          zslc(ic,jc)  = zslc(ic,jc)  + hs(i,j)*zda
          rullc(ic,jc) = rullc(ic,jc) + rhor*(runliqavo(i,j)/ztott)*zda
          ruflc(ic,jc) = ruflc(ic,jc) + rhor*(runfrozavo(i,j)/ztott)*zda
          arealc(ic,jc) = arealc(ic,jc) + zda
        enddo
      enddo

c        Set final LC arrays (and set falc to null if arealc = 0)

      do jc=1,nlatg 
        do ic=1,nlong
          if (arealc(ic,jc).gt.0.) then
            zslc(ic,jc) =  zslc(ic,jc) / arealc(ic,jc)
            zilc(ic,jc) =  zilc(ic,jc) / max (areaicelc(ic,jc), 1.e-6)
            hilc(ic,jc) =  hilc(ic,jc) / max (areaicelc(ic,jc), 1.e-6)
            filc(ic,jc) = (areaicelc(ic,jc)/arealc(ic,jc)) * falc(ic,jc)
            rullc(ic,jc) =  rullc(ic,jc) / darlc(ic,jc)
            ruflc(ic,jc) =  ruflc(ic,jc) / darlc(ic,jc)
          else
            falc(ic,jc)  = vnull_lc
            zslc(ic,jc)  = vnull_lc
            zilc(ic,jc)  = vnull_lc
            hilc(ic,jc)  = vnull_lc
            filc(ic,jc)  = vnull_lc
            rullc(ic,jc) = vnull_lc
            ruflc(ic,jc) = vnull_lc
          endif
        enddo
      enddo

c        Set hydraulic potential for LC routing. Assume oceans are
c        defined by folc_old > 0.5 (duplicate that by hypotlc=-1.).
c        For non-ocean, don't try to include ice vs bedrock (as above
c        for ice grid), just use suface elevation.

      sealev_lc =  0.            ! assume Loveclim sea level is always 0
      call readloveclim (zslc_old, 'surfelev',  1)  
      call readloveclim (folc_old, 'oceanfrac', 1)  
      do jc=1,nlatg
        do ic=1,nlong
c.........................................
          if (folc_old(ic,jc).gt.0.5) then
c.........................................
            hypotlc(ic,jc) = sealev_lc - 1.                ! LC oceans
c.............
          else
c.............
            if (falc(ic,jc).ne.vnull_lc) then              ! LC land/ice
              zslc_rou = (1.-falc(ic,jc))*zslc_old(ic,jc)
     *                 +     falc(ic,jc) *zslc(ic,jc)
            else
              zslc_rou = zslc_old(ic,jc)
            endif
            hypotlc(ic,jc) = max (zslc_rou, sealev_lc + 1.)
c..............
          endif
c..............
        enddo
      enddo

c        Find ocn+lakes and route LC ice-model-derived runoff to LC ocn 
c        (moves runllc to runllco, ruflc to ruflco) 

      call routerun (hypotlc, indroutelc,
     *               rullc, ruflc, rullco, ruflco,
     *               darlc, sealev_lc, vnull_lc, .true., nlong, nlatg) 

c        Do simple multiplicative global correction for runoffs
c        to ensure global total on final Loveclim grid is equal 
c        to total on ice-model grid. (cf. equivalent global corrections 
c        to P and E at start of ice-model segment (Loveclim->ice), 
c        done in subr globcorrect called from climate_gcm).

c     ice grid totals:
      totrul = 0.
      totruf = 0.
      do j=1,ny 
        do i=1,nx 
          totrul = totrul + rhor*(runliqav(i,j) /ztott)*darea(i,j)
          totruf = totruf + rhor*(runfrozav(i,j)/ztott)*darea(i,j)
        enddo
      enddo

c     loveclim grid totals:
      totrulo_lc = 0.
      totrufo_lc = 0.
      do jc=1,nlatg
        do ic=1,nlong
          if (rullco(ic,jc).ne.vnull_lc) 
     *      totrulo_lc = totrulo_lc + rullco(ic,jc)*darlc(ic,jc)
          if (ruflco(ic,jc).ne.vnull_lc) 
     *      totrufo_lc = totrufo_lc + ruflco(ic,jc)*darlc(ic,jc)
        enddo
      enddo

      ratioliq  = totrul / max (totrulo_lc, 1.e-10)
      ratiofroz = totruf / max (totrufo_lc, 1.e-10)

      write ( ioterm, '(/a, 2(/a,2e12.4,f10.4))' )
     *  'Ice-to-Loveclim model runoff corrections:',
     *  '  runoff_liquid, m^3/y (IM,LC), ratio:',
     *       totrul, totrulo_lc, ratioliq,
     *  '  runoff_frozen, m^3/y (IM,LC), ratio:',
     *       totruf, totrufo_lc, ratiofroz

      if (totrulo_lc.gt.1.e-10) then
        do jc=1,nlatg
          do ic=1,nlong
            if (rullco(ic,jc).ne.vnull_lc)
     *        rullco(ic,jc) = rullco(ic,jc) * ratioliq
          enddo
        enddo
      endif
      if (totrufo_lc.gt.1.e-10) then
        do jc=1,nlatg
          do ic=1,nlong
            if (ruflco(ic,jc).ne.vnull_lc)
     *        ruflco(ic,jc) = ruflco(ic,jc) * ratiofroz
          enddo
        enddo
      endif

c        Open new nc file (overwriting old), write LC arrays.
c        (filename climfileout is in comiceloveclim.h, from namelist).

      rcode = nf_create (climfileout, nf_clobber4, ncid)
      write (6,*) 'Opening new output file for Loveclim'

c        Define title, grid type, etc.

      cruntitle = 'icesheet model output for Loveclim'
      lenattr = lenchr(cruntitle)
      rcode = nf_put_att_text (ncid, ncglobal4, 'title',
     *                           lenattr, cruntitle)

      rcode = nf_put_att_text (ncid, ncglobal4, 'horiz_grid',
     *                         n7, 'regular')
      rcode = nf_put_att_text (ncid, ncglobal4, 'projection',
     *                         n6, 'latlon')

c        Define and write x dimension (longitude)

      nx4 = nlong
      rcode = nf_def_dim (ncid, 'lon', nx4, xid1)
      rcode = nf_def_var (ncid, 'lon', nf_float4, n1, xid1, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                         n9, 'longitude')
      rcode = nf_put_att_text (ncid, varid, 'units',
     *                         n12, 'degrees_east')
      rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                         n4, 'f8.3')
      rcode = nf_enddef (ncid)
      do i=1,nlong
        ztmp(i) = along(i)*180./pi
      enddo
      rcode = nf_put_vara_double (ncid, varid, n1, nx4, ztmp)
      rcode = nf_redef (ncid)

c        Define and write y dimension (latitude)

      ny4 = nlatg
      rcode = nf_def_dim (ncid, 'lat', ny4, yid1)
      rcode = nf_def_var (ncid, 'lat', nf_float4, n1, yid1, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n8, 'latitude')
      rcode = nf_put_att_text (ncid, varid, 'units',
     *                         n13, 'degrees_north')
      rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                         n4, 'f8.3')
      rcode = nf_enddef (ncid)
      do j=1,nlatg
        ztmp(j) = alatg(j)*180./pi
      enddo
      rcode = nf_put_vara_double (ncid, varid, n1, ny4, ztmp)
      rcode = nf_redef (ncid)

c        Define variable dimensions and attributes

c---------------
c        For 1D:
c---------------

c     longitude and latitude were set above as dimensions(?)

c     idim(1) = xid1
c     cvar = 'lon'
c
c     cvarl = 'longitude'
c     cunits = 'degrees_east'
c     ival4 = lenchr(cvarl)
c     jval4 = lenchr(cunits)
c     rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
c     rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
c     rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
c     rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c     ndim    = 1
c     idim(1) = yid1
c
c     cvar = 'lat'
c     cvarl = 'latitude'
c     cunits = 'degrees_north'
c     ival4 = lenchr(cvarl)
c     jval4 = lenchr(cunits)
c     rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
c     rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
c     rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
c     rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c---------------
c        For 2D:
c---------------

      ndim    = 2
      idim(1) = xid1
      idim(2) = yid1

      cvar = 'darea_lc'
      cvarl = 'cell area'
      cunits = 'm^2'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'zsurf_im'
      cvarl = 'elevation of ice-model surface'
      cunits = 'm'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'farea_im'
      cvarl = 'fractional area of ice-model grid'
      cunits = '0 to 1'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'fice_im'
      cvarl = 'fractional area of ice-model ice'
      cunits = '0 to 1'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'hice_im'
      cvarl = 'thickness of ice-model ice'
      cunits = 'm'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'zice_im'
      cvarl = 'elevation of ice-model ice'
      cunits = 'm'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'runoff_liq'
      cvarl = 'liquid discharge, ice to ocean'
      cunits = 'm/y liquid equivalent'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'runoff_froz'
      cvarl = 'frozen discharge, ice to ocean'
      cunits = 'm/y liquid equivalent'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'index_route'
      cvarl = 'index to routing ocean point'
      cunits = 'ilon+nlon(ilat-1)'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_int, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n3,'i10')

c        Leave define mode

      rcode = nf_enddef (ncid)

c        Write variables

c~~~~~~~~~~~~~~~
c        For 1D:
c~~~~~~~~~~~~~~~

c     longitude and latitude were set above as dimensions(?)

c     start(1) = 1
c     count(1) = nlong 
c     do i=1,nlong
c       ztmp(i) = along(i)*180./pi
c     enddo
c     rcode = nf_inq_varid (ncid, 'lon', varid)
c     rcode = nf_put_vara_double (ncid, varid, start, count, along)

c     start(1) = 1
c     count(1) = nlatg
c     do j=1,nlatg
c       ztmp(j) = alatg(j)*180./pi
c     enddo
c     rcode = nf_inq_varid (ncid, 'lat', varid)
c     rcode = nf_put_vara_double (ncid, varid, start, count, ztmp)

c~~~~~~~~~~~~~~~
c        For 2D:
c~~~~~~~~~~~~~~~

      start(1) = 1
      count(1) = nlong
      start(2) = 1
      count(2) = nlatg

      rcode = nf_inq_varid (ncid, 'darea_lc', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, darlc)

      rcode = nf_inq_varid (ncid, 'zsurf_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, zslc)

      rcode = nf_inq_varid (ncid, 'farea_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, falc)

      rcode = nf_inq_varid (ncid, 'fice_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, filc)

      rcode = nf_inq_varid (ncid, 'hice_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hilc)

      rcode = nf_inq_varid (ncid, 'zice_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, zilc)

      rcode = nf_inq_varid (ncid, 'runoff_liq', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, rullco)

      rcode = nf_inq_varid (ncid, 'runoff_froz', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, ruflco)

      rcode = nf_inq_varid (ncid, 'index_route', varid)
      rcode = nf_put_vara_int (ncid, varid, start, count, indroutelc)

c        Return to define mode (not necessary?) and close

      rcode = nf_redef (ncid)
      rcode = nf_close (ncid)

c        Ascii output (diagnostic only)

      iu = 154

c     land-ocean mask for reverse-video outlines:
      do j=1,nlatg
        do i=1,nlong
          if (folc_old(i,j).gt.0.5) then
             masklc(i,j) = 1
          else    
             masklc(i,j) = 0
          endif   
          worklc(i,j) = masklc(i,j)
        enddo   
      enddo
      call printmaplc (timeice, worklc, 'masklc', 0.1, masklc, iu)

      call printmaplc (timeice, folc_old, 'folc_old', 0.1, masklc, iu)

      call printmaplc (timeice, zslc, 'zslc', 100., masklc, iu)

      call printmaplc (timeice, falc, 'falc', 0.1, masklc, iu)

      call readloveclim (worklc, 'icefrac', 1)  
      do j=1,nlatg
        do i=1,nlong
          if (falc(i,j).ne.vnull_lc) 
     *      worklc(i,j) = (1.-falc(i,j))*worklc(i,j)
     *                  +     falc(i,j) *filc(i,j)
        enddo
      enddo
      call printmaplc (timeice, worklc, 'filc_all', 0.1, masklc, iu)

      call zero (worklc, nlong*nlatg)
      do j=1,nlatg
        do i=1,nlong
          if (masklc(i,j).eq.0) then
            jo = (indroutelc(i,j)-1)/nlong + 1
            io = indroutelc(i,j) - nlong*(jo-1)
            worklc(io,jo) = 1.
          endif
        enddo
      enddo
      call printmaplc (timeice, worklc, 'indroutelc', 0.1, masklc, iu)
      call printmaplc (timeice, rullc,  'rullc',  0.1, masklc, iu)
      call printmaplc (timeice, rullco, 'rullco', 0.1, masklc, iu)
      call printmaplc (timeice, ruflc,  'ruflc',  0.1, masklc, iu)
      call printmaplc (timeice, ruflco, 'ruflco', 0.1, masklc, iu)

      call flush(iu)

      iu = 153

      call zero (work, nx,ny)
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.0) then
            jo = (indroute(i,j)-1)/nx + 1
            io = indroute(i,j) - nx*(jo-1)
            work(io,jo) = 1.
          endif
        enddo
      enddo
      call printmap (timeice, work, 'indroute', 0.1, iu, 0)
      call printmap (timeice, runliqav,   'runliqav',   0.1, iu, 0)
      call printmap (timeice, runliqavo,  'runliqavo',  0.1, iu, 0)
      call printmap (timeice, runfrozav,  'runfrozav',  0.1, iu, 0)
      call printmap (timeice, runfrozavo, 'runfrozavo', 0.1, iu, 0)

      call flush(iu)

      return
      end

c-----------------------------------------------------------------------

      subroutine routerun (hypot, indroute,
     *                     runliqav, runfrozav, runliqavo, runfrozavo,
     *                     darea, sealev, vnull, ifwrap, nx, ny) 

c        Routes runoff (liquid and frozen separately) to ocean.
c        i.e, sets runliqavo given runliqav, runfrozavo given
c        runfrozav. Also returns indroute = i + nx*(j-1), index of
c        routed ocean point for each source point). 

c        Uses hydraulic potential for routing (hypot), 
c        = bed elevation + (rhoice/rholiq)*ice thickness, which
c        automatically accounts for basal ice pressure and whether 
c        ice would be grounded or floating (for routing purposes only, 
c        not necessarily in ice model - see comments in locatewater).
c        Also allows runliqav from baseperc,basefrml to be routed
c        under icesheet.

#include <comicephys.h>

      dimension
     *  hypot(nx,ny),     indroute(nx,ny),
     *  runliqav(nx,ny),  runfrozav(nx,ny),
     *  runliqavo(nx,ny), runfrozavo(nx,ny),
     *  darea(nx,ny)
      logical ifwrap

      dimension 
     *  maskh2o(nx,ny), sill(nx,ny), indsill(nx*ny)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 8)   ! 4 for E-W/N-S, 8 for diagonals too

c        First, find lakes, sills, ocean (maskh2o, indsill, nwater)

      call locatewater (maskh2o, sill, indsill, nwater, hypot, 
     *                  sealev, ifwrap, nx, ny)

c        Do routing to ocean, by steepest downhill path (hypot).
c        If encounter a lake, skip to its sill point.

      call izero (indroute, nx*ny)

c>>>>>>>>>>>>>>>>
      do j=1,ny
        do i=1,nx
c>>>>>>>>>>>>>>>>
          ia  = i
          ja  = j
c~~~~~~~~~~~~~~~~~~~~~~~~
          do iter=1,nx*ny
c~~~~~~~~~~~~~~~~~~~~~~~~
            if (maskh2o(ia,ja).eq.1) then             ! already in ocean
              ifin = ia
              jfin = ja
              goto 100

            else if (maskh2o(ia,ja).ge.2) then        ! in lake -> sill
              iwater = maskh2o(ia,ja)
              jb = (indsill(iwater)-1)/nx + 1
              ib = indsill(iwater) - (jb-1)*nx

            else if (maskh2o(ia,ja).eq.0) then  ! land, search neighbors
              hypotmin = hypot(ia,ja)
              do look=1,nlook                        
                if (ifwrap) then
                  ii = mod (ia+ish(look)+nx-1, nx) + 1
                else
                  ii = max (1, min (nx, ia+ish(look)))
                endif
                jj = max (1, min (ny, ja+jsh(look)))
                if (maskh2o(ii,jj).eq.1) then       ! found ocean
                  ifin = ii
                  jfin = jj
                  goto 100
                else                                ! found lake or land
                  if (maskh2o(ii,jj).ge.2) then
                    zhyp = sill(ii,jj) 
                  else
                    zhyp = hypot(ii,jj)
                  endif
                  if (zhyp.lt.hypotmin) then        ! downhill lk or lnd
                    hypotmin = zhyp
                    ib = ii
                    jb = jj
                  endif
                endif
              enddo   

              if (hypotmin.eq.hypot(ia,ja)) then   
c               if not ifwrap and stuck at a (land) domain boundary, 
c               then ok - set as final routing point
                if ( .not.ifwrap .and.
     *               (ia.eq.1.or.ia.eq.nx.or.ja.eq.1.or.ja.eq.ny) ) then
                  ifin = ia
                  jfin = ja
                  goto 100
                else
                  write (6,'(a/a,2i6,5x,l1/a,3i6)')
     *              '*** routerun: should not get here A', 
     *              'nx,ny,ifwrap=', nx, ny, ifwrap,
     *              'ia,ja,iter  =', ia, ja, iter
                  do jo = min(ja+1,ny), max(ja-1,1), -1
                    do io = max(ia-1,1), min(ia+1,nx)
                      write (6,'(3i6, f10.3)')
     *                  io, jo, maskh2o(io,jo), hypot(io,jo)  
                    enddo
                  enddo
                  stop
                endif
              endif

  200         continue    ! ok - found a lower land pt or sill, continue
            endif

            ia = ib
            ja = jb
c~~~~~~~~~~~~~~
          enddo
c~~~~~~~~~~~~~~
          write (6,'(a/a,2i6,5x,l1/a,3i6/a,2i6)')
     *      '*** routerun: should not get here B', 
     *      'nx,ny,ifwrap=', nx, ny, ifwrap,
     *      'i,j,iter  =', i, j, iter-1,
     *      'ia,ja,    =', ia, ja
          do jo = min(ja+1,ny), max(ja-1,1), -1
            do io = max(ia-1,1), min(ia+1,nx)
              if (maskh2o(io,jo).ge.1) then
                indz = indsill(maskh2o(io,jo))
                jz = (indz-1)/nx + 1
                iz = indz - nx*(jz-1)
              else
                indz = 0
                iz = 0
                jz = 0
              endif
              write (6,'(3i6,f10.3,3i6)')
     *          io, jo, maskh2o(io,jo), hypot(io,jo), indz,iz,jz
            enddo
          enddo
          stop

  100     continue       ! ok - found ocean or non-ifwrap edge
          indroute(i,j) = ifin + (jfin-1)*nx
c>>>>>>>>>>>>
        enddo
      enddo
c>>>>>>>>>>>>

c        Transfer runoff (run*av) to ocean locations (run*avo),
c        simply using routing index points (indroute), and
c        allowing for differences in cell area

      call zero (runliqavo,  nx*ny)
      call zero (runfrozavo, nx*ny)
      do j=1,ny
        do i=1,nx
          jr = (indroute(i,j)-1) / nx + 1
          ir =  indroute(i,j) - (jr-1)*nx
          if (runliqav(i,j).ne.vnull) then
            runliqavo(ir,jr) = runliqavo(ir,jr) 
     *                       + runliqav(i,j)*darea(i,j)/darea(ir,jr) 
          endif
          if (runfrozav(i,j).ne.vnull) then
            runfrozavo(ir,jr) = runfrozavo(ir,jr) 
     *                        + runfrozav(i,j)*darea(i,j)/darea(ir,jr) 
          endif
        enddo
      enddo

c        For ice model, vnull=0 so no effect here. 
c        For LC, change zero (routed) runoff to LC null value.   

      do j=1,ny
        do i=1,nx
          if (runliqavo(i,j).eq.0.)  runliqavo(i,j)  = vnull
          if (runfrozavo(i,j).eq.0.) runfrozavo(i,j) = vnull
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine locatewater (maskh2o, sill, indsill, nwater,
     *                        hypot, sealev, ifwrap, nx, ny)

c     Adapted from (unused) subr initwater in icehydrol.F, for locating
c     ocean, lake and sills, which is then used for routing of runoff
c     (called from routerun, called from writeloveclim).

c     Uses hydraulic potential hypot = hb + (rhoice/rholiq)*h
c     (see comments in routerun).

c     Just for purposes of routing to ocean, don't need to continue the
c     search for lakes within islands within lakes within continents,
c     because any runoff from those 2nd-order lakes would go over the
c     island and into a 1st-order lake, which has its sill calculated 
c     below, and from there over continental land to ocean.

c     nb: lakes, floating vs. grounded lake ice, etc, deduced here
c     do not jive (for now) with main program's physics. (Latter does 
c     not have any lakes at all unless -DWATERLAND. And if it does, uses
c     indlake, npoilake, nlake, separate from here (see findwater)).

c     Generally performs each search by searching one grid box 
c     from the current set at each iteration. Maintains list(1:nlist)
c     of indices (i,j, packed) of points added in current iteration,
c     and use only those for the next iteration. 

c     Ocean and rivers can propagate N-S, E-W and along diagonals,
c     and lakes too.

c     Returned: 
c       maskh2o = 0 for land; 1 for ocean; 2,3,4... for lakes
c       sill = sill elevation for lakes, = hypot for non-lake land
c       indsill(n) = i + nx*(j-1) index of sill for lake # n
c       nwater = number of lakes (including ocean)
c     Supplied:
c       hypot= hydraulic potential, bed elev + (rhoice/rholiq)*ice thic.
c       sealev = current eustatic sea level
c       ifwrap = .true. if global (wraparound longitude) grid
c       nx,ny = domain grid size

#include <comicephys.h>
c     nb: no comicegrid.h here. Array sizes are passed, so locatewater
c     can be used for Loveclim global grid as well as ice model grid.

      dimension
     *  maskh2o(nx,ny), sill(nx,ny),
     *  indsill(nx*ny),                        ! nx*ny >= max # of lakes
     *  hypot(nx,ny)
      logical ifwrap

c     local:
      dimension 
     *  icefloat(nx,ny), list(nx*ny), listprev(nx*ny)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 8)   ! 4 for E-W/N-S, 8 for diagonals too

c        Perform search for oceans (submerged points that are in 
c        contact with any submerged points at domain edges)

      call izero (maskh2o,  nx*ny)
      call izero (icefloat, nx*ny)

c----------------------
      do iter = 1,nx*ny
c----------------------

        nlist = 0

        if (iter.eq.1) then

c         initial ocean set is all submerged pts around domain edges:
          do j=1,ny
            iskip = max(nx-1,1)
            if (j.eq.1.or.j.eq.ny) iskip = 1
            do i=1,nx,iskip
              if (hypot(i,j).le.sealev) then 
                maskh2o(i,j) = 1
                icefloat(i,j) = 1
                nlist = nlist + 1
                list(nlist) = i + (j-1)*nx
              endif
            enddo
          enddo

        else

c         normal iteration: look at adjacent points for new ocean pts:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook
              if (ifwrap) then
                ii = mod (i+ish(look)+nx-1, nx) + 1
              else
                ii = max (1, min (nx, i+ish(look)))
              endif
              jj = max (1, min (ny, j+jsh(look)))
              if (maskh2o(ii,jj).eq.0) then 
                if (hypot(ii,jj).le.sealev) then
                  maskh2o(ii,jj) = 1
                  icefloat(ii,jj) = 1
                  nlist = nlist + 1
                  list(nlist) = ii + (jj-1)*nx
                endif
              endif
            enddo
          enddo

        endif

c          Skip out if no new ocean points found

        if (nlist.eq.0) goto 1000

c          Copy current list to previous list

        nlistprev = nlist
        do n=1,nlist
          listprev(n) = list(n)
        enddo

c----------
      enddo
c----------
 1000 continue

c        Set initial nwater, sill, indsill just for ocean
      
      call zero (sill, nx*ny)
      call izero (indsill, nx*ny)
      nwater = 0

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.1) then
            sill(i,j) = sealev
            if (nwater.eq.0) then
              nwater = 1                           ! "lake # 1" is ocean
              indsill(1) = 1                       ! arbitrary value
            endif
          endif
        enddo
      enddo

c        Perform search for lakes, by eliminating non-lakes (land pts
c        that have a downhill path to a coast).

c        Right now maskh2o=1 (ocean) or 0(not).
c        Initialize non-ocean points (maskh2o=0) to *possible* 
c        lake pts (maskh2o=-1) temporarily. 

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.0) maskh2o(i,j) = -1
        enddo
      enddo

c======================
      do iter = 1,nx*ny
c======================

        nlist = 0

        if (iter.eq.1) then

c         Set initial non-lake set (maskh2o=0) to all coastal points,
c         and, if not ifwrap, to all land pts at edges of domain 
c         (we don't know about slopes outside):
          do j=1,ny
            do i=1,nx
              if (maskh2o(i,j).eq.-1) then
                iflag = 0
                if ( .not.ifwrap .and. 
     *               (i.eq.1.or.i.eq.nx.or.j.eq.1.or.j.eq.ny) ) then
                  iflag = 1
                else
                  do look=1,nlook
                    if (ifwrap) then
                      ii = mod (i+ish(look)+nx-1, nx) + 1
                    else
                      ii = max (1, min (nx, i+ish(look)))
                    endif
                    jj = max (1, min (ny, j+jsh(look)))
                    if (maskh2o(ii,jj).eq.1) iflag = 1
                  enddo
                endif
                if (iflag.eq.1) then
                  maskh2o(i,j) = 0
                  nlist = nlist + 1
                  list(nlist) = i + (j-1)*nx
                endif
              endif
            enddo
          enddo

        else

c         normal iteration: look uphill for new non-lake points:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook
              if (ifwrap) then
                ii = mod (i+ish(look)+nx-1, nx) + 1
              else
                ii = max (1, min (nx, i+ish(look)))
              endif
              jj = max (1, min (ny, j+jsh(look)))
              if ( hypot(ii,jj).gt.hypot(i,j) .and. 
     *             maskh2o(ii,jj).eq.-1 ) then 
                maskh2o(ii,jj) = 0
                nlist = nlist + 1
                list(nlist) = ii + (jj-1)*nx
              endif
            enddo
          enddo

        endif

c          Skip out if no new non-lake points found

        if (nlist.eq.0) goto 2000

c          Copy current list to previous list

        nlistprev = nlist
        do n=1,nlist
          listprev(n) = list(n)
        enddo

c==========
      enddo
c==========
 2000 continue

c       Propagate sill-depth constraints through lakes.
c       sill already set above for ocean (maskh2o=1). 

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.-1) then
            sill(i,j) = 1.e20
          else if (maskh2o(i,j).eq.0) then 
            sill(i,j) = 0. 
          endif
        enddo
      enddo

      do iter = 1,nx*ny
        nchange = 0
        do j=1,ny
          do i=1,nx
            if (maskh2o(i,j).eq.-1) then
              sillold = sill(i,j)
              do look = 1,nlook
                if (ifwrap) then
                  ii = mod (i+ish(look)+nx-1, nx) + 1
                else
                  ii = max (1, min (nx, i+ish(look)))
                endif
                jj = max (1, min (ny, j+jsh(look)))
                if (maskh2o(ii,jj).eq.0) then
c                 land pt (cannot be ocn, all coastal pts are non-lake):
                  sill(i,j) = min (sill(i,j), hypot(ii,jj))
                else
c                 lake point:
                  sill(i,j) = min (sill(i,j), sill(ii,jj))
                endif
              enddo
c             remove "emergent land" from lake status:
              if (sill(i,j).lt.hypot(i,j)) then
                maskh2o(i,j) = 0
                icefloat(i,j) = 0
                sill(i,j) = 0.
              endif
              if (sillold.ne.sill(i,j)) then
                nchange = nchange + 1
              endif
            endif
          enddo
        enddo
        if (nchange.eq.0) goto 3000
        if (iter.gt.nx*ny-10) then 
          write (6,*) '*** Warning (locatewater): iter=',iter
        endif
      enddo
      write (6,*)'*** Error (locatewater): exceeded iteration: ',iter
      stop
 3000 continue

c        Convert -1 maskh2o pts (currently 0=land, 1=ocean, -1=inland 
c        lakes) to labels (2,3,..) that identify contiguous inland lakes

      do j=1,ny
        do i=1,nx

          if (maskh2o(i,j).eq.-1) then
            nwater = nwater + 1
            maskh2o(i,j) = nwater
            indsill(nwater) = 0

            nlistprev = 1
            listprev(1) = i + (j-1)*nx

            do iter = 1,nx*ny

              nlist = 0
              do n=1,nlistprev
                ja = (listprev(n)-1) / nx + 1
                ia =  listprev(n) - (ja-1)*nx
                do look = 1,nlook
                  if (ifwrap) then
                    ii = mod (ia+ish(look)+nx-1, nx) + 1
                  else
                    ii = max (1, min (nx, ia+ish(look)))
                  endif
                  jj = max (1, min (ny, ja+jsh(look)))
                  if ( maskh2o(ii,jj).gt.0) then
                    if (maskh2o(ii,jj).ne.maskh2o(i,j)) then
                      write (6,'(a,6i4)')
     *                  '*** Error 1 (locatewater) ii,jj,i,j,maskh2o=',
     *                  ii, jj, i, j, maskh2o(ii,jj), maskh2o(i,j)
                      stop
                    endif
                  else if (maskh2o(ii,jj).eq.-1) then
                    if (sill(ii,jj).ne.sill(i,j)) then
                      write (6,'(a,6i4/a,2f22.15)')
     *                  '*** Error 2 (locatewater) ii,jj,i,j,maskh2o=',
     *                  ii, jj, i, j, maskh2o(ii,jj), maskh2o(i,j),
     *                  '    sills=',sill(ii,jj), sill(i,j)
                      stop
                    endif
                    maskh2o(ii,jj) = maskh2o(i,j)
                    nlist = nlist + 1
                    list(nlist) = ii + (jj-1)*nx

                  else if (maskh2o(ii,jj).eq.0) then
c                   set location of first-found lake sill point (land)
                    if (hypot(ii,jj).eq.sill(i,j) .and.
     *                  indsill(nwater).eq.0) then
                      indsill(nwater) = ii + (jj-1)*nx
                    endif

                  endif
                enddo
              enddo

c                Skip out if no new lake points found

              if (nlist.eq.0) goto 4000

c                Copy current list to previous list

              nlistprev = nlist
              do n=1,nlist
                listprev(n) = list(n)
              enddo

            enddo
 4000       continue

          endif

        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine printmaplc (time, arrin, cdescin, scale, mask, iu)

c        Displays array arr to unit iu in "a1" format.
c        Like printmap (iceshow.F), but for Loveclim (LC) variables.
c        mask = 0 for land, 1 for ocean (for reverse video oocean bdry).

#include <comicegrid.h>

      dimension arrin(nlong,nlatg), mask(nlong,nlatg)
      character*(*) cdescin

      dimension arr(nlong,nlatg)
      character*1 carr(nlong)
      character cdesc*80, cbrr*20000
      logical first
      data first /.true./
      save iac, jac

      character crev*4, csul*4, coff*3
      parameter (crev = char(27)//'[7m',          ! reverse video
     *           csul = char(27)//'[4m',          ! underline (not used)
     *           coff = char(27)//'[m' )          ! reset all
c     character cbra*1
c     integer*1 icbra
c     equivalence (icbra,cbra)
c     data icbra /O'154000000000000000000'/        ! "ctrl["
c     save crev, csul, coff

      cdesc = cdescin

      if (first) then
        iac = ichar('A') - 1
        jac = ichar('a') - 1
      endif

      call scopy (nlong*nlatg, arrin, 1, arr, 1)

      ia = 1
      ib = nlong
      ic = 1

      if (.not. first) write (iu,*) 
      write(iu,100)
     *  time, cdesc(1:lenchr(cdesc)), scale
      write(iu,101) 
     *  (nint(along(i)*180./pi), i=ia+4,ib,5*ic)
  100 format(
     *  'time:',f12.1,'  ',a,'/',f11.4,'  (A-Z pos, a-z neg):')
  101 format(/6x,2x,500i5)

c>>>>>>>>>>>>>>>>>>>>>>>>>
      do j = nlatg, 1, -ic
c>>>>>>>>>>>>>>>>>>>>>>>>>

        do i=ia,ib
          if (abs(arr(i,j)).le.1.e-8) then
            carr(i) = '.'
          else
            zs = max (-1.e6, min (1.e6, arr(i,j)/scale))
            n = nint(zs)
            if (n.eq.0) then
              if (arr(i,j).eq.0.) then
                carr(i) = '.'
              else if (arr(i,j).gt.0.) then
                carr(i) = '+'
              else
                carr(i) = '-'
              endif
            else if (n.ge.1..and. n.le.26) then
              carr(i) = char(iac + n)
            else if (n.le.-1..and. n.ge.-26) then
              carr(i) = char(jac - n)
            else if (n.gt.26) then
              carr(i) = '*'
            else if (n.lt.-26) then
              carr(i) = '#'
            endif
          endif
        enddo

c........................
        if (.false.) then
c........................

c            Show ocean (not lakes) as blanks

          lenb = 0
          do i=ia,ib,ic 
            lenb = lenb + 1
            if (mask(i,j).eq.1) then
              cbrr(lenb:lenb) = ' '
            else
              cbrr(lenb:lenb) = carr(i)
            endif
          enddo
 
c...........
        else
c...........

c            Reverse-video for ocean pts touching land, using mask.
c            (nb: unlike printmap which shows land pts touching ocean)

          istatprev = 0
          lenb = 0
          do i = ia,ib,ic
            istat = 0
c           im = max (i-ic,1)                    ! walls
c           ip = min (i+ic,nlong)                ! walls
            im = mod (i-ic+nlong -1 ,nlong) + 1  ! wraparound
            ip = mod (i+ic+nlong -1 ,nlong) + 1  ! wraparound
            jm = max (j-ic,1)
            jp = min (j+ic,nlatg)
            if( mask(i,j).eq.1 .and. 
     *          (mask(im,j).ne.1.or.mask(ip,j).ne.1.or.
     *           mask(i,jm).ne.1.or.mask(i,jp).ne.1) 
     *         )  istat = 1

            if (istat.ne.istatprev) then
              cbrr(lenb+1:lenb+3) = coff
              lenb = lenb + 3
              if (istat.eq.1) then
                cbrr(lenb+1:lenb+4) = crev
                lenb = lenb + 4
              endif
            endif

            cbrr(lenb+1:lenb+1) = carr(i)
            lenb = lenb + 1

            if (i.gt.ib-ic .and. istat.ne.0) then
              cbrr(lenb+1:lenb+3) = coff
              lenb = lenb + 3
            endif

            istatprev = istat
          enddo

c............
        endif
c............

        write(iu,'(f6.1,2x,a)') alatg(j)*180./pi, cbrr(1:lenb)

c>>>>>>>>>>
      enddo
c>>>>>>>>>>

      first = .false.

      return
      end

c-----------------------------------------------------------------------

      subroutine calcfracim (falc)

c     Calculates fraction of each LC cell occupied by IM domain.
c     Sub-divides each LC cell into a fine mesh (nfile x nfine), 
c     and tests each of those lon,lat locations. 
c     Returns falc = # of locns within IM domain / (nfine*nfine)

#include <comicephys.h>
#include <comicegrid.h>

      dimension falc(nlong,nlatg)

      parameter (nfine=50)

      dimension alonge(0:nlong+1), alatge(0:nlatg+1)

c        Set extended LC grid variables, convenient for
c        setting grid interfaces (zlata,zlatb and zlona,zlonb) below

      alonge(0) = along(nlong) - 2.*pi
      call scopy (nlong,  along, 1, alonge(1), 1)
      alonge(nlong+1) = along(1) + 2.*pi

      alatge(0) = -0.5*pi
      call scopy (nlatg,  alatg, 1, alatge(1), 1)
      alatge(nlatg+1) = 0.5*pi

c=================
      do j=1,nlatg
c=================
        zlata = 0.5*(alatge(j-1)+alatge(j))
        zlatb = 0.5*(alatge(j)+alatge(j+1))
c===================
        do i=1,nlong
c===================
          zlona = 0.5*(alonge(i-1)+alonge(i))
          zlonb = 0.5*(alonge(i)+alonge(i+1))

c            Test on fine subdivision (nfine x nfine) of this LC cell

          ncountin = 0 
c---------------------------
          do jfine = 1,nfine
c---------------------------
            zlat = zlata + (zlatb-zlata)*((jfine-0.5)/nfine)
c-----------------------------
            do ifine = 1,nfine
c-----------------------------
              zlon =  zlona + (zlonb-zlona)*((ifine-0.5)/nfine)

#if defined (LONLAT)
              zlatd = zlat*180./pi
              zlond = zlon*180./pi
c             assume ice domain longitudes (alon1,alon2, comicegrid.h)
c             are in range -180 to 180:
              if (zlond.lt.-180.) zlond = zlond + 360.
              if (zlond.gt. 180.) zlond = zlond - 360.
              if (zlond.ge.alon1 .and. zlond.le.alon2 .and.
     *            zlatd.ge.alat1 .and. zlatd.le.alat2) then
                ncountin = ncountin + 1
              endif

#elif defined (STEREO)
              if ( (stdparallel.ge.0. .and. zlat.ge.0.) .or.
     *             (stdparallel.le.0. .and. zlat.le.0.) ) then
                call ps_xy (zx, zy, zlat, zlon, xoffa, yoffa)

                if ( zx. ge. -0.5*nx*dx0 .and.
     *               zx. le.  0.5*nx*dx0 .and.
     *               zy. ge. -0.5*ny*dy0 .and.
     *               zy. le.  0.5*ny*dy0 ) then
                  ncountin = ncountin + 1
                endif
              endif
#endif

c----------------
            enddo
          enddo
c----------------
            
          falc(i,j) = float(ncountin) / (nfine*nfine)

c============
        enddo
      enddo
c============

      return
      end

c-----------------------------------------------------------------------

c@@@@@
#endif
c@@@@@
