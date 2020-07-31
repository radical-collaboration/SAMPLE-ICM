c@@@@@@@@@@@@@@@@@@@@@
#if defined (LOVECLIP)
c@@@@@@@@@@@@@@@@@@@@@

c-----------------------------------------------------------------------

      subroutine readloveclim (arr, cvarin, nlonin, nlatin, nmonin)

c     Reads nc file written by Loveclim before each ice model segment 
c     of coupled runs, for one field, returned in arr on LC grid.
c     Grid size is passed (nlonin,nlatin) - either atmos (nlong,nlatg) 
c     or ocean (nlongo,nlatgo), in comicegrid.h.
c     Netcdf field names (cvarin) are: 
c       surfelev  = surface elevation (m)
c       tair2m    = 2-m surface air temperature (C)
c       precip    = precipitation (m/y liquid equiv)
c       evap      = evaporation (m/y liquid equiv)
c       solarflux = incident downward solar flux (W/m2), sfc or TOA(?)
c       tocean    = ocean temperature at ~400 m depth (C)
c     nlonin = longitude grid dimension
c     nlatin = latitude grid dimension
c     nmonin = number of months to be read (1 or 12)
c     Null data value is vnull_lc (set in comiceloveclim.h)

c     Can be called for anomaly method, using entry "readlovectl"(below)

#include <comicephys.h>
#include <comicegrid.h>
#include <comiceloveclim.h>
#include <netcdf.inc>

      dimension arr(nlonin,nlatin,nmonin)
      character*(*) cvarin

      character cvar*24
      integer*4 ncid, rcode, varid, start(3), count(3)

      real*4 arrin(nlonin,nlatin,12), 
     *       var4(nlong+nlatg), 
     *       var4o(nlongo+nlatgo)

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

      entry readlovectl (arr, cvarin, nlonin, nlatin, nmonin)
      cfile = climfilectl
      clab = 'modern control'
   10 continue
c---

      cvar = cvarin

c        If fake filename given, use simple hard-coded settings

c???????????????????????????????????????????????????
      if (cfile.eq.'FAKE' .or. cfile.eq.'fake') then
c???????????????????????????????????????????????????

        write (ioterm,'(2a)') 'Setting fake Loveclim input for ', cvar
  
c........................
        if (firstlc) then
c........................
c         atmos (in comicegrid.h):
          do i=1,nlong
            along(i) = 2.*pi*(i-1)/nlong  ! 1st grid center at Greenwich
          enddo
          do j=1,nlatg
            alatg(j) = -0.5*pi + (j-.5)*pi/nlatg
          enddo

c         ocean (in comicegrid.h):
          do i=1,nlongo
            alongo(i) = 2.*pi*(i-1)/nlongo ! 1st grd center at Greenwich
          enddo
          do j=1,nlatgo
            alatgo(j) = -0.5*pi + (j-.5)*pi/nlatgo
          enddo
c............
        endif
c............
  
        do m=1,nmonin
          do j=1,nlatin
            do i=1,nlonin
c---
              zelev = -100.
              if (j.gt.nlatin/2 .and. i.gt.nlonin/2) zelev = 100.
              if (j.le.nlatin/2 .and. i.le.nlonin/2) zelev = 100.
c---
c             zelev = -100.
c             zy = (j -.3*nlatin)/(.7*nlatin)
c             zx = (i -.3*nlonin)/(.7*nlonin)
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
              else 
                write (ioterm,'(/2a)') 
     *            '*** unknown fake Loveclim variable ', cvar
                stop
            endif
          enddo
        enddo
      enddo

c?????????
      else
c?????????

c        Open existing netcdf data file and read field

      write (6,'(4a)') 
     *  'reading ', clab(1:lenchr(clab)), ' Loveclim input file for ', 
     *  cvar
      rcode = nf_open(cfile, nf_nowrite, ncid)
      call checkrcode (rcode, cfile, 'open', 6000)

c        First call, set Loveclim atmos and ocean lons, lats 
c        (in common in  comicegrid.h)

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
        call checkrcode (rcode, cfile, 'lat', 7752)
        rcode = nf_get_vara_real (ncid, varid, start, count, var4)
        call checkrcode (rcode, cfile, 'lat', 7753)
        do j=1,nlatg
          alatg(j) = var4(j)*pi/180.
        enddo

        start(1) = 1
        count(1) = nlongo
        start(2) = 1
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lono', varid)
        call checkrcode (rcode, cfile, 'lono', 7754)
        rcode = nf_get_vara_real (ncid, varid, start, count, var4o)
        call checkrcode (rcode, cfile, 'lono', 7755)
        do i=1,nlongo
          alongo(i) = var4o(i)*pi/180.
        enddo

        start(1) = 1
        count(1) = nlatgo
        start(2) = 1
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lato', varid)
        call checkrcode (rcode, cfile, 'lato', 7756)
        rcode = nf_get_vara_real (ncid, varid, start, count, var4o)
        call checkrcode (rcode, cfile, 'lato', 7757)
        do j=1,nlatgo
          alatgo(j) = var4o(j)*pi/180.
        enddo
c==========
      endif
c==========

c        Read data field

      start(1) = 1
      count(1) = nlonin
      start(2) = 1
      count(2) = nlatin
      start(3) = 1
      count(3) = nmonin

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfile, cvar, 7770)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfile, cvar, 7771)

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 8000)

      do im=1,nmonin
        do j=1,nlatin
          do i=1,nlonin
            arr(i,j,im) = arrin(i,j,im)
          enddo
        enddo
      enddo

c??????????
      endif  ! fake or not
c??????????

c       If first call, set grid cell areas (atmos and ocean, 
c       fake or not, in comicegrid.h)

c>>>>>>>>>>>>>>>>>>>>>>
      if (firstlc) then
c>>>>>>>>>>>>>>>>>>>>>>

c       atmos:
        ztot = 0.
        do j=1,nlatg
          jp = min (j+1,nlatg)
          jm = max (j-1,1)
          zda = 2.*pi*(radius**2) * cos(alatg(j))
     *          * 0.5*(alatg(jp)-alatg(jm))
          if (j.eq.1 .or. j.eq.nlatg) zda = 2.*zda  ! no grid pt at pole
          do i=1,nlong
            darlc(i,j) = zda/nlong
            ztot = ztot + darlc(i,j)
          enddo
        enddo

c       adjust global total:
        zadj = 4.*pi*(radius**2) / ztot
        do j=1,nlatg
          do i=1,nlong
            darlc(i,j) = darlc(i,j) * zadj
          enddo
        enddo

c       ocean:
        ztot = 0.
        do j=1,nlatgo
          jp = min (j+1,nlatgo)
          jm = max (j-1,1)
          zda = 2.*pi*(radius**2) * cos(alatgo(j))
     *          * 0.5*(alatgo(jp)-alatgo(jm))
          if (j.eq.1 .or. j.eq.nlatgo) zda = 2.*zda ! no grid pt at pole
          do i=1,nlongo
            darlco(i,j) = zda/nlongo
            ztot = ztot + darlco(i,j)
          enddo
        enddo

c       adjust global total:
        zadj = 4.*pi*(radius**2) / ztot
        do j=1,nlatgo
          do i=1,nlongo
            darlco(i,j) = darlco(i,j) * zadj
          enddo
        enddo

        firstlc = .false.
c>>>>>>>>>>
      endif
c>>>>>>>>>>

      return
      end

c-----------------------------------------------------------------------

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
c         filc  = fraction of LC atmos cell occupied by IM model cells
c         zilc  = mean sfc elevation of IM model (within filc) (m)
c         fili  = fraction of LC atmos cell with IM ice (m)
c         zili  = mean sfc elevation of IM ice (within fili) (m)

c         ficco  = fraction of LC ocean cell occupied by IM model cells
c         filco  = fraction of LC ocean cell occupied by IM ocean
c                  (with maskwater=1)
c         hilco  = mean thickness of IM ice (floating, within filco) (m)
c         hblco  = mean elevation of IM bed (within filco) (m)
c         rullco = liquid runoff from i.s.,av for LC ocn cell(m.w.e./yr)
c         ruflco = ice discharge from i.s.,av for LC ocn cell(m.w.e./yr)

c        For discharge from ice sheet ("runoff" runliqav, runfrozav):
c          (1) liq is routed to ice-model ocean or ice domain edge. 
c              on ice grid. froz is too, but should aready be in
c              or adjacent to IM ocean.
c          (2) liq and froz are tiled to LC grid.
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
     *  inearlco(nx,ny),     jnearlco(nx,ny)

c     local, LC atmosphere grid: 
      dimension
     *  ailc(nlong,nlatg),   filc(nlong,nlatg),  zilc(nlong,nlatg), 
     *  aili(nlong,nlatg),   fili(nlong,nlatg),  zili(nlong,nlatg), 
     *  worklc(nlong,nlatg)

c     local, LC ocean grid: 
      dimension
     *  aicco(nlongo,nlatgo),  ficco(nlongo,nlatgo),
     *  ailco(nlongo,nlatgo),  filco(nlongo,nlatgo),
     *  hilco(nlongo,nlatgo),  hblco(nlongo,nlatgo),
     *  rullco(nlongo,nlatgo), ruflco(nlongo,nlatgo),
     *  worklco(nlongo,nlatgo)

c     for printmaplc output only:
      dimension masklc(nlong,nlatg), masklco(nlongo,nlatgo) 

c++++++++++++++++++++++++
c     for Netcdf writing:
c++++++++++++++++++++++++
#include <netcdf.inc>
      integer*4
     * lenattr, nx4, ny4, nlev4, ival4, jval4,
     * varid, rcode, ncerr, ncid,
     * xida, yida, xido, yido, levid, timid, minmaxid,
     * minmax(2), ndim, idim(4), start(4), count(4),
     * n1, n2, n3, n4, n5, n6, n7, n8, n9,n10,
     * n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,
     * n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls.
     * ncglobal4, ncunlim4,
     * nf_float4, nf_short4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_short4 = nf_short,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)
      save varid, rcode, ncid,
     *     xida, yida, xido, yido, levid, timid, minmaxid,
     *     minmax, idim, start, count,
     *     n1, n2, n3, n4, n5, n6, n7, n8, n9,n10,
     *     n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,
     *     n21,n22,n23,n24,n25,n26,n27,n28,n29,n30
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
     *               runliqav, runliqavo, runfrozav, runfrozavo,  
     *               darea, sealev, 0., .false., nx, ny) 

c        Calculate LC box indices containing IM grid centers

c     atmos:
      call setintnear (along, alatg, nlong, nlatg, 
     *                 alond, alatd, nx, ny,
     *                 inearlc, jnearlc)

c     ocean:
      call setintnear (alongo, alatgo, nlongo, nlatgo, 
     *                 alond, alatd, nx, ny,
     *                 inearlco, jnearlco)

c        Simple tiling to LC grid, ignoring overlaps of ice boxes 
c        across LC borders. 

      call zero (ailc, nlong*nlatg)
      call zero (filc, nlong*nlatg)
      call zero (zilc, nlong*nlatg)
      call zero (aili, nlong*nlatg)
      call zero (fili, nlong*nlatg)
      call zero (zili, nlong*nlatg)

      call zero (aicco, nlongo*nlatgo)
      call zero (ficco, nlongo*nlatgo)
      call zero (ailco, nlongo*nlatgo)
      call zero (filco, nlongo*nlatgo)
      call zero (hilco, nlongo*nlatgo)
      call zero (hblco, nlongo*nlatgo)
      call zero (rullco, nlongo*nlatgo)
      call zero (ruflco, nlongo*nlatgo)

c     ztott = max (.001, timeice-nyearstart)
      ztott = max (.001, timeice-timefluxprev)
      do j=1,ny 
        do i=1,nx 
          ic = inearlc(i,j) 
          jc = jnearlc(i,j) 
          zda = darea(i,j)
          ailc(ic,jc) = ailc(ic,jc) + zda
          zilc(ic,jc) = zilc(ic,jc) + (hs(i,j)-sealev)*zda
          if (h(i,j).gt.0.) then
            aili(ic,jc) = aili(ic,jc) + zda
            zili(ic,jc) = zili(ic,jc) + (hs(i,j)-sealev)*zda
          endif

          ico = inearlco(i,j) 
          jco = jnearlco(i,j) 
          zda = darea(i,j)
          aicco(ico,jco) = aicco(ico,jco) + zda
          if (maskwater(i,j).eq.1) then
            ailco(ico,jco) = ailco(ico,jco) + zda
            hilco(ico,jco) = hilco(ico,jco) + h(i,j)*zda
            hblco(ico,jco) = hblco(ico,jco) + (hb(i,j)-sealev)*zda
          endif
          rullco(ico,jco) = rullco(ico,jco) 
     *                      + rhor*(runliqavo(i,j)/ztott)*zda
          ruflco(ico,jco) = ruflco(ico,jco) 
     *                      + rhor*(runfrozavo(i,j)/ztott)*zda
        enddo
      enddo

c        Set final LC arrays 

      do jc=1,nlatg 
        do ic=1,nlong
          if (ailc(ic,jc).gt.1.e-6) then
            filc(ic,jc) = min (1., ailc(ic,jc)/darlc(ic,jc))
            zilc(ic,jc) = zilc(ic,jc) / ailc(ic,jc)
          endif
          if (aili(ic,jc).gt.1.e-6) then
            fili(ic,jc) = min (1., aili(ic,jc)/darlc(ic,jc))
            zili(ic,jc) = zili(ic,jc) / aili(ic,jc)
          endif
        enddo
      enddo

      do jco=1,nlatgo
        do ico=1,nlongo
          if (aicco(ico,jco).gt.1.e-6) then
            ficco(ico,jco) = min (1., aicco(ico,jco)/darlco(ico,jco))
          endif
          if (ailco(ico,jco).gt.1.e-6) then
            filco(ico,jco) = min (1., ailco(ico,jco)/darlco(ico,jco))
            hilco(ico,jco) = hilco(ico,jco) / ailco(ico,jco)
            hblco(ico,jco) = hblco(ico,jco) / ailco(ico,jco)
          endif
          rullco(ico,jco) = rullco(ico,jco) / darlco(ico,jco)
          ruflco(ico,jco) = ruflco(ico,jco) / darlco(ico,jco)
        enddo
      enddo

c        Do simple multiplicative global correction for runoffs
c        to ensure global total on final Loveclim grid is equal 
c        to total on ice-model grid

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
      do jco=1,nlatgo
        do ico=1,nlongo
          totrulo_lc = totrulo_lc + rullco(ico,jco)*darlco(ico,jco)
          totrufo_lc = totrufo_lc + ruflco(ico,jco)*darlco(ico,jco)
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
        do jco=1,nlatgo
          do ico=1,nlongo
            rullco(ico,jco) = rullco(ico,jco) * ratioliq
          enddo
        enddo
      endif
      if (totrufo_lc.gt.1.e-10) then
        do jco=1,nlatgo
          do ico=1,nlongo
            ruflco(ico,jco) = ruflco(ico,jco) * ratiofroz
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

c        Define and write x dimension (longitude, atmosphere)

      nx4 = nlong
      rcode = nf_def_dim (ncid, 'lon', nx4, xida)
      rcode = nf_def_var (ncid, 'lon', nf_float4, n1, xida, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                         n16, 'longitude, atmos')
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

c        Define and write y dimension (latitude, atmosphere)

      ny4 = nlatg
      rcode = nf_def_dim (ncid, 'lat', ny4, yida)
      rcode = nf_def_var (ncid, 'lat', nf_float4, n1, yida, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                         n15, 'latitude, atmos')
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

c        Define and write x dimension (longitude, ocean)

      nx4 = nlongo
      rcode = nf_def_dim (ncid, 'lono', nx4, xido)
      rcode = nf_def_var (ncid, 'lono', nf_float4, n1, xido, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                         n16, 'longitude, ocean')
      rcode = nf_put_att_text (ncid, varid, 'units',
     *                         n12, 'degrees_east')
      rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                         n4, 'f8.3')
      rcode = nf_enddef (ncid)
      do i=1,nlongo
        ztmp(i) = alongo(i)*180./pi
      enddo
      rcode = nf_put_vara_double (ncid, varid, n1, nx4, ztmp)
      rcode = nf_redef (ncid)

c        Define and write y dimension (latitude, ocean)

      ny4 = nlatgo
      rcode = nf_def_dim (ncid, 'lato', ny4, yido)
      rcode = nf_def_var (ncid, 'lato', nf_float4, n1, yido, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                         n15, 'latitude, ocean')
      rcode = nf_put_att_text (ncid, varid, 'units',
     *                         n13, 'degrees_north')
      rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                         n4, 'f8.3')
      rcode = nf_enddef (ncid)
      do j=1,nlatgo
        ztmp(j) = alatgo(j)*180./pi
      enddo
      rcode = nf_put_vara_double (ncid, varid, n1, ny4, ztmp)
      rcode = nf_redef (ncid)

c        Define variable dimensions and attributes

c---------------
c        For 2D:
c---------------

      ndim    = 2

c     atmos:

      idim(1) = xida
      idim(2) = yida

      cvar = 'fice_im'
      cvarl = 'fractional area of ice model'
      cunits = '0 to 1'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'zice_im'
      cvarl = 'surface elevation of ice model'
      cunits = 'm'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'fici_im'
      cvarl = 'fractional area of ice-model ice'
      cunits = '0 to 1'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'zici_im'
      cvarl = 'surface elevation of ice-model ice'
      cunits = 'm'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c     ocean:

      idim(1) = xido
      idim(2) = yido

      cvar = 'focc_im'
      cvarl = 'fractional area of ice model'
      cunits = '0 to 1'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'focn_im'
      cvarl = 'fractional area of ice-model ocean'
      cunits = '0 to 1'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'hiceo_im'
      cvarl = 'thickness of ice-model floating ice'
      cunits = 'm'
      ival4 = lenchr(cvarl)
      jval4 = lenchr(cunits)
      rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
      rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
      rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
      rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

      cvar = 'hbedo_im'
      cvarl = 'elevation of ice-model ocean bed'
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

c        Leave define mode

      rcode = nf_enddef (ncid)

c        Write variables

c~~~~~~~~~~~~~~~
c        For 2D:
c~~~~~~~~~~~~~~~

c     atmos:

      start(1) = 1
      count(1) = nlong
      start(2) = 1
      count(2) = nlatg

      rcode = nf_inq_varid (ncid, 'fice_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, filc)

      rcode = nf_inq_varid (ncid, 'zice_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, zilc)

      rcode = nf_inq_varid (ncid, 'fici_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, fili)

      rcode = nf_inq_varid (ncid, 'zici_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, zili)

c     ocean:

      start(1) = 1
      count(1) = nlongo
      start(2) = 1
      count(2) = nlatgo

      rcode = nf_inq_varid (ncid, 'focc_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, ficco)

      rcode = nf_inq_varid (ncid, 'focn_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, filco)

      rcode = nf_inq_varid (ncid, 'hiceo_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hilco)

      rcode = nf_inq_varid (ncid, 'hbedo_im', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hblco)

      rcode = nf_inq_varid (ncid, 'runoff_liq', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, rullco)

      rcode = nf_inq_varid (ncid, 'runoff_froz', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, ruflco)

c        Return to define mode (not necessary?) and close

      rcode = nf_redef (ncid)
      rcode = nf_close (ncid)

c        Ascii output (diagnostic only)

      iu = 154

c     atmos:

c     land-ocean mask for reverse-video outlines:
      do j=1,nlatg
        do i=1,nlong
          if (filc(i,j).gt.0.5) then
             masklc(i,j) = 1
          else    
             masklc(i,j) = 0
          endif   
          worklc(i,j) = masklc(i,j)
        enddo   
      enddo
      call printmaplc (timeice, worklc, 'masklc', 0.1, 
     *                 masklc, along, alatg, nlong, nlatg, iu)
      call printmaplc (timeice, filc, 'filc', 0.1, 
     *                 masklc, along, alatg, nlong, nlatg, iu)
      call printmaplc (timeice, zilc, 'zilc', 150., 
     *                 masklc, along, alatg, nlong, nlatg, iu)
      call printmaplc (timeice, fili, 'fili', 0.1, 
     *                 masklc, along, alatg, nlong, nlatg, iu)
      call printmaplc (timeice, zili, 'zili', 150., 
     *                 masklc, along, alatg, nlong, nlatg, iu)

c     ocean:

      do j=1,nlatgo
        do i=1,nlongo
          if (ficco(i,j).gt.0.5) then
             masklco(i,j) = 1
          else    
             masklco(i,j) = 0
          endif   
          worklco(i,j) = masklco(i,j)
        enddo   
      enddo
      call printmaplc (timeice, worklco, 'masklco', 0.1,
     *                 masklco, alongo, alatgo, nlongo, nlatgo, iu)
      call printmaplc (timeice, ficco,  'ficco', 0.1,
     *                 masklco, alongo, alatgo, nlongo, nlatgo, iu)
      call printmaplc (timeice, filco,  'filco', 0.1,
     *                 masklco, alongo, alatgo, nlongo, nlatgo, iu)
      call printmaplc (timeice, hilco,  'hilco', 50.,
     *                 masklco, alongo, alatgo, nlongo, nlatgo, iu)
      call printmaplc (timeice, hblco,  'hblco', 150.,
     *                 masklco, alongo, alatgo, nlongo, nlatgo, iu)
      call printmaplc (timeice, rullco, 'rullco', 0.1,
     *                 masklco, alongo, alatgo, nlongo, nlatgo, iu)
      call printmaplc (timeice, ruflco, 'ruflco', 0.1,
     *                 masklco, alongo, alatgo, nlongo, nlatgo, iu)

c     ice:

      call printmap (timeice, runliqav,   'runliqav',   0.1, iu, 0)
      call printmap (timeice, runliqavo,  'runliqavo',  0.1, iu, 0)
      call printmap (timeice, runfrozav,  'runfrozav',  0.1, iu, 0)
      call printmap (timeice, runfrozavo, 'runfrozavo', 0.1, iu, 0)

      call flush(iu)

      return
      end

c-----------------------------------------------------------------------

      subroutine routerun (hypot, indroute,
     *                     runliqav, runliqavo, runfrozav, runfrozavo,  
     *                     darea, sealev, vnull, ifwrap, nx, ny) 

c        Routes runoff (liquid) to IM ocean, i.e, sets runliqavo given 
c        runliqav. Also returns indroute = i + nx*(j-1), index of
c        routed ocean point for each liquid source point. 
c        Does same for frozen discharge runfrozav -> runfrozavo,
c        (runfrozav should already be in or adjacent to IM ocean).

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
     *  runliqav(nx,ny),  runliqavo(nx,ny),
     *  runfrozav(nx,ny), runfrozavo(nx,ny),
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

      subroutine printmaplc (time, arrin, cdescin, scale, 
     *                       mask, alonin, alatin, nlonin, nlatin, iu)

c        Displays array arr to unit iu in "a1" format.
c        Like printmap (iceshow.F), but for Loveclim (LC) variables.
c        mask = 0 for land, 1 for ocean (for reverse video oocean bdry).

#include <comicegrid.h>

      dimension arrin(nlonin,nlatin), mask(nlonin,nlatin),
     *          alonin(nlonin), alatin(nlatin)
      character*(*) cdescin

      dimension arr(nlonin,nlatin)
      character*1 carr(nlonin)
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

      call scopy (nlonin*nlatin, arrin, 1, arr, 1)

      ia = 1
      ib = nlonin
      ic = 1

      if (.not. first) write (iu,*) 
      write(iu,100)
     *  time, cdesc(1:lenchr(cdesc)), scale
      write(iu,101) 
     *  (nint(alonin(i)*180./pi), i=ia+4,ib,5*ic)
  100 format(
     *  'time:',f12.1,'  ',a,'/',f11.4,'  (A-Z pos, a-z neg):')
  101 format(/6x,2x,500i5)

c>>>>>>>>>>>>>>>>>>>>>>>>>
      do j = nlatin, 1, -ic
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
c           ip = min (i+ic,nlonin)                ! walls
            im = mod (i-ic+nlonin -1 ,nlonin) + 1  ! wraparound
            ip = mod (i+ic+nlonin -1 ,nlonin) + 1  ! wraparound
            jm = max (j-ic,1)
            jp = min (j+ic,nlatin)
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

        write(iu,'(f6.1,2x,a)') alatin(j)*180./pi, cbrr(1:lenb)

c>>>>>>>>>>
      enddo
c>>>>>>>>>>

      first = .false.

      return
      end

c-----------------------------------------------------------------------

c@@@@@
#endif
c@@@@@
