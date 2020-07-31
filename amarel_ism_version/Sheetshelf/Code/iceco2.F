c*********************
#if defined (CO2INTER)
c*********************

c-----------------------------------------------------------------------

      subroutine doco2 (tmon, pmon, runofa, h, maskwater, 
     *                  toge, tage, prge, ruge,
     *                  rco2inter, ecc, obl, prec, 
     *                  timeice, dtimeco2, nyearstart)

c     Calls dosoil twice to calculate (i) Antarctic CO2 weathering (ice
c     grid), and (ii) global non-Antarctic CO2 weathering (GCM grid),
c     and to step atmos co2.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

      dimension 
     *  tmon(nmon,nx,ny),       pmon(nmon,nx,ny),    runofa(nx,ny),
     *  h(nx,ny),               maskwater(nx,ny),
     *  toge(nlong,nlatg),     
     *  tage(nmon,nlong,nlatg), prge(nmon,nlong,nlatg), 
     *  ruge(nlong,nlatg)

      character*200 cfile

      dimension 
     *  masktmp(nx,ny), 
     *  maskg(nlong,nlatg),  workg(nlong,nlatg),
     *  dareag(nlong,nlatg), alongd(nlong,nlatg), alatgd(nlong,nlatg)
      save maskg, dareag, alongd, alatgd, firstg
      logical firstg
      data firstg /.true./

      dimension alatt31(48), wgaust31(48)
      data alatt31 /
     *  -1.5212131731, -1.4569823009, -1.3923721142, -1.3276762287,
     *  -1.2629478409, -1.1982037311, -1.1334508412, -1.0686925583,
     *  -1.0039307308, -0.9391664533, -0.8744004155, -0.8096330743,
     *  -0.7448647446, -0.6800956510, -0.6153259588, -0.5505557923,
     *  -0.4857852479, -0.4210144016, -0.3562433151, -0.2914720392,
     *  -0.2267006173, -0.1619290868, -0.0971574813, -0.0323858319,
     *   0.0323858319,  0.0971574813,  0.1619290868,  0.2267006173,
     *   0.2914720392,  0.3562433151,  0.4210144016,  0.4857852479,
     *   0.5505557923,  0.6153259588,  0.6800956510,  0.7448647446,
     *   0.8096330743,  0.8744004155,  0.9391664533,  1.0039307308,
     *   1.0686925583,  1.1334508412,  1.1982037311,  1.2629478409,
     *   1.3276762287,  1.3923721142,  1.4569823009,  1.5212131731
     *   /

      data wgaust31 /
     *   0.0031533461,  0.0073275539,  0.0114772346,  0.0155793157,
     *   0.0196161605,  0.0235707608,  0.0274265097,  0.0311672278,
     *   0.0347772226,  0.0382413511,  0.0415450829,  0.0446745609,
     *   0.0476166585,  0.0503590356,  0.0528901895,  0.0551995037,
     *   0.0572772921,  0.0591148397,  0.0607044392,  0.0620394232,
     *   0.0631141923,  0.0639242386,  0.0644661644,  0.0647376968,
     *   0.0647376968,  0.0644661644,  0.0639242386,  0.0631141923,
     *   0.0620394232,  0.0607044392,  0.0591148397,  0.0572772921,
     *   0.0551995037,  0.0528901895,  0.0503590356,  0.0476166585,
     *   0.0446745609,  0.0415450829,  0.0382413511,  0.0347772226,
     *   0.0311672278,  0.0274265097,  0.0235707608,  0.0196161605,
     *   0.0155793157,  0.0114772346,  0.0073275539,  0.0031533461
     *   /

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c        Run soil bucket model + weathering, for Antarctic (ice grid) 
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.) then
             masktmp(i,j) = 2
          else if (maskwater(i,j).ge.1) then 
             masktmp(i,j) = 1
          else
             masktmp(i,j) = 0
          endif
        enddo
      enddo

      call dosoil (tmon, pmon, runofa, masktmp,
     *             darea, radius, alond, alatd, 
     *             rco2inter, ecc, obl, prec, 
     *             timeice, dtimeco2, nyearstart, nx, ny, nmon,
     *             0)

c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c        Run soil bucket model + weathering, for globe except Antarctic
c        (gcm grid) 
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c     if first call, set gcm grid axes, etc:
c~~~~~~~~~~~~~~~~~~~~~
      if (firstg) then
c~~~~~~~~~~~~~~~~~~~~~
        do j=1,nlatg
          do i=1,nlong
            alongd(i,j) = (-pi + 2.*pi*(i-.5)/nlong) * 180./pi
            if (nlatg.eq.48) then
              alatgd(i,j) = alatt31(j) * 180./pi
              dareag(i,j) = wgaust31(j)
            else
              zlat        = -0.5*pi + pi*(j-0.5)/nlatg
              alatgd(i,j) = zlat * 180./pi
              dareag(i,j) = cos(zlat)
            endif
          enddo
        enddo

        totglob = 4.*pi*radius**2
        totcheck = 0.
        do j=1,nlatg
          do i=1,nlong
            totcheck =  totcheck + dareag(i,j)
          enddo
        enddo
        totcheck2 = 0.
        do j=1,nlatg
          do i=1,nlong
            dareag(i,j) = dareag(i,j)*totglob/totcheck
            totcheck2 =  totcheck2 + dareag(i,j)
          enddo
        enddo

c       Get lsx land-ocean mask from Genesis data input file 
c       (from no-ice file...)
        cfile = cpathmap(1:lenchr(cpathmap))
     *          // '34ma/data_35ma_non_sur_2X2'
        open (iuclim, file=cfile, form='formatted', status='old')
        call readdata_ice (iuclim, workg, nlong,nlatg,
     *                     nlong,nlatg, 'SURT    ',  1, 1)
        close (iuclim)
        do j=1,nlatg
          do i=1,nlong
            if (nint(workg(i,j)).eq.0) then
              maskg(i,j) =  1     ! ocean
            else if (nint(workg(i,j)).eq.1) then
              maskg(i,j) =  0     ! land
            else if (nint(workg(i,j)).eq.2) then
              maskg(i,j) = 2      ! ice sheet
            endif
          enddo
        enddo

c       diagnostic printouts:
c       do j=nlatg,1,-1
c         write (175,'(180i1)') (maskg(i,j),i=1,nlong)
c       enddo
c       write (175,'(i4, f8.2)') (j, alatgd(nlong/2,j), j=1,nlatg)
c       write (175,'(f12.6)') totcheck2/totglob 
c       do j=nlatg,1,-1
c       enddo
c       do m=1,nmon
c         write (176,'(a,i4)') 'mon=',m
c         write (177,'(a,i4)') 'mon=',m
c         if (m.eq.1) write (179,'(a,i4)') 'mon=',m
c         do j=nlatg,1,-1
c           write (176,'(180i4)') (nint(tage(m,i,j)),i=1,nlong)
c           write (177,'(180i4)') (nint(10.*prge(m,i,j)),i=1,nlong)
c           if (m.eq.1) then
c             write (178,'(180i4)') (nint(10.*ruge(i,j)),i=1,nlong)
c             write (179,'(180i4)') (nint(0.1*toge(i,j)),i=1,nlong)
c           endif
c         enddo
c       enddo
 
        firstg = .false.
c~~~~~~~~~~
      endif
c~~~~~~~~~~

      call dosoil (tage, prge, ruge, maskg, 
     *             dareag, radius, alongd, alatgd, 
     *             rco2inter, ecc, obl, prec, 
     *             timeice, dtimeco2, nyearstart, nlong, nlatg, nmon, 
     *             1)

      return
      end

c-----------------------------------------------------------------------

      subroutine dosoil (tsurf, precip, runof, mask,
     *                   darea, radius, alond, alatd,
     *                   rco2, ecc, obliq, prec,
     *                   timeice, dtimeco2, nyearstart, nx, ny, nmon,
     *                   iflag)

c        Runs soil bucket model for 2-D grid (ice model or global-Ant),
c        calling bucket and siweather for each point, and
c        accumulating total Si weathering in weathering_tot.
c        Called twice, 1st time with iflag=0 (Ant ice grid), then
c        with iflag=1 (global minus Ant). Saves results (eg, runoff
c        and weathering) for each grid in local vars, and if iflag=1,
c        steps co2 due to combined weathering.

c        tsurf = monthly surface air temperature, deg K (supp)
c        precip = monthly precip rate, m/yr (supp)
c                 nb: units in dosoil/bucket are deg C, mm/day. 
c        runof= annual runoff rate, mm/day (supp) (used only iflag=1) 
c        mask = 0,1,2 for land, ocean, ice (supp)
c        darea = area of each grid box (m2) (supp)
c        alond, alatd = lon, lat (deg) (supp)
c        rco2 = co2 level (x PAL) (mod)
c        ecc, obliq, prec = orbit (supp)
c        timeice = current time (yr BP) (supp)
c        dtimeco2 = time step each call (yr) (supp)
c        nx,ny = domain dimension (ice grid, or global) (supp)
c        iflag = 0 for ice grid, 1 for global (supp)
c        weathering_tot = total weathering,this domain (mol CO2/yr)(ret)

      dimension
     *  tsurf(nmon,nx,ny),  precip(nmon,nx,ny),  
     *  runof(nx,ny),       mask(nx,ny),
     *  darea(nx,ny),       alond(nx,ny) ,       alatd(nx,ny)

c     avoid dynamic array allocation:
      parameter (nxmax=200, nymax=200, nmonmax=12)
      dimension
     *  precipa(nxmax,nymax),    evapa(nxmax,nymax),
     *  runofa(nxmax,nymax),     dstora(nxmax,nymax),
     *  ta(nmonmax),             pr(nmonmax),       ho(nmonmax),
     *  weathering(nxmax,nymax),
     *  tempa(nxmax,nymax)      ! diagnostic only

      dimension ndaypm(12)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm

      logical firsta, firstsiw
      data firsta, firstsiw /.true., .true./
      save firsta, firstsiw

      save 
     *  area_ant,       area_glo,       area_all, 
     *  temp_ant,       temp_glo,       temp_all, 
     *  precip_ant,     precip_glo,     precip_all, 
     *  evap_ant,       evap_glo,       evap_all, 
     *  runof_ant,      runof_glo,      runof_all, 
     *  dstor_ant,      dstor_glo,      dstor_all, 
     *  weathering_ant, weathering_glo, weathering_all

      parameter (tmelt = 273.15)
      parameter (pi=3.14159265358979)

c     ecc =  0.016706             ! eccentricity
c     obliq = 23.4377 * pi/180.   ! obliquity
c     prec =  77.0613 * pi/180.   ! precession, from peri. to NH V.E.
      zs = sin (0.5*prec) / sqrt ((1.+ecc)/(1.-ecc))
      zc = cos (0.5*prec)
      ze = 2. * atan2 (zs,zc)
      vern = ze - ecc * sin(ze)

c        Set silicate weatherability factor

c     siweatherability = min(1., 0.8 + 0.2*(timeice-nyearstart)/150000.)

c     siweatherability = 0.8  ! co2a
c     siweatherability = 0.85 ! co2b
c     siweatherability = 0.9  ! co2c
c     siweatherability = 0.95 ! co2d
c     siweatherability = 1.   ! co2e

cb4   siweatherability = 1.58 ! 888888
c1    siweatherability = 1.8  ! 888888
c2    siweatherability = 2.0  ! 888888
c3    siweatherability = 1.5  ! 888888
c4    siweatherability = 1.4  ! 888888
c     siweatherability = 1.25 ! 888888 ! no ice
c     siweatherability = 1.7  ! 888888 ! no ice

      siweatherability = 1.    ! d0,d1
c     siweatherability = 0.8   ! d2
c     siweatherability = 1.2   ! d3

      area_tot   = 0.
      temp_tot   = 0.
      precip_tot = 0.
      evap_tot   = 0.
      runof_tot  = 0.
      dstor_tot  = 0.
      weathering_tot = 0.

c         Do soil bucket model (annual cycle) for every land pt

c================
      do j=1,ny
        do i=1,nx
c================

c--------------------------------------------
c         for land only (not ice, not ocean):
          if ( mask(i,j).eq.0) then
c--------------------------------------------

c...............................
            if (iflag.eq.0) then
c...............................

c                Calculate for Antarc (ice grid)

              do m=1,nmon
                pr(m) = precip(m,i,j)*(1000./365.)      ! m/yr to mm/day
                ta(m) = tsurf(m,i,j) - tmelt            ! deg K to C
c               calculate number of daylight hours at mid-month days:
                iyear = 0     
                if (m.eq.1) then 
                  zday  = 0.5*ndaypm(1)
                else
                  zday = zday + 0.5*(ndaypm(m-1)+ndaypm(m))
                endif
                isecyr = 86400*nint(zday)
                isecdy = 86400/2
                dt = 86400.
                call zencal (iyear, isecyr, isecdy, dt,
     *                       ecc, obliq, prec, vern, dist, eccf, 
     *                       cosq, fraq, cosq24, fraq24, 
     *                       alatd(i,j)*pi/180., 1)
        
                ho(m) = 24.*fraq
              enddo

c                Calculate runoff

              call bucket (pr, ta, ho, 
     *                     alond(i,j), alatd(i,j), nmon, 
     *                     precipa(i,j), evapa(i,j), 
     *                     runofa(i,j),  dstora(i,j), iflag)

c....................................
            else if (iflag.eq.1) then
c....................................

c                Calculate for globe except Antarc:

              do m=1,nmon
                ta(m) = tsurf(m,i,j) - tmelt                ! deg K to C
              enddo

              runofa(i,j) = runof(i,j)

c             diagnostic only:
              precipa(i,j) = 0.
              do m=1,nmon
                precipa(i,j) = precipa(i,j) 
     *                       + precip(m,i,j)*(1000./365.) ! m/y to mm/dy
     *                         *ndaypm(m)/365.
              enddo
              evapa(i,j) = 0.
c................
            endif
c................

c           diagnostic only:
            tempa(i,j) = 0.
            do m=1,nmon
              tempa(i,j) = tempa(i,j) 
     *                   + (tsurf(m,i,j)-tmelt)*ndaypm(m)/365. !  K to C
            enddo

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            zlan  = -60.
            if ( (iflag.eq.0.and.alatd(i,j).le.zlan) .or.    ! ice model
     *           (iflag.eq.1.and.alatd(i,j).gt.zlan)         ! gcm 
cz   *           (iflag.eq.1.and.alatd(i,j).le.zlan)         ! test only
     *         ) then
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c                Calculate silicate weathering

              call siweather1 (ta, nmon, runofa(i,j), weathering(i,j),
     *                         siweatherability)
c             call siweather2 (ta, nmon, runofa(i,j), weathering(i,j),
c    *                         siweatherability)


              area_tot = area_tot      +              darea(i,j)
c             (mm/day) * m^2:
              temp_tot   = temp_tot    + tempa(i,j)  *darea(i,j)
              precip_tot = precip_tot  + precipa(i,j)*darea(i,j)
              evap_tot   = evap_tot    + evapa(i,j)  *darea(i,j)
              runof_tot  = runof_tot   + runofa(i,j) *darea(i,j)
              dstor_tot  = dstor_tot   + dstora(i,j) *darea(i,j)
c             moles of CO2/yr:
              weathering_tot=weathering_tot + weathering(i,j)*darea(i,j)

c++++++++++++++++
            endif
c++++++++++++++++

c-------------
          else
c-------------

            tempa(i,j) = 0.
            precipa(i,j) = 0.
            evapa(i,j) = 0.
            runofa(i,j) = 0.
            dstora(i,j) = 0.
            weathering(i,j) = 0.

c--------------
          endif
c--------------
c============
        enddo
      enddo
c============

      if (iflag.eq.0) then
         area_ant       = area_tot
         temp_ant       = temp_tot
         precip_ant     = precip_tot
         evap_ant       = evap_tot
         runof_ant      = runof_tot
         dstor_ant      = dstor_tot
         weathering_ant = weathering_tot 
       else
         area_glo       = area_tot
         temp_glo       = temp_tot
         precip_glo     = precip_tot
         evap_glo       = evap_tot
         runof_glo      = runof_tot
         dstor_glo      = dstor_tot
         weathering_glo = weathering_tot 

         area_all       = area_ant       + area_glo
         temp_all       = temp_ant       + temp_glo
         precip_all     = precip_ant     + precip_glo
         evap_all       = evap_ant       + evap_glo
         runof_all      = runof_ant      + runof_glo
         dstor_all      = dstor_ant      + dstor_glo
         weathering_all =  weathering_ant + weathering_glo         !8888
co2o_lastweathering_all = 3.*weathering_ant + 0.9*weathering_glo   !8888
       endif

      if (iflag.eq.0) then
        if (firsta) then
c         do j=ny,1,-1
c           write (190,'(200i5)') (nint(10.*tempa(i,j)),i=1,nx)
c           write (191,'(200i5)') (nint(10.*precipa(i,j)),i=1,nx)
c           write (192,'(200i5)') (nint(10.*evapa(i,j)),i=1,nx)
c           write (193,'(200i5)') (nint(10.*runofa(i,j)),i=1,nx)
c           write (194,'(200i5)') (nint(10.*dstora(i,j)),i=1,nx)
c           write (195,'(200i5)') (mask(i,j),i=1,nx)
c           write (196,'(200i5)') (nint(10.*runof(i,j)),i=1,nx)
c         enddo
          firsta = .false.
        endif
      endif

c         Step to new atmospheric co2 level

c)))))))))))))))))))))))))
      if (iflag.eq.1) then
c)))))))))))))))))))))))))
        if (firstsiw) then
c         do j=ny,1,-1
c           write (180,'(200i5)') (nint(10.*tempa(i,j)),i=1,nx)
c           write (181,'(200i5)') (nint(10.*precipa(i,j)),i=1,nx)
c           write (182,'(200i5)') (nint(10.*evapa(i,j)),i=1,nx)
c           write (183,'(200i5)') (nint(10.*runofa(i,j)),i=1,nx)
c           write (184,'(200i5)') (nint(10.*dstora(i,j)),i=1,nx)
c           write (185,'(200i5)') (mask(i,j),i=1,nx)
c         enddo

          write (82,'(9a)') '        time',
     *   '  area_ant  area_glo  area_all',
     *   '  temp_ant  temp_glo  temp_all',
     *   '  prec_ant  prec_glo  prec_all',
     *   '  evap_ant  evap_glo  evap_all',
     *   '  runo_ant  runo_glo  runo_all',
     *   ' dstor_all   err_all',
     *   ' weath_ant weath_glo weath_all   fluxorg      volc      fnet',
     *   '      rco2   siweath'
          firstsiw = .false.
        endif

c       zarea = 0.3 * 4.*pi*(radius**2)       ! ~all land area

c          Get current organic carbon sink rate (mol CO2/yr)

        fluxorg = 0.
c       call orgsink (fluxorg, timeice, iuorg)       ! 8888
c       fluxorg = fluxorg - 0.6e12                   ! 8888
cb1     fluxorg = fluxorg - 0.9e12                   ! 8888
c       fluxorg = fluxorg - 1.2e12                   ! 8888

c       volc    = 6.4e12                             ! 8888
c       volc    = 6.0e12                             ! 8888
        volc    = 4.35080e12                         ! 8888  ! co2=2.003
c       volc    = 3.5e12                             ! 8888  ! co2=
        fnet    = volc - weathering_all - fluxorg

        write (82,'(i12,
     *              3i10,
     *              3f10.3,
     *              11f10.6,
     *              8f10.5)')
c    *               f10.5,7e15.6)')
     *    nint(timeice),

     *    nint(area_ant*1.e-6),                          ! km2
     *    nint(area_glo*1.e-6),                          ! km2
     *    nint(area_all*1.e-6),                          ! km2

     *    temp_ant/max(area_ant,1.),                     ! deg C
     *    temp_glo/area_glo,                             ! deg C
     *    temp_all/area_all,                             ! deg C

     *    precip_ant/max(area_ant,1.),                   ! mm/day
     *    precip_glo/area_glo,                           ! mm/day
     *    precip_all/area_all,                           ! mm/day
     *    evap_ant/max(area_ant,1.),                     ! mm/day
     *    evap_glo/area_glo,                             ! mm/day
     *    evap_all/area_all,                             ! mm/day
     *    runof_ant/max(area_ant,1.),                    ! mm/day
     *    runof_glo/area_glo,                            ! mm/day
     *    runof_all/area_all,                            ! mm/day
     *    dstor_all/area_all,
     *    (precip_all - evap_all - runof_all - dstor_all) / area_all,

     *   -weathering_ant/1.e12,                     ! (mol CO2/yr)/1.e12
     *   -weathering_glo/1.e12,
     *   -weathering_all/1.e12,
     *   -fluxorg/1.e12,
     *    volc/1.e12,
     *    fnet/1.e12,
     *    rco2,
     *    siweatherability
        call flush (82)

        sco2 = sqrt(rco2)     ! lk 1/16/09, accounting for ocean storage

        sco2 = sco2 + fnet * dtimeco2 
co2o *                     / 2.e16   ! 8888
co2ob*                     / 2.e17   ! 8888
co2oa*                     / 6.e17   ! 8888
co2q *                     / 3.8e18  ! 8888
cb1  *                     / 1.0e18  ! 8888
cb2  *                     / 3.0e18  ! 8888
cb3  *                     / 2.e17   ! 8888
     *                     / 0.6e17  ! 8888

        sco2 = max (sco2, 0.)
        rco2 = sco2**2        ! lk 1/16/09

        rco2 = max (0.02, 4. - 0.1*(timeice +34.e6)/200.) ! 888888
 
c))))))))))
      endif
c))))))))))

      return
      end

c----------------------------------------------------------------------

      subroutine bucket (pr, ta, ho, 
     *                   alond, alatd, nmon, 
     *                   precip, evap, runof, dstor, iflag)

c        Solves soil bucket model for one point. Monthly precip and
c        temperature are supplied to drive the model through nyear
c        years to ~equilibrium, and results (precip,evap,runof,dstor)
c        are returned for the last year.

c     pr    = monthly precip rates for nmon months (mm/day) (supp)
c     ta    = monthly air temps (deg C) (supp)
c     ho    = mid-month daylight hours (0 to 24) (supp)
c     alatd = latitude (deg N) (supp)
c     alond = longitude (deg E) (supp)
c     precip= annual mean precip (mm/day) (ret) 
c     evap  = annual mean evap  (mm/day) (ret) 
c     runof = annual mean runoff, sfc+base (mm/day) (ret) 
c     dstor = annual mean d[soil moisture content + snow] (mm/day) (ret)
c     iflag = 0 for ice grid, 1 for global (supp, for messages only)

      dimension ta(nmon), pr(nmon), ho(nmon)

      dimension mday(0:13)
      data mday / -15.,  16.,  45.,  75., 105., 136., 166.,
     *                  197., 228., 258., 289., 319., 350.,381. /
      save mday
      parameter (daypm = 30.)

c     if (iflag.eq.1 .and. abs(alatd).lt.20.) then
c        write (6,*) 'top of bucket'
c     endif

      conste = 1.
      consts = 1.
      constb = 1.

c     tune to make total Antarc runoff = that of gcm
      conste = 0.7
      consts = 1.3
      constb = 1.3

c     nyear =  5.
      nyear = 10.

c        Calculate zit and a (Huang et al, 1996, J. Clim, 9, 1350)

      zit = 0.
      do n=1,nmon
        if (ta(n).gt.0.) then
          zit = zit + (ta(n)/5.)**1.514
        endif
      enddo
      a = (6.75e-7)*(zit**3) - (7.71e-5)*(zit**2) + (1.79e-2)*zit +.49

      wmax = 760.
      w = 300.
      snow = 0.
c+++++++++++++++++++++++
      do iyear = 1,nyear
        do iday = 1,365
c+++++++++++++++++++++++

          if (iday.eq.1) then
            wprev = w
            snowprev = snow
            precip = 0.
            evap = 0.
            runof = 0.
          endif

c            Find interpolation indices for iday (like intdrive)
c            and interpolate to current time

          do im=1,nmon+1
            if (iday.le.mday(im)) then
              im1 = im-1
              im2 = im
              wm1 = (mday(im2)-iday) / (mday(im2)-mday(im1))
              if (im1.eq.0) im1 = nmon
              if (im2.eq.nmon+1) im2 = 1
              goto 102
            endif
          enddo
  102     continue
          p =  wm1*pr(im1) + (1.-wm1)*pr(im2)
          t =  wm1*ta(im1) + (1.-wm1)*ta(im2)
          h =  wm1*ho(im1) + (1.-wm1)*ho(im2)

c            Apply snowfall to snow
              
          if (t.lt.0.) then
            snow = snow + p
            p2 = 0.
          else
            p2 = p
          endif

c            Account for snow melt

c         smelt = min (5.*max(t,0.), snow)           ! pdd, .005 m/day/C
          smelt = min (7.*max(t,0.), snow)           ! pdd, .007 m/day/C
          snow = snow - smelt
          p2 = p2 + smelt

c            Compute potential evaporation

          if (t.lt.0.) then
            ep = 0.
          else if (t.lt.26.5) then
            ep = 16. * (h/12.) * ((10.*t/zit)**a) / daypm
          else 
            ep = (-415.85 + 32.25*t - 0.43*t*t) / daypm
          endif

c             Do time step (predictor-corrector) 

          do istep = 1,2
            if (istep.eq.1) then
              wp = w
            else
              wp = 0.5*(w+wp)
            endif
            if (snow.gt.0.) then
              e = 0.
              s = 0.
              b = 0.
            else
              e = conste * ep * (wp/wmax)
              s = consts * p2 * ((wp/wmax)**4.886)
              b = constb * 0.093 * wp / daypm
            endif
            wp = w + p2 - e - s - b 
            if (istep.eq.1) then
              wp = max (0., min (wmax, wp))
            else
              w = wp
            endif
          enddo

c            Account for w overshoots

          if (w.gt.wmax) then
            s = s + w - wmax
            w = wmax
          else if (w.lt.0.) then
            s = s + w
            if (s.lt.0.) then
              b = b + s
              s = 0
              if (b.lt.0.) then
                e = e + b
                b = 0
              endif
            endif
            w = 0.
          endif
          r = s + b

c            Accumulate annual mean quantities 

          precip = precip + p
          evap   = evap   + e
          runof  = runof  + r

          if (iday.eq.365) then
            precip = precip/365.
            evap   = evap/365.
            runof  = runof/365.
            dstor = (w + snow - wprev - snowprev) / 365.
          endif

c++++++++++++++
        enddo
      enddo
c++++++++++++++

      return
      end

c-----------------------------------------------------------------------

      subroutine siweather1 (ta, nmon, run, silicate, siweatherability)

c        Calculates silicate weathering (mol CO2/yr/m2) for one land pt,
c        following Godderis et al., 2010, Geoderma.

c     ta = monthly temps for nmon months (deg C) (supp)
c     run = annual mean runoff (mm/day) (supp)
c     silicate = weathering rate (mol CO2/yr/m2) (ret)
c     siweatherability = Si weatherability constant (supp)

      dimension ta(nmon)
      dimension ndaypm(12)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm 

      annualtemp = 0.
      do n=1,nmon
        annualtemp = annualtemp + ta(n)*ndaypm(n)/365.
      enddo

      if (annualtemp.lt.-10.) then
        silicate = 0.                                      ! polar
      else if (annualtemp.lt.10.) then
        silicate = 0.1027 * run * (365./1000.)             ! cold
      else if (annualtemp.lt.21.) then
        silicate = 0.2703 * run * (365./1000.)             ! mild
      else
        if ( (run*365./1000.).lt.0.2 ) then
          silicate = 0.6877 * run * (365./1000.)           ! warm arid
        else
          silicate = 0.0497 * run * (365./1000.)           ! warm humid
        endif
      endif

c     silicate = silicate * (6.      /9.295296) * siweatherability
c     silicate = silicate * (6.      /9.041016) * siweatherability
c     silicate = silicate * (6.942435/9.041016) * siweatherability
c     silicate = silicate * (6.      /7.813698) * siweatherability
      silicate = silicate                       * siweatherability

      return
      end

c-----------------------------------------------------------------------

      subroutine siweather2 (ta, nmon, run, silicate, siweatherability)

c        Calculates silicate weathering (mol CO2/yr/m2) for one land pt,
c        following Donnadieu et al., 2006, G3.

c     ta = monthly temps for nmon months (deg C) (supp)
c     run = annual mean runoff (mm/day) (supp)
c     silicate = weathering rate (mol CO2/yr/m2) (ret)
c     siweatherability = Si weatherability constant (supp)

      dimension ta(nmon)
      dimension ndaypm(12)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm 

      annualtemp = 0.
      do n=1,nmon
        annualtemp = annualtemp + ta(n)*ndaypm(n)/365.
      enddo

      zt = annualtemp + 273.16
      zt0 = 273.16

      zcon = (4.34362/304.474) ! arbitrary  scaling
      rgas = 8.314472
      Ea_silw = 48200.

      silicate = zcon * run * exp ( -(Ea_silw/rgas)*( 1./zt - 1./zt0) )

      silicate = silicate * siweatherability

      return
      end

c----------------------------------------------------------------------

      subroutine orgsink (fluxorgcur, timeice)

c        If first call, read deep-sea-core data file for
c        rate of organic carbon sink vs time, ~30.080 to 34.837 Ma
c        (from Lee Kump's Site744ForgCalc.[xls,txt] file).
c        Interpolate to model time (timeice), return CO2 sink rate.

c        fluxorgcur = organic sink rate, mol CO2 per year (ret)
c        timeice = years from present (-ve for BP) (supp)   
c        iu = unit number for input file (supp)

#include <comicepath.h>

      character*200 cfilin
      parameter (norg=374)
      dimension timeorg(norg), fluxorg(norg)
      logical firstorg
      data firstorg /.true./
      save timeorg, fluxorg, firstorg

      iu = 51

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_orgsink/Site744ForgCalc.txt'

c        If first call, read data file (derived from deep-sea core,
c        from Lee Kump's xls file)
      
c-----------------------
      if (firstorg) then
c-----------------------
        write (6,'(a)') 'Reading organic sink file'
        open (iu, file=cfilin, status='old')
        read (iu,'(//)') 
        do i=norg,1,-1                        ! so timeorg runs forwards
          read (iu,*,err=200,end=200) 
     *      timeorg(i), dummy, fluxorg(i)     ! dummy is d18O
          timeorg(i) = -1.e6*timeorg(i)       ! times to (-) years BP
        enddo
        close (iu)

        do i=norg,1,-1
          write (7,'(i6, f15.1, e15.3)') i, timeorg(i), fluxorg(i)
        enddo

        firstorg = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.le.timeorg(1)) then
c       fluxorgcur = fluxorg(1)
        fluxorgcur = 0.
        return
      else if (timeice.ge.timeorg(norg)) then
c       fluxorgcur = fluxorg(norg)
        fluxorgcur = 0.
        return
      else
        do i=1,norg-1
          if (timeice.ge.timeorg(i) .and. timeice.le.timeorg(i+1))then
            zwei = (timeice-timeorg(i)) / (timeorg(i+1)-timeorg(i))
            fluxorgcur = (1.-zwei)*fluxorg(i) + zwei*fluxorg(i+1)
            return
          endif
        enddo
      endif
      write (6,'(a)') "Error in orgsink - shouldn't get here"
      stop

  200 write (6,'(a,i6)') 
     *  'Error 200 in reading organik sink input file: i=', i
      stop

      end

c----------------------------------------------------------------------

c****
#else
c****

      subroutine doco2 ()
      return
      end

c*****
#endif
c*****
