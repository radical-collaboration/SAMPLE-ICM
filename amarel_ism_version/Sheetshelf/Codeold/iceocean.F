      subroutine ocean (h, hs, hb, hw, hbal, maskwater,
     *                  budgsnow, budgevap, budgrain, 
     *                  budgmelt, baseperc, basefrml, 
     *                  oceanmelt, facemelt, msector,
     *                  budgall, budgneg2, ubal, vbal,
     *                  tracer, arhap, visc,
     *                  div, calvice, crevliq, arcocn,
     *                  ifnexo, fedge, hedge, tempocn,
#if defined MELANGE
     *                  hmel,
#endif
     *                  ua, va, fracgu, fracgv,
     *                  timeice, sealev, ro18, dtantjan, weirun, 
     *                  dtocn, ifrest, iloop, iffirse, iffirsr, itry)

c        Set basal melting (oceanmelt) for floating ice,
c        vertical face melt (facemelt) if FACEMELT, and 
c        calving rate (calvice) if CALVING.

#include <comicephys.h>
#include <comicegrid.h>

#if defined (ISMIP6BMB)
#  include <comicepath.h>
#  include <netcdf.inc>
      character*120 cfile, cvarin
      integer*4 ncid, rcode, varid, start(3), count(3)
      character*16 cvar
      real*4 bmbanom(nx,ny)
      logical firstismip
      data firstismip /.true./
      save firstismip, bmbanom
#endif

#if defined (LARMIP)
#  include <comicepath.h>
#  include <netcdf.inc>
      character*120 cfile, cvarin
      integer*4 ncid, rcode, varid, start(3), count(3)
      character*16 cvar
      real*8 reglarmip(nx,ny)
      logical firstlarmip
      data firstlarmip /.true./
      save firstlarmip, reglarmip
#endif

      dimension
     *  h(nx,ny),            hs(nx,ny),          
     *  hb(nx,ny),           hw(nx,ny),    
     *  hbal(nx,ny),         maskwater(nx,ny),
     *  budgsnow(nx,ny),     budgevap(nx,ny),     budgrain(nx,ny),     
     *  budgmelt(nx,ny),     baseperc(nx,ny),     basefrml(nx,ny),
     *  oceanmelt(nx,ny),    facemelt(nx,ny),     msector(nx,ny),      
     *  budgall(nx,ny),      budgneg2(nx,ny),
     *  ubal(0:nxp,0:nyp),   vbal(0:nxp,0:nyp),
     *  tracer(nx,ny,0:nlevp,ntrace), arhap(nx,ny), visc(nx,ny),
     *  div(nx,ny),          calvice(nx,ny),      crevliq(nx,ny),
     *  arcocn(nx,ny),       
     *  ifnexo(nx,ny),       fedge(nx,ny),        hedge(nx,ny),       
     *  tempocn(nx,ny), 
#if defined MELANGE
     *  hmel(nx,ny),
#endif
     *  ua(0:nxp,0:nyp),     va(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp), fracgv(0:nxp,0:nyp)

      parameter (nlook = 4)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
c     E-W/N-S:
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /

      dimension distocn(nx,ny), distgl(nx,ny)
      dimension work(nx,ny), work2(nx,ny), work3(nx,ny)
      save distocn, distgl
c     parameter (dtarc = 0.)   ! yr-interval to call oceanarc (0=always)
      parameter (dtarc = 10.)

#if defined (OCMARTIN)
c     statement fnct for ocean temperature vs depth (no longer used):
c     ftoc (z, ta, za, tb, zb) = 
c    *     max (ta, min (tb, ta + (tb-ta)*(z-za)/(zb-za) ))
      character*16 csec
#endif

#if defined (LIU2D)
c     for o18 -> time lookup in calling liu_read. Smoothed d18O vs time,
c     last 20 kyrs, Lisiecki+Raymo:
      parameter (ntab = 18)
      dimension otab(0:ntab), ttab(0:ntab)
      save otab, ttab
      data ttab / 
     *   0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
     *        11.,  12.,  13.,  14.,  15.,  16.,  17.,  18./
      data otab / 
     * 3.23, 3.25, 3.27, 3.28, 3.30, 3.31, 3.33, 3.37, 3.40, 3.45, 3.52,
     *       3.60, 3.92, 4.06, 4.28, 4.49, 4.75, 4.88, 5.02/
#endif

#if defined (LIUCORR) || defined (RCPCORR) || defined (LOVEANOM) || defined (CSMANOM)
c     for bias correction vs modern obs
      dimension dtcor(nx,ny)
      logical first_dtcor
      data first_dtcor /.true./
      save dtcor, first_dtcor
#endif

#if defined (EISANTA) && defined (CLEMCLIM)
      parameter (nep=9)
      dimension
     *  alonep(nep), alatep(nep), xep(nep), yep(nep),
     *  allep(nx,ny)
      data alatep
     * / -80., -75., -72.5, -67.5, -66.5,  -65., -61., -55., -89. /
      data alonep
     * / -90., -70.,  -65.,  -67.,  -65., -62.5, -50.,  -1.,  -1. /
      save allep 
      logical firstep
      save firstep
      data firstep /.true./
#endif

      logical firstosec
      data firstosec /.true./
      save firstosec

c        For Antarctica, set geographical sectors (msector), sector
c        names (csector), and total number (nsec).
c        (csector, nsecmax, nsec are in comicegrid.h).
c        Used only for DTOCNMATCO2ASE and PITFUTOCNASE.

      if (firstosec) then
        call osector (msector, maskwater)
        firstosec = .false.
      endif

c        Set ifnexo=1 for ocean shelf points (h>0) adjacent to 
c        open ocn pts, else = 0 (for setfhedge, and returned for
c        cliffmelt in icedyn

      call izero (ifnexo, nx*ny)
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1 .and. h(i,j).eq.0.) then
            ifnexo(i,j) = 1
          else if (maskwater(i,j).eq.1 .and. h(i,j).gt.0.) then
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskwater(ii,jj).eq.1 .and. h(ii,jj).eq.0.)
     *          ifnexo(i,j) = 1
            enddo
          endif
        enddo
      enddo

c        Set sub-grid ice fraction, thickness([f,h]edge) at shelf edges.
c        Do here, since [f,h]edge may be needed by facemelt, calving 
c        below, and cliffmelt in icedyn.

#if ! defined (NOSUBGRIDEDGESHELF)
      call setfhedge (h, hb, hw, maskwater, ifnexo, fedge, hedge, 
     *                budgneg2, fracgu, fracgv, sealev, dtocn)
#else
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then
            fedge(i,j)  = 1.
            hedge(i,j)  = h(i,j)
          endif
        enddo
      enddo
#endif

c      If NOSHELF, ensure no floating ice and return 
c      (unlike MAXSHELFMELT below, which continues with calving, etc)

#if defined (NOSHELF)
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then
c           oceanmelt(i,j) = 1000. 
            oceanmelt(i,j) = 10000. 
            fedge(i,j) = 0.
            hedge(i,j) = 0.
            h(i,j) = 0.
          else    
            oceanmelt(i,j) = 0.
          endif   
          facemelt(i,j) = 0.
          div(i,j) = 0.
          calvice(i,j) = 0.
        enddo   
      enddo   
      return  
#endif

#if defined (TEST2D) 
c     for mismipplus:
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then

            zd = -rhoice*h(i,j)/rholiq
            zm1 = 0.2 * tanh(hw(i,j)/75.) * max (-100.-zd, 0.)

            zm2 = 0.
            if (xh(i,j).gt.480.e3) zm2 = 100.

#  if defined (MISMIPPLUS_spin)
            zo = 0.
#  elif defined (MISMIPPLUS_ctl)
            zo = 0.
#  elif defined (MISMIPPLUS_1rr)
            zo = zm1
#  elif defined (MISMIPPLUS_1ra)
            zo = zm1
            if (timeice.gt.100.) zo = 0.
#  elif defined (MISMIPPLUS_2rr)
            zo = zm2
#  elif defined (MISMIPPLUS_2ra)
            zo = zm2
            if (timeice.gt.100.) zo = 0.
#  else
            write (6,*) 'Error (ocean): must specify MISMIPPLUS_...'
            stop
#  endif

          else
            zo = 0.
          endif

          oceanmelt(i,j) = zo
          facemelt(i,j) = 0.
          div(i,j) = 0.
          calvice(i,j) = 0.
        enddo
      enddo
      return
#endif

c        Calculate degrees of arc with straight-line access to deep ocn,
c        and distance to open ocean and distance to grounding line. 
c        Used below for OCMARTIN and/or CALVING.
c        For nested runs, arcocn is optionally not calculated (if
c        -DKEEPARCOCN, to avoid problems with arc to open ocean for 
c        small domains), and has been read from nestdrive file and 
c        passed here. 

#if ( defined (OCMARTIN) || defined (CALVING) )  && ! ( defined (EISLINE) || defined (TRANSECTA) || defined (TESTMEL2D) )
      n = max (1, nint (dtarc/dtocn))
      if (mod (iloop, n).eq.0 .or. iffirsr.eq.1) 
     *  call oceanarc (arcocn, distocn, distgl, h, maskwater
#  if defined (MELANGE)
     *                 ,hmel
#  endif
     *                )
#else
      call resetr (arcocn, nx*ny, 360.)
#endif

c        Set climate weights zw[lgm,mod,hot]
c        (used for OCMARTIN and CALVING)

      zclim = weirun  ! weirun was set in longterm (0:1:2 = lgm:mod:hot)
      zclim = max (0., min (2., zclim))
      if (zclim.lt.1.) then
        zwmod = zclim
        zwlgm = 1.-zwmod
        zwhot = 0.
      else
        zwhot = zclim - 1.
        zwmod = 1.-zwhot
        zwlgm = 0.
      endif


c        Begin oceanmelt calculations

c=======================================================================

c>>>>>>>>>>>>>>>>>>>>>>>>>
#if defined (MAXSHELFMELT)
c>>>>>>>>>>>>>>>>>>>>>>>>>

c        Set huge oceanmelt value, but continue below and set fedge, 
c        facemelt, calving, if any

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then

            if (MAXSHELFMELT.eq.1) then
              oceanmelt(i,j) = 400.
            else 
              if (MAXSHELFMELT.le.999) then
                oceanmelt(i,j) = MAXSHELFMELT
              else
                oceanmelt(i,j) = MAXSHELFMELT
     *                         * max (0., min (1., timeice/10.))
              endif
            endif

          else
            oceanmelt(i,j) = 0.
          endif
        enddo
      enddo

c=======================================================================

c>>>>>>>>>>>>>>>>>>>>>>
#elif defined (EISLINE)
c>>>>>>>>>>>>>>>>>>>>>>

      do j=1,ny
        do i=1,nx
          oceanmelt(i,j) = 0.
          if (maskwater(i,j).eq.1) then

#  if defined (LINEB)
c           oceanmelt(i,j) = (1.-weirun)*1.00 + weirun*0.
            oceanmelt(i,j) = 0.
#  elif defined (LINEC)
c           oceanmelt(i,j) = (1.-weirun)*0.45 + weirun*0.  ! ..*0.20 +..
            oceanmelt(i,j) = 0.
#  elif defined (LINED)
            oceanmelt(i,j) = (1.-weirun)*0.20 + weirun*0.
#  elif defined (LINEF)
            oceanmelt(i,j) = 0.
#  elif defined (LINEH1) || defined (LINEH2)
            oceanmelt(i,j) = 0.
#  elif defined (LINEH3) || defined (LINEH4)
            oceanmelt(i,j) = 0.                            ! 0., 4.
            zr = sqrt (xh(i,j)**2 + yh(i,j)**2)
            if (zr.gt.1920.e3) oceanmelt(i,j) = 100.
#  elif defined (SLMODEL) && defined (SLOCEANMELT)
            oceanmelt(i,j) = SLOCEANMELT
#  elif defined (SLMODEL) && ! defined (SLOCEANMELT)
            oceanmelt(i,j) = 0.
#  endif

          endif
        enddo
      enddo

c=======================================================================

c>>>>>>>>>>>>>>>>>>>>>>>>
#elif defined (TESTMEL2D)
c>>>>>>>>>>>>>>>>>>>>>>>>

      do j=1,ny
        do i=1,nx
          oceanmelt(i,j) = 0.
          if (maskwater(i,j).eq.1) then
            oceanmelt(i,j) = 1.
#  if defined (TESTMELPROF)
c           oceanmelt(i,j) = 0.                                   ! 1-D
c           if (float(i)/float(nx).gt.0.70) oceanmelt(i,j) = 200. ! 1-D
            oceanmelt(i,j) = 15.                                  ! 2-D
#  endif
          else
            oceanmelt(i,j) = 0.
          endif
        enddo
      enddo

c=======================================================================

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#elif defined (EISANTA) && defined (SHELFINV)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then
            if (hbal(i,j).gt.0.) then
c             don't reset following adaptive-timestep blowup (itry>1):
              if (iffirsr.eq.1 .and. itry.eq.1) oceanmelt(i,j) = 0.2
              zdm = (h(i,j)-hbal(i,j)) / 50.
              oceanmelt(i,j) = max (-2., min (100., oceanmelt(i,j)+zdm))
            else
              oceanmelt(i,j) = 100.
            endif
          else
            oceanmelt(i,j) = 0.
          endif
        enddo
      enddo

c=======================================================================

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#elif defined (EISANTA) && defined (CLEMCLIM)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c        If first call, set array (allep) using polygon defining 
c        eastern side of Antarctic Peninsula (allep = 1, else = 0.
c        Polygon vertices are alatep,alonep, set in data stmts above).

      if (firstep) then
        do m=1,nep
          call ps_xy (xep(m), yep(m),
     *                alatep(m)*pi/180, alonep(m)*pi/180.,
     *                0., 0.)
        enddo

c       set allep array (1=inside E. Ant. Peninsula polygon, 0=outside):
        call zero (allep, nx*ny)
        call polygon_wais (allep, xh,yh,nx,ny, xep,yep, nep,1.)

c       diagnostic dump:
        write (133,'(/a)') 'Eastern Antarctic Peninsula:'
c       do j=ny,1,-1
c         write (133,'(1000i1)') (nint(allep(i,j)),i=1,nx)
c       enddo
        call printmap (timeice,allep,'allep', 1., 133, 0)
        ier = flush (133)

        firstep = .false.
      endif

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then   ! ocean

            if (nint(allep(i,j)).eq.1) then
c             eastern side of peninsula:
              oceanmelt(i,j) = 10. 
c fut*                       + 40.*max (0., min (1., timeice/100.))

            else 
c             western side of peninsula:
              oceanmelt(i,j) = 10. 
c fut*                       + 40.*max (0., min (1., timeice/100.))

            endif

          else
            oceanmelt(i,j) = 0.
          endif
        enddo
      enddo

c=======================================================================

c>>>>>>>>>>>>>>>>>>>>>>>
#elif defined (OCMARTIN)
c>>>>>>>>>>>>>>>>>>>>>>>

c        Ocean melt rate based simply on ocean water temperature minus 
c        freezing pt. of seawater at ice base, following Martin et al.
c        (PISM-PIK, The Cryosphere, 2010).

c        Various options for getting ocean temperatures:
c        LIU2D, RCP*, LOVECLIP, CSMCLIM, or Levitus (default).
c        First, get tempocn = current ogcm temps, propagated under 
c        shelves, at one ocean level that is hard-coded in 
c        liu_read/rcp_read/levitus_read/loveclim_ocean_read/
c        csmocean_read, using nearest dataset temp. at a given depth.

c(((((((((((((((((((((((
#  if defined (LOVECLIP)
c(((((((((((((((((((((((

      call loveclim_ocean_read (tempocn, maskwater, h, hw, vnull, 
     *                          0, 0, timeice)

clllllllllllllllllllllllll
#    if defined (LOVEANOM)
clllllllllllllllllllllllll

c        If requested, set bias correction (anomaly method) dtcor
c        (modern obs minus modern Loveclip), but done every time 
c        to interp to current maskwater (which doesn't correspond
c        to modern continents, but ignore that issue).

cc    if (first_dtcor) then

        call levitus_read (work,  maskwater, vnulla, 1, 0.)
        call resetr (dtcor, nx*ny, vnull)

        call loveclim_ocean_read (work2, maskwater, h, hw, vnullb, 
     *                            1, 1, 0.)

        do j=1,ny
          do i=1,nx
            if (work(i,j).ne.vnulla .and. work2(i,j).ne.vnullb)
     *        dtcor(i,j) = work(i,j) - work2(i,j)
          enddo
        enddo

c       call printmap (0., work, 'lev temp, corr (C)', 0.2, 155, 0)
c       call printmap (0., work2,'LC ctl temp, corr (C)', 0.2, 155, 0)
c       call printmap (0., dtcor, 'dtcor (C)', 0.2, 155, 0)
cc      first_dtcor = .false.
cc    endif

c     apply bias correction (anomaly method):
      do j=1,ny
        do i=1,nx
          if (tempocn(i,j).ne.vnull .and. dtcor(i,j).ne.vnull) then
            tempocn(i,j) = tempocn(i,j) + dtcor(i,j)
          endif
        enddo
      enddo

clllllllll
#    endif
clllllllll

c((((((((((((((((((((((((
#  elif defined (CSMCLIM)
c((((((((((((((((((((((((

      call csmocean_read (tempocn, maskwater, vnull, 0, 0, timeice)

cssssssssssssssssssssssss
#    if defined (CSMANOM)
cssssssssssssssssssssssss

c        If requested, set bias correction (anomaly method) dtcor
c        (modern obs minus modern CSM), but done every time 
c        to interp to current maskwater.

cc    if (first_dtcor) then

        call levitus_read (work,  maskwater, vnulla, 1, 0.)
        call resetr (dtcor, nx*ny, vnull)

        call csmocean_read (work2, maskwater, vnullb, 1, 1, 0.)

        do j=1,ny
          do i=1,nx
            if (work(i,j).ne.vnulla .and. work2(i,j).ne.vnullb)
     *        dtcor(i,j) = work(i,j) - work2(i,j)
          enddo
        enddo

c       call printmap (0., work, 'lev temp, corr (C)', 0.2, 155, 0)
c       call printmap (0., work2,'CSM ctl temp, corr (C)', 0.2, 155, 0)
c       call printmap (0., dtcor, 'dtcor (C)', 0.2, 155, 0)
cc      first_dtcor = .false.
cc    endif

c     apply bias correction (anomaly method):
      do j=1,ny
        do i=1,nx
          if (tempocn(i,j).ne.vnull .and. dtcor(i,j).ne.vnull) then
            tempocn(i,j) = tempocn(i,j) + dtcor(i,j)
          endif
        enddo
      enddo

csssssssss
#    endif
csssssssss

c((((((((((((((((((((((((
#  elif defined (CSMORDO)
c((((((((((((((((((((((((

      call csmordo_ocean_read (tempocn, maskwater, vnull, 0, timeice)

c((((((((((((((((((((((((
#  elif defined (GLACBAY)
c((((((((((((((((((((((((

      do j=1,ny
        do i=1,nx
c          tempocn(i,j) = 10.  ! e
           tempocn(i,j) = 20.  ! d
        enddo
      enddo

c((((((
#  else
c((((((

c        For t <= 0, LIU2D (if specified) or Levitus (default).
c        For t >  0, RCP*. (if specified) or Levitus (default).

c++++++++++++++++++++++++++++
      if (timeice.le.0.) then    
c++++++++++++++++++++++++++++

c......................
#    if defined (LIU2D) 
c......................

        oceandtliu = 0.

c         if model time outside Liu time span (before 22 ka),
c         look up time within 0 to 18 ka, with same o18 as current
c         o18 from Lisiecki+Raymo's 18O vs. time over last 18 kyrs, 
c         smoothed to be monotonic. Recall that ro18 = 0 for ~modern, 
c         ro18 = 1 for ~glacial, linearly related to o18 (its range 
c         between o18(0) and o18(18 ka)) in subr lisiecki.

        if (nint(timeice).lt.-22000) then

c         In longterm(icectl), weirun was set to 
c            1 -1*ro18 + 0.2*max(dtantjan,0.)
c         (Amplitude of dtantjan's influence compared to ro18 can be 
c         justified very approximately by comparing time series of ro18,
c         dtantjan and proximal-Antarctic-Liu temperatures over last 
c         ~20 kyrs.) If weirun is inside [0,1], i.e. LGM to modern, then
c         set ro18c (used below to deduce delta-o18 value). If resulting
c         weirun is outside [0,1], then apply excess as a uniform 
c         temperature increment (oceandtliu), and reset ro18c to 0 or 1.
c         Note that ro18c implicitly includes effect of dtantjan. 
c         nb: ro18 = 0 for modern, = 1 for lgm.
c         nb: weirun >=1 for modern and warmer, <= 0 for lgm and cooler.

          if (weirun.gt.1.) then
            ro18c = 0.
            oceandtliu = 4.*(weirun - 1.)  ! based on MIS31 ~2 C warming
          else if (weirun.lt.0.) then
            ro18c = 1.
c           oceandtliu = -4.*(0. - weirun)
            oceandtliu = 0.              ! better no cooling beyond LGM?
          else
            ro18c = 1. - weirun          !recover ro18 from longterm
            oceandtliu = 0.
          endif

c         ro18 -> o18, inverse of o18 -> ro18 in subr lisiecki (otab,
c         ntab settings above correspond to those in subr lisiecki):
          zo = otab(ntab)*ro18c + otab(0)*(1.-ro18c)
          if (zo.le.otab(0)) then
            ztime = ttab(0) 
          else if (zo.ge.otab(ntab)) then
            ztime = ttab(ntab) 
          else    
            do m=1,ntab
              if (zo.le.otab(m)) then
                ztime =((otab(m)-zo)*ttab(m-1) + (zo-otab(m-1))*ttab(m))
     *                / ( otab(m) - otab(m-1) )
                goto 100
              endif   
            enddo   
  100       continue
          endif   
          ztime = -ztime*1.e3                              ! +kyr to -yr

        else

          ztime = timeice

        endif

        call liu_read (tempocn, maskwater, vnull, 0, ztime)

c..............................................................
#    elif defined (RCP26) || defined (RCP45) || defined (RCP85)
c..............................................................

        call rcp_read (tempocn, maskwater, vnull, 0, timeice)

c........
#    else
c........

        call levitus_read (tempocn, maskwater, vnull, 0, timeice)

c..........
#     endif
c..........

c+++++++++++++++++++++++++++++++++
      else if (timeice.gt.0.) then    
c+++++++++++++++++++++++++++++++++

        ifforce = 0
        if (timeice.lt.1.1*dtocn) ifforce = 1

#    if defined (RCP26) || defined (RCP45) || defined (RCP85)
     
        call rcp_read (tempocn, maskwater, vnull, ifforce, timeice)

#    elif defined (LIU2D)

        oceandtliu = 0.
        call liu_read (tempocn, maskwater, vnull, ifforce, timeice)

#    else

        call levitus_read (tempocn, maskwater, vnull, ifforce, timeice)

#    endif

c++++++++++
      endif
c++++++++++

cccccccccccccccccccccccccccccccccccccccccccccc
#    if defined (LIUCORR) || defined (RCPCORR)
cccccccccccccccccccccccccccccccccccccccccccccc
c        If requested, set bias correction dtcor (modern obs minus 
c        modern Liu or RCP), but done every time to interp to current
c        maskwater. If additionally LIUCORR2, apply correction
c        only for mod/hot climates, not lgm (weighted).

cc    if (first_dtcor) then

        call levitus_read (work,  maskwater, vnulla, 1, 0.)
        call resetr (dtcor, nx*ny, vnull)

c------------------------------
        if (timeice.le.0.) then
c------------------------------
#      if defined (LIUCORR)
          call liu_read (work2, maskwater, vnullb, 1, 0.)
#      elif defined (RCPCORR)
          call rcp_read (work2, maskwater, vnullb, 1, 0.)
#      endif
          do j=1,ny
            do i=1,nx
              if (work(i,j).ne.vnulla .and. work2(i,j).ne.vnullb)
     *          dtcor(i,j) = work(i,j) - work2(i,j)
            enddo
          enddo
c-----------------------------------
        else if (timeice.gt.0.) then 
c-----------------------------------
#      if defined (RCPCORR)
          call rcp_read (work2, maskwater, vnullb, 1, 0.)
#      elif defined (LIUCORR)
          call liu_read (work2, maskwater, vnullb, 1, 0.)
#      endif
          do j=1,ny
            do i=1,nx
              if (work(i,j).ne.vnulla .and. work2(i,j).ne.vnullb)
     *          dtcor(i,j) = work(i,j) - work2(i,j)
            enddo
          enddo
c------------
        endif
c------------

c       call printmap (0., work, 'lev temp, corr (C)', 0.2, 155, 0)
c       call printmap (0., work2,'liu/rcp temp, corr (C)', 0.2, 155, 0)
c       call printmap (0., dtcor, 'dtcor (C)', 0.2, 155, 0)
cc      first_dtcor = .false.
cc    endif

c     apply bias correction:
      do j=1,ny
        do i=1,nx
          if (tempocn(i,j).ne.vnull .and. dtcor(i,j).ne.vnull) then
#       if defined (LIUCORR2)
            tempocn(i,j) = tempocn(i,j) + (zwmod+zwhot)*dtcor(i,j)
#       else
            tempocn(i,j) = tempocn(i,j) + dtcor(i,j)
#       endif
          endif
        enddo
      enddo
cccccccccc
#    endif
cccccccccc

c(((((((
#  endif
c(((((((

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c        Various uniform incremements to ocean temperatures. 
c        (Usually will be none or only one of them, but could
c        be oceandtliu (for times before 22 ka) and another).   
c        nb: tempocn constrained to be >= fr.pt. 400m, next loop below.

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then

            zdt = 0.


#  if defined (LIU2D)
            zdt = zdt + oceandtliu
#  endif

#  if defined (OCUNIFSHIFT)
            if (timeice.lt.0.) zdt = zdt - 5.*zwlgm + 2.*zwhot
#  endif

#  if defined (KN_DTO)
c           increase ocean temps uniformly (Knut,Nick):
            zdt = zdt + KN_DTO
#  endif

#  if defined (OCEANDT)
            zdt = zdt + OCEANDT
#  elif defined (OCEANDTRAMP)
            zdt = zdt + max (0.,min(1.,(timeice-1.e3)/1.e3))*OCEANDTRAMP
#  endif

#  if defined (OCEANDT_LIG)
c           LIG time-cts forcing,as DP16(timlig,dtolig in comicephys.h):
            if (timeice.le.timlig(1)) then
              ma = 1  
              mb = 1  
              zwa = 1.
            else if (timeice.ge.timlig(nlig)) then
              ma = nlig  
              mb = nlig  
              zwa = 1.
            else    
              do m=2,nlig
                if (timlig(m).ge.timeice) then
                  ma = m-1
                  mb = m  
                  zwa = (timlig(m)-timeice) / (timlig(m)-timlig(m-1))
                  goto 200
                endif   
              enddo   
            endif   
  200       zdt = zdt + zwa*dtolig(ma) + (1.-zwa)*dtolig(mb)  
#  endif

#  if defined (OCEANBA)
            if (timeice.gt.-14700. and. timeice.lt.-14200.) 
     *        zdt = zdt + OCEANBA
#  endif

#  if defined (DTOCNMATCO2) 
#    if defined (DTOCNMATCO2ASE)
            isec = msector(i,j)
            csec = csector(isec)
            if (csec.eq.'West') zdt = zdt + DTOCNMATCO2*zwhot
#    else    
            zdt = zdt + DTOCNMATCO2*zwhot
#    endif   
#  endif

#  if defined (PITFUTOCN)
#    if defined (PITFUTOCNASE)
            isec = msector(i,j)
            csec = csector(isec)
            if (csec.eq.'West') zdt = zdt + zwhot*PITFUTOCN
#    else
            zdt = zdt + zwhot*PITFUTOCN
#    endif
#  endif

            tempocn(i,j) = tempocn(i,j) + zdt

          endif
        enddo
      enddo

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c        Now calculate ocean melt, given tempocn(i,j) (= zto below). 
c        Possibly modify zto, and set ocfac, depending on climate 
c        (lgm/mod/hot), arcocn, bathymetry (distal)...

c%%%%%%%%%%%%%%%%
      do j=1,ny 
        do i=1,nx 
c%%%%%%%%%%%%%%%%

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if (maskwater(i,j).eq.1) then
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
            zhz = hedge(i,j)*rhor ! depth of ice base, +ve downwards (m)
            ztf = .0939 - .057*34.5 - 7.64e-4*zhz ! fr.pt.seawater, base
            ztof= .0939 - .057*34.5 - 7.64e-4*400.! fr.pt.seawater,400 m
            tempocn(i,j) = max (tempocn(i,j), ztof) ! constrain ocn temp

c******************
#  if defined (NHA)
c******************

            zto = tempocn(i,j)
#    if ! defined (LOVECLIP)
            zto = zwlgm*(zto - 3.) + (zwmod + zwhot)*zto
#    endif
            ocfac = 1.

c           increased contact with warm water in distal (deep) ocean:
c           zwb = max (0., min (1., (sealev-hb(i,j)-2500.)/200.))
c           zws2 = max (0., min (1., (arcocn(i,j)-150.)/20.))
c           zwbs = zwb*zws2
c           ocfac = ocfac*(1.-zwbs) + 10.*zwbs

c or param zto = f(latitude, depth), not from dataset:
c           zflat = max (0., min (1., (50.-abs(alatd(i,j)))/20. ))
c           ztprox = -1.*(1.-zflat)  + 5.*zflat
c           zfdeep = max (0., min (1., (sealev-hb(i,j)-900.)/200. ))
c           zto = (1.-zfdeep)*ztprox + zfdeep*5.
c           ocfac = 1.

c************************
#  elif defined (EISANTA)
c************************

            zto = tempocn(i,j)

            ocfac = 3.

c           Increase ocfac for deep "distal" ocean (zwb).
c           Weight also using arcocn (zws) to avoid doing this 
c           for deep interior troughs (e.g. in Weddell Sea):
            zwb = max (0., min (1., (sealev-hb(i,j)-2500.)/200.))
            zws2 = max (0., min (1., (arcocn(i,j)-150.)/20.))
            zwbs = zwb*zws2
            ocfac = ocfac*(1.-zwbs) + 10.*zwbs

c           reduction due to shallow seas (helps Abbot Ice Shelf):
            zws = max (0., min (1., (sealev-hb(i,j)-300.)/200. ))
            ocfac = ocfac * zws

c******
#  else
c******

            zto = tempocn(i,j)

            ocfac = 1.

c*******
#  endif
c*******

c           decrease ocfac for confined seaways (default zarc=20 deg):
#  if defined (OCFACARC)
            if (OCFACARC.eq.1) then
              zarc = 20.
            else
              zarc = OCFACARC 
            endif
#  else
            zarc = 20.
#  endif
            zws = max (0., min (1., (arcocn(i,j)-zarc)/20.))
            ocfac = ocfac * (1.0*zws + 0.1*(1.-zws))

#  if defined (RWRAP)
            isec = msector(i,j)
            csec = csector(isec)
            if (csec.eq.'West') then
              ocfac = ocfac * rwrap_ocfacmult * rwrap_ocfacmultase
            else
              ocfac = ocfac * rwrap_ocfacmult 
            endif
#  elif defined (OCFACMULT)
            ocfac = ocfac * OCFACMULT
#  endif 

#  if defined (OCFAC_JOHNSON)
c           Hudson Mts ~99.5 W, Mt.Murphy ~111.3 W (for -DPTS_JOHNSON).
c           108 W is between Thwaites and Pope Glacs.
            zlonw = -alond(i,j)             ! -180 to 180 deg E to deg W
            zwei = max (-1., min (1., -(zlonw-108.)/1. ))
            ocfac = ocfac * (2.**(OCFAC_JOHNSON*zwei))
#  endif 

c              Set oceanmelt

            zch2o = 4.218e3       ! specific heat of water (J/kg/C) 
            zgamma = 1.e-4        ! exchange coeff (m/s) 
            zfmar = 5.e-3         ! another adjustable factor
            oceanmelt(i,j) = 365.*86400.                 ! m/s to m/yr
     *                     * ocfac                       ! my factor
     *                     * abs(zto-ztf)                ! quadratic 
c                                                        ! Martin et al:
     *                     * rholiq * zch2o * zgamma * zfmar * (zto-ztf)
     *                     / (hfus*rhoice)

            work(i,j)  = zto 
            work2(i,j) = zto - ztf
            work3(i,j) = ocfac

#  if defined (JAKOBOCN)
           oceanmelt(i,j) = JAKOBOCN
#  endif

c%%%%%%%%%%%%%
          else                                          ! maskwater != 1
c%%%%%%%%%%%%%

            oceanmelt(i,j) = 0.
            work(i,j)  = 100.
            work2(i,j) = 100.
            work3(i,j) = 0.

c%%%%%%%%%%%%%%
          endif
c%%%%%%%%%%%%%%

#  if defined (FJORD)
c           if (maskwater(i,j).eq.1 .and. hb(i,j).lt.-200.)        ! a
c    *         oceanmelt(i,j) = 1000.                              ! a
c           if (hb(i,j).lt.-200.) oceanmelt(i,j) = 1000.           ! b
c           if (hb(i,j).lt.-300.) oceanmelt(i,j) = 1000.           ! c
c           if (maskwater(i,j).eq.1) oceanmelt(i,j) = 200.         ! d
c
c           zoc = 200.                                             ! e
c           zoc = 1000.                                            ! f
            zoc = 500.                                             ! g
            if (maskwater(i,j).eq.1) oceanmelt(i,j) = zoc          !efg
            do look=1,nlook                                        !efg
              ii = max (1, min (nx, i + ish(look)))                !efg
              jj = max (1, min (ny, j + jsh(look)))                !efg
              if (maskwater(ii,jj).eq.1 .and. hb(i,j).lt.sealev)   !efg
     *          oceanmelt(i,j) = zoc                               !efg
            enddo
#  endif

c%%%%%%%%%%%%
        enddo
      enddo
c%%%%%%%%%%%%

      if (iffirsr.eq.1) then
c     if (mod(timeice+.001,500.).lt.0.05) then
        call printmap (timeice,work,  'zto (C)',     0.2, 197, 0)
        call printmap (timeice,work,  'zto (C)',     1.0, 197, 0)
        call printmap (timeice,work2, 'zto-ztf (C)', 0.2, 197, 0)
        call printmap (timeice,work3, 'ocfac',       0.2, 197, 0)
        call flush (197)
      endif

#  if defined (KN_GLOMELT)
c     set oceanmelt at grounding line (Knut,Nick):
      do j=1,ny
        do i=1,nx
          ifgl = 0
          if ( h(i,j).gt.0. .and.
     *           (maskwater(i,j).eq.0 .and. sealev.gt.hb(i,j))
     *       )  then
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskwater(ii,jj).eq.1) ifgl = 1
            enddo
          endif
          if (ifgl.eq.1) oceanmelt(i,j) = KN_GLOMELT
        enddo
      enddo
#  endif

c. . . . . . . . . . . 
#if defined (ISMIP6BMB)
c. . . . . . . . . . . 
c        Read in and apply ocean melt anomaly for ISMIP6 abmb experiment

      if (firstismip) then
        write (ioterm,'(a)') 'Reading ISMIP6 abmb anomaly file'
        cfile = cpathhome(1:lenchr(cpathhome)) //  
     *          'Data_ismip6/basal_melt_anomaly_16km.nc'
        rcode = nf_open(cfile, nf_nowrite, ncid)
        call checkrcode (rcode, cfile, 'open', 6000)
        cvar = 'abmb'
        start(1) = 1
        count(1) = nx 
        start(2) = 1
        count(2) = ny 
        rcode = nf_inq_varid (ncid, cvar, varid)
        call checkrcode (rcode, cfile, cvar, 7760)
        rcode = nf_get_vara_real (ncid, varid, start, count, bmbanom)
        call checkrcode (rcode, cfile, cvar, 7761)
        rcode = nf_close (ncid)
        call checkrcode (rcode, cfile, 'close', 7762)
        firstismip = .false.
      endif
       
      zf = max (0., min (1., timeice/40.))                ! 40 year ramp
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then
            oceanmelt(i,j) = oceanmelt(i,j) +zf*bmbanom(i,j)*910./rhoice
          endif
        enddo
      enddo
c. . .
#endif
c. . .

c!!!!!!!!!!!!!!!!!!!
#if defined (LARMIP)
c!!!!!!!!!!!!!!!!!!!
c        Read in and apply ocean melt anomaly for LARMIP experiments

      if (firstlarmip) then
        write (ioterm,'(a)') 'Reading LARMIP mask file'
        cfile = cpathhome(1:lenchr(cpathhome)) //  
     *          'Data_larmip/LARMIP_regions_initMIPgrid_16.nc'
        rcode = nf_open(cfile, nf_nowrite, ncid)
        call checkrcode (rcode, cfile, 'open', 6000)
        cvar = 'regions'
        start(1) = 1
        count(1) = nx 
        start(2) = 1
        count(2) = ny 
        rcode = nf_inq_varid (ncid, cvar, varid)
        call checkrcode (rcode, cfile, cvar, 7760)
        rcode = nf_get_vara_double (ncid, varid, start,count, reglarmip)
        call checkrcode (rcode, cfile, cvar, 7761)
        rcode = nf_close (ncid)
        call checkrcode (rcode, cfile, 'close', 7762)
        call printmap (timeice, reglarmip, 'larmip regions', 1., 19, 1)
        firstlarmip = .false.
      endif

#  if defined (LARMIPREGION) 
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then
            if ( (LARMIPREGION.eq.0 .and. nint(reglarmip(i,j)).ge.1)
     *           .or.
     *           (LARMIPREGION.ne.0 .and. nint(reglarmip(i,j)).eq.
     *                                    LARMIPREGION)
     *         ) oceanmelt(i,j) = oceanmelt(i,j) + LARMIP
          endif
        enddo
      enddo
#  endif
c!!!!!
#endif
c!!!!!

c  end of "if defined (OCMARTIN)" block

c=======================================================================

c>>>>
#else
c>>>>

      write (ioterm, '(4(/a))') '*** Error:',
     *  '     Need to specify an ocean melt treatment in make file:',
     *  '     MAXSHELFMELT, EISLINE, EISANTA+SHELFINV, OCMARTIN,'
     *  '     or TESTMEL[2D,1D]')
      stop

c>>>>>
#endif
c>>>>>

c=======================================================================

      call zero (facemelt,  nx*ny)

c+++++++++++++++++++++
#if defined (FACEMELT)
c+++++++++++++++++++++

#  if defined (RWRAP)
      zfacemeltrate = rwrap_facemeltrate
#  elif defined (FACEMELTRATE)
      zfacemeltrate = FACEMELTRATE
#  else
      if (FACEMELT.eq.1) then 
        zfacemeltrate = 10.
      else
        zfacemeltrate = FACEMELT
      endif
#  endif

      do j=1,ny
        do i=1,nx
c..........................................................
          if ( h(i,j).gt.0. .and. 
     *         (maskwater(i,j).eq.1 .or. sealev.gt.hb(i,j))
     *       )  then
c..........................................................

c           Set facemelt for melting on vertical ice faces physically 
c           in contact with ocean, using oceanmelt of neigboring cells.
c           nb: facemelt is applied at sides, not surface or base, so is
c           not multiplied by fedge in icectl, and normalized by cell
c           area here.

c           vertical face extent potentially exposed to ocean: 
            if (maskwater(i,j).eq.0) then
c             if grounded, no interior basal oceanmelt
              zdz = min ( max (sealev-hb(i,j),0.), h(i,j) )
            else
c             if floating, interior basal melt captured by oceanmelt
              zdz = hedge(i,j)*rhor                 ! hedge=h if fedge=1
            endif

            zfacem = 0.
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
c             vertical face melt if adjacent cell is open ocean: 
              if (.not.(ii.eq.i .and. jj.eq.j)) then      ! for flowline
                if (maskwater(ii,jj).eq.1 .and. h(ii,jj).eq.0.) then 
                  if (look.eq.1) zdx = dyu(i-1,j)
                  if (look.ge.2) zdx = dyu(i,j) 
                  if (look.eq.3) zdx = dxv(i,j-1)
                  if (look.ge.4) zdx = dxv(i,j) 
                  zfacem = zfacem 
     *                   + zfacemeltrate*oceanmelt(ii,jj)
     *                     *zdx*zdz/darea(i,j)
c                 set for open ocean in case small ice after icedyn:
                  if (fedge(i,j).eq.1.) then
cunif               zdz2 = min(zdz,200.) ! unif small thickness, 3 plcs
                    zdz2 = zdz
                    facemelt(ii,jj) = max ( facemelt(ii,jj),
     *                                zfacemeltrate*oceanmelt(ii,jj)
     *                                *zdx*zdz2/darea(ii,jj) )
                  endif
                endif
              endif
            enddo
            facemelt(i,j) = zfacem

c..............
          endif
c..............
        enddo
      enddo

c+++++
#endif
c+++++

c====================
#if defined (CALVING)
c====================

c        Ice shelf calving (currently including dep on divergence)

      if (iffirse.eq.1 .and. ifrest.eq.-1) then
c       if first timestep of cold start, don't have ua,va yet  
        call zero (div, nx*ny)
        call zero (calvice, nx*ny)
      else
        call calving (h, hs, fedge, hedge, maskwater, 
     *                budgsnow, budgevap, budgrain, 
     *                budgmelt, baseperc, basefrml, oceanmelt, 
     *                budgall,  ubal,     vbal, 
     *                tracer, arhap, visc,
     *                arcocn, ua, va, 
     *                div, calvice, crevliq, 
     *                timeice, dtocn, zwlgm, zwmod, zwhot)
      endif

c====
#else
c====

      call zero (div, nx*ny)
      call zero (calvice, nx*ny)

c=====
#endif
c=====

      return
      end

c-----------------------------------------------------------------------

      subroutine osector (msector, maskwater)

c     Set geographical sector numbers (msector), names (csector),
c     and total number (nsec), for OCMARTIN (only for Antarctica).
c     (csector, nsec are in comicegrid.h).

#include <comicegrid.h>

      dimension msector(nx,ny)

c==========================================
#if defined (OCMARTIN) && defined (EISANTA)
c==========================================

      nsec = nsecmax   ! 4 for OCMARTIN

      do j=1,ny
        do i=1,nx
          zlon = alond(i,j)
          if (zlon.gt.180.) zlon = zlon - 360.
          zlat = alatd(i,j)

c*****

          if (
     *      (zlon.le. -10. .and. zlon.gt. -65.)
     *      .or.
     *      (zlon.le. -65. .and. zlon.gt. -90. .and. zlat.le. -75.) 
     *      .or.
     *      (zlon.le. -90. .and. zlon.gt.-120. .and. zlat.le. -85.)
     *      ) then
            msector(i,j) = 1
            csector(msector(i,j)) = 'Weddell'

c*****

          else if (
     *      (zlon.le. -65. .and. zlon.gt. -90. .and. zlat.gt. -75.) 
     *      .or.
     *      (zlon.le. -90. .and. zlon.gt.-120. .and. zlat.gt. -85.)
     *      .or.
     *      (zlon.le.-120. .and. zlon.gt.-140. .and. zlat.gt. -77.) 
     *      ) then
            msector(i,j) = 2
c           W. Peninsula., Bellinghausen, Amundsen:
            csector(msector(i,j)) = 'West'  

c*****

          else if (
     *      (zlon.le.-120. .and. zlon.gt.-140. .and. zlat.le. -77.) 
     *      .or.
     *      (zlon.le.-140. .or. zlon.gt. 160.) 
     *      ) then
            msector(i,j) = 3
            csector(msector(i,j)) = 'Ross'

c*****

          else if (zlon.le. 160. .and. zlon.gt. -10) then
            msector(i,j) = 4
c           most of East Antarctica:
            csector(msector(i,j)) = 'East'
          endif

c*****
        enddo
      enddo

c====
#else
c====

      nsec = 1
      call reseti (msector, nx*ny, 1)
#  if defined (EISLINE) || defined (TRANSECTA)
      csector(1) = 'flowline'
#  else
      csector(1) = ' '
#  endif

c=====
#endif
c=====

      return
      end

c-----------------------------------------------------------------------

      subroutine oceanarc (arcocn, distocn, distgl, h, maskwater
#  if defined (MELANGE)
     *                     ,hmel
#  endif
     *                    ) 

c       For each ocean point, computes degrees of arc (0 to 360) with 
c       straight-line all-ocean paths to open ocean (domain edge). 
c       For nested runs, optionally don't calculate arcocn (-DNESTING
c       and -DKEEPARCOCN)...see comments above.

#include <comicephys.h>
#include <comicegrid.h>

      dimension arcocn(nx,ny), distocn(nx,ny), distgl(nx,ny),
     *          h(nx,ny), maskwater(nx,ny)
#if defined MELANGE
     *         ,hmel(nx,ny)
#endif
      parameter (narc=72)   ! number of directions (every  5 deg)
c     parameter (narc=36)   ! number of directions (every 10 deg)

      dimension angarc(narc), tanarc(narc) 
      logical firstarc
      data firstarc /.true./
      save angarc, tanarc, firstarc

      parameter (nlook = 9)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
      data ish /0, -1, -1, -1,  0,  1, 1, 1, 0/
      data jsh /0,  1,  0, -1, -1, -1, 0, 1, 1/

      if (firstarc) then
        do m=1,narc
          angarc(m) = -pi + (m-0.5)*(2.*pi)/narc 
          tanarc(m) = tan(angarc(m))
        enddo
        firstarc = .false.
      endif

      do jo=1,ny
        do io=1,nx

c========================================
          if (maskwater(io,jo).eq.1) then
c========================================

#if defined (JAKOBOCN)
            ifdo = 1
#else
c           don't do for ocean pts with no adjacent h or land, for speed
            ifdo = 0
            do look=1,nlook
              ii = max (1, min (nx, io + ish(look)))
              jj = max (1, min (ny, jo + jsh(look)))
              if (h(ii,jj).gt.0. 
#  if defined MELANGE
     *           .or. hmel(ii,jj).gt.0.
#  endif
     *           .or. maskwater(ii,jj).ne.1) ifdo = 1
            enddo
#endif

c------------------------------
            if (ifdo.eq.1) then 
c------------------------------

#if ! (defined (NESTING) && defined (KEEPARCOCN))
             arcocn(io,jo) = 0.
#endif
             distocn(io,jo) = 10000.e3
             distgl(io,jo)  = 10000.e3

c>>>>>>>>>>>>>>>>>>>>>>>
             do m=1,narc
c>>>>>>>>>>>>>>>>>>>>>>>
               ang = angarc(m)
               tang = tanarc(m)

c                 Either go along one cell in i direction, and use 
c                 nearest j to current direction (ifi=1),or v.v. 
c                 for j and i(ifi=0):

               if (abs(ang).ge.0.75*pi) then
                 ifi = 1 
                 idir = -1
               else if (ang.ge.-0.75*pi .and. ang.le.-0.25*pi) then
                 ifi = 0 
                 idir = -1
               else if (abs(ang).le.0.25*pi) then
                 ifi = 1 
                 idir = 1
               else 
                 ifi = 0 
                 idir = 1
               endif

               i = io
               j = jo
               ifallice = 1
               zdistocn = 10000.e3
               if (h(io,jo).eq.0.) zdistocn = 0.   

  100          if (ifi.eq.1) then
                 i = i + idir
                 if (i.lt.1.or.i.gt.nx) goto 200
                 j = nint ((i-io)*tang) + jo
                 if (j.lt.1.or.j.gt.ny) goto 200
               else
                 j = j + idir
                 if (j.lt.1.or.j.gt.ny) goto 200
                 i = nint ((j-jo)/tang) + io
                 if (i.lt.1.or.i.gt.nx) goto 200
               endif

               if (maskwater(i,j).ne.1) then 
c                dist to gl must have gone through all ice, no open ocn 
                 if (ifallice.eq.1) then
                   zdist = dx0 * sqrt(float((i-io)**2+(j-jo)**2)) !polst
                   distgl(io,jo) = min (distgl(io,jo), zdist)
                 endif
                 goto 300
               endif

               if (h(i,j).eq.0.) then
                 zdist = dx0 * sqrt(float((i-io)**2+(j-jo)**2))   !polst
c                distocn(io,jo) = min (distocn(io,jo), zdist)
                 zdistocn = min (zdistocn, zdist)
                 ifallice = 0
               else
c                dist to open ocn must not encounter shelf ice again
                 zdistocn = 10000.e3
               endif

               go to 100

  200          continue
#if ! (defined (NESTING) && defined (KEEPARCOCN))
               arcocn(io,jo) = arcocn(io,jo) + 360./narc
#endif
               distocn(io,jo) = min (distocn(io,jo), zdistocn)

  300          continue

c>>>>>>>>>>>>>>>>>
             enddo 
c>>>>>>>>>>>>>>>>>

c---------------
            else 
c---------------

c                Open ocean, no adjacent ice or land

#if ! (defined (NESTING) && defined (KEEPARCOCN))
              arcocn(io,jo) = 360.
#endif
              distocn(io,jo)= 0.
              distgl(io,jo) = 10000.e3

c----------------
            endif  
c----------------

c=============
          else  
c=============

c              Land (or land ice)

#if ! (defined (NESTING) && defined (KEEPARCOCN))
            arcocn(io,jo) = 0.
#endif
            distocn(io,jo) = 10000.e3
            distgl(io,jo)  = 10000.e3

c==============
          endif  
c==============

        enddo  
      enddo  

      return
      end

c-----------------------------------------------------------------------

c,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#if ! defined (NOSUBGRIDEDGESHELF)
c,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      subroutine setfhedge (h, hb, hw, maskwater, ifnexo, fedge, hedge,
     *                      budgneg2, fracgu, fracgv, sealev, dtocn)

c        Set sub-grid ice fraction (fedge), and expected thickness of 
c        shelf ice (hedge) , at shelf edges next to open ocean. 
c        hedge = max of ice thicknesses in surrounding cells (zh,
c        allowing for "typical" downstream thinning).
c        Similar to sub-grid setting (SUBGRIDEDGE) for grounded ice in 
c        climate_all, but only to reduce oceanmelt by fedge in icectl
c        (no sfc elev interps). Don't change fedge for grounded ice 
c        (-DSUBGRIDEDGE  uses same array).

c        Also sets hedge>0 for ocean h=0 pts next to h>0,
c        to set calving in case small ice after icedyn.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny),            hb(nx,ny),         
     *  hw(nx,ny),           maskwater(nx,ny),
     *  ifnexo(nx,ny),       fedge(nx,ny),        hedge(nx,ny),
     *  budgneg2(nx,ny),     fracgu(0:nxp,0:nyp), fracgv(0:nxp,0:nyp)

      parameter (nlook = 4)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
c     E-W/N-S:
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /

      call zero (budgneg2, nx*ny)

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if (h(i,j).gt.0. and. ifnexo(i,j).eq.0) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              hedge(i,j) = h(i,j)
              fedge(i,j) = 1.
c~~~~~~~~~~~~~~~
            else    ! do also for ocn h=0, to set hedge for calving,
                    ! in case small ice after icedyn
c~~~~~~~~~~~~~~~
              zhav = 0.
              nav = 0
c............................
              do look=1,nlook
c............................
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (.not.(ii.eq.i .and. jj.eq.j)) then    ! for flowline

                  zh = 0.
                  if (maskwater(ii,jj).eq.0 .and. h(ii,jj).gt.0.) then
                    if (j.eq.jj) then
                      zfgl = fracgu (min(i,ii), j)
                    else
                      zfgl = fracgv (i, min(j,jj))
                    endif
                    hbg = hb(ii,jj)*(1.-zfgl) + hb(i,j)*zfgl
                    zh = (sealev-hbg)/rhor
                  else if (maskwater(ii,jj).eq.1 .and. h(ii,jj).gt.0.
     *                     .and. ifnexo(ii,jj).eq.0) then
                    zh = h(ii,jj)
                  endif

                  if (zh.gt.0.) then
cunif               zh3 = zh  ! unif small thickness, 3 plcs
c                   "downstream thinning, typical shelf profiles"...:
                    zh2 = zh * max (0.30, exp(-dd0/100.e3))
c                   for small coverage, weight towards zh, not zh2:
                    zwei =  min (1., h(i,j)/zh2)
                    zh3 = zh*(1.-zwei) + zh2*zwei 
                    zhav = max(zhav,zh3)          ! max, not average
c                   zhav = zhav + zh3             ! max, not average
                    nav = nav + 1
                  endif

                endif
c..................
              enddo
c..................
              if (nav.gt.0) then
c               zhav = zhav/nav                   ! max, not average
cunif           zhav = min (zhav, 200.) ! unif small thickness, 3 plcs
                hedge(i,j) = max (zhav, h(i,j))
                fedge(i,j) = h(i,j)/hedge(i,j)       ! so fedge=0 if h=0
              else if (h(i,j).gt.0.) then
                hedge(i,j) = h(i,j)
                fedge(i,j) = 1.
              else if (h(i,j).eq.0.) then
                hedge(i,j) = 0.
                fedge(i,j) = 0.
              endif
c~~~~~~~~~~~~~~~~
            endif
c~~~~~~~~~~~~~~~~

c           impose minimum limit (1 m) on sub-grid ice thickness: 7/5/13
            if (hedge(i,j).lt.1.) then
              budgneg2(i,j) = h(i,j)/dtocn
              h(i,j) = 0.
              fedge(i,j) = 0.
              hedge(i,j) = 0.
            endif

          endif
        enddo
      enddo

      return
      end

c,,,,,
#endif
c,,,,,

c-----------------------------------------------------------------------

c++++++++++++++++++++
#if defined (CALVING)
c++++++++++++++++++++

      subroutine calving (h, hs, fedge, hedge, maskwater, 
     *                    budgsnow, budgevap, budgrain, 
     *                    budgmelt, baseperc, basefrml, oceanmelt, 
     *                    budgall,  ubal,     vbal,
     *                    tracer, arhap, visc,
     *                    arcocn, ua, va,
     *                    div, calvice, crevliq,
     *                    timeice, dtocn,zwlgm, zwmod, zwhot)

#include <comicephys.h>
#include <comicegrid.h>
#include <comiceplot.h>

c       Calving parameterization for ice shelf pts.
c       Mechanical effects based on divergence, + climatic effects.

c     passed:
      dimension
     *  h(nx,ny),          hs(nx,ny),           fedge(nx,ny),        
     *  hedge(nx,ny),      maskwater(nx,ny),  
     *  budgsnow(nx,ny),   budgevap(nx,ny),     budgrain(nx,ny),     
     *  budgmelt(nx,ny),   baseperc(nx,ny),     basefrml(nx,ny),
     *  oceanmelt(nx,ny),
     *  budgall(nx,ny),    ubal(0:nxp,0:nyp),   vbal(0:nxp,0:nyp),
     *  tracer(nx,ny,0:nlevp,ntrace), arhap(nx,ny), visc(nx,ny),
     *  arcocn(nx,ny),     ua(0:nxp,0:nyp),     va(0:nxp,0:nyp),
     *  div(nx,ny),        calvice(nx,ny),      crevliq(nx,ny)

c     local:
      dimension 
     *  dudx(nx,ny), dvdy(nx,ny), maskcalv(nx,ny), ifopoc(nx,ny),
     *  work(nx,ny)

c     local, for search:
      dimension
     *  list(nx*ny), listprev(nx*ny)
      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook4 = 4)   ! 4 for E-W/N-S
      parameter (nlook8 = 8)   ! 8 for diagonals too

      call zero (dudx,  nx*ny)
      call zero (dvdy,  nx*ny)
      call zero (div,   nx*ny)
      call izero (maskcalv, nx*ny)

c        Set maskcalv (calving or not) flag, and divergence (div), for 
c        shelf pts. For fedge < 1.0, use free-floating value (f[hedge]) 

c        First, find open ocean pts (h=0) contiguous to domain edge 
c        (not polynyas), used for maskcalv=2:

      call findopenocn (ifopoc, maskwater, h)

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then

c-------------------------------------
            if (hedge(i,j).gt.0.) then
c-------------------------------------
c             maskcalv = 1 for any floating ice (including h=0 ocean pts
c             next to h>0, for which hedge>0):
              maskcalv(i,j) = 1

c             change maskcalv to 2 if next to *contiguous* open ocean:
c             do look=1,nlook4
              do look=1,nlook8                       ! nb: nlook8, not 4
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (.not.(ii.eq.i .and. jj.eq.j)) then    ! for flowline
c old:            if (maskwater(ii,jj).eq.1 .and. h(ii,jj).eq.0.) then
                  if (ifopoc(ii,jj).eq.1) then            ! new 11/14/14
                    maskcalv(i,j) = 2
                    goto 500
                  endif
                endif
              enddo
  500         continue
c----------------
            endif
c----------------

            if (h(i,j).gt.0. .and. fedge(i,j).eq.1) then
              dudx(i,j) = (ua(i,j)-ua(i-1,j))/dx(i,j)
              dvdy(i,j) = (va(i,j)-va(i,j-1))/dy(i,j)
              zdiv1 = dudx(i,j) + dvdy(i,j)
            else
              zdiv1 = 0.
            endif

c           free-floating div (cf. thetacalc):
            zdiv2 = ( arhap(i,j)*rhoip*grav*hedge(i,j)/4. ) ** powi

c           impose physical limit on du/dx+dv/dy (<= free-floating):
            zdiv1 = min (zdiv1, zdiv2)

            if (fedge(i,j).lt.1.) then
              div(i,j) = zdiv2
            else
              div(i,j) = zdiv1
            endif

          endif
        enddo
      enddo

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c Set extra crevasse depth due to water in surface crevasses, 
c dep on rain+sfc melt. Do everywhere, also used in cliffcalc (icedyn.F)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#if defined (CALVLIQ) || defined (RWRAP)
#  if defined (ZMAR)
      zmar = ZMAR
#  else
c     zmar = 1.0                                         ! EA nibble
c     zmar = 2.0                                         ! no EA nibbles
      zmar = 1.5                                         ! new 11/14/14
#  endif
      do j=1,ny
        do i=1,nx

          zm = max (0., (budgrain(i,j)+budgmelt(i,j))/1.)

c         crevliq(i,j) = zm**2                                 ! old std
cc        crevliq(i,j) = min (zm, baseperc(i,j)/1.) ** 2

c         reduce calving for small surface water, linear ramp to parab:
          if (zm.ge.2.*zmar) then
            crevliq(i,j) = zm**2
          else if (zm.ge.zmar) then
            crevliq(i,j) = 4.*zmar*(zm-zmar)
          else
            crevliq(i,j) = 0.
          endif

#  if defined (RWRAP)
          crevliq(i,j) = crevliq(i,j) * rwrap_calvliq
#  else
          crevliq(i,j) = crevliq(i,j) * CALVLIQ
#  endif
#  if defined (CALVLIQZWHOT)
     *                   * zwhot
#  endif

        enddo
      enddo
#else
      call zero (crevliq, nx*ny) 
#endif

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c Set calvice for ice shelf points
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c{{{{{{{{{{{{{{{{{{{{{{{
#if defined (CALVDAMAGE)
c{{{{{{{{{{{{{{{{{{{{{{{

c     following Bassis and Ma (EPSL, 2015). 
c     Track damage "r" (advected as tracer # 2, in icetherm),
c     where r = crevasse depths/ice thickness.

      call zero (calvice, nx*ny) 
      call zero (calvdam_t1, nx*ny)    ! diagnostic (in comiceplot)
      call zero (calvdam_t2, nx*ny)    ! diagnostic
      call zero (calvdam_t3, nx*ny)    ! diagnostic

      do j=1,ny
        do i=1,nx
          if (maskcalv(i,j).ge.1) then ! all floating ice

c           all layers are the same...just set zr = layer # 1.
c           at end, reset all layers.
            zr = tracer(i,j,1,2)

            zhedge = hedge(i,j)

c           B+M Eq.(11), neglecting beta terms:
            zn = powi

c           B+M Eq.(24):
            zrhos = 0.                         ! could depend on crevliq
            zmu0 = visc(i,j)/arhap(i,j)
            zs0 = 0.25 *
     *          * (rhoice-zrhos) * (rholiq-rhoice)*grav*zhedge
     *          / (max(zmu0,1.e-6)*max(div(i,j),1.e-6)*(rholiq-zrhos))
            zs0 = max (0., min (1., zs0))

c           as -budgall in icectl, no zfedge, no facemelt:
            zmdot = - (   budgsnow(i,j) - budgevap(i,j)
     *                  + budgrain(i,j) - baseperc(i,j) + basefrml(i,j)
     *                  - oceanmelt(i,j) )

c           B+M  Eq.(27):
            zterm1 = zn * (1.-zs0) * max(div(i,j),0.) 
            zterm2 = zmdot/zhedge
            zr = zr * exp((zterm1 + zterm2)*dtocn)

c---
c           small source term for damage creation:
c---
c either:
            ifnexl = 0
            do look=1,nlook4
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (.not.(ii.eq.i .and. jj.eq.j)) then    ! for flowline
                if (maskwater(ii,jj).eq.0) ifnexl = 1
              endif
            enddo
c           zterm3 = .01*ifnexl
            zterm3 = .05*ifnexl
            zr = max (zr, zterm3)
c---
c or:       zr = zr + f(fpin(hbsd,hw,...))*dtocn
c---
c or:
c           depth of sfc + basal crevasses due to div:
c           zhdiv = 2. * ( (max(div(i,j),0.)**powir)/arhap(i,j) )
c    *              / (rhoip*grav)
c           zterm3 = zhdiv/zhedge
ca          zr = max (zr, zterm3)
cb          zr = max (zr, 0.1*zterm3, .01)
c---

            zr = max (0., min (1., zr))

c           reset all tracer(#2) layers to new zr:
            do k=0,nlevp
              tracer(i,j,k,2) = zr
            enddo

            calvdam_t1(i,j) = zterm1    ! diagnostic
            calvdam_t2(i,j) = zterm2    ! diagnostic
            calvdam_t3(i,j) = zterm3    ! diagnostic

c           like default calving, below:
            zvmax = 3.e3
            zfc = 0.75
            calvice(i,j) = max (0., min (1., (zr-zfc)/(1.-zfc) ))
     *                   * zvmax*max(zhedge,5.)/dd0
c           based on actual remining undamaged thickness:
c           zvmax = 3.e3
c           zrem = (1.-zr)*zhedge
c           calvice(i,j) = calvice(i,j)
c    *                   + max (0., min (1., (10.-zrem)/10. ))
c    *                   * zvmax*max(zhedge,5.)/dd0

          endif ! maskcalv
        enddo
      enddo

c{{{{{{{{{{{{{{{{{{{{{{{{{{
#elif defined (CALVBALANCE)
c{{{{{{{{{{{{{{{{{{{{{{{{{{

c     following Amundson (J.Glac, 2015, in review).
c     Compare balance velocity (ubal or vbal, at interfaces)
c     with actual ice velocity (ua or va). Set calving veloc (zvc),
c     apply to ocean cells with fedge<1 adjacent to (ocn or no-ocn) 
c     fedge=1.
c     nb: ua,va are set to zero in icedyn unless one adjacent cell has 
c     fedge=1 (medge[u,v]), so can only do this calc for 
c     fedge=1:fedge<1 interfaces.

      alphabal = 1.2

      call zero (calvice, nx*ny) 
      call zero (calvbal_vi, nx*ny)    ! diagnostic (in comiceplot)
      call zero (calvbal_vb, nx*ny)    ! diagnostic
      call zero (calvbal_vc, nx*ny)    ! diagnostic

      call balvel (h, hs, budgall, ubal, vbal, timeice, dtocn)

      do j=1,ny
        do i=1,nx
          if ( maskwater(i,j).eq.1 .and. fedge(i,j).lt.1. 
c    *                             .and. fedge(i,j).gt.0. ! not open ocn
     *       ) then
            zvimax = 0.                 ! diagnostic
            zvbmax = 0.                 ! diagnostic
            zvcmax = 0.                 ! diagnostic
            do look=1,nlook4
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (.not.(ii.eq.i .and. jj.eq.j)) then      ! for flowline
                if (fedge(ii,jj).eq.1.) then

                  if (look.eq.1) then
                    zv  = ua(i-1,j)
                    zvb = ubal(i-1,j)
                    zdd = dyu(i-1,j)
                  else if (look.eq.2) then
                    zv  = ua(i,j)
                    zvb = ubal(i,j)
                    zdd = dyu(i,j)
                  else if (look.eq.3) then
                    zv  = va(i,j-1)
                    zvb = vbal(i,j-1)
                    zdd = dxu(i,j-1)
                  else if (look.eq.4) then
                    zv  = va(i,j)
                    zvb = vbal(i,j)
                    zdd = dxu(i,j)
                  endif
c                 like Amundson Eq.(7):
                  zvc = alphabal*abs(zv) - (alphabal-1.)*abs(zvb)
                  zvc = max (zvc, 0.)

                  zvimax = max (zvimax, abs(zv))    ! diagnostic
                  zvbmax = max (zvbmax, abs(zvb))   ! diagnostic
                  zvcmax = max (zvcmax, zvc)        ! diagnostic

c                 like default calving, below:
                  calvice(i,j) = calvice(i,j) 
     *                         + zvc * max(hedge(i,j),5.)*zdd/darea(i,j)

                endif  ! fedge(ii,jj).eq.1
              endif    ! look
            enddo      ! look
            calvbal_vi(i,j) = zvimax    ! diagnostic
            calvbal_vb(i,j) = zvbmax    ! diagnostic
            calvbal_vc(i,j) = zvcmax    ! diagnostic
          endif        ! maskwater.eq.1 .and. fedge.lt.1
        enddo
      enddo

c{{{{
#else
c{{{{

c     Default calving. Combo of:
c       Nye crevasse depths, 
c       simple accumulated damage param based on ice velocity,
c       minimum ice thickness (incl. arc),
c       crevliq.

      call zero (calvice, nx*ny) 

      do j=1,ny
        do i=1,nx

cmmmmmmmmmmmmmmmmmmmmmmmmm
#if defined (CALVEDGEONLY)
          if (maskcalv(i,j).eq.2) then ! hedge>0 next to contig open ocn
#else 
          if (maskcalv(i,j).ge.1) then ! all floating ice
#endif
cmmmmmmmmmmmmmmmmmmmmmmmmm

c              Dependency on sfc+basal crevasses (Nick et al, 2010)

            zhedge = hedge(i,j)

c           depth of sfc + basal crevasses due to div:
            zhdiv = 2. * ( (max(div(i,j),0.)**powir)/arhap(i,j) )
     *              / (rhoip*grav)

c---

            zhacc = 0.

c              Dependency on ice velocity (accumulated damage)

#  if ! defined (NOCALVACCU)
            if (maskcalv(i,j).eq.2) then  !incl opn ocn next to ice edge
c           if (maskcalv(i,j).eq.2 .and. h(i,j).gt.0.) then
              if (h(i,j).gt.0.) then
                zu = sqrt (   (0.5*(ua(i-1,j)+ua(i,j)))**2
     *                      + (0.5*(va(i,j-1)+va(i,j)))**2 )
              else
                zu = sqrt (   max (ua(i-1,j)**2, ua(i,j)**2)
     *                      + max (va(i,j-1)**2, va(i,j)**2) )
              endif
              zu = max (zu, .001)

c             nb: log(zu/1300.) = log(zu/100.) - log(1300./100.)
              zhacc = zhacc 
     *              + zhedge * max (0., log(zu/1600.)) / log(1.2)
c                   + zhedge * max (0., log(max(1.,(zu-600.)/200.)))
            endif
#  endif

c---

c               Dependency on minimum ice thickness, only at true 
c               shelf edge (h>0 next to h=0 and contig open ocean).
c               (good for modern, but problems with LGM expansion)

            if (maskcalv(i,j).eq.2) then  !incl opn ocn next to ice edge
c           if (maskcalv(i,j).eq.2 .and. h(i,j).gt.0.)  then
#  if defined (CALVACCH)
              zhsmall = CALVACCH
#  else
              zhsmall = 150.
c             zhsmall = 200. ! 11/14
#  endif
#  if defined (CALVACCHARC)
              if (CALVACCHARC.eq.1) then 
                zarc = 60.
c               zarc = 40.  ! 11/14
              else
                zarc = CALVACCHARC
              endif
#  else 
              zarc = 60.
c             zarc = 60.    ! 11/14
#  endif

#  if defined (MELANGE)
c             sinuosity effects on calving captured explicitly 
c             by back pressure of melange component:
              zws = 1.
#  else
c             implicit melange effect on calving:
              zws = max (0., min (1., (arcocn(i,j)-zarc)/20.))
#  endif
              zhsmall = zhsmall*zws

              zhacc = zhacc
     *              + zhedge * max (0., min (1., (zhsmall-zhedge)/50.))
            endif

c---

#  if defined (RWRAP)
            znick = rwrap_calvnick
#  elif defined (CALVNICK)
            znick = CALVNICK 
#  else
            znick = 1.
#  endif

#  if defined (CALVVMAX)
            zvmax = CALVVMAX
#  else
#    if defined (PREAGU13)
            zvmax = 12.e3  ! pre AGU13
#    else    
            zvmax =  3.e3  ! post AGU13
#    endif   
#  endif

#  if defined (CALVFRAC)
            zfc = CALVFRAC
#  else
            zfc = 0.75
#  endif

            zf = znick * (zhdiv + crevliq(i,j) + zhacc) / zhedge

c           use max of 5 m for ice thickness, to allow wastage to zero: 
            calvice(i,j) = max (0., min (1., (zf-zfc)/(1.-zfc) ))
     *                   * zvmax*max(zhedge,5.)/dd0

c---

cmmmmmmmmmmmmmm
          endif ! maskcalv
cmmmmmmmmmmmmmm

        enddo
      enddo

c{{{{{
#endif
c{{{{{

c#if (defined (EISLINE) || defined (TRANSECTA))
c     j = 1
c     i = 22
c     write (166, '(2a, i10,2i6, 7f10.4)') 
c    *  'timeice, i, maskcalv, fedg, hedge,
c    *  ' div, zf, zfc, zvmax, calv=', 
c    *  nint(timeice), i, maskcalv(i,j),
c    *  fedge(i,j), hedge(i,j),
c    *  div(i,j), zf, zfc, .001*zvmax, calvice(i,j)
c     call flush(166)
c---
c     zx = 0.
c     j = 1
c     do i=nx,1,-1
c       if (h(i,j).gt.0.) then
c         zx = dd0 * (i - 1 + fedge(i,1))
c         goto 900
c       endif
c     enddo
c 900 continue
c#endif

c---------------------
c     if (.true.) then
      if (.false.) then
c     if (mod(timeice+.001,1.).lt.0.05) then
c     if (mod(timeice+.001,100.).lt.0.05) then
c---------------------
      iu = 99
      write (iu,'(/a,f10.1/)') 'calving: timeice=',  timeice

      do j=1,ny
        do i=1,nx
          work(i,j) = maskcalv(i,j)
        enddo
      enddo
      call printmap (timeice,work,   'maskcalv',  0.1, iu, 1)

      do j=1,ny
        do i=1,nx
          work(i,j) = ifopoc(i,j)
        enddo
      enddo
      call printmap (timeice,work,   'ifopoc',   0.1, iu, 1)

      do j=1,ny
        do i=1,nx
          work(i,j) = ua(i,j)
        enddo
      enddo
      call printmap (timeice,work,   'ua_500',  500.0, iu, 1)
      do j=1,ny
        do i=1,nx
          work(i,j) = va(i,j)
        enddo
      enddo
      call printmap (timeice,work,   'va_500',  500.0, iu, 1)

      call printmap (timeice,dudx,   'dudx (1/y)', 1.e-3, iu, 0)
      call printmap (timeice,dvdy,   'dvdy (1/y)', 1.e-3, iu, 0)

      call printmap (timeice,div,    'div  (1/y)',    1.e-3, iu, 0)

      call printmap (timeice,h    ,'ice thickness (m)',      10., iu, 0)
      call printmap (timeice,hedge,'hedge (m)',              10., iu, 0)
      call printmap (timeice,fedge,'fedge (m)',              0.1, iu, 0)
      call printmap (timeice,arcocn, 'arcocn (deg)',          5., iu, 0)
      call printmap (timeice,calvice,'calvice(m)',           0.5, iu, 0)
      call flush (iu)
c----------
      endif
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine findopenocn (ifopoc, maskwater, h)

c       Find open ocean (no ice shelf) points that are contiguously
c       connected to domain boundaries, i.e., no polynyas. 
c       Sets ifopoc = 1 if so, 0 if not. 

#include <comicegrid.h>

      dimension ifopoc(nx,ny), maskwater(nx,ny), h(nx,ny) 
    
      dimension list(nx*ny), listprev(nx*ny)
      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook4 = 4)   ! 4 for E-W/N-S 
      parameter (nlook8 = 8)   ! 8 for diagonals too

      call izero (ifopoc, nx*ny)

c----------------------
      do iter = 1,nx*ny
c----------------------

        nlist = 0

        if (iter.eq.1) then

c         first iteration, open-ocean points at domain edges:
          do j=1,ny
            iskip = nx-1
            if (j.eq.1.or.j.eq.ny) iskip = 1
            do i=1,nx,iskip
              if (maskwater(i,j).eq.1 .and. h(i,j).eq.0.) then
                ifopoc(i,j) = 1
                nlist = nlist + 1
                list(nlist) = i + (j-1)*nx
              endif
            enddo
          enddo

        else

c         normal iteration, look at adjacent pts for new open-ocean pts:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook4
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (ifopoc(ii,jj).eq.0) then
                if (maskwater(ii,jj).eq.1 .and. h(ii,jj).eq.0.) then
                  ifopoc(ii,jj) = 1
                  nlist = nlist + 1
                  list(nlist) = ii + (jj-1)*nx
                endif
              endif
            enddo
          enddo

        endif

c          Skip out if no new open-ocean points found

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

      return
      end

c+++++
#endif
c+++++

c-----------------------------------------------------------------------

c##################
#if defined (LIU2D)
c##################

      subroutine liu_read (tempocn, maskwater, vnull, ifforce, timeice)

c     Reads Liu et al. (2009, Science) Netcdf ogcm data file. Times are
c     every 10 years from -22000 BP to +30 AP (2204 time points).
c     Lons (evely spaced) are indep of latitude, lats (unevenly spaced)
c     are indep of longitude, in Southern Hemisphere (not in NH!).
c
c     Temperatures are potential temperatures, should convert (Gsw...?)
c     (difference is approx 0.1 deg C per 1000 m depth).
c
c     Reads current ogcm temp vs.x,y (tempo) from file. Bilin. 
c     interpolates to ice grid (tempi), and propagate under shelves 
c     and to any liu-null regions using contiguous neighbor iteration. 
c     Returns as tempocn(nx,ny).

c     Each returned quantity uses ogcm temperatures at a particular
c     ogcm level (kleva,b), hard-coded below.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension 
     *  tempocn(nx,ny), maskwater(nx,ny)

      character cfilin*200, cvar*16
      integer*4 ncid, rcode, varid, start(4), count(4)
      save ncid, varid, cfilin

c     truncate last 3 time slices of Liu file (> 0 for ntliu=2202-4):
      parameter (nxliu=100, nyliu=116, nzliu=25, ntliu=2201)
      real*4 arr1d(10000), arr2d(nxliu,nyliu), brr2d(nxliu,nyliu)
      dimension lonliud(nxliu),     latliud(nyliu), 
     *          lonliu(nxliu),      latliu(nyliu), 
     *          zliu(nzliu),        tliu(ntliu),
     *          tempo(nxliu,nyliu),
     *          tempi(nx,ny)
      real lonliud, latliud, lonliu, latliu
c     for ocn_int:
      save tliu, lonliud, latliud, lonliu, latliu, tempi 
      real*4 tnull4, snull4 
      data tnull4, snull4 /9.96921e36, 9.96921e33/

c---
c     choose ocean layer indices (kleva, klevb)
c     z_t(1:25), from ncdump, in cm:
c     400.,      1222.025,   2108.881,   3100.537,   4239.197, 
c    5577.873,   7187.428,   9166.115,  11653.91,   14854.85,
c   19072.1,    24762.71,   32618.93,   43673.66,   59384.84, 
c   81538.86,  111251.3,   147480.4,   188381.7,   232692.8, 
c  279278.8,   327286.1,   376143.2,   425494.6,   475128.8 

c     surrounding 400 m level:
      parameter (kleva  = 13)   ! 32618.93 (converted to m below)
      parameter (klevb  = 14)   ! 43673.66 (converted to m below)
      parameter (zlevab = 400.) 

c     surrounding 500 m level:
c     parameter (kleva  = 14)   ! 43673.66 (converted to m below)
c     parameter (klevb  = 15)   ! 59384.84 (converted to m below)
c     parameter (zlevab = 500.) 
c---

      parameter (dtliu=10.)   ! time interval of data (yrs)

      logical firstliu, firstliuext
      data firstliu, firstliuext /.true., .true./
      save firstliu, firstliuext

c     diagnostic check output
      logical ifdiag
      data ifdiag /.true./
c     data ifdiag /.false./
      save ifdiag
      parameter (iudiag = 155)

      data it2dprev /-999/
      save it2dprev

      save weilev

      vnull = 1.e20                               ! vnull is passed back

c=======================
      if (firstliu) then
c=======================

        write (ioterm,'(a)')'reading existing Netcdf Liu OGCM data file'

c          Open existing Liu netcdf data file

        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *             '/Data_liu/b30.00_4kaDVTd.pop.ncrcat.ann.TEMP.nc'
        rcode = nf_open(cfilin, nf_nowrite, ncid)
        call checkrcode (rcode, cfilin, 'open', 8001)

c          ogcm depths (from surface to mid layer):

        start(1) = 1
        count(1) = nzliu

        rcode = nf_inq_varid (ncid,  'z_t',  varid)
        call checkrcode (rcode, cfilin, 'z_t', 8045)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr1d)
        call checkrcode (rcode, cfilin, 'z_t', 8046)
        if (ifdiag) write (iudiag,'(/a)') 'zliu:'
        do k=1,nzliu
          zliu(k) = arr1d(k)*.01  ! cm to m
          if (ifdiag) write (iudiag,'(i6,f10.3)') k, zliu(k)
        enddo

c       interpolate to zlevab from 2 surrounding Liu levels:
        weilev = (zlevab-zliu(kleva)) / (zliu(klevb)-zliu(kleva))
        if (ifdiag) write(iudiag,'(/a,i8,f10.3/a,i8,f10.3/a,f8.3)')
     *              'kleva  :', kleva, zliu(kleva),
     *              'klevb  :', klevb, zliu(klevb),
     *              'weilev :', weilev

c          ogcm times:

        start(1) = 1
        count(1) = ntliu

        rcode = nf_inq_varid (ncid,  'time',  varid)
        call checkrcode (rcode, cfilin, 'time', 8055)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr1d)
        call checkrcode (rcode, cfilin, 'time', 8056)
        if (ifdiag) write (iudiag,'(/a)') 'tliu:'
        do m=1,ntliu
          tliu(m) = arr1d(m)*1000.  ! -kyr BP to -yr BP
          if (ifdiag) write (iudiag,'(i6,f10.2)') m, tliu(m)
        enddo

c          ogcm latitudes (deg N), longitudes (deg E):

        start(1) = 1
        count(1) = nxliu

        start(2) = 1
        count(2) = nyliu

        rcode = nf_inq_varid (ncid,  'TLAT',  varid)
        call checkrcode (rcode, cfilin, 'TLAT', 8065)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
        call checkrcode (rcode, cfilin, 'TLAT', 8066)
        if (ifdiag) write (iudiag,'(/a)') 'latliu(d):'
        do j=1,nyliu
          latliud(j) = arr2d(1,j)
          latliu(j)  = latliud(j)*pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') j, latliud(j)
        enddo

        rcode = nf_inq_varid (ncid,  'TLONG',  varid)
        call checkrcode (rcode, cfilin, 'TLONG', 8067)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
        call checkrcode (rcode, cfilin, 'TLONG', 8068)
        if (ifdiag) write (iudiag,'(/a)') 'lonliu(d):'
        do i=1,nxliu
          lonliud(i) = arr2d(i,1)
          lonliu(i)  = lonliud(i)*pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') i, lonliud(i)
        enddo

        firstliu = .false.

c==========
      endif
c==========
     
c        Read current ogcm temp (x,y), interp to 
c        ice grid, propagate under shelves, return in tempocn

      it2d = nint((timeice-tliu(1))/dtliu) + 1
      it2d = max (1, min (ntliu, it2d))

c-------------------------------
      if (it2d.ne.it2dprev) then
c-------------------------------
        start(1) = 1
        count(1) = nxliu
        start(2) = 1
        count(2) = nyliu
        start(3) = kleva
        count(3) = 1
        start(4) = it2d
        count(4) = 1
        rcode = nf_inq_varid (ncid,  'TEMP',  varid)
        call checkrcode (rcode, cfilin, 'TEMP', 8070)

        start(3) = kleva
        rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
        call checkrcode (rcode, cfilin, 'TEMP', 8071)

        start(3) = klevb
        rcode = nf_get_vara_real (ncid, varid, start, count, brr2d)
        call checkrcode (rcode, cfilin, 'TEMP', 8072)

        do j=1,nyliu
          do i=1,nxliu
            if (arr2d(i,j).ne.tnull4 .and. brr2d(i,j).ne.tnull4) then
              tempo(i,j) = (1.-weilev)*arr2d(i,j) + weilev*brr2d(i,j)
            else
              tempo(i,j) = vnull
            endif
          enddo
        enddo

c          Interpolate to ice grid

        call ocn_int (tempo, lonliu, latliu, nxliu, nyliu, vnull, 
     *                tempi, maskwater, ifforce)

        it2dprev = it2d
c----------
      endif
c----------

c        Propagate under shelves and any other Liu-null areas.
c        Do every call in case maskwater changes.

      call scopy (nx*ny, tempi, 1, tempocn, 1)
      call ocn_ext (tempocn, maskwater, vnull)

      if (firstliuext) then
        call printmap (timeice, tempi,  'liu temp_orig (C)',1.,iudiag,0)
        call printmap (timeice, tempocn,'liu temp (C)',     1.,iudiag,0)
        call printmap (timeice, tempocn,'liu temp (C)',    0.2,iudiag,0)
      endif

      firstliuext = .false.
  
      return
      end

c#####
#endif
c#####

c-----------------------------------------------------------------------

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#    if defined (RCP26) || defined (RCP45) || defined (RCP85)
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      subroutine rcp_read (tempocn, maskwater, vnull, ifforce, timeice)

c     Reads RCP(8.5) Netcdf ogcm data file. Times are decadal averages
c     from 2010-2020 AD to 2290-2300 AD (29 time points).
c     Lons (evely spaced) are indep of latitude, lats (unevenly spaced)
c     are indep of longitude and only span Southern Hemisphere.
c
c     Reads current ogcm temp vs.x,y (tempo) from file. Bilin. 
c     interpolates to ice grid (tempi), and propagate under shelves 
c     and to any rcp-null regions using contiguous neighbor iteration. 
c     Returns as tempocn(nx,ny).

c     Each returned quantity uses ogcm temperatures at a particular
c     ogcm level (kleva,klevb), hard-coded below.

c     Similar to liu_read, but here do linear interpolation in time 
c     between 2 surrounding mid-decade time points of the data. 

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension 
     *  tempocn(nx,ny), maskwater(nx,ny)

      character cfilin*200, cvar*16
      integer*4 ncid, rcode, varid, start(4), count(4)
      save ncid, varid, cfilin

      parameter (nxrcp=320, nyrcp=187, nzrcp=60, ntrcp=29)
      real*4 arr1d(10000), arr2d(nxrcp,nyrcp), brr2d(nxrcp,nyrcp)
      dimension lonrcpd(nxrcp), latrcpd(nyrcp), 
     *          lonrcp(nxrcp),  latrcp(nyrcp), 
     *          zrcp(nzrcp),    trcp(ntrcp),
     *          tempo(nxrcp,nyrcp),
     *          tempi_a(nx,ny), tempi_b(nx,ny),  tempi(nx,ny)
      real lonrcpd, latrcpd, lonrcp, latrcp
c     for ocn_int:
      save trcp, lonrcpd, latrcpd, lonrcp, latrcp, tempi_a, tempi_b 

c---
c     choose ocean layer indices (kleva, klevb)
c     z_t(1:60), from ncdump, in cm:
c     500,    1500,    2500,     3500,   4500, 
c    5500,    6500,    7500,     8500,   9500, 
c   10500,   11500,   12500,    13500,  14500, 
c   15500,  16509.,  17547.,   18629.,  19766., 
c  20971.,  22257.,  23640.,   25137.,  26765., 
c  28548.,  30511.,  32686.,   35109.,  37822., 
c  40878.,  44337.,  48273.,   52772.,  57937., 
c  63886.,  70756.,  78700.,   87882.,  98470., 
c 110620., 124456., 140049.,  157394., 176400., 
c 196894., 218645., 241397.,  264900., 288938., 
c 313340., 337979., 362767., 387645.5, 412576., 
c 437539., 462519., 487508.,  512502., 537500. ;

c     surrounding 400 m level:
      parameter (kleva = 30)   ! 37822. (converted to m below)
      parameter (klevb = 31)   ! 40878. (converted to m below)
      parameter (zlevab = 400.) 

c     surrounding 500 m level:
c     parameter (kleva  = 33)   ! 48273. (converted to m below)
c     parameter (klevb  = 34)   ! 52772. (converted to m below)
c     parameter (zlevab = 500.) 
c---

      parameter (dtrcp=10.)   ! time interval of data (yrs)

      logical firstrcp, firstrcpext
      data firstrcp, firstrcpext /.true., .true./
      save firstrcp, firstrcpext

c     diagnostic check output
      logical ifdiag
      data ifdiag /.true./
c     data ifdiag /.false./
      save ifdiag
      parameter (iudiag = 155)

      data itrcpaprev, itrcpbprev /-999, -999/
      save itrcpaprev, itrcpbprev

      save weilev

      vnull = 1.e20                               ! vnull is passed back

c=======================
      if (firstrcp) then
c=======================

#if defined (RCP26)
        write (ioterm,'(a)')
     *    'reading existing Netcdf RCP2.6 ocean data file'
#elif defined (RCP45)
        write (ioterm,'(a)')
     *    'reading existing Netcdf RCP4.5 ocean data file'
#elif defined (RCP85)
        write (ioterm,'(a)')
     *    'reading existing Netcdf RCP8.5 ocean data file'
#else
        write (ioterm,'(/a)')
     *    '*** Error: missing RCP specifier (RCP26, RCP45, RCP85)'
        stop
#endif

c          Open existing rcp netcdf data file

        cfilin = cpathhome(1:lenchr(cpathhome)) //
#if defined (RCP26)
     *           '/Data_rcp2.6/decadal_avg_rcp2.6.nc'
#elif defined (RCP45)
     *           '/Data_rcp4.5/decadal_avg_rcp4.5.nc'
#elif defined (RCP85)
     *           '/Data_rcp8.5/decadal_avg_rcp8.5.nc'
#endif
        rcode = nf_open(cfilin, nf_nowrite, ncid)
        call checkrcode (rcode, cfilin, 'open', 8001)

c          ogcm depths (from surface to mid layer):

        start(1) = 1
        count(1) = nzrcp

        rcode = nf_inq_varid (ncid,  'z_t',  varid)
        call checkrcode (rcode, cfilin, 'z_t', 8045)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr1d)
        call checkrcode (rcode, cfilin, 'z_t', 8046)
        if (ifdiag) write (iudiag,'(/a)') 'zrcp (m):'
        do k=1,nzrcp
          zrcp(k) = arr1d(k)*.01  ! cm to m
          if (ifdiag) write (iudiag,'(i6,f10.3)') k, zrcp(k)
        enddo

c       interpolate to zlevab from 2 surrounding RCP levels:
        weilev = (zlevab-zrcp(kleva)) / (zrcp(klevb)-zrcp(kleva))
        if (ifdiag) write(iudiag,'(/a,i8,f10.3/a,i8,f10.3/a,f8.3)')
     *              'kleva  :', kleva, zrcp(kleva),
     *              'klevb  :', klevb, zrcp(klevb),
     *              'weilev :', weilev

c          ogcm times: m=1 to 29 for decades 2010-2020 to 2290-2300.
c          Convert from yr AD to yr AP (i.e., after 1950 AD).

        start(1) = 1
        count(1) = ntrcp

        if (ifdiag) write (iudiag,'(/a)') 'trcp (yr):'
        do m=1,ntrcp
          trcp(m) = 2010. + (m-.5)*10. - 1950.   
          if (ifdiag) write (iudiag,'(i6,f10.2,a)') 
     *                m, trcp(m) + 1950.,' AD'
        enddo

c          ogcm latitudes (deg N), longitudes (deg E, 0 to 360):
c          (file latitudes run north from -79.xxx to -0.1xxx)
c          (file longitudes run east from 320.xxx to 319.xxx)

        start(1) = 1
        count(1) = nxrcp

        start(2) = 1
        count(2) = nyrcp

        rcode = nf_inq_varid (ncid,  'TLAT',  varid)
        call checkrcode (rcode, cfilin, 'TLAT', 8065)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
        call checkrcode (rcode, cfilin, 'TLAT', 8066)
        if (ifdiag) write (iudiag,'(/a)') 'latrcp (deg):'
        do j=1,nyrcp
          latrcpd(j) = arr2d(1,j)
          latrcp(j)  = latrcpd(j)*pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') j, latrcpd(j)
        enddo

        rcode = nf_inq_varid (ncid,  'TLONG',  varid)
        call checkrcode (rcode, cfilin, 'TLONG', 8067)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
        call checkrcode (rcode, cfilin, 'TLONG', 8068)
        if (ifdiag) write (iudiag,'(/a)') 'lonrcp (deg):'
        do i=1,nxrcp
          lonrcpd(i) = arr2d(i,1)
          lonrcp(i)  = lonrcpd(i)*pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') i, lonrcpd(i)
        enddo

        firstrcp = .false.

c==========
      endif
c==========

c         Calculate data time indices and weight for 2 surrounding
c         mid-decades (itrcpa, itrcpb, wtrcpa)

      zt = (timeice-trcp(1))/dtrcp 
      itrcpa = nint(zt + 0.5)
      wtrcpa = itrcpa - zt  
      itrcpb = itrcpa + 1
      itrcpa = max (1, min (ntrcp, itrcpa))
      itrcpb = max (1, min (ntrcp, itrcpb))
      wtrcpa = max (0., min (1., wtrcpa))

c-------------------------------------------------------------
      if (itrcpa.ne.itrcpaprev .or. itrcpb.ne.itrcpbprev) then
c-------------------------------------------------------------
     
c           If time indices changed, read 2 new mid-decade ogcm temps, 
c           interp to ice grid, save in tempi_a, tempi_b

        start(1) = 1
        count(1) = nxrcp
        start(2) = 1
        count(2) = nyrcp
        start(3) = kleva
        count(3) = 1
c+++
c for tempi_a:
c+++
        start(4) = itrcpa
        count(4) = 1
        rcode = nf_inq_varid (ncid,  'TEMP',  varid)
        call checkrcode (rcode, cfilin, 'TEMP', 8070)

        start(3) = kleva
        rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
        call checkrcode (rcode, cfilin, 'TEMP', 8071)

        start(3) = klevb
        rcode = nf_get_vara_real (ncid, varid, start, count, brr2d)
        call checkrcode (rcode, cfilin, 'TEMP', 8072)

        do j=1,nyrcp
          do i=1,nxrcp
            if ( (.not.isnan(arr2d(i,j))) .and. (.not.isnan(brr2d(i,j)))
     *         ) then
              tempo(i,j) = (1.-weilev)*arr2d(i,j) + weilev*brr2d(i,j)
            else
              tempo(i,j) = vnull
            endif
          enddo
        enddo

c          Interpolate to ice grid 

        call ocn_int (tempo, lonrcp, latrcp, nxrcp, nyrcp, vnull, 
     *                tempi_a, maskwater, ifforce)
c+++
c for tempi_b:
c+++
        start(4) = itrcpb
        count(4) = 1
        rcode = nf_inq_varid (ncid,  'TEMP',  varid)
        call checkrcode (rcode, cfilin, 'TEMP', 8073)

        start(3) = kleva
        rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
        call checkrcode (rcode, cfilin, 'TEMP', 8074)

        start(3) = klevb
        rcode = nf_get_vara_real (ncid, varid, start, count, brr2d)
        call checkrcode (rcode, cfilin, 'TEMP', 8075)

        do j=1,nyrcp
          do i=1,nxrcp
            if ( (.not.isnan(arr2d(i,j))) .and. (.not.isnan(brr2d(i,j)))
     *         ) then
              tempo(i,j) = (1.-weilev)*arr2d(i,j) + weilev*brr2d(i,j)
            else
              tempo(i,j) = vnull
            endif
          enddo
        enddo

c          Interpolate to ice grid 

        call ocn_int (tempo, lonrcp, latrcp, nxrcp, nyrcp, vnull, 
     *                tempi_b, maskwater, ifforce)
c+++

        itrcpaprev = itrcpa
        itrcpbprev = itrcpb

c----------
      endif
c----------

c        Remainder done every timestep.
c        Linearly interpolate in time between surrounding mid-decades,
c        tempi_[a,b] (saved) to tempi.

      do j=1,ny
        do i=1,nx
          if (tempi_a(i,j).ne.vnull .and. tempi_b(i,j).ne.vnull) then
            tempi(i,j) = wtrcpa*tempi_a(i,j) + (1.-wtrcpa)*tempi_b(i,j)
          else
            tempi(i,j) = vnull
          endif
        enddo
      enddo

c        Propagate under shelves and any other rcp-null areas. (Do every
c        call in case maskwater changes, regardless of time interp).

      call scopy (nx*ny, tempi, 1, tempocn, 1)
      call ocn_ext (tempocn, maskwater, vnull)

      if (firstrcpext) then
        call printmap (timeice, tempi,  'rcp temp_orig (C)',1.,iudiag,0)
        call printmap (timeice, tempocn,'rcp temp (C)',     1.,iudiag,0)
        call printmap (timeice, tempocn,'rcp temp (C)',    0.2,iudiag,0)
      endif

      firstrcpext = .false.
  
      return
      end

c$$$$$
#endif
c$$$$$

c-----------------------------------------------------------------------

      subroutine levitus_read (tempocn, maskwater,vnull,ifforce,timeice)

c     Reads Levitus94 Netcdf data file for modern temperature.
c     Global 1x1 deg grid.
c     Bilin. interpolates to ice grid (tempi), and
c     propagate under shelves and to any null regions using
c     contiguous neighbor iteration. Returns as tempocn(nx,ny).
c     Returned ogcm temperatures at a particular level 
c     (klev), hard-coded below.
c
c     Optionally, read in World Ocean Atlas 2005 or 2009 netcdf file,
c     (selected by hard-coded cf,cx,cy,cz,ct below).
c     Same grid, different variable names.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension 
     *  tempocn(nx,ny), maskwater(nx,ny)

      character cfilin*200, cvar*16
      integer*4 ncid, rcode, varid, start(4), count(4)
      save ncid, varid, cfilin

      parameter (nxlev=360, nylev=180, nzlev=33)
      real*4 arr1d(10000), arr2d(nxlev,nylev)
      dimension lonlevd(nxlev), latlevd(nylev), 
     *          lonlev(nxlev),  latlev(nylev),  
     *          zlev(nzlev),
     *          tempo(nxlev,nylev),
     *          tempi(nx,ny)
      real lonlevd, latlevd, lonlev, latlev
      save lonlevd, latlevd, lonlev, latlev, tempi 

c---
c     all files have same lon, lat, z coords.
      character*80 cf, cx, cy, cz, ct
c---
c     Levitus 1994:
c     parameter (cf = 'temp_levitus94.nc')
c     parameter (cx = 'X')
c     parameter (cy = 'Y')
c     parameter (cz = 'Z')
c     parameter (ct = 'temp')
c     parameter (vn = 1.e20)
c---
c     WOA 2005:
c     parameter (cf = 'woa05_t00an1.nc')
c     parameter (cx = 'lon')
c     parameter (cy = 'lat')
c     parameter (cz = 'depth')
c     parameter (ct = 't00ann1')
c     parameter (vn = 1.e20)
c---
c     WOA 2009:
      parameter (cf = 'woa09_temperature_annual_1deg.nc')
      parameter (cx = 'lon')
      parameter (cy = 'lat')
      parameter (cz = 'depth')
      parameter (ct = 't_an')
      parameter (vn = 9.96921e36)
c---
c     choose ocean layer index (klev)
c     zlev(1:33), from ncdump, in m:
c        0,   10,   20,   30,   50,   75,  100,  125,  150,  200, 
c      250,  300,  400,  500,  600,  700,  800,  900, 1000, 1100, 
c     1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 
c     4500, 5000, 5500 ;
      parameter (klev = 13)   ! 400 m
c     parameter (klev = 14)   ! 500 m
c     parameter (klev = 15)   ! 600 m
c---

      logical firstlev, firstlevext
      data firstlev, firstlevext /.true., .true./
      save firstlev, firstlevext

c     diagnostic check output
      logical ifdiag
      data ifdiag /.true./
c     data ifdiag /.false./
      save ifdiag
      parameter (iudiag = 155)

      vnull = vn                                  ! vnull is passed back

c=======================
      if (firstlev) then
c=======================

        write (ioterm,'(3a)') 
     *    'reading existing Netcdf ', cf(1:lenchr(cf)), ' data file'

c          Open existing lev netcdf data file

        cfilin = cpathhome(1:lenchr(cpathhome)) // 
     *           '/Data_levitus94/' // cf
        rcode = nf_open(cfilin, nf_nowrite, ncid)
        call checkrcode (rcode, cfilin, 'open', 8001)

c          ocean depth levels (0 = surface, ...):

        start(1) = 1
        count(1) = nzlev

        rcode = nf_inq_varid (ncid, cz, varid)
        call checkrcode (rcode, cfilin, cz, 8045)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr1d)
        call checkrcode (rcode, cfilin, cz, 8046)
        if (ifdiag) write (iudiag,'(/a)') 'zlev:'
        do k=1,nzlev
          zlev(k) = arr1d(k)
          if (ifdiag) write (iudiag,'(i6,f10.3)') k, zlev(k)
        enddo
        if (ifdiag) then
          write(iudiag,'(a/i6,f10.3)') 'klev :', klev, zlev(klev)
        endif

c          latitudes (deg N), longitudes (deg E):

        start(1) = 1
        count(1) = nylev

        rcode = nf_inq_varid (ncid,  cy,  varid)
        call checkrcode (rcode, cfilin, cy, 8065)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr1d)
        call checkrcode (rcode, cfilin, cy, 8066)
        if (ifdiag) write (iudiag,'(/a)') 'latlev(d):'
        do j=1,nylev
          latlevd(j) = arr1d(j)
          latlev(j)  = latlevd(j)*pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') j, latlevd(j)
        enddo

        start(1) = 1
        count(1) = nxlev

        rcode = nf_inq_varid (ncid,  cx,  varid)
        call checkrcode (rcode, cfilin, cx, 8067)
        rcode = nf_get_vara_real (ncid, varid, start, count, arr1d)
        call checkrcode (rcode, cfilin, cx, 8068)
        if (ifdiag) write (iudiag,'(/a)') 'lonlev(d):'
        do i=1,nxlev
          lonlevd(i) = arr1d(i)
          lonlev(i)  = lonlevd(i)*pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') i, lonlevd(i)
        enddo

        firstlev = .false.

c==========
      endif
c==========
     
c        Read ocean temp (x,y), interp to ice grid, propagate under 
c        shelves, return in tempocn

      start(1) = 1
      count(1) = nxlev
      start(2) = 1
      count(2) = nylev
      start(3) = klev
      count(3) = 1
      start(4) = 1  ! time (only 1 at most)
      count(4) = 1
      rcode = nf_inq_varid (ncid,  ct,  varid)
      call checkrcode (rcode, cfilin, ct, 8071)
      rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
      call checkrcode (rcode, cfilin, ct, 8072)
      do j=1,nylev
        do i=1,nxlev
          tempo(i,j) = vnull 
          if (.not. isnan(arr2d(i,j))) then    ! works for real*4 and *8
            if (arr2d(i,j).gt.-10. .and. arr2d(i,j).lt.50.) then
              tempo(i,j) = arr2d(i,j)
            endif
          endif
        enddo
      enddo

c        Interpolate to ice grid

      call ocn_int (tempo, lonlev, latlev, nxlev, nylev, vnull, 
     *              tempi, maskwater, ifforce)

c        Propagate under shelves and any other lev-null areas.
c        Do every call in case maskwater changes.

      call scopy (nx*ny, tempi, 1, tempocn, 1)
      call ocn_ext (tempocn, maskwater, vnull)

      if (firstlevext) then
        call printmap (timeice, tempi,  'lev temp_orig (C)',1.,iudiag,0)
        call printmap (timeice, tempocn,'lev temp (C)',     1.,iudiag,0)
        firstlevext = .false.
      endif
  
      return
      end

c-----------------------------------------------------------------------

c:::::::::::::::::::::
#if defined (LOVECLIP)
c:::::::::::::::::::::

      subroutine loveclim_ocean_read (tempocn, maskwater, h, hw, vnull,
     *                                ifforce, ifctl, timeice)

c     Reads current Loveclim climate file for ocean temperatures.
c     Bilinearly interpolates to ice grid (first call(s) only, saved in 
c     tempicur or tempictl, depending on ifctl), and propagates 
c     under shelves and to any null regions using contiguous neighbor 
c     iteration (every call). Returns in tempocn(nx,ny).
c     ifforce = 1 to force (re)calculation of indices (in ocn_int).
c     ifctl = 0 to use subr  readloveclim, current field, or
c           = 1 to use entry readlovectl, for anomaly method
c             (-DLOVEANOM)

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <comiceloveclim.h>
#include <netcdf.inc>

      dimension 
     *  tempocn(nx,ny), maskwater(nx,ny), h(nx,ny), hw(nx,ny)

      dimension 
     * tempo(nlong,nlatg), oceanfrac(nlong,nlatg),
     * maskwaterdum(nx,ny),
     * tempicur(nx,ny), tempictl(nx,ny)
      save tempicur, tempictl

      logical firstcur, firstctl, firstout
      data firstcur, firstctl, firstout /.true., .true., .true./
      save firstcur, firstctl, firstout

c     diagnostic check output
      logical ifdiag
      data ifdiag /.true./
c     data ifdiag /.false./
      save ifdiag
      parameter (iudiag = 155)

c        If first call, read ocean temps (tempo) on Loveclim grid.
c        Interpolate to ice grid (all pts, maskwaterdum), save in tempi.
c        (along,alatg,vnull_lc, etc, in comicegrid.h, comiceloveclim.h).

      vnull = vnull_lc                            ! vnull is passed back

c-----------------------------------------------------------------------
      if ((ifctl.eq.0.and.firstcur) .or. (ifctl.eq.1.and.firstctl)) then
c-----------------------------------------------------------------------
        if (ifctl.eq.0) then
          call readloveclim (tempo, 'tocean', 1)
        else
          call readlovectl  (tempo, 'tocean', 1)    ! for anomaly method
        endif

c       protect against rounnd-off error vs. vnull: 
        do j=1,nlatg
          do i=1,nlong
            if (abs(tempo(i,j)).gt.1.e10) tempo(i,j) = vnull
          enddo
        enddo
c or:
c       Use Loveclim's oceanfrac to set null values in tempo,
c       don't rely on null values in tempo itself:
c       if (ifctl.eq.0) then
c         call readloveclim (oceanfrac, 'oceanfrac', 1)
c       else
c         call readlovectl  (oceanfrac, 'oceanfrac', 1)    ! anom method
c       endif
c       do j=1,nlatg
c         do i=1,nlong
c           if (oceanfrac(i,j).eq.vnull .or. oceanfrac(i,j).lt..001)
c    *        tempo(i,j) = vnull
c         enddo
c       enddo

        call reseti (maskwaterdum, nx*ny, 1)
        if (ifctl.eq.0) then
          call ocn_int (tempo, along,alatg,nlong,nlatg, vnull, tempicur,
     *                  maskwaterdum, ifforce)
        else
          call ocn_int (tempo, along,alatg,nlong,nlatg, vnull, tempictl,
     *                  maskwaterdum, ifforce)
        endif

        if (ifctl.eq.0) then
          firstcur = .false.
        else
          firstctl = .false.
        endif
c----------
      endif
c----------
     
c        Propagate into ice-grid ocean areas with null values 
c        (under shelves, or not ocean in Loveclim).
c        Do every call in case maskwater changes (could avoid by
c        propagating everywhere first call).

      if (ifctl.eq.0) then
        call scopy (nx*ny, tempicur, 1, tempocn, 1)
      else
        call scopy (nx*ny, tempictl, 1, tempocn, 1)
      endif

      call love_ocn_ext (tempocn, maskwater, h, hw, vnull)
c     call ocn_ext      (tempocn, maskwater, vnull)

c         For niciety (with LOVEANOM settings in calling subr) 
c         set tempocn back to null for current non-ocean pts

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).ne.1) tempocn(i,j) = vnull
        enddo
      enddo

      if (firstout) then
        call printmap (timeice, tempicur,'lc temp_orig (C)',1.,iudiag,0)
        call printmap (timeice, tempocn,'lc temp (C)',     1.,iudiag,0)
        firstout = .false.
      endif
  
      return
      end

c:::::
#endif
c:::::

c-----------------------------------------------------------------------

cOOOOOOOOOOOOOOOOOOOO
#if defined (CSMCLIM)
cOOOOOOOOOOOOOOOOOOOO

      subroutine csmocean_read (tempocn, maskwater, vnull, 
     *                          ifforce, ifanom, timeice)

c     Reads CSM data file for modern temperature.
c     Global lon,lat deg grid.
c     Bilin. interpolates to ice grid (tempi), and
c     propagate under shelves and to any null regions using
c     contiguous neighbor iteration. Returns as tempocn(nx,ny).
c     Returned ogcm temperatures at a particular level (kcsm), 
c     hard-coded below.
c     ifanom = 0 for experiment, = 1 for anomaly control.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension 
     *  tempocn(nx,ny), maskwater(nx,ny)

      character cfilin*200, cvar*16
      integer*4 ncid, rcode, varid, start(4), count(4)
      save ncid, varid, cfilin

c     parameter (nxcsm=nlong, nycsm=nlatg, nzcsm=60)
      parameter (nxcsm=nlong, nycsm=nlatg, nzcsm=1)
      real*4 arr1d(10000), arr2d(nxcsm,nycsm),  brr2d(nxcsm,nycsm)
      dimension loncsmd(nxcsm), latcsmd(nycsm), 
     *          loncsm(nxcsm),  latcsm(nycsm),
     *          zcsm(nzcsm),
     *          tempo(nxcsm,nycsm),
     *          tempi(nx,ny),
     *          work(1000)
      real loncsmd, latcsmd, loncsm, latcsm
      save loncsmd, latcsmd, loncsm, latcsm

      dimension timcsm(ntcsm)              ! ntcsm in comicegrid.h
      save timcsm                          ! timcsm not used for J. Tsai

c J. Tsai:
      dimension tempo_reg(nxcsm,nycsm), tempo_anom(nxcsm,nycsm)
      save tempo_reg,   tempo_anom, 
     *     itaprev_reg, itaprev_anom
      data itaprev_reg, itaprev_anom /-999999, -999999/

c---
c     choose ocean layer index (kcsm)
c     zcsm(1:60), from ncdump, in cm:
c 500, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500, 9500, 10500, 
c 11500, 12500, 13500, 14500, 15500, 16509.84,17547.9,18629.13,19766.03,
c 20971.14, 22257.83, 23640.88, 25137.02, 26765.42, 28548.37, 30511.92, 
c 32686.8, 35109.35, 37822.76, 40878.46, 44337.77, 48273.67, 52772.8, 
c 57937.29, 63886.26, 70756.33, 78700.25, 87882.52, 98470.59, 110620.4, 
c 124456.7, 140049.7, 157394.6, 176400.3, 196894.4, 218645.7, 241397.2, 
c 264900.1, 288938.5, 313340.5, 337979.3, 362767, 387645.2, 412576.8, 
c 437539.2, 462519, 487508.3, 512502.8, 537500 ;

c     parameter (kcsm = 31)   ! 408.78 m
      parameter (kcsm = 1)    ! 408.78 m
      data zcsm(1) / 40878.46 /
c---

      logical firstcsm, firstcsmext
      data firstcsm, firstcsmext /.true., .true./
      save firstcsm, firstcsmext

c     diagnostic check output
      logical ifdiag
      data ifdiag /.true./
c     data ifdiag /.false./
      save ifdiag
      parameter (iudiag = 155)

      vnull = 9.96921e+36                         ! vnull is passed back

c J. Tsai (all years, nyearstart must = 0, use just one current yr):
      ita = max (1, min (ntcsm, nint(timeice + 0.5) ))
c###################################################
      if (ifanom.eq.0 .and. ita.eq.itaprev_reg) then 
c###################################################
        call scopy (nxcsm*nycsm, tempo_reg, 1, tempo, 1)

c#########################################################
      else if (ifanom.eq.1 .and. ita.eq.itaprev_anom) then 
c#########################################################
        call scopy (nxcsm*nycsm, tempo_anom, 1, tempo, 1)

c#########
      else
c#########

c---

c..........................
      if (ifanom.eq.0) then
c..........................
        cfilin = cpathhome(1:lenchr(cpathhome)) // 
#  if defined (CSMEXP)
     *         'Data_csm/' // CSMEXP // '/' //
#  else
     *         'Data_csm/' //
#  endif
     *         'inputdata_OceanT.nc'                           ! J. Tsai
c...............................
      else if (ifanom.eq.1) then
c...............................
        cfilin = cpathhome(1:lenchr(cpathhome)) // 
#    if defined (CSMHIRES)
     *           'Data_csm/CESM1/Ctl/TOCN_CTL_Clim.nc'         ! J. Tsai
#    else
     *           'Data_csm/Ctl/TOCN_CTL_Clim.nc'               ! J. Tsai
#    endif
c..........
      endif
c..........

c     write (ioterm,'(3a,2i6)') 
c    *  'reading existing Netcdf ', cfilin(1:lenchr(cfilin)), 
c    *  ' data file', ifanom, ita  

c        Open existing CSM netcdf data file

      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 8001)

c=======================
      if (firstcsm) then
c=======================

c          ocean depth levels (0 = surface, ...):

c       start(1) = 1
c       count(1) = nzcsm
c       rcode = nf_inq_varid (ncid, 'z_t', varid)
c       call checkrcode (rcode, cfilin, 'z_t', 8045)
c       rcode = nf_get_vara_real (ncid, varid, start, count, arr1d)
c       call checkrcode (rcode, cfilin, 'z_t', 8046)
c       if (ifdiag) write (iudiag,'(/a)') 'zcsm:'
c       do k=1,nzcsm
c         zcsm(k) = arr1d(k)
c         if (ifdiag) write (iudiag,'(i6,f10.3)') k, zcsm(k)
c       enddo
c       if (ifdiag) then
c         write(iudiag,'(a/i6,f10.3)') 'kcsm :', kcsm, zcsm(kcsm)
c       endif

c          latitudes (deg N):

        start(1) = 1
        count(1) = nycsm
        start(2) = 1
        count(2) = 1

        rcode = nf_inq_varid (ncid,  'lat',  varid)
        call checkrcode (rcode, cfilin, 'lat', 8065)
        rcode = nf_get_vara_double (ncid, varid, start, count, latcsmd)
        call checkrcode (rcode, cfilin, 'lat', 8066)

c       protect interpolation logic in setinterp (needed?):
        latcsmd(1)     = max (latcsmd(1),     -89.999) 
        latcsmd(nycsm) = min (latcsmd(nycsm),  89.999) 

        if (ifdiag) write (iudiag,'(/a)') 'latcsm(d):'
        do j=1,nycsm
          latcsm(j) = latcsmd(j) * pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') j, latcsmd(j)
        enddo

c          longitudes (deg E):

        start(1) = 1
        count(1) = nxcsm
        start(2) = 1
        count(2) = 1

        rcode = nf_inq_varid (ncid,  'lon',  varid)
        call checkrcode (rcode, cfilin, 'lon', 8067)
        rcode = nf_get_vara_double (ncid, varid, start, count, work)
        call checkrcode (rcode, cfilin, 'lon', 8068)
        do i=1,nxcsm
c         shift cut from greenwich to dateline for ice model:
          ii = mod (i + nxcsm/2 -1, nxcsm)  + 1
          loncsmd(i) = work(ii)
          if (loncsmd(i).gt.180.) loncsmd(i) = loncsmd(i) - 360.
        enddo

c       protect interpolation logic in setinterp (needed?):
        if (loncsmd(1).gt. 179.999 .or. loncsmd(1).lt.-179.999) 
     *    loncsmd(1) = -179.999 
        if (loncsmd(nxcsm).gt. 179.999 .or. loncsmd(nxcsm).lt.-179.999) 
     *    loncsmd(nxcsm) = 179.999

        if (ifdiag) write (iudiag,'(/a)') 'loncsm(d):'
        do i=1,nxcsm
          loncsm(i)  = loncsmd(i) * pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') i, loncsmd(i)
        enddo

c          times (years) (ntcsm is in comicegrid.h):

        if (ntcsm.eq.1) then
          timcsm(1) = 0.
        else    
          start(1) = 1
          count(1) = ntcsm
          start(2) = 1
          count(2) = 1
          rcode = nf_inq_varid (ncid,  'year',  varid)  
          call checkrcode (rcode, cfilin, 'year', 8070)
          rcode = nf_get_vara_double (ncid, varid, start, count, timcsm)
          call checkrcode (rcode, cfilin, 'year', 8071)
          if (ifdiag) then
            write (iudiag,'(/a)') 'timcsm:' 
            do m=1,ntcsm
              write (iudiag,'(i6,f10.2)') m, timcsm(m) 
            enddo   
          endif
        endif

        firstcsm = .false.

c==========
      endif
c==========

c        Find surrounding time indices and weight (ita, itb, wta)

c---
c not J. Tsai:
c     if (timeice.le.timcsm(1)) then
c       ita = 1
c       itb = 1
c       wta = 1.
c     else if (timeice.ge.timcsm(ntcsm)) then
c       ita = ntcsm
c       itb = ntcsm
c       wta = 1.
c     else 
c       do m=2,ntcsm
c         if (timeice.le.timcsm(m)) then
c           ita = m-1
c           itb = m
c           wta = (timcsm(itb) - timeice) / (timcsm(itb) - timecsm(ita))
c           goto 10
c         endif
c       enddo
c  10   continue
c     endif
c---

c        Read ocean temp (x,y), interp to ice grid, propagate under 
c        shelves, return in tempocn. Read for 2 surrounding times, 
c        linear interp vs. time.

      start(1) = 1
      count(1) = nxcsm
      start(2) = 1
      count(2) = nycsm
      if (ifanom.eq.0) then
        start(3) = ita
      else if (ifanom.eq.1) then
        start(3) = 1
      endif
      count(3) = 1
c     start(3) = kcsm
c     count(3) = 1
c     start(4) = ita
c     count(4) = 1

      if (ifanom.eq.0) then
        cvar = 'TEMP_O_annual'
      else if (ifanom.eq.1) then
        cvar = 'TOCN_CTL'
      endif
      rcode = nf_inq_varid (ncid,  cvar,  varid)
      call checkrcode (rcode, cfilin, cvar, 8072)
      rcode = nf_get_vara_real (ncid, varid, start, count, arr2d)
      call checkrcode (rcode, cfilin, cvar, 8073)

c not J. Tsai:
c     start(4) = itb
c     rcode = nf_get_vara_real (ncid, varid, start, count, brr2d)
c     call checkrcode (rcode, cfilin, cvar, 8074)

      do j=1,nycsm
        do i=1,nxcsm
c         shift cut from greenwich to dateline for ice model:
          ii = mod (i + nxcsm/2 -1, nxcsm)  + 1
          tempo(i,j) = vnull 
c not J. Tsai:
c         if ( (.not. isnan(arr2d(ii,j))) .and.
c    *         (.not. isnan(brr2d(ii,j))) ) then !works for real*4 and 8
c           if ( arr2d(ii,j).gt.-10. .and. arr2d(ii,j).lt.50. .and.
c    *           brr2d(ii,j).gt.-10. .and. brr2d(ii,j).lt.50. ) then
c             tempo(i,j) = wta*arr2d(ii,j) + (1.-wta)*brr2d(ii,j)
c           endif
c         endif
c J. Tsai:
          if ( (.not. isnan(arr2d(ii,j))) ) then !works for real*4 and 8
            if ( arr2d(ii,j).gt.-10. .and. arr2d(ii,j).lt.50. ) then
              tempo(i,j) = arr2d(ii,j)
            endif
          endif
        enddo
      enddo

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 8080)

      if (ifanom.eq.0) then
        itaprev_reg = ita
        call scopy (nxcsm*nycsm, tempo, 1, tempo_reg, 1)
      else if (ifanom.eq.1) then
        itaprev_anom = ita
        call scopy (nxcsm*nycsm, tempo, 1, tempo_anom, 1)
      endif

c##########
      endif
c##########

c        Interpolate to ice grid

      call ocn_int (tempo, loncsm, latcsm, nxcsm, nycsm, vnull, 
     *              tempi, maskwater, ifforce)

c        Propagate under shelves and any other csm-null areas.
c        Do every call in case maskwater changes.

      call scopy (nx*ny, tempi, 1, tempocn, 1)
      call ocn_ext (tempocn, maskwater, vnull)

      if (firstcsmext) then
        call printmap (timeice, tempi,  'csm temp_orig (C)',1.,iudiag,0)
        call printmap (timeice, tempocn,'csm temp (C)',     1.,iudiag,0)
        firstcsmext = .false.
      endif
  
      return
      end

cOOOOO
#endif
cOOOOO

c-----------------------------------------------------------------------

cPPPPPPPPPPPPPPPPPPPP
#if defined (CSMORDO)
cPPPPPPPPPPPPPPPPPPPP

      subroutine csmordo_ocean_read (tempocn, maskwater, vnull, 
     *                               ifforce, timeice)

c     Reads CSM Ordovician data file for 400 m ocean temperature.
c     Global lon,lat deg grid.
c     Bilin. interpolates to ice grid (tempi), and
c     propagate under shelves and to any null regions using
c     contiguous neighbor iteration. Returns as tempocn(nx,ny).
c     Returned ogcm temperatures at a particular level (kcsm), 
c     hard-coded below (this file has only 1 level, so nzcsm=kcsm=1).

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension 
     *  tempocn(nx,ny), maskwater(nx,ny)

      character cfilin*200, cvar*16
      integer*4 ncid, rcode, varid, start(4), count(4)
      save ncid, varid, cfilin

      parameter (nxcsm=nlong, nycsm=nlatg, nzcsm=1)
      dimension arrin(nxcsm,nycsm,12)
      real*4 arr1d(10000), arr2d(nxcsm,nycsm)
      dimension loncsmd(nxcsm), latcsmd(nycsm), 
     *          loncsm(nxcsm),  latcsm(nycsm),
     *          zcsm(nzcsm),
     *          tempo(nxcsm,nycsm),
     *          tempi(nx,ny),
     *          work(1000)
      real loncsmd, latcsmd, loncsm, latcsm
      save loncsmd, latcsmd, loncsm, latcsm
      save tempo

      dimension ndaypm(12)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm

c---
c     choose ocean layer index (kcsm)  ! not used here 
c     parameter (kcsm = ...)
c     data zcsm(1) / 40000.00 /
c---

      logical firstcsmordo, firstcsmprint
      data firstcsmordo, firstcsmprint 
     *  /.true., .true./
      save firstcsmordo, firstcsmprint

c     diagnostic check output
      logical ifdiag
      data ifdiag /.true./
c     data ifdiag /.false./
      save ifdiag
      parameter (iudiag = 155)

      vnull = 9.96921e+36                         ! vnull is passed back

c===========================
      if (firstcsmordo) then
c===========================

        cfilin = cpathhome(1:lenchr(cpathhome)) // 
     *           'Data_csmordo/' // 
#  if defined (ORDO15X)
     *           'b30.153_t31_400TEMP_1870-1889.AC.nc'
#  elif defined (ORDO4X)
     *           'b30.152_t31_400TEMP_1840-1859.AC.nc'
#  endif

        write (ioterm,'(3a)') 
     *    'reading existing Netcdf ', cfilin(1:lenchr(cfilin)), 
     *    ' data file'

c          Open existing CSM Ordovician netcdf data file

        rcode = nf_open(cfilin, nf_nowrite, ncid)
        call checkrcode (rcode, cfilin, 'open', 8001)

c          ocean depth levels (0 = surface, ...):

c       start(1) = 1
c       count(1) = nzcsm
c       rcode = nf_inq_varid (ncid, 'z_t', varid)
c       call checkrcode (rcode, cfilin, 'z_t', 8045)
c       rcode = nf_get_vara_real (ncid, varid, start, count, arr1d)
c       call checkrcode (rcode, cfilin, 'z_t', 8046)
c       if (ifdiag) write (iudiag,'(/a)') 'zcsm:'
c       do k=1,nzcsm
c         zcsm(k) = arr1d(k)
c         if (ifdiag) write (iudiag,'(i6,f10.3)') k, zcsm(k)
c       enddo
c       if (ifdiag) then
c         write(iudiag,'(a/i6,f10.3)') 'kcsm :', kcsm, zcsm(kcsm)
c       endif

c          latitudes (deg N):

        start(1) = 1
        count(1) = nycsm
        start(2) = 1
        count(2) = 1

        rcode = nf_inq_varid (ncid,  'lat',  varid)
        call checkrcode (rcode, cfilin, 'lat', 8065)
        rcode = nf_get_vara_double (ncid, varid, start, count, latcsmd)
        call checkrcode (rcode, cfilin, 'lat', 8066)

c       protect interpolation logic in setinterp (needed?):
        latcsmd(1)     = max (latcsmd(1),     -89.999) 
        latcsmd(nycsm) = min (latcsmd(nycsm),  89.999) 

        if (ifdiag) write (iudiag,'(/a)') 'latcsm(d):'
        do j=1,nycsm
          latcsm(j) = latcsmd(j) * pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') j, latcsmd(j)
        enddo

c          longitudes (deg E):

        start(1) = 1
        count(1) = nxcsm
        start(2) = 1
        count(2) = 1

        rcode = nf_inq_varid (ncid,  'lon',  varid)
        call checkrcode (rcode, cfilin, 'lon', 8067)
        rcode = nf_get_vara_double (ncid, varid, start, count, work)
        call checkrcode (rcode, cfilin, 'lon', 8068)
        do i=1,nxcsm
c         shift cut from greenwich to dateline for ice model:
          ii = mod (i + nxcsm/2 -1, nxcsm)  + 1
          loncsmd(i) = work(ii)
          if (loncsmd(i).gt.180.) loncsmd(i) = loncsmd(i) - 360.
        enddo

c       protect interpolation logic in setinterp (needed?):
        if (loncsmd(1).gt. 179.999 .or. loncsmd(1).lt.-179.999) 
     *    loncsmd(1) = -179.999 
        if (loncsmd(nxcsm).gt. 179.999 .or. loncsmd(nxcsm).lt.-179.999) 
     *    loncsmd(nxcsm) = 179.999

        if (ifdiag) write (iudiag,'(/a)') 'loncsm(d):'
        do i=1,nxcsm
          loncsm(i)  = loncsmd(i) * pi/180.
          if (ifdiag) write (iudiag,'(i6,f10.2)') i, loncsmd(i)
        enddo

c---

c        Read ocean temp (x,y), interp to ice grid, propagate under 
c        shelves, return in tempocn. Read for 12 months into arrin,
c        convert to annual (arr2d).

      start(1) = 1
      count(1) = nxcsm
      start(2) = 1
      count(2) = nycsm
      start(3) = 1
      count(3) = 12
c     start(4) = kcsm
c     count(4) = 1

      cvar = 'TEMP400'
      rcode = nf_inq_varid (ncid,  cvar,  varid)
      call checkrcode (rcode, cfilin, cvar, 8072)
      rcode = nf_get_vara_double (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfilin, cvar, 8073)

c        Calculate annual means

      do j=1,nycsm
        do i=1,nxcsm
          zav = 0.
          zat = 0.
          do m=1,12
            if ( (.not. isnan(arrin(i,j,m))) ) then ! works real*4 & 8
              if ( arrin(i,j,m).gt.-10. .and. arrin(i,j,m).lt.50. ) then
                zav = zav + arrin(i,j,m)*ndaypm(m)
                zat = zat + ndaypm(m)
              endif
            endif
          enddo
          if (zat.gt.0.) then
            arr2d(i,j) = zav/zat 
          else
            arr2d(i,j) = vnull
          endif
        enddo
      enddo

      do j=1,nycsm
        do i=1,nxcsm
c         shift cut from greenwich to dateline for ice model:
c         ii = mod (i + nxcsm/2 -1, nxcsm)  + 1
          ii = i ! no shift for Ordovician fields
          tempo(i,j) = vnull 
          if ( (.not. isnan(arr2d(ii,j))) ) then !works for real*4 and 8
            if ( arr2d(ii,j).gt.-10. .and. arr2d(ii,j).lt.50. ) then
              tempo(i,j) = arr2d(ii,j)
            endif
          endif
        enddo
      enddo

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 8080)


      firstcsmordo = .false.

c==========
      endif
c==========

c        Do the rest every call in case maskwater changes. 
c        (tempo is saved, tempi is not).

c        Interpolate to ice grid

      call ocn_int (tempo, loncsm, latcsm, nxcsm, nycsm, vnull, 
     *              tempi, maskwater, ifforce)

c        Propagate under shelves and any other csm-null areas.

      call scopy (nx*ny, tempi, 1, tempocn, 1)
      call ocn_ext (tempocn, maskwater, vnull)

      if (firstcsmprint) then
        call printmap (timeice, tempi,  'csm temp_orig (C)',1.,iudiag,0)
        call printmap (timeice, tempocn,'csm temp (C)',     1.,iudiag,0)
        firstcsmprint = .false.
      endif
  
      return
      end

cPPPPP
#endif
cPPPPP

c-----------------------------------------------------------------------

      subroutine ocn_int (tempo, alono, alato, nlono, nlato, vnull, 
     *                    tempi, maskwater, ifforce)

c        Interpolate Liu/Levitus/Loveclim/Rcp/Csm ogcm temps (tempo) 
c        to ice grid (tempi) where ocean (maskwater=1) and non-null.
c        ifforce = 1 to force calculation of indices (setinterp),
c        needed only for liucorr/rcpcorr/loveanom/csmanom bias 
c        correction, or time of crossing from past (liu) to future 
c        (rcp). Otherwise setinterp is called first time only.

#include <comicegrid.h>

      dimension 
     *  tempo(nlono,nlato), alono(nlono), alato(nlato), 
     *  tempi(nx,ny), maskwater(nx,ny)

      dimension
     *  indoa(nx,ny), weioa(nx,ny), indob(nx,ny), weiob(nx,ny)
      save indoa, weioa, indob, weiob

      logical firstocnint
      data firstocnint /.true./
      save firstocnint

c        Set indices for bilinear interp (indo*,weio* local, saved)
c        from ice grid into ogcm grid

      if (firstocnint .or. ifforce.eq.1) then  
        call setinterp (alono, alato, nlono, nlato,
     *                  alond, alatd, nx,    ny,
     *                  indoa, weioa, indob, weiob, 0.)
        firstocnint = .false.
      endif

c        Interpolate non-null temps (tempo) to ice model grid (tempi)
    
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then
            il = indoa(i,j)
            ir = mod (il,nlono) + 1
            jb = indob(i,j)
            jt = min (jb+1, nlato)
            wlb =     weioa(i,j) *    weiob(i,j)
            wrb = (1.-weioa(i,j))*    weiob(i,j)
            wlt =     weioa(i,j) *(1.-weiob(i,j))
            wrt = (1.-weioa(i,j))*(1.-weiob(i,j))
            zarr = 0.
            zwei = 0.
            if (tempo(il,jb).ne.vnull) then
              zarr = zarr + tempo(il,jb)*wlb
              zwei = zwei +               wlb
            endif
            if (tempo(ir,jb).ne.vnull) then
              zarr = zarr + tempo(ir,jb)*wrb
              zwei = zwei +               wrb
            endif
            if (tempo(il,jt).ne.vnull) then
              zarr = zarr + tempo(il,jt)*wlt
              zwei = zwei +               wlt
            endif
            if (tempo(ir,jt).ne.vnull) then
              zarr = zarr + tempo(ir,jt)*wrt
              zwei = zwei +               wrt
            endif

            if (zwei.gt.0.) then
              tempi(i,j) = zarr/zwei
            else
              tempi(i,j) = vnull
            endif

          else
            tempi(i,j) = vnull
          endif
        enddo   
      enddo   

      return
      end

c-----------------------------------------------------------------------

      subroutine ocn_ext (tempi, maskwater, vnull)

c      Extends tempi (Liu/Levitus/Loveclim/Rcp/Csm temps, on ice grid) 
c      into null ocean areas (vnull, mostly under shelves), by 
c      iteratively propagating contiguous neighbors.

#include <comicegrid.h>

      dimension 
     *  tempi(nx,ny), maskwater(nx,ny)

c     local:
      dimension tempsav(nx,ny)
      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook4 = 4)   ! 4 for E-W/N-S
      parameter (nlook8 = 8)   ! 8 for diagonals too

c     fr.pt.seawater, 200 m depth:
      parameter (tfreezoc = .0939 - .057*34.5 - 7.64e-4*200.)

c        Propagate "good" tempi pts to null ocean pts. 
c        Do by 8pt-neighbor averaging, and iterating till no change.

c,,,,,,,,,,,,,,,,,,,,,,,
      do iter = 1, nx*ny
c,,,,,,,,,,,,,,,,,,,,,,,
        call scopy (nx*ny, tempi, 1, tempsav, 1)
        nnew = 0
        do j=1,ny
          do i=1,nx
            if (maskwater(i,j).eq.1 .and. tempi(i,j).eq.vnull) then
              zt = 0.
              nt = 0
              do look=1,nlook8
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (tempsav(ii,jj).ne.vnull) then
                  zt = zt + tempsav(ii,jj)
                  nt = nt + 1
                endif
              enddo
              if (nt.gt.0) then
                tempi(i,j) = zt/nt
                nnew = nnew + 1
              endif
            endif
          enddo
        enddo
        if (nnew.eq.0) go to 500
c,,,,,,,,,,,,,
      enddo
      write (6,*) 'ocn_ext: should not get here (iter 500)'
      stop
  500 continue
c,,,,,,,,,,,,,

c        For closed basins (surrounded by land, not reached above) 

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1 .and. tempi(i,j).eq.vnull) then
            tempi(i,j) = tfreezoc
          endif
        enddo
      enddo

      return
      end

c---------

c-----------------------------------------------------------------------

c:::::::::::::::::::::
#if defined (LOVECLIP)
c:::::::::::::::::::::

      subroutine love_ocn_ext (tempi, maskwater, h, hw, vnull)

c      More complex version of ocn_ext, for LOVECLIP.

c      Extends tempi (Loveclim temps, on ice grid) 
c      into null ocean areas (vnull, mostly under shelves), by 
c      iteratively propagating contiguous neighbors.

#include <comicegrid.h>

      dimension 
     *  tempi(nx,ny), maskwater(nx,ny), h(nx,ny), hw(nx,ny)

c     local:
      dimension tempsav(nx,ny), rdecay(nx,ny)
     *         ,rpath(nx,ny),   rpathsav(nx,ny)     ! diagnostic for now

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook4 = 4)   ! 4 for E-W/N-S
      parameter (nlook8 = 8)   ! 8 for diagonals too

c     fr.pt.seawater, 200 m depth:
      parameter (tfreezoc = .0939 - .057*34.5 - 7.64e-4*200.)

c        Set decay rates

      do j=1,ny
        do i=1,nx
          rdecay(i,j) = 0.
c--------------------------------------
          if (maskwater(i,j).eq.1) then
c--------------------------------------
c           ocean column thickness:
            zdz = hw(i,j)
            if (zdz.lt.500.) rdecay(i,j) = rdecay(i,j) 
     *                                   + (1./100.e3) * exp (-zdz/100.)

c           constriction due to land (straight coast = no constr):
            zrad = 200.e3
            iwid = max ( 1, nint(zrad/dx(i,j)) )
            jwid = max ( 1, nint(zrad/dy(i,j)) )
            arocn = 0
            arall  = 0
            do jj = j-jwid,j+jwid
              jja = max (1, min (ny,jj))
              do ii = i-iwid,i+iwid
#if defined (WRAPAROUND)
                iia = mod (ii+nx-1, nx) + 1
#else
                iia = max (1, min (nx, ii))
#endif
                if (maskwater(iia,jja).eq.1) 
     *            arocn = arocn + darea(iia,jja)
                arall = arall + darea(iia,jja)
              enddo
            enddo
            zopen = 2.*min (0.5, arocn/arall)   ! 0 (closed) to 1 (open)
            zf = 1./ max (zopen, 0.01) - 1.   ! 100 (closed) to 0 (open)
            rdecay(i,j) = rdecay(i,j) 
     *                  + (1./zrad) * zf

c           existence of floating ice:
            zice = max (0., min (1., h(i,j)/10.)) 
            rdecay(i,j) = rdecay(i,j)
     *                  + (1./500.e3)*zice

c           rdecay(i,j) = 0.         ! no decay effect (same as ocn_ext)
c--------------
          endif
c--------------
        enddo
      enddo

c     call zero (rpath, nx*ny)                              ! diagnostic

c        Propagate "good" tempi pts to null ocean pts. 
c        Do by 8pt-neighbor averaging, and iterating till no change.

c,,,,,,,,,,,,,,,,,,,,,,,
      do iter = 1, nx*ny
c,,,,,,,,,,,,,,,,,,,,,,,
        call scopy (nx*ny, tempi, 1, tempsav, 1)
c       call scopy (nx*ny, rpath, 1, rpathsav, 1)           ! diagnostic
        nnew = 0
        do j=1,ny
          do i=1,nx
            if (maskwater(i,j).eq.1 .and. tempi(i,j).eq.vnull) then
              zt = 0.
c             zp = 0.                                       ! diagnostic
              nt = 0
              do look=1,nlook8
#if defined (WRAPAROUND)
                ii = mod (i+ish(look)+nx-1, nx) + 1
#else
                ii = max (1, min (nx, i+ish(look)))
#endif
                jj = max (1, min (ny, j+jsh(look)))
                if (maskwater(ii,jj).eq.1 .and. tempsav(ii,jj).ne.vnull)
     *            then
                  zdd = sqrt ( (dx(i,j)*(ii-i))**2 + 
     *                         (dy(i,j)*(jj-j))**2 )
                  zrdx = 0.5*(rdecay(i,j)+rdecay(ii,jj)) * zdd
                  zwei = exp (-zrdx)
                  zt = zt + tempsav(ii,jj)*zwei + tfreezoc*(1.-zwei)
c                 zp = zp + rpathsav(ii,jj) + zrdx          ! diagnostic
                  nt = nt + 1
                endif
              enddo
              if (nt.gt.0) then
                tempi(i,j) = zt/nt
c               rpath(i,j) = zp/nt                          ! diagnostic
                nnew = nnew + 1
              endif
            endif
          enddo
        enddo
        if (nnew.eq.0) go to 500
c,,,,,,,,,,,,,
      enddo
      write (6,*) 'ocn_ext: should not get here (iter 500)'
      stop
  500 continue
c,,,,,,,,,,,,,

c        For closed basins (surrounded by land, not reached above) 

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1 .and. tempi(i,j).eq.vnull) then
            tempi(i,j) = tfreezoc
          endif
        enddo
      enddo

      return
      end

c:::::
#endif
c:::::
