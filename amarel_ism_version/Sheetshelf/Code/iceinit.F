c-----------------------------------------------------------------------

      subroutine initgrid (iczoom)

c     Initializes horizontal and vertical grids

#include <comicephys.h>
#include <comicegrid.h>

c     nb: dzeta..., dzsed... have to add up to 1:

      dimension dzeta10(10)
      data dzeta10
     *  /.05, .10, .15, .15, .20, .15, .10, .05, .03, .02/

      dimension dzeta12(12)
      data dzeta12
     *  /.05, .10, .10, .15, .20, .15, .10, .05, .04, .03, .02, .01/

c     dimension dzsed3(3)
c     data dzsed3
c    *  /.10, .45, .45/

      dimension dzbed6(6)
      data dzbed6
     *  /10., 30., 100., 300., 600., 1000./

      logical firstzoom
      save firstzoom
      data firstzoom /.true./

c        Set horizontal grids

c===================
#if defined (LONLAT)
c===================
#  if defined (ZOOMGL)
      write (ioterm,*)
     *  '*** Error (initgrid): cannot use ZOOMGL with LONLAT'
      stop
#  endif

      totarea = 0.
      do j=1,ny
        do i=1,nx
          alond(i,j)= alon1 + (alon2-alon1)*(i-.5)/nx
          alatd(i,j)= alat1 + (alat2-alat1)*(j-.5)/ny
        enddo
        alat(j) = alatd(1,j)*pi/180.
        clat(j) = cos(alat(j))
        do i=1,nx
          alon(i) = alond(i,1)*pi/180.
          clon(i) = cos(alon(i))
#  if defined (EISLINE) 
          xh(i,j) = radius*dlon*(i-0.5)*clat(j)
          yh(i,j) = radius*dlat*(j-0.5)
#  else
          xh(i,j) = radius*dlon*(i-0.5-0.5*nx)*clat(j)
          yh(i,j) = radius*dlat*(j-0.5-0.5*ny)
#  endif
          dx(i,j) = radius*dlon*clat(j)
          dy(i,j) = radius*dlat
          darea(i,j) = dx(i,j)*dy(i,j)
          totarea = totarea + darea(i,j)
        enddo
      enddo

      do j=0,nyp
        alatu(j) = (alat1 + (alat2-alat1)*(j-.5)/ny)*pi/180.
        clatu(j) = cos(alatu(j))
        do i=0,nxp
          dxu(i,j) = radius * clatu(j) * dlon
          dyu(i,j) = radius * dlat
        enddo
      enddo

      do j=0,nyp
        alatv(j) = (alat1 + (alat2-alat1)*float(j)/ny)*pi/180.
        clatv(j) = cos(alatv(j))
        do i=0,nxp
          dxv(i,j) = radius * clatv(j) * dlon
          dyv(i,j) = radius * dlat
        enddo
      enddo

      do j=0,ny
        do i=0,nx
          dxc(i,j) = dxv(i,j)
          dyc(i,j) = dyv(i,j)
        enddo
      enddo

c=====================
#elif defined (STEREO)
c=====================

c*********************************************
#  if defined (EISANTA) && defined (TRANSECTA)
c*********************************************

     xoffa = 0.
     yoffa = 0.

c~~~~~~~~~~~~~~~~~~~~~~~
#    if defined (ZOOMGL)
c~~~~~~~~~~~~~~~~~~~~~~~

c        If first call, set zoom-resolution grid vars along entire
c        transect (in common, comicegrid.h). Transect vertices are 
c        hard coded in settransect.

      if (firstzoom) then
        call settransect (xhzoom, yhzoom, alondzoom, alatdzoom,
     *                    distzoom, totlen, dzoom, nzoom, 'TRANSECTA')
        firstzoom = .false.
      endif

c        Model grid has (2*iwidzoom+1) coarse-scale boxes ("bins") with
c        zoom-grid resolution. Middle one (with coarse-index iczoom) 
c        contains the grounding line. (iwidzoom is set in comicegrid.h)

c        Set xh, yh, alond, alatd, dx, dy, distlen (model grid, 1:nx)

c..............
      do j=1,ny
c..............
        i = 0

c       coarse-scale boxes upstream from g.l.:
        if (iczoom-iwidzoom-1.ge.1) then
          do ic = 1, iczoom-iwidzoom-1
            m = nint((ic+0.5)*npzoom)
            i = i + 1
            xh(i,j) = xhzoom(m)
            yh(i,j) = yhzoom(m)
            alond(i,j) = alondzoom(m)
            alatd(i,j) = alatdzoom(m)
            call ps_scale (zscalex, zscaley, xh(i,j), yh(i,j), 
     *                     0., 0., dd0)
            dx(i,j) = dx0 * zscalex
            dy(i,j) = dy0 * zscaley
            distlen(i) = distzoom(m)
          enddo
        endif

c       2*iwidzoom+1 zoomed boxes surrounding g.l.:
        do ic = max(iczoom-iwidzoom,1), min(iczoom+iwidzoom, nx0)
          do m = (ic-1)*npzoom + 1, i*npzoom
            i = i + 1
            xh(i,j) = xhzoom(m)
            yh(i,j) = yhzoom(m)
            alond(i,j) = alondzoom(m)
            alatd(i,j) = alatdzoom(m)
            call ps_scale (zscalex, zscaley, xh(i,j), yh(i,j), 
     *                     0., 0., dd0)
            dx(i,j) = dzoom * zscalex
            dy(i,j) = dy0   * zscaley
            distlen(i) = distzoom(m)
          enddo
        enddo

c       coarse-scale boxes downstream from g.l.:
        if (iczoom+iwidzoom+1.le.nx0) then
          do ic = iczoom+iwidzoom+1, nx0
            m = nint((ic+0.5)*npzoom)
            i = i + 1
            xh(i,j) = xhzoom(m)
            yh(i,j) = yhzoom(m)
            alond(i,j) = alondzoom(m)
            alatd(i,j) = alatdzoom(m)
            call ps_scale (zscalex, zscaley, xh(i,j), yh(i,j), 
     *                     0., 0., dd0)
            dx(i,j) = dx0 * zscalex
            dy(i,j) = dy0 * zscaley
            distlen(i) = distzoom(m)
          enddo
        endif
c..........
      enddo
c..........

      if (i.ne.nx) then
        write (ioterm,'(a,i6,a,i6)')
     *     '*** Error (initgrid,zoomgl): i=',i,'  nx=',nx
        stop
      endif

c        Set darea, totarea

      totarea = 0.
      do j=1,ny
        do i=1,nx
          darea(i,j) = dx(i,j)*dy(i,j)
          totarea = totarea + darea(i,j)
        enddo
      enddo

c        Set dxu,dyu, dxv,dyv, dxc,dyc

      do j=0,nyp
        jj = max (1, min (ny, j))
        jp = max (1, min (ny, j+1))
        do i=0,nxp
          ii = max (1, min (nx, i))
          ip = max (1, min (nx, i+1))
          dxu(i,j) = 0.5*(dx(ii,jj)+dx(ip,jj))
          dyu(i,j) = 0.5*(dy(ii,jj)+dy(ip,jj))
          dxv(i,j) = 0.5*(dx(ii,jj)+dx(ii,jp))
          dyv(i,j) = 0.5*(dy(ii,jj)+dy(ii,jp))
          dxc(i,j) = 0.25*( dx(ii,jj)+dx(ii,jp)
     *                     +dx(ip,jj)+dx(ip,jp) )
          dyc(i,j) = 0.25*( dy(ii,jj)+dy(ii,jp)
     *                     +dy(ip,jj)+dy(ip,jp) )
        enddo
      enddo

c~~~~~~~~
#    else
c~~~~~~~~

c        For eisanta, transecta, not zoomgl:

c        Compute xh,yh,alond,alatd at equispaced points along transect.
c        Transect vertices are hard coded in settransect.

      call settransect (xh, yh, alond, alatd, 
     *                  distlen, totlen, dx0, nx, 'TRANSECTA')

c        Set dx,dy, dxu,dyu, dxv,dyv, dxc,dyc (allowing for p.s. scale)

      totarea = 0.
      do j=1,ny
        do i=1,nx
          call ps_scale (zscalex, zscaley, xh(i,j), yh(i,j), 
     *                   0., 0., dd0)
          dx(i,j) = dx0 * zscalex
          dy(i,j) = dy0 * zscaley
          darea(i,j) = dx(i,j)*dy(i,j)
          totarea = totarea + darea(i,j)
        enddo
      enddo

      do j=0,nyp
        do i=0,nxp
          ii = max (1, min (nx, i))
          jj = max (1, min (ny, j))
          call ps_scale (zscalex, zscaley, xh(ii,jj), yh(ii,jj), 
     *                   0., 0., dd0)
          dxu(i,j) = dx0 * zscalex
          dyu(i,j) = dy0 * zscaley
          dxv(i,j) = dx0 * zscalex
          dyv(i,j) = dy0 * zscaley
        enddo
      enddo

      do j=0,ny
        do i=0,nx
          ii = max (1, min (nx, i))
          jj = max (1, min (ny, j))
          call ps_scale (zscalex, zscaley, xh(ii,jj), yh(ii,jj), 
     *                   0., 0., dd0)
          dxc(i,j) = dx0 * zscalex
          dyc(i,j) = dy0 * zscaley
        enddo
      enddo

c~~~~~~~~~
#    endif
c~~~~~~~~~

c******
#  else
c******

#    if defined (NESTING)
c        Nested grid uses same polar stereo projection plane
c        as continental grid. Domain center is a[lon,lat]nesta,
c        set in comicegrid.h. 
c 
c        Calculate offset of nested domain's center, on continental 
c        grid's projected plane (which is polar stereo, centered on SP).
c        a[lon,lat]nesta, [x,y]offa are in comicegrid.h.

      zlon = alondnest0*pi/180.
      zlat = alatdnest0*pi/180.
      call ps_xy (zx, zy, zlat, zlon, 0., 0.)
      xoffa = zx
      yoffa = zy
#    else
c     for continental (non-nested) grid, domain center is exactly at SP.
c     (so SP is at grid corners for even nx,ny, or at a grid center
c     for odd nx,ny):  
      xoffa = 0.
      yoffa = 0.
#    endif

      totarea = 0.

c----------------------------------------------
#    if defined (TEST2D) || defined (TESTMEL2D)
c----------------------------------------------
c    for mismipplus or testmel2d:

      do j=1,ny
        do i=1,nx
          xh(i,j) = (i-0.5)*dx0
        enddo
      enddo

c        Set other x-coord metrics

      do j=1,ny
        dxu(0,j) = 2.*(xh(1,j) - 0.)
        do i=1,nx-1
          dxu(i,j) = xh(i+1,j) - xh(i,j)
        enddo
        dxu(nx,j) =  2.*(xdomain - xh(nx,j))
        dxu(nx+1,j) =  dxu(nx,j)
      enddo
      do i=0,nx+1
        dxu(i,0)    = dxu(i,1)
        dxu(i,ny+1) = dxu(i,ny)
      enddo

      do j=1,ny
        dx(1,j) = 2.*(xh(1,j) - 0.) 
        do i=2,nx-1
          dx(i,j) = 0.5*(xh(i+1,j)-xh(i-1,j)) 
        enddo
        dx(nx,j) = 2.*(xdomain - xh(nx,j))
      enddo

      do j=0,ny+1
        do i=0,nx+1
          dxv(i,j) = dx ( max(1,min(nx,i)), max(1,min(ny,j)) )
        enddo
      enddo

      call scopy ((nxp+1)*(nyp+1), dxu, 1, dxc, 1)

c       Set y-coord metrics (regular spacing)

      do j=1,ny
        do i=1,nx
          yh(i,j) = (j-0.5)*dy0
          dy(i,j) = dy0
          darea(i,j) = dx(i,j)*dy(i,j)
          totarea = totarea + darea(i,j)
        enddo
      enddo

      do j=0,ny+1
        do i=0,nx+1
          dyu(i,j) = dy0
          dyv(i,j) = dy0
          dyc(i,j) = dy0
        enddo
      enddo

      call zero (alond, nx*ny)
      call zero (alatd, nx*ny)

c--------
#    else
c--------

      do j=1,ny
        do i=1,nx
#    if defined (EISLINE)
          dx(i,j) = dx0
          dy(i,j) = dy0
          if (ny.eq.1) then
            alond(i,j) =   0.
            alatd(i,j) = -90. + (i-.5)*(dx0/radius)*(180./pi)
          else if (nx.eq.1) then
            alond(i,j) =   0.
            alatd(i,j) = -90. + (j-.5)*(dy0/radius)*(180./pi)
          else
            alond(i,j) =   0. + (i-.5)*(dx0/radius)*(180./pi)
            alatd(i,j) = -90. + (j-.5)*(dy0/radius)*(180./pi)
          endif
#    else
          zx = dx0*(i-0.5-0.5*nx)
          zy = dy0*(j-0.5-0.5*ny)
          call ps_scale (zscalex, zscaley, zx, zy, xoffa, yoffa, dd0)
          dx(i,j) = dx0 * zscalex
          dy(i,j) = dy0 * zscaley
          call ps_latlon (zlat, zlon, zx, zy, xoffa, yoffa)
          alatd(i,j) = zlat*180./pi
          alond(i,j) = zlon*180./pi
#    endif
          darea(i,j) = dx(i,j)*dy(i,j)
          totarea = totarea + darea(i,j)
        enddo
      enddo
          
      do j=0,nyp
        do i=0,nxp
#    if defined (EISLINE)
          dxu(i,j) = dx0
          dyu(i,j) = dy0
#    else
          zx = dx0*(i    -0.5*nx)
          zy = dy0*(j-0.5-0.5*ny)
          call ps_scale (zscalex, zscaley, zx, zy, xoffa, yoffa, dd0)
          dxu(i,j) = dx0 * zscalex
          dyu(i,j) = dy0 * zscaley
#    endif

#    if defined (EISLINE)
          dxv(i,j) = dx0
          dyv(i,j) = dy0
#    else
          zx = dx0*(i-0.5-0.5*nx)
          zy = dy0*(j    -0.5*ny)
          call ps_scale (zscalex, zscaley, zx, zy, xoffa, yoffa, dd0)
          dxv(i,j) = dx0 * zscalex
          dyv(i,j) = dy0 * zscaley
#    endif
        enddo
      enddo

      do j=0,ny
        do i=0,nx
#    if defined (EISLINE)
          dxc(i,j) = dx0
          dyc(i,j) = dy0
#    else
          zx = dx0*(i    -0.5*nx)
          zy = dy0*(j    -0.5*ny)
          call ps_scale (zscalex, zscaley, zx, zy, xoffa, yoffa, dd0)
          dxc(i,j) = dx0 * zscalex
          dyc(i,j) = dy0 * zscaley
#    endif
        enddo
      enddo

c---------
#    endif
c---------

#    if defined (EISLINE) 
      do j=1,ny
        do i=1,nx
          if (ny.eq.1) then
            xh(i,j) = (i-0.5)*dx0
            yh(i,j) = 0.
          else if (nx.eq.1) then
            xh(i,j) = 0.
            yh(i,j) = (j-0.5)*dy0
          else
            xh(i,j) = (i-0.5-0.5*nx)*dx0
            yh(i,j) = (j-0.5-0.5*ny)*dy0
          endif
        enddo
      enddo
#    elif ! defined (TEST2D) && ! defined (TESTMEL2D)
      do j=1,ny
        do i=1,nx
          xh(i,j) = (i-0.5)*dx0 + xoffa - 0.5*nx*dx0
          yh(i,j) = (j-0.5)*dy0 + yoffa - 0.5*ny*dy0
        enddo
      enddo
#    endif

c*******
#  endif
c*******

c====
#else
c====
      write (6,*) 'Error (initgrid): must define a projection'
      stop
c=====
#endif   
c=====

c  end of lonlat, stereo blocks

c        Set vertical ice grid

      if (nlev.eq.10) then
        call scopy (nlev, dzeta10, 1, dzeta, 1)
      else if (nlev.eq.12) then
        call scopy (nlev, dzeta12, 1, dzeta, 1)
      else if (nlev.eq.30) then
        do k=1,10
          dzeta(k) = 0.1/10.
        enddo
        do k=11,20
          dzeta(k) = 0.8/10.
        enddo
        do k=21,30
          dzeta(k) = 0.1/10.
        enddo
      else
        do k=1,nlev
          dzeta(k) = 1./nlev
        enddo
c       write (ioterm,*) 'Error: unknown nlev=',nlev
c       stop
      endif

      zeta(0) = 0.
      zeta(1) = 0.5*dzeta(1)
      do k=2,nlev
        zeta(k) = zeta(k-1) + 0.5*(dzeta(k-1)+dzeta(k))
      enddo
      zeta(nlevp) = 1.
    
      do k=0,nlev
        dzetah(k) = zeta(k+1) - zeta(k)
      enddo

      zetah(0) = 0.
      do k=1,nlev
        zetah(k) = zetah(k-1) + dzeta(k)
      enddo

c        Set vertical sediment grid

c     if (nsed.eq.3) then
c       call scopy (nsed, dzsed3, 1, dzsed, 1)
c     else
        do k=1,nsed
          dzsed(k) = 1./nsed
        enddo
c     endif

      zsedm(0) = 0.
      zsed(0) = 0.
      do k=1,nsed
        zsed(k)  = zsedm(k-1) + 0.5*dzsed(k)
        zsedm(k) = zsedm(k-1) +     dzsed(k)
      enddo
      if (abs(zsedm(nsed)-1.).gt.1.e-8) then
        write(ioterm,*) '*** Error: sum(dzsed) = ',zsedm(nsed)
        stop
      else
        zsedm(nsed) = 1.
        zsed(nsedp) = 1.
      endif

c        Set vertical bedrock grid

      if (nbed.eq.6) then
        call scopy (nbed, dzbed6, 1, dzbed, 1)
      else if (nbed.eq.1) then
        dzbed(1) = 30.  !   10.,  1.,  0.3
      else
        do k=1,nbed
          dzbed(k) = 300.
        enddo
      endif

      zbedm(0) = 0.
      zbed(0) = 0.
      bedthick = 0.
      do k=1,nbed
        zbedm(k) = bedthick
        zbed(k)  = bedthick + dzbed(k)*0.5
        bedthick = bedthick + dzbed(k)
      enddo
      zbedm(nbed) = bedthick
      zbed(nbedp) = bedthick

      return
      end

c-----------------------------------------------------------------------

      subroutine initphys (h, hs, hb, hbsd, t, tracer,
     *                     sedim, tsed, wsed, 
     *                     topbed, topbedeq, deflect, 
     *                     tbed, equiload, geoflux,
     *                     hw, tw, maskwater, maskinit, hbinit, 
#if defined (SEDIMENT) && defined (FORCE34MYR) && (defined (QUARRYREAD) || defined(QUARRYWRITE))
     *                     hbinitprev, hbrebprev, hbreb,  
     *                     quarrysec, quarrysecxy,
#endif
     *                     sedpres, arcocn, 
     *                     hsbal, hbal, indlake, npoilake, nlake,
     *                     ua, va, sealev, 
     *                     timeice, ifrest)

c     Initializes prognostic model variables

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

c     passed:
      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny),
     *  hbsd(nx,ny),
     *  t(nx,ny,0:nlevp),    tracer(nx,ny,0:nlevp,ntrace),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     *  topbed(nx,ny),       topbedeq(nx,ny),    deflect(nx,ny),
     *  tbed(nx,ny,nbed),    equiload(nx,ny),    geoflux(nx,ny),
     *  hw(nx,ny),           tw(nx,ny),
     *  maskwater(nx,ny),    maskinit(nx,ny),    hbinit(nx,ny),
#if defined (SEDIMENT) && defined (FORCE34MYR) && (defined (QUARRYREAD) || defined(QUARRYWRITE))
     *  hbinitprev(nx,ny),   hbrebprev(nx,ny),   hbreb(nx,ny),
     *  quarrysec(nsedsec),  quarrysecxy(nx,ny),
#endif
     *  sedpres(nx,ny),      arcocn(nx,ny),      
     *  hsbal(nx,ny),        hbal(nx,ny),
     *  indlake(npoimax,nlakemax), npoilake(nlakemax),
     *  ua(0:nxp,0:nyp),     va(0:nxp,0:nyp)

      dimension htmp(nx,ny), hbtmp(nx,ny), masktmp(nx,ny),
     *          topbedtmp(nx,ny), topbedsav(nx,ny)

#if defined (EISLINE)      

#  if defined (LINEB) 
c     std lineb (s.l. var):
      parameter (nprof=5)
      dimension xprof(nprof), zprof(nprof)
      data xprof /   0.,          600.,  1400.,  1401.,  1600./
      data zprof / 300.,         -250.,  -250., -2000., -2000./
c     data zprof / 300.,         -240.,  -260., -2000., -2000./

c     pinning point:
c     parameter (nprof=9)
c     dimension xprof(nprof), zprof(nprof)
c     data xprof /   0, 600,  700, 990,1010, 1030, 1400, 1401, 1600/
c     data zprof / 300,-240, -800,-800,-800, -800, -800,-2000,-2000/ !  
cc    data zprof / 300,-240, -800,-800,-139, -800, -800,-2000,-2000/ ! i
cc    data zprof / 300,-240, -800,-800,-130, -800, -800,-2000,-2000/ ! j

c     reverse slope:
c     parameter (nprof=5)
c     dimension xprof(nprof), zprof(nprof)
c     data xprof /   0.,          600.,  1400.,  1401.,  1600./
c     data zprof / 300.,         -250.,  -200., -2000., -2000./

#  elif defined (LINEC)
c     parameter (nprof=6)
c     dimension xprof(nprof), zprof(nprof)
c     data xprof /   0.,  600.,   900.,  1800.,  1801.,  2000./
c     data zprof / 400.,    0., -1000., -1000., -2500., -2500./
cc    data zprof / 400.,    0.,  -500.,  -500., -2500., -2500./ !noinsed
#    if defined (FLATANT) 
       parameter (nprof=3)
       dimension xprof(nprof), zprof(nprof)
       data xprof /   0.,  2000.,   2100. /
       data zprof /   1.,     1.,   -100./
#     else
        parameter (nprof=2)
        dimension xprof(nprof), zprof(nprof)
        data xprof /   0.,  1600. /
        data zprof / 400., -1680. /
#     endif

#  elif defined (LINED)
c     parameter (nprof=4)  ! a,b
      parameter (nprof=5)  ! c
      dimension xprof(nprof), zprof(nprof)
      parameter (xdomain = nx*dx0/1.e3) 
      parameter (xrise = 0.7*xdomain)
      parameter (xrise2 = xrise + 200.)
c     data xprof /    0.,  xrise, xrise2, xdomain /  ! a,b
c     data zprof / 1800.,  -200., -2000.,  -2000. /  ! a,b

      parameter (xplain = 0.4*xdomain)
      data xprof /    0., xplain, xrise, xrise2, xdomain /
      data zprof / 1800.,   200., -200., -2000.,  -2000. /

#  elif defined (LINEE) 
      parameter (nprof=5) 
      dimension xprof(nprof), zprof(nprof)
      parameter (xdomain = nx*dx0*1.e-3)
      data xprof /  0.,      16.,    24.,     32., xdomain/ 
#    if defined (INITEE)
ca    data zprof / -320.,  -320.,  -320.,   -320.,   -320./   ! a
      data zprof / -500.,  -500.,  -500.,   -500.,   -500./ 
#    else
      data zprof / -320.,  -320.,  -320.,   -320.,   -320./   ! b,d
c     data zprof / -320.,  -320.,  -289.,   -320.,   -320./   ! c,e
#    endif

#  elif defined (LINEF) 
cc    parameter (nprof=5) 
cc    dimension xprof(nprof), zprof(nprof)
cc    data xprof /   0., 100., 400,  1400., 1500. / 
cc    data zprof / 100., 100.,-400., -400., -500. /
c     parameter (nprof=2) 
c     dimension xprof(nprof), zprof(nprof)
c     data xprof /   0.,                    1000. /    ! 0 to 1500
c     data zprof /   0.,                    -600. /    ! 100 to -500
      parameter (nprof=4) 
      dimension xprof(nprof), zprof(nprof)
c     data xprof /   0., 250.,   350.,  1000. / 
c     data zprof / 100., 100., -2000., -100.  / 
c     data xprof /   0., 250.,   350.,  1000. / 
c     data zprof / 100., 100., -100., -2000.  / 
      data xprof /   0.,  700.,  1000.,   1500.  / 
      data zprof / 100.,  100., -1200.,  -1500. / 

#  elif defined (LINEG)
      parameter (nprof=3)
      dimension xprof(nprof), zprof(nprof)
      data xprof /    0.,  800.,  1000./
      data zprof /  100.,  100.,  -500./
c     data zprof / -330., -330.,  -330./  ! 1pt

#  elif defined (LINEH)
c     from ~/Levy/fort.21_5mf_4.5, 
c     [67.75E, -64.08S] to [177.40E, -73.82S], ~Prydz to Ross,
c     at 4.70 Ma (insignif h), in /iguana/s1/pollard/Antarc/Run5mf/4.5/.
      parameter (nprof=72)
      dimension xprof(nprof),    zprof(nprof), 
     *          alonprof(nprof), alatprof(nprof)
      data xprof /
     *   20.00000,   76.56854,  133.13708,  189.70563,  246.27417,  
     *  302.84271,  342.84271,  399.41125,  455.97980,  512.54834,  
     *  569.11688,  625.68542,  682.25397,  738.82251,  795.39105,  
     *  851.95959,  891.95959,  948.52814, 1005.09668, 1061.66522, 
     * 1118.23376, 1174.80231, 1231.37085, 1287.93939, 1344.50793, 
     * 1401.07648, 1441.07648, 1497.64502, 1554.21356, 1610.78210, 
     * 1667.35065, 1723.91919, 1780.48773, 1837.05627, 1893.62482, 
     * 1950.19336, 1990.19336, 2046.76190, 2103.33044, 2159.89899, 
     * 2216.46753, 2273.03607, 2329.60461, 2386.17316, 2442.74170, 
     * 2499.31024, 2539.31024, 2595.87878, 2652.44733, 2709.01587, 
     * 2765.58441, 2822.15295, 2878.72150, 2935.29004, 2991.85858, 
     * 3048.42712, 3088.42712, 3144.99567, 3201.56421, 3258.13275, 
     * 3314.70129, 3371.26984, 3427.83838, 3484.40692, 3540.97546, 
     * 3597.54401, 3637.54401, 3694.11255, 3750.68109, 3807.24963, 
     * 3863.81818, 3920.38672 /
      data zprof /
     *-3478.48438,-3245.79688,-2954.81250,-2638.04688,-2359.51562,
     *-2122.96875,-1563.15625, -513.15625, -621.31250, -673.67188,
     * -727.09375, -694.45312, -661.75000, -206.87500,  417.50900,
     *  538.83325,  727.75238, 1033.18848, 1029.63074, 1004.90131,
     * 1063.61792, 1107.81067, 1152.08643, 1210.93994, 1272.58240,
     * 1317.61743, 1337.34338, 1364.37537, 1402.31323, 1426.95386,
     * 1354.54163, 1382.04651, 1464.26501, 1443.16650,  979.72717,
     *  764.65674,  555.71338,  832.47711, 1128.37231, 1141.41565,
     *  941.42834,  864.74878,  791.70844,  753.87646,  784.68585,
     *  809.13367,  842.67194,  817.61993,  685.10980,  388.76340,
     *  383.82883,  299.17053,  238.81046,  452.93805,  453.79053,
     *  486.21939,  754.72327,  618.54211,  956.54541, 1337.78198,
     * 1480.35486,  487.20312, -523.85938, -627.81250, -508.35938,
     * -545.14062, -452.43750, -500.62500, -470.32812, -508.12500,
     * -463.90625, -502.25000 /
      data alonprof /
     *      67.75,      68.20,      68.66,      69.15,      69.65,
     *      70.17,      71.00,      71.57,      72.15,      72.76,     
     *      73.39,      74.05,      74.74,      75.47,      76.22,     
     *      77.01,      78.06,      78.91,      79.80,      80.73,     
     *      81.70,      82.72,      83.80,      84.92,      86.10,    
     *      87.34,      88.67,      90.00,      91.40,      92.86,     
     *      94.40,      96.01,      97.70,      99.46,     101.31, 
     *     103.24,     104.83,     106.86,     108.97,     111.16,     
     *     113.43,     115.77,     118.18,     120.65,     123.18,     
     *     125.75,     127.23,     129.81,     132.40,     135.00,     
     *     137.60,     140.19,     142.77,     145.30,     147.80,    
     *     150.26,     151.11,     153.43,     155.70,     157.89,     
     *     160.02,     162.07,     164.05,     165.96,     167.80,     
     *     169.56,     169.82,     171.47,     173.05,     174.56,     
     *     176.01,     177.40 /
      data alatprof /
     *     -64.08,     -64.54,     -65.00,     -65.45,     -65.91,     
     *     -66.36,     -66.48,     -66.93,     -67.38,     -67.82,     
     *     -68.27,     -68.71,     -69.15,     -69.59,     -70.03,     
     *     -70.46,     -70.53,     -70.96,     -71.38,     -71.80,     
     *     -72.21,     -72.61,     -73.01,     -73.41,     -73.80,     
     *     -74.18,     -74.20,     -74.56,     -74.92,     -75.27,     
     *     -75.61,     -75.94,     -76.26,     -76.57,     -76.86,     
     *     -77.14,     -77.05,     -77.30,     -77.54,     -77.75,     
     *     -77.95,     -78.13,     -78.29,     -78.43,     -78.55,     
     *     -78.65,     -78.43,     -78.49,     -78.52,     -78.53,     
     *     -78.52,     -78.49,     -78.43,     -78.35,     -78.24,     
     *     -78.12,     -77.80,     -77.64,     -77.47,     -77.28,     
     *     -77.07,     -76.84,     -76.60,     -76.34,     -76.07,     
     *     -75.78,     -75.43,     -75.13,     -74.82,     -74.49,     
     *     -74.16,     -73.82 /
#  elif defined (LINEH2)
      parameter (nprof=6)
      dimension xprof(nprof),    zprof(nprof), 
     *          alonprof(nprof), alatprof(nprof)
      data xprof    /    0.,  400.,  800., 1000., 1800., 2000. /
      data zprof    /  100., 1500.,  100.,  100.,  100., -500. /
      data alonprof /    0.,    0.,    0.,    0.,    0.,    0. /
      data alatprof /    0.,    0.,    0.,    0.,    0.,    0. /
#  elif defined (LINEH3) || defined (LINEH4)
      parameter (nprof=3)
      dimension xprof(nprof),    zprof(nprof)
      data xprof    /    0.,  800.,  2000./
      data zprof    /  100.,  100., -1000./
#  elif defined (LINEM) 
c     use comicebed for linem 
#  elif defined (SLMODEL)
#    if defined (SLPROFILEA)
      parameter (nprof=2)
      dimension xprof(nprof),    zprof(nprof)
      data xprof    /    0.,  3000./
      data zprof    /  100., -800./
#    elif defined (SLPROFILEB)
      parameter (nprof=4)
      dimension xprof(nprof),    zprof(nprof)
#      if defined (SLSIZEB)
      data xprof    /         0., SLSIZEB, 2000.,  3000./
#      else
      data xprof    /         0.,   1500., 2000.,  3000./
#      endif
c     data zprof    /      -500.,  -350., -500., -2000./
c     data zprof    /     -1000.,  -350., -500., -2000./
c     data zprof    /-SLPROFILEB,  -350., -500., -2000./
      data zprof    /-SLPROFILEB,  -600., -750., -2000./
#    endif
#  endif

#endif

      dimension work(nx,ny), work2(nx,ny)

#if defined (BEDMAP1)
      dimension bathy(nx,ny), bedelev(nx,ny), groundbed(nx,ny),
     *          thicice(nx,ny), surface(nx,ny), water(nx,ny)
#elif defined (NHA) || defined (CARB) || defined (ORDO)
      dimension worksur(nx,ny),    worktop(nx,ny), 
     *          worktopbed(nx,ny), worktopbedsd(nx,ny),
     *          workthick(nx,ny),  worksurmod(nx,ny)
#endif

#if defined (NOISEHB)
      dimension gnoise(nx,ny)
#endif

#if defined (NHA)
      character*3 cyr_nha
#endif

      character*200 cpath, cname

#include <comicebed.h>

c    Set sealev just for initial settings below, may be changed in 
c    first timestep of run (icectl, subr longterm)

      sealev = 0. 
#if defined (EISLINE) && defined (FORCE34MYR)
      sealev = 0. 
#elif defined (EISLINE) && (defined (LINEC) || defined (LINED) || defined (LINEE))
      sealev = -10000.
#elif defined (EISANTA) && defined (FORCE34MYR)
c     56 m for zero inital ice (~modern AIS equiv slr, as in longterm):
      if (ifrest.eq.-1) sealev = 56.
# endif

      call zero (hbsd, nx*ny)    ! stays 0 unless SEARISEINIT or BEDMAP2
      call zero (deflect, nx*ny) ! stays 0 unless BEDROCK,not SLMODEL

c        Set ice thickness, bed elevation (h,hb)

c=============================================================
#if defined (EISLINE) && (defined (LINEC) && defined (SINBED))
c=============================================================

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.
          zx = xh(i,j) * 1.e-3
          hb(i,j) = 800. -  1.3*zx + 1000.*cos(2.*pi*zx/500.) 
        enddo
      enddo

c=========================================
#elif defined (EISLINE) && defined (LINEF_OLD)
c=========================================

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.
          if (nx.gt.ny) then
            hb(i,j) = bedprofile(xh(i,j))
          else
            hb(i,j) = bedprofile(yh(i,j))
          endif
        enddo
      enddo

c=========================================
#elif defined (EISLINE) && defined (LINEM)
c=========================================

      do j=1,ny
        do i=1,nx
          if (nx.gt.ny) then
            hb(i,j) = bedprofile(xh(i,j))
          else
            hb(i,j) = bedprofile(yh(i,j))
          endif
          h(i,j) = 10.
        enddo
      enddo

c======================
#elif defined (EISLINE)
c======================

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.
          if (ny.eq.1) then
            zx = xh(i,j) * 1.e-3
          else if (nx.eq.1) then
            zx = yh(i,j) * 1.e-3
          else
c#  if defined (LINEH4)
            zx = sqrt (xh(i,j)**2 + yh(i,j)**2) * 1.e-3
c#  endif
          endif
          if (zx.lt.xprof(1)) then 
            hb(i,j) = zprof(1) 
          else if (zx.gt.xprof(nprof)) then 
            hb(i,j) = zprof(nprof) 
          else
            do n=1,nprof-1         
              if (zx.ge.xprof(n) .and. zx.le.xprof(n+1)) then
                zwei = (xprof(n+1)-zx) / (xprof(n+1)-xprof(n))
                hb(i,j) = zwei*zprof(n) + (1.-zwei)*zprof(n+1)
                goto 10
              endif
            enddo
   10       continue
          endif
#  if defined (LINEH)
          if (zx.lt.xprof(1)) then 
            alond(i,j) = alonprof(1)
            alatd(i,j) = alatprof(1)
          else if (zx.gt.xprof(nprof)) then 
            alond(i,j) = alonprof(nprof)
            alatd(i,j) = alatprof(nprof)
          else
            do n=1,nprof-1         
              if (zx.ge.xprof(n) .and. zx.le.xprof(n+1)) then
                zwei = (xprof(n+1)-zx) / (xprof(n+1)-xprof(n))
                alond(i,j) = zwei*alonprof(n) + (1.-zwei)*alonprof(n+1)
                alatd(i,j) = zwei*alatprof(n) + (1.-zwei)*alatprof(n+1)
                goto 15
              endif
            enddo
   15       continue
          endif
#  elif defined (LINEH2)
          alond(i,j) = 0.
          alatd(i,j) = -90. + 20.*abs(xh(i,j))/(nx*dd0)
#  endif
        enddo
      enddo

c===========================================
#elif defined (EISANTA) && defined (BEDMAP1)
c===========================================

c     Bedmap(v1) Antarctica:

c        Read Bedmap v1 files, aggregate to model grid 
c        (assumes all have same nulval).

#if defined (TRANSECTA) || defined (NESTING) 
      call bedmap1_interp ('bathy',     bathy,     nulval)
      call bedmap1_interp ('bedelev',   bedelev,   nulval)
      call bedmap1_interp ('groundbed', groundbed, nulval)
      call bedmap1_interp ('icethic',   thicice,   nulval)
      call bedmap1_interp ('surface',   surface,   nulval)
      call bedmap1_interp ('water',     water,     nulval)
#else
      call bedmap1_aggreg ('bathy',     bathy,     nulval)
      call bedmap1_aggreg ('bedelev',   bedelev,   nulval)
      call bedmap1_aggreg ('groundbed', groundbed, nulval)  
      call bedmap1_aggreg ('icethic',   thicice,   nulval)  
      call bedmap1_aggreg ('surface',   surface,   nulval)  
      call bedmap1_aggreg ('water',     water,     nulval)  
#endif

      do j=1,ny
        do i=1,nx
c         set bed elevations everywhere:
          if (nint(bedelev(i,j)).eq.nulval) then
            hb(i,j) = -4000.                        ! lats > -60
          else
            hb(i,j) = bedelev(i,j)
          endif

c         set ice thicknesses everywhere:
          if (nint(thicice(i,j)).eq.nulval) then
            h(i,j) = 0.
          else
            h(i,j) = thicice(i,j)
          endif
       
c         set hb, allowing for current local equilibrated 
c         isostatic depression by grounded ice:
          if ( alatd(i,j).gt.-60.5        .or.             
     *         nint(bathy(i,j)).ne.nulval .or.
     *         nint(water(i,j)).ne.nulval ) then
c           bathymetry:
            hb(i,j) = min (sealev-2., hb(i,j))            ! ocn bed < sl
          endif
        enddo
      enddo

c===========================================
#elif defined (EISANTA) && defined (BEDMAP2)
c===========================================

c        Read Bedmap2 netcdf file, aggregate to model grid 

      call bedmap2_aggreg (h,       'thk')
      call bedmap2_aggreg (hb,      'topg')
      call bedmap2_aggreg (hbsd,    'topgsd')
c     call bedmap2_aggreg (work,    'usurf')   ! not used

c===============================================
#elif defined (EISANTA) && defined (SEARISEINIT)
c===============================================

c        Read Searise netcdf file, aggregate to model grid 

      call searise_aggreg (h, idum, 'thk')
      call searise_aggreg (hb, idum, 'topg')
      call searise_aggreg (hbsd, idum, 'topgsd')

#if (! defined (NOLEBROCQ)) && defined (LEBROCQ2)
      call searise_aggreg2 (work,  idum, 'usrf')
      call searise_aggreg2 (work2, idum, 'lsrf2')
      call searise_aggreg2 (hb, idum, 'topg2')
      do j=1,ny
        do i=1,nx
          h(i,j) = work(i,j) - work2(i,j) 
        enddo
      enddo
#endif

c============================================
#elif defined (EISANTA) && defined (EISMINT2)
c============================================

c     EISMINT II Antarctica: 
 
      call readeisant ('BED', hb) 
c     call readeisant ('SUR', hs)           ! hs set in adjustpres below
      call readeisant ('THI', h) 
      call readeisant ('MAS', work) 
      do j=1,ny
        do i=1,nx
          if (nint(work(i,j)).eq.0) then
            h(i,j) = 0.
            hb(i,j) = -1000.
          else
c           h(i,j) = 0.     ! no initial ice
            hb(i,j) = sealev + 1.
          endif
        enddo
      enddo

c============================================
#elif defined (EISANTA) && defined (SHALDRIL)
c============================================

#  if defined (SIMPSYMINIT)
      do j=1,ny
        do i=1,nx
          zlat = abs(alatd(i,j))
          hb(i,j) = 2500. * ( 1. - min (2.,(90.-zlat)/(90.-75.)) )
        enddo
      enddo
#  else
c     Wilson+Luyendyk early Oligocene: 
      call shaldril_eoc (hb)
#  endif
      call zero (h, nx*ny)
#  if defined (FORCE34MYR)
      call bedmap2_aggreg (hbsd, 'topgsd')     ! use stddev from BEDMAP2
#  endif

c======================
#elif defined (EISANTA)
c======================

      write(ioterm,'(a/a)') 
     *  '*** Error: unknown -D option with -DEISANTA (in makefile):',
     *  '    Must be BEDMAP[1,2], SEARISEINIT, EISMINT2, or SHALDRIL'
      stop

c==================
#elif defined (NHA)
c==================

c........................
#  if defined (GREENLAND)
c........................

      call readbam (workthick,    'IceThickness')
      call readbam (worktopbed,   'BedrockElevation')
      call readbam (worktopbedsd, 'BedrockElevation_sd')
      call readbam (worksur,      'LandMask')
c     call readbam (worktop,      'SurfaceElevation')

      do j=1,ny
        do i=1,nx
          h(i,j)    = workthick(i,j)
          hb(i,j)   = worktopbed(i,j)
          hbsd(i,j) = worktopbedsd(i,j)

c         sedpres is set below by no-ice rebound calc:  
c         LandMask is 0=ocn, 1=land, 2=sheet, 3=non-Gr-land, 4=shelf:
c         if (nint(worksur(i,j)).eq.0 .or. 
c    *        nint(worksur(i,j)).eq.4) then
cc        if (nint(worktopbed(i,j)).lt.sealev) then
c           sedpres(i,j) = 1
c         else
c           sedpres(i,j) = 0
c         endif

        enddo
      enddo

c........................
#  elif defined (GLACBAY)
c........................

      call readglacbay (work,  'thick')
      call readglacbay (work2, 'topog')

      do j=1,ny
        do i=1,nx

c         h(i,j)  = work(i,j)
c         hb(i,j) = work2(i,j) - work(i,j)

c         with lia thicknesses:
c         if (work2(i,j).gt.1.) then
c           hb(i,j) = max (10., work2(i,j) - 0.3*work(i,j))
c           h(i,j) = 0.3*work(i,j)
c         else 
c           hb(i,j) = min (-100., work2(i,j)- 0.5*work(i,j)) 
c           h(i,j) = 0.
c         endif

          if (work2(i,j).gt.0.1) then
            hb(i,j) = work2(i,j) - work(i,j)
            h(i,j) = work(i,j)
          else 
            hb(i,j) = work2(i,j)
            h(i,j) = 0.
          endif

        enddo
      enddo

c......
#  else
c......

c     Northern Hemisphere, not Greenland or Glacier Bay

c     ICE5G:
c     cyr_nha = '21k'
c     cyr_nha = '15k'
      cyr_nha = '0k'

      call readnha (worksur,    'SURT',    cyr_nha,  360, 180)
      call readnha (worktop,    'TOPO',    cyr_nha,  360, 180)
      call readnha (worktopbed, 'TOPOBED', cyr_nha,  360, 180)
      call readnha (worksurmod, 'SURTSED', '0k',  360, 180)

      do j=1,ny
        do i=1,nx
          h(i,j)  = 0.
          hb(i,j) = worktopbed(i,j)
          if (nint(worksur(i,j)).eq.2) then
            h(i,j) = max (worktop(i,j) - worktopbed(i,j), 0.)
          else if (nint(worksur(i,j)).eq.0) then
            hb(i,j) = min (worktopbed(i,j), sealev - 1.)
          endif

c         sedpres set below by no-ice rebound calc:  
c         if (nint(worksurmod(i,j)).eq.0) then
c           sedpres(i,j) = 1
c         else
c           sedpres(i,j) = 0
c         endif
        enddo
      enddo

c.......
#  endif
c.......

c===================
#elif defined (CARB) 
c===================

c     global 

      call readcarb (worksur, 'SURT')
      call readcarb (worktop, 'TOPO')

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.
          if (nint(worksur(i,j)).eq.0) then
            hb(i,j) = -1000.
          else
            hb(i,j) = max (worktop(i,j), sealev + 1.)
          endif
        enddo
      enddo

c===================
#elif defined (ORDO) 
c===================

c     global 

      call readordo (iunh, worksur, 'LANDFRAC', timeice)
      call readordo (iunh, worktop, 'PHIS', timeice)

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.
          if (worksur(i,j).le.0.5) then
            hb(i,j) = -500.
          else
            hb(i,j) = max (worktop(i,j), sealev + 1.)
          endif
        enddo
      enddo

c=====================
#elif defined (TEST2D)
c=====================
c     for mismipplus:

      do j=1,ny
        do i=1,nx
          h(i,j) = 100.

          zx = xh(i,j)/300.e3
          zbx = -150. - 728.8*(zx**2) + 343.91*(zx**4) - 50.57*(zx**6)
          zy = yh(i,j) - 0.5*80.e3
          zby = 500. / (1.+ exp(-2.*(zy-24.e3)/4.e3))
     *        + 500. / (1.+ exp( 2.*(zy+24.e3)/4.e3))
          hb(i,j)  = max (zbx + zby, -720.)
        enddo
      enddo

c========================
#elif defined (TESTMEL2D)
c========================

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.
          hb(i,j) = 500. - 1000.*min (1., xh(i,j)/500.e3)
c         if (hb(i,j).lt.0.) hb(i,j) = - 2000.
#  if defined (TESTMELPROF)
          hb(i,j) = -1000.
#  endif
        enddo
      enddo

c====
#else
c====

      write(ioterm,'(a/4x,2a)') 
     *  '*** Error: unknown -D... domain (in makefile):',
     *  'Must be one of EISLINE, EISANTA, NHA, CARB, ORDO, TEST2D[AXI],'
     * ,' TESTMEL[2D,1D]'
      stop

c=====
#endif
c=====

#if defined (SEDIMENT)
c.........................
#  if defined (FORCE34MYR) 
c.........................
#    if defined (QUARRYREAD)
c     Read in new initial topog and coeff info (written by writequarry 
c     at end of prev run if QUARRYWRITE) 
#      if defined (QUARRYNOHB)
      call readquarry (hbinitprev, hbrebprev, work, quarrysec, timeice)
#      else
      call readquarry (hbinitprev, hbrebprev, hb,   quarrysec, timeice)
#      endif
#    else
#     if defined (SEDQUARRY) 
      zq = SEDQUARRY
#     else
      zq = 0.8e-10
#     endif
      call resetr (quarrysec, nsedsec, zq)
#    endif

c     set Wilson (P3, 2012) basin info (in comicegrid.h: csedsec, 
c     alonsedsec, totsedsecobs, msedsec):

      call calcwilsonbasins

c     transfer coeffs from indiv sector values to 2-D map, then smooth:
      do j=1,ny
        do i=1,nx
          quarrysecxy(i,j) = quarrysec(msedsec(i,j))
        enddo
      enddo
      call smooth9pt (quarrysecxy, nx, ny)
c......
#  else
c......
#     if defined (SEDQUARRY) 
      zq = SEDQUARRY
#     else
      zq = 0.8e-10
#     endif
      call resetr (quarrysecxy, nx*ny, zq)
c.......
#  endif
c.......
#endif

c        Modifications to Antarctic topog 

c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
#if defined (EISANTA) && defined (AGASEA)
c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

c     override h,hb with Agasea data (nb: already in SEARISEINIT
c     unless NOLEBROCQ, and BEDMAP2):

      call agasea_aggreg (h, hb)

c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
#elif defined (EISANTA) && defined (HOLSCHUH)
c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

c     for Nick Holschuh's bed topography map (all Antarctic, with 
c     his changes in Marie Byrd Land vicinity, else bedmap2)

      call bedmap2_aggreg (hb, 'z')

c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
#elif defined (EISANTA) && defined (ADJUSTPIGTHW)
c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

c     elevate PIG/THW Canisteo Peninsula bed (else no grounded ice 
c     with ALBMAP, perhaps ok for BEDMAP2(?)):

      do j=1,ny
        do i=1,nx
          zxh = xh(i,j)*1.e-3
          zyh = yh(i,j)*1.e-3
          if (zyh.gt. -374. .and. zyh.lt. -335. .and. 
     *        zxh.gt.-1743. .and. zxh.lt.-1714.) then
            zdhb = max (0., -100. - hb(i,j))
            hb(i,j) = hb(i,j) + zdhb
            h(i,j) = max (0., h(i,j) - zdhb)
          endif
        enddo
      enddo

c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
#elif defined (EISANTA) && defined (WANTBATHY)
c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

c     weight topog (modern rebounded) with W-L "Shaldril" E-O topog.
c     Use with SEARISEINIT (unless NOLEBROCQ), ZEROICEINIT,ZEROICEBEDEQ,
c     SHALP311[MIN,MAX,AVE], WANTBATHY. Linear weighting factor is given
c     in WANTBATHY (WANTBATHY*W-L + (1-WANTBATHY)*mod_reb)

      call shaldril_eoc (work)
      do j=1,ny
        do i=1,nx
          hb(i,j) = WANTBATHY*work(i,j) + (1.-WANTBATHY)*hb(i,j)
          h(i,j) = 0.
        enddo
      enddo

c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
#elif defined (EISANTA) && defined (DYNAMTOPOG)
c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

c     Change in (past) bedrock topog due to dynamic topography

      call dyntop_read (work)
      do j=1,ny
        do i=1,nx
          hb(i,j) = hb(i,j) - work(i,j)
        enddo
      enddo

c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
#elif defined (EISANTA) && defined (CLEMCLIM) && defined (CLEMMASSCON)
c}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}

c     Override by mass-conserving dataset in a subdomain within the 
c     Antarctic Peninsula (Clemens Schannwell)

      cname = 'Data_clemens/' //
     *        'MassConInput.nc'
      cname = cpathhome(1:lenchr(cpathhome)) // cname
      call clem_masscon (iuto, cname, 'topob', hb)
      call clem_masscon (iuto, cname, 'thk', h)

c}}}}}
#endif
c}}}}}

#if defined (GREENLAND) && defined (STRAIGHTJAKOB) 
        do j=1,ny
          do i=1,nx
c........................................
            if (STRAIGHTJAKOB.eq.1) then 
c........................................

c             only deepen in existing fjord (preserve geometry):
c             if (alatd(i,j).gt.69.09 .and. alatd(i,j).lt.69.33 .and.
c    *            alond(i,j).gt.-51.2 .and. alond(i,j).lt.-49.5) then
              if (alatd(i,j).gt.69.10 .and. alatd(i,j).lt.69.24 .and.
     *            alond(i,j).gt.-51.2 .and. alond(i,j).lt.-49.7) then
                if (hb(i,j).lt.0.) then
                  h(i,j) = 0.
                  hb(i,j) = -1000.
                endif   
              endif

c............................................
            else if (STRAIGHTJAKOB.eq.2) then 
c............................................

c             deepen everywhere in rectangle:
              if (alatd(i,j).gt.69.13 .and. alatd(i,j).lt.69.21 .and.
     *            alond(i,j).gt.-51.2 .and. alond(i,j).lt.-49.70) then
                h(i,j) = 0.
                hb(i,j) = -1000.
              endif

c             change fjord to land outside rectange:
              if ( (alatd(i,j).gt.69.09 .and. alatd(i,j).le.69.13 .and.
     *             (alond(i,j).gt.-51.0 .and. alond(i,j).lt.-49.70)) 
     *             .or.
     *             (alatd(i,j).ge.69.21 .and. alatd(i,j).lt.69.33 .and.
     *             (alond(i,j).gt.-50.9 .and. alond(i,j).lt.-49.70)) 
     *           ) then  
c               h(i,j) = 0.
                if (hb(i,j).lt.0.) h(i,j) = 0.
                hb(i,j) = max (hb(i,j), 10.)
              endif

c................
            endif
c................
          enddo
        enddo
#endif

c        Set initial sediment amount

#  if defined (FORCE34MYR) && defined (SEDIMENT)
      call findwater (masktmp, indlake, npoilake, nlake,
     *                h, hb, hw, sealev, timeice)
#  endif

      do j=1,ny
        do i=1,nx
          sedim(i,j)   = 0.
#if defined (SEDIMENT)
#  if defined (FORCE34MYR)
          if (masktmp(i,j).eq.0) sedim(i,j) = 0.
#  else
          sedim(i,j)   = 50.
#  endif
#endif
        enddo
      enddo

c#if defined (FORCE34MYR) && defined (SHALDRIL)
c     Fill Lambert graben w. sediment pre 34 Ma (Wilson et al.,P3,2012):
c     (actually to sealev-1, so sedpres will be 1, crhmel = sed value,
c     in basecoef_init):
c     do j=1,ny
c       do i=1,nx
c         if ( (alatd(i,j).gt.-77. .and. alatd(i,j).lt.-68.) .and. 
c    *         (alond(i,j).gt. 60. .and. alond(i,j).lt. 80.) 
c    *       ) then
c           sedim(i,j) = max (sealev-1. - hb(i,j), 0.)
c           hb(i,j) = max (sealev-1., hb(i,j))
c         endif
c       enddo
c     enddo
c#endif

c        At this point, h and hb have been set from external datasets.
c        Set hw, hs, maskwater, *lake, given sealev.
c        Then save hsbal, hbal, maskbal, for inverse (CRHINV) calcs,
c        usually modern observed (hsbal and hbal ret, maskbal in common)

c     call initwater (h, hb, hw, sealev) !set initial hw, no longer used

      call findwater (maskwater, indlake, npoilake, nlake,
     *                h, hb, hw, sealev, timeice)
      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

      call scopy (nx*ny, hs, 1, hsbal, 1)             ! for crhinv
      call scopy (nx*ny, h,  1, hbal,  1)             ! for crhinv
      call scopy_i (nx*ny, maskwater, 1, maskbal, 1)  ! for printmap
      call scopy_i (nx*ny, maskwater, 1, maskinit, 1) ! for readcrhmel/
                                                      ! crhshelf,fillcrh
c     hbinit is on restart file, so overwritten later if restart:
      call scopy (nx*ny, hb, 1, hbinit, 1)           ! for SLMODELITER,
                                                     ! QUARRYREAD/WRITE,
                                                     ! icetabtrough

c        Apply gaussian noise to bed topography

#if defined (NOISEHB)
      sigin = 2000                              ! raw noise std.dev (m)
      xsmooth = 200.e3                          ! smoothing box size (m)
      nsmooth = nint (xsmooth/dd0)
      call gaussnoise (gnoise, sigout, sigin, nsmooth, work, nx, ny)
      do j=1,ny
        do i=1,nx
          if (hw(i,j).gt.0.) then
            hb(i,j) = hb(i,j) + min (hw(i,j),gnoise(i,j))
            hw(i,j) = max (0., hw(i,j)-gnoise(i,j))
          else 
            hb(i,j) = hb(i,j) + gnoise(i,j)
            h(i,j)  = max (0., h(i,j)-gnoise(i,j))
          endif
        enddo
      enddo
      write (ioterm,'(a,f8.3,a,f8.3)') 
     *  'Gaussian noise added to hb.  sigout(m)=', sigout, 
     *  '  xsmooth(km)=', xsmooth*1.e-3
#endif

c        Set equilibium bedrock topography and isostatic load 
c        assuming prescribed initial state is in equililbrium

      do j=1,ny
        do i=1,nx
          topbed(i,j) = hb(i,j) - sedim(i,j)
          topbedeq(i,j) = topbed(i,j)
          equiload(i,j) = rhoice*h(i,j) + rholiq*hw(i,j)
     *                  + rhosed*sedim(i,j)
        enddo
      enddo

c       Lower equilibrium topography in Wilkes Subglacial Basin  
c       (and Aurora, dep on lon range), lowering more where
c       current hb is below sea level. Nb: Don't change (initial) h,hb.

#if defined (EISANTA) && defined (WILKESBATHY)
      do j=1,ny
        do i=1,nx
          if (alatd(i,j).gt. -82. .and. alatd(i,j).lt.-60. .and.
     *        alond(i,j).gt.  90. .and. alond(i,j).lt.160.) then
            if (maskwater(i,j).eq.0) then
c             if (hb(i,j).lt.0.) then
c               zold = topbedeq(i,j)
c               topbedeq(i,j) = -(abs(hb(i,j)/1.5e3)**2)*2.5e3
c             endif
              zdz = 200. * max (0., min (1., (400.-hb(i,j))/200.))
c             if (hb(i,j).le.0.) zdz = max(zdz,1.5*abs(hb(i,j)))
              if (hb(i,j).le.0.) zdz = max(zdz,WILKESBATHY*abs(hb(i,j)))
              topbedeq(i,j) = topbedeq(i,j) - zdz
              equiload(i,j) = rhoice*h(i,j) + rholiq*hw(i,j)
     *                      + rhosed*sedim(i,j)
            endif
          endif
        enddo
      enddo
#endif

c>>>>>>>>>>>>>>>>>>>>
#if defined (NESTING)
c>>>>>>>>>>>>>>>>>>>>

c        For nested runs, read in and interpolate required variables 
c        from nestdrive file written at end of previous continental run
c        (first one in namelist cnestlist, or default 'nestdrivein').
c        Sets variables (overriding initial settings above) either:
c           over entire domain (iflagnest=0, default), or
c           around domain edges  only (iflagnest=1).
c        For topbedeq+equiload, same choice with iflagnesteq=0,1, and 
c        also 2 (default), setting nowhere at all (keeping higher-res 
c        initial settings from above over whole domain). 
c        Nb: If ifrest = 0 or 1, all these variables (except
c        topbedeq,equiload) inside the domain will be overriden by 
c        readres, so iflagnest setting has no effect.
c
c        (h,hb,t,ua,va values around the domain edges are preserved in 
c        icedyn and icetherm, as are hb,topbed,sedim in icebed, icesed).
c        These (redundantly) and other variables are reset around edges,
c        by call to getnest at top of main timestep loop in icectl.

#if defined (IFLAGNEST)
      iflagnest = IFLAGNEST
#else
      iflagnest = 0
#endif

#if defined (IFLAGNESTEQ)
      iflagnesteq = IFLAGNESTEQ
#else
      iflagnesteq = 2
#endif

      if (iflagnest  .lt.0 .or. iflagnest  .gt.1 .or. 
     *    iflagnesteq.lt.0 .or. iflagnesteq.gt.2) then
        write (ioterm,'(/a,2i6)') 
     *    '*** Error: bad iflagnest,iflagnesteq:', iflagnest,iflagnesteq
        stop
      endif

      call getnest (h, hb, t, tbed, sedim, tsed, wsed,
     *              topbedeq, equiload, sedpres, arcocn, maskwater,
     *              ua, va, iflagnest, iflagnesteq, 0, timeice)

c        Reset topbed, maskwater, hw, etc.

      do j=1,ny
        do i=1,nx
          topbed(i,j) = hb(i,j) - sedim(i,j)
        enddo
      enddo
      call findwater (maskwater, indlake, npoilake, nlake,
     *                h, hb, hw, sealev, timeice)
      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

      call printmap (timeice, sedpres, 'sedpres', 0.1, iuout2d, 0)

c>>>>
#else
c>>>>

c        Calculations to set sedpres with all initial ice removed 
c        and bedrock rebounded. Use temporary htmp,hbtmp so as not 
c        to change h,hb, and reset maskwater,hw,hs (findwater) at end.
c        If (MODELINITREB), run bedrock model for up to 50 kyrs till
c        convergence. Else (default), simple calc w. local isos.rebound.
c
c        Optionally reset topbedeq, equiload to correspond to the 
c        calculated ice-free state (ZEROICEBEDEQ). Optionally use 
c        zero-ice and rebounded bedrock as initial state (ZEROICEINIT).
c
c        Only do for cold start (ifrest=-1). For restarts (ifrest=0,1),
c        sedpres will be read from restart file.  
c
c        Also, only do for 2-D (not EISLINE), and not NESTING.
c        For nesting, sedpres has been set from getnest (and may
c        also be reset from restart file).
c
c        One issue (Lev Tarasov, 4/10): whether to use sealev=0 here,
c        or sealev= ~77 equiv. to Gr+Ant all removed. Perhaps best 
c        to use 0 for runs where sealev always fixed at 0, and +77 
c        for runs where sealev varies (~d18O, or ~model ice vol) (?)
c        Or avoid the issue (for topbedeq, equiload) entirely by not 
c        using ZEROICEBEDEQ or ZEROICEINIT.

c#  if ! defined (GREENLAND)
      call zero (sedpres, nx*ny)
c#  endif

c::::::::::::::::::::::::::::
#  if ! defined (EISLINE)
c::::::::::::::::::::::::::::
      if (ifrest.eq.-1)  then
c::::::::::::::::::::::::::::

      do j=1,ny
        do i=1,nx
          htmp(i,j) = 0.
          hbtmp(i,j) = hb(i,j)
          topbedtmp(i,j) = hbtmp(i,j) - sedim(i,j)
        enddo
      enddo

c.....................
#    if (MODELINITREB)
c.....................

      call findwater (maskwater, indlake, npoilake, nlake,
     *                htmp, hbtmp, hw, sealev, timeice)
      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 htmp, hbtmp, hw, hs, sealev)

c        Spin up bedrock response for several 10,000 yrs with no ice
c        (including lithospheric effect). Nb: subr bedrock modifies 
c        both hbtmp and topbedtmp.

      write (6,'(a)') 'initphys: topbedeq spinup:   n     zdmax'

      zdt = 200.
      do n = 1, nint(50000./zdt)
        call scopy (nx*ny, topbedtmp, 1, topbedsav, 1)
        call bedrock (htmp, hbtmp, hw, topbedtmp, topbedeq, deflect, 
     *                sedim, equiload, sealev, timeice, zdt)
        call adjustpres (maskwater, indlake, npoilake, nlake,
     *                   htmp, hbtmp, hw, hs, sealev)

        call findwater (maskwater, indlake, npoilake, nlake,
     *                  htmp, hbtmp, hw, sealev, timeice)
        call adjustpres (maskwater, indlake, npoilake, nlake,
     *                   htmp, hbtmp, hw, hs, sealev)

        zdmax = -1.e20
        do j=1,ny
          do i=1,nx
            zdmax = max (zdmax, abs(topbedsav(i,j) - topbedtmp(i,j)))
          enddo
        enddo
        
        if (mod (n,10).eq.1) then 
          write (6,'(10x,".",$)')
        else
          write (6,'(".",$)')
        endif
        if (mod (n,10).eq.0) write (6,'(4x,i6,f10.4)') n, zdmax

        if (zdmax.lt.0.1) goto 500
      enddo
  500 continue
      write (6,*) 

c........
#    else
c........

c     simple initial rebound estimate (local isostatic asthen):

      write (6,'(a,i6,f10.4)') 'initphys: simple topbedeq estimate'

      do j=1,ny
        do i=1,nx
          if (hw(i,j).eq.0. and. h(i,j).gt.0.) then
c           allow for marine ice replaced by ocean (see notes 5/20/10):
c           but problems with inland deep basins, and nesting edges!
            zov = rhoice*h(i,j) - rholiq*max(0.,sealev-hb(i,j))
            hbtmp(i,j)= hb(i,j) + max(zov,0.)/rhoman 
            topbedtmp(i,j) = hbtmp(i,j) - sedim(i,j)
          endif
        enddo
      enddo

      call findwater (maskwater, indlake, npoilake, nlake,
     *                htmp, hbtmp, hw, sealev, timeice)
      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 htmp, hbtmp, hw, hs, sealev)

c.........
#    endif
c.........

c        Set sedpres (used in basecoef) depending on submerged or not.
c        Optionally reset topbedeq, equiload corresp. to the ice-free, 
c        rebounded, equilibrated state.

      do j=1,ny
        do i=1,nx
c#    if ! defined (GREENLAND)
          sedpres(i,j) = 0.
          if (hw(i,j).gt.0.) sedpres(i,j) = 1.
c#    endif
#    if defined (ZEROICEBEDEQ) && ! defined (MODOBS)
          topbedeq(i,j) = hbtmp(i,j) - sedim(i,j)
          equiload(i,j) = rhoice*htmp(i,j) + rholiq*hw(i,j)
     *                  + rhosed*sedim(i,j)
#    endif
        enddo
      enddo


#    if defined (ZEROICEINIT) && ! defined (MODOBS)
c     nb: don't change sea level to correspond to zeroed ice volume
c     (see comments above)    
      do j=1,ny
        do i=1,nx
          h(i,j) = htmp(i,j)
          hb(i,j) = hbtmp(i,j)
          topbed(i,j) = hb(i,j) - sedim(i,j)
        enddo
      enddo
#    endif

c        Reset maskwater, hs, hw, etc.

      call findwater (maskwater, indlake, npoilake, nlake,
     *                h, hb, hw, sealev, timeice)
      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

      call printmap (timeice, sedpres, 'sedpres', 0.1, iuout2d, 0)

c::::::::::
      endif
c::::::::::
#  endif
c::::::::::

c>>>>>
#endif
c>>>>>

c        Set geothermal heat flux (geoflux)

ccccccccccccccccccccc
#if defined (EISANTA)
ccccccccccccccccccccc

c..........................
#  if defined (NYBLADEFLUX)
c..........................
      call setgeoflux_nyblade (geoflux)

c.........................
#  elif defined (FOXMAULE)
c.........................
#    if (! defined (SEARISECLIM)) || defined (NOLEBROCQ)
      write (ioterm,'(/2a)') 
     *       '*** Error: For FOXMAULE, must use SEARISECLIM and',
     *       ' not NOLEBROCQ'
      stop
#    endif
      call searise_aggreg (geoflux, idum, 'bheatflx_fox')
c     impose max limit on geothermal heat flux:
c     do j=1,ny
c       do i=1,nx
c         geoflux(i,j) = min (geoflux(i,j), geoflux_wais)
c       enddo
c     enddo

c........................
#  elif defined (SHAPIRO)
c........................
#    if ! defined (SEARISECLIM)
      write (ioterm,'(/a)') 
     *       '*** Error: For SHAPIRO, must use SEARISECLIM'
      stop
#    endif
#    if defined (NOLEBROCQ)
      call searise_aggreg (geoflux, idum, 'bheatflx')          ! shapiro
#    else
      call searise_aggreg (geoflux, idum, 'bheatflx_shapiro')
#    endif
c     impose max limit on geothermal heat flux:
c     do j=1,ny
c       do i=1,nx
c         geoflux(i,j) = min (geoflux(i,j), geoflux_wais)
c       enddo
c     enddo

c......
#  else
c......
c     "simpgeoflux", default for eisanta):
      do j=1,ny
        do i=1,nx
c         higher geothermal flux under all WAIS:
          if ( alatd(i,j).gt. -86. .and.
     *         (alond(i,j).gt.170. .or. alond(i,j).lt.-30.) ) then
            geoflux(i,j) = geoflux_wais
          else
            geoflux(i,j) = geoflux_eais
          endif
        enddo
      enddo
c.......
#  endif
c.......

ccccccccccccccccccccccccc
#elif defined (GREENLAND)
ccccccccccccccccccccccccc

c        Hot spot at source of NEGIS, gaussian vs distance (distnegis)

      call resetr (geoflux, nx*ny, .042*31556926)      ! W/m2 to  J/m2/a

      zlatnegis = 76.
      zlonnegis = -35.
      distnegis = 150.e3
c     fluxnegis = .150 * 31556926
      fluxnegis = 0.                 ! no hot spot

c..............................
      if (fluxnegis.ne.0.) then
c..............................

      xnegis = cos(zlatnegis*pi/180.)*cos(zlonnegis*pi/180.)
      ynegis = cos(zlatnegis*pi/180.)*sin(zlonnegis*pi/180.)
      znegis = sin(zlatnegis*pi/180.)

      do j=1,ny
        do i=1,nx
          xpt = cos(alatd(i,j)*pi/180.)*cos(alond(i,j)*pi/180.)
          ypt = cos(alatd(i,j)*pi/180.)*sin(alond(i,j)*pi/180.)
          zpt = sin(alatd(i,j)*pi/180.)

c         either great circle:
          zdist = acos ( max (-1., min (1.,   xpt*xnegis
     *                                      + ypt*ynegis
     *                                      + zpt*znegis)) )
     *          * radius
c         or 3-D straight line:
c         zdist = sqrt (   (xpt-xnegis)**2
c    *                   + (ypt-ynegis)**2
c    *                   + (zpt-znegis)**2 )
c    *          * radius

          if (zdist.lt.2.*distnegis) then
            zwei = exp(-(zdist/distnegis)**2)
            geoflux(i,j) = fluxnegis*zwei + (1.-zwei)*geoflux(i,j) 
          endif
        enddo
      enddo

c..........
      endif
c..........

c     call printmap (0., geoflux, 'geoflux', .001*31556926*10, 177, 0)
c     call flush (177)

ccccc
#else
ccccc

      call resetr (geoflux, nx*ny, geoflux_unif)

cccccc
#endif
cccccc

c        Ice, sed, bed temperatures, sed liquid content

c----------------
      do j=1,ny
        do i=1,nx
c----------------

c            Initialize ice and sediment temperatures 
c            where significant ice or sed, assuming vertical linear 
c            diffusive profiles and basal ice at pressure melting point

c            Leave ice temps unchanged (from getnest) if nested run

#if ! defined (NESTING)
          dtdz = geoflux(i,j) / (2.1*31556926)         !conduc ~in vdif2
          if (h(i,j).gt.0.) then
c           don't have tsurf yet - use crude param for modern Antarctic:
            zts = tmelt + 34.46 - .00914*hs(i,j)
     *                          - .68775*abs(alatd(i,j))
            do k=0,nlevp
              t(i,j,k) = min (zts + dtdz*h(i,j)*zeta(k),
     *                        tmelt - dtmdh*h(i,j)*zeta(k))
            enddo
          else
            do k=0,nlevp
              t(i,j,k) = tmelt
            enddo
          endif
#endif

c            Sediment temperatures and liquid content of pores (wsed,
c            0-1, rest is ice, assumed saturated, set in vdif):

          do k=1,nsed
            if (sedim(i,j).gt.0.01) then
              tsed(i,j,k) = tmelt-dtmdh*h(i,j)
              wsed(i,j,k) = 1.
            else
              tsed(i,j,k) = tmelt
              wsed(i,j,k) = 0.
            endif
          enddo

c            Bedrock temperatures:

#if ! defined (NESTING)
          dtdz = geoflux(i,j) / (3.3*31556926)       !conduc as in vdif2
          do k=1,nbed 
            tbed(i,j,k) = t(i,j,nlevp) + dtdz*zbed(k)
          enddo
#endif

c           Water temperature (set to pressure melt point)

          tw(i,j) = tmelt - dtmdh*h(i,j)

c------------
        enddo
      enddo
c------------

c        True tracers (3-D temperatures are equiv. to tracer(*,*,*,1))

#if defined (CALVDAMAGE)
c     ice-shelf damage ("r" in Bassis and Ma, 2015, EPSL): 
      call zero (tracer(1,1,0,2), nx*ny*(nlevp+1))
#else
      if (ntraca.lt.ntrace) then
        do n=ntraca+1,ntrace
          call zero (tracer(1,1,0,n), nx*ny*(nlevp+1))
        enddo
      endif
#endif

      return
      end

c-----------------------------------------------------------------------

c))))))))))))))))))))
#if defined (EISANTA)
c))))))))))))))))))))

      subroutine settransect (xh, yh, alond, alatd, 
     *                        dist, totlen, dx0, nx, cflag)

c     Sets polar stereo coords (xh, yh, alond,alatd) along a 
c     piece-wise linear transect given by vertices (alonvert,alatvert).
c     First converts vertices to polar stereo coords (xvert, yvert).
c     Computes total length of transect, and divides evenly by nx.
c     Steps along the transect in even increments of totlen/nx, keeping
c     track of which segment the current point is in, linearly
c     interpolating between vertices of current segment in terms of
c     polar stero coords, setting (xh,yh) and converting that to 
c     lon,lat (alond,alatd). Also returns total transect length totlen.
c
c     If TRANSECTA, sets 1-D grid (from subr initgrid).
c     If not (2-D Antarctic run), define flowline for diagnostics 
c     (from subr icetabwais), and flowlines for trough output
c     (from subr icetabtrough).
c
c     cflag = transect choice:PIG, MURPHY, ROSS, WEDDELL, WILKES, SIPLE.
c             (if TRANSECTA, set by -DTRANSECTA to one of the above).

      dimension xh(nx), yh(nx), alond(nx), alatd(nx), dist(nx)
      character*(*) cflag

      parameter (pi = 3.14159265358979)

c+++++
c     PIG, beyond divide, down modern centerline, along paleo trough to
c     beyond continental rise (Joughin,J.Glac,2009; Kirshner,QSR,2012):
      parameter (npig=11)
      dimension alonvertpig(npig+1), alatvertpig(npig+1)
      data alatvertpig /
     *   -75.2,  -75.3,  -75.3,  -74.9,  -74.5,  -74.0,
     *   -73.6,  -72.8,  -72.4,  -72.0,  -70.8,  -68.0 /
      data alonvertpig /
     *   -80.0,  -85.0,  -97.5, -102.0, -103.8, -105.5, 
     *  -106.5, -107.2, -107.0, -106.5, -104.7, -101.2 /
c+++++
c     Mt. Murphy (-75.387, -111.262) into ASE:
      parameter (nmur=2)
      dimension alonvertmur(nmur+1), alatvertmur(nmur+1)
      data alatvertmur /
     *    -77.0,   -70.0,   -68.0 /
      data alonvertmur /
     *   -111.26, -111.26, -111.26 /
c+++++
c     Ross Embayment: 
      parameter (nros=3)
      dimension alonvertros(nros+1), alatvertros(nros+1)
      data alatvertros /
     *   -82.5,  -83.,   -78.,   -70. /
      data alonvertros /
     *   -110., -154.,   178.,   178. /
c+++++
c     Weddell Embayment:
      parameter (nwed=2)
      dimension alonvertwed(nwed+1), alatvertwed(nwed+1)
      data alatvertwed /
     *   -83.,  -82., -70. /
      data alonvertwed /
     *   -100., -70., -38. /
c+++++
c     Wilkes:
      parameter (nwil=2)
      dimension alonvertwil(nwil+1), alatvertwil(nwil+1)
      data alatvertwil /
     *   -80.0, -68.,   -65. /
      data alonvertwil /
     *    130.,  153.,  160. /
c+++++
c     Siple to Ross:
      parameter (nsip=3)
      dimension alonvertsip(nsip+1), alatvertsip(nsip+1)
      data alatvertsip /-81.0,  -81.0,  -70.0,  -45./   
      data alonvertsip / -90.,  -180.,  -170., -170./
c+++++
c     data alatvert /-65., -82.5,  -82.5,  -65./          ! weddell+ross
c     data alonvert /-30., -100.,  -170., -180./

c     data alatvert /-70., -82.5,  -82.5,  -68./          ! weddell+ross
c     data alonvert /-35., -100.,  -170., -179./

c     data alatvert /-70., -81.5,  -82.5,  -82.5,  -68./  ! weddell+ross
c     data alonvert /-54.,  -57.,  -100.,  -170., -179./

c     data alatvert /-82.5,  -82.5,  -68./                ! ross only
c     data alonvert /-100.,  -170., -179./

c     data alatvert /-65., -90.,  -65./                   ! wais to eais
c     data alonvert /-90., -90.,   100./

c     data alatvert /-73.82, -64.08 /                     ! ross-> prydz
c     data alonvert /177.40,  67.75 /
c+++++

      parameter (nsegmax=50, nsegmaxp=nsegmax+1)
      dimension 
     *  alonvert(nsegmaxp), alatvert(nsegmaxp), 
     *  xvert(nsegmaxp),    yvert(nsegmaxp),      seglen(nsegmax)

      if (cflag.eq.'PIG') then 
        nseg=npig
        call scopy (nseg+1, alonvertpig, 1, alonvert, 1)
        call scopy (nseg+1, alatvertpig, 1, alatvert, 1)
      else if (cflag.eq.'MUR') then 
        nseg=nmur
        call scopy (nseg+1, alonvertmur, 1, alonvert, 1)
        call scopy (nseg+1, alatvertmur, 1, alatvert, 1)
      else if (cflag.eq.'ROSS') then
        nseg=nros
        call scopy (nseg+1, alonvertros, 1, alonvert, 1)
        call scopy (nseg+1, alatvertros, 1, alatvert, 1)
      else if (cflag.eq.'WEDDELL') then
        nseg=nwed
        call scopy (nseg+1, alonvertwed, 1, alonvert, 1)
        call scopy (nseg+1, alatvertwed, 1, alatvert, 1)
      else if (cflag.eq.'WILKES') then
        nseg=nwil
        call scopy (nseg+1, alonvertwil, 1, alonvert, 1)
        call scopy (nseg+1, alatvertwil, 1, alatvert, 1)
      else if (cflag.eq.'SIPLE') then
        nseg=nsip
        call scopy (nseg+1, alonvertsip, 1, alonvert, 1)
        call scopy (nseg+1, alatvertsip, 1, alatvert, 1)
      else
        write (6,'(/a,a)') 
     *    '*** Error (settransect): unknown transect name ',cflag
        stop
      endif
      nsegp = nseg + 1
    
c        Convert lon,lat of vertices to x,y polar stereo coords

      do m=1,nsegp
        call ps_xy (xvert(m), yvert(m), 
     *              alatvert(m)*pi/180.,  alonvert(m)*pi/180.,
     *              0., 0.)
      enddo

c        Compute total transect length (totlen) and indiv segment
c        lengths (seglen)    

      totlen = 0.
      do m=1,nseg
        seglen(m) = sqrt (   (xvert(m+1)-xvert(m))**2 
     *                     + (yvert(m+1)-yvert(m))**2 )
        totlen = totlen + seglen(m)
      enddo

#if defined (TRANSECTA)
      if (abs((totlen/dx0) - nx). gt. 0.5) then
        write (6,'(a/a,f8.2/a,f8.2/a,f8.2/a,i8)')
     *    '*** Error (transect): need to adjust dx0 or nx for transect',
     *    '    totlen(km) =', totlen*.001,
     *    '    totlen/dx0 =', totlen/dx0,
     *    '    dx0        =', dx0,
     *    '    nx         =', nx
        stop
      endif
#endif

c        Step through whole transect in equal increments of totlen/nx,
c        interpolating between vertices of current segment

      iseg = 1
      curlena = 0.
      curlenb = curlena + seglen(iseg)
      do i = 1,nx
#if defined (TRANSECTA)
        dist(i) = (dx0)*(i-.5)
#else
        dist(i) = (totlen/nx)*(i-.5)
#endif
        if (dist(i).gt.curlenb .and. .not.(iseg.eq.nseg .and.i.eq.nx))
     *    then
c         Jumped into next segment: reset beginning and end distance
c         (except if last segment, can go past end point (by 1 grid pt))
          iseg = iseg + 1  
          if (iseg.gt.nseg) then 
            write (6,'(a,4i6)')
     *        '*** Error: iseg, nseg, i, nx =',iseg, nseg, i, nx
            stop
          endif
          curlena = curlenb
          curlenb = curlena + seglen(iseg) 
        endif
        zwei =  (dist(i)-curlena) / (curlenb-curlena)
        xh(i) = (1.-zwei)*xvert(iseg) + zwei*xvert(iseg+1)
        yh(i) = (1.-zwei)*yvert(iseg) + zwei*yvert(iseg+1)

c       convert xh,yh to lon,lat (as in initgrid, and inverse of 
c       calc above for vertices):
        call ps_latlon (zlat, zlon, xh(i), yh(i), 0., 0.)
        alatd(i) = zlat*pi/180.
        alond(i) = zlon*pi/180.
      enddo

c        Diagnostic output

c     iu = 88
c     if (iu.ne.0) then
c       write (iu,'(4x,a)') '   lonvert   latvert     xvert     yvert'
c       do m=1,nsegp
c         write (iu,'(4x,4f10.2)') 
c    *      alonvert(m), alatvert(m), xvert(m), yvert(m)
c       enddo
c
c       write (iu,'(/a,f10.2)')  'totlen=', .001*totlen
c
c       write (iu,'(/a)') '   i      lond      latd        xh        yh'
c       do i=1,nx
c         write (iu,'(i4,4f10.2)')
c    *      i, alond(i), alatd(i), .001*xh(i), .001*yh(i)
c       enddo
c       call flush(iu)
c     endif

      return
      end

c-----------------------------------------------------------------------

      subroutine setgeoflux_nyblade (geoflux)

#include <comicegrid.h>

      dimension geoflux(nx,ny)
   
c        Set Antarctic geotherm heat flux map (P,D,Nyblade, GPC, 2005)

      dimension work(nx,ny)

c     For heat flux region "a": 41 W/m2:
      parameter (nheata=10)
      dimension xheata(nheata), yheata(nheata)
      data xheata / 105,120,130,138,128,120, 93, 86, 86, 91/ 
      data yheata /  75, 70, 50, 45, 33, 31, 50, 59, 66, 71/ 
      dimension heata(nx,ny)
   
c     For heat flux region "b": 60 W/m2:
      parameter (nheatb=17)
      dimension xheatb(nheatb), yheatb(nheatb)
      data xheatb /
     * 55, 65, 68, 72, 86, 93,105,101, 95, 90, 80, 74, 60, 52, 50,46,50/
      data yheatb /
     *110, 99, 80, 69, 59, 50, 24, 21, 22, 25, 35, 46, 65,82,95,105,109/
      dimension heatb(nx,ny)
   
c     For heat flux region "c": 55 W/m2:
      parameter (nheatc=7)
      dimension xheatc(nheatc), yheatc(nheatc)
      data xheatc / 25,150,150,115,101, 70, 51/
      data yheatc /150,150,  1,  1, 21, 70,110/ 
      dimension heatc(nx,ny)

c        Convert from old (yyy_cycloggeo4) grid indices (i,j=1 to 150,
c        polar stereo, 40 km, centered on SP), to meters from SP.

      do n=1,nheata
        xheata(n) = (xheata(n)-0.5)*40.e3 - 0.5*150*40.e3
        yheata(n) = (yheata(n)-0.5)*40.e3 - 0.5*150*40.e3
      enddo
      do n=1,nheatb
        xheatb(n) = (xheatb(n)-0.5)*40.e3 - 0.5*150*40.e3
        yheatb(n) = (yheatb(n)-0.5)*40.e3 - 0.5*150*40.e3
      enddo
      do n=1,nheatc
        xheatc(n) = (xheatc(n)-0.5)*40.e3 - 0.5*150*40.e3
        yheatc(n) = (yheatc(n)-0.5)*40.e3 - 0.5*150*40.e3
      enddo

      call zero (heata, nx*ny)
      call zero (heatb, nx*ny)
      call zero (heatc, nx*ny)
      call polygon_wais (heata, xh,yh, nx, ny, xheata,yheata, nheata,1.)
      call polygon_wais (heatb, xh,yh, nx, ny, xheatb,yheatb, nheatb,1.)
      call polygon_wais (heatc, xh,yh, nx, ny, xheatc,yheatc, nheatc,1.)

c     write (177,*) 'heata:'
c     write (177,'(150i1)') ((nint(heata(i,j)),i=1,nx), j=ny,1,-1)
c     write (177,*) 'heatb:'
c     write (177,'(150i1)') ((nint(heatb(i,j)),i=1,nx), j=ny,1,-1)
c     write (177,*) 'heatc:'
c     write (177,'(150i1)') ((nint(heatc(i,j)),i=1,nx), j=ny,1,-1)

c       Set raw heat flux map

      do j=1,ny
        do i=1,nx
          if (nint(heata(i,j)).eq.1) then
            work(i,j) = 41.    ! mW/m2
          else if (nint(heatb(i,j)).eq.1) then
            work(i,j) = 60.    ! mW/m2
          else if (nint(heatc(i,j)).eq.1) then
            work(i,j) = 55.    ! mW/m2
          else
            work(i,j) = 70.    ! mW/m2 (default)
          endif
        enddo   
      enddo   
   
c       "krieg" to smooth boundaries between regions 
   
      iwid = max ( 1, nint(2.*(40.e3/dd0)) )
      znorm = sqrt(2.*(iwid**2))
      do j=1,ny
        do i=1,nx
          zval = 0.
          zwto = 0.
          do jj=max(1,j-iwid),min(ny,j+iwid)
            do ii=max(1,i-iwid),min(nx,i+iwid)
              zdist = sqrt( float((i-ii)**2  + (j-jj)**2) )
              zwei = max (0., 1.-zdist/znorm)
              zval = zval + work(ii,jj)*zwei
              zwto = zwto +             zwei
            enddo
          enddo
          geoflux(i,j) = (zval/zwto) * .001*31556926   ! mW/m2 to J/m2/a
          write (177,'(2i4,e20.8)') i,j,geoflux(i,j)
        enddo
      enddo

      call printmap (0., work, 'work', 3., 177, 0)
      call printmap (0., geoflux, 'geoflux', .001*31556926*3, 177, 0)

      return
      end

c)))))
#endif
c)))))

c-----------------------------------------------------------------------

c((((((((((((((((((((
#if defined (BEDMAP1)
c((((((((((((((((((((

      subroutine bedmap1_interp (cname, arro, nulval)

c     Reads Bedmap (v1) data file cname.asc, and bilinearly interpolates
c     to model grid given by alond(nx,ny),alatd(nx,ny) in comicegrid.
c     Also returns "no-data" value nulval, and sets arro to nulval
c     if surrounded by 4 nulval Bedmap points.

#include <comicegrid.h>

      character*(*) cname
      dimension arro(nx,ny)

      parameter (nximax=1500, nyimax=1500)
      dimension arri(nximax,nyimax)

      character*80 cfilout

c        Read Bedmap input file into local array arri

      call bedmapread (cname, arri, nximax, nyimax, 
     *                 nxi, nyi, ddi, nulval, iubed)

c        Bilinearly interpolate to model coords. First convert 
c        model lon, lat to polar stereo x,y, then find surrounding 
c        i,j's of Bedmap grid.

      stdparallel_bedm = -71.
      zlambda_bedm = 1. + cos((stdparallel_bedm + 90.)*pi/180.)
      radius_bedm = 6.371220e6

      do j=1,ny
        do i=1,nx
          zlon = alond(i,j)*pi/180.
          zlat = alatd(i,j)*pi/180.
          zcos = cos(alatd(i,j)*pi/180.) 
          zr = ( (1./zcos) - sqrt ((1./zcos)**2 - 1.) ) 
     *         * zlambda_bedm * radius_bedm
c         bedmap arrays are oriented s.t. 0 deg E is along neg x axis:
          zxi = zr * cos (pi - alond(i,j)*pi/180.)
          zyi = zr * sin (pi - alond(i,j)*pi/180.)

c         shift origin from SP to bottom left corner of Bedmap grid
          zxi = zxi + 0.5*nxi*ddi
          zyi = zyi + 0.5*nyi*ddi

c         find surrounding 4 pts, weights (ia,ia+1,weia) (ja,ja+1,weja)
          ia = nint(zxi/ddi)
          weia = ia + 0.5 - (zxi/ddi)
          ia = max (1, min (nxi, ia))
          ib = max (1, min (nxi, ia+1))

          ja = nint(zyi/ddi)
          weja = ja + 0.5 - (zyi/ddi)
          ja = max (1, min (nyi, ja))
          jb = max (1, min (nyi, ja+1))

c         blinearly  interpolate to model grid point
          zval = 0.
          zwet = 0.
          do jj=ja,jb
            do ii=ia,ib
              if (nint(arri(ii,jj)) .ne.nulval) then
                zwi = (ib-ii)*weia + (ii-ia)*(1.-weia) + (1.-(ib-ia))
                zwj = (jb-jj)*weja + (jj-ja)*(1.-weja) + (1.-(jb-ja))
                zval = zval + zwi*zwj*arri(ii,jj)
                zwet = zwet + zwi*zwj
              endif
            enddo
          enddo
          if (zwet.gt.0.) then
            arro(i,j) = zval/zwet
          else
            arro(i,j) = nulval
          endif
        enddo
      enddo

c        Write scratch output file

c----------------------------
      if (iubedout.ne.0) then
c----------------------------

      cfilout = cname(1:lenchr(cname)) // '.out'
      open (iubedout, file=cfilout, status='unknown')
      close (iubedout, status='delete')
      open (iubedout, file=cfilout, status='new', form='formatted')

      write (iubedout,'(a)') cname
      write (iubedout,'(a)') '   i   j       lon       lat      arro'
      do j=1,ny
        do i=1,nx
          write (iubedout,'(2i4,2f10.2, i10)') 
     *       i, j, alond(i,j), alatd(i,j), nint(arro(i,j))
        enddo
      enddo

      close (iubedout)

c----------
      endif
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine bedmap1_aggreg (cname, arro, nulval)

c     Reads Bedmap (v1) data file cname.asc, aggregates to model resol,
c     returns field in arro. Also returns "no-data" value nulval.   

#include <comicegrid.h>

      character*(*) cname
      dimension arro(nx,ny)

      parameter (nximax=1500, nyimax=1500)
      dimension arri(nximax,nyimax)

      parameter (nxo=nx, nyo=ny)
      dimension narro(nxo,nyo), narroall(nxo,nyo)
      character*1 carro(nxo,nyo)

      character*80 cfilout

c        Read Bedmap input file into local array arri

      call bedmapread (cname, arri, nximax, nyimax, 
     *                 nxi, nyi, ddi, nulval, iubed)

c        Aggregate to output resol. Loop over Bedmap points, calculate
c        the model grid box it's in (via zx,zy, distances from SP),
c        assign mean value of all those points if there's >= 1/2 good
c        data points out of all in the box.

      do jo=1,nyo
        do io=1,nxo
          arro(io,jo) = 0.
          narro(io,jo) = 0
          narroall(io,jo) = 0
        enddo
      enddo

c        zx[nr], zy[nr] are +/- distances from origin (SP):
c        rotate input 90 deg clockwise: so set zx to zynr, zy to -zxnr

      do j=1,nyi
        do i=1,nxi
          zxnr = ddi*(i-0.5-0.5*nxi)
          zynr = ddi*(j-0.5-0.5*nyi)
  
          zx =  zynr 
          zy = -zxnr 

          io = nint ( (zx + 0.5*dx0*nxo)/dx0 + 0.5 )
          jo = nint ( (zy + 0.5*dy0*nyo)/dy0 + 0.5 )

          if (io.ge.1.and.io.le.nxo .and. jo.ge.1.and.jo.le.nyo) then
            if (nint(arri(i,j)).ne.nulval) then
              arro(io,jo) = arro(io,jo) + arri(i,j)
              narro(io,jo) = narro(io,jo) + 1
            else
              narroall(io,jo) = narroall(io,jo) + 1
            endif
          endif
        enddo
      enddo

      do jo=1,nyo
        do io=1,nxo
          if (narro(io,jo).eq.0) then
            arro(io,jo) = nulval
          else if (narro(io,jo).ge.(narroall(io,jo)+1)/2) then
            arro(io,jo) = arro(io,jo)/narro(io,jo)
          else
            arro(io,jo) = nulval
          endif
        enddo
      enddo

c        Write scratch output file

c----------------------------
      if (iubedout.ne.0) then
c----------------------------

      zscale = 100.
      if (cname.eq.'water') zscale = 200.

      do jo=1,nyo
        do io=1,nxo
          if (nint(arro(io,jo)).eq.nulval) then
            carro(io,jo) = ' ' 
          else
            zv = arro(io,jo)/zscale
            izv = nint(zv)
            if (zv.eq.0.) then
              carro(io,jo) = '.'

            else if (izv.eq.0 .and. zv.gt.0.) then
              carro(io,jo) = '+'
            else if (izv.ge.1 .and. izv.le.26) then
              carro(io,jo) = char (ichar('A') + izv-1)
            else if (izv.gt.26) then
              carro(io,jo) = '*'

            else if (izv.eq.0 .and. zv.lt.0.) then
              carro(io,jo) = '-'
            else if (izv.le.-1 .and. izv.ge.-26) then
              carro(io,jo) = char (ichar('a') - izv-1)
            else if (izv.lt.-26) then
              carro(io,jo) = '#'
            endif
          endif
        enddo
      enddo

      cfilout = cname(1:lenchr(cname)) // '.out'
      open (iubedout, file=cfilout, status='unknown')
      close (iubedout, status='delete')
      open (iubedout, file=cfilout, status='new', form='formatted')

      write (iubedout,'(a)') cname
      write (iubedout,'(6x,200i5)') (io,io=5,nxo,5)
      do jo=nyo,1,-1
        write (iubedout,'(i4,2x,1000a1)') jo, (carro(io,jo),io=1,nxo)
      enddo

      close (iubedout)

c----------
      endif
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine bedmapread (cnamin, arri, nximax, nyimax, 
     *                       nxi, nyi, ddi, nulval, iu)

c     Basic read of one Bedmap (v1) data file (cnamin.asc) into temp 
c     array arri, and also reads and returns nxi,nyi, ddi, nulval 

#include <comicepath.h>

      character*(*) cnamin
      dimension arri(nximax,nyimax)

      character*20 cname, cdum
      character*200 cfilin

      cname = cnamin

      cfilin = cpathbed(1:lenchr(cpathbed))
     *         // cname(1:lenchr(cname)) // '.asc'

      open (iu, file=cfilin, status='old', form='formatted', err=500)

c        Read header lines, input resol, etc

      read (iu,*) cdum, nyi
      read (iu,*) cdum, nxi
      read (iu,*)
      read (iu,*)
      read (iu,*) cdum,ddi
      read (iu,*) cdum,nulval

      write (6,'(2a)')
     *  ' bedmapread: reading ', cfilin(1:lenchr(cfilin))
      write (6,'(a,2i6,f10.3,i6)') 
     *  '   nxi, nyi, ddi, nulval =', nxi,nyi,ddi,nulval

c        Read Bedmap data

      do i=1,nxi
        read (iu,*) (arri(i,j),j=1,nyi)
      enddo

      close (iu)

      return

  500 write (6,'(/2a)') '*** Error (bedmapread) opening file ', 
     *  cfilin(1:lenchr(cfilin))
      stop
      end

c(((((
#endif
c(((((

c-----------------------------------------------------------------------

c{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
#if defined (BEDMAP2) || defined (FORCE34MYR)
c{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{

      subroutine bedmap2_aggreg (arr, cvarin)

c     Reads Bedmap2 1-km modern data file 
c     for field name cvarin, aggregates to ice model array arr.
c     Field name cvarin = thk (ice thickness), topg (bed elev), 
c     topgsd (sub-grid bed elev s.d.) (or usurf, surface elev,not used).
c     For topgsd, process input bed elevs twice, to calculate std. dev. 
c     within each model grid box.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension arr(nx,ny)
      character*(*) cvarin

      dimension wei(nx,ny), arrb(nx,ny)

      character cfilin*200, cvar*16
      integer*4 ncid, rcode, varid, start(3), count(3)

      parameter (nxb2=6667)
      parameter (nyb2=6667)
      parameter (ddb2=1.e3)
      dimension xb2(nxb2), yb2(nyb2)
      save xb2, yb2
      real arrin(nxb2,nyb2)

      logical firstbedmap2
      save firstbedmap2
      data firstbedmap2 /.true./

      cvar = cvarin

c        Open existing Bedmap2 netcdf data file and read hi-res fields

      write (6,'(2a)') 
     *  'reading existing Bedmap2 data file for ',cvar

#  if defined (ROSETTAMINUS)
      if (cvar.eq.'thk') then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmap2_thickness_ndh.nc'
      else if (cvar.eq.'topg' .or. cvar.eq.'topgsd') then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmapminus.nc'
      else if (cvar.eq.'usurf') then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmap2_surface_ndh.nc'
      else if (cvar.eq.'z') then
        write (6,'(a)') 'using modified Holschuh file'
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_holschuh/bedmap2_bed.nc'
      endif
      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 6000)

#  elif defined (ROSETTAPLUS)
      if (cvar.eq.'thk') then 
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmap2_thickness_ndh.nc'
      else if (cvar.eq.'topg' .or. cvar.eq.'topgsd') then 
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmapplus.nc'
      else if (cvar.eq.'usurf') then 
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmap2_surface_ndh.nc'
      else if (cvar.eq.'z') then 
        write (6,'(a)') 'using modified Holschuh file'
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_holschuh/bedmap2_bed.nc'
      endif
      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 6000)

#  else
      if (cvar.eq.'thk') then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmap2_thickness_ndh.nc'
      else if (cvar.eq.'topg' .or. cvar.eq.'topgsd') then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmap2_bed_ndh.nc'
      else if (cvar.eq.'usurf') then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_bedmap2/bedmap2_surface_ndh.nc'
      else if (cvar.eq.'z') then
        write (6,'(a)') 'using modified Holschuh file'
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_holschuh/bedmap2_bed.nc'
      endif
      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 6000)
#  endif

c        First call, set Bedmap2 x and y coords

c===========================
      if (firstbedmap2) then
c===========================
        do i=1,nxb2
          xb2(i) = (-3333 + i - 1)*ddb2
        enddo
        do j=1,nyb2
          yb2(j) = (-3333 + j - 1)*ddb2
        enddo
        firstbedmap2 = .false.
c==========
      endif
c==========

c        Read Bedmap2 data field

      start(1) = 1
      count(1) = nxb2
      start(2) = 1
      count(2) = nyb2
      start(3) = 1
      count(3) = 1

      rcode = nf_inq_varid (ncid, 'z', varid)
      call checkrcode (rcode, cfilin, cvar, 7770)
      rcode = nf_get_vara_double (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfilin, cvar, 7771)

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfilin, 'close', 8000)

c        Tile to mode grid. Find [ita,itb], [jta,jtb], Bedmap2 index 
c        ranges that are within each model h-cell (like searise_aggreg).
c        Searise grid has even # of h-grid cells (SP at 4 corners).
c        Bedmap2 edges are on model centers and vice versa, acccounted
c        for by weights of 1, 1/2 or 1/4.

c        For non zoomgl applics, should use >= ~5 km grid, and 
c        this is like tiling (>> 1 km BEDMAP2 resol). For zoomgl,
c        should use coarse grid = 1 km, and this is like 
c        nearest-neighbor; coarse grid will capture all info available
c        in Bedmap2 (and invariant topbedeq, etc, preserved by *sav 
c        in rezoomgl). First call sets iczoom = nx0, so zoomed sector
c        is initially in ocean (hopefully) and out of the way.      

      xcorner = xb2(1) - 0.5*ddb2
      ycorner = yb2(1) - 0.5*ddb2

      do j=1,ny
        do i=1,nx
          jta = int((yh(i,j)-0.5*dy0-ycorner)/ddb2 + .001) + 1
          jtb = int((yh(i,j)+0.5*dy0-ycorner)/ddb2 - .001) + 1
          jtb = max (jta,jtb)
          ita = int((xh(i,j)-0.5*dx0-xcorner)/ddb2 + .001) + 1
          itb = int((xh(i,j)+0.5*dx0-xcorner)/ddb2 - .001) + 1
          itb = max (ita,itb)

          arr(i,j) = 0.
          zweitot = 0.
          do jt = jta,jtb
            jj = max (1, min (nyb2, jt))
            do it = ita,itb
              ii = max (1, min (nxb2, it))
              zwei = 1.
              if (jt.eq.jta .or. jt.eq.jtb) zwei = zwei*0.5
              if (it.eq.ita .or. it.eq.itb) zwei = zwei*0.5
              if (.not. isnan(arrin(ii,jj))) then
                arr(i,j) = arr(i,j) + zwei*arrin(ii,jj)
                zweitot = zweitot + zwei
              endif
            enddo
          enddo
          if (zweitot.gt.0.) then 
            arr(i,j) = arr(i,j) / zweitot
          else
            arr(i,j) = -9999.
          endif

c            For std. dev. of subgrid bed topog

          if (cvar.eq.'topgsd' .and. arr(i,j).ne.-9999.) then
            zmean = arr(i,j)
            zms = 0.
            zweitot = 0.
            do jt = jta,jtb
              jj = max (1, min (nyb2, jt))
              do it = ita,itb
                ii = max (1, min (nxb2, it))
                zwei = 1.
                if (jt.eq.jta .or. jt.eq.jtb) zwei = zwei*0.5
                if (it.eq.ita .or. it.eq.itb) zwei = zwei*0.5
                if (.not. isnan(arrin(ii,jj))) then
                  zms = zms + zwei*(arrin(ii,jj)-zmean)**2
                  zweitot = zweitot + zwei
                endif
              enddo
            enddo
            if (zweitot.gt.0.) then
              arr(i,j) = sqrt(zms/zweitot)
            else
              arr(i,j) = -9999.
            endif
          endif

        enddo
      enddo

c        Reset any no-data model grid points to default values

      do j=1,ny
        do i=1,nx
          if (arr(i,j).eq.-9999.) then
            if (cvar.eq.'topg' .or. cvar.eq.'z') then
              arr(i,j) = -5000.
            else
              arr(i,j) = 0.
            endif
          endif
        enddo
      enddo

      return
      end

c{{{{{
#endif
c{{{{{

c-----------------------------------------------------------------------

c^^^^^^^^^^^^^^^^^^^^^^
#if defined (GREENLAND)
c^^^^^^^^^^^^^^^^^^^^^^

      subroutine readbam (arr, cvarin)

c     Reads Bamber et al. (2013) modern data file for Greenland 
c     (Netcdf) for field name cvarin, aggregates to ice model array arr.
c     NP polar stereo, 1 km, std parallel=71, lon.cut=-39 E.
c     Null value = -9999.
c     Field names (cvarin) are: 
c        IceThickness (ice thickness) 
c        BedrockElevation (bed elev)
c        SurfaceElevation (surface elev)
c        LandMask (mask: 0=ocn, 1=land, 2=sheet, 3=non-Gr-land, 4=shelf)
c        projection_x_coordinate (m)
c        projection_y_coordinate (m)
c     For s.d.(bedrock elev), process input bed elevs twice, to 
c     calculate std. dev. within each model grid box.

c     Bamber et al., 2013, Science (canyon).
c     Bamber, J.L. et al., 2013, A new bed elevation dataset for 
c     Greenland, The Cryo., 7, 499-510.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension arr(nx,ny)
      character*(*) cvarin

      character cfilin*200, cvar*24
      integer*4 ncid, rcode, varid, start(3), count(3)

      parameter (nxb2=2501)
      parameter (nyb2=3001)
      parameter (ddb2=1.e3)
      real xb2(nxb2), yb2(nyb2)
      dimension indb2(nx,ny), jndb2(nx,ny), arrb(nx,ny), wei(nx,ny)
      save xb2, yb2, indb2, jndb2

      real*4 arrin(nxb2,nyb2)

      logical firstbam
      save firstbam
      data firstbam /.true./

      if (cvarin.eq.'BedrockElevation_sd') then
        cvar = 'BedrockElevation'
      else
        cvar = cvarin
      endif

c        Open existing netcdf data file and read field

      write (6,'(2a)') 
     *  'reading existing Greenland (Bamber) data file for ', cvar
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *         '/Data_greenland/Greenland_bedrock_topography_V3.nc'
c    *         '/Data_greenland/Bedmachine.nc' ! Morlighem 2017 BM3
      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 6000)

c        First call, calculate x and y coords of model grid points

c=======================
      if (firstbam) then
c=======================
        start(1) = 1
        count(1) = nxb2
        start(2) = 1
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'projection_x_coordinate', varid)
        call checkrcode (rcode, cfilin, 'projection_x_coordinate', 7750)
        rcode = nf_get_vara_double (ncid, varid, start, count, xb2)
        call checkrcode (rcode, cfilin, 'projection_x_coordinate', 7751)

        start(1) = 1
        count(1) = nyb2
        start(2) = 1
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'projection_y_coordinate', varid)
        call checkrcode (rcode, cfilin, 'projection_y_coordinate', 7760)
        rcode = nf_get_vara_double (ncid, varid, start, count, yb2)
        call checkrcode (rcode, cfilin, 'projection_y_coordinate', 7761)

        do j=1,ny
          do i=1,nx
            call ps_xy (xbam, ybam, 
     *                  alatd(i,j)*pi/180., alond(i,j)*pi/180., 0., 0.)
            indb2(i,j) = nint ((xbam - xb2(1))/ddb2 + 1.)
            jndb2(i,j) = nint ((ybam - yb2(1))/ddb2 + 1.)
          enddo
        enddo

        firstbam = .false.
c==========
      endif
c==========

c        Read data field

      start(1) = 1
      count(1) = nxb2
      start(2) = 1
      count(2) = nyb2
      start(3) = 1
      count(3) = 1

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfilin, cvar, 7770)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfilin, cvar, 7771)

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfilin, 'close', 8000)

c        Aggregate to model grid. Simple tiling, ignoring any overlaps
c        of data cells spanning 2 or more model cells. 

      do j=1,ny
        do i=1,nx
          iia =  max (indb2(i,j) - nint(dx(i,j)/ddb2), 1)
          iib =  min (indb2(i,j) + nint(dx(i,j)/ddb2), nxb2)
          jja =  max (jndb2(i,j) - nint(dy(i,j)/ddb2), 1)
          jjb =  min (jndb2(i,j) + nint(dy(i,j)/ddb2), nyb2)

          za = 0.
          zwei = 0.
          do jj=jja,jjb
            do ii=iia,iib
              if (arrin(ii,jj).ne.-9999.) then 
                za = za + arrin(ii,jj)
                zwei = zwei + 1.
              endif
            enddo
          enddo
          if (zwei.gt.0.) then
            arr(i,j) = za/zwei
          else
            arr(i,j) = -9999.
          endif

          if (cvarin.eq.'BedrockElevation_sd' .and. arr(i,j).ne.-9999)
     *      then
            za2 = 0.
            zwei2 = 0.
            do jj=jja,jjb
              do ii=iia,iib
                if (arrin(ii,jj).ne.-9999.) then 
                  za2 = za2 + (arrin(ii,jj)-arr(i,j))**2
                  zwei2 = zwei2 + 1.
                endif
              enddo
            enddo
            if (zwei2.gt.0.) then
              arr(i,j) = sqrt(za2/zwei2)
            else
              arr(i,j) = -9999.
            endif
          endif

        enddo
      enddo


c        Reset any no-data model grid points to default values

      do j=1,ny
        do i=1,nx
          if (arr(i,j).eq.-9999.) then
            if (cvarin.eq.'BedrockElevation') then
              arr(i,j) = -5000.
            else
              arr(i,j) = 0.
            endif
          endif
        enddo
      enddo

      return
      end

c^^^^^
#endif
c^^^^^
c-----------------------------------------------------------------------

c++++++++++++++++++++
#if defined (GLACBAY)
c++++++++++++++++++++

      subroutine readglacbay (arr, cvarin)

c     Reads Glacier Bay modern topog/bathym file, or ice thickness
c     file (Netcdf). Former if cvarin='thick',latter if = 'topog'.
c     lon-lat grid, same as ice model, return in array arr.
c     Null value = 1.175494e-38.
c     Field name on files is Band1.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension arr(nx,ny)
      character*(*) cvarin

      character cfilin*200, cvar*24
      integer*4 ncid, rcode, varid, start(3), count(3)

      parameter (nxb=nx)
      parameter (nyb=ny)
      real*4 arrin(nxb,nyb)
      data vnull /1.175494e-38/

c        Open existing netcdf data file and read field

      write (6,'(2a)') 
     *  'reading existing Glacier Bay data file for ', cvarin
#if defined (RESOL05)
      if (cvarin.eq.'thick') then 
        cfilin = cpathhome(1:lenchr(cpathhome)) //
c    *       '/Data_glacbay/GlacierBay_LIA_thickness_0.01x0.005.nc'
     *       '/Data_glacbay/thickness_map_0.01x0.005.nc'
      else if (cvarin.eq.'topog') then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *       '/Data_glacbay/GlacierBay_NOAA-Huss_cropped_0.01x0.005.nc'
#elif defined (RESOL1)
      if (cvarin.eq.'thick') then 
        cfilin = cpathhome(1:lenchr(cpathhome)) //
c    *       '/Data_glacbay/GlacierBay_LIA_thickness_0.02x0.01.nc'
     *       '/Data_glacbay/thickness_map_0.02x0.01.nc'
      else if (cvarin.eq.'topog') then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *       '/Data_glacbay/GlacierBay_NOAA-Huss_cropped_0.02x0.01.nc'
#endif
      endif
      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 6000)

c        Read data field

      cvar = 'Band1'

      start(1) = 1
      count(1) = nxb
      start(2) = 1
      count(2) = nyb
      start(3) = 1
      count(3) = 1

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfilin, cvar, 7770)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfilin, cvar, 7771)

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfilin, 'close', 8000)

      do j=1,ny
        do i=1,nx
          arr(i,j) = arrin(i,j)
        enddo
      enddo

      return
      end

c+++++
#endif
c+++++

c-----------------------------------------------------------------------

c<<<<<<<<<<<<<<<<<<<
#if defined (AGASEA)
c<<<<<<<<<<<<<<<<<<<

      subroutine agasea_aggreg (h, hb)

c     Reads Agasea data file, aggregates to model resolution,
c     sets h, hb for those pts with data.

#include <comicegrid.h>
#include <comicepath.h>

      dimension h(nx,ny), hb(nx,ny)

      dimension hag(nx,ny), hbag(nx,ny), weiag(nx,ny)
      character*1 coutag(nx,ny)
      character*200 cfilin

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_agasea/ASE05.lonlatthkbed'

      write (6,'(/a)') 
     *  'Reading agasea data file' 
      open (iuaga, file=cfilin, status='old')

c        Aggregate to ice model grid. Loop over Bedmap points, calculate
c        the model grid box it's in (via zx,zy, distances from SP),
c        assign mean value of all those points if there's >= 1/2 good
c        data points out of all in the box.

      do j=1,ny
        do i=1,nx
          hag(i,j) = 0.
          hbag(i,j) = 0.
          weiag(i,j) = 0.
        enddo
      enddo

      stdparallel_aga = -71.
      zlambda_aga = 1. + cos((stdparallel_aga + 90.)*pi/180.)
      radius_aga = 6.371220e6

      iline = 0
c-------------
  100 continue
c-------------
      iline = iline + 1
      read (iuaga, *, end=200, err=300) alonag, alatag, hagin, hbagin 
      hagin = max (hagin, 0.)

      zlon = alonag*pi/180.
      zlat = alatag*pi/180.
      zcos = cos(zlat)
      zr = ( (1./zcos) - sqrt ((1./zcos)**2 - 1.) )
     *         * zlambda_aga * radius_aga
      zx = zr * cos (0.5*pi - zlon) - xoffa
      zy = zr * sin (0.5*pi - zlon) - yoffa

      io = nint ( (zx + 0.5*dx0*nx)/dx0 + 0.5 )
      jo = nint ( (zy + 0.5*dy0*ny)/dy0 + 0.5 )

      zwei = 1.        ! equal weights for all pts within model grid box
c     zwei = sqrt (   ((zx/dx0 + 0.5*nx) - (io-0.5))**2
c    *              + ((zy/dy0 + 0.5*ny) - (jo-0.5))**2 )     ! krieging
c     zwei = max (0., min (10., 10.*(1. - (zwei-0.1)/0.9)))

c     write (196,'(a,2f8.2, 6f8.1, 2i5)') 
c    *  'alonag, alatag, zr, zx, zy, xoffa, yoffa, dd0, io, jo=',
c    *   alonag, alatag, 
c    *   zr/1.e3, zx/1.e3, zy/1.e3, xoffa/1.e3, yoffa/1.e3, dd0/1.e3, 
c    *   io, jo

      if (io.ge.1.and.io.le.nx .and. jo.ge.1.and.jo.le.ny) then
        hag(io,jo)   = hag(io,jo)   + hagin *zwei
        hbag(io,jo)  = hbag(io,jo)  + hbagin*zwei
        weiag(io,jo) = weiag(io,jo) +        zwei
      endif
      goto 100

c-------------
  200 continue
c-------------
      close (iuaga)
      write (6,'(a,i6)') '  agasea lines=', iline-1

      do j=1,ny
        do i=1,nx
          if (weiag(i,j).eq.0.) then
            hag(i,j) = -999999.
            hbag(i,j) = -999999.
          else 
            hag(i,j) = hag(i,j)/weiag(i,j)
            hbag(i,j) = hbag(i,j)/weiag(i,j)
            h(i,j) = hag(i,j)
            hb(i,j) = hbag(i,j)
          endif
        enddo
      enddo

c        Write scratch output maps

      do j=1,ny
        do i=1,nx
          if (weiag(i,j).eq.0.) then
            coutag(i,j) = ' ' 
          else
            zv = hag(i,j)/150.
            izv = nint(zv)
            if (zv.eq.0.) then
              coutag(i,j) = '.'
            else if (izv.eq.0 .and. zv.gt.0.) then
              coutag(i,j) = '+'
            else if (izv.ge.1 .and. izv.le.26) then
              coutag(i,j) = char (ichar('A') + izv-1)
            else if (izv.gt.26) then
              coutag(i,j) = '*'
            else if (izv.eq.0 .and. zv.lt.0.) then
              coutag(i,j) = '-'
            else if (izv.le.-1 .and. izv.ge.-26) then
              coutag(i,j) = char (ichar('a') - izv-1)
            else if (izv.lt.-26) then
              coutag(i,j) = '#'
            endif
          endif
        enddo
      enddo
      write (197,'(a)') 'hag:'
      write (197,'(5x,1000i5)') (i,i=5,nx,5)
      do j=ny,1,-1
        write (197,'(i4,1x,1000a1)') j, (coutag(i,j),i=1,nx)
      enddo

      do j=1,ny
        do i=1,nx
          if (weiag(i,j).eq.0.) then
            coutag(i,j) = ' ' 
          else
            zv = hbag(i,j)/150.
            izv = nint(zv)
            if (zv.eq.0.) then
              coutag(i,j) = '.'
            else if (izv.eq.0 .and. zv.gt.0.) then
              coutag(i,j) = '+'
            else if (izv.ge.1 .and. izv.le.26) then
              coutag(i,j) = char (ichar('A') + izv-1)
            else if (izv.gt.26) then
              coutag(i,j) = '*'
            else if (izv.eq.0 .and. zv.lt.0.) then
              coutag(i,j) = '-'
            else if (izv.le.-1 .and. izv.ge.-26) then
              coutag(i,j) = char (ichar('a') - izv-1)
            else if (izv.lt.-26) then
              coutag(i,j) = '#'
            endif
          endif
        enddo
      enddo
      write (197,'(/a)') 'hbag:'
      write (197,'(5x,1000i5)') (i,i=5,nx,5)
      do j=ny,1,-1
        write (197,'(i4,1x,1000a1)') j, (coutag(i,j),i=1,nx)
      enddo

      do j=1,ny
        do i=1,nx
          if (weiag(i,j).eq.0.) then
            coutag(i,j) = ' ' 
          else
            zv = weiag(i,j)
            izv = nint(zv)
            if (zv.eq.0.) then
              coutag(i,j) = '.'
            else if (izv.eq.0 .and. zv.gt.0.) then
              coutag(i,j) = '+'
            else if (izv.ge.1 .and. izv.le.26) then
              coutag(i,j) = char (ichar('A') + izv-1)
            else if (izv.gt.26) then
              coutag(i,j) = '*'
            else if (izv.eq.0 .and. zv.lt.0.) then
              coutag(i,j) = '-'
            else if (izv.le.-1 .and. izv.ge.-26) then
              coutag(i,j) = char (ichar('a') - izv-1)
            else if (izv.lt.-26) then
              coutag(i,j) = '#'
            endif
          endif
        enddo
      enddo
      write (197,'(/a)') 'weiag:'
      write (197,'(5x,1000i5)') (i,i=5,nx,5)
      do j=ny,1,-1
        write (197,'(i4,1x,1000a1)') j, (coutag(i,j),i=1,nx)
      enddo

      return

c-------------
  300 continue
c-------------
      write (6,'(a,i6)') 'Error reading AGASEA data file, line=',iline
      stop

      end

c<<<<<
#endif
c<<<<<

c-----------------------------------------------------------------------

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if defined (SHALDRIL) || defined (WANTBATHY)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine shaldril_eoc (hb)

#include <comicegrid.h>
#include <comicepath.h>
      parameter (iushal = 60)

c     Reads Eocene topog/bathym from Doug Wilson's datafile (ascii,
c     polar stereo, 10 or 20-km, 71S true lat, WGS84 ellipsoid), and
c     tiles crudely to model bedrock elevation hb.

      dimension hb(nx,ny) 

#if defined SHALBEDMAP
c     bednew4 BEDMAP(v1):
c     10 km input x runs -3190 to 3200 km, y runs 3200 to -3190 km:
      parameter (nxi=640, nyi=640, ddin=10.e3)
#elif defined (SHALGRL09)
c     Wilson-Luyendyk, GRL 2009, Eocene-Oligocene: 
c     10 km input x runs -3190 to 3200 km, y runs 2500 to -2600 km:
      parameter (nxi=640, nyi=511, ddin=10.e3)
c     20 km input x runs -3180 to 3200 km, y runs 2500 to -2600 km:
c     parameter (nxi=320, nyi=256, ddin=20.e3)
#elif defined (SHALP311MAX) || defined (SHALP311MIN) || defined (SHALP311AVE)
c     Wilson et al, Palaeo3 2011, Eocene-Oligocene: 
      parameter (nxi=640, nyi=640, ddin=10.e3)
#elif defined (SHALFILMIN)
c     Wilson et al, Palaeo3 2011, Eocene-Oligocene:
c     with TAM troughs filled in
c     10 km input x runs -2800 to 2800 km, y runs 2800 to -2800 km:  
      parameter (nxi=561, nyi=561, ddin=10.e3)
#elif defined (SHALMIOC)
c     Wilson Early Miocene:
c     10 km input x runs -3190 to 3200 km, y runs 2500 to -2600 km:
      parameter (nxi=640, nyi=511, ddin=10.e3)
#elif defined (SHALMIOC15)
c     Gasson version of Wilson Early Miocene:
c     10 km input x runs -2800 to 3000 km, y runs 2500 to -2500 km:
      parameter (nxi=581, nyi=501, ddin=10.e3)
#elif defined (BEDMAP2RBD)
c     Wilson rebounded Bedmap2 used in Jamieson et al., 2013 
c     10 km input x runs -2800 to 3000 km, y runs 2500 to -2500 km:
      parameter (nxi=581, nyi=501, ddin=10.e3)
#endif

      dimension hbi(nxi,nyi), xi(nxi,nyi), yi(nxi,nyi),
     *          hbi2(nxi)
      character*1 carr(nxi,nyi)
      character cfilin*200, cdum*16

      cfilin = cpathhome(1:lenchr(cpathhome))
#if defined (SHALBEDMAP)
     *         // 'Data_shaldril/bednew4.xyz'
#elif defined (SHALGRL09)
     *         // 'Data_shaldril/EoceneDec08_10k.txt'
c    *         // 'Data_shaldril/EoceneDec08_20k.txt'
#elif defined (SHALP311MAX) 
     *         // 'Data_shaldril/ANTscape34MyrMaximumTopography.asc'
#elif defined (SHALP311MIN) || defined (SHALP311AVE)
     *         // 'Data_shaldril/ANTscape34MyrMinimumTopography.asc'
#elif defined (SHALMIOC)
     *         // 'Data_shaldril/18MaTopo-20120604.txt'
elif defined (SHALMIOC15)
     *         // 'Data_shaldril/miocene_min_10k'
c    *         // 'Data_shaldril/miocene_max_10k'
#elif defined (BEDMAP2RBD)
     *         // 'Data_bedmap2/BM2bedrbd10k.xyz'
#elif defined (SHALFILMIN)
     *         // 'Data_shaldril/bedEOminfil.xyz'
#endif

      write (6,'(2a)') 'Reading SHALDRIL elevation file ',
     *                 cfilin(1:lenchr(cfilin))

c        data's j runs "top to bottom", so convert to model's j
c        bottom to top

      open (iushal, file=cfilin,  form='formatted', status='old')
c-----------------------------------------------------------------------
#if defined (SHALBEDMAP) || defined (SHALGRL09) || defined (SHALFILMIN) || defined (SHALMIOC) || defined (SHALMIOC15) || defined (BEDMAP2RBD)
c-----------------------------------------------------------------------
      do j=1,nyi
        jj = nyi -j + 1
        do i=1,nxi
          read (iushal, *, end=200, err=200) 
     *      xi(i,jj), yi(i,jj), hbi(i,jj)                    ! 10 km 
c    *      xi(i,jj), yi(i,jj), hbi(i,jj), zlon, zlat        ! 20 km
        enddo
      enddo
c-----------------------------------------------------------------------
#elif defined (SHALP311MAX) || defined (SHALP311MIN) || defined (SHALP311AVE)
c-----------------------------------------------------------------------
      read (iushal,*) cdum, ncols
      read (iushal,*) cdum, nrows
      read (iushal,*) cdum, xllcorner  ! m
      read (iushal,*) cdum, yllcorner  ! m
      read (iushal,*) cdum, cellsize   ! m
      read (iushal,*) cdum, znodata    ! no nodata values in these files
      do j=1,nyi
        jj = nyi -j + 1
        read (iushal,*, end=200, err=200) (hbi(i,jj),i=1,nxi)
        do i=1,nxi
          xi(i,jj) = ( i-0.5 -0.5*nxi)*cellsize*.001  ! m to km
          yi(i,jj) = (jj-0.5 -0.5*nyi)*cellsize*.001  ! m to km
        enddo
      enddo
c-----
#endif
c-----
      close (iushal)

c        For SHALP311AVE, have read "min" (hbi) above. Now read "max" 
c        (hbi2), take average 0.5*(hbi+hbi2)

c~~~~~~~~~~~~~~~~~~~~~~~~
#if defined (SHALP311AVE) 
c~~~~~~~~~~~~~~~~~~~~~~~~
      open (iushal, file=cfilin,  form='formatted', status='old')
      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_shaldril/ANTscape34MyrMaximumTopography.asc'
      read (iushal,*) cdum, ncols
      read (iushal,*) cdum, nrows
      read (iushal,*) cdum, xllcorner  ! m
      read (iushal,*) cdum, yllcorner  ! m
      read (iushal,*) cdum, cellsize   ! m
      read (iushal,*) cdum, znodata    ! no nodata values in these files
      do j=1,nyi
        jj = nyi -j + 1
        read (iushal,*, end=200, err=200) (hbi2(i),i=1,nxi)
        do i=1,nxi
          hbi(i,jj) = 0.5*(hbi(i,jj) + hbi2(i))
        enddo
cdone   do i=1,nxi
cdone     xi(i,jj) = ( i-0.5 -0.5*nxi)*cellsize*.001  ! m to km
cdone     yi(i,jj) = (jj-0.5 -0.5*nyi)*cellsize*.001  ! m to km
cdone   enddo
      enddo
      close (iushal)
c~~~~~
#endif
c~~~~~

c        Simple "tile" average from data to model grid

c     write (6,'(a)') '  tiling data to model grid'

      do jo=1,ny
c       write (6,'(a,2i6)') '  jo,ny=',jo,ny
        do io=1,nx

c         find nearest data grid point:
c---
c Either:
c         zdmin = 1.e20
c         do j=1,nyi
c           do i=1,nxi
c             zd = (   (1000.*xi(i,j)-xh(io,jo)))**2
c    *               + (1000.*yi(i,j)-yh(io,jo)))**2 )
c             if (zd.lt.zdmin) then
c               ia = i
c               ja = j
c               zdmin = zd
c             endif
c           enddo
c         enddo
c---
c Or:
c         quicker, especially for very fine (~5 or 10 km) grids:
          ia = nint ( (xh(io,jo) - 1000.*xi(1,1))/ddin + 1. )
          ja = nint ( (yh(io,jo) - 1000.*yi(1,1))/ddin + 1. )
          ia = max (1, min (nxi, ia))
          ja = max (1, min (nyi, ja))
c---

c         average involving surrounding data grid points
c         (ignore variations in tile area overlaps)
          zav = 0.
          nav = 0.
          inc = nint(dd0/ddin)
          do j=ja-inc/2,ja+inc/2
            jj = max (1, min (nyi, j))
            do i=ia-inc/2,ia+inc/2
              ii = max (1, min (nxi, i))
              zav = zav + hbi(ii,jj)
              nav = nav + 1
            enddo
          enddo
          hb(io,jo) = zav/nav

        enddo
      enddo

c        Diagnostic ascii dump of model field

      do j=1,ny
        do i=1,nx
          icarr = nint(hb(i,j)/150.)
          if (icarr.eq.0 .and. hb(i,j).ge.0.) then
            carr(i,j) = '+'
          else if (icarr.eq.0 .and. hb(i,j).lt.0.) then
            carr(i,j) = '-'
          else if (icarr.gt.26) then
            carr(i,j) = '*'
          else if (icarr.lt.-26) then
            carr(i,j) = '#'
          else if (icarr.gt.0 .and. icarr.le.26) then
            carr(i,j) = char(ichar('A') + icarr - 1 )
          else if (icarr.lt.0 .and. icarr.ge.-26) then
            carr(i,j) = char(ichar('a') - icarr - 1)
          endif
        enddo
      enddo
      write (87,'(5x,500i5)') (i, i=5,nx,5)
      do j=ny,1,-1
        write (87,'(i3,2x,2000a1)') j, (carr(i,j),i=1,nx) 
      enddo
      write (87,*)
      do j=ny,1,-1
        write (87,'(500i3)') (nint(hb(i,j)/100.),i=5,nx,5) 
      enddo

      return

  200 write (6,'(2a,i4,a,i4)') 'Error reading SHALDRIL file',
     *                '  i=',i,'  j=',j

      stop
      end
c>>>>>
#endif
c>>>>>

c-----------------------------------------------------------------------

c>>>>>>>>>>>>>>>>>>>>>>>
#if defined (FORCE34MYR)
c>>>>>>>>>>>>>>>>>>>>>>>

      subroutine shaldril_mio (hb)

#include <comicegrid.h>
#include <comicepath.h>
      parameter (iushal = 60)

c     Reads Miocene (~15 Ma) topog/bathym, like shaldril_eo 
c     but just one choice (no -D's), identrical processing code. 
c     For time-varying topbedeq,equiload, called from subr bedrock.

      dimension hb(nx,ny)

c     Wilson Early Miocene.
c     10 km input x runs -3190 to 3200 km, y runs 2500 to -2600 km:
c     parameter (nxi=640, nyi=511, ddin=10.e3)                ! SHALMIOC
c     Gasson version of Wilson Early Miocene.
c     10 km input x runs -2800 to 3000 km, y runs 2500 to -2500 km:
      parameter (nxi=581, nyi=501, ddin=10.e3)              ! SHALMIOC15

      dimension hbi(nxi,nyi), xi(nxi,nyi), yi(nxi,nyi)
      character*1 carr(nxi,nyi)
      character cfilin*200, cdum*16

      cfilin = cpathhome(1:lenchr(cpathhome))
c    *         // 'Data_shaldril/18MaTopo-20120604.txt'       ! SHALMIOC
     *         // 'Data_shaldril/miocene_min_10k'           ! SHALMIOC15
c    *         // 'Data_shaldril/miocene_max_10k'           ! SHALMIOC15

      write (6,'(a/2a)') 'shaldril_mio:',
     *                   'Reading SHALDRIL elevation file ',
     *                   cfilin(1:lenchr(cfilin))

c        data's j runs "top to bottom", so convert to model's j
c        bottom to top

      open (iushal, file=cfilin,  form='formatted', status='old')

      do j=1,nyi
        jj = nyi -j + 1
        do i=1,nxi
          read (iushal, *, end=200, err=200) 
     *      xi(i,jj), yi(i,jj), hbi(i,jj)                    ! 10 km 
c    *      xi(i,jj), yi(i,jj), hbi(i,jj), zlon, zlat        ! 20 km
        enddo
      enddo

      close (iushal)

c        Simple "tile" average from data to model grid

c     write (6,'(a)') '  tiling data to model grid'

      do jo=1,ny
c       write (6,'(a,2i6)') '  jo,ny=',jo,ny
        do io=1,nx

c         find nearest data grid point:
c---
c Either:
c         zdmin = 1.e20
c         do j=1,nyi
c           do i=1,nxi
c             zd = (   (1000.*xi(i,j)-xh(io,jo)))**2
c    *               + (1000.*yi(i,j)-yh(io,jo)))**2 )
c             if (zd.lt.zdmin) then
c               ia = i
c               ja = j
c               zdmin = zd
c             endif
c           enddo
c         enddo
c---
c Or:
c         quicker, especially for very fine (~5 or 10 km) grids:
          ia = nint ( (xh(io,jo) - 1000.*xi(1,1))/ddin + 1. )
          ja = nint ( (yh(io,jo) - 1000.*yi(1,1))/ddin + 1. )
          ia = max (1, min (nxi, ia))
          ja = max (1, min (nyi, ja))
c---

c         average involving surrounding data grid points
c         (ignore variations in tile area overlaps)
          zav = 0.
          nav = 0.
          inc = nint(dd0/ddin)
          do j=ja-inc/2,ja+inc/2
            jj = max (1, min (nyi, j))
            do i=ia-inc/2,ia+inc/2
              ii = max (1, min (nxi, i))
              zav = zav + hbi(ii,jj)
              nav = nav + 1
            enddo
          enddo
          hb(io,jo) = zav/nav

        enddo
      enddo

c        Diagnostic ascii dump of model field

      do j=1,ny
        do i=1,nx
          icarr = nint(hb(i,j)/150.)
          if (icarr.eq.0 .and. hb(i,j).ge.0.) then
            carr(i,j) = '+'
          else if (icarr.eq.0 .and. hb(i,j).lt.0.) then
            carr(i,j) = '-'
          else if (icarr.gt.26) then
            carr(i,j) = '*'
          else if (icarr.lt.-26) then
            carr(i,j) = '#'
          else if (icarr.gt.0 .and. icarr.le.26) then
            carr(i,j) = char(ichar('A') + icarr - 1 )
          else if (icarr.lt.0 .and. icarr.ge.-26) then
            carr(i,j) = char(ichar('a') - icarr - 1)
          endif
        enddo
      enddo
      write (87,'(5x,500i5)') (i, i=5,nx,5)
      do j=ny,1,-1
        write (87,'(i3,2x,2000a1)') j, (carr(i,j),i=1,nx) 
      enddo
      write (87,*)
      do j=ny,1,-1
        write (87,'(500i3)') (nint(hb(i,j)/100.),i=5,nx,5) 
      enddo

      return

  200 write (6,'(2a,i4,a,i4)') 
     *  'shaldril_mio: Error reading SHALDRIL file',
     *  '  i=',i,'  j=',j

      stop
      end
c>>>>>
#endif
c>>>>>

c-----------------------------------------------------------------------

c^^^^^^^^^^^^^^^^^^^^^^^
#if defined (DYNAMTOPOG)
c^^^^^^^^^^^^^^^^^^^^^^^

      subroutine dyntop_read (delhb)

c     Reads dynamic topography file (Jackie Austermann/Jerry Mitrovica) 
c     for difference in bed topog due to dynamic topography
c     (modern minus 3 Ma). Interp to ice grid, returned in delhb. 

#include <comicegrid.h>
#include <comicepath.h>

      parameter (iudyn=60, iudiag=157)

      dimension delhb(nx,ny)

      character cfilin*200, cname*200

c     parameter (nxdyn=360*4, nydyn=180*4)                ! older
      parameter (nxdyn=(360*2 + 1), nydyn=(180*2+1))      ! newer
      dimension alondynd(nxdyn), alatdynd(nydyn), 
     *          alondyn(nxdyn),  alatdyn(nydyn),  
     *          deldyn(nxdyn,nydyn)

      dimension
     *  inddyna(nx,ny), weidyna(nx,ny), inddynb(nx,ny), weidynb(nx,ny)

      parameter (vnull=999999.)

      write (ioterm,'(a)')
     *  'reading existing dynamic topography data file'

c        Read existing dynamic topog data file

#if defined (DYNAMFILE)
      cname = DYNAMFILE
#else
      cname = 'Ant_delT'
#endif
      if (cname.eq.'nodyntop') then
        call zero (delhb, nx*ny)
        return
      endif
  
      cfilin = cpathhome(1:lenchr(cpathhome)) // 'Data_dyntop/' // 
     *         cname(1:lenchr(cname)) // '_nulls.txt'
      open (iudyn, file=cfilin, form='formatted', status='old')
c     do j=1,nydyn      ! older
c       do i=1,nxdyn    ! older
      do i=1,nxdyn      ! newer
        do j=1,nydyn    ! newer
          read (iudyn,*) alondynd(i), alatdynd(j), deldyn(i,j)
          alondyn(i) = alondynd(i)*pi/180. 
          alatdyn(j) = alatdynd(j)*pi/180. 
        enddo
      enddo
      close (iudyn)

c       Diagnostic output

      write (iudiag,'(/a)') 'latitudes:'
      do j=1,nydyn
        write (iudiag,'(i6,f10.2)') j, alatdynd(j)
      enddo
      write (iudiag,'(/a)') 'longitudes:'
      do i=1,nxdyn
        write (iudiag,'(i6,f10.2)') i, alondynd(i)
      enddo
      write (iudiag,'(/a)') '|deldyn|/50:'
      do j=nydyn,1,-4
        write (iudiag,'(f8.3,2x,360i1)') 
     *    alatdynd(j), (nint(abs(deldyn(i,j))/50.), i=1,nxdyn,4)
      enddo

c        Interpolate to ice grid (like ocn_int in iceocean).

c        First set indices for bilinear interp (inddyn*,weidyn*)
c        from ice grid into ogcm grid

      call setinterp (alondyn, alatdyn, nxdyn, nydyn,
     *                alond, alatd, nx,    ny,
     *                inddyna, weidyna, inddynb, weidynb, 0.)

c        Interpolate non-null (deldyn) to ice model grid (delhb)
    
      do j=1,ny
        do i=1,nx
          il = inddyna(i,j)
          ir = mod (il,nxdyn) + 1
          jb = inddynb(i,j)
          jt = min (jb+1, nydyn)
          wlb =     weidyna(i,j) *    weidynb(i,j)
          wrb = (1.-weidyna(i,j))*    weidynb(i,j)
          wlt =     weidyna(i,j) *(1.-weidynb(i,j))
          wrt = (1.-weidyna(i,j))*(1.-weidynb(i,j))
          zarr = 0.
          zwei = 0.
          if (deldyn(il,jb).ne.vnull) then
            zarr = zarr + deldyn(il,jb)*wlb
            zwei = zwei +               wlb
          endif
          if (deldyn(ir,jb).ne.vnull) then
            zarr = zarr + deldyn(ir,jb)*wrb
            zwei = zwei +               wrb
          endif
          if (deldyn(il,jt).ne.vnull) then
            zarr = zarr + deldyn(il,jt)*wlt
            zwei = zwei +               wlt
          endif
          if (deldyn(ir,jt).ne.vnull) then
            zarr = zarr + deldyn(ir,jt)*wrt
            zwei = zwei +               wrt
          endif

          if (zwei.gt.0.) then
            delhb(i,j) = zarr/zwei
          else
            delhb(i,j) = 0.        ! if surrounded by nulls, set to zero
          endif

        enddo   
      enddo   

      call printmap (0., delhb, 'delhb (dyntop)', 20., iudiag, 0)
      call flush (iudiag)

      return
      end

c^^^^^
#endif
c^^^^^

c-----------------------------------------------------------------------

#if defined (EISMINT2)

      subroutine readeisant (code, arr)

c     Reads a 2-D field (arr) from an EISMINT II Antarctic data file. 
c     code = BED, SUR, THI or MAS.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

      character*(*) code
      dimension arr(nx,ny)

      character*240 cpath, cfile, fmt

      cpath = cpathhome(1:lenchr(cpathhome)) 
     *        // 'Eismintantarc/' /

      if (nx.ne.141 .or. ny.ne.141) then
        write (ioterm,*)' Error (readeisant): nx and ny must equal 141'
        stop
      endif

c        Set full name

      if (code.eq.'BED') then
        cfile = 'bed40eis.dat'
      else if (code.eq.'SUR') then
        cfile = 'sur40eis.dat'
      else if (code.eq.'THI') then
        cfile = 'thi40eis.dat'
      else if (code.eq.'MAS') then
        cfile = 'mask40eis.dat'
      else 
        write (ioterm,*) '*** Error: Unknown code in readeisant: ',code
        stop
      endif

      write (ioterm, '(/2a)')
     *   'Reading EISMINT data file ', cfile(1:lenchr(cfile))

      cfile = cpath(1:lenchr(cpath)) // cfile

c        Read file

      open (iueis, file=cfile, status='old', err=1000, iostat=ios)  
      if (ios.ne.0) goto 1000

      read (iueis,'(/a)') fmt
      do j=1,ny
        read (iueis,fmt,err=2000,iostat=ios) kdum, (arr(i,j),i=1,nx)
        if (ios.ne.0) goto 2000
      enddo

      close(iueis)

c        Adjust values 

      if (code.eq.'SUR') then
        do j=1,ny
          do i=1,nx
            if (nint(arr(i,j)).eq.25) arr(i,j) = 0.
          enddo
        enddo
      else if (code.eq.'THI') then
        do j=1,ny
          do i=1,nx
            if (nint(arr(i,j)).eq.218) arr(i,j) = 0.
          enddo
        enddo
      endif

      return

 1000 write (ioterm, '(/2a/a,i6)')
     *   '*** Error opening EISMINT data file ', cfile(1:lenchr(cfile)),
     *   '    ios = ',ios
      stop

 2000 write (ioterm, '(/2a/a,i6/a,i6)')
     *   '*** Error reading EISMINT data file ', cfile(1:lenchr(cfile)),
     *   '    ios = ',ios,
     *   '    j   = ',j

      stop
      end

#endif

c-----------------------------------------------------------------------

#if defined (NHA)

      subroutine readnha (arr, codein, cperin, nlonin, nlatin)

c     Reads a field (arrin) from a global Genesis data input file,
c     transfers it to arr on ice grid for "NHA" northern hemis runs.
c     code = 'SURT', 'TOPO', or 'TOPOBED', 'SURTSED' for ICE5G.
c     cperin  = '21k', '15k', '0k' for ICE5G.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

      dimension arr(nx,ny)
      character*(*) codein, cperin

      dimension arrin (nlonin,nlatin)
      character cpath*240, cfile*240, cres*4, cper*8
      character*8 code

      dimension alatt106(160)
      save alatt106
      data alatt106 /
     *  -1.5558130142, -1.5364033724, -1.5168791062, -1.4973289414,
     *  -1.4777689527, -1.4582042102, -1.4386368110, -1.4190677780,
     *  -1.3994976690, -1.3799268140, -1.3603554207, -1.3407836263,
     *  -1.3212115249, -1.3016391835, -1.2820666508, -1.2624939632,
     *  -1.2429211486, -1.2233482282, -1.2037752191, -1.1842021348,
     *  -1.1646289861, -1.1450557820, -1.1254825297, -1.1059092353,
     *  -1.0863359041, -1.0667625402, -1.0471891474, -1.0276157288,
     *  -1.0080422871, -0.9884688248, -0.9688953437, -0.9493218457,
     *  -0.9297483325, -0.9101748052, -0.8906012653, -0.8710277137,
     *  -0.8514541515, -0.8318805795, -0.8123069985, -0.7927334092,
     *  -0.7731598122, -0.7535862082, -0.7340125976, -0.7144389809,
     *  -0.6948653585, -0.6752917309, -0.6557180984, -0.6361444614,
     *  -0.6165708201, -0.5969971748, -0.5774235259, -0.5578498735,
     *  -0.5382762178, -0.5187025592, -0.4991288978, -0.4795552337,
     *  -0.4599815671, -0.4404078983, -0.4208342273, -0.4012605543,
     *  -0.3816868795, -0.3621132029, -0.3425395247, -0.3229658451,
     *  -0.3033921640, -0.2838184817, -0.2642447981, -0.2446711136,
     *  -0.2250974280, -0.2055237415, -0.1859500542, -0.1663763662,
     *  -0.1468026775, -0.1272289883, -0.1076552986, -0.0880816085,
     *  -0.0685079181, -0.0489342274, -0.0293605365, -0.0097868455,
     *   0.0097868455,  0.0293605365,  0.0489342274,  0.0685079181,
     *   0.0880816085,  0.1076552986,  0.1272289883,  0.1468026775,
     *   0.1663763662,  0.1859500542,  0.2055237415,  0.2250974280,
     *   0.2446711136,  0.2642447981,  0.2838184817,  0.3033921640,
     *   0.3229658451,  0.3425395247,  0.3621132029,  0.3816868795,
     *   0.4012605543,  0.4208342273,  0.4404078983,  0.4599815671,
     *   0.4795552337,  0.4991288978,  0.5187025592,  0.5382762178,
     *   0.5578498735,  0.5774235259,  0.5969971748,  0.6165708201,
     *   0.6361444614,  0.6557180984,  0.6752917309,  0.6948653585,
     *   0.7144389809,  0.7340125976,  0.7535862082,  0.7731598122,
     *   0.7927334092,  0.8123069985,  0.8318805795,  0.8514541515,
     *   0.8710277137,  0.8906012653,  0.9101748052,  0.9297483325,
     *   0.9493218457,  0.9688953437,  0.9884688248,  1.0080422871,
     *   1.0276157288,  1.0471891474,  1.0667625402,  1.0863359041,
     *   1.1059092353,  1.1254825297,  1.1450557820,  1.1646289861,
     *   1.1842021348,  1.2037752191,  1.2233482282,  1.2429211486,
     *   1.2624939632,  1.2820666508,  1.3016391835,  1.3212115249,
     *   1.3407836263,  1.3603554207,  1.3799268140,  1.3994976690,
     *   1.4190677780,  1.4386368110,  1.4582042102,  1.4777689527,
     *   1.4973289414,  1.5168791062,  1.5364033724,  1.5558130142/

      cper = cperin

      if (nlonin.eq.360 .and. nlatin.eq.360) then
        cres = '1XH'
        dlonin = 1.
        dlatin = 0.5
      else if (nlonin.eq.360 .and. nlatin.eq.180) then
        cres = '1X1'
        dlonin = 1.
        dlatin = 1.
      else if (nlonin.eq.180 .and. nlatin.eq.90) then
        cres = '2X2'
        dlonin = 2.
        dlatin = 2.
      else if (nlonin.eq.320 .and. nlatin.eq.160) then
        cres = 'T106'
        dlonin = 320./nlonin
        dlatin = 160./nlatin
c_sk  else if (nlonin.eq.2161 .and. nlatin.eq.1081) then
c_sk    S. Koenig 10 minute:
c_sk    cres = '10m'
c_sk    dlonin = 360./nlonin       ! dlonin =  0.166512488           
c_sk    dlatin = 180./nlatin       ! dlatin =  0.166512488
      else
        write (ioterm,'(2a,2i6)') 
     *    'Error (readnha): unknown input file resolution.',
     *    ' nlonin, nlatin = ', nlonin, nlatin
        stop
      endif

      code = codein
      cpath = cpathmap(1:lenchr(cpathmap))

      if (codein.eq.'SURT') then
        cfile = cpath(1:lenchr(cpath)) 
c    *        // '34ma/data_34ma_sur_' // cres
c_sk *        // 'V2/data_surnoice_' // cres
c_sk *        // 'V2/data_surnoice_openhudson_' // cres
     *        // 'Ice5g/data_ice5g_' // cper(1:lenchr(cper))
     *        // '_sur_' // cres
        icode = 1

      else if (codein.eq.'TOPO') then
        cfile = cpath(1:lenchr(cpath)) 
c    *          // '34ma/data_34ma_top_' // cres
c_sk *          // 'V2/data_topnoice_'   // cres
     *          // 'Ice5g/data_ice5g_' // cper(1:lenchr(cper))
     *          // '_top_' // cres
        icode = 2

      else if (codein.eq.'TOPOBED') then
        cfile = cpath(1:lenchr(cpath)) 
     *          // 'Ice5g/data_ice5g_' // cper(1:lenchr(cper))
     *          // '_topbed_' // cres
        icode = 2

      else if (codein.eq.'SURTSED') then
        cfile = cpath(1:lenchr(cpath)) 
     *          // 'Ice5g/data_ice5g_' // cper(1:lenchr(cper))
     *          // '_sur_' // cres
        code = 'SURT'
        icode = 1

      else 
        write (ioterm,*) 'Error (readnha): unknown code ',code
        stop
      endif

      write (6,'(4a)') 
     *  ' readnha: reading ', cfile(1:lenchr(cfile)),'  code=',code
      open (iunh, file = cfile, form='formatted', status='old')
c     no long shift, keep -180 to 180 as in data file:
      call readdata_ice (iunh, arrin, nlonin, nlatin, nlonin, nlatin,
     *                   code, icode, 1)
      close(iunh)

c        Tile (LONLAT) or nearest neighbor (STEREO) to ice grid

      
      iwid = 0.5*(dlond/dlonin) + .001
      jwid = 0.5*(dlatd/dlatin) + .001

      do j=1,ny 
        do i=1,nx 

c         set nearest neighbor ig,jg:
          if (cres.eq.'T106') then
c           T106 grid is gaussian:
            zlat = alatd(i,j)*pi/180.
            do jj=1,nlatin-1
              if (zlat .le. 0.5*(alatt106(jj) + alatt106(jj+1))) then
                jg = jj
                goto 10
              endif
           enddo
           jg = nlatin
   10      continue
          else
            jg = nint ((alatd(i,j) +  90.)/dlatin + 0.5)
            jg = max (1, min (nlatin, jg))
          endif

          zlon = mod (alond(i,j) + 360.,360.) 
          if (zlon.gt.180.) zlon = zlon - 360.
          ig = nint ((zlon + 180.)/dlonin + 0.5) 
          ig = max (1, min (nlonin, ig))

#if defined (LONLAT) 
c         tile:      
          zval = 0.
          nval = 0
          do jg2 = jg-jwid,jg+jwid
            jg3 = max (1, min (nlatin,jg2))
            do ig2 = ig-iwid,ig+iwid
              ig3 = mod (ig2 - 1 + nlonin, nlonin) + 1
              zval = zval + arrin(ig3,jg3)
              nval = nval + 1
            enddo
          enddo
          arr(i,j) = zval/nval
          if (code.eq.'SURT') then
            arr(i,j) = nint(arr(i,j))
            arr(i,j) = max (0., min (2., arr(i,j)))
          endif
#elif defined (STEREO)
c         nearest neighbor:
          arr(i,j) = arrin(ig,jg)
#endif
        enddo
      enddo

c     diagnostic map dump:
c     if (.true.) then
      if (.false.) then
        if (code.eq. 'SURT') then
          zsc = 1.
        else 
          zsc = 100.
        endif
        write (177,*) codein, '  input grid', nlonin, nlatin
        do j=nlatin,1,-1
         write(177,'(i5,3x,5000i1)') j,(nint(arrin(i,j)/zsc),i=1,nlonin)
        enddo
        write (177,*) codein, '  ice grid', nx, ny
        do j=ny,1,-1
          write (177,'(i5,3x,1000i1)') j, (nint(arr(i,j)/zsc),i=1,nx)
        enddo
        write (177,*)
        call flush(177)
      endif

      return
      end

#endif

c-----------------------------------------------------------------------

#if defined (CARB)

      subroutine readcarb (arr, code)

c     Reads a field (arrin) from a global Genesis data input file,
c     transfers it to arr on ice grid for "CARB" global hemis runs.
c     code = 'SURT' or 'TOPO' (key in Genesis data input file header).

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

      character*(*) code
      dimension arr(nx,ny)

c     ddin is the grid size (lat and lon, degrees) of the data input
c     files. Can be same or different from the ice model's.
c     parameter (ddin = dlatd)
      parameter (ddin = 1.)

      parameter (nlonin = 360/ddin + .001,
     *           nlatin = 180/ddin + .001)
      dimension arrin (nlonin,nlatin)
      character cpath*240, cfile*240, cres*3, cform*80

      if (nint(ddin).eq.1) then
        cres = '1X1'
      else if (nint(ddin).eq.2) then
        cres = '2X2'
      else
        write (ioterm,*) 'Error (readnha): resol must be 1 or 2 deg'
        stop
      endif

      cpath = cpathmap(1:lenchr(cpathmap))

      if (code.eq.'SURT') then
        cfile = cpath(1:lenchr(cpath))// '280ma/data_sur_280ma_' // cres
        icode = 1
      else if (code.eq.'TOPO') then
        cfile = cpath(1:lenchr(cpath))// '280ma/data_top_280ma_' // cres
        icode = 2
      else 
        write (ioterm,*) 'Error (readnha): unknown code ',code
        stop
      endif

      write (6,'(4a)') 
     *  ' readnha: reading ', cfile(1:lenchr(cfile)),'  code=',code
      open (iunh, file = cfile, form='formatted', status='old')
c     no long shift, keep -180 to 180 as in data file:
      call readdata_ice (iunh, arrin, nlonin, nlatin, nlonin, nlatin,
     *                   code, icode, 1)
      close(iunh)

c        Use nearest neighbor for now

      do j=1,ny 
        do i=1,nx 
          jg = nint ((alatd(i,j) +  90.)/ddin + 0.5)
          jg = max (1, min (nlatin, jg))

          zlon = mod (alond(i,j) + 360.,360.) 
          if (zlon.gt.180.) zlon = zlon - 360.
c         if data input file long resol = 2* model's, tile:
          if (nlonin.eq.2*nint(360./dlond)) then 
            ig  = nint ((zlon + 180.)/ddin)
            ig  = max (1, min (nlonin, ig))
            igp = mod (ig,nlonin) + 1
            arr(i,j) = 0.5*(arrin(ig,jg) + arrin(igp,jg))
          else
            ig = nint ((zlon + 180.)/ddin + 0.5 ) 
            ig = max (1, min (nlonin, ig))
            arr(i,j) = arrin(ig,jg)
          endif
        
        enddo
      enddo

c     diagnostic map dump:
c     if (code.eq. 'SURT') then
c       cform = '(i5,3x,360i1)'
c     else
c       cform = '(i5,3x,360i5)'
c     endif
c     write (177,*) code, '  global', nlonin, nlatin
c     do j=nlatin,1,-1
c       write (177,cform) j, (nint(arrin(i,j)),i=1,nlonin)
c     enddo
c     write (177,*) code, '  ice grid', nx, ny
c     do j=ny,1,-1
c       write (177,cform) j, (nint(arr(i,j)),i=1,nx)
c     enddo
c     call flush(177)

      return
      end

#endif

c-----------------------------------------------------------------------

#if defined (CLEMCLIM)

      subroutine clem_masscon (iu, cfilin, cvarin, field)

c     Reads modern mass-conserving data file for Antarctic Peninsula
c     (Clemens Schannwell), overrides bedrock elevations and
c     ice thickness of Bedmap2 within this subdomain (with CLEMCLIM) 

#include <comicephys.h>
#include <comicegrid.h>
#include <netcdf.inc>

      character*(*) cfilin, cvarin
      dimension field (nx,ny)

      integer*4 ncid, rcode, varid, start(4), count(4)
      character cfile*240, cvar*16

      parameter (nxcor=1024, nycor=1024) 
      parameter (dxcor=1.e3, dycor=1.e3)
      dimension arrin(nxcor,nycor), brrin(nycor,nxcor),
     *          xcor(nxcor), ycor(nycor)
      save xcor, ycor

      dimension kount(nx,ny), arrav(nx,ny), fieldold(nx,ny)

      logical firstcor
      data firstcor /.true./
      save firstcor

      cfile = cfilin
      cvar = cvarin

      write (6,'(2a)')
     *  'reading existing Netcdf file ', cfile(1:lenchr(cfile))

c        Open existing netcdf file and read requested field

      rcode = nf_open(cfile, nf_nowrite, ncid)
      call checkrcode (rcode, cfile, 'open', 6000)

c        If first call, read x and y axes (polar stereo grid,
c        despitye being called lat and lon on the file)

c-----------------------
      if (firstcor) then
c-----------------------

c          x axis:

        start(1) = 1
        count(1) = nxcor
        rcode = nf_inq_varid (ncid, 'lon', varid)
        call checkrcode (rcode, cfile, cvar, 7720)
        rcode = nf_get_vara_double (ncid, varid, start, count, xcor)
        call checkrcode (rcode, cfile, cvar, 7721)

c          y axis:

        start(1) = 1
        count(1) = nycor
        rcode = nf_inq_varid (ncid, 'lat', varid)
        call checkrcode (rcode, cfile, cvar, 7730)
        rcode = nf_get_vara_double (ncid, varid, start, count, ycor)
        call checkrcode (rcode, cfile, cvar, 7731)

        firstcor = .false.

c----------
      endif
c----------

c        bed elevations or ice thicknesses:

      start(1) = 1
      count(1) = nycor
      start(2) = 1
      count(2) = nxcor

      if (cvar.eq.'topob') then
        rcode = nf_inq_varid (ncid, 'topob', varid)
        call checkrcode (rcode, cfile, 'topob', 7735)
        rcode = nf_get_vara_double (ncid, varid, start, count, brrin)
        call checkrcode (rcode, cfile, 'topob', 7736)
      else if (cvar.eq.'thk') then
        rcode = nf_inq_varid (ncid, 'thk', varid)
        call checkrcode (rcode, cfile, 'thk', 7735)
        rcode = nf_get_vara_double (ncid, varid, start, count, brrin)
        call checkrcode (rcode, cfile, 'thk', 7736)
      else
        write (6,'(/2a)') '*** Error: unknown field ',cvar
        stop
      endif

c        Switch x and y axes (lon and lat indices are switched from 
c        normal order on the file)

      do jj=1,nycor 
        do ii=1,nxcor 
c         arrin(ii,jj) = brrin(jj,ii)
          arrin(ii,jj) = brrin(ii,jj)  ! don't switch - ok on file(?)
        enddo
      enddo

c       "Tile" onto ice model grid, (assuming each 1x1 km mass-con
c       cell fits exactly into an ice model cell, no overlaps).
c       Use a temporary ice model array (arrav) to accumulate tiling. 

      call izero (kount, nx*ny)
      call zero (arrav, nx*ny)

      xorig = -0.5*nx*dx0 + xoffa
      yorig = -0.5*ny*dy0 + yoffa

      do jj=1,nycor
        j = int((ycor(jj)-yorig)/dy0) + 1
        if (j.ge.1 .and. j.le.ny) then
          do ii=1,nxcor
            i = int((xcor(ii)-xorig)/dx0) + 1
            if (i.ge.1 .and. i.le.nx) then
              kount(i,j) = kount(i,j) + 1
              arrav(i,j) = arrav(i,j) + arrin(ii,jj)
            endif
          enddo
        endif
      enddo

      do j=1,ny
        do i=1,nx
          arrav(i,j) = arrav(i,j) / max (kount(i,j), 1)
        enddo
      enddo

c       Transfer to ice model array, only for ice model cells
c       completely within the mass-con domain (i..e, fully tiled) 

      call scopy (nx*ny, field, 1, fieldold, 1)

      nk = max (1, nint( (dx0*dy0)/(dxcor*dycor) ))
      do j=1,ny
        do i=1,nx
          if (kount(i,j).eq.nk) field(i,j) = arrav(i,j)
        enddo
      enddo

      rcode = nf_close (ncid)
        call checkrcode (rcode, cfile, 'close', 7770)
  
c     diagnostic output:
c     if (.true.) then
      if (.false.) then
        write (333,'(/a, 2f14.1)') 'xoffa, xorig=', xoffa, xorig
        write (333,'(/a, 2f14.1)') 'yoffa, yorig=', yoffa, yorig
        write (333,'(/a)') 'xcor'
        write (333,'(i4,f14.1)') (i, xcor(i),i=1,nxcor,10)
        write (333,'(/a)') 'ycor'
        write (333,'(i4,f14.1)') (j, ycor(j),j=1,nycor,10)
  
        write (333,'(/2a)') 'arrin ', cvar
        do jj=nycor,1,-10
          write (333,'(500i4)') (nint(.01*arrin(ii,jj)),ii=1,nxcor,10)
        enddo
  
        write (333,'(/2a)') 'arrav ', cvar
        do j=ny,1,-2
          write (333,'(500i4)') (nint(.01*arrav(i,j)),i=1,nx,2)
        enddo
  
        write (333,'(/2a)') 'field ', cvar
        do j=ny,1,-2
          write (333,'(500i4)') (nint(.01*field(i,j)),i=1,nx,2)
        enddo
  
        write (333,'(/2a)') 'fieldold ', cvar
        do j=ny,1,-2
          write (333,'(500i4)') (nint(.01*fieldold(i,j)),i=1,nx,2)
        enddo
  
        write (333,'(/a)') 'kount'
        do j=ny,1,-2
          write (333,'(500i2)') (kount(i,j),i=1,nx,2)
        enddo
c       ier = flush(333)
      endif

      return
      end

#endif

c-----------------------------------------------------------------------

#if defined (ORDO)

      subroutine readordo (iu, field, cvarin, timeice)

c     Read variable cvarin (LANDFRAC or PHIS) from CSM Ordovician "grid"
c     file, interpolate to ice grid, return in field 

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

      character*(*) cvarin
      dimension field (nx,ny)

      character cfile*240, cvar*16, cpath*240

      dimension arrg(nlong, nlatg)
      dimension
     *  indl(nx,ny),     weil(nx,ny),
     *  indb(nx,ny),     weib(nx,ny)
      save indl, weil, indb, weib

      logical firstordo
      data firstordo /.true./
      save firstordo

      cvar = cvarin

      cpath = cpathhome(1:lenchr(cpathhome)) // 'Data_csmordo/'

      call readdrive_csmordo (iu, cpath, arrg, cvar, 1, timeice)

c------------------------
      if (firstordo) then
c------------------------
        call setinterp (along, alatg, nlong, nlatg,
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, 0.)
        firstordo = .false.
c----------
      endif
c----------

      call intg2i (arrg, nlong, nlatg, field,
     *             indl, weil, indb, weib, nx, ny, 1)


c     diagnostic output:
c     if (.true.) then
      if (.false.) then
        if (cvar.eq.'LANDFRAC') zscale = 10.
        if (cvar.eq.'PHIS') zscale = .01
        write (333,'(/2a)') 'arrg ', cvar
        do jj=nlatg,1,-1
          write (333,'(500i3)') (nint(zscale*arrg(ii,jj)),ii=1,nlong)
        enddo
  
        write (333,'(/2a)') 'field ', cvar
        do j=ny,1,-2
          write (333,'(500i3)') (nint(zscale*field(i,j)),i=1,nx,2)
        enddo
        ier = flush(333)
      endif

      return
      end

#endif

c-----------------------------------------------------------------------

      subroutine readres (h, hs, hb, t, tracer,
     *                    baseperc, basefrml, heatb, 
     *                    sedim, tsed, wsed, heats, quarryacc, hbinit,
     *                    topbed, topbedeq, equiload, tbed, hw, tw, 
     *                    sedpres, crhmel, enhanceinv,
     *                    ua, va, ub, vb, dfu, dfv, dbu, dbv,
#if defined (MELANGE)
     *                    hmel,
#endif
     *                    sealev, rco2inter, timeicein, iczoom,
     *                    crestin)

c     Read restart file named crestin

#include <comicephys.h>
#include <comicegrid.h>

c     h grid:
      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny),
     *  t(nx,ny,0:nlevp),    tracer(nx,ny,0:nlevp,ntrace),
     *  baseperc(nx,ny),     basefrml(nx,ny),    heatb(nx,ny), 
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),   
     *  heats(nx,ny,nsed),   quarryacc(nx,ny),   hbinit(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny),    equiload(nx,ny),
     *  tbed(nx,ny,nbed),    hw(nx,ny),          tw(nx,ny),
     *  sedpres(nx,ny),      crhmel(nx,ny),      enhanceinv(nx,ny)

c     u grid:
      dimension
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp),
     *  dbu(0:nxp,0:nyp),       dbv(0:nxp,0:nyp)

#if defined (MELANGE)
     * ,hmel(nx,ny)
#endif

      character*(*) crestin

      dimension work(nx,ny), iwork(nx,ny), worku(0:nxp,0:nyp)

      character crest*120, cvar*10

      crest = crestin

      write (ioterm,'(/a,a)') 
     *  'reading restart file ', crest(1:lenchr(crest))

c        Open restart file, check version number 

      open (iures, file=crest, status='old', form='formatted')
      read (iures,'(e25.15)') versresin
      write (ioterm,'(a,f8.2)') '  versresin =', versresin
      if (       versresin.ne.1.
     *     .and. versresin.ne.2.
     *     .and. versresin.ne.3.
     *     .and. versresin.ne.4.
     *     .and. versresin.ne.5.
     *     .and. versresin.ne.6.
     *     .and. versresin.ne.7.
     *     .and. versresin.ne.8.
     *     .and. versresin.ne.9.
     *     .and. versresin.ne.10.
     *     .and. versresin.ne.11.
     *     .and. versresin.ne.12.
     *     .and. versresin.ne.13.
     *     .and. versresin.ne.14.
     *     .and. versresin.ne.15.
     *     .and. versresin.ne.16.
     *     .and. versresin.ne.17.
     *     .and. versresin.ne.18.
     *     .and. versresin.ne.19.
     *     .and. versresin.ne.20.
     *     .and. versresin.ne.21.
     *     .and. versresin.ne.22.
     *     .and. versresin.ne.23.
     *     .and. versresin.ne.24.
     *   ) then
        write(ioterm,*)
     *    '*** Error: restart file version # must be 1 to 24'
        stop
      endif

c        Read dates, resolutions

      if (versresin.le.9) then
        read (iures,'(e25.15,6i10)') 
     *    timeicein, nxin, nyin, nlevin, nsedin, nbedin, nwrit_dum
        ntracein = 1
      else if (versresin.le.16) then
        read (iures,'(e25.15,6i10)') 
     *    timeicein, nxin, nyin, nlevin, nsedin, nbedin
        ntracein = 1
      else
        read (iures,'(e25.15,6i10)') 
     *    timeicein, nxin, nyin, nlevin, nsedin, nbedin, ntracein
      endif
      write (ioterm,'(a, i10)') '  timeicein  =', nint(timeicein)

      if (nint(versresin).ge.18) then
        read (iures,'(4e25.15)') dxin, dyin, xoffin, yoffin
      else if (nint(versresin).ge.8) then
        read (iures,'(2e25.15)') dxin, dyin
        xoffin = 0.
        yoffin = 0.
c       just a guess:
        if (nxin.eq.nx  .and. nyin.eq.ny .and. 
     *      dxin.eq.dx0 .and. dyin.eq.dy0) then
          xoffin = xoffa
          yoffin = yoffa
        endif
      else
        if (nxin.eq.280 .or. nxin.eq.281) then
          dxin = 20.e3
          dyin = 20.e3
        else 
          dxin = 40.e3  ! edit if necessary (for old restart files)
          dyin = 40.e3  ! edit if necessary (for old restart files)
        endif
        xoffin = 0.
        yoffin = 0.
      endif

      if (nint(versresin).ge.11) then
        read (iures,'(e25.15)') (xhin(m),m=1,nxin*nyin)
        read (iures,'(e25.15)') (yhin(m),m=1,nxin*nyin)
      else
        do j=1,nyin
          do i=1,nxin
            m = (j-1)*nxin + i
#if defined (EISLINE) || defined (TEST2D) || defined (TESTMEL2D)
            xhin(m) = dxin*(i-0.5)
            yhin(m) = dyin*(j-0.5)
#else
            xhin(m) = dxin*(i-0.5) - 0.5*dxin*nxin
            yhin(m) = dyin*(j-0.5) - 0.5*dyin*nxin
#endif
          enddo
        enddo
      endif

      write (ioterm,'(a,2i8,2f12.3)') 
     *  '  nxin, nyin, dxin, dyin =', nxin, nyin, dxin, dyin
      write (ioterm,'(a,2i8,2f12.3)') 
     *  '  nx,   ny,   dx,   dy   =', nx,   ny,   dx0,  dy0

      if ((nxin+2)*(nyin+2).gt.ninmx) then
         write (ioterm,'(/a/a,2i8/a,i8/a,i8)') 
     *     ' *** Error (readres): input grid too large:',
     *     '      nxin, nyin        =', nxin, nyin, 
     *     '      (nxin+2)*(nyin+2) =', (nxin+2)*(nyin+2),
     *     '      ninmx             =', ninmx
        stop
      endif

      if ( nlev.ne.nlevin .or. nsed.ne.nsedin .or. nbed.ne.nbedin ) then
        write (ioterm,'(/a,4(/a,i4,a,i4))')
     *     '*** Error (readres): dimension mismatch:',
     *     '    nlev=',nlev,'  nlevin=',nlevin,
     *     '    nsed=',nsed,'  nsedin=',nsedin,
     *     '    nbed=',nbed,'  nbedin=',nbedin
        stop
      endif

#if defined (EISANTA) && defined (TRANSECTA) && defined (ZOOMGL)
      if (nx.ne.nxin .or. ny.ne.nyin) then
        write(ioterm,'(/a)')'*** Error: cannot change nx,ny with ZOOMGL'
        stop
      else
        read (iures) iczoom
        call initgrid (iczoom) ! reset to current zoom grid
      endif
#else
      iczoom = 0
#endif

      if (nint(versresin).ge.16) then
        read (iures,'(e25.15)') (alondin(m),m=1,nxin*nyin)
        read (iures,'(e25.15)') (alatdin(m),m=1,nxin*nyin)
      else
#if defined (LONLAT)
#  if defined (SAMEDOMAIN)
        if (.true.) then
#  else
        if (nxin.eq.nx .and. nyin.eq.ny) then
#  endif
c         assume same domain size as model:
          do j=1,nyin
            do i=1,nxin
              m = (j-1)*nxin + i
              alondin(m)= (alon1 + (alon2-alon1)*(i-.5)/nxin) 
              alatdin(m)= (alat1 + (alat2-alat1)*(j-.5)/nyin) 
            enddo
          enddo
        else
          write(ioterm,'(/a)')
     *      '*** Error: old restart file, different lon-lat grid'
          stop
        endif
#else
c       not necessary: not used unless LONLAT:
        call zero (alondin, nxin*nyin)
        call zero (alatdin, nxin*nyin)
#endif
      endif

c        Read 2-D and 3-D fields

      if (nint(versresin).le.7) then
        read (iures,'(e25.15)') ((zdum,i=1,nxin),j=1,nyin) ! xhin
        read (iures,'(e25.15)') ((zdum,i=1,nxin),j=1,nyin) ! yhin
      endif

      call readin  (h,        'h',       1,nx,1,ny,   1)
      call readin  (hs,       'hs',      1,nx,1,ny,   1)
      call readin  (hb,       'hb',      1,nx,1,ny,   1)
      if (nint(versresin).eq.6 .or. nint(versresin).eq.7)
     *  call readin  (work,   'hbsd',    1,nx,1,ny,   1)

c     ice temperature t equiv to tracer(...,1):
      do n=1,ntracein
        if (nint(versresin).le.16) then
          cvar = 'tracer'
        else
          write (cvar,'(a,i2.2)') 'tracer',n
        endif
        call readin (tracer(1,1,0,n), cvar, 1,nx,1,ny, nlevp+1)
      enddo
      if (ntracein.lt.ntrace) then
        do n=ntracein+1,ntrace
          call zero (tracer(1,1,0,n), nx*ny*(nlevp+1))
        enddo
      else if (ntrace.lt.ntracein) then
        write (ioterm,'(/a/a,i4,a,i4)')
     *    '*** Error: program ntrace < restart-file ntracein.',
     *    '    ntrace=',ntrace, '   ntracein=',ntracein
        stop
      endif

      call readin  (baseperc, 'baseperc',1,nx,1,ny,   1)
      call readin  (basefrml, 'basefrml',1,nx,1,ny,   1)
      call readin  (heatb,    'heatb',   1,nx,1,ny,   1)
      if (nint(versresin).eq.14) then
        call readin (work,'budgruna',1,nx,1,ny,   1)
        read (iures,'(e25.15)') zdum ! timruna
      endif
      call readin  (sedim,    'sedim',   1,nx,1,ny,   1)
      if (nint(versresin).le.7) 
     *  call readin  (work,   'sedimeq', 1,nx,1,ny,   1)
      if (nint(versresin).le.20)
     *  call readin  (work,   'sedimun', 1,nx,1,ny,   1)
      call readin  (tsed,     'tsed',    1,nx,1,ny,nsed)
      call readin  (wsed,     'wsed',    1,nx,1,ny,nsed)
      call readin  (heats,    'heats',   1,nx,1,ny,nsed)
      if (nint(versresin).ge.22) 
     *  call readin  (quarryacc,'quarryacc',1,nx,1,ny,1)
      if (nint(versresin).ge.24) 
     *  call readin  (hbinit, 'hbinit',  1,nx,1,ny,   1)
      call readin  (topbed,   'topbed',  1,nx,1,ny,   1)
      if (nint(versresin).le.7) 
     *  call readin (work,    'topbedeq',1,nx,1,ny,   1)
      if (nint(versresin).ge.23)  then
        call readin (topbedeq,'topbedeq',1,nx,1,ny,   1)
        call readin (equiload,'equiload',1,nx,1,ny,   1)
      endif
      call readin  (tbed,     'tbed',    1,nx,1,ny,nbed)
      if (nint(versresin).le.7) 
     *  call readin  (work,   'equiload',1,nx,1,ny,   1)
      call readin  (hw,       'hw',      1,nx,1,ny,   1)
      call readin  (tw,       'tw',      1,nx,1,ny,   1)
      if (nint(versresin).le.4)
     *  call readini (iwork,  'maskpres',1,nx,1,ny,   1)
      if (nint(versresin).le.7 .or. nint(versresin).ge.12) 
     *  call readin  (sedpres,'sedpres', 1,nx,1,ny,   1)
      if (nint(versresin).ge.7)
     *  call readin (crhmel,  'crhmel',  1,nx,1,ny,   1)
      if (nint(versresin).ge.9)
     *  call readin (enhanceinv,'enhanceinv',1,nx,1,ny,1)
      if (nint(versresin).ge.13 .and. nint(versresin).le.18) then
        call readin (work,     'oceanmeltav',  1,nx,1,ny,1)
        call readin (work,     'facemeltav',   1,nx,1,ny,1)
        if (nint(versresin).ge.15) then
          call readin (work,   'clifmeltav',   1,nx,1,ny,1)
          call readin (work,   'calviceav',    1,nx,1,ny,1)
          call readin (work,   'runliqav',     1,nx,1,ny,1)
          call readin (work,   'runfrozav',    1,nx,1,ny,1)
        else
          call readin (work,   'clifmelt0av',  1,nx,1,ny,1)
          call readin (work,   'clifmelt1av',  1,nx,1,ny,1)
          call readin (work,   'calviceav',    1,nx,1,ny,1)
        endif
      endif
      call readin  (ua,       'ua',      0,nxp,0,nyp, 1)
      call readin  (va,       'va',      0,nxp,0,nyp, 1)
      call readin  (ub,       'ub',      0,nxp,0,nyp, 1)
      call readin  (vb,       'vb',      0,nxp,0,nyp, 1)
      call readin  (dfu,      'dfu',     0,nxp,0,nyp, 1)
      call readin  (dfv,      'dfv',     0,nxp,0,nyp, 1)
      if (nint(versresin).ge.4) then
        call readin  (dbu,    'dbu',     0,nxp,0,nyp, 1)
        call readin  (dbv,    'dbv',     0,nxp,0,nyp, 1)
      endif
      if (nint(versresin).eq.7) then
        call readin  (worku,  'ubal',    0,nxp,0,nyp, 1)
        call readin  (worku,  'vbal',    0,nxp,0,nyp, 1)
        call readin  (work,   'hsbal',   1,nx,1,ny,   1)
        call readin  (work,   'hbal',    1,nx,1,ny,   1)
      endif
#if defined (MELANGE)
      if (nint(versresin).ge.20) then
        call readin  (hmel,   'hmel',    1,nx,1,ny,   1)
      else
        call zero (hmel, nx*ny)
      endif
#else
      if (nint(versresin).ge.20) then
        call readin  (work,   'hmel',    1,nx,1,ny,   1)
      endif
#endif

      if (nint(versresin).ge.8) read (iures,*)
      read (iures, '(e25.15)') sealev

      if (nint(versresin).ge.3) then
        if (nint(versresin).ge.8) read (iures,*)
        read (iures, '(e25.15)') rco2inter
      endif

c     sedtrack,itrtop in common (used in sedtracking_eul):
      call readin  (sedtrack, 'sedtrack',1,nx,1,ny,ntrack+1)
      if (nint(versresin).ge.8) read (iures,*)
      read (iures, '(e25.15)') zdum              ! old timebot, not used
      call readini (itrtop,   'itrtop',  1,nx,1,ny,   1)

      if (nint(versresin).ge.8) read (iures,*)
      read (iures, '(e25.15)') versresend
      close (iures)

      if (abs(versresin-versresend).gt.1.e-6) then
        write (ioterm,'(a,f8.2)') '   versresend=', versresend
        write (ioterm,'(a)')'*** Error: versresin and versresend differ'
        stop
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine readin (arr,  cvarin, nx1, nx2, ny1, ny2, nvert)
      entry     readini (iarr, cvarin, nx1, nx2, ny1, ny2, nvert) 

c     Read in one field (all horiz slices) from restart file.
c     If different resolution, process one horiz slice at a time:
c       if h-grid (nx1=1, calls interph[i],
c         else set to zero (u,v grids, nx1=0).
c       if real (all except 'irtrop'), bilinear interp for arr,
c         else integer, nearest neighbor for iarr.
c       if 'crhmel', logarithmic interp for arr, else linear.

#include <comicegrid.h>

      dimension  arr(nx1:nx2, ny1:ny2, nvert),
     *          iarr(nx1:nx2, ny1:ny2, nvert)
      character*(*) cvarin
      character*10 cvar, cvarfile
      dimension work(ninmx), iwork(ninmx)

#if defined (LONLAT)
      logical firstlonlat
      data firstlonlat /.true./
      save firstlonlat
      dimension zlonin(5000), zlatin(5000)               ! for setinterp
      dimension                                            
     *  indl(nx,ny),     weil(nx,ny),
     *  indb(nx,ny),     weib(nx,ny)
      save indl, weil, indb, weib, ifwrap
#elif defined (STEREO)
      dimension indh(nx,ny), jndh(nx,ny), weih(nx,ny,2)    ! for interph
      save indh, jndh, weih
#endif

c     write (ioterm,'(2a)') 'readin: reading ',cvarin

      cvar = cvarin

      if (nx1.eq.1) then
        ifhgrid = 1
      else
        ifhgrid = 0
      endif

      if (cvar.eq.'itrtop') then
        ifreal = 0
      else
        ifreal = 1
      endif

      if (cvar.eq.'crhmel') then
        iflog = 1
      else
        iflog = 0
      endif

      if (cvar.eq.'h') then
        iffirst = 1
      else
        iffirst = 0
      endif

c        Read field header record, check synch (versresin in common)

      if (nint(versresin).ge.8) then
        read (iures,'(a10,5x,i5)') cvarfile, nvertfile 
        if (cvarfile.ne.cvar .or. nvertfile.ne.nvert) then
          write (ioterm,'(/a,2(/a,a10,5x,i5))') 
     *      '*** Error reading restart file: incorrect field header',
     *      '    expected: ', cvar,     nvert,  
     *      '    file    : ', cvarfile, nvertfile  
          stop
        endif
      endif

c        If same resolution, just read into arr or iarr, and return.
c        (nxin,nyin are in common, have been read from restart file 
c        header lines in readres). 

      if ( nxin.eq.nx  .and. nyin.eq.ny  .and.   
     *     abs(dxin-dx0).lt.1.e-3 .and. abs(dyin-dy0).lt.1.e-3 ) then
        do k=1,nvert
          if (ifreal.eq.1) then
            read (iures,'(e25.15)') (( arr(i,j,k),i=nx1,nx2), j=ny1,ny2)
          else 
            read (iures,'(i10)') ((iarr(i,j,k),i=nx1,nx2), j=ny1,ny2)
          endif
        enddo
        return
      endif

c        Read in each horizontal slice for this field from 
c        restart file into [i]work, transfer to [i]arr if h-grid.
c        (nxin, nyin are in common, used in interph).

c-----------------
      do k=1,nvert 
c-----------------

c%%%%%%%%%%%%%%%%%%%
#if defined (LONLAT)
c%%%%%%%%%%%%%%%%%%%

        if (firstlonlat) then
          do i=1,nxin
            zlonin(i) = alondin(i + nxin*((nyin+1)/2-1)) * pi/180.
          enddo
          do j=1,nyin
            zlatin(j) = alatdin((nxin+1)/2 + nxin*(j-1)) * pi/180.
          enddo
          call setinterp (zlonin, zlatin, nxin, nyin,
     *                    alond, alatd, nx, ny,
     *                    indl, weil, indb, weib, 0.)

c         decide if input grid is global in longitude (so wraparound):
          dlonin = (zlonin(nxin)-zlonin(1)) / (nxin-1)
          ifwrap = 0
          if (dlonin*nxin. ge. 0.999*2.*pi) ifwrap = 1

          firstlonlat = .false.
        endif

c.............................
        if (ifhgrid.eq.1) then
c.............................

          if (ifreal.eq.1) then
            read (iures,'(e25.15)') (work(m),m=1,nxin*nyin)
          else
            read (iures,'(i10)') (iwork(m),m=1,nxin*nyin)
          endif

          do j=1,ny
            do i=1,nx
              il = indl(i,j)
              if (ifwrap.eq.0) then
                ir = min (il+1,nxin)    ! walls
              else
                ir = mod (il,nxin) + 1  ! wraparound
              endif
              jb = indb(i,j)
              jt = min (jb+1, nyin)
              wlb =     weil(i,j) *    weib(i,j)
              wrb = (1.-weil(i,j))*    weib(i,j)
              wlt =     weil(i,j) *(1.-weib(i,j))
              wrt = (1.-weil(i,j))*(1.-weib(i,j))
  
              if (ifreal.eq.1 .and. iflog.eq.0) then

                arr(i,j,k) = wlb * work(il + (jb-1)*nxin)
     *                     + wrb * work(ir + (jb-1)*nxin)
     *                     + wlt * work(il + (jt-1)*nxin)
     *                     + wrt * work(ir + (jt-1)*nxin)

              else if (ifreal.eq.1 .and. iflog.eq.1) then

                arr(i,j,k) = 10. ** (wlb*log10(work(il + (jb-1)*nxin))
     *                             + wrb*log10(work(ir + (jb-1)*nxin))
     *                             + wlt*log10(work(il + (jt-1)*nxin))
     *                             + wrt*log10(work(ir + (jt-1)*nxin)))

              else 

c               nearest neighbor, integer:
                weimax = max (wlb,wrb,wlt,wrt)
                if (weimax.eq.wlb) then
                  iarr(i,j,k) = iwork(il + (jb-1)*nxin)
                else if (weimax.eq.wrb) then
                  iarr(i,j,k) = iwork(ir + (jb-1)*nxin)
                else if (weimax.eq.wlt) then
                  iarr(i,j,k) = iwork(il + (jt-1)*nxin)
                else if (weimax.eq.wrt) then
                  iarr(i,j,k) = iwork(ir + (jt-1)*nxin)
                endif   

              endif   

            enddo   
          enddo   

c...........
        else
c...........

          if (ifreal.eq.1) then
            read (iures,'(e25.15)') (work(m),m=1,(nxin+2)*(nyin+2))
            call zero (arr(0,0,k), (nxp+1)*(nyp+1))
          else
            read (iures,'(i10)') (iwork(m),m=1,(nxin+2)*(nyin+2))
            call izero (iarr(0,0,k), (nxp+1)*(nyp+1))
          endif

c............
        endif
c............

c%%%%%%%%%%%%%%%%%%%%%
#elif defined (STEREO) 
c%%%%%%%%%%%%%%%%%%%%%

        if (ifreal.eq.1) then

          if (ifhgrid.eq.1) then
            read (iures,'(e25.15)') (work(m),m=1,nxin*nyin)
c           bilinear interp:
            call interph (arr(1,1,k), xh, yh, nx, ny,
     *                    work, xhin, yhin, nxin, nyin, 
     *                    indh, jndh, weih,
     *                    ifreal, iflog, iffirst)
          else
            read (iures,'(e25.15)') (work(m),m=1,(nxin+2)*(nyin+2))
            call zero (arr(0,0,k), (nxp+1)*(nyp+1))
          endif

        else

          if (ifhgrid.eq.1) then
            read (iures,'(i10)') (iwork(m),m=1,nxin*nyin)
c           nearest neighbor:
            call interphi (iarr(1,1,k), xh, yh, nx, ny,
     *                     iwork, xhin, yhin, nxin, nyin,
     *                     indh, jndh, weih,
     *                     ifreal, iflog, iffirst)
          else
            read (iures,'(i10)') (iwork(m),m=1,(nxin+2)*(nyin+2))
            call izero (iarr(0,0,k), (nxp+1)*(nyp+1))
          endif

        endif
   
c%%%%%
#endif
c%%%%%

c----------
      enddo
c----------

c     if (cvarin.eq.'h') then
c       write (250,'(i6,f10.3)') (i,arr(i,1,1),i=1,nx)
c       write (251,'(i6,f10.3)') (i,work(i),i=1,nxin)
c     else if (cvarin.eq.'hb') then
c       write (252,'(i6,f10.3)') (i,arr(i,1,1),i=1,nx)
c       write (253,'(i6,f10.3)') (i,arr(i,1,1),i=1,nx)
c     endif
 
c     write (ioterm,'(2a)') 'readin: read ok ',cvar

      return
      end

c-----------------------------------------------------------------------

      subroutine interph  (arr,    xh,    yh,   nx,   ny,
     *                     arrin,  xhin,  yhin, nxin, nyin,
     *                     indh,   jndh,  weih,
     *                     ifreal, iflog, iffirst)

      entry      interphi (iarr,   xh,    yh,   nx,   ny,
     *                     iarrin, xhin,  yhin, nxin, nyin,
     *                     indh,   jndh,  weih,
     *                     ifreal, iflog, iffirst)

c        Horizontal interpolation or nearest neigbor for one horizontal
c        slice, on h-grid, from arrin(nxin,nyin) to arr(nx,ny).
c        Only works for pole-centered stereo grids, or for rectangular
c        grids. Latter can be unevely spaced.
c
c        If iffirst=1, compute and save bilinear interpolation 
c          indices indh,jndh and weights weih (for h-grid).
c        If ifreal=1, arrays are real, do bilinear interp.
c          else arrays are integer, do nearest neighbor.
c        If iflog=1 (and ifreal=1), do logarithmic weighting 
c          (for crhmel), else linear.

      dimension  
     *  arr(nx,ny),  arrin(nxin,nyin),
     *  iarr(nx,ny), iarrin(nxin,nyin),
     *  xh(nx,ny),   xhin(nxin*nyin),
     *  yh(nx,ny),   yhin(nxin*nyin),
     *  indh(nx,ny), jndh(nx,ny), weih(nx,ny,2)

c     statement fn for 2D (nxin,nyin) to 1D xhin, yhin:
      m1d(i,j) = i + (j-1)*nxin 

c        If iffirst=1, set bilinear-interp indices and weights, 
c        from model grid into restart-file grid 

c---------------------------
      if (iffirst.eq.1) then
c---------------------------

        do j=1,ny
          do i=1,nx

            jin = (nyin+1)/2
            if (xh(i,j).le.xhin(m1d(1,jin))) then
              indh(i,j) = 1
              weih(i,j,1) = 0.
            else if (xh(i,j).ge.xhin(m1d(nxin,jin))) then
              indh(i,j) = nxin
              weih(i,j,1) = 0.
            else 
              do m = 2,nxin
                if (xh(i,j).le.xhin(m1d(m,jin))) then
                  indh(i,j) = m-1
                  weih(i,j,1) = (xh(i,j)          - xhin(m1d(m-1,jin)))
     *                        / (xhin(m1d(m,jin)) - xhin(m1d(m-1,jin)))
                  goto 10
                endif
              enddo
   10         continue
            endif

            iin = (nxin+1)/2
            if (yh(i,j).le.yhin(m1d(iin,1))) then
              jndh(i,j) = 1
              weih(i,j,2) = 0.
            else if (yh(i,j).ge.yhin(m1d(iin,nyin))) then
              jndh(i,j) = nyin
              weih(i,j,2) = 0.
            else 
              do m = 2,nyin
                if (yh(i,j).le.yhin(m1d(iin,m))) then
                  jndh(i,j) = m-1
                  weih(i,j,2) = (yh(i,j)          - yhin(m1d(iin,m-1)))
     *                        / (yhin(m1d(iin,m)) - yhin(m1d(iin,m-1)))
                  goto 20
                endif
              enddo
   20         continue
            endif
          enddo
        enddo

c----------
      endif
c----------

c        Interpolate to model grid

c================
      do j=1,ny
        do i=1,nx
c================

          ja = jndh(i,j)
          jb = min (nyin, jndh(i,j) + 1)
          ia = indh(i,j)
          ib = min (nxin, indh(i,j) + 1)

          zw1 = (1.-weih(i,j,1)) * (1.-weih(i,j,2))
          zw2 = (   weih(i,j,1)) * (1.-weih(i,j,2))
          zw3 = (1.-weih(i,j,1)) * (   weih(i,j,2))
          zw4 = (   weih(i,j,1)) * (   weih(i,j,2))

          if (ifreal.eq.1) then

c           bilinear interp, real:
            if (iflog.eq.0) then
              arr(i,j) = zw1*arrin(ia,ja)
     *                 + zw2*arrin(ib,ja)
     *                 + zw3*arrin(ia,jb)
     *                 + zw4*arrin(ib,jb)
            else
              arr(i,j) = 10. ** (   zw1*log10(arrin(ia,ja))
     *                            + zw2*log10(arrin(ib,ja))
     *                            + zw3*log10(arrin(ia,jb))
     *                            + zw4*log10(arrin(ib,jb)) )
            endif

          else 

c           nearest neighbor, integer:
            weimax = max (zw1,zw2,zw3,zw4)
            if (weimax.eq.zw1) then
              iarr(i,j) = iarrin (ia,ja)
            else if (weimax.eq.zw2) then
              iarr(i,j) = iarrin (ib,ja)
            else if (weimax.eq.zw3) then
              iarr(i,j) = iarrin (ia,jb)
            else if (weimax.eq.zw4) then
              iarr(i,j) = iarrin (ib,jb)
            endif

          endif

c============
        enddo
      enddo
c============

      return
      end

c-----------------------------------------------------------------------

      subroutine writeres (h, hs, hb, t, tracer,
     *                     baseperc, basefrml, heatb, 
     *                     sedim, tsed, wsed, heats, quarryacc, hbinit,
     *                     topbed, topbedeq, equiload, tbed, hw, tw, 
     *                     sedpres, crhmel, enhanceinv, 
     *                     ua, va, ub, vb, dfu, dfv, dbu, dbv,
#if defined (MELANGE)
     *                     hmel,
#endif
     *                     sealev, rco2inter, timeice, iczoom,
     *                     crestin, iflag)

c     Write restart file to filename crestin. 
c     If iflag=1, append year number to name.

#include <comicephys.h>
#include <comicegrid.h>

c     h grid:
      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny),
     *  t(nx,ny,0:nlevp),    tracer(nx,ny,0:nlevp,ntrace),
     *  baseperc(nx,ny),     basefrml(nx,ny),    heatb(nx,ny),        
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     *  heats(nx,ny,nsed),   quarryacc(nx,ny),   hbinit(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny),    equiload(nx,ny),
     *  tbed(nx,ny,nbed),    hw(nx,ny),          tw(nx,ny),
     *  sedpres(nx,ny),      crhmel(nx,ny),      enhanceinv(nx,ny)

c     u grid:
      dimension
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp),
     *  dbu(0:nxp,0:nyp),       dbv(0:nxp,0:nyp)

#if defined (MELANGE)
     * ,hmel(nx,ny)
#endif

      character*(*) crestin

      dimension work(nx,ny)
      character crest*80, cyr*80, cvar*10

      crest = crestin

c        If requested (iflag=1) , append year number to filename
   
      if (iflag.eq.1) then
        write (cyr, '(i10)') nint(abs(timeice))
        do i=1,lenchr(cyr)
          if (cyr(1:1).eq.' ') cyr = cyr(2:) 
        enddo   
c       if (timeice.lt.0.) cyr = cyr(1:lenchr(cyr)) // 'm'
        if (timeice.gt.0.) cyr = cyr(1:lenchr(cyr)) // 'p'
        crest = crest(1:lenchr(crest)) //  '_' // cyr
      endif

      write (ioterm,'(/a,a)') 
     *  'Writing restart file ', crest(1:lenchr(crest))

      open (iures, file=crest,status='unknown',form='formatted')
      close (iures, status='delete')
      open (iures, file=crest,status='new',form='formatted')
c     versres = 3. !           added rco2inter
c     versres = 4. ! 10/25/08, added dbu,dbv
c     versres = 5. ! 5/19/10,  removed maskpres
c     versres = 6. ! 7/08/10,  added hbsd
c     versres = 7. ! 11/04/10, added crhmel, ubal,vbal, hsbal,hbal
c     versres = 8. ! 11/23/10, added subheaders,
c                  !           removed xhin,yhin,hbsd,sedimeq,topbedeq,
c                  !                   equiload,sedpres,
c                  !                   ubal,vbal, hsbal,hbal
c     versres = 9. ! 12/28/10, added enhanceinv
c     versres = 10.! 1/5/11,   removed nwrit
c     versres = 11.! 1/23/12,  added xhin, yhin
c     versres = 12.! 6/6/13,   (re)added sedpres
c     versres = 13.! 7/27/13,  added oceanmeltav,facemeltav, 
c                  !                 clifmelt[0,1]av, calviceav
c     versres = 14.! 9/19/13,  added budgruna, timruna
c     versres = 15.! 9/19/13,  removed budgruna,timruna,clifmelt[0,1]av,
c                  !           added runliqav, runfrozav
c     versres = 16.! 11/13/13, added alond, alatd
c     versres = 17.! 5/6/15,   added ntrace
c     versres = 18.! 10/17/15, added xoffa, yoffa (although not used)
c     versres = 19.! 6/14/16,  removed oceanmeltav,facemeltav,
c                  !                   cliffmeltav, calviceav,
c                  !                   runliqav, runfrozav
c     versres = 20.! 7/14/17   added hmel (for melange)
c     versres = 21.! 7/17/18   removed sedimun
c     versres = 22.! 9/2/18    added quarryacc
c     versres = 23.! 10/22/18  added topbedeq, equiload (for 34myr)
      versres = 24.! 11/3/18   added hbinit (for 34myr quarryread/write)

      write (iures,'(e25.15)') versres
      write (iures,'(e25.15,6i10)') timeice,nx,ny,nlev,nsed,nbed,ntrace
      write (iures,'(4e25.15)') dx0, dy0, xoffa, yoffa
      write (iures,'(e25.15)') xh
      write (iures,'(e25.15)') yh
      write (iures,'(e25.15)') alond
      write (iures,'(e25.15)') alatd
#if defined (EISANTA) && defined (TRANSECTA) && defined (ZOOMGL)
      write (iures) iczoom
#endif

c       Write 2-D and 3-D fields

      call writeout  (h,        'h',       1,nx,1,ny,   1)
      call writeout  (hs,       'hs',      1,nx,1,ny,   1)
      call writeout  (hb,       'hb',      1,nx,1,ny,   1)
c     ice temperature t equiv to tracer(...,1):
      do n=1,ntrace
        write (cvar,'(a,i2.2)') 'tracer',n
        call writeout (tracer(1,1,0,n), cvar, 1,nx,1,ny, nlevp+1)
      enddo
      call writeout  (baseperc, 'baseperc',1,nx,1,ny,   1)
      call writeout  (basefrml, 'basefrml',1,nx,1,ny,   1)
      call writeout  (heatb,    'heatb',   1,nx,1,ny,   1)
      call writeout  (sedim,    'sedim',   1,nx,1,ny,   1)
      call writeout  (tsed,     'tsed',    1,nx,1,ny,nsed)
      call writeout  (wsed,     'wsed',    1,nx,1,ny,nsed)
      call writeout  (heats,    'heats',   1,nx,1,ny,nsed)
      call writeout  (quarryacc, 'quarryacc',1,nx,1,ny, 1)
      call writeout  (hbinit,   'hbinit',  1,nx,1,ny,   1)
      call writeout  (topbed,   'topbed',  1,nx,1,ny,   1)
      call writeout  (topbedeq, 'topbedeq',1,nx,1,ny,   1)
      call writeout  (equiload, 'equiload',1,nx,1,ny,   1)
      call writeout  (tbed,     'tbed',    1,nx,1,ny,nbed)
      call writeout  (hw,       'hw',      1,nx,1,ny,   1)
      call writeout  (tw,       'tw',      1,nx,1,ny,   1)
      call writeout  (sedpres,  'sedpres', 1,nx,1,ny,   1)
      call writeout  (crhmel,   'crhmel',  1,nx,1,ny,   1)
      call writeout  (enhanceinv,'enhanceinv',1,nx,1,ny,1)
      call writeout  (ua,       'ua',      0,nxp,0,nyp, 1)
      call writeout  (va,       'va',      0,nxp,0,nyp, 1)
      call writeout  (ub,       'ub',      0,nxp,0,nyp, 1)
      call writeout  (vb,       'vb',      0,nxp,0,nyp, 1)
      call writeout  (dfu,      'dfu',     0,nxp,0,nyp, 1)
      call writeout  (dfv,      'dfv',     0,nxp,0,nyp, 1)
      call writeout  (dbu,      'dbu',     0,nxp,0,nyp, 1)
      call writeout  (dbv,      'dbv',     0,nxp,0,nyp, 1)
#if defined (MELANGE)
      call writeout  (hmel,     'hmel',    1,nx,1,ny,   1)
#else
      call zero (work, nx*ny)
      call writeout  (work,     'hmel',    1,nx,1,ny,   1)
#endif

      write (iures, '(a/e25.15)') 'sealev:', sealev
      write (iures, '(a/e25.15)') 'rco2inter:', rco2inter

c     sedtrack,itrtop in common (used in sedtracking_eul):
      call writeout  (sedtrack, 'sedtrack',1,nx,1,ny,ntrack+1)
      write (iures, '(a/e25.15)') 'timebot:', 0. ! old timebot, not used
      call writeouti (itrtop, 'itrtop',    1,nx,1,ny,   1)

      write (iures, '(a/e25.15)') 'versres:', versres
      close (iures)

      return
      end

c-----------------------------------------------------------------------

      subroutine writeout (arr,  cvarin, nx1,nx2,ny1,ny2, nvert)
      entry     writeouti (iarr, cvarin, nx1,nx2,ny1,ny2, nvert)

c     Writes one field (all horiz slices) to restart file

#include <comicegrid.h>

      dimension  arr(nx1:nx2, ny1:ny2, nvert),
     *          iarr(nx1:nx2, ny1:ny2, nvert)
      character*(*) cvarin
      character*10 cvar

      cvar = cvarin

      if (cvar.eq.'itrtop') then
        ifreal = 0
      else
        ifreal = 1
      endif

      write (iures,'(a10,5x,i5)') cvar, nvert 

      do k=1,nvert
        if (ifreal.eq.1) then
          write (iures,'(e25.15)') (( arr(i,j,k),i=nx1,nx2), j=ny1,ny2)
        else 
          write (iures,'(i10)') ((iarr(i,j,k),i=nx1,nx2), j=ny1,ny2)
        endif
      enddo

      return
      end

c-----------------------------------------------------------------------

#if ! defined (NESTING)

      subroutine writenest (h, hb, hw, t, tbed, sedim, tsed, wsed, 
     *                      topbedeq, equiload, sedpres, arcocn, 
     *                      maskwater, ua, va, timeice)

c     Writes fields needed for subsequent nested run, for initializing 
c     and specifying edge b.c's. Called at end of continental-scale run.
c     Appends year number to basic file name ('nestdrive'). File will 
c     be read by getnest (readnest(2)) as 'nestdrivein' in later nested 
c     run. Loops twice, with iflag=0 and 1 to write to generic file,
c     and to file with year number appended to name.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),         hb(nx,ny),         hw(nx,ny),
     *  t(nx,ny,0:nlevp), tbed(nx,ny,nbed),
     *  sedim(nx,ny),     tsed(nx,ny,nsed),  wsed(nx,ny,nsed),
     *  topbedeq(nx,ny),  equiload(nx,ny),
     *  sedpres(nx,ny),   arcocn(nx,ny),     maskwater(nx,ny),
     *  ua(0:nxp,0:nyp),  va(0:nxp,0:nyp)

      dimension work(nx,ny)

      character cnestroot*80
      parameter (cnestroot = 'nestdrive')
      character cnest*80, cyr*80

c        Propagate arcocn to every point (using nearest ocean pt
c        at which it was calculated in oceanarc) 

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1 .and. arcocn(i,j).ne.360.) then
            work(i,j) = arcocn(i,j)
          else
            zdistmin = 1.e20
            do jj=1,ny
              do ii=1,nx
                if (maskwater(ii,jj).eq.1 .and. arcocn(ii,jj).ne.360.) 
     *            then
                  zdist =   (xh(i,j)-xh(ii,jj))**2 
     *                    + (yh(i,j)-yh(ii,jj))**2
                  if (zdist.lt.zdistmin) then
                    iimin = ii
                    jjmin = jj
                    zdistmin = zdist
                  endif
                endif
              enddo
            enddo
            work(i,j) = arcocn(iimin,jjmin)
          endif
        enddo
      enddo

c=================
      do iflag=0,1
c=================

c        If iflag=1, append year number to name

      if (iflag.eq.0) then
        cnest = cnestroot
      else
        write (cyr, '(i10)') nint(abs(timeice))
        do i=1,lenchr(cyr)
          if (cyr(1:1).eq.' ') cyr = cyr(2:)
        enddo
c       if (timeice.lt.0.) cyr = cyr(1:lenchr(cyr)) // 'm'
        if (timeice.gt.0.) cyr = cyr(1:lenchr(cyr)) // 'p'
        cnest = cnestroot(1:lenchr(cnestroot)) //  '_' // cyr
      endif

      write (ioterm,'(/a,a)') 'Writing nesting file ',cnest

      open (iunest, file=cnest,status='unknown',form='formatted')
      close (iunest, status='delete')
      open (iunest, file=cnest,status='new',form='formatted')

#if defined (LONLAT)
      write (iunest,'(a)') 'LONLAT'
#elif defined (STEREO)
      write (iunest,'(a)') 'STEREO'
#endif

      write (iunest,'(f15.3)') timeice
#if defined (LONLAT)
      write (iunest,'(5i6,4f10.4)') nx,ny,nlev,nbed,nsed,
     *                              alon1,alon2,alat1,alat2
#else
      write (iunest,'(5i6)') nx,ny,nlev,nbed,nsed
#endif
      write (iunest,'(e25.15)') h
      write (iunest,'(e25.15)') hb
      do k=0,nlevp
        write (iunest,'(e25.15)') ((t(i,j,k),i=1,nx),j=1,ny)
      enddo
      do k=1,nbed
        write (iunest,'(e25.15)') ((tbed(i,j,k),i=1,nx),j=1,ny)
      enddo

      write (iunest,'(e25.15)') sedim
      do k=1,nsed
        write (iunest,'(e25.15)') ((tsed(i,j,k),i=1,nx),j=1,ny)
      enddo
      do k=1,nsed
        write (iunest,'(e25.15)') ((wsed(i,j,k),i=1,nx),j=1,ny)
      enddo

      write (iunest,'(e25.15)') topbedeq
      write (iunest,'(e25.15)') equiload
      write (iunest,'(e25.15)') sedpres
      write (iunest,'(e25.15)') work                 ! propagated arcocn
      write (iunest,'(e25.15)') ua
      write (iunest,'(e25.15)') va

      write (iunest,'(f15.3)') timeice

      close (iunest)

c==========
      enddo  ! iflag
c==========

      return
      end

#endif

c-----------------------------------------------------------------------

c@@@@@@@@@@@@@@@@@@@@
#if defined (NESTING)
c@@@@@@@@@@@@@@@@@@@@

      subroutine getnest (h, hb, t, tbed, sedim, tsed, wsed, 
     *                    topbedeq, equiload, sedpres, arcocn,maskwater,
     *                    ua, va, iflagnest, iflagnesteq, iflagnesttemp,
     *                    timeice)
 
c     Calls readnest to get current set of continental-grid arrays 
c     (*in, written at end of a previous continental-scale run
c     by writenest). Interpolates needed quantities to nested grid.
c     Uses spatial interpolation indices set in setinterpnest
c     (inest, jnest, wnest, in comicegrid).
c
c     Called after init[hys settings from input datasets,
c     and if ifrest=0 or 1, also called after restart file is read. 
c     (And called at intervals during run for edge-resets only, 
c     if time-varying nest driving files are used).  
c
c     If iflagnest = 0 or 1, sets h,hb,ua,va variables over entire 
c     nested domain (0), or just along edges (1).
c
c     Similarly for topbedeq,equiload with iflagnesteq = 0, 1, 2 
c     (if = 2, not reset anywhere).
c
c     Similarly for temperatures if iflagnesttemp = 0 or 1. Call from
c     initphys has to set it to 0 (setting everywhere), in case
c     ifrest = -1, since temps not set by initial datasets.
c
#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny),         hb(nx,ny),
     *  t(nx,ny,0:nlevp), tbed(nx,ny,nbed),
     *  sedim(nx,ny),     tsed(nx,ny,nsed),  wsed(nx,ny,nsed),
     *  topbedeq(nx,ny),  equiload(nx,ny),
     *  sedpres(nx,ny),   arcocn(nx,ny),     maskwater(nx,ny),
     *  ua(0:nxp,0:nyp),  va(0:nxp,0:nyp)

      dimension 
     *  hin(nxa,nya),        hbin(nxa,nya),
     *  tin(nxa,nya,0:nlevp),tbedin(nxa,nya,nbed),
     *  sedimin(nxa,nya),    tsedin(nxa,nya,nsed), wsedin(nxa,nya,nsed),
     *  topbedeqin(nxa,nya), equiloadin(nxa,nya),
     *  sedpresin(nxa,nya),  arcocnin(nxa,nya),
     *  uain(0:nxa+1,0:nya+1), vain(0:nxa+1,0:nya+1) 

c    for continental-to-nested grid interpolation:
      dimension
     *  inest(nx,ny),           jnest(nx,ny),   
     *  wnest(nx,ny,2),
     *  inestu(0:nx+1,0:ny+1),  jnestu(0:nx+1,0:ny+1), 
     *  wnestu(0:nx+1,0:ny+1,2),
     *  inestv(0:nx+1,0:ny+1),  jnestv(0:nx+1,0:ny+1),
     *  wnestv(0:nx+1,0:ny+1,2)
      save inest,jnest,wnest, inestu,jnestu,wnestu, inestv,jnestv,wnestv

      logical firstreadnest
      save firstreadnest
      data firstreadnest /.true./

c        Read (+ time interp) continental arrays from nestdrive file(s)

      call readnest (hin, hbin, tin, tbedin, sedimin, tsedin, wsedin,
     *               topbedeqin, equiloadin,
     *               sedpresin, arcocnin, uain, vain, 
#if defined (LONLAT)
     *               alon1_in, alon2_in, alat1_in, alat2_in,
#endif
     *               timeice)

c        Get indices for continental-to-nested grid interpolation

      if (firstreadnest) then
        call setinterpnest (inest,  jnest,  wnest, 
     *                      inestu, jnestu, wnestu,
     *                      inestv, jnestv, wnestv
#if defined (LONLAT)
     *                     ,alon1_in, alon2_in, alat1_in, alat2_in
#endif
     *                     )
        firstreadnest = .false.
      endif

c        Spatially interpolate to nest-grid arrays

c>>>>>>>>>>>>>
c h-grid vars:
c>>>>>>>>>>>>>

      do j=1,ny
        do i=1,nx
          ja = jnest(i,j)
          jb = min (nya, jnest(i,j) + 1)
          ia = inest(i,j)
          ib = min (nxa, inest(i,j) + 1)

          zw1 = (1.-wnest(i,j,1)) * (1.-wnest(i,j,2))
          zw2 = (   wnest(i,j,1)) * (1.-wnest(i,j,2))
          zw3 = (1.-wnest(i,j,1)) * (   wnest(i,j,2))
          zw4 = (   wnest(i,j,1)) * (   wnest(i,j,2))

c--------------------------------------------------------------
          if ( ( iflagnest.eq.0 ) .or.
     *         ( iflagnest.eq.1 .and. (i.eq.1 .or. i.eq.nx .or.
     *                                 j.eq.1 .or. j.eq.ny) ) 
     *       ) then 
c--------------------------------------------------------------
            h(i,j) = zw1*hin(ia,ja) + zw2*hin(ib,ja)
     *             + zw3*hin(ia,jb) + zw4*hin(ib,jb)
            hb(i,j) = zw1*hbin(ia,ja) + zw2*hbin(ib,ja)
     *              + zw3*hbin(ia,jb) + zw4*hbin(ib,jb)
c--------------
          endif
c--------------

c..................................................................
          if ( ( iflagnesttemp.eq.0 ) .or.
     *         ( iflagnesttemp.eq.1 .and. (i.eq.1 .or. i.eq.nx .or.
     *                                 j.eq.1 .or. j.eq.ny) ) 
     *       ) then 
c..................................................................
          do k=0,nlevp
            t(i,j,k) = zw1*tin(ia,ja,k) + zw2*tin(ib,ja,k)
     *               + zw3*tin(ia,jb,k) + zw4*tin(ib,jb,k)
          enddo

          do k=1,nbed
            tbed(i,j,k) = zw1*tbedin(ia,ja,k) + zw2*tbedin(ib,ja,k)
     *                  + zw3*tbedin(ia,jb,k) + zw4*tbedin(ib,jb,k)
          enddo

          sedim(i,j) = zw1*sedimin(ia,ja) + zw2*sedimin(ib,ja)
     *               + zw3*sedimin(ia,jb) + zw4*sedimin(ib,jb)

          do k=1,nsed
            tsed(i,j,k) = zw1*tsedin(ia,ja,k) + zw2*tsedin(ib,ja,k)
     *                  + zw3*tsedin(ia,jb,k) + zw4*tsedin(ib,jb,k)
            wsed(i,j,k) = zw1*wsedin(ia,ja,k) + zw2*wsedin(ib,ja,k)
     *                  + zw3*wsedin(ia,jb,k) + zw4*wsedin(ib,jb,k)
          enddo
c..............
          endif
c..............

          sedpres(i,j) = zw1*sedpresin(ia,ja) + zw2*sedpresin(ib,ja)
     *                 + zw3*sedpresin(ia,jb) + zw4*sedpresin(ib,jb)
          sedpres(i,j) = nint(sedpres(i,j))

#if defined (KEEPARCOCN)
          arcocn(i,j) = zw1*arcocnin(ia,ja) + zw2*arcocnin(ib,ja)
     *                + zw3*arcocnin(ia,jb) + zw4*arcocnin(ib,jb)
#endif

c----------------------------------------------------------------
          if ( ( iflagnesteq.eq.0 ) .or.
     *         ( iflagnesteq.eq.1 .and. (i.eq.1 .or. i.eq.nx .or.
     *                                   j.eq.1 .or. j.eq.ny) ) 
     *       ) then 
c----------------------------------------------------------------
            topbedeq(i,j)= zw1*topbedeqin(ia,ja) + zw2*topbedeqin(ib,ja)
     *                   + zw3*topbedeqin(ia,jb) + zw4*topbedeqin(ib,jb)

            equiload(i,j)= zw1*equiloadin(ia,ja) + zw2*equiloadin(ib,ja)
     *                   + zw3*equiloadin(ia,jb) + zw4*equiloadin(ib,jb)
c--------------
          endif
c--------------

        enddo
      enddo

c>>>>>>>>>>>>>>>
c u,v-grid vars:
c>>>>>>>>>>>>>>>

      do j=0,nyp
        do i=0,nxp
c-------~~-------------------------------------------------------
          if ( ( iflagnest.eq.0 ) .or.
     *         ( iflagnest.eq.1 .and. (i.eq.0 .or. i.ge.nx .or.
     *                                 j.le.1 .or. j.ge.ny) ) 
     *       ) then 
c---------~~-----------------------------------------------------
            ja = jnestu(i,j)
            jb = min (nya+1, jnestu(i,j) + 1)
            ia = inestu(i,j)
            ib = min (nxa+1, inestu(i,j) + 1)

            zw1 = (1.-wnestu(i,j,1)) * (1.-wnestu(i,j,2))
            zw2 = (   wnestu(i,j,1)) * (1.-wnestu(i,j,2))
            zw3 = (1.-wnestu(i,j,1)) * (   wnestu(i,j,2))
            zw4 = (   wnestu(i,j,1)) * (   wnestu(i,j,2))

            ua(i,j) = zw1*uain(ia,ja) + zw2*uain(ib,ja)
     *              + zw3*uain(ia,jb) + zw4*uain(ib,jb)
c~~~~~~~~~~~~~~
          endif
c~~~~~~~~~~~~~~

c--------------------------------~~----------------------------
          if ( ( iflagnest.eq.0 ) .or.
     *         ( iflagnest.eq.1 .and. (i.le.1 .or. i.ge.nx .or.
     *                                 j.eq.0 .or. j.ge.ny) ) 
     *       ) then 
c-----------------------------------~~-------------------------
            ja = jnestv(i,j)
            jb = min (nya+1, jnestv(i,j) + 1)
            ia = inestv(i,j)
            ib = min (nxa+1, inestv(i,j) + 1)

            zw1 = (1.-wnestv(i,j,1)) * (1.-wnestv(i,j,2))
            zw2 = (   wnestv(i,j,1)) * (1.-wnestv(i,j,2))
            zw3 = (1.-wnestv(i,j,1)) * (   wnestv(i,j,2))
            zw4 = (   wnestv(i,j,1)) * (   wnestv(i,j,2))

            va(i,j) = zw1*vain(ia,ja) + zw2*vain(ib,ja)
     *              + zw3*vain(ia,jb) + zw4*vain(ib,jb)
c~~~~~~~~~~~~~~
          endif
c~~~~~~~~~~~~~~
        enddo
      enddo

#if defined (JAKOBOCN)
c     zero any floating ice at domain edge (or 1 cell in) 
c     for Jakobshavn nested runs:
      do j=1,ny
        do i=1,nx
          if ( (i.le.2 .or. i.ge.nx-1 .or. j.le.2 .or. j.ge.ny-1)
     *          .and. (maskwater(i,j).ne.0) ) then
            h(i,j) = 0.
            do k=0,nlevp
              t(i,j,k) = tmelt
            enddo
            if (i.le.2)    ua(i-1,j) = 0.
            if (i.ge.nx-1) ua(i,j)   = 0.
            if (j.le.2)    va(i,j-1) = 0.
            if (j.ge.ny-1) va(i,j)   = 0.
          endif
        enddo
      enddo
#endif

      return
      end

c-----------------------------------------------------------------------

      subroutine readnest (hin, hbin, tin, tbedin,sedimin,tsedin,wsedin,
     *                     topbedeqin, equiloadin,
     *                     sedpresin, arcocnin, uain, vain, 
#if defined (LONLAT)
     *                     alon1_in, alon2_in, alat1_in, alat2_in,
#endif
     *                     timeice)

c     Controls reading of nestdrive files, maintaining
c     2 sets of variables (*x, *y) from the 2 files whose times
c     surround current time. Linearly weights these in time
c     to returned set of variables (*in).

#include <comicephys.h>
#include <comicegrid.h>

c     passed:
      dimension 
     *  hin(nxa,nya),        hbin(nxa,nya),
     *  tin(nxa,nya,0:nlevp),tbedin(nxa,nya,nbed),
     *  sedimin(nxa,nya),    tsedin(nxa,nya,nsed), wsedin(nxa,nya,nsed),
     *  topbedeqin(nxa,nya), equiloadin(nxa,nya),
     *  sedpresin(nxa,nya),  arcocnin(nxa,nya),
     *  uain(0:nxa+1,0:nya+1), vain(0:nxa+1,0:nya+1) 

c     local (saved):
      dimension 
     *  hx(nxa,nya),          hbx(nxa,nya),
     *  tx(nxa,nya,0:nlevp),  tbedx(nxa,nya,nbed),
     *  sedimx(nxa,nya), tsedx(nxa,nya,nsed), wsedx(nxa,nya,nsed),
     *  topbedeqx(nxa,nya),   equiloadx(nxa,nya),
     *  sedpresx(nxa,nya),    arcocnx(nxa,nya),
     *  uax(0:nxa+1,0:nya+1), vax(0:nxa+1,0:nya+1),

     *  hy(nxa,nya),          hby(nxa,nya),
     *  ty(nxa,nya,0:nlevp),  tbedy(nxa,nya,nbed),
     *  sedimy(nxa,nya), tsedy(nxa,nya,nsed), wsedy(nxa,nya,nsed),
     *  topbedeqy(nxa,nya),   equiloady(nxa,nya),
     *  sedpresy(nxa,nya),    arcocny(nxa,nya),
     *  uay(0:nxa+1,0:nya+1), vay(0:nxa+1,0:nya+1)

      save hx,        hbx,
     *     tx,        tbedx,
     *     sedimx,    tsedx,     wsedx,
     *     topbedeqx, equiloadx,
     *     sedpresx,  arcocnx,
     *     uax,       vax,

     *     hy,        hby,
     *     ty,        tbedy,
     *     sedimy,    tsedy,     wsedy,
     *     topbedeqy, equiloady,
     *     sedpresy,  arcocny,
     *     uay,       vay

      logical firstreadnest 
      save firstreadnest
      data firstreadnest /.true./

      dimension tnestlist(1000)
      save tnestlist, mnestaprev, mnestbprev 
      data mnestaprev, mnestbprev / -999, -999/

c       If first call, read first lines of all nestdrive files,
c       save a list of their times (tnestlist(nnestlist))

      if (firstreadnest) then
        do m=1,nnestlist
          call readnest2a (m, tnestlist(m), timeice, 0) 
        enddo
        firstreadnest = .false.
      endif

c        Find indices of the 2 nest files spanning current time, 
c        and set weight for current time

      if (timeice.le.tnestlist(1) .or. nnestlist.eq.1) then
        mnesta = 1
        mnestb = 1
        weia   = 1.
      else if (timeice.ge.tnestlist(nnestlist)) then
        mnesta = nnestlist
        mnestb = nnestlist
        weia   = 1.
      else
        do m=2,nnestlist
          if (timeice.le.tnestlist(m)) then
            mnesta = m-1
            mnestb = m
            weia= (tnestlist(m)-timeice) / (tnestlist(m)-tnestlist(m-1))
            goto 100
          endif
        enddo
  100   continue
      endif

c       If 2 surrounding nest files have changed, read them into
c       saved arrays (*x, *y)

      if (mnesta.ne.mnestaprev .or. mnestb.ne.mnestbprev) then
        call readnest2 (hx, hbx, tx, tbedx, sedimx, tsedx, wsedx,
     *                  topbedeqx, equiloadx,
     *                  sedpresx, arcocnx, uax, vax,
#if defined (LONLAT)
     *                  alon1_in, alon2_in, alat1_in, alat2_in,
#endif
     *                  mnesta, timenestdum, timeice, 1)

        call readnest2 (hy, hby, ty, tbedy, sedimy, tsedy, wsedy,
     *                  topbedeqy, equiloady,
     *                  sedpresy, arcocny, uay, vay,
#if defined (LONLAT)
     *                  alon1_in, alon2_in, alat1_in, alat2_in,
#endif
     *                  mnestb, timenestdum, timeice, 1)

        mnestaprev = mnesta
        mnestbprev = mnestb
      endif

c        Set returned arrays (*in), linearly interpolating in time
c        between the 2 sets of saved arrays

      do j=1,nya
        do i=1,nxa
          hin(i,j)        = weia*hx(i,j)       +(1.-weia)*hy(i,j)
          hbin(i,j)       = weia*hbx(i,j)      +(1.-weia)*hby(i,j)
          sedimin(i,j)    = weia*sedimx(i,j)   +(1.-weia)*sedimy(i,j)
          topbedeqin(i,j) = weia*topbedeqx(i,j)+(1.-weia)*topbedeqy(i,j)
          equiloadin(i,j) = weia*equiloadx(i,j)+(1.-weia)*equiloady(i,j)
          sedpresin(i,j)  = weia*sedpresx(i,j) +(1.-weia)*sedpresy(i,j)
          arcocnin(i,j)   = weia*arcocnx(i,j)  +(1.-weia)*arcocny(i,j)
          do k=0,nlevp
            tin(i,j,k)    = weia*tx(i,j,k)     +(1.-weia)*ty(i,j,k)
          enddo
          do k=1,nbed
            tbedin(i,j,k) = weia*tbedx(i,j,k)  +(1.-weia)*tbedy(i,j,k)
          enddo
          do k=1,nsed
            tsedin(i,j,k) = weia*tsedx(i,j,k)  +(1.-weia)*tsedy(i,j,k)
            wsedin(i,j,k) = weia*wsedx(i,j,k)  +(1.-weia)*wsedy(i,j,k)
          enddo
        enddo
      enddo

      do j=0,nya+1
        do i=0,nxa+1
          uain(i,j) = weia*uax(i,j) + (1.-weia)*uay(i,j)
          vain(i,j) = weia*vax(i,j) + (1.-weia)*vay(i,j)
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine readnest2 (hin, hbin, tin,tbedin,sedimin,tsedin,wsedin,
     *                     topbedeqin, equiloadin,
     *                     sedpresin, arcocnin, uain, vain,
#if defined (LONLAT)
     *                     alon1_in, alon2_in, alat1_in, alat2_in,
#endif
     *                     mnest, timenest, timeice, iflag)    ! iflag=1

      entry readnest2a (mnest, timenest, timeice, iflag)       ! iflag=0

c     Reads one nestdrive file, # mnest in Namelist list 
c     (cnestlist, in comicegrd). If iflag=0 (entry readnest2a),
c     only read first line for time (timenest) (and second line to check
c     dimensions). Else, read full set of variables. Opens and closes 
c     the file each call.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  hin(nxa,nya),        hbin(nxa,nya),
     *  tin(nxa,nya,0:nlevp),tbedin(nxa,nya,nbed),
     *  sedimin(nxa,nya),    tsedin(nxa,nya,nsed), wsedin(nxa,nya,nsed),
     *  topbedeqin(nxa,nya), equiloadin(nxa,nya),
     *  sedpresin(nxa,nya),  arcocnin(nxa,nya),
     *  uain(0:nxa+1,0:nya+1), vain(0:nxa+1,0:nya+1)
#if defined (OLDNESTFILE)
      dimension bmassin(nxa,nya)                         ! no more bmass
#endif

      character*80 cfilin, cgridin
      
      if (cnestpath.eq.' ') then
        cfilin = cnestlist(mnest)
      else
        cfilin = cnestpath(1:lenchr(cnestpath)) // cnestlist(mnest)
      endif

      write (ioterm,'(2a, a,i4, a,i2)') 
     *  'reading nestdrive file ', cfilin(1:lenchr(cfilin)),
     *  '   mnest=',mnest, '   iflag=',iflag

      open (iunest, file=cfilin, status='old', form='formatted',
     *      iostat=ios, err=1000)

      read (iunest,'(a)') cgridin
      if (cgridin.ne.'LONLAT' .and. cgridin.ne.'STEREO') then
        cgridin = 'STEREO'
        backspace iunest
      endif

#if defined (LONLAT)
      if (cgridin.eq.'STEREO') then
#elif defined (STEREO)
      if (cgridin.eq.'LONLAT') then
#endif
        write (ioterm,'(/a)')
     *    '*** Error (readnest2): grid type different from model'
        stop
      endif

      read (iunest,'(f15.3)') timenest
      write (ioterm,'(2(a,f12.2))') 
     *  '   timenest=',timenest, '   timeice=',timeice

c        Return if just want time for this file

      if (iflag.eq.0) then
        close (iunest)
        return
      endif

#if defined (LONLAT)
      read (iunest,'(5i6,4f10.4)') 
     *  nxa_in, nya_in, nlev_in, nbed_in, nsed_in,
     *  alon1_in, alon2_in, alat1_in, alat2_in
#elif defined (STEREO)
      read (iunest,'(5i6)') nxa_in, nya_in, nlev_in, nbed_in, nsed_in
#endif

      if (nxa.ne.nxa_in .or. nya.ne.nya_in .or. nlev.ne.nlev_in .or.
     *    nbed.ne.nbed_in .or. nsed.ne.nsed_in) then
        write (ioterm,'(/a,2(/a,i6,a,i6,a,i6,a,i6,a,i6))')
     *    '*** Error (readnest2): mismatched resolutions',
     *    '    Expected: nxa=',nxa,'  nya=',nya,'  nlev=',nlev,
     *                   '  nbed=',nbed,'  nsed=',nsed,
     *    '    File    : nxa=',nxa_in,'  nya=',nya_in,'  nlev=',nlev_in,
     *                   '  nbed=',nbed_in,'  nsed=',nsed_in
        stop
      endif
  
      read (iunest,'(e25.15)') hin
      read (iunest,'(e25.15)') hbin
      do k=0,nlevp
        read (iunest,'(e25.15)') ((tin(i,j,k),i=1,nxa),j=1,nya)
      enddo
      do k=1,nbed
        read (iunest,'(e25.15)') ((tbedin(i,j,k),i=1,nxa),j=1,nya)
      enddo

      read (iunest,'(e25.15)') sedimin
      do k=1,nsed
        read (iunest,'(e25.15)') ((tsedin(i,j,k),i=1,nxa),j=1,nya)
      enddo
      do k=1,nsed
        read (iunest,'(e25.15)') ((wsedin(i,j,k),i=1,nxa),j=1,nya)
      enddo

      read (iunest,'(e25.15)') topbedeqin
      read (iunest,'(e25.15)') equiloadin
#if defined (OLDNESTFILE)
      read (iunest,'(e25.15)') bmassin                   ! no more bmass
#endif
      read (iunest,'(e25.15)') sedpresin
      read (iunest,'(e25.15)') arcocnin
      read (iunest,'(e25.15)') uain
      read (iunest,'(e25.15)') vain

      read (iunest,'(f15.3)') timecheck

      if (timenest.ne.timecheck) then
        write (ioterm,'(/a/2(a,f12.2))')
     *    '*** Error: timenest, timecheck differ:',
     *    '   timenest=',timenest, '   timecheck=',timecheck
        stop
      endif

      close (iunest)
      return

 1000 write (ioterm,'(/a,i6)') '*** Error in open: iostat=',ios 
      stop

      end

c-----------------------------------------------------------------------

      subroutine setinterpnest (inest,  jnest,  wnest, 
     *                          inestu, jnestu, wnestu,
     *                          inestv, jnestv, wnestv
#if defined (LONLAT)
     *                         ,alon1_in, alon2_in, alat1_in, alat2_in
#endif
     *                         )

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  inest(nx,ny),           jnest(nx,ny),   
     *  wnest(nx,ny,2),
     *  inestu(0:nx+1,0:ny+1),  jnestu(0:nx+1,0:ny+1), 
     *  wnestu(0:nx+1,0:ny+1,2),
     *  inestv(0:nx+1,0:ny+1),  jnestv(0:nx+1,0:ny+1),
     *  wnestv(0:nx+1,0:ny+1,2)

c        Calculate interpolation indices and weights, nested-grid to 
c        continental-grid (used in reading initial data, in getnest). 

c+++++++++++++++++++
#if defined (LONLAT)
c+++++++++++++++++++

      dlon_in = (alon2_in - alon1_in) / nxa
      dlat_in = (alat2_in - alat1_in) / nya

c     h grids:
      do j=1,ny
        do i=1,nx
          zx = (alond(i,j)-alon1_in) / dlon_in
          inest(i,j) = max (1, min (nxa, nint (zx)))
          wnest(i,j,1) = max (0., min (1., zx-inest(i,j)+0.5))

          zy = (alatd(i,j)-alat1_in) / dlat_in
          jnest(i,j) = max (1, min (nya, nint (zy)))
          wnest(i,j,2) = max (0., min (1., zy-jnest(i,j)+0.5))
        enddo
      enddo

c     u grids:
      do j=0,nyp
        jj = max (1, min (ny, j))
        do i=0,nxp
          ii = max (1, min (nx, i))
          if (i.eq.0) then  
            zx = (alond(1,jj)  -0.5*dlond - alon1_in) / dlon_in
          else 
            zx = (alond(ii,jj) -0.5*dlond - alon1_in) / dlon_in
          endif
          inestu(i,j) = max (0, min (nxa+1, nint (zx - 0.5)))
          wnestu(i,j,1) = max (0., min (1., zx-inestu(i,j)))

          jnestu(i,j) = jnest(ii,jj)
          wnestu(i,j,2) = wnest(ii,jj,2)
        enddo
      enddo

c     v grids:
      do i=0,nxp
        ii = max (1, min (nx, i))
        do j=0,nyp
          jj = max (1, min (ny, j))
          if (j.eq.0) then  
            zy = (alatd(ii,1)  - 0.5*dlatd - alat1_in) / dlat_in
          else 
            zy = (alatd(ii,jj) - 0.5*dlatd - alat1_in) / dlat_in
          endif
          jnestv(i,j) = max (0, min (nya+1, nint (zy - 0.5)))
          wnestv(i,j,2) = max (0., min (1., zy-jnestv(i,j)))

          inestv(i,j) = inest(ii,jj)
          wnestv(i,j,1) = wnest(ii,jj,1)
        enddo
      enddo

c+++++++++++++++++++++
#elif defined (STEREO)
c+++++++++++++++++++++

c     h grids:
      do j=1,ny
        do i=1,nx
          zx = (xh(i,j) + dxa*(0.5*nxa)) / dxa
          inest(i,j) = max (1, min (nxa, nint (zx)))
          wnest(i,j,1) = max (0., min (1., zx-inest(i,j)+0.5))

          zy = (yh(i,j) + dya*(0.5*nya)) / dya
          jnest(i,j) = max (1, min (nya, nint (zy)))
          wnest(i,j,2) = max (0., min (1., zy-jnest(i,j)+0.5))
        enddo
      enddo

c     u grids:
      do j=0,nyp
        jj = max (1, min (ny, j))
        do i=0,nxp
          ii = max (1, min (nx, i))
          if (i.eq.0) then  
            zx = (xh(1,jj) - 0.5*dx0 + dxa*(0.5*nxa))  / dxa
          else 
            zx = (xh(min(i,nx),jj) + 0.5*dx0 + dxa*(0.5*nxa)) / dxa
          endif
          inestu(i,j) = max (0, min (nxa+1, nint (zx - 0.5)))
          wnestu(i,j,1) = max (0., min (1., zx-inestu(i,j)))

          jnestu(i,j) = jnest(ii,jj)
          wnestu(i,j,2) = wnest(ii,jj,2)
        enddo
      enddo

c     v grids:
      do i=0,nxp
        ii = max (1, min (nx, i))
        do j=0,nyp
          jj = max (1, min (ny, j))
          if (j.eq.0) then  
            zy = (yh(ii,1) - 0.5*dy0 + dya*(0.5*nya)) / dya
          else 
            zy = (yh(ii,min(j,ny)) + 0.5*dy0 + dya*(0.5*nya)) / dya
          endif
          jnestv(i,j) = max (0, min (nya+1, nint (zy - 0.5)))
          wnestv(i,j,2) = max (0., min (1., zy-jnestv(i,j)))

          inestv(i,j) = inest(ii,jj)
          wnestv(i,j,1) = wnest(ii,jj,1)
        enddo
      enddo

c+++++
#endif
c+++++

      return
      end

c@@@@@
#endif
c@@@@@

c-----------------------------------------------------------------------

c!!!!!!!!!!!!!!!!!!!!!!!
#if defined (CRHMELREAD)
c!!!!!!!!!!!!!!!!!!!!!!!

      subroutine readcrhmel (crhmel, enhanceinv, maskinit, sedpres)

c     Read crhmel and enhanceinv from file written by a previous 
c     CRHINV or ENHINV run, and renamed to 'crhmelfilein'. 
c     Horiz resolution can be different, using bilinear interp here, 
c     assuming both grids are pole-centered stereo (like in 
c     readres/readin/interph).
 
#include <comicephys.h>
#include <comicegrid.h>

      dimension crhmel(nx,ny), enhanceinv(nx,ny), maskinit(nx,ny),
     *          sedpres(nx,ny)

      parameter (nxmax=700, nymax=700)
      dimension crhmelin(nxmax,nymax), enhanceinvin(nxmax,nymax),
     *          alonmelin(nxmax),      alatmelin(nymax)

#if defined (LONLAT)
c     local 2-D ice-grid arrays, indexes+weights to crhmelfinein grid
      dimension
     *  indl(nx,ny),     weil(nx,ny),
     *  indb(nx,ny),     weib(nx,ny)
#endif

      dimension crhmeltmp(nx,ny)

      character cfilin*12, cgridin*80

      dimension icrhmel(nx,ny), jcrhmel(nx,ny), wcrhmel(nx,ny,2)
c     save icrhmel, jcrhmel, wcrhmel   ! not necessary, only called once

#if defined (CRHSHELFASEIN)
c     to define a geographical region for interior Amundsen Sea 
c     hard bed region, Graham etal,QSR,2009,Fig.1):
      parameter (ncreg=12)
      dimension
     *  aloncreg(ncreg), alatcreg(ncreg), xcreg(ncreg), ycreg(ncreg),
     *  allcreg(nx,ny)
      data alatcreg
     * / -85.0, -74.0, -73.6, -73.6, -73.4, -73.9,
     *   -73.3, -72.9, -72.7, -71.7, -72.5, -85.0 /
      data aloncreg
     * / -150., -150., -124., -116., -115., -111.,   
     *   -106.,-105.5, -104., -100.,  -70.,  -70. /
      save allcreg
      logical firstcreg
      save firstcreg
      data firstcreg /.true./
#endif

c        Read crhmelin, enhanceinvin from input file 

      cfilin = 'crhmelfilein'
      write (ioterm,'(a,a)') 
     *  'reading crhmel file ',cfilin(1:lenchr(cfilin))

      open (iucrhmel, file=cfilin, status='old', form='formatted')

      read (iucrhmel,'(a)') cgridin
      if (cgridin.ne.'LONLAT' .and. cgridin.ne.'STEREO') then
        cgridin = 'UNKNOWN'
        backspace iucrhmel
      endif
      write (ioterm,'(2a)')
     *  '  input grid type= ',cgridin(1:lenchr(cgridin))

#if defined (LONLAT)
      if (cgridin.eq.'STEREO') then
#elif defined (STEREO)
      if (cgridin.eq.'LONLAT') then
#endif
        write (ioterm,'(/a)')
     *    '*** Error (readcrhmel): grid type different from model'
        stop
      endif

      read (iucrhmel,'(2i6,4e25.15)') 
     *  nxina, nyina, dxina, dyina, xoffina, yoffina

      if (nxina.gt.nxmax .or. nyina.gt.nymax) then
        write (ioterm,'(/a,2(/a,2i6))')
     *    '*** Error (readcrhmel): input dimensions too large',
     *    '    nxina, nyina =', nxina, nyina,
     *    '    nxmax, nymax =', nxmax, nymax
        stop
      endif

#if defined (LONLAT)
      write (ioterm,'(a,2i5)')
     *  '  nx,ny (input)=', nxina, nyina
      write (ioterm,'(a,2i5)')
     *  '  nx,ny (model)=', nx, ny
c...................................
      if (cgridin.eq.'UNKNOWN') then
c...................................
#  if defined (SAMEDOMAIN)
c       if (nx.ne.nxina.or.ny.ne.nyina) then
c       assume same domain size as model,(can be different grid size):
        do i=1,nxina 
          alonmelin(i)= (alon1 + (alon2-alon1)*(i-.5)/nxina) * pi/180.
        enddo
        do j=1,nyina 
          alatmelin(j)= (alat1 + (alat2-alat1)*(j-.5)/nyina) * pi/180.
        enddo
#  else
        write (ioterm,'(/2a)')
     *    '*** Error (readcrhmel):',
     *    ' old crhmelfile file, different lon-lat grid'
        stop
#  endif
c.........
      else 
c.........
        read (iucrhmel, '(e25.15)') (alonmelin(i),i=1,nxina) 
        read (iucrhmel, '(e25.15)') (alatmelin(j),j=1,nyina) 
        do i=1,nxina
          alonmelin(i) = alonmelin(i)*pi/180.
        enddo
        do j=1,nyina
          alatmelin(j) = alatmelin(j)*pi/180.
        enddo
c..........
      endif
c..........
#else
      write (ioterm,'(a,2i5,4f10.1)')
     *  '  nx,ny,dx0,dy0,xoffa,yoffa (input)=',
     *     nxina, nyina, dxina, dyina, xoffina, yoffina
      write (ioterm,'(a,2i5,4f10.1)')
     *  '  nx,ny,dx0,dy0,xoffa,yoffa (model)=', 
     *     nx, ny, dx0, dy0, xoffa, yoffa
#endif

      read (iucrhmel,'(e25.15)')
     *  ((crhmelin(i,j), i=1,nxina),j=1,nyina) 
      read (iucrhmel,'(e25.15)')
     *  ((enhanceinvin(i,j),i=1,nxina),j=1,nyina)
      close (iucrhmel)

c          If grid dimensions are same, simply copy fields
c          (nb: don't use scopy - input array size isn't nxina,nyina)!

c||||||||||||||||||||||||||||||||||||||||||||
      if (nxina.eq.nx  .and. ny.eq.nyina .and.
     *    dx0.eq.dxina .and. dy0.eq.dyina .and.
     *    xoffa.eq.xoffina .and. yoffa.eq.yoffina) then
c||||||||||||||||||||||||||||||||||||||||||||

        do j=1,ny
          do i=1,nx
            crhmel(i,j) = crhmelin(i,j)
            enhanceinv(i,j) = enhanceinvin(i,j)
          enddo   
        enddo

c|||||||||
      else
c|||||||||

c%%%%%%%%%%%%%%%%%%%
#if defined (LONLAT) 
c%%%%%%%%%%%%%%%%%%%

        call setinterp (alonmelin, alatmelin, nxina, nyina, 
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, 0.)

c       decide if input grid is global in longitude (so wraparound):
        dlonin = (alonmelin(nxina)-alonmelin(1)) / (nxina-1)
        ifwrap = 0
        if (dlonin*nxina. ge. 0.999*2.*pi) ifwrap = 1

        do j=1,ny
          do i=1,nx
c           indices and weights into crhmelfilein's grid (setinterp):
            il = indl(i,j)
            if (ifwrap.eq.0) then
              ir = min (il+1,nxina)    ! walls
            else
              ir = mod (il,nxina) + 1  ! wraparound
            endif
            jb = indb(i,j)
            jt = min (jb+1, nyina)
            wlb =     weil(i,j) *    weib(i,j)
            wrb = (1.-weil(i,j))*    weib(i,j)
            wlt =     weil(i,j) *(1.-weib(i,j))
            wrt = (1.-weil(i,j))*(1.-weib(i,j))

            crhmel(i,j) = 10. ** (   wlb*log10(crhmelin(il,jb))
     *                             + wrb*log10(crhmelin(ir,jb))
     *                             + wlt*log10(crhmelin(il,jt))
     *                             + wrt*log10(crhmelin(ir,jt)) )
            enhanceinv(i,j) =        wlb*enhanceinvin(il,jb)
     *                             + wrb*enhanceinvin(ir,jb)
     *                             + wlt*enhanceinvin(il,jt)
     *                             + wrt*enhanceinvin(ir,jt)
          enddo
        enddo

c%%%%%%%%%%%%%%%%%%%%%
#elif defined (STEREO) 
c%%%%%%%%%%%%%%%%%%%%%

c        Set bilin interp indices,weights,from h-grid into crhmelin grid

        do j=1,ny
          do i=1,nx
            zx = (xh(i,j) - (xoffina-dxina*(0.5*nxina))) / dxina
            icrhmel(i,j) = max (1, min (nxina, nint (zx)))
            wcrhmel(i,j,1) = max (0., min (1., zx-icrhmel(i,j)+0.5))
  
            zy = (yh(i,j) - (yoffina-dyina*(0.5*nyina))) / dyina
            jcrhmel(i,j) = max (1, min (nyina, nint (zy)))
            wcrhmel(i,j,2) = max (0., min (1., zy-jcrhmel(i,j)+0.5))
          enddo
        enddo

c          Do bilinear interp from input-grid to h-grid, using 
c          logarithmic averaging for crhmel, linear for enhanceinv.
 
        do j=1,ny
          do i=1,nx
            ja = jcrhmel(i,j)
            jb = min (nyina, jcrhmel(i,j) + 1)
            ia = icrhmel(i,j)
            ib = min (nxina, icrhmel(i,j) + 1)
  
            zw1 = (1.-wcrhmel(i,j,1)) * (1.-wcrhmel(i,j,2))
            zw2 = (   wcrhmel(i,j,1)) * (1.-wcrhmel(i,j,2))
            zw3 = (1.-wcrhmel(i,j,1)) * (   wcrhmel(i,j,2))
            zw4 = (   wcrhmel(i,j,1)) * (   wcrhmel(i,j,2))

            crhmel(i,j) = 10. ** (   zw1*log10(crhmelin(ia,ja))
     *                             + zw2*log10(crhmelin(ib,ja))
     *                             + zw3*log10(crhmelin(ia,jb))
     *                             + zw4*log10(crhmelin(ib,jb)) )

            enhanceinv(i,j) =        zw1*enhanceinvin(ia,ja)
     *                             + zw2*enhanceinvin(ib,ja)
     *                             + zw3*enhanceinvin(ia,jb)
     *                             + zw4*enhanceinvin(ib,jb)

          enddo
        enddo

c%%%%%
#endif
c%%%%%

c||||||||||
      endif
c||||||||||

c        Propagate "true-bed" values with crhmel >= crhhard
c        into "frozen" areas  with crhmel < crhhard 
c        (crhhard set in comicephys.h, unfrozen hard bedrock)

      call scopy (nx*ny, crhmel, 1, crhmeltmp, 1)
      call fill_crh (crhmel, crhmeltmp, maskinit, crhhard, nx,ny)

#if defined (RWRAP)
c    uniform multiplier for all (modern-grounded) basal coeffs:
      do j=1,ny
        do i=1,nx
          if (maskinit(i,j).eq.0) then
            crhmel(i,j) = crhmel(i,j) * rwrap_crhfac
          endif
        enddo
      enddo
#endif

c        Reset sub-ocean (non-grounded) values to CRHSHELF/CRHSHELFASEIN

c----------------------------------------
#if defined (CRHSHELF) || defined (RWRAP)
c----------------------------------------
#  if defined (CRHSHELFASEIN)
      if (firstcreg) then
        do m=1,ncreg
          call ps_xy (xcreg(m), ycreg(m),
     *                alatcreg(m)*pi/180., aloncreg(m)*pi/180.,
     *                0., 0.)
        enddo

c       set allcreg array (1=inside polygon, 0=outside):
        call zero (allcreg, nx*ny)
        call polygon_wais (allcreg, xh,yh,nx,ny, xcreg,ycreg, ncreg,1.)

c       diagnostic dump:
        do j=ny,1,-1
          write (133,'(1000i1)') (nint(allcreg(i,j)),i=1,nx)
        enddo
        call printmap (0., allcreg, 'allcreg', 1., 133, 0)
        call flush (133)

        firstcreg = .false.
      endif
#  endif

      do j=1,ny
        do i=1,nx
          if (maskinit(i,j).eq.1) then 
#  if defined (RWRAP)
            crhmel(i,j) = rwrap_crhshelf
#  elif defined (CRHSHELF)
            crhmel(i,j) = CRHSHELF
#    if defined (CRHSHELFASEIN)
            if (nint(allcreg(i,j)).eq.1) crhmel(i,j) = 
     *        max ( crhhard, CRHSHELFASEIN * CRHSHELF )
#    endif
#  endif
          endif
        enddo
      enddo
c-----
#endif
c-----

#if defined (CRHGREENLAND) 
      do j=1,ny
        do i=1,nx
          if (maskinit(i,j).ne.1) then 
            crhmel(i,j) = max (crhmel(i,j), CRHGREENLAND)
          endif
        enddo
      enddo
#endif
#if defined (CRHGREENLANDNORTH) 
      do j=1,ny
        do i=1,nx
          if (maskinit(i,j).ne.1 .and. alatd(i,j).gt.70.) then 
            crhmel(i,j) = max (crhmel(i,j), CRHGREENLANDNORTH)
          endif
        enddo
      enddo
#endif

c#if defined (FORCE34MYR) && defined (SHALDRIL) && ! defined (NOSED)
c     for long-term (34 Myr) with explicit sediment, crhmel represents
c     bedrock below pre 40-Ma regolith (initial land) and below
c     pre-40 Ma ocean sediment (initial ocean). For land,
c     impose >= 0.9e-8 to avoid crhinv values representing modern
c     streams (seds). For ocean, use 0.9e-8. (This value is just below 
c     SSA cuttof in icedyn). See icehydrol/basecoef, and icesed/sedflow.
c     ***
c     *** Don't do - keep modern-sedim info, even for 34 Myr sedim runs 
c     ***
c     do j=1,ny
c       do i=1,nx
c         if (sedpres(i,j).eq.0.) then 
c           crhmel(i,j) = min (crhmel(i,j), 0.9e-8)      ! initial land
c         else 
c           crhmel(i,j) = 0.9e-8                         ! initial ocean
c         endif
c       enddo
c     enddo
c#endif

      return
      end

c-----------------------------------------------------------------------

      subroutine fill_crh (crhmel, crhmelfile, maskinit, crhhard,
     *                     nx, ny)

c       For grounded points, propagate "true-rock" values of basal 
c       sliding coeff (from crhmelfile) >= crhhard into areas with
c       < crhhard. Do by neighbor-averaging, one step at a time
c       (similarly to div calc in calving). Return result in crhmel.

c       Called from readcrhmel. Usually used after CRHINV or ENHINV
c       run has generated crhmelfilein with NOBASET (no effect of basal
c       temps/h2o, crhnos lower limit on crhmel). This subr will
c       then fill the "deduced-frozen" areas, crhmel=crhnos,
c       with hard-bedrock values crhmel=crhhard.
c
c       nb: For grounded vs floating, uses "initial" maskinit (in 
c       common, copied from maskwater after initphys), not maskwater 
c       from restart file. Presumably maskinit is ~modern observed, 
c       same as used in the CRHINV/ENHINV run (with NOMOVECOAST) 
c       that generated the crhmelfilein file.

      dimension crhmel(nx,ny), crhmelfile(nx,ny), maskinit(nx,ny)

c     local (dynamic):
      dimension crhtmp(nx,ny)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook4 = 4)   ! 4 for E-W/N-S
      parameter (nlook8 = 8)   ! 8 for diagonals too

      do j=1,ny
        do i=1,nx
          if (maskinit(i,j).eq.0) then 
            if (crhmelfile(i,j).ge.crhhard) then
              crhmel(i,j) = crhmelfile(i,j)
            else
              crhmel(i,j) = 0.
            endif
          endif
        enddo
      enddo

c,,,,,,,,,,,,,,,,,,,,,,,
      do iter = 1, nx*ny
c,,,,,,,,,,,,,,,,,,,,,,,
        do j=1,ny
          do i=1,nx
            crhtmp(i,j) = crhmel(i,j)
          enddo
        enddo
        nnew = 0

        do j=1,ny
          do i=1,nx
            if (maskinit(i,j).eq.0 .and. crhtmp(i,j).eq.0.) then
              zv = 0.
              nv = 0
              do look=1,nlook4
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (maskinit(ii,jj).eq.0 .and. crhtmp(ii,jj).ne.0.)then
                  zv = zv + crhtmp(ii,jj)
                  nv = nv + 1
                endif
              enddo
              if (nv.gt.0) then
                crhmel(i,j) = zv / nv
                nnew = nnew + 1
              endif
            endif
          enddo
        enddo
        if (nnew.eq.0) go to 500
c,,,,,,,,,,
      enddo
c,,,,,,,,,,
      write (6,*) 'should not get here (iter max: fill_crh)'
      stop
  500 continue

c         For safety in case isolated points not reached

      do j=1,ny
        do i=1,nx
          if (maskinit(i,j).eq.0 .and. crhmel(i,j).eq.0.) then
            crhmel(i,j) = crhhard
          endif
        enddo
      enddo

      return
      end

c!!!!!
#endif
c!!!!!

c-----------------------------------------------------------------------

#if defined (CRHINV) || defined (ENHINV)

      subroutine writecrhmel (crhmel, enhanceinv, timeice)

c     Writes crhmel file, which can be used to initialize crhmel 
c     and enhanceinv at cold start of a subsequent run with 
c     CRHMELREAD requested (by subr readcrhmel).

#include <comicephys.h>
#include <comicegrid.h>

      dimension crhmel(nx,ny), enhanceinv(nx,ny)

      character cfile*10
      parameter (cfile = 'crhmelfile')
      character cfile2*80, cyr*80, cmd*200

      write (ioterm,'(/a,a)') 'Writing crhmel file ', cfile

      open (iucrhmel, file=cfile,status='unknown',form='formatted')
      close (iucrhmel, status='delete')
      open (iucrhmel, file=cfile,status='new',form='formatted')
#if defined (LONLAT)
      write (iucrhmel,'(a)') 'LONLAT'
#elif defined (STEREO)
      write (iucrhmel,'(a)') 'STEREO'
#endif
      write (iucrhmel,'(2i6,4e25.15)') nx, ny, dx0, dy0, xoffa, yoffa
#if defined (LONLAT)
      write (iucrhmel,'(e25.15)') (alond(i,(ny+1)/2),i=1,nx)
      write (iucrhmel,'(e25.15)') (alatd((nx+1)/2,j),j=1,ny)
#endif
      write (iucrhmel,'(e25.15)') crhmel
      write (iucrhmel,'(e25.15)') enhanceinv
      close (iucrhmel)

c        Save 2nd copy of crhmel file to unique filename

c     write (cyr, '(i10)') nint(abs(timeice))
c     do i=1,lenchr(cyr)
c       if (cyr(1:1).eq.' ') cyr = cyr(2:)
c     enddo
c     if (timeice.lt.0.) cyr = cyr(1:lenchr(cyr)) // 'm'
c     if (timeice.gt.0.) cyr = cyr(1:lenchr(cyr)) // 'p'
c     cfile2 = cfile //  '_' // cyr

c     write (ioterm,'(a,a)')
c    *  'Copying to crhmel file ', cfile2(1:lenchr(cfile2))
c     cmd = 'cp -p ' // cfile // ' ' // cfile2
c     ier = ishell (cmd)

      return
      end

#endif

c-----------------------------------------------------------------------

#if defined (SEDIMENT) && defined (FORCE34MYR) && defined (QUARRYREAD)

      subroutine readquarry (hbinitprev, hbrebprev, hbinit, quarrysec,
     *                       timeice)

c     Reads (1) quarrying coefficient each sector (quarrysec) and 
c     (2) bedrock elevations for this iteration (hbinitprev, hbrebprev,
c     hbinit), from file quarryfilein, copied from quarryfile 
c     written at end of prev run in subr writequarry if QUARRYWRITE.
c     Similar to subr readcrhmel.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  hbinitprev(nx,ny), hbrebprev(nx,ny), 
     *  hbinit(nx,ny), quarrysec(nsedsec)

      parameter (nxmax=700, nymax=700)
      dimension 
     *  hbinitprevin(nxmax,nymax), hbrebprevin(nxmax,nymax), 
     *  hbinitin(nxmax,nymax), 
     *  alonquarin(nxmax),      alatquarin(nymax)

#if defined (LONLAT)
c     local 2-D ice-grid arrays, indexes+weights to quarryfinein grid
      dimension
     *  indl(nx,ny),     weil(nx,ny),
     *  indb(nx,ny),     weib(nx,ny)
#endif

      character cfilin*12, cgridin*80, cdum*8, chead*24

      dimension iquarry(nx,ny), jquarry(nx,ny), wquarry(nx,ny,2)
c     save iquarry, jquarry, wquarry   ! not necessary, only called once

c        Read quarrysec, hb... from input file 

      cfilin = 'quarryfilein'
      write (ioterm,'(a,a)') 
     *  'reading quarry file ',cfilin(1:lenchr(cfilin))

      open (iuquarry, file=cfilin, status='old', form='formatted')

      read (iuquarry,*) cdum, nsedsecin
      if (nsedsecin.ne.nsedsec) then
        write (ioterm,'(a,2i6)') 
     *    '*** Error (readquarry): nsedsec[in]=', nsedsec, nsedsecin 
        stop
      endif
      do m=1,nsedsec
        read (iuquarry,*) quarrysec(m)
      enddo
      write (ioterm,'(a/(2x,e20.10))') 
     *  '  quarrysec:', (quarrysec(m),m=1,nsedsec)

      read (iuquarry,'(a)') cgridin
      if (cgridin.ne.'LONLAT' .and. cgridin.ne.'STEREO') then
        cgridin = 'UNKNOWN'
        backspace iuquarry
      endif
      write (ioterm,'(2a)')
     *  '  input grid type= ',cgridin(1:lenchr(cgridin))

#if defined (LONLAT)
      if (cgridin.eq.'STEREO') then
#elif defined (STEREO)
      if (cgridin.eq.'LONLAT') then
#endif
        write (ioterm,'(/a)')
     *    '*** Error (readquarry): grid type different from model'
        stop
      endif

      read (iuquarry,'(2i6,4e25.15)') 
     *  nxina, nyina, dxina, dyina, xoffina, yoffina

      if (nxina.gt.nxmax .or. nyina.gt.nymax) then
        write (ioterm,'(/a,2(/a,2i6))')
     *    '*** Error (readquarry): input dimensions too large',
     *    '    nxina, nyina =', nxina, nyina,
     *    '    nxmax, nymax =', nxmax, nymax
        stop
      endif

#if defined (LONLAT)
      write (ioterm,'(a,2i5)')
     *  '  nx,ny (input)=', nxina, nyina
      write (ioterm,'(a,2i5)')
     *  '  nx,ny (model)=', nx, ny
c...................................
      if (cgridin.eq.'UNKNOWN') then
c...................................
#  if defined (SAMEDOMAIN)
c       if (nx.ne.nxina.or.ny.ne.nyina) then
c       assume same domain size as model,(can be different grid size):
        do i=1,nxina 
          alonquarin(i)= (alon1 + (alon2-alon1)*(i-.5)/nxina) * pi/180.
        enddo
        do j=1,nyina 
          alatquarin(j)= (alat1 + (alat2-alat1)*(j-.5)/nyina) * pi/180.
        enddo
#  else
        write (ioterm,'(/2a)')
     *    '*** Error (readquarry):',
     *    ' old quarryfile file, different lon-lat grid'
        stop
#  endif
c.........
      else 
c.........
        read (iuquarry, '(e25.15)') (alonquarin(i),i=1,nxina) 
        read (iuquarry, '(e25.15)') (alatquarin(j),j=1,nyina) 
        do i=1,nxina
          alonquarin(i) = alonquarin(i)*pi/180.
        enddo
        do j=1,nyina
          alatquarin(j) = alatquarin(j)*pi/180.
        enddo
c..........
      endif
c..........
#else
      write (ioterm,'(a,2i5,4f10.1)')
     *  '  nx,ny,dx0,dy0,xoffa,yoffa (input)=',
     *     nxina, nyina, dxina, dyina, xoffina, yoffina
      write (ioterm,'(a,2i5,4f10.1)')
     *  '  nx,ny,dx0,dy0,xoffa,yoffa (model)=', 
     *     nx, ny, dx0, dy0, xoffa, yoffa
#endif

      read (iuquarry,'(a)') chead
      if (chead.ne.'HBINITPREV:') then
        write (ioterm,'(2a)')
     *    '*** Error (readquarry), header HBINITPREV: ', chead
        stop
      endif
      read (iuquarry,'(e25.15)')
     *  ((hbinitprevin(i,j), i=1,nxina),j=1,nyina) 

      read (iuquarry,'(a)') chead
      if (chead.ne.'HBREBPREV:') then
        write (ioterm,'(2a)') 
     *    '*** Error (readquarry), header HBREBPREV: ', chead
        stop
      endif
      read (iuquarry,'(e25.15)')
     *  ((hbrebprevin(i,j), i=1,nxina),j=1,nyina) 

      read (iuquarry,'(a)') chead
      if (chead.ne.'HBINIT:') then
        write (ioterm,'(2a)')
     *    '*** Error (readquarry), header HBINIT: ', chead
        stop
      endif
      read (iuquarry,'(e25.15)')
     *  ((hbinitin(i,j), i=1,nxina),j=1,nyina) 

      close (iuquarry)

c          If grid dimensions are same, simply copy fields
c          (nb: don't use scopy - input array size isn't nxina,nyina)!

c||||||||||||||||||||||||||||||||||||||||||||
      if (nxina.eq.nx  .and. ny.eq.nyina .and.
     *    dx0.eq.dxina .and. dy0.eq.dyina .and.
     *    xoffa.eq.xoffina .and. yoffa.eq.yoffina) then
c||||||||||||||||||||||||||||||||||||||||||||

        do j=1,ny
          do i=1,nx
            hbinitprev(i,j) = hbinitprevin(i,j)
            hbrebprev(i,j)  = hbrebprevin(i,j)
            hbinit(i,j)     = hbinitin(i,j)
          enddo   
        enddo

c|||||||||
      else
c|||||||||

c%%%%%%%%%%%%%%%%%%%
#if defined (LONLAT) 
c%%%%%%%%%%%%%%%%%%%

        call setinterp (alonquarin, alatquarin, nxina, nyina, 
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, 0.)

c       decide if input grid is global in longitude (so wraparound):
        dlonin = (alonquarin(nxina)-alonquarin(1)) / (nxina-1)
        ifwrap = 0
        if (dlonin*nxina. ge. 0.999*2.*pi) ifwrap = 1

        do j=1,ny
          do i=1,nx
c           indices and weights into quarryfilein's grid (setinterp):
            il = indl(i,j)
            if (ifwrap.eq.0) then
              ir = min (il+1,nxina)    ! walls
            else
              ir = mod (il,nxina) + 1  ! wraparound
            endif
            jb = indb(i,j)
            jt = min (jb+1, nyina)
            wlb =     weil(i,j) *    weib(i,j)
            wrb = (1.-weil(i,j))*    weib(i,j)
            wlt =     weil(i,j) *(1.-weib(i,j))
            wrt = (1.-weil(i,j))*(1.-weib(i,j))

            hbinitprev(i,j) = (   wlb*(hbinitprevin(il,jb))
     *                          + wrb*(hbinitprevin(ir,jb))
     *                          + wlt*(hbinitprevin(il,jt))
     *                          + wrt*(hbinitprevin(ir,jt)) )

            hbrebprev(i,j)  = (   wlb*(hbrebprevin(il,jb))
     *                          + wrb*(hbrebprevin(ir,jb))
     *                          + wlt*(hbrebprevin(il,jt))
     *                          + wrt*(hbrebprevin(ir,jt)) )

            hbinit(i,j)     = (   wlb*(hbinitin(il,jb))
     *                          + wrb*(hbinitin(ir,jb))
     *                          + wlt*(hbinitin(il,jt))
     *                          + wrt*(hbinitin(ir,jt)) )
          enddo
        enddo

c%%%%%%%%%%%%%%%%%%%%%
#elif defined (STEREO) 
c%%%%%%%%%%%%%%%%%%%%%

c        Set bilin interp indices,weights,from h-grid into input file
c        grid

        do j=1,ny
          do i=1,nx
            zx = (xh(i,j) - (xoffina-dxina*(0.5*nxina))) / dxina
            iquarry(i,j) = max (1, min (nxina, nint (zx)))
            wquarry(i,j,1) = max (0., min (1., zx-iquarry(i,j)+0.5))
  
            zy = (yh(i,j) - (yoffina-dyina*(0.5*nyina))) / dyina
            jquarry(i,j) = max (1, min (nyina, nint (zy)))
            wquarry(i,j,2) = max (0., min (1., zy-jquarry(i,j)+0.5))
          enddo
        enddo

c          Do bilinear interp from input-grid to h-grid
 
        do j=1,ny
          do i=1,nx
            ja = jquarry(i,j)
            jb = min (nyina, jquarry(i,j) + 1)
            ia = iquarry(i,j)
            ib = min (nxina, iquarry(i,j) + 1)
  
            zw1 = (1.-wquarry(i,j,1)) * (1.-wquarry(i,j,2))
            zw2 = (   wquarry(i,j,1)) * (1.-wquarry(i,j,2))
            zw3 = (1.-wquarry(i,j,1)) * (   wquarry(i,j,2))
            zw4 = (   wquarry(i,j,1)) * (   wquarry(i,j,2))

            hbinitprev(i,j) = (   zw1*(hbinitprevin(ia,ja))
     *                          + zw2*(hbinitprevin(ib,ja))
     *                          + zw3*(hbinitprevin(ia,jb))
     *                          + zw4*(hbinitprevin(ib,jb)) )

            hbrebprev(i,j) =  (   zw1*(hbrebprevin(ia,ja))
     *                          + zw2*(hbrebprevin(ib,ja))
     *                          + zw3*(hbrebprevin(ia,jb))
     *                          + zw4*(hbrebprevin(ib,jb)) )

            hbinit(i,j) =     (   zw1*(hbinitin(ia,ja))
     *                          + zw2*(hbinitin(ib,ja))
     *                          + zw3*(hbinitin(ia,jb))
     *                          + zw4*(hbinitin(ib,jb)) )

          enddo
        enddo

c%%%%%
#endif
c%%%%%

c||||||||||
      endif
c||||||||||

      return
      end

#endif

c-----------------------------------------------------------------------

#if defined (SEDIMENT) && defined (FORCE34MYR) && defined (QUARRYWRITE)

      subroutine writequarry (hbinitprev, hbrebprev, hbinit, quarrysec,
     *                        timeice)

c     Writes (1) quarrying coefficient each sector (quarrysec) and 
c     (2) bedrock elevations, for next iteration (hbinitprev, hbrebprev,
c     hbinit), to file quarryfile (and quarryfile_<year>), to be 
c     copied to quarryfilein, read by next iteration in subr 
c     readquarry if QUARRYREAD.
c     Similar to subr writecrhmel.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  hbinitprev(nx,ny), hbrebprev(nx,ny), hbinit(nx,ny),
     *  quarrysec(nsedsec)

      character cfile*10
      parameter (cfile = 'quarryfile')
      character cfile2*80, cyr*80, cmd*200

      write (ioterm,'(/a,a)') 'Writing quarry file ', cfile

      open (iuquarry, file=cfile,status='unknown',form='formatted')
      close (iuquarry, status='delete')
      open (iuquarry, file=cfile,status='new',form='formatted')

      write (iuquarry,'(a,i4/(2x,e20.10))') 'quarrysec:',nsedsec,
     *                                      (quarrysec(m),m=1,nsedsec)

#if defined (LONLAT)
      write (iuquarry,'(a)') 'LONLAT'
#elif defined (STEREO)
      write (iuquarry,'(a)') 'STEREO'
#endif
      write (iuquarry,'(2i6,4e25.15)') nx, ny, dx0, dy0, xoffa, yoffa
#if defined (LONLAT)
      write (iuquarry,'(e25.15)') (alond(i,(ny+1)/2),i=1,nx)
      write (iuquarry,'(e25.15)') (alatd((nx+1)/2,j),j=1,ny)
#endif

      write (iuquarry,'(a)') 'HBINITPREV:'
      write (iuquarry,'(e25.15)') hbinitprev
      write (iuquarry,'(a)') 'HBREBPREV:'
      write (iuquarry,'(e25.15)') hbrebprev
      write (iuquarry,'(a)') 'HBINIT:'
      write (iuquarry,'(e25.15)') hbinit
      close (iuquarry)

c        Save 2nd copy of file to unique filename

c     write (cyr, '(i10)') nint(abs(timeice))
c     do i=1,lenchr(cyr)
c       if (cyr(1:1).eq.' ') cyr = cyr(2:)
c     enddo
c     if (timeice.lt.0.) cyr = cyr(1:lenchr(cyr)) // 'm'
c     if (timeice.gt.0.) cyr = cyr(1:lenchr(cyr)) // 'p'
c     cfile2 = cfile //  '_' // cyr

c     write (ioterm,'(a,a)')
c    *  'Copying to quarry file ', cfile2(1:lenchr(cfile2))
c     cmd = 'cp -p ' // cfile // ' ' // cfile2
c     ier = ishell (cmd)

      return
      end

#endif

c-----------------------------------------------------------------------

#if defined (NOISEHB)

      subroutine gaussnoise (arr, sigout, sigin, nsmooth, work, nx, ny)

c     Use Box-Muller transformation to generate gaussian (normal) noise 
c     (mean=0, std.dev=sigin), and apply it at each point independently 
c     on a 2-D array arr(nx,ny).
c
c     If nsmooth>=1, apply uniform moving average filter to arr 
c     (with rectangular box size nsmooth x nsmooth). 
c     Calculate new std.dev (sigout).

c     To check, also calculates pdf(npdf) of generated noise,
c     and writes to stdout.

      dimension arr(nx,ny), work(nx,ny)

      parameter (pi = 3.14159265358979)
      parameter (npdf=30, xmin=-3., xmax=3.)
      parameter (dx=(xmax-xmin)/npdf)
      dimension pdf(npdf)

      do m = 1,npdf
        pdf(m) = 0.
      enddo

      call random_seed()

      do j=1,ny
        do i=1,nx
   10     call random_number(x1)  !  uniform on (0,1)
          x1 = 2.*x1 -1.
          call random_number(x2)  !  uniform on (0,1)
          x2 = 2.*x2 -1.
          if (x1.lt.1.e-20 .or. x1.gt.1.-1.e-20 .or.
     *        x2.lt.1.e-20 .or. x2.gt.1.-1.e-20) goto 10

          xa = sqrt(-2.*log(x1)) * cos(2.*pi*x2)
c         xb = sqrt(-2.*log(x1)) * sin(2.*pi*x2)   ! don't use xb

          m = nint ((xa - xmin)/dx + 0.5)
          m = max (1, min (npdf, m))
          pdf(m) = pdf(m) + 1.

          arr(i,j) = xa * sigin
        enddo
      enddo

      pnorm = 0.5*(pdf(npdf/2) + pdf((npdf+1)/2))
      do m = 2,npdf-1
        write (6,'(2f8.3)') xmin + dx*(m-0.5), pdf(m)/pnorm
      enddo

c---------------------------
      if (nsmooth.gt.0) then 
c---------------------------
        nhalf = max (1, nsmooth/2)

        do j=1,ny  
          do i=1,nx
            work(i,j) = arr(i,j)
          enddo
        enddo

        do j=1,ny  
          ja = max (1,  j-nhalf)
          jb = min (ny, j+nhalf)
          do i=1,nx
            ia = max (1,  i-nhalf)
            ib = min (nx, i+nhalf)
            zav = 0.
            nav = 0
            do jj=ja,jb
              do ii=ia,ib
                zav = zav + work(ii,jj)
                nav = nav + 1
              enddo
            enddo
            arr(i,j) = zav/nav
          enddo
        enddo
c----------
      endif
c----------

      sigout = 0.
      do j=1,ny  
        do i=1,nx
          sigout = sigout + arr(i,j)**2
        enddo
      enddo
      sigout = sqrt (sigout/(nx*ny))

      return
      end

#endif

c-----------------------------------------------------------------------

c###############################################################
#if defined (EISANTA) && defined (TRANSECTA) && defined (ZOOMGL)
c###############################################################

      subroutine rezoomgl (iczoom, maskwater,
     *                     h, hs, hb, tracer,
     *                     sedim, tsed, wsed,
     *                     topbed, tbed, hw, tw,
     *                     baseperc, basefrml, 
     *                     heatb, heats, heatw,
     *                     topbedeq, equiload, hbsd, crhmel,
     *                     u, v, ub, vb, uia, via, uit, vit, ua, va,
     *                     sih, siu, siv, dbu, dbv)
    
c          If grounding line has migrated out of central zoomed bin,
c          reset iczoom and grid vars, and regrid model variables.

#include <comicegrid.h>

c     h grid:
      dimension
     *  maskwater(nx,ny),
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny),
     *  tracer(nx,ny,0:nlevp,ntrace),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     *  topbed(nx,ny),       tbed(nx,ny,nbed),    
     *  hw(nx,ny),           tw(nx,ny),
     *  baseperc(nx,ny),     basefrml(nx,ny),    
     *  heatb(nx,ny),        heats(nx,ny,nsed),  heatw(nx,ny),
     *  topbedeq(nx,ny),     equiload(nx,ny),     
     *  hbsd(nx,ny),         crhmel(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),       via(0:nxp,0:nyp),
     *  uit(0:nxp,0:nyp),       vit(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  sih(nx,ny),
     *  siu(0:nxp,0:nyp),       siv(0:nxp,0:nyp),
     *  dbu(0:nxp,0:nyp),       dbv(0:nxp,0:nyp)

      dimension
     *  xhold(nx,ny),           yhold(nx,ny)

      dimension
     *  xhsav(nx,ny),           yhsav(nx,ny),
     *  topbedeqsav(nx,ny),     equiloadsav(nx,ny),
     *  hbsdsav(nx,ny),         crhmelsav(nx,ny)
      logical first
      data first /.true./
      save first, xhsav, yhsav, topbedeqsav, equiloadsav, 
     *     hbsdsav, crhmelsav

c        If first call this run, save variables that should be invariant

      if (first) then
        call scopy (nx*ny, xh, 1, xhsav, 1)
        call scopy (nx*ny, yh, 1, yhsav, 1)
        call scopy (nx*ny, topbedeq, 1, topbedeqsav, 1)
        call scopy (nx*ny, equiload, 1, equiloadsav, 1)
        call scopy (nx*ny, hbsd,     1, hbsdsav,     1)
        call scopy (nx*ny, crhmel,   1, crhmelsav,   1)
        first = .false.
      endif

      j = (ny+1)/2

c        Find grounding line, searching from rhs (from i=nx to i=1).
c        Set model-grid index igl, and polar stereo map coords (xgl,ygl)

      igl = 0
      do i=nx-1,1,-1
        if (h(i,j).gt.0. .and. maskwater(i,j).eq.0) then
          igl = i
          xgl = xh(i,j)
          ygl = yh(i,j)
          goto 50
        endif
      enddo
   50 continue

c        If no g.l., set iczoom to rhs of domain, return 

      if (igl.eq.0) then  
        iczoom = nx0
        return
      endif

c        Test if g.l. is within current iczoom bin (iga to igb).
c        If so, return with unchanged iczoom.
c        If not, find all-zoom grid index (mgl) closest to g.l., and
c        set new iczoom to the coarse bin containing mgl.

      iga = (iczoom-1)*npzoom + 1
      igb = (iczoom  )*npzoom

      if (igl.ge.iga .and. igl.le.igb) then
        return
      else
        zdmin = 1.e20
        do m=1,nzoom
          zd = (xhzoom(m)-xgl)**2 + (yhzoom(i)-ygl)**2
          if (zd.lt.zdmin) then
            mgl = m 
            zdmin = zd
          endif
        enddo
        iczoomnew = max (1, min (nx0, (mgl-1)/npzoom + 1))
      endif

      write (iuout,'(/a/2i4,3f10.3/)') 
     *        '*** Regridding (zoomgl):',
     *        '    iczoom[old,new], distlen[old,new], distzoom(new) =',
     *             iczoom, iczoomnew, 
     *             distlen(iczoom)/1.e3, distlen(iczoomnew)/1.e3, 
     *             distzoom(mgl)/1.e3
      iczoom = iczoomnew

c        Store old xh,yh in xhold,yhold, then reset xh,yh 
c        and other grid variables 

      call scopy (nx*ny, xh, 1, xhold, 1)
      call scopy (nx*ny, yh, 1, yhold, 1)
      call initgrid (iczoom)

c        Regrid prognostic variables (and some diagnostic)

c                                            nvert  ifreal iflog iffirst
      call regrid  (h,         xhold, yhold, 1,     1,     0,    1)
      call regrid  (hs,        xhold, yhold, 1,     1,     0,    0)
      call regrid  (hb,        xhold, yhold, 1,     1,     0,    0)
      call regrid  (hbsd,      xhold, yhold, 1,     1,     0,    0)
c     t equiv to tracer(...,1):
      call regrid  (tracer,    xhold, yhold, (nlevp+1)*ntrace, 1, 0, 0)
      call regrid  (sedim,     xhold, yhold, 1,     1,     0,    0)
      call regrid  (tsed,      xhold, yhold, nsed,  1,     0,    0)
      call regrid  (wsed,      xhold, yhold, nsed,  1,     0,    0)
      call regrid  (topbed,    xhold, yhold, 1,     1,     0,    0)
      call regrid  (tbed,      xhold, yhold, nbed,  1,     0,    0)
      call regrid  (hw,        xhold, yhold, 1,     1,     0,    0)
      call regrid  (tw,        xhold, yhold, 1,     1,     0,    0)
      call regrid  (sedtrack,  xhold, yhold, ntrack+1,1,   0,    0)
      call regridi (itrtop,    xhold, yhold, 1,     0,     0,    0)

c        For vars that should be invariant, always interpolate from 
c        first-called state (*sav, saved locally, to avoid diffusion 
c        with repeated interps) 

      call scopy (nx*ny, topbedeqsav, 1, topbedeq, 1)
      call scopy (nx*ny, equiloadsav, 1, equiload, 1)
      call scopy (nx*ny, hbsdsav,     1, hbsd,     1)
      call scopy (nx*ny, crhmelsav,   1, crhmel,   1)
      call regrid  (topbedeq,  xhsav, yhsav, 1,     1,     0,    1)
      call regrid  (equiload,  xhsav, yhsav, 1,     1,     0,    0)
      call regrid  (hbsd,      xhsav, yhsav, 1,     1,     0,    0)
      call regrid  (crhmel,    xhsav, yhsav, 1,     1,     0,    0)

c        Zero diagnostic variables for safety

      call zero (baseperc, nx*ny)
      call zero (basefrml, nx*ny)
      call zero (heatb,    nx*ny)
      call zero (heats,    nx*ny*nsed)
      call zero (heatw,    nx*ny)

      call zero (u,   (nxp+1)*(nyp+1)*(nlevp+1))
      call zero (v,   (nxp+1)*(nyp+1)*(nlevp+1))
      call zero (ub,  (nxp+1)*(nyp+1))
      call zero (vb,  (nxp+1)*(nyp+1))
      call zero (uia, (nxp+1)*(nyp+1))
      call zero (via, (nxp+1)*(nyp+1))
      call zero (uit, (nxp+1)*(nyp+1))  
      call zero (vit, (nxp+1)*(nyp+1))
      call zero (ua,  (nxp+1)*(nyp+1))
      call zero (va,  (nxp+1)*(nyp+1))
      call zero (sih, nx*ny)
      call zero (siu, (nxp+1)*(nyp+1))
      call zero (siv, (nxp+1)*(nyp+1))
      call zero (dbu, (nxp+1)*(nyp+1))
      call zero (dbv, (nxp+1)*(nyp+1))

      return
      end

c-----------------------------------------------------------------------

      subroutine regrid (arr, xhold, yhold, nvert, ifreal,iflog,iffirst)
   
      entry     regridi (iarr,xhold, yhold, nvert, ifreal,iflog,iffirst)
   
c     Regrid prognostic variable [i]arr from [xh,yh]old to [xh,yh] grid,
c     for zoom capability tracking grounding line
   
#include <comicegrid.h>

      dimension 
     *  arr(nx,ny,nvert),   iarr(nx,ny,nvert),
     *  xhold(nx,ny,nvert), yhold(nx,ny,nvert)

      dimension arrold(nx,ny), iarrold(nx,ny)

      dimension indh(nx,ny), jndh(nx,ny), weih(nx,ny,2)    ! for interph
      save indh, jndh, weih
   
      do k=1,nvert
        if (ifreal.eq.1) then
          call scopy (nx*ny, arr(1,1,k), 1, arrold, 1)
          call interph  (arr(1,1,k),  xh,    yh,    nx, ny,  
     *                   arrold,      xhold, yhold, nx, ny, 
     *                   indh, jndh, weih,
     *                   ifreal, iflog, iffirst)
        else    
          call scopy_i (nx*ny, iarr(1,1,k), 1, iarro, 1)
          call interphi (iarr(1,1,k), xh,    yh,    nx, ny,  
     *                   iarrold,     xhold, yhold, nx, ny, 
     *                   indh, jndh, weih,
     *                   ifreal, iflog, iffirst)
        endif   
      enddo   
   
      return  
      end    

c#####
#endif
c#####

c-----------------------------------------------------------------------
