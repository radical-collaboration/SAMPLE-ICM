      subroutine bedrock (h, hb, hw, topbed, topbedeq, deflect,
     *                    sedim, equiload, sealev, timeice, dt)

c        Steps bedrock elevation topbed through one timestep, then 
c        adjusts hb. Computes deflection of elastic
c        lithosphere +fluid mantle (deflect), then computes 
c        time-dependent relaxation of topbed to that profile 
c        either locally (kasth=0) or for asthenospheric channel flow.

#include <comicephys.h>
#include <comicegrid.h>

c        h      = ice-sheet thickness (m) (supplied)
c        hb     = bedrock+sedim elevation (m) (modified)
c        deflect  = lithospheric equilibrated defection (m) returned
c        topbed = bedrock elevation, not including sedim (m) (modified)
c        topbedeq = equilibrium bedrock elevation (m) (supplied)
c        sedim = sediment thickness (m) (supplied)
c        equiload = equil. load on bedrock (kg/m2) (supplied)
c        dt = bedrock-model timestep (yr) (supplied)

      dimension 
     *  h(nx,ny),         hb(nx,ny),         hw(nx,ny), 
     *  topbed(nx,ny),    topbedeq(nx,ny),   deflect(nx,ny),
     *  sedim(nx,ny),     equiload(nx,ny)

c     kasth=0 for local asthenospheric relax, 1 for channel flow:
c     taulocal = relaxation e-folding time (for kasth=0):
c     asthflow = diffusive bedrock flow coeff (m2/y) (for kasth=1):
#if defined (TAUASTH)
      parameter (kasth = 0, taulocal = TAUASTH) 
#else
      parameter (kasth = 0, taulocal = 3000.)   ! 5000.
#endif
c     parameter (kasth = 1, asthflow = 1.e8)

      dimension zmass(nx,ny), za(nx,ny)
      
c>>>>>>>>>>>>>>>>>>>
#if defined (IFLITH)
c>>>>>>>>>>>>>>>>>>>

c        For elastic lithospheric flexure (2D: Brotchie and Silvester,
c        1969,JGR,74,22,5240, and 1D: Turcotte and Schubert,1982,p.125)

      parameter (cruststiff_e = 1.e25) !Huybrechts,cf. 3.76e23 B+S (N m)
#  if defined (EISANTA)
#    if defined (LITHSTIF_W) 
      parameter (cruststiff_w = LITHSTIF_W) 
#    else
      parameter (cruststiff_w = 1.e24)
c     parameter (cruststiff_w = 1.e25)
#    endif
#  else
      parameter (cruststiff_w = cruststiff_e)
#  endif

c     for parallelization:
#if defined (EISLINE) 
      parameter (nseg=1)
#else
      parameter (nseg=4)
#endif
      dimension deflectseg(nx,ny,nseg)

      logical first
      save first, 
     *  crustlen_e, crustcoef_e,
     *  crustlen_w, crustcoef_w
      data first /.true./

c      Kelvin Function kei (0th order, from IMSL, for 2D lithosphere): 
c      akeic is stored coarse-resol (dkeic), linearly interpolated to
c      finer-resol akei (dkei) in first call:
       parameter (xkei=15., dkeic = .05,  nkeic=xkei/dkeic + .001,
     *                      dkei  = .005, nkei =xkei/dkei  + .001) 
       dimension akeic(0:nkeic), akei(0:nkei) 
       save akei
       data akeic / 
     * -0.7854,-0.7828,-0.7769,-0.7684,-0.7581,-0.7463,-0.7331,-0.7189,
     * -0.7038,-0.6880,-0.6716,-0.6547,-0.6374,-0.6199,-0.6022,-0.5843,
     * -0.5664,-0.5484,-0.5305,-0.5127,-0.4950,-0.4775,-0.4601,-0.4430,
     * -0.4262,-0.4096,-0.3933,-0.3773,-0.3617,-0.3464,-0.3314,-0.3168,
     * -0.3026,-0.2887,-0.2752,-0.2621,-0.2494,-0.2371,-0.2251,-0.2136,
     * -0.2024,-0.1916,-0.1812,-0.1711,-0.1614,-0.1521,-0.1431,-0.1345,
     * -0.1262,-0.1183,-0.1107,-0.1034,-0.0964,-0.0898,-0.0834,-0.0774,
     * -0.0716,-0.0661,-0.0608,-0.0558,-0.0511,-0.0466,-0.0424,-0.0384,
     * -0.0346,-0.0310,-0.0276,-0.0244,-0.0214,-0.0186,-0.0160,-0.0135,
     * -0.0112,-0.0091,-0.0071,-0.0052,-0.0035,-0.0019,-0.0004, 0.0009,
     *  0.0022, 0.0033, 0.0044, 0.0053, 0.0062, 0.0070, 0.0077, 0.0083,
     *  0.0088, 0.0093, 0.0097, 0.0101, 0.0104, 0.0106, 0.0108, 0.0110,
     *  0.0111, 0.0112, 0.0112, 0.0112, 0.0112, 0.0111, 0.0111, 0.0109,
     *  0.0108, 0.0107, 0.0105, 0.0103, 0.0101, 0.0099, 0.0097, 0.0095,
     *  0.0093, 0.0090, 0.0088, 0.0085, 0.0083, 0.0080, 0.0077, 0.0075,
     *  0.0072, 0.0070, 0.0067, 0.0064, 0.0062, 0.0059, 0.0057, 0.0054,
     *  0.0052, 0.0050, 0.0047, 0.0045, 0.0043, 0.0041, 0.0038, 0.0036,
     *  0.0034, 0.0032, 0.0031, 0.0029, 0.0027, 0.0025, 0.0024, 0.0022,
     *  0.0021, 0.0019, 0.0018, 0.0016, 0.0015, 0.0014, 0.0013, 0.0012,
     *  0.0010, 0.0009, 0.0008, 0.0008, 0.0007, 0.0006, 0.0005, 0.0004,
     *  0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0.0001, 0.0000, 0.0000,
     *  0.0000,-0.0001,-0.0001,-0.0001,-0.0002,-0.0002,-0.0002,-0.0002,
     * -0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,
     * -0.0003,-0.0004,-0.0004,-0.0004,-0.0004,-0.0004,-0.0004,-0.0004,
     * -0.0004,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,
     * -0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,-0.0003,
     * -0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,-0.0002,
     * -0.0002,-0.0002,-0.0002,-0.0002,-0.0001,-0.0001,-0.0001,-0.0001,
     * -0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001,
     * -0.0001,-0.0001,-0.0001,-0.0001,-0.0001,-0.0001, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
     *  0.0000, 0.0000, 0.0000, 0.0000, 0.0000 /

#if defined (LONLAT)
      dimension xcart(nx,ny), ycart(nx,ny), zcart(nx,ny)
      save xcart, ycart, zcart
#elif defined (STEREO)
      parameter (nbox=1200)  ! 100
      dimension abox_e(-nbox:nbox,-nbox:nbox),
     *          abox_w(-nbox:nbox,-nbox:nbox)
      save abox_e, abox_w, ibox
#endif

c        If first call, set flexural length scale (crustlen) and
c        coefficient for central deflection (crustcoef). 
c        Also interpolate akei from akeic for 2D lithospheric calcs,
c        and set 3D cartesian (x,y,z) coordinates (*cart, saved) 
c        for computation of great-circle distances between points.

      if (first) then

        if (nx.gt.1 .and. ny.gt.1) then
          crustlen_e = ( cruststiff_e / (rhoman*grav) ) ** 0.25
          crustlen_w = ( cruststiff_w / (rhoman*grav) ) ** 0.25
          crustcoef_e = crustlen_e**2 / (2.*pi*cruststiff_e)
          crustcoef_w = crustlen_w**2 / (2.*pi*cruststiff_w)
        else
          crustlen_e = ( 4.*cruststiff_e / (rhoman*grav) ) ** 0.25
          crustlen_w = ( 4.*cruststiff_w / (rhoman*grav) ) ** 0.25
          crustcoef_e = crustlen_e**3 / (8.*cruststiff_e)
          crustcoef_w = crustlen_w**3 / (8.*cruststiff_w)
        endif

        do ix=0,nkei
          zx = ix*dkei
          ia = max (0, int(zx/dkeic))
          ib = min (nkeic, ia+1)
          zwei = (zx - ia*dkeic) / dkeic
          akei(ix) = (1.-zwei)*akeic(ia) + zwei*akeic(ib)
        enddo

c;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
c       xcart,etc, or abox,etc, not used below unless nx or ny > 1:
        if (nx.gt.1 .and. ny.gt.1) then
c;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#if defined (LONLAT)

        zz = pi/180.
        do j=1,ny
          do i=1,nx
            xcart(i,j) = cos(alatd(i,j)*zz)*cos(alond(i,j)*zz)
            ycart(i,j) = cos(alatd(i,j)*zz)*sin(alond(i,j)*zz)
            zcart(i,j) = sin(alatd(i,j)*zz)
          enddo
        enddo

#elif defined (STEREO)

c       set local effect array (abox)
c       (8*crustlen gets forebulge...see akeic above)       
        zdx = 0.5 * ( dx((nx+1)/2,(ny+1)/2) + dy((nx+1)/2,(ny+1)/2) )
        ibox = nint (8.*max(crustlen_e,crustlen_w)/zdx) + 1
c       write(6,*) 'bedrock: ibox = ',ibox
        if (ibox.gt.nbox) then
          write(6,'(a,i6,a,i6/a,2f8.3,a,f8.3)')
     *     '*** bedrock: ibox too large. ibox=',ibox, '  nbox=',nbox,
     *     '    crustlen_[e,w] (km)=', crustlen_e*.001, crustlen_w*.001,
     *     '    zdx (km)=', zdx*.001
          stop
        endif
        do ja = -ibox,ibox
          do ia = -ibox,ibox
            zx = sqrt(float(ia**2 + ja**2)) * zdx / crustlen_e
            ix = min (nint(zx/dkei), nkei)
            abox_e(ia,ja) = akei(ix)

            zx = sqrt(float(ia**2 + ja**2)) * zdx / crustlen_w
            ix = min (nint(zx/dkei), nkei)
            abox_w(ia,ja) = akei(ix)
          enddo
        enddo

#endif

c;;;;;;;;;;;;
        endif
c;;;;;;;;;;;;

        first = .false.

      endif
c>>>>>
#endif
c>>>>>

      do j=1,ny
        do i=1,nx
          zmass(i,j) =   rhoice*h(i,j)     + rholiq*hw(i,j)
     *                 + rhosed*sedim(i,j) - equiload(i,j)
        enddo
      enddo

c>>>>>>>>>>>>>>>>>>>>>
#if ! defined (IFLITH)
c>>>>>>>>>>>>>>>>>>>>>

c        No lithosphere:

      do j=1,ny
        do i=1,nx
          deflect(i,j) = -zmass(i,j) / rhoman
        enddo
      enddo

c>>>>
#else
c>>>>
 
c        Compute elastic lithospheric and equilibrium asthenospheric
c        displacement (deflect, negative downwards)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (nx.gt.1 .and. ny.gt.1) then
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c          2D: Brotchie and Sylvester

c          Loop over multiple latitude segments, saving deflections
c          due to each in deflectseg, for parallelization. 

        call zero (deflectseg, nx*ny*nseg)

#if defined (OMPLITH)
!$OMP PARALLEL DO
!$OMP* shared  (zmass, darea, alatd, alond,
!$OMP*          deflectseg, akei, dx, dy, xh, yh,
!$OMP*          crustlen_e, crustlen_w, 
!$OMP*          abox_e, abox_w, ibox,
!$OMP*          crustcoef_e, crustcoef_w
#  if defined (LONLAT)
!$OMP*         ,xcart, ycart, zcart
#  endif
!$OMP*         )
!$OMP* private (iseg, ja,jb, i,j, zload,zdefcen,zdx, ii,jj, zx,ix, zcen,
!$OMP*          wei_w, wei_e, crustcoef, crustlen, zabox,
!$OMP*          nbufd, j1, j2, nj12, jin, iin,
!$OMP*          jja, jjb, iia, iib, jbox, idist, idirb, iic
!$OMP*         )
!$OMP* schedule(dynamic,2)
#endif

c=====================
        do iseg=1,nseg
c=====================

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c~~~~~~~~~~~~~~~~~~~
#if defined (LONLAT)
c~~~~~~~~~~~~~~~~~~~
          ja = (iseg-1)*(ny/nseg) + 1
          jb =  iseg   *(ny/nseg)
          if (iseg.eq.nseg) jb = ny
          do j=ja,jb
            jin = j
            do i=1,nx
              iin = i
c~~~~
#else
c~~~~
c         Expand domain to include deflections outside nested domain.
c         Only works with STEREO, not with LONLAT. Used to do this for 
c         nesting only, via bmass, read from the  nestdrive file,
c         but it led to complications if driving run uses ZEROICEBEDEQ
c         and nested run doesn't, or v.v (different defns of equiload).
c         Here, expand indices by 500 km, and use load (zmass, area)
c         of closest point (iin,jin) within the domain (actually 
c         2nd closest, see below). 
c         Array abox and ii,jj logic can handle the out-of-domain i,j. 
c         For nesting, don't use loads at domain edges, because their 
c         h and hb are set from nestdrive file, inconsistent with 
c         topbedeq and equiload.

#  if defined (NESTING)
          nbufd = nint(50.e3/dd0)          ! changed for nesting 10/4/17
#  else
          nbufd = nint(500.e3/dd0)
#  endif

          j1 = 1-nbufd
          j2 = ny+nbufd
          nj12 = j2 - j1 + 1
          ja = (iseg-1)*(nj12/nseg) + j1
          jb =  iseg   *(nj12/nseg) + j1 - 1
          if (iseg.eq.nseg) jb = j2
          do j=ja,jb
            do i = 1-nbufd, nx+nbufd
#  if defined (NESTING)
              jin = max (2, min (ny-1, j))      ! don't use domain edges
              iin = max (2, min (nx-1, i))      ! don't use domain edges
#  else
              jin = max (1, min (ny, j))
              iin = max (1, min (nx, i))
#  endif
c~~~~~
#endif
c~~~~~

#if defined (EISANTA)
c             West vs East Antarctic, like in iceinit.F for geoflux:
              if ( alatd(iin,jin).gt. -86. .and.
     *             (alond(iin,jin).gt.170. .or. alond(iin,jin).lt.-30.) 
     *           ) then
                wei_w = 1.
                wei_e = 0.
              else
                wei_w = 0.
                wei_e = 1.
              endif
#else
              wei_w = 0.               ! if not EISANTA, *_w = *_e above
              wei_e = 1.
#endif
              crustcoef = wei_w*crustcoef_w + wei_e*crustcoef_e
              crustlen  = wei_w*crustlen_w  + wei_e*crustlen_e
 
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

              zload = zmass(iin,jin) * grav
c---
c             For diagnostic testing of deflect with dbx:
c             zcen = sqrt (   (2.*(jin-1.)/(ny-.9999) - 1.)**2
c    *                      + (2.*(iin-1.)/(nx-.9999) - 1.)**2 )
c             if (nx.eq.1) zcen = abs (2.*(jin-1.)/(ny-1) - 1.)
c             zload = 0.
c             if (zcen.lt.0.5) zload = rhoice*1000.*grav
cc            if (jin.gt.ny/2) zload = rhoice*1000.*grav
c--
c                 Add this point-load's deflection to deflect(i,j)

              if (abs(zload).gt.1000.) then    ! > ~10 cm of ice equiv

                zdefcen = zload*darea(iin,jin)*crustcoef !nb:times darea

#if defined (LONLAT)
c               lon-lat:
c               jbox = nint ((8.*crustlen)/dy((nx+1)/2,(ny+1)/2)) + 1
                jbox = nint ((4.*crustlen)/dy((nx+1)/2,(ny+1)/2)) + 1
                jja = max (1,  j-jbox)
                jjb = min (ny, j+jbox)
c               lat range equal to 8*crustlen (or pole):
                do jj=jja, jjb
c                 lon range outward from i until exceed 8*crustlen:
                  do idist=0,nx/2
                    idirb= 1
c                   only count i and i +/- nx/2 pts once:
                    if (idist.eq.0 .or. idist.eq.nx/2) idirb= -1
                    do idir=-1,idirb,2
                      iic = i + idir*idist
                      ii = mod (iic + nx - 1, nx) + 1
c                     either great circle (exact):
                      zx = acos ( min (1.,   xcart(i,j)*xcart(ii,jj)
     *                                     + ycart(i,j)*ycart(ii,jj)
     *                                     + zcart(i,j)*zcart(ii,jj)))
     *                     * radius / crustlen
c                     or 3-D straight line (approx):
c                     zx = sqrt(  (xcart(i,j)-xcart(ii,jj))**2
c    *                          + (ycart(i,j)-ycart(ii,jj))**2
c    *                          + (zcart(i,j)-zcart(ii,jj))**2 )
c    *                     * radius / crustlen
                      if (zx.gt.8.) goto 500
                      ix = min (nint(zx/dkei), nkei)
                      deflectseg(ii,jj,iseg) = deflectseg(ii,jj,iseg)
     *                                       + zdefcen*akei(ix)
                    enddo
                  enddo
  500             continue
                enddo
#elif defined (STEREO)
c               polar stereographic:
                jja = max (1,  j-ibox)
                jjb = min (ny, j+ibox)
                iia = max (1,  i-ibox)
                iib = min (nx, i+ibox)
                do jj=jja,jjb
                  do ii=iia,iib
                    zabox = wei_w * abox_w(ii-i,jj-j) 
     *                    + wei_e * abox_e(ii-i,jj-j) 
                    deflectseg(ii,jj,iseg) = deflectseg(ii,jj,iseg)
     *                                     + zdefcen*zabox
                  enddo
                enddo
#endif

              endif
c~~~~~~~~~~~~~~~~
            enddo
          enddo
c~~~~~~~~~~~~~~~~
c============
        enddo
c============

c          Sum up deflections from each segment (for parallelization)

        call zero (deflect, nx*ny)
        do iseg=1,nseg
          do j=1,ny
            do i=1,nx
              deflect(i,j) = deflect(i,j) + deflectseg(i,j,iseg)
            enddo
          enddo
        enddo

c%%%%%%%%%
      else
c%%%%%%%%%

c          1D: Turcotte and Schubert:

        crustcoef = crustcoef_e  ! if 1D, *w = *e above
        crustlen  = crustlen_e

        call zero (deflect, nx*ny)

        do j=1,ny
          do i=1,nx
            zload = zmass(i,j) * grav

#if defined (EISANTA) || defined (TRANSECTA)
c           West vs East Antarctic, like in iceinit.F for geoflux:
            if ( alatd(i,j).gt. -86. .and.
     *           (alond(i,j).gt.170. .or. alond(i,j).lt.-30.) 
     *         ) then
              wei_w = 1.
              wei_e = 0.
            else
              wei_w = 0.
              wei_e = 1.
            endif
#else
            wei_w = 0.               ! if not EISANTA, *_w = *_e above
            wei_e = 1.
#endif
            crustcoef = wei_w*crustcoef_w + wei_e*crustcoef_e
            crustlen  = wei_w*crustlen_w  + wei_e*crustlen_e

c              Add this point-load's deflection to deflect(i,j)

            if (abs(zload).gt.1000.) then      ! > ~10 cm of ice equiv
              zdefcen = zload*darea(i,j)*crustcoef   ! nb: times darea
              zdx = 0.5 * (dx(i,j) + dy(i,j))
              zdefcen = zload*zdx*crustcoef
              do ii=1,nx
                do jj=1,ny
                  zx =  sqrt (   (xh(i,j)-xh(ii,jj))**2
     *                         + (yh(i,j)-yh(ii,jj))**2 )
     *                  / crustlen
                  deflect(ii,jj) = deflect(ii,jj)
     *                           - zdefcen*exp(-zx)*(cos(zx)+sin(zx))
                enddo
              enddo
            endif

          enddo
        enddo

c%%%%%%%%%%
      endif
c%%%%%%%%%%

c>>>>>
#endif
c>>>>>

c        Local bedrock relaxation towards equilibrium

c-------------------------
      if (kasth.eq.0) then
c-------------------------

#if defined (NESTING)
        do j=2,ny-1
          do i=2,nx-1
#else
        do j=1,ny
          do i=1,nx
#endif
            za(i,j) = topbed(i,j) - topbedeq(i,j) - deflect(i,j)

#if defined (RWRAP)
            topbed(i,j) = topbed(i,j) - (dt/rwrap_tauasth)*za(i,j)
#elif defined (TAUASE)
c           change default asthenospheric time scale 
c--
c           over ASE: 
            if (alatd(i,j).gt. -78. .and. alatd(i,j).lt.-71. .and.
     *          alond(i,j).gt.-120. .and. alond(i,j).lt.-90.) then
c--
c           over WAIS: 
c           if ( alatd(i,j).gt. -86. .and.
c    *          (alond(i,j).gt.170. .or. alond(i,j).lt.-30.) ) then
c--
              topbed(i,j) = topbed(i,j) - (dt/TAUASE)*za(i,j)
            else
              topbed(i,j) = topbed(i,j) - (dt/taulocal)*za(i,j)
            endif
#else
            topbed(i,j) = topbed(i,j) - (dt/taulocal)*za(i,j)
#endif
          enddo
        enddo

c---------
      else
c---------

c          Diffusive (thin-channel) asthenospheric flow

        do j=1,ny
          do i=1,nx
            za(i,j) = topbed(i,j) - topbedeq(i,j) - deflect(i,j)
          enddo
        enddo

#if defined (NESTING)
        do j=2,ny-1
          do i=2,nx-1
#else
        do j=1,ny
          do i=1,nx
#endif
            jm = max(j-1,1)
            jp = min(j+1,ny)
            im = max(i-1,1)
            ip = min(i+1,nx)
              topbed(i,j) = topbed(i,j)
     *                 + dt * asthflow
     *                 * (    ((za(ip,j)-za(i,j))/dxu(i,j) ) * dyu(i,j)
     *                      + ((za(im,j)-za(i,j))/dxu(im,j)) * dyu(im,j)
     *                      + ((za(i,jp)-za(i,j))/dyu(i,j) ) * dxu(i,j)
     *                      + ((za(i,jm)-za(i,j))/dyu(i,jm)) * dxu(i,jm)
     *                   )  / darea(i,j)

          enddo
        enddo

c          Boundary conditions (equilibrium at boundaries)
    
#if ! defined (NESTING)
        if (nx.eq.1) then
          topbed(1,1)  = topbedeq(1,1)  + deflect(1,1)
          topbed(1,ny) = topbedeq(1,ny) + deflect(1,ny)
        else if (ny.eq.1) then
          topbed(1,1)  = topbedeq(1,1)  + deflect(1,1)
          topbed(nx,1) = topbedeq(nx,1) + deflect(nx,1)
        else
          do j=1,ny
            if (j.eq.1 .or. j.eq.ny) then
              iskip = 1
            else
              iskip = nx-1
            endif
            do i=1,nx,iskip
              topbed(i,j) = topbedeq(i,j) + deflect(i,j)
            enddo
          enddo
        endif
#endif

c----------
      endif
c----------

c        Apply long-term thermal subsidence, tectonic uplift

#if defined (FORCE34MYR) && defined (EISANTA) && defined (SEDIMENT)
      call subsidence (hb, topbed, topbedeq, equiload, timeice, dt)

      call uplift (topbed, topbedeq, timeice, dt)
#endif

c        Reset bed+sed elevation hb

      do j=1,ny
        do i=1,nx
          hb(i,j) = topbed(i,j) + sedim(i,j)
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

c{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
#if defined (FORCE34MYR) && defined (EISANTA)
c{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{

      subroutine subsidence (hb, topbed, topbedeq, equiload, timeice,dt)

c     Change topbed, topbedeq for long-term thermal subsidence
c     in West Antarctica. As in Fig. S2 of Wilson et al. (Palaeo3,2012)
c     for 34 Ma, linearly interpolating in time from files for 34, 22, 
c     11 Ma, from Doug Wilson pers. comm. 10/18.

c     Except: If NOSED, set topbedeq and equiload, interpolating in time
c     from Antarctic reconstructions for 34 Ma, 15 Ma, and modern
c     observed.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  hb(nx,ny), topbed(nx,ny), topbedeq(nx,ny), equiload(nx,ny)

#if defined (NOSED)

      dimension 
     *  topbedeq34(nx,ny), equiload34(nx,ny),
     *  topbedeq15(nx,ny), equiload15(nx,ny),
     *  topbedeqmo(nx,ny), equiloadmo(nx,ny)
      logical firstnosed
      data firstnosed /.true./
      save topbedeq34, equiload34, 
     *     topbedeq15, equiload15, 
     *     topbedeqmo, equiloadmo, 
     *     firstnosed
      dimension
     *  htmp(nx,ny), hbtmp(nx,ny), hwtmp(nx,ny), hstmp(nx,ny),
     *  maskwatertmp(nx,ny), maskcurtmp(nx,ny),
     *  indlaketmp(npoimax,nlakemax), npoilaketmp(nlakemax)

#else

#  include <comicepath.h>
#  include <netcdf.inc>
      character*120 cfile
      integer*4 ncid, rcode, varid, start(3), count(3)
      character*16 cvar
      parameter (nwilson=640, dwilson=10.e3,
     *           iwidwil=int(0.5*dx0/dwilson), 
     *           jwidwil=int(0.5*dy0/dwilson))
      real*4 subsidin(nwilson,nwilson)
      dimension subsid34(nx,ny), subsid22(nx,ny), subsid11(nx,ny),
     *          work(nx,ny), work2(nx,ny)
      logical firstsubsid
      data firstsubsid /.true./
      save firstsubsid, subsid34, subsid22, subsid11

#endif

c@@@@@@@@@@@@@@@@@@
#if defined (NOSED)
c@@@@@@@@@@@@@@@@@@

c=========================
      if (firstnosed) then
c=========================

c          If first call, read topogs for 34 Ma, 15 Ma, and modern, 
c          calculate equiloads, and save in topogeq[34,15,mo] and 
c          equiload[34,15,mo]    

c...................
        do iloop=1,3
c...................

          if (iloop.eq.1) then
c           Eocene-Oligocene (as selected in -D's):
            call shaldril_eoc (hbtmp)
            call zero (htmp, nx*ny)
            zslev = 56.

          else if (iloop.eq.2) then
c           15 Ma (as SHALMIOC15):
            call shaldril_mio (hbtmp)
            call zero (htmp, nx*ny)
            zslev = 56.

          else if (iloop.eq.3) then
c           Modern:
            call bedmap2_aggreg (htmp,  'thk')
            call bedmap2_aggreg (hbtmp, 'topg')
            zslev = 0.
          endif

          call scopy_i (nx*ny,maskcur,1,maskcurtmp,1) ! save (findwater)

          call findwater (maskwatertmp,indlaketmp,npoilaketmp,nlaketmp,
     *                    htmp, hbtmp, hwtmp, zslev, timeice)
          call adjustpres (maskwatertmp,indlaketmp,npoilaketmp,nlaketmp,
     *                     htmp, hbtmp, hwtmp, hstmp, zslev)

          call scopy_i (nx*ny,maskcurtmp,1,maskcur,1) ! save (findwater)

          do j=1,ny
            do i=1,nx
              if (iloop.eq.1) then
                topbedeq34(i,j) = hbtmp(i,j)
                equiload34(i,j) = rhoice*htmp(i,j) + rholiq*hwtmp(i,j)
              else if (iloop.eq.2) then
                topbedeq15(i,j) = hbtmp(i,j)
                equiload15(i,j) = rhoice*htmp(i,j) + rholiq*hwtmp(i,j)
              else if (iloop.eq.3) then
                topbedeqmo(i,j) = hbtmp(i,j)
                equiloadmo(i,j) = rhoice*htmp(i,j) + rholiq*hwtmp(i,j)
              endif
            enddo
          enddo

c............
        enddo
c............

        firstnosed = .false.

c==========
      endif
c==========

c        Linearly interpolate in time, 34 to 15 Ma or 15 Ma to 0 Ma,
c        to set current topbedeq,equiload:

      if (timeice.lt.-15.e6) then
        zwb = max (0., min (1., (timeice+34.e6) / (-15.e6 + 34.e6) ))
        zwa = 1.-zwb
        zwc = 0.
      else if (timeice.ge.-15.e6) then
        zwc = max (0., min (1., (timeice+15.e6) / (0. + 15.e6) ))
        zwb = 1.-zwc
        zwa = 0.
      endif

      do j=1,ny
        do i=1,nx
          topbedeq(i,j) = zwa*topbedeq34(i,j)
     *                  + zwb*topbedeq15(i,j)
     *                  + zwc*topbedeqmo(i,j)

          equiload(i,j) = zwa*equiload34(i,j)
     *                  + zwb*equiload15(i,j)
     *                  + zwc*equiloadmo(i,j)
        enddo
      enddo

c@@@@
#else
c@@@@

c        If first call, read in 3 Netcdf files for subsidence from
c        particular times to present, for 34, 22, and 11 Ma, from 
c        Doug Wilson, equiv. to Fig. S2 in Palaeo3 2012. Saved
c        in subsid[34,22,11].

      if (firstsubsid) then

        call zero (subsid34, nx*ny)
        call zero (subsid22, nx*ny)
        call zero (subsid11, nx*ny)

c-------------------
        do iwilr=1,3
c-------------------

        if (iwilr.eq.1) then
          cfile = cpathhome(1:lenchr(cpathhome)) //
     *            'Data_wilson/thsubs34max.grd'
        else if (iwilr.eq.2) then
          cfile = cpathhome(1:lenchr(cpathhome)) //
     *            'Data_wilson/thsubs22max.grd'
        else if (iwilr.eq.3) then
          cfile = cpathhome(1:lenchr(cpathhome)) //
     *            'Data_wilson/thsubs11max.grd'
        endif

        write (ioterm,'(2a)') 'Reading Wilson subsidence file ',
     *                        cfile(1:lenchr(cfile))
        rcode = nf_open(cfile, nf_nowrite, ncid)
        call checkrcode (rcode, cfile, 'open', 6000)
        cvar = 'z'
        start(1) = 1
        count(1) = nwilson 
        start(2) = 1
        count(2) = nwilson 
        rcode = nf_inq_varid (ncid, cvar, varid)
        call checkrcode (rcode, cfile, cvar, 7760)
        rcode = nf_get_vara_real (ncid, varid, start, count, subsidin)
        call checkrcode (rcode, cfile, cvar, 7761)
        rcode = nf_close (ncid)
        call checkrcode (rcode, cfile, 'close', 7762)

c          Transfer Wilson-grid subsidin to ice-grid subsid*, by ~tiling
c          (Wilson-grid x=0,y=0 is South Pole, as ice-grid xh=0,yh=0.
c           Wilson-grid (nwilson/2, nwilson/2) cell center is SP.)

        do j=1,ny 
          do i=1,nx 
            iw = nint(xh(i,j)/dwilson) + nwilson/2
            jw = nint(yh(i,j)/dwilson) + nwilson/2
            if (iw+iwidwil.ge.1 .and. iw-iwidwil.le.nwilson .and.
     *          jw+jwidwil.ge.1 .and. jw-jwidwil.le.nwilson) then
              zav = 0.
              nav = 0
              do jj = jw-jwidwil,jw+jwidwil 
                do ii = iw-iwidwil,iw+iwidwil 
                  if (ii.ge.1 .and. ii.le.nwilson .and.
     *                jj.ge.1 .and. jj.le.nwilson) 
     *              zav = zav + subsidin(ii,jj)
                  nav = nav + 1
                enddo
              enddo
              if (iwilr.eq.1) then
                subsid34(i,j) = zav/nav
              else if (iwilr.eq.2) then
                subsid22(i,j) = zav/nav
              else if (iwilr.eq.3) then
                subsid11(i,j) = zav/nav
              endif
            endif
          enddo
        enddo

c------------
        enddo
c------------

        call printmap (timeice, subsid34, 'subsid34', 50., iuout2d, 0)
        call printmap (timeice, subsid22, 'subsid22', 50., iuout2d, 0)
        call printmap (timeice, subsid11, 'subsid11', 50., iuout2d, 0)
        do j=1,ny
          do i=1,nx
            work(i,j)  = subsid34(i,j) - subsid22(i,j)
            work2(i,j) = subsid22(i,j) - subsid11(i,j)
          enddo
        enddo
        call printmap (timeice, work,  'subsid34-22', 50., iuout2d, 0)
        call printmap (timeice, work2, 'subsid22-11', 50., iuout2d, 0)

        firstsubsid = .false. 
      endif   

c old, based only on subsid34 and tau ~25.e6:
c     ztau = 25.e6
c     if (timeice.lt.-34.e6) then
c       zf = 0.             ! no subsid until 34 Ma (Wilson init. topog)
c     else if (timeice.lt.0.) then
c       zf = (exp(-(timeice+34.e6)/ztau)/(1.-exp(-34.e6/ztau))) / ztau
c     else 
c       zf = 0.
c     endif
c     do j=1,ny
c       do i=1,nx
c         topbed(i,j)   = topbed(i,j)   - dt*zf*subsid34(i,j)
c         topbedeq(i,j) = topbedeq(i,j) - dt*zf*subsid34(i,j)
c       enddo
c     enddo
c     return

c        Increment topbed, topedeq for current subsidence rate.
c        zr is current subsidence rate, assuming constant rates
c        between subsid<nn> file times. nb: file values are subsidence
c        (m, positive) from file time to present. 

      do j=1,ny
        do i=1,nx
          if (timeice.lt.-34.e6) then
            zr = 0.         ! no subsid until 34 Ma (Wilson init. topog)
          else if (timeice.lt.-22.e6) then
            zr = (subsid34(i,j)-subsid22(i,j)) / (34.e6-22.e6)
          else if (timeice.lt.-11.e6) then
            zr = (subsid22(i,j)-subsid11(i,j)) / (22.e6-11.e6)
          else if (timeice.lt.0.) then
            zr = (subsid11(i,j)-0.) / (11.e6-0.)
          else 
            zr = 0.
          endif

          topbed(i,j)   = topbed(i,j)   - dt*zr
          topbedeq(i,j) = topbedeq(i,j) - dt*zr
        enddo
      enddo

c@@@@@
#endif
c@@@@@

      return
      end

c-----------------------------------------------------------------------

      subroutine uplift (topbed, topbedeq, timeice, dt)

c     Change topbed, topbedeq for long-term tectonic uplift of 
c     Marie Bird Land, only for Shaldril "minimum" E-O reconstruction
c     (as in Wilson et al., P3, 2012, Fig. S5, based on LeMasurier and
c     Landis, 1996), or 0.5x for Shaldril "average". 

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  topbed(nx,ny), topbedeq(nx,ny)

      dimension upliftmbl(nx,ny)
      logical firstuplift
      data firstuplift /.true./
      save firstuplift, upliftmbl

      data alatup, alonup /-77., -133./
      data dlatup, dlonup /  2.,   19./

#if defined (SHALP311MIN) || defined (SHALP311AVE)

      if (firstuplift) then
c       elliptical area as in Fig. S5 of Wilson et al, P3, 2012:
        do j=1,ny
          do i=1,nx
            zr = sqrt ( ((alatd(i,j)-alatup)/dlatup)**2 + 
     *                  ((alond(i,j)-alonup)/dlonup)**2 )
#if defined (SHALP311MIN)            
            upliftmbl(i,j) = 1000. * max (1.-zr, 0.)
#elif defined (SHALP311AVE)            
            upliftmbl(i,j) =  500. * max (1.-zr, 0.)
#endif
          enddo
        enddo
        call printmap (timeice, upliftmbl, 'upliftmbl', 50., iuout2d, 0)
        firstuplift = .false.
      endif

      if (timeice.gt.-34.e6 .and. timeice.lt.0.) then
        do j=1,ny
          do i=1,nx
            topbed(i,j)   = topbed(i,j)   + dt*upliftmbl(i,j)/34.e6
            topbedeq(i,j) = topbedeq(i,j) + dt*upliftmbl(i,j)/34.e6
          enddo
        enddo
      endif

#endif

      return
      end

c{{{{{
#endif
c{{{{{

c-----------------------------------------------------------------------

c####################
#if defined (SLMODEL)
c####################

       subroutine bedrock_sl (h, hb, hw, topbed, sedim, dgeoid,
     *                        sealev, sealev_init, dt, dtimeice,timeice,
     *                        ncalls, itersl)

c        Steps bedrock(+sedim) elevation hb through one timestep, then 
c        adjusts topbed (below sedim). Uses Natalya Gomez/Jerry 
c        Mitrovica sea-level/earth model by calling script 
c        slwrapper (old), or slmodel.exe (new).

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

c        h      = ice-sheet thickness (m) (supplied)
c        hb     = bedrock+sedim elevation (m) (modified)
c        topbed = bedrock elevation, not including sedim (m) (modified)
c        sedim = sediment thickness (m) (supplied)
c        dt = bedrock-model timestep (yr) (supplied)
c        timeice = current time (yr) (supplied)
c        ncalls = number of calls this run (= calls to slmodel.exe),
c                 not counting first timestep call (iterbed=0)
c        itersl = outer-script iter # if -DSLMODELITER, 0 if not (supp)

      dimension 
     *  h(nx,ny),         hb(nx,ny),        hw(nx,ny),
     *  topbed(nx,ny),    sedim(nx,ny),     dgeoid(nx,ny)

#if defined (EISLINE) 
c     1-D flowline vs lat. nxsaxi= # of identical "longitude" values   
c     for each latitude that the sl model expects in its files:
c     parameter (nxs=1, nxsaxi=1024, nys=512)      ! sl model (lat) grid
      parameter (nxs=1, nxsaxi=   1, nys=512)      ! sl model (lat) grid
#else
c     2-D:
#  if defined (R512)
      parameter (nxs=1024, nys=512)            ! sl model (lon,lat) grid
#  elif defined (R256)
      parameter (nxs= 512, nys=256)            ! sl model (lon,lat) grid
#  endif
#endif
      dimension alons(nxs), alats(nys)
      save alons, alats

      dimension
     *  indl(nx,ny),   weil(nx,ny),  indb(nx,ny),   weib(nx,ny)
      save indl, weil, indb, weib

      dimension
     *  inds(nxs,nys), weis(nxs,nys), jnds(nxs,nys), wejs(nxs,nys)
      save inds, weis, jnds, wejs

      dimension
     *  hb_sl(nxs,nys), dhb_sl(nxs,nys), h_sl(nxs,nys),
     *  dgeo_sl(nxs,nys)

      logical ifexist

      character*200 cfile, cmd
      integer system

      data iterbed /-1/
      save iterbed

      iterbed = iterbed + 1

c        First call only, read sl grid lons and lats, and calculate
c        indices and weights between sl and ice grids. Also write
c        initial topography file for sl model

c===========================
      if (iterbed.eq.0) then
c===========================

c#if defined (EISLINE) 
c        do i=1,nx
c          write (191,*) alond(i,1)
c        enddo
c        do j=ny,1,-1
c          write (192,*) alatd(1,j)
c        enddo
c#endif
 
#if defined (EISLINE) 
        alons(1) = 0.                                  ! arbitrary value
#else
        write (cfile,'(2a,i0)') cpathslgrid (1:lenchr(cpathslgrid)),
     *      'Gridfiles/GLlon_', nys
        open (iuslg, file=cfile, form='formatted', status='old') 
        do i=1,nxs
          read (iuslg,*) alons(i)
          alons(i) = alons(i)*pi/180.
        enddo
        close (iuslg)
#endif

        write (cfile,'(2a,i0)') cpathslgrid (1:lenchr(cpathslgrid)),
     *      'Gridfiles/GLlat_', nys
        open (iuslg, file=cfile, form='formatted', status='old') 
        do j=nys,1,-1
          read (iuslg,*) alats(j)
          alats(j) = alats(j)*pi/180.
        enddo
        close (iuslg)

c          Calculate indices and weights to interpolate from the 
c          sl grid (nxs,nys) to the ice grid (nx,ny).
c          Returns [indl,weil] from the ice grid into the
c          sl longitude grid, and [indb,weib] from the ice grid into
c          the sl latitude grid.

        call setinterp (alons, alats, nxs, nys,
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, 0.)

c          Calculate indices and weights to interpolate from 
c          the ice grid (nx,ny) to the sl grid (nxs,nys).
c          Returns [inds,jnds,weis,wejs] from sl to ice grid.

        call setinterp2 (alons, alats, nxs, nys,
     *                   inds, weis, jnds, wejs)

c=====
c       do j=1,ny
c         do i=1,nx
c           write (175,'(2i4, 2f8.2, 2(i4,f6.3))') i, j, 
c    *                       alond(i,j), alatd(i,j),
c    *                       indl(i,j), weil(i,j),
c    *                       indb(i,j), weib(i,j)
c         enddo
c       enddo

c       do js=1,nys
c         do is=1,nxs
c           write (176,'(2i4, 2f8.2, 2(i4,f6.3))') is, js, 
c    *                       alons(is)*180./pi, alats(js)*180./pi,
c    *                       inds(is,js), weis(is,js),
c    *                       jnds(is,js), wejs(is,js)
c         enddo
c       enddo
c       stop
c=====

c          Assemble initial sl topography (null outside ice domain)

        do js = 1,nys
          do is = 1,nxs
            i  = inds(is,js)
#if defined (WRAPAROUND)
            ip1 = mod (i,nx) + 1
#else
            ip1 = min (i+1,nx)
#endif
            j  = jnds(is,js)
            jp1 = min (j+1,ny)
            wi = weis(is,js)
            wj = wejs(is,js)
            if (i.eq.0 .or. j.eq.0) then
              hb_sl(is,js) = 9999.
            else
              hb_sl(is,js) =
     *                wi *    wj  * hb(i  ,j  )
     *          + (1.-wi)*    wj  * hb(ip1,j  )
     *          +     wi *(1.-wj) * hb(i  ,jp1)
     *          + (1.-wi)*(1.-wj) * hb(ip1,jp1)
            endif
c           write (177,'(2i4, 2(i4,f6.3), 5f8.2)') is, js, 
c    *        i,wi, j, wj, 
c    *        hb(i,j), hb(ip1,j), hb(i,jp1), hb(ip1,jp1), hb_sl(is,js)
          enddo
        enddo

c          Write initial sl topography file

        cfile = 'AntarcticBedrock'
        open (iuslb, file=cfile, form='formatted', status='unknown') 
        close (iuslb, status='delete') 
        open (iuslb, file=cfile, form='formatted', status='new') 
        do j=nys,1,-1
#if defined (EISLINE)
          write (iuslb,*) (hb_sl(1,j), i=1,nxsaxi) 
#elif defined (EISANTA) || defined (NHA)
          write (iuslb,*) (hb_sl(i,j),i=1,nxs) 
c         write (178,'(1024i6)') (nint(hb_sl(i,j)), i=1,nxs)
#endif
        enddo
        close (iuslb)

c          Set initial dgeoid to zero
 
        call zero (dgeoid, nx*ny)

c==========
      endif
c==========

c        assemble current sl iceload 

      do js = 1,nys
        do is = 1,nxs
          i  = inds(is,js)
#if defined (WRAPAROUND)
          ip1 = mod (i,nx) + 1
#else
          ip1 = min (i+1,nx)
#endif
          j  = jnds(is,js)
          jp1 = min (j+1,ny)
          wi = weis(is,js)
          wj = wejs(is,js)

          if (i.eq.0 .or. j.eq.0) then
            h_sl(is,js) = 9999.
          else
            zfl          =
     *              wi *    wj  * hw(i  ,j  ) / max (hw(i  ,j  ),1.e-6)
     *        + (1.-wi)*    wj  * hw(ip1,j  ) / max (hw(ip1,j  ),1.e-6)
     *        +     wi *(1.-wj) * hw(i  ,jp1) / max (hw(i  ,jp1),1.e-6)
     *        + (1.-wi)*(1.-wj) * hw(ip1,jp1) / max (hw(ip1,jp1),1.e-6)

            if (zfl.lt.0.5-.001) then
c             grounded ice (>= 2 surrounding ice pts are grounded)
              h_sl(is,js)  =
     *              wi *    wj  * h (i  ,j  )
     *        + (1.-wi)*    wj  * h (ip1,j  )
     *        +     wi *(1.-wj) * h (i  ,jp1)
     *        + (1.-wi)*(1.-wj) * h (ip1,jp1)
            else
c             floating ice: set to zero for sl model
              h_sl(is,js) = 0.
            endif
          endif
        enddo
      enddo

c        write current sl iceload file 

      cfile = 'anticeload'
      open (iusli, file=cfile, form='formatted', status='unknown') 
      close (iusli, status='delete') 
      open (iusli, file=cfile, form='formatted', status='new') 
      do j=nys,1,-1
#if defined (EISLINE)
        write (iusli,*) (h_sl(1,j),i=1,nxsaxi)
#elif defined (EISANTA) || defined (NHA)
        write (iusli,*) (h_sl(i,j),i=1,nxs) 
c       write (179,'(1024i6)') (nint(h_sl(i,j)), i=1,nxs)
#endif
      enddo
      close (iusli)

c       call sea level-earth model

      ztime = timeice
      if (iterbed.gt.0) ztime = timeice + dtimeice

c     write (cmd,'(a, 4i6, 2f10.4, i6, i12)') 
c    *  'slwrapper.cli', itersl, iterbed, iterbed-1, ncalls,
c    *  -dt*iterbed*.001, -dt*(iterbed-1)*.001, nint(dt), nint(ztime)

      write (cmd,'(a, 4i6, 2f10.4, i6, i12, i6)') 
     *  '  slmodel.exe', itersl, iterbed, iterbed-1, ncalls,
     *  -dt*iterbed*.001, -dt*(iterbed-1)*.001, nint(dt), nint(ztime)

      write (ioterm,'(a/2x,a)') 
     *  'Calling slmodel.exe. cmd=', cmd(1:lenchr(cmd))

      ier = system (cmd) 

      if (ier.ne.0) then
        write (ioterm,'(a,i6)') 
     *    '*** Error: slmodel.exe returned ier =', ier
        call exit (-1)
      endif

      inquire (file='okend_sl', exist=ifexist)
      if (.not.ifexist) then
        write (ioterm,'(a)') 
     *    '*** Error: okend_sl not found on return from slmodel.exe'
        call exit (-1)
      endif

c     delete unused files in TPWV, SLtemp, output120p55:
      ier = system ('./docleansltempfiles') ! clean  

c        Read in new delta-bedrock topography (dhb_sl is a drop in hb),
c        but not first call (iterbed=0)...only exists for iterbed >= 1.
c        This is R-S in Gomez et al. 2012, GJI. 

c===========================
      if (iterbed.ge.1) then
c===========================
        cfile = 'bedrock'
        open (iuslb, file=cfile, form='formatted', status='old') 
        do j=nys,1,-1
#if defined (EISLINE)
          if (nxsaxi.eq.1) then
            read (iuslb,*) dhb_sl(1,j)
          else
            read (iuslb,*) dhb_sl(1,j), (zdum,i=2,nxsaxi) 
          endif
#elif defined (EISANTA) || defined (NHA)
          read (iuslb,*) (dhb_sl(i,j),i=1,nxs) 
#endif
        enddo
        close (iuslb)
c       call zero (dhb_sl, nxs*nys)                     ! for testing 1D

c          Also read in change in geoid (dgeo_sl, S in Gomez et al),
c          accumulated since start of run in dgeoid (written to history)

#if (defined (EISANTA) || defined (NHA)) && defined (DGEOID)
        cfile = 'dgeoid'
        open (iuslb, file=cfile, form='formatted', status='old') 
        do j=nys,1,-1
          read (iuslb,*) (dgeo_sl(i,j),i=1,nxs) 
        enddo
        close (iuslb)
#else
        call zero (dgeo_sl, nxs*nys)
#endif
 
c          interpolate and apply delta-bedrock topog to ice grid hb

        do j=1,ny
          do i=1,nx
            il = indl(i,j)
            ir = mod (il,nxs) + 1
            jb = indb(i,j)
            jt = min (jb+1, nys)
            wlb =     weil(i,j) *    weib(i,j)
            wrb = (1.-weil(i,j))*    weib(i,j)
            wlt =     weil(i,j) *(1.-weib(i,j))
            wrt = (1.-weil(i,j))*(1.-weib(i,j))
            hb(i,j) = hb(i,j) 
     *              - dhb_sl(il,jb)*wlb - dhb_sl(ir,jb)*wrb
     *              - dhb_sl(il,jt)*wlt - dhb_sl(ir,jt)*wrt
c           and same for delta-geoid:
            dgeoid(i,j) = dgeoid(i,j)
     *                  + dgeo_sl(il,jb)*wlb + dgeo_sl(ir,jb)*wrb
     *                  + dgeo_sl(il,jt)*wlt + dgeo_sl(ir,jt)*wrt
          enddo
        enddo

#if defined (EISANTA) || defined (NHA)
c       write diagnostic table of current global mean dgeo_sl, dgeoid 
        zgsl = 0.
        zasl = 0.
        do j=1,nys
          zwei = cos(alats(j))
          do i=1,nxs
            zgsl = zgsl + dgeo_sl(i,j)*zwei
            zasl = zasl +              zwei
          enddo
        enddo
        zgsl = zgsl/zasl
        zg = 0.
        za = 0.
        do j=1,ny
          do i=1,nx
            zwei = cos(alatd(i,j)*pi/180.)
            zg = zg + dgeoid(i,j)*zwei
            za = za +             zwei
          enddo
        enddo
        zg = zg/za
        write (156,'(a,i12,2f10.4)') 'time_dgeosl_dgeoid=',
     *    nint(timeice), zgsl, zg
#endif

#if defined (SLMODELRADIAL)
c       update sealev (change from init. passed by s.l. radial model):
        open (iusle, file = 'eustatic', form='formatted', status='old')
        read (iusle,'(f25.0)') dsealev
        close (iusle)
        sealev = sealev_init + dsealev
        write (6,'(a,i10,a,3f10.3)') 'bedrock_sl, yr:', nint(timeice),
     *    '  slinit, dsl, sealev=', sealev_init, dsealev,sealev
#endif

c==========
      endif
c==========

c        Reset bed (below sed) elevation topbed

      do j=1,ny
        do i=1,nx
          topbed(i,j) = hb(i,j) - sedim(i,j)
        enddo
      enddo

c        Write 1-D ice-grid variables to ascii file for plotting

#if defined (EISLINE)
      write(iuslout,'(f12.2)')    timeice + dtimeice
      write(iuslout,'(400f10.3)') (1.e-3*yh(1,j),j=1,ny)
      write(iuslout,'(400f10.3)') (h(1,j),j=1,ny)
      write(iuslout,'(400f10.3)') (hb(1,j)-sealev+hw(1,j)+h(1,j),j=1,ny)
      write(iuslout,'(400f10.3)') (hb(1,j)-sealev,j=1,ny)
      write(iuslout,'(400f10.3)') (hw(1,j),j=1,ny)
#endif

      return
      end

c-----------------------------------------------------------------------

      subroutine setinterp2 (alons, alats, nxs, nys,
     *                       inds, weis, jnds, wejs)

c          Calculate indices and weights to interpolate from the
c          ice grid (nx,ny) to the sea-level model grid (nxs,nys).
c          Returns [inds,weis,jnds,wejs], indices into ice grid.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  alons(nxs), alats(nys),
     *  inds(nxs,nys), weis(nxs,nys), jnds(nxs,nys), wejs(nxs,nys)


#if defined (EISLINE)

c        For 1-D, if grid sizes identical, assume grids are identical

      if (nxs.eq.nx .and. nys.eq.ny) then
        do js=1,nys
          do is=1,nxs
            inds(is,js) = is
            jnds(is,js) = js
            weis(is,js) = 1.
            wejs(is,js) = 1.
          enddo
        enddo
        return
      endif

c        For 1-D, if grid sizes not identical, lin interp vs latitude

      is = 1
      do js=1,nys 
        inds(is,js) = 0
        jnds(is,js) = 0
        weis(is,js) = 0.
        wejs(is,js) = 0.
        zlats = alats(js)*180./pi 
        if (zlats.ge.alat1 .and. zlats.le.alat2) then
          if (zlats.le.alatd(1,1)) then
            inds(is,js) = 1
            jnds(is,js) = 1
            weis(is,js) = 1.
            wejs(is,js) = 1.
          else if (zlats.ge.alatd(1,ny)) then
            inds(is,js) = 1
            jnds(is,js) = ny
            weis(is,js) = 1.
            wejs(is,js) = 1.
          else
            do j=2,ny
              if (zlats.le.alatd(1,j)) then
                inds(is,js) = 1
                jnds(is,js) = j-1
                weis(is,js) = 1.
                wejs(is,js) = (alatd(1,j)-zlats) / dlatd
                go to 500
              endif
            enddo
  500       continue
          endif
        endif
      enddo

#elif defined (EISANTA)

      do js=1,nys
        do is=1,nxs
          inds(is,js) = 0
          jnds(is,js) = 0
          weis(is,js) = 0.
          wejs(is,js) = 0.
          if (alats(js).lt.0.) then           ! southern hemisphere only
            call ps_xy (zx, zy, alats(js), alons(is), xoffa, yoffa)
            zxa = (zx + 0.5*nx*dx0)/dx0 + 0.5
            zya = (zy + 0.5*ny*dy0)/dy0 + 0.5
            ia = int(zxa)
            ja = int(zya)
            if (ia.ge.1 .and. ia.le.nx-1 .and. ja.ge.1 .and. ja.le.ny-1)
     *        then
              inds(is,js) = ia
              jnds(is,js) = ja
              weis(is,js) = 1. - (zxa-ia)
              wejs(is,js) = 1. - (zya-ja)
            endif
          endif
        enddo
      enddo

#elif defined (NHA)

      do js=1,nys
        do is=1,nxs
          inds(is,js) = 0
          jnds(is,js) = 0
          weis(is,js) = 0.
          wejs(is,js) = 0.

c            Bring global-grid longitude into -180 to 180 E range 
c            used in ice model 

          zlons = alons(is)
          if (zlons.gt.pi) zlons = zlons - 2.*pi
          zlats = alats(js)
c.....................................................................
          if ( ( (zlons.ge.alon1*pi/180. .and. zlons.le.alon2*pi/180.)
#if defined (WRAPAROUND)
     *           .or. .true.
#endif
     *         )
     *         .and. 
     *         (zlats.ge.alat1*pi/180. .and. zlats.le.alat2*pi/180.) 
     *       ) then
c.....................................................................

            if (zlons.le.alon(1)) then
#if defined (WRAPAROUND)
              inds(is,js) = nx
              weis(is,js) = (alon(1)-zlons) / dlond
#else
              inds(is,js) = 1
              weis(is,js) = 1.
#endif
            else if (zlons.ge.alon(nx)) then
#if defined (WRAPAROUND)
              inds(is,js) = nx
              weis(is,js) = 1. - (zlons-alon(nx)) / dlond
#else
              inds(is,js) = nx
              weis(is,js) = 1.
#endif
            else
              do i=2,nx
                if (zlons.le.alon(i)) then
                  inds(is,js) = i-1
                  weis(is,js) = (alon(i)-zlons) / dlond
                  go to 700
                endif
              enddo
  700         continue
            endif

            if (zlats.le.alat(1)) then
              jnds(is,js) = 1
              wejs(is,js) = 1.
            else if (zlats.ge.alat(ny)) then
              jnds(is,js) = ny
              wejs(is,js) = 1.
            else
              do j=2,ny
                if (zlats.le.alat(j)) then
                  jnds(is,js) = j-1
                  wejs(is,js) = (alat(j)-zlats) / dlatd
                  go to 750
                endif
              enddo
  750         continue
            endif

c..............
          endif
c..............

        enddo
      enddo

#endif

      return
      end

c#####
#endif
c#####

c-----------------------------------------------------------------------

c++++++++++++++++++++
#if defined (EISANTA)
c++++++++++++++++++++

      subroutine hypsom (hb, h, maskwater, sealev, timeice, 
     *                   iffirse, iflast)

c     Calculates bedrock hypsometry, for both current bed and 
c     isostatically rebounded no-ice bed (calcs done in calchypsom). 
c     If first call this run, saves results in hyp*firs. If last call 
c     this run, also calcs modern rebounded hypsoms, and writes all
c     hypsoms to unit iuhyp (cf. Wilson et al., Palaoe3, 2012, Fig. 6, 
c     for Antarctica).

#include <comicephys.h>
#include <comicegrid.h>

      dimension hb(nx,ny), h(nx,ny), maskwater(nx,ny) 

      parameter (dhyp=80., hypmin=-2500., hypmax=4000.)
      parameter (nhyp = nint((hypmax-hypmin)/dhyp))
      dimension hyp(nhyp),     hypnoi(nhyp),
     *          hypfirs(nhyp), hypnoifirs(nhyp),
     *          hypmod(nhyp),  hypnoimod(nhyp)
      save hypfirs, hypnoifirs, timefirs

c     for modern incl. rebounded (also done in subr bedrock):
      dimension
     *  htmp(nx,ny), hbtmp(nx,ny), hwtmp(nx,ny), hstmp(nx,ny),
     *  maskwatertmp(nx,ny), maskcurtmp(nx,ny),
     *  indlaketmp(npoimax,nlakemax), npoilaketmp(nlakemax)

c        Calculate current hypsometries (hyp, hypnoi)

      call calchypsom (hb, h, maskwater, sealev, 
     *                 hyp, hypnoi, nhyp, hypmin, hypmax, dhyp)

c        If first call this run, save hypsoms in hypfirs, hypnoifirs

      if (iffirse.eq.1) then
        do m=1,nhyp
          hypfirs(m) = hyp(m)
          hypnoifirs(m) = hypnoi(m)
          timefirs = timeice
        enddo
      endif

c        If last call this run, also calc modern hypsoms, 
c        and output saved (first) hypsoms, current (last) and moden

c--------------------------
      if (iflast.eq.1) then
c--------------------------
        call scopy_i (nx*ny,maskcur,1,maskcurtmp,1)   ! save (findwater)

c       modern:
        call bedmap2_aggreg (htmp,  'thk')
        call bedmap2_aggreg (hbtmp, 'topg')
        call findwater (maskwatertmp,indlaketmp,npoilaketmp,nlaketmp,
     *                  htmp, hbtmp, hwtmp, sealev, timeice)
        call adjustpres (maskwatertmp,indlaketmp,npoilaketmp,nlaketmp,
     *                   htmp, hbtmp, hwtmp, hstmp, sealev)
        call calchypsom (hbtmp, htmp, maskwatertmp, sealev,
     *                   hypmod, hypnoimod, nhyp, hypmin, hypmax, dhyp)

        call scopy_i (nx*ny,maskcurtmp,1,maskcur,1)   ! save (findwater)

        write (iuhyp,'(2(a,i12))') 
     *   'timefirs=', nint(timefirs), 'timelast=', nint(timeice)
        write (iuhyp,'(a,i6,a,f8.2)') 
     *    'nhyp=',nhyp, '  dhyp=', dhyp
        write (iuhyp,'(6x,4a)') 
     *    '   elev(m)',
     *    '     f_first  fnoi_first',
     *    '      f_last   fnoi_last',
     *    '    f_modern fnoi_modern'
        do m=1,nhyp
        write (iuhyp,'(i6,i10,6f12.6)') 
     *    m, nint(hypmin+(m-.5)*dhyp), 
     *    hypfirs(m), hypnoifirs(m), 
     *    hyp(m),     hypnoi(m), 
     *    hypmod(m),  hypnoimod(m)
        enddo

        call flush (iuhyp)
c----------
      endif
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine calchypsom (hb, h, maskwater, sealev,
     *                       hyp, hypnoi, nhyp, hypmin, hypmax, dhyp)

#include <comicephys.h>
#include <comicegrid.h>

      dimension hb(nx,ny), h(nx,ny), maskwater(nx,ny),
     *          hyp(nhyp), hypnoi(nhyp)

c        Calculate current hypsoms (hyp, hypnoi)

      call zero (hyp,    nhyp)
      call zero (hypnoi, nhyp)
      ztota    = 0.
      ztotanoi = 0.
      do j=1,ny
        do i=1,nx
#if defined (EISANTA)
          if (alatd(i,j).lt.-65.) then              ! only south of 65 S
#else
          if (.true.) then
#endif
            m = nint ((hb(i,j)-hypmin)/dhyp + 0.5)
            if (m.ge.1 .and. m.le.nhyp) then    ! ignore if out of range
              hyp(m) = hyp(m) + darea(i,j)
              ztota = ztota   + darea(i,j)
            endif

c           local isostatic rebound for removal of grounded ice:
            if (maskwater(i,j).eq.0.) then
              zov = rhoice*h(i,j) - rholiq*max(0.,sealev-hb(i,j))
              zhb = hb(i,j) + max(zov,0.)/rhoman
            else
              zhb = hb(i,j)
            endif
            m = nint ((zhb-hypmin)/dhyp + 0.5)
            if (m.ge.1 .and. m.le.nhyp) then    ! ignore if out of range
              hypnoi(m) = hypnoi(m) + darea(i,j)
              ztotanoi = ztotanoi   + darea(i,j)
            endif
          endif
        enddo
      enddo

      zmax    = -1.e20 
      zmaxnoi = -1.e20 
      do m=1,nhyp
        hyp(m)    = hyp(m)/ztota
        hypnoi(m) = hypnoi(m)/ztotanoi
        zmax      = max (zmax, hyp(m))
        zmaxnoi   = max (zmaxnoi, hypnoi(m))
      enddo

      do m=1,nhyp
        hyp(m) = (hyp(m)/zmax)*.02   ! so ~same y-scale as Wilson et al.
        hypnoi(m) = (hypnoi(m)/zmaxnoi)*.02  
      enddo

      return
      end

c+++++
#endif
c+++++

c-----------------------------------------------------------------------
