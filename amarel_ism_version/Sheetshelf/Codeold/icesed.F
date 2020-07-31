
c-----------------------------------------------------------------------

      subroutine sedflow (h, hb, hbsd, hw,
     *                    topbed, topbedeq, sedim, sedimun,
     *                    disturb, t, baseperc, basefrml,
     *                    dbu, dbv, crhu, crhv, fsedu, fsedv,
     *                    maskwater, quarryrate,
     *                    heats, heatb, dt)

c        Calculates new sediment thickness (sedim), and modifies basal
c        flow coeffs due to sediment (crh[u,v], used in icdyn), 
c        following Jensen et al.,1995,Geomorph.,14,157-166 and
c        1996,JGR,191,8717-8728. Uses simple time-explicit numerics 
c        for sediment change.
c
c        If pows+1 (sed) = powb (basal) (set in comicephys), can 
c        combine crhu and crhsu rigorously weighted by fsedu, 
c        as done below...see notes 1/08.

#include <comicephys.h>
#include <comicegrid.h>

c     h     = ice-sheet thickness (m) (supplied)
c     hb    = bedrock (+ sediment) elevation (m) (modified, at end)
c     hbsd  = std dev bedrock (+ sediment) elev (m) (supp, for framp)
c     topbed = bedrock elevation (not including sed) (m) (supplied)
c     topbedeq= equil topbed (m) (supplied)
c     sedim = sediment thickness (m) (modified)
c     sedimun = undisturbed sediment thickness (m) (modified)
c     disturb = extent of disturbed sfc layer of sed (returned)
c     t       = ice layer temperatures (for ice-sed interface) (K)(supp)
c     baseperc = basal liquid water percolation from above (m/y)
c     basefrml = basal in-situ freeze/melt rate (m/y)
c     df[u,v]= [e,n]ward basal stress (= driving) (N/m2) (u,v grds,supp)
c     [u,v]b = [e,n]ward basal(=sed) top veloc (m/yr) (u,v grids, supp)
c     crh[u,v] = overall basal coeff for icedyn (u,v grids, modified)
c     fsed[u,v]  = areal fraction of sed for icedyn iter (u,v grds, ret)
c     maskwater = mask (0=land, 1=ocean, +/-2,3...=lake #) (supplied)
c     quarryrate= quarrying rate (m/y) (returned)
c     heats  = frictional heating in sediment layers (J/m2/a)(returned)
c     heatb  = frictional heating due to basal sliding (J/m2/a) (ret)
c     dt = sediment-model timestep (yr) (supplied)

      dimension 
     *  h(nx,ny),             hb(nx,ny),            
     *  hbsd(nx,ny),          hw(nx,ny), 
     *  topbed(nx,ny),        topbedeq(nx,ny),
     *  sedim(nx,ny),         sedimun(nx,ny),       
     *  disturb(nx,ny),       t(nx,ny,0:nlevp),
     *  baseperc(nx,ny),      basefrml(nx,ny),
     *  dbu(0:nxp,0:nyp),     dbv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),    crhv(0:nxp,0:nyp),
     *  fsedu(0:nxp,0:nyp),   fsedv(0:nxp,0:nyp),

     *  maskwater(nx,ny),     quarryrate(nx,ny),
     *  heats(nx,ny,nsed),    heatb(nx,ny)

c     local:
c     fsed        = areal fraction of sed (h grid)
c     fmelt       = reduction factor for frozen sed (top) (h grid)
c     seds[u,v]   = [e,n]ward sed flux (m2/yr) (u,v grids)
c     crhs[u,v]   = sediment coeffs, modif crhu,v (u,v grids)
c     taus[u,v]   = [e,n]ward sed-top stress (Pa) (u,v grids)
c     heats[u,v]  = [e,n]ward sed-top stress*veloc (J/m2/yr) (u,v grids)
c     heatb[u,v]  = [e,n]ward non-sed basal stress*veloc (" ") (" ")

      dimension
     *  fsed(nx,ny),           fmelt(nx,ny),
     *  sedsu(0:nxp,0:nyp),    sedsv(0:nxp,0:nyp),
     *  tausu(0:nxp,0:nyp),    tausv(0:nxp,0:nyp),
     *  heatsu(0:nxp,0:nyp),   heatsv(0:nxp,0:nyp),
     *  heatbu(0:nxp,0:nyp),   heatbv(0:nxp,0:nyp)

      parameter (sedd0    = 7.9e-7,     ! ref. deformation rate (1/s)
     *           sedm0    = 1.e10,      ! ref. viscosity (Pa.s)
c    *           sedm0    = 3.e9,       ! ref. viscosity (Pa.s)
c    *           sedm0    = 3.e8,       ! ref. viscosity (Pa.s)
     *           sedc0    = 0 .,        ! sed cohesion (Pa) 
c                pows is a parameter in comicephys.h (= powb-1)
     *           sedphi   = 22.,        ! angle of internal friction
c    *           sedphi   = 2.,   !777  ! angle of internal friction
c    *           sedphi   = 0.5,  !777  ! angle of internal friction
     *           sedfraca = 0.5 )       ! scale for frac.areal cover (m)

      save sedfac, sedbeta, sedgamma, first
      logical first
      data first /.true./

      if (first) then
c       convert from 1/s to 1/yr:
        sedfac = (86400.*365.) / (((2.*sedd0)**(pows-1))*(sedm0**pows))
        sedbeta  = (rhosed-rholiq) * grav * tan(sedphi*pi/180.)
        sedgamma = sedfac / ((pows+1.)*sedbeta)
        first = .false.
      endif

c        Set factor to reduce deformation if ice-sed base is frozen
c        (using same subr as in basecoef, basecoef_inv)

      call resetr (fmelt, nx*ny, 1.)

      do j=1,ny
        do i=1,nx
          call calc_crh (zcrhdum, fmelt(i,j), zfracb,
     *                   1.e10, 
     *                   t(i,j,nlevp) - (tmelt-dtmdh*h(i,j)),
     *                   baseperc(i,j)-basefrml(i,j),
     *                   h(i,j),
     *                   hbsd(i,j),
     *                   0)
        enddo
      enddo

c        Set areal fraction of sed cover

      do j=1,ny
        do i=1,nx
          fsed(i,j) = min (1., sedim(i,j)/sedfraca)
        enddo
      enddo

      call zero (fsedu,(nxp+1)*(nyp+1))
      call zero (fsedv,(nxp+1)*(nyp+1))
      call zero (disturb,nx*ny)
      call zero (sedsu,(nxp+1)*(nyp+1))
      call zero (sedsv,(nxp+1)*(nyp+1))
      call zero (tausu,(nxp+1)*(nyp+1))
      call zero (tausv,(nxp+1)*(nyp+1))
      call zero (heatsu,(nxp+1)*(nyp+1))
      call zero (heatsv,(nxp+1)*(nyp+1))
      call zero (heats, nx*ny*nsed)

c        Compute sediment coefficients crhs[u,v], and modify overall
c        overall coeffs crh[u,v] for icedyn. (Need to have just called
c        basecoef to set crhu for bare ground). Also compute sediment 
c        fluxes seds[u,v], and taus[u,v], heats[u,v] for later sedimun, 
c        heats calcs. Do this only where non-floating (as in basecoef 
c        for crh[u,v] - see comments there regarding g.l. and hw)

      zcrhnop = crhnos**(1./(pows+1.))

c        Eastward (u-grid):

      do j=1,ny
        do i=1,nx-1

          if (.not.(hw(i,j).gt.hwcut .and. hw(i+1,j).gt.hwcut)) then
            zdir = sign (1., -dbu(i,j))

            zm = 0.5*((1.+zdir)*fmelt(i,j) + (1.-zdir)*fmelt(i+1,j))!ups
c           zm = 0.5*(          fmelt(i,j) +           fmelt(i+1,j))!cen

            zf = 0.5*((1.+zdir)*fsed(i,j)  + (1.-zdir)*fsed(i+1,j)) !ups
c           zf = 0.5*(          fsed(i,j)  +           fsed(i+1,j)) !cen
            fsedu(i,j) = zf

c           zcrhs = zm*sedgamma
            zcrhs = 10.**((1.-zm)*log10(crhnos) + zm*log10(sedgamma))

            z1 = crhu(i,j) ** (1./(pows+1.))
            z2 = zcrhs     ** (1./(pows+1.))
            zdenom = max ( (1.-fsedu(i,j))*z2 + fsedu(i,j)*z1, zcrhnop )
            zdbs = dbu(i,j) * z1 / zdenom
            zdbs = sign ( max(0.,abs(zdbs)-sedc0), zdbs )

            zcrhb = crhu(i,j)
            zdbb  = dbu(i,j) * z2 / zdenom 

            crhu(i,j) = max ( zcrhb*zcrhs / (zdenom**(pows+1)),
     *                        crhnos )

            sedsu(i,j) = fsedu(i,j)
     *                   * zcrhs * (abs(zdbs)**(pows+1)) * (-zdbs)
     *                   / (sedbeta*(pows+2.))

            tausu(i,j) = -zdbs

            heatsu(i,j) =     fsedu(i,j)  * zcrhs *(abs(zdbs)**(pows+2))

            heatbu(i,j) = (1.-fsedu(i,j)) * zcrhb *(abs(zdbb)**(pows+2))

          else

c           adjust crhu underwater (crudely weighted) for icedyn iters:
            zf = 0.5*(fsed(i,j)+fsed(i+1,j))
            crhu(i,j) = zf*sedgamma + (1.-zf)*crhu(i,j)

          endif
        enddo
      enddo

c        Northward (v-grid):

      do j=1,ny-1
        do i=1,nx
          if (.not.(hw(i,j).gt.hwcut .and. hw(i,j+1).gt.hwcut)) then
            zdir = sign (1., -dbv(i,j))

            zm = 0.5*((1.+zdir)*fmelt(i,j) + (1.-zdir)*fmelt(i,j+1))!ups
c           zm = 0.5*(          fmelt(i,j) +           fmelt(i,j+1))!cen

            zf = 0.5*((1.+zdir)*fsed(i,j)  + (1.-zdir)*fsed(i,j+1)) !ups
c           zf = 0.5*(          fsed(i,j)  +           fsed(i,j+1)) !cen
            fsedv(i,j) = zf

c           zcrhs = zm*sedgamma
            zcrhs = 10.**((1.-zm)*log10(crhnos) + zm*log10(sedgamma))

            z1 = crhv(i,j) ** (1./(pows+1.))
            z2 = zcrhs     ** (1./(pows+1.))
            zdenom = max ( (1.-fsedv(i,j))*z2 + fsedv(i,j)*z1, zcrhnop )
            zdbs = dbv(i,j) * z1 / zdenom 
            zdbs = sign ( max(0.,abs(zdbs)-sedc0), zdbs )

            zcrhb = crhv(i,j)
            zdbb  = dbv(i,j) * z2 / zdenom 

            crhv(i,j) = max ( zcrhb*zcrhs / (zdenom**(pows+1)),
     *                        crhnos )

            sedsv(i,j) = fsedv(i,j)
     *                   * zcrhs * (abs(zdbs)**(pows+1)) * (-zdbs)
     *                   / (sedbeta*(pows+2.))

            tausv(i,j) = -zdbs

            heatsv(i,j) =     fsedv(i,j)  * zcrhs *(abs(zdbs)**(pows+2))

            heatbv(i,j) = (1.-fsedv(i,j)) * zcrhb *(abs(zdbb)**(pows+2))

          else

            zf = 0.5*(fsed(i,j)+fsed(i,j+1))
            crhv(i,j) = zf*sedgamma + (1.-zf)*crhv(i,j)

          endif
        enddo
      enddo

      totquar = 0.
      totpelag = 0.

c        Set new *undisturbed* sediment thicknesses (sedimun), 
c        ie, amount below the cutoff depth for sediment deformation,
c        and extent of disturbed sed (disturb). Also set frictional 
c        heating vs layer (heats, on h-grid) for therm/vdif.

      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.10) then

c           undisturbed thickness:
            ztau = sqrt (   (0.5*(tausu(i,j)+tausu(i-1,j)))**2
     *                    + (0.5*(tausv(i,j)+tausv(i,j-1)))**2 )
            disturb(i,j) = min (sedim(i,j), fmelt(i,j)*ztau/sedbeta)
            sedimun(i,j) = max (0., min (sedimun(i,j),
     *                                   sedim(i,j)-disturb(i,j)))

c           frictional heating (heats, J/m2/y):
            zheats = 0.5 * (heatsu(i,j) + heatsu(i-1,j))
     *             + 0.5 * (heatsv(i,j) + heatsv(i,j-1))
            do k = 1,nsed
              zz1 = sedim(i,j) * zsedm(k-1)
              zz2 = sedim(i,j) * zsedm(k)
              zzf = (   max (0., ztau - sedbeta*zz1)**(powb+1)
     *                - max (0., ztau - sedbeta*zz2)**(powb+1) )
     *              / max (ztau**(powb+1), 1.e-20) 
              heats(i,j,k) = zzf * zheats
            enddo

            heatb(i,j) = 0.5 * (heatbu(i,j) + heatbu(i-1,j))
     *                 + 0.5 * (heatbv(i,j) + heatbv(i,j-1))

          endif
        enddo
      enddo


c        Under ice: step sediment thickness (sedim) due to advection

#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
          sedim(i,j)= sedim(i,j)
     *              + (   sedsu(i-1,j)*dyu(i-1,j) - sedsu(i,j)*dyu(i,j)
     *                  + sedsv(i,j-1)*dxv(i,j-1) - sedsv(i,j)*dxv(i,j)
     *                ) * dt / darea(i,j)
          sedimun(i,j) = min (sedim(i,j), sedimun(i,j))
        enddo
      enddo

c        Under ice: local quarrying of bedrock, generating new sed(till)

      call zero (quarryrate, nx*ny)

#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
c         quarrying of bedrock under ice by basal sliding where 
c         sediment (large-scale) is thin, depending on heatb:
          if (h(i,j).gt.0.) then
#if defined (EISLINE) && defined (LINED)
c           zcoef = 0.2 e-9
            zcoef = 0.6 e-9
#else
c           zcoef = .0005 / (1.e5*10.)     ! m/y per <tau(N/m2)*ub(m/y)>
c           zcoef = .0010 / (1.e5*10.)
            zcoef = 0.2 e-9                ! as in P+D 2003
c           zcoef = 0.2 e-10
c           zcoef = 0.                     ! no quarrying source 777
#endif
            quarryrate(i,j) = zcoef*heatb(i,j)
            sedim(i,j) = sedim(i,j) + quarryrate(i,j)*dt
#if defined (EISLINE) && defined (LINED)
c           reduce bedrock topog:
            topbed(i,j)   = topbed(i,j)
     *                    - (rhosed/rhobed)*quarryrate(i,j)*dt
            topbedeq(i,j) = topbedeq(i,j)
     *                    - (rhosed/rhobed)*quarryrate(i,j)*dt
#endif
            totquar = totquar + quarryrate(i,j)*dt*darea(i,j)
          endif
        enddo
      enddo

c        Ice free: pelagic sediment deposition in ocean and lakes.
c        Different rates if ice shelf or not, and none in deep ocean.

#if defined (EISLINE) && defined (LINED)
      pelagrate1 = 0.
      pelagrate2 = 0.
#else
      pelagrate1 = 0.
      pelagrate2 = 0.
c     pelagrate1 = .0001   ! m/yr
c     pelagrate2 = .0001   ! m/yr
c     pelagrate1 = .0005   ! m/yr ! 777
c     pelagrate2 = .0005   ! m/yr ! 777
#endif
#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
          if (maskwater(i,j).eq.1 .and. topbedeq(i,j).gt.-1999.) then
            if (h(i,j).eq.0.) then
              sedim(i,j) = sedim(i,j) + pelagrate1*dt
              totpelag = totpelag + pelagrate1*dt*darea(i,j)
            else
              sedim(i,j) = sedim(i,j) + pelagrate2*dt
              totpelag = totpelag + pelagrate2*dt*darea(i,j)
            endif
          endif
        enddo
      enddo

c        Set negligible sedim to zero

      do j=1,ny
        do i=1,nx
          if (sedim(i,j).lt. 1.e-10) then
            sedim(i,j)   = 0.
            sedimun(i,j) = 0.
          endif
        enddo
      enddo

c        Adjust bed+sed elevation hb

      do j=1,ny
        do i=1,nx
          hb(i,j) = topbed(i,j) + sedim(i,j)
        enddo
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine sedocean (sedim, sedimun, h, hb, topbed, topbedeq,
     *                     maskwater)

c        "Slumps" sediment wherever dwnd slope (d[hb]/dx) exceeds a max 
c        limit (slopecrit). Only do this for points that are underwater
c        (maskwater ne 0). 
c        For each such point, find steepest-descent neighbor point
c        (nlook), and instantaneously transfers sediment to that 
c        neighboring point to reduce slope to slopecrit (or all
c        sediment is transferred). Sweep simply through domain, 
c        then iterate (nsweep). Also, apply "b.c." of no sediment
c        at domain edge, where sed is dumped to deep ocean.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  sedim(nx,ny),     sedimun(nx,ny),
     *  h(nx,ny),         hb(nx,ny), 
     *  topbed(nx,ny),    topbedeq(nx,ny),
     *  maskwater(nx,ny)

      parameter (slopecrit = 200./40.e3)     ! m/m
      parameter (nsweep = 4)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 8)   ! 4 for E-W/N-S, 8 for diagonals too

c####################
#if defined (NOSHELF)
c####################

c        If no ice shelf dynamics, set all marine sed to 1 m, return

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.1) then
            sedim(i,j) = 1.
          else
            sedim(i,j) = min (sedim(i,j),200.)
          endif
        enddo
      enddo
      return
c#####
#endif
c#####

      totdump  = 0.
      totslump = 0.

c=========================
      do isweep = 1,nsweep
c=========================
      nmove = 0

#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
          if (sedim(i,j).gt.0. .and. maskwater(i,j).ne.0) then

c              Find steepest-descent neighbor

            zsmax = -1.e20
            do look = 1,nlook     
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskwater(ii,jj).ne.0) then 
                zx = sqrt((xh(i,j)-xh(ii,jj))**2+(yh(i,j)-yh(ii,jj))**2)
                zs = (hb(i,j) - hb(ii,jj)) / max (1.e-20,zx)
                if (zs.gt.zsmax) then
                  zsmax = zs
                  iimax = ii 
                  jjmax = jj 
                  zxmax = zx
                endif
              endif
            enddo

c              If steepest slope exceeds limit, transfer sed

c---------------------------------------
            if (zsmax.gt.slopecrit) then
c---------------------------------------
              dsed = min ( 0.5*(zsmax - slopecrit)*zxmax, 
     *                     sedim(i,j) )                 ! all sed moved
c    *                     sedim(i,j)-sedimun(i,j) )    ! disturbed only
              zar2 = darea(i,j)/darea(iimax,jjmax)
              sedim(i,j)         = sedim(i,j)         - dsed
              sedimun(i,j)       = min (sedim(i,j), sedimun(i,j))
              sedim(iimax,jjmax) = sedim(iimax,jjmax) + dsed*zar2
              totslump = totslump + dsed*darea(i,j)

c                Adjust bed+sed elevations hb

              hb(i,j)         = topbed(i,j)         + sedim(i,j)
              hb(iimax,jjmax) = topbed(iimax,jjmax) + sedim(iimax,jjmax)

              nmove = nmove + 1
c----------------
            endif
c----------------

          endif
        enddo
      enddo

c        Skip out of domain-sweep iterations if nothing moved

      if (nmove.eq.0) goto 100
c==========
      enddo ! end of isweep loop
c==========
  100 continue

c       Apply "b.c.": sediment dumped to deep ocean outside domain

c remove at domain edge:
c#if defined (EISLINE)
c      totdump = totdump + sedim(nx,1)*darea(nx,1)
c      sedim(nx,1) = 0.
c#elif defined (EISANTA) || defined (NHA) || defined (CARB) || defined (ORDO)
c      do j=1,ny
c        iskip = nx-1
c        if (j.eq.1 .or. j.eq.ny) iskip = 1
c        do i=1,nx,iskip
c          totdump = totdump + sedim(i,j)*darea(i,j)
c          sedim(i,j) = 0.
c        enddo
c      enddo
c#endif

c remove everywhere in ocean where equil bathym < -2000m:
#if defined (EISLINE) || defined (EISANTA) || defined (NHA) || defined (CARB) || defined (ORDO)
#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
          if (maskwater(i,j).eq.1 .and. topbedeq(i,j).le.-1999.) then
            totdump = totdump + sedim(i,j)*darea(i,j)
            sedim(i,j) = 0.
            sedimun(i,j) = 0.
          endif
        enddo
      enddo
#endif

c        Adjust bed+sed elevation hb

      do j=1,ny
        do i=1,nx
          hb(i,j) = topbed(i,j) + sedim(i,j)
        enddo
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine sedbudg (sedim, timeice, dt, 
     *                    iffirse, iflast, nyeartot, nyearsedbud)

c        Accumulates and resets domain-wide sediment budget quantities
c        (tot*, tot*a, timesedprev in common)

#include <comicephys.h>
#include <comicegrid.h>

      dimension sedim(nx,ny)

      logical first
      save first
      data first /.true./

      if (nyearsedbud.eq.0) return

      totsed = 0.
      do j=1,ny
        do i=1,nx
          totsed = totsed + sedim(i,j)*darea(i,j)
        enddo
      enddo

c       Zero/reset accumulators first call, else accumulate domain 
c       quantities set in subrs during this timestep. (Don't accumulate
c       first call, since totsed set at end of timestep, not start...
c       do budget accounting actually starts at 2nd timestep)    

      if (first) then
        totquara  = 0.
        totpelaga = 0.
        totdumpa  = 0.
        totslumpa = 0.
        totsedprev  = totsed
        timesedprev = timeice
        first = .false.
      else
        totquara  = totquara  + totquar
        totpelaga = totpelaga + totpelag
        totdumpa  = totdumpa  + totdump
        totslumpa = totslumpa + totslump
      endif

c       If not sed budget time, return

      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyearsedbud),dt)).lt.dt
     *     .or. iffirse.eq.1 .or. iflast.eq.1
     *     )
     *   ) return 

c        Normalize accumulators 

      dtsedbud = max (timeice - timesedprev, 1.e-6)
      totquara  = totquara   / (totarea*dtsedbud)
      totpelaga = totpelaga  / (totarea*dtsedbud)
      totdumpa  = totdumpa   / (totarea*dtsedbud)
      totslumpa = totslumpa  / (totarea*dtsedbud)
      totdsed   = (totsed - totsedprev) /  (totarea*dtsedbud)
      toterr = totdsed - (totquara + totpelaga - totdumpa)

c        Write to tabular output file

      iu = iusedbud

c        Write header line(s) first write

      if (iffirse.eq.1) then
        nwrite = nint(float(nyeartot)/float(nyearsedbud)) + 1
        write (iu,'(a,i8)') 'nwritesedbud=',nwrite
        write (iu,'(3a)') 
     *    '      time',
     *    '      sed     dsed   quarry  pelagic     dump    slump',
     *    '     sederr'
      endif

      zz = 1.e6                              ! m/y to mm/ka
      write (iu,'(i10,6f9.3,f12.6)')
     *  nint(timeice), 
     *  totsed/totarea, zz*totdsed, zz*totquara, zz*totpelaga, 
     *  zz*totdumpa, zz*totslumpa, 
     *  max (-999., min (9999., zz*toterr))

      call flush (iu)

c       Zero/reset accumulators 

      totquara  = 0.
      totpelaga = 0.
      totdumpa  = 0.
      totslumpa = 0.
      totsedprev= totsed
      timesedprev = timeice

      return
      end

c-----------------------------------------------------------------------

#if defined (SEDTRACK)

      subroutine sedtracking_eul (h, hs, hb, 
     *                            sedim, sedimold, 
     *                            timeice, dt, weirun, sealev, 
     *                            iffirse, iflast, nyearstart, nyeartot)

c        Tracks dates of original deposition (sedim > sedimold) in
c        existing sediment, using an Eulerian grid with uniform
c        vertical spacing (dtr meters), but with origin  at base of
c        sediment. Track the field of "time-of deposition" on
c        this grid, in array sedtrack(i,j,0:ntr). If sed thickness moves
c        upwards through a vertical grid pt in this timestep (depostion,
c        sedim > grd pt, sedimold < grd pt), set to current time. If 
c        sedim falls through a grd pt (erosion), set back to valnul.

c        At output time, do vertical interp using the sedtrack array
c        to obtain elevations (m) of particular time contours, and 
c        write to output file, to be plotted by ginsedtrack.

c        Vertical index of sedtrack (1 to ntr) runs upwards from 
c        bottom of sed (opposite to other vertical coords in program:
c        zeta, zsed, zbed, which run downwards).

c        ntr,dtr,sedtrack(),timebot,itrtop() in common (comicegrid.h),
c        for ease in having sedtrack, timebot, itrtop on restart file. 
c        This is needed for restarts across time segments with 
c        dorunsheetshelf.    

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny), hs(nx,ny), hb(nx,ny), sedim(nx,ny), sedimold(nx,ny)

      parameter (valnul=-9999.)

c     dtrack is the interval (years) between particular-time contours
c     written to output file. Could be 1 ma or 100 ka, the finest
c     interval plotted in ginsedtrack (making output file smaller) .. 
      parameter (dtrack = 10000.)

c     nyearsedtrack is the output write interval (years)
      parameter (nyearsedtrack = 1000000) 
c     parameter (nyearsedtrack = 100000) 

      dimension work(nx)
      character*80 cform, cformt

      parameter (zbot = -2500., ztop=3000.) 

c        Initialize all sedtrack to null, except bottom of sed to 
c        start time (which will remain fixed at that value through run)

      if (iffirse.eq.1) then
        call resetr (sedtrack, nx*ny*(ntr+1), valnul)
        timebot = timeice
        do j=1,ny
          do i=1,nx
            itrtop(i,j) = 0
            sedtrack(i,j,0) = timebot
          enddo
        enddo
      endif

c         Adjust sedtracks due to changed sed thickness in this 
c         timestep (sedimold to sedim, due to advect in sedflow,
c         pelag, quarrying, erode/deposit, slump, dump) 

      do j=1,ny
        do i=1,nx
          izold = sedimold(i,j)/dtr
          iznew = sedim   (i,j)/dtr
          if (iznew.gt.ntr) then
            write (ioterm,*) 
     *        '*** Error: sedim exceeds sedtrack limit of ',
     *        (ntr+1)*dtr, ' meters at i,j=',i,j
            stop
          endif

          if (iznew.gt.izold) then
c           deposition through at least one sedtrack vertical grid pt:
            do n = izold+1,iznew
              sedtrack(i,j,n) = timeice 
            enddo
            itrtop(i,j) = iznew 
          else if (iznew.lt.izold) then
c           erosion through at least one sedtrack verrtical grid pt:
            do n = izold,iznew+1,-1
              sedtrack(i,j,n) = valnul 
            enddo
            itrtop(i,j) = iznew 
          endif
        enddo
      enddo

c        If output time (nyearsedtrack), write to sedtrack output file

      if (nyearsedtrack.eq.0) return

      if ( .not. 
     *     (mod(abs(timeice)+0.5*dt,max(float(nyearsedtrack),dt)).lt.dt
     *     )
     *   ) return 

      j = (ny+1)/2
      iu = iusedtrack
      cform  = '(a/151f12.5)'
      cformt = '(a,i4,i8,a/151f12.5)'

      if (iffirse.eq.1) then
        nwrite = nint(float(nyeartot)/float(nyearsedtrack)) + 1
        write (iu,'(a,i8)') 'nwritesedtrack=',nwrite
      endif

      ntrack = nint((timeice-nyearstart)/dtrack) + 1

c        Write 0-D quantities
        
      write (iu,'(80("="))')
      write (iu,'(a,i9 / a,f9.4 , 3(/a,i9) , 5(/a,f9.3))')
     *  'time     =',  nint(timeice),
     *  'weirun   =',  weirun,
     *  'nx       =',  nx,
     *  'ntrack   =',  ntrack,
     *  'dtrack   =',  nint(dtrack),
     *  'xlef     =',  0.,
     *  'xrit     =',  xh(nx,j)*1.e-3,
     *  'zbot     =',  zbot,
     *  'ztop     =',  ztop,
     *  'sealev   =',  sealev

c        Write 1-D quantities (h-grid)

      write(iu,cform) 'xh       :', (xh(i,j)*1.e-3,          i=1,nx)
      write(iu,cform) 'h        :', (h(i,j),                 i=1,nx)
      write(iu,cform) 'hs       :', (hs(i,j),                i=1,nx)
      write(iu,cform) 'hb       :', (hb(i,j),                i=1,nx)
      write(iu,cform) 'sedim    :', (sedim(i,j),             i=1,nx)

c        Write in reverse order so appears "right way up" in output file

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do itrack = ntrack,1,-1
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        timtrack = (itrack-1)*dtrack

c          interp elevs of fixed times from sedtrack

        do i=1,nx
          if (timtrack.le.sedtrack(i,j,0)) then
            work(i) = hb(i,j) - sedim(i,j)
          else if (timtrack.gt.sedtrack(i,j,itrtop(i,j))) then
            work(i) = hb(i,j)
          else
            do k = 1,itrtop(i,j)
              if (timtrack.le.sedtrack(i,j,k)) then
                koa = k-1
                kob = k
                wob = (timtrack - sedtrack(i,j,koa))
     *              / max (sedtrack(i,j,kob) - sedtrack(i,j,koa), .001)
                work(i) = (koa + wob)*dtr + hb(i,j) - sedim(i,j)
                goto 100
              endif
            enddo
            write (ioterm,*) '*** sedtrack: should never get here'
  100       continue           
          endif
        enddo

        iyr = nint(timtrack)
        write(iu,cformt) 'sedtrack ', itrack, nint(timtrack), ' :', 
     *                   (work(i),i=1,nx)

c>>>>>>>>>>
      enddo
c>>>>>>>>>>

      call flush (iu)

      return
      end

#endif

c-----------------------------------------------------------------------
