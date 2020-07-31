
c-----------------------------------------------------------------------

      subroutine sedflow (h, hb, hbsd, hw,
     *                    topbed, topbedeq, sedim,
     *                    t, baseperc, basefrml,
     *                    dbu, dbv, crhu, crhv, fsedu, fsedv, ub, vb,
     *                    maskwater, quarryrate, quarrycoef,quarrysecxy,
     *                    quarryacc, sedpres,
     *                    pelagrate, heats, heatb, dt)

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
c     t       = ice layer temperatures (for ice-sed interface) (K)(supp)
c     baseperc = basal liquid water percolation from above (m/y)
c     basefrml = basal in-situ freeze/melt rate (m/y)
c     df[u,v]= [e,n]ward basal stress (= driving) (N/m2) (u,v grds,supp)
c     [u,v]b = [e,n]ward basal(=sed) top veloc (m/yr) (u,v grids, supp)
c     crh[u,v] = overall basal coeff for icedyn (u,v grids, modified)
c     fsed[u,v]  = areal fraction of sed for icedyn iter (u,v grds, ret)
c     maskwater = mask (0=land, 1=ocean, +/-2,3...=lake #) (supplied)
c     quarryrate = quarrying rate (m/y) (returned)
c     pelagrate = pelagic deposition rate (m/y) (returned)
c     heats  = frictional heating in sediment layers (J/m2/a)(returned)
c     heatb  = frictional heating due to basal sliding (J/m2/a) (ret)
c     dt = sediment-model timestep (yr) (supplied)

      dimension 
     *  h(nx,ny),             hb(nx,ny),            
     *  hbsd(nx,ny),          hw(nx,ny), 
     *  topbed(nx,ny),        topbedeq(nx,ny),
     *  sedim(nx,ny),         t(nx,ny,0:nlevp),
     *  baseperc(nx,ny),      basefrml(nx,ny),
     *  dbu(0:nxp,0:nyp),     dbv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),    crhv(0:nxp,0:nyp),
     *  fsedu(0:nxp,0:nyp),   fsedv(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),      vb(0:nxp,0:nyp),
     *  maskwater(nx,ny),     quarryrate(nx,ny),  
     *  quarrycoef(nx,ny),    quarrysecxy(nx,ny),
     *  quarryacc(nx,ny),     sedpres(nx,ny),
     *  pelagrate(nx,ny),     heats(nx,ny,nsed),  heatb(nx,ny)

c     local:
c     fsed        = areal fraction of sed (h grid)
c     crhtmp      = sediment coeffs, temporary (h grid) 
c     seds[u,v]   = [e,n]ward sed flux (m2/yr) (u,v grids)
c     crhs[u,v]   = sediment coeffs, modif crhu,v (u,v grids)
c     taus[u,v]   = [e,n]ward sed-top stress (Pa) (u,v grids)
c     heats[u,v]  = [e,n]ward sed-top stress*veloc (J/m2/yr) (u,v grids)
c     heatb[u,v]  = [e,n]ward non-sed basal stress*veloc (" ") (" ")

      dimension
     *  fsed(nx,ny),           crhtmp(nx,ny),
     *  sedsu(0:nxp,0:nyp),    sedsv(0:nxp,0:nyp),
     *  tausu(0:nxp,0:nyp),    tausv(0:nxp,0:nyp),
     *  heatsu(0:nxp,0:nyp),   heatsv(0:nxp,0:nyp),
     *  heatbu(0:nxp,0:nyp),   heatbv(0:nxp,0:nyp)

      parameter (
#if defined (SEDENHANCE)
     *           enhancefac = SEDENHANCE,
#else
     *           enhancefac = 1.,       ! transport factor      !PD07,f1
#endif
 
#if defined (SEDM0)
     *           sedm0    = SEDM0,
#else
     *           sedm0    = 1.e10,      ! ref. viscosity (Pa.s) !PD07,f1
#endif

     *           sedd0    = 7.9e-7,     ! ref. deformation rate (1/s)
     *           sedphi   = 22.,        ! angle of internal friction
     *           sedc0    = 0 .,        ! sed cohesion (Pa) 
     *           sedfraca = 0.5 )       ! scale for frac.areal cover (m)

c    *           enhancefac = 10.,                              !PD07,f3
c    *           sedm0    = 1.e9,                               !PD07,f3
c    *           sedm0    = 3.e9,       ! ref. viscosity (Pa.s)
c    *           sedm0    = 3.e8,       ! ref. viscosity (Pa.s)
c    *           sedphi   = 2.,         ! angle of internal friction
c    *           sedphi   = .5,         ! angle of internal friction
c    *           pows     = 1.25,       ! exponent in sed.flow law
c                pows is a parameter in comicephys.h (= powb-1, powb=2)

c sedgamma (equiv. to crh for non-sed sliding) vs. sedm0:
c        sedm0            sedgamma
c 0.100000E+11        0.339446E-06
c 0.300000E+10        0.113149E-05
c 0.100000E+10        0.339446E-05

      save sedfac, sedbeta, sedgamma, first
      logical first
      data first /.true./

c....................
      if (first) then
c....................
c       convert from 1/s to 1/yr:
        sedfac = (86400.*365.) / (((2.*sedd0)**(pows-1))*(sedm0**pows))
        sedbeta  = (rhosed-rholiq) * grav * tan(sedphi*pi/180.)
        sedgamma = sedfac / ((pows+1.)*sedbeta)
        first = .false.
c..........
      endif
c..........

c        Set factor to reduce deformation if ice-sed base is frozen
c        or no basal melt (using same subr as in basecoef, basecoef_inv)

      do j=1,ny
        do i=1,nx
          call calc_crh (crhtmp(i,j), fracm,
     *                   sedgamma, 
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
c        fluxes seds[u,v], and taus[u,v], heats[u,v] for later heats 
c        calcs. Do this only where non-floating (as in basecoef 
c        for crh[u,v] - see comments there regarding g.l. and hw)

      zcrhnop = crhnos**(1./(pows+1.))

c        Eastward (u-grid):

      do j=1,ny
        do i=1,nx-1

          if (.not.(hw(i,j).gt.hwcut .and. hw(i+1,j).gt.hwcut)) then
            zdir = sign (1., -dbu(i,j))
            zma = 0.5*(1.+zdir)  ! ups
            zmb = 0.5*(1.-zdir)  ! ups
c           zma = 0.5            ! cen
c           zmb = 0.5            ! cen

            zcrhs=10.**(zma*log10(crhtmp(i,j))+zmb*log10(crhtmp(i+1,j)))
            fsedu(i,j) = zma*fsed(i,j) + zmb*fsed(i+1,j)

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

            sedsu(i,j) = sedsu(i,j) * enhancefac

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
            zma = 0.5*(1.+zdir)  ! ups
            zmb = 0.5*(1.-zdir)  ! ups
c           zma = 0.5            ! cen
c           zmb = 0.5            ! cen

            zcrhs=10.**(zma*log10(crhtmp(i,j))+zmb*log10(crhtmp(i,j+1)))
            fsedv(i,j) = zma*fsed(i,j) + zmb*fsed(i,j+1)

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

            sedsv(i,j) = sedsv(i,j) * enhancefac

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
      totzerosed = 0.

c        Set frictional heating vs layer (heats, on h-grid) for 
c        therm/vdif.

      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.10) then

c           undisturbed thickness:
            ztau = sqrt (   (0.5*(tausu(i,j)+tausu(i-1,j)))**2
     *                    + (0.5*(tausv(i,j)+tausv(i,j-1)))**2 )

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

c        Impose "cfl" limit to avoid -ve sedim (tracked by totzerosed) 

      zvmax = 0.25 * dd0 / dt

      do j=1,ny
        jp = min (j+1,ny)
        do i=1,nx
          ip = min (i+1,nx)
c         if (sedsu(i,j).gt.0.) then
c           sedsu(i,j) = min (sedsu(i,j),  zvmax*sedim(i ,j))
c         else
c           sedsu(i,j) = max (sedsu(i,j), -zvmax*sedim(ip,j))
c         endif
          sedsu(i,j) =  max (-zvmax*sedim(ip,j), min (zvmax*sedim(i,j),
     *                       sedsu(i,j) ))

c         if (sedsv(i,j).gt.0.) then
c           sedsv(i,j) = min (sedsv(i,j),  zvmax*sedim(i,j ))
c         else
c           sedsv(i,j) = max (sedsv(i,j), -zvmax*sedim(i,jp))
c         endif
          sedsv(i,j) =  max (-zvmax*sedim(i,jp), min (zvmax*sedim(i,j),
     *                       sedsv(i,j) ))
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
        enddo
      enddo

c        Under ice: local quarrying of bedrock, generating new sediment
c        (till) by basal sliding where sediment (large-scale) is thin, 
c        depending on heatb. quarrycoef is m/y per [tau(N/m2)*ub(m/y)].
c        nb: (1.e3/3.e7)/(1.e5*10.) = ~0.3e-10.
c        P+D 2003: 0.2e-9; P+D 2007, old LINED: 0.6e-9.

      call scopy (nx*ny, quarrysecxy, 1, quarrycoef, 1)
      call zero (quarryrate, nx*ny)

#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif

#if defined (FORCE34MYR) && !defined (NOSED)
c         fast quarrying (1.e-9) for remaining pre-40 Ma regolith on
c         initial land, or pre-40 Ma ocean sed on initial ocean
c         (as in icehydrol.F, subr basecoef, for crhmel; also see
c         iceinit.F, subr readcrhmel, for crhmel(sav)):
          if (sedpres(i,j).eq.0.) then
            zw = max (0., min(1.,(quarryacc(i,j)-regolithinit)/5. ))
            zq = 1.e-9
          else
            zw = max (0., min(1.,(quarryacc(i,j)-sedoceaninit)/5. ))
            zq = 1.e-9
          endif
          if (quarrycoef(i,j).eq.0.) then
            quarrycoef(i,j) = zw*quarrycoef(i,j) + (1.-zw)*zq
          else
            quarrycoef(i,j) = 10.** (     zw *log10(quarrycoef(i,j))
     *                              + (1.-zw)*log10(zq) )
          endif
#endif

          if (h(i,j).gt.0.) then
#if defined (QUARRYRATEJAM)
c           Jamieson et al., Geomorph., 2008, Eq. (8):
            zuvb = sqrt(  (0.5*(ub(i,j)+ub(i-1,j)))**2
     *                  + (0.5*(vb(i,j)+vb(i,j-1)))**2 )
            quarryrate(i,j) = quarrycoef(i,j)*h(i,j)*zuvb 
     *                      * (1.-fsed(i,j))
#else
            quarryrate(i,j) = quarrycoef(i,j)*heatb(i,j)
#endif

            sedim(i,j) = sedim(i,j) + quarryrate(i,j)*dt
            totquar = totquar + quarryrate(i,j)*dt*darea(i,j)

c           reduce bedrock topog, and equil topog (see notes 6/12/18):
            zdq = (rhosed/rhobed)*quarryrate(i,j)*dt
            topbed(i,j)   = topbed(i,j) - zdq
            topbedeq(i,j) = topbedeq(i,j) - (1.-(rhobed/rhoman))*zdq
          endif

        enddo
      enddo

c        Ice free: pelagic sediment deposition in ocean and lakes.
c        Different rates if ice shelf or not, and none in deep ocean.

      call zero (pelagrate, nx*ny)

      pelagrate1 = 0.
      pelagrate2 = 0.
c     pelagrate1 = .0001   ! m/yr
c     pelagrate2 = .0001   ! m/yr
c     pelagrate1 = .0005   ! m/yr  ! P+D 2007 Fig. 5
c     pelagrate2 = .0005   ! m/yr  ! P+D 2007 Fig. 5

#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
          if (maskwater(i,j).eq.1 .and. topbedeq(i,j).gt.-1999.) then
            if (h(i,j).eq.0.) then
              pelagrate(i,j) = pelagrate1
            else
              pelagrate(i,j) = pelagrate2
            endif
            sedim(i,j) = sedim(i,j) + pelagrate(i,j)*dt
            totpelag = totpelag + pelagrate(i,j)*dt*darea(i,j)
          endif
        enddo
      enddo

c        Set negligible or negative sedim to zero

      do j=1,ny
        do i=1,nx
          if (sedim(i,j).lt. 1.e-10) then
            totzerosed = totzerosed - sedim(i,j)*darea(i,j)
            sedim(i,j)   = 0.
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

      subroutine sedocean (sedim, hb, topbed, maskwater, dt)

c        Diffuses sediment, mimicking turbidity currents.
c        Linear diffusion wrt surface slope (hb), only between points 
c        that are underwater (maskwater ne 0), and only downslope
c        out of points that have sediment (sedim > 0). 

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  sedim(nx,ny),     hb(nx,ny), 
     *  topbed(nx,ny),    maskwater(nx,ny)

      dimension sednew(nx,ny)

#if defined (COEFTURB)
      parameter (coefturb = COEFTURB)  ! m^2/yr
#else
      parameter (coefturb = 1.e4)      ! m^2/yr,  ~(100.e3)^2 / 1.e6
#endif

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 4)   ! 4 for E-W/N-S, 8 for diagonals too

      totdump = 0.
      totslump = 0.

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

      call scopy (nx*ny, sedim, 1, sednew, 1)

#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
c--------------------------------------
          if (maskwater(i,j).ne.0) then
c--------------------------------------
            zds = 0.
c........................
            do look = 1,4     
c........................
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskwater(ii,jj).ne.0) then 
              
                if (ii.ne.i) then
                  zdx = abs(0.5*(xh(ii,jj) - xh(i,j)))
                  zdy = 0.5*(dy(ii,jj) + dy(i,j))
                else 
                  zdx = abs(0.5*(yh(ii,jj) - yh(i,j)))
                  zdy = 0.5*(dx(ii,jj) + dx(i,j))
                endif
                zsl  = (hb(ii,jj)-hb(i,j)) / zdx

c                  Only diffuse if "source" sedim > .01 m.
c                  Prevent making "source" sedim negative ("cfl").

                if (zsl.gt.0. .and. sedim(ii,jj).gt.0.01) then 
                  zds = zds 
     *                + min ( coefturb*abs(zsl)*zdy*dt, 
     *                        0.25*sedim(ii,jj)*darea(ii,jj) )
     *                / darea(i,j)  

                else if (zsl.lt.0. .and. sedim(i,j).gt.0.01) then 
                  zds = zds 
     *                - min ( coefturb*abs(zsl)*zdy*dt, 
     *                        0.25*sedim(i,j)*darea(i,j) )
     *                / darea(i,j)  
                endif

              endif
c................
            enddo
c................
            sednew(i,j) = sedim(i,j) + zds
c--------------
          endif
c--------------

        enddo
      enddo

c        Set new sedim, diagnostic totslump, adjust bed+sed elevation hb

      do j=1,ny
        do i=1,nx
          zsedold = sedim(i,j)
          sedim(i,j) = max (0., sednew(i,j))
          totslump = totslump + (sedim(i,j)-zsedold)*darea(i,j)
          hb(i,j) = topbed(i,j) + sedim(i,j)
        enddo
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine sedoceanold (sedim, hb, topbed, topbedeq, maskwater)

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
     *  sedim(nx,ny),     hb(nx,ny), 
     *  topbed(nx,ny),    topbedeq(nx,ny),
     *  maskwater(nx,ny)

#if defined (SLOPECRIT)
      parameter (slopecrit = SLOPECRIT / 80.e3) ! m/m
#else
      parameter (slopecrit = 600./80.e3)        ! m/m
c     parameter (slopecrit = 200./40.e3)        ! m/m
#endif
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
     *                     sedim(i,j) )
              zar2 = darea(i,j)/darea(iimax,jjmax)
              sedim(i,j)         = sedim(i,j)         - dsed
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
c#if defined (EISLINE) || defined (EISANTA) || defined (NHA) || defined (CARB) || defined (ORDO)
c#if defined (NESTING)
c      do j=2,ny-1
c        do i=2,nx-1
c#else
c      do j=1,ny
c        do i=1,nx
c#endif
c          if (maskwater(i,j).eq.1 .and. topbedeq(i,j).le.-1999.) then
c           totdump = totdump + sedim(i,j)*darea(i,j)
c           sedim(i,j) = 0.
c          endif
c        enddo
c      enddo
c#endif

c        Adjust bed+sed elevation hb

      do j=1,ny
        do i=1,nx
          hb(i,j) = topbed(i,j) + sedim(i,j)
        enddo
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine sedbudg (sedim, quarryacc, topbed, msector, timeice,dt,
     *                    iffirse, iflast, nyeartot, nyearsedbud, iloop)

c        Accumulates and resets domain-wide sediment budget quantities
c        (tot*, tot*a, timesedprev in comicephys).

#include <comicephys.h>
#include <comicegrid.h>

      dimension sedim(nx,ny), quarryacc(nx,ny), 
     *          topbed(nx,ny), msector(nx,ny)

c     totquar, etc (this timestep, without a) are in comicephys.h:
      save totquara, totpelaga, totzeroseda, totdumpa, totslumpa,
     *     totbedprev, totsedprev, timesedprev

      logical first
      save first
      data first /.true./

      if (nyearsedbud.eq.0) return

      call zero (totsedsec, nsedsec)
      call zero (totqacsec, nsedsec)
      do j=1,ny
        do i=1,nx
          totsedsec(msedsec(i,j)) = totsedsec(msedsec(i,j)) 
     *                            + sedim(i,j)*darea(i,j)
          totqacsec(msedsec(i,j)) = totqacsec(msedsec(i,j)) 
     *                            + quarryacc(i,j)*darea(i,j)
        enddo
      enddo

      totbed = 0.
      totsed = 0.
      totqac = 0.
      do j=1,ny
        do i=1,nx
          totbed = totbed + topbed(i,j)*darea(i,j)
          totsed = totsed + sedim(i,j)*darea(i,j)
          totqac = totqac + quarryacc(i,j)*darea(i,j)
        enddo
      enddo

c       Zero/reset accumulators first call (in comicephys.h), else 
c       accumulate domain quantities set in subrs during this timestep.
c       (Don't accumulate first call, since totsed set at end of 
c       timestep, not start..so budget accounting actually starts at 
c       2nd timestep)    

      if (first .or. iloop.eq.1) then   ! also do start of adaptdt stage
        totquara  = 0.
        totpelaga = 0.
        totzeroseda = 0.
        totdumpa  = 0.
        totslumpa = 0.
        totbedprev  = totbed
        totsedprev  = totsed
        timesedprev = timeice
        first = .false.
      else
        totquara  = totquara  + totquar
        totpelaga = totpelaga + totpelag
        totzeroseda = totzeroseda + totzerosed
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
      totzeroseda = totzeroseda / (totarea*dtsedbud)
      totdumpa  = totdumpa   / (totarea*dtsedbud)
      totslumpa = totslumpa  / (totarea*dtsedbud)
      totdsed   = (totsed - totsedprev) /  (totarea*dtsedbud)
      toterrs = totdsed - (totquara + totpelaga + totzeroseda -totdumpa)

c        Write to tabular output file

      iu = iusedbud

c        Write header line(s) if first write 

      if (iffirse.eq.1) then
        nwrite = nint(float(nyeartot)/float(nyearsedbud)) + 1
        write (iu,'(a,i8)') 'nwritesedbud=',nwrite
        write (iu,'(22a)') 
     *    '      time',
     *    (csedsec(m), m=1,nsedsec),'         ALL',
     *    '         bed         sed',
     *    '        dsed      quarry     pelagic     zerosed',
     *    '        dump       slump         sederr',
     *    (csedsec(m), m=1,nsedsec),'         ALL'
        write (iu,"( 10x, 9(9x,'km3'), 
     *               2(11x,'m'), 
     *               6(7x,'mm/ka'), 1(10x,'mm/ka'),
     *               9(9x,'km3')
     *             )")
      endif

      zz = 1.e6                              ! m/y to mm/ka
      write (iu,'(i10, 9i12, 8f12.3, e15.5, 9i12)')
     *  nint(timeice), 
     *  (nint(totsedsec(m)/1.e9),m=1,nsedsec), nint(totsed/1.e9), 
     *  totbed/totarea, totsed/totarea, 
     *  zz*totdsed, zz*totquara, zz*totpelaga, zz*totzeroseda,
     *  zz*totdumpa, zz*totslumpa, zz*toterrs,
     *  (nint(totqacsec(m)/1.e9),m=1,nsedsec), nint(totqac/1.e9)

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

      subroutine calcwilsonbasins

c     Set Wilson (P3, 2012) basin info: csedsec, nsedsec, alonsedsec,
c     totsedsecobs, msedsec, in comicegrid.h. Called once at start of 
c     run, via initphys.

#include <comicegrid.h>
#include <comicepath.h>

      parameter (nbasin=10)       ! number of Fig. 4 basin polygon files
      character*80 cbasfile(nbasin)
      character cfile*80, cpath*240
      dimension xbas(1000), ybas(1000)
      dimension arrbas(nx,ny), work(nx,ny)

c     Table 1 basins, combining Ross & Bell-Amun to match Fig. 4 basins:
      data csedsec /
     *  '         DML',
     *  '          EL',
     *  '     MRL-PEL',
     *  '      QML-WL',
     *  '      G5L-Ol',
     *  ' Ross-Bel-Am',    
     *  '    W Penins',
     *  '     Weddell' /

c     Table 1 longitude boundaries and obs. total sed, combining
c     Ross and Bell-Amun as above (Ross: 165-210, Bell-Amun: 210-280)

      data alonsedsec / 0.,30.,60.,94.,124.,165.,280.,310.,360./
c     data alonsedsec / 0.,30.,60.,94.,124.,165.,280.,315.,360./  !sedi9

      data totsedsecobs /                                 ! units are m3
#if defined (WILSONSEDMAX)
c Max:
     *  0.66e15, 0.77e15, 1.43e15, 0.66e15, 0.72e15, 
     *  5.60e15, 2.00e15, 4.40e15  /                    ! 5.60=3.00+2.60
#elif defined (WILSONSEDAVE)
c Min,Max average:
     *  0.60e15, 0.70e15, 1.30e15, 0.60e15, 0.60e15, 
     *  4.35e15, 1.50e15, 3.40e15  /                    ! 4.35=2.25+2.10
#else
c Min:
     *  0.54e15, 0.63e15, 1.17e15, 0.54e15, 0.48e15, 
     *  3.10e15, 1.00e15, 2.40e15  /                    ! 3.10=1.50+1.60
cc   *  3.10e15, 2.00e15, 1.40e15  /                    ! sedi9
#endif

c     names of Fig. 4 polygon files:
      data cbasfile /
     *  'eros015e.xy',
     *  'eros045e.xy',
     *  'eros085e.xy',
     *  'eros100e.xy',
     *  'eros115e.xy',
     *  'eros145e.xy',
     *  'erosprydz.xy',
     *  'erosross.xy',
     *  'rstrpenins.xy',
     *  'rstrweddell.xy' /

      call zero (work, nx*ny)

      cpath = cpathhome(1:lenchr(cpathhome))
     *        // 'Data_wilson/'

      write (ioterm, '(a)')
     *   'Reading Wilson et al. basin polygon files'

c        For each Fig. 4 polygon file, set arrbas to 1 (outside = 0)

c===================
      do ib=1,nbasin 
c===================
        cfile = cpath(1:lenchr(cpath)) // cbasfile(ib)
        open (iuwilbas, file=cfile, status='old')

        nvbas = 0
   10   read (iuwilbas,*, err=20, end=30) zx, zy 
        nvbas = nvbas + 1
        xbas(nvbas) = zx*1.e3                                  ! km to m
        ybas(nvbas) = zy*1.e3
        goto 10

   20   if (nvbas.eq.0) goto 10                    ! skip header line(s)
        write (ioterm,'(/2a)') '*** Error reading Wilson basin file ',
     *                         cbasfile(ib)(1:lenchr(cbasfile(ib)))
        stop

   30   close (iuwilbas)

        call zero (arrbas, nx*ny)
        call polygon_wais (arrbas, xh,yh,nx,ny, xbas,ybas,nvbas, 1.)

c          Accumulate basin #'s in work, warning if previous overlaps 

        do j=1,ny
          do i=1,nx
            if (arrbas(i,j).eq.1.) then
              if (work(i,j).ne.0.) then
                write (ioterm,'(a,5i4)') 
     *            '*** Warning (calcwilsonbasins):', 
     *            ib, i, j, nint(work(i,j)), nint(arrbas(i,j))
              endif
              work(i,j) = float(ib)
            endif
          enddo
        enddo
c==========
      enddo
c==========

      call printmap (timeice, work, 'Wilson basins orig', 1.,iuout2d,0)

c        Set final sector indices (msedsec, in comicegrid.h), 
c        the number of the Table 1 basin...combining Fig. 4's 
c        #3 (eros085) and #7 (erosprydz) into Table 1's 3 (MRL-PEL),
c        and Fig. 4's #4 (eros100) and #5 (eros115) into Table 1's
c        # 4 (QML-WL).

      do j=1,ny
        do i=1,nx
          ib = nint(work(i,j))
          if (ib.eq.1 .or. ib.eq.2) then
            ibt = ib 
          else if (ib.eq.3 .or. ib.eq.7) then
            ibt = 3 
          else if (ib.eq.4 .or. ib.eq.5) then
            ibt = 4
          else if (ib.eq.6) then
            ibt = 5
          else if (ib.eq.8) then
            ibt = 6
          else if (ib.eq.9) then
            ibt = 7
          else if (ib.eq.10) then
            ibt = 8
          else
            ibt = -100
          endif
          msedsec(i,j) = ibt
        enddo
      enddo

c        Fig. 4 basins only encompass terrestrial Antarctica.
c        Extend into ocean (whole domain) using longitude-line 
c        boundaries in Table 1.

      do j=1,ny
        do i=1,nx
          if (msedsec(i,j).eq.-100) then
            zlon = mod (alond(i,j)+360., 360.)
            do m=1,nsedsec
              if ( zlon.ge.alonsedsec(m-1) .and. zlon.le.alonsedsec(m)
     *           ) then
                msedsec(i,j) = m
                go to 50 
              endif
            enddo
   50       continue
          endif
        enddo
      enddo

c        maskcur (used in printmap) not set yet. msedsec also written
c        to iuout2d in 1st call to iceshow2d to see grounding line.  

      call itor (msedsec, work, nx*ny)
      call printmap (timeice, work, 'Wilson basins 1', 1., iuout2d, 0)

      return
      end

c-----------------------------------------------------------------------

#if defined (SEDTRACK)

      subroutine sedtracking_eul (hb, sedim, sedimold, 
     *                            timeice, dt, weirun, sealev, ifrest,
     *                            iffirse, iflast, nyearstart, nyeartot)

c        Tracks dates of original deposition (itrtop increases) in
c        existing sediment, using an Eulerian grid with uniform
c        vertical spacing (dtrack meters), but with origin  at base of
c        sediment. Track the field of "time-of deposition" on this grid,
c        in array sedtrack(i,j,0:ntrack). If sed thickness moves upwards
c        through a vertical grid pt in this timestep (deposition,
c        sedim increases through a grid pt, so itrtop increases), 
c        set to current time. If sedim falls through a grid pt 
c        (erosion), set back to vnulsedtrack (in comicegrid.h).

c        Vertical index of sedtrack (0 to ntrack) runs upwards from 
c        bottom of sed (opposite to other vertical coords in program:
c        zeta, zsed, zbed, which run downwards).

c        ntrack,dtrack,sedtrack(),itrtop() in common (comicegrid.h),
c        for ease in having sedtrack,itrtop on restart file. 
c        This is needed for restarts across time segments with 
c        dorunsheetshelf.    

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  hb(nx,ny), sedim(nx,ny), 
     *  sedimold(nx,ny)                   ! sedimold only for diagnostic

      character*80 cform, cformt
      parameter (zbot = -2500., ztop=3000.) 
      dimension work(nx,ny)

c        Initialize all sedtrack to null, or any initial sedim
c        to start time. Only for cold start - sedim,sedtrack,itrtop 
c        are on restart file.

      if (ifrest.eq.-1 .and. iffirse.eq.1) then
        do j=1,ny
          do i=1,nx
            if (sedim(i,j).eq.0.) then 
              itrtop(i,j) = -1
              do n=0,ntrack
                sedtrack(i,j,n) = vnulsedtrack
              enddo
            else 
              iz = sedim(i,j)/dtrack
              iz = min (iz, ntrack)
              itrtop(i,j) = iz
              do n=0,iz
                sedtrack(i,j,n) = timeice
              enddo
            endif
          enddo
        enddo
      endif

c         Adjust sedtracks due to changed sed thickness in this 
c         timestep (sedimold to sedim, due to advect in sedflow,
c         pelag, quarrying, erode/deposit, slump, dump) 

      do j=1,ny
        do i=1,nx
          izold = itrtop(i,j)
          iznew = sedim(i,j)/dtrack
          if (iznew.gt.ntrack) then
c           write (ioterm,*) 
c    *        '*** Error: sedim exceeds sedtrack limit of ',
c    *        (ntrack+1)*dtrack, ' meters at i,j=',i,j
c           stop
            iznew = min (iznew, ntrack)
          endif

          if (sedim(i,j).eq.0.) then
            do n = 0,ntrack
              sedtrack(i,j,n) = vnulsedtrack 
            enddo
            itrtop(i,j) = -1 
          else if (iznew.gt.izold) then
c           deposition through at least one sedtrack vertical grid pt:
            do n = izold+1,iznew
              sedtrack(i,j,n) = timeice 
            enddo
            itrtop(i,j) = iznew 
          else if (iznew.lt.izold) then
c           erosion through at least one sedtrack vertical grid pt:
            do n = izold,iznew+1,-1
              sedtrack(i,j,n) = vnulsedtrack 
            enddo
            itrtop(i,j) = iznew 
          endif

c         if (i.eq.68 .and. j.eq.27) then
c           iu = 777
c           write (iu,'(/5f15.3)') 
c    *        timeice, sedimold(i,j), sedim(i,j), 
c    *        hb(i,j),hb(i,j)-sedim(i,j)
c           write (iu,'(3i6)') izold, iznew, itrtop(i,j)
c           do n=0,max(0,itrtop(i,j))
c             write (iu,'(i6,2f15.3)') n, n*dtrack, sedtrack(i,j,n)
c           enddo
c         endif
       
        enddo
      enddo
      
      return
      end

#endif

c-----------------------------------------------------------------------
