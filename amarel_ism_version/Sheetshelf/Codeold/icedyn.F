
c-----------------------------------------------------------------------

      subroutine icedyn (h, hs, hb, hbsd, hw, sedim, t,
     *                   arhap, s1a0, s1b0, s2a0, s2b0, heath, heatb, 
     *                   budgall, budgneg, cliffmelt, calvice, crevliq,
     *                   ifnexo, fedge, hedge, maskh, visc,
     *                   sih, siu, siv, dfu, dfv, dbu, dbv,
     *                   ua, va, uia, via, ub, vb, uit, vit, uadv, vadv,
     *                   hu, hv, uw, vw, masku, maskv, 
     *                   crhu, crhv, fsedu, fsedv,
     *                   fracgu, fracgv, 
     *                   muind, mvind, thetau, thetav, 
     *                   hgu, hgv, qgu, qgv,
     *                   uschoof, vschoof, fschoofu, fschoofv,
     *                   uagrid, vagrid,
     *                   fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                   angnormu, angnormv,
#if defined (MELANGE)
     *                   ucimel, vcimel, hucimel, hvcimel, 
     *                   bpumel, bpvmel, thetaumel, thetavmel, hmel,
#endif
     *                   numh, nuvtot, itera, iterc, totflow,
     *                   maskwater, indlake, npoilake, nlake, sealev,
     *                   timeice, dt, ifrest, iffirse, iffirsr)

c     Steps combined ice sheet-shelf dynamics.
c     Shelf dyanimcs (SSA) solves for ua,va (dotria, dosparse[sor]). 
c     Sheet dynamics (SIA) solves for uia,via, and also time-implicitly
c     steps ice thickness h including SSA ua,va (icestep).  
c     ua,va,uia,via are passed, saved (in main prog) between iterations.
c 
c     SSA solutions really solve for ua = ub+uia, given uia: 
c        LHS(ua) = driving stress + C*(ua-uia).
c     ub are generally just diagnostic, set to ua-uia.
c     Where the SSA solution itself is needed (schoofgl, icestep,
c     ua-uia or uagrid-uia is used). 
c
c     Imposed Schoof fluxes and velocities are assumed (imagined) to 
c     include both stretching and internal deformation across g.l.,
c     i.e., uschoof are equiv. to ua, not ub, and are imposed on
c     the ua solution in dotria/dosparse[sor].
c
c     See notes 1/28/12 for flow diagram.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

c     passed:
      dimension 
     *  h(nx,ny),              hs(nx,ny),
     *  hb(nx,ny),             hbsd(nx,ny),       hw(nx,ny),
     *  sedim(nx,ny),          t(nx,ny,0:nlevp),  arhap(nx,ny),
     *  s1a0(nx,ny),           s1b0(nx,ny),
     *  s2a0(nx,ny),           s2b0(nx,ny), 
     *  heath(nx,ny),          heatb(nx,ny),
     *  budgall(nx,ny),        budgneg(nx,ny),
     *  cliffmelt(nx,ny),      calvice(nx,ny),    crevliq(nx,ny),
     *  ifnexo(nx,ny),         fedge(nx,ny),      hedge(nx,ny),   
     *  maskh(nx,ny),          visc(nx,ny),
     *  sih(nx,ny),            siu(0:nxp,0:nyp),  siv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),
     *  dbu(0:nxp,0:nyp),      dbv(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),      via(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),       vb(0:nxp,0:nyp),
     *  uit(0:nxp,0:nyp),      vit(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),     vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  uw(0:nxp,0:nyp),       vw(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp), 
     *  crhu(0:nxp,0:nyp),     crhv(0:nxp,0:nyp),
     *  fsedu(0:nxp,0:nyp),    fsedv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),   fracgv(0:nxp,0:nyp), 

     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp), 
     *  thetau(nx,ny),         thetav(nx,ny),
     *  hgu(0:nxp,0:nyp),      hgv(0:nxp,0:nyp),
     *  qgu(0:nxp,0:nyp),      qgv(0:nxp,0:nyp),
     *  uschoof(0:nxp,0:nyp),  vschoof(0:nxp,0:nyp),
     *  fschoofu(0:nxp,0:nyp), fschoofv(0:nxp,0:nyp),
     *  uagrid(0:nxp,0:nyp),   vagrid(0:nxp,0:nyp),

     *  fluxgrdu(0:nxp,0:nyp), fluxgrdv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp), fluxschv(0:nxp,0:nyp),
     *  angnormu(0:nxp,0:nyp), angnormv(0:nxp,0:nyp),
     *  maskwater(nx,ny),      indlake(npoimax,nlakemax), 
     *  npoilake(nlakemax)
#if defined (MELANGE)
     * ,ucimel(0:nxp,0:nyp),    vcimel(0:nxp,0:nyp),
     *  hucimel(0:nxp,0:nyp),   hvcimel(0:nxp,0:nyp),
     *  bpumel(0:nxp,0:nyp),    bpvmel(0:nxp,0:nyp),
     *  thetaumel(nx,ny),       thetavmel(nx,ny),     hmel(nx,ny)
#endif

c     local:

      dimension 
     *  h0(nx,ny),             h1(nx,ny),
     *  coef1(nx,ny),
     *  dhsu(0:nxp,0:nyp),     dhsv(0:nxp,0:nyp),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  medgeu(0:nxp,0:nyp),   medgev(0:nxp,0:nyp),
     *  lhsu(0:nxp,0:nyp),     lhsv(0:nxp,0:nyp),
c       for icestep:
     *  hu0(0:nxp,0:nyp),      hv0(0:nxp,0:nyp),
     *  dhsu0(0:nxp,0:nyp),    dhsv0(0:nxp,0:nyp),
     *  dfu0(0:nxp,0:nyp),     dfv0(0:nxp,0:nyp)
      real lhsu, lhsv

c     for basal heating calculations:
      dimension 
     *  heatbu(0:nxp,0:nyp),   heatbv(0:nxp,0:nyp)

c     cross grid:
      dimension
     *  coef2(0:nx,0:ny),      maskc(0:nx,0:ny)

c     for iteration (a,c) convergence criteria:
      dimension
     *  uaprev (0:nxp,0:nyp),  vaprev (0:nxp,0:nyp),     ! loop a
     *  h1prev(nx,ny)                                    ! loop c

c     for Hindmarsh and Payne (1996) convergence, loop c:
c     dimension delh1(nx,ny), delh1prev(nx,ny)

      parameter (delcrita = .01)   ! m/a  ! .01,.03 
      parameter (delcritc = .01)   ! m    ! .01,.03

#if defined (NOSHELF)
      parameter (nitera = 1)
      parameter (niterc = 1)
#else
#  if defined (NITERA)
      parameter (nitera = NITERA)
#  else
      parameter (nitera = 2)
#  endif
#  if defined (NITERC)
      parameter (niterc = NITERC)
#  else
      parameter (niterc = 2)
#  endif
#endif
      parameter (niterasch = nitera)

#if defined (SUBGRIDEDGE) || ! defined (NOSUBGRIDEDGESHELF)
      parameter (nlook = 4)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /
#endif

c-------------------------
c Start of executable code
c-------------------------
    
#if defined (TESTMEL2D) && defined (TESTMELPROF)
c     for 1-D melange flowline, just set ice in 1st grid point, return: 
      call zero (h, nx*ny)
      call zero (ua, (nxp+1)*(nyp+1))
      do j=1,ny
c       h(1,j) = -hb(1,j)/rhor
        h(1,j)  = 500. ! 1000. Jakob
        ua(1,j) = 5.e3 ! 10.e3. Jakob
      enddo
      return
#endif

c       Set basal coefficient cutoff to select points for inclusion
c       in elliptic SSA (shelf) eqns, or sheet only (m[u,v]ind below)

#  if defined (CRHUCRIT)
      crhucrit = CRHUCRIT * (1.e5**2)/(1.e5**powb)
#  elif defined (EISLINE) || defined (TEST2D)
      crhucrit = 0.       ! all points included
#  else
      crhucrit = 1.e-8    * (1.e5**2)/(1.e5**powb) ! old SHELFCRHU8
c     1.e-6 is ~stream for powb=2 (see basecoef, 1.e5 is "taunorm"):
c     crhucrit = 1.e50    ! no grounded points included 
#  endif

c       Set |ua|,|va| speed limit (after dosparse[sor], and in schoofgl)

#if defined (SPEEDLIMIT)
      if (SPEEDLIMIT.gt.0.) then
        vmax = SPEEDLIMIT                         ! 20.e3 ! 10.e3 ! 5.e3
      else if (SPEEDLIMIT.eq.0.) then
        vmax = 0.75 * dd0 / dt
      else
        vmax = 1.e20
      endif
#else
      vmax = 0.75 * dd0 / dt
#endif

c       Set dynamic mask for no flow *out* of grid boxes with fedge < 1.
c       fedge has been set for grounded pts in climate_all, and for 
c       floating in subr ocean. Implemented here just in the advection 
c       of h, in icestep.    

      call reseti (medgeu, (nxp+1)*(nyp+1), 1)
      call reseti (medgev, (nxp+1)*(nyp+1), 1)
#if defined (SUBGRIDEDGE) || ! defined (NOSUBGRIDEDGESHELF)
      do j = 1,ny
        do i = 1,nx
          if (fedge(i,j).lt.1.) then
            do look=1,4
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (.not.(ii.eq.i .and. jj.eq.j)) then      ! for flowline
                if ( .not. ( fedge(ii,jj).eq.1. .and. h(ii,jj).gt.0.
     *                       .and. hs(ii,jj).gt.hs(i,j)
     *                     ) ) then
                  if (look.eq.1) medgeu(i-1,j  ) = 0
                  if (look.eq.2) medgeu(i  ,j  ) = 0
                  if (look.eq.3) medgev(i  ,j-1) = 0
                  if (look.eq.4) medgev(i  ,j  ) = 0
                endif
              endif
            enddo
          endif
        enddo
      enddo
#endif

      nuvtot = 0
      iterasch = 0
      itera = 0
      iterc = 0

c        Impose boundary constraints on h for various expts  

c.......................................
#if defined (EISLINE) && defined (LINEE)
c.......................................
c     ice stream entering from left, prescribed l.h. thickness:
      do j=1,ny
#  if defined (RESTARTE) || defined (SEALEVCHANGE)
ca      h(1,j) = 510.  ! a
        h(1,j) = 550. 
#  else
        if (ifrest.eq.-1) then
ca        h(1,j) = 510. * min (1., timeice/500.)  ! a
          h(1,j) = 550. * min (1., timeice/500.)
        else
ca        h(1,j) = 510. ! a
          h(1,j) = 550.
        endif
#  endif
        h(nx,j) = 0.
c       kill off tiny amounts:
        do i=1,nx
          if (h(i,j).lt.1.e-8) h(i,j) = 0.
        enddo
      enddo
c.......................................................................
#elif ( defined (EISLINE) && ! defined (LINEE) ) || ( defined (EISANTA) && defined (TRANSECTA) )
c.......................................................................
      if (nx.gt.ny) then
        do j=1,ny
          h(nx,j) = 0.
        enddo
      else
        do i=1,nx
          h(i,ny) = 0.
        enddo
      endif

c.......................................................................
#elif ( defined (EISANTA) && ! defined (TRANSECTA) ) || defined (NHA) || defined (CARB) || defined (ORDO)
c.......................................................................

#  if ! defined (NESTING)
      if (nx.gt.1) then
        do j=1,ny
          h(1, j) = 0.
          h(nx,j) = 0.
        enddo
      endif
      if (ny.gt.1) then
        do i=1,nx
          h(i, 1) = 0.
          h(i,ny) = 0.
        enddo
      endif
#  endif

c............................................
#elif defined (TEST2D) || defined (TESTMEL2D)
c............................................
      do j=1,ny
        h(nx,j) = 0.
      enddo
c.....
#endif
c.....

c        Set masks (maskh, masku, maskv). All masks are set 
c        outside iteration loops, ignoring changes of non-zero <-> 
c        zero ice (h) between iterations. 

c        Set maskh (0/1 if no ice/ice, h grid)

      call izero (maskh,  nx*ny)
      do j = 1,ny
        do i = 1,nx
          if (h(i,j).gt.0.) then
            maskh(i,j) = 1
          endif
        enddo
      enddo

c        Set maskc (cross grid)

      call izero (maskc, (nx+1)*(ny+1))
      do j = 1,ny-1
        do i = 1,nx-1
          maskc(i,j)=maskh(i,j)*maskh(i+1,j)*maskh(i,j+1)*maskh(i+1,j+1)
        enddo
      enddo 

c        Set masku,v (1 to allow flow, 0 for no flow if no ice to
c        either side, or lower ice than non-ice surface) (u,v grids).
c        (nb: outside loop C, so ignores hs changes within iters).

      call izero (masku, (nxp+1)*(nyp+1))
      do j = 1,ny
        do i = 1,nx-1
          if ( maskh(i,j).eq.1 .and. maskh(i+1,j).eq.1 ) then
            masku(i,j) = 1
          else if ( maskh(i,j).eq.1 .and. maskh(i+1,j).eq.0 .and. 
     *              hs(i,j).gt.hs(i+1,j) ) then
            masku(i,j) = 1
          else if ( maskh(i,j).eq.0 .and. maskh(i+1,j).eq.1 .and. 
     *              hs(i,j).lt.hs(i+1,j) ) then
            masku(i,j) = 1
          endif
#if defined (NESTING)
c         to leave u along top and bottom edges unchanged (u=h(j) grid)
          if (j.le.1 .or. j.ge.ny) masku(i,j) = 0
#endif
        enddo
      enddo

      call izero (maskv, (nxp+1)*(nyp+1))
      do j = 1,ny-1
        do i = 1,nx
          if ( maskh(i,j).eq.1 .and. maskh(i,j+1).eq.1 ) then
            maskv(i,j) = 1
          else if ( maskh(i,j).eq.1 .and. maskh(i,j+1).eq.0 .and. 
     *              hs(i,j).gt.hs(i,j+1) ) then
            maskv(i,j) = 1
          else if (maskh(i,j).eq.0 .and. maskh(i,j+1).eq.1 .and. 
     *             hs(i,j).lt.hs(i,j+1) ) then
            maskv(i,j) = 1
          endif
#if defined (NESTING)
c         to leave v along left and right edges unchanged (v=h(i) grid)
          if (i.le.1 .or. i.ge.nx) maskv(i,j) = 0
#endif
        enddo
      enddo

c        Zero diagnostic quantities (some necessary for finite diffs 
c        near ice edges)

      call zero (hu,   (nxp+1)*(nyp+1))
      call zero (hv,   (nxp+1)*(nyp+1))
      call zero (dhsu, (nxp+1)*(nyp+1))
      call zero (dhsv, (nxp+1)*(nyp+1))
      call zero (uadv, (nxp+1)*(nyp+1))
      call zero (vadv, (nxp+1)*(nyp+1))

      call zero (lhsu, (nxp+1)*(nyp+1))
      call zero (lhsv, (nxp+1)*(nyp+1))
      call zero (dfu,  (nxp+1)*(nyp+1))
      call zero (dfv,  (nxp+1)*(nyp+1))
      call zero (dbu,  (nxp+1)*(nyp+1))
      call zero (dbv,  (nxp+1)*(nyp+1))

#if ! defined (NESTING)
      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.0) then
            ua(i,j)  = 0.
            uia(i,j) = 0.
            uit(i,j) = 0.
            ub(i,j)  = 0.
            uagrid(i,j) = 0.
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.0) then
            va(i,j)  = 0.
            via(i,j) = 0.
            vit(i,j) = 0.
            vb(i,j)  = 0.
            vagrid(i,j) = 0.
          endif
        enddo
      enddo
#endif

c        If first timestep of run, zero all velocities. Else, keep 
c        for next timestep's initial guess (u,v grids, all passed)
c        (some necessary for finite diffs near ice edges)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (iffirsr.eq.1) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~

c       ua,va,ub,vb are on restart file, so only set here if cold start:
c.............................
        if (ifrest.eq.-1) then
c.............................
#if defined (NESTING)
c         zero only in interior of nested domain, leave edges unchanged:
c         u-grid edges for left and right ua, 
c         h-grid edges for top and bottom ua, 
c         v-grid edges for top and bottom va, 
c         h-grid edges for left and right va. 
          do j=0,nyp
            do i=0,nxp
              if ((i.ge.1 .and.i.le.nx-1).and.(j.gt.1.and.j.lt.ny)) then
                ua(i,j) = 0.
                ub(i,j) = 0.
              else
                ub(i,j) = ua(i,j)
              endif

              if ((j.ge.1 .and.j.le.ny-1).and.(i.gt.1.and.i.lt.nx)) then
                va(i,j) = 0.
                vb(i,j) = 0.
              else
                vb(i,j) = va(i,j)
              endif

            enddo
          enddo
#else
          call zero (ua, (nxp+1)*(nyp+1))
          call zero (va, (nxp+1)*(nyp+1))
          call zero (ub, (nxp+1)*(nyp+1))
          call zero (vb, (nxp+1)*(nyp+1))
#endif
c............
        endif
c............
c       none of following are on restart file (only ua,va,ub,vb are):
        call scopy ((nxp+1)*(nyp+1), ua, 1, uagrid, 1)
        call scopy ((nxp+1)*(nyp+1), va, 1, vagrid, 1)
        call zero (uia,       (nxp+1)*(nyp+1))
        call zero (via,       (nxp+1)*(nyp+1))
        call zero (uit,       (nxp+1)*(nyp+1))
        call zero (vit,       (nxp+1)*(nyp+1))
        call zero (uschoof,   (nxp+1)*(nyp+1))
        call zero (vschoof,   (nxp+1)*(nyp+1))
        call zero (fschoofu,  (nxp+1)*(nyp+1))
        call zero (fschoofv,  (nxp+1)*(nyp+1))
        call resetr (thetau, nx*ny, -1.)
        call resetr (thetav, nx*ny, -1.)
c~~~~~~~~~~
      endif
c~~~~~~~~~~

c        Test for no ice. If so, just apply local forcing 
c        and skip all velocity calcs

      numh = 0
      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            numh = numh + 1
          endif
        enddo
      enddo
      if (numh.eq.0) then
        totflow = 0.
        do j=1,ny
          do i=1,nx
            h(i,j) = h(i,j) + budgall(i,j)*dt
            if (h(i,j).lt.0.) then
              budgneg(i,j) = -h(i,j)/dt 
              h(i,j) = 0.
            else
              budgneg(i,j) = 0.
            endif
            cliffmelt(i,j) = 0.
          enddo
        enddo
        goto 3500
      endif

c============================
c     Top of loop C (~Picard)
      do iterc = 1,niterc
c============================

c        Adjust hw and hs (for change in h):

#if ! defined (NOMOVECOAST)
c     don't call findwater here: causes blowups at 10 km (8/15/13) 
c     call findwater (maskwater, indlake, npoilake, nlake,
c    *                h, hb, hw, sealev, timeice)
#endif
      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

c       Set hu, hv, dhsu, dhsv, dfu, dfv  (u,v grids)
c       Do inside loop C to capture hs changes within iters.

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.1) then
            hu(i,j) = 0.5*(h(i,j) + h(i+1,j))
            dhsu(i,j) = (hs(i+1,j) - hs(i,j)) / dxu(i,j)
            dfu(i,j)  = rhoice*grav*hu(i,j)*dhsu(i,j)
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then
            hv(i,j) = 0.5*(h(i,j) + h(i,j+1))
            dhsv(i,j) = (hs(i,j+1) - hs(i,j)) / dyv(i,j)
            dfv(i,j)  = rhoice*grav*hv(i,j)*dhsv(i,j)
          endif
        enddo
      enddo

c     store original h (for whole time-step) in h0, etc,
c     (and initialize h1):
      if (iterc.eq.1) then 
        call scopy (nx*ny, h, 1, h0, 1)
        call scopy (nx*ny, h, 1, h1, 1)
        call scopy ((nxp+1)*(nyp+1), hu,   1, hu0,   1) !used in icestep
        call scopy ((nxp+1)*(nyp+1), hv,   1, hv0,   1) !used in icestep
        call scopy ((nxp+1)*(nyp+1), dhsu, 1, dhsu0, 1) !used in icestep
        call scopy ((nxp+1)*(nyp+1), dhsv, 1, dhsv0, 1) !used in icestep
        call scopy ((nxp+1)*(nyp+1), dfu,  1, dfu0, 1)  !used in icestep
        call scopy ((nxp+1)*(nyp+1), dfv,  1, dfv0, 1)  !used in icestep
      endif

      call scopy (nx*ny, h1, 1, h1prev, 1)

c        Impose boundary constraints,etc, on velocs for various expts  

c.......................................
#if defined (EISLINE) && defined (LINEE)
c.......................................
c     ice stream entering from left
      do j=1,ny
#  if defined (RESTARTE) || defined (INITEE) || defined (SEALEVCHANGE)
c       ua(0,j) = 0.    ! a
        ua(0,j) = 500. 
#  else
        ua(0,j) = 0.
#  endif
        ub(0,j) = ua(0,j)
      enddo
c.......................................................................
#elif defined EISLINE && ( defined (LINEF) || defined (LINEG) || defined (LINEM) || defined (SLMODEL) )
c.......................................................................
      if (nx.gt.ny) then
        do j=1,ny
          ua(0,j) = 0.
          ub(0,j) = 0.
        enddo
      else
        do i=1,nx
          ua(i,0) = 0.
          ub(i,0) = 0.
        enddo
      endif
c.............................................
#elif defined (EISANTA) && defined (TRANSECTA)
c.............................................
      do j=1,ny
        ua(0,j) = 0.
        ua(nx,j) = 0.
        ub(0,j) = 0.
        ub(nx,j) = 0.
      enddo
c............................................
#elif defined (TEST2D) || defined (TESTMEL2D)
c............................................
      do j=1,ny
        ua(0,j) = 0.
        ub(0,j) = 0.
      enddo
c.....
#endif
c.....

c        Set nuvtot and muind,mvind, indices from 2-D arrays into 
c        linear sequence of selected points for tridia or sparse solns
c        below. Cutoff basal coefficient crhucrit set above. 
c        (Can do outside A loop, depends only on mask[u,v],hw).  

      call izero (muind, (nxp+1)*(nyp+1))
      call izero (mvind, (nxp+1)*(nyp+1))
      nuvtot = 0
c----------------------
#if ! defined (NOSHELF)
c----------------------
      do j = 1,ny
        do i = 1,nx
          if (masku(i,j).eq.1) then
c           include if both adjacent h-points have water, 
c           or if this u-point is slippery enough (crhucrit set above):
            if ( (hw(i,j).gt.hwcut. and. hw(min(i+1,nx),j).gt.hwcut)
     *           .or. (crhu(i,j).gt.crhucrit)
     *         ) then
              nuvtot = nuvtot + 1
              muind(i,j) = nuvtot
            endif
          endif

          if (maskv(i,j).eq.1) then
            if ( (hw(i,j).gt.hwcut. and. hw(i,min(j+1,ny)).gt.hwcut)
     *           .or. (crhv(i,j).gt.crhucrit)
     *         ) then
              nuvtot = nuvtot + 1
              mvind(i,j) = nuvtot
            endif
          endif
        enddo
      enddo
c-----
#endif
c-----

c@@@@@@@@@@@@@@@@@@@@@
#if defined (SCHOOFGL)
c@@@@@@@@@@@@@@@@@@@@@

      call zero (uschoof,   (nxp+1)*(nyp+1))
      call zero (vschoof,   (nxp+1)*(nyp+1))
      call zero (fschoofu,  (nxp+1)*(nyp+1))
      call zero (fschoofv,  (nxp+1)*(nyp+1))
      call resetr (thetau, nx*ny, -1.)
      call resetr (thetav, nx*ny, -1.)

c        Extra elliptic calculation for velocities [u,v]bgrid without
c        any Schoof imposed conditions, used (i) to compare fluxes 
c        at g.l. with Schoof fluxes, for ad-hoc rule in schoofgl,
c        and (ii) to set buttressing factor in thetacalc.
c        n.b.: outside A loop, inside C loop.

      do j=1,ny
        do i=1,nx-1
          if (muind(i,j).eq.0) uagrid(i,j) = ua(i,j)
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          if (mvind(i,j).eq.0) vagrid(i,j) = va(i,j)
        enddo
      enddo

c&&&&&&&&&&&&&&&&&&&&&&&&&&
      if (nuvtot.gt.0) then
c&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     Top of loop ASCH (shelf flow)
      do iterasch = 1,niterasch
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call scopy ((nxp+1)*(nyp+1), uagrid, 1, uaprev, 1)
      call scopy ((nxp+1)*(nyp+1), vagrid, 1, vaprev, 1)

c        Set various coefficients (coef1, coef2, coefb[u,v], and heath)

      call setcoefs (coef1, coef2, visc, coefbu, coefbv, heath,
     *               maskh, arhap, h, hbsd, hw,
     *               masku, maskv, uagrid,vagrid,uia,via, fracgu,fracgv,
     *               uit, vit, sih, siu, siv, hu, hv, crhu, crhv,
     *               muind, mvind, maskc, iterasch, iterc) 

#   if defined (MELANGE)
      call setbpumel (bpumel, bpvmel, thetaumel, thetavmel,
     *                hmel, h, hedge, maskwater, fedge)
#   endif

c        Tridia or Sparse matrix solution of ice-shelf eqns 
c        (for uagrid,vagrid)

c----------------------
      if (ny.eq.1) then
c----------------------
        call dotridia  (uagrid, uschoof, fschoofu, 
     *                  hu, dhsu, hw, coef1, coefbu, uia, uw,
#  if defined (MELANGE)
     *                  bpumel, bpvmel,
#  endif
     *                  muind, nuvtot)
c---------------------------
      else if (nx.eq.1) then
c---------------------------
        call dotridia2 (vagrid, vschoof, fschoofv, 
     *                  hv, dhsv, hw, coef1, coefbv, via, vw,
#  if defined (MELANGE)
     *                  bpumel, bpvmel,
#  endif
     *                  mvind, nuvtot)
c---------
      else
c---------

#if defined (DOSPARSE)

        call dosparse (uagrid,vagrid, uschoof,vschoof,fschoofu,fschoofv,
     *                 hu, hv, dhsu, dhsv,
     *                 coef1, coef2, coefbu, coefbv, uia,via, uw,vw,
     *                 muind, mvind, nuvtot, hw,
#  if defined (MELANGE)
     *                 bpumel, bpvmel,
#  endif
     *                 timeice, iterasch, iterc, 0)

#elif defined (DOSPARSESOR)

        call dosparse_sor 
     *                (uagrid,vagrid, uschoof,vschoof,fschoofu,fschoofv,
     *                 hu, hv, dhsu, dhsv,
     *                 coef1, coef2, coefbu, coefbv, uia,via, uw,vw,
     *                 muind, mvind, maskwater, h, hs, hb, hw,
#  if defined (MELANGE)
     *                 bpumel, bpvmel,
#  endif
     *                 timeice, iterasch, iterc, 0)

#else
        write (ioterm,'(/a)') 
     *    '*** Error (icedyn): must specify -DDOSPARSE or -DDOSPARSESOR'
        stop
#endif

c----------
      endif
c----------

c     speed limit:
      do j=1,ny
        do i=1,nx-1
          uagrid(i,j) = max (-vmax, min (vmax, uagrid(i,j)))
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          vagrid(i,j) = max (-vmax, min (vmax, vagrid(i,j)))
        enddo
      enddo

c        Test for loop ASCH convergence (on velocs uagrid,vagrid)

      delmaxa = 0.
      do j=0,nyp
        do i=0,nxp
          delmaxa = max (delmaxa, abs(uagrid(i,j)-uaprev(i,j)),
     *                            abs(vagrid(i,j)-vaprev(i,j)))
        enddo
      enddo
      if (delmaxa.lt.delcrita) goto 500
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     End loop ASCH (iterasch)
      enddo
      iterasch = niterasch
  500 continue
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&
      endif ! nuvtot > 0 test
c&&&&&&&&&&

c        Sets Schoof-imposed velocities ([u,v]schoof) and weightings
c        for shelf and sheet dynamics (fschoof[u,v]), on u,v grids.
c        Uses [u,v]bgrid from above, and theta[u,v] (buttressing frac)
c        from thetacalc. Also returns various diagnostics:
c        hg[u,v], qg[u,v], at sub-grid g.l. location,
c        and flux* shifted to grid points.
     
      call schoofgl (uschoof, vschoof, fschoofu, fschoofv,
     *               hgu, hgv, qgu, qgv,
     *               fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *               angnormu, angnormv,
     *               h, hs, hb, hw, fedge, maskwater, 
     *               arhap, s2a0, s2b0, siu, siv, lhsu, lhsv, budgall,
     *               crhu, crhv, dfu, dfv, hu, hv, masku, maskv,
     *               muind, mvind, fracgu, fracgv, 
     *               ua, va, uagrid, vagrid, uia, via, thetau, thetav,
#if defined (MELANGE)
     *               thetaumel, thetavmel, hmel,
#endif
     *               sealev, timeice, vmax)

c@@@@@
#endif
c@@@@@

c==========================
      if (nuvtot.gt.0) then ! if no shelf-like pts, skip shelf calcs
c==========================
c===============================
c     Top of loop A (shelf flow)
      do itera = 1,nitera
c===============================

      call scopy ((nxp+1)*(nyp+1), ua, 1, uaprev, 1)
      call scopy ((nxp+1)*(nyp+1), va, 1, vaprev, 1)

c        Set various coefficients (coef1, coef2, coefb[u,v], and heath)

      call setcoefs (coef1, coef2, visc, coefbu, coefbv, heath,
     *               maskh, arhap, h, hbsd, hw,
     *               masku, maskv, ua, va, uia, via, fracgu, fracgv, 
     *               uit, vit, sih, siu, siv, hu, hv, crhu, crhv,
     *               muind, mvind, maskc, itera, iterc) 

c        Tridia or Sparse matrix solution of ice-shelf eqns (for ua,va)

c----------------------
      if (ny.eq.1) then
c----------------------
        call dotridia  (ua, uschoof, fschoofu, 
     *                  hu, dhsu, hw, coef1, coefbu, uia, uw,
#  if defined (MELANGE)
     *                  bpumel, bpvmel,
#  endif
     *                  muind, nuvtot)
c---------------------------
      else if (nx.eq.1) then
c---------------------------
        call dotridia2 (va, vschoof, fschoofv, 
     *                  hv, dhsv, hw, coef1, coefbv, via, vw,
#  if defined (MELANGE)
     *                  bpumel, bpvmel,
#  endif
     *                  mvind, nuvtot)
c---------
      else
c---------

#if defined (DOSPARSE)

        call dosparse (ua,va, uschoof,vschoof, fschoofu,fschoofv,
     *                 hu, hv, dhsu, dhsv,
     *                 coef1, coef2, coefbu, coefbv, uia, via, uw, vw,
     *                 muind, mvind, nuvtot, hw,
#  if defined (MELANGE)
     *                 bpumel, bpvmel,
#  endif
     *                 timeice, itera, iterc, 1)

#elif defined (DOSPARSESOR)

        call dosparse_sor
     *                (ua,va, uschoof,vschoof, fschoofu,fschoofv,
     *                 hu, hv, dhsu, dhsv,
     *                 coef1, coef2, coefbu, coefbv, uia, via, uw, vw,
     *                 muind, mvind, maskwater, h, hs, hb, hw,
#  if defined (MELANGE)
     *                 bpumel, bpvmel,
#  endif
     *                 timeice, itera, iterc, 1)

#else
      write (ioterm,'(/a)') 
     *  '*** Error (icedyn): must specify -DDOSPARSE or -DDOSPARSESOR'
      stop
#endif

c----------
      endif
c----------

c     speed limit:
      do j=1,ny
        do i=1,nx-1
          ua(i,j) = max (-vmax, min (vmax, ua(i,j)))
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          va(i,j) = max (-vmax, min (vmax, va(i,j)))
        enddo
      enddo

c        Test for loop A convergence (on velocs ua,va)

      delmaxa = 0.
      do j=0,nyp
        do i=0,nxp
          delmaxa = max (delmaxa, abs(ua(i,j)-uaprev(i,j)),
     *                            abs(va(i,j)-vaprev(i,j)))
        enddo
      enddo
c     monitor loop a convergence:
c     do iu = ioterm, iuout1d, iuout1d-ioterm
c       write (iu,'(a,i3,2(a,f10.4))')
c    *    'itera=',itera, '  delmaxa=',delmaxa,'  delcrita=',delcrita
c     enddo
      if (delmaxa.lt.delcrita) goto 1000
c=======================
c     End loop A (itera)
      enddo
      itera = nitera
 1000 continue
c=======================
c==========
      endif ! nuvtot > 0 test
c==========

c        Set lhsu,lhsv for icestep. Also set dbu,dbv (diagnostic except
c        used for heatb[u,v] after C loop).

      do j = 1,ny
        do i = 1,nx-1
          if (muind(i,j).ge.1) then
            dbu(i,j) = -coefbu(i,j) * (ua(i,j) - uia(i,j) - uw(i,j))
#if ! defined (NOLHSCROSS)
            lhsu(i,j) = dfu(i,j) - dbu(i,j)
#endif
          else if (masku(i,j).eq.1) then
            dbu(i,j)  = dfu(i,j)
            lhsu(i,j) = 0.
          else
            dbu(i,j)  = 0.
            lhsu(i,j) = 0.
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (mvind(i,j).ge.1) then
            dbv(i,j) = -coefbv(i,j) * (va(i,j) - via(i,j) - vw(i,j))
#if ! defined (NOLHSCROSS)
            lhsv(i,j) = dfv(i,j) - dbv(i,j)
#endif
          else if (maskv(i,j).eq.1) then
            dbv(i,j)  = dfv(i,j)
            lhsv(i,j) = 0.
          else
            dbv(i,j)  = 0.
            lhsv(i,j) = 0.
          endif
        enddo
      enddo

c        Step through one timestep to get h1, with time-implicit 
c        internal flow, and implicit h * explicit ub . 
c        Also return internal flow quantities uia, via, uit, vit 

cm    call mascalc (h0, ztoth0, ztota0)

      call icestep (h0, h, h1, hw, maskwater,
     *              hu0, hv0, dhsu0, dhsv0, dfu0, dfv0,
     *              lhsu, lhsv,
     *              s1a0, s1b0, s2a0, s2b0, siu, siv,
     *              crhu, crhv,
     *              uia, via, uit, vit, ub, vb, ua, va, 
     *              uschoof, vschoof, fschoofu, fschoofv,  
     *              masku, maskv, muind, mvind, medgeu, medgev, dt)

c        Set uadv,vadv fluxes for cliffcalc, and diagnostics in iceshow 
c        (do before h is updated below) 

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.1) then
c           centered for sia uia and sia ub, upstream for ssa ub:
            uadv(i,j) = uia(i,j)*hu(i,j)*dyu(i,j)
            if (muind(i,j).eq.0) then
              uadv(i,j) = uadv(i,j) + ub(i,j)*hu(i,j)*dyu(i,j)
            else
              uadv(i,j) = uadv(i,j) + 
     *          (  max (ub(i,j),0.)*h(i,j) + min (ub(i,j),0.)*h(i+1,j)
     *          ) * dyu(i,j)
            endif
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then
c           centered for sia via and sia vb, upstream for ssa vb:
            vadv(i,j) = via(i,j)*hv(i,j)*dxv(i,j)
            if (mvind(i,j).eq.0) then
              vadv(i,j) = vadv(i,j) + vb(i,j)*hv(i,j)*dxv(i,j)
            else
              vadv(i,j) = vadv(i,j) + 
     *          (  max (vb(i,j),0.)*h(i,j) + min (vb(i,j),0.)*h(i,j+1)
     *          ) * dxv(i,j)
            endif
          endif
        enddo
      enddo

c        Calculate cliffmelt (depends on ua,va at grounding lines)

#if defined (CLIFFMELT)
      call cliffcalc (h, hb, fedge, ifnexo, maskwater, 
     *                cliffmelt, calvice, crevliq,
     *                uadv, vadv, fracgu, fracgv, thetau, thetav,
#  if defined (MELANGE)
     *                ucimel, vcimel, hucimel, hvcimel,
#  endif
     *                sealev, timeice, dt)
#else
      call zero (cliffmelt, nx*ny)
#  if defined (MELANGE)
      call zero (ucimel, (nxp+1)*(nyp+1)) 
      call zero (vcimel, (nxp+1)*(nyp+1)) 
      call zero (hucimel, (nxp+1)*(nyp+1)) 
      call zero (hvcimel, (nxp+1)*(nyp+1)) 
#  endif
#endif

c        Explicitly apply surface mass balance to ice thickness.
c        cliffmelt and calvice are not yet included in budgall
c        (set/modified in cliffcalc), so appear explicitly here.

#if defined (NESTING)
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
          h1(i,j) = h1(i,j) 
     *            + (budgall(i,j) - cliffmelt(i,j) - calvice(i,j))*dt
        enddo
      enddo

cm    call mascalc (h1, ztoth1a, ztota1a)

c        Correct for negative h1

      do j=1,ny
        do i=1,nx
          if (h1(i,j).lt.0.) then
            budgneg(i,j) = -h1(i,j)/dt 
            h1(i,j) = 0.
          else
            budgneg(i,j) = 0.
          endif
        enddo
      enddo

c        Re-impose boundary constraints,etc, on h for various expts  

      totflow = 0.
c....................................................................
#if defined (EISLINE) || ( defined (EISANTA) && defined (TRANSECTA) )
c....................................................................
      if (nx.gt.ny) then
        do j=1,ny
          zh1 = h1(nx,j)
          h1(nx,j) = h0(nx,j)
          totflow = totflow + (h1(nx,j)-zh1)*darea(nx,j)
        enddo
      else
        do i=1,nx
          zh1 = h1(i,ny)
          h1(i,ny) = h0(i,ny)
          totflow = totflow + (h1(i,ny)-zh1)*darea(i,ny)
        enddo
      endif

c.......................................................................
#elif ( defined (EISANTA) && ! defined (TRANSECTA) ) || defined (NHA) || defined (CARB) || defined (ORDO)
c.......................................................................
#  if ! defined (NESTING)
      if (nx.gt.1) then
        do j=1,ny
          totflow = totflow - h1(1,j) *darea(1,j) 
          totflow = totflow - h1(nx,j)*darea(nx,j) 
          h1(1, j) = 0.
          h1(nx,j) = 0.
        enddo
      endif
      if (ny.gt.1) then
        do i=1,nx
          totflow = totflow - h1(i,1) *darea(i,1) 
          totflow = totflow - h1(i,ny)*darea(i,ny) 
          h1(i, 1) = 0.
          h1(i,ny) = 0.
        enddo
      endif
#  elif defined (NESTING)
      if (nx.gt.1) then
        do j=1,ny
          totflow = totflow
     * + ( max (ua(1,j),0.)*h1(1,j) + min (ua(1,j),0.)*h1(2,j)
     *   ) * dyu(1,j) * dt
     * - ( max (ua(nx-1,j),0.)*h1(nx-1,j) + min (ua(nx-1,j),0.)*h1(nx,j)
     *   ) * dyu(nx-1,j) * dt
        enddo
      endif
      if (ny.gt.1) then
        do i=1,nx
           totflow = totflow
     * + ( max (va(i,1),0.)*h1(i,1) + min (va(i,1),0.)*h1(i,2)
     *   ) * dxu(i,1) * dt
     * - ( max (va(i,ny-1),0.)*h1(i,ny-1) + min (va(i,ny-1),0.)*h1(i,ny)
     *   ) * dxu(i,ny-1) * dt
         enddo
       endif
#  endif

c............................................
#elif defined (TEST2D) || defined (TESTMEL2D)
c............................................
      do j=1,ny
        totflow = totflow - h1(nx,j)*darea(nx,j) 
        h1(nx,j) = 0.
      enddo
c.....
#endif
c.....

cm    call mascalc (h1, ztoth1b, ztota1b)
cm    ztotneg = 0.
cm    do j=1,ny
cm      do i=1,nx
cm        ztotneg = ztotneg + budgneg(i,j)*darea(nx,j)*dt
cm      enddo
cm    enddo
cm    write (6,'(a,4f20.14)') 
cm   *    '   toth1b, 1b-1a, totflow, totneg:', 
cm   *     ztoth1b/totarea, (ztoth1b-ztoth1a)/totarea,
cm   *     totflow/totarea, ztotneg/totarea

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Hindmarsh and Payne(1996) convergence:
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     do j=1,ny
c       do i=1,nx
c         delh1(i,j) = h1(i,j) - h1prev(i,j)
c       enddo
c     enddo
c
c     zcos = 999.
c     zthed= 999.
c     zrat = 999.
c     if (iterc.ge.2) then
c       znorma = 0.
c       znormb = 0.
c       znormc = 0.
c       zdot   = 0.
c       do j=1,ny
c         do i=1,nx
c           if (h1(i,j).gt.0.) then
c             znorma = znorma + delh1prev(i,j)**2 
c             znormb = znormb + delh1(i,j)**2
c             znormc = znormc + (delh1(i,j)-delh1prev(i,j))**2
c             zdot   = zdot   + delh1(i,j)*delh1prev(i,j)
c           endif
c         enddo
c       enddo
c       znorma = sqrt(znorma)
c       znormb = sqrt(znormb)
c       znormc = sqrt(znormc)
c       if (znorma.gt.0. .and. znormb.gt.0. .and. znormc.gt.0.) then
c         zcos = zdot/(znorma*znormb)
c         if (abs(zcos).le.1.) then
c           zthe = acos (zcos)
c           zthed= zthe*180./pi
c           if (abs(zthed).lt.30. .or. abs(zthed).gt.180.-30.) then
cc          if (abs(zthed).lt.45. .or. abs(zthed).gt.180.-45.) then
c             zrat = znorma/znormc
c             do j=1,ny
c               do i=1,nx
c                 h1(i,j) = h1prev(i,j) + delh1(i,j)*zrat
c               enddo
c             enddo
c           endif
c         endif
c       endif
c     endif
c
c     call scopy (nx*ny, delh1, 1, delh1prev, 1)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     End of Hindmarsh and Payne convergence
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c        Update h at t+0.5*dt to be used for next iter's setcoefs

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.5*(h1(i,j) + h0(i,j))
        enddo
      enddo

c        Test for loop C convergence (on ice thickness h1)

      delmaxc = 0.
      do j=1,ny
        do i=1,nx
          delmaxc = max (delmaxc, abs(h1(i,j)-h1prev(i,j)))
        enddo
      enddo

c     monitor loop c convergence:
c     do iu = ioterm, iuout1d, iuout1d-ioterm
c       write (iu,'(a,i3,2(a,f10.4))')
c    *    'iterc=',iterc, '  delmaxc=',delmaxc,'  delcritc=',delcritc
cc      for Hindmarsh and Payne convergence:
cc      write (iu,'(9x,3(a,f20.8))') 
cc   *    '  zcos=',zcos, '  zthe=',zthed, '  zrat=',zrat
c     enddo

      if (delmaxc.lt.delcritc) goto 3000

c=======================
c     End loop C (iterc)
      enddo
      iterc = niterc
 3000 continue
c=======================

c        Set new h for this timestep (at t+dt)

      call scopy (nx*ny, h1, 1, h, 1)

c        Calculate basal heating (heatb, J/m2/y) due to basal sliding,
c        equal to taub_u*ub + taub_v*vb. First on u-grid (heatbu) 
c        and v-grid (heatbv), then transfer to h-grid (heatv).
c        heatb is passed to icetherm.
c
c        If defined sediment, heatb has been calculated in sedflow
c        (using last time step quantities), allowing for tau-split
c        between sediment and bare rock.

#if ! defined (SEDIMENT)
c     u-direction:
      call zero (heatbu, (nxp+1)*(nyp+1))
      do j = 1,ny
        do i = 1,nx-1
          if ( masku(i,j).eq.1 
c    *         .and. fracgu(i,j).ne.0.    ! but ignores sub-grid pinning
     *         .and..not. (hw(i,j).gt.hwcut. and. hw(i+1,j).gt.hwcut) 
     *       ) then
            heatbu(i,j) = -(1.-fsedu(i,j)) 
     *                  * dbu(i,j) * (ua(i,j)-uia(i,j)-uw(i,j))
          endif
        enddo
      enddo

c     v-direction:
      call zero (heatbv, (nxp+1)*(nyp+1))
      do j = 1,ny-1
        do i = 1,nx
          if ( maskv(i,j).eq.1
c    *         .and. fracgv(i,j).ne.0.
     *         .and..not. (hw(i,j).gt.hwcut. and. hw(i,j+1).gt.hwcut)
     *       ) then
            heatbv(i,j) = -(1.-fsedv(i,j))
     *                  * dbv(i,j)*(va(i,j)-via(i,j)-vw(i,j))
          endif
        enddo
      enddo

c     transfer to h-grid:
      call zero (heatb, nx*ny)
      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            heatb(i,j) = 0.5*(heatbu(i-1,j) + heatbu(i,j))
     *                 + 0.5*(heatbv(i,j-1) + heatbv(i,j))
c           heatb(i,j) = 0.
          endif
        enddo
      enddo
#endif

c        Skip point for no ice (numh=0, top of loop C)

 3500 continue    

c        Zero very small shelf amounts (h ~1.e-20) created by icedyn

      do j=1,ny
        do i=1,nx
          if (h(i,j).lt.1.e-10) h(i,j) = 0.
        enddo
      enddo

c          Final adjust of maskwater, hw and hs (for change in h):

#if ! defined (NOMOVECOAST)
      call findwater (maskwater, indlake, npoilake, nlake,
     *                h, hb, hw, sealev, timeice)
#endif
      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

      return
      end

c-----------------------------------------------------------------------

      subroutine setcoefs (coef1, coef2, visc, coefbu, coefbv, heath,
     *                     maskh, arhap, h, hbsd, hw,
     *                     masku, maskv, ua,va, uia,via, fracgu, fracgv,
     *                     uit, vit, sih, siu, siv, hu, hv, crhu, crhv,
     *                     muind, mvind, maskc, itera, iterc)

c        Set various coefficients, etc, at top of A loop: 
c        coef[1,2], coefb[u,v], heath

#include <comicephys.h>
#include <comicegrid.h>

c     passed:
      dimension
     *  coef1(nx,ny),
     *  coef2(0:nx,0:ny),                                   ! cross grid
     *  visc(nx,ny), 
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  heath(nx,ny),

     *  maskh(nx,ny),          arhap (nx,ny), 
     *  h(nx,ny),              hbsd(nx,ny),           hw(nx,ny), 
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),      via(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),   fracgv(0:nxp,0:nyp),
     *  uit(0:nxp,0:nyp),      vit(0:nxp,0:nyp),
     *  sih(nx,ny),            siu(0:nxp,0:nyp),  siv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),     crhv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp), 
     *  maskc(0:nx,0:ny)                                    ! cross grid

c     local:
      dimension
     *  g1(nx,ny), g2(nx,ny), g3(nx,ny),
     *  g2c(0:nx,0:ny),                                     ! cross grid
     *  duiz(0:nxp,0:nyp),     dviz(0:nxp,0:nyp)
      save g1, g2, g3, g2c, duiz, dviz

c     Statement function (fpin) for fraction of u,v grid box with 
c     sub-grid bedrock mounds scraping ice bottom:

#include <comicepin.h>

c        Zero initially (some necessary for finite diffs near ice edges)

#if defined (SCHOOFGL)
      if (itera.eq.1) then   ! don't mix between calls for uagrid vs. ua
#else
      if (itera.eq.1 .and. iterc.eq.1) then    ! ok, no calls for uagrid
#endif
        call zero (coef1, nx*ny)
        call zero (coef2,(nx+1)*(ny+1))
        call zero (heath, nx*ny)
        call zero (visc,  nx*ny) 
        call zero (g1,    nx*ny)
        call zero (g2,    nx*ny)
        call zero (g3,    nx*ny)
        call zero (duiz, (nxp+1)*(nyp+1))
        call zero (dviz, (nxp+1)*(nyp+1))
        call zero (g2c,  (nx+1)*(ny+1))
        call zero (sih,   nx*ny)
        call zero (siu,  (nxp+1)*(nyp+1))
        call zero (siv,  (nxp+1)*(nyp+1))
      endif

c        Set g1 ((du/dx)**2 + (dv/dy)**2 + (du/dx)*(dv/dy)) (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            zdudx = (ua(i,j)-ua(i-1,j))/dx(i,j)
            zdvdy = (va(i,j)-va(i,j-1))/dy(i,j)
            g1(i,j) = zdudx**2 + zdvdy**2 + zdudx*zdvdy
          endif
        enddo
      enddo

c        Set g2c (cross grid), then g2 = [du/dy+dv/dx]**2}/4 (h grid)

      do j = 1,ny-1
        do i = 1,nx-1
          if (maskc(i,j).ne.0) then
            zdudy = (ua(i,j+1)-ua(i,j))/dyc(i,j)
            zdvdx = (va(i+1,j)-va(i,j))/dxc(i,j)
            g2c(i,j) = 0.25*((zdudy + zdvdx)**2)
          endif
        enddo
      enddo

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            g2(i,j) = 0.25*(g2c(i-1,j-1)+g2c(i,j-1)+g2c(i-1,j)+g2c(i,j))
          endif
        enddo
      enddo

c       Set internal(d/dz,sheet)-stress contrib to shelf rheology (g3)

c       Set duiz,dviz (u,v grids)

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.1) then
            duiz(i,j) = uit(i,j)/hu(i,j)
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then
            dviz(i,j) = vit(i,j)/hv(i,j)
          endif
        enddo
      enddo

c       Set g3 ([duiz**2 + dviz**2]/4) (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            g3(i,j) = 0.25 * (   (0.5*(duiz(i-1,j)+duiz(i,j)))**2
     *                         + (0.5*(dviz(i,j-1)+dviz(i,j)))**2 )
          endif
        enddo
      enddo

c        Set visc (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            visc(i,j) = 0.5 
     *                / (max(g1(i,j)+g2(i,j)+g3(i,j), gmin) ** powiv)
            heath(i,j) = ((g1(i,j)+g2(i,j)+g3(i,j))**powih) / arhap(i,j)
          else
            heath(i,j) = 0.
          endif
        enddo
      enddo

c        Set coef1 (2*visc*h/(A**1/n)) for ice-shelf eqns (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
c           coef1(i,j) = 2.*visc(i,j)*h(i,j)          / arhap(i,j)
            coef1(i,j) = 2.*visc(i,j)*max(h(i,j),.01) / arhap(i,j)
          endif
        enddo
      enddo

c        Set coef2 (visc*h/(A**1/n)) for ice-shelf eqns (cross grid)

      do j = 1,ny-1
        do i = 1,nx-1
          if (maskc(i,j).ne.0) then
            coef2(i,j) = 0.5 * 0.25 * (   coef1(i,j)   + coef1(i+1,j)
     *                                  + coef1(i,j+1) + coef1(i+1,j+1))
          endif
        enddo
      enddo

c        Set linearized basal sliding coeffs coefb[u,v] for shelf eqns 
c        = (crhu**(-1/m)) (ub2**((1-m)/2m)) (u,v grids)

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.1.) then

c           calculate zfpin, fraction of u-grid cell actually in 
c           contact with bottom, scraped by sub-grid bedrock mounds 
c           (stmt fnct fpin):
            zhw  = 0.5*(hw(i,j) + hw(i+1,j))
            zhbsd  = 0.5*(hbsd(i,j) + hbsd(i+1,j))
            zfpin = fpin(zhw,zhbsd)

c           calculate coefbu (N/m2 per m/y), weighting by water vs.
c           non-water contact fractions: 
            if (fracgu(i,j).eq.0. .and. zfpin.eq.0.) then
              coefbu(i,j) = coefbwater
            else
              zc = crhu(i,j) ** (-1./powb)
              zpow = -0.5*(1. - 1./powb)
              zu = ua(i,j) - uia(i,j)
              zv = 0.25*( va(i,j) + va(i+1,j) + va(i,j-1) + va(i+1,j-1)
     *                   -via(i,j)-via(i+1,j) -via(i,j-1) -via(i+1,j-1))
              zcoef = zc * ( max(zu**2+zv**2, ubmin**2)**zpow )
              zcoef = max (zcoef, coefbwater)
              if (fracgu(i,j).eq.0.) then
                coefbu(i,j) = (1.-zfpin) * coefbwater
     *                      +     zfpin  * zcoef
              else
                coefbu(i,j) = zcoef
              endif
            endif
          else
            coefbu(i,j) = 0. ! not used
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then

            zhw   = 0.5*(hw(i,j) + hw(i,j+1))
            zhbsd  = 0.5*(hbsd(i,j) + hbsd(i,j+1))
            zfpin = fpin(zhw,zhbsd)

            if (fracgv(i,j).eq.0. .and. zfpin.eq.0.) then
              coefbv(i,j) = coefbwater
            else
              zc = crhv(i,j) ** (-1./powb)
              zpow = -0.5*(1. - 1./powb)
              zu = 0.25*( ua(i,j) + ua(i,j+1) + ua(i-1,j) + ua(i-1,j+1)
     *                   -uia(i,j)-uia(i,j+1) -uia(i-1,j) -uia(i-1,j+1))
              zv = va(i,j) - via(i,j)
              zcoef = zc * ( max(zu**2+zv**2, ubmin**2)**zpow )
              zcoef = max (zcoef, coefbwater)
              if (fracgv(i,j).eq.0.) then
                coefbv(i,j) = (1.-zfpin) * coefbwater
     *                      +     zfpin  * zcoef
              else
                coefbv(i,j) = zcoef
              endif
            endif
          else
            coefbv(i,j) = 0. ! not used
          endif
        enddo
      enddo

c        Set sih on h grid, then siu, siv on u,v grids (shelf-stress 
c        contrib to sheet rheol, used in icestep,icetherm).
c        Don't need to do within ASCH,A loops, only after A loop,
c        but "belongs" in setcoefs (and uses visc, g1, g2).

#if ! defined (NOSOFTCROSS)
      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            sih(i,j)= ((2.*visc(i,j)/arhap(i,j))**2) * (g1(i,j)+g2(i,j))
          endif
        enddo
      enddo

      do j = 1,ny
        do i = 1,nx-1
          if (muind(i,j).ge.1) then
            siu(i,j)  = 0.5 * (sih(i,j) + sih(i+1,j))
          else
            siu(i,j) = 0.
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (mvind(i,j).ge.1) then
            siv(i,j)  = 0.5 * (sih(i,j) + sih(i,j+1))
          else
            siv(i,j) = 0.
          endif
        enddo
      enddo
#endif

      return
      end

c-----------------------------------------------------------------------

      subroutine dotridia (ua, uschoof, fschoofu, 
     *                     hu, dhsu, hw, coef1, coefbu, uia, uw,
#  if defined (MELANGE)
     *                     bpumel, bpvmel,
#  endif
     *                     muind,nuvtot)

c     Assembles *tri arrays for 1-D tridiagonal solution.
c     Calls tridia_i, sets new velocities ua.
c     Only for 1-D in x with ny=1 (use dotridia2 for 1-D in y, nx=1).
c     Weights imposed Schoof velocities (uschoof) into model eqation,
c     with weights fschoofu...which are only potentially non-zero 
c     for 2 u-grid points closest and downstream of grounding line.
c     Do this weighting as follows:
c        (1-fschoofu)*(orig tridiag eqn) + fschoofu*atri*(u=uschoof),
c     i.e., multiply the imposed equation by atri from the orig eqn,
c     so they're both the same order of magn. for intermediate fschoofs.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  ua(0:nxp,0:nyp),
     *  uschoof(0:nxp,0:nyp),  fschoofu(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       dhsu(0:nxp,0:nyp),     
     *  hw(nx,ny),             coef1(nx,ny),        coefbu(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),      uw(0:nxp,0:nyp),     muind(0:nxp,0:nyp)
#if defined (MELANGE)
     * ,bpumel(0:nxp,0:nyp),   bpvmel(0:nxp,0:nyp)
#endif

      dimension atri(nx), btri(nx), ctri(nx), rtri(nx), vtri(nx)

c        Tridiagonal solution for ua in 1-D problems (in x, with ny=1).
c        Relies on muind running left to right (gaps in ice are ok). 

      j = 1

      call zero (atri, nx)
      call zero (btri, nx)
      call zero (ctri, nx)
      call zero (rtri, nx)
      call zero (vtri, nx)

      do i=1,nx-1
        n1 = muind(i,j)

        if (n1.ge.1) then

          rtri(n1) = rhoice*grav*hu(i,j)*dhsu(i,j) 
#if defined (MELANGE)
     *               * bpumel(i,j)
#endif
     *             - coefbu(i,j)*(uia(i,j)+uw(i,j))

          z1 = 2.*coef1(i+1,j) / (dx(i+1,j)*dxu(i,j))  ! ua(i+1,j)
          z2 = 2.*coef1(i  ,j) / (dx(i,j  )*dxu(i,j))  ! ua(i-1,j)

c         diagonal:
          atri(n1) = - (z1 + z2) - coefbu(i,j)
#if defined (SIDEDRAG)
          if (hw(i,j).gt.hwcut. and. hw(i+1,j).gt.hwcut)
     *      atri(n1) = atri(n1) - SIDEDRAG*hu(i,j)
#endif

c         superdiagonal:
          if (muind(i+1,j).ge.1) then
            btri(n1) = z1
          else
            rtri(n1) = rtri(n1) - z1*ua(i+1,j)
          endif

c         subdiagonal:
          if (muind(i-1,j).ge.1) then
            ctri(n1) = z2
          else
            rtri(n1) = rtri(n1) - z2*ua(i-1,j)
          endif

c         impose Schoof velocity: 
          zwei = fschoofu(i,j)
c         atri(n1) = (1.-zwei)*atri(n1) + zwei*atri(n1) ! atri unchanged
          btri(n1) = (1.-zwei)*btri(n1)
          ctri(n1) = (1.-zwei)*ctri(n1)
          rtri(n1) = (1.-zwei)*rtri(n1) + zwei*atri(n1)*uschoof(i,j)

        endif

      enddo

      call tridia_i (atri, btri, ctri, vtri, rtri, nuvtot, 0)

c        Set new ua

      do i=1,nx
        n1 = muind(i,j)
        if (n1.ge.1) then
          ua(i,j) = vtri(n1)
        endif
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine dotridia2 (va, vschoof, fschoofv,
     *                      hv, dhsv, hw, coef1, coefbv, via, vw,
#  if defined (MELANGE)
     *                      bpumel, bpvmel,
#  endif
     *                      mvind, nuvtot)

c     Assembles *tri arrays for 1-D tridiagonal solution.
c     Calls tridia_i, sets new velocities va.
c     Only for 1-D in y with nx=1. (Use dotridia for 1-D in x, ny=1).
c     See comments in dotridia.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  va(0:nxp,0:nyp),
     *  vschoof(0:nxp,0:nyp),  fschoofv(0:nxp,0:nyp),
     *  hv(0:nxp,0:nyp),       dhsv(0:nxp,0:nyp),
     *  hw(nx,ny),             coef1(nx,ny),        coefbv(0:nxp,0:nyp),
     *  via(0:nxp,0:nyp),      vw(0:nxp,0:nyp),     mvind(0:nxp,0:nyp)
#  if defined (MELANGE)
     * ,bpumel(0:nxp,0:nyp),   bpvmel(0:nxp,0:nyp)
#  endif

      dimension atri(ny), btri(ny), ctri(ny), rtri(ny), vtri(ny)

c        Tridiagonal solution for va in 1-D problems (in y, with nx=1).
c        Relies on mvind running "bottom to top" (gaps in ice are ok). 

      i = 1

      call zero (atri, ny)
      call zero (btri, ny)
      call zero (ctri, ny)
      call zero (rtri, ny)
      call zero (vtri, ny)

      do j=1,ny-1
        n1 = mvind(i,j)

        if (n1.ge.1) then

          rtri(n1) = rhoice*grav*hv(i,j)*dhsv(i,j) 
#if defined (MELANGE)
     *               * bpvmel(i,j)
#endif
     *             - coefbv(i,j)*(via(i,j)+vw(i,j))

          z1 = 2.*coef1(i,j+1) / (dy(i,j+1)*dyv(i,j))  ! va(i,j+1)
          z2 = 2.*coef1(i,j  ) / (dy(i,j  )*dyv(i,j))  ! va(i,j-1)

c         diagonal:
          atri(n1) = - (z1 + z2) - coefbv(i,j)
#if defined (SIDEDRAG)
          if (hw(i,j).gt.hwcut. and. hw(i,j+1).gt.hwcut)
     *      atri(n1) = atri(n1) - SIDEDRAG*hv(i,j)
#endif

c         superdiagonal:
          if (mvind(i,j+1).ge.1) then
            btri(n1) = z1
          else
            rtri(n1) = rtri(n1) - z1*va(i,j+1)
          endif

c         subdiagonal:
          if (mvind(i,j-1).ge.1) then
            ctri(n1) = z2
          else
            rtri(n1) = rtri(n1) - z2*va(i,j-1)
          endif

c         impose Schoof velocity:
          zwei = fschoofv(i,j)
c         atri(n1) = (1.-zwei)*atri(n1) + zwei*atri(n1) ! atri unchanged
          btri(n1) = (1.-zwei)*btri(n1)
          ctri(n1) = (1.-zwei)*ctri(n1)
          rtri(n1) = (1.-zwei)*rtri(n1) + zwei*atri(n1)*vschoof(i,j)

        endif

      enddo

      call tridia_i (atri, btri, ctri, vtri, rtri, nuvtot, 0)

c        Set new va

      do j=1,ny
        n1 = mvind(i,j)
        if (n1.ge.1) then
          va(i,j) = vtri(n1)
        endif
      enddo

      return
      end

c-----------------------------------------------------------------------

c|||||||||||||||||||||
#if defined (DOSPARSE)
c|||||||||||||||||||||

      subroutine dosparse (ua, va, uschoof, vschoof, fschoofu, fschoofv,
     *                     hu, hv, dhsu, dhsv,
     *                     coef1, coef2, coefbu, coefbv, uia,via, uw,vw,
     *                     muind, mvind, nuvtot, hw,
#  if defined (MELANGE)
     *                     bpumel, bpvmel,
#  endif
     *                     timeice, itera, iterc, iflag)

c     Solves elliptic shelf equations for ua,va, 
c     just for points with m[u,v]ind >= 1. 
c     Solves using choice of sparse solvers (or Gaussian elimination). 
c     Assembles elspa, ijspa, rhs arrays for sparse solver,
c     calls sparse solver, sets new velocities ua,va
c     Weights imposed Schoof velocities (u,vschoof) into model eqation,
c     with weights fschoofu,v - only potentially non-zero for 
c     2 u,v-grid points closest and downstream of grounding line.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

      dimension 
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  uschoof(0:nxp,0:nyp),  vschoof(0:nxp,0:nyp),
     *  fschoofu(0:nxp,0:nyp), fschoofv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  dhsu(0:nxp,0:nyp),     dhsv(0:nxp,0:nyp),
     *  coef1(nx,ny),          coef2(0:nx,0:ny),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),      via(0:nxp,0:nyp),
     *  uw(0:nxp,0:nyp),       vw(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  hw(nx,ny)
#if defined (MELANGE)
     * ,bpumel(0:nxp,0:nyp),   bpvmel(0:nxp,0:nyp)
#endif

      dimension vecu(nuvmax)
      save vecu

      logical firstvecu
      data firstvecu /.true./
      save firstvecu

c     elspa, ijspa, in comicesparse.h, shared with solver.
c     vecu saved for subsequent calls to linbcg.
      dimension rhs(nuvmax)

c        Set sparse arrays and rhs's needed to solve system of eqns.
c        (using elspa,ijspa,rhs). Can't parallelize for 
c        sparse method, since off-diag terms in ijspa,elspa are 
c        assembled using incremented index kof. Must loop through 
c        i,j,muind,mvind in same sequence as muind,mvind were assembled
c        above, so that diagonal elements are in same order as
c        off-diagonal column groups in elspa,ijspa.
c
c        If the neighbor's u or v are not being computed (n2 = 0), add
c        their (explicit) contributions to rhs(n1).

      call izero (ijspa, nspamax)
      call  zero (elspa, nspamax)
      call  zero (rhs,   nuvtot)

      kof = nuvtot + 1

      do j=1,ny
        do i=1,nx

c         u velocities:

          n1 = muind(i,j)
c<<<<<<<<<<<<<<<<<<<<<<<<<<
          if (n1.ge.1) then
c<<<<<<<<<<<<<<<<<<<<<<<<<<
            kofprev = kof

c              rhs for ua(i,j):

              rhs(n1) = rhoice*grav*hu(i,j)*dhsu(i,j)
#if defined (MELANGE)
     *                  * bpumel(i,j)
#endif
     *                - coefbu(i,j)*(uia(i,j)+uw(i,j))

              z1 = 2.*coef1(i+1,j) / (dx(i+1,j )*dxu(i,j))  ! ua(i+1,j)
              z2 = 2.*coef1(i  ,j) / (dx(i,j   )*dxu(i,j))  ! ua(i-1,j)
              z3 =    coef2(i  ,j) / (dyc(i,j  )*dyu(i,j))  ! ua(i,j+1)
              z4 =    coef2(i,j-1) / (dyc(i,j-1)*dyu(i,j))  ! ua(i,j-1)

c                diagonal for ua(i,j):

              elspa(n1) = elspa(n1) - (z1 + z2 + z3 + z4) - coefbu(i,j)

              ijspa(n1) = kofprev + 1 

c                u-neighbors for ua(i,j):

              n2 = muind(i+1,j)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + z1
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - z1*ua(i+1,j)
              endif

              n2 = muind(i-1,j)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + z2
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - z2*ua(i-1,j)
              endif
  
              n2 = muind(i,j+1)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + z3
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - z3*ua(i,j+1)
              endif
  
              n2 = muind(i,j-1)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + z4
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - z4*ua(i,j-1)
              endif
  
c                v-neighbors for ua(i,j):
    
              za = coef1(i+1,j) / (dy(i+1,j )*dxu(i,j)) ! i+1,j; i+1,j-1
              zb = coef1(i  ,j) / (dy(i,j   )*dxu(i,j)) ! i,j  ; i,j-1
              zc = coef2(i  ,j) / (dxc(i,j  )*dyu(i,j)) ! i+1,j; i,j
              zd = coef2(i,j-1) / (dxc(i,j-1)*dyu(i,j)) ! i,j-1; i+1,j-1
  
              n2 = mvind(i+1,j)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + za + zc
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - (za+zc)*va(i+1,j)
              endif
  
              n2 = mvind(i+1,j-1)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) - za - zd
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) + (za+zd)*va(i+1,j-1)
              endif
  
              n2 = mvind(i,j)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) - zb - zc
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) + (zb+zc)*va(i,j)
              endif
  
              n2 = mvind(i,j-1)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + zb + zd
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - (zb+zd)*va(i,j-1)
              endif

c           endif  ! uag != 0

            if (fschoofu(i,j).gt.0.) then
              zwei = fschoofu(i,j)
              elspa(n1) = (1.-zwei)*elspa(n1) + zwei*elspa(n1)
              rhs(n1) = (1.-zwei)*rhs(n1) + zwei*elspa(n1)*uschoof(i,j)
              if (kof.ge.kofprev+1) then
                do kofsch = kofprev+1, kof
                  elspa(kofsch) = (1.-zwei)*elspa(kofsch)
                enddo
              endif
            endif
               
c<<<<<<<<<<<<<<
          endif    ! n1 >= 1
c<<<<<<<<<<<<<<

c         v velocities:

          n1 = mvind(i,j)
c^^^^^^^^^^^^^^^^^^^^^^^^^^
          if (n1.ge.1) then
c^^^^^^^^^^^^^^^^^^^^^^^^^^
            kofprev = kof

c                rhs for va(i,j):

              rhs(n1) = rhoice*grav*hv(i,j)*dhsv(i,j)
#if defined (MELANGE)
     *                  * bpvmel(i,j)
#endif
     *                - coefbv(i,j)*(via(i,j)+vw(i,j))

              z1 = 2.*coef1(i,j+1) / (dy(i,j+1 )*dyv(i,j))  ! va(i,j+1)
              z2 = 2.*coef1(i  ,j) / (dy(i,j   )*dyv(i,j))  ! va(i,j-1)
              z3 =    coef2(i  ,j) / (dxc(i,j  )*dxv(i,j))  ! va(i+1,j)
              z4 =    coef2(i-1,j) / (dxc(i-1,j)*dxv(i,j))  ! va(i-1,j)
 
c                diagonal for va(i,j):

              elspa(n1) = elspa(n1) - (z1 + z2 + z3 + z4) - coefbv(i,j)

              ijspa(n1) = kofprev + 1 

c                v-neighbors for va(i,j):

              n2 = mvind(i,j+1)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + z1
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - z1*va(i,j+1)
              endif
  
              n2 = mvind(i,j-1)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + z2
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - z2*va(i,j-1)
              endif
  
              n2 = mvind(i+1,j)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + z3
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - z3*va(i+1,j)
              endif
  
              n2 = mvind(i-1,j)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + z4
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - z4*va(i-1,j)
              endif
  
c                u-neighbors for va(i,j):
  
              za = coef1(i,j+1) / (dx(i,j+1 )*dyv(i,j)) ! i,j+1; i-1,j+1
              zb = coef1(i  ,j) / (dx(i,j   )*dyv(i,j)) ! i,j  ; i-1,j
              zc = coef2(i  ,j) / (dyc(i,j  )*dxv(i,j)) ! i,j+1; i,j
              zd = coef2(i-1,j) / (dyc(i-1,j)*dxv(i,j)) ! i-1,j; i-1,j+1
  
              n2 = muind(i,j+1)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + za + zc
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - (za+zc)*ua(i,j+1)
              endif
  
              n2 = muind(i-1,j+1)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) - za - zd
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) + (za+zd)*ua(i-1,j+1)
              endif
  
              n2 = muind(i,j)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) - zb - zc
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) + (zb+zc)*ua(i,j)
              endif
  
              n2 = muind(i-1,j)
              if (n2.ge.1) then
                kof = kof + 1
                elspa(kof) = elspa(kof) + zb + zd
                ijspa(kof) = n2
              else
                rhs(n1) = rhs(n1) - (zb+zd)*ua(i-1,j)
              endif

c           endif ! vag != 0

            if (fschoofv(i,j).gt.0.) then
              zwei = fschoofv(i,j)
              elspa(n1) = (1.-zwei)*elspa(n1) + zwei*elspa(n1)
              rhs(n1) = (1.-zwei)*rhs(n1) + zwei*elspa(n1)*vschoof(i,j)
              if (kof.ge.kofprev+1) then
                do kofsch = kofprev+1, kof
                  elspa(kofsch) = (1.-zwei)*elspa(kofsch)
                enddo
              endif
            endif

c^^^^^^^^^^^^^^
          endif   ! n1 >= 1
c^^^^^^^^^^^^^^

        enddo
      enddo

      if (kof.gt.nspamax) then
        write (6,*) 'Error: kof > nspamax. kof=',kof,' nspamax=',nspamax
        stop
      endif 

      ijspa(nuvtot+1) = kof + 1

      if (nuvtot.ne.ijspa(1) - 2) then
        write (6,*) 'Error: nuvtot, ijspa(1)-2 =',nuvtot,ijspa(1)-2
        stop
      endif

c        Solve linear system arr*vec = rhs for vec.
c        If number of pts <= nuvsmall, use gaussian elimination 
c        (avoids sparse-solution problems for nuvtot <= ~2)

c+++++++++++++++++++++++++++++++++
      if (nuvtot.gt.nuvsmall) then
c+++++++++++++++++++++++++++++++++

#if defined (NUMREC) || defined (SLAP)
      if (firstvecu) then
        call zero (vecu, nuvtot)
        firstvecu = .false.
      endif
#endif
#if defined (NUMREC)
      itol = 2 ! 4 ! 1
      tol = 1.e-4  ! 1.e-4 ! 1.e-6 ! 1.e-11
      itmax = 10000
      call linbcg (nuvtot, rhs, vecu, itol, tol, itmax, itlin, errlin)
c     write (ioterm,'(a,i5,a,e15.5)')
c    *  'linbcg(u): itlin=',itlin,' errlin=',errlin
      if (abs(errlin).gt.tol) write (ioterm,'(a,i5,a,e15.5)')
     *  '*** Warning linbcg(u): itlin=',itlin,' errlin=',errlin
#elif defined (SLAP)
      itol = 2
      tol = 1.e-4
      itmax = 10000
      call slap (nuvtot, rhs, vecu, itol, tol, itmax, itlin, errlin)
c     write (ioterm,'(a,i5,a,e15.5)')
c    *  'slap(u): itlin=',itlin,' errlin=',errlin
      if (abs(errlin).gt.tol) write (ioterm,'(a,i5,a,e15.5)')
     *  '*** Warning slap(u): itlin=',itlin,' errlin=',errlin
#elif defined (SUPERLU)
      call superlu (nuvtot, rhs, vecu, info)
      if (info.ne.0)write (ioterm,*)'*** Warning superlu(u): info=',info
#elif defined (WATSON)
      call watson (nuvtot, rhs, vecu, info)
      if (info.ne.0)write (ioterm,*) '*** Warning watson(u): info=',info
#elif defined (MKL)
      call mkl (nuvtot, rhs, vecu, info)
      if (info.ne.0) write (ioterm,*)'*** Warning mkl(u): info=',info
#elif defined (GAUSS)
      call gaussdo (nuvtot, rhs, vecu, info)
      if (info.ne.0) write (ioterm,*)'*** Warning gauss(u): info=',info
#endif

c+++++++++
      else
c+++++++++
        call gaussdo (nuvtot, rhs, vecu, info)
        if (info.ne.0) write (ioterm,*)'*** Warning gauss(u):info=',info
c++++++++++
      endif
c++++++++++

c        Set new ua,va

      do j=1,ny
        do i=1,nx
          n1 = muind(i,j)
          if (n1.ge.1) then
            ua(i,j) = vecu(n1)
          endif

          n2 = mvind(i,j)
          if (n2.ge.1) then
            va(i,j) = vecu(n2)
          endif
        enddo
      enddo

      return
      end

c|||||
#endif
c|||||

c-----------------------------------------------------------------------

c========================
#if defined (DOSPARSESOR)
c========================

      subroutine dosparse_sor
     *                    (ua, va, uschoof, vschoof, fschoofu, fschoofv,
     *                     hu, hv, dhsu, dhsv,
     *                     coef1, coef2, coefbu, coefbv, uia,via, uw,vw,
     *                     muind, mvind, maskwater, h, hs, hb, hw,
#  if defined (MELANGE)
     *                     bpumel, bpvmel,
#  endif
     *                     timeice, itera, iterc, iflag)

c     Solevs elliptic shelf equations for ua,va, using SOR (Successive 
c     Over Relaxation), just for points with m[u,v]ind >= 1. 
c     Weights imposed Schoof velocities (u,vschoof) into model eqation,
c     with weights fschoofu,v - only potentially non-zero for 
c     2 u,v-grid points closest and downstream of grounding line.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

      dimension 
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  uschoof(0:nxp,0:nyp),  vschoof(0:nxp,0:nyp),
     *  fschoofu(0:nxp,0:nyp), fschoofv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  dhsu(0:nxp,0:nyp),     dhsv(0:nxp,0:nyp),
     *  coef1(nx,ny),          coef2(0:nx,0:ny),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),      via(0:nxp,0:nyp),
     *  uw(0:nxp,0:nyp),       vw(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  maskwater(nx,ny),
     *  h(nx,ny),              hs(nx,ny),
     *  hb(nx,ny),             hw(nx,ny)
#if defined (MELANGE)
     * ,bpumel(0:nxp,0:nyp),   bpvmel(0:nxp,0:nyp)
#endif

#if defined (TEST2D)
      parameter (nitersurmax =   30, nitersurmin = 10, facsor=2.0)
      parameter (ducrit      = .003)                
#else
c     parameter (nitersurmax =  500, nitersurmin = 20, facsor=0.9)
c     parameter (nitersurmax =  100, nitersurmin = 20, facsor=1.2) ! std
      parameter (nitersurmax =   20, nitersurmin = 10, facsor=1.4) ! new
      parameter (ducrit      = .050)                               ! new
#endif

c     for diagnostic dump only:
c     dimension 
c    *  uaorig(0:nxp,0:nyp),       vaorig(0:nxp,0:nyp)

#if ! defined (NOCHECKERBOARD)
c     for checkerboard sor:
      logical firstcheck
      data firstcheck /.true./
      dimension mcheck(0:nxp, 0:nyp)
      save firstcheck, mcheck

      if (firstcheck) then
        do j=0,ny+1
          do i=0,nx+1
            mcheck(i,j) = mod (i+j,2)
          enddo
        enddo
        firstcheck = .false.
      endif
#endif

c     call scopy ((nxp+1)*(nyp+1), ua, 1, uaorig, 1)
c     call scopy ((nxp+1)*(nyp+1), va, 1, vaorig, 1)

c        Calculate mean veloc (was used to test iter skip, now diag)

      zuav = 0.
      ndav = 0.
      do j=1,ny
        do i=1,nx
          if (muind(i,j).ge.1) then
            zuav = zuav + ua(i,j)**2
            ndav = ndav + 1
          endif
          if (mvind(i,j).ge.1) then
            zuav = zuav + va(i,j)**2
            ndav = ndav + 1
          endif
        enddo
      enddo
      zuav = sqrt (zuav/ max(ndav,1))

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do itersur = 1,nitersurmax
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      zdumax =  0.

#if ! defined (NOCHECKERBOARD)
      do icheck = 0,1    ! for checkerboard sor updating
      do j=1,ny
        do i=1,nx-1
          if (muind(i,j).ge.1 .and. mcheck(i,j).eq.icheck) then
#else
      do j=1,ny
        do i=1,nx
          if (muind(i,j).ge.1) then
#endif

            ip1 = min (i+1,nx)

c              u-neighbors for ua(i,j):

            z1 = 2.*coef1(ip1,j) / (dx(ip1,j )*dxu(i,j))  ! i+1,j; i,j
            z2 = 2.*coef1(i  ,j) / (dx(i,j   )*dxu(i,j))  ! i,j  ; i-1,j
            z3 =    coef2(i  ,j) / (dyc(i,j  )*dyu(i,j))  ! i,j+1; i,j
            z4 =    coef2(i,j-1) / (dyc(i,j-1)*dyu(i,j))  ! i,j  ; i,j-1

c              v-neighbors for ua(i,j):
    
            za = coef1(ip1,j) / (dy(ip1,j )*dxu(i,j)) ! i+1,j; i+1,j-1
            zb = coef1(i  ,j) / (dy(i,j   )*dxu(i,j)) ! i,j  ; i,j-1
            zc = coef2(i  ,j) / (dxc(i,j  )*dyu(i,j)) ! i+1,j; i,j
            zd = coef2(i,j-1) / (dxc(i,j-1)*dyu(i,j)) ! i,j-1; i+1,j-1

            zforcu = 
     *          z1*(ua(i+1,j)        )     - z2*(       -ua(i-1,j)) 
     *        + z3*(ua(i,j+1)        )     - z4*(       -ua(i,j-1)) 

     *        + za*(va(i+1,j)-va(i+1,j-1)) - zb*(va(i,j)-va(i,j-1))
     *        + zc*(va(i+1,j)-va(i,j))     - zd*(va(i+1,j-1)-va(i,j-1))

     *        - rhoice*grav*hu(i,j)*dhsu(i,j) 
#if defined (MELANGE)
     *          * bpumel(i,j)
#endif
     *        + coefbu(i,j)*(uia(i,j)+uw(i,j))

            zclhsu = z1 + z2 + z3 + z4 + coefbu(i,j)

c@@@ std:
c           zdu = zforcu/zclhsu - ua(i,j)
c           ua(i,j) = facsor*zforcu/zclhsu + (1.-facsor)*ua(i,j)
c           impose (weighted) Schoof velocities:
c           zwei = fschoofu(i,j)
c           ua(i,j) = (1.-zwei)*ua(i,j) + zwei*uschoof(i,j)
c           if (zwei.eq.0.) zdumax = max (zdumax, abs(zdu))
c@@@ new:
c           Set new ua (zforcu/zclhsu or uschoof, weighted by fschoofu).
c           Modify by facsor (std sor if fschoofu=0), but facsor -> 1 
c           as fschoofu ->1, so ua=uschoof exactly for fschoofu=1.
            zs = fschoofu(i,j)
            zf = (1.-zs)*facsor + zs*1.
            uaold = ua(i,j)
            ua(i,j) = zf * ((1.-zs)*(zforcu/zclhsu) + zs*uschoof(i,j))
     *              + (1.-zf) * ua(i,j)
c           for test to skip out of iteration:
            zdu = abs(ua(i,j)-uaold) / max (abs(uaold),1.)
            zdumax = max (zdumax, abs(zdu))
c@@@

          endif
#if ! defined (NOCHECKERBOARD)
        enddo
      enddo
      enddo

      do icheck = 0,1    ! for checkerboard sor updating
      do j=1,ny-1
        do i=1,nx
          if (mvind(i,j).ge.1 .and. mcheck(i,j).eq.icheck) then
#else
          if (mvind(i,j).ge.1) then
#endif

            jp1 = min (j+1,ny)

c              v-neighbors for va(i,j):

            z1 = 2.*coef1(i,jp1) / (dy(i,jp1 )*dyv(i,j))  ! va(i,j+1)
            z2 = 2.*coef1(i  ,j) / (dy(i,j   )*dyv(i,j))  ! va(i,j-1)
            z3 =    coef2(i  ,j) / (dxc(i,j  )*dxv(i,j))  ! va(i+1,j)
            z4 =    coef2(i-1,j) / (dxc(i-1,j)*dxv(i,j))  ! va(i-1,j)

c              u-neighbors for va(i,j):
 
            za = coef1(i,jp1) / (dx(i,jp1 )*dyv(i,j)) ! i,j+1; i-1,j+1
            zb = coef1(i  ,j) / (dx(i,j   )*dyv(i,j)) ! i,j  ; i-1,j
            zc = coef2(i  ,j) / (dyc(i,j  )*dxv(i,j)) ! i,j+1; i,j
            zd = coef2(i-1,j) / (dyc(i-1,j)*dxv(i,j)) ! i-1,j; i-1,j+1

            zforcv = 
     *          z1*(va(i,j+1)        )     - z2*(       -va(i,j-1)) 
     *        + z3*(va(i+1,j)        )     - z4*(       -va(i-1,j)) 

     *        + za*(ua(i,j+1)-ua(i-1,j+1)) - zb*(ua(i,j)-ua(i-1,j))
     *        + zc*(ua(i,j+1)-ua(i,j))     - zd*(ua(i-1,j+1)-ua(i-1,j))

     *        - rhoice*grav*hv(i,j)*dhsv(i,j) 
#if defined (MELANGE)
     *          * bpvmel(i,j)
#endif
     *        + coefbv(i,j)*(via(i,j)+vw(i,j))

            zclhsv = z1 + z2 + z3 + z4 + coefbv(i,j)

c@@@ std:
c           zdv = zforcv/zclhsv - va(i,j)
c           va(i,j) = facsor*zforcv/zclhsv + (1.-facsor)*va(i,j)
c           impose (weighted) Schoof velocities:
c           zwei = fschoofv(i,j)
c           va(i,j) = (1.-zwei)*va(i,j) + zwei*vschoof(i,j)
c           if (zwei.eq.0.) zdumax = max (zdumax, abs(zdv))
c@@@ new (see comments above for u):
            zs = fschoofv(i,j)
            zf = (1.-zs)*facsor + zs*1.
            vaold = va(i,j)
            va(i,j) = zf * ((1.-zs)*(zforcv/zclhsv) + zs*vschoof(i,j))
     *              + (1.-zf) * va(i,j)
            zdv = abs(va(i,j)-vaold) / max (abs(vaold),1.)
            zdumax = max (zdumax, abs(zdv))
c@@@

          endif
        enddo
      enddo
#if ! defined (NOCHECKERBOARD)
      enddo
#endif

c        Print out sor diagnostics

c     if (.true.) then
      if (.false.) then
        if (iterc.eq.1 .and. itera.eq.1 .and. itersur.eq.1 ) 
     *    write (6,'(/25("=")/a,i10/25("="))') 
     *       'dosparse_sor: time=', nint(timeice)

        if (itersur.eq.1) 
     *    write (6,'(/a,3i8,f12.3/a)') 
     *       'iterc, itera, ndav, zuav=', iterc, itera, ndav, zuav,
     *       '  itersur      zdumax      ducrit'

c       if (itersur.eq.1 .or. itersur.eq.nitersurmax .or. 
c    *     (zdumax.lt.ducrit .and. itersur.ge.nitersurmin)) then
          write (6,'(i9, 2f12.3, i9)') itersur, zdumax, ducrit, iflag
c       endif
      endif

c######
c     diagnostic dump of maps around a given point:
c
c     if (.true.) then
c     if (.false.) then
c
c       io1 = 38 ! 25
c       jo1 = 73 ! 85
c       io2 = io1 + 4
c       jo2 = jo1 + 4
c       iu = 110
c
c       if (iterc.eq.1 .and. itera.eq.1 .and. itersur.eq.1 ) then
c         write (iu,'(/12("=")/a,i8/12("="))') 'time=', nint(timeice)
c       endif
c
c       if (itersur.eq.1 ) then
c
c         write (iu,'(/31("-")/a,2i8)') 'iterc, itera = ', iterc, itera
c
c         write (iu,'(/a,i6)') 'maskwater, h:', itersur
c         do jo=jo2,jo1,-1  
c           write (iu,'(5i9, 10x, 5f9.2)') 
c    *        (maskwater(io,jo),io=io1,io2), (h(io,jo),io=io1,io2)
c         enddo
c
c         write (iu,'(/a,i6)') 'hb, hs, hw:', itersur
c         do jo=jo2,jo1,-1  
c           write (iu,'(5f9.2, 10x, 5f9.2, 10x, 5e11.3)') 
c    *        (hb(io,jo),io=io1,io2), (hs(io,jo),io=io1,io2),
c    *        (hw(io,jo),io=io1,io2)
c         enddo
c
c         write (iu,'(/a,i6)') 'muind, mvind:', itersur
c         do jo=jo2,jo1,-1  
c           write (iu,'(5i9, 10x, 5i9)') 
c    *        (muind(io,jo),io=io1,io2), (mvind(io,jo),io=io1,io2)
c         enddo
c
c         write (iu,'(/a,i6)') 'uaorig, vaorig:', itersur
c         do jo=jo2,jo1,-1  
c           write (iu,'(5f9.1, 10x, 5f9.1)') 
c    *        (uaorig(io,jo),io=io1,io2), (vaorig(io,jo),io=io1,io2)
c         enddo
c
c       endif
c
cc      if (itersur.eq.1 .or. itersur.eq.nitersurmax .or. 
cc   *     (zdumax.lt.ducrit .and. itersur.ge.nitersurmin)) then
cc        write (iu,'(/a,i6)') 'ua, va:', itersur
cc        do jo=jo2,jo1,-1  
cc          write (iu,'(5f9.1, 10x, 5f9.1)') 
cc   *        (ua(io,jo),io=io1,io2), (va(io,jo),io=io1,io2)
cc        enddo
cc      endif
c
c       if (zdumax.gt.100000.) then
c         write (6,'(/a)') 'STOP...large zdumax (dosparse_sor, icedyn)'
c         stop
c       endif
c
c     endif
c######

c     skip out of sor iteration if max change is small enough:
      if (zdumax.lt.ducrit .and. itersur.ge.nitersurmin) goto 1000

c>>>>>>>>>>
      enddo ! itersur
c>>>>>>>>>>
 1000 continue

      return
      end

c=====
#endif
c=====

c-----------------------------------------------------------------------

c)))))))))))))))))))))
#if defined (SCHOOFGL)
c)))))))))))))))))))))

      subroutine schoofgl (uschoof, vschoof, fschoofu, fschoofv,
     *                     hgu, hgv, qgu, qgv,
     *                     fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                     angnormu, angnormv,
     *                     h, hs, hb, hw, fedge, maskwater,
     *                     arhap, s2a0,s2b0,siu,siv,lhsu,lhsv, budgall,
     *                     crhu, crhv, dfu, dfv, hu, hv, masku, maskv,
     *                     muind, mvind, fracgu, fracgv, 
     *                     ua,va,uagrid,vagrid,uia,via, thetau, thetav,
#if defined (MELANGE)
     *                     thetaumel, thetavmel, hmel,
#endif
     *                     sealev, timeice, vmax)

c     Does Schoof-based calculations for grounding lines (Schoof, JGR,
c     2007, see 9/27/07 and 11/10/07 notes).
c     (1) Calculates hg[u,v], qg[u,v] at sub-grid g.l. location
c         (given by fracg[u,v], passed). Schoof flux qg[u,v (in calcgl) 
c         is modified by buttressing fraction theta[u,v] (thetacalc).
c     (2) Sets [u,v]schoof velocs, and weightings fschoof[u,v], 
c         for [u,v]-grid point at g.l. or 1 pt downstream,
c         depending on Schoof vs. non-Schoof flux at former pt. 
c         (Non-Schoof flux is explicit SIA calced here for m[u,v]ind=0,
c         or from [u,v]bgrid (passed) for m[u,v]ind >= 1). 
c         These schoof velocs and weightings will be used in next 
c         dotria/dosparse[sor], and icestep. 
c     (3) Also returns various diagnostics:
c         hg[u,v], qg[u,v], at sub-grid g.l. location,
c         and flux* shifted to grid points.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
c       returned:
     *  uschoof(0:nxp,0:nyp),  vschoof(0:nxp,0:nyp),
     *  fschoofu(0:nxp,0:nyp), fschoofv(0:nxp,0:nyp),
     *  hgu(0:nxp,0:nyp),      hgv(0:nxp,0:nyp),
     *  qgu(0:nxp,0:nyp),      qgv(0:nxp,0:nyp),
     *  fluxgrdu(0:nxp,0:nyp), fluxgrdv(0:nxp,0:nyp), 
     *  fluxschu(0:nxp,0:nyp), fluxschv(0:nxp,0:nyp), 
     *  angnormu(0:nxp,0:nyp), angnormv(0:nxp,0:nyp),

c       passed:
     *  h(nx,ny),              hs(nx,ny),
     *  hb(nx,ny),             hw(nx,ny),       
     *  fedge(nx,ny),          maskwater(nx,ny),    
     *  arhap(nx,ny),          
     *  s2a0(nx,ny),           s2b0(nx,ny),           
     *  siu(0:nxp,0:nyp),      siv(0:nxp,0:nyp),
     *  lhsu(0:nxp,0:nyp),     lhsv(0:nxp,0:nyp),
     *  budgall(nx,ny),
     *  crhu(0:nxp,0:nyp),     crhv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),   fracgv(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  uagrid(0:nxp,0:nyp),   vagrid(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),      via(0:nxp,0:nyp),
     *  thetau(nx,ny),         thetav(nx,ny)
#if defined (MELANGE)
     * ,thetaumel(nx,ny),      thetavmel(nx,ny),   hmel(nx,ny)
#endif
      real lhsu, lhsv

c     local:
      dimension ifglu(0:nxp), ifglv(0:nyp)

#if defined (MELANGE)
#include <comicepresmel.h>
#endif

c     statement fn for Schoof weighting, dep on Schoof-"free" g.l. flux:
c     weischoof (fluxdif) = max (0., min (1., fluxdif/1.e4 ))  ! m^2 /yr
c     2/13: try to reduce butterfly-effect differences in totice(t) 
c     between machines:
c     weischoof (fluxdif) = max (0., min (1.,  fluxdif/1.e2 + 0.5 )) 
      weischoof (fluxdif) = max (0., min (1.,  fluxdif/1.e5 + 0.5 )) 

      call zero (uschoof,  (nxp+1)*(nyp+1))
      call zero (vschoof,  (nxp+1)*(nyp+1))
      call zero (fschoofu, (nxp+1)*(nyp+1))
      call zero (fschoofv, (nxp+1)*(nyp+1))
      call zero (hgu,      (nxp+1)*(nyp+1))
      call zero (hgv,      (nxp+1)*(nyp+1))
      call zero (qgu,      (nxp+1)*(nyp+1))
      call zero (qgv,      (nxp+1)*(nyp+1))
      call zero (fluxgrdu, (nxp+1)*(nyp+1))
      call zero (fluxgrdv, (nxp+1)*(nyp+1))
      call zero (fluxschu, (nxp+1)*(nyp+1))
      call zero (fluxschv, (nxp+1)*(nyp+1))
      call resetr (angnormu, (nxp+1)*(nyp+1), -999.)
      call resetr (angnormv, (nxp+1)*(nyp+1), -999.)

c        Set thetau,v = 1 - fraction of g.l. longitudinal stress 
c        buttressed, used below (i) to decide whether to impose 
c        schoofgl velocities (schoofgl), and (ii) to affect the 
c        imposed velocity itself (calcgl). Use non-Schoof "grid"
c        velocities, calculated in ASCH loop above.

      call thetacalc (h, fedge, maskwater, arhap, 
     *                uagrid, vagrid, thetau, thetav) ! from non-Schoof
c    *                ua,     va,     thetau, thetav) ! from prev C-iter

      if (powi.eq.3) then
        zpowi = 1.
      else
        zpowi = (powi-1)/2.
      endif

      if (powb.eq.3) then
        zpowb = 1.
      else
        zpowb = (powb-1)/2.
      endif

cuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
c        u grid:
cuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu

c==============
      do j=1,ny
c==============

c          Use ifglu to require at least 3 h-grid floating/water pts
c          next to g.l. If only 2, nearby g.l. calcs with opposite flow
c          direction could impose fluxes twice at same u-pt (if both
c          in advancing case). 

c          If advancing, procedure imposes veloc (ua,va) at 1st u-pt
c          *downstream* of g.l. u-pt. So to test on buttressing,
c          use theta[u,v] from closest floating h cell to g.l.
c          depending on whether u-pts were included in prev 
c          elliptic soln or not (muind, mvind).

        call izero (ifglu, nxp+1)
        do i=2,nx-2
          ip3 = min (i+3,nx)
          im2 = max (i-2,1)
          if (masku(i,j).eq.1) then
            if      ( h(i  ,j).gt.0. .and. maskwater(i  ,j).eq.0
     *                               .and. maskwater(i+1,j).ne.0
     *                               .and. maskwater(i+2,j).ne.0
     *                               .and. maskwater(ip3,j).ne.0 ) then
              ifglu(i) = 1
            else if ( h(i+1,j).gt.0. .and. maskwater(i+1,j).eq.0
     *                               .and. maskwater(i  ,j).ne.0
     *                               .and. maskwater(i-1,j).ne.0
     *                               .and. maskwater(im2,j).ne.0 ) then 
              ifglu(i) = -1
            endif
          endif
        enddo

c>>>>>>>>>>>>>>>>>>
        do i=2,nx-2
c>>>>>>>>>>>>>>>>>>
          if (iabs(ifglu(i)).eq.1) then

c              Set indices for points around g.l. On h grid,
c              (ia,j) is grounded with ice, (ib,j) is downstream
c              ocean or *open* lake. On u grid, (iau,j) contains g.l.,
c              and (ibu,j) is next downstream pt.

            if (ifglu(i).eq.1) then
              ia  = i
              ib  = i+1
              iau = i
              ibu = i+1
              idir = 1
            else
              ia  = i+1
              ib  = i
              iau = i
              ibu = i-1
              idir = -1
            endif
#if defined (MELANGE)
            if (fedge(ib,j).eq.0. .and. h(ia,j).gt.0.) then
              zr = 
     *          (0.5*rhomp*grav*(hmel(ib,j)**2) + presmel(hmel(ib,j)))
     *        / (0.5*rhoip*grav*(h(ia,j)**2))
              thetau(ib,j) = 1. - min(zr,1.) * (1.-thetaumel(ib,j))
              thetau(ib,j) = max (0., min (1., thetau(ib,j)))
            endif
#endif
            ztheta = thetau(ib,j)

c              Only do if marine ice is not abutting land, 
c              and g.l. is  not completely buttressed

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           if (hs(ia,j).gt.hs(ib,j) .and. ztheta.gt.0.) then  ! was 0.3
            if (hs(ia,j).gt.hs(ib,j)                   ) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c                Set needed Schoof-like values for calcgl:

              zasch = 0.5*(arhap(ia,j)**powi+arhap(ib,j)**powi)
              zbsch = (1./crhu(iau,j))**(1./powb)
              zfrac = fracgu(iau,j)
              zdx = 0.5*(dx(ia,j) + dx(ib,j))

c                Calculates hg, ug (g.l. ice thickness, schoof veloc)


              call calcgl ( h(ia,j), hs(ia,j), hb(ia,j),
     *                      h(ib,j), hs(ib,j), hb(ib,j),
     *                      zasch, zbsch, zfrac, zdx,
     *                      sealev, ztheta, zhg, zfg)

c                Calculate angle between actual normal to grounding line
c                and x axis (angnorm), and reduce flux accordingly

              call normcalc (iau, j, maskwater, angnorm, 1)
#if defined (ANGNORM)
              zfg = zfg * abs(cos(angnorm))
#endif

c                If requested, correct Schoof flux for fraction of 
c                budgall on first-floating h-cell between sub-grid g.l.
c                location (fracg[u,v]) and u-grid interface where flux 
c                is applied (iau or ibu, zfluxsa or zfluxsb). This 
c                assumes that dynamic flux variations are close to 
c                equilibrium  with local budget on small (sub-grid) 
c                scales. It yields equilibrated grounding-line locations
c                closer to the analytic Schoof solns in idealized tests
c                (MISMIP, Docquier/Pattyn; see notes 12/8/10, 12/12/10).
c                nb: Not valid for huge unrealistic budgall values
c                (MAXSHELFMELT, NOSHELF).

#if defined (SCHOOFADJUST)
c             zfluxsa = max ( 0., 
c    *                        zfg + (0.5-zfrac)*budgall(ib,j)*dx(ib,j) )
c             zfluxsb = max ( 0., zfluxsa + budgall(ib,j)*dx(ib,j) )
c or more exact area-accounting for non-Cartesian grids: 
              dyg = (1.-zfrac)*dy(ia,j) + zfrac*dy(ib,j) 
              zdar = 0.5*(dyu(iau,j)+dyg) * (zfrac-0.5)*dxu(iau,j)  
              zfluxsa = max ( 0., zfg*dyg - budgall(ib,j)*zdar )
     *                      / dyu(iau,j)
              zfluxsb = max ( 0., zfluxsa*dyu(iau,j)
     *                            + budgall(ib,j)*darea(ib,j) )
     *                      / dyu(ibu,j)
#else
              zfluxsa = zfg
              zfluxsb = zfg
#endif

c                Calculate Non-schoof "grid" flux at closest 
c                u-grid pt to g.l. (iau), like explicit zf's in 
c                icestep. zfluxg has model's i-index sign 
c                (+ve if flowing to "right", as idir).

              zdfu  = dfu(iau,j) - lhsu(iau,j)
              zdfv  = 0.25
     *              * (  dfv(ia,j) +dfv(ib,j) +dfv(ia,j-1) +dfv(ib,j-1)
     *                 -lhsv(ia,j)-lhsv(ib,j)-lhsv(ia,j-1)-lhsv(ib,j-1))
              zdf2  = max (zdfu**2 + zdfv**2, .001)

c             SIA internal deformation:
              zhu    = hu(iau,j)
              zs2a   = 0.5*(s2a0(ia,j) + s2a0(ib,j)) 
              zs2b   = 0.5*(s2b0(ia,j) + s2b0(ib,j))
              zbra   = zs2a*(zdf2**zpowi) + zs2b*siu(iau,j)
              zfluxg = 2. * zbra * (-zdfu) * zhu**2
              if (muind(iau,j).eq.0) then
c               plus SIA basal sliding:
                zcrh   = crhu(iau,j)
                zbra   = zcrh * (zdf2**zpowb)
                zfluxg = zfluxg + zbra * (-zdfu) *  zhu
              else
c               plus shelf elliptic soln (minus uia), with upstream h:
                zfluxg = zfluxg + (uagrid(iau,j)-uia(iau,j))*h(ia,j)
              endif

c                Impose Schoof velocity, weighted at closest g.l. 
c                u-grid pt (iau) or 1 point downstream (ibu),
c                as function (weischoof) of difference between
c                Schoof |flux| at iau vs. Schoof-free model |flux| 
c                at iau (zfluxdif). More weighting for iau if Schoof 
c                |flux| > model |flux| (g.l. retreat), more for ibu
c                if <  (g.l. advance). Sum of the two weights = 1.

              zfluxdif = zfluxsa - abs(zfluxg)
              fschoofu(iau,j) = weischoof(zfluxdif)  
              fschoofu(ibu,j) = 1. - weischoof(zfluxdif)  

c                Convert Schoof flux to veloc using upstream h-grid h,
c                to match upstream advection for it in icestep.
c                Limit |uschoof| for stability, in case h is v. small.

              uschoof(iau,j) = 
     *          idir * min (vmax, zfluxsa/max(h(ia,j),.01))
              uschoof(ibu,j) = 
     *          idir * min (vmax, zfluxsb/max(h(ib,j),.01))

c                Set diagnostic g.l. quantities, returned:

              hgu(iau,j) = zhg
              qgu(iau,j) = zfg
              fluxschu(iau,j) = idir*zfluxsa
              fluxgrdu(iau,j) = zfluxg
              angnormu(iau,j) = angnorm*180./pi
c~~~~~~~~~~~~~~~~
            endif
c~~~~~~~~~~~~~~~~

          endif
c>>>>>>>>>>>>
        enddo
c>>>>>>>>>>>>
c==========
      enddo
c==========

cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c        v grid (same comments as above for u grid):
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

c==============
      do i=1,nx
c==============

        call izero (ifglv, nyp+1)
        do j=2,ny-2
          jp3 = min (j+3,ny)
          jm2 = max (j-2,1)
          if (maskv(i,j).eq.1) then
            if      ( h(i  ,j).gt.0. .and. maskwater(i,j  ).eq.0
     *                               .and. maskwater(i,j+1).ne.0
     *                               .and. maskwater(i,j+2).ne.0
     *                               .and. maskwater(i,jp3).ne.0 ) then
              ifglv(j) = 1
            else if ( h(i,j+1).gt.0. .and. maskwater(i,j+1).eq.0
     *                               .and. maskwater(i,j  ).ne.0
     *                               .and. maskwater(i,j-1).ne.0
     *                               .and. maskwater(i,jm2).ne.0 ) then
              ifglv(j) = -1
            endif
          endif
        enddo

c<<<<<<<<<<<<<<<<<<
        do j=2,ny-2
c<<<<<<<<<<<<<<<<<<
          if (iabs(ifglv(j)).eq.1) then

            if (ifglv(j).eq.1) then
              ja  = j
              jb  = j+1
              jau = j
              jbu = j+1
              idir = 1
            else
              ja  = j+1
              jb  = j
              jau = j
              jbu = j-1
              idir = -1
            endif
#if defined (MELANGE)
            if (fedge(i,jb).eq.0. .and. h(i,ja).gt.0.) then
              zr = 
     *          (0.5*rhomp*grav*(hmel(i,jb)**2) + presmel(hmel(i,jb)))
     *        / (0.5*rhoip*grav*(h(i,ja)**2))
              thetav(i,jb) = 1. - min(zr,1.) * (1.-thetavmel(i,jb))
              thetav(i,jb) = max (0., min (1., thetav(i,jb)))
            endif
#endif
            ztheta = thetav(i,jb)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c           if (hs(i,ja).gt.hs(i,jb) .and. ztheta.gt.0.) then  ! was 0.3
            if (hs(i,ja).gt.hs(i,jb)                   ) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              zasch = 0.5*(arhap(i,ja)**powi+arhap(i,jb)**powi)
              zbsch = (1./crhv(i,jau))**(1./powb)
              zfrac = fracgv(i,jau)
              zdx = 0.5*(dy(i,ja) + dy(i,jb))

              call calcgl ( h(i,ja), hs(i,ja), hb(i,ja),
     *                      h(i,jb), hs(i,jb), hb(i,jb),
     *                      zasch, zbsch, zfrac, zdx,
     *                      sealev, ztheta, zhg, zfg)

              call normcalc (i, jau, maskwater, angnorm, 0)
#if defined (ANGNORM)
              zfg = zfg * abs(sin(angnorm))
#endif

#if defined (SCHOOFADJUST)
c             zfluxsa = max ( 0., 
c    *                        zfg + (0.5-zfrac)*budgall(i,jb)*dy(i,jb) )
c             zfluxsb = max ( 0., zfluxsa + budgall(i,jb)*dy(i,jb) )
              dxg = (1.-zfrac)*dx(i,ja) + zfrac*dx(i,jb) 
              zdar = 0.5*(dxv(i,jau)+dxg) * (zfrac-0.5)*dyu(i,jau)  
              zfluxsa = max ( 0., zfg*dxg - budgall(i,jb)*zdar )
     *                      / dxv(i,jau)
              zfluxsb = max ( 0., zfluxsa*dxv(i,jau)
     *                            + budgall(i,jb)*darea(i,jb) )
     *                      / dxv(i,jbu)
#else
              zfluxsa = zfg
              zfluxsb = zfg
#endif

              zdfu  = 0.25
     *              * ( dfu(i,ja) +dfu(i,jb) +dfu(i-1,ja) +dfu(i-1,jb)
     *                -lhsu(i,ja)-lhsu(i,jb)-lhsu(i-1,ja)-lhsu(i-1,jb))
              zdfv  = dfv(i,jau) - lhsv(i,jau)
              zdf2  = max (zdfu**2 + zdfv**2, .001)

c             SIA internal deformation:
              zhv    = hv(i,jau)
              zs2a   = 0.5*(s2a0(i,ja) + s2a0(i,jb))
              zs2b   = 0.5*(s2b0(i,ja) + s2b0(i,jb))
              zbra   = zs2a*(zdf2**zpowi) + zs2b*siv(i,jau)
              zfluxg = 2. * zbra * (-zdfv) * zhv**2
              if (mvind(i,jau).eq.0) then
c               plus SIA basal sliding:
                zcrh   = crhv(i,jau)
                zbra   = zcrh * (zdf2**zpowb)
                zfluxg = zfluxg + zbra * (-zdfv) * zhv
              else
c               plus shelf elliptic soln (minus via), with upstream h:
c bug fix 07/17/15:
c               zfluxg = zfluxg + (vagrid(i,jau)-via(i,jau))*h(ia,j)!bug
                zfluxg = zfluxg + (vagrid(i,jau)-via(i,jau))*h(i,ja) !ok
              endif

              zfluxdif = zfluxsa - abs(zfluxg)
              fschoofv(i,jau) = weischoof(zfluxdif)
              fschoofv(i,jbu) = 1. - weischoof(zfluxdif)

              vschoof(i,jau) =
     *          idir * min (vmax, zfluxsa/max(h(i,ja),.01))
              vschoof(i,jbu) =
     *          idir * min (vmax, zfluxsb/max(h(i,jb),.01))

              hgv(i,jau) = zhg
              qgv(i,jau) = zfg
              fluxschv(i,jau) = idir*zfluxsa
              fluxgrdv(i,jau) = zfluxg
              angnormv(i,jau) = angnorm*180./pi
c~~~~~~~~~~~~~~~~
            endif
c~~~~~~~~~~~~~~~~

          endif 
c<<<<<<<<<<<<
        enddo
c<<<<<<<<<<<<
c==========
      enddo
c==========

      return
      end

c-----------------------------------------------------------------------

      subroutine normcalc (iuc, juc, maskwater, angnorm, ifu)

c     Calculates angle (angnorm) between actual normal to grounding 
c     line and the x axis, for a single u-grid (ifu=1) or v-grid 
c     (ifu=0) point given by iuc,juc. The actual normal is estimated
c     using Ben P.'s algorithm! (Boulder, 2/12): averages the vectors 
c     (zdx,zdy) from the point to the h-grid centers of nearby
c     ocean points (maskwater ne 0), looking only in a box of
c     half width nwin (few tens of km), and weighted by 1/dist. 
c     Then the normal is just the arc tangent of the average (zdx,zdy).
c
c     Called from subr schoofgl. The u-grid Schoof flux is reduced there
c     by abs(cos(angnorm)), and the v-grid flux by abs(sin(angnorm))
c     (if defined ANGNORM).
c     For a long straight grounding line with (shallow) angle A to 
c     x axis, pixelated by cells that step one grid box perpendicularly
c     for many cells along, the resulting mean flux across an actual 
c     length L of grounding line calculated by this method is:
c       (L*cos(A)/dx)*F*cos(A)*dx + (L*sin(A)/dx)*F*sin(A)*dx = L*F
c     which is exactly what it should be for the actual segment L.

#include <comicephys.h>
#include <comicegrid.h>

      dimension maskwater(nx,ny)

      nwin = max (3, nint(15.e3/dd0))         ! half-width of window box

c        Set (zic,zjc) = h-grid coords of u- or v-grid point (iuc,juc).
c        Also set (ia,ib,ja,jb) =  window box limits, h grid.

      if (ifu.eq.1) then
        zic = iuc + 0.5
        zjc = juc
        ia = iuc   - (nwin-1)
        ib = iuc+1 + (nwin-1)
        ja = juc - nwin
        jb = juc + nwin
      else 
        zic = iuc
        zjc = juc + 0.5
        ia = iuc - nwin
        ib = iuc + nwin
        ja = juc   - (nwin-1)
        jb = juc+1 + (nwin-1)
      endif

c       Average vectors (zdx,zdy) from u/v point to ocean (not land) 
c       h-grid centers within window box, weighted by 1/distance 
c       (1/zwei). If at edge of domain, assume land-ocean map is
c       extended from last edge row or column.  

      zavx = 0.
      zavy = 0.
      zava = 0.
      do j=ja,jb
        jj = max (1, min (ny, j))
        do i=ia,ib
          ii = max (1, min (nx, i))
          if (maskwater(ii,jj).ne.0) then
            zdx = (i - zic) * dd0
            zdy = (j - zjc) * dd0
            zwei = sqrt(zdx**2 + zdy**2) 
            zavx = zavx + zdx / zwei 
            zavy = zavy + zdy / zwei
            zava = zava + 1.  / zwei
          endif
        enddo
      enddo
      if (zava.eq.0.) then
        write (6,'(a,2i6)') 'Error (normcalc): zava=0. iuc,juc=',iuc,juc
        stop
      endif
      zavx = zavx/zava
      zavy = zavy/zava

c       Set angle (counterclockwise) from x-axis to the estimated normal

      angnorm = atan2(zavy,zavx)
      
      return
      end

c-----------------------------------------------------------------------

      subroutine calcgl ( h1, hs1, hb1,
     *                    h2, hs2, hb2,
     *                    acoef, bcoef, frac, dx,
     *                    sealev, theta, hg, fg)

c     Calculates hg (g.l. depth) and fg (Schoof flux = hg*ug).
c     h1,hs1,hb1 are upstream, and h2,hs2,hb2 are downstream (all supp),
c     regardless of actual grid direction. Accounts for buttressing
c     fraction theta.

#include <comicephys.h>
c     statement function for grid resol/non-resol term, to
c     eliminate explicit stretching perturbation if grid resolves
c     g.l. (see notes 9/27/07).
c     fresol(x) = 1.-exp(-(x**2)) 
      fresol(x) = 1.   ! or not

c     calc g.l. depth (linear interp to u-grid of bed elev, and s.l.)
c     hbg =       0.5*hb1 +  0.5*hb2
      hbg = (1.-frac)*hb1 + frac*hb2      ! (a bit better 100 vs 300 ka)
      hg = max ((sealev-hbg)/rhor, 0.)
      hsg = hbg + hg

c     1/e-folding length of upstream perturbation (stretching) veloc
      zlen = ( acoef * (0.5*bcoef)**powb  * hg**(powi - 2*powb)
     *         * (0.25*rhoip*grav*theta)**(powi-powb)
     *       ) ** (1./(powb+1.))
  
c     perturbation (stretching) veloc
c     ug = (0.25*rhoip*grav*hg*theta)**(powb*(powi+1.)/(powb+1.))
c     bug fix, power of theta:
      ug = (0.25*rhoip*grav*hg)**(powb*(powi+1.)/(powb+1.))
     *     *            (theta)**(powb*(powi   )/(powb+1.))
     *     * (2.*hg*acoef/bcoef)**(powb/(powb+1.))
     *     / (0.5*(1.-rhor))**(powb/(powb+1.)) ! ~Sch eq.16

      ug = ug * fresol(zlen*dx)

#if defined (USCHOOFFAC)
      ug = ug * USCHOOFFAC       ! 0.1 ! 0.3 ! 0.5 ! 2. ! 3. ! 10. ! 30.
#endif

      fg = hg*ug

      return
      end

c-----------------------------------------------------------------------

      subroutine thetacalc (h, fedge, maskwater, arhap, 
     *                      ua, va, thetau, thetav)

c     Sets theta[u,v] = 1 - fraction of g.l. longitudinal stress 
c     buttressed, used in schoofgl and calcgl. 
c
c     Nb: theta[u,v] are on h-grid, for E-W and N-S directions 
c     respectively, and ice thickness is h, so that if all shelf is 
c     free-floating, theta=1 exactly (property of elliptic solution).
c     (And current h is that used in previous ASCH-loop, because
c     it's inside the C-loop).       
c

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),              fedge(nx,ny),
     *  maskwater(nx,ny),      arhap(nx,ny),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  thetau(nx,ny),         thetav(nx,ny)

      parameter (hthetmin = .01)   ! have tried 0.1, 1., 10.

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).ne.0) then
c~~~~~~~~~~~~~~~~~~~~~~~
#if defined (NOBUTTRESS)
c~~~~~~~~~~~~~~~~~~~~~~~
            thetau(i,j) = 1.
            thetav(i,j) = 1.
c~~~~
#else
c~~~~
            if (h(i,j).lt.hthetmin) then
              thetau(i,j) = 1.
              thetav(i,j) = 1.
            else
              zdudx = (ua(i,j)-ua(i-1,j)) / dx(i,j)
              zdvdy = (va(i,j)-va(i,j-1)) / dy(i,j)

c             zdiv = max (0., zdudx+zdvdy)
c             thetau(i,j) = ( zdiv**(1./powi) / arhap(i,j) )
c    *                      * 4. / (rhoip*grav*h(i,j))

              zg1 =  max (zdudx**2 + zdvdy**2 + zdudx*zdvdy, gmin)
              zvisc = 0.5 / (zg1**powiv)

#if defined (COP21RD)
              thetau(i,j)  = ( 4.*zvisc*2.*zdudx/arhap(i,j) )
     *                       / (rhoip*grav*h(i,j))
              thetav(i,j)  = ( 4.*zvisc*2.*zdvdy/arhap(i,j) )
     *                       / (rhoip*grav*h(i,j))
#else
              thetau(i,j)  = ( 2.*zvisc*(2.*zdudx+zdvdy)/arhap(i,j) )
     *                       / (0.5*rhoip*grav*h(i,j))
              thetav(i,j)  = ( 2.*zvisc*(2.*zdvdy+zdudx)/arhap(i,j) )
     *                       / (0.5*rhoip*grav*h(i,j))
#endif

              thetau(i,j) = max (0., min (1., thetau(i,j)))
              thetav(i,j) = max (0., min (1., thetav(i,j)))

              im1 = max (i-1,1)
              ip1 = min (i+1,nx)
              zfed = fedge(i,j)
              if (maskwater(im1,j).ne.0) zfed = min(zfed,fedge(im1,j))
              if (maskwater(ip1,j).ne.0) zfed = min(zfed,fedge(ip1,j))
              thetau(i,j) = (1.-zfed)*1. + zfed*thetau(i,j)

              jm1 = max (j-1,1)
              jp1 = min (j+1,ny)
              zfed = fedge(i,j)
              if (maskwater(i,jm1).ne.0) zfed = min(zfed,fedge(i,jm1))
              if (maskwater(i,jp1).ne.0) zfed = min(zfed,fedge(i,jp1))
              thetav(i,j) = (1.-zfed)*1. + zfed*thetav(i,j)
            endif
c~~~~~
#endif
c~~~~~
          else
            thetau(i,j) = -1.
            thetav(i,j) = -1.
          endif
        enddo
      enddo

      return
      end

c)))))
#endif
c)))))

c-----------------------------------------------------------------------

c((((((((((((((((((((((
#if defined (CLIFFMELT)
c((((((((((((((((((((((

      subroutine cliffcalc (h, hb, fedge, ifnexo, maskwater, 
     *                      cliffmelt, calvice, crevliq, 
     *                      uadv, vadv, fracgu, fracgv, thetau, thetav,
#if defined (MELANGE)
     *                      ucimel, vcimel, hucimel, hvcimel,
#endif
     *                      sealev, timeice, dt)

c     Set cliffmelt, for collapse of > ~100 m vertical exposed tidewater
c     cliff walls, Applied as erosional velocity on entire column
c     face. Use interpolated grounding line to get column-face 
c     dimensions. Rate of cliff failure depends on theta[u,v], and
c     details below depend on [u,v]adv (zug), so cliffcalc is called
c     from icedyn, and is not in iceocean.  
c
c     Can modify calvice - so cliffmelt and calvice are not yet included
c     in budgall (added in it after icedyn is called, in icectl).

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),            hb(nx,ny), 
     *  fedge(nx,ny),        ifnexo(nx,ny),      maskwater(nx,ny),    
     *  cliffmelt(nx,ny),    calvice(nx,ny),     crevliq(nx,ny),
     *  uadv(0:nxp,0:nyp),   vadv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp), fracgv(0:nxp,0:nyp),
     *  thetau(nx,ny),       thetav(nx,ny)
#if defined (MELANGE)
      dimension
     *  ucimel(0:nxp,0:nyp),    vcimel(0:nxp,0:nyp),
     *  hucimel(0:nxp,0:nyp),   hvcimel(0:nxp,0:nyp)
#endif
      parameter (nlook = 4)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
c     E-W/N-S:
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /

      dimension work4(nx,ny)
      save work4

      logical firstcliff
      data firstcliff /.true./
      save firstcliff

#  if defined (RWRAP)
      cvmax = rwrap_cliffvmax
#  elif defined (CLIFFVMAX)
      cvmax = CLIFFVMAX
#  else
#    if defined (PREAGU13)
      cvmax = 12.e3  ! pre AGU13 
#    else
      cvmax =  3.e3  ! post AGU13 
#    endif
#  endif

#  if defined (CLIFFZMAX)
      czmax = CLIFFZMAX
#  else
#    if defined (PREAGU13)
      czmax =  90.  ! pre AGU13 
#    else
      czmax = 100.  ! post AGU13 
#    endif
#  endif

#  if defined (PELTOVMAX)
      pvmax = PELTOVMAX
#  else
#    if defined (PREAGU13)
      pvmax = cvmax  ! pre AGU13
#    else
      pvmax = 0.     ! post AGU13
#    endif
#  endif

#  if defined (PELTOSLOPE)
      pslope = PELTOSLOPE
#  else
      pslope = 8.33           ! v = 70 + 8.33*hw, Pelto and Warren, 1991
#  endif

      if (firstcliff) then
        call zero (work4, nx*ny)     ! diagnostic output
        firstcliff = .false.
      endif

      call zero (cliffmelt, nx*ny)
#if defined (MELANGE)
      call zero (ucimel, (nxp+1)*(nyp+1)) 
      call zero (vcimel, (nxp+1)*(nyp+1)) 
      call zero (hucimel, (nxp+1)*(nyp+1)) 
      call zero (hvcimel, (nxp+1)*(nyp+1)) 
#endif

      do j=1,ny
        do i=1,nx
c...................................................................
          if ( h(i,j).gt.0. .and.
     *           (maskwater(i,j).eq.0 .and. sealev.gt.hb(i,j))
     *       )  then
c...................................................................

            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (.not.(ii.eq.i .and. jj.eq.j)) then      ! for flowline
#  if defined (PREAGU13)
c               pre AGU13 (could be moved back to iceocean):
                if (maskwater(ii,jj).eq.1 .and. fedge(ii,jj).eq.0.) then
#  else
c               "zugx":
                if (maskwater(ii,jj).eq.1) then
#  endif

                  if (j.eq.jj) then
                    zfgl = fracgu (min(i,ii), j)
                    zug = abs (uadv (min(i,ii), j)) / darea(ii,jj)
     *                    + h(ii,jj)/dt
                    zdd = dyu (min(i,ii), j)
                    ztheta = thetau(ii,jj)
                  else
                    zfgl = fracgv (i, min(j,jj))
                    zug = abs (vadv (i,  min(j,jj))) / darea(ii,jj)
     *                    + h(ii,jj)/dt
                    zdd = dxv (i, min(j,jj))
                    ztheta = thetav(ii,jj)
                  endif
                  hbg = hb(i,j)*(1.-zfgl) + hb(ii,jj)*zfgl
                  zhw = max (0., sealev-hbg)
                  zhs = zhw * ((1./rhor)-1.)
                  zh  = zhw *  (1./rhor)

#  if defined (PREAGU13)
c                 pre AGU13:
                  zhw2 = zhw
                  zhs2 = zhs
                  zh2  = zh
#  else
c                 "zugx":
c                 increase zh2 due to (i) theta increasing back stress
c                 (ztheta in numerator), (ii) theta decreasing 
c                 divergence and so (dry) crevass depth (as in calvice),
c                 (-.5*ztheta in denominator), and (iii) extra crevasse
c                 depth due to surface liquid (crevliq passed from
c                 calving/ocean, -zfl in denominator):
                  zfl = min ( 1., 0.5*(crevliq(i,j)+crevliq(ii,jj))
     *                            / max (zh, .001) )
                  zhw2 = zhw * ztheta / max(.1, 2.*(1.-.5*ztheta - zfl))
                  zhs2 = zhs * ztheta / max(.1, 2.*(1.-.5*ztheta - zfl))
                  zh2  = zh  * ztheta / max(.1, 2.*(1.-.5*ztheta - zfl))
#  endif

                  zclifv = cvmax * max (0., min (1., (zhs2-czmax)/20.))
                  zpelto = min (70. + pslope*zhw2, pvmax)    ! P+W, 1991
c                 zpelto = min (      pslope*zhw2, pvmax)
                  zcombv = max (zclifv, zpelto)
#  if defined (CLIFFMULT)
                  zcombv = zcombv * CLIFFMULT
#  endif

c                 Ok to accumulate cliffmelt here because each 
c                 (i,j), (ii,jj) pair is only encountered once (gr,fl).
c                 cliffmelt = "edge veloc" *thickess*width / cell area:

#  if defined (PREAGU13)
c                 pre AGU13:
                  cliffmelt(i,j) = cliffmelt(i,j)
     *                           + zcombv*zh*zdd/darea(i,j)
#  else
c                 "zugx":
c                 cliffmelt erodes some or all of Schoof flux crossing
c                 gl (zug), applied to (ii,jj), but no more than would
c                 completely erode ice in (ii,jj). Remainder of
c                 cliffmelt if any is applied to (i,j). 
c
c                 This is consistent with imagining terminus is still 
c                 a tidewater face even if fractional floating ice in 
c                 (ii,jj), as if that ice is grounded. 
c
c                 Even if ice shelf extends further, still apply 
c                 cliffmelt at g.l. (+1), within existing ice shelf 
c                 (strange, but included in Bassis+Walker concept)

                  cliffmelt(i,j) = cliffmelt(i,j)
     *                      + max (zcombv*zh*zdd/darea(i,j) - zug, 0.)
                  cliffmelt(ii,jj) = cliffmelt(ii,jj)
     *                      + min (zcombv*zh*zdd/darea(ii,jj), zug)
#  endif

#  if defined (MELANGE)
c                 set cliff erosion rates, u,v grid, for melange module:
                  if (i.ne.ii) then
                    ucimel(min(i,ii),j) =  zcombv
                    hucimel(min(i,ii),j) =  zh
                  else if (j.ne.jj) then
                    vcimel(i,min(j,jj)) =  zcombv
                    hvcimel(i,min(j,jj)) =  zh
                  endif
#  endif
                  work4(i,j) = max (work4(i,j), zhs)

                endif
              endif
            enddo

c..............
          endif
c..............
        enddo
      enddo

c     if (mod(timeice+.001,100.).lt.0.05) then
c       call printmap (timeice,cliffmelt,'cliffmelt (m/y)',  20., 149,0)
c       call printmap (timeice,work4,   'zhsmax(m)',         10., 149,0)
c       call printmap (timeice,h,       'ice thickness (m)', 150.,149,0)
c       call printmap (timeice,fedge,   'fedge',             .05, 149,0)
c       call flush (149)
c       call zero (work4, nx*ny)
c     endif

      return
      end

c(((((
#endif
c(((((

c-----------------------------------------------------------------------

      subroutine icestep (hold, h, hnew, hw, maskwater,
     *                    hu,  hv,  dhsu,  dhsv, dfu, dfv,
     *                    lhsu, lhsv,
     *                    s1a0, s1b0, s2a0, s2b0, siu, siv,
     *                    crhu, crhv,
     *                    uia, via, uit, vit, ub, vb, ua, va, 
     *                    uschoof, vschoof, fschoofu, fschoofv,  
     *                    masku, maskv, muind, mvind, medgeu, medgev,dt)

c        Step ice thickness through one timestep, from hold (beginning 
c        of timestep) to hnew (end of timestep, returned).
c        For all pts:
c          Does SIA internal deformation advection (fully time-implicit
c          Newton-Raphson ADI + tridia, see notes 2/9/11). Sets uia,via
c          for SSA solns, uit,vit for cross-softening in shelf coeffs.
c        For non-SSA (non-elliptic eqn) pts, i.e., m[u,v]ind = 0:
c          Adds SIA basal sliding (fully time-implicit in h advection).
c          Sets diagnostic ub,vb.
c        For SSA (elliptic eqn) pts, i.e., m[u,v]ind >= 1:
c          Adds prev elliptic ua,va with upstream (but time-implicit) h,
c          minus uia,via to represent basal sliding. 
c          Sets diagnostic ub,vb.
c
c        For closest pt and 1 pt downstream of grounding line (where
c        fschoof[u,v] are non zero), weights bewtween the above 
c        solutions, and [u,v]schoof with upstream (but time-implicit)
c        h. This assumes Schoof(2007) analysis yields g.l. flux and
c        velocity that includes internal deformation, although it 
c        isn't in the Schoof model. (So for idealized 1-D Mismip-like 
c        tests, need to set enhancesheet=0).

c        Sheet-pts ADI method operates on each horizontal direction 
c        in turn, calling subr ice to step the 1-D ice flow equation. 
c        Units are meters and years. The 1-D strip arrays *(n1d) are 
c        used both for both directions in the domain.

c        For ua advection (sliding/sed def/shelf),
c        use an upstream scheme for h interp (or could use cubic SLT?), 
c        since any centered-interp scheme gives spurious buildup
c        at the last non-zero grid box of a flat uniform ice shelf.

#include <comicephys.h>
#include <comicegrid.h>

c        hold  = initial thickness (m) (h-grid, supplied)
c        h     = mid-timestep ice-sheet thickness (m) (supplied)
c        hnew  = new ice-sheet thickness (m) (h-grid, returned)
c        hw    = liquid water thickness (m) (h-grid, supplied
c        uia = depth-averaged internal-shear u veloc (u-grid, returned)
c        via = depth-averaged internal-shear v veloc (v-grid, returned)
c        dt = ice-model timestep (yr) (supplied)

      dimension
     *  hold(nx,ny),           h(nx,ny),          hnew(nx,ny),
     *  hw(nx,ny),             maskwater(nx,ny),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  dhsu(0:nxp,0:nyp),     dhsv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),
     *  lhsu(0:nxp,0:nyp),     lhsv(0:nxp,0:nyp),

     *  s1a0(nx,ny),           s1b0(nx,ny),
     *  s2a0(nx,ny),           s2b0(nx,ny),
     *  siu(0:nxp,0:nyp),      siv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),     crhv(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),      via(0:nxp,0:nyp),
     *  uit(0:nxp,0:nyp),      vit(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),       vb(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  uschoof(0:nxp,0:nyp),  vschoof(0:nxp,0:nyp),
     *  fschoofu(0:nxp,0:nyp), fschoofv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  medgeu(0:nxp,0:nyp),   medgev(0:nxp,0:nyp)
       real lhsu, lhsv

c     local:
      dimension 
     *  pconv(nx,ny)     ! tendency from first ADI stage (m/yr) (h-grid)

c        1-D strip arrays:

      parameter (n1d=nx+ny)
      dimension
     *  hold1(n1d),   hnew1(n1d),   darea1(n1d),
     *  pconv1(n1d),  facfloat(n1d),
     *  e(0:n1d),     f(0:n1d),     g(0:n1d),      flux(0:n1d),
     *  ubsp(0:n1d),  ubsm(0:n1d),  fluxub(0:n1d)

      call scopy (nx*ny, hold, 1, hnew, 1)

#if defined (NESTING)
c     for nesting, only set uia,uit,ub to zero in interior
      do j=0,nyp
        do i=0,nxp
          if ((i.ge.1 .and.i.le.nx-1).and.(j.gt.1.and.j.lt.ny))then
            uia(i,j) = 0.
            uit(i,j) = 0.
            ub(i,j) = 0.
          endif
          if ((j.ge.1 .and.j.le.ny-1).and.(i.gt.1.and.i.lt.nx))then
            via(i,j) = 0.
            vit(i,j) = 0.
            vb(i,j) = 0.
          endif
        enddo
      enddo
#else
      call zero (uia,  (nxp+1)*(nyp+1))
      call zero (via,  (nxp+1)*(nyp+1))
      call zero (uit,  (nxp+1)*(nyp+1))
      call zero (vit,  (nxp+1)*(nyp+1))
      call zero (ub,   (nxp+1)*(nyp+1))
      call zero (vb,   (nxp+1)*(nyp+1))
#endif

      if (powi.eq.3) then
        zpowi = 1.
        zpowim1 = 0.
      else
        zpowi = (powi-1)/2.
        zpowim1 = zpowi - 1.
      endif

      if (powb.eq.3) then
        zpowb = 1.
        zpowbm1 = 0.
      else
        zpowb = (powb-1)/2.
        zpowbm1 = zpowb - 1.
      endif

c        First ADI pass, doing flow for all E-W strips. 
c        Calculate tendency (pconv) due to E-W flow, 
c        don't change 2-D ice thickness yet.

      call zero (pconv, nx*ny)

c~~~~~~~~~~~~~~~~~~~~~~
      if (nx.gt.1) then
c~~~~~~~~~~~~~~~~~~~~~~

#if defined (OMPSTEP)
!$OMP PARALLEL DO
!$OMP* shared (hold, h, hnew, hw, maskwater,
!$OMP*         hu,hv, dhsu,dhsv, dfu,dfv, lhsu,lhsv,
!$OMP*         s1a0,s1b0,s2a0,s2b0, siu,siv,
!$OMP*         uia, uit, ua, ub, dt, pconv,
!$OMP*         uschoof, fschoofu,
!$OMP*         darea, dxu, dxv, dyu, dyv,
!$OMP*         zpowi, zpowim1, zpowb, zpowbm1)
!$OMP* private (hold1, hnew1, darea1, pconv1, facfloat,
!$OMP*          e, f, g, flux, ubsp, ubsm, fluxub,
!$OMP*          i, j, ig, ip, jg, jp, n1, n2, 
!$OMP*          zdfu, zdfv, zdf2, zdhsu, zdhsv, zhu, 
!$OMP*          zs1a, zs1b, zs2a, zs2b,
!$OMP*          zbra, zf, zdfdh, zdfda1, zdfda2, zdh)
!$OMP* schedule(dynamic,2)
#endif

c>>>>>>>>>>>>>>
#if defined (NESTING)
      do j=2,ny-1
#else
      do j=1,ny
#endif
c>>>>>>>>>>>>>>

c        Find low-index edge of ice sheet (n1 is adjacent box with h=0).
c        If non-zero ice exists in first box, set n1=1, and b.c. will
c        implicitly be no flow out of left-hand edge of box.

        n1 = -1
        do i = 1,nx
          if (hold(i,j) + dt*max(pconv(i,j),0.).gt.0.) then
            n1 = max (i-1, 1)
            goto 10
          endif
        enddo
   10   if (n1.eq.-1) n1 = 1

c        Find high-index edge of ice sheet(n2 is adjacent box with h=0).
c        If non-zero ice exists in last box, set n2=nx, and b.c. will
c        implicitly be no flow out of right-hand edge of box.

        n2 = -1
        do i = nx,1,-1
          if (hold(i,j) + dt*max(pconv(i,j),0.).gt.0.) then
            n2 = min (i+1, nx)
            goto 20
          endif
        enddo
   20   if (n2.eq.-1) n2 = 1
#if defined (NESTING)        
        n1 = max (n1, 2)
        n2 = min (n2, nx-1)
#endif

        do i=1,nx
          hold1(i) = hold(i,j)
          pconv1(i)= pconv(i,j)
          darea1(i)= darea(i,j)
          facfloat(i) = 1.
          if (maskwater(i,j).ne.0) facfloat(i) = 1.-rhor 
        enddo

        call zero (e, n1d)
        call zero (f, n1d)
        call zero (g, n1d)
        call zero (flux, n1d)
        call zero (ubsp, n1d)
        call zero (ubsm, n1d)
        call zero (fluxub, n1d)

c*************************
        if (n2.gt.n1) then
c*************************

c          Set 1-D strips (E-W)

        do i=n1-1,n2
          ig = max(i,1)
          ip = min(i+1,nx)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (masku(i,j).eq.1) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c              Use quantities at beginning of timestep for time-implicit
c              tridiag calcs (hu0, dhsu0, dfu0 passed to icestep):

c           longitudinal-stress force (lhs), strain softening (beta2),
c           set to zero currently:

            zdfu  = dfu(i,j) - lhsu(i,j)
            zdfv  = 0.25 
     *            * ( dfv(i,j) + dfv(i+1,j) + dfv(i,j-1) + dfv(i+1,j-1)
     *              -lhsv(i,j) -lhsv(i+1,j) -lhsv(i,j-1) -lhsv(i+1,j-1))
            zdf2  = max (zdfu**2 + zdfv**2, .001)
            zdhsu = dhsu(i,j) 
            zdhsv = 0.25 * (  dhsv(i,j)   + dhsv(i+1,j)
     *                      + dhsv(i,j-1) + dhsv(i+1,j-1) )
            zhu   = hu(i,j)
            zdh   = (facfloat(ip)*hold(ip,j) - facfloat(ig)*hold(ig,j))
     *            / dxu(i,j)

c              Internal deformation (see notes 2/9/11):

            zs1a = 0.5*(s1a0(ig,j) + s1a0(ip,j))
            zs1b = 0.5*(s1b0(ig,j) + s1b0(ip,j))
            zs2a = 0.5*(s2a0(ig,j) + s2a0(ip,j))
            zs2b = 0.5*(s2b0(ig,j) + s2b0(ip,j))

            zbra   = zs2a*(zdf2**zpowi) + zs2b*siu(i,j)
            zf     = 2. * zbra * zdfu * zhu**2
            zdfdh  = 4. * zbra * zdfu * zhu
            zdfda1 = 2.*zbra*zhu**2 
     *             + 4.*zs2a*zpowi*(zdfu**2)  *(zhu**2)*(zdf2**zpowim1)
            zdfda2 = 4.*zs2a*zpowi*(zdfu*zdfv)*(zhu**2)*(zdf2**zpowim1)

            e(i)    = (1.-fschoofu(i,j))
     *                * ( - zdfdh*zhu 
     *                    - zdfda1*rhoice*grav*zhu*(zdhsu+zdh)
     *                    - zdfda2*rhoice*grav*zhu*zdhsv
     *                  ) * dyu(i,j)
            f(i)    = (1.-fschoofu(i,j))
     *                * (   zdfdh 
     *                    + zdfda1*rhoice*grav*zdhsu
     *                    + zdfda2*rhoice*grav*zdhsv
     *                  ) *dyu(i,j)
            g(i)    = (1.-fschoofu(i,j))
     *                * ( zdfda1*rhoice*grav*zhu/dxu(i,j) 
     *                  ) * dyu(i,j)
            flux(i) = (1.-fschoofu(i,j))
     *                * ( -zf
     *                  ) * dyu(i,j)

c           set uia,uit, used in icedyn and setcoefs (softening) 
            zuia = (1.-fschoofu(i,j))*(-zf)/max(zhu,1.e-20)
#if defined (UIACROSS)
            uia(i,j) = zuia
#endif
#if ! defined (NOSOFTCROSS)
            uit(i,j) = -2.*zdfu*zhu * (zdf2*zs1a + siu(i,j)*zs1b)
#endif

c              Plus basal sliding or elliptic (see notes 2/9/11):

c------------------------------------
c           sheet pts, not shelf:
            if (muind(i,j).eq.0) then
c------------------------------------
  
              zcrh = crhu(i,j)
              zbra = zcrh * (zdf2**zpowb)
              zf   = zbra * zdfu * zhu
              zdfdh = zbra * zdfu
              zdfda1 = zbra*zhu
     *               + 2.*zcrh*zpowb*(zdfu**2)  *zhu*(zdf2**zpowbm1)
              zdfda2 = 2.*zcrh*zpowb*(zdfu*zdfv)*zhu*(zdf2**zpowbm1)

              e(i)    = e(i)
     *                + (1.-fschoofu(i,j))
     *                  * ( - zdfdh*zhu 
     *                      - zdfda1*rhoice*grav*zhu*(zdhsu+zdh)
     *                      - zdfda2*rhoice*grav*zhu*zdhsv
     *                    ) * dyu(i,j)
              f(i)    = f(i)
     *                + (1.-fschoofu(i,j))
     *                  * (   zdfdh 
     *                      + zdfda1*rhoice*grav*zdhsu
     *                      + zdfda2*rhoice*grav*zdhsv
     *                    ) *dyu(i,j)
              g(i)    = g(i)
     *                + (1.-fschoofu(i,j))
     *                  * ( zdfda1*rhoice*grav*zhu/dxu(i,j) 
     *                    ) * dyu(i,j)
              flux(i) = flux(i)
     *                  + (1.-fschoofu(i,j))
     *                  * ( -zf
     *                    ) * dyu(i,j)

              fluxub(i) = fschoofu(i,j) 
     *                  * (  max(uschoof(i,j),0.)*hold(ig,j) 
     *                     + min(uschoof(i,j),0.)*hold(ip,j)
     *                    ) * dyu(i,j)
              ubsp(i) = fschoofu(i,j) * max(uschoof(i,j),0.) * dyu(i,j)
              ubsm(i) = fschoofu(i,j) * min(uschoof(i,j),0.) * dyu(i,j)

c             diagnostic:
              ub(i,j) = (1.-fschoofu(i,j)) *(-zf)/max(hu(i,j),1.e-20)
     *                +     fschoofu(i,j)  * uschoof(i,j)
#if defined (COP21RD)
              ua(i,j) = uia(i,j) + ub(i,j)
#else
              ua(i,j) = zuia + ub(i,j)
#endif

c-----------------------------------------
c           elliptic pts:
            else if (muind(i,j).ge.1) then
c-----------------------------------------

c             nb: if fschoofu > 0, Schoof velocities are already 
c             imposed in elliptic soln ua.
              ub(i,j) = ua(i,j) - uia(i,j)
              ubsp(i) = max(ub(i,j),0.) * dyu(i,j)
              ubsm(i) = min(ub(i,j),0.) * dyu(i,j)
              fluxub(i) = ubsp(i)*hold(ig,j) + ubsm(i)*hold(ip,j)

c----------------
            endif
c----------------

#if defined (SUBGRIDEDGE) || ! defined (NOSUBGRIDEDGESHELF)
c           no flow *out* of boxes with fedge < 1 (medgeu set above):
            if (medgeu(i,j).eq.0) then
              e(i) = 0.
              f(i) = 0.
              g(i) = 0.
              flux(i) = 0.
              ubsp(i) = 0.
              ubsm(i) = 0.
              fluxub(i) = 0.
              ua(i,j) = 0.
              ub(i,j) = 0.
              uia(i,j) = 0.
              uit(i,j) = 0.
            endif
#endif

c~~~~~~~~~~~~~~
          endif
c~~~~~~~~~~~~~~

        enddo

c Don't do following if nesting (bug fix 5/25/10)
#if ! defined (NESTING)
        e(n1-1) = 0.
        f(n1-1) = 0.
        g(n1-1) = 0.
        flux(n1-1) = 0.
        ubsp  (n1-1) = 0.
        ubsm  (n1-1) = 0.
        fluxub(n1-1) = 0.
        e(n2) = 0.
        f(n2) = 0.
        g(n2) = 0.
        flux(n2) = 0.
        ubsp  (n2) = 0.
        ubsm  (n2) = 0.
        fluxub(n2) = 0.
#endif

c          Call iceflow to do E-W flow

        call iceflow (hold1, hnew1, facfloat,
     *                e, f, g, flux, ubsp, ubsm, fluxub, pconv1, 
     *                darea1, n1, n2, nx, dt, j, 0)

c          Store tendency (dh/dt) due to E-W flow only.
c          But if ny=1 (no N-S flow), won't do N-S loop below, 
c          so store final result in hnew.

        if (ny.gt.1) then
          do i=n1,n2
            pconv(i,j) = (hnew1(i)-hold1(i)) / dt
          enddo
        else
          do i=n1,n2
            hnew(i,j) = hnew1(i)
          enddo
        endif

c************
        endif   ! n2 > n1
c************
c>>>>>>>>>>
      enddo     ! j loop
c>>>>>>>>>>
c~~~~~~~~~~
      endif     ! nx > 1
c~~~~~~~~~~

c        Second ADI pass, doing flow for all N-S strips. 
c        Use tendency (pconv) saved from first ADI pass.
c        Set new ice thickness (hnew).

c~~~~~~~~~~~~~~~~~~~~~~
      if (ny.gt.1) then
c~~~~~~~~~~~~~~~~~~~~~~

#if defined (OMPSTEP)
!$OMP PARALLEL DO
!$OMP* shared (hold, h, hnew, hw, maskwater,
!$OMP*         hu,hv, dhsu,dhsv, dfu,dfv, lhsu,lhsv,
!$OMP*         s1a0,s1b0,s2a0,s2b0, siu,siv,
!$OMP*         via, vit, vb, va, dt, pconv, 
!$OMP*         vschoof, fschoofv,
!$OMP*         darea, dxu, dxv, dyu, dyv,
!$OMP*         zpowi, zpowim1, zpowb, zpowbm1)
!$OMP* private (hold1, hnew1, darea1, pconv1, facfloat,
!$OMP*          e, f, g, flux, ubsp, ubsm, fluxub,
!$OMP*          i, j, ig, ip, jg, jp, n1, n2, 
!$OMP*          zdfu, zdfv, zdf2, zdhsu, zdhsv, zhv,
!$OMP*          zs1a, zs1b, zs2a, zs2b,
!$OMP*          zbra, zf, zdfdh, zdfda1, zdfda2, zdh)
!$OMP* schedule(dynamic,2)
#endif

c<<<<<<<<<<<<<<
#if defined (NESTING)
      do i=2,nx-1
#else
      do i=1,nx
#endif
c<<<<<<<<<<<<<<

c        Find low-index edge of ice sheet (n1 is adjacent box with h=0).
c        If non-zero ice exists in first box, set n1=1, and b.c. will
c        implicitly be no flow out of bottom edge of box.

        n1 = -1
        do j = 1,ny
          if (hold(i,j) + dt*max(pconv(i,j),0.).gt.0.) then
            n1 = max (j-1, 1)
            goto 30
          endif
        enddo
   30   if (n1.eq.-1) n1 = 1

c        Find high-index edge of ice sheet(n2 is adjacent box with h=0).
c        If non-zero ice exists in last box, set n2=ny, and b.c. will
c        implicitly be no flow out of top edge of box.

        n2 = -1
        do j = ny,1,-1
          if (hold(i,j) + dt*max(pconv(i,j),0.).gt.0.) then
            n2 = min (j+1, ny)
            goto 40
          endif
        enddo
   40   if (n2.eq.-1) n2 = 1
#if defined (NESTING)
        n1 = max (n1, 2)
        n2 = min (n2, ny-1)
#endif

        do j=1,ny
          hold1(j) = hold(i,j)
          pconv1(j)= pconv(i,j)
          darea1(j)= darea(i,j)
          facfloat(j) = 1.
          if (maskwater(i,j).ne.0) facfloat(j) = 1.-rhor
        enddo

        call zero (e, n1d)
        call zero (f, n1d)
        call zero (g, n1d)
        call zero (flux, n1d)
        call zero (ubsp, n1d)
        call zero (ubsm, n1d)
        call zero (fluxub, n1d)

c*************************
        if (n2.gt.n1) then
c*************************

c          Set 1-D strips (N-S)

        do j=n1-1,n2
          jg = max(j,1)
          jp = min(j+1,ny)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (maskv(i,j).eq.1) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c              Use quantities at beginning of timestep for time-implicit
c              tridiag calcs (hv0, dhsv0, dfv0 passed to icestep):

c           longitudinal-stress force (lhs), strain softening (beta2),
c           set to zero currently:

            zdfv   = dfv(i,j) - lhsv(i,j)
            zdfu   = 0.25 
     *             * (dfu(i,j) + dfu(i,j+1) + dfu(i-1,j) + dfu(i-1,j+1) 
     *              -lhsu(i,j) -lhsu(i,j+1) -lhsu(i-1,j) -lhsu(i-1,j+1))
            zdf2   = max (zdfu**2 + zdfv**2, .001)
            zdhsv = dhsv(i,j)
            zdhsu = 0.25 * (  dhsu(i,j)   + dhsu(i,j+1)
     *                      + dhsu(i-1,j) + dhsu(i-1,j+1) )
            zhv   = hv(i,j)
            zdh   = (facfloat(jp)*hold(i,jp) - facfloat(jg)*hold(i,jg))
     *            / dyv(i,j)

c              Internal deformation (see notes 2/9/11):

            zs1a = 0.5*(s1a0(i,jg) + s1a0(i,jp))
            zs1b = 0.5*(s1b0(i,jg) + s1b0(i,jp))
            zs2a = 0.5*(s2a0(i,jg) + s2a0(i,jp))
            zs2b = 0.5*(s2b0(i,jg) + s2b0(i,jp))

            zbra   = zs2a*(zdf2**zpowi) + zs2b*siv(i,j)
            zf     = 2. * zbra * zdfv * zhv**2
            zdfdh  = 4. * zbra * zdfv * zhv
            zdfda1 = 2.*zbra*zhv**2
     *             + 4.*zs2a*zpowi*(zdfv**2)  *(zhv**2)*(zdf2**zpowim1)
            zdfda2 = 4.*zs2a*zpowi*(zdfv*zdfu)*(zhv**2)*(zdf2**zpowim1)

            e(j)   = (1.-fschoofv(i,j))
     *                * ( - zdfdh*zhv
     *                    - zdfda1*rhoice*grav*zhv*(zdhsv+zdh)
     *                    - zdfda2*rhoice*grav*zhv*zdhsu
     *                  ) * dxv(i,j)
            f(j)    = (1.-fschoofv(i,j))
     *                * (   zdfdh
     *                    + zdfda1*rhoice*grav*zdhsv
     *                    + zdfda2*rhoice*grav*zdhsu
     *                  ) *dxv(i,j)
            g(j)    = (1.-fschoofv(i,j))
     *                * ( zdfda1*rhoice*grav*zhv/dyv(i,j)
     *                  ) * dxv(i,j)
            flux(j) = (1.-fschoofv(i,j))
     *                * ( -zf
     *                  ) * dxv(i,j)

c           set via,vit, used in icedyn and setcoefs (softening)
            zvia = (1.-fschoofv(i,j))*(-zf)/max(zhv,1.e-20)
#if defined (UIACROSS)
            via(i,j) = zvia
#endif
#if ! defined (NOSOFTCROSS)
            vit(i,j) = -2.*zdfv*zhv * (zdf2*zs1a + siv(i,j)*zs1b)
#endif

c              Plus basal sliding or elliptic (see notes 2/9/11):

c------------------------------------
c           sheet pts, not shelf:
            if (mvind(i,j).eq.0) then
c------------------------------------

              zcrh = crhv(i,j)
              zbra = zcrh * (zdf2**zpowb)
              zf   = zbra * zdfv * zhv
              zdfdh = zbra * zdfv
              zdfda1 = zbra*zhv
     *               + 2.*zcrh*zpowb*(zdfv**2)  *zhv*(zdf2**zpowbm1)
              zdfda2 = 2.*zcrh*zpowb*(zdfv*zdfu)*zhv*(zdf2**zpowbm1)

              e(j)    = e(j)
     *                + (1.-fschoofv(i,j))
     *                  * ( - zdfdh*zhv
     *                      - zdfda1*rhoice*grav*zhv*(zdhsv+zdh)
     *                      - zdfda2*rhoice*grav*zhv*zdhsu
     *                    ) * dxv(i,j)
              f(j)    = f(j)
     *                + (1.-fschoofv(i,j))
     *                  * (   zdfdh
     *                      + zdfda1*rhoice*grav*zdhsv
     *                      + zdfda2*rhoice*grav*zdhsu
     *                    ) *dxv(i,j)
              g(j)    = g(j)
     *                + (1.-fschoofv(i,j))
     *                  * ( zdfda1*rhoice*grav*zhv/dyv(i,j)
     *                    ) * dxv(i,j)
              flux(j) = flux(j)
     *                + (1.-fschoofv(i,j))
     *                  * ( -zf
     *                    ) * dxv(i,j)

              fluxub(j) = fschoofv(i,j)
     *                  * (  max(vschoof(i,j),0.)*hold(i,jg)
     *                     + min(vschoof(i,j),0.)*hold(i,jp)
     *                    ) * dxv(i,j)
              ubsp(j) = fschoofv(i,j) * max(vschoof(i,j),0.) * dxv(i,j)
              ubsm(j) = fschoofv(i,j) * min(vschoof(i,j),0.) * dxv(i,j)

c             diagnostic:
              vb(i,j) = (1.-fschoofv(i,j)) *(-zf)/max(hv(i,j),1.e-20)
     *                +     fschoofv(i,j)  * vschoof(i,j)
#if defined (COP21RD)
              va(i,j) = via(i,j) + vb(i,j)
#else
              va(i,j) = zvia + vb(i,j)
#endif

c-----------------------------------------
c           elliptic pts:
            else if (mvind(i,j).ge.1) then
c-----------------------------------------

c             nb: if fschoofv > 0, Schoof velocities are already 
c             imposed in elliptic soln va.
              vb(i,j) = va(i,j) - via(i,j)
              ubsp(j) = max(vb(i,j),0.)*dxv(i,j)
              ubsm(j) = min(vb(i,j),0.)*dxv(i,j)
              fluxub(j) = ubsp(j)*hold(i,jg) + ubsm(j)*hold(i,jp)

c----------------
            endif
c----------------

#if defined (SUBGRIDEDGE) || ! defined (NOSUBGRIDEDGESHELF)     
c           no flow *out* of boxes with fedge < 1 (medgev set above):
            if (medgev(i,j).eq.0) then
              e(j) = 0.
              f(j) = 0.
              g(j) = 0.
              f(j) = 0.
              flux(j) = 0.
              fluxub(j) = 0.
              ubsp(j) = 0.
              ubsm(j) = 0.
              va(i,j) = 0.
              vb(i,j) = 0.
              via(i,j) = 0.
              vit(i,j) = 0.
            endif
#endif

c~~~~~~~~~~~~~~
          endif
c~~~~~~~~~~~~~~

        enddo

c Don't do following if nesting (bug fix 5/25/10)
#if ! defined (NESTING)
        e(n1-1) = 0.
        f(n1-1) = 0.
        g(n1-1) = 0.
        flux(n1-1) = 0.
        ubsp  (n1-1) = 0.
        ubsm  (n1-1) = 0.
        fluxub(n1-1) = 0.
        e(n2) = 0.
        f(n2) = 0.
        g(n2) = 0.
        flux(n2) = 0.
        ubsp  (n2) = 0.
        ubsm  (n2) = 0.
        fluxub(n2) = 0.
#endif

c          Call iceflow to do N-S flow

        call iceflow (hold1, hnew1, facfloat,
     *                e, f, g, flux, ubsp, ubsm, fluxub, pconv1,
     *                darea1, n1, n2, ny, dt, j, 1)

c          Store new 2-D ice thickness

        do j=n1,n2
          hnew(i,j) = hnew1(j)
        enddo

c************
        endif   ! n2 > n1
c************
c<<<<<<<<<<
      enddo     ! i loop
c<<<<<<<<<<
c~~~~~~~~~~
      endif     ! ny > 1
c~~~~~~~~~~

      return
      end

c-----------------------------------------------------------------------

      subroutine iceflow (h0, h1, facfloat,
     *                    e, f, g, flux, ubsp, ubsm, fluxub, yconv, 
     *                    dar, n1, n2, nxory, dt, j, ifv)

c        Steps ice-sheet 1-D flow equation through one time step. 
c        Time-implicit contributions lead to a tridiagonal matrix 
c        equation. Units are meters and years.

c        The current direction is called the "x-direction" (index i)
c        here, and the perpendicular direction is the "y-direction",
c        accounted for either by yconv.
c        Only the x-direction has time-implicit terms.
c        The correspondence with real E-W or N-S directions depends
c        on the current ADI stage in the calling program.

c        Ice-sheet thickness are defined on h-grid (h).
c        Other terms (...u) are defined on a staggered u-grid, with 
c        ...u(i) at the interface between h-grid boxes i and i+1.

c        h0   = actual ice thickness at start of time step (m)(supplied)
c        h1   = actual ice thickness at end of time step (m) (returned)
c        yconv= convergence rate due to perpendicular flow (m/y) (supp)
c        dar  = grid-box area (h-grid) (m**2) (supplied)
c        dt   = timestep (yr) (supplied)

#include <comicephys.h>
#include <comicegrid.h>

      parameter (n1d=nx+ny)
      dimension 
     *  h0(n1d), h1(n1d), facfloat(n1d),
     *  e(0:n1d), f(0:n1d), g(0:n1d), flux(0:n1d), 
     *  ubsp(0:n1d), ubsm(0:n1d), fluxub(0:n1d), yconv(n1d), 
     *  dar(n1d)

      dimension
     *  a(n1d), b(n1d), c(n1d), forc(n1d)

      parameter (rimp = 1.)         ! time-implicit vs explicit fraction
c     parameter (rimp = 0.5)        ! time-implicit vs explicit fraction

c        Set implicit vs explicit timesteps

      dti =     rimp *dt
      dte = (1.-rimp)*dt

c       Set the 3 diagonals of the tridiagonal matrix, and the rhs

      do i = n1,n2
        im = max(i-1,1)
        ip = min(i+1,nxory)

c       main diagonal:
        a(i) = 1.+ dti*(   (-0.5*f(i)   + g(i)*facfloat(i)) 
     *                   + ( 0.5*f(i-1) + g(i-1)*facfloat(i))
     *                   + ( ubsp (i)   - ubsm (i-1))    ! upstream
c    *                   + ( 0.5*( -ubsp(i-1) + ubsp(i)  ! centered
c    *                             -ubsm(i-1) + ubsm(i)))! centered
     *                 ) / dar(i)
       
c       superdiagonal:
        b(i) =     dti*(   (-0.5*f(i)   - g(i)*facfloat(ip))
     *                              + ubsm (i)                ! upstream
c    *                              + 0.5*(ubsp(i) + ubsm(i)) ! centered
     *                 ) / dar(i)

c       subdiagonal:
        c(i) =     dti*(   (0.5*f(i-1) - g(i-1)*facfloat(im))
     *                              - ubsp (i-1)               !upstream
c    *                              - 0.5*(ubsp(i-1)+ubsm(i-1))!centered
     *                 ) / dar(i)

c       rhs:
        forc(i)= h0(i)
     *           + dti*(  e(i)      - e(i-1)      ) / dar(i)
     *           + dt *( -flux(i)   + flux(i-1)   ) / dar(i)
     *           + dte*( -fluxub(i) + fluxub(i-1) ) / dar(i)
     *           + dt * yconv(i)

#if defined (NESTING)
        if (i.eq.n1 .and. n1.eq.2) then
c         if at domain edge (2), add terms for time-invariant h(1) b.c.,
c         associated with time-implicit terms in a(2), b(2):
          forc(i) = forc(i) - h0(i-1)*c(i)

        else if (i.eq.n2 .and. n2.eq.nxory-1) then
c         Similarly for other edge (nxory-1), for fixed h(nxory),
c         associated with a(nxory-1), b(nxory-1):
          forc(i) = forc(i) - h0(i+1)*b(i)
        endif
#endif

      enddo

#if defined (EISLINE) && defined (LINEE)
c     ice stream entering from left, prescribed l.h. thickness:
      a(n1) = 1.
      b(n1) = 0.
      forc(n1) = h0(n1)
      a(nx) = 1.
      c(nx) = 0.  ! was b(nx) = 0...error?
      forc(nx) = h0(nx)
#endif

c        Solve tridiagonal system (steps ice thickness from h0 to h1)

      call tridia_i (a(n1), b(n1), c(n1), h1(n1), forc(n1), n2-n1+1, 0)

      return
      end

c-----------------------------------------------------------------------

      subroutine tridia_i (a, b, c, h, forc, n, ktrid)

c        Solves tridiagonal matrix equation. 
c        a,b,c are diagonals of the matrix.
c        a is main diagonal (1:n)
c        b is superdiagonal (1:n-1),
c        c is subdiagonal (2:n),
c        forc is right-hand side, 
c        h is vector to be solved for.
c        n is the size of the matrix and vectors. 
c        If ktrid is non-zero, uses condition that resulting h(i) 
c        must be >=  0.

      dimension a(n), b(n), c(n), h(n), forc(n)
c     parameter (nmax = 25000)
      parameter (nmax = 50000)
      dimension r(nmax), s(nmax), y(nmax)

c        Check that local arrays are big enough

      if (n.gt.nmax) then
        write(6,*)'*** Error in tridia_i: Too many points.',
     *            '   nmax=',nmax,'   n=',n
        stop
      endif

c     amin = 1.e-3
      amin = 1.e-6
      do i = 1,n
        if (abs(a(i)).lt.amin) then
          write (6,9000) i,a(i)
 9000     format(' tridia_i warning: a(',i3,')=',f12.8)
          if (a(i).ge.0.) a(i) = amin
          if (a(i).lt.0.) a(i) = -amin
        endif
      enddo

      s(1) = a(1)
      r(1) = b(1)/s(1)
      if (n.gt.2) then
        do i = 2,n-1
          s(i) = a(i)-c(i)*r(i-1)
          r(i) = b(i)/s(i)
        enddo
      endif
      if (n.ge.2) then
        s(n) = a(n)-c(n)*r(n-1)
      endif

      y(1) = forc(1)/s(1)
      if (n.ge.2) then
        do i = 2,n
          y(i) = (forc(i)-c(i)*y(i-1)) / s(i)
        enddo
      endif

      if (ktrid.eq.0) then
        h(n) = y(n)
        if (n.ge.2) then
          do i = n-1,1,-1
            h(i) = y(i)-r(i)*h(i+1)
          enddo
        endif
      else
        h(n) = max (0., y(n))
        if (n.ge.2) then
          do i = n-1,1,-1
            h(i) = max (0., y(i)-r(i)*h(i+1))
          enddo
        endif
      endif

      return
      end

c-----------------------------------------------------------------------

c====================
#if defined (MELANGE)
c====================

      subroutine setbpumel (bpumel, bpvmel, thetaumel, thetavmel,
     *                      hmel, h, hedge, maskwater, fedge)

c        Set back pressure by melange on regular-calving ice faces 
c        (bp[u,v]mel, on u and v grids), for bc in ssa equations
c        in icedyn (dotria's, dosparse's).
c        Applied at u,v-grid interfaces in dosparse, at edges of ice 
c        shelf, with one h-cell having fedge=0 and the other h,fedge>0.
c        Melange interacts with ice (usually shelf) either:
c        (1) in the same  h-grid cell if 0<fedge<1 as it flows to 
c             adjacent cell with fedge=0, or:
c        (2) with adjacent cell if this cell has fedge=0 and adjacent 
c            cell has fedge=1. 

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  bpumel(0:nxp,0:nyp),  bpvmel(0:nxp,0:nyp),
     *  thetaumel(nx,ny),     thetavmel(nx,ny),
     *  hmel(nx,ny),          h(nx,ny),            hedge(nx,ny),
     *  maskwater(nx,ny),     fedge(nx,ny) 

#include <comicepresmel.h>

      call resetr (bpumel, (nxp+1)*(nyp+1), 1.)
      call resetr (bpvmel, (nxp+1)*(nyp+1), 1.)

c        bpumel, u-grid:

      do j = 1,ny
        ja = j
        jb = j
        do i = 1,nx-1
          do iloop=1,2
            if (iloop.eq.1) then
              ia = i
              ib = i+1
            else
              ia = i+1
              ib = i
            endif
               
            if (h(ia,ja).gt.0. .and. fedge(ib,jb).eq.0.) then 
              if (fedge(ia,ja).eq.1.) then
                ic = ib
                jc = jb
              else 
                ic = ia
                jc = ja
              endif
              zr = 
     *          (0.5*rhomp*grav*(hmel(ic,jc)**2) + presmel(hmel(ic,jc)))
     *        / (0.5*rhoip*grav*(hedge(ia,ja)**2))
              bpumel(i,j) = 1. - min(zr,1.) * (1.-thetaumel(ic,jc))
              bpumel(i,j) = max (0., min (1., bpumel(i,j)))
            endif
          enddo
        enddo
      enddo

c        bpvmel, v-grid:

      do i = 1,nx
        ia = i
        ib = i
        do j = 1,ny-1
          do iloop=1,2
            if (iloop.eq.1) then
              ja = j
              jb = j+1
            else
              ja = j+1
              jb = j
            endif
               
            if (h(ia,ja).gt.0. .and. fedge(ib,jb).eq.0.) then 
              if (fedge(ia,ja).eq.1.) then
                ic = ib
                jc = jb
              else 
                ic = ia
                jc = ja
              endif
              zr = 
     *          (0.5*rhomp*grav*(hmel(ic,jc)**2) + presmel(hmel(ic,jc)))
     *        / (0.5*rhoip*grav*(hedge(ia,ja)**2))
              bpvmel(i,j) = 1. - min(zr,1.) * (1.-thetavmel(ic,jc))
              bpvmel(i,j) = max (0., min (1., bpvmel(i,j)))
            endif
          enddo
        enddo
      enddo

      return
      end

c=====
#endif
c=====

c-----------------------------------------------------------------------
