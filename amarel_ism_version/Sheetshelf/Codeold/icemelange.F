c||||||||||||||||||||
#if defined (MELANGE)
c||||||||||||||||||||

c-----------------------------------------------------------------------

      subroutine icemel (h, hs, hw, fmel, maskmel,
     *                   budgallmel,  budgsupmel,
     *                   budgnegmel,  budgadjmel, 
     *                   budgredmel,  budgsmomel,
     *                   curbackmel,  curthetmel,
     *                   curhmel,     curlmel,
     *                   ua, va, ucimel, vcimel, hucimel, hvcimel, 
     *                   thetau, thetav, 
     *                   uabound, vabound, 
c                        above are *mel in main prog. below are not:
     *                   hb, hbsd, hice, hedge, fedge, 
     *                   maskwater, tsurfi, uaice, vaice, calvice, 
     *                   sealev, timeice, dt, ifrest, iloopmel,nloopmel,
     *                   iffirse, iffirsr, nyeartabmel)

c     Steps combined SIA-SSA for melange dynamics (like icedyn.F).
c     All variable names above are in main program with *mel added.    
c     All subr names are as in icedyn.F but with *mel added.
c
c     Shelf dyanimcs (SSA) solves for ua,va (dotriamel, dosparsemel). 
c     Melange thickness is advected by dosparsemelh.
c     ua,va, are passed, saved (in main prog) between iterations.
c 
c     No prognostic melange temperature, but arhap can depend
c     on tsurfi.
c
c     Unlike ice dynamics, no melange thicknesses or velocities
c     are on nest-driving files, so for nested runs, apply same
c     b.c.s for melange as continental runs - no melange in grid
c     boxes touching domain edges. So all "if defined (NESTING)"
c     lines in icedyn.F are removed here, and "if ! defined" are kept.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

c     passed:
      dimension 
     *  h(nx,ny),              hs(nx,ny),          hw(nx,ny),
     *  fmel(nx,ny),           maskmel(nx,ny),
     *  budgallmel(nx,ny),     budgsupmel(nx,ny), 
     *  budgnegmel(nx,ny),     budgadjmel(nx,ny),  
     *  budgredmel(nx,ny),     budgsmomel(nx,ny),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  ucimel(0:nxp,0:nyp),   vcimel(0:nxp,0:nyp),
     *  hucimel(0:nxp,0:nyp),  hvcimel(0:nxp,0:nyp),
     *  thetau(nx,ny),         thetav(nx,ny),

     *  hb(nx,ny),             hbsd(nx,ny),    
     *  hice(nx,ny),           hedge(nx,ny),      fedge(nx,ny),
     *  maskwater(nx,ny),      tsurfi(nx,ny),
     *  uaice(0:nxp,0:nyp),    vaice(0:nxp,0:nyp),
     *  calvice(nx,ny)

c     local:
      dimension
     *  maskh(nx,ny),          visc(nx,ny),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp), 
     *  crhu(0:nxp,0:nyp),     crhv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  h0(nx,ny),             h1(nx,ny),
     *  coef1(nx,ny),
     *  dhsu(0:nxp,0:nyp),     dhsv(0:nxp,0:nyp),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp)

      dimension
     *  arhap(nx,ny)
      save arhap

c     cross grid:
      dimension
     *  coef2(0:nx,0:ny),      maskc(0:nx,0:ny)

c     for iteration (a,c) convergence criteria:
      dimension
     *  uaprev (0:nxp,0:nyp),  vaprev (0:nxp,0:nyp),     ! loop a
     *  h1prev(nx,ny)                                    ! loop c

c     for budget tracking and melange supply:
      dimension 
     *  ucamel(0:nxp,0:nyp),  vcamel(0:nxp,0:nyp),
     *  hucamel(0:nxp,0:nyp), hvcamel(0:nxp,0:nyp),
     *  supcalv(nx,ny), 
     *  budgallmelsav(nx,ny),
     *  budgallmeltmp(nx,ny),  
     *  budgsupmeltmp(nx,ny),  
     *  budgnegmeltmp(nx,ny),
     *  budgsmomeltmp(nx,ny)
      save budgallmelsav

c     for melange veloc bc's at calving/cliff ice faces:
      dimension
     *  uabound(0:nxp,0:nyp),  vabound(0:nxp,0:nyp)

      dimension hsmoprev(nx,ny)

      parameter (delcrita = .01)   ! m/a  ! .01,.03 
      parameter (delcritc = .01)   ! m    ! .01,.03

#if defined (NITERAMEL)
      parameter (niteramel = NITERAMEL)
#else
      parameter (niteramel = 4) 
#endif

#if defined (NITERCMEL)
      parameter (nitercmel = NITERCMEL)
#else
      parameter (nitercmel = 2)
#endif

#if defined (TESTMEL2DDIAG3)
      dimension zz1(nx), zz2(nx), zz3(nx), zz4(nx) 
#endif

      parameter (nlook = 4)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /

c     for arhap, as in icetherm (arrhenius):
      parameter (qact1    = 9.545e4,                       ! warm, J/mol
     *           qact2    = 7.820e4,                       ! cold, J/mol
     *           tact     = -6.5,                          ! deg K
     *           rgas     = 8.314,                         ! J/mol/K
     *           zq1  = qact1/rgas,  zq2 = qact2/rgas)

#include <comicepresmel.h>

c-------------------------
c Start of executable code
c-------------------------

c        Like rheol coeff settings for arhap in icetherm (arrhenius)

#if defined (ENHANCESHELFMEL)
      enhanceshelfmel = ENHANCESHELFMEL
#else
c     enhanceshelfmel = 0.3                  ! as for shelf in arrhenius
      enhanceshelfmel = 3.0
#endif

      if (powimel.eq.1) then
        crheoli1 = 0.6e-8                  ! warm, Pa-1 a-1
        crheoli2 = 0.5e-8                  ! cold, Pa-1 a-1
      else if (powimel.eq.2) then 
        crheoli1 = 0.6e-11                 ! warm, Pa-2 a-1
        crheoli2 = 0.5e-11                 ! cold, Pa-2 a-1
      else if (powimel.eq.3) then 
c       not used, see below (as in arrhenius for EISANTA):
        crheoli1 = 2.0 e-16                ! warm, Pa-3 a-1
        crheoli2 = 1.66e-16                ! cold, Pa-3 a-1
      else
        crheoli1 = 0.6e-8/((1.e4)**(powimel-1))
        crheoli2 = crheoli1
      endif

      do j=1,ny
        do i=1,nx
          ztstar = tsurfi(i,j)                   ! taken as melange temp
c         ztstar = tmelt

          if (powimel.eq.3) then 
            if (ztstar. ge. 263.15) then
              za = 5.47e10 * exp(-13.9e4/(8.314*ztstar)) !warm, Pa-3 a-1
            else
              za = 1.14e-5 * exp(-6.0 e4/(8.314*ztstar)) !cold, Pa-3 a-1
            endif
          else
            if (ztstar-tmelt. ge. tact) then
              za = crheoli1*exp(zq1*(1./tmelt - 1./ztstar))
            else
              za = crheoli2*exp(zq2*(1./tmelt - 1./ztstar))
            endif
          endif

          arhap(i,j) = (enhanceshelfmel*za)**powirmel 
c         arhap(i,j) = enhanceshelfmel*(za**powirmel)
        enddo
      enddo

c        Like basal coeff settings for crh[u,v] in icehydrol (basecoef)

#if defined (CRHUMEL)
      zcrh = CRHUMEL                  ! (m/y)/Pa^powbmel
#else
      if (powbmel.eq.1) zcrh = 1.e-1  ! (m/y)/Pa
      if (powbmel.eq.2) zcrh = 1.e-5  ! (m/y)/Pa^2
      if (powbmel.eq.3) zcrh = 1.e-9  ! (m/y)/Pa^3
#endif

      do j=0,ny+1
        do i=0,nx+1
          crhu(i,j) = zcrh
          crhv(i,j) = zcrh
        enddo   
      enddo   

c        For budget tracking - nb: icemelange called multiple times
c        (nloopmel) for every main model timestep (dtimeice)

      if (iloopmel.eq.1) then
        call scopy (nx*ny, budgallmel, 1, budgallmelsav, 1)
        call zero (budgallmel, nx*ny)
        call zero (budgsupmel, nx*ny)
        call zero (budgnegmel, nx*ny)
        call zero (budgadjmel, nx*ny)
        call zero (budgredmel, nx*ny)
        call zero (budgsmomel, nx*ny)
      endif
       
c       Set basal coefficient cutoff to select points for inclusion
c       in elliptic SSA (shelf) eqns, or sheet only (m[u,v]ind below)

c       Set |ua|,|va| speed limit (used after dosparse)

#if defined (SPEEDLIMITMEL)
      if (SPEEDLIMITMEL.gt.0.) then
        vmax = SPEEDLIMITMEL                      ! 20.e3 ! 10.e3 ! 5.e3
      else if (SPEEDLIMITMEL.eq.0.) then
        vmax = 0.75 * dd0 / dt
      else
        vmax = 1.e20
      endif
#else
      vmax = 0.75 * dd0 / dt
#endif

c============================
c     Top of loop C (~Picard).                ! Moved here from below to
      do iterc = 1,nitercmel                  ! set mask[h,u] for safety
c============================

c        Adjust melange surface elevations, water depths, maskmel, etc.
c        Also zero melange for maskwater=0 or fedge=1 (same after first
c        iloopmel call). For first iloopmel call only, redistribute  
c        zero-ed melange into neighboring eligible cells if any,
c        or if not, track in budgadjmel.

      call adjustmel (maskmel, h, hs, hw, fmel,
     *                hb, hice, fedge, maskwater, sealev, 
     *                budgadjmel, dt, iloopmel.eq.1.and.iterc.eq.1)

c        Impose boundary constraints on h for various expts 

c....................................................................
#if (defined (EISANTA) && defined (TRANSECTA)) || defined (TESTMEL2D)
c....................................................................
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

c.....
#endif
c.....

c        Set masks (maskh, masku, maskv). All masks are set 
c        outside iteration loops, ignoring changes of non-zero <-> 
c        zero ice (h) between iterations. 

c        Set maskh (0/1 if no ice/ice, h grid)

      call izero (maskh, nx*ny)
      do j = 1,ny
        do i = 1,nx
          if (h(i,j).gt.0.) then     ! h can only be > 0 if maskmel != 0
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
     *              hs(i,j).gt.hs(i+1,j) .and. maskmel(i+1,j).ne.0 )then
            masku(i,j) = 1
          else if ( maskh(i,j).eq.0 .and. maskh(i+1,j).eq.1 .and. 
     *              hs(i,j).lt.hs(i+1,j) .and. maskmel(i,j).ne.0) then
            masku(i,j) = 1
          endif
        enddo
      enddo

      call izero (maskv, (nxp+1)*(nyp+1))
      do j = 1,ny-1
        do i = 1,nx
          if ( maskh(i,j).eq.1 .and. maskh(i,j+1).eq.1 ) then
            maskv(i,j) = 1
          else if ( maskh(i,j).eq.1 .and. maskh(i,j+1).eq.0 .and. 
     *              hs(i,j).gt.hs(i,j+1) .and. maskmel(i,j+1).ne.0)then
            maskv(i,j) = 1
          else if (maskh(i,j).eq.0 .and. maskh(i,j+1).eq.1 .and. 
     *             hs(i,j).lt.hs(i,j+1) .and. maskmel(i,j).ne.0) then
            maskv(i,j) = 1
          endif
        enddo
      enddo

c        Zero diagnostic quantities (some necessary for finite diffs 
c        near ice edges)

      call zero (hu,   (nxp+1)*(nyp+1))
      call zero (hv,   (nxp+1)*(nyp+1))
      call zero (dhsu, (nxp+1)*(nyp+1))
      call zero (dhsv, (nxp+1)*(nyp+1))

      call zero (dfu,  (nxp+1)*(nyp+1))
      call zero (dfv,  (nxp+1)*(nyp+1))

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.0) then
            ua(i,j)  = 0.
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.0) then
            va(i,j)  = 0.
          endif
        enddo
      enddo

c        If first timestep of run, zero all velocities. Else, keep 
c        for next timestep's initial guess (u,v grids, all passed)
c        (some necessary for finite diffs near ice edges)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (iffirsr.eq.1 .and. iloopmel.eq.1) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       ua,va are not on restart file:
        call zero (ua, (nxp+1)*(nyp+1))
        call zero (va, (nxp+1)*(nyp+1))
c~~~~~~~~~~
      endif
c~~~~~~~~~~

c        Test for no ice. If so, just apply local forcing (including
c        bc flow from ice shelves or sheet) and skip all velocity calcs.
c        Supply from cliff failure is [h][u,v]icemel (u,v interfaces).
c        Supply from calving is calced in calv2mel (supcalv).

      numh = 0
      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            numh = numh + 1
          endif
        enddo
      enddo
c . . . . . . . . . . . .
      if (numh.eq.0) then
c . . . . . . . . . . . .
        call calv2mel (fedge, hedge, hice, calvice, maskwater, maskmel,
     *                 ucamel, vcamel, hucamel, hvcamel, supcalv)

        do j=1,ny
          do i=1,nx
            if (maskmel(i,j).ne.0) then
              zsup = (  (ucimel(i-1,j)*hucimel(i-1,j)) * dyu(i-1,j)
     *                + (ucimel(i,j)  *hucimel(i,j)  ) * dyu(i,j)
     *                + (vcimel(i,j-1)*hvcimel(i,j-1)) * dxu(i,j-1)
     *                + (vcimel(i,j)  *hvcimel(i,j)  ) * dxu(i,j)
     *               ) / darea(i,j)
     *               + supcalv(i,j)
    
#  if defined (TESTMEL2D) && defined (TESTMELPROF)
c             For 1-D melange flowline, set supply in 2nd grid pt:
              if (i.eq.2) then
c               zsup = 0.25 * 1.e4*1.e3*dyu(i-1,j) / darea(i,j) ! Jakob
                zsup = hice(1,j)*uaice(1,j)*dyu(1,j)/darea(2,j)
                supcalv(2,j) = zsup 
              else
                zsup = 0.
                supcalv(i,j) = 0.
              endif
#  endif

              h(i,j) = zsup * dt
              budgsupmel(i,j) = budgsupmel(i,j) + zsup            ! diag

              if (h(i,j).gt.0.) then
                h(i,j) = h(i,j) + budgallmelsav(i,j)*dt
                budgallmel(i,j) = budgallmel(i,j)+budgallmelsav(i,j) !d.
              endif

              if (h(i,j).lt.0.) then
                budgnegmel(i,j) = budgnegmel(i,j) - h(i,j)/dt     ! diag
                h(i,j) = 0.
              endif
            endif
          enddo   
        enddo   

        goto 3500
c . . . . . 
      endif
c . . . . . 

c============================
c     Top of loop C (~Picard)  ! moved further up for safety (mask[h,u])
c     do iterc = 1,nitercmel  
c============================

c        Adjust melange surface elevations, water depths, maskmel, etc 
c        (for change in h)

      call adjustmel (maskmel, h, hs, hw, fmel,
     *                hb, hice, fedge, maskwater, sealev,
     *                budgadjmel, dt, .false.)

c       Set hu, hv, dhsu, dhsv, dfu, dfv  (u,v grids)
c       Do inside loop C to capture hs changes within iters.

      do j = 1,ny
        do i = 1,nx-1
          if (masku(i,j).eq.1) then
            hu(i,j) = 0.5*(h(i,j) + h(i+1,j))
            dhsu(i,j) = (hs(i+1,j) - hs(i,j)) / dxu(i,j)
            dfu(i,j)  = rhomel*grav*hu(i,j)*dhsu(i,j)
     *                + (presmel(h(i+1,j))-presmel(h(i,j))) / dxu(i,j)
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then
            hv(i,j) = 0.5*(h(i,j) + h(i,j+1))
            dhsv(i,j) = (hs(i,j+1) - hs(i,j)) / dyv(i,j)
            dfv(i,j)  = rhomel*grav*hv(i,j)*dhsv(i,j)
     *                + (presmel(h(i,j+1))-presmel(h(i,j))) / dyv(i,j)
          endif
        enddo
      enddo

c     store original h (for whole time-step) in h0, etc,
c     (and initialize h1):
      if (iterc.eq.1) then 
        call scopy (nx*ny, h, 1, h0, 1)
        call scopy (nx*ny, h, 1, h1, 1)
      endif

      call scopy (nx*ny, h1, 1, h1prev, 1)

c....................................................................
#if (defined (EISANTA) && defined (TRANSECTA)) || defined (TESTMEL2D)
c....................................................................
      do j=1,ny
        ua(0,j) = 0.
c       ua(nx,j) = 0.
      enddo
c.....
#endif
c.....

c        Set nuvtot and muind,mvind, indices from 2-D arrays into 
c        linear sequence of selected points for tridia or sparse solns
c        below. (Can do outside A loop, depends only on mask[u,v],hw).  

      call izero (muind, (nxp+1)*(nyp+1))
      call izero (mvind, (nxp+1)*(nyp+1))
      nuvtot = 0

      do j = 1,ny
        do i = 1,nx
          if (masku(i,j).eq.1) then
            nuvtot = nuvtot + 1
            muind(i,j) = nuvtot
          endif

          if (maskv(i,j).eq.1) then
            nuvtot = nuvtot + 1
            mvind(i,j) = nuvtot
          endif
        enddo
      enddo

c==========================
      if (nuvtot.gt.0) then ! if no shelf-like pts, skip shelf calcs
c==========================

c===============================
c     Top of loop A (shelf flow)
      do itera = 1,niteramel
c===============================

      call scopy ((nxp+1)*(nyp+1), ua, 1, uaprev, 1)
      call scopy ((nxp+1)*(nyp+1), va, 1, vaprev, 1)

c        Set various coefficients (coef1, coef2, coefb[u,v])

      call setcoefsmel (coef1, coef2, visc, coefbu, coefbv,
     *                  maskh, arhap, h, hbsd, hw,
     *                  masku, maskv, ua, va, 
     *                  hu, hv, crhu, crhv,
     *                  muind, mvind, maskc, itera, iterc) 

c        Set melange veloc bc's at calving/cliff ice faces

      call calv2mel (fedge, hedge, hice, calvice, maskwater, maskmel,
     *               ucamel, vcamel, hucamel, hvcamel, supcalv)
      call setvelbc (uabound, vabound, 
     *               h, maskmel, supcalv,
     *               ucimel, vcimel, hucimel, hvcimel, 
     *               ucamel, vcamel, hucamel, hvcamel, 
     *               uaice, vaice, hice, fedge)

      do j=0,nyp
        do i=0,nxp
          if (uabound(i,j).ne.0.) ua(i,j) = uabound(i,j)
          if (vabound(i,j).ne.0.) va(i,j) = vabound(i,j)
        enddo
      enddo

c        Tridia or Sparse matrix solution of ice-shelf eqns (for ua,va)

c----------------------
      if (ny.eq.1) then
c----------------------
        call dotridiamel  (ua, 
     *                     hu, dfu, hw, coef1, coefbu,   
     *                     muind, nuvtot)
c---------
      else
c---------

#if defined (DOSPARSEMELU)

        call dosparsemelu (ua, va, 
     *                     hu, hv, dfu, dfv,
     *                     coef1, coef2, coefbu, coefbv, 
     *                     muind, mvind, nuvtot, hw, maskmel,
     *                     timeice, itera, iterc, iloopmel
     *                     )

#elif defined (DOSPARSESORMELU)

        call dosparse_sormelu
     *                (ua, va, 
     *                 hu, hv, dfu, dfv,
     *                 coef1, coef2, coefbu, coefbv, 
     *                 muind, mvind, h, hs, hb, hw, maskmel,
     *                 timeice, itera, iterc, iloopmel)

#else
        write (ioterm,'(/2a)')
     *    '*** Error (icemel): must specify -DDOSPARSEMELU or',
     *    ' -DDOSPARSESORMELU'
        stop
#endif

#if defined (SMOOTHMELU)
        if (SMOOTHMELU.gt.0) then
          do m=1,SMOOTHMELU
            call smoothermelu (ua, maskmel, 1)
            call smoothermelu (va, maskmel, 2)
          enddo
        endif
#endif

c----------
      endif
c----------

c!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined (TESTMEL2DDIAG3)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Explicit verification of solution of ua-eqn 

      call zero (zz1, nx)
      call zero (zz2, nx)
      call zero (zz3, nx)
      call zero (zz4, nx)

      j = (ny+1)/2
      jm= max (j-1,1)
      iah = 0
      do i=1,nx
        im= max (i-1,1)
        if (h(i,j).gt.0.) then
          if (iah.eq.0) iah = i
          zdudx = (ua(i,j)-ua(im,j))/ dx(i,j)
          zdvdy = (va(i,j)-va(i,jm))/ dy(i,j)
          zz1(i) = coef1(i,j)*(2.*zdudx + zdvdy)
          if (zdudx.gt.0.) zz1(i) = coefdivmel*zz1(i)
          zz2(i) = ( 0.5*rhomp*grav*(h(i,j)**2) + presmel(h(i,j)) )
        endif
      enddo

      totzz3 = 0.
      totzz4 = 0.
      do i=1,nx
        if (h(i,j).gt.0.) then
          zz3(i) = coefbu(i,j)*ua(i,j)
          zz4(i) = ua(i,j)
     *           * sidedragmel * 0.5*(0.5*(coef1(i+1,j)+coef1(i,j)))
     *           / (0.5*dyu(i,j)*dyu(i,j))
          totzz3 = totzz3 + zz3(i)*dxu(i,j)
          totzz4 = totzz4 + zz4(i)*dxu(i,j)
        endif
      enddo

      if (totzz3.gt.0. .or. totzz4.gt.0.) then
        write (888,*)
        do i=1,nx
          ip = min(i+1,nx)
          if (h(i,j).gt.0.) then
            write (888,'(f10.2,3i4, 4i4,2f10.4,f15.4, 10x,5e15.6)')
     *        timeice, iloopmel, iterc, i, 
     *        muind(i,j), nuvtot, maskmel(i,j), maskmel(ip,j), 
     *        h(i,j), h(ip,j), ua(i,j),

     *        (zz1(ip)-zz1(i))/dxu(i,j), (zz2(ip)-zz2(i))/dxu(i,j),
     *        zz3(i), zz4(i),
     *        (zz1(ip)-zz1(i))/dxu(i,j) - (zz2(ip)-zz2(i))/dxu(i,j)
     *        - zz3(i) - zz4(i)
          endif
        enddo

        write (888,'(/60x, 5e15.6)')
     *    zz1(iah), zz2(iah),
     *    totzz3, totzz4,
     *    -zz1(iah) + zz2(iah) - totzz3 - totzz4
      endif
c!!!!!
#endif
c!!!!!
       
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
c    *    'itera=',itera,'  delmaxa=',delmaxa,'  delcrita=',delcrita
c     enddo

      if (delmaxa.lt.delcrita) goto 1000

c=======================
c     End loop A (itera)
      enddo
      itera = niteramel
 1000 continue
c=======================
c==========
      endif ! nuvtot > 0 test
c==========

c        Diagnostic check of domain-wide force balance (to fort.281).
c        Must do before h changes in loop C below.

#if defined (CHECKMELFORCE)
      call checkmelforce (h, hs, maskmel,
     *                    thetau, thetav, coef1, arhap,
     *                    ua, va, muind, mvind,
     *                    coefbu, coefbv,
     *                    hice, fedge,
     *                    timeice, dt, 
     *                    iterc, iloopmel, nloopmel, nyeartabmel)
#endif

c        Explicitly apply surface mass balance to ice thickness.
c        For melange, apply only if melange exists already (h1 > 0).
c        budgallmel = budgsnow-budgevap-budgmelt-oceanmelt, set in 
c        icectl. Not multiplied by fmel - for dyanmics and supply,
c        consider all melange cells 100% covered. At end, consider 
c        displaced by any ice in same cell (1-fedge) - doesn't 
c        conserve melange ice.
c        See comments above re. supply from cliff and calving.

      call calv2mel (fedge, hedge, hice, calvice, maskwater, maskmel,
     *               ucamel, vcamel, hucamel, hvcamel, supcalv)

      call zero (budgallmeltmp, nx*ny)
      call zero (budgsupmeltmp, nx*ny)

      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).ne.0) then
            zsup = (  (ucimel(i-1,j)*hucimel(i-1,j)) * dyu(i-1,j)
     *              + (ucimel(i,j)  *hucimel(i,j)  ) * dyu(i,j)
     *              + (vcimel(i,j-1)*hvcimel(i,j-1)) * dxu(i,j-1)
     *              + (vcimel(i,j)  *hvcimel(i,j)  ) * dxu(i,j)
     *             ) / darea(i,j)
     *             + supcalv(i,j)

#  if defined (TESTMEL2D) && defined (TESTMELPROF)
c             For 1-D melange flowline, set supply in 2nd grid pt:
              if (i.eq.2) then
c               zsup = 0.25 * 1.e4*1.e3*dyu(i-1,j) / darea(i,j) ! Jakob
                zsup = hice(1,j)*uaice(1,j)*dyu(1,j)/darea(2,j)
                supcalv(2,j) = zsup 
              else
                zsup = 0.
                supcalv(i,j) = 0.
              endif
#  endif

            budgsupmeltmp(i,j) = zsup
            budgallmeltmp(i,j) = budgallmelsav(i,j)
          endif
        enddo
      enddo
       
c        "Slump" new supply of melange onto existing melange surface 

      call slumpsup (h0, maskmel, budgsupmeltmp,
     *               dt, timeice, iterc, iloopmel)

c       Advect melange thickness through timestep (h0 to h1)

cmas  call mascalc (h0, ztoth0, ztota0)

      call dosparsemelh (ua, va, 
     *                   h0, h, h1, maskmel,
     *                   budgallmeltmp, budgsupmeltmp, budgnegmeltmp,
     *                   timeice, dt, iterc, iloopmel)

cmas  call mascalc (h1, ztoth1, ztota1)
cmas  write (413,'(f10.2, 5e14.4, 2i10)')
cmas *  timeice, ztoth0,ztota0, ztoth1,ztota1,
cmas *  (ztoth1-ztoth0)/max (ztota1, dd0**2), iterc, iloopmel

c        Re-impose boundary constraints,etc, on h for various expts  

c....................................................................
#if (defined (EISANTA) && defined (TRANSECTA)) || defined (TESTMEL2D)
c....................................................................
      do j=1,ny
        budgnegmeltmp(nx,j) = budgnegmeltmp(nx,j) -  h1(nx,j)/dt   !diag
        h1(nx,j) = 0.
      enddo
c.......................................................................
#elif ( defined (EISANTA) && ! defined (TRANSECTA) ) || defined (NHA) || defined (CARB) || defined (ORDO)
c.......................................................................
      if (nx.gt.1) then
        do j=1,ny
          budgnegmeltmp(1, j) = budgnegmeltmp(1, j) -  h1(1, j)/dt !diag
          budgnegmeltmp(nx,j) = budgnegmeltmp(nx,j) -  h1(nx,j)/dt !diag
          h1(1, j) = 0.
          h1(nx,j) = 0.
        enddo
      endif
      if (ny.gt.1) then
        do i=1,nx
          budgnegmeltmp(i, 1) = budgnegmeltmp(i, 1) -  h1(i, 1)/dt !diag
          budgnegmeltmp(i,ny) = budgnegmeltmp(i,ny) -  h1(i,ny)/dt !diag
          h1(i, 1) = 0.
          h1(i,ny) = 0.
        enddo
      endif
c.................................................
#elif defined (TESTMEL2D) && defined (TESTMELPROF)
c.................................................
      do j=1,ny
        h1(1,j) = 0.
        h1(nx,j) = 0.
      enddo
c.....
#endif
c.....

c       Apply 9-pt smoother, to reduce alternating-grid
c       oscillations in hmel and +/-/+/- velocity noise

      call zero (budgsmomeltmp, nx*ny)

#if defined (SMOOTHMEL)
      if (SMOOTHMEL.gt.0) then
        call scopy (nx*ny, h1, 1, hsmoprev, 1)
 
        do m=1,SMOOTHMEL
          call smoothermel (h1, maskmel)
        enddo
  
c       zo = 0.
c       zn = 0.
c       za = 0.
        do j=1,ny
          do i=1,nx
            if (maskmel(i,j).ne.0) then
              budgsmomeltmp(i,j) =  (h1(i,j) - hsmoprev(i,j)) / dt
c             conservation check:
c             zo = zo + hsmoprev(i,j)*darea(i,j)
c             zn = zn + h1(i,j)*darea(i,j)
c             za = za + darea(i,j)
            endif
          enddo
        enddo
c       conservation check:
c       zo = zo / max (za, .001)
c       zn = zn / max (za, .001)
c       write (555,'(a/e25.15/2e25.15)') 'smoother: ', zo, zn, zn-zo
      endif
#endif

c        Update h at t+0.5*dt to be used for next iter's setcoefs

      do j=1,ny
        do i=1,nx
          h(i,j) = 0.5*(h1(i,j) + h0(i,j))
        enddo
      enddo

c        Don't do this - inhibits growth in testmel2d:  
c        Zero small melange thicknesses (avoiding infinite veloc bc's
c        at open ocean). Done here to avoid small h's in u,v calcs i
c        in next C-iteration. Don't need to account in budgredmel...
c        also done on final h outside C-loop below.
c
c     do j=1,ny
c       do i=1,nx
c         if (h(i,j).lt.hminmel) then
c           h(i,j) = 0.
c         endif
c       enddo
c     enddo

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
      iterc = nitercmel
 3000 continue
c=======================

c        Set new h for this timestep (at t+dt)

      call scopy (nx*ny, h1, 1, h, 1)

c        Accumulate budget quantities (for meltabbud) over
c        multiple calls to icemel for each main model timestep.
c        (Need to do outside loop C above, just for last iteration).

      do j=1,ny
        do i=1,nx
          budgallmel(i,j) = budgallmel(i,j) +  budgallmeltmp(i,j) ! diag
          budgsupmel(i,j) = budgsupmel(i,j) +  budgsupmeltmp(i,j) ! diag
          budgnegmel(i,j) = budgnegmel(i,j) +  budgnegmeltmp(i,j) ! diag
          budgsmomel(i,j) = budgsmomel(i,j) +  budgsmomeltmp(i,j) ! diag
        enddo
      enddo

c        Skip point for no ice (numh=0, top of loop C)

 3500 continue    

c        Zero very small shelf amounts (h ~1.e-20) created by icedyn

      do j=1,ny
        do i=1,nx
          if (h(i,j).lt.1.e-10) h(i,j) = 0.
        enddo
      enddo

c        Zero small melange thicknesses (avoiding infinite veloc bc's
c        at open ocean). Conceptually should correspond to size of a
c        typical individual "chunk" (see notes 7/10/17). But can have 
c        smaller values, as long as not smaller than can be created
c        by cliff/calv in one timestep where no melange before.

      do j=1,ny
        do i=1,nx
          if (h(i,j).lt.hminmel) then
            budgredmel(i,j) = budgredmel(i,j) - h(i,j)/dt         ! diag
            h(i,j) = 0.
          endif
        enddo
      enddo

c       Apply uniform melange thickness limit if specified

      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).ne.0) then
            if (h(i,j).gt.hmaxmel) then
              budgredmel(i,j)= budgredmel(i,j)-(h(i,j)-hmaxmel)/dt !diag
              h(i,j) = hmaxmel
            endif
          endif
        enddo
      enddo

c        Set melange buttressing factor (theta[u,v], on h-grid), 
c        ratio of actual to free-floating divergence, as in thetacalc.
c        Used for buttressing of cliff faces in icedyn (schoofgl).  

      call thetacalcmel (h, fedge, maskmel, arhap, 
     *                   ua, va, thetau, thetav)

c        Final adjust melange surface elevations, water depths, maskmel,
c        etc (for change in h)

      call adjustmel (maskmel, h, hs, hw, fmel,
     *                hb, hice, fedge, maskwater, sealev,
     *                budgadjmel, dt, .false.)

c        Set average budget quantities (for meltabbud), after  
c        multiple calls to icemel for each main model timestep

      if (iloopmel.eq.nloopmel) then
        do j=1,ny
          do i=1,nx
            budgallmel(i,j) = budgallmel(i,j) / nloopmel          ! diag
            budgsupmel(i,j) = budgsupmel(i,j) / nloopmel          ! diag
            budgnegmel(i,j) = budgnegmel(i,j) / nloopmel          ! diag
            budgsmomel(i,j) = budgsmomel(i,j) / nloopmel          ! diag
            budgadjmel(i,j) = budgadjmel(i,j) / nloopmel          ! diag
            budgredmel(i,j) = budgredmel(i,j) / nloopmel          ! diag
          enddo
        enddo
      endif

c        Diagnostic calculations of total |force difference|
c        (curbackmel), and average curthetmel, etc, on ice faces

      if (iloopmel.eq.1) then
        curbackmel = 0.
        curthetmel = 0.
        curhmel = 0.
        curlmel = 0.
      endif

      call checkmelforce2 (curbackmel, curthetmel,
     *                     curhmel, curlmel,
     *                     h, hs, maskmel,
     *                     thetau, thetav,
     *                     ua, va, 
     *                     hice, fedge, hedge,
     *                     sealev, timeice, dt)

      if (iloopmel.eq.nloopmel) then
        curbackmel = curbackmel/nloopmel
        curthetmel = curthetmel/nloopmel
        curhmel = curhmel/nloopmel
        curlmel = curlmel/nloopmel

        curthetmel = curthetmel / max (curlmel, .001)
        curhmel = curhmel / max (curlmel, .001)
      endif

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if defined (TESTMEL2DDIAG1) || defined (TESTMEL2DDIAG2)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        iuzmel = 450
        j=(ny+1)/2
        write (iuzmel,'(/a,f10.2)') 'time=',timeice

c       find beginning of melange supply (assume flow left to right):
        zs = 0.
        i0 = 0
        do i=1,nx
          if (zs.eq.0. .and. budgsupmel(i,j).gt.0.) then
            zs = (budgsupmel(i,j)+budgallmel(i,j))*dd0
            zh0 = h(i,j)                   ! upstream
c           zh0 = 0.5*(h(i,j)+h(i+1,j))    ! centered
            zu0 = zs/zh0
            zx0 = xh(i,j) + 0.5*dd0
            i0 = i
          endif
        enddo

#  if defined (TESTMEL2DDIAG1)

c       diagnostic table to check h,u vs. analytic soln (notes 7/11/17)
c       for no visc (coefdivmel=0): Only forces are driving stress, 
c       water friction, oceanmelt (budgallmel), budgsupmel.
c       zx, zu, zh (and zx0, zh0, zs, za, zm) are for analytic soln.

        if (i0.ne.0) then
          za = coefbwatermel/(rhomp*grav)
          zm = budgallmel(i0,j)

          do i=i0,nx-1
            zx = xh(i,j) + 0.5*dd0 - zx0
            zh3 = zh0**3 - 3.*za*(zs*zx + 0.5*zm*(zx**2))
            if (zh3.ge.0.) then
              zh = zh3**(1./3.)
              zu = (zs + zm*zx) / zh
            else
              zh = 0.
              zu = 0.
            endif
          
            write (iuzmel,'(i4, 2f10.3, 5x, 4f10.3)') 
     *        i, 0.5*(h(i,j)+h(i+1,j)), ua(i,j), zh, zu, zs/dd0, zm
          enddo
        endif

#  elif defined (TESTMEL2DDIAG2)

c       diagnostic table to check h,u vs. analytic soln (notes 9/27/17,
c       like Oerlemans & V. der V, pg. 64), for visc (coefdivmel) only,
c       no drag. Only forces are driving stress, negligible water 
c       friction, budgallmel, budgsupmel.
c       zx, zu, zh (and zx0, zh0, zc, zm) are for analytic soln.

        if (i0.ne.0) then
          zc = ( (rhomp*grav*arhap(i0,j))/(4.*coefdivmel) ) ** powimel
          zm = budgallmel(i0,j)

          izp = powimel+1

          do i=i0,nx-1
            zx = xh(i,j) + 0.5*dd0 - zx0
            zz = (zc/zm) - (zu0**izp)*((zc/zm)*(zh0**izp)-1.)
     *                   / ((zm*zx + zh0*zu0)**izp)
            if (zz.gt.0.) then
              zh = 1./(zz**(1./izp))
              zu = (zm*zx + zh0*zu0) / zh
              if (zu.lt.0.) then
                zh = 0.
                zu = 0. 
              endif
            else
              zh = 0.
              zu = 0. 
            endif
          
            write (iuzmel,'(i4,f10.3, 3x, 2f10.3, 3x, 4f10.3,2f12.6)') 
     *        i, zx*1.e-3, 
     *        0.5*(h(i,j)+h(i+1,j)), ua(i,j)*1.e-3,
     *        zh, zu*1.e-3, zs/dd0, zm, zc, thetau(i,j)
          enddo
        endif

#  endif
c~~~~~
#endif
c~~~~~

      return
      end

c-----------------------------------------------------------------------

      subroutine setcoefsmel (coef1, coef2, visc, coefbu, coefbv, 
     *                        maskh, arhap, h, hbsd, hw,
     *                        masku, maskv, ua, va, 
     *                        hu, hv, crhu, crhv,
     *                        muind, mvind, maskc, itera, iterc)

c        Set various coefficients, etc, at top of A loop: 
c        coef[1,2], coefb[u,v]

#include <comicephys.h>
#include <comicegrid.h>

c     passed:
      dimension
     *  coef1(nx,ny),
     *  coef2(0:nx,0:ny),                                   ! cross grid
     *  visc(nx,ny), 
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  maskh(nx,ny),          arhap (nx,ny), 
     *  h(nx,ny),              hbsd(nx,ny),           hw(nx,ny), 
     *  masku(0:nxp,0:nyp),    maskv(0:nxp,0:nyp),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),     crhv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp), 
     *  maskc(0:nx,0:ny)                                    ! cross grid

c     local:
      dimension
     *  g1(nx,ny), g2(nx,ny), g3(nx,ny),
     *  g2c(0:nx,0:ny)                                      ! cross grid
      save g1, g2, g3, g2c

c     Statement function (fpin) for fraction of u,v grid box with 
c     sub-grid bedrock mounds scraping ice bottom:

#include <comicepin.h>

c        Zero initially (some necessary for finite diffs near ice edges)

      if (itera.eq.1 .and. iterc.eq.1) then
        call zero (coef1, nx*ny)
        call zero (coef2,(nx+1)*(ny+1))
        call zero (visc,  nx*ny) 
        call zero (g1,    nx*ny)
        call zero (g2,    nx*ny)
        call zero (g3,    nx*ny)
      endif

c        Set g1 ((du/dx)**2 + (dv/dy)**2 + (du/dx)*(dv/dy)) (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            zdudx = (ua(i,j)-ua(i-1,j))/dx(i,j)
            zdvdy = (va(i,j)-va(i,j-1))/dy(i,j)

c           adjust visc for melange rheol - not done 7/23/17:
c           if (zdudx.gt.0.) then
c             zdudx = zdudx*coefdivmel
c           else
c             zdudx = zdudx*coefconmel
c           endif
c           if (zdvdy.gt.0.) then
c             zdvdy = zdvdy*coefdivmel
c           else
c             zdvdy = zdvdy*coefconmel
c           endif

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

c        Set visc (h grid)

      do j = 1,ny
        do i = 1,nx
          if (maskh(i,j).eq.1) then
            visc(i,j) = 0.5 
     *                / (max(g1(i,j)+g2(i,j)+g3(i,j), gmin) ** powivmel)
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
          if (masku(i,j).eq.1) then

c           calculate zfpin, fraction of u-grid cell actually in 
c           contact with bottom, scraped by sub-grid bedrock mounds 
c           (stmt fnct fpin):

c           Different from icedyn (shelves) 10/17:
            zfpin1 = fpin(hw(i,j),  hbsd(i,j))
            zfpin2 = fpin(hw(i+1,j),hbsd(i+1,j))
            zfpin = 0.5*(zfpin1+zfpin2)

c           calculate coefbu (N/m2 per m/y), weighting by water vs.
c           non-water contact fractions: 
            if (zfpin.eq.0.) then
              coefbu(i,j) = coefbwatermel
            else 
              zc = crhu(i,j) ** (-1./powbmel)
              zpow = -0.5*(1. - 1./powb)
              zu = ua(i,j) 
              zv = 0.25*(va(i,j) + va(i+1,j) + va(i,j-1) + va(i+1,j-1))
              zcoef = zc * ( max(zu**2+zv**2, ubmin**2)**zpow )
              zcoef = max (zcoef, coefbwatermel)
              coefbu(i,j) = (1.-zfpin) * coefbwatermel
     *                    +     zfpin  * zcoef
            endif
          else
            coefbu(i,j) = 0. ! not used
          endif
        enddo
      enddo

      do j = 1,ny-1
        do i = 1,nx
          if (maskv(i,j).eq.1) then

            zfpin1 = fpin(hw(i,j),  hbsd(i,j))
            zfpin2 = fpin(hw(i,j+1),hbsd(i,j+1))
            zfpin = 0.5*(zfpin1+zfpin2)

            if (zfpin.eq.0.) then
              coefbv(i,j) = coefbwatermel
            else
              zc = crhv(i,j) ** (-1./powb)
              zpow = -0.5*(1. - 1./powb)
              zu = 0.25*(ua(i,j) + ua(i,j+1) + ua(i-1,j) + ua(i-1,j+1))
     *                   
              zv = va(i,j)
              zcoef = zc * ( max(zu**2+zv**2, ubmin**2)**zpow )
              zcoef = max (zcoef, coefbwatermel)
              coefbv(i,j) = (1.-zfpin) * coefbwatermel
     *                    +     zfpin  * zcoef
            endif
          else
            coefbv(i,j) = 0. ! not used
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine setvelbc (uabound, vabound, 
     *                     h, maskmel, supcalv,
     *                     ucimel, vcimel, hucimel, hvcimel, 
     *                     ucamel, vcamel, hucamel, hvcamel, 
     *                     uaice, vaice, hice, fedge)

c        Set melange veloc bc's at calving/cliff ice faces

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  uabound(0:nxp,0:nyp),  vabound(0:nxp,0:nyp),
     *  h(nx,ny),              maskmel(nx,ny),         supcalv(nx,ny),
     *  ucimel(0:nxp,0:nyp),   vcimel(0:nxp,0:nyp),
     *  hucimel(0:nxp,0:nyp),  hvcimel(0:nxp,0:nyp),
     *  ucamel(0:nxp,0:nyp),   vcamel(0:nxp,0:nyp),
     *  hucamel(0:nxp,0:nyp),  hvcamel(0:nxp,0:nyp),
     *  uaice(0:nxp,0:nyp),    vaice(0:nxp,0:nyp),
     *  hice(nx,ny),           fedge(nx,ny)

      call zero (uabound, (nxp+1)*(nyp+1))
      call zero (vabound, (nxp+1)*(nyp+1))

c     u grid:

      do j=1,ny
        do i=1,nx
          ip = min (i+1,nx)
          idir = 0
          if ( (hice(i,j).gt.0. .and. fedge(i,j).eq.1.
     *           .and. maskmel(i,j).eq.0) .and.
     *         (h(ip,j).gt.0. .and. maskmel(ip,j).ne.0)
     *       ) then
            idir = 1
            imel = ip
            iice = i
          else if ( (hice(ip,j).gt.0. .and. fedge(ip,j).eq.1.
     *                .and. maskmel(ip,j).eq.0) .and. 
     *              (h(i,j).gt.0. .and. maskmel(i,j).ne.0)
     *            ) then
            idir = -1
            imel = i
            iice = ip
          endif
          if (iabs(idir).eq.1) then
            zca = supcalv(imel,j)*dx(imel,j) / hice(iice,j)
            zha = hice(iice,j)
            zfa = zca*zha
            zci = ucimel(i,j)
            zhi = hucimel(i,j)
            zfi = zci*zhi
            zcer  = (zci*zfi + zca*zfa) / max (zfi+zfa,.001)
#if defined (HNEWMEL)
            zhm = min ( float(HNEWMEL), hice(iice,j)*rhoice/rhomel )
#else
            zhm = min (30., hice(iice,j))
c           zhm = max (20.,  min (150., min (h(imel,j),hice(iice,j))))
c           zhm = min (hice(iice,j), max (150., h(imel,j)))
c           zhm = min (hice(iice,j),      150.            )
#endif
            uabound(i,j) = uaice(i,j) 
     *                   + idir * zcer 
     *                     * ((rhoice*hice(iice,j))/(rhomel*zhm) - 1.)
          endif
        enddo
      enddo

c     v grid:

      do j=1,ny
        jp = min (j+1,ny)
        do i=1,nx
          idir = 0
          if ( (hice(i,j).gt.0. .and. fedge(i,j).eq.1. 
     *           .and. maskmel(i,j).eq.0) .and. 
     *         (h(i,jp).gt.0. .and. maskmel(i,jp).ne.0)
     *       ) then
            idir = 1
            jmel = jp
            jice = j
          else if ( (hice(i,jp).gt.0. .and. fedge(i,jp).eq.1.
     *                .and. maskmel(i,jp).eq.0) .and. 
     *              (h(i,j).gt.0. .and. maskmel(i,j).ne.0)
     *            ) then
            idir = -1
            jmel = j
            jice = jp
          endif
          if (iabs(idir).eq.1) then
            zca = supcalv(i,jmel)*dy(i,jmel) / hice(i,jice)
            zha = hice(i,jice)
            zfa = zca*zha
            zci = vcimel(i,j)
            zhi = hvcimel(i,j)
            zfi = zci*zhi
            zcer  = (zci*zfi + zca*zfa) / max (zfi+zfa,.001)
#if defined (HNEWMEL)
            zhm = min ( float(HNEWMEL), hice(i,jice)*rhoice/rhomel )
#else
            zhm = min (30., hice(i,jice))
c           zhm = min (100., hice(i,jice))
c           zhm = max (20.,  min (150., min (h(i,jmel),hice(i,jice))))
c           zhm = min (hice(i,jice), max (150., h(i,jmel)))
c           zhm = min (hice(i,jice),      150.            )
#endif
            vabound(i,j) = vaice(i,j)
     *                   + idir * zcer
     *                     * ((rhoice*hice(i,jice))/(rhomel*zhm) - 1.)
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine dotridiamel (ua, 
     *                        hu, dfu, hw, coef1, coefbu, 
     *                        muind,nuvtot)

c     Assembles *tri arrays for 1-D tridiagonal solution.
c     Calls tridia_i, sets new velocities ua.
c     Only for 1-D in x with ny=1 (use dotridiamel2 for 1-D in y, nx=1).

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  ua(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       dfu(0:nxp,0:nyp),    hw(nx,ny), 
     *  coef1(nx,ny),          coefbu(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp)

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

          rtri(n1) = dfu(i,j) 

          z1 = 2.*coef1(i+1,j) / (dx(i+1,j)*dxu(i,j))  ! ua(i+1,j)
          z2 = 2.*coef1(i  ,j) / (dx(i,j  )*dxu(i,j))  ! ua(i-1,j)

c            Adjust z1,z2 for melange rheol

          if (ua(i+1,j)-ua(i,j).gt.0.) then
            z1 = z1*coefdivmel
          else 
            z1 = z1*coefconmel
          endif
          if (ua(i,j)-ua(i-1,j).gt.0.) then
            z2 = z2*coefdivmel
          else 
            z2 = z2*coefconmel
          endif

c         diagonal:
          atri(n1) = - (z1 + z2) - coefbu(i,j)

c         for side drag (no slip, using visc, * sidedragmel):
          z5 = sidedragmel * 0.5*(0.5*(coef1(i+1,j)+coef1(i,j)))
     *         / (0.5*dyu(i,j)*dyu(i,j))
          atri(n1) = atri(n1) - z5 

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

        endif

      enddo


      call tridia_imel (atri, btri, ctri, vtri, rtri, nuvtot, 0)

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

c|||||||||||||||||||||||||
#if defined (DOSPARSEMELU)
c|||||||||||||||||||||||||

      subroutine dosparsemelu (ua, va, 
     *                         hu, hv, dfu, dfv,
     *                         coef1, coef2, coefbu, coefbv,
     *                         muind, mvind, nuvtot, hw, maskmel,
     *                         timeice, itera, iterc, iloopmel
     *                        )

c     Solves elliptic shelf equations for ua,va, 
c     just for points with m[u,v]ind >= 1. 
c     Solves using choice of sparse solvers (or Gaussian elimination). 
c     Assembles elspa, ijspa, rhs arrays for sparse solver,
c     calls sparse solver, sets new velocities ua,va.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

      dimension 
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),
     *  coef1(nx,ny),          coef2(0:nx,0:ny),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  hw(nx,ny),             maskmel(nx,ny)

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

              rhs(n1) = dfu(i,j)

              z1 = 2.*coef1(i+1,j) / (dx(i+1,j )*dxu(i,j))  ! ua(i+1,j)
              z2 = 2.*coef1(i  ,j) / (dx(i,j   )*dxu(i,j))  ! ua(i-1,j)
              z3 =    coef2(i  ,j) / (dyc(i,j  )*dyu(i,j))  ! ua(i,j+1)
              z4 =    coef2(i,j-1) / (dyc(i,j-1)*dyu(i,j))  ! ua(i,j-1)

c                Adjust z1,z2 for melange rheol

              if (ua(i+1,j)-ua(i,j).gt.0.) then
                z1 = z1*coefdivmel
              else 
                z1 = z1*coefconmel
              endif
              if (ua(i,j)-ua(i-1,j).gt.0.) then
                z2 = z2*coefdivmel
              else 
                z2 = z2*coefconmel
              endif

c                diagonal for ua(i,j):

              elspa(n1) = elspa(n1) - (z1 + z2 + z3 + z4) - coefbu(i,j)

c---
c                Side drag (no slip, using visc, * sidedragmel):

              zside = 0.
#if defined (TESTMEL2D)
              if (nx.ge.ny .and. (j.eq.1 .or. j.eq.ny)) zside = 1.
#else
              if (j.gt.1) then 
                if (maskmel(i,  j-1).eq.0) zside = zside + 0.5
                if (maskmel(i+1,j-1).eq.0) zside = zside + 0.5
              endif
              if (j.lt.ny) then 
                if (maskmel(i,  j+1).eq.0) zside = zside + 0.5
                if (maskmel(i+1,j+1).eq.0) zside = zside + 0.5
              endif
#endif
              z5 = sidedragmel * 0.5*(0.5*(coef1(i+1,j)+coef1(i,j)))
     *           / (0.5*dyu(i,j)*dyu(i,j))
              elspa(n1) = elspa(n1) - zside*z5
c---

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

c                Adjust za,zb for melange rheol

              if (va(i+1,j)-va(i+1,j-1).gt.0.) then
                za = za*coefdivmel
              else 
                za = za*coefconmel
              endif
              if (va(i,j)-va(i,j-1).gt.0.) then
                zb = zb*coefdivmel
              else 
                zb = zb*coefconmel
              endif

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

              rhs(n1) = dfv(i,j)

              z1 = 2.*coef1(i,j+1) / (dy(i,j+1 )*dyv(i,j))  ! va(i,j+1)
              z2 = 2.*coef1(i  ,j) / (dy(i,j   )*dyv(i,j))  ! va(i,j-1)
              z3 =    coef2(i  ,j) / (dxc(i,j  )*dxv(i,j))  ! va(i+1,j)
              z4 =    coef2(i-1,j) / (dxc(i-1,j)*dxv(i,j))  ! va(i-1,j)

c                Adjust z1,z2 for melange rheol

              if (va(i,j+1)-va(i,j).gt.0.) then
                z1 = z1*coefdivmel
              else 
                z1 = z1*coefconmel
              endif
              if (va(i,j)-va(i,j-1).gt.0.) then
                z2 = z2*coefdivmel
              else 
                z2 = z2*coefconmel
              endif
 
c                diagonal for va(i,j):

              elspa(n1) = elspa(n1) - (z1 + z2 + z3 + z4) - coefbv(i,j)

c---
c                Side drag (no slip, using visc, * sidedragmel):

              zside = 0.
#if defined (TESTMEL2D)
              if (ny.ge.nx .and. (i.eq.1 .or. i.eq.nx)) zside = 1.
#else
              if (i.gt.1) then 
                if (maskmel(i-1,j  ).eq.0) zside = zside + 0.5
                if (maskmel(i-1,j+1).eq.0) zside = zside + 0.5
              endif
              if (i.lt.nx) then 
                if (maskmel(i+1,j  ).eq.0) zside = zside + 0.5
                if (maskmel(i+1,j+1).eq.0) zside = zside + 0.5
              endif
#endif
              z5 = sidedragmel * 0.5*(0.5*(coef1(i,j)+coef1(i,j+1)))
     *           / (0.5*dxv(i,j)*dxv(i,j))
              elspa(n1) = elspa(n1) - zside*z5
c---

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

c                Adjust za,zb for melange rheol

              if (ua(i,j+1)-ua(i-1,j+1).gt.0.) then
                za = za*coefdivmel
              else 
                za = za*coefconmel
              endif
              if (ua(i,j)-ua(i-1,j).gt.0.) then
                zb = zb*coefdivmel
              else 
                zb = zb*coefconmel
              endif
  
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

c^^^^^^^^^^^^^^
          endif   ! n1 >= 1
c^^^^^^^^^^^^^^

        enddo
      enddo

      if (kof.gt.nspamax) then
        write (6,*) 'Error (dosparsemelu): kof > nspamax. kof=', kof,
     *              ' nspamax=', nspamax
        stop
      endif 

      ijspa(nuvtot+1) = kof + 1

      if (nuvtot.ne.ijspa(1) - 2) then
        write (6,*) 'Error (dosparsemel): nuvtot, ijspa(1)-2 =',
     *              nuvtot, ijspa(1)-2
        stop
      endif

c        Solve linear system arr*vec = rhs for vec.
c        If number of pts <= nuvsmall, use gaussian elimination 
c        (avoids sparse-solution problems for nuvtot <= ~2)

c+++++++++++++++++++++++++++++++++
      if (nuvtot.gt.nuvsmall) then
c+++++++++++++++++++++++++++++++++

#if defined (NUMRECMELU) || defined (SLAPMELU)
      if (firstvecu) then
        call zero (vecu, nuvtot)
        firstvecu = .false.
      endif
#endif
#if defined (NUMRECMELU)
      itol = 2 ! 4 ! 1
      tol = 1.e-4  ! 1.e-4 ! 1.e-6 ! 1.e-11
      itmax = 10000
      call linbcg (nuvtot, rhs, vecu, itol, tol, itmax, itlin, errlin)
c     write (ioterm,'(a,i5,a,e15.5)')
c    *  'linbcg(melu): itlin=',itlin,' errlin=',errlin
      if (abs(errlin).gt.tol) write (ioterm,'(a,f10.2,a,i5,a,e15.5)')
     *  '*** Warning linbcg(melu): time=',timeice,
     *  ' itlin=',itlin,' errlin=',errlin

#elif defined (SLAPMELU)
      itol = 2
      tol = 1.e-4
      itmax = 10000
      call slap (nuvtot, rhs, vecu, itol, tol, itmax, itlin, errlin)
c     write (ioterm,'(a,i5,a,e15.5)')
c    *  'slap(melu): itlin=',itlin,' errlin=',errlin
      if (abs(errlin).gt.tol) write (ioterm,'(a,i5,a,e15.5)')
     *  '*** Warning slap(melu): itlin=',itlin,' errlin=',errlin
#elif defined (SUPERLUMELU)
      call superlu (nuvtot, rhs, vecu, info)
      if (info.ne.0)write (ioterm,*)'*** Warn superlu(melu): info=',info
#elif defined (WATSONMELU)
      call watson (nuvtot, rhs, vecu, info)
      if (info.ne.0)write (ioterm,*) '*** Warn watson(melu): info=',info
#elif defined (MKLMELU)
      call mkl (nuvtot, rhs, vecu, info)
      if (info.ne.0) write (ioterm,*)'*** Warning mkl(melu): info=',info
#elif defined (GAUSSMELU)
      call gaussdo (nuvtot, rhs, vecu, info)
      if (info.ne.0) write (ioterm,*)'*** Warn gauss(melu): info=',info
#endif

c+++++++++
      else
c+++++++++
        call gaussdo (nuvtot, rhs, vecu, info)
        if (info.ne.0) write (ioterm,*)'*** Warn gauss(melu):info=',info
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

c============================
#if defined (DOSPARSESORMELU)
c============================

      subroutine dosparse_sormelu
     *                    (ua, va, 
     *                     hu, hv, dfu, dfv,
     *                     coef1, coef2, coefbu, coefbv, 
     *                     muind, mvind, h, hs, hb, hw, maskmel,
     *                     timeice, itera, iterc, iloopmel)

c     Solevs elliptic shelf equations for ua,va, using SOR (Successive 
c     Over Relaxation), just for points with m[u,v]ind >= 1. 

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

      dimension 
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),       hv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),      dfv(0:nxp,0:nyp),
     *  coef1(nx,ny),          coef2(0:nx,0:ny),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  h(nx,ny),              hs(nx,ny),
     *  hb(nx,ny),             hw(nx,ny),            maskmel(nx,ny)

c     parameter (nitersurmax =  500, nitersurmin = 20, facsor=0.9)
c     parameter (nitersurmax =  100, nitersurmin = 20, facsor=1.2) ! std
      parameter (nitersurmax =   20, nitersurmin = 10, facsor=1.4) ! new
      parameter (ducrit      = .050)                               ! new

#if ! defined (NOCHECKERBOARD)
c     for checkerboard sor:
      logical firstcheck
      data firstcheck /.true./
      dimension mcheck(0:nxp, 0:nxp)
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

c              Adjust z1,z2 for melange rheol

            if (ua(i+1,j)-ua(i,j).gt.0.) then
              z1 = z1*coefdivmel
            else 
              z1 = z1*coefconmel
            endif
            if (ua(i,j)-ua(i-1,j).gt.0.) then
              z2 = z2*coefdivmel
            else 
              z2 = z2*coefconmel
            endif

c              v-neighbors for ua(i,j):
    
            za = coef1(ip1,j) / (dy(ip1,j )*dxu(i,j)) ! i+1,j; i+1,j-1
            zb = coef1(i  ,j) / (dy(i,j   )*dxu(i,j)) ! i,j  ; i,j-1
            zc = coef2(i  ,j) / (dxc(i,j  )*dyu(i,j)) ! i+1,j; i,j
            zd = coef2(i,j-1) / (dxc(i,j-1)*dyu(i,j)) ! i,j-1; i+1,j-1

c                Adjust za,zb for melange rheol

            if (va(i+1,j)-va(i+1,j-1).gt.0.) then
              za = za*coefdivmel
            else 
              za = za*coefconmel
            endif
            if (va(i,j)-va(i,j-1).gt.0.) then
              zb = zb*coefdivmel
            else 
              zb = zb*coefconmel
            endif

            zforcu = 
     *          z1*(ua(i+1,j)        )     - z2*(       -ua(i-1,j)) 
     *        + z3*(ua(i,j+1)        )     - z4*(       -ua(i,j-1)) 

     *        + za*(va(i+1,j)-va(i+1,j-1)) - zb*(va(i,j)-va(i,j-1))
     *        + zc*(va(i+1,j)-va(i,j))     - zd*(va(i+1,j-1)-va(i,j-1))

     *        - dfu(i,j) 

            zclhsu = z1 + z2 + z3 + z4 + coefbu(i,j)

c---
c              Side drag (no slip, using visc, * sidedragmel):

            zside = 0.
#if defined (TESTMEL2D)
            if (nx.ge.ny .and. (j.eq.1 .or. j.eq.ny)) zside = 1.
#else
            if (j.gt.1) then 
              if (maskmel(i,  j-1).eq.0) zside = zside + 0.5
              if (maskmel(i+1,j-1).eq.0) zside = zside + 0.5
            endif
            if (j.lt.ny) then 
              if (maskmel(i,  j+1).eq.0) zside = zside + 0.5
              if (maskmel(i+1,j+1).eq.0) zside = zside + 0.5
            endif
#endif
            z5 = sidedragmel * 0.5*(0.5*(coef1(ip1,j)+coef1(i,j)))
     *           / (0.5*dyu(i,j)*dyu(i,j))
            zclhsu = zclhsu + zside*z5
c---

c@@@ std:
c           zdu = zforcu/zclhsu - ua(i,j)
c           ua(i,j) = facsor*zforcu/zclhsu + (1.-facsor)*ua(i,j)
c           zdumax = max (zdumax, abs(zdu))
c@@@ new:
c           Set new ua (zforcu/zclhsu).
c           Modify by facsor (std sor).
            zf = facsor
            uaold = ua(i,j)
            ua(i,j) = zf * (zforcu/zclhsu)
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

c              Adjust z1,z2 for melange rheol

            if (va(i,j+1)-va(i,j).gt.0.) then
              z1 = z1*coefdivmel
            else 
              z1 = z1*coefconmel
            endif
            if (va(i,j)-va(i,j-1).gt.0.) then
              z2 = z2*coefdivmel
            else 
              z2 = z2*coefconmel
            endif

c              u-neighbors for va(i,j):
 
            za = coef1(i,jp1) / (dx(i,jp1 )*dyv(i,j)) ! i,j+1; i-1,j+1
            zb = coef1(i  ,j) / (dx(i,j   )*dyv(i,j)) ! i,j  ; i-1,j
            zc = coef2(i  ,j) / (dyc(i,j  )*dxv(i,j)) ! i,j+1; i,j
            zd = coef2(i-1,j) / (dyc(i-1,j)*dxv(i,j)) ! i-1,j; i-1,j+1

c                Adjust za,zb for melange rheol

            if (ua(i,j+1)-ua(i-1,j+1).gt.0.) then
              za = za*coefdivmel
            else 
              za = za*coefconmel
            endif
            if (ua(i,j)-ua(i-1,j).gt.0.) then
              zb = zb*coefdivmel
            else 
              zb = zb*coefconmel
            endif

            zforcv = 
     *          z1*(va(i,j+1)        )     - z2*(       -va(i,j-1)) 
     *        + z3*(va(i+1,j)        )     - z4*(       -va(i-1,j)) 

     *        + za*(ua(i,j+1)-ua(i-1,j+1)) - zb*(ua(i,j)-ua(i-1,j))
     *        + zc*(ua(i,j+1)-ua(i,j))     - zd*(ua(i-1,j+1)-ua(i-1,j))

     *        - dfv(i,j)

            zclhsv = z1 + z2 + z3 + z4 + coefbv(i,j)

c---
c              Side drag (no slip, using visc, * sidedragmel):

            zside = 0.
#if defined (TESTMEL2D)
            if (ny.ge.nx .and. (i.eq.1 .or. i.eq.nx)) zside = 1.
#else
            if (i.gt.1) then 
              if (maskmel(i-1,j  ).eq.0) zside = zside + 0.5
              if (maskmel(i-1,j+1).eq.0) zside = zside + 0.5
            endif
            if (i.lt.nx) then 
              if (maskmel(i+1,j  ).eq.0) zside = zside + 0.5
              if (maskmel(i+1,j+1).eq.0) zside = zside + 0.5
            endif
#endif
            z5 = sidedragmel * 0.5*(0.5*(coef1(i,jp1)+coef1(i,j)))
     *           / (0.5*dxv(i,j)*dxv(i,j))
            zclhsu = zclhsu + zside*z5
c---

c@@@ std:
c           zdv = zforcv/zclhsv - va(i,j)
c           va(i,j) = facsor*zforcv/zclhsv + (1.-facsor)*va(i,j)
c           if (zwei.eq.0.) zdumax = max (zdumax, abs(zdv))
c@@@ new (see comments above for u):
            zf = facsor
            vaold = va(i,j)
            va(i,j) = zf * (zforcv/zclhsv)
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
     *       'dosparse_sormelu: time=', nint(timeice)

        if (itersur.eq.1) 
     *    write (6,'(/a,3i8,f12.3/a)') 
     *       'iterc, itera, ndav, zuav=', iterc, itera, ndav, zuav,
     *       '  itersur      zdumax      ducrit'

c       if (itersur.eq.1 .or. itersur.eq.nitersurmax .or. 
c    *     (zdumax.lt.ducrit .and. itersur.ge.nitersurmin)) then
          write (6,'(i9, 2f12.3, i9)') itersur, zdumax, ducrit, iloopmel
c       endif
      endif

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

      subroutine tridia_imel (a, b, c, h, forc, n, ktrid)

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
        write(6,*)'*** Error in tridia_imel: Too many points.',
     *            '   nmax=',nmax,'   n=',n
        stop
      endif

c     amin = 1.e-3
      amin = 1.e-6
      do i = 1,n
        if (abs(a(i)).lt.amin) then
          write (6,9000) i,a(i)
 9000     format(' tridia_imel warning: a(',i3,')=',f12.8)
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

      subroutine dosparsemelh (ua, va, 
     *                         h0, h, h1, maskmel,
     *                         budgall, budgsup, budgneg,
     *                         timeice, dt, iterc, iloopmel)

c     Steps elliptic advection equations for melange thickness
c     (h0 to h1) through dt, just for points with h>0 or adjacent,
c     with maskmel != 0.
c     Uses upstream thickness (not centered), time implicit,
c     calls sparse solver, sets new melange thicknesses h1.
c     Assembles elspa, ijspa, rhs arrays for sparse solver,
c     elspa, ijspa, are in comicesparse.h, shared with solver.
c     vech is saved for subsequent calls to linbcg.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

      dimension 
     *  ua(0:nxp,0:nyp),    va(0:nxp,0:nyp),
     *  h0(nx,ny),          h(nx,ny),         
     *  h1(nx,ny),          maskmel(nx,ny),
     *  budgall(nx,ny),     budgsup(nx,ny),    budgneg(nx,ny)

      dimension
     *  mhind(nx,ny), vech(nuvmax)
      save vech

c     for "dotria" 1-D only:
      dimension atri(nx), btri(nx), ctri(nx), rtri(nx), vtri(nx)

      parameter (nlook = 4)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /

      logical firstvech
      data firstvech /.true./ 
      save firstvech

      dimension rhs(nuvmax)

cmat  for forming 2-D matrix from elspa, ijspa (linbcg): 
cmat  parameter (nhmat = 5000) 
cmat  dimension asparse(nhmat,nhmat), colsum(nhmat), 
cmat *          indsav(nx,ny), jndsav(nx,ny)

c        Set nhtot and mhind indices from 2-D arrays into 
c        linear sequence of selected points for sparse soln
c        below 
    
      call izero (mhind, nx*ny)
      nhtot = 0
    
      do j = 1,ny
        do i = 1,nx
          if (maskmel(i,j).ne.0) then
c           if (h(i,j).gt.0.) then           ! h0 not h to ensure advect
            if (h0(i,j).gt.0.) then          ! only out of h0>0 below
              nhtot = nhtot + 1
              mhind(i,j) = nhtot
cmat          indsav(nhtot) = i
cmat          jndsav(nhtot) = i
            else
              do look=1,4
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (.not.(ii.eq.i .and. jj.eq.j)) then    ! for flowline
c                 if (h(ii,jj).gt.0.) then  ! h0 not h to  ensure advect
                  if (h0(ii,jj).gt.0.) then ! only out of h0>0 below
                    nhtot = nhtot + 1
                    mhind(i,j) = nhtot
cmat                indsav(nhtot) = i
cmat                jndsav(nhtot) = i
                    goto 10
                  endif
                endif
              enddo
   10         continue
            endif
          endif
        enddo
      enddo

c#######################
      if (ny.eq.1) then
c#######################

c     Assembles *tri arrays for 1-D tridiagonal solution.
c     Calls tridia_i, sets new velocities ua.
c     Only for 1-D in x with ny=1 (use dotridia2 for 1-D in y, nx=1).


c        Tridiagonal solution for ua in 1-D problems (in x, with ny=1).
c        Relies on indexing running left to right (gaps in ice are ok). 
     
      j = 1
     
      call zero (atri, nx)
      call zero (btri, nx)
      call zero (ctri, nx)
      call zero (rtri, nx)
      call zero (vtri, nx)
     
      do i=1,nx
        n1 = mhind(i,j)
     
        if (n1.ge.1) then
     
          rtri(n1) = h0(i,j) + (budgall(i,j)+budgsup(i,j)) * dt
     
          z1 = dt*dyu(i-1,j)/darea(i,j)
          z2 = dt*dyu(i,j)  /darea(i,j)
     
c         diagonal:
          atri(n1) = 1.

          if (i.lt.nx) then
            if (mhind(i+1,j).ge.1) then
c             upstream:
              atri(n1) = atri(n1) + max (ua(i,j),0.)  *z2
            endif
          endif

          if (i.gt.1) then
            if (mhind(i-1,j).ge.1) then
              atri(n1) = atri(n1) - min (ua(i-1,j),0.)*z1
            endif
          endif
     
c         superdiagonal:
          if (i.lt.nx) then
            if (mhind(i+1,j).ge.1) then
c             upstream:
              btri(n1) = min (ua(i,j),0.)*z2
            endif
          endif

c         subdiagonal:
          if (i.gt.1) then
            if (mhind(i-1,j).ge.1) then
c             upstream:
              ctri(n1) = - max (ua(i-1,j),0.)*z1
            endif
          endif

        endif

      enddo

      call tridia_imel (atri, btri, ctri, vtri, rtri, nhtot, 1)

c        Set new h

      do i=1,nx
        n1 = mhind(i,j)
        if (n1.ge.1) then
          h1(i,j) = vtri(n1)
        endif
      enddo

      goto 1000

c##########
      endif
c##########

c        Set sparse arrays and rhs's needed to solve system of eqns.
c        (using elspa,ijspa,rhs). Can't parallelize for 
c        sparse method, since off-diag terms in ijspa,elspa are 
c        assembled using incremented index kof. Must loop through 
c        i,j in same sequence as mhind was assembled
c        above, so that diagonal elements are in same order as
c        off-diagonal column groups in elspa,ijspa.
c
c        If the neighbor's h iw not being computed (n2 = 0), 
c        velocities across boundary will be zero.

      call izero (ijspa, nspamax)
      call  zero (elspa, nspamax)
      call  zero (rhs,   nhtot)

      kof = nhtot + 1

      do j=1,ny
        do i=1,nx

          n1 = mhind(i,j)
c<<<<<<<<<<<<<<<<<<<<<<<<<<
          if (n1.ge.1) then
c<<<<<<<<<<<<<<<<<<<<<<<<<<
            kofprev = kof

c           rhs for h(i,j):

            rhs(n1) = h0(i,j) + (budgall(i,j)+budgsup(i,j)) * dt

            z1 = dt*dyu(i-1,j)/darea(i,j)
            z2 = dt*dyu(i,j)  /darea(i,j)
            z3 = dt*dxv(i,j-1)/darea(i,j)
            z4 = dt*dxv(i,j)  /darea(i,j)

c              diagonal for h(i,j):

            elspa(n1) = 1.    

c                       terms added in if's below to avoid exchange 
c                       with pts not in matrix: (possible because 
c                       ua,va are set from h/masku, changing in loop C,
c                       but only h0>0 and adjacent pts are included in
c                       matrix here). 
c    *                  - min (ua(i-1,j),0.)*z1
c    *                  + max (ua(i,j),0.)  *z2 
c    *                  - min (va(i,j-1),0.)*z3
c    *                  + max (va(i,j),0.)  *z4 

            ijspa(n1) = kofprev + 1 

c              h-neighbors for h(i,j):

            if (i.lt.nx) then 
              n2 = mhind(i+1,j)
              if (n2.ge.1) then
c               upstream:
                elspa(n1) = elspa(n1)   + max (ua(i,j),0.)*z2
                if (ua(i,j).lt.0.) then
                  kof = kof + 1
                  elspa(kof) = elspa(kof) + min (ua(i,j),0.)*z2
                  ijspa(kof) = n2
                endif
              endif
            endif

            if (i.gt.1) then 
              n2 = mhind(i-1,j)
              if (n2.ge.1) then
c               upstream:
                elspa(n1) = elspa(n1)   - min (ua(i-1,j),0.)*z1
                if (ua(i-1,j).gt.0.) then
                  kof = kof + 1
                  elspa(kof) = elspa(kof) - max (ua(i-1,j),0.)*z1
                  ijspa(kof) = n2
                endif
              endif
            endif
  
            if (j.lt.ny) then 
              n2 = mhind(i,j+1)
              if (n2.ge.1) then
c               upstream:
                elspa(n1) = elspa(n1)   + max (va(i,j),0.)*z4
                if (va(i,j).lt.0.) then
                  kof = kof + 1
                  elspa(kof) = elspa(kof) + min (va(i,j),0.)*z4
                  ijspa(kof) = n2
                endif
              endif
            endif
 
            if (j.gt.1) then
              n2 = mhind(i,j-1)
              if (n2.ge.1) then
c               upstream:
                elspa(n1) = elspa(n1)   - min (va(i,j-1),0.)*z3
                if (va(i,j-1).gt.0.) then
                  kof = kof + 1
                  elspa(kof) = elspa(kof) - max (va(i,j-1),0.)*z3
                  ijspa(kof) = n2
                endif
              endif
            endif
c<<<<<<<<<<<<<<
          endif    ! n1 >= 1
c<<<<<<<<<<<<<<

        enddo
      enddo

      if (kof.gt.nspamax) then
        write (6,*) 'Error (dosparsemelh): kof > nspamax. kof=', kof,
     *              ' nspamax=',nspamax
        stop
      endif 

      ijspa(nhtot+1) = kof + 1

      if (nhtot.ne.ijspa(1) - 2) then
        write (6,*) 'Error (dosparsemelh): nhtot, ijspa(1)-2 =',
     *              nhtot, ijspa(1)-2
        stop
      endif

cmat  call zero (asparse, nhmat*nhmat)
cmat  if (nhtot.gt.nhmat) then
cmat    write (6,*) 'nhmat too small: nhtot, nhmat=', nhtot, nhmat
cmat    stop
cmat  endif
cmat  do m=1,nhtot
cmat    asparse(m,m) = elspa(m) - 1. 
cmat    ka = ijspa(m)
cmat    kb = ijspa(m+1)
cmat    if (kb.gt.ka) then
cmat      do k=ka,kb-1
cmat        asparse(m,ijspa(k)) = elspa(k)
cmat      enddo
cmat    endif
cmat  enddo
cmat
cmat  do m=1,nhtot
cmat    i = indsav(m)
cmat    j = jndsav(m)
cmat    colsum(m) = 0.
cmat    do n=1,nhtot
cmat      ii = indsav(n)
cmat      jj = jndsav(n)
cmat      zar = darea(ii,jj)
cmat      colsum(m) = colsum(m) + asparse(n,m)*zar  
cmat    enddo
cmat  enddo
cmat
cmat  allcolsum = 0.
cmat  allcolmax = 0.
cmat  do m=1,nhtot
cmat    allcolsum = allcolsum + abs(colsum(m))
cmat    allcolmax = max(allcolmax,abs(colsum(m)))
cmat  enddo
cmat    
cmat  write (417,'(/a,f10.2, 2i6)') 'time=', timeice, iterc, iloopmel
cmat  write (417,'(a,2e14.6/)') 'allcolsum,max=',allcolsum/dd0**2,
cmat *                                           allcolmax/dd0**2 
cmat  do m=1,nhtot
cmat    write (417,'(i6,e14.6,2i6,e14.6)') m, colsum(m)/dd0**2, 
cmat *    indsav(m), jndsav(m), h0(indsav(m),jndsav(m))
cmat  enddo

c        Solve linear system arr*vec = rhs for vec.
c        If number of pts <= nuvsmall, use gaussian elimination 
c        (avoids sparse-solution problems for nhtot <= ~2)

c++++++++++++++++++++++++++++++++
      if (nhtot.gt.nuvsmall) then
c++++++++++++++++++++++++++++++++

#if defined (NUMRECMELH) || defined (SLAPMELH)
      if (firstvech) then
        call zero (vech, nhtot)
        firstvech = .false.
      endif
#endif
#if defined (NUMRECMELH)
      itol = 2       ! 4 ! 1
      tol = 1.e-11   ! 1.e-4 ! 1.e-6 ! 1.e-11
      itmax = 10000  ! 30000
      call linbcg (nhtot, rhs, vech, itol, tol, itmax, itlin, errlin)
c     write (ioterm,'(a,i5,a,e15.5)')
c    *  'linbcg(melh): itlin=',itlin,' errlin=',errlin
      if (abs(errlin).gt.tol) write (ioterm,'(a,f10.2,a,i5,a,e15.5)')
     *  '*** Warning linbcg(melh): time=',timeice,
     *  ' itlin=',itlin,' errlin=',errlin
#elif defined (SLAPMELH)
      itol = 2
      tol = 1.e-4
      itmax = 10000
      call slap (nhtot, rhs, vech, itol, tol, itmax, itlin, errlin)
c     write (ioterm,'(a,i5,a,e15.5)')
c    *  'slap(melh): itlin=',itlin,' errlin=',errlin
      if (abs(errlin).gt.tol) write (ioterm,'(a,i5,a,e15.5)')
     *  '*** Warning slap(melh): itlin=',itlin,' errlin=',errlin
#elif defined (SUPERLUMELH)
      call superlu (nhtot, rhs, vech, info)
      if (info.ne.0)write (ioterm,*)'*** Warn superlu(melh): info=',info
#elif defined (WATSONMELH)
      call watson (nhtot, rhs, vech, info)
      if (info.ne.0)write (ioterm,*) '*** Warn watson(melh): info=',info
#elif defined (MKLMELH)
      call mkl (nhtot, rhs, vech, info)
      if (info.ne.0) write (ioterm,*)'*** Warning mkl(melh): info=',info
#elif defined (GAUSSMELH)
      call gaussdo (nhtot, rhs, vech, info)
      if (info.ne.0) write (ioterm,*)'*** Warn gauss(melh): info=',info
#endif

c+++++++++
      else
c+++++++++
        call gaussdo (nhtot, rhs, vech, info)
        if (info.ne.0) write (ioterm,*)'*** Warn gauss(melh):info=',info
c++++++++++
      endif
c++++++++++

c        Set new h

      do j=1,ny
        do i=1,nx
          n1 = mhind(i,j)
          if (n1.ge.1) then
            h1(i,j) = vech(n1)
          endif
        enddo
      enddo

c#############
 1000 continue
c#############

      do j=1,ny
        do i=1,nx
          if (mhind(i,j).le.0) then
            if (maskmel(i,j).ne.0) then 
              h1(i,j) = h0(i,j)
              h1(i,j) = h1(i,j) + budgsup(i,j)*dt
              if (h1(i,j).gt.0.) then 
                h1(i,j) = h1(i,j) + budgall(i,j)*dt
              else
                budgall(i,j) = 0.
              endif
            else
              budgall(i,j) = 0.
              budgsup(i,j) = 0.
            endif
          endif
        enddo
      enddo

      call zero (budgneg, nx*ny)                                  ! diag
      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).ne.0) then
            if (h1(i,j).lt.0.) then
              budgneg(i,j) = - h1(i,j)/dt                         ! diag
              h1(i,j) = 0.
            endif
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine calv2mel (fedge, hedge, hice, calvice,
     *                     maskwater, maskmel, 
     *                     ucamel, vcamel, hucamel, hvcamel, supcalv)

c     Calculate supply of melange from calving (supcalv, from calvice).
c     If fedge < 1, all calving supply is local (same cell).   
c     If fedge = 1, no melange exists in cell, calving supplies 
c     adjacent ocean cells (with fedge<1), so first need to count
c     # of those adjacent cells (kount) and distribute calvice evenly.
c
c     nb: calvice can be > 0 for fedge=1 (ocean pt., and rarely for 
c     last grounded pt.) - see setfhedge, calving in iceocean.F.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  fedge(nx,ny), hedge(nx,ny), hice(nx,ny), calvice(nx,ny), 
     *  maskwater(nx,ny), maskmel(nx,ny), 
     *  ucamel(0:nxp,0:nyp),  vcamel(0:nxp,0:nyp),
     *  hucamel(0:nxp,0:nyp), hvcamel(0:nxp,0:nyp),
     *  supcalv(nx,ny)

      parameter (nlook = 4)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /

      call zero (supcalv, nx*ny)

      do j=1,ny
        do i=1,nx
          if (calvice(i,j).gt.0.) then                
c.....................................
            if (fedge(i,j).eq.1.) then                ! maskwater=0 or 1
c.....................................

c             no melange can exist (fedge=1). Count neighbors (ocean) 
c             with fedge < 1, that can receive this cell's calvice:

              kount = 0
              do look=1,4
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (.not.(ii.eq.i .and. jj.eq.j)) then    ! for flowline
c                 if (maskwater(ii,jj).eq.1.and.fedge(ii,jj).lt.1) then
                  if (maskmel(ii,jj).ne.0) then ! equiv
                    kount = kount + 1
                  endif
                endif
              enddo

c             repeat neighbor search, distribute calvice evenly:

              if (kount.ge.1) then
                do look=1,4
                  ii = max (1, min (nx, i + ish(look)))
                  jj = max (1, min (ny, j + jsh(look)))
                  if (.not.(ii.eq.i .and. jj.eq.j)) then  ! for flowline
c                   if (maskwater(ii,jj).eq.1 .and. fedge(ii,jj).lt.1)
                    if (maskmel(ii,jj).ne.0) then ! equiv
                      supcalv(ii,jj) = supcalv(ii,jj) 
     *                               + (calvice(i,j) / kount)
                    endif
                  endif
                enddo
              endif

c...............
            else
c...............

c             fedge < 1 (so must be ocean). calvice remains local:

              supcalv(i,j) = supcalv(i,j) + calvice(i,j)

c................
            endif
c................

          endif
        enddo
      enddo

c==================

      call zero (ucamel, (nxp+1)*(nyp+1))
      call zero (vcamel, (nxp+1)*(nyp+1))
      call zero (hucamel, (nxp+1)*(nyp+1))
      call zero (hvcamel, (nxp+1)*(nyp+1))

      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).ne.0 .and. supcalv(i,j).gt.0.) then 
            kount = 0
            zper= 0.
            zhc = 0.
            zc = 0.
            do look=1,4
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (.not.(ii.eq.i .and. jj.eq.j)) then    ! for flowline
                if (maskmel(ii,jj).eq.0 .and. hice(ii,jj).gt.0.) then
                  kount = kount + 1
                  if (ii.ne.i) zdd = dy(i,j)
                  if (jj.ne.j) zdd = dx(i,j)
                  zper = zper + zdd
                  zhc = zhc + hedge(ii,jj)*calvice(ii,jj)
                  zc = zc + calvice(ii,jj)
                endif
              endif
            enddo
            zhc = zhc + hedge(i,j)*calvice(i,j)
            zc = zc + calvice(i,j)

            if (kount.ge.1) then
c             zu = supcalv(i,j)*darea(i,j) / zaper
c             zh = zaper/zper
              zh = zhc / zc
              zu = supcalv(i,j)*darea(i,j) / (zh*zper)

              do look=1,4
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (.not.(ii.eq.i .and. jj.eq.j)) then  ! for flowline
                  if (maskmel(ii,jj).eq.0 .and. hice(ii,jj).gt.0.) then
                    if (ii.gt.i) then
                      ucamel(i,j) = -zu
                      hucamel(i,j) = zh
                    else if (ii.lt.i) then
                      ucamel(ii,j) = zu
                      hucamel(ii,j) = zh
                    else if (jj.gt.j) then
                      vcamel(i,j) = -zu
                      hvcamel(i,j) = zh
                    else if (jj.lt.j) then
                      vcamel(i,jj) = zu
                      hvcamel(i,jj) = zh
                    endif
                  endif
                endif
              enddo
            endif

          endif

        enddo
      enddo

c        Just return veloc magnitudes. +/- directions will be set in 
c        setvelbc

      do j=0,nyp
        do i=0,nxp
          ucamel(i,j) = abs(ucamel(i,j))
          vcamel(i,j) = abs(vcamel(i,j))
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine slumpsup (h0, maskmel, budgsup,
     *                     dt, timeice, iterc, iloopmel)

c        "Slump" new supply of melange onto existing melange surface 

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h0(nx,ny), h1(nx,ny), maskmel(nx,ny), budgsup(nx,ny)

      dimension h2(nx,ny), zslump(nx,ny)

      if (slopemaxmel.le.0.) return

      call scopy (nx*ny, h0, 1, h2, 1)

c=============================================================
      do j=1,ny
        do i=1,nx
c         if (maskmel(i,j).ne.0 .and. budgsup(i,j).gt.0.) then
          if (maskmel(i,j).ne.0) then
c=============================================================
            call zero (zslump, nx*ny)

            zdsup = budgsup(i,j)*dt 

            iwid = max (1, nint (500.e3/dx(i,j)))
            jwid = max (1, nint (500.e3/dy(i,j)))

c>>>>>>>>>>>>>>>>>>>>>>>>>>
            do iterbin=1,50
c>>>>>>>>>>>>>>>>>>>>>>>>>>
              if (iterbin.eq.1) then
                zremain = zdsup + h2(i,j)
              else if (iterbin.eq.2) then
                zremain = 0.
              else 
                zremain = 0.5*(zremain1+zremain2)
              endif

              zspace = 0.
              do jj = max(j-jwid,1), min(j+jwid,ny)
                do ii = max(i-iwid,1), min(i+iwid,nx)
                  if (ii.ne.i.and.jj.ne.j .and. maskmel(ii,jj).ne.0)then
                    zdist = iabs(i-ii)*dx(i,j) + iabs(j-jj)*dy(i,j)
                    zslump(ii,jj) = 
     *                       max ( zremain - zdist*slopemaxmel
     *                           - h2(ii,jj), 0. )
                    zspace = zspace + zslump(ii,jj) * darea(ii,jj)
                  endif
                enddo
              enddo

              zneed = (zdsup+h2(i,j) - zremain)*darea(i,j)
              zexcess = zspace - zneed

c             if (iterbin.eq.3) write (6,*)
c             if (iterbin.ge.3) write (6,'(a,f10.3,5i6,2f8.3,4f10.3)')
c    *          'slumpsup: ', 
c    *          timeice, iterc, iloopmel, iterbin, 
c    *          i,j, alatd(i,j), alond(i,j),
c    *          zremain1, zremain2, zremain, zexcess/(darea(i,j)*dt)

c===================================
              if (iterbin.eq.1) then 
c===================================
                if (zspace.le.0.) goto 6000   ! no slopemaxmel violation
                zremain1 = zremain
                zexcess1 = zexcess
c========================================
              else if (iterbin.eq.2) then
c===================================
                zremain2 = zremain
                zexcess2 = zexcess
c=================
              else 
c=================

c               almost converged, skip out (test vs. m/yr):
                if ( abs(zexcess/(darea(i,j)*dt)).lt..001) goto 5000

c..............................................................
                if (sign(1.,zexcess).eq.sign(1.,zexcess1) .and.
     *              sign(1.,zexcess).ne.sign(1.,zexcess2)) then
c..............................................................
                  zremain1 = zremain
                  zexcess1 = zexcess
c                 zremain2 = zremain2
c                 zexcess2 = zexcess2
c...................................................................
                else if (sign(1.,zexcess).ne.sign(1.,zexcess1) .and.
     *                   sign(1.,zexcess).eq.sign(1.,zexcess2)) then
c...................................................................
c                 zremain1 = zremain1
c                 zexcess1 = zexcess1
                  zremain2 = zremain 
                  zexcess2 = zexcess
c...................
                else
c...................
                  write (6,'(/a)') 'slumpsup: binary search error'
                  write (6,'(f10.3,5i6,2f8.3/
     *                       f16.3/2f16.3/2f16.3/4f16.3)')
     *              timeice, iterc, iloopmel, iterbin, i,j, 
     *              alatd(i,j), alond(i,j),
     *              zdsup,  
     *              zremain1, zexcess1,
     *              zremain2, zexcess2,
     *              zremain,  zexcess,  zspace, zneed

                  write (6,*)
                  zzremain = 0.        ! iterbin = 2
                  zzremain = zremain   ! iterbin = curr
                  zzneed = (zdsup-zzremain)*darea(i,j)
                  zzspace = 0.
                  do jj = max(j-jwid,1), min(j+jwid,ny)
                    do ii = max(i-iwid,1), min(i+iwid,nx)
                       if (ii.ne.i.and.jj.ne.j .and. 
     *                   maskmel(ii,jj).ne.0) then
                        zdist = iabs(i-ii)*dx(i,j) + iabs(j-jj)*dy(i,j)
                        zzslump = 
     *                    max ( h2(i,j) + zzremain - zdist*slopemaxmel
     *                        - h2(ii,jj), 0. )
                        zzspace = zzspace + zzslump * darea(ii,jj)
                        if (zzslump.gt.0.) 
     *                    write (6,'(2i6, 4f16.3)') ii, jj, 
     *                      zdist, 
     *                      h2(i,j)+zzremain-zdist*slopemaxmel,
     *                      h2(ii,jj), 
     *                      zzslump
                      endif
                    enddo
                  enddo
                  zzexcess = zzspace - zzneed
                  write (6,'(/4f16.3/32x,3f16.3)') 
     *              zzremain, zzexcess, zzspace, zzneed,
     *              zzspace/darea(i,j), zzneed/darea(i,j), h2(i,j)
                  stop
c....................
                endif
c....................

c==================
              endif
c==================
c>>>>>>>>>>>>>>>>
            enddo
c>>>>>>>>>>>>>>>>

 5000      continue

            do jj = max(j-jwid,1), min(j+jwid,ny)
              do ii = max(i-iwid,1), min(i+iwid,nx)
                if (ii.ne.i.and.jj.ne.j .and. maskmel(ii,jj).ne.0) then
                  h2(ii,jj) = h2(ii,jj) + zslump(ii,jj) 
                endif
              enddo
            enddo

 6000       continue

            h2(i,j) = zremain

c==============
          endif
        enddo
      enddo
c==============

      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).ne.0) then
            budgsup(i,j) = (h2(i,j) - h0(i,j)) / dt
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine adjustmel (maskmel, h, hs, hw, fmel,
     *                      hb, hice, fedge, maskwater, sealev,
     *                      budgadjmel, dt, ifadj)

c     Adjusts melange thickness, surface elevation, water depth, 
c     maskmel, fmel. Zeros melange for land and fedge=1.
c
c     nb: maskmel =  0 if not accessible to melange (land or ice),
c                 =  1 if floating melange (ocean, no ice), 
c                 = -1 if grounded melange (ocean, no ice).
c
c     nb: fmel is only 0<fmel<1 if fedge is < 1 (ice shelf).
c         fmel is not used in dynamics - no fractional treatment 
c         like shelves. After dynamics, melange can be considered
c         displaced (reduced) by any ice in cell (fmel = 1-fedge)
c         for diagnostics - but for conservation of melange mass,
c         imagine melange occupies entire cell (meltabbud in iceshow.F).
c
c     nb: Mask of maskwater=0 or fedge=1 remains same through the 
c         overall model timestep, after first iloopmel call to icemel,
c         so zeroing of previous melange only occurs for first call.
c         So for first iloopmel call only (ifadj=.true.), redistribute  
c         zero-ed melange into neighboring eligible cells, or if 
c         no such cells, track in budgadjmel.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  maskmel(nx,ny), h(nx,ny), hs(nx,ny), hw(nx,ny), fmel(nx,ny),
     *  hb(nx,ny), hice(nx,ny), fedge(nx,ny), maskwater(nx,ny),
     *  budgadjmel(nx,ny)

      logical ifadj

      dimension hprev(nx,ny)

      parameter (nlook = 4)
      dimension ish(nlook), jsh(nlook)
      save ish, jsh
      data ish /-1, 1,  0, 0 /
      data jsh / 0, 0, -1, 1 /

      call scopy (nx*ny, h, 1, hprev, 1)

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.0) then
            h(i,j) = 0.
            hw(i,j) = 0.
            hs(i,j) = hb(i,j) + hice(i,j)
            maskmel(i,j) = 0
            fmel(i,j) = 0.
          else if (fedge(i,j).eq.1.) then
            h(i,j) = 0.
            hw(i,j) = 0.
            hs(i,j) = sealev + hice(i,j)*(1.-rhor)
            maskmel(i,j) = 0
            fmel(i,j) = 0.
          else
            hw(i,j) = max (0., sealev - hb(i,j) - rhormel*h(i,j))
            hs(i,j) = hb(i,j) + hw(i,j) + h(i,j)
            if (hw(i,j).gt.0.) then 
              maskmel(i,j) = 1
            else
              maskmel(i,j) = -1
            endif
            if (h(i,j).gt.0.) then 
              fmel(i,j) = 1. - fedge(i,j)
            else
              fmel(i,j) = 0.
            endif
          endif
        enddo
      enddo

c        If first iloopmel call, redistribute zero-ed previous  
c        melange (due to maskwater=0 or fedge=1) equally into
c        melange-accessible cells (conserving mass), or track in 
c        budgadjmel if no such cells.

c--------------------
      if (ifadj) then
c--------------------
      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).eq.0 .and. h(i,j).ne.hprev(i,j)) then
            zadj = h(i,j) - hprev(i,j)                     ! zadj is -ve
            nnex = 0
            do look=1,4
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (.not.(ii.eq.i .and. jj.eq.j)) then      ! for flowline
                if (maskmel(ii,jj).ne.0) nnex = nnex + 1
              endif
            enddo
            if (nnex.gt.0) then
              do look=1,4
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (.not.(ii.eq.i .and. jj.eq.j)) then    ! for flowline
                  if (maskmel(ii,jj).ne.0) 
     *              h(ii,jj) = h(ii,jj) 
     *                         - (zadj*darea(i,j)/nnex)/darea(ii,jj)
                endif
              enddo
            else
              budgadjmel(i,j) = budgadjmel(i,j) + zadj/dt
            endif
          endif
        enddo
      enddo
c----------
      endif
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine smoothermel (h, maskmel)

c     Apply 9-pt smoother, to avoid alternating-grid
c     "ripples" in hmel and +/-/+/- velocity noise. Conserves
c     total melange mass.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny), maskmel(nx,ny)

      dimension 
     *  hold(nx,ny)

      dimension weis(-1:1,-1:1)
      save weis
      data weis /.0625, .125, .0625,
     *           .125,  .25,  .125,
     *           .0625, .125, .0625 /

      call scopy (nx*ny, h, 1, hold, 1)

      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).ne.0) then
            zh = 0.
            do jinc=-1,1
              jj = max (1, min(ny, j+jinc))
              do iinc=-1,1
                ii = max (1, min(nx, i+iinc))
                if (maskmel(ii,jj).ne.0) then  
                  zh = zh + hold(ii,jj)*darea(ii,jj)*weis(iinc,jinc)
                else
                  zh = zh + hold(i,j)*darea(i,j)*weis(iinc,jinc)
                endif
              enddo
            enddo
            h(i,j) = zh/darea(i,j)
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine smoothermelu (u, maskmel, igrid)

c     Apply 9-pt smoother, to avoid alternating-grid "waves" in 
c     [u,v]amel and hmel. igrid=1 for u-grid, 2 for v-grid.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  u(0:nxp,0:nyp), maskmel(nx,ny)

      dimension 
     *  uold(0:nxp,0:nyp), masku(0:nxp,0:nyp)

c     conserves total "volume", and flattens both checker-grid noise 
c     and alternating strips:
      dimension weis(-1:1,-1:1)
      save weis
      data weis /.0625, .125, .0625,
     *           .125,  .25,  .125,
     *           .0625, .125, .0625 /

      call scopy ((nxp+1)*(nyp+1), u, 1, uold, 1)

      call izero (masku, (nxp+1)*(nyp+1))
      do j=1,ny
        do i=1,nx
          if (igrid.eq.1) then
            ii = min (i+1, nx)
            jj = j
          else if (igrid.eq.2) then
            ii = i
            jj = min (j+1, ny)
          endif
          if (maskmel(i,j).ne.0 .and. maskmel(ii,jj).ne.0) masku(i,j)=1
        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if (masku(i,j).ne.0) then
            zu = 0.
            do jinc=-1,1
              jj = max (0, min(nyp, j+jinc))
              do iinc=-1,1
                ii = max (0, min(nxp, i+iinc))
                if (masku(ii,jj).ne.0) then  
                  zu = zu + uold(ii,jj)*weis(iinc,jinc)
                else
                  zu = zu + uold(i,j)  *weis(iinc,jinc)
                endif
              enddo
            enddo
            u(i,j) = zu
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine thetacalcmel (h, fedge, maskmel, arhap, 
     *                         ua, va, thetau, thetav)

c     Sets theta[u,v] = 1 - fraction of g.l. longitudinal stress 
c     buttressed, due to melange. Like thetacalc in icedyn.F.
c
c     Nb: theta[u,v] are on h-grid, for E-W and N-S directions 
c     respectively, and melange thickness is h, so that if all melange 
c     is free-floating, theta=1 exactly (property of elliptic solution).
c

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),
     *  maskmel(nx,ny),        arhap(nx,ny),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  thetau(nx,ny),         thetav(nx,ny)

      parameter (hthetmin = .01)   ! have tried 0.1, 1., 10.

#include <comicepresmel.h>

      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).ne.0) then

            if (h(i,j).lt.hthetmin) then
              thetau(i,j) = 1.
              thetav(i,j) = 1.
            else
              zdudx = (ua(i,j)-ua(i-1,j)) / dx(i,j)
              zdvdy = (va(i,j)-va(i,j-1)) / dy(i,j)
              if (zdudx.gt.0.) then
                zcoefu = coefdivmel
              else
                zcoefu = coefconmel
              endif
              if (zdvdy.gt.0.) then
                zcoefv = coefdivmel
              else
                zcoefv = coefconmel
              endif
              zg1 =  max (zdudx**2 + zdvdy**2 + zdudx*zdvdy, gmin)
              zvisc = 0.5 / (zg1**powivmel)

c             thetau(i,j)  = ( 4.*zcoefu*zvisc*2.*zdudx/arhap(i,j) )
c    *                       / (rhomp*grav*h(i,j))
c             thetav(i,j)  = ( 4.*zcoefv*zvisc*2.*zdvdy/arhap(i,j) )
c    *                       / (rhomp*grav*h(i,j))

              thetau(i,j) = 
     *           ( 2.*zcoefu*zvisc*(2.*zdudx+zdvdy)*h(i,j)/arhap(i,j) )
     *           / ( 0.5*rhomp*grav*(h(i,j)**2) + presmel(h(i,j)) )
              thetav(i,j) = 
     *           ( 2.*zcoefv*zvisc*(2.*zdvdy+zdudx)*h(i,j)/arhap(i,j) )
     *           / ( 0.5*rhomp*grav*(h(i,j)**2) + presmel(h(i,j)) )

              thetau(i,j) = min (1., thetau(i,j)) ! can be -ve, converge
              thetav(i,j) = min (1., thetav(i,j)) ! can be -ve, converge
            endif

          else

            thetau(i,j) = -1.
            thetav(i,j) = -1.

          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine checkmelforce (h, hs, maskmel,
     *                          thetau, thetav, coef1, arhap,
     *                          ua, va, muind, mvind,
     *                          coefbu, coefbv,
     *                          hice, fedge,
     *                          timeice, dt, 
     *                          iterc, iloopmel, nloopmel, nyeartabmel)

c     Diagnostic check of domain-wide force balance of mel SSA eqns.
c     Basically an areal integral of x and y eqns separately,
c     using b.c. balance at open ocean. 
c     Computes cur[side,base,back]mel[x,y] (vector) quantities,
c     and tests balance between them.
c     (see notes 10/18/17 and paper Appendix).
c     Must be called right after A-loop (veloc solutions),
c     before h is changed further down in C-loop.
c     Writes tabular results to iutabmelforce, every nyeartabmel years,
c     (same freq as in meltabbud). 
c
c     Similar to checkmelforce2, but here called within C-loop
c     for exact correspondence between h and ua,va.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),              hs(nx,ny),          maskmel(nx,ny),
     *  thetau(nx,ny),         thetav(nx,ny),      
     *  coef1(nx,ny),          arhap(nx,ny),
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),    mvind(0:nxp,0:nyp),
     *  coefbu(0:nxp,0:nyp),   coefbv(0:nxp,0:nyp),
     *  hice(nx,ny),           fedge(nx,ny)

      character*1 cgrmel

      save kountw
      data kountw /0/

#include <comicepresmel.h>

c        First update theta[u,v] (for mel)

      call thetacalcmel (h, fedge, maskmel, arhap, 
     *                   ua, va, thetau, thetav)

      cursidemelx = 0.
      cursidemely = 0.
      curbasemelx = 0.
      curbasemely = 0.
      curbackmelx = 0.
      curbackmely = 0.

c        Compute domain-wide side drag, diagnostic only (Newtons)
    
      do j=1,ny
        jm = max (j-1,1) 
        jp = min (j+1,ny)
        do i=1,nx
          im = max (i-1,1) 
          ip = min (i+1,nx)

          if ( (h(i,j).gt.0 .or. h(ip,j).gt.0.) .and. 
     *        muind(i,j).gt.0 ) then
            zside = 0.
#if (TESTMEL2D)
            if (nx.gt.ny .and. (j.eq.1 .or. j.eq.ny)) zside = 1.
#else
            if (j.gt.1) then 
              if (maskmel(i, j-1).eq.0) zside = zside + 0.5
              if (maskmel(ip,j-1).eq.0) zside = zside + 0.5
            endif
            if (j.lt.ny) then 
              if (maskmel(i, j+1).eq.0) zside = zside + 0.5
              if (maskmel(ip,j+1).eq.0) zside = zside + 0.5
            endif
#endif
            zdudy = ua(i,j)/(0.5*dyu(i,j))
            cursidemelx = cursidemelx 
     *        + zside*sidedragmel * 0.5*(0.5*(coef1(ip,j)+coef1(i,j)))
     *          * zdudy * dxu(i,j)
          endif
    
          if ( (h(i,j).gt.0 .or. h(i,jp).gt.0.) .and. 
     *          mvind(i,j).gt.0 ) then
            zside = 0.
#if (TESTMEL2D)
            if (ny.gt.nx .and. (i.eq.1 .or. i.eq.nx)) zside = 1.
#else
            if (i.gt.1) then 
              if (maskmel(i-1,j ).eq.0) zside = zside + 0.5
              if (maskmel(i-1,jp).eq.0) zside = zside + 0.5
            endif
            if (i.lt.nx) then 
              if (maskmel(i+1,j ).eq.0) zside = zside + 0.5
              if (maskmel(i+1,jp).eq.0) zside = zside + 0.5
            endif
#endif
            zdvdx = va(i,j)/(0.5*dxv(i,j))
            cursidemely = cursidemely 
     *        + zside*sidedragmel * 0.5*(0.5*(coef1(i,jp)+coef1(i,j)))
     *          * zdvdx * dyv(i,j)
          endif

        enddo   
      enddo   
    
c        Compute domain-wide basal drag, diagnostic only (Newtons)
    
      do j=1,ny
        jp = min (j+1, ny)
        do i=1,nx
          ip = min (i+1, nx)

          if ( (h(i,j).gt.0 .or. h(ip,j).gt.0.) .and.
     *         muind(i,j).gt.0  )
     *      curbasemelx = curbasemelx
     *                  + coefbu(i,j)*ua(i,j)*dxu(i,j)*dy(i,j)

          if ( (h(i,j).gt.0 .or. h(i,jp).gt.0.) .and.
     *         mvind(i,j).gt.0 )
     *      curbasemely = curbasemely
     *                  + coefbv(i,j)*va(i,j)*dyv(i,j)*dx(i,j)
        enddo
      enddo

c        Compute domain-wide back stress of melange on all ice faces:
c        (i) cliff failure (theta[u,v], in schoofgl), and 
c        (ii) dynamical  b.c. for regular calving faces (bp[u,v]mel, 
c             in dosparse, etc). 
c        See notes ("delta",10/18-26/17). Diagnostic only (Newtons).

      do j=1,ny
        jm = max (j-1,1) 
        jp = min (j+1,ny)
        do i=1,nx
          im = max (i-1,1) 
          ip = min (i+1,nx)
c...............................
          if (h(i,j).gt.0.) then
c...............................

c===
c x:
c===
            if (maskmel(im,j).eq.0.and.maskmel(ip,j).ne.0) then
              idir = 1
            else if (maskmel(ip,j).eq.0.and.maskmel(im,j).ne.0) then
              idir = -1
            else
              idir = 0
            endif

            curbackmelx = curbackmelx 
     *        + ( (1.-thetau(i,j))
     *            * (0.5*rhomp*grav*(h(i,j)**2) + presmel(h(i,j)))
     *          ) *  dy(i,j) * idir

c===
c y:
c===
            if (maskmel(i,jm).eq.0.and.maskmel(i,jp).ne.0) then
              idir = 1
            else if (maskmel(i,jp).eq.0.and.maskmel(i,jm).ne.0) then
              idir = -1
            else
              idir = 0
            endif

            curbackmely = curbackmely 
     *        + (   (1.-thetav(i,j))
     *            * (0.5*rhomp*grav*(h(i,j)**2) + presmel(h(i,j)))
     *          ) *  dx(i,j) * idir

c..............
          endif   
c..............
        enddo   
      enddo   

      curbackerrx = curbackmelx - (cursidemelx + curbasemelx)
      curbackerry = curbackmely - (cursidemely + curbasemely)

      cgrmel = ' '
      do j=1,ny
        do i=1,nx
          if (maskmel(i,j).eq.-1) cgrmel = '*'
        enddo
      enddo

c---
c     as in meltabbud, iceshow.F:
      if (nyeartabmel.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartabmel),dt)).lt.dt
     *                        .or. nyeartabmel.eq.-1 )
c    *       .or. iflast.eq.1 .or. nyeartabmel.eq.-1 )
     *   ) return  
c---

c        Write output 

      iu = iutabmelforce

      kountw = kountw + 1

      if (mod(kountw,50).eq.1) then
        write (iu,'(/a, 4a14,1x, 4a14,1x)')
     *    '        time ilmel iterc',

     *    '      sidex(N)',
     *    '      basex(N)',
     *    '      backx(N)',
     *    '   backerrx(N)',

     *    '      sidey(N)',
     *    '      basey(N)',
     *    '      backy(N)',
     *    '   backerry(N)'
      endif

      write (iu,'(f12.2,2i6, 4e14.4,a1, 4e14.4,a1)')
     *  timeice, iloopmel, iterc,

     *  cursidemelx,
     *  curbasemelx,
     *  curbackmelx,
     *  curbackerrx, 
     *  cgrmel,

     *  cursidemely,
     *  curbasemely,
     *  curbackmely,
     *  curbackerry, 
     *  cgrmel

      call flush (iu)

      return
      end

c-----------------------------------------------------------------------

      subroutine checkmelforce2 (curbackmel, curthetmel,
     *                           curhmel, curlmel,
     *                           h, hs, maskmel,
     *                           thetau, thetav,
     *                           ua, va, 
     *                           hice, fedge, hedge,
     *                           sealev, timeice, dt) 

c     Diagnostic calculation of domain-wide curbackmel, scalar (not 
c     vector) total |force difference| on all ice faces from that 
c     with no melange (see notes 10/18/17 and paper Appendix).
c     Also returns curthetmel, curhmel, curlmel, averages over
c     interfaces with ice. (All passed to meltabbud). 
c
c     Similar to checkmelforce. Don't need to call thetacalcmel
c     here, because called outside loop C, after call in icemel.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),              hs(nx,ny),          maskmel(nx,ny),
     *  thetau(nx,ny),         thetav(nx,ny),      
     *  ua(0:nxp,0:nyp),       va(0:nxp,0:nyp),
     *  hice(nx,ny),           fedge(nx,ny),       hedge(nx,ny)

#include <comicepresmel.h>

c        Compute domain-wide back stress of melange on all ice faces:
c        (i) cliff failure (theta[u,v], in schoofgl), and 
c        (ii) dynamical  b.c. for regular calving faces (bp[u,v]mel, 
c             in dosparse, etc). 
c        See notes ("delta",10/18-26/17). Diagnostic only (Newtons).

      do j=1,ny
        jm = max (j-1,1) 
        jp = min (j+1,ny)
        do i=1,nx
          im = max (i-1,1) 
          ip = min (i+1,nx)
c...............................
          if (h(i,j).gt.0.) then
c...............................

c           Modification if melange grounded next to ice face
c           (see notes "delta" 11/1/17):
c           if (maskmel(i,j).eq.-1) then
c             zhw = max (0., sealev-(hs(i,j)-h(i,j)))
c             zginc = 1.
c           else
c             zhw = 0.
c             zginc = 0.
c           endif

c===
c x:
c===
            if (maskmel(im,j).eq.0.and.maskmel(ip,j).ne.0) then
              idir = 1
            else if (maskmel(ip,j).eq.0.and.maskmel(im,j).ne.0) then
              idir = -1
            else
              idir = 0
            endif

            if ( (fedge(im,j).eq.1..and.hice(im,j).gt.0.) .or. 
     *           (fedge(ip,j).eq.1..and.hice(ip,j).gt.0.) ) then
              if (hice(i,j).gt.0.) then 
                zhice = hedge(i,j)
              else
                zhice = max (hedge(im,j), hedge(ip,j))
              endif
              zh = min (h(i,j), zhice)
c             zh = h(i,j)
              dcurbackmelx = 
     *               + (    max (0., min (1., 1.-thetau(i,j)))
     *                    * (0.5*rhomp*grav*(zh**2) + presmel(zh))
c    *                    + zginc * (  0.5*rhomel*rhormel*grav*(zh**2)
c    *                               - 0.5*rholiq        *grav*(zhw**2))
     *                 ) * dy(i,j) * idir
              curthetmel = curthetmel
     *               +                  thetau(i,j)  *dx(i,j)*iabs(idir)
              curhmel = curhmel + zh*dy(i,j)*iabs(idir)
              curlmel = curlmel + dy(i,j)*iabs(idir)
            else
              dcurbackmelx = 0.
            endif

c===
c y:
c===
            if (maskmel(i,jm).eq.0.and.maskmel(i,jp).ne.0) then
              idir = 1
            else if (maskmel(i,jp).eq.0.and.maskmel(i,jm).ne.0) then
              idir = -1
            else
              idir = 0
            endif

            if ( (fedge(i,jm).eq.1..and.hice(i,jm).gt.0.) .or. 
     *           (fedge(i,jp).eq.1..and.hice(i,jp).gt.0.) ) then
              if (hice(i,j).gt.0.) then 
                zhice = hedge(i,j)
              else
                zhice = max (hedge(i,jm), hedge(i,jp))
              endif
              zh = min (h(i,j), zhice)
c             zh = h(i,j)
              dcurbackmely = 
     *               + (    max (0., min (1., 1.-thetav(i,j)))
     *                    * (0.5*rhomp*grav*(zh**2) + presmel(zh))
c    *                    + zginc * (  0.5*rhomel*rhormel*grav*(zh**2)
c    *                               - 0.5*rholiq        *grav*(zhw**2))
     *                 ) * dx(i,j) * idir
              curthetmel = curthetmel
     *               +                  thetav(i,j)  *dy(i,j)*iabs(idir)
              curhmel = curhmel + zh*dx(i,j)*iabs(idir)
              curlmel = curlmel + dx(i,j)*iabs(idir)
            else
              dcurbackmely = 0.
            endif
c===

            curbackmel = curbackmel
     *                  + sqrt (dcurbackmelx**2 +dcurbackmely**2)

c..............
          endif   
c..............
        enddo   
      enddo   

      return
      end

c-----------------------------------------------------------------------

c||||
#else
c||||

      subroutine dummy_melange
      return
      end

c|||||
#endif
c|||||
