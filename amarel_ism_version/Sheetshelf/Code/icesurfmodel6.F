c=======================
#if defined (SURFMODEL6)
c=======================

c=======================================================================

      subroutine surfmodel6 
     *                    (ta, pa, ea, qa, mmax,
     *                     budgsnow,budgrain,budgmelt,budgevap,
     *                     tsurf, tsurfi,
     *                     hs, ticein, zicein, sdt, 
     *                     timeice, nyearfirn, iloc, jloc)

c     Surface snow-firn model, with nlay lagrangian layers.
c     Each layers tracks a lagrangian parcel of snow+liquid,
c     and liquid can percolate down, refreeze, or run off laterally. 
c     Does not track the channel/lake/moulin liquid system after that.

c     For each layer, snowm(k),sliqm(k) are column masses of snow
c     and liquid, with temperature tfirn(k). Snow porosity (air + liquid
c     space) is related to snow density (poros = 1 - rhosno/rhoice). 
c     Within snow fraction (1-poros), snow is essentially solid ice
c     with density rhoice. Maximum liquid holding capacity of porosity
c     space is fraction fret.

c     Seasonal cycles are run through each year (params dtimed,nspy are 
c     in comgrid.h/comparams.h), and  model is run through nyearfirn 
c     years (passed).

c     Atmos melt is computed not as pdd, but as a linearized surface 
c     energy flux + vertical conduction equation (tridiagonal), 
c     including insol (qa). Diurnal cycles of insol can be handled if 
c     dtimed < 1 day, and if so, diurnal cycles of temperature 
c     (sinusoid) can be superimposed (in subr surfbudg, via 
c     -DDTIMED (comgrid/comparams.h -> dtimed) and -DTADIURN).

c     Evaporation (passed as total evap from the climate model)
c     is taken from top layer, from snow+liquid in proportion to 
c     their masses.

c     Snow density (rhosno*) and grain size (rgrain*) of each layer are
c     time-stepped prognostically, following Arthern (2010) eqns
c     (nb: for dry snow, not for wet!), due to local processes, and 
c     mixed with properties of added snowfall.

c     Vertical temperature eqn is stepped time implicitly, with 2 passes
c     of the tridiag solver. For the 2nd pass, l.h.s. layer temperatures
c     are set to tmelt if they exceeded tmelt in the 1st pass. 

c     Liquid movement is handled by section "Adjustments":
c     One downward pass (1 to nlay), propagating excess heat (zheat)
c     and liquid from above (perccur). 
c     For impermeable layers (poros=0):
c       zheat warms/melts snow(ice), perccur runs off. 
c     For permeable layers (poros>0), 
c       perccur -> liquid stored.
c       zheat warms/melt some snow.
c       if t < tmelt and liquid, freeze liquid.
c       excess liquid (h*poros*fret) -> next layer's perccur.
c     If reached base of a permeable zone, i.e., contiguous layers
c     with poros > 0, below which is an impermeable layer or base
c     of model (nlay), then do a reverse pass upwards through the
c     permeable zone:
c       bottom perccur -> percrev.  
c       percrev -> liquid stored.
c       lateral runoff (stmt fn runlat) reduces liquid stored -> runoff.
c       excess liquid (h*poros, saturation) -> next (above) percrev.
c       any percrev at top -> runoff.

c     In a "layer shifting" section at the end of each time step,
c     layers are processed by:
c     c(1) remove zero-snow layers. 
c     c(2) maintain total col. mass by adding or removing from  bottom. 
c     c(3) merge very thin layers.
c     c(4) iteratively merge or split layers, one at a time, to restore 
c          thickness back towards target thicknesses vs. depth, and 
c          number of layers to nlay.   
c     Or: 
c     Option -DVADAPT for Eulerian layers, regridded each timestep.

c     Prognostic variables, saved between calls (in 2-D arrays a*):
c     rhosno*, rgrain*, snowm*, sliqm*, t* for layers 1:nlay.

c     Passed: 
c     ta, pa, ea, qa, mmax = climate (annual cycle, 1:nspy)
c     hs = surface elevation above sea level (m)
c     tice(in) = ice-model's temp at base of snow-firn (K)
c     sdt = std dev of air temp. for degree-day, deg C 
c           (only used for rainfall vs snowfall split)

c     Returned: 
c     budgsnow is snowfall, budgrain is rain, budgmelt is 
c     overal surface-model mobile melt, budgevap is evaporation.

c     Papers: 
c     Arthern et al., 2010, JGR.
c     Ligtenberg et al, 2011, The Cryo.
c     Langen et al., 2015, J. Clim.
c     Reijmer et al., 2012, The Cryo.
c     Rogozhina and Rau, 2014, The Cryo.
c     Robinson et al., 2010, The Cryo.
c     Janssens and Huybrechts, 2000, Ann Glac.
c     Reeh, 1989, Polarforschung.
c     Machguth et al., 2016, Nature CC (ice lenses).
c     Meyer et al., 2017, The Cryo (snow-firn model).
c     Linow et al., 2012, J. Glac. (grain sizes, cold snow)
c     Bougamont et al., 2005, JGR (internal runoff time scale)
c     Christianson et al., 2015, GRL (internal runoff time scale)

#if defined (STANDALONE6)
#include <comparams.h>
#else
#include <comicephys.h>
#include <comicegrid.h>
#endif

      dimension ta(nspy), pa(nspy), ea(nspy), qa(nspy)

c===
c     albsnodry = dry snow surface albedo (0.8)
c     albsnowet = wet snow surface albedo (0.6)
c     albsnoden = dense snow surface albedo (0.5)
c     acoeff = coefficent * air temp, W m-2 K-1 (10.)
c     atrans = atmos. tramsissivity with clouds (0.5) (=1 w/o clouds)
c     fcloud = cloudiness fraction 
c     aconst = constant in atmos flux eqn, W m-2 (-50.)
c     fret = max liq retention fraction of pore space 
c            (.02, Langen p.3696; .04 to .09, Ligtenberg eq.3)

      parameter (albsnodry = 0.8)
      parameter (albsnowet = 0.6)
      parameter (albsnoden = 0.5)
      parameter (atrans = 0.5)
      parameter (fcloud = 0.5)

#if defined (ACOEFFVAL)
      parameter (acoeff = ACOEFFVAL )
#else
#  if defined (MEYERHEWITT)  
      parameter (acoeff = 14.8)
#  else
      parameter (acoeff = 10.)
#  endif
#endif

#if defined (ACONSTVAL)
      parameter (aconst = ACONSTVAL )
#else
      parameter (aconst = -50.)
#endif
c===
c     thermal conductivity of snow, as in lsx
c     (Loth,1993,JGR,98,10451,eq.5; cf. Williams+Smith,1989,p.110)
c     condsno  = condsnoa + condsnob*(rhosno**2)
      parameter (condsnoa = .02*31556926)                ! J/a/m/K 
      parameter (condsnob = (2.5e-6)*31556926)           ! J/a/m/K/rho^2
c===
      parameter (rhoh2o    = 1000.) ! kg/m3, fresh water
#if defined (MEYERHEWITT)  
      parameter (rhosnofr  = rhoice*(1.-.64))
#else
      parameter (rhosnofr  = 350.)  ! cf. Ligtenberg, 2011, eq.2
#endif
      parameter (rhosnomin = .001)
c===
c     from Arthern et al. JGR 2010:
c     snow density equation:
      parameter (cra1 = 9.2e-9)    ! kg-1 m3 s, for rhosno < 550
      parameter (cra2 = 3.7e-9)    ! kg-1 m3 s, for rhosno > 550
      parameter (crb = 60.e3)      ! J mol-1
c     grain size equation:
      parameter (crr = 8.3144598)  ! gas constant, J mol-1 K-1
      parameter (cga = 1.3e-7)     ! m2 s-1
      parameter (cgb = 42.4e3)     ! J mol-1
c     see Linow et al. J.Glac 2012:
      parameter (rgrafr = .0002)   ! m  fresh snow, also min (was.0005)
      parameter (rgraic = .003)    ! m  ice, also max (was .005)
      parameter (rgrafr2 = rgrafr**2)
      parameter (rgraic2 = rgraic**2)
c===
      logical ift1, ift2
c===
#if defined (IFPROG6)
      logical firstsurfmodel6 (nx,ny)
      parameter (nxny=nx*ny)
      data firstsurfmodel6 /nxny*.true./
      save firstsurfmodel6
#endif
      logical firstdiag5
      save firstdiag5
      data firstdiag5 /.true./
c===
#if defined (MEYERHEWITT)  
c     parameter (totmass = 15000.)  ! kg/m2, corr. ~20 m
      parameter (totmass = 150000.) ! kg/m2
#else
      parameter (totmass = 150000.) ! kg/m2, corr. ~200-250 m (firnmice)
#endif

#if defined (ENSVALNL)
      parameter (nlay = ENSVALNL)
#else
      parameter (nlay = 5)
c     parameter (nlay = 10)
c     parameter (nlay = 30)
c     parameter (nlay = 100)
#endif

      parameter (nlaymax = nlay + 20)
      parameter (dmtop = 200.)                                   ! kg/m2
      parameter (dmmin = 10.)                                    ! kg/m2
      parameter (dtotmass = totmass/nlay)                        ! kg/m2

c     for time-step stability, cheat and cond:
      parameter (dzmin = .005)                                   ! m

      parameter (nskip = max (1,nlay/5))        ! for iudiag* output (k)
c     parameter (nskip = 1)
c===
c     for prognostic model variables, 2 horiz dims, saved between calls:
c     (but not on restart file, yet)
#if defined (IFPROG6)
      dimension
     *  asnowm(nx,ny,nlay),
     *  asliqm(nx,ny,nlay),
     *  atfirn(nx,ny,nlay),
     *  arhosno(nx,ny,nlay),
     *  argrain(nx,ny,nlay)
      save asnowm, asliqm, atfirn, arhosno, argrain
#endif
c===
c     prognostic: 
      dimension
     *  snowm(nlaymax), sliqm(nlaymax),
     *  tfirn(nlaymax), rhosno(nlaymax), rgrain(nlaymax)

c     diagnostic/temporary: 
      dimension
     *  thick(nlaymax), depth(nlaymax),
     *  fliq(nlaymax), poros(nlaymax), 
     *  perccur(0:nlaymax), percrev(0:nlaymax), 
     *  runofcur(nlaymax), runofzer(nlaymax),
     *  cheat(nlaymax), conds(nlaymax), cond(nlaymax),
     *  tfirn2(nlaymax), firnmass(nlaymax), firnmid(nlaymax),
     *  atri(nlaymax), btri(nlaymax), ctri(nlaymax), rhs(nlaymax)
c===
c     for target layer initialization, and layer shifting:
c     all non dim (0-1) w.r.t. column mass:
      save dmtarg, zmtarg, weitarg
      dimension dmtarg(nlay), zmtarg(nlay), weitarg(nlay)
      dimension dmlay(nlaymax), zmlay(nlaymax)
      parameter (nfine=10000, dfine=1./nfine)
      dimension dmfine(nfine)

#if defined (VADAPT)
      dimension dmvadapt(nlay)
      save dmvadapt
#endif

      parameter (score_ok = 1.5, zalph = 0.5, niter5max = 20)
      logical diagshift        ! for layer-shifting diagnostic dump info

      character*80 cform

c     for output to plot cross-sections vz z,day (gp, load.crosssec):

#if defined (CROSSSEC)
#  if defined (CROSSSEC2YR)
      parameter (nspycross = 2*nspy)
#  else
      parameter (nspycross = nspy)
#  endif
      dimension 
     *  crossz   (nlay,nspycross), crosst   (nlay,nspycross),
     *  crosspor (nlay,nspycross), crosssat (nlay,nspycross),  
     *  crosszvel(nlay,nspycross), crosswvel(nlay,nspycross),
     *  crossd (nspycross),
     *  zvel(nlay), zvelprev(nlay), wvel(nlay) 
#endif
     
c===
c     statement function for max capillary liquid capacity (fret,
c     max fraction of pore space, optional dep. on porosity):
#if defined (FRETVALA)
      parameter (fretvala = FRETVALA) 
#else
      parameter (fretvala = .04) 
#endif
#if defined (FRETVALB)
      parameter (fretvalb = FRETVALB) 
#else
      parameter (fretvalb = 0.) ! ~indep of poros, Coleou, Ann.Glac,1998
#endif
      fret(po) = fretvala + fretvalb*(1.-po)                ! po = poros
c===
c     statement function for internal (lateral) runoff rate, one layer
c     (applied in upward pass of "Adjustments" below) (kg m-2 s-1):
c
      parameter (trunlat = 100.*365.*86400.) ! 100-yr timescale(Knut)(s)
c     parameter (trunlat = 10.*86400.) ! much higher ? (K.Munneke, 2014)
      runlat (sl,slm,po) =           ! sl=sliqm, slm=saturated, po=poros
     *    ( po**1 )
     *  * ( (min (1., max (0.,sl/max(slm,.001)))) **1 )
     *  * ( sl/trunlat )
c     runlat (sl,slm,po) = 0.                       ! no internal runoff
c===

c        For diagnostic output (to fort.<iudiag*>, or 0 for no output)

      iudiag1 = 0
      iudiag2 = 0
      iudiag3 = 0
      iudiag4 = 0
      iudiag5 = 0
      iudiag6 = 0
      iudiagsh = 0
#if defined (STANDALONE6)
      iudiag1 = 100
      iudiag2 = 200
      iudiag3 = 300
      iudiag4 = 400
      iudiag5 = 500
#  if defined (CROSSSEC)
      iudiag6 = 600
#  endif
      iudiagsh = 700
#else
      if ( (alatd(iloc,jloc).ge. 79.4.and.alatd(iloc,jloc).le. 79.6.and.
     *      alond(iloc,jloc).ge.-31.0.and.alond(iloc,jloc).le.-18.0) 
     *     .and. iuz.le.9 ) then
        iuz = nint ((-18.-alond(iloc,jloc))/dlond + .501)
        iudiag1 = 100 + iuz
        iudiag2 = 200 + iuz
        iudiag3 = 300 + iuz
        iudiag4 = 400 + iuz
        iudiag5 = 500 + iuz
#  if defined (CROSSSEC)
        iudiag6 = 600 + iuz
#  endif
        iudiagsh = 700 + iuz
      endif
#endif

      tice = min (ticein, tmelt)           ! bottom b.c., from ice model
      zice = max (zicein, 1.)              ! bottom b.c., from ice model
#if defined (MEYERHEWITT)
      tice = tmelt - 10.
      zice = 10.
#endif

c        Initialize prognostic model variables. For all calls if
c        equilibrated, for first call each pt. if prognostic (IFPROG6).

#if defined (IFPROG6)
c,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
      if (firstsurfmodel6(iloc,jloc)) then
c,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
#endif

        zt = 0.
        do m=1,nspy
          zt = zt + min (ta(m),tmelt)/nspy ! K
        enddo

c          Depth-varying initial layer col. mass (snowm)

c          Set target layer thicknesses (dmtarg), 0-1 normalized by
c          totmass, from dmtop/totmass at top to slightly more than 
c          1./nlay at bottom. These two thicknesses are weighted 
c          according to weitarg (fn of depth), and then the 1./nlay
c          component is re-scaled (x zscal) so that total column mass 
c          = 1 (= totmass/totmass).

        ztota = 0.
        ztotb = 0.
        do k=1,nlay
          weitarg(k) = ((k-1.)/(nlay-1.))**(2.)
          ztota = ztota + (dmtop/totmass)*(1.-weitarg(k))
          ztotb = ztotb + (1./nlay)*weitarg(k)
        enddo

        zscal = (1. - ztota) / ztotb
        do k=1,nlay
          dmtarg(k) = (dmtop/totmass)*(1.-weitarg(k)) 
     *              + zscal*(1./nlay)*    weitarg(k)
        enddo

c          For use in layer shifting:
c          Set dmfine(1:nfine) on very fine regular vertical grid 
c          (0 to 1, non dim by totmass), to the interpolated values 
c          of target layer non-dim thicknesses dmtarg. For convenience,
c          first set zmtarg(1:nlay), the non-dim depths at the 
c          mid-points of the nlay target layers.

        zmtarg(1) = 0.5*dmtarg(1)
        do k=2,nlay
          zmtarg(k) = zmtarg(k-1) + 0.5*(dmtarg(k-1)+dmtarg(k))
        enddo

        do m=1,nfine
          zm = (m-0.5)/nfine
          if (zm.le.zmtarg(1)) then 
            dmfine(m) = dmtarg(1) 
          else if (zm.ge.zmtarg(nlay)) then 
            dmfine(m) = dmtarg(nlay) 
          else
            do k=2,nlay
              if (zm.ge.zmtarg(k-1) .and. zm.le.zmtarg(k)) then
                zw = (zm-zmtarg(k-1))/(zmtarg(k)-zmtarg(k-1))
                dmfine(m) = (1.-zw)*dmtarg(k-1) + zw*dmtarg(k)
              endif
            enddo
          endif
        enddo

#if defined (VADAPT)
c       set fixed layer mass thicknesses (kg m-2) for vadapt: 
        do k=1,nlay
          dmvadapt(k) = dmtarg(k)*totmass
        enddo
#endif

        do k=1,nlay
          snowm(k) = dmtarg(k)*totmass     ! kg m-2
          sliqm(k) = 0.                    ! kg m-2
          tfirn(k) = zt                    ! K
          rhosno(k) = 500.               ! kg m-3
          rgrain(k) = .001               ! m
#if defined (ANALTEST)
          tfirn(k) = -40. + tmelt
          rgrain(k) = .0003
#endif
        enddo

#if defined (IFPROG6)
        firstsurfmodel6(iloc,jloc) = .false.
c,,,,,,,,,
      else
c,,,,,,,,,
c       Recover prognostic model variables from 2-D arrays:

        do k=1,nlay
          snowm(k)  = asnowm(iloc,jloc,k)
          sliqm(k)  = asliqm(iloc,jloc,k)
          tfirn(k)  = atfirn(iloc,jloc,k)
          rhosno(k) = arhosno(iloc,jloc,k)
          rgrain(k) = argrain(iloc,jloc,k)
        enddo
c,,,,,,,,,,
      endif
c,,,,,,,,,,
#endif

      call zerobot (snowm, sliqm, tfirn, rhosno, rgrain,
     *              nlay, nlaymax)

      runofall = 0.
      deltotall = 0.

c        Overall annual loop

c@@@@@@@@@@@@@@@@@@@@@@@@@
      do iyear=1,nyearfirn
c@@@@@@@@@@@@@@@@@@@@@@@@@

c          Zero annual accumulators

        araina = 0.   ! annual rainfall (kg/m2)
        asnowa = 0.   ! annual snowfall (kg/m2)
        avapa  = 0.   ! annual evaporation (kg/m2)
        taira  = 0.   ! annual air temp (K)
        ticea  = 0.   ! top boundary condition for ice model (K)

        runof = 0.
        rmel = 0.
        refr = 0.
        deltot = 0.

        nlaychange_annp = 0.
        nlaychange_annm = 0.

c          Loop over timesteps through one year

c++++++++++++++++++
        do m=1,nspy
c++++++++++++++++++

          do k=1,nlay
            runofzer(k) = 0.
          enddo

          do k=1,nlay
            firnmass(k) = snowm(k) + sliqm(k)
          enddo
          firnmid(1) = 0.5*firnmass(1)
          do k=2,nlay
            firnmid(k) = firnmid(k-1) + 0.5*(firnmass(k-1)+firnmass(k))
          enddo

#if defined (CROSSSEC)
c         for streamlines in cross-sections: save z at start of timestep
          do k=1,nlay
            thick(k) = snowm(k)/max(rhosno(k),.001)
          enddo
          zvelprev(1) = 0.5*thick(1)
          do k=2,nlay
            zvelprev(k) = zvelprev(k-1) + 0.5*(thick(k-1) + thick(k))
          enddo
#endif

c-------------------------------------------------------------
c            Give snowfall (sna) to layer 1. Rainfall (raa) is 
c            given later, in "Adjustments".
c-------------------------------------------------------------

c         diurnal+synoptic, degree_day (erf):
          call degree_day (dtimed, ta(m), pdd, qdd, tmelt, sdt)
          raa = (1.-qdd) * pa(m)*1000./(365.*86400.)   ! m/yr to kg/m2/s
          sna =     qdd  * pa(m)*1000./(365.*86400.)   ! m/yr to kg/m2/s
#if defined (MEYERHEWITT)  
          pa(m) = 1.7 
          raa = 0.
          sna = pa(m)*1000./(365.*86400.)   ! m/yr to kg/m2/s
#endif

c         effect on layer 1 temperature, conserving heat
c         (effect on snow density, grain size are in eqns below):
          cheata = (cheaticea + cheaticeb*(tfirn(1)-tmelt)) 
     *             * snowm(1)
     *             + cheatliq * sliqm(1)
          zts = min (ta(m), tmelt)
          ztr = max (ta(m), tmelt)
          cheatb = (cheaticea + cheaticeb*(zts-tmelt)) * sna*dtimed
          cheatc = cheatliq * raa*dtimed
          tfirn(1) = (cheata*tfirn(1) + cheatb*zts + cheatc*ztr) 
     *             / (cheata+cheatb+cheatc)
#if defined (ANALTEST)
          tfirn(1) = -40. + tmelt
#endif
          snowm(1) = snowm(1) + sna*dtimed

c---------------------------------------------------------------------
c            Remove evaporation (eva) from layer 1. From snow+liquid 
c            (in proportion to their rel. masses)
c---------------------------------------------------------------------


          eva = ea(m)*1000./(365.*86400.)              ! m/yr to kg/m2/s

          e1sl = min (eva*dtimed, snowm(1)+sliqm(1))
          e1s  = e1sl * snowm(1)/max(snowm(1)+sliqm(1),.001) ! from snow
          e1l  = e1sl * sliqm(1)/max(snowm(1)+sliqm(1),.001) ! from liq
          snowm(1) = snowm(1) - e1s
          sliqm(1) = sliqm(1) - e1l

          call zeromass (snowm, sliqm, poros, rhosno, thick, fliq, 
     *                   runofzer, nlay)

c----------------------------------------------
c              Step snow density and grain size 
c----------------------------------------------
c....................
          do k=1,nlay
c....................

c              Step grain size squared (time implicit). For layer 1, 
c              "mixing" term due to snowfall follows assuming grain size
c              is a conservative tracer (see notes 2/17/17). 
c              Local term from Arthern et al. 2010.
c              d(rsq)/dt = cga*exp(-cgb/(crr*t))
c                        + (rgrafr2-rsq)*sna/snowm (if layer 1)  

            rsq = rgrain(k)**2
            rsqold = rsq
            zae = dtimed * cga*exp(-cgb/(crr*tfirn(k)))
c           modif to Arthern, tuning to Firnmice exp 1 rho(z), Fig. 3:
            zae = zae*8.
            zas = 0.
            if (k.eq.1) zas = sna*dtimed / max(snowm(k),.001)
            rsq = (rsq + zae + rgrafr2*zas) / (1.+ zas)

            rsq = max (rgrafr2, min (rgraic2, rsq))
            rgrain(k) = sqrt(rsq)
#if defined (ANALTEST)
            rgrain(k) = .0003
            rsq = rgrain(k)**2
#endif

c              Step layer 1 snow density (time implicit). For layer 1,
c              "mixing" term due to snowfall follows by assuming 
c              individual package volumes unchanged (see notes 2/17/17).
c              Local term from Arthern et al. 2010.
c              d(rhosno)/dt = 
c              (rhoice  -rhosno)*cra*exp(-crb/(crr*t))*csig/rsq
c              + (rhosnofr-rhosno)*(rhosno/rhosnofr)*sna/snowm (if ly 1)

#if defined (ANALTEST)
            cra = 7.e-9
#else
c           if (rhosno(k).le.550.) then
c             cra = cra1
c           else
c             cra = cra2
c           endif
c           or continuous transition:
c           zw = max (0., min (1., (rhosno(k)-500.)/100. ))
            zw = max(0.,min(1., (rhosno(k)-rhosnofr)/(rhoice-rhosnofr)))
            cra = (1.-zw)*cra1 + zw*cra2 
c           modif to Arthern, tuning to Firnmice exp 1 rho(z), Fig. 3:
            cra = cra*6.
#endif

            csig = firnmid(k)*grav

            zae = dtimed 
     *            * cra*exp(-crb/(crr*tfirn(k)))*csig/(0.5*(rsqold+rsq))
            zas = 0.
            if (k.eq.1) zas = 
     *            (rhosno(k)/rhosnofr) * (sna*dtimed/max(snowm(k),.001))

c           rhosnomax limit is to still accomodate sliqm with fliq=1:
            rhosnomax = 1. / ( (sliqm(k)/max(snowm(k),.001))/rhoh2o 
     *                          + 1./rhoice )
            rhosno(k) = (rhosno(k) + min(rhoice,rhosnomax)*zae 
     *                     + rhosnofr*zas)
     *                  / (1. + zae + zas)

            rhosno(k) = max (rhosnomin, min (rhoice, rhosno(k)))
c..............
          enddo
c..............

c            Calculate specific heats

          do k=1,nlay
            cheat(k) =  max ( (cheaticea+cheaticeb*(tfirn(k)-tmelt))
     *                        *snowm(k) + cheatliq*sliqm(k),
     *                        cheaticea*rhoice*dzmin )        !J m-2 K-1
          enddo

c            Calculate conduction cofficients
c            nb: no dependence of snow conductivity on pore liquid

          do k=1,nlay
            conds(k) = (condsnoa+condsnob*(rhosno(k)**2)) 
     *               / (0.5*max(snowm(k)/rhosno(k), dzmin))
          enddo

          do k=1,nlay-1
            cond(k) = 1. / (1./conds(k) + 1./conds(k+1))
          enddo
          cond(nlay) = 1. / ( 1./conds(nlay) 
                              + zice/(condicea*exp(-condiceb*tice)) )

          do k=1,nlay
            cond(k) = cond(k) / (86400.*365.)            ! J m-2 s-1 K-1
          enddo

c            Solar radiation, albedo calculation, for layer 1 

          poros(1) = max (0., 1. - rhosno(1)/rhoice)
          sliqmax = rhoh2o*fret(poros(1))*poros(1)*snowm(1)/rhosno(1)
          zl1 = max (0.,min (1., sliqm(1)/max(.001,sliqmax)))
          zr1 = max (0.,min (1.,(rhosno(1)-rhosnofr)/(rhoice-rhosnofr)))

c            Weight for embedded liquid amount (albsnodry to albsnowet)
c            and snow density (denser -> albsnoden), to get snow albedo.

          zalb = (1.-zr1) * ((1.-zl1)*albsnodry + zl1*albsnowet)
     *          +     zr1  * albsnoden 

#if defined (FCLOUDMIN)
          zcloud =fcloud*exp(-hs/1000.) + FCLOUDMIN*(1.-exp(-hs/1000.))
#else
          zcloud = fcloud
#endif

c         no multiple reflections:
c         zabspla = (1.-zalb)*(1.-zcloud + zcloud*atrans)
c         with multiple reflections, assuming atrans and cloud presence
c         are 100% correlated in time at each pt. (see notes 3/24/17):
          zabspla =   ((1.-zalb)/(1.-zalb*(1.-atrans)))
     *              * (1.-(1.-atrans)*(zcloud+zalb*(1.-zcloud)))

          zsol = zabspla * qa(m)

#if defined (MEYERHEWITT)  
          zq0 = 200.
          if (MEYERHEWITT.eq.1) then   
            zqm = -1.*zq0
          else
            zqm = MEYERHEWITT*zq0
          endif
          zsol = zqm - zq0*cos(2.*3.1415927*(m-.5)/nspy) - aconst
#endif

c----------------------------------------------------
c            Step temperatures (tfirn), time implicit
c----------------------------------------------------
c            cheat*d(T)/dt = acoeff*(ta-T) + zsol + aconst (if layer 1)
c                          + d(cond*dT/tz)/dz

c           2 passes (itriloop) of the tridiag solver. For the 2nd pass,
c           l.h.s. layer temperatures are set to tmelt if they exceeded 
c           tmelt in the 1st pass. 

c&&&&&&&&&&&&&&&&&&&&&&&&
          do itriloop=1,2
c&&&&&&&&&&&&&&&&&&&&&&&&
            if (itriloop.eq.1) then
              ifdo = 1
              itri = 1
              jtri = 0
            else if (itriloop.eq.2) then
              if (tfirn2(1).gt.tmelt) then
                ifdo = 1 
                itri = 0
                jtri = 1
              else
                ifdo = 0
              endif
            endif

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            if (ifdo.eq.1) then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

              do k=1,nlay
                zdt = dtimed/cheat(k)

                if (k.eq.1) then 
                  atri(k) = 1. + zdt*cond(k)*itri + zdt*acoeff*itri
                  btri(k) = -zdt*cond(k)
                  ctri(1) = 0.
                  rhs(k) = tfirn(k) 
     *                   - zdt * cond(k)*tmelt*jtri
     *                   + zdt * (acoeff*(ta(m)-tmelt*jtri)+zsol+aconst)

                else if (k.eq.2) then
                  atri(k) = 1. + zdt*cond(k) + zdt*cond(k-1)
                  btri(k) = -zdt*cond(k)
                  ctri(k) = -zdt*cond(k-1)*itri
                  rhs(k) = tfirn(k) 
     *                   + zdt*cond(k-1)*tmelt*jtri

                else if (k.gt.2) then
                  atri(k) = 1. + zdt*cond(k) + zdt*cond(k-1)
                  if (k.lt.nlay) then
                    btri(k) = -zdt*cond(k)
                  else
                    btri(nlay) = 0.
                  endif
                  ctri(k) = -zdt*cond(k-1)
                  rhs(k) = tfirn(k) 
                  if (k.eq.nlay) rhs(k) = rhs(k) 
     *                                  + zdt*cond(k)*tice
                endif

              enddo

              call tridia_i (atri, btri, ctri, tfirn2, rhs, nlay, 0)

              if (itriloop.eq.1) zheat = 0.
              if (itriloop.eq.2) zheat = (   acoeff*(ta(m)-tmelt) 
     *                                     + zsol + aconst
     *                                     - cond(1)*(tmelt-tfirn2(2))
     *                                   ) * zdt 

c>>>>>>>>>>>>>>>>
            endif
c>>>>>>>>>>>>>>>>
c&&&&&&&&&&&&&&
          enddo
c&&&&&&&&&&&&&&

          do k=1,nlay
            tfirn(k) = tfirn2(k)
#if defined (ANALTEST)
            tfirn(k) = -40. + tmelt
#endif
          enddo

c-----------------------
c            Adjustments.
c
c            One pass downwards, 1 to nlay, propagating zheat, perccur.
c            For solid-ice (poros=0) layers, use zheat to warm tfirn 
c            to(wards) tmelt, if reached melt some ice (-> runoff).

c            For permeable layers (poros>0), add perccur from above. 
c            Use zheat to warm/melt ice. Freeze liquid if tfirn < tmelt.
c            Pass excess liquid (> fret*thick*poros) downwards.

c            At bottom of each contiguous zone of permeable layers
c            (i.e., if solid-ice layer reached, or k+1=nlay), do
c            reverse upwards pass through this contiguous zone, 
c            allocating liquid (perccur -> percrev) to saturate
c            each layer in turn (fliq=1). Any excess liquid at
c            top of this permeable zone -> runoff. Track current 
c            permeable zone  by [kcura:kcurb].

c-----------------------

          perccur(0) = raa*dtimed
          percrev(0) = 0.
          do k=1,nlay
            perccur(k) = 0.
            percrev(k) = 0.
            runofcur(k) = 0.
          enddo
          rmelcur = 0.
          refrcur = 0.

          kcura = 0
          kcurb = 0

          do k=1,nlay

c-----------------------------------
            if (poros(k).eq.0.) then
c-----------------------------------

              kcura = 0
              kcurb = 0

c             solid-ice layer: use zheat to warm tfirn to(wards) 
c             tmelt, if reached melt some ice (-> runoff).   

              cheat(k) =  max ( (cheaticea+cheaticeb*(tfirn(k)-tmelt))
     *                          *snowm(k),
     *                          cheaticea*rhoice*dzmin )

              tfirn(k) = tfirn(k) + zheat/cheat(k)

              if (tfirn(k).gt.tmelt) then
                zheat = cheat(k)*(tfirn(k)-tmelt)
                zmelt = min (zheat/hfus, snowm(k))
                zheat = zheat - zmelt*hfus
                snowm(k) = snowm(k) - zmelt 
                runofcur(k) = runofcur(k) + zmelt + perccur(k-1)
                perccur(k) = 0.
                tfirn(k) = tmelt 
                rmelcur = rmelcur + zmelt 
              else
                zheat = 0.
                perccur(k) = 0.
              endif

c---------------
            else
c---------------

              if (kcura.eq.0) kcura = k

c             permeable layer: add perccur from above. Use zheat to 
c             warm layer/melt ice, or freeze liquid if tfirn<tmelt. 
c             Pass excess liquid (> fret*thick*poros) downwards.

              sliqm(k) = sliqm(k) + perccur(k-1)

              cheat(k) =  max ( (cheaticea+cheaticeb*(tfirn(k)-tmelt))
     *                          *snowm(k) + cheatliq*sliqm(k),
     *                          cheaticea*rhoice*dzmin )

              tfirn(k) = tfirn(k) + zheat/cheat(k)

c             If t > tmelt, melt some snow to liquid. 
c             Preserves layer thickness (snowm/rhosno), unless min snow
c             density reached (rhosnomin, with snowm=0, to protect 
c             sliqmax calcs below)
    
              if (tfirn(k).gt.tmelt) then
                zheat = cheat(k)*(tfirn(k)-tmelt)
                zmelt = min (zheat/hfus, snowm(k))
                zheat = zheat - zmelt*hfus
    
                rhosno(k) = rhosno(k) - zmelt / (snowm(k)/rhosno(k)) 
                rhosno(k) = max (rhosnomin, min (rhoice, rhosno(k)))
    
                sliqm(k) = sliqm(k) + zmelt 
                snowm(k) = snowm(k) - zmelt 
    
                tfirn(k) = tmelt 
                rmelcur = rmelcur + zmelt 
              else
                zheat = 0.
              endif   

c             If t < tmelt, freeze some liquid (sliqm) to snow.
c             Preserves layer thickness (snowm/rhosno).

              if (tfirn(k).lt.tmelt .and. sliqm(k).gt.0.) then
                zfreez = min (cheat(k)*(tmelt-tfirn(k))/hfus, sliqm(k))

                rhosno(k) = rhosno(k) + zfreez / (snowm(k)/rhosno(k)) 
                rhosno(k) = max (rhosnomin, min (rhoice, rhosno(k)))

                sliqm(k)  = sliqm(k) - zfreez
                snowm(k)  = snowm(k) + zfreez

                tfirn(k) = tfirn(k) + zfreez*hfus/cheat(k)
                refrcur = refrcur + zfreez
              endif

c             Adjust max liquid retention, give excess to percolation
c             to layer below

              poros(k) = max (0., 1. - rhosno(k)/rhoice)
              sliqmax= rhoh2o*fret(poros(k))*poros(k)*snowm(k)/rhosno(k)
              if (sliqm(k).gt.sliqmax) then
                perccur(k) = sliqm(k)-sliqmax
                sliqm(k) = sliqmax
              else
                perccur(k) = 0.
              endif

              kcurb = k

c......................................................
              if (poros(k+1).eq.0. .or. k.eq.nlay) then
c......................................................

c               At bottom of a contiguous permeable zone (kcura:kcurb).
c               Do reverse upwards pass through this contiguous zone, 
c               allocating liquid (perccur -> percrev) to saturate
c               each layer in turn (fliq=1). For each layer,
c               allow lateral runoff (which reduces percrev), depending
c               on porosity and saturation (sliqm/sliqmax), 
c               set by statement function runlat.
c               Any excess liquid at top of this permeable zone->runoff.

                percrev(kcurb) = perccur(kcurb)
                perccur(kcurb) = 0.

                do kr=kcurb, kcura, -1
                  poros(kr) = max (0., 1. - rhosno(kr)/rhoice)
                  sliqmax = rhoh2o*poros(kr)*snowm(kr)/rhosno(kr)
                  sliqm(kr) = sliqm(kr)+percrev(kr)
                  zr = min ( sliqm(kr), 
     *                       runlat (sliqm(kr), sliqmax, poros(kr))
     *                       * dtimed )
                  runofcur(kr) = runofcur(kr) + zr
                  sliqm(kr) = sliqm(kr) - zr
                  if (sliqm(kr).gt.sliqmax) then
                    percrev(kr-1) = sliqm(kr) - sliqmax
                    sliqm(kr) = sliqmax
                  else
                    percrev(kr-1) = 0.
                  endif
                  if (kr.eq.kcura) 
     *              runofcur(kr) = runofcur(kr) + percrev(kr-1)
                enddo

c..................
              endif
c..................

c----------------
            endif
c----------------

          enddo

c            snowm could have become zero in Adjustments. Protect 
c            this for layer shifting or vadapt below.

          call zeromass (snowm, sliqm, poros, rhosno, thick, fliq, 
     *                   runofzer, nlay)

#if defined (CROSSSEC)
c         for streamlines in cross-sections: calc z at end of timestep
c         (but before layer merge/splitting below). Up to here,
c         layer indexing is unchanged from staret of timestep.
c         Save z-change in timestep (wvel) for crosssec output below.
          do k=1,nlay
            thick(k) = snowm(k)/max(rhosno(k),.001)
          enddo
          zvel(1) = 0.5*thick(1)
          do k=2,nlay
            zvel(k) = zvel(k-1) + 0.5*(thick(k-1) + thick(k))
          enddo
          do k=1,nlay
            wvel(k) = (zvel(k)-zvelprev(k))/dtimed
          enddo
#endif

c///////////////////
#if defined (VADAPT)
c///////////////////
   
c            Re-adapt layer column mass thicknesses ("semi-Lagrangian")
   
          ztot = 0.
          do k=1,nlay
            ztot = ztot + snowm(k) + sliqm(k)
          enddo

          deltot = deltot + (ztot-totmass)

          call vadapt3 (dmvadapt,
     *                  snowm, 1000.*rhoice,
     *                  sliqm,  0.,
     *                  tfirn,  tice,
     *                  rhosno, rhoice, 
     *                  rgrain, rgraic, 
     *                  nlay)

c////
#else
c////

c            Layer shifting (lagrangian, temporary nlaycur)

c         for layer-shifting diagnostic dump info to iudiagsh:    
          diagshift = .false.
#if defined (DIAGSHIFT)
          if ( iudiagsh.ne.0 
c    *         .and. iyear.eq.nyearfirn
c    *         .and. mod (m, max(1,nint(5.*86400./dtimed))).eq.0
     *       ) diagshift = .true. 
#endif

          nlaycur = nlay

          if (diagshift) write (iudiagsh,'(/a,15f10.3, i10)')
     *      'above c(1) snowm(1-15):',(snowm(ix),ix=1,15), nlaycur

c--------------------------------------------------
c(1)         Remove zero-mass (snowm) layers, 
c            reducing number of layers
c--------------------------------------------------

          kk = 0
          do k=1,nlaycur
            if (snowm(k).eq.0.) then
c             remove layer: no snow
            else
c             keep layer: non-zero snow
              kk = kk + 1
              if (kk.ne.k) then
                snowm(kk) = snowm(k)
                sliqm(kk) = sliqm(k)
                tfirn(kk) = tfirn(k)
                rhosno(kk) = rhosno(k)
                rgrain(kk) = rgrain(k)
              endif
            endif
          enddo
          nlaycur = kk

          call zerobot (snowm, sliqm, tfirn, rhosno, rgrain,
     *                  nlaycur, nlaymax)

          if (diagshift) write (iudiagsh,'(a,15f10.3, i10)')
     *      'below c(1) snowm(1-15):',(snowm(ix),ix=1,15), nlaycur
              
c----------------------------------------------------------------------
c(2)         Adjust layer(s) to restore total column mass (-> totmass),
c            can increase or decrease number of layers 
c----------------------------------------------------------------------

          ztot = 0.
          do k=1,nlaycur
            ztot = ztot + snowm(k) + sliqm(k) 
          enddo

          deltot = deltot + (ztot-totmass)

          if (ztot.lt.totmass) then

c           add one layer on bottom with ice-model properties:
            nlaycur = nlaycur + 1
            if (nlaycur.gt.nlaymax) then
              write (6,'(a,i6)') '*** Error(3): nlaycur(>max)=',nlaycur
              stop
            endif
            snowm(nlaycur) = totmass - ztot
            sliqm(nlaycur) = 0.
            tfirn(nlaycur) = tice
            rhosno(nlaycur) = rhoice
            rgrain(nlaycur) = rgraic

          else if (ztot.gt.totmass) then

c           remove partial or ful1+partial layer(s) from bottom, iterate
  300       zb = snowm(nlaycur) + sliqm(nlaycur)
            zf = min (1., (ztot-totmass)/zb)
            if (zf.lt.1.) then
              snowm(nlaycur) = snowm(nlaycur)*(1.-zf)
              sliqm(nlaycur) = sliqm(nlaycur)*(1.-zf)
            else
              nlaycur = nlaycur - 1
              ztot = ztot - zb
              if (ztot.gt.totmass) goto 300 
            endif   

          endif   

          call zerobot (snowm, sliqm, tfirn, rhosno, rgrain,
     *                  nlaycur, nlaymax)

          if (diagshift) write (iudiagsh,'(a,15f10.3, i10)')
     *      'below c(2) snowm(1-15):',(snowm(ix),ix=1,15), nlaycur

c-----------------------------------------------------------------------
c(3)         Ged rid of very thin layers, by merging with layer below
c            (iteratively, until none left), decreasing number of layers
c-----------------------------------------------------------------------

  400     continue

          do k=1,nlaycur
            zmass = snowm(k) + sliqm(k)
            if (zmass.lt.dmmin) then
              ka = min (k, nlaycur-1) 
              call layermerge (snowm,sliqm,tfirn,rhosno,rgrain,
     *                         ka, nlaycur)
              go to 400
            endif
          enddo

          call zerobot (snowm, sliqm, tfirn, rhosno, rgrain,
     *                  nlaycur, nlaymax)

          if (diagshift) write (iudiagsh,'(a,15f10.3, i10)')
     *      'below c(3) snowm(1-15):',(snowm(ix),ix=1,15), nlaycur

c---------------------------------------------------------------
c(4)         Maintain number of layers (restore nlaycur -> nlay),
c            also restore thicknesses back towards target thicknesses 
c            vs. depth (dmtarg(1:nlay). Do this iteratively (500 loop),
c            by computing scores score_[split,merge] favoring
c            splitting or merging respectively. All scores are 
c            calculated using non-dim depths, normalized by totmass.
c---------------------------------------------------------------

          nlaychange      = nlaycur - nlay  
          nlaychange_annp = nlaychange_annp + max (nlaychange,0)
          nlaychange_annm = nlaychange_annm + min (nlaychange,0)

          niter5 = 0

c            Calc preliminary non-dim depth vars (dmlay, zmlay)

  500     continue

          niter5 = niter5 + 1

c         (ztot should equal totmass, after # 3 above)
          ztot = 0. 
          do k=1,nlaycur
            ztot = ztot + snowm(k) + sliqm(k)
          enddo

          do k=1,nlaycur
            dmlay(k) = (snowm(k) + sliqm(k))/ztot
          enddo

          zmlay(1) = 0.5*dmlay(1)
          do k=2,nlaycur
            zmlay(k) = zmlay(k-1) + 0.5*(dmlay(k-1)+dmlay(k))
          enddo

c            Calculate score for splitting (score_split) and
c            its location (k_split). score_split is the max ratio of 
c            layer thickness (dmlay) to target thickness for that 
c            depth (dmfine(ifine)).

          zsmax = -1.e20
          do k=1,nlaycur
            ifine = max (1, min (nfine, nint (zmlay(k)/dfine + 0.5) ))
            zs = dmlay(k) / dmfine (ifine)
            if (zs.gt.zsmax) then
              zsmax = zs
              ksmax = k
            endif
          enddo
          score_split = zsmax  
          k_split = ksmax  

c            Calculate score for merging (score_merge) and layer-pair 
c            location (k_merge, k_merge+1). score_merge is 
c            based on ratios of the 2 layer thicknesses to target 
c            thickness for their depth (or inverse to make > 1, 
c            zra and zrb, large favored), and the ratio of the 
c            hypothetical merged layer to its target thickness for the 
c            hypoth. depth (or inverse, zrmerge, small favored).

          zsmax = -1.e20
          do k=1,nlaycur-1
            ka = k
            kb = k+1

            ifa = max (1, min (nfine, nint (zmlay(ka)/dfine + 0.5) ))
            zra = dmlay(ka)/dmfine(ifa)
            zra = max (zra, 1./zra)

            ifb = max (1, min (nfine, nint (zmlay(kb)/dfine + 0.5) ))
            zrb = dmlay(kb)/dmfine(ifb)
            zrb = max (zrb, 1./zrb)

            dmerge = dmlay(ka) + dmlay(kb)  
c           zmerge = z(ka)-.5*d(ka)+.5*(d(ka)+d(kb)) = z(ka)+0.5*d(kb)
            zmerge = zmlay(ka) + 0.5*dmlay(kb)
            ifmerge = max (1, min (nfine, nint (zmerge/dfine + 0.5) ))
            zrmerge = dmerge/dmfine(ifmerge)
            zrmerge = max (zrmerge, 1./zrmerge)

            zs = max (zra, zrb) / zrmerge
            if (zs.gt.zsmax) then
              zsmax = zs
              ksmax = ka
            endif
          enddo
          score_merge = zsmax  
          k_merge = ksmax  

c            Now modify score2_[splot, merge] by combining with the 
c            current difference between nlaycur and (desired) nlay: 
c            if nlay > nlaycur, split is favored more; 
c            if nlay < nlaycur, merge is. The nlaycur influence is 
c            weighted by zalph. 

          score2_split = max (score_split - score_ok, 0.) 
     *                   + zalph * (nlay-nlaycur)
          score2_merge = max (score_merge - score_ok, 0.) 
     *                   + zalph * (nlaycur-nlay)

          if (diagshift) then
            write (iudiagsh,'(a,15f10.3, i10, 4f10.3, 2i6)')
     *        '   in c(4) snowm(1-15):',(snowm(ix),ix=1,15), nlaycur,
     *        score_split, score_merge, score2_split, score2_merge, 
     *        niter5, niter5max
          endif

c            Decide what to to: skip out, split a layer, or merge
c            a layer pair. To skip out (done with iteration), require 
c            nlaycur = nlay, and score_[split,merge] both <= score_ok.
c            (Or if # of iters = nitermax). 
c            If not skip out, break tie (score2_split = score2_merge)
c            by nlay < or > nlaycur.

c..................................................................
          if ( nlaycur.eq.nlay .and.  
     *         ( (score_split.le.score_ok .and. score_merge.le.score_ok)
     *           .or. 
     *           (niter5.ge.niter5max) )
     *       ) then
c..................................................................

c            done:

c.............
          else
c.............

c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
            if (  score2_split.gt.score2_merge .or.
     *           (score2_split.eq.score2_merge .and. nlay.gt.nlaycur)
     *         ) then
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

c              split layer k_split:

              nlaycur = nlaycur + 1
              if (nlaycur.gt.nlaymax) then
                write(6,'(a,i6)') '*** Error(5): nlaycur(>max)=',nlaycur
                stop    
              endif   
              if (k_split+2.le.nlaycur) then
                do k = nlaycur, k_split+2, -1
                  snowm(k) = snowm(k-1)
                  sliqm(k) = sliqm(k-1)
                  tfirn(k) = tfirn(k-1)
                  rhosno(k) = rhosno(k-1)
                  rgrain(k) = rgrain(k-1)
                enddo   
              endif   
              do k = k_split+1, k_split, -1
                snowm(k) = 0.5*snowm(k_split)
                sliqm(k) = 0.5*sliqm(k_split)
                tfirn(k) = tfirn(k_split)
                rhosno(k) = rhosno(k_split)
                rgrain(k) = rgrain(k_split)
              enddo   

c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            else if (  score2_split.lt.score2_merge .or.
     *                (score2_split.eq.score2_merge.and.nlay.le.nlaycur)
     *              ) then
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

c              merge layers k_merge, k_merge+1:

               call layermerge (snowm,sliqm,tfirn,rhosno,rgrain,
     *                          k_merge, nlaycur)

c . . . . . . . .
            endif
c . . . . . . . .

            call zerobot (snowm, sliqm, tfirn, rhosno, rgrain,
     *                    nlaycur, nlaymax)

            goto 500

c..............
          endif
c..............

          if (diagshift) then
            write (iudiagsh,'(a,15f10.3, i10)')
     *        'below c(4) snowm(1-15):',(snowm(ix),ix=1,15), nlaycur
            write (iudiagsh,'(a,10f10.3)')
     *        '                       ',(snowm(ix)/ztot,ix=1,10)
            write (iudiagsh,'(a,10f10.3)')
     *        '                       ',(dmtarg(ix),ix=1,10)
          endif

c            end of layer shifting  ! now nlaycur=nlay, col mass~totmass

c/////
#endif
c/////

c            Reset [thick,fliq,poros] from [snowm,sliqm,rhosno] after 
c            all changes above, enforcing rhosnomax condition which
c            could have been exceeded by layer-shifting or vadapt
c            (ssr2tfp)
   
          call zeromass (snowm, sliqm, poros, rhosno, thick, fliq,
     *                   runofzer, nlay)

          call ssr2tfp (thick, fliq,
     *                  snowm, sliqm,
     *                  poros,
     *                  rhosno, rhoh2o, rhoice, nlay)

c            Increment annual snowfall, rainfall (kg/m2), air temp,
c            and surface b.c. temp for ice model top (ticea)  

          araina = araina + dtimed*raa
          asnowa = asnowa + dtimed*sna
          avapa  = avapa  + dtimed*eva
          taira  = taira  + ta(m)/nspy
          ticea  = ticea  + tfirn(nlay)/nspy

          do k=1,nlay
            runof = runof + runofcur(k) + runofzer(k)
          enddo
          rmel = rmel + rmelcur
          refr = refr + refrcur

#if defined (ANALTEST)
          do k=1,nlay
            tfirn(k) = -40. + tmelt
          enddo
#endif

c            Reset thick, depth (firnmass/mid, not used), for
c            diagnostic outputs below

          do k=1,nlay
            thick(k) = snowm(k)/rhosno(k)
          enddo

          depth(1) = 0.5*thick(1)
          do k=2,nlay
            depth(k) = depth(k-1) + 0.5*(thick(k-1) + thick(k))
          enddo

          do k=1,nlay
            firnmass(k) = snowm(k) + sliqm(k)
          enddo
          firnmid(1) = 0.5*firnmass(1)
          do k=2,nlay
            firnmid(k) = firnmid(k-1) + 0.5*(firnmass(k-1)+firnmass(k))
          enddo

c            Diagnostic output to iudiag1, seasonal time step

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#if defined (DIAGALLYR)
          if (iudiag1.ne.0) then
#else
          if (iudiag1.ne.0 .and. iyear.eq.nyearfirn) then
#endif
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            if (m.eq.1) then
              tav = 0.
              pav = 0.
              do mm=1,nspy
                tav = tav + ta(mm)/nspy
                pav = pav + pa(mm)/nspy
              enddo

              write (iudiag1,'(/2i4, 5(3x,a,f7.3))')
     *          iloc,jloc, 
     *          'lon= ',alond(iloc,jloc), 'lat= ',alatd(iloc,jloc),
     *          'tann= ',tav-tmelt, 'pann= ',pav, 'tamax= ',ta(1)-tmelt

              write (iudiag1,'(4a, 9(10(a8,i3.3)))')
     *          '    year    m     day',
     *          '        ta       raa       sna       eva',
     *          '     fluxa     fsola       alb    albeff',
     *          '    nlaych',
     *          ('   thick',k, k=1,nlay,nskip),
     *          ('   snowm',k, k=1,nlay,nskip),
     *          ('   sliqm',k, k=1,nlay,nskip),
     *          ('    fliq',k, k=1,nlay,nskip),
     *          ('  rhosno',k, k=1,nlay,nskip),
     *          ('  rgrain',k, k=1,nlay,nskip),
     *          ('   tfirn',k, k=1,nlay,nskip),
     *          (' perccur',k, k=1,nlay,nskip),
     *          (' runofcr',k, k=1,nlay,nskip)
            endif

c           if (mod ( m, max(1,nint(5.*86400./dtimed)) ).eq.0) then
            if (.true.) then

              zfa = acoeff*(ta(m)-tfirn(1)) + zsol + aconst
c             zz = 86400.*365. ; zz2 = nspy 
              zz = dtimed      ; zz2 = 1.

              write ( iudiag1,
     *          '(i8,i5,f8.2, 
     *            4f10.3,
     *            4f10.3, 
     *            i10,
     *            9(10f11.3))' )
     *          nint(iyear+timeice), m, 
     *          (mod(m-mmax+nspy,nspy)+1)*dtimed/86400.,
     *          ta(m)-tmelt, raa*zz, sna*zz, eva*zz, 
     *          zfa, zsol, zalb, 1.-zabspla,
     *          nlaychange,
     *          (thick(k), k=1,nlay,nskip),
     *          (snowm(k), k=1,nlay,nskip),
     *          (sliqm(k), k=1,nlay,nskip),
     *          (fliq(k), k=1,nlay,nskip),
     *          (rhosno(k), k=1,nlay,nskip),
     *          (rgrain(k)*1000., k=1,nlay,nskip),
     *          (tfirn(k)-tmelt, k=1,nlay,nskip),
     *          ((perccur(k)-percrev(k))*zz2, k=1,nlay,nskip),
     *          ((runofcur(k)+runofzer(k))*zz2, k=1,nlay,nskip)
            endif
c^^^^^^^^^^^^^^
          endif
c^^^^^^^^^^^^^^

c{{{{{{{{{{{{{{{{{{{{{
#if defined (CROSSSEC)
c{{{{{{{{{{{{{{{{{{{{{

c            For cross-section (z vs day) output - save fields
c            through days of the last year. If last timestep of year, 
c            call crosssec6 to transfer to z coords, and output 
c            to files for gp (load.crosssec).
c            nb: wvel has different vert coords (zvel, before final
c            adaptions) than other fields.  

#if defined (CROSSSEC2YR)
          if (iudiag1.ne.0 .and. iyear.ge.nyearfirn-1) then
            mc = m + nspy*(iyear-(nyearfirn-1))
#else
          if (iudiag1.ne.0 .and. iyear.eq.nyearfirn) then
            mc = m
#endif
            crossd(mc) = (mc-.5)*dtimed/86400.

            do k=1,nlay
              crossz(k,mc) = depth(k)
              crosst(k,mc) = tfirn(k)-tmelt
              crosspor(k,mc) = poros(k)
c             saturation:
c             crosssat(k,mc) = fliq(k)
c             liquid mass fraction:
              crosssat(k,mc) = sliqm(k) / max (sliqm(k)+snowm(k), .001)
              crosszvel(k,mc) = zvel(k) 
              crosswvel(k,mc) = wvel(k) 
            enddo

            if (m.eq.nspy) then
              call crosssec6 (crossd, crossz, crosst,    nlay,nspycross,
     *                       'temperature', 1,timeice,nyearfirn,iudiag6)
              call crosssec6 (crossd, crossz, crosspor,  nlay,nspycross,
     *                       'porosity',    0,timeice,nyearfirn,iudiag6)
              call crosssec6 (crossd, crossz, crosssat,  nlay,nspycross,
     *                        'liqfrac',    0,timeice,nyearfirn,iudiag6)
              call crosssec6 (crossd,crosszvel,crosswvel,nlay,nspycross,
     *                        'wvel',       0,timeice,nyearfirn,iudiag6)
            endif
          endif
c{{{{{
#endif
c{{{{{

c++++++++++++
        enddo ! m, seasonal layer intra-annual loop
c++++++++++++

c          Accumulate annual runoff, dmass over yrs (runofall,deltotall)

        runofall  = runofall  + runof/nyearfirn            ! kg m-2 yr-1
        deltotall = deltotall + deltot/nyearfirn           ! kg m-2 yr-1

c          Diagnostic output to iudiag*, end of year

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        if (iudiag2.ne.0 .or. iudiag3.ne.0 .or. iudiag4.ne.0 .or. 
     *      iudiag5.ne.0) then
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

          tav = 0.
          pav = 0.
          do mm=1,nspy
            tav = tav + ta(mm)/nspy
            pav = pav + pa(mm)/nspy
          enddo

          zt = 0.
          zs = 0.
          zl = 0.
          zm = 0.
          do k=1,nlay
            zt = zt + thick(k) 
            zs = zs + snowm(k)
            zl = zl + sliqm(k)
            zm = zm + snowm(k) + sliqm(k)
          enddo

          zdip = 0.
          do k=1,nlay
            zdip = zdip + thick(k)*poros(k)*(1.-fliq(k))
          enddo

          zdipa = 0.
          do k=1,nlay
            zdipa = zdipa + thick(k)*poros(k)
          enddo

          zbco = 0.
          if (815.le.rhosno(1)) then
            zbco = 0.
          else
            do k=1,nlay-1
              if (815.gt.rhosno(k) .and. 815.le.rhosno(k+1)) then
                zw = (815.-rhosno(k)) / (rhosno(k+1)-rhosno(k))
                zbco = (1.-zw)*depth(k) + zw*depth(k+1)
                goto 800
              endif
            enddo
            zbco = zt
  800       continue
          endif

c         dliq should = annual change in all sliqm*, which should
c         be ~zero if model has reached "equilibrium" (lagrangian 
c         transfer with ice-model ice below does not involve liquid).
          dliq = araina + rmel - refr - runof

          if (iyear.eq.1 .and. iudiag2.ne.0) then
            write (iudiag2,'(38("=")/2i4, 2(3x,a,f7.3))')
     *        iloc,jloc, 
     *        'lon= ',alond(iloc,jloc), 'lat= ',alatd(iloc,jloc)

            write (iudiag2,'(4a, 7(10(a8,i3.3)))' )
     *        '    year',
     *        '    asnowa    araina     avapa     taira     ticea',
     *        '      rmel      refr     runof    deltot',
     *        '  nlaych+-',
     *        ('   thick',k, k=1,nlay,nskip),
     *        ('   snowm',k, k=1,nlay,nskip),
     *        ('   sliqm',k, k=1,nlay,nskip),
     *        ('    fliq',k, k=1,nlay,nskip),
     *        ('  rhosno',k, k=1,nlay,nskip),
     *        ('  rgrain',k, k=1,nlay,nskip),
     *        ('   tfirn',k, k=1,nlay,nskip)
          endif

          if (mod(iyear,max(nyearfirn/5,1)).eq.0) then
c         if (iyear.eq.nyearfirn) then
c         if (.true.) then

            if (iudiag2.ne.0) write (iudiag2, 
     *        '(i8, 
     *          5f10.3, 
     *          4f10.3, 
     *          2i5, 
     *          7(10(f11.3)))') 
     *        nint(iyear+timeice), 
     *        asnowa, araina, avapa, taira-tmelt, ticea-tmelt,
     *        rmel, refr, runof, deltot,
     *        nlaychange_annp, nlaychange_annm,
     *        (thick(k), k=1,nlay,nskip),
     *        (snowm(k), k=1,nlay,nskip),
     *        (sliqm(k), k=1,nlay,nskip),
     *        (fliq(k), k=1,nlay,nskip),
     *        (rhosno(k), k=1,nlay,nskip),
     *        (rgrain(k)*1000., k=1,nlay,nskip),
     *        (tfirn(k)-tmelt, k=1,nlay,nskip)
          endif

c.....................................
          if (iyear.eq.nyearfirn) then
c.....................................

            if (iudiag3.ne.0) write (iudiag3,'(a,i10/a,i10)') 
     *          'year   =',nint(iyear+timeice), 'nlay   = ',nlay

            do iloopout=1,2
              if (iloopout.eq.1) then 
                iu = iudiag2
                cform = '(/5(4x,a,f10.4)' //
     *                  ' /5(4x,a,f10.4), 2(4x,a,e10.3) ' //
     *                  ' /7(4x,a,f10.2)/)'
              else
                iu = iudiag3   ! for ensembles (doenssurf6, doanalsurf6)
                cform = '( 5(a,f10.4/),' // 
     *                  '  5(a,f10.4/), 2(a,e10.3/),' //
     *                  '  7(a,f10.2/) )'
              endif

              if (iu.ne.0) write (iu,cform)
     *        'tann   = ',taira-tmelt, 
     *        'tamax  = ',ta(1)-tmelt,
     *        'pann   = ',(araina+asnowa)*.001,    ! kg m-2 yr-1 to m/yr
     *        'snowa  = ',asnowa*.001,
     *        'raina  = ',araina*.001, 

     *        'rmel   = ',rmel*.001, 
     *        'refr   = ',refr*.001, 
     *        'runof  = ',runof*.001, 
     *        'deltot = ',deltot*.001,
     *        'budtot = ',(asnowa+araina-avapa-runof)*.001,
     *        'errbud = ',(asnowa+araina-avapa-runof - deltot)*.001,
     *        'dliq   = ',dliq*.001,

     *        'tthic  = ',zt,
     *        'tsnom  = ',zs,
     *        'tliqm  = ',zl,
     *        'tmassm = ',zm,
     *        'totmas = ',totmass,
     *        'dip    = ',zdip,
     *        'bco    = ',zbco

            enddo

            if (iudiag4.ne.0) then 
              write (iudiag4,'(/a,i8)') 'year=',nint(iyear+timeice) 
              write (iudiag4,'(5a)') 
     *          '   k',
     *          '     depth     thick      temp',
     *          '     snowm     poros    rhosno    rgrain',
     *          '     sliqm      fliq',
     *          '     runof'
              do k=1,nlay
                write (iudiag4,'(i4,10f10.3)') k, 
     *            depth(k), thick(k), tfirn(k)-tmelt, 
     *            snowm(k), poros(k), rhosno(k), rgrain(k)*1000., 
     *            sliqm(k), fliq(k),
     *            (runofcur(k)+runofzer(k))
     *              *((86400.*365.)/dtimed)*.001                  ! m/yr
              enddo
c             extra line for interface with ice model:
              write (iudiag4,'(i4,10f10.3)') nlay+1, 
     *          depth(nlay)+0.5*thick(nlay), 0., tice-tmelt, 
     *          0., 0., rhoice, rgraic*1000., 
     *          0., 0.,
     *          0.
            endif
c..............
          endif
c..............

          if (iudiag5.ne.0) then 
c           if (firstdiag5) write (iudiag5,'(5a)')
            if (iyear.eq.1) write (iudiag5,'(5a)')
     *   '    year', 
     *   '        tann        pann       snowa       raina',
     *   '       tthic       tsnom       tliqm      tmassm',
     *   '         dip        fliq         bco',
     *   '        rmel        refr       runof'
            firstdiag5 = .false.

            zfliq = 1. - zdip/max(zdipa,.001)
            write (iudiag5,'(i8,14f12.3)') nint(timeice+iyear),
     *        taira-tmelt, (araina+asnowa)*.001,asnowa*.001,araina*.001,
     *        zt, zs, zl, zm,
     *        zdip, zfliq, zbco, 
     *        rmel*.001, refr*.001, runof*.001     ! kg m-2 yr-1 to m/yr
          endif

c^^^^^^^^^^^^
        endif
c^^^^^^^^^^^^

c@@@@@@@@@@
      enddo  ! iyear annual loop
c@@@@@@@@@@

c        Transfer to budg* for ice model, converting flux units 
          
      budgsnow = asnowa*.001                       ! kg m-2 yr-1 to m/yr
      budgrain = araina*.001                       ! kg m-2 yr-1 to m/yr
      budgevap = avapa*.001                        ! kg m-2 yr-1 to m/yr
#if defined (IFPROG6)
c     if prognostic call, return time average:
      budgmelt = runofall*.001                   ! kg m-2 yr-1 to m/yr
#else
c     if equilibrated call, return last year:
      budgmelt = runof*.001                      ! kg m-2 yr-1 to m/yr
#endif

c     also returned (deg K):
      tsurf = taira
      tsurfi= ticea

c       Save prognostic model variables to 2-D arrays (local, saved):

#if defined (IFPROG6)
      do k=1,nlay
        asnowm(iloc,jloc,k) = snowm(k)
        asliqm(iloc,jloc,k) = sliqm(k)
        atfirn(iloc,jloc,k) = tfirn(k)
        arhosno(iloc,jloc,k) = rhosno(k)
        argrain(iloc,jloc,k) = rgrain(k)
      enddo
#endif

      return
      end

c-----------------------------------------------------------------------

c///////////////////
#if defined (VADAPT)
c///////////////////

      subroutine vadapt3 (dmvadapt,
     *                    snow, snow4,
     *                    sliq, sliq4,
     *                    t,     t4,
     *                    rho,   rho4,
     *                    rgr,   rgr4,
     *                    nlay)
#if defined (STANDALONE6)
#include <comparams.h>
#else
#include <comicephys.h>
#include <comicegrid.h>
#endif

c        Like vadapt, vadapt2 in lsx, dice.
c        Re-adapt surface model layer "thicknesses" (column masses),
c        resetting to dmvadapt(1:nlay) (kg m-2).
c        Adjusts profiles of fields so vertical integral of
c        relevant quantity is conserved.
c
c        Only re-adapt layers 1 to nlay. Layer nlay+1 is passed ("*4")
c        representing interior model ice below layer nlay 
c        (always pure ice, top layer of ice model).
c
c        All "thicknesses" and "depths" here (h,hold,z,zold,etc) 
c        are column masses (kg/m2).

      dimension dmvadapt(nlay)

      dimension 
     *  snow(nlay), snowold(nlay+1),
     *  sliq(nlay), sliqold(nlay+1), 
     *  t(nlay),    told(nlay+1),
     *  rho(nlay),  rhoold(nlay+1),
     *  rgr(nlay),  rgrold(nlay+1) 

      dimension
     *  h(nlay),    hold(nlay+1),
     *  z(nlay+1),  zold(nlay+2)

      dimension
     *  sumsno(nlay), sumliq(nlay),
     *  sumhea(nlay), sumspe(nlay), sumrho(nlay), sumrgr(nlay)

      nlayp = nlay+1

      do k=1,nlay
        snowold(k) = snow(k) 
        sliqold(k) = sliq(k) 
        told(k)    = t(k) 
        rhoold(k)  = rho(k) 
        rgrold(k)  = rgr(k) 
      enddo
      snowold(nlayp) = snow4 
      sliqold(nlayp) = sliq4
      told(nlayp)    = t4
      rhoold(nlayp)  = rho4 
      rgrold(nlayp)  = rgr4

      do k=1,nlayp
        hold(k) = max (.001, snowold(k) + sliqold(k))
      enddo

      do k=1,nlay
        h(k) = dmvadapt(k)
      enddo

c        set distances from surface to old layer boundaries

      zold(1) = 0.
      do k=1,nlayp
        zold(k+1) = zold(k) + hold(k)
      enddo

c     (z not used)
      z(1) = 0.
      do k=1,nlay
        z(k+1) = z(k) + hold(k)
      enddo

c        Integrate old profiles over each new layer, to get new fields
c
      zb = 0.
c~~~~~~~~~~~~~~~~
      do k=1,nlay
c~~~~~~~~~~~~~~~~
        za = zb
        zb = za + h(k)

        sumsno(k) = 0.
        sumliq(k) = 0.
        sumhea(k) = 0.
        sumspe(k) = 0.
        sumrho(k) = 0.
        sumrgr(k) = 0.

c       For snow and liquid, mass is conserved, sum(m), 
c       i.e. sum [sliq,sliq]old(ko)*dz/hold(ko).
c       For temperature, heat is conserved, sum(c*m*T) /sum (c*m).
c       For snow density, grain size, treat as tracers, sum(r*m)/sum(m).
c       See notes 3/2/17.

        do ko=1,nlayp
          if (za.lt.zold(ko+1) .and. zb.gt.zold(ko)) then
            dz = min (zold(ko+1),zb) - max(zold(ko),za)
            sumsno(k) = sumsno(k) + snowold(ko)*dz/hold(ko)
            sumliq(k) = sumliq(k) + sliqold(ko)*dz/hold(ko)

            cheat = (cheaticea + cheaticeb*(told(ko)-tmelt))
     *                         * snowold(ko)
     *            +  cheatliq  * sliqold(ko)
            sumhea(k) = sumhea(k) + cheat*told(ko)*dz/hold(ko)
            sumspe(k) = sumspe(k) + cheat*dz/hold(ko)  ! approx, old-new

            sumrho(k) = sumrho(k) + rhoold(ko)*snowold(ko)*dz/hold(ko)
c           or with volumes conserved:  
c           sumrho(k) = sumrho(k) + (snowold(ko)*dz/hold(ko))/rhoold(ko)
            sumrgr(k) = sumrgr(k) + rgrold(ko)*snowold(ko)*dz/hold(ko)
          endif
        enddo

        snow(k) = sumsno(k)
        sliq(k) = sumliq(k)
        t(k)    = sumhea(k) / max (.001, sumspe(k))
        rho(k)  = sumrho(k) / max (.001, sumsno(k))
c       rho(k)  = sumsno(k) / max (.001, sumrho(k))     ! vols conserved
        rgr(k)  = sumrgr(k) / max (.001, sumsno(k))

        if (abs(snow(k)+sliq(k)-h(k)).gt..01) then
          write (6,'(/2a,i3, 2(3x,a,f12.3))') '*** Error (vadapt3): ',
     *      'k=',k, 's+l+i=',snow(k)+sliq(k), 'h=',h(k)
          write (6,'(a,91f12.3)') '  hold   =', (hold(m),     m=1,nlayp)
          write (6,'(a,90f12.3)') '  h      =', (h(m),        m=1,nlay)
          write (6,'(a,92f12.3)') '  zold   =', (zold(m),   m=1,nlayp+1)
          write (6,'(a,91f12.3)') '  z      =', (z(m),        m=1,nlayp)
          write (6,'(a,91f12.3)') '  told   =',(told(m)-tmelt,m=1,nlayp)
          write (6,'(a,91f12.3)') '  snowold=', (snowold(m),  m=1,nlayp)
          write (6,'(a,91f12.3)') '  sliqold=', (sliqold(m),  m=1,nlayp)
          write (6,'(a,90f12.3)') '  snow   =', (snow(m),     m=1,nlay)
          write (6,'(a,90f12.3)') '  sliq   =', (sliq(m),     m=1,nlay)
          stop
        endif
c~~~~~~~~~~
      enddo
c~~~~~~~~~~

      return
      end

c/////
#endif
c/////

c-----------------------------------------------------------------------

c\\\\\\\\\\\\\\\\\\\\\
#if ! defined (VADAPT)
c\\\\\\\\\\\\\\\\\\\\\

      subroutine layermerge (snowm, sliqm, tfirn, rhosno, rgrain,
     *                       ka, n)

c     Calculates properties of two merged layers (ka, ka+1),
c     replaces ka with merged layer, shift all lower layers up by 1,
c     reduces n by 1.

#if defined (STANDALONE6)
#include <comparams.h>
#else
#include <comicephys.h>
#endif

      dimension
     *  snowm(n), sliqm(n), tfirn(n), rhosno(n), rgrain(n)

      kb = ka + 1

c        Similar to vadapt3 in dosurfmodel5 and earlier

      snowmm = snowm(ka) + snowm(kb)
      sliqmm = sliqm(ka) + sliqm(kb)

      cheata = (cheaticea + cheaticeb*(tfirn(ka)-tmelt)) 
     *         * snowm(ka)
     *         + cheatliq * sliqm(ka)
      cheatb = (cheaticea + cheaticeb*(tfirn(kb)-tmelt)) 
     *         * snowm(kb)
     *         + cheatliq * sliqm(kb)
      tfirnm = (cheata*tfirn(ka) + cheatb*tfirn(kb)) / (cheata+cheatb)

      rhosnom = (rhosno(ka)*snowm(ka) + rhosno(kb)*snowm(kb)) 
     *          / (snowm(ka)+snowm(kb))

      rgrainm = (rgrain(ka)*snowm(ka) + rgrain(kb)*snowm(kb))
     *          / (snowm(ka)+snowm(kb))

      snowm(ka) = snowmm
      sliqm(ka) = sliqmm
      tfirn(ka) = tfirnm
      rhosno(ka) = rhosnom
      rgrain(ka) = rgrainm

      if (kb.lt.n) then
        do k=kb,n-1
          snowm(k) = snowm(k+1)
          sliqm(k) = sliqm(k+1)
          tfirn(k) = tfirn(k+1)
          rhosno(k) = rhosno(k+1)
          rgrain(k) = rgrain(k+1)
        enddo
      endif
      n = n -1 

      return
      end

c\\\\\
#endif
c\\\\\

c-----------------------------------------------------------------------

      subroutine zeromass (snowm, sliqm, poros, rhosno, thick, fliq,
     *                     runofzer, nlay)

c     Kill off very small amounts of snow, liquid (individually)

      dimension snowm(nlay), sliqm(nlay), poros(nlay), rhosno(nlay),
     *          thick(nlay), fliq(nlay), runofzer(nlay)

      parameter (eps=1.e-6)                                      ! kg/m2

      do k=1,nlay
        if (snowm(k).lt.eps) then 
          snowm(k) = 0.
c         if snowm -> 0 with finite sliqm, give water to runoff: 
          runofzer(k) = runofzer(k) + sliqm(k)
          sliqm(k) = 0.
          poros(k) = 0.
          rhosno(k) = 0.
          thick(k) = 0.
          fliq(k) = 0.
        endif
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine ssr2tfp (thick, fliq,
     *                    snowm, sliqm, 
     *                    poros,
     *                    rhosno, rhoh2o, rhoice, nlay)

c     Calculate [thick,fliq,poros] from [snowm,sliqm,rhosno]. 

      dimension thick(nlay), fliq(nlay),
     *          snowm(nlay), sliqm(nlay),
     *          poros(nlay),
     *          rhosno(nlay)

      do k=1,nlay

c       rhosnomax limit is to still accomodate sliqm with fliq=1:
        rhosnomax = 1. / ( (sliqm(k)/max(snowm(k),.001))/rhoh2o 
     *                          + 1./rhoice )
        rhosno(k) = min (rhosno(k), rhosnomax)

        poros(k) = max (0., 1. - rhosno(k)/rhoice)

        thick(k) = snowm(k)/rhosno(k)

        fliq(k) = sliqm(k) / max ( rhoh2o*poros(k)*thick(k), 1.e-8 )

        fliq(k) = min (1., fliq(k))

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine zerobot (snowm, sliqm, tfirn, rhosno, rgrain,
     *                    n, nmax)

c     For nicety, zero prognostic variables below current bottom (n)

      dimension
     *  snowm(nmax), sliqm(nmax),
     *  tfirn(nmax), rhosno(nmax), rgrain(nmax)

      if (n.lt.nmax) then
        do k=n+1,nmax
          snowm(k) = 0.
          sliqm(k) = 0.
          tfirn(k) = 0.
          rhosno(k) = 0.
          rgrain(k) = 0.
        enddo
      endif

      return
      end

c-----------------------------------------------------------------------

c}}}}}}}}}}}}}}}}}}}}}
#if defined (CROSSSEC)
c}}}}}}}}}}}}}}}}}}}}}

      subroutine crosssec6 (zday, zin, valin, nz, nt, clabin,
     *                      iflag, timeice, nyearfirn, iu)

c     For cross-section output (versus z and day). Processes
c     one variable (valin, named clabin). Converts from
c     column-mass vertical coord to depth(m) from surface,
c     and writes to output file (for gp load.crosssec).

      dimension zday(nt), zin(nz,nt), valin(nz,nt)
      character*(*) clabin
      character*80 clab, cfile
#if defined (MEYERHEWITT)
      parameter (nzo=500, tzo=20.)   ! tzo is output max z below sfc (m)
#else
c     parameter (nzo=500, tzo=200.)  ! tzo is output max z below sfc (m)
c     parameter (nzo=500, tzo=40.)   ! tzo is output max z below sfc (m)
      parameter (nzo=500, tzo=20.)   ! tzo is output max z below sfc (m)
#endif
      parameter (ntmax=365*2)

      dimension 
     *  zout(nzo), indz(nzo,ntmax), jndz(nzo,ntmax), weiz(nzo,ntmax),
     *  valout(nzo,ntmax) 
      save zout, indz, jndz, weiz

      logical firstcrosssec6
      data firstcrosssec6 /.true./
      save firstcrosssec6

      if (firstcrosssec6) then
        if (nt.gt.ntmax) then
          write (6,'(a,2(a,i6))') '*** Error (surfmodel6 crosssec6):',
     *                            '  nt=',nt, '  ntmax=',ntmax
          stop
        endif
        do ko=1,nzo 
          zout(ko) = (ko-.5) * (tzo/nzo)
        enddo
        firstcrosssec6 = .false.
      endif

      clab = clabin

c        Calculate indices from column-mass vertical coord (zin, nz)
c        to depth (m) from surface (zout, 1:nzo, 0 to tzo). Do for
c        first field called for this day of year (iflag=1, passed),
c        or for wvel, that has different vertical coords (zvel, passed).
 
c.............................................
      if (iflag.eq.1 .or. clab.eq.'wvel') then
c.............................................
        do m=1,nt
          do ko=1,nzo
            if (zout(ko).le.zin(1,m)) then
              indz(ko,m) = 1
              jndz(ko,m) = 1
              weiz(ko,m) = 1.
            else if (zout(ko).ge.zin(nz,m)) then
              indz(ko,m) = nz
              jndz(ko,m) = nz
              weiz(ko,m) = 1.
            else
              do k=2,nz 
                if (zout(ko).le.zin(k,m)) then 
                  indz(ko,m) = k-1
                  jndz(ko,m) = k
                  weiz(ko,m) = (zin(k,m)-zout(ko))
     *                       / max(zin(k,m)-zin(k-1,m), .001)
                  goto 100
                endif
              enddo
              write (6,*) '*** Error (surfmodel6 crosssec6): not here' 
              stop
  100         continue
            endif
          enddo
        enddo
c..........
      endif
c..........

c        Interpolate field

      do m=1,nt
        do ko=1,nzo 
          valout(ko,m) =     weiz(ko,m) *valin(indz(ko,m),m) 
     *                 + (1.-weiz(ko,m))*valin(jndz(ko,m),m) 
        enddo
      enddo

      cfile = clab(1:lenchr(clab)) // '.cross'
      open (iu, file=cfile, status='unknown')
      close (iu, status='delete')
      open (iu, file=cfile, status='new')

      write (iu,'(2a, a,f10.2, a,i6)') 
     *  'CROSSSEC: ',    clab(1:lenchr(clab)), 
     *  '  timeice= ',   timeice,
     *  '  nyearfirn= ', nyearfirn

c-----------------------------
      if (clab.eq.'wvel') then
c-----------------------------

c          For vertical velocity, calculate streamlines by lagrangian
c          integration through year, oputput to wvel.cross for gp
c          (load.crosssec)

        call streamline6 (zout, valout, tzo, nzo, nt, ntmax,
     *                    timeice, nyearfirn, iu)
c---------
      else
c---------
c          For other fields, 
c          output 2-D array to file (*.cross) for gp (load.crosssec). 
c          day 1:nt, (nt,passed,= nspycross), 0 to 365.
c          depth 1:nzo, 0 to tzo (m).

        write (iu,'(a8,2x,a8/i8,2x,730f8.1)') '    z(m)','     day',
     *    nt, (zday(m),m=1,nt)
        do ko=1,nzo
          write (iu,'(f8.2,2x,730f8.2)') zout(ko), (valout(ko,m),m=1,nt)
        enddo
c----------
      endif
c----------

      close (iu)

      if (clab.eq.'temperature') then
        cfile = 'gfinfo.cross'
        open (iu, file=cfile, status='unknown')
        close (iu, status='delete')
        open (iu, file=cfile, status='new')
        write (iu,'(a,i6)') 'tzo= ', nint(tzo)
        close(iu)
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine streamline6 (zo, wo, tzo, nzo, nto, ntmax,
     *                        timeice, nyearfirn, iu)
c          Calculate streamlines, by lagrangian integration through 
c          year, using 2-D fields of vertical velocity (wo) vs.
c          day and depth from surface (zo). Output to wvel.cross for 
c          gp (load.crosssec)

      dimension zo(nzo,ntmax), wo(nzo,ntmax)

      parameter (nplag=50)  ! lagrangian paths, with various depths
c                             at start of year. Also various times 
c                             during year, starting af sfc or bottom.

c     either:
c     parameter (ntlag=100) ! # of lagr time pts through year (<= nto!)
c     dimension zlag (0:ntlag, nplag, 3) ! lagrangian path depths
c     or:
      dimension zlag (0:nto, nplag, 3)   ! lagrangian path depths
      ntlag = nto

      dzo = tzo/nzo

#if defined (CROSSSEC2YR)
      tto = 730.
#else
      tto = 365.
#endif
      dto = tto/nto
      dtlag = tto/ntlag

c        ibound = 1 for various depths, starting at day 0. 
c        ibound = 2 for various days, starting at surface.
c        ibound = 3 for various days, starting at bottom.

c>>>>>>>>>>>>>>>>>>
      do ibound=1,3
c>>>>>>>>>>>>>>>>>>

c        Loop over paths 

c=====================
      do iplag=1,nplag
c=====================

        if (ibound.eq.1) then
c         beginning of year, range over depths:
          itlaga = 0
          zlag(itlaga,iplag,ibound) = tzo*(iplag-.5)/float(nplag)

        else if (ibound.eq.2) then
c         surface, range over times:
          itlaga = nint( tto*(iplag-.5)/float(nplag) / dtlag)
          zlag(itlaga,iplag,ibound) = 0.

        else if (ibound.eq.3) then
c         bottom, range over times:
          itlaga = nint( tto*(iplag-.5)/float(nplag) / dtlag)
          zlag(itlaga,iplag,ibound) = tzo
        endif

        if (itlaga.gt.0) then
          do m=0,itlaga-1
            if (ibound.eq.2) zlag(m,iplag,ibound) = -1.
            if (ibound.eq.3) zlag(m,iplag,ibound) = tzo + 1.
          enddo
        endif

c................................
        if (itlaga.lt.ntlag) then
c.................................

c            Loop over lagrangian points for this path.
c            At each point, bilinear interpolation in (z,t) space
c            for vertical velocity (wcur). 

          zcur = zlag(itlaga,iplag,ibound)
c--------------------------------
          do itlag=itlaga+1,ntlag 
c--------------------------------
            tmid = (itlag-.5)*dtlag                 ! mid-interp in time

            iza = nint(zcur/dzo)
            izb = iza + 1
            wza = izb - (zcur/dzo + .5)
            iza = max (1, min (nzo, iza))
            izb = max (1, min (nzo, izb))

            ita = nint(tmid/dto)
            itb = ita + 1
            wta = itb - (tmid/dto + .5)
            if (ita.eq.0) ita = nto
            if (itb.eq.nto+1) itb = 1

            wcur =     wza *   wta  * wo(iza,ita)
     *           + (1.-wza)*   wta  * wo(izb,ita)
     *           +     wza *(1-wta) * wo(iza,itb)
     *           + (1.-wza)*(1-wta) * wo(izb,itb)

            znext = zcur + wcur*dtlag*86400.            ! wcur is in m/s

            if (znext.lt.0.) then
c             if path escapes top, set remaining pts in the year to -1:
              do m = itlag,ntlag 
                zlag(m,iplag,ibound) = -1.
              enddo
              if (itlag.eq.itlaga+1) zlag(itlaga,iplag,ibound) = -1.  
              goto 100
            else if (znext.gt.tzo) then
c             if path escapes bottom, set rem pts in the year to tzo+1:
              do m = itlag,ntlag 
                zlag(m,iplag,ibound) = tzo + 1.
              enddo
              if (itlag.eq.itlaga+1) zlag(itlaga,iplag,ibound) = tzo +1.
              goto 100
            endif

            zlag(itlag,iplag,ibound) = znext

            zcur = znext 
c--------------
          enddo
c--------------
  100     continue

c............
        endif
c............

c==========
      enddo
c==========

c>>>>>>>>>>
      enddo
c>>>>>>>>>>

c       Write paths as columns (row is time through year,
c       each column is a path), grouped by ibound=1,2,3,
c       to wvel.cross, for plotting (gp, load.crosssec).

      do ibound=1,3
        write (iu,'(/a,i3)') 'ibound=',ibound
        write (iu,'(a)') '     day'
        do m=0,ntlag
          write (iu,'(f8.2, 4x, 100f8.2)') 
     *      m*dtlag, (zlag(m,iplag,ibound), iplag=1,nplag)
        enddo
      enddo

      return
      end

c}}}}}
#endif
c}}}}}

c=======================================================================

c=====
#endif
c=====
