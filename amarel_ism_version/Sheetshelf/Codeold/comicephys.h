c comicephys.h

c        Basic quantities:

      parameter (tmelt  = 273.15)
c     parameter (hfus   = 0.3336e6)
      parameter (hfus   = 0.335 e6)
c     air-temp cutoff for rainfall vs. snowfall (in iceclimall.F):
#if defined (TRAINSNOW)
      parameter (trainsnow = tmelt + TRAINSNOW)
#else
      parameter (trainsnow = tmelt)
#endif

      parameter (rhobed = 3370., rhosed = 2390.) ! kg/m3

#if defined (EISLINE) && defined (LINEM)
      parameter (rhoice =  900.)
      parameter (rholiq = 1000.)
      parameter (grav   = 9.80)
#elif defined (NHA) && defined (GREENLAND) 
      parameter (rhoice =  917.)
      parameter (rholiq = 1027.)
      parameter (grav   = 9.81)
#elif defined (NHA) && ! defined (GREENLAND) 
      parameter (rhoice =  910.)
      parameter (rholiq = 1028.)
      parameter (grav   = 9.81)
#elif defined (EISLINE) || defined (EISANTA) || defined (CARB) || defined (ORDO)
      parameter (rhoice =  910.)
      parameter (rholiq = 1028.)
      parameter (grav   = 9.81)
#elif defined (TEST2D)
c     for mismipplus:
      parameter (rhoice =  918.)
      parameter (rholiq = 1028.)
      parameter (grav   = 9.81)
#elif defined (TESTMEL2D)
      parameter (rhoice =  918.)
      parameter (rholiq = 1028.)
      parameter (grav   = 9.81)
#else
      parameter (rhoice =  910.)
      parameter (rholiq = 1000.)
      parameter (grav   = 9.80616)
#endif
      parameter (rhor = rhoice/rholiq)
      parameter (rhoip = (1.-rhoice/rholiq)*rhoice)

c        Ice rheology, internal deformation:

      integer powi
c     parameter (powi = 1)   ! can't use n=1, only n=3 for shear 
      parameter (powi = 3)   ! softening in icedyn, icetherm

#if defined (TEST2D)
c     for mismipplus:
#  if defined (CRHEOLPLUS)
      parameter (crheolmismipplus=CRHEOLPLUS)
#  else
c     parameter (crheolmismipplus=2.0e-17)
      parameter (crheolmismipplus=3.5e-17)
#  endif
#endif

c     for shelf:
#if defined (LINSHELF) 
      parameter (powiv = 0.)
      parameter (powir = 1.)
#else
      parameter (powiv = (powi-1.)/(2.*powi))
      parameter (powir = (1./powi))
#endif
      parameter (powih = (powi+1.)/(2.*powi))      ! horiz shear heating

c        Basal (non-sed) sliding rheology:

      integer powb
#if defined (EISLINE) && defined (LINEF)
      parameter (powb = 2)
c     parameter (powb = 3)   ! schoof
#elif defined (EISLINE) && defined (LINEM)
#  if defined (BASALa)
      parameter (powb = 3)
#  elif defined (BASALb)
      parameter (powb = 1)
#  endif
#elif defined (TEST2D) 
c     for mismipplus:
      parameter (powb = 3)
#else
c     parameter (powb = 1)
      parameter (powb = 2)
c     parameter (powb = 3)
#endif

c~~~~~~~~~~~~~~~~~~~
#if defined (MELANGE)
c~~~~~~~~~~~~~~~~~~~
c     parameter (rhomel = rhoice)
      parameter (rhomel = 930.)
      parameter (rhomp = (1.-rhomel/rholiq)*rhomel)
      parameter (rhormel = rhomel/rholiq)

      integer powimel
#  if defined (POWIMEL) 
      parameter (powimel = POWIMEL)
#  else
c     parameter (powimel = 3)
c     parameter (powimel = 2)
      parameter (powimel = 1)
#  endif

      integer powbmel
#  if defined (POWBMEL) 
      parameter (powbmel = POWBMEL)
#  else
c     parameter (powbmel = 3)
c     parameter (powbmel = 2)
      parameter (powbmel = 1)
#  endif

      parameter (powivmel = (powimel-1.)/(2.*powimel))
      parameter (powirmel = (1./powimel))
      parameter (powihmel = (powimel+1.)/(2.*powimel)) ! horiz shear htg

#  if defined (SIDEDRAGMEL) 
      parameter (sidedragmel = SIDEDRAGMEL)
#  else
      parameter (sidedragmel = 1.)
#  endif

#  if defined (COEFDIVMEL) 
      parameter (coefdivmel = COEFDIVMEL)            ! for melange rheol
#  else
      parameter (coefdivmel = 0.)                    ! for melange rheol
#  endif

#  if defined (COEFCONMEL) 
      parameter (coefconmel = COEFCONMEL)            ! for melange rheol
#  else
      parameter (coefconmel = .05)                   ! for melange rheol
#  endif

#  if defined (COEFBWATERMEL) 
      parameter (coefbwatermel = COEFBWATERMEL) ! Pa/(m/a)
#  else
      parameter (coefbwatermel = 1.e-7)         ! Pa/(m/a)
#  endif

      parameter (hminmel = .5)

#  if defined (SLOPEMAXMEL)
      parameter (slopemaxmel = SLOPEMAXMEL)
#  else
      parameter (slopemaxmel = .01)   ! 10 m/km
c     parameter (slopemaxmel = 1.e20)
#  endif

#  if defined (HMAXMEL)
      parameter (hmaxmel = HMAXMEL)
#  else
      parameter (hmaxmel = 1.e20) 
c     parameter (hmaxmel = 200.) 
#  endif

c~~~~~
#endif
c~~~~~

c        Sediment rheology (other parameters set in sedflow)

c     parameter (pows  = 1.25)
      integer pows
      parameter (pows = powb - 1) ! so sed and basal can combine crh's

c        Avoid singular matrix if floating areas.
c        Also used for water friction (coefbsu*(ub-uw),coefbsv*(vb-vw))
c        in icedyn, with uw,vw set in movewater, to mimic ice-dam
c        breaking in catastrophic floods.

#if defined (EISLINE)
      parameter (coefbwater = 0.)
#else
      parameter (coefbwater = .001)  ! Pa/(m/a)
#endif


c        Water thickness (hw, m) above which considered lake or ocean

#if defined (NOMOVET) && ! defined (MOVEW)
      parameter (hwcut = 0.)
#else
c     parameter (hwcut = 0.2)   ! can't be zero in movewater 
      parameter (hwcut = 0.)    ! can't be zero in movewater 
#endif

c        Geothermal heat flux (J/m2/a, ~0.9 ucal/cm2/s):

#if defined (EISANTA)
      parameter (geoflux_eais = .0546 * 31556926)
#  if defined (GEOFLUX_W)
      parameter (geoflux_wais = (GEOFLUX_W/1000.) * 31556926)
#  else
      parameter (geoflux_wais = .070  * 31556926)
#  endif
#elif defined (EISLINE)
      parameter (geoflux_unif = .0711 * 31556926)
#elif defined (NHA) || defined (CARB) || defined (ORDO)
      parameter (geoflux_unif = .042  * 31556926)
#else
      parameter (geoflux_unif = .0546 * 31556926)
c     parameter (geoflux_unif = (0.9*4.187*.01)*31556926)
#endif

c        Rate of decrease of ice pressure-melting point with depth (K/m)

      parameter (dtmdh = 8.66e-4)
c     parameter (dtmdh = .0074*1.e-5*rhoice*9.80616)

c        Thermal conductivity, heat capacity, of ice     
c        condice  = condicea * exp(-condiceb*t)
c        cheatice = cheaticea + cheaticeb*(t-tmelt)

#if defined (EISLINE) || defined (EISANTA) || defined (NHA) || defined (CARB) || defined (ORDO)
      parameter (condicea = 2.1*31556926)                     ! J/a/m/K
      parameter (condiceb = 0.)                               ! J/a/m/K
      parameter (cheaticea= 2009.)                            ! J/kg/K
      parameter (cheaticeb= 0.)                               ! J/kg/K^2
#else
      parameter (condicea = 3.10e8)                           ! J/a/m/K
      parameter (condiceb = .0057)                            ! 1/K
      parameter (cheaticea= 2115.3)                           ! J/kg/K
      parameter (cheaticeb= 7.79)                             ! J/kg/K^2
#endif

      parameter (condsed  = 3.3*31556926)                     ! J/a/m/K
      parameter (cheatsed = 1000.)                            ! J/kg/K

      parameter (condbed  = 3.3*31556926)                     ! J/a/m/K
      parameter (cheatbed = 1000.)                            ! J/kg/K

      parameter (condliq = 70.*31556926)                      ! J/a/m2/K
      parameter (cheatliq = 4218.)                            ! J/kg/K

      parameter (condair = 10.*31556926)                      ! J/a/m2/K

c        lapse rate for surface-air temperature corrections (+ve, K/m)
c        (in iceclim[all,gcm].F)

#if defined (LAPSEPARAM)
      parameter (rlapse = LAPSEPARAM)
#elif defined (RLAPSE)
      parameter (rlapse = RLAPSE)
#else
#  if defined (NHA)
      parameter (rlapse = .0050)
#  else
      parameter (rlapse = .0080)
#  endif
#endif

c        temperature scaling for "zdt" precip corrections (C)

#if defined (TLAPSEPRECIP)                   
c     for elevation differences:
      parameter (tlapseprecip = TLAPSEPRECIP)          ! 0 for no effect
#else
      parameter (tlapseprecip = 10.)
#endif

#if defined (TLAPSEPRECIPCLIM)           
c     for parameterized climate shifts:
      parameter (tlapseprecipclim = TLAPSEPRECIPCLIM)  ! 0 for no effect
#else
      parameter (tlapseprecipclim = 10.)
#endif

#if defined (NHA)
      parameter (alorb =  70.) ! latitude for insol-based weighting (oN)
#else
      parameter (alorb = -70.) ! latitude for insol-based weighting (oN)
#endif

c        For sediment domain-wide diagnostic budget

      common /sedbud/ 
     *  totquar,  totpelag,  totdump,  totslump, 
     *  totquara, totpelaga, totdumpa, totslumpa,
     *  totsed, totsedprev, timesedprev

c        Various basal/sed coefficients, units are (m/a) / Pa^powb

      parameter (crhnos = 1.e-15)                 ! ~no sliding (frozen)
      parameter (crhhard= 1.e-10)                 ! hard bedrock
c     parameter (crhsed = 1.e-5)                  ! deformable sediment 
      parameter (crhsed = 1.e-4)                  ! deformable sediment 

c        For Last Interglacial atmos and ocean perturb forcing (as DP16)

#if defined (ATMOSDT_LIG) || defined (OCEANDT_LIG)
      parameter (nlig=4)
      dimension timlig(nlig), dtalig(nlig), dtolig(nlig)
      save timlig, dtalig, dtolig
      data timlig /-130000., -125000., -115000., -110000. /
      data dtalig /    1.97,     1.41,     0.83,    -1.57 /
      data dtolig /    1.70,     1.51,     1.09,     0.31 /
#endif

