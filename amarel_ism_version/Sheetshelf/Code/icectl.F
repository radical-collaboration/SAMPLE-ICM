c-----------------------------------------------------------------------

c List of preprocessor "define" keywords (-DKEYWORD in makefile). 
c   If indented like this, need to also use the keyword above.
c Stars (*XXX*) indicates a numeric value is expected, e.g., -DXXX=1.23

c=================
c Grid projections (must specify appropriate one for selected domain):
c=================

c LONLAT = longitude vs. latitude
c   WRAPAROUND = longitude grid wraps around globally (-180 to 180 E)
c STEREO = polar stereographic (ellipsoids and projections set in
c          comiceprojection.h, transformations in iceprojection.F)

c=========
c Domains:
c=========

c 2-D domains:

c EISANTA    = Antarctic 2-D (or 1-D if TRANSECTA) (using STEREO)
c  *TRANSECTA* = 1-D Antarctic transect (initgrid/settransect),
c                PIG, ROSS, WEDDELL, WILKES, or SIPLE
c     ZOOMGL = zoomed-grid around grounding line
c   NESTING    = nested domain (more parameters below)
c NHA        = 2-D Northern Hemisphere (using LONLAT)
c   GREENLAND  = Greenland domain (using LONLAT)
c   GLACBAY    = Glacier Bay domain (using LONLAT)
c CARB       = 2-D global (Permo-Carb) (using LONLAT)
c ORDO       = 2-D global (Ordovician) (using LONLAT) (w. CSMORDO)
c RESOL[40,20,10,5] = resolutions for EISANTA (km) or GREENLAND (degrees
c                     latitude*100). Defaults=40 (Ant.) and 20 (Green.)
c TEST2D     = 2-D idealized plane, divide at lhs, h(nx)=0(using STEREO)
c              for mismipplus tests.
c   TEST2DOUT = write requested output files for mismipplus.
c TESTMEL2D = 2-D idealized rectangle, with MELANGE.

c 1-D domains:

c EISLINE    = 1-D flowlines (idealized tests), with one of the options
c              below. Uses STEREO grid, except uses LONLAT for SLMODEL.
c   LINEB      = Quaternary sawtooth cycles (weirun: 0=modern, 1=lgm)
c   LINEC      = basic linear terrestrial slope (no marine), resol tests
c                with SUBGRIDEDGE (cf. van den Berg, J.Glac, 2006)
c     SINBED     = sinusoidal bed   (cf. van den Berg, J.Glac, 2006)
c     FLATANT    = flat bed, flat ELA, bounded by ocean
c   LINED      = basic mountain/ocean profile, quarrying+transport
c   LINEE      = wedge simulations, time independent or dependent.
c     Either INITEE, RESTARTE, SEALEVCHANGE (time independent),
c     or none of those (time dependent):
c     INITEE     = no wedge, time indep.
c     RESTARTE   = with wedge, time indep.
c     WEDGEnnKM  = wedge extent, nn = 4,8,12,16 (km), time indep.
c     SEALEV[CHANGE,RISE,FALL] = instantaneous s.l. change, time indep.
c   LINEF      = grounding line behavior, with divide at l.h. edge
c   LINEG      = sub-glacial lakes and hydrology  
c   LINEH      = Cenozoic "E-O", no marine, EAIS cross section.
c                Use with PARAMCLIM and SIMPEO (lon,lats set in iceinit)
c   LINEH2     = like LINEH (no marine), except simple bedrock profile 
c   LINEH3     = like LINEH2 except w. marine, no PARAMCLIM, time=0-> 
c   LINEH4     = like LINEH3 except 2-D axisym.
c   LINEM      = for MISMIP-like tests (D.Docquier, F.Pattyn, 12/2010).
c                Also uses BASAL[a,b], RESOL[n] (and not yet EXPT[n]).

c====================================================================
c Long-term sea level, orbit/insol, CO2 (sealev and orbit mostly from 
c obs time series, in subr longterm), and ocean melt (in subr ocean):
c====================================================================

c FORCEPLEIST = long-term forcing for EISANTA Plio-Pleistocene runs. 
c               lisiecki, or spratt if -DSPRATT (sealev),       
c               laskar orbit (dtant[ann,jan], dtseas), 
c               or hard-coded in longterm.
c               Also sets weirun for subr ocean. Nb:
c               dtantann,dtseas,ro18,rco2 affect atm.T,P (iceclimparam)
c               dtantjan,rco2,ro18 affect weirun (longterm -> ocean), 
c               but weirun (-> zclim,zwhot/mod/lgm) has only minor 
c               effects in subr ocean (via LIUCORR2, DTOCNMATCO2,
c               PITFUTOCN, CALVLIQZWHOT).
c FORCEEO     = long-term forcing for EISANTA Eocene-Oligocene runs.
c                sealev=0, laskar orb (dtant[ann,jan], dtseas).
c FORCEFUTRAMP = simple linear-ramped warming for years > 0.
c *PITFUTCO2*  = overrides rco2 in longterm, future ramp from 1 to 
c                PITFUTCO2.
c *PITFUTOCN*  = overrides weirun in longterm, future (weirun->zwhot)
c                ramp of ocean warming from 0 to zwhot*PITFUTOCN (ocean)
c    PITFUTOCNASE = only apply in 'West' sector (csec, Amundsen Sea).
c *DTOCNMATCO2* = ocean warming ramp (C), 0 to zwhot*DTOCNMATCO2 (ocean)
c    DTOCNMATCO2ASE = only apply in 'West' sector (csec, Amundsen Sea).
c ASYNCH      = long-term asynchronous runs using runasyn... scripts.
c               (sealev=0). (Only with GCMMATRIX)

c ORBPRIME    = synthetic (not real) orbits returned by subr laskar,
c               with prime numbers for prec,obl,ecc periods
c ORBMED      = fixed "me" (medium) orbit imposed, only with GCMMATRIX
c AUSTRALSPRING = 8/21 to 11/20 insolation used in laskar for dtantjan
c               (default is 12/21 to 1/20)  
c*ALATLASK*   = latitude to compute insolation in laskar (SH<0, def=-80)

c CO2INTER    = interactive CO2 due to silicate weathering, computed 
c               using a soil bucket model (on both ice-model and global
c               GCM grids) and Si weathering = f(runoff), in iceco2.F 
c               (Lee Kump) (only with EISANTA, (PARAMCLIM or GCMMATRIX))

c SL_ICE5G    = override sealev with ICE5G eustatic for last 40 ka
c               (only with EISANTA, or NHA+LOVECLIP)

c=================
c Climate forcing:
c=================

c GCMMATRIX = surface budget from GCM meteo file(s) with matrix method.
c   NHAOBS  = use modern observed climatol (from Genesis data files)
c             (with NHA).
c   NHAGLAC = use Genesis data file matrix (gen2) for glac-interglac,
c             (with NHA).
c RCMMATRIX = surface budget from RCM meteo file(s) with matrix method.
c  *RCMFILE*   = file name of single RCM climate file (if none, modern).
c   RCMMATCO2 = or: 1-D matrix interp based on rco2 (btwn rcm files for
c               mod,2x,4x). rco2+weirun set vs. time in subr longterm.
c     Either: 
c     RCPSPIN*, RCP[26*,45*,85*][FIX], COP* = like ICCP scenarios 
c               RCP 2.5,4.5,8.5, or COP21-based scenarios.
c     Or:
c    *RAMPMATCO2* = linear ramp time span, default = step fn (yr)
c       RAMP[PLIO,2X,4X,8X,8XW]=ramp from RCM pre_40 to plio,2x,4x,8x[w]
c     These both ramp future CO2(t) via rco2. 
c     RCP[26,45,85] determines future ocean temps (see below).  
c     RAMPMATCO2 ramps ocean temps via weirun(->zwhot) and DTOCNMATCO2.

c  *RCMCO2*   = rco2 level set in subr longterm, if not RCMMATCO2
c               (just for display - climate is from RCMFILE).
c  *RCMWEI*   = weirun level set in subr longterm, if not RCMMATCO2
c               (weirun determines ocean forcing in iceocean).
c   RCMANOM = anomaly method for T,P (T:exp+obs-ctl, P:exp*obs/ctl) 
c RACMOCLIM = surface budget from RACMO2 modern meteo files 
c             for Greenland (w. GREENLAND) or Antarctica (w. EISANTA)
c LOVECLIP = coupling (via files, separate programs) with LOVECLIM EMIC 
c   LOVEANOM = anomaly method for T,P,To (T:exp+obs-ctl, P:exp*obs/ctl).
c              nb: LOVEANOM for Ant only (EISANTA) - no obs avail for NH
c CSMCLIM = surface budget and ocean temps from CSM.
c   CSMANOM = anomaly method for T,P (T:exp+obs-ctl, P:exp*obs/ctl) 
c CLEMCLIM = surface budget from Clemens Schannwell CMIP GCM.
c CSMORDO  = surface budget from C. Shields CSM Ordovician (w. ORDO).
c EBM       = surface budget from my EBM (climate_ebm). 
c             Not currently maintained - would need work to resurrect.
c PARAMCLIM = parameterized climate, with one of:
c   SEARISECLIM = SeaRise Antarc modern datasets + param. perturbs.
c                 By  default, uses new Lebrocq modern data. 
c     NOLEBROCQ = use older non-Lebrocq data (also affects SEARISEINIT).
c     BASACCUM    = use BAS (Arthern) dataset for modern accumulation
c     SEARISEAR4  = additional perturbs, AR4 future minus modern
c   SIMPANT     = simple Antarc param dep on hs,lat,dtantann,sealev,rco2
c   SIMPEO      = other simple Antarc params, with one of SIMPEO_[1-3],
c                 or SIMPEO_HUY93 (like Huybrechts, Geogr.Annal., 1993)
c   SIMPNHA       = simple Northern Hemispheric param
c   *DTPERT*     = unif (f(t)) air temp shift, w. SIMPEO_HUY93, SIMPNHA
c RACMOANTSEAS   = add seasonal cycle (departures) of RACMO2 monthly 
c                  air temps (T) to SeaRise annual mean. Only with 
c                  EISANTA, and SEARISECLIM or [RCM,CSM,LOVE]ANOM.
c---
c for UW 5 Myr climate forcing modifications (with PARAMCLIM):
c   *WARMAMPA* = multplies ro18 in atmos dT if ro18 < 0 (def=1)
c   *COLDAMPA* = multplies ro18 in atmos dT if ro18 > 0 (def=1)
c   *WARMAMPO* = multplies ro18 in weirun (ocean) if ro18 < 0 (def=1)
c   *COLDAMPO* = multplies ro18 in weirun (ocean) if ro18 > 0 (def=1)
c   *QAMPANN*  = coeff for dtantann in atmos dT or weirun (def=0.1) 
c   *QAMPSEAS* = coeff for dtantjan,dtseas " " " (def=0.1)
c   *QAMPAUS*  = coeff for dtantaus " " " (def=0.1)
c    CO2PLIO   = linear ramp of rco2, CO2PLIO(ppm)/280 to 1, 3 to 2 Ma
c---

c *LAPSEPARAM* = lapse rate for sfc elev air temp correction, |deg C|/m
c *RLAPSE*     = same as LAPSEPARAM
c TLAPSEPRECIP = scaling in precip correction P*(2**(zdt/tlapseprecip)),
c                where zdt (deg C) = lapse rate * elev difference,
c                due to topography resol. differences (iceclimall.F). 
c                (default = 10 C; set to 0 for no effect). 
c TLAPSEPRECIPCLIM = same for warming/cooling shifts in parameterized
c                    climates (iceclimparam.F).
c *ATMOSDT* = unif. atmos warming, all climates (C) (default = 0).
c
c ALBEFFTOT = ice albedo feedback on air temps, prop. to total ice area,
c             (calcalbeff), with SIMPEO, SIMPNHA, or GCMMATRIX+NHA.
c*ALBEFFGREEN* = alternative to ALBEFFTOT. Sum of cooling footprints of
c                all ice points. Value = e-folding length (m).
c
c Surface snow-ice mass budget treatment (in iceclimall), intra-annual 
c "model" for PDD/surface melt with refreezing correction. Choice of 
c particular surface snow-ice mass model (in iceclimall, default is 
c SURFMODEL3):
c SURFMODEL1 = annual cycle, explicit snow, ice, refreezing fraction.
c SURFMODEL2 = annual cycle, explicit snow and embedded liq. in pores.
c SURFMODEL3 = 2 steps, annual cycle for pdd,etc, then ann. refr corr.
c SURFMODEL6 = snow-firn column, lagr.layers, ice lens, perc, 
c              d(rhos)/dt, d(grain)/dt, etc.
c with SURFMODEL[1-6]:
c   SURFINSOL = use lin. sfc. energy eqn. with TOA insol., not PDD. 
c               Diurnal cycles of insolation are resolved if dtimed 
c               (can be set by DTIMED) is < 1 day.
c               SURFMODEL6 use SURFINSOL automatically.
c  *NYEARFIRN*= years to integrate sfc model (def=2 to 4, except for  
c               SURFMODEl6: def=2000 if equil., def=dtimebudg if progn.)
c  IFPROG6    = for SURFMODEL6 only, prognostic each call (def=equil.)
c  *DTIMED*   = intraseasonal timestep for surface bal. model, hours.
c               In comicegrid.h, used in iceclimall.F (default=5 days). 
c  *TADIURN*  = Diurnal "cycle" (sinusoid) superimposed on seasonal air 
c               temps, amplitude TADIURN reduced by diurnal insolation
c               range, lagged to 2 pm. Done only if SURFINSOL or 
c               SURFMODEL6, and dtimed < 1 day.
c  *TASYNOP*  = "synoptic" sinusoid added to seasonal air temps.
c               amplitude TASYNOP, period 10 days.
c  *SDT*      = std. dev. of temps in degree_day (default = 5).
c               nb: if sfc energy balance model used, eg, SURFMODEL6,
c               only affects split between snowfall and rainfall.
c               nb: should reduce SDT if TADIURN or TASYNOP are used. 
c  *TRAINSNOW*= air temperature (C) determining whether precip is 
c               rain or snow (default=0)
c  *TPDD*     = air temperature (C) for zero pdd cutoff (default=0)
c SIMPPDD    = alternate simpler basic PDD calculation for ice budget,
c              no refreezing (and SURFINSOL has no effect).
c
c MORESHELFMELT  = arbitrary increase in rain & sfc melt on ice shelves
c NOSURFACEWATER = arbitrary setting of no rain, no sfc melt, everywhere
c*MOREWILKESSNOW*= arbitrary increase in snowfall, ~Wilkes-Aurora 
c                  basin sectors, by factor MOREWILKESSNOW

c LAGRANGEO18 = Lagrangian tracing back to surface, report mean d18O of
c               ice sheet from d18O_precip

c==============================================================
c Initial ice, sed, bed, and ice-free rebounded bed topography:
c (for initial starts with ifrest=-1. If ifrest=0 or 1, initial 
c conditions (h,hb) read from restart file and overwritten).
c==============================================================

c EISMINT2    = Eismint II modern Antarc (with EISANTA)
c BEDMAP1     = BEDMAP  (orig) modern Antarc (with EISANTA)
c BEDMAP2     = BEDMAP2 (v2)   modern Antarc (with EISANTA)
c SEARISEINIT = SeaRise modern Antarc (with EISANTA)
c   NOLEBROCQ = use older non-Lebrocq data (also affects SEARISECLIM)
c   LEBROCQ2  = with deeper Recovery basin (ice thickness,bed depth)
c SHALDRIL    = Wilson+Luyendyk early Oligocene reconstruct (w. EISANTA)
c   SHALBEDMAP  = use old bednew4 (Bedmap1) input file
c   SHALGRL09   = use Wilson-Luyendyk GRL 2009 reconstruction
c   SHALP311MAX = use Wilson et al. P3 2012 maximum reconstruction
c   SHALP311MIN = use Wilson et al. P3 2012 minimum reconstruction
c   SHALP311AVE = use average of the above
c   SHALFILMIN  = same SHALP311MIN with TAM troughs filled in
c   DWSEDPOLY   = impose Wilson sed patch (crhmel) in Weddell Embayment
c   SHALMIOC    = use Wilson early Miocene (~18 Ma,RSU4a) reconstruction
c   SHALMIOC15  = use Gasson interpolated 15Ma recons. based on Wilson13
c AGASEA      = override above for Antarctica with AGASEA data in 
c               Amundsen Sea region (nb. already in SEARISEINIT (unless
c               NOLEBROCQ), and in BEDMAP2)
c               (With EISANTA.)
c*WANTBATHY*  = arbitrary weighting of modern rebounded and "SHALDRIL"
c               (W-L early Olig) topographies. Weighting is
c               WANTBATHY*W-L + (1-WANTBATHY)*modern_rebounded.
c               Use with SEARISEINIT (unless NOLEBROCQ), ZEROICEINIT,
c               ZEROICEBEDEQ, SHALP311[MIN,MAX,AVE]. May need cold start
c               (ifrest=-1), else hb overridden by restart file.
c               (With EISANTA).
c NOISEHB     = Gaussian noise added to bedrock topography
c*WILKESBATHY*= arbitrary lowering of equilibrium topogograhy (topbedeq,
c               equiload) in ~Wilkes-to-Aurora basin sectors. Downward
c               shift of topbedeq = WILKESBATHY * modern |hb|, where
c               hb < 0 (and v.v.,hb>0). (Nb: Doesn't change init h,hb).
c               (With EISANTA).
c MODOBS      = Few-timestep run just to output initial (=modern obs) 
c               state. (With EISANTA).
c
c Only for cold starts (ifrest=-1), 2-D (not EISLINE), not nesting:
c   ZEROICEBEDEQ = reset bedrock equilibrium state (topbedeq, equiload) 
c                  to no-ice rebounded state (default is to set to  
c                  modern observed state, assumed in equilibrium)
c   ZEROICEINIT  = set initial model state to no-ice rebounded state
c   MODELINITREB = calc. of no-ice rebounded state, running bedrock 
c                  module for up to 50 kyr till convergence. Default
c                  is simple local isostatic rebound.

c=================
c Various physics:
c=================

c WATERLAND   = enable open inland lakes or subglacial water,
c               default = only ocean water (old NOWATER).
c MOVEW       = enable advection of liquid water (subr movewater)
c  CHANNEL_SPARSE or CHANNEL_TRIDIA = solver for liq. water flow 
c NOARRHEN    = no dependence of ice rheol coeff on T (subr arrhenius)
c NOBASET     = no dep of basal sliding coeff on basal vars (temp, etc)
c NOMOVET     = no advection or diffusion of temps (icetherm,vdif,hdif)
c NOMOVEI     = no change in ice thickness (icedyn)
c SEDIMENT    = deforming sediment model (not yet avail. with NESTING)
c  *NYEARSEDBUD* = interval for writing sediment budget quantities to 
c                  fort.23 (iusedbud). 
c   SEDTRACK     = sediment tracking for sediment "strata".
c                  Tracks time of deposition vs. x,y,sed_z (in array
c                  sedtrack). sedtrack is used in writehis (iceshow)
c                  to generate to sedstrat array (interpolated elevs
c                  for given times (Ma), written to history file.
c BEDROCK     = bedrock model (locally relaxing asthenosphere)
c  *TAUASTH*    = asthenospheric relaxation-folding time (yr, def=3.e3)
c   IFLITH      = plus elastic lithospheric flexure (if BEDROCK)
c   SLMODEL     = bed model replaced by self-grav. earth-sea level model
c                 (separate program, N.Gomez/J.Mitrovica) (if BEDROCK).
c                 Can be used with EISANTA (2-D) or EISLINE (1-D). 
c                 With EISLINE, also uses SLPROFILE[A,*B*], *SLSIZEB*,
c                 *SLSNOWFALL*, *SLSEALEVEL*,  SLRESTART, *SLNYEAREND*,
c                 SLOLDBED.
c     SLMODELITER = use outer-script iterative running (doslmodeliter)
c                   to correct initial bedrock elevs (hb) for 
c                   previous final difference from modern observed.
c     SLMODELRADIAL = modified SLMODEL with no ocean-sea level gravit.
c                     physics, just eustatic (and no outer iteration). 
c SIMPGEOFLUX = uniform geothermal heat flux values for EAIS and WAIS,
c               with EISANTA (default if not any of following).
c NYBLADEFLUX = Pollard,DeConto,Nyblade (2005) Antarctic geothermal 
c               heat flux map, with EISANTA.
c FOXMAULE    = with EISANTA, SEARISECLIM (unless NOLEBROCQ),  uses 
c               Fox Maule et al.(2005) geothermal heat flux map. 
c SHAPIRO     = like FOXMAULE, except uses Shapiro and Ritzwoller(2004).
c*ENHANCESHEET*= enhancement factor for sheet flow.
c*ENHANCESHELF*= enhancement factor for shelf flow.
c MAXSHELFMELT = very large (~1000 m/y) sub-ice-shelf melt (oceanmelt).
c                Can still have tiny amounts of floating ice.
c OCMARTIN    = sub-ice-shelf oceanic melt rate param. based on Martin 
c               et al.(The Cryo.,2010), and specified ocean temps:
c               LIU2D,RCP[26,45,85],LOVECLIP,CSMCLIM,Levitus(default).
c   Default     = sub-ice oceanic melt uses Levitus94, WOA05 or WOA09 
c                 ocean temps, modern climatol., nearest grid point.
c                 Levitus94/WOA05/WOA09 is hard-coded in iceocean.F.
c   LIU2D       = Liu et al. (2009) ogcm temps (22 to 0 ka, or based on
c                 same d18O from that interval), nearest ogcm grid pt,
c                 nearest Liu time.
c     LIUCORR   = apply bias correction (Levitus-modern) to LIU2D.
c       LIUCORR2= apply only for modern/hot, not lgm climate (weighted)
c   RCP[26,45,85] = RCP scenarios RCP2.5, 4.5, or 8.5, ocean temps 
c                 next 300 yr, from NCAR transient runs, nearest ogcm 
c                 grid pt, interp between decadal time points.
c     RCPCORR   = apply bias correction (Levitus-modern) to RCP.
c               nb: LIU2D and RCP can both be specified (def. Levitus).
c               If both specified, LIU2D used for t<0, RCP for t>0.
c   CSMCLIM     = sub-ice oceanic melt uses CSM ocean temps.
c     CSMANOM   = anomaly method for temp (T:exp+obs-ctl)
c *OCEANDT*     = unif. ocean warming (C) (default=0).
c *DTOCNMATCO2* = unif. ocean warming (C) (default=0), same as OCEANDT,
c                 for use with RAMPMATCO2 above.
c    DTOCNMATCO2ASE = only apply in 'West' (~Amundsen Sea) sector
c *OCFACMULT*   = unif. multiplier for all ocean melt rates (default=1).
c SUBGRIDEDGE = sub-grid interpolation scheme for surface ice budgets
c               on steep marginal slopes of *grounded* ice,
c               with fractional ice cover in last ice grid box, and
c               surface slope interpolation and extra budg calcs
c               for adjacent cell. (In subr climate_all: adjust by 
c               nedge, ifedgenext, hedgenorm).
c               Default is no subgridedge.
c NOSUBGRIDEDGESHELF = don't use subgridedgeshelf, similar sub-grid 
c               scheme for *floating* shelf edges, with fractional ice 
c               in last grid box, and effective ice thickness (hedge). 
c               Default is to use subgridedgeshelf.  
c FACEMELT  = ocean  melt on vertical faces of floating ice, at edges
c             of ice shelves. Needs subgridedgeshelf (def, i.e., 
c             not NOSUGRIDEDGESHELF).
c CLIFFMELT = rapid collapse if vertical exposed (aerial) 
c             calving face is > ~ 100 m. Needs subgridedgeshelf (def).
c  *CLIFFVMAX* = max erosion speed (with CLIFFMELT), def=3.e3.
c  *CLIFFZMAX* = subaerial cliff (m) at start of ramp (w. CLIFFMELT),
c                default = 100.
c  *PELTOVMAX* = max Pelto+Warren erosion speed (with CLIFFMELT),
c                default = 0 (no Pelto+Warren).
c  *PELTOSLOPE*= Pelto+Warren speed vs cliff size (with CLIFFMELT),
c                default = 8.33.
c  *CLIFFMULT* = factor multiplying overall erosion speed (CLIFFMELT),
c                default = 1.
c  *ZUGX*      = new (post AGU13) numerics and physics (see icedyn)
c CALVING   = ice shelf calving parameterization. Includes mechanical
c             and climatic processes: [divergence (div), baseperc],
c             For CALVNICK, crevasse depths=f(div,budgrain+budgmelt).
c             For CALVOLD, crevasse depths=f(div, baseperc).
c             div calcs need subgridedgeshelf (def) for fedge,hedge.
c             Default (no CALVING) is no calving. 
c             If CALVING, default is CALVNICK.
c   *CALVNICK* = based on crevasse depths:ice thickness > a critical
c                value (Nick et al 2010). Nominal crev. depths are 
c                multiplied by factor CALVNICK (default=1).
c     *CALVLIQ* = factor in budgrain+budgmelt influence on liq depth
c                 in sfc crevasses (default = none) (w. NICK)
c       CALVLIQZWHOT = CALVLIQ weighted by zwhot (0 for lgm<->mod).
c     *CALVFRAC* = critical ratio tot.crev.depths:ice thick (w.NICK),
c                  default = 0.75.
c     *CALVVMAX* = max calving erosion speed (w. NICK),
c                  default = 3.e3.
c   Alternate (old) calving, veloc proportional to div:
c   *CALVDIV* = factor (* div) for calving, default = 1.5e5.
c   *CALVPERC* = factor (* baseperc) in calving increase due to 
c                basal liquid percolation, default = 0.
c MELANGE = melange physical component, supplied by calving or cliff.

c==============
c Ice dynamics:
c==============

c SCHOOFGL   = Schoof(2007,JGR) imposed velocities at grounding line 
c   NOBUTTRESS  = no buttressing effect at g.l. (subr thetacalc) 
c  *USCHOOFFAC* = multiplier for Schoof velocs, default=1 (subr calcgl)
c   SCHOOFADJUST= correction in subr schoofgl for sub-grid budgall
c                (not sensible with MAXSHELFMELT, NOSHELF)
c   ANGNORM     = adjust Schoof flux for estimated g.l. orientation
c*SPEEDLIMIT*= cfl-like limit on ice speeds (m/a, subr icedyn).
c              Default (or = 0) is 0.75*dd0/dtimeice. If < 0, no limit. 
c*CRHUCRIT*  = basal coefficent cutoff (* norm. factor for powb):
c              SSA+SIA combo for points with crhu >, SIA with crhu <.
c              Default is 1.e-8 for 3-D, 0 (all SSA+SIA) for flowlines.
c*NITERA*    = max iterations of SSA "A" and "ASCH" loops (default=2)
c*NITERC*    = max iterations of outer "C" loop (default=2)
c*SUBGRIDPINSD*= sub-grid pinning pt area is fn. of hw/(sd(hb)*SUB...)
c                (in fpin, icedyn). Default=2. Set = 0 for none.
c*SUBGRIDPIN*= water-thickness ramp (m), 0->1 pinning (in fpin,icedyn).
c              (Old, overrides SUBGRIDPINSD which is the default).
c SIDEDRAG   = simple param of side drag on floating shelves,
c              proportional to ice veloc and thickness,only with EISLINE
c LINSHELF   = linear SSA rheology (still non-linear SIA)
c NOSHELF    = no floating ice allowed at all. No SSA flow in icedyn, 
c              only grounded ice and SIA (cf. MAXSHELFMELT above).
c NOMOVECOAST = no call to findwater, so coastlines do not move


c NOSOFTCROSS = no strain-softening interactions between SSA and SIA 
c               (uit,vit, sih,siu,siv zero-ed)
c NOLHSCROSS  = no effect of SSA stretching stresses on driving stress 
c               in SIA equations (lhsu,v zero-ed)
c UIACROSS    = include effect of SIA average flow in SSA equations
c               (more rigorous). Default is not, i.e., SSA solution 
c               is for basal flow (uia,via zero-ed).

c DOSPARSE = sparse matrix solver for ice SSA velocities, with one 
c            of the following external solvers:
c   NUMREC,SLAP,SUPERLU,WATSON,MKL,GAUSS
c DOSPARSESOR = Successive OverRelaxation code for ice SSA velocities

c EULER      = simple Euler (upstream) advection of ice temperature
c              (and any other tracers) in icetherm. Default is 
c              parabolic upstream advection.
c WRAPAROUND = E-W (dateline) wraparound in icetherm. Not implemented
c              yet in icedyn(!)...currently LONLAT 360-degree longitude 
c              domains have a no-flow cut at lon boundary.
c NONEWEXTREMA = constrain advected tracers (incl. T) to be within
c                min, max of (old) surrounding pts used in advection.

c=====================================================================
c Inversion runs for basal sliding coefficient, enhanceinv, oceanmelt:
c=====================================================================

c CRHINV     = basal sliding coeffs (crhmel) adjusted at intervals,
c              locally as ad-hoc function of hs-hsbal (in basecoef_inv).
c              Done only for grounded points with ice. 
c              If not NOBASET, don't do for frozen-bed pts, and limit
c              crhmel >= crhhard (unfrozen hard bed).
c              If NOBASET, no restriction on temps, and limit
c              crhmel >= crhnos (~no sliding, subsequently to be 
c              filled by fill_crh).
c              Usually use with NOMOVECOAST,and MAXSHELFMELT or NOSHELF
c              or SHELFINV. hsbal is hs from external datasets read
c              in initphys, i.e., usually modern observed surface elevs.
c  *DTIMECRHINV* = intervals (yrs) to do inverse CRHINV adjustment
c                  to crhmel (overrides default value)
c  *HSCALE*      = scaling for delta[surface elevation] (hs-hsbal) (m)
c                  (overrides default value)
c  *TAUSCALE*    = time interval to gradually ramp to new crhmel after
c                  each adjustment (default = instantaneous adjustment) 
c ENHINV     = enhancement factor for sheet flow (enhanceinv,
c              additional to ENHANCESHEET). enhanceinv is changed,
c              aiming for equal areas of frozen bed and crhmel=crhnos.
c              Adjustments done at same intervals and in same subr
c              (basecoef_inv) as for CRHINV.
c              enhanceinv can be uniform, or depend on distance to
c              nearest ice dome. Generally used in conjunction with
c              CRHINV + NOBASET. 
c CRHMELREAD = at start, read crhmel and enhanceinv from file 
c              (crhmelfilein) written by previous CRHINV/ENHINV run, 
c              overwriting standard initial setting or restart file.
c  *CRHSHELF*  = replace ocean (non-grounded) crhmel with this value
c    *CRHSHELFASEIN* = increment (x) to CRHSHELF for inner Amundsen Sea
c   SAMEDOMAIN = for LONLAT grids, tells program that input-file   
c                domain has same boundaries as model (for older files
c                that don't have grid lon,lats in file headers).
c                Also applies for reading restartin files (IFREST=0,1).
c SHELFINV   = like CRHINV, but for floating ice only, adjusting 
c              oceanmelt locally as ad-hoc function of local h-hbal,
c              every time subr ocean is called (dtimeocn).
c              Usually use with NOMOVECOAST and CRHINV, and no calving.

c===================
c Only with NESTING:
c===================

c IFLAGNEST  = 0 to reset initial vars everywhere from nestdrive file,
c              overriding h,hb settings from input datasets (if 
c              IFREST=-1), or 1 to reset only at domain edges 
c              (default=0).
c              nb: If IFREST=0 or 1, IFLAGNEST is ignored, and initial 
c              h,hb in the interior are set from restart file, 
c              overwriting nestdrive values.
c IFLAGNESTEQ= same for topbedeq, equiload. Also 2 (default)=no reset
c              of input datasets anywhere, neither edges nor interior.
c KEEPARCOCN = use arcocn from nestdrive file(s), not calculated in 
c              oceanarc.
c NESTEVOL   = Use current timeice in subr longterm, allowing
c              long-term forcing to vary (else fix at start time)

c======
c Misc:
c======

c ADAPTDT   = adaptive time stepping (see below: dtimestage, ntry,...)
c SUN       = FPE handling for Suns (subr hand1)
c NETCDF    = write NetCDF format history files(not needed, always done)
c HISTORY3D = include 3-D fields (ice t,u,v) in history files
c OMPLITH   = shared-memory OMP parallelization of lithosphere calcs.
c             Also need: ifort -openmp, setenv OMP_NUM_THREADS 2,3,...
c OMPSTEP   = shared-memory OMP parallelization of SIA dynamics calcs
c             (ADI, subr icestep). Also need: ifort... as above.
c OMPALBEFF = shared-memory OMP parallelization of ALBEFFGREEN calcs 
c             (in calcalbeff). Also need: ifort... as above.
c TABTROUGH = record grounding-line and calving edge locations along  
c             transects (with EISANTA), written to fort.29. Transects
c             defined in subr settransect, written in icetabtrough.
c
c Overrides hard-code settings in icectl.F, and any Namelist settings
c in namelist.in:
c*IFREST*     = -1,0,1 for cold start, initial start, or restart, resp.
c*DTIMEICE*   = primary (dynamical) ice model timestep, years.
c*NYEARSTART* = start of run, years (ignored if IFREST=1).
c*NYEAREND*   = end of run, years.
c*NYEAROUT2D* = interval between 2-D ascii writes (fort.19), years.
c*NYEARHIS*   = interval between history-file writes (fort.92.nc),years.
c*NYEARTAB*   = interval between tabular writes (fort.22), years.
c*NYEARTABBUD*= interval between tabular (budget) writes (fort.28), yrs.
c*NYEARTABMEL*= interval between tab.bud.(melange)writes (fort.280),yrs.
c*NYEARTABTRO*= interval between TABTROUGH writes (fort.29), years.
c*NYEARRES*   = interval between writing restart files, years.
c*NYEARNEST*  = interval between writing nest-driving files, years.

c-----------------------------------------------------------------------

#if defined (RWRAP)
c     for calling as dll module from R-wrapper (R.Fuller, B.Lee, 4/18):
      subroutine rwrap_sheetshelf (rwrap_param, 
     *                             rwrap_plio,
     *                             rwrap_lig,
     *                             rwrap_lgm,
     *                             rwrap_mod,
     *                             rwrap_obs,
     *                             rwrap_2100,
     *                             rwrap_2200,
     *                             rwrap_2300,
     *                             rwrap_2400,
     *                             rwrap_2500,
     *                             rwrap_mod_h, rwrap_mod_mask, 
     *                             rwrap_obs_h, rwrap_obs_mask, 
     *                             rwrap_mod_rms)
#else
      program sheetshelf
#endif

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

c     h grid:
      dimension
     *  h(nx,ny),            hs(nx,ny),          
     *  hb(nx,ny),           hbsd(nx,ny),
     *  t(nx,ny,0:nlevp),    tracer(nx,ny,0:nlevp,ntrace),
     *  arhap(nx,ny),        
     *  s1a0(nx,ny),         s1b0(nx,ny),
     *  s2a0(nx,ny),         s2b0(nx,ny),
     *  heati(nx,ny,nlev),   heath(nx,ny),       heatb(nx,ny),
     *  budgsnow(nx,ny),     budgrain(nx,ny),    
     *  budgevap(nx,ny),     budgmelt(nx,ny),
     *  baseperc(nx,ny),     basefrml(nx,ny),  
     *  oceanmelt(nx,ny),    msector(nx,ny),     
     *  facemelt(nx,ny),     cliffmelt(nx,ny),
     *  calvice(nx,ny),      crevliq(nx,ny),    
     *  div(nx,ny),          visc(nx,ny),
     *  arcocn(nx,ny),       
     *  budgall(nx,ny),      
     *  budgneg(nx,ny),      budgneg2(nx,ny),   ! just for diagnos.
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  maskh(nx,ny),        geoflux(nx,ny),

     *  ifnexo(nx,ny),
     *  bedge(nx,ny),        fedge(nx,ny),       hedge(nx,ny),
     *  tempocn(nx,ny),

     *  sedim(nx,ny),        sedimold(nx,ny),
     *  tsed(nx,ny,nsed),    wsed(nx,ny,nsed),   heats(nx,ny,nsed),
     *  quarryrate(nx,ny),   pelagrate(nx,ny),
     *  quarryacc(nx,ny),    pelagacc(nx,ny),    quarrycoef(nx,ny),   
     *  quarrysec(nsedsec),  quarrysecxy(nx,ny), 
     *  topbed(nx,ny),       topbedeq(nx,ny),
     *  equiload(nx,ny),     deflect(nx,ny),     tbed(nx,ny,nbed),   

     *  hw(nx,ny),           tw(nx,ny),          heatw(nx,ny),
     *  maskwater(nx,ny),    sedpres(nx,ny),     
     *  crhmel(nx,ny),       enhanceinv(nx,ny),
     *  indlake(npoimax,nlakemax), npoilake(nlakemax),

     *  tog(nx,ny),
     *  tmon(nmon,nx,ny),    pmon(nmon,nx,ny),
#if defined (LOVECLIP)
     *  emon(nmon,nx,ny),    qmon(nmon,nx,ny),
#endif
     *  dgeoid(nx,ny)

      equivalence (tracer,t)

#if defined (CO2INTER) 
      dimension
     *  toge(nlong,nlatg),    
     *  tage(nmon,nlong,nlatg), prge(nmon,nlong,nlatg), 
     *  ruge(nlong,nlatg),      runofa(nx,ny)
#elif defined (RACMOCLIM)
      dimension
     *  ruge(nlong,nlatg),      runofa(nx,ny)
#endif

c     for accumulating ocean-related liquid and frozen fluxes,
c     accumulated in subr icetabbud, reset in writehis.
c     (runliqav, runfrozav are for passing to loveclim).

      dimension
     *  oceanmeltav(nx,ny), facemeltav(nx,ny),
     *  cliffmeltav(nx,ny),  calviceav(nx,ny),
     *  runliqav(nx,ny),    runfrozav(nx,ny),
     *  hprev(nx,ny),       maskwaterprev(nx,ny)
      save timefluxprev

c     u and v grids (extra 0 and nxp,nyp pts are just for finite diffs 
c     in icedyn at edges, always zero). (ua,va,uia,via,uit,vit,ub,vb
c     need to be saved between calls to icedyn, and on restart file).
      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),       via(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uit(0:nxp,0:nyp),       vit(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),     maskv(0:nxp,0:nyp),  
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),    fracgv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),     mvind(0:nxp,0:nyp),
     *  thetau(nx,ny),          thetav(nx,ny),
     *  hgu(0:nxp,0:nyp),       hgv(0:nxp,0:nyp),
     *  qgu(0:nxp,0:nyp),       qgv(0:nxp,0:nyp),
     *  uschoof(0:nxp,0:nyp),   vschoof(0:nxp,0:nyp),
     *  fschoofu(0:nxp,0:nyp),  fschoofv(0:nxp,0:nyp),
     *  uagrid(0:nxp,0:nyp),    vagrid(0:nxp,0:nyp),
     *  fluxgrdu(0:nxp,0:nyp),  fluxgrdv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp),  fluxschv(0:nxp,0:nyp),
     *  angnormu(0:nxp,0:nyp),  angnormv(0:nxp,0:nyp),
     *  sih(nx,ny),
     *  siu(0:nxp,0:nyp),       siv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp),
     *  dbu(0:nxp,0:nyp),       dbv(0:nxp,0:nyp),
     *  fsedu(0:nxp,0:nyp),     fsedv(0:nxp,0:nyp),
     *  uw(0:nxp,0:nyp),        vw(0:nxp,0:nyp)

      dimension
     *  ubal(0:nxp,0:nyp),      vbal(0:nxp,0:nyp),
     *  hsbal(nx,ny),           hbal(nx,ny),
     *  maskinit(nx,ny),        hbinit(nx,ny)

#if defined (LOVECLIP)
c for climfile[in,out,ctl] (in common), to be in namelist:
#include <comiceloveclim.h>
#endif

      dimension hislist(1000)
      namelist /runparams/ 
     *  dtimeice, dtimebed, ifrest, 
     *  nyearstart, nyearend, nyearres, nyearnest, 
     *  nyearout1d, nyearout2d, nyearhis, nyeartab, nyearplot1d, 
     *  nyeartabwais, nyeartabbud, nyeartabtro,
     *  hislist, cnestpath, cnestlist,
     *  itersl
#if defined (LINEM)
     * ,amismip, xgexact
#endif
#if defined (LOVECLIP)
     * ,climfilein, climfileout, climfilectl   ! in common, comicepath.h
#endif
#if defined (MELANGE)
     * ,nyeartabmel
#endif

      logical ifexist, ifblow
      character*80 cresin

#if defined (LINEE) && defined (RESTARTE)
#include <fort.26x>
#endif

#if defined (SUN)
      external hand1
c     external hand1, hand2, hand3
#endif

#if defined (SLMODEL)
      character*24 cform
#  if defined (SLMODELITER) || defined (SLMODELRADIAL)
      character*80 cfile
      dimension hbstart(nx,ny), dhb(nx,ny)
#  endif
#endif

#if defined (SEDIMENT) && defined (FORCE34MYR) && (defined (QUARRYREAD) || defined(QUARRYWRITE))
      character*80 cfile
      dimension hmod(nx,ny), hbmod(nx,ny), hbrebmod(nx,ny),
     *          hbinitprev(nx,ny),               ! hbinit declared above
     *          hbrebprev(nx,ny), hbreb(nx,ny),
     *          work(nx,ny),  work2(nx,ny),
     *          quarrysecold(nsedsec),quarrysecxyold(nx,ny)
#endif

#if defined (NOMOVEI)
      dimension hsav(nx,ny)
#endif

#if defined (MELANGE)
      dimension 
     *  hmel(nx,ny),           hsmel(nx,ny),        hwmel(nx,ny),
     *  fmel(nx,ny),           maskmel(nx,ny),
     *  budgallmel(nx,ny),     budgsupmel(nx,ny),   
     *  budgnegmel(nx,ny),     budgadjmel(nx,ny),   
     *  budgredmel(nx,ny),     budgsmomel(nx,ny),
     *  oceanmeltmel(nx,ny),
     *  uamel(0:nxp,0:nyp),    vamel(0:nxp,0:nyp),
     *  ucimel(0:nxp,0:nyp),   vcimel(0:nxp,0:nyp),
     *  hucimel(0:nxp,0:nyp),  hvcimel(0:nxp,0:nyp),
     *  bpumel(0:nxp,0:nyp),   bpvmel(0:nxp,0:nyp),
     *  thetaumel(nx,ny),      thetavmel(nx,ny),
     *  uaboundmel(0:nxp,0:nyp), vaboundmel(0:nxp,0:nyp)
#endif

#if defined (RWRAP)
c     passed arguments for calling from R-wrapper:
      parameter (nrwrap_param = 11)
      dimension rwrap_param(nrwrap_param),         ! input
     *          rwrap_plio(3),                     ! output
     *          rwrap_lig(3),
     *          rwrap_lgm(3),
     *          rwrap_mod(3),
     *          rwrap_obs(3),
     *          rwrap_2100(3),
     *          rwrap_2200(3),
     *          rwrap_2300(3),
     *          rwrap_2400(3),
     *          rwrap_2500(3),
     *          rwrap_mod_h(nx,ny), rwrap_mod_mask(nx,ny),
     *          rwrap_obs_h(nx,ny), rwrap_obs_mask(nx,ny)

      logical firstrwrap
      data firstrwrap /.true./
      save firstrwrap
#endif

c-------------------------
c Start of executable code
c-------------------------

c        Set floating-point exception handler for sun (solaris) 

#if defined (SUN)
      ier = ieee_handler ('set', 'division', hand1)
      ier = ieee_handler ('set', 'overflow', hand1)
      ier = ieee_handler ('set', 'invalid',  hand1)
c     ier = ieee_handler ('set', 'overflow', hand2)
c     ier = ieee_handler ('set', 'invalid',  hand3)
#endif

c        Delete stub file that will signify normal end (not abend)

      open (iuokend, file='okend', status='unknown') 
      close (iuokend, status='delete') 

#if defined (RWRAP)
c     set input variables from R-wrapper argumnent (rwrap_param),
c     to variables in common (comicegrid.h), for use throughout code:
      rwrap_ocfacmult    = rwrap_param(1)
      rwrap_ocfacmultase = rwrap_param(2)
      rwrap_calvnick     = rwrap_param(3)
      rwrap_crhshelf     = rwrap_param(4)
      rwrap_tauasth      = rwrap_param(5)
      rwrap_calvliq      = rwrap_param(6)
      rwrap_cliffvmax    = rwrap_param(7)
      rwrap_facemeltrate = rwrap_param(8)
      rwrap_enhancesheet = rwrap_param(9)
      rwrap_enhanceshelf = rwrap_param(10)
      rwrap_crhfac       = rwrap_param(11)

c     initialize output variables to null values, for R-wrapper:
c     or - leave alone, so as not to-reset previous R-wrapper results.
c     do m=1,3
c       rwrap_plio(m) = 0.
c       rwrap_lig(m) = 0.
c       rwrap_lgm(m) = 0.
c       rwrap_mod(m) = 0.
c       rwrap_obs(m) = 0.
c       rwrap_2100(m) = 0.
c       rwrap_2200(m) = 0.
c       rwrap_2300(m) = 0.
c       rwrap_2400(m) = 0.
c       rwrap_2500(m) = 0.
c     enddo
c     call zero (rwrap_mod_h,    nx*ny)
c     call zero (rwrap_mod_mask, nx*ny)
c     call zero (rwrap_obs_h,    nx*ny)
c     call zero (rwrap_obs_mask, nx*ny)
c     rwrap_mod_rms = 0.
#endif

c        Set run parameters (which are defaults for namelist variables -
c        namelist read later below) 

      cnestpath = ' '
      call resetr (hislist, 1000, -999999.)
      do m=1,1000
        cnestlist(m) = ' '
      enddo 
      itersl = 0   ! for no outer iterations (no SLMODELITER)

#if defined (SLMODEL)
      dtimebed = 200.
#elif defined (NHA) && defined (LOVECLIP)
      dtimebed = 50.
#else
      dtimebed = 1.                   ! changed 4/15 (from 50)
#endif
#if defined (DTIMEBED)
      dtimebed = DTIMEBED             ! deconto:  = 1 to eliminate jumps
#endif

c--------------------
#if defined (EISLINE)
c--------------------

      if (dd0.le.0.2e3) then
        dtimeice = 0.1
      else 
        dtimeice = 0.5
      endif

      ifrest = -1
c     ifrest =  0
c     ifrest =  1
#  if defined (LINEE) && (defined (RESTARTE) || defined (SEALEVCHANGE))
      ifrest =  0
#  endif

      nyearstart  =  0
      nyearout2d  =  0
      nyearnest   =  0
      nyearhis    =  0
      nyearres    =  0
      nyeartab    = 50

#  if defined (LINEC)
#    if defined (SINBED) 
      nyearend    = 250000
      nyearout1d  =  10000
      nyearplot1d = 250000
#    elif defined (FLATANT)
      nyearend    = 1000000
      nyearout1d  =   10000
      nyearplot1d =  250000
#    else
      nyearend    = 750000
      nyearout1d  =  10000
      nyearplot1d = 250000
#    endif
#  elif defined (LINEE)
#    if defined (SEALEVCHANGE)
      nyearend    = 300
      nyearout1d  = 10
      nyearplot1d = 10
#    elif defined (INITEE) || defined (RESTARTE) 
      nyearend    = 5000 
      nyearout1d  = 50
      nyearplot1d = 10
#    else
      nyearend    = 10000
      nyearout1d  = 50
      nyearplot1d = 10
#    endif
#  elif defined (LINEF)
      nyearstart = 0
      nyearend    = 100000
      dtimeice    =     .5
c     dtimeice    =     1.
      nyearout1d  =    100
c     nyearout1d  =   1000
      nyearplot1d =   1000
c     nyearres    =    500
      nyearres    = 100000
#  elif defined (LINEH) || defined (LINEH2)
      nyearstart  = -37000000
      nyearend    = -25000000
      nyearout1d  =    100000
      nyearplot1d =    100000
      nyearres    = nyearend - nyearstart
      nyeartab    =      1000
#  elif defined (LINEH3) || defined (LINEH4)
      nyearstart  =         0
      nyearend    =    100000
      nyearout1d  =       100
      nyearplot1d =       100
      nyearres    =     10000
      nyeartab    =         1
#    if defined (LINEH4)
      nyearout2d  =       100
#    endif
#  elif defined (LINEM)
      nyearstart =     0
c     nyearend   = 50000
      nyearend   = 70000
c     nyearend   = 100000
#    if defined (RESOL1)
      dtimeice   = 1.
#    elif defined (RESOL2)
      dtimeice   = 1.
#    elif defined (RESOL3)
      dtimeice   = 1.
#    elif defined (RESOL4)
      dtimeice   = 1.
#    elif defined (RESOL5)
      dtimeice   = 0.2
#    elif defined (RESOL6)
      dtimeice   = 0.01
#    endif
      nyearout1d =   1000
      nyearplot1d =   0
      nyearres   = nyearend 
      nyeartab   =  50
      ifrest     = -1
#  elif defined (SLMODEL)
      nyearstart =       0
      nyearend    = 200000
      nyearout1d  =    100
      nyearplot1d =    100
      nyeartab    =     50
      nyearres    =  50000
#    if defined (SLRESTART)
      ifrest = 0
      nyearend    = 100000
#    else
      ifrest = -1
#    endif
#    if defined (SLNYEAREND)
      nyearend = SLNYEAREND
#    endif
#  else
      nyearstart = 0
      nyearend    = 200000
      nyearout1d  =   1000
      nyearplot1d =   1000
#  endif

c---------------------
#elif defined (ASYNCH)
c---------------------

#  if defined (NOSHELF)
      dtimeice    = 1.
#  else
      dtimeice    = 0.1 
#  endif

c     nyearstart,nyearend,ifrest set in namelist
      nyearout1d  =   1000
      nyearplot1d =    nyearout1d
      nyearout2d  =    nyearout1d
      nyearres    =  20000
      nyearnest   =      0
      nyearhis    =   1000
      nyeartab    =     50
      nyeartabwais=      0
      nyeartabbud =      0
      nyeartabtro =      0

c-------------------------------------------
#elif defined (EISANTA) && defined (NESTING)
c-------------------------------------------

      if (nint(dd0*1.e-3).le.1) then
        dtimeice = 0.02
      else if (nint(dd0*1.e-3).le.2) then
        dtimeice = 0.05
      else if (nint(dd0*1.e-3).le.5) then
        dtimeice = 0.1
      else
        dtimeice = 0.5
      endif

      ifrest =-1
c     ifrest = 0
c     ifrest = 1
 
      nyearstart  = -20000
      nyearend    = 0
c     nyearstart = 0
c     nyearend    = 2
c     nyearend    = nyearstart + 20000
      nyearout1d  =  100
c     nyearout1d  =  10
      nyearplot1d =  nyearout1d
      nyearout2d  =  nyearout1d
      nyearres    =  500 
      nyearnest   =  0
      nyearhis    =  nyearout1d 
      nyeartab    =  2
      nyeartabbud = nyearout1d
      nyeartabtro = nyearout1d

c---------------------------------------------
#elif defined (EISANTA) && ! defined (NESTING)
c---------------------------------------------

# if defined (NOSHELF)
      if (nint(dd0*1.e-3).le.10) then
        dtimeice = 0.25
      else if (nint(dd0*1.e-3).le.20) then
        dtimeice = 0.50
      else
        dtimeice = 1.
      endif
# else
      if (nint(dd0*1.e-3).le.10) then
        dtimeice = 0.125
      else if (nint(dd0*1.e-3).le.20) then
        dtimeice = 0.25
      else
        dtimeice = 0.5
      endif
# endif

#  if defined (TRANSECTA) 
      ifrest = -1
#  else
c     ifrest = -1
      ifrest =  0
c     ifrest =  1
#  endif

c.....................
#  if defined (MODOBS)
c.....................

      ifrest = -1
      nyearstart  = 0
      nyearend    = 2
      nyearout1d  = 2
      nyearplot1d = 2 
      nyearout2d  = 2
      nyearres    = 0
      nyearnest   = 0
      nyearhis    = 2
      nyeartab    = 1
      nyeartabwais= 0
      nyeartabbud = 0
      nyeartabtro = 0

c...................................................
#  elif defined (FORCEPLEIST) && ! defined (SLMODEL)
c...................................................

c     nyearstart  = -50000000
c     nyearend    = -40000000
c     nyearstart  = 0
c     nyearend    = 1000000
c     nyearstart  =  50000
c     nyearend    =  65000
c     nyearstart  =  65000
c     nyearend    =  90000
c     nyearstart  =      0
c     nyearend    =  70000
c mis31:
c     nyearstart  =  -1120000
c     nyearend    =  -1020000
c last glac cycle:
c     nyearstart  =  -160000
c     nyearstart  =   -80000
c     nyearstart  =   -40000
c     nyearstart  =   -20000
c     nyearend    =        0
c     nyearend    =     5000
c want:
c     nyearstart  =        0
c     nyearend    =   300000
c mod:
      nyearstart  =        0
      nyearend    =   100000
c     nyearstart  =   -40000
c     nyearend    =        0

c     nyearout1d  = 5000
      nyearout1d  = 1000
c     nyearout1d  =  500

c     nyearplot1d = 10000
      nyearplot1d =     0

c     nyearout2d  = 10000
c     nyearout2d  =  5000
      nyearout2d  =  1000
c     nyearout2d  =   100
c     nyearout2d  =    10

c     nyearres    =100000
c     nyearres    = 50000
c     nyearres    = 20000
      nyearres    = 10000

c     nyearnest   = 50000
c     nyearnest   = 10000
c     nyearnest   = 5000
c     nyearnest   = 1000
c     nyearnest   =  500
      nyearnest   = 0

      nyearhis    = nyearout2d
c     nyearhis    =1000
c     nyearhis    = 500
c     nyeartab    =  50
c     nyeartab    =  10
      nyeartab    =   1
c     nyeartabwais= nyeartab
      nyeartabwais= 0
      nyeartabbud = 0
      nyeartabtro = 0

#    if defined (RCMMATRIX)
      ifrest      =      0
      nyearstart  =     0
      nyearend    = 10000
      nyearout2d  =   100
      nyearhis    = nyearout2d
#    elif defined (CRHINV) || defined (ENHINV)
      ifrest     =     -1
      nyearstart =      0
      nyearend   = 400000
      nyearout2d =  10000
      nyearhis    = nyearout2d
#    endif

c.................................................
#  elif defined (FORCEPLEIST) && defined (SLMODEL)
c.................................................

c     nyearstart  =        0
c     nyearend    =    15000
c     nyearstart  =   -20000
      nyearstart  =   -40000
      nyearend    =        0

c     ifrest      = -1
      ifrest      = 0

      nyearout1d  =   50
      nyearplot1d =    0
      nyearout2d  =  200
      nyearres    = 5000
      nyearnest   =    0
      nyearhis    = nyearout2d
      nyeartab    =   10 
      nyeartabwais=    0
      nyeartabbud =    0
      nyeartabtro =    0

c........................
#  elif defined (FORCEEO)
c........................

      ifrest = -1

c     hyst:
c     nyearstart  = -37000000
c     nyearend    = -25000000
c     nyearout1d  =    100000
c     nyearout2d  =  nyearout1d
c     nyearplot1d =  nyearout1d
c     nyearres    =   1000000
c     nyearnest   =         0
c     nyearhis    =    100000
c     nyeartab    =      1000
c     nyeartabwais=         0
c     nyeartabbud =         0
c     nyeartabtro =         0

c     co2:
#    if defined (CO2INTER)
      nyearstart  = -34000000
      nyearend    = -32000000
#    else
      nyearstart  = -36000000
      nyearend    = -25000000
#    endif
      nyearout1d  =     50000
      nyearout2d  =     50000
      nyearplot1d =     50000
      nyearres    =   1000000
      nyearnest   =         0
      nyearhis    =     50000
      nyeartab    =      1000
      nyeartabwais=         0
      nyeartabbud =         0
      nyeartabtro =         0

c......
#  else
c......

c     not MODOBS, FORCEPLEIST or FORCEEO:
      nyearstart  =        0
      nyearend    =    10000
      ifrest      = -1
      nyearout1d  =   50
      nyearplot1d =    0
      nyearout2d  =  200
      nyearres    = 10000
      nyearnest   =    0
      nyearhis    = nyearout2d
      nyeartab    =   10 
      nyeartabwais=    0
      nyeartabbud =    0
      nyeartabtro =    0

c.......
#  endif
c.......

c------------------
#elif defined (NHA)
c------------------

#  if defined (NOSHELF)

      dtimeice = 1.

#  elif defined (LONLAT)

      if (nint(dlatd*100.).le.5) then
         dtimeice = 0.025
      else if (nint(dlatd*100.).le.10) then
         dtimeice = 0.05
      else if (nint(dlatd*100.).le.20) then
         dtimeice = 0.1
      else if (nint(dlatd*100.).le.50) then
         dtimeice = 0.2
      else
c        dtimeice = 0.5
         dtimeice = 1.0
      endif 

#  endif

      ifrest = -1
c     ifrest = 0 

      nyearstart  = 0
      nyearend    = 200000 
      nyearout1d  = 500
      nyearplot1d = 0
      nyearout2d  = 2000
      nyearres    = 20000
      nyearnest   = 0
      nyearhis    = nyearout2d
      nyeartab    = 10

c-------------------------------------
#elif defined (CARB) || defined (ORDO)
c-------------------------------------

c     dtimeice = 10.
      dtimeice = 5.
c     dtimeice = 2.
c     dtimeice = 1.

      ifrest = -1

      nyearstart  = 0
      nyearend    = 50000 
      nyearout1d  = 100
      nyearplot1d = 0
      nyearout2d  = nyearout1d
      nyearres    = 5000
      nyearnest   = 0
      nyearhis    = 5000
      nyeartab    = 50

c---------------------
#elif defined (TEST2D)
c---------------------

      zdx = min (dx0, dy0)
      if (nint(zdx*1.e-3).le.0.5) then
        dtimeice = 0.1
      else if (nint(zdx*1.e-3).le.1) then
        dtimeice = 0.25
      else if (nint(zdx*1.e-3).le.20) then
        dtimeice = 0.5
      else
        dtimeice = 1.
      endif

#if defined (MISMIPPLUS_spin)
      ifrest     = -1
      nyearstart = 0
      nyearend   = 20000 
      nyearhis    = 100
      nyearout1d  = 100
      nyearout2d  = 100
      nyeartab    = 10
      nyeartabbud = 10
      nyearres    = 5000
#else
c     MISMIPPLUS_[ctl,1rr,1ra,2rr,r2a]:
      ifrest     = 0
      nyearstart = 0
      nyearend   = 1000 
      nyearhis    = 10
      nyearout1d  = 10
      nyearout2d  = 10
      nyeartab    = -1
      nyeartabbud = -1
      nyearres    = 1000
#endif
      nyearplot1d = 0
      nyearnest   = 0

c~~~~~~~~~~~~~~~~~~~~~~~~
#elif defined (TESTMEL2D)
c~~~~~~~~~~~~~~~~~~~~~~~~

      zdx = min (dx0, dy0)
      if (nint(zdx*1.e-3).le.10) then
        dtimeice = 0.25
      else if (nint(zdx*1.e-3).le.20) then
        dtimeice = 0.5
      else if (nint(zdx*1.e-3).le.40) then
        dtimeice = 1.0
      else
        dtimeice = 0.5
      endif
      ifrest     = -1
      nyearstart = 0
      nyearend   = 20000 
      nyearhis    = 100
      nyearout1d  = 100
      nyearout2d  = 100
      nyeartab    = 10
      nyeartabbud = 10
      nyearres    = 20000
      nyearplot1d = 0
      nyearnest   = 0

c----
#else
c----
      write (ioterm,'(a)') 
     *  '*** Error (sheetshelf): no known run-type defined in makefile'
      stop
c-----
#endif
c-----

      nyearlag = 0
#if defined (MELANGE)
      nyeartabmel = 0
#endif

c        Read namelist (if file exists), overwriting values set above 

      inquire (file='namelist.in', exist=ifexist)
      if (ifexist) then
        write (ioterm,'(a)') 'Reading namelist file'
        open (iunamel, file='namelist.in', status='old')
        read (iunamel, runparams)
        close (iunamel)
      endif

      nhislist = 0
      do m=1000,1,-1
        if (hislist(m).ne.-999999.) then
          nhislist = m
          goto 100  
        endif
      enddo
  100 continue

c     nestdrive filename list, for nesting only (in comicegrid):
      nnestlist = 0
      do m=1000,1,-1
        if (cnestlist(m).ne.' ') then
          nnestlist = m
          goto 120  
        endif
      enddo
  120 continue
      if (nnestlist.eq.0) then
        nnestlist = 1
        cnestlist(1) = 'nestdrivein'
      endif

      if (cnestpath.ne.' ') then
        n = lenchr(cnestpath)
        if (cnestpath(n:n).ne.'/') then
          cnestpath(n+1:n+1) = '/'
        endif
      endif

c        Override variables set by -D in makefile (even if in namelist)

#if defined (IFREST)
      ifrest = IFREST
#endif
      if (ifrest.eq.2) ifrest = -1     ! in case can't have neg namelist

#if defined (DTIMEICE)
      dtimeice = DTIMEICE
#endif

#if defined (NYEARSTART)
      nyearstart = NYEARSTART
#endif

#if defined (NYEAREND)
      nyearend = NYEAREND
#endif

#if defined (NYEAROUT2D)
      nyearout2d = NYEAROUT2D
#endif

#if defined (NYEAROUT1D)
      nyearout1d = NYEAROUT1D
#endif

#if defined (NYEARPLOT1D)
      nyearplot1d = NYEARPLOT1D
#endif

#if defined (NYEARTAB)
      nyeartab = NYEARTAB
#endif

#if defined (NYEARHIS)
      nyearhis = NYEARHIS
#endif

#if defined (NYEARRES)
      nyearres = NYEARRES
#endif

#if defined (NYEARNEST)
      nyearnest = NYEARNEST
#endif

#if defined (NYEARTABBUD)
      nyeartabbud = NYEARTABBUD
#endif

#if defined (NYEARTABTRO)
      nyeartabtro = NYEARTABTRO
#endif

#if defined (NYEARLAG)
      nyearlag = NYEARLAG
#endif

#if defined (LINEM) && defined (AMISMIP)
      amismip = AMISMIP
#endif

#if defined (MELANGE) && defined (NYEARTABMEL)
      nyeartabmel = NYEARTABMEL
#endif

c        Set remaining run variables 

#if defined (GCMMATRIX)
#  if defined (CO2INTER)
      dtimeclim = 1000.
#  else
      dtimeclim = 100.
#  endif
#elif defined (RCMMATRIX)
      dtimeclim = 20.
#elif defined (EBM)
      dtimeclim = 1000.
#elif defined (NHA) && defined (LOVECLIP) || defined (CARB) 
      dtimeclim = 50.
#elif defined (CSMORDO)
      dtimeclim = 50.
#elif defined (CSMCLIM) || defined (CLEMCLIM)
      dtimeclim = 1.                  ! J. Tsai (call once a year)
#else
      dtimeclim = 1.                  ! changed 4/15 (from 50)
#endif
#if defined (DTIMECLIM)
      dtimeclim = DTIMECLIM           ! deconto:  = 1 to eliminate jumps
#endif

      dtimebudg = 100.
#if defined (DTIMEBUDG)
      dtimebudg = DTIMEBUDG
#endif

#if defined (NESTING) || defined (TEST2D) || defined (TESTMEL2D)
      dtimetherm = dtimeice
#elif defined (NHA) && defined (LOVECLIP)
      dtimetherm = 50.
#else
      dtimetherm = 1.                 ! changed 4/15 (from 50)
#endif
#if defined (DTIMETHERM)
      dtimetherm = DTIMETHERM         ! deconto:  = 1 to eliminate jumps
#endif

      dtimehyd = dtimeice
      dtimesed = dtimeice
      dtimeocn = dtimeice    ! 888888
c     dtimeocn = 10.         ! 888888
#if defined (MELANGE)
      dtimemel = dtimeice
#  if defined (DTIMEMEL)
      dtimemel = DTIMEMEL
#  endif
#endif

      dtimegetnest = dtimeice
      nyearsedbud  = nyeartab

#if defined (CO2INTER)
      dtimeco2 = 200.
#endif

      nyeartot = nyearend - nyearstart 

c        Initialize horizontal and vertical grids, 
c        and physical constants in commons

#if defined (EISANTA) && defined (TRANSECTA) && defined (ZOOMGL)
      iczoom = nx0                  ! initialize iczoom to rhs of domain
#else
      iczoom = 0
#endif

      call initgrid (iczoom)

      timeice = nyearstart

c        Top of main timestep loop. If ADAPTDT, divide into outer 
c        "stages" of length dtimestage (yrs), and inner loops through 
c        each stage with timestep dtimeice. Then, if a blowup is 
c        detected at the end of any inner timestep (ifblow), the 
c        current stage is restarted with halved timesteps (intermediate
c        "itry" loop), reading in the restart file (restart, not
c        restartin) written at the end of the previous successful stage.
c        (If not ADAPTDT, dtimestage = length of run, nstage=1, ntry=1).
c
c        Use "-fpe1" compile option to continue execution if fpe occurs,
c        Blowups are detected by extreme-value tests below 
c        (setting ifblow) at end of each timestep.
c
c        nb: If blowup+restart, can have spurious output for SLMODEL.

#if defined (ADAPTDT)
#  if defined (DTIMESTAGE)
      dtimestage = DTIMESTAGE
#  elif defined (NESTING)
      dtimestage = min (100., nyearend-timeice)
#  elif defined (FORCEEO)
      dtimestage = min (100000., nyearend-timeice)
#  elif defined (FORCE34MYR)
      dtimestage = min (200000., nyearend-timeice)
#  else
      dtimestage = min (1000.,nyearend-timeice)
#  endif
      nstage = nint((nyearend-timeice)/dtimestage)
c     ntry = 4
      ntry = 5  ! first retry (itry=2) with unaltered timestep
#else
      dtimestage = nyearend-timeice
      nstage = 1
      ntry = 1
#endif

c        Save basic model physical timesteps, possibly reduced
c        within each stage if blowups occur.
c        nb: for timing control below, these should all divide
c        exactly into dtimestage.

      dtimeicesav   = dtimeice
      dtimethermsav = dtimetherm
      dtimehydsav   = dtimehyd
      dtimebedsav   = dtimebed
      dtimesedsav   = dtimesed
      dtimeocnsav   = dtimeocn

      ifblow = .false.

c        Overall outer "stage" loop

c>>>>>>>>>>>>>>>>>>>>>>>
      do istage=1,nstage
c>>>>>>>>>>>>>>>>>>>>>>>

c        Adaptive time-step attempts to get through current stage

c>>>>>>>>>>>>>>>>>>>
      do itry=1,ntry
c>>>>>>>>>>>>>>>>>>>

c     dtimeice   = dtimeicesav   / (2**(itry-1))
c     first retry (itry=2) with unaltered timestep:
      dtimeice   = dtimeicesav   / (2**max(0,(itry-2)))
      dtimetherm = dtimethermsav / (2**max(0,(itry-2)))
      dtimehyd   = dtimehydsav   / (2**max(0,(itry-2)))
      dtimebed   = dtimebedsav   / (2**max(0,(itry-2)))
      dtimesed   = dtimesedsav   / (2**max(0,(itry-2)))
      dtimeocn   = dtimeocnsav   / (2**max(0,(itry-2)))

c     just for writemark message below:
      dtimeicenext = dtimeicesav   / (2**max(0,((itry+1)-2)))

      nloop = nint(dtimestage/dtimeice)

c        Perform full initialization (initphys, sets to zero,
c        call readres if ifrest!=0), first time through only,
c        or at beginning of each itry > 1 (i.e., ifblow) for ADAPTDT

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      if ( (istage.eq.1 .and. itry.eq.1) .or. ifblow )  then
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

c        Initialize physical variables

      call initphys (h, hs, hb, hbsd, t, tracer, 
     *               sedim,  tsed, wsed, 
     *               topbed, topbedeq, deflect, 
     *               tbed, equiload, geoflux,
     *               hw, tw, maskwater, maskinit, hbinit,
#if defined (SEDIMENT) && defined (FORCE34MYR) && (defined (QUARRYREAD) || defined(QUARRYWRITE))
     *               hbinitprev, hbrebprev, hbreb,
     *               quarrysec, quarrysecxy,
#endif
     *               sedpres, arcocn, 
     *               hsbal, hbal, indlake, npoilake, nlake,
     *               ua, va, sealev,
     *               float(nyearstart), ifrest)

c     used in icetherm before first set in icedyn:
c     (should be on restart file, as should uia,via):
      call zero (dfu,(nxp+1)*(nyp+1))    ! set in icedyn
      call zero (dfv,(nxp+1)*(nyp+1))    ! set in icedyn
      call zero (sih, nx*ny)             ! set in icedyn
      call zero (siu, (nxp+1)*(nyp+1))   ! set in icedyn
      call zero (siv, (nxp+1)*(nyp+1))   ! set in icedyn
      call zero (heath,nx*ny)            ! set in icedyn
      call zero (heatb,nx*ny)            ! set in icedyn, sedflow
      call zero (heats,nx*ny*nsed)       ! set in sedflow
      call zero (dbu,(nxp+1)*(nyp+1))    ! set in icedyn
      call zero (dbv,(nxp+1)*(nyp+1))    ! set in icedyn
c     used in iceocean before first set in icedyn:
      call resetr (visc, nx*ny, 1.)      ! set in icedyn
c     diagnostic in icedyn:
      call zero (ub, (nxp+1)*(nyp+1))    ! set in icedyn
      call zero (vb, (nxp+1)*(nyp+1))    ! set in icedyn
c     for safety, used (budgall,icedyn) before first set in icetherm
c     (vdif), unless icetherm called first step.Should be on rstrt file:
      call zero (baseperc,nx*ny)         ! set in vdif
      call zero (basefrml,nx*ny)         ! set in vdif

c     just for iceshow...:
      call zero (hedge,nx*ny)
      call zero (quarryrate,nx*ny)
      call zero (pelagrate,nx*ny)
      call zero (heati,nx*ny*nlev)
      call zero (fsedu,(nxp+1)*(nyp+1))
      call zero (fsedv,(nxp+1)*(nyp+1))
      call zero (u,(nxp+1)*(nyp+1)*(nlevp+1))
      call zero (v,(nxp+1)*(nyp+1)*(nlevp+1))
      call zero (ubal,(nxp+1)*(nyp+1))
      call zero (vbal,(nxp+1)*(nyp+1))

#if defined (MELANGE)
      call zero (hmel,nx*ny)
      call resetr (thetaumel, nx*ny, 1.)
      call resetr (thetavmel, nx*ny, 1.)
#endif

      call zero (quarryacc,nx*ny)             ! read  from  restart file
      call zero (pelagacc,nx*ny)

      rco2inter = 4. ! for interactive co2 (CO2INTER)

c        If initial start or restart (ifrest=0 or 1), read restart file
c        (for prognostic variables, and needed "inter-timestep" 
c        diagnostic variables)

c~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (ifrest.ge.0) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (ifblow) then
          cresin = 'restartadapt'
        else
          cresin = 'restartin'
        endif
        call readres  (h, hs, hb, t, tracer, 
     *                 baseperc, basefrml, heatb, 
     *                 sedim, tsed, wsed, heats, quarryacc, hbinit,
     *                 topbed, topbedeq, equiload, tbed, hw, tw, 
     *                 sedpres, crhmel, enhanceinv,
     *                 ua, va, ub, vb, dfu, dfv, dbu, dbv,
#if defined (MELANGE)
     *                 hmel,
#endif
     *                 sealev, rco2inter, timeicein, iczoom,
     *                 cresin)

        call printmap (timeice, sedpres, 'sedpres', 0.1, iuout2d, 0)

        if (ifrest.eq.0) rco2inter = 4.

#if defined (LINEE) && defined (RESTARTE)
        xglwedge1 = 1.e20
        jxglwedge2 = 1.e20
        do j=1,ny
          do i=1,nx
            if (hw_wedge(i).gt.0.) then
              if (xglwedge1.eq. 1.e20) xglwedge1 = xh(i,j)-0.5*dx(i,j)
              if (xh(i,j)-xglwedge1 .gt. 0. .and.
#if defined (WEDGE4KM)
     *            xh(i,j)-xglwedge1 .lt. 4.e3) then     !  4 km
#elif defined (WEDGE8KM)
     *            xh(i,j)-xglwedge1 .lt. 8.e3) then     !  8 km
#elif defined (WEDGE12KM)
     *            xh(i,j)-xglwedge1 .lt. 12.e3) then    ! 12 km
#elif defined (WEDGE16KM)
     *            xh(i,j)-xglwedge1 .lt. 16.e3) then    ! 16 km
#endif
                topbedeq(i,j) = hswedge(i) - hwedge(i) + .01
                topbed(i,j) = topbedeq(i,j)
                hb(i,j) = topbedeq(i,j)
                hw(i,j) = 0.
              else
                if (xglwedge2.eq. 1.e20) then
                  xglwedge2 = xh(i-1,j)
                  zglwedge2 = topbedeq(i-1,j)
                endif
                topbedeq(i,j) = max (topbedeq(i,j),
     *                               zglwedge2 -.01*(xh(i,j)-xglwedge2))
                topbed(i,j) = topbedeq(i,j)
                hb(i,j) = topbedeq(i,j)
                hw(i,j) = hs(i,j) - h(i,j) - topbed(i,j)
              endif
            endif
          enddo
        enddo
#endif

        if (ifrest.eq.1) then
          nyearendin = nint(timeicein)
          if (nyearend.le.nyearendin) then
            write (6,'(a/a,i8,a,i8)')
     *      '*** Error: end date is earlier than restart-file end date',
     *      '    nyearend=',nyearend, '  nyearendin=',nyearendin
            stop
          endif
          timeice = nyearendin
        endif

c          Reset variables around nested domain edges

#if defined (NESTING)
#  if defined (IFLAGNESTEQ)
        iflagnesteq = max (IFLAGNESTEQ, 1)
#  else
        iflagnesteq = 2
#  endif
        call getnest (h, hb, t, tbed, sedim, tsed, wsed,
     *                topbedeq, equiload, sedpres, arcocn, maskwater,
     *                ua, va, 1, iflagnesteq, 1, timeice)
#endif

c          Set initial maskwater, indlake, etc (not on restart file),
c          and adjust hw and hs (adjustpres)

#if ! defined (NOMOVECOAST)
        call findwater (maskwater, indlake, npoilake, nlake,
     *                  h, hb, hw, sealev, timeice)
#endif

        call adjustpres (maskwater, indlake, npoilake, nlake,
     *                   h, hb, hw, hs, sealev)

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      else if (ifrest.eq.-1) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c          If cold start, set basal sliding coefficients (crhmel, on 
c          h-grid, ignoring ice melting and grounded criteria), based 
c          on sedpres (set in initphys only for cold start). If restart,
c          crhmel (and enhanceinv, sedpres) have been read above from
c          restart file.

        call basecoef_init (sedpres, maskinit, crhmel, timeice) 
        call resetr (enhanceinv, nx*ny, 1.)

c~~~~~~~~~~~~~~~~~~~~
      endif  ! ifrest  
c~~~~~~~~~~~~~~~~~~~~

      call zero (oceanmeltav,nx*ny)
      call zero (facemeltav,nx*ny)
      call zero (cliffmeltav,nx*ny)
      call zero (calviceav,nx*ny)
      call zero (runliqav,nx*ny)
      call zero (runfrozav,nx*ny)
      timefluxprev = timeice

#if defined (CRHMELREAD) || ( defined (FORCE34MYR) && defined (NOSED) && defined (EISANTA) )
c        If requested, read crhmel and enhanceinv from file written by 
c        previous CRHINV or ENHINV run, renamed to 'crhmelfilein'
c        (overwriting settings above)
      call readcrhmel (crhmel, enhanceinv, maskinit, sedpres)
#endif

#if defined (TEST2D) 
c        For mismipplus, reset crhmel to hard-coded settings, in case
c        overridden above by readres (ifrest <= 0)
      call basecoef_init (sedpres, maskinit, crhmel, timeice) 
#endif

c*******************************************
#if defined (SLMODEL) && ! defined (EISLINE)
c*******************************************

c............................
#  if defined (SLMODELRADIAL)
c............................

c        For radial sea-level model, initialize bedrock elevs to those 
c        at the start of the last outer iteration of a previous 
c        non-radial run (from its sl_hb_iter.<n> file).
c        Add current (initial) sealev, to convert from non-radial
c        (rel. to sea surface) to  radial (absolute) elevations.

      write (cfile,'(a,i3.3)') 'sl_hb_iter.final'
      open  (iuslit,file=cfile,status='old',form='unformatted')
      read (iuslit) hbstart
      close (iuslit)
      do j=1,ny
        do i=1,nx
          hb(i,j) = hbstart(i,j) + sealev
          topbed(i,j) = hb(i,j) - sedim(i,j) 
        enddo
      enddo
      sealev_init = sealev 
      write (iuslout2,'(a,f10.3/)') 'initial sealev = ',sealev_init
      call printmap (timeice, hb, 'initial hb (m)', 150.,
     *               iuslout2 ,0)
c......
#  else
c......
c        If 2-D sea-level (non-radial) model, set sealev to zero, and 
c        adjust all initial elevations accordingly. sealev is maintained
c        at zero in subr longterm throughout the run. (hb's returned in 
c        subr bedrock_sl are *relative* to current local sea level).  

#    if defined (SLMODEL)
      do j=1,ny
        do i=1,nx
          hb(i,j) = hb(i,j) - sealev
          hs(i,j) = hs(i,j) - sealev
          topbed(i,j) = topbed(i,j) - sealev
          topbedeq(i,j) = topbedeq(i,j) - sealev
        enddo
      enddo
      sealev = 0.
      sealev_init = 0.  ! not used - only used for radial s.l.
#    endif

#    if defined (SLMODELITER)
c     if sea-level model and outer-script iterating for initial hb, and 
c     if 2nd or later iteration, read hbstart written at end of previous
c     iter (see below), which overrides hb from initial or restart file
c     above. If 1st iter, just copy initial hb to hbstart.

      if (itersl.gt.1) then
        write (cfile,'(a,i3.3)') 'sl_hb_iter.',itersl-1
        open  (iuslit,file=cfile,status='old',form='unformatted')
        read (iuslit) hbstart
        close (iuslit)
        do j=1,ny
          do i=1,nx
            hb(i,j) = hbstart(i,j)
            topbed(i,j) = hb(i,j) - sedim(i,j) 
          enddo
        enddo
        call printmap (timeice, hbstart, 'hbstart (m)', 150.,
     *                 iuslout2 ,0)
      else
        call scopy (nx*ny, hb, 1, hbstart, 1)                 ! 1st iter
      endif
#    endif
c.......
#  endif
c.......

      call findwater (maskwater, indlake, npoilake, nlake,
     *                h, hb, hw, sealev, timeice)
      call adjustpres (maskwater, indlake, npoilake, nlake,
     *                 h, hb, hw, hs, sealev)

#  if defined (SLMODELITER)
c     write restart file at start, with slmodeliter adjustment to hb,etc
      call writeres (h, hs, hb, t, tracer, 
     *               baseperc, basefrml, heatb,
     *               sedim, tsed, wsed, heats, quarryacc, hbinit,
     *               topbed, topbedeq, equiload, tbed, hw, tw, 
     *               sedpres, crhmel, enhanceinv,
     *               ua, va, ub, vb, dfu, dfv, dbu, dbv,
#    if defined (MELANGE)
     *               hmel,
#    endif
     *               sealev, rco2inter, timeice, iczoom,
     *               'restart', 1)
#  endif

c*****
#endif
c*****

      facice = 0.
      facorb = 0.
      facco2 = 0.

c@@@@@@@@@@
      endif   ! ((istage.eq.1 .and. itry.eq.1) .or. ifblow)  (full init)
c@@@@@@@@@@

      ifblow = .false.

#if defined (ADAPTDT)
c@@@@@@@@@@@@@@@@@@@@@@@@
      if (itry.eq.1) then
c@@@@@@@@@@@@@@@@@@@@@@@@
c        For adaptive time stepping (with -DADAPTDT). Write
c        a restart file (restartadapt) at start of current stage,
c        to be read below for loop-back if blowup in inner loop.
c        Also write a marker line to all ascii output files, for use 
c        by post-processing programs to delete repeated output 
c        if blowups and re-running of stages occur (matching
c        with 'REWINDING' writemark call below).

      if (istage.eq.1) then
        call writemark ('INIT',    timeice, dtimeice, dtimeicesav)
      else
        call writemark ('REACHED', timeice, dtimeice, dtimeicesav)
      endif

c        Write initial restart file for start of first stage

      call writeres (h, hs, hb, t, tracer, 
     *               baseperc, basefrml, heatb,
     *               sedim, tsed, wsed, heats, quarryacc, hbinit,
     *               topbed, topbedeq, equiload, tbed, hw, tw, 
     *               sedpres, crhmel, enhanceinv,
     *               ua, va, ub, vb, dfu, dfv, dbu, dbv,
#  if defined (MELANGE)
     *               hmel,
#  endif
     *               sealev, rco2inter, timeice, iczoom,
     *               'restartadapt', 0)
      timeicesav = timeice
c@@@@@@@@@@
      endif  ! (itry.eq.1)  (write restartadapt)
c@@@@@@@@@@
#endif

c        Inner time loop (dtimeice) through current stage

c>>>>>>>>>>>>>>>>>>>>>
      do iloop=1,nloop
c>>>>>>>>>>>>>>>>>>>>>

c          Set timing flags used by subrs:
c          iffirse = first timestep whole experiment (only ifrest=-1,0)
c          iffirsr = first timestep this run (incl. ifrest=1), 
c                    or any stage after blowup (itry gt.1) 
c          iflast  = last timestep this run

        iffirse = 0
        if (ifrest.ne.1 .and. abs(timeice-nyearstart).lt.0.5*dtimeice) 
     *    iffirse = 1

        iffirsr = 0
        if ((istage.eq.1.or.itry.gt.1) .and. iloop.eq.1) iffirsr = 1

        iflast = 0
        if (istage.eq.nstage .and. iloop.eq.nloop) iflast = 1

c       for iceshow*, icetabbud:
        call mascalc (h, toth0, tota0)

        call scopy   (nx*ny, h, 1, hprev, 1)
        call scopy_i (nx*ny, maskwater, 1, maskwaterprev, 1)

        call zero (budgneg,  nx*ny)
        call zero (budgneg2, nx*ny)

#if defined (EISANTA) && defined (TRANSECTA) && defined (ZOOMGL)
c       Check if grounding line has migrated out of central zoomed bin.
c       If so, reset iczoom, grid vars, and regrid model vars:
        call rezoomgl (iczoom, maskwater,
     *                 h, hs, hb, tracer,
     *                 sedim, tsed, wsed,
     *                 topbed, tbed, hw, tw,
     *                 baseperc, basefrml, 
     *                 heatb, heats, heatw,
     *                 topbedeq, equiload, hbsd, crhmel,
     *                 u, v, ub, vb, uia, via, uit, vit, ua, va,
     *                 sih, siu, siv, dbu, dbv)

        call findwater (maskwater, indlake, npoilake, nlake,
     *                  h, hb, hw, sealev, timeice)
        call adjustpres (maskwater, indlake, npoilake, nlake,
     *                   h, hb, hw, hs, sealev)
#endif

#if defined (SEDIMENT) && defined (SEDTRACK)
        call scopy (nx*ny, sedim, 1, sedimold, 1)
#endif

c          Reset nested variables around domain edges, for 2 reasons:
c          (i) multiple nestdrive files, so time varying b.c.'s. 
c          (ii) reset tsed, wsed, tbed, etc, (although redundant for 
c               h,hb,t,ua,va, preserved around domain edges in icedyn 
c               and icetherm, as are hb,topbed,sedim in icebed,icesed).

#if defined (NESTING)
        n = max (1, nint (dtimegetnest/dtimeice))
        if (mod (iloop, n).eq.0) then
#  if defined (IFLAGNESTEQ)
          iflagnesteq = IFLAGNESTEQ
#  else
          iflagnesteq = 2
#  endif
          call getnest (h, hb, t, tbed, sedim, tsed, wsed,
     *                  topbedeq, equiload, sedpres, arcocn, maskwater,
     *                  ua, va, 1, iflagnesteq, 1, timeice)
        endif
#endif

c          Set long-term weighting/driving quantities

#if defined (FORCE34MYR) || defined (NHAGLAC)
c       calc totsl = esl fall (m) due to all modeled ice,as in icetab2d.
c       For 1st timestep, sealev=56 if ifrest=-1,else from restart file:
        totsl = 0.
        do j=1,ny  
          do i=1,nx
            if (h(i,j).gt.0. .and. maskwater(i,j).eq.0)
     *        totsl = totsl + darea(i,j)
     *              * (h(i,j)*(rhoice/rholiq) - max(0.,sealev-hb(i,j)))
          enddo
        enddo
        totsl = totsl / (0.7 * 4.*pi*(radius**2))          ! ~ocean area
#else
        totsl = 0.                                         ! not used
#endif
        call longterm (weirun, ro18, sealev, totsl, dtocn34m,
     *                 dtantann, dtantaus, dtantjan, dtseas,
     *                 rco2, ecc, obl, prec,
     *                 timeice, nyearstart, nyearend, ifrest, iutsin)

#if defined (LONGTERMONLY)
        go to 4440
#endif
#if defined (CO2INTER)
        rco2 = rco2inter
#endif

c          Set surface climate, budg[snow,rain,evap,melt], tsurf[i]

        call climate_all (budgsnow, budgrain, budgevap, budgmelt,
     *                    tsurf, tsurfi,
     *                    h, hs, hsbal, hb, maskwater, t,
     *                    tog, tmon, pmon, 
#if defined (LOVECLIP)
     *                    emon, qmon,
#endif
     *                    bedge, fedge, hedge,
     *                    ro18, sealev, weirun,
     *                    dtantann, dtantjan, dtseas,
     *                    rco2, ecc, obl, prec,
     *                    facice, facorb, facco2,
#if defined (CO2INTER)
     *                    toge, tage, prge, ruge, runofa,
#elif defined (RACMOCLIM)
     *                    ruge, runofa,
#endif
     *                    timeice, dtimeice, dtimeclim, dtimebudg,
     *                    iloop, iffirsr)

#if defined (CO2INTER)
        n = max (1, nint (dtimeco2/dtimeice))
        if (mod (iloop, n).eq.0) then
          call doco2 (tmon, pmon, runofa, h, maskwater,
     *                toge, tage, prge, ruge,
     *                rco2inter, ecc, obl, prec,
     *                timeice, dtimeco2, nyearstart)
        endif
#endif

c          Step lithospheric and asthenosheric bedrock response  
        
c--------------------
#if defined (BEDROCK) 
c--------------------
        n = max (1, nint (dtimebed/dtimeice))

#  if defined (SLMODEL) && ! defined (SLOLDBED) 
        if ( (mod (iloop, n).eq.0 .or. iffirse.eq.1) 
     *     ) then
          ncalls = nint(dtimestage*nstage/dtimebed)
          call bedrock_sl (h, hb, hw, topbed, sedim, dgeoid,
     *                 sealev, sealev_init, dtimebed, dtimeice, timeice,
     *                 ncalls, itersl)
#  else
        if (mod (iloop, n).eq.0) then
            call bedrock (h, hb, hw, topbed, topbedeq, deflect,
     *                    sedim, equiload, sealev, timeice, dtimebed)
#  endif

c         adjust hw and hs for change in bedrock elev (and so hb):
          call adjustpres (maskwater, indlake, npoilake, nlake,
     *                      h, hb, hw, hs, sealev)
        endif
c-----
#endif
c-----

#if defined (SLMODEL) && (! defined (BEDROCK) || defined (SLOLDBED))
        n = max (1, nint (dtimebed/dtimeice))
        if (mod (iloop, n).eq.0 .or. iffirse.eq.1) then
          write (iuslout,'(f12.2)') timeice + dtimeice
          write (cform,"( '(', i4, 'f10.3)' )") ny
          write (iuslout,cform) (1.e-3*yh(1,j),j=1,ny)
          write (iuslout,cform) (h(1,j),j=1,ny)
          write (iuslout,cform) (hb(1,j)-sealev+hw(1,j)+h(1,j),j=1,ny)
          write (iuslout,cform) (hb(1,j)-sealev,j=1,ny)
          write (iuslout,cform) (hw(1,j),j=1,ny)
        endif
#  endif

c          Step liquid water (ocean, open lakes, sub-ice hydrol/lakes)

        n = max (1, nint (dtimehyd/dtimeice))
        if (mod (iloop, n).eq.0) then

#if defined (MOVEW)
          call movewater (maskwater, h, hb, hw, tw, 
     *                    heatw, baseperc, basefrml, 
     *                    uw, vw, dtimehyd)
#else
          call zero (heatw,nx*ny)
          call zero (uw, (nxp+1)*(nyp+1))
          call zero (vw, (nxp+1)*(nyp+1))
#endif

c         reset ocean/lake locations and indices 

#if ! defined (NOMOVECOAST)
          call findwater (maskwater, indlake, npoilake, nlake,
     *                    h, hb, hw, sealev, timeice)
#endif

c         adjust hw and hs for change in hw:
          call adjustpres (maskwater, indlake, npoilake, nlake,
     *                     h, hb, hw, hs, sealev)
        endif

c          Step temperatures t,tw,tsed,tbed, set basefrml,wsed

        n = max (1, nint (dtimetherm/dtimeice))
        if (mod (iloop, n).eq.0 .or. iffirsr.eq.1) then
          call icetherm (h, hs, t, tracer, hw, tw, maskwater,
     *                   budgsnow, budgrain, budgevap, budgmelt, 
     *                   baseperc, basefrml, oceanmelt,
     *                   tsurf, tsurfi,
     *                   arhap, s1a0, s1b0, s2a0, s2b0, enhanceinv,
     *                   heati, heath, heatb, heatw, geoflux,
     *                   w, wa,
     *                   sedim, tsed, wsed, heats, tbed,
     *                   u, v, ub, vb, sih, siu, siv, dbu, dbv,
     *                   timeice, dtimetherm, 
     *                   mod(iloop,n).eq.0 .or. iffirse.eq.1
#  if defined (LINEM)
     *                  ,amismip
#  endif
     *                  )
        endif

#if defined (LAGRANGEO18)
c          Lagrangian traces 3-D ice cells back to surface, 
c          report mean d18O of ice sheet from d18O_precip
        n = max (1, nint (dtimetherm/dtimeice))
        if (mod (iloop, n).eq.0) then
          call lagrangeo18 (h, maskwater, u,v,w, timeice)
        endif
#endif

#if defined (CRHINV) || defined (ENHINV)
c          For inverse runs, adjust crhmel aiming towards hs -> hsbal
c          (if CRHINV), and/or adjust enhanceinv aiming towards equal 
c          fractional areas of frozen bed and crhmel=crhnos (if ENHINV)

        call basecoef_inv (h, hs, hsbal, hbsd, hw, maskwater, t,
     *                     baseperc, basefrml,
     *                     crhmel, enhanceinv, 
     *                     u, v, w,
     *                     timeice+dtimeice, dtimeice)

c          Set balance velocities ubal,vbal, for display only.
c          nb: not currently called or displayed.
c
c       call balvel (hbal, hsbal, budgall, ubal, vbal,
c    *               timeice+dtimeice, dtimeice)
#endif

c          Set basal coefficients, transferrring intrinsic crhmel 
c          (h-grid) to u,v grids (crh[u,v]), and accounting for 
c          current basal conditions (temp, h2o) if grounded.
c          Also set grounded fraction (fracg[u,v]) and interpolated 
c          ice thickness at grounding line (hg[u,v])

        call basecoef (h, hb, hbsd, hw, 
     *                 baseperc, basefrml, 
     *                 sedpres, quarryacc,
     *                 t, tsurfi, sealev,
     *                 crhmel, crhu, crhv, 
     *                 fracgu, fracgv, hgu, hgv,
     *                 timeice)

c          Sub-glacial sediment deformation. Step sediment thickness, 
c          calculate frictional heating for therm/vdif

c.....................
#if defined (SEDIMENT)
c.....................
        nloopsed = nint (dtimeice/min(dtimesed,dtimeice))
        dt = dtimeice/nloopsed
        do iloopsed = 1,nloopsed
          call sedflow (h, hb, hbsd, hw, 
     *                  topbed, topbedeq, sedim,
     *                  t, baseperc, basefrml,
     *                  dbu, dbv, crhu, crhv, fsedu, fsedv, ub, vb,
     *                  maskwater, quarryrate, quarrycoef, quarrysecxy,
     *                  quarryacc, sedpres,
     *                  pelagrate, heats, heatb, dt)
          do j=1,ny
            do i=1,nx
              quarryacc(i,j) = quarryacc(i,j) + quarryrate(i,j)*dt
              pelagacc(i,j) = pelagacc(i,j)   + pelagrate(i,j)*dt
            enddo
          enddo
        enddo

c          Move sediment within ocean and lakes

#  if ! defined (SLOPECRIT)
c       linear diffusion (coefturb):
        call sedocean (sedim, hb, topbed, maskwater, dtimeice)
#  else
c       slumping (slopecrit):
        call sedoceanold (sedim, hb, topbed, topbedeq, maskwater)
#  endif

c          Adjust hw and hs (for change in sedim, hence hb):

        call adjustpres (maskwater, indlake, npoilake, nlake, 
     *                   h, hb, hw, hs, sealev)
c.....
#endif
c.....

c          Set oceanmelt, facemelt, div, calvice 

        n = max (1, nint (dtimeocn/dtimeice))
        if (mod(iloop, n).eq.0 .or. iffirse.eq.1 .or. iffirsr.eq.1) then
          call ocean (h, hs, hb, hw, hbal, maskwater,
     *                budgsnow, budgevap, budgrain, 
     *                budgmelt, baseperc, basefrml, 
     *                oceanmelt, facemelt, msector,
     *                budgall, budgneg2, ubal, vbal,
     *                tracer, arhap, visc,
     *                div, calvice, crevliq, arcocn,
     *                ifnexo, fedge, hedge,  tempocn,
#if defined (MELANGE)
     *                hmel,
#endif
     *                ua, va, fracgu, fracgv,
     *                timeice, sealev, ro18, dtantjan, weirun, dtocn34m,
     *                dtimeocn, ifrest, iloop, iffirse, iffirsr, itry)
        endif

c          Set net ice budget for convenient passing to icedyn, output.
c          Nb: budgmelt and budgrain have contributed to baseperc
c          in icetherm/vdif (= budgmelt+budgrain if NOMOVET or if h=0).
c          Throughout, budg[snow,rain,...],base..., and oceanmelt are 
c          local "unit-column" values for ice, whereas facemelt is
c          a whole-grid-box average. So multiply budgsnow,etc, by 
c          ice areal fraction fedge (unless sfc atmos budget is +ve 
c          over land, then set fraction = 1).
c 
c          Don't include cliffmelt and calvice in budgall here.
c          calvice was set in ocean, but will be possibly modified in 
c          icedyn (cliffmelt). cliffmelt is set in icedyn. They are 
c          applied (to h) along with budgall in icedyn, and added below
c          to budgall for diagnostic output.  

        do j=1,ny
          do i=1,nx
c           if sfc atmos budget is +ve over land, set fraction to 1:
c           (also see comments in iceocean for shelf, fedge = 1)
            if (maskwater(i,j).eq.0 .and. bedge(i,j).gt.0.) then
              zfedge = 1.
            else
              zfedge = fedge(i,j)
            endif
            budgall(i,j) = 
     *                 (   budgsnow(i,j) - budgevap(i,j) 
     *                   + budgrain(i,j) - baseperc(i,j) + basefrml(i,j)
     *                   - oceanmelt(i,j)
     *                 ) * zfedge
     *                 - facemelt(i,j)
          enddo
        enddo

c          Step ice dynamics, update ice thickness h and sfc elev hs,
c          set heath,heatb,ub,vb for icetherm) 

#if defined (NOMOVEI)
        call scopy (nx*ny, h, 1, hsav, 1)
#endif

        call icedyn (h, hs, hb, hbsd, hw, sedim, t,
     *               arhap, s1a0, s1b0, s2a0, s2b0, heath, heatb, 
     *               budgall, budgneg, cliffmelt, calvice, crevliq,
     *               ifnexo, fedge, hedge, maskh, visc,
     *               sih, siu, siv, dfu, dfv, dbu, dbv,
     *               ua, va, uia, via, ub, vb, uit, vit, uadv, vadv, 
     *               hu, hv, uw, vw, masku, maskv, 
     *               crhu, crhv, fsedu, fsedv,
     *               fracgu, fracgv, 
     *               muind, mvind, thetau, thetav, 
     *               hgu, hgv, qgu, qgv,
     *               uschoof, vschoof, fschoofu, fschoofv,
     *               uagrid, vagrid,
     *               fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *               angnormu, angnormv,
#  if defined (MELANGE)
     *               ucimel, vcimel, hucimel, hvcimel, 
     *               bpumel, bpvmel, thetaumel, thetavmel, hmel,
#  endif
     *               numh, nuvtot, itera, iterc, totflow,
     *               maskwater, indlake, npoilake, nlake, sealev,
     *               timeice, dtimeice, ifrest, iffirse, iffirsr)

c       add ocean's budgneg2 to icedyn's budgneg:
        do j=1,ny
          do i=1,nx
            budgneg(i,j) = budgneg(i,j) - budgneg2(i,j)
          enddo
        enddo

#if defined (NOMOVEI)
        call scopy (nx*ny, hsav, 1, h, 1)
#endif

        call mascalc (h, toth, tota)

c       cliffmelt was set, and calvice possibly modified, in icedyn. 
c       Add to budgall for diagnostic output (whole-grid-box averages).
        do j=1,ny
          do i=1,nx
            budgall(i,j) = budgall(i,j) - cliffmelt(i,j) - calvice(i,j)
          enddo
        enddo

c::::::::::::::::::::
#if defined (MELANGE)
c::::::::::::::::::::

c          Step melange dynamics, update melange thickness hmel, etc.

# if defined (OCFACMEL)
        ocfacmel = OCFACMEL
#  else
        ocfacmel = 1.
#  endif

# if defined (ATFACMEL)
        atfacmel = ATFACMEL
#  else
        atfacmel = 1.
#  endif

        do j=1,ny
          do i=1,nx
#  if defined (JAKOBOCN)
c           zws = max (0., min (1., (arcocn(i,j)-20.)/20.))
            zws = max (0., min (1., (arcocn(i,j)-50.)/20.))
#    if defined (JAKOBOCNINT)
            zocint = JAKOBOCNINT
#    else
            zocint = 0.
#    endif
            oceanmeltmel(i,j) = zws*oceanmelt(i,j) + (1.-zws)*zocint
#  else
            oceanmeltmel(i,j) = oceanmelt(i,j)
#  endif

            budgallmel(i,j) = budgsnow(i,j) - budgevap(i,j)
     *                      - ocfacmel*oceanmeltmel(i,j)
     *                      - atfacmel*budgmelt(i,j)
          enddo
        enddo

        nloopmel = nint (dtimeice/min(dtimemel,dtimeice))
        dt = dtimeice/nloopmel
c...............................
        do iloopmel = 1,nloopmel
c...............................
          call icemel (hmel, hsmel, hwmel, fmel, maskmel,
     *                 budgallmel,  budgsupmel, 
     *                 budgnegmel,  budgadjmel, 
     *                 budgredmel,  budgsmomel,
     *                 curbackmel,  curthetmel,
     *                 curhmel,     curlmel, 
     *                 uamel, vamel, ucimel, vcimel, hucimel, hvcimel, 
     *                 thetaumel,  thetavmel,
     *                 uaboundmel, vaboundmel,
c
     *                 hb, hbsd, h, hedge, fedge, 
     *                 maskwater, tsurfi, ua, va, calvice,
     *                 sealev, timeice, dt, ifrest, iloopmel, nloopmel,
     *                 iffirse, iffirsr, nyeartabmel)
c............
        enddo
c............
c:::::
#endif
c:::::

        timeice = timeice + dtimeice

c          Write ascii output, 1-D profiles (1-D or 2-D domains)

        call iceshow1d (h, hs, hb, t,
     *                  basefrml, oceanmelt, budgall, budgneg,
     *                  facemelt, cliffmelt, div, arhap, calvice,
     *                  bedge, fedge, hedge,
     *                  tsurf, tsurfi, heati, heath, heatb,
     *                  w, wa, maskh,
     *                  sedim, tsed, wsed, heats, 
     *                  quarryrate, quarryacc, quarrycoef,
     *                  pelagrate, pelagacc,
     *                  topbed, topbedeq, tbed,
     *                  hw, tw, maskwater,
     *                  u, v, ua, va, uia, via, ub, vb,
     *                  uadv, vadv, hu, hv, dfu, dfv, dbu, dbv,
     *                  fracgu, fracgv, muind, mvind, thetau, thetav,
     *                  crhu, crhv, fsedu, fsedv, 
#if defined (MELANGE)
     *                  hmel, hwmel, fmel, maskmel, 
     *                  budgallmel, budgsupmel,
     *                  uamel, vamel,
     *                  bpumel, bpvmel, 
     *                  thetaumel, thetavmel, 
#endif
     *                  numh, nuvtot, itera, iterc,
     *                  toth0, tota0, toth, tota, totflow,
     *                  sealev, dtantann, dtantjan, dtseas, rco2,
     *                  timeice, dtimeice, weirun, 
     *                  iffirse, iflast, nyearout1d)

c          Write ascii output, 2-D maps  (2-D domains only - 
c          does nothing if nyearout2d = 0, as set for 1-D domains)

        call iceshow2d (h, hs, hb, hbsd, t, tracer, 
     *                  fedge, hedge, 
     *                  budgall, budgsnow, budgrain, budgmelt, 
     *                  budgevap, budgneg, 
     *                  baseperc, basefrml, oceanmelt, msector,
     *                  facemelt, cliffmelt,
#if defined (RACMOCLIM) 
     *                  runofa,
#endif
     *                  div, calvice, arcocn,
     *                  tsurf, tsurfi, tmon,
     *                  heati, heath, heatb, arhap,
     *                  w, wa, sedim, 
     *                  quarryrate, quarryacc, quarrycoef,
     *                  pelagrate, pelagacc,
     *                  topbed, topbedeq, deflect, geoflux,
     *                  hw, tw, maskwater, crhmel, enhanceinv,
     *                  hsbal, hbinit, tempocn,
     *                  u, v, ua, va, ub, vb,
     *                  uadv, vadv, hu, hv, dfu, dfv, dbu, dbv,
     *                  thetau, thetav, hgu, hgv,
     *                  fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                  angnormu, angnormv,
     *                  crhu, crhv, 
#if defined (MELANGE)
     *                  hmel, hwmel, fmel, maskmel, 
     *                  budgallmel, budgsupmel, budgadjmel, budgsmomel,
     *                  oceanmeltmel,
     *                  uamel, vamel,
     *                  bpumel, bpvmel, 
     *                  thetaumel, thetavmel,
     *                  uaboundmel, vaboundmel,
#endif
     *                  numh, nuvtot, itera, iterc,
     *                  toth0, tota0, toth, tota, totflow,
     *                  sealev, timeice, dtimeice, weirun,
     *                  iffirse, iflast, nyearout2d)


#if defined (SEDIMENT)
c       Accumulate and write (if time) sediment budget quantities:
        call sedbudg (sedim, quarryacc, topbed,msector,timeice,dtimeice,
     *                iffirse, iflast, nyeartot, nyearsedbud, iloop)

c       Track dated strata in sediment, write to output file (if time):
#  if defined (SEDTRACK)
        call sedtracking_eul (hb, sedim, sedimold, 
     *                        timeice, dtimeice, weirun, sealev, ifrest,
     *                        iffirse, iflast, nyearstart, nyeartot)
#  endif
#endif

c       1D-flowline plotting output (1-D or 2-D domains):

        call iceplot1d (h, hs, hb, t, tracer, 
     *                  budgall, basefrml, 
     *                  oceanmelt, facemelt, cliffmelt,
     *                  tsurf, tsurfi, heati, heath, heatb,
     *                  w, wa, maskh,
     *                  sedim, tsed, wsed, heats, 
     *                  quarryrate, quarryacc, quarrycoef,
     *                  pelagrate, pelagacc,
     *                  topbed, topbedeq, tbed,
     *                  hw, tw,
     *                  u, v, ua, va, uia, via, ub, vb, 
     *                  uadv, vadv, hu, hv, 
     *                  masku, maskv, crhu, crhv, 
     *                  fracgu, fracgv, thetau, thetav,
     *                  sealev, dtantann, dtantjan, dtseas, rco2, 
     *                  timeice, dtimeice, weirun, 
     *                  iffirse, iflast, nyeartot, nyearplot1d)

c       tabular output:

#if defined (EISLINE) && ! defined (FORCE34MYR)

        call icetab1d (h, hs, hb, t, tracer, 
     *                 budgall, basefrml, 
     *                 oceanmelt, facemelt, cliffmelt,
     *                 tsurf, tsurfi, heati, heath, heatb,
     *                 w, wa, maskh,
     *                 sedim, tsed, wsed, heats, 
     *                 topbed, topbedeq, tbed,
     *                 hw, tw, maskwater,
     *                 u, v, ua, va, uia, via, ub, vb, 
     *                 uadv, vadv, hu, hv, 
     *                 masku, maskv, muind, mvind, crhu, crhv, 
     *                 fracgu, fracgv, thetau, thetav, 
     *                 hgu, hgv, qgu, qgv,
     *                 uschoof, vschoof, fschoofu, fschoofv,
     *                 uagrid, vagrid,
     *                 fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                 sealev, dtantann, dtantjan, dtseas, rco2,
     *                 timeice, dtimeice, weirun,
     *                 iffirse, iflast, nyeartot, nyeartab)

#elif defined (EISANTA) || defined (NHA) || defined (CARB) || defined (ORDO) || defined (TEST2D) || defined (TESTMEL2D) || defined (FORCE34MYR)

#  if defined (LONGTERMONLY)
 4440   timeice = timeice + dtimeice
#  endif

        call icetab2d (h, hs, hb, hw, maskwater, topbedeq,
     *                 ro18, sealev,
     *                 dtantann, dtantaus, dtantjan, dtseas, rco2,
     *                 ecc, obl, prec, facice, facorb, facco2,
     *                 timeice, dtimeice, weirun, dtocn34m,
     *                 iffirse, iflast, nyeartot, nyeartab)

#  if defined (LONGTERMONLY)
        goto 4441
#  endif

        call icetabbud (h, hprev, maskwaterprev, bedge, fedge,
     *                  budgsnow, budgevap, budgrain, budgmelt,
     *                  baseperc, basefrml, 
     *                  oceanmelt, facemelt, cliffmelt, calvice,
     *                  budgneg, 
     *                  oceanmeltav, facemeltav, cliffmeltav,
     *                  calviceav, runliqav, runfrozav,
     *                  toth0, tota0, toth, tota, totflow, 
     *                  timeice,  dtimeice, iffirse, iflast,nyeartabbud)

#  if defined (MELANGE)
        call meltabbud (hmel, fmel, maskmel,
     *                  budgallmel,  budgsupmel, 
     *                  budgnegmel,  budgadjmel, 
     *                  budgredmel,  budgsmomel,
     *                  curbackmel,  curthetmel,
     *                  curhmel,     curlmel, 
     *                  timeice, dtimeice, iffirsr, iflast, nyeartabmel)
#  endif

#endif

#if defined (EISANTA) && ! defined (TRANSECTA)
        call icetabwais (h, hb, hw, maskwater,
     *                   ub, vb,
     *                   sealev, timeice, dtimeice, weirun,
     *                   iffirse, iflast, nyeartot, nyeartabwais)

#endif

#if defined (EISANTA) &&  defined (TABTROUGH)
        call icetabtrough (h, hb, hw, maskwater, fedge, hbinit,
     *                     sealev, timeice, dtimeice,
     *                     iffirse, iflast, nyeartot, nyeartabtro)
#endif

#if defined (TEST2D) && defined (TEST2DOUT)
      call icemismipplus (h, hs, hb, 
     *                    w, wa, maskh,
     *                    hw, maskwater, 
     *                    u, v, ua, va, uia, via, ub, vb,
     *                    uadv, vadv, hu, hv,
     *                    masku, maskv, muind, mvind, crhu, crhv,
     *                    fracgu, fracgv, thetau, thetav,
     *                    hgu, hgv, qgu, qgv,
     *                    uschoof, vschoof,
     *                    sealev, timeice, dtimeice, 
     *                    iffirse, iffirsr, iflast)
#endif

c          Write to  history file

        call writehis (h, hs, hb, hbsd, t, tracer, fedge,
     *                 budgall, budgsnow, budgrain, budgmelt,
     *                 baseperc, basefrml, 
     *                 oceanmelt, msector, facemelt, cliffmelt, calvice,
#if defined (RACMOCLIM)
     *                 runofa,
#endif
     *                 tsurf, tsurfi, tmon,
     *                 heati, heath, heatb,
     *                 w, wa, sedim, tsed, wsed, heats, dgeoid,
     *                 quarryrate, quarryacc, quarrycoef,
     *                 pelagrate, pelagacc,
     *                 topbed, topbedeq, deflect, 
     *                 tbed, geoflux,
     *                 hw, tw, maskwater, crhmel, enhanceinv,
     *                 oceanmeltav, facemeltav, cliffmeltav, 
     *                 calviceav, runliqav, runfrozav, timefluxprev,
     *                 u, v,  ua, va, ub, vb, uadv, vadv, 
     *                 fluxschu, fluxschv,
     *                 thetau,   thetav,
#if defined (MELANGE)
     *                 hmel, maskmel, budgallmel, budgsupmel,
     *                 uamel, vamel, 
     *                 bpumel, bpvmel,
     *                 thetaumel, thetavmel,
#endif
     *                 sealev, timeice, dtimeicesav, nyearstart,
     *                 iffirse, iflast, nyearhis, hislist, nhislist) 

c          Write to lagrange file (for tracing postprocessing)

#if defined (LAGRANGEWRITE)
        call writelag (h, hb, hs, hw, maskwater, u, v, w, timeice, 
     *                 dtimeice, nyearstart, iffirse, iflast, nyearlag)
#endif

c          Write restart fileun_O1_C1_B-5_E3
c          This is a user-requested restart file, with year appended
c          to the name ('restart_...'), not one written at end of each 
c          stage in case of blowups ('restart')

        if ( ( nyearres.ne.0 .and. 
#if defined (LOVECLIP)
     *         ( mod ( abs(timeice-nyearstart)+0.5*dtimeice,
#else
     *         ( mod ( abs(timeice)+0.5*dtimeice,
#endif
     *                 max(float(nyearres),dtimeice) ) .lt. dtimeice
     *           .or. iflast.eq.1
     *         )
     *       )
     *     ) then
          call writeres (h, hs, hb, t, tracer, 
     *                   baseperc, basefrml, heatb,
     *                   sedim, tsed, wsed, heats, quarryacc, hbinit, 
     *                   topbed, topbedeq, equiload, tbed, hw, tw, 
     *                   sedpres, crhmel, enhanceinv,
     *                   ua, va, ub, vb, dfu, dfv, dbu, dbv,
#if defined (MELANGE)
     *                   hmel,
#endif
     *                   sealev, rco2inter, timeice, iczoom,
#if defined (NORESTARTFILEDATES)
     *                   'restart', 0)
#else
     *                   'restart', 1)
#endif
        endif

c***************************************
#if defined (EISLINE) && defined (LINEM)
c***************************************
c          output mismip 50-yr file:

        if ( mod (abs(timeice)+0.5*dtimeice,50.) .lt. dtimeice ) then

          j=(ny+1)/2
          tvol = 0.
          do i=1,nx
            tvol = tvol + h(i,j)*dx(i,j)
          enddo

          j=(ny+1)/2
          xgl = 0.
          do i=nx-1,1,-1
            if (fracgu(i,j).gt.0.) then
              igl = i
              xgl = (1.-fracgu(i,j))*xh(i,j) + fracgu(i,j)*xh(i+1,j)
              goto 500
            endif
          enddo
          igl = 1
          xgl = 0.
  500     continue

          i1 = max (1, igl-1)
          i2 = igl
          i3 = min (nx, igl+1)

          iu0 = max (0, igl-2)
          iu1 = max (0, igl-1)
          iu2 = igl
          iu3 = min (nx, igl+1)

          hbg = (1.-fracgu(igl,j))*hb(igl,1) + fracgu(igl,j)*hb(igl+1,j)
          hg = (sealev-hbg)/rhor

          write (121, '( f11.1, f11.1, f15.1, 2f11.3,
     *                   3(7x, f11.1, 2f11.3, f11.1) )' )
     *      timeice, xgl, tvol, h(1,j), hg,
     *      xh(i1,j), h(i1,j), -hb(i1,j),
     *         0.5*(uadv(iu0,j)/dyu(iu0,j) + uadv(iu1,j)/dyu(iu1,j)),
     *      xh(i2,j), h(i2,j), -hb(i2,j),
     *         0.5*(uadv(iu1,j)/dyu(iu1,j) + uadv(iu2,j)/dyu(iu2,j)),
     *      xh(i3,j), h(i3,j), -hb(i3,j),
     *         0.5*(uadv(iu2,j)/dyu(iu2,j) + uadv(iu3,j)/dyu(iu3,j))
        endif

c          output mismip _ss and _f files, if last timestep:

        if (iflast.eq.1) then
          j=(ny+1)/2
          do i=1,nx
            write (122,'(f10.1, f10.3)') xh(i,j), h(i,j)
          enddo

c         write (123,'(f10.1, f10.1)') xgl, timeice
          write (123,'(3f10.3, i10)') 
     *      xgl/1.e3, xgexact/1.e3, (xgl-xgexact)/1.e3, nint(timeice)
        endif
c*****
#endif
c*****

c          Write nesting file (if not nested run)

#if ! defined (NESTING)
        if ( nyearnest.ne.0 .and. 
     *         ( mod (abs(timeice)+0.5*dtimeice,
     *                max(float(nyearnest),dtimeice)) .lt. dtimeice
     *           .or. iffirse.eq.1 
     *         )
     *     ) then
          call writenest (h, hb, hw, t, tbed, sedim, tsed, wsed, 
     *                    topbedeq, equiload, sedpres, arcocn, 
     *                    maskwater, ua, va, timeice)
        endif
#endif

c          Check for blowup (extreme ice values) setting ifblow, 
c          to re-run current stage with reduced timestep (itry loop)

        ifblow = .false.

        call checkblow (h, 1,nx, 1,ny, 1,1, -10., 10000., 0.,
     *                  'ice thickness (h)', timeice, ifblow)

        call checkblow (t, 1,nx, 1,ny, 0,nlevp,  -200., 20., tmelt,
     *                  'ice temperature (t)', timeice, ifblow)

c       ua,va sometimes exceeded (rarely, w. hi res), then recovers:
c       call checkblow (ua, 0,nxp, 0,nyp, 1,1, -100.e3, 100.e3, 0.,
c    *                  'ice veloc (ua)', timeice, ifblow)

c       call checkblow (va, 0,nxp, 0,nyp, 1,1, -100.e3, 100.e3, 0.,
c    *                  'ice veloc (va)', timeice, ifblow)

#if defined (MELANGE)
        call checkblow (hmel, 1,nx, 1,ny, 1,1, -10., 10000., 0.,
     *                  'melange thickness (hmel)', timeice, ifblow)

c       call checkblow (uamel, 0,nxp, 0,nyp, 1,1, -500.e3, 500.e3, 0.,
c    *                  'melange veloc (uamel)', timeice, ifblow)

c       call checkblow (vamel, 0,nxp, 0,nyp, 1,1, -500.e3, 500.e3, 0.,
c    *                  'melange veloc (vamel)', timeice, ifblow)
#endif

c@@@@@@@@@@@@@@@@@@@@@@@
        if (ifblow) then
c@@@@@@@@@@@@@@@@@@@@@@@
c          Blowup has occurred (detected by checkblow calls above).
c          If adaptdt, jump out of inner loop (goto 1000), to re-run 
c          current stage (next itry loop) with halved timesteps. First 
c          set ifrest=1, to include readres in first loop of re-run. 
c          If not adaptdt (ntry=1), or used all itry's, exit.    

c-------------------------------
          if (itry.lt.ntry) then
c-------------------------------
c              Write 'REWINDING' marker line to ascii output files

            call writemark ('REWINDING', timeicesav, dtimeice,
     *                      dtimeicenext)

            timeice = timeicesav
            ifrest = 1
            goto 1000
c-------------
          else
c-------------
c              Reached limit of halving timesteps, so abort run

            call writemark ('ABEND', timeice, dtimeice, 0.)
#if defined (ADAPTDT)
            write (ioterm,'(/a,i4,2(a,f8.2))') '*** ABEND: itry=', itry,
     *        '  dtimeice=',dtimeice, '  dtimeicesav=',dtimeicesav
#else
            write (ioterm,'(/a,f8.2)') '*** ABEND: dtimeice=',dtimeice
#endif
            call exit (-1)
c--------------
          endif
c--------------
c@@@@@@@@@@@@
        endif   ! (ifblow)
c@@@@@@@@@@@@

#if defined (LONGTERMONLY)
 4441   continue 
#endif

#if defined (RACMOCLIM)
c       stop
#endif

c******************
#if defined (RWRAP) 
c******************

c          for R-wrapper, set output variables for particular times

        zao = 0.7* 4.*pi*(radius**2)                ! ~global ocean area

c          Calculate total ice vol, grounded ice area, equiv. sea level.
c          If at right time,set arguments rwap_[plio,lig,mod,<yr>](1:3).
c          If first time through, do for modern obs (*bal), 
c          set arguments rwrap_obs(1:3).

        if (firstrwrap) then
          zvolobs = 0.
          zareobs = 0.
          zeslobs = 0.
          do j=1,ny
            do i=1,nx
              if (hbal(i,j).gt.0.) then
                zvolobs = zvolobs + hbal(i,j)*darea(i,j)
                if (maskbal(i,j).eq.0) then
                   zareobs = zareobs + darea(i,j)
                   zeslobs = zeslobs + darea(i,j)
     *               * (  hbal(i,j)*(rhoice/rholiq) 
     *                   - max(0.,sealev-(hsbal(i,j)-hbal(i,j))) )
                endif
              endif
            enddo
          enddo
          zeslobs = zeslobs / zao
          rwrap_obs(1) = zvolobs
          rwrap_obs(2) = zareobs
          rwrap_obs(3) = zeslobs
          do j=1,ny
            do i=1,nx
              rwrap_obs_h(i,j) = hbal(i,j)
              rwrap_obs_mask(i,j) = maskbal(i,j)
            enddo
          enddo
          firstrwrap = .false.
        endif

        zvol = 0.
        zare = 0.
        zesl = 0.
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              zvol = zvol + h(i,j)*darea(i,j)
              if (maskwater(i,j).eq.0) then
                 zare = zare + darea(i,j)
                 zesl = zesl + darea(i,j)
     *             * (h(i,j)*(rhoice/rholiq) - max(0.,sealev-hb(i,j)))
              endif
            endif
          enddo
         enddo
         zesl = zesl / zao

#  if defined (RWRAP_PLIO)

        if (abs(timeice-nyearend).lt.0.6*dtimeice) then
          rwrap_plio(1) = zvol
          rwrap_plio(2) = zare
          rwrap_plio(3) = zesl
        endif

#  elif defined (RWRAP_LIG)

c       if (abs(timeice-nyearend).lt.0.6*dtimeice) then 
        if (abs(timeice+120000.).lt.0.6*dtimeice) then 
          rwrap_lig(1) = zvol
          rwrap_lig(2) = zare
          rwrap_lig(3) = zesl
        endif

#  elif defined (RWRAP_LGMFUT)

        torig = 1950.

        if (abs(timeice + 21000.).lt.0.6*dtimeice) then
          rwrap_lgm(1) = zvol
          rwrap_lgm(2) = zare
          rwrap_lgm(3) = zesl
        endif
        if (abs(timeice - 0.).lt.0.6*dtimeice) then
          rwrap_mod(1) = zvol
          rwrap_mod(2) = zare
          rwrap_mod(3) = zesl
        endif
        if (abs(timeice + torig - 2100.).lt.0.6*dtimeice) then
          rwrap_2100(1) = zvol
          rwrap_2100(2) = zare
          rwrap_2100(3) = zesl
        endif
        if (abs(timeice + torig - 2200.).lt.0.6*dtimeice) then
          rwrap_2200(1) = zvol
          rwrap_2200(2) = zare
          rwrap_2200(3) = zesl
        endif
        if (abs(timeice + torig - 2300.).lt.0.6*dtimeice) then
          rwrap_2300(1) = zvol
          rwrap_2300(2) = zare
          rwrap_2300(3) = zesl
        endif
        if (abs(timeice + torig - 2400.).lt.0.6*dtimeice) then
          rwrap_2400(1) = zvol
          rwrap_2400(2) = zare
          rwrap_2400(3) = zesl
        endif
        if (abs(timeice + torig - 2500.).lt.0.6*dtimeice) then
          rwrap_2500(1) = zvol
          rwrap_2500(2) = zare
          rwrap_2500(3) = zesl
        endif

c       For modern, also set 2-D map arguments (rwrap_mod_[h,mask]),
c       and rms ice thickness minus obs, where any ice (rwrap_mod_rms):
        if (abs(timeice - 0.).lt.0.6*dtimeice) then
          zrms = 0.
          zarr = 0.
          do j=1,ny
            do i=1,nx
              rwrap_mod_h(i,j) = h(i,j)
              rwrap_mod_mask(i,j) = maskwater(i,j)
              if (h(i,j).gt.0. .or. hbal(i,j).gt.0.) then
                zrms = zrms + darea(i,j)*(h(i,j)-hbal(i,j))**2
                zarr = zarr + darea(i,j)
              endif
            enddo
          enddo
          zrms = sqrt (zrms/max(zarr,.001)) 
          rwrap_mod_rms = zrms
        endif
#  endif
c*****
#endif
c*****

c          Write diagnostic hysometry (to iuhyp), start+end of run only

#if defined (EISANTA)
        if (iffirse.eq.1 .or. iflast.eq.1) then
          call hypsom (hb, h, maskwater, sealev, timeice,iffirse,iflast)
        endif
#endif

c>>>>>>>>>>
      enddo  ! inner time loop
c>>>>>>>>>>

      goto 2000

 1000 continue  ! for failed jump out of inner time loop  (ifblow)

c>>>>>>>>>>
      enddo  ! adaptive time-step loop
c>>>>>>>>>>

 2000 continue ! for successful jump out of adaptive timestep loop

c       Successful completion of stage - go to next stage

c>>>>>>>>>>
      enddo  ! outer stage loop
c>>>>>>>>>>

c************************
#if defined (SLMODELITER)
c************************
c        If sea-level model and outer-script iterating,
c        compute rms and max |dhb| (dhbrms, dhbmax).
c        Als set new initial hb (hbstart) for next iteration. 

      zarea = 0.
      dhbrms = 0.
      dhbmax = 0.
      do j=1,ny
        do i=1,nx
          dhb(i,j) = hb(i,j) - hbinit(i,j)
          hbstart(i,j) = hbstart(i,j) - dhb(i,j) 
          zarea  = zarea  + darea(i,j)
          dhbrms = dhbrms + darea(i,j)*dhb(i,j)**2
          dhbmax = max (dhbmax, abs(dhb(i,j)))
        enddo
      enddo
      dhbrms = sqrt (dhbrms/zarea)
      write (iuslout2,'(/a,2f8.3/)') 'dhbrms, dhbmax= ', dhbrms, dhbmax
      write (ioterm,  '(/a,2f8.3/)') 'dhbrms, dhbmax= ', dhbrms, dhbmax

c        Display final hb, observed modern hb (hbinit), difference 
c        (dhb), new initial hb for next iter (hbstart)

      call printmap (timeice, hb,    'hb (m)',150.,
     *               iuslout2 ,0)
      call printmap (timeice, hbinit,'hbinit (m)',150.,
     *               iuslout2 ,0)
      call printmap (timeice, dhb,   'dhb (m)',20.,
     *               iuslout2 ,0)
      call printmap (timeice, hbstart, 'hbstart (new) (m)',150.,
     *               iuslout2 ,0)

c        Write next-iter's hbstart (read at start of next iter above, 
c        to set hb)

      write (cfile,'(a,i3.3)') 'sl_hb_iter.',itersl
      open  (iuslit,file=cfile,status='unknown')
      close (iuslit,status='delete')
      open  (iuslit,file=cfile,status='new',form='unformatted')
      write (iuslit) hbstart
      close (iuslit)

      open (iuokend, file='okend', status='new') 
      close (iuokend)

c     exit with status=200 for do[slmodel,sed]iter to detect convergence
      if (dhbrms.lt.1.0) then     ! convergence criterion using dhbrms
c     if (dhbmax.lt.3.0) then     ! convergence criterion using dhbmax
        call exit (200)
      else
        call exit (0)
      endif
c*****
#endif
c*****

c**********************************************************************
#if defined (SEDIMENT) && defined (FORCE34MYR) && defined (QUARRYWRITE)
c**********************************************************************

c        Compute new initial bed topography, attempting to adjust 
c        final (modern) ice-free rebounded bedrock to observed 
c        (hbrebmod) at end of next iter run. Use Newton's Method
c        over previous iter and this one, using (hbinitprev, hbrebprev)
c        and (hbinit, hbreb), to estimate next iter's hbinit at which 
c        the linearly interpolated hbreb(next) = hbrebmod - zqd,
c        where zqd is the antoicipated effect of changing quarrying
c        coefficient (quarrysec[xy]) throughout next run. zqd is
c        estimated for each point from this run's coecwffs and 
c        quarryacc, sedim. 
c
c        Writes next iter's hbinintprev, hbprev, hbinit, quarrysec 
c        to quarryfile (subr writequarry), to be read at start of 
c        next iter run if QUARRYREAD.

      call bedmap2_aggreg (hmod,   'thk')         ! modern obs ice thick
      call bedmap2_aggreg (hbmod, 'topg')         ! modern obs bed elev

      call scopy (nx*ny, hbmod, 1, hbrebmod, 1)
      call scopy (nx*ny, hb,    1, hbreb,    1)

c        Simple rebound to ice-free equil where grounded, for modern 
c        hbrebmod (s.l. = 0), and model hbreb (s.l. = sealev)

#  if ! defined (QUARRYINVNOREB)
      do j=1,ny
        do i=1,nx
          zovmod = rhoice*hmod(i,j) - rholiq*max(0.,0.-hbmod(i,j))
          if (zovmod.gt.0.) hbrebmod(i,j) = hbmod(i,j) + zovmod/rhoman

          zov = rhoice*h(i,j) - rholiq*max(0.,sealev-hb(i,j)) 
          if (zov.gt.0.) hbreb(i,j) = hb(i,j) + zov/rhoman
        enddo
      enddo
#  endif

c        Save this run's quarry coeff arrays before modified below

      call scopy (nsedsec, quarrysec, 1, quarrysecold, 1)
      call scopy (nx*ny, quarrysecxy, 1, quarrysecxyold, 1)

c        Calc domain totals

      ztotsedobs = 0.
      do m=1,nsedsec
        ztotsedobs = ztotsedobs
     *             + totsedsecobs(m)  ! (m3) sum of Wilson etal.,P3,2012
      enddo
      ztotsed = 0.
      ztotqac = 0.
      do j=1,ny
        do i=1,nx
          ztotsed = ztotsed + sedim(i,j)*darea(i,j)
          ztotqac = ztotqac + quarryacc(i,j)*darea(i,j)
        enddo
      enddo

c       Totals per sector (tot[sed,qac]sec, using msedsec in
c       comicegrid.h, set to indices of Wilson et al. (P3, 2012) 
c       basins in calcwilsonbasins):

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

c       Adjust quarrying coefficients (quarrysec[xy]) aiming to make
c       next run's total modern sedment vol. = obs modern sed vol 

#  if defined (QUARRYINVGLOB)
c     adjust continental quarrying coeff (all same value):  
      do m=1,nsedsec
        quarrysec(m) = quarrysec(m) * (ztotsedobs/max(ztotsed,1.))
      enddo
#  else
c     adjust indiv (Wilson et al.) sector quarrying coeffs:
      do m=1,nsedsec
        quarrysec(m) = quarrysec(m) 
     *               * (totsedsecobs(m) / max(totsedsec(m),1.))
      enddo
#  endif

      do m=1,nsedsec
        quarrysec(m) = max (1.e-16, min (1.e-7, quarrysec(m)))  ! safety
      enddo

c     transfer coeffs from indiv sector values to 2-D map, then smooth:
      do j=1,ny
        do i=1,nx
          quarrysecxy(i,j) = quarrysec(msedsec(i,j))
        enddo
      enddo
      call smooth9pt (quarrysecxy, nx, ny)

c        Set next run's initial bed elevs (hbinit), which will set next
c        run's topbedeq via readquarry (hb->set topbedeq) in iceinit.
c        The changes here aim to make next run's modern rebounded bed 
c        elev = obs,using Netwon's Method after 1st iter (1st iter
c        does not have QUARRYREAD), allowing both for local error at 
c        end of this run (hbrebmod-hbreb), and anticipated change in 
c        final elev (zqd) due to change in quarrysecxy above and its 
c        effect through the run on local quarryacc and sedim, allowing 
c        for isostatic compensation in mantle (1.- rho.../rhoman).

      do j=1,ny
        do i=1,nx
          zqd = (quarrysecxy(i,j)/quarrysecxyold(i,j) - 1.)
     *        * ( - quarryacc(i,j)*(rhosed/rhobed)*(1.-rhobed/rhoman)
     *            + sedim(i,j)*(1.-rhosed/rhoman) )

#  if defined (QUARRYREAD)       
c         Newton's method:
          zdhbi = hbinit(i,j)-hbinitprev(i,j)
          zdhbi = sign (max(1.,abs(zdhbi)), zdhbi)
          zslope = (hbreb(i,j)-hbrebprev(i,j)) / zdhbi
          zslope = sign (max(0.1,abs(zslope)), zslope)
          zd = (hbrebmod(i,j) - hbreb(i,j) - zqd) / zslope
          zd = max (-500., min (500., zd))
#  else
c         first iter only, assuming change in hbinit -> same in hbreb 
          zd = hbrebmod(i,j) - hbreb(i,j) - zqd
#  endif

c         transfer this run's hbinit, hbreb, to hbinitprev, hbrebprev,
c         to be written to quarryfile:
          hbinitprev(i,j) = hbinit(i,j)
          hbrebprev(i,j) = hbreb(i,j)

c         set next iter's hbinit:
          hbinit(i,j) = hbinit(i,j) + zd
        enddo
      enddo

c        Smooth next iter's hbinit

c777  call smooth9pt (hbinit, nx, ny)

c        Write diagnostics to stdout

      zdhbrms = 0.
      zdhbmax = 0.
      zdhbrmsinit = 0.
      zdhbmaxinit = 0.
      ztotare = 0.
      do j=1,ny
        do i=1,nx
          work(i,j)  = hbrebmod(i,j)  - hbreb(i,j)
          zdhbrms = zdhbrms + darea(i,j)*work(i,j)**2
          zdhbmax = max (zdhbmax, abs(work(i,j)))

          work2(i,j) = hbinit(i,j) - hbinitprev(i,j)
          zdhbrmsinit = zdhbrmsinit + darea(i,j)*work2(i,j)**2
          zdhbmaxinit = max (zdhbmaxinit, abs(work2(i,j)))
          ztotare = ztotare + darea(i,j)
        enddo
      enddo
      zdhbrms     = sqrt (zdhbrms    /ztotare)
      zdhbrmsinit = sqrt (zdhbrmsinit/ztotare)

      zz = 1.e-9                                     ! m3 to km3
      zz2 = 1.e-15                                   ! m3 to million km3
      write (ioterm,
     *  '(/a/a,2f10.2/a,2f10.2/a,3f15.1
     *    /a,8f7.3,/a,8f7.3/a,8f7.3/a,8f7.3/a,8f7.3)')
     * 'QUARRYWRITE:',
     * ' dhb[rms,max]end  (m)   =', zdhbrms,       zdhbmax,
     * ' dhb[rms,max]next (m)   =', zdhbrmsinit,   zdhbmaxinit,
     * ' qaco,sedo,sedobs (km3) =', ztotqac*zz,ztotsed*zz,ztotsedobs*zz,
     * ' totqaco (Mkm3):', (totqacsec   (m)*zz2, m=1,nsedsec),
     * ' totsedo (Mkm3):', (totsedsec   (m)*zz2, m=1,nsedsec),
     * ' totsedob(Mkm3):', (totsedsecobs(m)*zz2, m=1,nsedsec),
     * ' qseco   *1.e10:', (quarrysecold(m)*1.e10, m=1,nsedsec),
     * ' qsecn   *1.e10:', (quarrysec   (m)*1.e10, m=1,nsedsec)

      call printmap (timeice,hbrebmod,'hbrebmod (m)',  150., iuout2d ,0)
      call printmap (timeice,hbreb,   'hbreb (m)',     150., iuout2d ,0)
      call printmap (timeice,work,'hbrebmod-hbreb(m)',  50., iuout2d ,0)
      call printmap (timeice,work2,'hbinit-hbinitprev(m)',50.,iuout2d,0)
      call printmap (timeice,hbinit,  'hbinit (m)',    150., iuout2d ,0)
      call printmap (timeice,hbinit,  'hbinit (m)',    500., iuout2d ,0)

c        Write hb... fields and coeffs (read at start if -DQUARRYREAD)

      call writequarry (hbinitprev, hbrebprev, hbinit, quarrysec, 
     *                  timeice)
c*****
#endif
c*****

c*********************
#if defined (LOVECLIP)
c*********************
c        Write final nc file to be read by Loveclim next segment

      call writeloveclim (hedge, hs, hb, fedge, maskwater,
     *                    runliqav, runfrozav, timefluxprev,  
     *                    sealev, timeice, nyearstart)
c*****
#endif
c*****

c     create stub file to signify normal end of run:
      open (iuokend, file='okend', status='new') 

#if defined (RWRAP)
      return
#else
      stop
#endif
      end

c-----------------------------------------------------------------------

      subroutine longterm (weirun, ro18, sealev, totsl, dtocn34m, 
     *                     dtantann, dtantaus, dtantjan, dtseas,
     *                     rco2, ecc, obl, prec,
     *                     timeicein, nyearstart,nyearend,ifrest,iutsin)

c        Sets long-term variables (sealev, dtant[ann,jan], dtseas, 
c        rco2, weirun):
c        sealev   = sea level (m)
c        dtantann = annual mean insol at 80 S - modern * 0.1 deg C/W/m2 
c        dtantjan = January     insol at 80 S - modern * 0.1 deg C/W/m2 
c        dtseas   = July - January insol * 0.1 deg C/W/m2  
c                   (or * 0.07 for SIMPEO, not 0.1, in subr laskar)
c        rco2     = atmos CO2 level / preindustrial  
c        weirun   = run weight (used in iceclimall for LINE*, 
c                   in ocean for LINE* and EISANTA).

#include <comicepath.h>

      parameter (pi     = 3.14159265358979)

c     for cyclic EISANTA tests (mod-LGM-mod-ig-mod):
      parameter (ntim=10)
      dimension ztim(ntim), wtim(ntim)
c     data ztim /0., 30., 40., 70., 80., 110., 120., 150., 160., 190./
c     data wtim /1.,  1.,  0.,  0.,  1.,   1.,   2.,   2.,   1.,   1./
      data ztim /0.,  5., 15., 20., 30.,  35.,  45.,  50.,  60.,  65./
      data wtim /1.,  1.,  0.,  0.,  1.,   1.,   2.,   2.,   1.,   1./
      save ztim, wtim

#if defined (RCMMATRIX) && defined (RCMMATCO2)
c     stmt fn for natural atmos co2 decay after co2 emissions end, 
c     from Archer et al,1997,GRL,Eq.3 (used in Charbit et al,2008):
#  if defined (NOARCHER)
      farch (ztb) = 1.
#  elif defined (RAPIDARCHER)
      farch (ztb) = exp (-ztb/100.)
#  else
      farch (ztb) =   0.750*exp(-ztb/365.)  + .135*exp(-ztb/5.5e3)
     *              + 0.035*exp(-ztb/8.2e3) + .080*exp(-ztb/200.e3)
#  endif
#endif

#if defined (NESTING) && ! defined (NESTEVOL)
c     avoid slow changes during nested run, use start time throughout
c     (unless time dependence specified with NESTEVOL):
      timeice = nyearstart
#elif defined (TSHIFT_JOHNSON)
      timeice = timeicein - TSHIFT_JOHNSON
#else
      timeice = timeicein 
#endif


c        Set defaults

      dtantann = 0.
      dtantjan = 0.
      dtseas   = 0.
      ecc  = 0.016706             ! modern
      obl  = 23.4377 * pi/180.    ! modern
      prec = 77.0613 * pi/180.    ! modern
      ro18 =  0.
#if ! defined (SLMODEL)
      sealev = 0.
#endif
      rco2 = 1.
      weirun = 1.
      dtocn34m = 0.

c==============================================
#if defined (EISLINE) && ! defined (FORCE34MYR)
c==============================================

c---
c     constant:
      weirun = 0.
c     weirun = 0.5
c     weirun = 1.
      sealev = 0. - 130.*weirun
c---
c     alternating steps:
c     if ( mod ( int((timeice-nyearstart)*1.e-3/20.), 2 ).eq.1 ) then
c       weirun = 0.
c     else
c       weirun = 1.
c     endif
c     sealev = 0. - 130.*weirun
c---
c     Eismint g.l. steps:
c     sealev0 =    0.
c     sealev0 = -100.
c     sealev0 =  100.
c     dsealev = -125.
c     dsealev =    0.
c     dsealev =  125.
c     kyear = nint(timeice-nyearstart)
c     if (kyear.le.50000) then
c       weirun = 0. 
c     else if (kyear.le.100000) then
c       weirun = min((kyear-50000)/1000., 1.)
c     else 
c       weirun = 1.- min((kyear-100000)/1000.,1.)
c     endif
c     sealev = sealev0 + dsealev *weirun
c-----
# if defined (LINEB)
c-----
c     Quaternary sawtooth cycles (weirun: 0=modern, 1=lgm):
      zper = 120000.
      zsub = 105000.
      zt = mod (timeice-nyearstart, zper)
      if (zt.lt.zsub) then 
        weirun = zt/zsub
      else
        weirun = 1. - (zt-zsub)/(zper-zsub)
      endif

c     gradual s.l. lowering, then gradual rising:
c     ztm = timeice - nyearstart
c     if(ztm.le.125000.) then
c       weirun = max (0., min (1.,    (ztm         )/125000.))
c     else
c       weirun = max (0., min (1., 1.-(ztm -125000.)/125000.))
c     endif

c     sealev = 0. - 130.*weirun  ! old
c     sealev = 0.                ! old

#  if defined (SIDEDRAG)
c     sealev = 180. -   50.*weirun
      sealev = 250. -  125.*weirun
#  else
      sealev =  30. -   95.*weirun
#  endif

c     reverse slope:
c     ztm = timeice - nyearstart
c     weirun = max (0., min (1., (ztm - 50000.)/50000.))
c#  if defined (SIDEDRAG)
c     sealev = 180. - 50.*weirun  ! sl falling
c     sealev = 130. + 50.*weirun  ! sl rising
c     sealev = 185. - 60.*weirun  ! sl sawtooth
c#  else
c     sealev = -55. -  10.*weirun ! sl falling
c     sealev =  20. +  10.*weirun ! sl rising
c     sealev =  35. - 105.*weirun ! sl sawtooth
c#  endif

c-----
# elif defined (LINEC)
c-----
c     Cenozoic "orbital" cycles (weirun: 0=warm, 1=cold):

c     weirun = 0.5 
c     sealev = 0. - 130.*weirun

c     zt = timeice - nyearstart 
c     weirun = 0.5 + 0.5*sin(2.*3.1415927*zt/100000.)
c     sealev = 100. - 200.*weirun

c     ztm = mod (timeice-nyearstart, 200000.)
c     weirun = min ( (ztm - 0.)/5000., (105000.-ztm)/5000. )
c     weirun = max (0., min (1., weirun))
c     sealev = 100. - 200.*weirun

c     ztm = mod (timeice-nyearstart, 100000.)
c     weirun = max ( (90000.-ztm)/5000., (ztm-95000.)/5000. )
c     weirun = 0.5* max (0., min (1., weirun))
c     weirun = 0.5
c     sealev = 300.*(0.5-weirun) +  (-130.)*weirun

      weirun = 0.
#  if defined (FLATANT)
      sealev = 0.
#  else
      sealev = -5000.
#  endif

c     if (timeice-nyearstart.lt.9900000.) then
c       sealev = -65.
c     else    
c       sealev = 150.
c     endif   
c     weirun = 0.5

c-----
# elif defined (LINED) 
c-----
c     basic mountain/ocean profile, just quarrying+transport

      ztm = mod (timeice-nyearstart, 1000000.)
      if (ztm.lt.100000.) then
        weirun = 0.5 * (ztm/100000.)
      else if (ztm.gt.800000.) then
        weirun = 0.5 * (1. - min (1.,(ztm-800000.)/100000.))
      else
        weirun = 0.5
      endif
      sealev = 0.

c     weirun = 0.5
c     if (timeice-nyearstart.gt.9900000.) weirun = 0.
c     sealev = 0.

c----
# elif defined (LINEE)
c----
c     Richard Alley sediment wedge WAIS tests
      weirun = 0.5
#  if defined (SEALEVCHANGE)
#   if defined (SEALEVFALL)
      sealev = -40. - SEALEVFALL
#   elif defined (SEALEVRISE)
c     sealev = -40. + 30   ! + 0 5 10 20 30
      sealev = -40. + SEALEVRISE
#   endif
#  elif defined (RESTARTE)
c     sealev = 0.
      sealev = -40. ! -35
#  elif defined (INITEE)
      zwei =  max (0., min (1., (timeice-1000.)/1000.))
ca    sealev = -320.*(1.-zwei) - 0.*zwei   ! a
c     sealev = -500.*(1.-zwei) - 0.*zwei
      sealev = -500.*(1.-zwei) - 40.*zwei  ! -35
#  else
      zwei =  max (0., min (1., (timeice-1000.)/(nyearend-1000.)))
      sealev = -320.*(1.-zwei) + 100.*    zwei   ! b,c
c     sealev = -320.*(   zwei) + 100.*(1.-zwei)  ! d,e
c     sealev = min (sealev, -75.)
c     sealev = -75. + 00.*min (1., timeice/2000.) ! 49.9 km
c     sealev = -75. +  5.*min (1., timeice/2000.) ! 
c     sealev = -75. + 10.*min (1., timeice/2000.) ! 20.9 km
c     sealev = -75. + 25.*min (1., timeice/2000.) ! 10.3 km
#  endif

c-----
# elif defined (LINEF)
c-----
      weirun = 0.5

c     zwei =  max (0., min (1., (timeice-20000.)/100.))
c     sealev = -500.*(1.-zwei) - 0.*zwei
c     sealev = -501.*(1.-zwei) - 200.*zwei

c     sealev = -310.  ! on slope
c     sealev = -320.  ! on shelf
c     sealev = -330.  ! to rh edge

c     zwei =  max (0., min (1., (timeice - 0.)/200000.))
c     sealev = -300.*(1.-zwei) - 330.*zwei

c     if (timeice.lt.50000.) then
c       zwei =  max (0., min (1., (timeice - 0.)/50000.))
c       sealev = -220.
c     else if (timeice.lt.100000.) then
c       zwei =  max (0., min (1., (timeice - 50000.)/50000.))
c       sealev = -220.*(1.-zwei) - 233.*zwei ! 240
c     else if (timeice.lt.150000.) then
c       zwei =  max (0., min (1., (timeice - 100000.)/50000.))
c       sealev = -233.
c     else if (timeice.lt.250000.) then
c       zwei =  max (0., min (1., (timeice - 100000.)/50000.))
c       sealev = -233.
c     else
c       zwei =  max (0., min (1., (timeice - 150000.)/50000.))
c       sealev = -233.*(1.-zwei) - 130.*zwei
c     endif

c     sealev = -150.
c     sealev = -150. + 1.*(timeice-40000)/2000.

c     sealev = 0. 
ca    sealev = -100. 
cb    sealev = 0. 
cc    sealev = -200. 

c     sealev = -100.
c     if (timeice.gt.100000.) sealev = -150.

      sealev = 0. 
c     if (timeice.gt.100000. .and. timeice.le.200000.) sealev = -125.
c     if (timeice.gt.100000. .and. timeice.le.200000.) sealev =  125.

c     schoof a,b:
c     if (timeice.gt.100000.) sealev = -125.
c     if (timeice.gt.200000.) sealev =  100.

c     schoof c:
c     zk = timeice/1000.
c     if (zk.lt.200.) sealev = 150. - 350.*zk/200.
c     if (zk.gt.200. .and. zk.lt.300.) sealev = -200.
c     if (zk.gt.300.) sealev = -200. + 350.*(zk-300.)/200.

c     if  (timeice.lt.50000.) then
c       sealev = 100.   
c     else if  (timeice.lt.100000.) then
c       sealev = -125.
c     else
c       sealev = 100.
c     endif

c     sealev = 0.

c-----
# elif defined (LINEG)
c-----

      weirun = 0.5
      sealev = 0.

c-----
# elif defined (LINEH) || defined (LINEH2)
c-----

      weirun = 0.
      sealev = 0.

      timea = -37.0e6
      timeb = -33.0e6
      if (timeice.lt.timeb) then
        rco2 = 6. - 1.*((timeice-timea)/1.e6)
      else
        rco2 = 6. - 1.*((timeb  -timea)/1.e6)
     *            + 1.*((timeice-timeb)/1.e6)
      endif
c     rco2 = 3.
      call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
     *             timeice, iutsin)

c-----
# elif defined (LINEH3) || defined (LINEH4)
c-----

      weirun = 0.
      sealev = 0.
      ztimeice = 0.
      call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
     *             ztimeice, iutsin)

c-----
# elif defined (LINEM) 
c-----

      weirun = 0.5
      sealev = 0.
      ztimeice = 0.
      call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
     *             ztimeice, iutsin)

c-----
# elif defined (SLMODEL)
c-----

c     for slmodel, 1-d profile:

      weirun = 0.

#   if defined (SLSEALEVEL) 
      sealev = SLSEALEVEL
#   else
      sealev = 0.
#   endif

#   if defined (SLPROFILEB) && ! defined (SLRESTART)
c     with all hb < 0, allow a large ice sheet to grow initially:
      if (timeice.lt.100000.) then
        sealev = -500.
      else
        sealev = -500. 
     *         + min (1.,(timeice-100000.)/10000.) * (SLSEALEVEL+500.) 
      endif
#   endif

      ztimeice = 0.
      call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
     *             ztimeice, iutsin)

c-----
# else
c-----

      weirun = 0.5
      sealev = 0.

c-----
# endif
c-----

c=====================
#elif defined (ASYNCH)
c=====================

      weirun = 1.
      sealev = 0.
      dtantann = 0.
      dtantaus = 0.
      dtantjan = 0.
      dtseas   = 0. ! not used

c=====================================================================
#elif defined (EISANTA) || (defined (EISLINE) && defined (FORCE34MYR))
c=====================================================================

c     dtantann = annual mean insol at alatlask-modern*dtdf_ann,  C/W/m2 
c     dtantaus = austral spring,insol alatlask-modern*dtdf_aus,  C/W/m2 
c     dtantjan = January     insol at alatlask-modern*dtdf_jan,  C/W/m2 
c     dtseas   = [July-January insol]      - modern * dtdf_seas, C/W/m2 
c                (dtdf_* = 0.07 for SIMPEO, else 0.1, in subr laskar)
c     rco2     = atmos CO2 level / preindustrial  
c     ro18     = modern (0) to LGM (1) weight for atmos param 
c     sealev   = sea level (m)
c     weirun   = run weight (used in ocean, via zclim)
 
      weirun = -1.e20   ! set at end of EISANTA section, if not before

      sealevin = sealev ! used to recover sealev for SLMODELRADIAL below

c--------------------------
#  if defined (FORCEPLEIST)
c--------------------------

      timeicex = timeice

c---

c     timeicex = timeice    ! for mod, or hot (fut)
c     timeicex = -20000.    ! for lgm
c     timeicex = -1082000.  ! for hot (MIS31)
c     rco2 = 2.             ! for hot (fut)

c---

#    if defined (CO2PLIO) 
c     long-term trend, during last 5 Myr
      zt1 = -3.e6
      zt2 = -2.e6
      zw =max (0., min  (1., (timeicex - zt1) / (zt2-zt1) ))
      rco2 = (1.-zw)*(CO2PLIO/280.) + zw*1.
#    else  
      rco2 = 1.
#    endif

c---

c     laskar for dtantann, dtantaus, dtantjan, dtseas:

      call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
     *             timeicex, iutsin)

c     Increase positive dtantjan 1/2 cycle ~116 ka, decrease it ~93 ka: 

#    if defined (FUDGE116) || defined (FUDGE116P)
      if (timeice.ge.-122000. .and. timeice.le.-108000.) then
#      if defined (FUDGE116)
c       if (dtantjan.gt.-4.) dtantjan = (dtantjan + 4.)*(8.5/5.5) - 4.
        if (dtantjan.gt.-4.) dtantjan = (dtantjan + 4.)*(2.0)     - 4.
#      elif defined (FUDGE116P)
        if (dtantjan.gt.-4.) dtantjan = (dtantjan + 4.)*(3.0)     - 4.
#      endif
      endif
#    endif

#    if defined (FUDGE90)
      if (timeice.ge.-100000. .and. timeice.le.-85000.) then  
        if (dtantjan.gt.-2.) dtantjan = (dtantjan + 2.)*(3.5/6.5) - 2.
      endif
#    endif

c     call epicadeltad (dtantann, timeicex, iutsin)  

c---

c     lisiecki for ro18, sealev:

#    if defined (SPRATT)
      call spratt (ro18, sealev, timeicex, iutsin)
#    else
      call lisiecki (ro18, sealev, timeicex, iutsin)
#    endif

c     call specmap  (ro18, sealev,   timeicex, iutsin)
c     call bassinot (ro18, sealev,   timeicex, iutsin)

c~~~~~~~~~~~~~~~~~~~~~~~~~~
#    if defined (RCMMATRIX)
c~~~~~~~~~~~~~~~~~~~~~~~~~~

cbbbbbbbbbbbbbbbbbbbbbbbbbbbb
#      if defined (PARAMCLIM)   
      if (timeice.gt.0.) then
#      else
      if (.true.) then
#      endif
cbbbbbbbbbbbbbbbbbbbbbbbbbbbb

c****************************
#      if defined (RCMMATCO2) 
c****************************

c        Set weirun and rco2 for scenario given by one of RCPSPIN*, 
c        RCP26*, RCP45*, RCP85*, COP*, RAMP*

      call scenario (timeice, nyearstart, weirun, rco2)

c**********
#      else
c**********

      rco2 = 1.
      weirun = 1.
#        if defined (RCMCO2) 
      rco2 = RCMCO2
#        endif
#        if defined (RCMWEI) 
      weirun = RCMWEI
#        endif

c***********
#      endif
c***********

cbbbbbbbbbb
      endif
cbbbbbbbbbb

c~~~~~~~~~
#    endif
c~~~~~~~~~

c-----------------------------
#  elif defined (FORCEFUTRAMP)
c-----------------------------

      ro18 = 0.
      sealev = 0.
      rco2   = 1. + 6. * max (0., min (1., timeice/200.))
      weirun = 1. + 1. * max (0., min (1., timeice/200.))

      dtantann = 0.           ! modern
      dtantjan = 0.           ! modern
      dtseas = 0.             ! modern
      ecc = 0.                ! neutral orbit
      obl = 23.5 * pi/180.    ! neutral orbit
      prec = 0.               ! neutral orbit

c------------------------
#  elif defined (FORCEEO)
c------------------------

      ro18 = 0.
      sealev = 0.

c     hyst:
c     timea = -37.0e6
c     timeb = -33.0e6
c     if (timeice.lt.timeb) then
c       rco2 = 6. - 1.*((timeice-timea)/1.e6)
c     else
c       rco2 = 6. - 1.*((timeb  -timea)/1.e6)
c    *            + 1.*((timeice-timeb)/1.e6)
c     endif
c     call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
c    *             timeice, iutsin)

c     co2:
      call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
     *             timeice, iutsin)
c     or:
c     dtantann = 0.           ! modern
c     dtantjan = 0.           ! modern
c     dtseas = 0.             ! modern
c     ecc = 0.                ! neutral orbit
c     obl = 23.5 * pi/180.    ! neutral orbit
c     prec = 0.               ! neutral orbit
c
c     if CO2INTER, rco2 will be overridden by rco2inter
cc    rco2 = max (4. - 2.*(timeice-nyearstart)/1.e6, 1.)
cc    call paganico2 (rco2, timeice)
c     call pearsonco2 (rco2, timeice)
      timea = -34.e6
      if (timeice.lt.timea) then
        rco2 = 2. + 2.*min (1., (timea-timeice)/1.e6)
      else
        rco2 = 3. - 1.*cos(pi*(timeice-timea)/1.e6)
      endif

c     resol:
c     timeicex= 0.
c     call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
c    *            timeicex, iutsin)
c     rco2 = 3.0

      weirun = 1.

c---------------------------
#  elif defined (FORCE34MYR)
c---------------------------

c===
c     idealized sinusoidal:
c     zw = 0.5*(1. + cos(2.*pi*(timeice+34.e6)/1.e5))   ! 0 cold, 1 warm
c     ro18 = 1.*(1.-zw) - 0.3*zw                   ! dtatm(C) = -ro18*10
c     sealev = (-120.)*(1.-zw) + 10.*zw
c     dtantann = 0.
c     dtantjan = 0.
c     dtantaus = 0.
c     dtseas = 0. 
c     rco2 = 1.
c     for iceocean.F: 
c     dtocn34m = (-2.)*(1.-zw) + 2.*zw  
c     weirun = 1. - ro18      ! not used, prevent setting below (-1.e20)
c===

c co2 for atmos (iceclimparam):

c     start with Zhang co2:
      call zhangco2 (rco2, timeice, iutsin)

c.............................
#    if defined (PAGANI_EOCO2)
c.............................
c     For window around E-O (35.7 to 32.5 Ma), use linear regression
c     to "total data" in Pagani et al., Science, 2011, SOM Fig. S8
c     (zco2_pag). Blend into Zhang 0.5 Myr each side of window.

c     first, ensure earlier Zhang doesn't dip below start of Pagani:
      if (timeice.lt.-35.7e6) then
        rco2 = max (rco2, 930/280.)
      endif

      if (timeice.gt.-35.7e6-0.5e6 .and. timeice.lt.-32.5e6+0.5e6) then
         if (timeice.lt.-34.1e6) then
           zwp = (timeice+35.7e6)/(-34.1e6+35.7e6)
           zco2_pag = (930./280.)**(1.-zwp) + (740./280.)**zwp
         else 
           zwp = (timeice+34.1e6)/(-32.5e6+34.1e6)
           zco2_pag = (740./280.)**(1.-zwp) + (570./280.)**zwp
         endif

         if (timeice.lt.-35.7e6) then
           zw = max (0., min (1., (timeice+35.7e6+0.5e6)/0.5e6 ))
         else if (timeice.gt.-32.5e6) then
           zw = max (0., min (1., 1. - (timeice+32.5e6)/0.5e6 ))
         else
           zw = 1.
         endif
         rco2 = (1.-zw)*rco2 +  zw*zco2_pag
       endif
c.........
#    endif
c.........

c     use Zhang co2 before 5.3 Ma, Lisiecki+Raymo after 5 Ma:
      call lisiecki (zo18_lis, zsl_lis, timeice, iutsin)
      zco2_lis = ((1.-zo18_lis)*280. + zo18_lis*180.) / 280. ! from d18O
      zw = max (0., min (1., (timeice+5.3e6)/0.3e6 ))
      rco2 = (1.-zw)*rco2 +  zw*zco2_lis

c     pre-Zhang (38.62 Ma), use ramp to 1200 ppm at 50 Ma
c     (based roughly on McKay et al., Phil.Trans.Roy.Soc., 2016):
      zw = max (0., min (1., (timeice+50.e6)/(-38.62e6+50.e6) ))
      rco2 = zw*rco2 + (1.-zw)*(1200./280.)

c ro18 for atmos, influence of NH glaciations (iceclimparam):

c     use modern (ro18=0) before 2.8 Ma, L+R after 2.5 Ma:
      zw = max (0., min (1., (timeice+2.8e6)/0.3e6 ))
      ro18 = (1.-zw)*0. + zw*zo18_lis                  ! 0=modern, 1=LGM

c dtantjan,ann for atmos (iceclimparam), dtantann for ocn (below):

      call laskar (dtantann, dtantaus, dtantjan, dtseas, ecc, obl, prec,
     *             timeice, iutsin)

c dtocn34m for ocean (iceocean.F):

c===
c     use Mudelsee ro18 before 5.3 Ma, Lisiecki+Raymo after 5 Ma:
c     call mudelsee (zo18_mud, timeice, iu)
c     zw = max (0., min (1., (timeice+5.3e6)/0.3e6 ))
c     zo18 = (1.-zw)*zo18_mud + zw*zo18_lis            ! 0=modern, 1=LGM
c
c     recover absolute d18O, inverse of scaling in lisieki and mudelsee:
c     zao = 3.23 + zo18*(5.02-3.23)
c
c     deep-time Mg/Ca (DSTs) ~proportional to d18O, with ranges 
c     over 0-50 Ma of 1 to 9 deg, 3 to 0 d18O (Gasson 2012; Lear 2000,
c     2010; Dutton 2005; Shevenell 2008; Dawber+Tripati 2011):
c     if (zao.gt.3.) then 
c       zw = max (0., min (1., (zao-3.)/(5.-3.)))
c       dtocn34m = zw*(-2). + (1.-zw)*0.                 ! LGM to modern
c     else
c       zw = max (0., min (1., (zao-1.)/(3.-1.)))
c       dtocn34m = zw*(0).  + (1.-zw)*8.                 ! deep time
c     endif
c===
      zl1200 = log(1200./280.)   ! 50 Ma
      zl400  = log(400. /280.)   ! LIG
      zl280  = log(280. /280.)   ! modern
      zl180  = log(180 ./280.)   ! LGM
      zl     = log(rco2)
      if (zl.gt.zl400) then 
        zw = max (0., (zl-zl400) / (zl1200-zl400) )
        dtocn34m = zw*(8.) + (1.-zw)*2.                ! warmer than LIG
      else if (zl.gt.zl280) then 
        zw = max (0., min (1., (zl-zl280) / (zl400-zl280) ))
        dtocn34m = zw*(2.) + (1.-zw)*0.                ! modern to LIG
      else 
        zw = max (0., min (1., (zl-zl180) / (zl280-zl180) ))
        dtocn34m = zw*(0.) + (1.-zw)*(-2.)             ! LGM to modern
      endif

c     add smaller influence of orbit (annual insol):
c     dtocn34m = dtocn34m + dtantann

#if defined (DTOCN34M)
      dtocn34m = dtocn34m + DTOCN34M
#endif
#if defined (COEFFCO2OCN)
      dtocn34m = COEFFCO2OCN*dtocn34m
#endif

      weirun = 1. - ro18      ! not used, prevent setting below (-1.e20)

c===

c sealev:
      
c     totsl is glob.s.l. fall (m) due to all modeled ice, from icetab2d.
c     Use it before 2.8 Ma (rel. to +56 m with zero Ant ice, 
c     so sealev = ~0 for modern-like Ant.),
c     and L+R (proportional to d18O) after 2.5 Ma:
      zw = max (0., min (1., (timeice+2.8e6)/0.3e6 ))
      sealev = (1.-zw)*(-totsl + 56.) + zw*zsl_lis     
c------
#  else
c------

      ro18 = 0.  
      sealev = 0.  
      dtantann = 0.
      dtantjan = 0.
      dtantaus = 0.
      dtseas = 0. 
      rco2 = 1.
      weirun = 1.

c-------
#  endif
c-------

c---weirun setting---
c        Set weirun (used in iceocean, via zclim, 0:1:2 = lgm:mod:hot),
c        if not set above 

c. . . . . . . . . . . . . . . .
      if (weirun.eq.-1.e20) then
c. . . . . . . . . . . . . . . .
#  if defined (WARMAMPO)
        zwarm = WARMAMPO
#  else    
        zwarm = 1.
#  endif   
#  if defined (COLDAMPO)
        zcold = COLDAMPO
#  else    
        zcold = 1.
#  endif   

        weirun  = 1.
c    *          + max (0., dtantjan) / 3.  ! 666666 ! old
c    *          - (ro18*125/85.)           ! 666666 ! old
c   
c    *          + max (0., dtantjan) / 3.  ! 666666 ! cur (40km/40f)
c    *          - 1.*ro18                  ! 666666 ! cur (40km/40f)
c
c   *           + dtantaus / 6.            ! 666666 ! 40km/Aus/ausd
c   *           - 0.2 - 1.*ro18            ! 666666 ! 40km/Aus/ausd

c               cf. iceocean, LIU2D, pre 22 kyr:

c    *          + 0.1*(dtantjan+2.431)     ! av dtantjan 5ma = -2.431
     *          + 0.2*max(dtantjan,0.)     ! like 2009 (.25)

     *          - 1.*zcold*max(ro18,0.)
     *          - 1.*zwarm*min(ro18,0.)

     *          + 1.*(log(rco2)/log(2.))
c. . . . .
      endif
c. . . . .

c     weirun = max (0., min (2., weirun))
c---weirun setting end---

#  if defined (PITFUTCO2) || defined (PITFUTOCN)
      if (timeice.gt.0.) then
        zfac = max (0., min (1., (timeice-0.)/150. ))
#    if defined (PITFUTCO2) 
        rco2   = 1.*(1.-zfac) + PITFUTCO2 * zfac
#    endif
#    if defined (PITFUTOCN) 
        weirun = 1.*(1.-zfac) + 2. * zfac
#    endif
      endif
#  endif

#  if defined (SL_ICE5G)
c     override sealev from above with ICE5G-derived sl (last 40 kyr):
      call slice5g (zslb, timeice, iutsin)            ! only since 40 ka
      zwei = max (0., min (1., (timeice+40000.)/5000.))
      sealev = (1.-zwei)*sealev + zwei*zslb
#  endif

#  if defined (SLMODEL) 
#    if defined (SLMODELRADIAL)
c        If 2-D EISANTA s.l. radial model, reset sealev if changed above
      sealev = sealevin
#    else
c        If 2-D EISANTA s.l. non-radial model, maintain sealev at zero.
c        Initial elevations (hb, etc) have been adjusted in main prog 
c        above main time loop. (hb's passed from subr bedrock_sl are 
c        *relative* to current local sea level).
      sealev = 0.
#    endif
#  endif

# if defined (RWRAP_LIG) 
      ro18 = 0.  
      sealev = 0.
      rco2 = 1.
      weirun = 1.
      dtantann = 0.
      dtantaus = 0.
      dtantjan = 0.
      dtseas   = 0.
# endif

c======================================================
#elif defined (NHA) || defined (CARB) || defined (ORDO)
c======================================================

      ro18 = 0.
      sealev = 0.
      dtantann = 0.
      dtantjan = 0.
      dtantaus = 0.
      dtseas = 0.
      rco2 = 1.
      weirun = 1.

      sealevin = sealev ! used to recover sealev for SLMODELRADIAL below

#  if defined (LOVECLIP) || defined (NHAGLAC)
c     set past sea level:
#    if defined (SPRATT)
      call spratt (zro18, zsla, timeice, iutsin)
#    else
      call lisiecki (zro18, zsla, timeice, iutsin)
#    endif
      call slice5g (zslb, timeice, iutsin)            ! only since 40 ka
      zwei = max (0., min (1., (timeice+40000.)/5000.))
      sealev = (1.-zwei)*zsla + zwei*zslb
#  endif

#  if defined (NHAGLAC)
c     set past weirun (for iceocean):
      weirun = 1. - 1.*zro18                       ! zro18 from lisiecki

c     set past ecc,obl,prec (for fac_matrix):
      call laskar (dtantann,dtantaus,dtantjan,dtseas, ecc,obl,prec,
     *             timeice, iutsin)

c     set past rco2 (for iceclimgcm):
      call epicaco2 (rco2, timeice, iutsin)
#  endif

#  if defined (SL_ICE5G)
c     override sealev from above with ICE5G-derived sl (last 40 kyr):
      call slice5g (zslb, timeice, iutsin)            ! only since 40 ka
      zwei = max (0., min (1., (timeice+40000.)/5000.))
      sealev = (1.-zwei)*sealev + zwei*zslb
#  endif
   
#  if defined (SLMODEL) 
#    if defined (SLMODELRADIAL)
c        If 2-D EISANTA s.l. radial model, reset sealev if changed above
      sealev = sealevin
#    else    
c        If 2-D EISANTA s.l. non-radial model, maintain sealev at zero.
c        Initial elevations (hb, etc) have been adjusted in main prog 
c        above main time loop. (hb's passed from subr bedrock_sl are 
c        *relative* to current local sea level). 
      sealev = 0.
#    endif   
#  endif

#  if defined (NHAGLAC)
      rco2 = 1.     ! runf
      weirun = 1.   ! runf
      ro18 = 0.     ! runf

c     sealev = 0.   ! runf
c     totsl is glob.s.l. fall (m) due to all modeled ice, calculated
c     above before call to longterm, as in icetab2d.
c     Start with 7.5 m due to modern Greenland.
      sealev = 7.5 - totsl  ! runj
#   endif

c=====================
#elif defined (TEST2D)
c=====================

      sealev = 0.
      dtantann = 0.
      dtantjan = 0.
      dtantaus = 0.
      dtseas = 0.
      rco2 = 1.
      weirun = 1.

c========================
#elif defined (TESTMEL2D)
c========================

      sealev = 0.
      dtantann = 0.
      dtantjan = 0.
      dtantaus = 0.
      dtseas = 0.
      rco2 = 1.
      weirun = 1.

c====
#else
c====

c     sea level unchanged from initial value     
      dtantann = 0.
      dtantjan = 0.
      dtantaus = 0.
      dtseas = 0.
      rco2 = 1.
      weirun = 1.

c=====
#endif
c=====

      return
      end

c----------------------------------------------------------------------

c!!!!!!!!!!!!!!!!!!!!
#if defined (EISANTA)
c!!!!!!!!!!!!!!!!!!!!

      subroutine specmap (ro18, sealev, timeice, iu)

c        If first call, read SPECMAP d18O record for last 782 ka, 
c        and convert to ro18spec (saved).
c        Interpolate to model time (timeice), return ro18, sealev.
c          ro18 = 1 for LGM, 0 for modern.
c          sealev = -125 for LGM, 0 for modern (m).

#include <comicepath.h>

      character*200 cfilin
      parameter (nspec=782)
      dimension timespec(0:nspec), o18spec(0:nspec), ro18spec(0:nspec)
      logical firstspec
      data firstspec /.true./
      save timespec, o18spec, ro18spec, firstspec

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_specmap/specmap.017'

c        If first call, read data file
      
c------------------------
      if (firstspec) then
c------------------------
        write (6,'(a)') 'Reading Specmap file'
        open (iu, file=cfilin, status='old')
        read (iu,'(/)') 
        do i=0,nspec
          read (iu,*,err=200,end=200) timespec(i), o18spec(i)
        enddo
        close (iu)

c          Convert to sea level, and convert time to (-) years BP  

        o18mod = o18spec(0)
        o18lgm = o18spec(19)
        do i=0,nspec
          ro18spec(i) = (o18spec(i)-o18mod) / (o18lgm-o18mod) 
          timespec(i) = -1000.*timespec(i)
c         write (7,'(i6,3f20.2)') i, timespec(i), o18spec(i),ro18spec(i)
        enddo

        firstspec = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timespec(0)) then
        ro18 = ro18spec(0)
      else if (timeice.le.timespec(nspec)) then
        ro18 = ro18spec(nspec)
      else
        do i=0,nspec-1
          if (timeice.le.timespec(i) .and. timeice.ge.timespec(i+1))then
            zwei = (timeice-timespec(i+1)) / (timespec(i)-timespec(i+1))
            ro18   = zwei*ro18spec(i)   + (1.-zwei)*ro18spec(i+1)
            goto 100
          endif
        enddo
        write (6,'(a)') "Error in specmap - shouldn't get here"
        stop
  100   continue
      endif
      sealev = 0.*(1.-ro18) - 125.*ro18
      return

  200 write (6,'(a,i6)') 'Error 200 in reading Specmap input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine bassinot (ro18, sealev, timeice, iu)

c        If first call, read Bassinot(1994) d18O record for last 890 ka,
c        and convert to ro18bas (saved).
c        Interpolate to model time (timeice), return ro18 and sealev.
c          ro18 = 1 for LGM, 0 for modern.
c          sealev = -125 for LGM, 0 for modern (m).

#include <comicepath.h>

      character*200 cfilin
      parameter (nbas=443)
      dimension timebas(nbas), o18bas(nbas), ro18bas(nbas)
      logical firstbas
      data firstbas /.true./
      save timebas, o18bas, ro18bas, firstbas

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_bassinot/bassinot.final.dat'

c        If first call, read data file
      
c------------------------
      if (firstbas) then
c------------------------
        write (6,'(a)') 'Reading Bassinot file'
        open (iu, file=cfilin, status='old')
        do i=1,nbas
          read (iu,*,err=200,end=200) timebas(i), o18bas(i)
        enddo
        close (iu)

c          Convert to sea level, and convert time to (-) years BP  

        o18mod = o18bas(1)
        o18lgm = o18bas(6)
        do i=1,nbas
          ro18bas(i) = (o18bas(i)-o18mod) / (o18lgm-o18mod) 
          timebas(i) = -1000.*timebas(i)
c         write (7,'(i6,3f20.2)') i, timebas(i), o18bas(i), ro18bas(i)
        enddo

        firstbas = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timebas(1)) then
        ro18 = ro18bas(1)
      else if (timeice.le.timebas(nbas)) then
        ro18 = ro18bas(nbas)
      else
        do i=1,nbas-1
          if (timeice.le.timebas(i) .and. timeice.ge.timebas(i+1))then
            zwei = (timeice-timebas(i+1)) / (timebas(i)-timebas(i+1))
            ro18   = zwei*ro18bas(i)   + (1.-zwei)*ro18bas(i+1)
            goto 100
          endif
        enddo
        write (6,'(a)') "Error in basssinot - shouldn't get here"
        stop
  100   continue
      endif
      sealev = 0.*(1.-ro18) - 125.*ro18
      return

  200 write(6,'(a,i6)') 'Error 200 in reading Bassinot input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine epicadeltad (dtantann, timeice, iu)

c        If first call, read EPICA (Dome C) dD record for last ~800 ka, 
c        bin into 1000 yr intervals. 
c        Interpolate to model time (timeice), return dtantann.

#include <comicepath.h>

      character*200 cfilin
      parameter (nepin=5788)
      dimension timein(nepin), dtin(nepin)

      parameter (nepmax=1600)
      dimension timeepica(nepmax), dtepica(nepmax)

      logical firstepica
      data firstepica /.true./
      save timeepica, dtepica, nepica, firstepica

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_epica/edc3deuttemp2007.txt'

c        If first call, read data file, bin into uniform intervals
      
c-------------------------
      if (firstepica) then
c-------------------------
        write (6,'(a)') 'Reading Epica file'
        open (iu, file=cfilin, status='old')
        read (iu,'()') 
        do i=1,nepin
          read (iu,*,err=200,end=200) 
     *      ibag, skip1, timein(i), dum2, dtin(i)
          if (ibag.ne.i+12) then
            write (6,'(2a,2i6)') 'Error reading Epica input file.',
     *                          ' ibag, i =', ibag, i
            stop
          endif
        enddo
        close (iu)

c          Bin into 500 yr intervals, make timeepica -ve (yrs BP).
c          Also shift dt's so most recent (dtepica(1)) = 0. 

        ie = 0
        do m=1,1600
          timea = (m-1.)*500.
          timeb = m*500.

          dtbin = 0.
          nbin = 0  
          do i=1,nepin 
            if (timein(i).gt.timea .and. timein(i).le.timeb) then
              dtbin = dtbin + dtin(i)
              nbin = nbin + 1
            endif
          enddo
          if (nbin.gt.0) then
            ie = ie + 1
            timeepica(ie) = -0.5*(timea+timeb) 
            dtepica(ie) = dtbin/nbin
          endif
        enddo

        nepica = ie

        do i=1,nepica
          dtepica(i) = dtepica(i) - dtepica(1)
c         write (7,'(i6,2f20.2)') i, timeepica(i), dtepica(i)
        enddo

        firstepica = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timeepica(1)) then
        dtantann = dtepica(1)
        return
      else if (timeice.le.timeepica(nepica)) then
        dtantann = dtepica(nepica)
        return
      else
        do i=1,nepica
          if (timeice.le.timeepica(i) .and. timeice.ge.timeepica(i+1))
     *      then
            zwei = (timeice     -timeepica(i+1))
     *           / (timeepica(i)-timeepica(i+1))
            dtantann = zwei*dtepica(i) + (1.-zwei)*dtepica(i+1)
            return
          endif
        enddo
      endif
      write (6,'(a)') "Error in Epica - shouldn't get here"
      stop

  200 write (6,'(a,i6)') 'Error 200 in reading Epica input file: i=',i
      stop

      end

c!!!!!
#endif
c!!!!!

c----------------------------------------------------------------------

      subroutine epicaco2 (rco2, timeice, iu)

c        If first call, read EPICA Dome C CO2 record for last 800 kyr,
c        save in co2sav (ppm, saved).
c        Interpolate to model time (timeice), return rco2 (=co2/280).

#include <comicepath.h>

      character*200 cfilin
      parameter (nco2=1901)
      dimension timeco2(nco2), co2sav(nco2)
      logical firstco2
      data firstco2 /.true./
      save timeco2, co2sav, firstco2

c        If first call, read data file
      
c------------------------
      if (firstco2) then
c------------------------
        cfilin = cpathhome(1:lenchr(cpathhome))
     *           // 'Data_co2/antarctica2015co2composite.dat'
        write (6,'(a)') 'Reading EPICA Dome C CO2 file'

        open (iu, file=cfilin, status='old')
        read (iu,'()') 
        do i=1,nco2
          read (iu,*,err=200,end=200) timeco2(i), co2sav(i), skip
        enddo
        close (iu)

        do i=1,nco2
          timeco2(i) = -timeco2(i)
c         write (7,'(i6,2f20.2)') i, timeco2(i), co2sav(i)
        enddo

        firstco2 = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timeco2(1)) then
        rco2 = co2sav(1)
      else if (timeice.le.timeco2(nco2)) then
        rco2 = co2sav(nco2)
      else
        do i=1,nco2-1
          if (timeice.le.timeco2(i) .and. timeice.ge.timeco2(i+1))then
            zwei = (timeice-timeco2(i+1)) / (timeco2(i)-timeco2(i+1))
            rco2 = zwei*co2sav(i) + (1.-zwei)*co2sav(i+1)
            goto 100
          endif
        enddo
        write (6,'(a)') "Error in epicaco2 - shouldn't get here"
        stop
  100   continue
      endif

      rco2 = rco2 / 280.
      return

  200 write (6,'(a,i6)') 
     *  'Error 200 in reading EPICA Dome C CO2 input file: i=',i
      stop

      end

c----------------------------------------------------------------------

      subroutine slice5g (sealev, timeice, iu)

c        Set sealev to ICE5G-derived eustatic sea level over
c        last 40 kyr (from eustatic_all output, s.l. radial run slzz2)

#include <comicepath.h>

c for sl_ice5g (ice5g sea level):
      character*200 cfilin
      dimension timeus(10000), sleus(10000)
      logical firsteus 
      save firsteus, timeus, sleus, neus
      data firsteus /.true./

      if (firsteus) then
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *             '/Data_searise/eus_sl_ice5g.dat'
        open (iu, file=cfilin, status='old')
        read (iu,*)   ! one header line
        neus = 0
   30   read (iu,*,end=40) timeus(neus+1), sleus(neus+1)
        neus = neus + 1
        if (neus.gt.10000) then
          write (6,'(/a/4x,a)') 
     *      '*** Error (longterm): sea-level input file too large:',
     *      cfilin(1:lenchr(cfilin))
          stop
        endif
        goto 30
   40   close (iu)
        do m=1,neus
          sleus(m) = sleus(m) - sleus(neus)
        enddo
        firsteus = .false.
      endif

c        Linearly interpolate to model time (eus... input file can be
c        unequally spaced in time, see doeusslice5gdat_smooth).

      if (timeice.le.timeus(1)) then
        sealev = sleus(1) 
      else if (timeice.ge.timeus(neus)) then
        sealev = sleus(neus) 
      else
        do i=1,neus-1
          if (timeice.le.timeus(i+1)) then
            sealev = (   sleus(i)  *(timeus(i+1)-timeice) 
     *                 + sleus(i+1)*(timeice-timeus(i)) )
     *              /  (timeus(i+1)-timeus(i))
            goto 50
          endif
        enddo
   50   continue
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine lisiecki (ro18, sealev, timeice, iu)

c        If first call, read Lisiecki+Raymo d18O record for last 5.32Ma,
c        and convert to ro18lis (saved).
c        Interpolate to model time (timeice), return ro18, sealev.
c          ro18 = 1 for LGM, 0 for modern.
c          sealev = -125 for LGM, 0 for modern (m).

#include <comicepath.h>

      character*200 cfilin
      parameter (nlis=2114)
      dimension timelis(0:nlis), o18lis(0:nlis), ro18lis(0:nlis)
      logical firstlis
      data firstlis /.true./
      save timelis, o18lis, ro18lis, firstlis

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_lisiecki/LR04stack.txt'

c        If first call, read data file
      
c------------------------
      if (firstlis) then
c------------------------
        write (6,'(a)') 'Reading Lisiecki+Raymo file'
        open (iu, file=cfilin, status='old')
        read (iu,'(////)') 
        do i=0,nlis
          read (iu,*,err=200,end=200) timelis(i), o18lis(i), skip
        enddo
        close (iu)

c          Convert to ro18 and sea level, and convert time to (-) yr BP 

        o18mod = o18lis(0)
        o18lgm = o18lis(18)
        do i=0,nlis
          ro18lis(i) = (o18lis(i)-o18mod) / (o18lgm-o18mod) 
          timelis(i) = -1000.*timelis(i)
c         write (7,'(i6,3f20.2)') i, timelis(i), o18lis(i), ro18lis(i)
        enddo

        firstlis = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timelis(0)) then
        ro18 = ro18lis(0)
      else if (timeice.le.timelis(nlis)) then
        ro18 = ro18lis(nlis)
      else
        do i=0,nlis-1
          if (timeice.le.timelis(i) .and. timeice.ge.timelis(i+1))then
            zwei = (timeice-timelis(i+1)) / (timelis(i)-timelis(i+1))
            ro18   = zwei*ro18lis(i)   + (1.-zwei)*ro18lis(i+1)
            goto 100
          endif
        enddo
        write (6,'(a)') "Error in lisiecki - shouldn't get here"
        stop
  100   continue
      endif
      sealev = 0.*(1.-ro18) - 125.*ro18
      return

  200 write (6,'(a,i6)') 
     *  'Error 200 in reading Lisiecki+Raymo input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine spratt (ro18, sealev, timeice, iu)

c        If first call, read Spratt+Lisiecki (CP, 2016) sea level stack
c        dataset for last 798 kyrs (saved).
c        Interpolate to model time (timeice), return sealev, ro18.
c          sealev = -130 for LGM (24 ka), 0 for modern (5 to 0 ka) (m).
c          ro18 = 1 for LGM, 0 for modern, lin. prop. to sealev.

#include <comicepath.h>

      character*200 cfilin
      parameter (nspr=798)
      dimension timespr(0:nspr), slspr(0:nspr)
      logical firstspr
      data firstspr /.true./
      save timespr, slspr, slmod, sllgm

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_spratt/spratt_sealevel.txt'

c        If first call, read data file
      
c------------------------
      if (firstspr) then
c------------------------
        write (6,'(a)') 'Reading Spratt+Lisiecki file'
        open (iu, file=cfilin, status='old')
        read (iu,'()') 
        do i=0,nspr
          read (iu,*,err=200,end=200) timespr(i), slspr(i)
          if (timespr(i).lt.5.) slspr(i) = 0.        ! < 5 ka, see SL16.
        enddo
        close (iu)

c          Convert time to (-) yr BP 

        slmod = slspr(0)    ! mod
        sllgm = slspr(24)   ! lgm
        do i=0,nspr
          timespr(i) = -1000.*timespr(i)
          write (7,'(i6,2f20.2)') i, timespr(i), slspr(i)
        enddo

        firstspr = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timespr(0)) then
        sealev = slspr(0)
      else if (timeice.le.timespr(nspr)) then
        sealev = slspr(nspr)
      else
        do i=0,nspr-1
          if (timeice.le.timespr(i) .and. timeice.ge.timespr(i+1))then
            zwei = (timeice-timespr(i+1)) / (timespr(i)-timespr(i+1))
            sealev = zwei*slspr(i)   + (1.-zwei)*slspr(i+1)
            goto 100
          endif
        enddo
        write (6,'(a)') "Error in spratt - shouldn't get here"
        stop
  100   continue
      endif

c     ro18 = 1 for LGM, 0 for modern:
      ro18 = 1. - (sealev-sllgm) / (slmod-sllgm)

      return

  200 write (6,'(a,i6)') 
     *  'Error 200 in reading Spratt+Lisiecki input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine mudelsee (ro18, timeice, iu)

c        Returns ro18 for 4 to 61 Ma, from Mudelsee et al. (2014, 
c        Rev. Geophys., supp inf. table, curve for "high latitudes").

#include <comicepath.h>

      parameter (nmud = 57001)
      dimension tmud(nmud), omud(nmud)

      logical firstmud
      data firstmud /.true./
      save tmud, omud, firstmud

      character cfilin*200, cskip*8

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 
     *  'Data_d18o/Mudelsee-etal-2014-RoG-stack-data-high-latitudes.txt'

c        If first call, read data file, reversing time order 
      
c-----------------------
      if (firstmud) then
c-----------------------
        write (6,'(a)') 'Reading Mudelsee d18O file'
        open (iu, file=cfilin, status='old')
        do i=1,11
          read (iu,*)  ! skip 11 header lines
        enddo
        do i=1,nmud
          j = nmud - i + 1
          read (iu,*,err=200,end=200) tmud(j), omud(j)
        enddo
        close (iu)

c          Convert time to (-) yr BP, and do18 to ro18, using 
c          Lisiecki and Raymo d18O values (3.23 modern, 5.02 LGM) 
c          as in subr lisiecki

        do i=1,nmud
          tmud(i) = -1.e6 * tmud(i)                    ! Ma to -yr
          zmud = omud(i)
          omud(i) = (omud(i)-3.23) / (5.02-3.23)       ! 0=modern, 1=LGM
c         write (766,*) i, tmud(i), zmud, omud(i)
        enddo
        firstmud = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.le.tmud(1)) then
        ro18 = omud(1)
        return
      else if (timeice.ge.tmud(nmud)) then
        ro18 = omud(nmud)
        return
      else
        do i=1,nmud-1
          if (timeice.ge.tmud(i) .and. timeice.le.tmud(i+1)) then
            zwei = (timeice-tmud(i)) / (tmud(i+1)-tmud(i))
            ro18 = (1.-zwei)*omud(i) + zwei*omud(i+1)
            return
          endif
        enddo
      endif

      write (6,'(a)') "Error in mudelsee - shouldn't get here"
      stop

  200 write (6,'(a)') "Error in reading Mudelsee d18O file"
      stop

      end

c-----------------------------------------------------------------------

      subroutine paganico2 (rco2, timeice)

c        Interpolates within a small range of time for atmospheric
c        co2 level, based on Pagani's new data (vertices from L. Kump)

      parameter (nco2 = 6)
      dimension  tco2(nco2), pco2(nco2)
      data tco2 / -34.5, -34.0, -33.5, -33.0, -32.5, -32.0/
      data pco2 / 1120., 1120.,  560.,  560.,  625.,  625./
      logical firstco2
      data firstco2 /.true./
      save tco2, pco2, firstco2

      if (firstco2) then
        do i=1,nco2
          tco2(i) = 1.e6* tco2(i)    ! Myr to yr
          pco2(i) = pco2(i) / 280.   ! ppmv to PAL
        enddo
        firstco2 = .false.
      endif

c        Interpolate to model time

      if (timeice.le.tco2(1)) then
        rco2 = pco2(1)
        return
      else if (timeice.ge.tco2(nco2)) then
        rco2 = pco2(nco2)
        return
      else
        do i=1,nco2-1
          if (timeice.ge.tco2(i) .and. timeice.le.tco2(i+1)) then
            zwei = (timeice-tco2(i)) / (tco2(i+1)-tco2(i))
            rco2 = (1.-zwei)*pco2(i) + zwei*pco2(i+1)
            return
          endif
        enddo
      endif

      write (6,'(a)') "Error in paganico2 - shouldn't get here"
      stop

      end

c-----------------------------------------------------------------------

      subroutine pearsonco2 (rco2, timeice)

c        Interpolates within a small range of time for atmospheric
c        co2 level, based on Pearson's new data (Nature, 2009,
c        variable CO3_2-)

c nominal:
c     parameter (nco2 = 9)
c     dimension  tco2(nco2), pco2(nco2), dpco2(nco2)
c     data tco2 /-34.63, -34.42, -34.21, -33.96, -33.59, 
c    *           -33.47, -33.33, -33.17, -32.98/ 
c     data pco2 /  876.,  1046.,  1074.,   959.,   763.,
c    *            1160.,  1105.,   786.,   624./ 
c     lower envelope departure ("err +" in xls file for variable CO3_2-)
c     dimension  dpco2(nco2)
c     data dpco2/  187.,   217.,   155.,   134.,    94.,
c    *             162.,   148.,   100.,    73./ 
c imposed dip at 33.55 Ma:
      parameter (nco2 = 10)
      dimension  tco2(nco2), pco2(nco2)
      data tco2 /-34.63, -34.42, -34.21, -33.96, -33.70, -33.55, 
     *           -33.47, -33.33, -33.17, -32.98/ 
c     data pco2 /  876.,  1046.,  1074.,   959.,   800.,   700.,  ! a
c     data pco2 /  876.,  1046.,  1074.,   959.,   800.,   625.,  ! b
      data pco2 /  876.,  1046.,  1074.,   959.,   700.,   600., 
     *            1160.,  1105.,   786.,   624./ 

      logical firstco2
      data firstco2 /.true./
      save tco2, pco2, firstco2

      if (firstco2) then
        do i=1,nco2
          tco2(i) = 1.e6* tco2(i)               ! Myr to yr
          pco2(i) = (pco2(i)         ) / 280.   ! ppmv to PAL
c         pco2(i) = (pco2(i)-dpco2(i)) / 280.   ! ppmv to PAL
        enddo
        firstco2 = .false.
      endif

c        Interpolate to model time

      if (timeice.le.tco2(1)) then
        rco2 = pco2(1)
        return
      else if (timeice.ge.tco2(nco2)) then
        rco2 = pco2(nco2)
        return
      else
        do i=1,nco2-1
          if (timeice.ge.tco2(i) .and. timeice.le.tco2(i+1)) then
            zwei = (timeice-tco2(i)) / (tco2(i+1)-tco2(i))
            rco2 = (1.-zwei)*pco2(i) + zwei*pco2(i+1)
            return
          endif
        enddo
      endif

      write (6,'(a)') "Error in pearsonco2 - shouldn't get here"
      stop

      end

c-----------------------------------------------------------------------

      subroutine zhangco2 (rco2, timeice, iu)

c        Returns rco2 for last ~40 Myr, from Zhang et al (2013, Phil.
c        Trans. Roy. Soc., supp inf. table).

#include <comicepath.h>

      parameter (nco2 = 87)
      dimension  tco2(nco2), pco2(nco2)

      logical firstco2
      data firstco2 /.true./
      save tco2, pco2, firstco2

      character cfilin*200, cskip*8

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_co2/zhang_2013_co2data.prn'

c        If first call, read data file, reversing time order 
      
c-----------------------
      if (firstco2) then
c-----------------------
        write (6,'(a)') 'Reading Zhang co2 file'
        open (iu, file=cfilin, status='old')
        read (iu,*)
        do i=1,nco2
          j = nco2 - i + 1
          read (iu,*,err=200,end=200) 
     *      (cskip,m=1,4), tco2(j), (cskip,m=1,7), pco2(j)
        enddo
        close (iu)

c          Convert time to (-) yr BP, co2 to PAL

        do i=1,nco2
          tco2(i) = -1.e6 * tco2(i)            ! Ma to -yr
          pco2(i) = pco2(i) / 280.             ! ppm to PAL
c         write (765,*) i, tco2(i), pco2(i)
        enddo
        firstco2 = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.le.tco2(1)) then
        rco2 = pco2(1)
        return
      else if (timeice.ge.tco2(nco2)) then
        rco2 = pco2(nco2)
        return
      else
        do i=1,nco2-1
          if (timeice.ge.tco2(i) .and. timeice.le.tco2(i+1)) then
            zwei = (timeice-tco2(i)) / (tco2(i+1)-tco2(i))
            rco2 = (1.-zwei)*pco2(i) + zwei*pco2(i+1)
            return
          endif
        enddo
      endif

      write (6,'(a)') "Error in zhangco2 - shouldn't get here"
      stop

  200 write (6,'(a)') "Error in reading Zhang co2 file"
      stop

      end

c-----------------------------------------------------------------------

      subroutine laskar (dtantann, dtantaus, dtantjan, dtseas, 
     *                   ecc, obl, prec, timeice, iu)

c        If first call, read Laskar (2004) file for orbital elements
c        from 51 Ma to present, and calculates annual, spring and Jan
c        insolation for latitude alatlask (can be -D..., default = -80).
c        Interpolate to model time (timeice), return dtant[ann,aus,jan] 
c        and dtseas.

#include <comicepath.h>

      character*200 cfilin
#if defined (FORCEEO) || defined (SIMPEO) || defined (FORCE34MYR)
      parameter (nlask=51000) ! Laskar file has last 51 Myr
#else
      parameter (nlask=5500)  ! just last 5.5 Myr
#endif
      dimension timelask(0:nlask), 
     *           ecclask(0:nlask), obllask(0:nlask), preclask(0:nlask)
      logical firstlaskread, firstlask
      data firstlaskread, firstlask /.true., .true./
      save firstlaskread, firstlask, 
     *     timelask, ecclask,  obllask, preclask, 
     *     fluxann0, fluxaus0, fluxjan0, fluxjul0
      parameter (pi = 3.14159265358979)
      dimension cosq(1), fraq(1), cosq24(1), fraq24(1), alatq(1)
      parameter (solcon = 1367.)
      parameter (nflux=365)
#if defined (SIMPEO)
      parameter (dtdf_ann  = 0.07)        ! d[~sfc temperature]/d[insol]
      parameter (dtdf_aus  = 0.07)        ! d[~sfc temperature]/d[insol]
      parameter (dtdf_seas = 0.07)        ! d[~sfc temperature]/d[insol]
#else

#  if defined (QAMPANN)
      parameter (dtdf_ann  = QAMPANN)     ! d[~sfc temperature]/d[insol]
#  else
      parameter (dtdf_ann  = 0.1)         ! d[~sfc temperature]/d[insol]
#  endif

#  if defined (QAMPAUS)
      parameter (dtdf_aus  = QAMPAUS)     ! d[~sfc temperature]/d[insol]
#  else
      parameter (dtdf_aus  = 0.1)         ! d[~sfc temperature]/d[insol]
#  endif

#  if defined (QAMPSEAS)
      parameter (dtdf_seas = QAMPSEAS)    ! d[~sfc temperature]/d[insol]
#  else
      parameter (dtdf_seas = 0.1)         ! d[~sfc temperature]/d[insol]
#  endif

#endif

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_laskar/INSOLN.LA2004.BTL.ASC'

c        If first call, read data file
      
c----------------------------
      if (firstlaskread) then
c----------------------------
        write (6,'(a)') 'Reading Laskar orbit file'
        open (iu, file=cfilin, status='old')

        do i=0,nlask
          read (iu,*,err=200,end=200) 
     *      timelask(i), ecclask(i), obllask(i), preclask(i)

c         1000's yrs to yrs (-ve past):
          timelask(i) = 1.e3*timelask(i)  

c         convert laskar prec (NHAE->peri) to Genesis (peri->NHVE):
          preclask(i) = mod (pi - preclask(i) + 2.*pi, 2.*pi)
        enddo

        close (iu)
        firstlaskread = .false.
c----------
      endif
c----------

c        Interpolate to model time

      if (timeice.ge.timelask(0)) then
        ecc_cur  = ecclask(0)
        obl_cur  = obllask(0)
        prec_cur = preclask(0)
        goto 100
      else if (timeice.le.timelask(nlask)) then
        ecc_cur  = ecclask(nlask)
        obl_cur  = obllask(nlask)
        prec_cur = preclask(nlask)
        goto 100
      else
        do i=0,nlask-1
          if (timeice.le.timelask(i) .and. timeice.ge.timelask(i+1))then
            zwei = (timeice-timelask(i+1)) / (timelask(i)-timelask(i+1))
            ecc_cur  = zwei*ecclask(i) + (1.-zwei)*ecclask(i+1)
            obl_cur  = zwei*obllask(i) + (1.-zwei)*obllask(i+1)

            if (preclask(i+1)-preclask(i).lt.-pi) then
              zpm = preclask(i)
              zpp = preclask(i+1) + 2.*pi
            else if (preclask(i+1)-preclask(i).gt.pi) then
              zpm = preclask(i) + 2.*pi
              zpp = preclask(i+1)
            else
              zpm = preclask(i)
              zpp = preclask(i+1)
            endif
            if (abs(zpp-zpm).gt.pi) then
              write (6,'(/2a,2f8.2, 3x, 2f8.2, 3x, 2i6)') 
     *           ' Error (laskar):',
     *           ' preclask(i:i+1), zp[m,p], i,i+1 =', 
     *           preclask(i), preclask(i+1), zpm,zpp, i,i+1
              stop
            endif
            prec_cur = zwei*zpm + (1.-zwei)*zpp
            prec_cur = mod (prec_cur + 2.*pi, 2.*pi)

            goto 100
          endif
        enddo
      endif
      write (6,'(a)') "Error in laskar - shouldn't get here"
      stop
  100 continue

c     ecc_cur = .05 
c     prec_cur = mod ( 0.5*pi - (2.*pi*timeice/40.e3) + 2.*pi, 2.*pi )
c     obl_cur = (  0.5*(22.2+24.5) 
c    *           + 0.5*(24.5-22.2)*cos(2.*pi*timeice/40.e3) ) * pi/180.
c
c     ecc_cur = 0.
c     prec_cur = 0.
c     obl_cur = 23.5*pi/180.

#if defined ORBPRIME
c     ecc,prec,obl kyr periods with prime numbers:
      ecc_cur = .025 + .025 * cos(2.*pi*timeice/101.e3)
      prec_cur = mod (-2.*pi*timeice/19.e3, 2.*pi)  
      if (prec_cur.lt.0.) prec_cur = prec_cur + 2.*pi 
      obl_cur = (  0.5*(22.2+24.5) 
     *           + 0.5*(24.5-22.2)*cos(2.*pi*timeice/41.e3) ) * pi/180. 
#endif

c        If first call, do flux calculations (zencal) twice,
c        first time for modern to store flux[ann,jan,jul]0.

      if (firstlask) then
        nloopf = 2
        firstlask = .false.
      else
        nloopf = 1
      endif

c~~~~~~~~~~~~~~~~~~~~~~~
      do iloopf=1,nloopf
c~~~~~~~~~~~~~~~~~~~~~~~

        if (nloopf.eq.2 .and. iloopf.eq.1) then
          ecc  = ecclask(0)
          obl  = obllask(0)
          prec = preclask(0)
          time = 0.
        else
          ecc  = ecc_cur
          obl  = obl_cur
          prec = prec_cur
          time = timeice
        endif

        zs = sin (0.5*prec) / sqrt ((1.+ecc)/(1.-ecc))
        zc = cos (0.5*prec)
        ze = 2. * atan2 (zs,zc)
        vern = ze - ecc * sin(ze)

c          Calculate daily mean insols at latitude alatlask,
c          and save averages for annual, January and July

#if defined (ALATLASK)
        alatlask = ALATLASK
#else
        alatlask = -80.
#endif

        fluxann = 0.
        fluxaus = 0.
        fluxjan = 0.
        fluxjul = 0.
        nann = 0
        naus = 0
        njan = 0
        njul = 0

        isecdy = 86400/2
        dtq    = 86400.
        alatq(1) = alatlask*pi/180.
        do m=1,nflux
          isecyr = nint (86400.*365.*(m-.5)/nflux)
          call zencal (nint(time), isecyr, isecdy, dtq,
     *                 ecc, obl, prec, 
     *                 vern, dist, eccf,
     *                 cosq, fraq, cosq24, fraq24, alatq, 1)
          zflux = cosq(1)*fraq(1)*eccf*solcon

          fluxann = fluxann + zflux 
          nann = nann + 1

c         austral spring, Aug 21 to Nov 20
c         (Stott et al. Science 2007; Timmermann et al. J.Clim. 2009):
          if (isecyr.gt.86400*232 .and. isecyr.lt.86400*324) then
            fluxaus = fluxaus + zflux
            naus = naus + 1
          endif

c         either calendar January:
c         if (isecyr.lt.86400*31) then  
c         or Dec 21 to Jan 20, as in Laskar program:
          if (isecyr.lt.86400*20 .or. isecyr.gt.86400*(365-11)) then
            fluxjan = fluxjan + zflux
            njan = njan + 1
          endif

c         either calendar July:
c         if (isecyr.gt.86400*181 .and. isecyr.lt.86400*212) then  
c         or Jun 21 to Jul 20:
          if (isecyr.gt.86400*171 .and. isecyr.lt.86400*202) then
            fluxjul = fluxjul + zflux
            njul = njul + 1
          endif
        enddo
        fluxann = fluxann / nann
        fluxaus = fluxaus / naus
        fluxjan = fluxjan / njan
        fluxjul = fluxjul / njul

        if (nloopf.eq.2 .and. iloopf.eq.1) then
          fluxann0 = fluxann
          fluxaus0 = fluxaus
          fluxjan0 = fluxjan
          fluxjul0 = fluxjul
        endif

c~~~~~~~~~~
      enddo
c~~~~~~~~~~

c        Change to (differences from modern insolation)*(dT/dinsol)

      dtantann = ( fluxann - fluxann0 ) * dtdf_ann
      dtantaus = ( fluxaus - fluxaus0 ) * dtdf_aus
      dtantjan = ( fluxjan - fluxjan0 ) * dtdf_seas
      dtseas   = ( (fluxjan-fluxjul) - (fluxjan0-fluxjul0) ) * dtdf_seas

c     write (7,'(f10.1, 3f10.6, 8f10.2)') 
c    *  time*1.e-3, 
c    *  ecc, obl*180./pi, prec*180./pi,
c    *  fluxann, fluxaus, fluxjan, fluxjul, 
c    *  (fluxann-fluxann0)*dtdf_ann, 
c    *  (fluxaus-fluxaus0)*dtdf_aus, 
c    *  (fluxjan-fluxjan0)*dtdf_seas,
c    *  (fluxjul-fluxjul0)*dtdf_seas

      return

  200 write (6,'(a,i6)') 'Error 200 in reading Laskar input file: i=',i
      stop

      end

c-----------------------------------------------------------------------

      subroutine writemark (clabin, timeice, dtimeice1, dtimeice2)

c     Write a marker line to all ascii output files, so that
c     a post-processing program can delete repeated output
c     if blowups and re-running of stages occur. 
c     clabin is either:
c       'INIT' to signal first timestep (restart file written),
c       'REACHED' to signal successful completion of a stage,
c       'REWINDING' to signal blowup and re-running current stage.
c
c     The resulting output files should have only REACHED lines, 
c     with possible REWINDING lines whose times (timeice) correspond
c     to the previous REACHED line, i.e.,
c       MARK INIT A 
c       MARK REACHED B
c       MARK REWINDING B
c       MARK REWINDING B
c       MARK REWINDING B
c       MARK REACHED C
c       etc.
c
c     Also written to stdout (ioterm) as a log. 

#include <comicegrid.h>

      character*(*) clabin

      character*80 cmark

      dimension iumark(50)
      logical first
      data first /.true./
      save iumark, nmark, first

c        If first call, assemble list of ascii output file unit numbers
      
c--------------------
      if (first) then
c--------------------
        nmark = 6
        iumark(1) = ioterm
        iumark(2) = iuout2d
        iumark(3) = iuout1d
        iumark(4) = iuplot1d
        iumark(5) = iutab
        iumark(6) = iutabbud

#if defined (LINEM) && defined (RESOL6) 
        nmark = nmark + 1
        iumark(nmark) = iuout1dzoom
#endif

#if defined (EISANTA)
        nmark = nmark + 1
        iumark(nmark) = iupts
        nmark = nmark + 1
        iumark(nmark) = iutabwais
#endif

#if defined (SEDIMENT)
        nmark = nmark + 1
        iumark(nmark) = iusedbud
#endif

#if defined (TABTROUGH)
        nmark = nmark + 1
        iumark(nmark) = iutabtro
#endif

#if defined (MELANGE)
        nmark = nmark + 1
        iumark(nmark) = iutabmel
#endif

        first = .false.
c----------
      endif
c----------

      cmark = 'MARK ' // clabin 

      do m = 1,nmark
        write (iumark(m),'(a, a,f15.2, a,f8.3, a,f8.3)') 
     *    cmark(1:15), '  time=',timeice,
     *    '  dtimeice(pre)=',dtimeice1, '  dtimeice(post)=',dtimeice2
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine checkblow (arr, nx1,nx2, ny1,ny2, nz1,nz2,  
     *                      emin, emax, anorm,
     *                      cvar, timeice, ifblow) 

c        Check for blowup (extreme ice value or NaN), setting 
c        ifblow = .true. to re-run current stage with reduced timestep

#include <comicegrid.h>

      dimension arr(nx1:nx2, ny1:ny2, nz1:nz2)
      character*(*) cvar
      logical ifblow

c     znaneg = sqrt(-1.)
c     znaneg2 =  1./0.
c     znaneg3 = -1./0.
      znaneg =  1.e50
      znaneg2 = 1.e50
      znaneg3 = 1.e50

      zmin =  1.e20
      zmax = -1.e20
      inan = 0
      do k=nz1,nz2
        do j=ny1,ny2
          do i=nx1,nx2

            if (arr(i,j,k)-anorm.lt.zmin) then
              imin = i
              jmin = j
              kmin = k
              zmin = arr(i,j,k)-anorm
            endif
            if (arr(i,j,k)-anorm.gt.zmax) then
              imax = i
              jmax = j
              kmax = k
              zmax = arr(i,j,k)-anorm
            endif

            if (inan.eq.0) then
c             zv = arr(i,j,k)
c             if ( (zv.eq.znaneg.or.zv.eq.znaneg2.or.zv.eq.znaneg3) .or.
c    *             ( .not. (zv.le.0..or.zv.gt.0.)) .or.
c    *             (zv.gt.1.e50 .or. zv.lt.-1.e50)
c    *           ) then 
              if (isnan(arr(i,j,k))) then
                inan = i
                jnan = j
                knan = k
              endif
            endif

          enddo
        enddo
      enddo

      if (zmin.lt.emin .or. zmax.gt.emax) then
        iminh = max (1, min (nx, imin))
        jminh = max (1, min (ny, jmin))
        imaxh = max (1, min (nx, imax))
        jmaxh = max (1, min (ny, jmax))
        write (ioterm,
     *    '(/3a,2e12.4 / a,f12.1 / 4x,a / (a4,f10.1,3i5,2f9.2,2f8.2))')
     *    ' *** Error: ', cvar,' outside range', emin, emax,
     *    '     time(yr)=', timeice,
     *    '       val    i    j    k    x(km)    y(km)     lon     lat',
 
     *    'min:', zmin, imin, jmin, kmin,
     *    xh(iminh,jminh)*.001, yh(iminh,jminh)*.001,
     *    alond(iminh,jminh), alatd(iminh,jminh),  

     *    'max:', zmax, imax, jmax, kmax,
     *    xh(imaxh,jmaxh)*.001, yh(imaxh,jmaxh)*.001,
     *    alond(imaxh,jmaxh), alatd(imaxh,jmaxh)  
        ifblow = .true.
      endif

      if (inan.ne.0) then
        inanh = max (1, min (nx, inan))
        jnanh = max (1, min (ny, jnan))
        write (ioterm,
     *    '(/3a / a,f12.1 / a / (3i5,2f9.2,2f8.2))')
     *    ' *** Error: ', cvar,' has NaN(s)',
     *    '     time(yr)=', timeice,
     *    '     i    j    k    x(km)    y(km)     lon     lat',
     *    inan, jnan, knan,
     *    xh(inanh,jnanh)*.001, yh(inanh,jnanh)*.001,
     *    alond(inanh,jnanh), alatd(inanh,jnanh)
        ifblow = .true. 
      endif

      return
      end

c-----------------------------------------------------------------------

c@@@@@@@@@@@@@@@@@@@@@@@@
#if defined (LAGRANGEO18)
c@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine lagrangeo18 (h, maskwater, u, v, w, timeice)

c        Traces all 3-D ice cells back to surface, using Lagrangian
c        stepping and 3-D linear interp of u,v,w. Assigns d18O of
c        precip to emerged surface pt, read from GCM driving file 
c        read by geto18p and bilin interped to sfc pt locn. Weights by 
c        cell volume, sums, reports mean d18O of ice sheet 
c        (to stdout, diagnostic maps to fort.192)
c
c        Lagrangian tracking iteratively steps each point backwards
c        in time, with an adjustable timestep (dtlag) depending
c        on local velocities (so v*dt ~ 0.05*grid size).
c
c        Vertical coord is ice sheet's zeta (0 to 1, sfc to base),
c        and vert veloc is w (=d(zeta)/dt).

c        Same Lagrangian tracking is also used in path2surf (icehydrol)
c        if ENHINV, to get distance to nearest dome.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny),               maskwater(nx,ny),
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  w(nx,ny,0:nlevp)

      dimension
     *  o18precip(nx,ny),    o18all(nx,ny,nlev),
     *  xsurf(nx,ny,nlev),   ysurf(nx,ny,nlev),
     *  statifm(nx,ny,nlev),    
     *  statu(nx,ny,nlev),   statv(nx,ny,nlev),
     *  stattim(nx,ny,nlev), statnum(nx,ny,nlev), 
     *  statvol(nx,ny),
     *  worka(nx,ny),        workb(nx,ny)

      save o18precip

      parameter (nzfine = 100000)
      dimension izfinea(nzfine), izfineb(nzfine), wzfinea(nzfine)
      save izfinea, izfineb, wzfinea

      parameter (niterlag = 1000000)

      logical first
      data first /.true./
      save first

c        If first call, set indices from very fine uniform vertical
c        grid (0 to 1) to model vertical grid (saved). Only use
c        vertical levels between k=1 and nlev (not 0 and nlevp).    
c        Also read in d18O of precip map from GCM driving file.

c---------------------
      if (first) then
c---------------------
        do m=1,nzfine
          zf = (m-0.5)/nzfine 
          if (zf.le.zeta(1)) then
            izfinea(m) = 1
            izfineb(m) = 1
            wzfinea(m) = 1.
          else if (zf.ge.zeta(nlev)) then
            izfinea(m) = nlev
            izfineb(m) = nlev
            wzfinea(m) = 1.
          else
            do k=2,nlev
              if (zf.le.zeta(k)) then
                izfinea(m) = k-1
                izfineb(m) = k
                wzfinea(m) = (zeta(k)-zf) / (zeta(k)-zeta(k-1))
                goto 100
              endif
            enddo
            write (6,'(/a)') '*** Error: lagrangeo18 (init)'
            stop
          endif
  100     continue
        enddo

c         Read d18O of precip from GCM ice driving file, made by 
c         debplotgcmnet ('O18', ipole=-3, data_icedrive_O18)

        call geto18p (o18precip, 192, timeice)

        first = .false.
c----------
      endif
c----------

      call resetr (o18all, nx*ny*nlev, 1.e20) 
      call zero (xsurf,   nx*ny*nlev) 
      call zero (ysurf,   nx*ny*nlev) 
      call zero (statu,   nx*ny*nlev)
      call zero (statv,   nx*ny*nlev)
      call zero (statifm, nx*ny*nlev)
      call zero (stattim, nx*ny*nlev)
      call zero (statnum, nx*ny*nlev)
      call zero (statvol, nx*ny)

c        Main loop over each 3-D cell with ice.
c        Only do for grounded ice (not sure if backward Lagrangian
c        tracking will work with fast shelf velocities). 

c>>>>>>>>>>>>>>
      do j=1,ny
        do i=1,nx
c         if (h(i,j).gt.0.) then                            ! all ice
          if (h(i,j).gt.0. .and. maskwater(i,j).eq.0) then  ! grounded
            do k=1,nlev
c>>>>>>>>>>>>>>
              xp = xh(i,j)
              yp = yh(i,j)
              zp = zeta(k)

c                Loop over each backward veloc*dt segment for this path.
c                Current location is (xp,yp,zp), current velocities
c                interpolated to this point are (up,vp,wp). When reach
c                surface of ice, record location xsurf,ysurf, skip out. 

c++++++++++++++++++++++++++++++++++++
              do iterlag = 1,niterlag
c++++++++++++++++++++++++++++++++++++

c                Interpolate velocities to current location

              xprev = xp
              yprev = yp
              zprev = zp

              zxp = ( xp - xoffa + 0.5*nx*dx0 ) / dx0
              zyp = ( yp - yoffa + 0.5*ny*dy0 ) / dy0

c---
              iua = int(zxp)
              iub = iua + 1
              wui = iub - zxp
              iua = max (0, min (nxp, iua))
              iub = max (0, min (nxp, iub))

              jua = nint(zyp)
              jub = jua + 1
              wuj = jub - (zyp+0.5)
              jua = max (1, min (ny, jua))
              jub = max (1, min (ny, jub))
c---
              iva = nint(zxp)
              ivb = iva + 1
              wvi = ivb - (zxp+0.5)
              iva = max (1, min (nx, iva))
              ivb = max (1, min (nx, ivb))

              jva = int(zyp)
              jvb = jva + 1
              wvj = jvb - zyp
              jva = max (0, min (nyp, jva))
              jvb = max (0, min (nyp, jvb))
c---
              iwa = iva
              iwb = ivb
              wwi = wvi

              jwa = jua
              jwb = jub
              wwj = wuj
c---
              m = nint(zp*nzfine + 0.5)
              m = max (1, min (nzfine, m))
              ka = izfinea(m)
              kb = izfineb(m)
              wk = wzfinea(m)
c---

c             if (h(iwa,jwa).eq.0. .or. h(iwb,jwa).eq.0. .or.
c    *            h(iwa,jwb).eq.0. .or. h(iwb,jwb).eq.0.) then
c               write (6,"(/a,i6,'  iwa,iwb,jwa,jwb=',4i6/
c    *                     /'h= ',4f12.5)")
c    *            'Warning (lagrangeo18 h=0): iterlag=', iterlag,
c    *            iwa,iwb,jwa,jwb,
c    *            h(iwa,jwa), h(iwb,jwa), h(iwb,jwa), h(iwb,jwb)
c             endif

c===
              up =       wui *    wuj *    wk  * u(iua,jua,ka)
     *             + (1.-wui)*    wuj *    wk  * u(iub,jua,ka)
     *             +     wui *(1.-wuj)*    wk  * u(iua,jub,ka)
     *             + (1.-wui)*(1.-wuj)*    wk  * u(iub,jub,ka)

     *             +     wui *    wuj *(1.-wk) * u(iua,jua,kb)
     *             + (1.-wui)*    wuj *(1.-wk) * u(iub,jua,kb)
     *             +     wui *(1.-wuj)*(1.-wk) * u(iua,jub,kb)
     *             + (1.-wui)*(1.-wuj)*(1.-wk) * u(iub,jub,kb)
c===
              vp =       wvi *    wvj *    wk  * v(iva,jva,ka)
     *             + (1.-wvi)*    wvj *    wk  * v(ivb,jva,ka)
     *             +     wvi *(1.-wvj)*    wk  * v(iva,jvb,ka)
     *             + (1.-wvi)*(1.-wvj)*    wk  * v(ivb,jvb,ka)

     *             +     wvi *    wvj *(1.-wk) * v(iva,jva,kb)
     *             + (1.-wvi)*    wvj *(1.-wk) * v(ivb,jva,kb)
     *             +     wvi *(1.-wvj)*(1.-wk) * v(iva,jvb,kb)
     *             + (1.-wvi)*(1.-wvj)*(1.-wk) * v(ivb,jvb,kb)
c===
              wp =       wwi *    wwj *    wk  * w(iwa,jwa,ka)
     *             + (1.-wwi)*    wwj *    wk  * w(iwb,jwa,ka)
     *             +     wwi *(1.-wwj)*    wk  * w(iwa,jwb,ka)
     *             + (1.-wwi)*(1.-wwj)*    wk  * w(iwb,jwb,ka)

     *             +     wwi *    wwj *(1.-wk) * w(iwa,jwa,kb)
     *             + (1.-wwi)*    wwj *(1.-wk) * w(iwb,jwa,kb)
     *             +     wwi *(1.-wwj)*(1.-wk) * w(iwa,jwb,kb)
     *             + (1.-wwi)*(1.-wwj)*(1.-wk) * w(iwb,jwb,kb)
c===

c                Arbitrarily avoid very small velocities

              up = sign ( max(  .01, abs(up)), up )
              vp = sign ( max(  .01, abs(vp)), vp )
              wp = sign ( max(1.e-8, abs(wp)), wp )

c                Set current timestep (dtlag)

              zdtx = 0.05 * dx0       / max (1.e-4, abs(up))
              zdty = 0.05 * dy0       / max (1.e-4, abs(vp))
              zdtz = 0.05 * (1./nlev) / max (1.e-6, abs(wp))
              dtlag = max (0.1, min (zdtx, zdty, zdtz))

c             write (198,'(4i8,f15.5)') iterlag, i,j,k,dtlag

c               Do backward Lagrangian step (-veloc*dt)

              xp = xp - dtlag*up
              yp = yp - dtlag*vp
              zp = zp - dtlag*wp

              stattim(i,j,k) = stattim(i,j,k) + dtlag
              statnum(i,j,k) = statnum(i,j,k) + 1.
   
c                If gone out of bottom of ice sheet (associated with
c                large -ve basefrml), stop tracking, skip out (to 200) 
   
              if (zp.gt.1.) then
c               write (6,"(/a,i6,'  xp,yp,zp=',2f10.2,f10.5)")
c    *            'Warning (lagrangeo18 zp>1): iterlag=', iterlag,
c    *            xp*.001, yp*.001, zp
                xsurf(i,j,k) = xp
                ysurf(i,j,k) = yp
                statu(i,j,k) = up
                statv(i,j,k) = vp
                statifm(i,j,k) = -1.
                go to 200
              endif   

c                If gone out of "sides" of ice sheet (kind of
c                an error), stop tracking this path, skip out (to 200)

              ip = nint ((xp - xoffa + 0.5*nx*dx0)/dx0 + 0.5)
              jp = nint ((yp - yoffa + 0.5*ny*dy0)/dy0 + 0.5)
              ip = max (1, min(nx,ip))
              jp = max (1, min(ny,jp))
              if (h(ip,jp).eq.0.) then
                xsurf(i,j,k) = xp
                ysurf(i,j,k) = yp
                statu(i,j,k) = up
                statv(i,j,k) = vp
                statifm(i,j,k) = -2.
c               write (191,"(a/a,3i6,i8/a,2f9.0,f8.4,3f10.5/a,4f10.3)")
c    *          'Warning (lagrangeo18 -1):',
c    *          '  i,j,k,iterlag=', i,j,k, iterlag-1,
c    *          '  xp,yp,zp,up,vp,wp=', xp*.001,yp*.001,zp, up,vp,wp,
c    *          '  h[1-4]=', h(iwa,jwa),h(iwb,jwa),h(iwb,jwa),h(iwb,jwb)
                go to 200
              endif

c                If gone out of top surface of ice sheet, interp
c                back in this timestep for more exact surface location,
c                skip out.

              if (zp.le.0.) then
                wei = (-zp) / (-zp + zprev)  
                xsurf(i,j,k) = (1.-wei)*xp + wei*xprev 
                ysurf(i,j,k) = (1.-wei)*yp + wei*yprev 
                statu(i,j,k) = up
                statv(i,j,k) = vp
                statifm(i,j,k) = 1.
c               write(190,"(a/a,3i6,i8 /a,2i8,f8.4, 3e14.5 /a,4f10.3)")
c    *          'lagrangeo18 1 (ok):',
c    *          '  i,j,k,iterlag=', i,j,k, iterlag,
c    *          '  xs,ys,zp,up,vp,wp=', 
c    *             nint(xsurf(i,j,k)*.001),nint(ysurf(i,j,k)*.001),zp,
c    *             up,vp,wp,
c    *          '  h[1-4]=', h(iwa,jwa),h(iwb,jwa),h(iwa,jwb),h(iwb,jwb)
                goto 200
              endif

c++++++++++++++++++
              enddo ! iterlag
c ++++++++++++++++++

c                If exceeded max iter (error), just record current 
c                location and go on.

              xsurf(i,j,k) = xp
              ysurf(i,j,k) = yp
              statu(i,j,k) = up
              statv(i,j,k) = vp
              statifm(i,j,k) = -3.
c             write (191,"(a /a,3i6,i8 /a,2f9.0,f8.4,3f10.5 /a,4f10.3)")
c    *          'Warning (lagrangeo18 -2):',
c    *          '  i,j,k,iterlag=', i,j,k, iterlag-1,
c    *          '  xp,yp,zp,up,vp,wp=', xp*.001,yp*.001,zp, up,vp,wp,
c    *          '  h[1-4]=', h(iwa,jwa),h(iwb,jwa),h(iwb,jwa),h(iwb,jwb)

  200         continue

c                Bilin. interpolate d18O of precip (o18precip) from
c                ice grid to (xsurf,ysurf), record as o18all (3-D array)

              zxp = ( xsurf(i,j,k) - xoffa + 0.5*nx*dx0 ) / dx0
              zyp = ( ysurf(i,j,k) - yoffa + 0.5*ny*dy0 ) / dy0
              ipa = nint(zxp)
              ipb = ipa + 1 
              wpi = ipb - (zxp+0.5)
              ipa = max (1, min (nx, ipa))
              ipb = max (1, min (nx, ipb))
 
              jpa = nint(zyp)
              jpb = jpa + 1 
              wpj = jpb - (zyp+0.5)
              jpa = max (1, min (ny, jpa))
              jpb = max (1, min (ny, jpb))

              zo =     wpi *    wpj * o18precip(ipa,jpa)
     *           + (1.-wpi)*    wpj * o18precip(ipb,jpa)
     *           +     wpi *(1.-wpj)* o18precip(ipa,jpb)
     *           + (1.-wpi)*(1.-wpj)* o18precip(ipb,jpb)

              o18all(i,j,k) = zo

              if (mod(i+j+k,100).eq.0) 
     *          write (6,'(a,5i8)') 
     *            'lagrangeo18: done i,j,k,iterlag,stattim:',
     *            i, j, k, iterlag, nint(stattim(i,j,k))
c>>>>>>>>>>>>>>>>
            enddo     ! k
          endif       ! h > 0
        enddo         ! i
      enddo           ! j
c>>>>>>>>>>>>>>>>

c        Sum all o18all, weighted by volume of each 3-D ice cell,
c        normalize by summed ice volume, report as mean d18O of ice.

c!!!!!!!!!!!!!!!
      zotot = 0.
      zvtot = 0.

      do j=1,ny
        do i=1,nx
          do k=1,nlev
            if (o18all(i,j,k).ne.1.e20) then
              zdv = darea(i,j)*h(i,j)*dzeta(k)
              zotot = zotot + o18all(i,j,k)*zdv
              zvtot = zvtot +               zdv
            endif
          enddo
        enddo
      enddo
      zotot = zotot / max (zvtot,.001) 

      write (6,'(/a/a,f12.2/a,i12/a,f12.4)') 
     *      'lagrangeo18:',
     *      '  timeice (yr) =', timeice,
     *      '  toti (km3)   =', nint(zvtot/1.e9),
     *      '  o18 (permil) =', zotot
c!!!!!!!!!!!!!!!

c        Print out diagnostic maps 
c        2-D:
c        statvol = summed volume of all 3-D ice cells emerging 
c                  within each grid box 
c        3-D (loop over k):
c        o18all, o18all-o18precip, i_xsurf, j_ysurf, statu, statv, w. 
c        stattim = total time (yrs) taken along each path.
c        statnum = number of Lagrangian steps taken along path.
c        statifm: 1=top (normal) exit, -1=side exit, -2=exceed iter.

      iu = 192

      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.) then
            do k=1,nlev
              iso = nint ((xsurf(i,j,k) - xoffa + 0.5*nx*dx0)/dx0 + 0.5)
              jso = nint ((ysurf(i,j,k) - yoffa + 0.5*ny*dy0)/dy0 + 0.5)
              iso = max (1, min(nx,iso))
              jso = max (1, min(ny,jso))
              statvol(iso,jso) = statvol(iso,jso) 
     *                         + darea(i,j)*dzeta(k)*h(i,j)/1.e9 
            enddo
          endif
        enddo
      enddo
      call printmap (timeice, statvol, 'statvol (km3)', 500., iu ,0)

c----------------
      do k=1,nlev
c----------------

        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, o18all(1,1,k), 'o18all', 2., iu ,0)

        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              worka(i,j) = o18all(i,j,k) - o18precip(i,j)
            else
              worka(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice, worka, 'o18all-o18precip', 1., iu ,0)

        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              worka(i,j) = (xsurf(i,j,k) - xoffa) / dx0
              workb(i,j) = (ysurf(i,j,k) - yoffa) / dy0
            else
              worka(i,j) = 0.
              workb(i,j) = 0.
            endif
          enddo
        enddo
        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, worka, 'xsurf', 10., iu ,0)
        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, workb, 'ysurf', 10., iu ,0)
        write (iu,'(a,i4)') 'k=',k 

        call printmap (timeice, statu(1,1,k), 'statu', 1., iu ,0)
        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, statv(1,1,k), 'statv', 1., iu ,0)
        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, statifm(1,1,k), 'statifm', 1., iu ,0)
        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, stattim(1,1,k), 'stattim', 10000.,iu ,0)
        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, statnum(1,1,k), 'statnum', 100., iu ,0)
        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, w(1,1,k), 'w (1/yr)', .0001, iu ,0)
c----------
      enddo
c----------
      call flush (iu)

      return
      end

c-----------------------------------------------------------------------

      subroutine geto18p (o18precip, iup, timeice)

c         Read d18O of precip from GCM ice driving file, made by 
c         debplotgcmnet ('O18', ipole=-3, data_icedrive_O18 (local)).
c         Bilin. interp from GCM to ice grid. Like in climate_gcm.

#include <comicephys.h>
#include <comicegrid.h>

      dimension o18precip(nx,ny)

      dimension o18pdrive(nlong,nlatg)

c     local 2-D ice-grid arrays, indexes and weights to GCM grid:
      dimension
     *  indl(nx,ny),     weil(nx,ny),
     *  indb(nx,ny),     weib(nx,ny)
c     save indl, indb, weil, weib
     
      call setinterp (along, alatg, nlong, nlatg,
     *                alond, alatd, nx, ny,
     *                indl, weil, indb, weib, rotate_to_gcm)

      call readdrive (iuto, 'data_icedrive_O18', o18pdrive, 'O18', 1)

      do j=1,ny
        do i=1,nx
c         indices and weights into agcm-grid arrays (setinterp):
          il = indl(i,j)
          ir = mod (il,nlong) + 1
          jb = indb(i,j)
          jt = min (jb+1, nlatg)
          wlb =     weil(i,j) *    weib(i,j)
          wrb = (1.-weil(i,j))*    weib(i,j)
          wlt =     weil(i,j) *(1.-weib(i,j))
          wrt = (1.-weil(i,j))*(1.-weib(i,j))

c         zo1 = o18pdrive(il,jb)*wlb 
c         zo2 = o18pdrive(ir,jb)*wrb
c         zo3 = o18pdrive(il,jt)*wlt 
c         zo4 = o18pdrive(ir,jt)*wrt
c         o18precip(i,j) = zo1 + zo2 + zo3 + zo4

          o18precip(i,j) = o18pdrive(il,jb)*wlb + o18pdrive(ir,jb)*wrb
     *                   + o18pdrive(il,jt)*wlt + o18pdrive(ir,jt)*wrt

        enddo
      enddo

      call printmap (timeice, o18precip, 'o18precip (permil)', 2.,iup,0)

      return
      end

c@@@@@
#endif
c@@@@@

c-----------------------------------------------------------------------

#if defined (RWRAP_TEST)

c        Dummy R-wrapper (R.Fuller, B.Lee, 4/18). 
c        Sets input params (rwrap_param(1:11)),
c        calls subr rwrap_sheetshelf,
c        writes outputs (rwrap_[plio,...]) to iurwrap.

      program rwrapper

#include <comicegrid.h>

      parameter (nrwrap_param = 11)
      dimension rwrap_param(nrwrap_param),         ! input
     *          rwrap_plio(3),                     ! output
     *          rwrap_lig(3),
     *          rwrap_lgm(3),
     *          rwrap_mod(3),
     *          rwrap_obs(3),
     *          rwrap_2100(3),
     *          rwrap_2200(3),
     *          rwrap_2300(3),
     *          rwrap_2400(3),
     *          rwrap_2500(3),
     *          rwrap_mod_h(nx,ny), rwrap_mod_mask(nx,ny),
     *          rwrap_obs_h(nx,ny), rwrap_obs_mask(nx,ny)

      character*1 cmap(nx,ny)

      parameter (iurwrap=444)

c old:
c     rwrap_param(1) = 1.      ! OCFACMULT     [.3, 1, 3]
c     rwrap_param(2) = 1.      ! OCFACMULTASE  [1, 3, 10]
c     rwrap_param(3) = 1.      ! CALVNICK      [.5, 1, 2]
c     rwrap_param(4) = 1.e-5   ! CRHSHELF      [1.e-7, 1.e-5, 1.e-4]
c     rwrap_param(5) = 3000.   ! TAUASTH       [1000, 3000, 5000]
c     rwrap_param(6) = 100.    ! CALVLIQ       [0, 100, 150]
c     rwrap_param(7) = 12.e3   ! CLIFFVMAX     [0.e3, 3.e3, 12.e3]
c     rwrap_param(8) = 10.     ! FACEMELTRATE  [3, 10, 30]
c     rwrap_param(9) = 1.      ! ENHANCESHEET  [.1, 1, 10]
c     rwrap_param(10)= .3      ! ENHANCESHELF  [.03, .3, 3]
c     rwrap_param(11)= 1.      ! crhfac        [.01, 1., 100.]

      rwrap_param(1) = 1.      ! OCFACMULT     [1.e-.5, 1.e0, 1.e.5]
      rwrap_param(2) = 1.      ! OCFACMULTASE  [1.e0, 1.e.5, 1.e1]
      rwrap_param(3) = 1.      ! CALVNICK      [0, 1, 2]
      rwrap_param(4) = 1.e-5   ! CRHSHELF      [1.e-7, 1.e-5.5, 1.e-4]
      rwrap_param(5) = 3000.   ! TAUASTH       [1000, 3000, 5000]
      rwrap_param(6) = 100.    ! CALVLIQ       [0, 100, 200]
      rwrap_param(7) = 12.e3   ! CLIFFVMAX     [0, 6000, 12000]
      rwrap_param(8) = 10.     ! FACEMELTRATE  [0, 10, 20] 
      rwrap_param(9) = 1.      ! ENHANCESHEET  [1.e-1, 1.e0, 1.e1]
      rwrap_param(10)= .3      ! ENHANCESHELF  [.3e-1, .3e0, .3e1]
      rwrap_param(11)= 1.      ! crhfac        [1.e-2, 1.e0, 1.e2]

#if defined (RWRAP_PLIO)
      write (iurwrap,'(a)') 'R-wrapper test: PLIO'
#elif defined (RWRAP_LIG)
      write (iurwrap,'(a)') 'R-wrapper test: LIG'
#elif defined (RWRAP_LGMFUT)
      write (iurwrap,'(a)') 'R-wrapper test: LGMFUT'
#else
      write (6,'(/a)') '*** Error (rwrapper): unknown -DRWRAP_...'
      stop
#endif

      write (iurwrap,'(/a)') 'Inputs (rwrap_param):'
      do m=1,nrwrap_param
         write (iurwrap,'(i4,e15.5)') m, rwrap_param(m)
      enddo

      call flush (iurwrap)

      call rwrap_sheetshelf (rwrap_param, 
     *                       rwrap_plio,
     *                       rwrap_lig,
     *                       rwrap_lgm,
     *                       rwrap_mod,
     *                       rwrap_obs,
     *                       rwrap_2100,
     *                       rwrap_2200,
     *                       rwrap_2300,
     *                       rwrap_2400,
     *                       rwrap_2500,
     *                       rwrap_mod_h, rwrap_mod_mask,
     *                       rwrap_obs_h, rwrap_obs_mask,
     *                       rwrap_mod_rms)

      write (iurwrap,'(/a/)') 'Outputs:'

c        Overal volume, area, esl, desl

      write (iurwrap,'(6x,5a)') '    vol(km3)', '  areag(km2)',
     *          '      esl(m)', '     desl(m)', '     rmsh(m)'

      write (iurwrap,'(a,2i12,2f12.3)') 'PLIO: ', 
     *  nint(rwrap_plio(1)*1.e-9), nint(rwrap_plio(2)*1.e-6),
     *  rwrap_plio(3), rwrap_obs(3)-rwrap_plio(3)

      write (iurwrap,'(a,2i12,2f12.3)') 'LIG : ', 
     *  nint(rwrap_lig(1)*1.e-9), nint(rwrap_lig(2)*1.e-6), 
     *  rwrap_lig(3), rwrap_obs(3)-rwrap_lig(3)

      write (iurwrap,'(a,2i12,2f12.3)') 'LGM : ', 
     *  nint(rwrap_lgm(1)*1.e-9), nint(rwrap_lgm(2)*1.e-6), 
     *  rwrap_lgm(3), rwrap_obs(3)-rwrap_lgm(3)

      write (iurwrap,'(a,2i12,3f12.3)') 'MOD : ', 
     *  nint(rwrap_mod(1)*1.e-9), nint(rwrap_mod(2)*1.e-6), 
     *  rwrap_mod(3), rwrap_obs(3)-rwrap_mod(3), rwrap_mod_rms

      write (iurwrap,'(a,2i12,2f12.3)') 'OBS : ', 
     *  nint(rwrap_obs(1)*1.e-9), nint(rwrap_obs(2)*1.e-6), 
     *  rwrap_obs(3), 0.

      write (iurwrap,'(a,2i12,2f12.3)') '2100: ', 
     *  nint(rwrap_2100(1)*1.e-9),nint(rwrap_2100(2)*1.e-6),
     *  rwrap_2100(3), rwrap_obs(3)-rwrap_2100(3)

      write (iurwrap,'(a,2i12,2f12.3)') '2200: ', 
     *  nint(rwrap_2200(1)*1.e-9),nint(rwrap_2200(2)*1.e-6),
     *  rwrap_2200(3), rwrap_obs(3)-rwrap_2200(3)

      write (iurwrap,'(a,2i12,2f12.3)') '2300: ', 
     *  nint(rwrap_2300(1)*1.e-9),nint(rwrap_2300(2)*1.e-6),
     *  rwrap_2300(3), rwrap_obs(3)-rwrap_2300(3)

      write (iurwrap,'(a,2i12,2f12.3)') '2400: ', 
     *  nint(rwrap_2400(1)*1.e-9),nint(rwrap_2400(2)*1.e-6),
     *  rwrap_2400(3), rwrap_obs(3)-rwrap_2400(3)

      write (iurwrap,'(a,2i12,2f12.3)') '2500: ', 
     *  nint(rwrap_2500(1)*1.e-9),nint(rwrap_2500(2)*1.e-6),
     *  rwrap_2500(3), rwrap_obs(3)-rwrap_2500(3)

c        Modern h maps, model and observed:

      write (iurwrap,'(/a)') 'mod_h (10 m):' 
      write (iurwrap,'(6x,70i8)') (nint(.001*xh(i,(ny+1)/2)), i=1,nx,2)
      do j=ny,1,-1
        write (iurwrap,'(i8,2x,70i4)') nint(.001*yh((nx+1)/2,j)),  
     *                               (nint(0.1*rwrap_mod_h(i,j)),i=1,nx)
      enddo

      write (iurwrap,'(/a)') 'obs_h (10 m):' 
      write (iurwrap,'(6x,70i8)') (nint(.001*xh(i,(ny+1)/2)), i=1,nx,2)
      do j=ny,1,-1
        write (iurwrap,'(i8,2x,70i4)') nint(.001*yh((nx+1)/2,j)),  
     *                               (nint(0.1*rwrap_obs_h(i,j)),i=1,nx)
      enddo

c        Modern mask maps, model and observed:

      do j=1,ny
        do i=1,nx
          if (rwrap_mod_h(i,j).gt.0.) then
            if (nint(rwrap_mod_mask(i,j)).eq.0) then
              cmap(i,j) = 'x'
            else
              cmap(i,j) = '.'
            endif
          else
            cmap(i,j) = ' '
          endif
        enddo
      enddo
      write (iurwrap,'(/a)') 'mod_mask (x = grounded, . = floating):' 
      write (iurwrap,'(3x,70i8)') (nint(.001*xh(i,(ny+1)/2)),i=1,nx,8)
      do j=ny,1,-1
        write (iurwrap,'(i8,2x,70a1)') nint(.001*yh((nx+1)/2,j)),  
     *                                 (cmap(i,j),i=1,nx)
      enddo

      do j=1,ny
        do i=1,nx
          if (rwrap_obs_h(i,j).gt.0.) then
            if (nint(rwrap_obs_mask(i,j)).eq.0) then
              cmap(i,j) = 'x'
            else
              cmap(i,j) = '.'
            endif
          else
            cmap(i,j) = ' '
          endif
        enddo
      enddo
      write (iurwrap,'(/a)') 'obs_mask (x = grounded, . = floating):' 
      write (iurwrap,'(3x,70i8)') (nint(.001*xh(i,(ny+1)/2)),i=1,nx,8)
      do j=ny,1,-1
        write (iurwrap,'(i8,2x,70a1)') nint(.001*yh((nx+1)/2,j)),  
     *                                 (cmap(i,j),i=1,nx) 
      enddo

      call flush (iurwrap)

      stop
      end

#endif

c-----------------------------------------------------------------------
