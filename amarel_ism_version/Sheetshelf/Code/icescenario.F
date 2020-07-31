cccccccccccccccccccccccccccccccccccccccccccccc
#if defined (RCMMATRIX) && defined (RCMMATCO2)
cccccccccccccccccccccccccccccccccccccccccccccc

c-----------------------------------------------------------------------

      subroutine scenario (timeice, nyearstart, weirun, rco2)

c        Sets rco2, given timeice (yr) for various future scenarios,
c        selected by one of RCPSPIN*, RCP26*, RCP45*, RCP85*, RCPCOP*.
c        These are merged from old (dp) scenarios and new (rd), 6/18.
c
c        Use timeice-nyearstart for all (new) scenarios, for which
c        nyearstart=0 always(?), corresp. to 1950.
c        Not for RCP*OLD and RAMP* for which timeice=0 corresp. to 1950.
c        So timeice (or timeice-nyearstart for new scenarios) always 
c        corresponds to 1950.
c
c        For RCP*OLD and some new scenarios, use natural atmos co2 
c        decay back to rco2=1 after emissions end. 
c        For RCP*OLD, set by stmt fn farch (Archer).
c
c        For RCP*OLD, use time lag (rcplag, years) between CO2 amount 
c        and RCM climates.
c
c        -D's for merging (see notes_codediffs, msg_2rd_codediffs):
c          COP21RD for minor physic diffs (icedyn.F, iceclimall.F).
c          DPOLD (for old dp input files, not rd,in iceclimgcm,iceocean)
c          RCP[SPIN,26,45,85]OLD for old dp without "OLD" (iceclimgcm.F)
c          RCP[26,45,85]OCN in iceocean.F (was RCP[26,45,85] in old dp).
c
c        And redone -D's for RACMO climates,Gr+Ant (see notes 3/16/18):
c          Old:          New:
c          RCMRACMO      RACMOCLIM (with GREENLAND)
c          RACANTCLIM    RACMOCLIM (with EISANTA)
c          RACANTOBS     RACMOANTSEAS. Ises Tann_searise + dTseas_racmo,
c                        only for Antarctica (EISANTA) with SEARISECLIM,
c                        or for [RCM,CSM,LOVE]CLIM with [...]ANOM.

#if defined (NOARCHER)
      farch (ztb) = 1.
#elif defined (RAPIDARCHER)
      farch (ztb) = exp (-ztb/100.)
#else
      farch (ztb) =   0.750*exp(-ztb/365.)  + .135*exp(-ztb/5.5e3)
     *              + 0.035*exp(-ztb/8.2e3) + .080*exp(-ztb/200.e3)
#endif

      weirun = 1. ! default (for all except RAMP*)

      rcplag = 10. ! only for RCP*OLD

crrrrrrrrrrrrrrrrrrrrrrr
#if defined (RCPSPINOLD)
crrrrrrrrrrrrrrrrrrrrrrr

c     dp old:
   
c        For spinup from 1850 (preind, rco2=1.0) to 1950 (rco2=1.09)

c     set zta = time since 1850:
      zta = max (0., timeice + 100. - rcplag) 
   
      if (zta.le.100.) then
        zw = zta/100.
        rco2 = 1.*(1.-zw) + 1.09 * zw 
      else
        rco2 = 1.09
      endif

crrrrrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCPSPIN2000)
crrrrrrrrrrrrrrrrrrrrrrrrrr

c       For spinup searise run with CCSM ocean  using Run2000 climate 

      rco2 = 1.0 

crrrrrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCPSPIN1950)
crrrrrrrrrrrrrrrrrrrrrrrrrr

c       For spinup from 1950 (rco2=1.08) to 2000 (rco2=1.30)

      zta = timeice-nyearstart

      if (zta.le.25.) then
        zw = zta/25.
        rco2 = 1.08*(1.-zw) + 1.15 * zw

      else if (zta.le.50.) then
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.15*(1.-zw) + 1.30 * zw

      else
        rco2 = 1.30
      endif

crrrrrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCPSPIN1850)
crrrrrrrrrrrrrrrrrrrrrrrrrr

c       For spinup from 1850 (rco2=1.02) to 2000 (rco2=1.30)

c     zta = timeice-nyearstart ! rd (?)

c     set zta = time since 1850:
      zta = max (0., timeice-nyearstart + 100.)

      if (zta.le.25.) then 
        zw = zta/25.
        rco2 = 1.02*(1.-zw) + 1.01 * zw 

      else if (zta.le.50.) then 
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.01*(1.-zw) + 1.03 * zw 

      else if (zta.le.75.) then
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.03*(1.-zw) + 1.06 * zw

      else if (zta.le.100.) then
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.06*(1.-zw) + 1.08 * zw

      else if (zta.le.125.) then
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.08*(1.-zw) + 1.15 * zw

      else if (zta.le.150.) then
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.15*(1.-zw) + 1.30 * zw

      else 
        rco2 = 1.30 
      endif

crrrrrrrrrrrrrrrrrrrrr
#elif defined (RCPFIX) 
crrrrrrrrrrrrrrrrrrrrr

c        For future RCP runs with constant CO2: 

      rco2 = 1.30 

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCP26OLD)
crrrrrrrrrrrrrrrrrrrrrrr

c     dp old:

c     for future RCP runs with constant or varying CO2 after Xyrs:

      zta = max (0., timeice - rcplag)

      if (zta.le.25.) then 
        zw = zta/25.
        rco2 = 1.09*(1.-zw) + 1.09 * zw 

      else if (zta.le.50.) then 
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.09*(1.-zw) + 1.25 * zw 

      else if (zta.le.75.) then
        ztb = zta - 50.
        zw = ztb/25.
        rco2 = 1.25*(1.-zw) + 1.45 * zw

      else if (zta.le.100.) then
        ztb = zta - 75.
        zw = ztb/25.
        rco2 = 1.45*(1.-zw) + 1.63 * zw

      else if (zta.le.125.) then 
        ztb = zta - 100. 
        zw = ztb/25.
        rco2 = 1.63-1.63*(0.+zw) + 1.59 * zw 

      else if (zta.le.150.) then
        ztb = zta - 125. 
        zw = ztb/25.
        rco2 = 1.59-1.59*(0.+zw) + 1.54 * zw

      else if (zta.le.200.) then 
        ztb = zta - 150. 
        zw = ztb/50.
        rco2 = 1.54-1.54*(0.+zw) + 1.45 * zw 

      else if (zta.le.250.) then 
        ztb = zta - 200. 
        zw = ztb/50.
        rco2 = 1.45-1.45*(0.+zw) + 1.38 * zw 

      else if (zta.le.300.) then 
        ztb = zta - 250. 
        zw = ztb/50.
        rco2 = 1.38-1.38*(0.+zw) + 1.33* zw

      else if (zta.le.350.) then
        ztb = zta - 300.
        zw = ztb/50.
        rco2 = 1.33-1.33*(0.+zw) + 1.29* zw

      else if (zta.le.550.) then
        ztb = zta - 350.
        zw = ztb/200.
        rco2 = 1.29-1.29*(0.+zw) + 1.17* zw

      else
        ztb = zta - 550.
        rco2 = 1.17*farch(ztb) + 1.*(1.-farch(ztb))
      endif

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCP26ATM)
crrrrrrrrrrrrrrrrrrrrrrr

      zta = timeice-nyearstart
      if (zta.le.10.) then
        zw = zta/10.
        rco2 = 1.08*(1.-zw) + 1.08 * zw

      else if (zta.le.20.) then
        ztb = zta - 10.
        zw = ztb/10.
        rco2 = 1.08*(1.-zw) + 1.11 * zw

      else if (zta.le.30.) then
        ztb = zta - 20.
        zw = ztb/10.
        rco2 = 1.11*(1.-zw) + 1.18 * zw

      else if (zta.le.40.) then
        ztb = zta - 30.
        zw = ztb/10.
        rco2 = 1.18*(1.-zw) + 1.25 * zw

      else if (zta.le.50.) then
       ztb = zta - 40.
       zw = ztb/10.
       rco2 = 1.25*(1.-zw) + 1.30 * zw

      else if (zta.le.60.) then
        ztb = zta - 50.
        zw = ztb/10.
        rco2 = 1.30*(1.-zw) + 1.39 * zw

      else if (zta.le.70.) then
        ztb = zta - 60.
        zw = ztb/10.
        rco2 = 1.39*(1.-zw) + 1.51 * zw

      else if (zta.le.80.) then
        ztb = zta - 70.
        zw = ztb/10.
        rco2 = 1.51*(1.-zw) + 1.59 * zw

      else if (zta.le.90.) then
        ztb = zta - 80.
        zw = ztb/10.
        rco2 = 1.59*(1.-zw) + 1.63 * zw

      else if (zta.le.100.) then
        ztb = zta - 90.
        zw = ztb/10.
        rco2 = 1.63*(1.-zw) + 1.63 * zw

      else if (zta.le.125.) then
        ztb = zta - 100.
        zw = ztb/25.
        rco2 = 1.63-1.63*(0.+zw) + 1.57 * zw

      else if (zta.le.150.) then
        ztb = zta - 125.
        zw = ztb/25.
        rco2 = 1.57-1.57*(0.+zw) + 1.53 * zw

      else if (zta.le.200.) then
        ztb = zta - 150.
        zw = ztb/50.
        rco2 = 1.53-1.53*(0.+zw) + 1.43 * zw

      else if (zta.le.250.) then
        ztb = zta - 200.
        zw = ztb/50.
        rco2 = 1.43-1.43*(0.+zw) + 1.37 * zw

      else if (zta.le.300.) then
        ztb = zta - 250.
        zw = ztb/50.
        rco2 = 1.37-1.37*(0.+zw) + 1.32* zw

      else if (zta.le.350.) then
        ztb = zta - 300.
        zw = ztb/50.
        rco2 = 1.32-1.32*(0.+zw) + 1.29* zw

      else if (zta.le.550.) then
        ztb = zta - 350.
        zw = ztb/200.
        rco2 = 1.29-1.29*(0.+zw) + 1.16* zw

      else
        rco2 = 1.16
      endif

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCP45OLD)
crrrrrrrrrrrrrrrrrrrrrrr

c     dp old:

      zta = max (0., timeice - rcplag)

      if (zta.le.25.) then
        zw = zta/25.
        rco2 = 1.09*(1.-zw) + 1.09 * zw

      else if (zta.le.50.) then
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.09*(1.-zw) + 1.25 * zw

      else if (zta.le.75.) then
        ztb = zta - 50.
        zw = ztb/25.
        rco2 = 1.25*(1.-zw) + 1.45 * zw

      else if (zta.le.100.) then
        ztb = zta - 75.
        zw = ztb/25.
        rco2 = 1.45*(1.-zw) + 1.76 * zw

      else if (zta.le.125.) then
        ztb = zta - 100.
        zw = ztb/25.
        rco2 = 1.76*(1.-zw) + 2.01 * zw

      else if (zta.le.150.) then
        ztb = zta - 125.
        zw = ztb/25.
        rco2 = 2.01*(1.-zw) + 2.06 * zw

      else if (zta.le.200.) then
        ztb = zta - 150.
        zw = ztb/50.
        rco2 = 2.08*(1.-zw) + 2.08 * zw    

      else if (zta.le.550.) then        ! level till 2500 CE, as in DP16
        rco2 = 2.08

      else 
        ztb = zta - 550.
        rco2 = 2.08*farch(ztb) + 1.*(1.-farch(ztb))
      endif

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCP45ATM)
crrrrrrrrrrrrrrrrrrrrrrr

      zta = timeice-nyearstart
      if (zta.le.10.) then 
        zw = zta/10.
        rco2 = 1.08*(1.-zw) + 1.08 * zw 

      else if (zta.le.20.) then 
        ztb = zta - 10.
        zw = ztb/10.
        rco2 = 1.08*(1.-zw) + 1.11 * zw 

      else if (zta.le.30.) then 
        ztb = zta - 20.
        zw = ztb/10.
        rco2 = 1.11*(1.-zw) + 1.18 * zw 

      else if (zta.le.40.) then 
        ztb = zta - 30.
        zw = ztb/10.
        rco2 = 1.18*(1.-zw) + 1.25 * zw 

      else if (zta.le.50.) then 
        ztb = zta - 40.
        zw = ztb/10.
        rco2 = 1.25*(1.-zw) + 1.30 * zw 

      else if (zta.le.60.) then 
        ztb = zta - 50.
        zw = ztb/10.
        rco2 = 1.30*(1.-zw) + 1.39 * zw 

      else if (zta.le.70.) then
        ztb = zta - 60.
        zw = ztb/10.
        rco2 = 1.39*(1.-zw) + 1.51 * zw

      else if (zta.le.80.) then
        ztb = zta - 70.
        zw = ztb/10.
        rco2 = 1.51*(1.-zw) + 1.63 * zw

      else if (zta.le.90.) then
        ztb = zta - 80.
        zw = ztb/10.
        rco2 = 1.63*(1.-zw) + 1.76 * zw

      else if (zta.le.100.) then
        ztb = zta - 90.
        zw = ztb/10.
        rco2 = 1.76*(1.-zw) + 1.88 * zw

      else if (zta.le.125.) then
        ztb = zta - 100.
        zw = ztb/25.
        rco2 = 1.88*(1.-zw) + 2.05 * zw

      else if (zta.le.150.) then
        ztb = zta - 125.
        zw = ztb/25.
        rco2 = 2.05*(1.-zw) + 2.08 * zw

      else if (zta.le.200.) then
        ztb = zta - 150.
        zw = ztb/50.
        rco2 = 2.08*(1.-zw) + 2.08 * zw

      else if (zta.le.550.) then
        rco2 = 2.08

      else

c       rco2 = 2.08

c       fast co2 scrubbing with 100-yr efolding time:
        ztb = zta - 550.
        zw = exp(-ztb/100)
        rco2 = 1.*(1.-zw) + 2.08* zw

      endif

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCP85OLD)
crrrrrrrrrrrrrrrrrrrrrrr

c     dp old:

      zta = max (0., timeice - rcplag)

      if (zta.le.25.) then
        zw = zta/25.
        rco2 = 1.09*(1.-zw) + 1.09 * zw

      else if (zta.le.50.) then
        ztb = zta - 25.
        zw = ztb/25.
        rco2 = 1.09*(1.-zw) + 1.25 * zw

      else if (zta.le.75.) then
        ztb = zta - 50.
        zw = ztb/25.
        rco2 = 1.25*(1.-zw) + 1.46 * zw

      else if (zta.le.100.) then
        ztb = zta - 75.
        zw = ztb/25.
        rco2 = 1.46*(1.-zw) + 1.95 * zw

      else if (zta.le.125.) then
        ztb = zta - 100.
        zw = ztb/25.
        rco2 = 1.95*(1.-zw) + 2.78 * zw

      else if (zta.le.150.) then
        ztb = zta - 125.
        zw = ztb/25.
        rco2 = 2.78*(1.-zw) + 3.90 * zw

      else if (zta.le.200.) then
        ztb = zta - 150.
        zw = ztb/50.
        rco2 = 3.90*(1.-zw) + 6.34 * zw    

      else if (zta.le.250.) then
        ztb = zta - 200.
        zw = ztb/50.
        rco2 = 6.34*(1.-zw) + 8.0 * zw

      else if (zta.le.550.) then        ! level till 2500 CE, as in DP16
        rco2 = 8.0

      else
        ztb = zta - 550.
        rco2 = 8.0*farch(ztb) + 1.*(1.-farch(ztb))
      endif

crrrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCP85SATM)
crrrrrrrrrrrrrrrrrrrrrrrr

       zta = timeice-nyearstart
       if (zta.le.30.) then 
        zw = zta/30.
        rco2 = 2.995*(1.-zw) + 4.40 * zw 

      else if (zta.le.80.) then
        ztb = zta - 30.
        zw = ztb/50.
        rco2 = 4.40*(1.-zw) + 6.84 * zw

      else if (zta.le.130.) then
        ztb = zta - 80.
        zw = ztb/50.
        rco2 = 6.84*(1.-zw) + 8.79 * zw

      else if (zta.le.180.) then
        ztb = zta - 130.
        zw = ztb/50.
        rco2 = 8.79*(1.-zw) + 9.43 * zw

c     else if (zta.le.550.) then
      else
        rco2 = 9.43
      endif

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCP85ATM)
crrrrrrrrrrrrrrrrrrrrrrr

c        This new RCP sequence begins in 1950 as in DP16,
c        There is no longer an imposed lag in rCO2. 
c        nmat is now 5.
c        1950: rco2 = 1.08, 2000: rco2 1.3, 2x, 4x, 8x
c        equivalent CO2 from Meisnhausen et al., 2011

      zta = timeice-nyearstart
      if (zta.le.10.) then
        zw = zta/10.
        rco2 = 1.08*(1.-zw) + 1.08 * zw

      else if (zta.le.20.) then 
        ztb = zta - 10.
        zw = ztb/10.
        rco2 = 1.08*(1.-zw) + 1.11 * zw 

      else if (zta.le.30.) then 
        ztb = zta - 20.
        zw = ztb/10.
        rco2 = 1.11*(1.-zw) + 1.18 * zw 

      else if (zta.le.40.) then 
        ztb = zta - 30.
        zw = ztb/10.
        rco2 = 1.18*(1.-zw) + 1.25 * zw 

      else if (zta.le.50.) then 
        ztb = zta - 40.
        zw = ztb/10.
        rco2 = 1.25*(1.-zw) + 1.30 * zw 

      else if (zta.le.60.) then 
        ztb = zta - 50.
        zw = ztb/10.
        rco2 = 1.30*(1.-zw) + 1.40 * zw

      else if (zta.le.70.) then
        ztb = zta - 60.
        zw = ztb/10.
        rco2 = 1.40*(1.-zw) + 1.53 * zw

      else if (zta.le.80.) then
        ztb = zta - 70.
        zw = ztb/10.
        rco2 = 1.53*(1.-zw) + 1.71 * zw

      else if (zta.le.90.) then
        ztb = zta - 80.
        zw = ztb/10.
        rco2 = 1.71*(1.-zw) + 1.95 * zw

      else if (zta.le.100.) then
        ztb = zta - 90.
        zw = ztb/10.
        rco2 = 1.95*(1.-zw) + 2.24 * zw

      else if (zta.le.125.) then
        ztb = zta - 100.
        zw = ztb/25.
        rco2 = 2.24*(1.-zw) + 3.19 * zw

      else if (zta.le.150.) then
        ztb = zta - 125.
        zw = ztb/25.
        rco2 = 3.19*(1.-zw) + 4.40 * zw

      else if (zta.le.175.) then
        ztb = zta - 150.
        zw = ztb/25.
        rco2 = 4.40*(1.-zw) + 5.61 * zw

      else if (zta.le.200.) then 
        ztb = zta - 175. 
        zw = ztb/25.
        rco2 = 5.61*(1.-zw) + 6.84 * zw 

      else if (zta.le.225.) then 
        ztb = zta - 200. 
        zw = ztb/25.
        rco2 = 6.84*(1.-zw) + 7.97 * zw 

      else if (zta.le.250.) then
        ztb = zta - 225.
        zw = ztb/25.
        rco2 = 7.97*(1.-zw) + 8.79 * zw

      else if (zta.le.275.) then 
        ztb = zta - 250. 
        zw = ztb/25.
        rco2 = 8.79*(1.-zw) + 9.28 * zw 

      else if (zta.le.300.) then 
        ztb = zta - 275. 
        zw = ztb/25.
        rco2 = 9.28*(1.-zw) + 9.43 * zw 

      else if (zta.le.550.) then
        rco2 = 9.43

      else

c       fast co2 scrubbing with 100-yr e-folding time:
        ztb = zta - 550.
        zw = exp(-ztb/100)
        rco2 = 1.*(1.-zw) + 8.*zw
c       rco2 = 1.*(1.-zw) + 9.43*zw ?

      endif

crrrrrrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RCP85ATM1900)
crrrrrrrrrrrrrrrrrrrrrrrrrrr

c         Begins in 1900.
c         There is no longer an imposed lag in rCO2.
c         nmat is 5.
c         1950 rco2=1.08, 2000 rco2= 1.3, 2x, 4x, 8x
c         equivalent CO2 from Meisnhausen et al., 2011

      zta = timeice-nyearstart
      if (zta.le.10.) then
        zw = zta/10.
        rco2 = 1.03*(1.-zw) + 1.03 * zw

      else if (zta.le.20.) then
        ztb = zta - 10.
        zw = ztb/10.
        rco2 = 1.03*(1.-zw) + 1.05 * zw

      else if (zta.le.30.) then
        ztb = zta - 20.
        zw = ztb/10.
        rco2 = 1.05*(1.-zw) + 1.07 * zw

      else if (zta.le.40.) then
        ztb = zta - 30.
        zw = ztb/10.
        rco2 = 1.07*(1.-zw) + 1.09 * zw

      else if (zta.le.50.) then
        ztb = zta - 40.
        zw = ztb/10.
        rco2 = 1.09-1.09*(0.+zw) + 1.08 * zw

      else if (zta.le.60.) then
        ztb = zta - 50.
        zw = ztb/10.
        rco2 = 1.08*(1.-zw) + 1.08 * zw

      else if (zta.le.70.) then
        ztb = zta - 60.
        zw = ztb/10.
        rco2 = 1.08*(1.-zw) + 1.11 * zw

      else if (zta.le.80.) then
        ztb = zta - 70.
        zw = ztb/10.
        rco2 = 1.11*(1.-zw) + 1.18 * zw

      else if (zta.le.90.) then
        ztb = zta - 80.
        zw = ztb/10.
        rco2 = 1.18*(1.-zw) + 1.25 * zw

      else if (zta.le.100.) then
        ztb = zta - 90.
        zw = ztb/10.
        rco2 = 1.25*(1.-zw) + 1.30 * zw

       else if (zta.le.110.) then
        ztb = zta - 100.
        zw = ztb/10.
        rco2 = 1.30*(1.-zw) + 1.40 * zw

      else if (zta.le.120.) then
        ztb = zta - 110.
        zw = ztb/10.
        rco2 = 1.40*(1.-zw) + 1.53 * zw

      else if (zta.le.130.) then
        ztb = zta - 120.
        zw = ztb/10.
        rco2 = 1.53*(1.-zw) + 1.71 * zw

      else if (zta.le.140.) then
        ztb = zta - 130.
        zw = ztb/10.
        rco2 = 1.71*(1.-zw) + 1.95 * zw

      else if (zta.le.150.) then
        ztb = zta - 140.
        zw = ztb/10.
        rco2 = 1.95*(1.-zw) + 2.24 * zw

      else if (zta.le.175.) then
        ztb = zta - 150.
        zw = ztb/25.
        rco2 = 2.24*(1.-zw) + 3.19 * zw

      else if (zta.le.200.) then
        ztb = zta - 175.
        zw = ztb/25.
        rco2 = 3.19*(1.-zw) + 4.40 * zw

      else if (zta.le.225.) then
        ztb = zta - 200.
        zw = ztb/25.
        rco2 = 4.40*(1.-zw) + 5.61 * zw

      else if (zta.le.250.) then
        ztb = zta - 225.
        zw = ztb/25.
        rco2 = 5.61*(1.-zw) + 6.84 * zw

      else if (zta.le.275.) then
        ztb = zta - 250.
        zw = ztb/25.
        rco2 = 6.84*(1.-zw) + 7.97 * zw

      else if (zta.le.300.) then
        ztb = zta - 275.
        zw = ztb/25.
        rco2 = 7.97*(1.-zw) + 8.79 * zw

      else if (zta.le.325.) then
        ztb = zta - 300.
        zw = ztb/25.
        rco2 = 8.79*(1.-zw) + 9.28 * zw

      else if (zta.le.350.) then
        ztb = zta - 325.
        zw = ztb/25.
        rco2 = 9.28*(1.-zw) + 9.43 * zw

      else if (zta.le.600.) then
        rco2 = 9.43

      else

c       fast co2 scrubbing with 100-yr e-folding time:
        ztb = zta - 600.
        zw = exp(-ztb/100)
        rco2 = 1.*(1.-zw) + 8.* zw
c       rco2 = 1.*(1.-zw) + 9.43 * zw ?

      endif

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (COP30ATM)
crrrrrrrrrrrrrrrrrrrrrrr

c        Limits global mean temp to ~3.0 deg.
c        based on RCP8.5 but co2 held constant after 2070CE

      zta = timeice-nyearstart
      if (zta.le.10.) then
        zw = zta/10.
        rco2 = 1.08*(1.-zw) + 1.08 * zw

      else if (zta.le.20.) then
        ztb = zta - 10.
        zw = ztb/10.
        rco2 = 1.08*(1.-zw) + 1.11 * zw

      else if (zta.le.30.) then
        ztb = zta - 20.
        zw = ztb/10.
        rco2 = 1.11*(1.-zw) + 1.18 * zw

      else if (zta.le.40.) then
        ztb = zta - 30.
        zw = ztb/10.
        rco2 = 1.18*(1.-zw) + 1.25 * zw

      else if (zta.le.50.) then
        ztb = zta - 40.
        zw = ztb/10.
        rco2 = 1.25*(1.-zw) + 1.30 * zw

      else if (zta.le.60.) then
        ztb = zta - 50.
        zw = ztb/10.
        rco2 = 1.30*(1.-zw) + 1.40 * zw

      else if (zta.le.70.) then
        ztb = zta - 60.
        zw = ztb/10.
        rco2 = 1.40*(1.-zw) + 1.53 * zw

      else if (zta.le.80.) then
        ztb = zta - 70.
        zw = ztb/10.
        rco2 = 1.53*(1.-zw) + 1.71 * zw

      else if (zta.le.90.) then
        ztb = zta - 80.
        zw = ztb/10.
        rco2 = 1.71*(1.-zw) + 1.95 * zw

      else if (zta.le.100.) then
        ztb = zta - 90.
        zw = ztb/10.
        rco2 = 1.95*(1.-zw) + 2.24 * zw
 
c====
c icectl.F.default, no overshoot:
      else if (zta.le.120.) then
        ztb = zta - 100.
        zw = ztb/20.
        rco2 = 2.24*(1.-zw) + 2.98 * zw
      else
        rco2 = 2.98
      endif
c====
c icectl.F.over2050:
c     else if (zta.le.125.) then 
c       rco2 = 2.24
c     else
c       fast co2 scrubbing with 100-yr efolding time (all overshoots):
c       ztb = zta - 125.
c       zw = exp(-ztb/100)
c       rco2 = 1.*(1.-zw) + 2.94 * zw         ! rd has + 2.98 * zw
c     endif
c====
c icectl.F.over2070:
c     else if (zta.le.120.) then
c       ztb = zta - 100.
c       zw = ztb/20.
c       rco2 = 2.24*(1.-zw) + 2.98 * zw
c     else
c       ztb = zta - 120.                      ! rd has - 125.
c       zw = exp(-ztb/100)
c       rco2 = 1.*(1.-zw) + 2.98 * zw
c     endif
c====
c icectl.F.over2100:
c     else if (zta.le.120.) then
c       ztb = zta - 100.
c       zw = ztb/20.
c       rco2 = 2.24*(1.-zw) + 2.98 * zw
c     else if (zta.le.150.) then
c       rco2 = 2.98
c     else
c       ztb = zta - 150.
c       zw = exp(-ztb/100)
c       rco2 = 1.*(1.-zw) + 2.98 * zw
c     endif
c====
c icectl.F.over2150:
c     else if (zta.le.120.) then
c       ztb = zta - 100.
c       zw = ztb/20.
c       rco2 = 2.24*(1.-zw) + 2.98 * zw
c     else if (zta.le.200.) then
c       rco2 = 2.98
c     else
c       ztb = zta - 200.
c       zw = exp(-ztb/100)
c       rco2 = 1.*(1.-zw) + 2.98 * zw
c     endif
c====
c icectl.F.over2200:
c     else if (zta.le.120.) then
c       ztb = zta - 100.
c       zw = ztb/20.
c       rco2 = 2.24*(1.-zw) + 2.98 * zw
c     else if (zta.le.250.) then
c       rco2 = 2.98
c     else
c       ztb = zta - 250.
c       zw = exp(-ztb/100)
c       rco2 = 1.*(1.-zw) + 2.98 * zw
c     endif
c====

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (COP20ATM)
crrrrrrrrrrrrrrrrrrrrrrr

c        Limits global mean temp to ~2.0 deg.
c        based on RCP4.5, but no CO2 drop after ~2060CE

      zta = timeice-nyearstart
      if (zta.le.10.) then
        zw = zta/10.
        rco2 = 1.08*(1.-zw) + 1.08 * zw

      else if (zta.le.20.) then
        ztb = zta - 10.
        zw = ztb/10.
        rco2 = 1.08*(1.-zw) + 1.11 * zw

      else if (zta.le.30.) then
        ztb = zta - 20.
        zw = ztb/10.
        rco2 = 1.11*(1.-zw) + 1.18 * zw

      else if (zta.le.40.) then
        ztb = zta - 30.
        zw = ztb/10.
        rco2 = 1.18*(1.-zw) + 1.25 * zw

      else if (zta.le.50.) then
        ztb = zta - 40.
        zw = ztb/10.
        rco2 = 1.25*(1.-zw) + 1.30 * zw

      else if (zta.le.60.) then
        ztb = zta - 50.
        zw = ztb/10.
        rco2 = 1.30*(1.-zw) + 1.39 * zw

      else if (zta.le.70.) then
        ztb = zta - 60.
        zw = ztb/10.
        rco2 = 1.39*(1.-zw) + 1.51 * zw

      else if (zta.le.80.) then
        ztb = zta - 70.
        zw = ztb/10.
        rco2 = 1.51*(1.-zw) + 1.63 * zw

      else if (zta.le.90.) then
        ztb = zta - 80.
        zw = ztb/10.
        rco2 = 1.63*(1.-zw) + 1.76 * zw

      else if (zta.le.100.) then
        ztb = zta - 90.
        zw = ztb/10.
        rco2 = 1.76*(1.-zw) + 1.88 * zw

      else if (zta.le.110.) then
        ztb = zta - 100.
        zw = ztb/10.
        rco2 = 1.88*(1.-zw) + 1.97 * zw

c     else if (zta.le.550.) then
      else
        rco2 = 1.97
      endif

crrrrrrrrrrrrrrrrrrrrrrr
#elif defined (COP15ATM)
crrrrrrrrrrrrrrrrrrrrrrr

c        Limits global mean temp to ~1.5 deg.
c        based on RCP4.5 up to year 2040

      zta = timeice-nyearstart
      if (zta.le.10.) then
        zw = zta/10.
        rco2 = 1.08*(1.-zw) + 1.08 * zw

      else if (zta.le.20.) then
        ztb = zta - 10.
        zw = ztb/10.
        rco2 = 1.08*(1.-zw) + 1.11 * zw

      else if (zta.le.30.) then
        ztb = zta - 20.
        zw = ztb/10.
        rco2 = 1.11*(1.-zw) + 1.18 * zw

      else if (zta.le.40.) then
        ztb = zta - 30.
        zw = ztb/10.
        rco2 = 1.18*(1.-zw) + 1.25 * zw

      else if (zta.le.50.) then
        ztb = zta - 40.
        zw = ztb/10.
        rco2 = 1.25*(1.-zw) + 1.30 * zw

      else if (zta.le.60.) then
        ztb = zta - 50.
        zw = ztb/10.
        rco2 = 1.30*(1.-zw) + 1.39 * zw

      else if (zta.le.70.) then
        ztb = zta - 60.
        zw = ztb/10.
        rco2 = 1.39*(1.-zw) + 1.51 * zw

      else if (zta.le.80.) then
        ztb = zta - 70.
        zw = ztb/10.
        rco2 = 1.51*(1.-zw) + 1.63 * zw

      else if (zta.le.90.) then
        ztb = zta - 80.
        zw = ztb/10.
        rco2 = 1.63*(1.-zw) + 1.76 * zw

c     else if (zta.le.550.) then
      else
        rco2 = 1.76
      endif

crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr
#elif defined (RAMPPLIO) || defined (RAMP2X) || defined (RAMP4X) || defined (RAMP8X) || defined (RAMP8XW)
crrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr

#          if defined (RAMPMATCO2)
      zwt = max (0., min (1., timeice/RAMPMATCO2))
#          else
      zwt = 1.
#          endif

      weirun = 1.*(1.-zwt) + 2.*zwt                 ! +0 to +DTOCNMATCO2

#          if defined (RAMPPLIO)
      rco2 = exp((1.-zwt)*log(1.) + zwt*log(1.43))               !  plio
#          elif defined (RAMP2X)
      rco2 = exp((1.-zwt)*log(1.) + zwt*log(2.))                 !  2x
#          elif defined (RAMP4X)
      rco2 = exp((1.-zwt)*log(1.) + zwt*log(4.))                 !  4x
#          elif defined (RAMP8X) || defined (RAMP8XW)
      rco2 = exp((1.-zwt)*log(1.) + zwt*log(8.))                 !  8x
#          endif

crrrr
#else
crrrr

      write (6,'(/a/a)') 
     * '*** Error (icescenario.F):',
     * ' no recognized scenario switch (-D) with RCMMATRIX, RCMMATCO2'
      stop

crrrrr
#endif
crrrrr

      return
      end

c-----------------------------------------------------------------------

cccccc
#endif
cccccc
