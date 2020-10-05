c comicepath.h

c     For Suns, must be in common (no substrings allowed with params).
c     For linux, must be parameters (multiply-defined errors with vars).

#if defined (SUN)
      common /cpathnamec/ cpathhome, cpathbed, 
     *                    cpathdrive, cpathdrivercm, cpathdriveobs,
     *                    cpathmap, cpathslgrid
#endif
      character*80 cpathhome, cpathbed, 
     *             cpathdrive, cpathdrivercm, cpathdriveobs,
     *             cpathmap, cpathslgrid

c-------------------------------------------------------------------
c     for empirical data files (time series, Ant precip, orbits...):
c-------------------------------------------------------------------
#if defined (SUN)
      data cpathhome /
#else
      parameter (cpathhome =
#endif

#if defined (RWRAP_PATH)
     *  RWRAP_PATH // '/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/'
c    *  '/storage/work/dxp21/Data_ice/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Sheetshelf_data/'

#elif defined (TERRA)
     *  '/home/deconto/Sheetshelf/'

#elif defined (MAC)
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/'           ! natalya

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Sheetshelf/'

#elif defined (LOKI)
     *  '/home/trhille/PSUICE/Data_ice/'
c    *  '/home/pspec/PSUICE/Data_ice/'

#elif defined (AMAREL)
     *  '/home/dp980/Sheetshelf/Data_ice/'

#elif defined (BRIDGES)
     *  '/pylon5/mc3bggp/vhayot/Sheetshelf/Data_ice/'
     
#elif defined (COMET)
     *  '/home/vhayot/' //
     *  'amarel_ism_version/' //
     *  'Sheetshelf/Data_ice/'
     
#else
     *  '../'
#endif

#if defined (SUN)
     *  /
#else
     *  )
#endif

c----------------------------
c     for BEDMAP1 data files:
c----------------------------
#if defined (SUN)
      data cpathbed /
#else
      parameter (cpathbed =
#endif

#if defined (RWRAP_PATH)
     *  RWRAP_PATH // '/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/'
c    *  '/storage/work/dxp21/Data_ice/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Data_bedmap1/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Bedmap/'

#elif defined (TERRA)
     *  '/home/deconto/Bedmap/'

#elif defined (MAC)
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Bedmap/'    ! natalya

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Bedmap/'

#elif defined (BRIDGES)
     *  '/pylon5/mc3bggp/vhayot/Sheetshelf/Data_ice/'
     
#elif defined (COMET)
     *  '/oasis/scratch/comet/vhayot/temp_project/SAMPLE-ICM
     *  /4amarel_ism_version/Sheetshelf/Data_ice/'

#else
     *  '../Bedmap1/'
#endif

#if defined (SUN)
     *  /
#else
     *  )
#endif

c---------------------------
c     for GCM driving files:
c---------------------------
#if defined (SUN)
      data cpathdrive /
#else
      parameter (cpathdrive =
#endif

#if defined (RWRAP_PATH)
     *  RWRAP_PATH // '/Driveice/Gcm/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/Driveice/Gcm/'
c    *  '/storage/work/dxp21/Data_ice/Driveice/Gcm/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Driveice/Gcm/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Driveice/Gcm/'

#elif defined (TERRA)
     *  '/home/deconto/Driveice/Gcm/'

#elif defined (MAC)
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Gcm/'       ! natalya

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Driveice/Gcm/'

#elif defined (BRIDGES)
     *  '/pylon5/mc3bggp/vhayot/Sheetshelf/Data_ice/Gcm/'

#elif defined (COMET)
     *  '/oasis/scratch/comet/vhayot/temp_project/SAMPLE-ICM
     *  /3amarel_ism_version/Sheetshelf/Data_ice/Gcm/'

#else
     *  '../Driveice/Gcm/'
#endif

#if defined (SUN)
     *  /
#else
     *  )
#endif

c---------------------------
c     for OBS driving files:
c---------------------------
#if defined (SUN)
      data cpathdriveobs /
#else
      parameter (cpathdriveobs =
#endif

#if defined (RWRAP_PATH)
     *  RWRAP_PATH // '/Driveice/Obs/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/Driveice/Obs/'
c    *  '/storage/work/dxp21/Data_ice/Driveice/Obs/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Driveice/Obs/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Obs/'

#elif defined (TERRA)
     *  '/home/deconto/Driveice/Obs/'

#elif defined (MAC)
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Obs/'       ! natalya

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Driveice/Obs/'
 
#elif defined (BRIDGES)
     *  '/pylon5/mc3bggp/vhayot/Sheetshelf/Data_ice/Driveice/Obs/'

#elif defined (COMET)
     *  '/oasis/scratch/comet/vhayot/temp_project/SAMPLE-ICM
     *  /1amarel_ism_version/Sheetshelf/Data_ice/Driveice/Obs/'

#else
     *  '../Driveice/Obs/'
#endif

#if defined (SUN)
     *  /
#else
     *  )
#endif

c---------------------------
c     for RCM driving files:
c---------------------------
#if defined (SUN)
      data cpathdrivercm /
#else
      parameter (cpathdrivercm =
#endif

#if defined (RWRAP_PATH)
     *  RWRAP_PATH // '/Driveice/Rcm/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/Driveice/Rcm/'
c    *  '/storage/work/dxp21/Data_ice/Driveice/Rcm/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Driveice/Rcm/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Rcm/'

#elif defined (TERRA)
     *  '/home/deconto/Driveice/Rcm/'

#elif defined (MAC)
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Rcm/'       ! natalya

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Driveice/Rcm/'

#elif defined (BRIDGE)
     *  '/pylon5/mc3bggp/vhayot/Sheetshelf/Data_ice/Driveice/Rcm'

#elif defined (COMET)
     *  '/oasis/scratch/comet/vhayot/temp_project/SAMPLE-ICM
     *  /2amarel_ism_version/Sheetshelf/Data_ice/Driveice/Rcm'

#else
     *  '../Driveice/Rcm/'
#endif

#if defined (SUN)
     *  /
#else
     *  )
#endif

c---------------------------------------
c     for Genesis global map data files:
c---------------------------------------
#if defined (SUN)
      data cpathmap /
#else
      parameter (cpathmap =
#endif

#if defined (RWRAP_PATH)
     *  RWRAP_PATH // '/Genesis.3.0/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/Genesis.3.0/'
     *  // 'Datafiles/'
c    *  '/storage/work/dxp21/Data_ice/Genesis.3.0/Datafiles/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Genesis.3.0/Datafiles/'

#elif defined (ACI)
     *  '/storage/home/dxp21/Genesis.3.0/Datafiles/'

#elif defined (TERRA)
     *  '/home/deconto/Datafiles/'

#elif defined (MAC)
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Genesis.3.0/Datafiles/'

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Genesis.3.0/Datafiles/'

#elif defined (BRIDGES)
     *  '/pylon5/mc3bggp/vhayot/Sheetshelf/Data_ice/
     *  Genesis.3.0/Datafiles'

#elif defined (COMET)
     *  '/home/vhayot/amarel_ism_version/Sheetshelf/Data_ice/
     *  Genesis.3.0/Datafiles'

#else
     *  '../Genesis.3.0/Datafiles/'
#endif

#if defined (SUN)
     *  /
#else
     *  )
#endif

c------------------------------------------------------------
c     for Gomez-Mitrovica bedrock model files (permanent area, 
c     only for Gridfiles, in bedrock_sl). 
c------------------------------------------------------------

#if defined (SUN)
      data cpathslgrid /
#else
      parameter (cpathslgrid =
#endif


#if defined (RWRAP_PATH)
     *  '/iguana/s1/pollard/Gomez/SEANEW/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Gomez/SEANEW/'

#elif defined (ACI)
     *  '/storage/home/dxp21/Gomez/SEANEW/'

#elif defined (MAC)
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/SEANEW/'

#else

     *  '../Gomez/SEANEW/'
#endif

#if defined (SUN)
     *  /
#else
     *  )
#endif

