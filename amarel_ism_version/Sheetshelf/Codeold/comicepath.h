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
     *  '/home/dp980/Sheetshelf/Data_ice/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/'
c    *  '/storage/work/dxp21/Data_ice/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/'

#elif defined (LION)
     *  '/gpfs/work/dxp21/Sheetshelf_data/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Sheetshelf_data/'

#elif defined (LUNA)
#  if defined (DECONTO)
     *  '/home/data01/deconto/Sheetshelf/'
#  else
     *  '/home/data01/pollard/Sheetshelf/'
#  endif

#elif defined (TERRA)
#  if defined (DECONTO)
     *  '/home/deconto/Sheetshelf/'
#  else
     *  '/home/pollard/Sheetshelf/'
#  endif

#elif defined (KUROSHIO)
     *  '/home/malteh/models/loveclip/psuim/'

#elif defined (CONDOR)
     *  '/condor/data2/pollard/Sheetshelf/'

#elif defined (MAC)
c    *  '/Users/natalya/Desktop/Pollard_Ice_Model/'
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/'

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Sheetshelf/'

#elif defined (ENDURANCE)
     *  '/endurance/s0/aum34/modeling/'

#elif defined (SPROUT)
     *  '/home/truffle/Sheetshelf/'

#elif defined (LOKI)
     *  '/home/trhille/PSUICE/Data_ice/'
c    *  '/home/pspec/PSUICE/Data_ice/'

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
     *  '/home/dp980/Sheetshelf/Data_ice/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/'
c    *  '/storage/work/dxp21/Data_ice/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Data_bedmap1/'

#elif defined (LION)
     *  '/gpfs/work/dxp21/Bedmap/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Bedmap/'

#elif defined (LUNA)
#  if defined (DECONTO)
     *  '/home/data01/deconto/Bedmap/'
#  else
     *  '/home/data01/pollard/Bedmap/'
#  endif

#elif defined (TERRA)
#  if defined (DECONTO)
     *  '/home/deconto/Bedmap/'
#  else
     *  '/home/pollard/Bedmap/'
#  endif

#elif defined (CONDOR)
     *  '/condor/data2/pollard/Bedmap/'

#elif defined (MAC)
c    *  '/Users/natalya/Desktop/Pollard_Ice_Model/Bedmap/'
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Bedmap/'

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Bedmap/'

#elif defined (ENDURANCE)
     *  '/endurance/s0/aum34/modeling/Bedmap'

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
     *  '/home/dp980/Sheetshelf/Data_ice/Driveice/Gcm/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/Driveice/Gcm/'
c    *  '/storage/work/dxp21/Data_ice/Driveice/Gcm/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Driveice/Gcm/'

#elif defined (LION)
     *  '/gpfs/work/dxp21/Driveice/Gcm/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Driveice/Gcm/'

#elif defined (LUNA)
#  if defined (DECONTO)
     *  '/home/data01/deconto/Driveice/Gcm/'
#  else
     *  '/home/data01/pollard/Driveice/Gcm/'
#  endif

#elif defined (TERRA)
#  if defined (DECONTO)
     *  '/home/deconto/Driveice/Gcm/'
#  else
     *  '/home/pollard/Driveice/Gcm/'
#  endif

#elif defined (CONDOR)
     *  '/condor/data2/pollard/Driveice/Gcm/'

#elif defined (MAC)
c    *  '/Users/natalya/Desktop/Pollard_Ice_Model/Gcm/'
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Gcm/'

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Driveice/Gcm/'

#elif defined (ENDURANCE)
     *  '/endurance/s0/aum34/modeling/Driveice/Gcm'

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
     *  '/home/dp980/Sheetshelf/Data_ice/Driveice/Obs/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/Driveice/Obs/'
c    *  '/storage/work/dxp21/Data_ice/Driveice/Obs/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Driveice/Obs/'

#elif defined (LION)
     *  '/gpfs/work/dxp21/Driveice/Obs/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Obs/'

#elif defined (LUNA)
#  if defined (DECONTO)
     *  '/home/data01/deconto/Driveice/Obs/'
#  else
     *  '/home/data01/pollard/Driveice/Obs/'
#  endif

#elif defined (TERRA)
#  if defined (DECONTO)
     *  '/home/deconto/Driveice/Obs/'
#  else
     *  '/home/pollard/Driveice/Obs/'
#  endif

#elif defined (CONDOR)
     *  '/condor/data2/pollard/Driveice/Obs/'

#elif defined (MAC)
c    *  '/Users/natalya/Desktop/Pollard_Ice_Model/Obs/'
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Obs/'

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Driveice/Obs/'

#elif defined (ENDURANCE)
     *  '/endurance/s0/aum34/modeling/Driveice/Obs'

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
     *  '/home/dp980/Sheetshelf/Data_ice/Driveice/Rcm/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/Driveice/Rcm/'
c    *  '/storage/work/dxp21/Data_ice/Driveice/Rcm/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Ebm2d/Codesheetshelf/Driveice/Rcm/'

#elif defined (LION)
     *  '/gpfs/work/dxp21/Driveice/Rcm/'

#elif defined (ACI)
     *  '/storage/home/dxp21/work/Rcm/'

#elif defined (LUNA)
#  if defined (DECONTO)
     *  '/home/data01/deconto/Driveice/Rcm/'
#  else
     *  '/home/data01/pollard/Driveice/Rcm/'
#  endif

#elif defined (TERRA)
#  if defined (DECONTO)
     *  '/home/deconto/Driveice/Rcm/'
#  else
     *  '/home/pollard/Driveice/Rcm/'
#  endif

#elif defined (CONDOR)
     *  '/condor/data2/pollard/Driveice/Rcm/'

#elif defined (MAC)
c    *  '/Users/natalya/Desktop/Pollard_Ice_Model/Rcm/'
     *  '/home/ng50/coupled_code/Pollard_Ice_Model/Rcm/'

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Driveice/Rcm/'

#elif defined (ENDURANCE)
     *  '/endurance/s0/aum34/modeling/Driveice/Rcm'

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
     *  '/home/dp980/Sheetshelf/Data_ice/Genesis.3.0/'
c    *  '/iguana/s1/pollard/People/Fuller_lee/Data_ice/Genesis.3.0/'
c    *  // 'Datafiles/'
c    *  '/storage/work/dxp21/Data_ice/Genesis.3.0/Datafiles/''

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Genesis.3.0/Datafiles/'

#elif defined (LION)
     *  '/gpfs/home/dxp21/Genesis.3.0/Datafiles/'

#elif defined (ACI)
     *  '/storage/home/dxp21/Genesis.3.0/Datafiles/'

#elif defined (LUNA)
#  if defined (DECONTO)
     *  '/home/data01/deconto/Genesis.3.0/Datafiles/'
#  else
     *  '/home/data01/pollard/Genesis.3.0/Datafiles/'
#  endif

#elif defined (TERRA)
#  if defined (DECONTO)
     *  '/home/deconto/Genesis.3.0/Datafiles/'
#  else
     *  '/home/pollard/Genesis.3.0/Datafiles/'
#  endif

#elif defined (KUROSHIO)
     *  '/home/malteh/models/loveclip/psuim/Genesis.3.0/Datafiles/'

#elif defined (CONDOR)
     *  '/condor/data2/pollard/Genesis.3.0/Datafiles/'

#elif defined (MAC)
c    *'/Users/natalya/Desktop/Pollard_Ice_Model/Genesis.3.0/Datafiles/'
     *'/home/ng50/coupled_code/Pollard_Ice_Model/Genesis.3.0/Datafiles/'

#elif defined (WOJU)
     *  '/home/scrim/dxp21/Genesis.3.0/Datafiles/'

#elif defined (ENDURANCE)
     *  '/endurance/s0/aum34/modeling/Genesis.3.0/Datafiles'

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
     *  '/home/dp980/Gomez/SEANEW/'
c    *  '/iguana/s1/pollard/Gomez/SEANEW/'

#elif defined (SUN) || defined (IGUANA)
     *  '/iguana/s1/pollard/Gomez/SEANEW/'

#elif defined (LION)
     *  '/gpfs/home/dxp21/Gomez/SEANEW/'

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

