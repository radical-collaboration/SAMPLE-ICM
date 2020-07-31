c***********************************************************************
#if defined (GCMMATRIX) || defined (RCMMATRIX) || defined (RACMOCLIM) || defined (CO2INTER) || defined (LOVECLIP) || defined (CSMCLIM) || defined (CLEMCLIM) || defined (CSMORDO)
c***********************************************************************

c-----------------------------------------------------------------------

      subroutine climate_gcm (tog, tmon, pmon, h, hs, hsbal,
#if defined (LOVECLIP)
     *                        emon, qmon,
#endif
#if defined (CO2INTER)
     *                        toge, tage, prge, ruge, unofa,
#elif defined (RACMOCLIM)
     *                        ruge, runofa,      
#endif
     *                        timeice, rco2, ecc, obl, prec, 
     *                        facice, facorb, facco2)

c     Sets monthly surface climate (tmon,pmon) for [GCM,RCM]MATRIX, 
c     on ice grid, weighted according to facice, facorb, facco2 
c     (fac_matrix, domatrix) from matrix of stored Genesis GCM files 
c     (generated previously by debplotdriveice). GCM data files are 
c     read just once at start, into large arrays *drive(..,..,nmon,nmat)
c     by readdrive. 
c
c     Nb: tmon,pmon are *not* lapse-rate corrected (done in climate_all)

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

c     passed:
      dimension
     *  tog(nx,ny),  
     *  tmon(nmon,nx,ny),  pmon(nmon,nx,ny),
#if defined (LOVECLIP)
     *  emon(nmon,nx,ny),  qmon(nmon,nx,ny),
#endif
     *  h(nx,ny),          hs(nx,ny),      hsbal(nx,ny)
#if defined (CO2INTER)
     * ,toge(nlong,nlatg),
     *  tage(nmon,nlong,nlatg), prge(nmon,nlong,nlatg),
     *  ruge(nlong,nlatg),      runofa(nx,ny)
#elif defined (RACMOCLIM)    
     * ,ruge(nlong,nlatg),      runofa(nx,ny)
#endif

c     agcm surface forcing arrays (read from files by debplotdriveice):
c     to* = spectral topography (m)
c     pr* = precipitation (mm/day)
c     ta* = air temperature at lowest agcm layer (C)
#if ! defined (CO2INTER)
      dimension
     *  toge(nlong,nlatg),
     *  tage(nmon,nlong,nlatg),
     *  prge(nmon,nlong,nlatg)
      save toge, tage, prge
#  if defined (LOVECLIP)
      dimension
     *  evge(nmon,nlong,nlatg),
     *  qsge(nmon,nlong,nlatg)
      save evge, qsge
#  endif
#endif

      dimension
     *  todrive(nlong,nlatg,   1,nmat),
     *  tadrive(nlong,nlatg,nmon,nmat),
     *  prdrive(nlong,nlatg,nmon,nmat),
     *  tag(nlong,nlatg), 
     *  prg(nlong,nlatg)
      save todrive, tadrive, prdrive
#if defined (LOVECLIP)
      dimension
     *  evdrive(nlong,nlatg,nmon,nmat),
     *  qsdrive(nlong,nlatg,nmon,nmat)
      save evdrive, qsdrive
#endif
#if defined (CO2INTER) || defined (RACMOCLIM)
      dimension
     *  rudrive(nlong,nlatg,   1,nmat)
      save rudrive
#endif

#if ! defined (RACMOCLIM)
c     local 2-D ice-grid arrays, indexes and weights to GCM grid:
      dimension
     *  indl(nx,ny),     weil(nx,ny),
     *  indb(nx,ny),     weib(nx,ny)
      save indl, weil, indb, weib
#endif

#if defined (NHA) && ! defined (NHAGLAC) && ! defined (LOVECLIP)
      dimension albeff(nx,ny)
#endif

c     for RCM or Loveclip anomaly method (observed - control):
#if defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM)
      dimension
     *  togobs(nx,ny), tmonobs(nmon,nx,ny), pmonobs(nmon,nx,ny),
     *  todrivectl(nlong,nlatg),         togectl(nlong,nlatg),
     *  tadrivectl(nlong,nlatg,nmon),    tagectl(nmon,nlong,nlatg),
     *  prdrivectl(nlong,nlatg,nmon),    prgectl(nmon,nlong,nlatg),
     *  togctl(nx,ny),
     *  tmonctl(nmon,nx,ny),
     *  pmonctl(nmon,nx,ny)
      save 
     *  togobs, tmonobs, pmonobs, togctl, tmonctl, pmonctl
# if defined (LOVEANOM)   
      dimension
     *  evdrivectl(nlong,nlatg,nmon),    evgectl(nmon,nlong,nlatg),
     *  emonctl(nmon,nx,ny)
      save emonctl
#  endif
      logical firstanom
      data firstanom /.true./
      save firstanom
#endif

      dimension ndaypm(12)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm

      logical first
      save first
      data first /.true./

      character*200 cpath(nmat), cfile, cname, cline,
     *              cfilez, cfilet, cfilep
      save cpath

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if defined (EISANTA) && ! defined (ASYNCH) && ! defined (RCMMATRIX) && ! defined (RACMOCLIM) && ! defined (LOVECLIP) && ! defined (CSMCLIM) && ! defined (CLEMCLIM) && ! defined (CSMORDO)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*20 csuffant(nmat)
      data csuffant / 
c 9x2 matrix (CO2 2x,3x):
c    *                '2noco', '2nome', '2nowa',
c    *                '2meco', '2meme', '2mewa',
c    *                '2fuco', '2fume', '2fuwa',
c    *                '3noco', '3nome', '3nowa',
c    *                '3meco', '3meme', '3mewa',
c    *                '3fuco', '3fume', '3fuwa' 

c 9x2 matrix (CO2 4x,8x):
c    *                '4nome', '4nome', '4nowa',
c    *                '4meme', '4meme', '4mewa',
c    *                '4fume', '4fume', '4fuwa',
c    *                '8nome', '8nome', '8nowa',
c    *                '8meme', '8meme', '8mewa',
c    *                '8fume', '8fume', '8fuwa' 

c 9x3 matrix (CO2 1x,4x,8x):
c    *                '1nome', '1nome', '1nome',
c    *                '1meme', '1meme', '1meme',
c    *                '1fume', '1fume', '1fume',
c    *                '4nome', '4nome', '4nome',
c    *                '4meme', '4meme', '4meme',
c    *                '4fume', '4fume', '4fume',
c    *                '8nome', '8nome', '8nome',
c    *                '8meme', '8meme', '8meme',
c    *                '8fume', '8fume', '8fume'

c 9x2 matrix (CO2 2/4, q/_/r, no/me/fu, all me orbit) (for CO2INTER):
c    *                '2nomes', '2nome', '2nomer',
c    *                '2memes', '2meme', '2memer',
c    *                '2fumes', '2fume', '2fumer',
c    *                '4nomes', '4nome', '4nomer',
c    *                '4memes', '4meme', '4memer',
c    *                '4fumes', '4fume', '4fumer'

c 9x4 matrix (CO2 0.5x,1x,2x,4x or 1x,2x,4x,8x) (see c05 in domatrix):
c05  *                '0.5nome', '0.5nome', '0.5nome',   ! c05
c05  *                '0.5meme', '0.5meme', '0.5meme',   ! c05
c05  *                '0.5fume', '0.5fume', '0.5fume',   ! c05
     *                '1noco', '1nome', '1nowa',
     *                '1meco', '1meme', '1mewa',
     *                '1fuco', '1fume', '1fuwa',
     *                '2noco', '2nome', '2nowa',
     *                '2meco', '2meme', '2mewa',
     *                '2fuco', '2fume', '2fuwa',
     *                '4noco', '4nome', '4nowa',
     *                '4meco', '4meme', '4mewa',
     *                '4fuco', '4fume', '4fuwa',
     *                '8noco', '8nome', '8nowa',         ! c05
     *                '8meco', '8meme', '8mewa',         ! c05
     *                '8fuco', '8fume', '8fuwa'          ! c05

c mis31 gcm (1 = 1.070 Ma, 2 = 1.082 Ma, see domatrix:
c    *                  'Mis31g' , 'Mis31f' 

     *             /
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#elif defined (NHA) && defined (NHAGLAC)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*20 csuffant(nmat)
      data csuffant / 
     *                'Noco',  'Nowa',
     *                '11kco', '11kwa',
     *                'Fuco',  'Fuwa',
     *                'Obs'
     *              /
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#elif defined (NHA) && ! defined (NHAOBS) && ! defined (NHAGLAC) && ! defined (LOVECLIP) && ! defined (CSMCLIM) && ! defined (CLEMCLIM) && ! defined (CSMORDO)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*20 csuffant(nmat)
      data csuffant / 
     *                'Modco', 'Modme', 'Modwa'
     *              /
c~~~~~
#endif
c~~~~~

      dimension alatt31(48)
      data alatt31 /
     *  -1.5212131731, -1.4569823009, -1.3923721142, -1.3276762287,
     *  -1.2629478409, -1.1982037311, -1.1334508412, -1.0686925583,
     *  -1.0039307308, -0.9391664533, -0.8744004155, -0.8096330743,
     *  -0.7448647446, -0.6800956510, -0.6153259588, -0.5505557923,
     *  -0.4857852479, -0.4210144016, -0.3562433151, -0.2914720392,
     *  -0.2267006173, -0.1619290868, -0.0971574813, -0.0323858319,
     *   0.0323858319,  0.0971574813,  0.1619290868,  0.2267006173,
     *   0.2914720392,  0.3562433151,  0.4210144016,  0.4857852479,
     *   0.5505557923,  0.6153259588,  0.6800956510,  0.7448647446,
     *   0.8096330743,  0.8744004155,  0.9391664533,  1.0039307308,
     *   1.0686925583,  1.1334508412,  1.1982037311,  1.2629478409,
     *   1.3276762287,  1.3923721142,  1.4569823009,  1.5212131731
     *   /

c        Initialization section for gcm forcing fields

c====================
      if (first) then
c====================

#if defined (ASYNCH)
        do n=1,nmat
          cpath(n) = './'
        enddo
#elif defined (GCMMATRIX) || defined (CO2INTER)
#  if defined (NHAOBS)
        do n=1,nmat
          cpath(n) = cpathdriveobs(1:lenchr(cpathdriveobs)) 
        enddo
#  elif defined (NHAGLAC)
        do n=1,nmat
          cpath(n) = cpathdrive(1:lenchr(cpathdrive))
     *               // 'gen2/'
     *               // csuffant(n)(1:lenchr(csuffant(n))) // '/'
        enddo
#  else
        do n=1,nmat
          cpath(n) = cpathdrive(1:lenchr(cpathdrive))
     *               // 'gen3/'
     *               // csuffant(n)(1:lenchr(csuffant(n))) // '/'
        enddo
#  endif
#elif defined (RCMMATRIX) 
        do n=1,nmat
          cpath(n) = cpathdrivercm
        enddo
#elif defined (RACMOCLIM) || defined (CSMCLIM) || defined (CLEMCLIM) || defined (CSMORDO)
        do n=1,nmat
          cpath(n) = cpathhome
        enddo
#endif

c          Set gcm axes (along,alatg in comicegrid). For LOVECLIP,
c          CSMCLIM, CLEMCLIM, CSMORDO, along,alatg are set in first 
c          call to readloveclim, readdrive_csm, readdrive_clem,
c          readdrive_csmordo, resp.

#if ( defined (GCMMATRIX) || defined (CO2INTER) ) && ! defined(CSMCLIM) && ! defined (CLEMCLIM) && ! defined (CSMORDO)
        do i=1,nlong
c         cut at dateline: fields on gcm files have been shifted 180 deg
          along(i) = -pi + 2.*pi*(i-0.5)/nlong
        enddo
        if (nlatg.eq.48) then
          do j=1,nlatg
            alatg(j) = alatt31(j)
          enddo
        else
          do j=1,nlatg
            alatg(j) = -0.5*pi + pi*(j-0.5)/nlatg
          enddo
        endif
#endif

c          Open and read gcm forcing files into *drive arrays

c------------------
        do m=1,nmat
c------------------

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#if defined (GCMMATRIX) || defined (CO2INTER)
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#  if defined (NHA)
#    if defined (NHAOBS)
          cfile = cpath(m)(1:lenchr(cpath(m)))
     *          // 'data_icedrive_TO_OBS'                          ! T31
#    elif defined (NHAGLAC)
          cfile = cpath(m)(1:lenchr(cpath(m)))
     *          // 'data_icedrive_TO'                              ! T31
#    else
          cfile = cpathmap(1:lenchr(cpathmap))
     *          // 'V2/data_top_2X2'
#    endif
#  elif defined (ASYNCH) 
          cfile = 'data_asyn_top_2X2'
#  else
c         read from TOPOG2 files (to get exactly corresponding topog 
c         if matrix files were generated using TA and TOPOGUN. TOPOG2
c         was used in LSX for delevtype=1. TOPOGUN does not correspond 
c         exactly to LSX correction (it's the T31 aggregate, re-interped
c         back to 2X2).
          mtop = mod((m-1)/3,3)
          if (mtop.eq.0) then
            cfile = cpathmap(1:lenchr(cpathmap))
     *            // '34ma/data_35ma_non_top_2X2'
          else if (mtop.eq.1) then
            cfile = cpathmap(1:lenchr(cpathmap))
     *            // '34ma/data_35ma_med_top_2X2'
          else if (mtop.eq.2) then
            cfile = cpathmap(1:lenchr(cpathmap))
     *            // '34ma/data_35ma_full_top_2X2'
          endif
#  endif
#  if defined (NHA) && (defined (NHAOBS) || defined (NHAGLAC))
          call readdrive (iuto, cfile, todrive(1,1,1,m), 'TO', 1)
#  else
          open (iuto, file=cfile, form='formatted', status='old')
          call readdata_ice (iuto, todrive(1,1,1,m), nlong,nlatg,
     *                       nlong,nlatg, 'TO      ',  2, 1)
          close (iuto)
#  endif

c^^^^^^^^^^^^^^^^^^^^^^^^
#elif defined (RCMMATRIX)
c^^^^^^^^^^^^^^^^^^^^^^^^

#  if defined (RCMMATCO2)

          if (m.eq.1) then
            cname = 'Runpre_40/LSX.nc'
          endif

#    if defined (RAMPPLIO) || defined (RAMP2X) || defined (RAMP4X) || defined (RAMP8X) || defined (RAMP8XW)
          if (m.eq.2) then
#      if defined (RAMPPLIO)
            cname = 'Runmis31_400w_40/LSX_mis31_400w_40.nc'       ! plio
#      elif defined (RAMP2X)
            cname = 'Run2x/LSX.nc'                                ! 2x
#      elif defined (RAMP4X)
            cname = 'Run4x/LSX.nc'                                ! 4x
#      elif defined (RAMP8X)
            cname = 'Run8x/LSX.nc'                                ! 8x
#      elif defined (RAMP8XW)
            cname = 'Run8xw/LSX.nc'                               ! 8xw
#      endif
          endif
#    else
          if (m.eq.2) then
            cname = 'Run2x/LSX.nc'
          else if (m.eq.3) then
            cname = 'Run4x/LSX.nc'
c           cname = 'Run4xw/LSX.nc'
          else if (m.eq.4) then
            cname = 'Run8x/LSX.nc'
c           cname = 'Run8xw/LSX.nc'
          endif
#    endif

#  elif defined (RCMFILE)

          cname = RCMFILE

#  else

          cname = 'Runpre_40/LSX.nc'

#  endif

          cfile = cpath(m)(1:lenchr(cpath(m))) // cname
          call readdrive_rcm (iuto,cfile,todrive(1,1,1,m),'TOPOGNOD',1)

c^^^^^^^^^^^^^^^^^^^^^^^^
#elif defined (RACMOCLIM)
c^^^^^^^^^^^^^^^^^^^^^^^^
c         all done by calling getracmo_[ant,gre]_fields below

c^^^^^^^^^^^^^^^^^^^^^^^
#elif defined (LOVECLIP)
c^^^^^^^^^^^^^^^^^^^^^^^
          call readloveclim (todrive(1,1,1,m), 'surfelev', 1)  

c^^^^^^^^^^^^^^^^^^^^^^
#elif defined (CSMCLIM)
c^^^^^^^^^^^^^^^^^^^^^^
          cname = 
#  if defined (CSMEXP)
     *      'Data_csm/' // CSMEXP // '/' //
#  else
     *      'Data_csm/' //
#  endif
     *      'inputdata_Topog.nc'                               ! J. Tsai
          cfile = cpath(m)(1:lenchr(cpath(m))) // cname
          call readdrive_csm (iuto,cfile,todrive(1,1,1,m),'PHIS',1,
     *                        0, timeice)

c^^^^^^^^^^^^^^^^^^^^^^
#elif defined (CLEMCLIM)
c^^^^^^^^^^^^^^^^^^^^^^
          cname = 
     *      'Data_clemens/' //
     *      'MPI-EPS-LR_RCP45TempPrecOro2200Monthly.nc'
          cfile = cpath(m)(1:lenchr(cpath(m))) // cname
          call readdrive_clem (iuto,cfile,todrive(1,1,1,m),'orog',1,
     *                         0, timeice)

c^^^^^^^^^^^^^^^^^^^^^^
#elif defined (CSMORDO)
c^^^^^^^^^^^^^^^^^^^^^^
c         path to files only:
          cfile = cpath(m)(1:lenchr(cpath(m))) // 'Data_csmordo/' 
          call readdrive_csmordo (iuto,cfile,todrive(1,1,1,m),'PHIS',1,
     *                            timeice)

c^^^^^
#endif
c^^^^^

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if defined (GCMMATRIX) || defined (CO2INTER)  
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  if defined (NHAOBS)
          cfile = cpath(m)(1:lenchr(cpath(m)))//'data_icedrive_PR_SHEA'
          call readdrive (iupr, cfile, prdrive(1,1,1,m), 'PR', nmon)

          cfile = cpath(m)(1:lenchr(cpath(m)))//'data_icedrive_TA_CRUCD'
          call readdrive (iuta, cfile, tadrive(1,1,1,m), 'TA', nmon)
#  else
          cfile = cpath(m)(1:lenchr(cpath(m)))//'data_icedrive_PR'
          call readdrive (iupr, cfile, prdrive(1,1,1,m), 'PR', nmon)

          cfile = cpath(m)(1:lenchr(cpath(m)))//'data_icedrive_TA'
          call readdrive (iuta, cfile, tadrive(1,1,1,m), 'TA', nmon)
#  endif
#  if defined (CO2INTER)
          cfile = cpath(m)(1:lenchr(cpath(m)))//'data_icedrive_RU'
          call readdrive (iuru, cfile, rudrive(1,1,1,m), 'RU',    1)
#  endif
c%%%%%%%%%%%%%%%%%%%%%%%%
#elif defined (RCMMATRIX)
c%%%%%%%%%%%%%%%%%%%%%%%%
c         cfile has been set above for TOPOGNOD:
          call readdrive_rcm (iuto,cfile,prdrive(1,1,1,m),'PRECIP',nmon)
          call readdrive_rcm (iuto,cfile,tadrive(1,1,1,m),'TS2',   nmon)
c%%%%%%%%%%%%%%%%%%%%%%%%
#elif defined (RACMOCLIM)
c%%%%%%%%%%%%%%%%%%%%%%%%
c         all done by call to getracmo_[ant,gre]_fields below
c%%%%%%%%%%%%%%%%%%%%%%%
#elif defined (LOVECLIP)
c%%%%%%%%%%%%%%%%%%%%%%%
          call readloveclim (tadrive(1,1,1,m), 'tair2m',    nmon)  
          call readloveclim (prdrive(1,1,1,m), 'precip',    nmon)  
          call readloveclim (evdrive(1,1,1,m), 'evap',      nmon)  
          call readloveclim (qsdrive(1,1,1,m), 'solarflux', nmon)  
c%%%%%%%%%%%%%%%%%%%%%%
#elif defined (CSMCLIM)
c%%%%%%%%%%%%%%%%%%%%%%
c J. Tsai: precip and temp read below, every call. 
c          (topog already read above, for lon,lat, needed for setinterp)
c%%%%%%%%%%%%%%%%%%%%%%
#elif defined (CLEMCLIM)
c%%%%%%%%%%%%%%%%%%%%%%
c precip and temp read below, every call. 
c (topog already read above, for lon,lat, needed for setinterp)
c%%%%%%%%%%%%%%%%%%%%%%%
#elif defined (CSMORDO)
c%%%%%%%%%%%%%%%%%%%%%%%
c        cfile has been set above (path to files only):
         call readdrive_csmordo (iuto,cfile,tadrive(1,1,1,m),'Tlowest',
     *                           nmon, timeice)
         call readdrive_csmordo (iuto,cfile,prdrive(1,1,1,m),'PRECT',
     *                           nmon, timeice)
c%%%%%
#endif
c%%%%%

c------------
        enddo
c------------

c          Set indices for bilinear interp (ind*,wei* local,saved)

#if defined (GCMMATRIX) || defined (CO2INTER)
        call setinterp (along, alatg, nlong, nlatg, 
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, rotate_to_gcm)
#elif defined (LOVECLIP)
        call setinterp (along, alatg, nlong, nlatg, 
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, 0.)
#elif defined (RCMMATRIX)
        call setinterp_rcm (nlong, nlatg, dd0rcm,
     *                      xh, yh, nx, ny,
     *                      indl, weil, indb, weib)
#elif defined (RACMOCLIM)
c       all done by calling getracmo_[ant,gre]_fields below
#elif defined (CSMCLIM)
        call setinterp (along, alatg, nlong, nlatg, 
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, 0.)
#elif defined (CLEMCLIM)
        call setinterp (along, alatg, nlong, nlatg, 
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, 0.)
#elif defined (CSMORDO)
        call setinterp (along, alatg, nlong, nlatg, 
     *                  alond, alatd, nx, ny,
     *                  indl, weil, indb, weib, 0.)
#endif

        first = .false.

c==========
      endif  ! end of initialization section
c==========

c%%%%%%%%%%%%%%%%%%%%
#if defined (CSMCLIM)
c%%%%%%%%%%%%%%%%%%%%
c  J. Tsai: read data files every call (once per year, dtimeclim=1)

      do m=1,nmat

c        todrive already read in, saved above - assume unchanged.
c        cname = 
c#  if defined (CSMEXP)
c     *      'Data_csm/' // CSMEXP // '/' //
c#  else
c     *      'Data_csm/' //
c#  endif
c     *      'inputdata_Topog.nc'                              ! J. Tsai
c        cfile = cpath(m)(1:lenchr(cpath(m))) // cname
c        call readdrive_csm (iuto,cfile,todrive(1,1,1,m),'PHIS',1,
c     *                      0, timeice)

        cname = 
#  if defined (CSMEXP)
     *      'Data_csm/' // CSMEXP // '/' //
#  else
     *      'Data_csm/' //
#  endif
     *      'inputdata_PRECP.nc'                               ! J. Tsai
        cfile = cpath(m)(1:lenchr(cpath(m))) // cname
        call readdrive_csm (iuto,cfile,prdrive(1,1,1,m),'PRECT',nmon,
     *                      0, timeice)

        cname = 
#  if defined (CSMEXP)
     *      'Data_csm/' // CSMEXP // '/' //
#  else
     *      'Data_csm/' //
#  endif
     *      'inputdata_2mSAT.nc'                               ! J. Tsai
        cfile = cpath(m)(1:lenchr(cpath(m))) // cname
        call readdrive_csm (iuto,cfile,tadrive(1,1,1,m),'T',  nmon,
     *                      0, timeice)
      enddo
c%%%%%
#endif
c%%%%%

c|||||||||||||||||||||
#if defined (CLEMCLIM)
c|||||||||||||||||||||
c  Clemens Schannwell: read data files every call (once per year, 
c  dtimeclim=1)

      do m=1,nmat

c        todrive already read in, saved above - assume unchanged.

        cname = 
     *      'Data_clemens/' // 
     *      'MPI-EPS-LR_RCP45TempPrecOro2200Monthly.nc'
        cfile = cpath(m)(1:lenchr(cpath(m))) // cname
        call readdrive_clem (iuto,cfile,prdrive(1,1,1,m),'pr',nmon,
     *                       0, timeice)

        cname = 
     *      'Data_clemens/' // 
     *      'MPI-EPS-LR_RCP45TempPrecOro2200Monthly.nc'
        cfile = cpath(m)(1:lenchr(cpath(m))) // cname
        call readdrive_clem (iuto,cfile,tadrive(1,1,1,m),'tas',  nmon,
     *                       0, timeice)
      enddo
c|||||
#endif
c|||||

c        Compute matrix factors facice, facorb, facco2

c. . . . . . . . . . . .
#if defined (RACMOCLIM)
c. . . . . . . . . . . .
c     all done by calling getracmo_[ant,gre]_fields below
c. . .
#else
c. . .

      call fac_matrix (timeice, rco2, ecc, obl, prec,
     *                 facice, facorb, facco2, h)

c        Transfer gcm [to,ta,pr]drive to [to,ta,pr]ge with matrix method

      call domatrix (toge,todrive, facice,facorb,facco2, timeice,1,   1)
      call domatrix (tage,tadrive, facice,facorb,facco2, timeice,nmon,2)
      call domatrix (prge,prdrive, facice,facorb,facco2, timeice,nmon,3)
#  if defined (LOVECLIP)
      call domatrix (evge,evdrive, facice,facorb,facco2, timeice,nmon,5)
      call domatrix (qsge,qsdrive, facice,facorb,facco2, timeice,nmon,6)
#  endif
#  if defined (CO2INTER) 
      call domatrix (ruge,rudrive, facice,facorb,facco2, timeice,1,   4)
#  endif
c. . .
#endif
c. . .

      do jg=1,nlatg
        do ig=1,nlong
          do im=1,nmon
c.....................
#if defined (LOVECLIP)
c.....................
            tage(im,ig,jg) = tage(im,ig,jg) + tmelt ! deg C to K
c           prge(im,ig,jg) = prge(im,ig,jg)*1.      ! already m/y liq eq
c           evge(im,ig,jg) = evge(im,ig,jg)*1.      ! already m/y liq eq
c           qsge(im,ig,jg) = qsge(im,ig,jg)*1.      ! already W/m2
c....
#else
c....
            tage(im,ig,jg) = tage(im,ig,jg) + tmelt         ! deg C to K

#  if defined (RCMMATRIX) && defined (RCMTEMPSHIFT)
     *                     + RCMTEMPSHIFT             ! regcm3 cold bias
#  endif

#  if defined (GCMMATRIX) && defined (GCMTEMPSHIFT)
     *                     + GCMTEMPSHIFT     ! gcm (Genesis) temp shift
#  endif

            prge(im,ig,jg) = prge(im,ig,jg)*.001*365.   ! mm/day to m/yr

c.....
#endif
c.....
          enddo
        enddo
      enddo

c        Set monthly mean T,P on ice grid

c......................
#if  defined (RACMOCLIM)
c......................

#  if defined (GREENLAND)
      call getracmo_gre_fields (tog, tmon, pmon, runofa)
#  elif defined (EISANTA)
      call getracmo_ant_fields (tog, tmon, pmon, runofa, hsbal)
#  endif

c....
#else
c....

      call intg2i (toge, nlong, nlatg, tog, 
     *             indl, weil, indb, weib, nx, ny, 1)
      call intg2i (tage, nlong, nlatg, tmon, 
     *             indl, weil, indb, weib, nx, ny, nmon)
      call intg2i (prge, nlong, nlatg, pmon, 
     *             indl, weil, indb, weib, nx, ny, nmon)
#  if defined (LOVECLIP)
      call intg2i (evge, nlong, nlatg, emon, 
     *             indl, weil, indb, weib, nx, ny, nmon)
      call intg2i (qsge, nlong, nlatg, qmon, 
     *             indl, weil, indb, weib, nx, ny, nmon)
#  endif
#  if defined (CO2INTER)
      call intg2i (ruge, nlong, nlatg, runofa,             ! keep mm/day
     *             indl, weil, indb, weib, nx, ny, 1)
#  endif

c.....
#endif
c.....

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if ( defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM) ) && ( defined (GREENLAND) || defined (EISANTA) )  
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c       RCM, CSM or LOVECLIP anomaly method (observed - control) 
c       for air temps and precip. Only for Greenland or Antarctica.

c........................
      if (firstanom) then
c........................

c          Get observed climate (RACMO for Greenland, SeaRise for Ant.)

#  if defined (GREENLAND)
        call getracmo_gre_fields (togobs, tmonobs, pmonobs, zrudum)
#  elif defined (EISANTA)
        call climate_param (togobs, tmonobs, pmonobs, h, hs, hsbal,
     *                      0., 0., 0.,
     *                      1., 0.)
#  endif

c          Get model control (modern) climate, interp to ice grid

#  if defined (RCMANOM)
        cfile = cpathdrivercm(1:lenchr(cpathdrivercm)) //
c    *          'Runpre_40/LSX.nc'  ! preindustrial (CO2=280 ppmv)
     *          'Runcon40/LSX.nc'   ! present day   (CO2=355 ppmv)
        call readdrive_rcm (iuto,cfile,todrivectl,'TOPOGNOD', 1)
        call readdrive_rcm (iuto,cfile,prdrivectl,'PRECIP',nmon)
        call readdrive_rcm (iuto,cfile,tadrivectl,'TS2',   nmon)
#  elif defined (LOVEANOM)
        call readlovectl (todrivectl, 'surfelev',  1)
        call readlovectl (prdrivectl, 'precip', nmon)  
        call readlovectl (evdrivectl, 'evap',   nmon)  
        call readlovectl (tadrivectl, 'tair2m', nmon)  
#  elif defined (CSMANOM)
        call scopy (nlong*nlatg, todrive,1, todrivectl,1)

        cfile = cpathhome(1:lenchr(cpathhome)) // 
#    if defined (CSMHIRES)
     *    'Data_csm/CESM1/Ctl/TOPO_CTL_Clim.nc'                ! J. Tsai
#    else
     *    'Data_csm/Ctl/TOPO_CTL_Clim.nc'                      ! J. Tsai
#    endif
        call readdrive_csm (iuto,cfile,todrivectl,'PHIS', 1, 
     *                      1, 0.)

        cfile = cpathhome(1:lenchr(cpathhome)) // 
#    if defined (CSMHIRES)
     *    'Data_csm/CESM1/Ctl/PREC_CTL_Clim.nc'                ! J. Tsai
#    else
     *    'Data_csm/Ctl/PREC_CTL_Clim.nc'                      ! J. Tsai
#    endif
        call readdrive_csm (iuto,cfile,prdrivectl,'PREC_CTL', nmon, 
     *                      1, 0.)

        cfile = cpathhome(1:lenchr(cpathhome)) // 
#    if defined (CSMHIRES)
     *    'Data_csm/CESM1/Ctl/TEMP_CTL_Clim.nc'                ! J. Tsai
#    else
     *    'Data_csm/Ctl/TEMP_CTL_Clim.nc'                      ! J. Tsai
#    endif
        call readdrive_csm (iuto,cfile,tadrivectl,'SAT2m_CTL', nmon, 
     *                      1, 0.)
#  endif

        do j=1,nlatg
          do i=1,nlong
            togectl(i,j) = todrivectl(i,j)
            do im=1,nmon
              tagectl(im,i,j) = tadrivectl(i,j,im) 
     *                          + tmelt                     ! deg C to K
#  if defined (RCMTEMPSHIFT)
     *                          + RCMTEMPSHIFT        ! regcm3 cold bias
#  endif

#  if defined (RCMANOM) || defined (CSMANOM)
              prgectl(im,i,j) = prdrivectl(i,j,im)
     *                          * .001*365.             ! mm/day to m/yr
#  elif defined (LOVEANOM)
              prgectl(im,i,j) = prdrivectl(i,j,im)  ! LC is already m/yr
              evgectl(im,i,j) = evdrivectl(i,j,im)  ! LC is already m/yr
#  endif
            enddo
          enddo
        enddo

        call intg2i (togectl, nlong, nlatg, togctl, 
     *               indl, weil, indb, weib, nx, ny, 1)
        call intg2i (tagectl, nlong, nlatg, tmonctl, 
     *               indl, weil, indb, weib, nx, ny, nmon)
        call intg2i (prgectl, nlong, nlatg, pmonctl, 
     *               indl, weil, indb, weib, nx, ny, nmon)
#  if defined (LOVEANOM)
        call intg2i (evgectl, nlong, nlatg, emonctl, 
     *               indl, weil, indb, weib, nx, ny, nmon)
#  endif

        firstanom = .false.
c..........
      endif
c..........

c     diagnostic dump then stop:
c     do im=1,nmon
c       do j=1,ny
c         do i=1,nx
c           write (203,'(3i6,4f15.8)')
c    *        im,i,j,pmon(im,i,j), pmonctl(im,i,j), pmonobs(im,i,j),
c    *        pmon(im,i,j)-pmonctl(im,i,j)
c           if (im.eq.1) then
c             if (tog(i,j).ne.0. .or. togctl(i,j).ne.0. 
c    *            or. togobs(i,j).ne.0.)
c    *          write (210,'(2i6,3f15.4,  f15.4)')
c    *            i,j,tog(i,j), togctl(i,j),  togobs(i,j),
c    *            togobs(i,j)-togctl(i,j)
c             if (tog(i,j).ne.togctl(i,j))
c    *          write (210,'(a)') "TOG,  TOGCTL NOT EQUAL"
c           endif
c         enddo
c       enddo
c     enddo
c     stop 

c        Apply anomaly to air temp (additive) and precip (multiplic.)

      do j=1,ny
        do i=1,nx
          zdt = rlapse*(togobs(i,j)-togctl(i,j))
          do im=1,nmon
            tmon(im,i,j) = tmon(im,i,j) 
     *                   + tmonobs(im,i,j)-tmonctl(im,i,j)
     *                   + zdt
          enddo   
        enddo   
      enddo   

c---
c      new temp anomaly method, using annual means for obs and ctl
c      (avoiding crude cosine seas cycle added for obs in 
c      iceclimparam.F, retaining seasonal cycle of exp)             7/16
#  if ! ( defined (CSMANOM) && defined (GREENLAND) )
c     do j=1,ny
c       do i=1,nx
c         zdt = rlapse*(togobs(i,j)-togctl(i,j))
c         ztanno = 0. 
c         ztannc = 0. 
c         do im=1,nmon
c           ztanno =  ztanno + tmonobs(im,i,j)*ndaypm(im)/365. 
c           ztannc =  ztannc + tmonctl(im,i,j)*ndaypm(im)/365. 
c         enddo
c         do im=1,nmon
c           tmon(im,i,j) = tmon(im,i,j) 
c    *                   + ztanno - ztannc
c    *                   + zdt
c         enddo   
c       enddo   
c     enddo   
#  endif
c---

c!!!!!!!!!!

c!!!!!!!!!!
c     old precip anomaly method:
c     do j=1,ny
c       do i=1,nx
c         do im=1,nmon
c           zr = ( pmonobs(im,i,j) / max (pmonctl(im,i,j),1.e-6) )
c    *         * (2.**(zdt/tlapseprecip))
c           if (zr.ge.0.25 .and. zr.le.4.) then
c             pmon(im,i,j) = pmon(im,i,j) * zr
c           else
cc            if multiplic. factor is too big, try additive:
c             zp = pmon(im,i,j) + pmonobs(im,i,j) - pmonctl(im,i,j)
c             if (zp.gt.0.) then
c               pmon(im,i,j) = zp
c             else
cc              if additive yields neg. precip, simple setting:
c               pmon(im,i,j) = min (pmon(im,i,j), pmonobs(im,i,j))
c             endif
c           endif
c         enddo   
c       enddo   
c     enddo   
c!!!!!!!!!!
c     new precip anomaly method, using annual means of ctl and obs, 
c     to preserve seasonal cycle of model (3/15, for MT 800 kyr)
c     do j=1,ny
c       do i=1,nx
c         zavc = 0.
c         zavo = 0.
c         do im=1,nmon
c           zavc = zavc + pmonctl(im,i,j)
c           zavo = zavo + pmonobs(im,i,j)
c         enddo   
c         zavc = zavc/nmon
c         zavo = zavo/nmon
c         do im=1,nmon
c           zr = pmon(im,i,j) / max (zavc,1.e-6)
c           zr = max (.25, min (4., zr))
c           pmon(im,i,j) = zavo * zr
c         enddo   
c       enddo   
c     enddo   
c!!!!!!!!!!
c     new precip anomaly method, using annual means of ctl and obs,
c     to multiply by same ratio for all months and times at each pt
c     (3/28, for MT 800 kyr):
      do j=1,ny
        do i=1,nx
          zavc = 0.
          zavo = 0.
          do im=1,nmon
            zavc = zavc + pmonctl(im,i,j)
            zavo = zavo + pmonobs(im,i,j)
          enddo
          zavc = zavc/nmon
          zavo = zavo/nmon
          zr = zavo / max (zavc, 1.e-6)
          zr = max (.03, min (30., zr))
          do im=1,nmon
            pmon(im,i,j) = pmon(im,i,j) * zr
          enddo
        enddo
      enddo
c!!!!!!!!!!


c  for  evaporation (Loveclim only).  Not currently used -  
c  "emonobs" not available. Also exp's (non-anom) evap not currently
c  used - budgevap set to 0 in iceclimall, as for other models.     7/16
#  if defined (LOVEANOM)
c     do j=1,ny
c       do i=1,nx
c         do im=1,nmon
c           emon(im,i,j) = emon(im,i,j)
c    *                   + emonobs(im,i,j) - emonctl(im,i,j) ! not avail
c         enddo
c       enddo
c     enddo
#  endif


c&&&&&
#endif
c&&&&&

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#if defined (NHA) && ! defined (NHAGLAC) && ! defined (LOVECLIP)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c       For NHA, apply corrections for ice albedo, co2, 
c       and orographic precip effect dep. on ice-gcm surface slope.
c       Nb: orb variations done via gcm matrix, unless NHAOBS.
c       Nb: lapse rate from gcm to ice sfc elev is applied in surfbudg.
c
c     call calcalbeff (h, albeff, timeice)
c
c     lapse rate from gcm to ice sfc elev applied in surfbudg.
c     do j=1,ny
c       jm1=max (j-1,1)
c       jp1=min (j+1,ny)
c       do i=1,nx
c         im1=max (i-1,1)
c         ip1=min (i+1,nx)
c
c         zdt = 0.
c         zptop = 1.
c
c         zdt = - 10.*albeff(i,j)                          ! ice alb
c    *          - 15.*albeff(i,j)                          ! ice alb
c    *          - 20.*albeff(i,j)                          ! ice alb
c    *          + 4.*(log(rco2)/log(2.))                   ! co2
c
c         zix = (hs(ip1,j)-hs(im1,j)) / (2.*dx(i,j))       ! orog precip
c         ziy = (hs(i,jp1)-hs(i,jm1)) / (2.*dy(i,j))
c         zis = sqrt(zix**2+ziy**2)
c         zgx = (tog(ip1,j)-tog(im1,j)) / (2.*dx(i,j))
c         zgy = (tog(i,jp1)-tog(i,jm1)) / (2.*dy(i,j))
c         zgs = sqrt(zgx**2+zgy**2)
c         zptop = max (0.3, min (3., 2.**((zis-zgs)/.01) ))
c
c         do im=1,nmon
c           tmon(im,i,j) = tmon(im,i,j) + zdt
c           pmon(im,i,j) = pmon(im,i,j) 
c    *                   * (2.**(zdt/tlapseprecip)) * zptop
c         enddo
c       enddo
c     enddo

c>>>>>
#endif
c>>>>>

c^^^^^^^^^^^^^^^^^^^^^
#if defined (LOVECLIP)
c^^^^^^^^^^^^^^^^^^^^^
c        Do domain-wide correction to mean precip and evap on ice sfcs, 
c        to (approx) conserve water between Loveclim and ice model.
c        Not done for anomaly method (-DLOVEANOM), which violates
c        conservation O(1)...possible spurious effects with globcorrect.
#  if ! defined (LOVEANOM)
      call globcorrect (prge, evge, pmon, emon, h)
#  endif
c^^^^^
#endif
c^^^^^

c     check topog maps:
c     write (188,*) 'toge:'
c     isk = 10
c     do j=nlatg,1,-isk
c       write (188,'(500i4)') (nint(toge(i,j)/10.),i=1,nlong,isk)
c     enddo
c 
c     write (188,*)
c     write (188,*) 'tog:'
c     do j=ny,1,-isk
c       write (188,'(500i4)') (nint(tog(i,j)/10.),i=1,nx,isk)
c     enddo
c     call flush (188)

      return
      end

c-----------------------------------------------------------------------

      subroutine readdrive (iu, cfile, fielddrive, ckey, nmonin)

c     Reads mid-month driving fields from a Genesis (debplotdriveice)
c     data input file (iu,cfile) into array fielddrive, for months
c     1 to nmonin. nmonin = nmon for all calls, except 1 for topog 
c     fields (ckey='TO'). Reads entire file into memory.

#include <comicephys.h>
#include <comicegrid.h>

      character*(*) cfile, ckey
      dimension fielddrive (nlong,nlatg,nmonin)

c        Open existing data input file

      open (iu, file=cfile, form='formatted', status='old')

c        Read scaling factor

      read (iu,'(7x,f12.5)') scale

c        Read monthly-mean global fields

      do 210 imon=1,nmonin
        call readdata_ice (iu, fielddrive(1,1,imon), nlong,nlatg,
     *                     nlong,nlatg, ckey,2, imon)
  210 continue

      close (iu)

c          Divide by scaling factor

      do imon=1,nmonin
        do jg=1,nlatg
          do ig=1,nlong
            fielddrive(ig,jg,imon) = fielddrive(ig,jg,imon) / scale
          enddo
        enddo
      enddo

      call checkfield (fielddrive, nmonin, ckey, -10000., 100000.)

      return
      end

c----------------------------------------------------------------------

      subroutine checkfield (fielddrive,nmonin,ckey,checkmin,checkmax)

c Checks extreme values in array fielddrive 

#include <comicephys.h>
#include <comicegrid.h>

      character*(*) ckey

      dimension fielddrive (nlong,nlatg,nmonin)

      zmin =  1.e20
      zmax = -1.e20
      zav  = 0.
      zava = 0.
      do imon=1,nmonin
        do jg=1,nlatg
#if defined (RCMMATRIX) || ( defined (RACMOCLIM) && defined (EISANTA) )
          zwei = 1.
#else
          zwei = cos(alatg(jg))
#endif
          do ig=1,nlong
            zmin = min (zmin, fielddrive(ig,jg,imon))
            zmax = max (zmax, fielddrive(ig,jg,imon))
            zav  = zav  + fielddrive(ig,jg,imon)*zwei
            zava = zava +                        zwei
          enddo
        enddo
      enddo
      zav = zav/zava

      if (zmin.lt.checkmin .or. zmax.gt.checkmax) then
        write (ioterm,100) ckey, zav, zmin, zmax
      else
        write (ioterm,102) ckey, zav, zmin, zmax
      endif
  100 format('**OUT**  ', a8, '  av,min,max =',3e15.4)
  102 format('   (in)  ', a8, '  av,min,max =',3f10.3)

      return
      end

c-----------------------------------------------------------------------

      subroutine fac_matrix (time, rco2, ecc, obl, prec, 
     *                       facice, facorb, facco2, h)

c        Computes weights facice, facorb, facco2 for matrix method.
c
c        time = current time, -yr BP (eg, -21000) (supplied)
c        facice = ice extent interpolation matrix factor (0-1,returned)
c        facorb = orbit interpolation matrix factor (W/m2, returned)
c        facco2 = co2 interpolation matrix factor (returned)
c        h     = ice sheet thickness, m (supplied)
c        darea = ice grid box area, m**2 (supplied)

#include <comicephys.h>
#include <comicegrid.h>

      dimension h(nx,ny)

c#######################################################################
#if defined (LOVECLIP) || defined (CSMCLIM) || defined (CLEMCLIM) || defined (CSMORDO) || ( defined (NHA) && defined (NHAOBS) )
c#######################################################################

c         Could probably just test on nmat.eq.1 (also in domatrix)

      facice = 0.
      facorb = 0.
      facco2 = 0.

c########################
#elif defined (RCMMATRIX)
c########################

      facice = 0.
      facorb = 0.
      facco2 = log(rco2)

c##############################################
#elif defined (GCMMATRIX) || defined (CO2INTER)
c##############################################

c        Compute facorb, current summer-half-year insolation at alorb
      
      call insolhalf (nint(time), alorb*pi/180., zflux, ecc,obl, prec)
      facorb = zflux
c     write (196,'(/i12,4(a,f10.5))') nint(time), 'ecc=',ecc,
c    *  '  obl=',obl*180./pi, '  prec=',prec*180./pi,'  zflux=',zflux
c     call flush (196)

c        Compute facice, current ice-sheet surface area 
c        (albedo/topog feedbck)

      ztotari = 0.
      ztotice = 0.
      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.001) ztotari = ztotari + darea(i,j)
          ztotice = ztotice + h(i,j)*darea(i,j)
        enddo
      enddo
      ztotari = ztotari*1.e-6 ! m2 to km2
      ztotice = ztotice*1.e-9 ! m3 to km3
# if defined (EISANTA)
      facice = min (1., ztotari/13.e6)         ! modern Antarctica
# elif defined (NHA)
      facice = min (1., ztotari/17.e6)         ! LGM N.Amer. + Greenland
c     facice = min (1., ztotice/33.e6)         ! LGM N.Amer. + Greenland
# endif

c        Compute facco2 (radiative effect of co2 is logarithmic
c        wrt concentration):

      facco2 = log(rco2)

c     write (6,'(a, i10, 3f8.3, i10, f8.3)') 
c    *  'time dflux facorb facco2 totice facice: ',
c    *  nint(time), zflux, facorb, facco2, nint(ztotice), facice

c#####
#endif
c#####

      return
      end

c-----------------------------------------------------------------------

      subroutine domatrix (arr, arrm, facice, facorb, facco2, timeice,
     *                     nmonin, itype)

c       Does "matrix method" interpolation of a meteorological driving
c       field. 
c       Computes:
c         sigma[weimat(m)*arrm(..,..,m)]   for topog,temp.
c         sigma[weimat(m)*arrm(..,..,m)]   for precip.
c         where weimat(m) are matrix weights depending on facice,facorb,
c               facco2,and the corresp values of the nmat experiments.
c       All arguments supplied unless noted.
c         arr  = interpolated field (returned, dim: nmonin,lon,lat)
c         arrm = matrix of experiment fields (dim:lon,lat,nmonin,nmat)
c         facice = weighting factor for ice (computed in fac_matrix)
c                  (0 to 1)
c         facorb = weighting factor for orbit (" " " " " " " " " " )
c                  (summer half-yr insol at alorb (fac_matrix), W/m2)
c         facco2 = weighting factor for co2   (" " " " " " " " " " )
c         itype  = 1 for topog  (additive), 
c                  2 for temp   (additive),
c                  3 for precip (additive (not ratio?))
c                  4 for runoff (additive (not ratio?))
c
c         nb: note switch of order of nmonin in arr vs arrm dimensions

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  arr(nmonin,nlong,nlatg), arrm(nlong,nlatg,nmonin,nmat)

      dimension exporb(3), expice(3), expco2(4)
      save exporb, expice, expco2

c     currently set by insolhalf below:
c     data exporb             ! facorb values, half-yr summer
c    *  /339.0, 360.1, 380.6/ ! insol (70S) as fluxco/me/wa below
c     from solin, plot.in.ant on pangea:
c    *  /331.0, 355.7, 378.3/ ! insol (70S) in old co,me,wa GCM expts
c    *  /380.0, 395.0, 420.0/ ! insol (55N) in 0k,11k,21k GCM expts.

#if defined (NHA)
c     data expice           ! facice values, representing Laurentide
c    *  /0.,   .25,    1./  ! extents in NH 0k,11k,21k GCM expts !777777
      data expice           ! facice values, representing Laur+Greenland
     *  /.13,  .30,    1./  ! extents in NH 0k,11k,21k GCM expts !777777
#else
      data expice           ! Antarc area/max_area in co,me,wa GCM expts
     *  /0.,  .528,    1./  ! from mask, plot.in.ant on pangea.
#endif

      dimension weiice(3), weiorb(3), weico2(4), weimat(nmat)

      logical first
      data first /.true./
      save first

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#if defined (LOVECLIP) || defined (CSMCLIM) || defined (CLEMCLIM) || defined (CSMORDO) || ( defined (NHA) && defined (NHAOBS) )
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

c         For all these, nmat=1 (comicegrid.h), so just transfer field 
c         and return. Could probably just test on nmat.eq.1.

      do imon=1,nmonin
        do j=1,nlatg
          do i=1,nlong
            arr(imon,i,j) = arrm(i,j,imon,1)
            if (itype.eq.3 .or. itype.eq.6)         ! precip, dwnd solar
     *         arr(imon,i,j) = max (arr(imon,i,j), 0.)
          enddo
        enddo
      enddo
      return
c@@@@@
#endif
c@@@@@

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if defined (GCMMATRIX) || defined (CO2INTER)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c         Compute reference summer-half-year insolation values

      if (first) then

# if defined (EISANTA)
        call insolhalf (0, alorb*pi/180., fluxco,
     *                  .05, 22.5*pi/180., 270.*pi/180.)
        call insolhalf (0, alorb*pi/180., fluxme,
     *                  .0 , 23.5*pi/180.,   0.*pi/180.)
        call insolhalf (0, alorb*pi/180., fluxwa,
     *                  .05, 24.5*pi/180.,  90.*pi/180.)
# elif defined (NHA)
        call insolhalf (0, alorb*pi/180., fluxco,
     *                  .05, 22.5*pi/180.,  90.*pi/180.)
        call insolhalf (0, alorb*pi/180., fluxme,
     *                  .0 , 23.5*pi/180.,   0.*pi/180.)
        call insolhalf (0, alorb*pi/180., fluxwa,
     *                  .05, 24.5*pi/180., 270.*pi/180.)
# endif
        write (6,*) 'fluxco=',fluxco
        write (6,*) 'fluxme=',fluxme
        write (6,*) 'fluxwa=',fluxwa
        exporb(1) = fluxco
        exporb(2) = fluxme
        exporb(3) = fluxwa

# if defined (EISANTA)
c       for nmat=27:
cc      expco2(1) = log(1.)
cc      expco2(2) = log(2.)
cc      expco2(3) = log(4.)
c       expco2(1) = log(1.)
c       expco2(2) = log(4.)
c       expco2(3) = log(8.)

c       for nmat=36:
c       expco2(1) = log(0.5)   ! c05
c       expco2(2) = log(1.)    ! c05
c       expco2(3) = log(2.)    ! c05
c       expco2(4) = log(4.)    ! c05
        expco2(1) = log(1.)    ! c05
        expco2(2) = log(2.)    ! c05
        expco2(3) = log(4.)    ! c05
        expco2(4) = log(8.)    ! c05
# elif defined (NHA)
        expco2(1) = log(185./280.)  ! lgm        (not used)
        expco2(2) = log(1.)         ! preindust  (not used)
        expco2(3) = log(2.)         !            (not used)
# endif

        first = .false.
      endif

c~~~~~
#endif
c~~~~~

c######################
#if defined (RCMMATRIX)
c######################

c       Just use CO2 weighting between RCM climates.
c       nb: nmat is set in comicegrid.h. cpathdrivercm set in
c       comicepath.h. rcm filenames set above in subr climate_gcm.
c       rco2 vs time of run set in icectl.F, subr longterm.
c       facco2 = log(rco2) has been set in subr fac_matrix.

      expco2(1) = log(1.)

c........................
#  if defined (RCMMATCO2)
c........................

#    if defined (RAMPPLIO) || defined (RAMP2X) || defined (RAMP4X) || defined (RAMP8X) || defined (RAMP8XW)
#      if defined (RAMPPLIO)
        expco2(2) = log(1.43)                                     ! plio
#      elif defined (RAMP2X)
        expco2(2) = log(2.)                                       ! 2x
#      elif defined (RAMP4X)
        expco2(2) = log(4.)                                       ! 4x
#      elif defined (RAMP8X) || defined (RAMP8XW)
        expco2(2) = log(8.)                                       ! 8x
#      endif
#    else
        expco2(2) = log(2.)
        expco2(3) = log(4.)
        expco2(4) = log(8.)
#    endif

c......
#  else
c......

c       if not RCMMATCO2, nmat=2, and both rcm files are RCMFILE or 
c       pre_40 (set in climate_gcm). Just set expco2(2) arbitrarily 
c       to allow "interpolation" wrt rco2 below.
        expco2(2) = log(2.)      

c.......
#  endif
c.......

c        Interpolate with respect to log(rco2). Nb: extrapolate beyond 
c        envelope of climates if rco2 < 1 or > 4 (or 8).

      call zero (weico2, nmat)
      do k=2,nmat
        if (       (facco2.lt.expco2(k  ).or.k.eq.nmat)
     *       .and. (facco2.ge.expco2(k-1).or.k.eq.2) ) then
          weico2(k) = (facco2-expco2(k-1))/(expco2(k)-expco2(k-1))
          weico2(k-1) = 1.-weico2(k)
        endif     
      enddo 

      do m=1,nmat
        weimat(m) = weico2(m)
      enddo

c######################
#elif defined (EISANTA)
c######################
 
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c     Prescribed matrix conditions:
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c---
c     no facice var (no albedo fbk):
c     facice = expice(1)              ! no  ice
c     facice = expice(2)              ! med ice
c---
c     no facco2 var (CO2 2x,3x):
c     facco2 = log(2.)
c---
c     no facorb variation: 
c     facorb = exporb(1)              ! cold orb
c     facorb = exporb(2)              ! med  orb
c     facorb = exporb(3)              ! warm orb
c---
c     facice = 0.                     ! mis31 gcm
c     facorb = dtantjan/0.1 + 529.85  ! mis31 gcm (Jan insol, 80S, W/m2)
c                                     ! nb: dtantjan not passed! 
c     facco2 = 1.                     ! mis31 gcm (CO2 2x,3x)
c---
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c     constrain weiice's to be <= 1:
      weiice(1) = min( 1., max (0.,
     *                 (expice(2)-facice)/(expice(2)-expice(1)) ))
      weiice(3) = min (1., max (0.,
     *                  (facice-expice(2))/(expice(3)-expice(2)) ))
      weiice(2) = 1. - (weiice(1)+weiice(3))         !weiice's sum to 1

c     allow weiorb's to be > 1:
      weiorb(1) = max (0., (exporb(2)-facorb)/(exporb(2)-exporb(1)))
      weiorb(3) = max (0., (facorb-exporb(2))/(exporb(3)-exporb(2)))
      weiorb(2) = 1. - (weiorb(1)+weiorb(3))        ! weiorb's sum to 1

c     no orbital signal, use "me" orbit
#  if defined (ORBMED)
      weiorb(1) = 0.
      weiorb(2) = 1.
      weiorb(3) = 0.
#  endif

      call zero (weico2, 4)
      do k=2,4
        if (       (facco2.lt.expco2(k  ).or.k.eq.4) 
     *       .and. (facco2.ge.expco2(k-1).or.k.eq.2) ) then
          weico2(k) = (facco2-expco2(k-1))/(expco2(k)-expco2(k-1))
          weico2(k-1) = 1.-weico2(k)
        endif
      enddo

      if (nmat.eq.1) then
c       if nmat=1, weights set above are not used, just for diag output 
        weimat(1) = 1.
      else
c       i,j,k loop limits must correspond to nmat (=i*j*k) and csuffant
c       settings (orb,ice,co2 dimensions):
        zw = 0.
        nzw = 0
        do k=1,4
          do j=1,3
            do i=1,3
              m = i + (j-1)*3 + (k-1)*9
              weimat(m) = weiorb(i)*weiice(j)*weico2(k)
              zw = zw + weimat(m)
              nzw = nzw + 1
            enddo
          enddo
        enddo
        if (nzw.ne.nmat) then
          write (ioterm,'(a,2i6)') 
     *      '*** Error (domatrix EISANTA): nzw, nmat =', nzw, nmat
          stop
        endif
        if (abs(zw-1.).gt.1.e-6) then
          write (ioterm,'(a,e15.6)') 
     *      '*** Error (domatrix EISANTA): zw = ', zw
          stop
        endif

c       Or (mis31 gcm):
c       z1070ma = 450.37   ! Jan insol, 80S, 1.070 Ma (W/m2)
c       z1082ma = 585.28   ! Jan insol, 80S, 1.082 Ma (W/m2)
c       weimat(1) = (z1082ma - facorb) / (z1082ma - z1070ma) 
c       weimat(1) = max (0., min (1., weimat(1)))
c       weimat(2) = 1. - weimat(1) 
      endif


c########################################
#elif defined (NHA) && ! defined (NHAOBS)
c########################################
c NHA + NHAOBS was set above, and returned. 

c......................
#  if defined (NHAGLAC)
c......................

c     3 weiice (No,11k,Fu), 2 weiorb (co,wa), no weico2.

      weiice(1) = max ( 0., (expice(2)-facice)/(expice(2)-expice(1)) )
      weiice(3) = max ( 0., (facice-expice(2))/(expice(3)-expice(2)) )
#    if defined (ICEAMP)
      weiice(1) = ICEAMP * weiice(1)
      weiice(3) = ICEAMP * weiice(3)
#    else
      weiice(1) = min ( 1., weiice(1) )
      weiice(3) = min ( 1., weiice(3) )
#    endif
      weiice(2) = 1. - (weiice(1)+weiice(3))         !weiice's sum to 1

c     No "me" orbit forcing available, only "co" and "wa":

      weiorb(1) = (exporb(3)-facorb)/(exporb(3)-exporb(1))
#    if defined (ORBAMP)
      weiorb(1) = ORBAMP * (weiorb(1)-0.5) + 0.5
#    else
      weiorb(1) = min ( 1., max ( 0., weiorb(1) ))
#    endif
      weiorb(2) = 1. - weiorb(1)
      weiorb(3) = 0.  ! for niceity only

      zw = 0.
      nzw = 0
      do j=1,3
        do i=1,2
          m = i + (j-1)*2
          weimat(m) = weiorb(i)*weiice(j)
          if (i.eq.1 .and. j.eq.1) weimat(m) = weimat(m) - 0.5 ! new
          if (i.eq.2 .and. j.eq.1) weimat(m) = weimat(m) - 0.5 ! new
          zw = zw + weimat(m)
          nzw = nzw + 1
        enddo
      enddo
c     weimat(7) = 0.                                    ! old
      weimat(7) = 1.                                    ! new
      zw = zw + weimat(7)
      nzw = nzw + 1

      if (nzw.ne.nmat) then
        write (ioterm,'(a,2i6)') 
     *    '*** Error (domatrix NHAGLAC): nzw, nmat =', nzw, nmat
        stop
      endif
      if (abs(zw-1.).gt.1.e-6) then
        write (ioterm,'(a,e15.6)') 
     *    '*** Error (domatrix NHAGLAC): zw = ', zw
        stop
      endif

      write (222,
     * "(f10.2, '   orb:',3f8.3,'   ice:',3f8.3,'   mat:'5x,7f8.3)")
     *  timeice, weiorb, weiice, weimat

c......
#  else
c......

c     for now, NHA nmat=3 (comicegrid.h), co/me/wa orbit (else modern):
      facorb = 0.5*(exporb(1)+exporb(3)) 
     *       + 0.5*(exporb(3)-exporb(1)) * sin(2.*pi*timeice/40.e3)
      facice = 0.
      facco2 = 1.

c     allow weiorb's to be > 1:
      weiorb(1) = max (0., (exporb(2)-facorb)/(exporb(2)-exporb(1)))
      weiorb(3) = max (0., (facorb-exporb(2))/(exporb(3)-exporb(2)))
      weiorb(2) = 1. - (weiorb(1)+weiorb(3))        ! weiorb's sum to 1

      do m=1,nmat
        weimat(m) = weiorb(m)
      enddo

c.......
#  endif
c.......

c#####
#endif
c#####

c        Do matrix weighting of climate variables (topo,T,P, or runoff)

      do imon=1,nmonin
        do j=1,nlatg
          do i=1,nlong
            arr(imon,i,j) = 0.
            do m=1,nmat
              arr(imon,i,j) = arr(imon,i,j) + weimat(m)*arrm(i,j,imon,m)
            enddo
c           precip:
            if (itype.eq.3 .or. itype.eq.4)
     *         arr(imon,i,j) = max (arr(imon,i,j), 0.)
          enddo
        enddo
      enddo

c        Reset facorb,facco2 for output in fort.22 - must be last call 
c        this timestep

c~~~~~~~~~~~~~~~~~~~~~~
#if defined (RCMMATRIX)
c~~~~~~~~~~~~~~~~~~~~~~
      if (itype.eq.3) then
        facorb = -1.*weiorb(1) + 0.*weiorb(2) + 1.*weiorb(3)

        zfac = facco2
        do k=2,nmat
          if (       (zfac.lt.expco2(k  ).or.k.eq.nmat) 
     *         .and. (zfac.ge.expco2(k-1).or.k.eq.2) ) then
            facco2 = weico2(k) + (k-1)
          endif
        enddo
      endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#elif defined (EISANTA) || ( defined (NHA) && ! defined (NHAOBS) )
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  if defined (CO2INTER)
      if (itype.eq.4) then
#  else
      if (itype.eq.3) then
#  endif

#  if defined (NHAGLAC)
        facice = -1.*weiice(1) + 0.*weiice(2) + 1.*weiice(3)
        facorb = weiorb(1)
        facco2 = 0.
#  else
        facice = -1.*weiice(1) + 0.*weiice(2) + 1.*weiice(3)
        facorb = -1.*weiorb(1) + 0.*weiorb(2) + 1.*weiorb(3)
        zfac = facco2
        do k=2,4
          if (       (zfac.lt.expco2(k  ).or.k.eq.4) 
     *         .and. (zfac.ge.expco2(k-1).or.k.eq.2) ) then
            facco2 = weico2(k) + (k-1)
          endif
        enddo
#  endif

      endif
c~~~~~
#endif
c~~~~~

      return
      end

c-----------------------------------------------------------------------

      subroutine intg2i (arrg, nlong, nlatg, arri, 
     *                   indl, weil, indb, weib, nx, ny, nmon)

c     Interpolates arrg (global gcm/rcm grid) to arri (ice model grid),
c     for nmon time periods, using indices indl, weil, etc.

      dimension 
     *  arrg(nmon,nlong,nlatg), arri(nmon,nx,ny), 
     *  indl(nx,ny),            weil(nx,ny),
     *  indb(nx,ny),            weib(nx,ny)

      do j=1,ny
        do i=1,nx

c         indices and weights into agcm-grid arrays (setinterp):
          il = indl(i,j)
#if defined (RCMMATRIX)
          ir = min (il+1, nlong)
#else
          ir = mod (il,nlong) + 1
#endif
          jb = indb(i,j)
          jt = min (jb+1, nlatg)
c         only needed if switched N<->S (LOVECLIP, done in setinterp):
          jb = max (1, min (nlatg, jb))
          wlb =     weil(i,j) *    weib(i,j)
          wrb = (1.-weil(i,j))*    weib(i,j)
          wlt =     weil(i,j) *(1.-weib(i,j))
          wrt = (1.-weil(i,j))*(1.-weib(i,j))
          do im=1,nmon
            arri(im,i,j) = arrg(im,il,jb)*wlb + arrg(im,ir,jb)*wrb
     *                   + arrg(im,il,jt)*wlt + arrg(im,ir,jt)*wrt
          enddo

        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

c++++++++++++++++++++++++
#  if defined (RCMMATRIX)
c++++++++++++++++++++++++

      subroutine readdrive_rcm (iu, cfile, fielddrive, cvarin, nmonin)

c     Reads mid-month driving fields from RCM (RegCM3) netcdf history 
c     file output (iu,cfile) into array fielddrive, for months
c     1 to nmonin. nmonin = nmon for all calls, except 1 for topog
c     fields (ckey='TO'). 

#include <comicephys.h>
#include <comicegrid.h>
#include <netcdf.inc>

      character*(*) cfile, cvarin
      dimension fielddrive (nlong,nlatg,nmonin)

      integer*4 ncid, rcode, varid, start(3), count(3)
      character*16 cvar

      real*4 arrin(nlong,nlatg,12)

      cvar = cvarin

c        Open existing RegCM3 netcdf file and read requested field

      rcode = nf_open(cfile, nf_nowrite, ncid)
      call checkrcode (rcode, cfile, 'open', 6000)

c        Define starting/ending position eof 2-D data values in array

      start(1) = 1
      count(1) = nlong

      start(2) = 1
      count(2) = nlatg

      start(3) = 1        
      count(3) = nmonin

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfile, cvar, 7760)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfile, cvar, 7761)

      do m=1,nmonin
        do j=1,nlatg
          do i=1,nlong
            if (cvar.eq.'PRECIP') then
              fielddrive(i,j,m) = arrin(i,j,m)*1000.*86400.! m/s to mm/d
            else if (cvar.eq.'TS2' .or. cvar.eq.'TS10') then
              fielddrive(i,j,m) = arrin(i,j,m) - tmelt     ! deg K to C
            else 
              fielddrive(i,j,m) = arrin(i,j,m)
            endif
          enddo
        enddo
      enddo

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 7762)

      return
      end

c-----------------------------------------------------------------------

      subroutine setinterp_rcm (nxrcm, nyrcm, dd0rcm,
     *                          xh, yh, nx, ny,
     *                          indl, weil, indb, weib)

c        Calculate indices and weights to interpolate from RCM (RegCM3)
c        grid (nxg,nyg) to the ice grid (nx,ny), both polar stereo.
c        Returns [indl,weil] from the ice grid into the RegCM3 "e-w" 
c        axis, and [indb,weib] from the ice grid into RegCM3 "n-s" axis.

      dimension
     *  xh(nx,ny),     yh(nx,ny),
     *  indl(nx,ny),   weil(nx,ny),
     *  indb(nx,ny),   weib(nx,ny)

      do j=1,ny
        do i=1,nx
          zx = (xh(i,j) + dd0rcm*(0.5*nxrcm)) / dd0rcm
          indl(i,j) = max (1, min (nxrcm, nint (zx) ))
          weil(i,j) = max (0., min (1., 0.5-(zx-indl(i,j)) ))

          zy = (yh(i,j) + dd0rcm*(0.5*nyrcm)) / dd0rcm
          indb(i,j) = max (1, min (nyrcm, nint (zy) ))
          weib(i,j) = max (0., min (1., 0.5-(zy-indb(i,j)) ))
        enddo
      enddo

      return
      end 

c+++++++
#  endif
c+++++++

c-----------------------------------------------------------------------

c||||||||||||||||||||||
#  if defined (CSMCLIM)
c||||||||||||||||||||||

      subroutine readdrive_csm (iu, cfilin, fielddrive, cvarin, nmonin,
     *                          ifanom, timeice)

c     Reads mid-month driving fields from CSM climate file into array 
c     fielddrive, for months 1 to nmonin. nmonin = nmon for all calls, 
c     except 1 for topog (cvarin ='PHIS'). 
c     ifanom = 0 for experiment, = 1 for anomaly ctrl (diagnostic only).

#include <comicephys.h>
#include <comicegrid.h>
#include <netcdf.inc>

      character*(*) cfilin, cvarin
      dimension fielddrive (nlong,nlatg,nmonin)

      integer*4 ncid, rcode, varid, start(4), count(4)
      character cfile*240, cvar*16

      real*4 arrin(nlong,nlatg,12), brrin(nlong,nlatg,12)
      dimension work(nlong), work2(nlong,nlatg)

      dimension timcsm(ntcsm)              ! ntcsm in comicegrid.h
      save timcsm                          ! timcsm not used for J. Tsai

      logical firstcsm, firstyear
      data firstcsm, firstyear /.true., .true./
      save firstcsm, firstyear

      cfile = cfilin
      cvar = cvarin

c J. Tsai: (all years, nyearstart must = 0, use just one current yr):
      ita = max (1, min (ntcsm, nint(timeice + 0.5) ))

c     write (ioterm,'(3a,2i6)')
c    *  'reading existing Netcdf ', cfile(1:lenchr(cfile)),
c    *  ' data file', ifanom, ita

c        Open existing CSM netcdf file and read requested field

      rcode = nf_open(cfile, nf_nowrite, ncid)
      call checkrcode (rcode, cfile, 'open', 6000)

c        If first call, read latitude and longitude grids into 
c        alatg, along, in common (comicegrid.h). Convert to radians.
c        nb: lon and lat are doubles on nc file.

c-----------------------
      if (firstcsm) then
c-----------------------

c          latitudes:

        start(1) = 1
        count(1) = nlatg
        start(2) = 1        
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lat', varid)
        call checkrcode (rcode, cfile, cvar, 7720)
        rcode = nf_get_vara_double (ncid, varid, start, count, alatg)
        call checkrcode (rcode, cfile, cvar, 7721)
        do j=1,nlatg
          alatg(j) = alatg(j) * pi/180.
        enddo

c       protect interpolation logic in setinterp (needed?):
        alatg(1)     = max (alatg(1),     -89.999 *pi/180.) 
        alatg(nlatg) = min (alatg(nlatg),  89.999 *pi/180.) 

c          longitudes:

        start(1) = 1
        count(1) = nlong
        start(2) = 1        
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lon', varid)
        call checkrcode (rcode, cfile, cvar, 7730)
        rcode = nf_get_vara_double (ncid, varid, start, count, work)
        call checkrcode (rcode, cfile, cvar, 7731)
c       shift cut from greenwich to dateline for ice model:
        do i=1,nlong
          ii = mod (i + nlong/2 -1, nlong)  + 1
          along(i) = work(ii) * pi/180.
          if (along(i).gt.pi) along(i) = along(i) - 2.*pi
        enddo

c       protect interpolation logic in setinterp (needed?):
        if ( along(1).gt. 179.999*pi/180. .or. 
     *       along(1).lt.-179.999*pi/180. ) 
     *    along(1) = -179.999 * pi/180.
        if ( along(nlong).gt. 179.999*pi/180. .or. 
     *       along(nlong).lt.-179.999*pi/180. ) 
     *    along(nlong) = 179.999 * pi/180.

        firstcsm = .false.

c----------
      endif
c----------

c.............................
      if (cvar.eq.'PHIS') then
c.............................

c        topography:

      start(1) = 1
      count(1) = nlong
      start(2) = 1
      count(2) = nlatg
      start(3) = 1        
      count(3) = 1
      rcode = nf_inq_varid (ncid, 'Topography', varid)
      call checkrcode (rcode, cfile, 'Topography', 7735)
      rcode = nf_get_vara_double (ncid, varid, start, count, work2)
      call checkrcode (rcode, cfile, 'Topography', 7736)
      do j=1,nlatg
        do i=1,nlong
c         shift cut from greenwich to dateline for ice model:
          ii = mod (i + nlong/2 -1, nlong)  + 1
c         fielddrive(i,j,1) = work2(ii,j) / 9.80616     ! geopot to elev
          fielddrive(i,j,1) = work2(ii,j)               ! already elev
        enddo
      enddo

c.........
      else
c.........

c        times (years) (ntcsm is in comicegrid.h):
   
      if (firstyear) then
        if (ntcsm.eq.1) then
          timcsm(1) = 0.
        else
          start(1) = 1
          count(1) = ntcsm
          start(2) = 1        
          count(2) = 1
          rcode = nf_inq_varid (ncid,  'year',  varid)  
          call checkrcode (rcode, cfilin, 'year', 7740)
          rcode = nf_get_vara_double (ncid, varid, start, count, timcsm)
          call checkrcode (rcode, cfilin, 'year', 7741)
        endif
        firstyear = .false.
      endif

c---
c not J. Tsai:
c
c        Find surrounding year indices and weight (ita, itb, wta)
c
c     if (timeice.le.timcsm(1)) then
c       ita = 1
c       itb = 1
c       wta = 1.
c     else if (timeice.ge.timcsm(ntsm)) then
c       ita = ntsm
c       itb = ntsm
c       wta = 1.
c     else
c       do m=2,ntcsm
c         if (timeice.le.timcsm(m)) then
c           ita = m-1
c           itb = m
c           wta = (timcsm(itb) - timeice) / (timcsm(itb) - timecsm(ita))
c           goto 10
c         endif
c       enddo
c  10   continue
c     endif
c---

c        Read field 

      start(1) = 1
      count(1) = nlong
      start(2) = 1
      count(2) = nlatg
      start(3) = 1        
      count(3) = nmonin
      if (ifanom.eq.0) then
        start(4) = ita
      else if (ifanom.eq.1) then
        start(4) = 1
      endif
      count(4) = 1

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfile, cvar, 7760)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfile, cvar, 7761)

c not J. Tsai:
c     start(4) = itb        
c     rcode = nf_get_vara_real (ncid, varid, start, count, brrin)
c     call checkrcode (rcode, cfile, cvar, 7762)

      do m=1,nmonin
        do j=1,nlatg
          do i=1,nlong

c           shift cut from greenwich to dateline for ice model:
            ii = mod (i + nlong/2 -1, nlong)  + 1

c not J. Tsai:
c           zrrin = wta*arrin(ii,j,m) + (1.-wta)*brrin(ii,j,m)
c J. Tsai:
            zrrin = arrin(ii,j,m)

            if (cvar.eq.'T' .or. cvar.eq.'SAT2m_CTL') then
              fielddrive(i,j,m) = zrrin - tmelt             ! deg K to C
            else if (cvar.eq.'PRECT' .or. cvar.eq.'PREC_CTL') then
              fielddrive(i,j,m) = zrrin*1000.*86400.       ! m/s to mm/d
            else 
              write (ioterm,'(/2a)')
     *          '*** Error (readdrive_csm): unknown fieldname: ', cvar
              stop
            endif

          enddo
        enddo
      enddo

c..........
      endif
c..........

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 7770)

      return
      end

c|||||||
#  endif
c|||||||

c-----------------------------------------------------------------------

c(((((((((((((((((((((((
#  if defined (CLEMCLIM)
c(((((((((((((((((((((((

      subroutine readdrive_clem (iu, cfilin, fielddrive, cvarin, nmonin,
     *                           ifanom, timeice)

c     Reads mid-month driving fields from Clemens climate file to array
c     fielddrive, for months 1 to nmonin. nmonin = nmon for all calls, 
c     except 1 for topog (cvarin ='orog'). 
c     ifanom = 0 for experiment, = 1 for anomaly ctrl (diagnostic only).

#include <comicephys.h>
#include <comicegrid.h>
#include <netcdf.inc>

      character*(*) cfilin, cvarin
      dimension fielddrive (nlong,nlatg,nmonin)

      integer*4 ncid, rcode, varid, start(4), count(4)
      character cfile*240, cvar*16

      real*4 arrin(nlong,nlatg,12), brrin(nlong,nlatg,12)
      real*4 crrin(nlong,nlatg)
      dimension work(nlong), work2(nlong,nlatg)

      logical firstclem
      data firstclem /.true./
      save firstclem

      cfile = cfilin
      cvar = cvarin

c     calculate time indices (file has all months, all years
c     as the single time index):
      iyearcal = 2001 + int(timeice) 
      iyearcal = iyearcal + 1    ! because this subr is called only
c                                ! once a year, the last timestep 
c                                ! of the year, so most accurate if 
c                                ! we read the *next* year of the gcm 

c     iyearfil = max (1, min (295, iyearcal-2006+1))  ! 2006 to 2300
      iyearfil = max (1, min (195, iyearcal-2006+1))  ! 2006 to 2200

#if defined (CRHINV)
c     if inverse run, just use perpetual modern (GCM year 2006 AD):
      iyearfil = 1
#endif

      ima = 1 + (iyearfil-1)*12

c     write (ioterm,'(3a,2i6)')
c    *  'reading existing Netcdf ', cfile(1:lenchr(cfile)),
c    *  ' data file', ifanom, ima

c        Open existing netcdf file and read requested field

      rcode = nf_open(cfile, nf_nowrite, ncid)
      call checkrcode (rcode, cfile, 'open', 6000)

c        If first call, read latitude and longitude grids into 
c        alatg, along, in common (comicegrid.h). Convert to radians.
c        nb: lon and lat are doubles on nc file.

c------------------------
      if (firstclem) then
c------------------------

c          latitudes:

        start(1) = 1
        count(1) = nlatg
        start(2) = 1        
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lat', varid)
        call checkrcode (rcode, cfile, cvar, 7720)
        rcode = nf_get_vara_double (ncid, varid, start, count, alatg)
        call checkrcode (rcode, cfile, cvar, 7721)
        do j=1,nlatg
          alatg(j) = alatg(j) * pi/180.
        enddo

c       protect interpolation logic in setinterp (needed?):
        alatg(1)     = max (alatg(1),     -89.999 *pi/180.) 
        alatg(nlatg) = min (alatg(nlatg),  89.999 *pi/180.) 

c          longitudes:

        start(1) = 1
        count(1) = nlong
        start(2) = 1        
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lon', varid)
        call checkrcode (rcode, cfile, cvar, 7730)
        rcode = nf_get_vara_double (ncid, varid, start, count, work)
        call checkrcode (rcode, cfile, cvar, 7731)
c       shift cut from greenwich to dateline for ice model:
        do i=1,nlong
          ii = mod (i + nlong/2 -1, nlong)  + 1
          along(i) = work(ii) * pi/180.
          if (along(i).gt.pi) along(i) = along(i) - 2.*pi
        enddo

c       protect interpolation logic in setinterp (needed?):
        if ( along(1).gt. 179.999*pi/180. .or. 
     *       along(1).lt.-179.999*pi/180. ) 
     *    along(1) = -179.999 * pi/180.
        if ( along(nlong).gt. 179.999*pi/180. .or. 
     *       along(nlong).lt.-179.999*pi/180. ) 
     *    along(nlong) = 179.999 * pi/180.

        firstclem = .false.

c----------
      endif
c----------

c.............................
      if (cvar.eq.'orog') then
c.............................

c        topography:

      start(1) = 1
      count(1) = nlong
      start(2) = 1
      count(2) = nlatg
      start(3) = 1        
      count(3) = 1
      rcode = nf_inq_varid (ncid, 'orog', varid)
      call checkrcode (rcode, cfile, 'orog', 7735)
      rcode = nf_get_vara_real (ncid, varid, start, count, crrin)
      call checkrcode (rcode, cfile, 'orog', 7736)
      do j=1,nlatg
        do i=1,nlong
c         shift cut from greenwich to dateline for ice model:
          ii = mod (i + nlong/2 -1, nlong)  + 1
          fielddrive(i,j,1) = crrin(ii,j)               ! already elev
        enddo
      enddo

c.........
      else
c.........


c        Read field (12 months)  

      start(1) = 1
      count(1) = nlong
      start(2) = 1
      count(2) = nlatg
      start(3) = ima        
      count(3) = 12        

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfile, cvar, 7760)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfile, cvar, 7761)

      do m=1,nmonin
        do j=1,nlatg
          do i=1,nlong

c           shift cut from greenwich to dateline for ice model:
            ii = mod (i + nlong/2 -1, nlong)  + 1

            zrrin = arrin(ii,j,m)

            if (cvar.eq.'tas') then
              fielddrive(i,j,m) = zrrin - tmelt             ! deg K to C
            else if (cvar.eq.'pr') then
              fielddrive(i,j,m) = zrrin*86400.            ! mm/s to mm/d
            else 
              write (ioterm,'(/2a)')
     *          '*** Error (readdrive_clem): unknown fieldname: ', cvar
              stop
            endif

          enddo
        enddo
      enddo

c..........
      endif
c..........

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 7770)

      return
      end

c(((((((
#  endif
c(((((((

c-----------------------------------------------------------------------

c))))))))))))))))))))))
#  if defined (CSMORDO)
c))))))))))))))))))))))

      subroutine readdrive_csmordo (iu,cfilin,fielddrive,cvarin,nmonin,
     *                              timeice)

c     Reads mid-month driving fields from CSMORDO climate file into arr
c     fielddrive, for months 1 to nmonin. nmonin = nmon for all calls, 
c     except 1 for topog (cvarin ='PHIS'). 

#include <comicephys.h>
#include <comicegrid.h>
#include <netcdf.inc>

      character*(*) cfilin, cvarin
      dimension fielddrive (nlong,nlatg,nmonin)

      integer*4 ncid, rcode, varid, start(4), count(4)
      character cfile*240, cfilegrid*240, cfilevar*240, cvar*16

      real*4 arrin(nlong,nlatg,12)
      dimension work(nlong)

      logical firstcsmordo
      data firstcsmordo /.true./
      save firstcsmordo

      cfile = cfilin      ! just the path to datafiles
      cvar = cvarin

      cfilegrid = cfile(1:lenchr(cfile)) // 
     *            'b30.153.cam2.GRID.T31.nc'

      if (cvar.eq.'LANDFRAC' .or. cvar.eq. 'PHIS') then
        cfile = cfile(1:lenchr(cfile)) // 
     *          'b30.153.cam2.GRID.T31.nc'

      else if (cvar.eq.'Tlowest') then
        cfile = cfile(1:lenchr(cfile)) // 
#    if defined (ORDO15X)
     *          'b30.153.cam2.Tlowest.1870-1889.AC.nc'
#    elif defined (ORDO4X)
     *          'b30.152.cam2.Tlowest.1840-1859.AC.nc'
#    endif

      else if (cvar.eq.'PRECT') then
        cfile = cfile(1:lenchr(cfile)) //
#    if defined (ORDO15X)
     *          'b30.153.cam2.PRECT.1870-1889.AC.nc'
#    elif defined (ORDO4X)
     *          'b30.152.cam2.PRECT.1840-1859.AC.nc'
#    endif
      else
        write (ioterm,'(/3a)') '*** Error (readdrive_csmordo):',
     *                         ' unknown fieldname ', cvarin
      endif

c        If first call, read latitude and longitude grids into 
c        alatg, along, in common (comicegrid.h). Convert to radians.
c        nb: lon and lat are doubles on nc file.

c---------------------------
      if (firstcsmordo) then
c---------------------------

c        Open CSMORDO netcdf "grid" file and read requested field

        write (ioterm,'(3a)')
     *    'reading existing Netcdf ', cfilegrid(1:lenchr(cfilegrid)),
     *    ' data file'

        rcode = nf_open(cfilegrid, nf_nowrite, ncid)
        call checkrcode (rcode, cfilegrid, 'open', 6000)

c          latitudes:

        start(1) = 1
        count(1) = nlatg
        start(2) = 1        
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lat', varid)
        call checkrcode (rcode, cfilegrid, cvar, 7720)
        rcode = nf_get_vara_double (ncid, varid, start, count, alatg)
        call checkrcode (rcode, cfilegrid, cvar, 7721)
        do j=1,nlatg
          alatg(j) = alatg(j) * pi/180.
        enddo

c       protect interpolation logic in setinterp (needed?):
        alatg(1)     = max (alatg(1),     -89.999 *pi/180.) 
        alatg(nlatg) = min (alatg(nlatg),  89.999 *pi/180.) 

c          longitudes:

        start(1) = 1
        count(1) = nlong
        start(2) = 1        
        count(2) = 1
        rcode = nf_inq_varid (ncid, 'lon', varid)
        call checkrcode (rcode, cfilegrid, cvar, 7730)
        rcode = nf_get_vara_double (ncid, varid, start, count, work)
        call checkrcode (rcode, cfilegrid, cvar, 7731)
c       shift cut from greenwich to dateline for ice model:
        do i=1,nlong
          ii = mod (i + nlong/2 -1, nlong)  + 1
          along(i) = work(ii) * pi/180.
          if (along(i).gt.pi) along(i) = along(i) - 2.*pi
        enddo

c       protect interpolation logic in setinterp (needed?):
        if ( along(1).gt. 179.999*pi/180. .or. 
     *       along(1).lt.-179.999*pi/180. ) 
     *    along(1) = -179.999 * pi/180.
        if ( along(nlong).gt. 179.999*pi/180. .or. 
     *       along(nlong).lt.-179.999*pi/180. ) 
     *    along(nlong) = 179.999 * pi/180.

        rcode = nf_close (ncid)
        call checkrcode (rcode, cfilegrid, 'close', 7770)

        firstcsmordo = .false.

c----------
      endif
c----------

c        Open existing CSMORDO netcdf file and read requested field

      write (ioterm,'(3a,i6)')
     *  'reading existing Netcdf ', cfile(1:lenchr(cfile)),
     *  ' data file'

      rcode = nf_open(cfile, nf_nowrite, ncid)
      call checkrcode (rcode, cfile, 'open', 6000)

c        Read field 

      start(1) = 1
      count(1) = nlong
      start(2) = 1
      count(2) = nlatg
      start(3) = 1        
      count(3) = nmonin

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfile, cvar, 7760)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfile, cvar, 7761)

      do m=1,nmonin
        do j=1,nlatg
          do i=1,nlong

c           shift cut from greenwich to dateline for ice model:
c           ii = mod (i + nlong/2 -1, nlong)  + 1
            ii = i ! no shift for Ordovician fields

            zrrin = arrin(ii,j,m)

            if (cvar.eq.'Tlowest') then
              fielddrive(i,j,m) = zrrin - tmelt             ! deg K to C
            else if (cvar.eq.'PRECT') then
              fielddrive(i,j,m) = zrrin*1000.*86400.       ! m/s to mm/d
            else if (cvar.eq.'PHIS') then
              fielddrive(i,j,m) = zrrin/9.80616         ! geopot to elev
            else
              fielddrive(i,j,m) = zrrin
            endif

          enddo
        enddo
      enddo

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 7770)

      return
      end

c)))))))
#  endif
c)))))))

c-----------------------------------------------------------------------

c****
#else
c****
      subroutine climate_gcm ()
      return
      end
c*****
#endif
c*****

c----------------------------------------------------------------------

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#if defined (RACMOCLIM) || defined (RACMOANTSEAS) || ( ( defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM) ) && defined (GREENLAND) )
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      subroutine getracmo_gre_fields (togobs, tmonobs, pmonobs,
     *                                runofobs)

c     Separate subr to read RACMO fields and interpolate from RACMO to 
c     ice grid, for Greenland.
c     Uses nlonrac,nlatrac dimensions, so can be different from 
c     nlong,nlatg for CSMCLIM. 
c     (nb: nlong,nlatg,nlonrac,nlatrac are in comicegrid.h).

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

      dimension
     *  togobs(nx,ny), tmonobs(nmon,nx,ny), pmonobs(nmon,nx,ny)
#if defined (RACMOCLIM)
     * ,runofobs(nx,ny)
#endif

      dimension 
     *  tograc(nlonrac,nlatrac),
     *  tarac(nlonrac,nlatrac,nmon),
     *  prrac(nlonrac,nlatrac,nmon)
#if defined (RACMOCLIM)
     * ,rurac(nlonrac,nlatrac)
#endif

      dimension
     *  indv(nx,ny,4), jndv(nx,ny,4),  weiv(nx,ny,4)
      save indv, jndv, weiv

      character*240 cfile

c        Open existing RACMO netcdf file and read fields
c        (alondrac, alatdrac are in common)

      cfile = cpathhome(1:lenchr(cpathhome)) // 
c    *        'Data_greenland/RACMO2_mo_avg_climatology_1960_2012.nc'
     *        'Data_greenland/' //
     *        'RACMO2_mo_avg_climatology_1960_2012_withrunoff.nc'

      call readdrive_racmo (iuto,cfile,tograc(1,1),'zs',1,
     *                      nlonrac, nlatrac)
      call readdrive_racmo (iuto,cfile,alondrac(1,1),'lon',1,
     *                      nlonrac, nlatrac)
      call readdrive_racmo (iuto,cfile,alatdrac(1,1),'lat',1,
     *                      nlonrac, nlatrac)
      call readdrive_racmo (iuto,cfile,tarac(1,1,1),'t2m',  nmon,
     *                      nlonrac, nlatrac)
      call readdrive_racmo (iuto,cfile,prrac(1,1,1),'precip',nmon,
     *                      nlonrac, nlatrac)
#if defined (RACMOCLIM)
      call readdrive_racmo (iuto,cfile,rurac(1,1),'runoff',1,
     *                      nlonrac, nlatrac)
#endif

c        Change units

      do im=1,nmon
        do j=1,nlatrac
          do i=1,nlonrac
            tarac(i,j,im) = tarac(i,j,im) + tmelt           ! deg C to K
            prrac(i,j,im) = prrac(i,j,im) *.001*365.    ! mm/day to m/yr
          enddo
        enddo
      enddo

c        Set indices and weights (in common) from RACMO to ice grid

      call setinterp_racmo (alondrac, alatdrac, nlonrac, nlatrac,
     *                      alond, alatd, nx, ny, 
     *                      indv, jndv, weiv)

c        Bilinearly interpolate fields from RACMO to ice grid 

      do j=1,ny
        do i=1,nx
          togobs(i,j) = 0.
#if defined (RACMOCLIM)
          runofobs(i,j) = 0.
#endif
          do m=1,4
            togobs(i,j) = togobs(i,j) 
     *                  + tograc(indv(i,j,m),jndv(i,j,m)) * weiv(i,j,m)
#if defined (RACMOCLIM)
            runofobs(i,j) = runofobs(i,j) 
     *                  + rurac(indv(i,j,m),jndv(i,j,m)) * weiv(i,j,m)
#endif
          enddo   
          do im=1,nmon
            tmonobs(im,i,j) = 0.
            pmonobs(im,i,j) = 0.
            do m=1,4
              tmonobs(im,i,j) = tmonobs(im,i,j) 
     *                 + tarac(indv(i,j,m),jndv(i,j,m),im) * weiv(i,j,m)
              pmonobs(im,i,j) = pmonobs(im,i,j) 
     *                 + prrac(indv(i,j,m),jndv(i,j,m),im) * weiv(i,j,m)
            enddo   
          enddo   
        enddo   
      enddo   

      return
      end

c-----------------------------------------------------------------------

      subroutine getracmo_ant_fields (togobs, tmonobs, pmonobs, 
     *                                runofobs, hsbal)

c     Separate subr to read RACMO fields and interpolate from RACMO to 
c     ice grid, for Antarctica.
c     Uses nlonrac,nlatrac dimensions, so can be different from 
c     nlong,nlatg for CSMCLIM. 
c     (nb: nlong,nlatg,nlonrac,nlatrac are in comicegrid.h).
c     Topog is not on these RACMO Antarcitc files, so set to
c     modern observed (Bedmap2), passed as hsbal from initphys).

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>

      dimension
     *  togobs(nx,ny), tmonobs(nmon,nx,ny), pmonobs(nmon,nx,ny),
     *  hsbal(nx,ny)
#if defined (RACMOCLIM)
     * ,runofobs(nx,ny)
#endif

      dimension 
     *  tarac(nlonrac,nlatrac,nmon),
     *  prrac(nlonrac,nlatrac,nmon)

      dimension
     *  indv(nx,ny,4), jndv(nx,ny,4),  weiv(nx,ny,4)
      save indv, jndv, weiv

      character*240 cfile

      dimension ndaypm(nmon)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm

c        Open existing RACMO netcdf file and read fields
c        (alondrac, alatdrac are in common)

      cfile = cpathhome(1:lenchr(cpathhome)) // 
     *        'Data_racmo_ant/' //
     *        'RACMO2_TS2.nc'
      call readdrive_racmo (iuto,cfile,alondrac(1,1),'lon',1,
     *                      nlonrac, nlatrac)
      call readdrive_racmo (iuto,cfile,alatdrac(1,1),'lat',1,
     *                      nlonrac, nlatrac)
      call readdrive_racmo (iuto,cfile,tarac(1,1,1),'t2m',  nmon,
     *                      nlonrac, nlatrac)
      cfile = cpathhome(1:lenchr(cpathhome)) // 
     *        'Data_racmo_ant/' //
     *        'RACMO2_CUMPRECIP.nc'
      call readdrive_racmo (iuto,cfile,prrac(1,1,1),'precip',nmon,
     *                      nlonrac, nlatrac)

c        Change units

      do im=1,nmon
        do j=1,nlatrac
          do i=1,nlonrac
            tarac(i,j,im) = tarac(i,j,im) + tmelt     ! deg C to deg K
            prrac(i,j,im) = prrac(i,j,im) *.001
     *                      * (365./ndaypm(im))       ! mm/month to m/yr
          enddo
        enddo
      enddo

c        Set indices and weights (in common) from RACMO to ice grid

      call setinterp_racmo (alondrac, alatdrac, nlonrac, nlatrac,
     *                      alond, alatd, nx, ny, 
     *                      indv, jndv, weiv)

c        Bilinearly interpolate fields from RACMO to ice grid 

      do j=1,ny
        do i=1,nx
          togobs(i,j) = hsbal(i,j) ! from modern obs (subr initphys)
          do im=1,nmon
            tmonobs(im,i,j) = 0.
            pmonobs(im,i,j) = 0.
            do m=1,4
              tmonobs(im,i,j) = tmonobs(im,i,j) 
     *                 + tarac(indv(i,j,m),jndv(i,j,m),im) * weiv(i,j,m)
              pmonobs(im,i,j) = pmonobs(im,i,j) 
     *                 + prrac(indv(i,j,m),jndv(i,j,m),im) * weiv(i,j,m)
            enddo   
          enddo   
        enddo   
      enddo   

#if defined (RACMOCLIM)
c     RACMO2 runoff not available for Antarctica:
      call zero (runofobs, nx*ny)
#endif

      return
      end

c-----------------------------------------------------------------------

      subroutine readdrive_racmo (iu, cfile, fielddrive, cvarin, nmonin,
     *                            nlonin, nlatin)

c     Reads mid-month driving fields from RACMO modern climatology for
c     Greenland (netcdf file) into array fielddrive, for months
c     1 to nmonin. nmonin = nmon for all calls, except 1 for topog
c     fields (ckey='zs'). 
c
c     Dimensions (nlonin,nlatin) are passed, so can use for 
c     (-DCSMANOM && -DGREENLAND (for RACMO dims, different from CSM's
c     nlong, nlatg in common)

#include <comicephys.h>
#include <comicegrid.h>
#include <netcdf.inc>

      character*(*) cfile, cvarin
      dimension fielddrive (nlonin,nlatin,nmonin)

      integer*4 ncid, rcode, varid, start(4), count(4)
      character*16 cvar

      real*4 arrin(nlonin,nlatin,12)

      dimension ndaypm(12)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      save ndaypm

      cvar = cvarin

c        Open existing RACMO netcdf file and read requested field

      rcode = nf_open(cfile, nf_nowrite, ncid)
      call checkrcode (rcode, cfile, 'open', 6000)

c        Define starting/ending position eof 2-D data values in array

      start(1) = 1
      count(1) = nlonin

      start(2) = 1
      count(2) = nlatin

#if defined (EISANTA)
c     RACMO2 Antarctic files have additional dimension "height"
      start(3) = 1        
      count(3) = 1
      start(4) = 1        
      count(4) = nmonin
#else
      start(3) = 1        
      count(3) = nmonin
      if (cvar.eq.'runoff') count(3) = 12
#endif

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfile, cvar, 7760)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfile, cvar, 7761)

      if (cvar.eq.'runoff') then
        call zero (fielddrive, nlonin*nlatin*1)
        do j=1,nlatin
          do i=1,nlonin
            do m=1,12
              fielddrive(i,j,1) = fielddrive(i,j,1)
     *                        + arrin(i,j,m)*(1000./365.) ! m/yr to mm/d
     *                          *(ndaypm(m)/365.)         ! av. over mon
            enddo
          enddo
        enddo
      else
        do m=1,nmonin
          do j=1,nlatin
            do i=1,nlonin
              if (cvar.eq.'precip') then
c               fielddrive(i,j,m)=arrin(i,j,m)/ndaypm(m) !mm/mon to mm/d
                fielddrive(i,j,m)=arrin(i,j,m)*1000./365.  !m/yr to mm/d
              else if (cvar.eq.'t2m') then
                 fielddrive(i,j,m) = arrin(i,j,m) - tmelt   ! deg K to C
              else
                fielddrive(i,j,m) = arrin(i,j,m)
              endif
            enddo
          enddo
        enddo
      endif

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfile, 'close', 7762)

      return
      end

c-----------------------------------------------------------------------

      subroutine setinterp_racmo (alondg, alatdg, nxg, nyg,
     *                            alond, alatd, nx, ny,
     *                            indv, jndv, weiv)

c        Calculate indices and weights to interpolate from RCM (RACMO)
c        grid (nxg,nyg) to the ice grid (nx,ny).
c        Returns [indl,weil] from the ice grid into the RACMO "e-w" 
c        axis, and [indb,weib] from the ice grid into RACMO "n-s" axis.

      dimension
     *  alondg(nxg,nyg), alatdg(nxg,nyg),
     *  alond(nx,ny),    alatd(nx,ny),
     *  indv(nx,ny,4),   jndv(nx,ny,4),  weiv(nx,ny,4)

      dimension
     *  xcart(nx,ny),     ycart(nx,ny),    zcart(nx,ny),
     *  xcartg(nxg,nyg),  ycartg(nxg,nyg), zcartg(nxg,nyg),
     *  zdist(nxg,nyg),
     *  igv(4), jgv(4), dgv(4)

      parameter (pi = 3.14159265358979)

      zz = pi/180.
      do j=1,ny
        do i=1,nx
          xcart(i,j) = cos(alatd(i,j)*zz)*cos(alond(i,j)*zz)
          ycart(i,j) = cos(alatd(i,j)*zz)*sin(alond(i,j)*zz)
          zcart(i,j) = sin(alatd(i,j)*zz)
        enddo
      enddo

      do jg=1,nyg
        do ig=1,nxg
          xcartg(ig,jg) = cos(alatdg(ig,jg)*zz)*cos(alondg(ig,jg)*zz)
          ycartg(ig,jg) = cos(alatdg(ig,jg)*zz)*sin(alondg(ig,jg)*zz)
          zcartg(ig,jg) = sin(alatdg(ig,jg)*zz)
        enddo
      enddo

c================
      do j=1,ny
        do i=1,nx
c================

          do jg=1,nyg
            do ig=1,nxg
c             either great circle (exact):
c             zdist(ig,jg) = acos ( min (1.,  xcart(i,j)*xcartg(ig,jg)
c    *                                      + ycart(i,j)*ycartg(ig,jg)
c    *                                      + zcart(i,j)*zcartg(ig,jg)))
c             or [3-D straight line]^2 (approx):
              zdist(ig,jg) =   (xcart(i,j)-xcartg(ig,jg))**2
     *                       + (ycart(i,j)-ycartg(ig,jg))**2
     *                       + (zcart(i,j)-zcartg(ig,jg))**2
            enddo
          enddo

c-----------------
          do m=1,4
c-----------------
            zdistmin = 1.e20
            do jg=1,nyg
              do ig=1,nxg
                if (zdist(ig,jg).lt.zdistmin) then
                  zdistmin = zdist(ig,jg)
                  igv(m) = ig
                  jgv(m) = jg
                endif
              enddo
            enddo
            dgv(m) = sqrt(zdistmin)                     ! dist^2 to dist
            zdist(igv(m),jgv(m)) = 1.e20
c--------------
          enddo
c--------------

          dnorm = 0.1*(pi/180.)   ! RACMO grid size equiv to 0.1 deg lat

          sumdgv = 0.
          do m=1,4
            dgv(m) = exp (-dgv(m) / (0.3*dnorm))
            sumdgv = sumdgv + dgv(m)
          enddo

          do m=1,4
            indv(i,j,m) = igv(m) 
            jndv(i,j,m) = jgv(m)
            weiv(i,j,m) = dgv(m)/sumdgv
          enddo

c============
        enddo
      enddo
c============

      return
      end 

c^^^^^
#endif
c^^^^^

c-----------------------------------------------------------------------

      subroutine setinterp (along, alatg, nlong, nlatg,
     *                      alond, alatd, nx, ny,
     *                      indl, weil, indb, weib, rotate)

c        Calculate indices and weights to interpolate from the agcm 
c        grid (nlong,nlatg) to the ice grid (nx,ny). 
c        Returns [indl,weil] from the ice grid into the 
c        agcm longitude grid, and [indb,weib] from the ice grid into 
c        the agcm latitude grid.
c
c        along,alatg supplied in radians, alond,alatd in degrees.
c        Works for any range (cut) of along - just must increase 
c        eastward, and be uniformly spaced.
c
c        rotate (deg eastward, passed) is a uniform longitudinal
c        rotational shift (+ve eastward) from ice grid to Genesis gcm
c        grid. (Usually 0, but for calls for Genesis gcm, set in 
c        rotate_to_gcm, in comicegrid.h).
c
c        Outside "if defined (GCMMATRIX)" so can be used by bedrock_sl
c        and liu_proc.

      dimension
     *  along(nlong),  alatg(nlatg),
     *  alond(nx,ny),  alatd(nx,ny),
     *  indl(nx,ny),   weil(nx,ny),
     *  indb(nx,ny),   weib(nx,ny)

      parameter (pi=3.14159265358979, tpi=(2.*pi))

#if defined (EISLINE) && defined (SLMODEL)
c     For Gomez-Mitrovica sea-level model, 1-D domains, if grid sizes 
c     are identical, assume grids are identical (as in setinterp2):
      if (nlong.eq.nx .and. nlatg.eq.ny) then
        do j=1,ny
          do i=1,nx
            indl(i,j) = i
            weil(i,j) = 1.
            indb(i,j) = j
            weib(i,j) = 1.
          enddo
        enddo
        return
      endif
#endif

c        Set monotonically increasing longitudes (starting from 
c        along(1), regardless of range, handled below), assuming 
c        uniform dlong

      if (along(2).gt.along(1)) then 
        dlong = along(2) - along(1) 
      else
        dlong = along(2) + 2.*pi - along(1) 
      endif
      zlong1 = along(1) - 0.5*dlong

c     decide if input grid is global in longitude (so wraparound):
      ifwrap = 0
      if (dlong*nlong .ge. 0.999*tpi) ifwrap = 1

c        Bilinear interp...calculate indices into gcm grid

      do j=1,ny
        do i=1,nx
          zlat = alatd(i,j)*pi/180. 
          zlon = (alond(i,j) + rotate) * pi/180.

c         shift ice grid longitude into gcm range [zlong1,zlong1+2*pi]
          do ishif = 1,10
            if (zlon.lt.zlong1) then 
              zlon = zlon + 2.*pi
            else if (zlon.gt.zlong1+tpi) then 
              zlon = zlon - 2.*pi
            else
              goto 50
            endif
          enddo
          write (6,'(/a/a,3f10.2)') 
     *       '*** Should not get here: (setinterp, ishif loop)',
     *       '    zlon, zlong1=',
     *       zlon*180./pi, zlong1*180./pi
          stop
   50     continue

c         compute gcm longitude index assuming uniform spacing:
          val = (zlon-zlong1)/dlong
          ia = nint(val)
          wa = ia + 0.5 - val
          if (ifwrap.eq.1) then
            if (ia.eq.0) ia = nlong
          else
            ia = max (ia, 1)
          endif
          indl(i,j) = ia
          weil(i,j) = wa

          val = zlat
          if (val.le.alatg(1)) then
            indb(i,j) = 1
            weib(i,j) = 1.
          else if (val.ge.alatg(nlatg)) then
            indb(i,j) = nlatg
            weib(i,j) = 1.
          else
            do 32 jo = 2,nlatg
              if (val.lt.alatg(jo)) then
                indb(i,j) = jo-1
                weib(i,j) =(alatg(jo)-val) / (alatg(jo)-alatg(jo-1))
                goto 34
              endif
   32       continue
   34       continue
          endif

        enddo
      enddo

c: no longer used: LC now changed to S->N. (Also, used to cause a bug 
c  when called from readin for restartin file of different  resolution).
c#if defined (LOVECLIP)
cc    latitudes run N->S in Loveclim file, but S->N here (in alatg). 
cc    So switch indices accordingly. indb can be 0, need line in
cc    binlin interp code (climate_gcm) constraining jb. 
c     do j=1,ny
c       do i=1,nx 
c         indb(i,j) = nlatg - indb(i,j)            ! switch N->S to S->N
c         weib(i,j) = 1. - weib(i,j)               ! switch N->S to S->N
c       enddo
c     enddo
c#endif

      return
      end

c----------------------------------------------------------------------

      SUBROUTINE READDATA_ICE (IU, ARR, MLONA,MLATA, MLON,MLAT,
     *                         CKEY, ICODE, INUM)
C
C        Reads an agcm or lsx data input file, checking on label and
C        dimensions in header record. Data files are ascii, with 1
C        header record, 3 skip records (blank, longitude grid, blank),
C        and MLAT data records. Each data record contains a
C        latitude grid value, followed by MLON data values either
C        in A1 or F5.0 format (flagged by ICODE). The longitudes are
C        western hemisphere first, then eastern, ie, (I=MLON/2+1,MLON),
C        (I=1,MLON/2). The data records are in reverse latitude-grid 
C        order (J=MLAT,1,-1). Data values are returned in *real* array
C        ARR.
C
C        Permissible characters in A1 format are 0-9, blank(=0),
C        A-Z (=10-35), or a-z(=10-35).
C
C        If file has more than one field, each set of data records is
C        preceded by the same 4 records as above. Each group is
C        separated by an extra blank line, between the last data record
C        of the preceding field and the header record of the next field.
C        
C        IU    = unit number
C        ARR   = data array (returned)
C        MLONA = longitude dimension of ARR
C        MLATA = latitude  dimension of ARR
C        MLON  = longitude dimension of data
C        MLAT  = latitude  dimension of data
C        CKEY  = expected 8-byte label in word # 1, rec # 1
C        ICODE = 1 for A1 format, 2 for F5.0 format.
C        INUM  = field number for this file

      DIMENSION ARR(MLONA,MLATA)
      CHARACTER CKEY*(*)
      CHARACTER CRR(5000)*1, CKEYZ*8, COMMENT*100

      IF (MLON.GT.5000) THEN
        WRITE(6,'(/A)') 
     *    'Error (readdata_ice): data records too big (MLON < 5000)'
        STOP
      ENDIF
C
C        Initialize record number (just for error messages)
C
      IREC = (INUM-1)*(4+MLAT) + MAX(INUM-2,0)

      ckeyz = ' '
      mlonz = 0
      mlatz = 0
      comment = ' '
C
C        If past first field, skip blank line separator
C
      
      IF (INUM.GT.1) THEN
        IREC = IREC + 1
        READ(IU, *, ERR=900, IOSTAT=IOS)
      ENDIF
C
C        Read header record
C
      IREC = IREC + 1
      READ(IU, 10, ERR=900, IOSTAT=IOS) CKEYZ, MLONZ, MLATZ, COMMENT
   10 FORMAT(A8,2I8,8X,A)

      IF (CKEYZ.NE.CKEY .OR. MLONZ.NE.MLON .OR. MLATZ.NE.MLAT) THEN
         WRITE(6,12) ckey,         COMMENT(1:LENCHR(COMMENT)),
     *                  CKEYZ,MLONZ,MLATZ, CKEY//'        ',MLON,MLAT
   12    FORMAT(/' *** Fatal Error.',
     *           ' Incorrect attributes for data file ',A
     *          /'     File comment                        = ',A,
     *          /'     File     label, longitude, latitude = ',A8, 2I6
     *          /'     Expected label, longitude, latitude = ',A8, 2I6/)
         stop
      ENDIF
C
C        Skip 3 records (blank, longitude grid, blank)
C
      IREC = IREC + 3
      READ(IU, 50, ERR=900, IOSTAT=IOS)
   50 FORMAT(//)
C
      IF (ICODE.EQ.1) THEN
C
C          Read "A1-format" data records
C
        DO 100 J=MLAT,1,-1
          IREC = IREC + 1
c         READ (IU,102,ERR=900,IOSTAT=IOS) (CRR(I),I=MLON/2+1,MLON),
c    *                                     (CRR(I),I=1,MLON/2)
c         Don't change from dateline to greenwich for ice grid
          READ (IU,102,ERR=900,IOSTAT=IOS) (CRR(I),I=1,MLON)
  102     FORMAT(5X,3X,5000A1)
          DO 104 I=1,MLON
            IC = ICHAR(CRR(I)(1:1))
            IF (IC.EQ.ICHAR(' ')) THEN
              ARR(I,J) = 0
            ELSE IF (IC.GE.ICHAR('0') .AND. IC.LE.ICHAR('9')) THEN
              ARR(I,J) = IC - ICHAR('0')
            ELSE IF (IC.GE.ICHAR('A') .AND. IC.LE.ICHAR('Z')) THEN
              ARR(I,J) = IC - ICHAR('A') + 10
            ELSE IF (IC.GE.ICHAR('a') .AND. IC.LE.ICHAR('z')) THEN
              ARR(I,J) = IC - ICHAR('a') + 10
            ELSE
              GOTO 910
            ENDIF
  104     CONTINUE
  100   CONTINUE
C
      ELSE
C
C          Read "F5.0-format" data records
C
        DO 150 J=MLAT,1,-1
          IREC = IREC + 1
c         READ (IU,152,ERR=900,IOSTAT=IOS) (ARR(I,J),I=MLON/2+1,MLON),
c    *                                     (ARR(I,J),I=1,MLON/2)
c         Don't change from dateline to greenwich for ice grid
          READ (IU,152,ERR=900,IOSTAT=IOS) (ARR(I,J),I=1,MLON)
  152     FORMAT(5X,3X,5000F5.0)
  150   CONTINUE
C
      ENDIF

      RETURN

C        Fatal errors
C
  900 WRITE(6,902) ckey,         COMMENT, IREC, IOS
  902 FORMAT(/' *** Error reading data file ',A
     *       /'     File comment = ',A
     *       /'     record number=',I6
     *       /'     Iostat code  =',I6/)
      stop
C
  910 WRITE(6,912) ckey,         COMMENT, IREC, I, J, CRR(I), IC
  912 FORMAT(/' *** Bad character in data file ',A
     *       /'     File comment    = ',A
     *       /'     record number   =',I6
     *       /'     longitude index =',I6
     *       /'     latitude  index =',I6
     *       /'     bad character   = ',A1
     *       /'     ichar           = ',I6/)
      stop
C
      END
