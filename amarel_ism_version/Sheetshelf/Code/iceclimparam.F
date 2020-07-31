c-----------------------------------------------------------------------

#if defined (PARAMCLIM) || defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM)

      subroutine climate_param (tog, tmon, pmon, h, hs, hsbal,
     *                          ro18, dtantann, dtantjan, dtseas,
     *                          rco2, timeice)

c     Sets surface monthly surface climate (tmon,pmon) using non-GCM
c     parameterizations. For Antarctica, either SEARISECLIM 
c     (with modern datasets), SIMPANT, SIMPEO, or SIMPNHA.
 
#include <comicephys.h>
#include <comicegrid.h>

c     passed:
      dimension 
     *  tog(nx,ny),    
     *  tmon(nmon,nx,ny),    pmon(nmon,nx,ny),
     *  h(nx,ny),            hs(nx,ny),          hsbal(nx,ny)

c     local:
      dimension tann(nx,ny), pann(nx,ny), albeff(nx,ny)

#if defined (SEARISECLIM) || defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM) || ( defined (EISLINE) && defined FORCE34MYR )
      logical firstseariseclim
      data firstseariseclim /.true./
      dimension budgsearise(nx,ny), tasearise(nx,ny), hssearise(nx,ny)
      save firstseariseclim, budgsearise, tasearise, hssearise
#  if defined (BASACCUM)
      logical firstbas
      data firstbas /.true./
      dimension budgbas (nx,ny)
      save firstbas, budgbas
#  endif
#  if defined (SEARISEAR4)
      dimension tar4(nx,ny), tar4mod(nx,ny),
     *          par4(nx,ny), par4mod(nx,ny)
      logical firstar4
      data firstar4 /.true./
      save tar4mod, par4mod, firstar4
#  endif
#endif

#if  defined (RACMOANTSEAS)
      logical firstracant
      data firstracant /.true./
      dimension budgracant(nmon,nx,ny), taracant(nmon,nx,ny), 
     *  taracantav(nx,ny), hsracant(nx,ny)
      save firstracant, budgracant, taracant, taracantav, hsracant
#endif
   
      dimension ndaypm(nmon)
      data ndaypm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/ 
      save ndaypm

c=======================================================================
#if defined (SEARISECLIM) || defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM) || ( defined (EISLINE) && defined FORCE34MYR )
c=======================================================================

C++++++++++++++++++++++++++++++++++++++++++++++
#  if defined (EISLINE) && defined (FORCE34MYR)
C++++++++++++++++++++++++++++++++++++++++++++++
      do j=1,ny
        do i=1,nx
           budgsearise(i,j) = 0.15
           tasearise(i,j) = tmelt - 15.
           hssearise(i,j) = 0.
         enddo
       enddo
C++++++
#  else
C++++++

      if (firstseariseclim) then
#    if defined (NOLEBROCQ)
        call searise_aggreg (budgsearise, idum, 'presprcp')
        call searise_aggreg (tasearise,   idum, 'presartm')
        call searise_aggreg (hssearise,   idum, 'usrf')
#    else
c       call searise_aggreg (budgsearise, idum, 'acca')   ! Arthern
        call searise_aggreg (budgsearise, idum, 'accr')   ! van de Berg
        call searise_aggreg (tasearise,   idum, 'temp')
        call searise_aggreg (hssearise,   idum, 'usrf')
#    endif
        firstseariseclim = .false.
      endif

c+++++++
#  endif
c+++++++

#  if defined (BASACCUM)
      if (firstbas) then
        write (6,'(a)') 'Reading BAS modern surface accumulation data'  
        call readbas (alond, alatd, budgbas, nx, ny, iubas)
        firstbas = .false.
      endif
#  endif

#  if defined (SEARISEAR4)
      if (firstar4) then
        call searise_aggreg_ar4 (tar4mod, 'instant_artm', 0.)
        call searise_aggreg_ar4 (par4mod, 'instant_prcp', 0.)
        firstar4 = .false.
      endif

      call searise_aggreg_ar4 (tar4, 'instant_artm', timeice)
      call searise_aggreg_ar4 (par4, 'instant_prcp', timeice)
#  endif

#  if defined (WARMAMPA)
      zwarm = WARMAMPA
#  else
      zwarm = 1.
#  endif
#  if defined (COLDAMPA)
      zcold = COLDAMPA
#  else
      zcold = 1.
#  endif

c-------------------------
#  if defined (FORCE34MYR)
c-------------------------
c     zdtorb = dtantann                      ! c old and c jan.A
      zdtorb = dtantjan - 0.5*dtseas         ! = 0.5*(dtantjan+dtantjul)

      zari = 0.
      do j=1,ny
        do i=1,nx
          if (h(i,j).ne.0.) zari = zari + darea(i,j)
        enddo
      enddo

c     albedo fdbk for smaller-than-modern (vs. Bedmap2 tot ice area).  
c     For greater-than-modern area, assume handled by ro18 (Pleist). 
#    if defined (EISANTA) 
#      if defined (COEFFALB)
      zdtalb = COEFFALB * max (0., 1.-(zari/13.924e12)) 
#      else
      zdtalb = 10. * max (0., 1.-(zari/13.924e12)) 
#      endif
#    else
      zdtalb = 0
#    endif

      coeffo18 = -5.
c     coeffo18 = -10.

c     coeffco2 =  5.
      coeffco2 =  10.
c     coeffco2 =  12.
c     coeffco2 =  14.
c     coeffco2 =  16.
c------
#  else
c------
      zdtorb = dtantann
      zdtalb = 0.
      coeffo18 = -10.
      coeffco2 =  6.
c-------
#  endif
c-------

c  can override above:
#  if defined (COEFFO18)
      coeffo18 = COEFFO18
#  endif
#  if defined (COEFFCO2)
      coeffco2 = COEFFCO2
#  endif

      do j=1,ny
        do i=1,nx
c         nb: no topog correction (hs vs hssearise), done in iceclimall:
          zdt =   zdtorb
     *          + zdtalb
     *          + coeffo18*zcold*max(ro18,0.)
     *          + coeffo18*zwarm*min(ro18,0.)
     *          + coeffco2*(log(rco2)/log(2.))

c    *          + (-10.)*zcold*max(ro18,0.)
c    *          + (-10.)*zwarm*min(ro18,0.)
c    *          + 6.*(log(rco2)/log(2.))

#  if defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM) 
          if (timeice.ge.0.) zdt = 0.
#  endif

          tann(i,j) = tasearise(i,j) + zdt

#  if defined (KN_DTA)
c         increase air temperatures uniformly (Knut,Nick):
          tann(i,j) = tann(i,j) + KN_DTA
#  endif

#  if defined (BASACCUM)
          pann(i,j) = budgbas(i,j)
#  else
          pann(i,j) = budgsearise(i,j) 
#  endif

          if (tlapseprecipclim.ne.0.) 
     *      pann(i,j) = pann(i,j) * (2.**(zdt/tlapseprecipclim))

#  if defined (SEARISEAR4)
          tann(i,j) = tann(i,j) + tar4(i,j) - tar4mod(i,j)
          pann(i,j) = pann(i,j) 
     *                * max (-5., min(5., par4(i,j)/par4mod(i,j)))
#  endif
#  if defined (PRECIPMULT)
          pann(i,j) = pann(i,j)*PRECIPMULT
#  endif
#  if defined (PRECIPHOL)
          pann(i,j) = pann(i,j) *
     *    (1. + PRECIPHOL*max(0., min(1., 1.-abs(timeice+4000.)/4000.)))
#  endif

        enddo
      enddo

c======================
#elif defined (SIMPANT)
c======================

      do j=1,ny
        do i=1,nx
c---
c         EISMINT II:
          tann(i,j) = tmelt + 34.46 - .00914*hs(i,j)
     *                              - .68775*abs(alatd(i,j))
c---
c         Ritz et al (2001):
c         if (hs(i,j).le.200.) then
c           tann(i,j) = tmelt + 49.642 - .943*abs(alatd(i,j))
c
c         else if (hs(i,j).gt.200. .and. hs(i,j).le.1500) then
c           tann(i,j) = tmelt + 36.689 - .725*abs(alatd(i,j))
c    *                                 - .005102*hs(i,j)
c
c         else if (hs(i,j).ge.1500.) then
c           tann(i,j) = tmelt +  7.405 - .180*abs(alatd(i,j))
c    *                                 - .014285*hs(i,j)
c         endif
c---

          tann(i,j) = tann(i,j) + dtantann
c                               antarctic cooling for lgm:
     *                          + (-10.)*ro18
c                               antarc warming per CO2 doubling:
     *                          + 6.*(log(rco2)/log(2.))

          pann(i,j) = 1.5 * (2.**((tann(i,j)-tmelt)/10.))
        enddo
      enddo

c=====================
#elif defined (SIMPEO)
c=====================

c     ice-sheet albedo feedback
      call calcalbeff (h, albeff, timeice)

      do j=1,ny
        do i=1,nx

c>>>>>>>>>>>>>>>>>>>>>>>
#  if defined (SIMPEO_1)
c>>>>>>>>>>>>>>>>>>>>>>>

          tann(i,j) = tmelt  
     *                        - 4.
     *                        + 0.1*(90. - abs(alatd(i,j))) ! latitude
     *                        - .006*hs(i,j)                ! lapse rate
     *                        - 10.*albeff(i,j)             ! ice alb
     *                        + 1.*dtantann                 ! orbit
     *                        + 5.*(log(rco2)/log(2.))      ! co2

          pann(i,j) = 1.5 * (2.**((tann(i,j)-tmelt)/10.))   ! simpeo1

c>>>>>>>>>>>>>>>>>>>>>>>>>
#  elif defined (SIMPEO_2)
c>>>>>>>>>>>>>>>>>>>>>>>>>

          tann(i,j) = tmelt  
     *                        -  5.  ! -15, -13, -10, -8, -5

#if defined (LINEH2)
     *                        + 0.1*(90. - abs(alatd(i,j))) ! latitude
c    *                        + 0.5*(90. - abs(alatd(i,j))) ! latitude
#else
     *                        + 0.1*(90. - abs(alatd(i,j))) ! latitude
#endif

     *                        - .007*hs(i,j)                ! lapse rate

     *                        - 8.*albeff(i,j)              ! ice alb
c    *                        -16.*albeff(i,j)              ! ice alb
c    *                        - 4.                          ! no ice alb

     *                        + 1.*dtantann                 ! orbit

     *                        + 6.*(log(rco2)/log(2.))      ! co2

          pann(i,j) = 0.8 * (2.**((tann(i,j)-tmelt)/10.))   ! simpeo2,3

c>>>>>>>>>>>>>>>>>>>>>>>>>
#  elif defined (SIMPEO_3)
c>>>>>>>>>>>>>>>>>>>>>>>>>

          tann(i,j) = tmelt  
     *                        - 7.
     *                        + 0.1*(90. - abs(alatd(i,j))) ! latitude
     *                        - .008*hs(i,j)                ! lapse rate
     *                        - 10.*albeff(i,j)             ! ice alb
     *                        + 1.*dtantann                 ! orbit
     *                        + 6.*(log(rco2)/log(2.))      ! co2

          pann(i,j) = 0.8 * (2.**((tann(i,j)-tmelt)/10.))  ! simpeo2,3

c         orographic precip: increase over ice for steep slopes:
c         pann(i,j) = 0.6 * (2.**((tann(i,j)-tmelt)/10.))
c         if (h(i,j).gt.0.) then
c           ip = min (i+1,nx)
c           im = max (i-1,1)
c           jp = min (j+1,ny)
c           jm = max (j-1,1)
c           zslope = max ( abs((hs(ip,j)-hs(im,j))/(xh(ip,j)-xh(im,j))),
c    *                     abs((hs(i,jp)-hs(i,jm))/(yh(i,jp)-yh(i,jm)))
c    *                   )
c           pann(i,j) = pann(i,j)
c    *                * (1. + max (0.,min(3.,zslope/.01)) )
c         endif

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#  elif defined (SIMPEO_HUY93)
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

chuy      tann(i,j) = tmelt  
chuy *                        + 34.46
chuy *                        - 0.00914*hs(i,j)             ! lapse rate
chuy *                        - 0.68775*abs(alatd(i,j))     ! latitude

          tann(i,j) = tmelt  
csta:
c    *                        - 15.0 
c    *                        - 0.008*hs(i,j)               ! lapse rate
c    *                        + 0.10*(70.-abs(alatd(i,j)))  ! latitude

cstb:
c    *                        - 15.0 
c    *                        - 0.008*hs(i,j)               ! lapse rate
c    *                        + 0.70*(70.-abs(alatd(i,j)))  ! latitude

cstc,d:
     *                        -  5.0 
c    *                        -  10.0 
     *                        - 0.008*hs(i,j)               ! lapse rate
     *                        + 0.70*(70.-abs(alatd(i,j)))  ! latitude
c    *                        + 0.10*(70.-abs(alatd(i,j)))  ! latitude
     *                        - 10.*albeff(i,j)             ! ice alb

#    if defined (DTPERT)
     *                        + ( DTPERT )
#    endif

chuy      pann(i,j) =  0.78 
chuy *              + 2.525e-2*(tann(i,j)-tmelt)
chuy *              + 2.225e-4*(tann(i,j)-tmelt)**2

          pann(i,j) = 0.8 * (2.**((tann(i,j)-tmelt)/10.))

c>>>>>>>
#  else
c>>>>>>>

        write (ioterm,'(2a)') 
     *     'Error (climate_param): no SIMPEO type',
     *     ' (SIMPEO_[1,2,3,HUY93]) defined in makefile' 
        stop

c>>>>>>>
#  endif
c>>>>>>>

        enddo
      enddo

c======================
#elif defined (SIMPNHA)
c======================

c     ice-sheet albedo feedback
      call calcalbeff (h, albeff, timeice)

      do j=1,ny
        jm1=max (j-1,1)
        jp1=min (j+1,ny)
        do i=1,nx
          im1=max (i-1,1)
          ip1=min (i+1,nx)

          tann(i,j) = tmelt 
     *                    - 15.0 
c    *                    - .0080*hs(i,j)                   ! lapse rate
     *                    - .0050*hs(i,j)                   ! lapse rate

c    *                    + 0.50*(70.-abs(alatd(i,j))) + 5. ! latitude  
c
c    *                    + 0.80*(70.-abs(alatd(i,j))) - 3. ! cycb
c    *                    + 0.80*(70.-abs(alatd(i,j))) - 7. ! cycc
c    *                    + 0.80*(70.-abs(alatd(i,j))) - 5. ! cycd
c    *                    + 0.80*(70.-abs(alatd(i,j))) - 3. ! cyce
     *                    + 0.80*(70.-abs(alatd(i,j))) - 4. ! cycf,g,h

c    *                    - 0.25*(alond(i,j)-(-80.))        ! cycb-g
     *                    - 0.15*(alond(i,j)-(-90.)) + 4.   ! cych
c    *                    - 0.15*(alond(i,j)-(-90.)) + 6.   ! cyci
c    *                    - 0.15*(alond(i,j)-(-90.)) + 2.   ! cycj

     *                    - 15.*albeff(i,j)                 ! ice alb

c                         nha cooling for lgm:
c    *                    + (-10.)*ro18
c                         nha warming per CO2 doubling:
c    *                    + 4.*(log(rco2)/log(2.))

#    if defined (DTPERT)
c    *                    + (DTPERT) * sin(2.*pi*timeice/40.e3)
     *                    - (DTPERT) * sin(2.*pi*timeice/40.e3)
c    *                    - (DTPERT) * 0.5
#    endif

          zix = (hs(ip1,j)-hs(im1,j)) / (2.*dx(i,j))       ! orog precip
          ziy = (hs(i,jp1)-hs(i,jm1)) / (2.*dy(i,j))
          zis = sqrt(zix**2+ziy**2)
          zporog = min (3., zis/.01)

c         pann(i,j) = 0.5 * (2.**((tann(i,j)-tmelt)/10.)) * zporog
          pann(i,j) = 0.8 * (2.**((tann(i,j)-tmelt)/10.)) * zporog

        enddo
      enddo

c====
#else
c====

      write (ioterm,'(2a)') 
     *   'Error (climate_param): no climate type',
     *   ' [SEARISECLIM, SIMPANT, SIMPEO, SIMPNHA] defined in makefile' 
      stop

c=====
#endif
c=====

c        Set tmon, pmon (seasonal), tog

c**********************************************
#if defined (EISANTA) && defined (RACMOANTSEAS)
c**********************************************

c        Add seasonal temp departures from RACMO, keep annual from above

      if (firstracant) then
c       (hsracant is not used):
        call getracmo_ant_fields (hsracant, taracant, budgracant, 
     *                            zrudum, hsbal)
        do j=1,ny
          do i=1,nx
            taracantav(i,j) = 0.
            do m=1,nmon
              taracantav(i,j) = taracantav(i,j) 
     *                        + taracant(m,i,j)*ndaypm(m)/365.
            enddo
          enddo
        enddo
        firstracant = .false.
      endif

      do j=1,ny
        do i=1,nx
          do m=1,nmon
c           tmon(m,i,j) = taracant(m,i,j)
c           pmon(m,i,j) = budgracant(m,i,j)
            tmon(m,i,j) = taracant(m,i,j) 
     *                  - taracantav(i,j)  + tann(i,j)
            pmon(m,i,j) = pann(i,j)
          enddo

c         set tog to be consistent with tmon:
#  if defined (SEARISECLIM)
          tog(i,j) = hssearise(i,j)
#  else
          tog(i,j) = hs(i,j)
#  endif
        enddo
      enddo

c****
#else
c****

c        Add sinusoidal temp cycle

      do j=1,ny
        do i=1,nx
#if defined (SIMPEO_HUY93)
c         Huybrechts (1993) TMS-TMA:
chuy      thamp =   (16.81 - 34.46)
chuy *            - (0.00692 - 0.00914)*hs(i,j)
chuy *            - (0.27937 - 0.68775)*abs(alatd(i,j))
chuy *          + 0.5*dtseas
          zwei = max (0., min (1., hs(i,j)/3000.))
          thamp = 0.5 * ( (1.-zwei)*20. + zwei*30. + dtseas )
#elif defined (EISANTA)
c         zonal modern Jan-Jul TS2 diff. is ~20 C at 70S, ~30 C 80-90S
c         (Genesis.3.0/Runpresilion/ts2jmj.cgm). Param vs elevation.
c         Also add orbital effect minus modern (dtseas).
          zwei = max (0., min (1., hs(i,j)/3000.))
          thamp = 0.5 * ((1.-zwei)*20. + zwei*30.)
#  if defined (FORCE34MYR)
          thamp = thamp 
c    *          + 0.5*dtseas + dtantjan - dtantann ! c old
c    *          + dtantjan - dtantann              ! c jan.A
     *          + 0.5*dtseas                 ! = 0.5*(dtantjan-dtantjul)
#  else
          thamp = thamp + 0.5*dtseas
#  endif
#else
          thamp = 0.5 * (20. + dtseas)
#endif

          if (alatd(i,j).ge.0.) thamp  = -thamp    ! SH vs NH phase 7/16

          do m=1,nmon
            tmon(m,i,j) = tann(i,j) + thamp * cos(2.*pi*(m-1.)/nmon)
            pmon(m,i,j) = pann(i,j)
          enddo

c         set tog to be consistent with tmon:
#  if defined (SEARISECLIM) || ( defined (EISLINE) && defined FORCE34MYR)
          tog(i,j) = hssearise(i,j)
#  else
          tog(i,j) = hs(i,j)
#  endif
        enddo
      enddo

c*****
#endif
c*****

      return
      end

#endif

c-----------------------------------------------------------------------

      subroutine calcalbeff (h, deflect, timeice)

c     For climate_simp, calculates spatial effect of presence of ice
c     on temperature (albedo effect).  For ALBEFFGREEN, similar to 
c     lithospheric Green's function calc in subr bedrock: sum a local 
c     decaying pattern (exp(-|x|/L) or exp(-(x/L)**2)) for each point
c     with ice, with length scale albefflen. Returns net effect 
c     (max-= ~1) in array deflect. (Effectively ~same as the fraction
c     of distance-weighted ice area within a circle around each point. 
c     For ALBEFFTOT, simly a uniform cooling proprtional to total ice
c     area (like in iceclimgcm matrix).

#include <comicephys.h>
#include <comicegrid.h>

c     h = ice-sheet thickness (m) (supp)
c     deflect = cooling effect of ice surfaces (0 to ~1) (ret)

      dimension 
     *  h(nx,ny), deflect(nx,ny)

c########################
#if defined (ALBEFFGREEN)
c########################

c     for parallelization:
#  if defined (EISLINE)
      parameter (nseg=1)
#  else
      parameter (nseg=4)
#  endif
      dimension deflectseg(nx,ny,nseg)

      save deflectnorm

      parameter (iu=0)   ! for no diagnostic printout
c     parameter (iu=159) ! for diagnostic printout (fort.159)

c     set e-folding length scale (albeffbox, km), and maximum distance
c     for calcs (albeffbox, km):
      parameter ( albefflen = ALBEFFGREEN,  ! 50.e3, 200.e3, 500.e3
     *            albeffbox = 5.*albefflen)

#  if defined (LONLAT)
      dimension xcart(nx,ny), ycart(nx,ny), zcart(nx,ny)
      save xcart, ycart, zcart
#  elif defined (STEREO)
      parameter (nbox=albeffbox/dd0)
      dimension abox(-nbox:nbox,-nbox:nbox) 
      save abox
#  endif

      logical first
      save first
      data first /.true./

c     statement function for effect vs. distance of ice cover:
c     falbeff(x) = exp (-abs(x)/albefflen)
cc    falbeff(x) = exp (-(x/albefflen)**2)
c     1 for x < albefflen, 0 for x >= albefflen:
      falbeff(x) = 0.5 * ((albefflen-x) + abs(albefflen-x))
     *           / max ( .001, abs(albefflen-x) )

c        If first call, set 3D cartesian (x,y,z) coordinates (*cart, 
c        saved) for calc of great-circle distances between points.

c====================
      if (first) then
c====================
#if defined (LONLAT)

        zz = pi/180.
        do j=1,ny
          do i=1,nx
            xcart(i,j) = cos(alatd(i,j)*zz)*cos(alond(i,j)*zz)
            ycart(i,j) = cos(alatd(i,j)*zz)*sin(alond(i,j)*zz)
            zcart(i,j) = sin(alatd(i,j)*zz)
          enddo
        enddo

#elif defined (STEREO)

c       set local effect array (abox)
        do ia = -nbox,nbox
          do ja = -nbox,nbox
            zx = sqrt(float(ia**2 + ja**2))*dd0
            abox(ia,ja) = falbeff (zx)
          enddo
          if (iu.ne.0) write (iu,'(100f8.2)')(abox(ia,ja),ja=-nbox,nbox)
        enddo

#endif

c         Calculate integrated falbeff * darea at center of box 
c         (sides 2*albeffbox), if all ice covered (deflectnorm, saved)

        nnorm = 50
        zdd = albeffbox/nnorm
#if defined (EISLINE)
        zda = zdd*dd0
        deflectnorm = 0.
        do i=1,nnorm
          zdx = sqrt ((i-.5)**2) * zdd
          deflectnorm = deflectnorm + zda*falbeff(zdx)
        enddo
        deflectnorm = 2.*deflectnorm  ! for 2 directions
#else
        zda = zdd**2
        deflectnorm = 0.
        do j=1,nnorm
          do i=1,nnorm
            zdx = sqrt ((j-.5)**2 + (i-.5)**2) * zdd
            deflectnorm = deflectnorm + zda*falbeff(zdx)
          enddo
        enddo
        deflectnorm = 4.*deflectnorm  ! for 4 quadrants
#endif
        if (iu.ne.0) write (iu, '(a,f10.4)') 'deflectnorm/area = ',
     *    deflectnorm / ((2.*albeffbox)**2)
            
        first = .false.
c==========
      endif
c==========

c          Loop over multiple latitude segments, saving deflections
c          due to each in deflectseg, for parallelization. 

        call zero (deflectseg, nx*ny*nseg)

#if defined (OMPALBEFF)
!$OMP PARALLEL DO
#if defined (LONLAT)
!$OMP* shared  (deflectseg, darea, xcart, ycart, zcart, dy)
!$OMP* private (iseg, ja,jb, i,j, jja, jjb, jj, idist, idirb, idir,
!$OMP*          iic, ii, zx)
#else
!$OMP* shared  (deflectseg, darea, abox)
!$OMP* private (iseg, ja,jb, i,j, iia,iib,jja,jjb, ii,jj)
#endif
!$OMP* schedule(dynamic,2)
#endif
c=====================
        do iseg=1,nseg
c=====================
          ja = (iseg-1)*(ny/nseg) + 1
          jb =  iseg   *(ny/nseg)
          if (iseg.eq.nseg) jb = ny
c~~~~~~~~~~~~~~~~~~~~
          do j=ja,jb
            do i=1,nx
c~~~~~~~~~~~~~~~~~~~~

              if (h(i,j).gt.2.) then
                zhwei = min (1., h(i,j)/20.)

c                   Add this point-load's deflection to deflect(i,j)

#if defined (LONLAT)

c               lat range equal to albeffbox (or pole):
                jbox = nint (albeffbox/dy((nx+1)/2,(ny+1)/2)) + 1
                jja = max (1,  j-jbox)
                jjb = min (ny, j+jbox)
                do jj=jja, jjb
c                 lon range outward from i until exceed albeffbox:
                  do idist=0,nx/2
                    idirb= 1
c                   only count i and i +/- nx/2 pts once:
                    if (idist.eq.0 .or. idist.eq.nx/2) idirb= -1
                    do idir=-1,idirb,2
                      iic = i + idir*idist
                      ii = mod (iic + nx - 1, nx) + 1
c                     either great circle (exact):
                      zx = acos ( min (1.,   xcart(i,j)*xcart(ii,jj)
     *                                     + ycart(i,j)*ycart(ii,jj)
     *                                     + zcart(i,j)*zcart(ii,jj)))
     *                     * radius
c                     or 3-D straight line (approx):
c                     zx = sqrt(  (xcart(i,j)-xcart(ii,jj))**2
c    *                          + (ycart(i,j)-ycart(ii,jj))**2
c    *                          + (zcart(i,j)-zcart(ii,jj))**2 )
c    *                     * radius
                      if (zx.gt.albeffbox) goto 500
                      deflectseg(ii,jj,iseg) = deflectseg(ii,jj,iseg)
     *                                + zhwei * darea(i,j) * falbeff(zx)
                    enddo
                  enddo
  500             continue
                enddo
#elif defined (STEREO)
c               polar stereographic:
                jja = max (1,  j-nbox)
                jjb = min (ny, j+nbox)
                iia = max (1,  i-nbox)
                iib = min (nx, i+nbox)
                do jj=jja,jjb
                  do ii=iia,iib
                    deflectseg(ii,jj,iseg) = deflectseg(ii,jj,iseg)
     *                            + zhwei * darea(i,j) * abox(ii-i,jj-j)
                  enddo
                enddo
#endif

              endif
c~~~~~~~~~~~~~~~~
            enddo
          enddo
c~~~~~~~~~~~~~~~~
c============
        enddo
c============

c          Sum up deflections from each segment (for parallelization)

        call zero (deflect, nx*ny)
        do iseg=1,nseg
          do j=1,ny
            do i=1,nx
              deflect(i,j) = deflect(i,j) + deflectseg(i,j,iseg)
            enddo
          enddo
        enddo

c          Normalize so that deflect = 1 if surrounded by extensive ice
c          (yields up to ~1.08...then take min (1, ...))

        do j=1,ny
          do i=1,nx
            deflect(i,j) = min (1., deflect(i,j)/deflectnorm)
          enddo
        enddo

        if (iu.ne.0) then
c         call printmap (timeice, h, 'ice thickness',   1., iu, 0)
          call printmap (timeice, h, 'ice thickness', 150., iu, 0)
          call printmap (timeice, deflect, 'albeff', 0.1, iu, 0)
        endif
        zmax = -1.e20
        do j=1,ny
          do i=1,nx
            zmax = max (zmax, deflect(i,j))
          enddo
        enddo
        if (iu.ne.0) write (iu,'(a,f10.5)') 'max deflect=', zmax

c########################
#elif defined (ALBEFFTOT)
c########################

c        Simple alternative: uniform value, proportional to ice area 
c        (like weiice in iceclimgcm for matrix):

      ztotari = 0.
      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.001) ztotari = ztotari + darea(i,j)
        enddo
      enddo
      ztotari = ztotari*1.e-6 ! m2 to km2
#  if defined (EISANTA)
      ztotari = min (1., ztotari/13.e6)        ! modern Antarctica
#  elif defined (NHA)
      ztotari =          ztotari/17.e6         ! LGM N.Amer. + Greenland
#  else
      ztotari = 0.
#  endif
      call resetr (deflect, nx*ny, ztotari)

c####
#else
c####

c       no albeff effect:
        call zero (deflect, nx*ny)

c#####
#endif
c#####

      return
      end

c-----------------------------------------------------------------------

#if defined (BASACCUM)

      subroutine readbas (alond, alatd, accum, nx, ny, iu)

c        Reads Arthern et al (BAS) Antarctic snow accumulation
c        data file. Transfers to model grid by nearest neighbor
c        using lon and lat.

#include <comicepath.h>

      dimension alond(nx,ny), alatd(nx,ny), accum(nx,ny)

c     input data:
      parameter (ninmax = 150000)
      dimension accumin(ninmax), xin(ninmax), yin(ninmax), zin(ninmax)

      parameter (tmelt = 273.16)
      parameter (pi = 3.14159265358979)

      character*200 cfilin

      cfilin = cpathhome(1:lenchr(cpathhome))
     *         // 'Data_bas/amsr_accumulation_map.txt'

c        Read BAS input file

      open (iu, file=cfilin, status='old')
      do  iskip=1,21
        read (iu,*)
      enddo
      nin = 0
  100 continue
      read (iu,*,end=200,err=200) zlat, zlon, izx, zy, zacc, zerr
      nin = nin + 1  
      if (nin.gt.ninmax) then
        write (6,'(a,2i6)') 'Error readbas, nin,ninmax=',nin,ninmax
        stop
      endif
      accumin(nin) = zacc
      xin(nin) = cos(zlat*pi/180.)*cos(zlon*pi/180.)
      yin(nin) = cos(zlat*pi/180.)*sin(zlon*pi/180.)
      zin(nin) = sin(zlat*pi/180.)
c     write (6,'(i6,3f10.2)') nin, zlat, zlon, zacc
      goto 100
  200 close (iu)
c     write (6,'(a,i6)') 'nin=',nin

c--------

c         Interpolate to model grid (nearest neighbor)

      do j=1,ny
        write (6,'(a,i3)') 'doing j=',j
        do i=1,nx
           x1 = cos(alatd(i,j)*pi/180.)*cos(alond(i,j)*pi/180.)
           y1 = cos(alatd(i,j)*pi/180.)*sin(alond(i,j)*pi/180.)
           z1 = sin(alatd(i,j)*pi/180.)
           zdistmin = 1.e20

           do k=1,nin
             zdist = 
c    *               sqrt ( 
     *                  (x1-xin(k))**2 + (y1-yin(k))**2 + (z1-zin(k))**2
c    *                    )
c            zdist = acos (x1*xin(k) + y1*yin(k) + z1*zin(k)) ! slower
             if (zdist.lt.zdistmin) then
               zdistmin = zdist
               kmin = k
             endif
           enddo

           if (nint(accumin(kmin)).eq.-999) then 
c            if null bas value, use crude parameterization vs lat 
             zhs = 0.  ! sea level
             zts = tmelt + 34.46 - .00914*zhs
     *                           - .68775*abs(alatd(i,j))
             accum(i,j) = 1.5 * (2.**((zts-tmelt)/10.))           ! m/yr
           else
             accum(i,j) = .001*accumin(kmin)          ! kg/m2/yr to m/yr
           endif
        enddo
      enddo

      return
      end

#endif

c-----------------------------------------------------------------------

#if defined (SEARISECLIM)

      subroutine searise_timeseries (sealev, dtantann, timeice)

#include <comicepath.h>
#include <netcdf.inc>

c     If first call, reads Searise 5 km Netcdf modern data file for 
c     sea-level and temperature time series (last 405 kyr, saved). 
c     Interpolates to current timeice, returns sealev, dtantann.

      parameter (ntemptim=4051, nsltim=406)

      real*4 temptim(ntemptim),  temptimval(ntemptim), 
     *       sltim(nsltim),      sltimval(nsltim)
       save  temptim,            temptimval,
     *       sltim,              sltimval

      character*200 cfilin
      integer*4 ncid, rcode, varid

      logical first
      save first
      data first /.true./

c--------------------
      if (first) then
c--------------------

        write (6,'(a)') 
     *    'reading existing Netcdf Searise data file for time series'

        cfilin = cpathhome(1:lenchr(cpathhome)) //
#if defined (NOLEBROCQ) 
     *           '/Data_searise/Antarctica_5km_withshelves_v4.nc'
#else
     *           '/Data_searise/Antarctica_5km_dev1.0.nc'
#endif

        rcode = nf_open(cfilin, nf_nowrite, ncid)
        call checkrcode (rcode, cfilin, 'open', 5000)

c       temperature time series (yr, deg C):
        rcode = nf_inq_varid (ncid,  'temptimes',  varid)
        call checkrcode (rcode, cfilin, 'temptimes', 5571)
        rcode = nf_get_var_real (ncid, varid, temptim)
        call checkrcode (rcode, cfilin, 'temptimes', 5572)

        rcode = nf_inq_varid (ncid, 'temp_time_series',  varid)
        call checkrcode (rcode, cfilin, 'temp_time_series', 5573)
        rcode = nf_get_var_real (ncid, varid, temptimval)
        call checkrcode (rcode, cfilin, 'temp_time_series', 5574)

c---

c       sea level time series (yr, m):
        rcode = nf_inq_varid (ncid,  'sealeveltimes',  varid)
        call checkrcode (rcode, cfilin, 'sealevtimes', 5575)
        rcode = nf_get_var_real (ncid, varid, sltim)
        call checkrcode (rcode, cfilin, 'sealevtimes', 5576)

        rcode = nf_inq_varid (ncid, 'sealevel_time_series',  varid)
        call checkrcode (rcode, cfilin, 'sealevel_time_series', 5577)
        rcode = nf_get_var_real (ncid, varid, sltimval)
        call checkrcode (rcode, cfilin, 'sealevel_time_series', 5578)

c---

        rcode = nf_close (ncid)
        call checkrcode (rcode, cfilin, 'close', 5900)

c---

c          Convert times to yr BP (neg)

        do i=1,ntemptim 
          temptim(i) = -temptim(i)
          temptimval(i) = temptimval(i) + 273.15                ! C to K
        enddo

        do i=1,nsltim 
          sltim(i) = -sltim(i)
        enddo

        first = .false.
c----------
      endif
c----------

      if (timeice.ge.temptim(1)) then
        dtantann = temptimval(1) - temptimval(1)
      else if (timeice.le.temptim(ntemptim)) then
        dtantann = temptimval(ntemptim) - temptimval(1)
      else
        do i=1,ntemptim-1
          if (timeice.le.temptim(i) .and. timeice.ge.temptim(i+1)) then 
            zwei = (timeice-temptim(i+1))/(temptim(i)-temptim(i+1))
            dtantann = zwei*temptimval(i) + (1.-zwei)*temptimval(i+1)
     *                 - temptimval(1)
            goto 100
          endif
        enddo
        write(6,'(a)')"Error: searise_timeseries_A - shouldn't get here"
        stop
      endif
  100 continue

      if (timeice.ge.sltim(1)) then
        sealev = sltimval(1) - sltimval(1)
      else if (timeice.le.sltim(nsltim)) then
        sealev = sltimval(nsltim) - sltimval(1)
      else
        do i=1,nsltim-1
          if (timeice.le.sltim(i) .and. timeice.ge.sltim(i+1)) then 
            zwei = (timeice-sltim(i+1))/(sltim(i)-sltim(i+1))
            sealev = zwei*sltimval(i) + (1.-zwei)*sltimval(i+1)
     *               - sltimval(1)
            goto 200
          endif
        enddo
        write(6,'(a)')"Error: searise_timeseries_B - shouldn't get here"
        stop
      endif
  200 continue

      return
      end

#endif

c-----------------------------------------------------------------------

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#if defined (SEARISEINIT) || defined (SEARISECLIM) || defined (RCMANOM) || defined (LOVEANOM) || defined (CSMANOM)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine searise_aggreg (arr, iarr, cvarin)

c     Reads Searise 5 km Netcdf modern data file for field name cvar,
c     Reads Searise 5 km Netcdf modern data file, newer Lebrocq version
c     by default, for field name cvar. Aggregates to ice model array, 
c     returned in arr(real) or iarr(int).
c     cvar = thk, topg[sd], acca/accr, temp, usrf, bheatflx_fox/shapiro.
c     Or different for NOLEBROCQ:
c     cvar = presprcp,presartm,usrf,bheatflx (real), or
c            maskwater,maskpres (->thkmask/mask) (int).
c     For topg[sd,sa], read topg and return std. dev. of topg values
c     (topgsd) within each model grid box.
c
c     Or different names for Lebrocq developmental dataset (-DLEBROCQ):
c     cvar = thk, topg, acca/accr, temp, usrf, bheatflx_fox/shapiro (r).

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension arr(nx,ny), iarr(nx,ny)
      character*(*) cvarin

      character cfilin*200, cvar*16, cvarmask*16
      integer*4 ncid, rcode, varid, start(3), count(3)

c     Searise input file attributes. [xcorner,ycorner] are coords (m)
c     of lower-left corner of lower-left (1,1) cell, relative to SP.
      parameter (nxsea=1200, nysea=1200, ddin=5.e3,
     *           xcorner=-2800000., ycorner=-2800000.)
      real*4    arrin(nxsea,nysea), lat(nxsea,nysea), lon(nxsea,nysea)
      integer*4 iarrin(nxsea,nysea)
      save lat, lon, cvarmask

      logical firstsearise
      save firstsearise
      data firstsearise /.true./

      cvar = cvarin

      write (6,'(2a)') 
     *   'reading existing Netcdf Searise data file for ',cvar

c        Open existing Searise netcdf data file and read hi-res fields

#if defined (NOLEBROCQ)
      cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_searise/Antarctica_5km_withshelves_v4.nc'
      cvarmask = 'thkmask'
#else
      cfilin = cpathhome(1:lenchr(cpathhome)) //
     *           '/Data_searise/Antarctica_5km_dev1.0.nc'
      cvarmask = 'mask'
#endif
      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 6000)

c        Define starting/ending position eof 2-D data values in array

      start(1) = 1
      count(1) = nxsea

      start(2) = 1
      count(2) = nysea

      start(3) = 1        
      count(3) = 1

c===========================
      if (firstsearise) then
c===========================
c       latitude (deg N):
        rcode = nf_inq_varid (ncid,  'lat',  varid)
        call checkrcode (rcode, cfilin, 'lat', 7765)
        rcode = nf_get_vara_real (ncid, varid, start, count, lat)
        call checkrcode (rcode, cfilin, 'lat', 7766)

c       longitude (deg E):
        rcode = nf_inq_varid (ncid,  'lon',  varid)
        call checkrcode (rcode, cfilin, 'lon', 7767)
        rcode = nf_get_vara_real (ncid, varid, start, count, lon)
        call checkrcode (rcode, cfilin, 'lon', 7768)

        firstsearise = .false.
c==========
      endif
c==========

      if (cvar.eq.'maskwater' .or. cvar.eq.'maskpres') then
        rcode = nf_inq_varid (ncid, cvarmask, varid)
        call checkrcode (rcode, cfilin, cvarmask, 7772)
        rcode = nf_get_vara_int (ncid, varid, start, count, iarrin)
        call checkrcode (rcode, cfilin, cvarmask, 7773)
      else
        if (cvar.eq.'topgsd') then
          rcode = nf_inq_varid (ncid, 'topg', varid)
        else
          rcode = nf_inq_varid (ncid, cvar, varid)
        endif
        call checkrcode (rcode, cfilin, cvar, 7770)
        rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
        call checkrcode (rcode, cfilin, cvar, 7771)
      endif

c--

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfilin, 'close', 8000)

c--

c        Fill null values in Searise fields (topg: equatorward of 60S;
c        presprcp and presartm: all open ocean).

      if (cvar.eq.'topg' .or. cvar.eq.'topgsd') then
        do j=1,nysea
          do i=1,nxsea
            if (nint(arrin(i,j)).eq.-9999) arrin(i,j) = -5000.       ! m
          enddo
        enddo

      else if (cvar.eq.'presprcp' .or. cvar.eq.'presartm' .or.
     *         cvar.eq.'acca'     .or. cvar.eq.'accr'     .or.  !lebrocq
     *         cvar.eq.'temp') then                             !lebrocq
        do j=1,nysea
          do i=1,nxsea
            if (arrin(i,j).eq.0.) then
c             crude parameterization vs latitude (assume at sea level):
              zhs = 0.
              zts = 34.46 - .00914*zhs - .68775*abs(lat(i,j))   ! deg C
              zpr = 1.5 * (2.**(zts/10.))                       ! m/yr
              if (cvar.eq.'presprcp' .or. 
     *            cvar.eq.'acca'     .or. cvar.eq.'accr') arrin(i,j)=zpr
              if (cvar.eq.'presartm' .or. cvar.eq.'temp') arrin(i,j)=zts
            endif
          enddo
        enddo

      else if (cvar(1:8).eq.'bheatflx') then
        do j=1,nysea
          do i=1,nxsea
            if (arrin(i,j).eq.0.) then
c             as in icectl (higher geothermal flux under WAIS):
              if ( lat(i,j).gt.-86. .and.
     *             (lon(i,j).gt.170. .or. lon(i,j).lt.-30.) ) then
                arrin(i,j) = geoflux_wais / 31556926.   ! J/m2/a to W/m2
              else
                arrin(i,j) = geoflux_eais / 31556926.   ! J/m2/a to W/m2
              endif
            endif
          enddo
        enddo

      endif

c        Aggregate (tile) to model grid

c-----------------------------------------------------------
      if (cvar.ne.'maskwater' .and. cvar.ne.'maskpres') then
c-----------------------------------------------------------

c           Find [ita,itb], [jta,jtb] Searise index ranges that are 
c           within each model h-cell.
c           Searise grid has even # of h-grid cells (SP at 4 corners).
c           For model 10,20,40 km, and Searise 5 km, model centers
c           are on Searise corners, and model edges on Searise edges. 
c           For model 5 km, cells coincide if nx,ny even (and jta=jtb,
c           ita=itb), or model centers are on corners if nx,ny odd
c           (and jtb=jta+1, itb=ita+1). +/-.001 only needed for 
c           model 5 km and nx,ny odd.  

        do j=1,ny
          do i=1,nx

            jta = nint((yh(i,j)-0.5*dy0-ycorner)/ddin - .001) + 1
            jtb = nint((yh(i,j)+0.5*dy0-ycorner)/ddin + .001)
            ita = nint((xh(i,j)-0.5*dx0-xcorner)/ddin - .001) + 1
            itb = nint((xh(i,j)+0.5*dx0-xcorner)/ddin + .001)
            jtb = max (jta,jtb) ! needed for dx,dy < 5km (should bilin)
            itb = max (ita,itb) ! " " " " " " " " " " " " " " " " " " "

            arr(i,j) = 0.
            nav = 0
            do jt = jta,jtb
              jj = max (1, min (nysea, jt))
              do it = ita,itb
                ii = max (1, min (nxsea, it))
                arr(i,j) = arr(i,j) + arrin(ii,jj)
                nav = nav + 1
              enddo
            enddo
            arr(i,j) = arr(i,j) / nav

c              For std. dev. of subgrid bed topog

            if (cvar.eq.'topgsd') then
              zmean = arr(i,j)
              zms = 0.
              nav = 0
              do jt = jta,jtb
                jj = max (1, min (nysea, jt))
                do it = ita,itb
                  ii = max (1, min (nxsea, it))
                  zms = zms + (arrin(ii,jj)-zmean)**2
                  nav = nav + 1
                enddo
              enddo
              arr(i,j) = sqrt(zms/nav)
            endif

          enddo
        enddo

c          Ad-hoc adjustments to bed topography for specific basins  

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        if (cvar.eq.'topg') then
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c         Pine Island region:
c         do j=1,ny
c           do i=1,nx
c             if (alatd(i,j).gt. -77. .and. alatd(i,j).lt.-73. .and.
c    *            alond(i,j).gt.-115. .and. alond(i,j).lt.-95. ) then
c               jta = nint((yh(i,j)-0.5*dy0-ycorner)/ddin - .001) + 1
c               jtb = nint((yh(i,j)+0.5*dy0-ycorner)/ddin + .001)
c               ita = nint((xh(i,j)-0.5*dx0-xcorner)/ddin - .001) + 1
c               itb = nint((xh(i,j)+0.5*dx0-xcorner)/ddin + .001)
c               zmean = arr(i,j)
c               zms = 0.
c               nav = 0
c               do jt = jta,jtb
c                 jj = max (1, min (nysea, jt))
c                 do it = ita,itb
c                   ii = max (1, min (nxsea, it))
c                   zms = zms +  (arrin(ii,jj)-zmean)**2
c                   nav = nav + 1
c                 enddo
c               enddo
c               arr(i,j) = arr(i,j) + sqrt(zms/nav)
c             endif
c           enddo
c         enddo

c         Lambert region:
c         do j=1,ny
c           do i=1,nx
c             if (alatd(i,j).gt. -74. .and. alatd(i,j).lt.-67. .and.
c    *            alond(i,j).gt.  63. .and. alond(i,j).lt. 77. ) then
c               jta = nint((yh(i,j)-0.5*dy0-ycorner)/ddin - .001) + 1
c               jtb = nint((yh(i,j)+0.5*dy0-ycorner)/ddin + .001)
c               ita = nint((xh(i,j)-0.5*dx0-xcorner)/ddin - .001) + 1
c               itb = nint((xh(i,j)+0.5*dx0-xcorner)/ddin + .001)
c               zmean = arr(i,j)
c               zms = 0.
c               nav = 0
c               do jt = jta,jtb
c                 jj = max (1, min (nysea, jt))
c                 do it = ita,itb
c                   ii = max (1, min (nxsea, it))
c                   zms = zms + (arrin(ii,jj)-zmean)**2
c                   nav = nav + 1
c                 enddo
c               enddo
c               arr(i,j) = arr(i,j) - 2.*sqrt(zms/nav)
c             endif
c           enddo
c         enddo

c         Wilkes region:
c         do j=1,ny
c           do i=1,nx
c             if (alatd(i,j).gt. -82. .and. alatd(i,j).lt.-60. .and.
c    *            alond(i,j).gt. 130. .and. alond(i,j).lt.160. ) then
c               zx1 = (147. - alond(i,j))/12.
c               zy2 = (-73. - alatd(i,j))/8.
c               zr = sqrt (zx1**2 + zy2**2)  
c               zfac = 1. * max (0., 1.-zr) 
c               if (arr(i,j).lt.0.) then
c                 arr(i,j) = arr(i,j) - max (abs(arr(i,j))*zfac, 100.)
c               endif
c             endif
c           enddo
c         enddo
c>>>>>>>>>>>>
        endif
c>>>>>>>>>>>>

c---------
      else 
c---------

        do j=1,ny
          do i=1,nx
            jta = nint((yh(i,j)-0.5*dy0-ycorner)/ddin - .001) + 1
            jtb = nint((yh(i,j)+0.5*dy0-ycorner)/ddin + .001)
            ita = nint((xh(i,j)-0.5*dx0-xcorner)/ddin - .001) + 1
            itb = nint((xh(i,j)+0.5*dx0-xcorner)/ddin + .001)
            jtb = max (jta,jtb) ! needed for dx,dy < 5km (should bilin)
            itb = max (ita,itb) ! " " " " " " " " " " " " " " " " " " "
            zmask0 = 0.
            zmask1 = 0.
            zmask2 = 0.
            do jt = jta,jtb
              jj = max (1, min (nysea, jt))
              do it = ita,itb
                ii = max (1, min (nxsea, it))
                if (iarrin(ii,jj).eq.0) then
                  zmask0 = zmask0 + 1.
                else if (iarrin(ii,jj).eq.1) then
                  zmask1 = zmask1 + 1.
                else if (iarrin(ii,jj).eq.2) then
                  zmask2 = zmask2 + 1.
                endif
              enddo
            enddo
            zmaskwin = max (zmask0, zmask1, zmask2) 
            if (zmask0.eq.zmaskwin) then
              iarr(i,j) = 1     ! open ocn
            else if (zmask1.eq.zmaskwin) then 
              iarr(i,j) = 0     ! gd ice or land
            else if (zmask2.eq.zmaskwin) then 
              if (cvar.eq.'maskwater') then 
                iarr(i,j) = 1     ! fl ice
              else if (cvar.eq.'maskpres') then 
                iarr(i,j) = 2     ! fl ice
              endif
            endif
          enddo
        enddo

c----------
      endif
c----------

c        Convert to model units

      if (cvar.eq.'presartm' .or. cvar.eq.'temp') then
        do j=1,ny
          do i=1,nx
            arr(i,j) = arr(i,j) + 273.15                        ! C to K
          enddo
        enddo

      else if (cvar.eq.'presprcp' .or.
     *         cvar.eq.'acca'     .or. cvar.eq.'accr') then
        do j=1,ny
          do i=1,nx
c           SeaRise precip is m/y ice equiv. Convert to m/y liq. here,
c           converted back to ice equiv at end of surfbudg (iceclimall):
            arr(i,j) = arr(i,j) * (rhoice/1000.)     ! m ice/a to mesl/a
          enddo
        enddo

      else if (cvar.eq.'bheatflx' .or. 
     *         cvar.eq.'bheatflx_fox' .or. cvar.eq.'bheatflx_shapiro') 
     *     then
        do j=1,ny
          do i=1,nx
            arr(i,j) = arr(i,j)*31556926.               ! W/m2 to J/m2/a
          enddo
        enddo
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine searise_aggreg2 (arr, iarr, cvar)

c     As searise_aggreg, except for "LEBROCQ2" developmental dataset #2 
c     (deeper Recovery Basin bedrock). Has difference nxsea.
c     Field names (cvar) = usrf, lrf2, topg2 (real).

#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension arr(nx,ny), iarr(nx,ny)
      character*(*) cvar

      character*200 cfilin
      integer*4 ncid, rcode, varid, start(3), count(3)

      parameter (nxsea=1160, nysea=1120, ddin=5.e3,
     *           xcorner=-2800000., ycorner=-2800000.)
      real*4    arrin(nxsea,nysea)
      integer*4 iarrin(nxsea,nysea)

      write (6,'(2a)') 
     *   'reading existing Netcdf Searise data file for ',cvar

c        Open existing Searise netcdf data file and read hi-res fields

      if (cvar.eq.'z') then 
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *             '/Data_holschuh/NDH_ALBMAP.grid.nc'
      else
        cfilin = cpathhome(1:lenchr(cpathhome)) //
     *             '/Data_searise/ALBMAPv1.nc'
      endif
      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 6000)

c        Define starting/ending position eof 2-D data values in array

      start(1) = 1
      count(1) = nxsea

      start(2) = 1
      count(2) = nysea

      start(3) = 1        
      count(3) = 1

      rcode = nf_inq_varid (ncid, cvar, varid)
      call checkrcode (rcode, cfilin, cvar, 7770)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfilin, cvar, 7771)

c--

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfilin, 'close', 8000)

c--

c        Fill null values in Searise fields (topg2: equatorward of 60S;
c        presprcp and presartm: all open ocean).

      if (cvar.eq.'topg2') then
        do j=1,nysea
          do i=1,nxsea
            if (nint(arrin(i,j)).eq.-9999) arrin(i,j) = -5000.       ! m
          enddo
        enddo
      endif

c        Aggregate (tile) to model grid

      do j=1,ny
        do i=1,nx
          jta = nint((yh(i,j)-0.5*dy0-ycorner)/ddin - .001) + 1
          jtb = nint((yh(i,j)+0.5*dy0-ycorner)/ddin + .001)
          ita = nint((xh(i,j)-0.5*dx0-xcorner)/ddin - .001) + 1
          itb = nint((xh(i,j)+0.5*dx0-xcorner)/ddin + .001)
          arr(i,j) = 0.
          nav = 0
          do jt = jta,jtb
            jj = max (1, min (nysea, jt))
            do it = ita,itb
              ii = max (1, min (nxsea, it))
              arr(i,j) = arr(i,j) + arrin(ii,jj)
              nav = nav + 1
            enddo
          enddo
          arr(i,j) = arr(i,j) / nav
        enddo
      enddo

      return
      end

c+++++
#endif
c+++++

c-----------------------------------------------------------------------

c***********************
#if defined (SEARISEAR4)
c***********************

      subroutine searise_aggreg_ar4 (arr, cvar, time)

c     As searise_aggreg, except for Searise AR4 file.
c     Reads Searise AR4 Netcdf data file (10 km) for field name cvar
c     and time (year AP), aggregates to ice model array, returned 
c     in arr (real).
c     Field names cvar = 'instant_artm' (deg C, converted to deg K)
c                     or 'instant_prcp' (m/y)

#include <comicephys.h>
#include <comicegrid.h>
#include <comicepath.h>
#include <netcdf.inc>

      dimension arr(nx,ny)
      character*(*) cvar

      character*80 cfilin
      integer*4 ncid, rcode, varid, start(3), count(3)

      parameter (nxsea=600, nysea=600, ddin=10.e3, ntimin=96,
     *           xcorner=-2800000., ycorner=-2800000.)
      real*4 arrin(nxsea,nysea)

c     On ar4 netcdf file, time index 1 (time=0) has fields all zeros (?)
      itim = max (2, min (ntimin, int(time)+2))

      write (6,'(3a,i4)') 
     *   'reading existing Netcdf searise AR4 data file for ',cvar,
     *   '   itim=',itim

c        Open existing Searise AR4 netcdf data file and read 
c        hi-res fields

      cfilin = cpathhome(1:lenchr(cpathhome)) //
     *         'Data_searise/Ant_pddvars_2004_2098.nc'

      rcode = nf_open(cfilin, nf_nowrite, ncid)
      call checkrcode (rcode, cfilin, 'open', 6000)

c---

c     Define starting/ending position of 2-D data values in array
      start(1) = 1
      count(1) = nxsea

      start(2) = 1
      count(2) = nysea

      start(3) = itim
      count(3) = 1

c---

      rcode = nf_inq_varid (ncid,  cvar,  varid)
      call checkrcode (rcode, cfilin, cvar, 7774)
      rcode = nf_get_vara_real (ncid, varid, start, count, arrin)
      call checkrcode (rcode, cfilin, cvar, 7775)

c--

      rcode = nf_close (ncid)
      call checkrcode (rcode, cfilin, 'close', 8001)

c--

c        Aggregate (tile) to model grid.
c        Like searise_aggreg, but only works for dx0,dy0 >= 10 km. 

      do j=1,ny
        do i=1,nx
          jta = nint((yh(i,j)-0.5*dy0-ycorner)/ddin - .001) + 1
          jtb = nint((yh(i,j)+0.5*dy0-ycorner)/ddin + .001)
          ita = nint((xh(i,j)-0.5*dx0-xcorner)/ddin - .001) + 1
          itb = nint((xh(i,j)+0.5*dx0-xcorner)/ddin + .001)
          arr(i,j) = 0.
          nav = 0
          do jt = jta,jtb
            jj = max (1, min (nysea, jt))
            do it = ita,itb
              ii = max (1, min (nxsea, it))
              arr(i,j) = arr(i,j) + arrin(ii,jj)
              nav = nav + 1
            enddo
          enddo
          arr(i,j) = arr(i,j) / nav
        enddo
      enddo

c        Adjust units

      if (cvar.eq.'instant_artm') then
        do j=1,ny
          do i=1,nx
            arr(i,j) = arr(i,j) + 273.15                        ! C to K
          enddo
        enddo

      else if (cvar.eq.'instant_prcp') then
        do j=1,ny
          do i=1,nx
c           SeaRise precip is m/y ice equiv. Convert to m/y liq. here,
c           converted back to ice equiv at end of surfbudg (iceclimall):
            arr(i,j) = arr(i,j) * (rhoice/1000.)     ! m ice/a to mesl/a
          enddo
        enddo

      endif

      return
      end

c*****
#endif
c*****

c----------------------------------------------------------------------

      subroutine dummy_climparam
c     avoid compile "empty file" error
      return
      end      
