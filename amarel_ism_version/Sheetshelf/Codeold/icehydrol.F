      subroutine initwater (h, hb, hw, sealev)

c     Nb: no longer used (1/29/13). Now all done by subr findwater.

c     Sets initial water depths of ocean and inland lakes (hw).
c     No sub-ice lakes or thin film initially (hw=0).
c        
c     Generally performs each search by searching one grid box 
c     from the current set at each iteration. Maintains list(1:nlist)
c     of indices (i,j, packed) of points added in current iteration,
c     and use only those for the next iteration. 

c     Ocean and rivers can propagate N-S, E-W and along diagonals,
c     and lakes too.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny), hb(nx,ny), hw(nx,ny)

c     local:
      parameter (nwatermax = 1000)
      dimension 
     *  icefloat(nx,ny), maskh2o(nx,ny), sill(nx,ny),indsill(nwatermax),
     *  hsa(nx,ny),                      ! hs as if no floating lake ice
     *  list(nx*ny), listprev(nx*ny)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 8)   ! 4 for E-W/N-S, 8 for diagonals too

c#######################################################################
#if defined (EISLINE) || defined (EISANTA) || defined (NHA) || defined (CARB) || defined (ORDO) || defined (TEST2D) || defined (TESTMEL2D)
c#######################################################################

c        Full calculations, considering ocean access, flotation, lakes

c        Perform search for oceans (submerged points that are in 
c        contact with any submerged points at domain edges)

      call izero (icefloat, nx*ny)
      call izero (maskh2o, nx*ny)

c----------------------
      do iter = 1,nx*ny
c----------------------

        nlist = 0

        if (iter.eq.1) then

c         initial ocean set is all submerged pts around domain edges:
          do j=1,ny
            iskip = max(nx-1,1)
            if (j.eq.1.or.j.eq.ny) iskip = 1
            do i=1,nx,iskip
              if (rhoice*h(i,j).le.rholiq*(sealev-hb(i,j))) then 
                maskh2o(i,j) = 1
                icefloat(i,j) = 1
                nlist = nlist + 1
                list(nlist) = i + (j-1)*nx
              endif
            enddo
          enddo

        else


c         normal iteration: look at adjacent points for new ocean pts:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskh2o(ii,jj).eq.0) then 
                if (rhoice*h(ii,jj).le.rholiq*(sealev-hb(ii,jj))) then
                  maskh2o(ii,jj) = 1
                  icefloat(ii,jj) = 1
                  nlist = nlist + 1
                  list(nlist) = ii + (jj-1)*nx
                endif
              endif
            enddo
          enddo

        endif

c          Skip out if no new ocean points found

        if (nlist.eq.0) goto 1000

c          Copy current list to previous list

        nlistprev = nlist
        do n=1,nlist
          listprev(n) = list(n)
        enddo

c----------
      enddo
c----------
 1000 continue

c        Set initial nwater, sill, indsill just for ocean
      
      nwater = 0
      indsill(1) = 0
      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.1) then
            sill(i,j) = sealev
            if (nwater.eq.0) then
              nwater = 1                         ! 1 is ocean_embayments
              indsill(1) = i + (j-1)*nx
            endif
          else
c           for non-ocean, just in case skip out to 5000 below
            sill(i,j) = h(i,j) + hb(i,j)
          endif
        enddo
      enddo

c     skip out with just ocean, no inland lakes

      goto 5000 ! 777   

c        Perform search for lakes, by eliminating non-lakes (land pts
c        that have a downhill path to a coast)

c        Right now maskh2o=1 (ocean) or 0(not).
c        Initialize non-ocean points (maskh2o=0) to *possible* 
c        lake pts (maskh2o=-1), temporarily. Also set  hsa as 
c        if no floating lake ice.

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.1) then
            hsa(i,j) = sealev + (1.-rhoice/rholiq)*h(i,j)
          else
            maskh2o(i,j) = -1
            hsa(i,j) = hb(i,j) + h(i,j)
          endif
        enddo
      enddo

c======================
      do iter = 1,nx*ny
c======================

        nlist = 0

        if (iter.eq.1) then

c         Set initial non-lake set (maskh2o=0) to all coastal points,
c         and all land pts at edges of domain (we don't know about 
c         slopes outside):
          do j=1,ny
            do i=1,nx
              if (maskh2o(i,j).eq.-1) then
                iflag = 0
                if (i.eq.1.or.i.eq.nx.or.j.eq.1.or.j.eq.ny) then
                  iflag = 1
                else
                  do look=1,nlook
                    ii = max (1, min (nx, i + ish(look)))
                    jj = max (1, min (ny, j + jsh(look)))
                    if (maskh2o(ii,jj).eq.1) iflag = 1
                  enddo
                endif
                if (iflag.eq.1) then
                  maskh2o(i,j) = 0
                  nlist = nlist + 1
                  list(nlist) = i + (j-1)*nx
                endif
              endif
            enddo
          enddo

        else

c         normal iteration: look uphill for new non-lake points:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if ( hsa(ii,jj).ge.hsa(i,j) .and. 
     *             maskh2o(ii,jj).eq.-1 ) then 
                maskh2o(ii,jj) = 0
                nlist = nlist + 1
                list(nlist) = ii + (jj-1)*nx
              endif
            enddo
          enddo

        endif

c          Skip out if no new non-lake points found

        if (nlist.eq.0) goto 2000

c          Copy current list to previous list

        nlistprev = nlist
        do n=1,nlist
          listprev(n) = list(n)
        enddo

c==========
      enddo
c==========
 2000 continue

c       Propagate sill-depth constraints through lakes.
c       sill already set above for ocean (maskh2o=1). 

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.-1) then
            sill(i,j) = 1.e20
          else if (maskh2o(i,j).eq.0) then 
            sill(i,j) = 0. 
          endif
        enddo
      enddo

      do iter = 1,nx*ny
        nchange = 0
        do j=1,ny
          do i=1,nx
            if (maskh2o(i,j).eq.-1) then
              sillold = sill(i,j)
              do look = 1,nlook
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))
                if (maskh2o(ii,jj).eq.0) then
c                 land pt (cannot be ocn, all coastal pts are non-lake):
                  sill(i,j) = min (sill(i,j), hsa(ii,jj))
                else
c                 lake point:
                  sill(i,j) = min (sill(i,j), sill(ii,jj))
                endif
              enddo
c             remove "emergent land" from lake status:
              if (sill(i,j).le.hsa(i,j)) then
                maskh2o(i,j) = 0
                icefloat(i,j) = 0
                sill(i,j) = 0.
              endif
              if (sillold.ne.sill(i,j)) then
                nchange = nchange + 1
              endif
            endif
          enddo
        enddo
        if (nchange.eq.0) goto 3000
        if (iter.gt.nx*ny-10) then 
          write (ioterm,*) '*** Warning (initwater): iter=',iter
        endif
      enddo
      write (ioterm,*)'*** Error (initwater): exceeded iteration: ',iter
      stop
 3000 continue

c        Convert -1 maskh2o pts (currently 0=land, 1=ocean, -1=inland 
c        lakes) to labels (2,3,..) that identify contiguous inland lakes

      do j=1,ny
        do i=1,nx

          if (maskh2o(i,j).eq.-1) then
            nwater = nwater + 1
            if (nwater.gt.nwatermax) then
              write(ioterm,*)'*** Error (initwater): too many lakes'
              stop
            endif
            maskh2o(i,j) = nwater
            indsill(nwater) = 0

            nlistprev = 1
            listprev(1) = i + (j-1)*nx

            do iter = 1,nx*ny

              nlist = 0
              do n=1,nlistprev
                ja = (listprev(n)-1) / nx + 1
                ia =  listprev(n) - (ja-1)*nx
                do look = 1,nlook
                  ii = max (1, min (nx, ia + ish(look)))
                  jj = max (1, min (ny, ja + jsh(look)))
                  if ( maskh2o(ii,jj).gt.0) then
                    if (maskh2o(ii,jj).ne.maskh2o(i,j)) then
                      write (ioterm,'(a,6i4)')
     *                  '*** Error 1 (initwater) ii,jj,i,j,maskh2o=',
     *                  ii, jj, i, j, maskh2o(ii,jj), maskh2o(i,j)
                      stop
                    endif
                  else if (maskh2o(ii,jj).eq.-1) then
                    if (sill(ii,jj).ne.sill(i,j)) then
                      write (ioterm,'(a,6i4/a,2f22.15)')
     *                  '*** Error 2 (initwater) ii,jj,i,j,maskh2o=',
     *                  ii, jj, i, j, maskh2o(ii,jj), maskh2o(i,j),
     *                  '    sills=',sill(ii,jj), sill(i,j)
                      stop
                    endif
                    maskh2o(ii,jj) = maskh2o(i,j)
                    nlist = nlist + 1
                    list(nlist) = ii + (jj-1)*nx

                  else if (maskh2o(ii,jj).eq.0) then
c                   set location of first-found lake sill point (land)
                    if (hsa(ii,jj).eq.sill(i,j) .and.
     *                  indsill(nwater).eq.0) then
                      indsill(nwater) = ii + (jj-1)*nx
                    endif

                  endif
                enddo
              enddo

c                Skip out if no new lake points found

              if (nlist.eq.0) goto 4000

c                Copy current list to previous list

              nlistprev = nlist
              do n=1,nlist
                listprev(n) = list(n)
              enddo

            enddo
 4000       continue

          endif

        enddo
      enddo

c        At this point, inland lakes rest on top of either land or ice.
c        [The only floating ice is in ocean embayments, bordered by 
c        coastal land with no ice or non-floating ice.]
c        Now look to see if inland lakes float adjacent (non-ocean) ice.
c        This expands inland lakes and creates floating ice on lakes.
c        Assume that lake hydrostatic pressure penetrates horizontally
c        only along lake-land interface; for lakes with ice bottoms, 
c        cannot penetrate solid ice (ie, not into adjacent cells
c        with ice...assume continuous ice even if the adjacent ice 
c        base is higher than the current cell's ice top).

c~~~~~~~~~~~~~~~~~~~~~~
      do iter = 1,nx*ny
c~~~~~~~~~~~~~~~~~~~~~~

        nfloata = 0
        nfloatb = 0

        do j=1,ny
          do i=1,nx
            if ( maskh2o(i,j).ne.1 .and.
     *           h(i,j).gt.0. .and. icefloat(i,j).eq.0 ) then

              do look=1,nlook
                ii = max (1, min (nx, i + ish(look)))
                jj = max (1, min (ny, j + jsh(look)))

c               Two ways to float ice:

c           (a) local lake over ice, transmitted to base by adjacent
c               ice-free land (which must be higher than local sill):
                if (maskh2o(i,j).ge.2 .and.
     *              maskh2o(ii,jj).eq.0 .and. h(ii,jj).eq.0.) then
                  if (rhoice*h(i,j).le.rholiq*(sill(i,j)-hb(i,j))) then
                    icefloat(i,j) = 1
                    nfloata = nfloata + 1
                  endif

c           (b) land ice, adjacent to lake with land bottom  
                else if (maskh2o(i,j).eq.0 .and.
     *                   maskh2o(ii,jj).ge.2 .and. h(ii,jj).eq.0.)then
                  if (rhoice*h(i,j).le.rholiq*(sill(ii,jj)-hb(i,j)))then
                    sill(i,j) = sill(ii,jj)
                    maskh2o(i,j) = maskh2o(ii,jj)
                    icefloat(i,j) = 1
                    nfloatb = nfloatb + 1
                  endif
                endif
              enddo

            endif
          enddo
        enddo

c          Skip out if no new floating-ice points found

        if (nfloata + nfloatb.eq.0) goto 5000

c~~~~~~~~~~
      enddo
c~~~~~~~~~~

 5000 continue

c       Set water thickness hw   

      do j=1,ny
        do i=1,nx
          if (maskh2o(i,j).eq.0) then
            hw(i,j) = 0.
          else
            hw(i,j) = sill(i,j) - (rhoice/rholiq)*h(i,j) - hb(i,j)
            if (hw(i,j).lt.0.) then
              write (ioterm,*) '*** Error (initwater): -ve hw, i,j=',i,j
              stop
            endif
          endif
        enddo
      enddo

      return

c####
#else
c####

c       Simple setting of hw (ignoring access to ocean issues)  

      do j=1,ny
        do i=1,nx
          hw(i,j) = max (0., sealev - (rhoice/rholiq)*h(i,j) - hb(i,j))
        enddo
      enddo

      return

c#####
#endif
c#####

      end

c-----------------------------------------------------------------------

      subroutine findwater (maskwater,
     *                      indlake, npoilake, nlake,
     *                      h, hb, hw, sealev, timeice)

c        Locates oceans from scratch, resets maskwater and hw 
c        for ocean and adjacent land. Resets lakes
c        (maskwater,indlake,npoilake,nlake) 
c        depending on existing hw.

c     maskwater = 0=grounded ice or ice-free land, 
c                 1=ocean,
c                 +/- 2,3,4,...=lake number (+ open, - subice) (ret)
c                 Here, lakes all +ve, set to + or - in adjustpres
c     indlake   = (packed) i,j of each lake's points (returned)
c     npoilake  = number of points in each lake      (returned)
c     nlake     = # of lakes (including ocean)        (returned)
c     h         = ice thickness                       (supplied)
c     hb        = bed+sed elevation                   (supplied)
c     hw        = water thickness           (modified for ocean)
c     sealev    = sea level                           (supplied)
c        
#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  maskwater(nx,ny), indlake(npoimax,nlakemax), 
     *  npoilake(nlakemax),
     *  h(nx,ny), hb(nx,ny), hw(nx,ny)

c     local:
      dimension
     *  list(nx*ny),      listprev(nx*ny)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 8)   ! 4 for E-W/N-S, 8 for diagonals too

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if ( defined (EISLINE) && (defined (LINEB) || defined (LINEE)) ) || ( defined (EISANTA) && defined (TRANSECTA) )
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c     simple setting of hw (flotation by ocean at sealev everywhere,
c     ignoring access issues):
      call izero (maskwater, nx*ny)
      call zero (hw, nx*ny)
      do j=1,ny
        do i=nx,1,-1
          if (rhoice*h(i,j).le.rholiq*(sealev-hb(i,j))) then
            maskwater(i,j) = 1
            hw(i,j) = sealev - hb(i,j) - (rhoice/rholiq)*h(i,j)
          endif
        enddo
      enddo
      nlake = 1 ! ocean
      call izero (indlake,    npoimax*nlakemax)
      call izero (npoilake,   nlakemax)
      go to 5000 ! return

c~~~~~~~~~~~~~~~~~~~~~~
#elif defined (EISLINE) 
c~~~~~~~~~~~~~~~~~~~~~~

c     simple setting of ocean (scanning r-to-l to first grounded pt):
      call izero (maskwater, nx*ny)
      call zero (hw, nx*ny)
      if (nx.gt.ny) then
        do j=1,ny
          do i=nx,1,-1
            if (rhoice*h(i,j).le.rholiq*(sealev-hb(i,j))) then
              maskwater(i,j) = 1
              hw(i,j) = sealev - hb(i,j) - (rhoice/rholiq)*h(i,j)
            else
              goto 50
            endif
          enddo
   50     continue
        enddo
      else
        do i=1,nx
          do j=ny,1,-1
            if (rhoice*h(i,j).le.rholiq*(sealev-hb(i,j))) then
              maskwater(i,j) = 1
              hw(i,j) = sealev - hb(i,j) - (rhoice/rholiq)*h(i,j)
            else
              goto 52
            endif
          enddo
   52     continue
        enddo
      endif
      nlake = 1 ! ocean
      call izero (indlake,    npoimax*nlakemax)
      call izero (npoilake,   nlakemax)
      go to 5000 ! return

c~~~~~
#endif
c~~~~~

      call izero (maskwater, nx*ny)

c        Perform search for oceans (submerged points that are in
c        contact with any submerged points at domain edges),
c        setting maskwater to 1 and setting hw for ocean pts,
c        and also setting hw to 0 for non-floating ocean-adjacent pts.

c----------------------
      do iter = 1,nx*ny
c----------------------

        nlist = 0

        if (iter.eq.1) then

          do j=1,ny
            iskip = nx-1
            if (j.eq.1.or.j.eq.ny) iskip = 1
            do i=1,nx,iskip
              if (sealev-hb(i,j)-rhor*h(i,j).gt.hwcut) then
                maskwater(i,j) = 1
                nlist = nlist + 1
                list(nlist) = i + (j-1)*nx
                hw(i,j) = sealev - hb(i,j) - rhor*h(i,j)
              else
                hw(i,j) = 0.
              endif
            enddo
          enddo

        else

c         normal iteration, look at adjacent points for new ocean pts:
          do n=1,nlistprev
            j = (listprev(n)-1) / nx + 1
            i =  listprev(n) - (j-1)*nx
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (maskwater(ii,jj).eq.0) then
                if (sealev-hb(ii,jj)-rhor*h(ii,jj).gt.hwcut) then
                  maskwater(ii,jj) = 1
                  nlist = nlist + 1
                  list(nlist) = ii + (jj-1)*nx
                  hw(ii,jj) = sealev - hb(ii,jj) - rhor*h(ii,jj)
                else 
                  hw(ii,jj) = 0.
                endif
              endif
            enddo
          enddo

        endif

c          Skip out if no new ocean points found

        if (nlist.eq.0) goto 1000

c          Copy current list to previous list

        nlistprev = nlist
        do n=1,nlist
          listprev(n) = list(n)
        enddo

c----------
      enddo
c----------
 1000 continue

c        Search for inland lakes, based on water thickness hw vs hwcut.
c        When find a new, not-already-found lake point (outer loop),
c        increment nlake, then iteratively search for adjacent lake 
c        points, forming list indlake and setting maskwater, 
c        and set npoilake when done.

      nlake = 1 ! ocean
      call izero (indlake,  npoimax*nlakemax)
      call izero (npoilake, nlakemax)

#if ! defined (WATERLAND) 
c     don't allow any non-ocean water (open lakes or sub-ice)
      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).ne.1) then
c            convert liquid water beneath ice to ice (or not..why?)
c            if (h(i,j).gt.hwcut) h(i,j) = h(i,j) + hw(i,j)/rhor
             hw(i,j) = 0.
          endif
        enddo
      enddo
      go to 5000 ! return 
#endif

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do ja=1,ny
        do ia=1,nx
          if (maskwater(ia,ja).eq.0 .and. hw(ia,ja).gt.hwcut) then    
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            nlake = nlake + 1
            if (nlake.gt.nlakemax) then
               write (ioterm, *) 
     *           'Error (findwater): max # lakes exceeded, nlake=',nlake
               stop
             endif

c----------------------------
            do iter = 1,nx*ny
c----------------------------

              nlist = 0

              if (iter.eq.1) then

                nlist = 1
                list(1) = ia + (ja-1)*nx 
                maskwater(ia,ja) = nlake
                npoi = 1
                indlake(npoi,nlake) = list(1)

              else

c               normal iteration, look at adjacent pts for new lake pts
                do n=1,nlistprev
                  j = (listprev(n)-1) / nx + 1
                  i =  listprev(n) - (j-1)*nx
                  do look=1,nlook
                    ii = max (1, min (nx, i + ish(look)))
                    jj = max (1, min (ny, j + jsh(look)))
                    if (maskwater(ii,jj).eq.0 .and. hw(ii,jj).gt.hwcut)
     *                then    
                      maskwater(ii,jj) = nlake
                      nlist = nlist + 1
                      list(nlist) = ii + (jj-1)*nx
                      npoi = npoi + 1
                      if (npoi.gt.npoimax) then
                        write (ioterm,'(/a/a,i6,a,i6)') 
     *                    'Error (findwater): max lake size exceeded:',
     *                    '   npoi=',npoi,'  nlake=',nlake
                        stop
                      endif
                      indlake(npoi,nlake) = list(nlist)
                    else if (maskwater(ii,jj).ne.nlake) then
                      write (ioterm,'(/a/a,i6/a,i6/a,2i6/a,2i6)') 
     *                  'Error (findwater): lake iteration error:',
     *                  '    maskwater=',maskwater(ii,jj),
     *                  '    nlake    =',nlake,
     *                  '    ii,jj    =',ii,jj,
     *                  '    ia,ja    =',ia,ja
                      stop
                    endif
                  enddo
                enddo

              endif

c                Skip out if no new open-lake points found

              if (nlist.eq.0) goto 2000

c                Copy current list to previous list

              nlistprev = nlist
              do n=1,nlist
                listprev(n) = list(n)
              enddo

c----------------
            enddo
c----------------
 2000       continue
            npoilake(nlake) = npoi

c>>>>>>>>>>>>>>
          endif
        enddo
      enddo
c>>>>>>>>>>>>>>

c        Finally, set all maskwater's within each lake to negative if
c        lake is closed (else maskwater remains positive, open lake)

      do ilake = 2,nlake

        ifopen = 0
        do n=1,npoilake(ilake)
          j = (indlake(n,ilake)-1) / nx + 1
          i =  indlake(n,ilake) - (j-1)*nx
          if (h(i,j).eq.0.) ifopen = 1
        enddo
  
        if (ifopen.eq.0) then
          do n=1,npoilake(ilake)
            j = (indlake(n,ilake)-1) / nx + 1
            i =  indlake(n,ilake) - (j-1)*nx
            maskwater(i,j) = -maskwater(i,j)
          enddo
        endif

      enddo

 5000 continue
      call scopy_i (nx*ny, maskwater, 1, maskcur, 1)      ! for printmap

      return
      end

c-----------------------------------------------------------------------

      subroutine adjustpres (maskwater, indlake, npoilake, nlake,
     *                       h, hb, hw, hs, sealev)

c     Adjusts water thicknesses hw within each separate lake or ocean
c     to equalize hydrostatic water pressure on a commom reference
c     level. Then resets surface elevation hs. 
c     Ignore consequent changes to maskwater, indlake, etc (i.e., if 
c     change in hw would make ice ground or float, or if hw would 
c     overflow sill of an open lake...wait till next call to findwater 
c     and/or movewater).
c
c     maskwater = 0,1,+/-2... (grounded ice/ice-free land, ocean, lake#)
c     indlake   = (packed) i,j of each lake's points
c     npoilake  = number of points in each lake
c     nlake     = # of lakes (including ocean)
c     h         = ice thickness
c     hb        = bed+sed elevation
c     hw        = water thickness (modified)
c     hs        = surface elevation (modified)
c     sealev    = sea level

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  maskwater(nx,ny), indlake(npoimax,nlakemax), npoilake(nlakemax),
     *  h(nx,ny), hb(nx,ny), hw(nx,ny), hs(nx,ny)

      dimension pres(npoimax)

c        Set hw and hs for land (no lake) and ocean

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.0) then
            hs(i,j) = hb(i,j) + hw(i,j) + h(i,j)
          else if (maskwater(i,j).eq.1) then
            hw(i,j) = max (0., sealev - hb(i,j) - rhor*h(i,j))
            hs(i,j) = hb(i,j) + hw(i,j) + h(i,j)
          endif
        enddo
      enddo

c        Loop for each lake (open or sub-ice, but not ocean)

c-----------------------
      do ilake = 2,nlake
c-----------------------

c          Calculate pres = pressure/(rholiq*grav) at "zero" elevation 
c          for each lake point, and its areal average totp

        totp = 0.
        tota = 0.
        do n=1,npoilake(ilake)
          j = (indlake(n,ilake)-1) / nx + 1
          i =  indlake(n,ilake) - (j-1)*nx
          pres(n) = rhor*h(i,j) + (hb(i,j)+hw(i,j))
          totp = totp + pres(n)*darea(i,j)
          tota = tota +         darea(i,j)
        enddo
        totp = totp/tota
  
c          Subtract pres-totp from each point's water thickness,
c          to equalize new pres's for all points. This conserves
c          total water in lake, unless...if new hw would be negative 
c          (newly grounded ice), just set hw to zero, ignore water
c          non-conservation, and wait till next call to findwater 
c          to correct lake indices.
c          This workd both for closed lakes (ice rides above
c          hw, and open lakes (ice floats in hw).

        do n=1,npoilake(ilake)
          j = (indlake(n,ilake)-1) / nx + 1
          i =  indlake(n,ilake) - (j-1)*nx
          hw(i,j) = max (0., hw(i,j) - (pres(n)-totp))
          hs(i,j) = hb(i,j) + hw(i,j) + h(i,j)
        enddo

c----------
      enddo
c----------

      return
      end

c-----------------------------------------------------------------------

      subroutine basecoef_init (sedpres, maskinit, crhmel, timeice)

c     Sets basal flow coeff crhmel (h-grid), at start of run only, 
c     ignoring melting and grounded criteria. crhmel wll be used to set 
c     crh[h,u,v] (u,v grids) in basecoef, allowing for those criteria.
c     If ifrest= 0 or 1,  crhmel will be overridden from restart file,
c     and/or from crhmelfilein file if -DCRHMELREAD.
c     If CRHINV or ENHINV, crhmel will be adjusted during runtime 
c     in basecoef_inv.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  sedpres(nx,ny),   maskinit(nx,ny), crhmel(nx,ny)

#if defined (SHALDRIL) && defined (DWSEDPOLY)
c     for D. Wilson's E-O Weddell Sea sediment area (polygon outline)

      dimension dwsed (nx,ny)
      save dwsed
      parameter (ndwsed = 44)
      dimension xdwsedin(2,ndwsed), xdwsed(ndwsed), ydwsed(ndwsed) 
      data xdwsedin / 
     *  -776.250,  1498.750,
     *  -773.750,  1376.250,
     *  -773.750,  1263.750,
     *  -798.750,  1183.750,
     *  -786.250,  1076.250,
     *  -741.250,  978.750,
     *  -698.750,  913.750,
     *  -673.750,  858.750,
     *  -663.750,  778.750,
     *  -678.750,  726.250,
     *  -681.250,  643.750,
     *  -681.250,  551.250,
     *  -686.250,  473.750,
     *  -701.250,  311.250,
     *  -708.750,  216.250,
     *  -758.750,  131.250,
     *  -811.250,  108.750,
     *  -881.250,  121.250,
     *  -946.250,  136.250,
     *  -993.750,  168.750,
     *  -1031.250,  226.250,
     *  -1101.250,  271.250,
     *  -1168.750,  303.750,
     *  -1206.250,  343.750,
     *  -1251.250,  381.250,
     *  -1296.250,  426.250,
     *  -1336.250,  478.750,
     *  -1368.750,  538.750,
     *  -1383.750,  603.750,
     *  -1396.250,  668.750,
     *  -1406.250,  731.250,
     *  -1436.250,  798.750,
     *  -1491.250,  851.250,
     *  -1536.250,  888.750,
     *  -1576.250,  921.250,
     *  -1611.250,  948.750,
     *  -1691.250,  1011.250,
     *  -1818.750,  1066.250,
     *  -1952.500,  1092.500,
     *  -2087.500,  1127.500,
     *  -2177.500,  1207.500,
     *  -2202.500,  1317.500,
     *  -2212.500,  1442.500,
     *  -2197.500,  1537.500 /
#endif

#if defined (SHALDRIL) && defined (DWSEDPOLY)
c     set D. Wilson's E-O Weddell Sea sediment area (polygon outline)
      do m=1,ndwsed
        xdwsed(m) = xdwsedin(1,m) * 1000.    ! km to m
        ydwsed(m) = xdwsedin(2,m) * 1000.    ! km to m
      enddo
      call zero (dwsed, nx*ny)
      call polygon_wais (dwsed, xh,yh,nx,ny, xdwsed,ydwsed,ndwsed,1.)

      call printmap (timeice, dwsed,  'dwsed',   0.1, 134, 0)
      call flush (134)
c     stop
#endif

c for powb=2: zcrh=1.e-15 ~no sliding, 1.e-10 ~some, 1.e-5 ~stream.
c             See crhnos, crhhard, crhsed in comicephys.h.

c        Set basal sliding coefficient crhmel for various experiments

      call resetr (crhmel, nx*ny, crhnos)

c================
      do j=1,ny
        do i=1,nx
c================
c....................
#if defined (EISLINE)
c....................

#  if defined (LINEB) 
          zr = xh(i,j)
          if (zr.lt.300.e3) then
            zcrh = 1.e-8
          else
            zcrh = 1.e-4
          endif
#  elif defined (LINEC)
          zcrh = crhnos
#  elif defined (LINED)
          zcrh = 1.e-10
#  elif defined (LINEE)
          zcrh = 2.e-5    ! abc
#  elif defined (LINEF)
c         zcrh = 1.e-8 

c step function between 2 values:
          zwei = mod (int(xh(i,j)/100.e3), 2)
          zcrh = 1.e-10*(1.-zwei) + 1.e-6*zwei
c uniform:
c         zcrh = 1.e-8

cc        zcrh = 3.e-9  ! std
cc        zcrh = 1.e-9  ! Fa
cc        zcrh = 1.e-8  ! Fb (and Fb_oldaslip wit budgsnow=0.1)
cc        zcrh = 3.e-8  ! Fc
cc        zcrh = (86400.*365.) / (7.624e6**powb)        ! powb=3?,schoof
#  elif defined (LINEG)
          zcrh = 1.e-10
#  elif defined (LINEH) || defined (LINEH2) || defined (LINEH3) || defined (LINEH4)
          zcrh = 1.e-10  ! shama
#  elif defined (LINEM)
#    if defined (BASALa)
          zcrh = (86400.*365.) / (1.e7**powb)           ! powb=3,mismip
#    elif defined (BASALb)
          zcrh = (86400.*365.) / (7.2082e10**powb)      ! powb=1,mismip
#    endif
#  elif defined (SLMODEL)
c         zcrh = 3.e-9
          zcrh = 1.e-10
#  endif

c......................
#elif defined (EISANTA)
c......................

          if ( sedpres(i,j).gt.0.
#    if defined (SHALDRIL) && defined (DWSEDPOLY)
     *         .or. dwsed(i,j).gt.0.
#    endif
     *       ) then
c           slippery bed where modern rebounded topog is below sl:
c           zcrh = crhsed
            zcrh = 1.e-5
          else
c           hard bed where modern rebounded topog is above sl:
            zcrh = crhhard
          endif

c    for Runwant..., Runwil...:
#    if defined (BEDALL) 
          zcrh = crhhard
#    elif defined (SEDALL)
c         zcrh = crhsed
          zcrh = 1.e-5
#    endif

c    for shama:
c         zcrh = crhhard

c..................
#elif defined (NHA)
c..................

#  if defined (GLACBAY)
          zcrh = 1.e-10  ! e,f
          zcrh = 1.e-7   ! g
#  else
c post MPT?:
          zcrh = 1.e-10
          if (sedpres(i,j).gt.0.) zcrh = 1.e-6   ! 1.e-5 ?
c pre MPT?: 
cc        zcrh = 1.e-10 ! cyci
c         zcrh = 1.e-8  ! cyca-h
cc        zcrh = 1.e-7  ! cycj
#  endif

c.....................................
#elif defined (CARB) || defined (ORDO)
c.....................................

          zcrh = 1.e-10
          if (sedpres(i,j).gt.0.) zcrh = 1.e-5

c.....................
#elif defined (TEST2D)
c.....................

c         for mismipplus:

          zcrh = 1.e-12   ! (1.e-4)^3   m a-1 Pa-3

c........................
#elif defined (TESTMEL2D)
c........................

          zcrh = 1.e-10

c....
#else
c....

          zcrh = 1.e-10

c.....
#endif
c.....

          crhmel(i,j) = zcrh

c============
        enddo
      enddo
c============

      return
      end

c-----------------------------------------------------------------------

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#if defined (CRHINV) || defined (ENHINV)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine basecoef_inv (h, hs, hsbal, hbsd, hw, maskwater, t,
     *                         baseperc, basefrml,
     *                         crhmel, enhanceinv, 
     *                         u, v, w, 
     *                         timeice, dtimeice)

c      If CRHINV:
c      "Inverse" calculations to dynamically adjust basal sliding
c      coefficients (crhmel) (on h-grid during a forward run,
c      aiming for hs = hsbal). Only do for grounded ice points. 
c      Limit max crhmel to <= crhsed.
c      If NOBASET, do for all grounded ice pts, ignoring fract
c      (from calc_crh) and allow lower min limit crhmel >= crhnos.
c      If not NOBASET, don't do for pts where bed temps are frozen
c      preventing sliding (fract,calc_crh), and limit 
c      min crhmel >= crhhard.
c
c      If ENHINV:
c      Also adjust additional enhancement factor for sheet flow
c      (enhanceinv), applied in subr arrhenius (icetherm.F), 
c      aiming for same fractional area of frozen bed as area
c      with crhmel=crhmin. enhanceinv can depend spatially on
c      distance to closest ice dome (enha, enhb, enhx).      
c
c      Writes crhmelfilein file, to be read by later non-inverse run.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny),         hs(nx,ny),     
     *  hsbal(nx,ny),     hbsd(nx,ny), 
     *  hw(nx,ny),        maskwater(nx,ny),  t(nx,ny,0:nlevp),
     *  baseperc(nx,ny),  basefrml(nx,ny),
     *  crhmel(nx,ny),    enhanceinv(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  w(nx,ny,0:nlevp)

#if defined (CRHINV)
      dimension crhmelprev(nx,ny), crhmelnew(nx,ny)
      logical firstprev
      data firstprev /.true./
      save crhmelprev, crhmelnew, timeprev, firstprev
#endif

#if defined (ENHINV)
      dimension distdome(nx,ny)
      save enha, enhb, enhx
      data enha, enhb, enhx /1., 1., 1000.e3/ 
c     data enha, enhb, enhx /1., .3,  500.e3/ 
#endif

c===

      crhmax = crhsed   ! sediment (crhsed in comicephys)
#if defined (NOBASET)
      crhmin = crhnos   ! frozen   (crhnos in comicephys)
#else
      crhmin = crhhard  ! hard bed (crhhard in comicephys)
#endif

c    set interval between adjustments below (yrs):
#if defined (DTIMECRHINV)
      dtimecrhinv = DTIMECRHINV
#else
#  if defined (NESTING)
      dtimecrhinv = 300.
#  else
      dtimecrhinv = 2000.
#  endif
#endif

c     initialize quantities to do gradual ramp adjustment to crhmel:
      if (firstprev) then
        call scopy (nx*ny, crhmel, 1, crhmelprev,1) 
        call scopy (nx*ny, crhmel, 1, crhmelnew, 1) 
        timeprev = timeice
        firstprev = .false.
      endif

c       Only adjust at intervals of dtimecrhinv yrs, to allow 
c       ice thicknesses and temps to adjust in between.
c       If inbetween (and CRHINV),gradually ramp crhmel from crhmelprev
c       to crhmelnew, over interval taucrh (to avoid blowups)
   
      dt = dtimeice
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if ( mod(abs(timeice)+0.5*dt,max(dtimecrhinv,dt)).ge.dt) then
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if defined (CRHINV) 

c       set time interval over which to ramp crhmel:
#  if defined (TAUSCALE)
        taucrh = max (.001, min (dtimecrhinv*0.5, TAUSCALE))
#  else
        taucrh = .001         ! no ramp, change immediately to crhmelnew
#  endif

        zwt = max (0., min (1., (timeice-timeprev) / taucrh))

c       ramp crhmel from crhmelprev to crhmelnew. For points not
c       adjusted last time (ocean, no ice, frozen,...), 
c       crhmelnew=crhmelprev=previous crhmel, so no change:
#  if defined (NESTING)
        do j=2,ny-1
          do i=2,nx-1
#  else
        do j=1,ny
          do i=1,nx
#  endif
            zp = (1. - zwt) * log10(crhmelprev(i,j))
     *               + zwt  * log10(crhmelnew(i,j))
            crhmel(i,j) = max (crhmin, min (crhmax, 10.**zp))
          enddo
        enddo

#endif

        return
c++++++++++
      endif
c++++++++++

c>>>>>>>>>>>>>>>>>>>
#if defined (CRHINV)
c>>>>>>>>>>>>>>>>>>>

c     set scaling of hs error (m) to set multiplicative factor x crhmel:
#if defined (HSCALE) 
      hscale = HSCALE
#else
#  if defined (NESTING)
      hscale = 1000.
#  else
      hscale = 2000.
#  endif
#endif

c        Save next crhmelprev and timeprev, and also crhmelnew
c        (crhmelnew is changed below for adjusted cells) 

      call scopy (nx*ny, crhmel, 1, crhmelprev, 1) 
      call scopy (nx*ny, crhmel, 1, crhmelnew, 1) 
      timeprev = timeice

c        Adjustments based on ice elevations (hs) compared to modern
c        observed (hsbal). Do only for grounded points with ice.
c        Use ad-hoc notional relation between [hs minus hsbal]
c        and multiplicative factor to adjust crhmel by.

#if defined (NESTING)
c     for nesting, edges are set in getnest, don't try to invert there
      do j=2,ny-1
        do i=2,nx-1
#else
      do j=1,ny
        do i=1,nx
#endif
          if (h(i,j).gt.0. .and. maskwater(i,j).eq.0) then

c           Only adjust where bed sliding is not completely limited by 
c           temperature, i.e., fract > 0. (nb: If NOBASET, calc_crh
c           sets fract = 1)
            call calc_crh (zcrh, fract, fracb,
     *                     crhmel(i,j),
     *                     t(i,j,nlevp) - (tmelt-dtmdh*h(i,j)),
     *                     baseperc(i,j)-basefrml(i,j),
     *                     h(i,j),
     *                     hbsd(i,j),
     *                     0)
            if ( (fract.gt.0.) .and.
     *           (
     *           (hs(i,j).gt.hsbal(i,j) .and. crhmel(i,j).lt.crhmax).or.
     *           (hs(i,j).lt.hsbal(i,j) .and. crhmel(i,j).gt.crhmin) 
     *           )
     *         ) then

              zfc = 10.** max ( -1.5, min ( 1.5, 
     *                          (hs(i,j)-hsbal(i,j)) / hscale ) )
              crhmelnew(i,j) = max ( crhmin, min ( crhmax, 
     *                               crhmel(i,j)*zfc ) )
            endif

          endif
        enddo
      enddo

c>>>>>
#endif
c>>>>>

c===

c<<<<<<<<<<<<<<<<<<<
#if defined (ENHINV)
c<<<<<<<<<<<<<<<<<<<

c        Adjust additional sheet-flow enhancement coefficient
c        (enhanceinv). First, calculate continental fractional 
c        (grounded-ice) areas where crhmel=crhmin (fac) and 
c        basal temp=frozen (fat)

      call basal_t_c (fac, fat, crhmel, t, hbsd, h, maskwater, crhmin)

c        Calculate distance to nearest dome (m), representing amount
c        of fabric hardening (distdome, the independent variable
c        in parameterization of enhanceinv)

      call path2dome (distdome, h, maskwater, u, v, w, timeice) 
c     call path2domesimp (distdome, h, hs, maskwater)

c        Adjust parameters in enhaniceinv parameterization 
c        (in do loop below), aiming to make fac = fat

      enhaold = enha
      enha = max ( 0.1, min ( 10., enha * 2.**(fat-fac) ))
      enhb = enhb*(enha/enhaold)      ! change enha, enhb by same factor
                                      ! (keeps ratio enha:enhb constant)

c        Set new enhanceinv (for non-ocean areas)

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.0) then
c           zw = max (0., min (1., 1. - distdome(i,j)/enhx))
            zw =  exp (-distdome(i,j)/enhx)
            enhanceinv(i,j) = zw*enha + (1.-zw)*enhb
          else
            enhanceinv(i,j) = enhb
          endif
        enddo
      enddo

      write (6,'(/a,i10,4f8.3/)') 
     *  'basecoef_inv: timeice, fac, fat, enha, enhb=',
     *  nint(timeice), fac, fat, enha, enhb

c<<<<<
#endif
c<<<<<

c===

c        Write crhmel and enhanceinv to separate file (crhmelfile),
c        which can be used to initialize crhmel,enhanceinv for   
c        subsequent run (crhmelfilein, read by readcrhmel if CRHMELREAD)

      call writecrhmel (crhmel, enhanceinv, timeice+dtimeice)

      return
      end

c%%%%%
#endif
c%%%%%

c-----------------------------------------------------------------------

c(((((((((((((((((((
#if defined (ENHINV)
c(((((((((((((((((((

      subroutine basal_t_c (fac, fat, 
     *                      crhmel, t, hbsd, h, maskwater, crhmin)

c     Set fractional continental grounded-ice area with
c     basal sliding coefficient crhmel = crhmin (fac),
c     and with basal temperature below pressure melt point - 1 (fat)

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  crhmel(nx,ny), t(nx,ny,0:nlevp), hbsd(nx,ny),   
     *  h(nx,ny), maskwater(nx,ny)

      fac = 0.
      fat = 0.
      atot = 0.
      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0. .and. maskwater(i,j).eq.0) then
            if (crhmel(i,j).eq.crhmin) then
              fac = fac + darea(i,j)
            endif
            if (t(i,j,nlevp).lt.tmelt-dtmdh*h(i,j) - 1.) then
              fat = fat + darea(i,j)
            endif
            atot = atot + darea(i,j)
          endif
        enddo
      enddo

      fac = fac/atot
      fat = fat/atot

      return
      end

c-----------------------------------------------------------------------

      subroutine path2dome (statdist, h, maskwater, u, v, w, timeice)

c        Traces ice cells, using Lagrangian stepping and 3-D (here 2-D)
c        linear interp of u,v (not w). Do only for starting cells in 
c        one layer (k=levpath), and arbitrarily set w=0 to keep
c        paths at that level. So traces paths back to where (u,v)
c        are zero and diverging (a dome). Returns distance from each 
c        cell to final path location (statdist), i.e. distance to 
c        nearest dome (statdist), used for setting enhanceinv in 
c        basecoef_inv. statdist results are ~same for any choice of
c        levpath (1,5,...).
c
c        Lagrangian tracking iteratively steps each point backwards
c        in time, with an adjustable timestep (dtlag) depending
c        on local velocities (so v*dt ~ 0.05*grid size).
c
c        Vertical coord is ice sheet's zeta (0 to 1, sfc to base),
c        and vert veloc is w (=d(zeta)/dt).

c        Uses exactly same Lagrangian tracking as in lagrango18 (icectl)
c        if LAGRANGEO18, but here restricted to one level (k=levpath).

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  statdist(nx,ny), 
     *  h(nx,ny),               maskwater(nx,ny),
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  w(nx,ny,0:nlevp)

      dimension
     *  xsurf(nx,ny),   ysurf(nx,ny),
     *  statifm(nx,ny),    
     *  statu(nx,ny),   statv(nx,ny),
     *  stattim(nx,ny), statnum(nx,ny), 
     *  worka(nx,ny),   workb(nx,ny)

      parameter (nzfine = 100000)
      dimension izfinea(nzfine), izfineb(nzfine), wzfinea(nzfine)
      save izfinea, izfineb, wzfinea

      parameter (niterlag = 2000000)    ! max number of iterations

      parameter (levpath = 1)           ! restrict to one level
c     parameter (levpath = 5)           ! restrict to one level

      parameter (stattimmax = 200000.)  ! max time on path

      logical first, firstprint
      data first, firstprint /.true., .true./
      save first, firstprint

c        If first call, set indices from very fine uniform vertical
c        grid (0 to 1) to model vertical grid (saved). Only use
c        vertical levels between k=1 and nlev (not 0 and nlevp).

c--------------------
      if (first) then
c--------------------
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
            write (6,'(/a)') '*** Error: path2dome (init)'
            stop
          endif
  100     continue
        enddo

        first = .false.
c----------
      endif
c----------

      call resetr (statdist,nx*ny, 1.e20)
      call zero (xsurf,   nx*ny) 
      call zero (ysurf,   nx*ny) 
      call zero (statu,   nx*ny)
      call zero (statv,   nx*ny)
      call zero (statifm, nx*ny)
      call zero (stattim, nx*ny)
      call zero (statnum, nx*ny)

c        Main loop over each starting cell with ice.
c        Only do for grounded ice (not sure if backward Lagrangian
c        tracking will work with fast shelf velocities). 

c>>>>>>>>>>>>>>
      do j=1,ny
        do i=1,nx
c         if (h(i,j).gt.0.) then                               ! all ice
          if (h(i,j).gt.0. .and. maskwater(i,j).eq.0) then    ! grounded
            k=levpath                          ! only one starting layer
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
c    *            'Warning (path2dome h=0): iterlag=', iterlag,
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
c             wp =       wwi *    wwj *    wk  * w(iwa,jwa,ka)
c    *             + (1.-wwi)*    wwj *    wk  * w(iwb,jwa,ka)
c    *             +     wwi *(1.-wwj)*    wk  * w(iwa,jwb,ka)
c    *             + (1.-wwi)*(1.-wwj)*    wk  * w(iwb,jwb,ka)
c
c    *             +     wwi *    wwj *(1.-wk) * w(iwa,jwa,kb)
c    *             + (1.-wwi)*    wwj *(1.-wk) * w(iwb,jwa,kb)
c    *             +     wwi *(1.-wwj)*(1.-wk) * w(iwa,jwb,kb)
c    *             + (1.-wwi)*(1.-wwj)*(1.-wk) * w(iwb,jwb,kb)
c             set vertical veloc = 0 to restrict to one level (levpath)
              wp = 0.
c===

c                Set current timestep (dtlag)

              zdtx = 0.05 * dx0       / max (1.e-4, abs(up))
              zdty = 0.05 * dy0       / max (1.e-4, abs(vp))
              dtlag = max (0.1, min (zdtx, zdty, 1000.))
c             zdtz = 0.05 * (1./nlev) / max (1.e-6, abs(wp))
c             dtlag = max (0.1, min (zdtx, zdty, zdtz, 1000.))
c             write (198,'(4i8,f15.5)') iterlag, i,j,k,dtlag

c               Do backward Lagrangian step (-veloc*dt)

              xp = xp - dtlag*up
              yp = yp - dtlag*vp
              zp = zp - dtlag*wp

              stattim(i,j) = stattim(i,j) + dtlag
              statnum(i,j) = statnum(i,j) + 1.

c                If gone out of bottom of ice sheet (kind of an error,
c                or +ve basefrml?), stop tracking, skip out (to 200).
c                Nb: Won't happen here, restriced to one level (levpath)

              if (zp.gt.1.) then
c               write (6,"(/a,i6,'  xp,yp,zp=',2f10.2,f10.5)")
c    *            'Warning (path2dome bot, -1): iterlag=', iterlag,
c    *            xp*.001, yp*.001, zp
                xsurf(i,j) = xp
                ysurf(i,j) = yp
                statu(i,j) = up
                statv(i,j) = vp
                statifm(i,j) = -1.
                go to 200
              endif

c                If gone out of "sides" of ice sheet (kind of an error),
c                stop tracking, skip out (to 200)

              ip = nint ((xp - xoffa + 0.5*nx*dx0)/dx0 + 0.5)
              jp = nint ((yp - yoffa + 0.5*ny*dy0)/dy0 + 0.5)
              ip = max (1, min(nx,ip))
              jp = max (1, min(ny,jp))
              if (h(ip,jp).eq.0.) then
                xsurf(i,j) = xp
                ysurf(i,j) = yp
                statu(i,j) = up
                statv(i,j) = vp
                statifm(i,j) = -2.
c               write (191,"(a/a,3i6,i8/a,2f9.0,f8.4,3f10.5/a,4f10.3)")
c    *          'Warning (path2dome side, -2):',
c    *          '  i,j,k,iterlag=', i,j,k, iterlag-1,
c    *          '  xp,yp,zp,up,vp,wp=', xp*.001,yp*.001,zp, up,vp,wp,
c    *          '  h[1-4]=', h(iwa,jwa),h(iwb,jwa),h(iwb,jwa),h(iwb,jwb)
                go to 200
              endif

c                If gone out of top surface of ice sheet, interp back
c                within timestep for more exact surface locn, skip out.
c                Nb: Won't happen here, restriced to one level (levpath)

              if (zp.le.0.) then
                wei = (-zp) / (-zp + zprev)  
                xsurf(i,j) = (1.-wei)*xp + wei*xprev 
                ysurf(i,j) = (1.-wei)*yp + wei*yprev 
                statu(i,j) = up
                statv(i,j) = vp
                statifm(i,j) = 1.
c               write(190,"(a/a,3i6,i8 /a,2i8,f8.4, 3e14.5 /a,4f10.3)")
c    *          'path2dome ok, 1:',
c    *          '  i,j,k,iterlag=', i,j,k, iterlag,
c    *          '  xs,ys,zp,up,vp,wp=', 
c    *             nint(xsurf(i,j)*.001),nint(ysurf(i,j)*.001),zp,
c    *             up,vp,wp,
c    *          '  h[1-4]=', h(iwa,jwa),h(iwb,jwa),h(iwa,jwb),h(iwb,jwb)
                goto 200
              endif

c               Limit total time on path, arbitrarily, but long enough
c               for all paths to reach a stationary point, skip out.
c               Nb: Occurs due to restriction here to one level,levpath.

              if (stattim(i,j).gt.stattimmax) then
                xsurf(i,j) = xp
                ysurf(i,j) = yp
                statu(i,j) = up
                statv(i,j) = vp
                statifm(i,j) = 1.
                goto 200
              endif

c++++++++++++++++++
              enddo ! iterlag
c ++++++++++++++++++

c                If exceeded max iter (error), just record current 
c                location and go on.

              xsurf(i,j) = xp
              ysurf(i,j) = yp
              statu(i,j) = up
              statv(i,j) = vp
              statifm(i,j) = -3.
c             write (191,"(a /a,3i6,i8 /a,2f9.0,f8.4,3f10.5 /a,4f10.3)")
c    *          'Warning (path2dome -3):',
c    *          '  i,j,k,iterlag=', i,j,k, iterlag-1,
c    *          '  xp,yp,zp,up,vp,wp=', xp*.001,yp*.001,zp, up,vp,wp,
c    *          '  h[1-4]=', h(iwa,jwa),h(iwb,jwa),h(iwb,jwa),h(iwb,jwb)

  200         continue

              statdist(i,j) = sqrt (  (xh(i,j) - xsurf(i,j))**2
     *                              + (yh(i,j) - ysurf(i,j))**2 )

c             if (mod(i+j+k,100).eq.0) 
c    *          write (6,'(a,6i8)') 
c    *          'path2dome: i,j,iterlag,dtlag,stattim,statdist:',
c    *          i, j, iterlag, nint(dtlag), 
c    *          nint(stattim(i,j)), nint(.001*statdist(i,j))
c>>>>>>>>>>>>>>>>
          endif       ! h > 0
        enddo         ! i
      enddo           ! j
c>>>>>>>>>>>>>>>>

c        Print out diagnostic maps (first call only).
c        stattim = total time (yrs) taken along each path.
c        statnum = number of Lagrangian steps taken along path.
c        statifm: 1=normal exit; -1,-2,etc=error exits.

c~~~~~~~~~~~~~~~~~~~~~~~~~
      if (firstprint) then
c~~~~~~~~~~~~~~~~~~~~~~~~~
        iu = 192
        k=levpath                              ! only one starting layer

        write (iu,'(a,i4)') 'k=',k 
        call printmap (timeice, statdist, 'statdist',100.e3,iu,0)

        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              worka(i,j) = (xsurf(i,j) - xoffa) / dx0
              workb(i,j) = (ysurf(i,j) - yoffa) / dy0
            else
              worka(i,j) = 0.
              workb(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice, worka,   'xsurf',    10.,   iu ,0)
        call printmap (timeice, workb,   'ysurf',    10.,   iu ,0)
        call printmap (timeice, statu,   'statu',    1.,    iu ,0)
        call printmap (timeice, statv,   'statv',    1.,    iu ,0)
        call printmap (timeice, statifm, 'statifm',  1.,    iu ,0)
        call printmap (timeice, stattim, 'stattim',  10000.,iu ,0)
        call printmap (timeice, statnum, 'statnum',  100.,  iu ,0)
        call printmap (timeice, w,       'w (1/yr)', .0001, iu ,0)

        call flush (iu)
        firstprint = .false.
c~~~~~~~~~~
      endif
c~~~~~~~~~~

      return
      end

c-----------------------------------------------------------------------

      subroutine path2domesimp (distdome, h, hs, maskwater)

c     Not currently used:
c     Much simpler and approximate alternative to path2dome.
c     First finds "dome points", i.e., local maxima of hs,
c     and stores their locations in a list (list, nlist).  
c     Then sets distdome for each (non-floating) grid point, 
c     i.e., distance to closest dome - but "closest" is decided by
c     comparing distances weighted by 1/ice thickness (1/h) 
c     of each dome, because ice flowing from thick domes should 
c     radiate further than from small domes.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  distdome(nx,ny), h(nx,ny), hs(nx,ny), maskwater(nx,ny)

      dimension list(nx*ny)

      dimension ish(8), jsh(8)
      save ish, jsh
c     First 4 are E-W/N-S, last 4 are diagonals:
      data ish /-1, 1,  0, 0, -1,  1, -1, 1/
      data jsh / 0, 0, -1, 1, -1, -1,  1, 1/
      parameter (nlook = 8)   ! 4 for E-W/N-S, 8 for diagonals too

      call resetr (distdome, nx*ny, 1.e20)

c        Make list of "domes" (local maxima of surface elevation hs)

      nlist = 0
      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0. .and. maskwater(i,j).eq.0) then

            ifdome = 1
            do look=1,nlook
              ii = max (1, min (nx, i + ish(look)))
              jj = max (1, min (ny, j + jsh(look)))
              if (hs(ii,jj).gt.hs(i,j)) ifdome = 0
            enddo

            if (ifdome.eq.1) then 
              nlist = nlist + 1
              list(nlist) = i + (j-1)*nx
            endif

          endif
        enddo
      enddo

      if (nlist.eq.0) return

c        For each non-ocn point, find distance to closest dome, judging
c        "closest" by comparing distances weighted by 1/h of each dome

      do j=1,ny
        do i=1,nx
          if (maskwater(i,j).eq.0) then
            do n=1,nlist
              ja = (list(n)-1) / nx + 1
              ia =  list(n) - (ja-1)*nx
              zda = dd0 * sqrt(float((i-ia)**2 + (j-ja)**2))      !polst
              if (zda/h(ia,ja).lt.distdome(i,j)) then
                ldome = list(n)
                distdome(i,j) = zda/h(ia,ja)
              endif
            enddo

            jd = (ldome-1) / nx + 1
            id =  ldome - (jd-1)*nx
            distdome(i,j) = distdome(i,j)*h(id,jd)
          endif
        enddo
      enddo

      return
      end

c(((((
#endif
c(((((

c-----------------------------------------------------------------------

      subroutine balvel (h, hsin, budgall, balvelu, balvelv,
     *                   timeice, dtimeice)

c        Calculates balance velocities, following Budd and Warner,
c        1996, Ann. Glac., 23, 21-27. For display only, 
c        if CRHINV or ENHINV. Or for use at calving edges, 
c        if CALVBALANCE.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny),             hsin(nx,ny),          budgall(nx,ny), 
     *  balvelu(0:nxp,0:nyp), balvelv(0:nxp,0:nyp)

      dimension
     *  work(nx*ny),      hs(nx,ny),        indbal(nx*ny),
     *  dfu(0:nxp,0:nyp), dfv(0:nxp,0:nyp),
     *  fu(0:nxp,0:nyp),  fv(0:nxp,0:nyp),
     *  gu(0:nxp,0:nyp),  gv(0:nxp,0:nyp)

c     not used (now based on powb basal sliding):
c     integer powbalh, powbals
cc    parameter (powbalh=5, powbals=3)       ! internal deformation n=3
c     parameter (powbalh=3, powbals=2)       ! basal sliding, ub ~ tau^2

      parameter (dtimebalvel= 1000.)

c       Only do at intervals of dtimebalvel yrs, for CRNINV or ENHINV.
c       Nb: meaningless if first timestep, since budgall and df[u,v]
c       are not set yet (done later in icectl). Also, only accurate
c       after first ~1000's yrs, after base[perc,frml] (in budgall) 
c       are spun up. hbal,hsbal have been set earlier in icectl.
   
#if ! defined (CALVBALANCE)
      dt = dtimeice
      if (mod (abs(timeice)+0.5*dt, max(dtimebalvel,dt)).ge.dt) return
#endif

c        Remove closed depressions (else get "zero zgo" error below),
c        crudely, by iteratively raising grid points with no downhill 
c        path out, to "sill" elev + 1 m. Use local hs, not passed hsin.
 
      call scopy (nx*ny, hsin, 1, hs, 1)
      do iterfix = 1,500
        nfix = 0
        do j=1,ny
          jm = max (j-1,1)
          jp = min (j+1,ny)
          do i=1,nx
            im = max (i-1,1)
            ip = min (i+1,nx)
            if (h(i,j).gt.0.) then
              if ( hs(i,j).le.hs(im,j) .and. 
     *             hs(i,j).le.hs(ip,j) .and. 
     *             hs(i,j).le.hs(i,jm) .and. 
     *             hs(i,j).le.hs(i,jp) ) then
                zhs = hs(i,j)
                hs(i,j) = min (hs(im,j),hs(ip,j),hs(i,jm),hs(i,jp)) + 1.
                nfix = nfix + 1
              endif
            endif
          enddo
        enddo
        if (iterfix.eq.1) nfixa = nfix
        if (nfix.eq.0) goto 500
      enddo
      write (ioterm,'(a)')'*** Error (balvel): still closed depressions'
      stop
  500 continue

c         Sort h-grid ice points into descending-elevation order

      nh = 0
      do j=1,ny  
        do i=1,nx  
          if (h(i,j).gt.0.) then 
            nh = nh + 1
            work(nh) = hs(i,j)
            indbal(nh) = i + (j-1)*nx
          endif
        enddo
      enddo
      call bubblesort_bal (work, indbal, nh) 

c        Set u,v-grid flux weights (estimated from slopes and interped 
c        h). Sign is +ve if *flux* (~ -slope) is increasing in i or j.

      call zero (dfu,  (nxp+1)*(nyp+1))
      call zero (dfv,  (nxp+1)*(nyp+1))
      call zero (gu, (nxp+1)*(nyp+1))
   
      do j=1,ny
        do i=1,nx-1
          zh = 0.5*(h(i+1,j)+h(i,j))
          zdhsu = (hs(i+1,j) - hs(i,j)) / dxu(i,j)
          dfu(i,j) = rhoice*grav*zh*zdhsu
        enddo   
      enddo   
   
      do j=1,ny-1
        do i=1,nx
          zh = 0.5*(h(i,j) + h(i,j+1))
          zdhsv = (hs(i,j+1) - hs(i,j)) / dyv(i,j)
          dfv(i,j) = rhoice*grav*zh*zdhsv
        enddo   
      enddo   

      do j=1,ny
        do i=1,nx-1
          zh = 0.5*(h(i+1,j)+h(i,j))
c either:
c         same u,v-grid aver. of driving (~basal) stress as in icestep:
          zdfv = 0.25
     *         * ( dfv(i,j) +dfv(i+1,j) +dfv(i,j-1) +dfv(i+1,j-1))
          zdfu = dfu(i,j)
          zdf2 = max (zdfu**2 + zdfv**2, .001)
          gu(i,j) = zh * (-zdfu) * (zdf2**((powb-1)/2.))
c or:
c         zs = (hs(i+1,j)-hs(i,j))/dxu(i,j)
c         gu(i,j) = (zh**powbalh) * (abs(zs)**(powbals-1)) * (-zs)
        enddo
      enddo

      call zero (gv, (nxp+1)*(nyp+1))
      do j=1,ny-1
        do i=1,nx
          zh = 0.5*(h(i,j+1)+h(i,j))
c either:
c         same u,v-grid aver. of driving (~basal) stress as in icestep:
          zdfu = 0.25
     *         * ( dfu(i-1,j) +dfu(i,j) +dfu(i-1,j+1) +dfu(i,j+1))
          zdfv = dfv(i,j)
          zdf2 = max (zdfu**2 + zdfv**2, .001)
          gv(i,j) = zh * (-zdfv) * (zdf2**((powb-1)/2.))
c or:
c         zs = (hs(i,j+1)-hs(i,j))/dyv(i,j)
c         gv(i,j) = (zh**powbalh) * (abs(zs)**(powbals-1)) * (-zs)
        enddo
      enddo

c        Calculate u,v-grid balance fluxes (= h*u*interface length)

      call zero (fu, (nxp+1)*(nyp+1))
      call zero (fv, (nxp+1)*(nyp+1))

c--------------
      do k=1,nh     ! process h-grid points, highest to lowest elev
c--------------
        kk = indbal(k)
        j = (kk-1)/nx + 1
        i = kk - (j-1)*nx   

c          Sum incoming fluxes

        zin = 0.
        if (gu(i-1,j).gt.0.) zin = zin  + fu(i-1,j)
        if (gu(i,j)  .lt.0.) zin = zin  - fu(i,j)
        if (gv(i,j-1).gt.0.) zin = zin  + fv(i,j-1)
        if (gv(i,j)  .lt.0.) zin = zin  - fv(i,j)

c          Set total outward flux = incomng flux + local mass balance

        zout = zin + max (budgall(i,j),0.)*dx(i,j)*dy(i,j)

        if (zout.lt.0.) then 
          write (ioterm,'(/a,2i6/a,3e14.6)') 
     *      '*** Error (balvel): negative zout at i,j', i, j,
     *      '    zin, budgall, zout=', zin, budgall(i,j), zout
          stop
        endif

c          Set weights for outgoing fluxes

        zgo = 0.
        if (gu(i-1,j).lt.0.) zgo = zgo + abs(gu(i-1,j))
        if (gu(i,j)  .gt.0.) zgo = zgo + abs(gu(i,j))
        if (gv(i,j-1).lt.0.) zgo = zgo + abs(gv(i,j-1))
        if (gv(i,j)  .gt.0.) zgo = zgo + abs(gv(i,j))

        if (zgo.eq.0.) then 
          write (ioterm,'(/a,2i4)') 
     *      '*** Error (balvel): zero zgo at i,j=',i,j
          stop
        endif

c          Set outgoing fluxes

        if (gu(i-1,j).lt.0.) fu(i-1,j) = -zout * abs(gu(i-1,j))/zgo
        if (gu(i,j)  .gt.0.) fu(i,j)   =  zout * abs(gu(i,j))/zgo
        if (gv(i,j-1).lt.0.) fv(i,j-1) = -zout * abs(gv(i,j-1))/zgo
        if (gv(i,j)  .gt.0.) fv(i,j)   =  zout * abs(gv(i,j))/zgo

c       write (199,'(a,2i6,4f10.4)') 'balvel: i,j,hs,h,zin,zout=',
c    *    i, j, hs(i,j), h(i,j), zin/darea(i,j), zout/darea(i,j)
c----------
      enddo
c----------

c        Calculate u,v-grid balance velocities

      call zero (balvelu, (nxp+1)*(nyp+1))
      do j=1,ny
        do i=1,nx-1
          zh = 0.5*(h(i+1,j)+h(i,j))
          zdy = 0.5*(dy(i+1,j)+dy(i,j))
          if (zh.gt.0.) balvelu(i,j) = fu(i,j) / (zh*zdy)
        enddo
      enddo

      call zero (balvelv, (nxp+1)*(nyp+1))
      do j=1,ny-1
        do i=1,nx
          zh = 0.5*(h(i,j+1)+h(i,j))
          zdx = 0.5*(dx(i,j+1)+dx(i,j))
          if (zh.gt.0.) balvelv(i,j) = fv(i,j) / (zh*zdx)
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine bubblesort_bal (a, ind, n)

c        Bubble sorts a(n) into descending order, 
c        maintaining index ind (i.e., sorted a(i) = a_orig(ind(i)),
c        for balance velocity calc (in balvel).

      dimension a(n), ind(n)

      do ibub = n-1,1,-1
        do j=ibub,n-1
c         if (a(j).gt.a(j+1)) then   ! ascending
          if (a(j).lt.a(j+1)) then   ! descending
            z = a(j+1)
            a(j+1) = a(j)
            a(j) = z

            iz = ind(j+1)
            ind(j+1) = ind(j)
            ind(j) = iz
          else
            go to 10
          endif
        enddo
   10   continue
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine basecoef (h, hb, hbsd, hw, 
     *                     baseperc, basefrml,
     *                     t, tsurfi, sealev,
     *                     crhmel, crhu, crhv,
     *                     fracgu, fracgv, hgu, hgv)

c     Sets basal sliding coefficients crh[u,v] (u,v grids). Setting 
c     starts from intrinsic crhmel (previously set in basecoef_init or
c     read from restart file), and modifies it due to basal temp,
c     water supply, etc. For floating/g.l. pts, uses unmodified crhmel.
c     (Actual values for ice-water interface are set in icedyn).
c
c     These calcs will be expanded when add explicit basal hydrology 
c     (probably on h-grid), and effects on basal hydrologic flow. 

c     Also set fracg[u,v] and hg[u,v] at grounding line.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny),            hb(nx,ny),            
     *  hbsd(nx,ny),         hw(nx,ny),
     *  baseperc(nx,ny),     basefrml(nx,ny),
     *  t(nx,ny,0:nlevp),    tsurfi(nx,ny),
     *  crhmel(nx,ny), 
     *  crhu(0:nxp,0:nyp),   crhv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp), fracgv(0:nxp,0:nyp),
     *  hgu(0:nxp,0:nyp),    hgv(0:nxp,0:nyp)

c        Set basal sliding coefficients crhu (u grid, iloopuv=1), and 
c        crhv (v grid, iloopuv=2), looping over u and v grids to avoid 
c        duplicate code. Do for all points regardless of  water - 
c        values will be re-set for water within icedyn iterations.
c        Also crh[u,v] will be adjusted in sedflow for deforming sed.

      call resetr (crhu, (nxp+1)*(nyp+1), crhnos)
      call resetr (crhv, (nxp+1)*(nyp+1), crhnos)

      call resetr (fracgu, (nxp+1)*(nyp+1), 1.)
      call resetr (fracgv, (nxp+1)*(nyp+1), 1.)

c>>>>>>>>>>>>>>>>>>>>>
      do iloopuv = 1,2
c>>>>>>>>>>>>>>>>>>>>>

      if (iloopuv.eq.1) then
        ish = 1
        jsh = 0
        nxtmp = nx-1
        nytmp = ny
      else
        ish = 0
        jsh = 1
        nxtmp = nx
        nytmp = ny-1
      endif

      do j=1,nytmp
        do i=1,nxtmp
          ia = i + ish
          ja = j + jsh

          zcrhmel = 10. ** (0.5*log10(crhmel(i,j)*crhmel(ia,ja)))

c            If at grounding line, calculate areal fraction of this 
c            u or v-grid cell that's grounded (fracg[u,v]). i.e., 
c            estimated sub-grid position of grounding line. Calculate 
c            height-above-flotation for the two surrounding h-grid 
c            points (one floating with h.a.f. < 0, one grounded with 
c            h.a.f. > 0), and linearly interpolated zero-location
c            is the estimated g.l position. [This is equivalent to 
c            linearly interpolating ice surface, ice base, and bed
c            elevations between the 2 h-points, and seeing where 
c            flotation occurs]. fracg[u,v] is used in icedyn, 
c            multiplying the basal coefficient coefbs[u,v].

c            Also calc g.l. ice thickness hg[u,v] (using fracg[u,v]).
c            Currently fracg[u,v] (and hg[u,v]) are diagnostic only
c            if schoofgl not defined. For schoofgl, fracg[u,v] are
c            used in calcgl to calculate hg[u,v] for schoof g.l. flux
c            (repeating hg[u,v] setting here, but within C loop).

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (     (hw(i,j).gt.hwcut .and. hw(ia,ja).le.hwcut) 
     *        .or. (hw(i,j).le.hwcut .and. hw(ia,ja).gt.hwcut) ) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zpot  = hb(i,j)   + hw(i,j)   + rhor*h(i,j)
            zpota = hb(ia,ja) + hw(ia,ja) + rhor*h(ia,ja)
            if (hw(i,j).gt.hwcut) then
              zaf = zpota - zpot     ! ht above flot (liq m), grounded
              zbf = -hw(i,j)         ! ht above flot (liq m), floating
            else if (hw(ia,ja).gt.hwcut) then
              zaf = zpot - zpota
              zbf = -hw(ia,ja)
            endif
            zfgl = zaf / max (zaf-zbf, 0.1)    ! for safety...should
            zfgl = max (0., min (1., zfgl))  ! have zaf>0, zbf<0
c           zfgl = 0.5 ! 777

c           or as in Pattyn et al (2006) (doesn't work if h=0):  
c           if (hw(i,j).gt.hwcut) then
c             zf = (sealev - hb(i,j))   / (rhor*max (h(i,j),.1))
c             zg = (sealev - hb(ia,ja)) / (rhor*max (h(ia,ja),.1))
c           else
c             zg = (sealev - hb(i,j))   / (rhor*max (h(i,j),.1))
c             zf = (sealev - hb(ia,ja)) / (rhor*max (h(ia,ja),.1))
c           endif
c           zfgl = (1.-zg) / max (zf-zg, .01)

c           calc g.l. depth (lin.interp bed elev to u,v grid, vs. s.l.,
c           so only works for marine g.l.s, not inland lakes):
c           zhbg =       0.5*hb(i,j) +  0.5*hb(ia,ja) 
c           slightly better in 100/300ka tests(?):
            if (hw(i,j).gt.hwcut) then
              zhbg = (1.-zfgl)*hb(ia,ja) + zfgl*hb(i,j) 
            else
              zhbg = (1.-zfgl)*hb(i,j)   + zfgl*hb(ia,ja) 
            endif
            zhg  = (sealev-zhbg)/rhor
#if defined (KN_GLSLIDE)            
c            increase sliding coefficient at grounding line (Knut,Nick):
             zcrhmel = min ( 1.e-4, zcrhmel*(10.**KN_GLSLIDE) ) 
#endif
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          else if (hw(i,j).gt.hwcut .and. hw(ia,ja).gt.hwcut) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            zfgl = 0.
            zhg  = 0.
c~~~~~~~~~~~~~
          else
c~~~~~~~~~~~~~
            zfgl = 1.
            zhg  = 0.
c~~~~~~~~~~~~~~
          endif
c~~~~~~~~~~~~~~

          if (iloopuv.eq.1) then
            fracgu(i,j) = zfgl
            hgu(i,j)    = zhg
          else
            fracgv(i,j) = zfgl
            hgv(i,j)    = zhg
          endif

c....

c            more influence of adjacent floating base:

          zt = 0. 
          nzt = 0 
          if (h(i,j).gt.0.) then
            zt = zt + t(i,j,nlevp)-(tmelt-dtmdh*h(i,j))
            nzt = nzt + 1 
          else if (hw(i,j).gt.hwcut) then
            zt = zt + 0. 
            nzt = nzt + 1 
          endif   
   
          if (h(ia,ja).gt.0.) then
            zt = zt + t(ia,ja,nlevp)-(tmelt-dtmdh*h(ia,ja))
            nzt = nzt + 1 
          else if (hw(ia,ja).gt.hwcut) then
            zt = zt + 0. 
            nzt = nzt + 1 
          endif   
   
          if (nzt.ge.1) then
            zt = zt/nzt
          else    
            zt = 0.5*(tsurfi(i,j) + tsurfi(ia,ja)) - tmelt 
          endif   

c....

c         influence of basal temperature, h2o on sliding: 
          call calc_crh (zcrh, fract, fracb,
     *                   zcrhmel, 
     *                   zt, 
     *                   0.5*(  baseperc(i,j)  -basefrml(i,j)
     *                        + baseperc(ia,ja)-basefrml(ia,ja)),
     *                   0.5*(h(i,j)+h(ia,ja)),
     *                   0.5*(hbsd(i,j)+hbsd(ia,ja)), 
     *                   0)

c         if floating or at grounding line, no influence:
          if (zfgl.lt.1.) then
            fract = 1.
            fracb = 1.
            zcrh = zcrhmel
          endif

          if (iloopuv.eq.1) then
            crhu(i,j) = zcrh
          else    
            crhv(i,j) = zcrh
          endif   

        enddo  ! i=1,nxtmp
      enddo    ! j=1,nytmp

c>>>>>>>>>>
      enddo   !iloopuv loop
c>>>>>>>>>>

      return
      end

c-----------------------------------------------------------------------

      subroutine calc_crh (crh, fract, fracb,
     *                     crhmel, tb, bh2o, h, hbsd, iflag)

c     Sets current basal sliding coefficient (crh) and also 
c     fractional influences of temp (fract, 0-1) and h2o (fracb, 0-inf),
c     given: 
c       crhmel = intrinsic basal sliding coeff
c       tb = basal temperature relative to pressure melt pt (deg C)
c       bh2o = basal water supply (baseperc-basefrml, m/y)
c       h = ice thickness (m)
c       hbsd = sub-grid std. dev. of bed topography (m)
c       iflag = 0 to ignore tb, bh2o influences if NOBASET(basecoef...)
c               1 to calculate anyway (iceshow2d, writehis)

#include <comicephys.h>

      crh = crhmel
      fract = 1.
      fracb = 1.

#if defined (NOBASET) 
      if (iflag.eq.0) return
#endif

c     effect of basal temperature:
      ztramp = max (3., .02*hbsd)
      fract = max (0., min (1., (tb+ztramp) / ztramp ))
      crh_t = 10. ** (fract*log10(crhmel) +  (1.-fract)*log10(crhnos))

c     effect of basal h2o supply (bh2o, m/yr):
cc    fracb =  bh2o * 1.e3
cc    crh_b = max (fracb*crhmel, crhnos)
      fracb = 1. ! no bh2o effect
      crh_b = 0. ! no bh2o effect

c+++
c     Runbasrcm4xp3:
c     fracb =  bh2o * 1.e2
c     crh_b = max (fracb*crhmel, crhnos)
c+++

      crh = max (crh_t, crh_b)
      crh = min (crh, crhsed)

      return
      end

c-----------------------------------------------------------------------

c******************
#if defined (MOVEW)
c******************

      subroutine movewater (maskwater, h, hb, hw, tw, 
     *                      heatw, baseperc, basefrml,
     *                      uw, vw, dt)

c Advects liquid water thickness (hw) due to horizontal flux between 
c adjacent cells, and applies local source/sink baseperc-basefrml. 
c (For now, no P-E on open water surfaces).
c
c Also advects water temperature (tw) for large water, or sets 
c heating advection tendency (heatw) for small water,
c used in vdif (where small-water tw is set).
c
c Thickness advection: (adapted from channelcompute in hydrology model)
c ~~~~~~~~~~~~~~~~~~~~
c Flow can be any of 8 quadrants. Also, applies negative volume 
c correction, tracked by enegsub, and computes diagnostic totosub, 
c flow from land to ocean. Flow equation is:  
c
c   d(hw)/dt = d/dx [ C(hw) hw^powdepth d(hb+hw+rhor*h]/dx^powslope ]
c
c Flow can occur out of a box in all directions, to avoid flip-flops
c at river junctions when just allowing only 1 (steepest) direction.
c Do time implicitly just with respect to hw term in the water
c potential slope.   
c
c 2 methods (if defines):
c channel_tridia:
c Implicit with contributions from all surrounding points. Tridiagonal
c solution, only for 1-D flowline. (6 out of the 8 directions in 
c iq,jq,etc, are wasted).
c
c channel_sparse:
c As channel_tridia but using sparse solver linbcg (Numerical Recipes)
c and row-indexed sparse storage method.
c
c Time implicit terms use linearized contribution from the slope term 
c in the flow equation. The exponent (slope**powslope) is set to 1 for 
c small slopes (< slopecrit) to avoid infinities at zero slope for 
c powslope < 1).
c
c Use nuptot and mupind to do this, index from full arrays into 
c linear sequence of selected points. 

c Temperature advection: 
c ~~~~~~~~~~~~~~~~~~~~~~
c Done at the end, time explicitly w.r.t. tw, in volume-flux form 
c using stored fluxes of hw across each cell face from hw calcs.

#include <comicephys.h>
#include <comicegrid.h>
#include <comicesparse.h>

      dimension
     *  maskwater(nx,ny), h(nx,ny),  hb(nx,ny), hw(nx,ny), tw(nx,ny),
     *  heatw(nx,ny), baseperc(nx,ny), basefrml(nx,ny),
     *  uw(0:nxp,0:nyp), vw(0:nxp,0:nyp) 

c Local arrays:
      dimension 
     *  hwprev(nx,ny), potw(nx,ny), zdt(nx,ny), 
     *  zfstora(nx,ny,8), zfstorb(nx,ny,8), work(nx,ny),
c       diagnostic: 
     *  streamflux(nx,ny),  ! stream flux out of cell (m3/s)
     *  istreamflux(nx,ny)  ! quadrant for streamflux (1-8)

#if defined (CHANNEL_TRIDIA)
      dimension atri(nx), btri(nx), ctri(nx), vech(nx), rhs(nx)
#elif defined (CHANNEL_SPARSE)
      dimension vech(nuvmax)          ! elspa,ijspa,etc, in comicesparse
      dimension rhs(nx*ny)
#endif
      dimension mupind(nx,ny), rhsexp(nx,ny)

      dimension iq(8), jq(8), kop(8), dists(8), dista(8)
      save iq, jq, kop
c directions to adjacent cells, in common in hydgrid.h:
c     k:         1   2   3   4   5   6   7   8
c anticlockwise from E:
c     data iq /  1,  1,  0, -1, -1, -1,  0,  1/
c     data jq /  0,  1,  1,  1,  0, -1, -1, -1/
c     data kop/  5,  6,  7,  8,  1,  2,  3,  4/
c sparse order (same as in sparse matrix in channeclcompute,
c lowest to highest i + (j-1)*nx):
      data iq / -1,  0,  1, -1,  1, -1,  0,  1/
      data jq / -1, -1, -1,  0,  0,  1,  1,  1/
      data kop/  8,  7,  6,  5,  4,  3,  2,  1/

      parameter (hwmin = hwcut)
      parameter (slopecrit = 1.e-10)
      logical first
      data first /.true./

c     statement function for flow coefficient:
      flowcoeff (d, dmin) = 1.  ! need to adjust as fn of d       777   

      call scopy (nx*ny, hw, 1, hwprev, 1)

c        Set nuptot and index mupind 

      nuptot = 0
      do j=1,ny
        do i=1,nx
#if defined (CHANNEL_TRIDIA) 
          if ( .true. ) then
#elif defined (CHANNEL_SPARSE)
          if ( .true. ) then  ! perhaps elim hw=0 & deep ocn pts  777   
#endif
            nuptot = nuptot + 1
            mupind(i,j) = nuptot
          else
            mupind(i,j) = 0
          endif
        enddo
      enddo
c     if (first) then
c       write (ioterm,*)
c    *     'movewater: nuptot=',nuptot, ' nx*ny=',nx*ny
c       first = .false.
c     endif

      call zero   (zfstora, nx*ny*8)
      call zero   (zfstorb, nx*ny*8)
      call resetr (zdt,     nx*ny, 1.)
      call zero   (rhsexp,  nx*ny)

#if defined (CHANNEL_TRIDIA)
      call zero (rhs,  nx)
      call resetr (atri, nx, 1.)
      call zero (btri, nx)
      call zero (ctri, nx)
#elif defined (CHANNEL_SPARSE)
      call zero (rhs, nuptot)
      call izero (ijspa, nspamax)
      call zero  (elspa, nspamax)
      do n1=1,nuptot
        elspa(n1) = 1.
      enddo
#endif

      do j=1,ny
        do i=1,nx
          zdt(i,j) = dt / darea(i,j)
          zlocal = dt * (baseperc(i,j) - basefrml(i,j)) !negl P-E 777   
     *             * (rhoice/rholiq)
          rhsexp(i,j) = rhsexp(i,j) + zlocal
          n1 = mupind(i,j)
          if (n1.ge.1) then
            rhs(n1) = rhs(n1) + zlocal
          endif
          potw(i,j) = hb(i,j) + hw(i,j) + rhor*h(i,j)
        enddo
      enddo

#if defined (OMPHYD)
!$OMP PARALLEL DO
!$OMP* shared  (iq,jq,dx,dy,hw,potw,
!$OMP*          zfstora,zfstorb)
!$OMP* private (i,j,i2,j2,k,dists,dista,zslope, zsabs,
!$OMP*          zdep,zwid,zcros,
!$OMP*          zcoef,zf,powslope,powconst,powdepth)
!$OMP* schedule(dynamic,2)
#endif

c        Compute zfstora, time-explicit coefficient multiplied
c        by potw(i2,j2)-potw(i,j)), for (i,j) and all its 8 
c        neighbors. Also zfstorb, time-implicit term that will be
c        multiplied by hw(i2,j2)-hw(i,j). To save cpu, only compute 
c        for 4 directions at each i,j, which gets all pairs in domain.
c        In subsequent do loops, all 8 directions used, so in next do 
c        loop below, zfstor[a,b](i,j,[5-8]) are copied from ...[1-4].

      do j=1,ny
        do i=1,nx

c            Look at 4 quadrants (gets all connections over domain
c            just once each)...save in zfstora,b

c~~~~~~~~~~~~~~~~~
          do k=1,4
c~~~~~~~~~~~~~~~~~
            i2 = i + iq(k)
            j2 = j + jq(k)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (i2.ge.1.and.i2.le.nx.and.j2.ge.1.and.j2.le.ny) then
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              if (j.eq.j2) then
                dists(k) = 0.5*(dx(i,j)+dx(i2,j2))
                dista(k) = 0.5*(dy(i,j)+dy(i2,j2))
              else if (i.eq.i2) then
                dists(k) = 0.5*(dy(i,j)+dy(i2,j2))
                dista(k) = 0.5*(dx(i,j)+dx(i2,j2))
              else 
                dists(k) = sqrt (   (0.5*(dx(i,j)+dx(i2,j2)))**2 
     *                            + (0.5*(dy(i,j)+dy(i2,j2)))**2 )
                dista(k) = 0.25*(dx(i,j)+dx(i2,j2)+dy(i,j)+dy(i2,j2))
              endif

              zslope = potw(i2,j2)-potw(i,j) / dists(k)
              zsabs  = abs(zslope)

              if (zslope.gt.0.) then              ! use upwind thickenss
                zdep = hw(i2,j2) 
              else
                zdep = hw(i,j) 
              endif
              zwid = dista(k)
              zcros = zdep * zwid

c Either Manning:           ! use for ice-free hydrol, rivers ?   777   
c             if (h(i,j).eq.0. .and. h(i2,j2).eq.0.) then
c               zcoef = 2.5/.06    ! 2.5/.15 
c               powdepth = 2./3.
c               if (zsabs.gt.slopecrit) then
c                 powslope = 0.5
c                 powconst = 1.
c               else
c                 powslope = 1.
c                 powconst = slopecrit**(powslope - 0.5)
c               endif
c             else
c Or sub-ice, lakes:
                zcoef = flowcoeff (zdep, hwmin)
                powdepth = 0.
                powslope = 1.
                powconst = 1.
c             endif
              
              zf = ( zcoef * zcros * zdep**powdepth )
     *             * zsabs**(powslope-1.) / (powconst*dists(k))
              zfstora(i,j,k) = zf
              zfstorb(i,j,k) = zf*powslope

c++++++++++++++++
            endif
c++++++++++++++++

c~~~~~~~~~~~~~~
          enddo
c~~~~~~~~~~~~~~

        enddo
      enddo

c        Set arrays and rhs's needed to solve system of eqns for 
c        each method, using zfstor[a,b] from above. 
c        rhs, *tri, *spa are for subset of pts solved
c        time implictly (index n1,n2), and rhsexp(i,j) are for 
c        other points, solved time explicitly.
c        Can't parallelize for sparse method, since off-diag terms
c        in ijspa,elspa are assembled using incremented index kof.   

      kof = nuptot + 1
      do j=1,ny
        do i=1,nx
          kofprev = kof
          n1 = mupind(i,j)
          do k=1,8
            i2 = i + iq(k)
            j2 = j + jq(k)
            if ( (i2.ge.1.and.i2.le.nx.and.j2.ge.1.and.j2.le.ny) ) then
              n2 = mupind(i2,j2)

c             fill in other 4 directions of zfstora,b:
              if (k.gt.4) then
                zfstora(i,j,k) = zfstora(i2,j2,kop(k))
                zfstorb(i,j,k) = zfstorb(i2,j2,kop(k))
              endif

              rhsexp(i,j) = rhsexp(i,j)  + zdt(i,j)  * zfstora(i,j,k)
     *                                   * (potw(i2,j2)-potw(i,j))

              if (n1.ge.1) then
#if defined (CHANNEL_TRIDIA)
                rhs(n1)  = rhs(n1)       + zdt(i,j)  * zfstora(i,j,k)
     *                                   * (potw(i2,j2)-potw(i,j))
                atri(n1) = atri(n1) + zdt(i,j) * zfstorb(i,j,k)
c               all except k=2 and k=7 are not used (n2=0 for others): 
                if (i2.eq.i+1) then       ! superdiagonal
                  btri(n2) = btri(n2) - zdt(i,j) * zfstorb(i,j,k)
                else if (i2.eq.i-1) then  ! subdiagonal
                  ctri(n2) = ctri(n2) - zdt(i,j) * zfstorb(i,j,k)
                endif
#elif defined (CHANNEL_SPARSE)
                rhs(n1)  = rhs(n1)     + zdt(i,j) * zfstora(i,j,k)
     *                                 * (potw(i2,j2)-potw(i,j))
                elspa(n1) = elspa(n1)  + zdt(i,j) * zfstorb(i,j,k)
                ijspa(n1) = kofprev + 1 
                if (n2.ge.1) then
                  kof = kof + 1
                  elspa(kof) = elspa(kof) - zdt(i,j)*zfstorb(i,j,k)
                  ijspa(kof) = n2
                endif
#endif
              endif

            endif
          enddo
        enddo
      enddo

#if defined (CHANNEL_SPARSE)
      ijspa(nuptot+1) = kof + 1
#endif

c        Solve linear system arr*vech = rhs for vech (=hw-hwprev)

      call zero (vech, nuptot)

#if defined (CHANNEL_TRIDIA)

      call tridia_i (atri, btri, ctri, vech, rhs, nuptot, 0)
c     call tridia_i (atri, btri, ctri, vech, rhs, nuptot, 1) !777

#elif defined (CHANNEL_SPARSE)

      itol = 1
      tol = 1.e-11
      itmax = 10000   ! 1000
      call linbcg (nuptot, rhs, vech, itol, tol, itmax, itlin, errlin)
      if (abs(errlin).gt.2.*tol) write (ioterm,*)
     *  '*** Warning linbcg(hw): itlin=',itlin,' errlin=',errlin

#endif

c        Set new water thickness hw

      do j=1,ny
        do i=1,nx
          n1 = mupind(i,j)
          if (n1.ge.1) then
            hw(i,j) = hwprev(i,j) + vech(n1)
          else
            hw(i,j) = hwprev(i,j) + rhsexp(i,j)
          endif
        enddo
      enddo

c        Adjust time explicit pts for neighboring time implicit changes
c        (conserves water)

      do j=1,ny
        do i=1,nx
          n1 = mupind(i,j)
          if (n1.eq.0) then
            do k=1,8
              i2 = i + iq(k)
              j2 = j + jq(k)
              if (i2.ge.1.and.i2.le.nx.and.j2.ge.1.and.j2.le.ny) then
                n2 = mupind(i2,j2)
                if (n2.ge.1) then
                  hw(i,j)= hw(i,j) + zdt(i,j)*zfstorb(i,j,k)
     *                               *(hw(i2,j2)-hwprev(i2,j2))
                endif
              endif
            enddo
          endif
        enddo
      enddo

c        Crude fix for negative water thickness

      enegsub = 0.
      do j=1,ny
        do i=1,nx
          if (hw(i,j).lt.0.) then
c           write (ioterm,*) 'neg vol i,j=',i,j
            enegsub = enegsub + (-hw(i,j)) * darea(i,j) / dt
            hw(i,j) = 0.
          endif
        enddo
      enddo

c        Set totosub (diagnostic, total stream flux from land to ocean,
c        m3/y), and streamflux,istreamflux (diagnostic, total stream 
c        flux (m3/y) and average quadrant *out* of each grid cell).
c
c        Also compute heat advection out of the cell (zhadv), 
c        and step temperature tw due to advection, or set tendency
c        heatw for small water.

      totosub = 0.
      call zero (streamflux, nx*ny)
      call zero (work, nx*ny)
      call zero (uw, (nxp+1)*(nyp+1))
      call zero (vw, (nxp+1)*(nyp+1))

      do j=1,ny
        do i=1,nx
          n1 = mupind(i,j)

          zhadv = 0.
          do k=1,8
            i2 = i + iq(k)
            j2 = j + jq(k)
            if ( (i2.ge.1.and.i2.le.nx.and.j2.ge.1.and.j2.le.ny) ) then
              n2 = mupind(i2,j2)

c                zs is h(upstream)*width*veloc out of each face (m3/y)

              zs = zfstora(i,j,k) * (potw(i,j)-potw(i2,j2))

c             correction for time implicit:
              if (n1.ge.1) then
                zs = zs + zfstorb(i,j,k)*vech(n1)
                if (n2.ge.1) zs = zs - zfstorb(i,j,k)*vech(n2)
              endif

              if (zs.gt.0.) then
                streamflux(i,j) = streamflux(i,j) + zs
                work(i,j) = work(i,j) + k*zs
              else 
                streamflux(i2,j2) = streamflux(i2,j2) - zs
                work(i2,j2) = work(i2,j2) - kop(k)*zs
              endif

c             land/lake to ocean:
              if (maskwater(i,j).ne.1.and.maskwater(i2,j2).eq.1) then
                totosub = totosub + zs
              endif

c                zhadv is T(upstream)*h*width*veloc out of cell (K m3/y)

              zhadv = zhadv - max(zs,0.)*tw(i,j) - min(zs,0.)*tw(i2,j2) 

c                Deduce water velocities uw,vw on u and v-grids 
c                (just set velocs at the rh and top faces of each 
c                h-grid cell, this will get all u and v-grid points).
c                Do only if surrounding hw's are large (used in icedyn).
c                Could try to factor in diagonal velocities ?     777   

              if (i2.eq.i+1 .and. j2.eq.j) then
                if (hw(i,j).gt.hwmin .and. hw(i2,j2).gt.hwmin) then
                  zdist =  0.5*(dx(i,j)+dx(i2,j2))
                  if (zs.gt.0.) then ! upstream
                    zdep = hwprev(i,j)
                  else
                    zdep = hwprev(i2,j2)
                  endif
                  uw(i,j) = zs / ( max(zdep,hwmin) * zdist )
                endif

              else if (i2.eq.i .and. j2.eq.j+1) then
                if (hw(i,j).gt.hwmin .and. hw(i2,j2).gt.hwmin) then
                  zdist =  0.5*(dy(i,j)+dy(i2,j2))
                  if (zs.gt.0.) then ! upstream
                    zdep = hwprev(i,j)
                  else
                    zdep = hwprev(i2,j2)
                  endif
                  vw(i,j) = zs / ( max(zdep,hwmin) * zdist )
                endif
              endif

            endif
          enddo

c            Set new water temperature (large) or heat tendency (small)

          if (hw(i,j).gt.hwmin) then
            ztme = tmelt - dtmdh*h(i,j)
            ztnew =  ( tw(i,j)*hwprev(i,j) + zhadv*zdt(i,j) 
     *                 + ztme*(baseperc(i,j)-basefrml(i,j))*dt 
     *                       *(rhoice/rholiq)
     *               ) / hw(i,j)
            tw(i,j) = ztnew
            heatw(i,j) = 0.
          else if (hw(i,j).gt.1.e-3) then
c           small and no water tw will be set in vdif
            heatw(i,j) = rholiq*cheatliq*zhadv*zdt(i,j) / dt    ! J/m2/y
          else
            heatw(i,j) = 0.
          endif

        enddo
      enddo

      do j=1,ny
        do i=1,nx
          if (streamflux(i,j).gt.0.) then
            istreamflux(i,j) = nint(work(i,j) / streamflux(i,j))
          else
            istreamflux(i,j) = 0
          endif
        enddo
      enddo

      return
      end

c*****
#endif
c*****
