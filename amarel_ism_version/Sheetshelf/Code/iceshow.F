      subroutine iceshow1d (h, hs, hb, t,
     *                      basefrml, oceanmelt, budgall, budgneg, 
     *                      facemelt, cliffmelt, div, arhap, calvice,
     *                      bedge, fedge, hedge,
     *                      tsurf, tsurfi, heati, heath, heatb,
     *                      w, wa, maskh,
     *                      sedim, tsed, wsed, heats, 
     *                      quarryrate, quarryacc, quarrycoef,
     *                      pelagrate, pelagacc,
     *                      topbed, topbedeq, tbed,
     *                      hw, tw, maskwater,
     *                      u, v, ua, va, uia, via, ub, vb,
     *                      uadv, vadv, hu, hv, dfu, dfv, dbu, dbv,
     *                      fracgu, fracgv, muind, mvind, thetau,thetav,
     *                      crhu, crhv, fsedu, fsedv, 
#if defined (MELANGE)
     *                      hmel, hwmel, fmel, maskmel, 
     *                      budgallmel, budgsupmel,
     *                      uamel, vamel, 
     *                      bpumel, bpvmel, 
     *                      thetaumel, thetavmel,
#endif
     *                      numh, nuvtot, itera, iterc,
     *                      toth0, tota0, toth, tota, totflow,
     *                      sealev, dtantann, dtantjan, dtseas, rco2,
     *                      timeice, dt, weirun, 
     *                      iffirse, iflast, nyearout1d)

c     Ascii 1-D printout - flowlines, within 1-D or 2-D domains.
c     Flowline for 2-D is defined by i1s,j1s,n1s (in comicegrid),
c     set by setflowline in first call below.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  t(nx,ny,0:nlevp),
     *  basefrml(nx,ny),
     *  oceanmelt(nx,ny),    budgall(nx,ny),     budgneg(nx,ny),
     *  facemelt(nx,ny),     cliffmelt(nx,ny),   
     *  div(nx,ny),          arhap(nx,ny),       calvice(nx,ny),
     *  bedge(nx,ny),        fedge(nx,ny),       hedge(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  heati(nx,ny,nlev),   heath(nx,ny),       heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  maskh(nx,ny),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),   
     *  heats(nx,ny,nsed),   
     *  quarryrate(nx,ny),   quarryacc(nx,ny),   quarrycoef(nx,ny), 
     *  pelagrate(nx,ny),    pelagacc(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  hw(nx,ny),           tw(nx,ny),          maskwater(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),       via(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp),
     *  dbu(0:nxp,0:nyp),       dbv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),    fracgv(0:nxp,0:nyp),
     *  muind(0:nxp,0:nyp),     mvind(0:nxp,0:nyp), 
     *  thetau(nx,ny),          thetav(nx,ny),
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp),
     *  fsedu(0:nxp,0:nyp),     fsedv(0:nxp,0:nyp)
#if defined (MELANGE)
     * ,hmel(nx,ny),            hwmel(nx,ny),
     *  fmel(nx,ny),            maskmel(nx,ny),     
     *  budgallmel(nx,ny),      budgsupmel(nx,ny),
     *  uamel(0:nxp,0:nyp),     vamel(0:nxp,0:nyp),
     *  bpumel(0:nxp,0:nyp),    bpvmel(0:nxp,0:nyp),
     *  thetaumel(nx,ny),       thetavmel(nx,ny)
#endif

      logical first1s
      data first1s /.true./
      save first1s

      dimension 
     *  work1d(nx+ny),
     *  work1dv(0:nx+ny+1),   iwork1dv(0:nx+ny+1), 
     *  work2d(0:nxp,0:nyp),  work2da(0:nxp,0:nyp)
      character*5 cwork1d(nx+ny), cwork1da(nx+ny)


c       Skip out if not time

      if (nyearout1d.eq.0) return
      if ( .not.
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyearout1d),dt)).lt.dt
     *       .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyearout1d.eq.-1 )
     *   ) return 

c        Overall mass budget 

      totbud  = 0.
      totbudi = 0.
      totnega = 0.
      do j=1,ny
        do i=1,nx
          totbud = totbud + budgall(i,j) * darea(i,j) * dt
          if (maskh(i,j).eq.1)
     *      totbudi = totbudi + budgall(i,j) * darea(i,j) * dt
          totnega = totnega + budgneg(i,j) * darea(i,j) * dt
        enddo
      enddo
      toterr = (toth-toth0) - (totbud + totnega + totflow) 

c        Print general 0-D quantities  

#if defined (EISLINE) && defined (LINEM) && defined (RESOL6)
      do iprint=1,3
        if (iprint.eq.1) iu = 6
        if (iprint.eq.2) iu = iuout1d
        if (iprint.eq.3) iu = iuout1dzoom
#else
      do iprint=1,2
        if (iprint.eq.1) iu = 6
        if (iprint.eq.2) iu = iuout1d
#endif
        write (iu,'(/80("=")/
     *            a,i10, a,f7.4, a,i6, a,2i4, a,i5, 3(a,f8.2), a,f8.3)')
     *    'time(yr)=',     nint(timeice),
     *    '   weirun=',    weirun,
     *    '   nuvtot=',    nuvtot,
     *    '   iter[a-c]=', itera,iterc,
     *    '   sealev=',    nint(sealev),
     *    '   dtantann=',  dtantann,
     *    '   dtantjan=',  dtantjan,
     *    '   dtseas  =',  dtseas,
     *    '   rco2    =',  rco2
        zz = dt*totarea
        write (iu,'(a,5f14.8,f21.15)')
     *   'd(h), bud, budi, flow, neg, err =',
     *    (toth-toth0)/zz, totbud/zz, totbudi/zz, totflow/zz,totnega/zz,
     *    toterr/zz
        call flush (iu)
      enddo

c        Print 1-D or flowline output 

c=======================================================================
#if defined (EISLINE) || defined (EISANTA) || defined (TEST2D) || defined (TESTMEL2D)
c=======================================================================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if (defined (EISANTA) && ! defined (TRANSECTA))
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do iprint=1,1
        iu = iuout1d
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#elif defined (EISLINE) || (defined (EISANTA) && defined (TRANSECTA)) || defined (TEST2D) || defined (TESTMEL2D)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  if defined (LINEM) && defined (RESOL6)
      do iprint=1,3
        if (iprint.eq.1) iu = 6
        if (iprint.eq.2) iu = iuout1d
        if (iprint.eq.3) iu = iuout1dzoom
#  else
      do iprint=2,2
        if (iprint.eq.1) iu = 6
        if (iprint.eq.2) iu = iuout1d
#  endif
c~~~~~
#endif
c~~~~~

c          Set 1-D flowline indices within full 1-D or 2-D domain.
c          i1s,j1s,n1s are in common (comicegrid).

        if (first1s) then
c...............................................
#if (defined (EISANTA) && ! defined (TRANSECTA))
c...............................................
          call setflowline
          first1s = .false.
c............................................
#elif defined (TEST2D) || defined (TESTMEL2D)
c............................................
          n1s = 0 
          do i = 1,nx
            n1s = n1s + 1 
            i1s(n1s) = i
            j1s(n1s) = (ny+1)/2
            d1s(n1s) = xh(i1s(n1s),j1s(n1s))
          enddo
c....
#else
c....
          n1s = 0 
          if (nx.ge.ny) then
            do i = 1,nx
              n1s = n1s + 1 
              i1s(n1s) = i
              j1s(n1s) = (ny+1)/2
#  if defined (EISANTA) && defined (TRANSECTA)
              d1s(n1s) = (i-0.5)*(totlen/nx)          ! totlen in common
#  else
              d1s(n1s) = xh(i1s(n1s),j1s(n1s))
#  endif 
            enddo
          else
            do j = 1, ny 
              n1s = n1s + 1 
              i1s(n1s) = (nx+1)/2
              j1s(n1s) = j
              d1s(n1s) = yh(i1s(n1s),j1s(n1s))
            enddo
          endif
#  if defined (LINEM) && defined (RESOL6)
          if (iu.eq.iuout1dzoom) then
            j = (ny+1)/2 
            do i=1,nx
              if (fracgu(i,j).lt.1.) then 
                igl = i
                go to 50
              endif
            enddo
            igl = nx/2
   50       continue
            n1s = 0 
            n1sa = max ( 1, min(igl-50,nx-100))
            n1sb = min (nx, n1sa+100)
            do i = n1sa,n1sb
              n1s = n1s + 1 
              i1s(n1s) = i
              j1s(n1s) = (ny+1)/2
              d1s(n1s) = xh(i1s(n1s),j1s(n1s))
            enddo
          endif
#  else
          first1s = .false.
#  endif
c.....
#endif
c.....
        endif

        inc = max (n1s/100, 1)            ! skip points for ok lwy width

        write(iu,'(a,1001i5)') 'i1s :', (i1s(m),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'j1s :', (j1s(m),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'km  :', (nint(d1s(m)*1.e-3),m=1,n1s,inc)
        write(iu,*)
        write(iu,'(a,1001i5)') 'h   :',
     *                              (nint(h(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'hs  :',
     *                             (nint(hs(i1s(m),j1s(m))),m=1,n1s,inc)
c       write(iu,'(a,1001i5)') 'hbi :',
c    *            (nint(hs(i1s(m),j1s(m))-h(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'hb  :',
     *                             (nint(hb(i1s(m),j1s(m))),m=1,n1s,inc)
 
        do m=1,n1s
          if (sedim(i1s(m),j1s(m)).gt.0.949) then
            write (cwork1d (m),'(i5)') nint(sedim(i1s(m),j1s(m)))
c           place holder, old sedimun:
            write (cwork1da(m),'(i5)') nint(sedim(i1s(m),j1s(m)))
          else
            write (cwork1d (m),'(f5.1)') sedim(i1s(m),j1s(m))
c           place holder, old sedimun:
            write (cwork1da(m),'(f5.1)') sedim(i1s(m),j1s(m))
          endif
        enddo
        write(iu,'(a,1001a5)') 'sed :', (cwork1d(m),m=1,n1s,inc)
        write(iu,'(a,1001a5)') 'sedu:', (cwork1da(m),m=1,n1s,inc)

        write(iu,'(a,1001i5)') 'topb:',
     *                         (nint(topbed(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'topq:',
     *                       (nint(topbedeq(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'hw  :',
     *                             (nint(hw(i1s(m),j1s(m))),m=1,n1s,inc)
#if defined (LINEE)
        write(iu,'(a,1001i5)') 'hwX :',
     *                        (nint(100.*hw(i1s(m),j1s(m))),m=1,n1s,inc)
#endif
        write(iu,'(a,1001i5)') 'tw  :',
     *                 (nint(10.*(tw(i1s(m),j1s(m))-tmelt)),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'mask:',
     *                            (maskwater(i1s(m),j1s(m)),m=1,n1s,inc)

        call iflowline (muind, mvind, iwork1dv, 0, 0)
        write(iu,'(a,2x,1001i5)') 'muin:', (iwork1dv(m),m=1,n1s-1,inc)

        write(iu,'(a,1001i5)') 'zflo:',
     *    ( nint (max ( 0., h(i1s(m),j1s(m)) 
     *                      - (sealev-hb(i1s(m),j1s(m)))/rhor )),
     *      m=1,n1s,inc )
#if defined (LINEE)
        write(iu,'(a,1001i5)') 'zflX:',
     *    ( nint ( 100. * max (0., h(i1s(m),j1s(m))
     *                            - (sealev-hb(i1s(m),j1s(m)))/rhor)),
     *      m=1,n1s,inc )
#endif

c       like flowline, but for h-grid thetau,v:
        do m=1,n1s
          if (j1s(m).eq.j1s(min(n1s,m+1))) then
            work1d(m) = thetau(i1s(m),j1s(m))
          else
            work1d(m) = thetav(i1s(m),j1s(m))
          endif
        enddo
        write(iu,'(a,1001i5)') 'thet:', 
     *                                (nint(work1d(m)*1.e2),m=1,n1s,inc)

        write(iu,'(a,1001i5)') 'budg:',
     *     ( max(-9999,nint(1000.*budgall(i1s(m),j1s(m)))) ,m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'clif:',
     *   ( max(-9999,nint(1000.*cliffmelt(i1s(m),j1s(m)))) ,m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'calv:',
     *     ( max(-9999,nint(1000.*calvice(i1s(m),j1s(m)))) ,m=1,n1s,inc)

c       sfc+basal crevasse depth:ice thickness (Nick, as subr calving):
        do m=1,n1s
          i = i1s(m)
          j = j1s(m)
          if (maskwater(i,j).eq.1 .and. hedge(i,j).gt.0.) then
            zhdiv = 2. * ( (max(div(i,j),0.)**powir)/arhap(i,j) )
     *              / (rhoip*grav)
            work1d(m) = zhdiv / hedge(i,j)  ! c:h ratio
          else
            work1d(m) = 0.
          endif
        enddo
        write(iu,'(a,1001i5)') 'crevrat:',
     *                    ( nint(100.*work1d(m)) ,m=1,n1s,inc )

        write(iu,'(a,1001i5)') 'bedg:',
     *    ( max(99999,min(-9999, nint(1000.*bedge(i1s(m),j1s(m))) )),
     *    m=1,n1s,inc )
        write(iu,'(a,1001i5)') 'fedg:',
     *                    (nint(100.*fedge(i1s(m),j1s(m))),m=1,n1s,inc)
        do m=1,n1s
          i = i1s(m) 
          j = j1s(m) 
          if (fedge(i,j).eq.0.) then
            work1d(m) = 0.
          else
            work1d(m) = min (9999., max (1., 10000.*fedge(i,j) ))
          endif
        enddo
        write(iu,'(a,1001i5)') 'fe01:',
     *                          (nint(work1d(m)),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'hedg:',
     *                          (nint(hedge(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'basf:',
     *                 (nint(1000.*basefrml(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'quar:',
     *                (nint(1.e6*quarryrate(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'ocnm:',
     *      (min(9999,nint(1000.*oceanmelt(i1s(m),j1s(m)))),m=1,n1s,inc)

        do m=1,n1s
          work1d(m) = 0.
          do k=1,nlev
            work1d(m) = work1d(m)
     *                + heati(i1s(m),j1s(m),k)*dzeta(k)*h(i1s(m),j1s(m))
          enddo
        enddo
        write(iu,'(a,1001i5)') 'heai:',
     *                     (nint(1000.*work1d(m)/31556926), m=1,n1s,inc)

        write(iu,'(a,1001i5)') 'heab:',
     *           (nint(1000.*heatb(i1s(m),j1s(m))/31556926),m=1,n1s,inc)

        do m=1,n1s
          work1d(m) = 0.
          do k=1,nsed
            work1d(m) = work1d(m) + heats(i1s(m),j1s(m),k)
          enddo
        enddo
        write(iu,'(a,1001i5)') 'heas:',
     *                     (nint(1000.*work1d(m)/31556926), m=1,n1s,inc)

        call flowline (ua, va, work1dv, 1, 1)
        write(iu,'(a,2x,1001i5)') 'ua  :',
     *                                 (nint(work1dv(m)), m=1,n1s-1,inc)

        call flowline (uia, via, work1dv, 1, 1)
        write(iu,'(a,2x,1001i5)') 'uia :',
     *                                 (nint(work1dv(m)), m=1,n1s-1,inc)

        call flowline (ub, vb, work1dv, 1, 1)
        write(iu,'(a,2x,1001i5)') 'ub  :',
     *                                 (nint(work1dv(m)), m=1,n1s-1,inc)

        do j=0,ny
          do i=0,nx
            work2d(i,j)  = uadv(i,j)/dyu(i,j)
            work2da(i,j) = vadv(i,j)/dxv(i,j)
          enddo
        enddo
        call flowline (work2d, work2da, work1dv, 1, 1)
        write(iu,'(a,2x,1001i5)') 'uadv:',
     *  ( max(99999,min(-9999, nint(work1dv(m)*1.e-2) )), m=1,n1s-1,inc)

        call flowline (dfu, dfv, work1dv, 1, 1)
        write(iu,'(a,2x,1001i5)') 'dfu :',
     *                            (nint(work1dv(m)/100.), m=1,n1s-1,inc)

        do j=0,ny
          do i=0,nx
            work2d (i,j) = log10(crhu(i,j))
            work2da(i,j) = log10(crhv(i,j))
          enddo
        enddo
        call flowline (work2d, work2da, work1dv, 1, 0)
        write(iu,'(a,2x,1001i5)') 'crhu:',
     *                              (nint(10.*work1dv(m)),m=1,n1s-1,inc)

        call flowline (fsedu, fsedv, work1dv, 1, 0)
        write(iu,'(a,2x,1001i5)') 'fseu:',
     *                             (nint(100.*work1dv(m)),m=1,n1s-1,inc)

#if defined (MELANGE)
        write(iu,'(/a,1001i5)') 'hmel:',
     *                          (nint(hmel(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'hwme:',
     *                         (nint(hwmel(i1s(m),j1s(m))),m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'fmel:',
     *                     (nint(100.*fmel(i1s(m),j1s(m))),m=1,n1s,inc)
        do m=1,n1s
          i = i1s(m) 
          j = j1s(m) 
          if (fmel(i,j).eq.0.) then
            work1d(m) = 0.
          else
            work1d(m) = min (9999., max (1., 10000.*fmel(i,j) ))
          endif
        enddo
        write(iu,'(a,1001i5)') 'fm01:',
     *                         (nint(work1d(m)),m=1,n1s,inc)

        call flowline (uamel, vamel, work1dv, 1, 1)
        write(iu,'(a,2x,1001i5)') 'uame:',
     *                            (nint(work1dv(m)), m=1,n1s-1,inc)

        write(iu,'(a,1001i5)') 'budm:',
     *   (max(-9999,nint(1000.*budgallmel(i1s(m),j1s(m)))), m=1,n1s,inc)
        write(iu,'(a,1001i5)') 'busm:',
     *   (max(-9999,nint(1000.*budgsupmel(i1s(m),j1s(m)))), m=1,n1s,inc)
#endif

        write (iu,'(/a)') 't*10:'
        do k=0,nlevp
          write (iu,'(f5.3,1001i5)') 
     *    1.-zeta(k), (nint(10.*(t(i1s(m),j1s(m),k)-tmelt)),m=1,n1s,inc)
        enddo

        write (iu,'(a,1001i5)') 'tbme:', 
     *      (nint(10.*(tmelt-dtmdh*h(i1s(m),j1s(m))-tmelt)),m=1,n1s,inc)

        write (iu,'(/a)') 'tsed:'
        do k=1,nsed
          write (iu,'(f5.3,1001i5)') 1.-zsed(k), 
     *             (nint(10.*(tsed(i1s(m),j1s(m),k)-tmelt)),m=1,n1s,inc)
        enddo

        write (iu,'(/a)') 'wsed:'
        do k=1,nsed
          write (iu,'(f5.3,1001i5)') 
     *      1.-zsed(k), (nint(100.*(wsed(i1s(m),j1s(m),k))),m=1,n1s,inc)
        enddo

        write (iu,'(/a,1001i5)') 'tbed:',
     *             (nint(10.*(tbed(i1s(m),j1s(m),1)-tmelt)),m=1,n1s,inc)

        write (iu,'(/a)') 'u:'
        do k=1,nlev
          do j=0,ny
            do i=0,nx
              work2d(i,j)  = u(i,j,k)
              work2da(i,j) = v(i,j,k)
            enddo
          enddo
          call flowline (work2d, work2da, work1dv, 1, 1)
          write (iu,'(f5.3,2x,1001i5)') 1.-zeta(k), 
     *                                  (nint(work1dv(m)),m=1,n1s-1,inc)
        enddo

        call flush (iu)
c~~~~~~~~~~
      enddo
c~~~~~~~~~~

c=====
#endif
c=====

#if defined (LINEE) && defined (INITEE)  
      if (iflast.eq.1) then
        iu = 26  
        write (iu,
     *    "(6x,'dimension hwedge(nx+1),hswedge(nx+1),hw_wedge(nx+1)')" )
        write (iu,"( '      data hwedge /')")
        write (iu,"(('     * ',5(f11.5,',')))") h 
        write (iu,"( '     *   0. /')")

        write (iu,"('      data hswedge /')")
        write (iu,"(('     * ',5(f11.5,',')))") hs
        write (iu,"( '     *   0. /')")

        write (iu,"('      data hw_wedge/')")
        write (iu,"(('     * ',5(f11.5,',')))") hw
        write (iu,"( '     *   0. /')")
      endif
#endif

      return
      end

c-----------------------------------------------------------------------

#if defined (EISANTA) && ! defined (TRANSECTA)

      subroutine setflowline 

c     For 2D Antarctic runs, sets coords for 1D-flowline in iceshow1d: 
c     n1s, i1s, j1s, d1s (in comicegrid, returned). 
c     a[lon,lat][a,b] are lons,lats of 2 end points (supplied).

#include <comicephys.h>
#include <comicegrid.h>

#if defined (FORCE34MYR)
      parameter (nseg=1, nsegp=nseg+1)
      dimension alonvert(nsegp), alatvert(nsegp)
      data  alonvert(1),  alatvert(1),  alonvert(2), alatvert(2)
c       from Weddell to Ross Seas, ~through Siple Dome:
     *   /
     *       -40.,   -70.,   178.,   -70.
     *   /
      data c1s /'Weddell to Ross'/

#elif (MURPHYLINE)
c     through Mt. Murphy (-75.387, -111.262) northwards:
      parameter (nseg=2, nsegp=nseg+1)
      dimension alonvert(nsegp), alatvert(nsegp)
      data alatvert /
     *    -77.0,   -70.0,   -68.0 / 
      data alonvert /
     *   -111.26, -111.26, -111.26 /
      data c1s /'Mt. Murphy'/

#elif (HUDSONLINE)
c     through Hudson Mts (Maish Nunatack -74.6,-99.5) -> PIG trough:
      parameter (nseg=10, nsegp=nseg+1)
      dimension alonvert(nsegp), alatvert(nsegp)
      data alatvert /
     *           -74.5,  -74.6,  -74.75, -74.5,  -74.0,
     *   -73.6,  -72.8,  -72.4,  -72.0,  -70.8,  -68.0 /
      data alonvert /
     *           -98.0,  -99.5, -102.5, -103.8, -105.5,
     *  -106.5, -107.2, -107.0, -106.5, -104.7, -101.2 /
      data c1s /'Hudson Mts.'/

#else
c     down PIG centerline -> PIG trough:
      parameter (nseg=11, nsegp=nseg+1)
      dimension alonvert(nsegp), alatvert(nsegp)
      data alatvert /
     *   -75.2,  -75.3,  -75.3,  -74.9,  -74.5,  -74.0,
     *   -73.6,  -72.8,  -72.4,  -72.0,  -70.8,  -68.0 /
      data alonvert /
     *   -80.0,  -85.0,  -97.5, -102.0, -103.8, -105.5,
     *  -106.5, -107.2, -107.0, -106.5, -104.7, -101.2 /
      data c1s /'PIG trough'/

#endif
c=====
c     single line between 2 points:
c     parameter (nseg=1, nsegp=nseg+1)
c     dimension alonvert(nsegp), alatvert(nsegp)
c     data  alonvert(1),  alatvert(1),  alonvert(2), alatvert(2)
c    *   /
cc      from Weddell to Ross Seas, ~through Siple Dome:
cc   *       -40.,   -70.,   178.,   -70.
cc   *      -110.,  -83.0,   170.,   -70. ! (?)
cc      across Mac.Robertson Land shelf (A. Mackintosh):
cc   *        65.,   -70.,    65.,   -65.
cc   *        63.,   -70.,    63.,   -65.
cc      across EAIS and WAIS:
cc   *      -102.,   -70.,    78.,   -65.
c    *        20.,   -68.,   140.,   -65.
c    *  /
c=====
c     like in subr settransect (iceinit.F), for 2-D with TRANSECTA: 
c     parameter (nseg=1,  nsegp=nseg+1)
c     parameter (nseg=2,  nsegp=nseg+1)
c     parameter (nseg=3,  nsegp=nseg+1)
c     parameter (nseg=4,  nsegp=nseg+1)
c     parameter (nseg=11, nsegp=nseg+1)
c     dimension alonvert(nsegp), alatvert(nsegp)
c+++++
c     data alatvert /-65., -82.5,  -82.5,  -65./          ! weddell+ross
c     data alonvert /-30., -100.,  -170., -180./
c+++++
c     data alatvert /-70., -82.5,  -82.5,  -68./          ! weddell+ross
c     data alonvert /-35., -100.,  -170., -179./
c+++++
c     data alatvert /-70., -81.5,  -82.5,  -82.5,  -68./  ! weddell+ross
c     data alonvert /-54.,  -57.,  -100.,  -170., -179./
c+++++
c     data alatvert /-82.5,  -82.5,  -68./                ! ross only
c     data alonvert /-100.,  -170., -179./
c+++++
c     data alatvert /-65., -90.,  -65./
c     data alonvert /-90., -90.,   100./
c+++++
c     data alatvert /-73.82, -64.08 /                     ! ross-> prydz
c     data alonvert /177.40,  67.75 /
c+++++
c     data alatvert /-77.5, -75.7, -74.3, -70.0, -65.0/   ! pig -> ocean
c     data alonvert /-92.0, -98.0,-103.0,-107.0,-107.0/
c+++++
c     data alatvert /                                     ! pig trough
c    *   -75.2,  -75.3,  -75.3,  -74.9,  -74.5,  -74.0,
c    *   -73.6,  -72.8,  -72.4,  -72.0,  -70.8,  -68.0 /
c     data alonvert /
c    *   -80.0,  -85.0,  -97.5, -102.0, -103.8, -105.5,
c    *  -106.5, -107.2, -107.0, -106.5, -104.7, -101.2 /
c+++++
c     data c1s /'...'/
c=====

      dimension xvert(nsegp), yvert(nsegp), seglen(nseg)

c        Loop over each straight-line segment in turn.
c        Only include 1st point in segment (the vertex) for 1st segment,
c        and alwats include the last point (the 2nd vertex).

      n1s = 0
c>>>>>>>>>>>>>>>>>>>
      do iseg=1,nseg
c>>>>>>>>>>>>>>>>>>>

      alona = alonvert(iseg)
      alonb = alonvert(iseg+1)
      alata = alatvert(iseg)
      alatb = alatvert(iseg+1)

c        Convert lon,lat of 2 end pts to i,j polar stereo coords.

      call ps_xy (zxa, zya, alata*pi/180, alona*pi/180., xoffa, yoffa)
      ia = (zxa + 0.5*nx*dx0)/dx0 + 1.0001
      ja = (zya + 0.5*ny*dy0)/dy0 + 1.0001

      call ps_xy (zxb, zyb, alatb*pi/180, alonb*pi/180., xoffa, yoffa)
      ib = (zxb + 0.5*nx*dx0)/dx0 + 1.0001
      jb = (zyb + 0.5*ny*dy0)/dy0 + 1.0001

      ia = max (1, min (nx, ia))
      ja = max (1, min (ny, ja))
      ib = max (1, min (nx, ib))
      jb = max (1, min (ny, jb))

      if (ia.eq.ja .and. ib.eq.jb) then

c        Do nothing: segment is too short, within one grid cell 

c        Set transect coords (i1s,j1s). If transect is more in x 
c        direction, increment i uniformly by 1 and calculate nearest j,
c        or v.v if transect is more in y direction.

      else if (iabs(ia-ib).ge.iabs(ja-jb)) then

        npoi = iabs(ib-ia) + 1
        mpoi = 0
        do i = ia,ib,isign(1,ib-ia)
          mpoi = mpoi + 1
c         include 1st pt for 1st segment only:
          if (iseg.eq.1 .or. i.ne.ia) then
            n1s = n1s + 1
            i1s(n1s) = i
            j1s(n1s) = ja + nint((mpoi-1.)*(jb-ja)/max(npoi-1.,1.))
          endif
        enddo

      else

        npoi = iabs(jb-ja) + 1
        mpoi = 0
        do j = ja,jb,isign(1,jb-ja)
          mpoi = mpoi + 1
c         include 1st pt for 1st segment only:
          if (iseg.eq.1 .or. j.ne.ja) then
            n1s = n1s + 1
            j1s(n1s) = j
            i1s(n1s) = ia + nint((mpoi-1.)*(ib-ia)/max(npoi-1.,1.))
          endif
        enddo

      endif

c>>>>>>>>>>
      enddo
c>>>>>>>>>>

      if (n1s.gt.nx+ny) then
        write (ioterm, '(/a/a,i8,a,i8)') 
     *        '*** Error (setflowline): 1-D transect is too long.',
     *        '    n1s=',n1s, '   max (nx+ny) =',nx+ny
        stop
      endif

c       Set distance (m) along transect 

      d1s(1) = 0.5*dd0
      do m = 2,n1s
        d1s(m) = d1s(m-1)
     *         + sqrt (  (xh(i1s(m),j1s(m)) - xh(i1s(m-1),j1s(m-1)))**2
     *                 + (yh(i1s(m),j1s(m)) - yh(i1s(m-1),j1s(m-1)))**2)
      enddo

c        Diagnostic output

      iu = 89
      if (iu.ne.0) then
        write (iu,'(a,i6,/2(a,2f8.2/),a)') 
     *    'setflowline:  n1s=', n1s,
     *    'alona, alata = ', alona, alata,
     *    'alonb, alatb = ', alonb, alatb,
     *    '   m   alond   alatd       d1s'

        do m=1,n1s
          write (iu,'(i4,2f8.2, f10.1)')
     *      m, alond(i1s(m),j1s(m)), alatd(i1s(m),j1s(m)), d1s(m)*1.e-3
        enddo
        call flush (iu)
      endif

      return
      end

#endif

c-----------------------------------------------------------------------

      subroutine flowline (uarr,   varr,  work, ifreal, ifvec)
      entry     iflowline (iuarr, ivarr, iwork, ifreal, ifvec)

#include <comicegrid.h>

c        Assembles 1-D flowline for a [u,v]-grid quantity,
c        based on h-grid indices i1s,j1s,n1s (in comicegrid).
c        ifreal = 1 for real (vector/scalar), else int (always scalar).
c        ifvec  = 1 for vector (for real only).
c        Distinguish between u or v direction by j1s changing or not.  
c
c        Flowline is defined by i1s,j1s,n1s (in comicegrid, set by 
c        setflowline/iceshow1d for 2-D, or within iceshow1d for 1-D).
c        Called by iceshow1d and iceplot1d.

      dimension 
     *  uarr (0:nxp,0:nyp), varr (0:nxp,0:nyp),  
     *  iuarr(0:nxp,0:nyp), ivarr(0:nxp,0:nyp),  
     *  work(0:nx+ny+1),    iwork(0:nx+ny+1)  

      do m=0,n1s-1

        if (m.eq.0) then
c         first (zero index) point possibly on domain boundary 
c         (erroneous if i1s or j1s decrease away from nx,ny boundary)
          iu = i1s(1)-1
          ju = j1s(1)
          iv = i1s(1)
          jv = j1s(1)-1
          mm = 1
        else
          iu = min (i1s(m),i1s(m+1))
          ju = j1s(m)
          iv = i1s(m)
          jv = min (j1s(m),j1s(m+1))
          mm = m
        endif
        
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (ifreal.eq.1) then
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          if (ifvec.eq.1) then
            work(m) = sqrt (uarr(iu,ju)**2 + varr(iv,jv)**2)
#if defined (TESTMEL2D)
            work(m) = sign ( work(m), uarr(iu,ju) )        ! for nx > ny
#else
c           positive if vector is -45 to +135 deg relative to x axis:
            work(m) = sign ( work(m), 
     *                  cos(atan2(varr(iv,jv),uarr(iu,ju)) - 0.5*pi) )
#endif
          else
            if (j1s(mm).eq.j1s(mm+1)) then
              work(m) = uarr(iu,ju)                        ! u direction
            else
              work(m) = varr(iv,jv)                        ! v direction
            endif
          endif
c~~~~~~~~~~~
        else
c~~~~~~~~~~~
          if (j1s(mm).eq.j1s(mm+1)) then 
            iwork(m) = iuarr(iu,ju)                        ! u direction
          else
            iwork(m) = ivarr(iv,jv)                        ! v direction
          endif
c~~~~~~~~~~~~
        endif
c~~~~~~~~~~~~

      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine iceshow2d (h, hs, hb, hbsd, t, tracer, 
     *                    fedge, hedge,
     *                    budgall, budgsnow, budgrain, budgmelt,
     *                    budgevap, budgneg,
     *                    baseperc, basefrml, oceanmelt, msector,
     *                    facemelt, cliffmelt,
#if defined (RACMOCLIM)
     *                    runofa,
#endif
     *                    div, calvice, arcocn,
     *                    tsurf, tsurfi, tmon,
     *                    heati, heath, heatb, arhap,
     *                    w, wa, sedim,  
     *                    quarryrate, quarryacc, quarrycoef,
     *                    pelagrate, pelagacc,
     *                    topbed, topbedeq, deflect, geoflux,
     *                    hw, tw, maskwater, crhmel, enhanceinv,
     *                    hsbal, hbinit, tempocn,
     *                    u, v, ua, va, ub, vb,
     *                    uadv, vadv, hu, hv, dfu, dfv, dbu, dbv,
     *                    thetau, thetav, hgu, hgv,
     *                    fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                    angnormu, angnormv,
     *                    crhu, crhv,
#if defined (MELANGE)
     *                    hmel, hwmel, fmel, maskmel, 
     *                    budgallmel, budgsupmel, budgadjmel,budgsmomel,
     *                    oceanmeltmel,
     *                    uamel, vamel, 
     *                    bpumel, bpvmel, 
     *                    thetaumel, thetavmel,
     *                    uaboundmel, vaboundmel,
#endif
     *                    numh, nuvtot, itera, iterc,
     *                    toth0, tota0, toth, tota, totflow,
     *                    sealev, timeice, dt, weirun,
     *                    iffirse, iflast, nyearout2d)

c     Ascii 2-D printout - x,y maps

#include <comicephys.h>
#include <comicegrid.h>
#include <comiceplot.h>

      dimension 
     *  h(nx,ny),            hs(nx,ny),          
     *  hb(nx,ny),           hbsd(nx,ny),
     *  t(nx,ny,0:nlevp),    tracer(nx,ny,0:nlevp,ntrace),
     *  fedge(nx,ny),        hedge(nx,ny),
     *  budgall(nx,ny),      budgsnow(nx,ny),   
     *  budgrain(nx,ny),     budgmelt(nx,ny),   
     *  budgevap(nx,ny),     budgneg(nx,ny),
     *  baseperc(nx,ny),     basefrml(nx,ny),     
     *  oceanmelt(nx,ny),    msector(nx,ny),    
     *  facemelt(nx,ny),     cliffmelt(nx,ny),   
#if defined (RACMOCLIM)
     *  runofa(nx,ny),
#endif
     *  div(nx,ny),          calvice(nx,ny),    arcocn(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),     tmon(nmon,nx,ny),
     *  heati(nx,ny,nlev),   heath(nx,ny),      
     *  heatb(nx,ny),        arhap(nx,ny),
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp), sedim(nx,ny),    
     *  quarryrate(nx,ny),   quarryacc(nx,ny),   quarrycoef(nx,ny), 
     *  pelagrate(nx,ny),    pelagacc(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny),
     *  deflect(nx,ny),      geoflux(nx,ny),
     *  hw(nx,ny),           tw(nx,ny),         maskwater(nx,ny),
     *  crhmel(nx,ny),       enhanceinv(nx,ny),
     *  hsbal(nx,ny),        hbinit(nx,ny),     tempocn(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  dfu(0:nxp,0:nyp),       dfv(0:nxp,0:nyp),
     *  dbu(0:nxp,0:nyp),       dbv(0:nxp,0:nyp),
     *  thetau(nx,ny),          thetav(nx,ny),
     *  hgu(0:nxp,0:nyp),       hgv(0:nxp,0:nyp),
     *  fluxgrdu(0:nxp,0:nyp),  fluxgrdv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp),  fluxschv(0:nxp,0:nyp),
     *  angnormu(0:nxp,0:nyp),  angnormv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp)
#if defined (MELANGE)
     * ,hmel(nx,ny),            hwmel(nx,ny),
     *  fmel(nx,ny),            maskmel(nx,ny),     
     *  budgallmel(nx,ny),      budgsupmel(nx,ny),  
     *  budgadjmel(nx,ny),      budgsmomel(nx,ny),
     *  oceanmeltmel(nx,ny),
     *  uamel(0:nxp,0:nyp),     vamel(0:nxp,0:nyp),
     *  bpumel(0:nxp,0:nyp),    bpvmel(0:nxp,0:nyp),
     *  thetaumel(nx,ny),       thetavmel(nx,ny),
     *  uaboundmel(0:nxp,0:nyp), vaboundmel(0:nxp,0:nyp)
#endif

      dimension 
     *  work(nx,ny), work2(nx,ny), work3(nx,ny)
      logical firstshow2
      data firstshow2 /.true./
      save firstshow2

#if defined (RACMOCLIM)   
      parameter (nbrac=10)
      dimension racz(nbrac), racl(nbrac), 
     *          rmod(nbrac,nbrac), rrac(nbrac,nbrac)
      character*10 creg(5)
#endif

c     Statement function (fpin) for fraction of u,v grid box with 
c     sub-grid bedrock mounds scraping ice bottom:
#include <comicepin.h>

c       Skip out if not time

      if (nyearout2d.eq.0) return
      if ( .not.
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyearout2d),dt)).lt.dt
     *       .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyearout2d.eq.-1 )
     *   ) return 

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      nprint = 1  ! 1, 2
      do iprint=1,nprint
       if (iprint.eq.1) iu = iuout2d
       if (iprint.eq.2) iu = 6
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

c          Write 2-D ascii maps

        if (firstshow2) then
          call printmap (timeice,alatd,   ' lat (deg S)',    5., iu, 0)
          call printmap (timeice,alond,   ' lon (deg E)',   10., iu, 0)
          call printmap (timeice,topbedeq,'equil bdrk (m)',150., iu, 0)

          call zero (work, nx*ny)
          do m=1,n1s
            work(i1s(m),j1s(m)) = 1.
          enddo
          call printmap (timeice,work,'flowline',1., iu,0)

#if defined (FORCE34MYR) && defined (SEDIMENT)
          call itor (msedsec, work, nx*ny)
          call printmap (timeice,work,'Wilson basins 2', 1., iu, 0)
#endif

          if (iprint.eq.nprint) firstshow2 = .false.
        endif
  
        call itor (maskwater, work, nx*ny)
        call printmap (timeice,work,    'maskwater',      1., iu, 0)
        call printmap (timeice,hw,      'hw(m)         ',40., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = tw(i,j) - tmelt
            if (maskwater(i,j).eq.0) work(i,j) = 0.
          enddo
        enddo
        call printmap (timeice,work,     'tw(C)        ', 0.1, iu, 0)

        call printmap (timeice,baseperc, 'baseperc(m/y)', .01, iu, 0)
        call printmap (timeice,baseperc, 'baseperc(m/y)', 0.1, iu, 0)
        call printmap (timeice,basefrml, 'basefrml(m/y)', .01, iu, 0)
        call printmap (timeice,basefrml, 'basefrml(m/y)', 0.1, iu, 0)
        call printmap (timeice,oceanmelt,'oceanmelt(m/y)', 10.,iu, 0)
        call printmap (timeice,oceanmelt,'oceanmelt(m/y)', 1., iu, 0)
        call printmap (timeice,oceanmelt,'oceanmelt(m/y)',0.05,iu, 0)
#if defined (OCMARTIN)
        call itor (msector, work, nx*ny)
        do j=1,ny
          do i=1,nx
            if (msector(i,j).eq.0) then
              work(i,j) = 0.
            else if (msector(i,j).ge.1 .and. msector(i,j).le.26) then
              work(i,j) = msector(i,j)
            else if (msector(i,j).ge.27) then
              work(i,j) = -msector(i,j) + 26
            else
              work(i,j) = -1.e20
            endif
          enddo
        enddo
        call printmap (timeice, work, 'msector', 1., iu, 0)
#endif

        call printmap (timeice,facemelt, 'facemelt (m/y)',  20., iu, 0)
        call printmap (timeice,cliffmelt, 'cliffmelt (m/y)',10., iu, 0)
        call printmap (timeice,cliffmelt, 'cliffmelt (m/y)', 1., iu, 0)
        call printmap (timeice,div,      'div  (1/y)',     .001, iu, 0)
        call printmap (timeice,calvice,  'calvice(m/y)',    10., iu, 0)
        call printmap (timeice,calvice,  'calvice(m/y)',     1., iu, 0)

c       sfc+basal crevasse depth:ice thickness (Nick, as subr calving):
        do j=1,ny
          do i=1,nx
            if (maskwater(i,j).eq.1 .and. hedge(i,j).gt.0.) then
              zhdiv = 2. * ( (max(div(i,j),0.)**powir)/arhap(i,j) )
     *                / (rhoip*grav)
              work(i,j) = zhdiv / hedge(i,j)  ! c:h ratio
            else
              work(i,j) = -1.
            endif
          enddo
        enddo
        call printmap (timeice,work,  'crevrat',  0.1, iu, 0)

#if defined (CALVDAMAGE) 
        call printmap (timeice, tracer(1,1,1,2), 
     *                                     'calvdam r',        .05,iu,1)
        call printmap (timeice, calvdam_t1,'calvdam term1',    .01,iu,1)
        call printmap (timeice, calvdam_t2,'calvdam term2',    .01,iu,1)
        call printmap (timeice, calvdam_t3,'calvdam term3',    .01,iu,1)
#elif defined (CALVBALANCE) 
        call printmap (timeice, calvbal_vi,'calvbal vi (m/y)',100.,iu,1)
        call printmap (timeice, calvbal_vb,'calvbal vb (m/y)',100.,iu,1)
        call printmap (timeice, calvbal_vc,'calvbal vc (m/y)',100.,iu,1)
#endif

        call printmap (timeice,arcocn,   'arcocn (deg)',  10., iu, 0)
        call printmap (timeice,budgall,  'budgall (m/y)', 0.1, iu, 0)
        call printmap (timeice,budgsnow, 'budgsnow(m/y)', 0.1, iu, 0)
c       call printmap (timeice,budgsnow, 'budgsnow(m/y)', .01, iu, 0)
        call printmap (timeice,budgrain, 'budgrain(m/y)', 0.1, iu, 0)
        call printmap (timeice,budgmelt, 'budgmelt(m/y)', 0.1, iu, 0)

c>>>>>>>>>>>>>>>>>>>>>>
#if defined (RACMOCLIM)
c>>>>>>>>>>>>>>>>>>>>>>

c          Output to compare budgmelt with RACMO runoff for Greenland

        iuracmo  = 351
        iuracmo2 = 352
c       hracmin = 150.              ! only include (non-thin) ice points
        hracmin = 170.              ! only include (non-thin) ice points

c          iuracmo output

        call printmap (timeice,h,  'ice thickness (m)', 150., iuracmo,0)

        call zero (work, nx*ny)
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.hracmin)
     *        work(i,j) = budgmelt(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'budgmelt(m/y)', 0.1,  iuracmo, 0)

        call zero (work, nx*ny)
        do j=1,ny
          do i=1,nx
c           mm/day liq. equiv to m/yr ice equiv:
            if (h(i,j).gt.hracmin)
     *        work(i,j) = runofa(i,j) * .001*365.*(rholiq/rhoice)
          enddo
        enddo
        call printmap (timeice, work,'racmo runoff (m/y)',0.1,iuracmo,0)

        call zero (work2, nx*ny)
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.hracmin)
     *        work2(i,j) = budgmelt(i,j) - work(i,j)  
          enddo
        enddo
        call printmap (timeice, work2, 'budgmelt-racmo runoff (m/y)',
     *                 0.1, iuracmo, 0)

        call zero (work2, nx*ny)
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.hracmin) 
     *        work2(i,j) = budgmelt(i,j) - work(i,j)  
          enddo
        enddo
        call printmap (timeice, work2, 'budgmelt-racmo runoff (m/y)',
     *                 0.1, iuracmo, 0)

        call zero (work2, nx*ny)
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.hracmin) then
              work2(i,j) = -1.e20
              do m=1,nmon
                work2(i,j) = max (work2(i,j), tmon(m,i,j)-tmelt)
              enddo
            endif
          enddo
        enddo
        call printmap (timeice,work2,'tamax (C)', 2., iuracmo, 0)

c>>>>>

c          iuracmo2 output

        write (iuracmo2,'(/80("=")/a,f12.2/)') 'time=',timeice

c          Compute overall avg budgmelt and RACMO runoff (zmtot,zrtot),
c          rms diff (zrms), and zero/non-zero mismatch fraction (zmis)

        zmtot = 0.
        zrtot = 0.
        zrms = 0.
        zmis = 0.
        zar  = 0.
        zmm = .001          ! small value below which melt ~ zero (m/yr)
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.hracmin) then
              zmtot = zmtot + budgmelt(i,j)*darea(i,j)
              zrtot = zrtot + work(i,j)    *darea(i,j)

              zrms = zrms + ((budgmelt(i,j)-work(i,j))**2) *darea(i,j)

              if ( (budgmelt(i,j).le.zmm .and. work(i,j).gt.zmm) .or.
     *             (budgmelt(i,j).gt.zmm .and. work(i,j).le.zmm) ) then
                zmis = zmis + darea(i,j)
              endif
              zar = zar + darea(i,j)
            endif
          enddo
        enddo
        zmtot = zmtot / max(zar,.001)
        zrtot = zrtot / max(zar,.001)
        zrms = sqrt (zrms / max(zar,.001))
        zmis = zmis / max(zar,.001)
        write (iuracmo2,'(a,3(f10.4,1x))') 
     *   'budgmelt_RACMOrunoff_diff(m_y-1)=   ', zmtot,zrtot,zmtot-zrtot
        write (iuracmo2,'(a,f10.4)') 
     *   'rms_diff_budgmelt_vs_racmo_(m_y-1)= ', zrms
        write (iuracmo2,'(a,f10.4)') 
     *   'mismatch_fraction_(0_vs_non-0)=     ', zmis

        zdifnorm = max (0., min  (1., 1. - abs(zmtot-zrtot)/2. ))
        zrmsnorm = max (0., min  (1., 1. - zrms/2. ))
        zmisnorm = max (0., min  (1., 1. - zmis/.5 ))

        write (iuracmo2,'(a,f10.4)') 
     *   'score_(0-1_perfect_is_1)=           ',
     *    (zdifnorm**0.2) * (zrmsnorm**0.2) * (zmisnorm**0.6)

c          Compute averages of budgmelt and RACMO runoff
c          for 2-D bins of elevation and latitude (nbrac**2),
c          Write 2-D tables to iuracmo2.
c          (From above, work = runoff, work2 = tamax-tmelt).

c---------------------
        do ibz=1,nbrac
c---------------------
          racza = 0. + 300.*(ibz-1)
          raczb = racza + 300.
          racz(ibz) = 0.5*(racza+raczb)
          if (ibz.eq.1) then
            racza = -1000.
          else if (ibz.eq.nbrac) then
            raczb = 10000.
          endif
c---------------------
        do ibl=1,nbrac
c---------------------
          racla = 60. + 2.5*(ibl-1)
          raclb = racla + 2.5
          racl(ibl) = 0.5*(racla+raclb)
          if (ibl.eq.1) then
            racla = 0.
          else if (ibl.eq.nbrac) then
            raclb = 90.
          endif

          rmod(ibl,ibz) = 0.
          rrac(ibl,ibz) = 0.
          zar = 0.
          do j=1,ny
            do i=1,nx
              if (h(i,j).gt.hracmin) then
                if (hs(i,j).gt.racza .and. hs(i,j).le.raczb .and.
     *              alatd(i,j).gt.racla .and. alatd(i,j).le.raclb) then
                  rmod(ibl,ibz)=rmod(ibl,ibz) + budgmelt(i,j)*darea(i,j)
                  rrac(ibl,ibz)=rrac(ibl,ibz) + work(i,j)    *darea(i,j)
                  zar = zar + darea(i,j)
                endif
              endif
            enddo
          enddo
          rmod(ibl,ibz) = rmod(ibl,ibz) / max (.001,zar)
          rrac(ibl,ibz) = rrac(ibl,ibz) / max (.001,zar)
c------------
        enddo
        enddo
c------------

        write (iuracmo2,'(/a/8x,10f7.2)')
     *    'model budgmelt bins, lat vs elev, (m/y):', 
     *    (racl(ibl), ibl=1,nbrac)
        do ibz=nbrac,1,-1
            write  (iuracmo2,'(f6.1, 2x,10f7.3)') 
     *        racz(ibz), (rmod(ibl,ibz),ibl=1,nbrac)
        enddo

        write (iuracmo2,'(/a/8x,10f7.2)')
     *    'racmo runoff bins, lat vs elev, (m/y):', 
     *    (racl(ibl), ibl=1,nbrac)
        do ibz=nbrac,1,-1
            write  (iuracmo2,'(f6.1, 2x,10f7.3)') 
     *        racz(ibz), (rrac(ibl,ibz),ibl=1,nbrac)
        enddo

        write (iuracmo2,'(/a/8x,10f7.2)')
     *    'budgmelt-racmo runoff bins, lat vs elev, (m/y):', 
     *    (racl(ibl), ibl=1,nbrac)
        do ibz=nbrac,1,-1
            write  (iuracmo2,'(f6.1, 2x,10f7.3)') 
     *        racz(ibz), (rmod(ibl,ibz)-rrac(ibl,ibz),ibl=1,nbrac)
        enddo

c          Compute linear regression coefficients a,b, for
c          RACMO runoff (or budgmelt) = a * (tamax-tmelt) + b.

c~~~~~~~~~~~~~~~~~~~~
        do ifield=1,2
c~~~~~~~~~~~~~~~~~~~~
c       from above, work = runoff, work2 = tamax-tmelt:
        if (ifield.eq.1) call  scopy (nx*ny, work, 1, work3, 1)
        if (ifield.eq.2) call  scopy (nx*ny, budgmelt, 1, work3, 1)

c....................
        do ireg = 1,5
c....................
        if (ireg.eq.1) then
          creg(ireg) = 'Southwest'
          zlat1 =   0. 
          zlat2 =  72. 
          zlon1 = -90. 
          zlon2 = -45. 
        else if (ireg.eq.2) then
          creg(ireg) = 'Southeast'
          zlat1 =   0. 
          zlat2 =  72. 
          zlon1 = -45. 
          zlon2 =   0. 
        else  if (ireg.eq.3) then
          creg(ireg) = 'Northwest'
          zlat1 =  72. 
          zlat2 =  90. 
          zlon1 = -90. 
          zlon2 = -45. 
        else if (ireg.eq.4) then
          creg(ireg) = 'Northeast'
          zlat1 =  72. 
          zlat2 =  90. 
          zlon1 = -45. 
          zlon2 =   0. 
        else if (ireg.eq.5) then
          creg(ireg) = 'All'
          zlat1 =   0. 
          zlat2 =  90. 
          zlon1 = -90. 
          zlon2 =   0. 
        endif
       
        zsx2 = 0.
        zsxy = 0.
        zsx  = 0.
        zsy  = 0.
        zsn  = 0.
        do j=1,ny
          do i=1,nx
            if (alatd(i,j).gt.zlat1 .and.  alatd(i,j).lt.zlat2 .and.  
     *          alond(i,j).gt.zlon1 .and.  alond(i,j).lt.zlon2 .and.  
     *          h(i,j).gt.hracmin .and. work3(i,j).gt.0.)  then
              zsx2 = zsx2 + work2(i,j)**2
              zsxy = zsxy + work2(i,j)*work3(i,j)
              zsx  = zsx  + work2(i,j)
              zsy  = zsy  + work3(i,j)
              zsn  = zsn  + 1.
            endif
          enddo
        enddo
        zdet = zsx2*zsn - zsx*zsx
        za =  (zsn*zsxy - zsx *zsy) / zdet
        zb = (-zsx*zsxy + zsx2*zsy) / zdet
       
c       alternate "covarience" calc of a (gives same result):
c       zt  = 0.
c       zt2 = 0.
c       zr  = 0.
c       ztr = 0.
c       n = 0
c       do j=1,ny
c         do i=1,nx
c           if (alatd(i,j).gt.zlat1 .and.  alatd(i,j).lt.zlat2 .and.  
c    *          alond(i,j).gt.zlon1 .and.  alond(i,j).lt.zlon2 .and.  
c    *          h(i,j).gt.hracmin .and. work3(i,j).gt.0.)  then
c             zt  = zt  + work2(i,j)
c             zt2 = zt2 + work2(i,j)**2
c             zr  = zr  + work3(i,j)
c             ztr = ztr + work2(i,j)*work3(i,j)
c             n = n + 1
c           endif
c         enddo
c       enddo
c       za = (ztr-zt*zr) / (zt2-zt*zt)   ! factors of 1/n**2 cancel

        if (ireg.eq.1) then
          if (ifield.eq.1) write (iuracmo2,'(/a)')
     *      'regression: RACMO runoff = a*(tamax-tmelt) + b:'
          if (ifield.eq.2) write (iuracmo2,'(a)')
     *      'regression: model budgmelt = a*(tamax-tmelt) + b:'
        endif
        write (iuracmo2,'(2x, 2a,2f10.4)') 
     *    creg(ireg), '  a,b=', za,zb
c............
        enddo  ! ireg
c............
c~~~~~~~~~~~~
        enddo  ! ifield
c~~~~~~~~~~~~

c>>>>>
#endif
c>>>>>

        call printmap (timeice,budgevap, 'budgevap(m/y)', 0.1, iu, 0)
        call printmap (timeice,budgneg,  'budgneg(m/y)',  0.1, iu, 0)
        call printmap (timeice,hb,'bed+sed elev (m)',    150., iu, 0)
        call printmap (timeice,hbsd,'bed+sed stddev (m)', 50., iu, 0)
        call printmap (timeice,topbed,'topbed (m)',      150., iu, 0)
#if defined (SEDIMENT)
        call printmap (timeice,sedim,'sedim (m)',         .05, iu, 0)
        call printmap (timeice,sedim,'sedim (m)',         10., iu, 0)
        call printmap (timeice,sedim,'sedim (m)',        150., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = log10(quarrycoef(i,j))
          enddo
        enddo
c       m/y per [tau(N/m2)*ub(m/y)]
        call printmap (timeice,work,'log10(quarrycoef)', 1., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = quarryrate(i,j) * 1.e6
          enddo
        enddo
        call printmap (timeice,work,'quarryrate (m/Myr)',  1., iu, 0)

        call printmap (timeice,quarryacc,'quarryacc (m)', 10., iu, 0)
        call printmap (timeice,quarryacc,'quarryacc (m)',150., iu, 0)
c       call printmap (timeice,pelagacc,'pelagacc (m)',  150., iu, 0)
#endif
        call printmap (timeice,deflect,  'deflect (m)',   10., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = geoflux(i,j) / (.001*31556926)
          enddo
        enddo
        call printmap (timeice, work, 'geoflux (mW/m2) (m)', 3., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j)  = 0.
            work2(i,j) = -1000.
            do m=1,nmon
              work(i,j) = work(i,j) + (tmon(m,i,j)-tmelt)/nmon
              work2(i,j) = max (work2(i,j), tmon(m,i,j)-tmelt)
            enddo
          enddo
        enddo
        call printmap (timeice,work,'ta (C)', 2., iu, 0)
        call printmap (timeice,work2,'tamax (C)', 2., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = tsurf(i,j) - tmelt
          enddo
        enddo
        call printmap (timeice,work,'tsurf (C)', 2., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = t(i,j,0) - tmelt
          enddo
        enddo
        call printmap (timeice,work,'ice surf t (C)', 2., iu, 0)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = t(i,j,nlevp) - (tmelt-dtmdh*h(i,j))
c           work(i,j) = t(i,j,nlevp) -  tmelt
c           if (h(i,j).eq.0.) work(i,j) = 0. 
          enddo
        enddo
        call printmap (timeice,work, 'basal t (C)',  1., iu, 0)
        call printmap (timeice,work, 'basal t (C)',  .1, iu, 0)
  
        do j=1,ny
          do i=1,nx
            call calc_crh (zcrh, work(i,j),
     *                     crhmel(i,j), 
     *                     t(i,j,nlevp) - (tmelt-dtmdh*h(i,j)),
     *                     baseperc(i,j)-basefrml(i,j),
     *                     h(i,j), 
     *                     hbsd(i,j), 
     *                     1)
          enddo
        enddo
        call printmap (timeice,work,'basal non-froz/lub. fraction',
     *                 .1, iu, 0)

        call calcz0 (work, h, t)
        call printmap (timeice, work, 'thickwarm', 150., iu, 0)

        call zero (work, nx*ny)
        do j=1,ny
          do i=1,nx-1
            zhw  = 0.5*(hw(i,j) + hw(i+1,j))
            zhbsd  = 0.5*(hbsd(i,j) + hbsd(i+1,j))
            work(i,j) = fpin (zhw, zhbsd)
          enddo
        enddo
        call printmap (timeice, work, 'fpinu', .1, iu, 1)

        do j=1,ny
          do i=1,nx
            work(i,j) = log10(crhmel(i,j))
          enddo
        enddo
        call printmap (timeice, work, 'log10(crhmel)', 1., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = 0.25*log10( crhu(i,j)*crhu(i-1,j) 
     *                             *crhv(i,j)*crhv(i,j-1))
          enddo
        enddo
        call printmap (timeice, work, 'log10(crhuv)', 1., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = log10(crhu(i,j))
          enddo
        enddo
        call printmap (timeice, work, 'log10(crhu)', 1., iu, 1)

        do j=1,ny
          do i=1,nx
            work(i,j) = log10(crhv(i,j))
          enddo
        enddo
        call printmap (timeice, work, 'log10(crhv)', 1., iu, 2)

        do j=1,ny
          do i=1,nx
            work(i,j) = log10(enhanceinv(i,j))
          enddo
        enddo
        call printmap (timeice, work, 'log10(enhanceinv)', .1, iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = ua(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'ua',      10.0, iu, 1)
        call printmap (timeice,work, 'ua_100', 100.0, iu, 1)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = va(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'va',  10.0, iu, 2)
        call printmap (timeice,work, 'va_100',  100.0, iu, 2)

c---------------------
#if defined (SCHOOFGL)
c---------------------
        do j=1,ny
          do i=1,nx
            work(i,j) = fluxgrdu(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'fluxgrdu', 1.e3, iu, -1)

        do j=1,ny
          do i=1,nx
            work(i,j) = fluxgrdv(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'fluxgrdv', 1.e3, iu, -2)

        do j=1,ny
          do i=1,nx
            work(i,j) = fluxschu(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'fluxschu', 1.e3, iu, -1)

        do j=1,ny
          do i=1,nx
            work(i,j) = fluxschv(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'fluxschv', 1.e3, iu, -2)

        do j=1,ny
          do i=1,nx
            if (abs(fluxschu(i,j)).ne.0.) then
              work(i,j) = min (1.e3, max (-1.e3, 
     *                                    fluxgrdu(i,j)/fluxschu(i,j) ))
            else
              work(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice, work, 'fluxgrdu/schu', 1., iu, -1)

        do j=1,ny
          do i=1,nx
            if (abs(fluxschv(i,j)).ne.0.) then
              work(i,j) = min (1.e3, max (-1.e3, 
     *                                    fluxgrdv(i,j)/fluxschv(i,j) ))
            else
              work(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice, work, 'fluxgrdv/schv', .1, iu, -2)

        do iuvloop=1,2
          do j=1,ny
            do i=1,nx
              if (iuvloop.eq.1) zang = angnormu(i,j)
              if (iuvloop.eq.2) zang = angnormv(i,j)
              if (zang.eq.-999.) then      ! set in schoofgl for no g.l.
                work(i,j) = 0.
              else if (abs(zang).le.5.) then
                work(i,j) = 1000.          ! show small values by "*"
              else if (zang.gt.180.) then
                work(i,j) = zang - 360.
              else 
                work(i,j) = zang
              endif
            enddo
          enddo
          if (iuvloop.eq.1) then
            call printmap (timeice, work, 'angnormu', 10., iu, -1)
          else
            call printmap (timeice, work, 'angnormv', 10., iu, -2)
          endif
        enddo
c-----
#endif
c-----

        do j=1,ny
          do i=1,nx
            work(i,j) = hgu(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'hgu (m)', 40., iu, 1)

        do j=1,ny
          do i=1,nx
            work(i,j) = hgv(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'hgv (m)', 40., iu, 2)

        call printmap (timeice, thetau, 'thetau', .05, iu, 0)
        call printmap (timeice, thetav, 'thetav', .05, iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = ub(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'ub',  10.0, iu, 1)
  
        do j=1,ny
          do i=1,nx
            work(i,j) = vb(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'vb',  10.0, iu, 2)

        do j=1,ny
          do i=1,nx
            work(i,j) = sqrt(  (0.5*(ub(i,j)+ub(i-1,j)))**2 
     *                       + (0.5*(vb(i,j)+vb(i,j-1)))**2 )
          enddo
        enddo
        call printmap (timeice,work, '|uvb|',   0.1, iu, 0)
        call printmap (timeice,work, '|uvb|',    1., iu, 0)
        call printmap (timeice,work, '|uvb|',  100., iu, 0)
  
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              zua = sqrt(  (0.5*(ua(i,j)+ua(i-1,j)))**2 
     *                   + (0.5*(va(i,j)+va(i,j-1)))**2 )
              zub = sqrt(  (0.5*(ub(i,j)+ub(i-1,j)))**2 
     *                   + (0.5*(vb(i,j)+vb(i,j-1)))**2 )
              work(i,j) = min (1.e6, zub / max (zua,1.e-20))
            else
              work(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice,work, '|uvb/uva|',  0.01, iu, 0)
        call printmap (timeice,work, '|uvb/uva|',  0.1,  iu, 0)
  
        do j=1,ny
          do i=1,nx
            if (h(i,j).gt.0.) then
              work(i,j) = sqrt(  (0.5*(ua(i,j)+ua(i-1,j)))**2 
     *                         + (0.5*(va(i,j)+va(i,j-1)))**2 )
            else
              work(i,j) = 0.
            endif
          enddo
        enddo
        call printmap (timeice,work, '|uvabs|',       10., iu, 0)
        call printmap (timeice,work, '|uvabs|_100',  100., iu, 0)
  
        do j=1,ny
          do i=1,nx
            zua = 0.5*(ua(i,j)+ua(i-1,j))
            zva = 0.5*(va(i,j)+va(i,j-1))
            work(i,j) = 1.e-6*h(i,j)*sqrt(zua**2+zva**2)
          enddo
        enddo
        call printmap (timeice,work, '|uvhabs|', 0.01, iu, 0)
    
        call zero (work,nx*ny)
        do k=1,nlev
          do j=1,ny
            do i=1,nx
              work(i,j) = work(i,j) + heati(i,j,k)*dzeta(k)*h(i,j)
     *                                / (rhoice*hfus)
            enddo
          enddo
        enddo
        call printmap (timeice,work,'heati (m/y)', 0.001, iu, 0)
    
        do j=1,ny
          do i=1,nx
            work(i,j) = heath(i,j)*h(i,j) / (rhoice*hfus)
          enddo
        enddo
        call printmap (timeice,work,'heath (m/y)', 0.001, iu, 0)
    
        call zero (work,nx*ny)
        do j=1,ny
         do i=1,nx
            work(i,j) = heatb(i,j) / (rhoice*hfus)
          enddo
        enddo
        call printmap (timeice,work,'heatb (m/y)', 0.001, iu, 0)

        do j=1,ny
          do i=1,nx
c           if (maskwater(i,j).ne.1) then
              work(i,j) = h(i,j) - (rholiq/rhoice)*(sealev-hb(i,j))
c           else
c             work(i,j) = 0.
c           endif
          enddo
        enddo
        call printmap (timeice,work,'hflot (m)', 10., iu, 0)

        call printmap (timeice,hs,'hs (m)', 150., iu, 0)
        call printmap (timeice,hsbal,'hsbal (m)', 150., iu, 0)
        do j=1,ny
          do i=1,nx
            work(i,j) = hs(i,j) - hsbal(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'hs-hsbal (m)', 50., iu, 0)
        call itor (maskbal, work, nx*ny)
        call printmap (timeice,work, 'maskbal',       1., iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j) = hb(i,j) - hbinit(i,j)
          enddo   
        enddo   
        call printmap (timeice,work, 'hb-hbinit (m)', 20., iu, 0)

        call printmap (timeice, tempocn,'ocean temp (C)',   1., iu, 0)
        call printmap (timeice, tempocn,'ocean temp (C)',  0.2, iu, 0)

        call printmap (timeice,h,  'ice thickness (m)',   150., iu, 0)
        call printmap (timeice,h,  'ice thickness_10 (m)', 10., iu, 0)
        call printmap (timeice,hedge, 'hedge (m)',        150., iu, 0)
        call printmap (timeice,hedge, 'hedge_10 (m)',      10., iu, 0)
        call printmap (timeice,fedge, 'fedge',             0.1, iu, 0)

c--------------------
#if defined (MELANGE)
c--------------------
        call itor (maskmel, work, nx*ny)
        call printmap (timeice, work, 'maskmel (-1,0,1)', 1., iu, 0)
        call printmap (timeice, hmel,  'hmel (m)',       50., iu, 0)
        call printmap (timeice, hmel,  'hmel_10 (m)',    10., iu, 0)
        call printmap (timeice, fmel,  'fmel (0-1)',    0.1,  iu, 0)
        call printmap (timeice, hwmel, 'hwmel (m)',      40., iu, 0)
        call printmap (timeice, budgallmel,'budgallmel (m/y)',  .1,iu,0)
        call printmap (timeice, budgallmel,'budgallmel (m/y)',  1.,iu,0)
        call printmap(timeice,budgsupmel,'budgsupmel_10 (m/y)',10.,iu,0)
        call printmap(timeice,budgadjmel,'budgadjmel_10 (m/y)',10.,iu,0)
        call printmap(timeice,budgsmomel,'budgsmomel_100(m/y)',100.,
     *                iu, 0)
        call printmap(timeice,budgsmomel,'budgsmomel_1000(m/y)',1000.,
     *                iu, 0)
        call printmap(timeice,oceanmeltmel,'oceanmeltmel_10(m/y)',10.,
     *                iu, 0)
        call printmap(timeice,oceanmeltmel,'oceanmeltmel_1(m/y)',1.,
     *                iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j)  = uamel(i,j)
            work2(i,j) = vamel(i,j)
          enddo
        enddo
        call printmap (timeice,work, 'uamel (m/y)',         10.0, iu, 1)
        call printmap (timeice,work, 'uamel_100 (m/y)',    100.0, iu, 1)
        call printmap (timeice,work, 'uamel_1000 (m/y)',  1000.0, iu, 1)
        call printmap (timeice,work, 'uamel_10000 (m/y)',10000.0, iu, 1)
        call printmap (timeice,work2,'vamel (m/y)',         10.0, iu, 2)
        call printmap (timeice,work2,'vamel_100 (m/y)',    100.0, iu, 2)
        call printmap (timeice,work2,'vamel_1000 (m/y)',  1000.0, iu, 2)
        call printmap (timeice,work2,'vamel_10000 (m/y)',10000.0, iu, 2)

        call printmap (timeice, thetaumel, 'thetaumel', .05, iu, 0)
        call printmap (timeice, thetavmel, 'thetavmel', .05, iu, 0)

        do j=1,ny
          do i=1,nx
            work(i,j)  = bpumel(i,j)
            work2(i,j) = bpvmel(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'bpumel', .05, iu, 1)
        call printmap (timeice, work2, 'bpvmel', .05, iu, 2)

        do j=1,ny
          do i=1,nx
            work(i,j)  = uaboundmel(i,j)
            work2(i,j) = vaboundmel(i,j)
          enddo
        enddo
        call printmap (timeice, work, 'uaboundmel_100 (m/y)',100.,iu,1)
        call printmap (timeice, work, 'uaboundmel_1000(m/y)',1000.,iu,1)
        call printmap (timeice, work2,'vaboundmel_100 (m/y)',100.,iu,2)
        call printmap (timeice, work2,'vaboundmel_1000(m/y)',1000.,iu,2)
c-----
#endif
c-----

        call flush (iu)

c~~~~~~~~~~
      enddo
c~~~~~~~~~~
      return
      end

c-----------------------------------------------------------------------

      subroutine printmap (time, arrin, cdescin, scale, iu, igrid)

c        Displays array arr to unit iu in "a1" format.
c        igrid = 0 for h-grid variables, 1 for u-grid, 2 for v-grid.

#include <comicegrid.h>

      dimension arrin(nx,ny)
      character*(*) cdescin

      dimension arr(nx,ny), maskdisplay(nx,ny)
      character*1 carr(nx)
      character cdesc*80, chem*80, cbrr*20000
      logical first
      data first /.true./
      save iac, jac

      character crev*4, csul*4, coff*3
      parameter (crev = char(27)//'[7m',           ! reverse video
     *           csul = char(27)//'[4m',           ! underline
     *           coff = char(27)//'[m' )           ! reset all
c     character cbra*1
c     integer*1 icbra
c     equivalence (icbra,cbra)
c     data icbra /O'154000000000000000000'/        ! "ctrl["
c     save crev, csul, coff

      cdesc = cdescin

      if (first) then
c       crev = cbra // '[7m'
c       csul = cbra // '[4m'
c       coff = cbra // '[m'
        iac = ichar('A') - 1
        jac = ichar('a') - 1
      endif

c       Set mask for reversed-video land-ocean outline
c       (maskcur, maskbal are in common)

      if (cdesc(1:4).eq.'hbal' .or.  cdesc(1:5).eq.'hsbal') then
        call scopy_i (nx*ny, maskbal, 1, maskdisplay, 1)
      else
        call scopy_i (nx*ny, maskcur, 1, maskdisplay, 1)
      endif

      call scopy (nx*ny, arrin, 1, arr, 1)

c     if (nx.le.1500) then
        nhem = 1
c     else
c       nhem = 2
c     endif


#if defined (TEST2D)
c     for mismipplus:
      ic = max (1, nx/640)
#elif defined (TESTMEL2D)
      ic = max (1, nx/300)
#elif defined (GREENLAND)
      ic = max (1, nx/350)
#else
      ic = max (1, nx/300)
c     ic = max (1, nx/240)
c     ic = 1
#endif

c-------------------
      do ihem=1,nhem
c-------------------

        if (nhem.eq.1) then
          chem = ' '
          ia = 1
          ib = nx
        else if (ihem.eq.1) then
          chem = 'western: '
          ia = nx/2 + 1
          ib = nx
        else
          chem = 'eastern: '
          ia = 1
          ib = nx/2
        endif

        if (.not. first) write (iu,*) 
        write(iu,100)
     *    time, chem(1:lenchr(chem)), cdesc(1:lenchr(cdesc)), scale
#if defined (LONLAT)
#  if defined (NESTING)
        za = 10.
#  else
        za = 1.
#  endif
        if (iabs(igrid).eq.0 .or. iabs(igrid).eq.2) then
          write(iu,'(/6x,2x,500i5)') 
     *      (nint(  za*alond(i,(ny+1)/2)), i=ia+4*ic,ib,5*ic)
        else
          write(iu,'(/6x,2x,500i5)') 
     *    (nint(za*0.5*(alond(i,(ny+1)/2)+alond(min(i+1,nx),(ny+1)/2))),
     *      i=ia+4*ic,ib,5*ic)
        endif
#else
        if (iabs(igrid).eq.0 .or. iabs(igrid).eq.2) then
          write(iu,'(/6x,2x,500i5)') 
     *      (nint(.001*xh(i,(ny+1)/2)), i=ia+4*ic,ib,5*ic)
        else
          write(iu,'(/6x,2x,500i5)') 
     *      (nint(.001*0.5* (xh(i,(ny+1)/2)+xh(min(i+1,nx),(ny+1/2)))),
     *       i=ia+4*ic,ib,5*ic)
        endif
#endif
  100   format(
     *    'time:',f12.2,'  ',a,'  ',a,'  /',f11.4,
     *    '  (A-Z pos, a-z neg):')

c>>>>>>>>>>>>>>>>>>>>>>>>
        do j = ny, 1, -ic
c>>>>>>>>>>>>>>>>>>>>>>>>

          do i=ia,ib
            if (abs(arr(i,j)).le.1.e-8) then
              carr(i) = '.'
            else
              zs = max (-1.e6, min (1.e6, arr(i,j)/scale))
              n = nint(zs)
              if (n.eq.0) then
                if (arr(i,j).eq.0.) then
                  carr(i) = '.'
                else if (arr(i,j).gt.0.) then
                  carr(i) = '+'
                else
                  carr(i) = '-'
                endif
              else if (n.ge.1..and. n.le.26) then
                carr(i) = char(iac + n)
              else if (n.le.-1..and. n.ge.-26) then
                carr(i) = char(jac - n)
              else if (n.gt.26) then
                carr(i) = '*'
              else if (n.lt.-26) then
                carr(i) = '#'
              endif
            endif
          enddo

          if (.false.) then

c              Show ocean (not lakes) as blanks

            lenb = 0
            do i=ia,ib,ic 
              lenb = lenb + 1
              if (maskdisplay(i,j).eq.1) then    ! maskdisplay in common
                cbrr(lenb:lenb) = ' '
              else
                cbrr(lenb:lenb) = carr(i)
              endif
            enddo
 
          else

c              Reverse-video for ocean/lake edge pts
c              (using h-grid maskdisplay, in common). igrid (passed) is
c              0 for h-grid arrays, 1 for u-grid, 2 for v-grid.
c              Also underline where ice (h > 0, ifi)...not done now.

            istatprev = 0
            lenb = 0
            do i = ia,ib,ic
              ifi = 0
cold          if (h(i,j).gt.0.) ifi = 1          ! doesn't work with lwy
              ifo = 0
              im = max (i-ic,1)
              ip = min (i+ic,nx)
              jm = max (j-ic,1)
              jp = min (j+ic,ny)
c             if (maskdisplay(i,j).ne.0) ifo = 1     ! for any water pts

c             if (igrid.eq.0) then
                if(
     *           ( maskdisplay(i,j).eq.0 .and. 
     *             (maskdisplay(im,j).ne.0.or.maskdisplay(ip,j).ne.0.or.
     *              maskdisplay(i,jm).ne.0.or.maskdisplay(i,jp).ne.0) 
     *           ) 
     *            ) ifo = 1
c             else if (igrid.eq.1) then
c               if ( (maskdisplay(i,j).eq.0.and.maskdisplay(ip,j).ne.0.)
c    *               .or.
c    *               (maskdisplay(i,j).ne.0.and.maskdisplay(ip,j).eq.0.)
c    *             ) ifo = 1
c             else if (igrid.eq.2) then
c               if ( (maskdisplay(i,j).eq.0.and.maskdisplay(i,jp).ne.0.)
c    *               .or.
c    *               (maskdisplay(i,j).ne.0.and.maskdisplay(i,jp).eq.0.)
c    *             ) ifo = 1
c             endif

              if (igrid.lt.0) ifo = 0  ! for g.l. related quantities ?

              istat = ifi + 2*ifo

              if (istat.ne.istatprev) then
                cbrr(lenb+1:lenb+3) = coff
                lenb = lenb + 3
                if (ifi.eq.1) then
                  cbrr(lenb+1:lenb+4) = csul
                  lenb = lenb + 4
                endif
                if (ifo.eq.1) then
                  cbrr(lenb+1:lenb+4) = crev
                  lenb = lenb + 4
                endif
              endif

              cbrr(lenb+1:lenb+1) = carr(i)
              lenb = lenb + 1

              if (i.gt.ib-ic .and. istat.ne.0) then
                cbrr(lenb+1:lenb+3) = coff
                lenb = lenb + 3
              endif

              istatprev = istat
            enddo

          endif

#if defined (LONLAT)
          if (iabs(igrid).eq.0 .or. iabs(igrid).eq.1) then
            write(iu,'(f6.2,2x,a)') alatd((nx+1)/2,j), cbrr(1:lenb)
          else 
            write(iu,'(f6.2,2x,a)')
     *         0.5*( alatd((nx+1)/2,j) + alatd((nx+1)/2,min(j+1,ny)) ),
     *         cbrr(1:lenb)
          endif
#else
          if (iabs(igrid).eq.0 .or. iabs(igrid).eq.1) then
            write(iu,'(i6,2x,a)') nint(.001*yh((nx+1)/2,j)),cbrr(1:lenb)
          else
            write(iu,'(i6,2x,a)') 
     *       nint(.001*0.5*(yh((nx+1)/2,j) + yh((nx+1)/2,min(j+1,ny)))),
     *       cbrr(1:lenb)
          endif
#endif

c>>>>>>>>>>>>
        enddo
c>>>>>>>>>>>>

c----------
      enddo
c----------

      first = .false.

      return
      end

c-----------------------------------------------------------------------

      subroutine iceplot1d (h, hs, hb, t, tracer, 
     *                      budgall, basefrml, 
     *                      oceanmelt, facemelt, cliffmelt,
     *                      tsurf, tsurfi, heati, heath, heatb,
     *                      w, wa, maskh,
     *                      sedim, tsed, wsed,heats,
     *                      quarryrate, quarryacc, quarrycoef,
     *                      pelagrate, pelagacc,
     *                      topbed, topbedeq, tbed,
     *                      hw, tw,
     *                      u, v, ua, va, uia, via, ub, vb, 
     *                      uadv, vadv, hu, hv, 
     *                      masku, maskv, crhu, crhv,
     *                      fracgu, fracgv, thetau, thetav,
     *                      sealev, dtantann, dtantjan, dtseas, rco2,
     *                      timeice, dt, weirun,
     *                      iffirse, iflast, nyeartot, nyearplot1d)
  
c     Writes ascii output file for 1-D flowline plotting, within 
c     1-D or 2-D domains. Flowline is defined by i1s,j1s,n1s
c     in comicegrid, set in iceshow1d.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  t(nx,ny,0:nlevp),    tracer(nx,ny,0:nlevp,ntrace),
     *  budgall(nx,ny),      basefrml(nx,ny),    
     *  oceanmelt(nx,ny),    facemelt(nx,ny),    cliffmelt(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  heati(nx,ny,nlev),   heath(nx,ny),       heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  maskh(nx,ny),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     *  heats(nx,ny,nsed),   
     *  quarryrate(nx,ny),   quarryacc(nx,ny),   quarrycoef(nx,ny), 
     *  pelagrate(nx,ny),    pelagacc(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  hw(nx,ny),           tw(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),       via(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),     maskv(0:nxp,0:nyp), 
     *  fracgu(0:nxp,0:nyp),    fracgv(0:nxp,0:nyp),
     *  thetau(nx,ny),          thetav(nx,ny)

      dimension work(0:(nxp+nyp))
      character*80 cform, cformi

c     number of cross-sec vertical points, and bottom and top extents
#if defined (EISLINE) && defined (LINEB)
      parameter (nz=200, zbot = -1000., ztop=3000.)
#elif defined (EISLINE) && defined (LINED)
      parameter (nz=200, zbot = -2500., ztop=3000.) 
#elif defined (EISLINE) && defined (LINEE)
      parameter (nz=200, zbot =  -600., ztop= 100.) 
#elif defined (EISLINE) && defined (LINEF)
      parameter (nz=200, zbot = -1000., ztop=4000.)
#elif defined (EISLINE) && defined (LINEG)
      parameter (nz=200, zbot = -1000., ztop=4000.)
#elif defined (EISLINE) && defined (SLMODEL)
      parameter (nz=  0, zbot =    0., ztop=    0.)      ! no 2-D output
#elif defined (EISANTA) && defined (TRANSECTA)
      parameter (nz=200, zbot = -1000., ztop=4000.)
#else
c     parameter (nz=200, zbot = -1000., ztop=3000.) 
c     parameter (nz=200, zbot = -2500., ztop=3000.) 
      parameter (nz=200, zbot = -2500., ztop=4000.) 
#endif

      dimension 
     *  work1d(nx+ny),
     *  work1dv(0:nx+ny+1),   iwork1dv(0:nx+ny+1),
     *  work2d(0:nxp,0:nyp),  work2da(0:nxp,0:nyp)

c        Skip out if not time

      if (nyearplot1d.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyearplot1d),dt)).lt.dt
     *       .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyearplot1d.eq.-1 )
     *   ) return 

      iu = iuplot1d
c     cform  = '(a/1001f12.5)'
c     cformi = '(a/1001i12  )'
c     cform  = '(a/(300f12.5))'
c     cformi = '(a/(300i12  ))'
      cform  = '(a/(400f10.3))'  ! berg   ! nb: gf limit 4096 chars/line
      cformi = '(a/(400i10  ))'  ! berg   ! see linesper below

      if (iffirse.eq.1) then
        nwrite = nint(float(nyeartot)/float(nyearplot1d)) + 1
#if defined (EISLINE) && defined (LINEE)
        nwrite =  4000./nyearplot1d + 1    ! bb,cc
#endif
        write (iu,'(a,i8)') 'nwriteplot=',nwrite
      endif

#if defined (EISLINE) && defined (LINEE)
      if (nint(timeice).lt.6000)  return ! bb,cc
#endif

c        Write 0-D quantities

      igl = 0
      fgl = 0.
      xgl  = 0.
      if (nx.ge.ny) then
        do m=1,n1s-1
          if (igl.eq.0 .and. fracgu(i1s(m),j1s(m)).lt.1.) then 
            igl = i1s(m)
            fgl = fracgu(i1s(m),j1s(m))
            xgl = (1.-fgl)*xh(i1s(m),j1s(m))+fgl*xh(i1s(m+1),j1s(m+1))
          endif
        enddo
      else 
        do m=1,n1s-1
          if (igl.eq.0 .and. fracgv(i1s(m),j1s(m)).lt.1.) then 
            igl = j1s(m)
            fgl = fracgu(i1s(m),j1s(m))
            xgl = (1.-fgl)*yh(i1s(m),j1s(m))+fgl*yh(i1s(m+1),j1s(m+1))
          endif
        enddo
      endif 

      write (iu,'(80("="))')
c     write (iu,'( a,i9/a,f9.4/4(a,i9/),11(a,f9.3/),a,i9)')
c     write (iu,'( a,i9/a,f9.4/4(a,i9/),11(a,f9.3/),3(a,i9))')   ! berg
      write (iu,'( a,i9/a,3a/a,f9.4/4(a,i9/),11(a,f9.3/),a,i9)') ! 34m
     *  'time     =',  nint(timeice),
     *  'ctitle   =',  " '",c1s(1:lenchr(c1s)),"'",              ! 34m
     *  'weirun   =',  weirun,
     *  'nx       =',  n1s,
     *  'nlev     =',  nlev,
     *  'nbed     =',  nbed,
     *  'nz       =',  nz,
     *  'xlef     =',  (d1s(1)   - 0.5*dx0)*1.e-3,
     *  'xrit     =',  (d1s(n1s) + 0.5*dx0)*1.e-3,
     *  'zbot     =',  zbot,
     *  'ztop     =',  ztop,
     *  'sealev   =',  sealev,
     *  'dtantann =',  dtantann,
     *  'dtantjan =',  dtantjan,
     *  'dtseas   =',  dtseas,     ! new..change ginice...
     *  'rco2     =',  rco2,       ! new..change ginice...
     *  'xgl      =',  xgl*1.e-3,
     *  'fgl      =',  fgl,
     *  'igl      =',  igl
c    * ,'   linesper  =', (n1s-1)/400 + 1,     ! berg, not 34m
c    *  '   linesperu =', (n1s  )/400 + 1      ! berg, not 34m

c        Write 1-D quantities (h-grid)

c     write(iu,cform) 'xh       :', (xh(i1s(m),j1s(m)),      m=1,n1s)
      write(iu,cform) 'd1s      :', (d1s(m)*1.e-3,           m=1,n1s)
      write(iu,cform) 'h        :', (h(i1s(m),j1s(m)),       m=1,n1s)
      write(iu,cform) 'hs       :', (hs(i1s(m),j1s(m)),      m=1,n1s)
      write(iu,cform) 'hb       :', (hb(i1s(m),j1s(m)),      m=1,n1s)
      write(iu,cform) 'sedim    :', (sedim(i1s(m),j1s(m)),   m=1,n1s)
c     place holder, old sedimun:
      write(iu,cform) 'sedimun  :', (sedim(i1s(m),j1s(m)),   m=1,n1s)
      do m=1,n1s
        work(m) = quarryrate(i1s(m),j1s(m))*1.e6
      enddo
      write(iu,cform) 'quarryrat:', (work(m),                m=1,n1s)
      write(iu,cform) 'hw       :', (hw(i1s(m),j1s(m)),      m=1,n1s)
      write(iu,cform) 'topbedeq :', (topbedeq(i1s(m),j1s(m)),m=1,n1s)
      write(iu,cform) 'thetau   :', (thetau(i1s(m),j1s(m)),  m=1,n1s)
      write(iu,cform) 'budgall  :', (budgall(i1s(m),j1s(m)), m=1,n1s)
      write(iu,cform) 'basefrml :', (basefrml(i1s(m),j1s(m)),m=1,n1s)
      write(iu,cform) 'oceanmelt:',(oceanmelt(i1s(m),j1s(m)),m=1,n1s)
      write(iu,cform) 'tw       :', (tw(i1s(m),j1s(m))-tmelt,m=1,n1s)
      do m=1,n1s
        work(m) = 0.
        do k=1,nlev
          work(m) = work(m) 
     *            + heati(i1s(m),j1s(m),k) * dzeta(k)*h(i1s(m),j1s(m))
        enddo
        work(m) = work(m) + heath(i1s(m),j1s(m)) * h(i1s(m),j1s(m))
      enddo
      write(iu,cform) 'heati    :', (work(m)/31556926,       m=1,n1s)
      write(iu,cform) 'heatb    :', 
     *                       (heatb(i1s(m),j1s(m))/31556926,m=1,n1s)
      do m=1,n1s
        work(m) = t(i1s(m),j1s(m),nlevp) 
     *          - (tmelt-dtmdh*h(i1s(m),j1s(m)))
      enddo
      write (iu,cform) 'tbme    :', (work(m),                m=1,n1s)

c        Write 1-D quantities (u-grid)

      call flowline (ua, va, work1dv, 1, 1)
      write(iu,cform) 'ua       :', (work1dv(m),             m=0,n1s)
      call flowline (ub, vb, work1dv, 1, 1)
      write(iu,cform) 'ub       :', (work1dv(m),             m=0,n1s)

      do j=0,ny
        do i=0,nx
          work2d (i,j) = log10(crhu(i,j))
          work2da(i,j) = log10(crhv(i,j))
        enddo
      enddo
      call flowline (work2d, work2da, work1dv, 1, 0)
      write(iu,cform) 'crhu     :', (work1dv(m),             m=0,n1s)

c          Write 2-D quantities

      if (nz.gt.0) then

        call crossec ( t,          tmelt, h, hs, hb, topbed, sedim, 
     *                't',   iu, 1, 1,  0,nlevp,nz, zbot, ztop)

        call crossec ( tsed,       tmelt, h, hs, hb, topbed, sedim,
     *                'tsed',iu, 2, 1,   1,nsed,nz, zbot, ztop)

        call crossec ( tbed,       tmelt, h, hs, hb, topbed, sedim,
     *                'tbed',iu, 3, 1,   1,nbed,nz, zbot, ztop)

c       u -grid:
        call crossecuv (u,   v,    0.  , h, hs, hb, topbed, sedim, 
     *                  'u',  iu, 1, 2, 0, nlevp, nz, zbot, ztop)

      endif

      call flush (iu)

      return
      end

c-----------------------------------------------------------------------

      subroutine crossec (arr,       shift, h, hs, hb, topbed, sedim,
     *                    cname, iu, iflag, igrid, n1,n2,nz, zbot, ztop)

      entry      crossecuv(arru, arrv, shift, h, hs, hb, topbed, sedim,
     *                    cname, iu, iflag, igrid, n1,n2,nz, zbot, ztop)

c     Called by iceplot1d, for flowline profiles.
c
c     Writes a 2-d ascii cross-section (x vs z, at j, for plotting) 
c     of array arr or arru, converting from vertical within-ice,
c     within-sed, or within-bed "sigma" coord to physical z coord. 
c     Accumulate only those points in x,z space in the cross sec
c     (ice, sed or bedrock), and write only those to iu for gf.
c
c     iflag = 1 for ice sheet , 2 for sediment, 3 for bedrock
c     igrid = 1 for h-grid (arr), 2 for u-grid (arru, interp here to h)

c     Uses flowline definition i1s,j1s,n1s in comicegrid.

#include <comicephys.h>
#include <comicegrid.h>

      dimension arr(nx,ny,n1:n2), 
     *          arru(0:nxp,0:nyp,n1:n2), arrv(0:nxp,0:nyp,n1:n2),
     *          h(nx,ny), hs(nx,ny), hb(nx,ny),
     *          topbed(nx,ny), sedim(nx,ny)
      character*(*) cname

      parameter (ncrosslim = (nx+ny)*200)
      dimension cross(ncrosslim),  indcross(ncrosslim,2)

      parameter (nfine=10000)
      dimension indfine   (nfine,2), weifine(nfine), 
     *          indfinesed(nfine,2), weifinesed(nfine), 
     *          indfinebed(nfine,2), weifinebed(nfine) 
      save indfine,    weifine, 
     *     indfinesed, weifinesed,
     *     indfinebed, weifinebed  

      logical first
      save first
      data first /.true./

c        Set regularly spaced fine resolution index into zeta/zsed/zbed
c        ([ind,wei]fine for ice, [ind,wei]finesed for sedim, 
c         [ind,wei]finebed for bedrock) 

      if (first) then

c       for ice sheet:
        do m=1,nfine
          zfine = (m-.5)/nfine
          do k=0,nlev
            if (zeta(k).le.zfine .and. zeta(k+1).ge.zfine) then
              indfine(m,1) = k 
              indfine(m,2) = k+1 
              weifine(m) = (zeta(k+1)-zfine) / (zeta(k+1)-zeta(k))
              goto 10
            endif
          enddo
          write (6,*) 'Error (crossec): setting indfine, weifine'
          stop
   10     continue
        enddo

c       for sedim:
        do m=1,nfine
          zfine = (m-.5)/nfine
          if (zfine.le.zsed(1)) then
            indfinesed(m,1) = 1 
            indfinesed(m,2) = 1 
            weifinesed(m) = 1.
          else if (zfine.ge.zsed(nsed)) then
            indfinesed(m,1) = nsed 
            indfinesed(m,2) = nsed 
            weifinesed(m) = 1.
          else
            do k=1,nsed-1
              if (zsed(k).le.zfine .and. zsed(k+1).ge.zfine) then
                indfinesed(m,1) = k 
                indfinesed(m,2) = k+1
                weifinesed(m) = (zsed(k+1)-zfine) / (zsed(k+1)-zsed(k))
                goto 20
              endif
            enddo
            write(6,*) 'Error (crossec): setting indfinesed, weifinesed'
            stop
   20       continue
          endif
        enddo

c       for bedrock:
        do m=1,nfine
          zfine = (m-.5)/nfine
          if (zfine.le.zbed(1)) then
            indfinebed(m,1) = 1 
            indfinebed(m,2) = 1 
            weifinebed(m) = 1.
          else if (zfine.ge.zbed(nbed)) then
            indfinebed(m,1) = nbed 
            indfinebed(m,2) = nbed 
            weifinebed(m) = 1.
          else
            do k=1,nbed-1
              if (zbed(k).le.zfine .and. zbed(k+1).ge.zfine) then
                indfinebed(m,1) = k 
                indfinebed(m,2) = k+1
                weifinebed(m) = (zbed(k+1)-zfine) / (zbed(k+1)-zbed(k))
                goto 30
              endif
            enddo
            write(6,*) 'Error (crossec): setting indfinebed, weifinebed'
            stop
   30       continue
          endif
        enddo

        first = .false.
      endif

c        Fill cross and indcross,  array, vertically interpolating 
c        wrt zeta or zbed. npoi is counter for points in the cross-sec.

      npoi = 0
c================
      do m=1,n1s
c================
        i = i1s(m)
        j = j1s(m)
        if ( (h(i,j).gt.0.     .and. iflag.eq.1) .or.
     *       (sedim(i,j).gt.0. .and. iflag.eq.2) .or.
     *       (                       iflag.eq.3) ) then

          if (iflag.eq.1) then
            za = hs(i,j) - h(i,j)
            zb = hs(i,j)
          else if (iflag.eq.2) then
            za = hb(i,j) - sedim(i,j)
            zb = hb(i,j)
          else if (iflag.eq.3) then
            za = topbed(i,j) - bedthick
            zb = topbed(i,j)
          endif

c         if (iflag.eq.1) then
c           dz = (ztop-zbot)/nz
c         else
            dz = 0.
c         endif

          do k=1,nz
            z = zbot + ((k-.5)/nz)*(ztop-zbot)
            if (z.ge.za .and. z.le.zb) then
              zsig = (zb-z)/(zb-za)

c             use regularly spaced fine res zeta-index (for speed)
              ifine = max (1, min (nfine, nint(zsig*nfine + 0.5) ))
              if (iflag.eq.1) then
                ka  = indfine(ifine,1)
                kb  = indfine(ifine,2)
                zwa = weifine(ifine)
              else if (iflag.eq.2) then
                ka  = indfinesed(ifine,1)
                kb  = indfinesed(ifine,2)
                zwa = weifinesed(ifine)
              else if (iflag.eq.3) then
                ka  = indfinebed(ifine,1)
                kb  = indfinebed(ifine,2)
                zwa = weifinebed(ifine)
              endif
              npoi = npoi + 1
              if (npoi.gt.ncrosslim) then
                write (6,*) 'Error (crossec): too many points'
                stop
              endif
              if (igrid.eq.1) then
                cross(npoi) = zwa*arr(i,j,ka) + (1.-zwa)*arr(i,j,kb)
              else if (igrid.eq.2) then
c either:
c               if (j1s(m).eq.j1s(min(m+1,n1s))) then
c                 cross(npoi)=    zwa *0.5*(arru(i,j,ka)+arru(i-1,j,ka))
c    *                      + (1.-zwa)*0.5*(arru(i,j,kb)+arru(i-1,j,kb))
c    *                      - shift
c               else
c                 cross(npoi)=    zwa *0.5*(arrv(i,j,ka)+arrv(i,j-1,ka))
c    *                      + (1.-zwa)*0.5*(arrv(i,j,kb)+arrv(i,j-1,kb))
c    *                      - shift
c               endif
c or:
                zu =       zwa *0.5*(arru(i,j,ka)+arru(i-1,j,ka))
     *               + (1.-zwa)*0.5*(arru(i,j,kb)+arru(i-1,j,kb)) -shift
                zv =       zwa *0.5*(arrv(i,j,ka)+arrv(i,j-1,ka))
     *               + (1.-zwa)*0.5*(arrv(i,j,kb)+arrv(i,j-1,kb)) -shift
                cross(npoi) = sqrt(zu**2 + zv**2)
              endif
              indcross(npoi,1) = m
              indcross(npoi,2) = k
            endif
          enddo


        endif
c==========
      enddo
c==========

c        Write cross section to output file

      write (iu,"(a,':')") cname
      write (iu,'(a,i9)')
     *  'npoi     =',  npoi

      if (npoi.gt.0) then
        do m=1,npoi
          write (iu,"(a,'(',i4,',',i4,') = ',f12.5)")
     *      cname, indcross(m,1), indcross(m,2), cross(m)
        enddo
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine icetab1d (h, hs, hb, t, tracer, 
     *                     budgall, basefrml, 
     *                     oceanmelt, facemelt, cliffmelt,
     *                     tsurf, tsurfi, heati, heath, heatb,
     *                     w, wa, maskh,
     *                     sedim, tsed, wsed, heats, 
     *                     topbed, topbedeq, tbed,
     *                     hw, tw, maskwater, 
     *                     u, v, ua, va, uia, via, ub, vb, 
     *                     uadv, vadv, hu, hv, 
     *                     masku, maskv, muind, mvind, crhu, crhv,
     *                     fracgu, fracgv, thetau, thetav,
     *                     hgu, hgv, qgu, qgv,
     *                     uschoof, vschoof, fschoofu, fschoofv,
     *                     uagrid, vagrid,
     *                     fluxgrdu, fluxgrdv, fluxschu, fluxschv,
     *                     sealev, dtantann, dtantjan, dtseas, rco2,
     *                     timeice, dt, weirun,
     *                     iffirse, iflast, nyeartot, nyeartab)
 
c     Writes ascii tabular file, time series of particular quantities,
c     for 1-D runs

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  t(nx,ny,0:nlevp),    tracer(nx,ny,0:nlevp,ntrace),
     *  budgall(nx,ny),      basefrml(nx,ny),    
     *  oceanmelt(nx,ny),    facemelt(nx,ny),    cliffmelt(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  heati(nx,ny,nlev),   heath(nx,ny),       heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  maskh(nx,ny),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     *  heats(nx,ny,nsed), 
     *  topbed(nx,ny),       topbedeq(nx,ny),    tbed(nx,ny,nbed),
     *  hw(nx,ny),           tw(nx,ny),          maskwater(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),       via(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),     maskv(0:nxp,0:nyp), 
     *  muind(0:nxp,0:nyp),     mvind(0:nxp,0:nyp), 
     *  fracgu(0:nxp,0:nyp),    fracgv(0:nxp,0:nyp), 
     *  thetau(nx,ny),          thetav(nx,ny),
     *  hgu(0:nxp,0:nyp),       hgv(0:nxp,0:nyp),
     *  qgu(0:nxp,0:nyp),       qgv(0:nxp,0:nyp),
     *  uschoof(0:nxp,0:nyp),   vschoof(0:nxp,0:nyp),
     *  uagrid(0:nxp,0:nyp),    vagrid(0:nxp,0:nyp),
     *  fschoofu(0:nxp,0:nyp),  fschoofv(0:nxp,0:nyp),
     *  fluxgrdu(0:nxp,0:nyp),  fluxgrdv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp),  fluxschv(0:nxp,0:nyp)

      save nlineout
      data nlineout /0/ 

      if (nyeartab.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartab),dt)).lt.dt
     *       .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyeartab.eq.-1 )
     *   ) return 

      iu = iutab

c        Total ice volume, non-ocean ice volume (average thicknesses)

      tota    = 0.
      toth    = 0.
      tothnf  = 0.
      totsedz = 0.
      do j=1,ny
        do i=1,nx
          tota = tota +        darea(i,j)
          toth = toth + h(i,j)*darea(i,j)
          if (maskwater(i,j).ne.1) tothnf = tothnf + h(i,j)*darea(i,j)
          totsedz = totsedz + sedim(i,j)*darea(i,j)
        enddo
      enddo
      toth   = toth  /tota
      tothnf = tothnf/tota
      totsedz = totsedz  /tota

c        Position of rightmost grounding line (u-grid) (search for 
c        g.l. from right to left), and thicknesses and elevs each side

      fgl = 0.
      xgl  = 0.
      zhgl = 0.
      ztheta = 0.

      zfscha = 0.
      zfschb = 0.
      zuscha = 0.
      zuschb = 0.
      zhscha = 0.
      zhschb = 0.
      zfluxscha = 0.
      zfluxschb = 0.
      zfluxgrd = 0.
      zuagrid = 0.
      zuadva = 0.
      zuadvb = 0.
      zbgla = 0.
      zbglb = 0.
c!!!!!!!!!!!!!!!!!!!!!!!
      if (nx.ge.ny) then
c!!!!!!!!!!!!!!!!!!!!!!!
        j = (ny+1)/2
        do i=nx-1,1,-1
          if (fracgu(i,j).gt.0.) then
            igl = i
            goto 100
          endif
        enddo
        igl = 0
  100   continue
        kgl = igl
        if (igl.ne.0) then
          fgl = fracgu(igl,j)
          xgl = (1.-fgl)*xh(igl,j)+fgl*xh(igl+1,j)
c         zhgl = (sealev - ((1.-fgl)*hb(igl,j)+fgl*hb(igl+1,j)))/rhor
          zhgl = hgu(igl,j)
          ztheta = thetau(igl+1,j)

          zfscha = fschoofu(igl,j)
          zfschb = fschoofu(igl+1,j)
          zuscha = uschoof(igl,j)
          zuschb = uschoof(igl+1,j)
          zhscha = h(igl,j)
          zhschb = h(igl+1,j)
          zfluxscha = fluxschu(igl,j)
          zfluxschb = zfluxscha 
     *              + darea(igl+1,j)*budgall(igl+1,j)/dyu(igl+1,j)
          zfluxgrd = fluxgrdu(igl,j)
          zuagrid  = uagrid(igl,j)
          zuadva = uadv(igl,j) / dyu(igl,j) 
          zuadvb = uadv(igl+1,j) / dyu(igl+1,j) 
          bgla = 0.
          do i=1,igl
            bgla = bgla + darea(i,j)*budgall(i,j)
          enddo
          bglb = bgla + darea(igl+1,j)*budgall(igl+1,j)
          bgla = bgla / dyu(igl,j)
          bglb = bglb / dyu(igl+1,j)
        endif
c!!!!!!!!!
      else 
c!!!!!!!!!
        i = (nx+1)/2
        do j=ny-1,1,-1
          if (fracgv(i,j).gt.0.) then
            jgl = j
            goto 200
          endif
        enddo
        jgl = 0
  200   continue
        kgl = jgl
        if (jgl.ne.0) then
          fgl = fracgv(i,jgl)
          xgl = (1.-fgl)*yh(i,jgl)+fgl*yh(i,jgl+1)
c         zhgl = (sealev - ((1.-fgl)*hb(i,jgl)+fgl*hb(i,jgl+1)))/rhor
          zhgl = hgv(i,jgl)
          ztheta = thetav(i,jgl+1)

          zfscha = fschoofv(i,jgl)
          zfschb = fschoofv(i,jgl+1)
          zuscha = vschoof(i,jgl)
          zuschb = vschoof(i,jgl+1)
          zhscha = h(i,jgl)
          zhschb = h(i,jgl+1)
          zfluxscha = fluxschv(i,jgl)
          zfluxschb = zfluxscha 
     *              + darea(i,jgl+1)*budgall(i,jgl+1)/dxv(i,jgl+1)
          zfluxgrd = fluxgrdv(i,jgl)
          zuagrid  = vagrid(i,jgl)
          zuadva = vadv(i,jgl) / dxv(i,jgl) 
          zuadvb = vadv(i,jgl+1) / dxv(i,jgl+1) 
          bgla = 0.
          do j=1,jgl
            bgla = bgla + darea(i,j)*budgall(i,j)
          enddo
          bglb = bgla + darea(i,jgl+1)*budgall(i,jgl+1)
          bgla = bgla / dxv(i,jgl)
          bglb = bglb / dxv(i,jgl+1)
        endif
c!!!!!!!!!!
      endif
c!!!!!!!!!!

c        Write header line(s) first call this run


      if (iffirse.eq.1) then
        nwrite = nint(float(nyeartot)/float(nyeartab)) + 1
        write(iu,'(a,i8/2(a,i5))') 'nwritetab=',nwrite,
     *                             'nx=',nx, '  ny=',ny
      endif
#if defined (EISLINE)
      if (iffirse.eq.1) then
#else
      if (iffirse.eq.1 .or. mod(nlineout,50).eq.0) then
#endif
        write (iu,'(/7a)')
     *    '       time',
     *    '    sl dtanta dtantj dtseas   rco2   toth',
     *    '    h500   hs500   sedim',
     *    '  igl    fgl    xgl    hgl  theta',
     *    '  fluxgrd uagrid',
     *    '     #     fscha   uscha  hscha fluxscha    uadva      bgla',
     *    '     #     fschb   uschb  hschb fluxschb    uadvb      bglb'
      endif

      i500 = min ( nx, nint((500.e3/dx0)+0.5) )

      write (iu,'(f11.1, 
     *            i6, 3f7.2, f7.3, f7.1,
     *            3f8.2,
     *            i5, f7.3, f7.1, f7.1, f7.3,
     *            i9, f7.1,
     *            2(9x, f7.3, f8.1, f7.1, 2i9, i10)
     *            )' )
     *  timeice, 
     *  nint(sealev), dtantann, dtantjan, dtseas, rco2, toth, 
     *  h(i500,j), hs(i500,j), totsedz, 
     *  kgl, fgl, xgl*1.e-3, zhgl, ztheta, 
     *  nint(zfluxgrd), zuagrid,
     *  zfscha, zuscha,zhscha,nint(zfluxscha), nint(zuadva), nint(bgla),
     *  zfschb, zuschb,zhschb,nint(zfluxschb), nint(zuadvb), nint(bglb)

      nlineout = nlineout + 1

      call flush (iu)

      return
      end

c-----------------------------------------------------------------------

      subroutine icetab2d (h, hs, hb, hw, maskwater, topbedeq,
     *                     ro18, sealev, 
     *                     dtantann, dtantaus, dtantjan, dtseas, rco2,
     *                     ecc, obl, prec, facice, facorb, facco2,
     *                     timeice, dt, weirun, dtocn34m,
     *                     iffirse, iflast, nyeartot, nyeartab)
 
c     Writes ascii tabular file, time series of particular quantities,
c     for 2-D runs

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny), hs(nx,ny), hb(nx,ny), hw(nx,ny), maskwater(nx,ny),
     *  topbedeq(nx,ny)

c for "esl":
      save firstesl, totslprev, totslwprev, cumesl, cumeslw
      logical firstesl
      data firstesl /.true./

c####################
#if defined (EISANTA)
c####################

c----------------------------------------
c for "esl", wais vs eais region (iutab):
c----------------------------------------
      logical firstwais
      data firstwais /.true./
      save firstwais
      parameter (nwais=10)
      dimension
     *  alonwais(nwais), alatwais(nwais), xwais(nwais), ywais(nwais),
     *  allwais(nx,ny)
      data alatwais
     * / -45., -45., -82.,  -85.,  -86., -82., -75., -72., -45.,  -45./
      data alonwais
     *  /-80., -40., -40.,  -70., -140., 163., 164., 170., 170., -120./
      save allwais

c------------------------------------
c for indiv pts (hb,hw,h,hs) (iupts):
c------------------------------------
      logical firstpts
      data firstpts /.true./
      save firstpts
      parameter (nptsmax=30)
      dimension alonpts(nptsmax), alatpts(nptsmax),
     *          ipts(nptsmax), jpts(nptsmax)
      character*10 cnamepts(nptsmax)
      save alonpts, alatpts, ipts, jpts, cnamepts

#  if defined (PTS_EAIS)
c     For EAIS domes:
      parameter (npts=6)
      data alonpts(1:npts) /
     *       -116.000,     106.800,       124.000,       39.700,    
     *          0.000,      76.517 /
      data alatpts(1:npts) /
    *         -85.000,     -78.467,       -75.000,      -79.300,   
     *         -75.00,     -80.817 /
      data cnamepts(1:npts) /
     *   '       ISA', '    Vostok', '    Dome C', '    Dome F',
     *   '      EDML', '    Dome A' /

#  elif defined (PTS_ACKERT)
c     For R. Ackert:
      parameter (npts=5)
      data alonpts(1:npts) /
     *       -119.517,     -112.085,     -148.816,     -115.167, 
     *       -125.083 /
      data alatpts(1:npts) /
     *        -80.017,      -79.467,      -81.666,      -84.717,  
     *        -77.217 /
      data cnamepts(1:npts) /
     *   '      Byrd', '     WAISD', '     Siple', '      Ohio', 
     *   '    Waesch' /

#  elif defined (PTS_JOHNSON)
c     For J. Johnson (Maish Ntk, Hudson Mts.; Icefall Ntk, Mt. Murphy):
      parameter (npts=2)
      data alonpts(1:npts) /
     *        -99.451,     -111.262 /
      data alatpts(1:npts) /
     *        -74.589,      -75.387 /
      data cnamepts(1:npts) /
     *   'Hudson Mts', 'Mt. Murphy' / 

#  else
c     For RSL (default):
      parameter (npts=10)
      data alonpts(1:npts) /
     *        163.615,      163.778,       78.064,       76.255,    
     *         39.614,      -58.950,      -61.000,      -67.323,
     *        110.474,      -57.333 /
      data alatpts(1:npts) /
     *        -77.323,      -74.886,      -68.588,      -69.380,   
     *        -69.335,      -62.167,      -63.000,      -67.688,
     *        -66.354,      -63.600 /
      data cnamepts(1:npts) /
     *   '     Scott', '     Terra', '    Vestfo', '    Larsem',
     *   '     Syowa', '    Fildes', '     Byers', '    Margue',
     *   '    Windmi', '    J.Ross'  /
#  endif

c#####
#endif
c#####

      if (nyeartab.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartab),dt)).lt.dt
     *       .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyeartab.eq.-1 )
     *   ) return 

      iu = iutab

c        Convert vertices to polar stereo x,y coords, for definition
c        of WAIS vs EAIS used in "esl" calcs

#if defined (EISANTA)
      if (firstwais) then
        do m=1,nwais
          call ps_xy (xwais(m), ywais(m), 
     *                alatwais(m)*pi/180, alonwais(m)*pi/180., 
     *                0., 0.)
        enddo

c       set allwais array (1=inside WAIS polygon, 0=outside):
        call zero (allwais, nx*ny)
        call polygon_wais (allwais, xh,yh,nx,ny, xwais,ywais, nwais,1.)

c       diagnostic dump:
c       do j=ny,1,-1
c         write (133,'(1000i1)') (nint(allwais(i,j)),i=1,nx)
c       enddo
c       call printmap (timeice,allwais,'allwais', 1., 133, 0)

        firstwais = .false.
      endif
#endif

c        Write header line(s) first call, if not a restart

      if (iffirse.eq.1) then
        nwrite = nint(float(nyeartot)/float(nyeartab)) + 1
        write (iu,'(a,i8)') 'nwritetab=',nwrite
#if defined (ISMIP6CTL) || defined (ISMIP6SMB) || defined (ISMIP6BMB) || defined (ABUMIP)
        write (iu,'(a10,27a24)')
     *    '      time',
     *    '    weirun',
     *    '      ro18',
     *    '    sealev',
     *    '    dtanta',
     *    '    dtants',
     *    '    dtantj',
     *    '    dtseas',
     *    '      rco2',
     *    '       ecc',
     *    '       obl',
     *    '      prec',
     *    '    facice',
     *    '    facorb',
     *    '    facco2',
     *    ' toti(km3)',
     *    'totig(km3)',
     *    'totif(km3)',
     *    ' tota(km2)',
     *    'totag(km2)',
     *    'totaf(km2)',
     *    '      h(m)',
     *    '   eofe(m)',
     *    '   eofw(m)',
     *    '    eof(m)',
     *    '   esle(m)',
     *    '   eslw(m)',
     *    '    esl(m)'
#else
        write (iu,'(8a)')
#  if defined (FORCE34MYR)
     *    '      time dtocn34m ro18 sealev',
#  else
     *    '      time weirun   ro18 sealev',
#  endif
     *    ' dtanta dtants dtantj dtseas   rco2',
     *    '     ecc     obl    prec facice facorb facco2',
     *    '   toti(km3)  totig(km3)  totif(km3)',
     *    '   tota(km2)  totag(km2)  totaf(km2)',
     *    '    h(m)',
     *    '   eofe(m)   eofw(m)    eof(m)',
     *    '   esle(m)   eslw(m)    esl(m)'
#endif
      endif

      toti  = 0.
      tota  = 0.
      totif = 0.
      totaf = 0.
      do j=1,ny
        do i=1,nx
          if (h(i,j).ne.0.) then 
            toti = toti + darea(i,j)*h(i,j)
            tota = tota + darea(i,j)
#if defined (ISMIP6CTL) || defined (ISMIP6SMB) || defined (ISMIP6BMB) || defined (ABUMIP)
c           for ismip6, apply smoothing based on height over flotation,
c           to reduce fluctuations in total grounded vs floating ice:  
            zhof = h(i,j) - max(sealev-hb(i,j),0.)/rhor
            zfgr = max (0., min (1., (zhof+25.)/50. ))
            totif = totif + darea(i,j)*h(i,j)*(1.-zfgr)
            totaf = totaf + darea(i,j)*(1.-zfgr)
#else
            if (maskwater(i,j).eq.1) then
              totif = totif + darea(i,j)*h(i,j)
              totaf = totaf + darea(i,j)
            endif
#endif
          endif
        enddo
      enddo

c        "esl" calculations for effect on sea level.
c        Compute for total domain, and for WAIS only (allwais).

      totsl  = 0.
      totslw = 0.
      do j=1,ny
        do i=1,nx
c---
c old:
c         zeff is effect of ice and bedrock on sea level (m^3 liquid):
c         s.l. rise if all ice melted, including reduced contrib from
c         any ice below sea level:
c         zsl = sealev ! changed to sealev 11/12/12 (was = 0)
c         zeff = 0.
c         if (maskwater(i,j).eq.0) 
c    *      zeff = darea(i,j)
c    *             * (h(i,j)*(rhoice/rholiq) - max(0.,zsl-hb(i,j)))
c---
c new (04/10/14)...not used, because any change in hb inside domain
c due to asthenosphere must actually be compensated by opposite hb 
c change *outside* the domain ("conservation of asthen.", not in ELRA).
c         Based on total volume of any ice or water in domain. 
c         Global conserv:
c         dO1 + dO2 + dI*(rhoi/rhow)  = 0
c         where O1 = ocean volume outside domain, O2 = inside domain, 
c         and I = ice vol inside domain. Also: 
c         dO1 = A1*dSL, 
c         dO2 = dO2(ice+bed model, SL const) + A2*dSL
c         where A1 = ocean area outside domain, A2 = inside domain,
c         dO2(ice...) = change in ocean vol inside domain computed
c         from ice model alone (hw, accounting for ice and bed changes,
c         but with SL constant), and dSL is change in sea level.
c         nb: SL = sea-surface elevation from Earth center. Rearranging:
c         dSL = - (dO2(ice...) + dI*(rhoi/rhow)) / (A1+A2)
c         zeff = darea(i,j) * (h(i,j)*(rhoice/rholiq) + hw(i,j))
c---
c new (05/09/14)...simply ice over flotation, where grounded ice exists:
          zsl = sealev ! changed to sealev 11/12/12 (was = 0)
          zeff = 0.
          if (h(i,j).gt.0. .and. maskwater(i,j).eq.0)
     *      zeff = darea(i,j) 
     *           * (h(i,j)*(rhoice/rholiq) - max(0.,zsl-hb(i,j)))
c---
          totsl = totsl + zeff
#if defined (EISANTA)
          if (allwais(i,j).eq.1) totslw = totslw + zeff
#endif

        enddo
      enddo

c     totad = total area in domain, totao = tot ocn area in domain
c     totad = 0.
c     totao = 0.
c     do j=1,ny
c       do i=1,nx
c         totad = totad + darea(i,j)
c         if (maskwater(i,j).eq.1) totao = totao + darea(i,j)
c       enddo
c     enddo
c     zao = 0.7*(4.*pi*(radius**2)-totad) + totao   ! ~global ocean area
c     Or, simpler (and constant): 
      zao = 0.7* 4.*pi*(radius**2)                  ! ~global ocean area

      totsl  = totsl  / zao
      totslw = totslw / zao

      if (firstesl) then
        cumesl  = 0.
        cumeslw = 0.
        firstesl = .false.
      else
        cumesl  = cumesl  - (totsl-totslprev)
        cumeslw = cumeslw - (totslw-totslwprev)
      endif

      totslprev = totsl
      totslwprev = totslw

      write (iu, 
#if defined (ISMIP6CTL) || defined (ISMIP6SMB) || defined (ISMIP6BMB) || defined (ABUMIP)
     *       '(f10.3, 14e24.14, 6i24, 7e24.14)'  
#else
     *       '(i10, 
     *         2f7.3, f7.1, 
     *         4f7.2, f7.3, 
     *         f8.5,f8.3,f8.3,  f7.3,f7.3,f7.3, 
     *         6i12, 
     *         f8.1,
     *         3f10.4,
     *         3f10.4
     *        )'
#endif
     *      )
#if defined (ISMIP6CTL) || defined (ISMIP6SMB) || defined (ISMIP6BMB) || defined (ABUMIP)
     *  timeice, 
#else
     *  nint(timeice), 
#endif
#  if defined (FORCE34MYR)
     *  dtocn34m, ro18, sealev, 
#  else
     *  weirun, ro18, sealev, 
#  endif
     *  dtantann, dtantaus, dtantjan, dtseas, rco2,
     *  ecc, obl*180./pi, prec*180./pi, facice, facorb, facco2,
     *  nint(toti*1.e-9), nint ((toti-totif)*1.e-9), nint(totif*1.e-9),
     *  nint(tota*1.e-6), nint ((tota-totaf)*1.e-6), nint(totaf*1.e-6),
     *  toti/max(tota,.001),
     *  totsl-totslw, totslw, totsl, 
     *  cumesl-cumeslw, cumeslw, cumesl

      call flush (iu)

c====================
#if defined (EISANTA)
c====================
c     for indiv pts (hb,hw,h,hs):

c.......................
      if (firstpts) then
c.......................
        do k=1,npts
          zdistmin = 1.e20
          zxs = cos (alatpts(k)*pi/180.) * cos(alonpts(k)*pi/180.)
          zys = cos (alatpts(k)*pi/180.) * sin(alonpts(k)*pi/180.)
          zzs = sin (alatpts(k)*pi/180.)
          do j=1,ny
            do i=1,nx
              zx = cos (alatd(i,j)*pi/180.) * cos(alond(i,j)*pi/180.)
              zy = cos (alatd(i,j)*pi/180.) * sin(alond(i,j)*pi/180.)
              zz = sin (alatd(i,j)*pi/180.)
              zdist = sqrt ((zxs-zx)**2 + (zys-zy)**2 + (zzs-zz)**2)
              if (zdist.lt.zdistmin) then
                ipts(k) = i
                jpts(k) = j
                zdistmin = zdist
              endif
            enddo
          enddo
        enddo
        write (iupts,'(a10, a10, 20(6x,5x,a10,2f10.3,5x))')
     *     '      site',' (deg N,E)',
     *    ( cnamepts(k), 
     *      alatd(ipts(k),jpts(k)), alond(ipts(k),jpts(k)), 
     *      k=1,npts )
        write (iupts,'(/a10, a10, 20(6x,4a10)/)') 
     *     '      year','   s.l.(m)',
     *    ('     hb(m)','     hw(m)','      h(m)','     hs(m)',k=1,npts)

        firstpts = .false.
c..........
      endif
c..........

      write (iupts, '(i10, f10.2, 20(6x,f10.2,f10.3,2f10.2))')
     *  nint(timeice), sealev, 
     *  ( hb(ipts(k),jpts(k)),
     *    hw(ipts(k),jpts(k)),
     *    h (ipts(k),jpts(k)), 
     *    hs(ipts(k),jpts(k)), 
     *    k=1,npts )
      call flush (iupts)
c=====
#endif
c=====

      return
      end

c-----------------------------------------------------------------------

      subroutine icetabbud (h, hprev, maskwaterprev, bedge, fedge,
     *                      budgsnow, budgevap, budgrain, budgmelt,
     *                      baseperc, basefrml, 
     *                      oceanmelt, facemelt, cliffmelt, calvice, 
     *                      budgneg,
     *                      oceanmeltav, facemeltav, cliffmeltav, 
     *                      calviceav, runliqav, runfrozav,
     *                      toth0, tota0, toth, tota, totflow,
     *                      timeice, dt, iffirse, iflast, nyeartabbud)
 
c     Calculates total budget quantities, for 2-D runs.
c     Accumulates ocean related *av quantities (e.g. runliqav,runfrozav
c     for passing to Loveclim), reset in writehis.
c     If time, writes to ascii tabular file - toterr column tests 
c     ice mass conservation.

c     Some passed parameters:
c     tot[h0,a0]= ice vol, area at start of timestep dt (m3,m2,mascalc).
c     tot[h,a]  = ice vol, area at end of timestep dt (m3,m2,mascalc).
c     totflow   = total flow*dt across domain boundaries (m3, icedyn).
c     budgneg   = array(i,j) correction for negative ice (m/y, icedyn).
c     dt        = dynamics timestep (yr, = dtimeice in icectl).

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  h(nx,ny),         hprev(nx,ny),     maskwaterprev(nx,ny),  
     *  bedge(nx,ny),     fedge(nx,ny),
     *  budgsnow(nx,ny),  budgevap(nx,ny), 
     *  budgrain(nx,ny),  budgmelt(nx,ny),
     *  baseperc(nx,ny),  basefrml(nx,ny), 
     *  oceanmelt(nx,ny), facemelt(nx,ny), 
     *  cliffmelt(nx,ny), calvice(nx,ny),
     *  budgneg(nx,ny),
     *  oceanmeltav(nx,ny), facemeltav(nx,ny), cliffmeltav(nx,ny), 
     *  calviceav(nx,ny),   runliqav(nx,ny),   runfrozav(nx,ny)

c       Calculate total accumulation and ablation. Note that fractional
c       ice area fedge = 0 for ocean points with no ice. Also, 
c       facemelt,cliffmelt,calvice are zero for points with no ice at 
c       the beginning of timestep (before call to icedyn in icectl), 
c       not necessarily at end of timestep (after call to icedyn). 
c       So summed zacc and zabl are non-zero only over ice (and 
c       exposed land, if any), and after dividing by zz (=dt*ice area)
c       below, are the mean accumulation and ablation rates over ice.
c
c       Inclusion of budgneg in zabl cancels out the large 
c       overshoots of calvice,facemelt,cliffmelt at the edges of 
c       ice shelves/tidewater faces. (budgneg is set in icedyn, where 
c       any resulting negative ice thicknesses are reset to zero).

      totacc = 0.
      totabl = 0.
      totacc2 = 0.
      totabl2 = 0.

      totrain = 0.
      totmelt = 0.
      totrefr = 0.
      totfric = 0.
      totneg  = 0.

      totbasp = 0.
      totbasf = 0.
      totocme = 0.
      totcalv = 0.
      totface = 0.
      totclif = 0.

      totare2 = 0.

      do j=1,ny
        do i=1,nx

c         as in setting of budgall in icectl, i.e.:
          if (maskwaterprev(i,j).eq.0 .and. bedge(i,j).gt.0.) then
            zfedge = 1.
          else
            zfedge = fedge(i,j)
          endif
c         budgall(i,j) = 
c    *           (   budgsnow(i,j) - budgevap(i,j) 
c    *             + budgrain(i,j) - baseperc(i,j) + basefrml(i,j)
c    *             - oceanmelt(i,j)
c    *           ) * zfedge
c    *           - facemelt(i,j) - cliffmelt(i,j) - calvice(i,j)

          zacc =   (   budgsnow(i,j) 
     *               + budgrain(i,j) 
     *             ) * zfedge
          zabl =   (                 budgevap(i,j) 
     *                             + baseperc(i,j) - basefrml(i,j)
     *               + oceanmelt(i,j)
     *             ) * zfedge
     *             + facemelt(i,j) + cliffmelt(i,j) + calvice(i,j) 
     *             - budgneg(i,j)

          totacc  = totacc  + zacc*darea(i,j)*dt
          totabl  = totabl  + zabl*darea(i,j)*dt        

c         Adjust magnitude of [baseperc-basefrml+oceanmelt+calving
c         +facemelt+cliffmelt] for budgneg. budgneg was set in icedyn 
c         for overshoots producing -ve ice thickness, and combined 
c         with budgneg2 set in iceocean (setfhedge) for small h -> 0.
c         Reduce baseperc,etc, all in same proportion. Maintains 
c         conserv of mass (checked by zacc2, zabl2, totacc2, totabl2,
c         toterr2).
c
c         Accumulate these corrected fluxes in time, in *av variables, 
c         written to history, with net runliqav,runfrozav, for
c         forcing of ocean models (eg, Loveclim). 
c
c         Only do all this for points with some ice at beginning or end 
c         of timestep, or with non-zero ice flow into box or +ve 
c         snowfall (if zero ice flow in and no snowfall, zacc=zabl 
c         exactly, with zabl including budgneg). For other points, 
c         zero fluxes are accumulated in tot* and av* variables 
c         for this timestep.

c................................................
          if ( hprev(i,j).gt.0. .or. h(i,j).gt.0.
     *        .or. abs(zacc-zabl).gt.1.e-8 ) then
c................................................

            zrefr = budgrain(i,j) + budgmelt(i,j) - baseperc(i,j)
            zfric = max ( 0.,    baseperc(i,j) - basefrml(i,j)
     *                        - (budgrain(i,j) + budgmelt(i,j)) )
            totmelt = totmelt + budgmelt(i,j)*zfedge*darea(i,j)*dt
            totrefr = totrefr + zrefr        *zfedge*darea(i,j)*dt
            totfric = totfric + zfric        *zfedge*darea(i,j)*dt
            totneg  = totneg  + budgneg(i,j) *zfedge*darea(i,j)*dt

c           net ablative terms:
            zat = (   budgevap(i,j)
     *              + baseperc(i,j) - basefrml(i,j) + oceanmelt(i,j) )
     *            * zfedge 
     *            + facemelt(i,j) + cliffmelt(i,j) + calvice(i,j)

            if (zat.gt..001) then
c             account for budgneg by uniform factor x abl terms:
              zfac = 1. - budgneg(i,j)/zat
              zdsn = 0.
            else
c             if net abl terms are 0 or neg, account by + d(snowfall):
              zfac = 1.
              zdsn = budgneg(i,j)
            endif

c           apply corrections for budgneg: 
            zsnow = budgsnow(i,j) *zfedge + zdsn
            zrain = budgrain(i,j) *zfedge
            zevap = budgevap(i,j) *zfedge * zfac
            zbasp = baseperc(i,j) *zfedge * zfac
            zbfrm = basefrml(i,j) *zfedge * zfac
            zocme = oceanmelt(i,j)*zfedge * zfac
            zcalv = calvice(i,j)          * zfac
            zface = facemelt(i,j)         * zfac
            zclif = cliffmelt(i,j)        * zfac

            totrain = totrain + zrain*darea(i,j)*dt
            totbasp = totbasp - zbasp*darea(i,j)*dt
            totbasf = totbasf + zbfrm*darea(i,j)*dt
            totocme = totocme - zocme*darea(i,j)*dt 
            totcalv = totcalv - zcalv*darea(i,j)*dt 
            totface = totface - zface*darea(i,j)*dt 
            totclif = totclif - zclif*darea(i,j)*dt 

            zacc2 = zsnow + zrain
            zabl2 = zevap + zbasp - zbfrm + zocme
     *                    + zface + zclif + zcalv

            totacc2  = totacc2 + zacc2*darea(i,j)*dt
            totabl2  = totabl2 + zabl2*darea(i,j)*dt

            totare2  = totare2 + zfedge*darea(i,j)

            oceanmeltav(i,j) = oceanmeltav(i,j) + zocme*dt
            facemeltav(i,j)  = facemeltav(i,j)  + zface*dt
            cliffmeltav(i,j) = cliffmeltav(i,j) + zclif*dt
            calviceav(i,j)   = calviceav(i,j)   + zcalv*dt

            runliqav(i,j)    = runliqav(i,j)
     *                       + zbasp*dt
     *                       - zbfrm*dt
     *                       + zocme*dt
     *                       + zface*dt
            runfrozav(i,j)   = runfrozav(i,j)
     *                       + zcalv*dt
     *                       + zclif*dt
c..............
          endif
c..............

        enddo
      enddo

      if (nyeartabbud.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartabbud),dt)).lt.dt
     *       .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyeartabbud.eq.-1 )
     *   ) return 

      iu = iutabbud

c        Set total budget error (toterr) = rate of change of ice 
c        volume minus net local fluxes. toterr = zero if model 
c        conserves ice mass (written in last column "err" below).

      toterr  = (toth-toth0) - (totacc -totabl +totflow)
      toterr2 = (toth-toth0) - (totacc2-totabl2+totflow)

c        Normalize all flux quantities to m/y equiv over ice surfaces
c        (totare2, from here, does not include bare land. tota0, from 
c        mascalc, also does not, but doesn't account for fedge).

      zz2 = dt * totare2
      if (totare2.eq.0.) zz2 = totarea  ! prevent divide by 0 if no ice

      if (iffirse.eq.1) then
#if defined (ISMIP6CTL) || defined (ISMIP6SMB) || defined (ISMIP6BMB) || defined (ABUMIP)
        write (iu,'(a12, 25a24)')
#else
        write (iu,'(5a12, 9a12, 2a14, 10a12)')
#endif
     *    '        time',
     *    '  toti0(km3)', 
     *    '   toti(km3)', 
     *    '   tota(km2)',
     *    '  tota2(km2)',
     *    '  dtoti(m/y)', 
     *    '    acc(m/y)',
     *    '   acc2(m/y)',
     *    '    abl(m/y)',
     *    '   abl2(m/y)',
     *    '  ac-ab(m/y)', 
     *    ' ac2-ab2(m/y', 
     *    '    neg(m/y)',
     *    '   flow(m/y)',
     *    '      err(m/y)',
     *    '     err2(m/y)',
     *    '   rain(m/y)',
     *    '   melt(m/y)',
     *    '   basp(m/y)',
     *    '   basf(m/y)',
     *    '   refr(m/y)',
     *    '   fric(m/y)',
     *    '   ocme(m/y)',
     *    '   calv(m/y)',
     *    '   face(m/y)',
     *    '   clif(m/y)'
      endif
#if defined (ISMIP6CTL) || defined (ISMIP6SMB) || defined (ISMIP6BMB) || defined (ABUMIP)
      write (iu,'(f12.3, 25e24.14)')
     *  timeice, 
     *  toth0*1.e-9, 
     *  toth*1.e-9, 
     *  tota*1.e-6, 
     *  totare2*1.e-6, 
#else
      write (iu,'(5i12, 9f12.4, 2e14.4, 10f12.4)')
     *  nint(timeice), 
     *  nint(toth0*1.e-9), 
     *  nint(toth*1.e-9), 
     *  nint(tota*1.e-6), 
     *  nint(totare2*1.e-6), 
#endif
     *  (toth-toth0)/zz2,
     *  totacc/zz2, 
     *  totacc2/zz2, 
     *  totabl/zz2, 
     *  totabl2/zz2, 
     *  (totacc-totabl)/zz2, 
     *  (totacc2-totabl2)/zz2, 
     *  totneg/zz2, 
     *  totflow/zz2, 
     *  toterr/zz2,    ! e14.4
     *  toterr2/zz2,   ! e14.4
     *  totrain/zz2,
     *  totmelt/zz2,
     *  totbasp/zz2,
     *  totbasf/zz2,
     *  totrefr/zz2,
     *  totfric/zz2,
     *  totocme/zz2,
     *  totcalv/zz2,
     *  totface/zz2,
     *  totclif/zz2

      call flush (iu)

      return
      end

c-----------------------------------------------------------------------

c^^^^^^^^^^^^^^^^^^^^
#if defined (MELANGE)
c^^^^^^^^^^^^^^^^^^^^

      subroutine meltabbud (hmel, fmel, maskmel,
     *                      budgallmel,  budgsupmel, 
     *                      budgnegmel,  budgadjmel, 
     *                      budgredmel,  budgsmomel,
     *                      curbackmel,  curthetmel,
     *                      curhmel,     curlmel,
     *                      timeice, dt, iffirsr, iflast, nyeartabmel)
 
c     Calculates total budget quantities for melange (like icetabbud).
c     If time (nyeartabmel, writes to ascii tabular file.

#include <comicephys.h>
#include <comicegrid.h>

      dimension
     *  hmel(nx,ny),         fmel(nx,ny),         maskmel(nx,ny), 
     *  budgallmel(nx,ny),   budgsupmel(nx,ny),   
     *  budgnegmel(nx,ny),   budgadjmel(nx,ny),   
     *  budgredmel(nx,ny),   budgsmomel(nx,ny)

      logical firstmeltab
      data firstmeltab /.true./
      save firstmeltab

      save totbudall, totbudsup, totbudneg, totbudadj, 
     *     totbudred, totbudsmo,
     *     totbackmel, totthetmel, tothmel, totlmel,  
     *     totmprev,  totm2prev, timeprev

      dimension work(nx,ny), work2(nx,ny)

      if (iffirsr.eq.1) then
        totbudall = 0.
        totbudsup = 0.
        totbudneg = 0.
        totbudadj = 0.
        totbudred = 0.
        totbudsmo = 0.
        totbackmel = 0.
        totthetmel = 0.
        tothmel = 0.
        totlmel = 0.  
      endif

      totm = 0.
      totm2  = 0.
      tota   = 0.
      tota2  = 0.
      totag  = 0.
      totam  = 0.

      do j=1,ny
        do i=1,nx
          totbudall   = totbudall + budgallmel(i,j)*darea(i,j)*dt
          totbudsup   = totbudsup + budgsupmel(i,j)*darea(i,j)*dt 
          totbudneg   = totbudneg + budgnegmel(i,j)*darea(i,j)*dt 
          totbudadj   = totbudadj + budgadjmel(i,j)*darea(i,j)*dt 
          totbudred   = totbudred + budgredmel(i,j)*darea(i,j)*dt 
          totbudsmo   = totbudsmo + budgsmomel(i,j)*darea(i,j)*dt 

          totm  = totm  + hmel(i,j)*darea(i,j)
          totm2 = totm2 + fmel(i,j)*hmel(i,j)*darea(i,j)

          if (hmel(i,j).gt.0.) then 
            tota = tota  + darea(i,j)
            tota2 = tota2 + fmel(i,j)*darea(i,j)
            if (maskmel(i,j).eq.-1)  
     *        totag = totag + darea(i,j)          ! grounded melange
          endif
          if (maskmel(i,j).ne.0)
     *      totam = totam + darea(i,j)            ! for normalizing only
        enddo
      enddo

      totbackmel = totbackmel + curbackmel*dt
      totthetmel = totthetmel + curthetmel*dt
      tothmel = tothmel + curhmel*dt
      totlmel = totlmel + curlmel*dt

      if (iffirsr.eq.1) then
        totmprev  = totm
        totm2prev = totm2
        timeprev  = timeice
      endif

c---
      if (nyeartabmel.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartabmel),dt)).lt.dt
     *       .or. iflast.eq.1 .or. nyeartabmel.eq.-1 )
     *   ) return 
      iu = iutabmel
c---

c        Set total budget error (toterr) = rate of change of ice 
c        volume minus net local fluxes. toterr = zero if model 
c        conserves ice mass (written in last column "err" below).

      toterr  = (totm -totmprev)  
     *        - (   totbudall + totbudsup + totbudneg + totbudadj 
     *            + totbudred + totbudsmo )

c        Normalize all flux quantities to m/y equiv over melange 
c        surfaces (tota). Prevent divide by zero if no melange, etc.

      if (tota.gt.0.) then
        zz2 = tota                         ! area (current) with melange
      else if (totam.gt.0.) then 
        zz2 = totam                        ! area accessible for melange
      else 
        zz2 = dd0**2
      endif
      zz2 = zz2 * max(timeice-timeprev,.001)

      if ( firstmeltab .or. 
     *     abs(mod(timeice,100.)-nyeartabmel).lt.0.75*nyeartabmel ) then
        write (iu,'(/12a12, 5a14)')
     *    '        time',
     *    '   totm(km3)', 
     *    '   tota(km2)',
     *    '  totag(km2)',
     *    '       hm(m)',
     *    '    dhm(m/y)', 
     *    '   buda(m/y)',
     *    '    sup(m/y)',
     *    '    neg(m/y)',
     *    '    adj(m/y)',
     *    '    red(m/y)',
     *    '    smo(m/y)',
     *    '      err(m/y)',
     *    '    backmel(N)',
     *    '     1-thetmel',
     *    '       hmel(m)',
     *    '      lmel(km)'
        firstmeltab = .false.
      endif

      write (iu,'(f12.2,3i12,f12.3, 7f12.4, 5e14.4)')
     *  timeice, 
     *  nint(totm*1.e-9), 
     *  nint(tota*1.e-6), 
     *  nint(totag*1.e-6),
     *  totm/max(tota,.001),   ! hm

     *  (totm-totmprev)/zz2,   ! dhm 
     *  totbudall/zz2, 
     *  totbudsup/zz2, 
     *  totbudneg/zz2, 
     *  totbudadj/zz2, 
     *  totbudred/zz2, 
     *  totbudsmo/zz2,

     *  toterr/zz2,                                ! e14.4
     *  totbackmel/max(timeice-timeprev,.001),     ! e14.4
     *  1.-totthetmel/max(timeice-timeprev,.001),  ! e14.4
     *  tothmel/max(timeice-timeprev,.001),        ! e14.4
     *  .001*totlmel/max(timeice-timeprev,.001)    ! e14.4

      call flush (iu)

      totbudall = 0.
      totbudsup = 0.
      totbudneg = 0.
      totbudadj = 0.
      totbudred = 0.
      totbudsmo = 0.
      totbackmel = 0.
      totthetmel = 0.
      tothmel = 0.
      totlmel = 0.
      totmprev  = totm
      totm2prev = totm2
      timeprev  = timeice

      return
      end

c^^^^^
#endif
c^^^^^

c-----------------------------------------------------------------------


c{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
#if defined (EISANTA) && defined (TABTROUGH)
c{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{

      subroutine icetabtrough (h, hb, hw, maskwater, fedge, hbinit,
     *                         sealev, timeice, dt,
     *                         iffirse, iflast, nyeartot, nyeartabtro)

c     Similar to icetabwais.
 
#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny), hb(nx,ny), hw(nx,ny), maskwater(nx,ny), fedge(nx,ny),
     *  hbinit(nx,ny)

      dimension work(nx,ny)

#if defined (TRANSECTA) && defined (ZOOMGL)
      parameter (ntro = nzoom)   
#elif defined (TRANSECTA) && ! defined (ZOOMGL)
      parameter (ntro = nx)   
#else
      parameter (ntro = 10000)
#endif

#if defined (TRANSECTA)
      parameter (ntrans=1)
#else
      parameter (ntrans=4)
#endif

      dimension
     *  xtro(ntro,ntrans),     ytro(ntro,ntrans), 
     *  alondtro(ntro,ntrans), alatdtro(ntro,ntrans),
     *  disttro(ntro,ntrans),  tottro(ntrans),
     *  itro(ntro,2,ntrans),   witro(ntro,ntrans), 
     *  jtro(ntro,2,ntrans),   wjtro(ntro,ntrans),
     *  hbtro(ntro,ntrans),    dxtro(ntrans),
     *  distmodgl(ntrans),     mmodgl(ntrans)

      save xtro, ytro, alondtro, alatdtro, disttro,  tottro, 
     *     itro, witro, jtro, wjtro, hbtro, dxtro, distmodgl, mmodgl

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     transect names, and modern observed grounding-line lat,lon
c     (names are used to identify transect in settransect)   
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      character*8 cflag(6), cflagin
      data cflag 
     *  / 'PIG     ', 'MUR     ', 'ROSS    ', 'WEDDELL ', 'WILKES  ',  
     *    'SIPLE   ' /
      dimension alonmodgl(6),  alatmodgl(6) 
      data alonmodgl(1), alatmodgl(1) / -99.5, -75.3 /         ! PIG
      data alonmodgl(2), alatmodgl(2) /-111.3, -75.2 /         ! MURPHY
      data alonmodgl(3), alatmodgl(3) / -157., -83.0 /         ! ROSS
      data alonmodgl(4), alatmodgl(4) /  -60., -82.0 /         ! WEDDELL
      data alonmodgl(5), alatmodgl(5) /  153., -68.0 /         ! WILKES
      data alonmodgl(6), alatmodgl(6) / -154., -82.0 /         ! SIPLE
      save cflag, alonmodgl, alatmodgl, ntransa, ntransb
c+++++

      dimension 
     *  zhaf(ntro,ntrans), zfed(ntro,ntrans),
     *  distgl(ntrans),    distca(ntrans)

      logical firsttro
      data firsttro /.true./
      save firsttro

      if (nyeartabtro.eq.0) return
      if ( .not.
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartabtro),dt)).lt.dt
     *       .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyeartabtro.eq.-1)
     *   ) return

c=======================================================================
c        If first call, set zoom-resolution grid vars along entire 
c        trough transect, and indices ([w]itro, [w]jtro) into ice-model
c        grid. Transect vertices are hardcoded in settransect.
c        Also calculate distance to modern grounding line (distmodgl).
c=======================================================================

c~~~~~~~~~~~~~~~~~~~~~~~
      if (firsttro) then
c~~~~~~~~~~~~~~~~~~~~~~~

#if defined (TRANSECTA)
        cflagin = 'TRANSECTA'
        do m = 1,5
          if (cflagin.eq.cflag(m)) then
            ntransa = m
            ntransb = m
            goto 100
          endif
        enddo
        write (ioterm,'(/2a)') 
     *    '*** Error (icetabtrough): unknown transect name ', cflagin 
        stop
  100   continue
#else
       ntransa=1
       ntransb=4
#endif

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        do itrans=ntransa,ntransb
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#if defined (TRANSECTA)
#  if defined (ZOOMGL)
        dxtro(itrans) = dzoom
#  else
        dxtro(itrans) = dd0
#  endif
#endif

c         Set transect info (xtro,ytro,alondtro,alatdtro,disttro)
c         vs distance index (1:ntro)

        call settransect (xtro(1,itrans),     ytro(1,itrans), 
     *                    alondtro(1,itrans), alatdtro(1,itrans),
     *                    disttro(1,itrans),  tottro(itrans),
     *                    dxtro(itrans), ntro, cflag(itrans))

#if ! defined (TRANSECTA)
        dxtro(itrans) = tottro(itrans)/ntro
#endif

c          Compute ice-grid indices and weights along transect  

        do m=1,ntro
          zz = (xtro(m,itrans) - xoffa + 0.5*nx*dx0) / dx0
          itro(m,1,itrans) = nint(zz)
          witro(m,itrans)  = 1. - (zz - itro(m,1,itrans) + 0.5)
          itro(m,2,itrans) = itro(m,1,itrans) + 1
          itro(m,1,itrans) = max (1, min (nx, itro(m,1,itrans)))
          itro(m,2,itrans) = max (1, min (nx, itro(m,2,itrans)))
          witro(m,itrans)  = max (0., min (1., witro(m,itrans)))

          zz = (ytro(m,itrans) - yoffa + 0.5*ny*dy0) / dy0
          jtro(m,1,itrans) = nint(zz)
          wjtro(m,itrans)  = 1. - (zz - jtro(m,1,itrans) + 0.5)
          jtro(m,2,itrans) = jtro(m,1,itrans) + 1
          jtro(m,1,itrans) = max (1, min (ny, jtro(m,1,itrans)))
          jtro(m,2,itrans) = max (1, min (ny, jtro(m,2,itrans)))
          wjtro(m,itrans)  = max (0., min (1., wjtro(m,itrans)))
        enddo

c          Set transect distance and index to modern grounding line
c          (distmodgl, mmodgl)

        zlon = alonmodgl(itrans) 
        zlat = alatmodgl(itrans) 
        call ps_xy (zx, zy, zlat*pi/180., zlon*pi/180, 0., 0.)
        zdmin = 1.e20
        do m=1,ntro
          zd = (xtro(m,itrans)-zx)**2 + (ytro(m,itrans)-zy)**2
          if (zd.lt.zdmin) then
            mdmin = m 
            zdmin = zd
          endif
        enddo
        distmodgl(itrans) = disttro(mdmin,itrans)
        mmodgl(itrans) = mdmin

c          Set modern observed bedrock profile along transect 
c          (using hbinit, passed from initphys)  
 
        do m=1,ntro
          i1 = itro(m,1,itrans)
          i2 = itro(m,2,itrans)
          wi = witro(m,itrans)
          j1 = jtro(m,1,itrans)
          j2 = jtro(m,2,itrans)
          wj = wjtro(m,itrans)
          hbtro(m,itrans) =       wi *    wj  * hbinit(i1,j1)
     *                      +     wi *(1.-wj) * hbinit(i1,j2)
     *                      + (1.-wi)    *wj  * hbinit(i2,j1)
     *                      + (1.-wi)*(1.-wj) * hbinit(i2,j2)
        enddo

c          Diagnostic dump

        write (133,'(/2a,f12.3, a,i6/)') 
     *    cflag(itrans), '  distmodgl=', distmodgl(itrans)*.001,
     *    '  mmodgl=',mmodgl(itrans)

        call zero (work, nx*ny)
        do m=1,ntro
          zd = disttro(m,itrans)*.001
          work(itro(m,1,itrans),jtro(m,1,itrans)) = zd
          work(itro(m,1,itrans),jtro(m,2,itrans)) = zd
          work(itro(m,2,itrans),jtro(m,1,itrans)) = zd
          work(itro(m,2,itrans),jtro(m,2,itrans)) = zd
        enddo
c       show modern grounding line position as "*":
        m = mmodgl(itrans)
        work(itro(m,1,itrans),jtro(m,1,itrans)) = 1.e6
        work(itro(m,1,itrans),jtro(m,2,itrans)) = 1.e6
        work(itro(m,2,itrans),jtro(m,1,itrans)) = 1.e6
        work(itro(m,2,itrans),jtro(m,2,itrans)) = 1.e6
        call printmap (timeice, work, cflag(itrans) // ' disttro', 
     *                 100., 133, 0)

        call flush (133)

c>>>>>>>>>>>>
        enddo
c>>>>>>>>>>>>

c          Write modern bedrock profiles along transects

        iu = iutabtrohb
        write (iu,'(10(4x,a8,12x,5x))')  
     *    ( cflag(itrans),itrans=ntransa,ntransb )
        write (iu,'(10(a24,5x))')    
     *    ( '    dist(km)   hb_mod(m)', itrans=ntransa,ntransb )
        do m=1,ntro
          write (iu,'(10(f12.1,f12.2,5x))') 
     *       ( ( (disttro(m,itrans)-distmodgl(itrans))*.001,
     *           hbtro(m,itrans) ), 
     *         itrans=ntransa,ntransb  )
        enddo
        call flush (iu)

        firsttro = .false.
c~~~~~~~~~~
      endif
c~~~~~~~~~~


c============================================================
c       Find current intersections of model g.l. & shelf edge 
c       with trough transect (distgl, distca)
c============================================================

c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        do itrans=ntransa,ntransb
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c        Look along trough transect, bilin. interp. heigh-over-flotation
c        (zhaf, m liq.equiv.), and ice-shelf fractional cover (zfed)

      do m=1,ntro
        i1 = itro(m,1,itrans)
        i2 = itro(m,2,itrans)
        wi = witro(m,itrans)
        j1 = jtro(m,1,itrans)
        j2 = jtro(m,2,itrans)
        wj = wjtro(m,itrans)
        zhaf(m,itrans) =
     *                wi *    wj  * (rhor*h(i1,j1) - (sealev-hb(i1,j1)))
     *          +     wi *(1.-wj) * (rhor*h(i1,j2) - (sealev-hb(i1,j2)))
     *          + (1.-wi)    *wj  * (rhor*h(i2,j1) - (sealev-hb(i2,j1)))
     *          + (1.-wi)*(1.-wj) * (rhor*h(i2,j2) - (sealev-hb(i2,j2)))

        zfed(m,itrans) =
     *                wi *    wj  * fedge(i1,j1) 
     *          +     wi *(1.-wj) * fedge(i1,j2)
     *          + (1.-wi)*    wj  * fedge(i2,j1) 
     *          + (1.-wi)*(1.-wj) * fedge(i2,j2)
      enddo

      distgl(itrans) = -1.
      distca(itrans) = -1.
c     In settransect, trough runs from inland to ocean (1 to ntro).  
c     Here, do the same (1 to ntro), find first non-grounded pt 
c     (distgl) and first open-ocean point (distca)
      do m=1,ntro
        if (zhaf(m,itrans).lt.0. .and. distgl(itrans).eq.-1.) then
           distgl(itrans) = disttro(m,itrans)
        endif
        if (zhaf(m,itrans).lt.0. .and. zfed(m,itrans).lt.0.5 
     *      .and. distca(itrans).eq.-1.) then
           distca(itrans) = disttro(m,itrans)
        endif
      enddo
c     reached edge of domain wiothout finding gl or ocean, set to edge:
      if (distgl(itrans).eq.-1.) distgl(itrans) = disttro(ntro,itrans) 
      if (distca(itrans).eq.-1.) distca(itrans) = disttro(ntro,itrans) 

c<<<<<<<<<<
      enddo
c<<<<<<<<<<


c=====================================================================
c        Write tabular output (1 line per current time, all transects)
c=====================================================================

      iu = iutabtro

c        Write header lines first call

c...........................
      if (iffirse.eq.1) then
c...........................
        nwrite = nint(float(nyeartot)/float(nyeartabtro)) + 1
        write (iu,'(a,i8)') 'nwrite=', nwrite 
        do itrans = ntransa,ntransb
          write (iu,'(a8, 2(a,i8), a,f8.2, a,i8)') 
     *      cflag(itrans),
     *      '  ntro=',          ntro,
     *      '  tottro(km)=',    nint(tottro(itrans)*.001),
     *      '  dxtro(km)=',     dxtro(itrans)*.001,
     *      '  distmodgl(km)=', nint(distmodgl(itrans)*.001)
        enddo

        write (iu,'(12x,10(6x,a8,10x,5x))') 
     *    ( cflag(itrans), itrans=ntransa,ntransb ) 

        write (iu,'(a,10(a,5x))') 
     *      '        time',
     *    ( '      distgl      distca',
     *      itrans=ntransa,ntransb )
c..........
      endif
c..........

c        Distances of gl and calving edge are km downstream from modern
c        gl (1992 for PIG?)

      write (iu, '(i12, 10(2f12.3,5x))')
     *  nint(timeice), 
     *  ( (distgl(itrans)-distmodgl(itrans))*.001,  
     *    (distca(itrans)-distmodgl(itrans))*.001, 
     *    itrans=ntransa,ntransb )

      call flush (iu)

      return
      end

c{{{{{
#endif
c{{{{{

c-----------------------------------------------------------------------

c#############################################
#if defined (EISANTA) && ! defined (TRANSECTA)
c#############################################

      subroutine icetabwais (h, hb, hw, maskwater, 
     *                       ub, vb,
     *                       sealev, timeice, dt, weirun, 
     *                       iffirse, iflast, nyeartot, nyeartabwais)
 
c     Calculates WAIS-specific quantities:

c     local conditions at specified points:
c       Whether pts "a" (alona,alata) is gd/shelf/ocn (istata=1,2,3).
c       Distances from pts "a" (alona,alata) to nearest g.l.,shelf edge.
c       Provenance of pts "a", following trajectories back to shore/g.l.
c     all WAIS:
c       WAIS-only ice vols, areas.
c       Distances along Siple->Ross Sea transect to g.l., shelf edge 

c     Writes these (and h,hb,hw of pt "a") to ascii tabular file.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny), hb(nx,ny), hw(nx,ny), maskwater(nx,ny)

      dimension
     *  ub(0:nxp,0:nyp), vb(0:nxp,0:nyp)

c     For local points "a":
      logical firsta
      save firsta 
      data firsta /.true./
#if defined (ALLTABWAIS)
      parameter (npa=1)
#else
      parameter (npa=3)
#endif
      dimension alona(npa), alata(npa), ia(npa), ja(npa), istata(npa), 
     *          glmin(npa), edmin(npa), alonp(npa), alatp(npa)
      save alona, alata, ia, ja
#if defined (ALLTABWAIS)
      data alata / -77.8894417 /                                ! AND-1B
      data alona / 167.0893282 /                                ! AND-1B
#else
      data alata / -77.8894417,  -77.8894417,   -81.0 /
      data alona / 167.0893282,  167.0893282,  -177.0 /
#endif

c     For defining all-WAIS region: 
      logical firstwais
      save firstwais
      data firstwais /.true./
      parameter (nwais=10)
      dimension 
     *  alonwais(nwais), alatwais(nwais), xwais(nwais), ywais(nwais),
     *  allwais(nx,ny)
      data alatwais 
     * / -45., -45., -82.,  -85.,  -86., -82., -75., -72., -45.,  -45./
      data alonwais 
     *  /-80., -40., -40.,  -70., -140., 163., 164., 170., 170., -120./
      save allwais

c     For Siple-> Ross Sea transect:
      logical firstsip
      save firstsip
      data firstsip /.true./
      parameter (nsipmax = 1000) 
      dimension xsip(nsipmax),    ysip(nsipmax), 
     *          alonsip(nsipmax), alatsip(nsipmax),  
     *          distsipa(nsipmax),distsip(nsipmax),  
     *          isip(nsipmax),    jsip(nsipmax)
c    *         ,work(nx,ny)                        ! for diagnostic dump
      save isip, jsip, nsip, distsip

      if (nyeartabwais.eq.0) return
      if ( .not. 
     *     ( mod(abs(timeice)+0.5*dt,max(float(nyeartabwais),dt)).lt.dt
     *       .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyeartabwais.eq.-1)
     *   ) return 

c        Convert lon,lat of pts "a" to i,j polar coords

c=====================
      if (firsta) then
c=====================
        do mp=1,npa
          call ps_xy (zxa, zya, alata(mp)*pi/180., alona(mp)*pi/180, 
     *                xoffa, yoffa)
          ia(mp) = (zxa + 0.5*nx*dx0)/dx0 + 1.0001
#if ! defined (NESTING)
c         shift AND-1B location 1 grid point to east:
          if ( nint(alata(mp)).eq.-78 .and. nint(alona(mp)).eq.167
#  if ! defined (ALLTABWAIS)
     *         .and. mp.eq.2
#  endif
     *       ) ia(mp) = ia(mp) -1
#endif
          ja(mp) = (zya + 0.5*ny*dy0)/dy0 + 1.0001

          ia(mp) = max (1, min (nx, ia(mp)))
          ja(mp) = max (1, min (ny, ja(mp)))
        enddo
        firsta = .false.
c==========
      endif
c==========

c===================================================================
c        Determine if pts "a" are grounded, sub-shelf, or open ocean
c===================================================================

      do mp=1,npa
        if (maskwater(ia(mp),ja(mp)).eq.0) then
          istata(mp) = 1                        ! grounded
        else
          if (h(ia(mp),ja(mp)).gt.0.) then
            istata(mp) = 2                      ! sub-shelf
          else
            istata(mp) = 3                      ! open ocean
          endif
        endif
      enddo

c=====================================================================
c        Calculate distances from pt "a" to closest grounding line and
c        to closest shelf edge (glmin, edmin)
c=====================================================================

#if defined (ALLTABWAIS)
      do mp=1,npa
        call distgled (glmin(mp), edmin(mp), istata(mp), ia(mp), ja(mp),
     *                 h, maskwater)
      enddo
#else
      do mp=1,npa
        glmin(mp) = 0.
        edmin(mp) = 0.
      enddo
#endif

c=====================================================================
c        Calculate "provenance" of pt "a" (alonp,alatp), tracking 
c        lagranginan path upstream to first-encountered grounding line  
c=====================================================================

#if defined (ALLTABWAIS)
      do mp=1,npa
        call provenance (alonp(mp), alatp(mp), istata(mp),ia(mp),ja(mp),
     *                   h, maskwater, ub, vb)
      enddo
#else
      do mp=1,npa
        alonp(mp) = 0.
        alatp(mp) = 0.
      enddo
#endif

c=============================================================
c        Calculate all-Ant and all-WAIS ice volumes, ice areas
c=============================================================

      if (firstwais) then
#if ! defined (NESTING)
c       convert vertices to polar stereo x,y coords:
        do m=1,nwais
          call ps_xy (xwais(m), ywais(m), 
     *                alatwais(m)*pi/180., alonwais(m)*pi/180, 
     *                0., 0.)
        enddo

c       set allwais array (1=inside WAIS polygon, 0=outside):
        call zero (allwais, nx*ny)
        call polygon_wais (allwais, xh,yh,nx,ny, xwais,ywais, nwais,1.)

c       diagnostic dump:
c       call printmap (timeice,allwais,'allwais', 1., 133, 0)
#else
        call resetr (allwais, nx*ny, 1.)
#endif
        firstwais = .false.
      endif

      toti      = 0.
      tota      = 0.
      totif     = 0.
      totaf     = 0.
      totiwais  = 0.
      totawais  = 0.
      totifwais = 0.
      totafwais = 0.

      do j=1,ny
        do i=1,nx

c         all Antarctic:
          if (h(i,j).ne.0.) then  
            toti = toti + darea(i,j)*h(i,j)
            tota = tota + darea(i,j)
            if (maskwater(i,j).eq.1) then
              totif = totif + darea(i,j)*h(i,j)
              totaf = totaf + darea(i,j)
            endif
          endif

c         all WAIS:
#if defined (NESTING)
          if (.true.) then
#else
c         if ( (alond(i,j).lt.-40. .or. alond(i,j).gt.170.) .and. 
c    *          alatd(i,j).gt.-85. ) then
          if (allwais(i,j).eq.1.) then
#endif
            if (h(i,j).ne.0.) then  
              totiwais = totiwais + darea(i,j)*h(i,j)
              totawais = totawais + darea(i,j)
              if (maskwater(i,j).eq.1) then
                totifwais = totifwais + darea(i,j)*h(i,j)
                totafwais = totafwais + darea(i,j)
              endif
            endif
          endif

        enddo
      enddo

c======================================================================
c       Find intersections of RIS g.l. & shelf edge with Siple transect
c======================================================================

c         Set transect coords (transect vertices are hard-coded in 
c         settransect). Also, convert lon,lats to i,j polar st. indices

      if (firstsip) then
        call settransect (xsip, ysip, alonsip, alatsip,
     *                    distsipa, totlensip,
     *                    dd0, nsipmax, 'SIPLE')
        nsip = 0
        isipprev = 0
        jsipprev = 0
        do i = 1,nsipmax
          isipz = (xsip(i) - xoffa + 0.5*nx*dx0)/dx0 + 1.0001
          jsipz = (ysip(i) - yoffa + 0.5*ny*dy0)/dy0 + 1.0001
          if (isipz.ge.1 .and. isipz.le.nx .and. 
     *        jsipz.ge.1 .and. jsipz.le.ny .and. 
     *        (isipz.ne.isipprev .or. jsipz.ne.jsipprev) ) then
            nsip = nsip + 1
            isip(nsip) = isipz
            jsip(nsip) = jsipz
            distsip(nsip) = distsipa(i)
c           distsip(nsip) = (i-1.)*(totlensip/nsipmax)
            isipprev = isipz
            jsipprev = jsipz
          endif
        enddo

c       diagnostic dump:
c       call zero (work, nx*ny)
c       do is=1,nsip
c          work(isip(is),jsip(is)) =  distsip(is)*1.e-3
c       enddo
c       call printmap (timeice,work,'distsip', 100., 133, 0)

        firstsip = .false.
      endif

c        Look along transect (runs from WAIS interior into Ross Sea),
c        setting first-encountered ocean (gl) point and shelf-edge point

      distsip_gl = -1.
      distsip_ed = -1.
      do is=1,nsip
        ii = isip(is)
        jj = jsip(is)
        if (maskwater(ii,jj).eq.1 .and. distsip_gl.eq.-1.)
     *    distsip_gl = distsip(is)

        if (maskwater(ii,jj).eq.1 .and. h(ii,jj).eq.0. .and. 
     *      distsip_ed.eq.-1.) 
     *    distsip_ed = distsip(is)
      enddo

c============================
c        Write tabular output
c============================

      iut = iutabwais  ! iu used above as a variable

c        Write header line(s) first call, if not a restart

c...........................
      if (iffirse.eq.1) then
c...........................
        nwrite = nint(float(nyeartot)/float(nyeartabwais)) + 1
        write (iut,'(2(a,i8))') 
     *    'nwritetab =',nwrite, '    npa =',npa
        do mp=1,npa
          write (iut,'( 2 (2(a,i1,a,f8.3,2x),5x) )')
     *      'alona_', mp,' = ',alona(mp),  
     *      'alata_', mp,' = ',alata(mp),
     *      'alonda_',mp,' = ',alond(ia(mp),ja(mp)),
     *      'alatda_',mp,' = ',alatd(ia(mp),ja(mp))  
        enddo

        write (iut,'(6a,10(5a))')
     *    '      time weirun sealev',

     *    '   toti(km3)  totig(km3)  totif(km3)',
     *    '   tota(km2)  totag(km2)  totaf(km2)    h(m)',

     *    '   ti_w(km3)  tig_w(km3)  tif_w(km3)',
     *    '   ta_w(km2)  tag_w(km2)  taf_w(km2)  h_w(m)',

     *    '  dsipgl  dsiped',

     *   ('      ',
     *    '     h_a    hb_a    hw_a',
     *    '   g/s/o',
     *    '   glmin   edmin',
     *    '   alonp   alatp',
     *    mp=1,npa)
c..........
      endif
c..........

      write (iut,
     *       '(i10, f7.3, f7.1, 

     *         6i12, f8.1,

     *         6i12, f8.1,

     *         2f8.1,

     *         10(6x,
     *            3f8.1,
     *            i8,
     *            2f8.1,
     *            2f8.2)
     *        )')
     *  nint(timeice), weirun, sealev, 

     *  nint(toti*1.e-9), nint ((toti-totif)*1.e-9), nint(totif*1.e-9),
     *  nint(tota*1.e-6), nint ((tota-totaf)*1.e-6), nint(totaf*1.e-6),
     *  toti/max(tota,.001),

     *  nint(totiwais*1.e-9), nint ((totiwais-totifwais)*1.e-9), 
     *    nint(totifwais*1.e-9),
     *  nint(totawais*1.e-6), nint ((totawais-totafwais)*1.e-6), 
     *    nint(totafwais*1.e-6),
     *  totiwais/max(totawais,.001),

     *  distsip_gl*1.e-3, distsip_ed*1.e-3, ! transect dist-> gl,sh.edge

     * (h(ia(mp),ja(mp)), hb(ia(mp),ja(mp)), hw(ia(mp),ja(mp)),
     *  istata(mp),                         ! gl/sh/oc (1/2/3)
     *  glmin(mp)*1.e-3, edmin(mp)*1.e-3,   ! upstm to nearest gl,sh.ed
     *  alonp(mp), alatp(mp),               ! provenance (or 0)
     *  mp=1,npa)

      call flush (iut)

      return
      end

c-----------------------------------------------------------------------

      subroutine distgled (glmin, edmin, istata, ia, ja, h, maskwater)

c     Calculate distances (glmin, edmin) from pt (ia,ja) to closest 
c     grounding/shore-line, and to nearest shelf edge

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny), maskwater(nx,ny)

      glmin = 1.e20
      edmin = 1.e20

      do j=1,ny-1
        do i=1,nx-1
          if (maskwater(i,j).ne.maskwater(i+1,j)) then
            glmin = min ( glmin, (i+.5-ia)**2 + (j-ja)**2 ) 
          endif
          if (maskwater(i,j).ne.maskwater(i,j+1)) then
            glmin = min ( glmin, (i-ia)**2 + (j+.5-ja)**2 ) 
          endif

          if ( maskwater(i,j).eq.1 .and. maskwater(i+1,j).eq.1 .and.
     *         h(i,j)*h(i+1,j).eq.0. .and. h(i,j)+h(i+1,j).gt.0. ) then
            edmin = min ( edmin, (i+.5-ia)**2 + (j-ja)**2 ) 
          endif
          if ( maskwater(i,j).eq.1 .and. maskwater(i,j+1).eq.1 .and.
     *         h(i,j)*h(i,j+1).eq.0. .and. h(i,j)+h(i,j+1).gt.0. ) then
            edmin = min ( edmin, (i-ia)**2 + (j+.5-ja)**2 ) 
          endif
        enddo
      enddo

c        Convert to meters, and set negative if appropriate, 
c        i.e., *upstream* distances from point "a"

      if (glmin.ne.1.e20) then
        glmin = dd0*sqrt(glmin)
        if (istata.eq.1) glmin = -glmin                     ! grounded
      else
        glmin = 99999900.                   ! matched for output format
      endif

      if (edmin.ne.1.e20) then
        edmin = dd0*sqrt(edmin)
        if (istata.eq.1 .or. istata.eq.2) edmin = -edmin    ! gd or sh
      else
        edmin = 99999900.                   ! matched for output format
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine provenance (alonp, alatp, istata, ia, ja, 
     *                       h, maskwater, ub, vb)

c     Find "provenance" (alonp,alatp) of point (ia,ja),
c     by using <ub,vb> field, bilinearly interpolated, to find 
c     lagrangian path from (a,ja) (if within ice shelf,istata=2) 
c     to first-encountered "grounding line" (where path first enters 
c     a grounded h-grid box).
c
c     Lagrangian path coords are (xp,yp) (km's from domain corner,
c     not the same as x,y from South Pole used elsewhere).
c
c     Scale interpolated velocity at each iteration to move the
c     current location a fraction "fracprov" of a grid box.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),        maskwater(nx,ny),
     *  ub(0:nxp,0:nyp), vb(0:nxp,0:nyp)

      parameter (fracprov = 0.2)
      parameter (niterprov = 100000)

c++++++++++++++++++++++++++
      if (istata.eq.2) then
c++++++++++++++++++++++++++

        xp = dd0*(ia-.5)
        yp = dd0*(ja-.5)

c>>>>>>>>>>>>>>>>>>>>>>>>>>
        do iter=1,niterprov
c>>>>>>>>>>>>>>>>>>>>>>>>>>

c            Interpolate ub,vb to current lagrangian point (xp,yp) 
c---
          iu = xp/dd0
          ziu = xp/dd0 - iu 

          ju = nint(yp/dd0)
          zju = yp/dd0 + 0.5 - ju

          zu= (1.-ziu)*(1.-zju)*ub(iu,ju)   + ziu*(1.-zju)*ub(iu+1,ju) 
     *      + (1.-ziu)*(   zju)*ub(iu,ju+1) + ziu*(   zju)*ub(iu+1,ju+1)
c---
          iv = nint(xp/dd0)
          ziv = xp/dd0 + 0.5 - iv 

          jv = yp/dd0
          zjv = yp/dd0 - jv

          zv= (1.-ziv)*(1.-zjv)*vb(iv,jv)   + ziv*(1.-zjv)*vb(iv+1,jv) 
     *      + (1.-ziv)*(   zjv)*vb(iv,jv+1) + ziv*(   zjv)*vb(iv+1,jv+1)
c---

          zsp = sqrt (zu**2 + zv**2)
          if (zsp.lt.1.e-20) then
            write (6,'(a,i10,2i6,e15.8)')
     *        'provenance error 1: zero speed: time,iu,ju,zsp=', 
     *        nint(timeice), iu, ju, zsp
c           stop
            alonp = -1.
            alatp = -1.
            goto 1000
          endif

c            Move lagrangian point (scaling velocity so shift is always 
c            fracprov*dd0), and find nearest h-grid indices ip,jp
          
          xp = xp - (zu/zsp) * fracprov*dd0
          yp = yp - (zv/zsp) * fracprov*dd0

          ip = nint(xp/dd0)
          zip = xp/dd0 + 0.5 - ip
          jp = nint(yp/dd0)
          zjp = yp/dd0 + 0.5 - jp

          if (ip.le.1 .or. ip.ge.nx .or. jp.le.1 .or. jp.ge.ny) then
            write (6,'(a,i10,2i6)')
     *        'provenance error 2: reached domain edge: time,ip,jp=',
     *        nint(timeice), ip, jp
c           stop
            alonp = -2.
            alatp = -2.
            goto 1000
          endif

c........................................
          if (maskwater(ip,jp).eq.0) then
c........................................

c           If reached grounded ice or ice-free ground, 
c           find interpolated lon,lat of lagrangian locn, and skip out:

            if (abs(alond(ip,jp)).lt.90.) then
c             this works across greenwich, not across dateline:
              alonp =  (1.-zip)*(1.-zjp)*alond(ip,  jp)
     *              +      zip *(1.-zjp)*alond(ip+1,jp)
     *              +  (1.-zip)*(   zjp)*alond(ip  ,jp+1)
     *              +  (   zip)*(   zjp)*alond(ip+1,jp+1)
            else
c             this works across dateline, not across greenwich:
              alonp =  (1.-zip)*(1.-zjp)*mod(alond(ip,  jp)  +360.,360.)
     *              +      zip *(1.-zjp)*mod(alond(ip+1,jp)  +360.,360.)
     *              +  (1.-zip)*(   zjp)*mod(alond(ip  ,jp+1)+360.,360.)
     *              +  (   zip)*(   zjp)*mod(alond(ip+1,jp+1)+360.,360.)
              alonp = mod (alonp+180., 360.) - 180.
            endif

            alatp =  (1.-zip)*(1.-zjp)*alatd(ip,  jp)
     *            +      zip *(1.-zjp)*alatd(ip+1,jp)
     *            +  (1.-zip)*(   zjp)*alatd(ip  ,jp+1)
     *            +  (   zip)*(   zjp)*alatd(ip+1,jp+1)

            goto 1000

c..................................................................
          else if (maskwater(ip,jp).eq.1 .and. h(ip,jp).gt.0.) then
c..................................................................

c           If still in ice shelf, continue iteration (do nothing here)

c..................................................................
          else if (maskwater(ip,jp).eq.1 .and. h(ip,jp).eq.0.) then
c..................................................................

c           If reached edge of ice shelf, error

            write (6,'(a,i10,2i6)')
     *        'provenance error 3: reached shelf edge: time,ip,jp=', 
     *        nint(timeice), ip, jp
c           stop
            alonp = -3.
            alatp = -3.
            goto 1000

c..............
          endif
c..............

          if (iter.eq.niterprov) then
            write (6,'(a,i10,2i6)')
     *        'provenance error 4: g.l. not found: time,ip,jp=', 
     *        nint(timeice), ip, jp
c           stop
            alatp = -4.
            alonp = -4.
            goto 1000
          endif

c>>>>>>>>>>>>
        enddo
c>>>>>>>>>>>>
 1000   continue  ! successfully found grounding line

c+++++++++
      else
c+++++++++

        alatp = 0.
        alonp = 0.

c++++++++++
      endif
c++++++++++

      return
      end

c#####
#endif
c#####

c-----------------------------------------------------------------------

      subroutine polygon_wais (arr, xlon, xlat, mlon, mlat,
     *                         xvert, yvert, nvert, valin)
 
c Identifies points in a regular 2-D mlon-by-mlat grid that lie
c within a polygon defined by nvert vertices (xvert, yvert).
c The method is to test one regular grid point at a time, by calculating
c the angle swept out by the vector originating at the grid point
c and moving from one vertex to the next (anticlockwise = positive).
c Then the sum of all such arcs will be zero for points ouside the
c polygon, and +/- 360 degrees for points inside.
 
c arr = regular-grid array whose "inside" pts are set to valin (mod)
c xlon = coords of regular grid (supplied)
c xlat = coords of regular grid (supplied)
c mlon = x-dimension of regular grid (supplied)
c mlat = y-dimension of regular grid (supplied)
c xvert = vertex coords (supplied)
c yvert = vertex coords (supplied)
c nvert = number of vertices (supplied)
c valin = value for "inside" points of arr (supplied)
 
      dimension arr (mlon,mlat), xlon(mlon,mlat), xlat(mlon,mlat),
     *          xvert(nvert), yvert(nvert)
 
      parameter (pi = 3.14159265358979)
c     degrees per radian
      parameter (dpr = 180./pi)
 
c     eps should be negligible but greater than machine precision
c     parameter (eps = 1.e-5)  ! if single precision
      parameter (eps = 1.e-10)
 
c-----
 
c        Loop over all points in the regular grid
 
      do 100 j=1,mlat
        do 102 i=1,mlon
          ytest = xlat(i,j)
          xtest = xlon(i,j)
 
c            Loop over vertex pairs (including nvert-to-1),
c            compute arc increment, and accumulate in angle.
 
          angle = 0.
          do 110 k=1,nvert
            kp = k+1
            if (k.eq.nvert) kp = 1
            xv1 = xvert(k)
            yv1 = yvert(k)
            xv2 = xvert(kp)
            yv2 = yvert(kp)
 
            if ( abs(yv1-ytest) + abs(xv1-xtest) .gt. eps ) then
              ang1 = dpr * atan2 (yv1-ytest, xv1-xtest)
            else
              ang1 = 0.
            endif
 
            if ( abs(yv2-ytest) + abs(xv2-xtest) .gt. eps ) then
              ang2 = dpr * atan2 (yv2-ytest, xv2-xtest)
            else
              ang2 = 0.
            endif
 
            dang = ang2 - ang1
 
c           If grid pt. on a vertex line, set to "in" and skip out
            if (abs(abs(dang)-180.).lt.eps) then
              arr(i,j) = valin
              goto 102
            endif
 
c           Increment total angle (I *think* the following logic is
c           independent of the range returned by atan2)
            if      (dang.le. 360 .and. dang .ge.180.) then
              angle = angle + dang - 360.
            else if (dang.lt. 180 .and. dang .ge. -180.) then
              angle = angle + dang
            else if (dang.lt.-180 .and. dang .ge.-360.) then
              angle = angle + dang + 360.
            endif
  110     continue
 
c         write(*,120) xlat(i,j), xlon(i,j), angle
c 120     format('xlat,xlon,angle=',3f10.4)
 
c            If angle is (almost) +/- 360 degrees, set to "in"
 
          if (abs(abs(angle)-360.) .lt. eps) arr(i,j) = valin
  102   continue
  100 continue
 
      return
      end

c-----------------------------------------------------------------------

      subroutine writehis (h, hs, hb, hbsd, t, tracer, fedge,
     *                     budgall, budgsnow, budgrain, budgmelt,
     *                     baseperc, basefrml, 
     *                     oceanmelt,msector,facemelt,cliffmelt,calvice,
#if defined (RACMOCLIM)
     *                     runofa,
#endif
     *                     tsurf, tsurfi, tmon,
     *                     heati, heath, heatb,
     *                     w, wa, sedim, tsed, wsed, heats, dgeoid,
     *                     quarryrate, quarryacc, quarrycoef,
     *                     pelagrate, pelagacc,
     *                     topbed, topbedeq, deflect, 
     *                     tbed, geoflux,
     *                     hw, tw, maskwater, crhmel, enhanceinv,
     *                     oceanmeltav, facemeltav, cliffmeltav,
     *                     calviceav, runliqav, runfrozav, timefluxprev,
     *                     u, v, ua, va, ub, vb, uadv, vadv,
     *                     fluxschu, fluxschv,
     *                     thetau, thetav,
#if defined (MELANGE)
     *                     hmel, maskmel, budgallmel, budgsupmel,
     *                     uamel, vamel, 
     *                     bpumel, bpvmel, 
     *                     thetaumel, thetavmel,
#endif                 
     *                     sealev, timeice, dt, nyearstart,
     *                     iffirse, iflast, nyearhis, hislist, nhislist)

c     Write to history file

#include <comicephys.h>
#include <comicegrid.h>
#include <comiceplot.h>

      dimension 
     *  h(nx,ny),            hs(nx,ny),     
     *  hb(nx,ny),           hbsd(nx,ny),
     *  t(nx,ny,0:nlevp),    tracer(nx,ny,0:nlevp,ntrace),
     *  fedge(nx,ny),
     *  budgall(nx,ny),      budgsnow(nx,ny),    
     *  budgrain(nx,ny),     budgmelt(nx,ny),
     *  baseperc(nx,ny),     basefrml(nx,ny),     
     *  oceanmelt(nx,ny),    msector(nx,ny),
     *  facemelt(nx,ny),     cliffmelt(nx,ny),   calvice(nx,ny),  
#if defined (RACMOCLIM)
     *  runofa(nx,ny),
#endif
     *  tsurf(nx,ny),        tsurfi(nx,ny),      tmon(nmon,nx,ny),
     *  heati(nx,ny,nlev),   heath(nx,ny),       heatb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),
     *  sedim(nx,ny),        tsed(nx,ny,nsed),   wsed(nx,ny,nsed),
     *  quarryrate(nx,ny),   quarryacc(nx,ny),   quarrycoef(nx,ny), 
     *  pelagrate(nx,ny),    pelagacc(nx,ny),
     *  heats(nx,ny,nsed),   dgeoid(nx,ny),
     *  topbed(nx,ny),       topbedeq(nx,ny),
     *  deflect(nx,ny),      tbed(nx,ny,nbed),   geoflux(nx,ny),
     *  hw(nx,ny),           tw(nx,ny),          maskwater(nx,ny),    
     *  crhmel(nx,ny),       enhanceinv(nx,ny),
     *  oceanmeltav(nx,ny),  facemeltav(nx,ny),  cliffmeltav(nx,ny),
     *  calviceav(nx,ny),    runliqav(nx,ny),    runfrozav(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  fluxschu(0:nxp,0:nyp),  fluxschv(0:nxp,0:nyp),
     *  thetau(nx,ny),          thetav(nx,ny)

#if defined (MELANGE)
     * ,hmel(nx,ny),            maskmel(nx,ny),      
     *  budgallmel(nx,ny),      budgsupmel(nx,ny),
     *  uamel(0:nxp,0:nyp),     vamel(0:nxp,0:nyp),
     *  bpumel(0:nxp,0:nyp),    bpvmel(0:nxp,0:nyp),
     *  thetaumel(nx,ny),       thetavmel(nx,ny)
#endif

#if defined (SEDIMENT) && defined (SEDTRACK)
      parameter (dstrat = 1000000.)
#  if defined (NYEARSTART) && defined (NYEAREND)
      parameter (nstrat = (NYEAREND-NYEARSTART)/dstrat + 1.001)
      parameter (tstrata = NYEARSTART)
#  else
      parameter (nstrat = 34000000./dstrat + 1.001)
      parameter (tstrata = -34000000.)
#  endif
      dimension sedstrat(nx,ny,nstrat), tstrat(nstrat)
      save tstrat
#endif

      dimension hislist(1000)

      dimension 
     *  work(nx,ny),      work2(nx,ny),     work3(nx,ny,0:nlevp)

      character*20 chist

      logical firsthis
      save firsthis, mcurfirst
      data firsthis /.true./

#include <netcdf.inc>
      integer*4
     * lenattr, nx4, ny4, nlev4, ival4, jval4, nstrat4,
     * varid, rcode, ncerr, ncid, 
     * xid1, yid1, xid0, yid0, levid, timid, minmaxid, tstratid,
     * minmax(2), ndim, idim(4), start(4), count(4), 
     * n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls.
     * ncglobal4, ncunlim4,
     * nf_float4, nf_short4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_short4 = nf_short,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)
      save varid, rcode, ncid,
     *     xid1, yid1, xid0, yid0, levid, timid, minmaxid, tstratid,
     *     minmax, idim, start, count,
     *     n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13
      data minmax /1,2/
      data  n1, n2, n3, n4, n5, n6, n7, n8, n9,n10,
     *     n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,
     *     n21,n22,n23,n24,n25,n26,n27,n28,n29,n30
     *     / 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     *      11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     *      21, 22, 23, 24, 25, 26, 27, 28, 29, 30/
      character cruntitle*80, cvar*16, cvarl*80, cunits*16
      dimension ztmp(20000)

c     parameter (chist = 'history.nc')
      parameter (chist = 'fort.92.nc')

      parameter (vershisnet = 3.0)

c        hislist (namelist, if entered) has precedence over nyearhis

c............................
      if (nhislist.ne.0) then
c............................
        do m=1,nhislist
          dttol = 0.5001*dt
          if (m.eq.1) dttol = 1.0001*dt
          if (abs(timeice - hislist(m)).lt.dttol) then
            mcur = m
            goto 100
          endif
        enddo
        return
  100   continue
c       if (firsthis) mcurfirst = mcur  ! let ifrest=1 append nc file(s)
        mcurfirst = 1                   ! let ifrest=1 append nc file(s)
        nwrit = mcur - mcurfirst + 1
c.................................
      else if (nyearhis.gt.0) then
c.................................
        if ( 
     *     .not.( mod(abs(timeice)+0.5*dt,max(float(nyearhis),dt)).lt.dt
     *           .or. iffirse.eq.1 .or. iflast.eq.1 .or. nyearhis.eq.-1)
     *     ) return 
        nwrit = nint((timeice-nyearstart)/nyearhis) + 1
        if (iffirse.eq.1) nwrit = 1
c..................................
      else if (nyearhis.eq.-1) then
c..................................
        nwrit = nint((timeice-nyearstart)/dt) + 1
        if (iffirse.eq.1) nwrit = 1
c.........
      else
c.........
        return
c..........
      endif
c..........

#if defined (SEDIMENT) && defined (SEDTRACK)
      if (firsthis) then
        do n=1,nstrat
          tstrat(n) = tstrata + dstrat*(n-1)
        enddo
      endif
#endif

c-------------------------
      if (nwrit.gt.1) then
c-------------------------

        rcode = nf_open(chist, nf_write4, ncid)
c       next write needed to avoid ifort/netcdf bug (?)
        write (6,*) 'Opening existing Netcdf history file'

        if (firsthis .and. rcode.ne.0) then
          write (ioterm, '(/3a/)') 
     *      "*** Error in opening previous run's history file ", 
     *      chist(1:lenchr(chist)), " needed for IFREST=1"
        endif

        call checkrcode (rcode, chist, 'open', 5000)

        firsthis = .false.

c---------
      else
c---------

        firsthis = .false.

        rcode = nf_create (chist, nf_clobber4, ncid)
        write (6,*) 'Opening new Netcdf history file'
        call checkrcode (rcode, chist, 'open', 5002) ! stop if open err

c          define title, etc

        cruntitle = 'ice sheet-shelf'
        lenattr = lenchr(cruntitle)
        rcode = nf_put_att_text (ncid, ncglobal4, 'title',
     *                             lenattr, cruntitle)

#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, ncglobal4, 'horiz_grid',
     *                           n7, 'regular')
        rcode = nf_put_att_text (ncid, ncglobal4, 'projection',
     *                           n6, 'latlon')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, ncglobal4, 'projection',
     *                           n11, 'polarstereo')
#endif

c          x dimension, h and v grids

        nx4 = nx
#if defined (PANOPLY)
        rcode = nf_def_dim (ncid, 'lon', nx4, xid1)
        rcode = nf_def_var (ncid, 'lon', nf_float4, n1, xid1, varid)
#else
        rcode = nf_def_dim (ncid, 'x1', nx4, xid1)
        rcode = nf_def_var (ncid, 'x1', nf_float4, n1, xid1, varid)
#endif
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n23, 'x coordinate, h,v-grids')
#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n12, 'degrees_east')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
#endif
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
#if defined (LONLAT)
        do i=1,nx
          ztmp(i) = alond(i,max(ny/2,1))
        enddo
#elif defined (STEREO)
        do i=1,nx
          ztmp(i) = xh(i,max(ny/2,1))*.001
        enddo
#endif
        rcode = nf_put_vara_double (ncid, varid, n1, nx4, ztmp)
        rcode = nf_redef (ncid)

c          y dimension, h and u grids

        ny4 = ny
#if defined (PANOPLY)
        rcode = nf_def_dim (ncid, 'lat', ny4, yid1)
        rcode = nf_def_var (ncid, 'lat', nf_float4, n1, yid1, varid)
#else
        rcode = nf_def_dim (ncid, 'y1', ny4, yid1)
        rcode = nf_def_var (ncid, 'y1', nf_float4, n1, yid1, varid)
#endif
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n23, 'y coordinate, h,u-grids')
#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n13, 'degrees_north')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
#endif
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
#if defined (LONLAT)
        do j=1,ny
          ztmp(j) = alatd(nx/2,j)
        enddo
#elif defined (STEREO)
        do j=1,ny
          ztmp(j) = yh(nx/2,j)*.001
        enddo
#endif
        rcode = nf_put_vara_double (ncid, varid, n1, ny4, ztmp)
        rcode = nf_redef (ncid)

c          x dimension, u-grid (one extra at i=nx, where u=0)

        nx4 = nx
#if defined (PANOPLY)
        rcode = nf_def_dim (ncid, 'lon0', nx4, xid0)
        rcode = nf_def_var (ncid, 'lon0', nf_float4, n1, xid0, varid)
#else
        rcode = nf_def_dim (ncid, 'x0', nx4, xid0)
        rcode = nf_def_var (ncid, 'x0', nf_float4, n1, xid0, varid)
#endif
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n20, 'x coordinate, u-grid')
#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n12, 'degrees_east')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
#endif
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
#if defined (LONLAT)
        j = max (ny/2,1)
        do i=1,nx-1
          ztmp(i) = 0.5 * (alond(i,j)+alond(i+1,j))
        enddo
        ztmp(nx) = 1.5*alond(nx,j) - 0.5*alond(nx-1,j)
#elif defined (STEREO)
        j = max (ny/2,1)
        do i=1,nx-1
          ztmp(i) = 0.5 * (xh(i,j)+xh(i+1,j)) *.001
        enddo
        ztmp(nx) = (xh(nx,j) + 0.5*dx(nx,j)) * .001
#endif
        rcode = nf_put_vara_double (ncid, varid, n1, nx4, ztmp)
        rcode = nf_redef (ncid)

c          y dimension, v-grid (one extra at j=ny, where v=0)

        ny4 = ny
#if defined (PANOPLY)
        rcode = nf_def_dim (ncid, 'lat0', ny4, yid0)
        rcode = nf_def_var (ncid, 'lat0', nf_float4, n1, yid0, varid)
#else
        rcode = nf_def_dim (ncid, 'y0', ny4, yid0)
        rcode = nf_def_var (ncid, 'y0', nf_float4, n1, yid0, varid)
#endif
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n20, 'y coordinate, v-grid')
#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n13, 'degrees_north')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
#endif
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
#if defined (LONLAT)
        i = max (nx/2,1)
        do j=1,ny-1
          ztmp(j) = 0.5 * (alatd(i,j)+alatd(i,j+1))
        enddo
        ztmp(ny) = 1.5*alatd(i,ny) - 0.5*alatd(i,ny-1)
#elif defined (STEREO)
        i = max (nx/2,1)
        do j=1,ny-1
          ztmp(j) = 0.5 * (yh(i,j)+yh(i,j+1)) * .001
        enddo
        ztmp(ny) = (yh(i,ny) + 0.5*dy0) * .001
#endif
        rcode = nf_put_vara_double (ncid, varid, n1, ny4, ztmp)
        rcode = nf_redef (ncid)

c          z dimension

        nlev4 = nlev + 2                                   !  0 to nlevp
        rcode = nf_def_dim (ncid, 'z', nlev4, levid)
        rcode = nf_def_var (ncid, 'z', nf_float4, n1, levid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n23, 'z coordinate within ice')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n18, '(ztop-z)/h, 0 to 1')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
        rcode = nf_put_vara_double (ncid, varid, n1, nlev4, zeta)
        rcode = nf_redef (ncid)

c          age dimension for sedstrat

#if defined (SEDIMENT) && defined (SEDTRACK)
        nstrat4 = nstrat
        rcode = nf_def_dim (ncid, 'tstrat', nstrat4, tstratid)
        rcode = nf_def_var (ncid, 'tstrat', nf_float4, n1, tstratid, 
     *                      varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n27, 'time coordinate for sedstrat')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n5, 'yr BP')
c       rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
c    *                           n5, 'e15.7')
c    *                           n5, 'f12.3')
        rcode = nf_put_att_text (ncid, varid, 'C_format',
     *                           n6, '%15.7e')
        rcode = nf_enddef (ncid)
        rcode = nf_put_vara_double (ncid, varid, n1, nstrat4, tstrat)
        rcode = nf_redef (ncid)
#endif

c          Time dimension 

        rcode = nf_def_dim (ncid, 'time', ncunlim4, timid)
        rcode = nf_def_var (ncid, 'time', nf_float4, n1, timid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n17, 'time (-=BP, +=AP)')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n5, 'years')
c       rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
c    *                           n5, 'e15.7')
c    *                           n5, 'f12.3')
        rcode = nf_put_att_text (ncid, varid, 'C_format',
     *                           n6, '%15.7e')

c         Variable dimensions and attributes

c---

c       For 0D, no time or + time:

        ndim    = 1
        idim(1) = timid

        cvar = 'vershis'
        cvarl = 'history file version number'
        cunits = 'none'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'sealev'
        cvarl = 'sea level'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'toti'
        cvarl = 'ice volume'
        cunits = 'm3'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'totig'
        cvarl = 'grounded ice volume'
        cunits = 'm3'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'totif'
        cvarl = 'floating ice volume'
        cunits = 'm3'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'tota'
        cvarl = 'ice area'
        cunits = 'm2'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'totag'
        cvarl = 'grounded ice area'
        cunits = 'm2'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'totaf'
        cvarl = 'floating ice area'
        cunits = 'm2'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c---

c       For 1D, no time:

#if defined (LONLAT)
        ndim    = 1

        idim(1) = yid1

        cvar = 'lat'
        cvarl = 'Latitude, h-grid'
        cunits = 'degrees_north'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        idim(1) = xid1

        cvar = 'lon'
        cvarl = 'Longitude, h-grid'
        cunits = 'degrees_east'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif
c---

c       For 2D, no time:

        ndim    = 2
        idim(1) = xid1
        idim(2) = yid1

        cvar = 'alatd'
        cvarl = 'latitude, h-grid'
        cunits = 'deg N'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'alond'
        cvarl = 'longitude, h-grid'
        cunits = 'deg E'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'darea'
        cvarl = 'cell area, h-grid'
        cunits = 'm^2'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c---

c       For 2D + time:

        ndim    = 3
        idim(1) = xid1
        idim(2) = yid1
        idim(3) = timid

        cvar = 'h'
        cvarl = 'ice thickness'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'hb'
        cvarl = 'bed elevation'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'hbsd'
        cvarl = 'standard deviation bed elevation'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'topbed'
        cvarl = 'bedrock elevation'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'topbedeq'
        cvarl = 'equilibrium bedrock elevation'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'hs'
        cvarl = 'surface elevation'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'fedge'
        cvarl = 'ice fractional area'
        cunits = '0-1'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'deflect'
        cvarl = 'lithospheric deflection'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'geoflux'
        cvarl = 'geothermal heat flux'
        cunits = 'J/m2/a'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'dgeoid'
        cvarl = 'geoid change since start'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'budgall'
        cvarl = 'net budget'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'budgsnow'
        cvarl = 'snowfall'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'budgrain'
        cvarl = 'rainfall'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'budgmelt'
        cvarl = 'surface melt'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'baseperc'
        cvarl = 'percolation to base'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'basefrml'
        cvarl = 'basal freezing rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'oceanmelt'
        cvarl = 'sub-ice oceanic melt'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'facemelt'
        cvarl = 'facemelt rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'cliffmelt'
        cvarl = 'cliffmelt rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'calvice'
        cvarl = 'calvice rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

#if defined (RACMOCLIM)
        cvar = 'runoff_racmo'
        cvarl = 'RACMO runoff'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

        cvar = 'oceanmeltav'
        cvarl = 'average sub-ice ocean melt'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'facemeltav'
        cvarl = 'average facemelt rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'cliffmeltav'
        cvarl = 'average cliffmelt rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'calviceav'
        cvarl = 'average calvice rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'runliqav'
        cvarl = 'average runliq rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'runfrozav'
        cvarl = 'average runfroz rate'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

#if defined (SEDIMENT)
        cvar = 'sedim'
        cvarl = 'sediment thickness'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'logquarrycoef'
        cvarl = 'quarrying coefficent'
        cunits = '1.e-10 m/yr/(taub*ub)'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'quarryrate'
        cvarl = 'quarrying rate'
        cunits = 'm/Myr'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'quarryacc'
        cvarl = 'accumulated quarrying'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'pelagacc'
        cvarl = 'accumulated pelagic'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

        cvar = 'maskwater'
        cvarl = '0=gd.ice/land, 1=fl.ice/ocean'
        cunits = '0/1'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'logcrhmel'
        cvarl = 'basal sliding coefficient'
        write (cunits,'(a,i1)') 'm/y/Pa^', powb 
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'enhanceinv'
        cvarl = 'extra enhancement (inv) factor for sheet flow'
        cunits = ' '
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits) 
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'tamax'
        cvarl = 'max monthly air temperature'
        cunits = 'C'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'ts'
        cvarl = 'surface temperature'
        cunits = 'C'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'tb'
        cvarl = 'basal temperature'
        cunits = 'C'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'tbhomol'
        cvarl = 'basal homologous temperature'
        cunits = 'C'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'fracm'
        cvarl = 'basal non-froz/lub. fraction'
        cunits = '0-1'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'thickwarm'
        cvarl = 'thickness of basal temperate ice'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'was'
        cvarl = 'surface ice z-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'wab'
        cvarl = 'basal ice z-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'fluxsch'
        cvarl = 'grounding-line flux'
        cunits = 'm2/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

#if defined (CALVDAMAGE) 
        cvar = 'damage_r'
        cvarl = 'crevasse:ice thickness (calvdamage, r)'
        cunits = '0-1'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'damage_term1'
        cvarl = 'n(1-So)Exx (calvdamage, rdot/r)'
        cunits = '1/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'damage_term2'
        cvarl = 'mdot/H (calvdamage, rdot/r)'
        cunits = '1/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'damage_term3'
        cvarl = 'source (calvdamage, rdot)'
        cunits = '1/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

#if defined (CALVBALANCE)
        cvar = 'balance_vice'
        cvarl = 'ice velocity (calvbalance)'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'balance_vbal'
        cvarl = 'balance velocity (calvbalance)'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'balance_vcalv'
        cvarl = 'calving velocity (calvbalance)'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

#if defined (MELANGE)
        cvar = 'hmel'
        cvarl = 'melange thickness'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'maskmel'
        cvarl = '-1 = gd.mel, 0 = no melange access, 1 = fl.mel. access'
        cunits = '-1/0/1'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'budgallmel'
        cvarl = 'net melange budget'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'budgsupmel'
        cvarl = 'ice->melange supply'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

c---

        idim(1) = xid0
        idim(2) = yid1

        cvar = 'utop'
        cvarl = 'surface ice x-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'ubot'
        cvarl = 'basal ice x-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'ua'
        cvarl = 'average ice x-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

#if defined (MELANGE)
        cvar = 'uamel'
        cvarl = 'average melange x-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

c---

        idim(1) = xid1
        idim(2) = yid0

        cvar = 'vtop'
        cvarl = 'surface ice y-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'vbot'
        cvarl = 'basal ice y-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'va'
        cvarl = 'average ice y-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

#if defined (MELANGE)
        cvar = 'vamel'
        cvarl = 'average melange y-velocity'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

c---

        idim(1) = xid1
        idim(2) = yid1

        cvar = 'msector'
        cvarl = 'OH basin sector number'
        cunits = 'integer'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c       For 3D + time:

#if defined (HISTORY3D)
        ndim    = 4
        idim(1) = xid1
        idim(2) = yid1
        idim(3) = levid
        idim(4) = timid

        cvar = 't'
        cvarl = 'ice temperature (3-D)'
        cunits = 'deg K'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'wa'
        cvarl = 'ice vertical velocity (3-D)'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        idim(1) = xid0
        idim(2) = yid1

        cvar = 'u'
        cvarl = 'ice x-velocity (3-D)'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        idim(1) = xid1
        idim(2) = yid0

        cvar = 'v'
        cvarl = 'ice y-velocity (3-D)'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

#if defined (SEDIMENT) && defined (SEDTRACK)
        ndim    = 4
        idim(1) = xid1
        idim(2) = yid1
        idim(3) = tstratid
        idim(4) = timid

        cvar = 'sedstrat'
        cvarl = 'sediment age (Ma)'
        cunits = 'Ma'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

c         Leave define mode

        rcode = nf_enddef (ncid)

c----------
      endif
c----------

c        Add current time to time dimension

      ival4 = nwrit 
      rcode = nf_inq_varid (ncid, 'time', varid)
      rcode = nf_put_var1_double (ncid, varid, ival4, timeice)

c        Write fields

c     For 0D, no time or + time:

      start(1) = nwrit
      count(1) = 1

      if (nwrit.eq.1) then
        rcode = nf_inq_varid (ncid, 'vershis', varid)
        rcode = nf_put_vara_double (ncid, varid, start,count,vershisnet)
      endif

      rcode = nf_inq_varid (ncid, 'sealev', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, sealev)

      toti  = 0.
      totig = 0.
      totif = 0.
      tota  = 0.
      totag = 0.
      totaf = 0.
      do j=1,ny
        do i=1,nx
          if (h(i,j).ne.0.) then
            toti = toti + darea(i,j)*h(i,j)
            tota = tota + darea(i,j)
            if (maskwater(i,j).eq.0) then
              totig = totig + darea(i,j)*h(i,j)
              totag = totag + darea(i,j)
            else
              totif = totif + darea(i,j)*h(i,j)
              totaf = totaf + darea(i,j)
            endif
          endif
        enddo
      enddo

      rcode = nf_inq_varid (ncid, 'toti', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, toti)

      rcode = nf_inq_varid (ncid, 'totig', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, totig)

      rcode = nf_inq_varid (ncid, 'totif', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, totif)

      rcode = nf_inq_varid (ncid, 'tota', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, tota)

      rcode = nf_inq_varid (ncid, 'totag', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, totag)

      rcode = nf_inq_varid (ncid, 'totaf', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, totaf)

c     For 1D, no time or + time:

#if defined (LONLAT)
      if (nwrit.eq.1) then
        start(1) = 1
        count(1) = ny
        start(2) = nwrit
        count(2) = 1

        do j=1,ny
          ztmp(j) = alatd(nx/2,j)
        enddo
        rcode = nf_inq_varid (ncid, 'lat', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, ztmp)

        start(1) = 1
        count(1) = nx
        start(2) = nwrit
        count(2) = 1

        do i=1,nx
          ztmp(i) = alond(i,max(ny/2,1))
        enddo
        rcode = nf_inq_varid (ncid, 'lon', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, ztmp)
      endif
#endif

c     For 2D, no time or + time:

      start(1) = 1
      count(1) = nx
      start(2) = 1
      count(2) = ny
      start(3) = nwrit
      count(3) = 1

      if (nwrit.eq.1) then
        rcode = nf_inq_varid (ncid, 'alatd', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, alatd)

        rcode = nf_inq_varid (ncid, 'alond', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, alond)

        rcode = nf_inq_varid (ncid, 'darea', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, darea)
      endif

      rcode = nf_inq_varid (ncid, 'h', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, h)

      rcode = nf_inq_varid (ncid, 'hb', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hb)

      rcode = nf_inq_varid (ncid, 'hbsd', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hbsd)

      rcode = nf_inq_varid (ncid, 'topbed', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, topbed)

      rcode = nf_inq_varid (ncid, 'topbedeq', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, topbedeq)

      rcode = nf_inq_varid (ncid, 'hs', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hs)

      rcode = nf_inq_varid (ncid, 'fedge', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, fedge)

      rcode = nf_inq_varid (ncid, 'deflect', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, deflect)

      rcode = nf_inq_varid (ncid, 'geoflux', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, geoflux)

      rcode = nf_inq_varid (ncid, 'dgeoid', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, dgeoid)

      rcode = nf_inq_varid (ncid, 'budgall', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, budgall)

      rcode = nf_inq_varid (ncid, 'budgsnow', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, budgsnow)

      rcode = nf_inq_varid (ncid, 'budgrain', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, budgrain)

      rcode = nf_inq_varid (ncid, 'budgmelt', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, budgmelt)

      rcode = nf_inq_varid (ncid, 'baseperc', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, baseperc)

      rcode = nf_inq_varid (ncid, 'basefrml', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, basefrml)

      rcode = nf_inq_varid (ncid, 'oceanmelt', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, oceanmelt)

      rcode = nf_inq_varid (ncid, 'facemelt', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, facemelt)

      rcode = nf_inq_varid (ncid, 'cliffmelt', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, cliffmelt)

      rcode = nf_inq_varid (ncid, 'calvice', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, calvice)

#if defined (RACMOCLIM)
      rcode = nf_inq_varid (ncid, 'runoff_racmo', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, runofa)
#endif

      do j=1,ny
        do i=1,nx
          work(i,j) = oceanmeltav(i,j) / max(.001, timeice-timefluxprev)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'oceanmeltav', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = facemeltav(i,j) / max (.001, timeice-timefluxprev)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'facemeltav', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) =cliffmeltav(i,j) / max (.001, timeice-timefluxprev)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'cliffmeltav', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = calviceav(i,j) / max (.001, timeice-timefluxprev)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'calviceav', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = runliqav(i,j) / max (.001, timeice-timefluxprev)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'runliqav', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = runfrozav(i,j) / max (.001, timeice-timefluxprev)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'runfrozav', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

c        Reset ocean-related liquid and frozen flux accumulators
c        (accumulated in subr icetabbud). Except...for Loveclip,
c        runliqav and runfrozav (and other *av) are accumulated over 
c        whole run (chunk), and used in writeloveclim.    
   
#if ! defined (LOVECLIP)
      call zero (oceanmeltav,nx*ny)
      call zero (facemeltav,nx*ny)
      call zero (cliffmeltav,nx*ny)
      call zero (calviceav,nx*ny)
      call zero (runliqav,nx*ny)
      call zero (runfrozav,nx*ny)
      timefluxprev = timeice
#endif

#if defined (SEDIMENT)
      rcode = nf_inq_varid (ncid, 'sedim', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, sedim)

      do j=1,ny
        do i=1,nx
          work(i,j) = log10(quarrycoef(i,j))
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'logquarrycoef', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = quarryrate(i,j)*1.e6
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'quarryrate', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      rcode = nf_inq_varid (ncid, 'quarryacc', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, quarryacc)

      rcode = nf_inq_varid (ncid, 'pelagacc', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, pelagacc)
#endif

      call itor (maskwater, work, nx*ny)
      rcode = nf_inq_varid (ncid, 'maskwater', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = log10(crhmel(i,j))
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'logcrhmel', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      rcode = nf_inq_varid (ncid, 'enhanceinv', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, enhanceinv)

      do j=1,ny
        do i=1,nx
          work(i,j) = -1000.
          do m=1,nmon
            work(i,j) = max (work(i,j), tmon(m,i,j)-tmelt)
          enddo
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'tamax', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = t(i,j,0) - tmelt
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'ts', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = t(i,j,nlevp) - tmelt
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'tb', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = t(i,j,nlevp) + dtmdh*h(i,j) - tmelt
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'tbhomol', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          call calc_crh (zcrh, work(i,j),
     *                   crhmel(i,j), 
     *                   t(i,j,nlevp) - (tmelt-dtmdh*h(i,j)),
     *                   baseperc(i,j)-basefrml(i,j),
     *                   h(i,j), 
     *                   hbsd(i,j), 
     *                   1)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'fracm', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      call calcz0 (work, h, t)
      rcode = nf_inq_varid (ncid, 'thickwarm', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = wa(i,j,0)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'was', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = wa(i,j,nlevp)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'wab', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

c     grounding-line flux out of last-grounded h-grid boxes, (m2/y),
c     averaged only around sides adjacent to non-grounded h-grid boxes:
      do j=1,ny
        jm1 = max (1, j-1)
        jp1 = min (ny,j+1)
        do i=1,nx
          zflx = 0.
          zlen = 0.
          if (maskwater(i,j).eq.0) then 
            im1 = max (1, i-1)
            ip1 = min (nx,i+1)
            if (maskwater(ip1,j).ne.0) then
              zflx = zflx + fluxschu(i,j)*dyu(i,j)
              zlen = zlen + dyu(i,j)
            endif
            if (maskwater(ip1,j).ne.0) then
              zflx = zflx + fluxschv(i,j)*dxv(i,j)
              zlen = zlen + dxv(i,j)
            endif
            if (maskwater(im1,j).ne.0) then
              zflx = zflx - fluxschu(i-1,j)*dyu(i-1,j)
              zlen = zlen + dyu(i-1,j)
            endif
            if (maskwater(i,jm1).ne.0) then
              zflx = zflx - fluxschv(i,j-1)*dxv(i,j-1)
              zlen = zlen + dxv(i,j-1)
            endif
          endif
          if (zlen.gt.0) then
            work(i,j) = zflx/zlen
          else
            work(i,j) = 0.
          endif
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'fluxsch', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

#if defined (CALVDAMAGE)
      rcode = nf_inq_varid (ncid, 'damage_r', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, 
     *                            tracer(1,1,1,2))

      rcode = nf_inq_varid (ncid, 'damage_term1', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, calvdam_t1)

      rcode = nf_inq_varid (ncid, 'damage_term2', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, calvdam_t2)

      rcode = nf_inq_varid (ncid, 'damage_term3', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, calvdam_t3)
#endif

#if defined (CALVBALANCE)
      rcode = nf_inq_varid (ncid, 'balance_vice', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, calvbal_vi)

      rcode = nf_inq_varid (ncid, 'balance_vbal', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, calvbal_vb)

      rcode = nf_inq_varid (ncid, 'balance_vcalv', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, calvbal_vc)
#endif

#if defined (MELANGE)
      rcode = nf_inq_varid (ncid, 'hmel', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hmel)

      call itor (maskmel, work, nx*ny)
      rcode = nf_inq_varid (ncid, 'maskmel', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      rcode = nf_inq_varid (ncid, 'budgallmel', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, budgallmel)

      rcode = nf_inq_varid (ncid, 'budgsupmel', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, budgsupmel)

      do j=1,ny
        do i=1,nx
          work(i,j) = 0.5*(uamel(i,j)+uamel(i-1,j))
          if (hmel(i,j).eq.0.) work(i,j) = 0.
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'uamel', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = 0.5*(vamel(i,j)+vamel(i,j-1))
          if (hmel(i,j).eq.0.) work(i,j) = 0.
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'vamel', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)
#endif

      do j=1,ny
        do i=1,nx
          work(i,j) = 0.5*(u(i,j,0)+u(i-1,j,0))
          if (h(i,j).eq.0.) work(i,j) = 0.
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'utop', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = 0.5*(ub(i,j)+ub(i-1,j))
          if (h(i,j).eq.0.) work(i,j) = 0.
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'ubot', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = 0.5*(ua(i,j)+ua(i-1,j))
          if (h(i,j).eq.0.) work(i,j) = 0.
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'ua', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = 0.5*(v(i,j,0)+v(i,j-1,0))
          if (h(i,j).eq.0.) work(i,j) = 0.
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'vtop', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = 0.5*(vb(i,j)+vb(i,j-1))
          if (h(i,j).eq.0.) work(i,j) = 0.
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'vbot', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = 0.5*(va(i,j)+va(i,j-1))
          if (h(i,j).eq.0.) work(i,j) = 0.
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'va', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

      do j=1,ny
        do i=1,nx
          work(i,j) = msector(i,j)
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'msector', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

c     For 3D + time, additional:

#if defined (HISTORY3D)
      start(3) = 1
      count(3) = nlev + 2
      start(4) = nwrit
      count(4) = 1

      rcode = nf_inq_varid (ncid, 't', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, t)

      rcode = nf_inq_varid (ncid, 'wa', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, wa)

      do k=0,nlevp
        do j=1,ny
          do i=1,nx
            work3(i,j,k) = 0.5*(u(i,j,k)+u(i-1,j,k))
            if (h(i,j).eq.0.) work3(i,j,k) = 0.
          enddo
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'u', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work3)

      do k=0,nlevp
        do j=1,ny
          do i=1,nx
            work3(i,j,k) = 0.5*(v(i,j,k)+v(i,j-1,k))
            if (h(i,j).eq.0.) work3(i,j,k) = 0.
          enddo
        enddo
      enddo
      rcode = nf_inq_varid (ncid, 'v', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work3)
#endif

#if defined (SEDIMENT) && defined (SEDTRACK)
      start(3) = 1
      count(3) = nstrat
      start(4) = nwrit
      count(4) = 1

c        Check sedtrack is monotonic in time 

      do j=1,ny
        do i=1,nx
          if (sedim(i,j).gt.0. .and. itrtop(i,j).ge.1) then
            do k=0,itrtop(i,j)-1
              kp = min (k+1,ntrack)
              if (sedtrack(i,j,k) .eq.vnulsedtrack .or.
     *            sedtrack(i,j,kp).eq.vnulsedtrack .or.
     *            sedtrack(i,j,k).gt.sedtrack(i,j,kp)) then
                write (6,'(/a/2i4,3f12.2)') 
     *             'Error: iceshow, sedtrack not monotonic:', 
     *             i, j, sedim(i,j), hb(i,j), topbed(i,j)
                do l=0,itrtop(i,j)
                  write (6,'(i4,2f15.2)') l, l*dtrack, sedtrack(i,j,l)
                enddo
                stop
              endif
            enddo
          endif
        enddo
      enddo

c         Interpolate sedtrack (times vs z-topbed) to sedstrat 
c         (z vs. times)
   
      do j=1,ny
        do i=1,nx

c...................................
          if (sedim(i,j).gt.0.) then
c...................................

c               For tstrat within the time-range of existing sediment 

            if (itrtop(i,j).ge.1) then
              do n=1,nstrat
                zt = tstrat(n)
                do k=0,itrtop(i,j)-1
                  if (sedtrack(i,j,k)  .le.zt .and.
     *                sedtrack(i,j,k+1).ge.zt) then
                    zwei = (zt - sedtrack(i,j,k)) / 
     *                     (sedtrack(i,j,k+1) - sedtrack(i,j,k))
                    sedstrat(i,j,n) = topbed(i,j) + dtrack*(k+zwei) 
                    goto 70
                  endif
                enddo
   70           continue
              enddo
            endif

c               For 1st tstrat younger than all sedim, set to top sedim

            do n=nstrat,2,-1
              if ( tstrat(n-1).lt.sedtrack(i,j,itrtop(i,j)) .and.
     *             tstrat(n)  .gt.sedtrack(i,j,itrtop(i,j)) ) then
                sedstrat(i,j,n) = hb(i,j)
c               for any younger tstrat, also set to sedim top
c               (for plotting in ginsedtrack):
                if (n.lt.nstrat) then
                  do nn=n+1,nstrat
                    sedstrat(i,j,nn) = hb(i,j)
                  enddo
                endif
                goto 80
              endif
            enddo
   80       continue

c               For 1st tstrat older than all sedim, set to bottom sedim

            do n=1,nstrat-1
              if ( tstrat(n)  .lt.sedtrack(i,j,0) .and.
     *             tstrat(n+1).gt.sedtrack(i,j,0) ) then
                sedstrat(i,j,n) = topbed(i,j)
c               for any older tstrat, also set to sedim bottom
c               (for plotting in ginsedtrack):
                if (n.gt.1) then
                  do nn=1,n-1
                    sedstrat(i,j,nn) = topbed(i,j)
                  enddo
                endif
                goto 90
              endif
            enddo
   90       continue

c.............
          else
c.............

c             If no sedim, set all sedstrat = bed surface 
c             (for plotting in ginsedtrack)
 
            do n=1,nstrat
              sedstrat(i,j,n) = hb(i,j)
            enddo

c..............
          endif
c..............

        enddo
      enddo

      rcode = nf_inq_varid (ncid, 'sedstrat', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, sedstrat)
#endif

c        Return to define mode

      rcode = nf_redef (ncid)

      rcode = nf_close (ncid)                ! to flush, re-opened above

      return
      end

c-----------------------------------------------------------------------

c@@@@@@@@@@@@@@@@@@@@@@@@@@
#if defined (LAGRANGEWRITE)
c@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine writelag (h, hb, hs, hw, maskwater, u, v, w, timeice,
     *                     dt, nyearstart, iffirse, iflast, nyearlag)

c     Write fields to Netcdf file for post-processing Lagrangian 
c     tracing. Grids, etc, are defined as in writehis.
c
c     Unlike writehis, only for -DSTEREO (not LONLAT). 
c     Unlike writehis, u and v are written with full array dimensions.  
c     Unlike writehis, no write for first timestep.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),             hb(nx,ny),  
     *  hs(nx,ny),            hw(nx,ny),
     *  maskwater(nx,ny),     w(nx,ny,0:nlevp)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp)

      dimension 
     *  work(nx,ny)

      character*20 chist

      logical firsthis
      save firsthis
      data firsthis /.true./

#include <netcdf.inc>
      integer*4
     * lenattr, nx4, ny4, nlev4, ival4, jval4,
     * varid, rcode, ncerr, ncid, 
     * xid1, yid1, xid0, yid0, levid, timid, minmaxid,
     * minmax(2), ndim, idim(4), start(4), count(4), 
     * n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,
c      following are from netcdf.inc...need to convert to int*4,
c      and use *4 variables as arguments in netcdf calls.
     * ncglobal4, ncunlim4,
     * nf_float4, nf_short4, nf_clobber4, nf_write4, nf_noerr4
      parameter (ncglobal4   = ncglobal,   ncunlim4  = ncunlim,
     *           nf_float4   = nf_float,   nf_short4 = nf_short,
     *           nf_clobber4 = nf_clobber, nf_write4 = nf_write,
     *           nf_noerr4   = nf_noerr)
      save varid, rcode, ncid,
     *     xid1, yid1, xid0, yid0, levid, timid, minmaxid,
     *     minmax, idim, start, count,
     *     n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13
      data minmax /1,2/
      data  n1, n2, n3, n4, n5, n6, n7, n8, n9,n10,
     *     n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,
     *     n21,n22,n23,n24,n25,n26,n27,n28,n29,n30
     *     / 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     *      11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     *      21, 22, 23, 24, 25, 26, 27, 28, 29, 30/
      character cruntitle*80, cvar*16, cvarl*80, cunits*16
      dimension ztmp(0:20000)

      parameter (chist = 'fort.93.nc')

      parameter (vershisnet = 3.0)

      if (nyearlag.gt.0) then
        if ( 
     *     .not.( mod(abs(timeice)+0.5*dt,max(float(nyearlag),dt)).lt.dt
     *                             .or. iflast.eq.1 .or. nyearlag.eq.-1)
     *     ) return 
        nwrit = nint((timeice-nyearstart)/nyearlag)
      else
        return
      endif

      firsthis = .false.

c-------------------------
      if (nwrit.gt.1) then
c-------------------------

        rcode = nf_open(chist, nf_write4, ncid)
c       next write needed to avoid ifort/netcdf bug (?)
        write (6,*) 'Opening existing Netcdf lagrangian file'

c---------
      else
c---------

        rcode = nf_create (chist, nf_clobber4, ncid)
        write (6,*) 'Opening new Netcdf lagrangian file'

c          define title, etc

        cruntitle = 'ice lagrangian'
        lenattr = lenchr(cruntitle)
        rcode = nf_put_att_text (ncid, ncglobal4, 'title',
     *                             lenattr, cruntitle)

#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, ncglobal4, 'horiz_grid',
     *                           n7, 'regular')
        rcode = nf_put_att_text (ncid, ncglobal4, 'projection',
     *                           n6, 'latlon')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, ncglobal4, 'projection',
     *                           n11, 'polarstereo')
#endif

c          x dimension, h grid

        nx4 = nx
#if defined (PANOPLY)
        rcode = nf_def_dim (ncid, 'lon', nx4, xid1)
        rcode = nf_def_var (ncid, 'lon', nf_float4, n1, xid1, varid)
#else
        rcode = nf_def_dim (ncid, 'x1', nx4, xid1)
        rcode = nf_def_var (ncid, 'x1', nf_float4, n1, xid1, varid)
#endif
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n23, 'x coordinate, h-grid')
#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n12, 'degrees_east')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
#endif
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
#if defined (LONLAT)
        do i=1,nx
          ztmp(i) = alond(i,max(ny/2,1))
        enddo
#elif defined (STEREO)
        do i=1,nx
          ztmp(i) = xh(i,max(ny/2,1))*.001
        enddo
#endif
        rcode = nf_put_vara_double (ncid, varid, n1, nx4, ztmp(1))
        rcode = nf_redef (ncid)

c          y dimension, h grid

        ny4 = ny
#if defined (PANOPLY)
        rcode = nf_def_dim (ncid, 'lat', ny4, yid1)
        rcode = nf_def_var (ncid, 'lat', nf_float4, n1, yid1, varid)
#else
        rcode = nf_def_dim (ncid, 'y1', ny4, yid1)
        rcode = nf_def_var (ncid, 'y1', nf_float4, n1, yid1, varid)
#endif
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n23, 'y coordinate, h-grid')
#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n13, 'degrees_north')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
#endif
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
#if defined (LONLAT)
        do j=1,ny
          ztmp(j) = alatd(nx/2,j)
        enddo
#elif defined (STEREO)
        do j=1,ny
          ztmp(j) = yh(nx/2,j)*.001
        enddo
#endif
        rcode = nf_put_vara_double (ncid, varid, n1, ny4, ztmp(1))
        rcode = nf_redef (ncid)

c          x dimension, u and v (0 to nx+1)

        nx4 = nx + 2
#if defined (PANOPLY)
        rcode = nf_def_dim (ncid, 'lon0', nx4, xid0)
        rcode = nf_def_var (ncid, 'lon0', nf_float4, n1, xid0, varid)
#else
        rcode = nf_def_dim (ncid, 'x0', nx4, xid0)
        rcode = nf_def_var (ncid, 'x0', nf_float4, n1, xid0, varid)
#endif
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n20, 'x coordinate, u,v-grid')
#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n12, 'degrees_east')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
#endif
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
#if defined (LONLAT)
        j = max (ny/2,1)
        ztmp(0) = alond(1,j) - 0.5*(alond(2,j)-alond(1,j)) 
        do i=1,nx-1
          ztmp(i) = 0.5 * (alond(i,j)+alond(i+1,j))
        enddo
        ztmp(nx)   = alond(nx,j) + 0.5(alond(nx,j)-alond(nx-1,j))
        ztmp(nx+1) = alond(nx,j) + 1.5(alond(nx,j)-alond(nx-1,j))
#elif defined (STEREO)
        j = max (ny/2,1)
        ztmp(0) = (xh(1,j) - 0.5*dx0) *.001
        do i=1,nx-1
          ztmp(i) = 0.5 * (xh(i,j)+xh(i+1,j)) *.001
        enddo
        ztmp(nx)   = (xh(nx,j) + 0.5*dx0) * .001
        ztmp(nx+1) = (xh(nx,j) + 1.5*dx0) * .001
#endif
        rcode = nf_put_vara_double (ncid, varid, n1, nx4, ztmp(0))
        rcode = nf_redef (ncid)

c          y dimension, u and v (0 to ny+1)

        ny4 = ny + 2
#if defined (PANOPLY)
        rcode = nf_def_dim (ncid, 'lat0', ny4, yid0)
        rcode = nf_def_var (ncid, 'lat0', nf_float4, n1, yid0, varid)
#else
        rcode = nf_def_dim (ncid, 'y0', ny4, yid0)
        rcode = nf_def_var (ncid, 'y0', nf_float4, n1, yid0, varid)
#endif
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n20, 'y coordinate, u,v-grid')
#if defined (LONLAT)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n13, 'degrees_north')
#elif defined (STEREO)
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n2, 'km')
#endif
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
#if defined (LONLAT)
        i = max (nx/2,1)
        ztmp(0) = alatd(i,1) - 0.5*(alatd(i,2)-alatd(i,1)) 
        do j=1,ny-1
          ztmp(j) = 0.5 * (alatd(i,j)+alatd(i,j+1))
        enddo
        ztmp(ny) = 1.5*alatd(i,ny) - 0.5*alatd(i,ny-1)
        ztmp(ny)   = alatd(i,ny) + 0.5(alatd(i,ny)-alatd(i,ny-1))
        ztmp(ny+1) = alatd(i,ny) + 1.5(alatd(i,ny)-alatd(i,ny-1))
#elif defined (STEREO)
        i = max (nx/2,1)
        ztmp(0) = (yh(i,1) - 0.5*dy0) *.001
        do j=1,ny-1
          ztmp(j) = 0.5 * (yh(i,j)+yh(i,j+1)) * .001
        enddo
        ztmp(ny)   = (yh(i,ny) + 0.5*dy0) * .001
        ztmp(ny+1) = (yh(i,ny) + 1.5*dy0) * .001
#endif
        rcode = nf_put_vara_double (ncid, varid, n1, ny4, ztmp(0))
        rcode = nf_redef (ncid)

c          z dimension

        nlev4 = nlev + 2                                   !  0 to nlevp
        rcode = nf_def_dim (ncid, 'z', nlev4, levid)
        rcode = nf_def_var (ncid, 'z', nf_float4, n1, levid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n23, 'z coordinate within ice')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n18, '(ztop-z)/h, 0 to 1')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n4, 'f8.3')
        rcode = nf_enddef (ncid)
        rcode = nf_put_vara_double (ncid, varid, n1, nlev4, zeta)
        rcode = nf_redef (ncid)

c          Time dimension 

        rcode = nf_def_dim (ncid, 'time', ncunlim4, timid)
        rcode = nf_def_var (ncid, 'time', nf_float4, n1, timid, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name',
     *                           n17, 'time (-=BP, +=AP)')
        rcode = nf_put_att_text (ncid, varid, 'units',
     *                           n5, 'years')
        rcode = nf_put_att_text (ncid, varid, 'FORTRAN_format',
     *                           n5, 'f12.3')

c         Variable dimensions and attributes

c---

c       For 0D, no time or + time:

        ndim    = 1
        idim(1) = timid

        cvar = 'vershis'
        cvarl = 'history file version number'
        cunits = 'none'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c---

c       For 1D, no time:

#if defined (LONLAT)
        ndim    = 1

        idim(1) = yid1

        cvar = 'lat'
        cvarl = 'Latitude, h-grid'
        cunits = 'degrees_north'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        idim(1) = xid1

        cvar = 'lon'
        cvarl = 'Longitude, h-grid'
        cunits = 'degrees_east'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')
#endif

c---

c       For 2D, no time:

        ndim    = 2
        idim(1) = xid1
        idim(2) = yid1

        cvar = 'alatd'
        cvarl = 'latitude, h-grid'
        cunits = 'deg N'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'alond'
        cvarl = 'longitude, h-grid'
        cunits = 'deg E'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'darea'
        cvarl = 'cell area, h-grid'
        cunits = 'm^2'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'xh'
        cvarl = 'x coord (2-D), h-grid'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'yh'
        cvarl = 'y coord (2-D), h-grid'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c---

c       For 2D + time:

        ndim    = 3
        idim(1) = xid1
        idim(2) = yid1
        idim(3) = timid

        cvar = 'h'
        cvarl = 'ice thickness'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'hb'
        cvarl = 'bed elevation'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'hs'
        cvarl = 'surface elevation'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'hw'
        cvarl = 'water thickness'
        cunits = 'm'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        cvar = 'maskwater'
        cvarl = '0=gd.ice/land, 1=fl.ice/ocean'
        cunits = '0/1'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c---

c       For 3D + time:

        ndim    = 4
        idim(1) = xid1
        idim(2) = yid1
        idim(3) = levid
        idim(4) = timid

        cvar = 'w'
        cvarl = 'ice vertical velocity d(zeta)/dt (3-D)'
        cunits = '1/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        idim(1) = xid0
        idim(2) = yid0

        cvar = 'u'
        cvarl = 'ice x-velocity (3-D)'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

        idim(1) = xid0
        idim(2) = yid0

        cvar = 'v'
        cvarl = 'ice y-velocity (3-D)'
        cunits = 'm/y'
        ival4 = lenchr(cvarl)
        jval4 = lenchr(cunits)
        rcode = nf_def_var (ncid, cvar, nf_float4, ndim, idim, varid)
        rcode = nf_put_att_text (ncid, varid, 'long_name', ival4, cvarl)
        rcode = nf_put_att_text (ncid, varid, 'units', jval4, cunits)
        rcode = nf_put_att_text (ncid,varid,'FORTRAN_format',n5,'e13.5')

c         Leave define mode

        rcode = nf_enddef (ncid)
c----------
      endif
c----------


c        Add current time to time dimension

      ival4 = nwrit 
      rcode = nf_inq_varid (ncid, 'time', varid)
      rcode = nf_put_var1_double (ncid, varid, ival4, timeice)

c        Write fields

c     For 0D, no time or + time:

      start(1) = nwrit
      count(1) = 1

      if (nwrit.eq.1) then
        rcode = nf_inq_varid (ncid, 'vershis', varid)
        rcode = nf_put_vara_double (ncid, varid, start,count,vershisnet)
      endif

c     For 1D, no time or + time:

#if defined (LONLAT)
      if (nwrit.eq.1) then
        start(1) = 1
        count(1) = ny
        start(2) = nwrit
        count(2) = 1

        do j=1,ny
          ztmp(j) = alatd(nx/2,j)
        enddo
        rcode = nf_inq_varid (ncid, 'lat', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, ztmp)

        start(1) = 1
        count(1) = nx
        start(2) = nwrit
        count(2) = 1

        do i=1,nx
          ztmp(i) = alond(i,max(ny/2,1))
        enddo
        rcode = nf_inq_varid (ncid, 'lon', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, ztmp)
      endif
#endif

c     For 2D, no time or + time:

      start(1) = 1
      count(1) = nx
      start(2) = 1
      count(2) = ny
      start(3) = nwrit
      count(3) = 1

      if (nwrit.eq.1) then
        rcode = nf_inq_varid (ncid, 'alatd', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, alatd)

        rcode = nf_inq_varid (ncid, 'alond', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, alond)

        rcode = nf_inq_varid (ncid, 'darea', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, darea)

        rcode = nf_inq_varid (ncid, 'xh', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, xh)

        rcode = nf_inq_varid (ncid, 'yh', varid)
        rcode = nf_put_vara_double (ncid, varid, start, count, yh)
      endif

      rcode = nf_inq_varid (ncid, 'h', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, h)

      rcode = nf_inq_varid (ncid, 'hb', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hb)

      rcode = nf_inq_varid (ncid, 'hs', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hs)

      rcode = nf_inq_varid (ncid, 'hw', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, hw)

      call itor (maskwater, work, nx*ny)
      rcode = nf_inq_varid (ncid, 'maskwater', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, work)

c     For 3D + time, additional:

      start(1) = 1
      count(1) = nx
      start(2) = 1
      count(2) = ny
      start(3) = 1
      count(3) = nlev + 2
      start(4) = nwrit
      count(4) = 1

      rcode = nf_inq_varid (ncid, 'w', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, w)

      start(1) = 1
      count(1) = nx + 2
      start(2) = 1
      count(2) = ny + 2

      rcode = nf_inq_varid (ncid, 'u', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, u)

      rcode = nf_inq_varid (ncid, 'v', varid)
      rcode = nf_put_vara_double (ncid, varid, start, count, v)

c        Return to define mode

      rcode = nf_redef (ncid)

      rcode = nf_close (ncid)                ! to flush, re-opened above

      return
      end

c@@@@@
#endif
c@@@@@

c-----------------------------------------------------------------------

      subroutine mascalc (arr, totv, tota)

c     Calculates total (domain) damount of variable arr (totv),
c     and total area with non-zero arr (tota)  

#include <comicephys.h>
#include <comicegrid.h>

      dimension arr(nx,ny)

      totv  = 0.
      tota  = 0.
      do j=1,ny
        do i=1,nx
          if (arr(i,j).ne.0.) then 
            totv = totv + darea(i,j)*arr(i,j)
            tota = tota + darea(i,j)
          endif
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine calcz0 (thickz0, h, t)

c     Calculate thickness of lower ice from first-from-top layer with 
c     homologous temperature >= tcritz0 to the base (thickz0), 
c     for output (history field "thickwarm"), i.e., thickness of basal 
c     temperate layer. Assumes dT/dz always warming down into ice.
c     If no ice, thickz0 = 0. 
c     If all ice temps < tcritz0, thickz0 = 0 (n.b. not h).
c     If all (i.e. top-layer) ice temp >= tcritz0, thickz0 = h.

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  thickz0(nx,ny), h(nx,ny), t(nx,ny,0:nlevp)

c     parameter (tcritz0 =  0.)    ! deg C (temperate ice)
      parameter (tcritz0 = -0.1)   ! deg C (temperate ice)

      do j=1,ny
        do i=1,nx
c--------------------------------
          if (h(i,j).gt.0.1) then
c--------------------------------
c           look down from sfc for first ice layer temp >= tcritz0
            do k=1,nlev
              ztme = tmelt 
              if (t(i,j,k) - (tmelt-dtmdh*h(i,j)*zeta(k)) .ge. tcritz0)
     *          then 
c               distance from top of layer to base:
                thickz0(i,j) = (1.-zetah(k-1))*h(i,j) 
                goto 10
              endif
            enddo
            thickz0(i,j) = 0.                  ! all ice temps < tcritz0
   10       continue
c-------------
          else
c-------------
            thickz0(i,j) = 0.                                   ! no ice
c--------------
          endif
c--------------
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

c}}}}}}}}}}}}}}}}}}}}}}
#if defined (TEST2DOUT)
c}}}}}}}}}}}}}}}}}}}}}}

      subroutine icemismipplus (h, hs, hb, 
     *                          w, wa, maskh,
     *                          hw, maskwater, 
     *                          u, v, ua, va, uia, via, ub, vb,
     *                          uadv, vadv, hu, hv,
     *                          masku, maskv, muind, mvind, crhu, crhv,
     *                          fracgu, fracgv, thetau, thetav,
     *                          hgu, hgv, qgu, qgv,
     *                          uschoof, vschoof,
     *                          sealev, timeice, dt, 
     *                          iffirse, iffirsr, iflast)
 
c     Writes ascii tabular files for mismipplus, snapshots in time. 
c     Written every nyearmismipplus years (set below)

#include <comicephys.h>
#include <comicegrid.h>

      dimension 
     *  h(nx,ny),            hs(nx,ny),          hb(nx,ny), 
     *  w(nx,ny,0:nlevp),    wa(nx,ny,0:nlevp),  maskh(nx,ny),
     *  hw(nx,ny),           maskwater(nx,ny)

      dimension
     *  u(0:nxp,0:nyp,0:nlevp), v(0:nxp,0:nyp,0:nlevp),
     *  ua(0:nxp,0:nyp),        va(0:nxp,0:nyp),
     *  uia(0:nxp,0:nyp),       via(0:nxp,0:nyp),
     *  ub(0:nxp,0:nyp),        vb(0:nxp,0:nyp),
     *  uadv(0:nxp,0:nyp),      vadv(0:nxp,0:nyp),
     *  hu(0:nxp,0:nyp),        hv(0:nxp,0:nyp),
     *  masku(0:nxp,0:nyp),     maskv(0:nxp,0:nyp), 
     *  muind(0:nxp,0:nyp),     mvind(0:nxp,0:nyp), 
     *  crhu(0:nxp,0:nyp),      crhv(0:nxp,0:nyp),
     *  fracgu(0:nxp,0:nyp),    fracgv(0:nxp,0:nyp), 
     *  thetau(nx,ny),          thetav(nx,ny),
     *  hgu(0:nxp,0:nyp),       hgv(0:nxp,0:nyp),
     *  qgu(0:nxp,0:nyp),       qgv(0:nxp,0:nyp),
     *  uschoof(0:nxp,0:nyp),   vschoof(0:nxp,0:nyp)

      character cfile*80, runlab*8

      parameter (nglmax=4*(nx+ny))
      dimension
     *  xgl(nglmax),   ygl(nglmax),
     *  hgl(nglmax),
     *  ubgl(nglmax),  uagl(nglmax), usgl(nglmax),
     *  vbgl(nglmax),  vagl(nglmax), vsgl(nglmax),
     *  ifugl(nglmax),
     *  indgl(nglmax), ifsel(nglmax)
      parameter (vnull=9999.)

      dimension work(nx,ny), work2(nx,ny), work3(nx,ny)

      logical firstmismipplus
      save firstmismipplus
      data firstmismipplus /.true./

#if defined (MISMIPPLUS_spin)
      runlab = '_spin_'
      nyearmismipplus = 100
#else
#  if defined (MISMIPPLUS_ctl)
      runlab = '_ctl_'
#  elif defined (MISMIPPLUS_1rr)
      runlab = '_1rr_'
#  elif defined (MISMIPPLUS_1ra)
      runlab = '_1ra_'
#  elif defined (MISMIPPLUS_2rr)
      runlab = '_2rr_'         
#  elif defined (MISMIPPLUS_2ra)
      runlab = '_2ra_'         
#  else
      write (6,*) 'Error (icemismipplus): must specify MISMIPPLUS_...'
      stop
#  endif
      nyearmismipplus = 10
      if (nint(timeice).gt.200) nyearmismipplus = 100
#endif

      if ( .not. 
     *   ( mod(abs(timeice)+0.5*dt,max(float(nyearmismipplus),dt)).lt.dt
     *    .or. iffirse.eq.1 .or. iffirsr.eq.1 .or. iflast.eq.1)
     *   ) return 

      iua = iumismipplus_a
      iub = iumismipplus_b

      if (firstmismipplus) then
        write (cfile,'(3a)') 
     *    'mismip', runlab(1:lenchr(runlab)), 'gl.dat'
        open (iua, file=cfile, status='unknown')
        close (iua, status='delete')
        open (iua, file=cfile, form='formatted', status='new')

        write (cfile,'(3a)') 
     *    'mismip', runlab(1:lenchr(runlab)), 'vol.dat'
        open (iub, file=cfile, status='unknown')
        close (iub, status='delete')
        open (iub, file=cfile, form='formatted', status='new')
        write (iub,'(4a)') 
     *    ' timeice    xcengl', 
     *    '            vol(km3)', 
     *    '          volaf(km3)', 
     *    '          areag(km2)' 

        firstmismipplus = .false.
      endif

      do n=1,nglmax
        xgl(n) = vnull
        ygl(n) = vnull
        hgl(n) = vnull
        ubgl(n) = vnull
        uagl(n) = vnull
        usgl(n) = vnull
        vbgl(n) = vnull
        vagl(n) = vnull
        vsgl(n) = vnull
        ifugl(n) = nint(vnull)
      enddo

      ngl = 0

c>>>>>>>>>>>>>>>>
      do j=1,ny     ! u grid
      do i=1,nx-1   ! u grid
c>>>>>>>>>>>>>>>>
c--------------------------------------------------------
c#if defined (RESOL10) || defined (RESOL20)
        if (hgu(i,j).ne.0. .or. 
     *      maskwater(i,j).ne.maskwater(i+1,j)) then 

          if (hgu(i,j).ne.0.) then 
            zhgl = hgu(i,j)
            zugl = uschoof(i,j)
          else
            if (maskwater(i,j).eq.0) then
              zhgl = fracgu(i,j)*h(i,j)   + (1.-fracgu(i,j))*h(i+1,j)
            else
              zhgl = fracgu(i,j)*h(i+1,j) + (1.-fracgu(i,j))*h(i,j)
            endif
            zugl = ua(i,j)
          endif
c#else
cc      if (uschoof(i,j).ne.0. .and. hgu(i,j).ne.0.) then
c       if (hgu(i,j).ne.0.) then    ! uschoof may = 0 due to buttressing
c         zhgl = hgu(i,j)
c         zugl = uschoof(i,j)
c#endif
c--------------------------------------------------------

c         calc sia utop and uav (reduced by buttressing like ssa):
          zdfu = (rhoice*grav*zhgl*(hs(i+1,j)-hs(i,j))/dxu(i,j))**powi
          zui0 = crheolmismipplus * (-zdfu) * zhgl / (powi+1.)
          zuia = zui0 * (powi+1.)/(powi+2.)

          if (maskwater(i+1,j).eq.1) then 
            zthet = thetau(i+1,j)
          else
            zthet = thetau(i,j)
          endif
          if (abs(zugl).lt.0.5) zthet = 0. 

          zub = zugl - zuia*zthet
          zua = zugl
          zus = zugl + (zui0-zuia)*zthet

c--
c         zvb = 0.
c         zva = 0.
c         zvs = 0.
c--
c         i (on h-grid) presumed upstream of gl:
          jp = min (ny,j+1) 
          jm = max (1,j-1)
          zdfv =  (  rhoice*grav*zhgl
     *             * (hs(i,jp)-hs(i,jm))/(yh(i,jp)-yh(i,jm)) )**powi
          zvi0 = crheolmismipplus * (-zdfv) * zhgl / (powi+1.)
          zvia = zvi0 * (powi+1.)/(powi+2.)
          zvb = 0.5*(va(i,j-1)+va(i,j)) - zvia*zthet
          zva = 0.5*(va(i,j-1)+va(i,j))
          zvs = 0.5*(va(i,j-1)+va(i,j)) + (zvi0-zvia)*zthet
c--

          ngl = ngl + 1
          if (ngl.gt.nglmax) then
            write (6,'(a)') '*** Error (icemismipplus): ngl too large'
            stop
          endif
          xgl(ngl)  = 0.5*(xh(i,j)+xh(i+1,j))
          ygl(ngl)  = yh(i,j)
          hgl(ngl)  = zhgl
          ubgl(ngl) = zub
          uagl(ngl) = zua
          usgl(ngl) = zus
          vbgl(ngl) = zvb
          vagl(ngl) = zva
          vsgl(ngl) = zvs
          ifugl(ngl) = 1
c------------
        endif
c------------
c>>>>>>>>>>
      enddo
      enddo
c>>>>>>>>>>

c<<<<<<<<<<<<<<<<
      do i=1,nx       ! v grid
      do j=1,ny-1     ! v grid
c<<<<<<<<<<<<<<<<
c--------------------------------------------------------
c#if defined (RESOL10) || defined (RESOL20)
        if (hgv(i,j).ne.0. .or. 
     *      maskwater(i,j).ne.maskwater(i,j+1)) then 

          if (hgv(i,j).ne.0.) then 
            zhgl = hgv(i,j)
            zvgl = vschoof(i,j)
          else
            if (maskwater(i,j).eq.0) then
              zhgl = fracgv(i,j)*h(i,j)   + (1.-fracgv(i,j))*h(i,j+1)
            else
              zhgl = fracgv(i,j)*h(i,j+1) + (1.-fracgv(i,j))*h(i,j)
            endif
            zvgl = va(i,j)
          endif
c#else
cc      if (vschoof(i,j).ne.0. .and. hgv(i,j).ne.0.) then
c       if (hgv(i,j).ne.0.) then    ! vschoof may = 0 due to buttressing
c         zhgl = hgv(i,j)
c         zvgl = vschoof(i,j)
c#endif
c--------------------------------------------------------

c         calc sia vtop and vav (reduced by buttressing like ssa):
          zdfv = (rhoice*grav*zhgl*(hs(i,j+1)-hs(i,j))/dyv(i,j))**powi
          zvi0 = crheolmismipplus * (-zdfv) * zhgl / (powi+1.)
          zvia = zvi0 * (powi+1.)/(powi+2.)
          if (maskwater(i,j+1).eq.1) then 
            zthet = thetav(i,j+1)
          else
            zthet = thetav(i,j)
          endif
          if (abs(zvgl).lt.0.5) zthet = 0. 

          zvb = zvgl - zvia*zthet
          zva = zvgl
          zvs = zvgl + (zvi0-zvia)*zthet

c--
c         zub = 0.
c         zua = 0.
c         zus = 0.
c--
c         ju (on h-grid) presumed upstream of gl:
          if (j.le.ny/2) ju = j
          if (j.gt.ny/2) ju = j+1 
          ip = min (nx,i+1)
          im = max (1,i-1)
          zdfu =  (  rhoice*grav*zhgl
     *             * (hs(ip,ju)-hs(im,ju))/(xh(ip,ju)-xh(im,ju)) )**powi
          zui0 = crheolmismipplus * (-zdfu) * zhgl / (powi+1.)
          zuia = zui0 * (powi+1.)/(powi+2.)
          zub = 0.5*(ua(i-1,ju)+ua(i,ju)) - zuia*zthet
          zua = 0.5*(ua(i-1,ju)+ua(i,ju))
          zus = 0.5*(ua(i-1,ju)+ua(i,ju)) + (zui0-zuia)*zthet
c--

          ngl = ngl + 1
          if (ngl.gt.nglmax) then
            write (6,'(a)') '*** Error (icemismipplus): ngl too large'
            stop
          endif
          xgl(ngl)  = xh(i,j)
          ygl(ngl)  = 0.5*(yh(i,j)+yh(i,j+1))
          hgl(ngl)  = zhgl
          ubgl(ngl) = zub
          uagl(ngl) = zua
          usgl(ngl) = zus
          vbgl(ngl) = zvb
          vagl(ngl) = zva
          vsgl(ngl) = zvs
          ifugl(ngl) = 0
c------------
        endif
c------------
c>>>>>>>>>>
      enddo
      enddo
c>>>>>>>>>>

c        Sort ngl grounding line points into order

      call izero (indgl, ngl)
      call izero (ifsel, ngl)

c         Start at (one and only one) u-grid gl point with j=1:

      do n=1,ngl
        if (ifugl(n).eq.1 .and. ygl(n).le.0.6*dy0) then 
          ncur = n 
          goto 100
        endif
      enddo
      write (6,'(a)') '*** Error (icemismipplus): first pt not found'
      stop
  100 continue

      nind = 1
      indgl(nind) = ncur 
      ifsel(ncur) = 1 

c---------------------
      do iloop = 2,ngl
c---------------------

c        Find closest not-yet-selected gl point to current point

      zdistmin = 1.e20
      zymin = 9999999.
      do n=1,ngl
        if (ifsel(n).eq.0) then
          zdist = (xgl(n)-xgl(ncur))**2 + (ygl(n)-ygl(ncur))**2
c         give priority to points not yet selected with lesser y posn:
          if ( ( nint(ygl(n)).gt.nint(ygl(ncur)) .and. 
     *           nint(zymin).gt.nint(ygl(ncur)) .and.
     *           zdist.lt.zdistmin )
     *         .or.
     *         ( nint(ygl(n)).le.nint(ygl(ncur)) .and. 
     *           nint(zymin).gt.nint(ygl(ncur)) ) 
     *         .or.
     *         ( nint(ygl(n)).le.nint(ygl(ncur)) .and. 
     *           nint(zymin).le.nint(ygl(ncur)) .and.
     *           zdist.lt.zdistmin )
     *       ) then
            zdistmin = zdist
            zymin = ygl(n)
            nnext = n 
          endif
        endif
      enddo
      if (zdist.eq.1.e20) then  
        write (6,'(a)')'*** Error (icemismipplus): closest pt not found'
        stop
      endif

c        Add to sorted list

      nind = nind + 1
      indgl(nind) = nnext
      ifsel(nnext) = 1
      ncur = nnext

c----------
      enddo
c----------

c        Write all sorted gl points (u and v grid) to 1st output file

      write (iua,'(/a,f8.2, a,i6)') 
     *  'ALL: timeice=',timeice, '   ngl=',ngl

      write (iua,'(5a)') 
     *  '    m    n  ifu',
     *  '       xgl       ygl',
     *  '       hgl',
     *  '      ubgl      uagl      usgl',
     *  '      vbgl      vagl      vsgl'
 
      do m=1,ngl
        n = indgl(m)
        write (iua,'(3i5, 2f10.3, f10.3, 3f10.3, 3f10.3)')
     *    m, n, ifugl(n),
     *    xgl(n)*1.e-3, ygl(n)*1.e-3, 
     *    hgl(n), 
     *    ubgl(n), uagl(n), usgl(n),
     *    vbgl(n), vagl(n), vsgl(n)
      enddo

c        Write consecutive pair averages (combining staggered u,v's)

      write (iua,'(/a,f8.2, a,i6)') 
     *  'PAIRED: timeice=',timeice, '   ngl-1=',ngl-1

      write (iua,'(5a)') 
     *  '   m1   m2   n1   n2 ifu1 ifu2',
     *  '       xgl       ygl',
     *  '       hgl',
     *  '      ubgl      uagl      usgl',
     *  '      vbgl      vagl      vsgl'

      do m1=1,ngl-1
        m2 = m1 + 1
        n1 = indgl(m1)
        n2 = indgl(m2)

        zx = 0.5*(xgl(n1)+xgl(n2))
        zy = 0.5*(ygl(n1)+ygl(n2))
        zh = 0.5*(hgl(n1)+hgl(n2))
        zub = 0.5*(ubgl(n1)+ubgl(n2))
        zua = 0.5*(uagl(n1)+uagl(n2))
        zus = 0.5*(usgl(n1)+usgl(n2))
        zvb = 0.5*(vbgl(n1)+vbgl(n2))
        zva = 0.5*(vagl(n1)+vagl(n2))
        zvs = 0.5*(vsgl(n1)+vsgl(n2))

        write (iua,'(6i5, 2f10.3, f10.3, 3f10.3, 3f10.3)')
     *    m1, m2, n1, n2, ifugl(n1), ifugl(n2),
     *    zx*1.e-3, zy*1.e-3, 
     *    zh, 
     *    zub, zua, zus,
     *    zvb, zva, zvs
      enddo

c       Write 2-D printmap dumps (to 1st output file)

      do j=1,ny
        do i=1,nx
          work(i,j) = uschoof(i,j)
          work2(i,j) = hgu(i,j)
          work3(i,j) = 1. - thetau(i,j)
        enddo
      enddo
      call printmap (timeice,work, 'uschoof', 10.0, iua, 1)
      call printmap (timeice,work2, 'hgu', 150.0, iua, 1)
      call printmap (timeice, work3, '1-thetau', .05, iua, 0)

      do j=1,ny
        do i=1,nx
          work(i,j) = vschoof(i,j)
          work2(i,j) = hgv(i,j)
          work3(i,j) = 1. - thetav(i,j)
        enddo
      enddo
      call printmap (timeice,work, 'vschoof', 10.0, iua, 1)
      call printmap (timeice,work2, 'hgv', 150.0, iua, 1)
      call printmap (timeice, work3, '1-thetav', .05, iua, 0)

      call flush (iua)

c       Compute total quantities, write to 2nd output file 

c     find centerline gl position:
      xcengl = 0. 
      j = (ny+1)/2
      do i=1,nx-1
        if (maskwater(i,j).eq.0 .and. maskwater(i+1,j).eq.1) then 
          xcengl = 0.5*(xh(i,j) + xh(i+1,j))
          goto 200
        endif
      enddo
  200 continue

      zvol = 0.
      zvolaf = 0.
      zagr = 0.
      do j=1,ny
        do i=1,nx
          zvol = zvol + h(i,j)*darea(i,j)
          if (maskwater(i,j).eq.0) then
            zhf = max (0., sealev-hb(i,j)) / rhor
            zvolaf = zvolaf + darea(i,j) * max (0., h(i,j)-zhf)
            zagr = zagr + darea(i,j)
          endif
        enddo
      enddo
      write (iub,'(f8.2, f10.3, 3f20.6)') timeice, 
     *  xcengl*1.e-3, zvol*1.e-9, zvolaf*1.e-9, zagr*1.e-6  

      call flush (iub)

      return
      end

c}}}}}
#endif
c}}}}}
