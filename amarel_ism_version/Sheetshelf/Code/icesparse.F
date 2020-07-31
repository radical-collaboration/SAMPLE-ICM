c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#if defined (DOSPARSE) || defined (MELANGE) || ( defined (MOVEW) && defined (CHANNEL_SPARSE) )
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

c Melange always uses dosparsemelh.

c-----------------------------------------------------------------------

c=======================================================================
#if defined (NUMREC) ||defined (NUMRECMELU) ||defined (NUMRECMELH) || defined (CHANNEL_SPARSE)
c=======================================================================

c     Will need to fix duplicate linbcg, etc, with EBM's (-DEBM) 

      SUBROUTINE linbcg(n,b,x,itol,tol,itmax,iter,err)
#include <comicegrid.h>
#include <comicesparse.h>
      INTEGER iter,itmax,itol,n
c     REAL err,tol,b(*),x(*),EPS
      REAL err,tol,b(n),x(n),EPS
      PARAMETER (EPS=1.d-14)
CU    USES atimes,asolve,snrm
      INTEGER j
      REAL 
     *  ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm,
     *  p(nspamax),pp(nspamax),            ! nspamax from comicesparse.h
     *  r(nspamax),rr(nspamax),
     *  z(nspamax),zz(nspamax),
     *  snrm
      iter=0
      call atimes(n,x,r,0)
      do 11 j=1,n
        r(j)=b(j)-r(j)
        rr(j)=r(j)
11    continue
C     call atimes(n,r,rr,0)
      if(itol.eq.1) then
        bnrm=snrm(n,b,itol)
        call asolve(n,r,z,0)
      else if (itol.eq.2) then
        call asolve(n,b,z,0)
        bnrm=snrm(n,z,itol)
        call asolve(n,r,z,0)
      else if (itol.eq.3.or.itol.eq.4) then
        call asolve(n,b,z,0)
        bnrm=snrm(n,z,itol)
        call asolve(n,r,z,0)
        znrm=snrm(n,z,itol)
      else
        write (6,*) 'illegal itol in linbcg'
        stop
      endif
100   if (iter.le.itmax) then
        iter=iter+1
        call asolve(n,rr,zz,1)
        bknum=0.d0
        do 12 j=1,n
          bknum=bknum+z(j)*rr(j)
12      continue
        if(iter.eq.1) then
          do 13 j=1,n
            p(j)=z(j)
            pp(j)=zz(j)
13        continue
        else
          bk=bknum/bkden
          do 14 j=1,n
            p(j)=bk*p(j)+z(j)
            pp(j)=bk*pp(j)+zz(j)
14        continue
        endif
        bkden=bknum
        call atimes(n,p,z,0)
        akden=0.d0
        do 15 j=1,n
          akden=akden+z(j)*pp(j)
15      continue
        ak=bknum/akden
        call atimes(n,pp,zz,1)
        do 16 j=1,n
          x(j)=x(j)+ak*p(j)
          r(j)=r(j)-ak*z(j)
          rr(j)=rr(j)-ak*zz(j)
16      continue
        call asolve(n,r,z,0)
        if(itol.eq.1)then
          err=snrm(n,r,itol)/bnrm
        else if(itol.eq.2)then
          err=snrm(n,z,itol)/bnrm
        else if(itol.eq.3.or.itol.eq.4)then
          zm1nrm=znrm
          znrm=snrm(n,z,itol)
          if(abs(zm1nrm-znrm).gt.EPS*znrm) then
            dxnrm=abs(ak)*snrm(n,p,itol)
            err=znrm/abs(zm1nrm-znrm)*dxnrm
          else
            err=znrm/bnrm
            goto 100
          endif
          xnrm=snrm(n,x,itol)
          if(err.le.0.5d0*xnrm) then
            err=err/xnrm
          else
            err=znrm/bnrm
            goto 100
          endif
        endif
c       write (6,*) ' iter=',iter,' err=',err
      if(err.gt.tol) goto 100
      endif
      return
      END

c-----------------------------------------------------------------------

      SUBROUTINE atimes(n,x,r,itrnsp)
#include <comicegrid.h>
#include <comicesparse.h>
      INTEGER n,itrnsp                  
      REAL x(n),r(n)          
CU    USES dsprsax,dsprstx
      if (itrnsp.eq.0) then
        call dsprsax(elspa,ijspa,x,r,n) !elspa,ijspa from comicesparse.h
      else
        call dsprstx(elspa,ijspa,x,r,n)
      endif
      return
      END

c-----------------------------------------------------------------------

      SUBROUTINE asolve(n,b,x,itrnsp)
#include <comicegrid.h>
#include <comicesparse.h>         
      INTEGER n,itrnsp,i      
      REAL x(n),b(n)
      do 11 i=1,n
        x(i)=b(i)/elspa(i)             ! elspa,ijspa from comicesparse.h
11    continue
      return
      END

c-----------------------------------------------------------------------

      SUBROUTINE dsprsax(sa,ija,x,b,n)
      INTEGER n,ija(*)
      REAL b(n),sa(*),x(n)
      INTEGER i,k
      if (ija(1).ne.n+2) then
        write (6,*) 'mismatched vector and matrix in sprsax'
        stop
      endif
      do 12 i=1,n
        b(i)=sa(i)*x(i)
        do 11 k=ija(i),ija(i+1)-1
          b(i)=b(i)+sa(k)*x(ija(k))
11      continue
12    continue
      return
      END

c-----------------------------------------------------------------------

      SUBROUTINE dsprstx(sa,ija,x,b,n)
      INTEGER n,ija(*)
      REAL b(n),sa(*),x(n)
      INTEGER i,j,k
      if (ija(1).ne.n+2) then
        write (6,*) 'mismatched vector and matrix in sprstx'
        stop
      endif
      do 11 i=1,n
        b(i)=sa(i)*x(i)
11    continue
      do 13 i=1,n
        do 12 k=ija(i),ija(i+1)-1
          j=ija(k)
          b(j)=b(j)+sa(k)*x(i)
12      continue
13    continue
      return
      END

c-----------------------------------------------------------------------

      FUNCTION snrm(n,sx,itol)
      INTEGER n,itol,i,isamax
      REAL sx(n),snrm
      if (itol.le.3)then
        snrm=0.
        do 11 i=1,n
          snrm=snrm+sx(i)**2
11      continue
        snrm=sqrt(snrm)
      else
        isamax=1
        do 12 i=1,n
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i
12      continue
        snrm=abs(sx(isamax))
      endif
      return
      END

c=====
#endif
c=====

c-----------------------------------------------------------------------

c===========================================================
#if defined (SLAP) ||defined (SLAPMELU) ||defined (SLAPMELH)
c===========================================================

      subroutine slap (n,rhs,vec,itol,tol,itmax,iter,err)

#include <comicegrid.h>
#include <comicesparse.h>

      dimension rhs(nuvmax), vec(nuvmax)

      real*8 valnz(nspamax)
      dimension ianz(nspamax), janz(nspamax)

      parameter (nswork = nspamax + 8*nuvmax)
      dimension rswork(nswork), iswork(nswork)

c        Convert Numerical Recipes sparse storage to SLAP triad format

c        n       = matrix dimension  
c        nnz     = number of non-zero elements
c        valnz   = non-zero values (real*8), any order 
c        ianz    = row numbers of each valnz
c        janz    = column numbers of each valnz

      nnz = 0    

c     loop over rows ! new = corrected (w. Malte Heinemann 2/1/15)
      do i=1,n
        ka = ijspa(i)
c old   kb = ijspa(i+1) - 1 
        kb = ijspa(i+1)        ! new
        if (kb.ge.ka) then
c         diagonal:
          nnz = nnz + 1
          valnz(nnz) = elspa(i)
          ianz(nnz) = i
          janz(nnz) = i
c         off diagonals:
c old     do k=ka,kb
          if (kb-1.ge.ka) then   ! new
            do k=ka,kb-1         ! new
              nnz = nnz + 1
              valnz(nnz) = elspa(k)
              ianz(nnz) = i
              janz(nnz) = ijspa(k)
            enddo
          endif                  ! mew
        else
          write(6,'(/a,3i8)') 'ERROR (slap): i, ijspa(i), ijspa(i+1)=',
     *                        i, ijspa(i), ijspa(i+1)
          stop
        endif
      enddo

c     write (6,*) 'slap: n,nnz = ',n,nnz

c        Call SLAP solver 

      isym = 0
      iunit = 0

      call dslucs
     *  (n, rhs, vec, nnz, ianz, janz, valnz, isym, itol, tol, itmax, 
     *   iter, err, ierr, iunit, rswork, nswork, iswork, nswork)

      return
      end

c=====
#endif
c=====

c-----------------------------------------------------------------------

c====================================================================
#if defined (SUPERLU) ||defined (SUPERLUMELU) ||defined (SUPERLUMELH)
c====================================================================

      subroutine superlu (n, rhs, vec, info)

#include <comicegrid.h>
#include <comicesparse.h>

      dimension rhs(nuvmax), vec(nuvmax)

      real*8 valnz(nspamax)
      dimension icolind(nspamax), irowptr(nuvmax+1)   

c        Convert Numerical Recipes sparse storage to Superlu in
c        compressed *row* (not column) order. 

c        n       = matrix dimension  
c        nnz     = number of non-zero elements
c        valnz   = non-zero values (real*8), grouped by row (but within 
c                  each row-group, elements can be in any order)
c        icolind = column numbers of each valnz
c        irowptr = index of valnz and icolind where each row-group 
c                  starts (for irowptr(1-n)). irowptr(n+1) = nnz + 1.

      nnz = 0    

c     loop over rows
      do i=1,n
        ka = ijspa(i)
        kb = ijspa(i+1)
        if (kb.ge.ka) then
c         diagonal:
          nnz = nnz + 1
          valnz(nnz) = elspa(i)
          icolind(nnz) = i
          irowptr(i) = nnz
c         off diagonals:
          if (kb-1.ge.ka) then
            do k=ka,kb-1
              nnz = nnz + 1
              valnz(nnz) = elspa(k)
              icolind(nnz) = ijspa(k)
            enddo
          endif
        endif
      enddo
      irowptr(n+1) = nnz + 1

c     write (6,*) 'superlu: n,nnz = ',n,nnz

c        Call superlu solver 

      nrhs = 1
      ldb = n
      nprocs = 1
      do i=1,n
        vec(i) = rhs(i)
      enddo

      call c_bridge_pdgssv
     * (nprocs, n, nnz, nrhs, valnz, icolind, irowptr, vec, ldb, info)

      return
      end

c=====
#endif
c=====

c-----------------------------------------------------------------------

c=================================================================
#if defined (WATSON) ||defined (WATSONMELU) ||defined (WATSONMELH)
c=================================================================

      subroutine watson (n, rhs, vec, info)

#include <comicegrid.h>
#include <comicesparse.h>

      integer*4 n, icolind, irowptr, nnz, info, ldb, nrhs, niter

      dimension rhs(nuvmax), vec(nuvmax)

      real*8 valnz(nspamax)
      dimension icolind(nspamax), irowptr(nuvmax+1)   
      real*8 berr, thresh, opc

      real*8 rmisc (nuvmax,1)
      integer iparm(64)
      real*8 dparm(64)

      logical first
      data first /.true./
      save first
      save iparm, dparm, rmisc


c        Convert Numerical Recipes sparse storage to WSMP 
c        (compressed row order). 

c        n       = matrix dimension  
c        nnz     = number of non-zero elements
c        valnz   = non-zero values (real*8), grouped by row (but within 
c                  each row-group, elements can be in any order)
c        icolind = column numbers of each valnz
c        irowptr = index of valnz and icolind where each row-group 
c                  starts (for irowptr(1-n)). irowptr(n+1) = nnz + 1.

c     call wsetmaxthrds(1)

      nnz = 0    

c     loop over rows
      do i=1,n
        ka = ijspa(i)
        kb = ijspa(i+1)
        if (kb.ge.ka) then
c         diagonal:
          nnz = nnz + 1
          valnz(nnz) = elspa(i)
          icolind(nnz) = i
          irowptr(i) = nnz
c         off diagonals:
          if (kb-1.ge.ka) then
            do k=ka,kb-1
              nnz = nnz + 1
              valnz(nnz) = elspa(k)
              icolind(nnz) = ijspa(k)
            enddo
          endif
        endif
      enddo
      irowptr(n+1) = nnz + 1

      call ws_sortindices_d (n, n, irowptr, icolind, valnz, info)
      if (info.ne.0) then
        write (6,*) 'ws_sortindices: info=',info
        stop
      endif

      ldb = n
      nrhs = 1

c     if (first) then
        do i=1,64
          iparm(i) = 0
          dparm(i) = 0.
        enddo
c       call wgsmp (n, irowptr, icolind, valnz, vec, ldb, nrhs, rmisc,
c    *              iparm, dparm)
c       first = .false.
c     endif

      do i = 1, n
        vec(i) = rhs(i)
      end do
c     iparm(1) = 1
      iparm(2) = 1
      iparm(3) = 4
c     iparm(11) = 1
c     iparm(12) = 0
      call wgsmp (n, irowptr, icolind, valnz, vec, ldb, nrhs, rmisc,
     *            iparm, dparm)
 
      if (iparm(64) .ne. 0) then
        print *,'wgsmp: iparm(2), iparm(64)=',iparm(2), iparm(64)
        stop
      end if

      return


c     write (6,*) 'watson: n,nnz = ',n,nnz

c        Call WSMP solver 

      ldb = n
      nrhs = 1
      niter = 3
      thresh = .05

c.. Analysis.

c     waltime = rtc()
      call wgralz (n, irowptr, icolind, valnz, nnz, opc, info)      

c     print *,'Analysis complete in time - ',rtc()-waltime

      if (info .ne. 0) then
        print *,'The following analysis ERROR was detected: ',info
        stop
      end if

c     print *,'Number of nonzeros in LU factors = ',nnz
c     print *,'Number of FLOPS in factorization = ',opc

C.. Factorization.

c     waltime = rtc()
      call wgrluf (n, irowptr, icolind, valnz, thresh, info)

c     waltime = rtc() - waltime
c     print *,'Factorization complete in time - ',waltime
      if (info .ne. 0) then
        print *,'The following factorization ERROR was detected: ',info
        stop
      end if
c     print *,'Factorization MegaFlops = ',(opc*1.d-6)/waltime

C.. Back substitution with iterative refinement.

c     waltime = rtc()
      do i=1,n
c       vec(i) = rhs(i)
        vec(i) = 1.
      enddo
      call wgrslv (n, irowptr, icolind, valnz, vec, ldb, nrhs, niter,
     *             berr, info)

c     print *,'Back substitution complete in time - ',rtc()-waltime
      if (info .ne. 0) then
        print *,'The following back subst ERROR was detected: ',info
        stop
      end if

c     print *,'Maximum relative error = ',berr
c     print *,'The solution of the system is as follows:'
c     do i = 1, n
c       print *,i,' : ',vec(i)
c     end do

      return
      end

c=====
#endif
c=====

c-----------------------------------------------------------------------

c========================================================
#if defined (MKL) ||defined (MKLMELU) ||defined (MKLMELH)
c========================================================

      subroutine mkl (n, rhs, vec, info)

#include <comicegrid.h>
#include <comicesparse.h>

      dimension rhs(nuvmax), vec(nuvmax)

      real*8 valnz(nspamax)
      dimension icolind(nspamax), irowptr(nuvmax+1)   

      integer*8 memptr(64) 
      dimension iparm(64)

c        Convert Numerical Recipes sparse storage to MKL
c        (compressed row order). 

c        n       = matrix dimension  
c        nnz     = number of non-zero elements
c        valnz   = non-zero values (real*8), grouped by row (but within 
c                  each row-group, elements in any order, sorted below)
c        icolind = column numbers of each valnz
c        irowptr = index of valnz and icolind where each row-group 
c                  starts (for irowptr(1-n)). irowptr(n+1) = nnz + 1.

      logical first
      save first
      data first /.true./

      nnz = 0    

c     loop over rows
      do i=1,n
c       diagonal:
        nnz = nnz + 1
        valnz(nnz) = elspa(i)
        icolind(nnz) = i
        irowptr(i) = nnz

c       off diagonals:
        ka = ijspa(i)
        kb = ijspa(i+1)
        if (kb-1.ge.ka) then
          do k=ka,kb-1
            nnz = nnz + 1
            valnz(nnz) = elspa(k)
            icolind(nnz) = ijspa(k)
          enddo
        endif
      enddo
      irowptr(n+1) = nnz + 1

c        Sort columns into ascending order within each row

      do i = 1,n
        if (irowptr(i).gt.0) then
          i1 = irowptr(i)  
          i2 = irowptr(i+1)-1  
          call mkl_bubblesort (icolind(i1), valnz(i1), i2-i1+1)   
        endif
      enddo

c        Set up for mkl solver 

      if (first) then
        do i=1,64
          memptr(i) = 0
        enddo
        first = .false.
      endif

      maxfct = 1
      mnum = 1
      mtype = 11
      iphase = 13 
      nrhs = 1
      do i=1,64
        iparm(i) = 0
      enddo
      msglvl = 0
      ierr = 0

c        Call mkl solver 

      call pardiso (memptr, maxfct, mnum, mtype, iphase, 
     *              n, valnz, irowptr, icolind, 
     *              idum, nrhs, iparm, msglvl, rhs, vec, ierra)
      if (ierra.ne.0) then
        write (6,*) 'pardiso (solve) error: ierra=',ierra
        stop
      endif

c        Release memory

      iphase = -1
      call pardiso (memptr, maxfct, mnum, mtype, iphase, 
     *              n, ddum, idum, idum, 
     *              idum, nrhs, iparm, msglvl, ddum, ddum, ierrb)
      if (ierrb.ne.0) then
        write (6,*) 'pardiso (memory release) error: ierrb=',ierrb
        stop
      endif

      info = 10*ierrb + ierra
      
      return
      end
c
c----------------------------------------------------------------------
c
      subroutine mkl_bubblesort (ia, va, n)

c        Bubble sorts ia(1:n) into ascending order, with va following

      dimension ia(n), va(n)

      do ibub = n-1,1,-1
        do j=ibub,n-1
          if (ia(j).gt.ia(j+1)) then
            iatmp = ia(j+1)
            ia(j+1) = ia(j)
            ia(j) = iatmp

            vatmp = va(j+1)
            va(j+1) = va(j)
            va(j) = vatmp
          else
            go to 10
          endif
        enddo
   10   continue
      enddo

      return
      end

c=====
#endif
c=====

c-----------------------------------------------------------------------

      subroutine gaussdo (nuvtot, rhs, vec, info) 

c     Sets up to do Gaussian elimination. Always called (regardless
c     of defines) if nuvtot <= nuvsmall, to avoid sparse-solution
c     problems for nuvtot <= ~2.  

#include <comicegrid.h>
#include <comicesparse.h>

      dimension rhs(nuvmax), vec(nuvmax)

c     If not defined GAUSS, gaussdo only called for nuvtot <= nuvsmall

#if defined (GAUSS) ||defined (GAUSSMELU) ||defined (GAUSSMELH)
      parameter (nuvgauss = nuvmax)
#else
      parameter (nuvgauss = nuvsmall)
#endif
      dimension a(nuvgauss,nuvgauss), asav(nuvgauss,nuvgauss)

      info = 0

      call zero (a, nuvgauss*nuvgauss)

      do i=1,nuvtot 
c       diagonal:
        a(i,i) = elspa(i)
      enddo

      do i=1,nuvtot
c       off diagonals:
        ka = ijspa(i)
        kb = ijspa(i+1)
        if (kb-1.ge.ka) then
          do k=ka,kb-1
            a(i,ijspa(k)) = elspa(k)
          enddo
        endif
      enddo

c       print out in 2-D matrix form

#if defined (MATPRINT) 
      if (.true.) then
        do iloop=1,2
          if (iloop.eq.1) iu = 6
          if (iloop.eq.2) iu = 20
          write (iu,*) 'nuvtot=', nuvtot

          write (iu,*) 'a:'
          do i=1,nuvtot
            write (iu,'(i6, 151f15.5)') i, (a(i,j), j=1,nuvtot)
          enddo 

          write (iu,*) 'rhs :'
          write (iu,'(6x, 151f15.5)')  (rhs (i), i=1,nuvtot)
        enddo
      endif
#endif

c       solve for vec by Gaussian elimination 

      do i=1,nuvtot
        vec(i) = rhs(i)
      enddo

      call gaussj (a, asav, nuvtot, nuvgauss, vec, 1, 1)

#if defined (MATPRINT)
      do iloop=1,2
        if (iloop.eq.1) iu = 6
        if (iloop.eq.2) iu = 20
        write (iu,*) 'vec:'
        write (iu,'(6x, 151f15.5)')  (vec(i), i=1,nuvtot)
      enddo
#endif

      return
      end

c-----------------------------------------------------------------------

      SUBROUTINE gaussj (a,asav,n,np,b,m,mp)

c     Solves A x = b using Gaussian Elimination 
c     (from Numerical Recipes?). Returns solution in b.  

      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp), asav(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      character*11 cline(900)

      do j=1,n
        do i=1,n
          asav(i,j) = a(i,j)
        enddo
      enddo

      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
c               pause            'singular matrix (i) in gaussj'   ! old
                write (6,'(/a)') 'singular matrix (i) in gaussj'
c               dump matrix a and rhs b:
c               write (6,'(a/8x,500i10)') 'a(n,n), b(n):', (ii,ii=1,n)
c               do jj=1,n
c                 write (6,'(i5,3x,500f10.4)') 
c    *              jj, (asav(ii,jj),ii=1,n), b(jj,1) 
c               enddo
                call exit (-1)
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) then
          do iloop=1,2
            if (iloop.eq.1) iu = 6
            if (iloop.eq.2) iu = 20
            write (iu,*) 'a (singular ii) irow,icol=',irow,icol
            do io=1,n
              write (iu,'(a,i4,a,f15.8)')'i=',io,'  a(i,i)=',asav(io,io)
            enddo
            do iseg=1,(n-1)/100 + 1
             joa = (iseg-1)*100 + 1
             job = min (iseg*100, n)
             write (iu,'(/6x, 900i11)') (jo, jo=joa,job)
             do io=1,n
               do jo=joa,job
                 if (abs(asav(io,jo)).le.1.e-10) then
                   cline(jo) = '          Z'
                 else
                   write (cline(jo),'(f11.6)') asav(io,jo)
                 endif
                 if (io.eq.164) cline(jo)(2:2) = 'X'
                 if (io.eq.jo ) cline(jo)(2:2) = '*'
               enddo
               write (iu,'(i6, 900a11)') io, (cline(jo),jo=joa,job)
             enddo
             write (iu,'(6x, 900i11/)') (jo, jo=joa,job)
           enddo
          enddo
c         pause            'singular (ii) matrix in gaussj'        ! old
          write (6,'(/a)') 'singular (ii) matrix in gaussj'
          call exit (-1)
        endif
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

c-----------------------------------------------------------------------

c|||||
#endif
c|||||
