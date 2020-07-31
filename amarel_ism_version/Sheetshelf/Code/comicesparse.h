c comicesparse.h

c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#if defined (DOSPARSE) || defined (MELANGE) || (defined (MOVEW) &&  defined (CHANNEL_SPARSE))
c|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

c     For sparse matrix soln (used in sparse solver and icedyn)

      parameter (nuvmax=2*nx*ny, nspamax=nuvmax*9 + 2)
      common /mat/ elspa(nspamax), ijspa(nspamax)

c     small point-number limit to use Gaussian elimination in icedyn
c     

#if defined (EISLINE) || defined (TRANSECTA)
c     parameter (nuvsmall = 2*nx*ny)
      parameter (nuvsmall = 20)
#else
c     parameter (nuvsmall = 6)
      parameter (nuvsmall = 22) ! for OCOH (no longer used)
#endif

c|||||
#endif
c|||||
