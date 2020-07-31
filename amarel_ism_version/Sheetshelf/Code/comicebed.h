c-------------------------------------------
#if defined (EISLINE) &&  defined (LINEF_OLD)
c-------------------------------------------

c     Statement function for bedrock profile versus x (m)

cc    bedprofile(x) = -x*.0006
cc    bedprofile(x) = max (-x*.0006, -360.)
cc    bedprofile(x) = min (200., max (-500., (350.e3-x)*.005))
c     bedprofile(x) = min (200., max (-1500.,(350.e3-x)*.005))
      bedprofile(x) = min (100., max (-1500., 
     *                100. - 200.*(x-600.e3)/100.e3 ))

c     Schoof:
cc    bedprofile(x) = 729. - 2184.8  * (x/750.e3)**2
cc   *                     + 1031.72 * (x/750.e3)**4
cc   *                     - 151.72  * (x/750.e3)**6

c------------------------------------------
#elif defined (EISLINE) &&  defined (LINEM)
c------------------------------------------

c     Statement function for bedrock profile versus x (m)

#  if defined (EXPT1)
      bedprofile(x) = -100. - (x/1.e3)
#  elif defined (EXPT3)
      bedprofile(x) = 729. - 2184.8  * (x/750.e3)**2
     *                     + 1031.72 * (x/750.e3)**4
     *                     - 151.72  * (x/750.e3)**6
#  endif

c-----
#endif
c-----
