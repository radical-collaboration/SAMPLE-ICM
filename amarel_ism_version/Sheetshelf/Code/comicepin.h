c     Statement function (fpin) for fraction of u,v grid box with
c     sub-grid bedrock mounds scraping ice bottom:

c-----------------------
#if defined (SUBGRIDPIN)
c-----------------------

c     /200. ! /300. ! **2
      fpin(zhw,zhbsd) = max (0., min (1., 1.-(zhw/max(SUBGRIDPIN,1.)) ))

c----
#else
c----

#  if defined (SUBGRIDPINSD)

c     based on standard deviation of high-res bathymetry dataset
c     within eacn grid cell (only with SEARISEINIT, else hbsd=0):

      fpin (zhw,zhbsd) = 0.5 * max (0., min (1.,  
     *                         1. - zhw/(max(SUBGRIDPINSD*zhbsd,1.)) ))

#   else

      fpin (zhw,zhbsd) = 0.5 * max (0., min (1.,  
     *                         1. - zhw/(max(2.*zhbsd,1.)) ))

#   endif

c-----
#endif
c-----
