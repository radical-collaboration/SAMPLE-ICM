c comiceplot.h

c For communicating fields to be plotted, from various parts of 
c the code, to subr writehis (iceshow.F) which writes to
c Netcdf history file fort.92.nc, and to subr iceshow2d (iceshow.F)
c which writes to ascii file fort.19.

#if defined (CALVDAMAGE)

      common /ciceplot1/ 
     *  calvdam_t1(nx,ny), calvdam_t2(nx,ny), calvdam_t3(nx,ny)

#elif defined (CALVBALANCE)

      common /ciceplot1/ 
     *  calvbal_vi(nx,ny), calvbal_vb(nx,ny), calvbal_vc(nx,ny)

#endif
