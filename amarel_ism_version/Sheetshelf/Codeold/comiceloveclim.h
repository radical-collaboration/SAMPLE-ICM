c comiceloveclim.h

c+++++++++++++++++++++
#if defined (LOVECLIP)
c+++++++++++++++++++++

c     for Loveclim <-> ice model asynchronous coupling 
c     (separate programs, file passing, overarching script).

c     climfilein is read at start of ice run, with Loveclim climate,etc.

c     climfileout is written at end of ice run, with ice sheet info.

c     climfilectl is like climfilein, but for Loveclim control (modern)
c     climate, for anomaly method (-DLOVEANOM). 

      common /cloveclip/ climfilein, climfileout, climfilectl
      character*240 climfilein, climfileout, climfilectl

      parameter (vnull_lc = -1.e34)             ! nb: real*8, not real*4

c+++++
#endif
c+++++
