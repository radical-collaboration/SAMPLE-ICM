c--------------------
#if defined (MELANGE)
c--------------------

c     Statement function for melange internal pressure vs thickness
c     (interpreting thickness as density of big berg-like pieces)...
c     cf. Flato & Hibler P* h exp (-K(1-A)), with P* = 2.7e4 N m-2.

c     Used where "rhomp" is used in icedyn.F (schoofgl, setbpumel)
c     and icemelange.F (thetacalcmel, checkmelforce[2]). Change 
c     rhomp*grav*hm**2 to 0.5*rhomp*grav*hm**2 + presmel(hm).
c     Also, in icemelange.F for do... veloc solvers, change
c     dfu = rhomel*grav*hu*dhsu to rhomel*grav*hu*dhsu + d(presmel)/dx.
c     See notes 12/7/17.

#  if defined (HPRESMEL)
      presmel(h) = rhomp*grav*(HPRESMEL**2)    ! steep rise for h 
     *           * (max(h-HPRESMEL,0.)/10.)    ! > HPRESMEL (by ~10 m's)
c b  *           * (max(h-HPRESMEL,0.)/30.)    ! > HPRESMEL (by ~10 m's)
c a  *           * (max(h-HPRESMEL,0.)/10.)**3 ! > HPRESMEL (by ~10 m's)
#  else
      presmel(h) = 0.
#  endif

c-----
#endif
c-----
