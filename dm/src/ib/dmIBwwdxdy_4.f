*****************************************************************************
*** auxiliary routine called by dsIBwwdxdy
*** author: Torsten Bringmann, 2007-02-16
*****************************************************************************

      real*8 function dmIBwwdxdy_4(x,y,m0,mw,mc1,mc2)
      real*8 x,y,m0,mw,mc1,mc2

      dmIBwwdxdy_4 = 
     -   (-256*m0**2*mc1*mc2*(-mw**14 + 
     -      1024*m0**14*(-1 + x)*(x - y)*y*
     -       (1 + 2*x*(-1 + y) - 2*y**2)**2 + 
     -      m0**2*mw**8*(mc1**2*
     -          (-48*mc2**2*(-1 + x) + 4*mw**2*(-11 + 13*x)) + 
     -         mw**2*(mc2**2*(-44 + 52*x) + mw**2*(41 - 69*x + 24*y)))
     -       - 4*m0**4*mw**6*
     -       (2*mc1**2*(8*mc2**2*(1 - 4*x + x**2 + 6*y - 6*x*y) + 
     -            mw**2*(-13 + 45*x - 26*x**2 - 40*y + 56*x*y)) + 
     -         mw**2*(-2*mc2**2*(13 - 45*x + 26*x**2 + 40*y - 56*x*y) + 
     -            mw**2*(34 - 117*x + 108*x**2 + 70*y - 194*x*y + 
     -               60*y**2))) - 
     -      64*m0**10*(-16*mc1**2*mc2**2*(-1 + x)*(x - y)*y - 
     -         4*(mc1**2 + mc2**2)*mw**2*
     -          (4*x**5 - 2*x**4*(1 + 6*y) + 2*x**3*y*(5 + 8*y) + 
     -            2*y*(-1 - 4*y + 4*y**2 + 2*y**3) - 
     -            2*x**2*(1 + 8*y - y**2 + 4*y**3) + 
     -            x*(1 + 14*y + 4*y**2 - 16*y**3 + 4*y**4)) + 
     -         mw**4*(-1 + 64*x**5 - 24*y - 36*y**2 + 96*y**3 - 
     -            28*y**4 - 96*y**5 - 4*x**4*(10 + 69*y) + 
     -            4*x**3*(4 + 57*y + 108*y**2) - 
     -            4*x**2*(9 + 41*y + 16*y**2 + 114*y**3) + 
     -            x*(17 + 124*y - 44*y**2 - 40*y**3 + 332*y**4))) - 
     -      256*m0**12*(-4*(mc1**2 + mc2**2)*y*
     -          (2*x**3*(-1 + y) + y - 2*y**3 + x**2*(3 - 4*y**2) + 
     -            x*(-1 - 3*y + 4*y**2 + 2*y**3)) + 
     -         mw**2*(8*x**6 - 8*x**5*(1 + 4*y) + 
     -            x**4*(2 + 32*y + 64*y**2) - 
     -            4*x**3*(1 + 10*y - 3*y**2 + 26*y**3) + 
     -            4*x**2*(1 + 13*y - 10*y**2 - 8*y**3 + 32*y**4) + 
     -            2*y*(1 + 6*y - 8*y**2 - 12*y**3 + 12*y**4 + 8*y**5)- 
     -          x*(1 + 22*y + 20*y**2 - 96*y**3 + 36*y**4 + 72*y**5)))
     -        - 16*m0**6*mw**4*
     -       (2*mc1**2*(mw**2*
     -             (-18*x**3 + 5*x**2*(3 + 16*y) + 
     -               4*(1 + 6*y + 3*y**2) - 4*x*(5 + 15*y + 9*y**2)) + 
     -            2*mc2**2*(-1 + 4*x**3 - 8*y - 20*x**2*y - 12*y**2 + 
     -               x*(5 + 20*y + 12*y**2))) + 
     -         mw**2*(mc2**2*
     -             (-36*x**3 + 10*x**2*(3 + 16*y) + 
     -               8*(1 + 6*y + 3*y**2) - 8*x*(5 + 15*y + 9*y**2)) + 
     -            mw**2*(-14 + 89*x**3 - 56*y - 15*y**2 - 80*y**3 - 
     -               x**2*(83 + 323*y) + x*(68 + 167*y + 231*y**2))))+ 
     -      64*m0**8*mw**2*(mc1**2*
     -          (-4*mc2**2*(2*x**4 - 4*x**3*y + 4*x**2*y*(1 + y) + 
     -               2*y*(1 + 2*y) - x*(1 + 6*y + 4*y**2)) + 
     -            mw**2*(1 + 20*x**4 + 16*y + 16*y**2 - 16*y**3 - 
     -               4*x**3*(1 + 19*y) + x**2*(8 + 60*y + 68*y**2) - 
     -               x*(11 + 48*y + 16*y**2 + 16*y**3))) + 
     -         mw**2*(mc2**2*
     -             (1 + 20*x**4 + 16*y + 16*y**2 - 16*y**3 - 
     -               4*x**3*(1 + 19*y) + x**2*(8 + 60*y + 68*y**2) - 
     -               x*(11 + 48*y + 16*y**2 + 16*y**3)) + 
     -            mw**2*(-3 - 50*x**4 - 26*y - 4*y**2 + 12*y**3 - 
     -               60*y**4 + x**3*(29 + 222*y) - 
     -               2*x**2*(15 + 74*y + 126*y**2) + 
     -               2*x*(12 + 39*y + 17*y**2 + 86*y**3))))))/
     -  (mw**4*(-2*mc1**2 + 3*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (-2*mc2**2 + 3*mw**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (mw**2 + 4*m0**2*(x - y))**2*(mw**2 - 4*m0**2*y)**2*
     -    (-2*mc1**2 + mw**2 + m0**2*(-2 + 4*y))*
     -    (-2*mc2**2 + mw**2 + m0**2*(-2 + 4*y)))

      return
      end   ! dmIBwwdxdy_4
