***
      real*8 FUNCTION lzamsf(m)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
      real*8 m,mx,a(200)
      common /MSCFF/ a
*
* A function to evaluate Lzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      mx = SQRT(m)
      lzamsf = (a(1)*m**5*mx + a(2)*m**11)/
     &         (a(3) + m**3 + a(4)*m**5 + a(5)*m**7 +
     &          a(6)*m**8 + a(7)*m**9*mx)
*
      return
      end
***
      real*8 FUNCTION rzamsf(m)
      implicit none
      real*8 m,mx,a(200)
      common /MSCFF/ a
*
* A function to evaluate Rzams
* ( from Tout et al., 1996, MNRAS, 281, 257 ).
*
      mx = SQRT(m)
      rzamsf = ((a(8)*m**2 + a(9)*m**6)*mx + a(10)*m**11 +
     &          (a(11) + a(12)*mx)*m**19)/
     &         (a(13) + a(14)*m**2 +
     &          (a(15)*m**8 + m**18 + a(16)*m**19)*mx)
*
      return
      end
***
      real*8 FUNCTION tbgbf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the lifetime to the BGB or to
* Helium ignition if no FGB exists.
* (JH 24/11/97)
*
      tbgbf = (a(17) + a(18)*m**4 + a(19)*m**(11.d0/2.d0) + m**7)/
     &        (a(20)*m**2 + a(21)*m**7)
*
      return
      end
***
      real*8 FUNCTION tbgbdf(m)
      implicit none
      real*8 m,mx,f,df,g,dg,a(200)
      common /MSCFF/ a
*
* A function to evaluate the derivitive of the lifetime to the BGB
* (or to Helium ignition if no FGB exists) wrt mass.
* (JH 24/11/97)
*
      mx = SQRT(m)
      f = a(17) + a(18)*m**4 + a(19)*m**5*mx + m**7
      df = 4.d0*a(18)*m**3 + 5.5d0*a(19)*m**4*mx + 7.d0*m**6
      g = a(20)*m**2 + a(21)*m**7
      dg = 2.d0*a(20)*m + 7.d0*a(21)*m**6
      tbgbdf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION tbgdzf(m)
      implicit none
      real*8 m,mx,f,df,g,dg,a(200)
      common /MSCFF/ a
*
* A function to evaluate the derivitive of the lifetime to the BGB
* (or to Helium ignition if no FGB exists) wrt Z.
* (JH 14/12/98)
*
      mx = m**5*SQRT(m)
      f = a(17) + a(18)*m**4 + a(19)*mx + m**7
      df = a(117) + a(118)*m**4 + a(119)*mx
      g = a(20)*m**2 + a(21)*m**7
      dg = a(120)*m**2
      tbgdzf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION thookf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the lifetime to the end of the MS
* hook ( for those models that have one ) as a fraction of
* the lifetime to the BGB
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      thookf = 1.d0 - 0.01d0*MAX(a(22)/m**a(23),a(24)+a(25)/m**a(26))
      thookf = MAX(thookf,0.5d0)
*
      return
      end
***
      real*8 FUNCTION ltmsf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the luminosity at the end of the MS
* (JH 24/11/97)
*
      ltmsf = (a(27)*m**3 + a(28)*m**4 + a(29)*m**(a(32)+1.8d0))/
     &        (a(30) + a(31)*m**5 + m**a(32))
*
      return
      end
***
      real*8 FUNCTION lalphf(m)
      implicit none
      real*8 m,mcut,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity alpha coefficent.
* (JH 24/11/97)
*
      mcut = 2.d0
      if(m.ge.mcut)then
         lalphf = (a(33) + a(34)*m**a(36))/(m**0.4d0 + a(35)*m**1.9d0)
      else
         if(m.le.0.5d0)then
            lalphf = a(39)
         elseif(m.le.0.7d0)then
            lalphf = a(39) + ((0.3d0 - a(39))/0.2d0)*(m - 0.5d0)
         elseif(m.le.a(37))then
            lalphf = 0.3d0 + ((a(40)-0.3d0)/(a(37)-0.7d0))*(m - 0.7d0)
         elseif(m.le.a(38))then
            lalphf = a(40) + ((a(41)-a(40))/(a(38)-a(37)))*(m - a(37))
         else
            lalphf = a(41) + ((a(42)-a(41))/(mcut-a(38)))*(m - a(38))
         endif
      endif
*
      return
      end
***
      real*8 FUNCTION lbetaf(m)
      implicit none
      real*8 m,a1,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity beta coefficent.
* (JH 24/11/97)
*
      lbetaf = a(43) - a(44)*m**a(45)
      lbetaf = MAX(lbetaf,0.d0)
      if(m.gt.a(46).and.lbetaf.gt.0.d0)then
         a1 = a(43) - a(44)*a(46)**a(45)
         lbetaf = a1 - 10.d0*a1*(m - a(46))
         lbetaf = MAX(lbetaf,0.d0)
      endif
*
      return
      end
***
      real*8 FUNCTION lnetaf(m)
      implicit none
      real*8 m,a(200)
      common /MSCFF/ a
*
* A function to evaluate the Luminosity neta exponent.
* (JH 24/11/97)
*
      if(m.le.1.d0)then
         lnetaf = 10.d0
      elseif(m.ge.1.1d0)then
         lnetaf = 20.d0
      else
         lnetaf = 10.d0 + 100.d0*(m - 1.d0)
      endif
      lnetaf = MIN(lnetaf,a(97))
*
      return
      end
***
      real*8 FUNCTION lhookf(m,mhook)
      implicit none
      real*8 m,mhook,a2,a(200)
      common /MSCFF/ a
*
* A function to evalute the luminosity at the start of
* the MS hook ( for those stars that have one ).
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      if(m.le.mhook)then
         lhookf = 0.d0
      elseif(m.ge.a(51))then
         lhookf = MIN(a(47)/m**a(48),a(49)/m**a(50))
      else
         a2 = MIN(a(47)/a(51)**a(48),a(49)/a(51)**a(50))
         lhookf = a2*((m-mhook)/(a(51)-mhook))**0.4d0
      endif
*
      return
      end
***
      real*8 FUNCTION rtmsf(m)
      implicit none
      real*8 m,m2,rchk,a(200)
      common /MSCFF/ a
      real*8 rzamsf
      external rzamsf
*
* A function to evaluate the radius at the end of the MS
* Note that a safety check is added to ensure Rtms > Rzams
* when extrapolating the function to low masses.
* (JH 24/11/97)
*
      m2 = a(62) + 0.1d0
      if(m.le.a(62))then
         rchk = 1.5d0*rzamsf(m)
         rtmsf = MAX(rchk,(a(52) + a(53)*m**a(55))/(a(54) + m**a(56)))
      elseif(m.ge.m2)then
         rtmsf = (a(57)*m**3+a(58)*m**a(61)+a(59)*m**(a(61)+1.5d0))/
     &           (a(60) + m**5)
      else
         rtmsf = a(63) + ((a(64) - a(63))/0.1d0)*(m - a(62))
      endif
*
      return
      end
***
      real*8 FUNCTION ralphf(m)
      implicit none
      real*8 m,a5,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius alpha coefficent.
* (JH 24/11/97)
*
      if(m.le.0.5d0)then
         ralphf = a(73)
      elseif(m.le.0.65d0)then
         ralphf = a(73) + ((a(74) - a(73))/0.15d0)*(m - 0.5d0)
      elseif(m.le.a(70))then
         ralphf = a(74) + ((a(75)-a(74))/(a(70)-0.65d0))*(m - 0.65d0)
      elseif(m.le.a(71))then
         ralphf = a(75) + ((a(76) - a(75))/(a(71) - a(70)))*(m - a(70))
      elseif(m.le.a(72))then
         ralphf = (a(65)*m**a(67))/(a(66) + m**a(68))
      else
         a5 = (a(65)*a(72)**a(67))/(a(66) + a(72)**a(68))
         ralphf = a5 + a(69)*(m - a(72))
      endif
*
      return
      end
***
      real*8 FUNCTION rbetaf(m)
      implicit none
      real*8 m,m2,m3,b2,b3,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius beta coefficent.
* (JH 24/11/97)
*
      m2 = 2.d0
      m3 = 16.d0
      if(m.le.1.d0)then
         rbetaf = 1.06d0
      elseif(m.le.a(82))then
         rbetaf = 1.06d0 + ((a(81)-1.06d0)/(a(82)-1.d0))*(m-1.d0)
      elseif(m.le.m2)then
         b2 = (a(77)*m2**(7.d0/2.d0))/(a(78) + m2**a(79))
         rbetaf = a(81) + ((b2-a(81))/(m2-a(82)))*(m-a(82))
      elseif(m.le.m3)then
         rbetaf = (a(77)*m**(7.d0/2.d0))/(a(78) + m**a(79))
      else
         b3 = (a(77)*m3**(7.d0/2.d0))/(a(78) + m3**a(79))
         rbetaf = b3 + a(80)*(m - m3)
      endif
      rbetaf = rbetaf - 1.d0
*
      return
      end
***
      real*8 FUNCTION rgammf(m)
      implicit none
      real*8 m,m1,b1,a(200)
      common /MSCFF/ a
*
* A function to evaluate the radius gamma coefficent.
* (JH 24/11/97)
*
      m1 = 1.d0
      if(m.gt.(a(88)+0.1d0))then
         rgammf = 0.d0
      else
         b1 = MAX(0.d0,a(83) + a(84)*(m1-a(85))**a(86))
         if(m.le.m1)then
            rgammf = a(83) + a(84)*ABS(m-a(85))**a(86)
         elseif(m.le.a(88))then
            rgammf = b1 + (a(89) - b1)*((m - m1)/(a(88) - m1))**a(87)
         else
            if(a(88).gt.m1) b1 = a(89)
            rgammf = b1 - 10.d0*b1*(m - a(88))
         endif
         rgammf = MAX(rgammf,0.d0)
      endif
*
      return
      end
***
      real*8 FUNCTION rhookf(m,mhook)
      implicit none
      real*8 m,mhook,m2,b2,a(200)
      common /MSCFF/ a
*
* A function to evalute the radius at the start of
* the MS hook ( for those stars that have one ).
* Note that this function is only valid for M > Mhook.
* (JH 24/11/97)
*
      if(m.le.mhook)then
         rhookf = 0.d0
      elseif(m.le.a(94))then
         rhookf = a(95)*SQRT((m-mhook)/(a(94)-mhook))
      elseif(m.le.2.d0)then
         m2 = 2.d0
         b2 = (a(90) + a(91)*m2**(7.d0/2.d0))/
     &        (a(92)*m2**3 + m2**a(93)) - 1.d0
         rhookf = a(95) + (b2-a(95))*((m-a(94))/(m2-a(94)))**a(96)
      else
         rhookf = (a(90) + a(91)*m**(7.d0/2.d0))/
     &            (a(92)*m**3 + m**a(93)) - 1.d0
      endif
*
      return
      end
***
      real*8 FUNCTION lbgbf(m)
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate the luminosity at the end of the
* FGB ( for those models that have one )
* Note that this function is only valid for LM & IM stars
* (JH 24/11/97)
*
      lbgbf = (a(1)*m**a(5) + a(2)*m**a(8))/
     &        (a(3) + a(4)*m**a(7) + m**a(6))
*
      return
      end
***
      real*8 FUNCTION lbgbdf(m)
      real*8 m,a(200)
      real*8 f,df,g,dg
      common /GBCFF/ a
*
* A function to evaluate the derivitive of the Lbgb function.
* Note that this function is only valid for LM & IM stars
* (JH 24/11/97)
*
      f = a(1)*m**a(5) + a(2)*m**a(8)
      df = a(5)*a(1)*m**(a(5)-1.d0) + a(8)*a(2)*m**(a(8)-1.d0)
      g = a(3) + a(4)*m**a(7) + m**a(6)
      dg = a(7)*a(4)*m**(a(7)-1.d0) + a(6)*m**(a(6)-1.d0)
*
      lbgbdf = (df*g - f*dg)/(g*g)
*
      return
      end
***
      real*8 FUNCTION lbagbf(m,mhefl)
      implicit none
      real*8 m,mhefl,a4,a(200)
      common /GBCFF/ a
*
* A function to evaluate the BAGB luminosity. (OP 21/04/98)
* Continuity between LM and IM functions is ensured by setting
* gbp(16) = lbagbf(mhefl,0.0) with gbp(16) = 1.0.
*
      a4 = (a(9)*mhefl**a(10) - a(16))/(exp(mhefl*a(11))*a(16))
*
      if(m.lt.mhefl)then
         lbagbf = a(9)*m**a(10)/(1.d0 + a4*exp(m*a(11)))
      else
         lbagbf = (a(12) + a(13)*m**(a(15)+1.8d0))/(a(14) + m**a(15))
      endif
*
      return
      end
***
      real*8 FUNCTION rgbf(m,lum)
      implicit none
      real*8 m,lum,a1,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius on the GB.
* (JH 24/11/97)
*
      a1 = MIN(a(20)/m**a(21),a(22)/m**a(23))
      rgbf = a1*(lum**a(18) + a(17)*lum**a(19))
*
      return
      end
***
      real*8 FUNCTION rgbdf(m,lum)
      implicit none
      real*8 m,lum,a1,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius derivitive on the GB (as f(L)).
* (JH 24/11/97)
*
      a1 = MIN(a(20)/m**a(21),a(22)/m**a(23))
      rgbdf = a1*(a(18)*lum**(a(18)-1.d0) +
     &            a(17)*a(19)*lum**(a(19)-1.d0))
*
      return
      end
***
      real*8 FUNCTION ragbf(m,lum,mhelf)
      implicit none
      real*8 m,lum,mhelf,m1,a1,a4,xx,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius on the AGB.
* (JH 24/11/97)
*
      m1 = mhelf - 0.2d0
      if(m.ge.mhelf)then
         xx = a(24)
      elseif(m.ge.m1)then
         xx = 1.d0 + 5.d0*(a(24)-1.d0)*(m-m1)
      else
         xx = 1.d0
      endif
      a4 = xx*a(19)
      if(m.le.m1)then
         a1 = a(29) + a(30)*m
      elseif(m.ge.mhelf)then
         a1 = MIN(a(25)/m**a(26),a(27)/m**a(28))
      else
         a1 = a(31) + 5.d0*(a(32)-a(31))*(m-m1)
      endif
*
      ragbf = a1*(lum**a(18) + a(17)*lum**a4)
*
      return
      end
***
      real*8 FUNCTION ragbdf(m,lum,mhelf)
      implicit none
      real*8 m,lum,mhelf,m1,a1,a4,xx,a(200)
      common /GBCFF/ a
*
* A function to evaluate radius derivitive on the AGB (as f(L)).
* (JH 24/11/97)
*
      m1 = mhelf - 0.2d0
      if(m.ge.mhelf)then
         xx = a(24)
      elseif(m.ge.m1)then
         xx = 1.d0 + 5.d0*(a(24)-1.d0)*(m-m1)
      else
         xx = 1.d0
      endif
      a4 = xx*a(19)
      if(m.le.m1)then
         a1 = a(29) + a(30)*m
      elseif(m.ge.mhelf)then
         a1 = MIN(a(25)/m**a(26),a(27)/m**a(28))
      else
         a1 = a(31) + 5.d0*(a(32)-a(31))*(m-m1)
      endif
*
      ragbdf = a1*(a(18)*lum**(a(18)-1.d0) +
     &             a(17)*a4*lum**(a4-1.d0))
*
      return
      end
***
      real*8 FUNCTION mctmsf(m)
      implicit none
      real*8 m,m525
*
* A function to evaluate core mass at the end of the MS as a
* fraction of the BGB value, i.e. this must be multiplied by
* the BGB value (see below) to give the actual core mass (JH 5/9/99)
*
      m525 = m**(21.d0/4.d0)
      mctmsf = (1.586d0 + m525)/(2.434d0 + 1.02d0*m525)
*
      return
      end
***
      real*8 FUNCTION mcheif(m,mhefl,mchefl)
      implicit none
      real*8 m,mhefl,mchefl,mcbagb,a3,a(200)
      common /GBCFF/ a
      real*8 mcagbf
      external mcagbf
*
* A function to evaluate core mass at BGB or He ignition
* (depending on mchefl) for IM & HM stars  (OP 25/11/97)
*
      mcbagb = mcagbf(m)
      a3 = mchefl**4 - a(33)*mhefl**a(34)
      mcheif = MIN(0.95d0*mcbagb,(a3 + a(33)*m**a(34))**(1.d0/4.d0))
*
      return
      end
***
      real*8 FUNCTION mheif(mc,mhefl,mchefl)
      implicit none
      real*8 mc,mhefl,mchefl,m1,m2,a3,a(200)
      common /GBCFF/ a
      real*8 mbagbf
      external mbagbf
*
* A function to evaluate mass at BGB or He ignition
* (depending on mchefl) for IM & HM stars by inverting
* mcheif
*
      m1 = mbagbf(mc/0.95d0)
      a3 = mchefl**4 - a(33)*mhefl**a(34)
      m2 = ((mc**4 - a3)/a(33))**(1.d0/a(34))
      mheif = MAX(m1,m2)
*
      return
      end
***
      real*8 FUNCTION mcagbf(m)
      implicit none
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate core mass at the BAGB  (OP 25/11/97)
*
      mcagbf = (a(37) + a(35)*m**a(36))**(1.d0/4.d0)
*
      return
      end
***
      real*8 FUNCTION mbagbf(mc)
      implicit none
      real*8 mc,mc4,a(200)
      common /GBCFF/ a
*
* A function to evaluate mass at the BAGB by inverting mcagbf.
*
      mc4 = mc**4
      if(mc4.gt.a(37))then
         mbagbf = ((mc4 - a(37))/a(35))**(1.d0/a(36))
      else
         mbagbf = 0.d0
      endif
*
      return
      end
***
      real*8 FUNCTION mcgbtf(t,A,GB,tinf1,tinf2,tx)
      implicit none
      real*8 t,A,GB(10),tinf1,tinf2,tx
*
* A function to evaluate Mc given t for GB, AGB and NHe stars
*
      if(t.le.tx)then
         mcgbtf = ((GB(5)-1.d0)*A*GB(4)*(tinf1 - t))**
     &                              (1.d0/(1.d0-GB(5)))
      else
         mcgbtf = ((GB(6)-1.d0)*A*GB(3)*(tinf2 - t))**
     &                              (1.d0/(1.d0-GB(6)))
      endif
*
      return
      end
***
      real*8 FUNCTION lgbtf(t,A,GB,tinf1,tinf2,tx)
      implicit none
      real*8 t,A,GB(10),tinf1,tinf2,tx
*
* A function to evaluate L given t for GB, AGB and NHe stars
*
      if(t.le.tx)then
         lgbtf = GB(4)*(((GB(5)-1.d0)*A*GB(4)*(tinf1 - t))**
     &                              (GB(5)/(1.d0-GB(5))))
      else
         lgbtf = GB(3)*(((GB(6)-1.d0)*A*GB(3)*(tinf2 - t))**
     &                              (GB(6)/(1.d0-GB(6))))
      endif
*
      return
      end
***
      real*8 FUNCTION mcgbf(lum,GB,lx)
      implicit none
      real*8 lum,GB(10),lx
*
* A function to evaluate Mc given L for GB, AGB and NHe stars
*
      if(lum.le.lx)then
         mcgbf = (lum/GB(4))**(1.d0/GB(5))
      else
         mcgbf = (lum/GB(3))**(1.d0/GB(6))
      endif
*
      return
      end
***
      real*8 FUNCTION lmcgbf(mc,GB)
      implicit none
      real*8 mc,GB(10)
*
* A function to evaluate L given Mc for GB, AGB and NHe stars
*
      if(mc.le.GB(7))then
         lmcgbf = GB(4)*(mc**GB(5))
      else
         lmcgbf = GB(3)*(mc**GB(6))
      endif
*
      return
      end
***
      real*8 FUNCTION lHeIf(m,mhefl)
      implicit none
      real*8 m,mhefl,a(200)
      common /GBCFF/ a
*
* A function to evaluate He-ignition luminosity  (OP 24/11/97)
* Continuity between the LM and IM functions is ensured with a first
* call setting lhefl = lHeIf(mhefl,0.0)
*
      if(m.lt.mhefl)then
         lHeIf = a(38)*m**a(39)/(1.d0 + a(41)*EXP(m*a(40)))
      else
         lHeIf = (a(42) + a(43)*m**3.8d0)/(a(44) + m**2)
      endif
*
      return
      end
***
      real*8 FUNCTION lHef(m)
      implicit none
      real*8 m,a(200)
      common /GBCFF/ a
*
* A function to evaluate the ratio LHe,min/LHeI  (OP 20/11/97)
* Note that this function is everywhere <= 1, and is only valid
* for IM stars
*
      lHef = (a(45) + a(46)*m**(a(48)+0.1d0))/(a(47) + m**a(48))
*
      return
      end
***
      real*8 FUNCTION rminf(m)
      implicit none
      real*8 m,mx,a(200)
      common /GBCFF/ a
*
* A function to evaluate the minimum radius during He-burning
* for IM & HM stars  (OP 20/11/97)
*
      mx = m**a(53)
      rminf = (a(49)*m + (a(50)*m)**a(52)*mx)/(a(51) + mx)
*
      return
      end
***
      real*8 FUNCTION tHef(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,mm,a(200)
      common /GBCFF/ a
      real*8 themsf
      external themsf
*
* A function to evaluate the He-burning lifetime.  (OP 26/11/97)
* For IM & HM stars, tHef is relative to tBGB.
* Continuity between LM and IM stars is ensured by setting
* thefl = tHef(mhefl,0.0,,0.0), and the call to themsf ensures
* continuity between HB and NHe stars as Menv -> 0.
*
      if(m.le.mhefl)then
         mm = MAX((mhefl - m)/(mhefl - mc),1.0d-12)
         tHef = (a(54) + (themsf(mc) - a(54))*mm**a(55))*
     &          (1.d0 + a(57)*EXP(m*a(56)))
      else
         mm = m**5
         tHef = (a(58)*m**a(61) + a(59)*mm)/(a(60) + mm)
      endif
*
      return
      end
***
      real*8 FUNCTION tblf(m,mhefl,mfgb)
      implicit none
      real*8 m,mhefl,mfgb,mr,m1,m2,r1,a(200)
      common /GBCFF/ a
      real*8 lheif,rminf,ragbf
      external lheif,rminf,ragbf
*
* A function to evaluate the blue-loop fraction of the He-burning
* lifetime for IM & HM stars  (OP 28/01/98)
*
      mr = mhefl/mfgb
      if(m.le.mfgb) then
         m1 = m/mfgb
         m2 = log10(m1)/log10(mr)
         m2 = max(m2,1.0d-12)
         tblf = a(64)*m1**a(63) + a(65)*m2**a(62)
      else
         r1 = 1.d0 - rminf(m)/ragbf(m,lheif(m,mhefl),mhefl)
         r1 = max(r1,1.0d-12)
         tblf = a(66)*m**a(67)*r1**a(68)
      end if
      tblf = MIN(1.d0,MAX(0.d0,tblf))
      if(tblf.lt.1.0d-10) tblf = 0.d0
*
      return
      end
***
      real*8 FUNCTION lzahbf(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,mm,a4,a5,a(200)
      common /GBCFF/ a
      real*8 lzhef
      external lzhef
*
* A function to evaluate the ZAHB luminosity for LM stars. (OP 28/01/98)
* Continuity with LHe,min for IM stars is ensured by setting
* lx = lHeif(mhefl,z,0.0,1.0)*lHef(mhefl,z,mfgb), and the call to lzhef
* ensures continuity between the ZAHB and the NHe-ZAMS as Menv -> 0.
*
      a5 = lzhef(mc)
      a4 = (a(69) + a5 - a(74))/((a(74) - a5)*exp(a(71)*mhefl))
      mm = MAX((m-mc)/(mhefl - mc),1.0d-12)
      lzahbf = a5 + (1.d0 + a(72))*a(69)*mm**a(70)/
     &         ((1.d0 + a(72)*mm**a(73))*(1.d0 + a4*EXP(m*a(71))))
*
      return
      end
***
      real*8 FUNCTION rzahbf(m,mc,mhefl)
      implicit none
      real*8 m,mc,mhefl,rx,ry,mm,f,a(200)
      common /GBCFF/ a
      real*8 rzhef,rgbf,lzahbf
*
* A function to evaluate the ZAHB radius for LM stars. (OP 28/01/98)
* Continuity with R(LHe,min) for IM stars is ensured by setting
* lx = lHeif(mhefl,z,0.0,1.0)*lHef(mhefl,z,mfgb), and the call to rzhef
* ensures continuity between the ZAHB and the NHe-ZAMS as Menv -> 0.
*
      rx = rzhef(mc)
      ry = rgbf(m,lzahbf(m,mc,mhefl))
      mm = MAX((m-mc)/(mhefl - mc),1.0d-12)
      f = (1.d0 + a(76))*mm**a(75)/(1.d0 + a(76)*mm**a(77))
      rzahbf = (1.d0 - f)*rx + f*ry
*
      return
      end
***
      real*8 FUNCTION lzhef(m)
      implicit none
      real*8 m,m15
*
* A function to evaluate Naked Helium star 'ZAMS' luminosity
*
      m15 = m*SQRT(m)
      lzhef = 1.5262d+04*m**(41.d0/4.d0)/
     &        (0.0469d0 + m**6*(31.18d0 + m15*(29.54d0 + m15)))
*
      return
      end
***
      real*8 FUNCTION rzhef(m)
      implicit none
      real*8 m
*
* A function to evaluate Helium star 'ZAMS' radius
*
      rzhef = 0.2391d0*m**4.6d0/(0.0065d0 + (0.162d0 + m)*m**3)
*
      return
      end
***
      real*8 FUNCTION themsf(m)
      implicit none
      real*8 m
*
* A function to evaluate Helium star main sequence lifetime
*
      themsf = (0.4129d0 + 18.81d0*m**4 + 1.853d0*m**6)/m**(13.d0/2.d0)
*
      return
      end
***
      real*8 FUNCTION rhehgf(m,lum,rx,lx)
      implicit none
      real*8 m,lum,rx,lx,cm
*
* A function to evaluate Helium star radius on the Hertzsprung gap
* from its mass and luminosity.
*
      cm = 2.0d-03*m**(5.d0/2.d0)/(2.d0 + m**5)
      rhehgf = rx*(lum/lx)**0.2d0 + 0.02d0*(EXP(cm*lum) - EXP(cm*lx))
*
      return
      end
***
      real*8 FUNCTION rhegbf(lum)
      implicit none
      real*8 lum
*
* A function to evaluate Helium star radius on the giant branch.
*
      rhegbf = 0.08d0*lum**(3.d0/4.d0)
*
      return
      end
***
      real*8 FUNCTION lpertf(m,mew)
      implicit none
      real*8 m,mew
      real*8 b,c
*
* A function to obtain the exponent that perturbs luminosity.
*
      b = 0.002d0*MAX(1.d0,2.5d0/m)
      c = 3.d0
      lpertf = ((1.d0 + b**c)*((mew/b)**c))/(1.d0+(mew/b)**c)
*
      return
      end
***
      real*8 FUNCTION rpertf(m,mew,r,rc)
      implicit none
      real*8 m,mew,r,rc
      real*8 a,b,c,q,fac,facmax
*
* A function to obtain the exponent that perturbs radius.
*
      if(mew.le.0.d0)then
         rpertf = 0.d0
      else
         a = 0.1d0
         b = 0.006d0*MAX(1.d0,2.5d0/m)
         c = 3.d0
         q = log(r/rc)
         fac = a/q
         facmax = -14.d0/log10(mew)
         fac = MIN(fac,facmax)
         rpertf = ((1.d0 + b**c)*((mew/b)**c)*(mew**fac))/
     &            (1.d0+(mew/b)**c)
      endif
*
      return
      end
***
      real*8 FUNCTION vrotf(m,ST)
      implicit none
      integer ST
      real*8 m
*
      if(ST.gt.0)then
         if(m.gt.6.35d0)then
            vrotf = (10.d0*m**(-0.0354d0))/(0.0389d0+m**(-7.95d0))
         else
            vrotf = (13.4d0*m**(-0.12d0))/(0.0389d0+m**(-7.95d0))
         endif
      else
         vrotf = 330.d0*m**3.3d0/(15.d0 + m**3.45d0)
      endif
*
      return
      end
***
      real*8 FUNCTION celamf(kw,m,lum,rad,rzams,menv,fac)
      implicit none
      integer kw
      real*8 m,lum,rad,rzams,menv,fac
      real*8 lam1,lam2,m1,logm,logl
      real*8 aa,bb,cc,dd
*
* A function to estimate lambda for common-envelope.
*
      if(fac.ge.0.d0)then
*
* No fits yet for naked He stars...
*
         if(kw.gt.6)then
            celamf = 0.5d0
            goto 90
         endif
*
         if(menv.gt.0.d0)then
* Formulae for giant-like stars; also used for HG and CHeB stars close
* to the Hayashi track.
            logl = log10(lum)
            logm = log10(m)
            if(kw.le.5)then
               m1 = m
               if(kw.gt.3) m1 = 100.d0
               lam1 = 3.d0/(2.4d0 + 1.d0/m1**(3.d0/2.d0)) - 0.15d0*logl
               lam1 = MIN(lam1,0.8d0)
            else
               lam1 = -3.5d0 - 0.75d0*logm + logl
            endif
            if(kw.gt.3)then
               lam2 = MIN(0.9d0,0.58d0 + 0.75d0*logm) - 0.08d0*logl
               if(kw.lt.6)then
                  lam1 = MIN(lam2,lam1)
               else
                  lam1 = MAX(lam2,lam1)
                  lam1 = MIN(lam1,1.d0)
               endif
            endif
            lam1 = 2.d0*lam1
            if(fac.gt.0.d0)then
* Use a fraction FAC of the ionization energy in the energy balance.
               if(kw.le.3)then
                  aa = MIN(1.2d0*(logm - 0.25d0)**2 - 0.7d0,-0.5d0)
               else
                  aa = MAX(-0.2d0 - logm,-0.5d0)
               endif
               bb = MAX(3.d0 - 5.d0*logm,1.5d0)
               cc = MAX(3.7d0 + 1.6d0*logm,3.3d0 + 2.1d0*logm)
               lam2 = aa + ATAN(bb*(cc - logl))
               if(kw.le.3)then
                  dd = MAX(0.d0,MIN(0.15d0,0.15d0 - 0.25d0*logm))
                  lam2 = lam2 + dd*(logl - 2.d0)
               endif
               lam2 = MAX(lam2,1.d-2)
               lam2 = MAX(1.d0/lam2,lam1)
               if(fac.ge.1.d0)then
                  lam1 = lam2
               else
                  lam1 = lam1 + fac*(lam2 - lam1)
               endif
            endif
         endif
*
         if(menv.lt.1.d0)then
* Formula for HG stars; also reasonable for CHeB stars in blue loop.
            lam2 = 0.42d0*(rzams/rad)**0.4d0
* Alternatively:
*           lam2 = 0.3d0*(rtms/rad)**0.4d0
            lam2 = 2.d0*lam2
         endif
*
         if(menv.le.0.d0)then
            celamf = lam2
         elseif(menv.ge.1.d0)then
            celamf = lam1
         else
* Interpolate between HG and GB values depending on conv. envelope mass.
            celamf = lam2 + sqrt(menv)*(lam1 - lam2)
         endif
*
 90      continue
*
      endif
*
      return
      end
***



***
      real*8 FUNCTION celamf_nanjing(kw,mzams,menv,zpars,rad,rzams,lam)
      implicit none
      integer kw,ceidx,i,stage_num
      real*8 mzams,rad,rzams,menv,lam,zpars(20)
      real*8 met,x,lambda_b,lambda_g,popI_popII_Z
      integer stage(44)
      real*8 masses(44)
      real*8 a_popI_lambdab(44),a_popI_lambdag(44)
      real*8 a_popII_lambdab(44),a_popII_lambdag(44)
      double precision b1_popI_lambdab(44),b1_popI_lambdag(44)
      double precision b1_popII_lambdab(44),b1_popII_lambdag(44)
      double precision b2_popI_lambdab(44),b2_popI_lambdag(44)
      double precision b2_popII_lambdab(44),b2_popII_lambdag(44)
      double precision b3_popI_lambdab(44),b3_popI_lambdag(44)
      double precision b3_popII_lambdab(44),b3_popII_lambdag(44)
      double precision b4_popI_lambdab(44),b4_popI_lambdag(44)
      double precision b4_popII_lambdab(44),b4_popII_lambdag(44)
      double precision b5_popI_lambdab(44),b5_popI_lambdag(44)
      double precision b5_popII_lambdab(44),b5_popII_lambdag(44)
*
* A function that uses the fitting formulae from Xu & Li 2010
* to determine lambda for common-envelope.
* Written by Michael Zevin, 2019
*

* Values for all fitting parameters
      stage=(/ 1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     &      2,2,2,2,2,2,2,2,2,2,2,2,2,2,
     &      3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3 /)
      masses=(/ 1,2,3,4,5,6,7,8,9,10,12,14,16,20,
     &      1,2,3,4,5,6,7,8,9,10,12,14,16,20,
     &      1,2,3,4,5,6,7,8,9,10,12,14,16,20,50,100 /)
* Population I, lambda_g
      a_popI_lambdag=(/ 17.583280,1.076580,1.307050,1.021830,0.857230,
     &     0.784280,0.760090,0.738260,0.715710,0.692450,
     &     0.685440,0.782150,0.853840,0.653340,63.612590,
     &     13.660580,-6.738420,-3.804550,-5.332790,7.687680,
     &     2.189520,-2.243540,-0.729900,0.367420,0.461560,
     &     0.568000,0.360140,0.226930,394.595560,0.482710,
     &     0.448890,0.131530,-0.004560,0.230830,0.262940,
     &     0.269560,0.405870,0.255490,-15.116720,-39.930890,
     &     -65.396020,-480055.679910,0.159000,0.246600 /)
      b1_popI_lambdag=(/ -3.484355d+01,-1.041000d-02,
     &     -2.292400d-01,-1.024000d-01,
     &     -4.922000d-02,-2.959000d-02,
     &     -2.412000d-02,-1.995000d-02,
     &     -1.657000d-02,-1.398000d-02,
     &     -1.394000d-02,-2.326000d-02,
     &     -3.086000d-02,-1.346000d-02,
     &     -3.998949d+02,-2.480310d+00,
     &     1.066560d+00,2.930800d-01,
     &     2.272800d-01,-3.072300d-01,
     &     -6.892000d-02,1.091800d-01,
     &     3.910000d-02,-3.440000d-03,
     &     -6.600000d-03,-4.700000d-03,
     &     -2.540000d-03,-8.767800d-04,
     &     -3.781684d+03,5.840000d-03,
     &     1.102000d-02,9.840000d-03,
     &     4.260000d-03,-2.660000d-03,
     &     -2.530000d-03,-2.190000d-03,
     &     -5.100000d-03,-1.520000d-03,
     &     6.331000d-02,1.366700d-01,
     &     1.976300d-01,7.874840d+06,
     &     -3.944510d-04,-1.210000d-03 /)
      b2_popI_lambdag=(/ 1.070536d+01,-4.905530d-05,
     &     1.847000d-02,4.930000d-03,
     &     1.370000d-03,5.201300d-04,
     &     3.471040d-04,2.378420d-04,
     &     1.646070d-04,1.172560d-04,
     &     1.208450d-04,3.259840d-04,
     &     5.508780d-04,1.169530d-04,
     &     9.596205d+02,1.527500d-01,
     &     -5.344000d-02,-6.030000d-03,
     &     -2.850000d-03,4.450000d-03,
     &     8.009360d-04,-1.790000d-03,
     &     -5.781320d-04,1.278380d-05,
     &     3.962500d-05,1.578180d-05,
     &     7.496390d-06,1.288520d-06,
     &     1.206513d+04,-6.220510d-05,
     &     -6.466290d-05,-2.898320d-05,
     &     4.711170d-06,2.217880d-05,
     &     1.322720d-05,7.977430d-06,
     &     2.738660d-05,3.352390d-06,
     &     -8.815420d-05,-1.559580d-04,
     &     -1.990780d-04,-4.305460d+07,
     &     2.884520d-07,1.890290d-06 /)
      b3_popI_lambdag=(/ 8.490420d+00,1.135280d-06,
     &     -5.062160d-04,-8.163430d-05,
     &     -1.361630d-05,-3.451720d-06,
     &     -1.923470d-06,-1.098030d-06,
     &     -6.319350d-07,-3.814870d-07,
     &     -4.290710d-07,-1.949910d-06,
     &     -4.376710d-06,-4.213260d-07,
     &     -7.952070d+02,-3.030000d-03,
     &     1.160000d-03,4.004710d-05,
     &     1.164080d-05,-2.704490d-05,
     &     -3.780920d-06,1.332440d-05,
     &     3.707200d-06,-1.072200d-08,
     &     -9.986670d-08,-2.212070d-08,
     &     -9.201030d-09,-6.129120d-10,
     &     -1.239746d+04,2.415310d-07,
     &     5.668570d-09,2.635190d-08,
     &     -1.728580d-08,-2.356960d-08,
     &     -7.122050d-09,-1.532960d-09,
     &     -5.744760d-08,2.242240d-10,
     &     4.098200d-08,5.940760d-08,
     &     6.687660d-08,7.846990d+07,
     &     -6.351320d-11,-1.120660d-09 /)
      b4_popI_lambdag=(/ 0.000000d+00,-3.916090d-09,
     &     4.570980d-06,4.554260d-07,
     &     4.686830d-08,8.172480d-09,
     &     3.796090d-09,1.790440d-09,
     &     8.520820d-10,4.358180d-10,
     &     5.291690d-10,4.080440d-09,
     &     1.250750d-08,5.244250d-10,
     &     0.000000d+00,0.000000d+00,
     &     -9.344460d-06,0.000000d+00,
     &     0.000000d+00,5.897120d-08,
     &     6.348200d-09,-4.578290d-08,
     &     -1.070360d-08,0.000000d+00,
     &     -8.841340d-11,1.084720d-11,
     &     3.938280d-12,0.000000d+00,
     &     0.000000d+00,-3.187200d-10,
     &     7.218180d-10,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     4.902180d-11,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,2.225800d-13 /)
      b5_popI_lambdag=(/ 0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,5.903130d-11,
     &     1.148330d-11,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     -1.220100d-12,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00 /)
* Population I, lambda_b
      a_popI_lambdab=(/ 8.358970,2.053630,2.408310,1.818600,1.525810,
     &     1.416010,1.383440,1.355160,1.325490,1.293120,
     &     1.285930,1.393320,1.431770,1.225470,46.009780,
     &     34.418260,-42.985130,-7.309800,-9.936470,13.914650,
     &     4.123870,-3.891890,0.863690,0.742330,-11.995370,
     &     1.128890,0.841430,0.487240,245.382300,0.889540,
     &     -0.046690,-0.373220,-0.800110,-2.771400,-0.632660,
     &     -0.128800,1.198040,0.370700,-58.037320,-106.905530,
     &     -154.705590,-260484.857240,0.313210,0.376000 /)
      b1_popI_lambdab=(/ -1.889048d+01,-6.850000d-03,
     &     -4.245900d-01,-1.746400d-01,
     &     -8.125000d-02,-4.965000d-02,
     &     -4.093000d-02,-3.414000d-02,
     &     -2.845000d-02,-2.371000d-02,
     &     -2.209000d-02,-3.180000d-02,
     &     -3.533000d-02,-1.957000d-02,
     &     -2.986499d+02,-6.652590d+00,
     &     7.901340d+00,5.664700d-01,
     &     4.283100d-01,-5.557900d-01,
     &     -1.297900d-01,1.937800d-01,
     &     -9.950000d-03,-6.230000d-03,
     &     9.920000d-02,-9.010000d-03,
     &     -5.760000d-03,-1.770000d-03,
     &     -2.357511d+03,9.800000d-03,
     &     7.640000d-03,9.430000d-03,
     &     9.920000d-03,6.467000d-02,
     &     2.054000d-02,9.900000d-03,
     &     -1.961000d-02,2.672210d-04,
     &     2.363300d-01,3.646900d-01,
     &     4.671800d-01,4.267590d+06,
     &     -7.503840d-04,-1.800000d-03 /)
      b2_popI_lambdab=(/ 1.047651d+01,-3.427390d-04,
     &     3.431000d-02,8.280000d-03,
     &     2.190000d-03,8.515270d-04,
     &     5.789520d-04,4.020650d-04,
     &     2.790970d-04,1.937640d-04,
     &     1.797640d-04,3.959170d-04,
     &     5.111280d-04,1.503320d-04,
     &     7.274094d+02,4.382300d-01,
     &     -5.464600d-01,-1.176000d-02,
     &     -5.440000d-03,8.090000d-03,
     &     1.530000d-03,-3.200000d-03,
     &     4.808370d-05,2.041970d-05,
     &     -2.898100d-04,3.040770d-05,
     &     1.688540d-05,2.602540d-06,
     &     7.542037d+03,-3.141100d-05,
     &     -4.327260d-05,-3.260330d-05,
     &     -3.032470d-05,-4.015370d-04,
     &     -1.364600d-04,-6.714550d-05,
     &     1.282220d-04,-9.864640d-06,
     &     -3.205350d-04,-4.147200d-04,
     &     -4.701690d-04,-2.330160d+07,
     &     5.385450d-07,2.810830d-06 /)
      b3_popI_lambdab=(/ 9.935200d-01,3.939870d-06,
     &     -9.268790d-04,-1.317270d-04,
     &     -2.052700d-05,-5.543840d-06,
     &     -3.192270d-06,-1.859310d-06,
     &     -1.072540d-06,-6.195760d-07,
     &     -6.215560d-07,-2.231320d-06,
     &     -3.576330d-06,-5.071290d-07,
     &     -6.076680d+02,-9.530000d-03,
     &     1.863000d-02,7.901120d-05,
     &     2.258480d-05,-4.948720d-05,
     &     -7.432270d-06,2.395040d-05,
     &     -6.104540d-08,-1.303880d-08,
     &     3.627510d-07,-4.319640d-08,
     &     -2.082700d-08,-1.258240d-09,
     &     -7.768917d+03,7.669790d-08,
     &     9.319420d-08,5.378230d-08,
     &     5.262350d-08,7.984660d-07,
     &     2.866100d-07,1.335680d-07,
     &     -3.412780d-07,2.261850d-08,
     &     1.451290d-07,1.573490d-07,
     &     1.577730d-07,4.241020d+07,
     &     -1.169460d-10,-1.673860d-09 /)
      b4_popI_lambdab=(/ 0.000000d+00,-1.182370d-08,
     &     8.245220d-06,7.083290d-07,
     &     6.791690d-08,1.323360d-08,
     &     6.409020d-09,3.088320d-09,
     &     1.468010d-09,7.042270d-10,
     &     7.594440d-10,4.508310d-09,
     &     9.367780d-09,6.068000d-10,
     &     0.000000d+00,0.000000d+00,
     &     -3.131010d-04,0.000000d+00,
     &     0.000000d+00,1.088990d-07,
     &     1.294180d-08,-8.289590d-08,
     &     -2.795040d-12,0.000000d+00,
     &     -1.655850d-10,2.145450d-11,
     &     8.978130d-12,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     3.356140d-10,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,3.350560d-13 /)
      b5_popI_lambdab=(/ 0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     2.074680d-06,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,1.078430d-10,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00 /)
* Population II, lambda_g
      a_popII_lambdag=(/ 0.571920,1.418960,0.992180,0.921720,0.876470,
     &     0.836360,0.816880,0.793960,0.805000,1.027000,
     &     1.179340,1.195260,1.177310,1.074960,0.248160,
     &     -56.034780,-6.474760,-3.212990,-0.385610,0.145760,
     &     0.218730,0.247480,0.207960,0.353550,0.371880,
     &     0.361630,0.341960,0.266910,0.155040,0.305940,
     &     0.054340,0.367790,0.312520,0.268020,0.268830,
     &     0.252490,0.319510,0.814910,-13.440000,-44.196400,
     &     -118.137570,-144.534560,0.492870,0.817160 /)
      b1_popII_lambdag=(/ 3.937000d-02,-4.266000d-01,
     &     -1.008200d-01,-6.187000d-02,
     &     -4.103000d-02,-2.806000d-02,
     &     -2.324000d-02,-1.903000d-02,
     &     -2.000000d-02,-5.685000d-02,
     &     -8.481000d-02,-8.503000d-02,
     &     -7.834000d-02,-5.737000d-02,
     &     -4.102000d-02,1.367490d+01,
     &     8.328000d-01,2.958300d-01,
     &     4.270000d-02,5.620000d-03,
     &     1.540000d-03,-9.933800d-05,
     &     6.629210d-04,-3.880000d-03,
     &     -3.650000d-03,-3.280000d-03,
     &     -2.800000d-03,-1.610000d-03,
     &     -1.238000d-02,-9.588580d-04,
     &     3.900000d-03,-9.910000d-03,
     &     -5.270000d-03,-2.480000d-03,
     &     -2.190000d-03,-1.610000d-03,
     &     -3.920000d-03,-1.610000d-03,
     &     8.141000d-02,2.259200d-01,
     &     5.273700d-01,4.657900d-01,
     &     -4.390000d-03,-1.436000d-02 /)
      b2_popII_lambdag=(/ -2.210000d-03,5.792000d-02,
     &     4.510000d-03,1.770000d-03,
     &     7.914440d-04,3.733460d-04,
     &     2.580400d-04,1.775740d-04,
     &     2.018720d-04,1.640000d-03,
     &     3.290000d-03,3.240000d-03,
     &     2.750000d-03,1.530000d-03,
     &     2.800000d-03,-1.095330d+00,
     &     -3.412000d-02,-8.330000d-03,
     &     -9.694800d-04,-1.302730d-04,
     &     -5.188060d-05,-1.992720d-05,
     &     -1.846630d-05,1.565730d-05,
     &     1.249440d-05,1.031190d-05,
     &     7.828650d-06,3.337800d-06,
     &     3.966330d-04,1.121740d-04,
     &     9.446090d-06,1.194110d-04,
     &     3.603480d-05,6.452290d-06,
     &     4.129410d-06,8.354780d-07,
     &     2.318150d-05,-8.133520d-06,
     &     -1.641000d-04,-3.851240d-04,
     &     -7.847900d-04,-4.991970d-04,
     &     1.067660d-05,9.311430d-05 /)
      b3_popII_lambdag=(/ 3.491920d-05,-2.810000d-03,
     &     -5.536320d-05,-1.426770d-05,
     &     -4.416440d-06,-1.470160d-06,
     &     -8.546960d-07,-5.042620d-07,
     &     -6.429500d-07,-1.652060d-05,
     &     -4.690960d-05,-4.589190d-05,
     &     -3.581080d-05,-1.490050d-05,
     &     -6.204190d-05,2.925000d-02,
     &     4.583990d-04,7.556460d-05,
     &     6.644550d-06,7.064590d-07,
     &     2.602830d-07,9.475040d-08,
     &     6.589830d-08,-1.981730d-08,
     &     -1.323880d-08,-9.927120d-09,
     &     -6.666840d-09,-2.155500d-09,
     &     -5.332900d-06,-1.040790d-06,
     &     -3.872780d-08,-3.595740d-07,
     &     -3.224450d-08,1.696090d-08,
     &     1.331380d-08,1.259990d-08,
     &     -6.594180d-08,1.957750d-08,
     &     1.106000d-07,2.193240d-07,
     &     3.895850d-07,1.780270d-07,
     &     -9.220150d-09,-2.653900d-07 /)
      b4_popII_lambdag=(/ -1.730460d-07,4.610000d-05,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     3.160520d-08,3.456400d-09,
     &     0.000000d+00,3.339570d-10,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     7.995750d-11,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     2.589260d-12,3.307730d-10 /)
      b5_popII_lambdag=(/ 0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     -6.842880d-11,-3.915360d-12,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,-1.512070d-13 /)
* Population II, lambda_b
      a_popII_lambdab=(/ 1.060310,2.561080,1.781400,1.659140,1.587010,
     &     1.527000,1.499950,1.468260,1.491960,1.739660,
     &     1.636340,1.455730,1.333780,1.271380,0.372940,
     &     -103.925380,-12.408320,-5.892530,-0.671760,0.309410,
     &     0.448620,0.502210,0.393420,0.757460,0.852490,
     &     0.852710,0.832540,0.697460,0.240120,0.545200,
     &     -0.475000,-0.210600,-0.120270,0.265780,0.815800,
     &     0.749240,0.731470,-9.265190,-51.152520,-140.000000,
     &     -358.400000,-436.007770,0.821000,1.253320 /)
      b1_popII_lambdab=(/ 8.173000d-02,-7.556200d-01,
     &     -1.713800d-01,-1.039800d-01,
     &     -6.897000d-02,-4.738000d-02,
     &     -3.921000d-02,-3.184000d-02,
     &     -3.247000d-02,-7.074000d-02,
     &     -4.646000d-02,-9.370000d-03,
     &     1.274000d-02,5.380000d-03,
     &     -5.825000d-02,2.537325d+01,
     &     1.590210d+00,5.429600d-01,
     &     7.708000d-02,9.650000d-03,
     &     2.340000d-03,-3.190210d-04,
     &     2.590000d-03,-8.520000d-03,
     &     -8.610000d-03,-7.930000d-03,
     &     -6.960000d-03,-4.300000d-03,
     &     -1.907000d-02,2.120000d-03,
     &     -3.280000d-03,-1.574000d-02,
     &     1.981000d-02,4.940000d-03,
     &     -1.633000d-02,-1.233000d-02,
     &     -1.076000d-02,8.064000d-02,
     &     3.023800d-01,7.126000d-01,
     &     1.599000d+00,1.413750d+00,
     &     -6.690000d-03,-2.065000d-02 /)
      b2_popII_lambdab=(/ -4.360000d-03,1.027000d-01,
     &     7.540000d-03,2.900000d-03,
     &     1.290000d-03,6.137300d-04,
     &     4.232700d-04,2.856220d-04,
     &     3.080660d-04,1.780000d-03,
     &     7.493510d-04,-1.310000d-03,
     &     -2.340000d-03,-1.200000d-03,
     &     3.750000d-03,-2.032730d+00,
     &     -6.494000d-02,-1.527000d-02,
     &     -1.750000d-03,-2.319750d-04,
     &     -9.231520d-05,-3.817170d-05,
     &     -4.977780d-05,3.516460d-05,
     &     2.992460d-05,2.517400d-05,
     &     1.959700d-05,8.973120d-06,
     &     6.095290d-04,6.429410d-05,
     &     1.311010d-04,2.011070d-04,
     &     -2.279080d-04,-7.022030d-05,
     &     1.465520d-04,9.557150d-05,
     &     7.543080d-05,-2.309520d-04,
     &     -5.953970d-04,-1.210000d-03,
     &     -2.380000d-03,-1.530000d-03,
     &     1.576650d-05,1.310700d-04 /)
      b3_popII_lambdab=(/ 6.784700d-05,-4.950000d-03,
     &     -9.026520d-05,-2.248620d-05,
     &     -6.993990d-06,-2.368350d-06,
     &     -1.376460d-06,-7.912280d-07,
     &     -9.532470d-07,-1.730720d-05,
     &     -5.236220d-06,3.070040d-05,
     &     4.603600d-05,1.807760d-05,
     &     -7.591910d-05,5.430000d-02,
     &     8.695870d-04,1.383540d-04,
     &     1.199100d-05,1.262730d-06,
     &     4.677970d-07,1.807260d-07,
     &     1.695330d-07,-4.577250d-08,
     &     -3.214160d-08,-2.445600d-08,
     &     -1.679850d-08,-5.834020d-09,
     &     -8.178190d-06,-1.467830d-07,
     &     -6.036690d-07,-6.903340d-07,
     &     7.555560d-07,2.252890d-07,
     &     -5.753080d-07,-3.371170d-07,
     &     -2.411400d-07,2.219860d-07,
     &     3.917980d-07,6.846000d-07,
     &     1.178000d-06,5.475730d-07,
     &     -1.342700d-08,-3.670060d-07 /)
      b4_popII_lambdab=(/ -3.342900d-07,8.054360d-05,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     4.837890d-08,0.000000d+00,
     &     8.495490d-10,7.927130d-10,
     &     0.000000d+00,0.000000d+00,
     &     8.777110d-10,4.673670d-10,
     &     2.955430d-10,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     3.742040d-12,4.587920d-10 /)
      b5_popII_lambdab=(/ 0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     -1.045680d-10,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,0.000000d+00,
     &     0.000000d+00,-2.090690d-13 /)



* find the stage/mass that matches the input, save as ceidx
      if(kw.ge.2.and.kw.lt.4)then
         stage_num=1
      elseif(kw.eq.4.or.kw.eq.7)then
         stage_num=2
      elseif((kw.ge.4.and.kw.lt.7).or.kw.ge.8)then
         stage_num=3
      else
         celamf_nanjing = 0.5d0
         goto 90
      endif

      DO i=1, SIZE(stage)
         IF(stage(i).eq.stage_num)THEN
            if(masses(i).eq.1.and.mzams.lt.1)then
               ceidx=i
            elseif(stage(i).eq.1.and.masses(i).eq.20.and.
     &            mzams.ge.20)then
               ceidx=i
            elseif(stage(i).eq.2.and.masses(i).eq.20.and.
     &            mzams.ge.20)then
               ceidx=i
            elseif(stage(i).eq.3.and.masses(i).eq.100.and.
     &            mzams.ge.100)then
               ceidx=i
            else
               if(mzams.ge.masses(i).and.mzams.lt.masses(i+1))then
                  ceidx=i
               endif
            endif
         ENDIF
      ENDDO


* do something based on metallicity here to separate out the popI and popII
      met = 1d0 - (zpars(11)+zpars(12))
      popI_popII_Z = 0.01

      IF(met.ge.popI_popII_Z)THEN
         if(ceidx.eq.1.or.ceidx.eq.15.or.ceidx.eq.29.or.ceidx.eq.42)then
            x = menv/mzams
         else
            x = rad
         endif

         lambda_g = a_popI_lambdag(ceidx)+b1_popI_lambdag(ceidx)*x
     &     +b2_popI_lambdag(ceidx)*x**2d0
     &     +b3_popI_lambdag(ceidx)*x**3d0
     &     +b4_popI_lambdag(ceidx)*x**4d0
     &     +b5_popI_lambdag(ceidx)*x**5d0

         lambda_b = a_popI_lambdab(ceidx)+b1_popI_lambdab(ceidx)*x
     &     +b2_popI_lambdab(ceidx)*x**2d0
     &     +b3_popI_lambdab(ceidx)*x**3d0
     &     +b4_popI_lambdab(ceidx)*x**4d0
     &     +b5_popI_lambdab(ceidx)*x**5d0

         if(ceidx.eq.1.or.ceidx.eq.15.or.ceidx.eq.29.or.ceidx.eq.42)then
            lambda_g = 1d0/lambda_g
            lambda_b = 1d0/lambda_b
         elseif(ceidx.eq.31.or.ceidx.eq.32.or.ceidx.eq.33)then
            lambda_g = EXP(lambda_g)
            lambda_b = EXP(lambda_b)
         endif

      ELSEIF(met.lt.popI_popII_Z)THEN
         x = rad

         lambda_g = a_popII_lambdag(ceidx)+b1_popII_lambdag(ceidx)*x
     &     +b2_popII_lambdag(ceidx)*x**2d0
     &     +b3_popII_lambdag(ceidx)*x**3d0
     &     +b4_popII_lambdag(ceidx)*x**4d0
     &     +b5_popII_lambdag(ceidx)*x**5d0

         lambda_b = a_popII_lambdab(ceidx)+b1_popII_lambdab(ceidx)*x
     &     +b2_popII_lambdab(ceidx)*x**2d0
     &     +b3_popII_lambdab(ceidx)*x**3d0
     &     +b4_popII_lambdab(ceidx)*x**4d0
     &     +b5_popII_lambdab(ceidx)*x**5d0

         if(ceidx.eq.31.or.ceidx.eq.32)then
            lambda_g = EXP(lambda_g)
            lambda_b = EXP(lambda_b)
         endif

      ENDIF

* combine lambda_b and lambda_g according to the number lam
* lambda_g assumes the internal energy does not contribute
* lambda_b assumes 100% of the internal energy contributes
      celamf_nanjing = (lam*lambda_b)+((1d0-lam)*lambda_g)

*
 90      continue
*
      return
      end
***




