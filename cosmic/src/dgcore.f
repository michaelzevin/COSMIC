***
      SUBROUTINE dgcore(kw1,kw2,kw3,m1,m2,m3,ebinde)
      IMPLICIT NONE
      INCLUDE 'const_bse.h'
*
* A routine to determine the outcome of a collision or coalescence
* of two degenerate cores.
* Entered with kw1,kw2 = 2 or 3 with M <= Mflash, 6, 10, 11 or 12
*
*
      integer kw1,kw2,kw3
*
      real*8 m1,m2,m3,ebinde
      real*8 r1,r2,r3,mhe,mc,mne,ebindi,ebindf,deleb,de,enuc
      real*8 temp,x,y,m0,mflash
      real*8 cvhe,cvc,cvne
      parameter(cvhe=3.1d+07,cvc=8.27d+06,cvne=7.44d+06)
      real*8 ehe,ec,ene
      parameter(ehe=5.812d+17,ec=2.21d+17,ene=2.06d+17)
      real*8 the,tc,gmr,mch
      parameter(the=1.0d+08,tc=1.0d+09,gmr=1.906d+15,mch=1.44d0)
*
      real*8 corerd
      external corerd
*
* Calculate the core radii setting m0 < mflash using dummy values as we
* know it to be true if kw = 2 or 3.
      m0 = 1.d0
      mflash = 2.d0
      r1 = corerd(kw1,m1,m0,mflash)
      r2 = corerd(kw2,m2,m0,mflash)
      r3 = corerd(kw3,m3,m0,mflash)
* Calculate the initial binding energy of the two seperate cores and the
* difference between this and the final binding energy of the new core.
      ebindi = m1*m1/r1 + m2*m2/r2
      ebindf = m3*m3/r3
      deleb = ABS(ebindi - ebindf)
* If an envelope is present reduce its binding energy by the amount
* of energy liberated by the coalescence.
      ebinde = MAX(0.d0,ebinde - deleb)
      if(kw1.gt.3) goto 90
* Distribute the mass into core mass groups where mhe represents Helium
* core mass, mc represents Carbon core mass and mne represents a Neon
* core mass or any mass that is all converted Carbon.
      mhe = 0.d0
      mc = 0.d0
      mne = 0.d0
      if(kw1.le.3.or.kw1.eq.10)then
         mhe = mhe + m1
      elseif(kw1.eq.12)then
         mne = mne + m1
      else
         mc = mc + m1
      endif
      if(kw2.le.3.or.kw2.eq.10)then
         mhe = mhe + m2
      elseif(kw2.eq.12)then
         mne = mne + m2
      else
         mc = mc + m2
      endif
* Calculate the temperature generated by the merging of the cores.
      temp = (deleb/(cvhe*mhe+cvc*mc+cvne*mne))*gmr
*
* To decide if He is converted to C we use:
*    3He4 -> C12 , T > 10^8 K , 7.274 Mev released,
* to decide if C is converted further we use:
*    2C12 -> Ne20 + alpha , T > 10^9 K , 4.616 Mev released.
* and to decide if O is converted further we use:
*    2O16 -> P31 + p , T > 10^9 K , 7.677 Mev released.
* To obtain the heat capacity of an O/Ne WD and to gain an idea of the
* energy released from the further processing of an O/Ne WD we use:
*    2Ne20 + gamma -> O16 + Mg24 +gamma , T > 10^9 K , 4.583 Mev released.
* For a CO core the composition is assumed to be 20% C, 80% O and for
* an ONe core 80% O, 20% Ne.
*
* Decide if conversion of elements can take place.
*     if(temp.gt.the)then
         x = 1.d0
*     else
*        x = 0.d0
*     endif
*     if(temp.gt.tc)then
*        y = 1.d0
*     else
         y = 0.d0
*     endif
* Calculate the nuclear energy generated from any element conversion.
      enuc = (x*ehe*mhe + y*(ec*mc + ene*mne))/gmr
* Calculate the difference between the binding energy of the star
* (core + envelope) and the nuclear energy. The new star will be
* destroyed in a SN if dE < 0.
      de = (ebindf + ebinde) - enuc
* If the star survives and an envelope is present then reduce the
* envelope binding energy by the amount of liberated nuclear energy.
* The envelope will not survive if its binding energy is reduced to <= 0.
      if(de.ge.0.d0) ebinde = MAX(0.d0,ebinde - enuc)
* Now determine the final evolution state of the merged core.
      if(de.lt.0.d0) kw3 = 15
      if(kw3.eq.3)then
         if(x.gt.0.5d0)then
            kw3 = 6
         elseif(ebinde.le.0.d0)then
            kw3 = 10
         endif
      elseif(kw3.eq.4)then
         if(x.gt.0.5d0)then
            kw3 = 6
         elseif(ebinde.le.0.d0)then
            kw3 = 7
         endif
      endif
      if(kw3.eq.6.and.y.lt.0.5d0)then
         if(ebinde.le.0.d0) kw3 = 11
      elseif(kw3.eq.6.and.y.gt.0.5d0)then
         if(ebinde.le.0.d0) kw3 = 12
      endif
      if(kw3.eq.10.and.x.gt.0.5d0) kw3 = 11
      if(kw3.eq.11.and.y.gt.0.5d0) kw3 = 12
      if(kw3.ge.10.and.kw3.le.12.and.m3.ge.mch) kw3 = 15
*
      if(kw3.eq.15) m3 = 0.d0
 90   continue
*
      return
      end
***
