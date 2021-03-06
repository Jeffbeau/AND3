
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zprod.F,v 1.6 2005/11/14 16:42:41 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
CDIR$ LIST
      SUBROUTINE p4zprod
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC           ROUTINE p4zprod : PISCES MODEL
CCC           ******************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC         Compute the phytoplankton production depending on
CCC         light, temperature and nutrient availability
CCC
CC   INPUT :
CC   -----
CC      argument
CC              None
CC      common
CC              all the common defined in opa
CC
CC
CC   OUTPUT :                   : no
CC   ------
CC
CC   EXTERNAL :
CC   --------
CC             p4zday
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : O. Aumont (2004) 
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
#include "domzgr_substitute.h90"
CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj, jk
      REAL silfac,pislopen(jpi,jpj,jpk),pislope2n(jpi,jpj,jpk)
      REAL zmixnano(jpi,jpj),zmixdiat(jpi,jpj),zfact
      REAL prdiachl,prbiochl,silim,ztn,zadap,zadap2
      REAL ysopt(jpi,jpj,jpk),pislopead(jpi,jpj,jpk)
      REAL prdia(jpi,jpj,jpk),prbio(jpi,jpj,jpk)
      REAL etot2(jpi,jpj,jpk),pislopead2(jpi,jpj,jpk)
      REAL xlim,silfac2,siborn,zprod,zprod2
      REAL zmxltst,zmxlday
C
C     Computation of the optimal production
C     -------------------------------------
C
C
        prmax(:,:,:)=0.6/rjjss*tgfunc(:,:,:)
#    if defined key_off_degrad
     &  *facvol(:,:,:)
#    endif
C
C     Computation of the day length
C     -----------------------------
C
        call p4zday 

        DO jk = 1,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
C
C      Computation of the P-I slope for nanos and diatoms
C      --------------------------------------------------
C
        ztn=max(0.,tn(ji,jj,jk)-15.)
        zadap=1.+2.*ztn/(2.+ztn)
        zadap2=1.

        zfact=exp(-0.21*emoy(ji,jj,jk))

        pislopead(ji,jj,jk)=pislope*(1.+zadap*zfact)
        pislopead2(ji,jj,jk)=pislope2*(1.+zadap2*zfact)

        pislopen(ji,jj,jk)=pislopead(ji,jj,jk)
     &    *trn(ji,jj,jk,jpnch)/(rtrn+trn(ji,jj,jk,jpphy)*12.)
     &    /(prmax(ji,jj,jk)*rjjss*xlimphy(ji,jj,jk)+rtrn)

        pislope2n(ji,jj,jk)=pislopead2(ji,jj,jk)
     &    *trn(ji,jj,jk,jpdch)/(rtrn+trn(ji,jj,jk,jpdia)*12.)
     &    /(prmax(ji,jj,jk)*rjjss*xlimdia(ji,jj,jk)+rtrn)
C
            END DO
          END DO
        END DO

        DO  jk = 1,jpkm1
          DO  jj = 1,jpj
            DO  ji = 1,jpi
C
C     Computation of production function
C     ----------------------------------
C
        prbio(ji,jj,jk) = prmax(ji,jj,jk)
     &    *(1.-exp(-pislopen(ji,jj,jk)*etot(ji,jj,jk)))
        prdia(ji,jj,jk) = prmax(ji,jj,jk)
     &    *(1.-exp(-pislope2n(ji,jj,jk)*etot(ji,jj,jk)))

            END DO
          END DO
        END DO

        DO  jk = 1,jpkm1
          DO  jj = 1,jpj
            DO  ji = 1,jpi
C
C    Si/C of diatoms
C    ------------------------
C    Si/C increases with iron stress and silicate availability
C    Si/C is arbitrariliy increased for very high Si concentrations
C    to mimic the very high ratios observed in the Southern Ocean
c    (silpot2)
C
C
        xlim=xdiatno3(ji,jj,jk)+xdiatnh4(ji,jj,jk)
C
        silim=min(prdia(ji,jj,jk)/(rtrn+prmax(ji,jj,jk)),
     &    trn(ji,jj,jk,jpfer)/(concdfe(ji,jj,jk)+trn(ji,jj,jk,jpfer)),
     &    trn(ji,jj,jk,jppo4)/(concdnh4+trn(ji,jj,jk,jppo4)),
     &    xlim)
        silfac=5.4*exp(-4.23*silim)+1.13
        siborn=max(0.,(trn(ji,jj,jk,jpsil)-15.E-6))
        silfac2=1.+3.*siborn/(siborn+xksi2)
        silfac=min(7.6,silfac*silfac2)
C
        ysopt(ji,jj,jk)=grosip*trn(ji,jj,jk,jpsil)/(trn(ji,jj,jk,jpsil)
     $    +xksi1)*silfac
C
            END DO
          END DO
        END DO
C
C    Computation of the limitation term due to
C    A mixed layer deeper than the euphotic depth
C    --------------------------------------------
C
        DO jj=1,jpj
          DO ji=1,jpi
         zmxltst=max(0.,hmld(ji,jj)-zmeu(ji,jj))
         zmxlday=zmxltst**2/rjjss
         zmixnano(ji,jj)=1.-zmxlday/(12.+zmxlday)
         zmixdiat(ji,jj)=1.-zmxlday/(36.+zmxlday)
          END DO
        END DO
                                                                                
        DO  jk = 1,jpkm1
          DO  jj = 1,jpj
            DO  ji = 1,jpi
         if (fsdepw(ji,jj,jk+1).le.hmld(ji,jj)) then
C
C     Mixed-layer effect on production
C     --------------------------------
C
         prbio(ji,jj,jk)=prbio(ji,jj,jk)*zmixnano(ji,jj)
         prdia(ji,jj,jk)=prdia(ji,jj,jk)*zmixdiat(ji,jj)
         endif
            END DO
          END DO
        END DO
C
        DO jk = 1,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
C
C      Computation of the maximum light intensity
C      ------------------------------------------
C
        etot2(ji,jj,jk)=etot(ji,jj,jk)*24./(strn(ji,jj)+rtrn)
        IF (strn(ji,jj).lt.1.) etot2(ji,jj,jk)=etot(ji,jj,jk)
C
            END DO
          END DO
        END DO

        DO jk = 1,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
C
C     Computation of the various production terms for nanophyto.
C     ----------------------------------------------------------
C
        prbiochl = prmax(ji,jj,jk)
     &    *(1.-exp(-pislopen(ji,jj,jk)*etot2(ji,jj,jk)))

        prorca(ji,jj,jk) = prbio(ji,jj,jk)
     &    *xlimphy(ji,jj,jk)*trn(ji,jj,jk,jpphy)*rfact2

        pronew(ji,jj,jk)=prorca(ji,jj,jk)*xnanono3(ji,jj,jk)
     &    /(xnanono3(ji,jj,jk)+xnanonh4(ji,jj,jk)+rtrn)
        proreg(ji,jj,jk)=prorca(ji,jj,jk)-pronew(ji,jj,jk)
C
        zprod=rjjss*prorca(ji,jj,jk)*prbiochl*trn(ji,jj,jk,jpphy)
     &    *xlimphy(ji,jj,jk)

        zprod2=rjjss*prorca(ji,jj,jk)*prbiochl*trn(ji,jj,jk,jpphy)
     &    *max(0.1,xlimphy(ji,jj,jk))

        prorca5(ji,jj,jk) = (fecnm)**2*zprod/chlcnm
     &    /(pislopead(ji,jj,jk)*etot2(ji,jj,jk)*trn(ji,jj,jk,jpnfe)
     &    +rtrn)

        prorca6(ji,jj,jk) = chlcnm*144.*zprod2/(pislopead(ji,jj,jk)
     &    *etot2(ji,jj,jk)*max(trn(ji,jj,jk,jpnch),1.E-10)+rtrn)

            END DO
          END DO
        END DO

        DO  jk = 1,jpkm1
          DO  jj = 1,jpj
            DO  ji = 1,jpi
C
C       Computation of the various production terms for diatoms
C       -------------------------------------------------------
C
        prdiachl = prmax(ji,jj,jk)
     &    *(1.-exp(-etot2(ji,jj,jk)*pislope2n(ji,jj,jk)))
C
        prorca2(ji,jj,jk) = prdia(ji,jj,jk)
     &    *xlimdia(ji,jj,jk)*trn(ji,jj,jk,jpdia)*rfact2
C
        pronew2(ji,jj,jk)=prorca2(ji,jj,jk)*xdiatno3(ji,jj,jk)
     &    /(xdiatno3(ji,jj,jk)+xdiatnh4(ji,jj,jk)+rtrn)
        proreg2(ji,jj,jk)=prorca2(ji,jj,jk)-pronew2(ji,jj,jk)
        prorca3(ji,jj,jk) = prorca2(ji,jj,jk)*ysopt(ji,jj,jk)
C
        zprod=rjjss*prorca2(ji,jj,jk)*prdiachl*trn(ji,jj,jk,jpdia)
     &    *max(0.1,xlimdia(ji,jj,jk))

        zprod2=rjjss*prorca2(ji,jj,jk)*prdiachl*xlimdia2(ji,jj,jk)
     &    *trn(ji,jj,jk,jpdia)

C
        prorca4(ji,jj,jk) = (fecdm)**2*zprod2/chlcdm
     &    /(pislopead2(ji,jj,jk)*etot2(ji,jj,jk)*trn(ji,jj,jk,jpdfe)
     &    +rtrn)
C
        prorca7(ji,jj,jk) = chlcdm*144.*zprod/(pislopead2(ji,jj,jk)
     &    *etot2(ji,jj,jk)*max(trn(ji,jj,jk,jpdch),1.E-10)+rtrn)
C
            END DO
          END DO
        END DO
C
#endif
      RETURN
      END

