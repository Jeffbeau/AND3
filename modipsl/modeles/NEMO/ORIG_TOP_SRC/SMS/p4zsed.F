
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zsed.F,v 1.6 2005/11/14 16:42:41 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
CDIR$ LIST
      SUBROUTINE p4zsed
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC          ROUTINE p4zsed : PISCES MODEL
CCC          *****************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC         Compute loss of organic matter in the sediments. This
CCC         is by no way a sediment model. The loss is simply 
CCC         computed to balance the inout from rivers and dust
CCC
CC   INPUT :
CC   -----
CC      common
CC              all the common defined in opa
CC
CC
CC   OUTPUT :                   : no
CC   ------
CC
CC   EXTERNAL :
CC   --------
CC             None
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : 2004 - O. Aumont 
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      USE lib_mpp
      IMPLICIT NONE
#include "domzgr_substitute.h90"
CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj, jk, ikt
      REAL sumsedsi,sumsedpo4,sumsedcal
      REAL xconctmp,denitot,nitrpottot,nitrpot(jpi,jpj,jpk)
      REAL xlim,xconctmp2,zstep,zfact
      REAL irondep(jpi,jpj,jpk),sidep(jpi,jpj)
CC
C
C     Time step duration for the biology
C     ----------------------------------
C
        zstep=rfact2/rjjss
C
C
C     Initialisation of variables used to compute deposition
C     ------------------------------------------------------
C
      irondep     = 0.
      sidep       = 0.
C
C     Iron and Si deposition at the surface
C     -------------------------------------
C
       do jj=1,jpj
         do ji=1,jpi
         irondep(ji,jj,1)=(0.014*dust(ji,jj)/(55.85*rmoss)
     &      +3E-10/raass)*rfact2/fse3t(ji,jj,1)
         sidep(ji,jj)=8.8*0.075*dust(ji,jj)*rfact2
     &      /(fse3t(ji,jj,1)*28.1*rmoss)
         end do
       end do
C
C     Iron solubilization of particles in the water column
C     ----------------------------------------------------
C
      do jk=2,jpk-1
        do jj=1,jpj
          do ji=1,jpi
          irondep(ji,jj,jk)=dust(ji,jj)/(10.*55.85*rmoss)*rfact2
     &      *0.0001
          end do
        end do
      end do
C
C    Add the external input of nutrients, carbon and alkalinity
C    ----------------------------------------------------------
C
        DO jj = 1,jpj
          DO ji = 1,jpi
          trn(ji,jj,1,jppo4) = trn(ji,jj,1,jppo4)
     &      +rivinp(ji,jj)*rfact2
          trn(ji,jj,1,jpno3) = trn(ji,jj,1,jpno3)
     &      +(rivinp(ji,jj)+nitdep(ji,jj))*rfact2
          trn(ji,jj,1,jpfer) = trn(ji,jj,1,jpfer)
     &      +rivinp(ji,jj)*9E-5*rfact2
          trn(ji,jj,1,jpsil) = trn(ji,jj,1,jpsil)
     &      +sidep(ji,jj)+cotdep(ji,jj)*rfact2/6.
          trn(ji,jj,1,jpdic) = trn(ji,jj,1,jpdic)
     &      +rivinp(ji,jj)*rfact2*2.631
          trn(ji,jj,1,jptal) = trn(ji,jj,1,jptal)
     &      +(cotdep(ji,jj)-rno3*(rivinp(ji,jj)
     &      +nitdep(ji,jj)))*rfact2
          END DO
        END DO
C
C     Add the external input of iron which is 3D distributed
C     (dust, river and sediment mobilization)
C     ------------------------------------------------------
C
        DO jk=1,jpkm1
          DO jj=1,jpj
            DO ji=1,jpi
          trn(ji,jj,jk,jpfer) = trn(ji,jj,jk,jpfer)
     &      +irondep(ji,jj,jk)+ironsed(ji,jj,jk)*rfact2
            END DO
          END DO
        END DO
C
C     Initialisation of variables used to compute Sinking Speed
C     ---------------------------------------------------------
C
        sumsedsi = 0.
        sumsedpo4 = 0.
        sumsedcal = 0.
C
C    Loss of biogenic silicon, Caco3 organic carbon in the sediments. 
C    First, the total loss is computed.
C    The factor for calcite comes from the alkalinity effect
C    -------------------------------------------------------------
C
        DO jj=2,jpjm1
          DO ji=2,jpim1
        ikt=max(mbathy(ji,jj)-1,1)
        zfact=e1t(ji,jj)*e2t(ji,jj)/rjjss
        sumsedsi=sumsedsi+trn(ji,jj,ikt,jpdsi)*wsbio4(ji,jj,ikt)
     &    *zfact
        sumsedcal=sumsedcal+trn(ji,jj,ikt,jpcal)*wscal(ji,jj,ikt)
     &    *2.*zfact
        sumsedpo4=sumsedpo4+(trn(ji,jj,ikt,jpgoc)*wsbio4(ji,jj,ikt)
     &    +trn(ji,jj,ikt,jppoc)*wsbio3(ji,jj,ikt))*zfact
          END DO
        END DO

         IF( lk_mpp ) THEN
            CALL mpp_sum( sumsedsi )   ! sums over the global domain
            CALL mpp_sum( sumsedcal )   ! sums over the global domain
            CALL mpp_sum( sumsedpo4 )   ! sums over the global domain
         ENDIF
C
C    Then this loss is scaled at each bottom grid cell for
C    equilibrating the total budget of silica in the ocean.
C    Thus, the amount of silica lost in the sediments equal
C    the supply at the surface (dust+rivers)
C    ------------------------------------------------------
C
        DO jj=1,jpj
          DO ji=1,jpi
        ikt=max(mbathy(ji,jj)-1,1)
        xconctmp=trn(ji,jj,ikt,jpdsi)*wsbio4(ji,jj,ikt)*zstep
     &    /fse3t(ji,jj,ikt)
        trn(ji,jj,ikt,jpdsi)=trn(ji,jj,ikt,jpdsi)-xconctmp
        trn(ji,jj,ikt,jpsil)=trn(ji,jj,ikt,jpsil)+xconctmp
     &    *(1.-(sumdepsi+rivalkinput/raass/6.)/sumsedsi)
          END DO
        END DO

        DO jj=1,jpj
          DO ji=1,jpi
        ikt=max(mbathy(ji,jj)-1,1)
        xconctmp=trn(ji,jj,ikt,jpcal)*wscal(ji,jj,ikt)*zstep
     &    /fse3t(ji,jj,ikt)
        trn(ji,jj,ikt,jpcal)=trn(ji,jj,ikt,jpcal)-xconctmp
        trn(ji,jj,ikt,jptal)=trn(ji,jj,ikt,jptal)+xconctmp
     &    *(1.-(rivalkinput/raass)/sumsedcal)*2.
        trn(ji,jj,ikt,jpdic)=trn(ji,jj,ikt,jpdic)+xconctmp
     &    *(1.-(rivalkinput/raass)/sumsedcal)
         END DO
       END DO

        DO jj=1,jpj
          DO ji=1,jpi
        ikt=max(mbathy(ji,jj)-1,1)
        xconctmp=trn(ji,jj,ikt,jpgoc)
        xconctmp2=trn(ji,jj,ikt,jppoc)
        trn(ji,jj,ikt,jpgoc)=trn(ji,jj,ikt,jpgoc)
     &    -xconctmp*wsbio4(ji,jj,ikt)*zstep/fse3t(ji,jj,ikt)
        trn(ji,jj,ikt,jppoc)=trn(ji,jj,ikt,jppoc)
     &    -xconctmp2*wsbio3(ji,jj,ikt)*zstep/fse3t(ji,jj,ikt)
        trn(ji,jj,ikt,jpdoc)=trn(ji,jj,ikt,jpdoc)
     &    +(xconctmp*wsbio4(ji,jj,ikt)+xconctmp2*wsbio3(ji,jj,ikt))
     &    *zstep/fse3t(ji,jj,ikt)*(1.-rivpo4input
     &    /(raass*sumsedpo4))
        trn(ji,jj,ikt,jpbfe)=trn(ji,jj,ikt,jpbfe)
     &    -trn(ji,jj,ikt,jpbfe)*wsbio4(ji,jj,ikt)*zstep
     &    /fse3t(ji,jj,ikt)
        trn(ji,jj,ikt,jpsfe)=trn(ji,jj,ikt,jpsfe)
     &    -trn(ji,jj,ikt,jpsfe)*wsbio3(ji,jj,ikt)*zstep
     &    /fse3t(ji,jj,ikt)
          END DO
        END DO
C
C  Nitrogen fixation (simple parameterization). The total gain
C  from nitrogen fixation is scaled to balance the loss by 
C  denitrification
C  -------------------------------------------------------------
C
        denitot=0.
        DO jk=1,jpk-1
          DO jj=2,jpj-1
            DO ji=2,jpi-1
        denitot=denitot+denitr(ji,jj,jk)*rdenit*e1t(ji,jj)*e2t(ji,jj)
     &    *fse3t(ji,jj,jk)*tmask(ji,jj,jk)*znegtr(ji,jj,jk)
            END DO
          END DO
        END DO

        IF( lk_mpp )   CALL mpp_sum( denitot )  ! sum over the global domain
C
C  Potential nitrogen fication dependant on temperature
C  and iron
C  ----------------------------------------------------
C
       DO jk=1,jpk
        DO jj=1,jpj
          DO ji=1,jpi
        xlim=(1.-xnanono3(ji,jj,jk)-xnanonh4(ji,jj,jk))
        if (xlim.le.0.2) xlim=0.01
        nitrpot(ji,jj,jk)=max(0.,(prmax(ji,jj,jk)-2.15/rjjss))
     &    *xlim*rfact2*trn(ji,jj,jk,jpfer)/(conc3
     &    +trn(ji,jj,jk,jpfer))*(1.-exp(-etot(ji,jj,jk)/50.))
          END DO
        END DO 
       END DO
C
        nitrpottot=0.
      DO jk=1,jpkm1
        DO jj=2,jpj-1
          DO ji=2,jpi-1
        nitrpottot=nitrpottot+nitrpot(ji,jj,jk)*e1t(ji,jj)
     &    *e2t(ji,jj)*tmask(ji,jj,jk)*fse3t(ji,jj,jk)
          END DO
        END DO
      END DO

        IF( lk_mpp )   CALL mpp_sum( nitrpottot )  ! sum over the global domain
C
C  Nitrogen change due to nitrogen fixation
C  ----------------------------------------
C
       DO jk=1,jpk
        DO jj=1,jpj
          DO ji=1,jpi
        zfact=nitrpot(ji,jj,jk)*1.E-7
        trn(ji,jj,jk,jpnh4)=trn(ji,jj,jk,jpnh4)+zfact
        trn(ji,jj,jk,jpoxy)=trn(ji,jj,jk,jpoxy)+zfact*o2nit
        trn(ji,jj,jk,jppo4)=trn(ji,jj,jk,jppo4)+30./46.*zfact
          END DO
        END DO
       END DO
C
#    if defined key_trc_diaadd
        DO jj = 1,jpj
          DO ji = 1,jpi
        trc2d(ji,jj,13) = nitrpot(ji,jj,1)*1E-7*fse3t(ji,jj,1)*1E3
     &    /rfact2
        trc2d(ji,jj,12) = irondep(ji,jj,1)*1e3*rfact2r
     &    *fse3t(ji,jj,1)
          END DO
        END DO
#    endif
C
#endif
      RETURN
      END
