
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zsink.F,v 1.6 2005/11/14 16:42:42 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
CDIR$ LIST
      SUBROUTINE p4zsink
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC          ROUTINE p4zsink : PISCES MODEL
CCC          ******************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC         Compute vertical flux of particulate matter due to
CCC         gravitational sinking
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
CC            p4zsink2
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
      IMPLICIT NONE
#include "domzgr_substitute.h90"
CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER jksed, ji, jj, jk
      REAL xagg1,xagg2,xagg3,xagg4
      REAL zfact,zstep,wsmax
C
C    Time step duration for biology
C    ------------------------------
C
       zstep=rfact2/rjjss
C
C    Sinking speeds of detritus is increased with depth as shown
C    by data and from the coagulation theory
C    -----------------------------------------------------------
C
       jksed=10
C
       DO jk=1,jpk-1
         DO jj=1,jpj
           DO ji=1,jpi
       zfact=max(0.,fsdepw(ji,jj,jk+1)-hmld(ji,jj))/2000.
       wsbio4(ji,jj,jk)=wsbio2+(200.-wsbio2)*zfact
             END DO
           END DO
         END DO
C
C      LIMIT THE VALUES OF THE SINKING SPEEDS 
C      TO AVOID NUMERICAL INSTABILITIES
C
      wsbio3(:,:,:)=wsbio

      Do jk=1,jpk-1
        DO jj=1,jpj
          DO ji=1,jpi
       wsmax=0.8*fse3t(ji,jj,jk)/zstep
       wsbio4(ji,jj,jk)=min(wsbio4(ji,jj,jk),wsmax)
       wsbio3(ji,jj,jk)=min(wsbio3(ji,jj,jk),wsmax)
          END DO
        END DO
       END DO

      wscal(:,:,:)=wsbio4(:,:,:)
C
C
C   INITIALIZE TO ZERO ALL THE SINKING ARRAYS
C   -----------------------------------------
C
         sinking=0.
         sinking2=0.
         sinkcal=0.
         sinkfer=0.
         sinksil=0.
         sinkfer2=0.
C
C   Compute the sedimentation term using p4zsink2 for all
C   the sinking particles
C   -----------------------------------------------------
C
         CALL p4zsink2(wsbio3,sinking,jppoc)
         CALL p4zsink2(wsbio3,sinkfer,jpsfe)
         CALL p4zsink2(wsbio4,sinking2,jpgoc)
         CALL p4zsink2(wsbio4,sinkfer2,jpbfe)
         CALL p4zsink2(wsbio4,sinksil,jpdsi)
         CALL p4zsink2(wscal,sinkcal,jpcal)
C
C  Exchange between organic matter compartments due to
C  coagulation/disaggregation
C  ---------------------------------------------------
C
         DO jk = 1,jpkm1
           DO jj = 1,jpj
             DO ji = 1,jpi
C
        zfact=zstep*zdiss(ji,jj,jk)
C
C    Part I : Coagulation dependent on turbulence
C    ----------------------------------------------
C
         xagg1=940.*zfact*trn(ji,jj,jk,jppoc)*trn(ji,jj,jk,jppoc)
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif

         xagg2=1.054E4*zfact*trn(ji,jj,jk,jppoc)*trn(ji,jj,jk,jpgoc)
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif
C
C    Aggregation of small into large particles
C    Part II : Differential settling
C    ----------------------------------------------
C
         xagg3=0.66*zstep*trn(ji,jj,jk,jppoc)*trn(ji,jj,jk,jppoc)
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif

         xagg4=0.*zstep*trn(ji,jj,jk,jppoc)*trn(ji,jj,jk,jpgoc)
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif
C
         xagg(ji,jj,jk)=xagg1+xagg2+xagg3+xagg4
         xaggfe(ji,jj,jk)=xagg(ji,jj,jk)*trn(ji,jj,jk,jpsfe)/
     &     (trn(ji,jj,jk,jppoc)+rtrn)
C
C     Aggregation of DOC to small particles
C     --------------------------------------
C
         xaggdoc(ji,jj,jk)=(80*trn(ji,jj,jk,jpdoc)+698.
     &     *trn(ji,jj,jk,jppoc))*zfact*trn(ji,jj,jk,jpdoc)
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif

         xaggdoc2(ji,jj,jk)=1.05E4*zfact*trn(ji,jj,jk,jpgoc)
     &     *trn(ji,jj,jk,jpdoc)
#    if defined key_off_degrad
     &     *facvol(ji,jj,jk)
#    endif
C
             END DO
           END DO
         END DO

#    if defined key_trc_dia3d
          trc2d(:,:,5) = sinking(:,:,jksed+1)*1.e3*rfact2r
          trc2d(:,:,6) = sinking2(:,:,jksed+1)*1.e3*rfact2r
          trc2d(:,:,7) = sinkfer(:,:,jksed+1)*1.e3*rfact2r
          trc2d(:,:,8) = sinkfer2(:,:,jksed+1)*1.e3*rfact2r
          trc2d(:,:,9) = sinksil(:,:,jksed+1)*1.e3*rfact2r
          trc2d(:,:,10) = sinkcal(:,:,jksed+1)*1.e3*rfact2r
#    endif
C
#endif
      RETURN
      END
