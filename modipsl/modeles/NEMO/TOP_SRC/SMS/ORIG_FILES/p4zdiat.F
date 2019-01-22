
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zdiat.F,v 1.6 2005/11/14 16:42:40 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
!! CDIR$ LIST
      SUBROUTINE p4zdiat
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC           ROUTINE p4zdiat : PISCES MODEL
CCC           ******************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC         Compute the mortality terms for diatoms
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
CC           None
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : O. Aumont (2002)
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
!! CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
!! CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj, jk
      REAL zfact,zstep,compadi
C
C      Time step duration for biology
C      ------------------------------
C
        zstep=rfact2/rjjss
C
C    Aggregation term for diatoms is increased in case of nutrient
C    stress as observed in reality. The stressed cells become more
C    sticky and coagulate to sink quickly out of the euphotic zone
C     ------------------------------------------------------------
C
        DO jk = 1,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
C
        compadi = max((trn(ji,jj,jk,jpdia)-1E-8),0.)
        zfact=1./(trn(ji,jj,jk,jpdia)+rtrn)
C
C    Aggregation term for diatoms is increased in case of nutrient
C    stress as observed in reality. The stressed cells become more
C    sticky and coagulate to sink quickly out of the euphotic zone
C     ------------------------------------------------------------
C
        respp2(ji,jj,jk) = 1E6*zstep
     &    *(wchl+wchld*(1.-xlimdia(ji,jj,jk)))
     &    *zdiss(ji,jj,jk)*compadi*trn(ji,jj,jk,jpdia)
#    if defined key_off_degrad
     &    *facvol(ji,jj,jk)
#    endif
                                                                               
        respds(ji,jj,jk) = respp2(ji,jj,jk)
     &    *trn(ji,jj,jk,jpbsi)*zfact

        respdf(ji,jj,jk) = respp2(ji,jj,jk)
     &    *trn(ji,jj,jk,jpdfe)*zfact
                                                                               
        respdch(ji,jj,jk)=respp2(ji,jj,jk)
     &    *trn(ji,jj,jk,jpdch)*zfact
C
C     Phytoplankton mortality. 
C     ------------------------
C
        tortp2(ji,jj,jk) = mprat2*zstep*trn(ji,jj,jk,jpdia)
     &    /(xkmort+trn(ji,jj,jk,jpdia))*compadi
#    if defined key_off_degrad
     &    *facvol(ji,jj,jk)
#    endif

        tortds(ji,jj,jk) = tortp2(ji,jj,jk)
     &    *trn(ji,jj,jk,jpbsi)*zfact

        tortdf(ji,jj,jk)=tortp2(ji,jj,jk)
     &    *trn(ji,jj,jk,jpdfe)*zfact

        tortdch(ji,jj,jk)=tortp2(ji,jj,jk)
     &    *trn(ji,jj,jk,jpdch)*zfact
C
            END DO
          END DO
        END DO
C
#endif
      RETURN
      END
