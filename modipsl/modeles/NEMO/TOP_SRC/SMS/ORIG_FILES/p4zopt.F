
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zopt.F,v 1.6 2005/11/14 16:42:41 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
! CDIR$ LIST
      SUBROUTINE p4zopt
#if defined key_passivetrc && defined key_trc_pisces
CCC---------------------------------------------------------------------
CCC
CCC             ROUTINE p4zopt : PISCES MODEL
CCC             *****************************
CCC
CCC  PURPOSE :
CCC  ---------
CCC         Compute the light availability in the water column
CCC         depending on the depth and the chlorophyll concentration
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
CC   MODIFICATIONS:
CC   --------------
CC      original  : O. Aumont (2004)
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
! CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
#include "domzgr_substitute.h90"
! CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj, jk, mrgb
      REAL xchl,ekg(jpi,jpj,jpk),ekr(jpi,jpj,jpk)
      REAL ekb(jpi,jpj,jpk)
      REAL parlux,e1(jpi,jpj,jpk),e2(jpi,jpj,jpk),e3(jpi,jpj,jpk)
      REAL zdepmoy(jpi,jpj)
      REAL etmp(jpi,jpj)
      REAL zrlight,zblight,zglight
C
C     Initialisation of variables used to compute PAR
C     -----------------------------------------------
C
        e1     = 0.
        e2     = 0.
        e3     = 0.
        etot   = 0.
        parlux = 0.43/3.

        DO jk=1,jpkm1
          DO jj=1,jpj
            DO ji=1,jpi
C
C     Separation in three light bands: red, green, blue
C     -------------------------------------------------
C
        xchl=(trn(ji,jj,jk,jpnch)+trn(ji,jj,jk,jpdch)+rtrn)*1.E6
        xchl=max(0.03,xchl)
        xchl=min(10.,xchl)
                                                                                
        mrgb = int(41+20.*log10(xchl)+rtrn)
                                                                                
        ekb(ji,jj,jk)=xkrgb(1,mrgb)
        ekg(ji,jj,jk)=xkrgb(2,mrgb)
        ekr(ji,jj,jk)=xkrgb(3,mrgb)
C
            END DO
          END DO
        END DO
C
          DO jj = 1,jpj
            DO ji = 1,jpi
C
C     Separation in three light bands: red, green, blue
C     -------------------------------------------------
C
        zblight=0.5*ekb(ji,jj,1)*fse3t(ji,jj,1)
        zglight=0.5*ekg(ji,jj,1)*fse3t(ji,jj,1)
        zrlight=0.5*ekr(ji,jj,1)*fse3t(ji,jj,1)
C
        e1(ji,jj,1) = parlux*qsr(ji,jj)*exp(-zblight)
        e2(ji,jj,1) = parlux*qsr(ji,jj)*exp(-zglight)
        e3(ji,jj,1) = parlux*qsr(ji,jj)*exp(-zrlight)
C
            END DO
          END DO
                                                                                
                                                                                
        DO jk = 2,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
C
C     Separation in three light bands: red, green, blue
C     -------------------------------------------------
C
        zblight=0.5*(ekb(ji,jj,jk-1)*fse3t(ji,jj,jk-1)
     &    +ekb(ji,jj,jk)*fse3t(ji,jj,jk))
        zglight=0.5*(ekg(ji,jj,jk-1)*fse3t(ji,jj,jk-1)
     &    +ekg(ji,jj,jk)*fse3t(ji,jj,jk))
        zrlight=0.5*(ekr(ji,jj,jk-1)*fse3t(ji,jj,jk-1)
     &    +ekr(ji,jj,jk)*fse3t(ji,jj,jk))
C
        e1(ji,jj,jk) = e1(ji,jj,jk-1)*exp(-zblight)
        e2(ji,jj,jk) = e2(ji,jj,jk-1)*exp(-zglight)
        e3(ji,jj,jk) = e3(ji,jj,jk-1)*exp(-zrlight)
C
            END DO
          END DO
        END DO
C
        etot(:,:,:) = e1(:,:,:)+e2(:,:,:)+e3(:,:,:)
C    
C     Computation of the euphotic depth
C     ---------------------------------
C    
        zmeu(:,:) = 300.

        DO jk = 2,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
        IF (etot(ji,jj,jk).GE.0.0043*qsr(ji,jj)) THEN
           zmeu(ji,jj) = fsdepw(ji,jj,jk+1)
        ENDIF
            END DO
          END DO
        END DO
C
        zmeu(:,:)=min(300.,zmeu(:,:))
C
C    Computation of the mean light over the mixed layer depth
C    --------------------------------------------------------
C
        zdepmoy  = 0
        etmp  = 0.
        emoy  = 0.

        DO jk = 1,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
         if (fsdepw(ji,jj,jk+1).le.hmld(ji,jj)) then
       etmp(ji,jj) = etmp(ji,jj)+etot(ji,jj,jk)*fse3t(ji,jj,jk)
       zdepmoy(ji,jj)=zdepmoy(ji,jj)+fse3t(ji,jj,jk)
         endif
            END DO
          END DO
        END DO

        emoy(:,:,:) = etot(:,:,:)

        DO jk = 1,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
        IF (fsdepw(ji,jj,jk+1).LE.hmld(ji,jj)) THEN
          emoy(ji,jj,jk) = etmp(ji,jj)/(zdepmoy(ji,jj)+rtrn)
        ENDIF
            END DO
          END DO
        END DO

#   if defined key_trc_diaadd
        trc2d(:,:,11) = zmeu(:,:)
#    endif
C
#endif
      RETURN
      END
