CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/trcexp.F,v 1.5 2005/11/14 12:30:28 opalod Exp $ 
      SUBROUTINE trcexp(kt)
#if defined key_passivetrc && defined key_trc_lobster1 
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcexp
CCC                     ******************
CCC
CC
CC     PURPOSE.
CC     --------
CC          *TRCEXP* MODELS EXPORT OF BIOGENIC MATTER (POC ''SOFT
CC                   TISSUE'') AND ITS DISTRIBUTION IN WATER COLUMN
CC
CC     METHOD.
CC     -------
CC          IN THE SURFACE LAYER POC IS PRODUCED ACCORDING TO
CC     NURTRIENTS AVAILABLE AND GROWTH CONDITIONS. NUTRIENT UPTAKE
CC     KINETICS FOLLOW MICHAELIS-MENTON FORMULATION. 
CC     THE TOTAL PARTICLE AMOUNT PRODUCED, IS DISTRIBUTED IN THE WATER
CC     COLUMN BELOW THE SURFACE LAYER.
CC
CC     EXTERNALS.
CC     ----------
CC          NONE.
CC
CC     REFERENCE.
CC     ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 1999    O. Aumont
CC      modifications : 1999    C. Le Quere
CC      additions   : 01-05 (O. Aumont, E. Kestenare):
CC                           add sediment computations
CC                  :  05-06  (AS. Kremeur) new temporal integration for sedpoc
CC ---------------------------------------------------------------------
c ------
CC parameters and commons
CC ======================
! CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      USE lbclnk
      USE trc
      USE trctrp_lec

      IMPLICIT NONE
! CDIR$ LIST
CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER kt
      INTEGER ji, jj, jk, zkbot(jpi,jpj)
      REAL zwork(jpi,jpj), zgeolpoc, zfact
CC----------------------------------------------------------------------
CC statement functions
CC ===================
! CDIR$ NOLIST
#include "domzgr_substitute.h90"
! CDIR$ LIST
C
C VERTICAL DISTRIBUTION OF NEWLY PRODUCED BIOGENIC
C POC IN THE WATER COLUMN
C (PARTS OF NEWLY FORMED MATTER REMAINING IN THE DIFFERENT
C LAYERS IS DETERMINED BY DMIN3 DEFINED IN common.passivetrc.*.h
C ----------------------------------------------------------------------
C
C
      DO jk = 1,jpkm1
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
            tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3)+
     &          (1./fse3t(ji,jj,jk))*
     &          dmin3(ji,jj,jk) *fbod(ji,jj)
          ENDDO
        ENDDO
      ENDDO
C
C     Find the last level of the water column
C     Compute fluxes due to sinking particles (slow)
C   
      zkbot = jpk
      zwork = 0.
C
C
      DO jk = 1,jpkm1
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
   
             IF ( tmask(ji,jj,jk) .eq. 1 .and.
     .            tmask(ji,jj,jk+1). eq. 0 ) THEN
                  zkbot(ji,jj) = jk
                  zwork(ji,jj) = vsed * trn(ji,jj,jk,jpdet)
              ENDIF
    
             ENDDO
         ENDDO
      ENDDO
C
C     Initialization
      zgeolpoc = 0.

C     Release of nutrients from the "simple" sediment
C
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
             tra(ji,jj,zkbot(ji,jj),jpno3) = 
     .          tra(ji,jj,zkbot(ji,jj),jpno3) +
     .               sedlam*sedpocn(ji,jj)/fse3t(ji,jj,zkbot(ji,jj))

C     Deposition of organic matter in the sediment
C
             zgeolpoc = zgeolpoc + sedlostpoc*sedpocn(ji,jj)*
     .                             e1t(ji,jj)*e2t(ji,jj)

             sedpoca(ji,jj) = zwork(ji,jj)*rdt +
     .                       dminl(ji,jj)*fbod(ji,jj)*rdt -
     .                       sedlam*sedpocn(ji,jj)*rdt -
     .                       sedlostpoc*sedpocn(ji,jj)*rdt
C
             ENDDO
         ENDDO
C
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
             tra(ji,jj,1,jpno3) = tra(ji,jj,1,jpno3) + zgeolpoc*
     .                            cmask(ji,jj)/areacot/fse3t(ji,jj,1)
           ENDDO
         ENDDO

         CALL lbc_lnk( sedpocn, 'T', 1. )
 
C Oa & Ek: diagnostics depending on jpdia2d
C          left as example
#     if defined key_trc_diaadd
           do jj=1,jpj
             do ji=1,jpi
              trc2d(ji,jj,19)=sedpocn(ji,jj)
             end do
           end do
#     endif

c      ! 1. Leap-frog scheme (only in explicit case, otherwise the 
c      ! -------------------  time stepping is already done in trczdf)
       IF(l_trczdf_exp .AND. (ln_trcadv_cen2 .OR. ln_trcadv_tvd)) THEN
         zfact = 2. * rdttra(jk) * FLOAT(ndttrc) 
         IF( neuler == 0 .AND. kt == nittrc000 ) 
     .     zfact = rdttra(jk) * FLOAT(ndttrc) 
         sedpoca(:,:) = ( sedpocb(:,:) + zfact * sedpoca(:,:) )
      ENDIF

      
c      ! 2. Time filter and swap of arrays
c      ! ---------------------------------
      IF ( ln_trcadv_cen2 .OR. ln_trcadv_tvd  ) THEN         
          IF( neuler == 0 .AND. kt == nittrc000 ) THEN
              DO jj = 1, jpj
                DO ji = 1, jpi
                  sedpocb(ji,jj) = sedpocn(ji,jj)
                  sedpocn(ji,jj) = sedpoca(ji,jj)
                  sedpoca(ji,jj) = 0.
                END DO
              END DO
         ELSE
             DO jj = 1, jpj
               DO ji = 1, jpi
                 sedpocb(ji,jj) = atfp*(sedpocb(ji,jj)+sedpoca(ji,jj)) 
     .                          + atfp1 * sedpocn(ji,jj)
                 sedpocn(ji,jj) = sedpoca(ji,jj)
                 sedpoca(ji,jj) = 0.
               END DO
             END DO
         ENDIF
         
      ELSE
c         !  case of smolar scheme or muscl
         DO jj = 1, jpj
            DO ji = 1, jpi
               sedpocb(ji,jj) = sedpoca(ji,jj)
               sedpocn(ji,jj) = sedpoca(ji,jj)
               sedpoca(ji,jj) = 0.
            END DO
         END DO
         
      ENDIF

#endif
      RETURN
      END
