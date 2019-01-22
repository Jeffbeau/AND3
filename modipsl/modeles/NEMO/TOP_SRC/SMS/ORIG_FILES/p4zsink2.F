
CCC $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/p4zsink2.F,v 1.6 2005/11/14 16:42:42 opalod Exp $ 
CCC  TOP 1.0 , LOCEAN-IPSL (2005) 
C This software is governed by CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
C ---------------------------------------------------------------------------
      SUBROUTINE p4zsink2(wstmp,sinktemp,jn)
! CDIR$ LIST
#if defined key_passivetrc && defined key_trc_pisces
!!!
!!!       p4zsink2 : PISCES model
!!!       ***********************
!!!
!!
!!  PURPOSE :
!!  ---------
!!     Compute the sedimentation terms for the various sinking
!!     particles. The scheme used to compute the trends is based
!!     on MUSCL. 
!!
!!   METHOD :
!!   -------
!!      this ROUTINE compute not exactly the advection but the
!!      transport term, i.e.  div(u*tra).
!!
!!
!!   REFERENCES :                
!!   ----------                  
!!
!!   References :
!!      Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
!!      IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
!!
!!
!!   MODIFICATIONS:
!!   --------------
!!       original :  06-00 (A.Estublier)
!!       modifications : 2004 (O. Aumont)
!!       
!!----------------------------------------------------------------------
CC ----------------------------------------------------------------
CC parameters and commons
CC ======================
! CDIR$ NOLIST
      USE oce_trc
      USE trp_trc
      USE sms
      IMPLICIT NONE
#include "domzgr_substitute.h90"
! CDIR$ LIST
CC-----------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji,jj,jk,jn
      REAL ztraz(jpi,jpj,jpk),zakz(jpi,jpj,jpk)
      REAL zkz(jpi,jpj,jpk)
      REAL zigma,zew,zstep,zign
      REAL wstmp(jpi,jpj,jpk),sinktemp(jpi,jpj,jpk)
      REAL wstmp2(jpi,jpj,jpk)

!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (01/00)
!!!---------------------------------------------------------------------
! 1. Initialization
! --------------

        zstep  = rdt*ndttrc

        ztraz  = 0
        zkz    = 0
        zakz   = 0.

        do jk=1,jpk-1
         wstmp2(:,:,jk+1)=-wstmp(:,:,jk)/rjjss*tmask(:,:,jk+1)
#    if defined key_off_degrad
     &      *facvol(:,:,jk)
#    endif
        end do
 
        wstmp2(:,:,1)=0.
!
! 3. Vertical advective flux
!-------------------------------
! ... first guess of the slopes
!   ... interior values
!
        DO jk=2,jpkm1
              ztraz(:,:,jk) = (trn(:,:,jk-1,jn) - trn(:,:,jk,jn))
     $                          *tmask(:,:,jk)
        ENDDO
!
! slopes
!
        DO jk=2,jpkm1
          DO jj = 1,jpj
            DO ji = 1, jpi
            zign = 0.5*(sign(1.,ztraz(ji,jj,jk)*ztraz(ji,jj,jk+1))+1)
            zakz(ji,jj,jk) = 0.5*(ztraz(ji,jj,jk)
     $                          +ztraz(ji,jj,jk+1))*zign
            ENDDO
          ENDDO
        ENDDO        
!
! Slopes limitation
!
        DO jk=2,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
              zakz(ji,jj,jk) = sign(1.,zakz(ji,jj,jk)) * 
     $                        min(abs(zakz(ji,jj,jk)),
     $                        2.*abs(ztraz(ji,jj,jk+1)),
     $                        2.*abs(ztraz(ji,jj,jk)))
            ENDDO
          ENDDO
        ENDDO        

! vertical advective flux
        DO jk=1,jpkm1
          DO jj = 1,jpj      
            DO ji = 1, jpi    
              zigma = wstmp2(ji,jj,jk+1)*zstep/fse3w(ji,jj,jk+1)
              zew   = wstmp2(ji,jj,jk+1)
              sinktemp(ji,jj,jk+1) = -zew*(trn(ji,jj,jk,jn)
     $           -0.5*(1+zigma)*zakz(ji,jj,jk))*rfact2
            ENDDO
          ENDDO
        ENDDO 
!
! Boundary conditions
!
         sinktemp(:,:,1)=0.
         sinktemp(:,:,jpk)=0.
C
       DO jk=1,jpkm1
          DO jj = 1,jpj
            DO ji = 1, jpi
!
            trn(ji,jj,jk,jn) = trn(ji,jj,jk,jn)
     &        + (sinktemp(ji,jj,jk)-sinktemp(ji,jj,jk+1))
     &        /fse3t(ji,jj,jk)
!
            ENDDO
          ENDDO
        ENDDO
!
        trb(:,:,:,jn)=trn(:,:,:,jn)
!
#endif
C
      RETURN
      END
