MODULE ocfzpt
   !!======================================================================
   !!                       ***  MODULE  ocfzpt  ***
   !! Ocean active tracers:  freezing point computation and freezing area
   !!======================================================================
   
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain 

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC oc_fz_pt        ! called by opa.F90 and step.F90

   !! * Shared module variables   
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      freeze, freezn,  &  !: after and now ice mask (0 or 1)
      fzptb, fzptn        !: before and now freezing point
   !!----------------------------------------------------------------------
   !! OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/ocfzpt.F90,v 1.1.1.1 2005/11/14 10:41:07 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE oc_fz_pt
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE oc_fz_pt  ***
      !!
      !! ** Purpose : - Calculate ocean surface freezing temperature
      !!              - Calculate related boolean for now ice distribution 
      !!
      !! ** Method  :   Caution, freezing point only for jackett & McDougall eos
      !!
      !! ** Action  : - fzptn  : now freezing temperature at ocean surface
      !!              - fzptb  : before freezing temperature at ocean surface
      !!              - freezn : boolean indicating freezing conditions at the 
      !!                ocean surface at time step "now"
      !!
      !! History :
      !!        !  94-08  (M.-A. Filiberti)  Original code
      !!   8.5  !  02-08  (G. Madec)  F90: Free form
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::  ji, jj      ! dummy loop indices
      !!----------------------------------------------------------------------      

!CDIR NOVERRCHK
      DO jj = 1, jpj
!CDIR NOVERRCHK
         DO ji = 1, jpi
            fzptb(ji,jj) = fzptn(ji,jj)        ! swap freezing point array

            !                                  ! ocean local surface freezing temperature
            !                                  ! sn >= 0 : this is verified in stpctl at each time step
            fzptn (ji,jj) = ( -0.0575 + 1.710523e-3 * SQRT( sn(ji,jj,1) )   &
                                      - 2.154996e-4 *       sn(ji,jj,1)   ) * sn(ji,jj,1)   !!   &
            !!                        - 7.53e-4 * pressure

            !                                  ! Define boolean related to freezing conditions
            freezn(ji,jj) = tmask(ji,jj,1) * MAX( 0., SIGN( 1., fzptn(ji,jj) - tn(ji,jj,1) )  )
         END DO
      END DO

   END SUBROUTINE oc_fz_pt

   !!======================================================================
END MODULE ocfzpt
