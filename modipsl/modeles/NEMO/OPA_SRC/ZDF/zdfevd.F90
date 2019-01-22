MODULE zdfevd
   !!======================================================================
   !!                       ***  MODULE  zdfevd  ***
   !! Ocean physics: parameterization of convection through an enhancement
   !!                of vertical eddy mixing coefficient
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   zdf_evd      : update momentum and tracer Kz at the location of
   !!                  statically unstable portion of the water column
   !!                  (called if ln_zdfevd=T)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics variables
   USE zdfkpp          ! KPP vertical mixing
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zdf_evd      ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfevd.F90,v 1.4 2005/09/02 15:02:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE zdf_evd( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_evd  ***
      !!                   
      !! ** Purpose :   Local increased the vertical eddy viscosity and diffu-
      !!      sivity coefficients when a static instability is encountered.
      !!
      !! ** Method  :   avt, and the 4 neighbouring avmu, avmv coefficients
      !!      are set to avevd (namelist parameter) if the water column is 
      !!      statically unstable (i.e. if rn2 < -1.e-12 )
      !!
      !! ** Action  :   Update avt, avmu, avmv in statically instable cases
      !!                and avt_evd which is avt due to convection
      !! References :
      !!      Lazar, A., these de l'universite Paris VI, France, 1997
      !! History :
      !!   7.0  !  97-06  (G. Madec, A. Lazar)  Original code
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  05-06  (C. Ethe) KPP parameterization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step indexocean time step

      !! * Local declarations
      INTEGER ::   ji, jj, jk               ! dummy loop indices
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'zdf_evd : Enhanced Vertical Diffusion (evd)'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         IF(lwp) WRITE(numout,*)
      ENDIF

      ! Initialisation of avt_evd (vertical diffusion due to convection) to avt and avmu_evd to avmu
      avt_evd  (:,:,:) = avt(:,:,:) 
      avmu_evd (:,:,:) = avmu(:,:,:) 

      SELECT CASE ( nevdm )
 
      CASE ( 1 )           ! enhance vertical eddy viscosity and diffusivity (if rn2<-1.e-12)
         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
#   if defined key_vectopt_loop   &&   ! defined key_autotasking
!!!         WHERE( rn2(:,:,jk) <= -1.e-12 ) avt(:,:,jk) = tmask(:,:,jk) * avevd   ! agissant sur T SEUL!
            jj = 1                     ! big loop forced
            DO ji = jpi+2, jpij   
#   if defined key_zdfkpp
!! no implicit mixing in the boundary layer with KPP
               IF( ( rn2(ji,jj,jk) <= -1.e-12 ) .AND. ( fsdepw(ji,jj,jk) > hkpp(ji,jj) ) ) THEN
#   else
               IF( rn2(ji,jj,jk) <= -1.e-12 ) THEN
#   endif
                  avt (ji  ,jj  ,jk) = avevd * tmask(ji  ,jj  ,jk)
                  avmu(ji  ,jj  ,jk) = avevd * umask(ji  ,jj  ,jk)
                  avmu(ji-1,jj  ,jk) = avevd * umask(ji-1,jj  ,jk)
                  avmv(ji  ,jj  ,jk) = avevd * vmask(ji  ,jj  ,jk)
                  avmv(ji  ,jj-1,jk) = avevd * vmask(ji  ,jj-1,jk)
               ENDIF
            END DO
#   else
            DO jj = 2, jpj             ! no vector opt.
               DO ji = 2, jpi
#   if defined key_zdfkpp
!! no implicit mixing in the boundary layer with KPP
               IF( ( rn2(ji,jj,jk) <= -1.e-12 ) .AND. ( fsdepw(ji,jj,jk) > hkpp(ji,jj) ) ) THEN
#   else
               IF( rn2(ji,jj,jk) <= -1.e-12 ) THEN
#   endif
                     avt (ji  ,jj  ,jk) = avevd * tmask(ji  ,jj  ,jk)
                     avmu(ji  ,jj  ,jk) = avevd * umask(ji  ,jj  ,jk)
                     avmu(ji-1,jj  ,jk) = avevd * umask(ji-1,jj  ,jk)
                     avmv(ji  ,jj  ,jk) = avevd * vmask(ji  ,jj  ,jk)
                     avmv(ji  ,jj-1,jk) = avevd * vmask(ji  ,jj-1,jk)
                  ENDIF
               END DO
            END DO
#   endif
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============

         ! Lateral boundary conditions on ( avt, avmu, avmv )   (unchanged sign)
         ! -------------------------------===================
         CALL lbc_lnk( avt , 'W', 1. )
         CALL lbc_lnk( avmu, 'U', 1. )
         CALL lbc_lnk( avmv, 'V', 1. )

      CASE DEFAULT         ! enhance vertical eddy diffusivity only (if rn2<-1.e-12) 
         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
!!!         WHERE( rn2(:,:,jk) <= -1.e-12 ) avt(:,:,jk) = tmask(:,:,jk) * avevd   ! agissant sur T SEUL! 
#   if defined key_vectopt_loop   &&   ! defined key_autotasking
            jj = 1                     ! big loop forced
            DO ji = 1, jpij   
#   if defined key_zdfkpp
!! no implicit mixing in the boundary layer with KPP
               IF( ( rn2(ji,jj,jk) <= -1.e-12 ) .AND. ( fsdepw(ji,jj,jk) > hkpp(ji,jj) ) ) &              
                  avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
#   else
               IF( rn2(ji,jj,jk) <= -1.e-12 )   avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
#   endif
            END DO
#   else
            DO jj = 1, jpj             ! loop over the whole domain (no lbc_lnk call)
               DO ji = 1, jpi
#   if defined key_zdfkpp
!! no implicit mixing in the boundary layer with KPP
               IF( ( rn2(ji,jj,jk) <= -1.e-12 ) .AND. ( fsdepw(ji,jj,jk) > hkpp(ji,jj) ) ) &          
                  avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
#   else
                  IF( rn2(ji,jj,jk) <= -1.e-12 )   avt(ji,jj,jk) = avevd * tmask(ji,jj,jk)
#   endif
               END DO
            END DO
#   endif
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============
      END SELECT 

      ! update of avt_evd and avmu_evd
      avt_evd  (:,:,:) = avt (:,:,:)  - avt_evd  (:,:,:) 
      avmu_evd (:,:,:) = avmu(:,:,:)  - avmu_evd (:,:,:) 

   END SUBROUTINE zdf_evd

   !!======================================================================
END MODULE zdfevd
