MODULE zdfmxl
   !!======================================================================
   !!                       ***  MODULE  zdfmxl  ***
   !! Ocean physics: mixed layer depth 
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   zdf_mxl      : Compute the turbocline and mixed layer depths.
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/ZDF/zdfmxl.F90,v 1.2 2005/11/16 16:16:03 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zdf_mxl           ! called by step.F90

   !! * Shared module variables
   INTEGER, PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      nmln                  !: number of level in the mixed layer
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
      hmld ,             &  !: mixing layer depth (turbocline) (m)
      hmlp ,             &  !: mixed layer depth  (rho=rho0+zdcrit) (m)
      hmlpt                 !: mixed layer depth at t-points (m)

   !! * module variables
   REAL(wp) ::   &
      avt_c = 5.e-4_wp,  &  ! Kz criterion for the turbocline depth
      rho_c = 0.01_wp       ! density criterion for mixed layer depth

   !! * Substitutions
#  include "domzgr_substitute.h90"

CONTAINS

# if defined key_autotasking
   !!----------------------------------------------------------------------
   !!   'key_autotasking'                               j-k-i loop (j-slab)
   !!----------------------------------------------------------------------

   SUBROUTINE zdf_mxl( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zdfmxl  ***
      !!                   
      !! ** Purpose :   Compute the turbocline depth and the mixed layer depth
      !!      with a density criteria.
      !!
      !! ** Method  :   The turbocline depth is the depth at which the vertical 
      !!      eddy diffusivity coefficient (resulting from the vertical physics
      !!      alone, not the isopycnal part, see trazdf.F) fall below a given
      !!      value defined locally (avt_c here taken equal to 5 cm/s2)
      !!
      !! ** Action  :
      !!
      !! History :
      !!   9.0  !  03-08  (G. Madec)  autotasking optimization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk     ! dummy loop indices
      INTEGER ::   ik             ! temporary integer
      INTEGER, DIMENSION(jpi,jpj) ::   &
         imld                     ! temporary workspace
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
         IF(lwp) WRITE(numout,*) '~~~~~~~   auto-tasking case : j-k-i loop'
         IF(lwp) WRITE(numout,*)
      ENDIF

      !                                                ! ===============
      DO jj = 1, jpj                                   !  Vertical slab
         !                                             ! ===============

         ! 1. Turbocline depth
         ! -------------------
         ! last w-level at which avt<avt_c (starting from the bottom jk=jpk)
         ! (since avt(.,.,jpk)=0, we have jpk=< imld =< 2 )
         DO jk = jpk, 2, -1
            DO ji = 1, jpi
               IF( avt(ji,jj,jk) < avt_c ) imld(ji,jj) = jk 
            END DO
         END DO

         ! Turbocline depth and sub-turbocline temperature
         DO ji = 1, jpi
            ik = imld(ji,jj)
            hmld (ji,jj) = fsdepw(ji,jj,ik) * tmask(ji,jj,1)
         END DO

!!gm idea
!!   
!!gm     DO jk = jpk, 2, -1
!!gm        DO ji = 1, jpi
!!gm           IF( avt(ji,jj,jk) < avt_c ) hmld(ji,jj) = fsdepw(ji,jj,jk) * tmask(ji,jj,1)
!!gm        END DO
!!gm     END DO
!!gm

         ! 2. Mixed layer depth
         ! --------------------
         ! Initialization to the number of w ocean point mbathy
         nmln(:,jj) = mbathy(:,jj)

         ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
         ! (rhop defined at t-point, thus jk-1 for w-level just above)
         DO jk = jpkm1, 2, -1
            DO ji = 1, jpi
               IF( rhop(ji,jj,jk) > rhop(ji,jj,1) + rho_c )   nmln(ji,jj) = jk
            END DO
         END DO

         ! Mixed layer depth
         DO ji = 1, jpi
            ik = nmln(ji,jj)
            hmlp (ji,jj) = fsdepw(ji,jj,ik) * tmask(ji,jj,1)
            hmlpt(ji,jj) = fsdept(ji,jj,ik-1)
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

   END SUBROUTINE zdf_mxl

# else
   !!----------------------------------------------------------------------
   !!   Default option :                                         k-j-i loop
   !!----------------------------------------------------------------------

   SUBROUTINE zdf_mxl( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdfmxl  ***
      !!                   
      !! ** Purpose :   Compute the turbocline depth and the mixed layer depth
      !!      with density criteria.
      !!
      !! ** Method  :   The turbocline depth is the depth at which the vertical
      !!      eddy diffusivity coefficient (resulting from the vertical physics
      !!      alone, not the isopycnal part, see trazdf.F) fall below a given
      !!      value defined locally (avt_c here taken equal to 5 cm/s2)
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  94-11  (M. Imbard)  Original code
      !!   8.0  !  96-01  (E. Guilyardi)  sub mixed layer temp.
      !!   8.1  !  97-07  (G. Madec)  optimization
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk     ! dummy loop indices
      INTEGER ::   ik             ! temporary integer
      INTEGER, DIMENSION(jpi,jpj) ::   &
         imld                     ! temporary workspace
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF


      ! 1. Turbocline depth
      ! -------------------
      ! last w-level at which avt<avt_c (starting from the bottom jk=jpk)
      ! (since avt(.,.,jpk)=0, we have jpk=< imld =< 2 )
      DO jk = jpk, 2, -1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( avt(ji,jj,jk) < avt_c ) imld(ji,jj) = jk 
            END DO
         END DO
      END DO

      ! Turbocline depth and sub-turbocline temperature
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = imld(ji,jj)
            hmld (ji,jj) = fsdepw(ji,jj,ik) * tmask(ji,jj,1)
         END DO
      END DO

!!gm idea
!!   
!!gm  DO jk = jpk, 2, -1
!!gm     DO jj = 1, jpj
!!gm        DO ji = 1, jpi
!!gm           IF( avt(ji,jj,jk) < avt_c ) hmld(ji,jj) = fsdepw(ji,jj,jk) * tmask(ji,jj,1)
!!gm        END DO
!!gm     END DO
!!gm  END DO
!!gm

      ! 2. Mixed layer depth
      ! --------------------
      ! Initialization to the number of w ocean point mbathy
      nmln(:,:) = mbathy(:,:)

      ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
      ! (rhop defined at t-point, thus jk-1 for w-level just above)
      DO jk = jpkm1, 2, -1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( rhop(ji,jj,jk) > rhop(ji,jj,1) + rho_c )   nmln(ji,jj) = jk
            END DO
         END DO
      END DO

      ! Mixed layer depth
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = nmln(ji,jj)
            hmlp (ji,jj) = fsdepw(ji,jj,ik) * tmask(ji,jj,1)
            hmlpt(ji,jj) = fsdept(ji,jj,ik-1)
         END DO
      END DO

   END SUBROUTINE zdf_mxl
#endif

   !!======================================================================
END MODULE zdfmxl
