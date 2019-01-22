!!DB: 2009.09.02 -- disabled thie key
!!sends message and increments nstop flag if this key is on
MODULE dtasst
   !!======================================================================
   !!                       ***  MODULE  dtasst  ***
   !! Data : Sea Surface Temperature (SST)
   
   !!      BUG initialisation  nyearsst !!!!!!bug
   
   !!======================================================================
   
   !!----------------------------------------------------------------------
   !!   dta_sst      : Reynolds sst data
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE ocfzpt          ! ???
   USE daymod          ! calendar

   IMPLICIT NONE
   PRIVATE

   !! * Shared routine
   PUBLIC dta_sst

   !! * Shared module variables
#if defined key_dtasst
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtasst = .TRUE.   !: sst data flag
#else
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtasst = .FALSE.  !: sst data flag
#endif
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      sst             !: surface temperature
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2) ::   &  !:
      rclice          !: climatological ice index (0/1) (2 months)

CONTAINS

#if defined key_dtasst
   !!----------------------------------------------------------------------
   !!   'key_dtasst'                                               SST data
   !!----------------------------------------------------------------------

   SUBROUTINE dta_sst( kt )
      
      !! * Arguments
      INTEGER ::   kt

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'DB: Defunct key dtasst defined ----> stop'
         IF(lwp) WRITE(numout,*) '    If you want this, find old code'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         nstop = nstop + 1
      ENDIF
   END SUBROUTINE dta_sst

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                        NO SST data
   !!----------------------------------------------------------------------

   SUBROUTINE dta_sst( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dta_sst  ***
      !!                    
      !! ** Purpose :   sea surface temperature data and update it
      !!     at each time step.   ???
      !!
      !! ** Method  : - sst   = tn
      !!              - rclice = 1. IF tn =< ztgel
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean timestep
      
      !! * Local declarations
      INTEGER :: ji, jj
      !!---------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dta_sst : No SST data'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF
      
      ! 1. Update at each time step
      ! ---------------------------

      sst   (:,:)   = tn   (:,:,1)
      rclice(:,:,1) = tmask(:,:,1)
      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( tn(ji,jj,1) >= fzptn(ji,jj) ) rclice(ji,jj,1) = 0.e0
         END DO
      END DO
      rclice(:,:,2) = rclice(:,:,1)
      
   END SUBROUTINE dta_sst
#endif

   !!======================================================================
END MODULE dtasst
