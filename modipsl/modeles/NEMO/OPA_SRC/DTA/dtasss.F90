!!DB: 2009.09.02 -- disabled thie key
!!sends message and increments nstop flag if this key is on
MODULE dtasss
   !!======================================================================
   !!                       ***  MODULE  dtasss  ***
   !! Data : Sea Surface Salinity (SSS)
   !!======================================================================
   
   !!----------------------------------------------------------------------
   !!   dta_sss      : sss data
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
   PUBLIC dta_sss

   !! * Shared module variables
#if defined key_dtasss
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtasss = .TRUE.   !: sss data flag
#else
   LOGICAL , PUBLIC, PARAMETER ::   lk_dtasss = .FALSE.  !: sss data flag
#endif
!   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
!      sss             !: surface salinity

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , IPSL-LODYC  (2005)
   !!----------------------------------------------------------------------

CONTAINS

#if defined key_dtasss
   !!----------------------------------------------------------------------
   !!   'key_dtasss'                                               SSS data
   !!----------------------------------------------------------------------

   SUBROUTINE dta_sss( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dta_sss  ***
      !!                    
      !! ** Purpose :   Read surface salinity data 
      !!
      !! ** Method  : - Read a specific sss.
      !!
      !! ** Action  : - sss 
      !!
      !! History :
      !!        !  90-03  (O. Marti and Ph Dandin)  Original code
      !!        !  92-07  (M. Imbard)
      !!        !  96-11  (E. Guilyardi)  Daily AGCM input files
      !!        !  00-04  (M. Imbard)  NetCDF FORMAT
      !!        !  00-10  (J.-P. Boulanger)  passage ORCA a TDH
      !!        !  01-10  (A. Lazar)  Reynolds default
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!        !  02-11  (C. Levy)  MPP/MPI NetCDF read
      !!        !  05-03  (M. Levy) adapt SST to SSS
      !!----------------------------------------------------------------------
      
      !! * Arguments
      INTEGER ::   kt

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'DB: Defunct key dtasss defined ----> stop'
         IF(lwp) WRITE(numout,*) '    If you want this, find old code'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         nstop = nstop + 1
      ENDIF

   END SUBROUTINE dta_sss

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                        NO SSS data
   !!----------------------------------------------------------------------

   SUBROUTINE dta_sss( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dta_sss  ***
      !!                    
      !! ** Purpose :   sea surface salinity data and update it
      !!     at each time step.   ???
      !!
      !! ** Method  : - sss  
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean timestep
      !!---------------------------------------------------------------------
      
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dta_sss : No SSS data'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      
   END SUBROUTINE dta_sss
#endif

   !!======================================================================
END MODULE dtasss
