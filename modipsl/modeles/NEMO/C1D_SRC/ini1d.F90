MODULE ini1d
   !!======================================================================
   !!                     ***  MODULE  ini1D  ***
   !! Ocean state   :  1D initialization
   !!=====================================================================
#if defined key_cfg_1d
   !!----------------------------------------------------------------------
   !!   'key_cfg_1d'               1D Configuration
   !!----------------------------------------------------------------------
   !!   init_1d   : initial mask
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain 
   USE phycst
   USE in_out_manager

   IMPLICIT NONE
   PRIVATE

   !! * Share Module variables
   LOGICAL, PUBLIC, PARAMETER ::  lk_cfg_1d = .TRUE.       !: 1D flag

   !! * Routine accessibility
   PUBLIC init_1d   ! routine called by OPA.F90

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/ini1d.F90,v 1.2 2005/09/02 15:33:59 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE init_1d
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE init_1d  ***
      !! 
      !! ** Purpose :   Re-Initialization of masks on 1D configuration
      !!
      !! ** Method  :
      !!
      !! History :
      !!   9.0  !  04-09  (C. Ethe) 1D configuration
      !!----------------------------------------------------------------------
      !! * Local declarations
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'init_1d :  masks on 1D configuration'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'

      ! set umask and vmask equal tmask in 1D configuration
      umask(:,:,:) = tmask(:,:,:)
      vmask(:,:,:) = tmask(:,:,:)     

   END SUBROUTINE init_1d

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO 1D Config
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_cfg_1d = .FALSE.    !: internal damping flag
CONTAINS
   SUBROUTINE init_1d       ! Empty routine

   END SUBROUTINE init_1d
#endif

   !!=====================================================================
END MODULE ini1d
