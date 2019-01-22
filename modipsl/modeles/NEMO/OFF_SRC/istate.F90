MODULE istate
   !!======================================================================
   !!                     ***  MODULE  istate  ***
   !! Ocean state   :  initial state setting
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   istate_init   : initial state setting
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers 
   USE dom_oce         ! ocean space and time domain 
   USE daymod          ! 
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC istate_init   ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/istate.F90,v 1.2 2005/11/16 16:19:34 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE istate_init
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers.
      !!
      !! ** Method  :
      !!
      !! History :
      !!   4.0  !  91-03  ()  Original code
      !!        !  91-11  (G. Madec)
      !!   9.0  !  03-09  (G. Madec)  F90: Free form, modules, orthogonality
      !!----------------------------------------------------------------------
      !! * Local declarations
      !!----------------------------------------------------------------------


      ! Initialization to zero
      ! ----------------------

      !     before fields       !       now fields        !      after fields       !
      ;   un   (:,:,:) = 0.e0   ;   ua   (:,:,:) = 0.e0
      ;   vn   (:,:,:) = 0.e0   ;   va   (:,:,:) = 0.e0
      ;                         ;   wn   (:,:,:) = 0.e0   
      ;   hdivn(:,:,:) = 0.e0   ;

      ;   tn   (:,:,:) = 0.e0   ;   ta   (:,:,:) = 0.e0
      ;   sn   (:,:,:) = 0.e0   ;   sa   (:,:,:) = 0.e0

      rhd  (:,:,:) = 0.e0
      rhop (:,:,:) = 0.e0
      rn2  (:,:,:) = 0.e0 


   END SUBROUTINE istate_init

   !!=====================================================================
END MODULE istate
