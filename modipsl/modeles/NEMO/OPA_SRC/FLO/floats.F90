!!DB 2008.05.16
!!Modified routine for more specific DB particle tracking needs
!!For older code, see OLD_CODE/ or look elsewhere for an older version
!!Modifications: 
!!(1) Make jpnfl = #-floats = something that is determined from the input
!!file "init_float". This means that arrays are dynamically allocated,
!!allowing also for zero floats (if init_float does not exist),
!!and that significant code changes occur, starting in flo_oce ...
!!(2) Note the re-definition of ln_flork4:
!!if TRUE (see namelist file) then  CALL flo_RDM(kt) which is a DB-written
!!Random Displacement Model for particle tracking. Note that at this time
!!this model is written for constant z-level particle tracking. Also,
!!if ln_flork4 is FALSE then existing Blanke routine is called.
!!(3) DB has eliminated most of the original float restart stuff, and the
!!argo float routine. 
!!Also, flo4rk.F90 has not been re-coded (it did not work anyways)

MODULE floats
   !!======================================================================
   !!                       ***  MODULE  floats  ***
   !! Ocean floats : floats
   !!======================================================================
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                     float trajectories
   !!----------------------------------------------------------------------
   !!   flo_stp   : float trajectories computation
   !!   flo_init  : initialization of float trajectories computation
   !!----------------------------------------------------------------------
   !! * Modules used
   USE flo_oce         ! floats variables
   USE lib_mpp         ! distributed memory computing
   USE flodom          ! initialisation Module 
   USE flowri          ! float output                     (flo_wri routine)
   USE flo4rk          ! Trajectories, Runge Kutta scheme (flo_4rk routine)
   USE floblk          ! Trajectories, Blanke scheme      (flo_blk routine)

   IMPLICIT NONE
   PRIVATE  

   !! * Routine accessibility
   PUBLIC flo_stp    ! routine called by step.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/floats.F90,v 1.3 2005/03/27 18:35:05 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE flo_stp( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE flo_stp  ***
      !!                    
      !! ** Purpose :   Compute the geographical position (lat., long., depth)
      !!      of each float at each time step with one of the algorithm.
      !! 
      !! ** Method  :   The position of a float is computed with Bruno Blanke 
      !!        algorithm by default and with a 4th order Runge-Kutta scheme
      !!        if ln_flork4 =T
      !!      
      !! History :
      !!   8.5  !  02-06  (A. Bozec, G. Madec )  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'flo_stp : call floats routine '
         IF(lwp) WRITE(numout,*) '~~~~~~~'

         CALL flo_init           ! read the namelist of floats             

         CALL flo_dom            ! compute/read initial position of floats

         ! Initialisation of wb for computation of floats trajectories at the first time step
         wb(:,:,:) = wn(:,:,:)
      ENDIF

      if(jpnfl == 0) return 


!!DB 2008.03.19
      if( ln_flork4 ) THEN
         CALL flo_RDM( kt )       ! Trajectories using DB RDM 
      else
         CALL flo_blk( kt )        ! Trajectories using Blanke' algorithme
      endif

      IF( lk_mpp )   CALL mppsync   ! synchronization of all the processor


      ! Writing and restart      
      
      ! trajectories file 
      IF( kt == nit000 .OR. MOD( kt, nwritefl ) == 0 )   CALL flo_wri( kt )
      ! restart file 
      IF( kt == nitend .OR. MOD( kt, nstockfl ) == 0 )   CALL flo_wri( kt )

      ! Save the old vertical velocity field
      wb(:,:,:) = wn(:,:,:)

   END SUBROUTINE flo_stp


   SUBROUTINE flo_init
      !!----------------------------------------------------------------
      !!                 ***  ROUTINE flo_init  ***
      !!                   
      !! ** Purpose :   Read the namelist of floats
      !!      
      !! History :
      !!   8.0  !         (CLIPPER)   original Code
      !!   8.5  !  02-06  (A. Bozec)  F90, Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
!!DB: 2009.09.02
!      USE ioipsl

      !! * Local declarations
      NAMELIST/namflo/ ln_rstflo, nwritefl, nstockfl, ln_argo, ln_flork4
      !!---------------------------------------------------------------------
      ! Namelist namflo : floats
      
      ! default values
      ln_rstflo  = .FALSE.
      nwritefl  = 150
      nstockfl  = 450
      
      ! lecture of namflo
      REWIND( numnam )
      READ  ( numnam, namflo )

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) '         Namelist floats :'
         WRITE(numout,*) '            restart                          ln_rstflo = ', ln_rstflo
         WRITE(numout,*) '            frequency of float output file   nwritefl  = ', nwritefl
         WRITE(numout,*) '            frequency of float restart file  nstockfl  = ', nstockfl
         WRITE(numout,*) ' '
      ENDIF

   END SUBROUTINE flo_init

#  else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE flo_stp( kt )          ! Empty routine
      if(lwp) WRITE(numout,*) 'flo_stp: You should not have seen this print! error?', kt
   END SUBROUTINE flo_stp
#endif

   !!======================================================================
 END MODULE floats
