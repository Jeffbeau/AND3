MODULE stpctl
   !!==============================================================================
   !!                       ***  MODULE  stpctl  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   stp_ctl      : Control the run
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC stp_ctl           ! routine called by step.F90
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE stp_ctl( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl  ***
      !!                     
      !! ** Purpose :   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Print solver statistics in numsol 
      !!              - Stop the run IF problem for the solver ( indec < 0 )
      !!
      !! History :
      !!        !  91-03  ()
      !!        !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)
      !!        !  97-06  (A.M. Treguier)
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index

      !!----------------------------------------------------------------------
      !!   OPA 9.0 , LOCEAN-IPSL  (2005)
      !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/stpctl.F90,v 1.2 2005/11/16 16:19:34 opalod Exp $
      !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
      !!----------------------------------------------------------------------

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stp_ctl : time-stepping control'
         WRITE(numout,*) '~~~~~~~'
         ! open time.step file
         CALL ctlopn( numstp, 'time.step', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
      ENDIF

      ! save the current time step in numstp
      ! ------------------------------------
      IF(lwp) WRITE(numstp,9100) kt
      IF(lwp) REWIND(numstp)
9100  FORMAT(1x, i8)


   END SUBROUTINE stp_ctl

   !!======================================================================
END MODULE stpctl
