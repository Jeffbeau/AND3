MODULE step
   !!======================================================================
   !!                       ***  MODULE step  ***
   !! Time-stepping    : manager of the ocean, tracer and ice time stepping
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   stp            : OPA system time-stepping
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE ldftra_oce
   USE in_out_manager  ! I/O manager
   USE lbclnk
   USE trcdia
   USE daymod          ! calendar                         (day     routine)
   USE trcstp          ! passive tracer time-stepping      (trc_stp routine)
   USE dtadyn          ! Lecture and interpolation of the dynamical fields
   USE eosbn2          ! equation of state                (eos_bn2 routine)
   USE trcrst          ! restart for passive tracers
   USE stpctl          ! time stepping control            (stp_ctl routine)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC stp            ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/step.F90,v 1.2 2005/11/16 16:19:34 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE stp( kstp )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE stp  ***
      !!                      
      !! ** Purpose : - Time stepping of OPA (momentum and active tracer eqs.)
      !!              - Time stepping of LIM (dynamic and thermodynamic eqs.)
      !!              - Tme stepping  of TRC (passive tracer eqs.)
      !! 
      !! ** Method  : -1- Update forcings and data  
      !!              -2- Update ocean physics 
      !!              -3- Compute the t and s trends 
      !!              -4- Update t and s 
      !!              -5- Compute the momentum trends
      !!              -6- Update the horizontal velocity
      !!              -7- Compute the diagnostics variables (rd,N2, div,cur,w)
      !!              -8- Outputs and diagnostics
      !!
      !! History :
      !!        !  91-03  ()  Original code
      !!        !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)  add a first output record
      !!        !  96-04  (G. Madec)  introduction of dynspg
      !!        !  96-04  (M.A. Foujols)  introduction of passive tracer
      !!   8.0  !  97-06  (G. Madec)  new architecture of call
      !!   8.2  !  97-06  (G. Madec, M. Imbard, G. Roullet)  free surface
      !!   8.2  !  99-02  (G. Madec, N. Grima)  hpg implicit
      !!   8.2  !  00-07  (J-M Molines, M. Imbard)  Open Bondary Conditions
      !!   9.0  !  02-06  (G. Madec)  free form, suppress macro-tasking
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kstp   ! ocean time-step index

      !! * local declarations
      INTEGER ::   indic    ! error indicator if < 0
      !! ---------------------------------------------------------------------

      indic = 1                    ! reset to no error condition

      CALL day( kstp )             ! Calendar

      CALL dta_dyn( kstp )          ! Interpolation of the dynamical fields

#if defined key_passivetrc
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Passive Tracer Model
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      ! N.B. ua, va, ta, sa arrays are used as workspace in this section
      !-----------------------------------------------------------------------

                             CALL trc_stp( kstp, indic)           ! time-stepping


#endif

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Control, diagnostics and outputs
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !                                            ! Time loop: control and print
                       CALL stp_ctl( kstp )


   END SUBROUTINE stp

   !!======================================================================
END MODULE step
