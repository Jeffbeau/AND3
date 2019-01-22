MODULE traadv_ctl
   !!==============================================================================
   !!                       ***  MODULE  traadv_ctl  ***
   !! Ocean active tracers:  advection scheme control
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   tra_adv_ctl  : control the different options of advection scheme
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC tra_adv_ctl     ! routine called by step module
 
   !! * Share module variables
   LOGICAL, PUBLIC ::   &
      ln_traadv_cen2   = .TRUE.  ,   &  ! 2nd order centered scheme flag
      ln_traadv_tvd    = .FALSE. ,   &  ! TVD scheme flag
      ln_traadv_muscl  = .FALSE. ,   &  ! MUSCL scheme flag
      ln_traadv_muscl2 = .FALSE.        ! MUSCL2 scheme flag

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traadv_ctl.F90,v 1.3 2005/12/21 10:46:45 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_adv_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_ctl  ***
      !!                
      !! ** Purpose :   Control the consistency between cpp options for 
      !!      tracer advection schemes
      !!
      !! History :
      !!   8.5  !  02-11  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ioptio

      NAMELIST/nam_traadv/ ln_traadv_cen2 , ln_traadv_tvd,   &
         &                 ln_traadv_muscl, ln_traadv_muscl2
      !!----------------------------------------------------------------------

      ! Read Namelist nam_traadv : tracer advection scheme
      ! -------------------------
      REWIND ( numnam )
      READ   ( numnam, nam_traadv )

      ! Parameter control and print
      ! ---------------------------
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_adv_ctl : choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist nam_tra_adv : chose a advection scheme for tracers'
         WRITE(numout,*)
         WRITE(numout,*) '             2nd order advection scheme     ln_traadv_cen2   = ', ln_traadv_cen2
         WRITE(numout,*) '             TVD advection scheme           ln_traadv_tvd    = ', ln_traadv_tvd
         WRITE(numout,*) '             MUSCL  advection scheme        ln_traadv_muscl  = ', ln_traadv_muscl
         WRITE(numout,*) '             MUSCL2 advection scheme        ln_traadv_muscl2 = ', ln_traadv_muscl2
      ENDIF

      ! Control of Advection scheme options
      ! -----------------------------------
      ioptio = 0
      IF( ln_traadv_cen2   )   ioptio = ioptio + 1
      IF( ln_traadv_tvd    )   ioptio = ioptio + 1
      IF( ln_traadv_muscl  )   ioptio = ioptio + 1
      IF( ln_traadv_muscl2 )   ioptio = ioptio + 1

      IF( lk_esopa ) THEN
         IF(lwp) WRITE(numout,*) ' esopa control : the use of all scheme is forced'
         ln_traadv_cen2   = .TRUE.
         ln_traadv_tvd    = .TRUE.
         ln_traadv_muscl  = .TRUE.
         ln_traadv_muscl2 = .TRUE.
      ELSEIF( ioptio > 1 .OR. ioptio == 0 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) ' Choose one advection scheme in namelist nam_traadv'
         IF(lwp) WRITE(numout,*) '        ***                              ***********'
         nstop = nstop + 1
      ENDIF

      IF( n_cla == 1 .AND. .NOT. ln_traadv_cen2 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '     cross-land advection only with 2nd order advection scheme'
         nstop = nstop + 1
      ENDIF

   END SUBROUTINE tra_adv_ctl

  !!======================================================================
END MODULE traadv_ctl
