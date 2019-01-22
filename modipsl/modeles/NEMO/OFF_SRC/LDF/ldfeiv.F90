MODULE ldfeiv
   !!======================================================================
   !!                     ***  MODULE  ldfeiv  ***
   !! Ocean physics:  variable eddy induced velocity coefficients
   !!======================================================================
#if   defined key_traldf_eiv   &&   defined key_traldf_c2d
   !!----------------------------------------------------------------------
   !!   'key_traldf_eiv'      and                     eddy induced velocity
   !!   'key_traldf_c2d'                    2D tracer lateral  mixing coef.
   !!----------------------------------------------------------------------
   !!   ldf_eiv      : compute the eddy induced velocity coefficients
   !!                  Same results but not same routine if 'key_autotasking'
   !!                  is defined or not
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC ldf_eiv               ! routine called by step.F90
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/LDF/ldfeiv.F90,v 1.1 2006/04/26 09:32:32 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------

CONTAINS

# if defined key_autotasking
   !!----------------------------------------------------------------------
   !!   'key_autotasking' :                            autotasking (j-slab)
   !!----------------------------------------------------------------------

   SUBROUTINE ldf_eiv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv  ***
      !!
      !! ** Purpose :   Compute the eddy induced velocity coefficient from the
      !!      growth rate of baroclinic instability.
      !!
      !! ** Method : Specific to the offline model. Computes the horizontal
      !!             values from the vertical value
      !!
      !! History :
      !!   9.0  !  06-03  (O. Aumont)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step inedx

      !! * Local declarations
      INTEGER ::   ji, jj, jk           ! dummy loop indices
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ldf_eiv : eddy induced velocity coefficients'
         IF(lwp) WRITE(numout,*) '~~~~~~~   key_autotasking'
      ENDIF

      ! Average the diffusive coefficient at u- v- points
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            aeiu(ji,jj) = .5 * (aeiw(ji,jj) + aeiw(ji+1,jj  ))
            aeiv(ji,jj) = .5 * (aeiw(ji,jj) + aeiw(ji  ,jj+1))
         END DO
      END DO
      !,,,,,,,,,,,,,,,,,,,,,,,,,,,,,synchro,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

      ! lateral boundary condition on aeiu, aeiv
      CALL lbc_lnk( aeiu, 'U', 1. )
      CALL lbc_lnk( aeiv, 'V', 1. )

   END SUBROUTINE ldf_eiv

# else
   !!----------------------------------------------------------------------
   !!   Default key                                             k-j-i loops
   !!----------------------------------------------------------------------

   SUBROUTINE ldf_eiv( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_eiv  ***
      !!
      !! ** Purpose :   Compute the eddy induced velocity coefficient from the
      !!      growth rate of baroclinic instability.
      !!
      !! ** Method : Specific to the offline model. Computes the horizontal
      !!             values from the vertical value
      !!
      !! History :
      !!   9.0  !  06-03  (O. Aumont)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step inedx

      !! * Local declarations
      INTEGER ::   ji, jj, jk           ! dummy loop indices
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ldf_eiv : eddy induced velocity coefficients'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! Average the diffusive coefficient at u- v- points
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! vector opt.
            aeiu(ji,jj) = .5 * ( aeiw(ji,jj) + aeiw(ji+1,jj  ) )
            aeiv(ji,jj) = .5 * ( aeiw(ji,jj) + aeiw(ji  ,jj+1) )
         END DO
      END DO

      ! lateral boundary condition on aeiu, aeiv
      CALL lbc_lnk( aeiu, 'U', 1. )
      CALL lbc_lnk( aeiv, 'V', 1. )

   END SUBROUTINE ldf_eiv

# endif

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Dummy module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE ldf_eiv( kt )       ! Empty routine
      WRITE(*,*) 'ldf_eiv: You should not have seen this print! error?', kt
   END SUBROUTINE ldf_eiv
#endif

   !!======================================================================
END MODULE ldfeiv


