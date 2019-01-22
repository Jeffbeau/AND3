MODULE wzvmod
   !!==============================================================================
   !!                       ***  MODULE  wzvmod  ***
   !! Ocean diagnostic variable : vertical velocity
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   wzv        : Compute the vertical velocity
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC wzv       ! routine called by step.F90 and inidtr.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------

CONTAINS

#if defined key_autotasking
   !!----------------------------------------------------------------------
   !!   'key_autotasking'                               j-k-i loop (j-slab)
   !!----------------------------------------------------------------------

   SUBROUTINE wzv( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE wzv  ***
      !!                     
      !! ** Purpose :   Compute the now vertical velocity after the array swap
      !!
      !! ** Method  :   Using the incompressibility hypothesis, the vertical
      !!     velocity is computed by integrating the horizontal divergence 
      !!     from the bottom to the surface.
      !!       The boundary conditions are w=0 at the bottom (no flux) and,
      !!     in regid-lid case, w=0 at the sea surface.
      !!
      !! ** action  :    wn array : the now vertical velocity
      !!
      !! History :
      !!   5.0  !  90-10  (C. Levy, G. Madec)  Original code
      !!   7.0  !  96-01  (G. Madec)  Statement function for e3
      !!   8.5  !  02-07  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER ::   jj, jk      ! dummy loop indices
      !!----------------------------------------------------------------------
      !!  OPA 9.0 , LOCEAN-IPSL (2005) 
      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/wzvmod.F90,v 1.4 2005/09/02 15:45:24 opalod Exp $ 
      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'wzv     : vertical velocity from continuity eq.'
         IF(lwp) WRITE(numout,*) '~~~~~~~   auto-tasking case : j-k-i loop '

         ! bottom boundary condition: w=0 (set once for all)
         wn(:,:,jpk) = 0.e0
      ENDIF

      !                                                ! ===============
      DO jj = 1, jpj                                   !  Vertical slab
         !                                             ! ===============
         ! Computation from the bottom
         DO jk = jpkm1, 1, -1
            wn(:,jj,jk) = wn(:,jj,jk+1) - fse3t(:,jj,jk) * hdivn(:,jj,jk)
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      IF(ln_ctl)   CALL prt_ctl(tab3d_1=wn, clinfo1=' w**2 -   : ', mask1=wn)

   END SUBROUTINE wzv

#else
   !!----------------------------------------------------------------------
   !!   Default option                                           k-j-i loop
   !!----------------------------------------------------------------------

   SUBROUTINE wzv( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE wzv  ***
      !!
      !! ** Purpose :   Compute the now vertical velocity after the array swap
      !!
      !! ** Method  :   Using the incompressibility hypothesis, the vertical
      !!      velocity is computed by integrating the horizontal divergence 
      !!      from the bottom to the surface.
      !!        The boundary conditions are w=0 at the bottom (no flux) and,
      !!      in regid-lid case, w=0 at the sea surface.
      !!
      !! ** action  :   wn array : the now vertical velocity
      !!
      !! History :
      !!   9.0  !  02-07  (G. Madec)  Vector optimization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER ::   jk          ! dummy loop indices
      !!----------------------------------------------------------------------
      !!  OPA 8.5, LODYC-IPSL (2002)
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'wzv     : vertical velocity from continuity eq.'
         IF(lwp) WRITE(numout,*) '~~~~~~~ ' 

         ! bottom boundary condition: w=0 (set once for all)
         wn(:,:,jpk) = 0.e0
      ENDIF

      ! Computation from the bottom
      DO jk = jpkm1, 1, -1
         wn(:,:,jk) = wn(:,:,jk+1) - fse3t(:,:,jk) * hdivn(:,:,jk)
      END DO

      IF(ln_ctl)   CALL prt_ctl(tab3d_1=wn, clinfo1=' w**2 -   : ', mask1=wn)

   END SUBROUTINE wzv
#endif

   !!======================================================================
END MODULE wzvmod
