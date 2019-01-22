!!DB: 2009.08.31 -- eliminated most code
MODULE closea
   !!======================================================================
   !!                       ***  MODULE  closea  ***
   !! Closed Seas  : 
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dom_clo    : modification of the ocean domain for closed seas cases
   !!   flx_clo    : Special handling of closed seas
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE ocesbc          ! ocean surface boundary conditions (fluxes)
   USE flxrnf          ! runoffs
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ???

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dom_clo      ! routine called by dom_init
   PUBLIC flx_clo      ! routine called by step

   !! * Share module variables
   INTEGER, PUBLIC, PARAMETER ::   &  !:
!      jpncs   = 4               !: number of closed sea
      jpncs   = 0               !: number of closed sea
   INTEGER, PUBLIC ::          & !!: namclo : closed seas and lakes
      nclosea =  0                !: = 0 no closed sea or lake
      !                           !  = 1 closed sea or lake in the domain
   INTEGER, PUBLIC, DIMENSION (jpncs) ::   &  !:
      ncstt,           &  !: Type of closed sea
      ncsi1, ncsj1,    &  !: closed sea limits                                                                 
      ncsi2, ncsj2,    &  !: 
      ncsnr               !: number of point where run-off pours
   INTEGER, PUBLIC, DIMENSION (jpncs,4) ::   &
      ncsir, ncsjr        !: Location of run-off

   !! * Module variable
   REAL(wp), DIMENSION (jpncs+1) ::   &
      surf               ! closed sea surface

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/closea.F90,v 1.5 2005/09/22 14:23:57 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dom_clo
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dom_clo  ***
      !!        
      !! ** Purpose :   Closed sea domain initialization
      !!
      !! ** Method  :   if a closed sea is located only in a model grid point
      !!      just the thermodynamic processes are applied.
      !!
      !! ** Action :   ncsi1(), ncsj1() : south-west closed sea limits (i,j)
      !!               ncsi2(), ncsj2() : north-east Closed sea limits (i,j)
      !!               ncsir(), ncsjr() : Location of runoff
      !!               ncsnr            : number of point where run-off pours
      !!               ncstt            : Type of closed sea
      !!                                  =0 spread over the world ocean
      !!                                  =2 put at location runoff
      !!
      !! History :
      !!        !  01-04  (E. Durand)  Original code
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Local variables
      INTEGER ::   jc            ! dummy loop indices
      !!----------------------------------------------------------------------
      

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*)'DB: dom_clo : closed seas ---> blank routine'
      IF(lwp) WRITE(numout,*)'~~~~~~~'

   END SUBROUTINE dom_clo


   SUBROUTINE flx_clo( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE flx_clo  ***
      !!                    
      !! ** Purpose :   Special handling of closed seas
      !!
      !! ** Method  :   Water flux is forced to zero over closed sea
      !!      Excess is shared between remaining ocean, or
      !!      put as run-off in open ocean.
      !!
      !! ** Action :
      !!
      !! History :
      !!   8.2  !  00-05  (O. Marti)  Original code
      !!   8.5  !  02-07  (G. Madec)  Free form, F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT (in) :: kt

      !! * Local declarations
      REAL(wp), DIMENSION (jpncs) :: zemp
      INTEGER  :: ji, jj, jc, jn
      REAL(wp) :: zze2
      !!----------------------------------------------------------------------

      ! 1 - Initialisation
      ! ------------------

      IF( kt == nit000 ) THEN 
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'DB: flx_clo : closed seas --- blank routine'
         IF(lwp) WRITE(numout,*)'~~~~~~~'

         ! Total surface of ocean
!         surf(jpncs+1) = SUM( e1t(:,:) * e2t(:,:) * tmask_i(:,:) )
         surf(jpncs+1) = 0.e0
      endif


   END SUBROUTINE flx_clo

   !!======================================================================
END MODULE closea
