MODULE limhdf
   !!======================================================================
   !!                    ***  MODULE limhdf   ***
   !! LIM ice model : horizontal diffusion of sea-ice quantities
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim'                                     LIM sea-ice model
   !!----------------------------------------------------------------------
   !!   lim_hdf  : diffusion trend on sea-ice variable
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce
   USE ice_oce         ! ice variables
   USE in_out_manager
   USE ice
   USE lbclnk
   USE lib_mpp
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC lim_hdf    ! called by lim_tra

   !! * Module variables
   LOGICAL  ::   linit = .TRUE.              ! ???
   REAL(wp) ::   epsi04 = 1e-04              ! constant
   REAL(wp), DIMENSION(jpi,jpj) ::   zfact   ! ???

   !! * Substitution 
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limhdf.F90,v 1.8 2005/09/22 13:50:30 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE lim_hdf( ptab )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE lim_hdf  ***
      !!
      !! ** purpose :   Compute and add the diffusive trend on sea-ice
      !!      variables
      !!
      !! ** method  :   Second order diffusive operator evaluated using a
      !!      Cranck-Nicholson time Scheme.
      !!
      !! ** Action  :    update ptab with the diffusive contribution
      !!
      !! History :
      !!        !  00-01 (LIM) Original code
      !!        !  01-05 (G. Madec, R. Hordoir) opa norm
      !!        !  02-08 (C. Ethe)  F90, free form
      !!-------------------------------------------------------------------
      ! * Arguments
      REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
         ptab                 ! Field on which the diffusion is applied  
      REAL(wp), DIMENSION(jpi,jpj) ::   &
         ptab0                ! ???

      ! * Local variables
      INTEGER ::  ji, jj      ! dummy loop indices
      INTEGER ::  &
         its, iter            ! temporary integers
      CHARACTER (len=55) :: charout
      REAL(wp) ::  &
         zalfa, zrlxint, zconv, zeps   ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj) ::  & 
         zrlx, zflu, zflv, &  ! temporary workspaces
         zdiv0, zdiv          !    "           "
      !!-------------------------------------------------------------------

      ! Initialisation
      ! ---------------   
      ! Time integration parameters
      zalfa = 0.5       ! =1.0/0.5/0.0 = implicit/Cranck-Nicholson/explicit
      its   = 100       ! Maximum number of iteration
      zeps  =  2. * epsi04

      ! Arrays initialization
      ptab0 (:, : ) = ptab(:,:)
!bug  zflu (:,jpj) = 0.e0
!bug  zflv (:,jpj) = 0.e0
      zdiv0(:, 1 ) = 0.e0
      zdiv0(:,jpj) = 0.e0
      IF( .NOT.lk_vopt_loop ) THEN
         zflu (jpi,:) = 0.e0   
         zflv (jpi,:) = 0.e0
         zdiv0(1,  :) = 0.e0
         zdiv0(jpi,:) = 0.e0
      ENDIF

      ! Metric coefficient (compute at the first call and saved in
      IF( linit ) THEN
         DO jj = 2, jpjm1  
            DO ji = fs_2 , fs_jpim1   ! vector opt.
               zfact(ji,jj) = ( e2u(ji,jj) + e2u(ji-1,jj  ) + e1v(ji,jj) + e1v(ji,jj-1) ) &
                  &          / ( e1t(ji,jj) * e2t(ji,jj) )
            END DO
         END DO
         linit = .FALSE.
      ENDIF


      ! Sub-time step loop
      zconv = 1.e0
      iter  = 0

      !                                                   !===================
      DO WHILE ( ( zconv > zeps ) .AND. (iter <= its) )   ! Sub-time step loop
         !                                                !===================
         ! incrementation of the sub-time step number
         iter = iter + 1

         ! diffusive fluxes in U- and V- direction
         DO jj = 1, jpjm1
            DO ji = 1 , fs_jpim1   ! vector opt.
               zflu(ji,jj) = pahu(ji,jj) * e2u(ji,jj) / e1u(ji,jj) * ( ptab(ji+1,jj) - ptab(ji,jj) )
               zflv(ji,jj) = pahv(ji,jj) * e1v(ji,jj) / e2v(ji,jj) * ( ptab(ji,jj+1) - ptab(ji,jj) )
            END DO
         END DO

         ! diffusive trend : divergence of the fluxes
         DO jj= 2, jpjm1
            DO ji = fs_2 , fs_jpim1   ! vector opt. 
               zdiv (ji,jj) = (  zflu(ji,jj) - zflu(ji-1,jj  )   &
                  &            + zflv(ji,jj) - zflv(ji  ,jj-1)  ) / ( e1t (ji,jj) * e2t (ji,jj) )
            END DO
         END DO

         ! save the first evaluation of the diffusive trend in zdiv0
         IF( iter == 1 )   zdiv0(:,:) = zdiv(:,:)       

         ! XXXX iterative evaluation?????
         DO jj = 2, jpjm1
            DO ji = fs_2 , fs_jpim1   ! vector opt.
               zrlxint = (   ptab0(ji,jj)    &
                  &       +  rdt_ice * (           zalfa   * ( zdiv(ji,jj) + zfact(ji,jj) * ptab(ji,jj) )   &
                  &                      + ( 1.0 - zalfa ) *   zdiv0(ji,jj) )  )                             & 
                  &    / ( 1.0 + zalfa * rdt_ice * zfact(ji,jj) )
               zrlx(ji,jj) = ptab(ji,jj) + om * ( zrlxint - ptab(ji,jj) )
            END DO
         END DO

         ! lateral boundary condition on ptab
         CALL lbc_lnk( zrlx, 'T', 1. )

         ! convergence test
         zconv = 0.e0
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zconv = MAX( zconv, ABS( zrlx(ji,jj) - ptab(ji,jj) )  )
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_max( zconv )   ! max over the global domain

         ptab(:,:) = zrlx(:,:)

         !                                         !==========================
      END DO                                       ! end of sub-time step loop
      !                                            !==========================

      IF(ln_ctl)   THEN
         zrlx(:,:) = ptab(:,:) - ptab0(:,:)
         WRITE(charout,FMT="(' lim_hdf  : zconv =',D23.16, ' iter =',I4,2X)") zconv, iter
         CALL prt_ctl(tab2d_1=zrlx, clinfo1=charout)
      ENDIF

   END SUBROUTINE lim_hdf

#else
   !!----------------------------------------------------------------------
   !!   Default option          Dummy module           NO LIM sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_hdf         ! Empty routine
   END SUBROUTINE lim_hdf
#endif

   !!======================================================================
END MODULE limhdf
