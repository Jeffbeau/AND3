MODULE dynnxt1d
   !!======================================================================
   !!                       ***  MODULE  dynnxt1d  ***
   !! Ocean dynamics: time stepping in 1D configuration
   !!======================================================================
#if defined key_cfg_1d
   !!----------------------------------------------------------------------
   !!   'key_cfg_1d'               1D Configuration
   !!----------------------------------------------------------------------  
   !!----------------------------------------------------------------------
   !!   dyn_nxt_1d   : update the horizontal velocity from the momentum trend
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! lateral boundary condition (or mpp link)
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_nxt_1d                ! routine called by step.F90
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_nxt_1d ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_nxt_1d  ***
      !!                   
      !! ** Purpose :   Compute the after horizontal velocity from the 
      !!      momentum trend.
      !!
      !! ** Method  :   Apply lateral boundary conditions on the trends (ua,va) 
      !!      through calls to routine lbc_lnk.
      !!      After velocity is compute using a leap-frog scheme environment:
      !!         (ua,va) = (ub,vb) + 2 rdt (ua,va)
      !!      Time filter applied on now horizontal velocity to avoid the
      !!      divergence of two consecutive time-steps and swap of dynamics
      !!      arrays to start the next time step:
      !!         (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
      !!         (un,vn) = (ua,va) 
      !!
      !! ** Action : - Update ub,vb arrays, the before horizontal velocity
      !!             - Update un,vn arrays, the now horizontal velocity
      !!
      !! History :
      !!        !  87-02  (P. Andrich, D. L Hostis)  Original code
      !!        !  90-10  (C. Levy, G. Madec)
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  97-02  (G. Madec & M. Imbard)  opa, release 8.0
      !!        !  97-04  (A. Weaver)  Euler forward step
      !!        !  97-06  (G. Madec)  lateral boudary cond., lbc routine
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!        !  04-10  (C. Ethe) 1D configuration
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   z2dt         ! temporary scalar
      !!----------------------------------------------------------------------
      !!   OPA 9.0 , LOCEAN-IPSL (2005) 
      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/dynnxt1d.F90,v 1.3 2005/10/03 09:20:35 opalod Exp $ 
      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_nxt_1d : time stepping on 1D configuation'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! Local constant initialization
      z2dt = 2. * rdt
      IF( neuler == 0 .AND. kt == nit000 )  z2dt = rdt

      ! Lateral boundary conditions on ( ua, va )
      CALL lbc_lnk( ua, 'U', -1. )
      CALL lbc_lnk( va, 'V', -1. )

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Next velocity
         ! -------------
         DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
            DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
               ! Leap-frog time stepping
               ua(ji,jj,jk) = ( ub(ji,jj,jk) + z2dt * ua(ji,jj,jk) ) * umask(ji,jj,jk)
               va(ji,jj,jk) = ( vb(ji,jj,jk) + z2dt * va(ji,jj,jk) ) * vmask(ji,jj,jk)
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
 
     !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Time filter and swap of dynamics arrays
         ! ------------------------------------------
         IF( neuler == 0 .AND. kt == nit000 ) THEN
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  ! Euler (forward) time stepping
                  ub(ji,jj,jk) = un(ji,jj,jk)
                  vb(ji,jj,jk) = vn(ji,jj,jk)
                  un(ji,jj,jk) = ua(ji,jj,jk)
                  vn(ji,jj,jk) = va(ji,jj,jk)
               END DO
            END DO
         ELSE
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  ! Leap-frog time stepping
                  ub(ji,jj,jk) = atfp * ( ub(ji,jj,jk) + ua(ji,jj,jk) ) + atfp1 * un(ji,jj,jk)
                  vb(ji,jj,jk) = atfp * ( vb(ji,jj,jk) + va(ji,jj,jk) ) + atfp1 * vn(ji,jj,jk)
                  un(ji,jj,jk) = ua(ji,jj,jk)
                  vn(ji,jj,jk) = va(ji,jj,jk)
               END DO
            END DO
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      IF(ln_ctl)   THEN
         CALL prt_ctl(tab3d_1=un, clinfo1=' nxt_1d  - Un: ', mask1=umask, &
            &         tab3d_2=vn, clinfo2=' Vn: ', mask2=vmask)
      ENDIF

!     IF(l_ctl)   WRITE(numout,*) ' nxt  - Un: ', SUM(un(2:nictl,2:njctl,1:jpkm1)*umask(2:nictl,2:njctl,1:jpkm1)), &
!     &                                  ' Vn: ', SUM(vn(2:nictl,2:njctl,1:jpkm1)*vmask(2:nictl,2:njctl,1:jpkm1))

   END SUBROUTINE dyn_nxt_1d
#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO 1D Config
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE dyn_nxt_1d ( kt )
      if(lwp) WRITE(numout,*) 'dyn_nxt_1d: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_nxt_1d
#endif
   !!======================================================================
END MODULE dynnxt1d
