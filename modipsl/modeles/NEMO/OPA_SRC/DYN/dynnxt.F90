MODULE dynnxt
   !!======================================================================
   !!                       ***  MODULE  dynnxt  ***
   !! Ocean dynamics: time stepping
   !!======================================================================
   
   !!----------------------------------------------------------------------
   !!   dyn_nxt      : update the horizontal velocity from the momentum trend
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE obc_oce         ! ocean open boundary conditions
   USE obcdyn          ! open boundary condition for momentum (obc_dyn routine)
   USE obcdyn_bt       ! 2D open boundary condition for momentum (obc_dyn_bt routine)
   USE obcvol          ! ocean open boundary condition (obc_vol routines)
   USE dynspg_oce      ! type of surface pressure gradient
   USE lbclnk          ! lateral boundary condition (or mpp link)
   USE prtctl          ! Print control
   USE agrif_opa_update
   USE agrif_opa_interp
#ifdef key_RIVER_INPUT
   USE rivers !JC:
#endif

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_nxt                ! routine called by step.F90
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_nxt ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_nxt  ***
      !!                   
      !! ** Purpose :   Compute the after horizontal velocity from the 
      !!      momentum trend.
      !!
      !! ** Method  :   Apply lateral boundary conditions on the trends (ua,va) 
      !!      through calls to routine lbc_lnk.
      !!      After velocity is compute using a leap-frog scheme environment:
      !!         (ua,va) = (ub,vb) + 2 rdt (ua,va)
      !!      Note that if lk_dynspg_flt=T, the time stepping has already been
      !!      performed in dynspg module
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
      !!        !  02-10  (C. Talandier, A-M. Treguier) Open boundary cond.
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   z2dt         ! temporary scalar
      !!----------------------------------------------------------------------
      !!  OPA 9.0 , LOCEAN-IPSL (2005) 
      !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynnxt.F90,v 1.10 2006/03/10 10:55:41 opalod Exp $ 
      !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_nxt : time stepping'
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
#if defined key_dynspg_flt
         ! Leap-frog time stepping already done in dynspg.F routine
#else
         DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
            DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
               ! Leap-frog time stepping
               ua(ji,jj,jk) = ( ub(ji,jj,jk) + z2dt * ua(ji,jj,jk) ) * umask(ji,jj,jk)
               va(ji,jj,jk) = ( vb(ji,jj,jk) + z2dt * va(ji,jj,jk) ) * vmask(ji,jj,jk)
            END DO
         END DO
#ifdef key_RIVER_INPUT
!JC:
   call riv_dyna( kt ) !JC: ua and va now need to be updated before swapping into un/vn 
#endif
# if defined key_obc
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      ! Update (ua,va) along open boundaries (only in the rigid-lid case)
      CALL obc_dyn( kt )

      IF ( lk_dynspg_exp .OR. lk_dynspg_ts ) THEN
         !Flather boundary condition :
         !        - Update sea surface height on each open boundary
         !                 sshn (= after ssh) for explicit case
         !                 sshn_b (= after ssha_b) for time-splitting case
         !        - Correct the barotropic velocities
         CALL obc_dyn_bt( kt )

         !Boundary conditions on sshn ( after ssh)
         CALL lbc_lnk( sshn, 'T', 1. )

         IF(ln_ctl) THEN         ! print sum trends (used for debugging)
            CALL prt_ctl(tab2d_1=sshn, clinfo1=' ssh      : ', mask1=tmask)
         ENDIF

         IF ( ln_vol_cst ) CALL obc_vol( kt )

      ENDIF

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
# endif
# if defined key_agrif
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      ! Update (ua,va) along open boundaries (only in the rigid-lid case)
      CALL Agrif_dyn( kt )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
# endif
#endif
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
         CALL prt_ctl(tab3d_1=un, clinfo1=' nxt  - Un: ', mask1=umask, &
            &         tab3d_2=vn, clinfo2=' Vn: ', mask2=vmask)
      ENDIF

#if defined key_agrif
      IF (.NOT.Agrif_Root())    CALL Agrif_Update_Dyn( kt )
#endif      

   END SUBROUTINE dyn_nxt

   !!======================================================================
END MODULE dynnxt
