MODULE tranxt
   !!======================================================================
   !!                       ***  MODULE  tranxt  ***
   !! Ocean active tracers:  time stepping on temperature and salinity
   !!======================================================================
   
   !!----------------------------------------------------------------------
   !!   tra_nxt     : time stepping on temperature and salinity
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ???
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE obctra          ! open boundary condition (obc_tra routine)
   USE prtctl          ! Print control
   USE agrif_opa_update
   USE agrif_opa_interp

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_nxt          ! routine called by step.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tranxt.F90,v 1.7 2006/03/10 10:55:45 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_nxt( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tranxt  ***
      !!
      !! ** Purpose :   Compute the temperature and salinity fields at the 
      !!      next time-step from their temporal trends and swap the fields.
      !! 
      !! ** Method  :   Apply lateral boundary conditions on (ta,ta) through 
      !!      call to lbc_lnk routine
      !!      After t and s are compute using a leap-frog scheme environment:
      !!         ta = tb + 2 rdttra(k) * ta
      !!         sa = sb + 2 rdttra(k) * sa
      !!      Compute and save in (ta,sa) an average over three time levels
      !!      (before,now and after) of temperature and salinity which is
      !!      used to compute rhd in eos routine and thus the hydrostatic 
      !!      pressure gradient (ln_dynhpg_imp = T)
      !!      Apply an Asselin time filter on now tracers (tn,sn) to avoid
      !!      the divergence of two consecutive time-steps and swap tracer
      !!      arrays to prepare the next time_step:
      !!         (zt,zs) = (ta+2tn+tb,sa+2sn+sb)/4       (ln_dynhpg_imp = T)
      !!         (zt,zs) = (0,0)                            (default option)
      !!         (tb,sb) = (tn,vn) + atfp [ (tb,sb) + (ta,sa) - 2 (tn,sn) ]
      !!         (tn,sn) = (ta,sa) 
      !!         (ta,sa) = (zt,zs)  (NB: reset to 0 after use in eos.F)
      !!
      !! ** Action  : - update (tb,sb) and (tn,sn) 
      !!              - (ta,sa) time averaged (t,s)      (ln_dynhpg_imp = T)
      !!
      !! History :
      !!   7.0  !  91-11  (G. Madec)  Original code
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  96-02  (G. Madec & M. Imbard)  opa release 8.0
      !!   8.0  !  96-04  (A. Weaver)  Euler forward step
      !!   8.2  !  99-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!        !  02-11  (C. Talandier, A-M Treguier) Open boundaries
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index

      !! * Local declarations
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zt, zs      ! temporary scalars
      REAL(wp) ::   zfact       ! temporary scalar
      !!----------------------------------------------------------------------


      ! 0. Lateral boundary conditions on ( ta, sa )   (T-point, unchanged sign)
      ! ---------------------------------============
      CALL lbc_lnk( ta, 'T', 1. )   
      CALL lbc_lnk( sa, 'T', 1. )


      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! 1. Leap-frog scheme (only in explicit case, otherwise the 
         ! -------------------  time stepping is already done in trazdf)
         IF( l_trazdf_exp ) THEN
            zfact = 2. * rdttra(jk)
            IF( neuler == 0 .AND. kt == nit000 ) zfact = rdttra(jk)
            ta(:,:,jk) = ( tb(:,:,jk) + zfact * ta(:,:,jk) ) * tmask(:,:,jk)
            sa(:,:,jk) = ( sb(:,:,jk) + zfact * sa(:,:,jk) ) * tmask(:,:,jk)
!           sa(:,:,jk) = max(sa(:,:,jk),5d0)
         ENDIF

         DO ji=1,jpi
         DO jj=1,jpj

            IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0 !BIO
         END DO
         END DO


#if defined key_obc
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! Update tracers on open boundaries.
      CALL obc_tra( kt )

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
#endif
#if defined key_agrif
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! Update tracers on open boundaries.
      CALL Agrif_tra( kt )

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
#endif


         ! 2. Time filter and swap of arrays
         ! ---------------------------------
 
         IF( ln_dynhpg_imp ) THEN                       ! semi-implicite hpg
            IF( neuler == 0 .AND. kt == nit000 ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zt = ( ta(ji,jj,jk) + 2. * tn(ji,jj,jk) + tb(ji,jj,jk) ) * 0.25
                     zs = ( sa(ji,jj,jk) + 2. * sn(ji,jj,jk) + sb(ji,jj,jk) ) * 0.25
                     tb(ji,jj,jk) = tn(ji,jj,jk)
                     sb(ji,jj,jk) = sn(ji,jj,jk)
                     tn(ji,jj,jk) = ta(ji,jj,jk)
                     IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
                     sn(ji,jj,jk) = sa(ji,jj,jk)
!                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
                     ta(ji,jj,jk) = zt
                     IF(zs.lt.0.0) zs=0.0 !BIO
                     sa(ji,jj,jk) = zs
!                    sa(ji,jj,jk) = max(zs,5d0)
                  END DO
               END DO
            ELSE
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zt = ( ta(ji,jj,jk) + 2. * tn(ji,jj,jk) + tb(ji,jj,jk) ) * 0.25
                     zs = ( sa(ji,jj,jk) + 2. * sn(ji,jj,jk) + sb(ji,jj,jk) ) * 0.25
                     tb(ji,jj,jk) = atfp  * ( tb(ji,jj,jk) + ta(ji,jj,jk) ) + atfp1 * tn(ji,jj,jk)
                     sb(ji,jj,jk) = atfp  * ( sb(ji,jj,jk) + sa(ji,jj,jk) ) + atfp1 * sn(ji,jj,jk)
                     tn(ji,jj,jk) = ta(ji,jj,jk)
                    IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
                     sn(ji,jj,jk) = sa(ji,jj,jk)
!                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
                     ta(ji,jj,jk) = zt
                     IF(zs.lt.0.0) zs=0.0 !BIO
                     sa(ji,jj,jk) = zs
!                    sa(ji,jj,jk) = max(zs,5d0)
                  END DO
               END DO
            ENDIF
         ELSE                                          ! Default case
            IF( neuler == 0 .AND. kt == nit000 ) THEN
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     tb(ji,jj,jk) = tn(ji,jj,jk)
                     sb(ji,jj,jk) = sn(ji,jj,jk)
                     tn(ji,jj,jk) = ta(ji,jj,jk)
                     IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
                     sn(ji,jj,jk) = sa(ji,jj,jk)
!                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
                  END DO
               END DO
            ELSE
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     tb(ji,jj,jk) = atfp  * ( tb(ji,jj,jk) + ta(ji,jj,jk) ) + atfp1 * tn(ji,jj,jk)
                     sb(ji,jj,jk) = atfp  * ( sb(ji,jj,jk) + sa(ji,jj,jk) ) + atfp1 * sn(ji,jj,jk)
                     tn(ji,jj,jk) = ta(ji,jj,jk)
                     IF(sa(ji,jj,jk).lt.0.0) sa(ji,jj,jk)=0.0  !BIO
                     sn(ji,jj,jk) = sa(ji,jj,jk)
!                    sn(ji,jj,jk) = max(sa(ji,jj,jk),5d0)
                  END DO
               END DO
            ENDIF
         ENDIF
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      IF(ln_ctl) THEN         ! print mean field (used for debugging)
         CALL prt_ctl(tab3d_1=tn, clinfo1=' nxt  - Tn: ', mask1=tmask, &
            &         tab3d_2=sn, clinfo2=' Sn: ', mask2=tmask)
      ENDIF
      
#if defined key_agrif
      IF (.NOT.Agrif_Root())    CALL Agrif_Update_Tra( kt )
#endif      

   END SUBROUTINE tra_nxt

   !!======================================================================
END MODULE tranxt
