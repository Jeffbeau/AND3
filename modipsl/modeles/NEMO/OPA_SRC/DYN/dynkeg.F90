MODULE dynkeg
   !!======================================================================
   !!                       ***  MODULE  dynkeg  ***
   !! Ocean dynamics:  kinetic energy gradient trend
   !!======================================================================
   
   !!----------------------------------------------------------------------
   !!   dyn_keg      : update the momentum trend with the horizontal tke
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_keg                ! routine called by step.F90
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynkeg.F90,v 1.7 2005/09/02 15:45:23 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_keg( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the 
      !!      general momentum trend.
      !!
      !! ** Method  :   Compute the now horizontal kinetic energy:
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend (ua,va).
      !!         ua = ua - 1/e1u di[ zhke ]
      !!         va = va - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (ua, va) with the hor. ke gradient trend
      !!             - Save the trends in (utrd,vtrd) ('key_trddyn')
      !!
      !! History :
      !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
      !!   7.0  !  97-05  (G. Madec)  Split dynber into dynkeg and dynhpg
      !!   9.0  !  02-07  (G. Madec)  F90: Free form and module
      !!    "   !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdua => ta,   & ! use ta as 3D workspace   
                         ztdva => sa      ! use sa as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step index

      !! * Local declarations
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      REAL(wp) ::   zua, zva, zu, zv    ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         zhke                           ! temporary workspace
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF
      
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Horizontal kinetic energy at T-point
         DO jj = 2, jpj
            DO ji = fs_2, jpi   ! vector opt.
               zv = 0.25 * (  vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
                            + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)  )
               zu = 0.25 * (  un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
                            + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)  )
               zhke(ji,jj,jk) = zv + zu
            END DO  
         END DO  
         
         ! Horizontal gradient of Horizontal kinetic energy
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               ! gradient of kinetic energy
               zua = -( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
               zva = -( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
               ! add to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO 
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the Kinetic Energy trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)

         CALL trd_mod(ztdua, ztdva, jpdtdkeg, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' keg  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_keg

   !!======================================================================
END MODULE dynkeg
