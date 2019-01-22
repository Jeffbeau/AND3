MODULE dynhpg_atsk
   !!======================================================================
   !!                       ***  MODULE  dynhpg_atsk  ***
   !! Ocean dynamics:  hydrostatic pressure gradient trend
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   dyn_hpg_atsk : update the momentum trend with the horizontal
   !!                  gradient of the hydrostatic pressure
   !!
   !!   default case : use of 3D work arrays (vector opt. available)
   !!   key_s_coord       : s-coordinate
   !!   key_partial_steps : z-coordinate with partial steps
   !!   default key       : z-coordinate
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_hpg_atsk ! routine called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynhpg_atsk.F90,v 1.9 2005/09/02 15:45:23 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

#if defined key_s_coord
   !!---------------------------------------------------------------------
   !!                  ***  dynhpg_atsk.h90  ***
   !!---------------------------------------------------------------------
   !!   'key_s_coord'                                         s-coordinate
   !!---------------------------------------------------------------------

   SUBROUTINE dyn_hpg_atsk( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg_atsk  ***
      !!        
      !! ** Purpose :   Compute the now momentum trend due to the horizontal 
      !!     gradient of the hydrostatic pressure. Add it to the general
      !!     momentum trend.
      !!
      !! ** Method  :   The now hydrostatic pressure gradient at a given level
      !!      jk is computed by taking the vertical integral of the in-situ 
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinate case ('key_s_coord'): a corrective term is
      !!      added to the horizontal pressure gradient :
      !!         zhpi = grav .....   + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....   + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!      j-k-i loop (j-slab) ('key_autotasking')
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
      !!
      !! History :
      !!   7.0  !  96-01  (G. Madec)  s-coordinates
      !!        !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * modules used
      USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
         &              zhpj => sa      ! use sa as 3D workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step index
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   &
         zcoef0, zcoef1, zuap, zvap     ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztdua, ztdva                   ! temporary scalars
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_hpg_atsk : s-coordinate hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   autotasking case (j-k-i loop)'
      ENDIF

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      zcoef0 = - grav * 0.5
      zuap   = 0.e0
      zvap   = 0.e0

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============
         ! 1. Surface value
         ! ----------------
         DO ji = 2, jpim1
            ! hydrostatic pressure gradient along s-surfaces
            zhpi(ji,jj,1) = zcoef0 / e1u(ji,jj)   &
                       * ( fse3w(ji+1,jj,1) * rhd(ji+1,jj,1) - fse3w(ji,jj,1) * rhd(ji,jj,1)  )
            zhpj(ji,jj,1) = zcoef0 / e2v(ji,jj)   &
                       * ( fse3w(ji,jj+1,1) * rhd(ji,jj+1,1) - fse3w(ji,jj,1) * rhd(ji,jj,1)  )
            ! s-coordinate pressure gradient correction
            zuap = -zcoef0 * ( rhd(ji+1,jj,1) + rhd(ji,jj,1) )   &
                 * ( fsde3w(ji+1,jj,1) - fsde3w(ji,jj,1) ) / e1u(ji,jj)
            zvap = -zcoef0 * ( rhd(ji,jj+1,1) + rhd(ji,jj,1) )   &
                 * ( fsde3w(ji,jj+1,1) - fsde3w(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1) + zuap
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1) + zvap
         END DO  

         ! 2. interior value (2=<jk=<jpkm1)
         ! -----------------
         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               ! hydrostatic pressure gradient along s-surfaces
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 / e1u(ji,jj)   &
                  &           * ( fse3w(ji+1,jj,jk) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) )   &
                  &              -fse3w(ji  ,jj,jk) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) )  )
               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 / e2v(ji,jj)   &
                  &           * ( fse3w(ji,jj+1,jk) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) )   &
                  &              -fse3w(ji,jj  ,jk) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) )  )
               ! s-coordinate pressure gradient correction 
               zuap = -zcoef0 * ( rhd(ji+1,jj  ,jk) + rhd(ji,jj,jk) )   &
                    * ( fsde3w(ji+1,jj,jk) - fsde3w(ji,jj,jk) ) / e1u(ji,jj)
               zvap = -zcoef0 * ( rhd(ji  ,jj+1,jk) + rhd(ji,jj,jk) )   &
                    * ( fsde3w(ji,jj+1,jk) - fsde3w(ji,jj,jk) ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk) + zuap
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk) + zvap
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the hydrostatic pressure gradient trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)
         CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_hpg_atsk

#elif defined key_partial_steps
   !!---------------------------------------------------------------------
   !!   'key_partial_steps'                     z-coordinate partial steps
   !!---------------------------------------------------------------------

   SUBROUTINE dyn_hpg_atsk( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg_atsk  ***
      !!  
      !! ** Purpose :   Compute the now momentum trend due to the hor. gradient
      !!      of the hydrostatic pressure. Add it to the general momentum trend.
      !!
      !! ** Method  :   The now hydrostatic pressure gradient at a given level
      !!      jk is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that 
      !!      level:    zhpi = grav .....
      !!                zhpj = grav .....
      !!      add it to the general momentum trend (ua,va).
      !!            ua = ua - 1/e1u * zhpi
      !!            va = va - 1/e2v * zhpj
      !!      j-k-i loop (j-slab) ('key_autotasking')
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
      !!
      !! History :
      !!   8.5  !  02-08  (A. Bozec)  Original code
      !!----------------------------------------------------------------------
      !! * modules used
      USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
         &              zhpj => sa      ! use sa as 3D workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step index

      !! * local declarations
      INTEGER ::   ji, jj, jk           ! dummy loop indices
      INTEGER ::   iku, ikv             ! temporary integers
      REAL(wp) ::   &
         zcoef0, zcoef1, zuap,       &  ! temporary scalars
         zcoef2, zcoef3, zvap           !    "         "
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztdua, ztdva                   ! temporary scalars
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_hpg_atsk : z-coord. partial steps hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   autotasking case (j-k-i loop)'
      ENDIF

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      zcoef0 = - grav * 0.5
      zuap   = 0.e0
      zvap   = 0.e0
      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============
         ! 1. Surface value
         ! ----------------
         DO ji = 2, jpim1
            zcoef1 = zcoef0 * fse3w(ji,jj,1)
            ! hydrostatic pressure gradient
            zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
            zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
         END DO

         ! 2. interior value (2=<jk=<jpkm1)
         ! -----------------
         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               zcoef1 = zcoef0 * fse3w(ji,jj,jk)
               ! hydrostatic pressure gradient
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
                  &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)

               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
                  &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO 
         END DO

         ! partial steps correction at the last level  (new gradient with  intgrd.F)
         DO ji = 2, jpim1
            iku = MIN ( mbathy(ji,jj), mbathy(ji+1,jj) ) - 1
            ikv = MIN ( mbathy(ji,jj), mbathy(ji,jj+1) ) - 1
            zcoef2 = zcoef0 * MIN( fse3w(ji,jj,iku), fse3w(ji+1,jj  ,iku) )
            zcoef3 = zcoef0 * MIN( fse3w(ji,jj,ikv), fse3w(ji  ,jj+1,ikv) )
            ! on i-direction
            IF ( iku > 2 ) THEN
               ! subtract old value  
               ua(ji,jj,iku) = ua(ji,jj,iku) - zhpi(ji,jj,iku)
               ! compute the new one   
               zhpi (ji,jj,iku) = zhpi(ji,jj,iku-1)   &
                  + zcoef2 * ( rhd(ji+1,jj,iku-1) - rhd(ji,jj,iku-1) + gru(ji,jj) ) / e1u(ji,jj)
               ! add the new one to the general momentum trend
               ua(ji,jj,iku) = ua(ji,jj,iku) + zhpi(ji,jj,iku)
            ENDIF
            ! on j-direction
            IF ( ikv > 2 ) THEN
               ! subtract old value  
               va(ji,jj,ikv) = va(ji,jj,ikv) - zhpj(ji,jj,ikv)
               ! compute the new one   
               zhpj (ji,jj,ikv) = zhpj(ji,jj,ikv-1)   &
                  + zcoef3 * ( rhd(ji,jj+1,ikv-1) - rhd(ji,jj,ikv-1) + grv(ji,jj) ) / e2v(ji,jj)
               ! add the new one to the general momentum trend
               va(ji,jj,ikv) = va(ji,jj,ikv) + zhpj(ji,jj,ikv)
            ENDIF
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the hydrostatic pressure gradient trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)
         CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF   

   END SUBROUTINE dyn_hpg_atsk

#else
   !!---------------------------------------------------------------------
   !!   Default case :                                        z-coordinate
   !!---------------------------------------------------------------------

   SUBROUTINE dyn_hpg_atsk( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg_atsk  ***
      !!   
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the hydrostatic pressure. Add it to the general
      !!      momentum trend.
      !!
      !! ** Method  :   The now hydrostatic pressure gradient at a given level
      !!      jk is computed by taking the vertical integral of the in-situ  
      !!      density gradient along the model level from the suface to that
      !!      level:    zhpi = grav .....
      !!                zhpj = grav .....
      !!      add it to the general momentum trend (ua,va).
      !!            ua = ua - 1/e1u * zhpi
      !!            va = va - 1/e2v * zhpj
      !!      j-k-i loop (j-slab) ('key_autotasking')
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend in (utrd,vtrd) ('key_trddyn')
      !!
      !! History :
      !!   1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  96-01  (G. Madec)  s-coordinates
      !!        !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
      !!   8.5  !  02-07  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * modules used
      USE oce, ONLY :   zhpi => ta,  &  ! use ta as 3D workspace
         &              zhpj => sa      ! use sa as 3D workspace

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step index

      !! * local declarations
      INTEGER ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   &
         zcoef0, zcoef1, zuap, zvap     ! temporary scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
         ztdua, ztdva                   ! temporary scalars
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_hpg_atsk : z-coordinate hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   auto-tasking (j-k-i loop)'
      ENDIF

      ! Save ua and va trends
      IF( l_trddyn )   THEN
         ztdua(:,:,:) = ua(:,:,:) 
         ztdva(:,:,:) = va(:,:,:) 
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      zcoef0 = - grav * 0.5
      zuap   = 0.e0
      zvap   = 0.e0

      !                                                ! ===============
      DO jj = 2, jpjm1                                 !  Vertical slab
         !                                             ! ===============
         ! 1. Surface value
         ! ----------------
         

         DO ji = 2, jpim1
            zcoef1 = zcoef0 * fse3w(ji,jj,1)
            ! hydrostatic pressure gradient
            zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) / e1u(ji,jj)
            zhpj(ji,jj,1) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua(ji,jj,1) = ua(ji,jj,1) + zhpi(ji,jj,1)
            va(ji,jj,1) = va(ji,jj,1) + zhpj(ji,jj,1)
         END DO

         ! 2. interior value (2=<jk=<jpkm1)
         ! -----------------
         DO jk = 2, jpkm1
            DO ji = 2, jpim1
               zcoef1 = zcoef0 * fse3w(ji,jj,jk)
               ! hydrostatic pressure gradient
               zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )   &
                  &                       - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) / e1u(ji,jj)

               zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )   &
                  &                       - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua(ji,jj,jk) = ua(ji,jj,jk) + zhpi(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) + zhpj(ji,jj,jk)
            END DO 
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! save the hydrostatic pressure gradient trends for diagnostic
      ! momentum trends
      IF( l_trddyn )   THEN
         zhpi(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
         zhpj(:,:,:) = va(:,:,:) - ztdva(:,:,:)

         CALL trd_mod(zhpi, zhpj, jpdtdhpg, 'DYN', kt)
      ENDIF

      IF(ln_ctl) THEN         ! print sum trends (used for debugging)
         CALL prt_ctl(tab3d_1=ua, clinfo1=' hpg  - Ua: ', mask1=umask, &
            &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
      ENDIF

   END SUBROUTINE dyn_hpg_atsk

#endif

   !!======================================================================
END MODULE dynhpg_atsk
