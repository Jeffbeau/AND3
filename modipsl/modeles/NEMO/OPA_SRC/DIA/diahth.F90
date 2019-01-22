MODULE diahth
   !!======================================================================
   !!                       ***  MODULE  diahth  ***
   !! Ocean diagnostics: thermocline and 20 degree depth
   !!======================================================================
#if   defined key_diahth   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_diahth' :                              thermocline depth diag.
   !!----------------------------------------------------------------------
   !!   dia_hth      : Compute diagnostics associated with the thermocline
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dia_hth    ! routine called by step.F90

   !! * Shared module variables
   LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .TRUE.   !: thermocline-20d depths flag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      hth  ,      &  !: depth of the max vertical temperature gradient (m)
      hd20 ,      &  !: depth of 20 C isotherm (m)
      hd28 ,      &  !: depth of 28 C isotherm (m)
      htc3           !: heat content of first 300 m

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diahth.F90,v 1.3 2005/03/27 18:34:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dia_hth( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_hth  ***
      !!
      !! ** Purpose :
      !!      Computes the depth of strongest vertical temperature gradient
      !!      Computes the depth of the 20 degree isotherm
      !!      Computes the depth of the 28 degree isotherm
      !!      Computes the heat content of first 300 m
      !!
      !! ** Method : 
      !!
      !! History :
      !!        !  94-09  (J.-P. Boulanger)  Original code
      !!        !  96-11  (E. Guilyardi)  OPA8 
      !!        !  97-08  (G. Madec)  optimization
      !!        !  99-07  (E. Guilyardi)  hd28 + heat content 
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!-------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !! * Local declarations
      INTEGER :: ji, jj, jk         ! dummy loop arguments
      INTEGER :: iid, iif, ilevel   ! temporary integers
      INTEGER, DIMENSION(jpi) ::   idepth
      INTEGER, DIMENSION(jpi,jpj) ::   ikc

      REAL(wp) :: zd, zmoy              ! temporary scalars
      REAL(wp), DIMENSION(jpi) ::   zmax
      REAL(wp), DIMENSION(jpi,jpk) ::   zdzt
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dia_hth : diagnostics of the thermocline depth'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         IF(lwp) WRITE(numout,*)
      ENDIF


      ! -------------------------- !
      !  Depth of the thermocline  !
      ! -------------------------- !
      ! The depth of the thermocline is defined as the depth of the 
      ! strongest vertical temperature gradient
      
      DO jj = 1, jpj
         
         ! vertical gradient of temperature
         DO jk = 2, jpkm1
            zdzt(:,jk) = ( tn(:,jj,jk-1) - tn(:,jj,jk) ) / fse3w(:,jj,jk) * tmask(:,jj,jk)
         END DO
         
         ! search the level of maximum vertical temperature gradient
         zmax  (:) = 0.e0
         idepth(:) = 1
         DO jk = jpkm1, 2, -1
            DO ji = 1, jpi
               IF( zdzt(ji,jk) > zmax(ji) ) THEN
                  zmax  (ji) = zdzt(ji,jk)
                  idepth(ji) = jk
               ENDIF
            END DO
         END DO

         ! depth of the thermocline
         DO ji = 1, jpi
            hth(ji,jj) = fsdepw(ji,jj,idepth(ji))
         END DO
         
      END DO


      ! ----------------------- !
      !  Depth of 20C isotherm  !
      ! ----------------------- !

      ! initialization to the number of ocean w-point mbathy
      ! (cf dommsk, minimum value: 1)
      ikc(:,:) = 1

      ! search the depth of 20 degrees isotherm
      ! ( starting from the top, last level above 20C, if not exist, = 1)
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tn(ji,jj,jk) >= 20. ) ikc(ji,jj) = jk
            END DO
         END DO
      END DO
      
      ! Depth of 20C isotherm
      DO jj = 1, jpj
         DO ji = 1, jpi
            iid = ikc(ji,jj)
            iif = mbathy(ji,jj)
            IF( iid /= 1 ) THEN 
               ! linear interpolation
               zd =  fsdept(ji,jj,iid)   &
                  + (    fsdept(ji,jj,iid+1) - fsdept(ji,jj,iid) )   &
                  * ( 20.*tmask(ji,jj,iid+1) -     tn(ji,jj,iid) )   &
                  / (        tn(ji,jj,iid+1) -     tn(ji,jj,iid)    &
                  + (1.-tmask(ji,jj,1))                       )
               ! bound by the ocean depth, minimum value, first T-point depth
               hd20(ji,jj) = MIN( zd*tmask(ji,jj,1), fsdepw(ji,jj,iif))
            ELSE 
               hd20(ji,jj)=0.
            ENDIF
         END DO
      END DO

      ! ----------------------- !
      !  Depth of 28C isotherm  ! 
      ! ----------------------- !
      
      ! initialization to the number of ocean w-point mbathy
      ! (cf dommsk, minimum value: 1)
      ikc(:,:) = 1
      
      ! search the depth of 28 degrees isotherm
      ! ( starting from the top, last level above 28C, if not exist, = 1)
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( tn(ji,jj,jk) >= 28. ) ikc(ji,jj) = jk
            END DO
         END DO
      END DO
      
      ! Depth of 28C isotherm
      DO jj = 1, jpj
         DO ji = 1, jpi
            iid = ikc(ji,jj)
            iif = mbathy(ji,jj)
            IF( iid /= 1 ) THEN 
               ! linear interpolation
               zd =  fsdept(ji,jj,iid)   &
                  + (    fsdept(ji,jj,iid+1) - fsdept(ji,jj,iid) )   &
                  * ( 28.*tmask(ji,jj,iid+1) -     tn(ji,jj,iid) )   &
                  / (        tn(ji,jj,iid+1) -     tn(ji,jj,iid)    &
                  + ( 1. - tmask(ji,jj,1) )  )
               ! bound by the ocean depth, minimum value, first T-point depth
               hd28(ji,jj) = MIN( zd*tmask(ji,jj,1), fsdepw(ji,jj,iif) )
            ELSE 
               hd28(ji,jj) = 0.
            ENDIF
         END DO
      END DO

      ! ----------------------------------------- !
      !  Heat content of first 300 m (18 levels)  !
      ! ----------------------------------------- !

      htc3(:,:) = 0.e0
      ilevel = 18
      zmoy = rau0 * rcp * 0.5
      
      ! intregrate tn from surface to klevel

      DO jk = 1, ilevel
               htc3(:,:) = htc3(:,:)   &
                         + zmoy * ( tn(:,:,jk) + tn(:,:,jk+1) ) * fse3w(:,:,jk) * tmask(:,:,jk)
      END DO

   END SUBROUTINE dia_hth

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   USE in_out_manager
   LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .FALSE.  !: thermocline-20d depths flag
CONTAINS
   SUBROUTINE dia_hth( kt )         ! Empty routine
      if(lwp) WRITE(numout,*) 'dia_hth: You should not have seen this print! error?', kt
   END SUBROUTINE dia_hth
#endif

   !!======================================================================
END MODULE diahth
