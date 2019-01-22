MODULE ldfdyn
   !!======================================================================
   !!                       ***  MODULE  ldfdyn  ***
   !! Ocean physics:  lateral viscosity coefficient 
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   ldf_dyn_init : initialization, namelist read, and parameters control
   !!   ldf_dyn_c3d   : 3D eddy viscosity coefficient initialization
   !!   ldf_dyn_c2d   : 2D eddy viscosity coefficient initialization
   !!   ldf_dyn_c1d   : 1D eddy viscosity coefficient initialization
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE ldfdyn_oce      ! ocean dynamics lateral physics
   USE phycst          ! physical constants
   USE ldfslp          ! ???
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC ldf_dyn_init   ! called by opa.F90
!!DB: due to new code changes must make the relevant routine PUBLIC or step() cannot find it
#if defined key_dynldf_smag || defined key_traldf_smag
   PUBLIC ldf_smag   ! possibly called by step.F90
#endif
#if defined key_dynldf_c3d
   PUBLIC ldf_dyn_c3d
#elif defined key_dynldf_c2d
   PUBLIC ldf_dyn_c2d
#elif defined key_dynldf_smag
   PUBLIC ldf_dyn_smag
#elif defined key_dynldf_c1d
   PUBLIC ldf_dyn_c1d
#else
!! do nothing
#endif




  INTERFACE ldf_zpf
     MODULE PROCEDURE ldf_zpf_1d, ldf_zpf_1d_3d, ldf_zpf_3d
  END INTERFACE

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn.F90,v 1.5 2005/03/27 18:35:06 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE ldf_dyn_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_init  ***
      !!                   
      !! ** Purpose :   set the horizontal ocean dynamics physics
      !!
      !! ** Method  :  
      !!      Eddy viscosity coefficients:
      !!         default option   : constant coef. ahm0 (namelist)
      !!        'key_dynldf_c1d': depth dependent coef. defined in 
      !!                        in ldf_dyn_c1d routine
      !!        'key_dynldf_c2d': latitude and longitude dependent coef.
      !!                        defined in ldf_dyn_c2d routine
      !!        'key_dynldf_c3d': latitude, longitude, depth dependent coef.
      !!                        defined in ldf_dyn_c3d routine
      !!      N.B. User defined include files.  By default, 3d and 2d coef.
      !!      are set to a constant value given in the namelist and the 1d
      !!      coefficients are initialized to a hyperbolic tangent vertical
      !!      profile.
      !!
      !! Reference :
      !!      Madec, G. and M. Imbard, 1996, A global ocean mesh to overcome
      !!      the North Pole singularity, Climate Dynamics, 12, 381-388.
      !!
      !! History :
      !!        !  07-97  (G. Madec)  from inimix.F split in 2 routines
      !!        !  08-97  (G. Madec)  multi dimensional coefficients
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Local declarations
      INTEGER ::   ioptio         ! ???
      LOGICAL :: ll_print = .FALSE.    ! Logical flag for printing viscosity coef.

       
      NAMELIST/nam_dynldf/ ln_dynldf_lap  , ln_dynldf_bilap,                &
         &                 ln_dynldf_level, ln_dynldf_hor, ln_dynldf_iso,   &
         &                 ahm0, ahmb0
      !!----------------------------------------------------------------------


      ! Define the lateral physics parameters
      ! ======================================
    
      ! Read Namelist nam_dynldf : Lateral physics
      REWIND( numnam )
      READ  ( numnam, nam_dynldf )

      ! Parameter print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_dyn : lateral momentum physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '          Namelist nam_dynldf : set lateral mixing parameters'
         WRITE(numout,*) '             laplacian operator          ln_dynldf_lap   = ', ln_dynldf_lap
         WRITE(numout,*) '             bilaplacian operator        ln_dynldf_bilap = ', ln_dynldf_bilap
         WRITE(numout,*) '             iso-level                   ln_dynldf_level = ', ln_dynldf_level
         WRITE(numout,*) '             horizontal (geopotential)   ln_dynldf_hor   = ', ln_dynldf_hor
         WRITE(numout,*) '             iso-neutral                 ln_dynldf_iso   = ', ln_dynldf_iso
         WRITE(numout,*) '             horizontal eddy viscosity            ahm0   = ', ahm0
         WRITE(numout,*) '             background viscosity                 ahmb0  = ', ahmb0
         WRITE(numout,*)
      ENDIF

      ! Parameter control

      ! control the input
      ioptio = 0
      IF( ln_dynldf_lap   )   ioptio = ioptio + 1
      IF( ln_dynldf_bilap )   ioptio = ioptio + 1
      IF( ioptio /= 1 )   THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) '          use ONE of the 2 lap/bilap operator type on momentum'
          nstop = nstop + 1
      ENDIF
      ioptio = 0
      IF( ln_dynldf_level )   ioptio = ioptio + 1
      IF( ln_dynldf_hor   )   ioptio = ioptio + 1
      IF( ln_dynldf_iso   )   ioptio = ioptio + 1
      IF( ioptio /= 1 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          use only ONE direction (level/hor/iso)'
         nstop = nstop + 1
      ENDIF

      IF( lk_sco ) THEN          ! s-coordinates: rotation required for horizontal or isopycnal direction
         IF( ( ln_dynldf_iso .OR. ln_dynldf_hor ) .AND. .NOT.lk_ldfslp ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          the rotation of the viscous tensor require key_ldfslp'
            IF( .NOT.lk_esopa )   nstop = nstop + 1
         ENDIF
      ELSE                       ! z-coordinates with/without partial step:
         ln_dynldf_level = ln_dynldf_level .OR. ln_dynldf_hor      ! level mixing = horizontal mixing
         ln_dynldf_hor   = .FALSE.
         IF(lwp) WRITE(numout,*) '          horizontal mixing in z-coord or partial steps: force ln_dynldf_level = T'
         IF(lwp) WRITE(numout,*) '                                                  and    force ln_dynldf_hor   = F'
         IF( ln_dynldf_iso .AND. .NOT.lk_ldfslp ) THEN             ! rotation required for isopycnal mixing
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          the rotation of the viscous tensor require key_ldfslp'
            IF( .NOT.lk_esopa )   nstop = nstop + 1
         ENDIF
      ENDIF

      l_dynldf_lap     =       ln_dynldf_lap   .AND. ln_dynldf_level     ! iso-level   laplacian operator
      l_dynldf_bilap   =       ln_dynldf_bilap .AND. ln_dynldf_level     ! iso-level bilaplacian operator
      l_dynldf_bilapg  =       ln_dynldf_bilap .AND. ln_dynldf_hor       ! geopotential bilap. (s-coord)
      l_dynldf_iso     =       ln_dynldf_lap   .AND.                  &  ! laplacian operator
         &                   ( ln_dynldf_iso   .OR.  ln_dynldf_hor )     ! iso-neutral (z-coord) or horizontal (s-coord)

      l_dynzdf_iso    = .FALSE.
      IF( l_dynldf_iso )   l_dynzdf_iso = .TRUE.

      ioptio = 0
      IF( l_dynldf_lap     )   ioptio = ioptio + 1
      IF( l_dynldf_bilap   )   ioptio = ioptio + 1
      IF( l_dynldf_bilapg  )   ioptio = ioptio + 1
      IF( l_dynldf_iso     )   ioptio = ioptio + 1
      IF( ioptio /= 1 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          this combination of operator and direction has not been implemented'
         nstop = nstop + 1
      ENDIF
      IF( lk_esopa ) THEN
         l_dynldf_lap = .TRUE.   ;   l_dynldf_bilap   = .TRUE.   ;   l_dynldf_bilapg  = .TRUE.
         l_dynldf_iso = .TRUE.   ;   l_dynzdf_iso     = .TRUE.
         IF(lwp ) WRITE(numout,*) '          esopa test: use all lateral physics options'
      ENDIF

      ! control print
      IF( l_dynldf_lap    .AND. lwp ) WRITE(numout,*) '          iso-level laplacian momentum operator'
      IF( l_dynldf_bilap  .AND. lwp ) WRITE(numout,*) '          iso-level bilaplacian momentum operator'
      IF( l_dynldf_bilapg .AND. lwp ) WRITE(numout,*) '          geopotential bilaplacian momentum operator'
      IF( l_dynldf_iso    .AND. lwp ) WRITE(numout,*) '          iso-neutral laplacian momentum operator'

      ! ... Space variation of eddy coefficients
      ioptio = 0
#if defined key_dynldf_c3d
      IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude, depth)'
      ioptio = ioptio+1
#endif
#if defined key_dynldf_smag
      IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude, depth)'
      ioptio = ioptio+1
#endif
#if defined key_dynldf_c2d
      IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( latitude, longitude)'
      ioptio = ioptio+1
#endif
#if defined key_dynldf_c1d
      IF(lwp) WRITE(numout,*) '          momentum mixing coef. = F( depth )'
      ioptio = ioptio+1
      IF( lk_sco ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          key_dynldf_c1d cannot be used in s-coordinate (key_s_coord)'
         nstop = nstop + 1
      ENDIF
#endif
      IF( ioptio == 0 ) THEN
          IF(lwp) WRITE(numout,*) '          momentum mixing coef. = constant  (default option)'
        ELSEIF( ioptio > 1 ) THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) '          use only one of the following keys:',   &
                                  ' key_dynldf_c3d, key_dynldf_c2d, key_dynldf_c1d'
          nstop = nstop + 1
      ENDIF


      IF( l_dynldf_bilap .OR. l_dynldf_bilapg ) THEN
         IF(lwp) WRITE(numout,*) '          biharmonic momentum diffusion'
         IF( ahm0 > 0 .AND. .NOT. lk_esopa ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) 'The horizontal viscosity coef. ahm0 must be negative'
            nstop = nstop + 1
         ENDIF
      ELSE
         IF(lwp) WRITE(numout,*) '          harmonic momentum diff. (default)'
         IF( ahm0 < 0 .AND. .NOT. lk_esopa ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          The horizontal viscosity coef. ahm0 must be positive'
            nstop = nstop + 1
         ENDIF
      ENDIF


      ! Lateral eddy viscosity
      ! ======================

#if defined key_dynldf_c3d
      CALL ldf_dyn_c3d( ll_print )   ! ahm = 3D coef. = F( longitude, latitude, depth )
#elif defined key_dynldf_c2d
      CALL ldf_dyn_c2d( ll_print )   ! ahm = 1D coef. = F( longitude, latitude )
#elif defined key_dynldf_smag
      CALL ldf_dyn_smag( ll_print )   ! ahm = 1D coef. = F( longitude, latitude )
#elif defined key_dynldf_c1d
      CALL ldf_dyn_c1d( ll_print )   ! ahm = 1D coef. = F( depth )
#else
                                     ! Constant coefficients
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'inildf: constant eddy viscosity coef. '
      IF(lwp) WRITE(numout,*) '~~~~~~'
      IF(lwp) WRITE(numout,*) '        ahm1 = ahm2 = ahm0 =  ',ahm0
#endif

   END SUBROUTINE ldf_dyn_init

#if defined key_dynldf_c3d
#   include "ldfdyn_c3d.h90"
#elif defined key_dynldf_c2d
#   include "ldfdyn_c2d.h90"
#elif defined key_dynldf_smag
#   include "ldfdyn_smag.h90"
#elif defined key_dynldf_c1d
#   include "ldfdyn_c1d.h90"
#endif

!!DB: 2008.04.09 -- changes to smag procedures
#if defined key_dynldf_smag || defined key_traldf_smag
#   include "ldfsmag.h90"
#endif 

   SUBROUTINE ldf_zpf_1d( ld_print, pdam, pwam, pbot, pdep, pah )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_zpf  ***
      !!                   
      !! ** Purpose :   vertical adimensional profile for eddy coefficient
      !!
      !! ** Method  :   1D eddy viscosity coefficients ( depth )
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT (in   ) ::   &
          pdam,     &  ! depth of the inflection point
          pwam,     &  ! width of inflection
          pbot         ! battom value (0<pbot<= 1)
      REAL(wp), INTENT (in   ), DIMENSION(jpk) ::   &
          pdep         ! depth of the gridpoint (T, U, V, F)
      REAL(wp), INTENT (inout), DIMENSION(jpk) ::   &
          pah          ! adimensional vertical profile

      !! * Local variables
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
      !!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )
      zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         pah(jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
      END DO

      ! Control print
      IF(lwp .AND. ld_print ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(jk), pdep(jk)
         END DO
      ENDIF

   END SUBROUTINE ldf_zpf_1d


   SUBROUTINE ldf_zpf_1d_3d( ld_print, pdam, pwam, pbot, pdep, pah )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_zpf  ***
      !!
      !! ** Purpose :   vertical adimensional profile for eddy coefficient
      !!
      !! ** Method  :   1D eddy viscosity coefficients ( depth )
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT (in   ) ::   &
          pdam,     &  ! depth of the inflection point
          pwam,     &  ! width of inflection
          pbot         ! battom value (0<pbot<= 1)
      REAL(wp), INTENT (in   ), DIMENSION(jpk) ::   &
          pdep         ! depth of the gridpoint (T, U, V, F)
      REAL(wp), INTENT (inout), DIMENSION(jpi,jpj,jpk) ::   &
          pah          ! adimensional vertical profile

      !! * Local variables
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs, zcf  ! temporary scalars
      !!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )
      zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         zcf = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
         pah(:,:,jk) = zcf
      END DO

      ! Control print
      IF(lwp .AND. ld_print ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(jk)
         END DO
      ENDIF

   END SUBROUTINE ldf_zpf_1d_3d


   SUBROUTINE ldf_zpf_3d( ld_print, pdam, pwam, pbot, pdep, pah )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_zpf  ***
      !!
      !! ** Purpose :   vertical adimensional profile for eddy coefficient
      !!
      !! ** Method  :   3D for partial step or s-coordinate
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      LOGICAL , INTENT (in   ) :: ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT (in   ) ::   &
          pdam,     &  ! depth of the inflection point
          pwam,     &  ! width of inflection
          pbot         ! reduction factor (surface value / bottom value)
      REAL(wp), INTENT (in   ), DIMENSION(jpi,jpj,jpk) ::   &
          pdep         ! dep of the gridpoint (T, U, V, F)
      REAL(wp), INTENT (inout), DIMENSION(jpi,jpj,jpk) ::   &
          pah          ! adimensional vertical profile

      !! * Local variables
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
      !!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept(1    ) ) / pwam )   
      zm01 = TANH( ( pdam - gdept(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         pah(:,:,jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(:,:,jk) ) / pwam )  )
      END DO

      ! Control print
      IF(lwp .AND. ld_print ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(1,1,jk)
         END DO
      ENDIF

   END SUBROUTINE ldf_zpf_3d
   !!======================================================================
END MODULE ldfdyn
