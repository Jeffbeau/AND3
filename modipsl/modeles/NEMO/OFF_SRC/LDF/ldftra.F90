MODULE ldftra
   !!======================================================================
   !!                       ***  MODULE  ldftra  ***
   !! Ocean physics:  lateral diffusivity coefficient 
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   ldf_tra_init : initialization, namelist read, and parameters control
   !!   ldf_tra_c3d   : 3D eddy viscosity coefficient initialization
   !!   ldf_tra_c2d   : 2D eddy viscosity coefficient initialization
   !!   ldf_tra_c1d   : 1D eddy viscosity coefficient initialization
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE ldfslp          ! ???
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC ldf_tra_init   ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/LDF/ldftra.F90,v 1.2 2005/11/16 16:13:32 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE ldf_tra_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra_init  ***
      !! 
      !! ** Purpose :   initializations of the horizontal ocean tracer physics
      !!
      !! ** Method :
      !!      Direction of lateral diffusion (tracers and/or momentum)
      !!        ln_traldf_iso  = T : initialize the slope arrays to zero
      !!        ln_traldf_geop = T : initialise the slope arrays to the i- and
      !!                            j-slopes of s-surfaces
      !!      Eddy diffusivity and eddy induced velocity cefficients:
      !!         default option   : constant coef. aht0, aeiv0 (namelist)
      !!        'key_traldf_c1d': depth dependent coef. defined in 
      !!                            in ldf_tra_c1d routine
      !!        'key_traldf_c2d': latitude and longitude dependent coef.
      !!                            defined in ldf_tra_c2d routine
      !!        'key_traldf_c3d': latitude, longitude, depth dependent coef.
      !!                            defined in ldf_tra_c3d routine
      !!
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
      INTEGER ::   ioptio               ! ???
      LOGICAL ::   ll_print = .FALSE.   ! =T print eddy coef. in numout
       
      NAMELIST/nam_traldf/ ln_traldf_lap  , ln_traldf_bilap,                &
         &                 ln_traldf_level, ln_traldf_hor, ln_traldf_iso,   &
         &                 aht0, ahtb0, aeiv0
      !!----------------------------------------------------------------------

      !  Define the lateral tracer physics parameters
      ! =============================================
    
      ! Read Namelist nam_traldf : Lateral physics on tracers
      REWIND( numnam )
      READ  ( numnam, nam_traldf )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_tra : lateral tracer physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '          Namelist nam_traldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '             laplacian operator          ln_traldf_lap   = ', ln_traldf_lap
         WRITE(numout,*) '             bilaplacian operator        ln_traldf_bilap = ', ln_traldf_bilap
         WRITE(numout,*) '             iso-level                   ln_traldf_level = ', ln_traldf_level
         WRITE(numout,*) '             horizontal (geopotential)   ln_traldf_hor   = ', ln_traldf_hor
         WRITE(numout,*) '             iso-neutral                 ln_traldf_iso   = ', ln_traldf_iso
         WRITE(numout,*) '             lateral eddy diffusivity             aht0   = ', aht0
         WRITE(numout,*) '             background hor. diffusivity          ahtb0  = ', ahtb0
         WRITE(numout,*) '             eddy induced velocity coef.          aeiv0  = ', aeiv0
         WRITE(numout,*)
      ENDIF

      ! Parameter control

      ! control the input
      ioptio = 0
      IF( ln_traldf_lap   )   ioptio = ioptio + 1
      IF( ln_traldf_bilap )   ioptio = ioptio + 1
      IF( ioptio /= 1 )   THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) '          use ONE of the 2 lap/bilap operator type on tracer'
          nstop = nstop + 1
      ENDIF
      ioptio = 0
      IF( ln_traldf_level )   ioptio = ioptio + 1
      IF( ln_traldf_hor   )   ioptio = ioptio + 1
      IF( ln_traldf_iso   )   ioptio = ioptio + 1
      IF( ioptio /= 1 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          use only ONE direction (level/hor/iso)'
         nstop = nstop + 1
      ENDIF

      ! ... Choice of the lateral scheme used
      IF( lk_traldf_eiv ) THEN
         IF(lwp) WRITE(numout,*) '          eddy induced velocity on tracers'
            IF( .NOT.ln_traldf_iso .OR. ln_traldf_bilap ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion'
            nstop = nstop + 1
         ENDIF
      ENDIF

      IF( lk_sco ) THEN          ! s-coordinates: rotation required for horizontal or isopycnal mixing
         IF( ( ln_traldf_iso .OR. ln_traldf_hor ) .AND. .NOT.lk_ldfslp ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
            IF( .NOT.lk_esopa )   nstop = nstop + 1
         ENDIF
      ELSE                       ! z-coordinates with/without partial step:
         ln_traldf_level = ln_traldf_level .OR. ln_traldf_hor      ! level diffusion = horizontal diffusion
         ln_traldf_hor   = .FALSE.
         IF(lwp) WRITE(numout,*) '          horizontal mixing in z-coord or partial steps: force ln_traldf_level = T'
         IF(lwp) WRITE(numout,*) '                                                  and    force ln_traldf_hor   = F'
         IF( ln_traldf_iso .AND. .NOT.lk_ldfslp ) THEN             ! rotation required for isopycnal mixing
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
            IF( .NOT.lk_esopa )   nstop = nstop + 1
         ENDIF
      ENDIF

      l_traldf_lap     =       ln_traldf_lap   .AND. ln_traldf_level     ! iso-level   laplacian operator
      l_traldf_bilap   =       ln_traldf_bilap .AND. ln_traldf_level     ! iso-level bilaplacian operator
      l_traldf_bilapg  =       ln_traldf_bilap .AND. ln_traldf_hor       ! geopotential bilap. (s-coord)
      l_traldf_iso     =       ln_traldf_lap   .AND.                  &  ! laplacian operator
         &                   ( ln_traldf_iso   .OR.  ln_traldf_hor )  &  ! iso-neutral (z-coord) or horizontal (s-coord)
         &                                     .AND. .NOT.lk_zps
      l_traldf_iso_zps =       ln_traldf_lap   .AND.                  &  ! laplacian operator
         &                   ( ln_traldf_iso   .OR.  ln_traldf_hor )  &  ! iso-neutral (partial steps)
         &                                     .AND. lk_zps              ! or geopotential in mixed partial steps/s-coord
      l_trazdf_iso    = .FALSE.
      l_trazdf_iso_vo = .FALSE.
      IF( l_traldf_iso     )   l_trazdf_iso = .TRUE.
      IF( l_traldf_iso_zps )   l_trazdf_iso = .TRUE.
#if defined key_vectopt_memory
      IF( l_trazdf_iso ) THEN
         l_trazdf_iso    = .FALSE.
         l_trazdf_iso_vo = .TRUE.
      ENDIF
#endif

      ioptio = 0
      IF( l_traldf_lap     )   ioptio = ioptio + 1
      IF( l_traldf_bilap   )   ioptio = ioptio + 1
      IF( l_traldf_bilapg  )   ioptio = ioptio + 1
      IF( l_traldf_iso     )   ioptio = ioptio + 1
      IF( l_traldf_iso_zps )   ioptio = ioptio + 1
      IF( ioptio /= 1 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          this combination of operator and direction has not been implemented'
         nstop = nstop + 1
      ENDIF
      IF( lk_esopa ) THEN
         l_traldf_lap = .TRUE.   ;   l_traldf_bilap   = .TRUE.   ;   l_traldf_bilapg  = .TRUE.
         l_traldf_iso = .TRUE.   ;   l_traldf_iso_zps = .TRUE.
         l_trazdf_iso = .TRUE.   ;   l_trazdf_iso_vo  = .TRUE.
         IF(lwp ) WRITE(numout,*) '          esopa test: use all lateral physics options'
      ENDIF


      ! ... Space variation of eddy coefficients
      ioptio = 0
#if defined key_traldf_c3d
      IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude, depth)'
      ioptio = ioptio + 1
#endif
#if defined key_traldf_c2d
      IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude)'
      ioptio = ioptio + 1
#endif
#if defined key_traldf_c1d
      IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( depth )'
      ioptio = ioptio + 1
      IF( lk_sco ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          key_traldf_c1d cannot be used in s-coordinate (key_s_coord)'
         nstop = nstop + 1
      ENDIF
#endif
      IF( ioptio == 0 ) THEN
          IF(lwp) WRITE(numout,*) '          tracer mixing coef. = constant (default option)'
        ELSEIF( ioptio > 1 ) THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) '          use only one of the following keys:',   &
             &                    ' key_traldf_c3d, key_traldf_c2d, key_traldf_c1d'
          nstop = nstop + 1
      ENDIF

      IF( l_traldf_bilap .OR. l_traldf_bilapg ) THEN
         IF(lwp) WRITE(numout,*) '          biharmonic tracer diffusion'
         IF( aht0 > 0 .AND. .NOT. lk_esopa ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          The horizontal diffusivity coef. aht0 must be negative'
            nstop = nstop + 1
         ENDIF
      ELSE
         IF(lwp) WRITE(numout,*) '          harmonic tracer diffusion (default)'
         IF( aht0 < 0 .AND. .NOT. lk_esopa ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          The horizontal diffusivity coef. aht0 must be positive'
            nstop = nstop + 1
         ENDIF
      ENDIF


      !  Lateral eddy diffusivity and eddy induced velocity coefficients
      ! ================================================================

#if defined key_traldf_c3d
      CALL ldf_tra_c3d( ll_print )           ! aht = 3D coef. = F( longitude, latitude, depth )
#elif defined key_traldf_c2d
      CALL ldf_tra_c2d( ll_print )           ! aht = 2D coef. = F( longitude, latitude )
#elif defined key_traldf_c1d
      CALL ldf_tra_c1d( ll_print )           ! aht = 1D coef. = F( depth )
#else
                                     ! Constant coefficients
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' inildf: constant eddy diffusivity coef.'
      IF(lwp)WRITE(numout,*) ' ~~~~~~'
      IF(lwp)WRITE(numout,*) '        ahtu = ahtv = ahtw = aht0 = ', aht0
      IF( lk_traldf_eiv ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) ' inildf: constant eddy induced velocity coef.'
         IF(lwp)WRITE(numout,*) ' ~~~~~~  '
         IF(lwp)WRITE(numout,*) '         aeiu = aeiv = aeiw = aeiv0 = ', aeiv0
      ENDIF
#endif

   END SUBROUTINE ldf_tra_init

#if defined key_traldf_c3d
#   include "ldftra_c3d.h90"
#elif defined key_traldf_c2d
#   include "ldftra_c2d.h90"
#elif defined key_traldf_c1d
#   include "ldftra_c1d.h90"
#endif

   !!======================================================================
END MODULE ldftra
