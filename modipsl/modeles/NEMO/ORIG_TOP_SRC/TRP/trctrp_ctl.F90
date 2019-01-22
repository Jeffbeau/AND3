MODULE trctrp_ctl
   !!==============================================================================
   !!                       ***  MODULE  trctrp_ctl  ***
   !! Ocean passive tracers:  transport option control
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_trp_ctl  : control the different options of transport
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and active tracers variables
   USE trc                 ! ocean passive tracers variables
   USE trctrp_lec          ! passive tracers transport

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_trp_ctl   

   !! * Module variable
#if defined key_trcldf_eiv
      LOGICAL, PARAMETER ::   lk_trcldf_eiv   = .TRUE.   !: eddy induced velocity flag
#else   
      LOGICAL, PARAMETER ::   lk_trcldf_eiv   = .FALSE.  !: eddy induced velocity flag
#endif

   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trctrp_ctl.F90,v 1.10 2006/04/11 13:49:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_trp_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trc_trp_ctl  ***
      !!                
      !! ** Purpose :   Control the consistency between cpp options for 
      !!      tracer transport
      !!
      !! History :
      !!   9.0  !  04-0.  (C. Ethe) 
      !!----------------------------------------------------------------------

      !!----------------------------------------------------------------------
      !!  TOP 1.0 , LOCEAN-IPSL (2005) 
      !!----------------------------------------------------------------------

      !! Control of Advection scheme options
      CALL trc_adv_ctl

      !! Control of Lateral diffusion scheme options
      CALL trc_ldf_ctl

      !! Control of Vertival diffusion scheme options
      CALL trc_zdf_ctl

      !! Control of Newtonian damping  options
      IF(lwp) THEN
         WRITE(numout,*) ' *** Tracer damping option'
         WRITE(numout,*)
      ENDIF

#if defined key_trcdmp
      IF(lwp) THEN 
         WRITE(numout,*)' key_trcdmp is defined'
         WRITE(numout,*)' Check trcdmp ROUTINE '
         WRITE(numout,*)'  '
      ENDIF 
      CALL trc_dmp_ctl
#else
      IF (lwp) WRITE(numout,*) ' No tracer damping'
#endif


   END SUBROUTINE trc_trp_ctl

   SUBROUTINE trc_adv_ctl
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trc_adv_ctl  ***
      !!                
      !! ** Purpose :   Control the consistency between cpp options for 
      !!      tracer advection schemes
      !!
      !! History :
      !!   8.5  !  02-11  (G. Madec)  Original code
      !!   9.0  !  04-0.  (C. Ethe)  adapted for passive tracers
      !!----------------------------------------------------------------------

      !! * Local declarations
      INTEGER ::   ioptio


      !!----------------------------------------------------------------------
      !!  TOP 1.0 , LOCEAN-IPSL (2005) 
      !!----------------------------------------------------------------------

      ! Control of Advection scheme options
      ! -----------------------------------
      ioptio = 0
      IF( ln_trcadv_cen2   )   ioptio = ioptio + 1
      IF( ln_trcadv_tvd    )   ioptio = ioptio + 1
      IF( ln_trcadv_muscl  )   ioptio = ioptio + 1
      IF( ln_trcadv_muscl2 )   ioptio = ioptio + 1
      IF( ln_trcadv_smolar )   ioptio = ioptio + 1

      IF( lk_esopa ) THEN
         IF(lwp) WRITE(numout,*) ' esopa control : the use of all scheme is forced'
         ln_trcadv_cen2   = .TRUE.
         ln_trcadv_tvd    = .TRUE.
         ln_trcadv_muscl  = .TRUE.
         ln_trcadv_muscl2 = .TRUE.
         ln_trcadv_smolar = .TRUE.
      ELSEIF( ioptio > 1 .OR. ioptio == 0 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) ' Choose one advection scheme in namelist nam_trcadv'
         IF(lwp) WRITE(numout,*) '        ***                              ***********'
         nstop = nstop + 1
      ENDIF

      IF( n_cla == 1 .AND. .NOT. ln_trcadv_cen2 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '     cross-land advection only with 2nd order advection scheme'
         nstop = nstop + 1
      ENDIF

      IF( lk_trccfg_1d ) THEN
         ln_trcadv_cen2   = .FALSE.    ;  ln_trcadv_tvd    = .FALSE. ; ln_trcadv_muscl  = .FALSE.
         ln_trcadv_muscl2 = .FALSE.    ;  ln_trcadv_smolar = .FALSE.
         IF(lwp) WRITE(numout,*) ' *******  1D configuration : No advection on passive tracers *******'
         IF(lwp) WRITE(numout,*) ' *******                                                     *******'
      ENDIF

   END SUBROUTINE trc_adv_ctl

   SUBROUTINE trc_ldf_ctl
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_ldf_ctl  ***
      !! 
      !! ** Purpose :   Control the consistency between cpp options for 
      !!      tracer lateral diffusion 
      !!
      !! History :
      !!   9.0  !  03-04  (C. Ethe) 
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ioptio               ! ???
      LOGICAL ::   ll_print = .FALSE.   ! =T print eddy coef. in numout      

      !!----------------------------------------------------------------------
      !!  TOP 1.0 , LOCEAN-IPSL (2005) 
      !!----------------------------------------------------------------------

      ! Parameter control

      ! control the input
      ioptio = 0
      IF( ln_trcldf_lap   )   ioptio = ioptio + 1
      IF( ln_trcldf_bilap )   ioptio = ioptio + 1
      IF( ioptio /= 1 )   THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) '          use ONE of the 2 lap/bilap operator type on tracer'
          nstop = nstop + 1
      ENDIF
      ioptio = 0
      IF( ln_trcldf_level )   ioptio = ioptio + 1
      IF( ln_trcldf_hor   )   ioptio = ioptio + 1
      IF( ln_trcldf_iso   )   ioptio = ioptio + 1
      IF( ioptio /= 1 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          use only ONE direction (level/hor/iso)'
         nstop = nstop + 1
      ENDIF

      ! ... Choice of the lateral scheme used
      IF( lk_trcldf_eiv ) THEN
         IF(lwp) WRITE(numout,*) '          eddy induced velocity on tracers'
            IF( .NOT.ln_trcldf_iso .OR. ln_trcldf_bilap ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion'
            nstop = nstop + 1
         ENDIF
      ENDIF

      IF( lk_sco ) THEN          ! s-coordinates: rotation required for horizontal or isopycnal mixing
         IF( ( ln_trcldf_iso .OR. ln_trcldf_hor ) .AND. .NOT.lk_ldfslp ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
            IF( .NOT.lk_esopa )   nstop = nstop + 1
         ENDIF
      ELSE                       ! z-coordinates with/without partial step:
         ln_trcldf_level = ln_trcldf_level .OR. ln_trcldf_hor      ! level diffusion = horizontal diffusion
         ln_trcldf_hor   = .FALSE.
         IF(lwp) WRITE(numout,*) '          horizontal mixing in z-coord or partial steps: force ln_trcldf_level = T'
         IF(lwp) WRITE(numout,*) '                                                  and    force ln_trcldf_hor   = F'
         IF( ln_trcldf_iso .AND. .NOT.lk_ldfslp ) THEN             ! rotation required for isopycnal mixing
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          the rotation of the diffusive tensor require key_ldfslp'
            IF( .NOT.lk_esopa )   nstop = nstop + 1
         ENDIF
      ENDIF

      l_trcldf_lap     =  ln_trcldf_lap   .AND. ln_trcldf_level     ! iso-level   laplacian operator
      l_trcldf_bilap   =  ln_trcldf_bilap .AND. ln_trcldf_level     ! iso-level bilaplacian operator
      l_trcldf_bilapg  =  ln_trcldf_bilap .AND. ln_trcldf_hor       ! geopotential bilap. (s-coord)
      l_trcldf_iso     =  ln_trcldf_lap   .AND.                  &  ! laplacian operator
         &                   ( ln_trcldf_iso   .OR.  ln_trcldf_hor )  &  ! iso-neutral (z-coord) or horizontal (s-coord)
         &                                     .AND. .NOT.lk_zps
      l_trcldf_iso_zps =       ln_trcldf_lap   .AND.                  &  ! laplacian operator
         &                   ( ln_trcldf_iso   .OR.  ln_trcldf_hor )  &  ! iso-neutral (partial steps)
         &                                     .AND. lk_zps              ! or geopotential in mixed partial steps/s-coord
      l_trczdf_iso    = .FALSE.
      l_trczdf_iso_vo = .FALSE.
      IF( l_trcldf_iso     )   l_trczdf_iso = .TRUE.
      IF( l_trcldf_iso_zps )   l_trczdf_iso = .TRUE.
#if defined key_vectopt_memory
      IF( l_trczdf_iso ) THEN
         l_trczdf_iso    = .FALSE.
         l_trczdf_iso_vo = .TRUE.
      ENDIF
#endif

 
      ioptio = 0
      IF( l_trcldf_lap     )   ioptio = ioptio + 1
      IF( l_trcldf_bilap   )   ioptio = ioptio + 1
      IF( l_trcldf_bilapg  )   ioptio = ioptio + 1
      IF( l_trcldf_iso     )   ioptio = ioptio + 1
      IF( l_trcldf_iso_zps )   ioptio = ioptio + 1
      IF( ioptio /= 1 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          this combination of operator and direction has not been implemented'
         nstop = nstop + 1
      ENDIF

      IF( lk_esopa ) THEN
         l_trcldf_lap = .TRUE.   ;   l_trcldf_bilap   = .TRUE.   ;   l_trcldf_bilapg  = .TRUE.
         l_trcldf_iso = .TRUE.   ;   l_trcldf_iso_zps = .TRUE.
         l_trczdf_iso = .TRUE.   ;   l_trczdf_iso_vo  = .TRUE.
         IF(lwp ) WRITE(numout,*) '          esopa test: use all lateral physics options'
      ENDIF

      IF( .NOT. ln_trcldf_diff .OR. lk_trccfg_1d ) THEN
         l_trcldf_lap = .FALSE.   ;   l_trcldf_bilap   = .FALSE.   ;   l_trcldf_bilapg  = .FALSE.
         l_trcldf_iso = .FALSE.   ;   l_trcldf_iso_zps = .FALSE.
         l_trczdf_iso = .FALSE.   ;   l_trczdf_iso_vo  = .FALSE.
         IF(lwp ) WRITE(numout,*) '************* No lateral physics on passive tracers *****************'
         IF(lwp ) WRITE(numout,*) '*************                                       *****************'
      ELSE
         ! ... Space variation of eddy coefficients
         ioptio = 0
#if defined key_traldf_c3d
         IF(lwp) WRITE(numout,*) 'tracer mixing coef. = F( latitude, longitude, depth)'
         ioptio = ioptio + 1
#endif
#if defined key_traldf_c2d
         IF(lwp) WRITE(numout,*) 'tracer mixing coef. = F( latitude, longitude)'
         ioptio = ioptio + 1
#endif
#if defined key_traldf_c1d
         IF(lwp) WRITE(numout,*) 'tracer mixing coef. = F( depth )'
         ioptio = ioptio + 1
         IF( lk_sco ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '          key_traldf_c1d cannot be used in s-coordinate (key_s_coord)'
            nstop = nstop + 1
         ENDIF
#endif
         IF( ioptio == 0 ) THEN
            IF(lwp) WRITE(numout,*) ' tracer mixing coef. = constant (default option)'
         ELSEIF( ioptio > 1 ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' use only one of the following keys:',   &
               &                    ' key_traldf_c3d, key_traldf_c2d, key_traldf_c1d'
            nstop = nstop + 1
         ENDIF
         
         IF( l_trcldf_bilap .OR. l_trcldf_bilapg ) THEN
            IF(lwp) WRITE(numout,*) '  biharmonic tracer diffusion'
            IF( ahtrc0 > 0 .AND. .NOT. lk_esopa ) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*) ' The horizontal diffusivity coef. aht0 must be negative'
               nstop = nstop + 1
            ENDIF
         ELSE
            IF(lwp) WRITE(numout,*) ' harmonic tracer diffusion (default)'
            IF( ahtrc0 < 0 .AND. .NOT. lk_esopa ) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*) 'The horizontal diffusivity coef. aht0 must be positive'
               nstop = nstop + 1
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE trc_ldf_ctl

   SUBROUTINE trc_zdf_ctl
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf_ctl  ***
      !! 
      !! ** Purpose :     Control the consistency between cpp options for 
      !!      tracer vertical diffusion
      !!
      !!   9.0  !  04-03  (C. Ethe)  
      !!----------------------------------------------------------------------
      !! * Local declarations

      !!----------------------------------------------------------------------
      !!  TOP 1.0 , LOCEAN-IPSL (2005) 
      !!----------------------------------------------------------------------

      ! Parameter & key controls
      ! ------------------------
      ! ... vertical mixing
      ! time stepping scheme (N.B. TKE scheme => force the use of implicit scheme)
#if defined key_zdftke
      l_trczdf_exp = .FALSE.          ! use implicit scheme
      l_trczdf_imp = .TRUE. 
#else
      IF( ln_zdfexp  ) THEN  
         l_trczdf_exp = .TRUE.           ! use explicit scheme
         l_trczdf_imp = .FALSE.
      ELSE
         l_trczdf_exp = .FALSE.          ! use implicit scheme
         l_trczdf_imp = .TRUE. 
      ENDIF
#endif

      IF( l_trczdf_iso .OR. l_trczdf_iso_vo ) THEN  
         l_trczdf_exp = .FALSE.          ! iso-neutral diffusion : 
         l_trczdf_imp = .FALSE.          ! implicit scheme included in iso-neutral routine
      ENDIF

#if defined key_esopa
      l_trczdf_exp = .TRUE.           ! esopa: use all options
      l_trczdf_imp = .TRUE.
#endif


   END SUBROUTINE trc_zdf_ctl

   SUBROUTINE trc_dmp_ctl
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_dmp_ctl  ***
      !! 
      !! ** Purpose :    Control the consistency between cpp options for 
      !!      tracer newtonian damping 
      !!
      !!
      !! History :
      !!   9.0  !  04-03  (C. Ethe) 
      !!----------------------------------------------------------------------
#if defined key_trcdmp

      SELECT CASE ( ndmptr )

      CASE ( -1 )               ! ORCA: damping in Red & Med Seas only
         IF(lwp) WRITE(numout,*) '          tracer damping in the Med & Red seas only'

      CASE ( 1:90 )             ! Damping poleward of 'ndmptr' degrees
         IF(lwp) WRITE(numout,*) '          tracer damping poleward of', ndmptr, ' degrees'

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for ndmptr = ', ndmptr
         nstop = nstop + 1

      END SELECT


      SELECT CASE ( nmldmptr )

      CASE ( 0 )                ! newtonian damping throughout the water column
         IF(lwp) WRITE(numout,*) '          tracer damping throughout the water column'

      CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
         IF(lwp) WRITE(numout,*) '          no tracer damping in the turbocline'

      CASE ( 2 )                ! no damping in the mixed layer 
         IF(lwp) WRITE(numout,*) '          no tracer damping in the mixed layer'

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for nmldmptr = ', nmldmptr
         nstop = nstop + 1

      END SELECT
#endif
 
   END SUBROUTINE trc_dmp_ctl

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_trp_ctl             ! Empty routine
   END SUBROUTINE trc_trp_ctl
#endif
   
  !!======================================================================
END MODULE trctrp_ctl
