MODULE ldftra_oce
   !!=====================================================================
   !!                      ***  MODULE  ldftra_oce  ***
   !! Ocean physics :  lateral tracer mixing coefficient defined in memory 
   !!=====================================================================
   !!
   !! ** Purpose : - Define in memory lateral tracer mixing coefficients
   !!
   !! History :
   !!   9.0  !  02-11  (G. Madec)  Original code (from common.h)
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_oce.F90,v 1.3 2005/03/27 18:35:07 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! ocean parameters

   IMPLICIT NONE
   PRIVATE

   !!----------------------------------------------------------------------
   !! Lateral eddy diffusivity coefficients (tracers)
   !!----------------------------------------------------------------------

   LOGICAL , PUBLIC ::              & !!: ** lateral mixing namelist (nam_traldf) **
      ln_traldf_lap   = .TRUE.  ,   &  !: laplacian operator
      ln_traldf_bilap = .FALSE. ,   &  !: bilaplacian operator
      ln_traldf_level = .FALSE. ,   &  !: iso-level direction
      ln_traldf_hor   = .FALSE. ,   &  !: horizontal (geopotential) direction
      ln_traldf_iso   = .TRUE.         !: iso-neutral direction

   REAL(wp), PUBLIC ::              & !!: ** lateral mixing namelist (namldf) **
      aht0  = 2000._wp     ,        &  !: lateral eddy diffusivity (m2/s)
      ahtb0 =    0._wp     ,        &  !: lateral background eddy diffusivity (m2/s)
      aeiv0 = 2000._wp                 !: eddy induced velocity coefficient (m2/s)

   LOGICAL , PUBLIC ::              &  !: flag of the lateral diff. scheme used 
      l_traldf_lap         ,        &  !: iso-level laplacian operator
      l_traldf_bilap       ,        &  !: iso-level bilaplacian operator
      l_traldf_bilapg      ,        &  !: geopotential bilap. (s-coord)
      l_traldf_iso         ,        &  !: iso-neutral laplacian or horizontal lapacian (s-coord)
      l_trazdf_iso         ,        &  !: idem for the vertical component
      l_trazdf_iso_vo      ,        &  !: idem with vectopt_memory
      l_traldf_iso_zps                 !: iso-neutral laplacian (partial steps)

#if defined key_traldf_c3d
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !: ** 3D coefficients **
#elif defined key_traldf_smag
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !: ** 3D coefficients **
#elif defined key_traldf_c2d
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj)     ::   &  !: ** 2D coefficients **
#elif defined key_traldf_c1d
   REAL(wp), PUBLIC, DIMENSION(jpk)         ::   &  !: ** 1D coefficients **
#else
   REAL(wp), PUBLIC                         ::   &  !: ** 0D coefficients **
#endif
      ahtt, ahtu, ahtv, ahtw                !: T-, U-, V-, W-points coefficients


#if defined key_traldf_eiv
   !!----------------------------------------------------------------------
   !!   'key_traldf_eiv'                              eddy induced velocity
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_traldf_eiv   = .TRUE.   !: eddy induced velocity flag
      
# if defined key_traldf_c3d
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !: ** 3D coefficients **
# elif defined key_traldf_c2d
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj)     ::   &  !: ** 2D coefficients **
# elif defined key_traldf_c1d
   REAL(wp), PUBLIC, DIMENSION(jpk)         ::   &  !: ** 1D coefficients **
# else
   REAL(wp), PUBLIC                         ::   &  !: ** 0D coefficients **
# endif
      aeiu, aeiv, aeiw                              !: U-, V-, W-points  induced velocity coef. (m2/s)

# if defined key_diaeiv
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::    &  !:
      u_eiv, v_eiv, w_eiv     !: The three component of the eddy induced velocity (m/s)
# endif

#else
   !!----------------------------------------------------------------------
   !!   Default option :                           NO eddy induced velocity
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_traldf_eiv   = .FALSE.   !: eddy induced velocity flag
   REAL(wp), PUBLIC ::   aeiu, aeiv, aeiw
#endif

   !!----------------------------------------------------------------------
END MODULE ldftra_oce
