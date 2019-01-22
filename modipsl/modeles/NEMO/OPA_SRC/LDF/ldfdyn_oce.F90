MODULE ldfdyn_oce
   !!======================================================================
   !!                  ***  MODULE  ldfdyn_oce  ***
   !! Ocean physics:  lateral momentum mixing coefficient defined in memory 
   !!======================================================================
   !!
   !! ** Purpose :
   !!       - Define in memory lateral momentum mixing coefficients
   !!
   !! History :
   !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldfdyn_oce.F90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce      ! ocean parameters

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !! Lateral eddy viscosity coefficients (dynamics)
   !!----------------------------------------------------------------------

   LOGICAL  ::                      & !!! ** lateral mixing namelist (nam_dynldf) **
      ln_dynldf_lap   = .TRUE.  ,   &  ! laplacian operator
      ln_dynldf_bilap = .FALSE. ,   &  ! bilaplacian operator
      ln_dynldf_level = .FALSE. ,   &  ! iso-level direction
      ln_dynldf_hor   = .TRUE.  ,   &  ! horizontal (geopotential) direction
      ln_dynldf_iso   = .FALSE.        ! iso-neutral direction

   REAL(wp) ::                      & !!! ** lateral mixing namelist (nam_dynldf) **
      ahm0  = 40000._wp ,   &  ! lateral eddy viscosity (m2/s)
      ahmb0 =     0._wp        ! lateral background eddy viscosity (m2/s)

   LOGICAL  ::                      &  ! flag of the lateral diff. scheme used
      l_dynldf_lap         ,        &  ! iso-level laplacian operator
      l_dynldf_bilap       ,        &  ! iso-level bilaplacian operator
      l_dynldf_bilapg      ,        &  ! geopotential bilap. (s-coord)
      l_dynldf_iso         ,        &  ! iso-neutral laplacian or horizontal lapacian (s-coord)
      l_dynzdf_iso                     ! iso-neutral laplacian or horizontal lapacian (s-coord)
   

#if defined key_dynldf_c3d
   REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &  ! ** 3D coefficients **
#elif defined key_dynldf_smag
   REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &  ! ** 3D coefficients **
#elif defined key_dynldf_c2d
   REAL(wp), DIMENSION(jpi,jpj)     ::   &  ! ** 2D coefficients **
#elif defined key_dynldf_c1d
   REAL(wp), DIMENSION(jpk)         ::   &  ! ** 2D coefficients **
#else
   REAL(wp)                         ::   &  ! ** 0D coefficients **
#endif
      ahm1, ahm2, ahm3, ahm4,tmph                ! ????

   !!----------------------------------------------------------------------
END MODULE ldfdyn_oce
