MODULE zdf_oce
   !!======================================================================
   !!              ***  MODULE  zdf_oce  ***
   !! Ocean physics : define vertical mixing variables
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   zdf_init    : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/ZDF/zdf_oce.F90,v 1.2 2005/11/16 16:16:03 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! mesh and scale factors

   IMPLICIT NONE
   PRIVATE

   !! * Share Module variables

   REAL(wp), PUBLIC ::   & !!: namzdf   vertical diffusion
      avt0  = 1.e-5_wp     !: vertical eddy diffusivity (m2/s)

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      avt                   !: vertical diffusivity coeff. at w-point
 
   REAL(wp), PUBLIC, DIMENSION(jpk) ::   &   !:
      avtb            !: background profile of avm and avt

   LOGICAL, PUBLIC ::    &   !:
      ln_zdfnpc        = .FALSE.        !: convection: non-penetrative convection flag
 
   !!======================================================================
END MODULE zdf_oce
