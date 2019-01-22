MODULE zdf_oce
   !!======================================================================
   !!              ***  MODULE  zdf_oce  ***
   !! Ocean physics : define vertical mixing variables
   !!=====================================================================
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf_oce.F90,v 1.4 2005/09/02 15:02:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   zdf_init    : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! mesh and scale factors

   IMPLICIT NONE
   PRIVATE

   !! * Share Module variables
   LOGICAL, PARAMETER, PUBLIC ::    &   !:
#if defined key_zdfcst   ||   defined key_esopa
      lk_zdfcst        = .TRUE.         !: constant vertical mixing flag
#else
      lk_zdfcst        = .FALSE.        !: constant vertical mixing flag
#endif
   LOGICAL, PUBLIC ::    &   !:
      ln_zdfevd        = .TRUE.  ,   &  !: convection: enhanced vertical diffusion flag
      ln_zdfnpc        = .FALSE.        !: convection: non-penetrative convection flag

   LOGICAL, PUBLIC ::    &   !:
      l_trazdf_exp     = .FALSE. ,   &  !: ???
      l_trazdf_imp     = .FALSE. ,   &  !: 
      l_dynzdf_exp     = .FALSE. ,   &  !: 
      l_dynzdf_imp     = .TRUE.  ,   &  !:
      l_dynzdf_imp_tsk = .FALSE.        !:

   INTEGER, PUBLIC ::    & !!: namzdf:  vertical diffusion
      n_zdfexp = 3    ,  &  !: number of sub-time step (explicit time stepping)
      nevdm    = 1          !: =0/1 flag to apply enhanced avm or not
 
   REAL(wp), PUBLIC ::   & !!: namzdf   vertical diffusion
      avm0  = 1.e-4_wp,  &  !: vertical eddy viscosity (m2/s)
      avt0  = 1.e-5_wp,  &  !: vertical eddy diffusivity (m2/s)
      avevd = 1._wp         !: vertical eddy coeff. for enhanced vert. diff. (m2/s)

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
      avmu,              &  !: vertical viscosity coeff. at uw-, vw-points
      avmv,              &  !: vertical viscosity coeff. at uw-, vw-points
      avt ,              &  !: vertical diffusivity coeff. at w-point
      avt_evd,           &  !: convection: enhanced vertical diffusivity coeff. at w-point
      avmu_evd              !: convection: enhanced vertical viscosity   coeff. at w-point
 
   REAL(wp), PUBLIC, DIMENSION(jpk) ::   &   !:
      avmb, avtb            !: background profile of avm and avt
 
   !!======================================================================
END MODULE zdf_oce
