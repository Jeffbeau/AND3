MODULE dynspg_oce
   !!----------------------------------------------------------------------
   !!                       ***  MODULE dynspg_oce  ***
   !!       
   !! ** Purpose :   Define in memory all the ocean space domain variables
   !!----------------------------------------------------------------------
   !! Modules used
   USE par_oce          ! ocean parameters

   IMPLICIT NONE
   PUBLIC           
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynspg_oce.F90,v 1.1 2005/12/28 09:25:06 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

   !! Surface pressure gradient logicals
   !! ----------------------------------
#if   defined key_dynspg_exp   ||  defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_exp = .TRUE.  !: Explicit free surface flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_exp = .FALSE. !: Explicit free surface flag
#endif
#if   defined key_dynspg_ts   ||  defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_ts  = .TRUE.  !: Free surface with time splitting flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_ts  = .FALSE. !: Free surface with time splitting flag
#endif
#if   defined key_dynspg_flt  ||  defined key_esopa
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_flt = .TRUE.  !: Filtered free surface cst volume flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_flt = .FALSE. !: Filtered free surface cst volume flag
#endif
#if   defined key_dynspg_rl
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_rl  = .TRUE.  !: Rigid-lid flag
#else
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_rl  = .FALSE. !: Rigid-lid flag
#endif

#if   defined key_dynspg_ts   ||  defined key_esopa
   !! Time splitting variables
   !! ------------------------
      REAL(wp), PUBLIC, DIMENSION(jpi,jpj) :: & ! variables averaged over the barotropic loop
         sshn_b, sshb_b,               &  ! sea surface heigth (now, before)
         un_b  , vn_b                     ! vertically integrated horizontal velocities (now)
      REAL(wp), PUBLIC, DIMENSION(jpi,jpj) :: & ! variables of the explicit barotropic loop
         sshn_e, ssha_e,               &  ! sea surface heigth (now,after)
         ua_e  , va_e                     ! vertically integrated horizontal velocities (after)
#endif

   !!----------------------------------------------------------------------

END MODULE dynspg_oce
