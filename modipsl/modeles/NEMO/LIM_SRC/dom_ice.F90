MODULE dom_ice
   !!======================================================================
   !!                   ***  MODULE  dom_ice  ***
   !! LIM Sea Ice :   Domain  variables
   !!======================================================================
   !! History :
   !!   2.0  !  03-08  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------
   !!   LIM 2.0, UCL-LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/dom_ice.F90,v 1.5 2006/03/21 08:42:22 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_ice

   IMPLICIT NONE
   PRIVATE

   !! * Share module variables
   LOGICAL, PUBLIC ::       &  !:
      l_jeq     = .TRUE. ,  &  !: Equator inside the domain flag
      ln_limini = .FALSE.,  &  !: Ice initialization state
      ln_limdmp = .FALSE.      !: Ice damping

   INTEGER, PUBLIC ::   &  !:
      njeq , njeqm1        !: j-index of the equator if it is inside the domain
      !                    !  (otherwise = jpj+10 (SH) or -10 (SH) )

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      fs2cor ,          &  !: coriolis factor
      fcor   ,          &  !: coriolis coefficient
      covrai ,          &  !: sine of geographic latitude
      area   ,          &  !: surface of grid cell 
      tms    , tmu         !: temperature and velocity points masks

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2,2) ::   &  !:
      wght   ,          &  !: weight of the 4 neighbours to compute averages
      akappa ,          &  !: first group of metric coefficients
      bkappa               !: third group of metric coefficients

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,2,2,2,2) ::   &  !:
      alambd               !: second group of metric coefficients

   !!======================================================================
END MODULE dom_ice
