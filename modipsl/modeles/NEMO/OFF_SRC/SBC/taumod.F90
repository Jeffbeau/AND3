MODULE taumod
   !!======================================================================
   !!                       ***  MODULE  taumod  ***
   !! Ocean forcing : stress at the the ocean surface
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   tau          : define the surface stress for the ocean
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/SBC/taumod.F90,v 1.2 2005/11/16 16:14:39 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce
   USE lib_mpp
   USE in_out_manager ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Share modules variables
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &
      taux, tauy          !: surface stress components 

   !!======================================================================
END MODULE taumod
