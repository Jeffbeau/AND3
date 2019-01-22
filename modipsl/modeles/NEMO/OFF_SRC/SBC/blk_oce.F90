MODULE blk_oce
   !!======================================================================
   !!                 ***  MODULE  blk_oce  ***
   !! Bulk   :  bulk parameter and  variables defined in memory 
   !!======================================================================
   !! History :
   !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/SBC/blk_oce.F90,v 1.2 2005/11/16 16:14:39 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
   !!----------------------------------------------------------------------
   !! ' key_flx_bulk_monthly or defined key_flx_bulk_daily             bulk
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce          ! ocean parameters

   IMPLICIT NONE

   LOGICAL, PUBLIC ::   l_bulk = .TRUE.   !: 
   
   !!----------------------------------------------------------------------
   !! bulk common variables
   !!----------------------------------------------------------------------

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &  !:
      vatm     ,      &  !: wind speed
      catm               !: percent of cloud cover

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC ::   l_bulk = .FALSE.  !:
#endif
   
   !!----------------------------------------------------------------------
END MODULE blk_oce
