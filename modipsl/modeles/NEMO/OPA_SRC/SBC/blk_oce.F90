MODULE blk_oce
   !!======================================================================
   !!                 ***  MODULE  blk_oce  ***
   !! Bulk   :  bulk parameter and  variables defined in memory 
   !!======================================================================
   !! History :
   !!   8.5  !  02-11  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SBC/blk_oce.F90,v 1.4 2005/09/22 10:58:15 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
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
      watm     ,      &  !: precipitation
      tatm     ,      &  !: atmospheric temperature
      hatm     ,      &  !: relative humidity
      vatm     ,      &  !: wind speed
      catm               !: percent of cloud cover

   REAL(wp), PUBLIC, DIMENSION(jpi,jpj)    ::   &  !:
      gsst               !: SST mean on nfbulk ocean time step

   REAL(wp) ::        &
      yearday  ,      &  !: number of days per year
      rdtbs2             !: bulk time step divide by 2
#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC ::   l_bulk = .FALSE.  !:
#endif
   
   INTEGER ::         & !!: namdom : space/time domain (namlist)
      nfbulk =  5        !: bulk computation frequency 
   !!----------------------------------------------------------------------
END MODULE blk_oce
