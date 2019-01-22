MODULE trc_oce
   !!======================================================================
   !!                      ***  MODULE  trc_oce  ***
   !! Ocean passive tracer  :  share SMS/Ocean variables
   !!======================================================================
   !! History :
   !!   9.0  !  04-03  (C. Ethe)  F90: Free form and module
   !!----------------------------------------------------------------------
#if defined key_passivetrc && defined key_trc_pisces
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/trc_oce.F90,v 1.3 2005/03/27 18:34:49 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   USE sms , ONLY :  &
      etot3    =>   etot3   !!:  Biological fluxes for light
   !! Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_qsr_sms = .TRUE. 
#else
   !!----------------------------------------------------------------------
   !! Default option                         No Biological fluxes for light          
   !!----------------------------------------------------------------------
   USE par_oce
   LOGICAL, PUBLIC, PARAMETER ::   lk_qsr_sms = .FALSE. 
   REAL(wp), PUBLIC , DIMENSION (jpi,jpj,jpk) :: &
      etot3
#endif

END MODULE trc_oce
