MODULE dynspg_oce
   !!======================================================================
   !!                   ***  MODULE  dynspg_oce  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/DYN/dynspg_oce.F90,v 1.1 2006/01/13 14:17:10 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------
   
   IMPLICIT NONE
   PRIVATE

#if defined key_dynspg_rl
   !!----------------------------------------------------------------------
   !!   'key_dynspg_rl'      rigid lid
   !!----------------------------------------------------------------------
   !! * Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_rl = .TRUE. !: rigid flag flag

#elif defined key_dynspg_flt || defined key_dynspg_ts || defined key_dynspg_exp
   !!----------------------------------------------------------------------
   !!   Default case :       free surface 
   !!   update the momentum trend with the surface pressure
   !!   gradient in the free surface case with vector
   !!   optimization
   !!----------------------------------------------------------------------
   !! * Shared module variables
   LOGICAL, PUBLIC, PARAMETER ::   lk_dynspg_rl = .FALSE.     !: rigid flag 
#endif
   
   !!======================================================================
END MODULE dynspg_oce
