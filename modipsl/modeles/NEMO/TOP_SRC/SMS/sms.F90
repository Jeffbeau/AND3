!!DB 2009.08.24 -- eliminate non BGCM code options 
MODULE sms
   !!======================================================================
   !!                        ***  sms  ***
   !! passive tracers :   set the passive tracers variables
   !!======================================================================
   !! History :
   !!   9.0  !  04-03  (C. Ethe)  Free form and module
   !!----------------------------------------------------------------------
   !!  TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/SMS/sms.F90,v 1.6 2005/11/14 16:42:42 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   !! * Modules used
#if defined key_passivetrc

   USE par_oce
   USE par_trc
!!DB
!   USE par_sms
   IMPLICIT NONE

#else
   !!======================================================================
   !!  Empty module : No passive tracer 
   !!======================================================================
#endif

END MODULE sms
