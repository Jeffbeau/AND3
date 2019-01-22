!!DB 2009.08.24 -- eliminate non BGCM code options 
!!DL 2012.09.27 -- added option for BGCM_02 
MODULE initrc
   !!================================================
   !!
   !!                       *** MODULE initrc ***
   !! Initialisation the tracer model
   !!================================================
                                                                                                                            
#if defined key_passivetrc

#if defined key_BGCM_01
!!lib_bgcm_01 contains ini_trc() and associated routines
  USE lib_bgcm_01

#elif defined key_BGCM_02
!!lib_bgcm_02 contains ini_trc() and associated routines
  USE lib_bgcm_02

#else
   !!======================================================================
   !!  Empty module : No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE ini_trc
      
   END SUBROUTINE ini_trc
#endif

#endif

END MODULE initrc 


