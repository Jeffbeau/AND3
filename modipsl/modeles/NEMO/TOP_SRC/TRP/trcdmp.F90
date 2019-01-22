!!DB 2009.08.24 -- eliminate non BGCM code options 
!!NB: Code is broken if tracer damping is desired 
!!     ====> must modify this code, or recover older version
MODULE trcdmp
   !!======================================================================
   !!                       ***  MODULE  trcdmp  ***
   !! Ocean physics: internal restoring trend on passive tracers
   !!======================================================================
#if  defined key_passivetrc && defined key_trcdmp 
  USE in_out_manager
!!DB
   LOGICAL , PUBLIC, PARAMETER ::   lk_trcdmp = .FALSE.    !: internal damping flag
CONTAINS
   SUBROUTINE trc_dmp( kt )        ! Empty routine
      INTEGER, INTENT(in) :: kt

      if(lwp) WRITE(numout,*) 'DBG -- trc_dmp: key_trcdmp inactivated ===> code will not work', kt

   END SUBROUTINE trc_dmp

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO internal damping
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_trcdmp = .FALSE.    !: internal damping flag
CONTAINS
   SUBROUTINE trc_dmp( kt )        ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_dmp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_dmp
#endif

   !!======================================================================
END MODULE trcdmp
