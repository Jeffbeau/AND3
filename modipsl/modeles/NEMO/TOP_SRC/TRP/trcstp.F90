!!DB 2010.08
!!Add more options
!!DB 2009.08 
!!Eliminate non-BGCM code options
MODULE trcstp
   !!======================================================================
   !!                       ***  MODULE trcstp  ***
   !! Time-stepping    : time loop of opa for passive tracer
   !!======================================================================


#if defined key_passivetrc

!!DB
#if defined key_BGCM_01

  USE bgcm_01_model

#elif defined key_PASSIVE_TRACER_01

  USE ptracer_01_model

#elif defined key_BGCM_02

  USE bgcm_02_model

#else
   !!----------------------------------------------------------------------
   !!   Default key                                     NO passive tracers
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_stp( kt )        ! Empty routine
      INTEGER, INTENT(in) :: kt
!      WRITE(*,*) 'trc_stp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_stp
#endif
#endif
   !!======================================================================
END MODULE trcstp
