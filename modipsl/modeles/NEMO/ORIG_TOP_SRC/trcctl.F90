MODULE trcctl
   !!==========================================================================
   !!
   !!                       *** MODULE trcctl ***
   !!
   !! Only for passive tracer
   !! control the cpp options for the run and IF files are availables
   !! control also consistancy between options and namelist values
   !!  O.Aumont and A.El Moussaoui 03/05 F90 
   !!=========================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trcctl.F90,v 1.4 2005/11/14 15:44:40 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   USE sms
   USE trctrp_ctl

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_ctl

CONTAINS

   SUBROUTINE trc_ctl
      !!===========================================================================================
      !!
      !!
      !!                       ROUTINE trcctl
      !!                     ******************
      !!
      !!      we use IF/ENDIF inside #IF defined option-cpp
      !!      FILE name must not exceed 21 characters
      !!      
      !!===========================================================================================
      
      !!----------------------------------------------------------------------
      !! local declarations
      !! ==================
      INTEGER  :: istop, jn
      
      !!---------------------------------------------------------------------
      !!  OPA.9    03/2005  
      !!---------------------------------------------------------------------

      ! 0. Parameter
      ! ------------
      istop = 0

      ! 1. LOGICAL UNIT initialization for specifi! files for passive tracer
      ! --------------------------------------------------------------------
      !     nutwrs : OUTPUT for passive tracer restart UNIT (always used)
      !     nutrst : restart FILE  INPUT  UNIT (lrsttr=.TRUE.)
      !     nutini(jptra) : UNIT for initial FILE for tracer

      nutwrs = 72
      nutrst = 73

      ! 2. restart for passive tracer (input)
      ! -----------------------------

      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) ' *** PASSIVE TRACER MODEL OPTIONS'
      IF(lwp) WRITE(numout,*) ' *** CONTROL'
      IF(lwp) WRITE(numout,*) ' '

      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) ' *** restart option for passive tracer'
      IF(lwp) WRITE(numout,*) ' '

      IF(lrsttr) THEN
         IF(lwp) WRITE(numout,*) ' READ a restart FILE for passive tracer'
         IF(lwp) WRITE(numout,*) ' '
      ELSE
         IF(lwp) WRITE(numout,*) ' no restart FILE'
         IF(lwp) WRITE(numout,*) ' '

         ! 3. OPEN FILES for initial tracer value
         ! --------------------------------------
         DO jn=1,jptra

            ! OPEN input FILE only IF lutini(jn) is true
            ! ------------------------------------------
            IF (lutini(jn)) THEN  

               ! prepare input FILE name a
               ! -------------------------                        
               IF(lwp) WRITE(numout,*)  &
                  ' READ an initial FILE  for passive tracer number :',jn        &
                  ,' traceur : ',ctrcnm(jn) 
               IF(lwp) WRITE(numout,*) ' '
            END IF
         END DO
      ENDIF

      ! 4. Don't USE non penetrative convective mixing option
      !     it's not implemented for passive tracer
      ! -----------------------------------------------------

      IF( ln_zdfnpc) THEN
         IF(lwp) WRITE (numout,*) ' ===>>>> : w a r n i n g '
         IF(lwp) WRITE (numout,*) ' =======   ============= '
         IF(lwp) WRITE (numout,*) ' STOP, this sheme is not implemented'
         IF(lwp) WRITE (numout,*) ' in passive tracer model:'
         IF(lwp) WRITE (numout,*) ' non penetrative convect. mixing scheme'
         istop = istop + 1
      ENDIF

      ! 5. transport scheme option
      ! --------------------------

      IF(lwp) WRITE(numout,*) '  '
      CALL trc_trp_ctl


      ! 6. SMS model
      ! ---------------------------------------------

      IF(lwp) WRITE(numout,*) '  '
      IF(lwp) WRITE(numout,*) ' *** Source/Sink model option'
      IF(lwp) WRITE(numout,*) '  '


#if defined key_trc_lobster1
#   include "trcctl.lobster1.h90"
#elif defined key_trc_pisces
#   include "trcctl.pisces.h90"
#elif defined key_cfc
#   include "trcctl.cfc.h90"
#else

      IF(lwp) WRITE (numout,*) ' No Source/Sink model '
      IF(lwp) WRITE (numout,*) ' '
#endif

      ! E r r o r  control
      ! ------------------

      IF ( istop > 0  ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) istop,' E R R O R found : we stop'
         IF(lwp)WRITE(numout,*) '**************************'
         IF(lwp)WRITE(numout,*)
         STOP 'trcctl'
      ENDIF

   END SUBROUTINE trc_ctl

#else
   !!======================================================================
   !!  Empty module : No passive tracer
   !!======================================================================
CONTAINS
   SUBROUTINE trc_ctl

   END SUBROUTINE trc_ctl
   
#endif

END MODULE trcctl
