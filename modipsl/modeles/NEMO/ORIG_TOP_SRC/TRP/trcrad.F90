MODULE trcrad
   !!======================================================================
   !!                       ***  MODULE  trcrad  ***
   !! Ocean passive tracers:  correction of negative concentrations
   !!======================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_rad    : correction of negative concentrations
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc             ! ocean dynamics and tracers variables
   USE trc                 ! ocean passive tracers variables
   USE lib_mpp
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_rad        ! routine called by trcstp.F90
   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcrad.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_rad( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_rad  ***
      !!
      !! ** Purpose : "crappy" routine to correct artificial negative
      !!      concentrations due to isopycnal scheme
      !!
      !! ** Method  : Set negative concentrations to zero
      !!              compute the corresponding mass added to the tracers
      !!              and remove it when possible 
      !!
      !! History :
      !!   8.2  !  01-01  (O. Aumont & E. Kestenare)  Original code
      !!   9.0  !  04-03  (C. Ethe)  free form F90
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      
      !! * Local declarations
      INTEGER ::  ji, jj, jk, jn             ! dummy loop indices
#if defined key_trc_pisces || defined key_trc_lobster1
      REAL(wp) :: zvolk, trcorb, trmasb ,trcorn, trmasn  
#endif
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_rad : Correct artificial negative concentrations '
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF


#if defined key_cfc
      DO jn = 1, jptra
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  trn(ji,jj,jk,jn) = MAX( 0. , trn(ji,jj,jk,jn) )
                  trb(ji,jj,jk,jn) = MAX( 0. , trb(ji,jj,jk,jn) )
               END DO
            END DO
         END DO
      END DO
      
#elif defined key_trc_pisces || defined key_trc_lobster1

      DO jn = 1, jptra
         trcorb = 0.
         trmasb = 0.
         trcorn = 0.
         trmasn = 0.
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zvolk = e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) &
#if defined key_off_degrad
                  &  * facvol(ji,jj,jk) &
#endif
                  &  * tmask(ji,jj,jk) * tmask_i(ji,jj)

                  trcorb = trcorb + MIN( 0., trb(ji,jj,jk,jn) )  * zvolk
                  trcorn = trcorn + MIN( 0., trn(ji,jj,jk,jn) )  * zvolk

                  trb(ji,jj,jk,jn) = MAX( 0. , trb(ji,jj,jk,jn) )
                  trn(ji,jj,jk,jn) = MAX( 0. , trn(ji,jj,jk,jn) )

                  trmasb = trmasb + trb(ji,jj,jk,jn) * zvolk
                  trmasn = trmasn + trn(ji,jj,jk,jn) * zvolk
               END DO
            END DO
         END DO

         IF( lk_mpp ) THEN
           CALL mpp_sum( trcorb )   ! sum over the global domain
           CALL mpp_sum( trcorn )   ! sum over the global domain
           CALL mpp_sum( trmasb )   ! sum over the global domain
           CALL mpp_sum( trmasn )   ! sum over the global domain
         ENDIF

         IF( trcorb /= 0 ) THEN
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     trb(ji,jj,jk,jn) = MAX( 0. , trb(ji,jj,jk,jn) )
                     trb(ji,jj,jk,jn) = trb(ji,jj,jk,jn) * ( 1. + trcorb/trmasb ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF

         IF( trcorn /= 0) THEN
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     trn(ji,jj,jk,jn) = MAX( 0. , trn(ji,jj,jk,jn) )
                     trn(ji,jj,jk,jn) = trn(ji,jj,jk,jn) * ( 1. + trcorn/trmasn ) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF

      END DO
     
#endif

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rad')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
      ENDIF

      
   END SUBROUTINE trc_rad

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_rad (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_rad: You should not have seen this print! error?', kt
   END SUBROUTINE trc_rad
#endif
   
   !!======================================================================
END MODULE trcrad
