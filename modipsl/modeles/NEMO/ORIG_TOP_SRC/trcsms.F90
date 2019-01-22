MODULE trcsms
   !!===========================================================================================
   !!
   !!                       *** MODULE trcsms ***
   !!
   !!  Time  loop of opa for passive tracer
   !!
   !!===========================================================================================
   !!  TOP 1.0,  LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/trcsms.F90,v 1.5 2006/04/10 15:40:29 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
#if defined key_passivetrc   
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   USE trcfreons
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC trc_sms

CONTAINS

   SUBROUTINE trc_sms( kt )
      !!===========================================================================================
      !!
      !!                       ROUTINE trcsms
      !!                     *****************
      !!
      !!  PURPOSE :
      !!  ---------
      !!          time loop of opa for passive tracer
      !!
      !!   METHOD :
      !!   -------
      !!      compute the well/spring evolution
      !!
      !!   INPUT :
      !!   -----
      !!      argument
      !!              ktask           : task identificator
      !!              kt              : time step
      !!      COMMON
      !!            all the COMMON defined in opa
      !!
      !!
      !!   OUTPUT :			: no
      !!   ------
      !!
      !!   WORKSPACE :
      !!   ---------
      !!
      !!   EXTERNAL :
      !!   --------
      !!      trcbio, trcsed, trcopt for NPZD or LOBSTER1 models
      !!
      !!      h3cprg for HAMOC3 and P3ZD models
      !!
      !!
      !!   History:
      !!   --------
      !!      original  : 96-11
      !!      additions : 99-07 (M. Levy)
      !!                  04-00 (O. Aumont, M.A. Foujols) HAMOCC3 and P3ZD
      !!                  12-00 (O. Aumont, E. Kestenare) add trcexp for instantaneous export 
      !!   05-03 (O. Aumont and A. El Moussaoui) F90
      !! -------------------------------------------------------------------------------------

      !! * Arguments
      !! -----------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index      

      !! * Local variables
      !! -----------------

      CHARACTER (len=25) :: charout

      !! this ROUTINE is called only every ndttrc time step
      !! --------------------------------------------------

      IF ( MOD(kt,ndttrc) /= 0) RETURN

      !! this first routines are parallelized on vertical slab
      !! ------------------------------------------------------

#if defined key_trc_lobster1

      !! tracers: optical model
      !! ----------------------

      CALL trcopt( kt)

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('OPT')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      !! tracers: biological model
      !! -------------------------

      CALL trcbio( kt)

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('BIO')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

      !! tracers: sedimentation model
      !! ----------------------------

      CALL trcsed(kt)
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('SED')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
 
      CALL trcexp(kt)

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('EXP')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF

#elif defined key_trc_pisces

      !! p4zprg: main PROGRAM for PISCES 
      !! -------------------------------
      CALL p4zprg(kt)

      !! SMS to DO

#elif defined key_cfc

      !! CFC's code taken from K. Rodgers

      !! This part is still experimental
      !! -------------------------------

      CALL trc_freons(kt)

#endif



   END SUBROUTINE trc_sms

#else
   !!======================================================================
   !!  Empty module : No passive tracer
   !!======================================================================
CONTAINS

   SUBROUTINE trc_sms( kt )

      ! no passive tracers
      INTEGER, INTENT( in ) ::   kt
      WRITE(*,*) 'trc_sms: You should not have seen this print! error?', kt
   END SUBROUTINE trc_sms

#endif 


END MODULE  trcsms
