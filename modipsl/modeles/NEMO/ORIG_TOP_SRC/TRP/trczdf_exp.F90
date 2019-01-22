MODULE trczdf_exp
   !!==============================================================================
   !!                    ***  MODULE  trczdf_exp  ***
   !! Ocean passive tracers:  vertical component of the tracer mixing trend using
   !!                        an explicit time-stepping (time spllitting scheme)
   !!==============================================================================
#if defined key_passivetrc
   !!----------------------------------------------------------------------
   !!   trc_zdf_exp  : update the tracer trend with the vertical diffusion
   !!                  using an explicit time stepping
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc          ! ocean dynamics and active tracers variables
   USE trc              ! ocean passive tracers variables
   USE trctrp_lec       ! passive tracers transport
   USE prtctl_trc          ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_zdf_exp          ! routine called by step.F90

   !! * Module variable
   REAL(wp), DIMENSION(jpk) ::   &
      rdttrc                     ! vertical profile of 2 x tracer time-step

   !! * Substitutions
#  include "passivetrc_substitute.h90"
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trczdf_exp.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_zdf_exp( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_zdf_exp  ***
      !!                   
      !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
      !!      using an explicit time stepping and add it to the general trend 
      !!      of the tracer equations.
      !!
      !! ** Method  :   The vertical diffusion of tracers  is given by:
      !!         difft = dz( avt dz(trb) ) = 1/e3t dk+1( avt/e3w dk(trb) )
      !!      It is evaluated with an Euler scheme, using a time splitting
      !!      technique.
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend tra :
      !!          tra = tra + dz( avt dz(t) ) if lk_zdfddm= T)
      !!
      !! ** Action : - Update tra with the before vertical diffusion trend
      !!             - Save the trends  in trtrd ('key_trc_diatrd')
      !!
      !! History :
      !!   6.0  !  90-10  (B. Blanke)  Original code
      !!   7.0  !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)  correction on tracer trend loops
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  97-05  (G. Madec)  vertical component of isopycnal
      !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
      !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
      !!        !  00-05  (MA Foujols) add lbc for tracer trends
      !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
      !!        !                     avt multiple correction
      !!        !  00-08  (G. Madec)  double diffusive mixing
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
      
      !! * Local declarations
      INTEGER ::   ji, jj, jk, jl, jn             ! dummy loop indices
      REAL(wp) ::   &
         zlavmr,                 &  ! ???
         zave3r, ze3tr,          &  ! ???
         ztra                  !
      REAL(wp), DIMENSION(jpi,jpk) ::   &
         zwx, zwy
      CHARACTER (len=22) :: charout
      !!---------------------------------------------------------------------

      IF( kt == nittrc000 ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_zdf_exp : vertical tracer mixing'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      ! 0. Local constant initialization
      ! --------------------------------
      IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
         ! time step = 2 rdttra with Arakawa or TVD advection scheme
         IF( neuler == 0 .AND. kt == nittrc000 ) THEN
            rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)             ! restarting with Euler time stepping
         ELSEIF( kt <= nittrc000 + 1 ) THEN
            rdttrc(:) = 2. * rdttra(:) * FLOAT(ndttrc)         ! leapfrog
         ENDIF
      ELSE
         rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)      
      ENDIF


      zlavmr = 1. / FLOAT( n_trczdf_exp )

      DO jn = 1, jptra

         !                                                ! ===============
         DO jj = 2, jpjm1                                 !  Vertical slab
            !                                             ! ===============
            ! 1. Initializations
            ! ------------------

            ! Surface & bottom boundary conditions: no flux
            DO ji = 2, jpim1
               zwy(ji, 1 ) = 0.e0
               zwy(ji,jpk) = 0.e0
            END DO

            ! zwx and zwz arrays set to before tracer values
            DO jk = 1, jpk
               DO ji = 2, jpim1
                  zwx(ji,jk) = trb(ji,jj,jk,jn)
               END DO
            END DO


            ! 2. Time splitting loop
            ! ----------------------

            DO jl = 1, n_trczdf_exp

               ! first vertical derivative
               ! double diffusion: fstravs(ji,jj,jk) = avt(ji,jj,jk) /= avs (key_trc_zdfddm) 
               !                   fstravs(ji,jj,jk) = avs(ji,jj,jk) = avt
               DO jk = 2, jpk
                  DO ji = 2, jpim1
                     zave3r = 1.e0 / fse3w(ji,jj,jk) 
                     zwy(ji,jk) = fstravs(ji,jj,jk) * ( zwx(ji,jk-1) - zwx(ji,jk) ) * zave3r
                  END DO
               END DO


               ! trend estimation at kt+l*2*rdt/n_zdfexp
               DO jk = 1, jpkm1
                  DO ji = 2, jpim1
                     ze3tr = zlavmr / fse3t(ji,jj,jk)
                     ! 2nd vertical derivative
                     ztra = ( zwy(ji,jk) - zwy(ji,jk+1) ) * ze3tr
                     ! update the tracer trends
                     tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
                     ! update tracer fields at kt+l*2*rdt/n_trczdf_exp
                     zwx(ji,jk) = zwx(ji,jk) + rdttrc(jk) * ztra * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('zdf - exp')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
      ENDIF

   END SUBROUTINE trc_zdf_exp

#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_zdf_exp (kt )              ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_zdf_exp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_zdf_exp
#endif
   
   !!==============================================================================
END MODULE trczdf_exp
