MODULE trcnxt
   !!======================================================================
   !!                       ***  MODULE  trcnxt  ***
   !! Ocean passive tracers:  time stepping on passives tracers
   !!======================================================================
#if defined key_passivetrc   
   !!----------------------------------------------------------------------
   !!   trc_nxt     : time stepping on passive tracers
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce_trc         ! ocean dynamics and tracers variables
   USE trc             ! ocean passive tracers variables
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE trctrp_lec      ! pasive tracers transport
   USE prtctl_trc      ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC trc_nxt          ! routine called by step.F90
   !!----------------------------------------------------------------------
   !!   TOP 1.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcnxt.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE trc_nxt( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trcnxt  ***
      !!
      !! ** Purpose :   Compute the passive tracers fields at the 
      !!      next time-step from their temporal trends and swap the fields.
      !! 
      !! ** Method  :   Apply lateral boundary conditions on (ua,va) through 
      !!      call to lbc_lnk routine
      !!   default:
      !!      arrays swap
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!         (trb) = (trn) 
      !!
      !!   For Arakawa or TVD Scheme : 
      !!      A Asselin time filter applied on now tracers (trn) to avoid
      !!      the divergence of two consecutive time-steps and tr arrays
      !!      to prepare the next time_step:
      !!         (trb) = (trn) + atfp [ (trb) + (tra) - 2 (trn) ]
      !!         (trn) = (tra) ; (tra) = (0,0)
      !!
      !!
      !! ** Action  : - update trb, trn
      !!
      !! History :
      !!   7.0  !  91-11  (G. Madec)  Original code
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  95-02  (M. Levy)   passive tracers 
      !!        !  96-02  (G. Madec & M. Imbard)  opa release 8.0
      !!   8.0  !  96-04  (A. Weaver)  Euler forward step
      !!   8.2  !  99-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!        !  02-11  (C. Talandier, A-M Treguier) Open boundaries
      !!   9.0  !  04-03  (C. Ethe) passive tracers 
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      !! * Local declarations
      INTEGER ::   ji, jj, jk,jn   ! dummy loop indices
      REAL(wp) ::   zfact, ztra    ! temporary scalar
      CHARACTER (len=22) :: charout
      !!----------------------------------------------------------------------

      IF( kt == nittrc000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_nxt : time stepping on passive tracers'
      ENDIF


      DO jn = 1, jptra

         ! 0. Lateral boundary conditions on tra (T-point, unchanged sign)
         ! ---------------------------------============
         CALL lbc_lnk( tra(:,:,:,jn), 'T', 1. )   
         
         !                                                ! ===============
         DO jk = 1, jpk                                   ! Horizontal slab
            !                                             ! ===============
            ! 1. Leap-frog scheme (only in explicit case, otherwise the 
            ! -------------------  time stepping is already done in trczdf)
            IF( l_trczdf_exp .AND. ( ln_trcadv_cen2 .OR. ln_trcadv_tvd) ) THEN
               zfact = 2. * rdttra(jk) * FLOAT(ndttrc) 
               IF( neuler == 0 .AND. kt == nittrc000 ) zfact = rdttra(jk) * FLOAT(ndttrc) 
               tra(:,:,jk,jn) = ( trb(:,:,jk,jn) + zfact * tra(:,:,jk,jn) ) * tmask(:,:,jk)
            ENDIF

         END DO

#if defined key_obc
        IF(lwp) WRITE(numout,cform_err)
        IF(lwp) WRITE(numout,*) '          Passive tracers and Open Boundary condition can not be used together '
        IF(lwp) WRITE(numout,*) '          Check in trc_nxt routine'
        nstop = nstop + 1
#endif

         DO jk = 1, jpk  

            ! 2. Time filter and swap of arrays
            ! ---------------------------------
            IF ( ln_trcadv_cen2 .OR. ln_trcadv_tvd  ) THEN

               IF( neuler == 0 .AND. kt == nittrc000 ) THEN
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        trb(ji,jj,jk,jn) = trn(ji,jj,jk,jn)
                        trn(ji,jj,jk,jn) = tra(ji,jj,jk,jn)
                        tra(ji,jj,jk,jn) = 0.
                     END DO
                  END DO
               ELSE
                  DO jj = 1, jpj
                     DO ji = 1, jpi
                        trb(ji,jj,jk,jn) = atfp  * ( trb(ji,jj,jk,jn) + tra(ji,jj,jk,jn) ) + atfp1 * trn(ji,jj,jk,jn)
                        trn(ji,jj,jk,jn) = tra(ji,jj,jk,jn)
                        tra(ji,jj,jk,jn) = 0.
                     END DO
                  END DO
               ENDIF

            ELSE
!  case of smolar scheme or muscl
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     trb(ji,jj,jk,jn) = tra(ji,jj,jk,jn)
                     trn(ji,jj,jk,jn) = tra(ji,jj,jk,jn)
                     tra(ji,jj,jk,jn) = 0.
                  END DO
               END DO

            ENDIF
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nxt')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=trn, mask=tmask, clinfo=ctrcnm)
      ENDIF


   END SUBROUTINE trc_nxt

#else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trc_nxt( kt )  
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'trc_nxt: You should not have seen this print! error?', kt
   END SUBROUTINE trc_nxt
#endif
   !!======================================================================
END MODULE trcnxt
