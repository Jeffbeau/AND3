MODULE traqsr
   !!======================================================================
   !!                       ***  MODULE  traqsr  ***
   !! Ocean physics: solar radiation penetration in the top ocean levels
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   tra_qsr      : trend due to the solar radiation penetration
   !!   tra_qsr_init : solar radiation penetration initialization
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE trdmod          ! ocean active tracers trends 
   USE trdmod_oce      ! ocean variables trends
   USE in_out_manager  ! I/O manager
   USE trc_oce         ! share SMS/Ocean variables
   USE ocesbc          ! thermohaline fluxes
   USE phycst          ! physical constants
   USE prtctl          ! Print control

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_qsr      ! routine called by step.F90 (ln_traqsr=T)
   PUBLIC tra_qsr_init ! routine called by opa.F90

   !! * Shared module variables
   LOGICAL, PUBLIC ::   ln_traqsr = .TRUE.   !: qsr flag (Default=T)

   !! * Module variables
   REAL(wp), PUBLIC ::  & !!! * penetrative solar radiation namelist *
      rabs = 0.58_wp,   &  ! fraction associated with xsi1
      xsi1 = 0.35_wp,   &  ! first depth of extinction 
      xsi2 = 23.0_wp       ! second depth of extinction 
      !                    ! (default values: water type Ib)
   LOGICAL ::           & 
      ln_qsr_sms = .false. ! flag to use or not the biological 
      !                    ! fluxes for light 
   
   INTEGER ::    &
      nksr                 ! number of levels
   REAL(wp), DIMENSION(jpk) ::   &
      gdsr                 ! profile of the solar flux penetration

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traqsr.F90,v 1.10 2005/09/22 10:55:46 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_qsr( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr  ***
      !!
      !! ** Purpose :   Compute the temperature trend due to the solar radiation
      !!      penetration and add it to the general temperature trend.
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is
      !!      defined through two penetration length scale (xsr1,xsr2) and a 
      !!      ratio (rabs) as :
      !!         I(k) = Qsr*( rabs*EXP(z(k)/xsr1) + (1.-rabs)*EXP(z(k)/xsr2) )
      !!         The temperature trend associated with the solar radiation
      !!      penetration is given by :
      !!            zta = 1/e3t dk[ I ] / (rau0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         In z-coordinate case, the computation is only done down to the
      !!      level where I(k) < 1.e-15 W/m2. In addition, the coefficients 
      !!      used for the computation are calculated one for once as they
      !!      depends on k only.
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - save the trend in ttrd ('key_trdtra')
      !!
      !! History :
      !!   6.0  !  90-10  (B. Blanke)  Original code
      !!   7.0  !  91-11  (G. Madec)
      !!        !  96-01  (G. Madec)  s-coordinates
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!----------------------------------------------------------------------
      !! * Modules used     
      USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
                         ztdsa => va      ! use va as 3D workspace   

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step

      !! * Local declarations
      INTEGER ::    ji, jj, jk            ! dummy loop indexes
      REAL(wp) ::   zc0, zta              ! temporary scalars
      REAL(wp) ::   zc1 , zc2 ,        &  ! temporary scalars
                    zdp1, zdp2            !
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF ( lwp )  WRITE(numout,*)
         IF ( lwp )  WRITE(numout,*) 'tra_qsr : penetration of the surface solar radiation'
         IF ( lwp )  WRITE(numout,*) '~~~~~~~'
      ENDIF

      ! Save ta and sa trends
      IF( l_trdtra )   THEN
         ztdta(:,:,:) = ta(:,:,:) 
         ztdsa(:,:,:) = 0.e0
      ENDIF

      IF( lk_qsr_sms .AND. ln_qsr_sms ) THEN    !  Biological fluxes  !
         !                                      ! =================== !
         !
         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  
                  zc0 = ro0cpr  / fse3t(ji,jj,jk)      ! compute the qsr trend
                  zta = zc0 * ( etot3(ji,jj,jk  ) * tmask(ji,jj,jk)     &
                     &        - etot3(ji,jj,jk+1) * tmask(ji,jj,jk+1) )
                  
                  ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
                  
               END DO
            END DO
            !                                             ! ===============
         END DO                                           !   End of slab
         !                                                ! ===============
         ! save the trends for diagnostic
         ! qsr tracers trends
         IF( l_trdtra )   THEN
            ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
            CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
         ENDIF

      ELSE
         !                                                ! =================== !
         IF( lk_sco ) THEN                                !    s-coordinate     !
         !                                                ! =================== !
         !
         !                                                   ! ===============
            DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                                ! ===============
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.

                     zdp1 = -fsdepw(ji,jj,jk  )              ! compute the qsr trend
                     zdp2 = -fsdepw(ji,jj,jk+1)
                     zc0  = qsr(ji,jj) * ro0cpr / fse3t(ji,jj,jk)
                     zc1  =   (  rabs * EXP(zdp1/xsi1) + (1.-rabs) * EXP(zdp1/xsi2)  )
                     zc2  = - (  rabs * EXP(zdp2/xsi1) + (1.-rabs) * EXP(zdp2/xsi2)  )
                     zta  = zc0 * (  zc1 * tmask(ji,jj,jk) + zc2 * tmask(ji,jj,jk+1)  )
                     
                     ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
                     
                  END DO
               END DO
               !                                                ! ===============
            END DO                                              !   End of slab
            !                                                   ! ===============
            ! save the trends for diagnostic
            ! qsr tracers trends
            IF( l_trdtra )   THEN
               ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
               CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
            ENDIF
            !
         ENDIF
         !                                                ! =================== !
         IF( lk_zps ) THEN                                !    partial steps    !
            !                                             ! =================== !
            !
            !                                                ! ===============
            DO jk = 1, nksr                                  ! Horizontal slab
               !                                             ! ===============
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     
                     zc0 = qsr(ji,jj) / fse3t(ji,jj,jk)      ! compute the qsr trend
                     zta = zc0 * ( gdsr(jk) * tmask(ji,jj,jk) - gdsr(jk+1) * tmask(ji,jj,jk+1) )
                     
                     ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
                     
                  END DO
               END DO
               !                                             ! ===============
            END DO                                           !   End of slab
            !                                                ! ===============
            ! save the trends for diagnostic
            ! qsr tracers trends
            IF( l_trdtra )   THEN
               ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
               CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
            ENDIF
            !
         ENDIF
         !                                                ! =================== !
         IF( lk_zco ) THEN                                !     z-coordinate    !
            !                                             ! =================== !
            !
            !                                                ! ===============
            DO jk = 1, nksr                                  ! Horizontal slab
               !                                             ! ===============
               zc0 = 1. / fse3t(1,1,jk)
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1   ! vector opt.
                     !                                       ! compute qsr forcing trend
                     zta = qsr(ji,jj) * zc0 * ( gdsr(jk)*tmask(ji,jj,jk) - gdsr(jk+1)*tmask(ji,jj,jk+1) )
                     
                     ta(ji,jj,jk) = ta(ji,jj,jk) + zta       ! add qsr trend to the temperature trend
                     
                  END DO
               END DO
               !                                             ! ===============
            END DO                                           !   End of slab
            !                                                ! ===============
            ! save the trends for diagnostic
            ! qsr tracers trends
            IF( l_trdtra )   THEN
               ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
               CALL trd_mod(ztdta, ztdsa, jpttdqsr, 'TRA', kt)
            ENDIF
            !
         ENDIF
         !
      ENDIF


      IF(ln_ctl) THEN         ! print mean trends (used for debugging)
         CALL prt_ctl(tab3d_1=ta, clinfo1=' qsr  - Ta: ', mask1=tmask, clinfo3='tra-ta')
      ENDIF

   END SUBROUTINE tra_qsr


   SUBROUTINE tra_qsr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr_init  ***
      !!
      !! ** Purpose :   Initialization for the penetrative solar radiation
      !!
      !! ** Method  :   The profile of solar radiation within the ocean is set
      !!      from two length scale of penetration (xsr1,xsr2) and a ratio
      !!      (rabs). These parameters are read in the namqsr namelist. The
      !!      default values correspond to clear water (type I in Jerlov' 
      !!      (1968) classification.
      !!         called by tra_qsr at the first timestep (nit000)
      !!
      !! ** Action  : - initialize xsr1, xsr2 and rabs
      !!
      !! Reference :
      !!   Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!
      !! History :
      !!   8.5  !  02-06  (G. Madec) Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::    ji,jj,jk, &  ! dummy loop index
                    indic        ! temporary integer
      REAL(wp) ::   zdp1         ! temporary scalar

      NAMELIST/namqsr/ ln_traqsr, rabs, xsi1, xsi2, ln_qsr_sms
      !!----------------------------------------------------------------------

      ! Read Namelist namqsr : ratio and length of penetration
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, namqsr )

      ! Parameter control and print
      ! ---------------------------
      IF( ln_traqsr  ) THEN
        IF ( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_qsr_init : penetration of the surface solar radiation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '    Namelist namqsr : set the parameter of penetration'
         WRITE(numout,*) '        fraction associated with xsi     rabs        = ',rabs
         WRITE(numout,*) '        first depth of extinction        xsi1        = ',xsi1
         WRITE(numout,*) '        second depth of extinction       xsi2        = ',xsi2
         IF( lk_qsr_sms ) THEN
            WRITE(numout,*) '     Biological fluxes for light(Y/N) ln_qsr_sms  = ',ln_qsr_sms
         ENDIF
         WRITE(numout,*) ' '
        END IF
      ELSE
        IF ( lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_qsr_init : NO solar flux penetration'
         WRITE(numout,*) '~~~~~~~~~~~~'
        END IF
      ENDIF

      IF( rabs > 1.e0 .OR. rabs < 0.e0 .OR. xsi1 < 0.e0 .OR. xsi2 < 0.e0 ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '             0<rabs<1, 0<xsi1, or 0<xsi2 not satisfied'
         nstop = nstop + 1
      ENDIF


      ! Initialization
      ! --------------
      IF( .NOT. lk_sco ) THEN
         ! z-coordinate with or without partial step : same before last ocean w-level everywhere
         gdsr(:) = 0.e0
         DO jk = 1, jpk
            zdp1 = -fsdepw(1,1,jk)
            gdsr(jk) = ro0cpr * (  rabs  * EXP( zdp1/xsi1 ) + (1.-rabs) * EXP( zdp1/xsi2 )  )
            IF ( gdsr(jk) <= 1.e-10 ) EXIT
         END DO
         indic = 0
         DO jk = 1, jpk
            IF( gdsr(jk) <= 1.e-15 .AND. indic == 0 ) THEN
               gdsr(jk) = 0.e0
               nksr = jk
               !!bug Edmee chg res   nksr = jk - 1
               indic = 1
            ENDIF
         END DO
         nksr = MIN( nksr, jpkm1 )
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) '        - z-coordinate, level max of computation =', nksr
            WRITE(numout,*) '             profile of coef. of penetration:'
            WRITE(numout,"('              ',7e11.2)") ( gdsr(jk), jk = 1, nksr )
            WRITE(numout,*)
         ENDIF
         ! Initialisation of Biological fluxes for light here because
         ! the optical biological model is call after the dynamical one
         IF( lk_qsr_sms .AND. ln_qsr_sms ) THEN
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     etot3(ji,jj,jk) = qsr(ji,jj) * gdsr(jk) * tmask(ji,jj,jk) / ro0cpr
                  END DO
               END DO
            END DO
         ENDIF

      ENDIF

   END SUBROUTINE tra_qsr_init

   !!======================================================================
END MODULE traqsr
