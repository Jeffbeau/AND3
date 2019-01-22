MODULE icestp1d
   !!======================================================================
   !!                       ***  MODULE icestp1d   ***
   !!   Sea-Ice model : 1D LIM Sea ice model time-stepping
   !!======================================================================
#if defined key_cfg_1d && defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_cfg_1d'                                       1D Configuration
   !!   'key_ice_lim' :                                   Lim sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_stp_1d       : sea-ice model time-stepping
   !!----------------------------------------------------------------------
   USE dom_oce
   USE oce  ! dynamics and tracers variables
   USE in_out_manager
   USE ice_oce         ! ice variables
   USE flx_oce         ! forcings variables
   USE dom_ice
   USE cpl_oce
   USE blk_oce
   USE daymod
   USE phycst          ! Define parameters for the routines
   USE taumod
   USE ice
   USE iceini
   USE lbclnk
   USE limdyn
   USE limtrp
   USE limthd
   USE limflx
   USE limdia
   USE limwri
   USE limrst

   USE ocesbc
   USE flxmod
   USE flxrnf
   USE tradmp         ! damping salinity trend
   USE dtatem
   USE dtasal
   USE ocfzpt
   USE prtctl          ! Print control


   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC ice_stp_1d  ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!-----------------------------------------------------
   !!   LIM 2.0 , UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/C1D_SRC/icestp1d.F90,v 1.6 2006/04/19 14:43:13 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!-----------------------------------------------------

CONTAINS

   SUBROUTINE ice_stp_1d ( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE ice_stp_1d  ***
      !!                   
      !! ** Purpose :   Louvain la Neuve Sea Ice Model time stepping 
      !!
      !! ** Action  : - call the ice dynamics routine 
      !!              - call the ice advection/diffusion routine 
      !!              - call the ice thermodynamics routine 
      !!              - call the routine that computes mass and 
      !!                heat fluxes at the ice/ocean interface
      !!              - save the outputs 
      !!              - save the outputs for restart when necessary
      !!
      !! History :
      !!   1.0  !  99-11  (M. Imbard)  Original code
      !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
      !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
      !!   9.0  !  04-10  (C. Ethe) 1D configuration
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index

      !! * Local declarations
      INTEGER   ::   ji, jj   ! dummy loop indices

      REAL(wp) , DIMENSION(jpi,jpj)    :: &
         zsss_io, zsss2_io, zsss3_io          ! tempory workspaces
      REAL(wp)  :: ztair2
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ice_stp_1d : Louvain la Neuve Ice Model (LIM)' 
         IF(lwp) WRITE(numout,*) '~~~~~~~   forced case using bulk formulea'
         !  Initialize fluxes fields
         gtaux(:,:) = 0.e0
         gtauy(:,:) = 0.e0
      ENDIF

      ! Temperature , salinity and horizonta wind
      ! sst_io  and sss_io, u_io and v_io  are initialized at nit000 in limistate.F90 (or limrst.F90) with :
      !              sst_io = sst_io + (nfice - 1) * (tn(:,:,1)+rt0 )
      !              sss_io = sss_io + (nfice - 1) * sn(:,:,1)
      !              u_io   = u_io   + (nfice - 1) * un(:,:,1)
      !              v_io   = v_io   + (nfice - 1) * vn(:,:,1)
      !    cumulate fields
      !
      sst_io(:,:) = sst_io(:,:) + tn(:,:,1) + rt0
      sss_io(:,:) = sss_io(:,:) + sn(:,:,1)
      u_io  (:,:) = u_io  (:,:) + un(:,:,1)
      v_io  (:,:) = v_io  (:,:) + vn(:,:,1)
 

      
      IF( MOD( kt-1, nfice ) == 0 ) THEN
         
         ! The LIM model is going to be call
         sst_io(:,:) = sst_io(:,:) / FLOAT( nfice ) * tmask(:,:,1)
         sss_io(:,:) = sss_io(:,:) / FLOAT( nfice )
         u_io  (:,:) = u_io  (:,:) / FLOAT( nfice )
         v_io  (:,:) = v_io  (:,:) / FLOAT( nfice )
         gtaux (:,:) = taux  (:,:)
         gtauy (:,:) = tauy  (:,:)

         zsss_io (:,:) = SQRT( sss_io(:,:) ) 
         zsss2_io(:,:) =  sss_io(:,:) *  sss_io(:,:)
         zsss3_io(:,:) = zsss_io(:,:) * zsss_io(:,:) * zsss_io(:,:)

         DO jj = 1, jpj
            DO ji = 1, jpi
               tfu(ji,jj)  = ABS ( rt0 - 0.0575       *   sss_io(ji,jj)   &
                  &                    + 1.710523e-03 * zsss3_io(ji,jj)   &
                  &                    - 2.154996e-04 * zsss2_io(ji,jj) ) * tms(ji,jj)
            END DO
         END DO
         
         
         IF(ln_ctl) THEN         ! print mean trends (used for debugging)
            CALL prt_ctl_info('Ice Forcings ')
            CALL prt_ctl(tab2d_1=qsr_oce ,clinfo1=' qsr_oce  : ', tab2d_2=qsr_ice , clinfo2=' qsr_ice   : ')
            CALL prt_ctl(tab2d_1=qnsr_oce,clinfo1=' qnsr_oce : ', tab2d_2=qnsr_ice, clinfo2=' qnsr_ice  : ')
            CALL prt_ctl(tab2d_1=evap    ,clinfo1=' evap     : ')
            CALL prt_ctl(tab2d_1=tprecip ,clinfo1=' precip   : ', tab2d_2=sprecip , clinfo2=' Snow      : ')
            CALL prt_ctl(tab2d_1=gtaux   ,clinfo1=' u-stress : ', tab2d_2=gtauy   , clinfo2=' v-stress  : ')
            CALL prt_ctl(tab2d_1=sst_io  ,clinfo1=' sst      : ', tab2d_2=sss_io  , clinfo2=' sss       : ')
            CALL prt_ctl(tab2d_1=u_io    ,clinfo1=' u_io     : ', tab2d_2=v_io    , clinfo2=' v_io      : ')
            CALL prt_ctl(tab2d_1=hsnif   ,clinfo1=' hsnif  1 : ', tab2d_2=hicif   , clinfo2=' hicif     : ')
            CALL prt_ctl(tab2d_1=frld    ,clinfo1=' frld   1 : ', tab2d_2=sist    , clinfo2=' sist      : ')
         ENDIF

         ! Ice model call
         numit = numit + nfice                                       ! Friction velocity
                                                                      
         DO jj = 1, jpj
            DO ji = 1, jpi
               tio_u(ji,jj) = - gtaux(ji,jj) / rau0
               tio_v(ji,jj) = - gtauy(ji,jj) / rau0               
               ztair2       = gtaux(ji,jj) * gtaux(ji,jj) + gtauy(ji,jj) * gtauy(ji,jj)           
               ust2s(ji,jj) = ( SQRT( ztair2  )  / rau0 )  * tms(ji,jj)  
            END DO
         END DO

         !                                                           !--------------------!
         CALL lim_thd                                                ! Ice thermodynamics !
         !                                                           !--------------------!
         IF(ln_ctl) THEN
            CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  2 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
            CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   2 : ', tab2d_2=sist  , clinfo2=' sist      : ')
            CALL prt_ctl(tab2d_1=u_io  ,clinfo1=' u_io   4 : ', tab2d_2=v_io  , clinfo2=' v_io      : ')
            CALL prt_ctl(tab2d_1=tio_u  ,clinfo1=' tio_u  4 : ', tab2d_2=tio_v  , clinfo2=' tio_v     : ')
         ENDIF



         ! Mass and heat fluxes from ice to ocean
         !                                                           !------------------------------!
         CALL lim_flx                                                ! Ice/Ocean Mass & Heat fluxes !
         !                                                           !------------------------------!

         IF(ln_ctl) THEN
            CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  7 : ', tab2d_2=hicif , clinfo2=' hicif   : ')
            CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   7 : ', tab2d_2=sist  , clinfo2=' sist      : ')
            CALL prt_ctl(tab2d_1=tio_u  ,clinfo1=' tio_u  7 : ', tab2d_2=tio_v  , clinfo2=' tio_v     : ')
         ENDIF
         !                                                           !-------------!
         CALL lim_wri                                                ! Ice outputs !
         !                                                           !-------------!

         IF( MOD( numit, nstock ) == 0 .OR. numit == nlast ) THEN
            !                                                        !------------------!
            CALL lim_rst_write( numit )                              ! Ice restart file !
            !                                                        !------------------!
         ENDIF

         ! Re-initialization of forcings
         qsr_oce (:,:) = 0.e0
         qsr_ice (:,:) = 0.e0
         qnsr_oce(:,:) = 0.e0
         qnsr_ice(:,:) = 0.e0 
         dqns_ice(:,:) = 0.e0 
         tprecip (:,:) = 0.e0 
         sprecip (:,:) = 0.e0
         qla_ice (:,:) = 0.e0
         dqla_ice(:,:) = 0.e0
         fr1_i0  (:,:) = 0.e0
         fr2_i0  (:,:) = 0.e0
         evap    (:,:) = 0.e0

        CALL oce_sbc_1d ( kt )

      ENDIF

   END SUBROUTINE ice_stp_1d
   
   
   SUBROUTINE oce_sbc_1d( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE oce_sbc_1d  ***
      !!                    
      !! ** Purpose : - Ocean surface boundary conditions with LIM sea-ice
      !!        model in forced mode using bulk formulea
      !!
      !! History :
      !!   1.0  !  99-11  (M. Imbard)  Original code
      !!        !  01-03  (D. Ludicone, E. Durand, G. Madec) free surf.
      !!   2.0  !  02-09  (G. Madec, C. Ethe)  F90: Free form and module
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * arguments
      INTEGER, INTENT( in  ) ::   kt   ! ocean time step

      !! * Local declarations
      INTEGER  ::   ji, jj                   ! dummy loop indices
      REAL(wp) ::   ztxy
      !!----------------------------------------------------------------------

      ! 1. initialization to zero at kt = nit000
      ! ---------------------------------------
      
      IF( kt == nit000 ) THEN     
         qsr    (:,:) = 0.e0
         qt     (:,:) = 0.e0
         qrp    (:,:) = 0.e0
         emp    (:,:) = 0.e0
         emps   (:,:) = 0.e0
         erp    (:,:) = 0.e0
#if ! defined key_dynspg_rl 
         dmp    (:,:) = 0.e0
#endif
      ENDIF

      CALL oce_sbc_dmp       ! Computation of internal and evaporation damping terms       

      ! Surface Ocean fluxes
      ! ====================
      
      ! Surface heat flux (W/m2)
      ! -----------------
      
      qt  (:,:) = fnsolar(:,:) + fsolar(:,:)     ! non solar heat flux + solar flux
      qsr (:,:) = fsolar(:,:)                     ! solar flux
      
#if ! defined key_dynspg_rl     
      ! total concentration/dilution effect (use on SSS)
      emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:) + empold
      
      ! total volume flux (use on sea-surface height)
      emp (:,:) = fmass(:,:) -   dmp(:,:) + runoff(:,:) + erp(:,:) + empold      
#else
      ! Rigid-lid (emp=emps=E-P-R+Erp)
      emps(:,:) = fmass(:,:) + fsalt(:,:) + runoff(:,:) + erp(:,:)     ! freshwater flux
      emp (:,:) = emps(:,:)
      
#endif
      
      ! Surface stress
      ! --------------
      
      ! update the stress beloww sea-ice area
      DO jj = 1, jpjm1
         DO ji = 1, fs_jpim1   ! vertor opt.
            ztxy        = freezn(ji,jj)             ! ice/ocean indicator at T-points
            taux(ji,jj) = (1.-ztxy) * taux(ji,jj) + ztxy * ftaux(ji,jj)    ! stress at the ocean surface
            tauy(ji,jj) = (1.-ztxy) * tauy(ji,jj) + ztxy * ftauy(ji,jj)
         END DO
      END DO
      
      ! boundary condition on the stress (taux,tauy)
      CALL lbc_lnk( taux, 'U', -1. )
      CALL lbc_lnk( tauy, 'V', -1. )
      
      ! Re-initialization of fluxes
      sst_io(:,:) = 0.e0
      sss_io(:,:) = 0.e0
      u_io  (:,:) = 0.e0
      v_io  (:,:) = 0.e0
      
      
   END SUBROUTINE oce_sbc_1d
   
#if defined key_dtasal
   !!----------------------------------------------------------------------
   !!   'key_dtasal'                                          salinity data
   !!----------------------------------------------------------------------
   SUBROUTINE oce_sbc_dmp
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE oce_sbc_dmp  ***
      !!                    
      !! ** Purpose : Computation of internal and evaporation damping terms 
      !!        for ocean surface boundary conditions 
      !!
      !! History :
      !!   9.0  !  04-01  (G. Madec, C. Ethe)  Original code
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj                   ! dummy loop indices
      REAL(wp) ::   zerp, zsrp
#if ! defined key_dynspg_rl
      REAL(wp) ::   zwei
      REAL(wp) ::   zerpplus(jpi,jpj), zerpminus(jpi,jpj)
      REAL(wp) ::   zplus, zminus, zadefi
# if defined key_tradmp
      INTEGER jk
      REAL(wp), DIMENSION(jpi,jpj) ::   zstrdmp
# endif
#endif
      !!----------------------------------------------------------------------


      ! sea ice indicator (1 or 0)
      DO jj = 1, jpj
         DO ji = 1, jpi
            freezn(ji,jj) = MAX(0., SIGN(1., freeze(ji,jj)-rsmall) )
         END DO
      END DO

      ! Initialisation
      ! --------------
      ! Restoring coefficients on SST and SSS   
      zsrp = dqdt0 * ro0cpr * rauw   ! (Kg/m2/s) 

#if ! defined key_dynspg_rl 
      ! Free-surface
         
      ! Internal damping
# if defined key_tradmp
      ! Vertical mean of dampind trend (computed in tradmp module)
      zstrdmp(:,:) = 0.e0
      DO jk = 1, jpk
         zstrdmp(:,:) = zstrdmp(:,:) + strdmp(:,:,jk) * fse3t(:,:,jk)
      END DO
      ! volume flux associated to internal damping to climatology
      dmp(:,:) = zstrdmp(:,:) * rauw / ( sss_io(:,:) + 1.e-20 )
# else
      dmp(:,:) = 0.e0            ! No internal damping
# endif
      
      !   evaporation damping term ( Surface restoring )
      zerpplus (:,:) = 0.e0
      zerpminus(:,:) = 0.e0
      zplus          =  15. / rday
      zminus         = -15. / rday
      
      DO jj = 1, jpj
         DO ji = 1, jpi
            zerp = ( 1. - 2.*upsrnfh(ji,jj) ) * zsrp   &
               & * ( sss_io(ji,jj) - s_dta(ji,jj,1) )     &
               & / ( sss_io(ji,jj) + 1.e-20        )
            erp(ji,jj) = zerp
            zerpplus (ji,jj) = MAX( erp(ji,jj), 0.e0 )
            zerpminus(ji,jj) = MIN( erp(ji,jj), 0.e0 )
         END DO
      END DO

      aplus  = 0.e0
      aminus = 0.e0
      DO jj = 1, jpj
         DO ji = 1, jpi
            zwei   = e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            aplus  = aplus  + zerpplus (ji,jj) * zwei
            aminus = aminus - zerpminus(ji,jj) * zwei
         END DO
      END DO

      IF(ln_ctl .AND. lwp)   WRITE(numout,*) ' oce_sbc_dmp : a+ = ', aplus, ' a- = ', aminus
#else
      ! Rigid-lid (emp=emps=E-P-R+Erp)
      
      erp(:,:) = ( 1. - freezn(:,:) ) * zsrp    &   ! surface restoring term
         &     * ( sss_io(:,:) - s_dta(:,:,1) )     &
         &     / ( sss_io(:,:) + 1.e-20      )
#endif

   END SUBROUTINE oce_sbc_dmp
   
#else
   !!----------------------------------------------------------------------
   !!   Dummy routine                                      NO salinity data
   !!----------------------------------------------------------------------
   USE in_out_manager
   SUBROUTINE oce_sbc_dmp         ! Dummy routine
      if(lwp) WRITE(numout,*) 'oce_sbc_dmp: you should not have seen that print! error?'
   END SUBROUTINE oce_sbc_dmp
#endif
#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module   NO 1D && NO LIM sea-ice model
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE ice_stp_1d ( kt )     ! Dummy routine
      if(lwp) WRITE(numout,*) 'ice_stp_1d: You should not have seen this print! error?', kt
   END SUBROUTINE ice_stp_1d
   SUBROUTINE oce_sbc_1d ( kt )     ! Dummy routine
      if(lwp)WRITE(numout,*) 'oce_sbc_1d: You should not have seen this print! error?', kt
   END SUBROUTINE oce_sbc_1d
#endif   
   !!======================================================================
END MODULE icestp1d
