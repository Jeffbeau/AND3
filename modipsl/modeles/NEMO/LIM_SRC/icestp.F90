MODULE icestp
   !!======================================================================
   !!                       ***  MODULE icestp   ***
   !!   Sea-Ice model : LIM Sea ice model time-stepping
   !!======================================================================
#if defined key_ice_lim
   !!----------------------------------------------------------------------
   !!   'key_ice_lim' :                                   Lim sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_stp       : sea-ice model time-stepping
   !!----------------------------------------------------------------------
   USE dom_oce
   USE oce  ! dynamics and tracers variables
   USE in_out_manager
   USE ice_oce         ! ice variables
   USE flx_oce         ! forcings variables
   USE dom_ice
   USE cpl_oce
   USE daymod
   USE phycst          ! Define parameters for the routines
   USE taumod
   USE ice
   USE iceini
   USE ocesbc
   USE lbclnk
   USE limdyn
   USE limtrp
   USE limthd
   USE limflx
   USE limdia
   USE limwri
   USE limrst
   USE limdmp          ! Ice damping
   USE prtctl          ! Print control

!!DB
   USE lib_ncdf        ! netCDF I/O library
   USE restart
!DBG
#if defined key_BGCM_01
!   USE bgcm_01_model, ONLY : sw_rad
#endif


   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC ice_stp  ! called by step.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!-----------------------------------------------------
   !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/icestp.F90,v 1.7 2006/03/21 08:42:22 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!-----------------------------------------------------

CONTAINS

  SUBROUTINE ice_stp ( kt )
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ice_stp  ***
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
    !!----------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT( in ) ::   kt         ! ocean time-step index

    !! * Local declarations
    INTEGER   ::   ji, jj   ! dummy loop indices

    REAL(wp) , DIMENSION(jpi,jpj)    :: &
         zsss_io, zsss2_io, zsss3_io          ! tempory workspaces
    !!----------------------------------------------------------------------

    IF( kt == nit000 ) THEN
       IF( lk_cpl ) THEN
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'ice_stp : Louvain la Neuve Ice Model (LIM)'
          IF(lwp) WRITE(numout,*) '~~~~~~~   coupled case'
       ELSE
          IF(lwp) WRITE(numout,*)
          IF(lwp) WRITE(numout,*) 'ice_stp : Louvain la Neuve Ice Model (LIM)' 
          IF(lwp) WRITE(numout,*) '~~~~~~~   forced case using bulk formulea'
       ENDIF
       !  Initialize fluxes fields
       gtaux(:,:) = 0.e0
       gtauy(:,:) = 0.e0
    ENDIF

    ! Temperature , salinity and horizonta wind
    ! sst_io  and sss_io, u_io and v_io  are initialized at nit000 in limistate.F90 (or limrst.F90) with :
    !              sst_io = sst_io + (nfice - 1) * (tn(:,:,1)+rt0 )
    !              sss_io = sss_io + (nfice - 1) * sn(:,:,1)
    !              u_io  = u_io  + (nfice - 1) * 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
    !              v_io  = v_io  + (nfice - 1) * 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
    !    cumulate fields
    !
    sst_io(:,:) = sst_io(:,:) + tn(:,:,1) + rt0
    sss_io(:,:) = sss_io(:,:) + sn(:,:,1)


    ! vectors at F-point
    DO jj = 2, jpj
       DO ji = fs_2, jpi   ! vector opt.
          u_io(ji,jj) = u_io(ji,jj) + 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
          v_io(ji,jj) = v_io(ji,jj) + 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
       END DO
    END DO

!!DB: 2009.09.30
!    IF( MOD( kt-1, nfice ) == 0 ) THEN
    IF( MOD( kt-nit000, nfice ) == 0 ) THEN

       ! The LIM model is going to be call
       sst_io(:,:) = sst_io(:,:) / FLOAT( nfice ) * tmask(:,:,1)
       sss_io(:,:) = sss_io(:,:) / FLOAT( nfice )

       ! stress from ocean U- and V-points to ice U,V point
       DO jj = 2, jpj
          DO ji = fs_2, jpi   ! vector opt.
             gtaux(ji,jj) = 0.5 * ( taux(ji-1,jj  ) + taux(ji-1,jj-1) )
             gtauy(ji,jj) = 0.5 * ( tauy(ji  ,jj-1) + tauy(ji-1,jj-1) )
             u_io  (ji,jj) = u_io(ji,jj) / FLOAT( nfice )
             v_io  (ji,jj) = v_io(ji,jj) / FLOAT( nfice )
          END DO
       END DO

       ! lateral boundary condition
       CALL lbc_lnk( gtaux(:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
       CALL lbc_lnk( gtauy(:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
       CALL lbc_lnk( u_io (:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)
       CALL lbc_lnk( v_io (:,:), 'I', -1. )   ! I-point (i.e. ice U-V point)

       !!gmbug  in the ocean freezing point computed as :
       !!gm           fzptn (ji,jj) = ( -0.0575 + 1.710523e-3 * SQRT( sn(ji,jj,1) )   &
       !!gm                                     - 2.154996e-4 *       sn(ji,jj,1)   ) * sn(ji,jj,1)   !!   &
       !!gm           !!                        - 7.53e-4 * pressure
       !!gm
       !!gm!bug this is much more accurate and efficient computation
       !!gm       **************************************************
       !!gm freezing point from unesco:
       !!gm     real function tf(s,p)
       !!gm   function to compute the freezing point of seawater
       !!gm
       !!gm   reference: unesco tech. papers in the marine science no. 28. 1978
       !!gm   eighth report jpots
       !!gm   annex 6 freezing point of seawater f.j. millero pp.29-35.
       !!gm
       !!gm  units:
       !!gm         pressure      p          decibars
       !!gm         salinity      s          pss-78
       !!gm         temperature   tf         degrees celsius
       !!gm         freezing pt.
       !!gm************************************************************
       !!gm  checkvalue: tf= -2.588567 deg. c for s=40.0, p=500. decibars
       !!gm     tf=(-.0575+1.710523e-3*sqrt(abs(s))-2.154996e-4*s)*s-7.53e-4*p
       !!gm     return
       !!gm     end
       !!gm!bug


       !!gm      DO jj = 1, jpj
       !!gm         DO ji = 1, jpi
       !!gm            tfu(ji,jj)  = (  rt0 + ( - 0.0575                              &
       !!gm               &                     + 1.710523e-3 * SQRT( sss_io(ji,jj) )   &
       !!gm               &                     - 2.154996e-4 *       sss_io(ji,jj)   ) * sss_io(ji,jj)  ) * tms(ji,jj)
       !!gm         END DO
       !!gm      END DO
       !!gm
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
       numit = numit + nfice 

       !                                                           !--------------!
       CALL lim_dyn                                                ! Ice dynamics !   ( rheology/dynamics )
       !                                                           !--------------!
       IF(ln_ctl) THEN
          CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  2 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
          CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   2 : ', tab2d_2=sist  , clinfo2=' sist      : ')
       ENDIF


       !                                                           !---------------!
       CALL lim_trp                                                ! Ice transport !  ( Advection/diffusion )
       !                                                           !---------------!
       IF(ln_ctl) THEN
          CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  3 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
          CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   3 : ', tab2d_2=sist  , clinfo2=' sist      : ')
       ENDIF

       !                                                           !-------------!
       IF( ln_limdmp ) CALL lim_dmp(kt)                            ! Ice damping !
       !                                                           !-------------!

       !                                                           !--------------------!
       CALL lim_thd                                                ! Ice thermodynamics !
       !                                                           !--------------------!
       IF(ln_ctl) THEN
          CALL prt_ctl(tab2d_1=hsnif ,clinfo1=' hsnif  4 : ', tab2d_2=hicif , clinfo2=' hicif     : ')
          CALL prt_ctl(tab2d_1=frld  ,clinfo1=' frld   4 : ', tab2d_2=sist  , clinfo2=' sist      : ')
       ENDIF


       ! Mass and heat fluxes from ice to ocean
       !                                                           !------------------------------!
       CALL lim_flx                                                ! Ice/Ocean Mass & Heat fluxes !
       !                                                           !------------------------------!

       IF( MOD( numit, ninfo ) == 0 .OR. ntmoy == 1 )  THEN        !-----------------!
          CALL lim_dia                                             ! Ice Diagnostics !
       ENDIF                                                       !-----------------!

       !                                                           !-------------!
       CALL lim_wri                                                ! Ice outputs !
       !                                                           !-------------!
!!DB: 2009.09.30: Moved outside of if-block
!       IF( MOD( numit, nstock ) == 0 .OR. numit == nlast ) THEN
!          !                                                        !------------------!
!          !!DB: old IOIPSL call
!          !            CALL lim_rst_write( numit )                              ! Ice restart file !
!          !                                                        !------------------!
!          !!DB
!          call rst_ice_write(numit)
!       ENDIF


#if defined key_passivetrc || defined key_BGCM_01
       !!DB 2009.08.17 -- need sw_rad variable for BGCM. For most situations this is 
       !! the variable qsr_oce, which is zeroed below ====> save a copy of it
!!DBG
!       sw_rad(:,:) = qsr_oce(:,:)
#endif

       ! Re-initialization of forcings
!DL       qsr_oce (:,:) = 0.e0
       qsr_ice (:,:) = 0.e0
       qnsr_oce(:,:) = 0.e0
       qnsr_ice(:,:) = 0.e0 
       dqns_ice(:,:) = 0.e0 
       tprecip (:,:) = 0.e0 
       sprecip (:,:) = 0.e0
       fr1_i0  (:,:) = 0.e0
       fr2_i0  (:,:) = 0.e0
       evap    (:,:) = 0.e0
#if defined key_coupled 
       rrunoff (:,:) = 0.e0
       calving (:,:) = 0.e0
#else
       qla_ice (:,:) = 0.e0
       dqla_ice(:,:) = 0.e0
#endif

    ENDIF   !!  IF( MOD( kt-1, nfice ) == 0 ) 

!!DB: 2009.09.30: Force ice restarts at same kt as restart()
    if( mod( kt, nstock ) == 0 .OR. kt == nitend ) THEN
       call rst_ice_write(kt)
    endif



  END SUBROUTINE ice_stp

#else
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module          NO LIM sea-ice model
   !!----------------------------------------------------------------------
   USE in_out_manager
CONTAINS
   SUBROUTINE ice_stp ( kt )     ! Dummy routine
      if(lpw)WRITE(numout,*) 'ice_stp: You should not have seen this print! error?', kt
   END SUBROUTINE ice_stp
#endif

   !!======================================================================
END MODULE icestp
