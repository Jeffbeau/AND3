MODULE opa
   !!==============================================================================
   !!                       ***  MODULE opa   ***
   !! Ocean system   : OPA ocean dynamics (including on-line tracers and sea-ice)
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   opa_model      : solve ocean dynamics, tracer and/or sea-ice
   !!   opa_flg        : initialisation of algorithm flag 
   !!   opa_closefile  : close remaining files
   !!----------------------------------------------------------------------
   !! * Modules used
   USE cpl_oce         ! ocean-atmosphere-sea ice coupled exchanges
   USE dom_oce         ! ocean space domain variables
   USE oce             ! dynamics and tracers variables
   USE trdmod_oce      ! ocean variables trends
   USE daymod          ! calendar
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing

   USE domcfg          ! domain configuration               (dom_cfg routine)
   USE mppini          ! shared/distributed memory setting (mpp_init routine)
   USE domain          ! domain initialization             (dom_init routine)
   USE obc_par         ! open boundary cond. parameters
   USE obcini          ! open boundary cond. initialization (obc_ini routine)
   USE solver          ! solver initialization          (solver_init routine)
   USE istate          ! initial state setting          (istate_init routine)
   USE eosbn2          ! equation of state            (eos bn2 routine)
   USE zpshde          ! partial step: hor. derivative (zps_hde routine)

   ! ocean physics
   USE traqsr          ! solar radiation penetration   (tra_qsr_init routine)
   USE ldfdyn          ! lateral viscosity setting      (ldfdyn_init routine)
   USE ldftra          ! lateral diffusivity setting    (ldftra_init routine)
   USE zdfini

   USE phycst          ! physical constant                  (par_cst routine)
   USE iceini          ! initialization of sea-ice         (ice_init routine)
   USE cpl             ! coupled ocean/atmos.              (cpl_init routine)
   USE ocfzpt          ! ocean freezing point              (oc_fz_pt routine)
   USE trdicp          ! momentum/tracers trends       (trd_icp_init routine)
   USE trdvor          ! vorticity trends              (trd_vor_init routine)
   USE trdmld          ! tracer mixed layer trends     (trd_mld_init routine)
   USE flxfwb          ! freshwater budget correction  (flx_fwb_init routine)
   USE flxmod          ! thermohaline forcing of the ocean (flx_init routine)

   USE diaptr          ! poleward transports           (dia_ptr_init routine)

   USE step            ! OPA time-stepping                  (stp     routine)
   USE dynspg_oce      ! Control choice of surface pressure gradient schemes
   USE prtctl          ! Print control                 (prt_ctl_init routine)
   USE ini1d           ! re-initialization of u-v mask for the 1D configuration
   USE dyncor1d        ! Coriolis factor at T-point
   USE step1d          ! Time stepping loop for the 1D configuration

   USE initrc          ! Initialization of the passive tracers

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC opa_model      ! called by model.F90
   PUBLIC opa_init
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/opa.F90,v 1.23 2006/04/19 14:43:14 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE opa_model
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE opa  ***
      !!
      !! ** Purpose :   opa solves the primitive equations on an orthogonal 
      !!      curvilinear mesh on the sphere.
      !!
      !! ** Method  : - model general initialization
      !!              - launch the time-stepping (stp routine)
      !!
      !! References :
      !!      Madec, Delecluse,Imbard, and Levy, 1997: reference manual.
      !!              internal report, IPSL.
      !!
      !! History :
      !!   4.0  !  90-10  (C. Levy, G. Madec)  Original code
      !!   7.0  !  91-11  (M. Imbard, C. Levy, G. Madec)
      !!   7.1  !  93-03  (M. Imbard, C. Levy, G. Madec, O. Marti,
      !!                   M. Guyon, A. Lazar, P. Delecluse, C. Perigaud,
      !!                   G. Caniaux, B. Colot, C. Maes ) release 7.1 
      !!        !  92-06  (L.Terray) coupling implementation
      !!        !  93-11  (M.A. Filiberti) IGLOO sea-ice 
      !!   8.0  !  96-03  (M. Imbard, C. Levy, G. Madec, O. Marti,
      !!                   M. Guyon, A. Lazar, P. Delecluse, L.Terray,
      !!                   M.A. Filiberti, J. Vialar, A.M. Treguier,
      !!                   M. Levy)  release 8.0
      !!   8.1  !  97-06  (M. Imbard, G. Madec)
      !!   8.2  !  99-11  (M. Imbard, H. Goosse)  LIM sea-ice model 
      !!        !  99-12  (V. Thierry, A-M. Treguier, M. Imbard, M-A. Foujols)  OPEN-MP 
      !!        !  00-07  (J-M Molines, M. Imbard)  Open Boundary Conditions  (CLIPPER)
      !!   9.0  !  02-08  (G. Madec)  F90: Free form and modules
      !!    "   !  04-08  (C. Talandier) New trends organization
      !!    "   !  05-06  (C. Ethe) Add the 1D configuration possibility
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   istp       ! time step index
      CHARACTER (len=64) ::        &
         cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing
      !!----------------------------------------------------------------------

#if defined key_agrif

      Call Agrif_Init_Grids()
#endif
      
      Call opa_init  ! Initializations

      IF( lk_cfg_1d  )  THEN 
         istp = nit000
         DO WHILE ( istp <= nitend .AND. nstop == 0 )
#if defined key_agrif
            CALL Agrif_Step(stp_1d)
#else
            CALL stp_1d( istp )
#endif
            istp = istp + 1
         END DO
      ELSE
         istp = nit000
         DO WHILE ( istp <= nitend .AND. nstop == 0 )
#if defined key_agrif
            CALL Agrif_Step(stp)
#else
            CALL stp( istp )
#endif
            istp = istp + 1
         END DO
      ENDIF
      !                                     ! ========= !
      !                                     !  Job end  !
      !                                     ! ========= !

      IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA

      IF( nstop /= 0 ) THEN                 ! error print
      IF(lwp) WRITE(numout,cform_err)
      IF(lwp) WRITE(numout,*) nstop, ' error have been found' 
      ENDIF

      CALL opa_closefile
      IF( lk_mpp )   CALL mppstop                          ! Close all files (mpp)

   END SUBROUTINE opa_model


   SUBROUTINE opa_flg
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE opa  ***
      !!
      !! ** Purpose :   Initialize logical flags that control the choice of
      !!      some algorithm or control print
      !!
      !! ** Method  :    Read in namilist namflg logical flags
      !!
      !! History :
      !!   9.0  !  03-11  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations

      NAMELIST/namflg/ ln_dynhpg_imp
      !!----------------------------------------------------------------------

      ! Read Namelist namflg : algorithm FLaG 
      ! --------------------
      REWIND ( numnam )
      READ   ( numnam, namflg )

      ! Parameter control and print
      ! ---------------------------
      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'opa_flg : algorithm flag initialization'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '          Namelist namflg : set algorithm flags'
         WRITE(numout,*)
         WRITE(numout,*) '             centered (F) or semi-implicit (T)   ln_dynhpg_imp = ', ln_dynhpg_imp
         WRITE(numout,*) '             hydrostatic pressure gradient'
      ENDIF

   END SUBROUTINE opa_flg

   SUBROUTINE opa_closefile
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE opa_closefile  ***
      !!
      !! ** Purpose :   Close the files
      !!           
      !! ** Method  : 
      !!
      !! History :
      !!   9.0  !  05-01  (O. Le Galloudec)  Original code
      !!----------------------------------------------------------------------
      !! * Modules used
      USE dtatem        ! temperature data
      USE dtasal        ! salinity data
      USE dtasst        ! sea surface temperature data
      !!----------------------------------------------------------------------

      IF ( lk_mpp ) CALL mppsync

      ! 1. Unit close
      ! -------------

      CLOSE( numnam )       ! namelist
      CLOSE( numout )       ! standard model output file
      CLOSE( numstp )       ! time-step file
      CLOSE( numwrs )       ! ocean restart file

      IF( lk_dtatem )   CLOSE( numtdt )
      IF( lk_dtasal )   CLOSE( numsdt )
      IF( lk_dtasst )   CLOSE( numsst )

      IF(lwp) CLOSE( numsol )

      IF( lk_cpl ) THEN
         CLOSE( numlhf )
         CLOSE( numlts )
      ENDIF

      CLOSE( numwri )

   END SUBROUTINE opa_closefile

   !!======================================================================
   SUBROUTINE opa_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE opa_init  ***
      !!
      !! ** Purpose :   initialization of the opa model
      !!
      !! ** Method  : 
      !!
      !! References :
      !!----------------------------------------------------------------------
      !! * Local declarations

#if defined key_coupled
      INTEGER ::   itro, istp0        ! ???
#endif
      CHARACTER (len=64) ::        &
         cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing
      CHARACTER (len=20) :: namelistname
      CHARACTER (len=28) :: file_out
      !!----------------------------------------------------------------------

      ! Initializations
      ! ===============

      file_out = 'ocean.output'
      
      ! open listing and namelist units
      IF ( numout /= 0 .AND. numout /= 6 ) THEN 
         CALL ctlopn(numout,file_out,'UNKNOWN', 'FORMATTED',   &
                      'SEQUENTIAL',1,numout,.FALSE.,1)
!         OPEN( UNIT=numout, FILE=TRIM(file_out), FORM='FORMATTED' )
      ENDIF

      namelistname = 'namelist'
      CALL ctlopn(numnam,namelistname,'OLD', 'FORMATTED', 'SEQUENTIAL',   &
                     1,numout,.FALSE.,1)
!!!!      OPEN( UNIT=numnam, FILE='namelist', FORM='FORMATTED', STATUS='OLD' )

!!DB 2009.09.09 -- restrict write to lwp
!      WRITE(numout,*)
!      WRITE(numout,*) '                 L O D Y C - I P S L'
!      WRITE(numout,*) '                     O P A model'
!      WRITE(numout,*) '            Ocean General Circulation Model'
!      WRITE(numout,*) '               version OPA 9.0  (2005) '
!      WRITE(numout,*)
!      WRITE(numout,*)

      ! Nodes selection
      narea = mynode()
      narea = narea + 1    ! mynode return the rank of proc (0 --> jpnij -1 )
      lwp   = narea == 1

      !                                     ! ============================== !
      !                                     !  Model general initialization  !
      !                                     ! ============================== !

      IF(lwp) then
         WRITE(numout,*)
         WRITE(numout,*) '                 L O D Y C - I P S L'
         WRITE(numout,*) '                     O P A model'
         WRITE(numout,*) '            Ocean General Circulation Model'
         WRITE(numout,*) '               version OPA 9.0  (2005) '
         WRITE(numout,*)
         WRITE(numout,*)
         WRITE(numout,cform_aaa)       ! Flag AAAAAAA
      endif
      
!!DB 2009.10.01 open run_info.out file as numout2 (in_out_manager)
      if(lwp) then 
         open(numout2, file = 'sopa_run_info.out',STATUS='replace')
         write(numout2,*) ' SOPA MC runtime information file '
         write(numout2,*) ' -------------------------------- '
         write(numout2,*) '                                  '
         write(numout2,*) '                 L O D Y C - I P S L'
         write(numout2,*) '                     O P A model'
         write(numout2,*) '            Ocean General Circulation Model'
         write(numout2,*) '               version OPA 9.0  (2005) '
         write(numout2,*) '                                  '
         write(numout2,*) '                                  '
      endif

                                            ! Domain decomposition
      IF( jpni*jpnj == jpnij ) THEN
         CALL mpp_init                          ! standard cutting out
      ELSE
         CALL mpp_init2                         ! eliminate land processors
      ENDIF
      
      CALL phy_cst                          ! Physical constants

      CALL dom_cfg                          ! Domain configuration
      
      CALL dom_init                         ! Domain
      IF( ln_ctl )      CALL prt_ctl_init   ! Print control

      IF( lk_cfg_1d )   CALL fcorio_1d      ! redefine Coriolis at T-point

      IF( lk_obc    )   CALL obc_init       ! Open boundaries 

      CALL day( nit000 )                    ! Calendar

      CALL istate_init                      ! ocean initial state (Dynamics and tracers)

      IF( lk_dynspg_flt .OR. lk_dynspg_rl ) THEN
         CALL solver_init( nit000 )         ! Elliptic solver
      ENDIF

!!add
                       CALL eos( tb, sb, rhd, rhop )        ! before potential and in situ densities
                       
                       CALL bn2( tb, sb, rn2 )              ! before Brunt-Vaisala frequency

      IF( lk_zps .AND. .NOT. lk_cfg_1d )   &
         &             CALL zps_hde( nit000, tb, sb, rhd,  &  ! Partial steps: before Horizontal DErivative
                                            gtu, gsu, gru, &  ! of t, s, rd at the bottom ocean level
                                            gtv, gsv, grv )

!!add

      CALL oc_fz_pt                         ! Surface freezing point

#if defined key_ice_lim
      CALL ice_init                         ! Sea ice model
#endif

      !                                     ! Ocean scheme

      CALL opa_flg                              ! Choice of algorithms

      !                                     ! Ocean physics

      CALL tra_qsr_init                         ! Solar radiation penetration

      CALL ldf_dyn_init                         ! Lateral ocean momentum physics

      CALL ldf_tra_init                         ! Lateral ocean tracer physics

      CALL zdf_init                             ! Vertical ocean physics

      !                                     ! Ocean trends
      ! Control parameters 
      IF( lk_trdtra .OR. lk_trdmld )   l_trdtra = .TRUE.
      IF( lk_trddyn .OR. lk_trdvor )   l_trddyn = .TRUE.

      IF( lk_trddyn .OR. lk_trdtra )   &
         &            CALL trd_icp_init         ! active tracers and/or momentum

      IF( lk_trdmld ) CALL trd_mld_init         ! mixed layer

      IF( lk_trdvor ) CALL trd_vor_init         ! vorticity

#if defined key_passivetrc
      CALL ini_trc                           ! Passive tracers
#endif

#if defined key_coupled
      itro  = nitend - nit000 + 1           ! Coupled
      istp0 = NINT( rdt )
      CALL cpl_init( itro, nexco, istp0 )   ! Signal processing and process id exchange
#endif

      CALL flx_init                         ! Thermohaline forcing initialization

      CALL flx_fwb_init                     ! FreshWater Budget correction

      CALL dia_ptr_init                     ! Poleward TRansports initialization

      !                                     ! =============== !
      !                                     !  time stepping  !
      !                                     ! =============== !


!!DB: 2007.12.06
!!Code for additional control parameters
!!Look for a file called run_params and get run params
!!This is a HARDWIRED list. At this time, the default values
!!should be found at the bottom of oce.F90
!!Exit gracefully if file not found
      if(lwp)write(numout2,*)'DBG Default perpetual_forcing,M2_ave,ioutput_ave: ', &
           perpetual_forcing,M2_ave,ioutput_ave
      open(66,file='run_params',STATUS='OLD',ERR=666)
      read(66,*,ERR=666)perpetual_forcing
      read(66,*,ERR=666)M2_ave
      read(66,*,ERR=666)ioutput_ave
      close(66)
!!DBG
      if(lwp) then
         write(numout2,*)'DBG: run_params --> perpetual_forcing,M2_ave,ioutput_ave: ', &
              perpetual_forcing,M2_ave, ioutput_ave
      endif
666   continue

!!DB 2009.10.01
      if(ioutput_ave /= 0) then
         ioutput_ave = nwrite
         if(lwp) write(numout2,*)'DBG: ioutput_ave /= 0 ; setting to nwrite  = ', nwrite
      endif


      IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA

      IF( lk_cfg_1d  )  THEN 
         CALL init_1d
      ENDIF
   END SUBROUTINE opa_init
   !!======================================================================
END MODULE opa
