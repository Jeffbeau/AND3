MODULE opa
   !!==============================================================================
   !!                       ***  MODULE opa   ***
   !! Ocean system   : OPA ocean dynamics (including on-line tracers and sea-ice)
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   opa_model      : solve ocean dynamics, tracer and/or sea-ice
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space domain variables
   USE oce             ! dynamics and tracers variables
   USE daymod          ! calendar
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing

   USE domcfg          ! domain configuration               (dom_cfg routine)
   USE mppini          ! shared/distributed memory setting (mpp_init routine)
   USE domain          ! domain initialization             (dom_init routine)
   USE istate          ! initial state setting          (istate_init routine)
   USE eosbn2          ! equation of state            (eos bn2 routine)
   USE zpshde          ! partial step: hor. derivative (zps_hde routine)

   ! ocean physics
   USE ldftra          ! lateral diffusivity setting    (ldftra_init routine)
   USE zdfini
   USE traqsr          ! solar radiation penetration   (tra_qsr_init routine)

   USE phycst          ! physical constant                  (par_cst routine)
   USE dtadyn          ! Lecture and Interpolation of the dynamical fields
   USE initrc          ! Initilization of the passive tracers
   USE step            ! OPA time-stepping                  (stp     routine)

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC opa_model      ! called by model.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/opa.F90,v 1.2 2005/11/16 16:19:34 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
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
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   istp       ! time step index
      CHARACTER (len=64) ::        &
         cform_aaa="( /, 'AAAAAAAA', / ) "     ! flag for output listing
      !!----------------------------------------------------------------------
      
      
      ! Initializations
      ! ===============
      
      ! open listing and namelist units
      IF ( numout /= 0 .AND. numout /= 6 ) THEN 
         OPEN( UNIT=numout, FILE='ocean.output', FORM='FORMATTED' )
      ENDIF

      ! Nodes selection
      narea = mynode()
      narea = narea + 1    ! mynode return the rank of proc (0 --> jpnij -1 )
      lwp   = narea == 1

      OPEN( UNIT=numnam, FILE='namelist', FORM='FORMATTED', STATUS='OLD' )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '                 L O D Y C - I P S L'
         WRITE(numout,*) '                     O P A model'
         WRITE(numout,*) '            Ocean General Circulation Model'
         WRITE(numout,*) '               version OPA 9.0  (2003)'
         WRITE(numout,*)
      ENDIF

      !                                     ! ============================== !
      !                                     !  Model general initialization  !
      !                                     ! ============================== !

      IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA

                                            ! Domain decomposition
      IF( jpni*jpnj == jpnij ) THEN
         CALL mpp_init                          ! standard cutting out
      ELSE
         CALL mpp_init2                         ! eliminate land processors
      ENDIF
      
      CALL phy_cst                          ! Physical constants

      CALL dom_cfg                          ! Domain configuration
      
      CALL dom_init                         ! Domain

      CALL day( nit000 )                    ! Calendar

      CALL istate_init                      ! ocean initial state (Dynamics and tracers)
!!add
                       CALL eos( tn, sn, rhd, rhop )        ! before potential and in situ densities

                       CALL bn2( tn, sn, rn2 )              ! before Brunt-Vaisala frequency

      IF( lk_zps    )   CALL zps_hde( nit000, tn, sn, rhd,  &  ! Partial steps: before Horizontal DErivative
                                          gtu, gsu, gru, &  ! of t, s, rd at the bottom ocean level
                                          gtv, gsv, grv )

!!add

      ! Initialization for the dynamics
      !
      CALL dta_dyn(nit000)     

#if defined key_passivetrc
      CALL ini_trc                           ! Passive tracers
#endif

      !                                     ! Ocean physics
      CALL tra_qsr_init                         ! Solar radiation penetration

      CALL ldf_tra_init                         ! Lateral ocean tracer physics

      CALL zdf_init                             ! Vertical ocean physics

      !                                     ! Ocean trends

      !                                     ! =============== !
      !                                     !  time stepping  !
      !                                     ! =============== !

      IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA

      istp = nit000
      DO WHILE ( istp <= nitend .AND. nstop == 0 )
         CALL stp( istp )
         istp = istp + 1
      END DO
      !                                     ! ========= !
      !                                     !  Job end  !
      !                                     ! ========= !

      IF(lwp) WRITE(numout,cform_aaa)       ! Flag AAAAAAA

      IF( nstop /= 0 ) THEN                 ! error print
      IF(lwp) WRITE(numout,cform_err)
      IF(lwp) WRITE(numout,*) nstop, ' error have been found' 
      ENDIF

      IF( lk_mpp )   CALL mppstop                          ! Close all files (mpp)

   END SUBROUTINE opa_model

   !!======================================================================
END MODULE opa
