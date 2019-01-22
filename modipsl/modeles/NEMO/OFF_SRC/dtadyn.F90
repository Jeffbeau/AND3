MODULE dtadyn
   !!======================================================================
   !!                       ***  MODULE  dtadyn  ***
   !! OFFLINE : interpolation of the physical fields
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   dta_dyn_init : initialization, namelist read, and parameters control
   !!   dta_dyn      : Interpolation of the fields
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE ocesbc
   USE ldfslp
   USE blk_oce
   USE ldfeiv          ! eddy induced velocity coef.      (ldf_eiv routine)
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE zdfmxl
   USE trabbl          ! tracers: bottom boundary layer
   USE ocfzpt
   USE zdfddm          ! vertical  physics: double diffusion
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing library

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC dta_dyn_init   ! called by opa.F90
   PUBLIC dta_dyn        ! called by step.F90

   !! * Module variables
   INTEGER , PUBLIC, PARAMETER :: jpflx = 13
   INTEGER , PUBLIC, PARAMETER :: &
      jptaux = 1 , & ! indice in flux for taux
      jptauy = 2 , & ! indice in flux for tauy
      jpwind = 3 , & ! indice in flux for wspd
      jpemp = 4  , & ! indice in flux for E-P
      jpice = 5  , & ! indice in flux for ice concentration
      jpqsr = 6      ! indice in flux for shortwave heat flux 

   LOGICAL , PUBLIC :: &
      lperdyn = .TRUE. , & ! boolean for periodic fields or not
      lfirdyn = .TRUE.     ! boolean for the first call or not

   INTEGER , PUBLIC :: &
      ndtadyn = 12 ,  & ! Number of dat in one year
      ndtatot = 12 ,  & ! Number of data in the input field
      nsptint = 1 ,   & ! type of spatial interpolation
      nficdyn = 2       ! number of dynamical fields 

   INTEGER :: ndyn1, ndyn2 , &
      nlecoff = 0  , & ! switch for the first read
      numfl_t, numfl_u, &
      numfl_v, numfl_w, numfl_s
      

   REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
      tdta   ,   & ! temperature at two consecutive times
      sdta   ,   & ! salinity at two consecutive times
      udta   ,   & ! zonal velocity at two consecutive times
      vdta   ,   & ! meridional velocity at two consecutive times
      wdta   ,   & ! vertical velocity at two consecutive times
      avtdta       ! vertical diffusivity coefficient

#if defined key_ldfslp
   REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
      uslpdta ,  & ! zonal isopycnal slopes
      vslpdta ,  & ! meridional isopycnal slopes
      wslpidta , & ! zonal diapycnal slopes
      wslpjdta     ! meridional diapycnal slopes
#endif

#if defined key_traldf_eiv   &&   defined key_traldf_c2d
   REAL(wp), DIMENSION(jpi,jpj,2) ::   &
      ahtwdta ,  & ! Lateral diffusivity
      eivwdta      ! G&M coefficient
#endif

   REAL(wp), DIMENSION(jpi,jpj,jpflx,2) ::   &
      flxdta       ! auxiliary 2-D forcing fields at two consecutive times
   REAL(wp), DIMENSION(jpi,jpj,2) ::       &
      zmxldta      ! mixed layer depth at two consecutive times

#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
   REAL(wp), DIMENSION(jpi,jpj,2) ::       &
      bblxdta ,  & ! frequency of bbl in the x direction at 2 consecutive times
      bblydta      ! frequency of bbl in the y direction at 2 consecutive times
#endif

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/dtadyn.F90,v 1.5 2006/04/26 09:32:54 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dta_dyn_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dta_dyn_init  ***
      !!
      !! ** Purpose :   initializations of parameters for the interpolation
      !!
      !! ** Method :
      !!
      !! History :
      !!    ! original  : 92-01 (M. Imbard: sub domain)
      !!    ! 98-04 (L.Bopp MA Foujols: slopes for isopyc.)
      !!    ! 98-05 (L. Bopp read output of coupled run)
      !!    ! 05-03 (O. Aumont and A. El Moussaoui) F90
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Local declarations


      NAMELIST/nam_offdyn/ ndtadyn, ndtatot, nsptint,            & 
          &                nficdyn, lperdyn
      !!----------------------------------------------------------------------

      !  Define the dynamical input parameters
      ! ======================================

      ! Read Namelist nam_offdyn : Lateral physics on tracers
      REWIND( numnam )
      READ  ( numnam, nam_offdyn )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'nam_offdyn : offline dynamical selection'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '  Namelist nam_offdyn : set parameters for the lecture of the dynamical fields'
         WRITE(numout,*) 
         WRITE(numout,*) ' number of elements in the FILE for a year  ndtadyn = ' , ndtadyn
         WRITE(numout,*) ' total number of elements in the FILE       ndtatot = ' , ndtatot
         WRITE(numout,*) ' type of interpolation                      nsptint = ' , nsptint
         WRITE(numout,*) ' number of dynamics FILE                    nficdyn = ' , nficdyn
         WRITE(numout,*) ' loop on the same FILE                      lperdyn = ' , lperdyn
         WRITE(numout,*) ' '
      ENDIF

   END SUBROUTINE dta_dyn_init

   SUBROUTINE dta_dyn(kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dta_dyn  ***
      !!
      !! ** Purpose : Prepares dynamics and physics fields from an 
      !!              OPA9 simulation  for an off-line simulation
      !!               for passive tracer
      !!
      !! ** Method : calculates the position of DATA to read READ DATA 
      !!             (example month changement) computes slopes IF needed
      !!             interpolates DATA IF needed
      !!
      !! ** History :
      !!   ! original  : 92-01 (M. Imbard: sub domain)
      !!   ! addition  : 98-04 (L.Bopp MA Foujols: slopes for isopyc.) 
      !!   ! addition  : 98-05 (L. Bopp read output of coupled run)
      !!   ! addition  : 05-03 (O. Aumont and A. El Moussaoui) F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE eosbn2

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index

      !! * Local declarations
      INTEGER ::   iper, iperm1, iswap   

      REAL(wp) :: zpdtan, zpdtpe, zdemi, zt
      REAL(wp) :: zweigh, zweighm1

      REAL(wp), DIMENSION(jpi,jpj,jpflx) ::   &
         flx  ! auxiliary field for 2-D surface boundary conditions


      ! 0. Initialization
      ! -----------------

      IF (lfirdyn) THEN
      !
      ! time step MUST BE nint000
      !
          IF (kt.ne.nit000) THEN
              IF (lwp) THEN
                  WRITE (numout,*) ' kt MUST BE EQUAL to nit000. kt=',kt  &
                     ,' nit000=',nit000
              END IF
              STOP 'dtadyn'
          END if
      ! Initialize the parameters of the interpolation
      CALL dta_dyn_init
      ENDIF


      zpdtan = raass / rdt
      zpdtpe = ((zpdtan / FLOAT (ndtadyn)))
      zdemi  = zpdtpe * 0.5
      zt     = (FLOAT (kt) + zdemi ) / zpdtpe

      zweigh   = zt - FLOAT(INT(zt))
      zweighm1 = 1. - zweigh

      IF (lperdyn) THEN
         iperm1 = MOD(INT(zt),ndtadyn)
      ELSE
         iperm1 = MOD(INT(zt),(ndtatot-1))
      ENDIF
      iper = iperm1 + 1
      IF (iperm1 == 0) THEN
          IF (lperdyn) THEN
              iperm1 = ndtadyn
          ELSE 
              IF (lfirdyn) THEN
                  IF (lwp) THEN 
                      WRITE (numout,*) ' dynamic file is not periodic '
                      WRITE (numout,*) ' with or without interpolation, '
                      WRITE (numout,*) ' we take the first value'
                      WRITE (numout,*) ' for the previous period '
                      WRITE (numout,*) ' iperm1 = 0  '
                  END IF 
              END IF
          END IF 
      END IF 

      iswap  = 0

      ! 1. First call lfirdyn = true
      ! ----------------------------

      IF (lfirdyn) THEN
      !
      ! store the information of the period read
      !
          ndyn1 = iperm1
          ndyn2 = iper

          IF (lwp) THEN
              WRITE (numout,*)         &
                 ' dynamics DATA READ for the period ndyn1 =',ndyn1, &
              & ' and for the period ndyn2 = ',ndyn2
              WRITE (numout,*) ' time step is :',kt
              WRITE (numout,*) ' we have ndtadyn = ',ndtadyn,&
                 &         ' records in the dynamic FILE for one year'
          END IF 
      !
      ! DATA READ for the iperm1 period
      !
          IF( iperm1 .NE. 0 ) THEN
             CALL dynrea( kt, iperm1 ) 
          ELSE 
             CALL dynrea( kt, 1 )
          ENDIF
      !
      ! Computes dynamical fields
      !
                tn(:,:,:)=tdta(:,:,:,2)
                sn(:,:,:)=sdta(:,:,:,2)
                avt(:,:,:)=avtdta(:,:,:,2)

         IF(lwp) THEN
            WRITE(numout,*)' temperature '
            WRITE(numout,*)
            CALL prihre(tn(1,1,1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)
            WRITE(numout,*) '  level = ',jpk/2
            CALL prihre(tn(1,1,jpk/2),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout)  
            WRITE(numout,*) '  level = ',jpkm1
            CALL prihre(tn(1,1,jpkm1),jpi,jpj,1,jpi,20,1,jpj,20,1.,numout) 
        ENDIF

#if defined key_ldfslp
            CALL eos( tn, sn, rhd, rhop )   ! Time-filtered in situ density 
            CALL bn2( tn, sn, rn2 )         ! before Brunt-Vaisala frequency
            CALL zdf_mxl( kt )              ! mixed layer depth
            CALL ldf_slp( kt, rhd, rn2 )

            uslpdta(:,:,:,2)=uslp(:,:,:)
            vslpdta(:,:,:,2)=vslp(:,:,:)
            wslpidta(:,:,:,2)=wslpi(:,:,:)
            wslpjdta(:,:,:,2)=wslpj(:,:,:)
#endif
       !
       ! swap from record 2 to 1
       !
                udta(:,:,:,1)=udta(:,:,:,2)
                vdta(:,:,:,1)=vdta(:,:,:,2)
                wdta(:,:,:,1)=wdta(:,:,:,2)
                avtdta(:,:,:,1)=avtdta(:,:,:,2)
                tdta(:,:,:,1)=tdta(:,:,:,2)
                sdta(:,:,:,1)=sdta(:,:,:,2)
#if defined key_ldfslp
                uslpdta(:,:,:,1)=uslpdta(:,:,:,2)
                vslpdta(:,:,:,1)=vslpdta(:,:,:,2)
                wslpidta(:,:,:,1)=wslpidta(:,:,:,2)
                wslpjdta(:,:,:,1)=wslpjdta(:,:,:,2)
#endif
                flxdta(:,:,:,1) = flxdta(:,:,:,2)
                zmxldta(:,:,1)=zmxldta(:,:,2)
#if defined key_traldf_eiv   &&   defined key_traldf_c2d
                ahtwdta(:,:,1)=ahtwdta(:,:,2)
                eivwdta(:,:,1)=eivwdta(:,:,2)
#endif
#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
                bblxdta(:,:,1)=bblxdta(:,:,2)
                bblydta(:,:,1)=bblydta(:,:,2)
#endif
      !
      ! indicates a swap
      !
          iswap = 1
      !
      ! DATA READ for the iper period
      !
          CALL dynrea(kt,iper)
      !
      ! Computes wdta (and slopes if key_trahdfiso)
      !
                tn(:,:,:)=tdta(:,:,:,2)
                sn(:,:,:)=sdta(:,:,:,2)
                avt(:,:,:)=avtdta(:,:,:,2)


#if defined key_ldfslp
            CALL eos( tn, sn, rhd, rhop )   ! Time-filtered in situ density
            CALL bn2( tn, sn, rn2 )         ! before Brunt-Vaisala frequency
            CALL zdf_mxl( kt )              ! mixed layer depth
            CALL ldf_slp( kt, rhd, rn2 )

            uslpdta(:,:,:,2)=uslp(:,:,:)
            vslpdta(:,:,:,2)=vslp(:,:,:)
            wslpidta(:,:,:,2)=wslpi(:,:,:)
            wslpjdta(:,:,:,2)=wslpj(:,:,:)
#endif
      !
      ! trace the first CALL
      !
          lfirdyn=.FALSE. 
      ENDIF
      !
      ! and now what we have to DO at every time step
      !
      ! check the validity of the period in memory
      !
      IF (iperm1.NE.ndyn1) THEN 
          IF (iperm1.EQ.0.) THEN
              IF (lwp) THEN
                  WRITE (numout,*) ' dynamic file is not periodic '
                  WRITE (numout,*) ' with or without interpolation, '
                  WRITE (numout,*) ' we take the last value'
                  WRITE (numout,*) ' for the last period '
                  WRITE (numout,*) ' iperm1 = 12  '
                  WRITE (numout,*) ' iper = 13'
              ENDIF
              iperm1 = 12
              iper =13
          ENDIF
      !
      ! we have to prepare a NEW READ of DATA
      !
      ! swap from record 2 to 1
      !
                udta(:,:,:,1)=udta(:,:,:,2)
                vdta(:,:,:,1)=vdta(:,:,:,2)
                wdta(:,:,:,1)=wdta(:,:,:,2)
                avtdta(:,:,:,1)=avtdta(:,:,:,2)
                tdta(:,:,:,1)=tdta(:,:,:,2)
                sdta(:,:,:,1)=sdta(:,:,:,2)
#if defined key_ldfslp
                uslpdta(:,:,:,1)=uslpdta(:,:,:,2)
                vslpdta(:,:,:,1)=vslpdta(:,:,:,2)
                wslpidta(:,:,:,1)=wslpidta(:,:,:,2)
                wslpjdta(:,:,:,1)=wslpjdta(:,:,:,2)
#endif
                flxdta(:,:,:,1) = flxdta(:,:,:,2)
                zmxldta(:,:,1)=zmxldta(:,:,2)
#if defined key_traldf_eiv   &&   defined key_traldf_c2d
                ahtwdta(:,:,1)=ahtwdta(:,:,2)
                eivwdta(:,:,1)=eivwdta(:,:,2)
#endif
#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
                bblxdta(:,:,1)=bblxdta(:,:,2)
                bblydta(:,:,1)=bblydta(:,:,2)
#endif
      !
      ! indicates a swap
      !
          iswap = 1
      !
      ! READ DATA for the iper period
      !
          CALL dynrea(kt,iper)
      !
      ! Computes wdta (and slopes if key_trahdfiso)
      !
                tn(:,:,:)=tdta(:,:,:,2)
                sn(:,:,:)=sdta(:,:,:,2)
                avt(:,:,:)=avtdta(:,:,:,2)

#if defined key_ldfslp
            CALL eos( tn, sn, rhd, rhop )   ! Time-filtered in situ density
            CALL bn2( tn, sn, rn2 )         ! before Brunt-Vaisala frequency
            CALL zdf_mxl( kt )              ! mixed layer depth
            CALL ldf_slp( kt, rhd, rn2 )

            uslpdta(:,:,:,2)=uslp(:,:,:)
            vslpdta(:,:,:,2)=vslp(:,:,:)
            wslpidta(:,:,:,2)=wslpi(:,:,:)
            wslpjdta(:,:,:,2)=wslpj(:,:,:)
#endif
       !
       ! store the information of the period read
       !
          ndyn1 = ndyn2
          ndyn2 = iper
       !
       ! we have READ another period of DATA
       !
          IF (lwp) THEN
              WRITE (numout,*) ' dynamics DATA READ for the period ndyn1 =',ndyn1
              WRITE (numout,*) ' and the period ndyn2 = ',ndyn2
              WRITE (numout,*) ' time step is :',kt
          END IF 

      END IF 

      !
      ! compute the DATA at the given time step
      !
      IF (nsptint.eq.0) THEN
      !
      ! no spatial interpolation
      !
      ! DATA are probably correct 
      ! we have to initialize DATA IF we have changed the period
      !
          IF (iswap.eq.1) THEN
      !
      ! initialize now fields with the NEW DATA READ
      !
                    un(:,:,:)=udta(:,:,:,2)
                    vn(:,:,:)=vdta(:,:,:,2)
                    wn(:,:,:)=wdta(:,:,:,2)
#if defined key_trc_zdfddm
                    avs(:,:,:)=avtdta(:,:,:,2)
#endif
                    avt(:,:,:)=avtdta(:,:,:,2)
                    tn(:,:,:)=tdta(:,:,:,2)
                    sn(:,:,:)=sdta(:,:,:,2)
#if defined key_ldfslp
                    uslp(:,:,:)=uslpdta(:,:,:,2)
                    vslp(:,:,:)=vslpdta(:,:,:,2)
                    wslpi(:,:,:)=wslpidta(:,:,:,2)
                    wslpj(:,:,:)=wslpjdta(:,:,:,2)
#endif
                    flx(:,:,:) = flxdta(:,:,:,2)
                    hmld(:,:)=zmxldta(:,:,2)
#if defined key_traldf_eiv   &&   defined key_traldf_c2d
                    ahtw(:,:)=ahtwdta(:,:,2)
                    aeiw(:,:)=eivwdta(:,:,2)
#endif
#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
                    bblx(:,:)=bblxdta(:,:,2)
                    bbly(:,:)=bblydta(:,:,2)
#endif
       !
       ! keep needed fluxes
       !
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
                    vatm(:,:) = flx(:,:,jpwind)
#endif
                    freeze(:,:) = flx(:,:,jpice)
                    emp(:,:) = flx(:,:,jpemp)
                    emps(:,:) = emp(:,:)
                    qsr(:,:) = flx(:,:,jpqsr)

          END IF 

      ELSE 
          IF (nsptint.eq.1) THEN
      !
      ! linear interpolation
      !
      ! initialize now fields with a linear interpolation
      !
                    un(:,:,:) = zweighm1 * udta(:,:,:,1) + zweigh * udta(:,:,:,2) 
                    vn(:,:,:) = zweighm1 * vdta(:,:,:,1) + zweigh * vdta(:,:,:,2)
                    wn(:,:,:) = zweighm1 * wdta(:,:,:,1) + zweigh * wdta(:,:,:,2)
#if defined key_zdfddm
                    avs(:,:,:)= zweighm1 * avtdta(:,:,:,1) + zweigh * avtdta(:,:,:,2)
#endif
                    avt(:,:,:)= zweighm1 * avtdta(:,:,:,1) + zweigh * avtdta(:,:,:,2)
                    tn(:,:,:) = zweighm1 * tdta(:,:,:,1) + zweigh * tdta(:,:,:,2)
                    sn(:,:,:) = zweighm1 * sdta(:,:,:,1) + zweigh * sdta(:,:,:,2)
   
         
#if defined key_ldfslp
                    uslp(:,:,:) = zweighm1 * uslpdta(:,:,:,1) + zweigh * uslpdta(:,:,:,2) 
                    vslp(:,:,:) = zweighm1 * vslpdta(:,:,:,1) + zweigh * vslpdta(:,:,:,2) 
                    wslpi(:,:,:) = zweighm1 * wslpidta(:,:,:,1) + zweigh * wslpidta(:,:,:,2) 
                    wslpj(:,:,:) = zweighm1 * wslpjdta(:,:,:,1) + zweigh * wslpjdta(:,:,:,2) 
#endif
                    flx(:,:,:) = zweighm1 * flxdta(:,:,:,1) + zweigh * flxdta(:,:,:,2) 
                    hmld(:,:) = zweighm1 * zmxldta(:,:,1) + zweigh  * zmxldta(:,:,2) 
#if defined key_traldf_eiv   &&   defined key_traldf_c2d
                    ahtw(:,:) =  zweighm1 * ahtwdta(:,:,1) + zweigh * ahtwdta(:,:,2)
                    aeiw(:,:) =  zweighm1 * eivwdta(:,:,1) + zweigh * eivwdta(:,:,2)
#endif
#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
                    bblx(:,:)= zweighm1 * bblxdta(:,:,1) + zweigh * bblxdta(:,:,2)
                    bbly(:,:)= zweighm1 * bblydta(:,:,1) + zweigh * bblydta(:,:,2)
#endif
       !
       ! keep needed fluxes
       !
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
                  vatm(:,:) = flx(:,:,jpwind)
#endif
                  freeze(:,:) = flx(:,:,jpice)
                  emp(:,:) = flx(:,:,jpemp)
                  emps(:,:) = emp(:,:)
                  qsr(:,:) = flx(:,:,jpqsr)
       !
       ! other interpolation
       !
          ELSE 

              WRITE (numout,*) ' this kind of interpolation don t EXIST'
              WRITE (numout,*) ' at the moment. we STOP '
              STOP 'dtadyn'

          END IF 

      END IF
      !
      ! lb in any case, we need rhop
      !
      CALL eos( tn, sn, rhd, rhop ) 

#if defined key_traldf_c2d
      ! In case of 2D varying coefficients, we need aeiv and aeiu
      IF( lk_traldf_eiv )   CALL ldf_eiv( kt )      ! eddy induced velocity coefficient
#endif

   END SUBROUTINE dta_dyn

   SUBROUTINE dynrea( kt, kenr )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dynrea  ***
      !!
      !! ** Purpose : READ dynamics fiels from OPA9 netcdf output
      !! 
      !! ** Method : READ the kenr records of DATA and store in
      !!             in udta(...,2), ....  
      !! 
      !! ** History : additions : M. Levy et M. Benjelloul jan 2001 
      !!              (netcdf FORMAT) 
      !!              05-03 (O. Aumont and A. El Moussaoui) F90
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt, kenr       ! time index
      !! * Local declarations
      INTEGER ::   ji, jj
      INTEGER ::   ipi,ipj,ipk,itime,jkenr,idtatot
      INTEGER , DIMENSION(ndtatot) :: istep

      REAL(wp) ::  zdate0

      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
        zu, zv, zw, zt, zs, zavt ! 3-D dynamical fields

# if defined key_traldf_eiv
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
        zaeiu, zaeiv, zaeiw
# endif

# if defined key_traldf_eiv   &&   defined key_traldf_c2d
      REAL(wp), DIMENSION(jpi,jpj) ::   &
        zeivw, zahtw
# endif

      REAL(wp), DIMENSION(jpi,jpj) :: &
        zlon, zlat, zemp, zqsr, zmld, zice, zwind 
#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
      REAL(wp), DIMENSION(jpi,jpj) :: &
        zbblx, zbbly
#endif
      REAL(wp), DIMENSION(jpk) :: zlev

      CHARACTER(len=45)  ::  &
         clname_t = 'dyna_grid_T.nc', &
         clname_u = 'dyna_grid_U.nc', &
         clname_v = 'dyna_grid_V.nc', &
         clname_w = 'dyna_grid_W.nc', &
         clname_s = 'dyna_wspd.nc'
      !
      ! 0. Initialization
      ! cas d'un fichier non periodique : on utilise deux fois le premier et
      ! le dernier champ temporel

      jkenr = kenr

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'Dynrea : reading dynamical fields, kenr = ', jkenr
         WRITE(numout,*) ' ~~~~~~~'
         WRITE(numout,*)
      ENDIF


      
      idtatot = ndtatot

      IF( kt == nit000 .AND. nlecoff == 0 ) THEN

         nlecoff = 1

         CALL flinopen(clname_t,mig(1),nlci,mjg(1),nlcj,.FALSE.,ipi,ipj, &
            &  ipk,zlon,zlat,zlev,itime,istep,zdate0,rdt,numfl_t)

         IF( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'problem with dimensions'
               WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               WRITE(numout,*) ' ipk ',ipk,' jpk    ',jpk
            ENDIF
            STOP 'dynrea  '
         ENDIF

         CALL flinopen(clname_u,mig(1),nlci,mjg(1),nlcj,.FALSE.,ipi,ipj, &
            &  ipk,zlon,zlat,zlev,itime,istep,zdate0,rdt,numfl_u)

         IF( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'problem with dimensions'
               WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               WRITE(numout,*) ' ipk ',ipk,' jpk    ',jpk
            ENDIF
            STOP 'dynrea  '
         ENDIF

         CALL flinopen(clname_v,mig(1),nlci,mjg(1),nlcj,.FALSE.,ipi,ipj, &
            &  ipk,zlon,zlat,zlev,itime,istep,zdate0,rdt,numfl_v)

         IF( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'problem with dimensions'
               WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               WRITE(numout,*) ' ipk ',ipk,' jpk    ',jpk
            ENDIF
            STOP 'dynrea '
         ENDIF

         CALL flinopen(clname_w,mig(1),nlci,mjg(1),nlcj,.FALSE.,ipi,ipj, &
            &  ipk,zlon,zlat,zlev,itime,istep,zdate0,rdt,numfl_w)

         IF( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'problem with dimensions'
               WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               WRITE(numout,*) ' ipk ',ipk,' jpk    ',jpk
            ENDIF
            STOP 'dynrea '
         ENDIF

         CALL flinopen(clname_s,mig(1),nlci,mjg(1),nlcj,.FALSE.,ipi,ipj, &
            &  ipk,zlon,zlat,zlev,itime,istep,zdate0,rdt,numfl_s)

         IF( ipi /= jpidta .OR. ipj /= jpjdta  ) THEN
            IF(lwp) THEN
               WRITE(numout,*)
               WRITE(numout,*) 'problem with dimensions'
               WRITE(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               WRITE(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
            ENDIF
            STOP 'dynrea'
         ENDIF

      ENDIF

      CALL flinget(numfl_u,'vozocrtx',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &         jkenr,mig(1),nlci,mjg(1),nlcj,zu(1:nlci,1:nlcj,1:jpk))

#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
      CALL flinget(numfl_u,'sobblcox',jpidta,jpjdta,1,idtatot,jkenr,  &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zbblx(1:nlci,1:nlcj))
#endif

# if defined key_traldf_eiv
      CALL flinget(numfl_u,'vozoeivu',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zaeiu(1:nlci,1:nlcj,1:jpk))
#endif

      CALL flinget(numfl_v,'vomecrty',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zv(1:nlci,1:nlcj,1:jpk))

#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
      CALL flinget(numfl_v,'sobblcoy',jpidta,jpjdta,1,idtatot,jkenr,  &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zbbly(1:nlci,1:nlcj))
#endif

# if defined key_traldf_eiv
      CALL flinget(numfl_v,'vomeeivv',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zaeiv(1:nlci,1:nlcj,1:jpk))
#endif

      CALL flinget(numfl_w,'vovecrtz',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zw(1:nlci,1:nlcj,1:jpk))

# if defined key_traldf_eiv
      CALL flinget(numfl_w,'voveeivw',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zaeiw(1:nlci,1:nlcj,1:jpk))
#endif


#if defined key_zdfddm
      CALL flinget(numfl_w,'voddmavs',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zavt(1:nlci,1:nlcj,1:jpk))
#else
      CALL flinget(numfl_w,'votkeavt',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zavt(1:nlci,1:nlcj,1:jpk))
#endif

#if   defined key_traldf_eiv   &&   defined key_traldf_c2d
      CALL flinget(numfl_w,'soleahtw',jpidta,jpjdta,1,idtatot,jkenr,   &
                  jkenr,mig(1),nlci,mjg(1),nlcj,zahtw(1:nlci,1:nlcj))

      CALL flinget(numfl_w,'soleaeiw',jpidta,jpjdta,1,idtatot,jkenr,   &
                  jkenr,mig(1),nlci,mjg(1),nlcj,zeivw(1:nlci,1:nlcj))
#endif

      CALL flinget(numfl_t,'votemper',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zt(1:nlci,1:nlcj,1:jpk))

      CALL flinget(numfl_t,'vosaline',jpidta,jpjdta,jpk,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zs(1:nlci,1:nlcj,1:jpk))

      CALL flinget(numfl_t,'somixhgt',jpidta,jpjdta,1,idtatot,jkenr,  &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zmld(1:nlci,1:nlcj))


      CALL flinget(numfl_t,'sowaflup',jpidta,jpjdta,1,idtatot,jkenr,  &
         &         jkenr,mig(1),nlci,mjg(1),nlcj,zemp(1:nlci,1:nlcj))

      CALL flinget(numfl_t,'soshfldo',jpidta,jpjdta,1,idtatot,jkenr,  &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zqsr(1:nlci,1:nlcj))

      CALL flinget(numfl_t,'soicecov',jpidta,jpjdta,1,idtatot,jkenr,  &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zice(1:nlci,1:nlcj))

      CALL flinget(numfl_s,'wspd',    jpidta,jpjdta,1,idtatot,jkenr,   &
         &        jkenr,mig(1),nlci,mjg(1),nlcj,zwind(1:nlci,1:nlcj))
 

        ! Extra-halo initialization in MPP
         IF( lk_mpp ) THEN
            DO ji = nlci+1, jpi
               zu(ji,:,:) = zu(1,:,:)   
               zv(ji,:,:) = zv(1,:,:)   
               zw(ji,:,:) = zw(1,:,:)   
               zavt(ji,:,:)=zavt(1,:,:)
               zt(ji,:,:)=zt(1,:,:)
               zs(ji,:,:)=zs(1,:,:)
               zmld(ji,:)=zmld(1,:)
               zwind(ji,:)=zwind(1,:)
               zemp(ji,:)=zemp(1,:)
               zqsr(ji,:)=zqsr(1,:)
               zice(ji,:)=zice(1,:)
#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
               zbblx(ji,:)=zbblx(1,:)
               zbbly(ji,:)=zbbly(1,:)
#endif
#if defined key_traldf_eiv
               zaeiu(ji,:,:)=zaeiu(1,:,:)
               zaeiv(ji,:,:)=zaeiv(1,:,:)
               zaeiw(ji,:,:)=zaeiw(1,:,:)
#endif
#if defined key_traldf_eiv   &&   defined key_traldf_c2d
               zahtw(ji,:)=zahtw(1,:)
               zeivw(ji,:)=zeivw(1,:)
#endif
            ENDDO
            DO jj = nlcj+1, jpj
               zu(:,jj,:) = zu(:,1,:)
               zv(:,jj,:) = zv(:,1,:)
               zw(:,jj,:) = zw(:,1,:)
               zavt(:,jj,:)=zavt(:,1,:)
               zt(:,jj,:)=zt(:,1,:)
               zs(:,jj,:)=zs(:,1,:)
               zmld(:,jj)=zmld(:,1)
               zwind(:,jj)=zwind(:,1)
               zemp(:,jj)=zemp(:,1)
               zqsr(:,jj)=zqsr(:,1)
               zice(:,jj)=zice(:,1)
#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
               zbblx(:,jj)=zbblx(:,1)
               zbbly(:,jj)=zbbly(:,1)
#endif
#if defined key_traldf_eiv
               zaeiu(:,jj,:)=zaeiu(:,1,:)
               zaeiv(:,jj,:)=zaeiv(:,1,:)
               zaeiw(:,jj,:)=zaeiw(:,1,:)
#endif
#if defined key_traldf_eiv   &&   defined key_traldf_c2d
               zahtw(:,jj)=zahtw(:,1)
               zeivw(:,jj)=zeivw(:,1)
#endif
            ENDDO
         ENDIF


            udta(:,:,:,2)=zu(:,:,:)*umask(:,:,:)
            vdta(:,:,:,2)=zv(:,:,:)*vmask(:,:,:)
            wdta(:,:,:,2)=zw(:,:,:)*tmask(:,:,:)
            tdta(:,:,:,2)=zt(:,:,:)*tmask(:,:,:)
            sdta(:,:,:,2)=zs(:,:,:)*tmask(:,:,:)
            avtdta(:,:,:,2)=zavt(:,:,:)*tmask(:,:,:)
#if defined key_traldf_eiv   &&   defined key_traldf_c2d
            ahtwdta(:,:,2)=zahtw(:,:)*tmask(:,:,1)
            eivwdta(:,:,2)=zeivw(:,:)*tmask(:,:,1)
#endif
      !
      !
      ! flux :
      !
            flxdta(:,:,jpwind,2)=zwind(:,:)*tmask(:,:,1)
            flxdta(:,:,jpice,2)=min(1.,zice(:,:))*tmask(:,:,1)
            flxdta(:,:,jpemp,2)=zemp(:,:)*tmask(:,:,1)
            flxdta(:,:,jpqsr,2)=zqsr(:,:)*tmask(:,:,1)
            zmxldta(:,:,2)=zmld(:,:)*tmask(:,:,1)

#if defined key_trcbbl_dif   ||   defined key_trcbbl_adv
            bblxdta(:,:,2)=max(0.,zbblx(:,:))
            bblydta(:,:,2)=max(0.,zbbly(:,:))

        DO ji=1,jpi
          DO jj=1,jpj
            if (bblxdta(ji,jj,2).gt.2.) bblxdta(ji,jj,2)=0.
            if (bblydta(ji,jj,2).gt.2.) bblydta(ji,jj,2)=0.
          END DO
        END DO
#endif

   END SUBROUTINE dynrea



END MODULE dtadyn
