!!DB/CN: 2008.11.14
!!Calls to IOIPSL eliminated and tested
!!More code cleanup likely possible, but not done

!!Key changes from standard obcdta.F90
!!(0) Look for  DB  to see various modifications

!!(1) obc_dta() has HARDWIRED areas where transports at certain boundary
!!locations are adjusted ---> key_ADJ_TRANSPORT; SBI_trans; SB_east; ...
!!- Note that at this time the code is not completely consistent with
!!respect to key_ADJ_TRANSPORT as some variables and operations exist that
!!are only required for that key but they are done even if that key is not defined.
!!These exceptions require minimal memory and execution time.

!!(2) Because OBC transports are adjusted, the eta-arrays need to be adjusted as
!!well to be consistent with the OBC vels (thermal wind is used). This is a major
!!change that requires new variables (search eta_). 
!!Note that 1 & 2 require definitions and computations on global variables defined
!!for these new OBC modifications. By global I mean that all relevant CPUs must
!!do the same computations and have the same values.

!!(3) The barotropic velocity in obc_dta_bt() is now computed as the vertical
!!integral of the baroclinic velocities. Thus the BT vels in the input files 
!!are not used anymore (they are ignored). Search (e.g.) ubtfoe0
!!(4) 5 tidal constituents are possible (search: nntide)
!!(5) tidal and non-tidal forcing are separated (search (e.g.) ubtfoe0, sshfoe0) 
!!(6) A routine is called that zeros the non-tidal net transport around the
!! boundaries (search: obc_ctl) 
!!(7) A ramp function is used for the forcing (search ramp). It is computed in
!!step.F90 and defined in oce.F90, It is based on kt as
!!opposed to (kt-nit000) so that by default it does not ramp a restart.  Also
!!it is complicated to follow when/where certain variables are ramped up -- but no
!!apologies for this.
!!(8) The routines are HARDWIRED for 3 open boundaries only. If north is also open, then
!!all the lp_obc_north areas should be checked and the code added as necessary. 
!!(9) Uses global masks: emaskg2, smaskg2, wmaskg2 (see obcini, obc_oce) 
!!(10) NB: use of e2v_e, e2u_s etc seem incorrect, but diff should be v.small so
!!     modification ... TO DO ... 


MODULE obcdta 
  !!==============================================================================
  !!                            ***  MODULE obcdta  ***
  !! Open boundary data : read the data for the open boundaries.
  !!==============================================================================
#if defined key_obc
  !!------------------------------------------------------------------------------
  !!   'key_obc'         :                                Open Boundary Conditions
  !!------------------------------------------------------------------------------
  !!   obc_dta           : read u, v, t, s data along each open boundary
  !!   obc_dta_psi       : read psi data along each open boundary (rigid lid only)
  !!------------------------------------------------------------------------------
  !! * Modules used
  USE oce             ! ocean dynamics and tracers
  USE dom_oce         ! ocean space and time domain
  USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
  USE phycst          ! physical constants
  USE obc_oce         ! ocean open boundary conditions
  USE daymod          ! calendar
  USE in_out_manager  ! I/O logical units
  USE lib_mpp         ! distributed memory computing
  USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
  USE ioipsl
  USE obcctl
!!CN
  USE lib_ncdf
!!DBG
#  if defined key_dynspg_rl
  USE obccli
#  endif

  IMPLICIT NONE
  PRIVATE

  !! * Accessibility
  PUBLIC obc_dta        ! routines called by step.F90
  PUBLIC obc_dta_bt     ! routines called by dynspg_ts.F90

  !! * Shared module variables
  INTEGER ::   &
     nlecto,   &  ! switch for the first read
     ntobc1,   &  ! first record used
     ntobc2,   &  ! second record used
     ntobc3,   &  ! last record used  for dynspg_ts
     itobc        ! number of time steps in OBC files


  REAL(wp), DIMENSION(:), ALLOCATABLE :: ztcobc      ! time_counter variable of BCs

  !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "obc_vectopt_loop_substitute.h90"

  !!---------------------------------------------------------------------------------
  !!   OPA 9.0 , LODYC-IPSL  (2003)
  !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdta.F90,v 1.9 2006/03/21 08:25:09 opalod Exp $
  !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
  !!---------------------------------------------------------------------------------

CONTAINS

  SUBROUTINE obc_dta (kt)
     !!--------------------------------------------------------------------
     !!              ***  SUBROUTINE obc_dta  ***
     !!
     !! ** Purpose :
     !!   Find the climatological boundary arrays for the specified date,
     !!   The boundary arrays are netcdf files. Three possible cases:
     !!   - one time frame only in the file (time dimension = 1).
     !!     in that case the boundary data does not change in time.
     !!   - many time frames. In that case,  if we have 12 frames
     !!     we assume monthly fields.
     !!     Else, we assume that time_counter is in seconds
     !!     since the beginning of either the current year or a reference
     !!     year given in the namelist.
     !!     (no check is done so far but one would have to check the "unit"
     !!     attribute of variable time_counter).
     !!
     !! History :
     !!        !  98-05 (J.M. Molines) Original code
     !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
     !!   9.0  !  04-06 (F. Durand, A-M. Treguier) Netcdf BC files on input
     !!--------------------------------------------------------------------
     !! * Arguments


     INTEGER, INTENT( in ) ::   kt          ! ocean time-step index

     !! * Local declarations
     INTEGER ::   ji, jj, jk, ii, ij   ! dummy loop indices
     INTEGER ::   itimo, iman, imois
     INTEGER ::   i15
     REAL(wp) ::   zxy
     !! * Ajouts FD
     INTEGER ::  isrel              ! number of seconds since 1/1/1992
     INTEGER, SAVE ::  itobce, itobcw,  & ! number of time steps in OBC files
                       itobcs, itobcn     !    "       "       "       "
     INTEGER ::  ikprint        ! frequency for printouts.
     INTEGER :: fid_e, fid_w, fid_n, fid_s       ! file identifiers
     LOGICAL :: l_exv
     INTEGER, DIMENSION(flio_max_dims) ::   f_d  ! dimensions lenght
     CHARACTER(LEN=25) :: v_name
!!CN
     INTEGER :: f_stat
     INTEGER :: len
     !!LOGICAL :: use_ioipsl = .FALSE.

!!DB
     LOGICAL :: DBG
     INTEGER :: ncells
     REAL(wp), DIMENSION(jpjdta,jpk) :: U_EW
!     REAL(wp), DIMENSION(jpidta,jpk) :: V_NS     !!for future use

!!DB: For boundary transport adjustments
!!NB: Could be inside of #ifdef key_ADJ_TRANSPORT
     REAL(wp) ::   fac_SBI, fac_SB_east, trans_in, trans_diff, trans_check, & 
                    vfreq, off1,a1,ph1, off2,a2,ph2  !!for time-dependent fac_*
     INTEGER :: jwest1, jwest2, jeast1, jeast2                !!indices where transport is adjusted 


     !!--------------------------------------------------------------------

!!DB
     DBG = .false. 

     IF( lk_dynspg_rl )  THEN
        CALL obc_dta_psi (kt)     ! update bsf data at open boundaries
        IF( nobc_dta == 1 .AND. kt == nit000 )   THEN
           IF(lwp) WRITE(numout,*) ' time-variable psi boundary data not allowed yet'
           STOP
        ENDIF
     ENDIF

     CALL ipslnlf (new_number=numout)

!!DB initialize time-dependent valve forcing parameters
!!Note that these are only used if key_ADJ_TRANSPORT is defined, but I
!!and lazy and do not #ifdef this

     vfreq = 2.0*rpi/365.0   ! seasonal cycle in days
!SBI params
!DL     off1 = 3.0
!DL     off1 = 2.375
!     a1= 1.625
!     off1 = 2.75
!     a1 = 2.0
     off1 = 1.975
     a1 = 1.225
!     off1 = 0.25
!     a1 = 0.1
     ph1 = 30.0   ! in days
!SB_east params: CLM0 values
! DL offset value reduced to 2.5 to reduce entry of warm water in the Gulf
! Jan 2013
!     off2 = 3.0
     off2 = 1.5

    if (kt.eq.nit000) then
      if (lwp) print*, '       '
      if (lwp) print*, 'off2 = ', off2
    endif 
     a2 = 0.5
     ph2 = -60.0   ! in days

     ! 1.   First call: check time frames available in files.
     ! -------------------------------------------------------

     IF( kt == nit000 )   THEN

        nlecto =  0

        IF(lwp) WRITE(numout,*)
        IF(lwp) WRITE(numout,*)     'obc_dta : find boundary data'
        IF(lwp) WRITE(numout,*)     '~~~~~~~'

        IF( nobc_dta == 0 )   THEN
           IF(lwp) WRITE(numout,*)  '  OBC data taken from initial conditions.'
           ntobc1 = 1
           ntobc2 = 1
        ELSE
           IF(lwp) WRITE(numout,*)  '  OBC data taken from netcdf files.'
           IF(lwp) WRITE(numout,*)  '  climatology (T/F):',ln_obc_clim
           ! check the number of time steps in the files.
           itobce =0 ; itobcw = 0; itobcn = 0; itobcs = 0
           v_name = 'time_counter'
           IF( lp_obc_east )   THEN
              !!CN: Replacing IOIPSL calls with lib_ncdf
              !CALL flioopfd ('obceast_TS.nc',fid_e)
              !CALL flioinqv (fid_e,TRIM(v_name),l_exv,len_dims=f_d)
              CALL ncdf_get_dim_size('obceast_TS.nc', 'time_counter', len, f_stat)
              IF( f_stat == 0 )   THEN
                 f_d(1) = len
                 itobce = f_d(1)
              ELSE
                 if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obceast_TS.nc'
              ENDIF
           ENDIF
           IF( lp_obc_west )   THEN
              !!CN: Replacing IOIPSL calls with lib_ncdf
              !CALL flioopfd ('obcwest_TS.nc',fid_w)
              !CALL flioinqv (fid_w,TRIM(v_name),l_exv,len_dims=f_d)
              CALL ncdf_get_dim_size('obcwest_TS.nc', 'time_counter', len, f_stat)
              IF( f_stat == 0 )   THEN
                 f_d(1) = len
                 itobcw = f_d(1)
              ELSE
                if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcwest_TS.nc'
              ENDIF
           ENDIF
           IF( lp_obc_north )   THEN
              !!CN: Replacing IOIPSL calls with lib_ncdf
              !CALL flioopfd ('obcnorth_TS.nc',fid_n)
              !CALL flioinqv (fid_n,TRIM(v_name),l_exv,len_dims=f_d)
              CALL ncdf_get_dim_size('obcnorth_TS.nc', 'time_counter', len, f_stat)
              IF( f_stat == 0 )   THEN
                 f_d(1) = len
                 itobcn = f_d(1)
              ELSE
                 if(lwp)WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcnorth_TS.nc'
              ENDIF
           ENDIF
           IF( lp_obc_south )   THEN
              !!CN: Replacing IOIPSL calls with lib_ncdf
              !CALL flioopfd ('obcsouth_TS.nc',fid_s)
              !CALL flioinqv (fid_s,TRIM(v_name),l_exv,len_dims=f_d)
              CALL ncdf_get_dim_size('obcsouth_TS.nc', 'time_counter', len, f_stat)
              IF( f_stat == 0 )   THEN
                 f_d(1) = len
                 itobcs = f_d(1)
              ELSE
                 if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcsouth_TS.nc'
              ENDIF
           ENDIF

           itobc = MAX(itobce,itobcw,itobcn,itobcs)
           nstop = 0
           IF( lp_obc_east  .AND. itobce /= itobc ) nstop = nstop+1
           IF( lp_obc_west  .AND. itobcw /= itobc ) nstop = nstop+1
           IF( lp_obc_north .AND. itobcn /= itobc ) nstop = nstop+1
           IF( lp_obc_south .AND. itobcs /= itobc ) nstop = nstop+1
           IF( nstop /= 0 )  THEN
              IF( lwp )   THEN
                 WRITE(numout,*) ' obcdta : all files must have the same number of time steps'
                 WRITE(numout,*) ' east, west, north, south: ', itobce, itobcw, itobcn, itobcs
              ENDIF
              STOP
           ENDIF
           IF( itobc == 1 )   THEN
              IF( lwp ) WRITE(numout,*) ' obcdta found one time step only in the OBC files'
           ELSE
              ALLOCATE (ztcobc(itobc))
              l_exv = .TRUE.
              IF( lp_obc_east )   THEN
                 !!CN: Replacing IOIPSL calls with lib_ncdf
                 IF( l_exv )   THEN
                    !CALL fliogetv (fid_e,TRIM(v_name),ztcobc)
                    CALL ncdf_read('obceast_TS.nc', 'time_counter', ztcobc, f_stat)
                    l_exv = .FALSE.
                 ENDIF
                 !CALL flioclo (fid_e)
              ENDIF
              IF( lp_obc_west )   THEN
                 !!CN: Replacing IOIPSL calls with lib_ncdf
                IF( l_exv )   THEN
                   !CALL fliogetv (fid_w,TRIM(v_name),ztcobc)
                   CALL ncdf_read('obcwest_TS.nc', 'time_counter', ztcobc, f_stat)
                   l_exv = .FALSE.
                ENDIF
                !CALL flioclo (fid_w)
              ENDIF
              IF( lp_obc_north )   THEN
                 !!CN: Replacing IOIPSL calls with lib_ncdf
                IF( l_exv )   THEN
                   !CALL fliogetv (fid_n,TRIM(v_name),ztcobc)
                   CALL ncdf_read('obcnorth_TS.nc', 'time_counter', ztcobc, f_stat)
                   l_exv = .FALSE.
                ENDIF
                !CALL flioclo (fid_n)
              ENDIF
              IF( lp_obc_south )   THEN
                 !!CN: replacing IOIPSL calls with lib_ncdf
                IF( l_exv )   THEN
                   !CALL fliogetv (fid_s,TRIM(v_name),ztcobc)
                   CALL ncdf_read('obcsouth_TS.nc', 'time_counter', ztcobc, f_stat)
                   l_exv = .FALSE.
                ENDIF
                !CALL flioclo (fid_s)
              ENDIF
              IF( lwp ) WRITE(numout,*) ' obcdta found', itobc,' time steps in the OBC files'
              IF( .NOT. ln_obc_clim .AND. itobc == 12 )   THEN
                 IF ( lwp ) WRITE(numout,*) '  WARNING: With monthly data we assume climatology'
                 ln_obc_clim = .true.
              ENDIF
           ENDIF
        ENDIF

      ! 1.1  Tangential velocities set to zero
      ! --------------------------------------
        IF( lp_obc_east  ) vfoe = 0.e0
        IF( lp_obc_west  ) vfow = 0.e0
        IF( lp_obc_south ) ufos = 0.e0
        IF( lp_obc_north ) ufon = 0.e0

      ! 1.2  Data temperature, salinity, normal velocities set to zero
      !                        or initial conditions if nobc_dta == 0
      ! --------------------------------------------------------------

        IF( lp_obc_east )   THEN
           ! initialisation to zero
!byoung           
           sedta(:,:,:,:) = 0.e0
           tedta(:,:,:,:) = 0.e0
           uedta(:,:,:) = 0.e0
           !                                    ! ================== !
           IF( nobc_dta == 0 )   THEN           ! initial state used
              !                                 ! ================== !
              !  Fills sedta, tedta, uedta (global arrays)
              !  Remark: this works for njzoom = 1.
              !          Should the definition of ij include njzoom?
              DO ji = nie0, nie1
                 DO jk = 1, jpkm1
                    DO jj = nje0p1, nje1m1
                       ij = jj -1 + njmpp
!byoung                       
                       sedta(ij,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                       tedta(ij,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                       uedta(ij,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
                    END DO
                 END DO
              END DO
!byoung
              DO ji = nie0-4, nie0-1
                 DO jk = 1, jpkm1
                    DO jj = nje0p1, nje1m1
                       ij = jj -1 + njmpp
                       sedta(ij,jk,nie0-ji+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                       tedta(ij,jk,nie0-ji+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                    END DO
                 END DO
              END DO              
           ENDIF
        ENDIF

        IF( lp_obc_west )   THEN
           ! initialisation to zero
!byoung           
           swdta(:,:,:,:) = 0.e0
           twdta(:,:,:,:) = 0.e0
           uwdta(:,:,:) = 0.e0
           !                                    ! ================== !
           IF( nobc_dta == 0 )   THEN           ! initial state used !
              !                                 ! ================== !
              !  Fills swdta, twdta, uwdta (global arrays)
              !  Remark: this works for njzoom = 1.
              !          Should the definition of ij include njzoom?
              DO ji = niw0, niw1
                 DO jk = 1, jpkm1
                    DO jj = njw0p1, njw1m1
                       ij = jj -1 + njmpp
!byoung                       
                       swdta(ij,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                       twdta(ij,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                       uwdta(ij,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
                    END DO
                 END DO
              END DO
!byoung
              DO ji = niw0+1, niw0+4
                 DO jk = 1, jpkm1
                    DO jj = njw0p1, njw1m1
                       ij = jj -1 + nimpp
                       swdta(ij,jk,ji-niw0+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                       twdta(ij,jk,ji-niw0+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                    END DO
                 END DO
              END DO
           ENDIF
        ENDIF

        IF( lp_obc_north)   THEN
           ! initialisation to zero

!byoung
           sndta(:,:,:,:) = 0.e0
           tndta(:,:,:,:) = 0.e0
           vndta(:,:,:) = 0.e0
           !                                    ! ================== !
           IF( nobc_dta == 0 )   THEN           ! initial state used
              !                                 ! ================== !
              !  Fills sndta, tndta, vndta (global arrays)
              !  Remark: this works for njzoom = 1.
              !          Should the definition of ij include njzoom?
              DO jj = njn0, njn1
                 DO jk = 1, jpkm1
                    DO ji = nin0p1, nin1m1
                       ii = ji -1 + nimpp
                       sndta(ii,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                       tndta(ii,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                       vndta(ii,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
                    END DO
                 END DO
              END DO
!sujie
              DO jj = njn0-4, njn0-1
                 DO jk = 1, jpkm1
                    DO ji = nin0p1, nin1m1
                       ii = ji -1 + nimpp
                       sndta(ii,jk,njn0-jj+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                       tndta(ii,jk,njn0-jj+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                    END DO
                 END DO
              END DO
           ENDIF
        ENDIF

        IF( lp_obc_south )   THEN
           ! initialisation to zero

           ssdta(:,:,:,:) = 0.e0
           tsdta(:,:,:,:) = 0.e0
           vsdta(:,:,:) = 0.e0
           !                                    ! ================== !
           IF( nobc_dta == 0 )   THEN           ! initial state used
              !                                 ! ================== !
              !  Fills ssdta, tsdta, vsdta (global arrays)
              !  Remark: this works for njzoom = 1.
              !          Should the definition of ij include njzoom?
              DO jj = njs0, njs1
                 DO jk = 1, jpkm1
                    DO ji = nis0p1, nis1m1
                       ii = ji -1 + nimpp

                       ssdta(ii,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                       tsdta(ii,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                       vsdta(ii,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
                    END DO
                 END DO
              END DO
!sujie
              DO jj = njs0+1, njs0+4
                 DO jk = 1, jpkm1
                    DO ji = nis0p1, nis1m1
                       ii = ji -1 + nimpp
                       ssdta(ii,jk,jj-njs0+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                       tsdta(ii,jk,jj-njs0+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                    END DO
                 END DO
              END DO
           ENDIF
        ENDIF


     ENDIF        !       end if kt == nit000

     ! 2.  Initialize the time we are at.
     !     Does this every time the routine is called,
     !     excepted when nobc_dta = 0
     !---------------------------------------------------------------------
     IF( nobc_dta == 0 )   THEN
        itimo = 1
        zxy   = 0
     ELSE
        IF( itobc == 1 )   THEN
           itimo = 1
        ELSE IF( itobc == 12 )   THEN      !   BC are monthly
           ! we assume we have climatology in that case
           iman  = 12
           i15   = nday / 16
           imois = nmonth + i15 - 1
           IF( imois == 0 )   imois = iman
           itimo = imois
        ELSE
           IF(lwp) WRITE(numout,*) 'data other than constant or monthly',kt
           iman  = itobc
           itimo = FLOOR( kt*rdt / (ztcobc(2)-ztcobc(1)) )
           isrel = kt*rdt
        ENDIF
     ENDIF

     ! 2.1 Read two records in the file if necessary
     ! ---------------------------------------------
     IF( ( nobc_dta == 1 ) .AND. ( ( kt == nit000 .AND. nlecto == 0 ) .OR. itimo  /= ntobc1 ) )   THEN
        nlecto = 1

        ! Calendar computation
        IF( itobc == 1 )   THEN            !  BC are constant in time
           ntobc1 = 1
           ntobc2 = 1
        ELSE IF( itobc == 12 )   THEN      !   BC are monthly
           ntobc1 = itimo         ! first file record used
           ntobc2 = ntobc1 + 1    ! last  file record used
           ntobc1 = MOD( ntobc1, iman )
           IF( ntobc1 == 0 )   ntobc1 = iman
           ntobc2 = MOD( ntobc2, iman )
           IF( ntobc2 == 0 )   ntobc2 = iman
           IF( lwp )   THEN
              WRITE(numout,*) ' read monthly obc first record file used ntobc1 ', ntobc1
              WRITE(numout,*) ' read monthly obc last  record file used ntobc2 ', ntobc2
           ENDIF
        ELSE
           isrel=kt*rdt
           ntobc1 = itimo         ! first file record used
           ntobc2 = ntobc1 + 1    ! last  file record used
           ntobc1 = MOD( ntobc1, iman )
           IF( ntobc1 == 0 )   ntobc1 = iman
           ntobc2 = MOD( ntobc2, iman )
           IF( ntobc2 == 0 )   ntobc2 = iman
           IF(lwp) WRITE(numout,*) ' read obc first record file used ntobc1 ', ntobc1
           IF(lwp) WRITE(numout,*) ' read obc last  record file used ntobc2 ', ntobc2
        ENDIF
                              ! ======================= !
                              !  BCs read               !
        !                     ! ======================= !


        IF( lp_obc_east )   THEN
           ! ... Read datafile and set temperature, salinity and normal velocity
           ! ... initialise the sedta, tedta, uedta arrays
!!DB/CN: Replace IOIPSL calls with lib_ncdf
           CALL obc_dta_gv ('y','vosaline',jpjef-jpjed+1,ntobc1,pdta_4D=sedta(jpjed:jpjef,:,1:5,1),fn='obceast_TS.nc')
           CALL obc_dta_gv ('y','vosaline',jpjef-jpjed+1,ntobc2,pdta_4D=sedta(jpjed:jpjef,:,1:5,2),fn='obceast_TS.nc')
           CALL obc_dta_gv ('y','votemper',jpjef-jpjed+1,ntobc1,pdta_4D=tedta(jpjed:jpjef,:,1:5,1),fn='obceast_TS.nc')
           CALL obc_dta_gv ('y','votemper',jpjef-jpjed+1,ntobc2,pdta_4D=tedta(jpjed:jpjef,:,1:5,2),fn='obceast_TS.nc')
           CALL obc_dta_gv ('y','vozocrtx',jpjef-jpjed+1,ntobc1,pdta_3D=uedta(jpjed:jpjef,:,1),fn='obceast_U.nc')
           CALL obc_dta_gv ('y','vozocrtx',jpjef-jpjed+1,ntobc2,pdta_3D=uedta(jpjed:jpjef,:,2),fn='obceast_U.nc')

           !  Usually printout is done only once at kt = nit000,
           !  unless nprint (namelist) > 1 !!DB -- code deleted, see older version

        ENDIF

        IF( lp_obc_west )   THEN
           ! ... Read datafile and set temperature, salinity and normal velocity
           ! ... initialise the swdta, twdta, uwdta arrays
!!DB/CN: Replace IOIPSL calls with lib_ncdf
           CALL obc_dta_gv ('y','vosaline',jpjwf-jpjwd+1,ntobc1,pdta_4D=swdta(jpjwd:jpjwf,:,1:5,1),fn='obcwest_TS.nc')
           CALL obc_dta_gv ('y','vosaline',jpjwf-jpjwd+1,ntobc2,pdta_4D=swdta(jpjwd:jpjwf,:,1:5,2),fn='obcwest_TS.nc')
           CALL obc_dta_gv ('y','votemper',jpjwf-jpjwd+1,ntobc1,pdta_4D=twdta(jpjwd:jpjwf,:,1:5,1),fn='obcwest_TS.nc')
           CALL obc_dta_gv ('y','votemper',jpjwf-jpjwd+1,ntobc2,pdta_4D=twdta(jpjwd:jpjwf,:,1:5,2),fn='obcwest_TS.nc')
           CALL obc_dta_gv ('y','vozocrtx',jpjwf-jpjwd+1,ntobc1,pdta_3D=uwdta(jpjwd:jpjwf,:,1),fn='obcwest_U.nc')
           CALL obc_dta_gv ('y','vozocrtx',jpjwf-jpjwd+1,ntobc2,pdta_3D=uwdta(jpjwd:jpjwf,:,2),fn='obcwest_U.nc')
!!DB: printout code deleted
        ENDIF

        IF( lp_obc_north )   THEN
!!DB/CN: Replace IOIPSL calls with lib_ncdf
           CALL obc_dta_gv ('x','vosaline',jpinf-jpind+1,ntobc1,pdta_4D=sndta(jpind:jpinf,:,1:5,1),fn='obcnorth_TS.nc')
           CALL obc_dta_gv ('x','vosaline',jpinf-jpind+1,ntobc2,pdta_4D=sndta(jpind:jpinf,:,1:5,2),fn='obcnorth_TS.nc')
           CALL obc_dta_gv ('x','votemper',jpinf-jpind+1,ntobc1,pdta_4D=tndta(jpind:jpinf,:,1:5,1),fn='obcnorth_TS.nc')
           CALL obc_dta_gv ('x','votemper',jpinf-jpind+1,ntobc2,pdta_4D=tndta(jpind:jpinf,:,1:5,2),fn='obcnorth_TS.nc')
           CALL obc_dta_gv ('x','vomecrty',jpinf-jpind+1,ntobc1,pdta_3D=vndta(jpind:jpinf,:,1),fn='obcnorth_V.nc')
           CALL obc_dta_gv ('x','vomecrty',jpinf-jpind+1,ntobc2,pdta_3D=vndta(jpind:jpinf,:,2),fn='obcnorth_V.nc')
        ENDIF

        IF( lp_obc_south )   THEN
!!DB/CN: Replace IOIPSL calls with lib_ncdf
           CALL obc_dta_gv ('x','vosaline',jpisf-jpisd+1,ntobc1,pdta_4D=ssdta(jpisd:jpisf,:,1:5,1),fn='obcsouth_TS.nc')
           CALL obc_dta_gv ('x','vosaline',jpisf-jpisd+1,ntobc2,pdta_4D=ssdta(jpisd:jpisf,:,1:5,2),fn='obcsouth_TS.nc')
           CALL obc_dta_gv ('x','votemper',jpisf-jpisd+1,ntobc1,pdta_4D=tsdta(jpisd:jpisf,:,1:5,1),fn='obcsouth_TS.nc')
           CALL obc_dta_gv ('x','votemper',jpisf-jpisd+1,ntobc2,pdta_4D=tsdta(jpisd:jpisf,:,1:5,2),fn='obcsouth_TS.nc')
           CALL obc_dta_gv ('x','vomecrty',jpisf-jpisd+1,ntobc1,pdta_3D=vsdta(jpisd:jpisf,:,1),fn='obcsouth_V.nc')
           CALL obc_dta_gv ('x','vomecrty',jpisf-jpisd+1,ntobc2,pdta_3D=vsdta(jpisd:jpisf,:,2),fn='obcsouth_V.nc')
        ENDIF

     ELSE

        nlecto = 0        !      no reading of OBC barotropic data

     ENDIF                !      end of the test on the condition to read or not the files

     ! 3.  Call at every time step :
     !     Linear interpolation of BCs to current time step
     ! ----------------------------------------------------

     IF( itobc == 1 .OR. nobc_dta == 0 )   THEN
        zxy = 0.
     ELSE IF( itobc == 12 )   THEN
        zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
     ELSE
        zxy = (ztcobc(ntobc1)-FLOAT(isrel))/(ztcobc(ntobc1)-ztcobc(ntobc2))
     ENDIF

 
!!DB: ramp for forcing. NB choice is kt rather than (kt-nit000) so the implicit 
!!assumption is that ramped forcing is not desirable for restarts.
!!NB: ramp now computed in step.F90
!     ramp=tanh(kt*rdt/(2.0*86400.0))
 
     IF( lp_obc_east )   THEN


        DO jk = 1, jpkm1
           DO jj = nje0p1, nje1m1
              ij = jj -1 + njmpp
!!DB
              sfoe(jj,jk,:) =  ( zxy * sedta(ij,jk,:,2) + &
                 &           (1.-zxy) * sedta(ij,jk,:,1) ) * temsk5(jj,jk,:) 
              tfoe(jj,jk,:) =  ( zxy * tedta(ij,jk,:,2) + &
                 &           (1.-zxy) * tedta(ij,jk,:,1) ) * temsk5(jj,jk,:)

           END DO
        END DO
!AD   interpolate to global
        DO jk = 1, jpkm1
           DO jj = jpjed,jpjef
              uedta1(jj,jk) =  ( zxy * uedta(jj,jk,2) + &
                 &           (1.-zxy) * uedta(jj,jk,1) ) !* global mask not * uemsk(jj,jk)
           END DO
        END DO

#ifdef key_ADJ_TRANSPORT
!!DB: adjust additional transport here, as applicable
!!DBG
        if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_east: ADJUSTING TRANSPORT, narea= ',narea
 

!!SBI Region
!!DB07.30: reverse part of SBI inflow as the Nfld coastal current is in wrong direction
        do jk = 1, jpkm1
           do jj = 222,224
              uedta1(jj,jk) = -uedta1(jj,jk)
           enddo
           uedta1(225,jk) = 0.0
        enddo

!!DB 2007.11.26
        fac_SBI = off1 + a1*cos(vfreq*(nday_year+ph1))
        if(lwp .AND. kt-nit000 == 0) write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SBI = ', &
                    nday_year,fac_SBI

        do jk = 1, jpkm1
           do jj = 220,230
              uedta1(jj,jk) = fac_SBI*uedta1(jj,jk)
           enddo
        enddo

!DB  East Inflow: Shelf break region ---------------------------------
        jeast1 = 102
        jeast2 = 108
        fac_SB_east = off2 + a2*cos(vfreq*(nday_year+ph2))
        if(lwp .AND. kt-nit000 == 0) &
             & write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SB_east = ', nday_year, fac_SB_east
        if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0)write(numout2,'(A51,2x,2(i7,1x),2(f8.3,1x))') &
          & 'DBG: obcdta -- kt,nday_year, fac_SBI, fac_SB_east ', kt,nday_year,fac_SBI,fac_SB_east

!!DB: globalize this calc
        ji = nie0
        ncells = 0      
        trans_in = 0.0
        do jk = 1, jpkm1  
           do jj = jeast1, jeast2
              trans_in = trans_in + uedta1(jj,jk)*e2v_e(jj)*e3t(jk)*emaskg2(jj,jk)
              ncells = ncells + emaskg2(jj,jk)
           enddo
        enddo

!!DB adjust SB_east inflow 
!!Note that I do not ramp U_EW as we need a correct uedta1 and correct ramping is done elsewhere
        trans_diff = (-fac_SB_east*1.e6) - trans_in
        do jk = 1, jpkm1  
           do jj = jeast1, jeast2
              U_EW(jj,jk) = (trans_diff)/ (e2v_e(jj) * e3t(jk)) * emaskg2(jj,jk)
           enddo
        enddo

	if (ncells .gt. 0) then
!AD/DB
           trans_check = 0.0
           do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
              do jj = jeast1, jeast2
                 uedta1(jj,jk) = uedta1(jj,jk) + U_EW(jj,jk)/float(ncells)
                 trans_check = trans_check + uedta1(jj,jk)*e2v_e(jj)*e3t(jk)*emaskg2(jj,jk)
              enddo
           enddo

!!DBG: check above calcs
           if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0) then
              write(numout2,'(A77,2x,2(i7,1x),4(f8.3,1x))') &
          & 'DBG: obcdta east -- kt,nday_year,fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
          &            fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
           endif
           if(DBG) then
              write(200+narea,'(A77,2x,3(i5,1x),4(f8.3,1x))') &
          & 'obcdta east -- kt,nday_year,ncells, fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
          &            ncells, fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
           endif

        end if
        
!DB0731:
!!modify the input at that Nfld bay that looks like a false open boundary
!!NB: for Barotropic tide cancel these mean transports
!AD

        do jk = 1, jpkm1  
           do jj = 140,143
              uedta1(jj,jk) = 0.0
           enddo
           do jj = 144, 147
!!DB 2009.06.26 -- reduce this transport to see if/how it affect transport into SW GSL
              uedta1(jj,jk) = -0.10 * emaskg2(jj,jk)
!              uedta1(jj,jk) = -0.02 * emaskg2(jj,jk)
           enddo
        enddo

#endif   !ADJ_TRANSPORTS EAST

!AD   Map global to local
        DO jk = 1, jpkm1
           DO jj = nje0p1, nje1m1
              ij = jj -1 + njmpp
              ufoe(jj,jk) =  uedta1(ij,jk)* uemsk(jj,jk)
           END DO
        END DO

!AD >>>: recompute elevation
        eta_e(jpjed)=0. !zero eta at SW corner
	DO ij=jpjed+1,jpjef
      	   eta_e(ij) = eta_e(ij-1) - e2v_e(ij-1)*ff_e(ij)/grav * &
           &    uedta1(ij-1,1)
	END DO
 !AD <<<
!!DBG
        if(DBG .AND. kt == nit000) then
           write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
                'DBG: kt,narea,max,min ufoe,uedta1,eta_e ', &
                kt,narea,MAXVAL(ufoe(:,:)),MINVAL(ufoe(:,:)), &
                MAXVAL(uedta1(:,:)),MINVAL(uedta1(:,:)), &
                MAXVAL(eta_e(:)),MINVAL(eta_e(:))
           DO jk = 1, jpkm1
              DO jj = nje0p1, nje1m1
                 ij = jj -1 + njmpp
                 write(3000+narea,'(3(i3,1x),3(e12.6,1x))') &
                      jk,jj,ij,ufoe(jj,jk),uedta1(ij,jk),uemsk(jj,jk)
              END DO
           END DO
           do ij = jpjed,jpjef
              write(8000+narea,'(i3,2x,2(e12.6,2x))')ij,eta_e(ij),uedta1(ij,1)
           enddo

        endif


        ufoe(:,:) = ramp*ufoe(:,:)

     ENDIF


     IF( lp_obc_west )   THEN

        DO jk = 1, jpkm1
           DO jj = njw0p1, njw1m1
              ij = jj -1 + njmpp
!!DB
              sfow(jj,jk,:) =  ( zxy * swdta(ij,jk,:,2) + &
                 &           (1.-zxy) * swdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
              tfow(jj,jk,:) =  ( zxy * twdta(ij,jk,:,2) + &
                 &           (1.-zxy) * twdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
           END DO
        END DO
!AD   interpolate to global
        DO jk = 1, jpkm1
           DO jj = jpjwd,jpjwf
              uwdta1(jj,jk) =  ( zxy * uwdta(jj,jk,2) + &
                 &           (1.-zxy) * uwdta(jj,jk,1) ) !* global mask not * uemsk(jj,jk)
           END DO
        END DO


#ifdef key_ADJ_TRANSPORT
!!DB: add additional transport here, as applicable
!!DBG
        if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_west ADJUSTING TRANSPORT, narea= ',narea

!DB  west Inflow---------------------------------
!!DBG
        if(lwp .AND. kt-nit000 == 0) &
             &  write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SB_west = fac_SB_east = ',nday_year, fac_SB_east

        jwest1 = 20
        jwest2 = 35
        ji = niw0
        ncells = 0     
        trans_in = 0.0
        do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
           do jj = jwest1,jwest2
              trans_in = trans_in + uwdta1(jj,jk)*e2v_w(jj)*e3t(jk)*wmaskg2(jj,jk)
              ncells = ncells + wmaskg2(jj,jk)
           enddo
        enddo

!!DB adjust SB_east inflow (NB:  -fac_SB_east*1.e6 == inflow in Sv)
!!Note that I do not ramp U_EW as we need a correct uedta1 and correct ramping is done elsewhere
        trans_diff = (-fac_SB_east*1.e6) - trans_in
        do jk = 1, jpkm1  
           do jj = jwest1, jwest2
              U_EW(jj,jk) = (trans_diff)/ (e2v_w(jj) * e3t(jk)) * wmaskg2(jj,jk)
           enddo
        enddo

	
	if (ncells .gt. 0) then
!AD/DB
           trans_check = 0.0
           do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
              do jj = jwest1, jwest2
                 uwdta1(jj,jk) = uwdta1(jj,jk) + U_EW(jj,jk)/float(ncells)
                 trans_check = trans_check + uwdta1(jj,jk)*e2v_w(jj)*e3t(jk)*wmaskg2(jj,jk)
              enddo
           enddo
!!DBG: check above calcs
           if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0) then
              write(numout2,'(A77,2x,2(i7,1x),4(f8.3,1x))') &
          & 'DBG: obcdta west -- kt,nday_year, fac_SB_east,trans_in,trans_diff,trans_check ', kt,nday_year, &
          &            fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6

           endif
           if(DBG) then
              write(200+narea,'(A77,2x,3(i5,1x),4(f8.3,1x))') &
          & 'obcdta west -- kt,nday_year,ncells, fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
          &            ncells, fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
           endif

        end if
#endif     !ADJ_TRANSPORT WEST

!AD   Map global to local
        DO jk = 1, jpkm1
           DO jj = njw0p1, njw1m1
              ij = jj -1 + njmpp
              ufow(jj,jk) =  uwdta1(ij,jk)* uwmsk(jj,jk)
           END DO
        END DO

!AD >>>: recompute elevation
        eta_w(jpjwd)=0. !zero eta at SW corner
	DO ij=jpjwd+1,jpjwf
      	   eta_w(ij) = eta_w(ij-1) - e2v_w(ij-1)*ff_w(ij)/grav * &
           &    uwdta1(ij-1,1)
	END DO
 !AD <<<

!!DBG
        if(DBG .AND. kt == nit000) then
           write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
                'DBG: kt,narea,max,min ufow,uwdta1,eta_w ', &
                kt,narea,MAXVAL(ufow(:,:)),MINVAL(ufow(:,:)), &
                MAXVAL(uwdta1(:,:)),MINVAL(uwdta1(:,:)), &
                MAXVAL(eta_w(:)),MINVAL(eta_w(:))
           DO jk = 1, jpkm1
              DO jj = njw0p1, njw1m1
                 ij = jj -1 + njmpp
                 write(4000+narea,'(3(i3,1x),3(e12.6,1x))') &
                      jk,jj,ij,ufow(jj,jk),uwdta1(ij,jk),uwmsk(jj,jk)
              END DO
           END DO
           do ij = jpjwd,jpjwf
              write(8100+narea,'(i3,2x,2(e12.6,2x))')ij,eta_w(ij),uwdta1(ij,1)
           enddo

        endif

        ufow(:,:) = ramp * ufow(:,:)

     ENDIF
!------------------------------------

!!DB: north is closed so do nothing. Note that this could lead to future probs
     IF( lp_obc_north )   THEN
        !  fills sfon, tfon, vfon (local to each processor)
        DO jk = 1, jpkm1
           DO ji = nin0p1, nin1m1
              ii = ji -1 + nimpp
!sujie
              sfon(ji,jk,:) =  ( zxy * sndta(ii,jk,:,2) + &
                 &           (1.-zxy) * sndta(ii,jk,:,1) )* tnmsk5(ji,jk,:)
              tfon(ji,jk,:) =  ( zxy * tndta(ii,jk,:,2) + &
                 &           (1.-zxy) * tndta(ii,jk,:,1) )* tnmsk5(ji,jk,:)
              vfon(ji,jk) =  ( zxy * vndta(ii,jk,2) + &
                 &           (1.-zxy) * vndta(ii,jk,1) ) * vnmsk(ji,jk)
           END DO
        END DO
     ENDIF

     IF( lp_obc_south )   THEN

        DO jk = 1, jpkm1
          DO ji = nis0p1, nis1m1
             ii = ji -1 + nimpp

             sfos(ji,jk,:) = ( zxy * ssdta(ii,jk,:,2) + &
                &          (1.-zxy) * ssdta(ii,jk,:,1) )* tsmsk5(ji,jk,:)
             tfos(ji,jk,:) = ( zxy * tsdta(ii,jk,:,2) + &
                &          (1.-zxy) * tsdta(ii,jk,:,1) )* tsmsk5(ji,jk,:)
   
            END DO
        END DO

!AD >>> interpolate to global
        DO jk = 1, jpkm1
          DO ii = jpisd,jpisf
             vsdta1(ii,jk) = ( zxy * vsdta(ii,jk,2) + &
                &          (1.-zxy) * vsdta(ii,jk,1) ) !* global mask not * vsmsk(ji,jk)       
  
            END DO
        END DO
!AD <<<

#ifdef key_ADJ_TRANSPORT
!!DB: Adjust transport here if desired
!        if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_south ADJUSTING TRANSPORT'
#endif    !ADJ_TRANSPORT SOUTH

!   Map global to local
        DO jk = 1, jpkm1
          DO ji = nis0p1, nis1m1
             ii = ji -1 + nimpp
             vfos(ji,jk) = vsdta1(ii,jk) * vsmsk(ji,jk)  
            END DO
        END DO

!  recompute elevation
        eta_s(jpisd)=0. !zero eta at SW corner
	DO ij=jpisd+1,jpisf
      	   eta_s(ij) = eta_s(ij-1) + e2u_s(ij-1)*ff_s(ij)/grav * &
           &    vsdta1(ij-1,1)
	END DO

!!DBG
        if(DBG .AND. kt == nit000) then
           write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
                'DBG: kt,narea,max,min vfos,vsdta1,eta_s ', &
                kt,narea,MAXVAL(vfos(:,:)),MINVAL(vfos(:,:)), &
                MAXVAL(vsdta1(:,:)),MINVAL(vsdta1(:,:)), &
                MAXVAL(eta_s(:)),MINVAL(eta_s(:))
           DO jk = 1, jpkm1
              DO ji = nis0p1, nis1m1
                 ii = ji -1 + nimpp
                 write(5000+narea,'(3(i3,1x),3(e12.6,1x))') &
                      jk,ji,ii,vfos(ji,jk),vsdta1(ii,jk), vsmsk(ji,jk)  
              END DO
           END DO
           do ij = jpisd,jpisf
              write(8200+narea,'(i3,2x,2(e12.6,2x))')ij,eta_s(ij),vsdta1(ij,1)
           enddo



        endif


        vfos(:,:) = ramp * vfos(:,:)

     ENDIF

!!DB: Add offset to East boundary SSH -- HARDWIRED
     do ij = jpjed+1, jpjef
        eta_e(ij) = eta_e(ij) + eta_s(jpisf)
     enddo


!!DB: call ZW routine that zeros non-tidal transports along obc's
     call obc_ctl(kt)


  END SUBROUTINE obc_dta


# if defined key_dynspg_rl
  !!-----------------------------------------------------------------------------
  !!   Rigid-lid
  !!-----------------------------------------------------------------------------

  SUBROUTINE obc_dta_psi ( kt )
     !!-----------------------------------------------------------------------------
     !!                       ***  SUBROUTINE obc_dta_psi  ***
     !!
     !! ** Purpose :
     !!      Update the climatological streamfunction OBC at each time step.
     !!      Depends on the user's configuration.  Here data are read only once
     !!      at the beginning of the run.
     !!
     !! ** Method :
     !!      1. initialization
     !!         kbsfstart: number of time steps over which increase bsf
     !!         during initialization. This is provided for a smooth start
     !!         in cases where the transport is large (like on the Antarctic
     !!         continent). also note that when kbfstart=1, the transport
     !!         increases a lot in one time step and the precision usually
     !!         required for the solver may not be enough.
     !!      2. set the time evolution of the climatological barotropic streamfunction
     !!         along the isolated coastlines ( gcbic(jnic) ).
     !!      3. set the climatological barotropic streamfunction at the boundary.
     !!
     !!      The last two steps are done only at first step (nit000) or if kt <= kbfstart
     !!
     !! History :
     !!        ! 97-08 (G. Madec, J.M. Molines)
     !!   8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
     !!   9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
     !!----------------------------------------------------------------------------
     !! * Arguments
     INTEGER, INTENT( in ) ::   kt          ! ocean time-step index

     !! * Local declarations
     INTEGER ::   ji, jj, jnic, jip         ! dummy loop indices
     INTEGER ::   inum = 11                 ! temporary logical unit
     INTEGER ::   ip, ii, ij, iii, ijj
     INTEGER ::   kbsfstart
     REAL(wp) ::   zsver1, zsver2, zsver3, z2dtr, zcoef
     !!----------------------------------------------------------------------------

     ! 1. initialisation
     ! -----------------

     kbsfstart =  1
     zsver1 =  bsfic0(1)
     zsver2 =  zsver1
     IF( kt <= kbsfstart )   THEN
        zcoef = float(kt)/float(kbsfstart)
     ELSE
        zcoef = 1.
     END IF
     bsfic(1) = zsver1*zcoef
     IF( lwp .AND. ( kt <= kbsfstart ) )   THEN
        IF(lwp) WRITE(numout,*)'            '
        IF(lwp) WRITE(numout,*)'obcdta: spinup phase in obc_dta_psi routine'
        IF(lwp) WRITE(numout,*)'~~~~~~  it=',kt,'  OBC: spinup coef: ', &
                                          zcoef, ' and transport: ',bsfic(1)
     END IF

     zsver2 =  bsfic(1)-bsfic(2)
     zsver3 =  bsfic(2)

     ! 2. Right hand side of the barotropic elliptic equation (isolated coastlines)
     ! ----------------------------------------------------------------------------

     IF( ( neuler == 0 ) .AND. ( kt == nit000 ) )   THEN
        z2dtr = 1./rdt
     ELSE
        z2dtr = 1./2./rdt
     END IF
     ! ... bsfb(ii,ij) should be constant but due to the Asselin filter it
     ! ... converges asymptotically towards bsfic(jnic)
     ! ... However, bsfb(ii,ij) is constant along the same coastlines
     ! ... ---> can be improved using an extra array for storing bsficb (before)
     IF( nbobc > 1 )   THEN
        DO jnic = 1,nbobc - 1
           gcbic(jnic) = 0.e0
           ip=mnic(0,jnic)
           DO jip = 1,ip
              ii = miic(jip,0,jnic)
              ij = mjic(jip,0,jnic)
              IF( ii >= nldi+ nimpp - 1 .AND. ii <= nlei+ nimpp - 1 .AND. &
                  ij >= nldj+ njmpp - 1 .AND. ij <= nlej+ njmpp - 1 )   THEN
                 iii=ii-nimpp+1
                 ijj=ij-njmpp+1
                 gcbic(jnic) = ( bsfic(jnic) - bsfb(iii,ijj) ) * z2dtr
              END IF
           END DO
        END DO
     END IF

     IF( lk_mpp )   CALL mpp_isl( gcbic, 3 )

     ! 3. Update the climatological barotropic function at the boundary
     ! ----------------------------------------------------------------

     IF( lpeastobc )   THEN

        IF( kt == nit000 .OR. kt <= kbsfstart )   THEN
           OPEN(inum,file='obceastbsf.dta')
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*) (bfoe(jj),jj=jpjed, jpjef)
           CLOSE(inum)
        END IF
        DO jj=jpjed, jpjefm1
           bfoe(jj)=bfoe(jj)*zcoef
        END DO

     END IF

     IF( lpwestobc)   THEN

        IF( kt == nit000 .OR. kt <= kbsfstart ) then
           OPEN(inum,file='obcwestbsf.dta')
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*) (bfow(jj),jj=jpjwd, jpjwf)
           CLOSE(inum)
        END IF
        DO jj=jpjwd, jpjwfm1
           bfow(jj)=bfow(jj)*zcoef
        END DO

     END IF

     IF( lpsouthobc)   THEN

        IF( kt == nit000 .OR. kt <= kbsfstart )   THEN
           OPEN(inum,file='obcsouthbsf.dta')
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*) (bfos(jj),jj=jpisd, jpisf)
           CLOSE(inum)
        END IF
        DO ji=jpisd, jpisfm1
           bfos(ji)=bfos(ji)*zcoef
        END DO

     END IF

     IF( lpnorthobc)   THEN
        IF( kt == nit000 .OR. kt <= kbsfstart )   THEN
           OPEN(inum,file='obcnorthbsf.dta')
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*)
           READ(inum,*) (bfon(jj),jj=jpind, jpinf)
           CLOSE(inum)
        END IF
        DO ji=jpind, jpinfm1
           bfon(ji)=bfon(ji)*zcoef
        END DO

     END IF

  END SUBROUTINE obc_dta_psi
#else
  !!-----------------------------------------------------------------------------
  !!   Default option
  !!-----------------------------------------------------------------------------
  SUBROUTINE obc_dta_psi ( kt )       ! Empty routine
     !! * Arguments
     INTEGER,INTENT(in) :: kt
     if(lwp) WRITE(numout,*) 'obc_dta_psi: You should not have seen this print! error?', kt
  END SUBROUTINE obc_dta_psi
# endif


#if defined key_dynspg_ts || defined key_dynspg_exp
  SUBROUTINE obc_dta_bt( kt, kbt )
     !!---------------------------------------------------------------------------
     !!                      ***  SUBROUTINE obc_dta  ***
     !!
     !! ** Purpose :   time interpolation of barotropic data for time-splitting scheme
     !!                Data at the boundary must be in m2/s
     !!
     !! History :
     !!   9.0  !  05-11 (V. garnier) Original code
     !!---------------------------------------------------------------------------
     !! * Arguments
     INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
     INTEGER, INTENT( in ) ::   kbt         ! barotropic ocean time-step index

     !! * Local declarations
     INTEGER ::   ji, jj, jk, ii, ij   ! dummy loop indices
     INTEGER ::   fid_e, fid_w, fid_n, fid_s, fid  ! file identifiers
     INTEGER ::   itimo, iman, imois, i15
     INTEGER ::   ntobcm, ntobcp, itimom, itimop
     REAL(wp) ::  zxy
     INTEGER ::   isrel, ikt           ! number of seconds since 1/1/1992
     INTEGER ::   ikprint              ! frequency for printouts.
!ylu/DB
     INTEGER ::   ntide, nntide   
     REAL(wp), DIMENSION(5) :: tideperiod    ! in hours
!!DB
!!Local arrays for BT transport, dimensioned as in obc_oce
     REAL(wp), DIMENSION(1:jpisf) ::    &   !:
     vbtfos0, sshfos0       !: south boundary barotropic transport without the tide
     REAL(wp), DIMENSION(1:jpjwf) ::   &  !:
     ubtfow0, sshfow0       !: west boundary barotropic transport without the tide
     REAL(wp), DIMENSION(1:jpjef) ::   &  !:
     ubtfoe0, sshfoe0       !: east boundary barotropic transport without the tide
     INTEGER :: status

!!DB: 1 tide= M2 = good for climatology runs; 5 tides ~ operational runs
     nntide =5
     !IF( kt == nit000 )print*,"Warning: nntide=5"
     !nntide =1
     tideperiod(1) = 12.4206  !M2
     tideperiod(2) = 12.0000  !S2
     tideperiod(3) = 25.8193  !O1
     tideperiod(4) = 23.9345  !K1
     tideperiod(5) = 12.6583  !N2


     !!---------------------------------------------------------------------------

     ! 1.   First call: check time frames available in files.
     ! -------------------------------------------------------

     IF( kt == nit000 ) THEN

        ! 1.1  Barotropic tangential velocities set to zero
        ! -------------------------------------------------
        IF( lp_obc_east  ) vbtfoe(:) = 0.e0
        IF( lp_obc_west  ) vbtfow(:) = 0.e0
        IF( lp_obc_south ) ubtfos(:) = 0.e0
        IF( lp_obc_north ) ubtfon(:) = 0.e0

        ! 1.2  Sea surface height and normal barotropic velocities set to zero
        !                               or initial conditions if nobc_dta == 0
        ! --------------------------------------------------------------------

         IF( lp_obc_east ) THEN
            ! initialisation to zero
            sshedta(:,:) = 0.e0
            ubtedta(:,:) = 0.e0
            !                                        ! ================== !
            IF( nobc_dta == 0 )   THEN               ! initial state used !
               !                                     ! ================== !
               !  Fills sedta, tedta, uedta (global arrays)
               !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
               DO ji = nie0, nie1
                  DO jj = nje0p1, nje1m1
                     ij = jj -1 + njmpp
                     sshedta(ij,1) = sshn(ji+1,jj) * tmask(ji+1,jj,1)
                  END DO
               END DO
            ENDIF
!ylu
           tidesshemag(:,:)=0.0
           tidesshepha(:,:)=0.0
           tidevbtemag(:,:)=0.0
           tidevbtepha(:,:)=0.0
!!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
           call ncdf_read_global('obceast_tide.nc','tidesshmag',tidesshemag(jpjed:jpjef,1:nntide),-1,status)
           call ncdf_read_global('obceast_tide.nc','tidesshpha',tidesshepha(jpjed:jpjef,1:nntide),-1,status)
           call ncdf_read_global('obceast_tide.nc','tidevbtmag',tidevbtemag(jpjed:jpjef,1:nntide),-1,status)
           call ncdf_read_global('obceast_tide.nc','tidevbtpha',tidevbtepha(jpjed:jpjef,1:nntide),-1,status)

         ENDIF

         IF( lp_obc_west) THEN
            ! initialisation to zero
            sshwdta(:,:) = 0.e0
            ubtwdta(:,:) = 0.e0
            !                                        ! ================== !
            IF( nobc_dta == 0 )   THEN               ! initial state used !
               !                                     ! ================== !
               !  Fills swdta, twdta, uwdta (global arrays)
               !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
               DO ji = niw0, niw1
                  DO jj = njw0p1, njw1m1
                     ij = jj -1 + njmpp
                     sshwdta(ij,1) = sshn(ji,jj) * tmask(ji,jj,1)
                  END DO
               END DO
            ENDIF
!ylu
           tidesshwmag(:,:)=0.0
           tidesshwpha(:,:)=0.0
           tidevbtwmag(:,:)=0.0
           tidevbtwpha(:,:)=0.0
!!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
           call ncdf_read_global('obcwest_tide.nc','tidesshmag',tidesshwmag(jpjwd:jpjwf,1:nntide),-1,status)
           call ncdf_read_global('obcwest_tide.nc','tidesshpha',tidesshwpha(jpjwd:jpjwf,1:nntide),-1,status)
           call ncdf_read_global('obcwest_tide.nc','tidevbtmag',tidevbtwmag(jpjwd:jpjwf,1:nntide),-1,status)
           call ncdf_read_global('obcwest_tide.nc','tidevbtpha',tidevbtwpha(jpjwd:jpjwf,1:nntide),-1,status)


           !CALL flioclo (fid_w)
         ENDIF

         IF( lp_obc_north) THEN
            ! initialisation to zero
            sshndta(:,:) = 0.e0
            vbtndta(:,:) = 0.e0
            !                                        ! ================== !
            IF( nobc_dta == 0 )   THEN               ! initial state used !
               !                                     ! ================== !
               !  Fills sndta, tndta, vndta (global arrays)
               !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
               DO jj = njn0, njn1
                  DO ji = nin0p1, nin1m1
                     DO jk = 1, jpkm1
                        ii = ji -1 + nimpp
                        vbtndta(ii,1) = vbtndta(ii,1) + vndta(ii,jk,1)*fse3v(ji,jj,jk)
                     END DO
                     sshndta(ii,1) = sshn(ii,jj+1) * tmask(ji,jj+1,1)
                  END DO
               END DO
            ENDIF
!ylu
           tidesshnmag(:,:)=0.0
           tidesshnpha(:,:)=0.0
           tidevbtnmag(:,:)=0.0
           tidevbtnpha(:,:)=0.0
!!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
           call ncdf_read_global('obcnorth_tide.nc','tidesshmag',tidesshnmag(jpind:jpinf,1:nntide),-1,status)
           call ncdf_read_global('obcnorth_tide.nc','tidesshpha',tidesshnpha(jpind:jpinf,1:nntide),-1,status)
           call ncdf_read_global('obcnorth_tide.nc','tidevbtmag',tidevbtnmag(jpind:jpinf,1:nntide),-1,status)
           call ncdf_read_global('obcnorth_tide.nc','tidevbtpha',tidevbtnpha(jpind:jpinf,1:nntide),-1,status)

         ENDIF

         IF( lp_obc_south) THEN
            ! initialisation to zero
            sshsdta(:,:) = 0.e0
            vbtsdta(:,:) = 0.e0
            !                                        ! ================== !
            IF( nobc_dta == 0 )   THEN               ! initial state used !
               !                                     ! ================== !
               !  Fills ssdta, tsdta, vsdta (global arrays)
               !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
               DO jj = njs0, njs1
                  DO ji = nis0p1, nis1m1
                     DO jk = 1, jpkm1
                        ii = ji -1 + nimpp
                        vbtsdta(ii,1) = vbtsdta(ii,1) + vsdta(ii,jk,1)*fse3v(ji,jj,jk)
                     END DO
                     sshsdta(ii,1) = sshn(ji,jj) * tmask(ii,jj,1)
                  END DO
               END DO
            ENDIF
!ylu
           tidesshsmag(:,:)=0.0
           tidesshspha(:,:)=0.0
           tidevbtsmag(:,:)=0.0
           tidevbtspha(:,:)=0.0
!!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
           call ncdf_read_global('obcsouth_tide.nc','tidesshmag',tidesshsmag(jpisd:jpisf,1:nntide),-1,status)
           call ncdf_read_global('obcsouth_tide.nc','tidesshpha',tidesshspha(jpisd:jpisf,1:nntide),-1,status)
           call ncdf_read_global('obcsouth_tide.nc','tidevbtmag',tidevbtsmag(jpisd:jpisf,1:nntide),-1,status)
           call ncdf_read_global('obcsouth_tide.nc','tidevbtpha',tidevbtspha(jpisd:jpisf,1:nntide),-1,status)

           !CALL flioclo (fid_s)
         ENDIF

      ENDIF        !       END IF kt == nit000

     
!!------------------------------------------------------------------------------------
     ! 2.      Initialize the time we are at. Does this every time the routine is called,
     !         excepted when nobc_dta = 0
     !
     IF( nobc_dta == 0) THEN
        itimo = 1
        zxy   = 0
     ELSE
        IF(itobc == 1) THEN
           itimo = 1
        ELSE IF (itobc == 12) THEN      !   BC are monthly
           ! we assume we have climatology in that case
           iman  = 12
           i15   = nday / 16
           imois = nmonth + i15 - 1
           IF( imois == 0 )   imois = iman
           itimo = imois
        ELSE
           IF(lwp) WRITE(numout,*) 'data other than constant or monthly',kt
           iman  = itobc
           itimo = FLOOR( kt*rdt / ztcobc(1))
           isrel=kt*rdt
        ENDIF
     ENDIF

     ! 2. Read two records in the file if necessary
     ! ---------------------------------------------

     IF( nobc_dta == 1 .AND. nlecto == 1 ) THEN
 !sujie add  --------------
      ! Calendar computation
         IF( itobc == 1 )   THEN            !  BC are constant in time
            ntobc1 = 1
            ntobc2 = 1
            ntobc3 = 1
         ELSE IF( itobc == 12 )   THEN      !   BC are monthly
            ntobc1 = itimo         ! first file record used
            ntobc2 = ntobc1 + 1    ! second  file record used
            ntobc3 = ntobc1 + 2    ! last  file record used
            ntobc1 = MOD( ntobc1, iman )
            IF( ntobc1 == 0 )   ntobc1 = iman
            ntobc2 = MOD( ntobc2, iman )
            IF( ntobc2 == 0 )   ntobc2 = iman
            ntobc3 = MOD( ntobc3, iman )
            IF( ntobc3 == 0 )   ntobc3 = iman

         ELSE
            isrel=kt*rdt
            ntobc1 = itimo         ! first file record used
            ntobc2 = ntobc1 + 1    ! second  file record used
            ntobc3 = ntobc1 + 2    ! last  file record used
            ntobc1 = MOD( ntobc1, iman )
            IF( ntobc1 == 0 )   ntobc1 = iman
            ntobc2 = MOD( ntobc2, iman )
            IF( ntobc2 == 0 )   ntobc2 = iman
            ntobc3 = MOD( ntobc3, iman )
            IF( ntobc3 == 0 )   ntobc3 = iman
         ENDIF
! ------------------

        IF( lp_obc_east ) THEN
           ! ... Read datafile and set sea surface height and barotropic velocity
           ! ... initialise the sshedta, ubtedta arrays
           sshedta(:,0) = sshedta(:,1)
           ubtedta(:,0) = ubtedta(:,1)
!!DB/CN: Replace IOIPSL with lib_ncdf
           CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc1,pdta_2D=sshedta(jpjed:jpjef,1),fn='obceast_TS.nc')
           CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc2,pdta_2D=sshedta(jpjed:jpjef,2),fn='obceast_TS.nc')
           IF( lk_dynspg_ts ) THEN
              CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc3,pdta_2D=sshedta(jpjed:jpjef,3),fn='obceast_TS.nc')
           ENDIF
           CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc1,pdta_2D=ubtedta(jpjed:jpjef,1),fn='obceast_U.nc')
           CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc2,pdta_2D=ubtedta(jpjed:jpjef,2),fn='obceast_U.nc')
           IF( lk_dynspg_ts ) THEN
              CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc3,pdta_2D=ubtedta(jpjed:jpjef,3),fn='obceast_U.nc')
           ENDIF

        ENDIF

        IF( lp_obc_west ) THEN
           ! ... Read datafile and set temperature, salinity and normal velocity
           ! ... initialise the swdta, twdta, uwdta arrays
           sshwdta(:,0) = sshwdta(:,1)
           ubtwdta(:,0) = ubtwdta(:,1)
!!DB/CN: Replace IOIPSL with lib_ncdf
           CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc1,pdta_2D=sshwdta(jpjwd:jpjwf,1),fn='obcwest_TS.nc')
           CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc2,pdta_2D=sshwdta(jpjwd:jpjwf,2),fn='obcwest_TS.nc')
           IF( lk_dynspg_ts ) THEN
              CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc3,pdta_2D=sshwdta(jpjwd:jpjwf,3),fn='obcwest_TS.nc')
           ENDIF
           CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc1,pdta_2D=ubtwdta(jpjwd:jpjwf,1),fn='obcwest_U.nc')
           CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc2,pdta_2D=ubtwdta(jpjwd:jpjwf,2),fn='obcwest_U.nc')
           IF( lk_dynspg_ts ) THEN
              CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc3,pdta_2D=ubtwdta(jpjwd:jpjwf,3),fn='obcwest_U.nc')
           ENDIF

        ENDIF

        IF( lp_obc_north) THEN
           ! ... Read datafile and set sea surface height and barotropic velocity
           ! ... initialise the sshndta, ubtndta arrays
           sshndta(:,0) = sshndta(:,1)
           vbtndta(:,0) = vbtndta(:,1)
!!DB/CN: Replace IOIPSL with lib_ncdf
           CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc1,pdta_2D=sshndta(jpind:jpinf,1),fn='obcnorth_TS.nc')
           CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc2,pdta_2D=sshndta(jpind:jpinf,2),fn='obcnorth_TS.nc')
           IF( lk_dynspg_ts ) THEN
               CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc3,pdta_2D=sshndta(jpind:jpinf,3),fn='obcnorth_TS.nc')
           ENDIF
           CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc1,pdta_2D=vbtndta(jpind:jpinf,1),fn='obcnorth_V.nc')
           CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc2,pdta_2D=vbtndta(jpind:jpinf,2),fn='obcnorth_V.nc')
           IF( lk_dynspg_ts ) THEN
              CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc3,pdta_2D=vbtndta(jpind:jpinf,3),fn='obcnorth_V.nc')
           ENDIF

        ENDIF

        IF( lp_obc_south) THEN
           ! ... Read datafile and set sea surface height and barotropic velocity
           ! ... initialise the sshsdta, ubtsdta arrays
           sshsdta(:,0) = sshsdta(:,1)
           vbtsdta(:,0) = vbtsdta(:,1)
!!DB/CN: Replace IOIPSL with lib_ncdf
           CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc1,pdta_2D=sshsdta(jpisd:jpisf,1),fn='obcsouth_TS.nc')
           CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc2,pdta_2D=sshsdta(jpisd:jpisf,2),fn='obcsouth_TS.nc')
           IF( lk_dynspg_ts ) THEN
               CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc3,pdta_2D=sshsdta(jpisd:jpisf,3),fn='obcsouth_TS.nc')
           ENDIF
           CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc1,pdta_2D=vbtsdta(jpisd:jpisf,1),fn='obcsouth_V.nc')
           CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc2,pdta_2D=vbtsdta(jpisd:jpisf,2),fn='obcsouth_V.nc')
           IF( lk_dynspg_ts ) THEN
               CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc3,pdta_2D=vbtsdta(jpisd:jpisf,3),fn='obcsouth_V.nc')
           ENDIF

        ENDIF

      ENDIF        !      end of the test on the condition to read or not the files

     ! 3.  Call at every time step : Linear interpolation of BCs to current time step
     ! ----------------------------------------------------------------------

      IF( lk_dynspg_ts ) THEN
         isrel = (kt-1)*rdt + kbt*rdtbt

         IF( nobc_dta == 1 ) THEN
            isrel = (kt-1)*rdt + kbt*rdtbt
            itimo  = FLOOR(  kt*rdt    / (ztcobc(2)-ztcobc(1)) )
            itimom = FLOOR( (kt-1)*rdt / (ztcobc(2)-ztcobc(1)) )
            itimop = FLOOR( (kt+1)*rdt / (ztcobc(2)-ztcobc(1)) )
            IF( itimom == itimo .AND. itimop == itimo ) THEN
               ntobcm = ntobc1
               ntobcp = ntobc2

            ELSEIF ( itimom <= itimo .AND. itimop == itimo ) THEN
               IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimo ) THEN
                  ntobcm = ntobc1-1
                  ntobcp = ntobc2-1
               ELSE
                  ntobcm = ntobc1
                  ntobcp = ntobc2
               ENDIF

            ELSEIF ( itimom == itimo .AND. itimop >= itimo ) THEN
               IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimop ) THEN
                  ntobcm = ntobc1
                  ntobcp = ntobc2
               ELSE
                  ntobcm = ntobc1+1
                  ntobcp = ntobc2+1
               ENDIF

            ELSEIF ( itimom == itimo-1 .AND. itimop == itimo+1 ) THEN
               IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimo ) THEN
                  ntobcm = ntobc1-1
                  ntobcp = ntobc2-1
               ELSEIF (  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimop ) THEN
                  ntobcm = ntobc1
                  ntobcp = ntobc2
               ELSEIF (  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) == itimop ) THEN
                  ntobcm = ntobc1+1
                  ntobcp = ntobc2+2
               ELSE
                  IF(lwp) WRITE(numout, *) 'obc_dta_bt: You should not have seen this print! error 1?'
               ENDIF
            ELSE
               IF(lwp) WRITE(numout, *) 'obc_dta_bt: You should not have seen this print! error 2?'
            ENDIF

         ENDIF

      ELSE IF( lk_dynspg_exp ) THEN
         isrel=kt*rdt
         ntobcm = ntobc1
         ntobcp = ntobc2
      ENDIF

      IF( itobc == 1 .OR. nobc_dta == 0 ) THEN
         zxy = 0.e0
      ELSE IF( itobc == 12 ) THEN
         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
      ELSE
         zxy = (ztcobc(ntobcm)-FLOAT(isrel)) / (ztcobc(ntobcm)-ztcobc(ntobcp))
      ENDIF


!DB: 2007.12.27 -- Note that I now keep the mean components separate
!    Formulae should work for any number of tide components
!    REM: velocity components already zxy-weighted & multiplied by ramp
      IF( lp_obc_east ) THEN    

         do ji = nie0, nie1
            do jj = nje0p1, nje1m1
!!DB: REM that ufoe has been multiplied by ramp
               ubtfoe0(jj) = 0.0
               do jk = 1, jpkm1
                  ubtfoe0(jj) = ubtfoe0(jj) + ufoe(jj,jk) * e3t(jk)
               enddo
               ij = jj -1 + njmpp
!!DB 
               sshfoe0(jj) = eta_e(ij) * temsk(jj,1)
               ubtfoe0(jj) = ubtfoe0(jj) * uemsk(jj,1) 
               ubtfoe(jj) = 0.0
               sshfoe(jj) = 0.0

#ifdef key_ADJ_TRANSPORT
!!DB: Keep this ifdef open

#endif
               do ntide = 1, nntide
                  sshfoe(jj) = sshfoe(jj) + tidesshemag(ij,ntide)  &
                       * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
                       - tidesshepha(ij,ntide) ) * temsk(jj,1)

                  ubtfoe(jj) = ubtfoe(jj) &
                       + tidevbtemag(ij,ntide)  &
                       * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
                       - tidevbtepha(ij,ntide) ) * uemsk(jj,1)
               enddo
!!DB 
               sshfoe(jj) = ramp*(sshfoe0(jj) + sshfoe(jj))
               ubtfoe(jj) = ramp*ubtfoe(jj) + ubtfoe0(jj)

            enddo
         enddo

!-----------------------------------------------

      ENDIF

!DB: 2007.12.27 -- Note that I now keep the mean components separate
!    Formulae should work for any number of tide components
!    REM: velocity components already zxy-weighted & multiplied by ramp
      IF( lp_obc_west) THEN   

         do ji = niw0, niw1
            do jj = njw0p1, njw1m1
               ubtfow0(jj) = 0.0 
               do jk = 1, jpkm1
                  ubtfow0(jj) = ubtfow0(jj) + ufow(jj,jk) * e3t(jk)
               enddo
               ij = jj -1 + njmpp

               ubtfow0(jj) = ubtfow0(jj) * uwmsk(jj,1) 
               ubtfow(jj) = 0.0
               sshfow0(jj) = eta_w(ij)*twmsk(jj,1)
               sshfow(jj) = 0.0

               do ntide = 1, nntide
                  sshfow(jj) =  sshfow(jj) + tidesshwmag(ij,ntide)  &
                       * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
                       - tidesshwpha(ij,ntide) ) * twmsk(jj,1) 
                  
                  ubtfow(jj) = ubtfow(jj) &
                       + tidevbtwmag(ij,ntide)  &
                       * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
                       - tidevbtwpha(ij,ntide) ) * uwmsk(jj,1)
               enddo
               sshfow(jj) =  ramp*(sshfow(jj) + sshfow0(jj))
               ubtfow(jj) = ramp*ubtfow(jj)+ ubtfow0(jj) 
            enddo
         enddo

!------------------------------------------------
      ENDIF

!---------------------------------------------------------------------------
!byoung for St. Lawrence River runoff: DELETED -- see old code

!----------------------------------------------------------------------------

!!DB: north should be closed. If not, modify this section
      IF( lp_obc_north) THEN           !  fills sshfon, vbtfon (local to each processor)
         DO ji = nin0p1, nin1m1
            ii = ji -1 + nimpp
            sshfon(ji) = ( zxy * sshndta(ii,2) + (1.-zxy) * sshndta(ii,1) ) * tnmsk(ji,1)
            vbtfon(ji) = ( zxy * vbtndta(ii,2) + (1.-zxy) * vbtndta(ii,1) ) * vnmsk(ji,1)
!ylu
          do ntide = 1,4  
            sshfon(ji) = sshfon(ji) + ramp*tidesshnmag(ii,ntide)  &
                         * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
                         - tidesshnpha(ii,ntide) ) * tnmsk(ji,1) 
            vbtfon(ji) = vbtfon(ji) + ramp*tidevbtnmag(ii,ntide)  &
                         * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
                         - tidevbtnpha(ii,ntide) ) * vnmsk(ji,1) 
          enddo
         END DO
      ENDIF


!DB: 2007.12.27 -- Note that I now keep the mean components separate
!    Formulae should work for any number of tide components
!    REM: velocity components already zxy-weighted & multiplied by ramp
      IF( lp_obc_south) THEN     

         do jj = njs0, njs1
            do ji = nis0p1, nis1m1
               vbtfos0(ji) = 0.0 
               do jk = 1, jpkm1
                  vbtfos0(ji) = vbtfos0(ji) + vfos(ji,jk) * e3t(jk)
               enddo
               ii = ji -1 + nimpp

               vbtfos0(ji) = vbtfos0(ji)*vsmsk(ji,1)
               vbtfos(ji) = 0.0
               sshfos0(ji) = eta_s(ii)* tsmsk(ji,1)
               sshfos(ji) = 0.0

               do ntide = 1, nntide
                  sshfos(ji) = sshfos(ji) + tidesshsmag(ii,ntide)  &
                       * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
                       - tidesshspha(ii,ntide) ) * tsmsk(ji,1)

                  vbtfos(ji) =  vbtfos(ji) + & 
                       tidevbtsmag(ii,ntide)  &
                       * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
                       - tidevbtspha(ii,ntide) ) * vsmsk(ji,1) 
               enddo
               sshfos(ji) = ramp*(sshfos(ji) + sshfos0(ji))
               vbtfos(ji) =  ramp*vbtfos(ji) +  vbtfos0(ji)
            enddo
         enddo

      ENDIF

  END SUBROUTINE obc_dta_bt


#else
  !!-----------------------------------------------------------------------------
  !!   Default option
  !!-----------------------------------------------------------------------------
  SUBROUTINE obc_dta_bt ( kt, kbt )       ! Empty routine
     !! * Arguments
     INTEGER,INTENT(in) :: kt
     INTEGER, INTENT( in ) ::   kbt         ! barotropic ocean time-step index
     if(lwp) WRITE(numout,*) 'obc_dta_bt: You should not have seen this print! error?', kt
     if(lwp) WRITE(numout,*) 'obc_dta_bt: You should not have seen this print! error?', kbt
  END SUBROUTINE obc_dta_bt
#endif

!!DB/CN 2008.11
  SUBROUTINE obc_dta_gv (cldim,clobc,kobcij,ktobc,pdta_2D,pdta_3D,pdta_4D,fn)
     !!-----------------------------------------------------------------------------
     !!                       ***  SUBROUTINE obc_dta_gv  ***
     !!
     !! ** Purpose :   Read an OBC forcing field from netcdf file
     !!                Input file are supposed to be 3D e.g.
     !!                - for a South or North OB : longitude x depth x time
     !!		- for a West or East OB : latitude x depth x time
     !!
     !! History :
     !!        !  04-06 (A.-M. Treguier, F. Durand) Original code
     !!        !  05-02 (J. Bellier, C. Talandier) use fliocom CALL
     !!----------------------------------------------------------------------------
     !! * Arguments
     INTEGER, INTENT(IN) ::   &
!        ifid  ,               & ! netcdf file name identifier
        kobcij,               & ! Horizontal (i or j) dimension of the array
        ktobc                   ! starting time index read
     CHARACTER(LEN=*), INTENT(IN)    ::   &
        cldim,                & ! dimension along which is the open boundary ('x' or 'y')
        clobc                   ! name of the netcdf variable read
     REAL, DIMENSION(kobcij,jpk,1), INTENT(OUT), OPTIONAL ::   &
        pdta_3D                 ! 3D array of OBC forcing field

     REAL, DIMENSION(kobcij,1), INTENT(OUT), OPTIONAL ::   &
        pdta_2D                 ! 2D array of OBC forcing field
     REAL, DIMENSION(kobcij,jpk,5,1), INTENT(OUT), OPTIONAL ::   &
        pdta_4D                 ! 4D array of OBC forcing field
!!CN
     CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fn ! name of netCDF file
     LOGICAL :: fnswitch ! logical indicating whether to use filename & lib_ncdf or not
                          ! to be removed once IOIPSL is entirely gone from obcdta
     REAL, DIMENSION(kobcij) :: buf_1D             ! Needed for lib_ncdf calls
     REAL, DIMENSION(kobcij,jpk) :: buf_2D         ! Needed for lib_ncdf calls
     REAL, DIMENSION(kobcij,jpk,5) :: buf_3D       ! Needed for lib_ncdf calls
     INTEGER :: f_stat
     LOGICAL :: tide

     !! * Local declarations
     INTEGER ::   indim
     LOGICAL ::   l_exv
     INTEGER,DIMENSION(4) ::   f_d, istart, icount
     REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   v_tmp_4
     !----------------------------------------------------------------------

     !!CN: Replacing IOIPSL calls with lib_ncdf
     fnswitch = .TRUE.
     CALL ncdf_get_num_dims(fn, TRIM(clobc), indim, f_stat);

     l_exv = .TRUE.
     IF( l_exv )   THEN
        ! checks the number of dimensions
        IF( indim == 2 )   THEN
           istart(1:2) = (/ 1     , ktobc /)
           icount(1:2) = (/ kobcij, 1     /)
           CALL ncdf_read(fn, TRIM(clobc), buf_1D, f_stat)
           pdta_2D(:,1) = buf_1D
        ELSE IF( indim == 3 )   THEN
           istart(1:3) = (/ 1     , 1    , ktobc /)
           icount(1:3) = (/ kobcij, jpk  , 1     /)
           CALL ncdf_read_global(fn, TRIM(clobc), buf_2D, -ktobc, f_stat)
           pdta_3D(:,:,1) = buf_2D
        ELSE IF( indim == 4 )   THEN
           istart(1:4) = (/ 1, 1, 1, ktobc /)
           icount(1:4) = (/ kobcij, jpk  , 5, 1 /)
           CALL ncdf_read_global(fn, TRIM(clobc), buf_3D, -ktobc, f_stat)
           pdta_4D(:,:,:,1) = buf_3d
        ELSE
           IF( lwp )   THEN
              WRITE(numout,*) ' Problem in OBC file for ',TRIM(clobc),' :'
              WRITE(numout,*) ' number of dimensions (not 3 or 4) =',indim
           ENDIF
           STOP
        ENDIF
     ELSE
        if(lwp) WRITE(numout,*) ' Variable ',TRIM(clobc),' not found'
     ENDIF

  END SUBROUTINE obc_dta_gv

#else
  !!--------------------------------------------------------------------
  !!  default option  :  Dummy module    NO Open Boundary Conditions
  !!--------------------------------------------------------------------
CONTAINS
  SUBROUTINE obc_dta( kt )             ! Dummy routine
    INTEGER, INTENT (in) :: kt
    if(lwp) WRITE(numout,*) 'obc_dta: You should not have seen this print! error?', kt
  END SUBROUTINE obc_dta
  SUBROUTINE obc_dta_bt( kt, jn)             ! Dummy routine
    INTEGER, INTENT (in) :: kt, jn
    if(lwp) WRITE(numout,*) 'obc_dta_bt: You should not have seen this print! error?', kt
    if(lwp) WRITE(numout,*) 'obc_dta_bt: You should not have seen this print! error?', jn
  END SUBROUTINE obc_dta_bt
#endif


  !!=====================================================================
END MODULE obcdta
