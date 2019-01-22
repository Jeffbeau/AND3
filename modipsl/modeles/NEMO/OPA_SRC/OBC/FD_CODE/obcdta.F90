MODULE obcdta
   !!==============================================================================
   !!                            ***  MODULE obcdta  ***
   !! Open boundary data : read the data for the open boundaries.
   !!==============================================================================
   !! History :
!! 20071221, FD, add correction by FD which reads the barotropic parameters within the
!!    same baroclinic routine.
   !!        !  98-05 (J.M. Molines) Original code
   !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
   !!   9.0  !  04-06 (F. Durand, A-M. Treguier) Netcdf BC files on input
!! FD add tangential velocity 03/08
   !!------------------------------------------------------------------------------
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
! FD: needed for AGRIF
   USE obc_par         ! ocean open boundary conditions
   USE obc_oce         ! ocean open boundary conditions
   USE daymod          ! calendar
   USE in_out_manager  ! I/O logical units
   USE lib_mpp         ! distributed memory computing
   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
   USE iom
#  if defined key_dynspg_rl
   USE obccli          ! climatological obc, use only in rigid-lid case
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
      ntobc        ! number of time steps in OBC files 

   REAL(wp), DIMENSION(:), ALLOCATABLE :: tcobc      ! time_counter variable of BCs

   !! * Substitutions
#  include "domzgr_substitute.h90"                                            
#  include "obc_vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !!   OPA 9.0 , LODYC-IPSL  (2003)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdta.F90,v 1.13 2007/02/26 17:26:05 opalod Exp $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE obc_dta( kt )
      !!--------------------------------------------------------------------
      !!              ***  SUBROUTINE obc_dta  ***
      !!                   
      !! ** Purpose :   Find the climatological boundary arrays for the specified date,
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
      INTEGER, DIMENSION(1) ::  itobce, itobcw,  & ! number of time steps in OBC files
                                itobcs, itobcn     !    "       "       "       "
      INTEGER ::  istop        
      INTEGER ::  iprint                              ! frequency for printouts.
      INTEGER ::  idvar, id_e, id_w, id_n, id_s       ! file identifiers
      LOGICAL :: llnot_done
      CHARACTER(LEN=25) :: cl_vname
      !!--------------------------------------------------------------------

      IF( lk_dynspg_rl )  THEN
         CALL obc_dta_psi (kt)     ! update bsf data at open boundaries
         IF ( nobc_dta == 1 .AND. kt == nit000 ) CALL ctl_stop( 'obcdta : time-variable psi boundary data not allowed yet' )
      ENDIF
           
      ! 1.   First call: check time frames available in files.
      ! -------------------------------------------------------

      IF ( kt == nit000 ) THEN
      
         nlecto = 0
         itobce(1) = 0   ;    itobcw(1) = 0
         itobcn(1) = 0   ;    itobcs(1) = 0

         IF (lwp) WRITE(numout,*)
         IF (lwp) WRITE(numout,*)     'obc_dta : find boundary data'
         IF (lwp) WRITE(numout,*)     '~~~~~~~'
             
         IF ( nobc_dta == 0 ) THEN
            IF(lwp) WRITE(numout,*)  '  OBC data taken from initial conditions.'
            ntobc1 = 1
            ntobc2 = 1
            ntobc  = 0
         ELSE    
            IF (lwp) WRITE(numout,*)  '  OBC data taken from netcdf files.'
            IF (lwp) WRITE(numout,*)  '  climatology (T/F):',ln_obc_clim
            ! check the number of time steps in the files.
            cl_vname = 'time_counter'
            IF ( lp_obc_east ) THEN
               CALL iom_open ( 'obceast_TS.nc' , id_e )
               idvar = iom_varid( id_e, cl_vname, kdimsz = itobce )
            ENDIF
            IF ( lp_obc_west ) THEN
               CALL iom_open ( 'obcwest_TS.nc' , id_w )
               idvar = iom_varid( id_w, cl_vname, kdimsz = itobcw )
            ENDIF
            IF ( lp_obc_north ) THEN
               CALL iom_open ( 'obcnorth_TS.nc', id_n )
               idvar = iom_varid( id_n, cl_vname, kdimsz = itobcn )
            ENDIF
            IF ( lp_obc_south ) THEN
               CALL iom_open ( 'obcsouth_TS.nc', id_s )
               idvar = iom_varid( id_s, cl_vname, kdimsz = itobcs )
            ENDIF

! DW
            ntobc = MAX(itobce(1),itobcw(1),itobcn(1),itobcs(1))
            istop = 0
            IF ( lp_obc_east  .AND. itobce(1) /= ntobc ) istop = 1 
            IF ( lp_obc_west  .AND. itobcw(1) /= ntobc ) istop = 1      
            IF ( lp_obc_north .AND. itobcn(1) /= ntobc ) istop = 1
            IF ( lp_obc_south .AND. itobcs(1) /= ntobc ) istop = 1 
            IF ( istop /= 0 )  THEN
               WRITE(ctmp1,*) ' east, west, north, south: ', itobce(1), itobcw(1), itobcn(1), itobcs(1)
               CALL ctl_stop( 'obcdta : all files must have the same number of time steps', ctmp1 )
            ENDIF
            IF ( ntobc == 1 ) THEN
               IF ( lwp ) WRITE(numout,*) ' obcdta found one time step only in the OBC files'
            ELSE
               ALLOCATE (tcobc(ntobc))
               llnot_done = .TRUE.
               IF ( lp_obc_east ) THEN
                  IF ( llnot_done ) THEN
                     CALL iom_gettime ( id_e, TRIM(cl_vname), tcobc )
                     llnot_done = .FALSE.
                  ENDIF
                  CALL iom_close (id_e)
               ENDIF
               IF ( lp_obc_west ) THEN
                  IF ( llnot_done ) THEN
                     CALL iom_gettime ( id_w, TRIM(cl_vname), tcobc )
                     llnot_done = .FALSE.
                 ENDIF
                 CALL iom_close (id_w)
               ENDIF
               IF ( lp_obc_north ) THEN
                  IF ( llnot_done ) THEN
                     CALL iom_gettime ( id_n, TRIM(cl_vname), tcobc )
                     llnot_done = .FALSE.
                  ENDIF
                  CALL iom_close (id_n)
               ENDIF
               IF ( lp_obc_south ) THEN
                  IF ( llnot_done ) THEN
                     CALL iom_gettime ( id_s, TRIM(cl_vname), tcobc )
                     llnot_done = .FALSE.
                  ENDIF
                  CALL iom_close (id_s)
               ENDIF
               IF ( lwp ) WRITE(numout,*) ' obcdta found', ntobc,' time steps in the OBC files'
               IF ( .NOT. ln_obc_clim .AND. ntobc == 12 ) THEN
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
      
! FD: barotropic init moved hereafter
         IF( lp_obc_east  ) vbtfoe(:) = 0.e0
         IF( lp_obc_west  ) vbtfow(:) = 0.e0
         IF( lp_obc_south ) ubtfos(:) = 0.e0
         IF( lp_obc_north ) ubtfon(:) = 0.e0

       ! 1.2  Data temperature, salinity, normal velocities set to zero
       !                        or initial conditions if nobc_dta == 0
       ! --------------------------------------------------------------

         IF( lp_obc_east )   THEN
            ! initialisation to zero
            sedta(:,:,:) = 0.e0
            tedta(:,:,:) = 0.e0
            uedta(:,:,:) = 0.e0
! FD: add tangential velocity
            vedta(:,:,:) = 0.e0
            !                                    ! ================== !
            IF( nobc_dta == 0 )   THEN           ! initial state used
               !                                 ! ================== !
               !  Fills sedta, tedta, uedta (global arrays)
               !  Remark: this works for njzoom = 1.
               !          Should the definition of ij include njzoom?
               DO ji = nie0, nie1
                  DO jk = 1, jpkm1
! FD: increase the range
! FD                     DO jj = nje0p1, nje1m1
                     DO jj = nje0, nje1
                        ij = jj -1 + njmpp
! FD                    sedta(ij,jk,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
! FD                    tedta(ij,jk,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                        sedta(ij,jk,1) = sn(ji+1,jj,jk)*tmask(ji+1,jj,jk)
                        tedta(ij,jk,1) = tn(ji+1,jj,jk)*tmask(ji+1,jj,jk)
                        uedta(ij,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
! FD: add tangential velocity
                        vedta(ij,jk,1) = vn(ji+1,jj,jk)*vmask(ji+1,jj,jk)
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( lp_obc_west )   THEN
            ! initialisation to zero
            swdta(:,:,:) = 0.e0
            twdta(:,:,:) = 0.e0
            uwdta(:,:,:) = 0.e0
! FD: add tangential velocity
            vwdta(:,:,:) = 0.e0
            !                                    ! ================== !
            IF( nobc_dta == 0 )   THEN           ! initial state used !
               !                                 ! ================== !
               !  Fills swdta, twdta, uwdta (global arrays)
               !  Remark: this works for njzoom = 1. 
               !          Should the definition of ij include njzoom?
               DO ji = niw0, niw1
                  DO jk = 1, jpkm1
! FD: increase the range
! FD                     DO jj = njw0p1, njw1m1
                     DO jj = njw0, njw1
                        ij = jj -1 + njmpp
                        swdta(ij,jk,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                        twdta(ij,jk,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                        uwdta(ij,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
! FD: add tangential velocity
                        vwdta(ij,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( lp_obc_north)   THEN
            ! initialisation to zero
            sndta(:,:,:) = 0.e0
            tndta(:,:,:) = 0.e0
! FD: add tangential velocity
            undta(:,:,:) = 0.e0
            vndta(:,:,:) = 0.e0
            !                                    ! ================== !
            IF( nobc_dta == 0 )   THEN           ! initial state used
               !                                 ! ================== !
               !  Fills sndta, tndta, vndta (global arrays)
               !  Remark: this works for njzoom = 1. 
               !          Should the definition of ij include njzoom?
               DO jj = njn0, njn1
                  DO jk = 1, jpkm1
! FD: increase the range
! FD                     DO ji = nin0p1, nin1m1
                     DO ji = nin0, nin1
                        ii = ji -1 + nimpp
! FD                    sndta(ii,jk,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
! FD                    tndta(ii,jk,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
                        sndta(ii,jk,1) = sn(ji,jj+1,jk)*tmask(ji,jj+1,jk)
                        tndta(ii,jk,1) = tn(ji,jj+1,jk)*tmask(ji,jj+1,jk)
! FD: add tangential velocity
                        undta(ii,jk,1) = un(ji,jj+1,jk)*umask(ji,jj+1,jk)
                        vndta(ii,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF

         IF( lp_obc_south )   THEN
            ! initialisation to zero
            ssdta(:,:,:) = 0.e0
            tsdta(:,:,:) = 0.e0
! FD: add tangential velocity
            usdta(:,:,:) = 0.e0
            vsdta(:,:,:) = 0.e0
            !                                    ! ================== !
            IF( nobc_dta == 0 )   THEN           ! initial state used
               !                                 ! ================== !
               !  Fills ssdta, tsdta, vsdta (global arrays)
               !  Remark: this works for njzoom = 1. 
               !          Should the definition of ij include njzoom?
               DO jj = njs0, njs1
                  DO jk = 1, jpkm1
! FD: increase the range
! FD                     DO ji = nis0p1, nis1m1
                     DO ji = nis0, nis1
                        ii = ji -1 + nimpp
                        ssdta(ii,jk,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
                        tsdta(ii,jk,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
! FD: add tangential velocity
                        usdta(ii,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
                        vsdta(ii,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
                     END DO
                  END DO
               END DO
            ENDIF
         ENDIF

! FD
! FD the barotropic init is done hereafter in obc_dta
! FD
# if defined key_dynspg_ts || defined key_dynspg_exp
! FD
          IF( lp_obc_east ) THEN
             ! initialisation to zero
             sshedta(:,:) = 0.e0
             !                                        ! ================== !
             IF( nobc_dta == 0 )   THEN               ! initial state used !
                !                                     ! ================== !
                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
                DO ji = nie0, nie1
                   DO jj = nje0, nje1
                      ij = jj -1 + njmpp
                      sshedta(ij,:) = sshn(ji+1,jj) * tmask(ji+1,jj,1)
                   END DO
                END DO
             ENDIF
          ENDIF

          IF( lp_obc_west) THEN
             ! initialisation to zero
             sshwdta(:,:) = 0.e0
             !                                        ! ================== !
             IF( nobc_dta == 0 )   THEN               ! initial state used !
                !                                     ! ================== !
                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
                DO ji = niw0, niw1
                   DO jj = njw0, njw1
                      ij = jj -1 + njmpp
                      sshwdta(ij,:) = sshn(ji,jj) * tmask(ji,jj,1)
                   END DO
                END DO
             ENDIF
          ENDIF

          IF( lp_obc_north) THEN
             ! initialisation to zero
             sshndta(:,:) = 0.e0
             !                                        ! ================== !
             IF( nobc_dta == 0 )   THEN               ! initial state used !
                !                                     ! ================== !
                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
                DO jj = njn0, njn1
                   DO ji = nin0, nin1
                      ii = ji -1 + nimpp
                      sshndta(ii,:) = sshn(ji,jj+1) * tmask(ji,jj+1,1)
                   END DO
                END DO
             ENDIF
          ENDIF

          IF( lp_obc_south) THEN
             ! initialisation to zero
             sshsdta(:,:) = 0.e0
             !                                        ! ================== !
             IF( nobc_dta == 0 )   THEN               ! initial state used !
                !                                     ! ================== !
                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
                DO jj = njs0, njs1
                   DO ji = nis0, nis1
                      ii = ji -1 + nimpp
                      sshsdta(ii,:) = sshn(ji,jj) * tmask(ji,jj,1)
                   END DO
                END DO
             ENDIF
          ENDIF


         ! 1.2  Compute depth averaged
         ! --------------------------------------------------------------------

          IF( lp_obc_east ) THEN
             ! initialisation to zero
             ubtedta(:,:) = 0.e0
! FD: add tangential velocity
             vbtedta(:,:) = 0.e0
                DO ji = nie0, nie1
                   DO jj = nje0, nje1
                      ij = jj -1 + njmpp
                      DO jk = 1, jpkm1
                         ubtedta(ij,1) = ubtedta(ij,1) + uedta(ij,jk,1)*fse3u(ji,jj,jk)
! FD: add tangential velocity
                         vbtedta(ij,1) = vbtedta(ij,1) + vedta(ij,jk,1)*fse3v(ji+1,jj,jk)
                      END DO
                      ubtedta(ij,1) = ubtedta(ij,1) * umask(ji,jj,1)
! FD: add tangential velocity
                      vbtedta(ij,1) = vbtedta(ij,1) * vmask(ji+1,jj,1)
                   END DO
                END DO
          ENDIF

          IF( lp_obc_west) THEN
             ! initialisation to zero
             ubtwdta(:,:) = 0.e0
! FD: add tangential velocity
             vbtwdta(:,:) = 0.e0
                DO ji = niw0, niw1
                   DO jj = njw0, njw1
                      ij = jj -1 + njmpp
                      DO jk = 1, jpkm1
                         ubtwdta(ij,1) = ubtwdta(ij,1) + uwdta(ij,jk,1)*fse3u(ji,jj,jk)
! FD: add tangential velocity
                         vbtwdta(ij,1) = vbtwdta(ij,1) + vwdta(ij,jk,1)*fse3v(ji,jj,jk)
                      END DO
                      ubtwdta(ij,1) = ubtwdta(ij,1) * umask(ji,jj,1)
! FD: add tangential velocity
                      vbtwdta(ij,1) = vbtwdta(ij,1) * vmask(ji,jj,1)
                   END DO
                END DO
          ENDIF

          IF( lp_obc_north) THEN
             ! initialisation to zero
! FD: add tangential velocity
             ubtndta(:,:) = 0.e0
             vbtndta(:,:) = 0.e0
                DO jj = njn0, njn1
                   DO ji = nin0, nin1
                      ii = ji -1 + nimpp
                      DO jk = 1, jpkm1
! FD: add tangential velocity
                         ubtndta(ii,1) = ubtndta(ii,1) + undta(ii,jk,1)*fse3u(ji,jj+1,jk)
                         vbtndta(ii,1) = vbtndta(ii,1) + vndta(ii,jk,1)*fse3v(ji,jj,jk)
                      END DO
! FD: add tangential velocity
                      ubtndta(ii,1) = ubtndta(ii,1) * umask(ji,jj+1,1)
                      vbtndta(ii,1) = vbtndta(ii,1) * vmask(ji,jj,1)
                   END DO
                END DO
          ENDIF

          IF( lp_obc_south) THEN
             ! initialisation to zero
! FD: add tangential velocity
             ubtsdta(:,:) = 0.e0
             vbtsdta(:,:) = 0.e0
                DO jj = njs0, njs1
                   DO ji = nis0, nis1
                      ii = ji -1 + nimpp
                      DO jk = 1, jpkm1
! FD: add tangential velocity
                         ubtsdta(ii,1) = ubtsdta(ii,1) + usdta(ii,jk,1)*fse3u(ji,jj,jk)
                         vbtsdta(ii,1) = vbtsdta(ii,1) + vsdta(ii,jk,1)*fse3v(ji,jj,jk)
                      END DO
! FD: add tangential velocity
                      ubtsdta(ii,1) = ubtsdta(ii,1) * umask(ji,jj,1)
                      vbtsdta(ii,1) = vbtsdta(ii,1) * vmask(ji,jj,1)
                   END DO
                END DO
          ENDIF
#endif

      ENDIF        !       end if kt == nit000
      
      ! 2.  Initialize the time we are at.
      !     Does this every time the routine is called,
      !     excepted when nobc_dta = 0
      !---------------------------------------------------------------------
      IF( nobc_dta == 0 )   THEN
         itimo = 1
         zxy   = 0
      ELSE
         IF( ntobc == 1 )   THEN
            itimo = 1
         ELSE IF( ntobc == 12 )   THEN      !   BC are monthly   
            ! we assume we have climatology in that case
            iman  = 12
            i15   = nday / 16
            imois = nmonth + i15 - 1
            IF( imois == 0 )   imois = iman
            itimo = imois   
         ELSE
            IF(lwp) WRITE(numout,*) 'data other than constant or monthly', kt
            iman  = ntobc
! FD            itimo = FLOOR( kt*rdt / (tcobc(2)-tcobc(1)) )+ 1		!! WAR
            itimo = MOD( FLOOR( kt*rdt / (tcobc(2)-tcobc(1)) ), iman ) + 1		!! WAR
            isrel = kt*rdt
         ENDIF
      ENDIF
      
      ! 2.1 Read two records in the file if necessary
      ! ---------------------------------------------
      IF( ( nobc_dta == 1 ) .AND. ( ( kt == nit000 .AND. nlecto == 0 ) .OR. itimo  /= ntobc1 ) )   THEN
         nlecto = 1
      
         ! Calendar computation
         IF( ntobc == 1 )   THEN            !  BC are constant in time
            ntobc1 = 1
            ntobc2 = 1  
         ELSE IF( ntobc == 12 )   THEN      !   BC are monthly
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
! FD            ntobc2 = ntobc1 + 1    ! last  file record used
            ntobc2 = MOD( ntobc1, iman ) + 1    ! last  file record used
!!WAR            ntobc1 = MOD( ntobc1, iman )
!!WAR            IF( ntobc1 == 0 )   ntobc1 = iman
!!WAR            ntobc2 = MOD( ntobc2, iman )
!!WAR            IF( ntobc2 == 0 )   ntobc2 = iman
            IF(lwp) WRITE(numout,*) ' read obc first record file used ntobc1 ', ntobc1
            IF(lwp) WRITE(numout,*) ' read obc last  record file used ntobc2 ', ntobc2
         ENDIF
                               ! ======================= !
                               !  BCs read               !
         !                     ! ======================= !
         IF( lp_obc_east )   THEN
            ! ... Read datafile and set temperature, salinity and normal velocity
            ! ... initialise the sedta, tedta, uedta arrays
            CALL iom_open ( 'obceast_TS.nc' , id_e )
            CALL iom_get ( id_e, jpdom_unknown, 'votemper', tedta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_e, jpdom_unknown, 'votemper', tedta(:,:,2), ktime=ntobc2 )
            CALL iom_get ( id_e, jpdom_unknown, 'vosaline', sedta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_e, jpdom_unknown, 'vosaline', sedta(:,:,2), ktime=ntobc2 )
! FD move
# if defined key_dynspg_ts || defined key_dynspg_exp
	     IF( kt /= nit000 ) THEN
	          sshedta(:,0) = sshedta(:,1)
	     ENDIF
            CALL iom_get ( id_e, jpdom_unknown, 'vossurfh', sshedta(:,1), ktime=ntobc1 )
            CALL iom_get ( id_e, jpdom_unknown, 'vossurfh', sshedta(:,2), ktime=ntobc2 )
# endif
            CALL iom_close (id_e)
            !
            CALL iom_open ( 'obceast_U.nc' , id_e )
            CALL iom_get ( id_e, jpdom_unknown, 'vozocrtx', uedta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_e, jpdom_unknown, 'vozocrtx', uedta(:,:,2), ktime=ntobc2 )
            CALL iom_close ( id_e )
# if defined key_obc_tangent_vel
! FD: add tangential velocity
            CALL iom_open ( 'obceast_V.nc' , id_e )
            CALL iom_get ( id_e, jpdom_unknown, 'vomecrty', vedta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_e, jpdom_unknown, 'vomecrty', vedta(:,:,2), ktime=ntobc2 )
            CALL iom_close ( id_e )
# endif
            !  Usually printout is done only once at kt = nit000,
            !  unless nprint (namelist) > 1      
            IF( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 ) )   THEN
               WRITE(numout,*)
               WRITE(numout,*) ' Read East OBC data records ', ntobc1, ntobc2
               iprint = (jpjef-jpjed+1)/20 +1
               WRITE(numout,*) ' Temperature record 1 - printout every 3 level'
               CALL prihre( tedta(:,:,1),jpjef-jpjed+1,jpk,1,jpjef-jpjed+1,iprint, &
                  &         jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Salinity  record 1 - printout every 3 level'
               CALL prihre( sedta(:,:,1), jpjef-jpjed+1, jpk, 1, jpjef-jpjed+1, iprint, &
                  &        jpk, 1, -3, 1., numout )
! FD move
# if defined key_dynspg_ts || defined key_dynspg_exp
               WRITE(numout,*)
               WRITE(numout,*) ' Sea surface height record 1'
               CALL prihre( sshedta(:,1), jpjef-jpjed+1, 1, 1, jpjef-jpjed+1, iprint, 1, 1, -3, 1., numout )
# endif
               WRITE(numout,*)
               WRITE(numout,*) ' Normal velocity U  record 1  - printout every 3 level'
               CALL prihre( uedta(:,:,1), jpjef-jpjed+1, jpk, 1, jpjef-jpjed+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
# if defined key_obc_tangent_vel
! FD: add tangential velocity
               WRITE(numout,*)
               WRITE(numout,*) ' Tangent velocity V record 1  - printout every 3 level'
               CALL prihre( vedta(:,:,1), jpjef-jpjed+1, jpk, 1, jpjef-jpjed+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
# endif
            ENDIF
         ENDIF

         IF( lp_obc_west )   THEN
            ! ... Read datafile and set temperature, salinity and normal velocity
            ! ... initialise the swdta, twdta, uwdta arrays
            CALL iom_open ( 'obcwest_TS.nc' , id_w )
            CALL iom_get ( id_w, jpdom_unknown, 'votemper', twdta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_w, jpdom_unknown, 'votemper', twdta(:,:,2), ktime=ntobc2 )
            CALL iom_get ( id_w, jpdom_unknown, 'vosaline', swdta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_w, jpdom_unknown, 'vosaline', swdta(:,:,2), ktime=ntobc2 )
# if defined key_dynspg_ts || defined key_dynspg_exp
	     IF( kt /= nit000 ) THEN
	          sshwdta(:,0) = sshwdta(:,1)
	     ENDIF
            CALL iom_get ( id_w, jpdom_unknown, 'vossurfh', sshwdta(:,1), ktime=ntobc1 )
            CALL iom_get ( id_w, jpdom_unknown, 'vossurfh', sshwdta(:,2), ktime=ntobc2 )
# endif
            CALL iom_close (id_w)
            !
            CALL iom_open ( 'obcwest_U.nc' , id_w )
            CALL iom_get ( id_w, jpdom_unknown, 'vozocrtx', uwdta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_w, jpdom_unknown, 'vozocrtx', uwdta(:,:,2), ktime=ntobc2 )
            CALL iom_close ( id_w )
# if defined key_obc_tangent_vel
! FD: add tangential velocity
            CALL iom_open ( 'obcwest_V.nc' , id_w )
            CALL iom_get ( id_e, jpdom_unknown, 'vomecrty', vwdta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_e, jpdom_unknown, 'vomecrty', vwdta(:,:,2), ktime=ntobc2 )
            CALL iom_close ( id_w )
#endif
            !
            IF( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 ) )   THEN
               WRITE(numout,*)
               WRITE(numout,*) ' Read West OBC data records ', ntobc1, ntobc2
               iprint = (jpjwf-jpjwd+1)/20 +1
               WRITE(numout,*) ' Temperature  record 1 - printout every 3 level'
               CALL prihre( twdta(:,:,1), jpjwf-jpjwd+1, jpk, 1, jpjwf-jpjwd+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Salinity  record 1 - printout every 3 level'
               CALL prihre( swdta(:,:,1), jpjwf-jpjwd+1, jpk, 1, jpjwf-jpjwd+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
! FD move
# if defined key_dynspg_ts || defined key_dynspg_exp
               WRITE(numout,*)
               WRITE(numout,*) ' Sea surface height record 1  - printout surface level'
               CALL prihre( sshwdta(:,1), jpjwf-jpjwd+1, 1, 1, jpjwf-jpjwd+1, iprint, 1, 1, -3, 1., numout )
# endif
               WRITE(numout,*)
               WRITE(numout,*) ' Normal velocity U  record 1  - printout every 3 level'
               CALL prihre( uwdta(:,:,1), jpjwf-jpjwd+1, jpk, 1, jpjwf-jpjwd+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
# if defined key_obc_tangent_vel
! FD: add tangential velocity
               WRITE(numout,*)
               WRITE(numout,*) ' Tangent velocity V record 1  - printout every 3 level'
               CALL prihre( vwdta(:,:,1), jpjwf-jpjwd+1, jpk, 1, jpjwf-jpjwd+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
# endif
            ENDIF
         ENDIF

         IF( lp_obc_north )   THEN      
            CALL iom_open ( 'obcnorth_TS.nc', id_n )
            CALL iom_get ( id_n, jpdom_unknown, 'votemper', tndta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_n, jpdom_unknown, 'votemper', tndta(:,:,2), ktime=ntobc2 )
            CALL iom_get ( id_n, jpdom_unknown, 'vosaline', sndta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_n, jpdom_unknown, 'vosaline', sndta(:,:,2), ktime=ntobc2 )
! FD move
# if defined key_dynspg_ts || defined key_dynspg_exp
	     IF( kt /= nit000 ) THEN
	          sshndta(:,0) = sshndta(:,1)
	     ENDIF
            CALL iom_get ( id_n, jpdom_unknown, 'vossurfh', sshndta(:,1), ktime=ntobc1 )
            CALL iom_get ( id_n, jpdom_unknown, 'vossurfh', sshndta(:,2), ktime=ntobc2 )
# endif
            CALL iom_close ( id_n )                                                           
            !
# if defined key_obc_tangent_vel
! FD: add tangential velocity
            CALL iom_open ( 'obcnorth_U.nc', id_n )                                          
            CALL iom_get ( id_n, jpdom_unknown, 'vozocrtx', undta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_n, jpdom_unknown ,'vozocrtx', undta(:,:,2), ktime=ntobc2 )
            CALL iom_close ( id_n )
# endif
            CALL iom_open ( 'obcnorth_V.nc', id_n )                                          
            CALL iom_get ( id_n, jpdom_unknown, 'vomecrty', vndta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_n, jpdom_unknown ,'vomecrty', vndta(:,:,2), ktime=ntobc2 )
            CALL iom_close ( id_n )
            !
            IF( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 ) )   THEN
               WRITE(numout,*)
               WRITE(numout,*) ' Read North OBC data records ', ntobc1, ntobc2
               iprint = (jpinf-jpind+1)/20 +1
               WRITE(numout,*) ' Temperature  record 1 - printout every 3 level'
               CALL prihre( tndta(:,:,1), jpinf-jpind+1, jpk, 1, jpinf-jpind+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Salinity  record 1 - printout every 3 level'
               CALL prihre( sndta(:,:,1), jpinf-jpind+1, jpk, 1, jpinf-jpind+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
! FD move
# if defined key_dynspg_ts || defined key_dynspg_exp
               WRITE(numout,*)
               WRITE(numout,*) ' Sea surface height record 1  - printout surface level'
               CALL prihre( sshndta(:,1), jpinf-jpind+1, 1, 1, jpinf-jpind+1, iprint, 1, 1, -3, 1., numout )
# endif
# if defined key_obc_tangent_vel
! FD: add tangential velocity
               WRITE(numout,*)
               WRITE(numout,*) ' Tangent velocity U record 1  - printout every 3 level'
               CALL prihre( undta(:,:,1), jpinf-jpind+1, jpk, 1, jpinf-jpind+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
# endif
               WRITE(numout,*)
               WRITE(numout,*) ' Normal velocity V  record 1  - printout every 3 level'
               CALL prihre( vndta(:,:,1), jpinf-jpind+1, jpk, 1, jpinf-jpind+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
            ENDIF
         ENDIF

         IF( lp_obc_south )   THEN      
            CALL iom_open ( 'obcsouth_TS.nc', id_s )
            CALL iom_get ( id_s, jpdom_unknown, 'votemper', tsdta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_s, jpdom_unknown, 'votemper', tsdta(:,:,2), ktime=ntobc2 )
            CALL iom_get ( id_s, jpdom_unknown, 'vosaline', ssdta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_s, jpdom_unknown, 'vosaline', ssdta(:,:,2), ktime=ntobc2 )
! FD move
# if defined key_dynspg_ts || defined key_dynspg_exp
	     IF( kt /= nit000 ) THEN
	          sshsdta(:,0) = sshsdta(:,1)
	     ENDIF
            CALL iom_get ( id_s, jpdom_unknown, 'vossurfh', sshsdta(:,1), ktime=ntobc1 )
            CALL iom_get ( id_s, jpdom_unknown, 'vossurfh', sshsdta(:,2), ktime=ntobc2 )
# endif
            CALL iom_close ( id_s )                                                           
            !
# if defined key_obc_tangent_vel
! FD: add tangential velocity
            CALL iom_open ( 'obcsouth_U.nc', id_s )                                          
            CALL iom_get ( id_s, jpdom_unknown, 'vozocrtx', usdta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_s, jpdom_unknown ,'vozocrtx', usdta(:,:,2), ktime=ntobc2 )
            CALL iom_close ( id_s )
# endif
            CALL iom_open ( 'obcsouth_V.nc', id_s )                                          
            CALL iom_get ( id_s, jpdom_unknown, 'vomecrty', vsdta(:,:,1), ktime=ntobc1 )
            CALL iom_get ( id_s, jpdom_unknown ,'vomecrty', vsdta(:,:,2), ktime=ntobc2 )
            CALL iom_close ( id_s )
            !
            IF( lwp .AND.  ( kt == nit000 .OR. nprint /= 0 ) )   THEN
               WRITE(numout,*)
               WRITE(numout,*) ' Read South OBC data records ', ntobc1, ntobc2
               iprint = (jpisf-jpisd+1)/20 +1
               WRITE(numout,*) ' Temperature  record 1 - printout every 3 level'
               CALL prihre( tsdta(:,:,1), jpisf-jpisd+1, jpk, 1, jpisf-jpisd+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
               WRITE(numout,*)
               WRITE(numout,*) ' Salinity  record 1 - printout every 3 level'
               CALL prihre( ssdta(:,:,1), jpisf-jpisd+1, jpk, 1, jpisf-jpisd+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
! FD move
# if defined key_dynspg_ts || defined key_dynspg_exp
               WRITE(numout,*)
               WRITE(numout,*) ' Sea surface height record 1  - printout surface level'
               CALL prihre( sshsdta(:,1), jpisf-jpisd+1, 1, 1, jpisf-jpisd+1, iprint, 1, 1, -3, 1., numout )
# endif
# if defined key_obc_tangent_vel
! FD: add tangential velocity
               WRITE(numout,*)
               WRITE(numout,*) ' Tangent velocity U record 1  - printout every 3 level'
               CALL prihre( usdta(:,:,1), jpisf-jpisd+1, jpk, 1, jpisf-jpisd+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
# endif
               WRITE(numout,*)
               WRITE(numout,*) ' Normal velocity V  record 1  - printout every 3 level'
               CALL prihre( vsdta(:,:,1), jpisf-jpisd+1, jpk, 1, jpisf-jpisd+1, iprint, &
                  &         jpk, 1, -3, 1., numout )
            ENDIF
         ENDIF

! FD move
# if defined key_dynspg_ts || defined key_dynspg_exp
         ! 2.2  After reading the velocity variables, compute depth averaged
         ! --------------------------------------------------------------------

          IF( lp_obc_east ) THEN
             ! initialisation to zero
	     IF( kt /= nit000 ) THEN 
	          ubtedta(:,0) = ubtedta(:,1)
	          vbtedta(:,0) = vbtedta(:,1)
	     ENDIF
             ubtedta(:,1:2) = 0.e0
! FD: add tangential velocity
             vbtedta(:,1:2) = 0.e0
                DO ji = nie0, nie1
                   DO jj = nje0, nje1
                      ij = jj -1 + njmpp
                      DO jk = 1, jpkm1
                         ubtedta(ij,1:2) = ubtedta(ij,1:2) + uedta(ij,jk,1:2)*fse3u(ji,jj,jk)
! FD: add tangential velocity
                         vbtedta(ij,1:2) = vbtedta(ij,1:2) + vedta(ij,jk,1:2)*fse3v(ji+1,jj,jk)
                      END DO
!                      ubtedta(ij,:) = ubtedta(ij,:) * umask(ji,jj,1)
! FD: add tangential velocity
!                      vbtedta(ij,:) = vbtedta(ij,:) * vmask(ji+1,jj,1)
                   END DO
                END DO
          ENDIF

          IF( lp_obc_west) THEN
             ! initialisation to zero
	     IF( kt /= nit000 ) THEN 
	          ubtwdta(:,0) = ubtwdta(:,1)
	          vbtwdta(:,0) = vbtwdta(:,1)
	     ENDIF
             ubtwdta(:,1:2) = 0.e0
! FD: add tangential velocity
             vbtwdta(:,1:2) = 0.e0
                DO ji = niw0, niw1
                   DO jj = njw0, njw1
                      ij = jj -1 + njmpp
                      DO jk = 1, jpkm1
                         ubtwdta(ij,1:2) = ubtwdta(ij,1:2) + uwdta(ij,jk,1:2)*fse3u(ji,jj,jk)
! FD: add tangential velocity
                         vbtwdta(ij,1:2) = vbtwdta(ij,1:2) + vwdta(ij,jk,1:2)*fse3v(ji,jj,jk)
                      END DO
!                      ubtwdta(ij,:) = ubtwdta(ij,:) * umask(ji,jj,1)
! FD: add tangential velocity
!                      vbtwdta(ij,:) = vbtwdta(ij,:) * vmask(ji,jj,1)
                   END DO
                END DO
          ENDIF

          IF( lp_obc_north) THEN
             ! initialisation to zero
! FD: add tangential velocity
	     IF( kt /= nit000 ) THEN 
	          ubtndta(:,0) = ubtndta(:,1)
	          vbtndta(:,0) = vbtndta(:,1)
	     ENDIF
             ubtndta(:,1:2) = 0.e0
             vbtndta(:,1:2) = 0.e0
                DO jj = njn0, njn1
                   DO ji = nin0, nin1
                      ii = ji -1 + nimpp
                      DO jk = 1, jpkm1
! FD: add tangential velocity
                         ubtndta(ii,1:2) = ubtndta(ii,1:2) + undta(ii,jk,1:2)*fse3u(ji,jj+1,jk)
                         vbtndta(ii,1:2) = vbtndta(ii,1:2) + vndta(ii,jk,1:2)*fse3v(ji,jj,jk)
                      END DO
! FD: add tangential velocity
!                      ubtndta(ii,:) = ubtndta(ii,:) * umask(ji,jj+1,1)
!                      vbtndta(ii,:) = vbtndta(ii,:) * vmask(ji,jj,1)
                   END DO
                END DO
          ENDIF

          IF( lp_obc_south) THEN
             ! initialisation to zero
	     IF( kt /= nit000 ) THEN
	          ubtsdta(:,0) = ubtsdta(:,1)
	          vbtsdta(:,0) = vbtsdta(:,1)
	     ENDIF
! FD: add tangential velocity
             ubtsdta(:,1:2) = 0.e0
             vbtsdta(:,1:2) = 0.e0
                DO jj = njs0, njs1
                   DO ji = nis0, nis1
                      ii = ji -1 + nimpp
                      DO jk = 1, jpkm1
! FD: add tangential velocity
                         ubtsdta(ii,1:2) = ubtsdta(ii,1:2) + usdta(ii,jk,1:2)*fse3u(ji,jj,jk)
                         vbtsdta(ii,1:2) = vbtsdta(ii,1:2) + vsdta(ii,jk,1:2)*fse3v(ji,jj,jk)
                      END DO
! FD: add tangential velocity
!                      ubtsdta(ii,:) = ubtsdta(ii,:) * umask(ji,jj,1)
!                      vbtsdta(ii,:) = vbtsdta(ii,:) * vmask(ji,jj,1)
                   END DO
                END DO
          ENDIF
# endif

      ELSE
         
         nlecto = 0        !      no reading of OBC barotropic data                         

      ENDIF                !      end of the test on the condition to read or not the files 
      
      ! 3.  Call at every time step :
      !     Linear interpolation of BCs to current time step 
      ! ----------------------------------------------------

      IF( ntobc == 1 .OR. nobc_dta == 0 )   THEN 
         zxy = 0.
      ELSE IF( ntobc == 12 )   THEN         
         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
      ELSE
#if defined key_print_WAR
WRITE(6,"(a, i5,i7,2i3,$)") "(obcdta, obc_dta) kt,isrel,ntobc1,ntobc2, zxy",kt,isrel,ntobc1,ntobc2	!! WAR
CALL FLUSH(6)									!! WAR
#endif
         zxy = (tcobc(ntobc1)-FLOAT(isrel)- tcobc(1))/(tcobc(ntobc1)-tcobc(ntobc2))
#if defined key_print_WAR
WRITE(6,*) zxy									!! WAR
CALL FLUSH(6)									!! WAR
#endif
      ENDIF

      IF( lp_obc_east )   THEN
         !  fills sfoe, tfoe, ufoe (local to each processor)
         DO jk = 1, jpkm1
! FD: increase the range
! FD            DO jj = nje0p1, nje1m1
            DO jj = nje0, nje1
               ij = jj -1 + njmpp
!!!			WRITE(6,*) "(obcdta,obc_dta) narea,kt,jpkm1,nje0,nje1,ij", narea, kt, jpkm1, nje0, nje1, ij				!! WAR
!!!			CALL FLUSH(6)				!! WAR
               sfoe(jj,jk) =  ( zxy * sedta(ij,jk,2) + &
                  &           (1.-zxy) * sedta(ij,jk,1) ) * temsk(jj,jk)
               tfoe(jj,jk) =  ( zxy * tedta(ij,jk,2) + &
                  &           (1.-zxy) * tedta(ij,jk,1) ) * temsk(jj,jk)
               ufoe(jj,jk) =  ( zxy * uedta(ij,jk,2) + &
                  &           (1.-zxy) * uedta(ij,jk,1) ) * uemsk(jj,jk)
! FD: add tangential velocity
               vfoe(jj,jk) =  ( zxy * vedta(ij,jk,2) + &
                  &           (1.-zxy) * vedta(ij,jk,1) ) * vemsk(jj,jk)
            END DO
         END DO
      ENDIF
      
      IF( lp_obc_west )   THEN
         !  fills sfow, tfow, ufow (local to each processor)
         DO jk = 1, jpkm1
! FD: increase the range
! FD            DO jj = njw0p1, njw1m1
            DO jj = njw0, njw1
               ij = jj -1 + njmpp
               sfow(jj,jk) =  ( zxy * swdta(ij,jk,2) + &
                  &           (1.-zxy) * swdta(ij,jk,1) ) * twmsk(jj,jk)
               tfow(jj,jk) =  ( zxy * twdta(ij,jk,2) + &
                  &           (1.-zxy) * twdta(ij,jk,1) ) * twmsk(jj,jk)
               ufow(jj,jk) =  ( zxy * uwdta(ij,jk,2) + &
                  &           (1.-zxy) * uwdta(ij,jk,1) ) * uwmsk(jj,jk)
! FD: add tangential velocity
               vfow(jj,jk) =  ( zxy * vwdta(ij,jk,2) + &
                  &           (1.-zxy) * vwdta(ij,jk,1) ) * vwmsk(jj,jk)
            END DO
         END DO
      ENDIF
      
      IF( lp_obc_north )   THEN
         !  fills sfon, tfon, vfon (local to each processor)
         DO jk = 1, jpkm1
! FD: increase the range
! FD            DO ji = nin0p1, nin1m1
            DO ji = nin0, nin1
               ii = ji -1 + nimpp
               sfon(ji,jk) =  ( zxy * sndta(ii,jk,2) + &
                  &           (1.-zxy) * sndta(ii,jk,1) ) * tnmsk(ji,jk)
               tfon(ji,jk) =  ( zxy * tndta(ii,jk,2) + &
                  &           (1.-zxy) * tndta(ii,jk,1) ) * tnmsk(ji,jk)
! FD: add tangential velocity
               ufon(ji,jk) =  ( zxy * undta(ii,jk,2) + &
                  &           (1.-zxy) * undta(ii,jk,1) ) * unmsk(ji,jk)
               vfon(ji,jk) =  ( zxy * vndta(ii,jk,2) + &
                  &           (1.-zxy) * vndta(ii,jk,1) ) * vnmsk(ji,jk)
            END DO
         END DO
      ENDIF
      
      IF( lp_obc_south )   THEN
         !  fills sfos, tfos, vfos (local to each processor)
         DO jk = 1, jpkm1
! FD: increase the range
! FD           DO ji = nis0p1, nis1m1
           DO ji = nis0, nis1
              ii = ji -1 + nimpp
              sfos(ji,jk) = ( zxy * ssdta(ii,jk,2) + &
                 &          (1.-zxy) * ssdta(ii,jk,1) ) * tsmsk(ji,jk)
              tfos(ji,jk) = ( zxy * tsdta(ii,jk,2) + &
                 &          (1.-zxy) * tsdta(ii,jk,1) ) * tsmsk(ji,jk)
! FD: add tangential velocity
              ufos(ji,jk) = ( zxy * usdta(ii,jk,2) + &
                 &          (1.-zxy) * usdta(ii,jk,1) ) * usmsk(ji,jk)
              vfos(ji,jk) = ( zxy * vsdta(ii,jk,2) + &
                 &          (1.-zxy) * vsdta(ii,jk,1) ) * vsmsk(ji,jk)
           END DO
         END DO             
      ENDIF

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
      INTEGER ::   inum                      ! temporary logical unit
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
            CALL ctlopn( inum, 'obceastbsf.dta', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
               &         1, numout, .TRUE., 1 )
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
            CALL ctlopn( inum, 'obcwestbsf.dta', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
               &         1, numout, .TRUE., 1 )
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
            CALL ctlopn( inum, 'obcsouthbsf.dta', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
               &         1, numout, .TRUE., 1 )
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
            CALL ctlopn( inum, 'obcnorthbsf.dta', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
               &         1, numout, .TRUE., 1 )
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
      WRITE(*,*) 'obc_dta_psi: You should not have seen this print! error?', kt
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
      INTEGER ::   id_e, id_w, id_n, id_s, fid  ! file identifiers
      INTEGER ::   itimo, iman, imois, i15
      INTEGER ::   itobcm, itobcp, itimom, itimop
      REAL(wp) ::  zxy
      INTEGER ::   isrel, ikt           ! number of seconds since 1/1/1992
      INTEGER ::   iprint              ! frequency for printouts.
! FD for swapping old time frame
      INTEGER ::   ibt1, ibt2
CHARACTER(LEN=3) :: scase						!! WAR

      !!---------------------------------------------------------------------------

! FD: this routine only interpolates in time during the barotropic subcycle
!     the barotropic OBC variables
      !!------------------------------------------------------------------------------------
      ! 1.      Initialize the time we are at. Does this every time the routine is called,
      !         excepted when nobc_dta = 0
      !

      IF( nobc_dta == 0) THEN
         itimo = 1
         zxy   = 0
	 ibt1  = 1
	 ibt2  = 1
      ELSE
         IF(ntobc == 1) THEN
            itimo = 1
            ibt1  = 1
            ibt2  = 1
         ELSE IF (ntobc == 12) THEN      !   BC are monthly
            ! we assume we have climatology in that case
            iman  = 12
            i15   = nday / 16
            imois = nmonth + i15 - 1
            IF( imois == 0 )   imois = iman
            itimo = imois
            ibt1  = 1
            ibt2  = 2
         ELSE
! FD too many printouts            IF(lwp) WRITE(numout,*) 'data other than constant or monthly',kt
            iman  = ntobc
            itimo = FLOOR( kt*rdt / (tcobc(2)-tcobc(1)))
            isrel=kt*rdt
            ibt1  = 1
            ibt2  = 2
         ENDIF
      ENDIF

      ! 3.  Call at every time step : Linear interpolation of BCs to current time step
      ! ----------------------------------------------------------------------

       IF( lk_dynspg_ts ) THEN
          isrel = (kt-1)*rdt + kbt*rdtbt

          IF( nobc_dta == 1 ) THEN
             itimo  = FLOOR(  kt*rdt    / (tcobc(2)-tcobc(1)) )
             itimom = FLOOR( (kt-1)*rdt / (tcobc(2)-tcobc(1)) )
             itimop = FLOOR( (kt+1)*rdt / (tcobc(2)-tcobc(1)) )
             itobcm = ntobc1
             itobcp = ntobc2
          ENDIF

       ELSE IF( lk_dynspg_exp ) THEN
          isrel=kt*rdt
          itobcm = ntobc1
          itobcp = ntobc2
       ENDIF

       IF( ntobc == 1 .OR. nobc_dta == 0 ) THEN
          zxy = 0.e0
       ELSE IF( ntobc == 12 ) THEN
          zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
       ELSE
#if defined key_print_WAR
WRITE(6,"(a, i5,i7,2i3,a5,$)") "(obcdta, obc_dta_bt) kt,isrel,itobcm,itobcp,scase, zxy",kt,isrel,itobcm,itobcp, scase	!! WAR
CALL FLUSH(6)									!! WAR
#endif
          zxy = (tcobc(itobcm)-FLOAT(isrel)-tcobc(1)) / (tcobc(itobcm)-tcobc(itobcp)) !! WAR
! FD case of extrapolation, switch to old time frames
	  IF( zxy < 0. ) THEN
	    ibt1 = 0
	    ibt2 = 1
	    itobcm = itobcm - 1
	    itobcp = itobcp - 1
            zxy = (tcobc(itobcm)-FLOAT(isrel)-tcobc(1)) / (tcobc(itobcm)-tcobc(itobcp)) !! WAR
	  ENDIF
#if defined key_print_WAR
WRITE(6,*) zxy									!! WAR
CALL FLUSH(6)									!! WAR
#endif
       ENDIF

       IF( lp_obc_east ) THEN           !  fills sshfoe, ubtfoe (local to each processor)
! FD: increase the range
! FD          DO jj = nje0p1, nje1m1
          DO jj = nje0, nje1
             ij = jj -1 + njmpp
             sshfoe(jj) = ( zxy * sshedta(ij,ibt2) + (1.-zxy) * sshedta(ij,ibt1) ) * temsk(jj,1)
             ubtfoe(jj) = ( zxy * ubtedta(ij,ibt2) + (1.-zxy) * ubtedta(ij,ibt1) ) * uemsk(jj,1)
! FD: add tangential velocity
             vbtfoe(jj) = ( zxy * vbtedta(ij,ibt2) + (1.-zxy) * vbtedta(ij,ibt1) ) * vemsk(jj,1)
          END DO
       ENDIF

       IF( lp_obc_west) THEN            !  fills sshfow, ubtfow (local to each processor)
! FD: increase the range
! FD          DO jj = njw0p1, njw1m1
          DO jj = njw0, njw1
             ij = jj -1 + njmpp
             sshfow(jj) = ( zxy * sshwdta(ij,ibt2) + (1.-zxy) * sshwdta(ij,ibt1) ) * twmsk(jj,1)
             ubtfow(jj) = ( zxy * ubtwdta(ij,ibt2) + (1.-zxy) * ubtwdta(ij,ibt1) ) * uwmsk(jj,1)
! FD: add tangential velocity
             vbtfow(jj) = ( zxy * vbtwdta(ij,ibt2) + (1.-zxy) * vbtwdta(ij,ibt1) ) * vwmsk(jj,1)
          END DO
       ENDIF

       IF( lp_obc_north) THEN           !  fills sshfon, vbtfon (local to each processor)
! FD: increase the range
! FD          DO ji = nin0p1, nin1m1
          DO ji = nin0, nin1
             ii = ji -1 + nimpp
             sshfon(ji) = ( zxy * sshndta(ii,ibt2) + (1.-zxy) * sshndta(ii,ibt1) ) * tnmsk(ji,1)
! FD: add tangential velocity
             ubtfon(ji) = ( zxy * ubtndta(ii,ibt2) + (1.-zxy) * ubtndta(ii,ibt1) ) * unmsk(ji,1)
             vbtfon(ji) = ( zxy * vbtndta(ii,ibt2) + (1.-zxy) * vbtndta(ii,ibt1) ) * vnmsk(ji,1)
          END DO
       ENDIF

       IF( lp_obc_south) THEN           !  fills sshfos, vbtfos (local to each processor)
! FD: increase the range
! FD          DO ji = nis0p1, nis1m1
          DO ji = nis0, nis1
             ii = ji -1 + nimpp
             sshfos(ji) = ( zxy * sshsdta(ii,ibt2) + (1.-zxy) * sshsdta(ii,ibt1) ) * tsmsk(ji,1)
! FD: add tangential velocity
             ubtfos(ji) = ( zxy * ubtsdta(ii,ibt2) + (1.-zxy) * ubtsdta(ii,ibt1) ) * usmsk(ji,1)
             vbtfos(ji) = ( zxy * vbtsdta(ii,ibt2) + (1.-zxy) * vbtsdta(ii,ibt1) ) * vsmsk(ji,1)
          END DO
       ENDIF
CALL FLUSH(6)
Call FLUSH(numout)
CALL MPPSYNC

   END SUBROUTINE obc_dta_bt

#else
   !!-----------------------------------------------------------------------------
   !!   Default option
   !!-----------------------------------------------------------------------------
   SUBROUTINE obc_dta_bt ( kt, kbt )       ! Empty routine
      !! * Arguments
      INTEGER,INTENT(in) :: kt
      INTEGER, INTENT( in ) ::   kbt         ! barotropic ocean time-step index
      WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kt
      WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kbt
   END SUBROUTINE obc_dta_bt
#endif

#else
   !!--------------------------------------------------------------------
   !!  default option  :  Dummy module    NO Open Boundary Conditions
   !!--------------------------------------------------------------------
CONTAINS
   SUBROUTINE obc_dta( kt )             ! Dummy routine
     INTEGER, INTENT (in) :: kt
     WRITE(*,*) 'obc_dta: You should not have seen this print! error?', kt
   END SUBROUTINE obc_dta
   SUBROUTINE obc_dta_bt( kt, jn)             ! Dummy routine
     INTEGER, INTENT (in) :: kt, jn
     WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kt
     WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', jn
   END SUBROUTINE obc_dta_bt
#endif

   !!=====================================================================
END MODULE obcdta
