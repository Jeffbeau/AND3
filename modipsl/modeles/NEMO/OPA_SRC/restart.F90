!!DB - 2007.12.04
!!This version has eliminated the IOIPSL code

MODULE restart
   !!======================================================================
   !!                     ***  MODULE  restart  ***
   !! Ocean restart :  write the ocean restart file
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   rst_write  : write of the restart file
   !!   rst_read   : read the restart file
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE oce             ! ocean dynamics and tracers 
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE daymod          ! calendar
   USE sol_oce         ! ocean elliptic solver
   USE zdf_oce         ! ???
   USE zdftke          ! turbulent kinetic energy scheme
   USE ice_oce         ! ice variables
   USE blk_oce         ! bulk variables
   USE flx_oce         ! sea-ice/ocean forcings variables
   USE dynspg_oce      ! free surface time splitting scheme variables
   USE cpl_oce,         ONLY : lk_cpl              !

   USE lib_ncdf          ! netCDF I/O library
#if defined key_ice_lim
   USE ice               ! for ice related I/O
#endif
   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC rst_write  ! routine called by step.F90
   PUBLIC rst_read   ! routine called by inidtr.F90
#if defined key_ice_lim
!!DB for ice
   PUBLIC rst_ice_write  ! routine called by LIM_SRC/icestp.F90
   PUBLIC rst_ice_read   ! routine called by LIM_SRC/iceini
#endif 



   !! * Module variables
   CHARACTER (len=48) ::   &
      crestart = 'initial.nc'   ! restart file name
   INTEGER :: write_count
#if defined key_ice_lim
   INTEGER ::  ice_write_count = 0
#endif
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/restart.F90,v 1.15 2006/03/10 10:55:35 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------


CONTAINS


   !!----------------------------------------------------------------------
   !!   Default option                                          NetCDF file
   !!----------------------------------------------------------------------

   SUBROUTINE rst_write( kt )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rstwrite  ***
      !!                     
      !! ** Purpose :   Write restart fields in NetCDF format
      !!
      !! ** Method  :   Write in numwrs file each nstock time step in NetCDF
      !!      file, save fields which are necessary for restart
      !!
      !! History :
      !!        !  99-11  (M. Imbard)  Original code
      !!   8.5  !  02-08  (G. Madec)  F90: Free form
      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Arguments 
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      LOGICAL ::   llbon
      CHARACTER (len=50) ::   clname, cln
      INTEGER ::   ic, jc, itime
      REAL(wp) ::   zdate0
      REAL(wp), DIMENSION( 1) ::   zfice, zfblk   ! used only in case of ice & bulk
      REAL(wp), DIMENSION(10) ::   zinfo(10)
      REAL(wp), DIMENSION(jpi,jpj) :: ztab 
      
      LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.  ! Use IOIPSL subroutines for restart output
      LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE. ! Use lib_ncdf for restart output
      INTEGER :: lncdf_stat                    ! lib_ncdf call return status
      CHARACTER(LEN=80) :: restfilename
      CHARACTER(LEN=4) :: chnum

#if defined key_agrif
       Integer :: knum
#endif
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'rst_wri : write restart.output NetCDF file'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         zfice(1) = 1.e0   ;   zfblk(1) = 1.e0
      ENDIF
      IF (kt == nit000) THEN
         write_count = 0
      END IF


      IF( MOD( kt, nstock ) == 0 .OR. kt == nitend ) THEN
         
         ! 0. Initializations
         ! ------------------

         IF(lwp) WRITE(numout,*) ' '
         IF(lwp) WRITE(numout,*) 'rst_write : write the restart file in NetCDF format ',   &
                                              'at it= ',kt,' date= ',ndastp
         IF(lwp) WRITE(numout,*) '~~~~~~~~~'

         ! Job informations
         zinfo(1) = FLOAT( no        )   ! job number
         zinfo(2) = FLOAT( kt        )   ! time-step
         zinfo(3) = FLOAT( 2 - nsolv )   ! pcg solver
         zinfo(4) = FLOAT( nsolv - 1 )   ! sor solver
         IF( lk_zdftke ) THEN
            zinfo(5) = 1.e0              ! TKE 
         ELSE
            zinfo(5) = 0.e0              ! no TKE 
         ENDIF
         zinfo(6) = FLOAT( ndastp )      ! date
         zinfo(7) = adatrj               ! number of elapsed days since the begining of the run


!!DB - code deleted ...
!         IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Are we using IOIPSL?
      ! LIB_NCDF CALLS

!DB 
!      IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Are we using lib_ncdf?
!         WRITE(UNIT=restfilename,FMT="('rstfile_',I4.4,'.nc')") write_count
!!DB
         write(chnum,'(I4.4)')write_count
         restfilename = trim(cexper)//'_rstfile_'//chnum//'.nc'

         CALL ncdf_create_restart(restfilename, lncdf_stat)
         CALL ncdf_write(restfilename, 'info', zinfo, lncdf_stat)   
         CALL ncdf_write(restfilename, 'ub', ub, nwrite, lncdf_stat) 
         CALL ncdf_write(restfilename, 'vb', vb, nwrite, lncdf_stat)  
         CALL ncdf_write(restfilename, 'tb', tb, nwrite, lncdf_stat) 
         CALL ncdf_write(restfilename, 'sb', sb, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'rotb', rotb, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'hdivb', hdivb, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'un', un, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'vn', vn, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'tn', tn, nwrite, lncdf_stat) 
         CALL ncdf_write(restfilename, 'sn', sn, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'rotn', rotn, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'hdivn', hdivn, nwrite, lncdf_stat)
         ztab(:,:) = gcx(1:jpi,1:jpj)
         CALL ncdf_write(restfilename, 'gcx', ztab, nwrite, lncdf_stat)
         ztab(:,:) = gcxb(1:jpi,1:jpj)
         CALL ncdf_write(restfilename, 'gcxb', ztab, nwrite, lncdf_stat)       
#if defined key_dynspg_rl
         CALL ncdf_write(restfilename, 'bsfb', bsfb, nwrite, lncdf_stat)! Rigid-lid formulation (bsf)
         CALL ncdf_write(restfilename, 'bsfn', bsfn, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'bsfd', bsfd, nwrite, lncdf_stat)
#else
         CALL ncdf_write(restfilename, 'sshb', sshb, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'sshn', sshn, nwrite, lncdf_stat)
# if defined key_dynspg_ts      
         CALL ncdf_write(restfilename, 'sshb_b', sshb_b, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'sshn_b', sshn_b, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'un_b', un_b, nwrite, lncdf_stat)  
         CALL ncdf_write(restfilename, 'vn_b', vn_b, nwrite, lncdf_stat)
# endif
#endif
#if defined key_zdftke   ||   defined key_esopa
         IF( lk_zdftke ) THEN
            CALL ncdf_write(restfilename, 'en', en, nwrite, lncdf_stat)   ! TKE arrays
         ENDIF
#endif
#if defined key_ice_lim
          zfice(1) = FLOAT( nfice )                                      ! Louvain La Neuve Sea Ice Model
!ADbug          CALL ncdf_write(restfilename, 'nfice', zfice, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'nfice', zfice, lncdf_stat)
         CALL ncdf_write(restfilename, 'sst_io', sst_io, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'sss_io', sss_io, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'u_io', u_io, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'v_io', v_io, nwrite, lncdf_stat)
# if defined key_coupled
         CALL ncdf_write(restfilename, 'alb_ice', alb_ice, nwrite, lncdf_stat)
# endif
#endif
#if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
         zfblk(1) = FLOAT( nfbulk )                                 ! Bulk
! ADbug         CALL ncdf_write(restfilename, 'nfbulk', zfblk, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'nfbulk', zfblk, lncdf_stat)
         CALL ncdf_write(restfilename, 'gsst', gsst, nwrite, lncdf_stat)
#endif

         CALL ncdf_write(restfilename, 'nav_lat', gphit, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'nav_lon', glamt, nwrite, lncdf_stat)
         CALL ncdf_write(restfilename, 'nav_lev', gdept, lncdf_stat)
         CALL ncdf_write(restfilename, 'time', REAL(0), 1, lncdf_stat)
         CALL ncdf_write(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)
         write_count = write_count + 1

!      END IF ! End lib_ncdf block

      ! END LIB_NCDF CALLS

   ENDIF

   END SUBROUTINE rst_write



   SUBROUTINE rst_read
     !!---------------------------------------------------------------------- 
     !!                   ***  ROUTINE rst_read  ***
     !! 
     !! ** Purpose :   Read files for restart
     !! 
     !! ** Method  :   Read the previous fields on the NetCDF file 
     !!      the first record indicates previous characterics
     !!      after control with the present run, we read :
     !!      - prognostic variables on the second record
     !!      - elliptic solver arrays 
     !!      - barotropic stream function arrays ("key_dynspg_rl" defined)
     !!        or free surface arrays 
     !!      - tke arrays (lk_zdftke=T)
     !!      for this last three records,  the previous characteristics 
     !!      could be different with those used in the present run. 
     !!
     !!   According to namelist parameter nrstdt,
     !!       nrstdt = 0  no control on the date (nit000 is  arbitrary).
     !!       nrstdt = 1  we verify that nit000 is equal to the last
     !!                   time step of previous run + 1.
     !!       In both those options, the  exact duration of the experiment
     !!       since the beginning (cumulated duration of all previous restart runs)
     !!       is not stored in the restart and is assumed to be (nit000-1)*rdt.
     !!       This is valid is the time step has remained constant.
     !!
     !!       nrstdt = 2  the duration of the experiment in days (adatrj)
     !!                    has been stored in the restart file.
     !!
     !! History :
     !!        !  99-05  (M. Imbard)  Original code
     !!   8.5  !  02-09  (G. Madec)  F90: Free form
     !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
     !!----------------------------------------------------------------------
     !! * Modules used
     USE ioipsl
     
     !! * Local declarations
     LOGICAL ::   llog
     CHARACTER (len=8 ) ::   clvnames(50)
     CHARACTER (len=32) ::   clname
     INTEGER  ::   &
          itime, ibvar,     &  !
          inum                 ! temporary logical unit
     REAL(wp) ::   zdate0, zdt, zinfo(10)
     REAL(wp) ::   zdept(jpk), zlamt(jpi,jpj), zphit(jpi,jpj)
     REAL(wp), DIMENSION(jpi,jpj) :: ztab 
#   if defined key_ice_lim
     INTEGER  ::   ios1, ji, jj, jn
     REAL(wp) ::   zfice(1)
#   endif
# if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
     INTEGER  ::   ios2, jk
     REAL(wp) ::   zfblk(1)
#   endif
     
     LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.   ! Use IOIPSL subroutines for restart input
     LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE.  ! Use lib_ncdf subroutines for restart input
     INTEGER :: lncdf_stat                     ! lib_ncdf routine return status
     REAL :: tstp
     CHARACTER(LEN=40) :: restfilename="restart.nc" !Default filename 
     
     
     !!----------------------------------------------------------------------
     !!  OPA 8.5, LODYC-IPSL (2002)
     !!----------------------------------------------------------------------
     clname = 'restart'
#if defined key_agrif       
     inum = Agrif_Get_Unit()
     If(.NOT. Agrif_root() ) clname = TRIM(Agrif_CFixed())//'_'//TRIM(clname)
#endif 
     
     IF(lwp) WRITE(numout,*)
     IF(lwp) WRITE(numout,*) 'rst_read : read the NetCDF restart file'
     IF(lwp) WRITE(numout,*) '~~~~~~~~'
     
     IF(lwp) WRITE(numout,*) ' Info on the present job : '
     IF(lwp) WRITE(numout,*) '   job number          : ', no
     IF(lwp) WRITE(numout,*) '   time-step           : ', nit000
     IF(lwp) WRITE(numout,*) '   solver type         : ', nsolv
     IF( lk_zdftke ) THEN
        IF(lwp) WRITE(numout,*) '   tke option          : 1 '
     ELSE
        IF(lwp) WRITE(numout,*) '   tke option          : 0 '
     ENDIF
     IF(lwp) WRITE(numout,*) '   date ndastp         : ', ndastp
     IF(lwp) WRITE(numout,*)
     
     ! Time domain : restart
     ! -------------------------
     
     IF(lwp) WRITE(numout,*)
     IF(lwp) WRITE(numout,*)
     IF(lwp) WRITE(numout,*) ' *** restart option'
     SELECT CASE ( nrstdt )
     CASE ( 0 ) 
        IF(lwp) WRITE(numout,*) ' nrstdt = 0 no control of nit000'
     CASE ( 1 ) 
        IF(lwp) WRITE(numout,*) ' nrstdt = 1 we control the date of nit000'
     CASE ( 2 )
        IF(lwp) WRITE(numout,*) ' nrstdt = 2 the date adatrj is read in restart file'
     CASE DEFAULT
        IF(lwp) WRITE(numout,*) '  ===>>>> nrstdt not equal 0, 1 or 2 : no control of the date'
        IF(lwp) WRITE(numout,*) ' =======                   ========='
     END SELECT
     
     itime = 0
     llog  = .FALSE.
     zlamt(:,:) = 0.e0
     zphit(:,:) = 0.e0
     zdept(:)   = 0.e0

!!DB -- code deleted ...
!     IF(USE_IOIPSL .EQV. .TRUE.) THEN ! Do we want to use IOIPSL?

     ! LIB_NCDF CALLS
!!DB
!     IF(USE_LIB_NCDF .EQV. .TRUE.) THEN ! Do we want to use lib_ncdf?
!        restfilename="restart.nc"
     CALL ncdf_read(restfilename, 'info', zinfo, lncdf_stat)
!!DB -- check if found file
     if(lwp .AND. lncdf_stat /= 0)then
        write(numout,*)'                '
        write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(numout,*)'STOP: Problem reading restart file ', restfilename
        write(numout,*)'narea = ', narea
        write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(numout,*)'                '
        stop
     endif

     
     ! NON LIB_NCDF RELATED CODE
     IF(lwp) WRITE(numout,*)
     IF(lwp) WRITE(numout,*) ' Info on the restart file read : '
     IF(lwp) WRITE(numout,*) '   job number          : ', NINT( zinfo(1) )
     IF(lwp) WRITE(numout,*) '   time-step           : ', NINT( zinfo(2) )
     IF(lwp) WRITE(numout,*) '   solver type         : ', NINT( zinfo(4) ) + 1
     IF(lwp) WRITE(numout,*) '   tke option          : ', NINT( zinfo(5) )
     IF(lwp) WRITE(numout,*) '   date ndastp         : ', NINT( zinfo(6) )
     IF(lwp) WRITE(numout,*)
     
     ! Control of date
     IF( nit000 - NINT( zinfo(2) )  /= 1 .AND. nrstdt /= 0 ) THEN
        IF(lwp) WRITE(numout,cform_err)
        IF(lwp) WRITE(numout,*) ' ===>>>> : problem with nit000 for the restart'
        IF(lwp) WRITE(numout,*) ' verify the restart file or rerun with nrstdt = 0 (namelist)'
        nstop = nstop + 1
     ENDIF
     
     ! re-initialisation of  adatrj0
     adatrj0 =  ( FLOAT( nit000-1 ) * rdttra(1) ) / rday
     
     IF ( nrstdt == 2 ) THEN
        !                             by default ndatsp has been set to ndate0 in dom_nam
        !                             ndate0 has been read in the namelist (standard OPA 8)
        !                             here when nrstdt=2 we keep the  final date of previous run
        ndastp = NINT( zinfo(6) )
        adatrj0 =  zinfo(7)
     ENDIF
     ! END NON LIB_NCDF RELATED CODE
     
     ! End of job info stuff
     CALL ncdf_read(restfilename, 'ub', ub, nwrite, lncdf_stat)   
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'vb', vb, nwrite, lncdf_stat)  
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'tb', tb, nwrite, lncdf_stat) 
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'sb', sb, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'rotb', rotb, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'hdivb', hdivb, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'un', un, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'vn', vn, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'tn', tn, nwrite, lncdf_stat) 
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'sn', sn, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'rotn', rotn, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'hdivn', hdivn, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'gcxb', ztab, nwrite, lncdf_stat)    
     if(lncdf_stat /= 0) stop
     gcxb(1:jpi,1:jpj) = ztab(:,:) 
     CALL ncdf_read(restfilename, 'gcx', ztab, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     gcx(1:jpi,1:jpj) = ztab(:,:)   
#if defined key_dynspg_rl
     CALL ncdf_read(restfilename, 'bsfb', bsfb, nwrite, lncdf_stat)! Rigid-lid formulation (bsf)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'bsfn', bsfn, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'bsfd', bsfd, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
#else
     CALL ncdf_read(restfilename, 'sshb', sshb, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'sshn', sshn, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
# if defined key_dynspg_ts      
     CALL ncdf_read(restfilename, 'sshb_b', sshb_b, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'sshn_b', sshn_b, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'un_b', un_b, nwrite, lncdf_stat)  
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'vn_b', vn_b, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
# endif
#endif
#if defined key_zdftke   ||   defined key_esopa
     IF( lk_zdftke ) THEN
        IF( NINT( zinfo(5) ) == 1 ) THEN                                ! Read tke arrays
           CALL ncdf_read(restfilename, 'en', en, nwrite, lncdf_stat)
           ln_rstke = .FALSE.
        ELSE
           IF(lwp) WRITE(numout,*) ' ===>>>> : the previous restart file didnot used  tke scheme'
           IF(lwp) WRITE(numout,*) ' =======                ======='
           nrstdt = 2
           ln_rstke = .TRUE.
        ENDIF
     ENDIF
#endif
# if defined key_ice_lim
      ! Louvain La Neuve Sea Ice Model
      ios1 = 0
      DO jn = 1, 30
         IF( clvnames(jn) == 'nfice' )  ios1 = 1
      END DO
      IF( ios1 == 1 ) THEN
!ADbug          CALL ncdf_read(restfilename, 'nfice', zfice, nwrite, lncdf_stat)
         CALL ncdf_read(restfilename, 'nfice', zfice, lncdf_stat)
         CALL ncdf_read(restfilename, 'sst_io', sst_io, nwrite, lncdf_stat)
         CALL ncdf_read(restfilename, 'sss_io', sss_io, nwrite, lncdf_stat)
         CALL ncdf_read(restfilename, 'u_io', u_io, nwrite, lncdf_stat)
         CALL ncdf_read(restfilename, 'v_io', v_io, nwrite, lncdf_stat)
# if defined key_coupled
         CALL ncdf_read(restfilename, 'alb_ice', alb_ice, nwrite, lncdf_stat)
# endif
      ENDIF
      IF( zfice(1) /= FLOAT(nfice) .OR. ios1 == 0 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'rst_read :  LLN sea Ice Model => Ice initialization'
         IF(lwp) WRITE(numout,*)
         sst_io(:,:) = ( nfice-1 )*( tn(:,:,1) + rt0 )          !!bug a explanation is needed here!
         sss_io(:,:) = ( nfice-1 )*  sn(:,:,1)
         DO jj = 2, jpj
            DO ji = 2, jpi
               u_io(ji,jj) = ( nfice-1 ) * 0.5 * ( un(ji-1,jj  ,1) + un(ji-1,jj-1,1) )
               v_io(ji,jj) = ( nfice-1 ) * 0.5 * ( vn(ji  ,jj-1,1) + vn(ji-1,jj-1,1) )
            END DO
         END DO
#    if defined key_coupled
         alb_ice(:,:) = 0.8 * tmask(:,:,1)
#    endif
      ENDIF
# endif
# if defined key_flx_bulk_monthly || defined key_flx_bulk_daily
      ! Louvain La Neuve Sea Ice Model
      ios2 = 0
      DO jk = 1, 30
         IF( clvnames(jk) == 'nfbulk' )  ios2 = 1
      END DO
      IF( ios2 == 1 ) THEN
!ADbug          CALL ncdf_read(restfilename, 'nfbulk', zfblk, nwrite, lncdf_stat)
         CALL ncdf_read(restfilename, 'nfbulk', zfblk, lncdf_stat)
         CALL ncdf_read(restfilename, 'gsst', gsst, nwrite, lncdf_stat)
      ENDIF
      IF( zfblk(1) /= FLOAT(nfbulk) .OR. ios2 == 0 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'rst_read :  LLN sea Ice Model => Ice initialization'
         IF(lwp) WRITE(numout,*)
         gsst(:,:) = 0.
         gsst(:,:) = gsst(:,:) + ( nfbulk-1 )*( tn(:,:,1) + rt0 )
      ENDIF
# endif

     CALL ncdf_read(restfilename, 'nav_lat', zlamt, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'nav_lon', zphit, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'nav_lev', gdept, lncdf_stat)
     if(lncdf_stat /= 0) stop
     tstp = REAL(itime)
     CALL ncdf_read(restfilename, 'time', tstp, 1, lncdf_stat)
     if(lncdf_stat /= 0) stop
     itime = NINT(tstp)
     CALL ncdf_read(restfilename, 'time_steps', zdate0, 1, lncdf_stat)
     if(lncdf_stat /= 0) stop
  ! END LIB_NCDF CALLS
  
  ! In case of restart with neuler = 0 then put all before fields = to now fields
  IF ( neuler == 0 ) THEN
     tb(:,:,:)=tn(:,:,:)
     sb(:,:,:)=sn(:,:,:)
     ub(:,:,:)=un(:,:,:)
     vb(:,:,:)=vn(:,:,:)
     rotb(:,:,:)=rotn(:,:,:)
     hdivb(:,:,:)=hdivn(:,:,:)
#if defined key_dynspg_rl
     ! rigid lid
     bsfb(:,:)=bsfn(:,:)
#else
     ! free surface formulation (eta)
     sshb(:,:)=sshn(:,:)
#endif
  ENDIF
  
END SUBROUTINE rst_read





#ifdef key_ice_lim
!!DB 2008.05.22 
!!NB: called when needed by icestp
!!         IF( MOD( numit, nstock ) == 0 .OR. numit == nlast ) THEN
!            CALL lim_rst_write( numit )                              ! Ice restart file !
!         ENDIF

   SUBROUTINE rst_ice_write( niter )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rst ice write  ***
      !!                     
      !! ** Purpose :   Write lim ice restart fields in NetCDF format
      !!
      !! ** Method  :   Write in numwrs file each nstock time step in NetCDF
      !!      file, save fields which are necessary for restart
      !!
      !! History : DB  2008.05.22
      !!
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl

      !! * Arguments 
      INTEGER, INTENT( in ) ::   niter         ! ice iteration 

      !! * Local declarations
      LOGICAL ::   llbon
      CHARACTER (len=50) ::   clname, cln
      INTEGER ::   ic, jc, itime, it0
      REAL(wp), DIMENSION(2) ::   zinfo(2)
      REAL(wp), DIMENSION(jpi,jpj) :: ztab 
      REAL(wp) ::  zsec, zdate0, zdt
      
      LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.  ! Use IOIPSL subroutines for restart output
      LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE. ! Use lib_ncdf for restart output
      INTEGER :: lncdf_stat                    ! lib_ncdf call return status
      CHARACTER(LEN=80) :: restfilename
      CHARACTER(LEN=4) :: chnum

      REAL(wp),DIMENSION(jpi,jpj,35) :: zmoment
      INTEGER :: ji, jj

#if defined key_agrif
       Integer :: knum
#endif
      !!----------------------------------------------------------------------


! Job informations
      it0      = niter
      zinfo(1) = FLOAT( nfice  )  ! coupling frequency OPA ICELLN  nfice
      zinfo(2) = FLOAT( it0   )   ! iteration number

      zsec     = 0.e0
      itime    = 0
!      zdept(1) = 0.e0   !! this is default z-level -- handled in ncdf_create_ice_restart()
      zdt      = rdt_ice * nstock  !!does not seem to enter restart file

      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) 'ice_rst_write : write the ice restart file in NetCDF format ',   &
           'at it= ',niter,' date= ',ndastp
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'


!!DB
      write(chnum,'(I4.4)')ice_write_count
      restfilename = trim(cexper)//'_ice_rstfile_'//chnum//'.nc'

      CALL ncdf_create_ice_restart(restfilename, lncdf_stat)
      CALL ncdf_write(restfilename, 'info', zinfo, lncdf_stat)   

      CALL ncdf_write(restfilename, 'nav_lat', gphit, nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'nav_lon', glamt, nwrite, lncdf_stat)
!      CALL ncdf_write(restfilename, 'nav_lev', real(0), nwrite, lncdf_stat)   !!nav_lev treated as a scalar
      CALL ncdf_write(restfilename, 'nav_lev', REAL(0), 1, lncdf_stat)   !!nav_lev treated as a scalar
      CALL ncdf_write(restfilename, 'time', REAL(0), 1, lncdf_stat)
      CALL ncdf_write(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)

      CALL ncdf_write(restfilename, 'hicif', hicif, nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'hsnif', hsnif   , nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'frld' , frld    , nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'sist' , sist    , nwrite, lncdf_stat)
# if defined key_coupled
      CALL ncdf_write(restfilename, 'albege', albege   , nwrite, lncdf_stat)
# endif, nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'tbif'  , tbif     , nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'u_ice' , u_ice, nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'v_ice' , v_ice , nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'gtaux' , gtaux , nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'gtauy' , gtauy , nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'qstoif', qstoif, nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'fsbbq' , fsbbq , nwrite, lncdf_stat)


!!DB load zmoment
      DO jj = 1, jpj  
         DO ji = 1, jpi
            zmoment(ji,jj,1)  = sxice(ji,jj)
            zmoment(ji,jj,2)  = syice(ji,jj)
            zmoment(ji,jj,3)  = sxxice(ji,jj)
            zmoment(ji,jj,4)  = syyice(ji,jj)
            zmoment(ji,jj,5)  = sxyice(ji,jj)
            zmoment(ji,jj,6)  = sxsn(ji,jj)
            zmoment(ji,jj,7)  = sysn(ji,jj)
            zmoment(ji,jj,8)  = sxxsn(ji,jj)
            zmoment(ji,jj,9)  = syysn(ji,jj)
            zmoment(ji,jj,10) = sxysn(ji,jj)
            zmoment(ji,jj,11) = sxa(ji,jj)
            zmoment(ji,jj,12) = sya(ji,jj)
            zmoment(ji,jj,13) = sxxa(ji,jj)
            zmoment(ji,jj,14) = syya(ji,jj)
            zmoment(ji,jj,15) = sxya(ji,jj)
            zmoment(ji,jj,16) = sxc0(ji,jj)
            zmoment(ji,jj,17) = syc0(ji,jj)
            zmoment(ji,jj,18) = sxxc0(ji,jj)
            zmoment(ji,jj,19) = syyc0(ji,jj)
            zmoment(ji,jj,20) = sxyc0(ji,jj)
            zmoment(ji,jj,21) = sxc1(ji,jj)
            zmoment(ji,jj,22) = syc1(ji,jj)
            zmoment(ji,jj,23) = sxxc1(ji,jj)
            zmoment(ji,jj,24) = syyc1(ji,jj)
            zmoment(ji,jj,25) = sxyc1(ji,jj)
            zmoment(ji,jj,26) = sxc2(ji,jj)
            zmoment(ji,jj,27) = syc2(ji,jj)
            zmoment(ji,jj,28) = sxxc2(ji,jj)
            zmoment(ji,jj,29) = syyc2(ji,jj)
            zmoment(ji,jj,30) = sxyc2(ji,jj)
            zmoment(ji,jj,31) = sxst(ji,jj)
            zmoment(ji,jj,32) = syst(ji,jj)
            zmoment(ji,jj,33) = sxxst(ji,jj)
            zmoment(ji,jj,34) = syyst(ji,jj)
            zmoment(ji,jj,35) = sxyst(ji,jj)
         END DO
      END DO
!!DB probably should call mppsync here ??? NB: Works if I do not
!      CALL MPI_BARRIER(MPI_COMM_WORLD, lncdf_stat)
      CALL ncdf_write(restfilename, 'moment', zmoment, nwrite, lncdf_stat)

      ice_write_count = ice_write_count + 1


   END SUBROUTINE rst_ice_write


!!DB: 2008.05.27 read ice restart 
!!Tested to run OK, but there is no proof that it did the right thing
!!NB; In general, there is no proof that ice restarts are perfect
!!Similar to old lim_rst_read()
!!NB: default restfilename = restart_ice.nc
!!NB: Agrif stuff not done yet
   SUBROUTINE rst_ice_read( niter )
      !-----------------------------------------------------------------------
      ! Arguments
      INTEGER  ::   niter        ! number of iteration

      !- dummy variables :
      CHARACTER(len=45)  ::  restfilename = 'restart_ice.nc'
      INTEGER :: &
        ji, jj, lncdf_stat
      INTEGER :: &
         inumrst, it0, it1, itime, ibvar, ifice
      LOGICAL :: &
         llog
      REAL(wp),DIMENSION(jpi,jpj) :: &
         zlamt, zphit
      REAL(wp),DIMENSION(jpi,jpj,35) :: &
         zmoment
      REAL(wp),DIMENSION(1) :: &
         zdept
      REAL(wp),DIMENSION(2) :: &
         zinfo
      REAL(wp) :: &
         zdate0, zdt
      CHARACTER ( len = 10 ) ::  &
         clvnames(60)       



#if defined key_agrif
      if ( .NOT. Agrif_Root() ) then
         restfilename= TRIM(Agrif_CFixed())//'_'//TRIM(restfilename)
      endif
#endif

     CALL ncdf_read(restfilename, 'info', zinfo, lncdf_stat)
!!DB -- check if found file
     if(lwp .AND. lncdf_stat /= 0)then
        write(numout,*)'                '
        write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(numout,*)'STOP: Problem reading ice restart file ', restfilename
        write(numout,*)'narea = ', narea
        write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(numout,*)'                '
        stop
     endif



      !Initialisations
      inumrst    = 71
      it0        = nit000
      itime      = 0
      llog       = .FALSE.
      zlamt(:,:) = 0.
      zphit(:,:) = 0.
      zdept(1)   = 0.

      ifice   = INT( zinfo(1) )
      it1     = INT( zinfo(2) )

!!DB Note that old code gets varnames but otherwise does not use this info 
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'lim_rst_read : READ restart file name ', restfilename,' at time step : ', it1
      ENDIF
      
      !Control of date
!!DB: 2009.09.30
      IF( ( it0 - it1 ) /= 1 .AND. ABS( nrstdt ) == 1 ) THEN
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) 'DBG: lim_rst_read nrstdt=1 and nit000 /= zinfo(2)'
            WRITE(numout,*) 'DBG: This case covered in restart() ====> Continue'
!            WRITE(numout,*) 'lim_rst_read ===>>>> : problem with nit000 for the restart'
!            WRITE(numout,*) '   we stop. verify the file or rerun with the value  0 for the'
!            WRITE(numout,*) '   control of time parameter  nrstdt'
!            nstop = nstop + 1
         ENDIF
      ENDIF

!!DBG: writes from rst_ice_write()
     CALL ncdf_read(restfilename, 'nav_lat', zlamt, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop
     CALL ncdf_read(restfilename, 'nav_lon', zphit, nwrite, lncdf_stat)
     if(lncdf_stat /= 0) stop

!!DBG ignore these vars to start
!      CALL ncdf_read(restfilename, 'nav_lev', zdept, 1, lncdf_stat)   !!nav_lev treated as a scalar
!      CALL ncdf_read(restfilename, 'time', REAL(0), 1, lncdf_stat)
!      CALL ncdf_read(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)

      CALL ncdf_read(restfilename, 'hicif', hicif, nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'hsnif', hsnif   , nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'frld' , frld    , nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'sist' , sist    , nwrite, lncdf_stat)
# if defined key_coupled
      CALL ncdf_read(restfilename, 'albege', albege   , nwrite, lncdf_stat)
# endif, nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'tbif'  , tbif     , nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'u_ice' , u_ice, nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'v_ice' , v_ice , nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'gtaux' , gtaux , nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'gtauy' , gtauy , nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'qstoif', qstoif, nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'fsbbq' , fsbbq , nwrite, lncdf_stat)
      CALL ncdf_read(restfilename, 'moment', zmoment, nwrite, lncdf_stat)

!!DB: 2009.09.30
!      niter = it1
      niter = nit000

      DO jj = 1, jpj
         DO ji = 1, jpi
            sxice(ji,jj)  = zmoment(ji,jj,1)
            syice(ji,jj)  = zmoment(ji,jj,2)
            sxxice(ji,jj) = zmoment(ji,jj,3)
            syyice(ji,jj) = zmoment(ji,jj,4)
            sxyice(ji,jj) = zmoment(ji,jj,5)
            sxsn(ji,jj)   = zmoment(ji,jj,6)
            sysn(ji,jj)   = zmoment(ji,jj,7)
            sxxsn(ji,jj)  = zmoment(ji,jj,8)
            syysn(ji,jj)  = zmoment(ji,jj,9)
            sxysn(ji,jj)  = zmoment(ji,jj,10)
            sxa(ji,jj)    = zmoment(ji,jj,11)
            sya(ji,jj)    = zmoment(ji,jj,12)
            sxxa(ji,jj)   = zmoment(ji,jj,13)
            syya(ji,jj)   = zmoment(ji,jj,14)
            sxya(ji,jj)   = zmoment(ji,jj,15)
            sxc0(ji,jj)   = zmoment(ji,jj,16)
            syc0(ji,jj)   = zmoment(ji,jj,17)
            sxxc0(ji,jj)  = zmoment(ji,jj,18)
            syyc0(ji,jj)  = zmoment(ji,jj,19)
            sxyc0(ji,jj)  = zmoment(ji,jj,20)
            sxc1(ji,jj)   = zmoment(ji,jj,21)
            syc1(ji,jj)   = zmoment(ji,jj,22)
            sxxc1(ji,jj)  = zmoment(ji,jj,23)
            syyc1(ji,jj)  = zmoment(ji,jj,24)
            sxyc1(ji,jj)  = zmoment(ji,jj,25)
            sxc2(ji,jj)   = zmoment(ji,jj,26)
            syc2(ji,jj)   = zmoment(ji,jj,27)
            sxxc2(ji,jj)  = zmoment(ji,jj,28)
            syyc2(ji,jj)  = zmoment(ji,jj,29)
            sxyc2(ji,jj)  = zmoment(ji,jj,30)
            sxst(ji,jj)   = zmoment(ji,jj,31)
            syst(ji,jj)   = zmoment(ji,jj,32)
            sxxst(ji,jj)  = zmoment(ji,jj,33)
            syyst(ji,jj)  = zmoment(ji,jj,34)
            sxyst(ji,jj)  = zmoment(ji,jj,35)
         END DO
      END DO

      
   END SUBROUTINE rst_ice_read



#endif   !!key_ice_lim


!!=====================================================================
END MODULE restart

