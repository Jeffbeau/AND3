!!DB: 2010.08 -- Module for simple passive tracer: key_PASSIVE_TRACER_01 
!! Contains most of the needed routines.

!!Note the sequence: (local) initrc.F90 USES this module if 
!! key_PASSIVE_TRACER_01 is defined. A consequence is that the model will use 
!! the subroutines in this module instead of other versions. A key example
!! is ini_trc() which initializes the tracer module. 

!!NB: I may not have added the north open boundary to all relevant routines so
!!     BEWARE !!!! 
!!(Also, note that must have north boundary for these code frags to work)


!!See lib_bgcm_01.F90 for more comments

#if defined key_passivetrc

MODULE lib_ptracer_01
   !!================================================
   !!
   !! Initialization of the BGCM 
   !!================================================
   !!--------------------------------------------------------------
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   USE lib_ncdf

   
   IMPLICIT NONE
   PRIVATE

!!DB: 2010.08 : kept wmask; computed in ptracer_ini
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk+1) ::  wmask

!!DB
!!Shared module variables, default values. NB: they are reset below
   !! Advection
   LOGICAL, PUBLIC ::   &
      ln_trcadv_cen2   = .FALSE. ,  & !!: 2nd order centered scheme flag
      ln_trcadv_tvd    = .TRUE. ,  &  !: TVD scheme flag
      ln_trcadv_muscl  = .FALSE. ,  &  !: MUSCL scheme flag
      ln_trcadv_muscl2 = .FALSE. ,  &  !: MUSCL2 scheme flag
      ln_trcadv_smolar = .FALSE.        !: Smolarkiewicz scheme flag (DB: must be .FALSE.)

   !! Lateral diffusion
   LOGICAL , PUBLIC ::              & !!: ** lateral mixing namelist (nam_trcldf) **
      ln_trcldf_diff  = .FALSE. ,   &  !: flag of perform or not the lateral diff.
      ln_trcldf_lap   = .TRUE.  ,   &  !: laplacian operator
      ln_trcldf_bilap = .FALSE. ,   &  !: bilaplacian operator
      ln_trcldf_level = .FALSE. ,   &  !: iso-level direction
      ln_trcldf_hor   = .TRUE. ,   &  !: horizontal (geopotential) direction
      ln_trcldf_iso   = .FALSE.         !: iso-neutral direction

   LOGICAL , PUBLIC ::              & !!: flag of the lateral diff. scheme used
      l_trcldf_lap         ,        &  !: iso-level laplacian operator
      l_trcldf_bilap       ,        &  !: iso-level bilaplacian operator
      l_trcldf_bilapg      ,        &  !: geopotential bilap. (s-coord)
      l_trcldf_iso         ,        &  !: iso-neutral laplacian or horizontal lapacian (s-coord)
      l_trczdf_iso         ,        &  !: idem for the vertical component
      l_trczdf_iso_vo      ,        &  !: idem with vectopt_memory
      l_trcldf_iso_zps                 !: iso-neutral laplacian (partial steps)

   !! Vertical diffusion
   LOGICAL , PUBLIC ::              & !!: nam_trczdf: vertical diffusion
      ln_trczdf_exp = .FALSE.          !: explicit vertical diffusion scheme flag

   INTEGER, PUBLIC ::               & !!: namzdf:  vertical diffusion
      n_trczdf_exp = 3                 !: number of sub-time step (explicit time stepping)

   LOGICAL, PUBLIC ::               &  !:
      l_trczdf_exp     = .FALSE. ,  &  !: explicit vertical diffusion
      l_trczdf_imp     = .FALSE.       !: implicit vertical diffusion

!!DB: not accounted for yet
#if defined key_trcdmp
   !! Newtonian damping
   INTEGER  , PUBLIC ::             & !!: * newtonian damping namelist (nam_trcdmp) *
      ndmptr   =   -1 ,             &  !: = 0/-1/'latitude' for damping over tracers
      ndmpftr  =    2 ,             &  !: = 1 create a damping.coeff NetCDF file 
      nmldmptr =    0                  !: = 0/1/2 flag for damping in the mixed layer

   REAL(wp) , PUBLIC ::             & !!:  * newtonian damping namelist *
      sdmptr   =   50.,             &  !: surface time scale for internal damping (days)
      bdmptr   =  360.,             &  !: bottom time scale for internal damping (days)
      hdmptr   =  800.                 !: depth of transition between sdmp and bdmp (meters)
#endif

!!DB -- arrays of boundary trn values
   REAL(wp), PUBLIC :: e_trn(jpj,jpk,jptra), w_trn(jpj,jpk,jptra), s_trn(jpi,jpk,jptra), &
        n_trn(jpi,jpk,jptra)
   REAL(wp), PUBLIC :: e_trb(jpj,jpk,jptra), w_trb(jpj,jpk,jptra), s_trb(jpi,jpk,jptra), &
        n_trb(jpi,jpk,jptra)
!
   CHARACTER (len=80), PUBLIC :: ptracer_fname
   
   !! * Accessibility
   PUBLIC ini_trc, ptracer_setup, ptracer_ini, extract_boundary_vals, &
        assign_boundary_vals, update_boundary_vals

   
CONTAINS

!!DB: routine ini_trc replaces regular ini_trc when key_PTRACER_01 is defined
!!
  SUBROUTINE ini_trc
      !!---------------------------------------------------------------------
      !!
      !!                       ROUTINE ini_trc
      !!                     ******************
      !!
      !!  PURPOSE :  initialize the PTRACER_01 model
      !!  DB: Replaces regular ini_trc when key_PTRACER_01 is defined
      !!---------------------------------------------------------------------


      !! 0.b PRINT the number of tracer
      !! ------------------------------

    IF(lwp) WRITE(numout,*) ' '
    IF(lwp) WRITE(numout,*) ' *** number of passive tracer jptra = ',jptra
    IF(lwp) WRITE(numout,*) ' '
    
    call ptracer_setup 
    
    call ptracer_ini
    
!!If restarting (not thought about yet) assume that the following call
!!will overwrite anything wrong that ptracer_ini might have done
    if( lrsttr ) THEN
!         CALL ptracer_rst       !!...TO DO ...
    endif

  END SUBROUTINE ini_trc


!!DB: Replace trc_lec, trc_trp_lec, trc_trp_ctl ... with 1 routine
!!One of the things I do in this routine is force the same adv/diff schemes
!!as are used by the TS module.
!!Initialized netcdf output file as well
!!2010.08 -- NB: Because I force the use of certain schemes (e.g. adv/diff), much of the
!!           below may not be necessary. Cleaning this up is left for a future exercise.

   SUBROUTINE ptracer_setup
     USE traadv_ctl
     USE ldftra_oce
     USE zdf_oce, ONLY : n_zdfexp
     USE lib_ncdf

      !!                     ******************
      !! local declarations
      !! ==================

      INTEGER ::  ji, status
      CHARACTER (len=32) :: clname

!!DB: These are all that are needed for PTRACER
      namelist/nattrc/nwritetrc,lrsttr,nrsttr
      namelist/natnum/ndttrc,rsc,rtrn,ncortrc,crosster
!!DB: Modified
      NAMELIST/namtrcldf/ ahtrc0, trcrat

!!DB: not yet
#if defined key_trcdmp
      NAMELIST/namtrcdmp/ ndmptr, ndmpftr, nmldmptr, sdmptr, bdmptr, hdmptr
#endif
      !!----------------------------------------------------------------------


!!DB: must do this as this routine is not called until step() is called
!!    Another example of the OPA logic that initializes various params and
!!    routines at unpredictable areas/times of the code. 
      call tra_adv_ctl 

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' ROUTINE PTRACER_SETUP'
         WRITE(numout,*) ' **************'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist for passive tracer'
         WRITE(numout,*) ' ***************************'
         WRITE(numout,*) ' '
      ENDIF

      numnat=80
!!DB: NB
      clname='namelist.PTRACER_01'
      OPEN( numnat, FILE= clname, FORM='formatted', STATUS = 'old')

!! initialization from namelist file
      !! ----------------------------------------------
      !! 1.0 namelist nattrc :

      nwritetrc = 10
      lrsttr=.FALSE.
      nrsttr = 0
      REWIND(numnat)
      READ(numnat,nattrc)

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) 'nattrc'
         WRITE(numout,*) ' '
         WRITE(numout,*)          &
            ' frequency of outputs for passive tracers nwritetrc = '    &
            ,nwritetrc  
         WRITE(numout,*) ' restart LOGICAL for passive tr. lrsttr = ',   &
            &         lrsttr
         WRITE(numout,*) ' control of time step for p. tr. nrsttr = ',   & 
            &         nrsttr
         WRITE(numout,*) ' '
      ENDIF

      !! 1.1 namelist natnum :
      !! ---------------------
      rsc=1.
      rtrn=1.e-15
      ncortrc=1
      ndttrc=4
      crosster=.FALSE.

      REWIND(numnat)
      READ(numnat,natnum)

!!Compute the first time step of tracer model
      nittrc000 = nit000 + ndttrc - 1

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) 'natnum'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' tuning coefficient              rsc     = ',    &
            rsc
         WRITE(numout,*) ' truncation value                rtrn    = ',    &
            rtrn
         WRITE(numout,*) ' number of corrective phase      ncortrc = ',    &
            ncortrc
         WRITE(numout,*) ' time step freq. for pass. trac. ndttrc  = ',    &
            ndttrc
         WRITE(numout,*) ' 1st time step for pass. trac. nittrc000 = ',    &
            nittrc000
         WRITE(numout,*) ' computes or not crossterms    crosster  = ',    &
            crosster
      ENDIF


      !! namelist of transport
      !! ---------------------
!!DB: Force tracer code to use the same adv/diff schemes as used for TS
      ln_trcadv_cen2   = ln_traadv_cen2  
      ln_trcadv_tvd    = ln_traadv_tvd   
      ln_trcadv_muscl  = ln_traadv_muscl
      ln_trcadv_muscl2 = ln_traadv_muscl2 

      l_trcldf_bilapg  = l_traldf_bilapg 
      l_trcldf_bilap   = l_traldf_bilap  
      l_trcldf_iso     = l_traldf_iso    
      l_trcldf_iso_zps = l_traldf_iso_zps
      l_trcldf_lap     = l_traldf_lap    

      l_trczdf_exp     = l_trazdf_exp
      if(l_trczdf_exp) then
         l_trczdf_imp = .false.
!!DBG: This is not debugged, so if using explicit scheme then you should check
!!     the following assignment
         n_trczdf_exp = n_zdfexp         
      else
         l_trczdf_imp = .true.
      endif
!!DB: I doubt that this will ever be .TRUE.; if so then there may
!!    be a problem with the above  ----->  refer to original trctrp_ctl.F90
      l_trczdf_iso_vo  = l_trazdf_iso_vo

      !!DB: keep this output, at least for a while
      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'DB: choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*)
         WRITE(numout,*) '             2nd order advection scheme     ln_trcadv_cen2   = ', ln_trcadv_cen2
         WRITE(numout,*) '             TVD advection scheme           ln_trcadv_tvd    = ', ln_trcadv_tvd
         WRITE(numout,*) '             MUSCL  advection scheme        ln_trcadv_muscl  = ', ln_trcadv_muscl
         WRITE(numout,*) '             MUSCL2 advection scheme        ln_trcadv_muscl2 = ', ln_trcadv_muscl2
      ENDIF

      !  Define the lateral tracer physics parameters
      ! =============================================

!!DB: I only need this to get ahtrc0 and trcrat
!!DB: ahtrc0 is needed #if NOT defined key_traldf_smag, 
!!    which likely would cause problems elsewhere 
!!DB: Note that above call to tra_adv_ctl and subsequent assignments
!!    means that most of this namelist is not needed and has been eliminated.
!!DB: Also, (eliminated param) ln_trcldf_diff controls whether or not to perform lateral diffusion
!!    so if you want NO lateral diffusion you must go to trctrp.F90
! Read Namelist namtrcldf : Lateral physics on tracers
      REWIND( numnat )
      READ  ( numnat, namtrcldf )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'DB: lateral passive tracer physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist namtrcldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '     laplacian operator                             ln_trcldf_lap   = ', l_trcldf_lap
         WRITE(numout,*) '     bilaplacian operator                           ln_trcldf_bilap = ', l_trcldf_bilap
         WRITE(numout,*) '     lateral eddy diffusivity                              ahtrc0   = ', ahtrc0
         WRITE(numout,*) '     ratio between passive and active tracer diffusion coef  trcrat = ', trcrat
      ENDIF

!!DB: REM, assigned above
      ! Parameter print
      ! ---------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'DB: vertical physics'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '          Namelist namtrczdf : set vertical diffusion parameters'
         WRITE(numout,*) '             time splitting / backward scheme ln_trczdf_exp = ', l_trczdf_exp
         WRITE(numout,*) '             number of time step               n_trczdf_exp = ', n_trczdf_exp
      ENDIF

!!DB: not accounted for yet
# if defined key_trcdmp
      !!DB: I have not covered this key yet so be careful
      ! Read Namelist namtdp : passive tracers damping term
      ! --------------------
      REWIND ( numnat )
      READ   ( numnat, namtrcdmp )
      IF( lzoom ) THEN
         nmldmptr = 0           ! restoring to climatology at closed north or south boundaries
      ENDIF

      ! Parameter control and print
      ! ---------------------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'newtonian damping'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '          Namelist namtrcdmp : set damping parameter'
         WRITE(numout,*)
         WRITE(numout,*) '             tracers damping option         ndmptr   = ', ndmptr
         WRITE(numout,*) '             create a damping.coeff file    ndmpftr  = ', ndmpftr
         WRITE(numout,*) '             mixed layer damping option     nmldmptr = ', nmldmptr, '(zoom: forced to 0)'
         WRITE(numout,*) '             surface time scale (days)      sdmptr   = ', sdmptr
         WRITE(numout,*) '             bottom time scale (days)       bdmptr   = ', bdmptr
         WRITE(numout,*) '             depth of transition (meters)   hdmptr   = ', hdmptr
         WRITE(numout,*)
      ENDIF

#endif


!!DB: create ncdf file 
      ptracer_fname = trim(cexper)//'_ptracer_01.nc'
      call ncdf_create_ptracer_file(status)


   END SUBROUTINE ptracer_setup


!!DB: Initialize model here
!!NB: For example of using pointers to trn(:,:,:,X) see lib_bgcm_01.F90
   SUBROUTINE ptracer_ini
      !!---------------------------------------------------------------------
      !!              
      !! ** Purpose : initialization of PTRACER_01 variables
      !!

     USE lbclnk

     INTEGER ::                   & 
          ji ,jj ,jk ,jn, jl        ! dummy loop indices  
     !!---------------------------------------------------------------------

     integer :: i1,i2,j1,j2

!!DB: -- compute wmask (
     wmask(:,:,1) = 0.0;   wmask(:,:,jpk+1) = 0.0
     do ji = 1, jpi
        do jj = 1, jpj
           do jk = 2, jpk
              wmask(ji,jj,jk) = min(tmask(ji,jj,jk-1),tmask(ji,jj,jk))
           enddo
        enddo
     enddo

     
     !! 1. initialization of passive tracer fields
     !! -------------------------------------------
     trn(:,:,:,:)=0.0
     tra(:,:,:,:)=0.0

!!DB: Arbitrary initialization for testing purposes -- HARDWIRED
!start one near boundary to test OBC code
     ji = 10
     jn = 1
     trn(mi0(ji):mi1(ji),mj0(10):mj1(30),1,jn) = 10.0
!     do jn = 1, jptra 
     do jn = 2, jptra 
        ji = 80 + 10*(jn-1)
        ji = min(ji,jpidta-10)   !!keep inside domain 
        trn(mi0(ji):mi1(ji),mj0(60):mj1(80),1,jn) = 10.0
     enddo


!!Look for and read restart_ptracer.nc 
!!NB: defaults to first record 
    i1 = 0
    open(125,file='restart_ptracer.nc',STATUS='OLD',ERR=666)
    i1 = 1
666 close(125)
    if(i1 > 0) then
       if(lwp) write(numout2,*)'PTRACER: Using restart_ptracer.nc'
       jj = 1
       CALL ncdf_read('restart_ptracer.nc', 'trn', trn, jptra, -jj, ji) 
    endif

!!To be safe
    do jn = 1, jptra
       CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )   
    enddo

!!DB: extract boundary values 
!!The below retained for possible future use
    call extract_boundary_vals
!! set before field
    e_trb(:,:,:) = e_trn(:,:,:)
    w_trb(:,:,:) = w_trn(:,:,:)
    s_trb(:,:,:) = s_trn(:,:,:)
    n_trb(:,:,:) = n_trn(:,:,:)


!! before field :
!! -------------
     trb(:,:,:,:) = trn(:,:,:,:)

     if( lwp ) then
        write(numout,*) ' '
        write(numout,*) 'PTRACER_01 initialisation done '
        write(numout,*) ' '
     endif

!write a file with the initial fields in it
      ptracer_fname = trim(cexper)//'_ptracer_init.nc'
      call ncdf_create_ptracer_file(ji)
      jj = 1
      CALL ncdf_write(ptracer_fname, 'time_counter', REAL(1), jj, ji)
      CALL ncdf_write(ptracer_fname, 'trn', trn, jptra, -jj, ji)
      CALL ncdf_write(ptracer_fname, 'ndastp',REAL(ndastp), jj, ji)
      CALL ncdf_write(ptracer_fname, 'model_time_step',REAL(1), jj, ji)
      CALL ncdf_write(ptracer_fname, 'model_time',REAL(1), jj, ji)

!Restore correct filename
      ptracer_fname = trim(cexper)//'_ptracer_01.nc'
      
 END SUBROUTINE ptracer_ini


!!DB: extract boundary values from trn and store them in 
!!the ?_trn boundary arrays
SUBROUTINE extract_boundary_vals
  USE obc_oce
#  include "obc_vectopt_loop_substitute.h90"
  IMPLICIT NONE
  INTEGER :: ji,jj,jk,jn

  do ji = fs_nie0+1, fs_nie1+1 ! isolate processor
     e_trn(:,:,:) = trn(ji,:,:,:)
  enddo
  do ji = niw0, niw1 ! isolate processor
     w_trn(:,:,:) = trn(ji,:,:,:)
  enddo
  do jj = njs0, njs1
     s_trn(:,:,:) = trn(:,jj,:,:)
  enddo
!!DB: no north OB
!  do jj = njn0, njn1
!     n_trn(:,:,:) = trn(:,jj,:,:)
!  enddo
END SUBROUTINE extract_boundary_vals

!!DB: Put the ?_trn boundary arrays values into tra at the boundary indices
!!Routine retained for possible future use
!!For example, if using data along Open Boundaries then a version of this
!!might be called
SUBROUTINE assign_boundary_vals
  USE obc_oce
#  include "obc_vectopt_loop_substitute.h90"
  IMPLICIT NONE
  INTEGER :: ji,jj,jk,jn

  do ji = fs_nie0+1, fs_nie1+1 !as in obctra
     tra(ji,:,:,:) = e_trn(:,:,:)
  enddo
  do ji = niw0, niw1 ! isolate processor
     tra(ji,:,:,:) = w_trn(:,:,:) 
  enddo
  do jj = njs0, njs1
     tra(:,jj,:,:) = s_trn(:,:,:) 
  enddo
!!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,:) = n_trn(:,:,:) 
!  enddo

END SUBROUTINE assign_boundary_vals

!!DB: 2009.04.20 -- copy values 1 cell in to OB positions
!!use "after" arrays
!!This routine is effectively the OBC code. It is called in trcnxt.F90
SUBROUTINE update_boundary_vals
  USE obc_oce
#  include "obc_vectopt_loop_substitute.h90"
  IMPLICIT NONE
  INTEGER :: ji,jj,jk,jn

!!DB: assign to outer(inner) values as well
  do ji = fs_nie0+1, fs_nie1+1 !as in obctra
     tra(ji,:,:,:) = tra(ji-1,:,:,:) 
  enddo
  do ji = niw0, niw1 ! isolate processor
     tra(ji,:,:,:) = tra(ji+1,:,:,:) 
  enddo
  do jj = njs0, njs1
     tra(:,jj,:,:) = tra(:,jj+1,:,:)
  enddo
!!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,:) = tra(:,jj-1,:,:)
!  enddo

END SUBROUTINE update_boundary_vals

!!DB NB: deleted the below but might want to include an example of assigning OBC values
!!Update N OBC vals and assign to tra()
!!To be called in trcnxt()
!  SUBROUTINE bgcm_N_obc ( kt )
!    USE obc_oce         !open boundary condition variables
!    USE lbclnk
!#  include "obc_vectopt_loop_substitute.h90"
    
    !!------------------------------------------------------------------------------
    !! * Arguments
!    INTEGER, INTENT( in ) ::   kt
!! the rest deleted
!! ...
!!
!!  END SUBROUTINE bgcm_N_obc

!!DB: 2008.08.28 ...
!!Create output file for PTRACER_01
!!The main variable is the trn(:,:,:,jptra) array which is what is integrated
!!by the tracer module. 
 SUBROUTINE ncdf_create_ptracer_file(status)

   USE lib_ncdf
   USE trc

    IMPLICIT NONE
    ! Subroutine argument declarations

    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
               varid,   &  ! ID of netCDF variable to be written to
               nfstat,  &  ! netCDF library call return status
               mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:5) :: dimids
!!DB -- unsure how to dimension varids, so choose a large number 
    INTEGER,DIMENSION(1:50) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since    
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    CHARACTER(LEN=20) :: op_type, varname
    INTEGER :: int_opp, &              ! Operation interval
               int_wri                 ! Write interval
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    INTEGER :: jf, varnum
!!DB NEW
    REAL(wp), DIMENSION(jpk) :: param_val = 0.0

    
    ! Initializations
    op_type = 'instantaneous'       !!default
!    op_type = 'ave(x)'       !!requires routine to do this  
    op_type = TRIM(op_type)


    status = NCDF_NOERR
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0

!!DB
    int_opp = nwrite * rdt
    int_wri = nwrite * rdt

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating output file:', ptracer_fname
       CALL FLUSH(100)
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ptracer_file - Creating file:', ptracer_fname
          CALL FLUSH(100)
       END IF
       ! Create the file
       nfstat = nf90_create(ptracer_fname, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ptracer_file - Defining dimensions in file:', ptracer_fname
          CALL FLUSH(100)
       END IF

       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'jptra', jptra, dimids(2))  !!CN: Changed dim creation order
       nfstat = nf90_def_dim(ncid, 'z', jpk, dimids(3))          !!ptracer model
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(4))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(5))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ptracer_file - Defining variables in file:', ptracer_fname
          CALL FLUSH(100)
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(5), dimids(4) /), &
            varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(5), dimids(4) /), &
            varids(2))
       nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &
            (/ dimids(3) /), &
            varids(3))
       nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
            (/ dimids(1) /), &
            varids(4))

       nfstat = nf90_def_var(ncid, 'trn', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(5))
       varnum = 5

!AD/DB: add new time-related variables
       varnum = varnum + 1
       ! ndate (ndastp)
       nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
            (/ dimids(1) /), &
            varids(varnum))
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'time step date in year/month/day aammjj')

       varnum = varnum + 1
       ! ndate (model_time)
       nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
            (/ dimids(1) /), &
            varids(varnum))
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', &
            'time step date (when output is written) in year/month/day aammjj (decimal day)')
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula1', 'nyear  =   model_time / 10000')       
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula2', & 
            'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula3', & 
            'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           

       varnum = varnum + 1
       ! kt 
       nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
            (/ dimids(1) /),  varids(varnum))

       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ptracer_file - Writing attributes in file:', ptracer_fname
          CALL FLUSH(100)
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', minval(glamt))
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', maxval(glamt))
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', minval(gphit))
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', maxval(gphit))
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! deptht
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', minval(gdept))      !for ptracer model
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', maxval(gdept))      !for ptracer model
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical Tracer levels')

!!DB: arbitrary attributes 
       ! trn 
       nfstat = nf90_put_att(ncid, varids(5), 'units', 'XXX')
       nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'PTRACER_01 tracers')
       nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'trn')
       nfstat = nf90_put_att(ncid, varids(5), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TNZYX')
       nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(5), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       ! time_counter
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(ptracer_fname))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'MISC', 'DB-created file')


       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_ptracer_file - Closing file:', ptracer_fname
          CALL FLUSH(100)
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF
    END IF


!!Write grid info to file
    CALL ncdf_write(ptracer_fname, 'nav_lat', gphit, -1, status)
    CALL ncdf_write(ptracer_fname, 'nav_lon', glamt, -1, status)
    CALL ncdf_write(ptracer_fname, 'deptht', gdept, status)

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_ptracer_file


END MODULE lib_ptracer_01

#endif 
