!!DB: 2010.08 -- Skeleton Module for : key_BGCM_02 = IML BioGeoChemical Model
!! Contains most of the needed routines.

!!NB: The intention is to include the basic routines that would
!!need to be accessed (mostly by bgcm_02_model) for the new model to run.
!!Some of these routines would have to be modified for the actual variables 
!!names, values, etc. 
!!DB has assumed 13 variables need to be time-stepped 
!!(====> trn(:,:,:,jptra) is the "now" array and jptra = 13 in par_trc_trp.F90),
!!and he has used pointers to them (v1, v2, ..., v13)


!!Note the sequence: (local) initrc.F90 USES this module if 
!! key_BGCM_02 is defined. A consequence is that the model will use 
!! the subroutines in this module instead of other versions. A key example
!! is ini_trc() which initializes the tracer module. 

!!NB: I may not have added the north open boundary to all relevant routines so
!!     BEWARE !!!! 
!!(Also, note that must have north boundary for these code frags to work)

!!See lib_bgcm_01.F90 for more comments

#if defined key_passivetrc

MODULE lib_bgcm_02
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
   !USE rivers, only:iriv_a, jriv_a, iriv_b, jriv_b, nriver, nrivermax, iside

   !! From dtatem for no3tem, DL
!DL   USE oce             ! ocean dynamics and tracers
!DL   USE dom_oce         ! ocean space and time domain
!DL   USE in_out_manager  ! I/O manager
!DL   USE daymod          ! calendar

   IMPLICIT NONE
   PRIVATE

!!DB: 2010.08 : kept wmask; computed in bgcm_ini
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk+1) ::  wmask

!!DB
!!Shared module variables, default values. NB: they are reset below
!!NB: These are required due to the overall structure of OPA. Eventually the
!!code could be simplified and much of the below would not be necessary
   !! Advection
   LOGICAL, PUBLIC ::   &
      ln_trcadv_cen2   = .TRUE. ,  & !!: 2nd order centered scheme flag
      ln_trcadv_tvd    = .FALSE. ,  &  !: TVD scheme flag
      ln_trcadv_muscl  = .FALSE. ,  &  !: MUSCL scheme flag
      ln_trcadv_muscl2 = .FALSE. ,  &  !: MUSCL2 scheme flag
      ln_trcadv_smolar = .FALSE.        !: Smolarkiewicz scheme flag (DB: must be .FALSE.)

   !! Lateral diffusion
   LOGICAL , PUBLIC ::              & !!: ** lateral mixing namelist (nam_trcldf) **
      ln_trcldf_diff  = .FALSE. ,   &  !: flag of perform or not the lateral diff.
      ln_trcldf_lap   = .TRUE.  ,   &  !: laplacian operator
      ln_trcldf_bilap = .FALSE. ,   &  !: bilaplacian operator
      ln_trcldf_level = .FALSE. ,   &  !: iso-level direction
      ln_trcldf_hor   = .FALSE. ,   &  !: horizontal (geopotential) direction
      ln_trcldf_iso   = .TRUE.         !: iso-neutral direction

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
   CHARACTER (len=80), PUBLIC :: bgcm_fname

#if defined DIAG_NPZD_GROWTH || defined DIAG_NPZD_flux 
   CHARACTER (len=80), PUBLIC :: bgcm_diag_fname   !NL#12
#endif
  
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
      n_dta             !: nitrate data at given time-step

   !! * dta_no3 related variables
   CHARACTER (len=45) ::   &
      cl_tdata
   INTEGER ::   &
      nlecte =  0,   &  ! switch for the first read
      ntem1      ,   &  ! first record used
      ntem2             ! second record used
   REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
      no3dta            ! nitrate data at two consecutive times, obc
!NL june2013 to reduce input of monthly obc to once a month
   INTEGER, DIMENSION(3) :: ntobc2_old_NO3(3)

#if defined OXYGEN
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
      o_dta             !: oxygen data at given time-step
   REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
      oxydta            ! oxygen data at two consecutive times
#endif

#if defined key_carbon
   REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &  !:
      c_dta             !: carbon data at given time-step
   REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   &
      dicdta            ! carbon data at two consecutive times
#endif
!DL ******************************************************************


   !! * Accessibility
   PUBLIC ini_trc, bgcm_setup, bgcm_ini, extract_boundary_vals, &
        assign_boundary_vals, update_boundary_vals, bgcm_obc_02, river_no3, &
        rst_bgcm_write
#if defined (NPZD_INT_PROD)
!    PUBLIC   output_intprod
#endif

CONTAINS

!!DB: routine ini_trc replaces regular ini_trc when key_BGCM_02 is defined
!!
  SUBROUTINE ini_trc
      !!---------------------------------------------------------------------
      !!
      !!                       ROUTINE ini_trc
      !!                     ******************
      !!
      !!  PURPOSE :  initialize the BGCM_02 model
      !!  DB: Replaces regular ini_trc when key_BGCM_02 is defined
      !!---------------------------------------------------------------------


      !! 0.b PRINT the number of tracer
      !! ------------------------------

    IF(lwp) WRITE(numout,*) ' '
    IF(lwp) WRITE(numout,*) ' *** number of passive tracer jptra = ',jptra
    IF(lwp) WRITE(numout,*) ' '
    
    call bgcm_setup 
    call bgcm_ini

! =====  (NL#2 : the reading of restart is done in  bgcm_ini (no need to re-do)
!!DB If restarting (not thought about yet) assume that the following call
!!DB will overwrite anything wrong that bgcm_ini might have done
!NL    if( lrsttr ) THEN
!         CALL bgcm_rst       !!...TO DO ...
!NL    endif

  END SUBROUTINE ini_trc


!!DB: Replace trc_lec, trc_trp_lec, trc_trp_ctl ... with 1 routine
!!One of the things I do in this routine is force the same adv/diff schemes
!!as are used by the TS module.
!!Initialized netcdf output file as well
!!2010.08 -- NB: Because I force the use of certain schemes (e.g. adv/diff), much of the
!!           below may not be necessary. Cleaning this up is left for a future exercise.

   SUBROUTINE bgcm_setup
     USE traadv_ctl
     USE ldftra_oce
     USE zdf_oce, ONLY : n_zdfexp
     USE lib_ncdf

      !!                     ******************
      !! local declarations
      !! ==================

      INTEGER ::  ji, status
      CHARACTER (len=32) :: clname

!!DB: These are all that are needed for BGCM
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
         WRITE(numout,*) ' ROUTINE BGCM_02_SETUP'
         WRITE(numout,*) ' **************'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist for passive tracer'
         WRITE(numout,*) ' ***************************'
         WRITE(numout,*) ' '
      ENDIF

      numnat=80
!!DB: NB
      clname='namelist.BGCM_02'
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
         WRITE(numout,*) 'Namelist namtrcldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) 'laplacian operator                             ln_trcldf_lap   = ', l_trcldf_lap
         WRITE(numout,*) 'bilaplacian operator                           ln_trcldf_bilap = ', l_trcldf_bilap
         WRITE(numout,*) 'lateral eddy diffusivity                              ahtrc0   = ', ahtrc0
         WRITE(numout,*) 'ratio between passive and active tracer diffusion coef  trcrat = ', trcrat
      ENDIF

!!DB: REM, assigned above
      ! Parameter print
      ! ---------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'DB: vertical physics'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) 'Namelist namtrczdf : set vertical diffusion parameters'
         WRITE(numout,*) ' time splitting / backward scheme ln_trczdf_exp = ', l_trczdf_exp
         WRITE(numout,*) ' number of time step               n_trczdf_exp = ', n_trczdf_exp
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
      bgcm_fname = trim(cexper)//'_bgcm_02.nc'
#if defined DIAG_NPZD_GROWTH || defined DIAG_NPZD_flux || defined DIAG_NPZD_PROD
     bgcm_diag_fname = trim(cexper)//'_bgcm_diag.nc'  !NL#12
#endif

      call ncdf_create_bgcm_file(status)


   END SUBROUTINE bgcm_setup


!!DB: Initialize model here
!!NB: Note example of using pointers to trn(:,:,:,?) (see bgcm_02_pointers.h90)
   SUBROUTINE bgcm_ini
      !!---------------------------------------------------------------------
      !!              
      !! ** Purpose : initialization of BGCM_02 variables
      !!

     USE lbclnk
     INTEGER :: ji ,jj ,jk ,jn, jl        ! dummy loop indices  
     integer :: i1,i2,j1,j2
!!DB: The below must go here, i.e. after last local declaration and before first assignment statement
!!    The file contains the assignment of variable names as pointers to trn() locations
!AD: inside lifemaker now #   include "bgcm_02_pointers.h90"

!!DB: -- compute wmask (may be needed)
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
     trb(:,:,:,:)=0.0

!DL jan 2012 Call new routine to initialise 3D fields
     call dta_no3(nit000)     ! use 3D nitrate field
!     call river_no3(nit000)     ! use quebec_nitrate
     trn(:,:,:,3)=n_dta(:,:,:)
     if (lwp) print*, 'called dta_no3 in bgcm_ini'

!DL november 2012, addition of oxygen
#if defined OXYGEN
     call dta_oxy(nit000)     ! use 3D oxygen field
     trn(:,:,:,9)=o_dta(:,:,:)
     if (lwp) print*, 'called dta_oxy in bgcm_ini'
# endif

!DL april 2013, addition of carbon
#if defined key_carbon
     call dta_dic(nit000)     ! use 3D DIC field
     trn(:,:,:,10)=c_dta(:,:,:)
     if (lwp) print*, 'called dta_dic in bgcm_ini'
# endif

!DL !!Look for and read restart_bgcm.nc 
!DL !!NB: defaults to first record 
!DL     i1 = 0
!DL     open(125,file='restart_bgcm.nc',STATUS='OLD',ERR=666)
!DL     i1 = 1
!DL 666 close(125)
!DL     if(i1 > 0) then
!DL        if(lwp) write(numout2,*)'BGCM: Using restart_bgcm.nc'
!DL        jj = 1
!DL        CALL ncdf_read('restart_bgcm.nc', 'trn', trn, jptra, -jj, ji) 
!DL    endif
! =========================!! NL#2
      IF( lrsttr ) THEN                    ! Restart from a file
        CALL rst_bgcm_read                       ! Read the restart file
      ENDIF
! =========================!! NL#2

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
        write(numout,*) 'BGCM_02 initialisation done '
        write(numout,*) ' '
     endif

!write a file with the initial fields in it
      bgcm_fname = trim(cexper)//'_bgcm_init.nc'
      call ncdf_create_bgcm_file(ji)
      jj = 1
      CALL ncdf_write(bgcm_fname, 'time_counter', REAL(1), jj, ji)
      CALL ncdf_write(bgcm_fname, 'trn', trn, jptra, -jj, ji)
      CALL ncdf_write(bgcm_fname, 'ndastp',REAL(ndastp), jj, ji)
      CALL ncdf_write(bgcm_fname, 'model_time_step',REAL(1), jj, ji)
      CALL ncdf_write(bgcm_fname, 'model_time',REAL(1), jj, ji)

!Restore correct filename
      bgcm_fname = trim(cexper)//'_bgcm_02.nc'
      
 END SUBROUTINE bgcm_ini


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
SUBROUTINE update_boundary_vals(jn)
  USE obc_oce
#  include "obc_vectopt_loop_substitute.h90"
  IMPLICIT NONE
  !! * Arguments
  INTEGER, INTENT( in ) ::   jn
  !DL, insertion of jn indice, may 2012 
  INTEGER :: ji,jj,jk
  REAL, DIMENSION(jpi,jpk) :: tra_smsk   !NL#6
  REAL, DIMENSION(jpj,jpk) :: tra_emsk,tra_wmsk   !NL#6

!   tra_smsk=0.

! NL#6 build the mask of the inflow (1 = in, 0= out), to put the value only there
      ! south border :  the method is copied from obcta.F90 => subroutine obc_tra_south
      !            the phase velocity taken is the phase velocity (here v)
      tra_smsk = sign( 1.,vfos)  ! now is ±1
!      tra_smsk = 1.   !NL#8
      tra_smsk= tsmsk * (tra_smsk + abs(tra_smsk))/2         ! now is 0 or 1
      ! east border :  the method is copied from obcta.F90 => subroutine obc_tra_east
      !            the phase velocity taken is the phase velocity (here -u)
      tra_emsk = sign(1.,-ufoe )  ! now is ±1
!      tra_emsk = 1.   !NL#8
      tra_emsk= temsk * (tra_emsk + abs(tra_emsk))/2     ! now is 0 or 1
      ! west border :  the method is copied from obcta.F90 => subroutine obc_tra_east
      !                 the phase velocity taken is the phase velocity (here u)
      tra_wmsk = sign(1., ufow )  ! now is ±1
!      tra_wmsk = 1.   !NL#8
      tra_wmsk= twmsk * (tra_wmsk + abs(tra_wmsk))/2        ! now is 0 or 1
   !NL#6
!DL need to put value from inside boundary if current are going out
!******************************************************************

  if (jn==1) then !diatoms
! a fixed value will tmporarily be assigned at the boundary for diatoms
! need to assign a climatology at the boundaries and applied only for mixed layer
      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
           DO jk = 1, jpkm1
              DO jj = 1, jpj
!        tra(ji,:,:,jn) = 0.16 
             tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &
                               tra_emsk(jj,jk) * 0.4
          enddo
        enddo 
      enddo
      do ji = fs_niw0, fs_niw1 ! isolate processor
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &
                                tra_wmsk(jj,jk) * 0.16 
!         tra(ji,:,:,jn) = 0.16
          enddo
        enddo 
      enddo
      do jj = fs_njs0, fs_njs1
             DO jk = 1, jpkm1
                DO ji = 1, jpi
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &
                                  tra_smsk(ji,jk) * 0.16 
        ! tra(:,jj,:,jn) = 0.16 
           enddo
        enddo
      enddo
!!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,jn) = tra(:,jj-1,:,jn)
!  enddo
  else if (jn==2) then !flagellates
! a fixed value will tmporarily be assigned at the boundary for diatoms
! need to assign a climatology at the boundaries and applied only for mixed layer
!      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
!        tra(ji,:,:,jn) = 0.1  
!      enddo
!      do ji = fs_niw0, fs_niw1 ! isolate processor
!         tra(ji,:,:,jn) = 0.1 
!      enddo
!      do jj = fs_njs0, fs_njs1
!         tra(:,jj,:,jn) = 0.1 
!      enddo
      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
           DO jk = 1, jpkm1
              DO jj = 1, jpj
             tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &
                               tra_emsk(jj,jk) * 0.15
          enddo
        enddo
      enddo
      do ji = fs_niw0, fs_niw1 ! isolate processor
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &
                                tra_wmsk(jj,jk) * 0.1
          enddo
        enddo
      enddo
      do jj = fs_njs0, fs_njs1
             DO jk = 1, jpkm1
                DO ji = 1, jpi
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &
                                  tra_smsk(ji,jk) * 0.1
           enddo
        enddo
      enddo

!!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,jn) = tra(:,jj-1,:,jn)
!  enddo
  else if (jn==4) then !ammonium
! a fixed value will tmporarily be assigned at the boundary for ammonium
! need to assign a climatology at the boundaries and applied only for mixed layer
  !    do ji = fs_nie0+1, fs_nie1+1 !as in obctra
  !      tra(ji,:,:,jn) = 0.5  
  !    enddo
  !    do ji = fs_niw0, fs_niw1 ! isolate processor
  !       tra(ji,:,:,jn) = 0.5 
  !    enddo
  !    do jj = fs_njs0, fs_njs1
  !       tra(:,jj,:,jn) = 0.2 
  !    enddo
      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
           DO jk = 1, jpkm1
              DO jj = 1, jpj
             tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &
                               tra_emsk(jj,jk) * 0.6
          enddo
        enddo
      enddo
      do ji = fs_niw0, fs_niw1 ! isolate processor
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &
                                tra_wmsk(jj,jk) * 0.1
          enddo
        enddo
      enddo
      do jj = fs_njs0, fs_njs1
             DO jk = 1, jpkm1
                DO ji = 1, jpi
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &
                                  tra_smsk(ji,jk) * 0.1
           enddo
        enddo
      enddo
!
!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,jn) = tra(:,jj-1,:,jn)
!  enddo
  else if (jn==5) then !large zoo
! a fixed value will tmporarily be assigned at the boundary for ammonium
! need to assign a climatology at the boundaries and applied only for mixed layer
      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
           DO jk = 1, jpkm1
              DO jj = 1, jpj
             tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &
                               tra_emsk(jj,jk) * 0.6
          enddo
        enddo
      enddo
      do ji = fs_niw0, fs_niw1 ! isolate processor
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &
                                tra_wmsk(jj,jk) * 0.2
          enddo
        enddo
      enddo
      do jj = fs_njs0, fs_njs1
             DO jk = 1, jpkm1
                DO ji = 1, jpi
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &
                                  tra_smsk(ji,jk) * 0.2
           enddo
        enddo
      enddo
!
  else if (jn==6) then ! small zoo
! a fixed value will tmporarily be assigned at the boundary for ammonium
! need to assign a climatology at the boundaries and applied only for mixed layer
      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
           DO jk = 1, jpkm1
              DO jj = 1, jpj
             tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &
                               tra_emsk(jj,jk) * 0.3
          enddo
        enddo
      enddo
      do ji = fs_niw0, fs_niw1 ! isolate processor
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &
                                tra_wmsk(jj,jk) * 0.2
          enddo
        enddo
      enddo
      do jj = fs_njs0, fs_njs1
             DO jk = 1, jpkm1
                DO ji = 1, jpi
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &
                                  tra_smsk(ji,jk) * 0.2
           enddo
        enddo
      enddo
!
!!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,jn) = tra(:,jj-1,:,jn)
!  enddo
  else if (jn==7) then !POM
! a fixed value will tmporarily be assigned at the boundary for ammonium
! need to assign a climatology at the boundaries and applied only for mixed layer
!      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
!        tra(ji,:,:,jn) = 0.002
!      enddo
!      do ji = fs_niw0, fs_niw1 ! isolate processor
!         tra(ji,:,:,jn) = 0.002
!      enddo
!      do jj = fs_njs0, fs_njs1
!         tra(:,jj,:,jn) = 0.002
!      enddo
      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
           DO jk = 1, jpkm1
              DO jj = 1, jpj
             tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &
                               tra_emsk(jj,jk) * 0.4
          enddo
        enddo
      enddo
      do ji = fs_niw0, fs_niw1 ! isolate processor
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &
                                tra_wmsk(jj,jk) * 0.04
          enddo
        enddo
      enddo
      do jj = fs_njs0, fs_njs1
             DO jk = 1, jpkm1
                DO ji = 1, jpi
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &
                                  tra_smsk(ji,jk) * 0.04
           enddo
        enddo
      enddo
!
!!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,jn) = tra(:,jj-1,:,jn)
!  enddo
  else if (jn==8) then !DOM
! a fixed value will tmporarily be assigned at the boundary for ammonium
! need to assign a climatology at the boundaries and applied only for mixed layer
!      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
!        tra(ji,:,:,jn) = 0.01
!      enddo
!      do ji = fs_niw0, fs_niw1 ! isolate processor
!         tra(ji,:,:,jn) = 0.01
!      enddo
!      do jj = fs_njs0, fs_njs1
!         tra(:,jj,:,jn) = 0.01
!      enddo
      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
           DO jk = 1, jpkm1
              DO jj = 1, jpj
             tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &
                               tra_emsk(jj,jk) * 0.8
          enddo
        enddo
      enddo
      do ji = fs_niw0, fs_niw1 ! isolate processor
           DO jk = 1, jpkm1
              DO jj = 1, jpj
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &
                                tra_wmsk(jj,jk) * 0.04
          enddo
        enddo
      enddo
      do jj = fs_njs0, fs_njs1
             DO jk = 1, jpkm1
                DO ji = 1, jpi
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &
                                  tra_smsk(ji,jk) * 0.04
           enddo
        enddo
      enddo
!
!!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,jn) = tra(:,jj-1,:,jn)
!  enddo
!    else
!!DB: assign to outer(inner) values as well
!      do ji = fs_nie0+1, fs_nie1+1 !as in obctra
!         tra(ji,:,:,jn) = tra(ji-1,:,:,jn) 
!      enddo
!      do ji = fs_niw0, fs_niw1 ! isolate processor
!         tra(ji,:,:,jn) = tra(ji+1,:,:,jn) 
!      enddo
!      do jj = fs_njs0, fs_njs1
!         tra(:,jj,:,jn) = tra(:,jj+1,:,jn)
!      enddo
!!DB: no north OB
!  do jj = njn0, njn1
!     tra(:,jj,:,jn) = tra(:,jj-1,:,jn)
!  enddo
    endif 
END SUBROUTINE update_boundary_vals

!!DB: example of assigning OBC values to a specific variable
!!To be called in trcnxt()
SUBROUTINE bgcm_obc_02( kt, jn )
!     USE par_oce
!     USE obc_par
     USE obc_oce
     USE lib_ncdf
#  include "obc_vectopt_loop_substitute.h90"
    
 !!------------------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT( in ) ::   kt, jn
    !! * Local declaration
    INTEGER ::   ji, jj, jk, jeast1, jeast2, ij, ii
    INTEGER ::   i, j, k
    INTEGER ::   ntobc1, ntobc2            !DL
    INTEGER ::   itimo, iman, imois        !DL
    INTEGER ::   i15                       !DL

    REAL(wp), DIMENSION(jpi,jpj) :: array1  !:
    REAL(wp) :: area1
    REAL(wp) :: zxy                         !DL

!DL Additions---------------------------
!: nitrate climatology at the boundaries

    REAL(wp), DIMENSION(jpj,jpk,5) ::  no3foe, no3fow
    REAL(wp), DIMENSION(jpi,jpk,5) ::  no3fos
!, no3fon              

    !: array used for interpolating monthly data on the east boundary
    REAL(wp), DIMENSION(1:jpjef,jpk,5,jptobc), SAVE ::  no3edta 

    !: array used for interpolating monthly data on the west boundary
    REAL(wp), DIMENSION(1:jpjwf,jpk,5,jptobc), SAVE ::  no3wdta 

    !: array used for interpolating monthly data on the south boundary
    REAL(wp), DIMENSION(1:jpisf,jpk,5,jptobc), SAVE ::  no3sdta
! no3ndta

#if defined OXYGEN
    REAL(wp), DIMENSION(jpj,jpk,5) ::  oxyfoe, oxyfow
    REAL(wp), DIMENSION(jpi,jpk,5) ::  oxyfos
!, oxyfon              

    !: array used for interpolating monthly data on the east boundary
    REAL(wp), DIMENSION(1:jpjef,jpk,5,jptobc) ::  oxyedta 

    !: array used for interpolating monthly data on the west boundary
    REAL(wp), DIMENSION(1:jpjwf,jpk,5,jptobc) ::  oxywdta 

    !: array used for interpolating monthly data on the south boundary
    REAL(wp), DIMENSION(1:jpisf,jpk,5,jptobc) ::  oxysdta
! oxyndta
#endif

#if defined key_carbon
    REAL(wp), DIMENSION(jpj,jpk,5) ::  dicfoe, dicfow
    REAL(wp), DIMENSION(jpi,jpk,5) ::  dicfos
!, dicfon              

    !: array used for interpolating monthly data on the east boundary
    REAL(wp), DIMENSION(1:jpjef,jpk,5,jptobc) ::  dicedta 

    !: array used for interpolating monthly data on the west boundary
    REAL(wp), DIMENSION(1:jpjwf,jpk,5,jptobc) ::  dicwdta 

    !: array used for interpolating monthly data on the south boundary
    REAL(wp), DIMENSION(1:jpisf,jpk,5,jptobc) ::  dicsdta
! dicndta
#endif

    REAL, DIMENSION(jpjef-jpjed+1,jpk,5) :: buf_3De       ! Needed for lib_ncdf calls
    REAL, DIMENSION(jpjwf-jpjwd+1,jpk,5) :: buf_3Dw       ! Needed for lib_ncdf calls
    REAL, DIMENSION(jpisf-jpisd+1,jpk,5) :: buf_3Ds       ! Needed for lib_ncdf calls
    INTEGER :: f_stat
    REAL, DIMENSION(jpi,jpk) :: tra_smsk   !NL#6
    REAL, DIMENSION(jpj,jpk) :: tra_emsk,tra_wmsk   !NL#6

!DL --------------------------------------
 
!DL    if (mod(kt-nit000,180) .eq. 0) then
! NL#6 build the mask of the inflow (1 = in, 0= out), to put the value only there
      ! south border :  the method is copied from obcta.F90 => subroutine obc_tra_south
      !            the phase velocity taken is the phase velocity (here v)
      tra_smsk = sign( 1.,vfos )  ! now is ±1
      tra_smsk= tsmsk * (tra_smsk + abs(tra_smsk))/2         ! now is 0 or 1
      ! east border :  the method is copied from obcta.F90 => subroutine obc_tra_east
      !            the phase velocity taken is the phase velocity (here -u)
      tra_emsk = sign(1.,-ufoe )  ! now is ±1
      tra_emsk= temsk * (tra_emsk + abs(tra_emsk))/2     ! now is 0 or 1
      ! west border :  the method is copied from obcta.F90 => subroutine obc_tra_east
      !                 the phase velocity taken is the phase velocity (here u)
      tra_wmsk = sign(1., ufow )  ! now is ±1
      tra_wmsk= twmsk * (tra_wmsk + abs(tra_wmsk))/2        ! now is 0 or 1
   !NL#6

! 2.  Initialize the time we are at.
     !     Does this every time the routine is called,
     !     excepted when nobc_dta = 0
     !---------------------------------------------------------------------
           iman  = 12
           i15   = nday / 16
           imois = nmonth + i15 - 1
           IF( imois == 0 )   imois = iman
           itimo = imois

     ! For linear interpolation of BCs to current time step
     ! ----------------------------------------------------
           zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.


! 2.1 Read two records in the file 
     ! ---------------------------------------------

        ! Calendar computation
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

       if (jn==3) then
!DL *** Reading boundary conditions for no3, as for T and S in obcdta.F90  
!DL *** Eastern boundary conditions

       IF( lp_obc_east )   THEN

!!!!!!!!!!!!!!!!!!!!!!!!
!! NL : read the 2 month at the first time, and read only the next when the month change.
          IF( kt == nit000 ) THEN  ! first step, read both month (past and post)
!DL test, augmention NO3 fronti'ere est (*1.1)
!************************************* 
            CALL ncdf_read_global('obceast_no3.nc', 'vonitrate', buf_3De, -ntobc1, f_stat)
            no3edta(jpjed:jpjef,:,1:5,1)=buf_3De*1.15
            CALL ncdf_read_global('obceast_no3.nc', 'vonitrate', buf_3De, -ntobc2, f_stat)
            no3edta(jpjed:jpjef,:,1:5,2)=buf_3De*1.15
            ntobc2_old_NO3(1)=ntobc2

            if( lwp ) print*,'First reading of the obceast_no3.nc'

          ELSEIF (ntobc2_old_NO3(1).ne.ntobc2) THEN ! when the month change, load the next month
             if( lwp ) print*,'testing Reading of the obceast_no3.nc for ',ntobc2
             if( lwp ) print*,'ntobc2_old_NO3, ntobc2 ',ntobc2_old_NO3(1),ntobc2, kt

             no3edta(jpjed:jpjef,:,1:5,1)=no3edta(jpjed:jpjef,:,1:5,2)
             CALL ncdf_read_global('obceast_no3.nc', 'vonitrate', buf_3De, - ntobc2, f_stat)
             no3edta(jpjed:jpjef,:,1:5,2)=buf_3De*1.15
             ntobc2_old_NO3(1)=ntobc2

             if( lwp ) print*,'Reading of the obceast_no3.nc for ',ntobc2

          ENDIF


           DO jk = 1, jpkm1
              DO jj = nje0p1, nje1m1
                 ij = jj -1 + njmpp
                 no3foe(jj,jk,:) =  ( zxy * no3edta(ij,jk,:,2) + &
                    &           (1.-zxy) * no3edta(ij,jk,:,1) ) * temsk5(jj,jk,:)
              END DO
           END DO

!DL What about the anomaly that is added to T and S, see obc_dta

        e_trn(:,:,jn)=no3foe(:,:,1) ! as in obctra.F90

        DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 1, jpj
!                 tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-temsk(jj,jk)) + &
!                               temsk(jj,jk) * no3foe(jj,jk,1)
                 tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &! NL#6
                               tra_emsk(jj,jk) * no3foe(jj,jk,1)! NL#6

              END DO
           END DO
        END DO

!DL    do ji = fs_nie0+1, fs_nie1+1 !as in obctra
!DL        tra(ji,:,:,jn) = e_trn(:,:,jn)
!DL    enddo

       ENDIF

!DL *** Western boundary conditions

       IF( lp_obc_west )   THEN

!! NL : read the 2 month at the first time, and read only the next when the month change.

          IF( kt == nit000 ) THEN  ! first step only

              CALL ncdf_read_global('obcwest_no3.nc', 'vonitrate', buf_3Dw, -ntobc1, f_stat)
              no3wdta(jpjwd:jpjwf,:,1:5,1)=buf_3Dw
              CALL ncdf_read_global('obcwest_no3.nc', 'vonitrate', buf_3Dw, -ntobc2, f_stat)
              no3wdta(jpjwd:jpjwf,:,1:5,2)=buf_3Dw
              ntobc2_old_NO3(2)=ntobc2

              if( lwp ) print*,'First reading of the obcwest_no3.nc'

           ELSEIF (ntobc2_old_NO3(2).ne.ntobc2) THEN ! when the month change, load the next month
              no3wdta(jpjwd:jpjwf,:,1:5,1)=no3wdta(jpjwd:jpjwf,:,1:5,2)
              CALL ncdf_read_global('obcwest_no3.nc', 'vonitrate', buf_3Dw, - ntobc2, f_stat)
              no3wdta(jpjwd:jpjwf,:,1:5,2)=buf_3Dw
              ntobc2_old_NO3(2)=ntobc2

              IF( lwp ) print*,'Reading of the obcwest_no3.nc for ',ntobc2

          ENDIF

           DO jk = 1, jpkm1
              DO jj = njw0p1, njw1m1
                 ij = jj -1 + njmpp
                 no3fow(jj,jk,:) =  ( zxy * no3wdta(ij,jk,:,2) + &
                    &           (1.-zxy) * no3wdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
              END DO
           END DO

         w_trn(:,:,jn)=no3fow(:,:,1) ! as in obctra.F90

        DO ji = fs_niw0, fs_niw1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 1, jpj
!DL                  tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-twmsk(jj,jk)) + &
!DL                                twmsk(jj,jk) * no3fow(jj,jk,1)
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &! NL#6
                                tra_wmsk(jj,jk) * no3fow(jj,jk,1)! NL#6
               END DO
            END DO
         END DO
       ENDIF

       IF( lp_obc_south )   THEN

!! NL : read the 12 month at the first time, take data otherwise
           IF( kt == nit000 ) THEN  ! first step only

              CALL ncdf_read_global('obcsouth_no3.nc', 'vonitrate', buf_3Ds, -ntobc1, f_stat)
              no3sdta(jpisd:jpisf,:,1:5,1)=buf_3Ds
              CALL ncdf_read_global('obcsouth_no3.nc', 'vonitrate', buf_3Ds, -ntobc2, f_stat)
              no3sdta(jpisd:jpisf,:,1:5,2)=buf_3Ds
              ntobc2_old_NO3(3)=ntobc2

              if( lwp ) print*,'First reading of the obcsouth_no3.nc'

           ELSEIF (ntobc2_old_NO3(3).ne.ntobc2) THEN ! when the month change, load the next month
              no3sdta(jpisd:jpisf,:,1:5,1)=no3sdta(jpisd:jpisf,:,1:5,2)
              CALL ncdf_read_global('obcsouth_no3.nc', 'vonitrate', buf_3Ds, -ntobc2, f_stat)
              no3sdta(jpisd:jpisf,:,1:5,2)=buf_3Ds
              ntobc2_old_NO3(3)=ntobc2

              if( lwp ) print*,'Reading of the obcsouth_no3 for ',ntobc2

           ENDIF

           DO jk = 1, jpkm1
             DO ji = nis0p1, nis1m1
                ii = ji -1 + nimpp
                no3fos(ji,jk,:) = ( zxy * no3sdta(ii,jk,:,2) + &
                   &          (1.-zxy) * no3sdta(ii,jk,:,1) ) * tsmsk5(ji,jk,:)
             END DO
           END DO

          s_trn(:,:,jn)=no3fos(:,:,1) ! as in obctra.F90

          DO jj = fs_njs0, fs_njs1  ! Vector opt.
             DO jk = 1, jpkm1
                DO ji = 1, jpi
!DL                    tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-tsmsk(ji,jk)) + &
!DL                                  tsmsk(ji,jk) * no3fos(ji,jk,1)
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &! NL#6
                                  tra_smsk(ji,jk) * no3fos(ji,jk,1)! NL#6

                 END DO
              END DO
           END DO

        ENDIF
!
!       IF( lp_obc_north )   THEN  !No northern boundary
!         CALL obc_dta_gv ('x','vonitrate',jpinf-jpind+1,ntobc1,pdta_4D=no3ndta(jpind:jpinf,:,1:5,1),fn='obcnorth_no3.nc')
!         CALL obc_dta_gv ('x','vonitrate',jpinf-jpind+1,ntobc2,pdta_4D=no3ndta(jpind:jpinf,:,1:5,2),fn='obcnorth_no3.nc')
!       ENDIF
!     !!DB: north is closed so do nothing.
!        IF( lp_obc_north )   THEN
!           DO jk = 1, jpkm1
!              DO ji = nin0p1, nin1m1
!                 ii = ji -1 + nimpp
!                 no3fon(ji,jk,:) =  ( zxy * no3ndta(ii,jk,:,2) + &
!                    &           (1.-zxy) * no3ndta(ii,jk,:,1) ) * tnmsk5(ji,jk,:)
!              END DO
!           END DO
!        ENDIF
!
    
!DL      endif !modulo

#if defined OXYGEN
       elseif (jn==9) then

!DL *** Reading boundary conditions for oxy  
!DL *** Eastern boundary conditions

       IF( lp_obc_east )   THEN
           
           CALL ncdf_read_global('obceast_OXY.nc', 'vo_OXY', buf_3De, -ntobc1, f_stat)
           oxyedta(jpjed:jpjef,:,1:5,1)=buf_3De*0.92
           CALL ncdf_read_global('obceast_OXY.nc', 'vo_OXY', buf_3De, -ntobc2, f_stat)
           oxyedta(jpjed:jpjef,:,1:5,2)=buf_3De*0.92

           DO jk = 1, jpkm1
              DO jj = nje0p1, nje1m1
                 ij = jj -1 + njmpp
                 oxyfoe(jj,jk,:) =  ( zxy * oxyedta(ij,jk,:,2) + &
                    &           (1.-zxy) * oxyedta(ij,jk,:,1) ) * temsk5(jj,jk,:)
              END DO
           END DO

        e_trn(:,:,jn)=oxyfoe(:,:,1) ! as in obctra.F90

        DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 1, jpj
!                 tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-temsk(jj,jk)) + &
!                               temsk(jj,jk) * oxyfoe(jj,jk,1)
                 tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &! NL#6
                               tra_emsk(jj,jk) * oxyfoe(jj,jk,1)! NL#6

              END DO
           END DO
        END DO

       ENDIF

!DL *** Western boundary conditions

       IF( lp_obc_west )   THEN

           CALL ncdf_read_global('obcwest_OXY.nc', 'vo_OXY', buf_3Dw, -ntobc1, f_stat)
           oxywdta(jpjwd:jpjwf,:,1:5,1)=buf_3Dw
           CALL ncdf_read_global('obcwest_OXY.nc', 'vo_OXY', buf_3Dw, -ntobc2, f_stat)
           oxywdta(jpjwd:jpjwf,:,1:5,2)=buf_3Dw


           DO jk = 1, jpkm1
              DO jj = njw0p1, njw1m1
                 ij = jj -1 + njmpp
                 oxyfow(jj,jk,:) =  ( zxy * oxywdta(ij,jk,:,2) + &
                    &           (1.-zxy) * oxywdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
              END DO
           END DO

         w_trn(:,:,jn)=oxyfow(:,:,1) ! as in obctra.F90

        DO ji = fs_niw0, fs_niw1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 1, jpj
!                  tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-twmsk(jj,jk)) + &
!                                twmsk(jj,jk) * oxyfow(jj,jk,1)
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &! NL#6
                                tra_wmsk(jj,jk) * oxyfow(jj,jk,1)! NL#6

               END DO
            END DO
         END DO

       ENDIF

!DL Southern boundary

       IF( lp_obc_south )   THEN

           CALL ncdf_read_global('obcsouth_OXY.nc', 'vo_OXY', buf_3Ds, -ntobc1, f_stat)
           oxysdta(jpisd:jpisf,:,1:5,1)=buf_3Ds
           CALL ncdf_read_global('obcsouth_OXY.nc', 'vo_OXY', buf_3Ds, -ntobc2, f_stat)
           oxysdta(jpisd:jpisf,:,1:5,2)=buf_3Ds

           DO jk = 1, jpkm1
             DO ji = nis0p1, nis1m1
                ii = ji -1 + nimpp
                oxyfos(ji,jk,:) = ( zxy * oxysdta(ii,jk,:,2) + &
                   &          (1.-zxy) * oxysdta(ii,jk,:,1) ) * tsmsk5(ji,jk,:)
             END DO
           END DO

          s_trn(:,:,jn)=oxyfos(:,:,1) ! as in obctra.F90

          DO jj = fs_njs0, fs_njs1  ! Vector opt.
             DO jk = 1, jpkm1
                DO ji = 1, jpi
!                    tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-tsmsk(ji,jk)) + &
!                                  tsmsk(ji,jk) * oxyfos(ji,jk,1)
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &! NL#6
                                  tra_smsk(ji,jk) * oxyfos(ji,jk,1)! NL#6
                 END DO
              END DO
           END DO

        ENDIF

# endif
!c
#if defined key_carbon
       elseif (jn==10) then

!DL *** Reading boundary conditions for dic  
!DL *** Eastern boundary conditions

       IF( lp_obc_east )   THEN
           
           CALL ncdf_read_global('obceast_DIC.nc', 'vo_DIC', buf_3De, -ntobc1, f_stat)
           dicedta(jpjed:jpjef,:,1:5,1)=buf_3De
           CALL ncdf_read_global('obceast_DIC.nc', 'vo_DIC', buf_3De, -ntobc2, f_stat)
           dicedta(jpjed:jpjef,:,1:5,2)=buf_3De

           DO jk = 1, jpkm1
              DO jj = nje0p1, nje1m1
                 ij = jj -1 + njmpp
                 dicfoe(jj,jk,:) =  ( zxy * dicedta(ij,jk,:,2) + &
                    &           (1.-zxy) * dicedta(ij,jk,:,1) ) * temsk5(jj,jk,:)
              END DO
           END DO

        e_trn(:,:,jn)=dicfoe(:,:,1) ! as in obctra.F90

        DO ji = fs_nie0+1, fs_nie1+1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 1, jpj
!                 tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-temsk(jj,jk)) + &
!                               temsk(jj,jk) * dicfoe(jj,jk,1)
                 tra(ji,jj,jk,jn)= tra(ji-1,jj,jk,jn) * (1.-tra_emsk(jj,jk)) + &! NL#6
                               tra_emsk(jj,jk) * dicfoe(jj,jk,1)! NL#6

              END DO
           END DO
        END DO

       ENDIF

!DL *** Western boundary conditions

       IF( lp_obc_west )   THEN

           CALL ncdf_read_global('obcwest_DIC.nc', 'vo_DIC', buf_3Dw, -ntobc1, f_stat)
           dicwdta(jpjwd:jpjwf,:,1:5,1)=buf_3Dw
           CALL ncdf_read_global('obcwest_DIC.nc', 'vo_DIC', buf_3Dw, -ntobc2, f_stat)
           dicwdta(jpjwd:jpjwf,:,1:5,2)=buf_3Dw


           DO jk = 1, jpkm1
              DO jj = njw0p1, njw1m1
                 ij = jj -1 + njmpp
                 dicfow(jj,jk,:) =  ( zxy * dicwdta(ij,jk,:,2) + &
                    &           (1.-zxy) * dicwdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
              END DO
           END DO

         w_trn(:,:,jn)=dicfow(:,:,1) ! as in obctra.F90

        DO ji = fs_niw0, fs_niw1 ! Vector opt.
           DO jk = 1, jpkm1
              DO jj = 1, jpj
!                  tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-twmsk(jj,jk)) + &
!                                twmsk(jj,jk) * dicfow(jj,jk,1)
                  tra(ji,jj,jk,jn)= tra(ji+1,jj,jk,jn) * (1.-tra_wmsk(jj,jk)) + &! NL#6
                                tra_wmsk(jj,jk) * dicfow(jj,jk,1)! NL#6
               END DO
            END DO
         END DO

       ENDIF

!DL Southern boundary

       IF( lp_obc_south )   THEN

           CALL ncdf_read_global('obcsouth_DIC.nc', 'vo_DIC', buf_3Ds, -ntobc1, f_stat)
           dicsdta(jpisd:jpisf,:,1:5,1)=buf_3Ds
           CALL ncdf_read_global('obcsouth_DIC.nc', 'vo_DIC', buf_3Ds, -ntobc2, f_stat)
           dicsdta(jpisd:jpisf,:,1:5,2)=buf_3Ds

           DO jk = 1, jpkm1
             DO ji = nis0p1, nis1m1
                ii = ji -1 + nimpp
                dicfos(ji,jk,:) = ( zxy * dicsdta(ii,jk,:,2) + &
                   &          (1.-zxy) * dicsdta(ii,jk,:,1) ) * tsmsk5(ji,jk,:)
             END DO
           END DO

          s_trn(:,:,jn)=dicfos(:,:,1) ! as in obctra.F90

          DO jj = fs_njs0, fs_njs1  ! Vector opt.
             DO jk = 1, jpkm1
                DO ji = 1, jpi
!                    tra(ji,jj,jk,jn)= tra(ji,jj,jk,jn) * (1.-tsmsk(ji,jk)) + &
!                                  tsmsk(ji,jk) * dicfos(ji,jk,1)
                    tra(ji,jj,jk,jn)= tra(ji,jj+1,jk,jn) * (1.-tra_smsk(ji,jk)) + &! NL#6
                                  tra_smsk(ji,jk) * dicfos(ji,jk,1)! NL#6

                 END DO
              END DO
           END DO

        ENDIF

# endif
  endif
 END SUBROUTINE bgcm_obc_02


!!DB: 2008.08.28 ...
!!Create output file for BGCM_02
!!The main variable is the trn(:,:,:,jptra) array which is what is integrated
!!by the tracer module. 
 SUBROUTINE ncdf_create_bgcm_file(status)

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
    INTEGER,DIMENSION(1:6) :: dimids
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
       WRITE(100,*) 'NCDF DEBUG: Creating output file:', bgcm_fname
       CALL FLUSH(100)
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Creating file:', bgcm_fname
          CALL FLUSH(100)
       END IF
       ! Create the file
       nfstat = nf90_create(bgcm_fname, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF

       ! global Attribut
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(bgcm_fname))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'MISC', 'DB-created file')
     
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Defining dimensions in file:', bgcm_fname
          CALL FLUSH(100)
       END IF

       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'jptra', jptra, dimids(2))  !!CN: Changed dim creation order
       nfstat = nf90_def_dim(ncid, 'z', jpk, dimids(3))          !!bgcm model
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(4))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(5))
       nfstat = nf90_def_dim(ncid, 'hour', 12, dimids(6))
       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Defining variables in file:', bgcm_fname
          CALL FLUSH(100)
       END IF

       ! longitude =================
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(5), dimids(4) /), &
            varids(1))
       ! attributs
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', minval(glamt))
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', maxval(glamt))
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! Latitude =================
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(5), dimids(4) /), &
            varids(2))
       ! attributs
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', minval(gphit))
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', maxval(gphit))
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! deptht =================
       nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &
            (/ dimids(3) /), &
            varids(3))
       ! Attributs
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', minval(gdept))      !for bgcm model
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', maxval(gdept))      !for bgcm model
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical Tracer levels')

       ! time_counter =================
       nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
            (/ dimids(1) /), &
            varids(4))
       ! Attribut
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! trn  =================
       nfstat = nf90_def_var(ncid, 'trn', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(2), dimids(1) /), &
            varids(5))
       !!DB: arbitrary attributes 
       nfstat = nf90_put_att(ncid, varids(5), 'units', 'XXX')
       nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'BGCM_02 tracers')
       nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'trn')
       nfstat = nf90_put_att(ncid, varids(5), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TNZYX')
       nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(5), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       varnum = 5

!       ! transx_no3  =================    ! NL#3
!       varnum = varnum+1     ! NL  add +1 the number of variable for the parz
!       nfstat = nf90_def_var(ncid, 'transx_no3', nf90_float, &
!            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
!            varids(varnum))   ! NL
!       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmol/m2/s')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'transx_no3')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Transport of nitrate through the cell in x')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter deptht nav_lat nav_lon')
!
!       ! transy_no3  =================    ! NL#3
!       varnum = varnum+1     ! NL  add +1 the number of variable for the parz
!       nfstat = nf90_def_var(ncid, 'transy_no3', nf90_float, &
!            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
!            varids(varnum))   ! NL
!       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmol/m2/s')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'transy_no3')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Transport of nitrate through the cell in y')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter deptht nav_lat nav_lon')
!
!       ! transz_no3  =================    ! NL#3
!       varnum = varnum+1     ! NL  add +1 the number of variable for the parz
!       nfstat = nf90_def_var(ncid, 'transz_no3', nf90_float, &
!            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
!            varids(varnum))   ! NL
!       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmol/m2/s')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'transz_no3')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Transport of nitrate through the cell in z')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter deptht nav_lat nav_lon')
!
!       ! parz  =================    ! NL
!       varnum = varnum+1     ! NL  add +1 the number of variables
!       nfstat = nf90_def_var(ncid, 'parz', nf90_float, &
!            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
!            varids(varnum))   ! NL
!       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'W/m2')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'parz')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'PAR at depth z (3D)')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter deptht nav_lat nav_lon')
!
       ! par_surf  =================    ! NL
       varnum = varnum+1     ! NL  add +1 the number of variables
       nfstat = nf90_def_var(ncid, 'par_surf', nf90_float, &
            (/ dimids(5), dimids(4), dimids(6), dimids(1) /), &
!            (/ dimids(5), dimids(4), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'W/m2')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'par_surf')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'PAR at the surface')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

!       ! kcdom  =================    ! NL
!       varnum = varnum+1     ! NL  add +1 the number of variables
!       nfstat = nf90_def_var(ncid, 'kcdom', nf90_float, &
!            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
!            varids(varnum))   ! NL
!       nfstat = nf90_put_att(ncid, varids(varnum), 'units', '/m')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'kcdom')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'attenuation coef. from CDOM')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter deptht nav_lat nav_lon')

!       ! euphotique layer  =================    ! NL#4
!       varnum = varnum+1     !  add +1 the number of variables
!       nfstat = nf90_def_var(ncid, 'ze01', nf90_float, &
!            (/ dimids(5), dimids(4), dimids(1) /), &
!            varids(varnum))   !
!       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'm')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
!       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'ze01')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Euphotic layer thickness')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
!       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TYX')
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
!       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')
!
!
       ! Production nouvelle des Diatoms  =================    ! NL#5
#if defined (NPZD_INT_PROD)
       varnum = varnum+1     !  add +1 the number of variable for pnwd
       nfstat = nf90_def_var(ncid, 'pnwd', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   !
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mgC/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pnwd')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'New Diatoms production')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

       ! Regenerated Diatoms production  =================    ! NL#5
       varnum = varnum+1     !  add +1 the number of variable for pregd
       nfstat = nf90_def_var(ncid, 'pregd', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   !
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mgC/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pregd')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Regenerated Diatoms production')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

       ! New Flagellate production  =================    ! NL#5
       varnum = varnum+1     !  add +1 to the number of variables
       nfstat = nf90_def_var(ncid, 'pnwf', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   !
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mgC/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pnwf')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'New Flagellate production')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

       ! Regenerated flagellated production  =================    ! NL#5
       varnum = varnum+1     !  add +1 the number of variable
       nfstat = nf90_def_var(ncid, 'pregf', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   !
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mgC/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pregf')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Regenerated Flagellated production')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

       ! Mesozoo production from Diatoms  =================    ! NL#5
       varnum = varnum+1     !  add +1 to the number of variable 
       nfstat = nf90_def_var(ncid, 'pmesd', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   !
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mgC/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pmesd')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Secondary Production => mesozoo diatoms')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

       ! Mesozoo production from microzoo  =================    ! NL#5
       varnum = varnum+1     !  add +1 to the number of variable 
       nfstat = nf90_def_var(ncid, 'pmesm', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   !
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mgC/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pmesm')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Secondary Production => mesozoo on microzoo')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

       ! Microzoo production  =================    ! NL#5
       varnum = varnum+1     !  add +1 to the number of variable 
       nfstat = nf90_def_var(ncid, 'pmicf', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   !
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mgC/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pmicf')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Secondary Production => micro')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')
#endif


#if defined (key_carbon)
       ! pH  =================    ! NL
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH 
       nfstat = nf90_def_var(ncid, 'pH', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       ! pH attributs , DL   + NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'none')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pH')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       ! satcalc  =================    ! NL
       varnum = varnum+1     ! NL  add +1 the number of variable for the satcalc
       nfstat = nf90_def_var(ncid, 'satca', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'none')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'satcalc')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'saturation state of calcite')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       ! satarg  =================    ! NL
       varnum = varnum+1     ! NL  add +1 the number of variable for the satarg
       nfstat = nf90_def_var(ncid, 'satar', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'none')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'satarg')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'saturation state of aragonite')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       ! Surface CO2 flux  =================    
       varnum = varnum+1     !  add +1 the number of variables
       nfstat = nf90_def_var(ncid, 'flxCO2', nf90_float, &
            (/ dimids(5), dimids(4), dimids(1) /), &
            varids(varnum))   !
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmol/m2/d')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'flxCO2')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Surface CO2 flux')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', op_type)
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter nav_lat nav_lon')

#endif

!AD/DB: add new time-related variables
       ! ndate (ndastp) =================
       varnum = varnum + 1
       nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
            (/ dimids(1) /), &
            varids(varnum))
       ! attribut
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'time step date in year/month/day aammjj')

       ! ndate (model_time) =================
       varnum = varnum + 1
       nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
            (/ dimids(1) /), &
            varids(varnum))
       ! Attribut
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', &
            'time step date (when output is written) in year/month/day aammjj (decimal day)')
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula1', 'nyear  =   model_time / 10000')       
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula2', & 
            'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
       nfstat = nf90_put_att(ncid, varids(varnum), 'formula3', & 
            'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           

       ! kt =================
       varnum = varnum + 1
       nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
            (/ dimids(1) /),  varids(varnum))

       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Writing attributes in file:', bgcm_fname
          CALL FLUSH(100)
       END IF

       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Closing file:', bgcm_fname
          CALL FLUSH(100)
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF
!    END IF
!NL#12 & NL#10 Create the file for the diagnostic output
#if defined DIAG_NPZD_GROWTH || defined DIAG_NPZD_flux || defined DIAG_NPZD_PROD
       ! Create the file
       nfstat = nf90_create(bgcm_diag_fname, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF

       ! global Attribut
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(bgcm_diag_fname))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'MISC', 'DB-created file')

       ! Define dimensions
       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'jptra', jptra, dimids(2))  !!CN: Changed dim creation order
       nfstat = nf90_def_dim(ncid, 'z', jpk, dimids(3))          !!bgcm model
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(4))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(5))

!  =================
       ! Define variables
       varnum =4 ! reset de varnum and put 4 for the first variable (lon,lat,depth, time)
       ! longitude
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
            (/ dimids(5), dimids(4) /), &
            varids(1))
       ! attributs
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', minval(glamt))
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', maxval(glamt))
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
       nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')

       ! Latitude =================
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
            (/ dimids(5), dimids(4) /), &
            varids(2))
       ! attributs
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', minval(gphit))
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', maxval(gphit))
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
       nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')

       ! deptht =================
       nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &
            (/ dimids(3) /), &
            varids(3))
       ! Attributs
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
       nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', minval(gdept))      !for bgcm model
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', maxval(gdept))      !for bgcm model
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical Tracer levels')

       ! time_counter =================
       nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
            (/ dimids(1) /), &
            varids(4))
       ! Attribut
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
#if defined DIAG_NPZD_GROWTH
       ! limitating factor  =================    ! NL#10
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'tcd', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'none')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'tcd')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Diatoms C doubling time')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_mean')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'tsd', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'none')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'tcd')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Diatoms N doubling time')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_mean')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'tcf', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'none')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'tcf')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flagellates C doubling time')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_mean')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'tsf', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'none')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'tcf')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flagellates N doubling time')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_mean')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       ! ^^^^ limitating factor  =================    ! NL#10

       varnum = varnum+1     ! sediment trap
       nfstat = nf90_def_var(ncid, 'trap', nf90_float, &
            (/ dimids(5), dimids(4), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'XX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'tcf')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Sediment trap')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_accumulation')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! bottom sediment
       nfstat = nf90_def_var(ncid, 'sediment', nf90_float, &
            (/ dimids(5), dimids(4), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'XX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'tcf')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Bottom sediment')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_accumulation')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
#endif
#if defined DIAG_NPZD_flux
       ! flux in modis  =================    ! NL#12
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'dmizdon', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'dmizdon')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Dead miz flow to don')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #G')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'dspdon', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'dspdon')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Dead flagelate (sp) to don')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #E')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'fpondon', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'fpondon')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Fragmentation of PON to DON')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #K')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'redonnh4', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'redonnh4')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Remineralization from DON to NH4')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #F')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'ppsp', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'ppsp')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : primary production of the flagelates (NO3 + NH$ to sp)')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #A2')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'pplp', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'pplp')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : primary production of the diat (NO3 + NH4 to lp)')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #A1')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'nh4no3', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'nh4no3')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : nitrification (NH4 to NO3)')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #L')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'dlppon', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'dlppon')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Dead diatom cells (LP to PON)')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #B')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'grlpmez', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'grlpmez')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Grazing by mesozoo (MEZ) from LP')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #C')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'grspmiz', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'grspmiz')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Grazing by microzoo (MIZ) from SP')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #D')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'grmezmiz', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'grmezmiz')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Grazing by mesozoo (MEZ) from MIZ')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #I')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'dmezpon', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'dmezpon')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Dead mesozoo (MEZ) to PON')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #J')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'grponmiz', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'grponmiz')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Grazing by microzoo (MIZ) from PON')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #M')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'emeznh4', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'emeznh4')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Excretion by MEZ to NH4')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #N')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'unmicdom', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'unmidom')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : un-assimilited part of microzoo (MIZ) routed to DON')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #M + flow #D (lost), add to flow #G')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'unmezpon', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'unmezpon')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Un-assimilated part of mesozoo (MEZ) diet routed to PON')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #C + flow #I  (lost), add to flow #J')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'unmiznh4', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'unmiznh4')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Un-assimilited part of microzoo diet routed to NH4')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #M + flow #D (lost), add to flow #H')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'sedlp', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'sedlp')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Sedimentation of diat (LP) to below layer')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #Z2')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       varnum = varnum+1     ! NL  add +1 the number of variable for the pH
       nfstat = nf90_def_var(ncid, 'sedpon', nf90_float, &
            (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
            varids(varnum))   ! NL
       nfstat = nf90_put_att(ncid, varids(varnum), 'units', 'mmolN/m3/day')
       nfstat = nf90_put_att(ncid, varids(varnum), 'missing_value', 1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_min', 1.000000E20 )
       nfstat = nf90_put_att(ncid, varids(varnum), 'valid_max', -1.000000E20)
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'sedpon')
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Flow Nit. : Sedimentation of PON to below layer')
       nfstat = nf90_put_att(ncid, varids(varnum), 'flow_code', 'flow #Z1')
       nfstat = nf90_put_att(ncid, varids(varnum), 'online_operation', 'Daily_total')
       nfstat = nf90_put_att(ncid, varids(varnum), 'axis', 'TZYX')
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_operation', float(int_opp))
       nfstat = nf90_put_att(ncid, varids(varnum), 'interval_write', float(int_wri))
       nfstat = nf90_put_att(ncid, varids(varnum), 'associate', 'time_counter jptra deptht nav_lat nav_lon')

       ! ^^^^ flux of NPZD =================    ! NL#12
#endif
       ! Close file
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF

#endif
    END IF  ! nproc ==0


!!Write grid info to file
    CALL ncdf_write(bgcm_fname, 'nav_lat', gphit, -1, status)
    CALL ncdf_write(bgcm_fname, 'nav_lon', glamt, -1, status)
    CALL ncdf_write(bgcm_fname, 'deptht', gdept, status)

!NL#12 & NL#10 Write grid info to file
#if defined DIAG_NPZD_GROWTH || defined DIAG_NPZD_flux || defined DIAG_NPZD_PROD
    CALL ncdf_write(bgcm_diag_fname, 'nav_lat', gphit, -1, status)
    CALL ncdf_write(bgcm_diag_fname, 'nav_lon', glamt, -1, status)
    CALL ncdf_write(bgcm_diag_fname, 'deptht', gdept, status)
#endif

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_bgcm_file

!******************************************
! New subroutine for nitrate initialisation
! Diane Lavoie, Januray 2012
! ******************************************

!DL Check is needed in step.F90
!
!   !! * Routine accessibility
!   PUBLIC dta_tem   ! called by step.F90 and inidta.F90

  SUBROUTINE dta_no3( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_no3  ***
      !!                    
      !! ** Purpose :   Reads monthly nitrate data 
      !!
      !! ** Method  :   Read on unit numtdt the interpolated nitrate 
      !!      onto the model grid.
      !!      Data begin at january. 
      !!      The value is centered at the middle of month.
      !!      In the opa model, kt=1 agree with january 1.
      !!      At each time step, a linear interpolation is applied between 
      !!      two monthly values.
      !!      Read on unit numtdt
      !!
      !! ** Action  :   define n_dta array at time-step kt
      !!
      !!----------------------------------------------------------------------
      !! * Modules used
!!DB
      USE lib_ncdf

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step

      !! * Local declarations
      INTEGER, PARAMETER ::   &
         jpmois = 12                    ! number of month
      INTEGER ::   ji, jj, jl           ! dummy loop indicies
      INTEGER ::   &
         imois, iman, itime, ik ,    &  ! temporary integers
         i15, ipi, ipj, ipk             !    "          "
#  if defined key_tradmp
      INTEGER ::   &
         il0, il1, ii0, ii1, ij0, ij1   ! temporary integers
# endif

      INTEGER, DIMENSION(jpmois) ::   istep
      REAL(wp) ::   zxy, zl, zdate0
      REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
      REAL(wp), DIMENSION(jpk) ::   zlev
!!DB
      INTEGER :: len, status 
      !!----------------------------------------------------------------------

      ! 0. Initialization
      ! -----------------

      iman  = jpmois
      i15   = nday / 16
      imois = nmonth + i15 - 1
      IF( imois == 0 )   imois = iman

      itime = jpmois
      ipi = jpiglo
      ipj = jpjglo
      ipk = jpk

      ! 1. First call kt=nit000
      ! -----------------------

      IF( kt == nit000 .AND. nlecte == 0 ) THEN
         ntem1 = 0
         cl_tdata = 'data_1m_no3_nomask.nc'
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' dta_no3 : Monthly Nitrate fields'
         IF(lwp) WRITE(numout,*) ' ~~~~~~'
         IF(lwp) WRITE(numout,*) '             NetCDF File'
         IF(lwp) WRITE(numout,*) cl_tdata
         IF(lwp) WRITE(numout,*)
         
#if defined key_agrif
         if ( .NOT. Agrif_Root() ) then
            cl_tdata = TRIM(Agrif_CFixed())//'_'//TRIM(cl_tdata)
         endif
#endif         
!!DB: 
         call ncdf_get_dim_size(cl_tdata, 'Time', itime, status)
         call ncdf_get_dim_size(cl_tdata, 'x', ipi, status)
         call ncdf_get_dim_size(cl_tdata, 'y', ipj, status)
         call ncdf_get_dim_size(cl_tdata, 'z', ipk, status)
         if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with dimensions'
               write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               write(numout,*) ' ipk ',ipk,' jpk ',jpk
            endif
            stop 'dta_no3'
         endif
!!DB: 
         if( itime /= jpmois ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with time coordinates'
               write(numout,*) ' itime ',itime,' jpmois ',jpmois
            endif
            stop 'dta_no3'
         endif

      ENDIF

      ! 2. Read monthly file
      ! -------------------

      IF( ( kt == nit000 .AND. nlecte == 0 ) .OR. imois /= ntem1 ) THEN
         nlecte = 1

!DL This variable needs to remain 0 if other variables are to be read
#if defined OXYGEN
         nlecte = 0
#endif
#if defined key_carbon
         nlecte = 0
#endif
         ! Calendar computation
         
         ntem1 = imois        ! first file record used 
         ntem2 = ntem1 + 1    ! last  file record used
         ntem1 = MOD( ntem1, iman )
         IF( ntem1 == 0 )   ntem1 = iman
         ntem2 = MOD( ntem2, iman )
         IF( ntem2 == 0 )   ntem2 = iman
         if(lwp) write(numout,*) 'dtatem reading records ',ntem1, ntem2
         
         ! Read monthly nitrate data
!!DB
         call ncdf_read(cl_tdata,'vonitrate',no3dta(:,:,:,1),-ntem1,status)
         call ncdf_read(cl_tdata,'vonitrate',no3dta(:,:,:,2),-ntem2,status)


         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' Monthly nitrate input OK'
         IF(lwp) WRITE(numout,*)

!DL Do I need to do something about the mask??
         !                                  ! Mask
         DO jl = 1, 2
            no3dta(:,:,:,jl) = no3dta(:,:,:,jl) * tmask(:,:,:)
            no3dta(:,:,jpk,jl) = 0.
!DL            IF( lk_zps ) THEN                ! z-coord. with partial steps
!DL               DO jj = 1, jpj                  ! interpolation of temperature at the last level
!DL                  DO ji = 1, jpi
!DL                     ik = mbathy(ji,jj) - 1
!DL                     IF( ik > 2 ) THEN
!DL                        zl = ( gdept(ik) - fsdept(ji,jj,ik) ) / ( gdept(ik) - gdept(ik-1) )
!DL                        no3dta(ji,jj,ik,jl) = (1.-zl) * no3dta(ji,jj,ik,jl) + zl * no3dta(ji,jj,ik-1,jl) 
!DL                     ENDIF
!DL                  END DO
!DL               END DO
!DL            ENDIF
         END DO

!         IF(lwp) THEN
!            WRITE(numout,*) ' temperature month ', ntem1, ntem2
!            WRITE(numout,*)
!            WRITE(numout,*) '  month = ', ntem1, '  level = 1'
!            CALL prihre( temdta(:,:,1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
!            WRITE(numout,*) '  month = ', ntem1, '  level = ', jpk/2
!            CALL prihre( temdta(:,:,jpk/2,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
!            WRITE(numout,*) '  month = ',ntem1,'  level = ', jpkm1
!            CALL prihre( temdta(:,:,jpkm1,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1., numout )
!         ENDIF
      ENDIF

      ! 2. At every time step compute nitrate data
      ! ----------------------------------------------
      zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!!DB: Note that for monthly data it is not necessary to do this every kt.
!!    So in the future use a code fragment like: 
!      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
!         t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)
!      endif
!DL Need to clarify this one
      n_dta(:,:,:) = (1.-zxy) * no3dta(:,:,:,1) + zxy * no3dta(:,:,:,2)


  END SUBROUTINE dta_no3

!DL #else
!DL    !!----------------------------------------------------------------------
!DL    !!   Default case                           NO 3D nitrate data field
!DL    !!----------------------------------------------------------------------
!DL    USE in_out_manager
!DL    LOGICAL , PUBLIC, PARAMETER ::   lk_dtatem = .FALSE.   !: temperature data flag
!DL CONTAINS
!DL    SUBROUTINE dta_no3( kt )        ! Empty routine
!DL       if(lwp) WRITE(numout,*) 'dta_no3: You should not have seen this print! error?', kt
!DL    END SUBROUTINE dta_no3

#if defined OXYGEN
 
 SUBROUTINE dta_oxy( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_oxy  ***
      !!                    
      !! ** Purpose :   Reads monthly oxygen data 
      !!
      !! ** Method  :   Read on unit numtdt the interpolated oxygen 
      !!      onto the model grid.
      !!      Data begin at january. 
      !!      The value is centered at the middle of month.
      !!      In the opa model, kt=1 agree with january 1.
      !!      At each time step, a linear interpolation is applied between 
      !!      two monthly values.
      !!      Read on unit numtdt
      !!
      !! ** Action  :   define o_dta array at time-step kt
      !!
      !!----------------------------------------------------------------------
      !! * Modules used
!!DB
      USE lib_ncdf

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step

      !! * Local declarations
      INTEGER, PARAMETER ::   &
         jpmois = 12                    ! number of month
      INTEGER ::   ji, jj, jl           ! dummy loop indicies
      INTEGER ::   &
         imois, iman, itime, ik ,    &  ! temporary integers
         i15, ipi, ipj, ipk             !    "          "
#  if defined key_tradmp
      INTEGER ::   &
         il0, il1, ii0, ii1, ij0, ij1   ! temporary integers
# endif

      INTEGER, DIMENSION(jpmois) ::   istep
      REAL(wp) ::   zxy, zl, zdate0
      REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
      REAL(wp), DIMENSION(jpk) ::   zlev
!!DB
      INTEGER :: len, status 
      !!----------------------------------------------------------------------

      ! 0. Initialization
      ! -----------------

      iman  = jpmois
      i15   = nday / 16
      imois = nmonth + i15 - 1
      IF( imois == 0 )   imois = iman

      itime = jpmois
      ipi = jpiglo
      ipj = jpjglo
      ipk = jpk

      ! 1. First call kt=nit000
      ! -----------------------

      IF( kt == nit000 .AND. nlecte == 0 ) THEN
         ntem1 = 0
         cl_tdata = 'data_1m_oxy_nomask.nc'
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' dta_oxy : Monthly oxygen fields'
         IF(lwp) WRITE(numout,*) ' ~~~~~~'
         IF(lwp) WRITE(numout,*) '             NetCDF File'
         IF(lwp) WRITE(numout,*) cl_tdata
         IF(lwp) WRITE(numout,*)
         
#if defined key_agrif
         if ( .NOT. Agrif_Root() ) then
            cl_tdata = TRIM(Agrif_CFixed())//'_'//TRIM(cl_tdata)
         endif
#endif         
!!DB: 
         call ncdf_get_dim_size(cl_tdata, 'Time', itime, status)
         call ncdf_get_dim_size(cl_tdata, 'x', ipi, status)
         call ncdf_get_dim_size(cl_tdata, 'y', ipj, status)
         call ncdf_get_dim_size(cl_tdata, 'z', ipk, status)
         if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with dimensions'
               write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               write(numout,*) ' ipk ',ipk,' jpk ',jpk
            endif
            stop 'dta_oxy'
         endif
!!DB: 
         if( itime /= jpmois ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with time coordinates'
               write(numout,*) ' itime ',itime,' jpmois ',jpmois
            endif
            stop 'dta_oxy'
         endif

      ENDIF

      ! 2. Read monthly file
      ! -------------------

      IF( ( kt == nit000 .AND. nlecte == 0 ) .OR. imois /= ntem1 ) THEN
         nlecte = 1

!NL This variable needs to remain 0 if other variables are to be read
#if defined key_carbon
         nlecte = 0
#endif
         ! Calendar computation
         
         ntem1 = imois        ! first file record used 
         ntem2 = ntem1 + 1    ! last  file record used
         ntem1 = MOD( ntem1, iman )
         IF( ntem1 == 0 )   ntem1 = iman
         ntem2 = MOD( ntem2, iman )
         IF( ntem2 == 0 )   ntem2 = iman
         if(lwp) write(numout,*) 'dta_oxy reading records ',ntem1, ntem2
         
         ! Read monthly nitrate data
!!DB
         call ncdf_read(cl_tdata,'vo_OXY',oxydta(:,:,:,1),-ntem1,status)
         call ncdf_read(cl_tdata,'vo_OXY',oxydta(:,:,:,2),-ntem2,status)


         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' Monthly oxygen input OK'
         IF(lwp) WRITE(numout,*)

!DL Do I need to do something about the mask??
         !                                  ! Mask
         DO jl = 1, 2
            oxydta(:,:,:,jl) = oxydta(:,:,:,jl) * tmask(:,:,:)
            oxydta(:,:,jpk,jl) = 0.
         END DO

      ENDIF

      ! 2. At every time step compute oxygen data
      ! ----------------------------------------------
      zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!!DB: Note that for monthly data it is not necessary to do this every kt.
!!    So in the future use a code fragment like: 
!      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
!         t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)
!      endif
!DL Need to clarify this one
      o_dta(:,:,:) = (1.-zxy) * oxydta(:,:,:,1) + zxy * oxydta(:,:,:,2)

  END SUBROUTINE dta_oxy
#endif 

#if defined key_carbon
 
 SUBROUTINE dta_dic( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dta_dic  ***
      !!                    
      !! ** Purpose :   Reads monthly dic data 
      !!
      !! ** Method  :   Read on unit numtdt the interpolated oxygen 
      !!      onto the model grid.
      !!      Data begin at january. 
      !!      The value is centered at the middle of month.
      !!      In the opa model, kt=1 agree with january 1.
      !!      At each time step, a linear interpolation is applied between 
      !!      two monthly values.
      !!      Read on unit numtdt
      !!
      !! ** Action  :   define c_dta array at time-step kt
      !!
      !!----------------------------------------------------------------------
      !! * Modules used
!!DB
      USE lib_ncdf

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt     ! ocean time-step

      !! * Local declarations
      INTEGER, PARAMETER ::   &
         jpmois = 12                    ! number of month
      INTEGER ::   ji, jj, jl           ! dummy loop indicies
      INTEGER ::   &
         imois, iman, itime, ik ,    &  ! temporary integers
         i15, ipi, ipj, ipk             !    "          "
#  if defined key_tradmp
      INTEGER ::   &
         il0, il1, ii0, ii1, ij0, ij1   ! temporary integers
# endif

      INTEGER, DIMENSION(jpmois) ::   istep
      REAL(wp) ::   zxy, zl, zdate0
      REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
      REAL(wp), DIMENSION(jpk) ::   zlev
!!DB
      INTEGER :: len, status 
      !!----------------------------------------------------------------------

      ! 0. Initialization
      ! -----------------
      iman  = jpmois
      i15   = nday / 16
      imois = nmonth + i15 - 1
      IF( imois == 0 )   imois = iman

      itime = jpmois
      ipi = jpiglo
      ipj = jpjglo
      ipk = jpk

      ! 1. First call kt=nit000
      ! -----------------------

      IF( kt == nit000 .AND. nlecte == 0 ) THEN
         ntem1 = 0
         cl_tdata = 'data_1m_DIC_nomask.nc'
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' dta_dic : Monthly DIC fields'
         IF(lwp) WRITE(numout,*) ' ~~~~~~'
         IF(lwp) WRITE(numout,*) '             NetCDF File'
         IF(lwp) WRITE(numout,*) cl_tdata
         IF(lwp) WRITE(numout,*)
         
#if defined key_agrif
         if ( .NOT. Agrif_Root() ) then
            cl_tdata = TRIM(Agrif_CFixed())//'_'//TRIM(cl_tdata)
         endif
#endif         
!!DB: 
         call ncdf_get_dim_size(cl_tdata, 'Time', itime, status)
         call ncdf_get_dim_size(cl_tdata, 'x', ipi, status)
         call ncdf_get_dim_size(cl_tdata, 'y', ipj, status)
         call ncdf_get_dim_size(cl_tdata, 'z', ipk, status)
         if( ipi /= jpidta .OR. ipj /= jpjdta .OR. ipk /= jpk ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with dimensions'
               write(numout,*) ' ipi ',ipi,' jpidta ',jpidta
               write(numout,*) ' ipj ',ipj,' jpjdta ',jpjdta
               write(numout,*) ' ipk ',ipk,' jpk ',jpk
            endif
            stop 'dta_dic'
         endif
!!DB: 
         if( itime /= jpmois ) then
            if(lwp) then
               write(numout,*)
               write(numout,*) 'problem with time coordinates'
               write(numout,*) ' itime ',itime,' jpmois ',jpmois
            endif
            stop 'dta_dic'
         endif

      ENDIF

      ! 2. Read monthly file
      ! -------------------

      IF( ( kt == nit000 .AND. nlecte == 0 ) .OR. imois /= ntem1 ) THEN
         nlecte = 1

         ! Calendar computation
         
         ntem1 = imois        ! first file record used 
         ntem2 = ntem1 + 1    ! last  file record used
         ntem1 = MOD( ntem1, iman )
         IF( ntem1 == 0 )   ntem1 = iman
         ntem2 = MOD( ntem2, iman )
         IF( ntem2 == 0 )   ntem2 = iman
         
         ! Read monthly nitrate data
!!DB
         call ncdf_read(cl_tdata,'vo_DIC',dicdta(:,:,:,1),-ntem1,status)
         call ncdf_read(cl_tdata,'vo_DIC',dicdta(:,:,:,2),-ntem2,status)

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' Monthly DIC input OK'
         IF(lwp) WRITE(numout,*)

!DL Do I need to do something about the mask??
         !                                  ! Mask
         DO jl = 1, 2
            dicdta(:,:,:,jl) = dicdta(:,:,:,jl) * tmask(:,:,:)
            dicdta(:,:,jpk,jl) = 0.
         END DO

      ENDIF

      ! 2. At every time step compute DIC data
      ! ----------------------------------------------
      zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
!!DB: Note that for monthly data it is not necessary to do this every kt.
!!    So in the future use a code fragment like: 
!      if(mod(kt,int(rday/rdt)) == 0) then !do interpolation ~ once per day
!         t_dta(:,:,:) = (1.-zxy) * temdta(:,:,:,1) + zxy * temdta(:,:,:,2)
!      endif
!DL Need to clarify this one
      c_dta(:,:,:) = (1.-zxy) * dicdta(:,:,:,1) + zxy * dicdta(:,:,:,2)

  END SUBROUTINE dta_dic
#endif
 
  SUBROUTINE river_no3( kt, jn )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE river_no3  ***
      !!                    
      !! ** Purpose :   New routine to read monthly nitrate data at Quebec City
      !!                Diane Lavoie, December 2012
      !!
      !! ** Method  :   Read on unit numtdt the nitrate data 
      !!      onto the model grid.
      !!      Data begin at january. 
      !!      The value is centered at the middle of month.
      !!      At each time step, a linear interpolation is applied between 
      !!      two monthly values.
      !!
      !!----------------------------------------------------------------------

      !! * Modules used
      
      USE rivers

      !! * Arguments
      INTEGER, INTENT( in ) ::  kt, jn     ! ocean time-step
          
!Attention, faire de facon identique 'a river
!*******************************************

     !! * Local declarations
     INTEGER, PARAMETER ::   &
        jpmois = 12                    ! number of month
     INTEGER ::   ji, jj, jk            ! dummy loop indicies
!DL     INTEGER ::   ir_sl, i, ii            ! dummy loop indicies
     INTEGER ::   ir, i, ii            ! dummy loop indicies
     INTEGER ::   &
        imois, iman, itime, ik ,    &  ! temporary integers
        i15, ipi, ipj, ipk             !    "          "
#  if defined key_tradmp
     INTEGER ::   &
        il0, il1, ii0, ii1, ij0, ij1   ! temporary integers
# endif
     INTEGER read_year,read_month,nriver_bgcm

!DL     INTEGER, DIMENSION(jpmois) ::   istep
!DL     REAL(wp), DIMENSION(2) ::   tmp_qnitr
!DL     REAL(wp) ::   zxy, zl, zdate0, qnitr
     REAL(wp), DIMENSION(nrivermax,2) ::   tmp_bgcm !            !NL#1
     REAL(wp), DIMENSION(nrivermax,12) ::  riv_bgcm !            !NL#1
     REAL(wp)   zxy, zl, zdate0
!     REAL(wp), DIMENSION(jpi,jpj) ::   zlon,zlat
!     REAL(wp), DIMENSION(jpk) ::   zlev
  
     REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
     rivndta            ! nitrate data at two consecutive times, river
#if defined (OXYGEN)
     REAL tk,s,saturated_02
# endif

!!----------------------------------------------------------------------

! 0. Initialization
! -----------------
      iman  = jpmois
      i15   = nday / 16
      imois = nmonth + i15 - 1
      IF( imois == 0 )   imois = iman

      itime = jpmois
      ipi = jpiglo
      ipj = jpjglo
      ipk = jpk

! Need to determine local and glocal indices
!c------------------------------------------
      if (kt==nit000) call Riv_init(kt)
!DL      ir_sl=1
!      print*
!      print*, "In river_no3 subroutine"
!      print*, "river at:",iriv_a(1),jriv_a(1)," to ",iriv_b(1),jriv_b(1)
!      print*, 'do_river(ir_sl)=',do_river(ir_sl)
!      print*,'TOTOTOTOTO'
! --------------------------------------------


     ! Read monthly nitrate data
     ! -------------------------

      ! Calendar computation
        
        ntem1 = imois        ! first file record used 
        ntem2 = ntem1 + 1    ! last  file record used
        ntem1 = MOD( ntem1, iman )
        IF( ntem1 == 0 )   ntem1 = iman
        ntem2 = MOD( ntem2, iman )
        IF( ntem2 == 0 )   ntem2 = iman
        

!DL      open(9, file = 'quebec_nitr.dat', status = 'old')
!DL      read(9,*) ii,qnitr
!DL      if (ntem1.eq.1) then
!DL          tmp_qnitr(1)=qnitr
!DL          read(9,*) ii,qnitr
!DL          tmp_qnitr(2)=qnitr
!DL      elseif (ntem1.eq.12) then
!DL          tmp_qnitr(2)=qnitr
!DL          do i=2,ntem1
!DL           read(9,*) ii,qnitr
!DL          enddo
!DL          tmp_qnitr(1)=qnitr
!DL      else          
!DL          do i=2,ntem1
!DL           read(9,*) ii,qnitr
!DL          enddo
!DL          tmp_qnitr(1)=qnitr
!DL          read(9,*) ii,qnitr
!DL          tmp_qnitr(2)=qnitr
!DL      endif
!DL      close(9)

!       if(lwp) print*, tmp_qnitr(1), tmp_qnitr(2)

        ! ===============   !NL#1
        ! load all the nitrate (#riv, 1 to ntem1+1 months)
        if (jn.eq.3) then
            open(9, file = 'river_nox.dat', status = 'old')
#if defined OXYGEN
        elseif (jn.eq.9) then
            IF(lwp.and.kt==nit000) write(*,*) 'Water from rivers saturated in oxygen'
            open(9, file = 'river_nox.dat', status = 'old') ! only read the nriver_bgcm
            read(9,*) nriver_bgcm  !number of rivers read in
            if(nriver_bgcm > nrivermax) stop 'nrivermax not big enough'
            close(9)
#endif
#if defined key_carbon
        elseif (jn.eq.10) then
            open(9, file = 'river_DIC.dat', status = 'old')
#endif
        else; stop
!        else; stop('File for jn='//char(jn+48)//' is not define')
        endif

        if (jn.ne.9) then        !NL#1 :to not read the oxygen (not from a file)
            read(9,*) nriver_bgcm  !number of rivers read in
            if(nriver_bgcm > nrivermax) stop 'nrivermax not big enough'

            do i=1,min(12,ntem1+1) ! no need to load futher than the next month
            read(9,*) read_year,read_month,(riv_bgcm(ir,i),ir=1,nriver_bgcm)
!            if (i.ne.read_month) stop('Error in the readin of ''river_nox/dic.dat''')
            if (i.ne.read_month) stop
            enddo
            close(9)

            if (ntem1.eq.12) then
            tmp_bgcm(:,2)=riv_bgcm(:,1)
            tmp_bgcm(:,1)=riv_bgcm(:,ntem1)
            else
            tmp_bgcm(:,1)=riv_bgcm(:,ntem1)
            tmp_bgcm(:,2)=riv_bgcm(:,ntem1+1)
            endif
        endif

        !NL#1
        ! ===============

!!c ------------------------------------------
!!c Compute nitrate data at local domain level
!!c ------------------------------------------

      zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.

      do ir=1,nriver_bgcm !nriver_bgcm            !NL#1
        if(do_river(ir)) then
          do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))
            do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
              do jk = 1, mbathy_riv(ir)
                 tra(ji+ishift(ir),jj+jshift(ir),jk,jn) = (1.-zxy) * tmp_bgcm(ir,1) + zxy * tmp_bgcm(ir,2)
!temporary additon of diat to rivers
                 tra(ji+ishift(ir),jj+jshift(ir),jk,1) = 0.05
!temporary additon of flag to rivers
                 tra(ji+ishift(ir),jj+jshift(ir),jk,2) = 0.025
!temporary additon of NH4 to rivers
                 tra(ji+ishift(ir),jj+jshift(ir),jk,4) = 1. 
!temporary additon of meso to rivers
                 tra(ji+ishift(ir),jj+jshift(ir),jk,5) = 0.075
!temporary additon of micro to rivers
                 tra(ji+ishift(ir),jj+jshift(ir),jk,6) = 0.05
!temporary additon of PON to rivers
                 tra(ji+ishift(ir),jj+jshift(ir),jk,7) = 1. 
!temporary additon of DON to rivers
                 tra(ji+ishift(ir),jj+jshift(ir),jk,8) = 3. 
!DL                  tra(ji,jj,jk,jn) = (1.-zxy) * tmp_qnitr(1) + zxy * tmp_qnitr(2)
!temporary, test for O2 at QC
#if defined OXYGEN
                 ! tra(ji,jj,jk,9)= 300. 
!temporary, test for oxy in rivers (put the saturation(t,s))
!c Solubility calculation from (Benson and Krauss, 1984, L&O) in [micromol/l] or [mmmol/m3]

                    tk = tn(ji,jj,jk) + 273.15
                    s = sn(ji,jj,jk)
                    saturated_02 = EXP( -135.90205 &
                    + (1.575701*1E5)/tk + (-6.642308*1E7)/(tk**2) &
                    + (1.243800*1E10)/(tk**3) + (-8.621949*1E11)/(tk**4) &
                    - s * ( 0.017674 + -10.754/tk + 2140.7/(tk**2) ) )

                    tra(ji+ishift(ir),jj+jshift(ir),jk,9) = saturated_02

#endif
!temporary, test for DIC at QC
#if defined key_carbon
                 ! tra(ji,jj,jk,10) =1500. 
                 tra(ji+ishift(ir),jj+jshift(ir),jk,10) = 1300. 
#endif

              enddo
            enddo
          enddo
        endif
      enddo

!          IF(lwp) WRITE(*,*)
!          IF(lwp) WRITE(*,*) ' Monthly nitrate at Quebec OK'
!          IF(lwp) WRITE(*,*), 'ntem1,ntem2 ' ,ntem1, ntem2
!          IF(lwp) WRITE(*,*)


  END SUBROUTINE river_no3


! =================================
! V V V V V       NL#2     V V V V V V

   SUBROUTINE rst_bgcm_write( niter )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE rst bgcm write  ***
      !!                     
      !! ** Purpose :   Write bgcm restart fields in NetCDF format
      !!
      !! ** Method  :   Write in numwrs file each nstock time step in NetCDF
      !!      file, save fields which are necessary for restart
      !!
      !! History : NL  2013.08.21
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
      REAL(wp), DIMENSION(1) ::   zinfo(1)
      REAL(wp) ::  zsec, zdate0, zdt
      
      LOGICAL,PARAMETER :: USE_IOIPSL=.FALSE.  ! Use IOIPSL subroutines for restart output
      LOGICAL,PARAMETER :: USE_LIB_NCDF=.TRUE. ! Use lib_ncdf for restart output
      INTEGER :: lncdf_stat                    ! lib_ncdf call return status
      CHARACTER(LEN=80) :: restfilename
      CHARACTER(LEN=4) :: chnum

      REAL(wp),DIMENSION(jpi,jpj,35) :: zmoment
      INTEGER :: ji, jj,jn
      INTEGER :: bgcm_write_count = 0

      !!----------------------------------------------------------------------

! Job informations
      it0      = niter
      zinfo(1) = it0   ! iteration number
      !zinfo(2) = adatrj  ! cumulated duration of all previous restart runs
      !zinfo(3) = bgcm_write_count  ! count write


      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) 'bgcm_rst_write : write the bgcm restart file in NetCDF format ',   &
           'at it= ',niter,' date= ',ndastp
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'

!!DB
      write(chnum,'(I4.4)')bgcm_write_count
      restfilename = trim(cexper)//'_bgcm_rstfile_'//chnum//'.nc'

      CALL ncdf_create_bgcm_restart(restfilename, lncdf_stat)
      CALL ncdf_write(restfilename, 'info', zinfo, lncdf_stat)   

      CALL ncdf_write(restfilename, 'nav_lat', gphit, nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'nav_lon', glamt, nwrite, lncdf_stat)
      CALL ncdf_write(restfilename, 'nav_lev', gdept, lncdf_stat)
      CALL ncdf_write(restfilename, 'time', REAL(0), 1, lncdf_stat)
      CALL ncdf_write(restfilename, 'time_steps', REAL(0), 1, lncdf_stat)

      CALL ncdf_write(restfilename, 'trn', trn, jptra, nwrite, lncdf_stat)

      bgcm_write_count = bgcm_write_count + 1


   END SUBROUTINE rst_bgcm_write

  SUBROUTINE ncdf_create_bgcm_restart(filename, status)
  ! NL 22.08.2013
  ! ncdf_create_bgcm_restart builds a standard OPA bgcm restart file with all the default
  ! dimensions, variables and attributes.
     IMPLICIT NONE
    ! Subroutine argument declarations
    CHARACTER(LEN=*),INTENT(IN) :: filename
    INTEGER,INTENT(OUT) :: status

    ! Local declarations
    INTEGER :: ncid,    &  ! netCDF file ID
         varid,   &  ! ID of netCDF variable to be written to
         nfstat,  &  ! netCDF library call return status
         mpistat     ! MPI library call return status
    INTEGER,DIMENSION(1:9) :: dimids
    INTEGER,DIMENSION(1:36) :: varids
    CHARACTER(LEN=20) :: cal_type      ! Calendar type
    CHARACTER(LEN=30) :: timestamp     ! File timestamp
    CHARACTER(LEN=100) :: sec_since
    CHARACTER(LEN=100) :: tstp_since
    CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
    CHARACTER(LEN=3),PARAMETER :: &
         &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
         &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
    INTEGER :: ji

    
    ! Initializations
    status = NCDF_NOERR
    
    CALL ioget_calendar(cal_type)
    CALL ioget_timestamp(timestamp)
    WRITE (UNIT=sec_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE (UNIT=tstp_since, &
         FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
         &  nyear,nmonth,nday,0, 0, 0
    WRITE(t_origin, &
         &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
         &   nyear,months(nmonth),nday,0,0,0

    IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
       WRITE(100,*) 'NCDF DEBUG: Creating default restart file:', filename
    END IF
    
    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Creating file:', filename
       END IF
       ! Create the file
       nfstat = nf90_create(filename, nf90_clobber, ncid)
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining dimensions in file:', filename
       END IF
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(1))
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(2))
       nfstat = nf90_def_dim(ncid, 'z', jpkdta, dimids(3))
       nfstat = nf90_def_dim(ncid, 'jptra', jptra, dimids(8))
       nfstat = nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimids(4))
       nfstat = nf90_def_dim(ncid, 'x_a', 1, dimids(5))
       nfstat = nf90_def_dim(ncid, 'y_a', 1, dimids(6))
       nfstat = nf90_def_dim(ncid, 'z_a', 1, dimids(7))

       
       ! Define variables
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Defining variables in file:', filename
       END IF
       nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, (/ dimids(1), dimids(2) /), varids(1))
       nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, (/ dimids(1), dimids(2) /), varids(2))
       nfstat = nf90_def_var(ncid, 'nav_lev', nf90_float, (/ dimids(3) /), varids(3))
       nfstat = nf90_def_var(ncid, 'time',    nf90_float, (/ dimids(4) /), varids(4))
       nfstat = nf90_def_var(ncid, 'time_steps', nf90_float,(/ dimids(4) /), varids(5))
       nfstat = nf90_def_var(ncid, 'info', nf90_float,   (/ dimids(7) /),  varids(6))
       nfstat = nf90_def_var(ncid, 'trn', nf90_float, & 
            & (/ dimids(1), dimids(2),dimids(3),dimids(8),dimids(4) /),  varids(7))

       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Writing attributes in file:', filename
       END IF
       ! nav_lon
       nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
       nfstat = nf90_put_att(ncid, varids(1), 'valid_min', -1.800000E2)
       nfstat = nf90_put_att(ncid, varids(1), 'valid_max', 1.800000E2)
       nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')

       ! nav_lat
       nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
       nfstat = nf90_put_att(ncid, varids(2), 'valid_min', -9.000000E1)
       nfstat = nf90_put_att(ncid, varids(2), 'valid_max', 9.000000E1)
       nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')

       ! nav_lev
       nfstat = nf90_put_att(ncid, varids(3), 'units', 'model_levels')
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', 0.0)
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', 0.0)
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Model levels')

       ! time
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'Time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', '0001-JUL-01 00:00:00')

       ! time_steps
       nfstat = nf90_put_att(ncid, varids(5), 'units', TRIM(tstp_since))
       nfstat = nf90_put_att(ncid, varids(5), 'title', 'Time steps')
       nfstat = nf90_put_att(ncid, varids(5), 'tstep_sec', 1.877760E7)
       nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'Time step axis')
       nfstat = nf90_put_att(ncid, varids(5), 'time_origin', TRIM(t_origin))

       ! info
       nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
!DB as per OLD code
       do ji = 7, 7
          nfstat = nf90_put_att(ncid, varids(ji), 'missing_value', 1.000000E20)
       enddo

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL , 'Conventions', 'GDT 1.2')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(filename))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number_total', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_number', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_dimensions_ids', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_global', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_size_local', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_first', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_position_last', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_start', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_halo_size_end', 0)
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'DOMAIN_type','box')

       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_restart - Closing file:', filename
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF
    END IF

    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF

  END SUBROUTINE ncdf_create_bgcm_restart


  SUBROUTINE rst_bgcm_read
     !! NL 22.08.2013
     !! rst_bgcm_read READ THE OPA bgcm restart file with all the default
     !! dimensions, variables and attributes. This will restart the simulation
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
     !! * Modules used
     USE ioipsl

     CHARACTER(len=45)  ::  restfilename = 'restart_bgcm.nc'
     INTEGER :: lncdf_stat                     ! lib_ncdf routine return status
     REAL(wp), DIMENSION(1) ::   zinfo(1)
     REAL(wp), DIMENSION(1) ::   zinfo2(10)

      ! print info on the present and past job
     IF(lwp) WRITE(numout,*)
     IF(lwp) WRITE(numout,*) 'rst_bgcm_read : read the NetCDF restart file for bgcm'
     IF(lwp) WRITE(numout,*) '~~~~~~~~'
     
     IF(lwp) WRITE(numout,*) ' Info on the present job : '
     IF(lwp) WRITE(numout,*) '   job number          : ', no
     IF(lwp) WRITE(numout,*) '   time-step           : ', nit000
     IF(lwp) WRITE(numout,*) '   date ndastp         : ', ndastp    
     IF(lwp) WRITE(numout,*) '   *** restart option'
     SELECT CASE ( nrstdt )
     CASE ( 0 ) 
        IF(lwp) WRITE(numout,*) '   nrstdt = 0 no control of nit000'
     CASE ( 1 ) 
        IF(lwp) WRITE(numout,*) '   nrstdt = 1 we control the date of nit000'
     CASE ( 2 )
        IF(lwp) WRITE(numout,*) '   nrstdt = 2 the date adatrj is read in restart file'
     CASE DEFAULT
        IF(lwp) WRITE(numout,*) '  ===>>>> nrstdt not equal 0, 1 or 2 : no control of the date'
        IF(lwp) WRITE(numout,*) ' =======                   ========='
     END SELECT

     CALL ncdf_read(restfilename, 'info', zinfo, lncdf_stat)
!!DB -- check if found file
     if(lwp .AND. lncdf_stat /= 0)then
        write(numout,*)'                '
        write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(numout,*)'STOP: Problem reading restart bgcm file ', restfilename
        write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(numout,*)'                '
        stop
     endif

     ! NON LIB_NCDF RELATED CODE
     IF(lwp) WRITE(numout,*) ' Info on the restart file read : '
     IF(lwp) WRITE(numout,*) '   time-step           : ', int(zinfo(1))
     !IF(lwp) WRITE(numout,*) '   date ndastp         : ', zinfo(2)
     !IF(lwp) WRITE(numout,*) '   rst count           : ', zinfo(3)
     IF(lwp) WRITE(numout,*)
     
     ! Control of date
     IF( nit000 - int(zinfo(1))  /= 1 .AND. nrstdt /= 0 ) THEN
        IF(lwp) WRITE(numout,cform_err)
        IF(lwp) WRITE(numout,*) ' ===>>>> : problem with nit000 for the restart'
        IF(lwp) WRITE(numout,*) ' verify the restart file or rerun with nrstdt = 0 (namelist)'
        WRITE(*,*) ' ===>>>> : problem with nit000 for the restart ',nit000,int(zinfo(1))+1
        nstop = nstop + 1
     ENDIF

     ! look if this "time step" is the same than in the restart of the ocean (restart.nc)
     CALL ncdf_read('restart.nc', 'info', zinfo2, lncdf_stat)
!!DB -- check if found file
     if(lwp .AND. lncdf_stat /= 0)then
        write(numout,*)'                '
        write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(numout,*)'STOP: Problem reading restart bgcm file ', restfilename
        write(numout,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(numout,*)'                '
        stop
     endif

     if (zinfo(1).ne.zinfo2(2) ) then
        IF(lwp) WRITE(numout,cform_err)
        IF(lwp) WRITE(numout,*) ' ===>>>> : problem with tstep for the restart'
        IF(lwp) WRITE(numout,*) ' restart.nc and restart_bgcm.nc don''t stop at the same tstep' 
     endif
     ! End of job info stuff

     ! Reading of the restart data
     print*, 'nwrite =', nwrite
     CALL ncdf_read(restfilename, 'trn', trn, jptra, nwrite, lncdf_stat)
     if(lncdf_stat /= 0)then
        write(*,*)'                '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'STOP: Problem reading trn variable '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'                '
        stop
     endif


  END SUBROUTINE rst_bgcm_read ! with the last rstfile saved


! ^ ^ ^ ^ ^       NL#2   ^ ^ ^ ^ ^   

END MODULE lib_bgcm_02

#endif 
