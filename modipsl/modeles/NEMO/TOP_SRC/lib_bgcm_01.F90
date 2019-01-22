!!DB 2009.10.09 -- code restructured to put most routines in
!! lib_bgcm_01.F90



!!DB 2008.08.18 ... 

!!This module does all the basic initialization for BGCM_01
!!It replaces the complicated series of calls to various *trc* subroutines
!!The calling sequence for this module is:
!! opa.F90: call ini_trc ----> these routines.
!!To keep some backward compatibility with generic tracer and biology code
!!this module is USE'd by initrc.F90 only:
!!                  #ifdef key_BGCM_01 (see initrc.F90)

!!NB: this is only used if key_BGCM_01 is defined
!!   (and key_passivetrc is defined)

!!Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!-1- I force the tracer model to use the same diffusion scheme as the dynamic model. 
!!    I found that the TVD advection scheme for tracers does not work, while the
!!    other 3 (cen2, muscl, muscl2) do. I hardwire the cen2 scheme which seems 
!!    slightly faster (see trctrp.F90)
!!-2- Related to -1- is the variable trcrat which is read from the
!!    namelist.passivetrc file. One must be careful with this variable. 
!!    This variable controls the tracer diffusivity:
!!    From passivetrc_substitute.h90:
!!    #elif defined key_traldf_smag
!!      SMAG scheme                    aht: 3D coefficient
!!    #  define fsahtrt(i,j,k)  trcrat * ahtt(i,j,k)
!!    DB's default for TS diffusion is the Smagorinsky scheme, in which a Pr-# is used that
!!    scales Ah versus Am: i.e. Ah(:,:,:) = Pr*Am(:,:,:), where Pr = 0.1 as
!!    recommended by Z.Wang. Therefore the assignment: trcrat*ahtt really is
!!    trcrat*Pr*Am. Currently trcrat = 15 -----> A_tracer = 1.5*A_momentum, which
!!    is considered reasonable to start. 
!!-3- I do not perform "consistency" controls (e.g. call *_ctl()). Some consistency is
!!    achieved by virtue of -1- which effectively uses the consistency from TS adv/diff.
!!-4- Some subroutines do not yet exist -- they are noted, or are obvious. 
!!-5- Some keys are verbotten -- they are noted, or are obvious. 
!!

!!END Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!To see some of the older calls and structure see (e.g.)
!!(drakes:dbrick) OLD_OPA_CODE_VERSIONS/

#if defined key_passivetrc

MODULE lib_bgcm_01
   !!================================================
   !!
   !! Initialization of the BGCM 
   !!================================================
   !!--------------------------------------------------------------
   !! * Modules used
   !! ==============
   USE oce_trc
   USE trc
   
   IMPLICIT NONE
   PRIVATE

!!DB: NPZ params 
   REAL(wp), PUBLIC  :: I_sfce, Ibar
!  These date to pre 04.01.20
   REAL(wp), PUBLIC  ::    & !Vm, ks, Rm, lambda, M_PHY, gamma, M_ZOO, kext
        Vm = 2.0/86400.0,  & !! convert per day to per second
        Rm = 0.5/86400.0,  &
        ks = 1.0,          &!! \mu Mol N per liter
        lambda = 0.3,      & !! grazing (OLD=0.2)
        gamma = 0.7,       &
        M_ZOO = 0.04/86400.0,  &!! Zoo mort
        M_PHY = 0.05/86400.0,  & !! M_PHY
        kext = 0.1,        & !!in m
        aaa = 0.025,       &!! initial slope of P-I curve (< 0.08)
        !       aaa = 0.080,       &!! initial slope of P-I curve (< 0.08)
        !       aaa = 0.050,       &!! initial slope of P-I curve (< 0.08)
        PAR = 0.43,        & !!from Oschlies
        w_sink = -1.5/86400.0 !!sinking rate in m/day (OLD=-2.5)
   
   !!DB: 2008.11.19: wmask for step_npz; computed in bgcm_ini
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
   REAL(wp), PUBLIC :: e_trn(jpj,jpk,jptra), w_trn(jpj,jpk,jptra), s_trn(jpi,jpk,jptra)
   REAL(wp), PUBLIC :: e_trb(jpj,jpk,jptra), w_trb(jpj,jpk,jptra), s_trb(jpi,jpk,jptra)
!!NEW 2009.10.07
   CHARACTER (len=80), PUBLIC :: bgcm_fname
   
   
   !! * Accessibility
   PUBLIC ini_trc, bgcm_setup, bgcm_ini, extract_boundary_vals, &
        assign_boundary_vals, update_boundary_vals, bgcm_N_obc

   
CONTAINS

!!DB: Replace trc_lec, trc_trp_lec, trc_trp_ctl ... with 1 routine
!!One of the things I do in this routine is force the same adv/diff schemes
!!as are used by the TS module.
!!Initialized netcdf output file as well

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
!!DBG
!      if(lwp) write(numout,*)'DBG: BGCM init: called tra_adv_ctl'
      call tra_adv_ctl 

      IF(lwp) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' ROUTINE BGCM_SETUP'
         WRITE(numout,*) ' **************'
         WRITE(numout,*) ' '
         WRITE(numout,*) ' namelist for passive tracer'
         WRITE(numout,*) ' ***************************'
         WRITE(numout,*) ' '
      ENDIF

      numnat=80
!!DB: new
      clname='namelist.BGCM_01'
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
!!     the following assingment
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


      !! namelist of SMS
!!DB
!!input of params for Source-Minus-Sink (SMS) part of code
!! ... TO DO ...
!!DBG
      if(lwp) write(numout,*)'DBG: bgcm_setup: no SMS setup yet'
!      call  some_routine_to_be_written(put at bottom of this module)

!!DB: create ncdf file 
      bgcm_fname = trim(cexper)//'_bgcm_01.nc'
      call ncdf_create_bgcm_file(status)


   END SUBROUTINE bgcm_setup


!!DB: Initialize model here
   SUBROUTINE bgcm_ini
      !!---------------------------------------------------------------------
      !!              
      !! ** Purpose : Specific initialization for BGCM_01
      !!
      !! * local declarations

     USE lbclnk

     INTEGER ::                   & 
          ji ,jj ,jk ,jn, jl        ! dummy loop indices  
     !!---------------------------------------------------------------------
!!DB -- example of using pointers to make NPZ code more explicit
     real(wp), pointer :: N(:,:,:), P(:,:,:), Z(:,:,:)
     real(wp), pointer :: N0(:,:,:), P0(:,:,:), Z0(:,:,:)
     integer :: i1,i2,j1,j2
     real(wp) :: c0,c1,s0,s1

!!DB -- assign pointers to default trn, trb arrays
     N => trn(:,:,:,1)
     P => trn(:,:,:,2)
     Z => trn(:,:,:,3)
     N0 => trb(:,:,:,1)
     P0 => trb(:,:,:,2)
     Z0 => trb(:,:,:,3)


!!DB: 2008.11.19 -- compute wmask
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


!!DB 2009.06.03 -- try a formula that decays offshelf
!!DB - from 1D model 
!!NB:  different profile in deep water
!     do ji = 1, jpi
!        do jj = 1, jpj
!           do jk = 1, jpk
!              N(ji,jj,jk) = min(5.0 + (10.-5.)/200.0 * gdept(jk),10.0)
!           enddo
!           if(mbathy(ji,jj) > 2) then
!              if(gdept(mbathy(ji,jj)-1) >= 500.0) N(ji,jj,:)=0.1
!!              if(gdept(mbathy(ji,jj)-1) >= 1500.0) N(ji,jj,:)=0.1
!           endif
!        enddo
!     enddo
!!DB 2009.06.03 -- try a formula that decays offshelf
!     c0=5.0
!     s0=(10.-5.)/250.
!     do ji = 1, jpi
!        do jj = 1, jpj
!           do jk = 1, jpk
!              N(ji,jj,jk) = min(c0 + s0 * gdept(jk),15.0)
!           enddo
!           if(mbathy(ji,jj) > 2) then
!              if(gdept(mbathy(ji,jj)-1) >= 500.0) then
!                 c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/2500.)) + 0.1
!                 s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/2500.))
!                 do jk = 1, jpk
!                    N(ji,jj,jk) = min(c1 + s1 * gdept(jk),15.0)
!                 enddo
!              endif
!           endif
!        enddo
!     enddo

!!DB 2009.06.08 -- try a linear N(z) formula that decays offshelf
!!with time-dependent params, which will also be applied as OBCs
!!NB: lots of HARDWIRED stuff at this time
!06.10: started with 2500 as the decay scale which was too broad
!replaced 2500 with 1500
!c0 == N(z=0) = c1+c2*cos(pi*((yd-c3)/180));
!c1=3.0;c2=2.25;c3=25.0;  
     c0 = 3.0 + 2.25*cos((nday_year-25.0)*rpi/180.)
     s0=(15.-c0)/150.
     do ji = 1, jpi
        do jj = 1, jpj
           do jk = 1, jpk
              N(ji,jj,jk) = min(c0 + s0 * gdept(jk),15.0)
           enddo
           if(mbathy(ji,jj) > 2) then
              if(gdept(mbathy(ji,jj)-1) >= 500.0) then
                 c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.)) + 0.1
                 s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.))
                 do jk = 1, jpk
                    N(ji,jj,jk) = min(c1 + s1 * gdept(jk),15.0)
                 enddo
              endif
           endif
        enddo
     enddo
!!DBG
     if(lwp) write(221+narea,*)'DBG: BGCM initial N(z=0), yd: ',c0,nday_year


    P(:,:,:) = 1.0e-2
    Z(:,:,:) = 1.0e-2

!!DBG
!    P(:,:,:) = 5.0
!    Z(:,:,:) = 2.0
!    N(:,:,:) = 10.0


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


!!DB: If some BGCM data are to be read from a file, do it here. 
!#if defined key_trc_dta
!!   Initialization of tracer from a file
!      CALL dta_trc( nit000 )
!      DO  jk = 1, jptra
!        IF( lutini(jk) ) THEN 
!           trn(:,:,:,jk) = trdta(:,:,:,jk)*tmask(:,:,:)
!        ENDIF
!      END DO
!#endif

!! before field :
!! -------------
     trb(:,:,:,:) = trn(:,:,:,:)

     if( lwp ) then
        write(numout,*) ' '
        write(numout,*) 'BGCM_01 initialisation done '
!!DBG
!        write(numout,*) 'DBG: Fields constant everywhere '
        write(numout,*) ' '
     endif
      
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
END SUBROUTINE extract_boundary_vals

!!DB: Put the ?_trn boundary arrays values into tra at the boundary indices
!!Routine retained for possible future use
!!For example, if using N data along Open Boundaries then a version of this
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

END SUBROUTINE assign_boundary_vals

!!DB: 2009.04.20 -- copy values 1 cell in to OB positions
!!use "after" arrays
!!This routine is effectively the OBC code. It is called in trcnxt.F90
SUBROUTINE update_boundary_vals
  USE obc_oce
#  include "obc_vectopt_loop_substitute.h90"
  IMPLICIT NONE
  INTEGER :: ji,jj,jk,jn

!!DBG: assign to outer(inner) values as well
  do ji = fs_nie0+1, fs_nie1+1 !as in obctra
     tra(ji,:,:,:) = tra(ji-1,:,:,:) 
  enddo
  do ji = niw0, niw1 ! isolate processor
     tra(ji,:,:,:) = tra(ji+1,:,:,:) 
  enddo
  do jj = njs0, njs1
     tra(:,jj,:,:) = tra(:,jj+1,:,:)
  enddo


END SUBROUTINE update_boundary_vals

!!2009.08
!!Update N OBC vals and assign to tra()
!!To be called in trcnxt()
  SUBROUTINE bgcm_N_obc ( kt )
    USE obc_oce         !open boundary condition variables
    USE lbclnk
#  include "obc_vectopt_loop_substitute.h90"
    
    !!------------------------------------------------------------------------------
    !! * Arguments
    INTEGER, INTENT( in ) ::   kt
    
    !! * Local declaration
    INTEGER ::   ji, jj, jk, jn
    REAL :: c0,c1,s0,s1


!!NB: eventually do this once/day
!!DB 2009.06.08 -- try a linear N(z) formula that decays offshelf
!!with time-dependent params, which will also be applied as OBCs
!!NB: lots of HARDWIRED stuff at this time
!06.10: started with 2500 as the decay scale which was too broad
!replaced 2500 with 1500
!c0 == N(z=0) = c1+c2*cos(pi*((yd-c3)/180));
!c1=3.0;c2=2.25;c3=25.0;  
     c0 = 3.0 + 2.25*cos((nday_year-25.0)*rpi/180.)
     s0=(15.-c0)/150.
!!isolate processors for each OB
     do ji = fs_nie0, fs_nie1
        do jj = 1, jpj
           do jk = 1, jpk
              e_trn(jj,jk,1) = min(c0 + s0 * gdept(jk),15.0)
           enddo
           if(mbathy(ji,jj) > 2) then
              if(gdept(mbathy(ji,jj)-1) >= 500.0) then
                 c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.)) + 0.1
                 s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.))
                 do jk = 1, jpk
                    e_trn(jj,jk,1) = min(c1 + s1 * gdept(jk),15.0)
                 enddo
              endif
           endif
        enddo
        tra(ji+1,:,:,1) = e_trn(:,:,1)
     enddo
     do ji = niw0, niw1
        do jj = 1, jpj
           do jk = 1, jpk
              w_trn(jj,jk,1) = min(c0 + s0 * gdept(jk),15.0)
           enddo
           if(mbathy(ji,jj) > 2) then
              if(gdept(mbathy(ji,jj)-1) >= 500.0) then
                 c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.)) + 0.1
                 s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.))
                 do jk = 1, jpk
                    w_trn(jj,jk,1) = min(c1 + s1 * gdept(jk),15.0)
                 enddo
              endif
           endif
        enddo
        tra(ji,:,:,1) = w_trn(:,:,1)
     enddo
     do jj = njs0, njs1
        do ji = 1, jpi
           do jk = 1, jpk
              s_trn(ji,jk,1) = min(c0 + s0 * gdept(jk),15.0)
           enddo
           if(mbathy(ji,jj) > 2) then
              if(gdept(mbathy(ji,jj)-1) >= 500.0) then
                 c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.)) + 0.1
                 s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.))
                 do jk = 1, jpk
                    s_trn(ji,jk,1) = min(c1 + s1 * gdept(jk),15.0)
                 enddo
              endif
           endif
        enddo
        tra(:,jj,:,1) = s_trn(:,:,1)
     enddo
!!DBG
     if(lwp) write(221+narea,*)'DBG: BGCM kt N(z=0), yd: ',kt, c0,nday_year

!!DBG: every day (hardwired) output a few vertical profiles to monitor output
     if(mod(kt-1,180) == 0) then
        
!!West
        do ji = niw0, niw1
           do jj = 2, jpj, int(jpj/4)
              do jk = 1, jpk
                 write(5000+narea,'(3(i8,1x),10(e10.4,2x))')kt,jj,jk, (w_trn(jj,jk,jn),jn=1,jptra), &
                      wmask(ji,jj,jk), gdept(mbathy(ji,jj)-1),float(mbathy(ji,jj))
              enddo
           enddo
        enddo
!!South
        do jj = fs_njs0, fs_njs1
           do ji = 2, jpi, int(jpi/4)
              do jk = 1, jpk
                 write(6000+narea,'(3(i8,1x),10(e10.4,2x))')kt,ji,jk, (s_trn(ji,jk,jn),jn=1,jptra), &
                      wmask(ji,jj,jk), gdept(mbathy(ji,jj)-1),float(mbathy(ji,jj))
              enddo
           enddo
        enddo
!!East
        do ji = fs_nie0, fs_nie1
           do jj = 2, jpj, int(jpj/4)
              do jk = 1, jpk
                 write(7000+narea,'(3(i8,1x),10(e10.4,2x))')kt,jj,jk, (e_trn(jj,jk,jn),jn=1,jptra), &
                      wmask(ji,jj,jk), gdept(mbathy(ji,jj)-1),float(mbathy(ji,jj))
              enddo
           enddo
        enddo
           
     endif


  END SUBROUTINE bgcm_N_obc

!!DB: 2008.08.28 ...
!!Create output file for BGCM_01
!!The main variable is the trn(:,:,:,jptra) array which is what is integrated
!!by the tracer module. 
!!N,P,Z can be used instead (see below and bgcm_01_model.F90) but using
!!both takes up double the disk space (even if no output is done),
!!so define and use only one of these possibilities.
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
       WRITE(100,*) 'NCDF DEBUG: Creating output file:', bgcm_fname
       CALL FLUSH
    END IF

    ! Only processor 0 does anything
    IF(nproc == 0) THEN
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Creating file:', bgcm_fname
          CALL FLUSH
       END IF
       ! Create the file
       nfstat = nf90_create(bgcm_fname, nf90_clobber, ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF
       
       ! Define dimensions
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Defining dimensions in file:', bgcm_fname
          CALL FLUSH
       END IF
       !nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       !nfstat = nf90_def_dim(ncid, 'jpkdta', jpkdta, dimids(2))
       !nfstat = nf90_def_dim(ncid, 'jpjdta', jpjdta, dimids(3))
       !nfstat = nf90_def_dim(ncid, 'jpidta', jpidta, dimids(4))
       nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
       nfstat = nf90_def_dim(ncid, 'jptra', jptra, dimids(2))  !!CN: Changed dim creation order
       nfstat = nf90_def_dim(ncid, 'z', jpk, dimids(3))          !!bgcm model
       nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(4))
       nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(5))
       !nfstat = nf90_def_dim(ncid, 'jptra', jptra, dimids(5)) !!CN: Changed dim creation order
       
       ! Define variables
       ! CN: changed order of dims in dim lists for var creation
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Defining variables in file:', bgcm_fname
          CALL FLUSH
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
!!DBG: Choose between NPZ or trn variable declarations 
!!NB: must be consistent with output in bgcm_01_model
       jf = 1
       if(jf == 0) then
          nfstat = nf90_def_var(ncid, 'trn', nf90_float, &
               (/ dimids(5), dimids(4), dimids(3), dimids(2), dimids(1) /), &
               varids(5))
          varnum = 5
       else
          nfstat = nf90_def_var(ncid, 'N', nf90_float, &
               (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
               varids(6))
          nfstat = nf90_def_var(ncid, 'P', nf90_float, &
               (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
               varids(7))
          nfstat = nf90_def_var(ncid, 'Z', nf90_float, &
               (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
               varids(8))
          varnum = 8
       endif

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
            'time step date (when output is writen) in year/month/day aammjj (decimal day)')
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

!!DB: HARDWIRED: add model params to .nc file
!!Use z=jpk as dimension so no new dimension need be defined.
!!If jpk becomes too small, use x or y
!!Use long_name for param names (other attribute names did not work)
       varnum = varnum + 1
       nfstat = nf90_def_var(ncid, 'model_params', nf90_float, &
            (/ dimids(3) /), varids(varnum))
       nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Hardwired -- 11 params: &
            Vm,Rm,ks,lambda,gamma,M_ZOO,M_PHY,kext,aaa,PAR,w_sink')
       nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'vals in model_params(1:11)')
       
       ! Add attributes
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Writing attributes in file:', bgcm_fname
          CALL FLUSH
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
       nfstat = nf90_put_att(ncid, varids(3), 'valid_min', minval(gdept))      !for bgcm model
       nfstat = nf90_put_att(ncid, varids(3), 'valid_max', maxval(gdept))      !for bgcm model
       nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
       nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical Tracer levels')

!!DB: temporary attributes (REM: jf set above)
       ! trn 
       if(jf == 0) then
          nfstat = nf90_put_att(ncid, varids(5), 'units', 'XXX')
          nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
          nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
          nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
          nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'BGCM_01 tracers N-P-Z')
          nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'trn')
          nfstat = nf90_put_att(ncid, varids(5), 'online_operation', op_type)
          nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TNZYX')
          nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', float(int_opp))
          nfstat = nf90_put_att(ncid, varids(5), 'interval_write', float(int_wri))
          nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
       else
       ! N
          nfstat = nf90_put_att(ncid, varids(6), 'units', 'XXX')
          nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
          nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
          nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
          nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'BGCM_01-- N')
          nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'N')
          nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
          nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TZYX')
          nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', float(int_opp))
          nfstat = nf90_put_att(ncid, varids(6), 'interval_write', float(int_wri))
          nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter deptht nav_lat nav_lon')
          !!DB: 
          ! P
          nfstat = nf90_put_att(ncid, varids(7), 'units', 'XXX')
          nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
          nfstat = nf90_put_att(ncid, varids(7), 'valid_min', 1.000000E20 )
          nfstat = nf90_put_att(ncid, varids(7), 'valid_max', -1.000000E20)
          nfstat = nf90_put_att(ncid, varids(7), 'long_name', 'BGCM_01-- P')
          nfstat = nf90_put_att(ncid, varids(7), 'short_name', 'P')
          nfstat = nf90_put_att(ncid, varids(7), 'online_operation', TRIM(op_type))
          nfstat = nf90_put_att(ncid, varids(7), 'axis', 'TZYX')
          nfstat = nf90_put_att(ncid, varids(7), 'interval_operation', float(int_opp))
          nfstat = nf90_put_att(ncid, varids(7), 'interval_write', float(int_wri))
          nfstat = nf90_put_att(ncid, varids(7), 'associate', 'time_counter deptht nav_lat nav_lon')
          !!DB: 
          ! Z
          nfstat = nf90_put_att(ncid, varids(8), 'units', 'XXX')
          nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
          nfstat = nf90_put_att(ncid, varids(8), 'valid_min', 1.000000E20 )
          nfstat = nf90_put_att(ncid, varids(8), 'valid_max', -1.000000E20)
          nfstat = nf90_put_att(ncid, varids(8), 'long_name', 'BGCM_01-- Z')
          nfstat = nf90_put_att(ncid, varids(8), 'short_name', 'Z')
          nfstat = nf90_put_att(ncid, varids(8), 'online_operation', TRIM(op_type))
          nfstat = nf90_put_att(ncid, varids(8), 'axis', 'TZYX')
          nfstat = nf90_put_att(ncid, varids(8), 'interval_operation', float(int_opp))
          nfstat = nf90_put_att(ncid, varids(8), 'interval_write', float(int_wri))
          nfstat = nf90_put_att(ncid, varids(8), 'associate', 'time_counter deptht nav_lat nav_lon')
       endif

       ! time_counter
       nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
       nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
       nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
       nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
       nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))

       ! global
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(bgcm_fname))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
       nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'MISC', 'DB-created file')


       ! Close file
       IF(DEBUG_OUT .EQV. .TRUE.) THEN
          WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Closing file:', bgcm_fname
          CALL FLUSH
       END IF
       nfstat = nf90_close(ncid)
       IF(nfstat /= nf90_noerr) THEN
          status = NCDF_NFERR
          RETURN
       END IF
    END IF


!!Write grid info to file
    CALL ncdf_write(bgcm_fname, 'nav_lat', gphit, -1, status)
    CALL ncdf_write(bgcm_fname, 'nav_lon', glamt, -1, status)
    CALL ncdf_write(bgcm_fname, 'deptht', gdept, status)
!!Write param_vals to model_params
    param_val(1) = Vm
    param_val(2) = Rm
    param_val(3) = ks
    param_val(4) = lambda
    param_val(5) = gamma
    param_val(6) = M_ZOO
    param_val(7) = M_PHY
    param_val(8) = kext
    param_val(9) = aaa
    param_val(10) = PAR
    param_val(11) = w_sink
    CALL ncdf_write(bgcm_fname, 'model_params', param_val, status)


    ! Sync up processors before returning from subroutine
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
    IF(mpistat /= 0) THEN
       status = NCDF_MPERR
       RETURN
    END IF
    
  END SUBROUTINE ncdf_create_bgcm_file

   
   SUBROUTINE ini_trc
      !!---------------------------------------------------------------------
      !!
      !!                       ROUTINE ini_trc
      !!                     ******************
      !!
      !!  PURPOSE :  initialize the BGCM_01 model
      !!  DB: Replaces regular ini_trc when key_BGCM_01 is defined
      !!---------------------------------------------------------------------


      !! 0.b PRINT the number of tracer
      !! ------------------------------

      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) ' *** number of passive tracer jptra = ',jptra
      IF(lwp) WRITE(numout,*) ' '

      call bgcm_setup 

      call bgcm_ini

!!If restarting (not thought about yet) assume that the following call
!!will overwrite anything wrong that bgcm_ini might have done
      if( lrsttr ) THEN
!         CALL bgcm_rst       !!...TO DO ...
      endif

   END SUBROUTINE ini_trc

END MODULE lib_bgcm_01

#endif 
