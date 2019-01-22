


PGF90 (Version     10.5)          09/23/2010  09:49:17      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir /usr/local/include
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /usr/local/include
          -idir /usr/include/mpich2-x86_64
          -idir /usr/include/mpich2-x86_64
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: lib_bgcm_01.F90

(    1) !!DB 2009.10.09 -- code restructured to put most routines in
(    2) !! lib_bgcm_01.F90
(    3) 
(    4) 
(    5) 
(    6) !!DB 2008.08.18 ... 
(    7) 
(    8) !!This module does all the basic initialization for BGCM_01
(    9) !!It replaces the complicated series of calls to various *trc* subroutines
(   10) !!The calling sequence for this module is:
(   11) !! opa.F90: call ini_trc ----> these routines.
(   12) !!To keep some backward compatibility with generic tracer and biology code
(   13) !!this module is USE'd by initrc.F90 only:
(   14) !!                  #ifdef key_BGCM_01 (see initrc.F90)
(   15) 
(   16) !!NB: this is only used if key_BGCM_01 is defined
(   17) !!   (and key_passivetrc is defined)
(   18) 
(   19) !!Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   20) !!-1- I force the tracer model to use the same diffusion scheme as the dynamic model. 
(   21) !!    I found that the TVD advection scheme for tracers does not work, while the
(   22) !!    other 3 (cen2, muscl, muscl2) do. I hardwire the cen2 scheme which seems 
(   23) !!    slightly faster (see trctrp.F90)
(   24) !!-2- Related to -1- is the variable trcrat which is read from the
(   25) !!    namelist.passivetrc file. One must be careful with this variable. 
(   26) !!    This variable controls the tracer diffusivity:
(   27) !!    From passivetrc_substitute.h90:
(   28) !!    #elif defined key_traldf_smag
(   29) !!      SMAG scheme                    aht: 3D coefficient
(   30) !!    #  define fsahtrt(i,j,k)  trcrat * ahtt(i,j,k)
(   31) !!    DB's default for TS diffusion is the Smagorinsky scheme, in which a Pr-# is used that
(   32) !!    scales Ah versus Am: i.e. Ah(:,:,:) = Pr*Am(:,:,:), where Pr = 0.1 as
(   33) !!    recommended by Z.Wang. Therefore the assignment: trcrat*ahtt really is
(   34) !!    trcrat*Pr*Am. Currently trcrat = 15 -----> A_tracer = 1.5*A_momentum, which
(   35) !!    is considered reasonable to start. 
(   36) !!-3- I do not perform "consistency" controls (e.g. call *_ctl()). Some consistency is
(   37) !!    achieved by virtue of -1- which effectively uses the consistency from TS adv/diff.
(   38) !!-4- Some subroutines do not yet exist -- they are noted, or are obvious. 
(   39) !!-5- Some keys are verbotten -- they are noted, or are obvious. 
(   40) !!
(   41) 
(   42) !!END Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   43) 
(   44) !!To see some of the older calls and structure see (e.g.)
(   45) !!(drakes:dbrick) OLD_OPA_CODE_VERSIONS/
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 2

(   47) #if defined key_passivetrc
(   48) 
(   49) MODULE lib_bgcm_01
(   50)    !!================================================
(   51)    !!
(   52)    !! Initialization of the BGCM 
(   53)    !!================================================
(   54)    !!--------------------------------------------------------------
(   55)    !! * Modules used
(   56)    !! ==============
(   57)    USE oce_trc
(   58)    USE trc
(   59)    
(   60)    IMPLICIT NONE
(   61)    PRIVATE
(   62) 
(   63) !!DB: NPZ params 
(   64)    REAL(wp), PUBLIC  :: I_sfce, Ibar
(   65) !  These date to pre 04.01.20
(   66)    REAL(wp), PUBLIC  ::    & !Vm, ks, Rm, lambda, M_PHY, gamma, M_ZOO, kext
(   67)         Vm = 2.0/86400.0,  & !! convert per day to per second
(   68)         Rm = 0.5/86400.0,  &
(   69)         ks = 1.0,          &!! \mu Mol N per liter
(   70)         lambda = 0.3,      & !! grazing (OLD=0.2)
(   71)         gamma = 0.7,       &
(   72)         M_ZOO = 0.04/86400.0,  &!! Zoo mort
(   73)         M_PHY = 0.05/86400.0,  & !! M_PHY
(   74)         kext = 0.1,        & !!in m
(   75)         aaa = 0.025,       &!! initial slope of P-I curve (< 0.08)
(   76)         !       aaa = 0.080,       &!! initial slope of P-I curve (< 0.08)
(   77)         !       aaa = 0.050,       &!! initial slope of P-I curve (< 0.08)
(   78)         PAR = 0.43,        & !!from Oschlies
(   79)         w_sink = -1.5/86400.0 !!sinking rate in m/day (OLD=-2.5)
(   80)    
(   81)    !!DB: 2008.11.19: wmask for step_npz; computed in bgcm_ini
(   82)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk+1) ::  wmask
(   83) 
(   84) !!DB
(   85) !!Shared module variables, default values. NB: they are reset below
(   86)    !! Advection
(   87)    LOGICAL, PUBLIC ::   &
(   88)       ln_trcadv_cen2   = .FALSE. ,  & !!: 2nd order centered scheme flag
(   89)       ln_trcadv_tvd    = .TRUE. ,  &  !: TVD scheme flag
(   90)       ln_trcadv_muscl  = .FALSE. ,  &  !: MUSCL scheme flag
(   91)       ln_trcadv_muscl2 = .FALSE. ,  &  !: MUSCL2 scheme flag
(   92)       ln_trcadv_smolar = .FALSE.        !: Smolarkiewicz scheme flag (DB: must be .FALSE.)
(   93) 
(   94)    !! Lateral diffusion
(   95)    LOGICAL , PUBLIC ::              & !!: ** lateral mixing namelist (nam_trcldf) **
(   96)       ln_trcldf_diff  = .FALSE. ,   &  !: flag of perform or not the lateral diff.
(   97)       ln_trcldf_lap   = .TRUE.  ,   &  !: laplacian operator
(   98)       ln_trcldf_bilap = .FALSE. ,   &  !: bilaplacian operator
(   99)       ln_trcldf_level = .FALSE. ,   &  !: iso-level direction
(  100)       ln_trcldf_hor   = .TRUE. ,   &  !: horizontal (geopotential) direction
(  101)       ln_trcldf_iso   = .FALSE.         !: iso-neutral direction
(  102) 
(  103)    LOGICAL , PUBLIC ::              & !!: flag of the lateral diff. scheme used
(  104)       l_trcldf_lap         ,        &  !: iso-level laplacian operator






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 3

(  105)       l_trcldf_bilap       ,        &  !: iso-level bilaplacian operator
(  106)       l_trcldf_bilapg      ,        &  !: geopotential bilap. (s-coord)
(  107)       l_trcldf_iso         ,        &  !: iso-neutral laplacian or horizontal lapacian (s-coord)
(  108)       l_trczdf_iso         ,        &  !: idem for the vertical component
(  109)       l_trczdf_iso_vo      ,        &  !: idem with vectopt_memory
(  110)       l_trcldf_iso_zps                 !: iso-neutral laplacian (partial steps)
(  111) 
(  112)    !! Vertical diffusion
(  113)    LOGICAL , PUBLIC ::              & !!: nam_trczdf: vertical diffusion
(  114)       ln_trczdf_exp = .FALSE.          !: explicit vertical diffusion scheme flag
(  115) 
(  116)    INTEGER, PUBLIC ::               & !!: namzdf:  vertical diffusion
(  117)       n_trczdf_exp = 3                 !: number of sub-time step (explicit time stepping)
(  118) 
(  119)    LOGICAL, PUBLIC ::               &  !:
(  120)       l_trczdf_exp     = .FALSE. ,  &  !: explicit vertical diffusion
(  121)       l_trczdf_imp     = .FALSE.       !: implicit vertical diffusion
(  122) 
(  123) !!DB: not accounted for yet
(  124) #if defined key_trcdmp
(  125)    !! Newtonian damping
(  126)    INTEGER  , PUBLIC ::             & !!: * newtonian damping namelist (nam_trcdmp) *
(  127)       ndmptr   =   -1 ,             &  !: = 0/-1/'latitude' for damping over tracers
(  128)       ndmpftr  =    2 ,             &  !: = 1 create a damping.coeff NetCDF file 
(  129)       nmldmptr =    0                  !: = 0/1/2 flag for damping in the mixed layer
(  130) 
(  131)    REAL(wp) , PUBLIC ::             & !!:  * newtonian damping namelist *
(  132)       sdmptr   =   50.,             &  !: surface time scale for internal damping (days)
(  133)       bdmptr   =  360.,             &  !: bottom time scale for internal damping (days)
(  134)       hdmptr   =  800.                 !: depth of transition between sdmp and bdmp (meters)
(  135) #endif
(  136) 
(  137) !!DB -- arrays of boundary trn values
(  138)    REAL(wp), PUBLIC :: e_trn(jpj,jpk,jptra), w_trn(jpj,jpk,jptra), s_trn(jpi,jpk,jptra)
(  139)    REAL(wp), PUBLIC :: e_trb(jpj,jpk,jptra), w_trb(jpj,jpk,jptra), s_trb(jpi,jpk,jptra)
(  140) !!NEW 2009.10.07
(  141)    CHARACTER (len=80), PUBLIC :: bgcm_fname
(  142)    
(  143)    
(  144)    !! * Accessibility
(  145)    PUBLIC ini_trc, bgcm_setup, bgcm_ini, extract_boundary_vals, &
(  146)         assign_boundary_vals, update_boundary_vals, bgcm_N_obc
(  147) 
(  148)    
(  149) CONTAINS
(  150) 
(  151) !!DB: Replace trc_lec, trc_trp_lec, trc_trp_ctl ... with 1 routine
(  152) !!One of the things I do in this routine is force the same adv/diff schemes
(  153) !!as are used by the TS module.
(  154) !!Initialized netcdf output file as well
(  155) 
(  156)    SUBROUTINE bgcm_setup
(  157)      USE traadv_ctl
(  158)      USE ldftra_oce
(  159)      USE zdf_oce, ONLY : n_zdfexp
(  160)      USE lib_ncdf
(  161) 
(  162)       !!                     ******************






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 4

(  163)       !! local declarations
(  164)       !! ==================
(  165) 
(  166)       INTEGER ::  ji, status
(  167)       CHARACTER (len=32) :: clname
(  168) 
(  169) !!DB: These are all that are needed for BGCM
(  170)       namelist/nattrc/nwritetrc,lrsttr,nrsttr
(  171)       namelist/natnum/ndttrc,rsc,rtrn,ncortrc,crosster
(  172) !!DB: Modified
(  173)       NAMELIST/namtrcldf/ ahtrc0, trcrat
(  174) 
(  175) !!DB: not yet
(  176) #if defined key_trcdmp
(  177)       NAMELIST/namtrcdmp/ ndmptr, ndmpftr, nmldmptr, sdmptr, bdmptr, hdmptr
(  178) #endif
(  179)       !!----------------------------------------------------------------------
(  180) 
(  181) 
(  182) !!DB: must do this as this routine is not called until step() is called
(  183) !!    Another example of the OPA logic that initializes various params and
(  184) !!    routines at unpredictable areas/times of the code. 
(  185) !!DBG
(  186) !      if(lwp) write(numout,*)'DBG: BGCM init: called tra_adv_ctl'
(  187)       call tra_adv_ctl 
(  188) 
(  189)       IF(lwp) THEN
(  190)          WRITE(numout,*) ' '
(  191)          WRITE(numout,*) ' ROUTINE BGCM_SETUP'
(  192)          WRITE(numout,*) ' **************'
(  193)          WRITE(numout,*) ' '
(  194)          WRITE(numout,*) ' namelist for passive tracer'
(  195)          WRITE(numout,*) ' ***************************'
(  196)          WRITE(numout,*) ' '
(  197)       ENDIF
(  198) 
(  199)       numnat=80
(  200) !!DB: new
(  201)       clname='namelist.BGCM_01'
(  202)       OPEN( numnat, FILE= clname, FORM='formatted', STATUS = 'old')
(  203) 
(  204) !! initialization from namelist file
(  205)       !! ----------------------------------------------
(  206)       !! 1.0 namelist nattrc :
(  207) 
(  208)       nwritetrc = 10
(  209)       lrsttr=.FALSE.
(  210)       nrsttr = 0
(  211)       REWIND(numnat)
(  212)       READ(numnat,nattrc)
(  213) 
(  214)       IF(lwp) THEN
(  215)          WRITE(numout,*) ' '
(  216)          WRITE(numout,*) 'nattrc'
(  217)          WRITE(numout,*) ' '
(  218)          WRITE(numout,*)          &
(  219)             ' frequency of outputs for passive tracers nwritetrc = '    &
(  220)             ,nwritetrc  






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 5

(  221)          WRITE(numout,*) ' restart LOGICAL for passive tr. lrsttr = ',   &
(  222)             &         lrsttr
(  223)          WRITE(numout,*) ' control of time step for p. tr. nrsttr = ',   & 
(  224)             &         nrsttr
(  225)          WRITE(numout,*) ' '
(  226)       ENDIF
(  227) 
(  228)       !! 1.1 namelist natnum :
(  229)       !! ---------------------
(  230)       rsc=1.
(  231)       rtrn=1.e-15
(  232)       ncortrc=1
(  233)       ndttrc=4
(  234)       crosster=.FALSE.
(  235) 
(  236)       REWIND(numnat)
(  237)       READ(numnat,natnum)
(  238) 
(  239) !!Compute the first time step of tracer model
(  240)       nittrc000 = nit000 + ndttrc - 1
(  241) 
(  242)       IF(lwp) THEN
(  243)          WRITE(numout,*) ' '
(  244)          WRITE(numout,*) 'natnum'
(  245)          WRITE(numout,*) ' '
(  246)          WRITE(numout,*) ' tuning coefficient              rsc     = ',    &
(  247)             rsc
(  248)          WRITE(numout,*) ' truncation value                rtrn    = ',    &
(  249)             rtrn
(  250)          WRITE(numout,*) ' number of corrective phase      ncortrc = ',    &
(  251)             ncortrc
(  252)          WRITE(numout,*) ' time step freq. for pass. trac. ndttrc  = ',    &
(  253)             ndttrc
(  254)          WRITE(numout,*) ' 1st time step for pass. trac. nittrc000 = ',    &
(  255)             nittrc000
(  256)          WRITE(numout,*) ' computes or not crossterms    crosster  = ',    &
(  257)             crosster
(  258)       ENDIF
(  259) 
(  260) 
(  261)       !! namelist of transport
(  262)       !! ---------------------
(  263) !!DB: Force tracer code to use the same adv/diff schemes as used for TS
(  264)       ln_trcadv_cen2   = ln_traadv_cen2  
(  265)       ln_trcadv_tvd    = ln_traadv_tvd   
(  266)       ln_trcadv_muscl  = ln_traadv_muscl
(  267)       ln_trcadv_muscl2 = ln_traadv_muscl2 
(  268) 
(  269)       l_trcldf_bilapg  = l_traldf_bilapg 
(  270)       l_trcldf_bilap   = l_traldf_bilap  
(  271)       l_trcldf_iso     = l_traldf_iso    
(  272)       l_trcldf_iso_zps = l_traldf_iso_zps
(  273)       l_trcldf_lap     = l_traldf_lap    
(  274) 
(  275)       l_trczdf_exp     = l_trazdf_exp
(  276)       if(l_trczdf_exp) then
(  277)          l_trczdf_imp = .false.
(  278) !!DBG: This is not debugged, so if using explicit scheme then you should check






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 6

(  279) !!     the following assingment
(  280)          n_trczdf_exp = n_zdfexp         
(  281)       else
(  282)          l_trczdf_imp = .true.
(  283)       endif
(  284) !!DB: I doubt that this will ever be .TRUE.; if so then there may
(  285) !!    be a problem with the above  ----->  refer to original trctrp_ctl.F90
(  286)       l_trczdf_iso_vo  = l_trazdf_iso_vo
(  287) 
(  288)       !!DB: keep this output, at least for a while
(  289)       ! Parameter control and print
(  290)       ! ---------------------------
(  291)       IF(lwp) THEN
(  292)          WRITE(numout,*)
(  293)          WRITE(numout,*) 'DB: choice/control of the tracer advection scheme'
(  294)          WRITE(numout,*) '~~~~~~~~~~~'
(  295)          WRITE(numout,*)
(  296)          WRITE(numout,*) '             2nd order advection scheme     ln_trcadv_cen2   = ', ln_trcadv_cen2
(  297)          WRITE(numout,*) '             TVD advection scheme           ln_trcadv_tvd    = ', ln_trcadv_tvd
(  298)          WRITE(numout,*) '             MUSCL  advection scheme        ln_trcadv_muscl  = ', ln_trcadv_muscl
(  299)          WRITE(numout,*) '             MUSCL2 advection scheme        ln_trcadv_muscl2 = ', ln_trcadv_muscl2
(  300)       ENDIF
(  301) 
(  302)       !  Define the lateral tracer physics parameters
(  303)       ! =============================================
(  304) 
(  305) !!DB: I only need this to get ahtrc0 and trcrat
(  306) !!DB: ahtrc0 is needed #if NOT defined key_traldf_smag, 
(  307) !!    which likely would cause problems elsewhere 
(  308) !!DB: Note that above call to tra_adv_ctl and subsequent assignments
(  309) !!    means that most of this namelist is not needed and has been eliminated.
(  310) !!DB: Also, (eliminated param) ln_trcldf_diff controls whether or not to perform lateral diffusion
(  311) !!    so if you want NO lateral diffusion you must go to trctrp.F90
(  312) ! Read Namelist namtrcldf : Lateral physics on tracers
(  313)       REWIND( numnat )
(  314)       READ  ( numnat, namtrcldf )
(  315) 
(  316)       IF(lwp) THEN
(  317)          WRITE(numout,*)
(  318)          WRITE(numout,*) 'DB: lateral passive tracer physics'
(  319)          WRITE(numout,*) '~~~~~~~'
(  320)          WRITE(numout,*) '   Namelist namtrcldf : set lateral mixing parameters (type, direction, coefficients)'
(  321)          WRITE(numout,*) '     laplacian operator                             ln_trcldf_lap   = ', l_trcldf_lap
(  322)          WRITE(numout,*) '     bilaplacian operator                           ln_trcldf_bilap = ', l_trcldf_bilap
(  323)          WRITE(numout,*) '     lateral eddy diffusivity                              ahtrc0   = ', ahtrc0
(  324)          WRITE(numout,*) '     ratio between passive and active tracer diffusion coef  trcrat = ', trcrat
(  325)       ENDIF
(  326) 
(  327) !!DB: REM, assigned above
(  328)       ! Parameter print
(  329)       ! ---------------
(  330)       IF(lwp) THEN
(  331)          WRITE(numout,*)
(  332)          WRITE(numout,*) 'DB: vertical physics'
(  333)          WRITE(numout,*) '~~~~~~~~'
(  334)          WRITE(numout,*) '          Namelist namtrczdf : set vertical diffusion parameters'
(  335)          WRITE(numout,*) '             time splitting / backward scheme ln_trczdf_exp = ', l_trczdf_exp
(  336)          WRITE(numout,*) '             number of time step               n_trczdf_exp = ', n_trczdf_exp






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 7

(  337)       ENDIF
(  338) 
(  339) !!DB: not accounted for yet
(  340) # if defined key_trcdmp
(  341)       !!DB: I have not covered this key yet so be careful
(  342)       ! Read Namelist namtdp : passive tracers damping term
(  343)       ! --------------------
(  344)       REWIND ( numnat )
(  345)       READ   ( numnat, namtrcdmp )
(  346)       IF( lzoom ) THEN
(  347)          nmldmptr = 0           ! restoring to climatology at closed north or south boundaries
(  348)       ENDIF
(  349) 
(  350)       ! Parameter control and print
(  351)       ! ---------------------------
(  352)       IF(lwp) THEN
(  353)          WRITE(numout,*)
(  354)          WRITE(numout,*) 'newtonian damping'
(  355)          WRITE(numout,*) '~~~~~~~'
(  356)          WRITE(numout,*) '          Namelist namtrcdmp : set damping parameter'
(  357)          WRITE(numout,*)
(  358)          WRITE(numout,*) '             tracers damping option         ndmptr   = ', ndmptr
(  359)          WRITE(numout,*) '             create a damping.coeff file    ndmpftr  = ', ndmpftr
(  360)          WRITE(numout,*) '             mixed layer damping option     nmldmptr = ', nmldmptr, '(zoom: forced to 0)'
(  361)          WRITE(numout,*) '             surface time scale (days)      sdmptr   = ', sdmptr
(  362)          WRITE(numout,*) '             bottom time scale (days)       bdmptr   = ', bdmptr
(  363)          WRITE(numout,*) '             depth of transition (meters)   hdmptr   = ', hdmptr
(  364)          WRITE(numout,*)
(  365)       ENDIF
(  366) 
(  367) #endif
(  368) 
(  369) 
(  370)       !! namelist of SMS
(  371) !!DB
(  372) !!input of params for Source-Minus-Sink (SMS) part of code
(  373) !! ... TO DO ...
(  374) !!DBG
(  375)       if(lwp) write(numout,*)'DBG: bgcm_setup: no SMS setup yet'
(  376) !      call  some_routine_to_be_written(put at bottom of this module)
(  377) 
(  378) !!DB: create ncdf file 
(  379)       bgcm_fname = trim(cexper)//'_bgcm_01.nc'
(  380)       call ncdf_create_bgcm_file(status)
(  381) 
(  382) 
(  383)    END SUBROUTINE bgcm_setup
(  384) 
(  385) 
(  386) !!DB: Initialize model here
(  387)    SUBROUTINE bgcm_ini
(  388)       !!---------------------------------------------------------------------
(  389)       !!              
(  390)       !! ** Purpose : Specific initialization for BGCM_01
(  391)       !!
(  392)       !! * local declarations
(  393) 
(  394)      USE lbclnk






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 8

(  395) 
(  396)      INTEGER ::                   & 
(  397)           ji ,jj ,jk ,jn, jl        ! dummy loop indices  
(  398)      !!---------------------------------------------------------------------
(  399) !!DB -- example of using pointers to make NPZ code more explicit
(  400)      real(wp), pointer :: N(:,:,:), P(:,:,:), Z(:,:,:)
(  401)      real(wp), pointer :: N0(:,:,:), P0(:,:,:), Z0(:,:,:)
(  402)      integer :: i1,i2,j1,j2
(  403)      real(wp) :: c0,c1,s0,s1
(  404) 
(  405) !!DB -- assign pointers to default trn, trb arrays
(  406)      N => trn(:,:,:,1)
(  407)      P => trn(:,:,:,2)
(  408)      Z => trn(:,:,:,3)
(  409)      N0 => trb(:,:,:,1)
(  410)      P0 => trb(:,:,:,2)
(  411)      Z0 => trb(:,:,:,3)
(  412) 
(  413) 
(  414) !!DB: 2008.11.19 -- compute wmask
(  415)      wmask(:,:,1) = 0.0;   wmask(:,:,jpk+1) = 0.0
(  416)      do ji = 1, jpi
(  417)         do jj = 1, jpj
(  418)            do jk = 2, jpk
(  419)               wmask(ji,jj,jk) = min(tmask(ji,jj,jk-1),tmask(ji,jj,jk))
(  420)            enddo
(  421)         enddo
(  422)      enddo
(  423) 
(  424)      
(  425)      !! 1. initialization of passive tracer fields
(  426)      !! -------------------------------------------
(  427)      trn(:,:,:,:)=0.0
(  428)      tra(:,:,:,:)=0.0
(  429) 
(  430) 
(  431) !!DB 2009.06.03 -- try a formula that decays offshelf
(  432) !!DB - from 1D model 
(  433) !!NB:  different profile in deep water
(  434) !     do ji = 1, jpi
(  435) !        do jj = 1, jpj
(  436) !           do jk = 1, jpk
(  437) !              N(ji,jj,jk) = min(5.0 + (10.-5.)/200.0 * gdept(jk),10.0)
(  438) !           enddo
(  439) !           if(mbathy(ji,jj) > 2) then
(  440) !              if(gdept(mbathy(ji,jj)-1) >= 500.0) N(ji,jj,:)=0.1
(  441) !!              if(gdept(mbathy(ji,jj)-1) >= 1500.0) N(ji,jj,:)=0.1
(  442) !           endif
(  443) !        enddo
(  444) !     enddo
(  445) !!DB 2009.06.03 -- try a formula that decays offshelf
(  446) !     c0=5.0
(  447) !     s0=(10.-5.)/250.
(  448) !     do ji = 1, jpi
(  449) !        do jj = 1, jpj
(  450) !           do jk = 1, jpk
(  451) !              N(ji,jj,jk) = min(c0 + s0 * gdept(jk),15.0)
(  452) !           enddo






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 9

(  453) !           if(mbathy(ji,jj) > 2) then
(  454) !              if(gdept(mbathy(ji,jj)-1) >= 500.0) then
(  455) !                 c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/2500.)) + 0.1
(  456) !                 s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/2500.))
(  457) !                 do jk = 1, jpk
(  458) !                    N(ji,jj,jk) = min(c1 + s1 * gdept(jk),15.0)
(  459) !                 enddo
(  460) !              endif
(  461) !           endif
(  462) !        enddo
(  463) !     enddo
(  464) 
(  465) !!DB 2009.06.08 -- try a linear N(z) formula that decays offshelf
(  466) !!with time-dependent params, which will also be applied as OBCs
(  467) !!NB: lots of HARDWIRED stuff at this time
(  468) !06.10: started with 2500 as the decay scale which was too broad
(  469) !replaced 2500 with 1500
(  470) !c0 == N(z=0) = c1+c2*cos(pi*((yd-c3)/180));
(  471) !c1=3.0;c2=2.25;c3=25.0;  
(  472)      c0 = 3.0 + 2.25*cos((nday_year-25.0)*rpi/180.)
(  473)      s0=(15.-c0)/150.
(  474)      do ji = 1, jpi
(  475)         do jj = 1, jpj
(  476)            do jk = 1, jpk
(  477)               N(ji,jj,jk) = min(c0 + s0 * gdept(jk),15.0)
(  478)            enddo
(  479)            if(mbathy(ji,jj) > 2) then
(  480)               if(gdept(mbathy(ji,jj)-1) >= 500.0) then
(  481)                  c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.)) + 0.1
(  482)                  s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.))
(  483)                  do jk = 1, jpk
(  484)                     N(ji,jj,jk) = min(c1 + s1 * gdept(jk),15.0)
(  485)                  enddo
(  486)               endif
(  487)            endif
(  488)         enddo
(  489)      enddo
(  490) !!DBG
(  491)      if(lwp) write(221+narea,*)'DBG: BGCM initial N(z=0), yd: ',c0,nday_year
(  492) 
(  493) 
(  494)     P(:,:,:) = 1.0e-2
(  495)     Z(:,:,:) = 1.0e-2
(  496) 
(  497) !!DBG
(  498) !    P(:,:,:) = 5.0
(  499) !    Z(:,:,:) = 2.0
(  500) !    N(:,:,:) = 10.0
(  501) 
(  502) 
(  503) !!To be safe
(  504)     do jn = 1, jptra
(  505)        CALL lbc_lnk( trn(:,:,:,jn), 'T', 1. )   
(  506)     enddo
(  507) 
(  508) !!DB: extract boundary values 
(  509) !!The below retained for possible future use
(  510)     call extract_boundary_vals






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 10

(  511) !! set before field
(  512)     e_trb(:,:,:) = e_trn(:,:,:)
(  513)     w_trb(:,:,:) = w_trn(:,:,:)
(  514)     s_trb(:,:,:) = s_trn(:,:,:)
(  515) 
(  516) 
(  517) !!DB: If some BGCM data are to be read from a file, do it here. 
(  518) !#if defined key_trc_dta
(  519) !!   Initialization of tracer from a file
(  520) !      CALL dta_trc( nit000 )
(  521) !      DO  jk = 1, jptra
(  522) !        IF( lutini(jk) ) THEN 
(  523) !           trn(:,:,:,jk) = trdta(:,:,:,jk)*tmask(:,:,:)
(  524) !        ENDIF
(  525) !      END DO
(  526) !#endif
(  527) 
(  528) !! before field :
(  529) !! -------------
(  530)      trb(:,:,:,:) = trn(:,:,:,:)
(  531) 
(  532)      if( lwp ) then
(  533)         write(numout,*) ' '
(  534)         write(numout,*) 'BGCM_01 initialisation done '
(  535) !!DBG
(  536) !        write(numout,*) 'DBG: Fields constant everywhere '
(  537)         write(numout,*) ' '
(  538)      endif
(  539)       
(  540)  END SUBROUTINE bgcm_ini
(  541) 
(  542) 
(  543) !!DB: extract boundary values from trn and store them in 
(  544) !!the ?_trn boundary arrays
(  545) SUBROUTINE extract_boundary_vals
(  546)   USE obc_oce
(  547) #  include "obc_vectopt_loop_substitute.h90"
(  548)   IMPLICIT NONE
(  549)   INTEGER :: ji,jj,jk,jn
(  550) 
(  551)   do ji = fs_nie0+1, fs_nie1+1 ! isolate processor
(  552)      e_trn(:,:,:) = trn(ji,:,:,:)
(  553)   enddo
(  554)   do ji = niw0, niw1 ! isolate processor
(  555)      w_trn(:,:,:) = trn(ji,:,:,:)
(  556)   enddo
(  557)   do jj = njs0, njs1
(  558)      s_trn(:,:,:) = trn(:,jj,:,:)
(  559)   enddo
(  560) END SUBROUTINE extract_boundary_vals
(  561) 
(  562) !!DB: Put the ?_trn boundary arrays values into tra at the boundary indices
(  563) !!Routine retained for possible future use
(  564) !!For example, if using N data along Open Boundaries then a version of this
(  565) !!might be called
(  566) SUBROUTINE assign_boundary_vals
(  567)   USE obc_oce
(  568) #  include "obc_vectopt_loop_substitute.h90"






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 11

(  569)   IMPLICIT NONE
(  570)   INTEGER :: ji,jj,jk,jn
(  571) 
(  572)   do ji = fs_nie0+1, fs_nie1+1 !as in obctra
(  573)      tra(ji,:,:,:) = e_trn(:,:,:)
(  574)   enddo
(  575)   do ji = niw0, niw1 ! isolate processor
(  576)      tra(ji,:,:,:) = w_trn(:,:,:) 
(  577)   enddo
(  578)   do jj = njs0, njs1
(  579)      tra(:,jj,:,:) = s_trn(:,:,:) 
(  580)   enddo
(  581) 
(  582) END SUBROUTINE assign_boundary_vals
(  583) 
(  584) !!DB: 2009.04.20 -- copy values 1 cell in to OB positions
(  585) !!use "after" arrays
(  586) !!This routine is effectively the OBC code. It is called in trcnxt.F90
(  587) SUBROUTINE update_boundary_vals
(  588)   USE obc_oce
(  589) #  include "obc_vectopt_loop_substitute.h90"
(  590)   IMPLICIT NONE
(  591)   INTEGER :: ji,jj,jk,jn
(  592) 
(  593) !!DBG: assign to outer(inner) values as well
(  594)   do ji = fs_nie0+1, fs_nie1+1 !as in obctra
(  595)      tra(ji,:,:,:) = tra(ji-1,:,:,:) 
(  596)   enddo
(  597)   do ji = niw0, niw1 ! isolate processor
(  598)      tra(ji,:,:,:) = tra(ji+1,:,:,:) 
(  599)   enddo
(  600)   do jj = njs0, njs1
(  601)      tra(:,jj,:,:) = tra(:,jj+1,:,:)
(  602)   enddo
(  603) 
(  604) 
(  605) END SUBROUTINE update_boundary_vals
(  606) 
(  607) !!2009.08
(  608) !!Update N OBC vals and assign to tra()
(  609) !!To be called in trcnxt()
(  610)   SUBROUTINE bgcm_N_obc ( kt )
(  611)     USE obc_oce         !open boundary condition variables
(  612)     USE lbclnk
(  613) #  include "obc_vectopt_loop_substitute.h90"
(  614)     
(  615)     !!------------------------------------------------------------------------------
(  616)     !! * Arguments
(  617)     INTEGER, INTENT( in ) ::   kt
(  618)     
(  619)     !! * Local declaration
(  620)     INTEGER ::   ji, jj, jk, jn
(  621)     REAL :: c0,c1,s0,s1
(  622) 
(  623) 
(  624) !!NB: eventually do this once/day
(  625) !!DB 2009.06.08 -- try a linear N(z) formula that decays offshelf
(  626) !!with time-dependent params, which will also be applied as OBCs






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 12

(  627) !!NB: lots of HARDWIRED stuff at this time
(  628) !06.10: started with 2500 as the decay scale which was too broad
(  629) !replaced 2500 with 1500
(  630) !c0 == N(z=0) = c1+c2*cos(pi*((yd-c3)/180));
(  631) !c1=3.0;c2=2.25;c3=25.0;  
(  632)      c0 = 3.0 + 2.25*cos((nday_year-25.0)*rpi/180.)
(  633)      s0=(15.-c0)/150.
(  634) !!isolate processors for each OB
(  635)      do ji = fs_nie0, fs_nie1
(  636)         do jj = 1, jpj
(  637)            do jk = 1, jpk
(  638)               e_trn(jj,jk,1) = min(c0 + s0 * gdept(jk),15.0)
(  639)            enddo
(  640)            if(mbathy(ji,jj) > 2) then
(  641)               if(gdept(mbathy(ji,jj)-1) >= 500.0) then
(  642)                  c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.)) + 0.1
(  643)                  s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.))
(  644)                  do jk = 1, jpk
(  645)                     e_trn(jj,jk,1) = min(c1 + s1 * gdept(jk),15.0)
(  646)                  enddo
(  647)               endif
(  648)            endif
(  649)         enddo
(  650)         tra(ji+1,:,:,1) = e_trn(:,:,1)
(  651)      enddo
(  652)      do ji = niw0, niw1
(  653)         do jj = 1, jpj
(  654)            do jk = 1, jpk
(  655)               w_trn(jj,jk,1) = min(c0 + s0 * gdept(jk),15.0)
(  656)            enddo
(  657)            if(mbathy(ji,jj) > 2) then
(  658)               if(gdept(mbathy(ji,jj)-1) >= 500.0) then
(  659)                  c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.)) + 0.1
(  660)                  s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.))
(  661)                  do jk = 1, jpk
(  662)                     w_trn(jj,jk,1) = min(c1 + s1 * gdept(jk),15.0)
(  663)                  enddo
(  664)               endif
(  665)            endif
(  666)         enddo
(  667)         tra(ji,:,:,1) = w_trn(:,:,1)
(  668)      enddo
(  669)      do jj = njs0, njs1
(  670)         do ji = 1, jpi
(  671)            do jk = 1, jpk
(  672)               s_trn(ji,jk,1) = min(c0 + s0 * gdept(jk),15.0)
(  673)            enddo
(  674)            if(mbathy(ji,jj) > 2) then
(  675)               if(gdept(mbathy(ji,jj)-1) >= 500.0) then
(  676)                  c1 = c0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.)) + 0.1
(  677)                  s1 = s0*(1.-tanh(gdept(mbathy(ji,jj)-1)/1500.))
(  678)                  do jk = 1, jpk
(  679)                     s_trn(ji,jk,1) = min(c1 + s1 * gdept(jk),15.0)
(  680)                  enddo
(  681)               endif
(  682)            endif
(  683)         enddo
(  684)         tra(:,jj,:,1) = s_trn(:,:,1)






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 13

(  685)      enddo
(  686) !!DBG
(  687)      if(lwp) write(221+narea,*)'DBG: BGCM kt N(z=0), yd: ',kt, c0,nday_year
(  688) 
(  689) !!DBG: every day (hardwired) output a few vertical profiles to monitor output
(  690)      if(mod(kt-1,180) == 0) then
(  691)         
(  692) !!West
(  693)         do ji = niw0, niw1
(  694)            do jj = 2, jpj, int(jpj/4)
(  695)               do jk = 1, jpk
(  696)                  write(5000+narea,'(3(i8,1x),10(e10.4,2x))')kt,jj,jk, (w_trn(jj,jk,jn),jn=1,jptra), &
(  697)                       wmask(ji,jj,jk), gdept(mbathy(ji,jj)-1),float(mbathy(ji,jj))
(  698)               enddo
(  699)            enddo
(  700)         enddo
(  701) !!South
(  702)         do jj = fs_njs0, fs_njs1
(  703)            do ji = 2, jpi, int(jpi/4)
(  704)               do jk = 1, jpk
(  705)                  write(6000+narea,'(3(i8,1x),10(e10.4,2x))')kt,ji,jk, (s_trn(ji,jk,jn),jn=1,jptra), &
(  706)                       wmask(ji,jj,jk), gdept(mbathy(ji,jj)-1),float(mbathy(ji,jj))
(  707)               enddo
(  708)            enddo
(  709)         enddo
(  710) !!East
(  711)         do ji = fs_nie0, fs_nie1
(  712)            do jj = 2, jpj, int(jpj/4)
(  713)               do jk = 1, jpk
(  714)                  write(7000+narea,'(3(i8,1x),10(e10.4,2x))')kt,jj,jk, (e_trn(jj,jk,jn),jn=1,jptra), &
(  715)                       wmask(ji,jj,jk), gdept(mbathy(ji,jj)-1),float(mbathy(ji,jj))
(  716)               enddo
(  717)            enddo
(  718)         enddo
(  719)            
(  720)      endif
(  721) 
(  722) 
(  723)   END SUBROUTINE bgcm_N_obc
(  724) 
(  725) !!DB: 2008.08.28 ...
(  726) !!Create output file for BGCM_01
(  727) !!The main variable is the trn(:,:,:,jptra) array which is what is integrated
(  728) !!by the tracer module. 
(  729) !!N,P,Z can be used instead (see below and bgcm_01_model.F90) but using
(  730) !!both takes up double the disk space (even if no output is done),
(  731) !!so define and use only one of these possibilities.
(  732)  SUBROUTINE ncdf_create_bgcm_file(status)
(  733) 
(  734)    USE lib_ncdf
(  735)    USE trc
(  736) 
(  737)     IMPLICIT NONE
(  738)     ! Subroutine argument declarations
(  739) 
(  740)     INTEGER,INTENT(OUT) :: status
(  741) 
(  742)     ! Local declarations






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 14

(  743)     INTEGER :: ncid,    &  ! netCDF file ID
(  744)                varid,   &  ! ID of netCDF variable to be written to
(  745)                nfstat,  &  ! netCDF library call return status
(  746)                mpistat     ! MPI library call return status
(  747)     INTEGER,DIMENSION(1:5) :: dimids
(  748) !!DB -- unsure how to dimension varids, so choose a large number 
(  749)     INTEGER,DIMENSION(1:50) :: varids
(  750)     CHARACTER(LEN=20) :: cal_type      ! Calendar type
(  751)     CHARACTER(LEN=30) :: timestamp     ! File timestamp
(  752)     CHARACTER(LEN=100) :: sec_since    
(  753)     CHARACTER(LEN=100) :: t_origin     ! Time origin of this run
(  754)     CHARACTER(LEN=20) :: op_type, varname
(  755)     INTEGER :: int_opp, &              ! Operation interval
(  756)                int_wri                 ! Write interval
(  757)     CHARACTER(LEN=3),PARAMETER :: &
(  758)          &  months(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
(  759)          &                 'JUL','AUG','SEP','OCT','NOV','DEC'/)
(  760)     INTEGER :: jf, varnum
(  761) !!DB NEW
(  762)     REAL(wp), DIMENSION(jpk) :: param_val = 0.0
(  763) 
(  764)     
(  765)     ! Initializations
(  766)     op_type = 'instantaneous'       !!default
(  767) !    op_type = 'ave(x)'       !!requires routine to do this  
(  768)     op_type = TRIM(op_type)
(  769) 
(  770) 
(  771)     status = NCDF_NOERR
(  772)     CALL ioget_calendar(cal_type)
(  773)     CALL ioget_timestamp(timestamp)
(  774)     WRITE (UNIT=sec_since, &
(  775)          FMT='("seconds since ",I4.4,2("-",I2.2)," ",I2.2,2(":",I2.2))') &
(  776)          &  nyear,nmonth,nday,0, 0, 0
(  777)     WRITE(t_origin, &
(  778)          &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)") &
(  779)          &   nyear,months(nmonth),nday,0,0,0
(  780) 
(  781) !!DB
(  782)     int_opp = nwrite * rdt
(  783)     int_wri = nwrite * rdt
(  784) 
(  785)     IF((DEBUG_OUT .EQV. .TRUE.) .AND. (nproc == 0)) THEN
(  786)        WRITE(100,*) 'NCDF DEBUG: Creating output file:', bgcm_fname
(  787)        CALL FLUSH
(  788)     END IF
(  789) 
(  790)     ! Only processor 0 does anything
(  791)     IF(nproc == 0) THEN
(  792)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  793)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Creating file:', bgcm_fname
(  794)           CALL FLUSH
(  795)        END IF
(  796)        ! Create the file
(  797)        nfstat = nf90_create(bgcm_fname, nf90_clobber, ncid)
(  798)        IF(nfstat /= nf90_noerr) THEN
(  799)           status = NCDF_NFERR
(  800)           RETURN






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 15

(  801)        END IF
(  802)        
(  803)        ! Define dimensions
(  804)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  805)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Defining dimensions in file:', bgcm_fname
(  806)           CALL FLUSH
(  807)        END IF
(  808)        !nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  809)        !nfstat = nf90_def_dim(ncid, 'jpkdta', jpkdta, dimids(2))
(  810)        !nfstat = nf90_def_dim(ncid, 'jpjdta', jpjdta, dimids(3))
(  811)        !nfstat = nf90_def_dim(ncid, 'jpidta', jpidta, dimids(4))
(  812)        nfstat = nf90_def_dim(ncid, 'time_counter', nf90_unlimited, dimids(1))
(  813)        nfstat = nf90_def_dim(ncid, 'jptra', jptra, dimids(2))  !!CN: Changed dim creation order
(  814)        nfstat = nf90_def_dim(ncid, 'z', jpk, dimids(3))          !!bgcm model
(  815)        nfstat = nf90_def_dim(ncid, 'y', jpjdta, dimids(4))
(  816)        nfstat = nf90_def_dim(ncid, 'x', jpidta, dimids(5))
(  817)        !nfstat = nf90_def_dim(ncid, 'jptra', jptra, dimids(5)) !!CN: Changed dim creation order
(  818)        
(  819)        ! Define variables
(  820)        ! CN: changed order of dims in dim lists for var creation
(  821)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  822)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Defining variables in file:', bgcm_fname
(  823)           CALL FLUSH
(  824)        END IF
(  825)        nfstat = nf90_def_var(ncid, 'nav_lon', nf90_float, &
(  826)             (/ dimids(5), dimids(4) /), &
(  827)             varids(1))
(  828)        nfstat = nf90_def_var(ncid, 'nav_lat', nf90_float, &
(  829)             (/ dimids(5), dimids(4) /), &
(  830)             varids(2))
(  831)        nfstat = nf90_def_var(ncid, 'deptht', nf90_float, &
(  832)             (/ dimids(3) /), &
(  833)             varids(3))
(  834)        nfstat = nf90_def_var(ncid, 'time_counter', nf90_float, &
(  835)             (/ dimids(1) /), &
(  836)             varids(4))
(  837) !!DBG: Choose between NPZ or trn variable declarations 
(  838) !!NB: must be consistent with output in bgcm_01_model
(  839)        jf = 1
(  840)        if(jf == 0) then
(  841)           nfstat = nf90_def_var(ncid, 'trn', nf90_float, &
(  842)                (/ dimids(5), dimids(4), dimids(3), dimids(2), dimids(1) /), &
(  843)                varids(5))
(  844)           varnum = 5
(  845)        else
(  846)           nfstat = nf90_def_var(ncid, 'N', nf90_float, &
(  847)                (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
(  848)                varids(6))
(  849)           nfstat = nf90_def_var(ncid, 'P', nf90_float, &
(  850)                (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
(  851)                varids(7))
(  852)           nfstat = nf90_def_var(ncid, 'Z', nf90_float, &
(  853)                (/ dimids(5), dimids(4), dimids(3), dimids(1) /), &
(  854)                varids(8))
(  855)           varnum = 8
(  856)        endif
(  857) 
(  858) !AD/DB: add new time-related variables






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 16

(  859)        varnum = varnum + 1
(  860)        ! ndate (ndastp)
(  861)        nfstat = nf90_def_var(ncid, 'ndastp', nf90_float, &
(  862)             (/ dimids(1) /), &
(  863)             varids(varnum))
(  864)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
(  865)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'time step date in year/month/day aammjj')
(  866) 
(  867)        varnum = varnum + 1
(  868)        ! ndate (model_time)
(  869)        nfstat = nf90_def_var(ncid, 'model_time', nf90_float, &
(  870)             (/ dimids(1) /), &
(  871)             varids(varnum))
(  872)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', &
(  873)             'time step date (when output is writen) in year/month/day aammjj (decimal day)')
(  874)        nfstat = nf90_put_att(ncid, varids(varnum), 'units', '=nyear*10000+nmonth*100+nday')
(  875)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula1', 'nyear  =   model_time / 10000')       
(  876)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula2', & 
(  877)             'nmonth = ( pmodel_time - (nyear * 10000) ) / 100')       
(  878)        nfstat = nf90_put_att(ncid, varids(varnum), 'formula3', & 
(  879)             'nday   =   model_time - (nyear * 10000) - ( nmonth * 100 )')                           
(  880) 
(  881)        varnum = varnum + 1
(  882)        ! kt 
(  883)        nfstat = nf90_def_var(ncid, 'model_time_step', nf90_float, &
(  884)             (/ dimids(1) /),  varids(varnum))
(  885) 
(  886) !!DB: HARDWIRED: add model params to .nc file
(  887) !!Use z=jpk as dimension so no new dimension need be defined.
(  888) !!If jpk becomes too small, use x or y
(  889) !!Use long_name for param names (other attribute names did not work)
(  890)        varnum = varnum + 1
(  891)        nfstat = nf90_def_var(ncid, 'model_params', nf90_float, &
(  892)             (/ dimids(3) /), varids(varnum))
(  893)        nfstat = nf90_put_att(ncid, varids(varnum), 'long_name', 'Hardwired -- 11 params: &
(  894)             Vm,Rm,ks,lambda,gamma,M_ZOO,M_PHY,kext,aaa,PAR,w_sink')
(  895)        nfstat = nf90_put_att(ncid, varids(varnum), 'short_name', 'vals in model_params(1:11)')
(  896)        
(  897)        ! Add attributes
(  898)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  899)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Writing attributes in file:', bgcm_fname
(  900)           CALL FLUSH
(  901)        END IF
(  902)        ! nav_lon
(  903)        nfstat = nf90_put_att(ncid, varids(1), 'units', 'degrees_east')
(  904)        nfstat = nf90_put_att(ncid, varids(1), 'valid_min', minval(glamt))
(  905)        nfstat = nf90_put_att(ncid, varids(1), 'valid_max', maxval(glamt))
(  906)        nfstat = nf90_put_att(ncid, varids(1), 'long_name', 'Longitude')
(  907)        nfstat = nf90_put_att(ncid, varids(1), 'nav_model', 'Default grid')
(  908) 
(  909)        ! nav_lat
(  910)        nfstat = nf90_put_att(ncid, varids(2), 'units', 'degrees_north')
(  911)        nfstat = nf90_put_att(ncid, varids(2), 'valid_min', minval(gphit))
(  912)        nfstat = nf90_put_att(ncid, varids(2), 'valid_max', maxval(gphit))
(  913)        nfstat = nf90_put_att(ncid, varids(2), 'long_name', 'Latitude')
(  914)        nfstat = nf90_put_att(ncid, varids(2), 'nav_model', 'Default grid')
(  915) 
(  916)        ! deptht






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 17

(  917)        nfstat = nf90_put_att(ncid, varids(3), 'units', 'm')
(  918)        nfstat = nf90_put_att(ncid, varids(3), 'positive', 'unknown')
(  919)        nfstat = nf90_put_att(ncid, varids(3), 'valid_min', minval(gdept))      !for bgcm model
(  920)        nfstat = nf90_put_att(ncid, varids(3), 'valid_max', maxval(gdept))      !for bgcm model
(  921)        nfstat = nf90_put_att(ncid, varids(3), 'title', 'deptht')
(  922)        nfstat = nf90_put_att(ncid, varids(3), 'long_name', 'Vertical Tracer levels')
(  923) 
(  924) !!DB: temporary attributes (REM: jf set above)
(  925)        ! trn 
(  926)        if(jf == 0) then
(  927)           nfstat = nf90_put_att(ncid, varids(5), 'units', 'XXX')
(  928)           nfstat = nf90_put_att(ncid, varids(5), 'missing_value', 1.000000E20)
(  929)           nfstat = nf90_put_att(ncid, varids(5), 'valid_min', 1.000000E20 )
(  930)           nfstat = nf90_put_att(ncid, varids(5), 'valid_max', -1.000000E20)
(  931)           nfstat = nf90_put_att(ncid, varids(5), 'long_name', 'BGCM_01 tracers N-P-Z')
(  932)           nfstat = nf90_put_att(ncid, varids(5), 'short_name', 'trn')
(  933)           nfstat = nf90_put_att(ncid, varids(5), 'online_operation', op_type)
(  934)           nfstat = nf90_put_att(ncid, varids(5), 'axis', 'TNZYX')
(  935)           nfstat = nf90_put_att(ncid, varids(5), 'interval_operation', float(int_opp))
(  936)           nfstat = nf90_put_att(ncid, varids(5), 'interval_write', float(int_wri))
(  937)           nfstat = nf90_put_att(ncid, varids(5), 'associate', 'time_counter jptra deptht nav_lat nav_lon')
(  938)        else
(  939)        ! N
(  940)           nfstat = nf90_put_att(ncid, varids(6), 'units', 'XXX')
(  941)           nfstat = nf90_put_att(ncid, varids(6), 'missing_value', 1.000000E20)
(  942)           nfstat = nf90_put_att(ncid, varids(6), 'valid_min', 1.000000E20 )
(  943)           nfstat = nf90_put_att(ncid, varids(6), 'valid_max', -1.000000E20)
(  944)           nfstat = nf90_put_att(ncid, varids(6), 'long_name', 'BGCM_01-- N')
(  945)           nfstat = nf90_put_att(ncid, varids(6), 'short_name', 'N')
(  946)           nfstat = nf90_put_att(ncid, varids(6), 'online_operation', TRIM(op_type))
(  947)           nfstat = nf90_put_att(ncid, varids(6), 'axis', 'TZYX')
(  948)           nfstat = nf90_put_att(ncid, varids(6), 'interval_operation', float(int_opp))
(  949)           nfstat = nf90_put_att(ncid, varids(6), 'interval_write', float(int_wri))
(  950)           nfstat = nf90_put_att(ncid, varids(6), 'associate', 'time_counter deptht nav_lat nav_lon')
(  951)           !!DB: 
(  952)           ! P
(  953)           nfstat = nf90_put_att(ncid, varids(7), 'units', 'XXX')
(  954)           nfstat = nf90_put_att(ncid, varids(7), 'missing_value', 1.000000E20)
(  955)           nfstat = nf90_put_att(ncid, varids(7), 'valid_min', 1.000000E20 )
(  956)           nfstat = nf90_put_att(ncid, varids(7), 'valid_max', -1.000000E20)
(  957)           nfstat = nf90_put_att(ncid, varids(7), 'long_name', 'BGCM_01-- P')
(  958)           nfstat = nf90_put_att(ncid, varids(7), 'short_name', 'P')
(  959)           nfstat = nf90_put_att(ncid, varids(7), 'online_operation', TRIM(op_type))
(  960)           nfstat = nf90_put_att(ncid, varids(7), 'axis', 'TZYX')
(  961)           nfstat = nf90_put_att(ncid, varids(7), 'interval_operation', float(int_opp))
(  962)           nfstat = nf90_put_att(ncid, varids(7), 'interval_write', float(int_wri))
(  963)           nfstat = nf90_put_att(ncid, varids(7), 'associate', 'time_counter deptht nav_lat nav_lon')
(  964)           !!DB: 
(  965)           ! Z
(  966)           nfstat = nf90_put_att(ncid, varids(8), 'units', 'XXX')
(  967)           nfstat = nf90_put_att(ncid, varids(8), 'missing_value', 1.000000E20)
(  968)           nfstat = nf90_put_att(ncid, varids(8), 'valid_min', 1.000000E20 )
(  969)           nfstat = nf90_put_att(ncid, varids(8), 'valid_max', -1.000000E20)
(  970)           nfstat = nf90_put_att(ncid, varids(8), 'long_name', 'BGCM_01-- Z')
(  971)           nfstat = nf90_put_att(ncid, varids(8), 'short_name', 'Z')
(  972)           nfstat = nf90_put_att(ncid, varids(8), 'online_operation', TRIM(op_type))
(  973)           nfstat = nf90_put_att(ncid, varids(8), 'axis', 'TZYX')
(  974)           nfstat = nf90_put_att(ncid, varids(8), 'interval_operation', float(int_opp))






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 18

(  975)           nfstat = nf90_put_att(ncid, varids(8), 'interval_write', float(int_wri))
(  976)           nfstat = nf90_put_att(ncid, varids(8), 'associate', 'time_counter deptht nav_lat nav_lon')
(  977)        endif
(  978) 
(  979)        ! time_counter
(  980)        nfstat = nf90_put_att(ncid, varids(4), 'units', TRIM(sec_since))
(  981)        nfstat = nf90_put_att(ncid, varids(4), 'calendar', TRIM(cal_type))
(  982)        nfstat = nf90_put_att(ncid, varids(4), 'title', 'Time')
(  983)        nfstat = nf90_put_att(ncid, varids(4), 'long_name', 'time axis')
(  984)        nfstat = nf90_put_att(ncid, varids(4), 'time_origin', TRIM(t_origin))
(  985) 
(  986)        ! global
(  987)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions', 'GDT 1.3')
(  988)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'file_name', TRIM(bgcm_fname))
(  989)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'TimeStamp', TRIM(timestamp))
(  990)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'associate_file', 'none')
(  991)        nfstat = nf90_put_att(ncid, NF90_GLOBAL, 'MISC', 'DB-created file')
(  992) 
(  993) 
(  994)        ! Close file
(  995)        IF(DEBUG_OUT .EQV. .TRUE.) THEN
(  996)           WRITE(100,*) 'NCDF DEBUG: ncdf_create_bgcm_file - Closing file:', bgcm_fname
(  997)           CALL FLUSH
(  998)        END IF
(  999)        nfstat = nf90_close(ncid)
( 1000)        IF(nfstat /= nf90_noerr) THEN
( 1001)           status = NCDF_NFERR
( 1002)           RETURN
( 1003)        END IF
( 1004)     END IF
( 1005) 
( 1006) 
( 1007) !!Write grid info to file
( 1008)     CALL ncdf_write(bgcm_fname, 'nav_lat', gphit, -1, status)
( 1009)     CALL ncdf_write(bgcm_fname, 'nav_lon', glamt, -1, status)
( 1010)     CALL ncdf_write(bgcm_fname, 'deptht', gdept, status)
( 1011) !!Write param_vals to model_params
( 1012)     param_val(1) = Vm
( 1013)     param_val(2) = Rm
( 1014)     param_val(3) = ks
( 1015)     param_val(4) = lambda
( 1016)     param_val(5) = gamma
( 1017)     param_val(6) = M_ZOO
( 1018)     param_val(7) = M_PHY
( 1019)     param_val(8) = kext
( 1020)     param_val(9) = aaa
( 1021)     param_val(10) = PAR
( 1022)     param_val(11) = w_sink
( 1023)     CALL ncdf_write(bgcm_fname, 'model_params', param_val, status)
( 1024) 
( 1025) 
( 1026)     ! Sync up processors before returning from subroutine
( 1027)     CALL MPI_BARRIER(MPI_COMM_WORLD, mpistat)
( 1028)     IF(mpistat /= 0) THEN
( 1029)        status = NCDF_MPERR
( 1030)        RETURN
( 1031)     END IF
( 1032)     






PGF90 (Version     10.5)          09/23/2010  09:49:17      page 19

( 1033)   END SUBROUTINE ncdf_create_bgcm_file
( 1034) 
( 1035)    
( 1036)    SUBROUTINE ini_trc
( 1037)       !!---------------------------------------------------------------------
( 1038)       !!
( 1039)       !!                       ROUTINE ini_trc
( 1040)       !!                     ******************
( 1041)       !!
( 1042)       !!  PURPOSE :  initialize the BGCM_01 model
( 1043)       !!  DB: Replaces regular ini_trc when key_BGCM_01 is defined
( 1044)       !!---------------------------------------------------------------------
( 1045) 
( 1046) 
( 1047)       !! 0.b PRINT the number of tracer
( 1048)       !! ------------------------------
( 1049) 
( 1050)       IF(lwp) WRITE(numout,*) ' '
( 1051)       IF(lwp) WRITE(numout,*) ' *** number of passive tracer jptra = ',jptra
( 1052)       IF(lwp) WRITE(numout,*) ' '
( 1053) 
( 1054)       call bgcm_setup 
( 1055) 
( 1056)       call bgcm_ini
( 1057) 
( 1058) !!If restarting (not thought about yet) assume that the following call
( 1059) !!will overwrite anything wrong that bgcm_ini might have done
( 1060)       if( lrsttr ) THEN
( 1061) !         CALL bgcm_rst       !!...TO DO ...
( 1062)       endif
( 1063) 
( 1064)    END SUBROUTINE ini_trc
( 1065) 
( 1066) END MODULE lib_bgcm_01
( 1067) 
( 1068) #endif 




























PGF90 (Version     10.5)          09/23/2010  09:49:17      page 20

(    1) # 1 "lib_bgcm_01.F90"
(    1) !!DB 2009.10.09 -- code restructured to put most routines in
(    2) !! lib_bgcm_01.F90
(    3) 
(    4) 
(    5) 
(    6) !!DB 2008.08.18 ... 
(    7) 
(    8) !!This module does all the basic initialization for BGCM_01
(    9) !!It replaces the complicated series of calls to various *trc* subroutines
(   10) !!The calling sequence for this module is:
(   11) !! opa.F90: call ini_trc ----> these routines.
(   12) !!To keep some backward compatibility with generic tracer and biology code
(   13) !!this module is USE'd by initrc.F90 only:
(   14) !!                  #ifdef key_BGCM_01 (see initrc.F90)
(   15) 
(   16) !!NB: this is only used if key_BGCM_01 is defined
(   17) !!   (and key_passivetrc is defined)
(   18) 
(   19) !!Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   20) !!-1- I force the tracer model to use the same diffusion scheme as the dynamic model. 
(   21) !!    I found that the TVD advection scheme for tracers does not work, while the
(   22) !!    other 3 (cen2, muscl, muscl2) do. I hardwire the cen2 scheme which seems 
(   23) !!    slightly faster (see trctrp.F90)
(   24) !!-2- Related to -1- is the variable trcrat which is read from the
(   25) !!    namelist.passivetrc file. One must be careful with this variable. 
(   26) !!    This variable controls the tracer diffusivity:
(   27) !!    From passivetrc_substitute.h90:
(   28) !!    #elif defined 1
(   29) !!      SMAG scheme                    aht: 3D coefficient
(   30) !!    #  define fsahtrt(i,j,k)  trcrat * ahtt(i,j,k)
(   31) !!    DB's default for TS diffusion is the Smagorinsky scheme, in which a Pr-# is used that
(   32) !!    scales Ah versus Am: i.e. Ah(:,:,:) = Pr*Am(:,:,:), where Pr = 0.1 as
(   33) !!    recommended by Z.Wang. Therefore the assignment: trcrat*ahtt really is
(   34) !!    trcrat*Pr*Am. Currently trcrat = 15 -----> A_tracer = 1.5*A_momentum, which
(   35) !!    is considered reasonable to start. 
(   36) !!-3- I do not perform "consistency" controls (e.g. call *_ctl()). Some consistency is
(   37) !!    achieved by virtue of -1- which effectively uses the consistency from TS adv/diff.
(   38) !!-4- Some subroutines do not yet exist -- they are noted, or are obvious. 
(   39) !!-5- Some keys are verbotten -- they are noted, or are obvious. 
(   40) !!
(   41) 
(   42) !!END Notes: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
(   43) 
(   44) !!To see some of the older calls and structure see (e.g.)
(   45) !!(drakes:dbrick) OLD_OPA_CODE_VERSIONS/
(   46) 
PGF90-W-0006-Input file empty (lib_bgcm_01.F90)
