


PGF90 (Version     10.5)          09/23/2010  09:49:39      page 1

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

Filename: obcdta.F90

(    1) !!DB/CN: 2008.11.14
(    2) !!Calls to IOIPSL eliminated and tested
(    3) !!More code cleanup likely possible, but not done
(    4) 
(    5) !!Key changes from standard obcdta.F90
(    6) !!(0) Look for  DB  to see various modifications
(    7) 
(    8) !!(1) obc_dta() has HARDWIRED areas where transports at certain boundary
(    9) !!locations are adjusted ---> key_ADJ_TRANSPORT; SBI_trans; SB_east; ...
(   10) !!- Note that at this time the code is not completely consistent with
(   11) !!respect to key_ADJ_TRANSPORT as some variables and operations exist that
(   12) !!are only required for that key but they are done even if that key is not defined.
(   13) !!These exceptions require minimal memory and execution time.
(   14) 
(   15) !!(2) Because OBC transports are adjusted, the eta-arrays need to be adjusted as
(   16) !!well to be consistent with the OBC vels (thermal wind is used). This is a major
(   17) !!change that requires new variables (search eta_). 
(   18) !!Note that 1 & 2 require definitions and computations on global variables defined
(   19) !!for these new OBC modifications. By global I mean that all relevant CPUs must
(   20) !!do the same computations and have the same values.
(   21) 
(   22) !!(3) The barotropic velocity in obc_dta_bt() is now computed as the vertical
(   23) !!integral of the baroclinic velocities. Thus the BT vels in the input files 
(   24) !!are not used anymore (they are ignored). Search (e.g.) ubtfoe0
(   25) !!(4) 5 tidal constituents are possible (search: nntide)
(   26) !!(5) tidal and non-tidal forcing are separated (search (e.g.) ubtfoe0, sshfoe0) 
(   27) !!(6) A routine is called that zeros the non-tidal net transport around the
(   28) !! boundaries (search: obc_ctl) 
(   29) !!(7) A ramp function is used for the forcing (search ramp). It is computed in
(   30) !!step.F90 and defined in oce.F90, It is based on kt as
(   31) !!opposed to (kt-nit000) so that by default it does not ramp a restart.  Also
(   32) !!it is complicated to follow when/where certain variables are ramped up -- but no
(   33) !!apologies for this.
(   34) !!(8) The routines are HARDWIRED for 3 open boundaries only. If north is also open, then
(   35) !!all the lp_obc_north areas should be checked and the code added as necessary. 
(   36) !!(9) Uses global masks: emaskg2, smaskg2, wmaskg2 (see obcini, obc_oce) 
(   37) !!(10) NB: use of e2v_e, e2u_s etc seem incorrect, but diff should be v.small so
(   38) !!     modification ... TO DO ... 
(   39) 
(   40) 
(   41) MODULE obcdta 
(   42)   !!==============================================================================
(   43)   !!                            ***  MODULE obcdta  ***
(   44)   !! Open boundary data : read the data for the open boundaries.
(   45)   !!==============================================================================
(   46) #if defined key_obc






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 2

(   47)   !!------------------------------------------------------------------------------
(   48)   !!   'key_obc'         :                                Open Boundary Conditions
(   49)   !!------------------------------------------------------------------------------
(   50)   !!   obc_dta           : read u, v, t, s data along each open boundary
(   51)   !!   obc_dta_psi       : read psi data along each open boundary (rigid lid only)
(   52)   !!------------------------------------------------------------------------------
(   53)   !! * Modules used
(   54)   USE oce             ! ocean dynamics and tracers
(   55)   USE dom_oce         ! ocean space and time domain
(   56)   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   57)   USE phycst          ! physical constants
(   58)   USE obc_oce         ! ocean open boundary conditions
(   59)   USE daymod          ! calendar
(   60)   USE in_out_manager  ! I/O logical units
(   61)   USE lib_mpp         ! distributed memory computing
(   62)   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   63)   USE ioipsl
(   64)   USE obcctl
(   65) !!CN
(   66)   USE lib_ncdf
(   67) !!DBG
(   68) #  if defined key_dynspg_rl
(   69)   USE obccli
(   70) #  endif
(   71) 
(   72)   IMPLICIT NONE
(   73)   PRIVATE
(   74) 
(   75)   !! * Accessibility
(   76)   PUBLIC obc_dta        ! routines called by step.F90
(   77)   PUBLIC obc_dta_bt     ! routines called by dynspg_ts.F90
(   78) 
(   79)   !! * Shared module variables
(   80)   INTEGER ::   &
(   81)      nlecto,   &  ! switch for the first read
(   82)      ntobc1,   &  ! first record used
(   83)      ntobc2,   &  ! second record used
(   84)      ntobc3,   &  ! last record used  for dynspg_ts
(   85)      itobc        ! number of time steps in OBC files
(   86) 
(   87) 
(   88)   REAL(wp), DIMENSION(:), ALLOCATABLE :: ztcobc      ! time_counter variable of BCs
(   89) 
(   90)   !! * Substitutions
(   91) #  include "domzgr_substitute.h90"
(   92) #  include "obc_vectopt_loop_substitute.h90"
(   93) 
(   94)   !!---------------------------------------------------------------------------------
(   95)   !!   OPA 9.0 , LODYC-IPSL  (2003)
(   96)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdta.F90,v 1.9 2006/03/21 08:25:09 opalod Exp $
(   97)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   98)   !!---------------------------------------------------------------------------------
(   99) 
(  100) CONTAINS
(  101) 
(  102)   SUBROUTINE obc_dta (kt)
(  103)      !!--------------------------------------------------------------------
(  104)      !!              ***  SUBROUTINE obc_dta  ***






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 3

(  105)      !!
(  106)      !! ** Purpose :
(  107)      !!   Find the climatological boundary arrays for the specified date,
(  108)      !!   The boundary arrays are netcdf files. Three possible cases:
(  109)      !!   - one time frame only in the file (time dimension = 1).
(  110)      !!     in that case the boundary data does not change in time.
(  111)      !!   - many time frames. In that case,  if we have 12 frames
(  112)      !!     we assume monthly fields.
(  113)      !!     Else, we assume that time_counter is in seconds
(  114)      !!     since the beginning of either the current year or a reference
(  115)      !!     year given in the namelist.
(  116)      !!     (no check is done so far but one would have to check the "unit"
(  117)      !!     attribute of variable time_counter).
(  118)      !!
(  119)      !! History :
(  120)      !!        !  98-05 (J.M. Molines) Original code
(  121)      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  122)      !!   9.0  !  04-06 (F. Durand, A-M. Treguier) Netcdf BC files on input
(  123)      !!--------------------------------------------------------------------
(  124)      !! * Arguments
(  125) 
(  126) 
(  127)      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(  128) 
(  129)      !! * Local declarations
(  130)      INTEGER ::   ji, jj, jk, ii, ij   ! dummy loop indices
(  131)      INTEGER ::   itimo, iman, imois
(  132)      INTEGER ::   i15
(  133)      REAL(wp) ::   zxy
(  134)      !! * Ajouts FD
(  135)      INTEGER ::  isrel              ! number of seconds since 1/1/1992
(  136)      INTEGER, SAVE ::  itobce, itobcw,  & ! number of time steps in OBC files
(  137)                        itobcs, itobcn     !    "       "       "       "
(  138)      INTEGER ::  ikprint        ! frequency for printouts.
(  139)      INTEGER :: fid_e, fid_w, fid_n, fid_s       ! file identifiers
(  140)      LOGICAL :: l_exv
(  141)      INTEGER, DIMENSION(flio_max_dims) ::   f_d  ! dimensions lenght
(  142)      CHARACTER(LEN=25) :: v_name
(  143) !!CN
(  144)      INTEGER :: f_stat
(  145)      INTEGER :: len
(  146)      !!LOGICAL :: use_ioipsl = .FALSE.
(  147) 
(  148) !!DB
(  149)      LOGICAL :: DBG
(  150)      INTEGER :: ncells
(  151)      REAL(wp), DIMENSION(jpjdta,jpk) :: U_EW
(  152) !     REAL(wp), DIMENSION(jpidta,jpk) :: V_NS     !!for future use
(  153) 
(  154) !!DB: For boundary transport adjustments
(  155) !!NB: Could be inside of #ifdef key_ADJ_TRANSPORT
(  156)      REAL(wp) ::   fac_SBI, fac_SB_east, trans_in, trans_diff, trans_check, & 
(  157)                     vfreq, off1,a1,ph1, off2,a2,ph2  !!for time-dependent fac_*
(  158)      INTEGER :: jwest1, jwest2, jeast1, jeast2                !!indices where transport is adjusted 
(  159) 
(  160) 
(  161)      !!--------------------------------------------------------------------
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 4

(  163) !!DB
(  164)      DBG = .false. 
(  165) 
(  166)      IF( lk_dynspg_rl )  THEN
(  167)         CALL obc_dta_psi (kt)     ! update bsf data at open boundaries
(  168)         IF( nobc_dta == 1 .AND. kt == nit000 )   THEN
(  169)            IF(lwp) WRITE(numout,*) ' time-variable psi boundary data not allowed yet'
(  170)            STOP
(  171)         ENDIF
(  172)      ENDIF
(  173) 
(  174)      CALL ipslnlf (new_number=numout)
(  175) 
(  176) !!DB initialize time-dependent valve forcing parameters
(  177) !!Note that these are only used if key_ADJ_TRANSPORT is defined, but I
(  178) !!and lazy and do not #ifdef this
(  179) 
(  180)      vfreq = 2.0*rpi/365.0   ! seasonal cycle in days
(  181) !SBI params
(  182)      off1 = 3.0
(  183)      a1 = 2.0
(  184)      ph1 = 30.0   ! in days
(  185) !SB_east params: CLM0 values
(  186)      off2 = 3.0
(  187)      a2 = 0.5
(  188)      ph2 = -60.0   ! in days
(  189) 
(  190)      ! 1.   First call: check time frames available in files.
(  191)      ! -------------------------------------------------------
(  192) 
(  193)      IF( kt == nit000 )   THEN
(  194) 
(  195)         nlecto =  0
(  196) 
(  197)         IF(lwp) WRITE(numout,*)
(  198)         IF(lwp) WRITE(numout,*)     'obc_dta : find boundary data'
(  199)         IF(lwp) WRITE(numout,*)     '~~~~~~~'
(  200) 
(  201)         IF( nobc_dta == 0 )   THEN
(  202)            IF(lwp) WRITE(numout,*)  '  OBC data taken from initial conditions.'
(  203)            ntobc1 = 1
(  204)            ntobc2 = 1
(  205)         ELSE
(  206)            IF(lwp) WRITE(numout,*)  '  OBC data taken from netcdf files.'
(  207)            IF(lwp) WRITE(numout,*)  '  climatology (T/F):',ln_obc_clim
(  208)            ! check the number of time steps in the files.
(  209)            itobce =0 ; itobcw = 0; itobcn = 0; itobcs = 0
(  210)            v_name = 'time_counter'
(  211)            IF( lp_obc_east )   THEN
(  212)               !!CN: Replacing IOIPSL calls with lib_ncdf
(  213)               !CALL flioopfd ('obceast_TS.nc',fid_e)
(  214)               !CALL flioinqv (fid_e,TRIM(v_name),l_exv,len_dims=f_d)
(  215)               CALL ncdf_get_dim_size('obceast_TS.nc', 'time_counter', len, f_stat)
(  216)               IF( f_stat == 0 )   THEN
(  217)                  f_d(1) = len
(  218)                  itobce = f_d(1)
(  219)               ELSE
(  220)                  if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obceast_TS.nc'






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 5

(  221)               ENDIF
(  222)            ENDIF
(  223)            IF( lp_obc_west )   THEN
(  224)               !!CN: Replacing IOIPSL calls with lib_ncdf
(  225)               !CALL flioopfd ('obcwest_TS.nc',fid_w)
(  226)               !CALL flioinqv (fid_w,TRIM(v_name),l_exv,len_dims=f_d)
(  227)               CALL ncdf_get_dim_size('obcwest_TS.nc', 'time_counter', len, f_stat)
(  228)               IF( f_stat == 0 )   THEN
(  229)                  f_d(1) = len
(  230)                  itobcw = f_d(1)
(  231)               ELSE
(  232)                 if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcwest_TS.nc'
(  233)               ENDIF
(  234)            ENDIF
(  235)            IF( lp_obc_north )   THEN
(  236)               !!CN: Replacing IOIPSL calls with lib_ncdf
(  237)               !CALL flioopfd ('obcnorth_TS.nc',fid_n)
(  238)               !CALL flioinqv (fid_n,TRIM(v_name),l_exv,len_dims=f_d)
(  239)               CALL ncdf_get_dim_size('obcnorth_TS.nc', 'time_counter', len, f_stat)
(  240)               IF( f_stat == 0 )   THEN
(  241)                  f_d(1) = len
(  242)                  itobcn = f_d(1)
(  243)               ELSE
(  244)                  if(lwp)WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcnorth_TS.nc'
(  245)               ENDIF
(  246)            ENDIF
(  247)            IF( lp_obc_south )   THEN
(  248)               !!CN: Replacing IOIPSL calls with lib_ncdf
(  249)               !CALL flioopfd ('obcsouth_TS.nc',fid_s)
(  250)               !CALL flioinqv (fid_s,TRIM(v_name),l_exv,len_dims=f_d)
(  251)               CALL ncdf_get_dim_size('obcsouth_TS.nc', 'time_counter', len, f_stat)
(  252)               IF( f_stat == 0 )   THEN
(  253)                  f_d(1) = len
(  254)                  itobcs = f_d(1)
(  255)               ELSE
(  256)                  if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcsouth_TS.nc'
(  257)               ENDIF
(  258)            ENDIF
(  259) 
(  260)            itobc = MAX(itobce,itobcw,itobcn,itobcs)
(  261)            nstop = 0
(  262)            IF( lp_obc_east  .AND. itobce /= itobc ) nstop = nstop+1
(  263)            IF( lp_obc_west  .AND. itobcw /= itobc ) nstop = nstop+1
(  264)            IF( lp_obc_north .AND. itobcn /= itobc ) nstop = nstop+1
(  265)            IF( lp_obc_south .AND. itobcs /= itobc ) nstop = nstop+1
(  266)            IF( nstop /= 0 )  THEN
(  267)               IF( lwp )   THEN
(  268)                  WRITE(numout,*) ' obcdta : all files must have the same number of time steps'
(  269)                  WRITE(numout,*) ' east, west, north, south: ', itobce, itobcw, itobcn, itobcs
(  270)               ENDIF
(  271)               STOP
(  272)            ENDIF
(  273)            IF( itobc == 1 )   THEN
(  274)               IF( lwp ) WRITE(numout,*) ' obcdta found one time step only in the OBC files'
(  275)            ELSE
(  276)               ALLOCATE (ztcobc(itobc))
(  277)               l_exv = .TRUE.
(  278)               IF( lp_obc_east )   THEN






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 6

(  279)                  !!CN: Replacing IOIPSL calls with lib_ncdf
(  280)                  IF( l_exv )   THEN
(  281)                     !CALL fliogetv (fid_e,TRIM(v_name),ztcobc)
(  282)                     CALL ncdf_read('obceast_TS.nc', 'time_counter', ztcobc, f_stat)
(  283)                     l_exv = .FALSE.
(  284)                  ENDIF
(  285)                  !CALL flioclo (fid_e)
(  286)               ENDIF
(  287)               IF( lp_obc_west )   THEN
(  288)                  !!CN: Replacing IOIPSL calls with lib_ncdf
(  289)                 IF( l_exv )   THEN
(  290)                    !CALL fliogetv (fid_w,TRIM(v_name),ztcobc)
(  291)                    CALL ncdf_read('obcwest_TS.nc', 'time_counter', ztcobc, f_stat)
(  292)                    l_exv = .FALSE.
(  293)                 ENDIF
(  294)                 !CALL flioclo (fid_w)
(  295)               ENDIF
(  296)               IF( lp_obc_north )   THEN
(  297)                  !!CN: Replacing IOIPSL calls with lib_ncdf
(  298)                 IF( l_exv )   THEN
(  299)                    !CALL fliogetv (fid_n,TRIM(v_name),ztcobc)
(  300)                    CALL ncdf_read('obcnorth_TS.nc', 'time_counter', ztcobc, f_stat)
(  301)                    l_exv = .FALSE.
(  302)                 ENDIF
(  303)                 !CALL flioclo (fid_n)
(  304)               ENDIF
(  305)               IF( lp_obc_south )   THEN
(  306)                  !!CN: replacing IOIPSL calls with lib_ncdf
(  307)                 IF( l_exv )   THEN
(  308)                    !CALL fliogetv (fid_s,TRIM(v_name),ztcobc)
(  309)                    CALL ncdf_read('obcsouth_TS.nc', 'time_counter', ztcobc, f_stat)
(  310)                    l_exv = .FALSE.
(  311)                 ENDIF
(  312)                 !CALL flioclo (fid_s)
(  313)               ENDIF
(  314)               IF( lwp ) WRITE(numout,*) ' obcdta found', itobc,' time steps in the OBC files'
(  315)               IF( .NOT. ln_obc_clim .AND. itobc == 12 )   THEN
(  316)                  IF ( lwp ) WRITE(numout,*) '  WARNING: With monthly data we assume climatology'
(  317)                  ln_obc_clim = .true.
(  318)               ENDIF
(  319)            ENDIF
(  320)         ENDIF
(  321) 
(  322)       ! 1.1  Tangential velocities set to zero
(  323)       ! --------------------------------------
(  324)         IF( lp_obc_east  ) vfoe = 0.e0
(  325)         IF( lp_obc_west  ) vfow = 0.e0
(  326)         IF( lp_obc_south ) ufos = 0.e0
(  327)         IF( lp_obc_north ) ufon = 0.e0
(  328) 
(  329)       ! 1.2  Data temperature, salinity, normal velocities set to zero
(  330)       !                        or initial conditions if nobc_dta == 0
(  331)       ! --------------------------------------------------------------
(  332) 
(  333)         IF( lp_obc_east )   THEN
(  334)            ! initialisation to zero
(  335) !byoung           
(  336)            sedta(:,:,:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 7

(  337)            tedta(:,:,:,:) = 0.e0
(  338)            uedta(:,:,:) = 0.e0
(  339)            !                                    ! ================== !
(  340)            IF( nobc_dta == 0 )   THEN           ! initial state used
(  341)               !                                 ! ================== !
(  342)               !  Fills sedta, tedta, uedta (global arrays)
(  343)               !  Remark: this works for njzoom = 1.
(  344)               !          Should the definition of ij include njzoom?
(  345)               DO ji = nie0, nie1
(  346)                  DO jk = 1, jpkm1
(  347)                     DO jj = nje0p1, nje1m1
(  348)                        ij = jj -1 + njmpp
(  349) !byoung                       
(  350)                        sedta(ij,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  351)                        tedta(ij,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  352)                        uedta(ij,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
(  353)                     END DO
(  354)                  END DO
(  355)               END DO
(  356) !byoung
(  357)               DO ji = nie0-4, nie0-1
(  358)                  DO jk = 1, jpkm1
(  359)                     DO jj = nje0p1, nje1m1
(  360)                        ij = jj -1 + njmpp
(  361)                        sedta(ij,jk,nie0-ji+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  362)                        tedta(ij,jk,nie0-ji+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  363)                     END DO
(  364)                  END DO
(  365)               END DO              
(  366)            ENDIF
(  367)         ENDIF
(  368) 
(  369)         IF( lp_obc_west )   THEN
(  370)            ! initialisation to zero
(  371) !byoung           
(  372)            swdta(:,:,:,:) = 0.e0
(  373)            twdta(:,:,:,:) = 0.e0
(  374)            uwdta(:,:,:) = 0.e0
(  375)            !                                    ! ================== !
(  376)            IF( nobc_dta == 0 )   THEN           ! initial state used !
(  377)               !                                 ! ================== !
(  378)               !  Fills swdta, twdta, uwdta (global arrays)
(  379)               !  Remark: this works for njzoom = 1.
(  380)               !          Should the definition of ij include njzoom?
(  381)               DO ji = niw0, niw1
(  382)                  DO jk = 1, jpkm1
(  383)                     DO jj = njw0p1, njw1m1
(  384)                        ij = jj -1 + njmpp
(  385) !byoung                       
(  386)                        swdta(ij,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  387)                        twdta(ij,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  388)                        uwdta(ij,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
(  389)                     END DO
(  390)                  END DO
(  391)               END DO
(  392) !byoung
(  393)               DO ji = niw0+1, niw0+4
(  394)                  DO jk = 1, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 8

(  395)                     DO jj = njw0p1, njw1m1
(  396)                        ij = jj -1 + nimpp
(  397)                        swdta(ij,jk,ji-niw0+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  398)                        twdta(ij,jk,ji-niw0+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  399)                     END DO
(  400)                  END DO
(  401)               END DO
(  402)            ENDIF
(  403)         ENDIF
(  404) 
(  405)         IF( lp_obc_north)   THEN
(  406)            ! initialisation to zero
(  407) 
(  408) !byoung
(  409)            sndta(:,:,:,:) = 0.e0
(  410)            tndta(:,:,:,:) = 0.e0
(  411)            vndta(:,:,:) = 0.e0
(  412)            !                                    ! ================== !
(  413)            IF( nobc_dta == 0 )   THEN           ! initial state used
(  414)               !                                 ! ================== !
(  415)               !  Fills sndta, tndta, vndta (global arrays)
(  416)               !  Remark: this works for njzoom = 1.
(  417)               !          Should the definition of ij include njzoom?
(  418)               DO jj = njn0, njn1
(  419)                  DO jk = 1, jpkm1
(  420)                     DO ji = nin0p1, nin1m1
(  421)                        ii = ji -1 + nimpp
(  422)                        sndta(ii,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  423)                        tndta(ii,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  424)                        vndta(ii,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
(  425)                     END DO
(  426)                  END DO
(  427)               END DO
(  428) !sujie
(  429)               DO jj = njn0-4, njn0-1
(  430)                  DO jk = 1, jpkm1
(  431)                     DO ji = nin0p1, nin1m1
(  432)                        ii = ji -1 + nimpp
(  433)                        sndta(ii,jk,njn0-jj+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  434)                        tndta(ii,jk,njn0-jj+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  435)                     END DO
(  436)                  END DO
(  437)               END DO
(  438)            ENDIF
(  439)         ENDIF
(  440) 
(  441)         IF( lp_obc_south )   THEN
(  442)            ! initialisation to zero
(  443) 
(  444)            ssdta(:,:,:,:) = 0.e0
(  445)            tsdta(:,:,:,:) = 0.e0
(  446)            vsdta(:,:,:) = 0.e0
(  447)            !                                    ! ================== !
(  448)            IF( nobc_dta == 0 )   THEN           ! initial state used
(  449)               !                                 ! ================== !
(  450)               !  Fills ssdta, tsdta, vsdta (global arrays)
(  451)               !  Remark: this works for njzoom = 1.
(  452)               !          Should the definition of ij include njzoom?






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 9

(  453)               DO jj = njs0, njs1
(  454)                  DO jk = 1, jpkm1
(  455)                     DO ji = nis0p1, nis1m1
(  456)                        ii = ji -1 + nimpp
(  457) 
(  458)                        ssdta(ii,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  459)                        tsdta(ii,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  460)                        vsdta(ii,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
(  461)                     END DO
(  462)                  END DO
(  463)               END DO
(  464) !sujie
(  465)               DO jj = njs0+1, njs0+4
(  466)                  DO jk = 1, jpkm1
(  467)                     DO ji = nis0p1, nis1m1
(  468)                        ii = ji -1 + nimpp
(  469)                        ssdta(ii,jk,jj-njs0+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  470)                        tsdta(ii,jk,jj-njs0+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  471)                     END DO
(  472)                  END DO
(  473)               END DO
(  474)            ENDIF
(  475)         ENDIF
(  476) 
(  477) 
(  478)      ENDIF        !       end if kt == nit000
(  479) 
(  480)      ! 2.  Initialize the time we are at.
(  481)      !     Does this every time the routine is called,
(  482)      !     excepted when nobc_dta = 0
(  483)      !---------------------------------------------------------------------
(  484)      IF( nobc_dta == 0 )   THEN
(  485)         itimo = 1
(  486)         zxy   = 0
(  487)      ELSE
(  488)         IF( itobc == 1 )   THEN
(  489)            itimo = 1
(  490)         ELSE IF( itobc == 12 )   THEN      !   BC are monthly
(  491)            ! we assume we have climatology in that case
(  492)            iman  = 12
(  493)            i15   = nday / 16
(  494)            imois = nmonth + i15 - 1
(  495)            IF( imois == 0 )   imois = iman
(  496)            itimo = imois
(  497)         ELSE
(  498)            IF(lwp) WRITE(numout,*) 'data other than constant or monthly',kt
(  499)            iman  = itobc
(  500)            itimo = FLOOR( kt*rdt / (ztcobc(2)-ztcobc(1)) )
(  501)            isrel = kt*rdt
(  502)         ENDIF
(  503)      ENDIF
(  504) 
(  505)      ! 2.1 Read two records in the file if necessary
(  506)      ! ---------------------------------------------
(  507)      IF( ( nobc_dta == 1 ) .AND. ( ( kt == nit000 .AND. nlecto == 0 ) .OR. itimo  /= ntobc1 ) )   THEN
(  508)         nlecto = 1
(  509) 
(  510)         ! Calendar computation






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 10

(  511)         IF( itobc == 1 )   THEN            !  BC are constant in time
(  512)            ntobc1 = 1
(  513)            ntobc2 = 1
(  514)         ELSE IF( itobc == 12 )   THEN      !   BC are monthly
(  515)            ntobc1 = itimo         ! first file record used
(  516)            ntobc2 = ntobc1 + 1    ! last  file record used
(  517)            ntobc1 = MOD( ntobc1, iman )
(  518)            IF( ntobc1 == 0 )   ntobc1 = iman
(  519)            ntobc2 = MOD( ntobc2, iman )
(  520)            IF( ntobc2 == 0 )   ntobc2 = iman
(  521)            IF( lwp )   THEN
(  522)               WRITE(numout,*) ' read monthly obc first record file used ntobc1 ', ntobc1
(  523)               WRITE(numout,*) ' read monthly obc last  record file used ntobc2 ', ntobc2
(  524)            ENDIF
(  525)         ELSE
(  526)            isrel=kt*rdt
(  527)            ntobc1 = itimo         ! first file record used
(  528)            ntobc2 = ntobc1 + 1    ! last  file record used
(  529)            ntobc1 = MOD( ntobc1, iman )
(  530)            IF( ntobc1 == 0 )   ntobc1 = iman
(  531)            ntobc2 = MOD( ntobc2, iman )
(  532)            IF( ntobc2 == 0 )   ntobc2 = iman
(  533)            IF(lwp) WRITE(numout,*) ' read obc first record file used ntobc1 ', ntobc1
(  534)            IF(lwp) WRITE(numout,*) ' read obc last  record file used ntobc2 ', ntobc2
(  535)         ENDIF
(  536)                               ! ======================= !
(  537)                               !  BCs read               !
(  538)         !                     ! ======================= !
(  539) 
(  540) 
(  541)         IF( lp_obc_east )   THEN
(  542)            ! ... Read datafile and set temperature, salinity and normal velocity
(  543)            ! ... initialise the sedta, tedta, uedta arrays
(  544) !!DB/CN: Replace IOIPSL calls with lib_ncdf
(  545)            CALL obc_dta_gv ('y','vosaline',jpjef-jpjed+1,ntobc1,pdta_4D=sedta(jpjed:jpjef,:,1:5,1),fn='obceast_TS.nc')
(  546)            CALL obc_dta_gv ('y','vosaline',jpjef-jpjed+1,ntobc2,pdta_4D=sedta(jpjed:jpjef,:,1:5,2),fn='obceast_TS.nc')
(  547)            CALL obc_dta_gv ('y','votemper',jpjef-jpjed+1,ntobc1,pdta_4D=tedta(jpjed:jpjef,:,1:5,1),fn='obceast_TS.nc')
(  548)            CALL obc_dta_gv ('y','votemper',jpjef-jpjed+1,ntobc2,pdta_4D=tedta(jpjed:jpjef,:,1:5,2),fn='obceast_TS.nc')
(  549)            CALL obc_dta_gv ('y','vozocrtx',jpjef-jpjed+1,ntobc1,pdta_3D=uedta(jpjed:jpjef,:,1),fn='obceast_U.nc')
(  550)            CALL obc_dta_gv ('y','vozocrtx',jpjef-jpjed+1,ntobc2,pdta_3D=uedta(jpjed:jpjef,:,2),fn='obceast_U.nc')
(  551) 
(  552)            !  Usually printout is done only once at kt = nit000,
(  553)            !  unless nprint (namelist) > 1 !!DB -- code deleted, see older version
(  554) 
(  555)         ENDIF
(  556) 
(  557)         IF( lp_obc_west )   THEN
(  558)            ! ... Read datafile and set temperature, salinity and normal velocity
(  559)            ! ... initialise the swdta, twdta, uwdta arrays
(  560) !!DB/CN: Replace IOIPSL calls with lib_ncdf
(  561)            CALL obc_dta_gv ('y','vosaline',jpjwf-jpjwd+1,ntobc1,pdta_4D=swdta(jpjwd:jpjwf,:,1:5,1),fn='obcwest_TS.nc')
(  562)            CALL obc_dta_gv ('y','vosaline',jpjwf-jpjwd+1,ntobc2,pdta_4D=swdta(jpjwd:jpjwf,:,1:5,2),fn='obcwest_TS.nc')
(  563)            CALL obc_dta_gv ('y','votemper',jpjwf-jpjwd+1,ntobc1,pdta_4D=twdta(jpjwd:jpjwf,:,1:5,1),fn='obcwest_TS.nc')
(  564)            CALL obc_dta_gv ('y','votemper',jpjwf-jpjwd+1,ntobc2,pdta_4D=twdta(jpjwd:jpjwf,:,1:5,2),fn='obcwest_TS.nc')
(  565)            CALL obc_dta_gv ('y','vozocrtx',jpjwf-jpjwd+1,ntobc1,pdta_3D=uwdta(jpjwd:jpjwf,:,1),fn='obcwest_U.nc')
(  566)            CALL obc_dta_gv ('y','vozocrtx',jpjwf-jpjwd+1,ntobc2,pdta_3D=uwdta(jpjwd:jpjwf,:,2),fn='obcwest_U.nc')
(  567) !!DB: printout code deleted
(  568)         ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 11

(  569) 
(  570)         IF( lp_obc_north )   THEN
(  571) !!DB/CN: Replace IOIPSL calls with lib_ncdf
(  572)            CALL obc_dta_gv ('x','vosaline',jpinf-jpind+1,ntobc1,pdta_4D=sndta(jpind:jpinf,:,1:5,1),fn='obcnorth_TS.nc')
(  573)            CALL obc_dta_gv ('x','vosaline',jpinf-jpind+1,ntobc2,pdta_4D=sndta(jpind:jpinf,:,1:5,2),fn='obcnorth_TS.nc')
(  574)            CALL obc_dta_gv ('x','votemper',jpinf-jpind+1,ntobc1,pdta_4D=tndta(jpind:jpinf,:,1:5,1),fn='obcnorth_TS.nc')
(  575)            CALL obc_dta_gv ('x','votemper',jpinf-jpind+1,ntobc2,pdta_4D=tndta(jpind:jpinf,:,1:5,2),fn='obcnorth_TS.nc')
(  576)            CALL obc_dta_gv ('x','vomecrty',jpinf-jpind+1,ntobc1,pdta_3D=vndta(jpind:jpinf,:,1),fn='obcnorth_V.nc')
(  577)            CALL obc_dta_gv ('x','vomecrty',jpinf-jpind+1,ntobc2,pdta_3D=vndta(jpind:jpinf,:,2),fn='obcnorth_V.nc')
(  578)         ENDIF
(  579) 
(  580)         IF( lp_obc_south )   THEN
(  581) !!DB/CN: Replace IOIPSL calls with lib_ncdf
(  582)            CALL obc_dta_gv ('x','vosaline',jpisf-jpisd+1,ntobc1,pdta_4D=ssdta(jpisd:jpisf,:,1:5,1),fn='obcsouth_TS.nc')
(  583)            CALL obc_dta_gv ('x','vosaline',jpisf-jpisd+1,ntobc2,pdta_4D=ssdta(jpisd:jpisf,:,1:5,2),fn='obcsouth_TS.nc')
(  584)            CALL obc_dta_gv ('x','votemper',jpisf-jpisd+1,ntobc1,pdta_4D=tsdta(jpisd:jpisf,:,1:5,1),fn='obcsouth_TS.nc')
(  585)            CALL obc_dta_gv ('x','votemper',jpisf-jpisd+1,ntobc2,pdta_4D=tsdta(jpisd:jpisf,:,1:5,2),fn='obcsouth_TS.nc')
(  586)            CALL obc_dta_gv ('x','vomecrty',jpisf-jpisd+1,ntobc1,pdta_3D=vsdta(jpisd:jpisf,:,1),fn='obcsouth_V.nc')
(  587)            CALL obc_dta_gv ('x','vomecrty',jpisf-jpisd+1,ntobc2,pdta_3D=vsdta(jpisd:jpisf,:,2),fn='obcsouth_V.nc')
(  588)         ENDIF
(  589) 
(  590)      ELSE
(  591) 
(  592)         nlecto = 0        !      no reading of OBC barotropic data
(  593) 
(  594)      ENDIF                !      end of the test on the condition to read or not the files
(  595) 
(  596)      ! 3.  Call at every time step :
(  597)      !     Linear interpolation of BCs to current time step
(  598)      ! ----------------------------------------------------
(  599) 
(  600)      IF( itobc == 1 .OR. nobc_dta == 0 )   THEN
(  601)         zxy = 0.
(  602)      ELSE IF( itobc == 12 )   THEN
(  603)         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
(  604)      ELSE
(  605)         zxy = (ztcobc(ntobc1)-FLOAT(isrel))/(ztcobc(ntobc1)-ztcobc(ntobc2))
(  606)      ENDIF
(  607) 
(  608)  
(  609) !!DB: ramp for forcing. NB choice is kt rather than (kt-nit000) so the implicit 
(  610) !!assumption is that ramped forcing is not desirable for restarts.
(  611) !!NB: ramp now computed in step.F90
(  612) !     ramp=tanh(kt*rdt/(2.0*86400.0))
(  613)  
(  614)      IF( lp_obc_east )   THEN
(  615) 
(  616) 
(  617)         DO jk = 1, jpkm1
(  618)            DO jj = nje0p1, nje1m1
(  619)               ij = jj -1 + njmpp
(  620) !!DB
(  621)               sfoe(jj,jk,:) =  ( zxy * sedta(ij,jk,:,2) + &
(  622)                  &           (1.-zxy) * sedta(ij,jk,:,1) ) * temsk5(jj,jk,:) 
(  623)               tfoe(jj,jk,:) =  ( zxy * tedta(ij,jk,:,2) + &
(  624)                  &           (1.-zxy) * tedta(ij,jk,:,1) ) * temsk5(jj,jk,:)
(  625) 
(  626)            END DO






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 12

(  627)         END DO
(  628) !AD   interpolate to global
(  629)         DO jk = 1, jpkm1
(  630)            DO jj = jpjed,jpjef
(  631)               uedta1(jj,jk) =  ( zxy * uedta(jj,jk,2) + &
(  632)                  &           (1.-zxy) * uedta(jj,jk,1) ) !* global mask not * uemsk(jj,jk)
(  633)            END DO
(  634)         END DO
(  635) 
(  636) #ifdef key_ADJ_TRANSPORT
(  637) !!DB: adjust additional transport here, as applicable
(  638) !!DBG
(  639)         if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_east: ADJUSTING TRANSPORT, narea= ',narea
(  640)  
(  641) 
(  642) !!SBI Region
(  643) !!DB07.30: reverse part of SBI inflow as the Nfld coastal current is in wrong direction
(  644)         do jk = 1, jpkm1
(  645)            do jj = 222,224
(  646)               uedta1(jj,jk) = -uedta1(jj,jk)
(  647)            enddo
(  648)            uedta1(225,jk) = 0.0
(  649)         enddo
(  650) 
(  651) !!DB 2007.11.26
(  652)         fac_SBI = off1 + a1*cos(vfreq*(nday_year+ph1))
(  653)         if(lwp .AND. kt-nit000 == 0) write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SBI = ', &
(  654)                     nday_year,fac_SBI
(  655) 
(  656)         do jk = 1, jpkm1
(  657)            do jj = 220,230
(  658)               uedta1(jj,jk) = fac_SBI*uedta1(jj,jk)
(  659)            enddo
(  660)         enddo
(  661) 
(  662) !DB  East Inflow: Shelf break region ---------------------------------
(  663)         jeast1 = 102
(  664)         jeast2 = 108
(  665)         fac_SB_east = off2 + a2*cos(vfreq*(nday_year+ph2))
(  666)         if(lwp .AND. kt-nit000 == 0) &
(  667)              & write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SB_east = ', nday_year, fac_SB_east
(  668)         if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0)write(numout2,'(A51,2x,2(i7,1x),2(f8.3,1x))') &
(  669)           & 'DBG: obcdta -- kt,nday_year, fac_SBI, fac_SB_east ', kt,nday_year,fac_SBI,fac_SB_east
(  670) 
(  671) !!DB: globalize this calc
(  672)         ji = nie0
(  673)         ncells = 0      
(  674)         trans_in = 0.0
(  675)         do jk = 1, jpkm1  
(  676)            do jj = jeast1, jeast2
(  677)               trans_in = trans_in + uedta1(jj,jk)*e2v_e(jj)*e3t(jk)*emaskg2(jj,jk)
(  678)               ncells = ncells + emaskg2(jj,jk)
(  679)            enddo
(  680)         enddo
(  681) 
(  682) !!DB adjust SB_east inflow 
(  683) !!Note that I do not ramp U_EW as we need a correct uedta1 and correct ramping is done elsewhere
(  684)         trans_diff = (-fac_SB_east*1.e6) - trans_in






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 13

(  685)         do jk = 1, jpkm1  
(  686)            do jj = jeast1, jeast2
(  687)               U_EW(jj,jk) = (trans_diff)/ (e2v_e(jj) * e3t(jk)) * emaskg2(jj,jk)
(  688)            enddo
(  689)         enddo
(  690) 
(  691) 	if (ncells .gt. 0) then
(  692) !AD/DB
(  693)            trans_check = 0.0
(  694)            do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
(  695)               do jj = jeast1, jeast2
(  696)                  uedta1(jj,jk) = uedta1(jj,jk) + U_EW(jj,jk)/float(ncells)
(  697)                  trans_check = trans_check + uedta1(jj,jk)*e2v_e(jj)*e3t(jk)*emaskg2(jj,jk)
(  698)               enddo
(  699)            enddo
(  700) 
(  701) !!DBG: check above calcs
(  702)            if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0) then
(  703)               write(numout2,'(A77,2x,2(i7,1x),4(f8.3,1x))') &
(  704)           & 'DBG: obcdta east -- kt,nday_year,fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
(  705)           &            fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
(  706)            endif
(  707)            if(DBG) then
(  708)               write(200+narea,'(A77,2x,3(i5,1x),4(f8.3,1x))') &
(  709)           & 'obcdta east -- kt,nday_year,ncells, fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
(  710)           &            ncells, fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
(  711)            endif
(  712) 
(  713)         end if
(  714)         
(  715) !DB0731:
(  716) !!modify the input at that Nfld bay that looks like a false open boundary
(  717) !!NB: for Barotropic tide cancel these mean transports
(  718) !AD
(  719) 
(  720)         do jk = 1, jpkm1  
(  721)            do jj = 140,143
(  722)               uedta1(jj,jk) = 0.0
(  723)            enddo
(  724)            do jj = 144, 147
(  725) !!DB 2009.06.26 -- reduce this transport to see if/how it affect transport into SW GSL
(  726)               uedta1(jj,jk) = -0.10 * emaskg2(jj,jk)
(  727) !              uedta1(jj,jk) = -0.02 * emaskg2(jj,jk)
(  728)            enddo
(  729)         enddo
(  730) 
(  731) #endif   !ADJ_TRANSPORTS EAST
(  732) 
(  733) !AD   Map global to local
(  734)         DO jk = 1, jpkm1
(  735)            DO jj = nje0p1, nje1m1
(  736)               ij = jj -1 + njmpp
(  737)               ufoe(jj,jk) =  uedta1(ij,jk)* uemsk(jj,jk)
(  738)            END DO
(  739)         END DO
(  740) 
(  741) !AD >>>: recompute elevation
(  742)         eta_e(jpjed)=0. !zero eta at SW corner






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 14

(  743) 	DO ij=jpjed+1,jpjef
(  744)       	   eta_e(ij) = eta_e(ij-1) - e2v_e(ij-1)*ff_e(ij)/grav * &
(  745)            &    uedta1(ij-1,1)
(  746) 	END DO
(  747)  !AD <<<
(  748) !!DBG
(  749)         if(DBG .AND. kt == nit000) then
(  750)            write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
(  751)                 'DBG: kt,narea,max,min ufoe,uedta1,eta_e ', &
(  752)                 kt,narea,MAXVAL(ufoe(:,:)),MINVAL(ufoe(:,:)), &
(  753)                 MAXVAL(uedta1(:,:)),MINVAL(uedta1(:,:)), &
(  754)                 MAXVAL(eta_e(:)),MINVAL(eta_e(:))
(  755)            DO jk = 1, jpkm1
(  756)               DO jj = nje0p1, nje1m1
(  757)                  ij = jj -1 + njmpp
(  758)                  write(3000+narea,'(3(i3,1x),3(e12.6,1x))') &
(  759)                       jk,jj,ij,ufoe(jj,jk),uedta1(ij,jk),uemsk(jj,jk)
(  760)               END DO
(  761)            END DO
(  762)            do ij = jpjed,jpjef
(  763)               write(8000+narea,'(i3,2x,2(e12.6,2x))')ij,eta_e(ij),uedta1(ij,1)
(  764)            enddo
(  765) 
(  766)         endif
(  767) 
(  768) 
(  769)         ufoe(:,:) = ramp*ufoe(:,:)
(  770) 
(  771)      ENDIF
(  772) 
(  773) 
(  774)      IF( lp_obc_west )   THEN
(  775) 
(  776)         DO jk = 1, jpkm1
(  777)            DO jj = njw0p1, njw1m1
(  778)               ij = jj -1 + njmpp
(  779) !!DB
(  780)               sfow(jj,jk,:) =  ( zxy * swdta(ij,jk,:,2) + &
(  781)                  &           (1.-zxy) * swdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
(  782)               tfow(jj,jk,:) =  ( zxy * twdta(ij,jk,:,2) + &
(  783)                  &           (1.-zxy) * twdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
(  784)            END DO
(  785)         END DO
(  786) !AD   interpolate to global
(  787)         DO jk = 1, jpkm1
(  788)            DO jj = jpjwd,jpjwf
(  789)               uwdta1(jj,jk) =  ( zxy * uwdta(jj,jk,2) + &
(  790)                  &           (1.-zxy) * uwdta(jj,jk,1) ) !* global mask not * uemsk(jj,jk)
(  791)            END DO
(  792)         END DO
(  793) 
(  794) 
(  795) #ifdef key_ADJ_TRANSPORT
(  796) !!DB: add additional transport here, as applicable
(  797) !!DBG
(  798)         if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_west ADJUSTING TRANSPORT, narea= ',narea
(  799) 
(  800) !DB  west Inflow---------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 15

(  801) !!DBG
(  802)         if(lwp .AND. kt-nit000 == 0) &
(  803)              &  write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SB_west = fac_SB_east = ',nday_year, fac_SB_east
(  804) 
(  805)         jwest1 = 20
(  806)         jwest2 = 35
(  807)         ji = niw0
(  808)         ncells = 0     
(  809)         trans_in = 0.0
(  810)         do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
(  811)            do jj = jwest1,jwest2
(  812)               trans_in = trans_in + uwdta1(jj,jk)*e2v_w(jj)*e3t(jk)*wmaskg2(jj,jk)
(  813)               ncells = ncells + wmaskg2(jj,jk)
(  814)            enddo
(  815)         enddo
(  816) 
(  817) !!DB adjust SB_east inflow (NB:  -fac_SB_east*1.e6 == inflow in Sv)
(  818) !!Note that I do not ramp U_EW as we need a correct uedta1 and correct ramping is done elsewhere
(  819)         trans_diff = (-fac_SB_east*1.e6) - trans_in
(  820)         do jk = 1, jpkm1  
(  821)            do jj = jwest1, jwest2
(  822)               U_EW(jj,jk) = (trans_diff)/ (e2v_w(jj) * e3t(jk)) * wmaskg2(jj,jk)
(  823)            enddo
(  824)         enddo
(  825) 
(  826) 	
(  827) 	if (ncells .gt. 0) then
(  828) !AD/DB
(  829)            trans_check = 0.0
(  830)            do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
(  831)               do jj = jwest1, jwest2
(  832)                  uwdta1(jj,jk) = uwdta1(jj,jk) + U_EW(jj,jk)/float(ncells)
(  833)                  trans_check = trans_check + uwdta1(jj,jk)*e2v_w(jj)*e3t(jk)*wmaskg2(jj,jk)
(  834)               enddo
(  835)            enddo
(  836) !!DBG: check above calcs
(  837)            if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0) then
(  838)               write(numout2,'(A77,2x,2(i7,1x),4(f8.3,1x))') &
(  839)           & 'DBG: obcdta west -- kt,nday_year, fac_SB_east,trans_in,trans_diff,trans_check ', kt,nday_year, &
(  840)           &            fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
(  841) 
(  842)            endif
(  843)            if(DBG) then
(  844)               write(200+narea,'(A77,2x,3(i5,1x),4(f8.3,1x))') &
(  845)           & 'obcdta west -- kt,nday_year,ncells, fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
(  846)           &            ncells, fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
(  847)            endif
(  848) 
(  849)         end if
(  850) #endif     !ADJ_TRANSPORT WEST
(  851) 
(  852) !AD   Map global to local
(  853)         DO jk = 1, jpkm1
(  854)            DO jj = njw0p1, njw1m1
(  855)               ij = jj -1 + njmpp
(  856)               ufow(jj,jk) =  uwdta1(ij,jk)* uwmsk(jj,jk)
(  857)            END DO
(  858)         END DO






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 16

(  859) 
(  860) !AD >>>: recompute elevation
(  861)         eta_w(jpjwd)=0. !zero eta at SW corner
(  862) 	DO ij=jpjwd+1,jpjwf
(  863)       	   eta_w(ij) = eta_w(ij-1) - e2v_w(ij-1)*ff_w(ij)/grav * &
(  864)            &    uwdta1(ij-1,1)
(  865) 	END DO
(  866)  !AD <<<
(  867) 
(  868) !!DBG
(  869)         if(DBG .AND. kt == nit000) then
(  870)            write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
(  871)                 'DBG: kt,narea,max,min ufow,uwdta1,eta_w ', &
(  872)                 kt,narea,MAXVAL(ufow(:,:)),MINVAL(ufow(:,:)), &
(  873)                 MAXVAL(uwdta1(:,:)),MINVAL(uwdta1(:,:)), &
(  874)                 MAXVAL(eta_w(:)),MINVAL(eta_w(:))
(  875)            DO jk = 1, jpkm1
(  876)               DO jj = njw0p1, njw1m1
(  877)                  ij = jj -1 + njmpp
(  878)                  write(4000+narea,'(3(i3,1x),3(e12.6,1x))') &
(  879)                       jk,jj,ij,ufow(jj,jk),uwdta1(ij,jk),uwmsk(jj,jk)
(  880)               END DO
(  881)            END DO
(  882)            do ij = jpjwd,jpjwf
(  883)               write(8100+narea,'(i3,2x,2(e12.6,2x))')ij,eta_w(ij),uwdta1(ij,1)
(  884)            enddo
(  885) 
(  886)         endif
(  887) 
(  888)         ufow(:,:) = ramp * ufow(:,:)
(  889) 
(  890)      ENDIF
(  891) !------------------------------------
(  892) 
(  893) !!DB: north is closed so do nothing. Note that this could lead to future probs
(  894)      IF( lp_obc_north )   THEN
(  895)         !  fills sfon, tfon, vfon (local to each processor)
(  896)         DO jk = 1, jpkm1
(  897)            DO ji = nin0p1, nin1m1
(  898)               ii = ji -1 + nimpp
(  899) !sujie
(  900)               sfon(ji,jk,:) =  ( zxy * sndta(ii,jk,:,2) + &
(  901)                  &           (1.-zxy) * sndta(ii,jk,:,1) )* tnmsk5(ji,jk,:)
(  902)               tfon(ji,jk,:) =  ( zxy * tndta(ii,jk,:,2) + &
(  903)                  &           (1.-zxy) * tndta(ii,jk,:,1) )* tnmsk5(ji,jk,:)
(  904)               vfon(ji,jk) =  ( zxy * vndta(ii,jk,2) + &
(  905)                  &           (1.-zxy) * vndta(ii,jk,1) ) * vnmsk(ji,jk)
(  906)            END DO
(  907)         END DO
(  908)      ENDIF
(  909) 
(  910)      IF( lp_obc_south )   THEN
(  911) 
(  912)         DO jk = 1, jpkm1
(  913)           DO ji = nis0p1, nis1m1
(  914)              ii = ji -1 + nimpp
(  915) 
(  916)              sfos(ji,jk,:) = ( zxy * ssdta(ii,jk,:,2) + &






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 17

(  917)                 &          (1.-zxy) * ssdta(ii,jk,:,1) )* tsmsk5(ji,jk,:)
(  918)              tfos(ji,jk,:) = ( zxy * tsdta(ii,jk,:,2) + &
(  919)                 &          (1.-zxy) * tsdta(ii,jk,:,1) )* tsmsk5(ji,jk,:)
(  920)    
(  921)             END DO
(  922)         END DO
(  923) 
(  924) !AD >>> interpolate to global
(  925)         DO jk = 1, jpkm1
(  926)           DO ii = jpisd,jpisf
(  927)              vsdta1(ii,jk) = ( zxy * vsdta(ii,jk,2) + &
(  928)                 &          (1.-zxy) * vsdta(ii,jk,1) ) !* global mask not * vsmsk(ji,jk)       
(  929)   
(  930)             END DO
(  931)         END DO
(  932) !AD <<<
(  933) 
(  934) #ifdef key_ADJ_TRANSPORT
(  935) !!DB: Adjust transport here if desired
(  936) !        if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_south ADJUSTING TRANSPORT'
(  937) #endif    !ADJ_TRANSPORT SOUTH
(  938) 
(  939) !   Map global to local
(  940)         DO jk = 1, jpkm1
(  941)           DO ji = nis0p1, nis1m1
(  942)              ii = ji -1 + nimpp
(  943)              vfos(ji,jk) = vsdta1(ii,jk) * vsmsk(ji,jk)  
(  944)             END DO
(  945)         END DO
(  946) 
(  947) !  recompute elevation
(  948)         eta_s(jpisd)=0. !zero eta at SW corner
(  949) 	DO ij=jpisd+1,jpisf
(  950)       	   eta_s(ij) = eta_s(ij-1) + e2u_s(ij-1)*ff_s(ij)/grav * &
(  951)            &    vsdta1(ij-1,1)
(  952) 	END DO
(  953) 
(  954) !!DBG
(  955)         if(DBG .AND. kt == nit000) then
(  956)            write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
(  957)                 'DBG: kt,narea,max,min vfos,vsdta1,eta_s ', &
(  958)                 kt,narea,MAXVAL(vfos(:,:)),MINVAL(vfos(:,:)), &
(  959)                 MAXVAL(vsdta1(:,:)),MINVAL(vsdta1(:,:)), &
(  960)                 MAXVAL(eta_s(:)),MINVAL(eta_s(:))
(  961)            DO jk = 1, jpkm1
(  962)               DO ji = nis0p1, nis1m1
(  963)                  ii = ji -1 + nimpp
(  964)                  write(5000+narea,'(3(i3,1x),3(e12.6,1x))') &
(  965)                       jk,ji,ii,vfos(ji,jk),vsdta1(ii,jk), vsmsk(ji,jk)  
(  966)               END DO
(  967)            END DO
(  968)            do ij = jpisd,jpisf
(  969)               write(8200+narea,'(i3,2x,2(e12.6,2x))')ij,eta_s(ij),vsdta1(ij,1)
(  970)            enddo
(  971) 
(  972) 
(  973) 
(  974)         endif






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 18

(  975) 
(  976) 
(  977)         vfos(:,:) = ramp * vfos(:,:)
(  978) 
(  979)      ENDIF
(  980) 
(  981) !!DB: Add offset to East boundary SSH -- HARDWIRED
(  982)      do ij = jpjed+1, jpjef
(  983)         eta_e(ij) = eta_e(ij) + eta_s(jpisf)
(  984)      enddo
(  985) 
(  986) 
(  987) !!DB: call ZW routine that zeros non-tidal transports along obc's
(  988)      call obc_ctl(kt)
(  989) 
(  990) 
(  991)   END SUBROUTINE obc_dta
(  992) 
(  993) 
(  994) # if defined key_dynspg_rl
(  995)   !!-----------------------------------------------------------------------------
(  996)   !!   Rigid-lid
(  997)   !!-----------------------------------------------------------------------------
(  998) 
(  999)   SUBROUTINE obc_dta_psi ( kt )
( 1000)      !!-----------------------------------------------------------------------------
( 1001)      !!                       ***  SUBROUTINE obc_dta_psi  ***
( 1002)      !!
( 1003)      !! ** Purpose :
( 1004)      !!      Update the climatological streamfunction OBC at each time step.
( 1005)      !!      Depends on the user's configuration.  Here data are read only once
( 1006)      !!      at the beginning of the run.
( 1007)      !!
( 1008)      !! ** Method :
( 1009)      !!      1. initialization
( 1010)      !!         kbsfstart: number of time steps over which increase bsf
( 1011)      !!         during initialization. This is provided for a smooth start
( 1012)      !!         in cases where the transport is large (like on the Antarctic
( 1013)      !!         continent). also note that when kbfstart=1, the transport
( 1014)      !!         increases a lot in one time step and the precision usually
( 1015)      !!         required for the solver may not be enough.
( 1016)      !!      2. set the time evolution of the climatological barotropic streamfunction
( 1017)      !!         along the isolated coastlines ( gcbic(jnic) ).
( 1018)      !!      3. set the climatological barotropic streamfunction at the boundary.
( 1019)      !!
( 1020)      !!      The last two steps are done only at first step (nit000) or if kt <= kbfstart
( 1021)      !!
( 1022)      !! History :
( 1023)      !!        ! 97-08 (G. Madec, J.M. Molines)
( 1024)      !!   8.5  ! 02-10 (C. Talandier, A-M. Treguier) Free surface, F90
( 1025)      !!   9.0  ! 05-11  (V. Garnier) Surface pressure gradient organization
( 1026)      !!----------------------------------------------------------------------------
( 1027)      !! * Arguments
( 1028)      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
( 1029) 
( 1030)      !! * Local declarations
( 1031)      INTEGER ::   ji, jj, jnic, jip         ! dummy loop indices
( 1032)      INTEGER ::   inum = 11                 ! temporary logical unit






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 19

( 1033)      INTEGER ::   ip, ii, ij, iii, ijj
( 1034)      INTEGER ::   kbsfstart
( 1035)      REAL(wp) ::   zsver1, zsver2, zsver3, z2dtr, zcoef
( 1036)      !!----------------------------------------------------------------------------
( 1037) 
( 1038)      ! 1. initialisation
( 1039)      ! -----------------
( 1040) 
( 1041)      kbsfstart =  1
( 1042)      zsver1 =  bsfic0(1)
( 1043)      zsver2 =  zsver1
( 1044)      IF( kt <= kbsfstart )   THEN
( 1045)         zcoef = float(kt)/float(kbsfstart)
( 1046)      ELSE
( 1047)         zcoef = 1.
( 1048)      END IF
( 1049)      bsfic(1) = zsver1*zcoef
( 1050)      IF( lwp .AND. ( kt <= kbsfstart ) )   THEN
( 1051)         IF(lwp) WRITE(numout,*)'            '
( 1052)         IF(lwp) WRITE(numout,*)'obcdta: spinup phase in obc_dta_psi routine'
( 1053)         IF(lwp) WRITE(numout,*)'~~~~~~  it=',kt,'  OBC: spinup coef: ', &
( 1054)                                           zcoef, ' and transport: ',bsfic(1)
( 1055)      END IF
( 1056) 
( 1057)      zsver2 =  bsfic(1)-bsfic(2)
( 1058)      zsver3 =  bsfic(2)
( 1059) 
( 1060)      ! 2. Right hand side of the barotropic elliptic equation (isolated coastlines)
( 1061)      ! ----------------------------------------------------------------------------
( 1062) 
( 1063)      IF( ( neuler == 0 ) .AND. ( kt == nit000 ) )   THEN
( 1064)         z2dtr = 1./rdt
( 1065)      ELSE
( 1066)         z2dtr = 1./2./rdt
( 1067)      END IF
( 1068)      ! ... bsfb(ii,ij) should be constant but due to the Asselin filter it
( 1069)      ! ... converges asymptotically towards bsfic(jnic)
( 1070)      ! ... However, bsfb(ii,ij) is constant along the same coastlines
( 1071)      ! ... ---> can be improved using an extra array for storing bsficb (before)
( 1072)      IF( nbobc > 1 )   THEN
( 1073)         DO jnic = 1,nbobc - 1
( 1074)            gcbic(jnic) = 0.e0
( 1075)            ip=mnic(0,jnic)
( 1076)            DO jip = 1,ip
( 1077)               ii = miic(jip,0,jnic)
( 1078)               ij = mjic(jip,0,jnic)
( 1079)               IF( ii >= nldi+ nimpp - 1 .AND. ii <= nlei+ nimpp - 1 .AND. &
( 1080)                   ij >= nldj+ njmpp - 1 .AND. ij <= nlej+ njmpp - 1 )   THEN
( 1081)                  iii=ii-nimpp+1
( 1082)                  ijj=ij-njmpp+1
( 1083)                  gcbic(jnic) = ( bsfic(jnic) - bsfb(iii,ijj) ) * z2dtr
( 1084)               END IF
( 1085)            END DO
( 1086)         END DO
( 1087)      END IF
( 1088) 
( 1089)      IF( lk_mpp )   CALL mpp_isl( gcbic, 3 )
( 1090) 






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 20

( 1091)      ! 3. Update the climatological barotropic function at the boundary
( 1092)      ! ----------------------------------------------------------------
( 1093) 
( 1094)      IF( lpeastobc )   THEN
( 1095) 
( 1096)         IF( kt == nit000 .OR. kt <= kbsfstart )   THEN
( 1097)            OPEN(inum,file='obceastbsf.dta')
( 1098)            READ(inum,*)
( 1099)            READ(inum,*)
( 1100)            READ(inum,*)
( 1101)            READ(inum,*)
( 1102)            READ(inum,*)
( 1103)            READ(inum,*) (bfoe(jj),jj=jpjed, jpjef)
( 1104)            CLOSE(inum)
( 1105)         END IF
( 1106)         DO jj=jpjed, jpjefm1
( 1107)            bfoe(jj)=bfoe(jj)*zcoef
( 1108)         END DO
( 1109) 
( 1110)      END IF
( 1111) 
( 1112)      IF( lpwestobc)   THEN
( 1113) 
( 1114)         IF( kt == nit000 .OR. kt <= kbsfstart ) then
( 1115)            OPEN(inum,file='obcwestbsf.dta')
( 1116)            READ(inum,*)
( 1117)            READ(inum,*)
( 1118)            READ(inum,*)
( 1119)            READ(inum,*)
( 1120)            READ(inum,*)
( 1121)            READ(inum,*) (bfow(jj),jj=jpjwd, jpjwf)
( 1122)            CLOSE(inum)
( 1123)         END IF
( 1124)         DO jj=jpjwd, jpjwfm1
( 1125)            bfow(jj)=bfow(jj)*zcoef
( 1126)         END DO
( 1127) 
( 1128)      END IF
( 1129) 
( 1130)      IF( lpsouthobc)   THEN
( 1131) 
( 1132)         IF( kt == nit000 .OR. kt <= kbsfstart )   THEN
( 1133)            OPEN(inum,file='obcsouthbsf.dta')
( 1134)            READ(inum,*)
( 1135)            READ(inum,*)
( 1136)            READ(inum,*)
( 1137)            READ(inum,*)
( 1138)            READ(inum,*)
( 1139)            READ(inum,*) (bfos(jj),jj=jpisd, jpisf)
( 1140)            CLOSE(inum)
( 1141)         END IF
( 1142)         DO ji=jpisd, jpisfm1
( 1143)            bfos(ji)=bfos(ji)*zcoef
( 1144)         END DO
( 1145) 
( 1146)      END IF
( 1147) 
( 1148)      IF( lpnorthobc)   THEN






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 21

( 1149)         IF( kt == nit000 .OR. kt <= kbsfstart )   THEN
( 1150)            OPEN(inum,file='obcnorthbsf.dta')
( 1151)            READ(inum,*)
( 1152)            READ(inum,*)
( 1153)            READ(inum,*)
( 1154)            READ(inum,*)
( 1155)            READ(inum,*)
( 1156)            READ(inum,*) (bfon(jj),jj=jpind, jpinf)
( 1157)            CLOSE(inum)
( 1158)         END IF
( 1159)         DO ji=jpind, jpinfm1
( 1160)            bfon(ji)=bfon(ji)*zcoef
( 1161)         END DO
( 1162) 
( 1163)      END IF
( 1164) 
( 1165)   END SUBROUTINE obc_dta_psi
( 1166) #else
( 1167)   !!-----------------------------------------------------------------------------
( 1168)   !!   Default option
( 1169)   !!-----------------------------------------------------------------------------
( 1170)   SUBROUTINE obc_dta_psi ( kt )       ! Empty routine
( 1171)      !! * Arguments
( 1172)      INTEGER,INTENT(in) :: kt
( 1173)      if(lwp) WRITE(numout,*) 'obc_dta_psi: You should not have seen this print! error?', kt
( 1174)   END SUBROUTINE obc_dta_psi
( 1175) # endif
( 1176) 
( 1177) 
( 1178) #if defined key_dynspg_ts || defined key_dynspg_exp
( 1179)   SUBROUTINE obc_dta_bt( kt, kbt )
( 1180)      !!---------------------------------------------------------------------------
( 1181)      !!                      ***  SUBROUTINE obc_dta  ***
( 1182)      !!
( 1183)      !! ** Purpose :   time interpolation of barotropic data for time-splitting scheme
( 1184)      !!                Data at the boundary must be in m2/s
( 1185)      !!
( 1186)      !! History :
( 1187)      !!   9.0  !  05-11 (V. garnier) Original code
( 1188)      !!---------------------------------------------------------------------------
( 1189)      !! * Arguments
( 1190)      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
( 1191)      INTEGER, INTENT( in ) ::   kbt         ! barotropic ocean time-step index
( 1192) 
( 1193)      !! * Local declarations
( 1194)      INTEGER ::   ji, jj, jk, ii, ij   ! dummy loop indices
( 1195)      INTEGER ::   fid_e, fid_w, fid_n, fid_s, fid  ! file identifiers
( 1196)      INTEGER ::   itimo, iman, imois, i15
( 1197)      INTEGER ::   ntobcm, ntobcp, itimom, itimop
( 1198)      REAL(wp) ::  zxy
( 1199)      INTEGER ::   isrel, ikt           ! number of seconds since 1/1/1992
( 1200)      INTEGER ::   ikprint              ! frequency for printouts.
( 1201) !ylu/DB
( 1202)      INTEGER ::   ntide, nntide   
( 1203)      REAL(wp), DIMENSION(5) :: tideperiod    ! in hours
( 1204) !!DB
( 1205) !!Local arrays for BT transport, dimensioned as in obc_oce
( 1206)      REAL(wp), DIMENSION(1:jpisf) ::    &   !:






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 22

( 1207)      vbtfos0, sshfos0       !: south boundary barotropic transport without the tide
( 1208)      REAL(wp), DIMENSION(1:jpjwf) ::   &  !:
( 1209)      ubtfow0, sshfow0       !: west boundary barotropic transport without the tide
( 1210)      REAL(wp), DIMENSION(1:jpjef) ::   &  !:
( 1211)      ubtfoe0, sshfoe0       !: east boundary barotropic transport without the tide
( 1212)      INTEGER :: status
( 1213) 
( 1214) !!DB: 1 tide= M2 = good for climatology runs; 5 tides ~ operational runs
( 1215) !     nntide =5
( 1216)      nntide =1
( 1217)      tideperiod(1) = 12.4206  !M2
( 1218)      tideperiod(2) = 12.0000  !S2
( 1219)      tideperiod(3) = 25.8193  !O1
( 1220)      tideperiod(4) = 23.9345  !K1
( 1221)      tideperiod(5) = 12.6583  !N2
( 1222) 
( 1223) 
( 1224)      !!---------------------------------------------------------------------------
( 1225) 
( 1226)      ! 1.   First call: check time frames available in files.
( 1227)      ! -------------------------------------------------------
( 1228) 
( 1229)      IF( kt == nit000 ) THEN
( 1230) 
( 1231)         ! 1.1  Barotropic tangential velocities set to zero
( 1232)         ! -------------------------------------------------
( 1233)         IF( lp_obc_east  ) vbtfoe(:) = 0.e0
( 1234)         IF( lp_obc_west  ) vbtfow(:) = 0.e0
( 1235)         IF( lp_obc_south ) ubtfos(:) = 0.e0
( 1236)         IF( lp_obc_north ) ubtfon(:) = 0.e0
( 1237) 
( 1238)         ! 1.2  Sea surface height and normal barotropic velocities set to zero
( 1239)         !                               or initial conditions if nobc_dta == 0
( 1240)         ! --------------------------------------------------------------------
( 1241) 
( 1242)          IF( lp_obc_east ) THEN
( 1243)             ! initialisation to zero
( 1244)             sshedta(:,:) = 0.e0
( 1245)             ubtedta(:,:) = 0.e0
( 1246)             !                                        ! ================== !
( 1247)             IF( nobc_dta == 0 )   THEN               ! initial state used !
( 1248)                !                                     ! ================== !
( 1249)                !  Fills sedta, tedta, uedta (global arrays)
( 1250)                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
( 1251)                DO ji = nie0, nie1
( 1252)                   DO jj = nje0p1, nje1m1
( 1253)                      ij = jj -1 + njmpp
( 1254)                      sshedta(ij,1) = sshn(ji+1,jj) * tmask(ji+1,jj,1)
( 1255)                   END DO
( 1256)                END DO
( 1257)             ENDIF
( 1258) !ylu
( 1259)            tidesshemag(:,:)=0.0
( 1260)            tidesshepha(:,:)=0.0
( 1261)            tidevbtemag(:,:)=0.0
( 1262)            tidevbtepha(:,:)=0.0
( 1263) !!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
( 1264)            call ncdf_read_global('obceast_tide.nc','tidesshmag',tidesshemag(jpjed:jpjef,1:nntide),-1,status)






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 23

( 1265)            call ncdf_read_global('obceast_tide.nc','tidesshpha',tidesshepha(jpjed:jpjef,1:nntide),-1,status)
( 1266)            call ncdf_read_global('obceast_tide.nc','tidevbtmag',tidevbtemag(jpjed:jpjef,1:nntide),-1,status)
( 1267)            call ncdf_read_global('obceast_tide.nc','tidevbtpha',tidevbtepha(jpjed:jpjef,1:nntide),-1,status)
( 1268) 
( 1269)          ENDIF
( 1270) 
( 1271)          IF( lp_obc_west) THEN
( 1272)             ! initialisation to zero
( 1273)             sshwdta(:,:) = 0.e0
( 1274)             ubtwdta(:,:) = 0.e0
( 1275)             !                                        ! ================== !
( 1276)             IF( nobc_dta == 0 )   THEN               ! initial state used !
( 1277)                !                                     ! ================== !
( 1278)                !  Fills swdta, twdta, uwdta (global arrays)
( 1279)                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
( 1280)                DO ji = niw0, niw1
( 1281)                   DO jj = njw0p1, njw1m1
( 1282)                      ij = jj -1 + njmpp
( 1283)                      sshwdta(ij,1) = sshn(ji,jj) * tmask(ji,jj,1)
( 1284)                   END DO
( 1285)                END DO
( 1286)             ENDIF
( 1287) !ylu
( 1288)            tidesshwmag(:,:)=0.0
( 1289)            tidesshwpha(:,:)=0.0
( 1290)            tidevbtwmag(:,:)=0.0
( 1291)            tidevbtwpha(:,:)=0.0
( 1292) !!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
( 1293)            call ncdf_read_global('obcwest_tide.nc','tidesshmag',tidesshwmag(jpjwd:jpjwf,1:nntide),-1,status)
( 1294)            call ncdf_read_global('obcwest_tide.nc','tidesshpha',tidesshwpha(jpjwd:jpjwf,1:nntide),-1,status)
( 1295)            call ncdf_read_global('obcwest_tide.nc','tidevbtmag',tidevbtwmag(jpjwd:jpjwf,1:nntide),-1,status)
( 1296)            call ncdf_read_global('obcwest_tide.nc','tidevbtpha',tidevbtwpha(jpjwd:jpjwf,1:nntide),-1,status)
( 1297) 
( 1298) 
( 1299)            !CALL flioclo (fid_w)
( 1300)          ENDIF
( 1301) 
( 1302)          IF( lp_obc_north) THEN
( 1303)             ! initialisation to zero
( 1304)             sshndta(:,:) = 0.e0
( 1305)             vbtndta(:,:) = 0.e0
( 1306)             !                                        ! ================== !
( 1307)             IF( nobc_dta == 0 )   THEN               ! initial state used !
( 1308)                !                                     ! ================== !
( 1309)                !  Fills sndta, tndta, vndta (global arrays)
( 1310)                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
( 1311)                DO jj = njn0, njn1
( 1312)                   DO ji = nin0p1, nin1m1
( 1313)                      DO jk = 1, jpkm1
( 1314)                         ii = ji -1 + nimpp
( 1315)                         vbtndta(ii,1) = vbtndta(ii,1) + vndta(ii,jk,1)*fse3v(ji,jj,jk)
( 1316)                      END DO
( 1317)                      sshndta(ii,1) = sshn(ii,jj+1) * tmask(ji,jj+1,1)
( 1318)                   END DO
( 1319)                END DO
( 1320)             ENDIF
( 1321) !ylu
( 1322)            tidesshnmag(:,:)=0.0






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 24

( 1323)            tidesshnpha(:,:)=0.0
( 1324)            tidevbtnmag(:,:)=0.0
( 1325)            tidevbtnpha(:,:)=0.0
( 1326) !!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
( 1327)            call ncdf_read_global('obcnorth_tide.nc','tidesshmag',tidesshnmag(jpind:jpinf,1:nntide),-1,status)
( 1328)            call ncdf_read_global('obcnorth_tide.nc','tidesshpha',tidesshnpha(jpind:jpinf,1:nntide),-1,status)
( 1329)            call ncdf_read_global('obcnorth_tide.nc','tidevbtmag',tidevbtnmag(jpind:jpinf,1:nntide),-1,status)
( 1330)            call ncdf_read_global('obcnorth_tide.nc','tidevbtpha',tidevbtnpha(jpind:jpinf,1:nntide),-1,status)
( 1331) 
( 1332)          ENDIF
( 1333) 
( 1334)          IF( lp_obc_south) THEN
( 1335)             ! initialisation to zero
( 1336)             sshsdta(:,:) = 0.e0
( 1337)             vbtsdta(:,:) = 0.e0
( 1338)             !                                        ! ================== !
( 1339)             IF( nobc_dta == 0 )   THEN               ! initial state used !
( 1340)                !                                     ! ================== !
( 1341)                !  Fills ssdta, tsdta, vsdta (global arrays)
( 1342)                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
( 1343)                DO jj = njs0, njs1
( 1344)                   DO ji = nis0p1, nis1m1
( 1345)                      DO jk = 1, jpkm1
( 1346)                         ii = ji -1 + nimpp
( 1347)                         vbtsdta(ii,1) = vbtsdta(ii,1) + vsdta(ii,jk,1)*fse3v(ji,jj,jk)
( 1348)                      END DO
( 1349)                      sshsdta(ii,1) = sshn(ji,jj) * tmask(ii,jj,1)
( 1350)                   END DO
( 1351)                END DO
( 1352)             ENDIF
( 1353) !ylu
( 1354)            tidesshsmag(:,:)=0.0
( 1355)            tidesshspha(:,:)=0.0
( 1356)            tidevbtsmag(:,:)=0.0
( 1357)            tidevbtspha(:,:)=0.0
( 1358) !!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
( 1359)            call ncdf_read_global('obcsouth_tide.nc','tidesshmag',tidesshsmag(jpisd:jpisf,1:nntide),-1,status)
( 1360)            call ncdf_read_global('obcsouth_tide.nc','tidesshpha',tidesshspha(jpisd:jpisf,1:nntide),-1,status)
( 1361)            call ncdf_read_global('obcsouth_tide.nc','tidevbtmag',tidevbtsmag(jpisd:jpisf,1:nntide),-1,status)
( 1362)            call ncdf_read_global('obcsouth_tide.nc','tidevbtpha',tidevbtspha(jpisd:jpisf,1:nntide),-1,status)
( 1363) 
( 1364)            !CALL flioclo (fid_s)
( 1365)          ENDIF
( 1366) 
( 1367)       ENDIF        !       END IF kt == nit000
( 1368) 
( 1369)      
( 1370) !!------------------------------------------------------------------------------------
( 1371)      ! 2.      Initialize the time we are at. Does this every time the routine is called,
( 1372)      !         excepted when nobc_dta = 0
( 1373)      !
( 1374)      IF( nobc_dta == 0) THEN
( 1375)         itimo = 1
( 1376)         zxy   = 0
( 1377)      ELSE
( 1378)         IF(itobc == 1) THEN
( 1379)            itimo = 1
( 1380)         ELSE IF (itobc == 12) THEN      !   BC are monthly






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 25

( 1381)            ! we assume we have climatology in that case
( 1382)            iman  = 12
( 1383)            i15   = nday / 16
( 1384)            imois = nmonth + i15 - 1
( 1385)            IF( imois == 0 )   imois = iman
( 1386)            itimo = imois
( 1387)         ELSE
( 1388)            IF(lwp) WRITE(numout,*) 'data other than constant or monthly',kt
( 1389)            iman  = itobc
( 1390)            itimo = FLOOR( kt*rdt / ztcobc(1))
( 1391)            isrel=kt*rdt
( 1392)         ENDIF
( 1393)      ENDIF
( 1394) 
( 1395)      ! 2. Read two records in the file if necessary
( 1396)      ! ---------------------------------------------
( 1397) 
( 1398)      IF( nobc_dta == 1 .AND. nlecto == 1 ) THEN
( 1399)  !sujie add  --------------
( 1400)       ! Calendar computation
( 1401)          IF( itobc == 1 )   THEN            !  BC are constant in time
( 1402)             ntobc1 = 1
( 1403)             ntobc2 = 1
( 1404)             ntobc3 = 1
( 1405)          ELSE IF( itobc == 12 )   THEN      !   BC are monthly
( 1406)             ntobc1 = itimo         ! first file record used
( 1407)             ntobc2 = ntobc1 + 1    ! second  file record used
( 1408)             ntobc3 = ntobc1 + 2    ! last  file record used
( 1409)             ntobc1 = MOD( ntobc1, iman )
( 1410)             IF( ntobc1 == 0 )   ntobc1 = iman
( 1411)             ntobc2 = MOD( ntobc2, iman )
( 1412)             IF( ntobc2 == 0 )   ntobc2 = iman
( 1413)             ntobc3 = MOD( ntobc3, iman )
( 1414)             IF( ntobc3 == 0 )   ntobc3 = iman
( 1415) 
( 1416)          ELSE
( 1417)             isrel=kt*rdt
( 1418)             ntobc1 = itimo         ! first file record used
( 1419)             ntobc2 = ntobc1 + 1    ! second  file record used
( 1420)             ntobc3 = ntobc1 + 2    ! last  file record used
( 1421)             ntobc1 = MOD( ntobc1, iman )
( 1422)             IF( ntobc1 == 0 )   ntobc1 = iman
( 1423)             ntobc2 = MOD( ntobc2, iman )
( 1424)             IF( ntobc2 == 0 )   ntobc2 = iman
( 1425)             ntobc3 = MOD( ntobc3, iman )
( 1426)             IF( ntobc3 == 0 )   ntobc3 = iman
( 1427)          ENDIF
( 1428) ! ------------------
( 1429) 
( 1430)         IF( lp_obc_east ) THEN
( 1431)            ! ... Read datafile and set sea surface height and barotropic velocity
( 1432)            ! ... initialise the sshedta, ubtedta arrays
( 1433)            sshedta(:,0) = sshedta(:,1)
( 1434)            ubtedta(:,0) = ubtedta(:,1)
( 1435) !!DB/CN: Replace IOIPSL with lib_ncdf
( 1436)            CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc1,pdta_2D=sshedta(jpjed:jpjef,1),fn='obceast_TS.nc')
( 1437)            CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc2,pdta_2D=sshedta(jpjed:jpjef,2),fn='obceast_TS.nc')
( 1438)            IF( lk_dynspg_ts ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 26

( 1439)               CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc3,pdta_2D=sshedta(jpjed:jpjef,3),fn='obceast_TS.nc')
( 1440)            ENDIF
( 1441)            CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc1,pdta_2D=ubtedta(jpjed:jpjef,1),fn='obceast_U.nc')
( 1442)            CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc2,pdta_2D=ubtedta(jpjed:jpjef,2),fn='obceast_U.nc')
( 1443)            IF( lk_dynspg_ts ) THEN
( 1444)               CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc3,pdta_2D=ubtedta(jpjed:jpjef,3),fn='obceast_U.nc')
( 1445)            ENDIF
( 1446) 
( 1447)         ENDIF
( 1448) 
( 1449)         IF( lp_obc_west ) THEN
( 1450)            ! ... Read datafile and set temperature, salinity and normal velocity
( 1451)            ! ... initialise the swdta, twdta, uwdta arrays
( 1452)            sshwdta(:,0) = sshwdta(:,1)
( 1453)            ubtwdta(:,0) = ubtwdta(:,1)
( 1454) !!DB/CN: Replace IOIPSL with lib_ncdf
( 1455)            CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc1,pdta_2D=sshwdta(jpjwd:jpjwf,1),fn='obcwest_TS.nc')
( 1456)            CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc2,pdta_2D=sshwdta(jpjwd:jpjwf,2),fn='obcwest_TS.nc')
( 1457)            IF( lk_dynspg_ts ) THEN
( 1458)               CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc3,pdta_2D=sshwdta(jpjwd:jpjwf,3),fn='obcwest_TS.nc')
( 1459)            ENDIF
( 1460)            CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc1,pdta_2D=ubtwdta(jpjwd:jpjwf,1),fn='obcwest_U.nc')
( 1461)            CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc2,pdta_2D=ubtwdta(jpjwd:jpjwf,2),fn='obcwest_U.nc')
( 1462)            IF( lk_dynspg_ts ) THEN
( 1463)               CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc3,pdta_2D=ubtwdta(jpjwd:jpjwf,3),fn='obcwest_U.nc')
( 1464)            ENDIF
( 1465) 
( 1466)         ENDIF
( 1467) 
( 1468)         IF( lp_obc_north) THEN
( 1469)            ! ... Read datafile and set sea surface height and barotropic velocity
( 1470)            ! ... initialise the sshndta, ubtndta arrays
( 1471)            sshndta(:,0) = sshndta(:,1)
( 1472)            vbtndta(:,0) = vbtndta(:,1)
( 1473) !!DB/CN: Replace IOIPSL with lib_ncdf
( 1474)            CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc1,pdta_2D=sshndta(jpind:jpinf,1),fn='obcnorth_TS.nc')
( 1475)            CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc2,pdta_2D=sshndta(jpind:jpinf,2),fn='obcnorth_TS.nc')
( 1476)            IF( lk_dynspg_ts ) THEN
( 1477)                CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc3,pdta_2D=sshndta(jpind:jpinf,3),fn='obcnorth_TS.nc')
( 1478)            ENDIF
( 1479)            CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc1,pdta_2D=vbtndta(jpind:jpinf,1),fn='obcnorth_V.nc')
( 1480)            CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc2,pdta_2D=vbtndta(jpind:jpinf,2),fn='obcnorth_V.nc')
( 1481)            IF( lk_dynspg_ts ) THEN
( 1482)               CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc3,pdta_2D=vbtndta(jpind:jpinf,3),fn='obcnorth_V.nc')
( 1483)            ENDIF
( 1484) 
( 1485)         ENDIF
( 1486) 
( 1487)         IF( lp_obc_south) THEN
( 1488)            ! ... Read datafile and set sea surface height and barotropic velocity
( 1489)            ! ... initialise the sshsdta, ubtsdta arrays
( 1490)            sshsdta(:,0) = sshsdta(:,1)
( 1491)            vbtsdta(:,0) = vbtsdta(:,1)
( 1492) !!DB/CN: Replace IOIPSL with lib_ncdf
( 1493)            CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc1,pdta_2D=sshsdta(jpisd:jpisf,1),fn='obcsouth_TS.nc')
( 1494)            CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc2,pdta_2D=sshsdta(jpisd:jpisf,2),fn='obcsouth_TS.nc')
( 1495)            IF( lk_dynspg_ts ) THEN
( 1496)                CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc3,pdta_2D=sshsdta(jpisd:jpisf,3),fn='obcsouth_TS.nc')






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 27

( 1497)            ENDIF
( 1498)            CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc1,pdta_2D=vbtsdta(jpisd:jpisf,1),fn='obcsouth_V.nc')
( 1499)            CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc2,pdta_2D=vbtsdta(jpisd:jpisf,2),fn='obcsouth_V.nc')
( 1500)            IF( lk_dynspg_ts ) THEN
( 1501)                CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc3,pdta_2D=vbtsdta(jpisd:jpisf,3),fn='obcsouth_V.nc')
( 1502)            ENDIF
( 1503) 
( 1504)         ENDIF
( 1505) 
( 1506)       ENDIF        !      end of the test on the condition to read or not the files
( 1507) 
( 1508)      ! 3.  Call at every time step : Linear interpolation of BCs to current time step
( 1509)      ! ----------------------------------------------------------------------
( 1510) 
( 1511)       IF( lk_dynspg_ts ) THEN
( 1512)          isrel = (kt-1)*rdt + kbt*rdtbt
( 1513) 
( 1514)          IF( nobc_dta == 1 ) THEN
( 1515)             isrel = (kt-1)*rdt + kbt*rdtbt
( 1516)             itimo  = FLOOR(  kt*rdt    / (ztcobc(2)-ztcobc(1)) )
( 1517)             itimom = FLOOR( (kt-1)*rdt / (ztcobc(2)-ztcobc(1)) )
( 1518)             itimop = FLOOR( (kt+1)*rdt / (ztcobc(2)-ztcobc(1)) )
( 1519)             IF( itimom == itimo .AND. itimop == itimo ) THEN
( 1520)                ntobcm = ntobc1
( 1521)                ntobcp = ntobc2
( 1522) 
( 1523)             ELSEIF ( itimom <= itimo .AND. itimop == itimo ) THEN
( 1524)                IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimo ) THEN
( 1525)                   ntobcm = ntobc1-1
( 1526)                   ntobcp = ntobc2-1
( 1527)                ELSE
( 1528)                   ntobcm = ntobc1
( 1529)                   ntobcp = ntobc2
( 1530)                ENDIF
( 1531) 
( 1532)             ELSEIF ( itimom == itimo .AND. itimop >= itimo ) THEN
( 1533)                IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimop ) THEN
( 1534)                   ntobcm = ntobc1
( 1535)                   ntobcp = ntobc2
( 1536)                ELSE
( 1537)                   ntobcm = ntobc1+1
( 1538)                   ntobcp = ntobc2+1
( 1539)                ENDIF
( 1540) 
( 1541)             ELSEIF ( itimom == itimo-1 .AND. itimop == itimo+1 ) THEN
( 1542)                IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimo ) THEN
( 1543)                   ntobcm = ntobc1-1
( 1544)                   ntobcp = ntobc2-1
( 1545)                ELSEIF (  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimop ) THEN
( 1546)                   ntobcm = ntobc1
( 1547)                   ntobcp = ntobc2
( 1548)                ELSEIF (  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) == itimop ) THEN
( 1549)                   ntobcm = ntobc1+1
( 1550)                   ntobcp = ntobc2+2
( 1551)                ELSE
( 1552)                   IF(lwp) WRITE(numout, *) 'obc_dta_bt: You should not have seen this print! error 1?'
( 1553)                ENDIF
( 1554)             ELSE






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 28

( 1555)                IF(lwp) WRITE(numout, *) 'obc_dta_bt: You should not have seen this print! error 2?'
( 1556)             ENDIF
( 1557) 
( 1558)          ENDIF
( 1559) 
( 1560)       ELSE IF( lk_dynspg_exp ) THEN
( 1561)          isrel=kt*rdt
( 1562)          ntobcm = ntobc1
( 1563)          ntobcp = ntobc2
( 1564)       ENDIF
( 1565) 
( 1566)       IF( itobc == 1 .OR. nobc_dta == 0 ) THEN
( 1567)          zxy = 0.e0
( 1568)       ELSE IF( itobc == 12 ) THEN
( 1569)          zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
( 1570)       ELSE
( 1571)          zxy = (ztcobc(ntobcm)-FLOAT(isrel)) / (ztcobc(ntobcm)-ztcobc(ntobcp))
( 1572)       ENDIF
( 1573) 
( 1574) 
( 1575) !DB: 2007.12.27 -- Note that I now keep the mean components separate
( 1576) !    Formulae should work for any number of tide components
( 1577) !    REM: velocity components already zxy-weighted & multiplied by ramp
( 1578)       IF( lp_obc_east ) THEN    
( 1579) 
( 1580)          do ji = nie0, nie1
( 1581)             do jj = nje0p1, nje1m1
( 1582) !!DB: REM that ufoe has been multiplied by ramp
( 1583)                ubtfoe0(jj) = 0.0
( 1584)                do jk = 1, jpkm1
( 1585)                   ubtfoe0(jj) = ubtfoe0(jj) + ufoe(jj,jk) * e3t(jk)
( 1586)                enddo
( 1587)                ij = jj -1 + njmpp
( 1588) !!DB 
( 1589)                sshfoe0(jj) = eta_e(ij) * temsk(jj,1)
( 1590)                ubtfoe0(jj) = ubtfoe0(jj) * uemsk(jj,1) 
( 1591)                ubtfoe(jj) = 0.0
( 1592)                sshfoe(jj) = 0.0
( 1593) 
( 1594) #ifdef key_ADJ_TRANSPORT
( 1595) !!DB: Keep this ifdef open
( 1596) 
( 1597) #endif
( 1598)                do ntide = 1, nntide
( 1599)                   sshfoe(jj) = sshfoe(jj) + tidesshemag(ij,ntide)  &
( 1600)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1601)                        - tidesshepha(ij,ntide) ) * temsk(jj,1)
( 1602) 
( 1603)                   ubtfoe(jj) = ubtfoe(jj) &
( 1604)                        + tidevbtemag(ij,ntide)  &
( 1605)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1606)                        - tidevbtepha(ij,ntide) ) * uemsk(jj,1)
( 1607)                enddo
( 1608) !!DB 
( 1609)                sshfoe(jj) = ramp*(sshfoe0(jj) + sshfoe(jj))
( 1610)                ubtfoe(jj) = ramp*ubtfoe(jj) + ubtfoe0(jj)
( 1611) 
( 1612)             enddo






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 29

( 1613)          enddo
( 1614) 
( 1615) !-----------------------------------------------
( 1616) 
( 1617)       ENDIF
( 1618) 
( 1619) !DB: 2007.12.27 -- Note that I now keep the mean components separate
( 1620) !    Formulae should work for any number of tide components
( 1621) !    REM: velocity components already zxy-weighted & multiplied by ramp
( 1622)       IF( lp_obc_west) THEN   
( 1623) 
( 1624)          do ji = niw0, niw1
( 1625)             do jj = njw0p1, njw1m1
( 1626)                ubtfow0(jj) = 0.0 
( 1627)                do jk = 1, jpkm1
( 1628)                   ubtfow0(jj) = ubtfow0(jj) + ufow(jj,jk) * e3t(jk)
( 1629)                enddo
( 1630)                ij = jj -1 + njmpp
( 1631) 
( 1632)                ubtfow0(jj) = ubtfow0(jj) * uwmsk(jj,1) 
( 1633)                ubtfow(jj) = 0.0
( 1634)                sshfow0(jj) = eta_w(ij)*twmsk(jj,1)
( 1635)                sshfow(jj) = 0.0
( 1636) 
( 1637)                do ntide = 1, nntide
( 1638)                   sshfow(jj) =  sshfow(jj) + tidesshwmag(ij,ntide)  &
( 1639)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1640)                        - tidesshwpha(ij,ntide) ) * twmsk(jj,1) 
( 1641)                   
( 1642)                   ubtfow(jj) = ubtfow(jj) &
( 1643)                        + tidevbtwmag(ij,ntide)  &
( 1644)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1645)                        - tidevbtwpha(ij,ntide) ) * uwmsk(jj,1)
( 1646)                enddo
( 1647)                sshfow(jj) =  ramp*(sshfow(jj) + sshfow0(jj))
( 1648)                ubtfow(jj) = ramp*ubtfow(jj)+ ubtfow0(jj) 
( 1649)             enddo
( 1650)          enddo
( 1651) 
( 1652) !------------------------------------------------
( 1653)       ENDIF
( 1654) 
( 1655) !---------------------------------------------------------------------------
( 1656) !byoung for St. Lawrence River runoff: DELETED -- see old code
( 1657) 
( 1658) !----------------------------------------------------------------------------
( 1659) 
( 1660) !!DB: north should be closed. If not, modify this section
( 1661)       IF( lp_obc_north) THEN           !  fills sshfon, vbtfon (local to each processor)
( 1662)          DO ji = nin0p1, nin1m1
( 1663)             ii = ji -1 + nimpp
( 1664)             sshfon(ji) = ( zxy * sshndta(ii,2) + (1.-zxy) * sshndta(ii,1) ) * tnmsk(ji,1)
( 1665)             vbtfon(ji) = ( zxy * vbtndta(ii,2) + (1.-zxy) * vbtndta(ii,1) ) * vnmsk(ji,1)
( 1666) !ylu
( 1667)           do ntide = 1,4  
( 1668)             sshfon(ji) = sshfon(ji) + ramp*tidesshnmag(ii,ntide)  &
( 1669)                          * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1670)                          - tidesshnpha(ii,ntide) ) * tnmsk(ji,1) 






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 30

( 1671)             vbtfon(ji) = vbtfon(ji) + ramp*tidevbtnmag(ii,ntide)  &
( 1672)                          * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1673)                          - tidevbtnpha(ii,ntide) ) * vnmsk(ji,1) 
( 1674)           enddo
( 1675)          END DO
( 1676)       ENDIF
( 1677) 
( 1678) 
( 1679) !DB: 2007.12.27 -- Note that I now keep the mean components separate
( 1680) !    Formulae should work for any number of tide components
( 1681) !    REM: velocity components already zxy-weighted & multiplied by ramp
( 1682)       IF( lp_obc_south) THEN     
( 1683) 
( 1684)          do jj = njs0, njs1
( 1685)             do ji = nis0p1, nis1m1
( 1686)                vbtfos0(ji) = 0.0 
( 1687)                do jk = 1, jpkm1
( 1688)                   vbtfos0(ji) = vbtfos0(ji) + vfos(ji,jk) * e3t(jk)
( 1689)                enddo
( 1690)                ii = ji -1 + nimpp
( 1691) 
( 1692)                vbtfos0(ji) = vbtfos0(ji)*vsmsk(ji,1)
( 1693)                vbtfos(ji) = 0.0
( 1694)                sshfos0(ji) = eta_s(ii)* tsmsk(ji,1)
( 1695)                sshfos(ji) = 0.0
( 1696) 
( 1697)                do ntide = 1, nntide
( 1698)                   sshfos(ji) = sshfos(ji) + tidesshsmag(ii,ntide)  &
( 1699)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1700)                        - tidesshspha(ii,ntide) ) * tsmsk(ji,1)
( 1701) 
( 1702)                   vbtfos(ji) =  vbtfos(ji) + & 
( 1703)                        tidevbtsmag(ii,ntide)  &
( 1704)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1705)                        - tidevbtspha(ii,ntide) ) * vsmsk(ji,1) 
( 1706)                enddo
( 1707)                sshfos(ji) = ramp*(sshfos(ji) + sshfos0(ji))
( 1708)                vbtfos(ji) =  ramp*vbtfos(ji) +  vbtfos0(ji)
( 1709)             enddo
( 1710)          enddo
( 1711) 
( 1712)       ENDIF
( 1713) 
( 1714)   END SUBROUTINE obc_dta_bt
( 1715) 
( 1716) 
( 1717) #else
( 1718)   !!-----------------------------------------------------------------------------
( 1719)   !!   Default option
( 1720)   !!-----------------------------------------------------------------------------
( 1721)   SUBROUTINE obc_dta_bt ( kt, kbt )       ! Empty routine
( 1722)      !! * Arguments
( 1723)      INTEGER,INTENT(in) :: kt
( 1724)      INTEGER, INTENT( in ) ::   kbt         ! barotropic ocean time-step index
( 1725)      if(lwp) WRITE(numout,*) 'obc_dta_bt: You should not have seen this print! error?', kt
( 1726)      if(lwp) WRITE(numout,*) 'obc_dta_bt: You should not have seen this print! error?', kbt
( 1727)   END SUBROUTINE obc_dta_bt
( 1728) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 31

( 1729) 
( 1730) !!DB/CN 2008.11
( 1731)   SUBROUTINE obc_dta_gv (cldim,clobc,kobcij,ktobc,pdta_2D,pdta_3D,pdta_4D,fn)
( 1732)      !!-----------------------------------------------------------------------------
( 1733)      !!                       ***  SUBROUTINE obc_dta_gv  ***
( 1734)      !!
( 1735)      !! ** Purpose :   Read an OBC forcing field from netcdf file
( 1736)      !!                Input file are supposed to be 3D e.g.
( 1737)      !!                - for a South or North OB : longitude x depth x time
( 1738)      !!		- for a West or East OB : latitude x depth x time
( 1739)      !!
( 1740)      !! History :
( 1741)      !!        !  04-06 (A.-M. Treguier, F. Durand) Original code
( 1742)      !!        !  05-02 (J. Bellier, C. Talandier) use fliocom CALL
( 1743)      !!----------------------------------------------------------------------------
( 1744)      !! * Arguments
( 1745)      INTEGER, INTENT(IN) ::   &
( 1746) !        ifid  ,               & ! netcdf file name identifier
( 1747)         kobcij,               & ! Horizontal (i or j) dimension of the array
( 1748)         ktobc                   ! starting time index read
( 1749)      CHARACTER(LEN=*), INTENT(IN)    ::   &
( 1750)         cldim,                & ! dimension along which is the open boundary ('x' or 'y')
( 1751)         clobc                   ! name of the netcdf variable read
( 1752)      REAL, DIMENSION(kobcij,jpk,1), INTENT(OUT), OPTIONAL ::   &
( 1753)         pdta_3D                 ! 3D array of OBC forcing field
( 1754) 
( 1755)      REAL, DIMENSION(kobcij,1), INTENT(OUT), OPTIONAL ::   &
( 1756)         pdta_2D                 ! 2D array of OBC forcing field
( 1757)      REAL, DIMENSION(kobcij,jpk,5,1), INTENT(OUT), OPTIONAL ::   &
( 1758)         pdta_4D                 ! 4D array of OBC forcing field
( 1759) !!CN
( 1760)      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fn ! name of netCDF file
( 1761)      LOGICAL :: fnswitch ! logical indicating whether to use filename & lib_ncdf or not
( 1762)                           ! to be removed once IOIPSL is entirely gone from obcdta
( 1763)      REAL, DIMENSION(kobcij) :: buf_1D             ! Needed for lib_ncdf calls
( 1764)      REAL, DIMENSION(kobcij,jpk) :: buf_2D         ! Needed for lib_ncdf calls
( 1765)      REAL, DIMENSION(kobcij,jpk,5) :: buf_3D       ! Needed for lib_ncdf calls
( 1766)      INTEGER :: f_stat
( 1767)      LOGICAL :: tide
( 1768) 
( 1769)      !! * Local declarations
( 1770)      INTEGER ::   indim
( 1771)      LOGICAL ::   l_exv
( 1772)      INTEGER,DIMENSION(4) ::   f_d, istart, icount
( 1773)      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   v_tmp_4
( 1774)      !----------------------------------------------------------------------
( 1775) 
( 1776)      !!CN: Replacing IOIPSL calls with lib_ncdf
( 1777)      fnswitch = .TRUE.
( 1778)      CALL ncdf_get_num_dims(fn, TRIM(clobc), indim, f_stat);
( 1779) 
( 1780)      l_exv = .TRUE.
( 1781)      IF( l_exv )   THEN
( 1782)         ! checks the number of dimensions
( 1783)         IF( indim == 2 )   THEN
( 1784)            istart(1:2) = (/ 1     , ktobc /)
( 1785)            icount(1:2) = (/ kobcij, 1     /)
( 1786)            CALL ncdf_read(fn, TRIM(clobc), buf_1D, f_stat)






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 32

( 1787)            pdta_2D(:,1) = buf_1D
( 1788)         ELSE IF( indim == 3 )   THEN
( 1789)            istart(1:3) = (/ 1     , 1    , ktobc /)
( 1790)            icount(1:3) = (/ kobcij, jpk  , 1     /)
( 1791)            CALL ncdf_read_global(fn, TRIM(clobc), buf_2D, -ktobc, f_stat)
( 1792)            pdta_3D(:,:,1) = buf_2D
( 1793)         ELSE IF( indim == 4 )   THEN
( 1794)            istart(1:4) = (/ 1, 1, 1, ktobc /)
( 1795)            icount(1:4) = (/ kobcij, jpk  , 5, 1 /)
( 1796)            CALL ncdf_read_global(fn, TRIM(clobc), buf_3D, -ktobc, f_stat)
( 1797)            pdta_4D(:,:,:,1) = buf_3d
( 1798)         ELSE
( 1799)            IF( lwp )   THEN
( 1800)               WRITE(numout,*) ' Problem in OBC file for ',TRIM(clobc),' :'
( 1801)               WRITE(numout,*) ' number of dimensions (not 3 or 4) =',indim
( 1802)            ENDIF
( 1803)            STOP
( 1804)         ENDIF
( 1805)      ELSE
( 1806)         if(lwp) WRITE(numout,*) ' Variable ',TRIM(clobc),' not found'
( 1807)      ENDIF
( 1808) 
( 1809)   END SUBROUTINE obc_dta_gv
( 1810) 
( 1811) #else
( 1812)   !!--------------------------------------------------------------------
( 1813)   !!  default option  :  Dummy module    NO Open Boundary Conditions
( 1814)   !!--------------------------------------------------------------------
( 1815) CONTAINS
( 1816)   SUBROUTINE obc_dta( kt )             ! Dummy routine
( 1817)     INTEGER, INTENT (in) :: kt
( 1818)     if(lwp) WRITE(numout,*) 'obc_dta: You should not have seen this print! error?', kt
( 1819)   END SUBROUTINE obc_dta
( 1820)   SUBROUTINE obc_dta_bt( kt, jn)             ! Dummy routine
( 1821)     INTEGER, INTENT (in) :: kt, jn
( 1822)     if(lwp) WRITE(numout,*) 'obc_dta_bt: You should not have seen this print! error?', kt
( 1823)     if(lwp) WRITE(numout,*) 'obc_dta_bt: You should not have seen this print! error?', jn
( 1824)   END SUBROUTINE obc_dta_bt
( 1825) #endif
( 1826) 
( 1827) 
( 1828)   !!=====================================================================
( 1829) END MODULE obcdta





















PGF90 (Version     10.5)          09/23/2010  09:49:39      page 33

(    1) # 1 "obcdta.F90"
(    1) !!DB/CN: 2008.11.14
(    2) !!Calls to IOIPSL eliminated and tested
(    3) !!More code cleanup likely possible, but not done
(    4) 
(    5) !!Key changes from standard obcdta.F90
(    6) !!(0) Look for  DB  to see various modifications
(    7) 
(    8) !!(1) obc_dta() has HARDWIRED areas where transports at certain boundary
(    9) !!locations are adjusted ---> 1; SBI_trans; SB_east; ...
(   10) !!- Note that at this time the code is not completely consistent with
(   11) !!respect to 1 as some variables and operations exist that
(   12) !!are only required for that key but they are done even if that key is not defined.
(   13) !!These exceptions require minimal memory and execution time.
(   14) 
(   15) !!(2) Because OBC transports are adjusted, the eta-arrays need to be adjusted as
(   16) !!well to be consistent with the OBC vels (thermal wind is used). This is a major
(   17) !!change that requires new variables (search eta_). 
(   18) !!Note that 1 & 2 require definitions and computations on global variables defined
(   19) !!for these new OBC modifications. By global I mean that all relevant CPUs must
(   20) !!do the same computations and have the same values.
(   21) 
(   22) !!(3) The barotropic velocity in obc_dta_bt() is now computed as the vertical
(   23) !!integral of the baroclinic velocities. Thus the BT vels in the input files 
(   24) !!are not used anymore (they are ignored). Search (e.g.) ubtfoe0
(   25) !!(4) 5 tidal constituents are possible (search: nntide)
(   26) !!(5) tidal and non-tidal forcing are separated (search (e.g.) ubtfoe0, sshfoe0) 
(   27) !!(6) A routine is called that zeros the non-tidal net transport around the
(   28) !! boundaries (search: obc_ctl) 
(   29) !!(7) A ramp function is used for the forcing (search ramp). It is computed in
(   30) !!step.F90 and defined in oce.F90, It is based on kt as
(   31) !!opposed to (kt-nit000) so that by default it does not ramp a restart.  Also
(   32) !!it is complicated to follow when/where certain variables are ramped up -- but no
(   33) !!apologies for this.
(   34) !!(8) The routines are HARDWIRED for 3 open boundaries only. If north is also open, then
(   35) !!all the lp_obc_north areas should be checked and the code added as necessary. 
(   36) !!(9) Uses global masks: emaskg2, smaskg2, wmaskg2 (see obcini, obc_oce) 
(   37) !!(10) NB: use of e2v_e, e2u_s etc seem incorrect, but diff should be v.small so
(   38) !!     modification ... TO DO ... 
(   39) 
(   40) 
(   41) MODULE obcdta 
(   42)   !!==============================================================================
(   43)   !!                            ***  MODULE obcdta  ***
(   44)   !! Open boundary data : read the data for the open boundaries.
(   45)   !!==============================================================================
(   46) # 47
(   47)   !!------------------------------------------------------------------------------
(   48)   !!   'key_obc'         :                                Open Boundary Conditions
(   49)   !!------------------------------------------------------------------------------
(   50)   !!   obc_dta           : read u, v, t, s data along each open boundary
(   51)   !!   obc_dta_psi       : read psi data along each open boundary (rigid lid only)
(   52)   !!------------------------------------------------------------------------------
(   53)   !! * Modules used
(   54)   USE oce             ! ocean dynamics and tracers
(   55)   USE dom_oce         ! ocean space and time domain
(   56)   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   57)   USE phycst          ! physical constants






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 34

(   58)   USE obc_oce         ! ocean open boundary conditions
(   59)   USE daymod          ! calendar
(   60)   USE in_out_manager  ! I/O logical units
(   61)   USE lib_mpp         ! distributed memory computing
(   62)   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   63)   USE ioipsl
(   64)   USE obcctl
(   65) !!CN
(   66)   USE lib_ncdf
(   67) !!DBG
(   68) 
(   69) # 72
(   72)   IMPLICIT NONE
(   73)   PRIVATE
(   74) 
(   75)   !! * Accessibility
(   76)   PUBLIC obc_dta        ! routines called by step.F90
(   77)   PUBLIC obc_dta_bt     ! routines called by dynspg_ts.F90
(   78) 
(   79)   !! * Shared module variables
(   80)   INTEGER ::   &
(   81)      nlecto,   &  ! switch for the first read
(   82)      ntobc1,   &  ! first record used
(   83)      ntobc2,   &  ! second record used
(   84)      ntobc3,   &  ! last record used  for dynspg_ts
(   85)      itobc        ! number of time steps in OBC files
(   86) 
(   87) 
(   88)   REAL(wp), DIMENSION(:), ALLOCATABLE :: ztcobc      ! time_counter variable of BCs
(   89) 
(   90)   !! * Substitutions
(   91) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 92 "obcdta.F90"
(   92) # 1 "./obc_vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                 ***  obc_vectopt_loop_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obc_vectopt_loop_substitute.h90,v 1.3 2005/12/12 14:20:26 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 35

(   11)    !!----------------------------------------------------------------------
(   12) # 93 "obcdta.F90"
(   93) 
(   94) # 94
(   94)   !!---------------------------------------------------------------------------------
(   95)   !!   OPA 9.0 , LODYC-IPSL  (2003)
(   96)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcdta.F90,v 1.9 2006/03/21 08:25:09 opalod Exp $
(   97)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   98)   !!---------------------------------------------------------------------------------
(   99) 
(  100) CONTAINS
(  101) 
(  102)   SUBROUTINE obc_dta (kt)
(  103)      !!--------------------------------------------------------------------
(  104)      !!              ***  SUBROUTINE obc_dta  ***
(  105)      !!
(  106)      !! ** Purpose :
(  107)      !!   Find the climatological boundary arrays for the specified date,
(  108)      !!   The boundary arrays are netcdf files. Three possible cases:
(  109)      !!   - one time frame only in the file (time dimension = 1).
(  110)      !!     in that case the boundary data does not change in time.
(  111)      !!   - many time frames. In that case,  if we have 12 frames
(  112)      !!     we assume monthly fields.
(  113)      !!     Else, we assume that time_counter is in seconds
(  114)      !!     since the beginning of either the current year or a reference
(  115)      !!     year given in the namelist.
(  116)      !!     (no check is done so far but one would have to check the "unit"
(  117)      !!     attribute of variable time_counter).
(  118)      !!
(  119)      !! History :
(  120)      !!        !  98-05 (J.M. Molines) Original code
(  121)      !!   8.5  !  02-10 (C. Talandier, A-M. Treguier) Free surface, F90
(  122)      !!   9.0  !  04-06 (F. Durand, A-M. Treguier) Netcdf BC files on input
(  123)      !!--------------------------------------------------------------------
(  124)      !! * Arguments
(  125) 
(  126) 
(  127)      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(  128) 
(  129)      !! * Local declarations
(  130)      INTEGER ::   ji, jj, jk, ii, ij   ! dummy loop indices
(  131)      INTEGER ::   itimo, iman, imois
(  132)      INTEGER ::   i15
(  133)      REAL(wp) ::   zxy
(  134)      !! * Ajouts FD
(  135)      INTEGER ::  isrel              ! number of seconds since 1/1/1992
(  136)      INTEGER, SAVE ::  itobce, itobcw,  & ! number of time steps in OBC files
(  137)                        itobcs, itobcn     !    "       "       "       "
(  138)      INTEGER ::  ikprint        ! frequency for printouts.
(  139)      INTEGER :: fid_e, fid_w, fid_n, fid_s       ! file identifiers
(  140)      LOGICAL :: l_exv
(  141)      INTEGER, DIMENSION(flio_max_dims) ::   f_d  ! dimensions lenght
(  142)      CHARACTER(LEN=25) :: v_name
(  143) !!CN
(  144)      INTEGER :: f_stat
(  145)      INTEGER :: len
(  146)      !!LOGICAL :: use_ioipsl = .FALSE.
(  147) 






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 36

(  148) !!DB
(  149)      LOGICAL :: DBG
(  150)      INTEGER :: ncells
(  151)      REAL(wp), DIMENSION(jpjdta,jpk) :: U_EW
(  152) !     REAL(wp), DIMENSION(jpidta,jpk) :: V_NS     !!for future use
(  153) 
(  154) !!DB: For boundary transport adjustments
(  155) !!NB: Could be inside of #ifdef 1
(  156)      REAL(wp) ::   fac_SBI, fac_SB_east, trans_in, trans_diff, trans_check, & 
(  157)                     vfreq, off1,a1,ph1, off2,a2,ph2  !!for time-dependent fac_*
(  158)      INTEGER :: jwest1, jwest2, jeast1, jeast2                !!indices where transport is adjusted 
(  159) 
(  160) 
(  161)      !!--------------------------------------------------------------------
(  162) 
(  163) !!DB
(  164)      DBG = .false. 
(  165) 
(  166)      IF( lk_dynspg_rl )  THEN
(  167)         CALL obc_dta_psi (kt)     ! update bsf data at open boundaries
(  168)         IF( nobc_dta == 1 .AND. kt == nit000 )   THEN
(  169)            IF(lwp) WRITE(numout,*) ' time-variable psi boundary data not allowed yet'
(  170)            STOP
(  171)         ENDIF
(  172)      ENDIF
(  173) 
(  174)      CALL ipslnlf (new_number=numout)
(  175) 
(  176) !!DB initialize time-dependent valve forcing parameters
(  177) !!Note that these are only used if 1 is defined, but I
(  178) !!and lazy and do not #ifdef this
(  179) 
(  180)      vfreq = 2.0*rpi/365.0   ! seasonal cycle in days
(  181) !SBI params
(  182)      off1 = 3.0
(  183)      a1 = 2.0
(  184)      ph1 = 30.0   ! in days
(  185) !SB_east params: CLM0 values
(  186)      off2 = 3.0
(  187)      a2 = 0.5
(  188)      ph2 = -60.0   ! in days
(  189) 
(  190)      ! 1.   First call: check time frames available in files.
(  191)      ! -------------------------------------------------------
(  192) 
(  193)      IF( kt == nit000 )   THEN
(  194) 
(  195)         nlecto =  0
(  196) 
(  197)         IF(lwp) WRITE(numout,*)
(  198)         IF(lwp) WRITE(numout,*)     'obc_dta : find boundary data'
(  199)         IF(lwp) WRITE(numout,*)     '~~~~~~~'
(  200) 
(  201)         IF( nobc_dta == 0 )   THEN
(  202)            IF(lwp) WRITE(numout,*)  '  OBC data taken from initial conditions.'
(  203)            ntobc1 = 1
(  204)            ntobc2 = 1
(  205)         ELSE






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 37

(  206)            IF(lwp) WRITE(numout,*)  '  OBC data taken from netcdf files.'
(  207)            IF(lwp) WRITE(numout,*)  '  climatology (T/F):',ln_obc_clim
(  208)            ! check the number of time steps in the files.
(  209)            itobce =0 ; itobcw = 0; itobcn = 0; itobcs = 0
(  210)            v_name = 'time_counter'
(  211)            IF( lp_obc_east )   THEN
(  212)               !!CN: Replacing IOIPSL calls with lib_ncdf
(  213)               !CALL flioopfd ('obceast_TS.nc',fid_e)
(  214)               !CALL flioinqv (fid_e,TRIM(v_name),l_exv,len_dims=f_d)
(  215)               CALL ncdf_get_dim_size('obceast_TS.nc', 'time_counter', len, f_stat)
(  216)               IF( f_stat == 0 )   THEN
(  217)                  f_d(1) = len
(  218)                  itobce = f_d(1)
(  219)               ELSE
(  220)                  if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obceast_TS.nc'
(  221)               ENDIF
(  222)            ENDIF
(  223)            IF( lp_obc_west )   THEN
(  224)               !!CN: Replacing IOIPSL calls with lib_ncdf
(  225)               !CALL flioopfd ('obcwest_TS.nc',fid_w)
(  226)               !CALL flioinqv (fid_w,TRIM(v_name),l_exv,len_dims=f_d)
(  227)               CALL ncdf_get_dim_size('obcwest_TS.nc', 'time_counter', len, f_stat)
(  228)               IF( f_stat == 0 )   THEN
(  229)                  f_d(1) = len
(  230)                  itobcw = f_d(1)
(  231)               ELSE
(  232)                 if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcwest_TS.nc'
(  233)               ENDIF
(  234)            ENDIF
(  235)            IF( lp_obc_north )   THEN
(  236)               !!CN: Replacing IOIPSL calls with lib_ncdf
(  237)               !CALL flioopfd ('obcnorth_TS.nc',fid_n)
(  238)               !CALL flioinqv (fid_n,TRIM(v_name),l_exv,len_dims=f_d)
(  239)               CALL ncdf_get_dim_size('obcnorth_TS.nc', 'time_counter', len, f_stat)
(  240)               IF( f_stat == 0 )   THEN
(  241)                  f_d(1) = len
(  242)                  itobcn = f_d(1)
(  243)               ELSE
(  244)                  if(lwp)WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcnorth_TS.nc'
(  245)               ENDIF
(  246)            ENDIF
(  247)            IF( lp_obc_south )   THEN
(  248)               !!CN: Replacing IOIPSL calls with lib_ncdf
(  249)               !CALL flioopfd ('obcsouth_TS.nc',fid_s)
(  250)               !CALL flioinqv (fid_s,TRIM(v_name),l_exv,len_dims=f_d)
(  251)               CALL ncdf_get_dim_size('obcsouth_TS.nc', 'time_counter', len, f_stat)
(  252)               IF( f_stat == 0 )   THEN
(  253)                  f_d(1) = len
(  254)                  itobcs = f_d(1)
(  255)               ELSE
(  256)                  if(lwp) WRITE(numout,*) ' Variable ',TRIM(v_name),' not found in file ','obcsouth_TS.nc'
(  257)               ENDIF
(  258)            ENDIF
(  259) 
(  260)            itobc = MAX(itobce,itobcw,itobcn,itobcs)
(  261)            nstop = 0
(  262)            IF( lp_obc_east  .AND. itobce /= itobc ) nstop = nstop+1
(  263)            IF( lp_obc_west  .AND. itobcw /= itobc ) nstop = nstop+1






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 38

(  264)            IF( lp_obc_north .AND. itobcn /= itobc ) nstop = nstop+1
(  265)            IF( lp_obc_south .AND. itobcs /= itobc ) nstop = nstop+1
(  266)            IF( nstop /= 0 )  THEN
(  267)               IF( lwp )   THEN
(  268)                  WRITE(numout,*) ' obcdta : all files must have the same number of time steps'
(  269)                  WRITE(numout,*) ' east, west, north, south: ', itobce, itobcw, itobcn, itobcs
(  270)               ENDIF
(  271)               STOP
(  272)            ENDIF
(  273)            IF( itobc == 1 )   THEN
(  274)               IF( lwp ) WRITE(numout,*) ' obcdta found one time step only in the OBC files'
(  275)            ELSE
(  276)               ALLOCATE (ztcobc(itobc))
(  277)               l_exv = .TRUE.
(  278)               IF( lp_obc_east )   THEN
(  279)                  !!CN: Replacing IOIPSL calls with lib_ncdf
(  280)                  IF( l_exv )   THEN
(  281)                     !CALL fliogetv (fid_e,TRIM(v_name),ztcobc)
(  282)                     CALL ncdf_read('obceast_TS.nc', 'time_counter', ztcobc, f_stat)
(  283)                     l_exv = .FALSE.
(  284)                  ENDIF
(  285)                  !CALL flioclo (fid_e)
(  286)               ENDIF
(  287)               IF( lp_obc_west )   THEN
(  288)                  !!CN: Replacing IOIPSL calls with lib_ncdf
(  289)                 IF( l_exv )   THEN
(  290)                    !CALL fliogetv (fid_w,TRIM(v_name),ztcobc)
(  291)                    CALL ncdf_read('obcwest_TS.nc', 'time_counter', ztcobc, f_stat)
(  292)                    l_exv = .FALSE.
(  293)                 ENDIF
(  294)                 !CALL flioclo (fid_w)
(  295)               ENDIF
(  296)               IF( lp_obc_north )   THEN
(  297)                  !!CN: Replacing IOIPSL calls with lib_ncdf
(  298)                 IF( l_exv )   THEN
(  299)                    !CALL fliogetv (fid_n,TRIM(v_name),ztcobc)
(  300)                    CALL ncdf_read('obcnorth_TS.nc', 'time_counter', ztcobc, f_stat)
(  301)                    l_exv = .FALSE.
(  302)                 ENDIF
(  303)                 !CALL flioclo (fid_n)
(  304)               ENDIF
(  305)               IF( lp_obc_south )   THEN
(  306)                  !!CN: replacing IOIPSL calls with lib_ncdf
(  307)                 IF( l_exv )   THEN
(  308)                    !CALL fliogetv (fid_s,TRIM(v_name),ztcobc)
(  309)                    CALL ncdf_read('obcsouth_TS.nc', 'time_counter', ztcobc, f_stat)
(  310)                    l_exv = .FALSE.
(  311)                 ENDIF
(  312)                 !CALL flioclo (fid_s)
(  313)               ENDIF
(  314)               IF( lwp ) WRITE(numout,*) ' obcdta found', itobc,' time steps in the OBC files'
(  315)               IF( .NOT. ln_obc_clim .AND. itobc == 12 )   THEN
(  316)                  IF ( lwp ) WRITE(numout,*) '  WARNING: With monthly data we assume climatology'
(  317)                  ln_obc_clim = .true.
(  318)               ENDIF
(  319)            ENDIF
(  320)         ENDIF
(  321) 






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 39

(  322)       ! 1.1  Tangential velocities set to zero
(  323)       ! --------------------------------------
(  324)         IF( lp_obc_east  ) vfoe = 0.e0
(  325)         IF( lp_obc_west  ) vfow = 0.e0
(  326)         IF( lp_obc_south ) ufos = 0.e0
(  327)         IF( lp_obc_north ) ufon = 0.e0
(  328) 
(  329)       ! 1.2  Data temperature, salinity, normal velocities set to zero
(  330)       !                        or initial conditions if nobc_dta == 0
(  331)       ! --------------------------------------------------------------
(  332) 
(  333)         IF( lp_obc_east )   THEN
(  334)            ! initialisation to zero
(  335) !byoung           
(  336)            sedta(:,:,:,:) = 0.e0
(  337)            tedta(:,:,:,:) = 0.e0
(  338)            uedta(:,:,:) = 0.e0
(  339)            !                                    ! ================== !
(  340)            IF( nobc_dta == 0 )   THEN           ! initial state used
(  341)               !                                 ! ================== !
(  342)               !  Fills sedta, tedta, uedta (global arrays)
(  343)               !  Remark: this works for njzoom = 1.
(  344)               !          Should the definition of ij include njzoom?
(  345)               DO ji = nie0, nie1
(  346)                  DO jk = 1, jpkm1
(  347)                     DO jj = nje0p1, nje1m1
(  348)                        ij = jj -1 + njmpp
(  349) !byoung                       
(  350)                        sedta(ij,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  351)                        tedta(ij,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  352)                        uedta(ij,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
(  353)                     END DO
(  354)                  END DO
(  355)               END DO
(  356) !byoung
(  357)               DO ji = nie0-4, nie0-1
(  358)                  DO jk = 1, jpkm1
(  359)                     DO jj = nje0p1, nje1m1
(  360)                        ij = jj -1 + njmpp
(  361)                        sedta(ij,jk,nie0-ji+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  362)                        tedta(ij,jk,nie0-ji+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  363)                     END DO
(  364)                  END DO
(  365)               END DO              
(  366)            ENDIF
(  367)         ENDIF
(  368) 
(  369)         IF( lp_obc_west )   THEN
(  370)            ! initialisation to zero
(  371) !byoung           
(  372)            swdta(:,:,:,:) = 0.e0
(  373)            twdta(:,:,:,:) = 0.e0
(  374)            uwdta(:,:,:) = 0.e0
(  375)            !                                    ! ================== !
(  376)            IF( nobc_dta == 0 )   THEN           ! initial state used !
(  377)               !                                 ! ================== !
(  378)               !  Fills swdta, twdta, uwdta (global arrays)
(  379)               !  Remark: this works for njzoom = 1.






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 40

(  380)               !          Should the definition of ij include njzoom?
(  381)               DO ji = niw0, niw1
(  382)                  DO jk = 1, jpkm1
(  383)                     DO jj = njw0p1, njw1m1
(  384)                        ij = jj -1 + njmpp
(  385) !byoung                       
(  386)                        swdta(ij,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  387)                        twdta(ij,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  388)                        uwdta(ij,jk,1) = un(ji,jj,jk)*umask(ji,jj,jk)
(  389)                     END DO
(  390)                  END DO
(  391)               END DO
(  392) !byoung
(  393)               DO ji = niw0+1, niw0+4
(  394)                  DO jk = 1, jpkm1
(  395)                     DO jj = njw0p1, njw1m1
(  396)                        ij = jj -1 + nimpp
(  397)                        swdta(ij,jk,ji-niw0+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  398)                        twdta(ij,jk,ji-niw0+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  399)                     END DO
(  400)                  END DO
(  401)               END DO
(  402)            ENDIF
(  403)         ENDIF
(  404) 
(  405)         IF( lp_obc_north)   THEN
(  406)            ! initialisation to zero
(  407) 
(  408) !byoung
(  409)            sndta(:,:,:,:) = 0.e0
(  410)            tndta(:,:,:,:) = 0.e0
(  411)            vndta(:,:,:) = 0.e0
(  412)            !                                    ! ================== !
(  413)            IF( nobc_dta == 0 )   THEN           ! initial state used
(  414)               !                                 ! ================== !
(  415)               !  Fills sndta, tndta, vndta (global arrays)
(  416)               !  Remark: this works for njzoom = 1.
(  417)               !          Should the definition of ij include njzoom?
(  418)               DO jj = njn0, njn1
(  419)                  DO jk = 1, jpkm1
(  420)                     DO ji = nin0p1, nin1m1
(  421)                        ii = ji -1 + nimpp
(  422)                        sndta(ii,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  423)                        tndta(ii,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  424)                        vndta(ii,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
(  425)                     END DO
(  426)                  END DO
(  427)               END DO
(  428) !sujie
(  429)               DO jj = njn0-4, njn0-1
(  430)                  DO jk = 1, jpkm1
(  431)                     DO ji = nin0p1, nin1m1
(  432)                        ii = ji -1 + nimpp
(  433)                        sndta(ii,jk,njn0-jj+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  434)                        tndta(ii,jk,njn0-jj+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  435)                     END DO
(  436)                  END DO
(  437)               END DO






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 41

(  438)            ENDIF
(  439)         ENDIF
(  440) 
(  441)         IF( lp_obc_south )   THEN
(  442)            ! initialisation to zero
(  443) 
(  444)            ssdta(:,:,:,:) = 0.e0
(  445)            tsdta(:,:,:,:) = 0.e0
(  446)            vsdta(:,:,:) = 0.e0
(  447)            !                                    ! ================== !
(  448)            IF( nobc_dta == 0 )   THEN           ! initial state used
(  449)               !                                 ! ================== !
(  450)               !  Fills ssdta, tsdta, vsdta (global arrays)
(  451)               !  Remark: this works for njzoom = 1.
(  452)               !          Should the definition of ij include njzoom?
(  453)               DO jj = njs0, njs1
(  454)                  DO jk = 1, jpkm1
(  455)                     DO ji = nis0p1, nis1m1
(  456)                        ii = ji -1 + nimpp
(  457) 
(  458)                        ssdta(ii,jk,1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  459)                        tsdta(ii,jk,1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  460)                        vsdta(ii,jk,1) = vn(ji,jj,jk)*vmask(ji,jj,jk)
(  461)                     END DO
(  462)                  END DO
(  463)               END DO
(  464) !sujie
(  465)               DO jj = njs0+1, njs0+4
(  466)                  DO jk = 1, jpkm1
(  467)                     DO ji = nis0p1, nis1m1
(  468)                        ii = ji -1 + nimpp
(  469)                        ssdta(ii,jk,jj-njs0+1,1) = sn(ji,jj,jk)*tmask(ji,jj,jk)
(  470)                        tsdta(ii,jk,jj-njs0+1,1) = tn(ji,jj,jk)*tmask(ji,jj,jk)
(  471)                     END DO
(  472)                  END DO
(  473)               END DO
(  474)            ENDIF
(  475)         ENDIF
(  476) 
(  477) 
(  478)      ENDIF        !       end if kt == nit000
(  479) 
(  480)      ! 2.  Initialize the time we are at.
(  481)      !     Does this every time the routine is called,
(  482)      !     excepted when nobc_dta = 0
(  483)      !---------------------------------------------------------------------
(  484)      IF( nobc_dta == 0 )   THEN
(  485)         itimo = 1
(  486)         zxy   = 0
(  487)      ELSE
(  488)         IF( itobc == 1 )   THEN
(  489)            itimo = 1
(  490)         ELSE IF( itobc == 12 )   THEN      !   BC are monthly
(  491)            ! we assume we have climatology in that case
(  492)            iman  = 12
(  493)            i15   = nday / 16
(  494)            imois = nmonth + i15 - 1
(  495)            IF( imois == 0 )   imois = iman






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 42

(  496)            itimo = imois
(  497)         ELSE
(  498)            IF(lwp) WRITE(numout,*) 'data other than constant or monthly',kt
(  499)            iman  = itobc
(  500)            itimo = FLOOR( kt*rdt / (ztcobc(2)-ztcobc(1)) )
(  501)            isrel = kt*rdt
(  502)         ENDIF
(  503)      ENDIF
(  504) 
(  505)      ! 2.1 Read two records in the file if necessary
(  506)      ! ---------------------------------------------
(  507)      IF( ( nobc_dta == 1 ) .AND. ( ( kt == nit000 .AND. nlecto == 0 ) .OR. itimo  /= ntobc1 ) )   THEN
(  508)         nlecto = 1
(  509) 
(  510)         ! Calendar computation
(  511)         IF( itobc == 1 )   THEN            !  BC are constant in time
(  512)            ntobc1 = 1
(  513)            ntobc2 = 1
(  514)         ELSE IF( itobc == 12 )   THEN      !   BC are monthly
(  515)            ntobc1 = itimo         ! first file record used
(  516)            ntobc2 = ntobc1 + 1    ! last  file record used
(  517)            ntobc1 = MOD( ntobc1, iman )
(  518)            IF( ntobc1 == 0 )   ntobc1 = iman
(  519)            ntobc2 = MOD( ntobc2, iman )
(  520)            IF( ntobc2 == 0 )   ntobc2 = iman
(  521)            IF( lwp )   THEN
(  522)               WRITE(numout,*) ' read monthly obc first record file used ntobc1 ', ntobc1
(  523)               WRITE(numout,*) ' read monthly obc last  record file used ntobc2 ', ntobc2
(  524)            ENDIF
(  525)         ELSE
(  526)            isrel=kt*rdt
(  527)            ntobc1 = itimo         ! first file record used
(  528)            ntobc2 = ntobc1 + 1    ! last  file record used
(  529)            ntobc1 = MOD( ntobc1, iman )
(  530)            IF( ntobc1 == 0 )   ntobc1 = iman
(  531)            ntobc2 = MOD( ntobc2, iman )
(  532)            IF( ntobc2 == 0 )   ntobc2 = iman
(  533)            IF(lwp) WRITE(numout,*) ' read obc first record file used ntobc1 ', ntobc1
(  534)            IF(lwp) WRITE(numout,*) ' read obc last  record file used ntobc2 ', ntobc2
(  535)         ENDIF
(  536)                               ! ======================= !
(  537)                               !  BCs read               !
(  538)         !                     ! ======================= !
(  539) 
(  540) 
(  541)         IF( lp_obc_east )   THEN
(  542)            ! ... Read datafile and set temperature, salinity and normal velocity
(  543)            ! ... initialise the sedta, tedta, uedta arrays
(  544) !!DB/CN: Replace IOIPSL calls with lib_ncdf
(  545)            CALL obc_dta_gv ('y','vosaline',jpjef-jpjed+1,ntobc1,pdta_4D=sedta(jpjed:jpjef,:,1:5,1),fn='obceast_TS.nc')
(  546)            CALL obc_dta_gv ('y','vosaline',jpjef-jpjed+1,ntobc2,pdta_4D=sedta(jpjed:jpjef,:,1:5,2),fn='obceast_TS.nc')
(  547)            CALL obc_dta_gv ('y','votemper',jpjef-jpjed+1,ntobc1,pdta_4D=tedta(jpjed:jpjef,:,1:5,1),fn='obceast_TS.nc')
(  548)            CALL obc_dta_gv ('y','votemper',jpjef-jpjed+1,ntobc2,pdta_4D=tedta(jpjed:jpjef,:,1:5,2),fn='obceast_TS.nc')
(  549)            CALL obc_dta_gv ('y','vozocrtx',jpjef-jpjed+1,ntobc1,pdta_3D=uedta(jpjed:jpjef,:,1),fn='obceast_U.nc')
(  550)            CALL obc_dta_gv ('y','vozocrtx',jpjef-jpjed+1,ntobc2,pdta_3D=uedta(jpjed:jpjef,:,2),fn='obceast_U.nc')
(  551) 
(  552)            !  Usually printout is done only once at kt = nit000,
(  553)            !  unless nprint (namelist) > 1 !!DB -- code deleted, see older version






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 43

(  554) 
(  555)         ENDIF
(  556) 
(  557)         IF( lp_obc_west )   THEN
(  558)            ! ... Read datafile and set temperature, salinity and normal velocity
(  559)            ! ... initialise the swdta, twdta, uwdta arrays
(  560) !!DB/CN: Replace IOIPSL calls with lib_ncdf
(  561)            CALL obc_dta_gv ('y','vosaline',jpjwf-jpjwd+1,ntobc1,pdta_4D=swdta(jpjwd:jpjwf,:,1:5,1),fn='obcwest_TS.nc')
(  562)            CALL obc_dta_gv ('y','vosaline',jpjwf-jpjwd+1,ntobc2,pdta_4D=swdta(jpjwd:jpjwf,:,1:5,2),fn='obcwest_TS.nc')
(  563)            CALL obc_dta_gv ('y','votemper',jpjwf-jpjwd+1,ntobc1,pdta_4D=twdta(jpjwd:jpjwf,:,1:5,1),fn='obcwest_TS.nc')
(  564)            CALL obc_dta_gv ('y','votemper',jpjwf-jpjwd+1,ntobc2,pdta_4D=twdta(jpjwd:jpjwf,:,1:5,2),fn='obcwest_TS.nc')
(  565)            CALL obc_dta_gv ('y','vozocrtx',jpjwf-jpjwd+1,ntobc1,pdta_3D=uwdta(jpjwd:jpjwf,:,1),fn='obcwest_U.nc')
(  566)            CALL obc_dta_gv ('y','vozocrtx',jpjwf-jpjwd+1,ntobc2,pdta_3D=uwdta(jpjwd:jpjwf,:,2),fn='obcwest_U.nc')
(  567) !!DB: printout code deleted
(  568)         ENDIF
(  569) 
(  570)         IF( lp_obc_north )   THEN
(  571) !!DB/CN: Replace IOIPSL calls with lib_ncdf
(  572)            CALL obc_dta_gv ('x','vosaline',jpinf-jpind+1,ntobc1,pdta_4D=sndta(jpind:jpinf,:,1:5,1),fn='obcnorth_TS.nc')
(  573)            CALL obc_dta_gv ('x','vosaline',jpinf-jpind+1,ntobc2,pdta_4D=sndta(jpind:jpinf,:,1:5,2),fn='obcnorth_TS.nc')
(  574)            CALL obc_dta_gv ('x','votemper',jpinf-jpind+1,ntobc1,pdta_4D=tndta(jpind:jpinf,:,1:5,1),fn='obcnorth_TS.nc')
(  575)            CALL obc_dta_gv ('x','votemper',jpinf-jpind+1,ntobc2,pdta_4D=tndta(jpind:jpinf,:,1:5,2),fn='obcnorth_TS.nc')
(  576)            CALL obc_dta_gv ('x','vomecrty',jpinf-jpind+1,ntobc1,pdta_3D=vndta(jpind:jpinf,:,1),fn='obcnorth_V.nc')
(  577)            CALL obc_dta_gv ('x','vomecrty',jpinf-jpind+1,ntobc2,pdta_3D=vndta(jpind:jpinf,:,2),fn='obcnorth_V.nc')
(  578)         ENDIF
(  579) 
(  580)         IF( lp_obc_south )   THEN
(  581) !!DB/CN: Replace IOIPSL calls with lib_ncdf
(  582)            CALL obc_dta_gv ('x','vosaline',jpisf-jpisd+1,ntobc1,pdta_4D=ssdta(jpisd:jpisf,:,1:5,1),fn='obcsouth_TS.nc')
(  583)            CALL obc_dta_gv ('x','vosaline',jpisf-jpisd+1,ntobc2,pdta_4D=ssdta(jpisd:jpisf,:,1:5,2),fn='obcsouth_TS.nc')
(  584)            CALL obc_dta_gv ('x','votemper',jpisf-jpisd+1,ntobc1,pdta_4D=tsdta(jpisd:jpisf,:,1:5,1),fn='obcsouth_TS.nc')
(  585)            CALL obc_dta_gv ('x','votemper',jpisf-jpisd+1,ntobc2,pdta_4D=tsdta(jpisd:jpisf,:,1:5,2),fn='obcsouth_TS.nc')
(  586)            CALL obc_dta_gv ('x','vomecrty',jpisf-jpisd+1,ntobc1,pdta_3D=vsdta(jpisd:jpisf,:,1),fn='obcsouth_V.nc')
(  587)            CALL obc_dta_gv ('x','vomecrty',jpisf-jpisd+1,ntobc2,pdta_3D=vsdta(jpisd:jpisf,:,2),fn='obcsouth_V.nc')
(  588)         ENDIF
(  589) 
(  590)      ELSE
(  591) 
(  592)         nlecto = 0        !      no reading of OBC barotropic data
(  593) 
(  594)      ENDIF                !      end of the test on the condition to read or not the files
(  595) 
(  596)      ! 3.  Call at every time step :
(  597)      !     Linear interpolation of BCs to current time step
(  598)      ! ----------------------------------------------------
(  599) 
(  600)      IF( itobc == 1 .OR. nobc_dta == 0 )   THEN
(  601)         zxy = 0.
(  602)      ELSE IF( itobc == 12 )   THEN
(  603)         zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
(  604)      ELSE
(  605)         zxy = (ztcobc(ntobc1)-FLOAT(isrel))/(ztcobc(ntobc1)-ztcobc(ntobc2))
(  606)      ENDIF
(  607) 
(  608)  
(  609) !!DB: ramp for forcing. NB choice is kt rather than (kt-nit000) so the implicit 
(  610) !!assumption is that ramped forcing is not desirable for restarts.
(  611) !!NB: ramp now computed in step.F90






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 44

(  612) !     ramp=tanh(kt*rdt/(2.0*86400.0))
(  613)  
(  614)      IF( lp_obc_east )   THEN
(  615) 
(  616) 
(  617)         DO jk = 1, jpkm1
(  618)            DO jj = nje0p1, nje1m1
(  619)               ij = jj -1 + njmpp
(  620) !!DB
(  621)               sfoe(jj,jk,:) =  ( zxy * sedta(ij,jk,:,2) + &
(  622)                  &           (1.-zxy) * sedta(ij,jk,:,1) ) * temsk5(jj,jk,:) 
(  623)               tfoe(jj,jk,:) =  ( zxy * tedta(ij,jk,:,2) + &
(  624)                  &           (1.-zxy) * tedta(ij,jk,:,1) ) * temsk5(jj,jk,:)
(  625) 
(  626)            END DO
(  627)         END DO
(  628) !AD   interpolate to global
(  629)         DO jk = 1, jpkm1
(  630)            DO jj = jpjed,jpjef
(  631)               uedta1(jj,jk) =  ( zxy * uedta(jj,jk,2) + &
(  632)                  &           (1.-zxy) * uedta(jj,jk,1) ) !* global mask not * uemsk(jj,jk)
(  633)            END DO
(  634)         END DO
(  635) 
(  636) # 637
(  637) !!DB: adjust additional transport here, as applicable
(  638) !!DBG
(  639)         if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_east: ADJUSTING TRANSPORT, narea= ',narea
(  640)  
(  641) 
(  642) !!SBI Region
(  643) !!DB07.30: reverse part of SBI inflow as the Nfld coastal current is in wrong direction
(  644)         do jk = 1, jpkm1
(  645)            do jj = 222,224
(  646)               uedta1(jj,jk) = -uedta1(jj,jk)
(  647)            enddo
(  648)            uedta1(225,jk) = 0.0
(  649)         enddo
(  650) 
(  651) !!DB 2007.11.26
(  652)         fac_SBI = off1 + a1*cos(vfreq*(nday_year+ph1))
(  653)         if(lwp .AND. kt-nit000 == 0) write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SBI = ', &
(  654)                     nday_year,fac_SBI
(  655) 
(  656)         do jk = 1, jpkm1
(  657)            do jj = 220,230
(  658)               uedta1(jj,jk) = fac_SBI*uedta1(jj,jk)
(  659)            enddo
(  660)         enddo
(  661) 
(  662) !DB  East Inflow: Shelf break region ---------------------------------
(  663)         jeast1 = 102
(  664)         jeast2 = 108
(  665)         fac_SB_east = off2 + a2*cos(vfreq*(nday_year+ph2))
(  666)         if(lwp .AND. kt-nit000 == 0) &
(  667)              & write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SB_east = ', nday_year, fac_SB_east
(  668)         if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0)write(numout2,'(A51,2x,2(i7,1x),2(f8.3,1x))') &
(  669)           & 'DBG: obcdta -- kt,nday_year, fac_SBI, fac_SB_east ', kt,nday_year,fac_SBI,fac_SB_east






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 45

(  670) 
(  671) !!DB: globalize this calc
(  672)         ji = nie0
(  673)         ncells = 0      
(  674)         trans_in = 0.0
(  675)         do jk = 1, jpkm1  
(  676)            do jj = jeast1, jeast2
(  677)               trans_in = trans_in + uedta1(jj,jk)*e2v_e(jj)*e3t(jk)*emaskg2(jj,jk)
(  678)               ncells = ncells + emaskg2(jj,jk)
(  679)            enddo
(  680)         enddo
(  681) 
(  682) !!DB adjust SB_east inflow 
(  683) !!Note that I do not ramp U_EW as we need a correct uedta1 and correct ramping is done elsewhere
(  684)         trans_diff = (-fac_SB_east*1.e6) - trans_in
(  685)         do jk = 1, jpkm1  
(  686)            do jj = jeast1, jeast2
(  687)               U_EW(jj,jk) = (trans_diff)/ (e2v_e(jj) * e3t(jk)) * emaskg2(jj,jk)
(  688)            enddo
(  689)         enddo
(  690) 
(  691) 	if (ncells .gt. 0) then
(  692) !AD/DB
(  693)            trans_check = 0.0
(  694)            do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
(  695)               do jj = jeast1, jeast2
(  696)                  uedta1(jj,jk) = uedta1(jj,jk) + U_EW(jj,jk)/float(ncells)
(  697)                  trans_check = trans_check + uedta1(jj,jk)*e2v_e(jj)*e3t(jk)*emaskg2(jj,jk)
(  698)               enddo
(  699)            enddo
(  700) 
(  701) !!DBG: check above calcs
(  702)            if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0) then
(  703)               write(numout2,'(A77,2x,2(i7,1x),4(f8.3,1x))') &
(  704)           & 'DBG: obcdta east -- kt,nday_year,fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
(  705)           &            fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
(  706)            endif
(  707)            if(DBG) then
(  708)               write(200+narea,'(A77,2x,3(i5,1x),4(f8.3,1x))') &
(  709)           & 'obcdta east -- kt,nday_year,ncells, fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
(  710)           &            ncells, fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
(  711)            endif
(  712) 
(  713)         end if
(  714)         
(  715) !DB0731:
(  716) !!modify the input at that Nfld bay that looks like a false open boundary
(  717) !!NB: for Barotropic tide cancel these mean transports
(  718) !AD
(  719) 
(  720)         do jk = 1, jpkm1  
(  721)            do jj = 140,143
(  722)               uedta1(jj,jk) = 0.0
(  723)            enddo
(  724)            do jj = 144, 147
(  725) !!DB 2009.06.26 -- reduce this transport to see if/how it affect transport into SW GSL
(  726)               uedta1(jj,jk) = -0.10 * emaskg2(jj,jk)
(  727) !              uedta1(jj,jk) = -0.02 * emaskg2(jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 46

(  728)            enddo
(  729)         enddo
(  730) 
(  731) 
(  732) # 733
(  733) !AD   Map global to local
(  734)         DO jk = 1, jpkm1
(  735)            DO jj = nje0p1, nje1m1
(  736)               ij = jj -1 + njmpp
(  737)               ufoe(jj,jk) =  uedta1(ij,jk)* uemsk(jj,jk)
(  738)            END DO
(  739)         END DO
(  740) 
(  741) !AD >>>: recompute elevation
(  742)         eta_e(jpjed)=0. !zero eta at SW corner
(  743) 	DO ij=jpjed+1,jpjef
(  744)       	   eta_e(ij) = eta_e(ij-1) - e2v_e(ij-1)*ff_e(ij)/grav * &
(  745)            &    uedta1(ij-1,1)
(  746) 	END DO
(  747)  !AD <<<
(  748) !!DBG
(  749)         if(DBG .AND. kt == nit000) then
(  750)            write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
(  751)                 'DBG: kt,narea,max,min ufoe,uedta1,eta_e ', &
(  752)                 kt,narea,MAXVAL(ufoe(:,:)),MINVAL(ufoe(:,:)), &
(  753)                 MAXVAL(uedta1(:,:)),MINVAL(uedta1(:,:)), &
(  754)                 MAXVAL(eta_e(:)),MINVAL(eta_e(:))
(  755)            DO jk = 1, jpkm1
(  756)               DO jj = nje0p1, nje1m1
(  757)                  ij = jj -1 + njmpp
(  758)                  write(3000+narea,'(3(i3,1x),3(e12.6,1x))') &
(  759)                       jk,jj,ij,ufoe(jj,jk),uedta1(ij,jk),uemsk(jj,jk)
(  760)               END DO
(  761)            END DO
(  762)            do ij = jpjed,jpjef
(  763)               write(8000+narea,'(i3,2x,2(e12.6,2x))')ij,eta_e(ij),uedta1(ij,1)
(  764)            enddo
(  765) 
(  766)         endif
(  767) 
(  768) 
(  769)         ufoe(:,:) = ramp*ufoe(:,:)
(  770) 
(  771)      ENDIF
(  772) 
(  773) 
(  774)      IF( lp_obc_west )   THEN
(  775) 
(  776)         DO jk = 1, jpkm1
(  777)            DO jj = njw0p1, njw1m1
(  778)               ij = jj -1 + njmpp
(  779) !!DB
(  780)               sfow(jj,jk,:) =  ( zxy * swdta(ij,jk,:,2) + &
(  781)                  &           (1.-zxy) * swdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
(  782)               tfow(jj,jk,:) =  ( zxy * twdta(ij,jk,:,2) + &
(  783)                  &           (1.-zxy) * twdta(ij,jk,:,1) ) * twmsk5(jj,jk,:)
(  784)            END DO
(  785)         END DO






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 47

(  786) !AD   interpolate to global
(  787)         DO jk = 1, jpkm1
(  788)            DO jj = jpjwd,jpjwf
(  789)               uwdta1(jj,jk) =  ( zxy * uwdta(jj,jk,2) + &
(  790)                  &           (1.-zxy) * uwdta(jj,jk,1) ) !* global mask not * uemsk(jj,jk)
(  791)            END DO
(  792)         END DO
(  793) 
(  794) 
(  795) # 796
(  796) !!DB: add additional transport here, as applicable
(  797) !!DBG
(  798)         if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_west ADJUSTING TRANSPORT, narea= ',narea
(  799) 
(  800) !DB  west Inflow---------------------------------
(  801) !!DBG
(  802)         if(lwp .AND. kt-nit000 == 0) &
(  803)              &  write(numout2,*)'DBG: (obcdta, dt=1) nday_year, fac_SB_west = fac_SB_east = ',nday_year, fac_SB_east
(  804) 
(  805)         jwest1 = 20
(  806)         jwest2 = 35
(  807)         ji = niw0
(  808)         ncells = 0     
(  809)         trans_in = 0.0
(  810)         do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
(  811)            do jj = jwest1,jwest2
(  812)               trans_in = trans_in + uwdta1(jj,jk)*e2v_w(jj)*e3t(jk)*wmaskg2(jj,jk)
(  813)               ncells = ncells + wmaskg2(jj,jk)
(  814)            enddo
(  815)         enddo
(  816) 
(  817) !!DB adjust SB_east inflow (NB:  -fac_SB_east*1.e6 == inflow in Sv)
(  818) !!Note that I do not ramp U_EW as we need a correct uedta1 and correct ramping is done elsewhere
(  819)         trans_diff = (-fac_SB_east*1.e6) - trans_in
(  820)         do jk = 1, jpkm1  
(  821)            do jj = jwest1, jwest2
(  822)               U_EW(jj,jk) = (trans_diff)/ (e2v_w(jj) * e3t(jk)) * wmaskg2(jj,jk)
(  823)            enddo
(  824)         enddo
(  825) 
(  826) 	
(  827) 	if (ncells .gt. 0) then
(  828) !AD/DB
(  829)            trans_check = 0.0
(  830)            do jk = 1, jpkm1   !!jpkm1 could be changed to some other level
(  831)               do jj = jwest1, jwest2
(  832)                  uwdta1(jj,jk) = uwdta1(jj,jk) + U_EW(jj,jk)/float(ncells)
(  833)                  trans_check = trans_check + uwdta1(jj,jk)*e2v_w(jj)*e3t(jk)*wmaskg2(jj,jk)
(  834)               enddo
(  835)            enddo
(  836) !!DBG: check above calcs
(  837)            if(lwp .AND. mod(kt-nit000,int(86400/rdt)) == 0) then
(  838)               write(numout2,'(A77,2x,2(i7,1x),4(f8.3,1x))') &
(  839)           & 'DBG: obcdta west -- kt,nday_year, fac_SB_east,trans_in,trans_diff,trans_check ', kt,nday_year, &
(  840)           &            fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
(  841) 
(  842)            endif
(  843)            if(DBG) then






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 48

(  844)               write(200+narea,'(A77,2x,3(i5,1x),4(f8.3,1x))') &
(  845)           & 'obcdta west -- kt,nday_year,ncells, fac_SB_east, trans_in,trans_diff,trans_check ', kt,nday_year, &
(  846)           &            ncells, fac_SB_east, trans_in/1.e6, trans_diff/1.e6, trans_check/1.e6
(  847)            endif
(  848) 
(  849)         end if
(  850) 
(  851) # 852
(  852) !AD   Map global to local
(  853)         DO jk = 1, jpkm1
(  854)            DO jj = njw0p1, njw1m1
(  855)               ij = jj -1 + njmpp
(  856)               ufow(jj,jk) =  uwdta1(ij,jk)* uwmsk(jj,jk)
(  857)            END DO
(  858)         END DO
(  859) 
(  860) !AD >>>: recompute elevation
(  861)         eta_w(jpjwd)=0. !zero eta at SW corner
(  862) 	DO ij=jpjwd+1,jpjwf
(  863)       	   eta_w(ij) = eta_w(ij-1) - e2v_w(ij-1)*ff_w(ij)/grav * &
(  864)            &    uwdta1(ij-1,1)
(  865) 	END DO
(  866)  !AD <<<
(  867) 
(  868) !!DBG
(  869)         if(DBG .AND. kt == nit000) then
(  870)            write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
(  871)                 'DBG: kt,narea,max,min ufow,uwdta1,eta_w ', &
(  872)                 kt,narea,MAXVAL(ufow(:,:)),MINVAL(ufow(:,:)), &
(  873)                 MAXVAL(uwdta1(:,:)),MINVAL(uwdta1(:,:)), &
(  874)                 MAXVAL(eta_w(:)),MINVAL(eta_w(:))
(  875)            DO jk = 1, jpkm1
(  876)               DO jj = njw0p1, njw1m1
(  877)                  ij = jj -1 + njmpp
(  878)                  write(4000+narea,'(3(i3,1x),3(e12.6,1x))') &
(  879)                       jk,jj,ij,ufow(jj,jk),uwdta1(ij,jk),uwmsk(jj,jk)
(  880)               END DO
(  881)            END DO
(  882)            do ij = jpjwd,jpjwf
(  883)               write(8100+narea,'(i3,2x,2(e12.6,2x))')ij,eta_w(ij),uwdta1(ij,1)
(  884)            enddo
(  885) 
(  886)         endif
(  887) 
(  888)         ufow(:,:) = ramp * ufow(:,:)
(  889) 
(  890)      ENDIF
(  891) !------------------------------------
(  892) 
(  893) !!DB: north is closed so do nothing. Note that this could lead to future probs
(  894)      IF( lp_obc_north )   THEN
(  895)         !  fills sfon, tfon, vfon (local to each processor)
(  896)         DO jk = 1, jpkm1
(  897)            DO ji = nin0p1, nin1m1
(  898)               ii = ji -1 + nimpp
(  899) !sujie
(  900)               sfon(ji,jk,:) =  ( zxy * sndta(ii,jk,:,2) + &
(  901)                  &           (1.-zxy) * sndta(ii,jk,:,1) )* tnmsk5(ji,jk,:)






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 49

(  902)               tfon(ji,jk,:) =  ( zxy * tndta(ii,jk,:,2) + &
(  903)                  &           (1.-zxy) * tndta(ii,jk,:,1) )* tnmsk5(ji,jk,:)
(  904)               vfon(ji,jk) =  ( zxy * vndta(ii,jk,2) + &
(  905)                  &           (1.-zxy) * vndta(ii,jk,1) ) * vnmsk(ji,jk)
(  906)            END DO
(  907)         END DO
(  908)      ENDIF
(  909) 
(  910)      IF( lp_obc_south )   THEN
(  911) 
(  912)         DO jk = 1, jpkm1
(  913)           DO ji = nis0p1, nis1m1
(  914)              ii = ji -1 + nimpp
(  915) 
(  916)              sfos(ji,jk,:) = ( zxy * ssdta(ii,jk,:,2) + &
(  917)                 &          (1.-zxy) * ssdta(ii,jk,:,1) )* tsmsk5(ji,jk,:)
(  918)              tfos(ji,jk,:) = ( zxy * tsdta(ii,jk,:,2) + &
(  919)                 &          (1.-zxy) * tsdta(ii,jk,:,1) )* tsmsk5(ji,jk,:)
(  920)    
(  921)             END DO
(  922)         END DO
(  923) 
(  924) !AD >>> interpolate to global
(  925)         DO jk = 1, jpkm1
(  926)           DO ii = jpisd,jpisf
(  927)              vsdta1(ii,jk) = ( zxy * vsdta(ii,jk,2) + &
(  928)                 &          (1.-zxy) * vsdta(ii,jk,1) ) !* global mask not * vsmsk(ji,jk)       
(  929)   
(  930)             END DO
(  931)         END DO
(  932) !AD <<<
(  933) 
(  934) # 935
(  935) !!DB: Adjust transport here if desired
(  936) !        if(lwp .AND. kt == nit000) write(numout2,*)'IN lp_obc_south ADJUSTING TRANSPORT'
(  937) 
(  938) # 939
(  939) !   Map global to local
(  940)         DO jk = 1, jpkm1
(  941)           DO ji = nis0p1, nis1m1
(  942)              ii = ji -1 + nimpp
(  943)              vfos(ji,jk) = vsdta1(ii,jk) * vsmsk(ji,jk)  
(  944)             END DO
(  945)         END DO
(  946) 
(  947) !  recompute elevation
(  948)         eta_s(jpisd)=0. !zero eta at SW corner
(  949) 	DO ij=jpisd+1,jpisf
(  950)       	   eta_s(ij) = eta_s(ij-1) + e2u_s(ij-1)*ff_s(ij)/grav * &
(  951)            &    vsdta1(ij-1,1)
(  952) 	END DO
(  953) 
(  954) !!DBG
(  955)         if(DBG .AND. kt == nit000) then
(  956)            write(2000+narea,'(A45,2(i3,1x),6(e12.6,1x))') &
(  957)                 'DBG: kt,narea,max,min vfos,vsdta1,eta_s ', &
(  958)                 kt,narea,MAXVAL(vfos(:,:)),MINVAL(vfos(:,:)), &
(  959)                 MAXVAL(vsdta1(:,:)),MINVAL(vsdta1(:,:)), &






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 50

(  960)                 MAXVAL(eta_s(:)),MINVAL(eta_s(:))
(  961)            DO jk = 1, jpkm1
(  962)               DO ji = nis0p1, nis1m1
(  963)                  ii = ji -1 + nimpp
(  964)                  write(5000+narea,'(3(i3,1x),3(e12.6,1x))') &
(  965)                       jk,ji,ii,vfos(ji,jk),vsdta1(ii,jk), vsmsk(ji,jk)  
(  966)               END DO
(  967)            END DO
(  968)            do ij = jpisd,jpisf
(  969)               write(8200+narea,'(i3,2x,2(e12.6,2x))')ij,eta_s(ij),vsdta1(ij,1)
(  970)            enddo
(  971) 
(  972) 
(  973) 
(  974)         endif
(  975) 
(  976) 
(  977)         vfos(:,:) = ramp * vfos(:,:)
(  978) 
(  979)      ENDIF
(  980) 
(  981) !!DB: Add offset to East boundary SSH -- HARDWIRED
(  982)      do ij = jpjed+1, jpjef
(  983)         eta_e(ij) = eta_e(ij) + eta_s(jpisf)
(  984)      enddo
(  985) 
(  986) 
(  987) !!DB: call ZW routine that zeros non-tidal transports along obc's
(  988)      call obc_ctl(kt)
(  989) 
(  990) 
(  991)   END SUBROUTINE obc_dta
(  992) 
(  993) 
(  994) # 1167
( 1167)   !!-----------------------------------------------------------------------------
( 1168)   !!   Default option
( 1169)   !!-----------------------------------------------------------------------------
( 1170)   SUBROUTINE obc_dta_psi ( kt )       ! Empty routine
( 1171)      !! * Arguments
( 1172)      INTEGER,INTENT(in) :: kt
( 1173)      if(lwp) WRITE(numout,*) 'obc_dta_psi: You should not have seen this print! error?', kt
( 1174)   END SUBROUTINE obc_dta_psi
( 1175) 
( 1176) 
( 1177) # 1179
( 1179)   SUBROUTINE obc_dta_bt( kt, kbt )
( 1180)      !!---------------------------------------------------------------------------
( 1181)      !!                      ***  SUBROUTINE obc_dta  ***
( 1182)      !!
( 1183)      !! ** Purpose :   time interpolation of barotropic data for time-splitting scheme
( 1184)      !!                Data at the boundary must be in m2/s
( 1185)      !!
( 1186)      !! History :
( 1187)      !!   9.0  !  05-11 (V. garnier) Original code
( 1188)      !!---------------------------------------------------------------------------
( 1189)      !! * Arguments
( 1190)      INTEGER, INTENT( in ) ::   kt          ! ocean time-step index






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 51

( 1191)      INTEGER, INTENT( in ) ::   kbt         ! barotropic ocean time-step index
( 1192) 
( 1193)      !! * Local declarations
( 1194)      INTEGER ::   ji, jj, jk, ii, ij   ! dummy loop indices
( 1195)      INTEGER ::   fid_e, fid_w, fid_n, fid_s, fid  ! file identifiers
( 1196)      INTEGER ::   itimo, iman, imois, i15
( 1197)      INTEGER ::   ntobcm, ntobcp, itimom, itimop
( 1198)      REAL(wp) ::  zxy
( 1199)      INTEGER ::   isrel, ikt           ! number of seconds since 1/1/1992
( 1200)      INTEGER ::   ikprint              ! frequency for printouts.
( 1201) !ylu/DB
( 1202)      INTEGER ::   ntide, nntide   
( 1203)      REAL(wp), DIMENSION(5) :: tideperiod    ! in hours
( 1204) !!DB
( 1205) !!Local arrays for BT transport, dimensioned as in obc_oce
( 1206)      REAL(wp), DIMENSION(1:jpisf) ::    &   !:
( 1207)      vbtfos0, sshfos0       !: south boundary barotropic transport without the tide
( 1208)      REAL(wp), DIMENSION(1:jpjwf) ::   &  !:
( 1209)      ubtfow0, sshfow0       !: west boundary barotropic transport without the tide
( 1210)      REAL(wp), DIMENSION(1:jpjef) ::   &  !:
( 1211)      ubtfoe0, sshfoe0       !: east boundary barotropic transport without the tide
( 1212)      INTEGER :: status
( 1213) 
( 1214) !!DB: 1 tide= M2 = good for climatology runs; 5 tides ~ operational runs
( 1215) !     nntide =5
( 1216)      nntide =1
( 1217)      tideperiod(1) = 12.4206  !M2
( 1218)      tideperiod(2) = 12.0000  !S2
( 1219)      tideperiod(3) = 25.8193  !O1
( 1220)      tideperiod(4) = 23.9345  !K1
( 1221)      tideperiod(5) = 12.6583  !N2
( 1222) 
( 1223) 
( 1224)      !!---------------------------------------------------------------------------
( 1225) 
( 1226)      ! 1.   First call: check time frames available in files.
( 1227)      ! -------------------------------------------------------
( 1228) 
( 1229)      IF( kt == nit000 ) THEN
( 1230) 
( 1231)         ! 1.1  Barotropic tangential velocities set to zero
( 1232)         ! -------------------------------------------------
( 1233)         IF( lp_obc_east  ) vbtfoe(:) = 0.e0
( 1234)         IF( lp_obc_west  ) vbtfow(:) = 0.e0
( 1235)         IF( lp_obc_south ) ubtfos(:) = 0.e0
( 1236)         IF( lp_obc_north ) ubtfon(:) = 0.e0
( 1237) 
( 1238)         ! 1.2  Sea surface height and normal barotropic velocities set to zero
( 1239)         !                               or initial conditions if nobc_dta == 0
( 1240)         ! --------------------------------------------------------------------
( 1241) 
( 1242)          IF( lp_obc_east ) THEN
( 1243)             ! initialisation to zero
( 1244)             sshedta(:,:) = 0.e0
( 1245)             ubtedta(:,:) = 0.e0
( 1246)             !                                        ! ================== !
( 1247)             IF( nobc_dta == 0 )   THEN               ! initial state used !
( 1248)                !                                     ! ================== !






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 52

( 1249)                !  Fills sedta, tedta, uedta (global arrays)
( 1250)                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
( 1251)                DO ji = nie0, nie1
( 1252)                   DO jj = nje0p1, nje1m1
( 1253)                      ij = jj -1 + njmpp
( 1254)                      sshedta(ij,1) = sshn(ji+1,jj) * tmask(ji+1,jj,1)
( 1255)                   END DO
( 1256)                END DO
( 1257)             ENDIF
( 1258) !ylu
( 1259)            tidesshemag(:,:)=0.0
( 1260)            tidesshepha(:,:)=0.0
( 1261)            tidevbtemag(:,:)=0.0
( 1262)            tidevbtepha(:,:)=0.0
( 1263) !!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
( 1264)            call ncdf_read_global('obceast_tide.nc','tidesshmag',tidesshemag(jpjed:jpjef,1:nntide),-1,status)
( 1265)            call ncdf_read_global('obceast_tide.nc','tidesshpha',tidesshepha(jpjed:jpjef,1:nntide),-1,status)
( 1266)            call ncdf_read_global('obceast_tide.nc','tidevbtmag',tidevbtemag(jpjed:jpjef,1:nntide),-1,status)
( 1267)            call ncdf_read_global('obceast_tide.nc','tidevbtpha',tidevbtepha(jpjed:jpjef,1:nntide),-1,status)
( 1268) 
( 1269)          ENDIF
( 1270) 
( 1271)          IF( lp_obc_west) THEN
( 1272)             ! initialisation to zero
( 1273)             sshwdta(:,:) = 0.e0
( 1274)             ubtwdta(:,:) = 0.e0
( 1275)             !                                        ! ================== !
( 1276)             IF( nobc_dta == 0 )   THEN               ! initial state used !
( 1277)                !                                     ! ================== !
( 1278)                !  Fills swdta, twdta, uwdta (global arrays)
( 1279)                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
( 1280)                DO ji = niw0, niw1
( 1281)                   DO jj = njw0p1, njw1m1
( 1282)                      ij = jj -1 + njmpp
( 1283)                      sshwdta(ij,1) = sshn(ji,jj) * tmask(ji,jj,1)
( 1284)                   END DO
( 1285)                END DO
( 1286)             ENDIF
( 1287) !ylu
( 1288)            tidesshwmag(:,:)=0.0
( 1289)            tidesshwpha(:,:)=0.0
( 1290)            tidevbtwmag(:,:)=0.0
( 1291)            tidevbtwpha(:,:)=0.0
( 1292) !!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
( 1293)            call ncdf_read_global('obcwest_tide.nc','tidesshmag',tidesshwmag(jpjwd:jpjwf,1:nntide),-1,status)
( 1294)            call ncdf_read_global('obcwest_tide.nc','tidesshpha',tidesshwpha(jpjwd:jpjwf,1:nntide),-1,status)
( 1295)            call ncdf_read_global('obcwest_tide.nc','tidevbtmag',tidevbtwmag(jpjwd:jpjwf,1:nntide),-1,status)
( 1296)            call ncdf_read_global('obcwest_tide.nc','tidevbtpha',tidevbtwpha(jpjwd:jpjwf,1:nntide),-1,status)
( 1297) 
( 1298) 
( 1299)            !CALL flioclo (fid_w)
( 1300)          ENDIF
( 1301) 
( 1302)          IF( lp_obc_north) THEN
( 1303)             ! initialisation to zero
( 1304)             sshndta(:,:) = 0.e0
( 1305)             vbtndta(:,:) = 0.e0
( 1306)             !                                        ! ================== !






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 53

( 1307)             IF( nobc_dta == 0 )   THEN               ! initial state used !
( 1308)                !                                     ! ================== !
( 1309)                !  Fills sndta, tndta, vndta (global arrays)
( 1310)                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
( 1311)                DO jj = njn0, njn1
( 1312)                   DO ji = nin0p1, nin1m1
( 1313)                      DO jk = 1, jpkm1
( 1314)                         ii = ji -1 + nimpp
( 1315)                         vbtndta(ii,1) = vbtndta(ii,1) + vndta(ii,jk,1)*e3t(jk)
( 1316)                      END DO
( 1317)                      sshndta(ii,1) = sshn(ii,jj+1) * tmask(ji,jj+1,1)
( 1318)                   END DO
( 1319)                END DO
( 1320)             ENDIF
( 1321) !ylu
( 1322)            tidesshnmag(:,:)=0.0
( 1323)            tidesshnpha(:,:)=0.0
( 1324)            tidevbtnmag(:,:)=0.0
( 1325)            tidevbtnpha(:,:)=0.0
( 1326) !!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
( 1327)            call ncdf_read_global('obcnorth_tide.nc','tidesshmag',tidesshnmag(jpind:jpinf,1:nntide),-1,status)
( 1328)            call ncdf_read_global('obcnorth_tide.nc','tidesshpha',tidesshnpha(jpind:jpinf,1:nntide),-1,status)
( 1329)            call ncdf_read_global('obcnorth_tide.nc','tidevbtmag',tidevbtnmag(jpind:jpinf,1:nntide),-1,status)
( 1330)            call ncdf_read_global('obcnorth_tide.nc','tidevbtpha',tidevbtnpha(jpind:jpinf,1:nntide),-1,status)
( 1331) 
( 1332)          ENDIF
( 1333) 
( 1334)          IF( lp_obc_south) THEN
( 1335)             ! initialisation to zero
( 1336)             sshsdta(:,:) = 0.e0
( 1337)             vbtsdta(:,:) = 0.e0
( 1338)             !                                        ! ================== !
( 1339)             IF( nobc_dta == 0 )   THEN               ! initial state used !
( 1340)                !                                     ! ================== !
( 1341)                !  Fills ssdta, tsdta, vsdta (global arrays)
( 1342)                !  Remark: this works for njzoom = 1. Should the definition of ij include njzoom?
( 1343)                DO jj = njs0, njs1
( 1344)                   DO ji = nis0p1, nis1m1
( 1345)                      DO jk = 1, jpkm1
( 1346)                         ii = ji -1 + nimpp
( 1347)                         vbtsdta(ii,1) = vbtsdta(ii,1) + vsdta(ii,jk,1)*e3t(jk)
( 1348)                      END DO
( 1349)                      sshsdta(ii,1) = sshn(ji,jj) * tmask(ii,jj,1)
( 1350)                   END DO
( 1351)                END DO
( 1352)             ENDIF
( 1353) !ylu
( 1354)            tidesshsmag(:,:)=0.0
( 1355)            tidesshspha(:,:)=0.0
( 1356)            tidevbtsmag(:,:)=0.0
( 1357)            tidevbtspha(:,:)=0.0
( 1358) !!DB: direct call to read_global(); NB time index should be irrelevant as variables have no time axis
( 1359)            call ncdf_read_global('obcsouth_tide.nc','tidesshmag',tidesshsmag(jpisd:jpisf,1:nntide),-1,status)
( 1360)            call ncdf_read_global('obcsouth_tide.nc','tidesshpha',tidesshspha(jpisd:jpisf,1:nntide),-1,status)
( 1361)            call ncdf_read_global('obcsouth_tide.nc','tidevbtmag',tidevbtsmag(jpisd:jpisf,1:nntide),-1,status)
( 1362)            call ncdf_read_global('obcsouth_tide.nc','tidevbtpha',tidevbtspha(jpisd:jpisf,1:nntide),-1,status)
( 1363) 
( 1364)            !CALL flioclo (fid_s)






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 54

( 1365)          ENDIF
( 1366) 
( 1367)       ENDIF        !       END IF kt == nit000
( 1368) 
( 1369)      
( 1370) !!------------------------------------------------------------------------------------
( 1371)      ! 2.      Initialize the time we are at. Does this every time the routine is called,
( 1372)      !         excepted when nobc_dta = 0
( 1373)      !
( 1374)      IF( nobc_dta == 0) THEN
( 1375)         itimo = 1
( 1376)         zxy   = 0
( 1377)      ELSE
( 1378)         IF(itobc == 1) THEN
( 1379)            itimo = 1
( 1380)         ELSE IF (itobc == 12) THEN      !   BC are monthly
( 1381)            ! we assume we have climatology in that case
( 1382)            iman  = 12
( 1383)            i15   = nday / 16
( 1384)            imois = nmonth + i15 - 1
( 1385)            IF( imois == 0 )   imois = iman
( 1386)            itimo = imois
( 1387)         ELSE
( 1388)            IF(lwp) WRITE(numout,*) 'data other than constant or monthly',kt
( 1389)            iman  = itobc
( 1390)            itimo = FLOOR( kt*rdt / ztcobc(1))
( 1391)            isrel=kt*rdt
( 1392)         ENDIF
( 1393)      ENDIF
( 1394) 
( 1395)      ! 2. Read two records in the file if necessary
( 1396)      ! ---------------------------------------------
( 1397) 
( 1398)      IF( nobc_dta == 1 .AND. nlecto == 1 ) THEN
( 1399)  !sujie add  --------------
( 1400)       ! Calendar computation
( 1401)          IF( itobc == 1 )   THEN            !  BC are constant in time
( 1402)             ntobc1 = 1
( 1403)             ntobc2 = 1
( 1404)             ntobc3 = 1
( 1405)          ELSE IF( itobc == 12 )   THEN      !   BC are monthly
( 1406)             ntobc1 = itimo         ! first file record used
( 1407)             ntobc2 = ntobc1 + 1    ! second  file record used
( 1408)             ntobc3 = ntobc1 + 2    ! last  file record used
( 1409)             ntobc1 = MOD( ntobc1, iman )
( 1410)             IF( ntobc1 == 0 )   ntobc1 = iman
( 1411)             ntobc2 = MOD( ntobc2, iman )
( 1412)             IF( ntobc2 == 0 )   ntobc2 = iman
( 1413)             ntobc3 = MOD( ntobc3, iman )
( 1414)             IF( ntobc3 == 0 )   ntobc3 = iman
( 1415) 
( 1416)          ELSE
( 1417)             isrel=kt*rdt
( 1418)             ntobc1 = itimo         ! first file record used
( 1419)             ntobc2 = ntobc1 + 1    ! second  file record used
( 1420)             ntobc3 = ntobc1 + 2    ! last  file record used
( 1421)             ntobc1 = MOD( ntobc1, iman )
( 1422)             IF( ntobc1 == 0 )   ntobc1 = iman






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 55

( 1423)             ntobc2 = MOD( ntobc2, iman )
( 1424)             IF( ntobc2 == 0 )   ntobc2 = iman
( 1425)             ntobc3 = MOD( ntobc3, iman )
( 1426)             IF( ntobc3 == 0 )   ntobc3 = iman
( 1427)          ENDIF
( 1428) ! ------------------
( 1429) 
( 1430)         IF( lp_obc_east ) THEN
( 1431)            ! ... Read datafile and set sea surface height and barotropic velocity
( 1432)            ! ... initialise the sshedta, ubtedta arrays
( 1433)            sshedta(:,0) = sshedta(:,1)
( 1434)            ubtedta(:,0) = ubtedta(:,1)
( 1435) !!DB/CN: Replace IOIPSL with lib_ncdf
( 1436)            CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc1,pdta_2D=sshedta(jpjed:jpjef,1),fn='obceast_TS.nc')
( 1437)            CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc2,pdta_2D=sshedta(jpjed:jpjef,2),fn='obceast_TS.nc')
( 1438)            IF( lk_dynspg_ts ) THEN
( 1439)               CALL obc_dta_gv ('y','vossurfh',jpjef-jpjed+1,ntobc3,pdta_2D=sshedta(jpjed:jpjef,3),fn='obceast_TS.nc')
( 1440)            ENDIF
( 1441)            CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc1,pdta_2D=ubtedta(jpjed:jpjef,1),fn='obceast_U.nc')
( 1442)            CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc2,pdta_2D=ubtedta(jpjed:jpjef,2),fn='obceast_U.nc')
( 1443)            IF( lk_dynspg_ts ) THEN
( 1444)               CALL obc_dta_gv ('y','vozoubt',jpjef-jpjed+1,ntobc3,pdta_2D=ubtedta(jpjed:jpjef,3),fn='obceast_U.nc')
( 1445)            ENDIF
( 1446) 
( 1447)         ENDIF
( 1448) 
( 1449)         IF( lp_obc_west ) THEN
( 1450)            ! ... Read datafile and set temperature, salinity and normal velocity
( 1451)            ! ... initialise the swdta, twdta, uwdta arrays
( 1452)            sshwdta(:,0) = sshwdta(:,1)
( 1453)            ubtwdta(:,0) = ubtwdta(:,1)
( 1454) !!DB/CN: Replace IOIPSL with lib_ncdf
( 1455)            CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc1,pdta_2D=sshwdta(jpjwd:jpjwf,1),fn='obcwest_TS.nc')
( 1456)            CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc2,pdta_2D=sshwdta(jpjwd:jpjwf,2),fn='obcwest_TS.nc')
( 1457)            IF( lk_dynspg_ts ) THEN
( 1458)               CALL obc_dta_gv ('y','vossurfh',jpjwf-jpjwd+1,ntobc3,pdta_2D=sshwdta(jpjwd:jpjwf,3),fn='obcwest_TS.nc')
( 1459)            ENDIF
( 1460)            CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc1,pdta_2D=ubtwdta(jpjwd:jpjwf,1),fn='obcwest_U.nc')
( 1461)            CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc2,pdta_2D=ubtwdta(jpjwd:jpjwf,2),fn='obcwest_U.nc')
( 1462)            IF( lk_dynspg_ts ) THEN
( 1463)               CALL obc_dta_gv ('y','vozoubt',jpjwf-jpjwd+1,ntobc3,pdta_2D=ubtwdta(jpjwd:jpjwf,3),fn='obcwest_U.nc')
( 1464)            ENDIF
( 1465) 
( 1466)         ENDIF
( 1467) 
( 1468)         IF( lp_obc_north) THEN
( 1469)            ! ... Read datafile and set sea surface height and barotropic velocity
( 1470)            ! ... initialise the sshndta, ubtndta arrays
( 1471)            sshndta(:,0) = sshndta(:,1)
( 1472)            vbtndta(:,0) = vbtndta(:,1)
( 1473) !!DB/CN: Replace IOIPSL with lib_ncdf
( 1474)            CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc1,pdta_2D=sshndta(jpind:jpinf,1),fn='obcnorth_TS.nc')
( 1475)            CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc2,pdta_2D=sshndta(jpind:jpinf,2),fn='obcnorth_TS.nc')
( 1476)            IF( lk_dynspg_ts ) THEN
( 1477)                CALL obc_dta_gv ('x','vossurfh',jpinf-jpind+1,ntobc3,pdta_2D=sshndta(jpind:jpinf,3),fn='obcnorth_TS.nc')
( 1478)            ENDIF
( 1479)            CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc1,pdta_2D=vbtndta(jpind:jpinf,1),fn='obcnorth_V.nc')
( 1480)            CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc2,pdta_2D=vbtndta(jpind:jpinf,2),fn='obcnorth_V.nc')






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 56

( 1481)            IF( lk_dynspg_ts ) THEN
( 1482)               CALL obc_dta_gv ('x','vomevbt',jpinf-jpind+1,ntobc3,pdta_2D=vbtndta(jpind:jpinf,3),fn='obcnorth_V.nc')
( 1483)            ENDIF
( 1484) 
( 1485)         ENDIF
( 1486) 
( 1487)         IF( lp_obc_south) THEN
( 1488)            ! ... Read datafile and set sea surface height and barotropic velocity
( 1489)            ! ... initialise the sshsdta, ubtsdta arrays
( 1490)            sshsdta(:,0) = sshsdta(:,1)
( 1491)            vbtsdta(:,0) = vbtsdta(:,1)
( 1492) !!DB/CN: Replace IOIPSL with lib_ncdf
( 1493)            CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc1,pdta_2D=sshsdta(jpisd:jpisf,1),fn='obcsouth_TS.nc')
( 1494)            CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc2,pdta_2D=sshsdta(jpisd:jpisf,2),fn='obcsouth_TS.nc')
( 1495)            IF( lk_dynspg_ts ) THEN
( 1496)                CALL obc_dta_gv ('x','vossurfh',jpisf-jpisd+1,ntobc3,pdta_2D=sshsdta(jpisd:jpisf,3),fn='obcsouth_TS.nc')
( 1497)            ENDIF
( 1498)            CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc1,pdta_2D=vbtsdta(jpisd:jpisf,1),fn='obcsouth_V.nc')
( 1499)            CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc2,pdta_2D=vbtsdta(jpisd:jpisf,2),fn='obcsouth_V.nc')
( 1500)            IF( lk_dynspg_ts ) THEN
( 1501)                CALL obc_dta_gv ('x','vomevbt',jpisf-jpisd+1,ntobc3,pdta_2D=vbtsdta(jpisd:jpisf,3),fn='obcsouth_V.nc')
( 1502)            ENDIF
( 1503) 
( 1504)         ENDIF
( 1505) 
( 1506)       ENDIF        !      end of the test on the condition to read or not the files
( 1507) 
( 1508)      ! 3.  Call at every time step : Linear interpolation of BCs to current time step
( 1509)      ! ----------------------------------------------------------------------
( 1510) 
( 1511)       IF( lk_dynspg_ts ) THEN
( 1512)          isrel = (kt-1)*rdt + kbt*rdtbt
( 1513) 
( 1514)          IF( nobc_dta == 1 ) THEN
( 1515)             isrel = (kt-1)*rdt + kbt*rdtbt
( 1516)             itimo  = FLOOR(  kt*rdt    / (ztcobc(2)-ztcobc(1)) )
( 1517)             itimom = FLOOR( (kt-1)*rdt / (ztcobc(2)-ztcobc(1)) )
( 1518)             itimop = FLOOR( (kt+1)*rdt / (ztcobc(2)-ztcobc(1)) )
( 1519)             IF( itimom == itimo .AND. itimop == itimo ) THEN
( 1520)                ntobcm = ntobc1
( 1521)                ntobcp = ntobc2
( 1522) 
( 1523)             ELSEIF ( itimom <= itimo .AND. itimop == itimo ) THEN
( 1524)                IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimo ) THEN
( 1525)                   ntobcm = ntobc1-1
( 1526)                   ntobcp = ntobc2-1
( 1527)                ELSE
( 1528)                   ntobcm = ntobc1
( 1529)                   ntobcp = ntobc2
( 1530)                ENDIF
( 1531) 
( 1532)             ELSEIF ( itimom == itimo .AND. itimop >= itimo ) THEN
( 1533)                IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimop ) THEN
( 1534)                   ntobcm = ntobc1
( 1535)                   ntobcp = ntobc2
( 1536)                ELSE
( 1537)                   ntobcm = ntobc1+1
( 1538)                   ntobcp = ntobc2+1






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 57

( 1539)                ENDIF
( 1540) 
( 1541)             ELSEIF ( itimom == itimo-1 .AND. itimop == itimo+1 ) THEN
( 1542)                IF(  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimo ) THEN
( 1543)                   ntobcm = ntobc1-1
( 1544)                   ntobcp = ntobc2-1
( 1545)                ELSEIF (  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) < itimop ) THEN
( 1546)                   ntobcm = ntobc1
( 1547)                   ntobcp = ntobc2
( 1548)                ELSEIF (  FLOOR( isrel / (ztcobc(2)-ztcobc(1)) ) == itimop ) THEN
( 1549)                   ntobcm = ntobc1+1
( 1550)                   ntobcp = ntobc2+2
( 1551)                ELSE
( 1552)                   IF(lwp) WRITE(numout, *) 'obc_dta_bt: You should not have seen this print! error 1?'
( 1553)                ENDIF
( 1554)             ELSE
( 1555)                IF(lwp) WRITE(numout, *) 'obc_dta_bt: You should not have seen this print! error 2?'
( 1556)             ENDIF
( 1557) 
( 1558)          ENDIF
( 1559) 
( 1560)       ELSE IF( lk_dynspg_exp ) THEN
( 1561)          isrel=kt*rdt
( 1562)          ntobcm = ntobc1
( 1563)          ntobcp = ntobc2
( 1564)       ENDIF
( 1565) 
( 1566)       IF( itobc == 1 .OR. nobc_dta == 0 ) THEN
( 1567)          zxy = 0.e0
( 1568)       ELSE IF( itobc == 12 ) THEN
( 1569)          zxy = FLOAT( nday + 15 - 30 * i15 ) / 30.
( 1570)       ELSE
( 1571)          zxy = (ztcobc(ntobcm)-FLOAT(isrel)) / (ztcobc(ntobcm)-ztcobc(ntobcp))
( 1572)       ENDIF
( 1573) 
( 1574) 
( 1575) !DB: 2007.12.27 -- Note that I now keep the mean components separate
( 1576) !    Formulae should work for any number of tide components
( 1577) !    REM: velocity components already zxy-weighted & multiplied by ramp
( 1578)       IF( lp_obc_east ) THEN    
( 1579) 
( 1580)          do ji = nie0, nie1
( 1581)             do jj = nje0p1, nje1m1
( 1582) !!DB: REM that ufoe has been multiplied by ramp
( 1583)                ubtfoe0(jj) = 0.0
( 1584)                do jk = 1, jpkm1
( 1585)                   ubtfoe0(jj) = ubtfoe0(jj) + ufoe(jj,jk) * e3t(jk)
( 1586)                enddo
( 1587)                ij = jj -1 + njmpp
( 1588) !!DB 
( 1589)                sshfoe0(jj) = eta_e(ij) * temsk(jj,1)
( 1590)                ubtfoe0(jj) = ubtfoe0(jj) * uemsk(jj,1) 
( 1591)                ubtfoe(jj) = 0.0
( 1592)                sshfoe(jj) = 0.0
( 1593) 
( 1594) # 1595
( 1595) !!DB: Keep this ifdef open
( 1596) 






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 58

( 1597) # 1598
( 1598)                do ntide = 1, nntide
( 1599)                   sshfoe(jj) = sshfoe(jj) + tidesshemag(ij,ntide)  &
( 1600)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1601)                        - tidesshepha(ij,ntide) ) * temsk(jj,1)
( 1602) 
( 1603)                   ubtfoe(jj) = ubtfoe(jj) &
( 1604)                        + tidevbtemag(ij,ntide)  &
( 1605)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1606)                        - tidevbtepha(ij,ntide) ) * uemsk(jj,1)
( 1607)                enddo
( 1608) !!DB 
( 1609)                sshfoe(jj) = ramp*(sshfoe0(jj) + sshfoe(jj))
( 1610)                ubtfoe(jj) = ramp*ubtfoe(jj) + ubtfoe0(jj)
( 1611) 
( 1612)             enddo
( 1613)          enddo
( 1614) 
( 1615) !-----------------------------------------------
( 1616) 
( 1617)       ENDIF
( 1618) 
( 1619) !DB: 2007.12.27 -- Note that I now keep the mean components separate
( 1620) !    Formulae should work for any number of tide components
( 1621) !    REM: velocity components already zxy-weighted & multiplied by ramp
( 1622)       IF( lp_obc_west) THEN   
( 1623) 
( 1624)          do ji = niw0, niw1
( 1625)             do jj = njw0p1, njw1m1
( 1626)                ubtfow0(jj) = 0.0 
( 1627)                do jk = 1, jpkm1
( 1628)                   ubtfow0(jj) = ubtfow0(jj) + ufow(jj,jk) * e3t(jk)
( 1629)                enddo
( 1630)                ij = jj -1 + njmpp
( 1631) 
( 1632)                ubtfow0(jj) = ubtfow0(jj) * uwmsk(jj,1) 
( 1633)                ubtfow(jj) = 0.0
( 1634)                sshfow0(jj) = eta_w(ij)*twmsk(jj,1)
( 1635)                sshfow(jj) = 0.0
( 1636) 
( 1637)                do ntide = 1, nntide
( 1638)                   sshfow(jj) =  sshfow(jj) + tidesshwmag(ij,ntide)  &
( 1639)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1640)                        - tidesshwpha(ij,ntide) ) * twmsk(jj,1) 
( 1641)                   
( 1642)                   ubtfow(jj) = ubtfow(jj) &
( 1643)                        + tidevbtwmag(ij,ntide)  &
( 1644)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1645)                        - tidevbtwpha(ij,ntide) ) * uwmsk(jj,1)
( 1646)                enddo
( 1647)                sshfow(jj) =  ramp*(sshfow(jj) + sshfow0(jj))
( 1648)                ubtfow(jj) = ramp*ubtfow(jj)+ ubtfow0(jj) 
( 1649)             enddo
( 1650)          enddo
( 1651) 
( 1652) !------------------------------------------------
( 1653)       ENDIF
( 1654) 






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 59

( 1655) !---------------------------------------------------------------------------
( 1656) !byoung for St. Lawrence River runoff: DELETED -- see old code
( 1657) 
( 1658) !----------------------------------------------------------------------------
( 1659) 
( 1660) !!DB: north should be closed. If not, modify this section
( 1661)       IF( lp_obc_north) THEN           !  fills sshfon, vbtfon (local to each processor)
( 1662)          DO ji = nin0p1, nin1m1
( 1663)             ii = ji -1 + nimpp
( 1664)             sshfon(ji) = ( zxy * sshndta(ii,2) + (1.-zxy) * sshndta(ii,1) ) * tnmsk(ji,1)
( 1665)             vbtfon(ji) = ( zxy * vbtndta(ii,2) + (1.-zxy) * vbtndta(ii,1) ) * vnmsk(ji,1)
( 1666) !ylu
( 1667)           do ntide = 1,4  
( 1668)             sshfon(ji) = sshfon(ji) + ramp*tidesshnmag(ii,ntide)  &
( 1669)                          * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1670)                          - tidesshnpha(ii,ntide) ) * tnmsk(ji,1) 
( 1671)             vbtfon(ji) = vbtfon(ji) + ramp*tidevbtnmag(ii,ntide)  &
( 1672)                          * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1673)                          - tidevbtnpha(ii,ntide) ) * vnmsk(ji,1) 
( 1674)           enddo
( 1675)          END DO
( 1676)       ENDIF
( 1677) 
( 1678) 
( 1679) !DB: 2007.12.27 -- Note that I now keep the mean components separate
( 1680) !    Formulae should work for any number of tide components
( 1681) !    REM: velocity components already zxy-weighted & multiplied by ramp
( 1682)       IF( lp_obc_south) THEN     
( 1683) 
( 1684)          do jj = njs0, njs1
( 1685)             do ji = nis0p1, nis1m1
( 1686)                vbtfos0(ji) = 0.0 
( 1687)                do jk = 1, jpkm1
( 1688)                   vbtfos0(ji) = vbtfos0(ji) + vfos(ji,jk) * e3t(jk)
( 1689)                enddo
( 1690)                ii = ji -1 + nimpp
( 1691) 
( 1692)                vbtfos0(ji) = vbtfos0(ji)*vsmsk(ji,1)
( 1693)                vbtfos(ji) = 0.0
( 1694)                sshfos0(ji) = eta_s(ii)* tsmsk(ji,1)
( 1695)                sshfos(ji) = 0.0
( 1696) 
( 1697)                do ntide = 1, nntide
( 1698)                   sshfos(ji) = sshfos(ji) + tidesshsmag(ii,ntide)  &
( 1699)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1700)                        - tidesshspha(ii,ntide) ) * tsmsk(ji,1)
( 1701) 
( 1702)                   vbtfos(ji) =  vbtfos(ji) + & 
( 1703)                        tidevbtsmag(ii,ntide)  &
( 1704)                        * cos(2.0*3.1416/(tideperiod(ntide)*3600.)*kt*rdt  &
( 1705)                        - tidevbtspha(ii,ntide) ) * vsmsk(ji,1) 
( 1706)                enddo
( 1707)                sshfos(ji) = ramp*(sshfos(ji) + sshfos0(ji))
( 1708)                vbtfos(ji) =  ramp*vbtfos(ji) +  vbtfos0(ji)
( 1709)             enddo
( 1710)          enddo
( 1711) 
( 1712)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 60

( 1713) 
( 1714)   END SUBROUTINE obc_dta_bt
( 1715) 
( 1716) 
( 1717) 
( 1718) # 1730
( 1730) !!DB/CN 2008.11
( 1731)   SUBROUTINE obc_dta_gv (cldim,clobc,kobcij,ktobc,pdta_2D,pdta_3D,pdta_4D,fn)
( 1732)      !!-----------------------------------------------------------------------------
( 1733)      !!                       ***  SUBROUTINE obc_dta_gv  ***
( 1734)      !!
( 1735)      !! ** Purpose :   Read an OBC forcing field from netcdf file
( 1736)      !!                Input file are supposed to be 3D e.g.
( 1737)      !!                - for a South or North OB : longitude x depth x time
( 1738)      !!		- for a West or East OB : latitude x depth x time
( 1739)      !!
( 1740)      !! History :
( 1741)      !!        !  04-06 (A.-M. Treguier, F. Durand) Original code
( 1742)      !!        !  05-02 (J. Bellier, C. Talandier) use fliocom CALL
( 1743)      !!----------------------------------------------------------------------------
( 1744)      !! * Arguments
( 1745)      INTEGER, INTENT(IN) ::   &
( 1746) !        ifid  ,               & ! netcdf file name identifier
( 1747)         kobcij,               & ! Horizontal (i or j) dimension of the array
( 1748)         ktobc                   ! starting time index read
( 1749)      CHARACTER(LEN=*), INTENT(IN)    ::   &
( 1750)         cldim,                & ! dimension along which is the open boundary ('x' or 'y')
( 1751)         clobc                   ! name of the netcdf variable read
( 1752)      REAL, DIMENSION(kobcij,jpk,1), INTENT(OUT), OPTIONAL ::   &
( 1753)         pdta_3D                 ! 3D array of OBC forcing field
( 1754) 
( 1755)      REAL, DIMENSION(kobcij,1), INTENT(OUT), OPTIONAL ::   &
( 1756)         pdta_2D                 ! 2D array of OBC forcing field
( 1757)      REAL, DIMENSION(kobcij,jpk,5,1), INTENT(OUT), OPTIONAL ::   &
( 1758)         pdta_4D                 ! 4D array of OBC forcing field
( 1759) !!CN
( 1760)      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fn ! name of netCDF file
( 1761)      LOGICAL :: fnswitch ! logical indicating whether to use filename & lib_ncdf or not
( 1762)                           ! to be removed once IOIPSL is entirely gone from obcdta
( 1763)      REAL, DIMENSION(kobcij) :: buf_1D             ! Needed for lib_ncdf calls
( 1764)      REAL, DIMENSION(kobcij,jpk) :: buf_2D         ! Needed for lib_ncdf calls
( 1765)      REAL, DIMENSION(kobcij,jpk,5) :: buf_3D       ! Needed for lib_ncdf calls
( 1766)      INTEGER :: f_stat
( 1767)      LOGICAL :: tide
( 1768) 
( 1769)      !! * Local declarations
( 1770)      INTEGER ::   indim
( 1771)      LOGICAL ::   l_exv
( 1772)      INTEGER,DIMENSION(4) ::   f_d, istart, icount
( 1773)      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   v_tmp_4
( 1774)      !----------------------------------------------------------------------
( 1775) 
( 1776)      !!CN: Replacing IOIPSL calls with lib_ncdf
( 1777)      fnswitch = .TRUE.
( 1778)      CALL ncdf_get_num_dims(fn, TRIM(clobc), indim, f_stat);
( 1779) 
( 1780)      l_exv = .TRUE.
( 1781)      IF( l_exv )   THEN






PGF90 (Version     10.5)          09/23/2010  09:49:39      page 61

( 1782)         ! checks the number of dimensions
( 1783)         IF( indim == 2 )   THEN
( 1784)            istart(1:2) = (/ 1     , ktobc /)
( 1785)            icount(1:2) = (/ kobcij, 1     /)
( 1786)            CALL ncdf_read(fn, TRIM(clobc), buf_1D, f_stat)
( 1787)            pdta_2D(:,1) = buf_1D
( 1788)         ELSE IF( indim == 3 )   THEN
( 1789)            istart(1:3) = (/ 1     , 1    , ktobc /)
( 1790)            icount(1:3) = (/ kobcij, jpk  , 1     /)
( 1791)            CALL ncdf_read_global(fn, TRIM(clobc), buf_2D, -ktobc, f_stat)
( 1792)            pdta_3D(:,:,1) = buf_2D
( 1793)         ELSE IF( indim == 4 )   THEN
( 1794)            istart(1:4) = (/ 1, 1, 1, ktobc /)
( 1795)            icount(1:4) = (/ kobcij, jpk  , 5, 1 /)
( 1796)            CALL ncdf_read_global(fn, TRIM(clobc), buf_3D, -ktobc, f_stat)
( 1797)            pdta_4D(:,:,:,1) = buf_3d
( 1798)         ELSE
( 1799)            IF( lwp )   THEN
( 1800)               WRITE(numout,*) ' Problem in OBC file for ',TRIM(clobc),' :'
( 1801)               WRITE(numout,*) ' number of dimensions (not 3 or 4) =',indim
( 1802)            ENDIF
( 1803)            STOP
( 1804)         ENDIF
( 1805)      ELSE
( 1806)         if(lwp) WRITE(numout,*) ' Variable ',TRIM(clobc),' not found'
( 1807)      ENDIF
( 1808) 
( 1809)   END SUBROUTINE obc_dta_gv
( 1810) 
( 1811) 
( 1812) 
( 1813) # 1828
( 1828)   !!=====================================================================
( 1829) END MODULE obcdta
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (obcdta.F90: 603)
