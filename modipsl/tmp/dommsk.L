


PGF90 (Version     10.5)          09/23/2010  09:49:32      page 1

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

Filename: dommsk.F90

(    1) MODULE dommsk
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE dommsk   ***
(    4)    !! Ocean initialization : domain land/sea mask 
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_msk        : compute land/ocean mask
(    9)    !!   dom_msk_nsa    : update land/ocean mask when no-slip accurate
(   10)    !!                    option is used.
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE obc_oce         ! ocean open boundary conditions
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18)    USE lib_mpp
(   19)    USE solisl          ! ???
(   20)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC dom_msk        ! routine called by inidom.F90
(   27) 
(   28)    !! * Module variables
(   29)    REAL(wp) ::   &
(   30)       shlat = 2.   ! type of lateral boundary condition on velocity (namelist namlbc)
(   31)    
(   32)    !! * Substitutions
(   33) #  include "vectopt_loop_substitute.h90"
(   34)    !!---------------------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/dommsk.F90,v 1.10 2006/03/21 08:25:42 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!---------------------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41)    
(   42)    SUBROUTINE dom_msk
(   43)       !!---------------------------------------------------------------------
(   44)       !!                 ***  ROUTINE dom_msk  ***
(   45)       !!
(   46)       !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 2

(   47)       !!      zontal velocity points (u & v), vorticity points (f) and baro-
(   48)       !!      tropic stream function  points (b).
(   49)       !!        Set mbathy to the number of non-zero w-levels of a water column
(   50)       !!      (if island in the domain (lk_isl=T), this is done latter in
(   51)       !!      routine solver_init)
(   52)       !!
(   53)       !! ** Method  :   The ocean/land mask is computed from the basin bathy-
(   54)       !!      metry in level (mbathy) which is defined or read in dommba.
(   55)       !!      mbathy equals 0 over continental T-point, -n over the nth 
(   56)       !!      island T-point, and the number of ocean level over the ocean.
(   57)       !!
(   58)       !!      At a given position (ji,jj,jk) the ocean/land mask is given by:
(   59)       !!      t-point : 0. IF mbathy( ji ,jj) =< 0
(   60)       !!                1. IF mbathy( ji ,jj) >= jk
(   61)       !!      u-point : 0. IF mbathy( ji ,jj)  or mbathy(ji+1, jj ) =< 0
(   62)       !!                1. IF mbathy( ji ,jj) and mbathy(ji+1, jj ) >= jk.
(   63)       !!      v-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1) =< 0
(   64)       !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1) >= jk.
(   65)       !!      f-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1)
(   66)       !!                   or mbathy(ji+1,jj)  or mbathy(ji+1,jj+1) =< 0
(   67)       !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1)
(   68)       !!                and mbathy(ji+1,jj) and mbathy(ji+1,jj+1) >= jk.
(   69)       !!      b-point : the same definition as for f-point of the first ocean
(   70)       !!                level (surface level) but with 0 along coastlines.
(   71)       !!
(   72)       !!        The lateral friction is set through the value of fmask along
(   73)       !!      the coast and topography. This value is defined by shlat, a
(   74)       !!      namelist parameter:
(   75)       !!         shlat = 0, free slip  (no shear along the coast)
(   76)       !!         shlat = 2, no slip  (specified zero velocity at the coast)
(   77)       !!         0 < shlat < 2, partial slip   | non-linear velocity profile
(   78)       !!         2 < shlat, strong slip        | in the lateral boundary layer
(   79)       !!
(   80)       !!      N.B. If nperio not equal to 0, the land/ocean mask arrays
(   81)       !!      are defined with the proper value at lateral domain boundaries,
(   82)       !!      but bmask. indeed, bmask defined the domain over which the
(   83)       !!      barotropic stream function is computed. this domain cannot
(   84)       !!      contain identical columns because the matrix associated with
(   85)       !!      the barotropic stream function equation is then no more inverti-
(   86)       !!      ble. therefore bmask is set to 0 along lateral domain boundaries
(   87)       !!      even IF nperio is not zero.
(   88)       !!
(   89)       !!      In case of open boundaries (lk_obc=T):
(   90)       !!        - tmask is set to 1 on the points to be computed bay the open
(   91)       !!          boundaries routines.
(   92)       !!        - bmask is  set to 0 on the open boundaries.
(   93)       !!
(   94)       !!      Set mbathy to the number of non-zero w-levels of a water column
(   95)       !!                  mbathy = min( mbathy, 1 ) + 1
(   96)       !!      (note that the minimum value of mbathy is 2).
(   97)       !!
(   98)       !! ** Action :
(   99)       !!                     tmask    : land/ocean mask at t-point (=0. or 1.)
(  100)       !!                     umask    : land/ocean mask at u-point (=0. or 1.)
(  101)       !!                     vmask    : land/ocean mask at v-point (=0. or 1.)
(  102)       !!                     fmask    : land/ocean mask at f-point (=0. or 1.)
(  103)       !!                          =shlat along lateral boundaries
(  104)       !!                     bmask    : land/ocean mask at barotropic stream






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 3

(  105)       !!                          function point (=0. or 1.) and set to
(  106)       !!                          0 along lateral boundaries
(  107)       !!                   mbathy   : number of non-zero w-levels 
(  108)       !!
(  109)       !! History :
(  110)       !!        !  87-07  (G. Madec)  Original code
(  111)       !!        !  91-12  (G. Madec)
(  112)       !!        !  92-06  (M. Imbard)
(  113)       !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
(  114)       !!        !  96-01  (G. Madec)  suppression of common work arrays
(  115)       !!        !  96-05  (G. Madec)  mask computed from tmask and sup-
(  116)       !!                 pression of the double computation of bmask
(  117)       !!        !  97-02  (G. Madec)  mesh information put in domhgr.F
(  118)       !!        !  97-07  (G. Madec)  modification of mbathy and fmask
(  119)       !!        !  98-05  (G. Roullet)  free surface
(  120)       !!        !  00-03  (G. Madec)  no slip accurate
(  121)       !!        !  01-09  (J.-M. Molines)  Open boundaries
(  122)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  123)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  124)       !!----------------------------------------------------------------------
(  125)       !! *Local declarations
(  126)       INTEGER  ::   ji, jj, jk, ii     ! dummy loop indices
(  127)       INTEGER  ::   iif, iil, ijf, ijl
(  128)       INTEGER  ::   ii0, ii1, ij0, ij1
(  129)       INTEGER, DIMENSION(jpi,jpj) ::  imsk
(  130) 
(  131)       REAL(wp), DIMENSION(jpi,jpj) ::   zwf
(  132) 
(  133)       NAMELIST/namlbc/ shlat
(  134)       !!---------------------------------------------------------------------
(  135)       
(  136) 
(  137)       ! Namelist namlbc : lateral momentum boundary condition
(  138)       REWIND( numnam )
(  139)       READ  ( numnam, namlbc )
(  140)       IF(lwp) THEN
(  141)          WRITE(numout,*)
(  142)          WRITE(numout,*) 'dommsk : ocean mask '
(  143)          WRITE(numout,*) '~~~~~~'
(  144)          WRITE(numout,*) '         Namelist namlbc'
(  145)          WRITE(numout,*) '            lateral momentum boundary cond. shlat = ',shlat
(  146)       ENDIF
(  147) 
(  148)       IF ( shlat == 0. ) THEN
(  149)           IF(lwp) WRITE(numout,*) '         ocean lateral free-slip '
(  150)         ELSEIF ( shlat  ==  2. ) THEN
(  151)           IF(lwp) WRITE(numout,*) '         ocean lateral  no-slip '
(  152)         ELSEIF ( 0. < shlat .AND. shlat < 2. ) THEN
(  153)           IF(lwp) WRITE(numout,*) '         ocean lateral  partial-slip '
(  154)         ELSEIF ( 2. < shlat ) THEN
(  155)           IF(lwp) WRITE(numout,*) '         ocean lateral  strong-slip '
(  156)         ELSE
(  157)           IF(lwp) WRITE(numout,cform_err)
(  158)           IF(lwp) WRITE(numout,*) ' shlat is negative = ', shlat
(  159)           nstop = nstop + 1
(  160)       ENDIF
(  161) 
(  162)       ! 1. Ocean/land mask at t-point (computed from mbathy)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 4

(  163)       ! -----------------------------
(  164)       ! Tmask has already the right boundary conditions since mbathy is ok
(  165) 
(  166)       tmask(:,:,:) = 0.e0
(  167)       DO jk = 1, jpk
(  168)          DO jj = 1, jpj
(  169)             DO ji = 1, jpi
(  170)                IF( FLOAT( mbathy(ji,jj)-jk )+.1 >= 0.e0 ) tmask(ji,jj,jk) = 1.e0
(  171)             END DO  
(  172)          END DO  
(  173)       END DO  
(  174) 
(  175) !ylu
(  176) !!DB -- 2008.11.24 -- Is this the cause of my probs in SE/SW corners ???
(  177) !SPINUP Run ---> garbage ~10d
(  178) !      write(2200+narea,*)'DBG: not setting tmask in SE/SW corners to zero'
(  179)       tmask(mi0(1):mi1(2),mj0(1):mj1(2),1:jpk)=0.0
(  180)       tmask(mi0(196):mi1(197),mj0(1):mj1(2),1:jpk)=0.0
(  181) 
(  182) 
(  183) #if defined key_zdfkpp
(  184) !!DB: delete ORCA
(  185) !      IF( cp_cfg == 'orca' )   THEN
(  186) #endif
(  187) 
(  188)       ! Interior domain mask (used for global sum)
(  189)       ! --------------------
(  190) 
(  191)       tmask_i(:,:) = tmask(:,:,1)
(  192)       iif = jpreci                         ! ???
(  193)       iil = nlci - jpreci + 1
(  194)       ijf = jprecj                         ! ???
(  195)       ijl = nlcj - jprecj + 1
(  196) 
(  197)       tmask_i( 1 :iif,   :   ) = 0.e0      ! first columns
(  198)       tmask_i(iil:jpi,   :   ) = 0.e0      ! last  columns (including mpp extra columns)
(  199)       tmask_i(   :   , 1 :ijf) = 0.e0      ! first rows
(  200)       tmask_i(   :   ,ijl:jpj) = 0.e0      ! last  rows (including mpp extra rows)
(  201) 
(  202)       ! north fold mask
(  203)       tpol(1:jpiglo) = 1.e0 
(  204)       fpol(1:jpiglo) = 1.e0
(  205)       IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
(  206)          tpol(jpiglo/2+1:jpiglo) = 0.e0
(  207)          fpol(     1    :jpiglo) = 0.e0
(  208)          ! T-point pivot: only half of the nlcj-1 row
(  209)          IF( mjg(nlej) == jpjglo )   THEN
(  210)             DO ji = iif+1, iil-1
(  211)                tmask_i(ji,nlej-1) = tmask_i(ji,nlej-1) * tpol(mig(ji))
(  212)             END DO
(  213)          ENDIF
(  214)       ENDIF
(  215)       IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
(  216)          tpol(     1    :jpiglo) = 0.e0
(  217)          fpol(jpiglo/2+1:jpiglo) = 0.e0
(  218)       ENDIF
(  219) 
(  220)       ! 2. Ocean/land mask at u-,  v-, and z-points (computed from tmask)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 5

(  221)       ! -------------------------------------------
(  222)       
(  223)       ! Computation
(  224)       DO jk = 1, jpk
(  225)          DO jj = 1, jpjm1
(  226)             DO ji = 1, fs_jpim1   ! vector loop
(  227)                umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
(  228)                vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
(  229)                fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
(  230)                   &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
(  231)             END DO
(  232)          END DO
(  233)       END DO
(  234) 
(  235) !!DB: delete ORCA
(  236) !      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN
(  237) 
(  238)       ! Lateral boundary conditions
(  239)       CALL lbc_lnk( umask, 'U', 1. )
(  240)       CALL lbc_lnk( vmask, 'V', 1. )
(  241)       CALL lbc_lnk( fmask, 'F', 1. )
(  242) 
(  243)       nmask(:,:,:)=fmask(:,:,:)
(  244)       ! 4. ocean/land mask for the elliptic equation
(  245)       ! --------------------------------------------
(  246)       
(  247)       ! Computation
(  248)       IF( lk_dynspg_rl ) THEN
(  249)          bmask(:,:) = fmask(:,:,1)       ! elliptic equation is written at f-point
(  250)       ELSE
(  251)          bmask(:,:) = tmask(:,:,1)       ! elliptic equation is written at t-point
(  252)       ENDIF
(  253)       
(  254)       ! Boundary conditions
(  255)       !   cyclic east-west : bmask must be set to 0. on rows 1 and jpi
(  256)       IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
(  257)          bmask( 1 ,:) = 0.e0
(  258)          bmask(jpi,:) = 0.e0
(  259)       ENDIF
(  260)       
(  261)       !   south symmetric :  bmask must be set to 0. on row 1
(  262)       IF( nperio == 2 ) THEN
(  263)          bmask(:, 1 ) = 0.e0
(  264)       ENDIF
(  265)       
(  266)       !   north fold : 
(  267)       IF( nperio == 3 .OR. nperio == 4 ) THEN
(  268)          IF( lk_dynspg_rl ) THEN
(  269)             ! T-pt pivot and F-pt elliptic eq. : bmask set to 0. on rows jpj-1 and jpj
(  270)             bmask(:,jpj-1) = 0.e0
(  271)             bmask(:,jpj  ) = 0.e0
(  272)          ELSE
(  273)             ! T-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj and on half jpjglo-1 row
(  274)             DO ji = 1, jpi
(  275)                ii = ji + nimpp - 1
(  276)                bmask(ji,jpj-1) = bmask(ji,jpj-1) * tpol(ii)
(  277)                bmask(ji,jpj  ) = 0.e0
(  278)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 6

(  279)          ENDIF
(  280)       ENDIF
(  281)       IF( nperio == 5 .OR. nperio == 6 ) THEN
(  282)          IF( lk_dynspg_rl ) THEN
(  283)             ! F-pt pivot and F-pt elliptic eq. : bmask set to 0. on row jpj and on half jpjglo-1 row
(  284)             DO ji = 1, jpi
(  285)                ii = ji + nimpp - 1
(  286)                bmask(ji,jpj-1) = bmask(ji,jpj-1) * fpol(ii)
(  287)                bmask(ji,jpj  ) = 0.e0
(  288)             END DO
(  289)          ELSE
(  290)             ! F-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj
(  291)             bmask(:,jpj) = 0.e0
(  292)          ENDIF
(  293)       ENDIF
(  294) 
(  295)       ! Mpp boundary conditions: bmask is set to zero on the overlap
(  296)       ! region for all elliptic solvers
(  297) 
(  298)       IF( lk_mpp ) THEN
(  299)          IF( nbondi /= -1 .AND. nbondi /= 2 )   bmask(  1 :jpreci,:) = 0.e0
(  300)          IF( nbondi /=  1 .AND. nbondi /= 2 )   bmask(nlci:jpi   ,:) = 0.e0
(  301)          IF( nbondj /= -1 .AND. nbondj /= 2 )   bmask(:,  1 :jprecj) = 0.e0
(  302)          IF( nbondj /=  1 .AND. nbondj /= 2 )   bmask(:,nlcj:jpj   ) = 0.e0
(  303)       
(  304)          ! north fold : bmask must be set to 0. on rows jpj-1 and jpj 
(  305)          IF( npolj == 3 .OR. npolj == 4 ) THEN
(  306)             IF( lk_dynspg_rl ) THEN
(  307)                DO ji = 1, nlci
(  308)                   bmask(ji,nlcj-1) = 0.e0
(  309)                   bmask(ji,nlcj  ) = 0.e0
(  310)                END DO
(  311)             ELSE
(  312)                DO ji = 1, nlci
(  313)                   ii = ji + nimpp - 1
(  314)                   bmask(ji,nlcj-1) = bmask(ji,nlcj-1) * tpol(ii)
(  315)                   bmask(ji,nlcj  ) = 0.e0
(  316)                END DO
(  317)             ENDIF
(  318)          ENDIF
(  319)          IF( npolj == 5 .OR. npolj == 6 ) THEN
(  320)             IF( lk_dynspg_rl ) THEN
(  321)                DO ji = 1, nlci
(  322)                   ii = ji + nimpp - 1
(  323)                   bmask(ji,nlcj-1) = bmask(ji,nlcj-1) * fpol(ii)
(  324)                   bmask(ji,nlcj  ) = 0.e0
(  325)                END DO
(  326)             ELSE
(  327)                DO ji = 1, nlci
(  328)                   bmask(ji,nlcj  ) = 0.e0
(  329)                END DO
(  330)             ENDIF
(  331)          ENDIF
(  332)       ENDIF
(  333) 
(  334) 
(  335)       ! mask for second order calculation of vorticity
(  336)       ! ----------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 7

(  337)       
(  338)       CALL dom_msk_nsa
(  339) 
(  340)       
(  341)       ! Lateral boundary conditions on velocity (modify fmask)
(  342)       ! ---------------------------------------
(  343)       
(  344)       DO jk = 1, jpk
(  345) 
(  346)          zwf(:,:) = fmask(:,:,jk)
(  347)          
(  348)          DO jj = 2, jpjm1
(  349)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  350)                IF( fmask(ji,jj,jk) == 0. ) THEN
(  351)                   fmask(ji,jj,jk) = shlat * MIN( 1., MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
(  352)                      &                                    zwf(ji-1,jj), zwf(ji,jj-1)  )  )
(  353)                ENDIF
(  354)             END DO
(  355)          END DO
(  356)          
(  357)          DO jj = 2, jpjm1
(  358)             IF( fmask(1,jj,jk) == 0. ) THEN
(  359)                fmask(1  ,jj,jk) = shlat * MIN( 1., MAX( zwf(2,jj), zwf(1,jj+1), zwf(1,jj-1) ) )
(  360)             ENDIF
(  361)             IF( fmask(jpi,jj,jk) == 0. ) THEN
(  362)                fmask(jpi,jj,jk) = shlat * MIN( 1., MAX( zwf(jpi,jj+1), zwf(jpim1,jj), zwf(jpi,jj-1) ) )
(  363)             ENDIF
(  364)          END DO
(  365)          
(  366)          DO ji = 2, jpim1
(  367)             IF( fmask(ji,1,jk) == 0. ) THEN
(  368)                fmask(ji, 1 ,jk) = shlat * MIN( 1., MAX( zwf(ji+1,1), zwf(ji,2), zwf(ji-1,1) ) )
(  369)             ENDIF
(  370)             IF( fmask(ji,jpj,jk) == 0. ) THEN
(  371)                fmask(ji,jpj,jk) = shlat * MIN( 1., MAX( zwf(ji+1,jpj), zwf(ji-1,jpj), zwf(ji,jpjm1) ) )
(  372)             ENDIF
(  373)          END DO
(  374)       END DO
(  375)       
(  376) 
(  377) !!DB: delete ORCA
(  378) !      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN
(  379) 
(  380)       
(  381)       ! Lateral boundary conditions on fmask
(  382)       CALL lbc_lnk( fmask, 'F', 1. )
(  383)       
(  384)       ! Mbathy set to the number of w-level (minimum value 2)
(  385)       ! -----------------------------------
(  386)       IF( lk_isl ) THEN
(  387)          ! this is done at the end of solver_init routine
(  388)       ELSE
(  389)          DO jj = 1, jpj
(  390)             DO ji = 1, jpi
(  391)                mbathy(ji,jj) = MAX( 1, mbathy(ji,jj) ) + 1
(  392)             END DO
(  393)          END DO
(  394)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 8

(  395)       
(  396)       ! Control print
(  397)       ! -------------
(  398)       IF( nprint == 1 .AND. lwp ) THEN
(  399)          imsk(:,:) = INT( tmask_i(:,:) )
(  400)          WRITE(numout,*) ' tmask_i : '
(  401)          CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  402)                &                           1, jpj, 1, 1, numout)
(  403)          WRITE (numout,*)
(  404)          WRITE (numout,*) ' dommsk: tmask for each level'
(  405)          WRITE (numout,*) ' ----------------------------'
(  406)          DO jk = 1, jpk
(  407)             imsk(:,:) = INT( tmask(:,:,jk) )
(  408) 
(  409)             WRITE(numout,*)
(  410)             WRITE(numout,*) ' level = ',jk
(  411)             CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  412)                &                              1, jpj, 1, 1, numout)
(  413)          END DO
(  414)          WRITE(numout,*)
(  415)          WRITE(numout,*) ' dom_msk: vmask for each level'
(  416)          WRITE(numout,*) ' -----------------------------'
(  417)          DO jk = 1, jpk
(  418)             imsk(:,:) = INT( vmask(:,:,jk) )
(  419)             WRITE(numout,*)
(  420)             WRITE(numout,*) ' level = ',jk
(  421)             CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  422)                &                              1, jpj, 1, 1, numout)
(  423)          END DO
(  424)          WRITE(numout,*)
(  425)          WRITE(numout,*) ' dom_msk: fmask for each level'
(  426)          WRITE(numout,*) ' -----------------------------'
(  427)          DO jk = 1, jpk
(  428)             imsk(:,:) = INT( fmask(:,:,jk) )
(  429)             WRITE(numout,*)
(  430)             WRITE(numout,*) ' level = ',jk
(  431)             CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  432)                &                              1, jpj, 1, 1, numout )
(  433)          END DO
(  434)          WRITE(numout,*)
(  435)          WRITE(numout,*) ' dom_msk: bmask '
(  436)          WRITE(numout,*) ' ---------------'
(  437)          WRITE(numout,*)
(  438)          imsk(:,:) = INT( bmask(:,:) )
(  439)          CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  440)                &                           1, jpj, 1, 1, numout )
(  441)       ENDIF
(  442) 
(  443)    END SUBROUTINE dom_msk
(  444) 
(  445) #if defined key_noslip_accurate
(  446)    !!----------------------------------------------------------------------
(  447)    !!   'key_noslip_accurate' :         accurate no-slip boundary condition
(  448)    !!----------------------------------------------------------------------
(  449)    
(  450)    SUBROUTINE dom_msk_nsa
(  451)       !!---------------------------------------------------------------------
(  452)       !!                 ***  ROUTINE dom_msk_nsa  ***






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 9

(  453)       !! 
(  454)       !! ** Purpose :
(  455)       !!
(  456)       !! ** Method  :
(  457)       !!
(  458)       !! ** Action :
(  459)       !!
(  460)       !! History :
(  461)       !!        !  00-03  (G. Madec)  no slip accurate
(  462)       !!----------------------------------------------------------------------
(  463)       !! *Local declarations
(  464)       INTEGER  :: ji, jj, jk, ii, jl ! dummy loop indices
(  465)       INTEGER ::   ine, inw, ins, inn, itest, ierror, iind, ijnd
(  466)       INTEGER, DIMENSION(jpi*jpj*jpk,3) ::  icoord
(  467)       REAL(wp) ::   zaa
(  468)       !!---------------------------------------------------------------------
(  469)       !!  OPA 9.0 , LOCEAN-IPSL (2005)
(  470)       !!---------------------------------------------------------------------
(  471)       
(  472) 
(  473)       IF(lwp)WRITE(numout,*)
(  474)       IF(lwp)WRITE(numout,*) 'dom_msk_nsa : noslip accurate boundary condition'
(  475)       IF(lwp)WRITE(numout,*) '~~~~~~~~~~~   using Schchepetkin and O Brian scheme'
(  476)       IF( lk_mpp ) THEN
(  477)          IF(lwp)WRITE(numout,cform_err)
(  478)          IF(lwp)WRITE(numout,*) ' mpp version is not yet implemented'
(  479)          nstop = nstop + 1
(  480)       ENDIF
(  481) 
(  482)       ! mask for second order calculation of vorticity
(  483)       ! ----------------------------------------------
(  484)       ! noslip boundary condition: fmask=1  at convex corner, store
(  485)       ! index of straight coast meshes ( 'west', refering to a coast,
(  486)       ! means west of the ocean, aso)
(  487)       
(  488)       DO jk = 1, jpk
(  489)          DO jl = 1, 4
(  490)             npcoa(jl,jk) = 0
(  491)             DO ji = 1, 2*(jpi+jpj)
(  492)                nicoa(ji,jl,jk) = 0
(  493)                njcoa(ji,jl,jk) = 0
(  494)             END DO
(  495)          END DO
(  496)       END DO
(  497)       
(  498) !      IF( jperio == 2 ) THEN
(  499)       IF( jperio == 2 .AND. lwp) THEN
(  500)          WRITE(numout,*) ' '
(  501)          WRITE(numout,*) ' symetric boundary conditions need special'
(  502)          WRITE(numout,*) ' treatment not implemented. we stop.'
(  503)          STOP
(  504)       ENDIF
(  505)       
(  506)       ! convex corners
(  507)       
(  508)       DO jk = 1, jpkm1
(  509)          DO jj = 1, jpjm1
(  510)             DO ji = 1, jpim1






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 10

(  511)                zaa = tmask(ji  ,jj,jk) + tmask(ji  ,jj+1,jk)   &
(  512)                   &+ tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk)
(  513)                IF( ABS(zaa-3.) <= 0.1 )   fmask(ji,jj,jk) = 1.
(  514)             END DO
(  515)          END DO
(  516)       END DO
(  517) 
(  518)       ! north-south straight coast
(  519) 
(  520)       DO jk = 1, jpkm1
(  521)          inw = 0
(  522)          ine = 0
(  523)          DO jj = 2, jpjm1
(  524)             DO ji = 2, jpim1
(  525)                zaa = tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk)
(  526)                IF( ABS(zaa-2.) <= 0.1 .AND. fmask(ji,jj,jk) == 0 ) THEN
(  527)                   inw = inw + 1
(  528)                   nicoa(inw,1,jk) = ji
(  529)                   njcoa(inw,1,jk) = jj
(  530) !                  IF( nprint == 1 ) WRITE(numout,*) ' west  : ', jk,inw, ji, jj
(  531)                   IF( nprint == 1 .AND. lwp) WRITE(numout,*) ' west  : ', jk,inw, ji, jj
(  532) 
(  533)                ENDIF
(  534)                zaa = tmask(ji,jj,jk) + tmask(ji,jj+1,jk)
(  535)                IF( ABS(zaa-2.) <= 0.1 .AND. fmask(ji,jj,jk) == 0 ) THEN
(  536)                   ine = ine + 1
(  537)                   nicoa(ine,2,jk) = ji
(  538)                   njcoa(ine,2,jk) = jj
(  539) !                  IF( nprint == 1 ) WRITE(numout,*) ' east  : ', jk, ine, ji, jj
(  540)                   IF( nprint == 1 .AND. lwp) WRITE(numout,*) ' east  : ', jk, ine, ji, jj
(  541)                ENDIF
(  542)             END DO
(  543)          END DO
(  544)          npcoa(1,jk) = inw
(  545)          npcoa(2,jk) = ine
(  546)       END DO
(  547) 
(  548)       ! west-east straight coast
(  549) 
(  550)       DO jk = 1, jpkm1
(  551)          ins = 0
(  552)          inn = 0
(  553)          DO jj = 2, jpjm1
(  554)             DO ji =2, jpim1
(  555)                zaa = tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk)
(  556)                IF( ABS(zaa-2.) <= 0.1 .AND. fmask(ji,jj,jk) == 0 ) THEN
(  557)                   ins = ins + 1
(  558)                   nicoa(ins,3,jk) = ji
(  559)                   njcoa(ins,3,jk) = jj
(  560) !                  IF( nprint == 1 ) WRITE(numout,*) ' south : ', jk, ins, ji, jj
(  561)                   IF( nprint == 1 .AND. lwp) WRITE(numout,*) ' south : ', jk, ins, ji, jj
(  562)                ENDIF
(  563)                zaa = tmask(ji+1,jj,jk) + tmask(ji,jj,jk)
(  564)                IF( ABS(zaa-2.) <= 0.1 .AND. fmask(ji,jj,jk) == 0 ) THEN
(  565)                   inn = inn + 1
(  566)                   nicoa(inn,4,jk) = ji
(  567)                   njcoa(inn,4,jk) = jj
(  568)                   IF( nprint == 1 .AND. lwp) WRITE(numout,*) ' north : ', jk, inn, ji, jj






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 11

(  569)                ENDIF
(  570)             END DO
(  571)          END DO
(  572)          npcoa(3,jk) = ins
(  573)          npcoa(4,jk) = inn
(  574)       END DO
(  575) 
(  576)       itest = 2 * ( jpi + jpj )
(  577)       DO jk = 1, jpk
(  578)          IF( npcoa(1,jk) > itest .OR. npcoa(2,jk) > itest .OR.   &
(  579)              npcoa(3,jk) > itest .OR. npcoa(4,jk) > itest ) THEN
(  580)             if(lwp) then
(  581)                WRITE(numout,*)
(  582)                WRITE(numout,*) ' level jk = ',jk
(  583)                WRITE(numout,*) ' straight coast index arraies are too small.:'
(  584)                WRITE(numout,*) ' npe, npw, nps, npn = ', npcoa(1,jk), npcoa(2,jk),   &
(  585)                     &                                     npcoa(3,jk), npcoa(4,jk)
(  586)                WRITE(numout,*) ' 2*(jpi+jpj) = ',itest,'. we stop.'
(  587)             endif
(  588)             STOP   !!bug nstop to be used
(  589)         ENDIF
(  590)       END DO
(  591) 
(  592)       ierror = 0
(  593)       iind = 0
(  594)       ijnd = 0
(  595)       IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) iind = 2
(  596)       IF( nperio == 3 .OR. nperio == 4 .OR. nperio == 5 .OR. nperio == 6 ) ijnd = 2
(  597)       DO jk = 1, jpk
(  598)          DO jl = 1, npcoa(1,jk)
(  599)             IF( nicoa(jl,1,jk)+3 > jpi+iind ) THEN
(  600)                ierror = ierror+1
(  601)                icoord(ierror,1) = nicoa(jl,1,jk)
(  602)                icoord(ierror,2) = njcoa(jl,1,jk)
(  603)                icoord(ierror,3) = jk
(  604)             ENDIF
(  605)          END DO
(  606)          DO jl = 1, npcoa(2,jk)
(  607)             IF(nicoa(jl,2,jk)-2 < 1-iind ) THEN
(  608)                ierror = ierror + 1
(  609)                icoord(ierror,1) = nicoa(jl,2,jk)
(  610)                icoord(ierror,2) = njcoa(jl,2,jk)
(  611)                icoord(ierror,3) = jk
(  612)             ENDIF
(  613)          END DO
(  614)          DO jl = 1, npcoa(3,jk)
(  615)             IF( njcoa(jl,3,jk)+3 > jpj+ijnd ) THEN
(  616)                ierror = ierror + 1
(  617)                icoord(ierror,1) = nicoa(jl,3,jk)
(  618)                icoord(ierror,2) = njcoa(jl,3,jk)
(  619)                icoord(ierror,3) = jk
(  620)             ENDIF
(  621)          END DO
(  622)          DO jl=1,npcoa(4,jk)
(  623)             IF( njcoa(jl,4,jk)-2 < 1) THEN
(  624)                ierror=ierror+1
(  625)                icoord(ierror,1)=nicoa(jl,4,jk)
(  626)                icoord(ierror,2)=njcoa(jl,4,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 12

(  627)                icoord(ierror,3)=jk
(  628)             ENDIF
(  629)          END DO
(  630)       END DO
(  631)       
(  632)       IF( ierror > 0 ) THEN
(  633)          IF(lwp) WRITE(numout,*)
(  634)          IF(lwp) WRITE(numout,*) '              Problem on lateral conditions'
(  635)          IF(lwp) WRITE(numout,*) '                 Bad marking off at points:'
(  636)          DO jl = 1, ierror
(  637)             IF(lwp) WRITE(numout,*) 'Level:',icoord(jl,3),   &
(  638)                &                  '  Point(',icoord(jl,1),',',icoord(jl,2),')'
(  639)          END DO
(  640)          IF(lwp) WRITE(numout,*) 'We stop...'   !!cr print format to be used
(  641)          nstop = nstop + 1
(  642)       ENDIF
(  643) 
(  644)    END SUBROUTINE dom_msk_nsa
(  645) 
(  646) #else
(  647)    !!----------------------------------------------------------------------
(  648)    !!   Default option :                                      Empty routine
(  649)    !!----------------------------------------------------------------------
(  650)    SUBROUTINE dom_msk_nsa       
(  651)    END SUBROUTINE dom_msk_nsa
(  652) #endif
(  653)    
(  654)    !!======================================================================
(  655) END MODULE dommsk



































PGF90 (Version     10.5)          09/23/2010  09:49:32      page 13

(    1) # 1 "dommsk.F90"
(    1) MODULE dommsk
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE dommsk   ***
(    4)    !! Ocean initialization : domain land/sea mask 
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dom_msk        : compute land/ocean mask
(    9)    !!   dom_msk_nsa    : update land/ocean mask when no-slip accurate
(   10)    !!                    option is used.
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE obc_oce         ! ocean open boundary conditions
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18)    USE lib_mpp
(   19)    USE solisl          ! ???
(   20)    USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC dom_msk        ! routine called by inidom.F90
(   27) 
(   28)    !! * Module variables
(   29)    REAL(wp) ::   &
(   30)       shlat = 2.   ! type of lateral boundary condition on velocity (namelist namlbc)
(   31)    
(   32)    !! * Substitutions
(   33) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 34 "dommsk.F90"
(   34) # 34
(   34)    !!---------------------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/dommsk.F90,v 1.10 2006/03/21 08:25:42 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!---------------------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41)    
(   42)    SUBROUTINE dom_msk
(   43)       !!---------------------------------------------------------------------
(   44)       !!                 ***  ROUTINE dom_msk  ***






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 14

(   45)       !!
(   46)       !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
(   47)       !!      zontal velocity points (u & v), vorticity points (f) and baro-
(   48)       !!      tropic stream function  points (b).
(   49)       !!        Set mbathy to the number of non-zero w-levels of a water column
(   50)       !!      (if island in the domain (lk_isl=T), this is done latter in
(   51)       !!      routine solver_init)
(   52)       !!
(   53)       !! ** Method  :   The ocean/land mask is computed from the basin bathy-
(   54)       !!      metry in level (mbathy) which is defined or read in dommba.
(   55)       !!      mbathy equals 0 over continental T-point, -n over the nth 
(   56)       !!      island T-point, and the number of ocean level over the ocean.
(   57)       !!
(   58)       !!      At a given position (ji,jj,jk) the ocean/land mask is given by:
(   59)       !!      t-point : 0. IF mbathy( ji ,jj) =< 0
(   60)       !!                1. IF mbathy( ji ,jj) >= jk
(   61)       !!      u-point : 0. IF mbathy( ji ,jj)  or mbathy(ji+1, jj ) =< 0
(   62)       !!                1. IF mbathy( ji ,jj) and mbathy(ji+1, jj ) >= jk.
(   63)       !!      v-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1) =< 0
(   64)       !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1) >= jk.
(   65)       !!      f-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1)
(   66)       !!                   or mbathy(ji+1,jj)  or mbathy(ji+1,jj+1) =< 0
(   67)       !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1)
(   68)       !!                and mbathy(ji+1,jj) and mbathy(ji+1,jj+1) >= jk.
(   69)       !!      b-point : the same definition as for f-point of the first ocean
(   70)       !!                level (surface level) but with 0 along coastlines.
(   71)       !!
(   72)       !!        The lateral friction is set through the value of fmask along
(   73)       !!      the coast and topography. This value is defined by shlat, a
(   74)       !!      namelist parameter:
(   75)       !!         shlat = 0, free slip  (no shear along the coast)
(   76)       !!         shlat = 2, no slip  (specified zero velocity at the coast)
(   77)       !!         0 < shlat < 2, partial slip   | non-linear velocity profile
(   78)       !!         2 < shlat, strong slip        | in the lateral boundary layer
(   79)       !!
(   80)       !!      N.B. If nperio not equal to 0, the land/ocean mask arrays
(   81)       !!      are defined with the proper value at lateral domain boundaries,
(   82)       !!      but bmask. indeed, bmask defined the domain over which the
(   83)       !!      barotropic stream function is computed. this domain cannot
(   84)       !!      contain identical columns because the matrix associated with
(   85)       !!      the barotropic stream function equation is then no more inverti-
(   86)       !!      ble. therefore bmask is set to 0 along lateral domain boundaries
(   87)       !!      even IF nperio is not zero.
(   88)       !!
(   89)       !!      In case of open boundaries (lk_obc=T):
(   90)       !!        - tmask is set to 1 on the points to be computed bay the open
(   91)       !!          boundaries routines.
(   92)       !!        - bmask is  set to 0 on the open boundaries.
(   93)       !!
(   94)       !!      Set mbathy to the number of non-zero w-levels of a water column
(   95)       !!                  mbathy = min( mbathy, 1 ) + 1
(   96)       !!      (note that the minimum value of mbathy is 2).
(   97)       !!
(   98)       !! ** Action :
(   99)       !!                     tmask    : land/ocean mask at t-point (=0. or 1.)
(  100)       !!                     umask    : land/ocean mask at u-point (=0. or 1.)
(  101)       !!                     vmask    : land/ocean mask at v-point (=0. or 1.)
(  102)       !!                     fmask    : land/ocean mask at f-point (=0. or 1.)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 15

(  103)       !!                          =shlat along lateral boundaries
(  104)       !!                     bmask    : land/ocean mask at barotropic stream
(  105)       !!                          function point (=0. or 1.) and set to
(  106)       !!                          0 along lateral boundaries
(  107)       !!                   mbathy   : number of non-zero w-levels 
(  108)       !!
(  109)       !! History :
(  110)       !!        !  87-07  (G. Madec)  Original code
(  111)       !!        !  91-12  (G. Madec)
(  112)       !!        !  92-06  (M. Imbard)
(  113)       !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
(  114)       !!        !  96-01  (G. Madec)  suppression of common work arrays
(  115)       !!        !  96-05  (G. Madec)  mask computed from tmask and sup-
(  116)       !!                 pression of the double computation of bmask
(  117)       !!        !  97-02  (G. Madec)  mesh information put in domhgr.F
(  118)       !!        !  97-07  (G. Madec)  modification of mbathy and fmask
(  119)       !!        !  98-05  (G. Roullet)  free surface
(  120)       !!        !  00-03  (G. Madec)  no slip accurate
(  121)       !!        !  01-09  (J.-M. Molines)  Open boundaries
(  122)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  123)       !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(  124)       !!----------------------------------------------------------------------
(  125)       !! *Local declarations
(  126)       INTEGER  ::   ji, jj, jk, ii     ! dummy loop indices
(  127)       INTEGER  ::   iif, iil, ijf, ijl
(  128)       INTEGER  ::   ii0, ii1, ij0, ij1
(  129)       INTEGER, DIMENSION(jpi,jpj) ::  imsk
(  130) 
(  131)       REAL(wp), DIMENSION(jpi,jpj) ::   zwf
(  132) 
(  133)       NAMELIST/namlbc/ shlat
(  134)       !!---------------------------------------------------------------------
(  135)       
(  136) 
(  137)       ! Namelist namlbc : lateral momentum boundary condition
(  138)       REWIND( numnam )
(  139)       READ  ( numnam, namlbc )
(  140)       IF(lwp) THEN
(  141)          WRITE(numout,*)
(  142)          WRITE(numout,*) 'dommsk : ocean mask '
(  143)          WRITE(numout,*) '~~~~~~'
(  144)          WRITE(numout,*) '         Namelist namlbc'
(  145)          WRITE(numout,*) '            lateral momentum boundary cond. shlat = ',shlat
(  146)       ENDIF
(  147) 
(  148)       IF ( shlat == 0. ) THEN
(  149)           IF(lwp) WRITE(numout,*) '         ocean lateral free-slip '
(  150)         ELSEIF ( shlat  ==  2. ) THEN
(  151)           IF(lwp) WRITE(numout,*) '         ocean lateral  no-slip '
(  152)         ELSEIF ( 0. < shlat .AND. shlat < 2. ) THEN
(  153)           IF(lwp) WRITE(numout,*) '         ocean lateral  partial-slip '
(  154)         ELSEIF ( 2. < shlat ) THEN
(  155)           IF(lwp) WRITE(numout,*) '         ocean lateral  strong-slip '
(  156)         ELSE
(  157)           IF(lwp) WRITE(numout,cform_err)
(  158)           IF(lwp) WRITE(numout,*) ' shlat is negative = ', shlat
(  159)           nstop = nstop + 1
(  160)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 16

(  161) 
(  162)       ! 1. Ocean/land mask at t-point (computed from mbathy)
(  163)       ! -----------------------------
(  164)       ! Tmask has already the right boundary conditions since mbathy is ok
(  165) 
(  166)       tmask(:,:,:) = 0.e0
(  167)       DO jk = 1, jpk
(  168)          DO jj = 1, jpj
(  169)             DO ji = 1, jpi
(  170)                IF( FLOAT( mbathy(ji,jj)-jk )+.1 >= 0.e0 ) tmask(ji,jj,jk) = 1.e0
(  171)             END DO  
(  172)          END DO  
(  173)       END DO  
(  174) 
(  175) !ylu
(  176) !!DB -- 2008.11.24 -- Is this the cause of my probs in SE/SW corners ???
(  177) !SPINUP Run ---> garbage ~10d
(  178) !      write(2200+narea,*)'DBG: not setting tmask in SE/SW corners to zero'
(  179)       tmask(mi0(1):mi1(2),mj0(1):mj1(2),1:jpk)=0.0
(  180)       tmask(mi0(196):mi1(197),mj0(1):mj1(2),1:jpk)=0.0
(  181) 
(  182) 
(  183) 
(  184) # 188
(  188)       ! Interior domain mask (used for global sum)
(  189)       ! --------------------
(  190) 
(  191)       tmask_i(:,:) = tmask(:,:,1)
(  192)       iif = jpreci                         ! ???
(  193)       iil = nlci - jpreci + 1
(  194)       ijf = jprecj                         ! ???
(  195)       ijl = nlcj - jprecj + 1
(  196) 
(  197)       tmask_i( 1 :iif,   :   ) = 0.e0      ! first columns
(  198)       tmask_i(iil:jpi,   :   ) = 0.e0      ! last  columns (including mpp extra columns)
(  199)       tmask_i(   :   , 1 :ijf) = 0.e0      ! first rows
(  200)       tmask_i(   :   ,ijl:jpj) = 0.e0      ! last  rows (including mpp extra rows)
(  201) 
(  202)       ! north fold mask
(  203)       tpol(1:jpiglo) = 1.e0 
(  204)       fpol(1:jpiglo) = 1.e0
(  205)       IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
(  206)          tpol(jpiglo/2+1:jpiglo) = 0.e0
(  207)          fpol(     1    :jpiglo) = 0.e0
(  208)          ! T-point pivot: only half of the nlcj-1 row
(  209)          IF( mjg(nlej) == jpjglo )   THEN
(  210)             DO ji = iif+1, iil-1
(  211)                tmask_i(ji,nlej-1) = tmask_i(ji,nlej-1) * tpol(mig(ji))
(  212)             END DO
(  213)          ENDIF
(  214)       ENDIF
(  215)       IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
(  216)          tpol(     1    :jpiglo) = 0.e0
(  217)          fpol(jpiglo/2+1:jpiglo) = 0.e0
(  218)       ENDIF
(  219) 
(  220)       ! 2. Ocean/land mask at u-,  v-, and z-points (computed from tmask)
(  221)       ! -------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 17

(  222)       
(  223)       ! Computation
(  224)       DO jk = 1, jpk
(  225)          DO jj = 1, jpjm1
(  226)             DO ji = 1, jpim1   ! vector loop
(  227)                umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
(  228)                vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
(  229)                fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
(  230)                   &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
(  231)             END DO
(  232)          END DO
(  233)       END DO
(  234) 
(  235) !!DB: delete ORCA
(  236) !      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN
(  237) 
(  238)       ! Lateral boundary conditions
(  239)       CALL lbc_lnk( umask, 'U', 1. )
(  240)       CALL lbc_lnk( vmask, 'V', 1. )
(  241)       CALL lbc_lnk( fmask, 'F', 1. )
(  242) 
(  243)       nmask(:,:,:)=fmask(:,:,:)
(  244)       ! 4. ocean/land mask for the elliptic equation
(  245)       ! --------------------------------------------
(  246)       
(  247)       ! Computation
(  248)       IF( lk_dynspg_rl ) THEN
(  249)          bmask(:,:) = fmask(:,:,1)       ! elliptic equation is written at f-point
(  250)       ELSE
(  251)          bmask(:,:) = tmask(:,:,1)       ! elliptic equation is written at t-point
(  252)       ENDIF
(  253)       
(  254)       ! Boundary conditions
(  255)       !   cyclic east-west : bmask must be set to 0. on rows 1 and jpi
(  256)       IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
(  257)          bmask( 1 ,:) = 0.e0
(  258)          bmask(jpi,:) = 0.e0
(  259)       ENDIF
(  260)       
(  261)       !   south symmetric :  bmask must be set to 0. on row 1
(  262)       IF( nperio == 2 ) THEN
(  263)          bmask(:, 1 ) = 0.e0
(  264)       ENDIF
(  265)       
(  266)       !   north fold : 
(  267)       IF( nperio == 3 .OR. nperio == 4 ) THEN
(  268)          IF( lk_dynspg_rl ) THEN
(  269)             ! T-pt pivot and F-pt elliptic eq. : bmask set to 0. on rows jpj-1 and jpj
(  270)             bmask(:,jpj-1) = 0.e0
(  271)             bmask(:,jpj  ) = 0.e0
(  272)          ELSE
(  273)             ! T-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj and on half jpjglo-1 row
(  274)             DO ji = 1, jpi
(  275)                ii = ji + nimpp - 1
(  276)                bmask(ji,jpj-1) = bmask(ji,jpj-1) * tpol(ii)
(  277)                bmask(ji,jpj  ) = 0.e0
(  278)             END DO
(  279)          ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 18

(  280)       ENDIF
(  281)       IF( nperio == 5 .OR. nperio == 6 ) THEN
(  282)          IF( lk_dynspg_rl ) THEN
(  283)             ! F-pt pivot and F-pt elliptic eq. : bmask set to 0. on row jpj and on half jpjglo-1 row
(  284)             DO ji = 1, jpi
(  285)                ii = ji + nimpp - 1
(  286)                bmask(ji,jpj-1) = bmask(ji,jpj-1) * fpol(ii)
(  287)                bmask(ji,jpj  ) = 0.e0
(  288)             END DO
(  289)          ELSE
(  290)             ! F-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj
(  291)             bmask(:,jpj) = 0.e0
(  292)          ENDIF
(  293)       ENDIF
(  294) 
(  295)       ! Mpp boundary conditions: bmask is set to zero on the overlap
(  296)       ! region for all elliptic solvers
(  297) 
(  298)       IF( lk_mpp ) THEN
(  299)          IF( nbondi /= -1 .AND. nbondi /= 2 )   bmask(  1 :jpreci,:) = 0.e0
(  300)          IF( nbondi /=  1 .AND. nbondi /= 2 )   bmask(nlci:jpi   ,:) = 0.e0
(  301)          IF( nbondj /= -1 .AND. nbondj /= 2 )   bmask(:,  1 :jprecj) = 0.e0
(  302)          IF( nbondj /=  1 .AND. nbondj /= 2 )   bmask(:,nlcj:jpj   ) = 0.e0
(  303)       
(  304)          ! north fold : bmask must be set to 0. on rows jpj-1 and jpj 
(  305)          IF( npolj == 3 .OR. npolj == 4 ) THEN
(  306)             IF( lk_dynspg_rl ) THEN
(  307)                DO ji = 1, nlci
(  308)                   bmask(ji,nlcj-1) = 0.e0
(  309)                   bmask(ji,nlcj  ) = 0.e0
(  310)                END DO
(  311)             ELSE
(  312)                DO ji = 1, nlci
(  313)                   ii = ji + nimpp - 1
(  314)                   bmask(ji,nlcj-1) = bmask(ji,nlcj-1) * tpol(ii)
(  315)                   bmask(ji,nlcj  ) = 0.e0
(  316)                END DO
(  317)             ENDIF
(  318)          ENDIF
(  319)          IF( npolj == 5 .OR. npolj == 6 ) THEN
(  320)             IF( lk_dynspg_rl ) THEN
(  321)                DO ji = 1, nlci
(  322)                   ii = ji + nimpp - 1
(  323)                   bmask(ji,nlcj-1) = bmask(ji,nlcj-1) * fpol(ii)
(  324)                   bmask(ji,nlcj  ) = 0.e0
(  325)                END DO
(  326)             ELSE
(  327)                DO ji = 1, nlci
(  328)                   bmask(ji,nlcj  ) = 0.e0
(  329)                END DO
(  330)             ENDIF
(  331)          ENDIF
(  332)       ENDIF
(  333) 
(  334) 
(  335)       ! mask for second order calculation of vorticity
(  336)       ! ----------------------------------------------
(  337)       






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 19

(  338)       CALL dom_msk_nsa
(  339) 
(  340)       
(  341)       ! Lateral boundary conditions on velocity (modify fmask)
(  342)       ! ---------------------------------------
(  343)       
(  344)       DO jk = 1, jpk
(  345) 
(  346)          zwf(:,:) = fmask(:,:,jk)
(  347)          
(  348)          DO jj = 2, jpjm1
(  349)             DO ji = 2, jpim1   ! vector opt.
(  350)                IF( fmask(ji,jj,jk) == 0. ) THEN
(  351)                   fmask(ji,jj,jk) = shlat * MIN( 1., MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
(  352)                      &                                    zwf(ji-1,jj), zwf(ji,jj-1)  )  )
(  353)                ENDIF
(  354)             END DO
(  355)          END DO
(  356)          
(  357)          DO jj = 2, jpjm1
(  358)             IF( fmask(1,jj,jk) == 0. ) THEN
(  359)                fmask(1  ,jj,jk) = shlat * MIN( 1., MAX( zwf(2,jj), zwf(1,jj+1), zwf(1,jj-1) ) )
(  360)             ENDIF
(  361)             IF( fmask(jpi,jj,jk) == 0. ) THEN
(  362)                fmask(jpi,jj,jk) = shlat * MIN( 1., MAX( zwf(jpi,jj+1), zwf(jpim1,jj), zwf(jpi,jj-1) ) )
(  363)             ENDIF
(  364)          END DO
(  365)          
(  366)          DO ji = 2, jpim1
(  367)             IF( fmask(ji,1,jk) == 0. ) THEN
(  368)                fmask(ji, 1 ,jk) = shlat * MIN( 1., MAX( zwf(ji+1,1), zwf(ji,2), zwf(ji-1,1) ) )
(  369)             ENDIF
(  370)             IF( fmask(ji,jpj,jk) == 0. ) THEN
(  371)                fmask(ji,jpj,jk) = shlat * MIN( 1., MAX( zwf(ji+1,jpj), zwf(ji-1,jpj), zwf(ji,jpjm1) ) )
(  372)             ENDIF
(  373)          END DO
(  374)       END DO
(  375)       
(  376) 
(  377) !!DB: delete ORCA
(  378) !      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN
(  379) 
(  380)       
(  381)       ! Lateral boundary conditions on fmask
(  382)       CALL lbc_lnk( fmask, 'F', 1. )
(  383)       
(  384)       ! Mbathy set to the number of w-level (minimum value 2)
(  385)       ! -----------------------------------
(  386)       IF( lk_isl ) THEN
(  387)          ! this is done at the end of solver_init routine
(  388)       ELSE
(  389)          DO jj = 1, jpj
(  390)             DO ji = 1, jpi
(  391)                mbathy(ji,jj) = MAX( 1, mbathy(ji,jj) ) + 1
(  392)             END DO
(  393)          END DO
(  394)       ENDIF
(  395)       






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 20

(  396)       ! Control print
(  397)       ! -------------
(  398)       IF( nprint == 1 .AND. lwp ) THEN
(  399)          imsk(:,:) = INT( tmask_i(:,:) )
(  400)          WRITE(numout,*) ' tmask_i : '
(  401)          CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  402)                &                           1, jpj, 1, 1, numout)
(  403)          WRITE (numout,*)
(  404)          WRITE (numout,*) ' dommsk: tmask for each level'
(  405)          WRITE (numout,*) ' ----------------------------'
(  406)          DO jk = 1, jpk
(  407)             imsk(:,:) = INT( tmask(:,:,jk) )
(  408) 
(  409)             WRITE(numout,*)
(  410)             WRITE(numout,*) ' level = ',jk
(  411)             CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  412)                &                              1, jpj, 1, 1, numout)
(  413)          END DO
(  414)          WRITE(numout,*)
(  415)          WRITE(numout,*) ' dom_msk: vmask for each level'
(  416)          WRITE(numout,*) ' -----------------------------'
(  417)          DO jk = 1, jpk
(  418)             imsk(:,:) = INT( vmask(:,:,jk) )
(  419)             WRITE(numout,*)
(  420)             WRITE(numout,*) ' level = ',jk
(  421)             CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  422)                &                              1, jpj, 1, 1, numout)
(  423)          END DO
(  424)          WRITE(numout,*)
(  425)          WRITE(numout,*) ' dom_msk: fmask for each level'
(  426)          WRITE(numout,*) ' -----------------------------'
(  427)          DO jk = 1, jpk
(  428)             imsk(:,:) = INT( fmask(:,:,jk) )
(  429)             WRITE(numout,*)
(  430)             WRITE(numout,*) ' level = ',jk
(  431)             CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  432)                &                              1, jpj, 1, 1, numout )
(  433)          END DO
(  434)          WRITE(numout,*)
(  435)          WRITE(numout,*) ' dom_msk: bmask '
(  436)          WRITE(numout,*) ' ---------------'
(  437)          WRITE(numout,*)
(  438)          imsk(:,:) = INT( bmask(:,:) )
(  439)          CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
(  440)                &                           1, jpj, 1, 1, numout )
(  441)       ENDIF
(  442) 
(  443)    END SUBROUTINE dom_msk
(  444) 
(  445) # 647
(  647)    !!----------------------------------------------------------------------
(  648)    !!   Default option :                                      Empty routine
(  649)    !!----------------------------------------------------------------------
(  650)    SUBROUTINE dom_msk_nsa       
(  651)    END SUBROUTINE dom_msk_nsa
(  652) # 653
(  653)    
(  654)    !!======================================================================






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 21

(  655) END MODULE dommsk
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (dommsk.F90: 170)
