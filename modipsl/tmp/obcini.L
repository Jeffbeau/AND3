


PGF90 (Version     10.5)          09/23/2010  09:50:07      page 1

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

Filename: obcini.F90

(    1) !!DB 2008.05.08
(    2) !!(0) In general search for DB and DBG to see my modifications
(    3) !!(1) Regarding the fact that certain CPU tilings COULD result in overlaps
(    4) !!between tiles of greater than 1 column: I have written and retained code
(    5) !!that flags this occurence, and outputs information. Search for:
(    6) !! write(1X00,...)  where X=5,6,7 == E,W,S ----> fort.1500, fort.1600, fort.1700 files
(    7) !!(2) I create 2D global OBC masking files for use in obcdta. They are called:
(    8) !!emaskg2, wmaskg2, smaskg2 (where e.g.: e=east, g=global dimension, 2=2 dimensional).
(    9) !!They are output to: fort.6X00 where X=5,6,7 == E,W,S.
(   10) !!Note the use of mpp_sum() to ensure that all CPUs share the same global values.
(   11) !!(3) If you want output of a lot of debugging info on OBC variable values for
(   12) !!various CPUs (nareas) set DBG=.true., and search DBG and write(####+narea,...)
(   13) 
(   14) 
(   15) MODULE obcini
(   16)   !!=================================================================================
(   17)   !!                       ***  MODULE  obcini  ***
(   18)   !! OBC initial state :  Open boundary initial state
(   19)   !!=================================================================================
(   20) #if defined key_obc
(   21)   !!---------------------------------------------------------------------------------
(   22)   !!   'key_obc'                                             Open Boundary Conditions
(   23)   !!---------------------------------------------------------------------------------
(   24)   !!   obc_init       : initialization for the open boundary condition
(   25)   !!---------------------------------------------------------------------------------
(   26)   !! * Modules used
(   27)   USE oce             ! ocean dynamics and tracers variables
(   28)   USE dom_oce         ! ocean space and time domain variables
(   29)   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   30)   USE phycst          ! physical constants
(   31)   USE obc_oce         ! ocean open boundary conditions
(   32)   USE lib_mpp         ! for mpp_sum
(   33)   USE in_out_manager  ! I/O units
(   34) 
(   35)   IMPLICIT NONE
(   36)   PRIVATE
(   37) 
(   38)   !! * Routine accessibility
(   39)   PUBLIC obc_init        ! routine called by opa.F90
(   40) 
(   41)   !! * Substitutions
(   42) #  include "obc_vectopt_loop_substitute.h90"
(   43)   !!---------------------------------------------------------------------------------
(   44)   !!   OPA 9.0 , LOCEAN-IPSL (2005)
(   45)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcini.F90,v 1.8 2006/03/21 07:52:26 opalod Exp $
(   46)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 2

(   47)   !!---------------------------------------------------------------------------------
(   48) 
(   49) CONTAINS
(   50) 
(   51)   SUBROUTINE obc_init
(   52)      !!----------------------------------------------------------------------
(   53)      !!                 ***  ROUTINE obc_init  ***
(   54)      !!
(   55)      !! ** Purpose :   Initialization of the dynamics and tracer fields at
(   56)      !!      the open boundaries.
(   57)      !!
(   58)      !! ** Method  :   initialization of open boundary variables
(   59)      !!      (u, v, bsf) over 3 time step and 3 rows
(   60)      !!      (t, s)      over 2 time step and 2 rows
(   61)      !!      if ln_rstart = .FALSE. : no restart, fields set to zero
(   62)      !!      if ln_rstart = .TRUE.  : restart, fields are read in a file
(   63)      !!      if rdpxxx = 0 then lfbc is set true for this boundary.
(   64)      !!
(   65)      !! ** Input   :   restart.obc file, restart file for open boundaries
(   66)      !!
(   67)      !! History :
(   68)      !!   8.0  !  97-07  (G. Madec)  Original code
(   69)      !!        !  97-11  (J.M. Molines)
(   70)      !!   8.5  !  02-11  (C. Talandier, A-M. Treguier) Free surface, F90
(   71)      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   72)      !!----------------------------------------------------------------------
(   73)      !! * Modules used
(   74)      USE obcrst,   ONLY :   obc_rst_lec   ! Make obc_rst_lec routine available
(   75)      USE obcdom,   ONLY :   obc_dom       ! Make obc_dom routine available
(   76) 
(   77)      !! * Local declarations
(   78) !sujie      INTEGER  ::   ji, jj, istop , inumfbc
(   79)      INTEGER  ::   ji, jj, istop , inumfbc, j1, ii, ij
(   80)      INTEGER, DIMENSION(4) ::   icorner
(   81)      REAL(wp) ::   zbsic1, zbsic2, zbsic3
(   82)      REAL(wp), DIMENSION(2) ::   ztestmask
(   83) !!DB: 2008.05.01
(   84)      INTEGER :: s_bnd1(jpni) = 0, s_bnd2(jpni) = 0
(   85)      INTEGER :: w_bnd1(jpnj) = 0, w_bnd2(jpnj) = 0
(   86)      INTEGER :: e_bnd1(jpnj) = 0, e_bnd2(jpnj) = 0
(   87)      LOGICAL :: DBG
(   88) 
(   89)      NAMELIST/namobc/ rdpein, rdpwin, rdpnin, rdpsin,   &
(   90)         &             rdpeob, rdpwob, rdpnob, rdpsob,   &
(   91)         &             zbsic1, zbsic2, zbsic3,           &
(   92)         &             nbic, volemp, nobc_dta,           &
(   93)         &             ln_obc_clim, ln_vol_cst, ln_obc_fla
(   94)      !!----------------------------------------------------------------------
(   95) 
(   96) 
(   97) !!DB
(   98) !     DBG = .true.
(   99)      DBG = .false.
(  100) 
(  101) 
(  102)      IF(lwp) WRITE(numout,*)
(  103)      IF(lwp) WRITE(numout,*) 'obc_init : initialization of open boundaries'
(  104)      IF(lwp) WRITE(numout,*) '~~~~~~~~'






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 3

(  105) 
(  106) 
(  107)      ! 0. read namelist parameters
(  108)      ! ---------------------------
(  109)      ! default values already set except:
(  110)      zbsic1 = 0.e0
(  111)      zbsic2 = 0.e0
(  112)      zbsic3 = 0.e0
(  113) 
(  114)      ! Namelist namobc : open boundaries
(  115)      REWIND( numnam )
(  116)      READ  ( numnam, namobc )
(  117) 
(  118)      bsfic0(1) = zbsic1
(  119)      bsfic (2) = zbsic2
(  120)      bsfic (3) = zbsic3
(  121) 
(  122)      ! By security we set rdpxin and rdpxob respectively
(  123)      ! to 1. and 15. if the corresponding OBC is not activated
(  124)      IF( .NOT.lp_obc_east ) THEN
(  125)         rdpein = 1.
(  126)         rdpeob = 15.
(  127)      END IF
(  128)      IF( .NOT.lp_obc_west ) THEN
(  129)         rdpwin = 1.
(  130)         rdpwob = 15.
(  131)      END IF
(  132)      IF( .NOT.lp_obc_north ) THEN
(  133)         rdpnin = 1.
(  134)         rdpnob = 15.
(  135)      END IF
(  136)      IF( .NOT.lp_obc_south ) THEN
(  137)         rdpsin = 1.
(  138)         rdpsob = 15.
(  139)      END IF
(  140) 
(  141)      ! number of open boudaries and open boundary indicators
(  142)      nbobc = 0
(  143)      IF( lp_obc_east  )   nbobc = nbobc + 1
(  144)      IF( lp_obc_west  )   nbobc = nbobc + 1
(  145)      IF( lp_obc_north )   nbobc = nbobc + 1
(  146)      IF( lp_obc_south )   nbobc = nbobc + 1
(  147) 
(  148)      IF(lwp) WRITE(numout,*) '         Number of open boundaries    nbobc = ',nbobc
(  149)      IF(lwp) WRITE(numout,*)
(  150)      IF( nbobc /= 0 .AND. jperio /= 0 ) THEN
(  151)         IF(lwp) WRITE(numout,*)
(  152)         IF(lwp) WRITE(numout,*) ' E R R O R : Cyclic or symmetric,',   &
(  153)            ' and open boundary condition are not compatible'
(  154)         IF(lwp) WRITE(numout,*) ' ========== '
(  155)         IF(lwp) WRITE(numout,*)
(  156)         nstop = nstop + 1
(  157)      END IF
(  158) 
(  159)      ! control prints
(  160)      IF(lwp) WRITE(numout,*) '         namobc'
(  161)      IF(lwp) WRITE(numout,*) ' '
(  162)      IF(lwp) WRITE(numout,*) '         data in file (=1) or     nobc_dta = ', nobc_dta






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 4

(  163)      IF(lwp) WRITE(numout,*) '         initial state used (=0)             '
(  164)      IF(lwp) WRITE(numout,*) '         climatology (true) or not:', ln_obc_clim
(  165)      IF(lwp) WRITE(numout,*) ' '
(  166)      IF(lwp) WRITE(numout,*) '                                 WARNING                     '
(  167)      IF(lwp) WRITE(numout,*) '         Flather"s algorithm is applied with explicit free surface scheme                 '
(  168)      IF(lwp) WRITE(numout,*) '         or with free surface time-splitting scheme          '
(  169)      IF(lwp) WRITE(numout,*) '         Nor radiation neither relaxation is allowed with explicit free surface scheme:   '
(  170)      IF(lwp) WRITE(numout,*) '         Radiation and/or relaxation is allowed with free surface time-splitting scheme '
(  171)      IF(lwp) WRITE(numout,*) '         depending of the choice of rdpXin = rdpXob  = 0. for open boundaries             '
(  172)      IF(lwp) WRITE(numout,*) ' '
(  173)      IF(lwp) WRITE(numout,*) '         For the rigid-lid case or the filtered free surface case,                        '
(  174)      IF(lwp) WRITE(numout,*) '         radiation, relaxation or presciption of data can be applied                      '
(  175)      IF( lwp.AND.lp_obc_east ) THEN
(  176)         WRITE(numout,*) '         East open boundary :'
(  177)         WRITE(numout,*) '              i index                    jpieob = ', jpieob
(  178)         WRITE(numout,*) '              damping time scale (days)  rdpeob = ', rdpeob
(  179)         WRITE(numout,*) '              damping time scale (days)  rdpein = ', rdpein
(  180)      ENDIF
(  181) 
(  182)      IF( lwp.AND.lp_obc_west ) THEN
(  183)         WRITE(numout,*) '         West open boundary :'
(  184)         WRITE(numout,*) '              i index                    jpiwob = ', jpiwob
(  185)         WRITE(numout,*) '              damping time scale (days)  rdpwob = ', rdpwob
(  186)         WRITE(numout,*) '              damping time scale (days)  rdpwin = ', rdpwin
(  187)      ENDIF
(  188) 
(  189)      IF( lwp.AND.lp_obc_north ) THEN
(  190)         WRITE(numout,*) '         North open boundary :'
(  191)         WRITE(numout,*) '               j index                    jpjnob = ', jpjnob
(  192)         WRITE(numout,*) '               damping time scale (days)  rdpnob = ', rdpnob
(  193)         WRITE(numout,*) '               damping time scale (days)  rdpnin = ', rdpnin
(  194)      ENDIF
(  195) 
(  196)      IF( lwp.AND.lp_obc_south ) THEN
(  197)         WRITE(numout,*) '         South open boundary :'
(  198)         WRITE(numout,*) '               j index                    jpjsob = ', jpjsob
(  199)         WRITE(numout,*) '               damping time scale (days)  rdpsob = ', rdpsob
(  200)         WRITE(numout,*) '               damping time scale (days)  rdpsin = ', rdpsin
(  201)         WRITE(numout,*) ' '
(  202)      ENDIF
(  203) 
(  204)      ! 1. Initialisation of constants
(  205)      ! ------------------------------
(  206) 
(  207)      ! ... convert rdp$ob in seconds
(  208)      rdpein = rdpein * rday
(  209)      rdpwin = rdpwin * rday
(  210)      rdpnin = rdpnin * rday
(  211)      rdpsin = rdpsin * rday
(  212)      rdpeob = rdpeob * rday
(  213)      rdpwob = rdpwob * rday
(  214)      rdpnob = rdpnob * rday
(  215)      rdpsob = rdpsob * rday
(  216)      lfbceast  = .FALSE.
(  217)      lfbcwest  = .FALSE.
(  218)      lfbcnorth = .FALSE.
(  219)      lfbcsouth = .FALSE.
(  220)      inumfbc = 0






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 5

(  221)      ! ... look for Fixed Boundaries (rdp = 0 )
(  222)      ! ... When specified, lbcxxx flags are set to TRUE and rdpxxx are set to
(  223)      ! ...  a small arbitrary value, (to avoid division by zero further on).
(  224)      ! ...  rdpxxx is not used anymore.
(  225)      IF( lp_obc_east )  THEN
(  226)         IF( (rdpein+rdpeob) == 0 )  THEN
(  227)            lfbceast = .TRUE.
(  228)            rdpein = 1e-3
(  229)            rdpeob = 1e-3
(  230)            inumfbc = inumfbc+1
(  231)         ELSEIF ( (rdpein*rdpeob) == 0 )  THEN
(  232)            IF(lwp) THEN
(  233)               WRITE(numout,cform_err)
(  234)               WRITE(numout,*) 'obc_init : rdpein & rdpeob must be both zero or non zero'
(  235)               nstop = nstop + 1
(  236)            ENDIF
(  237)         END IF
(  238)      END IF
(  239)      IF( lp_obc_west )  THEN
(  240)         IF( (rdpwin + rdpwob) == 0 )  THEN
(  241)            lfbcwest = .TRUE.
(  242)            rdpwin = 1e-3
(  243)            rdpwob = 1e-3
(  244)            inumfbc = inumfbc+1
(  245)         ELSEIF ( (rdpwin*rdpwob) == 0 )  THEN
(  246)            IF(lwp) THEN
(  247)               WRITE(numout,cform_err)
(  248)               WRITE(numout,*) 'obc_init : rdpwin & rdpwob must be both zero or non zero'
(  249)               nstop = nstop + 1
(  250)            ENDIF
(  251)         END IF
(  252)      END IF
(  253)      IF( lp_obc_north )  THEN
(  254)         IF( (rdpnin + rdpnob) == 0 )  THEN
(  255)            lfbcnorth = .TRUE.
(  256)            rdpnin = 1e-3
(  257)            rdpnob = 1e-3
(  258)            inumfbc = inumfbc+1
(  259)         ELSEIF ( (rdpnin*rdpnob) == 0 )  THEN
(  260)            IF(lwp) THEN
(  261)               WRITE(numout,cform_err)
(  262)               WRITE(numout,*) 'obc_init : rdpnin & rdpnob must be both zero or non zero'
(  263)               nstop = nstop + 1
(  264)            ENDIF
(  265)         END IF
(  266)      END IF
(  267)      IF( lp_obc_south )  THEN
(  268)         IF( (rdpsin + rdpsob) == 0 )  THEN
(  269)            lfbcsouth = .TRUE.
(  270)            rdpsin = 1e-3
(  271)            rdpsob = 1e-3
(  272)            inumfbc = inumfbc+1
(  273)         ELSEIF ( (rdpsin*rdpsob) == 0 )  THEN
(  274)            IF(lwp) THEN
(  275)               WRITE(numout,cform_err)
(  276)               WRITE(numout,*) 'obc_init : rdpsin & rdpsob must be both zero or non zero'
(  277)               nstop = nstop + 1
(  278)            ENDIF






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 6

(  279)         END IF
(  280)      END IF
(  281) 
(  282)      ! 2.  Clever mpp indices for loops on the open boundaries.
(  283)      !     The loops will be performed only on the processors
(  284)      !     that contain a given open boundary.
(  285)      ! --------------------------------------------------------
(  286) 
(  287)      IF( lp_obc_east ) THEN
(  288)         ! ...   mpp initialization
(  289)         nie0   = max( 1, min(jpieob   - nimpp+1, jpi     ) )
(  290)         nie1   = max( 0, min(jpieob   - nimpp+1, jpi - 1 ) )
(  291)         nie0p1 = max( 1, min(jpieob+1 - nimpp+1, jpi     ) )
(  292)         nie1p1 = max( 0, min(jpieob+1 - nimpp+1, jpi - 1 ) )
(  293)         nie0m1 = max( 1, min(jpieob-1 - nimpp+1, jpi     ) )
(  294)         nie1m1 = max( 0, min(jpieob-1 - nimpp+1, jpi - 1 ) )
(  295)         nje0   = max( 2, min(jpjed    - njmpp+1, jpj     ) )
(  296)         nje1   = max( 0, min(jpjef    - njmpp+1, jpj - 1 ) )
(  297)         nje0p1 = max( 1, min(jpjedp1  - njmpp+1, jpj     ) )
(  298)         nje0m1 = max( 1, min(jpjed    - njmpp+1, jpj     ) )
(  299)         nje1m1 = max( 0, min(jpjefm1  - njmpp+1, jpj - 1 ) )
(  300)         nje1m2 = max( 0, min(jpjefm1-1- njmpp+1, jpj - 1 ) )
(  301)         IF(lwp) THEN
(  302)            IF( lfbceast ) THEN
(  303)               WRITE(numout,*)'     '
(  304)               WRITE(numout,*)'         Specified East Open Boundary'
(  305)            ELSE
(  306)               WRITE(numout,*)'     '
(  307)               WRITE(numout,*)'         Radiative East Open Boundary'
(  308)            END IF
(  309)         END IF
(  310) 
(  311) !DB: output above plus nimpp njmpp in a readable form:
(  312)         if(DBG) then
(  313)            write(2100+narea,'(a35,6(i3,1x))')'OBC east variables for narea ',narea
(  314)            write(2100+narea,'(a35,6(i3,1x))')'nimpp njmpp ', nimpp, njmpp
(  315)            write(2100+narea,'(a35,6(i3,1x))')'jpi jpj jpieob jpjed jpjef ',jpi,jpj, jpieob, jpjed, jpjef 
(  316)            write(2100+narea,'(a35,6(i3,1x))')'nie0, nie1, nje0, nje1: ',nie0, nie1, nje0, nje1
(  317)            write(2100+narea,'(a35,6(i3,1x))')'nje0p1, nje1m1, nie0p1, nie1m1: ',nje0p1, nje1m1, nie0p1, nie1m1
(  318)            write(2100+narea,'(a35,6(i3,1x))')'nje0m1, nje1m2, nie0m1, nie1p1: ',nje0m1, nje1m2, nie0m1, nie1p1
(  319)         endif
(  320) 
(  321) !!DBG
(  322)         if(mod(narea,jpni) == 0) then        !!Assume tiling goes from left-to-right btm-to-top
(  323)            ji = (narea-1)/jpni + 1
(  324)            e_bnd1(ji) = nje0p1-1+njmpp
(  325)            e_bnd2(ji) = nje1m1-1+njmpp
(  326)         endif
(  327) 
(  328)      END IF
(  329)      if(lk_mpp) then
(  330)         call mpp_sum(e_bnd1,jpnj)
(  331)         call mpp_sum(e_bnd2,jpnj)
(  332)      endif
(  333) !!DB: Check if any tiles overlap by more than 1 and report
(  334) !!NB: Assumes halo zone = 1
(  335)      if(lwp) then
(  336) !!DBG






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 7

(  337)         write(1500,'(A30,100(i4,1x))')'(EAST) 1st index in tile:  ',(e_bnd1(ji),ji=1,jpnj)
(  338)         write(1500,'(A30,100(i4,1x))')'(EAST) Last index in tile: ',(e_bnd2(ji),ji=1,jpnj)
(  339)         do ji = 2, jpnj
(  340)            if(e_bnd1(ji) < e_bnd2(ji-1)) then
(  341)               write(1500,*)'Problem at East boundary between tiles ',ji, ji-1
(  342)            endif
(  343)         enddo
(  344)      endif
(  345) 
(  346) 
(  347) 
(  348)      IF( lp_obc_west ) THEN
(  349)         ! ...   mpp initialization
(  350)         niw0   = max( 1, min(jpiwob   - nimpp+1, jpi     ) )
(  351)         niw1   = max( 0, min(jpiwob   - nimpp+1, jpi - 1 ) )
(  352)         niw0p1 = max( 1, min(jpiwob+1 - nimpp+1, jpi     ) )
(  353)         niw1p1 = max( 0, min(jpiwob+1 - nimpp+1, jpi - 1 ) )
(  354)         njw0   = max( 2, min(jpjwd    - njmpp+1, jpj     ) )
(  355)         njw1   = max( 0, min(jpjwf    - njmpp+1, jpj - 1 ) )
(  356)         njw0p1 = max( 1, min(jpjwdp1  - njmpp+1, jpj     ) )
(  357)         njw0m1 = max( 1, min(jpjwd    - njmpp+1, jpj     ) )
(  358)         njw1m1 = max( 0, min(jpjwfm1  - njmpp+1, jpj - 1 ) )
(  359)         njw1m2 = max( 0, min(jpjwfm1-1- njmpp+1, jpj - 1 ) )
(  360)         IF(lwp) THEN
(  361)            IF( lfbcwest ) THEN
(  362)               WRITE(numout,*)'     '
(  363)               WRITE(numout,*)'         Specified West Open Boundary'
(  364)            ELSE
(  365)               WRITE(numout,*)'     '
(  366)               WRITE(numout,*)'         Radiative West Open Boundary'
(  367)            END IF
(  368)         END IF
(  369) !DB: output above plus nimpp njmpp in a readable form:
(  370)         if(DBG) then
(  371)            write(2200+narea,'(a35,6(i3,1x))')'OBC west variables for narea ',narea
(  372)            write(2200+narea,'(a35,6(i3,1x))')'nimpp njmpp ', nimpp, njmpp
(  373)            write(2200+narea,'(a35,6(i3,1x))')'jpi jpj jpiwob jpjwd jpjwf ',jpi,jpj, jpiwob, jpjwd, jpjwf 
(  374)            write(2200+narea,'(a35,6(i3,1x))')'niw0, niw1, njw0, njw1: ',niw0, niw1, njw0, njw1
(  375)            write(2200+narea,'(a35,6(i3,1x))')'njw0p1, njw1m1, niw0p1, niw1p1: ',njw0p1, njw1m1, niw0p1, niw1p1
(  376)            write(2200+narea,'(a35,6(i3,1x))')'njw0m1, njw1m2, niw0p1, niw1p1: ',njw0m1, njw1m2, niw0p1, niw1p1
(  377)         endif
(  378) 
(  379) !!DBG
(  380)         if(mod(narea,jpni) == 1) then        !!Assume tiling goes from left-to-right btm-to-top
(  381)            ji = (narea-1)/jpni + 1
(  382)            w_bnd1(ji) = njw0p1-1+njmpp
(  383)            w_bnd2(ji) = njw1m1-1+njmpp
(  384)         endif
(  385) 
(  386)      END IF
(  387)      if(lk_mpp) then
(  388)         call mpp_sum(w_bnd1,jpnj)
(  389)         call mpp_sum(w_bnd2,jpnj)
(  390)      endif
(  391) !!DB: Check if any tiles overlap by more than 1 and report
(  392) !!NB: Assumes halo zone = 1
(  393)      if(lwp) then
(  394) !!DBG






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 8

(  395)         write(1600,'(A30,100(i4,1x))')'(WEST) 1st index in tile:  ',(w_bnd1(ji),ji=1,jpnj)
(  396)         write(1600,'(A30,100(i4,1x))')'(WEST) Last index in tile: ',(w_bnd2(ji),ji=1,jpnj)
(  397)         do ji = 2, jpnj
(  398)            if(w_bnd1(ji) < w_bnd2(ji-1)) then
(  399)               write(1600,*)'Problem at West boundary between tiles ',ji, ji-1
(  400)            endif
(  401)         enddo
(  402)      endif
(  403) 
(  404) 
(  405)      IF( lp_obc_north ) THEN
(  406)         ! ...   mpp initialization
(  407)         nin0   = max( 2, min(jpind    - nimpp+1, jpi     ) )
(  408)         nin1   = max( 0, min(jpinf    - nimpp+1, jpi - 1 ) )
(  409)         nin0p1 = max( 1, min(jpindp1  - nimpp+1, jpi     ) )
(  410)         nin0m1 = max( 1, min(jpind    - nimpp+1, jpi     ) )
(  411)         nin1m1 = max( 0, min(jpinfm1  - nimpp+1, jpi - 1 ) )
(  412)         nin1m2 = max( 0, min(jpinfm1-1- nimpp+1, jpi - 1 ) )
(  413)         njn0   = max( 1, min(jpjnob   - njmpp+1, jpj     ) )
(  414)         njn1   = max( 0, min(jpjnob   - njmpp+1, jpj - 1 ) )
(  415)         njn0p1 = max( 1, min(jpjnob+1 - njmpp+1, jpj     ) )
(  416)         njn1p1 = max( 0, min(jpjnob+1 - njmpp+1, jpj - 1 ) )
(  417)         njn0m1 = max( 1, min(jpjnob-1 - njmpp+1, jpj     ) )
(  418)         njn1m1 = max( 0, min(jpjnob-1 - njmpp+1, jpj - 1 ) )
(  419)         IF(lwp) THEN
(  420)            IF( lfbcnorth ) THEN
(  421)               WRITE(numout,*)'     '
(  422)               WRITE(numout,*)'         Specified North Open Boundary'
(  423)            ELSE
(  424)               WRITE(numout,*)'     '
(  425)               WRITE(numout,*)'         Radiative North Open Boundary'
(  426)            END IF
(  427)         END IF
(  428)      END IF
(  429) 
(  430)      IF( lp_obc_south ) THEN
(  431)         ! ...   mpp initialization
(  432)         nis0   = max( 2, min(jpisd    - nimpp+1, jpi     ) )
(  433)         nis1   = max( 0, min(jpisf    - nimpp+1, jpi - 1 ) )
(  434)         nis0p1 = max( 1, min(jpisdp1  - nimpp+1, jpi     ) )
(  435)         nis0m1 = max( 1, min(jpisd    - nimpp+1, jpi     ) )
(  436)         nis1m1 = max( 0, min(jpisfm1  - nimpp+1, jpi - 1 ) )
(  437)         nis1m2 = max( 0, min(jpisfm1-1- nimpp+1, jpi - 1 ) )
(  438)         njs0   = max( 1, min(jpjsob   - njmpp+1, jpj     ) )
(  439)         njs1   = max( 0, min(jpjsob   - njmpp+1, jpj - 1 ) )
(  440)         njs0p1 = max( 1, min(jpjsob+1 - njmpp+1, jpj     ) )
(  441)         njs1p1 = max( 0, min(jpjsob+1 - njmpp+1, jpj - 1 ) )
(  442)         IF(lwp) THEN
(  443)            IF( lfbcsouth ) THEN
(  444)               WRITE(numout,*)'     '
(  445)               WRITE(numout,*)'         Specified South Open Boundary'
(  446)            ELSE
(  447)               WRITE(numout,*)'     '
(  448)               WRITE(numout,*)'         Radiative South Open Boundary'
(  449)            END IF
(  450)         END IF
(  451) !DB: output above plus nimpp njmpp in a readable form:
(  452)         if(DBG) then






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 9

(  453)            write(2300+narea,'(a35,6(i3,1x))')'OBC south variables for narea ',narea
(  454)            write(2300+narea,'(a35,6(i3,1x))')'nimpp njmpp ', nimpp, njmpp
(  455)            write(2300+narea,'(a35,6(i3,1x))')'jpi jpj jpjsob jpisd jpisf ',jpi,jpj, jpjsob, jpisd, jpisf 
(  456)            write(2300+narea,'(a35,6(i3,1x))')'nis0, nis1, njs0, njs1: ',nis0, nis1, njs0, njs1
(  457)            write(2300+narea,'(a35,6(i3,1x))')'nis0p1, nis1m1, nis0m1, nis1m2: ',nis0p1, nis1m1, nis0m1, nis1m2
(  458)            write(2300+narea,'(a35,6(i3,1x))')'njs0p1, njs1p1, njs0p1, njs1p1: ',njs0p1, njs1p1, njs0p1, njs1p1
(  459)         endif
(  460) 
(  461) !!DBG
(  462)         if(narea <= jpni) then        !!Assume tiling goes from left-to-right btm-to-top
(  463)            ji = mod(narea,jpni)
(  464)            if(ji==0) ji=jpni
(  465)            s_bnd1(ji) = nis0p1-1+nimpp
(  466)            s_bnd2(ji) = nis1m1-1+nimpp
(  467)         endif
(  468) 
(  469)      END IF
(  470) 
(  471)      if(lk_mpp) then
(  472)         call mpp_sum(s_bnd1,jpni)
(  473)         call mpp_sum(s_bnd2,jpni)
(  474)      endif
(  475) !!DB: Check if any tiles overlap by more than 1 and report
(  476) !!NB: Assumes halo zone = 1
(  477)      if(lwp) then
(  478) 
(  479) !!DBG
(  480)         write(1700,'(A30,100(i4,1x))')'(SOUTH) 1st index in tile:  ',(s_bnd1(ji),ji=1,jpni)
(  481)         write(1700,'(A30,100(i4,1x))')'(SOUTH) Last index in tile: ',(s_bnd2(ji),ji=1,jpni)
(  482)         do ji = 2, jpni
(  483)            if(s_bnd1(ji) < s_bnd2(ji-1)) then
(  484)               write(1700,*)'Problem at South boundary between tiles ',ji, ji-1
(  485)            endif
(  486)         enddo
(  487)      endif
(  488) 
(  489) 
(  490) 
(  491)      ! 3. mask correction for OBCs
(  492)      ! ---------------------------
(  493) 
(  494)      IF( lp_obc_east ) THEN
(  495)         !... (jpjed,jpjefm1),jpieob
(  496)         DO jj = nje0, nje1m1
(  497) # if defined key_dynspg_rl
(  498)            DO ji = nie0, nie1
(  499) # else
(  500)            DO ji = nie0p1, nie1p1
(  501) # endif
(  502)               bmask(ji,jj) = 0.e0
(  503)            END DO
(  504)         END DO
(  505) 
(  506)         ! ... initilization to zero
(  507)         uemsk(:,:) = 0.e0
(  508)         vemsk(:,:) = 0.e0
(  509)         temsk(:,:) = 0.e0
(  510) !byoung






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 10

(  511)         temsk5(:,:,:) = 0.e0
(  512)         ! ... set 2D mask on East OBC,  Vopt
(  513)         DO ji = fs_nie0, fs_nie1
(  514)            DO jj = nje0, nje1
(  515)               uemsk(jj,:) = umask(ji,  jj,:)
(  516)               vemsk(jj,:) = vmask(ji+1,jj,:)
(  517)               temsk(jj,:) = tmask(ji+1,jj,:)
(  518)               !byoung
(  519)               DO j1 = 1,5
(  520)                  temsk5(jj,:,j1) = tmask(ji+2-j1,jj,:)
(  521)               END DO
(  522)            END DO
(  523)         END DO
(  524) 
(  525) !!DB: create global east mask
(  526)         do ij = fs_nie0, fs_nie1
(  527)            do jj = nje0, nje1
(  528)               ii = jj-1+njmpp
(  529)               emaskg(ii) = uemsk(jj,1)
(  530)               emaskg2(ii,:) = uemsk(jj,:)
(  531)            enddo
(  532)         enddo
(  533)         call mpp_sum(emaskg,jpjdta)
(  534) !!DB: correct multiple overlaps 
(  535)         do ii = 1, jpjdta
(  536)            if(emaskg(ii) > 1) emaskg(ii) = 1
(  537)         enddo
(  538) !!DB create global mask using mpp_sum (Check=OK)
(  539)         do ii = 1, jpjdta
(  540)            call mpp_sum(emaskg2(ii,:),jpk)
(  541)         enddo
(  542) !!DB: correct multiple overlaps 
(  543)         do ii = 1, jpjdta
(  544)            do jj = 1, jpk
(  545)               if(emaskg2(ii,jj) > 1)emaskg2(ii,jj) = 1
(  546)            enddo
(  547)         enddo
(  548) !!DBG: output global mask. Checked that all CPUs are the same so only 1 output required
(  549) !        if(lwp) then
(  550) !           do ii = 1, jpjdta
(  551) !              write(6500,'(i3,3x,200(f3.0,1x))')ii, (emaskg2(ii,jj),jj=1,jpk)
(  552) !           enddo
(  553) !        endif
(  554)         
(  555)      END IF
(  556) 
(  557)      IF( lp_obc_west ) THEN
(  558)         ! ... (jpjwd,jpjwfm1),jpiwob
(  559)         DO jj = njw0, njw1m1
(  560)            DO ji = niw0, niw1
(  561)               bmask(ji,jj) = 0.e0
(  562)            END DO
(  563)         END DO
(  564) 
(  565)         ! ... initilization to zero
(  566)         uwmsk(:,:) = 0.e0
(  567)         vwmsk(:,:) = 0.e0
(  568)         twmsk(:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 11

(  569) !byoung
(  570)         twmsk5(:,:,:) = 0.e0
(  571) 
(  572)         ! ... set 2D mask on West OBC,  Vopt
(  573)         DO ji = fs_niw0, fs_niw1
(  574)            DO jj = njw0, njw1
(  575)               uwmsk(jj,:) = umask(ji,jj,:)
(  576)               vwmsk(jj,:) = vmask(ji,jj,:)
(  577)               twmsk(jj,:) = tmask(ji,jj,:)
(  578)               !byoung
(  579)               DO j1 = 1,5
(  580)                  twmsk5(jj,:,j1) = tmask(ji+j1-1,jj,:)
(  581)               END DO
(  582)            END DO
(  583)         END DO
(  584) !!DB: create global west mask
(  585)         do ij = fs_niw0, fs_niw1
(  586)            do jj = njw0, njw1
(  587)               ii = jj-1+njmpp
(  588)               wmaskg(ii) = uwmsk(jj,1)
(  589)               wmaskg2(ii,:) = uwmsk(jj,:)
(  590)            enddo
(  591)         enddo
(  592)         call mpp_sum(wmaskg,jpjdta)
(  593) !DB: correct multiple overlaps here
(  594)         do ii = 1, jpjdta
(  595)            if(wmaskg(ii) > 1) wmaskg(ii) = 1
(  596)         enddo
(  597) !!DB create global mask using mpp_sum (Check=OK)
(  598)         do ii = 1, jpjdta
(  599)            call mpp_sum(wmaskg2(ii,:),jpk)
(  600)         enddo
(  601) !!DB: correct multiple overlaps 
(  602)         do ii = 1, jpjdta
(  603)            do jj = 1, jpk
(  604)               if(wmaskg2(ii,jj) > 1)wmaskg2(ii,jj) = 1
(  605)            enddo
(  606)         enddo
(  607) !!DBG: output global mask. Checked that all CPUs are the same so only 1 output required
(  608) !        if(lwp) then
(  609) !           do ii = 1, jpjdta
(  610) !              write(6600,'(i3,3x,200(f3.0,1x))')ii, (wmaskg2(ii,jj),jj=1,jpk)
(  611) !           enddo
(  612) !        endif
(  613) 
(  614)      END IF
(  615) 
(  616)      IF( lp_obc_north ) THEN
(  617)         ! ... jpjnob,(jpind,jpisfm1)
(  618) # if defined key_dynspg_rl
(  619)         DO jj = njn0, njn1
(  620) # else
(  621)         DO jj = njn0p1, njn1p1
(  622) # endif
(  623)            DO ji = nin0, nin1m1
(  624)               bmask(ji,jj) = 0.e0
(  625)            END DO
(  626)         END DO






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 12

(  627) 
(  628)         ! ... initilization to zero
(  629)         unmsk(:,:) = 0.e0
(  630)         vnmsk(:,:) = 0.e0
(  631)         tnmsk(:,:) = 0.e0
(  632) !sujie
(  633)         tnmsk5(:,:,:) = 0.e0
(  634) 
(  635)         ! ... set 2D mask on North OBC,  Vopt
(  636)         DO jj = fs_njn0, fs_njn1
(  637)            DO ji = nin0, nin1
(  638)               unmsk(ji,:) = umask(ji,jj+1,:)
(  639)               vnmsk(ji,:) = vmask(ji,jj  ,:)
(  640)               tnmsk(ji,:) = tmask(ji,jj+1,:)
(  641) !sujie
(  642)            DO j1 = 1,5
(  643)               tnmsk5(ji,:,j1) = tmask(ji,jj+2-j1,:)
(  644)            END DO
(  645)            END DO
(  646)         END DO
(  647) !!DB: create global north mask
(  648)         do ij = fs_njn0, fs_njn1
(  649)            do jj = nin0, nin1
(  650)               ii = jj-1+nimpp
(  651)               nmaskg(ii) = vnmsk(jj,1)
(  652)            enddo
(  653)         enddo
(  654)         call mpp_sum(nmaskg,jpidta)
(  655) !!DB: correct multiple overlaps here
(  656)         do ii = 1, jpidta
(  657)            if(nmaskg(ii) > 1) nmaskg(ii) = 1
(  658)         enddo
(  659) 
(  660)      END IF
(  661) 
(  662)      IF( lp_obc_south ) THEN
(  663)         ! ... jpjsob,(jpisd,jpisfm1)
(  664)         DO jj = njs0, njs1
(  665)            DO ji = nis0, nis1m1
(  666)               bmask(ji,jj) = 0.e0
(  667)            END DO
(  668)         END DO
(  669) 
(  670)         ! ... initilization to zero
(  671)         usmsk(:,:) = 0.e0
(  672)         vsmsk(:,:) = 0.e0
(  673)         tsmsk(:,:) = 0.e0
(  674) !sujie
(  675)         tsmsk5(:,:,:) = 0.e0
(  676) 
(  677)         ! ... set 2D mask on South OBC,  Vopt
(  678) !        DO jj = njs0, njs1
(  679)         DO jj = fs_njs0, fs_njs1
(  680)            DO ji = nis0, nis1
(  681)               usmsk(ji,:) = umask(ji,jj,:)
(  682)               vsmsk(ji,:) = vmask(ji,jj,:)
(  683)               tsmsk(ji,:) = tmask(ji,jj,:)
(  684)               !sujie






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 13

(  685)               DO j1 = 1,5
(  686)                  tsmsk5(ji,:,j1) = tmask(ji,jj+j1-1,:)
(  687)               END DO
(  688)            END DO
(  689)         END DO
(  690) !!DB: create global south mask
(  691)         do ij = fs_njs0, fs_njs1
(  692)            do jj = nis0, nis1
(  693)               ii = jj-1+nimpp
(  694)               smaskg(ii) = vsmsk(jj,1)
(  695)               smaskg2(ii,:) = vsmsk(jj,:)
(  696)            enddo
(  697)         enddo
(  698)         call mpp_sum(smaskg,jpidta)
(  699) !!DB: correct multiple overlaps here
(  700)         do ii = 1, jpidta
(  701)            if(smaskg(ii) > 1) smaskg(ii) = 1
(  702)         enddo
(  703) !!DB create global mask using mpp_sum (Check=OK)
(  704)         do ii = 1, jpidta
(  705)            call mpp_sum(smaskg2(ii,:),jpk)
(  706)         enddo
(  707) !!DB: correct multiple overlaps 
(  708)         do ii = 1, jpidta
(  709)            do jj = 1, jpk
(  710)               if(smaskg2(ii,jj) > 1)smaskg2(ii,jj) = 1
(  711)            enddo
(  712)         enddo
(  713) !!DBG: output global mask. Checked that all CPUs are the same so only 1 output required
(  714) !        if(lwp) then
(  715) !           do ii = 1, jpidta
(  716) !              write(6700,'(i3,3x,200(f3.0,1x))')ii, (smaskg2(ii,jj),jj=1,jpk)
(  717) !           enddo
(  718) !        endif
(  719) 
(  720)      END IF
(  721) 
(  722) # if defined key_dynspg_flt
(  723) 
(  724)      ! ... Initialize obcumask and obcvmask for the Force filtering
(  725)      !     boundary condition in dynspg_flt
(  726)      obcumask(:,:) = umask(:,:,1)
(  727)      obcvmask(:,:) = vmask(:,:,1)
(  728) 
(  729)      ! ... Initialize obctmsk on overlap region and obcs. This mask
(  730)      !     is used in obcvol.F90 to calculate cumulate flux E-P.
(  731)      !     - no flux E-P on obcs and overlap region (jpereci = jprecj = 1)
(  732)      obctmsk(:,:) = tmask(:,:,1)
(  733)      obctmsk(1  ,:) = 0.e0
(  734)      obctmsk(jpi,:) = 0.e0
(  735)      obctmsk(:  ,1) = 0.e0
(  736)      obctmsk(:,jpj) = 0.e0
(  737) 
(  738)      IF( lp_obc_east ) THEN
(  739)         ! ... East obc Force filtering mask for the grad D
(  740)         DO ji = nie0, nie1
(  741)            DO jj = nje0p1, nje1m1
(  742)               obcumask(ji  ,jj)=0.e0






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 14

(  743)               obcvmask(ji+1,jj)=0.e0
(  744)            END DO
(  745)         END DO
(  746) 
(  747)         ! ... set to 0 on East OBC
(  748)         DO jj = nje0p1, nje1m1
(  749)            DO ji = nie0p1, nie1p1
(  750)               obctmsk(ji,jj) = 0.e0
(  751)            END DO
(  752)         END DO
(  753)      END IF
(  754) 
(  755)      IF( lp_obc_west ) THEN
(  756)         ! ... West obc Force filtering mask for the grad D
(  757)         DO ji = niw0, niw1
(  758)            DO jj = njw0p1, njw1m1
(  759)               obcumask(ji,jj)=0.e0
(  760)               obcvmask(ji,jj)=0.e0
(  761)            END DO
(  762)         END DO
(  763) 
(  764)         ! ... set to 0 on West OBC
(  765)         DO jj = njw0p1, njw1m1
(  766)            DO ji = niw0, niw1
(  767)               obctmsk(ji,jj) = 0.e0
(  768)            END DO
(  769)         END DO
(  770)      END IF
(  771) 
(  772)      IF( lp_obc_north ) THEN
(  773)         ! ... North obc Force filtering mask for the grad D
(  774)         DO jj = njn0, njn1
(  775)            DO ji = nin0p1, nin1m1
(  776)               obcvmask(ji,jj  )=0.e0
(  777)               obcumask(ji,jj+1)=0.e0
(  778)            END DO
(  779)         END DO
(  780) 
(  781)         ! ... set to 0 on North OBC
(  782)         DO jj = njn0p1, njn1p1
(  783)            DO ji = nin0p1, nin1m1
(  784)               obctmsk(ji,jj) = 0.e0
(  785)            END DO
(  786)         END DO
(  787)      END IF
(  788) 
(  789)      IF( lp_obc_south ) THEN
(  790)         ! ... South obc Force filtering mask for the grad D
(  791)         DO jj = njs0, njs1
(  792)            DO ji = nis0p1, nis1m1
(  793)               obcumask(ji,jj)=0.e0
(  794)               obcvmask(ji,jj)=0.e0
(  795)            END DO
(  796)         END DO
(  797) 
(  798)         ! ... set to 0 on South OBC
(  799)         DO jj = njs0, njs1
(  800)            DO ji = nis0p1, nis1m1






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 15

(  801)               obctmsk(ji,jj) = 0.e0
(  802)            END DO
(  803)         END DO
(  804)      END IF
(  805) 
(  806) # endif
(  807) 
(  808) # if ! defined key_dynspg_rl
(  809) 
(  810)      IF ( ln_vol_cst ) THEN
(  811) 
(  812)         ! 3.1 Total lateral surface for each open boundary
(  813)         ! ------------------------------------------------
(  814) 
(  815)         ! ... West open boundary surface
(  816)         IF( lp_obc_west ) THEN
(  817)            DO ji = niw0, niw1
(  818)               DO jj = 1, jpj
(  819)                  obcsurftot = obcsurftot+hu(ji,jj)*e2u(ji,jj)*uwmsk(jj,1)
(  820)               END DO
(  821)            END DO
(  822)         END IF
(  823) 
(  824)         ! ... East open boundary surface
(  825)         IF( lp_obc_east ) THEN
(  826)            DO ji = nie0, nie1
(  827)               DO jj = 1, jpj
(  828)                  obcsurftot = obcsurftot+hu(ji,jj)*e2u(ji,jj)*uemsk(jj,1)
(  829)               END DO
(  830)            END DO
(  831)         END IF
(  832) 
(  833)         ! ... North open boundary vertical surface
(  834)         IF( lp_obc_north ) THEN
(  835)            DO jj = njn0, njn1
(  836)               DO ji = 1, jpi
(  837)                  obcsurftot = obcsurftot+hv(ji,jj)*e1v(ji,jj)*vnmsk(ji,1)
(  838)               END DO
(  839)            END DO
(  840)         END IF
(  841) 
(  842)         ! ... South open boundary vertical surface
(  843)         IF( lp_obc_south ) THEN
(  844)            DO jj = njs0, njs1
(  845)               DO ji = 1, jpi
(  846)                  obcsurftot = obcsurftot+hv(ji,jj)*e1v(ji,jj)*vsmsk(ji,1)
(  847)               END DO
(  848)            END DO
(  849)         END IF
(  850)         IF( lk_mpp )   CALL mpp_sum( obcsurftot )   ! sum over the global domain
(  851)      ENDIF
(  852) # endif
(  853) 
(  854)      ! 5. Control print on mask
(  855)      !    The extremities of the open boundaries must be in land
(  856)      !    or else correspond to an "ocean corner" between two open boundaries.
(  857)      !    corner 1 is southwest, 2 is south east, 3 is northeast, 4 is northwest.
(  858)      ! --------------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 16

(  859) 
(  860)      icorner(:)=0
(  861) 
(  862)      ! ... control of the west boundary
(  863)      IF( lp_obc_west ) THEN
(  864)         IF( jpiwob < 2 .OR.  jpiwob >= jpiglo-2 ) THEN
(  865)            IF(lwp) WRITE(numout,*)
(  866)            IF(lwp) WRITE(numout,*) ' E R R O R : jpiwob exceed ', jpiglo-2, 'or less than 2'
(  867)            IF(lwp) WRITE(numout,*) ' ========== '
(  868)            IF(lwp) WRITE(numout,*)
(  869)            nstop = nstop + 1
(  870)         END IF
(  871)         ztestmask(:)=0.
(  872)         DO ji=niw0,niw1
(  873)            IF( (njw0 + njmpp - 1) == jpjwd ) ztestmask(1)=ztestmask(1)+ tmask(ji,njw0,1)
(  874)            IF( (njw1 + njmpp - 1) == jpjwf ) ztestmask(2)=ztestmask(2)+ tmask(ji,njw1,1)
(  875)         END DO
(  876)         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain
(  877) 
(  878)         IF( ztestmask(1) /= 0. ) icorner(1)=icorner(1)+1
(  879)         IF( ztestmask(2) /= 0. ) icorner(4)=icorner(4)+1
(  880)      END IF
(  881) 
(  882)      ! ... control of the east boundary
(  883)      IF( lp_obc_east ) THEN
(  884)         IF( jpieob < 4 .OR.  jpieob >= jpiglo ) THEN
(  885)            IF(lwp) WRITE(numout,*)
(  886)            IF(lwp) WRITE(numout,*) ' E R R O R : jpieob exceed ', jpiglo, ' or less than 4'
(  887)            IF(lwp) WRITE(numout,*) ' ========== '
(  888)            IF(lwp) WRITE(numout,*)
(  889)            nstop = nstop + 1
(  890)         END IF
(  891)         ztestmask(:)=0.
(  892)         DO ji=nie0p1,nie1p1
(  893)            IF( (nje0 + njmpp - 1) == jpjed ) ztestmask(1)=ztestmask(1)+ tmask(ji,nje0,1)
(  894)            IF( (nje1 + njmpp - 1) == jpjef ) ztestmask(2)=ztestmask(2)+ tmask(ji,nje1,1)
(  895)         END DO
(  896)         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain
(  897) 
(  898)        IF( ztestmask(1) /= 0. ) icorner(2)=icorner(2)+1
(  899)        IF( ztestmask(2) /= 0. ) icorner(3)=icorner(3)+1
(  900)      END IF
(  901) 
(  902)      ! ... control of the north boundary
(  903)      IF( lp_obc_north ) THEN
(  904)         IF( jpjnob < 4 .OR.  jpjnob >= jpjglo ) THEN
(  905)            IF(lwp) WRITE(numout,*)
(  906)            IF(lwp) WRITE(numout,*) ' E R R O R : jpjnob exceed ', jpjglo, ' or less than 4'
(  907)            IF(lwp) WRITE(numout,*) ' ========== '
(  908)            IF(lwp) WRITE(numout,*)
(  909)            nstop = nstop + 1
(  910)         END IF
(  911)         ztestmask(:)=0.
(  912)         DO jj=njn0p1,njn1p1
(  913)            IF( (nin0 + nimpp - 1) == jpind ) ztestmask(1)=ztestmask(1)+ tmask(nin0,jj,1)
(  914)            IF( (nin1 + nimpp - 1) == jpinf ) ztestmask(2)=ztestmask(2)+ tmask(nin1,jj,1)
(  915)         END DO
(  916)         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 17

(  917) 
(  918)         IF( ztestmask(1) /= 0. ) icorner(4)=icorner(4)+1
(  919)         IF( ztestmask(2) /= 0. ) icorner(3)=icorner(3)+1
(  920)      END IF
(  921) 
(  922)      ! ... control of the south boundary
(  923)      IF( lp_obc_south ) THEN
(  924)         IF( jpjsob < 2 .OR.  jpjsob >= jpjglo-2 ) THEN
(  925)            IF(lwp) WRITE(numout,*)
(  926)            IF(lwp) WRITE(numout,*) ' E R R O R : jpjsob exceed ', jpjglo-2, ' or less than 2'
(  927)            IF(lwp) WRITE(numout,*) ' ========== '
(  928)            IF(lwp) WRITE(numout,*)
(  929)            nstop = nstop + 1
(  930)         END IF
(  931)         ztestmask(:)=0.
(  932)         DO jj=njs0,njs1
(  933)            IF( (nis0 + nimpp - 1) == jpisd ) ztestmask(1)=ztestmask(1)+ tmask(nis0,jj,1)
(  934)            IF( (nis1 + nimpp - 1) == jpisf ) ztestmask(2)=ztestmask(2)+ tmask(nis1,jj,1)
(  935)         END DO
(  936)         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain
(  937) 
(  938)         IF( ztestmask(1) /= 0. ) icorner(1)=icorner(1)+1
(  939)         IF( ztestmask(2) /= 0. ) icorner(2)=icorner(2)+1
(  940)      END IF
(  941) 
(  942)      IF( icorner(1) == 2 ) THEN
(  943)         IF(lwp) WRITE(numout,*)
(  944)         IF(lwp) WRITE(numout,*) ' South West ocean corner, two open boudaries'
(  945)         IF(lwp) WRITE(numout,*) ' ========== '
(  946)         IF(lwp) WRITE(numout,*)
(  947)         IF( jpisd /= jpiwob.OR.jpjsob /= jpjwd ) THEN
(  948)            IF(lwp) WRITE(numout,*) ' Open boundaries do not fit, we stop'
(  949)            nstop = nstop + 1
(  950)         END IF
(  951)      ELSE IF( icorner(1) == 1 ) THEN
(  952)              IF(lwp) WRITE(numout,*) ' Open boundaries do not fit at SW corner, we stop'
(  953)              nstop = nstop + 1
(  954)      END IF
(  955) 
(  956)      IF( icorner(2) == 2 ) THEN
(  957)          IF(lwp) WRITE(numout,*)
(  958)          IF(lwp) WRITE(numout,*) ' South East ocean corner, two open boudaries'
(  959)          IF(lwp) WRITE(numout,*) ' ========== '
(  960)          IF(lwp) WRITE(numout,*)
(  961)          IF( jpisf /= jpieob+1.OR.jpjsob /= jpjed ) THEN
(  962)             IF(lwp) WRITE(numout,*) ' Open boundaries do not fit, we stop'
(  963)             nstop = nstop + 1
(  964)          END IF
(  965)      ELSE IF( icorner(2) == 1 ) THEN
(  966)              IF(lwp) WRITE(numout,*) ' Open boundaries do not fit at SE corner, we stop'
(  967)              nstop = nstop + 1
(  968)      END IF
(  969) 
(  970)      IF( icorner(3) == 2 ) THEN
(  971)         IF(lwp) WRITE(numout,*)
(  972)         IF(lwp) WRITE(numout,*) ' North East ocean corner, two open boudaries'
(  973)         IF(lwp) WRITE(numout,*) ' ========== '
(  974)         IF(lwp) WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 18

(  975)         IF( jpinf /= jpieob+1 .OR. jpjnob+1 /= jpjef ) THEN
(  976)            IF(lwp) WRITE(numout,*) ' Open boundaries do not fit, we stop'
(  977)            nstop = nstop + 1
(  978)         END IF
(  979)       ELSE IF( icorner(3) == 1 ) THEN
(  980)               IF(lwp) WRITE(numout,*) ' Open boundaries do not fit at NE corner, we stop'
(  981)               nstop = nstop + 1
(  982)       END IF
(  983) 
(  984)      IF( icorner(4) == 2 ) THEN
(  985)         IF(lwp) WRITE(numout,*)
(  986)         IF(lwp) WRITE(numout,*) ' North West ocean corner, two open boudaries'
(  987)         IF(lwp) WRITE(numout,*) ' ========== '
(  988)         IF(lwp) WRITE(numout,*)
(  989)         IF( jpind /= jpiwob.OR.jpjnob+1 /= jpjwf ) THEN
(  990)            IF(lwp) WRITE(numout,*) ' Open boundaries do not fit, we stop'
(  991)            nstop = nstop + 1
(  992)         END IF
(  993)       ELSE IF( icorner(4) == 1 ) THEN
(  994)               IF(lwp) WRITE(numout,*) ' Open boundaries do not fit at NW corner, we stop'
(  995)               nstop = nstop + 1
(  996)       END IF
(  997) 
(  998)      ! 6. Initialization of open boundary variables (u, v, bsf, t, s)
(  999)      ! --------------------------------------------------------------
( 1000)      !   only if at least one boundary is  radiative
( 1001) 
( 1002)      ! ... Restart from restart.obc
( 1003)      IF ( inumfbc < nbobc .AND.  ln_rstart ) THEN
( 1004)         CALL obc_rst_lec
( 1005)      ELSE
( 1006) 
( 1007)          ! ... Initialization to zero of radiation arrays.
( 1008)          !     Those have dimensions of local subdomains
( 1009) 
( 1010)          bebnd(:,:,:)   = 0.e0   ;   bnbnd(:,:,:)   = 0.e0
( 1011)          uebnd(:,:,:,:) = 0.e0   ;   unbnd(:,:,:,:) = 0.e0
( 1012)          vebnd(:,:,:,:) = 0.e0   ;   vnbnd(:,:,:,:) = 0.e0
( 1013)          tebnd(:,:,:,:) = 0.e0   ;   tnbnd(:,:,:,:) = 0.e0
( 1014)          sebnd(:,:,:,:) = 0.e0   ;   snbnd(:,:,:,:) = 0.e0
( 1015) 
( 1016)          bwbnd(:,:,:)   = 0.e0   ;   bsbnd(:,:,:)   = 0.e0
( 1017)          uwbnd(:,:,:,:) = 0.e0   ;   usbnd(:,:,:,:) = 0.e0
( 1018)          vwbnd(:,:,:,:) = 0.e0   ;   vsbnd(:,:,:,:) = 0.e0
( 1019)          twbnd(:,:,:,:) = 0.e0   ;   tsbnd(:,:,:,:) = 0.e0
( 1020)          swbnd(:,:,:,:) = 0.e0   ;   ssbnd(:,:,:,:) = 0.e0
( 1021) 
( 1022)      END IF
( 1023) 
( 1024) # if defined key_dynspg_rl
( 1025)      ! 7. Isolated coastline arrays initialization (rigid lid case only)
( 1026)      ! -----------------------------------------------------------------
( 1027)      CALL obc_dom
( 1028) # endif
( 1029) 
( 1030)      ! 8. Control print
( 1031)      ! ... control of the east boundary
( 1032)      IF( lp_obc_east ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 19

( 1033)         istop = 0
( 1034)         IF( jpieob < 4 .OR.  jpieob >= jpiglo ) THEN
( 1035)            IF(lwp) WRITE(numout,cform_err)
( 1036)            IF(lwp) WRITE(numout,*) '            jpieob exceed ', jpim1, ' or less than 4'
( 1037)            istop = istop + 1
( 1038)         END IF
( 1039) 
( 1040)         IF( lk_mpp ) THEN
( 1041)            ! ...
( 1042)            IF( nimpp > jpieob-5) THEN
( 1043)               IF(lwp) WRITE(numout,cform_err)
( 1044)               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the East OBC'
( 1045)               IF(lwp) WRITE(numout,*) '        nimpp must be < jpieob-5'
( 1046)               istop = istop + 1
( 1047)            ENDIF
( 1048)         ELSE
( 1049) 
( 1050)            ! ... stop if  e r r o r (s)   detected
( 1051)            IF( istop /= 0 ) THEN
( 1052)               IF(lwp)WRITE(numout,*)
( 1053)               IF(lwp)WRITE(numout,*) istop,' E R R O R (S) detected : stop'
( 1054)               IF(lwp)WRITE(numout,*) ' =============== '
( 1055)               IF(lwp)WRITE(numout,*)
( 1056)               nstop = nstop + 1
( 1057)            ENDIF
( 1058)         ENDIF
( 1059)      ENDIF
( 1060) 
( 1061)      ! ... control of the west boundary
( 1062)      IF( lp_obc_west ) THEN
( 1063)         istop = 0
( 1064)         IF( jpiwob < 2 .OR.  jpiwob >= jpiglo ) THEN
( 1065)            IF(lwp) WRITE(numout,cform_err)
( 1066)            IF(lwp) WRITE(numout,*) '            jpiwob exceed ', jpim1, ' or less than 2'
( 1067)            istop = istop + 1
( 1068)         END IF
( 1069) 
( 1070)         IF( lk_mpp ) THEN
( 1071)            IF( (nimpp < jpiwob+5) .AND. (nimpp > 1) ) THEN
( 1072)               IF(lwp) WRITE(numout,cform_err)
( 1073)               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the West OBC'
( 1074)               IF(lwp) WRITE(numout,*) '        nimpp must be > jpiwob-5 or =1'
( 1075)               istop = istop + 1
( 1076)            ENDIF
( 1077)         ELSE
( 1078) 
( 1079)            ! ... stop if  e r r o r (s)   detected
( 1080)            IF( istop /= 0 ) THEN
( 1081)               IF(lwp)WRITE(numout,*)
( 1082)               IF(lwp)WRITE(numout,*) istop,' E R R O R (S) detected : stop'
( 1083)               IF(lwp)WRITE(numout,*) ' =============== '
( 1084)               IF(lwp)WRITE(numout,*)
( 1085)               nstop = nstop + 1
( 1086)            ENDIF
( 1087)         ENDIF
( 1088)      ENDIF
( 1089) 
( 1090)      ! control of the north boundary






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 20

( 1091)      IF( lp_obc_north ) THEN
( 1092)         istop = 0
( 1093)         IF( jpjnob < 4 .OR.  jpjnob >= jpjglo ) THEN
( 1094)            IF(lwp) WRITE(numout,cform_err)
( 1095)            IF(lwp) WRITE(numout,*) '          jpjnob exceed ', jpjm1,' or less than 4'
( 1096)            istop = istop + 1
( 1097)         END IF
( 1098) 
( 1099)         IF( lk_mpp ) THEN
( 1100)            IF( njmpp > jpjnob-5) THEN
( 1101)               IF(lwp) WRITE(numout,cform_err)
( 1102)               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the North OBC'
( 1103)               IF(lwp) WRITE(numout,*) '        njmpp must be < jpjnob-5'
( 1104)               istop = istop + 1
( 1105)            ENDIF
( 1106)         ELSE
( 1107) 
( 1108)            ! ... stop if  e r r o r (s)   detected
( 1109)            IF( istop /= 0 ) THEN
( 1110)               IF(lwp)WRITE(numout,*)
( 1111)               IF(lwp)WRITE(numout,*) istop,' E R R O R (S) detected : stop'
( 1112)               IF(lwp)WRITE(numout,*) ' =============== '
( 1113)               IF(lwp)WRITE(numout,*)
( 1114)               nstop = nstop + 1
( 1115)            ENDIF
( 1116)         ENDIF
( 1117)      ENDIF
( 1118) 
( 1119)      ! control of the south boundary
( 1120)      IF( lp_obc_south ) THEN
( 1121)         istop = 0
( 1122)         IF( jpjsob < 2 .OR. jpjsob >= jpjglo ) THEN
( 1123)            IF(lwp) WRITE(numout,cform_err)
( 1124)            IF(lwp) WRITE(numout,*) '          jpjsob exceed ', jpjm1,' or less than 2'
( 1125)            istop = istop + 1
( 1126)         END IF
( 1127) 
( 1128)         IF( lk_mpp ) THEN
( 1129)            IF( (njmpp < jpjsob+5) .AND. (njmpp > 1) ) THEN
( 1130)               IF(lwp) WRITE(numout,cform_err)
( 1131)               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the South OBC'
( 1132)               IF(lwp) WRITE(numout,*) '        njmpp must be > jpjsob+5 or =1'
( 1133)               istop = istop + 1
( 1134)            ENDIF
( 1135)         ELSE
( 1136) 
( 1137)            ! ... stop if  e r r o r (s)   detected
( 1138)            IF( istop /= 0 ) THEN
( 1139)               IF(lwp)WRITE(numout,*)
( 1140)               IF(lwp)WRITE(numout,*) istop,' E R R O R (S) detected : stop'
( 1141)               IF(lwp)WRITE(numout,*) ' =============== '
( 1142)               IF(lwp)WRITE(numout,*)
( 1143)               nstop = nstop + 1
( 1144)            ENDIF
( 1145)         ENDIF
( 1146)      ENDIF
( 1147) 
( 1148)   END SUBROUTINE obc_init






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 21

( 1149) 
( 1150) #else
( 1151)   !!---------------------------------------------------------------------------------
( 1152)   !!   Dummy module                                                NO open boundaries
( 1153)   !!---------------------------------------------------------------------------------
( 1154) CONTAINS
( 1155)   SUBROUTINE obc_init      ! Dummy routine
( 1156)   END SUBROUTINE obc_init
( 1157) #endif
( 1158) 
( 1159)   !!=================================================================================
( 1160) END MODULE obcini




















































PGF90 (Version     10.5)          09/23/2010  09:50:07      page 22

(    1) # 1 "obcini.F90"
(    1) !!DB 2008.05.08
(    2) !!(0) In general search for DB and DBG to see my modifications
(    3) !!(1) Regarding the fact that certain CPU tilings COULD result in overlaps
(    4) !!between tiles of greater than 1 column: I have written and retained code
(    5) !!that flags this occurence, and outputs information. Search for:
(    6) !! write(1X00,...)  where X=5,6,7 == E,W,S ----> fort.1500, fort.1600, fort.1700 files
(    7) !!(2) I create 2D global OBC masking files for use in obcdta. They are called:
(    8) !!emaskg2, wmaskg2, smaskg2 (where e.g.: e=east, g=global dimension, 2=2 dimensional).
(    9) !!They are output to: fort.6X00 where X=5,6,7 == E,W,S.
(   10) !!Note the use of mpp_sum() to ensure that all CPUs share the same global values.
(   11) !!(3) If you want output of a lot of debugging info on OBC variable values for
(   12) !!various CPUs (nareas) set DBG=.true., and search DBG and write(####+narea,...)
(   13) 
(   14) 
(   15) MODULE obcini
(   16)   !!=================================================================================
(   17)   !!                       ***  MODULE  obcini  ***
(   18)   !! OBC initial state :  Open boundary initial state
(   19)   !!=================================================================================
(   20) # 21
(   21)   !!---------------------------------------------------------------------------------
(   22)   !!   'key_obc'                                             Open Boundary Conditions
(   23)   !!---------------------------------------------------------------------------------
(   24)   !!   obc_init       : initialization for the open boundary condition
(   25)   !!---------------------------------------------------------------------------------
(   26)   !! * Modules used
(   27)   USE oce             ! ocean dynamics and tracers variables
(   28)   USE dom_oce         ! ocean space and time domain variables
(   29)   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   30)   USE phycst          ! physical constants
(   31)   USE obc_oce         ! ocean open boundary conditions
(   32)   USE lib_mpp         ! for mpp_sum
(   33)   USE in_out_manager  ! I/O units
(   34) 
(   35)   IMPLICIT NONE
(   36)   PRIVATE
(   37) 
(   38)   !! * Routine accessibility
(   39)   PUBLIC obc_init        ! routine called by opa.F90
(   40) 
(   41)   !! * Substitutions
(   42) # 1 "./obc_vectopt_loop_substitute.h90"
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
(   11)    !!----------------------------------------------------------------------
(   12) # 43 "obcini.F90"
(   43) # 43
(   43)   !!---------------------------------------------------------------------------------
(   44)   !!   OPA 9.0 , LOCEAN-IPSL (2005)






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 23

(   45)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcini.F90,v 1.8 2006/03/21 07:52:26 opalod Exp $
(   46)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   47)   !!---------------------------------------------------------------------------------
(   48) 
(   49) CONTAINS
(   50) 
(   51)   SUBROUTINE obc_init
(   52)      !!----------------------------------------------------------------------
(   53)      !!                 ***  ROUTINE obc_init  ***
(   54)      !!
(   55)      !! ** Purpose :   Initialization of the dynamics and tracer fields at
(   56)      !!      the open boundaries.
(   57)      !!
(   58)      !! ** Method  :   initialization of open boundary variables
(   59)      !!      (u, v, bsf) over 3 time step and 3 rows
(   60)      !!      (t, s)      over 2 time step and 2 rows
(   61)      !!      if ln_rstart = .FALSE. : no restart, fields set to zero
(   62)      !!      if ln_rstart = .TRUE.  : restart, fields are read in a file
(   63)      !!      if rdpxxx = 0 then lfbc is set true for this boundary.
(   64)      !!
(   65)      !! ** Input   :   restart.obc file, restart file for open boundaries
(   66)      !!
(   67)      !! History :
(   68)      !!   8.0  !  97-07  (G. Madec)  Original code
(   69)      !!        !  97-11  (J.M. Molines)
(   70)      !!   8.5  !  02-11  (C. Talandier, A-M. Treguier) Free surface, F90
(   71)      !!   9.0  !  05-11  (V. Garnier) Surface pressure gradient organization
(   72)      !!----------------------------------------------------------------------
(   73)      !! * Modules used
(   74)      USE obcrst,   ONLY :   obc_rst_lec   ! Make obc_rst_lec routine available
(   75)      USE obcdom,   ONLY :   obc_dom       ! Make obc_dom routine available
(   76) 
(   77)      !! * Local declarations
(   78) !sujie      INTEGER  ::   ji, jj, istop , inumfbc
(   79)      INTEGER  ::   ji, jj, istop , inumfbc, j1, ii, ij
(   80)      INTEGER, DIMENSION(4) ::   icorner
(   81)      REAL(wp) ::   zbsic1, zbsic2, zbsic3
(   82)      REAL(wp), DIMENSION(2) ::   ztestmask
(   83) !!DB: 2008.05.01
(   84)      INTEGER :: s_bnd1(jpni) = 0, s_bnd2(jpni) = 0
(   85)      INTEGER :: w_bnd1(jpnj) = 0, w_bnd2(jpnj) = 0
(   86)      INTEGER :: e_bnd1(jpnj) = 0, e_bnd2(jpnj) = 0
(   87)      LOGICAL :: DBG
(   88) 
(   89)      NAMELIST/namobc/ rdpein, rdpwin, rdpnin, rdpsin,   &
(   90)         &             rdpeob, rdpwob, rdpnob, rdpsob,   &
(   91)         &             zbsic1, zbsic2, zbsic3,           &
(   92)         &             nbic, volemp, nobc_dta,           &
(   93)         &             ln_obc_clim, ln_vol_cst, ln_obc_fla
(   94)      !!----------------------------------------------------------------------
(   95) 
(   96) 
(   97) !!DB
(   98) !     DBG = .true.
(   99)      DBG = .false.
(  100) 
(  101) 
(  102)      IF(lwp) WRITE(numout,*)






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 24

(  103)      IF(lwp) WRITE(numout,*) 'obc_init : initialization of open boundaries'
(  104)      IF(lwp) WRITE(numout,*) '~~~~~~~~'
(  105) 
(  106) 
(  107)      ! 0. read namelist parameters
(  108)      ! ---------------------------
(  109)      ! default values already set except:
(  110)      zbsic1 = 0.e0
(  111)      zbsic2 = 0.e0
(  112)      zbsic3 = 0.e0
(  113) 
(  114)      ! Namelist namobc : open boundaries
(  115)      REWIND( numnam )
(  116)      READ  ( numnam, namobc )
(  117) 
(  118)      bsfic0(1) = zbsic1
(  119)      bsfic (2) = zbsic2
(  120)      bsfic (3) = zbsic3
(  121) 
(  122)      ! By security we set rdpxin and rdpxob respectively
(  123)      ! to 1. and 15. if the corresponding OBC is not activated
(  124)      IF( .NOT.lp_obc_east ) THEN
(  125)         rdpein = 1.
(  126)         rdpeob = 15.
(  127)      END IF
(  128)      IF( .NOT.lp_obc_west ) THEN
(  129)         rdpwin = 1.
(  130)         rdpwob = 15.
(  131)      END IF
(  132)      IF( .NOT.lp_obc_north ) THEN
(  133)         rdpnin = 1.
(  134)         rdpnob = 15.
(  135)      END IF
(  136)      IF( .NOT.lp_obc_south ) THEN
(  137)         rdpsin = 1.
(  138)         rdpsob = 15.
(  139)      END IF
(  140) 
(  141)      ! number of open boudaries and open boundary indicators
(  142)      nbobc = 0
(  143)      IF( lp_obc_east  )   nbobc = nbobc + 1
(  144)      IF( lp_obc_west  )   nbobc = nbobc + 1
(  145)      IF( lp_obc_north )   nbobc = nbobc + 1
(  146)      IF( lp_obc_south )   nbobc = nbobc + 1
(  147) 
(  148)      IF(lwp) WRITE(numout,*) '         Number of open boundaries    nbobc = ',nbobc
(  149)      IF(lwp) WRITE(numout,*)
(  150)      IF( nbobc /= 0 .AND. jperio /= 0 ) THEN
(  151)         IF(lwp) WRITE(numout,*)
(  152)         IF(lwp) WRITE(numout,*) ' E R R O R : Cyclic or symmetric,',   &
(  153)            ' and open boundary condition are not compatible'
(  154)         IF(lwp) WRITE(numout,*) ' ========== '
(  155)         IF(lwp) WRITE(numout,*)
(  156)         nstop = nstop + 1
(  157)      END IF
(  158) 
(  159)      ! control prints
(  160)      IF(lwp) WRITE(numout,*) '         namobc'






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 25

(  161)      IF(lwp) WRITE(numout,*) ' '
(  162)      IF(lwp) WRITE(numout,*) '         data in file (=1) or     nobc_dta = ', nobc_dta
(  163)      IF(lwp) WRITE(numout,*) '         initial state used (=0)             '
(  164)      IF(lwp) WRITE(numout,*) '         climatology (true) or not:', ln_obc_clim
(  165)      IF(lwp) WRITE(numout,*) ' '
(  166)      IF(lwp) WRITE(numout,*) '                                 WARNING                     '
(  167)      IF(lwp) WRITE(numout,*) '         Flather"s algorithm is applied with explicit free surface scheme                 '
(  168)      IF(lwp) WRITE(numout,*) '         or with free surface time-splitting scheme          '
(  169)      IF(lwp) WRITE(numout,*) '         Nor radiation neither relaxation is allowed with explicit free surface scheme:   '
(  170)      IF(lwp) WRITE(numout,*) '         Radiation and/or relaxation is allowed with free surface time-splitting scheme '
(  171)      IF(lwp) WRITE(numout,*) '         depending of the choice of rdpXin = rdpXob  = 0. for open boundaries             '
(  172)      IF(lwp) WRITE(numout,*) ' '
(  173)      IF(lwp) WRITE(numout,*) '         For the rigid-lid case or the filtered free surface case,                        '
(  174)      IF(lwp) WRITE(numout,*) '         radiation, relaxation or presciption of data can be applied                      '
(  175)      IF( lwp.AND.lp_obc_east ) THEN
(  176)         WRITE(numout,*) '         East open boundary :'
(  177)         WRITE(numout,*) '              i index                    jpieob = ', jpieob
(  178)         WRITE(numout,*) '              damping time scale (days)  rdpeob = ', rdpeob
(  179)         WRITE(numout,*) '              damping time scale (days)  rdpein = ', rdpein
(  180)      ENDIF
(  181) 
(  182)      IF( lwp.AND.lp_obc_west ) THEN
(  183)         WRITE(numout,*) '         West open boundary :'
(  184)         WRITE(numout,*) '              i index                    jpiwob = ', jpiwob
(  185)         WRITE(numout,*) '              damping time scale (days)  rdpwob = ', rdpwob
(  186)         WRITE(numout,*) '              damping time scale (days)  rdpwin = ', rdpwin
(  187)      ENDIF
(  188) 
(  189)      IF( lwp.AND.lp_obc_north ) THEN
(  190)         WRITE(numout,*) '         North open boundary :'
(  191)         WRITE(numout,*) '               j index                    jpjnob = ', jpjnob
(  192)         WRITE(numout,*) '               damping time scale (days)  rdpnob = ', rdpnob
(  193)         WRITE(numout,*) '               damping time scale (days)  rdpnin = ', rdpnin
(  194)      ENDIF
(  195) 
(  196)      IF( lwp.AND.lp_obc_south ) THEN
(  197)         WRITE(numout,*) '         South open boundary :'
(  198)         WRITE(numout,*) '               j index                    jpjsob = ', jpjsob
(  199)         WRITE(numout,*) '               damping time scale (days)  rdpsob = ', rdpsob
(  200)         WRITE(numout,*) '               damping time scale (days)  rdpsin = ', rdpsin
(  201)         WRITE(numout,*) ' '
(  202)      ENDIF
(  203) 
(  204)      ! 1. Initialisation of constants
(  205)      ! ------------------------------
(  206) 
(  207)      ! ... convert rdp$ob in seconds
(  208)      rdpein = rdpein * rday
(  209)      rdpwin = rdpwin * rday
(  210)      rdpnin = rdpnin * rday
(  211)      rdpsin = rdpsin * rday
(  212)      rdpeob = rdpeob * rday
(  213)      rdpwob = rdpwob * rday
(  214)      rdpnob = rdpnob * rday
(  215)      rdpsob = rdpsob * rday
(  216)      lfbceast  = .FALSE.
(  217)      lfbcwest  = .FALSE.
(  218)      lfbcnorth = .FALSE.






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 26

(  219)      lfbcsouth = .FALSE.
(  220)      inumfbc = 0
(  221)      ! ... look for Fixed Boundaries (rdp = 0 )
(  222)      ! ... When specified, lbcxxx flags are set to TRUE and rdpxxx are set to
(  223)      ! ...  a small arbitrary value, (to avoid division by zero further on).
(  224)      ! ...  rdpxxx is not used anymore.
(  225)      IF( lp_obc_east )  THEN
(  226)         IF( (rdpein+rdpeob) == 0 )  THEN
(  227)            lfbceast = .TRUE.
(  228)            rdpein = 1e-3
(  229)            rdpeob = 1e-3
(  230)            inumfbc = inumfbc+1
(  231)         ELSEIF ( (rdpein*rdpeob) == 0 )  THEN
(  232)            IF(lwp) THEN
(  233)               WRITE(numout,cform_err)
(  234)               WRITE(numout,*) 'obc_init : rdpein & rdpeob must be both zero or non zero'
(  235)               nstop = nstop + 1
(  236)            ENDIF
(  237)         END IF
(  238)      END IF
(  239)      IF( lp_obc_west )  THEN
(  240)         IF( (rdpwin + rdpwob) == 0 )  THEN
(  241)            lfbcwest = .TRUE.
(  242)            rdpwin = 1e-3
(  243)            rdpwob = 1e-3
(  244)            inumfbc = inumfbc+1
(  245)         ELSEIF ( (rdpwin*rdpwob) == 0 )  THEN
(  246)            IF(lwp) THEN
(  247)               WRITE(numout,cform_err)
(  248)               WRITE(numout,*) 'obc_init : rdpwin & rdpwob must be both zero or non zero'
(  249)               nstop = nstop + 1
(  250)            ENDIF
(  251)         END IF
(  252)      END IF
(  253)      IF( lp_obc_north )  THEN
(  254)         IF( (rdpnin + rdpnob) == 0 )  THEN
(  255)            lfbcnorth = .TRUE.
(  256)            rdpnin = 1e-3
(  257)            rdpnob = 1e-3
(  258)            inumfbc = inumfbc+1
(  259)         ELSEIF ( (rdpnin*rdpnob) == 0 )  THEN
(  260)            IF(lwp) THEN
(  261)               WRITE(numout,cform_err)
(  262)               WRITE(numout,*) 'obc_init : rdpnin & rdpnob must be both zero or non zero'
(  263)               nstop = nstop + 1
(  264)            ENDIF
(  265)         END IF
(  266)      END IF
(  267)      IF( lp_obc_south )  THEN
(  268)         IF( (rdpsin + rdpsob) == 0 )  THEN
(  269)            lfbcsouth = .TRUE.
(  270)            rdpsin = 1e-3
(  271)            rdpsob = 1e-3
(  272)            inumfbc = inumfbc+1
(  273)         ELSEIF ( (rdpsin*rdpsob) == 0 )  THEN
(  274)            IF(lwp) THEN
(  275)               WRITE(numout,cform_err)
(  276)               WRITE(numout,*) 'obc_init : rdpsin & rdpsob must be both zero or non zero'






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 27

(  277)               nstop = nstop + 1
(  278)            ENDIF
(  279)         END IF
(  280)      END IF
(  281) 
(  282)      ! 2.  Clever mpp indices for loops on the open boundaries.
(  283)      !     The loops will be performed only on the processors
(  284)      !     that contain a given open boundary.
(  285)      ! --------------------------------------------------------
(  286) 
(  287)      IF( lp_obc_east ) THEN
(  288)         ! ...   mpp initialization
(  289)         nie0   = max( 1, min(jpieob   - nimpp+1, jpi     ) )
(  290)         nie1   = max( 0, min(jpieob   - nimpp+1, jpi - 1 ) )
(  291)         nie0p1 = max( 1, min(jpieob+1 - nimpp+1, jpi     ) )
(  292)         nie1p1 = max( 0, min(jpieob+1 - nimpp+1, jpi - 1 ) )
(  293)         nie0m1 = max( 1, min(jpieob-1 - nimpp+1, jpi     ) )
(  294)         nie1m1 = max( 0, min(jpieob-1 - nimpp+1, jpi - 1 ) )
(  295)         nje0   = max( 2, min(jpjed    - njmpp+1, jpj     ) )
(  296)         nje1   = max( 0, min(jpjef    - njmpp+1, jpj - 1 ) )
(  297)         nje0p1 = max( 1, min(jpjedp1  - njmpp+1, jpj     ) )
(  298)         nje0m1 = max( 1, min(jpjed    - njmpp+1, jpj     ) )
(  299)         nje1m1 = max( 0, min(jpjefm1  - njmpp+1, jpj - 1 ) )
(  300)         nje1m2 = max( 0, min(jpjefm1-1- njmpp+1, jpj - 1 ) )
(  301)         IF(lwp) THEN
(  302)            IF( lfbceast ) THEN
(  303)               WRITE(numout,*)'     '
(  304)               WRITE(numout,*)'         Specified East Open Boundary'
(  305)            ELSE
(  306)               WRITE(numout,*)'     '
(  307)               WRITE(numout,*)'         Radiative East Open Boundary'
(  308)            END IF
(  309)         END IF
(  310) 
(  311) !DB: output above plus nimpp njmpp in a readable form:
(  312)         if(DBG) then
(  313)            write(2100+narea,'(a35,6(i3,1x))')'OBC east variables for narea ',narea
(  314)            write(2100+narea,'(a35,6(i3,1x))')'nimpp njmpp ', nimpp, njmpp
(  315)            write(2100+narea,'(a35,6(i3,1x))')'jpi jpj jpieob jpjed jpjef ',jpi,jpj, jpieob, jpjed, jpjef 
(  316)            write(2100+narea,'(a35,6(i3,1x))')'nie0, nie1, nje0, nje1: ',nie0, nie1, nje0, nje1
(  317)            write(2100+narea,'(a35,6(i3,1x))')'nje0p1, nje1m1, nie0p1, nie1m1: ',nje0p1, nje1m1, nie0p1, nie1m1
(  318)            write(2100+narea,'(a35,6(i3,1x))')'nje0m1, nje1m2, nie0m1, nie1p1: ',nje0m1, nje1m2, nie0m1, nie1p1
(  319)         endif
(  320) 
(  321) !!DBG
(  322)         if(mod(narea,jpni) == 0) then        !!Assume tiling goes from left-to-right btm-to-top
(  323)            ji = (narea-1)/jpni + 1
(  324)            e_bnd1(ji) = nje0p1-1+njmpp
(  325)            e_bnd2(ji) = nje1m1-1+njmpp
(  326)         endif
(  327) 
(  328)      END IF
(  329)      if(lk_mpp) then
(  330)         call mpp_sum(e_bnd1,jpnj)
(  331)         call mpp_sum(e_bnd2,jpnj)
(  332)      endif
(  333) !!DB: Check if any tiles overlap by more than 1 and report
(  334) !!NB: Assumes halo zone = 1






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 28

(  335)      if(lwp) then
(  336) !!DBG
(  337)         write(1500,'(A30,100(i4,1x))')'(EAST) 1st index in tile:  ',(e_bnd1(ji),ji=1,jpnj)
(  338)         write(1500,'(A30,100(i4,1x))')'(EAST) Last index in tile: ',(e_bnd2(ji),ji=1,jpnj)
(  339)         do ji = 2, jpnj
(  340)            if(e_bnd1(ji) < e_bnd2(ji-1)) then
(  341)               write(1500,*)'Problem at East boundary between tiles ',ji, ji-1
(  342)            endif
(  343)         enddo
(  344)      endif
(  345) 
(  346) 
(  347) 
(  348)      IF( lp_obc_west ) THEN
(  349)         ! ...   mpp initialization
(  350)         niw0   = max( 1, min(jpiwob   - nimpp+1, jpi     ) )
(  351)         niw1   = max( 0, min(jpiwob   - nimpp+1, jpi - 1 ) )
(  352)         niw0p1 = max( 1, min(jpiwob+1 - nimpp+1, jpi     ) )
(  353)         niw1p1 = max( 0, min(jpiwob+1 - nimpp+1, jpi - 1 ) )
(  354)         njw0   = max( 2, min(jpjwd    - njmpp+1, jpj     ) )
(  355)         njw1   = max( 0, min(jpjwf    - njmpp+1, jpj - 1 ) )
(  356)         njw0p1 = max( 1, min(jpjwdp1  - njmpp+1, jpj     ) )
(  357)         njw0m1 = max( 1, min(jpjwd    - njmpp+1, jpj     ) )
(  358)         njw1m1 = max( 0, min(jpjwfm1  - njmpp+1, jpj - 1 ) )
(  359)         njw1m2 = max( 0, min(jpjwfm1-1- njmpp+1, jpj - 1 ) )
(  360)         IF(lwp) THEN
(  361)            IF( lfbcwest ) THEN
(  362)               WRITE(numout,*)'     '
(  363)               WRITE(numout,*)'         Specified West Open Boundary'
(  364)            ELSE
(  365)               WRITE(numout,*)'     '
(  366)               WRITE(numout,*)'         Radiative West Open Boundary'
(  367)            END IF
(  368)         END IF
(  369) !DB: output above plus nimpp njmpp in a readable form:
(  370)         if(DBG) then
(  371)            write(2200+narea,'(a35,6(i3,1x))')'OBC west variables for narea ',narea
(  372)            write(2200+narea,'(a35,6(i3,1x))')'nimpp njmpp ', nimpp, njmpp
(  373)            write(2200+narea,'(a35,6(i3,1x))')'jpi jpj jpiwob jpjwd jpjwf ',jpi,jpj, jpiwob, jpjwd, jpjwf 
(  374)            write(2200+narea,'(a35,6(i3,1x))')'niw0, niw1, njw0, njw1: ',niw0, niw1, njw0, njw1
(  375)            write(2200+narea,'(a35,6(i3,1x))')'njw0p1, njw1m1, niw0p1, niw1p1: ',njw0p1, njw1m1, niw0p1, niw1p1
(  376)            write(2200+narea,'(a35,6(i3,1x))')'njw0m1, njw1m2, niw0p1, niw1p1: ',njw0m1, njw1m2, niw0p1, niw1p1
(  377)         endif
(  378) 
(  379) !!DBG
(  380)         if(mod(narea,jpni) == 1) then        !!Assume tiling goes from left-to-right btm-to-top
(  381)            ji = (narea-1)/jpni + 1
(  382)            w_bnd1(ji) = njw0p1-1+njmpp
(  383)            w_bnd2(ji) = njw1m1-1+njmpp
(  384)         endif
(  385) 
(  386)      END IF
(  387)      if(lk_mpp) then
(  388)         call mpp_sum(w_bnd1,jpnj)
(  389)         call mpp_sum(w_bnd2,jpnj)
(  390)      endif
(  391) !!DB: Check if any tiles overlap by more than 1 and report
(  392) !!NB: Assumes halo zone = 1






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 29

(  393)      if(lwp) then
(  394) !!DBG
(  395)         write(1600,'(A30,100(i4,1x))')'(WEST) 1st index in tile:  ',(w_bnd1(ji),ji=1,jpnj)
(  396)         write(1600,'(A30,100(i4,1x))')'(WEST) Last index in tile: ',(w_bnd2(ji),ji=1,jpnj)
(  397)         do ji = 2, jpnj
(  398)            if(w_bnd1(ji) < w_bnd2(ji-1)) then
(  399)               write(1600,*)'Problem at West boundary between tiles ',ji, ji-1
(  400)            endif
(  401)         enddo
(  402)      endif
(  403) 
(  404) 
(  405)      IF( lp_obc_north ) THEN
(  406)         ! ...   mpp initialization
(  407)         nin0   = max( 2, min(jpind    - nimpp+1, jpi     ) )
(  408)         nin1   = max( 0, min(jpinf    - nimpp+1, jpi - 1 ) )
(  409)         nin0p1 = max( 1, min(jpindp1  - nimpp+1, jpi     ) )
(  410)         nin0m1 = max( 1, min(jpind    - nimpp+1, jpi     ) )
(  411)         nin1m1 = max( 0, min(jpinfm1  - nimpp+1, jpi - 1 ) )
(  412)         nin1m2 = max( 0, min(jpinfm1-1- nimpp+1, jpi - 1 ) )
(  413)         njn0   = max( 1, min(jpjnob   - njmpp+1, jpj     ) )
(  414)         njn1   = max( 0, min(jpjnob   - njmpp+1, jpj - 1 ) )
(  415)         njn0p1 = max( 1, min(jpjnob+1 - njmpp+1, jpj     ) )
(  416)         njn1p1 = max( 0, min(jpjnob+1 - njmpp+1, jpj - 1 ) )
(  417)         njn0m1 = max( 1, min(jpjnob-1 - njmpp+1, jpj     ) )
(  418)         njn1m1 = max( 0, min(jpjnob-1 - njmpp+1, jpj - 1 ) )
(  419)         IF(lwp) THEN
(  420)            IF( lfbcnorth ) THEN
(  421)               WRITE(numout,*)'     '
(  422)               WRITE(numout,*)'         Specified North Open Boundary'
(  423)            ELSE
(  424)               WRITE(numout,*)'     '
(  425)               WRITE(numout,*)'         Radiative North Open Boundary'
(  426)            END IF
(  427)         END IF
(  428)      END IF
(  429) 
(  430)      IF( lp_obc_south ) THEN
(  431)         ! ...   mpp initialization
(  432)         nis0   = max( 2, min(jpisd    - nimpp+1, jpi     ) )
(  433)         nis1   = max( 0, min(jpisf    - nimpp+1, jpi - 1 ) )
(  434)         nis0p1 = max( 1, min(jpisdp1  - nimpp+1, jpi     ) )
(  435)         nis0m1 = max( 1, min(jpisd    - nimpp+1, jpi     ) )
(  436)         nis1m1 = max( 0, min(jpisfm1  - nimpp+1, jpi - 1 ) )
(  437)         nis1m2 = max( 0, min(jpisfm1-1- nimpp+1, jpi - 1 ) )
(  438)         njs0   = max( 1, min(jpjsob   - njmpp+1, jpj     ) )
(  439)         njs1   = max( 0, min(jpjsob   - njmpp+1, jpj - 1 ) )
(  440)         njs0p1 = max( 1, min(jpjsob+1 - njmpp+1, jpj     ) )
(  441)         njs1p1 = max( 0, min(jpjsob+1 - njmpp+1, jpj - 1 ) )
(  442)         IF(lwp) THEN
(  443)            IF( lfbcsouth ) THEN
(  444)               WRITE(numout,*)'     '
(  445)               WRITE(numout,*)'         Specified South Open Boundary'
(  446)            ELSE
(  447)               WRITE(numout,*)'     '
(  448)               WRITE(numout,*)'         Radiative South Open Boundary'
(  449)            END IF
(  450)         END IF






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 30

(  451) !DB: output above plus nimpp njmpp in a readable form:
(  452)         if(DBG) then
(  453)            write(2300+narea,'(a35,6(i3,1x))')'OBC south variables for narea ',narea
(  454)            write(2300+narea,'(a35,6(i3,1x))')'nimpp njmpp ', nimpp, njmpp
(  455)            write(2300+narea,'(a35,6(i3,1x))')'jpi jpj jpjsob jpisd jpisf ',jpi,jpj, jpjsob, jpisd, jpisf 
(  456)            write(2300+narea,'(a35,6(i3,1x))')'nis0, nis1, njs0, njs1: ',nis0, nis1, njs0, njs1
(  457)            write(2300+narea,'(a35,6(i3,1x))')'nis0p1, nis1m1, nis0m1, nis1m2: ',nis0p1, nis1m1, nis0m1, nis1m2
(  458)            write(2300+narea,'(a35,6(i3,1x))')'njs0p1, njs1p1, njs0p1, njs1p1: ',njs0p1, njs1p1, njs0p1, njs1p1
(  459)         endif
(  460) 
(  461) !!DBG
(  462)         if(narea <= jpni) then        !!Assume tiling goes from left-to-right btm-to-top
(  463)            ji = mod(narea,jpni)
(  464)            if(ji==0) ji=jpni
(  465)            s_bnd1(ji) = nis0p1-1+nimpp
(  466)            s_bnd2(ji) = nis1m1-1+nimpp
(  467)         endif
(  468) 
(  469)      END IF
(  470) 
(  471)      if(lk_mpp) then
(  472)         call mpp_sum(s_bnd1,jpni)
(  473)         call mpp_sum(s_bnd2,jpni)
(  474)      endif
(  475) !!DB: Check if any tiles overlap by more than 1 and report
(  476) !!NB: Assumes halo zone = 1
(  477)      if(lwp) then
(  478) 
(  479) !!DBG
(  480)         write(1700,'(A30,100(i4,1x))')'(SOUTH) 1st index in tile:  ',(s_bnd1(ji),ji=1,jpni)
(  481)         write(1700,'(A30,100(i4,1x))')'(SOUTH) Last index in tile: ',(s_bnd2(ji),ji=1,jpni)
(  482)         do ji = 2, jpni
(  483)            if(s_bnd1(ji) < s_bnd2(ji-1)) then
(  484)               write(1700,*)'Problem at South boundary between tiles ',ji, ji-1
(  485)            endif
(  486)         enddo
(  487)      endif
(  488) 
(  489) 
(  490) 
(  491)      ! 3. mask correction for OBCs
(  492)      ! ---------------------------
(  493) 
(  494)      IF( lp_obc_east ) THEN
(  495)         !... (jpjed,jpjefm1),jpieob
(  496)         DO jj = nje0, nje1m1
(  497) # 500
(  500)            DO ji = nie0p1, nie1p1
(  501) # 502
(  502)               bmask(ji,jj) = 0.e0
(  503)            END DO
(  504)         END DO
(  505) 
(  506)         ! ... initilization to zero
(  507)         uemsk(:,:) = 0.e0
(  508)         vemsk(:,:) = 0.e0
(  509)         temsk(:,:) = 0.e0
(  510) !byoung






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 31

(  511)         temsk5(:,:,:) = 0.e0
(  512)         ! ... set 2D mask on East OBC,  Vopt
(  513)         DO ji = nie0 , nie1
(  514)            DO jj = nje0, nje1
(  515)               uemsk(jj,:) = umask(ji,  jj,:)
(  516)               vemsk(jj,:) = vmask(ji+1,jj,:)
(  517)               temsk(jj,:) = tmask(ji+1,jj,:)
(  518)               !byoung
(  519)               DO j1 = 1,5
(  520)                  temsk5(jj,:,j1) = tmask(ji+2-j1,jj,:)
(  521)               END DO
(  522)            END DO
(  523)         END DO
(  524) 
(  525) !!DB: create global east mask
(  526)         do ij = nie0 , nie1
(  527)            do jj = nje0, nje1
(  528)               ii = jj-1+njmpp
(  529)               emaskg(ii) = uemsk(jj,1)
(  530)               emaskg2(ii,:) = uemsk(jj,:)
(  531)            enddo
(  532)         enddo
(  533)         call mpp_sum(emaskg,jpjdta)
(  534) !!DB: correct multiple overlaps 
(  535)         do ii = 1, jpjdta
(  536)            if(emaskg(ii) > 1) emaskg(ii) = 1
(  537)         enddo
(  538) !!DB create global mask using mpp_sum (Check=OK)
(  539)         do ii = 1, jpjdta
(  540)            call mpp_sum(emaskg2(ii,:),jpk)
(  541)         enddo
(  542) !!DB: correct multiple overlaps 
(  543)         do ii = 1, jpjdta
(  544)            do jj = 1, jpk
(  545)               if(emaskg2(ii,jj) > 1)emaskg2(ii,jj) = 1
(  546)            enddo
(  547)         enddo
(  548) !!DBG: output global mask. Checked that all CPUs are the same so only 1 output required
(  549) !        if(lwp) then
(  550) !           do ii = 1, jpjdta
(  551) !              write(6500,'(i3,3x,200(f3.0,1x))')ii, (emaskg2(ii,jj),jj=1,jpk)
(  552) !           enddo
(  553) !        endif
(  554)         
(  555)      END IF
(  556) 
(  557)      IF( lp_obc_west ) THEN
(  558)         ! ... (jpjwd,jpjwfm1),jpiwob
(  559)         DO jj = njw0, njw1m1
(  560)            DO ji = niw0, niw1
(  561)               bmask(ji,jj) = 0.e0
(  562)            END DO
(  563)         END DO
(  564) 
(  565)         ! ... initilization to zero
(  566)         uwmsk(:,:) = 0.e0
(  567)         vwmsk(:,:) = 0.e0
(  568)         twmsk(:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 32

(  569) !byoung
(  570)         twmsk5(:,:,:) = 0.e0
(  571) 
(  572)         ! ... set 2D mask on West OBC,  Vopt
(  573)         DO ji = niw0 , niw1
(  574)            DO jj = njw0, njw1
(  575)               uwmsk(jj,:) = umask(ji,jj,:)
(  576)               vwmsk(jj,:) = vmask(ji,jj,:)
(  577)               twmsk(jj,:) = tmask(ji,jj,:)
(  578)               !byoung
(  579)               DO j1 = 1,5
(  580)                  twmsk5(jj,:,j1) = tmask(ji+j1-1,jj,:)
(  581)               END DO
(  582)            END DO
(  583)         END DO
(  584) !!DB: create global west mask
(  585)         do ij = niw0 , niw1
(  586)            do jj = njw0, njw1
(  587)               ii = jj-1+njmpp
(  588)               wmaskg(ii) = uwmsk(jj,1)
(  589)               wmaskg2(ii,:) = uwmsk(jj,:)
(  590)            enddo
(  591)         enddo
(  592)         call mpp_sum(wmaskg,jpjdta)
(  593) !DB: correct multiple overlaps here
(  594)         do ii = 1, jpjdta
(  595)            if(wmaskg(ii) > 1) wmaskg(ii) = 1
(  596)         enddo
(  597) !!DB create global mask using mpp_sum (Check=OK)
(  598)         do ii = 1, jpjdta
(  599)            call mpp_sum(wmaskg2(ii,:),jpk)
(  600)         enddo
(  601) !!DB: correct multiple overlaps 
(  602)         do ii = 1, jpjdta
(  603)            do jj = 1, jpk
(  604)               if(wmaskg2(ii,jj) > 1)wmaskg2(ii,jj) = 1
(  605)            enddo
(  606)         enddo
(  607) !!DBG: output global mask. Checked that all CPUs are the same so only 1 output required
(  608) !        if(lwp) then
(  609) !           do ii = 1, jpjdta
(  610) !              write(6600,'(i3,3x,200(f3.0,1x))')ii, (wmaskg2(ii,jj),jj=1,jpk)
(  611) !           enddo
(  612) !        endif
(  613) 
(  614)      END IF
(  615) 
(  616)      IF( lp_obc_north ) THEN
(  617)         ! ... jpjnob,(jpind,jpisfm1)
(  618) # 621
(  621)         DO jj = njn0p1, njn1p1
(  622) # 623
(  623)            DO ji = nin0, nin1m1
(  624)               bmask(ji,jj) = 0.e0
(  625)            END DO
(  626)         END DO
(  627) 
(  628)         ! ... initilization to zero






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 33

(  629)         unmsk(:,:) = 0.e0
(  630)         vnmsk(:,:) = 0.e0
(  631)         tnmsk(:,:) = 0.e0
(  632) !sujie
(  633)         tnmsk5(:,:,:) = 0.e0
(  634) 
(  635)         ! ... set 2D mask on North OBC,  Vopt
(  636)         DO jj = njn0, njn1
(  637)            DO ji = nin0, nin1
(  638)               unmsk(ji,:) = umask(ji,jj+1,:)
(  639)               vnmsk(ji,:) = vmask(ji,jj  ,:)
(  640)               tnmsk(ji,:) = tmask(ji,jj+1,:)
(  641) !sujie
(  642)            DO j1 = 1,5
(  643)               tnmsk5(ji,:,j1) = tmask(ji,jj+2-j1,:)
(  644)            END DO
(  645)            END DO
(  646)         END DO
(  647) !!DB: create global north mask
(  648)         do ij = njn0, njn1
(  649)            do jj = nin0, nin1
(  650)               ii = jj-1+nimpp
(  651)               nmaskg(ii) = vnmsk(jj,1)
(  652)            enddo
(  653)         enddo
(  654)         call mpp_sum(nmaskg,jpidta)
(  655) !!DB: correct multiple overlaps here
(  656)         do ii = 1, jpidta
(  657)            if(nmaskg(ii) > 1) nmaskg(ii) = 1
(  658)         enddo
(  659) 
(  660)      END IF
(  661) 
(  662)      IF( lp_obc_south ) THEN
(  663)         ! ... jpjsob,(jpisd,jpisfm1)
(  664)         DO jj = njs0, njs1
(  665)            DO ji = nis0, nis1m1
(  666)               bmask(ji,jj) = 0.e0
(  667)            END DO
(  668)         END DO
(  669) 
(  670)         ! ... initilization to zero
(  671)         usmsk(:,:) = 0.e0
(  672)         vsmsk(:,:) = 0.e0
(  673)         tsmsk(:,:) = 0.e0
(  674) !sujie
(  675)         tsmsk5(:,:,:) = 0.e0
(  676) 
(  677)         ! ... set 2D mask on South OBC,  Vopt
(  678) !        DO jj = njs0, njs1
(  679)         DO jj = njs0, njs1
(  680)            DO ji = nis0, nis1
(  681)               usmsk(ji,:) = umask(ji,jj,:)
(  682)               vsmsk(ji,:) = vmask(ji,jj,:)
(  683)               tsmsk(ji,:) = tmask(ji,jj,:)
(  684)               !sujie
(  685)               DO j1 = 1,5
(  686)                  tsmsk5(ji,:,j1) = tmask(ji,jj+j1-1,:)






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 34

(  687)               END DO
(  688)            END DO
(  689)         END DO
(  690) !!DB: create global south mask
(  691)         do ij = njs0, njs1
(  692)            do jj = nis0, nis1
(  693)               ii = jj-1+nimpp
(  694)               smaskg(ii) = vsmsk(jj,1)
(  695)               smaskg2(ii,:) = vsmsk(jj,:)
(  696)            enddo
(  697)         enddo
(  698)         call mpp_sum(smaskg,jpidta)
(  699) !!DB: correct multiple overlaps here
(  700)         do ii = 1, jpidta
(  701)            if(smaskg(ii) > 1) smaskg(ii) = 1
(  702)         enddo
(  703) !!DB create global mask using mpp_sum (Check=OK)
(  704)         do ii = 1, jpidta
(  705)            call mpp_sum(smaskg2(ii,:),jpk)
(  706)         enddo
(  707) !!DB: correct multiple overlaps 
(  708)         do ii = 1, jpidta
(  709)            do jj = 1, jpk
(  710)               if(smaskg2(ii,jj) > 1)smaskg2(ii,jj) = 1
(  711)            enddo
(  712)         enddo
(  713) !!DBG: output global mask. Checked that all CPUs are the same so only 1 output required
(  714) !        if(lwp) then
(  715) !           do ii = 1, jpidta
(  716) !              write(6700,'(i3,3x,200(f3.0,1x))')ii, (smaskg2(ii,jj),jj=1,jpk)
(  717) !           enddo
(  718) !        endif
(  719) 
(  720)      END IF
(  721) 
(  722) 
(  723) 
(  724) # 810
(  810)      IF ( ln_vol_cst ) THEN
(  811) 
(  812)         ! 3.1 Total lateral surface for each open boundary
(  813)         ! ------------------------------------------------
(  814) 
(  815)         ! ... West open boundary surface
(  816)         IF( lp_obc_west ) THEN
(  817)            DO ji = niw0, niw1
(  818)               DO jj = 1, jpj
(  819)                  obcsurftot = obcsurftot+hu(ji,jj)*e2u(ji,jj)*uwmsk(jj,1)
(  820)               END DO
(  821)            END DO
(  822)         END IF
(  823) 
(  824)         ! ... East open boundary surface
(  825)         IF( lp_obc_east ) THEN
(  826)            DO ji = nie0, nie1
(  827)               DO jj = 1, jpj
(  828)                  obcsurftot = obcsurftot+hu(ji,jj)*e2u(ji,jj)*uemsk(jj,1)
(  829)               END DO






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 35

(  830)            END DO
(  831)         END IF
(  832) 
(  833)         ! ... North open boundary vertical surface
(  834)         IF( lp_obc_north ) THEN
(  835)            DO jj = njn0, njn1
(  836)               DO ji = 1, jpi
(  837)                  obcsurftot = obcsurftot+hv(ji,jj)*e1v(ji,jj)*vnmsk(ji,1)
(  838)               END DO
(  839)            END DO
(  840)         END IF
(  841) 
(  842)         ! ... South open boundary vertical surface
(  843)         IF( lp_obc_south ) THEN
(  844)            DO jj = njs0, njs1
(  845)               DO ji = 1, jpi
(  846)                  obcsurftot = obcsurftot+hv(ji,jj)*e1v(ji,jj)*vsmsk(ji,1)
(  847)               END DO
(  848)            END DO
(  849)         END IF
(  850)         IF( lk_mpp )   CALL mpp_sum( obcsurftot )   ! sum over the global domain
(  851)      ENDIF
(  852) 
(  853) # 854
(  854)      ! 5. Control print on mask
(  855)      !    The extremities of the open boundaries must be in land
(  856)      !    or else correspond to an "ocean corner" between two open boundaries.
(  857)      !    corner 1 is southwest, 2 is south east, 3 is northeast, 4 is northwest.
(  858)      ! --------------------------------------------------------------------------
(  859) 
(  860)      icorner(:)=0
(  861) 
(  862)      ! ... control of the west boundary
(  863)      IF( lp_obc_west ) THEN
(  864)         IF( jpiwob < 2 .OR.  jpiwob >= jpiglo-2 ) THEN
(  865)            IF(lwp) WRITE(numout,*)
(  866)            IF(lwp) WRITE(numout,*) ' E R R O R : jpiwob exceed ', jpiglo-2, 'or less than 2'
(  867)            IF(lwp) WRITE(numout,*) ' ========== '
(  868)            IF(lwp) WRITE(numout,*)
(  869)            nstop = nstop + 1
(  870)         END IF
(  871)         ztestmask(:)=0.
(  872)         DO ji=niw0,niw1
(  873)            IF( (njw0 + njmpp - 1) == jpjwd ) ztestmask(1)=ztestmask(1)+ tmask(ji,njw0,1)
(  874)            IF( (njw1 + njmpp - 1) == jpjwf ) ztestmask(2)=ztestmask(2)+ tmask(ji,njw1,1)
(  875)         END DO
(  876)         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain
(  877) 
(  878)         IF( ztestmask(1) /= 0. ) icorner(1)=icorner(1)+1
(  879)         IF( ztestmask(2) /= 0. ) icorner(4)=icorner(4)+1
(  880)      END IF
(  881) 
(  882)      ! ... control of the east boundary
(  883)      IF( lp_obc_east ) THEN
(  884)         IF( jpieob < 4 .OR.  jpieob >= jpiglo ) THEN
(  885)            IF(lwp) WRITE(numout,*)
(  886)            IF(lwp) WRITE(numout,*) ' E R R O R : jpieob exceed ', jpiglo, ' or less than 4'
(  887)            IF(lwp) WRITE(numout,*) ' ========== '






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 36

(  888)            IF(lwp) WRITE(numout,*)
(  889)            nstop = nstop + 1
(  890)         END IF
(  891)         ztestmask(:)=0.
(  892)         DO ji=nie0p1,nie1p1
(  893)            IF( (nje0 + njmpp - 1) == jpjed ) ztestmask(1)=ztestmask(1)+ tmask(ji,nje0,1)
(  894)            IF( (nje1 + njmpp - 1) == jpjef ) ztestmask(2)=ztestmask(2)+ tmask(ji,nje1,1)
(  895)         END DO
(  896)         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain
(  897) 
(  898)        IF( ztestmask(1) /= 0. ) icorner(2)=icorner(2)+1
(  899)        IF( ztestmask(2) /= 0. ) icorner(3)=icorner(3)+1
(  900)      END IF
(  901) 
(  902)      ! ... control of the north boundary
(  903)      IF( lp_obc_north ) THEN
(  904)         IF( jpjnob < 4 .OR.  jpjnob >= jpjglo ) THEN
(  905)            IF(lwp) WRITE(numout,*)
(  906)            IF(lwp) WRITE(numout,*) ' E R R O R : jpjnob exceed ', jpjglo, ' or less than 4'
(  907)            IF(lwp) WRITE(numout,*) ' ========== '
(  908)            IF(lwp) WRITE(numout,*)
(  909)            nstop = nstop + 1
(  910)         END IF
(  911)         ztestmask(:)=0.
(  912)         DO jj=njn0p1,njn1p1
(  913)            IF( (nin0 + nimpp - 1) == jpind ) ztestmask(1)=ztestmask(1)+ tmask(nin0,jj,1)
(  914)            IF( (nin1 + nimpp - 1) == jpinf ) ztestmask(2)=ztestmask(2)+ tmask(nin1,jj,1)
(  915)         END DO
(  916)         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain
(  917) 
(  918)         IF( ztestmask(1) /= 0. ) icorner(4)=icorner(4)+1
(  919)         IF( ztestmask(2) /= 0. ) icorner(3)=icorner(3)+1
(  920)      END IF
(  921) 
(  922)      ! ... control of the south boundary
(  923)      IF( lp_obc_south ) THEN
(  924)         IF( jpjsob < 2 .OR.  jpjsob >= jpjglo-2 ) THEN
(  925)            IF(lwp) WRITE(numout,*)
(  926)            IF(lwp) WRITE(numout,*) ' E R R O R : jpjsob exceed ', jpjglo-2, ' or less than 2'
(  927)            IF(lwp) WRITE(numout,*) ' ========== '
(  928)            IF(lwp) WRITE(numout,*)
(  929)            nstop = nstop + 1
(  930)         END IF
(  931)         ztestmask(:)=0.
(  932)         DO jj=njs0,njs1
(  933)            IF( (nis0 + nimpp - 1) == jpisd ) ztestmask(1)=ztestmask(1)+ tmask(nis0,jj,1)
(  934)            IF( (nis1 + nimpp - 1) == jpisf ) ztestmask(2)=ztestmask(2)+ tmask(nis1,jj,1)
(  935)         END DO
(  936)         IF( lk_mpp )   CALL mpp_sum( ztestmask, 2 )   ! sum over the global domain
(  937) 
(  938)         IF( ztestmask(1) /= 0. ) icorner(1)=icorner(1)+1
(  939)         IF( ztestmask(2) /= 0. ) icorner(2)=icorner(2)+1
(  940)      END IF
(  941) 
(  942)      IF( icorner(1) == 2 ) THEN
(  943)         IF(lwp) WRITE(numout,*)
(  944)         IF(lwp) WRITE(numout,*) ' South West ocean corner, two open boudaries'
(  945)         IF(lwp) WRITE(numout,*) ' ========== '






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 37

(  946)         IF(lwp) WRITE(numout,*)
(  947)         IF( jpisd /= jpiwob.OR.jpjsob /= jpjwd ) THEN
(  948)            IF(lwp) WRITE(numout,*) ' Open boundaries do not fit, we stop'
(  949)            nstop = nstop + 1
(  950)         END IF
(  951)      ELSE IF( icorner(1) == 1 ) THEN
(  952)              IF(lwp) WRITE(numout,*) ' Open boundaries do not fit at SW corner, we stop'
(  953)              nstop = nstop + 1
(  954)      END IF
(  955) 
(  956)      IF( icorner(2) == 2 ) THEN
(  957)          IF(lwp) WRITE(numout,*)
(  958)          IF(lwp) WRITE(numout,*) ' South East ocean corner, two open boudaries'
(  959)          IF(lwp) WRITE(numout,*) ' ========== '
(  960)          IF(lwp) WRITE(numout,*)
(  961)          IF( jpisf /= jpieob+1.OR.jpjsob /= jpjed ) THEN
(  962)             IF(lwp) WRITE(numout,*) ' Open boundaries do not fit, we stop'
(  963)             nstop = nstop + 1
(  964)          END IF
(  965)      ELSE IF( icorner(2) == 1 ) THEN
(  966)              IF(lwp) WRITE(numout,*) ' Open boundaries do not fit at SE corner, we stop'
(  967)              nstop = nstop + 1
(  968)      END IF
(  969) 
(  970)      IF( icorner(3) == 2 ) THEN
(  971)         IF(lwp) WRITE(numout,*)
(  972)         IF(lwp) WRITE(numout,*) ' North East ocean corner, two open boudaries'
(  973)         IF(lwp) WRITE(numout,*) ' ========== '
(  974)         IF(lwp) WRITE(numout,*)
(  975)         IF( jpinf /= jpieob+1 .OR. jpjnob+1 /= jpjef ) THEN
(  976)            IF(lwp) WRITE(numout,*) ' Open boundaries do not fit, we stop'
(  977)            nstop = nstop + 1
(  978)         END IF
(  979)       ELSE IF( icorner(3) == 1 ) THEN
(  980)               IF(lwp) WRITE(numout,*) ' Open boundaries do not fit at NE corner, we stop'
(  981)               nstop = nstop + 1
(  982)       END IF
(  983) 
(  984)      IF( icorner(4) == 2 ) THEN
(  985)         IF(lwp) WRITE(numout,*)
(  986)         IF(lwp) WRITE(numout,*) ' North West ocean corner, two open boudaries'
(  987)         IF(lwp) WRITE(numout,*) ' ========== '
(  988)         IF(lwp) WRITE(numout,*)
(  989)         IF( jpind /= jpiwob.OR.jpjnob+1 /= jpjwf ) THEN
(  990)            IF(lwp) WRITE(numout,*) ' Open boundaries do not fit, we stop'
(  991)            nstop = nstop + 1
(  992)         END IF
(  993)       ELSE IF( icorner(4) == 1 ) THEN
(  994)               IF(lwp) WRITE(numout,*) ' Open boundaries do not fit at NW corner, we stop'
(  995)               nstop = nstop + 1
(  996)       END IF
(  997) 
(  998)      ! 6. Initialization of open boundary variables (u, v, bsf, t, s)
(  999)      ! --------------------------------------------------------------
( 1000)      !   only if at least one boundary is  radiative
( 1001) 
( 1002)      ! ... Restart from restart.obc
( 1003)      IF ( inumfbc < nbobc .AND.  ln_rstart ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 38

( 1004)         CALL obc_rst_lec
( 1005)      ELSE
( 1006) 
( 1007)          ! ... Initialization to zero of radiation arrays.
( 1008)          !     Those have dimensions of local subdomains
( 1009) 
( 1010)          bebnd(:,:,:)   = 0.e0   ;   bnbnd(:,:,:)   = 0.e0
( 1011)          uebnd(:,:,:,:) = 0.e0   ;   unbnd(:,:,:,:) = 0.e0
( 1012)          vebnd(:,:,:,:) = 0.e0   ;   vnbnd(:,:,:,:) = 0.e0
( 1013)          tebnd(:,:,:,:) = 0.e0   ;   tnbnd(:,:,:,:) = 0.e0
( 1014)          sebnd(:,:,:,:) = 0.e0   ;   snbnd(:,:,:,:) = 0.e0
( 1015) 
( 1016)          bwbnd(:,:,:)   = 0.e0   ;   bsbnd(:,:,:)   = 0.e0
( 1017)          uwbnd(:,:,:,:) = 0.e0   ;   usbnd(:,:,:,:) = 0.e0
( 1018)          vwbnd(:,:,:,:) = 0.e0   ;   vsbnd(:,:,:,:) = 0.e0
( 1019)          twbnd(:,:,:,:) = 0.e0   ;   tsbnd(:,:,:,:) = 0.e0
( 1020)          swbnd(:,:,:,:) = 0.e0   ;   ssbnd(:,:,:,:) = 0.e0
( 1021) 
( 1022)      END IF
( 1023) 
( 1024) 
( 1025) # 1030
( 1030)      ! 8. Control print
( 1031)      ! ... control of the east boundary
( 1032)      IF( lp_obc_east ) THEN
( 1033)         istop = 0
( 1034)         IF( jpieob < 4 .OR.  jpieob >= jpiglo ) THEN
( 1035)            IF(lwp) WRITE(numout,cform_err)
( 1036)            IF(lwp) WRITE(numout,*) '            jpieob exceed ', jpim1, ' or less than 4'
( 1037)            istop = istop + 1
( 1038)         END IF
( 1039) 
( 1040)         IF( lk_mpp ) THEN
( 1041)            ! ...
( 1042)            IF( nimpp > jpieob-5) THEN
( 1043)               IF(lwp) WRITE(numout,cform_err)
( 1044)               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the East OBC'
( 1045)               IF(lwp) WRITE(numout,*) '        nimpp must be < jpieob-5'
( 1046)               istop = istop + 1
( 1047)            ENDIF
( 1048)         ELSE
( 1049) 
( 1050)            ! ... stop if  e r r o r (s)   detected
( 1051)            IF( istop /= 0 ) THEN
( 1052)               IF(lwp)WRITE(numout,*)
( 1053)               IF(lwp)WRITE(numout,*) istop,' E R R O R (S) detected : stop'
( 1054)               IF(lwp)WRITE(numout,*) ' =============== '
( 1055)               IF(lwp)WRITE(numout,*)
( 1056)               nstop = nstop + 1
( 1057)            ENDIF
( 1058)         ENDIF
( 1059)      ENDIF
( 1060) 
( 1061)      ! ... control of the west boundary
( 1062)      IF( lp_obc_west ) THEN
( 1063)         istop = 0
( 1064)         IF( jpiwob < 2 .OR.  jpiwob >= jpiglo ) THEN
( 1065)            IF(lwp) WRITE(numout,cform_err)






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 39

( 1066)            IF(lwp) WRITE(numout,*) '            jpiwob exceed ', jpim1, ' or less than 2'
( 1067)            istop = istop + 1
( 1068)         END IF
( 1069) 
( 1070)         IF( lk_mpp ) THEN
( 1071)            IF( (nimpp < jpiwob+5) .AND. (nimpp > 1) ) THEN
( 1072)               IF(lwp) WRITE(numout,cform_err)
( 1073)               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the West OBC'
( 1074)               IF(lwp) WRITE(numout,*) '        nimpp must be > jpiwob-5 or =1'
( 1075)               istop = istop + 1
( 1076)            ENDIF
( 1077)         ELSE
( 1078) 
( 1079)            ! ... stop if  e r r o r (s)   detected
( 1080)            IF( istop /= 0 ) THEN
( 1081)               IF(lwp)WRITE(numout,*)
( 1082)               IF(lwp)WRITE(numout,*) istop,' E R R O R (S) detected : stop'
( 1083)               IF(lwp)WRITE(numout,*) ' =============== '
( 1084)               IF(lwp)WRITE(numout,*)
( 1085)               nstop = nstop + 1
( 1086)            ENDIF
( 1087)         ENDIF
( 1088)      ENDIF
( 1089) 
( 1090)      ! control of the north boundary
( 1091)      IF( lp_obc_north ) THEN
( 1092)         istop = 0
( 1093)         IF( jpjnob < 4 .OR.  jpjnob >= jpjglo ) THEN
( 1094)            IF(lwp) WRITE(numout,cform_err)
( 1095)            IF(lwp) WRITE(numout,*) '          jpjnob exceed ', jpjm1,' or less than 4'
( 1096)            istop = istop + 1
( 1097)         END IF
( 1098) 
( 1099)         IF( lk_mpp ) THEN
( 1100)            IF( njmpp > jpjnob-5) THEN
( 1101)               IF(lwp) WRITE(numout,cform_err)
( 1102)               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the North OBC'
( 1103)               IF(lwp) WRITE(numout,*) '        njmpp must be < jpjnob-5'
( 1104)               istop = istop + 1
( 1105)            ENDIF
( 1106)         ELSE
( 1107) 
( 1108)            ! ... stop if  e r r o r (s)   detected
( 1109)            IF( istop /= 0 ) THEN
( 1110)               IF(lwp)WRITE(numout,*)
( 1111)               IF(lwp)WRITE(numout,*) istop,' E R R O R (S) detected : stop'
( 1112)               IF(lwp)WRITE(numout,*) ' =============== '
( 1113)               IF(lwp)WRITE(numout,*)
( 1114)               nstop = nstop + 1
( 1115)            ENDIF
( 1116)         ENDIF
( 1117)      ENDIF
( 1118) 
( 1119)      ! control of the south boundary
( 1120)      IF( lp_obc_south ) THEN
( 1121)         istop = 0
( 1122)         IF( jpjsob < 2 .OR. jpjsob >= jpjglo ) THEN
( 1123)            IF(lwp) WRITE(numout,cform_err)






PGF90 (Version     10.5)          09/23/2010  09:50:07      page 40

( 1124)            IF(lwp) WRITE(numout,*) '          jpjsob exceed ', jpjm1,' or less than 2'
( 1125)            istop = istop + 1
( 1126)         END IF
( 1127) 
( 1128)         IF( lk_mpp ) THEN
( 1129)            IF( (njmpp < jpjsob+5) .AND. (njmpp > 1) ) THEN
( 1130)               IF(lwp) WRITE(numout,cform_err)
( 1131)               IF(lwp) WRITE(numout,*) '        A sub-domain is too close to the South OBC'
( 1132)               IF(lwp) WRITE(numout,*) '        njmpp must be > jpjsob+5 or =1'
( 1133)               istop = istop + 1
( 1134)            ENDIF
( 1135)         ELSE
( 1136) 
( 1137)            ! ... stop if  e r r o r (s)   detected
( 1138)            IF( istop /= 0 ) THEN
( 1139)               IF(lwp)WRITE(numout,*)
( 1140)               IF(lwp)WRITE(numout,*) istop,' E R R O R (S) detected : stop'
( 1141)               IF(lwp)WRITE(numout,*) ' =============== '
( 1142)               IF(lwp)WRITE(numout,*)
( 1143)               nstop = nstop + 1
( 1144)            ENDIF
( 1145)         ENDIF
( 1146)      ENDIF
( 1147) 
( 1148)   END SUBROUTINE obc_init
( 1149) 
( 1150) 
( 1151) # 1159
( 1159)   !!=================================================================================
( 1160) END MODULE obcini
