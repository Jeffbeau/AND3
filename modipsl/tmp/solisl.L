


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

Filename: solisl.F90

(    1) MODULE solisl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  solisl  ***
(    4)    !! Ocean island : specific treatment of island in rigid-lid case
(    5)    !!==============================================================================
(    6) #if defined key_islands
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_islands' :                           islands in rigid-lid case
(    9)    !!----------------------------------------------------------------------
(   10)    !!   isl_dom     : locate islands in the domain
(   11)    !!   isl_pri     : control print of island gridpoint position
(   12)    !!   isl_pth     : Compute coeff. associated with path round each island
(   13)    !!   isl_mat     : Compute the matrix associated with islands
(   14)    !!   isl_bsf     : Compute the barotropic streamfunction of each island
(   15)    !!   isl_dyn_spg : Update the barotropic streamfunction trend with the
(   16)    !!                 Island contribution (call by dyn_spg)
(   17)    !!   isl_stp_ctl : print island information (call by stp_ctl)
(   18)    !!----------------------------------------------------------------------
(   19)    !! * Modules used
(   20)    USE oce             ! ocean dynamics and tracers variables
(   21)    USE dom_oce         ! ocean space and time domain variables 
(   22)    USE in_out_manager  ! I/O manager
(   23)    USE sol_oce         ! ocean solver 
(   24)    USE obc_oce         ! ocean open boundary condition
(   25)    USE lib_mpp         ! distributed memory computing
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Routine accessibility
(   31)    PUBLIC isl_dom        ! routine called by solver_init
(   32)    PUBLIC isl_mat        ! routine called by solver_init
(   33)    PUBLIC isl_bsf        ! routine called by solver_init
(   34)    PUBLIC isl_dyn_spg    ! routine called by dyn_spg
(   35)    PUBLIC isl_stp_ctl    ! routine called by stp_ctl
(   36) 
(   37)    !! * Shared module variables
(   38)    LOGICAL, PUBLIC, PARAMETER ::   lk_isl = .TRUE.    !: 'key_islands' flag
(   39) 
(   40)    !! * module variable
(   41)    INTEGER :: numisl = 11      ! logical unit for island file only used
(   42)    !                           ! here during the initialization phase
(   43)    INTEGER ::   &
(   44)       nimlu, njmlu, nkmlu,  &  ! i-j-k-dimensions read
(   45)       nlmlu, nmmlu, nnmlu      ! read islands number
(   46)    INTEGER, DIMENSION(jpnisl,0:4,jpisl) ::   &






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 2

(   47)       miisl, mjisl             ! position of island grid-points
(   48)    INTEGER, DIMENSION(0:4,jpisl) ::   &
(   49)       mnisl                    ! number of grid-points for each island
(   50)  
(   51)    REAL(wp) ::    &
(   52)       replu                    ! read absolute precision
(   53)    REAL(wp), DIMENSION(jpisl,jpisl) ::   &
(   54)       aisl, aislm1             ! island matrix and its inverse
(   55)    REAL(wp), DIMENSION(jpi,jpj,jpisl) ::   &
(   56)       bsfisl                   ! barotropic streamfunction of island
(   57)    REAL(wp), DIMENSION(jpi,jpj,2,jpisl) ::   &
(   58)       acisl1, acisl2           ! coef. to compute circulations round islands
(   59)    REAL(wp), DIMENSION(jpisl) ::   &
(   60)       bisl,                 &  ! second member of island linear system
(   61)       visl                     ! trend of island stream function
(   62)    !!----------------------------------------------------------------------
(   63)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   64)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solisl.F90,v 1.6 2005/12/12 14:18:06 opalod Exp $ 
(   65)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   66)    !!----------------------------------------------------------------------
(   67) 
(   68) CONTAINS
(   69) 
(   70)    SUBROUTINE isl_dom
(   71)       !!----------------------------------------------------------------------
(   72)       !!                  ***  ROUTINE isl_dom  ***
(   73)       !!                    
(   74)       !! ** Purpose :   Locate island grid-points from mbathy array
(   75)       !!
(   76)       !! ** Method  :   The coordinates of ocean grid-points round an island 
(   77)       !!      are found from mbathy array read or computed in dommba.
(   78)       !!      we first compute zwb, an ocean/land mask defined as follows:
(   79)       !!              zwb(i,j)  =  0. over the main land and the ocean
(   80)       !!                        = -n  over the nth island
(   81)       !!      With the proper boundary conditions (defined by nperio)
(   82)       !!        Note: IF i,j are the coordinates of an ocean grid-point west
(   83)       !!      or east of a island, the corresponding coordinates miisl(ip,4,n)
(   84)       !!      or mjisl(ip,2,n) are those of the western or southern side of
(   85)       !!      the island (i.e. i+1 ou j+1, respectively)
(   86)       !!
(   87)       !! ** Action  :   compute mnisl, miisl, mjisl arrays defined as:
(   88)       !!      mnisl(i,n) : nb of grid-points along (i=0), north (i=1), north
(   89)       !!      (i=1), south (i=2), east (i=3), or west (i=4) of the nth island
(   90)       !!      miisl(ip,i,n), mjisl(ip,i,n) : (i,j) index of ipth u- or v-point
(   91)       !!      along (i=0), north (i=1),south (i=2), east (i=3), or west (i=4)
(   92)       !!      of the nth island
(   93)       !!
(   94)       !! History :
(   95)       !!    4.0  !  88-03  (G. Madec)  Original code
(   96)       !!    6.0  !  96-01  (G. Madec)  Suppress common workspace, use of bmask
(   97)       !!    8.5  !  96-01  (G. Madec)  Free form, F90
(   98)       !!----------------------------------------------------------------------
(   99)       !! * Local declarations
(  100)       INTEGER ::   ji, jj, jn, jnil   ! dummy loop indices
(  101)       INTEGER ::   ind, inilt, ip, ipn, ips, ipe, ipw, ii, ij, iju
(  102)       INTEGER ::   iista, iiend, ijsta, ijend, ijstm1, ijenm1
(  103)       INTEGER ::   isrchne
(  104)       INTEGER, DIMENSION(jpj) ::   indil






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 3

(  105)       INTEGER, DIMENSION(jpi,jpj) ::   idil
(  106) 
(  107)       REAL(wp) znil
(  108)       REAL(wp), DIMENSION(jpi,jpj) ::   zwb
(  109)       !!----------------------------------------------------------------------
(  110) 
(  111) 
(  112)       ! 0. Islands index computed from mbathy
(  113)       ! -------------------------------------
(  114)       ! zwb=0 over the continent and ocean, =-n over the nth island
(  115)       
(  116)       ! Computation
(  117)       DO jj = 1, jpjm1
(  118)          DO ji = 1, jpim1
(  119)             zwb(ji,jj) = MIN( 0 , mbathy(ji,jj  ), mbathy(ji+1,jj  ),   &
(  120)                                   mbathy(ji,jj+1), mbathy(ji+1,jj+1)  )
(  121)          END DO
(  122)       END DO
(  123)       
(  124)       ! Lateral boundary conditions on zwb
(  125) 
(  126)       !  mono- or macro-tasking environnement:
(  127)       IF( nperio == 2 ) THEN
(  128)          zwb(:, 1 ) = zwb(:, 2 )
(  129)       ELSEIF( nperio == 3 .OR. nperio == 4) THEN
(  130) ! om 
(  131) !$$$          DO ji = 1, jpim1
(  132) ! On ne peut pas partir de ji=1, car alors jiu=jpi, et cette valeur
(  133) ! n'est pas encore initialiee. 
(  134) ! Le cas ji=1 est de toute facon traite a partir de la ligne 135
(  135)          DO ji = 2, jpim1
(  136)             iju = jpi-ji+1
(  137)             zwb(ji,jpj  ) = zwb(iju,jpj-3)
(  138)             zwb(ji,jpjm1) = zwb(iju,jpj-2)
(  139)          END DO
(  140)       ELSEIF( nperio == 5 .OR. nperio == 6) THEN
(  141)          DO ji = 1, jpim1
(  142)             iju = jpi-ji
(  143)             zwb(ji,jpj  ) = zwb(iju,jpj-2)
(  144)          END DO
(  145)          DO ji = jpi/2+1, jpi-1
(  146)             iju = jpi-ji
(  147)             zwb(ji,jpjm1) = zwb(iju,jpjm1)
(  148)          END DO
(  149)       ELSE
(  150)          zwb(:,jpj) = 0.e0
(  151)       ENDIF
(  152) 
(  153)       IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
(  154)          zwb( 1 ,:) = zwb(jpim1,:)
(  155)          zwb(jpi,:) = zwb(  2  ,:)
(  156)       ELSE
(  157)          zwb( 1 ,:) = 0.e0
(  158)          zwb(jpi,:) = 0.e0
(  159)       ENDIF
(  160)       IF( lk_mpp )   CALL lbc_lnk( zwb, 'G', 1. )
(  161) 
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 4

(  163)       ! 1. Initialization for the search of island grid-points
(  164)       ! ------------------------------------------------------
(  165) 
(  166)       IF( lk_mpp ) THEN
(  167) 
(  168)          ! Mpp : The overlap region are not taken into account
(  169)          ! (islands bondaries are searched over subdomain only)
(  170)          iista =  1   + jpreci
(  171)          iiend = nlci - jpreci
(  172)          ijsta =  1   + jprecj
(  173)          ijend = nlcj - jprecj
(  174)          ijstm1=  1   + jprecj
(  175)          ijenm1= nlcj - jprecj
(  176)          IF( nbondi == -1 .OR. nbondi == 2 ) THEN
(  177)             iista  = 1
(  178)          ENDIF
(  179)          IF( nbondi ==  1 .OR. nbondi == 2 ) THEN
(  180)             iiend  = nlci
(  181)          ENDIF
(  182)          IF( nbondj == -1 .OR. nbondj == 2 ) THEN
(  183)             ijsta  = 1
(  184)             ijstm1 = 2
(  185)          ENDIF
(  186)          IF( nbondj ==  1 .OR. nbondj == 2 ) THEN
(  187)             ijend  = nlcj
(  188)             ijenm1 = nlcj-1
(  189)          ENDIF
(  190)          IF( npolj == 3 .OR. npolj == 4 ) THEN
(  191)             ijend  = nlcj-2
(  192)             ijenm1 = nlcj-2
(  193)          ENDIF 
(  194)       ELSE
(  195)          ! mono- or macro-tasking environnement: full domain scan
(  196)          iista  = 1
(  197)          iiend  = jpi
(  198)          ijsta  = 1
(  199)          ijstm1 = 2
(  200)          IF( nperio == 3 .OR. nperio == 4 ) THEN
(  201)             ijend  = jpj-2
(  202)             ijenm1 = jpj-2
(  203)          ELSEIF( nperio == 5 .OR. nperio == 6 ) THEN
(  204)             ijend  = jpj-1
(  205)             ijenm1 = jpj-1
(  206)          ELSE
(  207)             ijend  = jpj
(  208)             ijenm1 = jpj-1
(  209)          ENDIF
(  210)       ENDIF
(  211) 
(  212) 
(  213)       ! 2. Loop over island
(  214)       ! -------------------
(  215) 
(  216)       DO jnil = 1, jpisl
(  217) 
(  218)          ! 2.1 Initialization to zero of miisl, mjisl of the jnil island
(  219) 
(  220)          miisl(:,:,jnil) = 0






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 5

(  221)          mjisl(:,:,jnil) = 0
(  222)          
(  223)          ! 2.2 Search grid-points of island jnil
(  224)  
(  225)          indil(:)   = 0
(  226)          idil (:,:) = 0
(  227)          
(  228)          znil = - FLOAT( jnil )
(  229)          DO jj = ijsta, ijend
(  230)             indil(jj) = 0
(  231)             ind = 0
(  232)             DO ji = iista, iiend
(  233)                IF( zwb(ji,jj) == znil ) THEN
(  234)                   ind = ind + 1
(  235)                   idil(ind,jj) = ji
(  236)                   indil(jj) = indil(jj) + 1
(  237)                ENDIF
(  238)             END DO
(  239)          END DO
(  240)          
(  241)          ! 2.3 Check the number of island
(  242)          
(  243)          inilt = 0
(  244)          DO jj = ijsta, ijend
(  245)             inilt = inilt + indil(jj)
(  246)          END DO
(  247)          IF( lk_mpp )   CALL mpp_sum( inilt )   ! sum over the global domain
(  248) 
(  249)          IF( inilt == 0 ) THEN
(  250)             IF(lwp) THEN
(  251)                WRITE(numout,*) ' isldom: there is not island number: ', jnil,' while jpisl= ', jpisl
(  252)                WRITE(numout,*) ' change parameter.h'
(  253)             ENDIF
(  254)             STOP 'isldom'      !cr replace by nstop
(  255)          ENDIF
(  256)          
(  257)          ! 2.4 Coastal island grid-points (miisl,mjisl)
(  258)          
(  259)          ip  = 0
(  260)          ipn = 0
(  261)          ips = 0
(  262)          ipe = 0
(  263)          ipw = 0
(  264) 
(  265)          ! South line (ij=1)
(  266)          ij = 1
(  267)          DO jn = 1, indil(ij)
(  268)             ii = idil(jn,ij)
(  269)             IF( (ij+njmpp-1) == 1 .AND. ii > jpreci .AND. ii < (nlci-jpreci+1) ) THEN
(  270)                IF(  zwb(ii-1,ij) * zwb(ii+1,ij) * zwb(ii,ij+1)  == 0. ) THEN
(  271)                   ip = ip + 1
(  272)                   miisl(ip,0,jnil) = ii
(  273)                   mjisl(ip,0,jnil) = ij
(  274)                   IF( zwb(ii-1,ij) == 0. ) THEN
(  275)                      ipw = ipw + 1
(  276)                      miisl(ipw,4,jnil) = ii
(  277)                      mjisl(ipw,4,jnil) = ij
(  278)                   ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 6

(  279)                   IF( zwb(ii+1,ij) == 0. ) THEN
(  280)                      ipe = ipe+1
(  281)                      IF( (nperio == 1 .OR. nperio == 4.OR. nperio == 6)  &
(  282)                         .AND. ii == (nlci-jpreci) ) THEN
(  283)                         miisl(ipe,3,jnil) = 1+jpreci
(  284)                      ELSE
(  285)                         miisl(ipe,3,jnil) = ii + 1
(  286)                      ENDIF
(  287)                      mjisl(ipe,3,jnil) = ij
(  288)                   ENDIF
(  289)                   IF( zwb(ii,ij+1) == 0. ) THEN
(  290)                      ipn = ipn+1
(  291)                      miisl(ipn,1,jnil) = ii
(  292)                      mjisl(ipn,1,jnil) = ij + 1
(  293)                   ENDIF
(  294)                ENDIF
(  295)             ENDIF
(  296)          END DO
(  297)          
(  298)          ! Middle lines (2=<jj=<jpjm1 or jpj-2 if north fold b.c.)
(  299)          DO jj = ijstm1, ijenm1
(  300)             DO jn = 1, indil(jj)
(  301)                ii = idil(jn,jj)
(  302)                IF( ii > jpreci .AND. ii < (nlci-jpreci+1) ) THEN
(  303)                   IF( (zwb(ii-1, jj )*zwb(ii+1, jj )*   &
(  304)                        zwb( ii ,jj-1)*zwb( ii ,jj+1)) == 0. ) THEN
(  305)                      ip = ip + 1
(  306)                      miisl(ip,0,jnil) = ii
(  307)                      mjisl(ip,0,jnil) = jj
(  308)                      IF( zwb(ii-1,jj) == 0. ) THEN
(  309)                         ipw = ipw + 1
(  310)                         miisl(ipw,4,jnil) = ii
(  311)                         mjisl(ipw,4,jnil) = jj
(  312)                      ENDIF
(  313)                      IF( zwb(ii+1,jj) == 0. ) THEN
(  314)                         ipe = ipe + 1
(  315)                         IF((nperio == 1.OR.nperio == 4.OR.nperio == 6)   &
(  316)                            .AND.ii == (nlci-jpreci) ) THEN
(  317)                            miisl(ipe,3,jnil) = 1 + jpreci
(  318)                         ELSE
(  319)                            miisl(ipe,3,jnil) = ii + 1
(  320)                         ENDIF
(  321)                         mjisl(ipe,3,jnil) = jj
(  322)                      ENDIF
(  323)                      IF( zwb(ii,jj-1) == 0. ) THEN
(  324)                         ips = ips + 1
(  325)                         miisl(ips,2,jnil) = ii
(  326)                         mjisl(ips,2,jnil) = jj
(  327)                      ENDIF
(  328)                      IF( zwb(ii,jj+1) == 0. ) THEN
(  329)                         ipn = ipn + 1
(  330)                         miisl(ipn,1,jnil) = ii
(  331)                         mjisl(ipn,1,jnil) = jj + 1
(  332)                      ENDIF
(  333)                   ENDIF
(  334)                ENDIF
(  335)             END DO
(  336)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 7

(  337)          
(  338)          ! North line (jj=jpj) only if not north fold b.c.
(  339)          IF( nperio /= 3 .AND. nperio /= 4 .AND. nperio /= 5 .AND. nperio /= 6 ) THEN
(  340)             ij = jpj
(  341)             DO jn = 1, indil(ij)
(  342)                ii = idil(jn,ij)
(  343)                IF( (ij+njmpp-1) == jpjglo .AND. ii > jpreci .AND. ii < (nlci-jpreci+1) ) THEN
(  344)                   IF( (zwb(ii-1, ij )*zwb(ii+1, ij )* zwb( ii ,ij-1) ) == 0. ) THEN
(  345)                      ip = ip+1
(  346)                      miisl(ip,0,jnil) = ii
(  347)                      mjisl(ip,0,jnil) = ij
(  348)                      IF( zwb(ii-1,ij) == 0. ) THEN
(  349)                         ipw = ipw+1
(  350)                         miisl(ipw,4,jnil) = ii
(  351)                         mjisl(ipw,4,jnil) = ij
(  352)                      ENDIF
(  353)                      IF( zwb(ii+1,ij) == 0. ) THEN
(  354)                         ipe = ipe+1
(  355)                         IF( (nperio == 1) .AND. ii == (nlci-jpreci) ) THEN
(  356)                            miisl(ipe,3,jnil) = 1+jpreci
(  357)                         ELSE
(  358)                            miisl(ipe,3,jnil) = ii+1
(  359)                         ENDIF
(  360)                         mjisl(ipe,3,jnil) = ij
(  361)                      ENDIF
(  362)                      IF( zwb(ii,ij-1) == 0. ) THEN
(  363)                         ips = ips+1
(  364)                         miisl(ips,2,jnil) = ii
(  365)                         mjisl(ips,2,jnil) = ij
(  366)                      ENDIF
(  367)                   ENDIF
(  368)                ENDIF
(  369)             END DO
(  370)          ENDIF
(  371)          
(  372)          mnisl(0,jnil) = ip
(  373)          mnisl(1,jnil) = ipn
(  374)          mnisl(2,jnil) = ips
(  375)          mnisl(3,jnil) = ipe
(  376)          mnisl(4,jnil) = ipw
(  377)          
(  378)          ! Take account of redundant points
(  379)          
(  380)          IF( lk_mpp )   CALL mpp_sum( ip )   ! sum over the global domain
(  381)          
(  382)          IF( ip > jpnisl ) THEN
(  383)             IF(lwp) THEN
(  384)                WRITE(numout,*) ' isldom: the island ',jnil,' has ',   &
(  385)                   mnisl(0,jnil),' grid-points, while jpnisl= ', jpnisl,ip
(  386)                WRITE(numout,*) ' change parameter.h'
(  387)             ENDIF
(  388)             STOP 'isldom'    !cr => nstop
(  389)          ENDIF
(  390)          
(  391)          ! 2.5 Set to zero the grid-points of the jnil island in x
(  392)          
(  393)          DO jj = 1, jpj
(  394)             DO ji = 1, jpi






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 8

(  395)                IF( zwb(ji,jj)+FLOAT(jnil) == 0. ) zwb(ji,jj) = 0. 
(  396)             END DO
(  397)          END DO
(  398)          
(  399)       END DO
(  400)       
(  401) 
(  402)       ! 3. Check the number of island
(  403)       ! -----------------------------
(  404) 
(  405)       inilt = isrchne( jpij, zwb(1,1), 1, 0. )
(  406)       IF( lk_mpp )   CALL mpp_min( inilt )   ! min over the global domain
(  407) 
(  408)       IF( inilt /= jpij+1 ) THEN
(  409)          IF(lwp) THEN
(  410)             WRITE(numout,*) ' isldom: there is at least one more ',   &
(  411)                   'island in the domain and jpisl=', jpisl
(  412)             WRITE(numout,*) ' change parameter.h'
(  413)          ENDIF
(  414)          STOP 'isldom'
(  415)       ENDIF
(  416) 
(  417) 
(  418)       ! 4. Print of island parametres and arrays
(  419)       ! ----------------------------------------
(  420) 
(  421)       CALL isl_pri
(  422) 
(  423) 
(  424)       ! 5. Array for computation of circulation arround islands
(  425)       ! -------------------------------------------------------
(  426) 
(  427)       CALL isl_pth
(  428) 
(  429)    END SUBROUTINE isl_dom
(  430) 
(  431) 
(  432)    SUBROUTINE isl_pri
(  433)       !!----------------------------------------------------------------------
(  434)       !!                  ***  ROUTINE isl_pri  ***
(  435)       !!              
(  436)       !! ** Purpose :   Print islands variables and islands arrays
(  437)       !!
(  438)       !! ** Method  :
(  439)       !!
(  440)       !! History :
(  441)       !!   1.0  !  88-03  (G. Madec)  Original code
(  442)       !!   8.5  !  02-08  (G. Madec)  Free form, F90
(  443)       !!----------------------------------------------------------------------
(  444)       !! * Local declarations
(  445)       INTEGER ::   ji, jni, jnp      ! dummy loop variables
(  446)       INTEGER ::   &
(  447)          ip, ipn, ips, ipe, ipw      ! temporary integers
(  448)       !!----------------------------------------------------------------------
(  449)       
(  450)       IF(lwp) THEN
(  451)          WRITE(numout,*) ' '
(  452)          WRITE(numout,*) '*** islpri number of islands : jpisl =',jpisl






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 9

(  453)       ENDIF
(  454) 
(  455)       DO jni = 1, jpisl
(  456)          ip  = mnisl(0,jni)
(  457)          ipn = mnisl(1,jni)
(  458)          ips = mnisl(2,jni)
(  459)          ipe = mnisl(3,jni)
(  460)          ipw = mnisl(4,jni)
(  461)          IF( lk_mpp ) THEN
(  462)             CALL mpp_sum( ip  )   ! sums over the global domain
(  463)             CALL mpp_sum( ipn )
(  464)             CALL mpp_sum( ips )
(  465)             CALL mpp_sum( ipe )
(  466)             CALL mpp_sum( ipw )
(  467)          ENDIF
(  468)          IF(lwp) THEN
(  469)             WRITE(numout,9000) jni
(  470)             WRITE(numout,9010) ip, ipn, ips, ipe, ipw
(  471)             WRITE(numout,9020)
(  472)             DO jnp = 1, mnisl(0,jni)
(  473)                WRITE(numout,9030) jnp, ( miisl(jnp,ji,jni)+nimpp-1,   &
(  474)                                          mjisl(jnp,ji,jni)+njmpp-1, ji = 0, 4 )
(  475)             END DO
(  476)          ENDIF
(  477)       END DO
(  478) 
(  479)       ! FORMAT   !!cr => no more format
(  480)  9000 FORMAT(/, /, 'island number= ', i2 )
(  481)  9010 FORMAT(/, 'npil=',i4,' npn=',i3,' nps=',i3,' npe=',i3,' npw=',i3 )
(  482)  9020 FORMAT(/,'     * isl point *  point n  *  point s  *  point e ', '*  point w  *')
(  483)  9030 FORMAT(i4,' * (',i3,',',i3,') * (',i3,',',i3,') * (',i3,',',i3,   &
(  484)           ') * (',i3,',',i3,') * (',i3,',',i3,') *' )
(  485) 
(  486)    END SUBROUTINE isl_pri
(  487) 
(  488) 
(  489)    SUBROUTINE isl_pth
(  490)       !!----------------------------------------------------------------------
(  491)       !!                  ***  ROUTINE isl_pth  ***
(  492)       !!             
(  493)       !! ** Purpose :   intialize arrays for the computation of streamfunction
(  494)       !!      around islands
(  495)       !!
(  496)       !! ** Method  :
(  497)       !!
(  498)       !! ** Action  : - acisl1(i,j,ii,n): coefficient n-s (ii=1) and e-w (ii=2)
(  499)       !!                for the calculation of mu and mv around the island n      
(  500)       !!              - acisl2(i,j,ii,n): coefficient n-s (ii=1) and e-w (ii=2)
(  501)       !!                for the calculation of bsfd around the island n
(  502)       !!
(  503)       !! History :
(  504)       !!   1.0  !  88-03  (G. Madec)  Original code
(  505)       !!   8.5  !  02-08  (G. Madec)  Free form, F90
(  506)       !!----------------------------------------------------------------------
(  507)       !! * Local declarations
(  508)       INTEGER ::   jni, jii, jnp    ! dummy loop indices
(  509)       INTEGER ::   ii, ij           ! temporary integers
(  510)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 10

(  511)       
(  512)       ! 1. Initialisation
(  513)       ! -----------------
(  514)       acisl1(:,:,:,:) = 0.e0
(  515)       acisl2(:,:,:,:) = 0.e0
(  516)       bsfisl(:,:,:)   = 0.e0
(  517) 
(  518)       ! 2. Coefficient arrays
(  519)       ! ---------------------
(  520)       DO jni = 1, jpisl
(  521)          DO jii = 1, 4
(  522)             DO jnp = 1,mnisl(jii,jni)
(  523)                ii = miisl(jnp,jii,jni)
(  524)                ij = mjisl(jnp,jii,jni)
(  525)                IF( jii <= 2 ) THEN
(  526)                   ! north and south points
(  527)                   acisl1(ii,ij,1,jni) = float( 2*jii-3) * e1u(ii,ij)
(  528)                   acisl2(ii,ij,1,jni) = float( 2*jii-3) * e1u(ii,ij) * hur(ii,ij) / e2u(ii,ij)
(  529)                ELSE
(  530)                   ! east and west points
(  531)                   acisl1(ii,ij,2,jni) = float(-2*jii+7) * e2v(ii,ij)
(  532)                   acisl2(ii,ij,2,jni) = float(-2*jii+7) * e2v(ii,ij) * hvr(ii,ij) / e1v(ii,ij)
(  533)                ENDIF
(  534)             END DO
(  535)          END DO
(  536)       END DO
(  537)       
(  538)    END SUBROUTINE isl_pth
(  539) 
(  540)    !!----------------------------------------------------------------------
(  541)    !!   Default option :                                        NetCDF file
(  542)    !!----------------------------------------------------------------------
(  543) 
(  544)    SUBROUTINE isl_mat
(  545)       !!----------------------------------------------------------------------
(  546)       !!                  ***  ROUTINE isl_mat  ***
(  547)       !!                
(  548)       !! ** Purpose :   Compute and invert the island matrix 
(  549)       !!
(  550)       !! ** Method  :   aisl(jni,jnj) matrix constituted by bsfisl circulation
(  551)       !!
(  552)       !! ** Action  : - aisl     : island matrix
(  553)       !!              - aislm1   : invert of the island matrix
(  554)       !!      file
(  555)       !!              - islands  : contain bsfisl, aisl and aislm1
(  556)       !!
(  557)       !! History :
(  558)       !!   1.0  !  88-10  (G. Madec)  Original code
(  559)       !!   7.0  !  96-01  (G. Madec)  suppression of common work arrays
(  560)       !!   8.5  !  02-08  (G. Madec)  Free form, F90
(  561)       !!----------------------------------------------------------------------
(  562)       !! * Modules used
(  563)       USE ioipsl
(  564)       
(  565)       !! * Local declarations
(  566)       INTEGER ::   ji, jj, jni, jnj, jn, jl   ! dummy loop indices
(  567)       INTEGER ::   itime, ibvar, ios          ! temporary integers
(  568)       LOGICAL ::   llog






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 11

(  569)       CHARACTER (len=32) ::   clname
(  570)       CHARACTER (len=8 ) ::   clvnames(100)
(  571)       REAL(wp), DIMENSION(1) ::   zdept
(  572)       REAL(wp), DIMENSION(jpi,jpj) ::   zlamt, zphit
(  573)       REAL(wp), DIMENSION(jpi,jpj,2) ::   zwx
(  574)       REAL(wp), DIMENSION(jpisl*jpisl) ::   ztab
(  575)       !!----------------------------------------------------------------------
(  576) 
(  577) 
(  578)       ! I. Island matrix lecture in numisl (if it exists)
(  579)       ! ==================================
(  580) 
(  581)       ! Lecture
(  582)       zlamt(:,:) = 0.
(  583)       zphit(:,:) = 0.
(  584)       zdept(1)   = 0.
(  585)       itime = 0
(  586)       clvnames="        "
(  587)       clname = 'islands'
(  588)       CALL ioget_vname(numisl, ibvar, clvnames)
(  589)       IF(lwp) WRITE(numout,*) clvnames
(  590)       ios=0
(  591)       DO jn=1,100
(  592)         IF(clvnames(jn) == 'aisl') ios=1
(  593)       END DO
(  594)       IF( ios == 0 ) go to 110 
(  595) 
(  596)       CALL restget( numisl, 'aisl'  , jpisl, jpisl, 1, 0, llog, aisl   )
(  597)       CALL restget( numisl, 'aislm1', jpisl, jpisl, 1, 0, llog, aislm1 )
(  598)       CALL restclo( numisl )
(  599)       ! Control print
(  600)       IF(lwp) THEN
(  601)          WRITE(numout,*)
(  602)          WRITE(numout,*)' islmat: lecture aisl/aislm1 in numisl done'
(  603)          WRITE(numout,*)' ~~~~~~'
(  604)          WRITE(numout,*)
(  605)          WRITE(numout,*) '        island matrix : '
(  606)          WRITE(numout,*)
(  607)          
(  608)          DO jnj = 1, jpisl
(  609)             WRITE(numout,'(8e12.4)') ( aisl(jni,jnj), jni = 1, jpisl )
(  610)          END DO
(  611) 
(  612)          WRITE(numout,*)
(  613)          WRITE(numout,*) '       inverse of the island matrix'
(  614)          WRITE(numout,*)
(  615) 
(  616)          DO jnj = 1, jpisl
(  617)             WRITE(numout,'(12e11.3)') ( aislm1(jni,jnj), jni=1,jpisl )
(  618)          END DO
(  619)       ENDIF
(  620)       
(  621)       RETURN
(  622) 
(  623)  110  CONTINUE
(  624) 
(  625) 
(  626)       ! II. Island matrix computation






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 12

(  627)       ! =============================
(  628) 
(  629)       DO jnj = 1, jpisl
(  630) 
(  631)          ! Circulation of bsf(jnj) around island jni
(  632)          
(  633)          DO jj = 2, jpj
(  634)                zwx(:,jj,1) = -( bsfisl(:,jj,jnj) - bsfisl(:,jj-1,jnj) )
(  635)          END DO
(  636)          zwx(:,1,1) = 0.e0
(  637)          
(  638)          DO jj = 1, jpj
(  639)             DO ji = 2, jpi
(  640)                zwx(ji,jj,2) = ( bsfisl(ji,jj,jnj) - bsfisl(ji-1,jj,jnj) )
(  641)             END DO
(  642)          END DO
(  643)          zwx(1,:,2) = 0.e0
(  644)          
(  645)          ! Island matrix
(  646)          
(  647)          DO jni = 1, jpisl
(  648)             aisl(jni,jnj) = 0.e0
(  649)             DO jl = 1, 2
(  650)                DO jj=1,jpj
(  651)                   DO ji=1,jpi
(  652)                      aisl(jni,jnj) = aisl(jni,jnj) + acisl2(ji,jj,jl,jni)*zwx(ji,jj,jl)
(  653)                   END DO
(  654)                END DO
(  655)             END DO
(  656)          END DO
(  657)          
(  658)       END DO
(  659)       IF( lk_mpp ) THEN
(  660)          DO jnj = 1, jpisl
(  661)             DO jni = 1, jpisl
(  662)                ztab(jni+(jnj-1)*jpisl) = aisl(jni,jnj)
(  663)             END DO
(  664)          END DO
(  665)          CALL mpp_sum( ztab, jpisl*jpisl )   ! sum over the global domain
(  666)       ENDIF
(  667) 
(  668)       ! 1.3 Control print
(  669) 
(  670)       IF(lwp) THEN
(  671)          WRITE(numout,*)
(  672)          WRITE(numout,*) 'islmat : island matrix'
(  673)          WRITE(numout,*) '~~~~~~'
(  674)          WRITE(numout,*)
(  675)          
(  676)          DO jnj = 1, jpisl
(  677)             WRITE(numout,'(8e12.4)') ( aisl(jni,jnj), jni = 1, jpisl )
(  678)          END DO
(  679)       ENDIF
(  680)       
(  681) 
(  682)       ! 2. Invertion of the island matrix
(  683)       ! ---------------------------------
(  684)       






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 13

(  685)       ! 2.1 Call of an imsl routine for the matrix invertion
(  686) 
(  687)       CALL linrg( jpisl, aisl, jpisl, aislm1, jpisl )
(  688) 
(  689)       ! 2.2 Control print
(  690)       
(  691)       IF(lwp) THEN
(  692)          WRITE(numout,*)
(  693)          WRITE(numout,*) 'islmat : inverse of the island matrix'
(  694)          WRITE(numout,*) '~~~~~~'
(  695)          WRITE(numout,*)
(  696)          
(  697)          DO jnj = 1, jpisl
(  698)             WRITE(numout, '(12e11.3)')  '        ', ( aislm1(jni,jnj), jni=1, jpisl )
(  699)          END DO
(  700)       ENDIF
(  701)       
(  702) 
(  703)       ! 3. Output of aisl and aislm1 in numisl
(  704)       ! --------------------------------------
(  705) 
(  706)       CALL restput( numisl, 'aisl'  , jpisl, jpisl, 1, 0, aisl   )
(  707)       CALL restput( numisl, 'aislm1', jpisl, jpisl, 1, 0, aislm1 )
(  708)       CALL restclo( numisl )
(  709) 
(  710)    END SUBROUTINE isl_mat
(  711) 
(  712) 
(  713)    SUBROUTINE isl_bsf
(  714)       !!----------------------------------------------------------------------
(  715)       !!                  ***  ROUTINE isl_bsf  ***
(  716)       !!          
(  717)       !! ** Purpose :
(  718)       !!         Compute the barotropic stream function associated with each
(  719)       !!      island using FETI or preconditioned conjugate gradient method
(  720)       !!      or read them in numisl.
(  721)       !!
(  722)       !! ** Method :
(  723)       !!
(  724)       !! ** input/output file :
(  725)       !!            numisl            : barotropic stream function associated
(  726)       !!                          with each island of the domain
(  727)       !!
(  728)       !! ** Action :
(  729)       !!       bsfisl, the streamfunction which takes the value 1 over island ni,
(  730)       !!    and 0 over the others islands
(  731)       !!       file 'numisl' barotropic stream function associated
(  732)       !!                              with each island of the domain
(  733)       !!
(  734)       !! History :
(  735)       !!        !  87-10  (G. Madec)  Original code
(  736)       !!        !  91-11  (G. Madec)
(  737)       !!        !  93-03  (G. Madec)  release 7.1
(  738)       !!        !  93-04  (M. Guyon)  loops and suppress pointers
(  739)       !!        !  96-11  (A. Weaver)  correction to preconditioning
(  740)       !!        !  98-02  (M. Guyon)  FETI method
(  741)       !!        !  99-11  (M. Imbard)  NetCDF FORMAT with IOIPSL
(  742)       !!   8.5  !  02-08  (G. Madec)  Free form, F90






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 14

(  743)       !!----------------------------------------------------------------------
(  744)       !! * Modules used
(  745)       USE ioipsl
(  746)       USE solpcg
(  747)       USE solfet
(  748)       USE solsor
(  749) 
(  750)       !! * Local declarations
(  751)       LOGICAL  ::   llog, llbon
(  752)       CHARACTER (len=10) ::  clisl
(  753)       CHARACTER (len=32) ::  clname, clname2
(  754)       INTEGER  ::   ji, jj, jni, jii, jnp, je   ! dummy loop indices
(  755)       INTEGER  ::   iimlu, ijmlu, inmlu, iju
(  756)       INTEGER  ::   ii, ij, icile, icut, inmax, indic
(  757)       INTEGER  ::   itime, ie
(  758)       REAL(wp) ::   zepsr, zeplu, zgwgt
(  759)       REAL(wp) ::   zep(jpisl), zlamt(jpi,jpj), zphit(jpi,jpj), zdept(1), zprec(4)
(  760)       REAL(wp) ::   zdate0, zdt
(  761)       REAL(wp) ::   t2p1(jpi,1,1)
(  762)       INTEGER  ::   iloc
(  763)       !!----------------------------------------------------------------------
(  764) 
(  765) 
(  766)       ! 0. Initializations
(  767)       ! ==================
(  768) 
(  769)       icile         = 0      ! set to zero the convergence indicator
(  770)       bsfisl(:,:,:) = 0.e0   ! set to zero of bsfisl
(  771) 
(  772) 
(  773)       ! I. Lecture of bsfisl in numisl (if it exists)
(  774)       ! =============================================
(  775) 
(  776)       icut  = 0
(  777)       iimlu = 0
(  778)       ijmlu = 0
(  779)       inmlu = 0
(  780)       zeplu = 0.
(  781)       zlamt(:,:) = 0.
(  782)       zphit(:,:) = 0.
(  783)       zdept(1)   = 0.
(  784)       itime = 0
(  785)       clname = 'islands'
(  786)       ie=1
(  787)       DO je = 1, 32
(  788)         IF( clname(je:je) /= ' ' ) ie = je
(  789)       END DO
(  790)       clname2 = clname(1:ie)//".nc"
(  791)       INQUIRE( FILE=clname2, EXIST=llbon )
(  792) ! islands FILE does not EXIST : icut=999
(  793)       IF( llbon ) THEN 
(  794)          ! island FILE is present 
(  795)          CALL restini(clname,jpi,jpj,zlamt,zphit,1,zdept,  &
(  796)             &         'NONE',itime,zdate0,zdt,numisl,domain_id=nidom)
(  797)          CALL restget(numisl,'PRECISION',1,1,4,0,llog,zprec)
(  798)          iimlu = NINT( zprec(1) )
(  799)          ijmlu = NINT( zprec(2) )
(  800)          inmlu = NINT( zprec(3) )






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 15

(  801)          zeplu = zprec(4)
(  802)          ! the read domain does not correspond to the model one : icut=999
(  803)          IF( iimlu /= jpi .OR. ijmlu /= jpj .OR. inmlu /= jpisl ) THEN
(  804)             icut = 999
(  805)             CALL restclo(numisl)
(  806)          ELSE 
(  807)             DO jni = 1, jpisl
(  808)                IF( jni < 10 ) THEN
(  809)                   WRITE(clisl,'("island",I1)') jni
(  810)                ELSEIF( jni < 100 ) THEN
(  811)                   WRITE(clisl,'("island",I2)') jni
(  812)                ELSE 
(  813)                   WRITE(clisl,'("island",I3)') jni
(  814)                ENDIF
(  815)                CALL restget(numisl,clisl,jpi,jpj,1,0,llog, bsfisl(:,:,jni))
(  816)             END DO
(  817)          ENDIF
(  818)       ELSE 
(  819)          ! islands FILE does not EXIST : icut=999 
(  820)          icut = 999
(  821)          CALL restclo(numisl)
(  822)       ENDIF
(  823)       
(  824)       ! the read precision is not the required one : icut=888
(  825)       IF( zeplu > epsisl ) THEN 
(  826)          icut = 888
(  827)          CALL restclo(numisl)
(  828)       ENDIF
(  829) 
(  830)       ! Control print
(  831)       IF( icut == 999 ) THEN
(  832)           IF(lwp) THEN
(  833)               WRITE(numout,*)
(  834)               WRITE(numout,*) 'islbsf : lecture bsfisl in numisl failed'
(  835)               WRITE(numout,*) '~~~~~~'
(  836)               WRITE(numout,*) '         icut= ', icut
(  837)               WRITE(numout,*) '         imlu= ', iimlu, ' jmlu= ', ijmlu,   &
(  838)                                       ' nilu= ', inmlu, ' epsisl lu= ', zeplu
(  839)               WRITE(numout,*) '         the bsfisl are computed from zero'
(  840)           ENDIF
(  841)         ELSEIF( icut == 888 ) THEN
(  842)           IF(lwp) THEN
(  843)               WRITE(numout,*)
(  844)               WRITE(numout,*) 'islbsf : lecture bsfisl in numisl done'
(  845)               WRITE(numout,*) '~~~~~~'
(  846)               WRITE(numout,*) '         the required accuracy is not reached'
(  847)               WRITE(numout,*) '         epsisl lu= ', zeplu,' epsisl required ', epsisl
(  848)               WRITE(numout,*) '         the bsfisl are computed from the read values'
(  849)           ENDIF
(  850)         ELSE
(  851)           IF(lwp) THEN
(  852)               WRITE(numout,*)
(  853)               WRITE(numout,*) 'islbsf : lecture bsfisl in numisl done'
(  854)               WRITE(numout,*) '~~~~~~'
(  855)           ENDIF
(  856)           RETURN
(  857)       ENDIF
(  858)       






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 16

(  859)       
(  860)       ! II. Compute the bsfisl (if icut=888 or 999)
(  861)       ! ============================================
(  862) 
(  863)       ! save nmax
(  864)       inmax = nmax
(  865) 
(  866)       ! set the number of iteration of island computation
(  867)       nmax = nmisl
(  868) 
(  869)       ! Loop over islands
(  870)       ! -----------------
(  871)       
(  872)       DO jni = 1, jpisl
(  873) 
(  874) 
(  875)          ! 1. Initalizations of island computation
(  876)          ! ---------------------------------------
(  877)          
(  878)          ! Set the pcg solution gcb to zero
(  879)          gcb(:,:) = 0.e0
(  880) 
(  881)          ! Set first guess gcx either to zero or to the read bsfisl
(  882)          IF( icut == 999 ) THEN
(  883)             gcx(:,:) = 0.e0
(  884)          ELSEIF( icut == 888 ) THEN
(  885)             IF(lwp) WRITE(numout,*) ' islbsf: bsfisl read are used as first guess'
(  886)             ! c a u t i o n: bsfisl masked because it contains 1 along island seaside
(  887)             gcx(:,:) = bsfisl(:,:,jni) * bmask(:,:)
(  888)          ENDIF
(  889)          
(  890)          ! Right hand side of the streamfunction equation
(  891)          
(  892)          IF( lk_mpp ) THEN
(  893) 
(  894)             ! north fold treatment
(  895)             IF( npolj == 3 .OR. npolj == 5)   iloc=jpiglo-(nimpp-1+nimppt(nono+1)-1)
(  896)             IF( npolj == 4 .OR. npolj == 6)   iloc=jpiglo-2*(nimpp-1)
(  897)             t2p1(:,1,1) = 0.e0
(  898)             ! north and south grid-points
(  899)             DO jii = 1, 2
(  900)                DO jnp = 1, mnisl(jii,jni)
(  901)                   ii = miisl(jnp,jii,jni)
(  902)                   ij = mjisl(jnp,jii,jni)
(  903)                   IF( ( npolj == 3 .OR. npolj == 4 ) .AND.   &
(  904)                      ( ij == nlcj-1 .AND. jii == 1) ) THEN 
(  905)                      iju=iloc-ii+1
(  906)                      t2p1(iju,1,1) =  t2p1(iju,1,1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
(  907)                   ELSEIF( ( npolj == 5 .OR. npolj == 6 ) .AND.   &
(  908)                      ( ij == nlcj-1 .AND. jii == 1) ) THEN 
(  909)                      iju=iloc-ii
(  910)                      gcb(ii,ij) =  gcb(ii,ij) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
(  911)                      t2p1(iju,1,1) =  t2p1(iju,1,1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
(  912)                   ELSE  
(  913)                      gcb(ii,ij-jii+1) =  gcb(ii,ij-jii+1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
(  914)                   ENDIF
(  915)                END DO
(  916)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 17

(  917)          
(  918)             ! east and west grid-points
(  919)          
(  920)             DO jii = 3, 4
(  921)                DO jnp = 1, mnisl(jii,jni)
(  922)                   ii = miisl(jnp,jii,jni)
(  923)                   ij = mjisl(jnp,jii,jni)
(  924)                   gcb(ii-jii+3,ij) = gcb(ii-jii+3,ij) + hvr(ii,ij) * e2v(ii,ij) / e1v(ii,ij)
(  925)                END DO
(  926)             END DO
(  927) 
(  928)             IF( lk_mpp )   CALL mpplnks( gcb )   !!bug ? should use an lbclnk ? is it possible?
(  929) 
(  930)          ELSE
(  931) 
(  932)             ! north and south grid-points
(  933)             DO jii = 1, 2
(  934)                DO jnp = 1, mnisl(jii,jni)
(  935)                   ii = miisl(jnp,jii,jni)
(  936)                   ij = mjisl(jnp,jii,jni)
(  937)                   IF( ( nperio == 3 .OR. nperio == 4 ) .AND.   &
(  938)                      ( ij == jpj-1 .AND. jii == 1) ) THEN 
(  939)                      gcb(jpi-ii+1,ij-1) = gcb(jpi-ii+1,ij-1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
(  940)                   ELSEIF( ( nperio == 5 .OR. nperio == 6 ) .AND.   &
(  941)                      ( ij == jpj-1 .AND. jii == 1) ) THEN 
(  942)                      gcb(ii,ij) =  gcb(ii,ij) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij)
(  943)                      gcb(jpi-ii,ij) = gcb(jpi-ii,ij) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij) 
(  944)                   ELSE  
(  945)                      gcb(ii,ij-jii+1) =  gcb(ii,ij-jii+1) + hur(ii,ij) * e1u(ii,ij) / e2u(ii,ij)
(  946)                   ENDIF
(  947)                END DO
(  948)             END DO
(  949) 
(  950)             ! east and west grid-points
(  951)             DO jii = 3, 4
(  952)                DO jnp = 1, mnisl(jii,jni)
(  953)                   ii = miisl(jnp,jii,jni)
(  954)                   ij = mjisl(jnp,jii,jni)
(  955)                   IF( bmask(ii-jii+3,ij) /= 0. ) THEN
(  956)                      gcb(ii-jii+3,ij) = gcb(ii-jii+3,ij) + hvr(ii,ij) * e2v(ii,ij) / e1v(ii,ij)
(  957)                   ELSE
(  958)                      ! east-west cyclic boundary conditions
(  959)                      IF( ii-jii+3 == 1 ) THEN
(  960)                         gcb(jpim1,ij) = gcb(jpim1,ij) + hvr(ii,ij) * e2v(ii,ij) / e1v(ii,ij)
(  961)                      ENDIF
(  962)                   ENDIF
(  963)                END DO
(  964)             END DO
(  965)          ENDIF
(  966) 
(  967)          ! Preconditioned right hand side and absolute precision
(  968) 
(  969)          IF( nsolv == 3 ) THEN 
(  970)             ! FETI method
(  971)             ncut    = 0
(  972)             rnorme  = 0.e0
(  973)             gcb(:,:) = bmask(:,:) * gcb(:,:)
(  974)             DO jj = 1, jpj






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 18

(  975)                DO ji = 1, jpi
(  976)                   rnorme = rnorme + gcb(ji,jj) * gcb(ji,jj)
(  977)                END DO
(  978)             END DO
(  979)             
(  980)             IF( lk_mpp )   CALL mpp_sum( rnorme )
(  981) 
(  982)             IF(lwp) WRITE(numout,*) 'rnorme ', rnorme
(  983)             epsr  = epsisl * epsisl * rnorme
(  984)             indic = 0
(  985)          ELSE 
(  986)             ncut   = 0
(  987)             rnorme = 0.e0
(  988)             DO jj = 1, jpj
(  989)                DO ji = 1, jpi
(  990)                   gcb  (ji,jj) = gcdprc(ji,jj) * gcb(ji,jj)
(  991)                   zgwgt  = gcdmat(ji,jj) * gcb(ji,jj)
(  992)                   rnorme = rnorme + gcb(ji,jj) * zgwgt
(  993)                END DO
(  994)             END DO
(  995)             IF( lk_mpp )   CALL mpp_sum( rnorme )   ! sum over the global domain
(  996) 
(  997)             IF(lwp) WRITE(numout,*) 'rnorme ', rnorme
(  998)             epsr = epsisl * epsisl * rnorme
(  999)             indic = 0
( 1000)          ENDIF
( 1001) 
( 1002) 
( 1003)          ! 3. PCG solver for gcp.gcx=gcb in monotask
( 1004)          ! -------------------------------------------
( 1005) 
( 1006)          IF( nsolv == 3 ) THEN 
( 1007)             epsilo = epsisl       ! precision to compute Islands matrix A
( 1008)             CALL sol_fet( indic )  ! FETI method
( 1009)             epsilo = eps          ! precision to compute grad PS
( 1010)          ELSE 
( 1011)             CALL sol_pcg( indic )  ! pcg method
( 1012)          ENDIF
( 1013) 
( 1014) 
( 1015)          ! 4. Save the solution in bsfisl
( 1016)          ! ------------------------------
( 1017) 
( 1018)          bsfisl(:,:,jni) = gcx(:,:)
( 1019) 
( 1020) 
( 1021)          ! 5. Boundary conditions
( 1022)          ! ----------------------
( 1023) 
( 1024)          ! set to 1. coastal gridpoints of the island
( 1025)          DO jnp = 1, mnisl(0,jni)
( 1026)             ii = miisl(jnp,0,jni)
( 1027)             ij = mjisl(jnp,0,jni)
( 1028)             bsfisl(ii,ij,jni) = 1.
( 1029)          END DO
( 1030) 
( 1031)          ! cyclic boundary conditions
( 1032)          IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 19

( 1033)             bsfisl( 1 ,:,jni) = bsfisl(jpim1,:,jni)
( 1034)             bsfisl(jpi,:,jni) = bsfisl(  2  ,:,jni)
( 1035)          ENDIF
( 1036)          IF( nperio == 3 .OR. nperio == 4 ) THEN
( 1037)             DO ji = 1, jpim1
( 1038)                iju = jpi-ji+1
( 1039)                bsfisl(ji,jpj-1,jni) = bsfisl(iju,jpj-2,jni)
( 1040)                bsfisl(ji, jpj ,jni) = bsfisl(iju,jpj-3,jni)
( 1041)             END DO
( 1042)          ENDIF
( 1043)          IF( nperio == 5 .OR. nperio == 6 ) THEN
( 1044)             DO ji = 1, jpi-1
( 1045)                iju=jpi-ji
( 1046)                bsfisl(ji,jpj,jni) = bsfisl(iju,jpj-2,jni)
( 1047)             END DO
( 1048)             DO ji = jpi/2+1, jpi-1
( 1049)                iju=jpi-ji
( 1050)                bsfisl (ji,jpjm1,jni) = bsfisl(iju,jpjm1,jni)
( 1051)             END DO
( 1052)          ENDIF
( 1053)          IF( lk_mpp )   CALL lbc_lnk( bsfisl(:,:,jni), 'G', 1. )   ! link at G-point
( 1054) 
( 1055) 
( 1056)          ! 6. Control print
( 1057)          ! ----------------
( 1058) 
( 1059)          IF(lwp) WRITE(numout,*)
( 1060)          IF(lwp) WRITE(numout,*) ' islbsf: island number: ', jni
( 1061)          IF(lwp) WRITE (numout,9290) niter, res, SQRT(epsr)/epsisl
( 1062)          zep(jni) =  MAX(epsisl, res/(SQRT(epsr)/epsisl))
( 1063)          IF( indic <  0 ) THEN
( 1064)             icile = icile-1
( 1065)             IF(lwp) WRITE(numout,*) '  pcg do not converge for island: ', jni
( 1066)             IF(lwp) WRITE(numout,*) '      Precision reached: ',zep(jni)
( 1067)          ENDIF
( 1068) 
( 1069) 9290     FORMAT('        niter :',i4,' , res :',e20.10,' , gcb :',e20.10)
( 1070) 
( 1071)          !                                          !====================
( 1072)       END DO                                        !  End Loop islands
( 1073)       !                                             !====================
( 1074) 
( 1075)       ! 7. Reset PCG
( 1076)       ! ------------
( 1077) 
( 1078)       ! reset the number of iteration for pcg
( 1079)       nmax = inmax
( 1080) 
( 1081)       ! reset to zero pcg arrays
( 1082)       gcx  (:,:) = 0.e0
( 1083)       gcxb (:,:) = 0.e0
( 1084)       gcb  (:,:) = 0.e0
( 1085)       gcr  (:,:) = 0.e0
( 1086)       gcdes(:,:) = 0.e0
( 1087)       gccd (:,:) = 0.e0
( 1088) 
( 1089) 
( 1090)       ! III. Output of bsfisl in numisl






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 20

( 1091)       ! ===============================
( 1092) 
( 1093)       CALL ymds2ju( 0, 1, 1, 0.e0, zdate0 )
( 1094)       zprec(1) = FLOAT(jpi)
( 1095)       zprec(2) = FLOAT(jpj)
( 1096)       zprec(3) = FLOAT(jpisl)
( 1097)       IF(lwp) WRITE(numout,*) clname
( 1098)       CALL restini( 'NONE', jpi, jpj, glamt, gphit, 1, zdept,  &
( 1099)          &          clname, itime, zdate0, rdt, numisl, domain_id=nidom )
( 1100)       IF( icile == 0 .AND. icut /= 0 ) THEN
( 1101)          IF(lwp) THEN
( 1102)             WRITE(numout,*)
( 1103)             WRITE(numout,*)' islbsf: write bsfisl in numisl ', numisl
( 1104)             WRITE(numout,*)' -------------------------------------'
( 1105)          ENDIF
( 1106)          zprec(4) = epsisl
( 1107)          CALL restput(numisl,'PRECISION',1,1,4,0,zprec)
( 1108)          DO jni = 1, jpisl
( 1109)             IF(jni < 10) THEN
( 1110)                WRITE(clisl,'("island",I1)') jni
( 1111)             ELSE IF(jni < 100) THEN
( 1112)                WRITE(clisl,'("island",I2)') jni
( 1113)             ELSE 
( 1114)                WRITE(clisl,'("island",I3)') jni
( 1115)             ENDIF
( 1116)             CALL restput( numisl, clisl, jpi, jpj, 1, 0, bsfisl(:,:,jni) )
( 1117)          END DO
( 1118)       ENDIF
( 1119) 
( 1120)       IF( icile < 0 ) THEN
( 1121)          IF(lwp) THEN
( 1122)             WRITE(numout,*)
( 1123)             WRITE(numout,*) ' islbsf: number of island without convergence : ',ABS(icile)
( 1124)             WRITE(numout,*) ' ---------------------------------------------'
( 1125)          ENDIF
( 1126)          zepsr = epsisl
( 1127)          DO jni = 1, jpisl
( 1128)             IF(lwp) WRITE(numout,*) '    isl ',jni,' precision reached ', zep(jni)
( 1129)             zepsr = MAX( zep(jni), zepsr )
( 1130)          END DO
( 1131)          IF( zepsr == 0. ) zepsr = epsisl
( 1132)          IF(lwp) THEN
( 1133)             WRITE(numout,*) ' save value of precision reached: ',zepsr
( 1134)             WRITE(numout,*)
( 1135)             WRITE(numout,*)' islbsf: save bsfisl in numisl ',numisl
( 1136)             WRITE(numout,*)' -------------------------------------'
( 1137)          ENDIF
( 1138) 
( 1139)          zprec(4) = zepsr
( 1140)          CALL restput( numisl, 'PRECISION', 1, 1, 1, 0, zprec )
( 1141)          DO jni = 1, jpisl
( 1142)             IF( jni < 10 ) THEN
( 1143)                WRITE(clisl,'("island",I1)') jni
( 1144)             ELSE IF( jni < 100 ) THEN
( 1145)                WRITE(clisl,'("island",I2)') jni
( 1146)             ELSE 
( 1147)                WRITE(clisl,'("island",I3)') jni
( 1148)             ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 21

( 1149)             CALL restput( numisl, clisl, jpi, jpj, 1, 0, bsfisl(:,:,jni) )
( 1150)          END DO
( 1151)          CALL restclo(numisl)
( 1152)          nstop = nstop + 1
( 1153)       ENDIF
( 1154) 
( 1155)    END SUBROUTINE isl_bsf
( 1156) 
( 1157) 
( 1158)    SUBROUTINE isl_dyn_spg
( 1159)       !!----------------------------------------------------------------------
( 1160)       !!                  ***  routine isl_dyn_spg  ***
( 1161)       !!
( 1162)       !! ** Purpose :   Compute and add the island contribution to the 
( 1163)       !!      barotropic stream function trend.
( 1164)       !!
( 1165)       !!
( 1166)       !! ** Method  :   Rigid-lid appromimation: ...????
( 1167)       !!
( 1168)       !! ** Action : - Update bsfd with the island contribution
( 1169)       !!
( 1170)       !! History :
( 1171)       !!   9.0  !  03-09  (G. Madec)  isolate island computation
( 1172)       !!---------------------------------------------------------------------
( 1173) 
( 1174)       !! * Local declarations
( 1175)       INTEGER ::   ji, jj, jni, jnj    ! dummy loop indices
( 1176)       !!----------------------------------------------------------------------
( 1177)       
( 1178) 
( 1179)       ! compute the island potential
( 1180)       ! ----------------------------
( 1181)       DO jni = 1, jpisl                      ! second member
( 1182)          bisl(jni) =  0.e0
( 1183)          DO jj = 2, jpj
( 1184)             DO ji = 2, jpi
( 1185)                bisl(jni) =  bisl(jni) + acisl1(ji,jj,1,jni) *  spgu(ji,jj)                  &
( 1186)                   &                   + acisl1(ji,jj,2,jni) *  spgv(ji,jj)                  &
( 1187)                   &                   + acisl2(ji,jj,1,jni) * ( gcx(ji,jj)-gcx(ji,jj-1) )   &
( 1188)                   &                   - acisl2(ji,jj,2,jni) * ( gcx(ji,jj)-gcx(ji-1,jj) )
( 1189)             END DO
( 1190)          END DO
( 1191)       END DO
( 1192)       IF( lk_mpp )   CALL mpp_sum( bisl, jpisl )   ! sum over the global domain
( 1193) 
( 1194)       DO jni = 1, jpisl                     ! Island stream function trend
( 1195)          visl(jni) = 0.e0
( 1196)          DO jnj = 1, jpisl
( 1197)             visl(jni) = visl(jni) + aislm1(jni,jnj) * bisl(jnj)
( 1198)          END DO
( 1199)       END DO
( 1200) 
( 1201)       ! update the bsf trend ( caution : bsfd is not zero along island coastlines, dont mask it ! )
( 1202)       ! --------------------
( 1203)       DO jj = 1, jpj
( 1204)          DO jni = 1, jpisl
( 1205)             DO ji = 1, jpi
( 1206)                bsfd(ji,jj) = bsfd(ji,jj) + visl(jni) * bsfisl(ji,jj,jni)






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 22

( 1207)             END DO
( 1208)          END DO
( 1209)       END DO
( 1210) 
( 1211)    END SUBROUTINE isl_dyn_spg
( 1212) 
( 1213) 
( 1214)    SUBROUTINE isl_stp_ctl( kt, kindic )
( 1215)       !!----------------------------------------------------------------------
( 1216)       !!                    ***  ROUTINE isl_stp_ctl  ***
( 1217)       !!                     
( 1218)       !! ** Purpose :   ???
( 1219)       !!
( 1220)       !! ** Method  : - print island potential
( 1221)       !!
( 1222)       !! History :
( 1223)       !!   9.0  !  03-09  (G. Madec)  isolated from stp_ctl
( 1224)       !!----------------------------------------------------------------------
( 1225)       !! * Arguments
( 1226)       INTEGER, INTENT(in   ) ::   kt        ! ocean time-step index
( 1227)       INTEGER, INTENT(inout) ::   kindic    ! indicator of solver convergence
( 1228) 
( 1229)       !! * local declarations
( 1230)       INTEGER  ::   jni                     ! dummy loop indice
( 1231)       REAL(wp) ::   zfact                   ! temporary scalar
( 1232)       !!----------------------------------------------------------------------
( 1233)       !!  OPA 8.5, LODYC-IPSL (2002)
( 1234)       !!----------------------------------------------------------------------
( 1235) 
( 1236)       IF( kt == nit000 .AND. lwp ) THEN
( 1237)          WRITE(numout,*)
( 1238)          WRITE(numout,*) 'isl_stp_ctl : time-stepping control'
( 1239)          WRITE(numout,*) '~~~~~~~~~~~'
( 1240)       ENDIF
( 1241) 
( 1242)       ! Island trends
( 1243)       DO jni = 1, jpisl
( 1244)          zfact = 0.
( 1245)          IF( miisl(1,0,jni) /= 0 .AND. mjisl(1,0,jni) /= 0 ) THEN
( 1246)             zfact = 1.e-6 * bsfn(miisl(1,0,jni),mjisl(1,0,jni))
( 1247)          ENDIF
( 1248)          IF( lk_mpp )   CALL mpp_isl( zfact )
( 1249) 
( 1250)          IF(lwp) WRITE(numisp,9300) kt, jni, zfact, visl(jni)
( 1251)          IF( MOD( kt, nwrite ) == 0 .OR. kindic < 0     &
( 1252)             .OR. ( kt == nit000 .AND. kindic > 0 )      &
( 1253)             .OR.   kt == nitend                         ) THEN
( 1254)             IF( jni == 1 .AND. lwp ) THEN
( 1255)                WRITE(numout,*)
( 1256)                WRITE(numout,*) 'isl_stp_ctl : island bsf'
( 1257)                WRITE(numout,*) '~~~~~~~~~~~'
( 1258)             ENDIF
( 1259)             IF(lwp) WRITE(numout,9300) kt, jni, zfact, visl(jni)
( 1260)          ENDIF
( 1261)       END DO
( 1262) 9300  FORMAT(' it : ',i8,' island  :',i4,'   BSF (Sverdrup) : ',f7.2, '   visl : ',e15.6)
( 1263) 
( 1264)    END SUBROUTINE isl_stp_ctl






PGF90 (Version     10.5)          09/23/2010  09:49:32      page 23

( 1265) 
( 1266) #else
( 1267)    !!----------------------------------------------------------------------
( 1268)    !!   Default option                                         Empty module
( 1269)    !!----------------------------------------------------------------------
( 1270)    LOGICAL, PUBLIC, PARAMETER ::   lk_isl = .FALSE.    !: 'key_islands' flag
( 1271) CONTAINS
( 1272)    SUBROUTINE isl_dom                        ! Empty routine
( 1273)    END SUBROUTINE isl_dom
( 1274)    SUBROUTINE isl_bsf                        ! Empty routine
( 1275)    END SUBROUTINE isl_bsf
( 1276)    SUBROUTINE isl_mat                        ! Empty routine
( 1277)    END SUBROUTINE isl_mat
( 1278)    SUBROUTINE isl_dyn_spg                    ! Empty routine
( 1279)    END SUBROUTINE isl_dyn_spg
( 1280)    SUBROUTINE isl_stp_ctl( kt, kindic )      ! Empty routine
( 1281) !      WRITE(*,*) 'isl_stp_ctl: You should not have seen this print! error?', kt, kindic
( 1282)    END SUBROUTINE isl_stp_ctl
( 1283) #endif
( 1284) 
( 1285)    !!======================================================================
( 1286) END MODULE solisl










































PGF90 (Version     10.5)          09/23/2010  09:49:32      page 24

(    1) # 1 "solisl.F90"
(    1) MODULE solisl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  solisl  ***
(    4)    !! Ocean island : specific treatment of island in rigid-lid case
(    5)    !!==============================================================================
(    6) # 1267
( 1267)    !!----------------------------------------------------------------------
( 1268)    !!   Default option                                         Empty module
( 1269)    !!----------------------------------------------------------------------
( 1270)    LOGICAL, PUBLIC, PARAMETER ::   lk_isl = .FALSE.    !: 'key_islands' flag
( 1271) CONTAINS
( 1272)    SUBROUTINE isl_dom                        ! Empty routine
( 1273)    END SUBROUTINE isl_dom
( 1274)    SUBROUTINE isl_bsf                        ! Empty routine
( 1275)    END SUBROUTINE isl_bsf
( 1276)    SUBROUTINE isl_mat                        ! Empty routine
( 1277)    END SUBROUTINE isl_mat
( 1278)    SUBROUTINE isl_dyn_spg                    ! Empty routine
( 1279)    END SUBROUTINE isl_dyn_spg
( 1280)    SUBROUTINE isl_stp_ctl( kt, kindic )      ! Empty routine
( 1281) !      WRITE(*,*) 'isl_stp_ctl: You should not have seen this print! error?', kt, kindic
( 1282)    END SUBROUTINE isl_stp_ctl
( 1283) 
( 1284) # 1285
( 1285)    !!======================================================================
( 1286) END MODULE solisl
