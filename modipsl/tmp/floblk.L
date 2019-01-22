


PGF90 (Version     10.5)          09/23/2010  09:49:47      page 1

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

Filename: floblk.F90

(    1) !!DB 2008.03.19 ...
(    2) !!Added a Random Displacement Model (RDM) = particle tracking plus turbulence
(    3) !!Routine is called flo_RDM. It replaces the flo_4rk() routine if(ln_flork4)
(    4) !!NB: the Rk4 routine did not work for me.
(    5) !!To use RDM scheme: in namelist set ln_flork4 = .true. 
(    6) !!Requires a init_float file with float positions which I have also modified. 
(    7) !!See flodom.F90 and look for DB
(    8) !!NB: Routine works but it is not really complete. For example:
(    9) !!   (a) it ignores w-direction, i.e. floats remain at initial level
(   10) !!   (b) following from (a): there is no vertical turbulence
(   11) !!   (c) Add more notes ... TO DO ...
(   12) 
(   13) 
(   14) MODULE floblk
(   15)    !!======================================================================
(   16)    !!                     ***  MODULE  floblk  ***
(   17)    !! Ocean floats :   trajectory computation
(   18)    !!======================================================================
(   19) #if   defined key_floats   ||   defined key_esopa
(   20)    !!----------------------------------------------------------------------
(   21)    !!   'key_floats'                                     float trajectories
(   22)    !!----------------------------------------------------------------------
(   23) 
(   24)    !!----------------------------------------------------------------------
(   25)    !!    flotblk     : compute float trajectories with Blanke algorithme
(   26)    !!----------------------------------------------------------------------
(   27)    !! * Modules used
(   28)    USE flo_oce         ! ocean drifting floats
(   29)    USE oce             ! ocean dynamics and tracers
(   30)    USE dom_oce         ! ocean space and time domain
(   31)    USE phycst          ! physical constants
(   32)    USE in_out_manager  ! I/O manager
(   33)    USE lib_mpp         ! distribued memory computing library
(   34) !!DB
(   35)    USE obc_oce         ! ocean open boundary conditions and other definitions
(   36)    USE ldfdyn_oce          ! to get momentum diffusivities ??? 
(   37)    USE ldftra_oce          ! to get tracer diffusivities ??? 
(   38) 
(   39)    IMPLICIT NONE
(   40)    PRIVATE
(   41) 
(   42) ! CN: For Mersenne twister RNG
(   43) ! Default seed
(   44)     integer, parameter :: defaultsd = 4357
(   45) ! Period parameters
(   46)     integer, parameter :: N = 624, N1 = N + 1






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 2

(   47) 
(   48) ! the array for the state vector
(   49)     integer, save, dimension(0:N-1) :: mt
(   50)     integer, save                   :: mti = N1
(   51) 
(   52)    !! * Accessibility
(   53) !!DB
(   54)    PRIVATE gasdev      ! gauss function 
(   55)    PUBLIC flo_blk      ! routine called by floats.F90
(   56)    PUBLIC flo_RDM      ! routine called by floats.F90
(   57)    !!----------------------------------------------------------------------
(   58)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   59)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/floblk.F90,v 1.4 2005/03/27 18:35:05 opalod Exp $ 
(   60)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   61)    !!----------------------------------------------------------------------
(   62) 
(   63) CONTAINS
(   64) 
(   65)    SUBROUTINE flo_blk( kt )
(   66)       !!---------------------------------------------------------------------
(   67)       !!                  ***  ROUTINE flo_blk  ***
(   68)       !!           
(   69)       !! ** Purpose :   Compute the geographical position,latitude, longitude
(   70)       !!      and depth of each float at each time step.
(   71)       !! 
(   72)       !! ** Method  :   The position of a float is computed with Bruno Blanke
(   73)       !!      algorithm. We need to know the velocity field, the old positions
(   74)       !!      of the floats and the grid defined on the domain.
(   75)       !!
(   76)       !!----------------------------------------------------------------------
(   77)       !! * arguments
(   78)       INTEGER, INTENT( in  ) ::   kt ! ocean time step
(   79) 
(   80)       !! * Local declarations
(   81)       INTEGER :: jfl              ! dummy loop arguments
(   82)       INTEGER :: ind, ifin, iloop
(   83) !!DB
(   84) !      INTEGER , DIMENSION ( jpnfl )  ::   &
(   85)       INTEGER , DIMENSION (:),ALLOCATABLE, SAVE  ::   &
(   86)          iil, ijl, ikl,             &     ! index of nearest mesh
(   87)          iiloc , ijloc,             &
(   88)          iiinfl, ijinfl, ikinfl,    &     ! index of input mesh of the float.
(   89)          iioutfl, ijoutfl, ikoutfl        ! index of output mesh of the float.
(   90) !!DB
(   91) !      REAL(wp) , DIMENSION ( jpnfl )  ::    &
(   92)       REAL(wp) , DIMENSION (:),ALLOCATABLE, SAVE  ::    &
(   93)          zgifl, zgjfl, zgkfl,       &     ! position of floats, index on 
(   94)                                           ! velocity mesh.
(   95)          ztxfl, ztyfl, ztzfl,       &     ! time for a float to quit the mesh
(   96)                                           ! across one of the face x,y and z 
(   97)          zttfl,                     &     ! time for a float to quit the mesh 
(   98)          zagefl,                    &     ! time during which, trajectorie of 
(   99)                                           ! the float has been computed
(  100)          zagenewfl,                 &     ! new age of float after calculation 
(  101)                                           ! of new position
(  102)          zufl, zvfl, zwfl,          &     ! interpolated vel. at float position
(  103)          zudfl, zvdfl, zwdfl,       &     ! velocity diff input/output of mesh
(  104)          zgidfl, zgjdfl, zgkdfl           ! direction index of float 






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 3

(  105)       REAL(wp)   ::       &
(  106)          zuinfl,zvinfl,zwinfl,      &     ! transport across the input face
(  107)          zuoutfl,zvoutfl,zwoutfl,   &     ! transport across the ouput face
(  108)          zvol,                      &     ! volume of the mesh
(  109)          zsurfz,                    &     ! surface of the face of the mesh 
(  110)          zind
(  111)       REAL(wp), DIMENSION ( 2 )  ::   &
(  112)          zsurfx, zsurfy                   ! surface of the face of the mesh
(  113) 
(  114) 
(  115)       !!---------------------------------------------------------------------
(  116)       
(  117)       IF( kt == nit000 ) THEN
(  118)          IF(lwp) WRITE(numout,*)
(  119)          IF(lwp) WRITE(numout,*) 'flo_blk : compute Blanke trajectories for floats '
(  120)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(  121) !!DB Allocate arrays
(  122)          ALLOCATE(iil(jpnfl));ALLOCATE(ijl(jpnfl));ALLOCATE(ikl(jpnfl))
(  123)          ALLOCATE(iiloc(jpnfl));ALLOCATE(ijloc(jpnfl))
(  124)          ALLOCATE(iiinfl(jpnfl));ALLOCATE(ijinfl(jpnfl));ALLOCATE(ikinfl(jpnfl))
(  125)          ALLOCATE(iioutfl(jpnfl));ALLOCATE(ijoutfl(jpnfl));ALLOCATE(ikoutfl(jpnfl))
(  126) 
(  127)          ALLOCATE(zgifl(jpnfl));ALLOCATE(zgjfl(jpnfl));ALLOCATE(zgkfl(jpnfl))
(  128)          ALLOCATE(zgidfl(jpnfl));ALLOCATE(zgjdfl(jpnfl));ALLOCATE(zgkdfl(jpnfl))
(  129)          ALLOCATE(ztxfl(jpnfl));ALLOCATE(ztyfl(jpnfl));ALLOCATE(ztzfl(jpnfl))
(  130)          ALLOCATE(zudfl(jpnfl));ALLOCATE(zvdfl(jpnfl));ALLOCATE(zwdfl(jpnfl))
(  131)          ALLOCATE(zufl(jpnfl));ALLOCATE(zvfl(jpnfl));ALLOCATE(zwfl(jpnfl))
(  132)          ALLOCATE(zttfl(jpnfl));ALLOCATE(zagefl(jpnfl));ALLOCATE(zagenewfl(jpnfl))
(  133)       ENDIF
(  134) 
(  135) 
(  136)       ! Initialisation of parameters
(  137)       
(  138)       DO jfl = 1, jpnfl
(  139)          ! ages of floats are put at zero
(  140)          zagefl(jfl) = 0.
(  141)          ! index on the velocity grid 
(  142)          ! We considere k coordinate negative, with this transformation 
(  143)          ! the computation in the 3 direction is the same. 
(  144)          zgifl(jfl) = tpifl(jfl) - 0.5
(  145)          zgjfl(jfl) = tpjfl(jfl) - 0.5
(  146)          zgkfl(jfl) = MIN(-1.,-(tpkfl(jfl)))
(  147)          ! surface drift every 10 days 
(  148)          IF( ln_argo ) THEN
(  149)             IF( MOD(kt,150) >= 146 .OR. MOD(kt,150) == 0 )  zgkfl(jfl) = -1.
(  150)          ENDIF
(  151)          ! index of T mesh
(  152)          iil(jfl) = 1 + INT(zgifl(jfl))
(  153)          ijl(jfl) = 1 + INT(zgjfl(jfl))
(  154)          ikl(jfl) =     INT(zgkfl(jfl))
(  155)       END DO
(  156)        
(  157)       iloop = 0
(  158) 222   DO jfl = 1, jpnfl
(  159) # if   defined key_mpp_mpi   ||   defined key_mpp_shmem
(  160)          IF( (iil(jfl) >= (mig(nldi)-jpizoom+1)) .AND. (iil(jfl) <= (mig(nlei)-jpizoom+1)) .AND.   &
(  161)              (ijl(jfl) >= (mjg(nldj)-jpjzoom+1)) .AND. (ijl(jfl) <= (mjg(nlej)-jpjzoom+1)) ) THEN
(  162)             iiloc(jfl) = iil(jfl) - (mig(1)-jpizoom+1) + 1






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 4

(  163)             ijloc(jfl) = ijl(jfl) - (mjg(1)-jpjzoom+1) + 1
(  164) # else 
(  165)             iiloc(jfl) = iil(jfl)
(  166)             ijloc(jfl) = ijl(jfl)
(  167) # endif
(  168)             
(  169)             ! compute the transport across the mesh where the float is.            
(  170)             zsurfx(1) = e2u(iiloc(jfl)-1,ijloc(jfl)  ) * e3t(-ikl(jfl))
(  171)             zsurfx(2) = e2u(iiloc(jfl)  ,ijloc(jfl)  ) * e3t(-ikl(jfl))
(  172)             zsurfy(1) = e1v(iiloc(jfl)  ,ijloc(jfl)-1) * e3t(-ikl(jfl))
(  173)             zsurfy(2) = e1v(iiloc(jfl)  ,ijloc(jfl)  ) * e3t(-ikl(jfl))
(  174) 
(  175)             ! for a isobar float zsurfz is put to zero. The vertical velocity will be zero too.
(  176)             zsurfz=  e1t(iiloc(jfl),ijloc(jfl)) * e2t(iiloc(jfl),ijloc(jfl))
(  177)             zvol  =( e1t(iiloc(jfl),ijloc(jfl)) * e2t(iiloc(jfl),ijloc(jfl)) * e3t(-ikl(jfl)) )
(  178) 
(  179)             !
(  180)             zuinfl =( ub(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) )/2.*zsurfx(1)
(  181)             zuoutfl=( ub(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) )/2.*zsurfx(2)
(  182)             zvinfl =( vb(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) )/2.*zsurfy(1)
(  183)             zvoutfl=( vb(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) )/2.*zsurfy(2)
(  184)             zwinfl =-(wb(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1))    &
(  185)                &   +  wn(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1)) )/2. *  zsurfz*nisobfl(jfl)
(  186)             zwoutfl=-(wb(iiloc(jfl),ijloc(jfl),- ikl(jfl)   )   &
(  187)                &   +  wn(iiloc(jfl),ijloc(jfl),- ikl(jfl)   ) )/2. *  zsurfz*nisobfl(jfl)
(  188)             
(  189)             ! interpolation of velocity field on the float initial position            
(  190)             zufl(jfl)=  zuinfl  + ( zgifl(jfl) - float(iil(jfl)-1) ) * ( zuoutfl - zuinfl)
(  191)             zvfl(jfl)=  zvinfl  + ( zgjfl(jfl) - float(ijl(jfl)-1) ) * ( zvoutfl - zvinfl)
(  192)             zwfl(jfl)=  zwinfl  + ( zgkfl(jfl) - float(ikl(jfl)-1) ) * ( zwoutfl - zwinfl)
(  193)             
(  194)             ! faces of input and output
(  195)             ! u-direction
(  196)             IF( zufl(jfl) < 0. ) THEN
(  197)                iioutfl(jfl) = iil(jfl) - 1.
(  198)                iiinfl (jfl) = iil(jfl)
(  199)                zind   = zuinfl
(  200)                zuinfl = zuoutfl
(  201)                zuoutfl= zind
(  202)             ELSE
(  203)                iioutfl(jfl) = iil(jfl)
(  204)                iiinfl (jfl) = iil(jfl) - 1
(  205)             ENDIF
(  206)             ! v-direction       
(  207)             IF( zvfl(jfl) < 0. ) THEN
(  208)                ijoutfl(jfl) = ijl(jfl) - 1.
(  209)                ijinfl (jfl) = ijl(jfl)
(  210)                zind    = zvinfl
(  211)                zvinfl  = zvoutfl
(  212)                zvoutfl = zind
(  213)             ELSE
(  214)                ijoutfl(jfl) = ijl(jfl)
(  215)                ijinfl (jfl) = ijl(jfl) - 1.
(  216)             ENDIF
(  217)             ! w-direction
(  218)             IF( zwfl(jfl) < 0. ) THEN
(  219)                ikoutfl(jfl) = ikl(jfl) - 1.
(  220)                ikinfl (jfl) = ikl(jfl)






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 5

(  221)                zind    = zwinfl
(  222)                zwinfl  = zwoutfl
(  223)                zwoutfl = zind
(  224)             ELSE
(  225)                ikoutfl(jfl) = ikl(jfl)
(  226)                ikinfl (jfl) = ikl(jfl) - 1.
(  227)             ENDIF
(  228)             
(  229)             ! compute the time to go out the mesh across a face
(  230)             ! u-direction
(  231)             zudfl (jfl) = zuoutfl - zuinfl
(  232)             zgidfl(jfl) = float(iioutfl(jfl) - iiinfl(jfl))
(  233)             IF( zufl(jfl)*zuoutfl <= 0. ) THEN
(  234)                ztxfl(jfl) = 1.E99
(  235)             ELSE
(  236)                IF( ABS(zudfl(jfl)) >= 1.E-5 ) THEN
(  237)                   ztxfl(jfl)= zgidfl(jfl)/zudfl(jfl) * LOG(zuoutfl/zufl (jfl))
(  238)                ELSE
(  239)                   ztxfl(jfl)=(float(iioutfl(jfl))-zgifl(jfl))/zufl(jfl)
(  240)                ENDIF
(  241)                IF( (ABS(zgifl(jfl)-float(iiinfl (jfl))) <=  1.E-7) .OR.   &
(  242)                    (ABS(zgifl(jfl)-float(iioutfl(jfl))) <=  1.E-7) ) THEN
(  243)                   ztxfl(jfl)=(zgidfl(jfl))/zufl(jfl)
(  244)                ENDIF
(  245)             ENDIF
(  246)             ! v-direction
(  247)             zvdfl (jfl) = zvoutfl - zvinfl
(  248)             zgjdfl(jfl) = float(ijoutfl(jfl)-ijinfl(jfl))
(  249)             IF( zvfl(jfl)*zvoutfl <= 0. ) THEN
(  250)                ztyfl(jfl) = 1.E99
(  251)             ELSE
(  252)                IF( ABS(zvdfl(jfl)) >= 1.E-5 ) THEN
(  253)                   ztyfl(jfl) = zgjdfl(jfl)/zvdfl(jfl) * LOG(zvoutfl/zvfl (jfl))
(  254)                ELSE
(  255)                   ztyfl(jfl) = (float(ijoutfl(jfl)) - zgjfl(jfl))/zvfl(jfl)
(  256)                ENDIF
(  257)                IF( (ABS(zgjfl(jfl)-float(ijinfl (jfl))) <= 1.E-7) .OR.   &
(  258)                    (ABS(zgjfl(jfl)-float(ijoutfl(jfl))) <=  1.E-7) ) THEN
(  259)                   ztyfl(jfl) = (zgjdfl(jfl)) / zvfl(jfl)
(  260)                ENDIF
(  261)             ENDIF
(  262)             ! w-direction        
(  263)             IF( nisobfl(jfl) == 1. ) THEN 
(  264)                zwdfl (jfl) = zwoutfl - zwinfl
(  265)                zgkdfl(jfl) = float(ikoutfl(jfl) - ikinfl(jfl))
(  266)                IF( zwfl(jfl)*zwoutfl <= 0. ) THEN
(  267)                   ztzfl(jfl) = 1.E99
(  268)                ELSE
(  269)                   IF( ABS(zwdfl(jfl)) >= 1.E-5 ) THEN
(  270)                      ztzfl(jfl) = zgkdfl(jfl)/zwdfl(jfl) * LOG(zwoutfl/zwfl (jfl))
(  271)                   ELSE
(  272)                      ztzfl(jfl) = (float(ikoutfl(jfl)) - zgkfl(jfl))/zwfl(jfl)
(  273)                   ENDIF
(  274)                   IF( (ABS(zgkfl(jfl)-float(ikinfl (jfl))) <=  1.E-7) .OR.   &
(  275)                       (ABS(zgkfl(jfl)-float(ikoutfl(jfl))) <= 1.E-7) ) THEN
(  276)                      ztzfl(jfl) = (zgkdfl(jfl)) / zwfl(jfl)
(  277)                   ENDIF
(  278)                ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 6

(  279)             ENDIF
(  280)             
(  281)             ! the time to go leave the mesh is the smallest time
(  282)                    
(  283)             IF( nisobfl(jfl) == 1. ) THEN 
(  284)                zttfl(jfl) = MIN(ztxfl(jfl),ztyfl(jfl),ztzfl(jfl))
(  285)             ELSE
(  286)                zttfl(jfl) = MIN(ztxfl(jfl),ztyfl(jfl))
(  287)             ENDIF
(  288)             ! new age of the FLOAT
(  289)             zagenewfl(jfl) = zagefl(jfl) + zttfl(jfl)*zvol
(  290)             ! test to know if the "age" of the float is not bigger than the 
(  291)             ! time step
(  292)             IF( zagenewfl(jfl) > rdt ) THEN
(  293)                zttfl(jfl) = (rdt-zagefl(jfl)) / zvol
(  294)                zagenewfl(jfl) = rdt
(  295)             ENDIF
(  296)             
(  297)             ! In the "minimal" direction we compute the index of new mesh
(  298)             ! on i-direction
(  299)             IF( ztxfl(jfl) <=  zttfl(jfl) ) THEN
(  300)                zgifl(jfl) = float(iioutfl(jfl))
(  301)                ind = iioutfl(jfl)
(  302)                IF( iioutfl(jfl) >= iiinfl(jfl) ) THEN
(  303)                   iioutfl(jfl) = iioutfl(jfl) + 1
(  304)                ELSE
(  305)                   iioutfl(jfl) = iioutfl(jfl) - 1
(  306)                ENDIF
(  307)                iiinfl(jfl) = ind
(  308)             ELSE
(  309)                IF( ABS(zudfl(jfl)) >= 1.E-5 ) THEN 
(  310)                   zgifl(jfl) = zgifl(jfl) + zgidfl(jfl)*zufl(jfl)    &
(  311)                      &       * ( EXP( zudfl(jfl)/zgidfl(jfl)*zttfl(jfl) ) - 1. ) /  zudfl(jfl)
(  312)                ELSE
(  313)                   zgifl(jfl) = zgifl(jfl) + zufl(jfl) * zttfl(jfl)
(  314)                ENDIF
(  315)             ENDIF
(  316)             ! on j-direction
(  317)             IF( ztyfl(jfl) <= zttfl(jfl) ) THEN
(  318)                zgjfl(jfl) = float(ijoutfl(jfl))
(  319)                ind = ijoutfl(jfl)
(  320)                IF( ijoutfl(jfl) >= ijinfl(jfl) ) THEN
(  321)                   ijoutfl(jfl) = ijoutfl(jfl) + 1
(  322)                ELSE
(  323)                   ijoutfl(jfl) = ijoutfl(jfl) - 1
(  324)                ENDIF
(  325)                ijinfl(jfl) = ind
(  326)             ELSE
(  327)                IF( ABS(zvdfl(jfl)) >= 1.E-5 ) THEN 
(  328)                   zgjfl(jfl) = zgjfl(jfl)+zgjdfl(jfl)*zvfl(jfl)   &
(  329)                      &       * ( EXP(zvdfl(jfl)/zgjdfl(jfl)*zttfl(jfl)) - 1. ) /  zvdfl(jfl)
(  330)                ELSE
(  331)                   zgjfl(jfl) = zgjfl(jfl)+zvfl(jfl)*zttfl(jfl)
(  332)                ENDIF
(  333)             ENDIF
(  334)             ! on k-direction
(  335)             IF( nisobfl(jfl) == 1. ) THEN
(  336)                IF( ztzfl(jfl) <= zttfl(jfl) ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 7

(  337)                   zgkfl(jfl) = float(ikoutfl(jfl))
(  338)                   ind = ikoutfl(jfl)
(  339)                   IF( ikoutfl(jfl) >= ikinfl(jfl) ) THEN
(  340)                      ikoutfl(jfl) = ikoutfl(jfl)+1
(  341)                   ELSE
(  342)                      ikoutfl(jfl) = ikoutfl(jfl)-1
(  343)                   ENDIF
(  344)                   ikinfl(jfl) = ind
(  345)                ELSE
(  346)                   IF( ABS(zwdfl(jfl)) >= 1.E-5 ) THEN 
(  347)                      zgkfl(jfl) = zgkfl(jfl)+zgkdfl(jfl)*zwfl(jfl)    &
(  348)                         &       * ( EXP(zwdfl(jfl)/zgkdfl(jfl)*zttfl(jfl)) - 1. ) /  zwdfl(jfl)
(  349)                   ELSE
(  350)                      zgkfl(jfl) = zgkfl(jfl)+zwfl(jfl)*zttfl(jfl)
(  351)                   ENDIF
(  352)                ENDIF
(  353)             ENDIF
(  354)             
(  355)             ! coordinate of the new point on the temperature grid
(  356)             
(  357)             iil(jfl) = MAX(iiinfl(jfl),iioutfl(jfl))
(  358)             ijl(jfl) = MAX(ijinfl(jfl),ijoutfl(jfl))
(  359)             IF( nisobfl(jfl) ==  1 ) ikl(jfl) = MAX(ikinfl(jfl),ikoutfl(jfl))
(  360)             ! reinitialisation of the age of FLOAT
(  361)             zagefl(jfl) = zagenewfl(jfl)
(  362) # if   defined key_mpp_mpi   ||   defined key_mpp_shmem
(  363)          ELSE
(  364)             ! we put zgifl, zgjfl, zgkfl, zagefl
(  365)             zgifl (jfl) = 0.
(  366)             zgjfl (jfl) = 0.
(  367)             zgkfl (jfl) = 0.
(  368)             zagefl(jfl) = 0.
(  369)             iil(jfl) = 0
(  370)             ijl(jfl) = 0
(  371)          ENDIF
(  372) # endif
(  373)       END DO
(  374)       
(  375)       ! synchronisation
(  376)       IF( lk_mpp )   CALL mpp_sum( zgifl , jpnfl )   ! sums over the global domain
(  377)       IF( lk_mpp )   CALL mpp_sum( zgjfl , jpnfl )
(  378)       IF( lk_mpp )   CALL mpp_sum( zgkfl , jpnfl )
(  379)       IF( lk_mpp )   CALL mpp_sum( zagefl, jpnfl )
(  380)       IF( lk_mpp )   CALL mpp_sum( iil   , jpnfl )
(  381)       IF( lk_mpp )   CALL mpp_sum( ijl   , jpnfl )
(  382)       
(  383)       ! in the case of open boundaries we need to test if the floats don't
(  384)       ! go out of the domain. If it goes out, the float is put at the 
(  385)       ! middle of the mesh in the domain but the trajectory isn't compute 
(  386)       ! more time.      
(  387) # if defined key_obc
(  388) !!DB - 2008.03.10 -- errors in variable names if free-surface ----> fix 
(  389)       DO jfl = 1, jpnfl
(  390)          IF( lp_obc_east ) THEN
(  391) !            IF( jped <=  zgjfl(jfl) .AND. zgjfl(jfl) <= jpef .AND. nieob-1 <=  zgifl(jfl) ) THEN
(  392)             IF( jpjed <=  zgjfl(jfl) .AND. zgjfl(jfl) <= jpjef .AND. jpieob-1 <=  zgifl(jfl) ) THEN
(  393)                zgifl (jfl) = INT(zgifl(jfl)) + 0.5
(  394)                zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 8

(  395)                zagefl(jfl) = rdt
(  396)             END IF
(  397)          END IF
(  398)          IF( lp_obc_west ) THEN
(  399) !            IF( jpwd <= zgjfl(jfl) .AND. zgjfl(jfl) <= jpwf .AND. niwob >=  zgifl(jfl) ) THEN
(  400)             IF( jpjwd <= zgjfl(jfl) .AND. zgjfl(jfl) <= jpjwf .AND. jpiwob >=  zgifl(jfl) ) THEN
(  401)                zgifl (jfl) = INT(zgifl(jfl)) + 0.5
(  402)                zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
(  403)                zagefl(jfl) = rdt
(  404)             END IF
(  405)          END IF
(  406)          IF( lp_obc_north ) THEN
(  407) !            IF( jpnd <=  zgifl(jfl) .AND. zgifl(jfl) <= jpnf .AND. njnob-1 >=  zgjfl(jfl) ) THEN
(  408)             IF( jpind <=  zgifl(jfl) .AND. zgifl(jfl) <= jpinf .AND. jpjnob-1 >=  zgjfl(jfl) ) THEN
(  409)                zgifl (jfl) = INT(zgifl(jfl)) + 0.5
(  410)                zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
(  411)                zagefl(jfl) = rdt
(  412)             END IF
(  413)          END IF
(  414)          IF( lp_obc_south ) THEN
(  415) !            IF( jpsd <=  zgifl(jfl) .AND. zgifl(jfl) <= jpsf .AND.  njsob >= zgjfl(jfl) ) THEN
(  416)             IF( jpisd <=  zgifl(jfl) .AND. zgifl(jfl) <= jpisf .AND.  jpjsob >= zgjfl(jfl) ) THEN
(  417)                zgifl (jfl) = INT(zgifl(jfl)) + 0.5
(  418)                zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
(  419)                zagefl(jfl) = rdt
(  420)             END IF
(  421)          END IF
(  422)       END DO
(  423) #endif
(  424) 
(  425)       ! Test to know if a  float hasn't integrated enought time
(  426)       IF( ln_argo ) THEN
(  427)          ifin = 1
(  428)          DO jfl = 1, jpnfl
(  429)             IF( zagefl(jfl) < rdt )   ifin = 0
(  430)             tpifl(jfl) = zgifl(jfl) + 0.5
(  431)             tpjfl(jfl) = zgjfl(jfl) + 0.5
(  432)          END DO
(  433)       ELSE
(  434)          ifin = 1
(  435)          DO jfl = 1, jpnfl
(  436)             IF( zagefl(jfl) < rdt )   ifin = 0
(  437)             tpifl(jfl) = zgifl(jfl) + 0.5
(  438)             tpjfl(jfl) = zgjfl(jfl) + 0.5
(  439)             IF( nisobfl(jfl) == 1 ) tpkfl(jfl) = -(zgkfl(jfl))
(  440)          END DO
(  441)       ENDIF
(  442)       IF( ifin == 0 ) THEN
(  443)          iloop = iloop + 1 
(  444)          GO TO 222
(  445)       ENDIF
(  446) 
(  447)    END SUBROUTINE flo_blk
(  448) 
(  449) !!DB
(  450)    SUBROUTINE flo_RDM( kt )
(  451)       !!---------------------------------------------------------------------
(  452)       !!                  ***  ROUTINE flo_RDM  ***






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 9

(  453)       !!           
(  454)       !! ** Purpose :   Compute the geographical position,latitude, longitude
(  455)       !!      and depth of each float at each time step.
(  456)       !! 
(  457)       !! ** Method  :   The position of a float is updated using a simple Euler 
(  458)       !!      algorithm. We need to know the velocity field, the old positions
(  459)       !!      of the floats and the grid defined on the domain.
(  460)       !!
(  461)       !!----------------------------------------------------------------------
(  462)       !! * arguments
(  463)       INTEGER, INTENT( in  ) ::   kt ! ocean time step
(  464) 
(  465)       !! * Local declarations
(  466)       INTEGER :: jfl              ! dummy loop arguments
(  467)       INTEGER :: ind, ifin, iloop
(  468) !!DB
(  469) !      INTEGER , DIMENSION ( jpnfl )  ::   &
(  470)       INTEGER , DIMENSION (:),ALLOCATABLE, SAVE  ::   &
(  471)          iil, ijl, ikl,             &     ! index of nearest mesh
(  472)          iiloc , ijloc 
(  473) !!DB
(  474) !      REAL(wp) , DIMENSION ( jpnfl )  ::    &
(  475)       REAL(wp) , DIMENSION (:),ALLOCATABLE, SAVE  ::    &
(  476)          zgifl, zgjfl, zgkfl,       &     ! position of floats, index on velocity mesh.
(  477)          zufl, zvfl, zwfl
(  478)       REAL(wp)   ::       &
(  479)          zuinfl,zvinfl,zwinfl,      &     ! transport across the input face
(  480)          zuoutfl,zvoutfl,zwoutfl
(  481) !!DB
(  482)       REAL(wp), DIMENSION(:), ALLOCATABLE, SAVE :: dxfl, dyfl, dzfl, dxifl, dyjfl, dzkfl
(  483)       REAL(wp) ::   &
(  484)            val, hor_diff, dx_rdm, dy_rdm, Ah_lo, Ah_hi, dAh_dx, dAh_dy, &
(  485)            Ahx_top, Ahx_btm, Ahy_left, Ahy_right, Ah_bl, Ah_br, Ah_tl, Ah_tr
(  486)       integer :: iseed, i1, i2, i3, j1, j2, k3
(  487)       LOGICAL :: RDM = .TRUE.
(  488) 
(  489)       !!---------------------------------------------------------------------
(  490)       !!---------------------------------------------------------------------
(  491) 
(  492)       IF( kt == nit000 ) THEN
(  493)          IF(lwp) WRITE(numout,*)
(  494)          IF(lwp) WRITE(numout,*) 'DBG: flo_RDM : RDM trajectories for floats '
(  495)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(  496) !!DB
(  497) !!NB: iseed = 5 means that the same random numbers will be generated every run
(  498) !!    For a different series must "randomize" the assignment of iseed
(  499)          iseed = 5
(  500)          call sgrnd(iseed)   !!initialize Mersenne Twister
(  501)          val = gasdev(iseed)
(  502)          if(lwp) write(numout,*)'DBG: kt, gasdev = ', kt, val
(  503) !!DB Allocate arrays
(  504)          ALLOCATE(iil(jpnfl));ALLOCATE(ijl(jpnfl));ALLOCATE(ikl(jpnfl))
(  505)          ALLOCATE(iiloc(jpnfl));ALLOCATE(ijloc(jpnfl))
(  506) 
(  507)          ALLOCATE(zgifl(jpnfl));ALLOCATE(zgjfl(jpnfl));ALLOCATE(zgkfl(jpnfl))
(  508)          ALLOCATE(zufl(jpnfl));ALLOCATE(zvfl(jpnfl));ALLOCATE(zwfl(jpnfl))
(  509) 
(  510)          ALLOCATE(dxfl(jpnfl)); ALLOCATE(dyfl(jpnfl)); ALLOCATE(dzfl(jpnfl))






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 10

(  511)          ALLOCATE(dxifl(jpnfl)); ALLOCATE(dyjfl(jpnfl)); ALLOCATE(dzkfl(jpnfl))
(  512) 
(  513)       ENDIF
(  514) 
(  515)       ! Initialisation of parameters
(  516)       
(  517)       do jfl = 1, jpnfl
(  518)          ! index on the velocity grid 
(  519)          ! We considere k coordinate negative, with this transformation 
(  520)          ! the computation in the 3 direction is the same. 
(  521)          zgifl(jfl) = tpifl(jfl) - 0.5
(  522)          zgjfl(jfl) = tpjfl(jfl) - 0.5
(  523)          zgkfl(jfl) = MIN(-1.,-(tpkfl(jfl)))
(  524)          ! index of T mesh
(  525)          iil(jfl) = 1 + INT(zgifl(jfl))
(  526)          ijl(jfl) = 1 + INT(zgjfl(jfl))
(  527)          ikl(jfl) =     INT(zgkfl(jfl))
(  528)       enddo
(  529)        
(  530)       do jfl = 1, jpnfl
(  531) # if   defined key_mpp_mpi   ||   defined key_mpp_shmem
(  532)          IF( (iil(jfl) >= (mig(nldi)-jpizoom+1)) .AND. (iil(jfl) <= (mig(nlei)-jpizoom+1)) .AND.   &
(  533)              (ijl(jfl) >= (mjg(nldj)-jpjzoom+1)) .AND. (ijl(jfl) <= (mjg(nlej)-jpjzoom+1)) ) THEN
(  534)             iiloc(jfl) = iil(jfl) - (mig(1)-jpizoom+1) + 1
(  535)             ijloc(jfl) = ijl(jfl) - (mjg(1)-jpjzoom+1) + 1
(  536) # else 
(  537)             iiloc(jfl) = iil(jfl)
(  538)             ijloc(jfl) = ijl(jfl)
(  539) # endif
(  540)             
(  541)             ! compute the transport across the mesh where the float is.            
(  542) !            zsurfx(1) = e2u(iiloc(jfl)-1,ijloc(jfl)  ) * e3t(-ikl(jfl))
(  543) !            zsurfx(2) = e2u(iiloc(jfl)  ,ijloc(jfl)  ) * e3t(-ikl(jfl))
(  544) !            zsurfy(1) = e1v(iiloc(jfl)  ,ijloc(jfl)-1) * e3t(-ikl(jfl))
(  545) !            zsurfy(2) = e1v(iiloc(jfl)  ,ijloc(jfl)  ) * e3t(-ikl(jfl))
(  546) 
(  547)             ! for a isobar float zsurfz is put to zero. The vertical velocity will be zero too.
(  548) !            zsurfz=  e1t(iiloc(jfl),ijloc(jfl)) * e2t(iiloc(jfl),ijloc(jfl))
(  549) !            zvol  =( e1t(iiloc(jfl),ijloc(jfl)) * e2t(iiloc(jfl),ijloc(jfl)) * e3t(-ikl(jfl)) )
(  550) 
(  551) 
(  552) !!DB  Modify to compute velocity only
(  553)             !
(  554)             zuinfl =( ub(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl)) )/2.
(  555)             zuoutfl=( ub(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) + un(iiloc(jfl)  ,ijloc(jfl),-ikl(jfl)) )/2.
(  556)             zvinfl =( vb(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl)) )/2.
(  557)             zvoutfl=( vb(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) + vn(iiloc(jfl),ijloc(jfl)  ,-ikl(jfl)) )/2.
(  558)             zwinfl =-(wb(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1))    &
(  559)                &   +  wn(iiloc(jfl),ijloc(jfl),-(ikl(jfl)-1)) )/2. * nisobfl(jfl)
(  560)             zwoutfl=-(wb(iiloc(jfl),ijloc(jfl),- ikl(jfl)   )   &
(  561)                &   +  wn(iiloc(jfl),ijloc(jfl),- ikl(jfl)   ) )/2. * nisobfl(jfl)
(  562)             
(  563)             ! interpolation of velocity field on the float initial position            
(  564)             zufl(jfl)=  zuinfl  + ( zgifl(jfl) - float(iil(jfl)-1) ) * ( zuoutfl - zuinfl)
(  565)             zvfl(jfl)=  zvinfl  + ( zgjfl(jfl) - float(ijl(jfl)-1) ) * ( zvoutfl - zvinfl)
(  566)             zwfl(jfl)=  zwinfl  + ( zgkfl(jfl) - float(ikl(jfl)-1) ) * ( zwoutfl - zwinfl)
(  567) 
(  568) !!DB update position






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 11

(  569) !!the zgifl index holds the global U index == fractional position on the vel-grid
(  570) !!tpifl holds the global T index == fractional position on the T-grid
(  571) !!These are displacements in meters (NB: could be -ve)
(  572)             dxfl(jfl) = zufl(jfl)*rdt
(  573)             dyfl(jfl) = zvfl(jfl)*rdt
(  574) ! ignore z-direction for now
(  575) !            dzfl(jfl) = zwfl(jfl)*rdt
(  576) 
(  577) 
(  578) !!DB: correct Random Displacement Model
(  579) !! ie it uses the local Ah plus its spatial derivative
(  580) !!Algorithm: 
(  581) !!Part 1: Ah interpolated to float location; use F-pt diffusivity and bilinear interp
(  582) !! Displacement ~ N(0,Ah_local)
(  583) !!Part 2: Grad of Ah
(  584) !! Add correction to displacement to account for horizontal gradients in diffusivity
(  585) !!NB: z-direction ignored
(  586)             if(RDM) then
(  587) 
(  588) !!Code for tracer diffusivity
(  589) #ifdef  key_traldf_smag
(  590)                i1 = max(iiloc(jfl)-1,1)
(  591)                i2 = iiloc(jfl)
(  592)                i3 = min(iiloc(jfl)+1,jpi)
(  593)                j1 = max(ijloc(jfl)-1,1)
(  594)                j2 = ijloc(jfl)
(  595)                k3 = -ikl(jfl)
(  596)                Ah_tl = 0.5*(ahtv(i1,j2,k3)+ahtv(i2,j2,k3))
(  597)                Ah_tr = 0.5*(ahtv(i2,j2,k3)+ahtv(i3,j2,k3))
(  598)                Ah_br = 0.5*(ahtv(i2,j1,k3)+ahtv(i3,j1,k3))
(  599)                Ah_bl = 0.5*(ahtv(i1,j1,k3)+ahtv(i2,j1,k3))
(  600) #else
(  601) !Use momentum diffusivity 
(  602) !-0- Ah at 4 corners of cell: tr=top-right; bl=btm-left, etc
(  603)                Ah_tl = ahm2(iiloc(jfl)-1,ijloc(jfl),-ikl(jfl))
(  604)                Ah_tr = ahm2(iiloc(jfl),ijloc(jfl),-ikl(jfl))
(  605)                Ah_br = ahm2(iiloc(jfl),ijloc(jfl)-1,-ikl(jfl))
(  606)                Ah_bl = ahm2(iiloc(jfl)-1,ijloc(jfl)-1,-ikl(jfl))
(  607) #endif
(  608) 
(  609) !!Part 1:
(  610) !-1- Ah at x-position at top of cell
(  611)                Ahx_top = Ah_tl + (zgifl(jfl) - float(iil(jfl)-1))*(Ah_tr-Ah_tl) 
(  612) !-2- Ah at x-position at btm of cell
(  613)                Ahx_btm = Ah_bl + (zgifl(jfl) - float(iil(jfl)-1))*(Ah_br-Ah_bl) 
(  614) !-3- Interp to y-position
(  615)                hor_diff = Ahx_btm + (zgjfl(jfl) - float(ijl(jfl)-1))*(Ahx_top-Ahx_btm)
(  616)                dx_rdm = sqrt(2.0*hor_diff*rdt)*gasdev(iseed)
(  617)                dy_rdm = sqrt(2.0*hor_diff*rdt)*gasdev(iseed)
(  618)                dxfl(jfl) = dxfl(jfl) + dx_rdm
(  619)                dyfl(jfl) = dyfl(jfl) + dy_rdm
(  620) 
(  621) !!Part 2: Compute drift correction term ~ d(Ah)/dx * rdt
(  622) !!DB: need local gradients of Ah, which is more complicated as local Ah indices +/- 2 are needed
(  623) !!which can overindex local domain. A solution is to compute grad(Ah) everywhere, call lbc_lnk,
(  624) !!and then interpolate the halo-updated field. This is too computationally expensive, so I will
(  625) !!use a simpler method. 
(  626) !!Note that the bilinear interpolation used above is equivalent to assuming constant gradients






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 12

(  627) !!within the cell ====>
(  628)                dAh_dy = (Ahx_top-Ahx_btm)/e2v(iiloc(jfl),ijloc(jfl)-1)
(  629)                Ahy_left = Ah_bl + (zgjfl(jfl) - float(ijl(jfl)-1))*(Ah_tl-Ah_bl)
(  630)                Ahy_right = Ah_br + (zgjfl(jfl) - float(ijl(jfl)-1))*(Ah_tr-Ah_br)
(  631)                dAh_dx = (Ahy_right-Ahy_left)/e1u(iiloc(jfl)-1,ijloc(jfl))
(  632) 
(  633)                hor_diff = dAh_dx * rdt
(  634)                dxfl(jfl) = dxfl(jfl) + hor_diff
(  635)                hor_diff = dAh_dy * rdt
(  636)                dyfl(jfl) = dyfl(jfl) + hor_diff
(  637) 
(  638)             endif   !!RDM true
(  639) 
(  640) !displacement relative to cell size
(  641)             dxifl(jfl) = dxfl(jfl)/e1u(iiloc(jfl)-1,ijloc(jfl))
(  642)             dyjfl(jfl) = dyfl(jfl)/e2v(iiloc(jfl),ijloc(jfl)-1)
(  643) 
(  644) 
(  645) 
(  646) !!DB: 04.04: I'll ignore 2 possibilities (for now):
(  647) !-1- If particle leaves cell then relative displacement should be adjusted to account
(  648) !for part of its travel in the new cell where e1u may be different. As this is likely a
(  649) !small difference, I'll ignore it for now and fix it at some future time (likely never!)
(  650) !-2- displacement is > cell size; ie dx = U*rdt > e1u: This is like a CFL condition
(  651) !which I'll assume is satisfied as dt = rdt (ie routine called every kt). If the routine
(  652) !changes so that it is called every N*rdt seconds then we would have to check if displacement
(  653) !is greater than 1 cell size and adjust as necessary.
(  654) 
(  655)             zgifl(jfl) = zgifl(jfl) + dxifl(jfl)
(  656)             zgjfl(jfl) = zgjfl(jfl) + dyjfl(jfl)
(  657) ! ignore w-direction for now
(  658) 
(  659) 
(  660) # if   defined key_mpp_mpi   ||   defined key_mpp_shmem
(  661)          ELSE
(  662)             ! we put zgifl, zgjfl, zgkfl
(  663)             zgifl (jfl) = 0.
(  664)             zgjfl (jfl) = 0.
(  665)             zgkfl (jfl) = 0.
(  666)             iil(jfl) = 0
(  667)             ijl(jfl) = 0
(  668)          ENDIF
(  669) # endif
(  670)       END DO
(  671)       
(  672)       ! synchronisation
(  673)       IF( lk_mpp )   CALL mpp_sum( zgifl , jpnfl )   ! sums over the global domain
(  674)       IF( lk_mpp )   CALL mpp_sum( zgjfl , jpnfl )
(  675)       IF( lk_mpp )   CALL mpp_sum( zgkfl , jpnfl )
(  676)       IF( lk_mpp )   CALL mpp_sum( iil   , jpnfl )
(  677)       IF( lk_mpp )   CALL mpp_sum( ijl   , jpnfl )
(  678)       
(  679) !!DB: Make sure float is inside domain. 
(  680) !!At this point all processors hold the same zg*fl variables, in global indices 
(  681) !!Cavalier at this stage but DBG 
(  682)       do jfl = 1, jpnfl      
(  683)          if(zgifl(jfl) > jpidta-1) zgifl(jfl) = jpidta-1
(  684)          if(zgjfl(jfl) > jpjdta-1) zgjfl(jfl) = jpjdta-1






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 13

(  685)          if(zgifl(jfl) < 1) zgifl(jfl) = 1
(  686)          if(zgjfl(jfl) < 1) zgjfl(jfl) = 1
(  687)       enddo
(  688) 
(  689) 
(  690) !!DB: the below (comment and code) is old and I do not think it is necessary 
(  691) !!I think what I did above is OK
(  692)       ! in the case of open boundaries we need to test if the floats don't
(  693)       ! go out of the domain. If it goes out, the float is put at the 
(  694)       ! middle of the mesh in the domain 
(  695) # if defined key_obc
(  696)       DO jfl = 1, jpnfl
(  697)          IF( lp_obc_east ) THEN
(  698)             IF( jpjed <=  zgjfl(jfl) .AND. zgjfl(jfl) <= jpjef .AND. jpieob-1 <=  zgifl(jfl) ) THEN
(  699)                zgifl (jfl) = INT(zgifl(jfl)) + 0.5
(  700)                zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
(  701)             END IF
(  702)          END IF
(  703)          IF( lp_obc_west ) THEN
(  704)             IF( jpjwd <= zgjfl(jfl) .AND. zgjfl(jfl) <= jpjwf .AND. jpiwob >=  zgifl(jfl) ) THEN
(  705)                zgifl (jfl) = INT(zgifl(jfl)) + 0.5
(  706)                zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
(  707)             END IF
(  708)          END IF
(  709)          IF( lp_obc_north ) THEN
(  710)             IF( jpind <=  zgifl(jfl) .AND. zgifl(jfl) <= jpinf .AND. jpjnob-1 >=  zgjfl(jfl) ) THEN
(  711)                zgifl (jfl) = INT(zgifl(jfl)) + 0.5
(  712)                zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
(  713)             END IF
(  714)          END IF
(  715)          IF( lp_obc_south ) THEN
(  716)             IF( jpisd <=  zgifl(jfl) .AND. zgifl(jfl) <= jpisf .AND.  jpjsob >= zgjfl(jfl) ) THEN
(  717)                zgifl (jfl) = INT(zgifl(jfl)) + 0.5
(  718)                zgjfl (jfl) = INT(zgjfl(jfl)) + 0.5
(  719)             END IF
(  720)          END IF
(  721)       END DO
(  722) #endif
(  723) 
(  724) !!DB
(  725)       do jfl = 1, jpnfl
(  726)          tpifl(jfl) = zgifl(jfl) + 0.5
(  727)          tpjfl(jfl) = zgjfl(jfl) + 0.5
(  728) !!Not sure what to do here
(  729) !         IF( nisobfl(jfl) == 1 ) tpkfl(jfl) = -(zgkfl(jfl))
(  730)       enddo
(  731) 
(  732) 
(  733) 
(  734)    END SUBROUTINE flo_RDM
(  735) 
(  736) 
(  737) !-----------------------------------------------------------------------
(  738) !
(  739)       FUNCTION GASDEV(IDUM)
(  740) !
(  741) ! GASDEV - GASDEV RETURNS A NORMALLY DISTRIBUTED DEVIATE
(  742) !          WITH ZERO MEAN AND UNIT VARIAN!E, USING RAN1(IDUM) AS






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 14

(  743) !          THE SOURCE FOR THE UNIFORM DEVIATES.
(  744) !
(  745)       INTEGER :: IDUM
(  746)       REAL :: GASDEV
(  747)       REAL(8) :: n1, n2
(  748)       INTEGER :: ISET
(  749)       REAL :: FAC,GSET,RSQ,V1,V2,RAN1,RAN2
(  750)       !CN: no longer needed since grnd is a subroutine
(  751)       !double precision :: grnd
(  752)       SAVE ISET,GSET
(  753)       DATA ISET/0/
(  754)       IF (ISET.EQ.0) THEN
(  755)         !CN: Changed grnd to subroutine to fix linker errors
(  756)         CALL grnd(n1)
(  757)         CALL grnd(n2)
(  758) 1       V1=2.*n1-1.
(  759)         V2=2.*n2-1.
(  760)         !V1=2.*grnd()-1.
(  761)         !V2=2.*grnd()-1. 
(  762)         RSQ=V1**2+V2**2
(  763)         IF(RSQ.GE.1..OR.RSQ.EQ.0.)GOTO 1
(  764)         FAC=SQRT(-2.*LOG(RSQ)/RSQ)
(  765)         GSET=V1*FAC
(  766)         GASDEV=V2*FAC
(  767)         ISET=1
(  768)       ELSE
(  769)         GASDEV=GSET
(  770)         ISET=0
(  771)       ENDIF
(  772)       RETURN
(  773) 
(  774)       end function GASDEV
(  775) 
(  776) 
(  777) 
(  778) !  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE &OL`.
(  779) !
(  780) 
(  781) !     MERSENNE TWISTER CODE
(  782) ! A C-program for MT19937: Real number version
(  783) !   genrand() generates one pseudorandom real number (double)
(  784) ! which is uniformly distributed on [0,1]-interval, for each
(  785) ! call. sgenrand(seed) set initial values to the working area
(  786) ! of 624 words. Before genrand(), sgenrand(seed) must be
(  787) ! called once. (seed is any 32-bit integer except for 0).
(  788) ! Integer generator is obtained by modifying two lines.
(  789) !   Coded by Takuji Nishimura, considering the suggestions by
(  790) ! Topher Cooper and Marc Rieffel in July-Aug. 1997.
(  791) !
(  792) ! This library is free software; you can redistribute it and/or
(  793) ! modify it under the terms of the GNU Library General Public
(  794) ! License as published by the Free Software Foundation; either
(  795) ! version 2 of the License, or (at your option) any later
(  796) ! version.
(  797) ! This library is distributed in the hope that it will be useful,
(  798) ! but WITHOUT ANY WARRANTY; without even the implied warranty of
(  799) ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
(  800) ! See the GNU Library General Public License for more details.






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 15

(  801) ! You should have received a copy of the GNU Library General
(  802) ! Public License along with this library; if not, write to the
(  803) ! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
(  804) ! 02111-1307  USA
(  805) !
(  806) ! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
(  807) ! When you use this, send an email to: matumoto@math.keio.ac.jp
(  808) ! with an appropriate reference to your work.
(  809) !
(  810) !***********************************************************************
(  811) ! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
(  812) !
(  813) !   genrand()      -> double precision function grnd()
(  814) !   sgenrand(seed) -> subroutine sgrnd(seed)
(  815) !                     integer seed
(  816) !
(  817) ! This program uses the following non-standard intrinsics.
(  818) !   ishft(i,n): If n>0, shifts bits in i by n positions to left.
(  819) !               If n<0, shifts bits in i by n positions to right.
(  820) !   iand (i,j): Performs logical AND on corresponding bits of i and j.
(  821) !   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
(  822) !   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
(  823) !
(  824) !***********************************************************************
(  825) ! Fortran version rewritten as an F90 module and mt state saving and getting
(  826) ! subroutines added by Richard Woloshyn. (rwww@triumf.ca). June 30, 1999
(  827) 
(  828) !Initialization subroutine
(  829)   subroutine sgrnd(seed)
(  830)     implicit none
(  831) !
(  832) !      setting initial seeds to mt[N] using
(  833) !      the generator Line 25 of Table 1 in
(  834) !      [KNUTH 1981, The Art of Computer Programming
(  835) !         Vol. 2 (2nd Ed.), pp102]
(  836) !
(  837)     integer, intent(in) :: seed
(  838) 
(  839)     mt(0) = iand(seed,-1)
(  840)     do mti=1,N-1
(  841)       mt(mti) = iand(69069 * mt(mti-1),-1)
(  842)     enddo
(  843) !
(  844)     return
(  845)   end subroutine sgrnd
(  846) 
(  847) !Random number generator
(  848) !CN: Changed grnd to a subroutine to fix linker errors
(  849)   subroutine grnd(rnd_num)
(  850)     implicit integer(a-z)
(  851) 
(  852)     REAL(8),INTENT(OUT) :: rnd_num
(  853) 
(  854) ! Period parameters
(  855)     integer, parameter :: M = 397, MATA  = -1727483681
(  856) !                                    constant vector a
(  857)     integer, parameter :: LMASK =  2147483647
(  858) !                                    least significant r bits






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 16

(  859)     integer, parameter :: UMASK = -LMASK - 1
(  860) !                                    most significant w-r bits
(  861) ! Tempering parameters
(  862)     integer, parameter :: TMASKB= -1658038656, TMASKC= -272236544
(  863) 
(  864)     dimension mag01(0:1)
(  865)     data mag01/0, MATA/
(  866)     save mag01
(  867) !                        mag01(x) = x * MATA for x=0,1
(  868) 
(  869)     TSHFTU(y)=ishft(y,-11)
(  870)     TSHFTS(y)=ishft(y,7)
(  871)     TSHFTT(y)=ishft(y,15)
(  872)     TSHFTL(y)=ishft(y,-18)
(  873) 
(  874)     if(mti.ge.N) then
(  875) !                       generate N words at one time
(  876)       if(mti.eq.N+1) then
(  877) !                            if sgrnd() has not been called,
(  878)         call sgrnd( defaultsd )
(  879) !                              a default initial seed is used
(  880)       endif
(  881) 
(  882)       do kk=0,N-M-1
(  883)           y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
(  884)           mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
(  885)       enddo
(  886)       do kk=N-M,N-2
(  887)           y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
(  888)           mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
(  889)       enddo
(  890)       y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
(  891)       mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
(  892)       mti = 0
(  893)     endif
(  894) 
(  895)     y=mt(mti)
(  896)     mti = mti + 1 
(  897)     y=ieor(y,TSHFTU(y))
(  898)     y=ieor(y,iand(TSHFTS(y),TMASKB))
(  899)     y=ieor(y,iand(TSHFTT(y),TMASKC))
(  900)     y=ieor(y,TSHFTL(y))
(  901) 
(  902)     if(y .lt. 0) then
(  903)       rnd_num=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
(  904)     else
(  905)       rnd_num=dble(y)/(2.0d0**32-1.0d0)
(  906)     endif
(  907) 
(  908)     return
(  909)   end subroutine grnd
(  910) 
(  911) #  else
(  912)    !!----------------------------------------------------------------------
(  913)    !!   Default option                                         Empty module
(  914)    !!----------------------------------------------------------------------
(  915) CONTAINS
(  916)    SUBROUTINE flo_blk                  ! Empty routine






PGF90 (Version     10.5)          09/23/2010  09:49:47      page 17

(  917)    END SUBROUTINE flo_blk 
(  918) #endif
(  919)    
(  920)    !!======================================================================
(  921) END MODULE floblk 



























































PGF90 (Version     10.5)          09/23/2010  09:49:47      page 18

(    1) # 1 "floblk.F90"
(    1) !!DB 2008.03.19 ...
(    2) !!Added a Random Displacement Model (RDM) = particle tracking plus turbulence
(    3) !!Routine is called flo_RDM. It replaces the flo_4rk() routine if(ln_flork4)
(    4) !!NB: the Rk4 routine did not work for me.
(    5) !!To use RDM scheme: in namelist set ln_flork4 = .true. 
(    6) !!Requires a init_float file with float positions which I have also modified. 
(    7) !!See flodom.F90 and look for DB
(    8) !!NB: Routine works but it is not really complete. For example:
(    9) !!   (a) it ignores w-direction, i.e. floats remain at initial level
(   10) !!   (b) following from (a): there is no vertical turbulence
(   11) !!   (c) Add more notes ... TO DO ...
(   12) 
(   13) 
(   14) MODULE floblk
(   15)    !!======================================================================
(   16)    !!                     ***  MODULE  floblk  ***
(   17)    !! Ocean floats :   trajectory computation
(   18)    !!======================================================================
(   19) # 912
(  912)    !!----------------------------------------------------------------------
(  913)    !!   Default option                                         Empty module
(  914)    !!----------------------------------------------------------------------
(  915) CONTAINS
(  916)    SUBROUTINE flo_blk                  ! Empty routine
(  917)    END SUBROUTINE flo_blk 
(  918) # 919
(  919)    
(  920)    !!======================================================================
(  921) END MODULE floblk 
