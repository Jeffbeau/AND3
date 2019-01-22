


PGF90 (Version     10.5)          09/23/2010  09:49:02      page 1

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

Filename: lbclnk.F90

(    1) MODULE lbclnk
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  lbclnk  ***
(    4)    !! Ocean        : lateral boundary conditions
(    5)    !!=====================================================================
(    6) #if   defined key_mpp_mpi   ||   defined key_mpp_shmem
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_mpp_mpi'     OR      MPI massively parallel processing library
(    9)    !!   'key_mpp_shmem'         SHMEM massively parallel processing library
(   10)    !!----------------------------------------------------------------------
(   11)    !!----------------------------------------------------------------------
(   12)    !!   lbc_lnk      : generic interface for mpp_lnk_3d and mpp_lnk_2d
(   13)    !!                  routines defined in lib_mpp
(   14)    !!   lbc_lnk_e    : generic interface for mpp_lnk_2d_e
(   15)    !!                   routinee defined in lib_mpp
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE lib_mpp          ! distributed memory computing library
(   19) 
(   20)    INTERFACE lbc_lnk
(   21)       MODULE PROCEDURE mpp_lnk_3d, mpp_lnk_2d
(   22)    END INTERFACE
(   23) 
(   24)    INTERFACE lbc_lnk_e
(   25)       MODULE PROCEDURE mpp_lnk_2d_e
(   26)    END INTERFACE
(   27) 
(   28)    PUBLIC lbc_lnk       ! ocean lateral boundary conditions
(   29)    PUBLIC lbc_lnk_e
(   30)    !!----------------------------------------------------------------------
(   31) 
(   32) #else
(   33)    !!----------------------------------------------------------------------
(   34)    !!   Default option                              shared memory computing
(   35)    !!----------------------------------------------------------------------
(   36)    !!   lbc_lnk      : generic interface for lbc_lnk_3d and lbc_lnk_2d
(   37)    !!   lbc_lnk_3d   : set the lateral boundary condition on a 3D variable
(   38)    !!                  on OPA ocean mesh
(   39)    !!   lbc_lnk_2d   : set the lateral boundary condition on a 2D variable
(   40)    !!                  on OPA ocean mesh
(   41)    !!----------------------------------------------------------------------
(   42)    !! * Modules used
(   43)    USE oce             ! ocean dynamics and tracers   
(   44)    USE dom_oce         ! ocean space and time domain 
(   45)    USE in_out_manager  ! I/O manager
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 2

(   47)    IMPLICIT NONE
(   48)    PRIVATE
(   49) 
(   50)    INTERFACE lbc_lnk
(   51)       MODULE PROCEDURE lbc_lnk_3d, lbc_lnk_2d
(   52)    END INTERFACE
(   53) 
(   54)    INTERFACE lbc_lnk_e
(   55)       MODULE PROCEDURE lbc_lnk_2d
(   56)    END INTERFACE
(   57) 
(   58)    PUBLIC lbc_lnk       ! ocean/ice  lateral boundary conditions
(   59)    PUBLIC  lbc_lnk_e 
(   60)    !!----------------------------------------------------------------------
(   61) 
(   62) CONTAINS
(   63) 
(   64)    SUBROUTINE lbc_lnk_3d( pt3d, cd_type, psgn )
(   65)       !!---------------------------------------------------------------------
(   66)       !!                  ***  ROUTINE lbc_lnk_3d  ***
(   67)       !!
(   68)       !! ** Purpose :   set lateral boundary conditions (non mpp case)
(   69)       !!
(   70)       !! ** Method  :
(   71)       !!
(   72)       !! History :
(   73)       !!        !  97-06  (G. Madec)  Original code
(   74)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(   75)       !!----------------------------------------------------------------------
(   76)       !! * Arguments
(   77)       CHARACTER(len=1), INTENT( in ) ::   &
(   78)          cd_type       ! nature of pt3d grid-points
(   79)          !             !   = T ,  U , V , F or W  gridpoints
(   80)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ) ::   &
(   81)          pt3d          ! 3D array on which the boundary condition is applied
(   82)       REAL(wp), INTENT( in ) ::   &
(   83)          psgn          ! control of the sign change
(   84)          !             !   =-1 , the sign is changed if north fold boundary
(   85)          !             !   = 1 , no sign change
(   86)          !             !   = 0 , no sign change and > 0 required (use the inner
(   87)          !             !         row/column if closed boundary)
(   88) 
(   89)       !! * Local declarations
(   90)       INTEGER  ::   ji, jk
(   91)       INTEGER  ::   ijt, iju
(   92)       !!----------------------------------------------------------------------
(   93)       !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   94)       !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lbclnk.F90,v 1.7 2005/09/30 10:18:51 opalod Exp $ 
(   95)       !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   96)       !!----------------------------------------------------------------------
(   97)       
(   98)       !                                                      ! ===============
(   99)       DO jk = 1, jpk                                         ! Horizontal slab
(  100)          !                                                   ! ===============
(  101) 
(  102)          !                                     ! East-West boundaries
(  103)          !                                     ! ====================
(  104)          SELECT CASE ( nperio )






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 3

(  105) 
(  106)          CASE ( 1 , 4 , 6 )                    ! * cyclic east-west
(  107)             pt3d( 1 ,:,jk) = pt3d(jpim1,:,jk)          ! all points
(  108)             pt3d(jpi,:,jk) = pt3d(  2  ,:,jk)
(  109) 
(  110)          CASE DEFAULT                          ! * closed
(  111)             SELECT CASE ( cd_type )
(  112)             CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
(  113)                pt3d( 1 ,:,jk) = 0.e0
(  114)                pt3d(jpi,:,jk) = 0.e0
(  115)             CASE ( 'F' )                               ! F-point
(  116)                pt3d(jpi,:,jk) = 0.e0
(  117)             END SELECT
(  118) 
(  119)          END SELECT
(  120) 
(  121)          !                                     ! North-South boundaries
(  122)          !                                     ! ======================
(  123)          SELECT CASE ( nperio )
(  124) 
(  125)          CASE ( 2 )                            ! *  south symmetric
(  126) 
(  127)             SELECT CASE ( cd_type )
(  128)             CASE ( 'T' , 'U' , 'W' )                   ! T-, U-, W-points
(  129)                pt3d(:, 1 ,jk) = pt3d(:,3,jk)
(  130)                pt3d(:,jpj,jk) = 0.e0
(  131)             CASE ( 'V' , 'F' )                         ! V-, F-points
(  132)                pt3d(:, 1 ,jk) = psgn * pt3d(:,2,jk)
(  133)                pt3d(:,jpj,jk) = 0.e0
(  134)             END SELECT
(  135) 
(  136)          CASE ( 3 , 4 )                        ! *  North fold  T-point pivot
(  137) 
(  138)             pt3d( 1 ,jpj,jk) = 0.e0
(  139)             pt3d(jpi,jpj,jk) = 0.e0
(  140) 
(  141)             SELECT CASE ( cd_type )
(  142)             CASE ( 'T' , 'W' )                         ! T-, W-point
(  143)                DO ji = 2, jpi
(  144)                   ijt = jpi-ji+2
(  145)                   pt3d(ji, 1 ,jk) = 0.e0
(  146)                   pt3d(ji,jpj,jk) = psgn * pt3d(ijt,jpj-2,jk)
(  147)                END DO
(  148)                DO ji = jpi/2+1, jpi
(  149)                   ijt = jpi-ji+2
(  150)                   pt3d(ji,jpjm1,jk) = psgn * pt3d(ijt,jpjm1,jk)
(  151)                END DO
(  152)             CASE ( 'U' )                               ! U-point
(  153)                DO ji = 1, jpi-1
(  154)                   iju = jpi-ji+1
(  155)                   pt3d(ji, 1 ,jk) = 0.e0
(  156)                   pt3d(ji,jpj,jk) = psgn * pt3d(iju,jpj-2,jk)
(  157)                END DO
(  158)                DO ji = jpi/2, jpi-1
(  159)                   iju = jpi-ji+1
(  160)                   pt3d(ji,jpjm1,jk) = psgn * pt3d(iju,jpjm1,jk)
(  161)                END DO
(  162)             CASE ( 'V' )                               ! V-point






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 4

(  163)                   DO ji = 2, jpi
(  164)                      ijt = jpi-ji+2
(  165)                      pt3d(ji,  1  ,jk) = 0.e0
(  166)                      pt3d(ji,jpj-1,jk) = psgn * pt3d(ijt,jpj-2,jk)
(  167)                      pt3d(ji,jpj  ,jk) = psgn * pt3d(ijt,jpj-3,jk)
(  168)                   END DO
(  169)             CASE ( 'F' )                               ! F-point
(  170)                   DO ji = 1, jpi-1
(  171)                      iju = jpi-ji+1
(  172)                      pt3d(ji,jpj-1,jk) = psgn * pt3d(iju,jpj-2,jk)
(  173)                      pt3d(ji,jpj  ,jk) = psgn * pt3d(iju,jpj-3,jk)
(  174)                   END DO
(  175)             END SELECT
(  176) 
(  177)          CASE ( 5 , 6 )                        ! *  North fold  F-point pivot
(  178) 
(  179)             pt3d( 1 ,jpj,jk) = 0.e0
(  180)             pt3d(jpi,jpj,jk) = 0.e0
(  181) 
(  182)             SELECT CASE ( cd_type )
(  183)             CASE ( 'T' , 'W' )                         ! T-, W-point
(  184)                DO ji = 1, jpi
(  185)                   ijt = jpi-ji+1
(  186)                   pt3d(ji, 1 ,jk) = 0.e0
(  187)                   pt3d(ji,jpj,jk) = psgn * pt3d(ijt,jpj-1,jk)
(  188)                END DO
(  189)             CASE ( 'U' )                               ! U-point
(  190)                   DO ji = 1, jpi-1
(  191)                      iju = jpi-ji
(  192)                      pt3d(ji, 1 ,jk) = 0.e0
(  193)                      pt3d(ji,jpj,jk) = psgn * pt3d(iju,jpj-1,jk)
(  194)                   END DO
(  195)             CASE ( 'V' )                               ! V-point
(  196)                   DO ji = 1, jpi
(  197)                      ijt = jpi-ji+1
(  198)                      pt3d(ji, 1 ,jk) = 0.e0
(  199)                      pt3d(ji,jpj,jk) = psgn * pt3d(ijt,jpj-2,jk)
(  200)                   END DO
(  201)                   DO ji = jpi/2+1, jpi
(  202)                      ijt = jpi-ji+1
(  203)                      pt3d(ji,jpjm1,jk) = psgn * pt3d(ijt,jpjm1,jk)
(  204)                   END DO
(  205)             CASE ( 'F' )                               ! F-point
(  206)                   DO ji = 1, jpi-1
(  207)                      iju = jpi-ji
(  208)                      pt3d(ji,jpj  ,jk) = psgn * pt3d(iju,jpj-2,jk)
(  209)                   END DO
(  210)                   DO ji = jpi/2+1, jpi-1
(  211)                      iju = jpi-ji
(  212)                      pt3d(ji,jpjm1,jk) = psgn * pt3d(iju,jpjm1,jk)
(  213)                   END DO
(  214)             END SELECT
(  215) 
(  216)          CASE DEFAULT                          ! *  closed
(  217) 
(  218)             SELECT CASE ( cd_type )
(  219)             CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
(  220)                pt3d(:, 1 ,jk) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 5

(  221)                pt3d(:,jpj,jk) = 0.e0
(  222)             CASE ( 'F' )                               ! F-point
(  223)                pt3d(:,jpj,jk) = 0.e0
(  224)             END SELECT
(  225) 
(  226)          END SELECT
(  227)          !                                                   ! ===============
(  228)       END DO                                                 !   End of slab
(  229)       !                                                      ! ===============
(  230)    END SUBROUTINE lbc_lnk_3d
(  231) 
(  232) 
(  233)    SUBROUTINE lbc_lnk_2d( pt2d, cd_type, psgn )
(  234)       !!---------------------------------------------------------------------
(  235)       !!                 ***  ROUTINE lbc_lnk_2d  ***
(  236)       !!
(  237)       !! ** Purpose :   set lateral boundary conditions (non mpp case)
(  238)       !!
(  239)       !! ** Method  :
(  240)       !!
(  241)       !! History :
(  242)       !!        !  97-06  (G. Madec)  Original code
(  243)       !!        !  01-05  (E. Durand)  correction
(  244)       !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
(  245)       !!----------------------------------------------------------------------
(  246)       !! * Arguments
(  247)       CHARACTER(len=1), INTENT( in ) ::   &
(  248)          cd_type       ! nature of pt2d grid-point
(  249)          !             !   = T , U , V , F or W  gridpoints
(  250)          !             !   = I sea-ice U-V gridpoint (= F ocean grid point with indice shift)
(  251)       REAL(wp), INTENT( in ) ::   &
(  252)          psgn          ! control of the sign change
(  253)          !             !   =-1 , the sign is modified following the type of b.c. used
(  254)          !             !   = 1 , no sign change
(  255)       REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ) ::   &
(  256)          pt2d          ! 2D array on which the boundary condition is applied
(  257) 
(  258)       !! * Local declarations
(  259)       INTEGER  ::   ji
(  260)       INTEGER  ::   ijt, iju
(  261)       !!----------------------------------------------------------------------
(  262)       !!  OPA 8.5, LODYC-IPSL (2002)
(  263)       !!----------------------------------------------------------------------
(  264)       
(  265)       
(  266)       !                                        ! East-West boundaries
(  267)       !                                        ! ====================
(  268)       SELECT CASE ( nperio )
(  269) 
(  270)       CASE ( 1 , 4 , 6 )                       ! * cyclic east-west
(  271)          pt2d( 1 ,:) = pt2d(jpim1,:)
(  272)          pt2d(jpi,:) = pt2d(  2  ,:)
(  273) 
(  274)       CASE DEFAULT                             ! * closed 
(  275)          SELECT CASE ( cd_type )
(  276)          CASE ( 'T' , 'U' , 'V' , 'W' )                ! T-, U-, V-, W-points
(  277)             pt2d( 1 ,:) = 0.e0
(  278)             pt2d(jpi,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 6

(  279)          CASE ( 'F' )                                  ! F-point, ice U-V point
(  280)             pt2d(jpi,:) = 0.e0 
(  281)          CASE ( 'I' )                                  ! F-point, ice U-V point
(  282)             pt2d( 1 ,:) = 0.e0 
(  283)             pt2d(jpi,:) = 0.e0 
(  284)          END SELECT
(  285) 
(  286)       END SELECT
(  287) 
(  288)       !                                        ! North-South boundaries
(  289)       !                                        ! ======================
(  290)       SELECT CASE ( nperio )
(  291) 
(  292)       CASE ( 2 )                               ! * South symmetric
(  293) 
(  294)          SELECT CASE ( cd_type )
(  295)          CASE ( 'T' , 'U' , 'W' )                      ! T-, U-, W-points
(  296)             pt2d(:, 1 ) = pt2d(:,3)
(  297)             pt2d(:,jpj) = 0.e0
(  298)          CASE ( 'V' , 'F' , 'I' )                      ! V-, F-points, ice U-V point
(  299)             pt2d(:, 1 ) = psgn * pt2d(:,2)
(  300)             pt2d(:,jpj) = 0.e0
(  301)          END SELECT
(  302) 
(  303)       CASE ( 3 , 4 )                           ! * North fold  T-point pivot
(  304) 
(  305)          pt2d( 1 , 1 ) = 0.e0        !!!!!  bug gm ??? !Edmee
(  306)          pt2d( 1 ,jpj) = 0.e0
(  307)          pt2d(jpi,jpj) = 0.e0
(  308) 
(  309)          SELECT CASE ( cd_type )
(  310) 
(  311)          CASE ( 'T' , 'W' )                            ! T-, W-point
(  312)             DO ji = 2, jpi
(  313)                ijt = jpi-ji+2
(  314)                pt2d(ji, 1 ) = 0.e0
(  315)                pt2d(ji,jpj) = psgn * pt2d(ijt,jpj-2)
(  316)             END DO
(  317)             DO ji = jpi/2+1, jpi
(  318)                ijt = jpi-ji+2
(  319)                pt2d(ji,jpjm1) = psgn * pt2d(ijt,jpjm1)
(  320)             END DO
(  321) 
(  322)          CASE ( 'U' )                                  ! U-point
(  323)             DO ji = 1, jpi-1
(  324)                iju = jpi-ji+1
(  325)                pt2d(ji, 1 ) = 0.e0
(  326)                pt2d(ji,jpj) = psgn * pt2d(iju,jpj-2)
(  327)             END DO
(  328)             DO ji = jpi/2, jpi-1
(  329)                iju = jpi-ji+1
(  330)                pt2d(ji,jpjm1) = psgn * pt2d(iju,jpjm1)
(  331)             END DO
(  332) 
(  333)          CASE ( 'V' )                                  ! V-point
(  334)             DO ji = 2, jpi
(  335)                ijt = jpi-ji+2
(  336)                pt2d(ji, 1   ) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 7

(  337)                pt2d(ji,jpj-1) = psgn * pt2d(ijt,jpj-2)
(  338)                pt2d(ji,jpj  ) = psgn * pt2d(ijt,jpj-3)
(  339)             END DO
(  340) 
(  341)          CASE ( 'F' )                                  ! F-point
(  342)             DO ji = 1, jpi-1
(  343)                iju = jpi - ji + 1
(  344)                pt2d(ji,jpj-1) = psgn * pt2d(iju,jpj-2)
(  345)                pt2d(ji,jpj  ) = psgn * pt2d(iju,jpj-3)
(  346)             END DO
(  347) 
(  348)          CASE ( 'I' )                                  ! ice U-V point
(  349)             pt2d(:, 1 ) = 0.e0
(  350)             pt2d(2,jpj) = psgn * pt2d(3,jpj-1)
(  351)             DO ji = 3, jpi
(  352)                iju = jpi - ji + 3
(  353)                pt2d(ji,jpj) = psgn * pt2d(iju,jpj-1)
(  354)             END DO
(  355) 
(  356)          END SELECT
(  357) 
(  358)       CASE ( 5 , 6 )                           ! * North fold  F-point pivot
(  359) 
(  360)          pt2d( 1 , 1 ) = 0.e0           !!bug  ???
(  361)          pt2d( 1 ,jpj) = 0.e0
(  362)          pt2d(jpi,jpj) = 0.e0
(  363) 
(  364)          SELECT CASE ( cd_type )
(  365) 
(  366)          CASE ( 'T' , 'W' )                            ! T-, W-point
(  367)             DO ji = 1, jpi
(  368)                ijt = jpi-ji+1
(  369)                pt2d(ji, 1 ) = 0.e0
(  370)                pt2d(ji,jpj) = psgn * pt2d(ijt,jpj-1)
(  371)             END DO
(  372) 
(  373)          CASE ( 'U' )                                  ! U-point
(  374)             DO ji = 1, jpi-1
(  375)                iju = jpi-ji
(  376)                pt2d(ji, 1 ) = 0.e0
(  377)                pt2d(ji,jpj) = psgn * pt2d(iju,jpj-1)
(  378)             END DO
(  379) 
(  380)          CASE ( 'V' )                                  ! V-point
(  381)             DO ji = 1, jpi
(  382)                ijt = jpi-ji+1
(  383)                pt2d(ji, 1 ) = 0.e0
(  384)                pt2d(ji,jpj) = psgn * pt2d(ijt,jpj-2)
(  385)             END DO
(  386)             DO ji = jpi/2+1, jpi
(  387)                ijt = jpi-ji+1
(  388)                pt2d(ji,jpjm1) = psgn * pt2d(ijt,jpjm1)
(  389)             END DO
(  390) 
(  391)          CASE ( 'F' )                                  ! F-point
(  392)             DO ji = 1, jpi-1
(  393)                iju = jpi-ji
(  394)                pt2d(ji,jpj  ) = psgn * pt2d(iju,jpj-2)






PGF90 (Version     10.5)          09/23/2010  09:49:02      page 8

(  395)             END DO
(  396)             DO ji = jpi/2+1, jpi-1
(  397)                iju = jpi-ji
(  398)                pt2d(ji,jpjm1) = psgn * pt2d(iju,jpjm1)
(  399)             END DO
(  400) 
(  401)          CASE ( 'I' )                                  ! ice U-V point
(  402)             pt2d( : , 1 ) = 0.e0
(  403)             pt2d( 2 ,jpj) = 0.e0
(  404)             DO ji = 2 , jpim1
(  405)                ijt = jpi - ji + 2
(  406)                pt2d(ji,jpj)= 0.5 * ( pt2d(ji,jpjm1) + psgn * pt2d(ijt,jpjm1) )
(  407)             END DO
(  408) 
(  409)          END SELECT
(  410) 
(  411)       CASE DEFAULT                             ! * closed
(  412) 
(  413)          SELECT CASE ( cd_type )
(  414)          CASE ( 'T' , 'U' , 'V' , 'W' )                ! T-, U-, V-, W-points
(  415)             pt2d(:, 1 ) = 0.e0
(  416)             pt2d(:,jpj) = 0.e0
(  417)          CASE ( 'F' )                                  ! F-point
(  418)             pt2d(:,jpj) = 0.e0
(  419)          CASE ( 'I' )                                  ! ice U-V point
(  420)             pt2d(:, 1 ) = 0.e0
(  421)             pt2d(:,jpj) = 0.e0
(  422)          END SELECT
(  423) 
(  424)       END SELECT
(  425) 
(  426)    END SUBROUTINE lbc_lnk_2d
(  427) 
(  428) #endif
(  429) 
(  430)    !!======================================================================
(  431) END MODULE lbclnk



























PGF90 (Version     10.5)          09/23/2010  09:49:02      page 9

(    1) # 1 "lbclnk.F90"
(    1) MODULE lbclnk
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  lbclnk  ***
(    4)    !! Ocean        : lateral boundary conditions
(    5)    !!=====================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_mpp_mpi'     OR      MPI massively parallel processing library
(    9)    !!   'key_mpp_shmem'         SHMEM massively parallel processing library
(   10)    !!----------------------------------------------------------------------
(   11)    !!----------------------------------------------------------------------
(   12)    !!   lbc_lnk      : generic interface for mpp_lnk_3d and mpp_lnk_2d
(   13)    !!                  routines defined in lib_mpp
(   14)    !!   lbc_lnk_e    : generic interface for mpp_lnk_2d_e
(   15)    !!                   routinee defined in lib_mpp
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE lib_mpp          ! distributed memory computing library
(   19) 
(   20)    INTERFACE lbc_lnk
(   21)       MODULE PROCEDURE mpp_lnk_3d, mpp_lnk_2d
(   22)    END INTERFACE
(   23) 
(   24)    INTERFACE lbc_lnk_e
(   25)       MODULE PROCEDURE mpp_lnk_2d_e
(   26)    END INTERFACE
(   27) 
(   28)    PUBLIC lbc_lnk       ! ocean lateral boundary conditions
(   29)    PUBLIC lbc_lnk_e
(   30)    !!----------------------------------------------------------------------
(   31) 
(   32) 
(   33) # 430
(  430)    !!======================================================================
(  431) END MODULE lbclnk
