


PGF90 (Version     10.5)          09/23/2010  09:50:16      page 1

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

Filename: tranpc.F90

(    1) MODULE tranpc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  tranpc  ***
(    4)    !! Ocean active tracers:  non penetrative convection scheme
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_npc      : apply the non penetrative convection scheme
(    9)    !!   tra_npc_init : initialization and control of the scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers 
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE trdmod          ! ocean active tracer trends
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE eosbn2          ! equation of state (eos routine) 
(   17)    USE lbclnk          ! lateral boundary conditions (or mpp link)
(   18)    USE in_out_manager  ! I/O manager
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Routine accessibility
(   24)    PUBLIC tra_npc      ! routine called by step.F90
(   25) 
(   26)    !! * Module variable
(   27)    INTEGER ::       &
(   28)       nnpc1 =   1,  &  ! nnpc1   non penetrative convective scheme frequency
(   29)       nnpc2 =  15      ! nnpc2   non penetrative convective scheme print frequency
(   30) 
(   31)    !! * Substitutions
(   32) #  include "domzgr_substitute.h90"
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tranpc.F90,v 1.3 2005/03/27 18:35:20 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE tra_npc( kt )
(   42)       !!----------------------------------------------------------------------
(   43)       !!                  ***  ROUTINE tranpc  ***
(   44)       !!
(   45)       !! ** Purpose :   Non penetrative convective adjustment scheme. solve 
(   46)       !!      the static instability of the water column (now, after the swap) 






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 2

(   47)       !!      while conserving heat and salt contents.
(   48)       !!
(   49)       !! ** Method  :   The algorithm used converges in a maximium of jpk 
(   50)       !!      iterations. instabilities are treated when the vertical density
(   51)       !!      gradient is less than 1.e-5.
(   52)       !!
(   53)       !!      'key_trdtra' defined: the trend associated with this
(   54)       !!                               algorithm is saved.
(   55)       !!
(   56)       !!      macro-tasked on vertical slab (jj-loop)
(   57)       !!
(   58)       !! ** Action  : - (tn,sn) after the application od the npc scheme
(   59)       !!              - save the associated trends (ttrd,strd) ('key_trdtra')
(   60)       !!
(   61)       !! References :
(   62)       !!      Madec, et al., 1991, JPO, 21, 9, 1349-1371.
(   63)       !!
(   64)       !! History :
(   65)       !!   1.0  !  90-09  (G. Madec)  Original code
(   66)       !!        !  91-11  (G. Madec)
(   67)       !!        !  92-06  (M. Imbard)  periodic conditions on t and s
(   68)       !!        !  93-03  (M. Guyon)  symetrical conditions 
(   69)       !!        !  96-01  (G. Madec)  statement function for e3
(   70)       !!                                  suppression of common work arrays
(   71)       !!   8.5  !  02-06  (G. Madec)  free form F90
(   72)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   73)       !!----------------------------------------------------------------------
(   74)       !! * Modules used     
(   75)       USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
(   76)                          ztdsa => va      ! use va as 3D workspace   
(   77) 
(   78)       !! * Arguments
(   79)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   80) 
(   81)       !! * Local declarations
(   82)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(   83)       INTEGER ::   &
(   84)          inpcc ,                        & ! number of statically instable water column
(   85)          inpci ,                        & ! number of iteration for npc scheme
(   86)          jiter, jkdown, jkp,            & ! ???
(   87)          ikbot, ik, ikup, ikdown          ! ???
(   88)       REAL(wp) ::   &                     ! temporary arrays
(   89)          ze3tot, zta, zsa, zraua, ze3dwn
(   90)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(   91)          zwx, zwy, zwz                    ! temporary arrays
(   92)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   93)          zrhop                            ! temporary arrays
(   94)       !!----------------------------------------------------------------------
(   95) 
(   96)       IF( kt == nit000  )   CALL tra_npc_init
(   97) 
(   98) 
(   99)       IF( MOD( kt, nnpc1 ) == 0 ) THEN
(  100) 
(  101)          inpcc = 0
(  102)          inpci = 0
(  103) 
(  104)          ! 0. Potential density






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 3

(  105)          ! --------------------
(  106) 
(  107)          CALL eos( tn, sn, rhd, zrhop )
(  108) 
(  109)          ! Save tn and sn trends
(  110)          IF( l_trdtra )   THEN
(  111)             ztdta(:,:,:) = tn(:,:,:) 
(  112)             ztdsa(:,:,:) = sn(:,:,:) 
(  113)          ENDIF
(  114) 
(  115)          !                                                ! ===============
(  116)          DO jj = 1, jpj                                   !  Vertical slab
(  117)             !                                             ! ===============
(  118) 
(  119)             ! 1. Static instability pointer 
(  120)             ! -----------------------------
(  121) 
(  122)             DO jk = 1, jpkm1
(  123)                DO ji = 1, jpi
(  124)                   zwx(ji,jk) = ( zrhop(ji,jj,jk) - zrhop(ji,jj,jk+1) ) * tmask(ji,jj,jk+1)
(  125)                END DO
(  126)             END DO
(  127) 
(  128)             ! 1.1 do not consider the boundary points
(  129) 
(  130)             ! even if east-west cyclic b. c. do not considere ji=1 or jpi
(  131)             DO jk = 1, jpkm1
(  132)                zwx( 1 ,jk) = 0.e0
(  133)                zwx(jpi,jk) = 0.e0
(  134)             END DO
(  135)             ! even if south-symmetric b. c. used, do not considere jj=1
(  136)             IF( jj == 1 ) zwx(:,:) = 0.e0
(  137) 
(  138)             DO jk = 1, jpkm1
(  139)                DO ji = 1, jpi
(  140)                   zwx(ji,jk) = 1.
(  141)                   IF( zwx(ji,jk) < 1.e-5 ) zwx(ji,jk)=0.
(  142)                END DO
(  143)             END DO
(  144) 
(  145)             zwy(:,1) = 0.
(  146)             DO ji = 1, jpi
(  147)                DO jk = 1, jpkm1
(  148)                   zwy(ji,1) = zwy(ji,1) + zwx(ji,jk)
(  149)                END DO
(  150)             END DO
(  151) 
(  152)             zwz(1,1) = 0.
(  153)             DO ji = 1, jpi
(  154)                zwz(1,1) = zwz(1,1) + zwy(ji,1)
(  155)             END DO
(  156) 
(  157)             inpcc = inpcc + NINT( zwz(1,1) )
(  158) 
(  159) 
(  160)             ! 2. Vertical mixing for each instable portion of the density profil
(  161)             ! ------------------------------------------------------------------
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 4

(  163)             IF (zwz(1,1) /= 0.) THEN
(  164) 
(  165)                ! -->> the density profil is statically instable :
(  166) 
(  167)                DO ji = 1, jpi
(  168)                   IF( zwy(ji,1) /= 0. ) THEN
(  169) 
(  170)                      ! ikbot: ocean bottom level
(  171) 
(  172)                      ikbot = mbathy(ji,jj)
(  173) 
(  174)                      ! vertical iteration
(  175) 
(  176)                      DO jiter = 1, jpk
(  177) 
(  178)                         ! search of ikup : the first static instability from the sea surface
(  179) 
(  180)                         ik = 0
(  181) 220                     CONTINUE
(  182)                         ik = ik + 1
(  183)                         IF( ik >= ikbot-1 ) GO TO 200
(  184)                         zwx(ji,ik) = zrhop(ji,jj,ik) - zrhop(ji,jj,ik+1)
(  185)                         IF( zwx(ji,ik) <= 0. ) GO TO 220
(  186)                         ikup = ik
(  187)                         ! the density profil is instable below ikup
(  188)    
(  189)                         ! ikdown : bottom of the instable portion of the density profil
(  190) 
(  191)                         ! search of ikdown and vertical mixing from ikup to ikdown
(  192) 
(  193)                         ze3tot= fse3t(ji,jj,ikup)
(  194)                         zta   = tn   (ji,jj,ikup)
(  195)                         zsa   = sn   (ji,jj,ikup)
(  196)                         zraua = zrhop(ji,jj,ikup)
(  197) 
(  198)                         DO jkdown = ikup+1, ikbot-1
(  199)                            IF( zraua <= zrhop(ji,jj,jkdown) ) THEN
(  200)                               ikdown = jkdown
(  201)                               GO TO 240
(  202)                            ENDIF
(  203)                            ze3dwn =  fse3t(ji,jj,jkdown)
(  204)                            ze3tot =  ze3tot + ze3dwn
(  205)                            zta   = ( zta*(ze3tot-ze3dwn) + tn(ji,jj,jkdown)*ze3dwn )/ze3tot
(  206)                            zsa   = ( zsa*(ze3tot-ze3dwn) + sn(ji,jj,jkdown)*ze3dwn )/ze3tot
(  207)                            zraua = ( zraua*(ze3tot-ze3dwn) + zrhop(ji,jj,jkdown)*ze3dwn )/ze3tot
(  208)                            inpci = inpci+1
(  209)                         END DO
(  210)                         ikdown = ikbot-1
(  211) 240                     CONTINUE
(  212) 
(  213)                         DO jkp = ikup, ikdown-1
(  214)                            tn(ji,jj,jkp) = zta
(  215)                            sn(ji,jj,jkp) = zsa
(  216)                            zrhop(ji,jj,jkp) = zraua
(  217)                         END DO
(  218)                         IF (ikdown == ikbot-1 .AND. zraua >= zrhop(ji,jj,ikdown) ) THEN
(  219)                            tn(ji,jj,ikdown) = zta
(  220)                            sn(ji,jj,ikdown) = zsa






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 5

(  221)                            zrhop(ji,jj,ikdown) = zraua
(  222)                         ENDIF
(  223) 
(  224)                      END DO
(  225)                   ENDIF
(  226) 200               CONTINUE
(  227)                END DO
(  228) 
(  229)                ! <<-- no more static instability on slab jj
(  230) 
(  231)             ENDIF
(  232)             !                                             ! ===============
(  233)          END DO                                           !   End of slab
(  234)          !                                                ! ===============
(  235) 
(  236) 
(  237)          ! save the trends for diagnostic
(  238)          ! Non penetrative mixing trends
(  239)          IF( l_trdtra )   THEN
(  240)             ztdta(:,:,:) = tn(:,:,:) - ztdta(:,:,:)
(  241)             ztdsa(:,:,:) = sn(:,:,:) - ztdsa(:,:,:)
(  242) 
(  243)             CALL trd_mod(ztdta, ztdsa, jpttdnpc, 'TRA', kt)
(  244)          ENDIF
(  245)       
(  246)          ! Lateral boundary conditions on ( tn, sn )   ( Unchanged sign)
(  247)          ! ------------------------------============
(  248)          CALL lbc_lnk( tn, 'T', 1. )
(  249)          CALL lbc_lnk( sn, 'T', 1. )
(  250)       
(  251) 
(  252)          !  2. non penetrative convective scheme statistics
(  253)          !  -----------------------------------------------
(  254) 
(  255)          IF( nnpc2 /= 0 .AND. MOD( kt, nnpc2 ) == 0 ) THEN
(  256)             IF(lwp) WRITE(numout,*)' kt=',kt, ' number of statically instable',   &
(  257)                ' water column : ',inpcc, ' number of iteration : ',inpci
(  258)          ENDIF
(  259) 
(  260)       ENDIF
(  261)       
(  262)    END SUBROUTINE tra_npc
(  263) 
(  264) 
(  265)    SUBROUTINE tra_npc_init
(  266)       !!----------------------------------------------------------------------
(  267)       !!                  ***  ROUTINE tra_npc_init  ***
(  268)       !!                   
(  269)       !! ** Purpose :   initializations of the non-penetrative adjustment scheme
(  270)       !!
(  271)       !! History :
(  272)       !!   8.5  !  02-12  (G. Madec)  F90 : free form
(  273)       !!----------------------------------------------------------------------
(  274)       !! * Namelist
(  275)       NAMELIST/namnpc/ nnpc1, nnpc2
(  276)       !!----------------------------------------------------------------------
(  277) 
(  278)       ! Namelist namzdf : vertical diffusion






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 6

(  279)       REWIND( numnam )
(  280)       READ  ( numnam, namnpc )
(  281) 
(  282)       ! Parameter print
(  283)       ! ---------------
(  284)       IF(lwp) THEN
(  285)          WRITE(numout,*)
(  286)          WRITE(numout,*) 'tra_npc_init : Non Penetrative Convection (npc) scheme'
(  287)          WRITE(numout,*) '~~~~~~~~~~~~'
(  288)          WRITE(numout,*) '          Namelist namnpc : set npc scheme parameters'
(  289)          WRITE(numout,*)
(  290)          WRITE(numout,*) '             npc scheme frequency           nnpc1  = ', nnpc1
(  291)          WRITE(numout,*) '             npc scheme print frequency     nnpc2  = ', nnpc2
(  292)          WRITE(numout,*)
(  293)       ENDIF
(  294) 
(  295) 
(  296)       ! Parameter controls
(  297)       ! ------------------
(  298)       IF ( nnpc1 == 0 ) THEN
(  299)           IF(lwp) WRITE(numout,cform_war)
(  300)           IF(lwp) WRITE(numout,*) '             nnpc1 = ', nnpc1, ' is forced to 1'
(  301)           nnpc1 = 1
(  302)           nwarn = nwarn + 1
(  303)       ENDIF
(  304)       
(  305)    END SUBROUTINE tra_npc_init
(  306) 
(  307)    !!======================================================================
(  308) END MODULE tranpc


































PGF90 (Version     10.5)          09/23/2010  09:50:16      page 7

(    1) # 1 "tranpc.F90"
(    1) MODULE tranpc
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  tranpc  ***
(    4)    !! Ocean active tracers:  non penetrative convection scheme
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_npc      : apply the non penetrative convection scheme
(    9)    !!   tra_npc_init : initialization and control of the scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers 
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE trdmod          ! ocean active tracer trends
(   15)    USE trdmod_oce      ! ocean variables trends
(   16)    USE eosbn2          ! equation of state (eos routine) 
(   17)    USE lbclnk          ! lateral boundary conditions (or mpp link)
(   18)    USE in_out_manager  ! I/O manager
(   19) 
(   20)    IMPLICIT NONE
(   21)    PRIVATE
(   22) 
(   23)    !! * Routine accessibility
(   24)    PUBLIC tra_npc      ! routine called by step.F90
(   25) 
(   26)    !! * Module variable
(   27)    INTEGER ::       &
(   28)       nnpc1 =   1,  &  ! nnpc1   non penetrative convective scheme frequency
(   29)       nnpc2 =  15      ! nnpc2   non penetrative convective scheme print frequency
(   30) 
(   31)    !! * Substitutions
(   32) # 1 "./domzgr_substitute.h90"
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
(   48) # 33 "tranpc.F90"
(   33) # 33
(   33)    !!----------------------------------------------------------------------
(   34)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/tranpc.F90,v 1.3 2005/03/27 18:35:20 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE tra_npc( kt )






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 8

(   42)       !!----------------------------------------------------------------------
(   43)       !!                  ***  ROUTINE tranpc  ***
(   44)       !!
(   45)       !! ** Purpose :   Non penetrative convective adjustment scheme. solve 
(   46)       !!      the static instability of the water column (now, after the swap) 
(   47)       !!      while conserving heat and salt contents.
(   48)       !!
(   49)       !! ** Method  :   The algorithm used converges in a maximium of jpk 
(   50)       !!      iterations. instabilities are treated when the vertical density
(   51)       !!      gradient is less than 1.e-5.
(   52)       !!
(   53)       !!      'key_trdtra' defined: the trend associated with this
(   54)       !!                               algorithm is saved.
(   55)       !!
(   56)       !!      macro-tasked on vertical slab (jj-loop)
(   57)       !!
(   58)       !! ** Action  : - (tn,sn) after the application od the npc scheme
(   59)       !!              - save the associated trends (ttrd,strd) ('key_trdtra')
(   60)       !!
(   61)       !! References :
(   62)       !!      Madec, et al., 1991, JPO, 21, 9, 1349-1371.
(   63)       !!
(   64)       !! History :
(   65)       !!   1.0  !  90-09  (G. Madec)  Original code
(   66)       !!        !  91-11  (G. Madec)
(   67)       !!        !  92-06  (M. Imbard)  periodic conditions on t and s
(   68)       !!        !  93-03  (M. Guyon)  symetrical conditions 
(   69)       !!        !  96-01  (G. Madec)  statement function for e3
(   70)       !!                                  suppression of common work arrays
(   71)       !!   8.5  !  02-06  (G. Madec)  free form F90
(   72)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   73)       !!----------------------------------------------------------------------
(   74)       !! * Modules used     
(   75)       USE oce, ONLY :    ztdta => ua,   & ! use ua as 3D workspace   
(   76)                          ztdsa => va      ! use va as 3D workspace   
(   77) 
(   78)       !! * Arguments
(   79)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   80) 
(   81)       !! * Local declarations
(   82)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(   83)       INTEGER ::   &
(   84)          inpcc ,                        & ! number of statically instable water column
(   85)          inpci ,                        & ! number of iteration for npc scheme
(   86)          jiter, jkdown, jkp,            & ! ???
(   87)          ikbot, ik, ikup, ikdown          ! ???
(   88)       REAL(wp) ::   &                     ! temporary arrays
(   89)          ze3tot, zta, zsa, zraua, ze3dwn
(   90)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(   91)          zwx, zwy, zwz                    ! temporary arrays
(   92)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   93)          zrhop                            ! temporary arrays
(   94)       !!----------------------------------------------------------------------
(   95) 
(   96)       IF( kt == nit000  )   CALL tra_npc_init
(   97) 
(   98) 
(   99)       IF( MOD( kt, nnpc1 ) == 0 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 9

(  100) 
(  101)          inpcc = 0
(  102)          inpci = 0
(  103) 
(  104)          ! 0. Potential density
(  105)          ! --------------------
(  106) 
(  107)          CALL eos( tn, sn, rhd, zrhop )
(  108) 
(  109)          ! Save tn and sn trends
(  110)          IF( l_trdtra )   THEN
(  111)             ztdta(:,:,:) = tn(:,:,:) 
(  112)             ztdsa(:,:,:) = sn(:,:,:) 
(  113)          ENDIF
(  114) 
(  115)          !                                                ! ===============
(  116)          DO jj = 1, jpj                                   !  Vertical slab
(  117)             !                                             ! ===============
(  118) 
(  119)             ! 1. Static instability pointer 
(  120)             ! -----------------------------
(  121) 
(  122)             DO jk = 1, jpkm1
(  123)                DO ji = 1, jpi
(  124)                   zwx(ji,jk) = ( zrhop(ji,jj,jk) - zrhop(ji,jj,jk+1) ) * tmask(ji,jj,jk+1)
(  125)                END DO
(  126)             END DO
(  127) 
(  128)             ! 1.1 do not consider the boundary points
(  129) 
(  130)             ! even if east-west cyclic b. c. do not considere ji=1 or jpi
(  131)             DO jk = 1, jpkm1
(  132)                zwx( 1 ,jk) = 0.e0
(  133)                zwx(jpi,jk) = 0.e0
(  134)             END DO
(  135)             ! even if south-symmetric b. c. used, do not considere jj=1
(  136)             IF( jj == 1 ) zwx(:,:) = 0.e0
(  137) 
(  138)             DO jk = 1, jpkm1
(  139)                DO ji = 1, jpi
(  140)                   zwx(ji,jk) = 1.
(  141)                   IF( zwx(ji,jk) < 1.e-5 ) zwx(ji,jk)=0.
(  142)                END DO
(  143)             END DO
(  144) 
(  145)             zwy(:,1) = 0.
(  146)             DO ji = 1, jpi
(  147)                DO jk = 1, jpkm1
(  148)                   zwy(ji,1) = zwy(ji,1) + zwx(ji,jk)
(  149)                END DO
(  150)             END DO
(  151) 
(  152)             zwz(1,1) = 0.
(  153)             DO ji = 1, jpi
(  154)                zwz(1,1) = zwz(1,1) + zwy(ji,1)
(  155)             END DO
(  156) 
(  157)             inpcc = inpcc + NINT( zwz(1,1) )






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 10

(  158) 
(  159) 
(  160)             ! 2. Vertical mixing for each instable portion of the density profil
(  161)             ! ------------------------------------------------------------------
(  162) 
(  163)             IF (zwz(1,1) /= 0.) THEN
(  164) 
(  165)                ! -->> the density profil is statically instable :
(  166) 
(  167)                DO ji = 1, jpi
(  168)                   IF( zwy(ji,1) /= 0. ) THEN
(  169) 
(  170)                      ! ikbot: ocean bottom level
(  171) 
(  172)                      ikbot = mbathy(ji,jj)
(  173) 
(  174)                      ! vertical iteration
(  175) 
(  176)                      DO jiter = 1, jpk
(  177) 
(  178)                         ! search of ikup : the first static instability from the sea surface
(  179) 
(  180)                         ik = 0
(  181) 220                     CONTINUE
(  182)                         ik = ik + 1
(  183)                         IF( ik >= ikbot-1 ) GO TO 200
(  184)                         zwx(ji,ik) = zrhop(ji,jj,ik) - zrhop(ji,jj,ik+1)
(  185)                         IF( zwx(ji,ik) <= 0. ) GO TO 220
(  186)                         ikup = ik
(  187)                         ! the density profil is instable below ikup
(  188)    
(  189)                         ! ikdown : bottom of the instable portion of the density profil
(  190) 
(  191)                         ! search of ikdown and vertical mixing from ikup to ikdown
(  192) 
(  193)                         ze3tot= e3t(ikup)
(  194)                         zta   = tn   (ji,jj,ikup)
(  195)                         zsa   = sn   (ji,jj,ikup)
(  196)                         zraua = zrhop(ji,jj,ikup)
(  197) 
(  198)                         DO jkdown = ikup+1, ikbot-1
(  199)                            IF( zraua <= zrhop(ji,jj,jkdown) ) THEN
(  200)                               ikdown = jkdown
(  201)                               GO TO 240
(  202)                            ENDIF
(  203)                            ze3dwn =  e3t(jkdown)
(  204)                            ze3tot =  ze3tot + ze3dwn
(  205)                            zta   = ( zta*(ze3tot-ze3dwn) + tn(ji,jj,jkdown)*ze3dwn )/ze3tot
(  206)                            zsa   = ( zsa*(ze3tot-ze3dwn) + sn(ji,jj,jkdown)*ze3dwn )/ze3tot
(  207)                            zraua = ( zraua*(ze3tot-ze3dwn) + zrhop(ji,jj,jkdown)*ze3dwn )/ze3tot
(  208)                            inpci = inpci+1
(  209)                         END DO
(  210)                         ikdown = ikbot-1
(  211) 240                     CONTINUE
(  212) 
(  213)                         DO jkp = ikup, ikdown-1
(  214)                            tn(ji,jj,jkp) = zta
(  215)                            sn(ji,jj,jkp) = zsa






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 11

(  216)                            zrhop(ji,jj,jkp) = zraua
(  217)                         END DO
(  218)                         IF (ikdown == ikbot-1 .AND. zraua >= zrhop(ji,jj,ikdown) ) THEN
(  219)                            tn(ji,jj,ikdown) = zta
(  220)                            sn(ji,jj,ikdown) = zsa
(  221)                            zrhop(ji,jj,ikdown) = zraua
(  222)                         ENDIF
(  223) 
(  224)                      END DO
(  225)                   ENDIF
(  226) 200               CONTINUE
(  227)                END DO
(  228) 
(  229)                ! <<-- no more static instability on slab jj
(  230) 
(  231)             ENDIF
(  232)             !                                             ! ===============
(  233)          END DO                                           !   End of slab
(  234)          !                                                ! ===============
(  235) 
(  236) 
(  237)          ! save the trends for diagnostic
(  238)          ! Non penetrative mixing trends
(  239)          IF( l_trdtra )   THEN
(  240)             ztdta(:,:,:) = tn(:,:,:) - ztdta(:,:,:)
(  241)             ztdsa(:,:,:) = sn(:,:,:) - ztdsa(:,:,:)
(  242) 
(  243)             CALL trd_mod(ztdta, ztdsa, jpttdnpc, 'TRA', kt)
(  244)          ENDIF
(  245)       
(  246)          ! Lateral boundary conditions on ( tn, sn )   ( Unchanged sign)
(  247)          ! ------------------------------============
(  248)          CALL lbc_lnk( tn, 'T', 1. )
(  249)          CALL lbc_lnk( sn, 'T', 1. )
(  250)       
(  251) 
(  252)          !  2. non penetrative convective scheme statistics
(  253)          !  -----------------------------------------------
(  254) 
(  255)          IF( nnpc2 /= 0 .AND. MOD( kt, nnpc2 ) == 0 ) THEN
(  256)             IF(lwp) WRITE(numout,*)' kt=',kt, ' number of statically instable',   &
(  257)                ' water column : ',inpcc, ' number of iteration : ',inpci
(  258)          ENDIF
(  259) 
(  260)       ENDIF
(  261)       
(  262)    END SUBROUTINE tra_npc
(  263) 
(  264) 
(  265)    SUBROUTINE tra_npc_init
(  266)       !!----------------------------------------------------------------------
(  267)       !!                  ***  ROUTINE tra_npc_init  ***
(  268)       !!                   
(  269)       !! ** Purpose :   initializations of the non-penetrative adjustment scheme
(  270)       !!
(  271)       !! History :
(  272)       !!   8.5  !  02-12  (G. Madec)  F90 : free form
(  273)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:16      page 12

(  274)       !! * Namelist
(  275)       NAMELIST/namnpc/ nnpc1, nnpc2
(  276)       !!----------------------------------------------------------------------
(  277) 
(  278)       ! Namelist namzdf : vertical diffusion
(  279)       REWIND( numnam )
(  280)       READ  ( numnam, namnpc )
(  281) 
(  282)       ! Parameter print
(  283)       ! ---------------
(  284)       IF(lwp) THEN
(  285)          WRITE(numout,*)
(  286)          WRITE(numout,*) 'tra_npc_init : Non Penetrative Convection (npc) scheme'
(  287)          WRITE(numout,*) '~~~~~~~~~~~~'
(  288)          WRITE(numout,*) '          Namelist namnpc : set npc scheme parameters'
(  289)          WRITE(numout,*)
(  290)          WRITE(numout,*) '             npc scheme frequency           nnpc1  = ', nnpc1
(  291)          WRITE(numout,*) '             npc scheme print frequency     nnpc2  = ', nnpc2
(  292)          WRITE(numout,*)
(  293)       ENDIF
(  294) 
(  295) 
(  296)       ! Parameter controls
(  297)       ! ------------------
(  298)       IF ( nnpc1 == 0 ) THEN
(  299)           IF(lwp) WRITE(numout,cform_war)
(  300)           IF(lwp) WRITE(numout,*) '             nnpc1 = ', nnpc1, ' is forced to 1'
(  301)           nnpc1 = 1
(  302)           nwarn = nwarn + 1
(  303)       ENDIF
(  304)       
(  305)    END SUBROUTINE tra_npc_init
(  306) 
(  307)    !!======================================================================
(  308) END MODULE tranpc
