


PGF90 (Version     10.5)          09/23/2010  09:49:11      page 1

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

Filename: geo2ocean.F90

(    1) MODULE geo2ocean
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  geo2ocean  ***
(    4)    !! Ocean mesh    :  ???
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   repcmo      : 
(    9)    !!   angle       :
(   10)    !!   geo2oce     :
(   11)    !!   repere      :   old routine suppress it ???
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE dom_oce         ! mesh and scale factors
(   15)    USE phycst          ! physical constants
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18) 
(   19)    IMPLICIT NONE
(   20) 
(   21)    !! * Accessibility
(   22)    PRIVATE
(   23)    PUBLIC repcmo   ! routine called by ???.F90
(   24)    PUBLIC geo2oce  ! routine called by ???.F90
(   25)    PUBLIC repere   ! routine called by ???.F90
(   26) 
(   27)    !! * Module variables
(   28)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   29)       gsinu , gcosu ,   &  ! matrix element for change grid u (repcmo.F)
(   30)       gsinv , gcosv ,   &  ! matrix element for change grid v (repcmo.F)
(   31)       gsinus, gcosin      ! matrix element for change grid (repere.F)
(   32) 
(   33)   !! * Substitutions
(   34) #  include "vectopt_loop_substitute.h90"
(   35)    !!---------------------------------------------------------------------------------
(   36)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/geo2ocean.F90,v 1.4 2005/03/27 18:34:46 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!---------------------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE repcmo ( pxu1, pyu1, pxv1, pyv1,   &
(   44)                        px2 , py2 , kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE repcmo  ***






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 2

(   47)       !!
(   48)       !! ** Purpose :   Change vector componantes from a geographic grid to a
(   49)       !!      stretched coordinates grid.
(   50)       !!
(   51)       !! ** Method  :   Initialization of arrays at the first call.
(   52)       !!
(   53)       !! ** Action  : - px2 : first componante (defined at u point)
(   54)       !!              - py2 : second componante (defined at v point)
(   55)       !!
(   56)       !! History :
(   57)       !!   7.0  !  07-96  (O. Marti)  Original code
(   58)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   59)       !!----------------------------------------------------------------------
(   60)       !! * Arguments 
(   61)       INTEGER,  INTENT( in ) ::   &
(   62)          kt                ! ocean time-step
(   63)       REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   & 
(   64)          pxu1, pyu1,     & ! geographic vector componantes at u-point
(   65)          pxv1, pyv1        ! geographic vector componantes at v-point
(   66)       REAL(wp), INTENT( out ), DIMENSION(jpi,jpj) ::   & 
(   67)          px2,            & ! i-componante (defined at u-point)
(   68)          py2               ! j-componante (defined at v-point)
(   69)       !!----------------------------------------------------------------------
(   70) 
(   71) 
(   72)       ! Initialization of gsin* and gcos* at first call
(   73)       ! -----------------------------------------------
(   74) 
(   75)       IF( kt <= nit000 + 1 ) THEN
(   76)          IF(lwp) WRITE(numout,*)
(   77)          IF(lwp) WRITE(numout,*) 'repcmo : use the geographic to stretched'
(   78)          IF(lwp) WRITE(numout,*) ' ~~~~~   coordinate transformation'
(   79) 
(   80)          CALL angle       ! initialization of the transformation
(   81)       ENDIF
(   82)       
(   83)       ! Change from geographic to stretched coordinate
(   84)       ! ----------------------------------------------
(   85)       
(   86)       px2(:,:) = pxu1(:,:) * gcosu(:,:) + pyu1(:,:) * gsinu(:,:)
(   87)       py2(:,:) = pyv1(:,:) * gcosv(:,:) - pxv1(:,:) * gsinv(:,:)   
(   88)       
(   89)    END SUBROUTINE repcmo
(   90) 
(   91) 
(   92)    SUBROUTINE angle
(   93)       !!----------------------------------------------------------------------
(   94)       !!                  ***  ROUTINE angle  ***
(   95)       !! 
(   96)       !! ** Purpose :   Compute angles between model grid lines and the 
(   97)       !!      direction of the North
(   98)       !!
(   99)       !! ** Method  :
(  100)       !!
(  101)       !! ** Action  :   Compute (gsinu, gcosu, gsinv, gcosv) arrays: sinus and 
(  102)       !!      cosinus of the angle between the north-south axe and the 
(  103)       !!      j-direction at u and v-points
(  104)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 3

(  105)       !! History :
(  106)       !!   7.0  !  96-07  (O. Marti)  Original code
(  107)       !!   8.0  !  98-06  (G. Madec)
(  108)       !!   8.5  !  98-06  (G. Madec)  Free form, F90 + opt.
(  109)       !!----------------------------------------------------------------------
(  110)       !! * local declarations
(  111)       INTEGER ::   ji, jj      ! dummy loop indices
(  112) 
(  113)       REAL(wp) ::   &
(  114)          zlam, zphi,             &  ! temporary scalars
(  115)          zlan, zphh,             &  !    "         "
(  116)          zxnpu, zxnpv , znnpu,   &  !    "         "
(  117)          zynpu, zynpv , znnpv,   &  !    "         "
(  118)          zxffu, zmnpfu, zxffv,   &  !    "         "
(  119)          zyffu, zmnpfv, zyffv       !    "         "
(  120)       !!----------------------------------------------------------------------
(  121) 
(  122)       ! ============================= !
(  123)       ! Compute the cosinus and sinus !
(  124)       ! ============================= !
(  125)       ! (computation done on the north stereographic polar plane)
(  126) 
(  127)       DO jj = 2, jpj
(  128) !CDIR NOVERRCHK
(  129)          DO ji = fs_2, jpi   ! vector opt.
(  130) 
(  131)             ! north pole direction & modulous (at u-point)
(  132)             zlam = glamu(ji,jj)
(  133)             zphi = gphiu(ji,jj)
(  134)             zxnpu = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
(  135)             zynpu = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
(  136)             znnpu = zxnpu*zxnpu + zynpu*zynpu
(  137) 
(  138)             ! north pole direction & modulous (at v-point)
(  139)             zlam = glamv(ji,jj)
(  140)             zphi = gphiv(ji,jj)
(  141)             zxnpv = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
(  142)             zynpv = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
(  143)             znnpv = zxnpv*zxnpv + zynpv*zynpv
(  144) 
(  145)             ! j-direction: f-point segment direction (u-point)
(  146)             zlam = glamf(ji,jj  )
(  147)             zphi = gphif(ji,jj  )
(  148)             zlan = glamf(ji,jj-1)
(  149)             zphh = gphif(ji,jj-1)
(  150)             zxffu =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
(  151)                &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
(  152)             zyffu =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
(  153)                &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
(  154)             zmnpfu = SQRT ( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
(  155)             zmnpfu = MAX( zmnpfu, 1.e-14 )
(  156) 
(  157)             ! i-direction: f-point segment direction (v-point)
(  158)             zlam = glamf(ji  ,jj)
(  159)             zphi = gphif(ji  ,jj)
(  160)             zlan = glamf(ji-1,jj)
(  161)             zphh = gphif(ji-1,jj)
(  162)             zxffv =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 4

(  163)                &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
(  164)             zyffv =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
(  165)                &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
(  166)             zmnpfv = SQRT ( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
(  167)             zmnpfv = MAX( zmnpfv, 1.e-14 )
(  168) 
(  169)             ! cosinus and sinus using scalar and vectorial products
(  170)             gsinu(ji,jj) = ( zxnpu*zyffu - zynpu*zxffu ) / zmnpfu
(  171)             gcosu(ji,jj) = ( zxnpu*zxffu + zynpu*zyffu ) / zmnpfu
(  172) 
(  173)             ! cosinus and sinus using scalar and vectorial products
(  174)             ! (caution, rotation of 90 degres)
(  175)             gsinv(ji,jj) = ( zxnpv*zxffv + zynpv*zyffv ) / zmnpfv
(  176)             gcosv(ji,jj) =-( zxnpv*zyffv - zynpv*zxffv ) / zmnpfv
(  177) 
(  178)          END DO
(  179)       END DO
(  180) 
(  181)       ! =============== !
(  182)       ! Geographic mesh !
(  183)       ! =============== !
(  184) 
(  185)       DO jj = 2, jpj
(  186)          DO ji = fs_2, jpi   ! vector opt.
(  187)             IF( ABS( glamf(ji,jj) - glamf(ji,jj-1) ) < 1.e-8 ) THEN
(  188)                gsinu(ji,jj) = 0.
(  189)                gcosu(ji,jj) = 1.
(  190)             ENDIF
(  191)             IF( ABS( gphif(ji,jj) - gphif(ji-1,jj) ) < 1.e-8 ) THEN
(  192)                gsinv(ji,jj) = 0.
(  193)                gcosv(ji,jj) = 1.
(  194)             ENDIF
(  195)          END DO
(  196)       END DO
(  197) 
(  198)       ! =========================== !
(  199)       ! Lateral boundary conditions !
(  200)       ! =========================== !
(  201) 
(  202)       ! lateral boundary cond.: U-, V-pts, sgn
(  203)       CALL lbc_lnk ( gsinu, 'U', -1. )   ;   CALL lbc_lnk( gsinv, 'V', -1. )
(  204)       CALL lbc_lnk ( gcosu, 'U', -1. )   ;   CALL lbc_lnk( gcosv, 'V', -1. )
(  205) 
(  206)    END SUBROUTINE angle
(  207) 
(  208) 
(  209)    SUBROUTINE geo2oce ( pxx , pyy , pzz, cgrid,     &
(  210)                         plon, plat, pte, ptn  , ptv )
(  211)       !!----------------------------------------------------------------------
(  212)       !!                    ***  ROUTINE geo2oce  ***
(  213)       !!      
(  214)       !! ** Purpose :
(  215)       !!
(  216)       !! ** Method  :   Change wind stress from geocentric to east/north
(  217)       !!
(  218)       !! History :
(  219)       !!        !         (O. Marti)  Original code
(  220)       !!        !  91-03  (G. Madec)






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 5

(  221)       !!        !  92-07  (M. Imbard)
(  222)       !!        !  99-11  (M. Imbard) NetCDF format with IOIPSL
(  223)       !!        !  00-08  (D. Ludicone) Reduced section at Bab el Mandeb
(  224)       !!   8.5  !  02-06  (G. Madec)  F90: Free form
(  225)       !!----------------------------------------------------------------------
(  226)       !! * Local declarations
(  227)       REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
(  228)          pxx, pyy, pzz
(  229)       CHARACTER (len=1), INTENT( in) ::   &
(  230)          cgrid
(  231)       REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
(  232)          plon, plat
(  233)       REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::    &
(  234)          pte, ptn, ptv
(  235)       REAL(wp), PARAMETER :: rpi = 3.141592653E0
(  236)       REAL(wp), PARAMETER :: rad = rpi / 180.e0
(  237) 
(  238)       !! * Local variables
(  239)       INTEGER ::   ig     !
(  240) 
(  241)       !! * Local save
(  242)       REAL(wp), SAVE, DIMENSION(jpi,jpj,4) ::   &
(  243)          zsinlon, zcoslon,   &
(  244)          zsinlat, zcoslat
(  245)       LOGICAL, SAVE, DIMENSION (4) ::   &
(  246)          linit = .FALSE.
(  247)       !!----------------------------------------------------------------------
(  248) 
(  249)       SELECT CASE( cgrid)
(  250) 
(  251)          CASE ( 't' ) ;; ig = 1
(  252)          CASE ( 'u' ) ;; ig = 2
(  253)          CASE ( 'v' ) ;; ig = 3
(  254)          CASE ( 'f' ) ;; ig = 4
(  255) 
(  256)          CASE default
(  257)             IF(lwp) WRITE(numout,cform_err)
(  258)             IF(lwp) WRITE(numout,*) 'geo2oce : bad grid argument : ', cgrid
(  259)             nstop = nstop + 1
(  260)        END SELECT
(  261)       
(  262)       IF( .NOT. linit(ig) ) THEN 
(  263)          zsinlon (:,:,ig) = SIN (rad * plon)
(  264)          zcoslon (:,:,ig) = COS (rad * plon)
(  265)          zsinlat (:,:,ig) = SIN (rad * plat)
(  266)          zcoslat (:,:,ig) = COS (rad * plat)
(  267)          linit (ig) = .TRUE.
(  268)       ENDIF
(  269)       
(  270)       pte = - zsinlon (:,:,ig) * pxx + zcoslon (:,:,ig) * pyy
(  271)       ptn = - zcoslon (:,:,ig) * zsinlat (:,:,ig) * pxx    &
(  272)             - zsinlon (:,:,ig) * zsinlat (:,:,ig) * pyy    &
(  273)             + zcoslat (:,:,ig) * pzz
(  274)       ptv =   zcoslon (:,:,ig) * zcoslat (:,:,ig) * pxx    &
(  275)             + zsinlon (:,:,ig) * zcoslat (:,:,ig) * pyy    &
(  276)             + zsinlat (:,:,ig) * pzz
(  277) 
(  278)    END SUBROUTINE geo2oce






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 6

(  279) 
(  280) 
(  281)    SUBROUTINE repere ( px1, py1, px2, py2, kchoix )
(  282)       !!----------------------------------------------------------------------
(  283)       !!                 ***  ROUTINE repere  ***
(  284)       !!        
(  285)       !! ** Purpose :   Change vector componantes between a geopgraphic grid 
(  286)       !!      and a stretched coordinates grid.
(  287)       !!
(  288)       !! ** Method  :   initialization of arrays at the first call.
(  289)       !!
(  290)       !! ** Action  :
(  291)       !!
(  292)       !! History :
(  293)       !!        !  89-03  (O. Marti)  original code
(  294)       !!        !  92-02  (M. Imbard)
(  295)       !!        !  93-03  (M. Guyon)  symetrical conditions
(  296)       !!        !  98-05  (B. Blanke)
(  297)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(  298)       !!----------------------------------------------------------------------
(  299)       !! * Arguments
(  300)       REAL(wp), INTENT( in   ), DIMENSION(jpi,jpj) ::   &
(  301)          px1, py1          ! two horizontal components to be rotated
(  302)       REAL(wp), INTENT( out  ), DIMENSION(jpi,jpj) ::   &
(  303)          px2, py2          ! the two horizontal components in the model repere
(  304)       INTEGER, INTENT( inout ) ::   &
(  305)          kchoix   ! type of transformation
(  306)                   ! = 1 change from geographic to model grid.
(  307)                   ! =-1 change from model to geographic grid
(  308)                   ! = 0 same as the previous call
(  309)       !! * Local declarations
(  310)       INTEGER, SAVE :: nmem
(  311) 
(  312)       INTEGER ::   ji, jj                    ! dummy loop indices
(  313) 
(  314)       REAL(wp) :: zxx, zcof1, zcof2,    &
(  315)          ze1t, ze2t
(  316)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  317)          zlamdu, zphiu,   &
(  318)          zlamdv, zphiv
(  319)       !!----------------------------------------------------------------------
(  320) 
(  321) 
(  322)       ! 0. Initialization of gsinus and gcosin IF first call
(  323)       ! ----------------------------------------------------
(  324)       
(  325)       ! 0.1 control argument
(  326)       
(  327)       IF( kchoix == 0 ) THEN
(  328)          IF( nmem == 0 ) THEN
(  329)             IF(lwp) WRITE(numout,cform_err)
(  330)             IF(lwp) WRITE(numout,*) 'repere : e r r o r  in kchoix : ', kchoix
(  331)             IF(lwp) WRITE(numout,*) ' for the first call , you must indicate '
(  332)             IF(lwp) WRITE(numout,*) ' the direction of change '
(  333)             IF(lwp) WRITE(numout,*) ' kchoix = 1 geo       --> stretched '
(  334)             IF(lwp) WRITE(numout,*) ' kchoix =-1 stretched --> geo '
(  335)             nstop = nstop + 1
(  336)          ELSE






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 7

(  337)             kchoix = nmem
(  338)          ENDIF
(  339)       ELSEIF( kchoix == 1 .OR. kchoix == -1 ) THEN
(  340)          nmem = kchoix
(  341)       ELSE
(  342)          IF(lwp) WRITE(numout,cform_err)
(  343)          IF(lwp) WRITE(numout,*) 'repere : e r r o r  in kchoix : ', kchoix
(  344)          IF(lwp) WRITE(numout,*) ' kchoix must be equal to -1, 0 or 1 '
(  345)          nstop = nstop + 1
(  346)       ENDIF
(  347) 
(  348)       ! 0.2 Initialization
(  349) 
(  350)       zxx = gsinus(jpi/2,jpj/2)**2+gcosin(jpi/2,jpj/2)**2
(  351)       IF( zxx <= 0.1 ) THEN
(  352)          IF(lwp) WRITE(numout,*) 'repere : initialization '
(  353)          DO jj = 1, jpj
(  354)             DO ji = 2, jpi
(  355)                zlamdu(ji,jj) = glamu(ji,jj) - glamu(ji-1,jj)
(  356)                zlamdu(ji,jj) = ASIN( SIN( rad*zlamdu(ji,jj) ) )/rad
(  357)                zphiu (ji,jj) = gphiu(ji,jj) - gphiu(ji-1,jj)
(  358)             END DO
(  359)          END DO
(  360)          DO jj = 2, jpj
(  361)             zlamdv(:,jj) = glamv(:,jj)-glamv(:,jj-1)
(  362)             zlamdv(:,jj) = ASIN(SIN(rad*zlamdv(:,jj)))/rad
(  363)             zphiv (:,jj) = gphiv(:,jj)-gphiv(:,jj-1)
(  364)          END DO
(  365)          
(  366)          ! 0.3 Boudary conditions and periodicity
(  367)          
(  368)          IF( nperio == 1 .OR.nperio == 4.OR.nperio == 6 ) THEN
(  369)             DO jj = 1, jpj
(  370)                zlamdu(1,jj) = zlamdu(jpim1,jj)
(  371)                zphiu (1,jj) = zphiu (jpim1,jj)
(  372)             END DO
(  373)          ELSE
(  374)             DO jj = 1, jpj
(  375)                zlamdu(1,jj) = zlamdu(2,jj)
(  376)                zphiu (1,jj) = zphiu (2,jj)
(  377)             END DO
(  378)          ENDIF
(  379)          
(  380)          DO ji = 1, jpi
(  381)             zlamdv(ji,1) = zlamdv(ji,2)
(  382)             zphiv (ji,1) = zphiv (ji,2)
(  383)          END DO
(  384)          
(  385)          IF( nperio == 2 ) THEN
(  386)             DO ji = 1, jpi
(  387)                zphiv (ji,1) = zphiv (ji,3)
(  388)             END DO
(  389)          ENDIF
(  390)          
(  391)          ! 0.4 gsinus gcosin
(  392)          
(  393) !CDIR NOVERRCHK
(  394)          DO jj = 1, jpj






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 8

(  395) !CDIR NOVERRCHK
(  396)             DO ji = 1, jpi
(  397)                zcof1 = rad * ra * COS( rad * gphit(ji,jj) )
(  398)                zcof2 = rad * ra
(  399)                zlamdu(ji,jj) = zlamdu(ji,jj) * zcof1
(  400)                zlamdv(ji,jj) = zlamdv(ji,jj) * zcof1
(  401)                zphiu (ji,jj) = zphiu (ji,jj) * zcof2
(  402)                zphiv (ji,jj) = zphiv (ji,jj) * zcof2
(  403)             END DO
(  404)          END DO
(  405) 
(  406) !CDIR NOVERRCHK
(  407)          DO jj = 1, jpj
(  408) !CDIR NOVERRCHK
(  409)             DO ji = 1, jpi
(  410)                ze1t = SQRT( zlamdu(ji,jj)*zlamdu(ji,jj) + zphiu(ji,jj)*zphiu(ji,jj) )
(  411)                ze2t = SQRT( zlamdv(ji,jj)*zlamdv(ji,jj) + zphiv(ji,jj)*zphiv(ji,jj) )
(  412)                gsinus(ji,jj) = 0.5*( zphiu(ji,jj)/ze1t - zlamdv(ji,jj)/ze2t )
(  413)                gcosin(ji,jj) = 0.5*( zphiv(ji,jj)/ze2t + zlamdu(ji,jj)/ze1t )
(  414)             END DO
(  415)          END DO
(  416)          CALL lbc_lnk( gsinus, 'T', -1. )
(  417)          CALL lbc_lnk( gcosin, 'T', -1. )
(  418)          
(  419)       ENDIF
(  420) 
(  421) 
(  422)       ! 1. Change from geographic to tretched
(  423)       ! -------------------------------------
(  424)       
(  425)       IF( kchoix == 1 ) THEN
(  426)           px2(:,:) =  px1(:,:) * gcosin(:,:) + py1(:,:) * gsinus(:,:)
(  427)           py2(:,:) = -px1(:,:) * gsinus(:,:) + py1(:,:) * gcosin(:,:)
(  428)       ENDIF
(  429)       
(  430) 
(  431)       ! 2. Change from tretched to geographic
(  432)       ! -------------------------------------
(  433)       
(  434)       IF( kchoix == -1 ) THEN
(  435)           px2(:,:) =  px1(:,:) * gcosin(:,:) - py1(:,:) * gsinus(:,:)
(  436)           py2(:,:) =  px1(:,:) * gsinus(:,:) + py1(:,:) * gcosin(:,:)
(  437)       ENDIF
(  438)       
(  439)    END SUBROUTINE repere
(  440) 
(  441)   !!======================================================================
(  442) END MODULE geo2ocean
















PGF90 (Version     10.5)          09/23/2010  09:49:11      page 9

(    1) # 1 "geo2ocean.F90"
(    1) MODULE geo2ocean
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  geo2ocean  ***
(    4)    !! Ocean mesh    :  ???
(    5)    !!=====================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   repcmo      : 
(    9)    !!   angle       :
(   10)    !!   geo2oce     :
(   11)    !!   repere      :   old routine suppress it ???
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE dom_oce         ! mesh and scale factors
(   15)    USE phycst          ! physical constants
(   16)    USE in_out_manager  ! I/O manager
(   17)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   18) 
(   19)    IMPLICIT NONE
(   20) 
(   21)    !! * Accessibility
(   22)    PRIVATE
(   23)    PUBLIC repcmo   ! routine called by ???.F90
(   24)    PUBLIC geo2oce  ! routine called by ???.F90
(   25)    PUBLIC repere   ! routine called by ???.F90
(   26) 
(   27)    !! * Module variables
(   28)    REAL(wp), DIMENSION(jpi,jpj) ::   &
(   29)       gsinu , gcosu ,   &  ! matrix element for change grid u (repcmo.F)
(   30)       gsinv , gcosv ,   &  ! matrix element for change grid v (repcmo.F)
(   31)       gsinus, gcosin      ! matrix element for change grid (repere.F)
(   32) 
(   33)   !! * Substitutions
(   34) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 35 "geo2ocean.F90"
(   35) # 35
(   35)    !!---------------------------------------------------------------------------------
(   36)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/geo2ocean.F90,v 1.4 2005/03/27 18:34:46 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!---------------------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE repcmo ( pxu1, pyu1, pxv1, pyv1,   &
(   44)                        px2 , py2 , kt )






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 10

(   45)       !!----------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE repcmo  ***
(   47)       !!
(   48)       !! ** Purpose :   Change vector componantes from a geographic grid to a
(   49)       !!      stretched coordinates grid.
(   50)       !!
(   51)       !! ** Method  :   Initialization of arrays at the first call.
(   52)       !!
(   53)       !! ** Action  : - px2 : first componante (defined at u point)
(   54)       !!              - py2 : second componante (defined at v point)
(   55)       !!
(   56)       !! History :
(   57)       !!   7.0  !  07-96  (O. Marti)  Original code
(   58)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(   59)       !!----------------------------------------------------------------------
(   60)       !! * Arguments 
(   61)       INTEGER,  INTENT( in ) ::   &
(   62)          kt                ! ocean time-step
(   63)       REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   & 
(   64)          pxu1, pyu1,     & ! geographic vector componantes at u-point
(   65)          pxv1, pyv1        ! geographic vector componantes at v-point
(   66)       REAL(wp), INTENT( out ), DIMENSION(jpi,jpj) ::   & 
(   67)          px2,            & ! i-componante (defined at u-point)
(   68)          py2               ! j-componante (defined at v-point)
(   69)       !!----------------------------------------------------------------------
(   70) 
(   71) 
(   72)       ! Initialization of gsin* and gcos* at first call
(   73)       ! -----------------------------------------------
(   74) 
(   75)       IF( kt <= nit000 + 1 ) THEN
(   76)          IF(lwp) WRITE(numout,*)
(   77)          IF(lwp) WRITE(numout,*) 'repcmo : use the geographic to stretched'
(   78)          IF(lwp) WRITE(numout,*) ' ~~~~~   coordinate transformation'
(   79) 
(   80)          CALL angle       ! initialization of the transformation
(   81)       ENDIF
(   82)       
(   83)       ! Change from geographic to stretched coordinate
(   84)       ! ----------------------------------------------
(   85)       
(   86)       px2(:,:) = pxu1(:,:) * gcosu(:,:) + pyu1(:,:) * gsinu(:,:)
(   87)       py2(:,:) = pyv1(:,:) * gcosv(:,:) - pxv1(:,:) * gsinv(:,:)   
(   88)       
(   89)    END SUBROUTINE repcmo
(   90) 
(   91) 
(   92)    SUBROUTINE angle
(   93)       !!----------------------------------------------------------------------
(   94)       !!                  ***  ROUTINE angle  ***
(   95)       !! 
(   96)       !! ** Purpose :   Compute angles between model grid lines and the 
(   97)       !!      direction of the North
(   98)       !!
(   99)       !! ** Method  :
(  100)       !!
(  101)       !! ** Action  :   Compute (gsinu, gcosu, gsinv, gcosv) arrays: sinus and 
(  102)       !!      cosinus of the angle between the north-south axe and the 






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 11

(  103)       !!      j-direction at u and v-points
(  104)       !!
(  105)       !! History :
(  106)       !!   7.0  !  96-07  (O. Marti)  Original code
(  107)       !!   8.0  !  98-06  (G. Madec)
(  108)       !!   8.5  !  98-06  (G. Madec)  Free form, F90 + opt.
(  109)       !!----------------------------------------------------------------------
(  110)       !! * local declarations
(  111)       INTEGER ::   ji, jj      ! dummy loop indices
(  112) 
(  113)       REAL(wp) ::   &
(  114)          zlam, zphi,             &  ! temporary scalars
(  115)          zlan, zphh,             &  !    "         "
(  116)          zxnpu, zxnpv , znnpu,   &  !    "         "
(  117)          zynpu, zynpv , znnpv,   &  !    "         "
(  118)          zxffu, zmnpfu, zxffv,   &  !    "         "
(  119)          zyffu, zmnpfv, zyffv       !    "         "
(  120)       !!----------------------------------------------------------------------
(  121) 
(  122)       ! ============================= !
(  123)       ! Compute the cosinus and sinus !
(  124)       ! ============================= !
(  125)       ! (computation done on the north stereographic polar plane)
(  126) 
(  127)       DO jj = 2, jpj
(  128) !CDIR NOVERRCHK
(  129)          DO ji = 2, jpi   ! vector opt.
(  130) 
(  131)             ! north pole direction & modulous (at u-point)
(  132)             zlam = glamu(ji,jj)
(  133)             zphi = gphiu(ji,jj)
(  134)             zxnpu = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
(  135)             zynpu = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
(  136)             znnpu = zxnpu*zxnpu + zynpu*zynpu
(  137) 
(  138)             ! north pole direction & modulous (at v-point)
(  139)             zlam = glamv(ji,jj)
(  140)             zphi = gphiv(ji,jj)
(  141)             zxnpv = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
(  142)             zynpv = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
(  143)             znnpv = zxnpv*zxnpv + zynpv*zynpv
(  144) 
(  145)             ! j-direction: f-point segment direction (u-point)
(  146)             zlam = glamf(ji,jj  )
(  147)             zphi = gphif(ji,jj  )
(  148)             zlan = glamf(ji,jj-1)
(  149)             zphh = gphif(ji,jj-1)
(  150)             zxffu =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
(  151)                &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
(  152)             zyffu =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
(  153)                &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
(  154)             zmnpfu = SQRT ( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
(  155)             zmnpfu = MAX( zmnpfu, 1.e-14 )
(  156) 
(  157)             ! i-direction: f-point segment direction (v-point)
(  158)             zlam = glamf(ji  ,jj)
(  159)             zphi = gphif(ji  ,jj)
(  160)             zlan = glamf(ji-1,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 12

(  161)             zphh = gphif(ji-1,jj)
(  162)             zxffv =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
(  163)                &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
(  164)             zyffv =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
(  165)                &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
(  166)             zmnpfv = SQRT ( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
(  167)             zmnpfv = MAX( zmnpfv, 1.e-14 )
(  168) 
(  169)             ! cosinus and sinus using scalar and vectorial products
(  170)             gsinu(ji,jj) = ( zxnpu*zyffu - zynpu*zxffu ) / zmnpfu
(  171)             gcosu(ji,jj) = ( zxnpu*zxffu + zynpu*zyffu ) / zmnpfu
(  172) 
(  173)             ! cosinus and sinus using scalar and vectorial products
(  174)             ! (caution, rotation of 90 degres)
(  175)             gsinv(ji,jj) = ( zxnpv*zxffv + zynpv*zyffv ) / zmnpfv
(  176)             gcosv(ji,jj) =-( zxnpv*zyffv - zynpv*zxffv ) / zmnpfv
(  177) 
(  178)          END DO
(  179)       END DO
(  180) 
(  181)       ! =============== !
(  182)       ! Geographic mesh !
(  183)       ! =============== !
(  184) 
(  185)       DO jj = 2, jpj
(  186)          DO ji = 2, jpi   ! vector opt.
(  187)             IF( ABS( glamf(ji,jj) - glamf(ji,jj-1) ) < 1.e-8 ) THEN
(  188)                gsinu(ji,jj) = 0.
(  189)                gcosu(ji,jj) = 1.
(  190)             ENDIF
(  191)             IF( ABS( gphif(ji,jj) - gphif(ji-1,jj) ) < 1.e-8 ) THEN
(  192)                gsinv(ji,jj) = 0.
(  193)                gcosv(ji,jj) = 1.
(  194)             ENDIF
(  195)          END DO
(  196)       END DO
(  197) 
(  198)       ! =========================== !
(  199)       ! Lateral boundary conditions !
(  200)       ! =========================== !
(  201) 
(  202)       ! lateral boundary cond.: U-, V-pts, sgn
(  203)       CALL lbc_lnk ( gsinu, 'U', -1. )   ;   CALL lbc_lnk( gsinv, 'V', -1. )
(  204)       CALL lbc_lnk ( gcosu, 'U', -1. )   ;   CALL lbc_lnk( gcosv, 'V', -1. )
(  205) 
(  206)    END SUBROUTINE angle
(  207) 
(  208) 
(  209)    SUBROUTINE geo2oce ( pxx , pyy , pzz, cgrid,     &
(  210)                         plon, plat, pte, ptn  , ptv )
(  211)       !!----------------------------------------------------------------------
(  212)       !!                    ***  ROUTINE geo2oce  ***
(  213)       !!      
(  214)       !! ** Purpose :
(  215)       !!
(  216)       !! ** Method  :   Change wind stress from geocentric to east/north
(  217)       !!
(  218)       !! History :






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 13

(  219)       !!        !         (O. Marti)  Original code
(  220)       !!        !  91-03  (G. Madec)
(  221)       !!        !  92-07  (M. Imbard)
(  222)       !!        !  99-11  (M. Imbard) NetCDF format with IOIPSL
(  223)       !!        !  00-08  (D. Ludicone) Reduced section at Bab el Mandeb
(  224)       !!   8.5  !  02-06  (G. Madec)  F90: Free form
(  225)       !!----------------------------------------------------------------------
(  226)       !! * Local declarations
(  227)       REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
(  228)          pxx, pyy, pzz
(  229)       CHARACTER (len=1), INTENT( in) ::   &
(  230)          cgrid
(  231)       REAL(wp), INTENT( in ), DIMENSION(jpi,jpj) ::   &
(  232)          plon, plat
(  233)       REAL(wp), INTENT(out), DIMENSION(jpi,jpj) ::    &
(  234)          pte, ptn, ptv
(  235)       REAL(wp), PARAMETER :: rpi = 3.141592653E0
(  236)       REAL(wp), PARAMETER :: rad = rpi / 180.e0
(  237) 
(  238)       !! * Local variables
(  239)       INTEGER ::   ig     !
(  240) 
(  241)       !! * Local save
(  242)       REAL(wp), SAVE, DIMENSION(jpi,jpj,4) ::   &
(  243)          zsinlon, zcoslon,   &
(  244)          zsinlat, zcoslat
(  245)       LOGICAL, SAVE, DIMENSION (4) ::   &
(  246)          linit = .FALSE.
(  247)       !!----------------------------------------------------------------------
(  248) 
(  249)       SELECT CASE( cgrid)
(  250) 
(  251)          CASE ( 't' ) ;; ig = 1
(  252)          CASE ( 'u' ) ;; ig = 2
(  253)          CASE ( 'v' ) ;; ig = 3
(  254)          CASE ( 'f' ) ;; ig = 4
(  255) 
(  256)          CASE default
(  257)             IF(lwp) WRITE(numout,cform_err)
(  258)             IF(lwp) WRITE(numout,*) 'geo2oce : bad grid argument : ', cgrid
(  259)             nstop = nstop + 1
(  260)        END SELECT
(  261)       
(  262)       IF( .NOT. linit(ig) ) THEN 
(  263)          zsinlon (:,:,ig) = SIN (rad * plon)
(  264)          zcoslon (:,:,ig) = COS (rad * plon)
(  265)          zsinlat (:,:,ig) = SIN (rad * plat)
(  266)          zcoslat (:,:,ig) = COS (rad * plat)
(  267)          linit (ig) = .TRUE.
(  268)       ENDIF
(  269)       
(  270)       pte = - zsinlon (:,:,ig) * pxx + zcoslon (:,:,ig) * pyy
(  271)       ptn = - zcoslon (:,:,ig) * zsinlat (:,:,ig) * pxx    &
(  272)             - zsinlon (:,:,ig) * zsinlat (:,:,ig) * pyy    &
(  273)             + zcoslat (:,:,ig) * pzz
(  274)       ptv =   zcoslon (:,:,ig) * zcoslat (:,:,ig) * pxx    &
(  275)             + zsinlon (:,:,ig) * zcoslat (:,:,ig) * pyy    &
(  276)             + zsinlat (:,:,ig) * pzz






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 14

(  277) 
(  278)    END SUBROUTINE geo2oce
(  279) 
(  280) 
(  281)    SUBROUTINE repere ( px1, py1, px2, py2, kchoix )
(  282)       !!----------------------------------------------------------------------
(  283)       !!                 ***  ROUTINE repere  ***
(  284)       !!        
(  285)       !! ** Purpose :   Change vector componantes between a geopgraphic grid 
(  286)       !!      and a stretched coordinates grid.
(  287)       !!
(  288)       !! ** Method  :   initialization of arrays at the first call.
(  289)       !!
(  290)       !! ** Action  :
(  291)       !!
(  292)       !! History :
(  293)       !!        !  89-03  (O. Marti)  original code
(  294)       !!        !  92-02  (M. Imbard)
(  295)       !!        !  93-03  (M. Guyon)  symetrical conditions
(  296)       !!        !  98-05  (B. Blanke)
(  297)       !!   8.5  !  02-08  (G. Madec)  F90: Free form
(  298)       !!----------------------------------------------------------------------
(  299)       !! * Arguments
(  300)       REAL(wp), INTENT( in   ), DIMENSION(jpi,jpj) ::   &
(  301)          px1, py1          ! two horizontal components to be rotated
(  302)       REAL(wp), INTENT( out  ), DIMENSION(jpi,jpj) ::   &
(  303)          px2, py2          ! the two horizontal components in the model repere
(  304)       INTEGER, INTENT( inout ) ::   &
(  305)          kchoix   ! type of transformation
(  306)                   ! = 1 change from geographic to model grid.
(  307)                   ! =-1 change from model to geographic grid
(  308)                   ! = 0 same as the previous call
(  309)       !! * Local declarations
(  310)       INTEGER, SAVE :: nmem
(  311) 
(  312)       INTEGER ::   ji, jj                    ! dummy loop indices
(  313) 
(  314)       REAL(wp) :: zxx, zcof1, zcof2,    &
(  315)          ze1t, ze2t
(  316)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(  317)          zlamdu, zphiu,   &
(  318)          zlamdv, zphiv
(  319)       !!----------------------------------------------------------------------
(  320) 
(  321) 
(  322)       ! 0. Initialization of gsinus and gcosin IF first call
(  323)       ! ----------------------------------------------------
(  324)       
(  325)       ! 0.1 control argument
(  326)       
(  327)       IF( kchoix == 0 ) THEN
(  328)          IF( nmem == 0 ) THEN
(  329)             IF(lwp) WRITE(numout,cform_err)
(  330)             IF(lwp) WRITE(numout,*) 'repere : e r r o r  in kchoix : ', kchoix
(  331)             IF(lwp) WRITE(numout,*) ' for the first call , you must indicate '
(  332)             IF(lwp) WRITE(numout,*) ' the direction of change '
(  333)             IF(lwp) WRITE(numout,*) ' kchoix = 1 geo       --> stretched '
(  334)             IF(lwp) WRITE(numout,*) ' kchoix =-1 stretched --> geo '






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 15

(  335)             nstop = nstop + 1
(  336)          ELSE
(  337)             kchoix = nmem
(  338)          ENDIF
(  339)       ELSEIF( kchoix == 1 .OR. kchoix == -1 ) THEN
(  340)          nmem = kchoix
(  341)       ELSE
(  342)          IF(lwp) WRITE(numout,cform_err)
(  343)          IF(lwp) WRITE(numout,*) 'repere : e r r o r  in kchoix : ', kchoix
(  344)          IF(lwp) WRITE(numout,*) ' kchoix must be equal to -1, 0 or 1 '
(  345)          nstop = nstop + 1
(  346)       ENDIF
(  347) 
(  348)       ! 0.2 Initialization
(  349) 
(  350)       zxx = gsinus(jpi/2,jpj/2)**2+gcosin(jpi/2,jpj/2)**2
(  351)       IF( zxx <= 0.1 ) THEN
(  352)          IF(lwp) WRITE(numout,*) 'repere : initialization '
(  353)          DO jj = 1, jpj
(  354)             DO ji = 2, jpi
(  355)                zlamdu(ji,jj) = glamu(ji,jj) - glamu(ji-1,jj)
(  356)                zlamdu(ji,jj) = ASIN( SIN( rad*zlamdu(ji,jj) ) )/rad
(  357)                zphiu (ji,jj) = gphiu(ji,jj) - gphiu(ji-1,jj)
(  358)             END DO
(  359)          END DO
(  360)          DO jj = 2, jpj
(  361)             zlamdv(:,jj) = glamv(:,jj)-glamv(:,jj-1)
(  362)             zlamdv(:,jj) = ASIN(SIN(rad*zlamdv(:,jj)))/rad
(  363)             zphiv (:,jj) = gphiv(:,jj)-gphiv(:,jj-1)
(  364)          END DO
(  365)          
(  366)          ! 0.3 Boudary conditions and periodicity
(  367)          
(  368)          IF( nperio == 1 .OR.nperio == 4.OR.nperio == 6 ) THEN
(  369)             DO jj = 1, jpj
(  370)                zlamdu(1,jj) = zlamdu(jpim1,jj)
(  371)                zphiu (1,jj) = zphiu (jpim1,jj)
(  372)             END DO
(  373)          ELSE
(  374)             DO jj = 1, jpj
(  375)                zlamdu(1,jj) = zlamdu(2,jj)
(  376)                zphiu (1,jj) = zphiu (2,jj)
(  377)             END DO
(  378)          ENDIF
(  379)          
(  380)          DO ji = 1, jpi
(  381)             zlamdv(ji,1) = zlamdv(ji,2)
(  382)             zphiv (ji,1) = zphiv (ji,2)
(  383)          END DO
(  384)          
(  385)          IF( nperio == 2 ) THEN
(  386)             DO ji = 1, jpi
(  387)                zphiv (ji,1) = zphiv (ji,3)
(  388)             END DO
(  389)          ENDIF
(  390)          
(  391)          ! 0.4 gsinus gcosin
(  392)          






PGF90 (Version     10.5)          09/23/2010  09:49:11      page 16

(  393) !CDIR NOVERRCHK
(  394)          DO jj = 1, jpj
(  395) !CDIR NOVERRCHK
(  396)             DO ji = 1, jpi
(  397)                zcof1 = rad * ra * COS( rad * gphit(ji,jj) )
(  398)                zcof2 = rad * ra
(  399)                zlamdu(ji,jj) = zlamdu(ji,jj) * zcof1
(  400)                zlamdv(ji,jj) = zlamdv(ji,jj) * zcof1
(  401)                zphiu (ji,jj) = zphiu (ji,jj) * zcof2
(  402)                zphiv (ji,jj) = zphiv (ji,jj) * zcof2
(  403)             END DO
(  404)          END DO
(  405) 
(  406) !CDIR NOVERRCHK
(  407)          DO jj = 1, jpj
(  408) !CDIR NOVERRCHK
(  409)             DO ji = 1, jpi
(  410)                ze1t = SQRT( zlamdu(ji,jj)*zlamdu(ji,jj) + zphiu(ji,jj)*zphiu(ji,jj) )
(  411)                ze2t = SQRT( zlamdv(ji,jj)*zlamdv(ji,jj) + zphiv(ji,jj)*zphiv(ji,jj) )
(  412)                gsinus(ji,jj) = 0.5*( zphiu(ji,jj)/ze1t - zlamdv(ji,jj)/ze2t )
(  413)                gcosin(ji,jj) = 0.5*( zphiv(ji,jj)/ze2t + zlamdu(ji,jj)/ze1t )
(  414)             END DO
(  415)          END DO
(  416)          CALL lbc_lnk( gsinus, 'T', -1. )
(  417)          CALL lbc_lnk( gcosin, 'T', -1. )
(  418)          
(  419)       ENDIF
(  420) 
(  421) 
(  422)       ! 1. Change from geographic to tretched
(  423)       ! -------------------------------------
(  424)       
(  425)       IF( kchoix == 1 ) THEN
(  426)           px2(:,:) =  px1(:,:) * gcosin(:,:) + py1(:,:) * gsinus(:,:)
(  427)           py2(:,:) = -px1(:,:) * gsinus(:,:) + py1(:,:) * gcosin(:,:)
(  428)       ENDIF
(  429)       
(  430) 
(  431)       ! 2. Change from tretched to geographic
(  432)       ! -------------------------------------
(  433)       
(  434)       IF( kchoix == -1 ) THEN
(  435)           px2(:,:) =  px1(:,:) * gcosin(:,:) - py1(:,:) * gsinus(:,:)
(  436)           py2(:,:) =  px1(:,:) * gsinus(:,:) + py1(:,:) * gcosin(:,:)
(  437)       ENDIF
(  438)       
(  439)    END SUBROUTINE repere
(  440) 
(  441)   !!======================================================================
(  442) END MODULE geo2ocean
