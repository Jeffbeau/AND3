


PGF90 (Version     10.5)          09/23/2010  09:49:54      page 1

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

Filename: limadv.F90

(    1) MODULE limadv 
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limadv   ***
(    4)    !! LIM sea-ice model : sea-ice advection
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_adv_x  : advection of sea ice on x axis
(   11)    !!   lim_adv_y  : advection of sea ice on y axis
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE dom_oce
(   15)    USE dom_ice
(   16)    USE ice_oce         ! ice variables
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE lbclnk
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC lim_adv_x    ! called by lim_trp
(   26)    PUBLIC lim_adv_y    ! called by lim_trp
(   27) 
(   28)    !! * Module variables
(   29)    REAL(wp)  ::            &  ! constant values
(   30)       epsi20 = 1e-20    ,  &
(   31)       rzero  = 0.e0     ,  &
(   32)       rone   = 1.e0
(   33)    !!----------------------------------------------------------------------
(   34)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limadv.F90,v 1.5 2005/09/02 15:45:16 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE lim_adv_x( pdf, put , pcrh, psm , ps0 ,   &
(   42)       &                  psx, psxx, psy , psyy, psxy )
(   43)       !!---------------------------------------------------------------------
(   44)       !!                **  routine lim_adv_x  **
(   45)       !!  
(   46)       !! ** purpose :   Computes and adds the advection trend to sea-ice






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 2

(   47)       !!      variable on x axis
(   48)       !!
(   49)       !! ** method  :   Uses Prather second order scheme that advects
(   50)       !!      tracers but also theirquadratic forms. The method preserves
(   51)       !!      tracer structures by conserving second order moments.
(   52)       !!      Ref.: "Numerical Advection by Conservation of Second Order
(   53)       !!      Moments", JGR, VOL. 91. NO. D6. PAGES 6671-6681. MAY 20, 1986
(   54)       !!     
(   55)       !! History :
(   56)       !!        !  00-01 (LIM)
(   57)       !!        !  01-05 (G. Madec, R. Hordoir) opa norm
(   58)       !!        !  03-10 (C. Ethe) F90, module
(   59)       !!        !  03-12 (R. Hordoir, G. Madec) mpp
(   60)       !!--------------------------------------------------------------------
(   61)       !! * Arguments
(   62)       REAL(wp)                    , INTENT(in)  ::  &
(   63)          pdf ,       &  ! ???
(   64)          pcrh           ! = 1. : lim_adv_x is called before lim_adv_y
(   65)          !              ! = 0. : lim_adv_x is called after  lim_adv_y
(   66)       REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::  &
(   67)          put            ! i-direction ice velocity at ocean U-point (m/s)
(   68)       REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::  & 
(   69)          ps0 , psm , &  ! ???
(   70)          psx , psy , &  ! ??? 
(   71)          psxx, psyy, psxy
(   72) 
(   73)       !! * Local declarations
(   74)       INTEGER ::   ji, jj      ! dummy loop indices
(   75)       REAL(wp)       ::  &
(   76)          zrdt, zslpmax, ztemp, zin0,     &  ! temporary scalars
(   77)          zs1max, zs1new, zs2new,         &  !    "         "
(   78)          zalf, zalfq, zalf1, zalf1q,     &  !    "         "
(   79)          zbt , zbt1                         !    "         "
(   80)       REAL(wp), DIMENSION(jpi,jpj)  ::   &  ! temporary workspace
(   81)          zf0 , zfx , zfy , zbet,         &  !    "           "
(   82)          zfxx, zfyy, zfxy,               &  !    "           "
(   83)          zfm, zalg, zalg1, zalg1q           !    "           "
(   84)       !---------------------------------------------------------------------
(   85) 
(   86)       ! Limitation of moments.                                           
(   87) 
(   88)       zrdt      = rdt_ice * pdf   ! If ice drift field is too fast, use an appropriate time step for advection.
(   89) 
(   90)       DO jj = 1, jpj
(   91)          DO ji = 1, jpi
(   92)             zslpmax = MAX( rzero, ps0(ji,jj) )
(   93)             zs1max  = 1.5 * zslpmax
(   94)             zs1new  = MIN( zs1max, MAX( -zs1max, psx(ji,jj) ) )
(   95)             zs2new  = MIN(  2.0 * zslpmax - 0.3334 * ABS( zs1new ),      &
(   96)                &            MAX( ABS( zs1new ) - zslpmax, psxx(ji,jj) )  )
(   97)             zin0    = ( 1.0 - MAX( rzero, sign ( rone, -zslpmax) ) ) * tms(ji,jj)   ! Case of empty boxes & Apply mask
(   98) 
(   99)             ps0 (ji,jj) = zslpmax  
(  100)             psx (ji,jj) = zs1new      * zin0
(  101)             psxx(ji,jj) = zs2new      * zin0
(  102)             psy (ji,jj) = psy (ji,jj) * zin0
(  103)             psyy(ji,jj) = psyy(ji,jj) * zin0
(  104)             psxy(ji,jj) = MIN( zslpmax, MAX( -zslpmax, psxy(ji,jj) ) ) * zin0






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 3

(  105)          END DO
(  106)       END DO
(  107) 
(  108)       !  Initialize volumes of boxes  (=area if adv_x first called, =psm otherwise)                                     
(  109)       psm (:,:)  = MAX( pcrh * area(:,:) + ( 1.0 - pcrh ) * psm(:,:) , epsi20 )
(  110) 
(  111)       !  Calculate fluxes and moments between boxes i<-->i+1              
(  112)       DO jj = 2, jpjm1                    !  Flux from i to i+1 WHEN u GT 0 
(  113) !i bug   DO ji = 1, jpim1 
(  114) !i    DO jj = 1, jpj                      !  Flux from i to i+1 WHEN u GT 0 
(  115)          DO ji = 1, jpi
(  116)             zbet(ji,jj)  =  MAX( rzero, SIGN( rone, put(ji,jj) ) )
(  117)             zalf         =  MAX( rzero, put(ji,jj) ) * zrdt * e2u(ji,jj) / psm(ji,jj)
(  118)             zalfq        =  zalf * zalf
(  119)             zalf1        =  1.0 - zalf
(  120)             zalf1q       =  zalf1 * zalf1
(  121)             zfm (ji,jj)  =  zalf  * psm(ji,jj)
(  122)             zf0 (ji,jj)  =  zalf  * ( ps0(ji,jj) + zalf1 * ( psx(ji,jj) + (zalf1 - zalf) * psxx(ji,jj)  ) )
(  123)             zfx (ji,jj)  =  zalfq * ( psx(ji,jj) + 3.0 * zalf1 * psxx(ji,jj) )
(  124)             zfxx(ji,jj)  =  zalf  * zalfq * psxx(ji,jj)
(  125)             zfy (ji,jj)  =  zalf  * ( psy(ji,jj) + zalf1 * psxy(ji,jj) )
(  126)             zfxy(ji,jj)  =  zalfq * psxy(ji,jj)
(  127)             zfyy(ji,jj)  =  zalf  * psyy(ji,jj)
(  128) 
(  129)             !  Readjust moments remaining in the box.
(  130)             psm (ji,jj)  =  psm (ji,jj) - zfm(ji,jj)
(  131)             ps0 (ji,jj)  =  ps0 (ji,jj) - zf0(ji,jj)
(  132)             psx (ji,jj)  =  zalf1q * ( psx(ji,jj) - 3.0 * zalf * psxx(ji,jj) )
(  133)             psxx(ji,jj)  =  zalf1  * zalf1q * psxx(ji,jj)
(  134)             psy (ji,jj)  =  psy (ji,jj) - zfy(ji,jj)
(  135)             psyy(ji,jj)  =  psyy(ji,jj) - zfyy(ji,jj)
(  136)             psxy(ji,jj)  =  zalf1q * psxy(ji,jj)
(  137)          END DO
(  138)       END DO
(  139) 
(  140)       DO jj = 2, jpjm1                      !  Flux from i+1 to i when u LT 0.
(  141) !i    DO jj = 1, jpjm1                      !  Flux from i+1 to i when u LT 0.
(  142)          DO ji = 1, jpim1
(  143)             zalf          = MAX( rzero, -put(ji,jj) ) * zrdt * e2u(ji,jj) / psm(ji+1,jj) 
(  144)             zalg  (ji,jj) = zalf
(  145)             zalfq         = zalf * zalf
(  146)             zalf1         = 1.0 - zalf
(  147)             zalg1 (ji,jj) = zalf1
(  148)             zalf1q        = zalf1 * zalf1
(  149)             zalg1q(ji,jj) = zalf1q
(  150)             zfm   (ji,jj) = zfm (ji,jj) + zalf  * psm(ji+1,jj)
(  151)             zf0   (ji,jj) = zf0 (ji,jj) + zalf  * ( ps0(ji+1,jj) - zalf1 * ( psx(ji+1,jj) - (zalf1 - zalf ) * psxx(ji+1,jj) 
(  152)             zfx   (ji,jj) = zfx (ji,jj) + zalfq * ( psx(ji+1,jj) - 3.0 * zalf1 * psxx(ji+1,jj) )
(  153)             zfxx  (ji,jj) = zfxx(ji,jj) + zalf  * zalfq * psxx(ji+1,jj)
(  154)             zfy   (ji,jj) = zfy (ji,jj) + zalf  * ( psy(ji+1,jj) - zalf1 * psxy(ji+1,jj) )
(  155)             zfxy  (ji,jj) = zfxy(ji,jj) + zalfq * psxy(ji+1,jj)
(  156)             zfyy  (ji,jj) = zfyy(ji,jj) + zalf  * psyy(ji+1,jj)
(  157)          END DO
(  158)       END DO
(  159) 
(  160)       DO jj = 2, jpjm1                     !  Readjust moments remaining in the box. 
(  161)          DO ji = 2, jpim1
(  162)             zbt  =       zbet(ji-1,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 4

(  163)             zbt1 = 1.0 - zbet(ji-1,jj)
(  164)             psm (ji,jj) = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) - zfm(ji-1,jj) )
(  165)             ps0 (ji,jj) = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) - zf0(ji-1,jj) )
(  166)             psx (ji,jj) = zalg1q(ji-1,jj) * ( psx(ji,jj) + 3.0 * zalg(ji-1,jj) * psxx(ji,jj) )
(  167)             psxx(ji,jj) = zalg1 (ji-1,jj) * zalg1q(ji-1,jj) * psxx(ji,jj)
(  168)             psy (ji,jj) = zbt * psy (ji,jj) + zbt1 * ( psy (ji,jj) - zfy (ji-1,jj) )
(  169)             psyy(ji,jj) = zbt * psyy(ji,jj) + zbt1 * ( psyy(ji,jj) - zfyy(ji-1,jj) )
(  170)             psxy(ji,jj) = zalg1q(ji-1,jj) * psxy(ji,jj)
(  171)          END DO
(  172)       END DO
(  173) 
(  174)       !   Put the temporary moments into appropriate neighboring boxes.    
(  175)       DO jj = 2, jpjm1                     !   Flux from i to i+1 IF u GT 0.
(  176)          DO ji = 2, jpim1
(  177)             zbt  =       zbet(ji-1,jj)
(  178)             zbt1 = 1.0 - zbet(ji-1,jj)
(  179)             psm(ji,jj)  = zbt * ( psm(ji,jj) + zfm(ji-1,jj) ) + zbt1 * psm(ji,jj)
(  180)             zalf        = zbt * zfm(ji-1,jj) / psm(ji,jj)
(  181)             zalf1       = 1.0 - zalf
(  182)             ztemp       = zalf * ps0(ji,jj) - zalf1 * zf0(ji-1,jj)
(  183)             ps0(ji,jj)  = zbt * (ps0(ji,jj) + zf0(ji-1,jj)) + zbt1 * ps0(ji,jj)
(  184)             psx(ji,jj)  = zbt * ( zalf * zfx(ji-1,jj) + zalf1 * psx(ji,jj) + 3.0 * ztemp ) + zbt1 * psx(ji,jj)
(  185)             psxx(ji,jj) = zbt * ( zalf * zalf * zfxx(ji-1,jj) + zalf1 * zalf1 * psxx(ji,jj)   &
(  186)                &                + 5.0 * ( zalf * zalf1 * ( psx (ji,jj) - zfx(ji-1,jj) ) - ( zalf1 - zalf ) * ztemp )  )   &
(  187)                &        + zbt1 * psxx(ji,jj)
(  188)             psxy(ji,jj) = zbt * ( zalf * zfxy(ji-1,jj) + zalf1 * psxy(ji,jj)             &
(  189)                &                + 3.0 * (- zalf1*zfy(ji-1,jj)  + zalf * psy(ji,jj) ) )   &
(  190)                &         + zbt1 * psxy(ji,jj)
(  191)             psy (ji,jj) = zbt * ( psy (ji,jj) + zfy (ji-1,jj) ) + zbt1 * psy (ji,jj)
(  192)             psyy(ji,jj) = zbt * ( psyy(ji,jj) + zfyy(ji-1,jj) ) + zbt1 * psyy(ji,jj)
(  193)          END DO
(  194)       END DO
(  195) 
(  196)       DO jj = 2, jpjm1                     !  Flux from i+1 to i IF u LT 0.
(  197)          DO ji = 2, jpim1
(  198)             zbt  =       zbet(ji,jj)
(  199)             zbt1 = 1.0 - zbet(ji,jj)
(  200)             psm(ji,jj)  = zbt * psm(ji,jj)  + zbt1 * ( psm(ji,jj) + zfm(ji,jj) )
(  201)             zalf        = zbt1 * zfm(ji,jj) / psm(ji,jj)
(  202)             zalf1       = 1.0 - zalf
(  203)             ztemp       = -zalf * ps0(ji,jj) + zalf1 * zf0(ji,jj)
(  204)             ps0(ji,jj)  = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) + zf0(ji,jj) )
(  205)             psx(ji,jj)  = zbt  * psx(ji,jj)   &
(  206)                &        + zbt1 * ( zalf * zfx(ji,jj) + zalf1 * psx(ji,jj) + 3.0 * ztemp )
(  207)             psxx(ji,jj) = zbt  * psxx(ji,jj)   &
(  208)                &        + zbt1 * (  zalf * zalf * zfxx(ji,jj)  + zalf1 * zalf1 * psxx(ji,jj)  &
(  209)                &                 + 5.0 *( zalf * zalf1 * ( - psx(ji,jj) + zfx(ji,jj) ) + ( zalf1 - zalf ) * ztemp ) )
(  210)             psxy(ji,jj) = zbt  * psxy(ji,jj)   &
(  211)                &        + zbt1 * (  zalf * zfxy(ji,jj) + zalf1 * psxy(ji,jj)  &
(  212)                &                 + 3.0 * ( zalf1 * zfy(ji,jj) - zalf * psy(ji,jj) )  )
(  213)             psy(ji,jj)  = zbt * psy (ji,jj)  + zbt1 * ( psy (ji,jj) + zfy (ji,jj) )
(  214)             psyy(ji,jj) = zbt * psyy(ji,jj)  + zbt1 * ( psyy(ji,jj) + zfyy(ji,jj) )
(  215)          END DO
(  216)       END DO
(  217) 
(  218)       !-- Lateral boundary conditions
(  219)       CALL lbc_lnk( psm , 'T', 1. )
(  220)       CALL lbc_lnk( ps0 , 'T', 1. )






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 5

(  221)       CALL lbc_lnk( psx , 'T', 1. )
(  222)       CALL lbc_lnk( psxx, 'T', 1. )
(  223)       CALL lbc_lnk( psy , 'T', 1. )
(  224)       CALL lbc_lnk( psyy, 'T', 1. )
(  225)       CALL lbc_lnk( psxy, 'T', 1. )
(  226) 
(  227)       IF(ln_ctl)   THEN
(  228)          CALL prt_ctl(tab2d_1=psm  , clinfo1=' lim_adv_x: psm  :', tab2d_2=ps0 , clinfo2=' ps0  : ')
(  229)          CALL prt_ctl(tab2d_1=psx  , clinfo1=' lim_adv_x: psx  :', tab2d_2=psxx, clinfo2=' psxx : ')
(  230)          CALL prt_ctl(tab2d_1=psy  , clinfo1=' lim_adv_x: psy  :', tab2d_2=psyy, clinfo2=' psyy : ')
(  231)          CALL prt_ctl(tab2d_1=psxy , clinfo1=' lim_adv_x: psxy :')
(  232)       ENDIF
(  233) 
(  234)    END SUBROUTINE lim_adv_x
(  235) 
(  236) 
(  237)    SUBROUTINE lim_adv_y( pdf, pvt , pcrh, psm , ps0 ,   &
(  238)       &                  psx, psxx, psy , psyy, psxy )
(  239)       !!---------------------------------------------------------------------
(  240)       !!                **  routine lim_adv_y  **
(  241)       !!            
(  242)       !! ** purpose :   Computes and adds the advection trend to sea-ice 
(  243)       !!      variable on y axis
(  244)       !!
(  245)       !! ** method  :   Uses Prather second order scheme that advects tracers
(  246)       !!      but also their quadratic forms. The method preserves tracer
(  247)       !!      structures by conserving second order moments.
(  248)       !!
(  249)       !! History :
(  250)       !!   1.0  !  00-01 (LIM)
(  251)       !!        !  01-05 (G. Madec, R. Hordoir) opa norm
(  252)       !!   2.0  !  03-10 (C. Ethe) F90, module
(  253)       !!        !  03-12 (R. Hordoir, G. Madec) mpp
(  254)       !!---------------------------------------------------------------------
(  255)       !! * Arguments
(  256)       REAL(wp),                     INTENT(in)  :: &
(  257)          pdf,        &  ! ???
(  258)          pcrh           ! = 1. : lim_adv_x is called before lim_adv_y
(  259)          !              ! = 0. : lim_adv_x is called after  lim_adv_y
(  260)       REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: &
(  261)          pvt            ! j-direction ice velocity at ocean V-point (m/s)
(  262)       REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: &
(  263)          psm , ps0 , psx , psy,   &
(  264)          psxx, psyy, psxy
(  265) 
(  266)       !! * Local Variables
(  267)       INTEGER  ::   ji, jj                ! dummy loop indices
(  268)       REAL(wp) ::   &
(  269)          zrdt, zslpmax, zin0, ztemp,  &   ! temporary scalars
(  270)          zs1max, zs1new, zs2new,      &   !    "         "
(  271)          zalf, zalfq, zalf1, zalf1q,  &   !    "         "
(  272)          zbt , zbt1                       !
(  273)       REAL(wp), DIMENSION(jpi,jpj)  :: &
(  274)          zf0 , zfx , zfy ,             &  ! temporary workspace
(  275)          zfxx, zfyy, zfxy,             &  !    "           "
(  276)          zfm , zbet,                   &  !    "           "
(  277)          zalg, zalg1, zalg1q              !    "           "
(  278)       !---------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 6

(  279) 
(  280)       ! Limitation of moments.
(  281) 
(  282)       zrdt = rdt_ice * pdf ! If ice drift field is too fast, use an appropriate time step for advection.
(  283) 
(  284)        DO jj = 1, jpj
(  285)           DO ji = 1, jpi
(  286)              zslpmax = MAX( rzero, ps0(ji,jj) )
(  287)              zs1max  = 1.5 * zslpmax
(  288)              zs1new  = MIN( zs1max, MAX( -zs1max, psy(ji,jj) ) )
(  289)              zs2new  = MIN(  ( 2.0 * zslpmax - 0.3334 * ABS( zs1new ) ),   &
(  290)                 &             MAX( ABS( zs1new )-zslpmax, psyy(ji,jj) )  )
(  291)              zin0    = ( 1.0 - MAX( rzero, sign ( rone, -zslpmax) ) ) * tms(ji,jj)   ! Case of empty boxes & Apply mask
(  292)              ps0 (ji,jj) = zslpmax  
(  293)              psx (ji,jj)  = psx (ji,jj) * zin0
(  294)              psxx(ji,jj)  = psxx(ji,jj) * zin0
(  295)              psy (ji,jj) = zs1new * zin0
(  296)              psyy(ji,jj) = zs2new * zin0
(  297)              psxy(ji,jj) = MIN( zslpmax, MAX( -zslpmax, psxy(ji,jj) ) ) * zin0
(  298)           END DO
(  299)        END DO
(  300) 
(  301)        !  Initialize volumes of boxes (=area if adv_x first called, =psm otherwise)
(  302)        psm (:,:)  = MAX( pcrh * area(:,:) + ( 1.0 - pcrh ) * psm(:,:) , epsi20 )
(  303) 
(  304)        !  Calculate fluxes and moments between boxes j<-->j+1              
(  305) !!bug  DO jj = 2, jpjm1
(  306)        DO jj = 1, jpj
(  307)           DO ji = 2, jpim1
(  308) !!bug     DO ji = 1, jpim1
(  309)              !  Flux from j to j+1 WHEN v GT 0   
(  310)              zbet(ji,jj)  =  MAX( rzero, SIGN( rone, pvt(ji,jj) ) )
(  311)              zalf         =  MAX( rzero, pvt(ji,jj) ) * zrdt * e1v(ji,jj) / psm(ji,jj)
(  312)              zalfq        =  zalf * zalf
(  313)              zalf1        =  1.0 - zalf
(  314)              zalf1q       =  zalf1 * zalf1
(  315)              zfm (ji,jj)  =  zalf  * psm(ji,jj)
(  316)              zf0 (ji,jj)  =  zalf  * ( ps0(ji,jj) + zalf1 * ( psy(ji,jj)  + (zalf1-zalf) * psyy(ji,jj)  ) ) 
(  317)              zfy (ji,jj)  =  zalfq *( psy(ji,jj) + 3.0*zalf1*psyy(ji,jj) )
(  318)              zfyy(ji,jj)  =  zalf  * zalfq * psyy(ji,jj)
(  319)              zfx (ji,jj)  =  zalf  * ( psx(ji,jj) + zalf1 * psxy(ji,jj) )
(  320)              zfxy(ji,jj)  =  zalfq * psxy(ji,jj)
(  321)              zfxx(ji,jj)  =  zalf  * psxx(ji,jj)
(  322) 
(  323)              !  Readjust moments remaining in the box.
(  324)              psm (ji,jj)  =  psm (ji,jj) - zfm(ji,jj)
(  325)              ps0 (ji,jj)  =  ps0 (ji,jj) - zf0(ji,jj)
(  326)              psy (ji,jj)  =  zalf1q * ( psy(ji,jj) -3.0 * zalf * psyy(ji,jj) )
(  327)              psyy(ji,jj)  =  zalf1 * zalf1q * psyy(ji,jj)
(  328)              psx (ji,jj)  =  psx (ji,jj) - zfx(ji,jj)
(  329)              psxx(ji,jj)  =  psxx(ji,jj) - zfxx(ji,jj)
(  330)              psxy(ji,jj)  =  zalf1q * psxy(ji,jj)
(  331)           END DO
(  332)        END DO
(  333) 
(  334)        DO jj = 1, jpjm1                   !  Flux from j+1 to j when v LT 0.
(  335)           DO ji = 2, jpim1
(  336) !i     DO jj = 1, jpjm1                   !  Flux from j+1 to j when v LT 0.






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 7

(  337) !i        DO ji = 2, jpim1
(  338)              zalf          = ( MAX(rzero, -pvt(ji,jj) ) * zrdt * e1v(ji,jj) ) / psm(ji,jj+1) 
(  339)              zalg  (ji,jj) = zalf
(  340)              zalfq         = zalf * zalf
(  341)              zalf1         = 1.0 - zalf
(  342)              zalg1 (ji,jj) = zalf1
(  343)              zalf1q        = zalf1 * zalf1
(  344)              zalg1q(ji,jj) = zalf1q
(  345)              zfm   (ji,jj) = zfm (ji,jj) + zalf  * psm(ji,jj+1)
(  346)              zf0   (ji,jj) = zf0 (ji,jj) + zalf  * ( ps0(ji,jj+1) - zalf1 * (psy(ji,jj+1) - (zalf1 - zalf ) * psyy(ji,jj+1) 
(  347)              zfy   (ji,jj) = zfy (ji,jj) + zalfq * ( psy(ji,jj+1) - 3.0 * zalf1 * psyy(ji,jj+1) )
(  348)              zfyy  (ji,jj) = zfyy(ji,jj) + zalf  * zalfq * psyy(ji,jj+1)
(  349)              zfx   (ji,jj) = zfx (ji,jj) + zalf  * ( psx(ji,jj+1) - zalf1 * psxy(ji,jj+1) )
(  350)              zfxy  (ji,jj) = zfxy(ji,jj) + zalfq * psxy(ji,jj+1)
(  351)              zfxx  (ji,jj) = zfxx(ji,jj) + zalf  * psxx(ji,jj+1)
(  352)           END DO
(  353)        END DO
(  354)  
(  355)        !  Readjust moments remaining in the box. 
(  356)        DO jj = 2, jpjm1
(  357)           DO ji = 2, jpim1
(  358)              zbt  =         zbet(ji,jj-1)
(  359)              zbt1 = ( 1.0 - zbet(ji,jj-1) )
(  360)              psm (ji,jj) = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) - zfm(ji,jj-1) )
(  361)              ps0 (ji,jj) = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) - zf0(ji,jj-1) )
(  362)              psy (ji,jj) = zalg1q(ji,jj-1) * ( psy(ji,jj) + 3.0 * zalg(ji,jj-1) * psyy(ji,jj) )
(  363)              psyy(ji,jj) = zalg1 (ji,jj-1)  * zalg1q(ji,jj-1) * psyy(ji,jj)
(  364)              psx (ji,jj) = zbt * psx (ji,jj) + zbt1 * ( psx (ji,jj) - zfx (ji,jj-1) )
(  365)              psxx(ji,jj) = zbt * psxx(ji,jj) + zbt1 * ( psxx(ji,jj) - zfxx(ji,jj-1) )
(  366)              psxy(ji,jj) = zalg1q(ji,jj-1) * psxy(ji,jj)
(  367)           END DO
(  368)        END DO
(  369) 
(  370)        !   Put the temporary moments into appropriate neighboring boxes.    
(  371)        DO jj = 2, jpjm1                    !   Flux from j to j+1 IF v GT 0.
(  372)           DO ji = 2, jpim1
(  373)              zbt  =         zbet(ji,jj-1)
(  374)              zbt1 = ( 1.0 - zbet(ji,jj-1) )
(  375)              psm(ji,jj)  = zbt * ( psm(ji,jj) + zfm(ji,jj-1) ) + zbt1 * psm(ji,jj) 
(  376)              zalf        = zbt * zfm(ji,jj-1) / psm(ji,jj) 
(  377)              zalf1       = 1.0 - zalf
(  378)              ztemp       = zalf * ps0(ji,jj) - zalf1 * zf0(ji,jj-1)
(  379)              ps0(ji,jj)  = zbt * (ps0(ji,jj) + zf0(ji,jj-1)) + zbt1 * ps0(ji,jj)
(  380) 
(  381)              psy(ji,jj)  = zbt  * ( zalf * zfy(ji,jj-1) + zalf1 * psy(ji,jj) + 3.0 * ztemp )   &
(  382)                 &        + zbt1 * psy(ji,jj)  
(  383) 
(  384)              psyy(ji,jj) = zbt  * ( zalf * zalf * zfyy(ji,jj-1) + zalf1 * zalf1 * psyy(ji,jj)                             &
(  385)                 &                 + 5.0 * ( zalf * zalf1 * ( psy(ji,jj) - zfy(ji,jj-1) ) - ( zalf1 - zalf ) * ztemp ) )   & 
(  386)                 &        + zbt1 * psyy(ji,jj)
(  387) 
(  388)              psxy(ji,jj) = zbt  * ( zalf * zfxy(ji,jj-1) + zalf1 * psxy(ji,jj)              &
(  389)                                   + 3.0 * (- zalf1 * zfx(ji,jj-1) + zalf * psx(ji,jj) ) )   &
(  390)                          + zbt1 * psxy(ji,jj)
(  391)              psx (ji,jj) = zbt * ( psx (ji,jj) + zfx (ji,jj-1) ) + zbt1 * psx (ji,jj)
(  392)              psxx(ji,jj) = zbt * ( psxx(ji,jj) + zfxx(ji,jj-1) ) + zbt1 * psxx(ji,jj)
(  393)           END DO
(  394)        END DO






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 8

(  395) 
(  396)        DO jj = 2, jpjm1                   !  Flux from j+1 to j IF v LT 0.
(  397)           DO ji = 2, jpim1
(  398)              zbt  =         zbet(ji,jj)
(  399)              zbt1 = ( 1.0 - zbet(ji,jj) )
(  400)              psm(ji,jj)  = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) + zfm(ji,jj) )
(  401)              zalf        = zbt1 * zfm(ji,jj) / psm(ji,jj)
(  402)              zalf1       = 1.0 - zalf
(  403)              ztemp       = -zalf * ps0(ji,jj) + zalf1 * zf0(ji,jj)
(  404)              ps0(ji,jj)  = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) + zf0(ji,jj) )
(  405)              psy(ji,jj)  = zbt  * psy(ji,jj)  &
(  406)                 &        + zbt1 * ( zalf*zfy(ji,jj) + zalf1 * psy(ji,jj) + 3.0 * ztemp )
(  407)              psyy(ji,jj) = zbt  * psyy(ji,jj)  &
(  408)                 &        + zbt1 * ( zalf * zalf * zfyy(ji,jj) + zalf1 * zalf1 * psyy(ji,jj)   &
(  409)                 &                 + 5.0 *( zalf *zalf1 *( -psy(ji,jj) + zfy(ji,jj) ) + ( zalf1 - zalf ) * ztemp ) )
(  410)              psxy(ji,jj) = zbt  * psxy(ji,jj)   &
(  411)                 &        + zbt1 * ( zalf * zfxy(ji,jj) + zalf1 * psxy(ji,jj)   &
(  412)                 &                 + 3.0 * ( zalf1 * zfx(ji,jj) - zalf * psx(ji,jj) ) )
(  413)              psx(ji,jj)  = zbt * psx (ji,jj) + zbt1 * ( psx (ji,jj) + zfx (ji,jj) )
(  414)              psxx(ji,jj) = zbt * psxx(ji,jj) + zbt1 * ( psxx(ji,jj) + zfxx(ji,jj) )
(  415)           END DO
(  416)        END DO
(  417) 
(  418)       !-- Lateral boundary conditions
(  419)       CALL lbc_lnk( psm , 'T', 1. )
(  420)       CALL lbc_lnk( ps0 , 'T', 1. )
(  421)       CALL lbc_lnk( psx , 'T', 1. )
(  422)       CALL lbc_lnk( psxx, 'T', 1. )
(  423)       CALL lbc_lnk( psy , 'T', 1. )
(  424)       CALL lbc_lnk( psyy, 'T', 1. )
(  425)       CALL lbc_lnk( psxy, 'T', 1. )
(  426) 
(  427)       IF(ln_ctl) THEN
(  428)          CALL prt_ctl(tab2d_1=psm  , clinfo1=' lim_adv_y: psm  :', tab2d_2=ps0 , clinfo2=' ps0  : ')
(  429)          CALL prt_ctl(tab2d_1=psx  , clinfo1=' lim_adv_y: psx  :', tab2d_2=psxx, clinfo2=' psxx : ')
(  430)          CALL prt_ctl(tab2d_1=psy  , clinfo1=' lim_adv_y: psy  :', tab2d_2=psyy, clinfo2=' psyy : ')
(  431)          CALL prt_ctl(tab2d_1=psxy , clinfo1=' lim_adv_y: psxy :')
(  432)       ENDIF
(  433) 
(  434)    END SUBROUTINE lim_adv_y
(  435) 
(  436) #else
(  437)    !!----------------------------------------------------------------------
(  438)    !!   Default option            Dummy module         NO LIM sea-ice model
(  439)    !!----------------------------------------------------------------------
(  440) CONTAINS
(  441)    SUBROUTINE lim_adv_x         ! Empty routine
(  442)    END SUBROUTINE lim_adv_x
(  443)    SUBROUTINE lim_adv_y           ! Empty routine
(  444)    END SUBROUTINE lim_adv_y
(  445) 
(  446) #endif
(  447) 
(  448) END MODULE limadv










PGF90 (Version     10.5)          09/23/2010  09:49:54      page 9

(    1) # 1 "limadv.F90"
(    1) MODULE limadv 
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limadv   ***
(    4)    !! LIM sea-ice model : sea-ice advection
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_adv_x  : advection of sea ice on x axis
(   11)    !!   lim_adv_y  : advection of sea ice on y axis
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE dom_oce
(   15)    USE dom_ice
(   16)    USE ice_oce         ! ice variables
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE lbclnk
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC lim_adv_x    ! called by lim_trp
(   26)    PUBLIC lim_adv_y    ! called by lim_trp
(   27) 
(   28)    !! * Module variables
(   29)    REAL(wp)  ::            &  ! constant values
(   30)       epsi20 = 1e-20    ,  &
(   31)       rzero  = 0.e0     ,  &
(   32)       rone   = 1.e0
(   33)    !!----------------------------------------------------------------------
(   34)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limadv.F90,v 1.5 2005/09/02 15:45:16 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE lim_adv_x( pdf, put , pcrh, psm , ps0 ,   &
(   42)       &                  psx, psxx, psy , psyy, psxy )
(   43)       !!---------------------------------------------------------------------
(   44)       !!                **  routine lim_adv_x  **
(   45)       !!  
(   46)       !! ** purpose :   Computes and adds the advection trend to sea-ice
(   47)       !!      variable on x axis
(   48)       !!
(   49)       !! ** method  :   Uses Prather second order scheme that advects
(   50)       !!      tracers but also theirquadratic forms. The method preserves
(   51)       !!      tracer structures by conserving second order moments.
(   52)       !!      Ref.: "Numerical Advection by Conservation of Second Order
(   53)       !!      Moments", JGR, VOL. 91. NO. D6. PAGES 6671-6681. MAY 20, 1986
(   54)       !!     
(   55)       !! History :
(   56)       !!        !  00-01 (LIM)
(   57)       !!        !  01-05 (G. Madec, R. Hordoir) opa norm






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 10

(   58)       !!        !  03-10 (C. Ethe) F90, module
(   59)       !!        !  03-12 (R. Hordoir, G. Madec) mpp
(   60)       !!--------------------------------------------------------------------
(   61)       !! * Arguments
(   62)       REAL(wp)                    , INTENT(in)  ::  &
(   63)          pdf ,       &  ! ???
(   64)          pcrh           ! = 1. : lim_adv_x is called before lim_adv_y
(   65)          !              ! = 0. : lim_adv_x is called after  lim_adv_y
(   66)       REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  ::  &
(   67)          put            ! i-direction ice velocity at ocean U-point (m/s)
(   68)       REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::  & 
(   69)          ps0 , psm , &  ! ???
(   70)          psx , psy , &  ! ??? 
(   71)          psxx, psyy, psxy
(   72) 
(   73)       !! * Local declarations
(   74)       INTEGER ::   ji, jj      ! dummy loop indices
(   75)       REAL(wp)       ::  &
(   76)          zrdt, zslpmax, ztemp, zin0,     &  ! temporary scalars
(   77)          zs1max, zs1new, zs2new,         &  !    "         "
(   78)          zalf, zalfq, zalf1, zalf1q,     &  !    "         "
(   79)          zbt , zbt1                         !    "         "
(   80)       REAL(wp), DIMENSION(jpi,jpj)  ::   &  ! temporary workspace
(   81)          zf0 , zfx , zfy , zbet,         &  !    "           "
(   82)          zfxx, zfyy, zfxy,               &  !    "           "
(   83)          zfm, zalg, zalg1, zalg1q           !    "           "
(   84)       !---------------------------------------------------------------------
(   85) 
(   86)       ! Limitation of moments.                                           
(   87) 
(   88)       zrdt      = rdt_ice * pdf   ! If ice drift field is too fast, use an appropriate time step for advection.
(   89) 
(   90)       DO jj = 1, jpj
(   91)          DO ji = 1, jpi
(   92)             zslpmax = MAX( rzero, ps0(ji,jj) )
(   93)             zs1max  = 1.5 * zslpmax
(   94)             zs1new  = MIN( zs1max, MAX( -zs1max, psx(ji,jj) ) )
(   95)             zs2new  = MIN(  2.0 * zslpmax - 0.3334 * ABS( zs1new ),      &
(   96)                &            MAX( ABS( zs1new ) - zslpmax, psxx(ji,jj) )  )
(   97)             zin0    = ( 1.0 - MAX( rzero, sign ( rone, -zslpmax) ) ) * tms(ji,jj)   ! Case of empty boxes & Apply mask
(   98) 
(   99)             ps0 (ji,jj) = zslpmax  
(  100)             psx (ji,jj) = zs1new      * zin0
(  101)             psxx(ji,jj) = zs2new      * zin0
(  102)             psy (ji,jj) = psy (ji,jj) * zin0
(  103)             psyy(ji,jj) = psyy(ji,jj) * zin0
(  104)             psxy(ji,jj) = MIN( zslpmax, MAX( -zslpmax, psxy(ji,jj) ) ) * zin0
(  105)          END DO
(  106)       END DO
(  107) 
(  108)       !  Initialize volumes of boxes  (=area if adv_x first called, =psm otherwise)                                     
(  109)       psm (:,:)  = MAX( pcrh * area(:,:) + ( 1.0 - pcrh ) * psm(:,:) , epsi20 )
(  110) 
(  111)       !  Calculate fluxes and moments between boxes i<-->i+1              
(  112)       DO jj = 2, jpjm1                    !  Flux from i to i+1 WHEN u GT 0 
(  113) !i bug   DO ji = 1, jpim1 
(  114) !i    DO jj = 1, jpj                      !  Flux from i to i+1 WHEN u GT 0 
(  115)          DO ji = 1, jpi






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 11

(  116)             zbet(ji,jj)  =  MAX( rzero, SIGN( rone, put(ji,jj) ) )
(  117)             zalf         =  MAX( rzero, put(ji,jj) ) * zrdt * e2u(ji,jj) / psm(ji,jj)
(  118)             zalfq        =  zalf * zalf
(  119)             zalf1        =  1.0 - zalf
(  120)             zalf1q       =  zalf1 * zalf1
(  121)             zfm (ji,jj)  =  zalf  * psm(ji,jj)
(  122)             zf0 (ji,jj)  =  zalf  * ( ps0(ji,jj) + zalf1 * ( psx(ji,jj) + (zalf1 - zalf) * psxx(ji,jj)  ) )
(  123)             zfx (ji,jj)  =  zalfq * ( psx(ji,jj) + 3.0 * zalf1 * psxx(ji,jj) )
(  124)             zfxx(ji,jj)  =  zalf  * zalfq * psxx(ji,jj)
(  125)             zfy (ji,jj)  =  zalf  * ( psy(ji,jj) + zalf1 * psxy(ji,jj) )
(  126)             zfxy(ji,jj)  =  zalfq * psxy(ji,jj)
(  127)             zfyy(ji,jj)  =  zalf  * psyy(ji,jj)
(  128) 
(  129)             !  Readjust moments remaining in the box.
(  130)             psm (ji,jj)  =  psm (ji,jj) - zfm(ji,jj)
(  131)             ps0 (ji,jj)  =  ps0 (ji,jj) - zf0(ji,jj)
(  132)             psx (ji,jj)  =  zalf1q * ( psx(ji,jj) - 3.0 * zalf * psxx(ji,jj) )
(  133)             psxx(ji,jj)  =  zalf1  * zalf1q * psxx(ji,jj)
(  134)             psy (ji,jj)  =  psy (ji,jj) - zfy(ji,jj)
(  135)             psyy(ji,jj)  =  psyy(ji,jj) - zfyy(ji,jj)
(  136)             psxy(ji,jj)  =  zalf1q * psxy(ji,jj)
(  137)          END DO
(  138)       END DO
(  139) 
(  140)       DO jj = 2, jpjm1                      !  Flux from i+1 to i when u LT 0.
(  141) !i    DO jj = 1, jpjm1                      !  Flux from i+1 to i when u LT 0.
(  142)          DO ji = 1, jpim1
(  143)             zalf          = MAX( rzero, -put(ji,jj) ) * zrdt * e2u(ji,jj) / psm(ji+1,jj) 
(  144)             zalg  (ji,jj) = zalf
(  145)             zalfq         = zalf * zalf
(  146)             zalf1         = 1.0 - zalf
(  147)             zalg1 (ji,jj) = zalf1
(  148)             zalf1q        = zalf1 * zalf1
(  149)             zalg1q(ji,jj) = zalf1q
(  150)             zfm   (ji,jj) = zfm (ji,jj) + zalf  * psm(ji+1,jj)
(  151)             zf0   (ji,jj) = zf0 (ji,jj) + zalf  * ( ps0(ji+1,jj) - zalf1 * ( psx(ji+1,jj) - (zalf1 - zalf ) * psxx(ji+1,jj) ) )
(  152)             zfx   (ji,jj) = zfx (ji,jj) + zalfq * ( psx(ji+1,jj) - 3.0 * zalf1 * psxx(ji+1,jj) )
(  153)             zfxx  (ji,jj) = zfxx(ji,jj) + zalf  * zalfq * psxx(ji+1,jj)
(  154)             zfy   (ji,jj) = zfy (ji,jj) + zalf  * ( psy(ji+1,jj) - zalf1 * psxy(ji+1,jj) )
(  155)             zfxy  (ji,jj) = zfxy(ji,jj) + zalfq * psxy(ji+1,jj)
(  156)             zfyy  (ji,jj) = zfyy(ji,jj) + zalf  * psyy(ji+1,jj)
(  157)          END DO
(  158)       END DO
(  159) 
(  160)       DO jj = 2, jpjm1                     !  Readjust moments remaining in the box. 
(  161)          DO ji = 2, jpim1
(  162)             zbt  =       zbet(ji-1,jj)
(  163)             zbt1 = 1.0 - zbet(ji-1,jj)
(  164)             psm (ji,jj) = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) - zfm(ji-1,jj) )
(  165)             ps0 (ji,jj) = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) - zf0(ji-1,jj) )
(  166)             psx (ji,jj) = zalg1q(ji-1,jj) * ( psx(ji,jj) + 3.0 * zalg(ji-1,jj) * psxx(ji,jj) )
(  167)             psxx(ji,jj) = zalg1 (ji-1,jj) * zalg1q(ji-1,jj) * psxx(ji,jj)
(  168)             psy (ji,jj) = zbt * psy (ji,jj) + zbt1 * ( psy (ji,jj) - zfy (ji-1,jj) )
(  169)             psyy(ji,jj) = zbt * psyy(ji,jj) + zbt1 * ( psyy(ji,jj) - zfyy(ji-1,jj) )
(  170)             psxy(ji,jj) = zalg1q(ji-1,jj) * psxy(ji,jj)
(  171)          END DO
(  172)       END DO
(  173) 






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 12

(  174)       !   Put the temporary moments into appropriate neighboring boxes.    
(  175)       DO jj = 2, jpjm1                     !   Flux from i to i+1 IF u GT 0.
(  176)          DO ji = 2, jpim1
(  177)             zbt  =       zbet(ji-1,jj)
(  178)             zbt1 = 1.0 - zbet(ji-1,jj)
(  179)             psm(ji,jj)  = zbt * ( psm(ji,jj) + zfm(ji-1,jj) ) + zbt1 * psm(ji,jj)
(  180)             zalf        = zbt * zfm(ji-1,jj) / psm(ji,jj)
(  181)             zalf1       = 1.0 - zalf
(  182)             ztemp       = zalf * ps0(ji,jj) - zalf1 * zf0(ji-1,jj)
(  183)             ps0(ji,jj)  = zbt * (ps0(ji,jj) + zf0(ji-1,jj)) + zbt1 * ps0(ji,jj)
(  184)             psx(ji,jj)  = zbt * ( zalf * zfx(ji-1,jj) + zalf1 * psx(ji,jj) + 3.0 * ztemp ) + zbt1 * psx(ji,jj)
(  185)             psxx(ji,jj) = zbt * ( zalf * zalf * zfxx(ji-1,jj) + zalf1 * zalf1 * psxx(ji,jj)   &
(  186)                &                + 5.0 * ( zalf * zalf1 * ( psx (ji,jj) - zfx(ji-1,jj) ) - ( zalf1 - zalf ) * ztemp )  )   &
(  187)                &        + zbt1 * psxx(ji,jj)
(  188)             psxy(ji,jj) = zbt * ( zalf * zfxy(ji-1,jj) + zalf1 * psxy(ji,jj)             &
(  189)                &                + 3.0 * (- zalf1*zfy(ji-1,jj)  + zalf * psy(ji,jj) ) )   &
(  190)                &         + zbt1 * psxy(ji,jj)
(  191)             psy (ji,jj) = zbt * ( psy (ji,jj) + zfy (ji-1,jj) ) + zbt1 * psy (ji,jj)
(  192)             psyy(ji,jj) = zbt * ( psyy(ji,jj) + zfyy(ji-1,jj) ) + zbt1 * psyy(ji,jj)
(  193)          END DO
(  194)       END DO
(  195) 
(  196)       DO jj = 2, jpjm1                     !  Flux from i+1 to i IF u LT 0.
(  197)          DO ji = 2, jpim1
(  198)             zbt  =       zbet(ji,jj)
(  199)             zbt1 = 1.0 - zbet(ji,jj)
(  200)             psm(ji,jj)  = zbt * psm(ji,jj)  + zbt1 * ( psm(ji,jj) + zfm(ji,jj) )
(  201)             zalf        = zbt1 * zfm(ji,jj) / psm(ji,jj)
(  202)             zalf1       = 1.0 - zalf
(  203)             ztemp       = -zalf * ps0(ji,jj) + zalf1 * zf0(ji,jj)
(  204)             ps0(ji,jj)  = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) + zf0(ji,jj) )
(  205)             psx(ji,jj)  = zbt  * psx(ji,jj)   &
(  206)                &        + zbt1 * ( zalf * zfx(ji,jj) + zalf1 * psx(ji,jj) + 3.0 * ztemp )
(  207)             psxx(ji,jj) = zbt  * psxx(ji,jj)   &
(  208)                &        + zbt1 * (  zalf * zalf * zfxx(ji,jj)  + zalf1 * zalf1 * psxx(ji,jj)  &
(  209)                &                 + 5.0 *( zalf * zalf1 * ( - psx(ji,jj) + zfx(ji,jj) ) + ( zalf1 - zalf ) * ztemp ) )
(  210)             psxy(ji,jj) = zbt  * psxy(ji,jj)   &
(  211)                &        + zbt1 * (  zalf * zfxy(ji,jj) + zalf1 * psxy(ji,jj)  &
(  212)                &                 + 3.0 * ( zalf1 * zfy(ji,jj) - zalf * psy(ji,jj) )  )
(  213)             psy(ji,jj)  = zbt * psy (ji,jj)  + zbt1 * ( psy (ji,jj) + zfy (ji,jj) )
(  214)             psyy(ji,jj) = zbt * psyy(ji,jj)  + zbt1 * ( psyy(ji,jj) + zfyy(ji,jj) )
(  215)          END DO
(  216)       END DO
(  217) 
(  218)       !-- Lateral boundary conditions
(  219)       CALL lbc_lnk( psm , 'T', 1. )
(  220)       CALL lbc_lnk( ps0 , 'T', 1. )
(  221)       CALL lbc_lnk( psx , 'T', 1. )
(  222)       CALL lbc_lnk( psxx, 'T', 1. )
(  223)       CALL lbc_lnk( psy , 'T', 1. )
(  224)       CALL lbc_lnk( psyy, 'T', 1. )
(  225)       CALL lbc_lnk( psxy, 'T', 1. )
(  226) 
(  227)       IF(ln_ctl)   THEN
(  228)          CALL prt_ctl(tab2d_1=psm  , clinfo1=' lim_adv_x: psm  :', tab2d_2=ps0 , clinfo2=' ps0  : ')
(  229)          CALL prt_ctl(tab2d_1=psx  , clinfo1=' lim_adv_x: psx  :', tab2d_2=psxx, clinfo2=' psxx : ')
(  230)          CALL prt_ctl(tab2d_1=psy  , clinfo1=' lim_adv_x: psy  :', tab2d_2=psyy, clinfo2=' psyy : ')
(  231)          CALL prt_ctl(tab2d_1=psxy , clinfo1=' lim_adv_x: psxy :')






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 13

(  232)       ENDIF
(  233) 
(  234)    END SUBROUTINE lim_adv_x
(  235) 
(  236) 
(  237)    SUBROUTINE lim_adv_y( pdf, pvt , pcrh, psm , ps0 ,   &
(  238)       &                  psx, psxx, psy , psyy, psxy )
(  239)       !!---------------------------------------------------------------------
(  240)       !!                **  routine lim_adv_y  **
(  241)       !!            
(  242)       !! ** purpose :   Computes and adds the advection trend to sea-ice 
(  243)       !!      variable on y axis
(  244)       !!
(  245)       !! ** method  :   Uses Prather second order scheme that advects tracers
(  246)       !!      but also their quadratic forms. The method preserves tracer
(  247)       !!      structures by conserving second order moments.
(  248)       !!
(  249)       !! History :
(  250)       !!   1.0  !  00-01 (LIM)
(  251)       !!        !  01-05 (G. Madec, R. Hordoir) opa norm
(  252)       !!   2.0  !  03-10 (C. Ethe) F90, module
(  253)       !!        !  03-12 (R. Hordoir, G. Madec) mpp
(  254)       !!---------------------------------------------------------------------
(  255)       !! * Arguments
(  256)       REAL(wp),                     INTENT(in)  :: &
(  257)          pdf,        &  ! ???
(  258)          pcrh           ! = 1. : lim_adv_x is called before lim_adv_y
(  259)          !              ! = 0. : lim_adv_x is called after  lim_adv_y
(  260)       REAL(wp), DIMENSION(jpi,jpj), INTENT(in)  :: &
(  261)          pvt            ! j-direction ice velocity at ocean V-point (m/s)
(  262)       REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) :: &
(  263)          psm , ps0 , psx , psy,   &
(  264)          psxx, psyy, psxy
(  265) 
(  266)       !! * Local Variables
(  267)       INTEGER  ::   ji, jj                ! dummy loop indices
(  268)       REAL(wp) ::   &
(  269)          zrdt, zslpmax, zin0, ztemp,  &   ! temporary scalars
(  270)          zs1max, zs1new, zs2new,      &   !    "         "
(  271)          zalf, zalfq, zalf1, zalf1q,  &   !    "         "
(  272)          zbt , zbt1                       !
(  273)       REAL(wp), DIMENSION(jpi,jpj)  :: &
(  274)          zf0 , zfx , zfy ,             &  ! temporary workspace
(  275)          zfxx, zfyy, zfxy,             &  !    "           "
(  276)          zfm , zbet,                   &  !    "           "
(  277)          zalg, zalg1, zalg1q              !    "           "
(  278)       !---------------------------------------------------------------------
(  279) 
(  280)       ! Limitation of moments.
(  281) 
(  282)       zrdt = rdt_ice * pdf ! If ice drift field is too fast, use an appropriate time step for advection.
(  283) 
(  284)        DO jj = 1, jpj
(  285)           DO ji = 1, jpi
(  286)              zslpmax = MAX( rzero, ps0(ji,jj) )
(  287)              zs1max  = 1.5 * zslpmax
(  288)              zs1new  = MIN( zs1max, MAX( -zs1max, psy(ji,jj) ) )
(  289)              zs2new  = MIN(  ( 2.0 * zslpmax - 0.3334 * ABS( zs1new ) ),   &






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 14

(  290)                 &             MAX( ABS( zs1new )-zslpmax, psyy(ji,jj) )  )
(  291)              zin0    = ( 1.0 - MAX( rzero, sign ( rone, -zslpmax) ) ) * tms(ji,jj)   ! Case of empty boxes & Apply mask
(  292)              ps0 (ji,jj) = zslpmax  
(  293)              psx (ji,jj)  = psx (ji,jj) * zin0
(  294)              psxx(ji,jj)  = psxx(ji,jj) * zin0
(  295)              psy (ji,jj) = zs1new * zin0
(  296)              psyy(ji,jj) = zs2new * zin0
(  297)              psxy(ji,jj) = MIN( zslpmax, MAX( -zslpmax, psxy(ji,jj) ) ) * zin0
(  298)           END DO
(  299)        END DO
(  300) 
(  301)        !  Initialize volumes of boxes (=area if adv_x first called, =psm otherwise)
(  302)        psm (:,:)  = MAX( pcrh * area(:,:) + ( 1.0 - pcrh ) * psm(:,:) , epsi20 )
(  303) 
(  304)        !  Calculate fluxes and moments between boxes j<-->j+1              
(  305) !!bug  DO jj = 2, jpjm1
(  306)        DO jj = 1, jpj
(  307)           DO ji = 2, jpim1
(  308) !!bug     DO ji = 1, jpim1
(  309)              !  Flux from j to j+1 WHEN v GT 0   
(  310)              zbet(ji,jj)  =  MAX( rzero, SIGN( rone, pvt(ji,jj) ) )
(  311)              zalf         =  MAX( rzero, pvt(ji,jj) ) * zrdt * e1v(ji,jj) / psm(ji,jj)
(  312)              zalfq        =  zalf * zalf
(  313)              zalf1        =  1.0 - zalf
(  314)              zalf1q       =  zalf1 * zalf1
(  315)              zfm (ji,jj)  =  zalf  * psm(ji,jj)
(  316)              zf0 (ji,jj)  =  zalf  * ( ps0(ji,jj) + zalf1 * ( psy(ji,jj)  + (zalf1-zalf) * psyy(ji,jj)  ) ) 
(  317)              zfy (ji,jj)  =  zalfq *( psy(ji,jj) + 3.0*zalf1*psyy(ji,jj) )
(  318)              zfyy(ji,jj)  =  zalf  * zalfq * psyy(ji,jj)
(  319)              zfx (ji,jj)  =  zalf  * ( psx(ji,jj) + zalf1 * psxy(ji,jj) )
(  320)              zfxy(ji,jj)  =  zalfq * psxy(ji,jj)
(  321)              zfxx(ji,jj)  =  zalf  * psxx(ji,jj)
(  322) 
(  323)              !  Readjust moments remaining in the box.
(  324)              psm (ji,jj)  =  psm (ji,jj) - zfm(ji,jj)
(  325)              ps0 (ji,jj)  =  ps0 (ji,jj) - zf0(ji,jj)
(  326)              psy (ji,jj)  =  zalf1q * ( psy(ji,jj) -3.0 * zalf * psyy(ji,jj) )
(  327)              psyy(ji,jj)  =  zalf1 * zalf1q * psyy(ji,jj)
(  328)              psx (ji,jj)  =  psx (ji,jj) - zfx(ji,jj)
(  329)              psxx(ji,jj)  =  psxx(ji,jj) - zfxx(ji,jj)
(  330)              psxy(ji,jj)  =  zalf1q * psxy(ji,jj)
(  331)           END DO
(  332)        END DO
(  333) 
(  334)        DO jj = 1, jpjm1                   !  Flux from j+1 to j when v LT 0.
(  335)           DO ji = 2, jpim1
(  336) !i     DO jj = 1, jpjm1                   !  Flux from j+1 to j when v LT 0.
(  337) !i        DO ji = 2, jpim1
(  338)              zalf          = ( MAX(rzero, -pvt(ji,jj) ) * zrdt * e1v(ji,jj) ) / psm(ji,jj+1) 
(  339)              zalg  (ji,jj) = zalf
(  340)              zalfq         = zalf * zalf
(  341)              zalf1         = 1.0 - zalf
(  342)              zalg1 (ji,jj) = zalf1
(  343)              zalf1q        = zalf1 * zalf1
(  344)              zalg1q(ji,jj) = zalf1q
(  345)              zfm   (ji,jj) = zfm (ji,jj) + zalf  * psm(ji,jj+1)
(  346)              zf0   (ji,jj) = zf0 (ji,jj) + zalf  * ( ps0(ji,jj+1) - zalf1 * (psy(ji,jj+1) - (zalf1 - zalf ) * psyy(ji,jj+1) ) )
(  347)              zfy   (ji,jj) = zfy (ji,jj) + zalfq * ( psy(ji,jj+1) - 3.0 * zalf1 * psyy(ji,jj+1) )






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 15

(  348)              zfyy  (ji,jj) = zfyy(ji,jj) + zalf  * zalfq * psyy(ji,jj+1)
(  349)              zfx   (ji,jj) = zfx (ji,jj) + zalf  * ( psx(ji,jj+1) - zalf1 * psxy(ji,jj+1) )
(  350)              zfxy  (ji,jj) = zfxy(ji,jj) + zalfq * psxy(ji,jj+1)
(  351)              zfxx  (ji,jj) = zfxx(ji,jj) + zalf  * psxx(ji,jj+1)
(  352)           END DO
(  353)        END DO
(  354)  
(  355)        !  Readjust moments remaining in the box. 
(  356)        DO jj = 2, jpjm1
(  357)           DO ji = 2, jpim1
(  358)              zbt  =         zbet(ji,jj-1)
(  359)              zbt1 = ( 1.0 - zbet(ji,jj-1) )
(  360)              psm (ji,jj) = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) - zfm(ji,jj-1) )
(  361)              ps0 (ji,jj) = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) - zf0(ji,jj-1) )
(  362)              psy (ji,jj) = zalg1q(ji,jj-1) * ( psy(ji,jj) + 3.0 * zalg(ji,jj-1) * psyy(ji,jj) )
(  363)              psyy(ji,jj) = zalg1 (ji,jj-1)  * zalg1q(ji,jj-1) * psyy(ji,jj)
(  364)              psx (ji,jj) = zbt * psx (ji,jj) + zbt1 * ( psx (ji,jj) - zfx (ji,jj-1) )
(  365)              psxx(ji,jj) = zbt * psxx(ji,jj) + zbt1 * ( psxx(ji,jj) - zfxx(ji,jj-1) )
(  366)              psxy(ji,jj) = zalg1q(ji,jj-1) * psxy(ji,jj)
(  367)           END DO
(  368)        END DO
(  369) 
(  370)        !   Put the temporary moments into appropriate neighboring boxes.    
(  371)        DO jj = 2, jpjm1                    !   Flux from j to j+1 IF v GT 0.
(  372)           DO ji = 2, jpim1
(  373)              zbt  =         zbet(ji,jj-1)
(  374)              zbt1 = ( 1.0 - zbet(ji,jj-1) )
(  375)              psm(ji,jj)  = zbt * ( psm(ji,jj) + zfm(ji,jj-1) ) + zbt1 * psm(ji,jj) 
(  376)              zalf        = zbt * zfm(ji,jj-1) / psm(ji,jj) 
(  377)              zalf1       = 1.0 - zalf
(  378)              ztemp       = zalf * ps0(ji,jj) - zalf1 * zf0(ji,jj-1)
(  379)              ps0(ji,jj)  = zbt * (ps0(ji,jj) + zf0(ji,jj-1)) + zbt1 * ps0(ji,jj)
(  380) 
(  381)              psy(ji,jj)  = zbt  * ( zalf * zfy(ji,jj-1) + zalf1 * psy(ji,jj) + 3.0 * ztemp )   &
(  382)                 &        + zbt1 * psy(ji,jj)  
(  383) 
(  384)              psyy(ji,jj) = zbt  * ( zalf * zalf * zfyy(ji,jj-1) + zalf1 * zalf1 * psyy(ji,jj)                             &
(  385)                 &                 + 5.0 * ( zalf * zalf1 * ( psy(ji,jj) - zfy(ji,jj-1) ) - ( zalf1 - zalf ) * ztemp ) )   & 
(  386)                 &        + zbt1 * psyy(ji,jj)
(  387) 
(  388)              psxy(ji,jj) = zbt  * ( zalf * zfxy(ji,jj-1) + zalf1 * psxy(ji,jj)              &
(  389)                                   + 3.0 * (- zalf1 * zfx(ji,jj-1) + zalf * psx(ji,jj) ) )   &
(  390)                          + zbt1 * psxy(ji,jj)
(  391)              psx (ji,jj) = zbt * ( psx (ji,jj) + zfx (ji,jj-1) ) + zbt1 * psx (ji,jj)
(  392)              psxx(ji,jj) = zbt * ( psxx(ji,jj) + zfxx(ji,jj-1) ) + zbt1 * psxx(ji,jj)
(  393)           END DO
(  394)        END DO
(  395) 
(  396)        DO jj = 2, jpjm1                   !  Flux from j+1 to j IF v LT 0.
(  397)           DO ji = 2, jpim1
(  398)              zbt  =         zbet(ji,jj)
(  399)              zbt1 = ( 1.0 - zbet(ji,jj) )
(  400)              psm(ji,jj)  = zbt * psm(ji,jj) + zbt1 * ( psm(ji,jj) + zfm(ji,jj) )
(  401)              zalf        = zbt1 * zfm(ji,jj) / psm(ji,jj)
(  402)              zalf1       = 1.0 - zalf
(  403)              ztemp       = -zalf * ps0(ji,jj) + zalf1 * zf0(ji,jj)
(  404)              ps0(ji,jj)  = zbt * ps0(ji,jj) + zbt1 * ( ps0(ji,jj) + zf0(ji,jj) )
(  405)              psy(ji,jj)  = zbt  * psy(ji,jj)  &






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 16

(  406)                 &        + zbt1 * ( zalf*zfy(ji,jj) + zalf1 * psy(ji,jj) + 3.0 * ztemp )
(  407)              psyy(ji,jj) = zbt  * psyy(ji,jj)  &
(  408)                 &        + zbt1 * ( zalf * zalf * zfyy(ji,jj) + zalf1 * zalf1 * psyy(ji,jj)   &
(  409)                 &                 + 5.0 *( zalf *zalf1 *( -psy(ji,jj) + zfy(ji,jj) ) + ( zalf1 - zalf ) * ztemp ) )
(  410)              psxy(ji,jj) = zbt  * psxy(ji,jj)   &
(  411)                 &        + zbt1 * ( zalf * zfxy(ji,jj) + zalf1 * psxy(ji,jj)   &
(  412)                 &                 + 3.0 * ( zalf1 * zfx(ji,jj) - zalf * psx(ji,jj) ) )
(  413)              psx(ji,jj)  = zbt * psx (ji,jj) + zbt1 * ( psx (ji,jj) + zfx (ji,jj) )
(  414)              psxx(ji,jj) = zbt * psxx(ji,jj) + zbt1 * ( psxx(ji,jj) + zfxx(ji,jj) )
(  415)           END DO
(  416)        END DO
(  417) 
(  418)       !-- Lateral boundary conditions
(  419)       CALL lbc_lnk( psm , 'T', 1. )
(  420)       CALL lbc_lnk( ps0 , 'T', 1. )
(  421)       CALL lbc_lnk( psx , 'T', 1. )
(  422)       CALL lbc_lnk( psxx, 'T', 1. )
(  423)       CALL lbc_lnk( psy , 'T', 1. )
(  424)       CALL lbc_lnk( psyy, 'T', 1. )
(  425)       CALL lbc_lnk( psxy, 'T', 1. )
(  426) 
(  427)       IF(ln_ctl) THEN
(  428)          CALL prt_ctl(tab2d_1=psm  , clinfo1=' lim_adv_y: psm  :', tab2d_2=ps0 , clinfo2=' ps0  : ')
(  429)          CALL prt_ctl(tab2d_1=psx  , clinfo1=' lim_adv_y: psx  :', tab2d_2=psxx, clinfo2=' psxx : ')
(  430)          CALL prt_ctl(tab2d_1=psy  , clinfo1=' lim_adv_y: psy  :', tab2d_2=psyy, clinfo2=' psyy : ')
(  431)          CALL prt_ctl(tab2d_1=psxy , clinfo1=' lim_adv_y: psxy :')
(  432)       ENDIF
(  433) 
(  434)    END SUBROUTINE lim_adv_y
(  435) 
(  436) 
(  437) # 448
(  448) END MODULE limadv
