


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

Filename: limrhg.F90

(    1) MODULE limrhg
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  limrhg  ***
(    4)    !!   Ice rheology :  performs sea ice rheology
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_rhg   : computes ice velocities
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE phycst
(   14)    USE par_oce
(   15)    USE ice_oce         ! ice variables
(   16)    USE dom_ice
(   17)    USE ice
(   18)    USE lbclnk
(   19)    USE lib_mpp
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE prtctl          ! Print control
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC lim_rhg  ! routine called by lim_dyn
(   28) 
(   29)    !! * Module variables
(   30)    REAL(wp)  ::           &  ! constant values
(   31)       rzero   = 0.e0   ,  &
(   32)       rone    = 1.e0
(   33)    !!----------------------------------------------------------------------
(   34)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limrhg.F90,v 1.7 2006/03/10 10:35:42 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE lim_rhg( k_j1, k_jpj )
(   42)       !!-------------------------------------------------------------------
(   43)       !!                 ***  SUBROUTINR lim_rhg  ***
(   44)       !!
(   45)       !! ** purpose :   determines the velocity field of sea ice by using
(   46)       !!  atmospheric (wind stress) and oceanic (water stress and surface






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 2

(   47)       !!  tilt) forcings. Ice-ice interaction is described by a non-linear
(   48)       !!  viscous-plastic law including shear strength and a bulk rheology.
(   49)       !!
(   50)       !! ** Action  : - compute u_ice, v_ice the sea-ice velocity
(   51)       !!
(   52)       !! History :
(   53)       !!   0.0  !  93-12  (M.A. Morales Maqueda.)  Original code
(   54)       !!   1.0  !  94-12  (H. Goosse) 
(   55)       !!   2.0  !  03-12  (C. Ethe, G. Madec)  F90, mpp
(   56)       !!-------------------------------------------------------------------
(   57)       ! * Arguments
(   58)       INTEGER, INTENT(in) ::   k_j1 ,  &  ! southern j-index for ice computation
(   59)          &                     k_jpj      ! northern j-index for ice computation
(   60) 
(   61)       ! * Local variables
(   62)       INTEGER ::   ji, jj              ! dummy loop indices
(   63) 
(   64)       INTEGER  :: &
(   65)          iim1, ijm1, iip1 , ijp1   , & ! temporary integers
(   66)          iter, jter                    !    "          "
(   67) 
(   68)       CHARACTER (len=50) :: charout
(   69) 
(   70)       REAL(wp) :: &
(   71)          ze11  , ze12  , ze22  , ze21  ,   &  ! temporary scalars
(   72)          zt11  , zt12  , zt21  , zt22  ,   &  !    "         "
(   73)          zvis11, zvis21, zvis12, zvis22,   &  !    "         "
(   74)          zgphsx, ztagnx, zusw  ,           &  !    "         "
(   75)          zgphsy, ztagny                       !    "         "
(   76)       REAL(wp) :: &
(   77)          zresm, zunw, zvnw, zur, zvr, zmod, za, zac, &
(   78)          zmpzas, zstms, zindu, zindu1, zusdtp, zmassdt, zcorlal,  &
(   79)          ztrace2, zdeter, zdelta, zsang, zmask, zdgp, zdgi, zdiag
(   80)       REAL(wp),DIMENSION(jpi,jpj) ::   &
(   81)          zpresh, zfrld, zmass, zcorl,     &
(   82)          zu0, zv0, zviszeta, zviseta,     &
(   83)          zc1u, zc1v, zc2u, zc2v, za1ct, za2ct, za1, za2, zb1, zb2,  &
(   84)          zc1, zc2, zd1, zd2, zden, zu_ice, zv_ice, zresr
(   85)       REAL(wp),DIMENSION(jpi,jpj,2,2) :: &
(   86)          zs11, zs12, zs22, zs21
(   87)       !!-------------------------------------------------------------------
(   88)       
(   89)       !  Store initial velocities
(   90)       !  ------------------------
(   91)       zu0(:,:) = u_ice(:,:)
(   92)       zv0(:,:) = v_ice(:,:)
(   93) 
(   94)       ! Ice mass, ice strength, and wind stress at the center            |
(   95)       ! of the grid squares.                                             |
(   96)       !-------------------------------------------------------------------
(   97) 
(   98)       DO jj = k_j1 , k_jpj-1
(   99)          DO ji = 1 , jpi
(  100)             za1(ji,jj)    = tms(ji,jj) * ( rhosn * hsnm(ji,jj) + rhoic * hicm(ji,jj) )
(  101)             zpresh(ji,jj) = tms(ji,jj) *  pstarh * hicm(ji,jj) * EXP( -c_rhg * frld(ji,jj) )
(  102) #if defined key_lim_cp1 && defined key_coupled
(  103)             zb1(ji,jj)    = tms(ji,jj) * gtaux(ji,jj) * ( 1.0 - frld(ji,jj) )
(  104)             zb2(ji,jj)    = tms(ji,jj) * gtauy(ji,jj) * ( 1.0 - frld(ji,jj) )






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 3

(  105) #else
(  106)             zb1(ji,jj)    = tms(ji,jj) * ( 1.0 - frld(ji,jj) )
(  107)             zb2(ji,jj)    = tms(ji,jj) * ( 1.0 - frld(ji,jj) )
(  108) #endif
(  109)          END DO
(  110)       END DO
(  111) 
(  112) 
(  113)       !---------------------------------------------------------------------------
(  114)       !  Wind stress, coriolis and mass terms at the corners of the grid squares |
(  115)       !  Gradient of ice strenght.                                               |
(  116)       !---------------------------------------------------------------------------
(  117)          
(  118)       DO jj = k_j1+1, k_jpj-1
(  119)          DO ji = 2, jpi
(  120)             zstms = tms(ji,jj  ) * wght(ji,jj,2,2) + tms(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  121)                &  + tms(ji,jj-1) * wght(ji,jj,2,1) + tms(ji-1,jj-1) * wght(ji,jj,1,1)
(  122)             zusw  = 1.0 / MAX( zstms, epsd )
(  123) 
(  124)             zt11 = tms(ji  ,jj  ) * frld(ji  ,jj  ) 
(  125)             zt12 = tms(ji-1,jj  ) * frld(ji-1,jj  ) 
(  126)             zt21 = tms(ji  ,jj-1) * frld(ji  ,jj-1) 
(  127)             zt22 = tms(ji-1,jj-1) * frld(ji-1,jj-1)
(  128) 
(  129)             ! Leads area.
(  130)             zfrld(ji,jj) =  (  zt11 * wght(ji,jj,2,2) + zt12 * wght(ji,jj,1,2)   &
(  131)                &             + zt21 * wght(ji,jj,2,1) + zt22 * wght(ji,jj,1,1) ) * zusw
(  132) 
(  133)             ! Mass and coriolis coeff.
(  134)             zmass(ji,jj) = ( za1(ji,jj  ) * wght(ji,jj,2,2) + za1(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  135)                &           + za1(ji,jj-1) * wght(ji,jj,2,1) + za1(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw
(  136)             zcorl(ji,jj) = zmass(ji,jj) * fcor(ji,jj)
(  137) 
(  138)             ! Wind stress.
(  139) #if defined key_lim_cp1 && defined key_coupled
(  140)             ztagnx = ( zb1(ji,jj  ) * wght(ji,jj,2,2) + zb1(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  141)                &     + zb1(ji,jj-1) * wght(ji,jj,2,1) + zb1(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw
(  142)             ztagny = ( zb2(ji,jj  ) * wght(ji,jj,2,2) + zb2(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  143)                &     + zb2(ji,jj-1) * wght(ji,jj,2,1) + zb2(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw
(  144) #else
(  145)             ztagnx = ( zb1(ji,jj  ) * wght(ji,jj,2,2) + zb1(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  146)                &     + zb1(ji,jj-1) * wght(ji,jj,2,1) + zb1(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw * gtaux(ji,jj)
(  147)             ztagny = ( zb2(ji,jj  ) * wght(ji,jj,2,2) + zb2(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  148)                &     + zb2(ji,jj-1) * wght(ji,jj,2,1) + zb2(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw * gtauy(ji,jj)
(  149) #endif
(  150) 
(  151)             ! Gradient of ice strength
(  152)             zgphsx =   ( alambd(ji,jj,2,2,2,1) - alambd(ji,jj,2,1,2,1) ) * zpresh(ji  ,jj-1)   &
(  153)                &     + ( alambd(ji,jj,2,2,2,2) - alambd(ji,jj,2,1,2,2) ) * zpresh(ji  ,jj  )   &
(  154)                &     - ( alambd(ji,jj,2,2,1,1) + alambd(ji,jj,2,1,1,1) ) * zpresh(ji-1,jj-1)   &
(  155)                &     - ( alambd(ji,jj,2,2,1,2) + alambd(ji,jj,2,1,1,2) ) * zpresh(ji-1,jj  )
(  156) 
(  157)             zgphsy = - ( alambd(ji,jj,1,1,2,1) + alambd(ji,jj,1,2,2,1) ) * zpresh(ji  ,jj-1)   &
(  158)                &     - ( alambd(ji,jj,1,1,1,1) + alambd(ji,jj,1,2,1,1) ) * zpresh(ji-1,jj-1)   &
(  159)                &     + ( alambd(ji,jj,1,1,2,2) - alambd(ji,jj,1,2,2,2) ) * zpresh(ji  ,jj  )   &
(  160)                &     + ( alambd(ji,jj,1,1,1,2) - alambd(ji,jj,1,2,1,2) ) * zpresh(ji-1,jj  )
(  161) 
(  162)             ! Computation of the velocity field taking into account the ice-ice interaction.                                






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 4

(  163)             ! Terms that are independent of the velocity field.
(  164)             za1ct(ji,jj) = ztagnx - zcorl(ji,jj) * v_oce(ji,jj) - zgphsx
(  165)             za2ct(ji,jj) = ztagny + zcorl(ji,jj) * u_oce(ji,jj) - zgphsy
(  166)          END DO
(  167)       END DO
(  168) 
(  169) !! inutile!!
(  170) !!??    CALL lbc_lnk( za1ct, 'I', -1. )
(  171) !!??    CALL lbc_lnk( za2ct, 'I', -1. )
(  172) 
(  173) 
(  174)       ! SOLUTION OF THE MOMENTUM EQUATION.
(  175)       !------------------------------------------
(  176)       !                                                   ! ==================== !
(  177)       DO iter = 1 , 2 * nbiter                            !    loop over iter    !
(  178)          !                                                ! ==================== !        
(  179)          zindu = MOD( iter , 2 )
(  180)          zusdtp = ( zindu * 2.0 + ( 1.0 - zindu ) * 1.0 )  * REAL( nbiter ) / rdt_ice
(  181) 
(  182)          ! Computation of free drift field for free slip boundary conditions.
(  183) 
(  184)            DO jj = k_j1, k_jpj-1
(  185)               DO ji = 1, jpim1
(  186)                  !- Rate of strain tensor.
(  187)                  zt11 =   akappa(ji,jj,1,1) * ( u_ice(ji+1,jj) + u_ice(ji+1,jj+1) - u_ice(ji,jj  ) - u_ice(ji  ,jj+1) )  &
(  188)                     &   + akappa(ji,jj,1,2) * ( v_ice(ji+1,jj) + v_ice(ji+1,jj+1) + v_ice(ji,jj  ) + v_ice(ji  ,jj+1) )
(  189)                  zt12 = - akappa(ji,jj,2,2) * ( u_ice(ji  ,jj) + u_ice(ji+1,jj  ) - u_ice(ji,jj+1) - u_ice(ji+1,jj+1) )  &
(  190)                     &   - akappa(ji,jj,2,1) * ( v_ice(ji  ,jj) + v_ice(ji+1,jj  ) + v_ice(ji,jj+1) + v_ice(ji+1,jj+1) )
(  191)                  zt22 = - akappa(ji,jj,2,2) * ( v_ice(ji  ,jj) + v_ice(ji+1,jj  ) - v_ice(ji,jj+1) - v_ice(ji+1,jj+1) )  &
(  192)                     &   + akappa(ji,jj,2,1) * ( u_ice(ji  ,jj) + u_ice(ji+1,jj  ) + u_ice(ji,jj+1) + u_ice(ji+1,jj+1) )
(  193)                  zt21 =   akappa(ji,jj,1,1) * ( v_ice(ji+1,jj) + v_ice(ji+1,jj+1) - v_ice(ji,jj  ) - v_ice(ji  ,jj+1) )  &
(  194)                     &   - akappa(ji,jj,1,2) * ( u_ice(ji+1,jj) + u_ice(ji+1,jj+1) + u_ice(ji,jj  ) + u_ice(ji  ,jj+1) )
(  195) 
(  196)                  !- Rate of strain tensor. 
(  197)                  zdgp = zt11 + zt22
(  198)                  zdgi = zt12 + zt21
(  199)                  ztrace2 = zdgp * zdgp 
(  200)                  zdeter  = zt11 * zt22 - 0.25 * zdgi * zdgi
(  201) 
(  202)                  !  Creep limit depends on the size of the grid.
(  203)                  zdelta = MAX( SQRT( ztrace2 + ( ztrace2 - 4.0 * zdeter ) * usecc2),  creepl)
(  204) 
(  205)                  !-  Computation of viscosities.
(  206)                  zviszeta(ji,jj) = MAX( zpresh(ji,jj) / zdelta, etamn )
(  207)                  zviseta (ji,jj) = zviszeta(ji,jj) * usecc2
(  208)               END DO
(  209)            END DO
(  210) !!??       CALL lbc_lnk( zviszeta, 'I', -1. )  ! or T point???   semble reellement inutile
(  211) !!??       CALL lbc_lnk( zviseta , 'I', -1. )
(  212) 
(  213) 
(  214)            !-  Determination of zc1u, zc2u, zc1v and zc2v.
(  215)            DO jj = k_j1+1, k_jpj-1
(  216)               DO ji = 2, jpim1
(  217)                  ze11   =  akappa(ji-1,jj-1,1,1)
(  218)                  ze12   = +akappa(ji-1,jj-1,2,2)
(  219)                  ze22   =  akappa(ji-1,jj-1,2,1)
(  220)                  ze21   = -akappa(ji-1,jj-1,1,2)






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 5

(  221)                  zvis11 = 2.0 * zviseta (ji-1,jj-1) + dm
(  222)                  zvis22 =       zviszeta(ji-1,jj-1) - zviseta(ji-1,jj-1)
(  223)                  zvis12 =       zviseta (ji-1,jj-1) + dm
(  224)                  zvis21 =       zviseta (ji-1,jj-1)
(  225) 
(  226)                  zdiag = zvis22 * ( ze11 + ze22 )
(  227)                  zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
(  228)                  zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
(  229)                  zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag
(  230)                  zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12
(  231) 
(  232)                  ze11   = -akappa(ji,jj-1,1,1)
(  233)                  ze12   = +akappa(ji,jj-1,2,2)
(  234)                  ze22   =  akappa(ji,jj-1,2,1)
(  235)                  ze21   = -akappa(ji,jj-1,1,2)
(  236)                  zvis11 = 2.0 * zviseta (ji,jj-1) + dm
(  237)                  zvis22 =       zviszeta(ji,jj-1) - zviseta(ji,jj-1)
(  238)                  zvis12 =       zviseta (ji,jj-1) + dm
(  239)                  zvis21 =       zviseta (ji,jj-1)
(  240) 
(  241)                  zdiag = zvis22 * ( ze11 + ze22 )
(  242)                  zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
(  243)                  zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
(  244)                  zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
(  245)                  zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12
(  246) 
(  247)                  ze11   =  akappa(ji-1,jj,1,1)
(  248)                  ze12   = -akappa(ji-1,jj,2,2)
(  249)                  ze22   =  akappa(ji-1,jj,2,1)
(  250)                  ze21   = -akappa(ji-1,jj,1,2)
(  251)                  zvis11 = 2.0 * zviseta (ji-1,jj) + dm
(  252)                  zvis22 =       zviszeta(ji-1,jj) - zviseta(ji-1,jj)
(  253)                  zvis12 =       zviseta (ji-1,jj) + dm
(  254)                  zvis21 =       zviseta (ji-1,jj)
(  255) 
(  256)                  zdiag = zvis22 * ( ze11 + ze22 ) 
(  257)                  zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
(  258)                  zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
(  259)                  zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
(  260)                  zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12
(  261) 
(  262)                  ze11   = -akappa(ji,jj,1,1)
(  263)                  ze12   = -akappa(ji,jj,2,2)
(  264)                  ze22   =  akappa(ji,jj,2,1)
(  265)                  ze21   = -akappa(ji,jj,1,2)
(  266)                  zvis11 = 2.0 * zviseta (ji,jj) + dm
(  267)                  zvis22 =       zviszeta(ji,jj) - zviseta(ji,jj)
(  268)                  zvis12 =       zviseta (ji,jj) + dm
(  269)                  zvis21 =       zviseta (ji,jj)
(  270) 
(  271)                  zdiag = zvis22 * ( ze11 + ze22 )
(  272)                  zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
(  273)                  zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
(  274)                  zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag 
(  275)                  zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12
(  276)               END DO
(  277)            END DO
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 6

(  279)            DO jj = k_j1+1, k_jpj-1
(  280)               DO ji = 2, jpim1
(  281)                  zc1u(ji,jj) =   &
(  282)                     + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)   &
(  283)                     - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)   &
(  284)                     - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)   &
(  285)                     + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)   &
(  286)                     + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)   &
(  287)                     + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)   &
(  288)                     - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)   &
(  289)                     - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
(  290)                  
(  291)                  zc2u(ji,jj) =   &
(  292)                     + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)   &
(  293)                     - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)   &
(  294)                     - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)   &
(  295)                     + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)   &
(  296)                     - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)   &
(  297)                     - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)   &
(  298)                     + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)   &
(  299)                     + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
(  300)              END DO
(  301)            END DO
(  302) 
(  303)            DO jj = k_j1+1, k_jpj-1
(  304)               DO ji = 2, jpim1
(  305)                  !  zc1v , zc2v.
(  306)                  ze11   =  akappa(ji-1,jj-1,1,2)
(  307)                  ze12   = -akappa(ji-1,jj-1,2,1)
(  308)                  ze22   = +akappa(ji-1,jj-1,2,2)
(  309)                  ze21   =  akappa(ji-1,jj-1,1,1)
(  310)                  zvis11 = 2.0 * zviseta (ji-1,jj-1) + dm
(  311)                  zvis22 =       zviszeta(ji-1,jj-1) - zviseta(ji-1,jj-1)
(  312)                  zvis12 =       zviseta (ji-1,jj-1) + dm
(  313)                  zvis21 =       zviseta (ji-1,jj-1)
(  314) 
(  315)                  zdiag = zvis22 * ( ze11 + ze22 )
(  316)                  zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
(  317)                  zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
(  318)                  zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag
(  319)                  zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12
(  320)  
(  321)                  ze11   =  akappa(ji,jj-1,1,2)
(  322)                  ze12   = -akappa(ji,jj-1,2,1)
(  323)                  ze22   = +akappa(ji,jj-1,2,2)
(  324)                  ze21   = -akappa(ji,jj-1,1,1)
(  325)                  zvis11 = 2.0 * zviseta (ji,jj-1) + dm
(  326)                  zvis22 =       zviszeta(ji,jj-1) - zviseta(ji,jj-1)
(  327)                  zvis12 =       zviseta (ji,jj-1) + dm
(  328)                  zvis21 =       zviseta (ji,jj-1)
(  329) 
(  330)                  zdiag = zvis22 * ( ze11 + ze22 )
(  331)                  zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
(  332)                  zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
(  333)                  zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
(  334)                  zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12
(  335) 
(  336)                  ze11   =  akappa(ji-1,jj,1,2)






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 7

(  337)                  ze12   = -akappa(ji-1,jj,2,1)
(  338)                  ze22   = -akappa(ji-1,jj,2,2)
(  339)                  ze21   =  akappa(ji-1,jj,1,1)
(  340)                  zvis11 = 2.0 * zviseta (ji-1,jj) + dm
(  341)                  zvis22 =       zviszeta(ji-1,jj) - zviseta(ji-1,jj)
(  342)                  zvis12 =       zviseta (ji-1,jj) + dm
(  343)                  zvis21 =       zviseta (ji-1,jj)
(  344) 
(  345)                  zdiag = zvis22 * ( ze11 + ze22 )
(  346)                  zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
(  347)                  zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
(  348)                  zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
(  349)                  zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12
(  350) 
(  351)                  ze11   =  akappa(ji,jj,1,2)
(  352)                  ze12   = -akappa(ji,jj,2,1)
(  353)                  ze22   = -akappa(ji,jj,2,2)
(  354)                  ze21   = -akappa(ji,jj,1,1)
(  355)                  zvis11 = 2.0 * zviseta (ji,jj) + dm
(  356)                  zvis22 =       zviszeta(ji,jj) - zviseta(ji,jj)
(  357)                  zvis12 =       zviseta (ji,jj) + dm
(  358)                  zvis21 =       zviseta (ji,jj)
(  359) 
(  360)                  zdiag = zvis22 * ( ze11 + ze22 )
(  361)                  zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
(  362)                  zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
(  363)                  zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag
(  364)                  zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12
(  365) 
(  366)               END DO
(  367)            END DO
(  368) 
(  369)            DO jj = k_j1+1, k_jpj-1
(  370)               DO ji = 2, jpim1
(  371)                  zc1v(ji,jj) =   &
(  372)                     + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)   &
(  373)                     - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)   &
(  374)                     - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)   &
(  375)                     + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)   &
(  376)                     + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)   &
(  377)                     + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)   &
(  378)                     - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)   &
(  379)                     - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
(  380)                  zc2v(ji,jj) =   &
(  381)                     + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)   &
(  382)                     - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)   &
(  383)                     - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)   &
(  384)                     + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)   &
(  385)                     - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)   &
(  386)                     - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)   &
(  387)                     + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)   &
(  388)                     + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
(  389)               END DO
(  390)            END DO
(  391) 
(  392)          ! Relaxation.
(  393)            
(  394) iflag:   DO jter = 1 , nbitdr






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 8

(  395) 
(  396)             !  Store previous drift field.   
(  397)             DO jj = k_j1, k_jpj-1
(  398)                zu_ice(:,jj) = u_ice(:,jj)
(  399)                zv_ice(:,jj) = v_ice(:,jj)
(  400)             END DO
(  401) 
(  402)             DO jj = k_j1+1, k_jpj-1
(  403)                zsang  = SIGN( 1.e0, fcor(1,jj) ) * sangvg   ! only the sinus changes its sign with the hemisphere
(  404)                DO ji = 2, jpim1
(  405)                  zur     = u_ice(ji,jj) - u_oce(ji,jj)
(  406)                  zvr     = v_ice(ji,jj) - v_oce(ji,jj)
(  407)                  zmod    = SQRT( zur * zur + zvr * zvr) * ( 1.0 - zfrld(ji,jj) )
(  408)                  za      = rhoco * zmod
(  409)                  zac     = za * cangvg
(  410)                   zmpzas  = alpha * zcorl(ji,jj) + za * zsang
(  411)                   zmassdt = zusdtp * zmass(ji,jj)
(  412)                   zcorlal = ( 1.0 - alpha ) * zcorl(ji,jj)
(  413) 
(  414)                   za1(ji,jj) =  zmassdt * zu0(ji,jj) + zcorlal * zv0(ji,jj) + za1ct(ji,jj)   &
(  415)                      &        + za * ( cangvg * u_oce(ji,jj) - zsang * v_oce(ji,jj) )
(  416) 
(  417)                   za2(ji,jj) =  zmassdt * zv0(ji,jj) - zcorlal * zu0(ji,jj) + za2ct(ji,jj)   &
(  418)                      &        + za * ( cangvg * v_oce(ji,jj) + zsang * u_oce(ji,jj) )
(  419) 
(  420)                   zb1(ji,jj)  = zmassdt + zac - zc1u(ji,jj)
(  421)                   zb2(ji,jj)  = zmpzas  - zc2u(ji,jj)
(  422)                   zc1(ji,jj)  = zmpzas  + zc1v(ji,jj)
(  423)                   zc2(ji,jj)  = zmassdt + zac  - zc2v(ji,jj) 
(  424)                   zdeter      = zc1(ji,jj) * zb2(ji,jj) + zc2(ji,jj) * zb1(ji,jj)
(  425)                   zden(ji,jj) = SIGN( rone, zdeter) / MAX( epsd , ABS( zdeter ) )
(  426)                END DO
(  427)             END DO
(  428) 
(  429)             ! The computation of ice interaction term is splitted into two parts
(  430)             !-------------------------------------------------------------------------
(  431) 
(  432)             ! Terms that do not involve already up-dated velocities.
(  433)          
(  434)             DO jj = k_j1+1, k_jpj-1
(  435)                DO ji = 2, jpim1
(  436)                   iim1 = ji
(  437)                   ijm1 = jj - 1
(  438)                   iip1 = ji + 1
(  439)                   ijp1 = jj
(  440)                   ze11 =   akappa(iim1,ijm1,1,1) * u_ice(iip1,ijp1) + akappa(iim1,ijm1,1,2) * v_ice(iip1,ijp1)
(  441)                   ze12 = + akappa(iim1,ijm1,2,2) * u_ice(iip1,ijp1) - akappa(iim1,ijm1,2,1) * v_ice(iip1,ijp1)
(  442)                   ze22 = + akappa(iim1,ijm1,2,2) * v_ice(iip1,ijp1) + akappa(iim1,ijm1,2,1) * u_ice(iip1,ijp1)
(  443)                   ze21 =   akappa(iim1,ijm1,1,1) * v_ice(iip1,ijp1) - akappa(iim1,ijm1,1,2) * u_ice(iip1,ijp1)
(  444)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  445)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  446)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  447)                   zvis21 =       zviseta (iim1,ijm1)
(  448)                   zdiag = zvis22 * ( ze11 + ze22 )
(  449)                   zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
(  450)                   zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
(  451)                   zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
(  452)                   zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 9

(  453) 
(  454) 
(  455)                   iim1 = ji - 1
(  456)                   ijm1 = jj
(  457)                   iip1 = ji
(  458)                   ijp1 = jj + 1                   
(  459)                   ze11 =   akappa(iim1,ijm1,1,1) * ( u_ice(iip1,ijp1) - u_ice(iim1,ijp1) )   &
(  460)                      &   + akappa(iim1,ijm1,1,2) * ( v_ice(iip1,ijp1) + v_ice(iim1,ijp1) )
(  461)                   ze12 = + akappa(iim1,ijm1,2,2) * ( u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )   &
(  462)                      &   - akappa(iim1,ijm1,2,1) * ( v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )
(  463)                   ze22 = + akappa(iim1,ijm1,2,2) * ( v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )   &
(  464)                      &   + akappa(iim1,ijm1,2,1) * ( u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )
(  465)                   ze21 =   akappa(iim1,ijm1,1,1) * ( v_ice(iip1,ijp1) - v_ice(iim1,ijp1) )   &
(  466)                      &   - akappa(iim1,ijm1,1,2) * ( u_ice(iip1,ijp1) + u_ice(iim1,ijp1) )
(  467)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  468)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  469)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  470)                   zvis21 =       zviseta (iim1,ijm1)
(  471)                   zdiag = zvis22 * ( ze11 + ze22 )
(  472)                   zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
(  473)                   zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
(  474)                   zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
(  475)                   zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12
(  476) 
(  477)                   iim1 = ji
(  478)                   ijm1 = jj
(  479)                   iip1 = ji + 1
(  480)                   ijp1 = jj + 1
(  481)                   ze11 =   akappa(iim1,ijm1,1,1) * ( u_ice(iip1,ijm1) + u_ice(iip1,ijp1) - u_ice(iim1,ijp1) )   &
(  482)                      &   + akappa(iim1,ijm1,1,2) * ( v_ice(iip1,ijm1) + v_ice(iip1,ijp1) + v_ice(iim1,ijp1) )
(  483)                   ze12 = - akappa(iim1,ijm1,2,2) * ( u_ice(iip1,ijm1) - u_ice(iim1,ijp1) - u_ice(iip1,ijp1) )   &
(  484)                      &   - akappa(iim1,ijm1,2,1) * ( v_ice(iip1,ijm1) + v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )
(  485)                   ze22 = - akappa(iim1,ijm1,2,2) * ( v_ice(iip1,ijm1) - v_ice(iim1,ijp1) - v_ice(iip1,ijp1) )   &
(  486)                      &   + akappa(iim1,ijm1,2,1) * ( u_ice(iip1,ijm1) + u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )
(  487)                   ze21 =   akappa(iim1,ijm1,1,1) * ( v_ice(iip1,ijm1) + v_ice(iip1,ijp1) - v_ice(iim1,ijp1) )   &
(  488)                      &   - akappa(iim1,ijm1,1,2) * ( u_ice(iip1,ijm1) + u_ice(iip1,ijp1) + u_ice(iim1,ijp1) ) 
(  489)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  490)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  491)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  492)                   zvis21 =       zviseta (iim1,ijm1)
(  493) 
(  494)                   zdiag = zvis22 * ( ze11 + ze22 )
(  495)                   zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
(  496)                   zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
(  497)                   zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag
(  498)                   zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12
(  499) 
(  500)                END DO
(  501)             END DO
(  502) 
(  503)             ! Terms involving already up-dated velocities.
(  504)             !-Using the arrays zu_ice and zv_ice in the computation of the terms ze leads to JACOBI's method; 
(  505)             ! Using arrays u and v in the computation of the terms ze leads to GAUSS-SEIDEL method.
(  506)              
(  507)             DO jj = k_j1+1, k_jpj-1
(  508)                DO ji = 2, jpim1
(  509)                   iim1 = ji - 1
(  510)                   ijm1 = jj - 1






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 10

(  511)                   iip1 = ji
(  512)                   ijp1 = jj
(  513)                   ze11 =   akappa(iim1,ijm1,1,1) * ( zu_ice(iip1,ijm1) - zu_ice(iim1,ijm1) - zu_ice(iim1,ijp1) )   &
(  514)                      &   + akappa(iim1,ijm1,1,2) * ( zv_ice(iip1,ijm1) + zv_ice(iim1,ijm1) + zv_ice(iim1,ijp1) )
(  515)                   ze12 = - akappa(iim1,ijm1,2,2) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) - zu_ice(iim1,ijp1) )   &
(  516)                      &   - akappa(iim1,ijm1,2,1) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) + zv_ice(iim1,ijp1) )
(  517)                   ze22 = - akappa(iim1,ijm1,2,2) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) - zv_ice(iim1,ijp1) )   &
(  518)                      &   + akappa(iim1,ijm1,2,1) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) + zu_ice(iim1,ijp1) )
(  519)                   ze21 =   akappa(iim1,ijm1,1,1) * ( zv_ice(iip1,ijm1) - zv_ice(iim1,ijm1) - zv_ice(iim1,ijp1) )   &
(  520)                      &  -  akappa(iim1,ijm1,1,2) * ( zu_ice(iip1,ijm1) + zu_ice(iim1,ijm1) + zu_ice(iim1,ijp1) )
(  521)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  522)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  523)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  524)                   zvis21 =       zviseta (iim1,ijm1)
(  525) 
(  526)                   zdiag = zvis22 * ( ze11 + ze22 )
(  527)                   zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
(  528)                   zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
(  529)                   zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag
(  530)                   zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12
(  531) 
(  532) #if defined key_agrif
(  533)              END DO
(  534)           END DO
(  535) 
(  536)           DO jj = k_j1+1, k_jpj-1
(  537)              DO ji = 2, jpim1
(  538) #endif
(  539) 
(  540)                   iim1 = ji
(  541)                   ijm1 = jj - 1
(  542)                   iip1 = ji + 1
(  543)                   ze11 =   akappa(iim1,ijm1,1,1) * ( zu_ice(iip1,ijm1) - zu_ice(iim1,ijm1) )   &
(  544)                      &   + akappa(iim1,ijm1,1,2) * ( zv_ice(iip1,ijm1) + zv_ice(iim1,ijm1) )
(  545)                   ze12 = - akappa(iim1,ijm1,2,2) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) )   &
(  546)                      &   - akappa(iim1,ijm1,2,1) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) )
(  547)                   ze22 = - akappa(iim1,ijm1,2,2) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) )   &
(  548)                      &   + akappa(iim1,ijm1,2,1) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) )
(  549)                   ze21 =   akappa(iim1,ijm1,1,1) * ( zv_ice(iip1,ijm1) - zv_ice(iim1,ijm1) )   &
(  550)                      &   - akappa(iim1,ijm1,1,2) * ( zu_ice(iip1,ijm1) + zu_ice(iim1,ijm1) )
(  551)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  552)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  553)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  554)                   zvis21 =       zviseta (iim1,ijm1)
(  555) 
(  556)                   zdiag = zvis22 * ( ze11 + ze22 )
(  557)                   zs11(ji,jj,2,1) =  zs11(ji,jj,2,1) + zvis11 * ze11 + zdiag
(  558)                   zs12(ji,jj,2,1) =  zs12(ji,jj,2,1) + zvis12 * ze12 + zvis21 * ze21
(  559)                   zs22(ji,jj,2,1) =  zs22(ji,jj,2,1) + zvis11 * ze22 + zdiag
(  560)                   zs21(ji,jj,2,1) =  zs21(ji,jj,2,1) + zvis12 * ze21 + zvis21 * ze12
(  561) 
(  562) 
(  563)                   iim1 = ji - 1
(  564)                   ijm1 = jj 
(  565)                   ze11 = - akappa(iim1,ijm1,1,1) * zu_ice(iim1,ijm1) + akappa(iim1,ijm1,1,2) * zv_ice(iim1,ijm1)
(  566)                   ze12 = - akappa(iim1,ijm1,2,2) * zu_ice(iim1,ijm1) - akappa(iim1,ijm1,2,1) * zv_ice(iim1,ijm1)
(  567)                   ze22 = - akappa(iim1,ijm1,2,2) * zv_ice(iim1,ijm1) + akappa(iim1,ijm1,2,1) * zu_ice(iim1,ijm1)
(  568)                   ze21 = - akappa(iim1,ijm1,1,1) * zv_ice(iim1,ijm1) - akappa(iim1,ijm1,1,2) * zu_ice(iim1,ijm1)






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 11

(  569)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  570)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  571)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  572)                   zvis21 =       zviseta (iim1,ijm1)
(  573) 
(  574)                   zdiag = zvis22 * ( ze11 + ze22 )
(  575)                   zs11(ji,jj,1,2) =  zs11(ji,jj,1,2) + zvis11 * ze11 + zdiag 
(  576)                   zs12(ji,jj,1,2) =  zs12(ji,jj,1,2) + zvis12 * ze12 + zvis21 * ze21
(  577)                   zs22(ji,jj,1,2) =  zs22(ji,jj,1,2) + zvis11 * ze22 + zdiag
(  578)                   zs21(ji,jj,1,2) =  zs21(ji,jj,1,2) + zvis12 * ze21 + zvis21 * ze12
(  579) 
(  580) #if defined key_agrif
(  581)              END DO
(  582)           END DO
(  583) 
(  584)           DO jj = k_j1+1, k_jpj-1
(  585)              DO ji = 2, jpim1
(  586) #endif
(  587)                   zd1(ji,jj) =   &
(  588)                      + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)  &
(  589)                      - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)  &
(  590)                      - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)  &
(  591)                      + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)  &
(  592)                      + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)  &
(  593)                      + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)  &
(  594)                      - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)  &
(  595)                      - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
(  596)                   zd2(ji,jj) =   &
(  597)                      + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)  &
(  598)                      - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)  &
(  599)                      - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)  &
(  600)                      + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)  &
(  601)                      - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)  &
(  602)                      - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)  &
(  603)                      + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)  &
(  604)                      + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
(  605)                END DO
(  606)             END DO
(  607) 
(  608)             DO jj = k_j1+1, k_jpj-1
(  609)                DO ji = 2, jpim1
(  610)                   zunw = (  ( za1(ji,jj) + zd1(ji,jj) ) * zc2(ji,jj)        &
(  611)                      &    + ( za2(ji,jj) + zd2(ji,jj) ) * zc1(ji,jj) ) * zden(ji,jj)
(  612) 
(  613)                   zvnw = (  ( za2(ji,jj) + zd2(ji,jj) ) * zb1(ji,jj)        &
(  614)                      &    - ( za1(ji,jj) + zd1(ji,jj) ) * zb2(ji,jj) ) * zden(ji,jj)
(  615) 
(  616)                   zmask = ( 1.0 - MAX( rzero, SIGN( rone , 1.0 - zmass(ji,jj) ) ) ) * tmu(ji,jj)
(  617) 
(  618)                   u_ice(ji,jj) = ( u_ice(ji,jj) + om * ( zunw - u_ice(ji,jj) ) * tmu(ji,jj) ) * zmask
(  619)                   v_ice(ji,jj) = ( v_ice(ji,jj) + om * ( zvnw - v_ice(ji,jj) ) * tmu(ji,jj) ) * zmask
(  620)                END DO
(  621)             END DO
(  622) 
(  623)             CALL lbc_lnk( u_ice, 'I', -1. )
(  624)             CALL lbc_lnk( v_ice, 'I', -1. )
(  625) 
(  626)             !---  5.2.5.4. Convergence test.






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 12

(  627)             DO jj = k_j1+1 , k_jpj-1
(  628)                zresr(:,jj) = MAX( ABS( u_ice(:,jj) - zu_ice(:,jj) ) , ABS( v_ice(:,jj) - zv_ice(:,jj) ) )
(  629)             END DO
(  630)             zresm = MAXVAL( zresr( 1:jpi , k_j1+1:k_jpj-1 ) )
(  631)             IF( lk_mpp )   CALL mpp_max( zresm )   ! max over the global domain
(  632) 
(  633)             IF ( zresm <= resl) EXIT iflag
(  634) 
(  635)          END DO iflag
(  636) 
(  637)          zindu1 = 1.0 - zindu
(  638)          DO jj = k_j1 , k_jpj-1
(  639)             zu0(:,jj) = zindu * zu0(:,jj) + zindu1 * u_ice(:,jj)
(  640)             zv0(:,jj) = zindu * zv0(:,jj) + zindu1 * v_ice(:,jj)
(  641)          END DO
(  642)       !                                                   ! ==================== !
(  643)       END DO                                              !  end loop over iter  !
(  644)       !                                                   ! ==================== !
(  645) 
(  646)       IF(ln_ctl) THEN
(  647)          WRITE(charout,FMT="('lim_rhg  : res =',D23.16, ' iter =',I4)") zresm, jter
(  648)          CALL prt_ctl_info(charout)
(  649)          CALL prt_ctl(tab2d_1=u_ice, clinfo1=' lim_rhg  : u_ice :', tab2d_2=v_ice, clinfo2=' v_ice :')
(  650)       ENDIF
(  651) 
(  652)    END SUBROUTINE lim_rhg
(  653) 
(  654) #else
(  655)    !!----------------------------------------------------------------------
(  656)    !!   Default option          Dummy module           NO LIM sea-ice model
(  657)    !!----------------------------------------------------------------------
(  658)    USE in_out_manager
(  659) CONTAINS
(  660)    SUBROUTINE lim_rhg( k1 , k2 )         ! Dummy routine
(  661)       if(lwp) WRITE(numout,*) 'lim_rhg: You should not have seen this print! error?', k1, k2
(  662)    END SUBROUTINE lim_rhg
(  663) #endif
(  664) 
(  665)    !!==============================================================================
(  666) END MODULE limrhg
























PGF90 (Version     10.5)          09/23/2010  09:49:54      page 13

(    1) # 1 "limrhg.F90"
(    1) MODULE limrhg
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  limrhg  ***
(    4)    !!   Ice rheology :  performs sea ice rheology
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_rhg   : computes ice velocities
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE phycst
(   14)    USE par_oce
(   15)    USE ice_oce         ! ice variables
(   16)    USE dom_ice
(   17)    USE ice
(   18)    USE lbclnk
(   19)    USE lib_mpp
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE prtctl          ! Print control
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC lim_rhg  ! routine called by lim_dyn
(   28) 
(   29)    !! * Module variables
(   30)    REAL(wp)  ::           &  ! constant values
(   31)       rzero   = 0.e0   ,  &
(   32)       rone    = 1.e0
(   33)    !!----------------------------------------------------------------------
(   34)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   35)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limrhg.F90,v 1.7 2006/03/10 10:35:42 opalod Exp $ 
(   36)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   37)    !!----------------------------------------------------------------------
(   38) 
(   39) CONTAINS
(   40) 
(   41)    SUBROUTINE lim_rhg( k_j1, k_jpj )
(   42)       !!-------------------------------------------------------------------
(   43)       !!                 ***  SUBROUTINR lim_rhg  ***
(   44)       !!
(   45)       !! ** purpose :   determines the velocity field of sea ice by using
(   46)       !!  atmospheric (wind stress) and oceanic (water stress and surface
(   47)       !!  tilt) forcings. Ice-ice interaction is described by a non-linear
(   48)       !!  viscous-plastic law including shear strength and a bulk rheology.
(   49)       !!
(   50)       !! ** Action  : - compute u_ice, v_ice the sea-ice velocity
(   51)       !!
(   52)       !! History :
(   53)       !!   0.0  !  93-12  (M.A. Morales Maqueda.)  Original code
(   54)       !!   1.0  !  94-12  (H. Goosse) 
(   55)       !!   2.0  !  03-12  (C. Ethe, G. Madec)  F90, mpp
(   56)       !!-------------------------------------------------------------------
(   57)       ! * Arguments






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 14

(   58)       INTEGER, INTENT(in) ::   k_j1 ,  &  ! southern j-index for ice computation
(   59)          &                     k_jpj      ! northern j-index for ice computation
(   60) 
(   61)       ! * Local variables
(   62)       INTEGER ::   ji, jj              ! dummy loop indices
(   63) 
(   64)       INTEGER  :: &
(   65)          iim1, ijm1, iip1 , ijp1   , & ! temporary integers
(   66)          iter, jter                    !    "          "
(   67) 
(   68)       CHARACTER (len=50) :: charout
(   69) 
(   70)       REAL(wp) :: &
(   71)          ze11  , ze12  , ze22  , ze21  ,   &  ! temporary scalars
(   72)          zt11  , zt12  , zt21  , zt22  ,   &  !    "         "
(   73)          zvis11, zvis21, zvis12, zvis22,   &  !    "         "
(   74)          zgphsx, ztagnx, zusw  ,           &  !    "         "
(   75)          zgphsy, ztagny                       !    "         "
(   76)       REAL(wp) :: &
(   77)          zresm, zunw, zvnw, zur, zvr, zmod, za, zac, &
(   78)          zmpzas, zstms, zindu, zindu1, zusdtp, zmassdt, zcorlal,  &
(   79)          ztrace2, zdeter, zdelta, zsang, zmask, zdgp, zdgi, zdiag
(   80)       REAL(wp),DIMENSION(jpi,jpj) ::   &
(   81)          zpresh, zfrld, zmass, zcorl,     &
(   82)          zu0, zv0, zviszeta, zviseta,     &
(   83)          zc1u, zc1v, zc2u, zc2v, za1ct, za2ct, za1, za2, zb1, zb2,  &
(   84)          zc1, zc2, zd1, zd2, zden, zu_ice, zv_ice, zresr
(   85)       REAL(wp),DIMENSION(jpi,jpj,2,2) :: &
(   86)          zs11, zs12, zs22, zs21
(   87)       !!-------------------------------------------------------------------
(   88)       
(   89)       !  Store initial velocities
(   90)       !  ------------------------
(   91)       zu0(:,:) = u_ice(:,:)
(   92)       zv0(:,:) = v_ice(:,:)
(   93) 
(   94)       ! Ice mass, ice strength, and wind stress at the center            |
(   95)       ! of the grid squares.                                             |
(   96)       !-------------------------------------------------------------------
(   97) 
(   98)       DO jj = k_j1 , k_jpj-1
(   99)          DO ji = 1 , jpi
(  100)             za1(ji,jj)    = tms(ji,jj) * ( rhosn * hsnm(ji,jj) + rhoic * hicm(ji,jj) )
(  101)             zpresh(ji,jj) = tms(ji,jj) *  pstarh * hicm(ji,jj) * EXP( -c_rhg * frld(ji,jj) )
(  102) # 106
(  106)             zb1(ji,jj)    = tms(ji,jj) * ( 1.0 - frld(ji,jj) )
(  107)             zb2(ji,jj)    = tms(ji,jj) * ( 1.0 - frld(ji,jj) )
(  108) # 109
(  109)          END DO
(  110)       END DO
(  111) 
(  112) 
(  113)       !---------------------------------------------------------------------------
(  114)       !  Wind stress, coriolis and mass terms at the corners of the grid squares |
(  115)       !  Gradient of ice strenght.                                               |
(  116)       !---------------------------------------------------------------------------
(  117)          
(  118)       DO jj = k_j1+1, k_jpj-1






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 15

(  119)          DO ji = 2, jpi
(  120)             zstms = tms(ji,jj  ) * wght(ji,jj,2,2) + tms(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  121)                &  + tms(ji,jj-1) * wght(ji,jj,2,1) + tms(ji-1,jj-1) * wght(ji,jj,1,1)
(  122)             zusw  = 1.0 / MAX( zstms, epsd )
(  123) 
(  124)             zt11 = tms(ji  ,jj  ) * frld(ji  ,jj  ) 
(  125)             zt12 = tms(ji-1,jj  ) * frld(ji-1,jj  ) 
(  126)             zt21 = tms(ji  ,jj-1) * frld(ji  ,jj-1) 
(  127)             zt22 = tms(ji-1,jj-1) * frld(ji-1,jj-1)
(  128) 
(  129)             ! Leads area.
(  130)             zfrld(ji,jj) =  (  zt11 * wght(ji,jj,2,2) + zt12 * wght(ji,jj,1,2)   &
(  131)                &             + zt21 * wght(ji,jj,2,1) + zt22 * wght(ji,jj,1,1) ) * zusw
(  132) 
(  133)             ! Mass and coriolis coeff.
(  134)             zmass(ji,jj) = ( za1(ji,jj  ) * wght(ji,jj,2,2) + za1(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  135)                &           + za1(ji,jj-1) * wght(ji,jj,2,1) + za1(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw
(  136)             zcorl(ji,jj) = zmass(ji,jj) * fcor(ji,jj)
(  137) 
(  138)             ! Wind stress.
(  139) # 145
(  145)             ztagnx = ( zb1(ji,jj  ) * wght(ji,jj,2,2) + zb1(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  146)                &     + zb1(ji,jj-1) * wght(ji,jj,2,1) + zb1(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw * gtaux(ji,jj)
(  147)             ztagny = ( zb2(ji,jj  ) * wght(ji,jj,2,2) + zb2(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  148)                &     + zb2(ji,jj-1) * wght(ji,jj,2,1) + zb2(ji-1,jj-1) * wght(ji,jj,1,1) ) * zusw * gtauy(ji,jj)
(  149) 
(  150) # 151
(  151)             ! Gradient of ice strength
(  152)             zgphsx =   ( alambd(ji,jj,2,2,2,1) - alambd(ji,jj,2,1,2,1) ) * zpresh(ji  ,jj-1)   &
(  153)                &     + ( alambd(ji,jj,2,2,2,2) - alambd(ji,jj,2,1,2,2) ) * zpresh(ji  ,jj  )   &
(  154)                &     - ( alambd(ji,jj,2,2,1,1) + alambd(ji,jj,2,1,1,1) ) * zpresh(ji-1,jj-1)   &
(  155)                &     - ( alambd(ji,jj,2,2,1,2) + alambd(ji,jj,2,1,1,2) ) * zpresh(ji-1,jj  )
(  156) 
(  157)             zgphsy = - ( alambd(ji,jj,1,1,2,1) + alambd(ji,jj,1,2,2,1) ) * zpresh(ji  ,jj-1)   &
(  158)                &     - ( alambd(ji,jj,1,1,1,1) + alambd(ji,jj,1,2,1,1) ) * zpresh(ji-1,jj-1)   &
(  159)                &     + ( alambd(ji,jj,1,1,2,2) - alambd(ji,jj,1,2,2,2) ) * zpresh(ji  ,jj  )   &
(  160)                &     + ( alambd(ji,jj,1,1,1,2) - alambd(ji,jj,1,2,1,2) ) * zpresh(ji-1,jj  )
(  161) 
(  162)             ! Computation of the velocity field taking into account the ice-ice interaction.                                 
(  163)             ! Terms that are independent of the velocity field.
(  164)             za1ct(ji,jj) = ztagnx - zcorl(ji,jj) * v_oce(ji,jj) - zgphsx
(  165)             za2ct(ji,jj) = ztagny + zcorl(ji,jj) * u_oce(ji,jj) - zgphsy
(  166)          END DO
(  167)       END DO
(  168) 
(  169) !! inutile!!
(  170) !!??    CALL lbc_lnk( za1ct, 'I', -1. )
(  171) !!??    CALL lbc_lnk( za2ct, 'I', -1. )
(  172) 
(  173) 
(  174)       ! SOLUTION OF THE MOMENTUM EQUATION.
(  175)       !------------------------------------------
(  176)       !                                                   ! ==================== !
(  177)       DO iter = 1 , 2 * nbiter                            !    loop over iter    !
(  178)          !                                                ! ==================== !        
(  179)          zindu = MOD( iter , 2 )
(  180)          zusdtp = ( zindu * 2.0 + ( 1.0 - zindu ) * 1.0 )  * REAL( nbiter ) / rdt_ice
(  181) 






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 16

(  182)          ! Computation of free drift field for free slip boundary conditions.
(  183) 
(  184)            DO jj = k_j1, k_jpj-1
(  185)               DO ji = 1, jpim1
(  186)                  !- Rate of strain tensor.
(  187)                  zt11 =   akappa(ji,jj,1,1) * ( u_ice(ji+1,jj) + u_ice(ji+1,jj+1) - u_ice(ji,jj  ) - u_ice(ji  ,jj+1) )  &
(  188)                     &   + akappa(ji,jj,1,2) * ( v_ice(ji+1,jj) + v_ice(ji+1,jj+1) + v_ice(ji,jj  ) + v_ice(ji  ,jj+1) )
(  189)                  zt12 = - akappa(ji,jj,2,2) * ( u_ice(ji  ,jj) + u_ice(ji+1,jj  ) - u_ice(ji,jj+1) - u_ice(ji+1,jj+1) )  &
(  190)                     &   - akappa(ji,jj,2,1) * ( v_ice(ji  ,jj) + v_ice(ji+1,jj  ) + v_ice(ji,jj+1) + v_ice(ji+1,jj+1) )
(  191)                  zt22 = - akappa(ji,jj,2,2) * ( v_ice(ji  ,jj) + v_ice(ji+1,jj  ) - v_ice(ji,jj+1) - v_ice(ji+1,jj+1) )  &
(  192)                     &   + akappa(ji,jj,2,1) * ( u_ice(ji  ,jj) + u_ice(ji+1,jj  ) + u_ice(ji,jj+1) + u_ice(ji+1,jj+1) )
(  193)                  zt21 =   akappa(ji,jj,1,1) * ( v_ice(ji+1,jj) + v_ice(ji+1,jj+1) - v_ice(ji,jj  ) - v_ice(ji  ,jj+1) )  &
(  194)                     &   - akappa(ji,jj,1,2) * ( u_ice(ji+1,jj) + u_ice(ji+1,jj+1) + u_ice(ji,jj  ) + u_ice(ji  ,jj+1) )
(  195) 
(  196)                  !- Rate of strain tensor. 
(  197)                  zdgp = zt11 + zt22
(  198)                  zdgi = zt12 + zt21
(  199)                  ztrace2 = zdgp * zdgp 
(  200)                  zdeter  = zt11 * zt22 - 0.25 * zdgi * zdgi
(  201) 
(  202)                  !  Creep limit depends on the size of the grid.
(  203)                  zdelta = MAX( SQRT( ztrace2 + ( ztrace2 - 4.0 * zdeter ) * usecc2),  creepl)
(  204) 
(  205)                  !-  Computation of viscosities.
(  206)                  zviszeta(ji,jj) = MAX( zpresh(ji,jj) / zdelta, etamn )
(  207)                  zviseta (ji,jj) = zviszeta(ji,jj) * usecc2
(  208)               END DO
(  209)            END DO
(  210) !!??       CALL lbc_lnk( zviszeta, 'I', -1. )  ! or T point???   semble reellement inutile
(  211) !!??       CALL lbc_lnk( zviseta , 'I', -1. )
(  212) 
(  213) 
(  214)            !-  Determination of zc1u, zc2u, zc1v and zc2v.
(  215)            DO jj = k_j1+1, k_jpj-1
(  216)               DO ji = 2, jpim1
(  217)                  ze11   =  akappa(ji-1,jj-1,1,1)
(  218)                  ze12   = +akappa(ji-1,jj-1,2,2)
(  219)                  ze22   =  akappa(ji-1,jj-1,2,1)
(  220)                  ze21   = -akappa(ji-1,jj-1,1,2)
(  221)                  zvis11 = 2.0 * zviseta (ji-1,jj-1) + dm
(  222)                  zvis22 =       zviszeta(ji-1,jj-1) - zviseta(ji-1,jj-1)
(  223)                  zvis12 =       zviseta (ji-1,jj-1) + dm
(  224)                  zvis21 =       zviseta (ji-1,jj-1)
(  225) 
(  226)                  zdiag = zvis22 * ( ze11 + ze22 )
(  227)                  zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
(  228)                  zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
(  229)                  zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag
(  230)                  zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12
(  231) 
(  232)                  ze11   = -akappa(ji,jj-1,1,1)
(  233)                  ze12   = +akappa(ji,jj-1,2,2)
(  234)                  ze22   =  akappa(ji,jj-1,2,1)
(  235)                  ze21   = -akappa(ji,jj-1,1,2)
(  236)                  zvis11 = 2.0 * zviseta (ji,jj-1) + dm
(  237)                  zvis22 =       zviszeta(ji,jj-1) - zviseta(ji,jj-1)
(  238)                  zvis12 =       zviseta (ji,jj-1) + dm
(  239)                  zvis21 =       zviseta (ji,jj-1)






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 17

(  240) 
(  241)                  zdiag = zvis22 * ( ze11 + ze22 )
(  242)                  zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
(  243)                  zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
(  244)                  zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
(  245)                  zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12
(  246) 
(  247)                  ze11   =  akappa(ji-1,jj,1,1)
(  248)                  ze12   = -akappa(ji-1,jj,2,2)
(  249)                  ze22   =  akappa(ji-1,jj,2,1)
(  250)                  ze21   = -akappa(ji-1,jj,1,2)
(  251)                  zvis11 = 2.0 * zviseta (ji-1,jj) + dm
(  252)                  zvis22 =       zviszeta(ji-1,jj) - zviseta(ji-1,jj)
(  253)                  zvis12 =       zviseta (ji-1,jj) + dm
(  254)                  zvis21 =       zviseta (ji-1,jj)
(  255) 
(  256)                  zdiag = zvis22 * ( ze11 + ze22 ) 
(  257)                  zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
(  258)                  zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
(  259)                  zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
(  260)                  zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12
(  261) 
(  262)                  ze11   = -akappa(ji,jj,1,1)
(  263)                  ze12   = -akappa(ji,jj,2,2)
(  264)                  ze22   =  akappa(ji,jj,2,1)
(  265)                  ze21   = -akappa(ji,jj,1,2)
(  266)                  zvis11 = 2.0 * zviseta (ji,jj) + dm
(  267)                  zvis22 =       zviszeta(ji,jj) - zviseta(ji,jj)
(  268)                  zvis12 =       zviseta (ji,jj) + dm
(  269)                  zvis21 =       zviseta (ji,jj)
(  270) 
(  271)                  zdiag = zvis22 * ( ze11 + ze22 )
(  272)                  zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
(  273)                  zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
(  274)                  zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag 
(  275)                  zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12
(  276)               END DO
(  277)            END DO
(  278) 
(  279)            DO jj = k_j1+1, k_jpj-1
(  280)               DO ji = 2, jpim1
(  281)                  zc1u(ji,jj) =   &
(  282)                     + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)   &
(  283)                     - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)   &
(  284)                     - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)   &
(  285)                     + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)   &
(  286)                     + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)   &
(  287)                     + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)   &
(  288)                     - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)   &
(  289)                     - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
(  290)                  
(  291)                  zc2u(ji,jj) =   &
(  292)                     + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)   &
(  293)                     - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)   &
(  294)                     - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)   &
(  295)                     + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)   &
(  296)                     - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)   &
(  297)                     - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)   &






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 18

(  298)                     + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)   &
(  299)                     + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
(  300)              END DO
(  301)            END DO
(  302) 
(  303)            DO jj = k_j1+1, k_jpj-1
(  304)               DO ji = 2, jpim1
(  305)                  !  zc1v , zc2v.
(  306)                  ze11   =  akappa(ji-1,jj-1,1,2)
(  307)                  ze12   = -akappa(ji-1,jj-1,2,1)
(  308)                  ze22   = +akappa(ji-1,jj-1,2,2)
(  309)                  ze21   =  akappa(ji-1,jj-1,1,1)
(  310)                  zvis11 = 2.0 * zviseta (ji-1,jj-1) + dm
(  311)                  zvis22 =       zviszeta(ji-1,jj-1) - zviseta(ji-1,jj-1)
(  312)                  zvis12 =       zviseta (ji-1,jj-1) + dm
(  313)                  zvis21 =       zviseta (ji-1,jj-1)
(  314) 
(  315)                  zdiag = zvis22 * ( ze11 + ze22 )
(  316)                  zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
(  317)                  zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
(  318)                  zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag
(  319)                  zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12
(  320)  
(  321)                  ze11   =  akappa(ji,jj-1,1,2)
(  322)                  ze12   = -akappa(ji,jj-1,2,1)
(  323)                  ze22   = +akappa(ji,jj-1,2,2)
(  324)                  ze21   = -akappa(ji,jj-1,1,1)
(  325)                  zvis11 = 2.0 * zviseta (ji,jj-1) + dm
(  326)                  zvis22 =       zviszeta(ji,jj-1) - zviseta(ji,jj-1)
(  327)                  zvis12 =       zviseta (ji,jj-1) + dm
(  328)                  zvis21 =       zviseta (ji,jj-1)
(  329) 
(  330)                  zdiag = zvis22 * ( ze11 + ze22 )
(  331)                  zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
(  332)                  zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
(  333)                  zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
(  334)                  zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12
(  335) 
(  336)                  ze11   =  akappa(ji-1,jj,1,2)
(  337)                  ze12   = -akappa(ji-1,jj,2,1)
(  338)                  ze22   = -akappa(ji-1,jj,2,2)
(  339)                  ze21   =  akappa(ji-1,jj,1,1)
(  340)                  zvis11 = 2.0 * zviseta (ji-1,jj) + dm
(  341)                  zvis22 =       zviszeta(ji-1,jj) - zviseta(ji-1,jj)
(  342)                  zvis12 =       zviseta (ji-1,jj) + dm
(  343)                  zvis21 =       zviseta (ji-1,jj)
(  344) 
(  345)                  zdiag = zvis22 * ( ze11 + ze22 )
(  346)                  zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
(  347)                  zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
(  348)                  zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
(  349)                  zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12
(  350) 
(  351)                  ze11   =  akappa(ji,jj,1,2)
(  352)                  ze12   = -akappa(ji,jj,2,1)
(  353)                  ze22   = -akappa(ji,jj,2,2)
(  354)                  ze21   = -akappa(ji,jj,1,1)
(  355)                  zvis11 = 2.0 * zviseta (ji,jj) + dm






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 19

(  356)                  zvis22 =       zviszeta(ji,jj) - zviseta(ji,jj)
(  357)                  zvis12 =       zviseta (ji,jj) + dm
(  358)                  zvis21 =       zviseta (ji,jj)
(  359) 
(  360)                  zdiag = zvis22 * ( ze11 + ze22 )
(  361)                  zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
(  362)                  zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
(  363)                  zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag
(  364)                  zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12
(  365) 
(  366)               END DO
(  367)            END DO
(  368) 
(  369)            DO jj = k_j1+1, k_jpj-1
(  370)               DO ji = 2, jpim1
(  371)                  zc1v(ji,jj) =   &
(  372)                     + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)   &
(  373)                     - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)   &
(  374)                     - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)   &
(  375)                     + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)   &
(  376)                     + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)   &
(  377)                     + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)   &
(  378)                     - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)   &
(  379)                     - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
(  380)                  zc2v(ji,jj) =   &
(  381)                     + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)   &
(  382)                     - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)   &
(  383)                     - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)   &
(  384)                     + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)   &
(  385)                     - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)   &
(  386)                     - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)   &
(  387)                     + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)   &
(  388)                     + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
(  389)               END DO
(  390)            END DO
(  391) 
(  392)          ! Relaxation.
(  393)            
(  394) iflag:   DO jter = 1 , nbitdr
(  395) 
(  396)             !  Store previous drift field.   
(  397)             DO jj = k_j1, k_jpj-1
(  398)                zu_ice(:,jj) = u_ice(:,jj)
(  399)                zv_ice(:,jj) = v_ice(:,jj)
(  400)             END DO
(  401) 
(  402)             DO jj = k_j1+1, k_jpj-1
(  403)                zsang  = SIGN( 1.e0, fcor(1,jj) ) * sangvg   ! only the sinus changes its sign with the hemisphere
(  404)                DO ji = 2, jpim1
(  405)                  zur     = u_ice(ji,jj) - u_oce(ji,jj)
(  406)                  zvr     = v_ice(ji,jj) - v_oce(ji,jj)
(  407)                  zmod    = SQRT( zur * zur + zvr * zvr) * ( 1.0 - zfrld(ji,jj) )
(  408)                  za      = rhoco * zmod
(  409)                  zac     = za * cangvg
(  410)                   zmpzas  = alpha * zcorl(ji,jj) + za * zsang
(  411)                   zmassdt = zusdtp * zmass(ji,jj)
(  412)                   zcorlal = ( 1.0 - alpha ) * zcorl(ji,jj)
(  413) 






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 20

(  414)                   za1(ji,jj) =  zmassdt * zu0(ji,jj) + zcorlal * zv0(ji,jj) + za1ct(ji,jj)   &
(  415)                      &        + za * ( cangvg * u_oce(ji,jj) - zsang * v_oce(ji,jj) )
(  416) 
(  417)                   za2(ji,jj) =  zmassdt * zv0(ji,jj) - zcorlal * zu0(ji,jj) + za2ct(ji,jj)   &
(  418)                      &        + za * ( cangvg * v_oce(ji,jj) + zsang * u_oce(ji,jj) )
(  419) 
(  420)                   zb1(ji,jj)  = zmassdt + zac - zc1u(ji,jj)
(  421)                   zb2(ji,jj)  = zmpzas  - zc2u(ji,jj)
(  422)                   zc1(ji,jj)  = zmpzas  + zc1v(ji,jj)
(  423)                   zc2(ji,jj)  = zmassdt + zac  - zc2v(ji,jj) 
(  424)                   zdeter      = zc1(ji,jj) * zb2(ji,jj) + zc2(ji,jj) * zb1(ji,jj)
(  425)                   zden(ji,jj) = SIGN( rone, zdeter) / MAX( epsd , ABS( zdeter ) )
(  426)                END DO
(  427)             END DO
(  428) 
(  429)             ! The computation of ice interaction term is splitted into two parts
(  430)             !-------------------------------------------------------------------------
(  431) 
(  432)             ! Terms that do not involve already up-dated velocities.
(  433)          
(  434)             DO jj = k_j1+1, k_jpj-1
(  435)                DO ji = 2, jpim1
(  436)                   iim1 = ji
(  437)                   ijm1 = jj - 1
(  438)                   iip1 = ji + 1
(  439)                   ijp1 = jj
(  440)                   ze11 =   akappa(iim1,ijm1,1,1) * u_ice(iip1,ijp1) + akappa(iim1,ijm1,1,2) * v_ice(iip1,ijp1)
(  441)                   ze12 = + akappa(iim1,ijm1,2,2) * u_ice(iip1,ijp1) - akappa(iim1,ijm1,2,1) * v_ice(iip1,ijp1)
(  442)                   ze22 = + akappa(iim1,ijm1,2,2) * v_ice(iip1,ijp1) + akappa(iim1,ijm1,2,1) * u_ice(iip1,ijp1)
(  443)                   ze21 =   akappa(iim1,ijm1,1,1) * v_ice(iip1,ijp1) - akappa(iim1,ijm1,1,2) * u_ice(iip1,ijp1)
(  444)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  445)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  446)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  447)                   zvis21 =       zviseta (iim1,ijm1)
(  448)                   zdiag = zvis22 * ( ze11 + ze22 )
(  449)                   zs11(ji,jj,2,1) =  zvis11 * ze11 + zdiag
(  450)                   zs12(ji,jj,2,1) =  zvis12 * ze12 + zvis21 * ze21
(  451)                   zs22(ji,jj,2,1) =  zvis11 * ze22 + zdiag
(  452)                   zs21(ji,jj,2,1) =  zvis12 * ze21 + zvis21 * ze12
(  453) 
(  454) 
(  455)                   iim1 = ji - 1
(  456)                   ijm1 = jj
(  457)                   iip1 = ji
(  458)                   ijp1 = jj + 1                   
(  459)                   ze11 =   akappa(iim1,ijm1,1,1) * ( u_ice(iip1,ijp1) - u_ice(iim1,ijp1) )   &
(  460)                      &   + akappa(iim1,ijm1,1,2) * ( v_ice(iip1,ijp1) + v_ice(iim1,ijp1) )
(  461)                   ze12 = + akappa(iim1,ijm1,2,2) * ( u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )   &
(  462)                      &   - akappa(iim1,ijm1,2,1) * ( v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )
(  463)                   ze22 = + akappa(iim1,ijm1,2,2) * ( v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )   &
(  464)                      &   + akappa(iim1,ijm1,2,1) * ( u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )
(  465)                   ze21 =   akappa(iim1,ijm1,1,1) * ( v_ice(iip1,ijp1) - v_ice(iim1,ijp1) )   &
(  466)                      &   - akappa(iim1,ijm1,1,2) * ( u_ice(iip1,ijp1) + u_ice(iim1,ijp1) )
(  467)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  468)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  469)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  470)                   zvis21 =       zviseta (iim1,ijm1)
(  471)                   zdiag = zvis22 * ( ze11 + ze22 )






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 21

(  472)                   zs11(ji,jj,1,2) =  zvis11 * ze11 + zdiag
(  473)                   zs12(ji,jj,1,2) =  zvis12 * ze12 + zvis21 * ze21
(  474)                   zs22(ji,jj,1,2) =  zvis11 * ze22 + zdiag
(  475)                   zs21(ji,jj,1,2) =  zvis12 * ze21 + zvis21 * ze12
(  476) 
(  477)                   iim1 = ji
(  478)                   ijm1 = jj
(  479)                   iip1 = ji + 1
(  480)                   ijp1 = jj + 1
(  481)                   ze11 =   akappa(iim1,ijm1,1,1) * ( u_ice(iip1,ijm1) + u_ice(iip1,ijp1) - u_ice(iim1,ijp1) )   &
(  482)                      &   + akappa(iim1,ijm1,1,2) * ( v_ice(iip1,ijm1) + v_ice(iip1,ijp1) + v_ice(iim1,ijp1) )
(  483)                   ze12 = - akappa(iim1,ijm1,2,2) * ( u_ice(iip1,ijm1) - u_ice(iim1,ijp1) - u_ice(iip1,ijp1) )   &
(  484)                      &   - akappa(iim1,ijm1,2,1) * ( v_ice(iip1,ijm1) + v_ice(iim1,ijp1) + v_ice(iip1,ijp1) )
(  485)                   ze22 = - akappa(iim1,ijm1,2,2) * ( v_ice(iip1,ijm1) - v_ice(iim1,ijp1) - v_ice(iip1,ijp1) )   &
(  486)                      &   + akappa(iim1,ijm1,2,1) * ( u_ice(iip1,ijm1) + u_ice(iim1,ijp1) + u_ice(iip1,ijp1) )
(  487)                   ze21 =   akappa(iim1,ijm1,1,1) * ( v_ice(iip1,ijm1) + v_ice(iip1,ijp1) - v_ice(iim1,ijp1) )   &
(  488)                      &   - akappa(iim1,ijm1,1,2) * ( u_ice(iip1,ijm1) + u_ice(iip1,ijp1) + u_ice(iim1,ijp1) ) 
(  489)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  490)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  491)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  492)                   zvis21 =       zviseta (iim1,ijm1)
(  493) 
(  494)                   zdiag = zvis22 * ( ze11 + ze22 )
(  495)                   zs11(ji,jj,2,2) =  zvis11 * ze11 + zdiag
(  496)                   zs12(ji,jj,2,2) =  zvis12 * ze12 + zvis21 * ze21
(  497)                   zs22(ji,jj,2,2) =  zvis11 * ze22 + zdiag
(  498)                   zs21(ji,jj,2,2) =  zvis12 * ze21 + zvis21 * ze12
(  499) 
(  500)                END DO
(  501)             END DO
(  502) 
(  503)             ! Terms involving already up-dated velocities.
(  504)             !-Using the arrays zu_ice and zv_ice in the computation of the terms ze leads to JACOBI's method; 
(  505)             ! Using arrays u and v in the computation of the terms ze leads to GAUSS-SEIDEL method.
(  506)              
(  507)             DO jj = k_j1+1, k_jpj-1
(  508)                DO ji = 2, jpim1
(  509)                   iim1 = ji - 1
(  510)                   ijm1 = jj - 1
(  511)                   iip1 = ji
(  512)                   ijp1 = jj
(  513)                   ze11 =   akappa(iim1,ijm1,1,1) * ( zu_ice(iip1,ijm1) - zu_ice(iim1,ijm1) - zu_ice(iim1,ijp1) )   &
(  514)                      &   + akappa(iim1,ijm1,1,2) * ( zv_ice(iip1,ijm1) + zv_ice(iim1,ijm1) + zv_ice(iim1,ijp1) )
(  515)                   ze12 = - akappa(iim1,ijm1,2,2) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) - zu_ice(iim1,ijp1) )   &
(  516)                      &   - akappa(iim1,ijm1,2,1) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) + zv_ice(iim1,ijp1) )
(  517)                   ze22 = - akappa(iim1,ijm1,2,2) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) - zv_ice(iim1,ijp1) )   &
(  518)                      &   + akappa(iim1,ijm1,2,1) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) + zu_ice(iim1,ijp1) )
(  519)                   ze21 =   akappa(iim1,ijm1,1,1) * ( zv_ice(iip1,ijm1) - zv_ice(iim1,ijm1) - zv_ice(iim1,ijp1) )   &
(  520)                      &  -  akappa(iim1,ijm1,1,2) * ( zu_ice(iip1,ijm1) + zu_ice(iim1,ijm1) + zu_ice(iim1,ijp1) )
(  521)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  522)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  523)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  524)                   zvis21 =       zviseta (iim1,ijm1)
(  525) 
(  526)                   zdiag = zvis22 * ( ze11 + ze22 )
(  527)                   zs11(ji,jj,1,1) =  zvis11 * ze11 + zdiag
(  528)                   zs12(ji,jj,1,1) =  zvis12 * ze12 + zvis21 * ze21
(  529)                   zs22(ji,jj,1,1) =  zvis11 * ze22 + zdiag






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 22

(  530)                   zs21(ji,jj,1,1) =  zvis12 * ze21 + zvis21 * ze12
(  531) 
(  532) 
(  533) # 540
(  540)                   iim1 = ji
(  541)                   ijm1 = jj - 1
(  542)                   iip1 = ji + 1
(  543)                   ze11 =   akappa(iim1,ijm1,1,1) * ( zu_ice(iip1,ijm1) - zu_ice(iim1,ijm1) )   &
(  544)                      &   + akappa(iim1,ijm1,1,2) * ( zv_ice(iip1,ijm1) + zv_ice(iim1,ijm1) )
(  545)                   ze12 = - akappa(iim1,ijm1,2,2) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) )   &
(  546)                      &   - akappa(iim1,ijm1,2,1) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) )
(  547)                   ze22 = - akappa(iim1,ijm1,2,2) * ( zv_ice(iim1,ijm1) + zv_ice(iip1,ijm1) )   &
(  548)                      &   + akappa(iim1,ijm1,2,1) * ( zu_ice(iim1,ijm1) + zu_ice(iip1,ijm1) )
(  549)                   ze21 =   akappa(iim1,ijm1,1,1) * ( zv_ice(iip1,ijm1) - zv_ice(iim1,ijm1) )   &
(  550)                      &   - akappa(iim1,ijm1,1,2) * ( zu_ice(iip1,ijm1) + zu_ice(iim1,ijm1) )
(  551)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  552)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  553)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  554)                   zvis21 =       zviseta (iim1,ijm1)
(  555) 
(  556)                   zdiag = zvis22 * ( ze11 + ze22 )
(  557)                   zs11(ji,jj,2,1) =  zs11(ji,jj,2,1) + zvis11 * ze11 + zdiag
(  558)                   zs12(ji,jj,2,1) =  zs12(ji,jj,2,1) + zvis12 * ze12 + zvis21 * ze21
(  559)                   zs22(ji,jj,2,1) =  zs22(ji,jj,2,1) + zvis11 * ze22 + zdiag
(  560)                   zs21(ji,jj,2,1) =  zs21(ji,jj,2,1) + zvis12 * ze21 + zvis21 * ze12
(  561) 
(  562) 
(  563)                   iim1 = ji - 1
(  564)                   ijm1 = jj 
(  565)                   ze11 = - akappa(iim1,ijm1,1,1) * zu_ice(iim1,ijm1) + akappa(iim1,ijm1,1,2) * zv_ice(iim1,ijm1)
(  566)                   ze12 = - akappa(iim1,ijm1,2,2) * zu_ice(iim1,ijm1) - akappa(iim1,ijm1,2,1) * zv_ice(iim1,ijm1)
(  567)                   ze22 = - akappa(iim1,ijm1,2,2) * zv_ice(iim1,ijm1) + akappa(iim1,ijm1,2,1) * zu_ice(iim1,ijm1)
(  568)                   ze21 = - akappa(iim1,ijm1,1,1) * zv_ice(iim1,ijm1) - akappa(iim1,ijm1,1,2) * zu_ice(iim1,ijm1)
(  569)                   zvis11 = 2.0 * zviseta (iim1,ijm1) + dm
(  570)                   zvis22 =       zviszeta(iim1,ijm1) - zviseta(iim1,ijm1)
(  571)                   zvis12 =       zviseta (iim1,ijm1) + dm
(  572)                   zvis21 =       zviseta (iim1,ijm1)
(  573) 
(  574)                   zdiag = zvis22 * ( ze11 + ze22 )
(  575)                   zs11(ji,jj,1,2) =  zs11(ji,jj,1,2) + zvis11 * ze11 + zdiag 
(  576)                   zs12(ji,jj,1,2) =  zs12(ji,jj,1,2) + zvis12 * ze12 + zvis21 * ze21
(  577)                   zs22(ji,jj,1,2) =  zs22(ji,jj,1,2) + zvis11 * ze22 + zdiag
(  578)                   zs21(ji,jj,1,2) =  zs21(ji,jj,1,2) + zvis12 * ze21 + zvis21 * ze12
(  579) 
(  580) # 587
(  587)                   zd1(ji,jj) =   &
(  588)                      + alambd(ji,jj,2,2,2,1) * zs11(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs11(ji,jj,2,2)  &
(  589)                      - alambd(ji,jj,2,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs11(ji,jj,1,2)  &
(  590)                      - alambd(ji,jj,1,1,2,1) * zs12(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs12(ji,jj,1,1)  &
(  591)                      + alambd(ji,jj,1,1,2,2) * zs12(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs12(ji,jj,1,2)  &
(  592)                      + alambd(ji,jj,1,2,1,1) * zs21(ji,jj,1,1) + alambd(ji,jj,1,2,2,1) * zs21(ji,jj,2,1)  &
(  593)                      + alambd(ji,jj,1,2,1,2) * zs21(ji,jj,1,2) + alambd(ji,jj,1,2,2,2) * zs21(ji,jj,2,2)  &
(  594)                      - alambd(ji,jj,2,1,1,1) * zs22(ji,jj,1,1) - alambd(ji,jj,2,1,2,1) * zs22(ji,jj,2,1)  &
(  595)                      - alambd(ji,jj,2,1,1,2) * zs22(ji,jj,1,2) - alambd(ji,jj,2,1,2,2) * zs22(ji,jj,2,2)
(  596)                   zd2(ji,jj) =   &
(  597)                      + alambd(ji,jj,2,2,2,1) * zs21(ji,jj,2,1) + alambd(ji,jj,2,2,2,2) * zs21(ji,jj,2,2)  &
(  598)                      - alambd(ji,jj,2,2,1,1) * zs21(ji,jj,1,1) - alambd(ji,jj,2,2,1,2) * zs21(ji,jj,1,2)  &
(  599)                      - alambd(ji,jj,1,1,2,1) * zs22(ji,jj,2,1) - alambd(ji,jj,1,1,1,1) * zs22(ji,jj,1,1)  &






PGF90 (Version     10.5)          09/23/2010  09:49:54      page 23

(  600)                      + alambd(ji,jj,1,1,2,2) * zs22(ji,jj,2,2) + alambd(ji,jj,1,1,1,2) * zs22(ji,jj,1,2)  &
(  601)                      - alambd(ji,jj,1,2,1,1) * zs11(ji,jj,1,1) - alambd(ji,jj,1,2,2,1) * zs11(ji,jj,2,1)  &
(  602)                      - alambd(ji,jj,1,2,1,2) * zs11(ji,jj,1,2) - alambd(ji,jj,1,2,2,2) * zs11(ji,jj,2,2)  &
(  603)                      + alambd(ji,jj,2,1,1,1) * zs12(ji,jj,1,1) + alambd(ji,jj,2,1,2,1) * zs12(ji,jj,2,1)  &
(  604)                      + alambd(ji,jj,2,1,1,2) * zs12(ji,jj,1,2) + alambd(ji,jj,2,1,2,2) * zs12(ji,jj,2,2)
(  605)                END DO
(  606)             END DO
(  607) 
(  608)             DO jj = k_j1+1, k_jpj-1
(  609)                DO ji = 2, jpim1
(  610)                   zunw = (  ( za1(ji,jj) + zd1(ji,jj) ) * zc2(ji,jj)        &
(  611)                      &    + ( za2(ji,jj) + zd2(ji,jj) ) * zc1(ji,jj) ) * zden(ji,jj)
(  612) 
(  613)                   zvnw = (  ( za2(ji,jj) + zd2(ji,jj) ) * zb1(ji,jj)        &
(  614)                      &    - ( za1(ji,jj) + zd1(ji,jj) ) * zb2(ji,jj) ) * zden(ji,jj)
(  615) 
(  616)                   zmask = ( 1.0 - MAX( rzero, SIGN( rone , 1.0 - zmass(ji,jj) ) ) ) * tmu(ji,jj)
(  617) 
(  618)                   u_ice(ji,jj) = ( u_ice(ji,jj) + om * ( zunw - u_ice(ji,jj) ) * tmu(ji,jj) ) * zmask
(  619)                   v_ice(ji,jj) = ( v_ice(ji,jj) + om * ( zvnw - v_ice(ji,jj) ) * tmu(ji,jj) ) * zmask
(  620)                END DO
(  621)             END DO
(  622) 
(  623)             CALL lbc_lnk( u_ice, 'I', -1. )
(  624)             CALL lbc_lnk( v_ice, 'I', -1. )
(  625) 
(  626)             !---  5.2.5.4. Convergence test.
(  627)             DO jj = k_j1+1 , k_jpj-1
(  628)                zresr(:,jj) = MAX( ABS( u_ice(:,jj) - zu_ice(:,jj) ) , ABS( v_ice(:,jj) - zv_ice(:,jj) ) )
(  629)             END DO
(  630)             zresm = MAXVAL( zresr( 1:jpi , k_j1+1:k_jpj-1 ) )
(  631)             IF( lk_mpp )   CALL mpp_max( zresm )   ! max over the global domain
(  632) 
(  633)             IF ( zresm <= resl) EXIT iflag
(  634) 
(  635)          END DO iflag
(  636) 
(  637)          zindu1 = 1.0 - zindu
(  638)          DO jj = k_j1 , k_jpj-1
(  639)             zu0(:,jj) = zindu * zu0(:,jj) + zindu1 * u_ice(:,jj)
(  640)             zv0(:,jj) = zindu * zv0(:,jj) + zindu1 * v_ice(:,jj)
(  641)          END DO
(  642)       !                                                   ! ==================== !
(  643)       END DO                                              !  end loop over iter  !
(  644)       !                                                   ! ==================== !
(  645) 
(  646)       IF(ln_ctl) THEN
(  647)          WRITE(charout,FMT="('lim_rhg  : res =',D23.16, ' iter =',I4)") zresm, jter
(  648)          CALL prt_ctl_info(charout)
(  649)          CALL prt_ctl(tab2d_1=u_ice, clinfo1=' lim_rhg  : u_ice :', tab2d_2=v_ice, clinfo2=' v_ice :')
(  650)       ENDIF
(  651) 
(  652)    END SUBROUTINE lim_rhg
(  653) 
(  654) 
(  655) # 665
(  665)    !!==============================================================================
(  666) END MODULE limrhg



