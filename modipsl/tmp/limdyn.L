


PGF90 (Version     10.5)          09/23/2010  09:49:58      page 1

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

Filename: limdyn.F90

(    1) MODULE limdyn
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  limdyn  ***
(    4)    !!   Sea-Ice dynamics :  
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!    lim_dyn      : computes ice velocities
(   11)    !!    lim_dyn_init : initialization and namelist read
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE dom_ice
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE ice
(   19)    USE ice_oce
(   20)    USE iceini
(   21)    USE limistate
(   22)    USE limrhg          ! ice rheology
(   23)    USE lbclnk
(   24)    USE lib_mpp
(   25)    USE prtctl          ! Print control
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Accessibility
(   31)    PUBLIC lim_dyn  ! routine called by ice_step
(   32) 
(   33)    !! * Module variables
(   34)    REAL(wp)  ::  rone    = 1.e0   ! constant value
(   35) 
(   36)    !!----------------------------------------------------------------------
(   37)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdyn.F90,v 1.7 2005/09/22 10:43:39 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE lim_dyn
(   45)       !!-------------------------------------------------------------------
(   46)       !!               ***  ROUTINE lim_dyn  ***






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 2

(   47)       !!               
(   48)       !! ** Purpose :   compute ice velocity and ocean-ice stress
(   49)       !!                
(   50)       !! ** Method  : 
(   51)       !!
(   52)       !! ** Action  : - Initialisation
(   53)       !!              - Call of the dynamic routine for each hemisphere
(   54)       !!              - computation of the stress at the ocean surface         
(   55)       !!              - treatment of the case if no ice dynamic
(   56)       !! History :
(   57)       !!   1.0  !  01-04  (LIM)  Original code
(   58)       !!   2.0  !  02-08  (C. Ethe, G. Madec)  F90, mpp
(   59)       !!---------------------------------------------------------------------
(   60)       !! * Loal variables
(   61)       INTEGER ::   ji, jj             ! dummy loop indices
(   62)       INTEGER ::   i_j1, i_jpj        ! Starting/ending j-indices for rheology
(   63)       REAL(wp) ::   &
(   64)          ztairx, ztairy,           &  ! tempory scalars
(   65)          zsang , zmod,             &
(   66)          ztglx , ztgly ,           &
(   67)          zt11, zt12, zt21, zt22 ,  &
(   68)          zustm, zsfrld, zsfrldm4,  &
(   69)          zu_ice, zv_ice, ztair2
(   70)       REAL(wp),DIMENSION(jpj) ::   &
(   71)          zind,                     &  ! i-averaged indicator of sea-ice
(   72)          zmsk                         ! i-averaged of tmask
(   73)       !!---------------------------------------------------------------------
(   74) 
(   75)       IF( numit == nstart  )   CALL lim_dyn_init   ! Initialization (first time-step only)
(   76)       
(   77)       IF ( ln_limdyn ) THEN
(   78) 
(   79)          ! Mean ice and snow thicknesses.          
(   80)          hsnm(:,:)  = ( 1.0 - frld(:,:) ) * hsnif(:,:)
(   81)          hicm(:,:)  = ( 1.0 - frld(:,:) ) * hicif(:,:)
(   82) 
(   83)          u_oce(:,:)  = u_io(:,:) * tmu(:,:)
(   84)          v_oce(:,:)  = v_io(:,:) * tmu(:,:)
(   85)        
(   86)          !                                         ! Rheology (ice dynamics)
(   87)          !                                         ! ========
(   88)          
(   89)          !  Define the j-limits where ice rheology is computed
(   90)          ! ---------------------------------------------------
(   91)          
(   92)          IF( lk_mpp ) THEN                    ! mpp: compute over the whole domain
(   93)             i_j1 = 1   
(   94)             i_jpj = jpj
(   95)             IF(ln_ctl)    THEN
(   96)                CALL prt_ctl_info('lim_dyn  :    i_j1 = ', ivar1=i_j1, clinfo2=' ij_jpj = ', ivar2=i_jpj)
(   97)             ENDIF
(   98)             CALL lim_rhg( i_j1, i_jpj )
(   99) 
(  100)          ELSE                                 ! optimization of the computational area
(  101) 
(  102)             DO jj = 1, jpj
(  103)                zind(jj) = SUM( frld (:,jj  ) )   ! = FLOAT(jpj) if ocean everywhere on a j-line
(  104)                zmsk(jj) = SUM( tmask(:,jj,1) )   ! = 0          if land  everywhere on a j-line






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 3

(  105)             END DO
(  106) 
(  107)             IF( l_jeq ) THEN                     ! local domain include both hemisphere
(  108)                !                                 ! Rheology is computed in each hemisphere
(  109)                !                                 ! only over the ice cover latitude strip
(  110)                ! Northern hemisphere
(  111)                i_j1  = njeq
(  112)                i_jpj = jpj
(  113)                DO WHILE ( i_j1 <= jpj .AND. zind(i_j1) == FLOAT(jpi) .AND. zmsk(i_j1) /=0 )
(  114)                   i_j1 = i_j1 + 1
(  115)                END DO
(  116)                i_j1 = MAX( 1, i_j1-1 )
(  117)                IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : NH i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
(  118)     
(  119)                CALL lim_rhg( i_j1, i_jpj )
(  120)     
(  121)                ! Southern hemisphere
(  122)                i_j1  =  1 
(  123)                i_jpj = njeq
(  124)                DO WHILE ( i_jpj >= 1 .AND. zind(i_jpj) == FLOAT(jpi) .AND. zmsk(i_jpj) /=0 )
(  125)                   i_jpj = i_jpj - 1
(  126)                END DO
(  127)                i_jpj = MIN( jpj, i_jpj+2 )
(  128)                IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : SH i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
(  129)     
(  130)                CALL lim_rhg( i_j1, i_jpj )
(  131)     
(  132)             ELSE                                 ! local domain extends over one hemisphere only
(  133)                !                                 ! Rheology is computed only over the ice cover
(  134)                !                                 ! latitude strip
(  135)                i_j1  = 1
(  136)                DO WHILE ( i_j1 <= jpj .AND. zind(i_j1) == FLOAT(jpi) .AND. zmsk(i_j1) /=0 )
(  137)                   i_j1 = i_j1 + 1
(  138)                END DO
(  139)                i_j1 = MAX( 1, i_j1-1 )
(  140)     
(  141)                i_jpj  = jpj
(  142)                DO WHILE ( i_jpj >= 1  .AND. zind(i_jpj) == FLOAT(jpi) .AND. zmsk(i_jpj) /=0 )
(  143)                   i_jpj = i_jpj - 1
(  144)                END DO
(  145)                i_jpj = MIN( jpj, i_jpj+2)
(  146)     
(  147)                IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : one hemisphere: i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
(  148)     
(  149)                CALL lim_rhg( i_j1, i_jpj )
(  150) 
(  151)             ENDIF
(  152) 
(  153)          ENDIF
(  154) 
(  155)          IF(ln_ctl)   THEN 
(  156)             CALL prt_ctl(tab2d_1=u_oce , clinfo1=' lim_dyn  : u_oce :', tab2d_2=v_oce , clinfo2=' v_oce :')
(  157)             CALL prt_ctl(tab2d_1=u_ice , clinfo1=' lim_dyn  : u_ice :', tab2d_2=v_ice , clinfo2=' v_ice :')
(  158)          ENDIF
(  159)          
(  160)          !                                         ! Ice-Ocean stress
(  161)          !                                         ! ================
(  162)          DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 4

(  163)             zsang  = SIGN(1.e0, gphif(1,jj-1) ) * sangvg
(  164)             DO ji = 2, jpim1
(  165)                ! computation of wind stress over ocean in X and Y direction
(  166) #if defined key_coupled && defined key_lim_cp1
(  167)                ztairx =  frld(ji-1,jj  ) * gtaux(ji-1,jj  ) + frld(ji,jj  ) * gtaux(ji,jj  )      &
(  168)                   &    + frld(ji-1,jj-1) * gtaux(ji-1,jj-1) + frld(ji,jj-1) * gtaux(ji,jj-1)
(  169) 
(  170)                ztairy =  frld(ji-1,jj  ) * gtauy(ji-1,jj  ) + frld(ji,jj  ) * gtauy(ji,jj  )      &
(  171)                   &    + frld(ji-1,jj-1) * gtauy(ji-1,jj-1) + frld(ji,jj-1) * gtauy(ji,jj-1)
(  172) #else
(  173)                zsfrld  = frld(ji,jj) + frld(ji-1,jj) + frld(ji-1,jj-1) + frld(ji,jj-1)
(  174)                ztairx  = zsfrld * gtaux(ji,jj)
(  175)                ztairy  = zsfrld * gtauy(ji,jj)
(  176) #endif
(  177)                zsfrldm4 = 4 - frld(ji,jj) - frld(ji-1,jj) - frld(ji-1,jj-1) - frld(ji,jj-1)
(  178)                zu_ice   = u_ice(ji,jj) - u_oce(ji,jj)
(  179)                zv_ice   = v_ice(ji,jj) - v_oce(ji,jj)
(  180)                zmod     = SQRT( zu_ice * zu_ice + zv_ice * zv_ice ) 
(  181)                ztglx   = zsfrldm4 * rhoco * zmod * ( cangvg * zu_ice - zsang * zv_ice ) 
(  182)                ztgly   = zsfrldm4 * rhoco * zmod * ( cangvg * zv_ice + zsang * zu_ice ) 
(  183) 
(  184)                tio_u(ji,jj) = - ( ztairx + 1.0 * ztglx ) / ( 4 * rau0 )
(  185)                tio_v(ji,jj) = - ( ztairy + 1.0 * ztgly ) / ( 4 * rau0 )
(  186)             END DO
(  187)          END DO
(  188)          
(  189)          ! computation of friction velocity
(  190)          DO jj = 2, jpjm1
(  191)             DO ji = 2, jpim1
(  192) 
(  193)                zu_ice   = u_ice(ji-1,jj-1) - u_oce(ji-1,jj-1)
(  194)                zv_ice   = v_ice(ji-1,jj-1) - v_oce(ji-1,jj-1)
(  195)                zt11  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice )
(  196) 
(  197)                zu_ice   = u_ice(ji-1,jj) - u_oce(ji-1,jj)
(  198)                zv_ice   = v_ice(ji-1,jj) - v_oce(ji-1,jj)
(  199)                zt12  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 
(  200) 
(  201)                zu_ice   = u_ice(ji,jj-1) - u_oce(ji,jj-1)
(  202)                zv_ice   = v_ice(ji,jj-1) - v_oce(ji,jj-1)
(  203)                zt21  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 
(  204) 
(  205)                zu_ice   = u_ice(ji,jj) - u_oce(ji,jj)
(  206)                zv_ice   = v_ice(ji,jj) - v_oce(ji,jj)
(  207)                zt22  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 
(  208) 
(  209)                ztair2 = gtaux(ji,jj) * gtaux(ji,jj) + gtauy(ji,jj) * gtauy(ji,jj)
(  210) 
(  211)                zustm =  ( 1 - frld(ji,jj) ) * 0.25 * ( zt11 + zt12 + zt21 + zt22 )        &
(  212)                   &  +        frld(ji,jj)   * SQRT( ztair2 )
(  213) 
(  214)                ust2s(ji,jj) = ( zustm / rau0 ) * ( rone + sdvt(ji,jj) ) * tms(ji,jj)
(  215)             END DO
(  216)          END DO
(  217) 
(  218)        ELSE      ! no ice dynamics : transmit directly the atmospheric stress to the ocean
(  219)                     
(  220)           DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 5

(  221)              DO ji = 2, jpim1
(  222) #if defined key_coupled && defined key_lim_cp1
(  223)                 tio_u(ji,jj) = - (  gtaux(ji  ,jj  ) + gtaux(ji-1,jj  )       &
(  224)                    &              + gtaux(ji-1,jj-1) + gtaux(ji  ,jj-1) ) / ( 4 * rau0 )
(  225) 
(  226)                 tio_v(ji,jj) = - (  gtauy(ji  ,jj )  + gtauy(ji-1,jj  )       &
(  227)                    &              + gtauy(ji-1,jj-1) + gtauy(ji  ,jj-1) ) / ( 4 * rau0 )
(  228) #else
(  229)                 tio_u(ji,jj) = - gtaux(ji,jj) / rau0
(  230)                 tio_v(ji,jj) = - gtauy(ji,jj) / rau0 
(  231) #endif
(  232)                 ztair2       = gtaux(ji,jj) * gtaux(ji,jj) + gtauy(ji,jj) * gtauy(ji,jj)
(  233)                 zustm        = SQRT( ztair2  )
(  234) 
(  235)                 ust2s(ji,jj) = ( zustm / rau0 ) * ( rone + sdvt(ji,jj) ) * tms(ji,jj)
(  236)             END DO
(  237)          END DO
(  238) 
(  239)       ENDIF
(  240) 
(  241)       CALL lbc_lnk( ust2s, 'T',  1. )   ! T-point
(  242)       CALL lbc_lnk( tio_u, 'I', -1. )   ! I-point (i.e. ice U-V point)
(  243)       CALL lbc_lnk( tio_v, 'I', -1. )   ! I-point (i.e. ice U-V point)
(  244) 
(  245)       IF(ln_ctl) THEN 
(  246)             CALL prt_ctl(tab2d_1=tio_u , clinfo1=' lim_dyn  : tio_u :', tab2d_2=tio_v , clinfo2=' tio_v :')
(  247)             CALL prt_ctl(tab2d_1=ust2s , clinfo1=' lim_dyn  : ust2s :')
(  248)       ENDIF
(  249) 
(  250)    END SUBROUTINE lim_dyn
(  251) 
(  252) 
(  253)    SUBROUTINE lim_dyn_init
(  254)       !!-------------------------------------------------------------------
(  255)       !!                  ***  ROUTINE lim_dyn_init  ***
(  256)       !!
(  257)       !! ** Purpose : Physical constants and parameters linked to the ice
(  258)       !!      dynamics
(  259)       !!
(  260)       !! ** Method  :  Read the namicedyn namelist and check the ice-dynamic
(  261)       !!       parameter values called at the first timestep (nit000)
(  262)       !!
(  263)       !! ** input   :   Namelist namicedyn
(  264)       !!
(  265)       !! history :
(  266)       !!  8.5  ! 03-08 (C. Ethe) original code
(  267)       !!-------------------------------------------------------------------
(  268)       NAMELIST/namicedyn/ epsd, alpha,     &
(  269)          &                dm, nbiter, nbitdr, om, resl, cw, angvg, pstar,   &
(  270)          &                c_rhg, etamn, creepl, ecc, ahi0
(  271)       !!-------------------------------------------------------------------
(  272) 
(  273)       ! Define the initial parameters
(  274)       ! -------------------------
(  275) 
(  276)       ! Read Namelist namicedyn
(  277)       REWIND ( numnam_ice )
(  278)       READ   ( numnam_ice  , namicedyn )






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 6

(  279)       IF(lwp) THEN
(  280)          WRITE(numout,*)
(  281)          WRITE(numout,*) 'lim_dyn_init : ice parameters for ice dynamics '
(  282)          WRITE(numout,*) '~~~~~~~~~~~~'
(  283)          WRITE(numout,*) '       tolerance parameter                              epsd   = ', epsd
(  284)          WRITE(numout,*) '       coefficient for semi-implicit coriolis           alpha  = ', alpha
(  285)          WRITE(numout,*) '       diffusion constant for dynamics                  dm     = ', dm
(  286)          WRITE(numout,*) '       number of sub-time steps for relaxation          nbiter = ', nbiter
(  287)          WRITE(numout,*) '       maximum number of iterations for relaxation      nbitdr = ', nbitdr
(  288)          WRITE(numout,*) '       relaxation constant                              om     = ', om
(  289)          WRITE(numout,*) '       maximum value for the residual of relaxation     resl   = ', resl
(  290)          WRITE(numout,*) '       drag coefficient for oceanic stress              cw     = ', cw
(  291)          WRITE(numout,*) '       turning angle for oceanic stress                 angvg  = ', angvg
(  292)          WRITE(numout,*) '       first bulk-rheology parameter                    pstar  = ', pstar
(  293)          WRITE(numout,*) '       second bulk-rhelogy parameter                    c_rhg  = ', c_rhg
(  294)          WRITE(numout,*) '       minimun value for viscosity                      etamn  = ', etamn
(  295)          WRITE(numout,*) '       creep limit                                      creepl = ', creepl
(  296)          WRITE(numout,*) '       eccentricity of the elliptical yield curve       ecc    = ', ecc
(  297)          WRITE(numout,*) '       horizontal diffusivity coeff. for sea-ice        ahi0   = ', ahi0
(  298)       ENDIF
(  299) 
(  300)       usecc2 = 1.0 / ( ecc * ecc )
(  301)       rhoco  = rau0 * cw
(  302)       angvg  = angvg * rad
(  303)       sangvg = SIN( angvg )
(  304)       cangvg = COS( angvg )
(  305)       pstarh = pstar / 2.0
(  306) 
(  307)       !  Diffusion coefficients.
(  308)       ahiu(:,:) = ahi0 * umask(:,:,1)
(  309)       ahiv(:,:) = ahi0 * vmask(:,:,1)
(  310) 
(  311)    END SUBROUTINE lim_dyn_init
(  312) 
(  313) #else
(  314)    !!----------------------------------------------------------------------
(  315)    !!   Default option          Empty module           NO LIM sea-ice model
(  316)    !!----------------------------------------------------------------------
(  317) CONTAINS
(  318)    SUBROUTINE lim_dyn         ! Empty routine
(  319)    END SUBROUTINE lim_dyn
(  320) #endif 
(  321) 
(  322)    !!======================================================================
(  323) END MODULE limdyn



















PGF90 (Version     10.5)          09/23/2010  09:49:58      page 7

(    1) # 1 "limdyn.F90"
(    1) MODULE limdyn
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  limdyn  ***
(    4)    !!   Sea-Ice dynamics :  
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!    lim_dyn      : computes ice velocities
(   11)    !!    lim_dyn_init : initialization and namelist read
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst
(   15)    USE in_out_manager  ! I/O manager
(   16)    USE dom_ice
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE ice
(   19)    USE ice_oce
(   20)    USE iceini
(   21)    USE limistate
(   22)    USE limrhg          ! ice rheology
(   23)    USE lbclnk
(   24)    USE lib_mpp
(   25)    USE prtctl          ! Print control
(   26) 
(   27)    IMPLICIT NONE
(   28)    PRIVATE
(   29) 
(   30)    !! * Accessibility
(   31)    PUBLIC lim_dyn  ! routine called by ice_step
(   32) 
(   33)    !! * Module variables
(   34)    REAL(wp)  ::  rone    = 1.e0   ! constant value
(   35) 
(   36)    !!----------------------------------------------------------------------
(   37)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limdyn.F90,v 1.7 2005/09/22 10:43:39 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43) 
(   44)    SUBROUTINE lim_dyn
(   45)       !!-------------------------------------------------------------------
(   46)       !!               ***  ROUTINE lim_dyn  ***
(   47)       !!               
(   48)       !! ** Purpose :   compute ice velocity and ocean-ice stress
(   49)       !!                
(   50)       !! ** Method  : 
(   51)       !!
(   52)       !! ** Action  : - Initialisation
(   53)       !!              - Call of the dynamic routine for each hemisphere
(   54)       !!              - computation of the stress at the ocean surface         
(   55)       !!              - treatment of the case if no ice dynamic
(   56)       !! History :
(   57)       !!   1.0  !  01-04  (LIM)  Original code






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 8

(   58)       !!   2.0  !  02-08  (C. Ethe, G. Madec)  F90, mpp
(   59)       !!---------------------------------------------------------------------
(   60)       !! * Loal variables
(   61)       INTEGER ::   ji, jj             ! dummy loop indices
(   62)       INTEGER ::   i_j1, i_jpj        ! Starting/ending j-indices for rheology
(   63)       REAL(wp) ::   &
(   64)          ztairx, ztairy,           &  ! tempory scalars
(   65)          zsang , zmod,             &
(   66)          ztglx , ztgly ,           &
(   67)          zt11, zt12, zt21, zt22 ,  &
(   68)          zustm, zsfrld, zsfrldm4,  &
(   69)          zu_ice, zv_ice, ztair2
(   70)       REAL(wp),DIMENSION(jpj) ::   &
(   71)          zind,                     &  ! i-averaged indicator of sea-ice
(   72)          zmsk                         ! i-averaged of tmask
(   73)       !!---------------------------------------------------------------------
(   74) 
(   75)       IF( numit == nstart  )   CALL lim_dyn_init   ! Initialization (first time-step only)
(   76)       
(   77)       IF ( ln_limdyn ) THEN
(   78) 
(   79)          ! Mean ice and snow thicknesses.          
(   80)          hsnm(:,:)  = ( 1.0 - frld(:,:) ) * hsnif(:,:)
(   81)          hicm(:,:)  = ( 1.0 - frld(:,:) ) * hicif(:,:)
(   82) 
(   83)          u_oce(:,:)  = u_io(:,:) * tmu(:,:)
(   84)          v_oce(:,:)  = v_io(:,:) * tmu(:,:)
(   85)        
(   86)          !                                         ! Rheology (ice dynamics)
(   87)          !                                         ! ========
(   88)          
(   89)          !  Define the j-limits where ice rheology is computed
(   90)          ! ---------------------------------------------------
(   91)          
(   92)          IF( lk_mpp ) THEN                    ! mpp: compute over the whole domain
(   93)             i_j1 = 1   
(   94)             i_jpj = jpj
(   95)             IF(ln_ctl)    THEN
(   96)                CALL prt_ctl_info('lim_dyn  :    i_j1 = ', ivar1=i_j1, clinfo2=' ij_jpj = ', ivar2=i_jpj)
(   97)             ENDIF
(   98)             CALL lim_rhg( i_j1, i_jpj )
(   99) 
(  100)          ELSE                                 ! optimization of the computational area
(  101) 
(  102)             DO jj = 1, jpj
(  103)                zind(jj) = SUM( frld (:,jj  ) )   ! = FLOAT(jpj) if ocean everywhere on a j-line
(  104)                zmsk(jj) = SUM( tmask(:,jj,1) )   ! = 0          if land  everywhere on a j-line
(  105)             END DO
(  106) 
(  107)             IF( l_jeq ) THEN                     ! local domain include both hemisphere
(  108)                !                                 ! Rheology is computed in each hemisphere
(  109)                !                                 ! only over the ice cover latitude strip
(  110)                ! Northern hemisphere
(  111)                i_j1  = njeq
(  112)                i_jpj = jpj
(  113)                DO WHILE ( i_j1 <= jpj .AND. zind(i_j1) == FLOAT(jpi) .AND. zmsk(i_j1) /=0 )
(  114)                   i_j1 = i_j1 + 1
(  115)                END DO






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 9

(  116)                i_j1 = MAX( 1, i_j1-1 )
(  117)                IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : NH i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
(  118)     
(  119)                CALL lim_rhg( i_j1, i_jpj )
(  120)     
(  121)                ! Southern hemisphere
(  122)                i_j1  =  1 
(  123)                i_jpj = njeq
(  124)                DO WHILE ( i_jpj >= 1 .AND. zind(i_jpj) == FLOAT(jpi) .AND. zmsk(i_jpj) /=0 )
(  125)                   i_jpj = i_jpj - 1
(  126)                END DO
(  127)                i_jpj = MIN( jpj, i_jpj+2 )
(  128)                IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : SH i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
(  129)     
(  130)                CALL lim_rhg( i_j1, i_jpj )
(  131)     
(  132)             ELSE                                 ! local domain extends over one hemisphere only
(  133)                !                                 ! Rheology is computed only over the ice cover
(  134)                !                                 ! latitude strip
(  135)                i_j1  = 1
(  136)                DO WHILE ( i_j1 <= jpj .AND. zind(i_j1) == FLOAT(jpi) .AND. zmsk(i_j1) /=0 )
(  137)                   i_j1 = i_j1 + 1
(  138)                END DO
(  139)                i_j1 = MAX( 1, i_j1-1 )
(  140)     
(  141)                i_jpj  = jpj
(  142)                DO WHILE ( i_jpj >= 1  .AND. zind(i_jpj) == FLOAT(jpi) .AND. zmsk(i_jpj) /=0 )
(  143)                   i_jpj = i_jpj - 1
(  144)                END DO
(  145)                i_jpj = MIN( jpj, i_jpj+2)
(  146)     
(  147)                IF(ln_ctl .AND. lwp)   WRITE(numout,*) 'lim_dyn : one hemisphere: i_j1 = ', i_j1, ' ij_jpj = ', i_jpj
(  148)     
(  149)                CALL lim_rhg( i_j1, i_jpj )
(  150) 
(  151)             ENDIF
(  152) 
(  153)          ENDIF
(  154) 
(  155)          IF(ln_ctl)   THEN 
(  156)             CALL prt_ctl(tab2d_1=u_oce , clinfo1=' lim_dyn  : u_oce :', tab2d_2=v_oce , clinfo2=' v_oce :')
(  157)             CALL prt_ctl(tab2d_1=u_ice , clinfo1=' lim_dyn  : u_ice :', tab2d_2=v_ice , clinfo2=' v_ice :')
(  158)          ENDIF
(  159)          
(  160)          !                                         ! Ice-Ocean stress
(  161)          !                                         ! ================
(  162)          DO jj = 2, jpjm1
(  163)             zsang  = SIGN(1.e0, gphif(1,jj-1) ) * sangvg
(  164)             DO ji = 2, jpim1
(  165)                ! computation of wind stress over ocean in X and Y direction
(  166) # 173
(  173)                zsfrld  = frld(ji,jj) + frld(ji-1,jj) + frld(ji-1,jj-1) + frld(ji,jj-1)
(  174)                ztairx  = zsfrld * gtaux(ji,jj)
(  175)                ztairy  = zsfrld * gtauy(ji,jj)
(  176) # 177
(  177)                zsfrldm4 = 4 - frld(ji,jj) - frld(ji-1,jj) - frld(ji-1,jj-1) - frld(ji,jj-1)
(  178)                zu_ice   = u_ice(ji,jj) - u_oce(ji,jj)
(  179)                zv_ice   = v_ice(ji,jj) - v_oce(ji,jj)






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 10

(  180)                zmod     = SQRT( zu_ice * zu_ice + zv_ice * zv_ice ) 
(  181)                ztglx   = zsfrldm4 * rhoco * zmod * ( cangvg * zu_ice - zsang * zv_ice ) 
(  182)                ztgly   = zsfrldm4 * rhoco * zmod * ( cangvg * zv_ice + zsang * zu_ice ) 
(  183) 
(  184)                tio_u(ji,jj) = - ( ztairx + 1.0 * ztglx ) / ( 4 * rau0 )
(  185)                tio_v(ji,jj) = - ( ztairy + 1.0 * ztgly ) / ( 4 * rau0 )
(  186)             END DO
(  187)          END DO
(  188)          
(  189)          ! computation of friction velocity
(  190)          DO jj = 2, jpjm1
(  191)             DO ji = 2, jpim1
(  192) 
(  193)                zu_ice   = u_ice(ji-1,jj-1) - u_oce(ji-1,jj-1)
(  194)                zv_ice   = v_ice(ji-1,jj-1) - v_oce(ji-1,jj-1)
(  195)                zt11  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice )
(  196) 
(  197)                zu_ice   = u_ice(ji-1,jj) - u_oce(ji-1,jj)
(  198)                zv_ice   = v_ice(ji-1,jj) - v_oce(ji-1,jj)
(  199)                zt12  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 
(  200) 
(  201)                zu_ice   = u_ice(ji,jj-1) - u_oce(ji,jj-1)
(  202)                zv_ice   = v_ice(ji,jj-1) - v_oce(ji,jj-1)
(  203)                zt21  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 
(  204) 
(  205)                zu_ice   = u_ice(ji,jj) - u_oce(ji,jj)
(  206)                zv_ice   = v_ice(ji,jj) - v_oce(ji,jj)
(  207)                zt22  = rhoco * ( zu_ice * zu_ice + zv_ice * zv_ice ) 
(  208) 
(  209)                ztair2 = gtaux(ji,jj) * gtaux(ji,jj) + gtauy(ji,jj) * gtauy(ji,jj)
(  210) 
(  211)                zustm =  ( 1 - frld(ji,jj) ) * 0.25 * ( zt11 + zt12 + zt21 + zt22 )        &
(  212)                   &  +        frld(ji,jj)   * SQRT( ztair2 )
(  213) 
(  214)                ust2s(ji,jj) = ( zustm / rau0 ) * ( rone + sdvt(ji,jj) ) * tms(ji,jj)
(  215)             END DO
(  216)          END DO
(  217) 
(  218)        ELSE      ! no ice dynamics : transmit directly the atmospheric stress to the ocean
(  219)                     
(  220)           DO jj = 2, jpjm1
(  221)              DO ji = 2, jpim1
(  222) # 229
(  229)                 tio_u(ji,jj) = - gtaux(ji,jj) / rau0
(  230)                 tio_v(ji,jj) = - gtauy(ji,jj) / rau0 
(  231) # 232
(  232)                 ztair2       = gtaux(ji,jj) * gtaux(ji,jj) + gtauy(ji,jj) * gtauy(ji,jj)
(  233)                 zustm        = SQRT( ztair2  )
(  234) 
(  235)                 ust2s(ji,jj) = ( zustm / rau0 ) * ( rone + sdvt(ji,jj) ) * tms(ji,jj)
(  236)             END DO
(  237)          END DO
(  238) 
(  239)       ENDIF
(  240) 
(  241)       CALL lbc_lnk( ust2s, 'T',  1. )   ! T-point
(  242)       CALL lbc_lnk( tio_u, 'I', -1. )   ! I-point (i.e. ice U-V point)
(  243)       CALL lbc_lnk( tio_v, 'I', -1. )   ! I-point (i.e. ice U-V point)






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 11

(  244) 
(  245)       IF(ln_ctl) THEN 
(  246)             CALL prt_ctl(tab2d_1=tio_u , clinfo1=' lim_dyn  : tio_u :', tab2d_2=tio_v , clinfo2=' tio_v :')
(  247)             CALL prt_ctl(tab2d_1=ust2s , clinfo1=' lim_dyn  : ust2s :')
(  248)       ENDIF
(  249) 
(  250)    END SUBROUTINE lim_dyn
(  251) 
(  252) 
(  253)    SUBROUTINE lim_dyn_init
(  254)       !!-------------------------------------------------------------------
(  255)       !!                  ***  ROUTINE lim_dyn_init  ***
(  256)       !!
(  257)       !! ** Purpose : Physical constants and parameters linked to the ice
(  258)       !!      dynamics
(  259)       !!
(  260)       !! ** Method  :  Read the namicedyn namelist and check the ice-dynamic
(  261)       !!       parameter values called at the first timestep (nit000)
(  262)       !!
(  263)       !! ** input   :   Namelist namicedyn
(  264)       !!
(  265)       !! history :
(  266)       !!  8.5  ! 03-08 (C. Ethe) original code
(  267)       !!-------------------------------------------------------------------
(  268)       NAMELIST/namicedyn/ epsd, alpha,     &
(  269)          &                dm, nbiter, nbitdr, om, resl, cw, angvg, pstar,   &
(  270)          &                c_rhg, etamn, creepl, ecc, ahi0
(  271)       !!-------------------------------------------------------------------
(  272) 
(  273)       ! Define the initial parameters
(  274)       ! -------------------------
(  275) 
(  276)       ! Read Namelist namicedyn
(  277)       REWIND ( numnam_ice )
(  278)       READ   ( numnam_ice  , namicedyn )
(  279)       IF(lwp) THEN
(  280)          WRITE(numout,*)
(  281)          WRITE(numout,*) 'lim_dyn_init : ice parameters for ice dynamics '
(  282)          WRITE(numout,*) '~~~~~~~~~~~~'
(  283)          WRITE(numout,*) '       tolerance parameter                              epsd   = ', epsd
(  284)          WRITE(numout,*) '       coefficient for semi-implicit coriolis           alpha  = ', alpha
(  285)          WRITE(numout,*) '       diffusion constant for dynamics                  dm     = ', dm
(  286)          WRITE(numout,*) '       number of sub-time steps for relaxation          nbiter = ', nbiter
(  287)          WRITE(numout,*) '       maximum number of iterations for relaxation      nbitdr = ', nbitdr
(  288)          WRITE(numout,*) '       relaxation constant                              om     = ', om
(  289)          WRITE(numout,*) '       maximum value for the residual of relaxation     resl   = ', resl
(  290)          WRITE(numout,*) '       drag coefficient for oceanic stress              cw     = ', cw
(  291)          WRITE(numout,*) '       turning angle for oceanic stress                 angvg  = ', angvg
(  292)          WRITE(numout,*) '       first bulk-rheology parameter                    pstar  = ', pstar
(  293)          WRITE(numout,*) '       second bulk-rhelogy parameter                    c_rhg  = ', c_rhg
(  294)          WRITE(numout,*) '       minimun value for viscosity                      etamn  = ', etamn
(  295)          WRITE(numout,*) '       creep limit                                      creepl = ', creepl
(  296)          WRITE(numout,*) '       eccentricity of the elliptical yield curve       ecc    = ', ecc
(  297)          WRITE(numout,*) '       horizontal diffusivity coeff. for sea-ice        ahi0   = ', ahi0
(  298)       ENDIF
(  299) 
(  300)       usecc2 = 1.0 / ( ecc * ecc )
(  301)       rhoco  = rau0 * cw






PGF90 (Version     10.5)          09/23/2010  09:49:58      page 12

(  302)       angvg  = angvg * rad
(  303)       sangvg = SIN( angvg )
(  304)       cangvg = COS( angvg )
(  305)       pstarh = pstar / 2.0
(  306) 
(  307)       !  Diffusion coefficients.
(  308)       ahiu(:,:) = ahi0 * umask(:,:,1)
(  309)       ahiv(:,:) = ahi0 * vmask(:,:,1)
(  310) 
(  311)    END SUBROUTINE lim_dyn_init
(  312) 
(  313) 
(  314) # 322
(  322)    !!======================================================================
(  323) END MODULE limdyn
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (limdyn.F90: 113)
