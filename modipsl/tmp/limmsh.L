


PGF90 (Version     10.5)          09/23/2010  09:49:50      page 1

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

Filename: limmsh.F90

(    1) MODULE limmsh
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  limmsh  ***
(    4)    !! LIM ice model :   definition of the ice mesh parameters
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_msh   : definition of the ice mesh
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE phycst
(   14)    USE dom_oce
(   15)    USE dom_ice
(   16)    USE lbclnk
(   17)    USE in_out_manager
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Accessibility
(   23)    PUBLIC lim_msh      ! routine called by ice_ini.F90
(   24) 
(   25)    !!----------------------------------------------------------------------
(   26)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limmsh.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32) 
(   33)    SUBROUTINE lim_msh
(   34)       !!-------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE lim_msh  ***
(   36)       !!              
(   37)       !! ** Purpose : Definition of the charact. of the numerical grid
(   38)       !!       
(   39)       !! ** Action  : - Initialisation of some variables
(   40)       !!              - Definition of some constants linked with the grid
(   41)       !!              - Definition of the metric coef. for the sea/ice
(   42)       !!              - Initialization of the ice masks (tmsk, umsk)
(   43)       !! 
(   44)       !! ** Refer.  : Deleersnijder et al. Ocean Modelling 100, 7-10 
(   45)       !!
(   46)       !! ** History :






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 2

(   47)       !!         original    : 01-04 (LIM)
(   48)       !!         addition    : 02-08 (C. Ethe, G. Madec)
(   49)       !!--------------------------------------------------------------------- 
(   50)       !! * Local variables
(   51)       INTEGER :: ji, jj      ! dummy loop indices
(   52) 
(   53)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   54)          zd2d1 , zd1d2       ! Derivative of zh2 (resp. zh1) in the x direction
(   55)          !                   ! (resp. y direction) (defined at the center)
(   56)       REAL(wp) ::         &
(   57)          zh1p  , zh2p   , &  ! Idem zh1, zh2 for the bottom left corner of the grid
(   58)          zd2d1p, zd1d2p , &  ! Idem zd2d1, zd1d2 for the bottom left corner of the grid
(   59)          zusden, zusden2     ! temporary scalars
(   60)       !!---------------------------------------------------------------------
(   61) 
(   62)       IF(lwp) THEN
(   63)          WRITE(numout,*)
(   64)          WRITE(numout,*) 'lim_msh : LIM sea-ice model, mesh initialization'
(   65)          WRITE(numout,*) '~~~~~~~'
(   66)       ENDIF
(   67)       
(   68)       !----------------------------------------------------------                          
(   69)       !    Initialization of local and some global (common) variables 
(   70)       !------------------------------------------------------------------ 
(   71)       
(   72)       njeq   = INT( jpj / 2 )   !i bug mpp potentiel
(   73)       njeqm1 = njeq - 1 
(   74) 
(   75)       fcor(:,:) = 2. * omega * SIN( gphit(:,:) * rad )   !  coriolis factor
(   76)  
(   77) !i    DO jj = 1, jpj
(   78) !i       zmsk(jj) = SUM( tmask(:,jj,:) )   ! = 0          if land  everywhere on a j-line
(   79) !!ii     write(numout,*) jj, zind(jj)
(   80) !i    END DO
(   81) 
(   82)       IF( fcor(1,1) * fcor(1,nlcj) < 0.e0 ) THEN   ! local domain include both hemisphere
(   83)          l_jeq = .TRUE.
(   84)          njeq  = 1
(   85)          DO WHILE ( njeq <= jpj .AND. fcor(1,njeq) < 0.e0 )
(   86)             njeq = njeq + 1
(   87)          END DO
(   88)          IF(lwp ) WRITE(numout,*) '          the equator is inside the domain at about njeq = ', njeq
(   89)       ELSEIF( fcor(1,1) < 0.e0 ) THEN
(   90)          l_jeq = .FALSE.
(   91)          njeq = jpj
(   92)          IF(lwp ) WRITE(numout,*) '          the model domain is entirely in the southern hemisphere: njeq = ', njeq
(   93)       ELSE
(   94)          l_jeq = .FALSE.
(   95)          njeq = 2
(   96)          IF(lwp ) WRITE(numout,*) '          the model domain is entirely in the northern hemisphere: njeq = ', njeq
(   97)       ENDIF
(   98) 
(   99)       njeqm1 = njeq - 1
(  100) 
(  101) 
(  102)       !   For each grid, definition of geometric tables 
(  103)       !------------------------------------------------------------------
(  104)       






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 3

(  105)       !-------------------
(  106)       ! Conventions :    !
(  107)       !-------------------
(  108)       !  indices 1 \ 2 <-> localisation in the 2 direction x \ y
(  109)       !  3rd indice <-> localisation on the mesh :
(  110)       !  0 = Centre ;  1 = corner W x(i-1/2) ; 2 = corner S y(j-1/2) ;
(  111)       !  3 = corner SW x(i-1/2),y(j-1/2)
(  112)       !-------------------
(  113) !!ibug ???
(  114)       akappa(:,:,:,:) = 0.e0
(  115)       wght(:,:,:,:) = 0.e0
(  116)       alambd(:,:,:,:,:,:) = 0.e0
(  117)       tmu(:,:) = 0.e0
(  118) !!i
(  119)       
(  120)       
(  121)       ! metric coefficients for sea ice dynamic
(  122)       !----------------------------------------
(  123)       !                                                       ! akappa
(  124)       DO jj = 2, jpj
(  125)          zd1d2(:,jj) = e1v(:,jj) - e1v(:,jj-1)
(  126)       END DO
(  127)       CALL lbc_lnk( zd1d2, 'T', -1. )
(  128) 
(  129)       DO ji = 2, jpi
(  130)          zd2d1(ji,:) = e2u(ji,:) - e2u(ji-1,:)
(  131)       END DO
(  132)       CALL lbc_lnk( zd2d1, 'T', -1. )
(  133) 
(  134)       akappa(:,:,1,1) =        1.0 / ( 2.0 * e1t(:,:) )
(  135)       akappa(:,:,1,2) = zd1d2(:,:) / ( 4.0 * e1t(:,:) * e2t(:,:) )
(  136)       akappa(:,:,2,1) = zd2d1(:,:) / ( 4.0 * e1t(:,:) * e2t(:,:) )
(  137)       akappa(:,:,2,2) =        1.0 / ( 2.0 * e2t(:,:) )
(  138)       
(  139)       !                                                      ! weights (wght)
(  140)       DO jj = 2, jpj
(  141)          DO ji = 2, jpi
(  142)             zusden = 1. / (  ( e1t(ji,jj) + e1t(ji-1,jj  ) )   &
(  143)                &           * ( e2t(ji,jj) + e2t(ji  ,jj-1) ) )
(  144)             wght(ji,jj,1,1) = zusden * e1t(ji  ,jj) * e2t(ji,jj  )
(  145)             wght(ji,jj,1,2) = zusden * e1t(ji  ,jj) * e2t(ji,jj-1)
(  146)             wght(ji,jj,2,1) = zusden * e1t(ji-1,jj) * e2t(ji,jj  )
(  147)             wght(ji,jj,2,2) = zusden * e1t(ji-1,jj) * e2t(ji,jj-1)
(  148)          END DO
(  149)       END DO
(  150)       CALL lbc_lnk( wght(:,:,1,1), 'I', 1. )      ! CAUTION: even with the lbc_lnk at ice U-V-point
(  151)       CALL lbc_lnk( wght(:,:,1,2), 'I', 1. )      ! the value of wght at jpj is wrong
(  152)       CALL lbc_lnk( wght(:,:,2,1), 'I', 1. )      ! but it is never used
(  153)       CALL lbc_lnk( wght(:,:,2,2), 'I', 1. )
(  154)     
(  155)       ! Coefficients for divergence of the stress tensor
(  156)       !-------------------------------------------------
(  157) 
(  158)       DO jj = 2, jpj
(  159)          DO ji = 2, jpi
(  160)             zh1p  =  e1t(ji  ,jj  ) * wght(ji,jj,2,2)   &
(  161)                &   + e1t(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  162)                &   + e1t(ji  ,jj-1) * wght(ji,jj,2,1)   &






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 4

(  163)                &   + e1t(ji-1,jj-1) * wght(ji,jj,1,1)
(  164) 
(  165)             zh2p  =  e2t(ji  ,jj  ) * wght(ji,jj,2,2)   &
(  166)                &   + e2t(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  167)                &   + e2t(ji  ,jj-1) * wght(ji,jj,2,1)   &
(  168)                &   + e2t(ji-1,jj-1) * wght(ji,jj,1,1)
(  169) 
(  170) ! better written but change the last digit and thus solver in less than 100 timestep
(  171) !           zh1p  = e1t(ji-1,jj  ) * wght(ji,jj,1,2) + e1t(ji,jj  ) * wght(ji,jj,2,2)   &
(  172) !              &  + e1t(ji-1,jj-1) * wght(ji,jj,1,1) + e1t(ji,jj-1) * wght(ji,jj,2,1) 
(  173) 
(  174) !           zh2p  = e2t(ji-1,jj  ) * wght(ji,jj,1,2) + e2t(ji,jj  ) * wght(ji,jj,2,2)   &
(  175) !              &  + e2t(ji-1,jj-1) * wght(ji,jj,1,1) + e2t(ji,jj-1) * wght(ji,jj,2,1)
(  176) 
(  177) !!ibug =0   zusden = 1.0 / ( zh1p * zh2p * 4.e0 )
(  178)             zusden = 1.0 / MAX( zh1p * zh2p * 4.e0 , 1.e-20 )
(  179)             zusden2 = zusden * 2.0 
(  180) 
(  181)             zd1d2p = zusden * 0.5 * ( -e1t(ji-1,jj-1) + e1t(ji-1,jj  ) - e1t(ji,jj-1) + e1t(ji  ,jj)   )
(  182)             zd2d1p = zusden * 0.5 * (  e2t(ji  ,jj-1) - e2t(ji-1,jj-1) + e2t(ji,jj  ) - e2t(ji-1,jj)   )
(  183) 
(  184)             alambd(ji,jj,2,2,2,1) = zusden2 * e2t(ji  ,jj-1)
(  185)             alambd(ji,jj,2,2,2,2) = zusden2 * e2t(ji  ,jj  )
(  186)             alambd(ji,jj,2,2,1,1) = zusden2 * e2t(ji-1,jj-1)
(  187)             alambd(ji,jj,2,2,1,2) = zusden2 * e2t(ji-1,jj  )
(  188) 
(  189)             alambd(ji,jj,1,1,2,1) = zusden2 * e1t(ji  ,jj-1)
(  190)             alambd(ji,jj,1,1,2,2) = zusden2 * e1t(ji  ,jj  )
(  191)             alambd(ji,jj,1,1,1,1) = zusden2 * e1t(ji-1,jj-1)
(  192)             alambd(ji,jj,1,1,1,2) = zusden2 * e1t(ji-1,jj  )
(  193) 
(  194)             alambd(ji,jj,1,2,2,1) = zd1d2p
(  195)             alambd(ji,jj,1,2,2,2) = zd1d2p
(  196)             alambd(ji,jj,1,2,1,1) = zd1d2p
(  197)             alambd(ji,jj,1,2,1,2) = zd1d2p
(  198) 
(  199)             alambd(ji,jj,2,1,2,1) = zd2d1p
(  200)             alambd(ji,jj,2,1,2,2) = zd2d1p
(  201)             alambd(ji,jj,2,1,1,1) = zd2d1p
(  202)             alambd(ji,jj,2,1,1,2) = zd2d1p
(  203)          END DO
(  204)       END DO
(  205) 
(  206)       CALL lbc_lnk( alambd(:,:,2,2,2,1), 'I', 1. )      ! CAUTION: even with the lbc_lnk at ice U-V point
(  207)       CALL lbc_lnk( alambd(:,:,2,2,2,2), 'I', 1. )      ! the value of wght at jpj is wrong
(  208)       CALL lbc_lnk( alambd(:,:,2,2,1,1), 'I', 1. )      ! but it is never used
(  209)       CALL lbc_lnk( alambd(:,:,2,2,1,2), 'I', 1. )      ! 
(  210) 
(  211)       CALL lbc_lnk( alambd(:,:,1,1,2,1), 'I', 1. )      ! CAUTION: idem
(  212)       CALL lbc_lnk( alambd(:,:,1,1,2,2), 'I', 1. )      ! 
(  213)       CALL lbc_lnk( alambd(:,:,1,1,1,1), 'I', 1. )      !
(  214)       CALL lbc_lnk( alambd(:,:,1,1,1,2), 'I', 1. )      !
(  215) 
(  216)       CALL lbc_lnk( alambd(:,:,1,2,2,1), 'I', 1. )      ! CAUTION: idem
(  217)       CALL lbc_lnk( alambd(:,:,1,2,2,2), 'I', 1. )      !
(  218)       CALL lbc_lnk( alambd(:,:,1,2,1,1), 'I', 1. )      !
(  219)       CALL lbc_lnk( alambd(:,:,1,2,1,2), 'I', 1. )      !
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 5

(  221)       CALL lbc_lnk( alambd(:,:,2,1,2,1), 'I', 1. )      ! CAUTION: idem
(  222)       CALL lbc_lnk( alambd(:,:,2,1,2,2), 'I', 1. )      !
(  223)       CALL lbc_lnk( alambd(:,:,2,1,1,1), 'I', 1. )      !
(  224)       CALL lbc_lnk( alambd(:,:,2,1,1,2), 'I', 1. )      !
(  225)             
(  226) 
(  227)       ! Initialization of ice masks
(  228)       !----------------------------
(  229)       
(  230)       tms(:,:) = tmask(:,:,1)      ! ice T-point  : use surface tmask
(  231) 
(  232) !i here we can use umask with a i and j shift of -1,-1
(  233)       tmu(:,1) = 0.e0
(  234)       tmu(1,:) = 0.e0
(  235)       DO jj = 2, jpj               ! ice U.V-point: computed from ice T-point mask
(  236)          DO ji = 2, jpim1
(  237)             tmu(ji,jj) =  tms(ji,jj) * tms(ji-1,jj) * tms(ji,jj-1) * tms(ji-1,jj-1)            
(  238)          END DO
(  239)       END DO
(  240)       
(  241)       !--lateral boundary conditions    
(  242)       CALL lbc_lnk( tmu(:,:), 'I', 1. )
(  243)       
(  244)       ! unmasked and masked area of T-grid cell
(  245)       area(:,:) = e1t(:,:) * e2t(:,:)
(  246)       
(  247)    END SUBROUTINE lim_msh
(  248) 
(  249) #else
(  250)    !!----------------------------------------------------------------------
(  251)    !!   Default option            Dummy Module         NO LIM sea-ice model
(  252)    !!----------------------------------------------------------------------
(  253) CONTAINS
(  254)    SUBROUTINE lim_msh           ! Dummy routine
(  255)    END SUBROUTINE lim_msh
(  256) #endif
(  257) 
(  258)    !!======================================================================
(  259) END MODULE limmsh

























PGF90 (Version     10.5)          09/23/2010  09:49:50      page 6

(    1) # 1 "limmsh.F90"
(    1) MODULE limmsh
(    2)    !!======================================================================
(    3)    !!                     ***  MODULE  limmsh  ***
(    4)    !! LIM ice model :   definition of the ice mesh parameters
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim'                                     LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_msh   : definition of the ice mesh
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE phycst
(   14)    USE dom_oce
(   15)    USE dom_ice
(   16)    USE lbclnk
(   17)    USE in_out_manager
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Accessibility
(   23)    PUBLIC lim_msh      ! routine called by ice_ini.F90
(   24) 
(   25)    !!----------------------------------------------------------------------
(   26)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limmsh.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32) 
(   33)    SUBROUTINE lim_msh
(   34)       !!-------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE lim_msh  ***
(   36)       !!              
(   37)       !! ** Purpose : Definition of the charact. of the numerical grid
(   38)       !!       
(   39)       !! ** Action  : - Initialisation of some variables
(   40)       !!              - Definition of some constants linked with the grid
(   41)       !!              - Definition of the metric coef. for the sea/ice
(   42)       !!              - Initialization of the ice masks (tmsk, umsk)
(   43)       !! 
(   44)       !! ** Refer.  : Deleersnijder et al. Ocean Modelling 100, 7-10 
(   45)       !!
(   46)       !! ** History :
(   47)       !!         original    : 01-04 (LIM)
(   48)       !!         addition    : 02-08 (C. Ethe, G. Madec)
(   49)       !!--------------------------------------------------------------------- 
(   50)       !! * Local variables
(   51)       INTEGER :: ji, jj      ! dummy loop indices
(   52) 
(   53)       REAL(wp), DIMENSION(jpi,jpj) ::  &
(   54)          zd2d1 , zd1d2       ! Derivative of zh2 (resp. zh1) in the x direction
(   55)          !                   ! (resp. y direction) (defined at the center)
(   56)       REAL(wp) ::         &
(   57)          zh1p  , zh2p   , &  ! Idem zh1, zh2 for the bottom left corner of the grid






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 7

(   58)          zd2d1p, zd1d2p , &  ! Idem zd2d1, zd1d2 for the bottom left corner of the grid
(   59)          zusden, zusden2     ! temporary scalars
(   60)       !!---------------------------------------------------------------------
(   61) 
(   62)       IF(lwp) THEN
(   63)          WRITE(numout,*)
(   64)          WRITE(numout,*) 'lim_msh : LIM sea-ice model, mesh initialization'
(   65)          WRITE(numout,*) '~~~~~~~'
(   66)       ENDIF
(   67)       
(   68)       !----------------------------------------------------------                          
(   69)       !    Initialization of local and some global (common) variables 
(   70)       !------------------------------------------------------------------ 
(   71)       
(   72)       njeq   = INT( jpj / 2 )   !i bug mpp potentiel
(   73)       njeqm1 = njeq - 1 
(   74) 
(   75)       fcor(:,:) = 2. * omega * SIN( gphit(:,:) * rad )   !  coriolis factor
(   76)  
(   77) !i    DO jj = 1, jpj
(   78) !i       zmsk(jj) = SUM( tmask(:,jj,:) )   ! = 0          if land  everywhere on a j-line
(   79) !!ii     write(numout,*) jj, zind(jj)
(   80) !i    END DO
(   81) 
(   82)       IF( fcor(1,1) * fcor(1,nlcj) < 0.e0 ) THEN   ! local domain include both hemisphere
(   83)          l_jeq = .TRUE.
(   84)          njeq  = 1
(   85)          DO WHILE ( njeq <= jpj .AND. fcor(1,njeq) < 0.e0 )
(   86)             njeq = njeq + 1
(   87)          END DO
(   88)          IF(lwp ) WRITE(numout,*) '          the equator is inside the domain at about njeq = ', njeq
(   89)       ELSEIF( fcor(1,1) < 0.e0 ) THEN
(   90)          l_jeq = .FALSE.
(   91)          njeq = jpj
(   92)          IF(lwp ) WRITE(numout,*) '          the model domain is entirely in the southern hemisphere: njeq = ', njeq
(   93)       ELSE
(   94)          l_jeq = .FALSE.
(   95)          njeq = 2
(   96)          IF(lwp ) WRITE(numout,*) '          the model domain is entirely in the northern hemisphere: njeq = ', njeq
(   97)       ENDIF
(   98) 
(   99)       njeqm1 = njeq - 1
(  100) 
(  101) 
(  102)       !   For each grid, definition of geometric tables 
(  103)       !------------------------------------------------------------------
(  104)       
(  105)       !-------------------
(  106)       ! Conventions :    !
(  107)       !-------------------
(  108)       !  indices 1 \ 2 <-> localisation in the 2 direction x \ y
(  109)       !  3rd indice <-> localisation on the mesh :
(  110)       !  0 = Centre ;  1 = corner W x(i-1/2) ; 2 = corner S y(j-1/2) ;
(  111)       !  3 = corner SW x(i-1/2),y(j-1/2)
(  112)       !-------------------
(  113) !!ibug ???
(  114)       akappa(:,:,:,:) = 0.e0
(  115)       wght(:,:,:,:) = 0.e0






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 8

(  116)       alambd(:,:,:,:,:,:) = 0.e0
(  117)       tmu(:,:) = 0.e0
(  118) !!i
(  119)       
(  120)       
(  121)       ! metric coefficients for sea ice dynamic
(  122)       !----------------------------------------
(  123)       !                                                       ! akappa
(  124)       DO jj = 2, jpj
(  125)          zd1d2(:,jj) = e1v(:,jj) - e1v(:,jj-1)
(  126)       END DO
(  127)       CALL lbc_lnk( zd1d2, 'T', -1. )
(  128) 
(  129)       DO ji = 2, jpi
(  130)          zd2d1(ji,:) = e2u(ji,:) - e2u(ji-1,:)
(  131)       END DO
(  132)       CALL lbc_lnk( zd2d1, 'T', -1. )
(  133) 
(  134)       akappa(:,:,1,1) =        1.0 / ( 2.0 * e1t(:,:) )
(  135)       akappa(:,:,1,2) = zd1d2(:,:) / ( 4.0 * e1t(:,:) * e2t(:,:) )
(  136)       akappa(:,:,2,1) = zd2d1(:,:) / ( 4.0 * e1t(:,:) * e2t(:,:) )
(  137)       akappa(:,:,2,2) =        1.0 / ( 2.0 * e2t(:,:) )
(  138)       
(  139)       !                                                      ! weights (wght)
(  140)       DO jj = 2, jpj
(  141)          DO ji = 2, jpi
(  142)             zusden = 1. / (  ( e1t(ji,jj) + e1t(ji-1,jj  ) )   &
(  143)                &           * ( e2t(ji,jj) + e2t(ji  ,jj-1) ) )
(  144)             wght(ji,jj,1,1) = zusden * e1t(ji  ,jj) * e2t(ji,jj  )
(  145)             wght(ji,jj,1,2) = zusden * e1t(ji  ,jj) * e2t(ji,jj-1)
(  146)             wght(ji,jj,2,1) = zusden * e1t(ji-1,jj) * e2t(ji,jj  )
(  147)             wght(ji,jj,2,2) = zusden * e1t(ji-1,jj) * e2t(ji,jj-1)
(  148)          END DO
(  149)       END DO
(  150)       CALL lbc_lnk( wght(:,:,1,1), 'I', 1. )      ! CAUTION: even with the lbc_lnk at ice U-V-point
(  151)       CALL lbc_lnk( wght(:,:,1,2), 'I', 1. )      ! the value of wght at jpj is wrong
(  152)       CALL lbc_lnk( wght(:,:,2,1), 'I', 1. )      ! but it is never used
(  153)       CALL lbc_lnk( wght(:,:,2,2), 'I', 1. )
(  154)     
(  155)       ! Coefficients for divergence of the stress tensor
(  156)       !-------------------------------------------------
(  157) 
(  158)       DO jj = 2, jpj
(  159)          DO ji = 2, jpi
(  160)             zh1p  =  e1t(ji  ,jj  ) * wght(ji,jj,2,2)   &
(  161)                &   + e1t(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  162)                &   + e1t(ji  ,jj-1) * wght(ji,jj,2,1)   &
(  163)                &   + e1t(ji-1,jj-1) * wght(ji,jj,1,1)
(  164) 
(  165)             zh2p  =  e2t(ji  ,jj  ) * wght(ji,jj,2,2)   &
(  166)                &   + e2t(ji-1,jj  ) * wght(ji,jj,1,2)   &
(  167)                &   + e2t(ji  ,jj-1) * wght(ji,jj,2,1)   &
(  168)                &   + e2t(ji-1,jj-1) * wght(ji,jj,1,1)
(  169) 
(  170) ! better written but change the last digit and thus solver in less than 100 timestep
(  171) !           zh1p  = e1t(ji-1,jj  ) * wght(ji,jj,1,2) + e1t(ji,jj  ) * wght(ji,jj,2,2)   &
(  172) !              &  + e1t(ji-1,jj-1) * wght(ji,jj,1,1) + e1t(ji,jj-1) * wght(ji,jj,2,1) 
(  173) 






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 9

(  174) !           zh2p  = e2t(ji-1,jj  ) * wght(ji,jj,1,2) + e2t(ji,jj  ) * wght(ji,jj,2,2)   &
(  175) !              &  + e2t(ji-1,jj-1) * wght(ji,jj,1,1) + e2t(ji,jj-1) * wght(ji,jj,2,1)
(  176) 
(  177) !!ibug =0   zusden = 1.0 / ( zh1p * zh2p * 4.e0 )
(  178)             zusden = 1.0 / MAX( zh1p * zh2p * 4.e0 , 1.e-20 )
(  179)             zusden2 = zusden * 2.0 
(  180) 
(  181)             zd1d2p = zusden * 0.5 * ( -e1t(ji-1,jj-1) + e1t(ji-1,jj  ) - e1t(ji,jj-1) + e1t(ji  ,jj)   )
(  182)             zd2d1p = zusden * 0.5 * (  e2t(ji  ,jj-1) - e2t(ji-1,jj-1) + e2t(ji,jj  ) - e2t(ji-1,jj)   )
(  183) 
(  184)             alambd(ji,jj,2,2,2,1) = zusden2 * e2t(ji  ,jj-1)
(  185)             alambd(ji,jj,2,2,2,2) = zusden2 * e2t(ji  ,jj  )
(  186)             alambd(ji,jj,2,2,1,1) = zusden2 * e2t(ji-1,jj-1)
(  187)             alambd(ji,jj,2,2,1,2) = zusden2 * e2t(ji-1,jj  )
(  188) 
(  189)             alambd(ji,jj,1,1,2,1) = zusden2 * e1t(ji  ,jj-1)
(  190)             alambd(ji,jj,1,1,2,2) = zusden2 * e1t(ji  ,jj  )
(  191)             alambd(ji,jj,1,1,1,1) = zusden2 * e1t(ji-1,jj-1)
(  192)             alambd(ji,jj,1,1,1,2) = zusden2 * e1t(ji-1,jj  )
(  193) 
(  194)             alambd(ji,jj,1,2,2,1) = zd1d2p
(  195)             alambd(ji,jj,1,2,2,2) = zd1d2p
(  196)             alambd(ji,jj,1,2,1,1) = zd1d2p
(  197)             alambd(ji,jj,1,2,1,2) = zd1d2p
(  198) 
(  199)             alambd(ji,jj,2,1,2,1) = zd2d1p
(  200)             alambd(ji,jj,2,1,2,2) = zd2d1p
(  201)             alambd(ji,jj,2,1,1,1) = zd2d1p
(  202)             alambd(ji,jj,2,1,1,2) = zd2d1p
(  203)          END DO
(  204)       END DO
(  205) 
(  206)       CALL lbc_lnk( alambd(:,:,2,2,2,1), 'I', 1. )      ! CAUTION: even with the lbc_lnk at ice U-V point
(  207)       CALL lbc_lnk( alambd(:,:,2,2,2,2), 'I', 1. )      ! the value of wght at jpj is wrong
(  208)       CALL lbc_lnk( alambd(:,:,2,2,1,1), 'I', 1. )      ! but it is never used
(  209)       CALL lbc_lnk( alambd(:,:,2,2,1,2), 'I', 1. )      ! 
(  210) 
(  211)       CALL lbc_lnk( alambd(:,:,1,1,2,1), 'I', 1. )      ! CAUTION: idem
(  212)       CALL lbc_lnk( alambd(:,:,1,1,2,2), 'I', 1. )      ! 
(  213)       CALL lbc_lnk( alambd(:,:,1,1,1,1), 'I', 1. )      !
(  214)       CALL lbc_lnk( alambd(:,:,1,1,1,2), 'I', 1. )      !
(  215) 
(  216)       CALL lbc_lnk( alambd(:,:,1,2,2,1), 'I', 1. )      ! CAUTION: idem
(  217)       CALL lbc_lnk( alambd(:,:,1,2,2,2), 'I', 1. )      !
(  218)       CALL lbc_lnk( alambd(:,:,1,2,1,1), 'I', 1. )      !
(  219)       CALL lbc_lnk( alambd(:,:,1,2,1,2), 'I', 1. )      !
(  220) 
(  221)       CALL lbc_lnk( alambd(:,:,2,1,2,1), 'I', 1. )      ! CAUTION: idem
(  222)       CALL lbc_lnk( alambd(:,:,2,1,2,2), 'I', 1. )      !
(  223)       CALL lbc_lnk( alambd(:,:,2,1,1,1), 'I', 1. )      !
(  224)       CALL lbc_lnk( alambd(:,:,2,1,1,2), 'I', 1. )      !
(  225)             
(  226) 
(  227)       ! Initialization of ice masks
(  228)       !----------------------------
(  229)       
(  230)       tms(:,:) = tmask(:,:,1)      ! ice T-point  : use surface tmask
(  231) 






PGF90 (Version     10.5)          09/23/2010  09:49:50      page 10

(  232) !i here we can use umask with a i and j shift of -1,-1
(  233)       tmu(:,1) = 0.e0
(  234)       tmu(1,:) = 0.e0
(  235)       DO jj = 2, jpj               ! ice U.V-point: computed from ice T-point mask
(  236)          DO ji = 2, jpim1
(  237)             tmu(ji,jj) =  tms(ji,jj) * tms(ji-1,jj) * tms(ji,jj-1) * tms(ji-1,jj-1)            
(  238)          END DO
(  239)       END DO
(  240)       
(  241)       !--lateral boundary conditions    
(  242)       CALL lbc_lnk( tmu(:,:), 'I', 1. )
(  243)       
(  244)       ! unmasked and masked area of T-grid cell
(  245)       area(:,:) = e1t(:,:) * e2t(:,:)
(  246)       
(  247)    END SUBROUTINE lim_msh
(  248) 
(  249) 
(  250) # 258
(  258)    !!======================================================================
(  259) END MODULE limmsh
