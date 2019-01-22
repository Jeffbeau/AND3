


PGF90 (Version     10.5)          09/23/2010  09:49:59      page 1

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

Filename: limtrp.F90

(    1) MODULE limtrp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limtrp   ***
(    4)    !! LIM transport ice model : sea-ice advection/diffusion
(    5)    !!======================================================================
(    6) #if defined key_ice_lim
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_trp      : advection/diffusion process of sea ice
(   11)    !!   lim_trp_init : initialization and namelist read
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst
(   15)    USE dom_oce
(   16)    USE daymod
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE ice_oce         ! ice variables
(   19)    USE dom_ice
(   20)    USE ice
(   21)    USE iceini
(   22)    USE limistate
(   23)    USE limadv
(   24)    USE limhdf
(   25)    USE lbclnk
(   26)    USE lib_mpp
(   27) 
(   28)    IMPLICIT NONE
(   29)    PRIVATE
(   30) 
(   31)    !! * Routine accessibility
(   32)    PUBLIC lim_trp       ! called by ice_step
(   33) 
(   34)    !! * Shared module variables
(   35)    REAL(wp), PUBLIC  ::   &  !:
(   36)       bound  = 0.e0          !: boundary condit. (0.0 no-slip, 1.0 free-slip)
(   37) 
(   38)    !! * Module variables
(   39)    REAL(wp)  ::           &  ! constant values
(   40)       epsi06 = 1.e-06  ,  &
(   41)       epsi03 = 1.e-03  ,  &
(   42)       epsi16 = 1.e-16  ,  &
(   43)       rzero  = 0.e0    ,  &
(   44)       rone   = 1.e0
(   45) 
(   46)    !! * Substitution






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 2

(   47) #  include "vectopt_loop_substitute.h90"
(   48)    !!----------------------------------------------------------------------
(   49)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   50)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limtrp.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
(   51)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   52)    !!----------------------------------------------------------------------
(   53) 
(   54) CONTAINS
(   55) 
(   56)    SUBROUTINE lim_trp
(   57)       !!-------------------------------------------------------------------
(   58)       !!                   ***  ROUTINE lim_trp ***
(   59)       !!                    
(   60)       !! ** purpose : advection/diffusion process of sea ice
(   61)       !!
(   62)       !! ** method  : variables included in the process are scalar,   
(   63)       !!     other values are considered as second order. 
(   64)       !!     For advection, a second order Prather scheme is used.  
(   65)       !!
(   66)       !! ** action :
(   67)       !!
(   68)       !! History :
(   69)       !!   1.0  !  00-01 (LIM)  Original code
(   70)       !!        !  01-05 (G. Madec, R. Hordoir) opa norm
(   71)       !!   2.0  !  04-01 (G. Madec, C. Ethe)  F90, mpp
(   72)       !!---------------------------------------------------------------------
(   73)       !! * Local Variables
(   74)       INTEGER  ::   ji, jj, jk,   &  ! dummy loop indices
(   75)          &          initad           ! number of sub-timestep for the advection
(   76) 
(   77)       REAL(wp) ::  &                              
(   78)          zindb  ,  &
(   79)          zacrith, &
(   80)          zindsn , &
(   81)          zindic , &
(   82)          zusvosn, &
(   83)          zusvoic, &
(   84)          zignm  , &
(   85)          zindhe , &
(   86)          zvbord , &
(   87)          zcfl   , &
(   88)          zusnit , &
(   89)          zrtt, ztsn, ztic1, ztic2
(   90) 
(   91)       REAL(wp), DIMENSION(jpi,jpj)  ::   &  ! temporary workspace
(   92)          zui_u , zvi_v , zsm   ,         &
(   93)          zs0ice, zs0sn , zs0a  ,         &
(   94)          zs0c0 , zs0c1 , zs0c2 ,         &
(   95)          zs0st
(   96)       !---------------------------------------------------------------------
(   97) 
(   98)       IF( numit == nstart  )   CALL lim_trp_init      ! Initialization (first time-step only)
(   99) 
(  100)       zsm(:,:) = area(:,:)
(  101)       
(  102)       IF( ln_limdyn ) THEN
(  103)          !-------------------------------------!
(  104)          !   Advection of sea ice properties   !






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 3

(  105)          !-------------------------------------!
(  106) 
(  107)          ! ice velocities at ocean U- and V-points (zui_u,zvi_v)
(  108)          ! ---------------------------------------
(  109)          ! zvbord factor between 1 and 2 to take into account slip or no-slip boundary conditions.        
(  110)          zvbord = 1.0 + ( 1.0 - bound )
(  111)          DO jj = 1, jpjm1
(  112)             DO ji = 1, jpim1
(  113)                zui_u(ji,jj) = ( u_ice(ji+1,jj  ) + u_ice(ji+1,jj+1) ) / ( MAX( tmu(ji+1,jj  ) + tmu(ji+1,jj+1), zvbord ) )
(  114)                zvi_v(ji,jj) = ( v_ice(ji  ,jj+1) + v_ice(ji+1,jj+1) ) / ( MAX( tmu(ji  ,jj+1) + tmu(ji+1,jj+1), zvbord ) )
(  115)             END DO
(  116)          END DO
(  117)          ! Lateral boundary conditions on zui_u, zvi_v
(  118)          CALL lbc_lnk( zui_u, 'U', -1. )
(  119)          CALL lbc_lnk( zvi_v, 'V', -1. )
(  120) 
(  121)          ! CFL test for stability
(  122)          ! ----------------------
(  123)          zcfl  = 0.e0
(  124)          zcfl  = MAX( zcfl, MAXVAL( ABS( zui_u(1:jpim1, :     ) ) * rdt_ice / e1u(1:jpim1, :     ) ) )
(  125)          zcfl  = MAX( zcfl, MAXVAL( ABS( zvi_v( :     ,1:jpjm1) ) * rdt_ice / e2v( :     ,1:jpjm1) ) )
(  126) 
(  127)          IF (lk_mpp ) CALL mpp_max(zcfl)
(  128) 
(  129)          IF ( zcfl > 0.5 .AND. lwp )   WRITE(numout,*) 'lim_trp : violation of cfl criterion the ',nday,'th day, cfl = ',zcf
(  130) 
(  131)          ! content of properties
(  132)          ! ---------------------
(  133)          zs0sn (:,:) =  hsnm(:,:) * area(:,:)                 ! Snow volume.
(  134)          zs0ice(:,:) =  hicm (:,:) * area(:,:)                ! Ice volume.
(  135)          zs0a  (:,:) =  ( 1.0 - frld(:,:) ) * area(:,:)       ! Surface covered by ice.
(  136)          zs0c0 (:,:) =  tbif(:,:,1) / rt0_snow * zs0sn(:,:)   ! Heat content of the snow layer.
(  137)          zs0c1 (:,:) =  tbif(:,:,2) / rt0_ice  * zs0ice(:,:)  ! Heat content of the first ice layer.
(  138)          zs0c2 (:,:) =  tbif(:,:,3) / rt0_ice  * zs0ice(:,:)  ! Heat content of the second ice layer.
(  139)          zs0st (:,:) =  qstoif(:,:) / xlic     * zs0a(:,:)    ! Heat reservoir for brine pockets.
(  140)          
(  141)  
(  142)          ! Advection 
(  143)          ! ---------
(  144)          ! If ice drift field is too fast, use an appropriate time step for advection.         
(  145)          initad = 1 + INT( MAX( rzero, SIGN( rone, zcfl-0.5 ) ) )
(  146)          zusnit = 1.0 / REAL( initad ) 
(  147)          
(  148)          IF ( MOD( nday , 2 ) == 0) THEN
(  149)             DO jk = 1,initad
(  150)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
(  151)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
(  152)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
(  153)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
(  154)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
(  155)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
(  156)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
(  157)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
(  158)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
(  159)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
(  160)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
(  161)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
(  162)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 4

(  163)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
(  164)             END DO
(  165)          ELSE
(  166)             DO jk = 1, initad
(  167)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
(  168)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
(  169)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
(  170)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
(  171)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
(  172)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
(  173)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
(  174)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
(  175)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
(  176)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
(  177)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
(  178)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
(  179)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
(  180)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
(  181)             END DO
(  182)          ENDIF
(  183)                         
(  184)          ! recover the properties from their contents
(  185)          ! ------------------------------------------
(  186)          zs0ice(:,:) = zs0ice(:,:) / area(:,:)
(  187)          zs0sn (:,:) = zs0sn (:,:) / area(:,:)
(  188)          zs0a  (:,:) = zs0a  (:,:) / area(:,:)
(  189)          zs0c0 (:,:) = zs0c0 (:,:) / area(:,:)
(  190)          zs0c1 (:,:) = zs0c1 (:,:) / area(:,:)
(  191)          zs0c2 (:,:) = zs0c2 (:,:) / area(:,:)
(  192)          zs0st (:,:) = zs0st (:,:) / area(:,:)
(  193) 
(  194) 
(  195)          !-------------------------------------!
(  196)          !   Diffusion of sea ice properties   !
(  197)          !-------------------------------------!
(  198) 
(  199)          ! Masked eddy diffusivity coefficient at ocean U- and V-points
(  200)          ! ------------------------------------------------------------
(  201)          DO jj = 1, jpjm1          ! NB: has not to be defined on jpj line and jpi row
(  202)             DO ji = 1 , fs_jpim1   ! vector opt.
(  203)                pahu(ji,jj) = ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji  ,jj) ) ) )   &
(  204)                   &        * ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji+1,jj) ) ) ) * ahiu(ji,jj)
(  205)                pahv(ji,jj) = ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji,jj  ) ) ) )   &
(  206)                   &        * ( 1.0 - MAX( rzero, SIGN( rone,- zs0a(ji,jj+1) ) ) ) * ahiv(ji,jj)
(  207)             END DO
(  208)          END DO
(  209) 
(  210)          ! diffusion
(  211)          ! ---------
(  212)          CALL lim_hdf( zs0ice )
(  213)          CALL lim_hdf( zs0sn  )
(  214)          CALL lim_hdf( zs0a   )
(  215)          CALL lim_hdf( zs0c0  )
(  216)          CALL lim_hdf( zs0c1  )
(  217)          CALL lim_hdf( zs0c2  )
(  218)          CALL lim_hdf( zs0st  )
(  219) 
(  220)          zs0ice(:,:) = MAX( rzero, zs0ice(:,:) * area(:,:) )    !!bug:  est-ce utile






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 5

(  221)          zs0sn (:,:) = MAX( rzero, zs0sn (:,:) * area(:,:) )    !!bug:  cf /area  juste apres
(  222)          zs0a  (:,:) = MAX( rzero, zs0a  (:,:) * area(:,:) )    !! suppression des 2 change le resultat...
(  223)          zs0c0 (:,:) = MAX( rzero, zs0c0 (:,:) * area(:,:) )
(  224)          zs0c1 (:,:) = MAX( rzero, zs0c1 (:,:) * area(:,:) )
(  225)          zs0c2 (:,:) = MAX( rzero, zs0c2 (:,:) * area(:,:) )
(  226)          zs0st (:,:) = MAX( rzero, zs0st (:,:) * area(:,:) )
(  227) 
(  228) 
(  229)          ! -------------------------------------------------------------------!
(  230)          !   Up-dating and limitation of sea ice properties after transport   !
(  231)          ! -------------------------------------------------------------------!
(  232) 
(  233)          ! Up-dating and limitation of sea ice properties after transport.
(  234)          DO jj = 1, jpj
(  235) !!!iii      zindhe = REAL( MAX( 0, isign(1, jj - njeq ) ) )              !ibug mpp  !!bugmpp  njeq!
(  236)             zindhe = MAX( 0.e0, SIGN( 1.e0, fcor(1,jj) ) )              ! = 0 for SH, =1 for NH
(  237)             DO ji = 1, jpi
(  238) 
(  239)                ! Recover mean values over the grid squares.
(  240)                zs0sn (ji,jj) = MAX( rzero, zs0sn (ji,jj)/area(ji,jj) )
(  241)                zs0ice(ji,jj) = MAX( rzero, zs0ice(ji,jj)/area(ji,jj) )
(  242)                zs0a  (ji,jj) = MAX( rzero, zs0a  (ji,jj)/area(ji,jj) )
(  243)                zs0c0 (ji,jj) = MAX( rzero, zs0c0 (ji,jj)/area(ji,jj) )
(  244)                zs0c1 (ji,jj) = MAX( rzero, zs0c1 (ji,jj)/area(ji,jj) )
(  245)                zs0c2 (ji,jj) = MAX( rzero, zs0c2 (ji,jj)/area(ji,jj) )
(  246)                zs0st (ji,jj) = MAX( rzero, zs0st (ji,jj)/area(ji,jj) )
(  247) 
(  248)                ! Recover in situ values.
(  249)                zindb         = MAX( rzero, SIGN( rone, zs0a(ji,jj) - epsi06 ) )
(  250)                zacrith       = 1.0 - ( zindhe * acrit(1) + ( 1.0 - zindhe ) * acrit(2) )
(  251)                zs0a (ji,jj)  = zindb * MIN( zs0a(ji,jj), zacrith )
(  252)                hsnif(ji,jj)  = zindb * ( zs0sn(ji,jj) /MAX( zs0a(ji,jj), epsi16 ) )
(  253)                hicif(ji,jj)  = zindb * ( zs0ice(ji,jj)/MAX( zs0a(ji,jj), epsi16 ) )
(  254)                zindsn        = MAX( rzero, SIGN( rone, hsnif(ji,jj) - epsi06 ) )
(  255)                zindic        = MAX( rzero, SIGN( rone, hicif(ji,jj) - epsi03 ) )
(  256)                zindb         = MAX( zindsn, zindic )
(  257)                zs0a (ji,jj)  = zindb * zs0a(ji,jj)
(  258)                frld (ji,jj)  = 1.0 - zs0a(ji,jj)
(  259)                hsnif(ji,jj)  = zindsn * hsnif(ji,jj)
(  260)                hicif(ji,jj)  = zindic * hicif(ji,jj)
(  261)                zusvosn       = 1.0/MAX( hsnif(ji,jj) * zs0a(ji,jj), epsi16 )
(  262)                zusvoic       = 1.0/MAX( hicif(ji,jj) * zs0a(ji,jj), epsi16 )
(  263)                zignm         = MAX( rzero,  SIGN( rone, hsndif - hsnif(ji,jj) ) )
(  264)                zrtt          = 173.15 * rone 
(  265)                ztsn          =          zignm   * tbif(ji,jj,1)  &
(  266)                               + ( 1.0 - zignm ) * MIN( MAX( zrtt, rt0_snow * zusvosn * zs0c0(ji,jj)) , tfu(ji,jj) ) 
(  267)                ztic1          = MIN( MAX( zrtt, rt0_ice * zusvoic * zs0c1(ji,jj) ) , tfu(ji,jj) )
(  268)                ztic2          = MIN( MAX( zrtt, rt0_ice * zusvoic * zs0c2(ji,jj) ) , tfu(ji,jj) )
(  269)  
(  270)                tbif(ji,jj,1) = zindsn * ztsn  + ( 1.0 - zindsn ) * tfu(ji,jj)               
(  271)                tbif(ji,jj,2) = zindic * ztic1 + ( 1.0 - zindic ) * tfu(ji,jj)
(  272)                tbif(ji,jj,3) = zindic * ztic2 + ( 1.0 - zindic ) * tfu(ji,jj)
(  273)                qstoif(ji,jj) = zindb  * xlic * zs0st(ji,jj) /  MAX( zs0a(ji,jj), epsi16 )
(  274)             END DO
(  275)          END DO
(  276)          
(  277)       ENDIF
(  278)       






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 6

(  279)    END SUBROUTINE lim_trp
(  280) 
(  281) 
(  282)    SUBROUTINE lim_trp_init
(  283)       !!-------------------------------------------------------------------
(  284)       !!                  ***  ROUTINE lim_trp_init  ***
(  285)       !!
(  286)       !! ** Purpose :   initialization of ice advection parameters
(  287)       !!
(  288)       !! ** Method  : Read the namicetrp namelist and check the parameter 
(  289)       !!       values called at the first timestep (nit000)
(  290)       !!
(  291)       !! ** input   :   Namelist namicetrp
(  292)       !!
(  293)       !! history :
(  294)       !!   2.0  !  03-08 (C. Ethe)  Original code
(  295)       !!-------------------------------------------------------------------
(  296)       NAMELIST/namicetrp/ bound
(  297)       !!-------------------------------------------------------------------
(  298) 
(  299)       ! Read Namelist namicetrp
(  300)       REWIND ( numnam_ice )
(  301)       READ   ( numnam_ice  , namicetrp )
(  302)       IF(lwp) THEN
(  303)          WRITE(numout,*)
(  304)          WRITE(numout,*) 'lim_trp_init : Ice parameters for advection '
(  305)          WRITE(numout,*) '~~~~~~~~~~~~'
(  306)          WRITE(numout,*) '   boundary conditions (0. no-slip, 1. free-slip) bound  = ', bound
(  307)       ENDIF
(  308)             
(  309)    END SUBROUTINE lim_trp_init
(  310) 
(  311) #else
(  312)    !!----------------------------------------------------------------------
(  313)    !!   Default option         Empty Module                No sea-ice model
(  314)    !!----------------------------------------------------------------------
(  315) CONTAINS
(  316)    SUBROUTINE lim_trp        ! Empty routine
(  317)    END SUBROUTINE lim_trp
(  318) #endif
(  319) 
(  320)    !!======================================================================
(  321) END MODULE limtrp





















PGF90 (Version     10.5)          09/23/2010  09:49:59      page 7

(    1) # 1 "limtrp.F90"
(    1) MODULE limtrp
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE limtrp   ***
(    4)    !! LIM transport ice model : sea-ice advection/diffusion
(    5)    !!======================================================================
(    6) # 7
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ice_lim' :                                   LIM sea-ice model
(    9)    !!----------------------------------------------------------------------
(   10)    !!   lim_trp      : advection/diffusion process of sea ice
(   11)    !!   lim_trp_init : initialization and namelist read
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE phycst
(   15)    USE dom_oce
(   16)    USE daymod
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE ice_oce         ! ice variables
(   19)    USE dom_ice
(   20)    USE ice
(   21)    USE iceini
(   22)    USE limistate
(   23)    USE limadv
(   24)    USE limhdf
(   25)    USE lbclnk
(   26)    USE lib_mpp
(   27) 
(   28)    IMPLICIT NONE
(   29)    PRIVATE
(   30) 
(   31)    !! * Routine accessibility
(   32)    PUBLIC lim_trp       ! called by ice_step
(   33) 
(   34)    !! * Shared module variables
(   35)    REAL(wp), PUBLIC  ::   &  !:
(   36)       bound  = 0.e0          !: boundary condit. (0.0 no-slip, 1.0 free-slip)
(   37) 
(   38)    !! * Module variables
(   39)    REAL(wp)  ::           &  ! constant values
(   40)       epsi06 = 1.e-06  ,  &
(   41)       epsi03 = 1.e-03  ,  &
(   42)       epsi16 = 1.e-16  ,  &
(   43)       rzero  = 0.e0    ,  &
(   44)       rone   = 1.e0
(   45) 
(   46)    !! * Substitution
(   47) # 1 "./vectopt_loop_substitute.h90"
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






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 8

(   11)    !!----------------------------------------------------------------------
(   12) # 48 "limtrp.F90"
(   48) # 48
(   48)    !!----------------------------------------------------------------------
(   49)    !!   LIM 2.0,  UCL-LOCEAN-IPSL (2005) 
(   50)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/LIM_SRC/limtrp.F90,v 1.5 2005/03/27 18:34:42 opalod Exp $ 
(   51)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   52)    !!----------------------------------------------------------------------
(   53) 
(   54) CONTAINS
(   55) 
(   56)    SUBROUTINE lim_trp
(   57)       !!-------------------------------------------------------------------
(   58)       !!                   ***  ROUTINE lim_trp ***
(   59)       !!                    
(   60)       !! ** purpose : advection/diffusion process of sea ice
(   61)       !!
(   62)       !! ** method  : variables included in the process are scalar,   
(   63)       !!     other values are considered as second order. 
(   64)       !!     For advection, a second order Prather scheme is used.  
(   65)       !!
(   66)       !! ** action :
(   67)       !!
(   68)       !! History :
(   69)       !!   1.0  !  00-01 (LIM)  Original code
(   70)       !!        !  01-05 (G. Madec, R. Hordoir) opa norm
(   71)       !!   2.0  !  04-01 (G. Madec, C. Ethe)  F90, mpp
(   72)       !!---------------------------------------------------------------------
(   73)       !! * Local Variables
(   74)       INTEGER  ::   ji, jj, jk,   &  ! dummy loop indices
(   75)          &          initad           ! number of sub-timestep for the advection
(   76) 
(   77)       REAL(wp) ::  &                              
(   78)          zindb  ,  &
(   79)          zacrith, &
(   80)          zindsn , &
(   81)          zindic , &
(   82)          zusvosn, &
(   83)          zusvoic, &
(   84)          zignm  , &
(   85)          zindhe , &
(   86)          zvbord , &
(   87)          zcfl   , &
(   88)          zusnit , &
(   89)          zrtt, ztsn, ztic1, ztic2
(   90) 
(   91)       REAL(wp), DIMENSION(jpi,jpj)  ::   &  ! temporary workspace
(   92)          zui_u , zvi_v , zsm   ,         &
(   93)          zs0ice, zs0sn , zs0a  ,         &
(   94)          zs0c0 , zs0c1 , zs0c2 ,         &
(   95)          zs0st
(   96)       !---------------------------------------------------------------------
(   97) 
(   98)       IF( numit == nstart  )   CALL lim_trp_init      ! Initialization (first time-step only)
(   99) 
(  100)       zsm(:,:) = area(:,:)
(  101)       
(  102)       IF( ln_limdyn ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 9

(  103)          !-------------------------------------!
(  104)          !   Advection of sea ice properties   !
(  105)          !-------------------------------------!
(  106) 
(  107)          ! ice velocities at ocean U- and V-points (zui_u,zvi_v)
(  108)          ! ---------------------------------------
(  109)          ! zvbord factor between 1 and 2 to take into account slip or no-slip boundary conditions.        
(  110)          zvbord = 1.0 + ( 1.0 - bound )
(  111)          DO jj = 1, jpjm1
(  112)             DO ji = 1, jpim1
(  113)                zui_u(ji,jj) = ( u_ice(ji+1,jj  ) + u_ice(ji+1,jj+1) ) / ( MAX( tmu(ji+1,jj  ) + tmu(ji+1,jj+1), zvbord ) )
(  114)                zvi_v(ji,jj) = ( v_ice(ji  ,jj+1) + v_ice(ji+1,jj+1) ) / ( MAX( tmu(ji  ,jj+1) + tmu(ji+1,jj+1), zvbord ) )
(  115)             END DO
(  116)          END DO
(  117)          ! Lateral boundary conditions on zui_u, zvi_v
(  118)          CALL lbc_lnk( zui_u, 'U', -1. )
(  119)          CALL lbc_lnk( zvi_v, 'V', -1. )
(  120) 
(  121)          ! CFL test for stability
(  122)          ! ----------------------
(  123)          zcfl  = 0.e0
(  124)          zcfl  = MAX( zcfl, MAXVAL( ABS( zui_u(1:jpim1, :     ) ) * rdt_ice / e1u(1:jpim1, :     ) ) )
(  125)          zcfl  = MAX( zcfl, MAXVAL( ABS( zvi_v( :     ,1:jpjm1) ) * rdt_ice / e2v( :     ,1:jpjm1) ) )
(  126) 
(  127)          IF (lk_mpp ) CALL mpp_max(zcfl)
(  128) 
(  129)          IF ( zcfl > 0.5 .AND. lwp )   WRITE(numout,*) 'lim_trp : violation of cfl criterion the ',nday,'th day, cfl = ',zcfl
(  130) 
(  131)          ! content of properties
(  132)          ! ---------------------
(  133)          zs0sn (:,:) =  hsnm(:,:) * area(:,:)                 ! Snow volume.
(  134)          zs0ice(:,:) =  hicm (:,:) * area(:,:)                ! Ice volume.
(  135)          zs0a  (:,:) =  ( 1.0 - frld(:,:) ) * area(:,:)       ! Surface covered by ice.
(  136)          zs0c0 (:,:) =  tbif(:,:,1) / rt0_snow * zs0sn(:,:)   ! Heat content of the snow layer.
(  137)          zs0c1 (:,:) =  tbif(:,:,2) / rt0_ice  * zs0ice(:,:)  ! Heat content of the first ice layer.
(  138)          zs0c2 (:,:) =  tbif(:,:,3) / rt0_ice  * zs0ice(:,:)  ! Heat content of the second ice layer.
(  139)          zs0st (:,:) =  qstoif(:,:) / xlic     * zs0a(:,:)    ! Heat reservoir for brine pockets.
(  140)          
(  141)  
(  142)          ! Advection 
(  143)          ! ---------
(  144)          ! If ice drift field is too fast, use an appropriate time step for advection.         
(  145)          initad = 1 + INT( MAX( rzero, SIGN( rone, zcfl-0.5 ) ) )
(  146)          zusnit = 1.0 / REAL( initad ) 
(  147)          
(  148)          IF ( MOD( nday , 2 ) == 0) THEN
(  149)             DO jk = 1,initad
(  150)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
(  151)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
(  152)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
(  153)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
(  154)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
(  155)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
(  156)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
(  157)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
(  158)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
(  159)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
(  160)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 10

(  161)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
(  162)                CALL lim_adv_x( zusnit, zui_u, rone , zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
(  163)                CALL lim_adv_y( zusnit, zvi_v, rzero, zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
(  164)             END DO
(  165)          ELSE
(  166)             DO jk = 1, initad
(  167)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
(  168)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0ice, sxice, sxxice, syice, syyice, sxyice )
(  169)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
(  170)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0sn , sxsn , sxxsn , sysn , syysn , sxysn  )
(  171)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
(  172)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0a  , sxa  , sxxa  , sya  , syya  , sxya   )
(  173)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
(  174)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c0 , sxc0 , sxxc0 , syc0 , syyc0 , sxyc0  )
(  175)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
(  176)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c1 , sxc1 , sxxc1 , syc1 , syyc1 , sxyc1  )
(  177)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
(  178)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0c2 , sxc2 , sxxc2 , syc2 , syyc2 , sxyc2  )
(  179)                CALL lim_adv_y( zusnit, zvi_v, rone , zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
(  180)                CALL lim_adv_x( zusnit, zui_u, rzero, zsm, zs0st , sxst , sxxst , syst , syyst , sxyst  )
(  181)             END DO
(  182)          ENDIF
(  183)                         
(  184)          ! recover the properties from their contents
(  185)          ! ------------------------------------------
(  186)          zs0ice(:,:) = zs0ice(:,:) / area(:,:)
(  187)          zs0sn (:,:) = zs0sn (:,:) / area(:,:)
(  188)          zs0a  (:,:) = zs0a  (:,:) / area(:,:)
(  189)          zs0c0 (:,:) = zs0c0 (:,:) / area(:,:)
(  190)          zs0c1 (:,:) = zs0c1 (:,:) / area(:,:)
(  191)          zs0c2 (:,:) = zs0c2 (:,:) / area(:,:)
(  192)          zs0st (:,:) = zs0st (:,:) / area(:,:)
(  193) 
(  194) 
(  195)          !-------------------------------------!
(  196)          !   Diffusion of sea ice properties   !
(  197)          !-------------------------------------!
(  198) 
(  199)          ! Masked eddy diffusivity coefficient at ocean U- and V-points
(  200)          ! ------------------------------------------------------------
(  201)          DO jj = 1, jpjm1          ! NB: has not to be defined on jpj line and jpi row
(  202)             DO ji = 1 , jpim1   ! vector opt.
(  203)                pahu(ji,jj) = ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji  ,jj) ) ) )   &
(  204)                   &        * ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji+1,jj) ) ) ) * ahiu(ji,jj)
(  205)                pahv(ji,jj) = ( 1.0 - MAX( rzero, SIGN( rone, -zs0a(ji,jj  ) ) ) )   &
(  206)                   &        * ( 1.0 - MAX( rzero, SIGN( rone,- zs0a(ji,jj+1) ) ) ) * ahiv(ji,jj)
(  207)             END DO
(  208)          END DO
(  209) 
(  210)          ! diffusion
(  211)          ! ---------
(  212)          CALL lim_hdf( zs0ice )
(  213)          CALL lim_hdf( zs0sn  )
(  214)          CALL lim_hdf( zs0a   )
(  215)          CALL lim_hdf( zs0c0  )
(  216)          CALL lim_hdf( zs0c1  )
(  217)          CALL lim_hdf( zs0c2  )
(  218)          CALL lim_hdf( zs0st  )






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 11

(  219) 
(  220)          zs0ice(:,:) = MAX( rzero, zs0ice(:,:) * area(:,:) )    !!bug:  est-ce utile
(  221)          zs0sn (:,:) = MAX( rzero, zs0sn (:,:) * area(:,:) )    !!bug:  cf /area  juste apres
(  222)          zs0a  (:,:) = MAX( rzero, zs0a  (:,:) * area(:,:) )    !! suppression des 2 change le resultat...
(  223)          zs0c0 (:,:) = MAX( rzero, zs0c0 (:,:) * area(:,:) )
(  224)          zs0c1 (:,:) = MAX( rzero, zs0c1 (:,:) * area(:,:) )
(  225)          zs0c2 (:,:) = MAX( rzero, zs0c2 (:,:) * area(:,:) )
(  226)          zs0st (:,:) = MAX( rzero, zs0st (:,:) * area(:,:) )
(  227) 
(  228) 
(  229)          ! -------------------------------------------------------------------!
(  230)          !   Up-dating and limitation of sea ice properties after transport   !
(  231)          ! -------------------------------------------------------------------!
(  232) 
(  233)          ! Up-dating and limitation of sea ice properties after transport.
(  234)          DO jj = 1, jpj
(  235) !!!iii      zindhe = REAL( MAX( 0, isign(1, jj - njeq ) ) )              !ibug mpp  !!bugmpp  njeq!
(  236)             zindhe = MAX( 0.e0, SIGN( 1.e0, fcor(1,jj) ) )              ! = 0 for SH, =1 for NH
(  237)             DO ji = 1, jpi
(  238) 
(  239)                ! Recover mean values over the grid squares.
(  240)                zs0sn (ji,jj) = MAX( rzero, zs0sn (ji,jj)/area(ji,jj) )
(  241)                zs0ice(ji,jj) = MAX( rzero, zs0ice(ji,jj)/area(ji,jj) )
(  242)                zs0a  (ji,jj) = MAX( rzero, zs0a  (ji,jj)/area(ji,jj) )
(  243)                zs0c0 (ji,jj) = MAX( rzero, zs0c0 (ji,jj)/area(ji,jj) )
(  244)                zs0c1 (ji,jj) = MAX( rzero, zs0c1 (ji,jj)/area(ji,jj) )
(  245)                zs0c2 (ji,jj) = MAX( rzero, zs0c2 (ji,jj)/area(ji,jj) )
(  246)                zs0st (ji,jj) = MAX( rzero, zs0st (ji,jj)/area(ji,jj) )
(  247) 
(  248)                ! Recover in situ values.
(  249)                zindb         = MAX( rzero, SIGN( rone, zs0a(ji,jj) - epsi06 ) )
(  250)                zacrith       = 1.0 - ( zindhe * acrit(1) + ( 1.0 - zindhe ) * acrit(2) )
(  251)                zs0a (ji,jj)  = zindb * MIN( zs0a(ji,jj), zacrith )
(  252)                hsnif(ji,jj)  = zindb * ( zs0sn(ji,jj) /MAX( zs0a(ji,jj), epsi16 ) )
(  253)                hicif(ji,jj)  = zindb * ( zs0ice(ji,jj)/MAX( zs0a(ji,jj), epsi16 ) )
(  254)                zindsn        = MAX( rzero, SIGN( rone, hsnif(ji,jj) - epsi06 ) )
(  255)                zindic        = MAX( rzero, SIGN( rone, hicif(ji,jj) - epsi03 ) )
(  256)                zindb         = MAX( zindsn, zindic )
(  257)                zs0a (ji,jj)  = zindb * zs0a(ji,jj)
(  258)                frld (ji,jj)  = 1.0 - zs0a(ji,jj)
(  259)                hsnif(ji,jj)  = zindsn * hsnif(ji,jj)
(  260)                hicif(ji,jj)  = zindic * hicif(ji,jj)
(  261)                zusvosn       = 1.0/MAX( hsnif(ji,jj) * zs0a(ji,jj), epsi16 )
(  262)                zusvoic       = 1.0/MAX( hicif(ji,jj) * zs0a(ji,jj), epsi16 )
(  263)                zignm         = MAX( rzero,  SIGN( rone, hsndif - hsnif(ji,jj) ) )
(  264)                zrtt          = 173.15 * rone 
(  265)                ztsn          =          zignm   * tbif(ji,jj,1)  &
(  266)                               + ( 1.0 - zignm ) * MIN( MAX( zrtt, rt0_snow * zusvosn * zs0c0(ji,jj)) , tfu(ji,jj) ) 
(  267)                ztic1          = MIN( MAX( zrtt, rt0_ice * zusvoic * zs0c1(ji,jj) ) , tfu(ji,jj) )
(  268)                ztic2          = MIN( MAX( zrtt, rt0_ice * zusvoic * zs0c2(ji,jj) ) , tfu(ji,jj) )
(  269)  
(  270)                tbif(ji,jj,1) = zindsn * ztsn  + ( 1.0 - zindsn ) * tfu(ji,jj)               
(  271)                tbif(ji,jj,2) = zindic * ztic1 + ( 1.0 - zindic ) * tfu(ji,jj)
(  272)                tbif(ji,jj,3) = zindic * ztic2 + ( 1.0 - zindic ) * tfu(ji,jj)
(  273)                qstoif(ji,jj) = zindb  * xlic * zs0st(ji,jj) /  MAX( zs0a(ji,jj), epsi16 )
(  274)             END DO
(  275)          END DO
(  276)          






PGF90 (Version     10.5)          09/23/2010  09:49:59      page 12

(  277)       ENDIF
(  278)       
(  279)    END SUBROUTINE lim_trp
(  280) 
(  281) 
(  282)    SUBROUTINE lim_trp_init
(  283)       !!-------------------------------------------------------------------
(  284)       !!                  ***  ROUTINE lim_trp_init  ***
(  285)       !!
(  286)       !! ** Purpose :   initialization of ice advection parameters
(  287)       !!
(  288)       !! ** Method  : Read the namicetrp namelist and check the parameter 
(  289)       !!       values called at the first timestep (nit000)
(  290)       !!
(  291)       !! ** input   :   Namelist namicetrp
(  292)       !!
(  293)       !! history :
(  294)       !!   2.0  !  03-08 (C. Ethe)  Original code
(  295)       !!-------------------------------------------------------------------
(  296)       NAMELIST/namicetrp/ bound
(  297)       !!-------------------------------------------------------------------
(  298) 
(  299)       ! Read Namelist namicetrp
(  300)       REWIND ( numnam_ice )
(  301)       READ   ( numnam_ice  , namicetrp )
(  302)       IF(lwp) THEN
(  303)          WRITE(numout,*)
(  304)          WRITE(numout,*) 'lim_trp_init : Ice parameters for advection '
(  305)          WRITE(numout,*) '~~~~~~~~~~~~'
(  306)          WRITE(numout,*) '   boundary conditions (0. no-slip, 1. free-slip) bound  = ', bound
(  307)       ENDIF
(  308)             
(  309)    END SUBROUTINE lim_trp_init
(  310) 
(  311) 
(  312) # 320
(  320)    !!======================================================================
(  321) END MODULE limtrp
