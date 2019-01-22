


PGF90 (Version     10.5)          09/23/2010  09:49:07      page 1

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

Filename: zdfddm.F90

(    1) MODULE zdfddm
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfddm  ***
(    4)    !! Ocean physics : double diffusion mixing parameterization
(    5)    !!======================================================================
(    6) #if defined key_zdfddm   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_zdfddm' :                                     double diffusion
(    9)    !!----------------------------------------------------------------------
(   10)    !!   zdf_ddm       : compute the Ks for salinity
(   11)    !!   zdf_ddm_init  : read namelist and control the parameters
(   12)    !!----------------------------------------------------------------------
(   13)    !! * Modules used
(   14)    USE oce             ! ocean dynamics and tracers variables
(   15)    USE dom_oce         ! ocean space and time domain variables 
(   16)    USE zdf_oce         ! ocean vertical physics variables
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   19)    USE prtctl          ! Print control
(   20) 
(   21)    IMPLICIT NONE
(   22)    PRIVATE
(   23) 
(   24)    !! * Routine accessibility
(   25)    PUBLIC zdf_ddm     ! called by step.F90
(   26) 
(   27)    !! * Shared module variables
(   28)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdfddm = .TRUE.    !: double diffusive mixing flag
(   29)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj,jpk) ::   &   !:
(   30)       avs ,               &  !: salinity vertical diffusivity coeff. at w-point
(   31)       rrau                   !: heat/salt buoyancy flux ratio
(   32) 
(   33)    !! * Module variables
(   34)    REAL(wp) ::            & !!! * double diffusive mixing namelist *
(   35)       avts  = 1.e-4_wp ,  &  ! maximum value of avs for salt fingering
(   36)       hsbfr = 1.6_wp         ! heat/salt buoyancy flux ratio
(   37) 
(   38)    !! * Substitutions
(   39) #  include "vectopt_loop_substitute.h90"
(   40)    !!----------------------------------------------------------------------
(   41)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   42)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfddm.F90,v 1.6 2005/09/02 15:45:43 opalod Exp $ 
(   43)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   44)    !!----------------------------------------------------------------------
(   45) 
(   46) CONTAINS






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 2

(   47) 
(   48)    SUBROUTINE zdf_ddm( kt )
(   49)       !!----------------------------------------------------------------------
(   50)       !!                  ***  ROUTINE zdf_ddm  ***
(   51)       !!                    
(   52)       !! ** Purpose :   Add to the vertical eddy diffusivity coefficient the 
(   53)       !!      effect of salt fingering and diffusive convection. 
(   54)       !!
(   55)       !! ** Method  :   Diapycnal mixing is increased in case of double
(   56)       !!      diffusive mixing (i.e. salt fingering and diffusive layering)
(   57)       !!      following Merryfield et al. (1999). The rate of double diffusive 
(   58)       !!      mixing depend on the buoyancy ratio: Rrau=alpha/beta dk[T]/dk[S]
(   59)       !!      which is computed in rn2.F
(   60)       !!         * salt fingering (Schmitt 1981):
(   61)       !!      for Rrau > 1 and rn2 > 0 : zavfs = avts / ( 1 + (Rrau/hsbfr)^6 )
(   62)       !!      for Rrau > 1 and rn2 > 0 : zavfs = O
(   63)       !!      otherwise                : zavft = 0.7 zavs / Rrau
(   64)       !!         * diffusive layering (Federov 1988):
(   65)       !!      for 0< Rrau < 1 and rn2 > 0 : zavdt = 1.3635e-6  
(   66)       !!                                 * exp( 4.6 exp(-0.54 (1/Rrau-1) ) )
(   67)       !!      otherwise                   : zavdt = 0 
(   68)       !!      for .5 < Rrau < 1 and rn2 > 0 : zavds = zavdt (1.885 Rrau -0.85)
(   69)       !!      for  0 < Rrau <.5 and rn2 > 0 : zavds = zavdt 0.15 Rrau      
(   70)       !!      otherwise                     : zavds = 0 
(   71)       !!         * update the eddy diffusivity:
(   72)       !!      avt = avt + zavft + zavdt
(   73)       !!      avs = avs + zavfs + zavds
(   74)       !!      avmu, avmv are required to remain at least above avt and avs.
(   75)       !!      
(   76)       !! ** Action  :   avt, avs : update vertical eddy diffusivity coef.
(   77)       !!                           for temperature and salinity
(   78)       !!
(   79)       !! References :
(   80)       !!      Merryfield et al., JPO, 29, 1124-1142, 1999.
(   81)       !! History :
(   82)       !!        !  00-08  (G. Madec)  double diffusive mixing
(   83)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   84)       !!----------------------------------------------------------------------
(   85)       !! * Arguments
(   86)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step indexocean time step
(   87) 
(   88)       !! * Local declarations
(   89)       INTEGER ::   ji, jj , jk              ! dummy loop indices
(   90)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   91)          zmsks, zmskf,                    & ! temporary workspace 
(   92)          zmskd1, zmskd2, zmskd3             !    "           "
(   93)       REAL(wp) ::   &
(   94)          zinr, zrr,                       & ! temporary scalars
(   95)          zavft, zavfs,                    & !    "         "
(   96)          zavdt, zavds                       !    "         "
(   97)       !!----------------------------------------------------------------------
(   98) 
(   99) 
(  100)       IF ( kt == nit000 )   CALL zdf_ddm_init          ! Initialization (first time-step only)
(  101) 
(  102) 
(  103)       ! Compute avs
(  104)       ! -----------






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 3

(  105)       !                                                ! ===============
(  106)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  107)          !                                             ! ===============
(  108)          ! Define the mask 
(  109)          ! ---------------
(  110)          ! only retains positive value of rrau
(  111)          rrau(:,:,jk) = MAX( 1.e-20, rrau(:,:,jk) )
(  112) 
(  113)          ! indicators:
(  114)          DO jj = 1, jpj
(  115)             DO ji = 1, jpi
(  116)                ! stability indicator: msks=1 if rn2>0; 0 elsewhere
(  117)                IF( rn2(ji,jj,jk) + 1.e-12  <= 0. ) THEN
(  118)                   zmsks(ji,jj) = 0.e0
(  119)                ELSE
(  120)                   zmsks(ji,jj) = 1.e0
(  121)                ENDIF
(  122)                ! salt fingering indicator: msksf=1 if rrau>1; 0 elsewhere            
(  123)                IF( rrau(ji,jj,jk) <= 1. ) THEN
(  124)                   zmskf(ji,jj) = 0.e0
(  125)                ELSE
(  126)                   zmskf(ji,jj) = 1.e0
(  127)                ENDIF
(  128)                ! diffusive layering indicators: 
(  129)                !   mskdl1=1 if 0<rrau<1; 0 elsewhere
(  130)                IF( rrau(ji,jj,jk) >= 1. ) THEN
(  131)                   zmskd1(ji,jj) = 0.e0
(  132)                ELSE
(  133)                   zmskd1(ji,jj) = 1.e0
(  134)                ENDIF
(  135)                !   mskdl2=1 if 0<rrau<0.5; 0 elsewhere
(  136)                IF( rrau(ji,jj,jk) >= 0.5 ) THEN
(  137)                   zmskd2(ji,jj) = 0.e0
(  138)                ELSE
(  139)                   zmskd2(ji,jj) = 1.e0
(  140)                ENDIF
(  141)                !   mskdl3=1 if 0.5<rrau<1; 0 elsewhere
(  142)                IF( rrau(ji,jj,jk) <= 0.5 .OR. rrau(ji,jj,jk) >= 1. ) THEN
(  143)                   zmskd3(ji,jj) = 0.e0
(  144)                ELSE
(  145)                   zmskd3(ji,jj) = 1.e0
(  146)                ENDIF
(  147)             END DO
(  148)          END DO
(  149)          ! mask zmsk in order to have avt and avs masked
(  150)          zmsks(:,:) = zmsks(:,:) * tmask(:,:,jk)
(  151) 
(  152) 
(  153)          ! Update avt and avs
(  154)          ! ------------------
(  155)          ! Constant eddy coefficient: reset to the background value
(  156) !CDIR NOVERRCHK
(  157)          DO jj = 1, jpj
(  158) !CDIR NOVERRCHK
(  159)             DO ji = 1, jpi
(  160)                zinr = 1./rrau(ji,jj,jk)
(  161)                ! salt fingering
(  162)                zrr = rrau(ji,jj,jk)/hsbfr






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 4

(  163)                zrr = zrr * zrr
(  164)                zavfs = avts / ( 1 + zrr*zrr*zrr ) * zmsks(ji,jj) *zmskf(ji,jj)
(  165)                zavft = 0.7 * zavfs / rrau(ji,jj,jk)
(  166)                ! diffusive layering
(  167)                zavdt = 1.3635e-6 * EXP(4.6*EXP(-0.54*(zinr-1.) ) )   &
(  168)                                  * zmsks(ji,jj) * zmskd1(ji,jj)
(  169)                zavds = zavdt * zmsks(ji,jj)   &
(  170)                      * ( (1.85 * rrau(ji,jj,jk) - 0.85 ) * zmskd3(ji,jj)   &
(  171)                          + zavdt * 0.15 * rrau(ji,jj,jk) * zmskd2(ji,jj)  )
(  172)                ! add to the eddy viscosity coef. previously computed
(  173)                avs (ji,jj,jk) = avt(ji,jj,jk) + zavfs + zavds
(  174)                avt (ji,jj,jk) = avt(ji,jj,jk) + zavft + zavdt
(  175)             END DO
(  176)          END DO
(  177) 
(  178) 
(  179)          ! Increase avmu, avmv if necessary
(  180)          ! --------------------------------
(  181)          DO jj = 1, jpjm1
(  182)             DO ji = 1, fs_jpim1   ! vector opt.
(  183)                avmu(ji,jj,jk) = MAX( avmu(ji,jj,jk),    &
(  184)                                      avt(ji,jj,jk), avt(ji+1,jj,jk),   &
(  185)                                      avs(ji,jj,jk), avs(ji+1,jj,jk) )   &
(  186)                               * umask(ji,jj,jk)
(  187)                avmv(ji,jj,jk) = MAX( avmv(ji,jj,jk),    &
(  188)                                      avt(ji,jj,jk), avt(ji,jj+1,jk),   &
(  189)                                      avs(ji,jj,jk), avs(ji,jj+1,jk) )   &
(  190)                               * vmask(ji,jj,jk)
(  191)             END DO
(  192)          END DO
(  193)          !                                                ! ===============
(  194)       END DO                                              !   End of slab
(  195)       !                                                   ! ===============
(  196)       
(  197)       ! Lateral boundary conditions on ( avt, avs, avmu, avmv )   (unchanged sign)
(  198)       ! -------------------------------========================
(  199)       CALL lbc_lnk( avt , 'W', 1. )
(  200)       CALL lbc_lnk( avs , 'W', 1. )
(  201)       CALL lbc_lnk( avmu, 'U', 1. ) 
(  202)       CALL lbc_lnk( avmv, 'V', 1. )
(  203) 
(  204)       IF(ln_ctl) THEN
(  205)          CALL prt_ctl(tab3d_1=avt , clinfo1=' ddm  - t: ', tab3d_2=avs , clinfo2=' s: ', ovlap=1, kdim=jpk)
(  206)          CALL prt_ctl(tab3d_1=avmu, clinfo1=' ddm  - u: ', tab3d_2=avmv, clinfo2=' v: ', ovlap=1, kdim=jpk)
(  207)       ENDIF
(  208)       
(  209)    END SUBROUTINE zdf_ddm
(  210)    
(  211)    
(  212)    SUBROUTINE zdf_ddm_init
(  213)       !!----------------------------------------------------------------------
(  214)       !!                  ***  ROUTINE zdf_ddm_init  ***
(  215)       !!
(  216)       !! ** Purpose :   Initialization of double diffusion mixing scheme
(  217)       !!
(  218)       !! ** Method  :   Read the nammbf namelist and check the parameter values
(  219)       !!      called by zdf_ddm at the first timestep (nit000)
(  220)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 5

(  221)       !! History :
(  222)       !!   8.5  !  02-08  (G. Madec)  Original code
(  223)       !!----------------------------------------------------------------------
(  224)       NAMELIST/namddm/ avts, hsbfr
(  225)       !!----------------------------------------------------------------------
(  226) 
(  227)       ! Read Namelist namddm : double diffusion mixing scheme
(  228)       ! --------------------
(  229)       REWIND ( numnam )
(  230)       READ   ( numnam, namddm )
(  231) 
(  232) 
(  233)       ! Parameter control and print
(  234)       ! ---------------------------
(  235)       IF(lwp) THEN
(  236)          WRITE(numout,*)
(  237)          WRITE(numout,*) 'zdf_ddm : double diffusive mixing'
(  238)          WRITE(numout,*) '~~~~~~~'
(  239)          WRITE(numout,*) '          Namelist namddm : set dd mixing parameter'
(  240)          WRITE(numout,*) '             maximum avs for dd mixing      avts   = ', avts
(  241)          WRITE(numout,*) '             heat/salt buoyancy flux ratio  hsbfr  = ', hsbfr
(  242)          WRITE(numout,*)
(  243)       ENDIF
(  244) 
(  245)    END SUBROUTINE zdf_ddm_init
(  246) 
(  247) #else
(  248)    !!----------------------------------------------------------------------
(  249)    !!   Default option :          Dummy module          No double diffusion
(  250)    !!----------------------------------------------------------------------
(  251)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdfddm = .FALSE.   !: double diffusion flag
(  252) CONTAINS
(  253)    SUBROUTINE zdf_ddm( kt )           ! Dummy routine
(  254) !      WRITE(*,*) 'zdf_ddm: You should not have seen this print! error?', kt
(  255)    END SUBROUTINE zdf_ddm
(  256) #endif
(  257) 
(  258)    !!======================================================================
(  259) END MODULE zdfddm

























PGF90 (Version     10.5)          09/23/2010  09:49:07      page 6

(    1) # 1 "zdfddm.F90"
(    1) MODULE zdfddm
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfddm  ***
(    4)    !! Ocean physics : double diffusion mixing parameterization
(    5)    !!======================================================================
(    6) # 248
(  248)    !!----------------------------------------------------------------------
(  249)    !!   Default option :          Dummy module          No double diffusion
(  250)    !!----------------------------------------------------------------------
(  251)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdfddm = .FALSE.   !: double diffusion flag
(  252) CONTAINS
(  253)    SUBROUTINE zdf_ddm( kt )           ! Dummy routine
(  254) !      WRITE(*,*) 'zdf_ddm: You should not have seen this print! error?', kt
(  255)    END SUBROUTINE zdf_ddm
(  256) 
(  257) # 258
(  258)    !!======================================================================
(  259) END MODULE zdfddm
