


PGF90 (Version     10.5)          09/23/2010  09:50:20      page 1

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

Filename: zdfric.F90

(    1) MODULE zdfric
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfric  ***
(    4)    !! Ocean physics:  vertical mixing coefficient compute from the local
(    5)    !!                 Richardson number dependent formulation
(    6)    !!======================================================================
(    7) #if defined key_zdfric   ||   defined key_esopa
(    8)    !!----------------------------------------------------------------------
(    9)    !!   'key_zdfric'                                             Kz = f(Ri)
(   10)    !!----------------------------------------------------------------------
(   11)    !!   zdf_ric      : update momentum and tracer Kz from the Richardson
(   12)    !!                  number computation
(   13)    !!   zdf_ric_init : initialization, namelist read, & parameters control
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers variables
(   17)    USE dom_oce         ! ocean space and time domain variables
(   18)    USE zdf_oce         ! ocean vertical physics
(   19) !  USE phycst          ! physical constants
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE lbclnk          ! ocean lateral boundary condition (or mpp link)
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC zdf_ric   ! called by step.F90
(   28) 
(   29)    !! * Shared module variables
(   30)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdfric = .TRUE.    !: Richardson vertical mixing flag
(   31) 
(   32)    !! * Module variables
(   33)    INTEGER ::               & !!! namric   richardson number dependent Kz
(   34)       nric  = 2                ! coefficient of the parameterization
(   35)    REAL(wp) ::              & !!! namric   richardson number dependent Kz
(   36)       avmri = 100.e-4_wp ,  &  ! maximum value of the vertical eddy viscosity
(   37)       alp   =   5._wp          ! coefficient of the parameterization
(   38)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   39)       tmric                    ! coef. for the horizontal mean at t-point
(   40) 
(   41)    !! * Substitutions
(   42) #  include "domzgr_substitute.h90"
(   43)    !!----------------------------------------------------------------------
(   44)    !!----------------------------------------------------------------------
(   45)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   46)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfric.F90,v 1.3 2005/03/27 18:35:26 opalod Exp $ 






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 2

(   47)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   48)    !!----------------------------------------------------------------------
(   49) CONTAINS
(   50) 
(   51)    SUBROUTINE zdf_ric( kt )
(   52)       !!----------------------------------------------------------------------
(   53)       !!                 ***  ROUTINE zdfric  ***
(   54)       !!                    
(   55)       !! ** Purpose :   Compute the before eddy viscosity and diffusivity as
(   56)       !!      a function of the local richardson number.
(   57)       !!
(   58)       !! ** Method  :   Local richardson number dependent formulation of the 
(   59)       !!      vertical eddy viscosity and diffusivity coefficients. the eddy
(   60)       !!      coefficients are given by:
(   61)       !!              avm = avm0 + avmb
(   62)       !!              avt = avm0 / (1 + alp*ri)
(   63)       !!      with    ri  = N^2 / dz(u)**2
(   64)       !!                  = e3w**2 * rn2/[ mi( dk(ub) )+mj( dk(vb) ) ]
(   65)       !!              avm0= avmri / (1 + alp*ri)**nric
(   66)       !!      Where ri is the before local Richardson number, avmri the maximum
(   67)       !!      value reaches by the vertical eddy coefficients, avmb and avtb
(   68)       !!      the background (or minimum) values of these coefficients for
(   69)       !!      momemtum and tracers, and alp, nric are adjustable parameters.
(   70)       !!      typical values used are : avm0=1.e-2 m2/s, avmb=1.e-6 m2/s
(   71)       !!      avtb=1.e-7 m2/s, alp=5. and nric=2.
(   72)       !!      this formulation needs ri>=0 : ri is set to zero if dz(rau)<0.
(   73)       !!      a numerical threshold is impose on the vertical shear (1.e-20)
(   74)       !!        N.B. the mask are required for implicit scheme, and surface
(   75)       !!      and bottom value already set in inimix.F
(   76)       !!
(   77)       !! References :
(   78)       !!      pacanowski & philander 1981, j. phys. oceanogr., 1441-1451.
(   79)       !! History :
(   80)       !!        !  87-09  (P. Andrich)  Original code
(   81)       !!        !  91-11  (G. Madec)
(   82)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   83)       !!        !  96-01  (G. Madec)  complet rewriting of multitasking
(   84)       !!                                  suppression of common work arrays
(   85)       !!        !  97-06 (G. Madec)  complete rewriting of zdfmix
(   86)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   87)       !!----------------------------------------------------------------------
(   88)       !! * Arguments
(   89)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step indexocean time step
(   90) 
(   91)       !! * Local declarations
(   92)       INTEGER ::   ji, jj, jk               ! dummy loop indices
(   93)       REAL(wp) ::   &
(   94)          zcoef, zdku, zdkv, zri, z05alp     ! temporary scalars
(   95)       REAL(wp), DIMENSION(jpi,jpj) ::   zwx ! temporary workspace
(   96) 
(   97)       IF( kt == nit000  ) CALL zdf_ric_init            ! Initialization (first time-step only)
(   98) 
(   99)       !                                                ! ===============
(  100)       DO jk = 2, jpkm1                                 ! Horizontal slab
(  101)          !                                             ! ===============
(  102)          ! Richardson number (put in zwx(ji,jj))
(  103)          ! -----------------
(  104)          ! minimum value set to zero






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 3

(  105)          DO jj = 2, jpjm1
(  106)             DO ji = 2, jpim1
(  107)                zcoef = 0.5 / fse3w(ji,jj,jk)
(  108)                ! shear of horizontal velocity
(  109)                zdku = zcoef * (  ub(ji-1,jj,jk-1) + ub(ji,jj,jk-1)   &
(  110)                                 -ub(ji-1,jj,jk  ) - ub(ji,jj,jk  )  )
(  111)                zdkv = zcoef * (  vb(ji,jj-1,jk-1) + vb(ji,jj,jk-1)   &
(  112)                                 -vb(ji,jj-1,jk  ) - vb(ji,jj,jk  )  )
(  113)                ! richardson number (minimum value set to zero)
(  114)                zri = rn2(ji,jj,jk) / ( zdku*zdku + zdkv*zdkv + 1.e-20 )
(  115)                zwx(ji,jj) = MAX( zri, 0.e0 )
(  116)             END DO
(  117)          END DO
(  118) 
(  119)          ! Boundary condition on zwx   (sign unchanged)
(  120)          CALL lbc_lnk( zwx, 'W', 1. )
(  121) 
(  122) 
(  123)          ! Vertical eddy viscosity and diffusivity coefficients
(  124)          ! -------------------------------------------------------
(  125)          ! Eddy viscosity coefficients
(  126)          z05alp = 0.5 * alp
(  127)          DO jj = 1, jpjm1
(  128)             DO ji = 1, jpim1
(  129)                avmu(ji,jj,jk) = umask(ji,jj,jk)   &
(  130)                               * avmri / ( 1. + z05alp*( zwx(ji+1,jj)+zwx(ji,jj) ) )**nric
(  131)                avmv(ji,jj,jk) = vmask(ji,jj,jk)   &
(  132)                               * avmri / ( 1. + z05alp*( zwx(ji,jj+1)+zwx(ji,jj) ) )**nric
(  133)             END DO
(  134)          END DO
(  135) 
(  136)          ! Eddy diffusivity coefficients
(  137)          DO jj = 2, jpjm1
(  138)             DO ji = 2, jpim1
(  139)                avt(ji,jj,jk) = tmric(ji,jj,jk) / ( 1. + alp * zwx(ji,jj) )   &
(  140)                              * (  avmu(ji,jj,jk) + avmu(ji-1, jj ,jk)        &
(  141)                                 + avmv(ji,jj,jk) + avmv( ji ,jj-1,jk)  )     &
(  142)                              + avtb(jk) * tmask(ji,jj,jk)
(  143)             END DO
(  144)          END DO
(  145) 
(  146)          ! Add the background coefficient on eddy viscosity
(  147)          DO jj = 2, jpjm1
(  148)             DO ji = 2, jpim1
(  149)                avmu(ji,jj,jk) = avmu(ji,jj,jk) + avmb(jk) * umask(ji,jj,jk)
(  150)                avmv(ji,jj,jk) = avmv(ji,jj,jk) + avmb(jk) * vmask(ji,jj,jk)
(  151)             END DO
(  152)          END DO
(  153)          !                                             ! ===============
(  154)       END DO                                           !   End of slab
(  155)       !                                                ! ===============
(  156) 
(  157)       ! Boundary conditions on (avt,avmu,avmv)   (unchanged sign)
(  158)       ! -----------------------===============
(  159)       CALL lbc_lnk( avt , 'W', 1. )
(  160)       CALL lbc_lnk( avmu, 'U', 1. )
(  161)       CALL lbc_lnk( avmv, 'V', 1. )
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 4

(  163)    END SUBROUTINE zdf_ric
(  164) 
(  165) 
(  166)    SUBROUTINE zdf_ric_init
(  167)       !!----------------------------------------------------------------------
(  168)       !!                 ***  ROUTINE zdfbfr_init  ***
(  169)       !!                    
(  170)       !! ** Purpose :   Initialization of the vertical eddy diffusivity and
(  171)       !!      viscosity coef. for the Richardson number dependent formulation.
(  172)       !!
(  173)       !! ** Method  :   Read the namric namelist and check the parameter values
(  174)       !!
(  175)       !! ** input   :   Namelist namric
(  176)       !!
(  177)       !! ** Action  :   increase by 1 the nstop flag is setting problem encounter
(  178)       !!
(  179)       !! history :
(  180)       !!  8.5  !  02-06  (G. Madec)  original code
(  181)       !!----------------------------------------------------------------------
(  182)       !! * local declarations
(  183)       INTEGER :: ji, jj, jk        ! dummy loop indices
(  184) 
(  185)       NAMELIST/namric/ avmri, alp, nric
(  186)       !!----------------------------------------------------------------------
(  187)       !!  OPA 8.5, LODYC-IPSL (2002)
(  188)       !!----------------------------------------------------------------------
(  189) 
(  190)       ! Read Namelist namric : richardson number dependent Kz
(  191)       ! --------------------
(  192)       REWIND ( numnam )
(  193)       READ   ( numnam, namric )
(  194) 
(  195) 
(  196)       ! Parameter control and print
(  197)       ! ---------------------------
(  198)       ! Control print
(  199)       IF(lwp) WRITE(numout,*)
(  200)       IF(lwp) WRITE(numout,*) 'zdf_ric : Ri depend vertical mixing scheme'
(  201)       IF(lwp) WRITE(numout,*) '======='
(  202)       IF(lwp) WRITE(numout,*) '          Namelist namric : set Kz(Ri) parameters'
(  203) 
(  204)       IF(lwp) THEN
(  205)          WRITE(numout,*)
(  206)          WRITE(numout,*) '             maximum vertical viscosity     avmri  = ', avmri
(  207)          WRITE(numout,*) '             coefficient                    alp    = ', alp
(  208)          WRITE(numout,*) '             coefficient                    nric   = ', nric
(  209)          WRITE(numout,*)
(  210)       ENDIF
(  211) 
(  212) 
(  213)       ! Work arrays for Ri number formulation
(  214)       ! -------------------------------------
(  215) 
(  216)       ! background eddy viscosity and diffusivity profiles
(  217)       avmb(:) = avm0
(  218)       avtb(:) = avt0
(  219) 
(  220)       ! background profile of avm (fit the theoretical/observational






PGF90 (Version     10.5)          09/23/2010  09:50:20      page 5

(  221)       !     profile shown by Krauss (1990) and avt
(  222) !!!   avtb(:) = 1.e-5 + 2.8e-8 * gdepw(:) ! m2/s
(  223) 
(  224)       ! Increase the background in the surface layers
(  225)       avmb(1) = 10.  * avmb(1)      ;      avtb(1) = 10.  * avtb(1)
(  226)       avmb(2) = 10.  * avmb(2)      ;      avtb(2) = 10.  * avtb(2)
(  227)       avmb(3) =  5.  * avmb(3)      ;      avtb(3) =  5.  * avtb(3)
(  228)       avmb(4) =  2.5 * avmb(4)      ;      avtb(4) =  2.5 * avtb(4)
(  229) 
(  230)       ! weighting mean array tmric for 4 T-points which accounts for coastal boundary conditions.
(  231)       DO jk = 1, jpk
(  232)          DO jj = 2, jpj
(  233)             DO ji = 2, jpi
(  234)                tmric(ji,jj,jk) =  tmask(ji,jj,jk)                                  &
(  235)                                / MAX( 1.,  umask(ji-1,jj  ,jk) + umask(ji,jj,jk)   &
(  236)                                          + vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)  )
(  237)             END DO
(  238)          END DO
(  239)       END DO
(  240) 
(  241)       tmric(:,1,:) = 0.e0
(  242) 
(  243)       ! Initialization of vertical eddy coef. to the background value
(  244)       DO jk = 1, jpk
(  245)          avt (:,:,jk) = avtb(jk) * tmask(:,:,jk)
(  246)          avmu(:,:,jk) = avmb(jk) * umask(:,:,jk)
(  247)          avmv(:,:,jk) = avmb(jk) * vmask(:,:,jk)
(  248)       END DO
(  249) 
(  250)    END SUBROUTINE zdf_ric_init
(  251) 
(  252) #else
(  253)    !!----------------------------------------------------------------------
(  254)    !!   Dummy module :              NO Richardson dependent vertical mixing
(  255)    !!----------------------------------------------------------------------
(  256)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdfric = .FALSE.   !: Richardson mixing flag
(  257) CONTAINS
(  258)    SUBROUTINE zdf_ric( kt )        ! Dummy routine
(  259) !      WRITE(*,*) 'zdf_ric: You should not have seen this print! error?', kt
(  260)    END SUBROUTINE zdf_ric
(  261) #endif
(  262) 
(  263)    !!======================================================================
(  264) END MODULE zdfric




















PGF90 (Version     10.5)          09/23/2010  09:50:20      page 6

(    1) # 1 "zdfric.F90"
(    1) MODULE zdfric
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfric  ***
(    4)    !! Ocean physics:  vertical mixing coefficient compute from the local
(    5)    !!                 Richardson number dependent formulation
(    6)    !!======================================================================
(    7) # 253
(  253)    !!----------------------------------------------------------------------
(  254)    !!   Dummy module :              NO Richardson dependent vertical mixing
(  255)    !!----------------------------------------------------------------------
(  256)    LOGICAL, PUBLIC, PARAMETER ::   lk_zdfric = .FALSE.   !: Richardson mixing flag
(  257) CONTAINS
(  258)    SUBROUTINE zdf_ric( kt )        ! Dummy routine
(  259) !      WRITE(*,*) 'zdf_ric: You should not have seen this print! error?', kt
(  260)    END SUBROUTINE zdf_ric
(  261) 
(  262) # 263
(  263)    !!======================================================================
(  264) END MODULE zdfric
