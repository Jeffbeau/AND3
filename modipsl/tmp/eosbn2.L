


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

Filename: eosbn2.F90

(    1) MODULE eosbn2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  eosbn2  ***
(    4)    !! Ocean diagnostic variable : equation of state - in situ and potential density
(    5)    !!                                               - Brunt-Vaisala frequency 
(    6)    !!==============================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   eos            : generic interface of the equation of state
(   10)    !!   eos_insitu     : Compute the in situ density
(   11)    !!   eos_insitu_pot : Compute the insitu and surface referenced potential
(   12)    !!                    volumic mass
(   13)    !!   eos_insitu_2d  : Compute the in situ density for 2d fields
(   14)    !!   eos_bn2        : Compute the Brunt-Vaisala frequency
(   15)    !!   eos_init       : set eos parameters (namelist)
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE phycst          ! physical constants
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE zdfddm          ! vertical physics: double diffusion
(   22)    USE prtctl          ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Interface 
(   28)    INTERFACE eos
(   29)       MODULE PROCEDURE eos_insitu, eos_insitu_pot, eos_insitu_2d
(   30)    END INTERFACE 
(   31)    INTERFACE bn2
(   32)       MODULE PROCEDURE eos_bn2
(   33)    END INTERFACE 
(   34) 
(   35)    !! * Routine accessibility
(   36)    PUBLIC eos        ! called by step.F90, inidtr.F90, tranpc.F90 and intgrd.F90
(   37)    PUBLIC bn2        ! called by step.F90
(   38) 
(   39)    !! * Share module variables
(   40)    INTEGER , PUBLIC ::   &  !: nameos : ocean physical parameters
(   41)       neos      = 0,     &  !: = 0/1/2 type of eq. of state and Brunt-Vaisala frequ.
(   42)       neos_init = 0         !: control flag for initialization
(   43) 
(   44)    REAL(wp), PUBLIC ::   &  !: nameos : ocean physical parameters
(   45)       ralpha = 2.0e-4,   &  !: thermal expension coeff. (linear equation of state)
(   46)       rbeta  = 7.7e-4       !: saline  expension coeff. (linear equation of state)






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 2

(   47)    
(   48)    !! * Substitutions
(   49) #  include "domzgr_substitute.h90"
(   50) #  include "vectopt_loop_substitute.h90"
(   51)    !!----------------------------------------------------------------------
(   52)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   53)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/eosbn2.F90,v 1.9 2005/09/02 15:45:19 opalod Exp $ 
(   54)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   55)    !!----------------------------------------------------------------------
(   56) 
(   57) CONTAINS
(   58) 
(   59)    SUBROUTINE eos_insitu ( ptem, psal, prd )
(   60)       !!----------------------------------------------------------------------
(   61)       !!                   ***  ROUTINE eos_insitu  ***
(   62)       !! 
(   63)       !! ** Purpose :   Compute the in situ density (ratio rho/rau0) from 
(   64)       !!       potential temperature and salinity using an equation of state
(   65)       !!       defined through the namelist parameter neos.
(   66)       !!
(   67)       !! ** Method  :   3 cases:
(   68)       !!      neos = 0 : Jackett and McDougall (1994) equation of state.
(   69)       !!         the in situ density is computed directly as a function of
(   70)       !!         potential temperature relative to the surface (the opa t
(   71)       !!         variable), salt and pressure (assuming no pressure variation
(   72)       !!         along geopotential surfaces, i.e. the pressure p in decibars
(   73)       !!         is approximated by the depth in meters.
(   74)       !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
(   75)       !!         with pressure                      p        decibars
(   76)       !!              potential temperature         t        deg celsius
(   77)       !!              salinity                      s        psu
(   78)       !!              reference volumic mass        rau0     kg/m**3
(   79)       !!              in situ volumic mass          rho      kg/m**3
(   80)       !!              in situ density anomalie      prd      no units
(   81)       !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
(   82)       !!          t = 40 deg celcius, s=40 psu
(   83)       !!      neos = 1 : linear equation of state function of temperature only
(   84)       !!              prd(t) = 0.0285 - ralpha * t
(   85)       !!      neos = 2 : linear equation of state function of temperature and
(   86)       !!               salinity
(   87)       !!              prd(t,s) = rbeta * s - ralpha * tn - 1.
(   88)       !!      Note that no boundary condition problem occurs in this routine
(   89)       !!      as (ptem,psal) are defined over the whole domain.
(   90)       !!
(   91)       !! ** Action  :   compute prd , the in situ density (no units)
(   92)       !!
(   93)       !! References :
(   94)       !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
(   95)       !!
(   96)       !! History :
(   97)       !!        !  89-03 (o. Marti)  Original code
(   98)       !!        ! 94-08 (G. Madec)
(   99)       !!        !  96-01 (G. Madec) statement function for e3
(  100)       !!        !  97-07 (G. Madec) introduction of neos, OPA8.1
(  101)       !!        !  97-07 (G. Madec) density instead of volumic mass
(  102)       !!        !  99-02 (G. Madec, N. Grima) semi-implicit pressure gradient
(  103)       !!        !  01-09 (M. Ben Jelloul) bugfix   
(  104)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 3

(  105)       !! * Arguments
(  106)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  107)          ptem,                 &  ! potential temperature
(  108)          psal                     ! salinity
(  109)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  110)          prd                      ! potential density (surface referenced)
(  111) 
(  112)       !! * Local declarations
(  113)       INTEGER ::  ji, jj, jk      ! dummy loop indices
(  114)       REAL(wp) ::   &          
(  115)          zt , zs , zh , zsr,   &  ! temporary scalars
(  116)          zr1, zr2, zr3, zr4,   &  !    "         " 
(  117)          zrhop, ze, zbw, zb,   &  !    "         "
(  118)          zd , zc , zaw, za ,   &  !    "         "
(  119)          zb1, za1, zkw, zk0       !    "         "
(  120)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  121)          zws                      ! temporary workspace
(  122)       !!----------------------------------------------------------------------
(  123) 
(  124) 
(  125)       ! initialization (in not already done)
(  126)       IF( neos_init == 0 ) CALL eos_init
(  127) 
(  128) 
(  129)       SELECT CASE ( neos )
(  130) 
(  131)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  132) 
(  133) !CDIR NOVERRCHK
(  134)          zws(:,:,:) = SQRT( ABS( psal(:,:,:) ) )
(  135) 
(  136)          !                                                ! ===============
(  137)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  138)             !                                             ! ===============
(  139)             DO jj = 1, jpj
(  140)                DO ji = 1, jpi
(  141)                   zt = ptem(ji,jj,jk)
(  142)                   zs = psal(ji,jj,jk)
(  143)                   ! depth
(  144)                   zh = fsdept(ji,jj,jk)
(  145)                   ! square root salinity
(  146)                   zsr= zws(ji,jj,jk)
(  147)                   ! compute volumic mass pure water at atm pressure
(  148)                   zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
(  149)                      -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
(  150)                   ! seawater volumic mass atm pressure
(  151)                   zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
(  152)                      -4.0899e-3 ) *zt+0.824493
(  153)                   zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
(  154)                   zr4= 4.8314e-4
(  155) 
(  156)                   ! potential volumic mass (reference to the surface)
(  157)                   zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
(  158) 
(  159)                   ! add the compression terms
(  160)                   ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
(  161)                   zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
(  162)                   zb = zbw + ze * zs






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 4

(  163) 
(  164)                   zd = -2.042967e-2
(  165)                   zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
(  166)                   zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
(  167)                   za = ( zd*zsr + zc ) *zs + zaw
(  168) 
(  169)                   zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
(  170)                   za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
(  171)                   zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
(  172)                   zk0= ( zb1*zsr + za1 )*zs + zkw
(  173) 
(  174)                   ! masked in situ density anomaly
(  175)                   prd(ji,jj,jk) = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )    &
(  176)                      - rau0 ) / rau0 * tmask(ji,jj,jk)
(  177)                END DO
(  178)             END DO
(  179)             !                                             ! ===============
(  180)          END DO                                           !   End of slab
(  181)          !                                                ! ===============
(  182) 
(  183) 
(  184)       CASE ( 1 )               ! Linear formulation function of temperature only
(  185) 
(  186)          !                                                ! ===============
(  187)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  188)             !                                             ! ===============
(  189)             DO jj = 1, jpj
(  190)                DO ji = 1, jpi
(  191)                   zt = ptem(ji,jj,jk)
(  192)                   zs = psal(ji,jj,jk)
(  193)                   !   ... density and potential volumic mass
(  194)                   prd(ji,jj,jk) = ( 0.0285 - ralpha * zt ) * tmask(ji,jj,jk)
(  195)                END DO
(  196)             END DO
(  197)             !                                             ! ===============
(  198)          END DO                                           !   End of slab
(  199)          !                                                ! ===============
(  200) 
(  201) 
(  202)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  203) 
(  204)          !                                                ! ===============
(  205)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  206)             !                                             ! ===============
(  207)             DO jj = 1, jpj
(  208)                DO ji = 1, jpi
(  209)                   zt = ptem(ji,jj,jk)
(  210)                   zs = psal(ji,jj,jk)
(  211)                   !   ... density and potential volumic mass
(  212)                   prd(ji,jj,jk) = (   rbeta  * zs - ralpha * zt ) * tmask(ji,jj,jk)
(  213)                END DO
(  214)             END DO
(  215)             !                                             ! ===============
(  216)          END DO                                           !   End of slab
(  217)          !                                                ! ===============
(  218) 
(  219)       CASE DEFAULT
(  220) 






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 5

(  221)          IF(lwp) WRITE(numout,cform_err)
(  222)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  223)          nstop = nstop + 1
(  224) 
(  225)       END SELECT
(  226) 
(  227)       IF(ln_ctl)   THEN
(  228)          CALL prt_ctl(tab3d_1=prd, clinfo1=' eos  : ', ovlap=1, kdim=jpk)
(  229)       ENDIF
(  230) 
(  231)    END SUBROUTINE eos_insitu
(  232) 
(  233) 
(  234)    SUBROUTINE eos_insitu_pot ( ptem, psal, prd, prhop)
(  235)       !!----------------------------------------------------------------------
(  236)       !!                  ***  ROUTINE eos_insitu_pot  ***
(  237)       !!           
(  238)       !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
(  239)       !!      potential volumic mass (Kg/m3) from potential temperature and
(  240)       !!      salinity fields using an equation of state defined through the 
(  241)       !!     namelist parameter neos.
(  242)       !!
(  243)       !! ** Method  :
(  244)       !!      neos = 0 : Jackett and McDougall (1994) equation of state.
(  245)       !!         the in situ density is computed directly as a function of
(  246)       !!         potential temperature relative to the surface (the opa t
(  247)       !!         variable), salt and pressure (assuming no pressure variation
(  248)       !!         along geopotential surfaces, i.e. the pressure p in decibars
(  249)       !!         is approximated by the depth in meters.
(  250)       !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
(  251)       !!              rhop(t,s)  = rho(t,s,0)
(  252)       !!         with pressure                      p        decibars
(  253)       !!              potential temperature         t        deg celsius
(  254)       !!              salinity                      s        psu
(  255)       !!              reference volumic mass        rau0     kg/m**3
(  256)       !!              in situ volumic mass          rho      kg/m**3
(  257)       !!              in situ density anomalie      prd      no units
(  258)       !!
(  259)       !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
(  260)       !!          t = 40 deg celcius, s=40 psu
(  261)       !!
(  262)       !!      neos = 1 : linear equation of state function of temperature only
(  263)       !!              prd(t) = ( rho(t) - rau0 ) / rau0 = 0.028 - ralpha * t
(  264)       !!              rhop(t,s)  = rho(t,s)
(  265)       !!
(  266)       !!      neos = 2 : linear equation of state function of temperature and
(  267)       !!               salinity
(  268)       !!              prd(t,s) = ( rho(t,s) - rau0 ) / rau0 
(  269)       !!                       = rbeta * s - ralpha * tn - 1.
(  270)       !!              rhop(t,s)  = rho(t,s)
(  271)       !!      Note that no boundary condition problem occurs in this routine
(  272)       !!      as (tn,sn) or (ta,sa) are defined over the whole domain.
(  273)       !!
(  274)       !! ** Action  : - prd  , the in situ density (no units)
(  275)       !!              - prhop, the potential volumic mass (Kg/m3)
(  276)       !!
(  277)       !! References :
(  278)       !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 6

(  279)       !!      Brown, J. A. and K. A. Campana. Mon. Weather Rev., 1978
(  280)       !!
(  281)       !! History :
(  282)       !!   4.0  !  89-03  (O. Marti)
(  283)       !!        !  94-08  (G. Madec)
(  284)       !!        !  96-01  (G. Madec) statement function for e3
(  285)       !!        !  97-07  (G. Madec) introduction of neos, OPA8.1
(  286)       !!        !  97-07  (G. Madec) density instead of volumic mass
(  287)       !!        !  99-02  (G. Madec, N. Grima) semi-implicit pressure gradient
(  288)       !!        !  01-09  (M. Ben Jelloul) bugfix   
(  289)       !!   9.0  !  03-08  (G. Madec)  F90, free form
(  290)       !!----------------------------------------------------------------------
(  291)       !! * Arguments
(  292)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  293)          ptem,   &  ! potential temperature
(  294)          psal       ! salinity
(  295)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  296)          prd,    &  ! potential density (surface referenced)
(  297)          prhop      ! potential density (surface referenced)
(  298) 
(  299)       !! * Local declarations
(  300)       INTEGER  ::  ji, jj, jk                ! dummy loop indices
(  301)       REAL(wp) ::   &             ! temporary scalars
(  302)          zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw,   &
(  303)          zb, zd, zc, zaw, za, zb1, za1, zkw, zk0
(  304)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zws
(  305)       !!----------------------------------------------------------------------
(  306) 
(  307)       ! initialization (in not already done)
(  308)       IF( neos_init == 0 ) CALL eos_init
(  309) 
(  310) 
(  311)       SELECT CASE ( neos )
(  312) 
(  313)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  314) 
(  315) !CDIR NOVERRCHK
(  316)          zws(:,:,:) = SQRT( ABS( psal(:,:,:) ) )
(  317) 
(  318)          !                                                ! ===============
(  319)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  320)             !                                             ! ===============
(  321)             DO jj = 1, jpj
(  322)                DO ji = 1, jpi
(  323)                   zt = ptem(ji,jj,jk)
(  324)                   zs = psal(ji,jj,jk)
(  325)                   ! depth
(  326)                   zh = fsdept(ji,jj,jk)
(  327)                   ! square root salinity
(  328) !!Edmee           zsr= SQRT( ABS( zs ) )
(  329)                   zsr= zws(ji,jj,jk)
(  330)                   ! compute volumic mass pure water at atm pressure
(  331)                   zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
(  332)                      -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
(  333)                   ! seawater volumic mass atm pressure
(  334)                   zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
(  335)                      -4.0899e-3 ) *zt+0.824493
(  336)                   zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 7

(  337)                   zr4= 4.8314e-4
(  338) 
(  339)                   ! potential volumic mass (reference to the surface)
(  340)                   zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
(  341) 
(  342)                   ! save potential volumic mass
(  343)                   prhop(ji,jj,jk) = zrhop * tmask(ji,jj,jk)
(  344) 
(  345)                   ! add the compression terms
(  346)                   ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
(  347)                   zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
(  348)                   zb = zbw + ze * zs
(  349) 
(  350)                   zd = -2.042967e-2
(  351)                   zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
(  352)                   zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
(  353)                   za = ( zd*zsr + zc ) *zs + zaw
(  354) 
(  355)                   zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
(  356)                   za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
(  357)                   zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
(  358)                   zk0= ( zb1*zsr + za1 )*zs + zkw
(  359) 
(  360)                   ! masked in situ density anomaly
(  361)                   prd(ji,jj,jk) = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )    &
(  362)                      - rau0 ) / rau0 * tmask(ji,jj,jk)
(  363)                END DO
(  364)             END DO
(  365)             !                                             ! ===============
(  366)          END DO                                           !   End of slab
(  367)          !                                                ! ===============
(  368) 
(  369) 
(  370)       CASE ( 1 )               ! Linear formulation function of temperature only
(  371) 
(  372)          !                                                ! ===============
(  373)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  374)             !                                             ! ===============
(  375)             DO jj = 1, jpj
(  376)                DO ji = 1, jpi
(  377)                   zt = ptem(ji,jj,jk)
(  378)                   zs = psal(ji,jj,jk)
(  379)                   !   ... density and potential volumic mass
(  380)                   prd  (ji,jj,jk) = ( 0.0285 - ralpha * zt )        * tmask(ji,jj,jk)
(  381)                   prhop(ji,jj,jk) = ( rau0 * prd(ji,jj,jk) + rau0 ) * tmask(ji,jj,jk)
(  382)                END DO
(  383)             END DO
(  384)             !                                             ! ===============
(  385)          END DO                                           !   End of slab
(  386)          !                                                ! ===============
(  387) 
(  388) 
(  389)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  390) 
(  391)          !                                                ! ===============
(  392)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  393)             !                                             ! ===============
(  394)             DO jj = 1, jpj






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 8

(  395)                DO ji = 1, jpi
(  396)                   zt = ptem(ji,jj,jk)
(  397)                   zs = psal(ji,jj,jk)
(  398)                   !   ... density and potential volumic mass
(  399)                   prd  (ji,jj,jk) = ( rbeta  * zs - ralpha * zt ) * tmask(ji,jj,jk)
(  400)                   prhop(ji,jj,jk) = ( rau0 * prd(ji,jj,jk) + rau0 )    * tmask(ji,jj,jk)
(  401)                END DO
(  402)             END DO
(  403)             !                                             ! ===============
(  404)          END DO                                           !   End of slab
(  405)          !                                                ! ===============
(  406) 
(  407)       CASE DEFAULT
(  408) 
(  409)          IF(lwp) WRITE(numout,cform_err)
(  410)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  411)          nstop = nstop + 1
(  412) 
(  413)       END SELECT
(  414) 
(  415)       IF(ln_ctl)   THEN
(  416)          CALL prt_ctl(tab3d_1=prd, clinfo1=' eos-p: ', tab3d_2=prhop, clinfo2=' pot : ', ovlap=1, kdim=jpk)
(  417)       ENDIF
(  418) 
(  419)    END SUBROUTINE eos_insitu_pot
(  420) 
(  421)    SUBROUTINE eos_insitu_2d ( ptem, psal, pdep, prd )
(  422)       !!----------------------------------------------------------------------
(  423)       !!                  ***  ROUTINE eos_insitu_2d  ***
(  424)       !!
(  425)       !! ** Purpose :   Compute the in situ density (ratio rho/rau0) from 
(  426)       !!      potential temperature and salinity using an equation of state
(  427)       !!      defined through the namelist parameter neos. * 2D field case
(  428)       !!
(  429)       !! ** Method :
(  430)       !!      neos = 0 : Jackett and McDougall (1994) equation of state.
(  431)       !!         the in situ density is computed directly as a function of
(  432)       !!         potential temperature relative to the surface (the opa t
(  433)       !!         variable), salt and pressure (assuming no pressure variation
(  434)       !!         along geopotential surfaces, i.e. the pressure p in decibars
(  435)       !!         is approximated by the depth in meters.
(  436)       !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
(  437)       !!         with pressure                      p        decibars
(  438)       !!              potential temperature         t        deg celsius
(  439)       !!              salinity                      s        psu
(  440)       !!              reference volumic mass        rau0     kg/m**3
(  441)       !!              in situ volumic mass          rho      kg/m**3
(  442)       !!              in situ density anomalie      prd      no units
(  443)       !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
(  444)       !!          t = 40 deg celcius, s=40 psu
(  445)       !!      neos = 1 : linear equation of state function of temperature only
(  446)       !!              prd(t) = 0.0285 - ralpha * t
(  447)       !!      neos = 2 : linear equation of state function of temperature and
(  448)       !!               salinity
(  449)       !!              prd(t,s) = rbeta * s - ralpha * tn - 1.
(  450)       !!      Note that no boundary condition problem occurs in this routine
(  451)       !!      as (ptem,psal) are defined over the whole domain.
(  452)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 9

(  453)       !! ** Action  : - prd , the in situ density (no units)
(  454)       !!
(  455)       !! References :
(  456)       !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
(  457)       !!
(  458)       !! History :
(  459)       !!   8.5  !  02-11  (G. Madec, A. Bozec)  partial step
(  460)       !!----------------------------------------------------------------------
(  461)       !! * Arguments
(  462)       REAL(wp), DIMENSION(jpi,jpj), INTENT( in ) ::   &
(  463)          ptem,                           &  ! potential temperature
(  464)          psal,                           &  ! salinity
(  465)          pdep                               ! depth
(  466)       REAL(wp), DIMENSION(jpi,jpj), INTENT( out ) ::   &
(  467)          prd                                ! potential density (surface referenced)
(  468) 
(  469)       !! * Local declarations
(  470)       INTEGER ::  ji, jj                    ! dummy loop indices
(  471)       REAL(wp) ::   &             ! temporary scalars
(  472)          zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw,   &
(  473)          zb, zd, zc, zaw, za, zb1, za1, zkw, zk0,               &
(  474)          zmask
(  475)       REAL(wp), DIMENSION(jpi,jpj) :: zws
(  476)       !!----------------------------------------------------------------------
(  477) 
(  478) 
(  479)       ! initialization (in not already done)
(  480)       IF( neos_init == 0 ) CALL eos_init
(  481) 
(  482)       prd(:,:) = 0.e0
(  483) 
(  484)       SELECT CASE ( neos )
(  485) 
(  486)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  487) 
(  488) !CDIR NOVERRCHK
(  489)          DO jj = 1, jpjm1
(  490) !CDIR NOVERRCHK
(  491) #if defined key_autotasking
(  492)             DO ji = 1, jpim1
(  493) #else
(  494)             DO ji = 1, fs_jpim1   ! vector opt.
(  495) #endif
(  496)                zws(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
(  497)             END DO
(  498)          END DO
(  499) 
(  500)          !                                                ! ===============
(  501)          DO jj = 1, jpjm1                                 ! Horizontal slab
(  502)             !                                             ! ===============
(  503) #if defined key_autotasking
(  504)             DO ji = 1, jpim1
(  505) #else
(  506)             DO ji = 1, fs_jpim1   ! vector opt.
(  507) #endif
(  508) 
(  509)                zmask = tmask(ji,jj,1)      ! land/sea bottom mask = surf. mask
(  510) 






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 10

(  511)                zt = ptem (ji,jj)            ! interpolated T
(  512)                zs = psal (ji,jj)            ! interpolated S
(  513)                zsr= zws(ji,jj)            ! square root of interpolated S
(  514)                zh = pdep(ji,jj)            ! depth at the partial step level
(  515) 
(  516)                ! compute volumic mass pure water at atm pressure
(  517)                zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
(  518)                          -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
(  519)                ! seawater volumic mass atm pressure
(  520)                zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 )*zt+7.6438e-5 ) *zt   &
(  521)                          -4.0899e-3 ) *zt+0.824493
(  522)                zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
(  523)                zr4= 4.8314e-4
(  524) 
(  525)                ! potential volumic mass (reference to the surface)
(  526)                zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
(  527) 
(  528)                ! add the compression terms
(  529)                ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
(  530)                zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
(  531)                zb = zbw + ze * zs
(  532) 
(  533)                zd = -2.042967e-2
(  534)                zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
(  535)                zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt -4.721788
(  536)                za = ( zd*zsr + zc ) *zs + zaw
(  537) 
(  538)                zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
(  539)                za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
(  540)                zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt   &
(  541)                          +2098.925 ) *zt+190925.6
(  542)                zk0= ( zb1*zsr + za1 )*zs + zkw
(  543) 
(  544)                ! masked in situ density anomaly
(  545)                prd(ji,jj) = ( zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rau0 )   &
(  546)                           / rau0 * zmask
(  547)             END DO
(  548)             !                                             ! ===============
(  549)          END DO                                           !   End of slab
(  550)          !                                                ! ===============
(  551) 
(  552) 
(  553)       CASE ( 1 )               ! Linear formulation function of temperature only
(  554) 
(  555)          !                                                ! ===============
(  556)          DO jj = 1, jpjm1                                 ! Horizontal slab
(  557)             !                                             ! ===============
(  558) #if defined key_autotasking
(  559)             DO ji = 1, jpim1
(  560) #else
(  561)             DO ji = 1, fs_jpim1   ! vector opt.
(  562) #endif
(  563)                prd(ji,jj) = ( 0.0285 - ralpha * ptem(ji,jj) ) * tmask(ji,jj,1)
(  564)             END DO
(  565)             !                                             ! ===============
(  566)          END DO                                           !   End of slab
(  567)          !                                                ! ===============
(  568) 






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 11

(  569) 
(  570)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  571) 
(  572)          !                                                ! ===============
(  573)          DO jj = 1, jpjm1                                 ! Horizontal slab
(  574)             !                                             ! ===============
(  575) #if defined key_autotasking
(  576)             DO ji = 1, jpim1
(  577) #else
(  578)             DO ji = 1, fs_jpim1   ! vector opt.
(  579) #endif
(  580)                prd(ji,jj) = ( rbeta * psal(ji,jj) - ralpha * ptem(ji,jj) ) * tmask(ji,jj,1) 
(  581)             END DO
(  582)             !                                             ! ===============
(  583)          END DO                                           !   End of slab
(  584)          !                                                ! ===============
(  585) 
(  586)       CASE DEFAULT
(  587) 
(  588)          IF(lwp) WRITE(numout,cform_err)
(  589)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  590)          nstop = nstop + 1
(  591) 
(  592)       END SELECT
(  593) 
(  594)       IF(ln_ctl)   CALL prt_ctl(tab2d_1=prd, clinfo1=' eos2d: ')
(  595) 
(  596)    END SUBROUTINE eos_insitu_2d
(  597) 
(  598) 
(  599)    SUBROUTINE eos_bn2( ptem, psal, pn2 )
(  600)       !!----------------------------------------------------------------------
(  601)       !!                  ***  ROUTINE eos_bn2  ***
(  602)       !!
(  603)       !! ** Purpose :   Compute the local Brunt-Vaisala frequency at the time-
(  604)       !!      step of the input arguments
(  605)       !!      
(  606)       !! ** Method :
(  607)       !!       * neos = 0  : UNESCO sea water properties
(  608)       !!         The brunt-vaisala frequency is computed using the polynomial
(  609)       !!      polynomial expression of McDougall (1987):
(  610)       !!            N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
(  611)       !!      If lk_zdfddm=T, the heat/salt buoyancy flux ratio Rrau is
(  612)       !!      computed and used in zdfddm module :
(  613)       !!              Rrau = alpha/beta * ( dk[ t ] / dk[ s ] )
(  614)       !!       * neos = 1  : linear equation of state (temperature only)
(  615)       !!            N^2 = grav * ralpha * dk[ t ]/e3w
(  616)       !!       * neos = 2  : linear equation of state (temperature & salinity)
(  617)       !!            N^2 = grav * (ralpha * dk[ t ] - rbeta * dk[ s ] ) / e3w
(  618)       !!      The use of potential density to compute N^2 introduces e r r o r
(  619)       !!      in the sign of N^2 at great depths. We recommand the use of 
(  620)       !!      neos = 0, except for academical studies.
(  621)       !!        Macro-tasked on horizontal slab (jk-loop)
(  622)       !!      N.B. N^2 is set to zero at the first level (JK=1) in inidtr
(  623)       !!      and is never used at this level.
(  624)       !!
(  625)       !! ** Action  : - pn2 : the brunt-vaisala frequency
(  626)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 12

(  627)       !! References :
(  628)       !!      McDougall, T. J., J. Phys. Oceanogr., 17, 1950-1964, 1987.
(  629)       !!
(  630)       !! History :
(  631)       !!   6.0  !  94-07  (G. Madec, M. Imbard)  Original code
(  632)       !!   8.0  !  97-07  (G. Madec) introduction of statement functions
(  633)       !!   8.5  !  02-07  (G. Madec) Free form, F90
(  634)       !!   8.5  !  02-08  (G. Madec) introduction of arguments
(  635)       !!----------------------------------------------------------------------
(  636)       !! * Arguments
(  637)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  638)          ptem,                           &  ! potential temperature
(  639)          psal                               ! salinity
(  640)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  641)          pn2                               ! Brunt-Vaisala frequency
(  642) 
(  643)       !! * Local declarations
(  644)       INTEGER  ::   ji, jj, jk   ! dummy loop indices
(  645)       REAL(wp) ::   &
(  646)          zgde3w, zt, zs, zh,  &  ! temporary scalars 
(  647)          zalbet, zbeta           !    "         "
(  648) #if defined key_zdfddm
(  649)       REAL(wp) ::   zds          ! temporary scalars
(  650) #endif
(  651)       !!----------------------------------------------------------------------
(  652)       !!  OPA8.5, LODYC-IPSL (2002)
(  653)       !!----------------------------------------------------------------------
(  654) 
(  655)       ! pn2 : first and last levels
(  656)       ! ---------------------------
(  657)       ! bn^2=0. at jk=1 and jpk set in inidtr.F : no computation
(  658) 
(  659) 
(  660)       ! pn2 : interior points only (2=< jk =< jpkm1 )
(  661)       ! -------------------------- 
(  662) 
(  663)       SELECT CASE ( neos )
(  664) 
(  665)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  666) 
(  667)          !                                                ! ===============
(  668)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  669)             !                                             ! ===============
(  670)             DO jj = 1, jpj
(  671)                DO ji = 1, jpi
(  672)                   zgde3w = grav / fse3w(ji,jj,jk)
(  673)                   zt = 0.5 * ( ptem(ji,jj,jk) + ptem(ji,jj,jk-1) )          ! potential temperature at w-point
(  674)                   zs = 0.5 * ( psal(ji,jj,jk) + psal(ji,jj,jk-1) ) - 35.0   ! salinity anomaly (s-35) at w-point
(  675)                   zh = fsdepw(ji,jj,jk)                                     ! depth in meters  at w-point
(  676) 
(  677)                   zalbet = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta
(  678)                      &                               - 0.203814e-03 ) * zt   &
(  679)                      &                               + 0.170907e-01 ) * zt   &
(  680)                      &   + 0.665157e-01                                      &
(  681)                      &   +     ( - 0.678662e-05 * zs                         &
(  682)                      &           - 0.846960e-04 * zt + 0.378110e-02 ) * zs   &
(  683)                      &   +   ( ( - 0.302285e-13 * zh                         &
(  684)                      &           - 0.251520e-11 * zs                         &






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 13

(  685)                      &           + 0.512857e-12 * zt * zt           ) * zh   &
(  686)                      &           - 0.164759e-06 * zs                         &
(  687)                      &        +(   0.791325e-08 * zt - 0.933746e-06 ) * zt   &
(  688)                      &                               + 0.380374e-04 ) * zh
(  689) 
(  690)                   zbeta  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt      &   ! beta
(  691)                      &                            - 0.301985e-05 ) * zt      &
(  692)                      &   + 0.785567e-03                                      &
(  693)                      &   + (     0.515032e-08 * zs                           &
(  694)                      &         + 0.788212e-08 * zt - 0.356603e-06 ) * zs     &
(  695)                      &   +(  (   0.121551e-17 * zh                           &
(  696)                      &         - 0.602281e-15 * zs                           &
(  697)                      &         - 0.175379e-14 * zt + 0.176621e-12 ) * zh     &
(  698)                      &                             + 0.408195e-10   * zs     &
(  699)                      &     + ( - 0.213127e-11 * zt + 0.192867e-09 ) * zt     &
(  700)                      &                             - 0.121555e-07 ) * zh
(  701) 
(  702)                   pn2(ji,jj,jk) = zgde3w * zbeta * tmask(ji,jj,jk)           &   ! N^2 
(  703)                      &          * ( zalbet * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )   &
(  704)                      &                     - ( psal(ji,jj,jk-1) - psal(ji,jj,jk) ) )
(  705) #if defined key_zdfddm
(  706)                   !                                                         !!bug **** caution a traiter zds=dk[S]= 0 !!!!
(  707)                   zds = ( psal(ji,jj,jk-1) - psal(ji,jj,jk) )                    ! Rrau = (alpha / beta) (dk[t] / dk[s])
(  708)                   IF ( ABS( zds) <= 1.e-20 ) zds = 1.e-20
(  709)                   rrau(ji,jj,jk) = zalbet * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) ) / zds
(  710) #endif
(  711)                END DO
(  712)             END DO
(  713)             !                                             ! ===============
(  714)          END DO                                           !   End of slab
(  715)          !                                                ! ===============
(  716) 
(  717) 
(  718)       CASE ( 1 )               ! Linear formulation function of temperature only
(  719) 
(  720)          !                                                ! ===============
(  721)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  722)             !                                             ! ===============
(  723)             DO jj = 1, jpj
(  724)                DO ji = 1, jpi
(  725)                   zgde3w = grav / fse3w(ji,jj,jk) * tmask(ji,jj,jk)
(  726)                   pn2(ji,jj,jk) = zgde3w * ralpha * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )
(  727)                END DO
(  728)             END DO
(  729)             !                                             ! ===============
(  730)          END DO                                           !   End of slab
(  731)          !                                                ! ===============
(  732) 
(  733) 
(  734)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  735) 
(  736)          !                                                ! ===============
(  737)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  738)             !                                             ! ===============
(  739)             DO jj = 1, jpj
(  740)                DO ji = 1, jpi
(  741)                   zgde3w = grav / fse3w(ji,jj,jk) * tmask(ji,jj,jk)
(  742)                   pn2(ji,jj,jk) = zgde3w * (  ralpha * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )   &






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 14

(  743)                      &                      - rbeta  * ( psal(ji,jj,jk-1) - psal(ji,jj,jk) )  )
(  744)                END DO
(  745)             END DO
(  746) #if defined key_zdfddm
(  747)             !                                                ! Rrau = (alpha / beta) (dk[t] / dk[s])
(  748)             zalbet = ralpha / rbeta
(  749)             DO jj = 1, jpj
(  750)                DO ji = 1, jpi
(  751)                   zds = ( psal(ji,jj,jk-1) - psal(ji,jj,jk) )
(  752)                   IF ( ABS( zds ) <= 1.e-20 ) zds = 1.e-20
(  753)                   rrau(ji,jj,jk) = zalbet * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) ) / zds
(  754)                END DO
(  755)             END DO
(  756) #endif
(  757)             !                                             ! ===============
(  758)          END DO                                           !   End of slab
(  759)          !                                                ! ===============
(  760) 
(  761)       CASE DEFAULT
(  762) 
(  763)          IF(lwp) WRITE(numout,cform_err)
(  764)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  765)          nstop = nstop + 1
(  766) 
(  767)       END SELECT
(  768) 
(  769)       IF(ln_ctl)   THEN
(  770)          CALL prt_ctl(tab3d_1=pn2, clinfo1=' bn2  : ', ovlap=1, kdim=jpk)
(  771) #if defined key_zdfddm
(  772)          CALL prt_ctl(tab3d_1=rrau, clinfo1=' rrau : ', ovlap=1, kdim=jpk)
(  773) #endif
(  774)       ENDIF
(  775) 
(  776)    END SUBROUTINE eos_bn2
(  777) 
(  778) 
(  779)    SUBROUTINE eos_init
(  780)       !!----------------------------------------------------------------------
(  781)       !!                 ***  ROUTINE eos_init  ***
(  782)       !!
(  783)       !! ** Purpose :   initializations for the equation of state
(  784)       !!
(  785)       !! ** Method  :   Read the namelist nameos
(  786)       !!
(  787)       !! ** Action  :   blahblah....
(  788)       !!
(  789)       !! History :
(  790)       !!   8.5  !  02-10  (G. Madec)  Original code
(  791)       !!----------------------------------------------------------------------
(  792)       NAMELIST/nameos/ neos, ralpha, rbeta
(  793)       !!----------------------------------------------------------------------
(  794)       !!  OPA 8.5, LODYC-IPSL (2002)
(  795)       !!----------------------------------------------------------------------
(  796) 
(  797)       ! set the initialization flag to 1
(  798)       neos_init = 1           ! indicate that the initialization has been done
(  799) 
(  800)       ! namelist nameos : ocean physical parameters






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 15

(  801) 
(  802)       ! Read Namelist nameos : equation of state
(  803)       REWIND( numnam )
(  804)       READ  ( numnam, nameos )
(  805) 
(  806)       ! Control print
(  807)       IF(lwp) THEN
(  808)          WRITE(numout,*)
(  809)          WRITE(numout,*) 'eos_init : equation of state'
(  810)          WRITE(numout,*) '~~~~~~~~'
(  811)          WRITE(numout,*) '          Namelist nameos : set eos parameters'
(  812)          WRITE(numout,*)
(  813)          WRITE(numout,*) '             flag for eq. of state and N^2  neos   = ', neos
(  814)          WRITE(numout,*) '             thermal exp. coef. (linear)    ralpha = ', ralpha
(  815)          WRITE(numout,*) '             saline  exp. coef. (linear)    rbeta  = ', rbeta
(  816)          WRITE(numout,*)
(  817)       ENDIF
(  818) 
(  819)       SELECT CASE ( neos )
(  820) 
(  821)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  822) 
(  823)          IF(lwp) WRITE(numout,*) '          use of Jackett & McDougall (1994) equation of state and'
(  824)          IF(lwp) WRITE(numout,*) '                 McDougall (1987) Brunt-Vaisala frequency'
(  825) 
(  826)       CASE ( 1 )               ! Linear formulation function of temperature only
(  827) 
(  828)          IF(lwp) WRITE(numout,*) '          use of linear eos rho(T) = rau0 * ( 1.0285 - ralpha * T )'
(  829)          IF( lk_zdfddm ) THEN
(  830)             IF(lwp) WRITE(numout,cform_err)
(  831)             IF(lwp) WRITE(numout,*) '          double diffusive mixing parameterization requires',   &
(  832)                                              ' that T and S are used as state variables'
(  833)             nstop = nstop + 1
(  834)          ENDIF
(  835) 
(  836)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  837) 
(  838)          IF(lwp) WRITE(numout,*) '          use of linear eos rho(T,S) = rau0 * ( rbeta * S - ralpha * T )'
(  839) 
(  840)       CASE DEFAULT
(  841) 
(  842)          IF(lwp) WRITE(numout,cform_err)
(  843)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  844)          nstop = nstop + 1
(  845) 
(  846)       END SELECT
(  847) 
(  848)    END SUBROUTINE eos_init
(  849) 
(  850)    !!======================================================================
(  851) END MODULE eosbn2  













PGF90 (Version     10.5)          09/23/2010  09:49:07      page 16

(    1) # 1 "eosbn2.F90"
(    1) MODULE eosbn2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  eosbn2  ***
(    4)    !! Ocean diagnostic variable : equation of state - in situ and potential density
(    5)    !!                                               - Brunt-Vaisala frequency 
(    6)    !!==============================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   eos            : generic interface of the equation of state
(   10)    !!   eos_insitu     : Compute the in situ density
(   11)    !!   eos_insitu_pot : Compute the insitu and surface referenced potential
(   12)    !!                    volumic mass
(   13)    !!   eos_insitu_2d  : Compute the in situ density for 2d fields
(   14)    !!   eos_bn2        : Compute the Brunt-Vaisala frequency
(   15)    !!   eos_init       : set eos parameters (namelist)
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE dom_oce         ! ocean space and time domain
(   19)    USE phycst          ! physical constants
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE zdfddm          ! vertical physics: double diffusion
(   22)    USE prtctl          ! Print control
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Interface 
(   28)    INTERFACE eos
(   29)       MODULE PROCEDURE eos_insitu, eos_insitu_pot, eos_insitu_2d
(   30)    END INTERFACE 
(   31)    INTERFACE bn2
(   32)       MODULE PROCEDURE eos_bn2
(   33)    END INTERFACE 
(   34) 
(   35)    !! * Routine accessibility
(   36)    PUBLIC eos        ! called by step.F90, inidtr.F90, tranpc.F90 and intgrd.F90
(   37)    PUBLIC bn2        ! called by step.F90
(   38) 
(   39)    !! * Share module variables
(   40)    INTEGER , PUBLIC ::   &  !: nameos : ocean physical parameters
(   41)       neos      = 0,     &  !: = 0/1/2 type of eq. of state and Brunt-Vaisala frequ.
(   42)       neos_init = 0         !: control flag for initialization
(   43) 
(   44)    REAL(wp), PUBLIC ::   &  !: nameos : ocean physical parameters
(   45)       ralpha = 2.0e-4,   &  !: thermal expension coeff. (linear equation of state)
(   46)       rbeta  = 7.7e-4       !: saline  expension coeff. (linear equation of state)
(   47)    
(   48)    !! * Substitutions
(   49) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 17

(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 50 "eosbn2.F90"
(   50) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 51 "eosbn2.F90"
(   51) # 51
(   51)    !!----------------------------------------------------------------------
(   52)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   53)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/eosbn2.F90,v 1.9 2005/09/02 15:45:19 opalod Exp $ 
(   54)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   55)    !!----------------------------------------------------------------------
(   56) 
(   57) CONTAINS
(   58) 
(   59)    SUBROUTINE eos_insitu ( ptem, psal, prd )
(   60)       !!----------------------------------------------------------------------
(   61)       !!                   ***  ROUTINE eos_insitu  ***
(   62)       !! 
(   63)       !! ** Purpose :   Compute the in situ density (ratio rho/rau0) from 
(   64)       !!       potential temperature and salinity using an equation of state
(   65)       !!       defined through the namelist parameter neos.
(   66)       !!
(   67)       !! ** Method  :   3 cases:
(   68)       !!      neos = 0 : Jackett and McDougall (1994) equation of state.
(   69)       !!         the in situ density is computed directly as a function of
(   70)       !!         potential temperature relative to the surface (the opa t
(   71)       !!         variable), salt and pressure (assuming no pressure variation
(   72)       !!         along geopotential surfaces, i.e. the pressure p in decibars
(   73)       !!         is approximated by the depth in meters.
(   74)       !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
(   75)       !!         with pressure                      p        decibars
(   76)       !!              potential temperature         t        deg celsius
(   77)       !!              salinity                      s        psu
(   78)       !!              reference volumic mass        rau0     kg/m**3
(   79)       !!              in situ volumic mass          rho      kg/m**3
(   80)       !!              in situ density anomalie      prd      no units
(   81)       !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
(   82)       !!          t = 40 deg celcius, s=40 psu
(   83)       !!      neos = 1 : linear equation of state function of temperature only
(   84)       !!              prd(t) = 0.0285 - ralpha * t
(   85)       !!      neos = 2 : linear equation of state function of temperature and
(   86)       !!               salinity
(   87)       !!              prd(t,s) = rbeta * s - ralpha * tn - 1.






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 18

(   88)       !!      Note that no boundary condition problem occurs in this routine
(   89)       !!      as (ptem,psal) are defined over the whole domain.
(   90)       !!
(   91)       !! ** Action  :   compute prd , the in situ density (no units)
(   92)       !!
(   93)       !! References :
(   94)       !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
(   95)       !!
(   96)       !! History :
(   97)       !!        !  89-03 (o. Marti)  Original code
(   98)       !!        ! 94-08 (G. Madec)
(   99)       !!        !  96-01 (G. Madec) statement function for e3
(  100)       !!        !  97-07 (G. Madec) introduction of neos, OPA8.1
(  101)       !!        !  97-07 (G. Madec) density instead of volumic mass
(  102)       !!        !  99-02 (G. Madec, N. Grima) semi-implicit pressure gradient
(  103)       !!        !  01-09 (M. Ben Jelloul) bugfix   
(  104)       !!----------------------------------------------------------------------
(  105)       !! * Arguments
(  106)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  107)          ptem,                 &  ! potential temperature
(  108)          psal                     ! salinity
(  109)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  110)          prd                      ! potential density (surface referenced)
(  111) 
(  112)       !! * Local declarations
(  113)       INTEGER ::  ji, jj, jk      ! dummy loop indices
(  114)       REAL(wp) ::   &          
(  115)          zt , zs , zh , zsr,   &  ! temporary scalars
(  116)          zr1, zr2, zr3, zr4,   &  !    "         " 
(  117)          zrhop, ze, zbw, zb,   &  !    "         "
(  118)          zd , zc , zaw, za ,   &  !    "         "
(  119)          zb1, za1, zkw, zk0       !    "         "
(  120)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  121)          zws                      ! temporary workspace
(  122)       !!----------------------------------------------------------------------
(  123) 
(  124) 
(  125)       ! initialization (in not already done)
(  126)       IF( neos_init == 0 ) CALL eos_init
(  127) 
(  128) 
(  129)       SELECT CASE ( neos )
(  130) 
(  131)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  132) 
(  133) !CDIR NOVERRCHK
(  134)          zws(:,:,:) = SQRT( ABS( psal(:,:,:) ) )
(  135) 
(  136)          !                                                ! ===============
(  137)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  138)             !                                             ! ===============
(  139)             DO jj = 1, jpj
(  140)                DO ji = 1, jpi
(  141)                   zt = ptem(ji,jj,jk)
(  142)                   zs = psal(ji,jj,jk)
(  143)                   ! depth
(  144)                   zh = gdept(jk)
(  145)                   ! square root salinity






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 19

(  146)                   zsr= zws(ji,jj,jk)
(  147)                   ! compute volumic mass pure water at atm pressure
(  148)                   zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
(  149)                      -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
(  150)                   ! seawater volumic mass atm pressure
(  151)                   zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
(  152)                      -4.0899e-3 ) *zt+0.824493
(  153)                   zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
(  154)                   zr4= 4.8314e-4
(  155) 
(  156)                   ! potential volumic mass (reference to the surface)
(  157)                   zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
(  158) 
(  159)                   ! add the compression terms
(  160)                   ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
(  161)                   zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
(  162)                   zb = zbw + ze * zs
(  163) 
(  164)                   zd = -2.042967e-2
(  165)                   zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
(  166)                   zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
(  167)                   za = ( zd*zsr + zc ) *zs + zaw
(  168) 
(  169)                   zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
(  170)                   za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
(  171)                   zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
(  172)                   zk0= ( zb1*zsr + za1 )*zs + zkw
(  173) 
(  174)                   ! masked in situ density anomaly
(  175)                   prd(ji,jj,jk) = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )    &
(  176)                      - rau0 ) / rau0 * tmask(ji,jj,jk)
(  177)                END DO
(  178)             END DO
(  179)             !                                             ! ===============
(  180)          END DO                                           !   End of slab
(  181)          !                                                ! ===============
(  182) 
(  183) 
(  184)       CASE ( 1 )               ! Linear formulation function of temperature only
(  185) 
(  186)          !                                                ! ===============
(  187)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  188)             !                                             ! ===============
(  189)             DO jj = 1, jpj
(  190)                DO ji = 1, jpi
(  191)                   zt = ptem(ji,jj,jk)
(  192)                   zs = psal(ji,jj,jk)
(  193)                   !   ... density and potential volumic mass
(  194)                   prd(ji,jj,jk) = ( 0.0285 - ralpha * zt ) * tmask(ji,jj,jk)
(  195)                END DO
(  196)             END DO
(  197)             !                                             ! ===============
(  198)          END DO                                           !   End of slab
(  199)          !                                                ! ===============
(  200) 
(  201) 
(  202)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  203) 






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 20

(  204)          !                                                ! ===============
(  205)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  206)             !                                             ! ===============
(  207)             DO jj = 1, jpj
(  208)                DO ji = 1, jpi
(  209)                   zt = ptem(ji,jj,jk)
(  210)                   zs = psal(ji,jj,jk)
(  211)                   !   ... density and potential volumic mass
(  212)                   prd(ji,jj,jk) = (   rbeta  * zs - ralpha * zt ) * tmask(ji,jj,jk)
(  213)                END DO
(  214)             END DO
(  215)             !                                             ! ===============
(  216)          END DO                                           !   End of slab
(  217)          !                                                ! ===============
(  218) 
(  219)       CASE DEFAULT
(  220) 
(  221)          IF(lwp) WRITE(numout,cform_err)
(  222)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  223)          nstop = nstop + 1
(  224) 
(  225)       END SELECT
(  226) 
(  227)       IF(ln_ctl)   THEN
(  228)          CALL prt_ctl(tab3d_1=prd, clinfo1=' eos  : ', ovlap=1, kdim=jpk)
(  229)       ENDIF
(  230) 
(  231)    END SUBROUTINE eos_insitu
(  232) 
(  233) 
(  234)    SUBROUTINE eos_insitu_pot ( ptem, psal, prd, prhop)
(  235)       !!----------------------------------------------------------------------
(  236)       !!                  ***  ROUTINE eos_insitu_pot  ***
(  237)       !!           
(  238)       !! ** Purpose :   Compute the in situ density (ratio rho/rau0) and the
(  239)       !!      potential volumic mass (Kg/m3) from potential temperature and
(  240)       !!      salinity fields using an equation of state defined through the 
(  241)       !!     namelist parameter neos.
(  242)       !!
(  243)       !! ** Method  :
(  244)       !!      neos = 0 : Jackett and McDougall (1994) equation of state.
(  245)       !!         the in situ density is computed directly as a function of
(  246)       !!         potential temperature relative to the surface (the opa t
(  247)       !!         variable), salt and pressure (assuming no pressure variation
(  248)       !!         along geopotential surfaces, i.e. the pressure p in decibars
(  249)       !!         is approximated by the depth in meters.
(  250)       !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
(  251)       !!              rhop(t,s)  = rho(t,s,0)
(  252)       !!         with pressure                      p        decibars
(  253)       !!              potential temperature         t        deg celsius
(  254)       !!              salinity                      s        psu
(  255)       !!              reference volumic mass        rau0     kg/m**3
(  256)       !!              in situ volumic mass          rho      kg/m**3
(  257)       !!              in situ density anomalie      prd      no units
(  258)       !!
(  259)       !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
(  260)       !!          t = 40 deg celcius, s=40 psu
(  261)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 21

(  262)       !!      neos = 1 : linear equation of state function of temperature only
(  263)       !!              prd(t) = ( rho(t) - rau0 ) / rau0 = 0.028 - ralpha * t
(  264)       !!              rhop(t,s)  = rho(t,s)
(  265)       !!
(  266)       !!      neos = 2 : linear equation of state function of temperature and
(  267)       !!               salinity
(  268)       !!              prd(t,s) = ( rho(t,s) - rau0 ) / rau0 
(  269)       !!                       = rbeta * s - ralpha * tn - 1.
(  270)       !!              rhop(t,s)  = rho(t,s)
(  271)       !!      Note that no boundary condition problem occurs in this routine
(  272)       !!      as (tn,sn) or (ta,sa) are defined over the whole domain.
(  273)       !!
(  274)       !! ** Action  : - prd  , the in situ density (no units)
(  275)       !!              - prhop, the potential volumic mass (Kg/m3)
(  276)       !!
(  277)       !! References :
(  278)       !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
(  279)       !!      Brown, J. A. and K. A. Campana. Mon. Weather Rev., 1978
(  280)       !!
(  281)       !! History :
(  282)       !!   4.0  !  89-03  (O. Marti)
(  283)       !!        !  94-08  (G. Madec)
(  284)       !!        !  96-01  (G. Madec) statement function for e3
(  285)       !!        !  97-07  (G. Madec) introduction of neos, OPA8.1
(  286)       !!        !  97-07  (G. Madec) density instead of volumic mass
(  287)       !!        !  99-02  (G. Madec, N. Grima) semi-implicit pressure gradient
(  288)       !!        !  01-09  (M. Ben Jelloul) bugfix   
(  289)       !!   9.0  !  03-08  (G. Madec)  F90, free form
(  290)       !!----------------------------------------------------------------------
(  291)       !! * Arguments
(  292)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  293)          ptem,   &  ! potential temperature
(  294)          psal       ! salinity
(  295)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  296)          prd,    &  ! potential density (surface referenced)
(  297)          prhop      ! potential density (surface referenced)
(  298) 
(  299)       !! * Local declarations
(  300)       INTEGER  ::  ji, jj, jk                ! dummy loop indices
(  301)       REAL(wp) ::   &             ! temporary scalars
(  302)          zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw,   &
(  303)          zb, zd, zc, zaw, za, zb1, za1, zkw, zk0
(  304)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zws
(  305)       !!----------------------------------------------------------------------
(  306) 
(  307)       ! initialization (in not already done)
(  308)       IF( neos_init == 0 ) CALL eos_init
(  309) 
(  310) 
(  311)       SELECT CASE ( neos )
(  312) 
(  313)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  314) 
(  315) !CDIR NOVERRCHK
(  316)          zws(:,:,:) = SQRT( ABS( psal(:,:,:) ) )
(  317) 
(  318)          !                                                ! ===============
(  319)          DO jk = 1, jpkm1                                 ! Horizontal slab






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 22

(  320)             !                                             ! ===============
(  321)             DO jj = 1, jpj
(  322)                DO ji = 1, jpi
(  323)                   zt = ptem(ji,jj,jk)
(  324)                   zs = psal(ji,jj,jk)
(  325)                   ! depth
(  326)                   zh = gdept(jk)
(  327)                   ! square root salinity
(  328) !!Edmee           zsr= SQRT( ABS( zs ) )
(  329)                   zsr= zws(ji,jj,jk)
(  330)                   ! compute volumic mass pure water at atm pressure
(  331)                   zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
(  332)                      -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
(  333)                   ! seawater volumic mass atm pressure
(  334)                   zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt   &
(  335)                      -4.0899e-3 ) *zt+0.824493
(  336)                   zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
(  337)                   zr4= 4.8314e-4
(  338) 
(  339)                   ! potential volumic mass (reference to the surface)
(  340)                   zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
(  341) 
(  342)                   ! save potential volumic mass
(  343)                   prhop(ji,jj,jk) = zrhop * tmask(ji,jj,jk)
(  344) 
(  345)                   ! add the compression terms
(  346)                   ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
(  347)                   zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
(  348)                   zb = zbw + ze * zs
(  349) 
(  350)                   zd = -2.042967e-2
(  351)                   zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
(  352)                   zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt - 4.721788
(  353)                   za = ( zd*zsr + zc ) *zs + zaw
(  354) 
(  355)                   zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
(  356)                   za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
(  357)                   zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt + 2098.925 ) *zt+190925.6
(  358)                   zk0= ( zb1*zsr + za1 )*zs + zkw
(  359) 
(  360)                   ! masked in situ density anomaly
(  361)                   prd(ji,jj,jk) = (  zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  )    &
(  362)                      - rau0 ) / rau0 * tmask(ji,jj,jk)
(  363)                END DO
(  364)             END DO
(  365)             !                                             ! ===============
(  366)          END DO                                           !   End of slab
(  367)          !                                                ! ===============
(  368) 
(  369) 
(  370)       CASE ( 1 )               ! Linear formulation function of temperature only
(  371) 
(  372)          !                                                ! ===============
(  373)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  374)             !                                             ! ===============
(  375)             DO jj = 1, jpj
(  376)                DO ji = 1, jpi
(  377)                   zt = ptem(ji,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 23

(  378)                   zs = psal(ji,jj,jk)
(  379)                   !   ... density and potential volumic mass
(  380)                   prd  (ji,jj,jk) = ( 0.0285 - ralpha * zt )        * tmask(ji,jj,jk)
(  381)                   prhop(ji,jj,jk) = ( rau0 * prd(ji,jj,jk) + rau0 ) * tmask(ji,jj,jk)
(  382)                END DO
(  383)             END DO
(  384)             !                                             ! ===============
(  385)          END DO                                           !   End of slab
(  386)          !                                                ! ===============
(  387) 
(  388) 
(  389)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  390) 
(  391)          !                                                ! ===============
(  392)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  393)             !                                             ! ===============
(  394)             DO jj = 1, jpj
(  395)                DO ji = 1, jpi
(  396)                   zt = ptem(ji,jj,jk)
(  397)                   zs = psal(ji,jj,jk)
(  398)                   !   ... density and potential volumic mass
(  399)                   prd  (ji,jj,jk) = ( rbeta  * zs - ralpha * zt ) * tmask(ji,jj,jk)
(  400)                   prhop(ji,jj,jk) = ( rau0 * prd(ji,jj,jk) + rau0 )    * tmask(ji,jj,jk)
(  401)                END DO
(  402)             END DO
(  403)             !                                             ! ===============
(  404)          END DO                                           !   End of slab
(  405)          !                                                ! ===============
(  406) 
(  407)       CASE DEFAULT
(  408) 
(  409)          IF(lwp) WRITE(numout,cform_err)
(  410)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  411)          nstop = nstop + 1
(  412) 
(  413)       END SELECT
(  414) 
(  415)       IF(ln_ctl)   THEN
(  416)          CALL prt_ctl(tab3d_1=prd, clinfo1=' eos-p: ', tab3d_2=prhop, clinfo2=' pot : ', ovlap=1, kdim=jpk)
(  417)       ENDIF
(  418) 
(  419)    END SUBROUTINE eos_insitu_pot
(  420) 
(  421)    SUBROUTINE eos_insitu_2d ( ptem, psal, pdep, prd )
(  422)       !!----------------------------------------------------------------------
(  423)       !!                  ***  ROUTINE eos_insitu_2d  ***
(  424)       !!
(  425)       !! ** Purpose :   Compute the in situ density (ratio rho/rau0) from 
(  426)       !!      potential temperature and salinity using an equation of state
(  427)       !!      defined through the namelist parameter neos. * 2D field case
(  428)       !!
(  429)       !! ** Method :
(  430)       !!      neos = 0 : Jackett and McDougall (1994) equation of state.
(  431)       !!         the in situ density is computed directly as a function of
(  432)       !!         potential temperature relative to the surface (the opa t
(  433)       !!         variable), salt and pressure (assuming no pressure variation
(  434)       !!         along geopotential surfaces, i.e. the pressure p in decibars
(  435)       !!         is approximated by the depth in meters.






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 24

(  436)       !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
(  437)       !!         with pressure                      p        decibars
(  438)       !!              potential temperature         t        deg celsius
(  439)       !!              salinity                      s        psu
(  440)       !!              reference volumic mass        rau0     kg/m**3
(  441)       !!              in situ volumic mass          rho      kg/m**3
(  442)       !!              in situ density anomalie      prd      no units
(  443)       !!         Check value: rho = 1060.93298 kg/m**3 for p=10000 dbar,
(  444)       !!          t = 40 deg celcius, s=40 psu
(  445)       !!      neos = 1 : linear equation of state function of temperature only
(  446)       !!              prd(t) = 0.0285 - ralpha * t
(  447)       !!      neos = 2 : linear equation of state function of temperature and
(  448)       !!               salinity
(  449)       !!              prd(t,s) = rbeta * s - ralpha * tn - 1.
(  450)       !!      Note that no boundary condition problem occurs in this routine
(  451)       !!      as (ptem,psal) are defined over the whole domain.
(  452)       !!
(  453)       !! ** Action  : - prd , the in situ density (no units)
(  454)       !!
(  455)       !! References :
(  456)       !!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
(  457)       !!
(  458)       !! History :
(  459)       !!   8.5  !  02-11  (G. Madec, A. Bozec)  partial step
(  460)       !!----------------------------------------------------------------------
(  461)       !! * Arguments
(  462)       REAL(wp), DIMENSION(jpi,jpj), INTENT( in ) ::   &
(  463)          ptem,                           &  ! potential temperature
(  464)          psal,                           &  ! salinity
(  465)          pdep                               ! depth
(  466)       REAL(wp), DIMENSION(jpi,jpj), INTENT( out ) ::   &
(  467)          prd                                ! potential density (surface referenced)
(  468) 
(  469)       !! * Local declarations
(  470)       INTEGER ::  ji, jj                    ! dummy loop indices
(  471)       REAL(wp) ::   &             ! temporary scalars
(  472)          zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw,   &
(  473)          zb, zd, zc, zaw, za, zb1, za1, zkw, zk0,               &
(  474)          zmask
(  475)       REAL(wp), DIMENSION(jpi,jpj) :: zws
(  476)       !!----------------------------------------------------------------------
(  477) 
(  478) 
(  479)       ! initialization (in not already done)
(  480)       IF( neos_init == 0 ) CALL eos_init
(  481) 
(  482)       prd(:,:) = 0.e0
(  483) 
(  484)       SELECT CASE ( neos )
(  485) 
(  486)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  487) 
(  488) !CDIR NOVERRCHK
(  489)          DO jj = 1, jpjm1
(  490) !CDIR NOVERRCHK
(  491) # 494
(  494)             DO ji = 1, jpim1   ! vector opt.
(  495) # 496






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 25

(  496)                zws(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
(  497)             END DO
(  498)          END DO
(  499) 
(  500)          !                                                ! ===============
(  501)          DO jj = 1, jpjm1                                 ! Horizontal slab
(  502)             !                                             ! ===============
(  503) # 506
(  506)             DO ji = 1, jpim1   ! vector opt.
(  507) 
(  508) # 509
(  509)                zmask = tmask(ji,jj,1)      ! land/sea bottom mask = surf. mask
(  510) 
(  511)                zt = ptem (ji,jj)            ! interpolated T
(  512)                zs = psal (ji,jj)            ! interpolated S
(  513)                zsr= zws(ji,jj)            ! square root of interpolated S
(  514)                zh = pdep(ji,jj)            ! depth at the partial step level
(  515) 
(  516)                ! compute volumic mass pure water at atm pressure
(  517)                zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt   &
(  518)                          -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
(  519)                ! seawater volumic mass atm pressure
(  520)                zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 )*zt+7.6438e-5 ) *zt   &
(  521)                          -4.0899e-3 ) *zt+0.824493
(  522)                zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
(  523)                zr4= 4.8314e-4
(  524) 
(  525)                ! potential volumic mass (reference to the surface)
(  526)                zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1
(  527) 
(  528)                ! add the compression terms
(  529)                ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
(  530)                zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
(  531)                zb = zbw + ze * zs
(  532) 
(  533)                zd = -2.042967e-2
(  534)                zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
(  535)                zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt -4.721788
(  536)                za = ( zd*zsr + zc ) *zs + zaw
(  537) 
(  538)                zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
(  539)                za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
(  540)                zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt   &
(  541)                          +2098.925 ) *zt+190925.6
(  542)                zk0= ( zb1*zsr + za1 )*zs + zkw
(  543) 
(  544)                ! masked in situ density anomaly
(  545)                prd(ji,jj) = ( zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) - rau0 )   &
(  546)                           / rau0 * zmask
(  547)             END DO
(  548)             !                                             ! ===============
(  549)          END DO                                           !   End of slab
(  550)          !                                                ! ===============
(  551) 
(  552) 
(  553)       CASE ( 1 )               ! Linear formulation function of temperature only
(  554) 
(  555)          !                                                ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 26

(  556)          DO jj = 1, jpjm1                                 ! Horizontal slab
(  557)             !                                             ! ===============
(  558) # 561
(  561)             DO ji = 1, jpim1   ! vector opt.
(  562) # 563
(  563)                prd(ji,jj) = ( 0.0285 - ralpha * ptem(ji,jj) ) * tmask(ji,jj,1)
(  564)             END DO
(  565)             !                                             ! ===============
(  566)          END DO                                           !   End of slab
(  567)          !                                                ! ===============
(  568) 
(  569) 
(  570)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  571) 
(  572)          !                                                ! ===============
(  573)          DO jj = 1, jpjm1                                 ! Horizontal slab
(  574)             !                                             ! ===============
(  575) # 578
(  578)             DO ji = 1, jpim1   ! vector opt.
(  579) # 580
(  580)                prd(ji,jj) = ( rbeta * psal(ji,jj) - ralpha * ptem(ji,jj) ) * tmask(ji,jj,1) 
(  581)             END DO
(  582)             !                                             ! ===============
(  583)          END DO                                           !   End of slab
(  584)          !                                                ! ===============
(  585) 
(  586)       CASE DEFAULT
(  587) 
(  588)          IF(lwp) WRITE(numout,cform_err)
(  589)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  590)          nstop = nstop + 1
(  591) 
(  592)       END SELECT
(  593) 
(  594)       IF(ln_ctl)   CALL prt_ctl(tab2d_1=prd, clinfo1=' eos2d: ')
(  595) 
(  596)    END SUBROUTINE eos_insitu_2d
(  597) 
(  598) 
(  599)    SUBROUTINE eos_bn2( ptem, psal, pn2 )
(  600)       !!----------------------------------------------------------------------
(  601)       !!                  ***  ROUTINE eos_bn2  ***
(  602)       !!
(  603)       !! ** Purpose :   Compute the local Brunt-Vaisala frequency at the time-
(  604)       !!      step of the input arguments
(  605)       !!      
(  606)       !! ** Method :
(  607)       !!       * neos = 0  : UNESCO sea water properties
(  608)       !!         The brunt-vaisala frequency is computed using the polynomial
(  609)       !!      polynomial expression of McDougall (1987):
(  610)       !!            N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
(  611)       !!      If lk_zdfddm=T, the heat/salt buoyancy flux ratio Rrau is
(  612)       !!      computed and used in zdfddm module :
(  613)       !!              Rrau = alpha/beta * ( dk[ t ] / dk[ s ] )
(  614)       !!       * neos = 1  : linear equation of state (temperature only)
(  615)       !!            N^2 = grav * ralpha * dk[ t ]/e3w
(  616)       !!       * neos = 2  : linear equation of state (temperature & salinity)
(  617)       !!            N^2 = grav * (ralpha * dk[ t ] - rbeta * dk[ s ] ) / e3w






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 27

(  618)       !!      The use of potential density to compute N^2 introduces e r r o r
(  619)       !!      in the sign of N^2 at great depths. We recommand the use of 
(  620)       !!      neos = 0, except for academical studies.
(  621)       !!        Macro-tasked on horizontal slab (jk-loop)
(  622)       !!      N.B. N^2 is set to zero at the first level (JK=1) in inidtr
(  623)       !!      and is never used at this level.
(  624)       !!
(  625)       !! ** Action  : - pn2 : the brunt-vaisala frequency
(  626)       !!
(  627)       !! References :
(  628)       !!      McDougall, T. J., J. Phys. Oceanogr., 17, 1950-1964, 1987.
(  629)       !!
(  630)       !! History :
(  631)       !!   6.0  !  94-07  (G. Madec, M. Imbard)  Original code
(  632)       !!   8.0  !  97-07  (G. Madec) introduction of statement functions
(  633)       !!   8.5  !  02-07  (G. Madec) Free form, F90
(  634)       !!   8.5  !  02-08  (G. Madec) introduction of arguments
(  635)       !!----------------------------------------------------------------------
(  636)       !! * Arguments
(  637)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ) ::   &
(  638)          ptem,                           &  ! potential temperature
(  639)          psal                               ! salinity
(  640)       REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ) ::   &
(  641)          pn2                               ! Brunt-Vaisala frequency
(  642) 
(  643)       !! * Local declarations
(  644)       INTEGER  ::   ji, jj, jk   ! dummy loop indices
(  645)       REAL(wp) ::   &
(  646)          zgde3w, zt, zs, zh,  &  ! temporary scalars 
(  647)          zalbet, zbeta           !    "         "
(  648) # 651
(  651)       !!----------------------------------------------------------------------
(  652)       !!  OPA8.5, LODYC-IPSL (2002)
(  653)       !!----------------------------------------------------------------------
(  654) 
(  655)       ! pn2 : first and last levels
(  656)       ! ---------------------------
(  657)       ! bn^2=0. at jk=1 and jpk set in inidtr.F : no computation
(  658) 
(  659) 
(  660)       ! pn2 : interior points only (2=< jk =< jpkm1 )
(  661)       ! -------------------------- 
(  662) 
(  663)       SELECT CASE ( neos )
(  664) 
(  665)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  666) 
(  667)          !                                                ! ===============
(  668)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  669)             !                                             ! ===============
(  670)             DO jj = 1, jpj
(  671)                DO ji = 1, jpi
(  672)                   zgde3w = grav / e3w(jk)
(  673)                   zt = 0.5 * ( ptem(ji,jj,jk) + ptem(ji,jj,jk-1) )          ! potential temperature at w-point
(  674)                   zs = 0.5 * ( psal(ji,jj,jk) + psal(ji,jj,jk-1) ) - 35.0   ! salinity anomaly (s-35) at w-point
(  675)                   zh = gdepw(jk)                                     ! depth in meters  at w-point
(  676) 
(  677)                   zalbet = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 28

(  678)                      &                               - 0.203814e-03 ) * zt   &
(  679)                      &                               + 0.170907e-01 ) * zt   &
(  680)                      &   + 0.665157e-01                                      &
(  681)                      &   +     ( - 0.678662e-05 * zs                         &
(  682)                      &           - 0.846960e-04 * zt + 0.378110e-02 ) * zs   &
(  683)                      &   +   ( ( - 0.302285e-13 * zh                         &
(  684)                      &           - 0.251520e-11 * zs                         &
(  685)                      &           + 0.512857e-12 * zt * zt           ) * zh   &
(  686)                      &           - 0.164759e-06 * zs                         &
(  687)                      &        +(   0.791325e-08 * zt - 0.933746e-06 ) * zt   &
(  688)                      &                               + 0.380374e-04 ) * zh
(  689) 
(  690)                   zbeta  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt      &   ! beta
(  691)                      &                            - 0.301985e-05 ) * zt      &
(  692)                      &   + 0.785567e-03                                      &
(  693)                      &   + (     0.515032e-08 * zs                           &
(  694)                      &         + 0.788212e-08 * zt - 0.356603e-06 ) * zs     &
(  695)                      &   +(  (   0.121551e-17 * zh                           &
(  696)                      &         - 0.602281e-15 * zs                           &
(  697)                      &         - 0.175379e-14 * zt + 0.176621e-12 ) * zh     &
(  698)                      &                             + 0.408195e-10   * zs     &
(  699)                      &     + ( - 0.213127e-11 * zt + 0.192867e-09 ) * zt     &
(  700)                      &                             - 0.121555e-07 ) * zh
(  701) 
(  702)                   pn2(ji,jj,jk) = zgde3w * zbeta * tmask(ji,jj,jk)           &   ! N^2 
(  703)                      &          * ( zalbet * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )   &
(  704)                      &                     - ( psal(ji,jj,jk-1) - psal(ji,jj,jk) ) )
(  705) # 711
(  711)                END DO
(  712)             END DO
(  713)             !                                             ! ===============
(  714)          END DO                                           !   End of slab
(  715)          !                                                ! ===============
(  716) 
(  717) 
(  718)       CASE ( 1 )               ! Linear formulation function of temperature only
(  719) 
(  720)          !                                                ! ===============
(  721)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  722)             !                                             ! ===============
(  723)             DO jj = 1, jpj
(  724)                DO ji = 1, jpi
(  725)                   zgde3w = grav / e3w(jk) * tmask(ji,jj,jk)
(  726)                   pn2(ji,jj,jk) = zgde3w * ralpha * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )
(  727)                END DO
(  728)             END DO
(  729)             !                                             ! ===============
(  730)          END DO                                           !   End of slab
(  731)          !                                                ! ===============
(  732) 
(  733) 
(  734)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  735) 
(  736)          !                                                ! ===============
(  737)          DO jk = 2, jpkm1                                 ! Horizontal slab
(  738)             !                                             ! ===============
(  739)             DO jj = 1, jpj
(  740)                DO ji = 1, jpi






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 29

(  741)                   zgde3w = grav / e3w(jk) * tmask(ji,jj,jk)
(  742)                   pn2(ji,jj,jk) = zgde3w * (  ralpha * ( ptem(ji,jj,jk-1) - ptem(ji,jj,jk) )   &
(  743)                      &                      - rbeta  * ( psal(ji,jj,jk-1) - psal(ji,jj,jk) )  )
(  744)                END DO
(  745)             END DO
(  746) # 757
(  757)             !                                             ! ===============
(  758)          END DO                                           !   End of slab
(  759)          !                                                ! ===============
(  760) 
(  761)       CASE DEFAULT
(  762) 
(  763)          IF(lwp) WRITE(numout,cform_err)
(  764)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  765)          nstop = nstop + 1
(  766) 
(  767)       END SELECT
(  768) 
(  769)       IF(ln_ctl)   THEN
(  770)          CALL prt_ctl(tab3d_1=pn2, clinfo1=' bn2  : ', ovlap=1, kdim=jpk)
(  771) # 774
(  774)       ENDIF
(  775) 
(  776)    END SUBROUTINE eos_bn2
(  777) 
(  778) 
(  779)    SUBROUTINE eos_init
(  780)       !!----------------------------------------------------------------------
(  781)       !!                 ***  ROUTINE eos_init  ***
(  782)       !!
(  783)       !! ** Purpose :   initializations for the equation of state
(  784)       !!
(  785)       !! ** Method  :   Read the namelist nameos
(  786)       !!
(  787)       !! ** Action  :   blahblah....
(  788)       !!
(  789)       !! History :
(  790)       !!   8.5  !  02-10  (G. Madec)  Original code
(  791)       !!----------------------------------------------------------------------
(  792)       NAMELIST/nameos/ neos, ralpha, rbeta
(  793)       !!----------------------------------------------------------------------
(  794)       !!  OPA 8.5, LODYC-IPSL (2002)
(  795)       !!----------------------------------------------------------------------
(  796) 
(  797)       ! set the initialization flag to 1
(  798)       neos_init = 1           ! indicate that the initialization has been done
(  799) 
(  800)       ! namelist nameos : ocean physical parameters
(  801) 
(  802)       ! Read Namelist nameos : equation of state
(  803)       REWIND( numnam )
(  804)       READ  ( numnam, nameos )
(  805) 
(  806)       ! Control print
(  807)       IF(lwp) THEN
(  808)          WRITE(numout,*)
(  809)          WRITE(numout,*) 'eos_init : equation of state'
(  810)          WRITE(numout,*) '~~~~~~~~'






PGF90 (Version     10.5)          09/23/2010  09:49:07      page 30

(  811)          WRITE(numout,*) '          Namelist nameos : set eos parameters'
(  812)          WRITE(numout,*)
(  813)          WRITE(numout,*) '             flag for eq. of state and N^2  neos   = ', neos
(  814)          WRITE(numout,*) '             thermal exp. coef. (linear)    ralpha = ', ralpha
(  815)          WRITE(numout,*) '             saline  exp. coef. (linear)    rbeta  = ', rbeta
(  816)          WRITE(numout,*)
(  817)       ENDIF
(  818) 
(  819)       SELECT CASE ( neos )
(  820) 
(  821)       CASE ( 0 )               ! Jackett and McDougall (1994) formulation
(  822) 
(  823)          IF(lwp) WRITE(numout,*) '          use of Jackett & McDougall (1994) equation of state and'
(  824)          IF(lwp) WRITE(numout,*) '                 McDougall (1987) Brunt-Vaisala frequency'
(  825) 
(  826)       CASE ( 1 )               ! Linear formulation function of temperature only
(  827) 
(  828)          IF(lwp) WRITE(numout,*) '          use of linear eos rho(T) = rau0 * ( 1.0285 - ralpha * T )'
(  829)          IF( lk_zdfddm ) THEN
(  830)             IF(lwp) WRITE(numout,cform_err)
(  831)             IF(lwp) WRITE(numout,*) '          double diffusive mixing parameterization requires',   &
(  832)                                              ' that T and S are used as state variables'
(  833)             nstop = nstop + 1
(  834)          ENDIF
(  835) 
(  836)       CASE ( 2 )               ! Linear formulation function of temperature and salinity
(  837) 
(  838)          IF(lwp) WRITE(numout,*) '          use of linear eos rho(T,S) = rau0 * ( rbeta * S - ralpha * T )'
(  839) 
(  840)       CASE DEFAULT
(  841) 
(  842)          IF(lwp) WRITE(numout,cform_err)
(  843)          IF(lwp) WRITE(numout,*) '          bad flag value for neos = ', neos
(  844)          nstop = nstop + 1
(  845) 
(  846)       END SELECT
(  847) 
(  848)    END SUBROUTINE eos_init
(  849) 
(  850)    !!======================================================================
(  851) END MODULE eosbn2  
