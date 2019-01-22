


PGF90 (Version     10.5)          09/23/2010  09:49:45      page 1

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

Filename: dynzdf_imp.F90

(    1) MODULE dynzdf_imp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  dynzdf_imp  ***
(    4)    !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_zdf_imp  : update the momentum trend with the vertical diffu-
(    9)    !!                  sion using a implicit time-stepping.
(   10)    !!----------------------------------------------------------------------
(   11)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   12)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_imp.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   13)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE phycst          ! physical constants
(   19)    USE zdf_oce         ! ocean vertical physics
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE taumod          ! surface ocean stress
(   22)    USE trdmod          ! ocean dynamics trends 
(   23)    USE trdmod_oce      ! ocean variables trends
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC dyn_zdf_imp    ! called by step.F90
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34) #  include "vectopt_loop_substitute.h90"
(   35)    !!----------------------------------------------------------------------
(   36)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_imp.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43) 
(   44)    SUBROUTINE dyn_zdf_imp( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE dyn_zdf_imp  ***






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 2

(   47)       !!                   
(   48)       !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
(   49)       !!      and the surface forcing, and add it to the general trend of 
(   50)       !!      the momentum equations.
(   51)       !!
(   52)       !! ** Method  :   The vertical momentum mixing trend is given by :
(   53)       !!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
(   54)       !!      backward time stepping
(   55)       !!      Surface boundary conditions: wind stress input
(   56)       !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
(   57)       !!      Add this trend to the general trend ua :
(   58)       !!         ua = ua + dz( avmu dz(u) )
(   59)       !!
(   60)       !! ** Action : - Update (ua,va) arrays with the after vertical diffusive
(   61)       !!               mixing trend.
(   62)       !!             - Save the trends in (ztdua,ztdva) ('l_trddyn')
(   63)       !!
(   64)       !! History :
(   65)       !!        !  90-10  (B. Blanke)  Original code
(   66)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   67)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   68)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   69)       !!---------------------------------------------------------------------
(   70)       !! * Modules used
(   71)       USE oce, ONLY :  zwd   => ta,   &  ! use ta as workspace
(   72)                        zws   => sa       ! use sa as workspace
(   73) 
(   74)       !! * Arguments
(   75)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   76) 
(   77)       !! * Local declarations
(   78)       INTEGER ::   &
(   79)          ji, jj, jk,                  &  ! dummy loop indices
(   80)          ikbu, ikbum1, ikbv, ikbvm1      ! temporary integers
(   81)       REAL(wp) ::   &
(   82)          zrau0r, z2dt,                &  ! temporary scalars
(   83)          z2dtf, zcoef, zzws, zrhs        !    "         "
(   84)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   85)          ztsx, ztsy, ztbx, ztby          ! temporary workspace arrays
(   86)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   87)          zwi, ztdua, ztdva               ! temporary workspace arrays
(   88)       !!----------------------------------------------------------------------
(   89) 
(   90)       IF( kt == nit000 ) THEN
(   91)          IF(lwp) WRITE(numout,*)
(   92)          IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
(   93)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(   94)       ENDIF
(   95) 
(   96)       ! 0. Local constant initialization
(   97)       ! --------------------------------
(   98)       zrau0r = 1. / rau0      ! inverse of the reference density
(   99)       z2dt   = 2. * rdt       ! Leap-frog environnement
(  100)       ztsx(:,:) = 0.e0
(  101)       ztsy(:,:) = 0.e0 
(  102)       ztbx(:,:) = 0.e0
(  103)       ztby(:,:) = 0.e0
(  104)       ! Euler time stepping when starting from rest






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 3

(  105)       IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt
(  106) 
(  107)       ! Save previous ua and va trends
(  108)       IF( l_trddyn )   THEN
(  109)          ztdua(:,:,:) = ua(:,:,:) 
(  110)          ztdva(:,:,:) = va(:,:,:) 
(  111)       ENDIF
(  112) 
(  113)       ! 1. Vertical diffusion on u
(  114)       ! ---------------------------
(  115)       ! Matrix and second member construction
(  116)       ! bottom boundary condition: only zws must be masked as avmu can take
(  117)       ! non zero value at the ocean bottom depending on the bottom friction
(  118)       ! used (see zdfmix.F)
(  119)       DO jk = 1, jpkm1
(  120)          DO jj = 2, jpjm1 
(  121)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  122)                zcoef = - z2dt / fse3u(ji,jj,jk)
(  123)                zwi(ji,jj,jk) = zcoef * avmu(ji,jj,jk  ) / fse3uw(ji,jj,jk  )
(  124)                zzws          = zcoef * avmu(ji,jj,jk+1) / fse3uw(ji,jj,jk+1)
(  125)                zws(ji,jj,jk) = zzws * umask(ji,jj,jk+1)
(  126)                zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zzws
(  127)             END DO
(  128)          END DO
(  129)       END DO
(  130) 
(  131)       ! Surface boudary conditions
(  132)       DO jj = 2, jpjm1 
(  133)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  134)             zwi(ji,jj,1) = 0.
(  135)             zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  136)          END DO
(  137)       END DO
(  138) 
(  139)       ! Matrix inversion starting from the first level
(  140)       !-----------------------------------------------------------------------
(  141)       !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  142)       !
(  143)       !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  144)       !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  145)       !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  146)       !        (        ...               )( ...  ) ( ...  )
(  147)       !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  148)       !
(  149)       !   m is decomposed in the product of an upper and a lower triangular matrix
(  150)       !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  151)       !   The solution (the after velocity) is in ua
(  152)       !-----------------------------------------------------------------------
(  153) 
(  154)       ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)
(  155)       DO jk = 2, jpkm1
(  156)          DO jj = 2, jpjm1   
(  157)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  158)                zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
(  159)             END DO
(  160)          END DO
(  161)       END DO
(  162) 






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 4

(  163)       ! second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1
(  164)       DO jj = 2, jpjm1   
(  165)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  166) !!! change les resultats (derniers digit, pas significativement + rapide 1* de moins)
(  167) !!!         ua(ji,jj,1) = ub(ji,jj,1)  &
(  168) !!!                      + z2dt * ( ua(ji,jj,1) + taux(ji,jj) / ( fse3u(ji,jj,1)*rau0 ) )
(  169)             z2dtf = z2dt / ( fse3u(ji,jj,1)*rau0 )
(  170)             ua(ji,jj,1) = ub(ji,jj,1)  &
(  171)                          + z2dt *  ua(ji,jj,1) + z2dtf * taux(ji,jj)
(  172)          END DO
(  173)       END DO
(  174)       DO jk = 2, jpkm1
(  175)          DO jj = 2, jpjm1   
(  176)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  177)                zrhs = ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)   ! zrhs=right hand side
(  178)                ua(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * ua(ji,jj,jk-1)
(  179)             END DO
(  180)          END DO
(  181)       END DO
(  182) 
(  183)       ! thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk
(  184)       DO jj = 2, jpjm1   
(  185)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  186)             ua(ji,jj,jpkm1) = ua(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
(  187)          END DO
(  188)       END DO
(  189)       DO jk = jpk-2, 1, -1
(  190)          DO jj = 2, jpjm1   
(  191)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  192)                ua(ji,jj,jk) =( ua(ji,jj,jk) - zws(ji,jj,jk) * ua(ji,jj,jk+1) ) / zwd(ji,jj,jk)
(  193)             END DO
(  194)          END DO
(  195)       END DO
(  196) 
(  197)       IF( l_trddyn )  THEN 
(  198)          ! diagnose surface and bottom momentum fluxes
(  199)          DO jj = 2, jpjm1   
(  200)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  201)                ! save the surface forcing momentum fluxes
(  202)                ztsx(ji,jj) = taux(ji,jj) / ( fse3u(ji,jj,1)*rau0 )
(  203)                ! save bottom friction momentum fluxes
(  204)                ikbu   = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
(  205)                ikbum1 = MAX( ikbu-1, 1 )
(  206)                ztbx(ji,jj) = - avmu(ji,jj,ikbu) * ua(ji,jj,ikbum1)   &
(  207)                   / ( fse3u(ji,jj,ikbum1)*fse3uw(ji,jj,ikbu) )
(  208)                ! subtract surface forcing and bottom friction trend from vertical
(  209)                ! diffusive momentum trend
(  210)                ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
(  211)                ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
(  212)             END DO
(  213)          END DO
(  214)       ENDIF
(  215) 
(  216)       ! Normalization to obtain the general momentum trend ua
(  217)       DO jk = 1, jpkm1
(  218)          DO jj = 2, jpjm1   
(  219)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  220)                ua(ji,jj,jk) = ( ua(ji,jj,jk) - ub(ji,jj,jk) ) / z2dt






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 5

(  221)             END DO
(  222)          END DO
(  223)       END DO
(  224) 
(  225) 
(  226)       ! 2. Vertical diffusion on v
(  227)       ! ---------------------------
(  228)       ! Matrix and second member construction
(  229)       ! bottom boundary condition: only zws must be masked as avmv can take
(  230)       ! non zero value at the ocean bottom depending on the bottom friction
(  231)       ! used (see zdfmix.F)
(  232)       DO jk = 1, jpkm1
(  233)          DO jj = 2, jpjm1   
(  234)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  235)                zcoef = -z2dt / fse3v(ji,jj,jk)
(  236)                zwi(ji,jj,jk) = zcoef * avmv(ji,jj,jk  ) / fse3vw(ji,jj,jk  )
(  237)                zzws       = zcoef * avmv(ji,jj,jk+1) / fse3vw(ji,jj,jk+1)
(  238)                zws(ji,jj,jk) =  zzws * vmask(ji,jj,jk+1)
(  239)                zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zzws
(  240)             END DO
(  241)          END DO
(  242)       END DO
(  243) 
(  244)       ! Surface boudary conditions
(  245)       DO jj = 2, jpjm1   
(  246)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  247)             zwi(ji,jj,1) = 0.e0
(  248)             zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  249)          END DO
(  250)       END DO
(  251) 
(  252)       ! Matrix inversion
(  253)       !-----------------------------------------------------------------------
(  254)       !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  255)       !
(  256)       !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  257)       !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  258)       !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  259)       !        (        ...               )( ...  ) ( ...  )
(  260)       !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  261)       !
(  262)       !   m is decomposed in the product of an upper and lower triangular
(  263)       !   matrix
(  264)       !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  265)       !   The solution (after velocity) is in 2d array va
(  266)       !-----------------------------------------------------------------------
(  267) 
(  268)       ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)
(  269)       DO jk = 2, jpkm1
(  270)          DO jj = 2, jpjm1   
(  271)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  272)                zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
(  273)             END DO
(  274)          END DO
(  275)       END DO
(  276) 
(  277)       ! second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1
(  278)       DO jj = 2, jpjm1






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 6

(  279)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  280) !!! change les resultats (derniers digit, pas significativement + rapide 1* de moins)
(  281) !!!         va(ji,jj,1) = vb(ji,jj,1)  &
(  282) !!!                      + z2dt * ( va(ji,jj,1) + tauy(ji,jj) / ( fse3v(ji,jj,1)*rau0 ) )
(  283)             z2dtf = z2dt / ( fse3v(ji,jj,1)*rau0 )
(  284)             va(ji,jj,1) = vb(ji,jj,1)  &
(  285)                          + z2dt * va(ji,jj,1) + z2dtf * tauy(ji,jj)
(  286)          END DO
(  287)       END DO
(  288)       DO jk = 2, jpkm1
(  289)          DO jj = 2, jpjm1
(  290)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  291)                zrhs = vb(ji,jj,jk) + z2dt * va(ji,jj,jk)   ! zrhs=right hand side
(  292)                va(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * va(ji,jj,jk-1)
(  293)             END DO
(  294)          END DO
(  295)       END DO
(  296) 
(  297)       ! thrid recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk
(  298)       DO jj = 2, jpjm1   
(  299)          DO ji = fs_2, fs_jpim1   ! vector opt.
(  300)             va(ji,jj,jpkm1) = va(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
(  301)          END DO
(  302)       END DO
(  303)       DO jk = jpk-2, 1, -1
(  304)          DO jj = 2, jpjm1   
(  305)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  306)                va(ji,jj,jk) =( va(ji,jj,jk) - zws(ji,jj,jk) * va(ji,jj,jk+1) ) / zwd(ji,jj,jk)
(  307)             END DO
(  308)          END DO
(  309)       END DO
(  310) 
(  311)       IF( l_trddyn )  THEN 
(  312)          ! diagnose surface and bottom momentum fluxes
(  313)          DO jj = 2, jpjm1   
(  314)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  315)                ! save the surface forcing momentum fluxes
(  316)                ztsy(ji,jj) = tauy(ji,jj) / ( fse3v(ji,jj,1)*rau0 )
(  317)                ! save bottom friction momentum fluxes
(  318)                ikbv   = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
(  319)                ikbvm1 = MAX( ikbv-1, 1 )
(  320)                ztby(ji,jj) = - avmv(ji,jj,ikbv) * va(ji,jj,ikbvm1)   &
(  321)                   / ( fse3v(ji,jj,ikbvm1)*fse3vw(ji,jj,ikbv) )
(  322)                ! subtract surface forcing and bottom friction trend from vertical
(  323)                ! diffusive momentum trend
(  324)                ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
(  325)                ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
(  326)             END DO
(  327)          END DO
(  328)       ENDIF
(  329) 
(  330)       ! Normalization to obtain the general momentum trend va
(  331)       DO jk = 1, jpkm1
(  332)          DO jj = 2, jpjm1   
(  333)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  334)                va(ji,jj,jk) = ( va(ji,jj,jk) - vb(ji,jj,jk) ) / z2dt
(  335)             END DO
(  336)          END DO






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 7

(  337)       END DO
(  338) 
(  339)       ! save the vertical diffusion trends for diagnostic
(  340)       ! momentum trends
(  341)       IF( l_trddyn )  THEN 
(  342)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  343)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  344) 
(  345)          CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
(  346)          ztdua(:,:,:) = 0.e0
(  347)          ztdva(:,:,:) = 0.e0
(  348)          ztdua(:,:,1) = ztsx(:,:)
(  349)          ztdva(:,:,1) = ztsy(:,:)
(  350)          CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
(  351)          ztdua(:,:,:) = 0.e0
(  352)          ztdva(:,:,:) = 0.e0
(  353)          ztdua(:,:,1) = ztbx(:,:)
(  354)          ztdva(:,:,1) = ztby(:,:)
(  355)          CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
(  356)       ENDIF
(  357) 
(  358)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  359)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
(  360)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  361)       ENDIF
(  362) 
(  363)    END SUBROUTINE dyn_zdf_imp
(  364) 
(  365)    !!==============================================================================
(  366) END MODULE dynzdf_imp


































PGF90 (Version     10.5)          09/23/2010  09:49:45      page 8

(    1) # 1 "dynzdf_imp.F90"
(    1) MODULE dynzdf_imp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  dynzdf_imp  ***
(    4)    !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   dyn_zdf_imp  : update the momentum trend with the vertical diffu-
(    9)    !!                  sion using a implicit time-stepping.
(   10)    !!----------------------------------------------------------------------
(   11)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   12)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_imp.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   13)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   14)    !!----------------------------------------------------------------------
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE phycst          ! physical constants
(   19)    USE zdf_oce         ! ocean vertical physics
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE taumod          ! surface ocean stress
(   22)    USE trdmod          ! ocean dynamics trends 
(   23)    USE trdmod_oce      ! ocean variables trends
(   24)    USE prtctl          ! Print control
(   25) 
(   26)    IMPLICIT NONE
(   27)    PRIVATE
(   28) 
(   29)    !! * Routine accessibility
(   30)    PUBLIC dyn_zdf_imp    ! called by step.F90
(   31) 
(   32)    !! * Substitutions
(   33) # 1 "./domzgr_substitute.h90"
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
(   48) # 34 "dynzdf_imp.F90"
(   34) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 
(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 9

(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 35 "dynzdf_imp.F90"
(   35) # 35
(   35)    !!----------------------------------------------------------------------
(   36)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DYN/dynzdf_imp.F90,v 1.7 2005/09/02 15:45:24 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43) 
(   44)    SUBROUTINE dyn_zdf_imp( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE dyn_zdf_imp  ***
(   47)       !!                   
(   48)       !! ** Purpose :   Compute the trend due to the vert. momentum diffusion
(   49)       !!      and the surface forcing, and add it to the general trend of 
(   50)       !!      the momentum equations.
(   51)       !!
(   52)       !! ** Method  :   The vertical momentum mixing trend is given by :
(   53)       !!             dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ua) )
(   54)       !!      backward time stepping
(   55)       !!      Surface boundary conditions: wind stress input
(   56)       !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F)
(   57)       !!      Add this trend to the general trend ua :
(   58)       !!         ua = ua + dz( avmu dz(u) )
(   59)       !!
(   60)       !! ** Action : - Update (ua,va) arrays with the after vertical diffusive
(   61)       !!               mixing trend.
(   62)       !!             - Save the trends in (ztdua,ztdva) ('l_trddyn')
(   63)       !!
(   64)       !! History :
(   65)       !!        !  90-10  (B. Blanke)  Original code
(   66)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   67)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   68)       !!   9.0  !  04-08  (C. Talandier)  New trends organization
(   69)       !!---------------------------------------------------------------------
(   70)       !! * Modules used
(   71)       USE oce, ONLY :  zwd   => ta,   &  ! use ta as workspace
(   72)                        zws   => sa       ! use sa as workspace
(   73) 
(   74)       !! * Arguments
(   75)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   76) 
(   77)       !! * Local declarations
(   78)       INTEGER ::   &
(   79)          ji, jj, jk,                  &  ! dummy loop indices
(   80)          ikbu, ikbum1, ikbv, ikbvm1      ! temporary integers
(   81)       REAL(wp) ::   &
(   82)          zrau0r, z2dt,                &  ! temporary scalars
(   83)          z2dtf, zcoef, zzws, zrhs        !    "         "
(   84)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   85)          ztsx, ztsy, ztbx, ztby          ! temporary workspace arrays
(   86)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   87)          zwi, ztdua, ztdva               ! temporary workspace arrays






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 10

(   88)       !!----------------------------------------------------------------------
(   89) 
(   90)       IF( kt == nit000 ) THEN
(   91)          IF(lwp) WRITE(numout,*)
(   92)          IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
(   93)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(   94)       ENDIF
(   95) 
(   96)       ! 0. Local constant initialization
(   97)       ! --------------------------------
(   98)       zrau0r = 1. / rau0      ! inverse of the reference density
(   99)       z2dt   = 2. * rdt       ! Leap-frog environnement
(  100)       ztsx(:,:) = 0.e0
(  101)       ztsy(:,:) = 0.e0 
(  102)       ztbx(:,:) = 0.e0
(  103)       ztby(:,:) = 0.e0
(  104)       ! Euler time stepping when starting from rest
(  105)       IF( neuler == 0 .AND. kt == nit000 )   z2dt = rdt
(  106) 
(  107)       ! Save previous ua and va trends
(  108)       IF( l_trddyn )   THEN
(  109)          ztdua(:,:,:) = ua(:,:,:) 
(  110)          ztdva(:,:,:) = va(:,:,:) 
(  111)       ENDIF
(  112) 
(  113)       ! 1. Vertical diffusion on u
(  114)       ! ---------------------------
(  115)       ! Matrix and second member construction
(  116)       ! bottom boundary condition: only zws must be masked as avmu can take
(  117)       ! non zero value at the ocean bottom depending on the bottom friction
(  118)       ! used (see zdfmix.F)
(  119)       DO jk = 1, jpkm1
(  120)          DO jj = 2, jpjm1 
(  121)             DO ji = 2, jpim1   ! vector opt.
(  122)                zcoef = - z2dt / e3t(jk)
(  123)                zwi(ji,jj,jk) = zcoef * avmu(ji,jj,jk  ) / e3w(jk  )
(  124)                zzws          = zcoef * avmu(ji,jj,jk+1) / e3w(jk+1)
(  125)                zws(ji,jj,jk) = zzws * umask(ji,jj,jk+1)
(  126)                zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zzws
(  127)             END DO
(  128)          END DO
(  129)       END DO
(  130) 
(  131)       ! Surface boudary conditions
(  132)       DO jj = 2, jpjm1 
(  133)          DO ji = 2, jpim1   ! vector opt.
(  134)             zwi(ji,jj,1) = 0.
(  135)             zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  136)          END DO
(  137)       END DO
(  138) 
(  139)       ! Matrix inversion starting from the first level
(  140)       !-----------------------------------------------------------------------
(  141)       !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  142)       !
(  143)       !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  144)       !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  145)       !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 11

(  146)       !        (        ...               )( ...  ) ( ...  )
(  147)       !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  148)       !
(  149)       !   m is decomposed in the product of an upper and a lower triangular matrix
(  150)       !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  151)       !   The solution (the after velocity) is in ua
(  152)       !-----------------------------------------------------------------------
(  153) 
(  154)       ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)
(  155)       DO jk = 2, jpkm1
(  156)          DO jj = 2, jpjm1   
(  157)             DO ji = 2, jpim1   ! vector opt.
(  158)                zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
(  159)             END DO
(  160)          END DO
(  161)       END DO
(  162) 
(  163)       ! second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1
(  164)       DO jj = 2, jpjm1   
(  165)          DO ji = 2, jpim1   ! vector opt.
(  166) !!! change les resultats (derniers digit, pas significativement + rapide 1* de moins)
(  167) !!!         ua(ji,jj,1) = ub(ji,jj,1)  &
(  168) !!!                      + z2dt * ( ua(ji,jj,1) + taux(ji,jj) / ( e3t(1)*rau0 ) )
(  169)             z2dtf = z2dt / ( e3t(1)*rau0 )
(  170)             ua(ji,jj,1) = ub(ji,jj,1)  &
(  171)                          + z2dt *  ua(ji,jj,1) + z2dtf * taux(ji,jj)
(  172)          END DO
(  173)       END DO
(  174)       DO jk = 2, jpkm1
(  175)          DO jj = 2, jpjm1   
(  176)             DO ji = 2, jpim1   ! vector opt.
(  177)                zrhs = ub(ji,jj,jk) + z2dt * ua(ji,jj,jk)   ! zrhs=right hand side
(  178)                ua(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * ua(ji,jj,jk-1)
(  179)             END DO
(  180)          END DO
(  181)       END DO
(  182) 
(  183)       ! thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk
(  184)       DO jj = 2, jpjm1   
(  185)          DO ji = 2, jpim1   ! vector opt.
(  186)             ua(ji,jj,jpkm1) = ua(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
(  187)          END DO
(  188)       END DO
(  189)       DO jk = jpk-2, 1, -1
(  190)          DO jj = 2, jpjm1   
(  191)             DO ji = 2, jpim1   ! vector opt.
(  192)                ua(ji,jj,jk) =( ua(ji,jj,jk) - zws(ji,jj,jk) * ua(ji,jj,jk+1) ) / zwd(ji,jj,jk)
(  193)             END DO
(  194)          END DO
(  195)       END DO
(  196) 
(  197)       IF( l_trddyn )  THEN 
(  198)          ! diagnose surface and bottom momentum fluxes
(  199)          DO jj = 2, jpjm1   
(  200)             DO ji = 2, jpim1   ! vector opt.
(  201)                ! save the surface forcing momentum fluxes
(  202)                ztsx(ji,jj) = taux(ji,jj) / ( e3t(1)*rau0 )
(  203)                ! save bottom friction momentum fluxes






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 12

(  204)                ikbu   = MIN( mbathy(ji+1,jj), mbathy(ji,jj) )
(  205)                ikbum1 = MAX( ikbu-1, 1 )
(  206)                ztbx(ji,jj) = - avmu(ji,jj,ikbu) * ua(ji,jj,ikbum1)   &
(  207)                   / ( e3t(ikbum1)*e3w(ikbu) )
(  208)                ! subtract surface forcing and bottom friction trend from vertical
(  209)                ! diffusive momentum trend
(  210)                ztdua(ji,jj,1     ) = ztdua(ji,jj,1     ) - ztsx(ji,jj)
(  211)                ztdua(ji,jj,ikbum1) = ztdua(ji,jj,ikbum1) - ztbx(ji,jj)
(  212)             END DO
(  213)          END DO
(  214)       ENDIF
(  215) 
(  216)       ! Normalization to obtain the general momentum trend ua
(  217)       DO jk = 1, jpkm1
(  218)          DO jj = 2, jpjm1   
(  219)             DO ji = 2, jpim1   ! vector opt.
(  220)                ua(ji,jj,jk) = ( ua(ji,jj,jk) - ub(ji,jj,jk) ) / z2dt
(  221)             END DO
(  222)          END DO
(  223)       END DO
(  224) 
(  225) 
(  226)       ! 2. Vertical diffusion on v
(  227)       ! ---------------------------
(  228)       ! Matrix and second member construction
(  229)       ! bottom boundary condition: only zws must be masked as avmv can take
(  230)       ! non zero value at the ocean bottom depending on the bottom friction
(  231)       ! used (see zdfmix.F)
(  232)       DO jk = 1, jpkm1
(  233)          DO jj = 2, jpjm1   
(  234)             DO ji = 2, jpim1   ! vector opt.
(  235)                zcoef = -z2dt / e3t(jk)
(  236)                zwi(ji,jj,jk) = zcoef * avmv(ji,jj,jk  ) / e3w(jk  )
(  237)                zzws       = zcoef * avmv(ji,jj,jk+1) / e3w(jk+1)
(  238)                zws(ji,jj,jk) =  zzws * vmask(ji,jj,jk+1)
(  239)                zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zzws
(  240)             END DO
(  241)          END DO
(  242)       END DO
(  243) 
(  244)       ! Surface boudary conditions
(  245)       DO jj = 2, jpjm1   
(  246)          DO ji = 2, jpim1   ! vector opt.
(  247)             zwi(ji,jj,1) = 0.e0
(  248)             zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  249)          END DO
(  250)       END DO
(  251) 
(  252)       ! Matrix inversion
(  253)       !-----------------------------------------------------------------------
(  254)       !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  255)       !
(  256)       !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  257)       !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  258)       !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  259)       !        (        ...               )( ...  ) ( ...  )
(  260)       !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  261)       !






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 13

(  262)       !   m is decomposed in the product of an upper and lower triangular
(  263)       !   matrix
(  264)       !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  265)       !   The solution (after velocity) is in 2d array va
(  266)       !-----------------------------------------------------------------------
(  267) 
(  268)       ! First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)
(  269)       DO jk = 2, jpkm1
(  270)          DO jj = 2, jpjm1   
(  271)             DO ji = 2, jpim1   ! vector opt.
(  272)                zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
(  273)             END DO
(  274)          END DO
(  275)       END DO
(  276) 
(  277)       ! second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1
(  278)       DO jj = 2, jpjm1
(  279)          DO ji = 2, jpim1   ! vector opt.
(  280) !!! change les resultats (derniers digit, pas significativement + rapide 1* de moins)
(  281) !!!         va(ji,jj,1) = vb(ji,jj,1)  &
(  282) !!!                      + z2dt * ( va(ji,jj,1) + tauy(ji,jj) / ( e3t(1)*rau0 ) )
(  283)             z2dtf = z2dt / ( e3t(1)*rau0 )
(  284)             va(ji,jj,1) = vb(ji,jj,1)  &
(  285)                          + z2dt * va(ji,jj,1) + z2dtf * tauy(ji,jj)
(  286)          END DO
(  287)       END DO
(  288)       DO jk = 2, jpkm1
(  289)          DO jj = 2, jpjm1
(  290)             DO ji = 2, jpim1   ! vector opt.
(  291)                zrhs = vb(ji,jj,jk) + z2dt * va(ji,jj,jk)   ! zrhs=right hand side
(  292)                va(ji,jj,jk) = zrhs - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * va(ji,jj,jk-1)
(  293)             END DO
(  294)          END DO
(  295)       END DO
(  296) 
(  297)       ! thrid recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk
(  298)       DO jj = 2, jpjm1   
(  299)          DO ji = 2, jpim1   ! vector opt.
(  300)             va(ji,jj,jpkm1) = va(ji,jj,jpkm1) / zwd(ji,jj,jpkm1)
(  301)          END DO
(  302)       END DO
(  303)       DO jk = jpk-2, 1, -1
(  304)          DO jj = 2, jpjm1   
(  305)             DO ji = 2, jpim1   ! vector opt.
(  306)                va(ji,jj,jk) =( va(ji,jj,jk) - zws(ji,jj,jk) * va(ji,jj,jk+1) ) / zwd(ji,jj,jk)
(  307)             END DO
(  308)          END DO
(  309)       END DO
(  310) 
(  311)       IF( l_trddyn )  THEN 
(  312)          ! diagnose surface and bottom momentum fluxes
(  313)          DO jj = 2, jpjm1   
(  314)             DO ji = 2, jpim1   ! vector opt.
(  315)                ! save the surface forcing momentum fluxes
(  316)                ztsy(ji,jj) = tauy(ji,jj) / ( e3t(1)*rau0 )
(  317)                ! save bottom friction momentum fluxes
(  318)                ikbv   = MIN( mbathy(ji,jj+1), mbathy(ji,jj) )
(  319)                ikbvm1 = MAX( ikbv-1, 1 )






PGF90 (Version     10.5)          09/23/2010  09:49:45      page 14

(  320)                ztby(ji,jj) = - avmv(ji,jj,ikbv) * va(ji,jj,ikbvm1)   &
(  321)                   / ( e3t(ikbvm1)*e3w(ikbv) )
(  322)                ! subtract surface forcing and bottom friction trend from vertical
(  323)                ! diffusive momentum trend
(  324)                ztdva(ji,jj,1     ) = ztdva(ji,jj,1     ) - ztsy(ji,jj)
(  325)                ztdva(ji,jj,ikbvm1) = ztdva(ji,jj,ikbvm1) - ztby(ji,jj)
(  326)             END DO
(  327)          END DO
(  328)       ENDIF
(  329) 
(  330)       ! Normalization to obtain the general momentum trend va
(  331)       DO jk = 1, jpkm1
(  332)          DO jj = 2, jpjm1   
(  333)             DO ji = 2, jpim1   ! vector opt.
(  334)                va(ji,jj,jk) = ( va(ji,jj,jk) - vb(ji,jj,jk) ) / z2dt
(  335)             END DO
(  336)          END DO
(  337)       END DO
(  338) 
(  339)       ! save the vertical diffusion trends for diagnostic
(  340)       ! momentum trends
(  341)       IF( l_trddyn )  THEN 
(  342)          ztdua(:,:,:) = ua(:,:,:) - ztdua(:,:,:)
(  343)          ztdva(:,:,:) = va(:,:,:) - ztdva(:,:,:)
(  344) 
(  345)          CALL trd_mod(ztdua, ztdva, jpdtdzdf, 'DYN', kt)
(  346)          ztdua(:,:,:) = 0.e0
(  347)          ztdva(:,:,:) = 0.e0
(  348)          ztdua(:,:,1) = ztsx(:,:)
(  349)          ztdva(:,:,1) = ztsy(:,:)
(  350)          CALL trd_mod(ztdua , ztdva , jpdtdswf, 'DYN', kt)
(  351)          ztdua(:,:,:) = 0.e0
(  352)          ztdva(:,:,:) = 0.e0
(  353)          ztdua(:,:,1) = ztbx(:,:)
(  354)          ztdva(:,:,1) = ztby(:,:)
(  355)          CALL trd_mod(ztdua , ztdva , jpdtdbfr, 'DYN', kt)
(  356)       ENDIF
(  357) 
(  358)       IF(ln_ctl) THEN         ! print sum trends (used for debugging)
(  359)          CALL prt_ctl(tab3d_1=ua, clinfo1=' zdf  - Ua: ', mask1=umask, &
(  360)             &         tab3d_2=va, clinfo2=' Va: ', mask2=vmask, clinfo3='dyn')
(  361)       ENDIF
(  362) 
(  363)    END SUBROUTINE dyn_zdf_imp
(  364) 
(  365)    !!==============================================================================
(  366) END MODULE dynzdf_imp
