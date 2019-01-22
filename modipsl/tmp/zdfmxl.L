


PGF90 (Version     10.5)          09/23/2010  09:49:04      page 1

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

Filename: zdfmxl.F90

(    1) MODULE zdfmxl
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfmxl  ***
(    4)    !! Ocean physics: mixed layer depth 
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   zdf_mxl      : Compute the turbocline and mixed layer depths.
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables
(   13)    USE zdf_oce         ! ocean vertical physics
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE prtctl          ! Print control
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! * Routine accessibility
(   21)    PUBLIC zdf_mxl           ! called by step.F90
(   22) 
(   23)    !! * Shared module variables
(   24)    INTEGER, PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   25)       nmln                  !: number of level in the mixed layer
(   26)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   27)       hmld ,             &  !: mixing layer depth (turbocline) (m)
(   28)       hmlp ,             &  !: mixed layer depth  (rho=rho0+zdcrit) (m)
(   29)       hmlpt                 !: mixed layer depth at t-points (m)
(   30) 
(   31)    !! * module variables
(   32)    REAL(wp) ::   &
(   33)       avt_c = 5.e-4_wp,  &  ! Kz criterion for the turbocline depth
(   34)       rho_c = 0.01_wp       ! density criterion for mixed layer depth
(   35) 
(   36)    !! * Substitutions
(   37) #  include "domzgr_substitute.h90"
(   38)    !!----------------------------------------------------------------------
(   39)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   40)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfmxl.F90,v 1.7 2006/03/09 16:47:16 opalod Exp $ 
(   41)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   42)    !!----------------------------------------------------------------------
(   43) 
(   44) CONTAINS
(   45) 
(   46) # if defined key_autotasking






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 2

(   47)    !!----------------------------------------------------------------------
(   48)    !!   'key_autotasking'                               j-k-i loop (j-slab)
(   49)    !!----------------------------------------------------------------------
(   50) 
(   51)    SUBROUTINE zdf_mxl( kt )
(   52)       !!----------------------------------------------------------------------
(   53)       !!                    ***  ROUTINE zdfmxl  ***
(   54)       !!                   
(   55)       !! ** Purpose :   Compute the turbocline depth and the mixed layer depth
(   56)       !!      with a density criteria.
(   57)       !!
(   58)       !! ** Method  :   The turbocline depth is the depth at which the vertical 
(   59)       !!      eddy diffusivity coefficient (resulting from the vertical physics
(   60)       !!      alone, not the isopycnal part, see trazdf.F) fall below a given
(   61)       !!      value defined locally (avt_c here taken equal to 5 cm/s2)
(   62)       !!
(   63)       !! ** Action  :
(   64)       !!
(   65)       !! History :
(   66)       !!   9.0  !  03-08  (G. Madec)  autotasking optimization
(   67)       !!----------------------------------------------------------------------
(   68)       !! * Arguments
(   69)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(   70) 
(   71)       !! * Local declarations
(   72)       INTEGER ::   ji, jj, jk     ! dummy loop indices
(   73)       INTEGER ::   ik             ! temporary integer
(   74)       INTEGER, DIMENSION(jpi,jpj) ::   &
(   75)          imld                     ! temporary workspace
(   76)       !!----------------------------------------------------------------------
(   77) 
(   78)       IF( kt == nit000 ) THEN
(   79)          IF(lwp) WRITE(numout,*)
(   80)          IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
(   81)          IF(lwp) WRITE(numout,*) '~~~~~~~   auto-tasking case : j-k-i loop'
(   82)          IF(lwp) WRITE(numout,*)
(   83)       ENDIF
(   84) 
(   85)       !                                                ! ===============
(   86)       DO jj = 1, jpj                                   !  Vertical slab
(   87)          !                                             ! ===============
(   88) 
(   89)          ! 1. Turbocline depth
(   90)          ! -------------------
(   91)          ! last w-level at which avt<avt_c (starting from the bottom jk=jpk)
(   92)          ! (since avt(.,.,jpk)=0, we have jpk=< imld =< 2 )
(   93)          DO jk = jpk, 2, -1
(   94)             DO ji = 1, jpi
(   95)                IF( avt(ji,jj,jk) < avt_c ) imld(ji,jj) = jk 
(   96)             END DO
(   97)          END DO
(   98) 
(   99)          ! Turbocline depth and sub-turbocline temperature
(  100)          DO ji = 1, jpi
(  101)             ik = imld(ji,jj)
(  102)             hmld (ji,jj) = fsdepw(ji,jj,ik) * tmask(ji,jj,1)
(  103)          END DO
(  104) 






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 3

(  105) !!gm idea
(  106) !!   
(  107) !!gm     DO jk = jpk, 2, -1
(  108) !!gm        DO ji = 1, jpi
(  109) !!gm           IF( avt(ji,jj,jk) < avt_c ) hmld(ji,jj) = fsdepw(ji,jj,jk) * tmask(ji,jj,1)
(  110) !!gm        END DO
(  111) !!gm     END DO
(  112) !!gm
(  113) 
(  114)          ! 2. Mixed layer depth
(  115)          ! --------------------
(  116)          ! Initialization to the number of w ocean point mbathy
(  117)          nmln(:,jj) = mbathy(:,jj)
(  118) 
(  119)          ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
(  120)          ! (rhop defined at t-point, thus jk-1 for w-level just above)
(  121)          DO jk = jpkm1, 2, -1
(  122)             DO ji = 1, jpi
(  123)                IF( rhop(ji,jj,jk) > rhop(ji,jj,1) + rho_c )   nmln(ji,jj) = jk
(  124)             END DO
(  125)          END DO
(  126) 
(  127)          ! Mixed layer depth
(  128)          DO ji = 1, jpi
(  129)             ik = nmln(ji,jj)
(  130)             hmlp (ji,jj) = fsdepw(ji,jj,ik) * tmask(ji,jj,1)
(  131)             hmlpt(ji,jj) = fsdept(ji,jj,ik-1)
(  132)          END DO
(  133)          !                                             ! ===============
(  134)       END DO                                           !   End of slab
(  135)       !                                                ! ===============
(  136) 
(  137)       IF(ln_ctl)   THEN
(  138)          CALL prt_ctl(tab2d_1=REAL(nmln,wp), clinfo1=' nmln : ', tab2d_2=hmld, clinfo2=' hmld : ', ovlap=1)
(  139)       ENDIF
(  140)       
(  141) 
(  142)    END SUBROUTINE zdf_mxl
(  143) 
(  144) # else
(  145)    !!----------------------------------------------------------------------
(  146)    !!   Default option :                                         k-j-i loop
(  147)    !!----------------------------------------------------------------------
(  148) 
(  149)    SUBROUTINE zdf_mxl( kt )
(  150)       !!----------------------------------------------------------------------
(  151)       !!                  ***  ROUTINE zdfmxl  ***
(  152)       !!                   
(  153)       !! ** Purpose :   Compute the turbocline depth and the mixed layer depth
(  154)       !!      with density criteria.
(  155)       !!
(  156)       !! ** Method  :   The turbocline depth is the depth at which the vertical
(  157)       !!      eddy diffusivity coefficient (resulting from the vertical physics
(  158)       !!      alone, not the isopycnal part, see trazdf.F) fall below a given
(  159)       !!      value defined locally (avt_c here taken equal to 5 cm/s2)
(  160)       !!
(  161)       !! ** Action  :
(  162)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 4

(  163)       !! History :
(  164)       !!        !  94-11  (M. Imbard)  Original code
(  165)       !!   8.0  !  96-01  (E. Guilyardi)  sub mixed layer temp.
(  166)       !!   8.1  !  97-07  (G. Madec)  optimization
(  167)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  168)       !!----------------------------------------------------------------------
(  169)       !! * Arguments
(  170)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(  171) 
(  172)       !! * Local declarations
(  173)       INTEGER ::   ji, jj, jk     ! dummy loop indices
(  174)       INTEGER ::   ik             ! temporary integer
(  175)       INTEGER, DIMENSION(jpi,jpj) ::   &
(  176)          imld                     ! temporary workspace
(  177)       !!----------------------------------------------------------------------
(  178) 
(  179)       IF( kt == nit000 ) THEN
(  180)          IF(lwp) WRITE(numout,*)
(  181)          IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
(  182)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(  183)       ENDIF
(  184) 
(  185) 
(  186)       ! 1. Turbocline depth
(  187)       ! -------------------
(  188)       ! last w-level at which avt<avt_c (starting from the bottom jk=jpk)
(  189)       ! (since avt(.,.,jpk)=0, we have jpk=< imld =< 2 )
(  190)       DO jk = jpk, 2, -1
(  191)          DO jj = 1, jpj
(  192)             DO ji = 1, jpi
(  193)                IF( avt(ji,jj,jk) < avt_c ) imld(ji,jj) = jk 
(  194)             END DO
(  195)          END DO
(  196)       END DO
(  197) 
(  198)       ! Turbocline depth and sub-turbocline temperature
(  199)       DO jj = 1, jpj
(  200)          DO ji = 1, jpi
(  201)             ik = imld(ji,jj)
(  202)             hmld (ji,jj) = fsdepw(ji,jj,ik) * tmask(ji,jj,1)
(  203)          END DO
(  204)       END DO
(  205) 
(  206) !!gm idea
(  207) !!   
(  208) !!gm  DO jk = jpk, 2, -1
(  209) !!gm     DO jj = 1, jpj
(  210) !!gm        DO ji = 1, jpi
(  211) !!gm           IF( avt(ji,jj,jk) < avt_c ) hmld(ji,jj) = fsdepw(ji,jj,jk) * tmask(ji,jj,1)
(  212) !!gm        END DO
(  213) !!gm     END DO
(  214) !!gm  END DO
(  215) !!gm
(  216) 
(  217)       ! 2. Mixed layer depth
(  218)       ! --------------------
(  219)       ! Initialization to the number of w ocean point mbathy
(  220)       nmln(:,:) = mbathy(:,:)






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 5

(  221) 
(  222)       ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
(  223)       ! (rhop defined at t-point, thus jk-1 for w-level just above)
(  224)       DO jk = jpkm1, 2, -1
(  225)          DO jj = 1, jpj
(  226)             DO ji = 1, jpi
(  227)                IF( rhop(ji,jj,jk) > rhop(ji,jj,1) + rho_c )   nmln(ji,jj) = jk
(  228)             END DO
(  229)          END DO
(  230)       END DO
(  231) 
(  232)       ! Mixed layer depth
(  233)       DO jj = 1, jpj
(  234)          DO ji = 1, jpi
(  235)             ik = nmln(ji,jj)
(  236)             hmlp (ji,jj) = fsdepw(ji,jj,ik) * tmask(ji,jj,1)
(  237)             hmlpt(ji,jj) = fsdept(ji,jj,ik-1)
(  238)          END DO
(  239)       END DO
(  240) 
(  241)       IF(ln_ctl)   THEN
(  242)          CALL prt_ctl(tab2d_1=REAL(nmln,wp), clinfo1=' nmln : ', tab2d_2=hmld, clinfo2=' hmld : ', ovlap=1)
(  243)       ENDIF
(  244) 
(  245)    END SUBROUTINE zdf_mxl
(  246) #endif
(  247) 
(  248)    !!======================================================================
(  249) END MODULE zdfmxl



































PGF90 (Version     10.5)          09/23/2010  09:49:04      page 6

(    1) # 1 "zdfmxl.F90"
(    1) MODULE zdfmxl
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  zdfmxl  ***
(    4)    !! Ocean physics: mixed layer depth 
(    5)    !!======================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   zdf_mxl      : Compute the turbocline and mixed layer depths.
(    9)    !!----------------------------------------------------------------------
(   10)    !! * Modules used
(   11)    USE oce             ! ocean dynamics and tracers variables
(   12)    USE dom_oce         ! ocean space and time domain variables
(   13)    USE zdf_oce         ! ocean vertical physics
(   14)    USE in_out_manager  ! I/O manager
(   15)    USE prtctl          ! Print control
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! * Routine accessibility
(   21)    PUBLIC zdf_mxl           ! called by step.F90
(   22) 
(   23)    !! * Shared module variables
(   24)    INTEGER, PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   25)       nmln                  !: number of level in the mixed layer
(   26)    REAL(wp), PUBLIC, DIMENSION(jpi,jpj) ::   &   !:
(   27)       hmld ,             &  !: mixing layer depth (turbocline) (m)
(   28)       hmlp ,             &  !: mixed layer depth  (rho=rho0+zdcrit) (m)
(   29)       hmlpt                 !: mixed layer depth at t-points (m)
(   30) 
(   31)    !! * module variables
(   32)    REAL(wp) ::   &
(   33)       avt_c = 5.e-4_wp,  &  ! Kz criterion for the turbocline depth
(   34)       rho_c = 0.01_wp       ! density criterion for mixed layer depth
(   35) 
(   36)    !! * Substitutions
(   37) # 1 "./domzgr_substitute.h90"
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
(   48) # 38 "zdfmxl.F90"
(   38) # 38
(   38)    !!----------------------------------------------------------------------
(   39)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   40)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfmxl.F90,v 1.7 2006/03/09 16:47:16 opalod Exp $ 
(   41)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 7

(   42)    !!----------------------------------------------------------------------
(   43) 
(   44) CONTAINS
(   45) 
(   46) # 145
(  145)    !!----------------------------------------------------------------------
(  146)    !!   Default option :                                         k-j-i loop
(  147)    !!----------------------------------------------------------------------
(  148) 
(  149)    SUBROUTINE zdf_mxl( kt )
(  150)       !!----------------------------------------------------------------------
(  151)       !!                  ***  ROUTINE zdfmxl  ***
(  152)       !!                   
(  153)       !! ** Purpose :   Compute the turbocline depth and the mixed layer depth
(  154)       !!      with density criteria.
(  155)       !!
(  156)       !! ** Method  :   The turbocline depth is the depth at which the vertical
(  157)       !!      eddy diffusivity coefficient (resulting from the vertical physics
(  158)       !!      alone, not the isopycnal part, see trazdf.F) fall below a given
(  159)       !!      value defined locally (avt_c here taken equal to 5 cm/s2)
(  160)       !!
(  161)       !! ** Action  :
(  162)       !!
(  163)       !! History :
(  164)       !!        !  94-11  (M. Imbard)  Original code
(  165)       !!   8.0  !  96-01  (E. Guilyardi)  sub mixed layer temp.
(  166)       !!   8.1  !  97-07  (G. Madec)  optimization
(  167)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(  168)       !!----------------------------------------------------------------------
(  169)       !! * Arguments
(  170)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
(  171) 
(  172)       !! * Local declarations
(  173)       INTEGER ::   ji, jj, jk     ! dummy loop indices
(  174)       INTEGER ::   ik             ! temporary integer
(  175)       INTEGER, DIMENSION(jpi,jpj) ::   &
(  176)          imld                     ! temporary workspace
(  177)       !!----------------------------------------------------------------------
(  178) 
(  179)       IF( kt == nit000 ) THEN
(  180)          IF(lwp) WRITE(numout,*)
(  181)          IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
(  182)          IF(lwp) WRITE(numout,*) '~~~~~~~ '
(  183)       ENDIF
(  184) 
(  185) 
(  186)       ! 1. Turbocline depth
(  187)       ! -------------------
(  188)       ! last w-level at which avt<avt_c (starting from the bottom jk=jpk)
(  189)       ! (since avt(.,.,jpk)=0, we have jpk=< imld =< 2 )
(  190)       DO jk = jpk, 2, -1
(  191)          DO jj = 1, jpj
(  192)             DO ji = 1, jpi
(  193)                IF( avt(ji,jj,jk) < avt_c ) imld(ji,jj) = jk 
(  194)             END DO
(  195)          END DO
(  196)       END DO
(  197) 






PGF90 (Version     10.5)          09/23/2010  09:49:04      page 8

(  198)       ! Turbocline depth and sub-turbocline temperature
(  199)       DO jj = 1, jpj
(  200)          DO ji = 1, jpi
(  201)             ik = imld(ji,jj)
(  202)             hmld (ji,jj) = gdepw(ik) * tmask(ji,jj,1)
(  203)          END DO
(  204)       END DO
(  205) 
(  206) !!gm idea
(  207) !!   
(  208) !!gm  DO jk = jpk, 2, -1
(  209) !!gm     DO jj = 1, jpj
(  210) !!gm        DO ji = 1, jpi
(  211) !!gm           IF( avt(ji,jj,jk) < avt_c ) hmld(ji,jj) = gdepw(jk) * tmask(ji,jj,1)
(  212) !!gm        END DO
(  213) !!gm     END DO
(  214) !!gm  END DO
(  215) !!gm
(  216) 
(  217)       ! 2. Mixed layer depth
(  218)       ! --------------------
(  219)       ! Initialization to the number of w ocean point mbathy
(  220)       nmln(:,:) = mbathy(:,:)
(  221) 
(  222)       ! Last w-level at which rhop>=rho surf+rho_c (starting from jpk-1)
(  223)       ! (rhop defined at t-point, thus jk-1 for w-level just above)
(  224)       DO jk = jpkm1, 2, -1
(  225)          DO jj = 1, jpj
(  226)             DO ji = 1, jpi
(  227)                IF( rhop(ji,jj,jk) > rhop(ji,jj,1) + rho_c )   nmln(ji,jj) = jk
(  228)             END DO
(  229)          END DO
(  230)       END DO
(  231) 
(  232)       ! Mixed layer depth
(  233)       DO jj = 1, jpj
(  234)          DO ji = 1, jpi
(  235)             ik = nmln(ji,jj)
(  236)             hmlp (ji,jj) = gdepw(ik) * tmask(ji,jj,1)
(  237)             hmlpt(ji,jj) = gdept(ik-1)
(  238)          END DO
(  239)       END DO
(  240) 
(  241)       IF(ln_ctl)   THEN
(  242)          CALL prt_ctl(tab2d_1=REAL(nmln,wp), clinfo1=' nmln : ', tab2d_2=hmld, clinfo2=' hmld : ', ovlap=1)
(  243)       ENDIF
(  244) 
(  245)    END SUBROUTINE zdf_mxl
(  246) 
(  247) # 248
(  248)    !!======================================================================
(  249) END MODULE zdfmxl
