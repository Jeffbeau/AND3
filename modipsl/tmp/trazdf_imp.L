


PGF90 (Version     10.5)          09/23/2010  09:50:17      page 1

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

Filename: trazdf_imp.F90

(    1) MODULE trazdf_imp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trazdf_imp  ***
(    4)    !! Ocean active tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_zdf_imp  : update the tracer trend with the vertical diffusion
(    9)    !!                  using an implicit time-stepping.
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE zdf_oce         ! ocean vertical physics
(   15)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   16)    USE zdfddm          ! ocean vertical physics: double diffusion
(   17)    USE zdfkpp          ! KPP parameterisation
(   18)    USE trdmod          ! ocean active tracers trends 
(   19)    USE trdmod_oce      ! ocean variables trends
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE prtctl          ! Print control
(   22) 
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC tra_zdf_imp          ! routine called by step.F90
(   29) 
(   30)    !! * Module variable
(   31)    REAL(wp), DIMENSION(jpk) ::   &
(   32)       r2dt                     ! vertical profile of 2 x tracer time-step
(   33) 
(   34)    !! * Substitutions
(   35) #  include "domzgr_substitute.h90"
(   36) #  include "zdfddm_substitute.h90"
(   37)    !!----------------------------------------------------------------------
(   38)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_imp.F90,v 1.6 2005/09/02 15:45:34 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------
(   42) 
(   43) CONTAINS
(   44) 
(   45)    SUBROUTINE tra_zdf_imp( kt )
(   46)       !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 2

(   47)       !!                  ***  ROUTINE tra_zdf_imp  ***
(   48)       !!
(   49)       !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
(   50)       !!      using an implicit time stepping and add it to the general trend
(   51)       !!      of the tracer equations.
(   52)       !!
(   53)       !! ** Method  :   The vertical diffusion of tracers (t & s) is given by:
(   54)       !!          difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(ta) )
(   55)       !!      It is thus evaluated using a backward time scheme
(   56)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   57)       !!      both tracers (bottom, applied through the masked field avt).
(   58)       !!      Add this trend to the general trend ta,sa :
(   59)       !!          ta = ta + dz( avt dz(t) )
(   60)       !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T)
(   61)       !!
(   62)       !! ** Action  : - Update (ta,sa) with the before vertical diffusion trend
(   63)       !!              - save the trends in (ttrd,strd) ('key_trdtra')
(   64)       !!
(   65)       !! History :
(   66)       !!   6.0  !  90-10  (B. Blanke)  Original code
(   67)       !!   7.0  !  91-11  (G. Madec)
(   68)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(   69)       !!        !  96-01  (G. Madec)  statement function for e3
(   70)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   71)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(   72)       !!        !  00-08  (G. Madec)  double diffusive mixing
(   73)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   74)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   75)       !!   9.0  !  05-01  (C. Ethe )  non-local flux in KPP vertical mixing scheme
(   76)       !!---------------------------------------------------------------------
(   77)       !! * Modules used     
(   78)       USE oce, ONLY :    ztdta => ua,       & ! use ua as 3D workspace   
(   79)                          ztdsa => va          ! use va as 3D workspace   
(   80) 
(   81)       !! * Arguments
(   82)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   83) 
(   84)       !! * Local declarations
(   85)       INTEGER ::   ji, jj, jk                 ! dummy loop indices
(   86)       INTEGER ::   ikst, ikenm2, ikstp1
(   87)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(   88)          zwd, zws, zwi,          &  ! ???
(   89)          zwx, zwy, zwz, zwt         ! ???
(   90)       !!---------------------------------------------------------------------
(   91) 
(   92) 
(   93)       ! 0. Local constant initialization
(   94)       ! --------------------------------
(   95) 
(   96)       ! time step = 2 rdttra ex
(   97)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(   98)          r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
(   99)       ELSEIF( kt <= nit000 + 1) THEN
(  100)          r2dt(:) = 2. * rdttra(:)          ! leapfrog
(  101)       ENDIF
(  102) 
(  103)       ! Save ta and sa trends
(  104)       IF( l_trdtra )   THEN






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 3

(  105)          ztdta(:,:,:) = ta(:,:,:) 
(  106)          ztdsa(:,:,:) = sa(:,:,:) 
(  107)       ENDIF
(  108) 
(  109)       !                                                ! ===============
(  110)       DO jj = 2, jpjm1                                 !  Vertical slab
(  111)          !                                             ! ===============
(  112)          ! 0. Matrix construction
(  113)          ! ----------------------
(  114) 
(  115)          ! Diagonal, inferior, superior (including the bottom boundary condition via avt masked)
(  116)          DO jk = 1, jpkm1
(  117)             DO ji = 2, jpim1
(  118)                zwi(ji,jk) = - r2dt(jk) * avt(ji,jj,jk  )   &
(  119)                                        / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  120)                zws(ji,jk) = - r2dt(jk) * avt(ji,jj,jk+1)   &
(  121)                                        / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  122)                zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
(  123)             END DO
(  124)          END DO
(  125) 
(  126)          ! Surface boudary conditions
(  127)          DO ji = 2, jpim1
(  128)             zwi(ji,1) = 0.e0
(  129)             zwd(ji,1) = 1. - zws(ji,1)
(  130)          END DO
(  131) 
(  132)          ! 1. Vertical diffusion on temperature
(  133)          ! -------------------------===========
(  134) 
(  135)          ! Second member construction
(  136) #if defined key_zdfkpp
(  137)          ! add non-local temperature flux ( in convective case only)
(  138)          DO jk = 1, jpkm1
(  139)             DO ji = 2, jpim1  
(  140)                zwy(ji,jk) = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)  &
(  141)                   &  - r2dt(jk) * ( ghats(ji,jj,jk) * avt(ji,jj,jk) - ghats(ji,jj,jk+1) * avt(ji,jj,jk+1) ) &
(  142)                   &               * wt0(ji,jj) / fse3t(ji,jj,jk) 
(  143)             END DO
(  144)          END DO
(  145) #else
(  146)          DO jk = 1, jpkm1
(  147)             DO ji = 2, jpim1             
(  148)                zwy(ji,jk) = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)
(  149)             END DO
(  150)          END DO
(  151) #endif
(  152) 
(  153)          ! Matrix inversion from the first level
(  154)          ikst = 1
(  155) 
(  156) #   include "zdf.matrixsolver.h90"
(  157) 
(  158)          ! Save the masked temperature after in ta
(  159)          ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done
(  160)          !                  it will not be done in tranxt)
(  161)          DO jk = 1, jpkm1
(  162)             DO ji = 2, jpim1






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 4

(  163)                ta(ji,jj,jk) = zwx(ji,jk) * tmask(ji,jj,jk)
(  164)             END DO
(  165)          END DO
(  166) 
(  167) 
(  168)          ! 2. Vertical diffusion on salinity
(  169)          ! -------------------------========
(  170) 
(  171) #if defined key_zdfddm
(  172)          ! Rebuild the Matrix as avt /= avs
(  173) 
(  174)          ! Diagonal, inferior, superior
(  175)          ! (including the bottom boundary condition via avs masked)
(  176)          DO jk = 1, jpkm1
(  177)             DO ji = 2, jpim1
(  178)                zwi(ji,jk) = - r2dt(jk) * fsavs(ji,jj,jk  )   &
(  179)                   /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  180)                zws(ji,jk) = - r2dt(jk) * fsavs(ji,jj,jk+1)   &
(  181)                   /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  182)                zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
(  183)             END DO
(  184)          END DO
(  185) 
(  186)          ! Surface boudary conditions
(  187)          DO ji = 2, jpim1
(  188)             zwi(ji,1) = 0.e0
(  189)             zwd(ji,1) = 1. - zws(ji,1)
(  190)          END DO
(  191) #endif
(  192)          ! Second member construction
(  193) #if defined key_zdfkpp
(  194)          ! add non-local salinity flux ( in convective case only)
(  195)          DO jk = 1, jpkm1
(  196)             DO ji = 2, jpim1  
(  197)                zwy(ji,jk) = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)  &
(  198)                   &  - r2dt(jk) * ( ghats(ji,jj,jk) * fsavs(ji,jj,jk) - ghats(ji,jj,jk+1) * fsavs(ji,jj,jk+1) ) &
(  199)                   &               * ws0(ji,jj) / fse3t(ji,jj,jk) 
(  200)             END DO
(  201)          END DO
(  202) #else
(  203)          DO jk = 1, jpkm1
(  204)             DO ji = 2, jpim1             
(  205)                zwy(ji,jk) = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)
(  206)             END DO
(  207)          END DO
(  208) #endif
(  209)  
(  210)          ! Matrix inversion from the first level
(  211)          ikst = 1
(  212) 
(  213) #   include "zdf.matrixsolver.h90"
(  214) 
(  215) 
(  216)          ! Save the masked salinity after in sa
(  217)          ! (c a u t i o n:  salinity not its trend, Leap-frog scheme done
(  218)          !                  it will not be done in tranxt)
(  219)          DO jk = 1, jpkm1
(  220)             DO ji = 2, jpim1






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 5

(  221)                sa(ji,jj,jk) = zwx(ji,jk)  * tmask(ji,jj,jk)
(  222)             END DO
(  223)          END DO
(  224) 
(  225)          !                                             ! ===============
(  226)       END DO                                           !   End of slab
(  227)       !                                                ! ===============
(  228) 
(  229)       ! save the trends for diagnostic
(  230)       ! Compute and save the vertical diffusive temperature & salinity trends
(  231)       IF( l_trdtra )   THEN
(  232)          ! compute the vertical diffusive trends in substracting the previous 
(  233)          ! trends ztdta()/ztdsa() to the new one computed (dT/dt or dS/dt) 
(  234)          ! with the new temperature/salinity ta/sa
(  235)          DO jk = 1, jpkm1
(  236)             ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) )   & ! new trend
(  237)                 &           - ztdta(:,:,jk)                                ! old trend
(  238)             ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) )   & ! new trend
(  239)                 &           - ztdsa(:,:,jk)                                ! old trend
(  240)          END DO
(  241) 
(  242)          CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
(  243)       ENDIF
(  244) 
(  245)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  246)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf  - Ta: ', mask1=tmask, &
(  247)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  248)       ENDIF
(  249) 
(  250)    END SUBROUTINE tra_zdf_imp
(  251) 
(  252)    !!==============================================================================
(  253) END MODULE trazdf_imp































PGF90 (Version     10.5)          09/23/2010  09:50:17      page 6

(    1) # 1 "trazdf_imp.F90"
(    1) MODULE trazdf_imp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trazdf_imp  ***
(    4)    !! Ocean active tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_zdf_imp  : update the tracer trend with the vertical diffusion
(    9)    !!                  using an implicit time-stepping.
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE zdf_oce         ! ocean vertical physics
(   15)    USE ldftra_oce      ! ocean active tracers: lateral physics
(   16)    USE zdfddm          ! ocean vertical physics: double diffusion
(   17)    USE zdfkpp          ! KPP parameterisation
(   18)    USE trdmod          ! ocean active tracers trends 
(   19)    USE trdmod_oce      ! ocean variables trends
(   20)    USE in_out_manager  ! I/O manager
(   21)    USE prtctl          ! Print control
(   22) 
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC tra_zdf_imp          ! routine called by step.F90
(   29) 
(   30)    !! * Module variable
(   31)    REAL(wp), DIMENSION(jpk) ::   &
(   32)       r2dt                     ! vertical profile of 2 x tracer time-step
(   33) 
(   34)    !! * Substitutions
(   35) # 1 "./domzgr_substitute.h90"
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
(   48) # 36 "trazdf_imp.F90"
(   36) # 1 "./zdfddm_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** zdfddm_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 7

(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfddm_substitute.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 16
(   16) !   Defautl option :                     avs = avt
(   17) # 37 "trazdf_imp.F90"
(   37) # 37
(   37)    !!----------------------------------------------------------------------
(   38)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   39)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_imp.F90,v 1.6 2005/09/02 15:45:34 opalod Exp $ 
(   40)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   41)    !!----------------------------------------------------------------------
(   42) 
(   43) CONTAINS
(   44) 
(   45)    SUBROUTINE tra_zdf_imp( kt )
(   46)       !!----------------------------------------------------------------------
(   47)       !!                  ***  ROUTINE tra_zdf_imp  ***
(   48)       !!
(   49)       !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
(   50)       !!      using an implicit time stepping and add it to the general trend
(   51)       !!      of the tracer equations.
(   52)       !!
(   53)       !! ** Method  :   The vertical diffusion of tracers (t & s) is given by:
(   54)       !!          difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(ta) )
(   55)       !!      It is thus evaluated using a backward time scheme
(   56)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   57)       !!      both tracers (bottom, applied through the masked field avt).
(   58)       !!      Add this trend to the general trend ta,sa :
(   59)       !!          ta = ta + dz( avt dz(t) )
(   60)       !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T)
(   61)       !!
(   62)       !! ** Action  : - Update (ta,sa) with the before vertical diffusion trend
(   63)       !!              - save the trends in (ttrd,strd) ('key_trdtra')
(   64)       !!
(   65)       !! History :
(   66)       !!   6.0  !  90-10  (B. Blanke)  Original code
(   67)       !!   7.0  !  91-11  (G. Madec)
(   68)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(   69)       !!        !  96-01  (G. Madec)  statement function for e3
(   70)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   71)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(   72)       !!        !  00-08  (G. Madec)  double diffusive mixing
(   73)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   74)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   75)       !!   9.0  !  05-01  (C. Ethe )  non-local flux in KPP vertical mixing scheme
(   76)       !!---------------------------------------------------------------------
(   77)       !! * Modules used     
(   78)       USE oce, ONLY :    ztdta => ua,       & ! use ua as 3D workspace   
(   79)                          ztdsa => va          ! use va as 3D workspace   
(   80) 
(   81)       !! * Arguments
(   82)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   83) 
(   84)       !! * Local declarations
(   85)       INTEGER ::   ji, jj, jk                 ! dummy loop indices






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 8

(   86)       INTEGER ::   ikst, ikenm2, ikstp1
(   87)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(   88)          zwd, zws, zwi,          &  ! ???
(   89)          zwx, zwy, zwz, zwt         ! ???
(   90)       !!---------------------------------------------------------------------
(   91) 
(   92) 
(   93)       ! 0. Local constant initialization
(   94)       ! --------------------------------
(   95) 
(   96)       ! time step = 2 rdttra ex
(   97)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(   98)          r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
(   99)       ELSEIF( kt <= nit000 + 1) THEN
(  100)          r2dt(:) = 2. * rdttra(:)          ! leapfrog
(  101)       ENDIF
(  102) 
(  103)       ! Save ta and sa trends
(  104)       IF( l_trdtra )   THEN
(  105)          ztdta(:,:,:) = ta(:,:,:) 
(  106)          ztdsa(:,:,:) = sa(:,:,:) 
(  107)       ENDIF
(  108) 
(  109)       !                                                ! ===============
(  110)       DO jj = 2, jpjm1                                 !  Vertical slab
(  111)          !                                             ! ===============
(  112)          ! 0. Matrix construction
(  113)          ! ----------------------
(  114) 
(  115)          ! Diagonal, inferior, superior (including the bottom boundary condition via avt masked)
(  116)          DO jk = 1, jpkm1
(  117)             DO ji = 2, jpim1
(  118)                zwi(ji,jk) = - r2dt(jk) * avt(ji,jj,jk  )   &
(  119)                                        / ( e3t(jk) * e3w(jk  ) )
(  120)                zws(ji,jk) = - r2dt(jk) * avt(ji,jj,jk+1)   &
(  121)                                        / ( e3t(jk) * e3w(jk+1) )
(  122)                zwd(ji,jk) = 1. - zwi(ji,jk) - zws(ji,jk)
(  123)             END DO
(  124)          END DO
(  125) 
(  126)          ! Surface boudary conditions
(  127)          DO ji = 2, jpim1
(  128)             zwi(ji,1) = 0.e0
(  129)             zwd(ji,1) = 1. - zws(ji,1)
(  130)          END DO
(  131) 
(  132)          ! 1. Vertical diffusion on temperature
(  133)          ! -------------------------===========
(  134) 
(  135)          ! Second member construction
(  136) # 146
(  146)          DO jk = 1, jpkm1
(  147)             DO ji = 2, jpim1             
(  148)                zwy(ji,jk) = tb(ji,jj,jk) + r2dt(jk) * ta(ji,jj,jk)
(  149)             END DO
(  150)          END DO
(  151) 
(  152) # 153






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 9

(  153)          ! Matrix inversion from the first level
(  154)          ikst = 1
(  155) 
(  156) # 1 "./zdf.matrixsolver.h90"
(    1) !! 
(    2) !!         ZDF.MATRIXSOLVER
(    3) !!       ********************
(    4) !! 
(    5) !! Matrix inversion
(    6) !!----------------------------------------------------------------------
(    7) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(    8) !
(    9) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(   10) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(   11) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(   12) !        (        ...               )( ...  ) ( ...  )
(   13) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(   14) !
(   15) !   m is decomposed in the product of an upper and lower triangular
(   16) !   matrix
(   17) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(   18) !   The second member is in 2d array zwy
(   19) !   The solution is in 2d array zwx
(   20) !   The 2d arry zwt and zwz are work space arrays
(   21) !
(   22) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(   23) !   the resolution of tke matrix where surface tke value is prescribed
(   24) !   so that ikstrt=2.
(   25) !!----------------------------------------------------------------------
(   26) !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   27) !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf.matrixsolver.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   28) !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   29) !!----------------------------------------------------------------------
(   30) 
(   31)       ikstp1 = ikst + 1
(   32)       ikenm2 = jpk - 2
(   33)       DO ji = 2, jpim1
(   34)          zwt(ji,ikst) = zwd(ji,ikst)
(   35)       END DO
(   36)       DO jk = ikstp1, jpkm1
(   37)          DO ji = 2, jpim1
(   38)             zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(   39)          END DO
(   40)       END DO
(   41)       DO ji = 2, jpim1
(   42)          zwz(ji,ikst) = zwy(ji,ikst)
(   43)       END DO
(   44)       DO jk = ikstp1, jpkm1
(   45)          DO ji = 2, jpim1
(   46)             zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(   47)          END DO
(   48)       END DO
(   49)       DO ji = 2, jpim1
(   50)          zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(   51)       END DO
(   52)       DO jk = ikenm2, ikst, -1
(   53)          DO ji = 2, jpim1
(   54)             zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 10

(   55)          END DO
(   56)       END DO
(   57)       
(   58) # 157 "trazdf_imp.F90"
(  157) 
(  158) # 158
(  158)          ! Save the masked temperature after in ta
(  159)          ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done
(  160)          !                  it will not be done in tranxt)
(  161)          DO jk = 1, jpkm1
(  162)             DO ji = 2, jpim1
(  163)                ta(ji,jj,jk) = zwx(ji,jk) * tmask(ji,jj,jk)
(  164)             END DO
(  165)          END DO
(  166) 
(  167) 
(  168)          ! 2. Vertical diffusion on salinity
(  169)          ! -------------------------========
(  170) 
(  171) # 192
(  192)          ! Second member construction
(  193) # 203
(  203)          DO jk = 1, jpkm1
(  204)             DO ji = 2, jpim1             
(  205)                zwy(ji,jk) = sb(ji,jj,jk) + r2dt(jk) * sa(ji,jj,jk)
(  206)             END DO
(  207)          END DO
(  208) # 209
(  209)  
(  210)          ! Matrix inversion from the first level
(  211)          ikst = 1
(  212) 
(  213) # 1 "./zdf.matrixsolver.h90"
(    1) !! 
(    2) !!         ZDF.MATRIXSOLVER
(    3) !!       ********************
(    4) !! 
(    5) !! Matrix inversion
(    6) !!----------------------------------------------------------------------
(    7) !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(    8) !
(    9) !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(   10) !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(   11) !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(   12) !        (        ...               )( ...  ) ( ...  )
(   13) !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(   14) !
(   15) !   m is decomposed in the product of an upper and lower triangular
(   16) !   matrix
(   17) !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(   18) !   The second member is in 2d array zwy
(   19) !   The solution is in 2d array zwx
(   20) !   The 2d arry zwt and zwz are work space arrays
(   21) !
(   22) !   N.B. the starting vertical index (ikst) is equal to 1 except for
(   23) !   the resolution of tke matrix where surface tke value is prescribed
(   24) !   so that ikstrt=2.
(   25) !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 11

(   26) !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   27) !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdf.matrixsolver.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   28) !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   29) !!----------------------------------------------------------------------
(   30) 
(   31)       ikstp1 = ikst + 1
(   32)       ikenm2 = jpk - 2
(   33)       DO ji = 2, jpim1
(   34)          zwt(ji,ikst) = zwd(ji,ikst)
(   35)       END DO
(   36)       DO jk = ikstp1, jpkm1
(   37)          DO ji = 2, jpim1
(   38)             zwt(ji,jk) = zwd(ji,jk) - zwi(ji,jk) * zws(ji,jk-1) / zwt(ji,jk-1)
(   39)          END DO
(   40)       END DO
(   41)       DO ji = 2, jpim1
(   42)          zwz(ji,ikst) = zwy(ji,ikst)
(   43)       END DO
(   44)       DO jk = ikstp1, jpkm1
(   45)          DO ji = 2, jpim1
(   46)             zwz(ji,jk) = zwy(ji,jk) - zwi(ji,jk) / zwt(ji,jk-1) * zwz(ji,jk-1)
(   47)          END DO
(   48)       END DO
(   49)       DO ji = 2, jpim1
(   50)          zwx(ji,jpkm1) = zwz(ji,jpkm1) / zwt(ji,jpkm1)
(   51)       END DO
(   52)       DO jk = ikenm2, ikst, -1
(   53)          DO ji = 2, jpim1
(   54)             zwx(ji,jk) =( zwz(ji,jk) - zws(ji,jk) * zwx(ji,jk+1) ) / zwt(ji,jk)
(   55)          END DO
(   56)       END DO
(   57)       
(   58) # 214 "trazdf_imp.F90"
(  214) 
(  215) 
(  216) # 216
(  216)          ! Save the masked salinity after in sa
(  217)          ! (c a u t i o n:  salinity not its trend, Leap-frog scheme done
(  218)          !                  it will not be done in tranxt)
(  219)          DO jk = 1, jpkm1
(  220)             DO ji = 2, jpim1
(  221)                sa(ji,jj,jk) = zwx(ji,jk)  * tmask(ji,jj,jk)
(  222)             END DO
(  223)          END DO
(  224) 
(  225)          !                                             ! ===============
(  226)       END DO                                           !   End of slab
(  227)       !                                                ! ===============
(  228) 
(  229)       ! save the trends for diagnostic
(  230)       ! Compute and save the vertical diffusive temperature & salinity trends
(  231)       IF( l_trdtra )   THEN
(  232)          ! compute the vertical diffusive trends in substracting the previous 
(  233)          ! trends ztdta()/ztdsa() to the new one computed (dT/dt or dS/dt) 
(  234)          ! with the new temperature/salinity ta/sa
(  235)          DO jk = 1, jpkm1
(  236)             ztdta(:,:,jk) = ( ( ta(:,:,jk) - tb(:,:,jk) ) / r2dt(jk) )   & ! new trend
(  237)                 &           - ztdta(:,:,jk)                                ! old trend






PGF90 (Version     10.5)          09/23/2010  09:50:17      page 12

(  238)             ztdsa(:,:,jk) = ( ( sa(:,:,jk) - sb(:,:,jk) ) / r2dt(jk) )   & ! new trend
(  239)                 &           - ztdsa(:,:,jk)                                ! old trend
(  240)          END DO
(  241) 
(  242)          CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
(  243)       ENDIF
(  244) 
(  245)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  246)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf  - Ta: ', mask1=tmask, &
(  247)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  248)       ENDIF
(  249) 
(  250)    END SUBROUTINE tra_zdf_imp
(  251) 
(  252)    !!==============================================================================
(  253) END MODULE trazdf_imp
