


PGF90 (Version     10.5)          09/23/2010  09:49:30      page 1

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

Filename: trazdf_exp.F90

(    1) MODULE trazdf_exp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trazdf_exp  ***
(    4)    !! Ocean active tracers:  vertical component of the tracer mixing trend using
(    5)    !!                        an explicit time-stepping (time spllitting scheme)
(    6)    !!==============================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   tra_zdf_exp  : update the tracer trend with the vertical diffusion
(   10)    !!                  using an explicit time stepping
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and active tracers 
(   14)    USE dom_oce         ! ocean space and time domain 
(   15)    USE trdmod          ! ocean active tracers trends 
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE zdf_oce         ! ocean vertical physics
(   18)    USE zdfddm          ! ocean vertical physics: double diffusion
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE prtctl          ! Print control
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC tra_zdf_exp          ! routine called by step.F90
(   27) 
(   28)    !! * Module variable
(   29)    REAL(wp), DIMENSION(jpk) ::   &
(   30)       r2dt                     ! vertical profile of 2 x tracer time-step
(   31) 
(   32)    !! * Substitutions
(   33) #  include "domzgr_substitute.h90"
(   34) #  include "zdfddm_substitute.h90"
(   35)    !!----------------------------------------------------------------------
(   36)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_exp.F90,v 1.4 2005/09/02 15:45:34 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE tra_zdf_exp( kt )
(   44)       !!----------------------------------------------------------------------
(   45)       !!                  ***  ROUTINE tra_zdf_exp  ***
(   46)       !!                   






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 2

(   47)       !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
(   48)       !!      using an explicit time stepping and add it to the general trend 
(   49)       !!      of the tracer equations.
(   50)       !!
(   51)       !! ** Method  :   The vertical diffusion of tracers (t & s) is given by:
(   52)       !!         difft = dz( avt dz(tb) ) = 1/e3t dk+1( avt/e3w dk(tb) )
(   53)       !!      It is evaluated with an Euler scheme, using a time splitting
(   54)       !!      technique.
(   55)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   56)       !!      both tracers (bottom, applied through the masked field avt).
(   57)       !!      Add this trend to the general trend ta,sa :
(   58)       !!          ta = ta + dz( avt dz(t) )
(   59)       !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm= T)
(   60)       !!
(   61)       !! ** Action : - Update (ta,sa) with the before vertical diffusion trend
(   62)       !!             - Save the trends  in (ztdta,ztdsa) ('key_trdtra')
(   63)       !!
(   64)       !! History :
(   65)       !!   6.0  !  90-10  (B. Blanke)  Original code
(   66)       !!   7.0  !  91-11  (G. Madec)
(   67)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(   68)       !!        !  96-01  (G. Madec)  statement function for e3
(   69)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   70)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(   71)       !!        !  00-08  (G. Madec)  double diffusive mixing
(   72)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   73)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   74)       !!---------------------------------------------------------------------
(   75)       !! * Modules used     
(   76)       USE oce, ONLY :    ztdta => ua,       & ! use ua as 3D workspace   
(   77)                          ztdsa => va          ! use va as 3D workspace   
(   78) 
(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   81)       
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk, jl             ! dummy loop indices
(   84)       REAL(wp) ::   &
(   85)          zlavmr,                            & ! temporary scalars
(   86)          zave3r, ze3tr,                     & !    "         "
(   87)          zta, zsa                             !    "         " 
(   88)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(   89)          zwx, zwy, zwz, zww
(   90)       !!---------------------------------------------------------------------
(   91) 
(   92) 
(   93)       ! 0. Local constant initialization
(   94)       ! --------------------------------
(   95)       ! time step = 2 rdttra 
(   96)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(   97)          r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
(   98)       ELSEIF( kt <= nit000 + 1) THEN
(   99)          r2dt(:) = 2. * rdttra(:)          ! leapfrog
(  100)       ENDIF
(  101)       zlavmr = 1. / float( n_zdfexp )
(  102) 
(  103)       ! Save ta and sa trends
(  104)       IF( l_trdtra )   THEN






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 3

(  105)          ztdta(:,:,:) = ta(:,:,:) 
(  106)          ztdsa(:,:,:) = sa(:,:,:) 
(  107)       ENDIF
(  108) 
(  109)       !                                                ! ===============
(  110)       DO jj = 2, jpjm1                                 !  Vertical slab
(  111)          !                                             ! ===============
(  112)          ! 1. Initializations
(  113)          ! ------------------
(  114) 
(  115)          ! Surface & bottom boundary conditions: no flux
(  116)          DO ji = 2, jpim1
(  117)             zwy(ji, 1 ) = 0.e0
(  118)             zwy(ji,jpk) = 0.e0
(  119)             zww(ji, 1 ) = 0.e0
(  120)             zww(ji,jpk) = 0.e0
(  121)          END DO
(  122) 
(  123)          ! zwx and zwz arrays set to before tracer values
(  124)          DO jk = 1, jpk
(  125)             DO ji = 2, jpim1
(  126)                zwx(ji,jk) = tb(ji,jj,jk)
(  127)                zwz(ji,jk) = sb(ji,jj,jk)
(  128)             END DO
(  129)          END DO
(  130) 
(  131) 
(  132)          ! 2. Time splitting loop
(  133)          ! ----------------------
(  134) 
(  135)          DO jl = 1, n_zdfexp
(  136) 
(  137)             ! first vertical derivative
(  138)             IF( lk_zdfddm ) THEN       ! double diffusion: avs /= avt
(  139)                DO jk = 2, jpk
(  140)                   DO ji = 2, jpim1
(  141)                      zave3r = 1.e0 / fse3w(ji,jj,jk) 
(  142)                      zwy(ji,jk) =   avt(ji,jj,jk) * ( zwx(ji,jk-1) - zwx(ji,jk) ) * zave3r
(  143)                      zww(ji,jk) = fsavs(ji,jj,jk) * ( zwz(ji,jk-1) - zwz(ji,jk) ) * zave3r
(  144)                   END DO
(  145)                END DO
(  146)             ELSE                      ! default : avs = avt
(  147)                DO jk = 2, jpk
(  148)                   DO ji = 2, jpim1
(  149)                      zave3r = avt(ji,jj,jk) / fse3w(ji,jj,jk)
(  150)                      zwy(ji,jk) = zave3r *(zwx(ji,jk-1) - zwx(ji,jk) )
(  151)                      zww(ji,jk) = zave3r *(zwz(ji,jk-1) - zwz(ji,jk) )
(  152)                   END DO
(  153)                END DO
(  154)             ENDIF
(  155) 
(  156)             ! trend estimation at kt+l*2*rdt/n_zdfexp
(  157)             DO jk = 1, jpkm1
(  158)                DO ji = 2, jpim1
(  159)                   ze3tr = zlavmr / fse3t(ji,jj,jk)
(  160)                   ! 2nd vertical derivative
(  161)                   zta = ( zwy(ji,jk) - zwy(ji,jk+1) ) * ze3tr
(  162)                   zsa = ( zww(ji,jk) - zww(ji,jk+1) ) * ze3tr






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 4

(  163)                   ! update the tracer trends
(  164)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  165)                   sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  166)                   ! update tracer fields at kt+l*2*rdt/n_zdfexp
(  167)                   zwx(ji,jk) = zwx(ji,jk) + r2dt(jk) * zta * tmask(ji,jj,jk)
(  168)                   zwz(ji,jk) = zwz(ji,jk) + r2dt(jk) * zsa * tmask(ji,jj,jk)
(  169)                END DO
(  170)             END DO
(  171)          END DO
(  172)          !                                             ! ===============
(  173)       END DO                                           !   End of slab
(  174)       !                                                ! ===============
(  175)       ! save the trends for diagnostic
(  176)       ! vertical diffusive tracers trends
(  177)       IF( l_trdtra )   THEN
(  178)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  179)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  180)          CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
(  181)       ENDIF
(  182) 
(  183)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  184)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf  - Ta: ', mask1=tmask, &
(  185)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  186)       ENDIF
(  187) 
(  188)    END SUBROUTINE tra_zdf_exp
(  189) 
(  190)    !!==============================================================================
(  191) END MODULE trazdf_exp



































PGF90 (Version     10.5)          09/23/2010  09:49:30      page 5

(    1) # 1 "trazdf_exp.F90"
(    1) MODULE trazdf_exp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trazdf_exp  ***
(    4)    !! Ocean active tracers:  vertical component of the tracer mixing trend using
(    5)    !!                        an explicit time-stepping (time spllitting scheme)
(    6)    !!==============================================================================
(    7) 
(    8)    !!----------------------------------------------------------------------
(    9)    !!   tra_zdf_exp  : update the tracer trend with the vertical diffusion
(   10)    !!                  using an explicit time stepping
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and active tracers 
(   14)    USE dom_oce         ! ocean space and time domain 
(   15)    USE trdmod          ! ocean active tracers trends 
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE zdf_oce         ! ocean vertical physics
(   18)    USE zdfddm          ! ocean vertical physics: double diffusion
(   19)    USE in_out_manager  ! I/O manager
(   20)    USE prtctl          ! Print control
(   21) 
(   22)    IMPLICIT NONE
(   23)    PRIVATE
(   24) 
(   25)    !! * Routine accessibility
(   26)    PUBLIC tra_zdf_exp          ! routine called by step.F90
(   27) 
(   28)    !! * Module variable
(   29)    REAL(wp), DIMENSION(jpk) ::   &
(   30)       r2dt                     ! vertical profile of 2 x tracer time-step
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
(   48) # 34 "trazdf_exp.F90"
(   34) # 1 "./zdfddm_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** zdfddm_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 6

(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfddm_substitute.h90,v 1.2 2005/03/27 18:35:25 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 16
(   16) !   Defautl option :                     avs = avt
(   17) # 35 "trazdf_exp.F90"
(   35) # 35
(   35)    !!----------------------------------------------------------------------
(   36)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/trazdf_exp.F90,v 1.4 2005/09/02 15:45:34 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE tra_zdf_exp( kt )
(   44)       !!----------------------------------------------------------------------
(   45)       !!                  ***  ROUTINE tra_zdf_exp  ***
(   46)       !!                   
(   47)       !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
(   48)       !!      using an explicit time stepping and add it to the general trend 
(   49)       !!      of the tracer equations.
(   50)       !!
(   51)       !! ** Method  :   The vertical diffusion of tracers (t & s) is given by:
(   52)       !!         difft = dz( avt dz(tb) ) = 1/e3t dk+1( avt/e3w dk(tb) )
(   53)       !!      It is evaluated with an Euler scheme, using a time splitting
(   54)       !!      technique.
(   55)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   56)       !!      both tracers (bottom, applied through the masked field avt).
(   57)       !!      Add this trend to the general trend ta,sa :
(   58)       !!          ta = ta + dz( avt dz(t) )
(   59)       !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm= T)
(   60)       !!
(   61)       !! ** Action : - Update (ta,sa) with the before vertical diffusion trend
(   62)       !!             - Save the trends  in (ztdta,ztdsa) ('key_trdtra')
(   63)       !!
(   64)       !! History :
(   65)       !!   6.0  !  90-10  (B. Blanke)  Original code
(   66)       !!   7.0  !  91-11  (G. Madec)
(   67)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(   68)       !!        !  96-01  (G. Madec)  statement function for e3
(   69)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   70)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(   71)       !!        !  00-08  (G. Madec)  double diffusive mixing
(   72)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   73)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   74)       !!---------------------------------------------------------------------
(   75)       !! * Modules used     
(   76)       USE oce, ONLY :    ztdta => ua,       & ! use ua as 3D workspace   
(   77)                          ztdsa => va          ! use va as 3D workspace   
(   78) 
(   79)       !! * Arguments
(   80)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   81)       
(   82)       !! * Local declarations
(   83)       INTEGER ::   ji, jj, jk, jl             ! dummy loop indices
(   84)       REAL(wp) ::   &
(   85)          zlavmr,                            & ! temporary scalars






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 7

(   86)          zave3r, ze3tr,                     & !    "         "
(   87)          zta, zsa                             !    "         " 
(   88)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(   89)          zwx, zwy, zwz, zww
(   90)       !!---------------------------------------------------------------------
(   91) 
(   92) 
(   93)       ! 0. Local constant initialization
(   94)       ! --------------------------------
(   95)       ! time step = 2 rdttra 
(   96)       IF( neuler == 0 .AND. kt == nit000 ) THEN
(   97)          r2dt(:) =  rdttra(:)              ! restarting with Euler time stepping
(   98)       ELSEIF( kt <= nit000 + 1) THEN
(   99)          r2dt(:) = 2. * rdttra(:)          ! leapfrog
(  100)       ENDIF
(  101)       zlavmr = 1. / float( n_zdfexp )
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
(  112)          ! 1. Initializations
(  113)          ! ------------------
(  114) 
(  115)          ! Surface & bottom boundary conditions: no flux
(  116)          DO ji = 2, jpim1
(  117)             zwy(ji, 1 ) = 0.e0
(  118)             zwy(ji,jpk) = 0.e0
(  119)             zww(ji, 1 ) = 0.e0
(  120)             zww(ji,jpk) = 0.e0
(  121)          END DO
(  122) 
(  123)          ! zwx and zwz arrays set to before tracer values
(  124)          DO jk = 1, jpk
(  125)             DO ji = 2, jpim1
(  126)                zwx(ji,jk) = tb(ji,jj,jk)
(  127)                zwz(ji,jk) = sb(ji,jj,jk)
(  128)             END DO
(  129)          END DO
(  130) 
(  131) 
(  132)          ! 2. Time splitting loop
(  133)          ! ----------------------
(  134) 
(  135)          DO jl = 1, n_zdfexp
(  136) 
(  137)             ! first vertical derivative
(  138)             IF( lk_zdfddm ) THEN       ! double diffusion: avs /= avt
(  139)                DO jk = 2, jpk
(  140)                   DO ji = 2, jpim1
(  141)                      zave3r = 1.e0 / e3w(jk) 
(  142)                      zwy(ji,jk) =   avt(ji,jj,jk) * ( zwx(ji,jk-1) - zwx(ji,jk) ) * zave3r
(  143)                      zww(ji,jk) = avt(ji,jj,jk) * ( zwz(ji,jk-1) - zwz(ji,jk) ) * zave3r






PGF90 (Version     10.5)          09/23/2010  09:49:30      page 8

(  144)                   END DO
(  145)                END DO
(  146)             ELSE                      ! default : avs = avt
(  147)                DO jk = 2, jpk
(  148)                   DO ji = 2, jpim1
(  149)                      zave3r = avt(ji,jj,jk) / e3w(jk)
(  150)                      zwy(ji,jk) = zave3r *(zwx(ji,jk-1) - zwx(ji,jk) )
(  151)                      zww(ji,jk) = zave3r *(zwz(ji,jk-1) - zwz(ji,jk) )
(  152)                   END DO
(  153)                END DO
(  154)             ENDIF
(  155) 
(  156)             ! trend estimation at kt+l*2*rdt/n_zdfexp
(  157)             DO jk = 1, jpkm1
(  158)                DO ji = 2, jpim1
(  159)                   ze3tr = zlavmr / e3t(jk)
(  160)                   ! 2nd vertical derivative
(  161)                   zta = ( zwy(ji,jk) - zwy(ji,jk+1) ) * ze3tr
(  162)                   zsa = ( zww(ji,jk) - zww(ji,jk+1) ) * ze3tr
(  163)                   ! update the tracer trends
(  164)                   ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  165)                   sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  166)                   ! update tracer fields at kt+l*2*rdt/n_zdfexp
(  167)                   zwx(ji,jk) = zwx(ji,jk) + r2dt(jk) * zta * tmask(ji,jj,jk)
(  168)                   zwz(ji,jk) = zwz(ji,jk) + r2dt(jk) * zsa * tmask(ji,jj,jk)
(  169)                END DO
(  170)             END DO
(  171)          END DO
(  172)          !                                             ! ===============
(  173)       END DO                                           !   End of slab
(  174)       !                                                ! ===============
(  175)       ! save the trends for diagnostic
(  176)       ! vertical diffusive tracers trends
(  177)       IF( l_trdtra )   THEN
(  178)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  179)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  180)          CALL trd_mod(ztdta, ztdsa, jpttdzdf, 'TRA', kt)
(  181)       ENDIF
(  182) 
(  183)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  184)          CALL prt_ctl(tab3d_1=ta, clinfo1=' zdf  - Ta: ', mask1=tmask, &
(  185)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  186)       ENDIF
(  187) 
(  188)    END SUBROUTINE tra_zdf_exp
(  189) 
(  190)    !!==============================================================================
(  191) END MODULE trazdf_exp
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (trazdf_exp.F90: 101)
