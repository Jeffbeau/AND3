


PGF90 (Version     10.5)          09/23/2010  09:50:13      page 1

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

Filename: traldf_bilap.F90

(    1) MODULE traldf_bilap
(    2)    !!==============================================================================
(    3)    !!                   ***  MODULE  traldf_bilap  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_ldf_bilap : update the tracer trend with the horizontal diffusion
(    9)    !!                   using a iso-level biharmonic operator
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE ldftra_oce      ! ocean tracer   lateral physics
(   15)    USE trdmod          ! ocean active tracers trends 
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE ldfslp          ! iso-neutral slopes 
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20)    USE diaptr          ! poleward transport diagnostics
(   21)    USE prtctl          ! Print control
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC tra_ldf_bilap   ! routine called by step.F90
(   28) 
(   29)    !! * Substitutions
(   30) #  include "domzgr_substitute.h90"
(   31) #  include "ldftra_substitute.h90"
(   32) #  include "ldfeiv_substitute.h90"
(   33) #  include "vectopt_loop_substitute.h90"
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_bilap.F90,v 1.6 2005/09/02 15:45:33 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41)    
(   42)    SUBROUTINE tra_ldf_bilap( kt )
(   43)       !!----------------------------------------------------------------------
(   44)       !!                  ***  ROUTINE tra_ldf_bilap  ***
(   45)       !!
(   46)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 2

(   47)       !!      trend and add it to the general trend of tracer equation.
(   48)       !!
(   49)       !! ** Method  :   4th order diffusive operator along model level surfaces 
(   50)       !!      evaluated using before fields (forward time scheme). The hor.
(   51)       !!      diffusive trends of temperature (idem for salinity) is given by:
(   52)       !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
(   53)       !!      factors e3. are inside the derivatives:
(   54)       !!      Laplacian of tb:
(   55)       !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(tb) ]
(   56)       !!                                  + dj-1[ e1v*e3v/e2v dj(tb) ]  }
(   57)       !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
(   58)       !!        zlt   = ahtt * zlt
(   59)       !!        call to lbc_lnk
(   60)       !!      Bilaplacian (laplacian of zlt):
(   61)       !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
(   62)       !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
(   63)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(   64)       !!      Laplacian of tb:
(   65)       !!         zlt   = 1/(e1t*e2t) {  di-1[ e2u/e1u di(tb) ]
(   66)       !!                              + dj-1[ e1v/e2v dj(tb) ] }
(   67)       !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
(   68)       !!        zlt   = ahtt * zlt
(   69)       !!        call to lbc_lnk
(   70)       !!      Bilaplacian (laplacian of zlt):
(   71)       !!         difft = 1/(e1t*e2t) {  di-1[ e2u/e1u di(zlt) ]
(   72)       !!                              + dj-1[ e1v/e2v dj(zlt) ]  }
(   73)       !!
(   74)       !!      Add this trend to the general trend (ta,sa):
(   75)       !!         (ta,sa) = (ta,sa) + ( difft , diffs )
(   76)       !!
(   77)       !! ** Action : - Update (ta,sa) arrays with the before iso-level
(   78)       !!               biharmonic mixing trend.
(   79)       !!             - Save the trends in (ztdta,ztdsa) ('key_trdtra')
(   80)       !!
(   81)       !! History :
(   82)       !!        !  91-11  (G. Madec)  Original code
(   83)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   84)       !!        !  95-11  (G. Madec)  suppress volumetric scale factors
(   85)       !!        !  96-01  (G. Madec)  statement function for e3
(   86)       !!        !  96-01  (M. Imbard)  mpp exchange
(   87)       !!        !  97-07  (G. Madec)  optimization, and ahtt
(   88)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   89)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   90)       !!----------------------------------------------------------------------
(   91)       !! * Modules used
(   92)       USE oce           , ztu => ua,  &  ! use ua as workspace
(   93)          &                ztv => va      ! use va as workspace
(   94) 
(   95)       !! * Arguments
(   96)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   97) 
(   98)       !! * Local declarations
(   99)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(  100) #if defined key_partial_steps
(  101)       INTEGER ::   iku, ikv               ! temporary integers
(  102) #endif
(  103)       REAL(wp) ::   zta, zsa              ! temporary scalars
(  104)       REAL(wp), DIMENSION(jpi,jpj) ::   & 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 3

(  105)          zeeu, zeev, zbtr,              & ! workspace
(  106)          zlt, zls
(  107)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   & 
(  108)          zsu, zsv,                          & ! workspace arrays
(  109)          ztdta, ztdsa
(  110)       !!----------------------------------------------------------------------
(  111) 
(  112)       IF( kt == nit000 ) THEN
(  113)          IF(lwp) WRITE(numout,*)
(  114)          IF(lwp) WRITE(numout,*) 'tra_ldf_bilap : iso-level biharmonic operator'
(  115)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
(  116)       ENDIF
(  117) 
(  118)       ! Save ta and sa trends
(  119)       IF( l_trdtra )   THEN
(  120)          ztdta(:,:,:) = ta(:,:,:) 
(  121)          ztdsa(:,:,:) = sa(:,:,:) 
(  122)       ENDIF
(  123) 
(  124)       !                                                ! ===============
(  125)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  126)          !                                             ! ===============
(  127) 
(  128)          ! 0. Initialization of metric arrays (for z- or s-coordinates)
(  129)          ! ----------------------------------
(  130) 
(  131)          DO jj = 1, jpjm1
(  132)             DO ji = 1, fs_jpim1   ! vector opt.
(  133) #if defined key_s_coord || defined key_partial_steps
(  134)                ! s-coordinates, vertical scale factor are used
(  135)                zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  136)                zeeu(ji,jj) = e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj) * umask(ji,jj,jk)
(  137)                zeev(ji,jj) = e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj) * vmask(ji,jj,jk)
(  138) #else
(  139)                ! z-coordinates, no vertical scale factors
(  140)                zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  141)                zeeu(ji,jj) = e2u(ji,jj) / e1u(ji,jj) * umask(ji,jj,jk)
(  142)                zeev(ji,jj) = e1v(ji,jj) / e2v(ji,jj) * vmask(ji,jj,jk)
(  143) #endif
(  144)             END DO
(  145)          END DO
(  146) 
(  147) 
(  148)          ! 1. Laplacian
(  149)          ! ------------
(  150) 
(  151)          ! First derivative (gradient)
(  152)          DO jj = 1, jpjm1
(  153)             DO ji = 1, fs_jpim1   ! vector opt.
(  154)                ztu(ji,jj,jk) = zeeu(ji,jj) * ( tb(ji+1,jj  ,jk) - tb(ji,jj,jk) )
(  155)                zsu(ji,jj,jk) = zeeu(ji,jj) * ( sb(ji+1,jj  ,jk) - sb(ji,jj,jk) )
(  156)                ztv(ji,jj,jk) = zeev(ji,jj) * ( tb(ji  ,jj+1,jk) - tb(ji,jj,jk) )
(  157)                zsv(ji,jj,jk) = zeev(ji,jj) * ( sb(ji  ,jj+1,jk) - sb(ji,jj,jk) )
(  158)             END DO
(  159)          END DO
(  160) #if defined key_partial_steps
(  161)          DO jj = 1, jpj-1
(  162)             DO ji = 1, jpi-1






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 4

(  163)                ! last level
(  164)                iku = MIN ( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
(  165)                ikv = MIN ( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
(  166)                IF( iku == jk ) THEN
(  167)                   ztu(ji,jj,jk) = zeeu(ji,jj) * gtu(ji,jj)
(  168)                   zsu(ji,jj,jk) = zeeu(ji,jj) * gsu(ji,jj)
(  169)                ENDIF
(  170)                IF( ikv == jk ) THEN
(  171)                   ztv(ji,jj,jk) = zeev(ji,jj) * gtv(ji,jj)
(  172)                   zsv(ji,jj,jk) = zeev(ji,jj) * gsv(ji,jj)
(  173)                ENDIF
(  174)             END DO
(  175)          END DO
(  176) #endif
(  177) 
(  178)          ! Second derivative (divergence)
(  179)          DO jj = 2, jpjm1
(  180)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  181)                zlt(ji,jj) = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
(  182)                zls(ji,jj) = zbtr(ji,jj) * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk) + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
(  183)             END DO
(  184)          END DO
(  185) 
(  186)          ! Multiply by the eddy diffusivity coefficient
(  187)          DO jj = 2, jpjm1
(  188)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  189)                zlt(ji,jj) = fsahtt(ji,jj,jk) * zlt(ji,jj)
(  190)                zls(ji,jj) = fsahtt(ji,jj,jk) * zls(ji,jj)
(  191)             END DO
(  192)          END DO
(  193) 
(  194)          ! Lateral boundary conditions on the laplacian (zlt,zls)   (unchanged sgn)
(  195)          CALL lbc_lnk( zlt, 'T', 1. )   ;    CALL lbc_lnk( zls, 'T', 1. )
(  196) 
(  197)          ! 2. Bilaplacian
(  198)          ! --------------
(  199) 
(  200)          ! third derivative (gradient)
(  201)          DO jj = 1, jpjm1
(  202)             DO ji = 1, fs_jpim1   ! vector opt.
(  203)                ztu(ji,jj,jk) = zeeu(ji,jj) * ( zlt(ji+1,jj  ) - zlt(ji,jj) )
(  204)                zsu(ji,jj,jk) = zeeu(ji,jj) * ( zls(ji+1,jj  ) - zls(ji,jj) )
(  205)                ztv(ji,jj,jk) = zeev(ji,jj) * ( zlt(ji  ,jj+1) - zlt(ji,jj) )
(  206)                zsv(ji,jj,jk) = zeev(ji,jj) * ( zls(ji  ,jj+1) - zls(ji,jj) )
(  207)             END DO
(  208)          END DO
(  209) 
(  210)          ! fourth derivative (divergence) and add to the general tracer trend
(  211)          DO jj = 2, jpjm1
(  212)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  213)                ! horizontal diffusive trends
(  214)                zta = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
(  215)                zsa = zbtr(ji,jj) * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk) + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
(  216)                ! add it to the general tracer trends
(  217)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  218)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  219)             END DO
(  220)          END DO






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 5

(  221)          !                                             ! ===============
(  222)       END DO                                           ! Horizontal slab
(  223)       !                                                ! ===============
(  224) 
(  225)       ! save the trends for diagnostic
(  226)       ! save the horizontal diffusive trends
(  227)       IF( l_trdtra )   THEN
(  228)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  229)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  230) 
(  231)          CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
(  232)       ENDIF
(  233) 
(  234)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  235)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  236)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  237)       ENDIF
(  238) 
(  239)       ! "zonal" mean lateral diffusive heat and salt transport
(  240)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  241) # if defined key_s_coord || defined key_partial_steps
(  242)          pht_ldf(:) = ptr_vj( ztv(:,:,:) )
(  243)          pst_ldf(:) = ptr_vj( zsv(:,:,:) )
(  244) # else
(  245)          DO jk = 1, jpkm1
(  246)             DO jj = 2, jpjm1
(  247)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  248)                   ztv(ji,jj,jk) = ztv(ji,jj,jk) * fse3v(ji,jj,jk)
(  249)                   zsv(ji,jj,jk) = zsv(ji,jj,jk) * fse3v(ji,jj,jk)
(  250)                END DO
(  251)             END DO
(  252)          END DO
(  253)          pht_ldf(:) = ptr_vj( ztv(:,:,:) )
(  254)          pst_ldf(:) = ptr_vj( zsv(:,:,:) )
(  255) # endif
(  256)       ENDIF
(  257) 
(  258)    END SUBROUTINE tra_ldf_bilap
(  259) 
(  260)    !!==============================================================================
(  261) END MODULE traldf_bilap























PGF90 (Version     10.5)          09/23/2010  09:50:13      page 6

(    1) # 1 "traldf_bilap.F90"
(    1) MODULE traldf_bilap
(    2)    !!==============================================================================
(    3)    !!                   ***  MODULE  traldf_bilap  ***
(    4)    !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_ldf_bilap : update the tracer trend with the horizontal diffusion
(    9)    !!                   using a iso-level biharmonic operator
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce             ! ocean dynamics and active tracers
(   13)    USE dom_oce         ! ocean space and time domain
(   14)    USE ldftra_oce      ! ocean tracer   lateral physics
(   15)    USE trdmod          ! ocean active tracers trends 
(   16)    USE trdmod_oce      ! ocean variables trends
(   17)    USE in_out_manager  ! I/O manager
(   18)    USE ldfslp          ! iso-neutral slopes 
(   19)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   20)    USE diaptr          ! poleward transport diagnostics
(   21)    USE prtctl          ! Print control
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC tra_ldf_bilap   ! routine called by step.F90
(   28) 
(   29)    !! * Substitutions
(   30) # 1 "./domzgr_substitute.h90"
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
(   48) # 31 "traldf_bilap.F90"
(   31) # 1 "./ldftra_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    *** ldftra_substitute.h90  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
(    5)    !!      with a constant or 1D or 2D or 3D array, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/LDF/ldftra_substitute.h90,v 1.2 2005/03/27 18:35:07 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 7

(   12) # 19
(   19) !   'key_traldf_c3d' :                 aht: 3D coefficient
(   20) # 32 "traldf_bilap.F90"
(   32) # 1 "./ldfeiv_substitute.h90"
(    1) # 40 "./ldfeiv_substitute.h90"
(   40) !   Default option :                         eiv: dummy variables
(   41) # 33 "traldf_bilap.F90"
(   33) # 1 "./vectopt_loop_substitute.h90"
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
(   12) # 34 "traldf_bilap.F90"
(   34) # 34
(   34)    !!----------------------------------------------------------------------
(   35)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   36)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/TRA/traldf_bilap.F90,v 1.6 2005/09/02 15:45:33 opalod Exp $ 
(   37)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   38)    !!----------------------------------------------------------------------
(   39) 
(   40) CONTAINS
(   41)    
(   42)    SUBROUTINE tra_ldf_bilap( kt )
(   43)       !!----------------------------------------------------------------------
(   44)       !!                  ***  ROUTINE tra_ldf_bilap  ***
(   45)       !!
(   46)       !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
(   47)       !!      trend and add it to the general trend of tracer equation.
(   48)       !!
(   49)       !! ** Method  :   4th order diffusive operator along model level surfaces 
(   50)       !!      evaluated using before fields (forward time scheme). The hor.
(   51)       !!      diffusive trends of temperature (idem for salinity) is given by:
(   52)       !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
(   53)       !!      factors e3. are inside the derivatives:
(   54)       !!      Laplacian of tb:
(   55)       !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(tb) ]
(   56)       !!                                  + dj-1[ e1v*e3v/e2v dj(tb) ]  }
(   57)       !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
(   58)       !!        zlt   = ahtt * zlt
(   59)       !!        call to lbc_lnk
(   60)       !!      Bilaplacian (laplacian of zlt):
(   61)       !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
(   62)       !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
(   63)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(   64)       !!      Laplacian of tb:
(   65)       !!         zlt   = 1/(e1t*e2t) {  di-1[ e2u/e1u di(tb) ]
(   66)       !!                              + dj-1[ e1v/e2v dj(tb) ] }
(   67)       !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
(   68)       !!        zlt   = ahtt * zlt
(   69)       !!        call to lbc_lnk
(   70)       !!      Bilaplacian (laplacian of zlt):






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 8

(   71)       !!         difft = 1/(e1t*e2t) {  di-1[ e2u/e1u di(zlt) ]
(   72)       !!                              + dj-1[ e1v/e2v dj(zlt) ]  }
(   73)       !!
(   74)       !!      Add this trend to the general trend (ta,sa):
(   75)       !!         (ta,sa) = (ta,sa) + ( difft , diffs )
(   76)       !!
(   77)       !! ** Action : - Update (ta,sa) arrays with the before iso-level
(   78)       !!               biharmonic mixing trend.
(   79)       !!             - Save the trends in (ztdta,ztdsa) ('key_trdtra')
(   80)       !!
(   81)       !! History :
(   82)       !!        !  91-11  (G. Madec)  Original code
(   83)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   84)       !!        !  95-11  (G. Madec)  suppress volumetric scale factors
(   85)       !!        !  96-01  (G. Madec)  statement function for e3
(   86)       !!        !  96-01  (M. Imbard)  mpp exchange
(   87)       !!        !  97-07  (G. Madec)  optimization, and ahtt
(   88)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   89)       !!   9.0  !  04-08  (C. Talandier) New trends organization
(   90)       !!----------------------------------------------------------------------
(   91)       !! * Modules used
(   92)       USE oce           , ztu => ua,  &  ! use ua as workspace
(   93)          &                ztv => va      ! use va as workspace
(   94) 
(   95)       !! * Arguments
(   96)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   97) 
(   98)       !! * Local declarations
(   99)       INTEGER ::   ji, jj, jk             ! dummy loop indices
(  100) # 103
(  103)       REAL(wp) ::   zta, zsa              ! temporary scalars
(  104)       REAL(wp), DIMENSION(jpi,jpj) ::   & 
(  105)          zeeu, zeev, zbtr,              & ! workspace
(  106)          zlt, zls
(  107)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   & 
(  108)          zsu, zsv,                          & ! workspace arrays
(  109)          ztdta, ztdsa
(  110)       !!----------------------------------------------------------------------
(  111) 
(  112)       IF( kt == nit000 ) THEN
(  113)          IF(lwp) WRITE(numout,*)
(  114)          IF(lwp) WRITE(numout,*) 'tra_ldf_bilap : iso-level biharmonic operator'
(  115)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
(  116)       ENDIF
(  117) 
(  118)       ! Save ta and sa trends
(  119)       IF( l_trdtra )   THEN
(  120)          ztdta(:,:,:) = ta(:,:,:) 
(  121)          ztdsa(:,:,:) = sa(:,:,:) 
(  122)       ENDIF
(  123) 
(  124)       !                                                ! ===============
(  125)       DO jk = 1, jpkm1                                 ! Horizontal slab
(  126)          !                                             ! ===============
(  127) 
(  128)          ! 0. Initialization of metric arrays (for z- or s-coordinates)
(  129)          ! ----------------------------------
(  130) 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 9

(  131)          DO jj = 1, jpjm1
(  132)             DO ji = 1, jpim1   ! vector opt.
(  133) # 139
(  139)                ! z-coordinates, no vertical scale factors
(  140)                zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  141)                zeeu(ji,jj) = e2u(ji,jj) / e1u(ji,jj) * umask(ji,jj,jk)
(  142)                zeev(ji,jj) = e1v(ji,jj) / e2v(ji,jj) * vmask(ji,jj,jk)
(  143) # 144
(  144)             END DO
(  145)          END DO
(  146) 
(  147) 
(  148)          ! 1. Laplacian
(  149)          ! ------------
(  150) 
(  151)          ! First derivative (gradient)
(  152)          DO jj = 1, jpjm1
(  153)             DO ji = 1, jpim1   ! vector opt.
(  154)                ztu(ji,jj,jk) = zeeu(ji,jj) * ( tb(ji+1,jj  ,jk) - tb(ji,jj,jk) )
(  155)                zsu(ji,jj,jk) = zeeu(ji,jj) * ( sb(ji+1,jj  ,jk) - sb(ji,jj,jk) )
(  156)                ztv(ji,jj,jk) = zeev(ji,jj) * ( tb(ji  ,jj+1,jk) - tb(ji,jj,jk) )
(  157)                zsv(ji,jj,jk) = zeev(ji,jj) * ( sb(ji  ,jj+1,jk) - sb(ji,jj,jk) )
(  158)             END DO
(  159)          END DO
(  160) 
(  161) # 178
(  178)          ! Second derivative (divergence)
(  179)          DO jj = 2, jpjm1
(  180)             DO ji = 2, jpim1   ! vector opt.
(  181)                zlt(ji,jj) = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
(  182)                zls(ji,jj) = zbtr(ji,jj) * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk) + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
(  183)             END DO
(  184)          END DO
(  185) 
(  186)          ! Multiply by the eddy diffusivity coefficient
(  187)          DO jj = 2, jpjm1
(  188)             DO ji = 2, jpim1   ! vector opt.
(  189)                zlt(ji,jj) = ahtt(ji,jj,jk) * zlt(ji,jj)
(  190)                zls(ji,jj) = ahtt(ji,jj,jk) * zls(ji,jj)
(  191)             END DO
(  192)          END DO
(  193) 
(  194)          ! Lateral boundary conditions on the laplacian (zlt,zls)   (unchanged sgn)
(  195)          CALL lbc_lnk( zlt, 'T', 1. )   ;    CALL lbc_lnk( zls, 'T', 1. )
(  196) 
(  197)          ! 2. Bilaplacian
(  198)          ! --------------
(  199) 
(  200)          ! third derivative (gradient)
(  201)          DO jj = 1, jpjm1
(  202)             DO ji = 1, jpim1   ! vector opt.
(  203)                ztu(ji,jj,jk) = zeeu(ji,jj) * ( zlt(ji+1,jj  ) - zlt(ji,jj) )
(  204)                zsu(ji,jj,jk) = zeeu(ji,jj) * ( zls(ji+1,jj  ) - zls(ji,jj) )
(  205)                ztv(ji,jj,jk) = zeev(ji,jj) * ( zlt(ji  ,jj+1) - zlt(ji,jj) )
(  206)                zsv(ji,jj,jk) = zeev(ji,jj) * ( zls(ji  ,jj+1) - zls(ji,jj) )
(  207)             END DO
(  208)          END DO
(  209) 






PGF90 (Version     10.5)          09/23/2010  09:50:13      page 10

(  210)          ! fourth derivative (divergence) and add to the general tracer trend
(  211)          DO jj = 2, jpjm1
(  212)             DO ji = 2, jpim1   ! vector opt.
(  213)                ! horizontal diffusive trends
(  214)                zta = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
(  215)                zsa = zbtr(ji,jj) * (  zsu(ji,jj,jk) - zsu(ji-1,jj,jk) + zsv(ji,jj,jk) - zsv(ji,jj-1,jk)  )
(  216)                ! add it to the general tracer trends
(  217)                ta(ji,jj,jk) = ta(ji,jj,jk) + zta
(  218)                sa(ji,jj,jk) = sa(ji,jj,jk) + zsa
(  219)             END DO
(  220)          END DO
(  221)          !                                             ! ===============
(  222)       END DO                                           ! Horizontal slab
(  223)       !                                                ! ===============
(  224) 
(  225)       ! save the trends for diagnostic
(  226)       ! save the horizontal diffusive trends
(  227)       IF( l_trdtra )   THEN
(  228)          ztdta(:,:,:) = ta(:,:,:) - ztdta(:,:,:)
(  229)          ztdsa(:,:,:) = sa(:,:,:) - ztdsa(:,:,:)
(  230) 
(  231)          CALL trd_mod(ztdta, ztdsa, jpttdldf, 'TRA', kt)
(  232)       ENDIF
(  233) 
(  234)       IF(ln_ctl) THEN         ! print mean trends (used for debugging)
(  235)          CALL prt_ctl(tab3d_1=ta, clinfo1=' ldf  - Ta: ', mask1=tmask, &
(  236)             &         tab3d_2=sa, clinfo2=' Sa: ', mask2=tmask, clinfo3='tra')
(  237)       ENDIF
(  238) 
(  239)       ! "zonal" mean lateral diffusive heat and salt transport
(  240)       IF( ln_diaptr .AND. ( MOD( kt, nf_ptr ) == 0 ) ) THEN
(  241) # 245
(  245)          DO jk = 1, jpkm1
(  246)             DO jj = 2, jpjm1
(  247)                DO ji = 2, jpim1   ! vector opt.
(  248)                   ztv(ji,jj,jk) = ztv(ji,jj,jk) * e3t(jk)
(  249)                   zsv(ji,jj,jk) = zsv(ji,jj,jk) * e3t(jk)
(  250)                END DO
(  251)             END DO
(  252)          END DO
(  253)          pht_ldf(:) = ptr_vj( ztv(:,:,:) )
(  254)          pst_ldf(:) = ptr_vj( zsv(:,:,:) )
(  255) # 256
(  256)       ENDIF
(  257) 
(  258)    END SUBROUTINE tra_ldf_bilap
(  259) 
(  260)    !!==============================================================================
(  261) END MODULE traldf_bilap
